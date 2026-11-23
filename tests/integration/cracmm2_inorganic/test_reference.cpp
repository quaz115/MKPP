#include "mkpp_host/cracmm2_inorganic.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>

extern "C" void mkpp_reference_lite(const double* totals, const double* organic, double rh,
                                    double temperature, double* gases, double* particles,
                                    double* sulfate, double* water, double* hydronium, int* status);

namespace {
struct Fixture {
  std::array<double, 8> component_mass;  // ug/m3 as reported by the original input deck
  double organic_mass;                  // ug/m3
  double rh;
  double temperature;
};

// First record of each named input file supplied with the original Lite distribution.
constexpr Fixture ammonium_nitrate{{0, 1, 2, 1, 0, 0, 0, 0}, 10, 0.80, 298.15};  // test1.inp
constexpr Fixture mixed_ions{{0.54, 1.088, 0.37, 0.087, 0.87, 0.01, 0.006, 0.004},
                             1.508, 0.83, 291.89};  // Partitioning_with_organics.INP

// Ammonia in excess of sulfate neutralization with an appreciable total nitrate pool.
// No organic mass is available to supply organic water in this cold case.
constexpr Fixture cold_inorganic{{0, 1, 2, 3.15, 0, 0, 0, 0}, 0, 0.70, 278.15};
// 0.060 umol/m3 total NHx and 0.040 umol/m3 sulfate, below full neutralization.
constexpr Fixture sulfate_rich{{0, 3.92, 1.02, 1.26, 0, 0, 0, 0}, 2, 0.55, 288.15};

bool compare(std::string_view field, std::string_view unit, double actual, double reference,
             double absolute_tolerance) {
  // Same-solver port comparison. Floors correspond to 5 pmol/m3 air or 0.1 ng water/m3.
  constexpr double relative_tolerance = 0.01;
  const double tolerance = absolute_tolerance + relative_tolerance * std::abs(reference);
  const double difference = std::abs(actual - reference);
  const bool passed = std::isfinite(actual) && std::isfinite(reference) && actual >= 0.0 &&
                      reference >= 0.0 && difference <= tolerance;
  std::cout << std::left << std::setw(16) << field << std::setw(10) << unit << std::right << std::scientific
            << std::setprecision(6) << std::setw(16) << actual << std::setw(16) << reference
            << std::setw(16) << difference << std::setw(16) << tolerance << (passed ? "  PASS\n" : "  FAIL\n");
  return passed;
}

int run(std::string_view name, const Fixture& fixture) {
  std::cout << "Case " << name << ", T=" << fixture.temperature << " K, RH=" << fixture.rh << '\n';
  // Input deck mass basis is Na, H2SO4, NH3, HNO3, HCl, Ca, K and Mg (g/mol).
  constexpr std::array<double, 8> molar_mass{23, 98, 17, 63, 36.5, 40.1, 39.1, 24.3};
  constexpr double molecules_cm3_per_mol_m3 = 6.02214076e17;
  std::array<double, 8> totals{};
  for (std::size_t i = 0; i < totals.size(); ++i) {
    totals[i] = fixture.component_mass[i] * 1.0e-6 / molar_mass[i];
  }
  const std::array<double, 3> organic{fixture.organic_mass * 1.0e-9, 0.15, 1000.0};
  std::array<double, 3> gases{}, particles{};
  double sulfate{};
  double hydronium{};
  std::array<double, 2> water{};
  int reference_status = -1;
  mkpp_reference_lite(totals.data(), organic.data(), fixture.rh, fixture.temperature, gases.data(),
                      particles.data(), &sulfate, water.data(), &hydronium, &reference_status);
  const bool expected_organic_water = fixture.organic_mass > 0.0 ? water[1] > 0.0 : water[1] == 0.0;
  if (reference_status != 0 || !(water[0] > 0.0) || !expected_organic_water ||
      !(particles[0] > 0.0) || !(sulfate > 0.0)) {
    std::cerr << "Legacy Lite reference failed to produce a wet equilibrium, status=" << reference_status << '\n';
    return 1;
  }

  using namespace mkpp::host;
  Cracmm2InorganicState host{};
  host.ana = totals[0] * molecules_cm3_per_mol_m3;
  host.aso4 = totals[1] * molecules_cm3_per_mol_m3;
  host.nh3 = totals[2] * molecules_cm3_per_mol_m3;
  host.hno3 = totals[3] * molecules_cm3_per_mol_m3;
  host.acl = totals[4] * molecules_cm3_per_mol_m3;
  host.aca = totals[5] * molecules_cm3_per_mol_m3;
  host.ak = totals[6] * molecules_cm3_per_mol_m3;
  host.amg = totals[7] * molecules_cm3_per_mol_m3;
  const Cracmm2ThermoEnvironment environment{fixture.temperature, fixture.rh, organic[0], organic[1], organic[2]};
  Cracmm2ThermoDiagnostics host_water{};
  if (partition_cracmm2_inorganic(host, environment, host_water) != Cracmm2ThermoStatus::success) {
    std::cerr << "Host thermodynamic solve failed\n";
    return 1;
  }

  constexpr double molar_tolerance = 5.0e-12;  // mol/m3
  constexpr double water_tolerance = 1.0e-13;   // kg/m3
  const std::array<double, 3> host_gases{host.nh3, host.hno3, host.hcl};
  const std::array<double, 3> host_particles{host.anh4, host.ano3, host.acl};
  constexpr std::array<std::string_view, 3> gas_names{"NH3", "HNO3", "HCl"};
  constexpr std::array<std::string_view, 3> particle_names{"ANH4", "ANO3", "ACL"};
  std::cout << std::left << std::setw(16) << "Field" << std::setw(10) << "Unit" << std::right
            << std::setw(16) << "MKPP" << std::setw(16) << "Legacy Lite" << std::setw(16) << "Abs. diff."
            << std::setw(16) << "Tolerance" << "  Result\n";
  bool ok = true;
  for (std::size_t i = 0; i < gases.size(); ++i) {
    ok &= compare(gas_names[i], "mol m-3", host_gases[i] / molecules_cm3_per_mol_m3, gases[i], molar_tolerance);
    ok &= compare(particle_names[i], "mol m-3", host_particles[i] / molecules_cm3_per_mol_m3,
                  particles[i], molar_tolerance);
  }
  ok &= compare("ASO4", "mol m-3", host.aso4 / molecules_cm3_per_mol_m3, sulfate, molar_tolerance);
  ok &= compare("AH3OP", "mol m-3", host_water.hydronium_cm3 / molecules_cm3_per_mol_m3, hydronium, 1.0e-16);
  ok &= compare("total water", "kg m-3", host_water.total_kg_m3, water[0], water_tolerance);
  ok &= compare("organic water", "kg m-3", host_water.organic_kg_m3, water[1], water_tolerance);
  return ok ? 0 : 1;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc == 2 && std::string_view(argv[1]) == "ammonium-nitrate") return run(argv[1], ammonium_nitrate);
  if (argc == 2 && std::string_view(argv[1]) == "mixed-ions") return run(argv[1], mixed_ions);
  if (argc == 2 && std::string_view(argv[1]) == "cold-inorganic") return run(argv[1], cold_inorganic);
  if (argc == 2 && std::string_view(argv[1]) == "sulfate-rich") return run(argv[1], sulfate_rich);
  std::cerr << "Expected a case name: ammonium-nitrate, mixed-ions, cold-inorganic, or sulfate-rich\n";
  return 2;
}
