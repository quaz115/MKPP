#include "mkpp_host/cracmm2_inorganic.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>

extern "C" void mkpp_reference_lite(const double* totals, const double* organic, double rh,
                                    double temperature, double* gases, double* particles,
                                    double* water, int* status);

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

bool compare(std::string_view field, double actual, double reference, double absolute_tolerance) {
  // Same-solver port comparison. Floors correspond to 5 pmol/m3 air or 0.1 ng water/m3.
  constexpr double relative_tolerance = 0.01;
  const double tolerance = absolute_tolerance + relative_tolerance * std::abs(reference);
  if (std::isfinite(actual) && std::isfinite(reference) && actual >= 0.0 && reference >= 0.0 &&
      std::abs(actual - reference) <= tolerance) {
    return true;
  }
  std::cerr << std::setprecision(12) << field << ": host=" << actual << " reference=" << reference
            << " tolerance=" << tolerance << '\n';
  return false;
}

int run(const Fixture& fixture) {
  // Input deck mass basis is Na, H2SO4, NH3, HNO3, HCl, Ca, K and Mg (g/mol).
  constexpr std::array<double, 8> molar_mass{23, 98, 17, 63, 36.5, 40.1, 39.1, 24.3};
  constexpr double molecules_cm3_per_mol_m3 = 6.02214076e17;
  std::array<double, 8> totals{};
  for (std::size_t i = 0; i < totals.size(); ++i) {
    totals[i] = fixture.component_mass[i] * 1.0e-6 / molar_mass[i];
  }
  const std::array<double, 3> organic{fixture.organic_mass * 1.0e-9, 0.15, 1000.0};
  std::array<double, 3> gases{}, particles{};
  std::array<double, 2> water{};
  int reference_status = -1;
  mkpp_reference_lite(totals.data(), organic.data(), fixture.rh, fixture.temperature, gases.data(),
                      particles.data(), water.data(), &reference_status);
  if (reference_status != 0 || !(water[0] > 0.0) || !(water[1] > 0.0) || !(particles[0] > 0.0)) {
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
  Cracmm2AerosolWater host_water{};
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
  bool ok = true;
  for (std::size_t i = 0; i < gases.size(); ++i) {
    ok &= compare(gas_names[i], host_gases[i] / molecules_cm3_per_mol_m3, gases[i], molar_tolerance);
    ok &= compare(particle_names[i], host_particles[i] / molecules_cm3_per_mol_m3, particles[i], molar_tolerance);
  }
  ok &= compare("total water", host_water.total_kg_m3, water[0], water_tolerance);
  ok &= compare("organic water", host_water.organic_kg_m3, water[1], water_tolerance);
  return ok ? 0 : 1;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc == 2 && std::string_view(argv[1]) == "ammonium-nitrate") return run(ammonium_nitrate);
  if (argc == 2 && std::string_view(argv[1]) == "mixed-ions") return run(mixed_ions);
  std::cerr << "Expected a case name: ammonium-nitrate or mixed-ions\n";
  return 2;
}
