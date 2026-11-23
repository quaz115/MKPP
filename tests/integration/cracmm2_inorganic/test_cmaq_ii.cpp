#include "mkpp_host/cracmm2_inorganic.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>

extern "C" void mkpp_reference_isorropia_ii(const double* totals, double rh, double temperature, double* gases,
                                       double* particles, double* sulfate, double* water, double* hydronium, int* status);

namespace {
constexpr std::string_view reference_name = MKPP_II_REFERENCE_NAME;
constexpr double to_host = 6.02214076e17; // mol/m3 air to molecules/cm3 air
constexpr double rh = 0.70;
// Same inorganic totals and RH in both cases. Only temperature changes.
// Original input-deck basis: 1, 2 and 3.15 ug/m3 of H2SO4, NH3 and HNO3 equivalents.
constexpr std::array<double, 8> totals{0.0, 1e-6/98, 2e-6/17, 3.15e-6/63, 0.0, 0.0, 0.0, 0.0};

struct Result {
  double nh3, hno3, nh4, no3, sulfate, water, hydronium;
};

bool physical(const Result& r, std::string_view model) {
  bool ok = true;
  for (double v : {r.nh3, r.hno3, r.nh4, r.no3, r.sulfate, r.water, r.hydronium}) {
    if (!std::isfinite(v) || v < 0.0) ok = false;
  }
  // Independent input-versus-sum checks, including every form returned by the II reference.
  const std::array<double, 3> final_totals{r.sulfate, r.nh3 + r.nh4, r.hno3 + r.no3};
  constexpr std::array<std::string_view, 3> labels{"sulfate", "NHx", "nitrate"};
  for (std::size_t i = 0; i < final_totals.size(); ++i) {
    const double residual = final_totals[i] - totals[i+1];
    const double tolerance = 5.0e-12 + 1.0e-3 * totals[i+1]; // mol/m3, plus 0.1% of total
    std::cout << model << ' ' << labels[i] << " closure residual=" << residual << " mol m-3, limit=" << tolerance << '\n';
    if (std::abs(residual) > tolerance) ok = false;
  }
  if (!(r.water > 0.0)) ok = false; // these tests select wet, metastable equilibria
  if (!ok) std::cerr << model << " returned an invalid or unbalanced wet equilibrium\n";
  return ok;
}

void report(std::string_view name, std::string_view unit, double lite, double ii) {
  std::cout << std::left << std::setw(22) << name << std::setw(12) << unit << std::right
            << std::setw(17) << lite << std::setw(17) << ii << std::setw(17) << (lite - ii) << '\n';
}

bool check_particle_fraction(std::string_view species, double lite, double ii, double total) {
  const double lite_fraction = lite / total;
  const double ii_fraction = ii / total;
  report(species, "1", lite_fraction, ii_fraction);
  // Regression limit for the two zero-organic NHx/sulfate/nitrate fixtures.
  constexpr double tolerance = 0.005; // absolute fraction; 0.5 percentage point
  if (std::abs(lite_fraction - ii_fraction) <= tolerance) return true;
  std::cerr << species << " exceeds the absolute fraction tolerance " << tolerance << '\n';
  return false;
}

int run(double temperature) {
  using namespace mkpp::host;
  std::cout << reference_name << " and MKPP Lite, forward/metastable, zero organics, T="
            << temperature << " K, RH=" << rh << '\n';
  Cracmm2InorganicState state{};
  state.aso4 = totals[1] * to_host;
  state.nh3 = totals[2] * to_host;
  state.hno3 = totals[3] * to_host;
  Cracmm2ThermoDiagnostics diagnostics{};
  if (partition_cracmm2_inorganic(state, {temperature, rh, 0.0, 0.0, 1000.0}, diagnostics)
      != Cracmm2ThermoStatus::success || diagnostics.organic_kg_m3 != 0.0) {
    std::cerr << "Invalid inorganic-only MKPP Lite result\n";
    return 1;
  }
  const Result lite{state.nh3/to_host, state.hno3/to_host, state.anh4/to_host, state.ano3/to_host,
                    state.aso4/to_host, diagnostics.total_kg_m3, diagnostics.hydronium_cm3/to_host};
  std::array<double, 3> gas{};
  std::array<double, 2> particle{};
  double sulfate{}, water{}, hydronium{};
  int status = -1;
  mkpp_reference_isorropia_ii(totals.data(), rh, temperature, gas.data(), particle.data(), &sulfate, &water, &hydronium, &status);
  if (status != 0) {
    std::cerr << reference_name << " reported solver diagnostics, status=" << status << '\n';
    return 1;
  }
  const Result ii{gas[0], gas[1], particle[0], particle[1], sulfate, water, hydronium};
  std::cout << std::scientific << std::setprecision(7) << std::left << std::setw(22) << "Field"
            << std::setw(12) << "Unit" << std::right << std::setw(17) << "MKPP Lite"
            << std::setw(17) << reference_name << std::setw(17) << "Lite minus II" << '\n';
  report("NH3", "mol m-3", lite.nh3, ii.nh3);
  report("HNO3", "mol m-3", lite.hno3, ii.hno3);
  report("ANH4", "mol m-3", lite.nh4, ii.nh4);
  report("ANO3", "mol m-3", lite.no3, ii.no3);
  report("ASO4", "mol m-3", lite.sulfate, ii.sulfate);
  report("AH3OP", "mol m-3", lite.hydronium, ii.hydronium);
  report("aerosol water", "kg m-3", lite.water, ii.water);
  const bool nhx_ok = check_particle_fraction("NHx particle fraction", lite.nh4, ii.nh4, totals[2]);
  const bool nitrate_ok = check_particle_fraction("NO3 particle fraction", lite.no3, ii.no3, totals[3]);
  const bool lite_ok = physical(lite, "MKPP Lite");
  const bool ii_ok = physical(ii, reference_name);
  return lite_ok && ii_ok && nhx_ok && nitrate_ok ? 0 : 1;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc == 2 && std::string_view(argv[1]) == "cold") return run(278.15);
  if (argc == 2 && std::string_view(argv[1]) == "warm") return run(298.15);
  std::cerr << "Expected cold or warm\n";
  return 2;
}
