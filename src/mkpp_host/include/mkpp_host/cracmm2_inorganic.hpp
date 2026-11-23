#pragma once

#include <Isorropia/Isorropia.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace mkpp::host {

// Bulk CRACMM2 molecular or ionic formula units per cm3 of air.
// ASO4 is particulate sulfate. The host owns gaseous SULF and its condensation.
struct Cracmm2InorganicState {
    double nh3{}, anh4{}, hno3{}, ano3{}, hcl{}, acl{};
    double aso4{}, ana{}, aca{}, ak{}, amg{};
};

struct Cracmm2ThermoEnvironment {
    double temperature_k{};
    double relative_humidity{}; // fraction
    double organic_mass_kg_m3{};
    double organic_kappa{};
    double organic_density_kg_m3{1000.0};
};

struct Cracmm2ThermoDiagnostics {
    // Organic-associated water is included in total aerosol water.
    double total_kg_m3{};
    double organic_kg_m3{};
    double hydronium_cm3{}; // aqueous H3O+ ions per cm3 of air (CRACMM2 AH3OP)
};

enum class Cracmm2ThermoStatus { success, invalid_input, unsupported_input, solver_diagnostic, invalid_result };

// Forward equilibrium on one host cell. On failure neither output is modified.
// If a second process advances inorganic partitioning, the caller must select
// a single owner. This call does not advance kinetics, allocate modes, or add
// condensed H2SO4, organic nitrates, or water to the host chemical state.
inline Cracmm2ThermoStatus partition_cracmm2_inorganic(Cracmm2InorganicState& state,
                                                      const Cracmm2ThermoEnvironment& env,
                                                      Cracmm2ThermoDiagnostics& diagnostics) {
    constexpr double to_mol_m3 = 1.0e6 / 6.02214076e23;
    auto nonnegative = [](double value) { return std::isfinite(value) && value >= 0.0; };

    const std::array<double, 11> species = {state.nh3,  state.anh4, state.hno3, state.ano3,
                                            state.hcl,  state.acl,  state.aso4, state.ana,
                                            state.aca,  state.ak,   state.amg};
    if (!std::all_of(species.begin(), species.end(), nonnegative) ||
        !std::isfinite(env.temperature_k) || env.temperature_k <= 0.0 ||
        env.temperature_k >= std::numeric_limits<int>::max() || // Lite casts T to int
        !nonnegative(env.relative_humidity) || env.relative_humidity > 1.0 ||
        !nonnegative(env.organic_mass_kg_m3) || !nonnegative(env.organic_kappa) ||
        !nonnegative(env.organic_density_kg_m3) ||
        (env.organic_mass_kg_m3 > 0.0 && env.organic_density_kg_m3 == 0.0)) {
        return Cracmm2ThermoStatus::invalid_input;
    }

    if (std::all_of(species.begin(), species.end(), [](double value) { return value == 0.0; })) {
        // The inorganic solver floors absent ions. Preserve an empty dry cell
        // exactly; organic-only aerosol water is outside this adapter's scope.
        if (env.organic_mass_kg_m3 > 0.0 && env.organic_kappa > 0.0) {
            return Cracmm2ThermoStatus::unsupported_input;
        }
        diagnostics = {};
        return Cracmm2ThermoStatus::success;
    }

    const std::array<double, 8> totals = {state.ana, state.aso4, state.nh3 + state.anh4,
                                          state.hno3 + state.ano3, state.hcl + state.acl,
                                          state.aca, state.ak, state.amg};
    if (!std::all_of(totals.begin(), totals.end(), nonnegative)) {
        return Cracmm2ThermoStatus::invalid_input;
    }

    // The pinned Lite CLI converts its input file to mol/m3 BEFORE the C++
    // call. The C ABI forwards w[] unchanged and expects the same molar basis.
    IsorropiaInput input{};
    for (std::size_t i = 0; i < totals.size(); ++i) {
        input.w[i] = totals[i] * to_mol_m3;
    }
    input.iprob = 0; // forward, total-component problem
    input.temp = env.temperature_k;
    input.rh = env.relative_humidity;
    input.org[0] = env.organic_mass_kg_m3;
    input.org[1] = env.organic_kappa;
    input.org[2] = env.organic_density_kg_m3;

    IsorropiaState solved{};
    isorropia_solve_c(&input, &solved);
    if (solved.num_errors != 0) {
        // The flat ABI supplies an error count without severity or case codes.
        return Cracmm2ThermoStatus::solver_diagnostic;
    }

    if (!nonnegative(solved.water) || !nonnegative(solved.watcmp[23]) ||
        solved.watcmp[23] > solved.water) {
        return Cracmm2ThermoStatus::invalid_result;
    }

    // Lite molal[1] corresponds to legacy AERLIQ(1), moles of aqueous H+
    // per m3 of AIR. Each ion contributes one AH3OP formula unit. A dry cell
    // or one at Lite's liquid-water floor has no resolved aqueous hydronium.
    constexpr double water_floor_kg_m3 = 1.0e-20;
    const bool dry = solved.dryf != 0 || solved.water <= water_floor_kg_m3;
    const double hydronium_cm3 = dry ? 0.0 : solved.molal[1] / to_mol_m3;
    if (!nonnegative(hydronium_cm3)) {
        return Cracmm2ThermoStatus::invalid_result;
    }

    const std::array<double, 3> gas_mol_m3 = {solved.gnh3, solved.ghno3, solved.ghcl};
    std::array<double, 3> gas_host{};
    for (std::size_t i = 0; i < gas_mol_m3.size(); ++i) {
        const double total_mol_m3 = input.w[2 + i];
        // Lite imposes an internal floor of 1e-20 mol/m3 on absent components.
        // Correct only that floor and floating-point roundoff at the endpoints.
        const double tolerance = 1.0e-20 + 64 * std::numeric_limits<double>::epsilon() * total_mol_m3;
        const double gas = gas_mol_m3[i];
        if (!std::isfinite(gas) || gas < -tolerance || gas > total_mol_m3 + tolerance) {
            return Cracmm2ThermoStatus::invalid_result;
        }
        gas_host[i] = std::min(totals[2 + i], std::max(0.0, gas) / to_mol_m3);
    }

    // Forward waer is an input copy. Bulk aerosol formula-unit totals follow
    // from the original conserved component minus its equilibrium gas.
    state.nh3 = gas_host[0];
    state.anh4 = totals[2] - state.nh3;
    state.hno3 = gas_host[1];
    state.ano3 = totals[3] - state.hno3;
    state.hcl = gas_host[2];
    state.acl = totals[4] - state.hcl;
    diagnostics = {solved.water, solved.watcmp[23], hydronium_cm3};
    return Cracmm2ThermoStatus::success;
}

} // namespace mkpp::host
