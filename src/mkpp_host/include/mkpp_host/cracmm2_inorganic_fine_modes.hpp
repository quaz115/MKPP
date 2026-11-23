#pragma once

#include "mkpp_host/cracmm2_inorganic.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace mkpp::host {

struct Cracmm2InorganicGas {
    double nh3{}, hno3{}, hcl{};
};

// Each field is the number of molecular or ionic formula units per cm3 of air.
struct Cracmm2InorganicMode {
    double anh4{}, ano3{}, acl{}, aso4{}, ana{}, aca{}, ak{}, amg{};
};

// Equilibrate shared gas against just the fine modes supplied by the host.
// The host owns the condensation sinks (common units) and any transfer of SULF.
// Water and hydronium are diagnostics for the entire selected fine-mode pool.
// Gas, modal amounts and diagnostics are committed together on success.
inline Cracmm2ThermoStatus partition_cracmm2_fine_modes(
    Cracmm2InorganicGas& gas, std::span<Cracmm2InorganicMode> fine_modes,
    std::span<const double> condensation_sinks, const Cracmm2ThermoEnvironment& env,
    Cracmm2ThermoDiagnostics& diagnostics) {
    if (fine_modes.empty() || fine_modes.size() != condensation_sinks.size()) {
        return Cracmm2ThermoStatus::invalid_input;
    }
    using Mode = Cracmm2InorganicMode;
    constexpr std::array<double Mode::*, 8> components{
        &Mode::anh4, &Mode::ano3, &Mode::acl, &Mode::aso4, &Mode::ana, &Mode::aca, &Mode::ak, &Mode::amg};
    std::array<long double, components.size()> totals{};
    std::size_t largest_sink = 0;
    for (std::size_t m = 0; m < fine_modes.size(); ++m) {
        const double sink = condensation_sinks[m];
        if (!std::isfinite(sink) || sink < 0.0) return Cracmm2ThermoStatus::invalid_input;
        if (sink > condensation_sinks[largest_sink]) largest_sink = m;
        for (std::size_t c = 0; c < components.size(); ++c) {
            const double value = fine_modes[m].*components[c];
            if (!std::isfinite(value) || value < 0.0) return Cracmm2ThermoStatus::invalid_input;
            totals[c] += value;
        }
    }
    for (const auto amount : totals) {
        if (amount > std::numeric_limits<double>::max()) return Cracmm2ThermoStatus::invalid_input;
    }

    Cracmm2InorganicState bulk{
        .nh3 = gas.nh3, .anh4 = static_cast<double>(totals[0]),
        .hno3 = gas.hno3, .ano3 = static_cast<double>(totals[1]),
        .hcl = gas.hcl, .acl = static_cast<double>(totals[2]),
        .aso4 = static_cast<double>(totals[3]), .ana = static_cast<double>(totals[4]),
        .aca = static_cast<double>(totals[5]), .ak = static_cast<double>(totals[6]),
        .amg = static_cast<double>(totals[7])};
    Cracmm2ThermoDiagnostics new_diagnostics{};
    const auto status = partition_cracmm2_inorganic(bulk, env, new_diagnostics);
    if (status != Cracmm2ThermoStatus::success) return status;

    const long double max_sink = condensation_sinks[largest_sink];
    long double scaled_sink_sum = 0.0L;
    if (max_sink > 0.0L) {
        for (const double sink : condensation_sinks) scaled_sink_sum += sink / max_sink;
    }

    // Keep the common Aitken/accumulation case allocation-free. Staging also
    // guarantees rollback if another family cannot be allocated.
    using Volatiles = std::array<double, 3>;
    std::array<Volatiles, 2> local{};
    std::vector<Volatiles> extra;
    std::span<Volatiles> updated;
    if (fine_modes.size() <= local.size()) {
        updated = std::span(local).first(fine_modes.size());
    } else {
        extra.resize(fine_modes.size());
        updated = extra;
    }
    const Volatiles old_gas{gas.nh3, gas.hno3, gas.hcl};
    Volatiles new_gas{bulk.nh3, bulk.hno3, bulk.hcl};
    const Volatiles target{bulk.anh4, bulk.ano3, bulk.acl};
    for (std::size_t c = 0; c < target.size(); ++c) {
        const auto member = components[c];
        const long double before = totals[c];
        const long double after = target[c];
        const bool condensing = after > before;
        bool retain_modes = target[c] == static_cast<double>(before);

        if (!retain_modes && condensing && max_sink == 0.0L) {
            // The scalar adapter obtains its particle amount by subtracting
            // the final gas from the family total. Ignore only changes within
            // subtraction roundoff when no condensation sink is available.
            const long double roundoff = 8.0L * std::numeric_limits<double>::epsilon() *
                                         std::max(before, static_cast<long double>(old_gas[c]));
            if (after - before > roundoff) return Cracmm2ThermoStatus::unsupported_input;
            new_gas[c] = old_gas[c];
            retain_modes = true;
        }
        if (retain_modes) {
            for (std::size_t m = 0; m < fine_modes.size(); ++m) updated[m][c] = fine_modes[m].*member;
            continue;
        }

        // Condensation distributes the gain by sink. Evaporation distributes
        // the surviving material by the initial inventory of THIS species.
        // Allocating the surviving amount avoids cancellation near exhaustion.
        const long double amount = condensing ? after - before : after;
        std::size_t pivot = largest_sink;
        if (!condensing) {
            for (std::size_t m = 0; m < fine_modes.size(); ++m) {
                if (fine_modes[m].*member > fine_modes[pivot].*member) pivot = m;
            }
        }
        long double assigned = 0.0L;
        for (std::size_t m = 0; m < fine_modes.size(); ++m) {
            if (m == pivot) continue;
            const long double old = fine_modes[m].*member;
            const long double fraction = condensing ? (condensation_sinks[m] / max_sink) / scaled_sink_sum
                                                    : old / before;
            const long double share = amount * fraction;
            updated[m][c] = static_cast<double>((condensing ? old : 0.0L) + share);
            assigned += share;
        }
        // The largest participating mode takes the residual so tiny modes can
        // retain their representable concentrations even in an asymmetric split.
        const long double remainder = amount - assigned;
        if (remainder < 0.0L) return Cracmm2ThermoStatus::invalid_result;
        updated[pivot][c] = static_cast<double>((condensing ? fine_modes[pivot].*member : 0.0L) + remainder);
    }

    for (const auto& mode : updated) {
        for (double amount : mode) {
            if (!std::isfinite(amount) || amount < 0.0) return Cracmm2ThermoStatus::invalid_result;
        }
    }
    gas = {new_gas[0], new_gas[1], new_gas[2]};
    for (std::size_t m = 0; m < fine_modes.size(); ++m) {
        fine_modes[m].anh4 = updated[m][0];
        fine_modes[m].ano3 = updated[m][1];
        fine_modes[m].acl = updated[m][2];
    }
    diagnostics = new_diagnostics;
    return Cracmm2ThermoStatus::success;
}

} // namespace mkpp::host
