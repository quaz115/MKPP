#pragma once

#include "mkpp_host/cracmm2_inorganic_fine_modes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>

namespace mkpp::host {

// Construct once using the generated CRACMM2 mkpp::MechanismContext. The
// lookup uses the indices of the live kinetics build, including any reordering.
class Cracmm2InorganicBinding {
public:
    template <class MechanismContext>
    explicit Cracmm2InorganicBinding(const MechanismContext& context)
        : species_count_(context.getSpeciesCount()) {
        for (std::size_t i = 0; i < names_.size(); ++i) {
            const int index = context.getSpeciesIndex(names_[i]);
            if (index < 0 || static_cast<std::size_t>(index) >= species_count_) {
                throw std::invalid_argument("CRACMM2 inorganic species index outside mechanism");
            }
            indices_[i] = static_cast<std::size_t>(index);
            for (std::size_t j = 0; j < i; ++j) {
                if (indices_[j] == indices_[i]) {
                    throw std::invalid_argument("CRACMM2 inorganic species indices must be distinct");
                }
            }
        }
    }

    // Use the same contiguous host memory passed to MechanismContext::setStatePointers.
    // LayoutLeft indexes (cell,species) as cell + cells * species. Concentrations
    // are molecules or ion formula units cm^-3. Call after a successful kinetics step.
    [[nodiscard]] Cracmm2ThermoStatus partition_cell(std::span<double> concentrations,
                                                     std::size_t cells, std::size_t cell,
                                                     const Cracmm2ThermoEnvironment& env) const {
        if (!valid_cell(concentrations, cells, cell)) {
            return Cracmm2ThermoStatus::invalid_input;
        }
        auto value = [&](Slot slot) -> double& {
            return concentrations[cell + cells * indices_[slot]];
        };
        Cracmm2InorganicState state{};
        state.nh3 = value(nh3);
        state.anh4 = value(anh4);
        state.hno3 = value(hno3);
        state.ano3 = value(ano3);
        state.hcl = value(hcl);
        state.acl = value(acl);
        state.aso4 = value(aso4);
        state.ana = value(ana);
        state.aca = value(aca);
        state.ak = value(ak);
        state.amg = value(amg);

        Cracmm2ThermoDiagnostics diagnostics{};
        const auto status = partition_cracmm2_inorganic(state, env, diagnostics);
        if (status != Cracmm2ThermoStatus::success) {
            return status;
        }
        // AH2O and AORGH2O are separate CRACMM2 tracers, each with declared
        // molar mass 18 g/mol. Lite reports total water including organic water.
        const double inorganic_water = (diagnostics.total_kg_m3 - diagnostics.organic_kg_m3) * water_per_kg_m3;
        const double organic_water = diagnostics.organic_kg_m3 * water_per_kg_m3;
        if (!std::isfinite(inorganic_water) || !std::isfinite(organic_water) ||
            inorganic_water < 0.0 || organic_water < 0.0) {
            return Cracmm2ThermoStatus::invalid_result;
        }

        value(nh3) = state.nh3;
        value(anh4) = state.anh4;
        value(hno3) = state.hno3;
        value(ano3) = state.ano3;
        value(hcl) = state.hcl;
        value(acl) = state.acl;
        value(ah2o) = inorganic_water;
        value(aorgh2o) = organic_water;
        value(ah3op) = diagnostics.hydronium_cm3;
        return Cracmm2ThermoStatus::success;
    }

    // The generated particle slots contain all modes. Pass the two fine-mode
    // inventories, their condensation sinks and their previous diagnostics.
    // On success, the three gases and the selected fine contribution to each
    // bulk tracer are updated once. Excluded reservoirs are conserved. The
    // caller retains the returned fine diagnostics across timesteps/restarts.
    [[nodiscard]] Cracmm2ThermoStatus partition_fine_modes_cell(
        std::span<double> concentrations, std::size_t cells, std::size_t cell,
        std::span<Cracmm2InorganicMode, 2> modes, std::span<const double, 2> condensation_sinks,
        const Cracmm2ThermoEnvironment& env, Cracmm2ThermoDiagnostics& fine_diagnostics) const {
        if (!valid_cell(concentrations, cells, cell)) {
            return Cracmm2ThermoStatus::invalid_input;
        }
        auto value = [&](Slot slot) -> double& {
            return concentrations[cell + cells * indices_[slot]];
        };

        // The first three families exchange with the gas phase. The remaining
        // five are checked because they enter the fine equilibrium calculation.
        constexpr std::array<Slot, 8> aerosol_slots{anh4, ano3, acl, aso4, ana, aca, ak, amg};
        using Mode = Cracmm2InorganicMode;
        constexpr std::array<double Mode::*, 8> members{
            &Mode::anh4, &Mode::ano3, &Mode::acl, &Mode::aso4, &Mode::ana, &Mode::aca, &Mode::ak, &Mode::amg};
        constexpr std::array<Slot, 3> diagnostic_slots{ah2o, aorgh2o, ah3op};

        std::array<double, 8> excluded_aerosol{};
        for (std::size_t i = 0; i < aerosol_slots.size(); ++i) {
            const double selected = modes[0].*members[i] + modes[1].*members[i];
            if (!excluded_amount(value(aerosol_slots[i]), selected, excluded_aerosol[i])) {
                return Cracmm2ThermoStatus::invalid_input;
            }
        }
        const auto prior = diagnostic_amounts(fine_diagnostics);
        std::array<double, 3> excluded_diagnostics{};
        for (std::size_t i = 0; i < diagnostic_slots.size(); ++i) {
            if (!excluded_amount(value(diagnostic_slots[i]), prior[i], excluded_diagnostics[i])) {
                return Cracmm2ThermoStatus::invalid_input;
            }
        }

        // Stage all outputs: a late overflow when combining the fine and
        // excluded amounts must leave both the modal and generated state intact.
        Cracmm2InorganicGas gas{value(nh3), value(hno3), value(hcl)};
        std::array<Mode, 2> updated_modes{modes[0], modes[1]};
        Cracmm2ThermoDiagnostics next{};
        const auto status = partition_cracmm2_fine_modes(gas, updated_modes, condensation_sinks, env, next);
        if (status != Cracmm2ThermoStatus::success) {
            return status;
        }
        std::array<double, 3> aerosol{}, diagnostics{};
        const auto next_amounts = diagnostic_amounts(next);
        for (std::size_t i = 0; i < aerosol.size(); ++i) {
            aerosol[i] = excluded_aerosol[i] + (updated_modes[0].*members[i] + updated_modes[1].*members[i]);
            diagnostics[i] = excluded_diagnostics[i] + next_amounts[i];
            if (!std::isfinite(aerosol[i]) || !std::isfinite(diagnostics[i])) {
                return Cracmm2ThermoStatus::invalid_result;
            }
        }

        value(nh3) = gas.nh3;
        value(hno3) = gas.hno3;
        value(hcl) = gas.hcl;
        for (std::size_t i = 0; i < aerosol.size(); ++i) {
            value(aerosol_slots[i]) = aerosol[i];
            value(diagnostic_slots[i]) = diagnostics[i];
        }
        std::copy(updated_modes.begin(), updated_modes.end(), modes.begin());
        fine_diagnostics = next;
        return Cracmm2ThermoStatus::success;
    }

private:
    static constexpr double water_per_kg_m3 = 6.02214076e23 / (0.018 * 1.0e6);

    bool valid_cell(std::span<const double> state, std::size_t cells, std::size_t cell) const {
        return cell < cells && state.size() % cells == 0 && state.size() / cells == species_count_;
    }

    static std::array<double, 3> diagnostic_amounts(const Cracmm2ThermoDiagnostics& values) {
        return {(values.total_kg_m3 - values.organic_kg_m3) * water_per_kg_m3,
                values.organic_kg_m3 * water_per_kg_m3, values.hydronium_cm3};
    }

    static bool excluded_amount(double total, double selected, double& excluded) {
        if (!std::isfinite(total) || !std::isfinite(selected) || total < 0.0 || selected < 0.0) {
            return false;
        }
        // All arguments are formula units per cm3 of air. Allow roundoff
        // from summation/unit conversion, scaled to each tracer's own amount.
        constexpr double roundoff = 16 * std::numeric_limits<double>::epsilon();
        const double difference = total - selected;
        if (difference < -roundoff * std::max(total, selected)) {
            return false;
        }
        excluded = std::max(0.0, difference);
        return true;
    }

    enum Slot { nh3, anh4, hno3, ano3, hcl, acl, aso4, ana, aca, ak, amg, ah2o, aorgh2o, ah3op, count };
    static constexpr std::array<const char*, count> names_ = {
        "NH3", "ANH4", "HNO3", "ANO3", "HCL", "ACL", "ASO4", "ANA", "ACA", "AK", "AMG", "AH2O", "AORGH2O", "AH3OP"};
    std::array<std::size_t, count> indices_{};
    std::size_t species_count_;
};

} // namespace mkpp::host
