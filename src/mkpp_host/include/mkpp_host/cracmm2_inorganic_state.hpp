#pragma once

#include "mkpp_host/cracmm2_inorganic.hpp"

#include <array>
#include <cmath>
#include <cstddef>
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
        if (cell >= cells || concentrations.size() % cells != 0 ||
            concentrations.size() / cells != species_count_) {
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

        Cracmm2AerosolWater water{};
        const auto status = partition_cracmm2_inorganic(state, env, water);
        if (status != Cracmm2ThermoStatus::success) {
            return status;
        }
        // AH2O and AORGH2O are separate CRACMM2 tracers, each with declared
        // molar mass 18 g/mol. Lite reports total water including organic water.
        constexpr double water_per_kg_m3 = 6.02214076e23 / (0.018 * 1.0e6);
        const double inorganic_water = (water.total_kg_m3 - water.organic_kg_m3) * water_per_kg_m3;
        const double organic_water = water.organic_kg_m3 * water_per_kg_m3;
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
        return Cracmm2ThermoStatus::success;
    }

private:
    enum Slot { nh3, anh4, hno3, ano3, hcl, acl, aso4, ana, aca, ak, amg, ah2o, aorgh2o, count };
    static constexpr std::array<const char*, count> names_ = {
        "NH3", "ANH4", "HNO3", "ANO3", "HCL", "ACL", "ASO4", "ANA", "ACA", "AK", "AMG", "AH2O", "AORGH2O"};
    std::array<std::size_t, count> indices_{};
    std::size_t species_count_;
};

} // namespace mkpp::host
