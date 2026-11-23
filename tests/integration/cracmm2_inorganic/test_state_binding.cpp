#include "mkpp_host/cracmm2_inorganic_state.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <string_view>
#include <vector>

namespace {

using namespace mkpp::host;

// Test double for the generated MechanismContext name lookup. Deliberately
// change species order and interleave unrelated species to detect positional
// assumptions independently of the order used by a particular kinetics build.
struct PermutedContext {
    std::array<std::string_view, 17> names = {"ASO4", "HNO3", "O3", "ANA", "ACL", "AH2O", "AH3OP",
        "NH3", "ACA", "AORGH2O", "AMG", "SULF", "ANO3", "AK", "ANH4", "HCL", "OH"};

    std::size_t getSpeciesCount() const { return names.size(); }

    int getSpeciesIndex(std::string_view name) const {
        const auto it = std::find(names.begin(), names.end(), name);
        return it == names.end() ? -1 : static_cast<int>(it - names.begin());
    }
};

TEST(Cracmm2InorganicBinding, EquilibratesOneCellOfLayoutLeftHostState) {
    const PermutedContext context;
    const Cracmm2InorganicBinding binding(context);
    constexpr std::size_t cells = 2;
    std::vector<double> state(cells * context.getSpeciesCount());
    std::iota(state.begin(), state.end(), 1.0);
    auto offset = [&](std::string_view name) { return cells * context.getSpeciesIndex(name) + 1; };

    // Forward with organics, first row of EPFL's Partitioning_with_organics.INP.
    // The totals in the file are acid/base equivalents in micrograms/m3.
    constexpr double umol_m3 = 6.02214076e11;
    Cracmm2InorganicState cell{};
    cell.nh3 = 0.37 / 17.0 * umol_m3;
    cell.hno3 = 0.087 / 63.0 * umol_m3;
    cell.hcl = 0.87 / 36.5 * umol_m3;
    cell.aso4 = 1.088 / 98.0 * umol_m3;
    cell.ana = 0.54 / 23.0 * umol_m3;
    cell.aca = 0.01 / 40.1 * umol_m3;
    cell.ak = 0.006 / 39.1 * umol_m3;
    cell.amg = 0.004 / 24.3 * umol_m3;

    state[offset("NH3")] = cell.nh3;
    state[offset("ANH4")] = cell.anh4;
    state[offset("HNO3")] = cell.hno3;
    state[offset("ANO3")] = cell.ano3;
    state[offset("HCL")] = cell.hcl;
    state[offset("ACL")] = cell.acl;
    state[offset("ASO4")] = cell.aso4;
    state[offset("ANA")] = cell.ana;
    state[offset("ACA")] = cell.aca;
    state[offset("AK")] = cell.ak;
    state[offset("AMG")] = cell.amg;
    auto expected = state;
    Cracmm2ThermoDiagnostics diagnostics{};
    const Cracmm2ThermoEnvironment env{291.89, 0.83, 1.508e-9, 0.15, 1000.0};
    ASSERT_EQ(partition_cracmm2_inorganic(cell, env, diagnostics), Cracmm2ThermoStatus::success);
    expected[offset("NH3")] = cell.nh3;
    expected[offset("ANH4")] = cell.anh4;
    expected[offset("HNO3")] = cell.hno3;
    expected[offset("ANO3")] = cell.ano3;
    expected[offset("HCL")] = cell.hcl;
    expected[offset("ACL")] = cell.acl;
    constexpr double molecules_per_kg_m3 = 6.02214076e23 / (0.018 * 1.0e6);
    expected[offset("AH2O")] = (diagnostics.total_kg_m3 - diagnostics.organic_kg_m3) * molecules_per_kg_m3;
    expected[offset("AORGH2O")] = diagnostics.organic_kg_m3 * molecules_per_kg_m3;
    expected[offset("AH3OP")] = diagnostics.hydronium_cm3;
    ASSERT_EQ(binding.partition_cell(state, cells, 1, env), Cracmm2ThermoStatus::success);
    EXPECT_EQ(state, expected); // Includes first cell, cations, SULF and unrelated species.

    EXPECT_EQ(binding.partition_cell(state, cells, 0, {298.15, 83.0}), Cracmm2ThermoStatus::invalid_input);
    EXPECT_EQ(binding.partition_cell(state, cells, 2, env), Cracmm2ThermoStatus::invalid_input);
    EXPECT_EQ(binding.partition_cell(state, 0, 0, env), Cracmm2ThermoStatus::invalid_input);
    EXPECT_EQ(binding.partition_cell(std::span(state).subspan(1), cells, 0, env), Cracmm2ThermoStatus::invalid_input);
    EXPECT_EQ(state, expected);
}

TEST(Cracmm2InorganicBinding, RejectsMissingSpeciesAtInitialization) {
    PermutedContext context;
    context.names[context.getSpeciesIndex("ASO4")] = "ASO4_ALT";
    EXPECT_THROW((void)Cracmm2InorganicBinding(context), std::invalid_argument);
}

} // namespace
