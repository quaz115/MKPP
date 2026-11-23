#include "mkpp_host/cracmm2_inorganic.hpp"

#include <gtest/gtest.h>

#include <limits>
#include <utility>

namespace {

using namespace mkpp::host;

// Conversion from micromoles/m3 air to molecular formula units/cm3 air.
constexpr double per_umol_m3 = 6.02214076e11;

Cracmm2InorganicState wet_nitrate_cell(double ammonia_ug_m3 = 2.0) {
    Cracmm2InorganicState cell{};
    // Input mass columns in the legacy fixture are H2SO4, NH3 and HNO3
    // equivalents. Their molecular weights differ from NH4+, NO3-, and SO4--.
    cell.aso4 = (1.0 / 98.0) * per_umol_m3;
    cell.nh3 = (ammonia_ug_m3 / 17.0) * per_umol_m3;
    cell.hno3 = (1.0 / 63.0) * per_umol_m3;
    return cell;
}

constexpr Cracmm2ThermoEnvironment wet_conditions{298.15, 0.8, 1.0e-8, 0.15, 1000.0};

TEST(Cracmm2Inorganic, WetReferenceReturnsGasParticleAndWater) {
    auto cell = wet_nitrate_cell();
    const auto initial = cell;
    Cracmm2ThermoDiagnostics diagnostics{};
    ASSERT_EQ(partition_cracmm2_inorganic(cell, wet_conditions, diagnostics), Cracmm2ThermoStatus::success);

    // Pinned Lite tests/test_c_api.cpp, Run 1 of isolite1_0_src/test1.inp.
    // The gas tolerances here are expressed in micromoles/m3, matching that fixture.
    EXPECT_NEAR(cell.nh3 / per_umol_m3, 0.093891, 1.0e-4);
    EXPECT_NEAR(cell.hno3 / per_umol_m3, 0.012439, 1.0e-4);
    EXPECT_NEAR(diagnostics.total_kg_m3, 8.083366e-9, 1.0e-11);
    // AERLIQ(1) = 2.170486e-11 mol H+ / m3 air in the pinned C API fixture.
    EXPECT_NEAR(diagnostics.hydronium_cm3 / per_umol_m3, 2.170486e-5, 1.0e-7);
    // kappa-Koehler uptake at RH=0.8 with organic/water density ratio of 1.
    EXPECT_NEAR(diagnostics.organic_kg_m3, 6.0e-9, 1.0e-14);
    EXPECT_GT(diagnostics.total_kg_m3, diagnostics.organic_kg_m3);
    EXPECT_GT(cell.anh4, 0.0);
    EXPECT_GT(cell.ano3, 0.0);
    EXPECT_NEAR(cell.nh3 + cell.anh4, initial.nh3, 1.0e-12 * initial.nh3);
    EXPECT_NEAR(cell.hno3 + cell.ano3, initial.hno3, 1.0e-12 * initial.hno3);
    EXPECT_DOUBLE_EQ(cell.aso4, initial.aso4);
    EXPECT_DOUBLE_EQ(cell.hcl, 0.0);
    EXPECT_DOUBLE_EQ(cell.acl, 0.0);
}

TEST(Cracmm2Inorganic, InitialPhasesAndCellOrderDoNotDetermineEquilibrium) {
    auto gas_start = wet_nitrate_cell();
    auto particle_start = gas_start;
    std::swap(particle_start.nh3, particle_start.anh4);
    std::swap(particle_start.hno3, particle_start.ano3);
    Cracmm2ThermoDiagnostics gas_result{}, particle_result{}, other_result{};
    ASSERT_EQ(partition_cracmm2_inorganic(gas_start, wet_conditions, gas_result), Cracmm2ThermoStatus::success);

    // A different cell between the two identical-total cells exercises reuse.
    auto other = wet_nitrate_cell(6.0);
    ASSERT_EQ(partition_cracmm2_inorganic(other, wet_conditions, other_result), Cracmm2ThermoStatus::success);
    EXPECT_GT(other.nh3, gas_start.nh3);
    ASSERT_EQ(partition_cracmm2_inorganic(particle_start, wet_conditions, particle_result),
              Cracmm2ThermoStatus::success);
    EXPECT_DOUBLE_EQ(particle_start.nh3, gas_start.nh3);
    EXPECT_DOUBLE_EQ(particle_start.hno3, gas_start.hno3);
    EXPECT_DOUBLE_EQ(particle_start.anh4, gas_start.anh4);
    EXPECT_DOUBLE_EQ(particle_start.ano3, gas_start.ano3);
    EXPECT_DOUBLE_EQ(particle_result.total_kg_m3, gas_result.total_kg_m3);
    EXPECT_DOUBLE_EQ(particle_result.hydronium_cm3, gas_result.hydronium_cm3);
}

TEST(Cracmm2Inorganic, EmptyInorganicCellPreservesExactZeros) {
    Cracmm2InorganicState empty{};
    Cracmm2ThermoDiagnostics diagnostics{1.0, 2.0, 3.0};
    ASSERT_EQ(partition_cracmm2_inorganic(empty, {280.0, 0.5}, diagnostics), Cracmm2ThermoStatus::success);
    EXPECT_DOUBLE_EQ(empty.nh3 + empty.anh4 + empty.hno3 + empty.ano3 + empty.hcl + empty.acl, 0.0);
    EXPECT_DOUBLE_EQ(diagnostics.total_kg_m3, 0.0);
    EXPECT_DOUBLE_EQ(diagnostics.organic_kg_m3, 0.0);
    EXPECT_DOUBLE_EQ(diagnostics.hydronium_cm3, 0.0);
    // Organic-only water is handled by the host when forcing is nonzero.
    EXPECT_EQ(partition_cracmm2_inorganic(empty, wet_conditions, diagnostics), Cracmm2ThermoStatus::unsupported_input);
}

TEST(Cracmm2Inorganic, InvalidStateOrMeteorologyLeavesOutputsUntouched) {
    auto cell = wet_nitrate_cell();
    const double old_nh3 = cell.nh3;
    Cracmm2ThermoDiagnostics diagnostics{2.0e-8, 1.0e-8, 7.0e6};
    EXPECT_EQ(partition_cracmm2_inorganic(cell, {298.15, 80.0}, diagnostics), Cracmm2ThermoStatus::invalid_input);
    EXPECT_EQ(partition_cracmm2_inorganic(cell, {0.0, 0.8}, diagnostics), Cracmm2ThermoStatus::invalid_input);
    cell.ano3 = -1.0;
    EXPECT_EQ(partition_cracmm2_inorganic(cell, wet_conditions, diagnostics), Cracmm2ThermoStatus::invalid_input);
    EXPECT_DOUBLE_EQ(cell.ano3, -1.0);
    cell.ano3 = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(partition_cracmm2_inorganic(cell, wet_conditions, diagnostics), Cracmm2ThermoStatus::invalid_input);
    EXPECT_DOUBLE_EQ(cell.nh3, old_nh3);
    EXPECT_DOUBLE_EQ(diagnostics.total_kg_m3, 2.0e-8);
    EXPECT_DOUBLE_EQ(diagnostics.organic_kg_m3, 1.0e-8);
    EXPECT_DOUBLE_EQ(diagnostics.hydronium_cm3, 7.0e6);
}

TEST(Cracmm2Inorganic, SolverCationAdjustmentDoesNotModifyHostState) {
    auto cell = wet_nitrate_cell();
    // Lite forward ISRP4F adjusts Ca if it exceeds all available anions and
    // records diagnostic 51. The bulk adapter has no reservoir for that loss.
    cell.aca = 1.0 * per_umol_m3;
    const auto initial = cell;
    Cracmm2ThermoDiagnostics diagnostics{1.0e-9, 0.0, 1.0e6};
    EXPECT_EQ(partition_cracmm2_inorganic(cell, wet_conditions, diagnostics), Cracmm2ThermoStatus::solver_diagnostic);
    EXPECT_DOUBLE_EQ(cell.aca, initial.aca);
    EXPECT_DOUBLE_EQ(cell.nh3, initial.nh3);
    EXPECT_DOUBLE_EQ(cell.anh4, initial.anh4);
    EXPECT_DOUBLE_EQ(diagnostics.total_kg_m3, 1.0e-9);
    EXPECT_DOUBLE_EQ(diagnostics.hydronium_cm3, 1.0e6);
}

} // namespace
