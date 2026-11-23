#include "mkpp_host/cracmm2_inorganic_fine_modes.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <span>

namespace {

using namespace mkpp::host;

constexpr double per_umol_m3 = 6.02214076e11; // formula units/cm3 per micromole/m3 air
constexpr double sulfate = per_umol_m3 / 98.0;
constexpr double ammonia = 2.0 * per_umol_m3 / 17.0;
constexpr double nitrate = per_umol_m3 / 63.0;
constexpr Cracmm2ThermoEnvironment wet{298.15, 0.80, 1.0e-8, 0.15, 1000.0};

void expect_concentration(double actual, double expected) {
    EXPECT_NEAR(actual, expected, 1.0e-7 + 1.0e-11 * std::abs(expected));
}

auto ions(const Cracmm2InorganicMode& m) {
    return std::array{m.anh4, m.ano3, m.acl, m.aso4, m.ana, m.aca, m.ak, m.amg};
}

auto nonvolatile(const Cracmm2InorganicMode& m) {
    return std::array{m.aso4, m.ana, m.aca, m.ak, m.amg};
}

void expect_diagnostics(const Cracmm2ThermoDiagnostics& actual, const Cracmm2ThermoDiagnostics& expected) {
    EXPECT_NEAR(actual.total_kg_m3, expected.total_kg_m3, 1.e-19 + 1.e-10 * expected.total_kg_m3);
    EXPECT_NEAR(actual.organic_kg_m3, expected.organic_kg_m3, 1.e-19 + 1.e-10 * expected.organic_kg_m3);
    expect_concentration(actual.hydronium_cm3, expected.hydronium_cm3);
}

TEST(Cracmm2FineModes, CondensesOnlyIntoSelectedModesWithPositiveSink) {
    // First legacy Lite test1.inp record. The receiving second mode has seed
    // surface supplied by its host and initially has no inorganic material.
    Cracmm2InorganicGas gas{ammonia, nitrate, 0.0};
    std::array<Cracmm2InorganicMode, 4> modes{};
    modes[0].aso4 = sulfate;
    modes[3] = {123.0, 456.0, 789.0, 101.0, 112.0, 131.0, 415.0, 161.0};
    const auto coarse_before = ions(modes[3]);
    const double large = std::numeric_limits<double>::max();
    const std::array sinks{large, large / 2, 0.0}; // finite 2:1 ratio, raw sum would overflow
    Cracmm2InorganicState bulk{};
    bulk.nh3 = ammonia;
    bulk.hno3 = nitrate;
    bulk.aso4 = sulfate;
    Cracmm2ThermoDiagnostics reference{}, result{};
    ASSERT_EQ(partition_cracmm2_inorganic(bulk, wet, reference), Cracmm2ThermoStatus::success);

    ASSERT_EQ(partition_cracmm2_fine_modes(gas, std::span(modes).first<3>(), sinks, wet, result),
              Cracmm2ThermoStatus::success);
    ASSERT_GT(bulk.anh4, 0.0);
    ASSERT_GT(bulk.ano3, 0.0);
    for (std::size_t i = 0; i < 3; ++i) {
        const double share = std::array{2./3, 1./3, 0.0}[i];
        expect_concentration(modes[i].anh4, share * bulk.anh4);
        expect_concentration(modes[i].ano3, share * bulk.ano3);
        EXPECT_EQ(modes[i].acl, 0.0);
    }
    EXPECT_EQ(modes[0].aso4, sulfate);
    EXPECT_EQ(modes[1].aso4, 0.0);
    EXPECT_EQ(modes[2].aso4, 0.0);
    EXPECT_EQ(ions(modes[3]), coarse_before);
    expect_concentration(gas.nh3, bulk.nh3);
    expect_concentration(gas.hno3, bulk.hno3);
    // Original Lite Fortran test1.inp, record 1: condensed NHx and nitrate
    // in micromoles/m3. Use the independent reference harness's 1% + 5 pmol/m3.
    EXPECT_NEAR((modes[0].anh4 + modes[1].anh4) / per_umol_m3, 0.02380612, 0.01 * 0.02380612 + 5.e-6);
    EXPECT_NEAR((modes[0].ano3 + modes[1].ano3) / per_umol_m3, 0.003484534, 0.01 * 0.003484534 + 5.e-6);
    expect_concentration(gas.nh3 + modes[0].anh4 + modes[1].anh4 + modes[2].anh4, ammonia);
    expect_concentration(gas.hno3 + modes[0].ano3 + modes[1].ano3 + modes[2].ano3, nitrate);
    expect_diagnostics(result, reference);
}

TEST(Cracmm2FineModes, MixedFluxesUseSpeciesSpecificEvaporationAndHostCondensationWeights) {
    // First Partitioning_with_organics.INP record. The same component totals
    // start with all NHx and Cl condensed and all nitrate in the gas.
    constexpr Cracmm2ThermoEnvironment environment{291.89, 0.83, 1.508e-9, 0.15, 1000.0};
    Cracmm2InorganicState bulk{};
    bulk.anh4 = 0.37 / 17.0 * per_umol_m3;
    bulk.hno3 = 0.087 / 63.0 * per_umol_m3;
    bulk.acl = 0.87 / 36.5 * per_umol_m3;
    bulk.aso4 = 1.088 / 98.0 * per_umol_m3;
    bulk.ana = 0.54 / 23.0 * per_umol_m3;
    bulk.aca = 0.01 / 40.1 * per_umol_m3;
    bulk.ak = 0.006 / 39.1 * per_umol_m3;
    bulk.amg = 0.004 / 24.3 * per_umol_m3;
    const auto initial = bulk;
    Cracmm2InorganicGas gas{0.0, bulk.hno3, 0.0};
    std::array<Cracmm2InorganicMode, 2> modes{};
    constexpr double first_ammonium_share = 1.e-9;
    modes[0].anh4 = first_ammonium_share * bulk.anh4;
    modes[1].anh4 = bulk.anh4 - modes[0].anh4;
    modes[0].acl = 0.75 * bulk.acl;
    modes[1].acl = bulk.acl - modes[0].acl;
    modes[0].aso4 = bulk.aso4;
    modes[1].ana = bulk.ana;
    modes[1].aca = bulk.aca;
    modes[0].ak = bulk.ak;
    modes[0].amg = bulk.amg;
    const auto before = modes;
    Cracmm2ThermoDiagnostics result{}, reference{};
    ASSERT_EQ(partition_cracmm2_inorganic(bulk, environment, reference), Cracmm2ThermoStatus::success);
    ASSERT_LT(bulk.anh4, initial.anh4);
    ASSERT_LT(bulk.acl, initial.acl);
    ASSERT_GT(bulk.ano3, 0.0);
    ASSERT_EQ(partition_cracmm2_fine_modes(gas, modes, std::array{1.0, 3.0}, environment, result),
              Cracmm2ThermoStatus::success);

    expect_concentration(modes[0].anh4, first_ammonium_share * bulk.anh4);
    expect_concentration(modes[1].anh4, (1.0 - first_ammonium_share) * bulk.anh4);
    expect_concentration(modes[0].acl, 0.75 * bulk.acl);
    expect_concentration(modes[1].acl, 0.25 * bulk.acl);
    expect_concentration(modes[0].ano3, 0.25 * bulk.ano3);
    expect_concentration(modes[1].ano3, 0.75 * bulk.ano3);
    for (std::size_t i = 0; i < modes.size(); ++i) {
        EXPECT_EQ(nonvolatile(modes[i]), nonvolatile(before[i]));
    }
    expect_concentration(gas.nh3 + modes[0].anh4 + modes[1].anh4, initial.anh4);
    expect_concentration(gas.hcl + modes[0].acl + modes[1].acl, initial.acl);
    expect_concentration(gas.hno3 + modes[0].ano3 + modes[1].ano3, initial.hno3);
    expect_diagnostics(result, reference);

    // Reversing storage together with its host-supplied sinks preserves each mode.
    std::array reversed{before[1], before[0]};
    Cracmm2InorganicGas other_gas{0.0, initial.hno3, 0.0};
    Cracmm2ThermoDiagnostics other{};
    ASSERT_EQ(partition_cracmm2_fine_modes(other_gas, reversed, std::array{3.0, 1.0}, environment, other),
              Cracmm2ThermoStatus::success);
    for (std::size_t i = 0; i < modes.size(); ++i) {
        expect_concentration(reversed[1-i].anh4, modes[i].anh4);
        expect_concentration(reversed[1-i].ano3, modes[i].ano3);
        expect_concentration(reversed[1-i].acl, modes[i].acl);
    }
}

TEST(Cracmm2FineModes, SingleEvaporatingReservoirNeedsNoCondensationSink) {
    Cracmm2InorganicGas gas{};
    std::array<Cracmm2InorganicMode, 1> modes{};
    modes[0].anh4 = ammonia;
    modes[0].ano3 = nitrate;
    modes[0].aso4 = sulfate;
    Cracmm2ThermoDiagnostics result{};
    ASSERT_EQ(partition_cracmm2_fine_modes(gas, modes, std::array{0.0}, wet, result),
              Cracmm2ThermoStatus::success);
    EXPECT_GT(gas.nh3, 0.0);
    EXPECT_GT(gas.hno3, 0.0);
    EXPECT_GT(modes[0].anh4, 0.0);
    EXPECT_GT(modes[0].ano3, 0.0);
    expect_concentration(gas.nh3 + modes[0].anh4, ammonia);
    expect_concentration(gas.hno3 + modes[0].ano3, nitrate);
}

TEST(Cracmm2FineModes, InvalidInputMissingSinkAndSolverDiagnosticLeaveOutputsUntouched) {
    Cracmm2InorganicGas gas{ammonia, nitrate, 0.0};
    std::array<Cracmm2InorganicMode, 2> modes{};
    modes[0].aso4 = sulfate;
    Cracmm2ThermoDiagnostics result{1.0, 2.0, 3.0};
    auto check_failure = [&](std::span<const double> sinks, Cracmm2ThermoStatus expected) {
        const auto before = modes;
        EXPECT_EQ(partition_cracmm2_fine_modes(gas, modes, sinks, wet, result), expected);
        EXPECT_EQ(gas.nh3, ammonia);
        EXPECT_EQ(gas.hno3, nitrate);
        EXPECT_EQ(gas.hcl, 0.0);
        EXPECT_EQ(ions(modes[0]), ions(before[0]));
        EXPECT_EQ(ions(modes[1]), ions(before[1]));
        EXPECT_EQ(result.total_kg_m3, 1.0);
        EXPECT_EQ(result.organic_kg_m3, 2.0);
        EXPECT_EQ(result.hydronium_cm3, 3.0);
    };
    check_failure(std::array{0.0, 0.0}, Cracmm2ThermoStatus::unsupported_input);
    check_failure(std::array{1.0}, Cracmm2ThermoStatus::invalid_input);
    check_failure(std::array{1.0, -1.0}, Cracmm2ThermoStatus::invalid_input);
    check_failure(std::array{1.0, std::numeric_limits<double>::infinity()}, Cracmm2ThermoStatus::invalid_input);
    modes[1].ano3 = -1.0;
    check_failure(std::array{1.0, 1.0}, Cracmm2ThermoStatus::invalid_input);
    modes[1].ano3 = 0.0;
    modes[1].aca = per_umol_m3; // excess Ca triggers Lite forward diagnostic 51
    check_failure(std::array{1.0, 1.0}, Cracmm2ThermoStatus::solver_diagnostic);
    modes[1].aca = 0.0;
    EXPECT_EQ(partition_cracmm2_fine_modes(gas, modes, std::array{1.0, 1.0}, wet, result),
              Cracmm2ThermoStatus::success);
}

} // namespace
