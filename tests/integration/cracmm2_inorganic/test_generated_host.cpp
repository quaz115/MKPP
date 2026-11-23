#include "mkpp.hpp"
#include "mkpp_host/cracmm2_inorganic_state.hpp"

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string_view>
#include <vector>

namespace {

TEST(Cracmm2GeneratedHost, AppliesDarkKineticsThenPartitionsOneOfTwoCells) {
    constexpr std::size_t cells = 2;
    mkpp::MechanismContext context(cells);
    context.initialize();
    const mkpp::host::Cracmm2InorganicBinding binding(context);
    const auto species = context.getSpeciesCount();
    ASSERT_EQ(species, static_cast<std::size_t>(MKPP_NUM_SPECIES));
    auto offset = [&](std::string_view name, std::size_t cell) {
        return cell + cells * static_cast<std::size_t>(context.getSpeciesIndex(name));
    };

    std::array<double, cells> temp{298.15, 291.89}, rh{0.65, 0.83};
    std::array<double, cells> pressure{101325., 101325.}, density{2.46e19, 2.51e19};
    std::vector<double> state(cells * species, 0.);
    for (std::size_t cell = 0; cell < cells; ++cell) {
        state[offset("AIR", cell)] = density[cell];
        state[offset("M", cell)] = density[cell];
        state[offset("O2", cell)] = 0.2095 * density[cell];
        state[offset("N2", cell)] = 0.7808 * density[cell];
        state[offset("H2O", cell)] = 0.015 * density[cell];
        state[offset("NO", cell)] = 1.e10 * (cell + 1);
        state[offset("O3", cell)] = 7.e11;
    }
    // EPFL's first forward organic-aerosol case. Initialize the volatile
    // inorganic components in the gas to exercise transfer on the second cell.
    constexpr double umol_m3 = 6.02214076e11;
    constexpr std::size_t target = 1;
    state[offset("NH3", target)] = 0.37 / 17.0 * umol_m3;
    state[offset("HNO3", target)] = 0.087 / 63.0 * umol_m3;
    state[offset("HCL", target)] = 0.87 / 36.5 * umol_m3;
    state[offset("ASO4", target)] = 1.088 / 98.0 * umol_m3;
    state[offset("ANA", target)] = 0.54 / 23.0 * umol_m3;
    state[offset("ACA", target)] = 0.01 / 40.1 * umol_m3;
    state[offset("AK", target)] = 0.006 / 39.1 * umol_m3;
    state[offset("AMG", target)] = 0.004 / 24.3 * umol_m3;

    std::vector<double> dark_forcing(cells * MKPP_NUM_RATE_PARAMETERS, 0.);
    context.setStatePointers(state.data(), temp.data(), pressure.data(), density.data());
    context.setRateParameterPointers(dark_forcing.data(), MKPP_NUM_RATE_PARAMETERS);
    context.setAerosolPointers(rh.data(), nullptr);
    // A short positive step makes NO2 via the dark NO + O3 reaction (R036).
    context.integrate(0.1);
    ASSERT_GT(state[offset("NO2", 0)], 0.);
    ASSERT_GT(state[offset("NO2", target)], 0.);
    const auto after_kinetics = state;

    const mkpp::host::Cracmm2ThermoEnvironment env{temp[target], rh[target], 1.508e-9, 0.15, 1000.};
    ASSERT_EQ(binding.partition_cell(state, cells, target, env), mkpp::host::Cracmm2ThermoStatus::success);
    for (std::size_t i = 0; i < species; ++i) {
        EXPECT_EQ(state[cells * i], after_kinetics[cells * i]) << "other cell, species " << i;
    }
    for (auto name : {"ASO4", "ANA", "ACA", "AK", "AMG", "SULF", "O3", "NO2"}) {
        EXPECT_EQ(state[offset(name, target)], after_kinetics[offset(name, target)]) << name;
    }
    auto check_total = [&](std::string_view gas, std::string_view ion) {
        const auto g = offset(gas, target);
        const auto p = offset(ion, target);
        EXPECT_GE(state[g], 0.);
        EXPECT_GE(state[p], 0.);
        EXPECT_DOUBLE_EQ(state[g] + state[p], after_kinetics[g] + after_kinetics[p]);
    };
    check_total("NH3", "ANH4");
    check_total("HNO3", "ANO3");
    check_total("HCL", "ACL");
    EXPECT_GT(state[offset("ANH4", target)], after_kinetics[offset("ANH4", target)]);
    EXPECT_GT(state[offset("AH2O", target)], 0.);
    EXPECT_GT(state[offset("AORGH2O", target)], 0.);
    EXPECT_GT(state[offset("AH3OP", target)], 0.);
}

TEST(Cracmm2GeneratedHost, UpdatesTwoFineModesAndRetainsExcludedReservoirs) {
    using namespace mkpp::host;
    constexpr std::size_t cells = 2, target = 1;
    constexpr double umol_m3 = 6.02214076e11;
    constexpr double water_per_kg_m3 = 6.02214076e23 / (0.018 * 1.e6);
    mkpp::MechanismContext context(cells);
    context.initialize();
    const Cracmm2InorganicBinding binding(context);
    std::vector<double> state(cells * context.getSpeciesCount(), 0.);
    auto offset = [&](std::string_view name, std::size_t cell) {
        return cell + cells * static_cast<std::size_t>(context.getSpeciesIndex(name));
    };
    auto value = [&](std::string_view name) -> double& { return state[offset(name, target)]; };
    std::array<double, cells> temperature{298.15, 291.89}, pressure{101325., 101325.}, density{2.46e19, 2.51e19};
    context.setStatePointers(state.data(), temperature.data(), pressure.data(), density.data());
    for (std::size_t i = 0; i < context.getSpeciesCount(); ++i) state[cells * i] = 10. + i;

    // The first EPFL mixed-ion fixture supplies gas and fine-mode totals.
    // The generated particle fields include an additional excluded reservoir.
    // Give that reservoir a different composition and separate water amounts.
    value("NH3") = 0.37 / 17. * umol_m3;
    value("HNO3") = 0.087 / 63. * umol_m3;
    value("HCL") = 0.87 / 36.5 * umol_m3;
    value("O3") = 7.e11;
    value("SULF") = 5.e8;
    using Mode = Cracmm2InorganicMode;
    std::array<Mode, 2> modes{};
    modes[0].aso4 = 1.088 / 98. * umol_m3;
    modes[1].ana = 0.54 / 23. * umol_m3;
    modes[1].aca = 0.01 / 40.1 * umol_m3;
    modes[1].ak = 0.006 / 39.1 * umol_m3;
    modes[1].amg = 0.004 / 24.3 * umol_m3;
    constexpr std::array<std::string_view, 8> particle{"ANH4","ANO3","ACL","ASO4","ANA","ACA","AK","AMG"};
    constexpr std::array<double Mode::*, 8> fields{&Mode::anh4, &Mode::ano3, &Mode::acl, &Mode::aso4,
                                                    &Mode::ana, &Mode::aca, &Mode::ak, &Mode::amg};
    constexpr std::array<double, 8> outside{1.e10, 2.e10, 3.e10, 4.e10, 5.e10, 6.e10, 7.e10, 8.e10};
    for (std::size_t i = 0; i < fields.size(); ++i) {
        value(particle[i]) = outside[i] + (modes[0].*fields[i] + modes[1].*fields[i]);
    }
    Cracmm2ThermoDiagnostics fine{1.4e-9, 0.4e-9, 3.e5};
    constexpr double outside_inorganic_water = 0.7e-9 * water_per_kg_m3;
    constexpr double outside_organic_water = 0.2e-9 * water_per_kg_m3;
    constexpr double outside_hydronium = 9.e6;
    value("AH2O") = outside_inorganic_water + (fine.total_kg_m3 - fine.organic_kg_m3) * water_per_kg_m3;
    value("AORGH2O") = outside_organic_water + fine.organic_kg_m3 * water_per_kg_m3;
    value("AH3OP") = outside_hydronium + fine.hydronium_cm3;
    const auto original = state;
    const std::array<double, 2> sinks{1., 3.};
    const Cracmm2ThermoEnvironment env{temperature[target], 0.83, 1.508e-9, 0.15, 1000.};

    // Repeat from the returned state to detect double counting on the next
    // timestep. The gas and fine-mode arrays must describe the same cell.
    for (int step = 0; step < 2; ++step) {
        ASSERT_EQ(binding.partition_fine_modes_cell(state, cells, target, modes, sinks, env, fine),
                  Cracmm2ThermoStatus::success);
        for (std::size_t i = 0; i < particle.size(); ++i) {
            EXPECT_DOUBLE_EQ(value(particle[i]), outside[i] + (modes[0].*fields[i] + modes[1].*fields[i])) << particle[i];
        }
        for (std::size_t i = 0; i < context.getSpeciesCount(); ++i) EXPECT_EQ(state[cells * i], original[cells * i]);
        for (const auto name : {"ASO4","ANA","ACA","AK","AMG","O3","SULF"}) {
            EXPECT_EQ(value(name), original[offset(name, target)]) << name;
        }
        constexpr std::array<std::string_view, 3> gases{"NH3","HNO3","HCL"};
        for (std::size_t i = 0; i < gases.size(); ++i) {
            EXPECT_DOUBLE_EQ(value(gases[i]) + value(particle[i]),
                             original[offset(gases[i], target)] + original[offset(particle[i], target)]) << gases[i];
        }
        EXPECT_DOUBLE_EQ(value("AH2O"), outside_inorganic_water +
                         (fine.total_kg_m3 - fine.organic_kg_m3) * water_per_kg_m3);
        EXPECT_DOUBLE_EQ(value("AORGH2O"), outside_organic_water + fine.organic_kg_m3 * water_per_kg_m3);
        EXPECT_DOUBLE_EQ(value("AH3OP"), outside_hydronium + fine.hydronium_cm3);
        EXPECT_GT(modes[0].ano3, 0.);
        EXPECT_GT(modes[1].ano3, 0.);
    }

    // Declare less total than the supplied fine-only contribution. Both
    // nonvolatile material and stale prior water diagnostics must be rejected.
    for (const auto name : {"ASO4", "AORGH2O"}) {
        const double saved = value(name);
        value(name) = 0.;
        const auto before = state;
        const auto modes_before = modes;
        const auto fine_before = fine;
        EXPECT_EQ(binding.partition_fine_modes_cell(state, cells, target, modes, sinks, env, fine),
                  Cracmm2ThermoStatus::invalid_input);
        EXPECT_EQ(state, before);
        for (std::size_t i = 0; i < modes.size(); ++i) {
            for (auto field : fields) EXPECT_EQ(modes[i].*field, modes_before[i].*field);
        }
        EXPECT_EQ(fine.total_kg_m3, fine_before.total_kg_m3);
        EXPECT_EQ(fine.organic_kg_m3, fine_before.organic_kg_m3);
        EXPECT_EQ(fine.hydronium_cm3, fine_before.hydronium_cm3);
        value(name) = saved;
    }
}

} // namespace

int main(int argc, char** argv) {
    Kokkos::ScopeGuard scope(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
