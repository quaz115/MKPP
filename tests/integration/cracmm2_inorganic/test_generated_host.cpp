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

} // namespace

int main(int argc, char** argv) {
    Kokkos::ScopeGuard scope(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
