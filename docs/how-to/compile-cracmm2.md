---
type: how-to
category: mechanisms
tags: [cracmm, epa, openatmos, compilation]
---

# Compile CRACMM2

MKPP includes the Community Regional Atmospheric Chemistry Multiphase Mechanism
(CRACMM2) published by the U.S. Environmental Protection Agency. The canonical
input is [mechanism.json](../../mechanisms/openatmos/cracmm2/mechanism.json),
generated from the pinned EPA reaction and species tables in the adjacent
`source` directory.

The source snapshot is pinned to commit
`8a71fdbfb2183a1e37a53eb237e265d12f09b245` from
[USEPA/CRACMM](https://github.com/USEPA/CRACMM). Regenerate the canonical input
after updating those source files with:

```bash
python3 scripts/import_cracmm2.py \
  --species mechanisms/openatmos/cracmm2/source/cracmm2_metadata.csv \
  --reactions mechanisms/openatmos/cracmm2/source/cracmm2_rxn_metadata.csv \
  --output mechanisms/openatmos/cracmm2/mechanism.json
```

CRACMM2 uses MKPP's default kinetic units (`molecule cm-3`), so no unit
declaration or rate-coefficient conversion is needed. Supply solver
concentrations and air number density on this same basis.

Compile the generated MKPP solver with the normal catalog command:

```bash
python3 -m mkpp.cli compile \
  mechanisms/openatmos/cracmm2/mechanism.json \
  --test-env example_env.yaml \
  --out mkpp-generated \
  --adjoint --emit-manifest --report --no-cache
```

Recognized multi-channel reaction rates are mapped to MKPP's `EP2` and `EP3`
forms; all 11 recognized pressure-dependent associations are mapped to the
standard `TROE` reaction type. The six remaining `UNKNOWN` rates are
expressions that cannot be represented safely and are listed in
`metadata.unconverted_rates` with their original expression and reason. The
current unresolved reactions are R051, R064, R067, R167, R169, and HAL_Ozone.
They need a separately implemented custom-rate path; the importer does not
substitute guessed constants.

## Test bulk inorganic equilibrium on the host

This optional step uses the C++ ISORROPIA-Lite revision pinned by MKPP. From the
MKPP checkout, with a C++23 compiler, CMake, and Kokkos available, run:

```sh
git submodule update --init src/isorropialite
cmake -S . -B build-lite -DBUILD_TESTING=ON -DMKPP_ENABLE_CRACMM2_ISORROPIA_LITE=ON -DMKPP_TEST_LITE_FORTRAN=ON
cmake --build build-lite --target cracmm2_inorganic_checks --parallel 2
ctest --test-dir build-lite -L cracmm2_inorganic --output-on-failure --no-tests=error
```

The reference option requires a Fortran compiler and builds the original Lite
routines from the pinned submodule. Four forward box cases compare the three gas
amounts, their particulate families, sulfate, particulate hydronium and organic and total water. The condensed
reference sums the corresponding aqueous and solid species independently.
The cases cover a wet ammonium–nitrate–sulfate mixture, a wet mixture with
sodium, chloride and crustal ions, colder inorganic-only aerosol and sulfate-rich aerosol. The comparison uses 1% of the reference
plus 5 pmol m⁻³ for each chemical family or 0.1 ng m⁻³ for water. These
tolerances apply to these fixtures only.

To exercise both kinetics and partitioning through the generated host API, use
an MKPP Python environment, emit a matching bundle, and repeat the build:

```sh
python3 -m mkpp.cli compile mechanisms/openatmos/cracmm2/mechanism.json --test-env example_env.yaml --host-api --no-cache --out build-host-api
cmake -S . -B build-lite -DMKPP_CRACMM2_HOST_API_DIR="$PWD/build-host-api"
cmake --build build-lite --target cracmm2_inorganic_checks e2e_cracmm2_runner --parallel 2
ctest --test-dir build-lite -L cracmm2_inorganic --output-on-failure
ctest --test-dir build-lite -R '^E2ESolver_cracmm2$' --output-on-failure
```

For a cross-model comparison, supply `-DMKPP_CMAQ_ISORROPIA_II_DIR=<cmaq>/CCTM/src/aero/aero6`
at configure time, build `cracmm2_inorganic_checks`, and run
`ctest --test-dir build-lite -R '^cracmm2_cmaq_ii_' --verbose`. CI uses the
[CMAQ 5.5.0.4 reference sources](https://github.com/USEPA/CMAQ/tree/1b70e0a395fc42cc34a4355807d39df1d3dbae84/CCTM/src/aero/aero6).
Compare forward metastable inorganic-only cases at the same temperature, RH
and component totals. Cross-model differences quantify sensitivity to the
thermodynamics implementation. The separate legacy Lite tests supply the
same-model numerical acceptance criteria.

CMAQ applies forward equilibrium to fine modes. In hybrid mode it also uses
a reverse equilibrium calculation at the coarse-particle surface to drive
time-dependent mass transfer. Validation of that coarse-mode path requires a size-resolved transfer test. See the [CMAQ fine-mode calculation](https://github.com/USEPA/CMAQ/blob/1b70e0a395fc42cc34a4355807d39df1d3dbae84/CCTM/src/aero/aero6/aero_subs.F#L1465-L1506)
and [coarse-mode calculation](https://github.com/USEPA/CMAQ/blob/1b70e0a395fc42cc34a4355807d39df1d3dbae84/CCTM/src/aero/aero6/aero_subs.F#L1198-L1211).

### Call from a CRACMM2 host

Link your host to `MKPP::cracmm2_inorganic` and include
`mkpp_host/cracmm2_inorganic_state.hpp`. Construct the binding once with the
generated `mkpp::MechanismContext`. After a successful kinetic step on the
host, pass the same writable concentration buffer registered with the
context. Its storage is contiguous, with cell index varying fastest.

```cpp
mkpp::host::Cracmm2InorganicBinding binding(context);
// At each thermodynamics step, after synchronizing state to host memory:
const auto status = binding.partition_cell(concentrations, n_cells, cell, environment);
if (status != mkpp::host::Cracmm2ThermoStatus::success) {
    // The cell is unchanged. Propagate the failure to the model driver.
    throw std::runtime_error("CRACMM2 inorganic equilibrium failed");
}
```

Use molecules or ionic formula units cm⁻³ for the concentration buffer,
temperature in K and RH as a fraction. The adapter converts the eight
component totals to mol m⁻³ of air before the direct solver call. Optional
organic forcing comprises dry organic mass in kg m⁻³, hygroscopicity and
material density in kg m⁻³. It writes back NH3/ANH4, HNO3/ANO3 and HCL/ACL.
AH2O receives inorganic water and AORGH2O the organic share. AH3OP receives
one hydronium ion per reported aqueous H+ ion. The direct solver reports the
H+ amount per air volume. A pH diagnostic also requires an explicit water basis and activity convention. Select one
process to own these water tracers. Provide the host organic forcing whenever
organic-associated water is intended. The typed single-cell interface
`partition_cracmm2_inorganic` can instead return water diagnostics to a
driver that owns those tracers.

Sulfate supplied to the adapter is particulate ASO4. The host must handle
SULF condensation and nucleation before supplying the particulate pool.
The bulk entry point returns cell totals. The fine-mode entry point below
also distributes the change in particulate ammonium, nitrate, and chloride.
The driving model schedules deposition, SOA chemistry and heterogeneous uptake. It must be scheduled once for each
chosen equilibrium step. Compiling the kinetic mechanism alone does not
schedule this process. Tests against legacy Lite and tests against a tagged
CMAQ configuration with ISORROPIA II answer different scientific questions.

### Partition gas against selected fine modes

For a size-resolved host, call `partition_fine_modes_cell` in place of the
bulk call for the same pool and timestep. The generated-host overload takes
the cell's Aitken and accumulation reservoirs. Supply their ion amounts in
formula units cm⁻³, the host's condensation sink for each mode, and organic
forcing for **only these modes**. Keep the gas concentrations in the generated
state. Each aggregate inorganic tracer in the generated state must include
the selected fine-mode amounts.

```cpp
const auto status = binding.partition_fine_modes_cell(
    concentrations, n_cells, cell, fine_modes, condensation_sinks, environment,
    fine_diagnostics);
if (status != mkpp::host::Cracmm2ThermoStatus::success) {
    throw std::runtime_error("CRACMM2 fine-mode equilibrium failed");
}
```

`fine_diagnostics` is per-cell state. Initialize it from the selected modes'
previous total water, organic-associated water (both kg m⁻³ of air), and
hydronium (ions cm⁻³ of air). Use zeros for an initially dry fine aerosol.
Pass its updated value on the next step and restore it with the modal
reservoirs on restart. This lets the binding replace only the fine contribution
to the all-mode AH2O, AORGH2O and AH3OP tracers. A host that owns its own
aggregate fields can call `partition_cracmm2_fine_modes` directly with a span
of selected modes and consume the returned fine-total diagnostics.

One forward Lite solve supplies the equilibrium gas amounts and particulate
totals. Each mode receives a share of **newly condensed** material proportional
to its supplied sink. The sinks use common units and at least one must be
positive when material condenses. For evaporation, each mode loses a share
proportional to its original amount of that particular ion. The concentrations
of sulfate, sodium, calcium, potassium and magnesium in each mode are held
fixed through this step. The coarse reservoir has its own exchange step.

The two allocation rules follow the
[CMAQ 5.5.0.4 fine-mode transfer calculation](https://github.com/USEPA/CMAQ/blob/1b70e0a395fc42cc34a4355807d39df1d3dbae84/CCTM/src/aero/aero6/aero_subs.F#L1542-L1581).
CMAQ computes its sinks from wet modal growth factors and also applies a
trace-chloride guard and evaporation floor. Here the aerosol driver supplies
the sinks. Water and hydronium are diagnosed for the selected fine pool as
a whole. Applications needing mode-specific water or acidity must calculate
them from modal composition. See the
[CMAQ modal water](https://github.com/USEPA/CMAQ/blob/1b70e0a395fc42cc34a4355807d39df1d3dbae84/CCTM/src/aero/aero6/aero_subs.F#L1608-L1625)
and [hydronium calculations](https://github.com/USEPA/CMAQ/blob/1b70e0a395fc42cc34a4355807d39df1d3dbae84/CCTM/src/aero/aero6/aero_subs.F#L1658-L1672).
