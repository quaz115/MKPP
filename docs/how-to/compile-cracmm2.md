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
routines from the pinned submodule. Two forward box cases compare the three gas
amounts, their particulate families, and organic and total water. The condensed
reference sums the corresponding aqueous and solid species independently.
The cases cover a wet ammonium–nitrate–sulfate mixture and a wet mixture with
sodium, chloride and crustal ions. The comparison uses 1% of the reference
plus 5 pmol m⁻³ for each chemical family or 0.1 ng m⁻³ for water. These
tolerances apply to these fixtures only.

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
AH2O receives the inorganic water and AORGH2O the organic share. Select one
process to own these water tracers. Provide the host organic forcing whenever
organic-associated water is intended. The typed single-cell interface
`partition_cracmm2_inorganic` can instead return water diagnostics to a
driver that owns those tracers.

Sulfate supplied to the adapter is particulate ASO4. The host must handle
SULF condensation and nucleation before supplying the particulate pool.
The adapter operates on bulk amounts and leaves modal redistribution, dry
deposition, SOA chemistry, heterogeneous rate updates and the choice of
chemistry timestep to the driving model. It must be scheduled once for each
chosen equilibrium step. Compiling the kinetic mechanism alone does not
schedule this process. Tests against legacy Lite and tests against a tagged
CMAQ configuration with ISORROPIA II answer different scientific questions.
