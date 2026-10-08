# PEXSI solver fixes and findings (ABACUS side)

These files record an investigation of ABACUS's PEXSI (`ks_solver pexsi`) support: why its
results and timings diverge from the dense solvers (`genelpa`, `scalapack_gvx`), a minimal
patch that fixes the parts on the ABACUS side, and the verification data.

They are **not** part of the aiida-abacus Python package — they are kept here so the ABACUS
patch and the evidence travel together with the plugin that drives those calculations.

## Files

| File | Content |
|---|---|
| `pexsi-energy-accounting.patch` | 3-file ABACUS patch (+19/−2) against **v3.10.1**: two compile fixes for the `#ifdef __PEXSI` code paths plus the energy-accounting fix |
| `VERIFICATION_REPORT.md` | What the patch changes and how it was verified (rank consistency, energy decomposition, `pexsi_method` sweep) |
| `SOLVER_COMPARISON.md` | Earlier local benchmark: `pexsi` vs `genelpa` vs `scalapack_gvx` on diamond carbon (accuracy + wall time), including the build recipe for a PEXSI-enabled ABACUS |

## The patch in one paragraph

`HSolverLCAO::solve` assigned PEXSI's **free energy** (`totalFreeEnergy`) to `f_en.eband` and
never set `f_en.demet` (−TS), and it did so only on the ranks that drive the pole expansion, so
every other MPI rank kept `eband = 0` and reported a different total energy (measured: rank 0
`−279.6757 eV`, ranks 1–3 `−358.9418 eV` — the difference is exactly one `E_band`). The patch
publishes the three PEXSI scalars on the whole communicator, assigns the true band energy
(`totalEnergyH = Tr[H·DM]`) to `eband`, and puts the remainder (`totalFreeEnergy − totalEnergyH`)
into `demet`. After the patch all ranks agree, and `E_band` matches the dense solver to 0.4 meV,
but the **total energy is unchanged** (the split is algebraically invariant).

## What actually fixes the accuracy (verified)

The accuracy switch is PEXSI's pole-expansion method, not the accounting:

| setting (diamond C, Γ-only, 2 atoms) | total energy | vs dense reference |
|---|---|---|
| `pexsi_method 1` + `npole 40` (**ABACUS default**) | −279.6757230584546 eV | **+333.3 meV** |
| `pexsi_method 2` + `npole 40` | −280.0089654571346 eV | **+0.013 meV** |

This reproduces, on a different PEXSI/ABACUS stack (PEXSI 1.2.0 + ABACUS 3.10.1), the pole-table
accuracy limits measured on the yeesuan cluster with PEXSI 2.0.0 — so `pexsi_method 2` should be
treated as required, not optional, whenever PEXSI results are compared with `genelpa`.

## Applying the patch

```bash
cd <abacus-develop checkout>          # tested against tag v3.10.1
git apply patches/abacus-pexsi/pexsi-energy-accounting.patch
```

Build with `-DENABLE_PEXSI=ON` (see `SOLVER_COMPARISON.md` §2 for a conda-forge recipe:
`pexsi`, `superlu_dist`, `parmetis`, `metis`, `scalapack`, `fftw`, `cereal`).

## Status

* Reported here as evidence for review; not yet proposed upstream.
* The two compile fixes show ABACUS's PEXSI glue regressed between v3.9.x (where it builds) and
  v3.10.x/develop (where `elecstate_lcao.cpp` and `edm.cpp` fail under `#ifdef __PEXSI`).
