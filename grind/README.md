# GRIND METHOD — qgraf-vs-Julia diagnostic infrastructure

Setup used to find the root cause of BUG 1 (qgen flavor-loop under-count,
HANDOFF.md Session 22). Workflow:

  1. Copy `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08` → `qgraf-grind.f08`.
  2. Add `WRITE(0, ...)` instrumentation at qgen entry/match/emit points
     (qgen:13889, 13900, 14455) and at the end of qrvi to dump dpntro/rotvpo.
  3. Build with debug symbols: `gfortran -o qgraf-grind -cpp -P -O0 -g -J fmodules qgraf-grind.f08`.
  4. Run `./qgraf-grind ctrl.dat 2> grind_trace.txt` for ground truth.
  5. Mirror instrumentation in `dump_julia_emissions.jl` for the same case.
  6. Diff topologies/emissions to find the first divergence.

## Files committed (small, our own work)

| File | Purpose |
|------|---------|
| `ctrl.dat` | qgraf control file for ee→μμ 1L (model=qed2, no opts) |
| `dump_julia_emissions.jl` | Julia diagnostic — dumps every (topo, ps1, pmap) emission |
| `inspect_dpntro.gdb` | gdb script to inspect rotvpo/nrot after qrvi |
| `inspect_qgen.gdb` | gdb script to break at qgen vfo/match/emit |
| `run_v2_tests.sh` | Sequential v2 test runner with incremental output |

## Files NOT committed (gitignored)

- `qgraf-grind`, `qgraf-grind.f08` — derived from qgraf source (externally licensed)
- `fmodules/`, `models/`, `styles/`, `out_tmp/` — qgraf build artefacts and copies
- `*_trace.txt`, `*_stdout.txt` — qgraf-output diagnostic logs (regenerable)
- `v2_test_results.txt` — ad-hoc test log (regenerable)

To rebuild, copy the source from refs/, instrument, compile.
