# `refs/qgraf/` — qgraf-4.0.6 reference

qgraf is Paulo Nogueira's Feynman diagram generator. We use it as the
**porting oracle** for `src/v2/qgraf/` (Strategy C port). The Fortran source,
PDF manuals, and distribution tarballs are © Paulo Nogueira and are
explicitly marked "not to be redistributed without explicit permission" —
these files are **gitignored**.

## Committed in this directory

- `ALGORITHM.md` — our cleanroom algorithm spec written from static source
  analysis. This is **our work product**, not derived from the qgraf
  documentation or PDFs. See the header of the file for the analysis
  methodology.
- `v4.0.6/qgraf-4.0.6.dir/generate_golden_masters.py` — our script that
  runs qgraf on a battery of (model, process, loops, options) configurations
  and writes `golden_masters/**/*.out` + `golden_masters/manifest.json` +
  `golden_masters/SUMMARY.md`.
- `v4.0.6/qgraf-4.0.6.dir/parse_golden_master.py` — our script that parses
  qgraf's `array.sty`-formatted `.out` files into structured JSON.
- `v4.0.6/qgraf-4.0.6.dir/golden_masters/` — 104 generated test fixtures:
  per-case `.out` files, `manifest.json` (5.5 MB), `SUMMARY.md`. These are
  **generated data**, used directly by `scripts/qgraf_golden_master_report.jl`
  and `test/v2/test_diagram_gen.jl`.

## Acquiring qgraf source (not committed)

The canonical citation is:
> P. Nogueira, *Automatic Feynman graph generation*, J. Comput. Phys.
> 105 (1993) 279–289, doi:10.1006/jcph.1993.1074.

qgraf is distributed from Paulo Nogueira's CeFEMA page at IST:

```bash
# Check latest pointer:
#   https://cfif.ist.utl.pt/~paulo/qgraf.html
# (mirrors may vary — contact Paulo Nogueira directly for a recent version)

cd refs/qgraf/
# Typical layout of a downloaded distribution:
#   qgraf-4.0.6.tgz      — distribution tarball
#   qgraf-4.0.6.pdf      — manual
# After extracting:
tar xf qgraf-4.0.6.tgz                        # creates v4.0.6/qgraf-4.0.6.dir/
cd v4.0.6/qgraf-4.0.6.dir
make                                           # builds ./qgraf Fortran executable
```

The related paper `Nogueira1993_JCompPhys105_279.pdf` lives in
`../papers/` (also gitignored; see `refs/papers/README.md`).

## Regenerating the golden master test suite

Requires a working qgraf binary in `./qgraf`:

```bash
cd refs/qgraf/v4.0.6/qgraf-4.0.6.dir
python3 generate_golden_masters.py
# writes golden_masters/{manifest.json, SUMMARY.md, <model>/<case>.out}
```

`scripts/qgraf_golden_master_report.jl` (in the repo root) consumes these
fixtures to report PASS/FAIL status of our Julia port against qgraf's
outputs without requiring qgraf to be installed locally at test time.

## Attribution

When citing qgraf in source code comments (Rule 2 in `CLAUDE.md`), use:

```julia
# Ref: qgraf-4.0.6.f08:<line_number> (<subroutine_name>)
# Ref: Nogueira1993_JCompPhys105_279.pdf §<section>
# Ref: refs/qgraf/ALGORITHM.md §<section>
```
