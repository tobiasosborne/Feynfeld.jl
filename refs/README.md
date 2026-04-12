# `refs/` — Reference Sources

Most of this directory is **gitignored**. Agents working on other devices
must re-acquire the ignored pieces using the instructions below.

## What IS committed

| Path | Purpose |
|------|---------|
| `refs/README.md` | this file |
| `refs/qgraf/README.md` | qgraf acquisition instructions |
| `refs/qgraf/ALGORITHM.md` | cleanroom qg21/qpg11/qgen spec (written from source analysis) |
| `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/generate_golden_masters.py` | our script; regenerates test fixtures |
| `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/parse_golden_master.py` | our script; parses `.out` → JSON |
| `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/` | 104 generated test fixtures (our outputs, used by `scripts/qgraf_golden_master_report.jl`) |
| `refs/papers/README.md` | paper acquisition instructions |
| `refs/papers/fetch_nogueira.mjs` | our Playwright fetch script for Nogueira (1993) |

## What is NOT committed (and why)

| Path | Reason |
|------|--------|
| `refs/qgraf/*.f08`, `*.pdf`, `*.tgz`, `v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08`, `qgraf-4.0.6.pdf` | qgraf source is © Paulo Nogueira, explicitly marked "not to be redistributed without explicit permission" |
| `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/models/`, `styles/`, `q4fc/`, module binaries | part of the qgraf distribution, same licence constraint |
| `refs/papers/*.pdf` | journal-copyrighted PDFs (Elsevier, Springer, APS, …) — fetched via TIB VPN for local reading only, not redistributable |
| `refs/FeynCalc/`, `refs/FeynArts/`, `refs/FeynRules/`, `refs/LoopTools/`, `refs/FormCalc/` | upstream reference codebases used as porting oracles; each has its own licence (mostly LGPL/GPL) and is mirrored locally via git clone |

## How to re-acquire

See `refs/qgraf/README.md` for qgraf and `refs/papers/README.md` for papers.
For FeynCalc / FeynArts / FeynRules / LoopTools / FormCalc:

```bash
cd refs/
git clone https://github.com/FeynCalc/feyncalc FeynCalc
wget https://wwwth.mpp.mpg.de/members/hahn/FeynArts.tar.gz && tar xf FeynArts.tar.gz
# FeynRules: http://feynrules.irmp.ucl.ac.be/ (registration)
# LoopTools: https://feynarts.de/looptools/  (tarball)
# FormCalc: https://feynarts.de/formcalc/   (tarball)
```
