# `refs/papers/` — cited literature

Journal-copyrighted PDFs are gitignored. Acquire via TIB VPN (institutional
access at LUH) or DoltHub/JSTOR if open-access. Our fetch scripts (committed)
automate retrieval via Playwright + authenticated browser session.

## Committed here

- `README.md` — this file
- `fetch_nogueira.mjs` — Playwright script to download Nogueira (1993) from
  ScienceDirect, reusing the persistent Chromium profile under
  `~/Projects/Sturm.jl/docs/literature/quantum_simulation/.browser-profile`.
  Usage: `node refs/papers/fetch_nogueira.mjs` (TIB VPN must be active).

## Minimum paper set (all gitignored)

| File | Reference | Acquisition |
|------|-----------|-------------|
| `Nogueira1993_JCompPhys105_279.pdf` | P. Nogueira, J. Comput. Phys. 105 (1993) 279 — qgraf algorithm | `node refs/papers/fetch_nogueira.mjs` |
| `Denner1993_FortschrPhys41.pdf` | A. Denner, Fortschr. Phys. 41 (1993) 307 — one-loop techniques, PV reduction | Wiley via TIB |
| `tHooftVeltman1979_NuclPhysB153.pdf` | 't Hooft, Veltman, Nucl. Phys. B153 (1979) 365 — scalar integrals | Elsevier via TIB |
| `PassarinoVeltman1979_NuclPhysB160.pdf` | Passarino, Veltman, Nucl. Phys. B160 (1979) 151 — PV decomposition | Elsevier via TIB |
| `MertigBohmDenner1991_FeynCalc_CPC64.pdf` | Mertig, Böhm, Denner, Comput. Phys. Commun. 64 (1991) 345 — FeynCalc | Elsevier via TIB |
| `Shtabovenko2024_FeynCalc10_2312.14089.pdf` | Shtabovenko et al., arXiv:2312.14089 — FeynCalc 10 | arXiv (open access) |
| `vanOldenborgh1990_ZPhysC46.pdf` | van Oldenborgh, Z. Phys. C46 (1990) 425 — LoopTools origin | Springer via TIB |

PDG 2024 review chapters and the P&S textbook (djvu) are also in this
directory, likewise gitignored.

## Writing a new fetch script

Pattern (see `fetch_nogueira.mjs`):

1. Import `chromium.launchPersistentContext` from an existing node_modules
   path (currently `~/Projects/qvls-sturm/viz/node_modules/playwright`).
2. Launch headless = false. User solves CAPTCHA / institutional login.
3. Wait on a DOM signal that institutional access is live
   (`meta[name="citation_pdf_url"]` on ScienceDirect is reliable).
4. Request the PDF via `page.request.get(pdfUrl, { Referer })` using
   session cookies.
5. Validate `body[0..5] == '%PDF-'` before saving.
