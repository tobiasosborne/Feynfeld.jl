#!/usr/bin/env python3
"""
Generate comprehensive golden master test suite from qgraf-4.0.6.

Runs qgraf on a systematic battery of (model, process, loops, options)
configurations and captures structured output for red-green TDD.

Output: golden_masters/<model>/<process>_<loops>L[_<options>].out
        golden_masters/manifest.json  (diagram counts + metadata)
"""

import subprocess, json, os, re, sys
from pathlib import Path

QGRAF = "./qgraf"
GM_DIR = Path("golden_masters")
GM_DIR.mkdir(exist_ok=True)

# ── Models ──────────────────────────────────────────────────────────

MODELS = {
    "phi3": {"file": "phi3", "submodel": None},
    "qed1": {"file": "qedx", "submodel": "qed1"},
    "qed2": {"file": "qedx", "submodel": "qed2"},
    "qed3": {"file": "qedx", "submodel": "qed3"},
    "qcd":  {"file": "qcd",  "submodel": None},
}

# ── Test battery ────────────────────────────────────────────────────
# Each entry: (model, in_particles, out_particles, max_loops, options_list)
# options_list is a list of option strings to try (empty string = no options)

BATTERY = [
    # ═══ phi3: pure graph theory, simplest model ═══
    # 1→2 decay at increasing loop order
    ("phi3", ["phi"],         ["phi", "phi"],                4, [""]),
    ("phi3", ["phi"],         ["phi", "phi"],                2, ["onepi"]),
    ("phi3", ["phi"],         ["phi", "phi"],                2, ["notadpole"]),
    ("phi3", ["phi"],         ["phi", "phi"],                2, ["noselfloop"]),
    ("phi3", ["phi"],         ["phi", "phi"],                2, ["nosnail"]),
    ("phi3", ["phi"],         ["phi", "phi"],                2, ["nosigma"]),
    ("phi3", ["phi"],         ["phi", "phi"],                2, ["onepi", "noselfloop"]),
    # 2→2 scattering
    ("phi3", ["phi", "phi"],  ["phi", "phi"],                3, [""]),
    ("phi3", ["phi", "phi"],  ["phi", "phi"],                2, ["onepi"]),
    ("phi3", ["phi", "phi"],  ["phi", "phi"],                1, ["notadpole"]),
    ("phi3", ["phi", "phi"],  ["phi", "phi"],                1, ["nosnail"]),
    ("phi3", ["phi", "phi"],  ["phi", "phi"],                1, ["nosigma"]),
    ("phi3", ["phi", "phi"],  ["phi", "phi"],                1, ["onevi"]),
    ("phi3", ["phi", "phi"],  ["phi", "phi"],                1, ["nodiloop"]),
    ("phi3", ["phi", "phi"],  ["phi", "phi"],                1, ["noparallel"]),
    # Higher multiplicity
    ("phi3", ["phi"],         ["phi", "phi", "phi"],         2, [""]),
    ("phi3", ["phi"],         ["phi", "phi", "phi"],         1, ["onepi"]),
    ("phi3", ["phi", "phi"],  ["phi", "phi", "phi"],         1, [""]),
    ("phi3", ["phi", "phi"],  ["phi", "phi", "phi"],         0, ["onepi"]),
    ("phi3", ["phi"],         ["phi", "phi", "phi", "phi"],  1, [""]),
    ("phi3", ["phi", "phi"],  ["phi", "phi", "phi", "phi"],  0, [""]),
    # 3→3 (tree only — combinatorial explosion)
    ("phi3", ["phi", "phi", "phi"], ["phi", "phi", "phi"],   0, [""]),

    # ═══ QED 1-gen: electrons + photons ═══
    # e+e- → γγ (pair annihilation)
    ("qed1", ["e_minus", "e_plus"], ["photon", "photon"],    1, [""]),
    ("qed1", ["e_minus", "e_plus"], ["photon", "photon"],    1, ["onepi"]),
    # Compton: eγ → eγ
    ("qed1", ["e_minus", "photon"], ["e_minus", "photon"],   1, [""]),
    ("qed1", ["e_minus", "photon"], ["e_minus", "photon"],   1, ["onepi"]),
    # Bhabha: e+e- → e+e-
    ("qed1", ["e_minus", "e_plus"], ["e_minus", "e_plus"],   1, [""]),
    ("qed1", ["e_minus", "e_plus"], ["e_minus", "e_plus"],   1, ["onepi"]),
    # γγ → e+e-
    ("qed1", ["photon", "photon"],  ["e_minus", "e_plus"],   1, [""]),
    # e- → e- γ (bremsstrahlung vertex)
    ("qed1", ["e_minus"],           ["e_minus", "photon"],   2, [""]),
    ("qed1", ["e_minus"],           ["e_minus", "photon"],   1, ["onepi"]),
    # γ → e+e- (vertex)
    ("qed1", ["photon"],            ["e_minus", "e_plus"],   2, [""]),
    ("qed1", ["photon"],            ["e_minus", "e_plus"],   1, ["onepi"]),
    # 5-point tree: e+e- → e+e- γ
    ("qed1", ["e_minus", "e_plus"], ["e_minus", "e_plus", "photon"], 0, [""]),
    # e+e- → 3γ
    ("qed1", ["e_minus", "e_plus"], ["photon", "photon", "photon"], 0, [""]),
    # γγ → γγ (light-by-light, purely loop)
    ("qed1", ["photon", "photon"],  ["photon", "photon"],    1, [""]),
    # Vacuum polarization: photon self-energy
    ("qed1", ["photon"],            ["photon"],              2, [""]),
    ("qed1", ["photon"],            ["photon"],              1, ["onepi"]),
    # Electron self-energy
    ("qed1", ["e_minus"],           ["e_minus"],             2, [""]),
    ("qed1", ["e_minus"],           ["e_minus"],             1, ["onepi"]),

    # ═══ QED 2-gen: e + μ (the Feynfeld benchmark process) ═══
    # e+e- → μ+μ- (THE benchmark)
    ("qed2", ["e_minus", "e_plus"], ["mu_minus", "mu_plus"], 1, [""]),
    ("qed2", ["e_minus", "e_plus"], ["mu_minus", "mu_plus"], 1, ["onepi"]),
    ("qed2", ["e_minus", "e_plus"], ["mu_minus", "mu_plus"], 1, ["onepi", "floop"]),
    ("qed2", ["e_minus", "e_plus"], ["mu_minus", "mu_plus"], 1, ["notadpole"]),
    ("qed2", ["e_minus", "e_plus"], ["mu_minus", "mu_plus"], 1, ["nosigma"]),
    ("qed2", ["e_minus", "e_plus"], ["mu_minus", "mu_plus"], 0, [""]),
    # e-μ- → e-μ- (elastic scattering, t-channel only)
    ("qed2", ["e_minus", "mu_minus"], ["e_minus", "mu_minus"], 0, [""]),
    ("qed2", ["e_minus", "mu_minus"], ["e_minus", "mu_minus"], 1, ["onepi"]),
    # μ+μ- → μ+μ-
    ("qed2", ["mu_minus", "mu_plus"], ["mu_minus", "mu_plus"], 0, [""]),
    # μ+μ- → e+e- (crossing of benchmark)
    ("qed2", ["mu_minus", "mu_plus"], ["e_minus", "e_plus"],   0, [""]),
    # Muon self-energy
    ("qed2", ["mu_minus"],            ["mu_minus"],            1, ["onepi"]),

    # ═══ QED 3-gen: e + μ + τ (VP with all flavors) ═══
    ("qed3", ["e_minus", "e_plus"],   ["mu_minus", "mu_plus"], 1, [""]),
    ("qed3", ["e_minus", "e_plus"],   ["tau_minus", "tau_plus"], 0, [""]),
    ("qed3", ["e_minus", "e_plus"],   ["tau_minus", "tau_plus"], 1, ["onepi"]),
    # Photon self-energy with 3 generations
    ("qed3", ["photon"],              ["photon"],              1, ["onepi"]),
    ("qed3", ["photon"],              ["photon"],              2, ["onepi"]),

    # ═══ QCD: quarks + gluons + ghosts ═══
    # qq̄ → gg
    ("qcd", ["quark", "antiquark"],   ["gluon", "gluon"],     0, [""]),
    ("qcd", ["quark", "antiquark"],   ["gluon", "gluon"],     1, ["onepi"]),
    ("qcd", ["quark", "antiquark"],   ["gluon", "gluon"],     1, ["onepi", "floop"]),
    # gg → gg
    ("qcd", ["gluon", "gluon"],       ["gluon", "gluon"],     0, [""]),
    ("qcd", ["gluon", "gluon"],       ["gluon", "gluon"],     1, ["onepi"]),
    # qq̄ → qq̄
    ("qcd", ["quark", "antiquark"],   ["quark", "antiquark"],  0, [""]),
    ("qcd", ["quark", "antiquark"],   ["quark", "antiquark"],  1, ["onepi"]),
    # qg → qg
    ("qcd", ["quark", "gluon"],       ["quark", "gluon"],     0, [""]),
    ("qcd", ["quark", "gluon"],       ["quark", "gluon"],     1, ["onepi"]),
    # gg → qq̄
    ("qcd", ["gluon", "gluon"],       ["quark", "antiquark"], 0, [""]),
    # qq̄ → ggg (5-point tree)
    ("qcd", ["quark", "antiquark"],   ["gluon", "gluon", "gluon"], 0, [""]),
    # gg → ggg (5-point pure glue)
    ("qcd", ["gluon", "gluon"],       ["gluon", "gluon", "gluon"], 0, [""]),
    # Gluon self-energy
    ("qcd", ["gluon"],                ["gluon"],              1, ["onepi"]),
    ("qcd", ["gluon"],                ["gluon"],              1, ["onepi", "floop"]),
    # Quark self-energy
    ("qcd", ["quark"],                ["quark"],              1, ["onepi"]),
    # Ghost contribution check
    ("qcd", ["ghost"],                ["ghost"],              1, ["onepi"]),
]


def make_process_name(in_parts, out_parts):
    """e.g. 'eminus_eplus_TO_muminus_muplus'"""
    def name(p):
        return p.replace("_", "").replace("+", "plus").replace("-", "minus")
    ins = "_".join(name(p) for p in in_parts)
    outs = "_".join(name(p) for p in out_parts)
    return f"{ins}_TO_{outs}"


def make_case_name(model, in_parts, out_parts, loops, options):
    proc = make_process_name(in_parts, out_parts)
    name = f"{proc}_{loops}L"
    if options:
        opts = "_".join(sorted(options))
        name += f"_{opts}"
    return name


def make_control_file(model_info, in_parts, out_parts, loops, options, outfile):
    """Generate a qgraf control file."""
    model_str = model_info["file"]
    if model_info["submodel"]:
        model_str = f'{model_info["submodel"]} // \'{model_info["file"]}\''
    else:
        model_str = f'\'{model_info["file"]}\''

    in_str = ", ".join(
        f"{p}[p{i+1}]" for i, p in enumerate(in_parts)
    )
    out_str = ", ".join(
        f"{p}[q{i+1}]" for i, p in enumerate(out_parts)
    )

    opts = ", ".join(options) if options else ""

    return f""" config = delete, flush ;
 output_dir = 'golden_masters/{outfile.parent.name}/' ;
 style_dir = 'styles/' ;
 model_dir = 'models/' ;
 style = 'array.sty' ;
 output = '{outfile.name}' ;
 model = {model_str} ;
 in = {in_str} ;
 out = {out_str} ;
 loops = {loops} ;
 loop_momentum = k ;
 options = {opts} ;
"""


def parse_diagram_count(stdout):
    """Extract total diagram count from qgraf stdout."""
    m = re.search(r'total\s+\.\.\.\s+(\d+)\s+connected', stdout)
    if m:
        return int(m.group(1))
    # Check for zero
    if "no diagrams" in stdout.lower() or "total   ...   0" in stdout:
        return 0
    # Parse per-loop counts
    counts = re.findall(r'\.\.\.\s+(\d+)\s*$', stdout, re.MULTILINE)
    return sum(int(c) for c in counts) if counts else -1


def parse_output_file(path):
    """Parse array.sty output to extract diagram structures."""
    if not path.exists():
        return []
    text = path.read_text()
    diagrams = []
    # Split on a(N):= pattern
    parts = re.split(r'\ba\((\d+)\)\s*:=', text)
    # parts[0] is header, then alternating index, body
    for i in range(1, len(parts), 2):
        idx = int(parts[i])
        body = parts[i + 1].strip().rstrip(";").strip()
        diag = {"index": idx, "raw": body}

        # Extract components
        diag["sign"] = "+" if body.startswith("(+") else "-"
        diag["symmetry_factor"] = re.search(r'\(([+-]\d+(?:/\d+)?)\)', body)
        if diag["symmetry_factor"]:
            diag["symmetry_factor"] = diag["symmetry_factor"].group(1)

        # Count propagators and vertices
        diag["n_propagators"] = len(re.findall(r'\bprop\(', body))
        diag["n_vertices"] = len(re.findall(r'\bvrtx\(', body))
        diag["n_external_in"] = len(re.findall(r'\bpol\(', body.split("prop")[0])) if "prop" in body else 0

        # Extract field types in propagators
        props = re.findall(r'prop\((\w+)\(', body)
        diag["propagator_fields"] = props

        # Extract vertex structures
        verts = re.findall(r'vrtx\(([^)]+(?:\([^)]*\))*[^)]*)\)', body)
        diag["vertices_raw"] = verts

        diagrams.append(diag)
    return diagrams


def run_case(model_name, in_parts, out_parts, loops, options):
    """Run a single golden master case."""
    model_info = MODELS[model_name]
    case_name = make_case_name(model_name, in_parts, out_parts, loops, options)

    out_dir = GM_DIR / model_name
    out_dir.mkdir(exist_ok=True)
    outfile = out_dir / f"{case_name}.out"

    ctrl = make_control_file(model_info, in_parts, out_parts, loops, options, outfile)
    ctrl_path = Path(f"ctrl.dat")
    ctrl_path.write_text(ctrl)

    try:
        result = subprocess.run(
            [QGRAF, str(ctrl_path)],
            capture_output=True, text=True, timeout=120
        )
        stdout = result.stdout + result.stderr
        count = parse_diagram_count(stdout)

        # Parse output
        diagrams = parse_output_file(outfile)

        return {
            "model": model_name,
            "in": in_parts,
            "out": out_parts,
            "loops": loops,
            "options": options,
            "case_name": case_name,
            "n_diagrams": count,
            "n_diagrams_parsed": len(diagrams),
            "output_file": str(outfile),
            "success": result.returncode == 0 and count >= 0,
            "diagrams_summary": [
                {
                    "index": d["index"],
                    "sign": d["sign"],
                    "symmetry_factor": d["symmetry_factor"],
                    "n_propagators": d["n_propagators"],
                    "n_vertices": d["n_vertices"],
                    "propagator_fields": d["propagator_fields"],
                }
                for d in diagrams
            ],
            "qgraf_stdout": stdout.strip(),
        }
    except subprocess.TimeoutExpired:
        return {
            "model": model_name, "in": in_parts, "out": out_parts,
            "loops": loops, "options": options, "case_name": case_name,
            "n_diagrams": -1, "success": False, "error": "TIMEOUT",
        }
    except Exception as e:
        return {
            "model": model_name, "in": in_parts, "out": out_parts,
            "loops": loops, "options": options, "case_name": case_name,
            "n_diagrams": -1, "success": False, "error": str(e),
        }


def main():
    manifest = {"generator": "qgraf-4.0.6", "cases": []}
    total = 0
    failed = 0

    # Expand battery: each entry generates cases for loops 0..max_loops
    cases = []
    for model, ins, outs, max_loops, opts_list in BATTERY:
        for opts in [opts_list]:  # opts_list is already the option set
            for L in range(max_loops + 1):
                # Skip redundant: options that only matter at loop level
                if L == 0 and any(o in ["onepi", "floop", "nosigma", "nosnail",
                                        "noselfloop", "notadpole", "nodiloop",
                                        "noparallel", "onevi"]
                                  for o in opts_list if o):
                    # Still run L=0 with no options, skip L=0 with loop-only opts
                    # unless it's the only loop level
                    if max_loops > 0:
                        continue
                cases.append((model, ins, outs, L, [o for o in opts_list if o]))

    # Deduplicate
    seen = set()
    unique_cases = []
    for c in cases:
        key = (c[0], tuple(c[1]), tuple(c[2]), c[3], tuple(c[4]))
        if key not in seen:
            seen.add(key)
            unique_cases.append(c)

    print(f"Running {len(unique_cases)} golden master cases...")

    for i, (model, ins, outs, loops, opts) in enumerate(unique_cases):
        case_name = make_case_name(model, ins, outs, loops, opts)
        print(f"  [{i+1}/{len(unique_cases)}] {case_name}", end=" ... ", flush=True)

        result = run_case(model, ins, outs, loops, opts)
        manifest["cases"].append(result)
        total += 1

        if result["success"]:
            print(f"{result['n_diagrams']} diagrams")
        else:
            print(f"FAILED: {result.get('error', 'unknown')}")
            failed += 1

    # Write manifest
    manifest_path = GM_DIR / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    # Write summary table
    summary_path = GM_DIR / "SUMMARY.md"
    with open(summary_path, "w") as f:
        f.write("# Golden Master Test Suite — qgraf-4.0.6\n\n")
        f.write(f"Generated: {len(unique_cases)} cases, {failed} failed\n\n")
        f.write("| Model | Process | Loops | Options | Diagrams |\n")
        f.write("|-------|---------|-------|---------|----------|\n")
        for c in manifest["cases"]:
            proc = f"{' '.join(c['in'])} → {' '.join(c['out'])}"
            opts = ", ".join(c["options"]) if c["options"] else "—"
            count = c["n_diagrams"] if c["success"] else "FAIL"
            f.write(f"| {c['model']} | {proc} | {c['loops']} | {opts} | {count} |\n")

    print(f"\nDone: {total} cases ({total - failed} passed, {failed} failed)")
    print(f"Manifest: {manifest_path}")
    print(f"Summary:  {summary_path}")


if __name__ == "__main__":
    main()
