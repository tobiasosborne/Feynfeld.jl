#!/usr/bin/env python3
"""
Parse qgraf array.sty output into structured JSON for Julia test generation.

Each diagram becomes a dict with:
- index, sign, symmetry_factor
- external particles (in/out with momenta)
- propagators (field, index, momentum)
- vertices (list of (field, index, momentum) at each vertex)
- topology signature (sorted vertex degrees, propagator count)

Usage: python3 parse_golden_master.py golden_masters/qed2/eminus_eplus_TO_muminus_muplus_1L_onepi.out
"""

import re, json, sys
from pathlib import Path


def parse_momentum(s):
    """Parse momentum string like 'k1-p1+q1' into structured form."""
    s = s.strip()
    terms = []
    # Split on +/- keeping the sign
    parts = re.findall(r'[+-]?[a-z]\d*', s)
    for p in parts:
        sign = -1 if p.startswith('-') else 1
        name = p.lstrip('+-')
        terms.append({"name": name, "sign": sign})
    return terms


def parse_field_entry(s):
    """Parse 'field(index,momentum)' into dict."""
    m = re.match(r'(\w+)\((-?\d+),([^)]+)\)', s.strip())
    if not m:
        return None
    return {
        "field": m.group(1),
        "index": int(m.group(2)),
        "momentum": m.group(3).strip(),
    }


def parse_vertex(s):
    """Parse vertex content: 'field1(idx1,mom1),field2(idx2,mom2),...'"""
    # Match individual field entries within the vertex
    entries = re.findall(r'(\w+\(-?\d+,[^)]+\))', s)
    return [parse_field_entry(e) for e in entries if parse_field_entry(e)]


def parse_diagram(index, body):
    """Parse a single diagram body from array.sty format."""
    diag = {
        "index": index,
        "external_in": [],
        "external_out": [],
        "propagators": [],
        "vertices": [],
    }

    # Sign and symmetry factor
    sf_match = re.search(r'\(([+-]\d+(?:/\d+)?)\)', body)
    if sf_match:
        sf_str = sf_match.group(1)
        diag["sign"] = 1 if sf_str.startswith('+') else -1
        diag["symmetry_factor"] = sf_str

    # External polarizations
    for m in re.finditer(r'pol\((\w+)\((-?\d+),([^)]+)\)\)', body):
        entry = {
            "field": m.group(1),
            "index": int(m.group(2)),
            "momentum": m.group(3).strip(),
        }
        # Negative odd indices = in, negative even = out (qgraf convention)
        # Actually: in particles get listed first in the output
        diag["external_in" if abs(entry["index"]) % 2 == 1 else "external_out"].append(entry)

    # Propagators
    for m in re.finditer(r'prop\((\w+)\((-?\d+),([^)]+)\)\)', body):
        diag["propagators"].append({
            "field": m.group(1),
            "index": int(m.group(2)),
            "momentum": m.group(3).strip(),
        })

    # Vertices
    for m in re.finditer(r'vrtx\(([^;]*?)\)\*', body):
        vraw = m.group(1)
        diag["vertices"].append(parse_vertex(vraw))

    # Also catch last vertex (before semicolon)
    last_v = re.search(r'vrtx\(([^;]*?)\)\s*;', body)
    if last_v:
        diag["vertices"].append(parse_vertex(last_v.group(1)))

    # Topology signature
    vertex_degrees = sorted(len(v) for v in diag["vertices"])
    diag["topology"] = {
        "n_external": len(diag["external_in"]) + len(diag["external_out"]),
        "n_propagators": len(diag["propagators"]),
        "n_vertices": len(diag["vertices"]),
        "vertex_degrees": vertex_degrees,
        "n_loops": len(diag["propagators"]) - len(diag["vertices"]) + 1,
        "propagator_fields": sorted(p["field"] for p in diag["propagators"]),
    }

    return diag


def parse_output_file(path):
    """Parse complete qgraf array.sty output file."""
    text = Path(path).read_text()

    # Extract metadata from header comments
    metadata = {}
    for line in text.split('\n'):
        if line.startswith('%') and '=' in line:
            key, _, val = line.lstrip('% ').partition('=')
            metadata[key.strip()] = val.strip().rstrip(';').strip().strip("'")

    # Split into diagrams
    parts = re.split(r'\ba\((\d+)\)\s*:=', text)
    diagrams = []
    for i in range(1, len(parts), 2):
        idx = int(parts[i])
        body = parts[i + 1].strip()
        diagrams.append(parse_diagram(idx, body))

    return {
        "metadata": metadata,
        "n_diagrams": len(diagrams),
        "diagrams": diagrams,
    }


def topology_fingerprint(diag):
    """Generate a canonical topology fingerprint for duplicate detection."""
    t = diag["topology"]
    return (
        t["n_external"],
        t["n_propagators"],
        t["n_vertices"],
        tuple(t["vertex_degrees"]),
        tuple(t["propagator_fields"]),
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: parse_golden_master.py <output_file>")
        sys.exit(1)

    result = parse_output_file(sys.argv[1])

    if "--json" in sys.argv:
        print(json.dumps(result, indent=2))
    else:
        print(f"File: {sys.argv[1]}")
        print(f"Model: {result['metadata'].get('model', '?')}")
        print(f"Process: {result['metadata'].get('in', '?')} → {result['metadata'].get('out', '?')}")
        print(f"Loops: {result['metadata'].get('loops', '?')}")
        print(f"Diagrams: {result['n_diagrams']}")
        print()

        # Topology summary
        from collections import Counter
        topo_counts = Counter()
        for d in result["diagrams"]:
            fp = topology_fingerprint(d)
            topo_counts[fp] += 1

        print(f"Distinct topologies: {len(topo_counts)}")
        for fp, count in topo_counts.most_common():
            n_ext, n_prop, n_vert, vdeg, pfields = fp
            n_loops = n_prop - n_vert + 1
            print(f"  {count}× : {n_loops}L, {n_prop} props {list(pfields)}, "
                  f"vertices {list(vdeg)}")
