#!/usr/bin/env julia
# Diagnostic: run every golden-master case through count_diagrams_qg21 (the
# Strategy C qg21 port) and report PASS / FAIL / SKIP. Sibling to
# qgraf_golden_master_report.jl which uses the legacy count_diagrams path.
#
# Usage: julia --project=. scripts/qgraf_golden_master_qg21_report.jl [max_loops]
#
# Env: QGRAF_MAX_SECONDS (default 60) — per-case time budget.

using Feynfeld
using Feynfeld.QgrafPort: count_diagrams_qg21

const SUMMARY = joinpath(@__DIR__, "..", "refs", "qgraf", "v4.0.6",
                         "qgraf-4.0.6.dir", "golden_masters", "SUMMARY.md")

const FIELD_MAP = Dict(
    "phi"       => :phi,
    "e_minus"   => :e,   "e_plus"      => :e,
    "photon"    => :gamma,
    "mu_minus"  => :mu,  "mu_plus"     => :mu,
    "tau_minus" => :tau, "tau_plus"    => :tau,
    "quark"     => :q,   "antiquark"   => :q,
    "gluon"     => :g,
    "ghost"     => :ghost,
)

const MODEL_CTOR = Dict(
    "phi3" => () -> phi3_model(),
    "qed1" => () -> qed1_model(),
    "qed2" => () -> qed_model(),
    "qed3" => () -> qed3_model(),
    "qcd"  => () -> qcd_model(),
)

parse_fields(side::AbstractString) = begin
    parts = [strip(s) for s in split(side) if !isempty(strip(s))]
    syms = Symbol[]
    for p in parts
        haskey(FIELD_MAP, p) || return nothing
        push!(syms, FIELD_MAP[p])
    end
    syms
end

parse_options(opt::AbstractString) =
    opt == "—" ? Symbol[] : [Symbol(strip(s)) for s in split(opt, ",")]

struct Case
    model::String
    in_fields::Vector{Symbol}
    out_fields::Vector{Symbol}
    loops::Int
    options::Vector{Symbol}
    expected::Int
end

complexity(c::Case) = c.loops * 1000 +
                      (length(c.in_fields) + length(c.out_fields)) * 10 +
                      length(c.options)

function parse_cases()
    cases = Case[]
    for line in eachline(SUMMARY)
        startswith(line, "| ") || continue
        cells = [strip(s) for s in split(line, "|")][2:end-1]
        length(cells) == 5 || continue
        cells[1] in ("Model", "-------") && continue

        model = cells[1]
        process = cells[2]
        loops = tryparse(Int, cells[3])
        loops === nothing && continue

        arrow_idx = findfirst("→", process)
        arrow_idx === nothing && continue
        in_str  = strip(process[1:prevind(process, first(arrow_idx))])
        out_str = strip(process[nextind(process, last(arrow_idx)):end])
        in_fields  = parse_fields(in_str)
        out_fields = parse_fields(out_str)
        (in_fields === nothing || out_fields === nothing) && continue

        opts = parse_options(cells[4])
        n = cells[5] == "FAIL" ? -1 : something(tryparse(Int, cells[5]), -2)
        n == -2 && continue

        push!(cases, Case(model, in_fields, out_fields, loops, opts, n))
    end
    cases
end

const SUPPORTED_OPTS = Set([:onepi, :nosbridge, :notadpole, :onshell,
                             :nosnail, :onevi, :noselfloop, :nodiloop,
                             :noparallel])

function run_case(c::Case, max_seconds::Float64)
    haskey(MODEL_CTOR, c.model) || return (:skip, "no model ctor for $(c.model)", 0.0)
    c.expected == -1 && return (:skip, "qgraf FAIL case", 0.0)

    unsupported = setdiff(c.options, SUPPORTED_OPTS)
    isempty(unsupported) || return (:skip, "unsupported options: $unsupported", 0.0)

    length(c.in_fields) + length(c.out_fields) == 2 && c.loops == 0 &&
        return (:skip, "1→1 0-loop has 0 diagrams trivially", 0.0)

    m = MODEL_CTOR[c.model]()
    kw = (; loops = c.loops,
           onepi      = :onepi      in c.options,
           nosbridge  = :nosbridge  in c.options,
           notadpole  = :notadpole  in c.options,
           onshell    = :onshell    in c.options,
           nosnail    = :nosnail    in c.options,
           onevi      = :onevi      in c.options,
           noselfloop = :noselfloop in c.options,
           nodiloop   = :nodiloop   in c.options,
           noparallel = :noparallel in c.options)
    local got
    t = 0.0
    try
        t = @elapsed got = count_diagrams_qg21(m, c.in_fields, c.out_fields; kw...)
    catch e
        return (:error, "$(typeof(e)): $(first(sprint(showerror, e), 120))", t)
    end
    if t > max_seconds
        return (:slow, "over budget ($(round(t, digits=1))s > $(max_seconds)s): got $got, want $(c.expected)", t)
    end
    got == c.expected ? (:pass, "$got", t) : (:fail, "got $got, want $(c.expected)", t)
end

_status_tag(s::Symbol) = (pass="PASS ", fail="FAIL ", error="ERROR", skip="SKIP ", slow="SLOW ")[s]

function main()
    max_loops = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 3
    max_seconds = parse(Float64, get(ENV, "QGRAF_MAX_SECONDS", "60"))

    cases = parse_cases()
    sort!(cases; by=complexity)

    println("qgraf golden-master report — QG21 PATH (count_diagrams_qg21)")
    println("$(length(cases)) total cases, filter: loops ≤ $max_loops, max $(max_seconds)s per case")
    println()
    println("[ idx/total] STATUS  label $(repeat(" ", 55)) time    result")
    println("-"^110)

    buckets = Dict{Symbol,Int}(k => 0 for k in (:pass, :fail, :skip, :error, :slow))
    running = Tuple{Case,Symbol,String,Float64}[]

    for (i, c) in enumerate(cases)
        if c.loops > max_loops
            buckets[:skip] += 1
            continue
        end
        opts = isempty(c.options) ? "-" : join(c.options, ",")
        proc = "$(join(c.in_fields, ' ')) → $(join(c.out_fields, ' '))"
        label = "$(c.model) | $proc | $(c.loops)L | $opts"
        flush(stdout)
        status, msg, t = run_case(c, max_seconds)
        buckets[status] += 1
        push!(running, (c, status, msg, t))
        tstr = t > 0 ? "$(lpad(round(t, digits=2), 5))s" : "    -"
        println("[$(lpad(i,3))/$(length(cases))] $(_status_tag(status)) $(rpad(label, 60)) $tstr  $msg")
        flush(stdout)
    end

    println()
    println("="^70)
    println("SUMMARY (QG21 PATH): PASS=$(buckets[:pass])  FAIL=$(buckets[:fail])  SKIP=$(buckets[:skip])  SLOW=$(buckets[:slow])  ERROR=$(buckets[:error])")
    println("="^70)

    for status in (:fail, :error, :slow)
        rows = [(c, s, m, t) for (c, s, m, t) in running if s == status]
        isempty(rows) && continue
        println("\n─── $(uppercase(string(status))) ($(length(rows))) ───")
        for (c, _, m, t) in rows
            opts = isempty(c.options) ? "-" : join(c.options, ",")
            proc = "$(join(c.in_fields, ' ')) → $(join(c.out_fields, ' '))"
            tstr = t > 0 ? " ($(round(t, digits=2))s)" : ""
            println("  $(c.model) | $proc | $(c.loops)L | $opts$tstr  →  $m")
        end
    end
end

main()
