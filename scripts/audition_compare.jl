# Compare audition (A/B/C) vs legacy on a battery of cases.
using Feynfeld
using Feynfeld.QgrafPort: count_dedup_burnside, count_dedup_canonical,
                              count_dedup_prefilter

cases = [
    ("phi3 tree φ→φφ",      phi3_model(), [:phi], [:phi, :phi], 0),
    ("phi3 tree φφ→φφ",     phi3_model(), [:phi, :phi], [:phi, :phi], 0),
    ("phi3 tree φ→φφφ",     phi3_model(), [:phi], [:phi, :phi, :phi], 0),
    ("phi3 1L φ→φφ",        phi3_model(), [:phi], [:phi, :phi], 1),
    ("phi3 1L φφ→φφ",       phi3_model(), [:phi, :phi], [:phi, :phi], 1),
    ("phi3 1L φ→φφφ",       phi3_model(), [:phi], [:phi, :phi, :phi], 1),
    ("QED ee→μμ tree",      qed_model(),  [:e, :e], [:mu, :mu], 0),
    ("QED ee→μμ 1L 2gen",   qed_model(),  [:e, :e], [:mu, :mu], 1),
    ("QED ee→ee tree (1g)", qed1_model(), [:e, :e], [:e, :e], 0),
    ("QED eγ→eγ tree (1g)", qed1_model(), [:e, :gamma], [:e, :gamma], 0),
]

println(rpad("CASE", 26), "  ", "LEG", "    A", "    B", "    C")
for (name, m, in_, out_, loops) in cases
    leg = count_diagrams(m, in_, out_; loops=loops)
    a   = count_dedup_burnside(m,  in_, out_; loops=loops)
    b   = count_dedup_canonical(m, in_, out_; loops=loops)
    c   = count_dedup_prefilter(m, in_, out_; loops=loops)
    flag = (leg == a == b == c) ? "✓" : "✗"
    println(rpad(name, 26), "  ", lpad(leg, 3), "  ", lpad(a, 3), "  ",
            lpad(b, 3), "  ", lpad(c, 3), "  ", flag)
end
