# Layer 3: External leg definition for Feynman diagrams.
#
# ExternalLeg is the shared type used by channels.jl, amplitude.jl,
# and cross_section.jl. It specifies a particle entering or leaving
# the scattering process.

struct ExternalLeg
    field_name::Symbol
    momentum::Momentum
    incoming::Bool
    antiparticle::Bool
    mass::Rational{Int}
end

# Backward-compatible constructor: mass defaults to 0 (massless)
ExternalLeg(name, mom, incoming, anti) = ExternalLeg(name, mom, incoming, anti, 0//1)

function Base.show(io::IO, l::ExternalLeg)
    dir = l.incoming ? "in" : "out"
    anti = l.antiparticle ? "bar" : ""
    print(io, "$(l.field_name)$(anti)($(dir), $(l.momentum))")
end
