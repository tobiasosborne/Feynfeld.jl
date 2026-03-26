# Feynfeld.jl — Dimensional regularisation algebra
#
# Implements the BMHV (Breitenlohner-Maison / 't Hooft-Veltman) dimensional
# projection rules. Every Lorentz object lives in one of three dimension slots:
#   4    — purely 4-dimensional
#   D    — full D-dimensional (D = 4 - 2ε)
#   D-4  — evanescent (D-4)-dimensional
#
# The projection algebra is:
#   4 ∩ 4     = 4        D ∩ D     = D        (D-4) ∩ (D-4) = D-4
#   4 ∩ D     = 4        D ∩ (D-4) = D-4      4 ∩ (D-4)     = 0
#
# Ported from FeynCalc SharedObjects.m Pair/DiracGamma dimensional rules.

export Dim4, DimD, DimDm4, DimSlot, dim_contract, dim_trace, to_dim

"""
    Dim4

Marker for 4-dimensional objects.
"""
struct Dim4 end

"""
    DimD

Marker for D-dimensional objects (dimensional regularisation).
"""
struct DimD end

"""
    DimDm4

Marker for (D-4)-dimensional evanescent objects.
"""
struct DimDm4 end

const DimSlot = Union{Dim4, DimD, DimDm4}

"""
    dim_contract(d1, d2) -> Union{DimSlot, Nothing}

Compute the BMHV dimensional projection of two dimension slots.
Returns `nothing` when the projection vanishes (4 ∩ (D-4) = 0).

This encodes the complete BMHV algebra for projecting between dimension spaces.
Ref: FeynCalc SharedObjects.m, PairContract.m dimEval rules.
"""
dim_contract(::Dim4,  ::Dim4)  = Dim4()
dim_contract(::DimD,  ::DimD)  = DimD()
dim_contract(::DimDm4, ::DimDm4) = DimDm4()
dim_contract(::Dim4,  ::DimD)  = Dim4()
dim_contract(::DimD,  ::Dim4)  = Dim4()
dim_contract(::DimD,  ::DimDm4) = DimDm4()
dim_contract(::DimDm4, ::DimD)  = DimDm4()
dim_contract(::Dim4,  ::DimDm4) = nothing  # vanishes
dim_contract(::DimDm4, ::Dim4)  = nothing  # vanishes

"""
    dim_trace(d) -> Union{Int, Symbol}

The metric trace g^μ_μ in the given dimension slot.
"""
dim_trace(::Dim4)  = 4
dim_trace(::DimD)  = :D
dim_trace(::DimDm4) = :(D - 4)

"""
    to_dim(x) -> DimSlot

Convert FeynCalc-style dimension specification to a DimSlot.
- `4` or omitted → Dim4()
- `:D` → DimD()
- `:(D-4)` → DimDm4()
"""
to_dim(x::Int) = x == 4 ? Dim4() : error("Unsupported integer dimension $x; use 4, :D, or :(D-4)")
to_dim(s::Symbol) = s === :D ? DimD() : error("Unsupported symbolic dimension :$s; use :D or :(D-4)")
function to_dim(ex::Expr)
    (ex == :(D - 4)) && return DimDm4()
    error("Unsupported dimension expression $ex; use 4, :D, or :(D-4)")
end
to_dim(::Dim4) = Dim4()
to_dim(::DimD) = DimD()
to_dim(::DimDm4) = DimDm4()
