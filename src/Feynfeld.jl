module Feynfeld

using LinearAlgebra
import TensorGR

# Layer 0: Foundation — dimensional regularisation and TensorGR bridge
include("algebra/dimensions.jl")

# Layer 1: Algebra — the type system
# Core types (LorentzIndex, Dirac, Colour, PaVe)
include("algebra/types.jl")
# Momentum (Pair slot argument)
include("algebra/momentum.jl")
# Pair (universal Lorentz bilinear) + convenience constructors
include("algebra/pair.jl")
# Scalar product context
include("algebra/sp_context.jl")
# Dirac algebra types (gamma matrices, spinors, chains)
include("algebra/dirac_types.jl")
# Non-commutative Dirac chain (DOT) and DotSimplify
include("algebra/dirac_chain.jl")
# DiracTrick: core Dirac simplification rules
include("algebra/dirac_trick.jl")
# DiracTrace: trace evaluation
include("algebra/dirac_trace.jl")
# DiracOrder: normal ordering
include("algebra/dirac_order.jl")
# Dirac equation at chain boundaries
include("algebra/dirac_equation.jl")
# DiracSimplify: master orchestrator
include("algebra/dirac_simplify.jl")
# Dirac gamma5 scheme system (NDR/BMHV/Larin)
include("algebra/dirac_scheme.jl")
# SU(N) colour algebra
include("algebra/colour_types.jl")
include("algebra/colour_trace.jl")
include("algebra/colour_simplify.jl")
# ExpandScalarProduct (bilinear expansion)
include("algebra/expand_sp.jl")
# Levi-Civita tensor
include("algebra/eps.jl")
# Algebraic expression tree (AlgTerm/AlgSum)
include("algebra/alg_expr.jl")
# Lorentz index contraction
include("algebra/contract.jl")
# AlgSum operations (contract, expand_sp, dirac_trace_alg)
include("algebra/alg_ops.jl")
# FermionSpinSum: completeness relations for spin-averaged |M|²
include("algebra/fermion_spin_sum.jl")
# TensorGR bridge (Minkowski registry)
include("algebra/minkowski.jl")

# Layer 2: Integrals — PaVe symbols and propagator types
include("integrals/pave.jl")
include("integrals/feynamp_denominator.jl")
include("integrals/tdec.jl")
include("integrals/pave_reduce.jl")

# Layer 3: Model
# Lagrangian, fields, parameters, symmetries
# include("model/model.jl")

# Layer 4: Rules
# Lagrangian → Feynman rules, vertex extraction
# include("rules/rules.jl")

# Layer 5: Diagrams
# Topology generation, field insertion, amplitude construction
# include("diagrams/diagrams.jl")

# Layer 6: Evaluate
# |M|², cross-sections, decay rates, coupling RGE, codegen
# include("evaluate/evaluate.jl")

end # module Feynfeld
