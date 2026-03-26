module Feynfeld

using LinearAlgebra
import TensorGR

# Layer 0: Foundation — dimensional regularisation and TensorGR bridge
include("algebra/dimensions.jl")

# Layer 1: Algebra — the type system (port first)
# Lorentz · Dirac · Colour · Tensor
# Contraction, traces, simplify, PV reduction
include("algebra/types.jl")
include("algebra/minkowski.jl")

# Layer 2: Integrals
# PV scalar functions A₀ B₀ C₀ D₀
# include("integrals/integrals.jl")

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
