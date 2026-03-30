# FeynfeldX: experimental Julia-idiomatic rebuild of Feynfeld.jl
# Key design choices: see DESIGN.md
# Six-layer pipeline: Model → Rules → Diagrams → Algebra → Integrals → Evaluate
module FeynfeldX

# ---- Layer 4: Algebra (core) ----
include("coeff.jl")
include("types.jl")
include("colour_types.jl")
include("pair.jl")
include("expr.jl")
include("sp_context.jl")
include("contract.jl")
include("eps_contract.jl")
include("expand_sp.jl")
include("dirac.jl")
include("dirac_trace.jl")
include("dirac_expr.jl")
include("spin_sum.jl")
include("interference.jl")
include("dirac_trick.jl")
include("colour_trace.jl")
include("colour_simplify.jl")

# ---- Layer 1: Model ----
include("model.jl")

# ---- Layer 2: Rules ----
include("rules.jl")

# ---- Layer 3: Diagrams / Channels / Amplitudes ----
include("diagrams.jl")
include("channels.jl")
include("amplitude.jl")

# ---- Layer 1b: QCD Model ----
include("qcd_model.jl")

# ---- Layer 1c: EW Model ----
include("ew_model.jl")

# ---- Layer 4b: Polarization sums ----
include("polarization_sum.jl")

# ---- Layer 5: Integrals ----
include("pave.jl")
include("c0_analytical.jl")
include("d0_collier.jl")
include("pave_eval.jl")
include("schwinger.jl")
include("vertex.jl")
include("running_alpha.jl")
include("ew_parameters.jl")
include("ew_cross_section.jl")

# ---- Layer 6: Evaluate ----
include("cross_section.jl")

# ---- Exports: Algebra ----
export DimPoly, DIM, DIM_MINUS_4, Coeff, evaluate_dim, normalise_coeff, mul_coeff, add_coeff
export PhysicsIndex, Dim4, DimD, DimDm4, LorentzIndex, Momentum, MomentumSum, momentum_sum
export MetricTensor, FourVector, ScalarProduct
export pair, SP, FV, MT, SPD, FVD, MTD
export Eps
export AdjointIndex, FundIndex, SUNT, SUNDelta, FundDelta, SUNF, SUND, ColourChain
export AlgSum, AlgFactor, FactorKey, alg, alg_from_factors
export SPContext, set_sp, get_sp, with_sp, sp_context, evaluate_sp, CURRENT_SP
export contract, eps_contract, expand_scalar_product, substitute_index
export DiracGamma, DiracSlot, LISlot, MomSlot, MomSumSlot, Gamma5Slot
export ProjPSlot, ProjMSlot, IdSlot
export DiracChain, DiracElement, Spinor, SpinorKind
export UKind, VKind, UBarKind, VBarKind
export GA, GAD, GS, GA5, GA6, GA7
export u, v, ubar, vbar
export dot, gammas, gamma_pair, lorentz_index
export dirac_trace
export fermion_spin_sum, spin_sum_amplitude_squared, spin_sum_interference
export DiracExpr, dirac_trick
export colour_trace, colour_delta_trace, contract_colour
export casimir_fundamental, casimir_adjoint, trace_normalization

# ---- Exports: Model (Layer 1) ----
export GaugeGroup, U1, SU
export FieldSpecies, Fermion, Boson, Scalar
export Field, fermion, vector_boson, scalar
export AbstractModel, QEDModel, qed_model, QCDModel, qcd_model, triple_gauge_vertex
export EWModel, ew_model
export model_name, model_fields, gauge_groups
export get_field, fermion_fields, boson_fields
export mass_trait, charge_trait, species, Massive, Massless, Charged, Neutral

# ---- Exports: Rules (Layer 2) ----
export VertexRule, PropagatorRule, FeynmanRules, feynman_rules
export vertex_factor, vertex_structure, propagator_num

# ---- Exports: Diagrams / Channels (Layer 3) ----
export ExternalLeg
export TreeChannel, tree_channels, vertex_legs, build_amplitude, propagator_momentum

# ---- Exports: Integrals (Layer 5) ----
export polarization_sum
export PaVe, A0, B0, B1, B00, B11, C0, C1, C2, D0
export evaluate
export schwinger_correction, vacuum_polarization, sigma_nlo_ee_mumu
export vertex_f2_zero, vertex_f2
export SM_LEPTONS, SM_QUARKS, SM_FERMIONS
export delta_alpha, running_alpha, sigma_improved_ee_mumu
export polarization_sum_massive
export EW_M_W, EW_M_Z, EW_SIN2_W, EW_COS2_W, EW_SIN_W, EW_COS_W, EW_ALPHA
export EW_GV_E, EW_GA_E, EW_GL_E, EW_GR_E
export sigma_ee_ww

# ---- Exports: Evaluate (Layer 6) ----
export Mandelstam, sp_context_from_mandelstam
export CrossSectionProblem, solve_tree, evaluate_m_squared
export dsigma_domega, sigma_total_tree_ee_mumu

end # module FeynfeldX
