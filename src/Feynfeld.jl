# Feynfeld: Julia-native, agent-facing, full-stack physics computation suite.
# Six-layer pipeline: Model → Rules → Diagrams → Algebra → Integrals → Evaluate.
# Design choices: see src/v2/DESIGN.md. Rules + vision: see CLAUDE.md + Feynfeld_PRD.md.
module Feynfeld

# ---- Layer 4: Algebra (core) ----
include("v2/coeff.jl")
include("v2/types.jl")
include("v2/colour_types.jl")
include("v2/pair.jl")
include("v2/expr.jl")
include("v2/sp_context.jl")
include("v2/contract.jl")
include("v2/eps_contract.jl")
include("v2/expand_sp.jl")
include("v2/dirac.jl")
include("v2/dirac_trace.jl")
include("v2/dirac_expr.jl")
include("v2/spin_sum.jl")
include("v2/interference.jl")
include("v2/dirac_trick.jl")
include("v2/colour_trace.jl")
include("v2/colour_simplify.jl")

# ---- Layer 1: Model ----
include("v2/model.jl")

# ---- Layer 2: Rules ----
include("v2/rules.jl")

# ---- Layer 3: Diagrams / Channels / Amplitudes ----
include("v2/diagrams.jl")
include("v2/channels.jl")
include("v2/amplitude.jl")
include("v2/gauge_exchange.jl")

# ---- Layer 1d: φ³ Model ----
include("v2/phi3_model.jl")

# ---- Layer 3d: qgraf-faithful topology generator (Strategy C port) ----
# Submodule included FIRST so the legacy canonicality filter can delegate
# to `QgrafPort.is_canonical_feynman` (full per-class lex-next-permutation).
include("v2/qgraf/QgrafPort.jl")

# ---- Layer 3c: Algorithmic diagram generation ----
include("v2/degree_partition.jl")
include("v2/topology_types.jl")
include("v2/topology_filter.jl")
include("v2/topology_enum.jl")
include("v2/vertex_check.jl")
include("v2/field_assign.jl")
include("v2/diagram_gen.jl")

# ---- Layer 3b: Loop diagrams ----
include("v2/loop_channels.jl")
include("v2/loop_amplitude.jl")
include("v2/loop_interference.jl")

# ---- Layer 1b: QCD Model ----
include("v2/qcd_model.jl")

# ---- Layer 1c: EW Model ----
include("v2/ew_model.jl")

# ---- Layer 4b: Polarization sums ----
include("v2/polarization_sum.jl")

# ---- Layer 5: Integrals ----
include("v2/pave.jl")
include("v2/c0_analytical.jl")
include("v2/d0_collier.jl")
include("v2/b0_eval.jl")
include("v2/pave_eval.jl")
include("v2/d_tensor.jl")
include("v2/sp_lookup.jl")
include("v2/tid.jl")
include("v2/schwinger.jl")
include("v2/vertex.jl")
include("v2/running_alpha.jl")
include("v2/ew_parameters.jl")
include("v2/ew_cross_section.jl")

# ---- Layer 6: Evaluate ----
include("v2/eps_evaluate.jl")
include("v2/cross_section.jl")

# ---- Layer 6b: NLO assembly (depends on cross_section.jl) ----
include("v2/nlo_box.jl")

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
export spin_sum_amplitude_squared, spin_sum_interference
export DiracExpr, dirac_trick
export colour_trace, colour_delta_trace, contract_colour
export casimir_fundamental, casimir_adjoint, trace_normalization

# ---- Exports: Model (Layer 1) ----
export GaugeGroup, U1, SU
export FieldSpecies, Fermion, Boson, Scalar
export Field, fermion, vector_boson, scalar
export AbstractModel, QEDModel, qed_model, qed1_model, qed3_model, QCDModel, qcd_model, triple_gauge_vertex
export EWModel, ew_model
export model_name, model_fields, gauge_groups
export get_field, fermion_fields, boson_fields
export mass_trait, charge_trait, species, Massive, Massless, Charged, Neutral

# ---- Exports: Rules (Layer 2) ----
export VertexRule, FeynmanRules, feynman_rules
export vertex_factor, vertex_structure, propagator_num, gauge_coupling_phase

# ---- Exports: Diagrams / Channels (Layer 3) ----
export ExternalLeg
export TreeChannel, tree_channels, vertex_legs, build_amplitude, propagator_momentum
export Phi3Model, phi3_model
export count_diagrams, generate_tree_channels
export LoopChannel, box_channels, build_loop_box_amplitude, BoxDenominators
export spin_sum_tree_loop_interference
export evaluate_box_integral

# ---- Exports: Integrals (Layer 5) ----
export polarization_sum
export PaVe, A0, B0, B1, B00, B11, C0, C1, C2, D0, D1, D2, D3
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
export evaluate_single_box_channel, evaluate_box_channels, born_virtual_box
export Mandelstam, sp_context_from_mandelstam
export CrossSectionProblem, solve_tree, solve_tree_pipeline, evaluate_m_squared
export evaluate_numeric, sp_values_2to2
export dsigma_domega, sigma_total_tree_ee_mumu

end # module Feynfeld
