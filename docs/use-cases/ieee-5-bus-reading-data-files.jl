#=
# [Example of Ordinary Power Flow with IEEE 5 Bus Test System](@id ieee-5-bus-ordinary-power-flow)
=#

using Revise
using ePowerSim


# using SciMLBase

# using SciMLNLSolve: NLSolveJL

using NLSolvers

using NLsolve: nlsolve, converged, OnceDifferentiable

using NLsolve

using NonlinearSolve

using NonlinearSolve: TrustRegion

#---------------------------------------------------

# using DifferentialEquations

# using OrdinaryDiffEq, Sundials, ODEInterfaceDiffEq

using FiniteDiff, LinearSolve

using ForwardDiff, Zygote

using DiffRules

# using PreallocationTools

#---------------------------------------------------

using LinearAlgebra, GenericSchur, Arblib

#---------------------------------------------------

using OrderedCollections: OrderedDict

using Permutations

using SparseArrays, StaticArrays, ComponentArrays

using DataFrames, DataFramesMeta

using JSONTables, JSON, JSON3

using Query, CSV, Tables, XLSX

using Graphs

#---------------------------------------------------

using StatsBase

using Plots

import StatsPlots

using Latexify

using LatexPrint

#---------------------------------------------------
using NamedTupleTools

using Accessors, AccessorsExtra 

#---------------------------------------------------
#---------------------------------------------------


package_dir = pkgdir(ePowerSim)


data =
    joinpath( package_dir, "data")

data_dir = joinpath(data,
             "converted-data")
src_dir =
    joinpath(package_dir, "src")

components_libs_dir =
    joinpath(src_dir,
             "components-lib")

script_dir = @__DIR__


#---------------------------------------------------
#---------------------------------------------------
# Reading network data
#---------------------------------------------------
#---------------------------------------------------

net_type = "dynamic-net-data"

gen_type = "gen-sauer"

avr_type = "avr-t1-cb-sauer"

gov_type = "gov-t1-cb-sauer"

data_ext = "json"

dynamic_net_data_by_components_file =
    "$(net_type)-"*
    "$(gen_type)-"*
    "$(avr_type)-"*
    "$(gov_type)"*
    ".$(data_ext)"

json_net_data_by_components_file =
    dynamic_net_data_by_components_file

#---------------------------------------------------

case_name = "case14"

#---------------------------------------------------

case_data_dir =
   joinpath( data_dir,
             case_name,)

json_case_dir = joinpath(case_data_dir, "json" )


if  (json_net_data_by_components_file == "" ||
    json_net_data_by_components_file == nothing) 

    net_data_by_components_file =
        joinpath(
            json_case_dir,
            "net_data_by_components_file.json")
else

    net_data_by_components_file =
        joinpath(
            json_case_dir,
            json_net_data_by_components_file)

end

#---------------------------------------------------


sim_type  = "sim-reading-data-files"

cd(script_dir)

results_dir =
    joinpath(
        script_dir,
        "results",
        sim_type)

if !(isdir(results_dir))
    
    mkpath(results_dir)
    
end

figure_dir =
    joinpath(
        script_dir,
        "figure",
        sim_type)

if !(isdir(figure_dir))
    
    mkpath(figure_dir)
    
end

tex_filename =
    joinpath(results_dir,
             "$(case_name)-" *
                 "$(sim_type)-" *
                 "sim-results.tex")

sd_dynamics_sim_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "$(sim_type)-" *
            "sim-results.csv")

#---------------------------------------------------
# base setting and some booleans 
#---------------------------------------------------

use_pu_in_PQ    = true

line_data_in_pu = true

use_init_u0     = false

use_nlsolve     = false

with_faults     = false

#---------------------------------------------------
## solvers and settings
#---------------------------------------------------

pf_alg      = NewtonRaphson()

ode_alg     = Rodas4()

dae_alg     = IDA()

# ode_alg   = ImplicitMidpoint()

dt          = 0.0001

Δt          = 1.0 / 2^(4)

abstol      = 1e-12

reltol      = 1e-12

#---------------------------------------------------
# Network pertubation (fault time and data)
#---------------------------------------------------

list_fault_point_from_node_a = [0.3]

list_fault_resistance        = [0.001]

list_no_line_circuit         = [1]

list_edges_to_have_fault     = [ 3 ] # [ 4 ]

clear_fault_selection_list   = [1]

#---------------------------------------------------

on_fault_time = 10.617 # 10.0

Δt_clear_time = 0.5

Δt_generation_adjustment_time =
    Δt_clear_time # 0.2

clear_fault_time =
    on_fault_time +
    Δt_clear_time

line_outage_time =
    on_fault_time

generation_adjustment_time =
    line_outage_time +
    Δt_generation_adjustment_time

#---------------------------------------------------
# sudden load change
#---------------------------------------------------

pertubation_factor = 1.10

restoration_factor = 1.0

pertubation_time   = on_fault_time

restoration_time   = clear_fault_time

Δt1 = 1.5

Δt2 = 1.5

#---------------------------------------------------
# Simulation Period
#---------------------------------------------------

sim_duration  = 20.0

timespan_min  = restoration_time + Δt1 + Δt2

timespan      =
    timespan_min < sim_duration ? sim_duration : timespan_min

time_start    = 0.0

time_final    = timespan

tspan         = (0.0, timespan)

sim_timespan  = (0.0, timespan)

plot_timespan = (0.0, timespan)

#---------------------------------------------------

#---------------------------------------------------
# Reading network data
#---------------------------------------------------

(;gens_govs_avrs_states_syms,
 gens_govs_avrs_types) =
     get_gens_govs_avrs_states_syms_by_json(
            net_data_by_components_file;
            components_libs_dir =
                components_libs_dir )

#---------------------------------------------------

(;plant_generators_data_from_json,
 plant_loads_data_from_json,
 plant_transmission_data_from_json,
 edge_data_from_json,
 shunt_data_from_json,
 baseMVA_data_from_json,
 gencost_data_from_json) =
    NamedTupleTools.select(
        get_net_data_by_components_from_json_file(
            net_data_by_components_file;
            in_components_type_sym =
                false ),
        (:plant_generators_data_from_json,
         :plant_loads_data_from_json,
         :plant_transmission_data_from_json,
         :edge_data_from_json,
         :shunt_data_from_json,
         :baseMVA_data_from_json,
         :gencost_data_from_json))

baseMVA = baseMVA_data_from_json

#----------------------------------------
# Dynamic components functions
#----------------------------------------

(;comps_callback_paras_funs,
 comps_init_funs,
 comps_output_funs,
 ode_comps_dyn_funs,
 dae_comps_dyn_funs,
 comps_dyn_funs) =
     get_dynamic_comps_init_out_dyn_callback_funcs(
    gens_govs_avrs_types)

#----------------------------------------
# Generators dynamic and static parameters
#----------------------------------------

# gens_para_sequence_order =
#     (:components_data, :gen)

gens_generic_sequence_order =
    (:components_data, :gen)


ode_gens_para_selections  =
    (:H, :D,
     :X_d, :X_q,                  
     :X_d_dash, :X_q_dash,
     :T_d_dash, :T_q_dash, :Sn )

ode_gens_generic_selections =
    (:H, :D,
     :ra, :xℓ,
     :X_d, :X_q,
     :X_d_dash,  :X_q_dash,
     :X_d_2dash, :X_q_2dash,
     :T_d_dash,  :T_q_dash, :Sn )

opf_gens_generic_selections =
   (
    :Sn,:vh,
    :P, :Q,
    :Pmin, :Pmax,
    :Qmin, :Qmax,
    :vmin, :vmax )

govs_and_avrs_sequence_order =
    ( :components_data,)

govs_and_avrs_selections =
    ( :gov, :avr )

#----------------------------------------

ode_gens_generic_para =
     get_ode_gens_generic_para(
         plant_generators_data_from_json;
         sequence_order =
             gens_generic_sequence_order,
         selections =
             ode_gens_generic_selections)

#----------------------------------------

(;generic_gens_para,
 generic_govs_para,
 generic_avrs_para) =
     get_generic_gens_avr_gov_para(
         plant_generators_data_from_json;
         gens_sequence_order =
             gens_generic_sequence_order,
        gens_selections =
            ode_gens_generic_selections,
        govs_and_avrs_sequence_order =
            govs_and_avrs_sequence_order,
        govs_and_avrs_selections =
            govs_and_avrs_selections)

#------------------------------------------------

 pf_generic_gens_para =
    NamedTupleTools.select(
        get_selected_vec_nt_to_vec_vec(
            generic_gens_para,
            nothing;
            selections =
                (:ra, :X_d, :X_q,     
                 :X_d_dash, :X_q_dash, :Sn),
            vec_datatype = Float64 ),
        (:ra, :X_d, :X_q,     
         :X_d_dash, :X_q_dash, :Sn) )

#------------------------------------------------

opf_generic_each_gen_para =
    get_components_properties_by_json(
        plant_generators_data_from_json;
        sequence_order =
            gens_generic_sequence_order,
         selections =
              opf_gens_generic_selections )


"To make sure Vector{NamedTuple} is returned
 instead of Vector{Any}"
opf_generic_each_gen_para =
    NamedTuple[
        item for item in
            opf_generic_each_gen_para]

#----------------------------------------

 opf_generic_gens_para =
    NamedTupleTools.select(
        get_selected_vec_nt_to_vec_vec(
            opf_generic_each_gen_para,
            nothing;
            selections =
                (:vh,  :P, :Q,
                 :Pmin, :Pmax, :Qmin, :Qmax,
                 :vmin, :vmax, :Sn ),
            vec_datatype = Float64 ),
        (:vh, :P, :Q,
         :Pmin, :Pmax, :Qmin, :Qmax,
         :vmin, :vmax, :Sn ) )

#---------------------------------------------------
# basic edges parameters parameters
#---------------------------------------------------

(;edges_orientation,
 edges_Ybr_cal,
 Ybr_cal_and_edges_orientation,
 Ynet_wt_nodes_idx_wt_adjacent_nodes) =
     NamedTupleTools.select(
         get_transmission_network_parameters_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json,
             edge_data_from_json,
             shunt_data_from_json;
             baseMVA =
                 baseMVA,
             basekV =
                 basekV,
             use_pu_in_PQ =
                 use_pu_in_PQ,
             line_data_in_pu =
                 line_data_in_pu ),
         (:edges_orientation,
          :edges_Ybr_cal,
          :Ybr_cal_and_edges_orientation,
          :Ynet_wt_nodes_idx_wt_adjacent_nodes))

#----------------------------------------

edges_ftbus_and_generic_data =
      get_edges_ftbus_and_generic_data_by_json(
         edge_data_from_json )

#----------------------------------------

edges_generic_data =
     get_edges_generic_data_by_json(
         edge_data_from_json )

#----------------------------------------
# basic nodes parameters parameters
#----------------------------------------

pf_PQ_param =
    get_pf_PQ_param_by_json(
        plant_generators_data_from_json,
        plant_loads_data_from_json,
        plant_transmission_data_from_json;
        baseMVA = baseMVA,
        use_pu_in_PQ =
            use_pu_in_PQ )
#----------------------------------------

sta_pf_PQ_para =
    get_pf_PQ_param_by_json(
        plant_generators_data_from_json,
        plant_loads_data_from_json,
        plant_transmission_data_from_json;
        baseMVA =
            baseMVA,
        use_pu_in_PQ =
            use_pu_in_PQ)

#----------------------------------------

gens_vh_slack_θh_para =
    get_gens_vh_slack_θh_para_by_json(
        plant_generators_data_from_json )

#----------------------------------------

sta_pf_vars_and_paras_idx =
    get_sta_pf_vars_and_paras_idx_by_json(
        plant_generators_data_from_json,
        plant_loads_data_from_json,
        plant_transmission_data_from_json )

#----------------------------------------

pf_sta_ΔPQ_mismatch_parameters =
    get_pf_sta_ΔPQ_mismatch_parameters_by_json(
        plant_generators_data_from_json,
        plant_loads_data_from_json,
        plant_transmission_data_from_json,
        edge_data_from_json,
        shunt_data_from_json;
        baseMVA = baseMVA,
        basekV = 1.0,
        use_pu_in_PQ = use_pu_in_PQ,
        line_data_in_pu = line_data_in_pu)

#----------------------------------------
#----------------------------------------
# Some basic indices, labels and symbols
#----------------------------------------
#----------------------------------------


#----------------------------------------
# Nodes type indices and translation dict
#----------------------------------------

net_nodes_type_idxs =
    get_net_nodes_type_idxs_by_json(
        plant_generators_data_from_json,
        plant_loads_data_from_json,
        plant_transmission_data_from_json )

loc_load_exist =
    getproperty(
        net_nodes_type_idxs,
        :loc_load_exist)

#----------------------------------------

dyn_pf_fun_kwd_net_idxs =
    NamedTupleTools.select(
        net_nodes_type_idxs,
        (:slack_gens_nodes_idx,
         :non_slack_gens_nodes_idx,
         :gens_nodes_idx,
         :non_gens_nodes_idx,
         :gens_with_loc_load_idx,
         :gens_nodes_with_loc_loads_idx,
         :all_nodes_idx))

#----------------------------------------

dyn_pf_fun_kwd_n2s_idxs =
    NamedTupleTools.select(
        get_dict_net_streamlined_idx_by_nodes_type_idxs(
            net_nodes_type_idxs ),
        (:n2s_slack_gens_idx,
         :n2s_non_slack_gens_idx,
         :n2s_gens_idx,
         :n2s_non_gens_idx,
         :n2s_gens_with_loc_load_idxs,
         :n2s_all_nodes_idx))

(gens_nodes_idx,
 non_gens_nodes_idx,
 gens_nodes_with_loc_loads_idx,
 all_nodes_idx) =
     NamedTupleTools.select(
         dyn_pf_fun_kwd_net_idxs,
         (:gens_nodes_idx,
          :non_gens_nodes_idx,
          :gens_nodes_with_loc_loads_idx,
          :all_nodes_idx))

(;n2s_gens_idx,
n2s_non_gens_idx,
n2s_gens_with_loc_load_idxs,
n2s_all_nodes_idx ) =
    NamedTupleTools.select(
        dyn_pf_fun_kwd_n2s_idxs,
        (:n2s_gens_idx,
         :n2s_non_gens_idx,
         :n2s_gens_with_loc_load_idxs,
         :n2s_all_nodes_idx))

#----------------------------------------
# States and algebraic variables labels
#----------------------------------------

plants_states_syms =
    get_plants_states_syms(
        gens_govs_avrs_states_syms)

generic_state_sym = 
    get_generic_state_sym(
        gens_govs_avrs_states_syms,
        gens_nodes_idx;
        label_prefix = "bus")

state_labels = 
    get_state_labels(
        gens_govs_avrs_states_syms,
        gens_nodes_idx;
        label_prefix = "bus",
        plants_states_by_per_comp = false,
        plants_states_by_per_plant = true )


algebraic_vars_labels = 
    get_algebraic_vars_labels(
        dyn_pf_fun_kwd_net_idxs;
        label_prefix = "bus" )


network_vars_labels = 
    get_network_vars_labels(
        gens_govs_avrs_states_syms,
        dyn_pf_fun_kwd_net_idxs;
        label_prefix = "bus",
        plants_states_by_per_comp = false,
        plants_states_by_per_plant = true )

state_vars_idx = 
   get_state_vars_idx(
       gens_govs_avrs_states_syms)

(;nodes_names,
 gens_nodes_names,
 non_gens_nodes_names) =
     get_model_nodes_types_names(
          dyn_pf_fun_kwd_net_idxs,
          dyn_pf_fun_kwd_n2s_idxs )

(;state_vars_idx,
 vec_comp_states_Idx,
 plants_states_syms,
 generic_state_sym,
 state_labels,
 algebraic_vars_labels,
 network_vars_labels) =
     get_plants_states_syms_and_labels(
         gens_govs_avrs_states_syms,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_n2s_idxs)

(;state_vars_idx,
 vec_comp_states_Idx,
 plants_states_syms,
 generic_state_sym,
 state_labels,
 algebraic_vars_labels,
 network_vars_labels,
 
 model_syms,
 nodes_names,
 gens_nodes_names,
 non_gens_nodes_names,
 SM_gens_nodes_names,
 SC_gens_nodes_names) = 
    get_plants_states_syms_wt_labels_wt_names(
        gens_govs_avrs_states_syms,
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs)

(;state_labels,
 algebraic_vars_labels,
 network_vars_labels) =
     get_generic_network_vars_labels(
         plants_states_syms,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_n2s_idxs
         ;label_prefix = "bus",
         plants_states_by_per_comp = false,
         plants_states_by_per_plant = true)

model_syms = 
get_model_syms(
    state_labels,
    dyn_pf_fun_kwd_net_idxs)


#----------------------------------------
# Mass matrix and bool dae vars
#----------------------------------------

model_mass_matrix =
    DAE_MassMatrix(
        length(state_labels),
        length(algebraic_vars_labels) )

model_bool_dae_vars =
    DAE_BoolVector(
        length(state_labels),
        length(algebraic_vars_labels) )

ode_gens_mass_matrix =
    DAE_MassMatrix(
        length(state_labels),
        0 )

ode_gens_bool_dae_vars =
    DAE_BoolVector(
        length(state_labels),
        0 )

#----------------------------------------

(;model_mass_matrix,
 model_bool_dae_vars,
 ode_gens_mass_matrix,
 ode_gens_bool_dae_vars) =
     get_mass_matrix_and_bool_dae_vars(
    state_labels,
    algebraic_vars_labels)

#----------------------------------------
# Parameter indices
#----------------------------------------

Png_Qng_Pll_Qll_Idx = 
get_generic_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)

Pg_Qg_Png_Qng_Pll_Qll_Idx = 
 get_generic_Pg_Qg_Png_Qng_Pll_Qll_Idx(
     dyn_pf_fun_kwd_net_idxs)

#----------------------------------------

scale_Pg_Qg_Png_Qng_Pll_Qll_Idx = 
 get_generic_scale_Pg_Qg_Png_Qng_Pll_Qll_Idx(
     dyn_pf_fun_kwd_net_idxs)

scale_Pg_Png_Qng_Idx =
get_generic_scale_Pg_Png_Qng_Idx(
    dyn_pf_fun_kwd_net_idxs)

#----------------------------------------

Png_Qng_Pll_Qll_Idx = 
get_generic_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)

#----------------------------------------

Pg_Png_Qng_Idx = 
get_generic_Pg_Png_Qng_Idx(
    dyn_pf_fun_kwd_net_idxs)

#----------------------------------------

dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx =
get_dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)


dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx =
get_dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)


dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx =
get_dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)


dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx =
get_dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)


dyn_vh_id_iq_V_ref_Tm_Idx = 
get_dyn_vh_id_iq_V_ref_Tm_Idx(
    gens_nodes_idx;
    reverse_idx = false)


dyn_Tm_V_ref_id_iq_vh_Idx =
get_dyn_vh_id_iq_V_ref_Tm_Idx(
    gens_nodes_idx;
    reverse_idx = true )


dyn_V_ref_Tm_id_iq_vh_Idx  =
get_dyn_V_ref_Tm_vh_id_iq_Idx(
    gens_nodes_idx)


dyn_vh_id_iq_ωref0_vref0_porder0_Idx =
get_dyn_vh_id_iq_ωref0_vref0_porder0_Idx(
    gens_nodes_idx;
    reverse_idx = false )


 ωref0_vref0_porder0_id_iq_vh_Idx = 
 get_ωref0_vref0_porder0_id_iq_vh_Idx(
     gens_nodes_idx)


dyn_ωref0_vref0_porder0_id_iq_vh_Idx =
get_dyn_ωref0_vref0_porder0_id_iq_vh_Idx(
    gens_nodes_idx)


id_iq_pg_vh_Idx = 
get_id_iq_pg_vh_Idx(
    gens_nodes_idx)

#--------------------------------------

dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx =
    get_dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)


dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx =
    get_dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)

#----------------------------------------
# States and algebraic variable  indices
#----------------------------------------


(;state_vars_and_i_dq_Idx_in_state,
 state_vars_and_i_dq_wt_fault_Idx_in_state,
 state_algebraic_vars_Idx_in_state,
 state_algebraic_vars_wt_fault_Idx_in_state) =
     get_state_and_algebraic_vars_Idx_in_state(
         state_labels,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_n2s_idxs;
         no_lines_fault = 1 )

(;state_vars_and_i_dq_Idx_in_state,
 state_vars_and_i_dq_wt_fault_Idx_in_state,
 state_algebraic_vars_Idx_in_state,
 state_algebraic_vars_wt_fault_Idx_in_state) =
     get_state_and_algebraic_vars_Idx_in_state(
         state_labels,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_n2s_idxs;
         no_lines_fault = 1 )

gens_state_vars_idx_in_state =
get_gens_state_vars_idx_in_state(
    network_vars_labels,
    # all_nodes_idx,
    dyn_pf_fun_kwd_net_idxs,
    n2s_all_nodes_idx;
    selected_gens_state_vars_syms =
        (:δ, :ed_dash, :eq_dash) )


state_vars_and_i_dq_Idx_in_state =
get_state_vars_and_i_dq_Idx_in_state(
    generic_state_sym, #generic_state_labels,
    gens_nodes_idx,
    all_nodes_idx )

state_vars_and_i_dq_wt_fault_Idx_in_state =
get_state_vars_and_i_dq_wt_fault_Idx_in_state(
    generic_state_sym, # generic_state_labels,
    gens_nodes_idx,
    all_nodes_idx;
    no_lines_fault = 1)

#----------------------------------------

state_algebraic_vars_Idx_in_state =
get_state_algebraic_vars_Idx_in_state(
    generic_state_sym, # generic_state_labels,
    gens_nodes_idx,
    all_nodes_idx )


state_algebraic_vars_wt_fault_Idx_in_state = 
get_state_algebraic_vars_wt_fault_Idx_in_state(
    generic_state_sym, # generic_state_labels,
    gens_nodes_idx,
    all_nodes_idx;
    no_lines_fault = 1)

#----------------------------------------
#  Algebraic variables  indices
#----------------------------------------

## state_vars_idx

pf_red_vh_θh_wt_slack_value_Idx = 
get_generic_red_vh_θh_wt_slack_value_Idx(
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)


pf_vh_θh_idx_and_idx2Idx = 
get_pf_vh_θh_idx_and_idx2Idx(
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)


# get_generic_flat_vh_flat_θh_Idx(
#     gens_nodes_idx,
#     all_nodes_idx)

dyn_pf_flat_vh_flat_θh_Idx  =
get_generic_flat_vh_flat_θh_Idx(
    all_nodes_idx)

#----------------------------------------

dyn_pf_flat_vh_flat_θh_id_iq_Idx =
get_generic_flat_vh_flat_θh_id_iq_Idx(
    gens_nodes_idx,
    all_nodes_idx)

#----------------------------------------

dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx =
get_generic_flat_vh_flat_θh_wt_slack_value_Idx(
    all_nodes_idx)

#----------------------------------------

dyn_pf_vh_θh_id_iq_vhf_θhf_Idx =
 get_generic_vh_θh_id_iq_vhf_θhf_Idx(
    gens_nodes_idx,
    all_nodes_idx;
     no_lines_fault = 1)


dyn_pf_vh_vhf_θh_θhf_id_iq_Idx =
get_generic_vh_vhf_θh_θhf_id_iq_Idx(
    gens_nodes_idx,
    all_nodes_idx;
    no_lines_fault = 1)

dyn_pf_vh_vhf_Idx =
 get_generic_vh_vhf_Idx(
    all_nodes_idx;
     no_lines_fault = 1)

dyn_pf_θh_θhf_Idx = 
get_generic_θh_θhf_Idx(
    all_nodes_idx;
    no_lines_fault = 1)

#----------------------------------------


get_generic_red_vh_θh_wt_slack_value_Idx(
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)


dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx =
get_generic_flat_vh_flat_θh_wt_slack_value_Idx(
    all_nodes_idx)


(;dyn_slack_value_Idxs,) =
    NamedTupleTools.select(
        dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
        (:dyn_slack_value_Idxs, ))

pf_vh_θh_idx_and_idx2Idx =
get_pf_vh_θh_idx_and_idx2Idx(
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)

#----------------------------------------
#----------------------------------------

get_model_states_comp_idxs_in_Idx(
        network_vars_labels,
        all_nodes_idx,
        n2s_all_nodes_idx;
    vars =
        [:δ, :ω, :ed_dash, :eq_dash] )



#----------------------------------------
# faults
#----------------------------------------

list_faulted_line_a_b_orientation =
    edges_orientation[
        list_edges_to_have_fault  ] 

#----------------------------------------

no_lines_fault =
    length(list_faulted_line_a_b_orientation)

no_current_lines_fault =
    no_lines_fault - length(
        clear_fault_selection_list)


no_cleared_lines_fault =
    length(clear_fault_selection_list)

#----------------------------------------

# (;faulty_Ynet,
# faulty_nodes_idx_with_adjacent_nodes_idx,

# faulty_all_nodes_idx,
# n2s_faulty_all_nodes_idx,

# fault_nodes_idx,
# n2s_fault_nodes_idx,            

# list_Ya_nkf,
# list_Ynkf_b,

# list_faulty_line_Yl,
# list_healthy_lines_Yl,

# list_node_b_idx_in_a_node_row,
# list_node_a_idx_in_b_node_row,

# list_faulted_line_a_b_orientation,
# list_fault_point_from_node_a,
# list_fault_resistance,
# list_no_line_circuit)

(Ynet,
nodes_idx_with_adjacent_nodes_idx) =
    Ynet_wt_nodes_idx_wt_adjacent_nodes

on_fault_net_para =
     make_lines_faults_data_set(
         Ynet,
         nodes_idx_with_adjacent_nodes_idx,

         all_nodes_idx,
         n2s_all_nodes_idx,

         list_faulted_line_a_b_orientation ,
         list_fault_point_from_node_a,
         list_fault_resistance,
         list_no_line_circuit )

(;fault_nodes_idx,) =
    NamedTupleTools.select(
        on_fault_net_para,
        (:fault_nodes_idx,))

#--------------------------------------

# (;pre_clear_fault_Ynet,
# pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,

# pre_clear_fault_all_nodes_idx,
# n2s_pre_clear_fault_all_nodes_idx,

# pre_clear_fault_nodes_idx,
# n2s_pre_clear_fault_nodes_idx,

# pre_clear_list_Ya_nkf, 
# pre_clear_list_Ynkf_b, 

# pre_clear_list_faulty_line_Yl, 
# pre_clear_list_healthy_lines_Yl, 

# pre_clear_list_node_b_idx_in_a_node_row, 
# pre_clear_list_node_a_idx_in_b_node_row, 

# pre_clear_list_faulted_line_a_b_orientation, 
# pre_clear_list_fault_point_from_node_a,

# pre_clear_list_fault_resistance,
# pre_clear_list_no_line_circuit,

# post_clear_fault_Ynet,    
# post_clear_fault_nodes_idx_with_adjacent_nodes_idx,

# post_clear_fault_all_nodes_idx,    
# n2s_post_clear_fault_all_nodes_idx,

# post_clear_fault_nodes_idx,    
# n2s_post_clear_fault_nodes_idx,

# faulty_Ynet,
# faulty_nodes_idx_with_adjacent_nodes_idx,

# faulty_all_nodes_idx,
# n2s_faulty_all_nodes_idx,

# fault_nodes_idx,
# n2s_fault_nodes_idx,

# list_faulted_line_a_b_orientation,
# list_fault_point_from_node_a,
# list_fault_resistance,
# list_no_line_circuit,

# list_Ya_nkf,
# list_Ynkf_b,

# list_faulty_line_Yl,
# list_healthy_lines_Yl,

# list_node_b_idx_in_a_node_row,
# list_node_a_idx_in_b_node_row)

cleared_selected_lines_faults_net_para =
    get_cleared_selected_lines_faults_data_set(
        clear_fault_selection_list;
        deepcopy(on_fault_net_para)...)

(;pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
post_clear_fault_nodes_idx_with_adjacent_nodes_idx) =
 NamedTupleTools.select(
     cleared_selected_lines_faults_net_para,
     (
      :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
      :post_clear_fault_nodes_idx_with_adjacent_nodes_idx))

#----------------------------------------
# Callbacks
#----------------------------------------

plants_cb_paras_kwd_para =
    (generic_gens_para ,
     generic_avrs_para,
     generic_govs_para,
     comps_callback_paras_funs )

#----------------------------------------

"""
generic_model_callback_paras =
    plants_generic_model_callback_paras_func(
        state_vars_idx,
        plants_states_syms;
        kwd_para =
            plants_cb_paras_kwd_para )

plants_cb_paras_switches =
    getproperty(
        generic_model_callback_paras,
        :plants_cb_paras_switches)


avrs_govs_cb_sw =
    getproperty(
        generic_model_callback_paras,
        :plants_avr_gov_cb_para_sw_in_plant)

avrs_govs_cb_sw_Idx =
    getproperty(
        generic_model_callback_paras,
        :plants_avr_gov_cb_para_sw_idx_in_plant )

cb = cb_fun_make_state_callbacks(
    generic_model_callback_paras)
"""

#----------------------------------------

(plants_cb_paras_switches,
 list_selected_plants_state_event_cb_paras,
 list_selected_plants_state_affect_cb_paras,

 avrs_govs_cb_sw,
 avrs_govs_cb_sw_Idx ) =
     NamedTupleTools.select(
         plants_generic_model_callback_paras_func(
             state_vars_idx,
             plants_states_syms;
             kwd_para =
                 plants_cb_paras_kwd_para ) ,
         (:plants_cb_paras_switches,
          :list_selected_plants_state_event_cb_paras,
          :list_selected_plants_state_affect_cb_paras,

          :plants_avr_gov_cb_para_sw_in_plant,
          :plants_avr_gov_cb_para_sw_idx_in_plant))

cb = cb_fun_make_state_callbacks(
    list_selected_plants_state_event_cb_paras,
    list_selected_plants_state_affect_cb_paras )


#----------------------------------------
#########################################
#----------------------------------------

#----------------------------------------
# Composition and aggregations
#----------------------------------------

(pf_kw_para,
 red_types_Idxs_etc,
 pf_PQ_param) =
    NamedTupleTools.select(
        pf_sta_ΔPQ_mismatch_parameters,
        (:pf_kw_para,
         :red_types_Idxs_etc,
         :pf_PQ_param) )

(red_vh_Idxs,
 red_θh_Idxs) =
    NamedTupleTools.select(
        red_types_Idxs_etc,
        (:red_vh_Idxs,
         :red_θh_Idxs) )

#----------------------------------------

sta_red_vh_θh_0 =
   [ ones(length(red_vh_Idxs));
     zeros(length(red_θh_Idxs))]

#----------------------------------------
# Powerflow func and prob
#----------------------------------------

kwd_sta_sta_ΔPQ_sol_by_json =
   (;
    pf_alg,
    pf_kw_para,
    red_vh_Idxs,
    red_θh_Idxs,
    sta_red_vh_θh_0) 


#----------------------------------------
# Powerflow func and prob
#----------------------------------------

pf_sol =
    get_pf_sta_ΔPQ_mismatch_sol_by_generic(
        pf_PQ_param;
        kwd_para =
            kwd_sta_sta_ΔPQ_sol_by_json )

#----------------------------------------
# Results    
#----------------------------------------

generic_red_sol_kwd_para =
    (;Ybr_cal_and_edges_orientation,
      Ynet_wt_nodes_idx_wt_adjacent_nodes,
      sta_pf_PQ_para,
      ode_gens_generic_para,
      pf_kw_para) 

# generic_dyn_sol_kwd_para =
#    (;loc_load_exist,
#     sta_pf_PQ_para,
#     # ode_gens_generic_para,
#     # dyn_pf_flat_vh_flat_θh_id_iq_Idx,
#     dyn_pf_fun_kwd_n2s_idxs,
#     dyn_pf_fun_kwd_net_idxs)


generic_results_pf_sta_red_sol =
    get_generic_results_pf_sta_red_sol_u(
        pf_sol;
        generic_red_sol_kwd_para =
            generic_red_sol_kwd_para,
        baseMVA =
            baseMVA,
        basekV =
            1.0 )

#----------------------------------------
    
(pf_P_gens,
 pf_Q_gens,
 vh,
 θh,
 gens_vh,
 gens_θh) =
    NamedTupleTools.select(
        generic_results_pf_sta_red_sol,
        (:pf_P_gens,
         :pf_Q_gens,
         :vh,
         :θh,
         :gens_vh,
         :gens_θh ) )

#----------------------------------------
# Init
#----------------------------------------

plants_init_kwd_para =
    (generic_gens_para ,
     generic_avrs_para,
     generic_govs_para ,
     comps_init_funs )

plants_generic_states_init_wt_ref =
    plants_generic_model_init_func(
        gens_vh,
        gens_θh,
        pf_P_gens,
        pf_Q_gens,
        ωs;
        kwd_para =
            plants_init_kwd_para )

(plants_states_init,
 plants_refs ) =
     NamedTupleTools.select(
         plants_generic_states_init_wt_ref,
         (:plants_states_init,
          :plants_refs))


( nt_vec_per_paras,
  vec_vec_per_paras ) =
      get_nt_vec_wt_vec_vec_per_paras(
    plants_refs ;
    nt_syms =
        (:ωs,
         :ω_ref,
         :v_ref,
         :p_order,
         :i_d,
         :i_q,
         :gen_vh ) )

(ω_ref,
 v_ref,
 p_order,
 gens_i_d,
 gens_i_q ) =
     NamedTupleTools.select(
         nt_vec_per_paras, (
             :ω_ref,
             :v_ref,
             :p_order,
             :i_d,
             :i_q))

#----------------------------------------
# System model init
#----------------------------------------

u0_model_states_init =
    Float64[plants_states_init;
            vh;
            θh;
            gens_i_d;
            gens_i_q]

du0_model_states_init =
    zeros(length(u0_model_states_init))

res = similar(u0_model_states_init)

#----------------------------------------
# System model para
#----------------------------------------

(P_non_gens,
 Q_non_gens, 
 P_g_loc_load,
 Q_g_loc_load) =
    NamedTupleTools.select(
        sta_pf_PQ_para,
        (:P_non_gens,
         :Q_non_gens,
         :P_g_loc_load,
         :Q_g_loc_load ) )

generic_model_dynamics_para =
    Float64[ω_ref;
            v_ref;
            p_order;
            P_non_gens;
            Q_non_gens;
            P_g_loc_load;
            Q_g_loc_load]


ω_ref_v_ref_p_order_Png_Qng_Pll_Qll =
    Float64[ω_ref;
            v_ref;
            p_order;
            P_non_gens;
            Q_non_gens;
            P_g_loc_load;
            Q_g_loc_load]        

ωref0_vref0_porder0_id_iq_vh =
    [ω_ref;
     v_ref;
     p_order;
     gens_i_d;
     gens_i_q;
     gens_vh]


flat_vh_flat_θh_id_iq_u0 =
    [vh;
     θh;
     gens_i_d;
     gens_i_q ]

dyn_pf_solver =
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol )

algebraic_generic_model_kwd_para =
    (;loc_load_exist,

     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes )

algebraic_generic_model_sol_kwd_para =
    (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_kwd_para
     )


algebraic_generic_model_wt_fault_kwd_para =
    (;loc_load_exist,

     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

     on_fault_net_para,

     clear_fault_selection_list,

     no_lines_fault,
     no_cleared_lines_fault,

     list_fault_resistance,

     with_faults,

     cleared_selected_lines_faults_net_para )

algebraic_generic_model_wt_fault_sol_kwd_para =
    (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_wt_fault_kwd_para,
     cleared_selected_lines_faults_net_para
     )

#----------------------------------------
#----------------------------------------
# generic model_dynamics_kwd_para
#----------------------------------------
#----------------------------------------

#----------------------------------------
# ODE system dyamanics kwd parameters
#----------------------------------------    


plants_kwd_para =
    (;state_vars_idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,

     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,

     vec_comp_states_Idx,

     avrs_govs_cb_sw,
     avrs_govs_cb_sw_Idx,

     comps_dyn_funs,
     comps_output_funs,

    ωs) 

# plants_kwd_para =
#     deepcopy(ode_plants_kwd_para)


ode_plants_kwd_para =
    (;state_vars_idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,

     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,

     vec_comp_states_Idx,

     avrs_govs_cb_sw,
     avrs_govs_cb_sw_Idx,

     ode_comps_dyn_funs,
     comps_output_funs,
     ωs)


ode_generic_model_dynamics_kwd_para =
    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     # pf_generic_gens_para,

     ode_plants_kwd_para,
     plants_kwd_para,

     algebraic_generic_model_sol_kwd_para )


ode_generic_model_kwd_para =
    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     ode_plants_kwd_para,

     algebraic_generic_model_sol_kwd_para )


#----------------------------------------
# DAE system dyamanics kwd parameters
#----------------------------------------    

dae_plants_kwd_para =
    (;state_vars_idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,

     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,

     vec_comp_states_Idx,

     avrs_govs_cb_sw,
     avrs_govs_cb_sw_Idx,

     dae_comps_dyn_funs,
     comps_output_funs,
     ωs)



dae_generic_model_dynamics_kwd_para =
    (;
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dae_plants_kwd_para )


dae_generic_model_kwd_para =
    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_net_idxs,

     dae_plants_kwd_para,
     
     algebraic_generic_model_kwd_para )


#----------------------------------------
# ODE mass matrix system dyamanics kwd parameters
#----------------------------------------    


mm_generic_model_kwd_para = 
    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, #

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs, #
     dyn_pf_fun_kwd_net_idxs, #

     ode_plants_kwd_para,

     algebraic_generic_model_kwd_para ) 


mm_generic_model_dynamics_kwd_para = 
    (;
     ωs,
     loc_load_exist,
     state_vars_idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, #

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs, #
     dyn_pf_fun_kwd_net_idxs, #

     ode_plants_kwd_para,

     algebraic_generic_model_sol_kwd_para ) 


#----------------------------------------
# System dyamanics kwd parameters with faults
#----------------------------------------    

generic_system_dynamics_wt_fault_kwd_para =
    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,

     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     ode_plants_kwd_para,
     dae_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,
     algebraic_generic_model_wt_fault_sol_kwd_para,

     algebraic_generic_model_kwd_para,
     algebraic_generic_model_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     cleared_selected_lines_faults_net_para,

     with_faults,
     generic_results_pf_sta_red_sol )

pre_post_generic_system_dynamics_wt_fault_kwd_para =
    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, 

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     # dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 
     # pf_generic_gens_para, 
     # Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     dae_plants_kwd_para,

     algebraic_generic_model_kwd_para,

     # cleared_selected_lines_faults_net_para,
     nodes_idx_with_adjacent_nodes_idx,
     pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
     post_clear_fault_nodes_idx_with_adjacent_nodes_idx)

#---------------------------------------------------
#---------------------------------------------------
