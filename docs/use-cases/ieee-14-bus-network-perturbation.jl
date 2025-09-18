#=
# [Example of Network Perturbation with IEEE 14 Bus Test System](@id ieee-14-bus-network-perturbation)
=#


#---------------------------------------------------

using Revise

# using Pkg

using ePowerSim

#---------------------------------------------------

using SciMLBase

using SciMLNLSolve: NLSolveJL

using NLSolvers

using NLsolve: nlsolve, converged, OnceDifferentiable

using NLsolve

using NonlinearSolve

using NonlinearSolve: TrustRegion

#---------------------------------------------------

using DifferentialEquations

using OrdinaryDiffEq, Sundials, ODEInterfaceDiffEq

using FiniteDiff, LinearSolve

using ForwardDiff, Zygote

using DiffRules

using PreallocationTools

#---------------------------------------------------

using LinearAlgebra, GenericSchur, Arblib

#---------------------------------------------------

using OrderedCollections: OrderedDict

using Permutations

using SparseArrays, StaticArrays, ComponentArrays

using DataFrames, DataFramesMeta

using JSONTables, JSON, JSON3

using Query, CSV, Tables, XLSX

using JLD2

#---------------------------------------------------

using Optimization

using OptimizationOptimJL

using JuMP

using HiGHS

import Clarabel

import Ipopt

import PATHSolver

#---------------------------------------------------

using Graphs

#---------------------------------------------------

using StatsBase

using Plots

import StatsPlots

using Latexify

using LatexPrint

#---------------------------------------------------

using BenchmarkTools, BenchmarkPlots

#---------------------------------------------------

using NamedTupleTools

using Accessors, AccessorsExtra 


#---------------------------------------------------
# global settings
#---------------------------------------------------

freq = 60

Ωb = 2 * pi * freq

ωs = Ωb 

ω_ref0 = ωs

basekV = 1.0

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

"""

sauer_net_wt_avr_string    = "net-static-data-avr-sauer-"

rtds_net_wt_avr_string     = "net-static-data-avr-rtds-"

sauer_gov_string           = "gov-sauer"

ieee_tgov2sauer_gov_string = "gov-ieee-tgov2sauer"

ieee_tgov1_gov_string      = "gov-ieee-tgov1"

#---------------------------------------------------

# net_wt_avr_string = rtds_net_wt_avr_string
# gov_string        = ieee_tgov1_gov_string

net_wt_avr_string = sauer_net_wt_avr_string
gov_string        = sauer_gov_string

dynamic_net_data_by_components_file =
    "$(net_wt_avr_string)"*
    "$(gov_string)" *
    ".json"

json_net_data_by_components_file =
    dynamic_net_data_by_components_file

"""

json_net_data_by_components_file =
    "net-static-data-avr-sauer-gov-sauer.json"

#---------------------------------------------------

# case_name = "case9"
case_name = "case14"

# This ensures a lowercase name

case_name = lowercase(case_name) 

#---------------------------------------------------
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
#---------------------------------------------------

# sim_type  = "$(case_name)-"*"network-pertubation-"*
#     "-$(net_wt_avr_string)-$(gov_string)"

sim_type  = "$(case_name)-"*"network-pertubation"


# sim_type  = "sim-network-perturbation"

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

basekV = 1.0

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

list_edges_to_have_fault     = [ 8 ] # [ 4 ]

clear_fault_selection_list   = [1]

#---------------------------------------------------

on_fault_time = 10.617 # 10.0

Δt_clear_time = 0.2

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
## ntuple_status_steady_state_data
#---------------------------------------------------

ntuple_status_steady_state_data =
    get_ntuple_status_steady_state_data(
    ;with_faults = with_faults,
    net_data_by_components_file =
        net_data_by_components_file,
    components_libs_dir =
        components_libs_dir,
    
        timespan =
            timespan,
        on_fault_time =
            on_fault_time,
        clear_fault_time =
            clear_fault_time,
    
        list_fault_point_from_node_a =
            list_fault_point_from_node_a,
        list_fault_resistance =
            list_fault_resistance,
        list_no_line_circuit =
            list_no_line_circuit,

        list_edges_to_have_fault =
            list_edges_to_have_fault,
        clear_fault_selection_list =
            clear_fault_selection_list,
    
        basekV =
            basekV,    
        use_pu_in_PQ =
            use_pu_in_PQ,
        line_data_in_pu =
            line_data_in_pu)

#---------------------------------------------------


(;state_labels,
 algebraic_vars_labels,

 dyn_pf_fun_kwd_n2s_idxs,
 dyn_pf_fun_kwd_net_idxs,

 system_fault_status,
 generic_system_dynamics_wt_fault_kwd_para,
 on_fault_net_para,
 cleared_selected_lines_faults_net_para,

 ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

 # model_bool_dae_vars_wt_fault,
 # model_syms_wt_fault,         
 # u0_model_states_init_wt_fault,

 model_bool_dae_vars,     
 model_syms,
 u0_model_states_init,
 model_mass_matrix,     

 cb_states,
 plants_cb_paras_switches,

 nodes_names,
 gens_nodes_names,
 non_gens_nodes_names,
 SM_gens_nodes_names,
 SC_gens_nodes_names,

 ωref0_vref0_porder0_id_iq_vh_Idx,
 dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

 dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

 Ybr_cal_and_edges_orientation,
 Ynet_wt_nodes_idx_wt_adjacent_nodes) =
    NamedTupleTools.select(
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :pre_fault_state),
            :static_prefault_paras),        
        (:state_labels,
         :algebraic_vars_labels,
         :dyn_pf_fun_kwd_n2s_idxs,
         :dyn_pf_fun_kwd_net_idxs,

         :system_fault_status,
         :generic_system_dynamics_wt_fault_kwd_para,
         :on_fault_net_para,
         :cleared_selected_lines_faults_net_para,

         :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

         # :model_bool_dae_vars_wt_fault,
         # :model_syms_wt_fault,         
         # :u0_model_states_init_wt_fault,

         :model_bool_dae_vars,     
         :model_syms,
         :u0_model_states_init,
         :model_mass_matrix,             

         :cb_states,
         :plants_cb_paras_switches,

         :nodes_names,
         :gens_nodes_names,
         :non_gens_nodes_names,
         :SM_gens_nodes_names,
         :SC_gens_nodes_names,

         :ωref0_vref0_porder0_id_iq_vh_Idx,
         :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

         :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

         :Ybr_cal_and_edges_orientation,
         :Ynet_wt_nodes_idx_wt_adjacent_nodes))


#---------------------------------------------------

# po := post_outage

(po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,) =
    NamedTupleTools.select(
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :post_fault_state),
            :dynamic_status_paras),
        (:ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,))

#---------------------------------------------------

(Ynet, ) =
     NamedTupleTools.select(
    Ynet_wt_nodes_idx_wt_adjacent_nodes,
         (:Ynet, ) )

#---------------------------------------------------

(fault_Ynet,
 post_fault_Ynet) =
    NamedTupleTools.select(
    cleared_selected_lines_faults_net_para,
        (:pre_clear_fault_Ynet,
         :post_clear_fault_Ynet, ))

#---------------------------------------------------

(;dyn_ω_ref_Idx,
 dyn_v_ref_Idx,
 dyn_p_order_Idx,
 dyn_Png_Idx,
 dyn_Qng_Idx,
 dyn_Pll_Idx,
 dyn_Qll_Idx ) =
     NamedTupleTools.select(
         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
         (:dyn_ω_ref_Idx,
          :dyn_v_ref_Idx,
          :dyn_p_order_Idx,
          :dyn_Png_Idx,
          :dyn_Qng_Idx,
          :dyn_Pll_Idx,
          :dyn_Qll_Idx))

#---------------------------------------------------

@show system_fault_status

if system_fault_status[1] != 0

    system_fault_status[1] = 0

end

#---------------------------------------
#---------------------------------------
# line loss with only porder_adj
#---------------------------------------
#---------------------------------------

cb_line_outage = DiscreteCallback(
    (u, t, integrator) ->
        on_line_outage_condition(
            u, t, integrator,
            line_outage_time),

   on_line_outage_affect!;
    save_positions=(true, true),
    initializealg =
        ShampineCollocationInit() )


cb_clear_outage = DiscreteCallback(
    (u, t, integrator) ->
        clear_fault_condition(
            u, t, integrator,
            clear_fault_time),

   clear_outage_affect!;
    save_positions=(true, true),
    initializealg =
        ShampineCollocationInit() )

#---------------------------------------------------

gens_porder_adj =
    po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]

cb_gens_porder_adjustment = DiscreteCallback(
    (u, t, integrator) ->
        on_generation_adjustment_condition(
            u, t, integrator,
            generation_adjustment_time),

    (integrator) ->
        on_generation_adjustment_affect!(
            integrator,
            gens_porder_adj,
            dyn_p_order_Idx );
    save_positions=(true, true),
    initializealg =
        ShampineCollocationInit() )

#---------------------------------------------------
#---------------------------------------------------

vref_and_porder_Idx =
    [dyn_v_ref_Idx;
     dyn_p_order_Idx]

vref_and_porder_adj =
    po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        vref_and_porder_Idx]

cb_vref_and_porder_adj = DiscreteCallback(
    (u, t, integrator) ->
        on_generation_adjustment_condition(
            u, t, integrator,
            generation_adjustment_time),

    (integrator) ->
        on_generation_adjustment_affect!(
        integrator,
        vref_and_porder_adj,
        vref_and_porder_Idx);
    save_positions=(true, true),
    initializealg =
        ShampineCollocationInit() )

#---------------------------------------------------
#---------------------------------------------------

cb_outage_set =
    CallbackSet(
        cb_line_outage,
        cb_clear_outage)

cb_outage_set =
    CallbackSet(cb_line_outage,
    cb_gens_porder_adjustment)

cb_outage_set =
    CallbackSet(cb_line_outage,
    cb_vref_and_porder_adj)

#---------------------------------------------------

tstop_outage =
    [line_outage_time,
     clear_fault_time]

tstop_outage =
    [line_outage_time,
     generation_adjustment_time]

#---------------------------------------------------
#---------------------------------------------------

# generic_model_dynamics_para =
#     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

model_dynamics_para =
    (;generic_model_dynamics_para =
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     Ynet,
     post_fault_Ynet,
     system_fault_status,
     plants_cb_paras_switches )

model_dynamics_kwd_para =
    generic_system_dynamics_wt_fault_kwd_para

#----------------------------------------

system_dynamics_fun! =
    mm_line_outage_generic_dynamics_wt_pre_post_fault_by_ode_pf_funcs!

#----------------------------------------

system_sol =
    DifferentialEquations.solve(
        ODEProblem(
    ODEFunction(
    (dx,x,p,t) ->
        system_dynamics_fun!(
            dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                model_dynamics_kwd_para);
        syms =
            model_syms,
        mass_matrix =
            model_mass_matrix ) ,
            u0_model_states_init,
            sim_timespan,
            model_dynamics_para,
            callback =
                cb_states),
            ode_alg,
            callback =
                cb_outage_set,
            tstops =
                tstop_outage,
            abstol = abstol,
            reltol = reltol )


#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
#---------------------------------------------------

pre_fault_paras =
    getproperty(
        getproperty(
            ntuple_status_steady_state_data,
            :pre_fault_state),
        :static_prefault_paras)


ω_ref_v_ref_p_order_Png_Qng_Pll_Qll =
    getproperty(
        pre_fault_paras,
        :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll)

fault_paras =
    getproperty(
        getproperty(
            ntuple_status_steady_state_data,
            :fault_state),
        :dynamic_status_paras)

post_fault_paras =
    getproperty(
        getproperty(
            ntuple_status_steady_state_data,
            :post_fault_state),
        :dynamic_status_paras)

(;system_fault_status,
 generic_system_dynamics_wt_fault_kwd_para,
 Ynet_wt_nodes_idx_wt_adjacent_nodes,
 on_fault_net_para,
 cleared_selected_lines_faults_net_para,

 model_bool_dae_vars_wt_fault,
 model_syms_wt_fault,         
 u0_model_states_init_wt_fault,

 nodes_names,
 gens_nodes_names,
 non_gens_nodes_names,
 SM_gens_nodes_names,
 SC_gens_nodes_names) =
    NamedTupleTools.select(
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :pre_fault_state),
            :static_prefault_paras),        
        (:system_fault_status,
         :generic_system_dynamics_wt_fault_kwd_para,
         :Ynet_wt_nodes_idx_wt_adjacent_nodes,
         :on_fault_net_para,
         :cleared_selected_lines_faults_net_para,

         :model_bool_dae_vars_wt_fault,
         :model_syms_wt_fault,         
         :u0_model_states_init_wt_fault,

         :nodes_names,
         :gens_nodes_names,
         :non_gens_nodes_names,
         :SM_gens_nodes_names,
         :SC_gens_nodes_names))

#----------------------------------------

(Ynet, ) =
     NamedTupleTools.select(
    Ynet_wt_nodes_idx_wt_adjacent_nodes,
         (:Ynet, ) )

#----------------------------------------

(faulty_Ynet,
 faulty_nodes_idx_with_adjacent_nodes_idx) =
    NamedTupleTools.select(
    on_fault_net_para,
        (:faulty_Ynet,
         :faulty_nodes_idx_with_adjacent_nodes_idx))

(fault_Ynet,
 post_fault_Ynet) =
    NamedTupleTools.select(
    cleared_selected_lines_faults_net_para,
        (:pre_clear_fault_Ynet,
         :post_clear_fault_Ynet,))

#----------------------------------------
# DAE system dyamanics simulation
#----------------------------------------

cb_on_fault = DiscreteCallback(
    (u, t, integrator) ->
        on_fault_condition(
            u, t, integrator,
            on_fault_time),

    on_fault_Ynet_affect!; 
    save_positions=(true,true),
    initializealg =
        ShampineCollocationInit() )

cb_clear_fault = DiscreteCallback(
    (u, t, integrator) ->
        clear_fault_condition(
            u, t, integrator,
            clear_fault_time),

   clear_fault_Ynet_affect!;
    save_positions=(true,true),
    initializealg =
        ShampineCollocationInit())

#--------------------------------------

cb_faults = CallbackSet(
    cb_on_fault,
    cb_clear_fault)

#----------------------------------------

model_dynamics_para =
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     Ynet,
     fault_Ynet,
     post_fault_Ynet,
     system_fault_status )          

model_dynamics_kwd_para =
    generic_system_dynamics_wt_fault_kwd_para

system_dynamics_fun! =
    Ynet_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!

#----------------------------------------

model_bool_dae_vars =
    model_bool_dae_vars_wt_fault

model_syms =
    model_syms_wt_fault

u0_model_states_init =
    u0_model_states_init_wt_fault

#----------------------------------------

du0_model_states_init =
    zeros(length(u0_model_states_init))

res = similar(u0_model_states_init)

#----------------------------------------

faults_and_clear_times =
    [on_fault_time,
     clear_fault_time]

#----------------------------------------

system_sol =
    DifferentialEquations.solve(
        DAEProblem(
    DAEFunction(
    (res, dx, x, p, t) ->
        system_dynamics_fun!(
            res, dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                model_dynamics_kwd_para);
    syms =
        model_syms),
    du0_model_states_init,
    u0_model_states_init,
    sim_timespan,
    model_dynamics_para,
    differential_vars =
        model_bool_dae_vars ),
        dae_alg,
        callback = cb_faults,
        tstops = faults_and_clear_times,
        abstol = abstol,
        reltol = reltol )

#---------------------------------------------------

nt_dae_results =
    (;system_sol,
     model_syms,
     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,
     sim_timespan )

#---------------------------------------------------

plot_dae_dynamics =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            nt_dae_results,
            (:system_sol,
             :model_syms,
             :gens_nodes_names,
             :SM_gens_nodes_names,
             :non_gens_nodes_names,
             :sim_timespan ) )... )

#---------------------------------------------------

# plot_dae_dynamics properties
# (:δ_a_plot, :ω_a_plot, :ed_dash_a_plot,
#  :eq_dash_a_plot, :vr1_a_plot, :vr2_a_plot,
#  :vf_tilade_a_plot, :xg1_a_plot, :xg2_a_plot,
#  :plot_gens_vh, :plot_gens_θh, :plot_non_gens_vh,
#  :plot_non_gens_θh)

δ_a_plot =
    getproperty(
        plot_dae_dynamics,
        :δ_a_plot)

ω_a_plot =
    getproperty(
        plot_dae_dynamics,
        :ω_a_plot)

plot_gens_vh =
    getproperty(
        plot_dae_dynamics,
        :plot_gens_vh)

plot_gens_θh =
    getproperty(
        plot_dae_dynamics,
        :plot_gens_θh)


plot_non_gens_vh =
    getproperty(
        plot_dae_dynamics,
        :plot_non_gens_vh)


plot_non_gens_θh =
    getproperty(
        plot_dae_dynamics,
        :plot_non_gens_θh)

#---------------------------------------------------
#---------------------------------------------------

(s_pf_P_gens,
 s_pf_Q_gens,
 s_vh,
 s_θh,
 s_gens_vh,
 s_gens_θh,
 s_gens_id,
 s_gens_iq,
 
 s_gens_mag_E,
 s_gens_ang_E,
 s_post_sta_PQ,
 s_Yred,
 s_Yint,

 s_ω_ref, s_v_ref, s_p_order,
 s_gens_i_d, s_gens_i_q,
 
 s_flat_vh_flat_θh_id_iq_u0,
 s_flat_vh_flat_θh_id_iq_vfh_θfh,
 s_gens_δ,
 s_gens_ed_dash,
 s_gens_eq_dash) =
     NamedTupleTools.select(
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :pre_fault_state),
            :static_prefault_paras),
        (:pf_P_gens,
         :pf_Q_gens,
         :vh,
         :θh,
         :gens_vh,
         :gens_θh,
         :gens_id,
         :gens_iq,
         
         :gens_mag_E,
         :gens_ang_E,
         :post_sta_PQ,
         :Yred,
         :Yint,
         
         :ω_ref, :v_ref, :p_order,
         :gens_i_d, :gens_i_q,
         
         :flat_vh_flat_θh_id_iq_u0,
         :flat_vh_flat_θh_id_iq_vfh_θfh,
         :gens_δ,
         :gens_ed_dash,
         :gens_eq_dash ) )


(f_dyn_pf_P_gens,
 f_dyn_pf_Q_gens,
 f_dyn_vh,
 f_dyn_θh,
 f_dyn_gens_vh,
 f_dyn_gens_θh,
 f_dyn_gens_id,
 f_dyn_gens_iq,
 f_dyn_gens_mag_E,
 f_dyn_gens_ang_E,
 f_post_dyn_PQ,
 f_dyn_Yred,
 f_dyn_Yint,
 f_flat_vh_flat_θh_id_iq_vfh_θfh,
 f_dyn_gens_δ,
 f_dyn_gens_ed_dash,
 f_dyn_gens_eq_dash) =
     NamedTupleTools.select(
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :fault_state),
            :dynamic_status_paras),
        (:dyn_pf_P_gens,
         :dyn_pf_Q_gens,
         :dyn_vh,
         :dyn_θh,
         :dyn_gens_vh,
         :dyn_gens_θh,
         :dyn_gens_id,
         :dyn_gens_iq,
         :dyn_gens_mag_E,
         :dyn_gens_ang_E,
         :post_dyn_PQ,
         :dyn_Yred, :dyn_Yint,
         :flat_vh_flat_θh_id_iq_vfh_θfh,
         :dyn_gens_δ,
         :dyn_gens_ed_dash,
         :dyn_gens_eq_dash ) )

#----------------------------------------

(p_dyn_pf_P_gens,
 p_dyn_pf_Q_gens,
 p_dyn_vh,
 p_dyn_θh,
 p_dyn_gens_vh,
 p_dyn_gens_θh,
 p_dyn_gens_id,
 p_dyn_gens_iq,
 p_dyn_gens_mag_E,
 p_dyn_gens_ang_E,
 p_post_dyn_PQ,
 p_dyn_Yred,
 p_dyn_Yint,
 p_flat_vh_flat_θh_id_iq_vfh_θfh,
 p_dyn_gens_δ,
 p_dyn_gens_ed_dash,
 p_dyn_gens_eq_dash) =
     NamedTupleTools.select(
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :post_fault_state),
            :dynamic_status_paras),
        (:dyn_pf_P_gens,
         :dyn_pf_Q_gens,
         :dyn_vh,
         :dyn_θh,
         :dyn_gens_vh,
         :dyn_gens_θh,
         :dyn_gens_id,
         :dyn_gens_iq,
         :dyn_gens_mag_E,
         :dyn_gens_ang_E,
         :post_dyn_PQ,
         :dyn_Yred, :dyn_Yint,
         :flat_vh_flat_θh_id_iq_vfh_θfh,
         :dyn_gens_δ,
         :dyn_gens_ed_dash,
         :dyn_gens_eq_dash) )

#---------------------------------------------------
#---------------------------------------------------


"""
# Possible system_status

system_status = :pre_fault_state,
system_status = :fault_state
system_status = :post_fault_state

"""


"""
# outage_type = :line_outage    
# outage_type = :line_outage_wt_pref_adjs
# outage_type = :line_outage_wt_vpref_adjs
 
"""

line_outage_pertubation_by_mm_ode = 
sim_line_outage_pertubation_by_mm_ode(
    :line_outage_wt_pref_adjs, # outage_type

    on_fault_time,
    clear_fault_time,

    line_outage_time,
    generation_adjustment_time;

    net_data_by_components_file,

    timespan,

    components_libs_dir,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit  = [1],

    list_edges_to_have_fault  = [ 8 ],
    clear_fault_selection_list = [1] )

#---------------------------------------------------
#---------------------------------------------------
