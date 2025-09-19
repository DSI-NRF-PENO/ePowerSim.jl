#=
# [Example of Network Perturbation with IEEE 5 Bus Test System](@id ieee-5-bus-network-perturbation)
=#

#=
# [Example of Network Perturbation with IEEE 5 Bus Test System](@id ieee-5-bus-network-perturbation)
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


sauer_net_wt_avr_string    = "net-static-data-avr-sauer"

rtds_net_wt_avr_string     = "net-static-data-avr-rtds"

sauer_gov_string           = "gov-sauer"

ieee_tgov2sauer_gov_string = "gov-ieee-tgov2sauer"

ieee_tgov1_gov_string      = "gov-ieee-tgov1"

#---------------------------------------------------

# # rtds

# net_wt_avr_string = rtds_net_wt_avr_string
# gov_string        = ieee_tgov1_gov_string


# sauer

net_wt_avr_string = sauer_net_wt_avr_string
gov_string        = sauer_gov_string


dynamic_net_data_by_components_file =
    "$(net_wt_avr_string)-$(gov_string).json"

json_net_data_by_components_file =
    dynamic_net_data_by_components_file


# net_wt_avr_string   = "net-static-data-avr-sauer"

# gov_string                 = "gov-sauer"

# dynamic_net_data_by_components_file =
#     "$(net_wt_avr_string)-$(gov_string).json"

# json_net_data_by_components_file =
#     dynamic_net_data_by_components_file

#---------------------------------------------------
#---------------------------------------------------

case_name = "case9"

# case_name = "case14"

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

sim_type  = "$(case_name)-"*"test-"*"network-pertubation"*
    "-$(net_wt_avr_string)-$(gov_string)"


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

sim_duration  = 30.0

timespan_min  = restoration_time + Δt1 + Δt2

timespan      =
    timespan_min < sim_duration ? sim_duration : timespan_min

time_start    = 0.0

time_final    = timespan

tspan         = (0.0, timespan)

sim_timespan  = (0.0, timespan)

plot_timespan = (0.0, timespan)


#---------------------------------------------------
# A generalised function that takes functions
# as a namedtuple
# 
# get_generic_line_loss_outage_wt_or_no_ref_adjs
# 
#---------------------------------------------------

list_network_status = 
    [:pre_fault_state,
     :fault_state,
     :post_fault_state ]

nt_system_dynamics_fun_type =
        (;system_dynamics_fun_type =
        :dae_line_loss_pre_fault_post_pf_funcs,
        system_dynamics_fun =
        line_loss_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs! )

 
generic_plot_line_outage =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            get_generic_line_loss_outage_wt_or_no_ref_adjs(
                :line_outage, # outage_type,    

                   on_fault_time,
                   clear_fault_time,

                   line_outage_time,
                   generation_adjustment_time,

                   net_data_by_components_file;

                   timespan,    

                   with_faults,   

                   components_libs_dir,

                   data_dir, 

                   list_fault_point_from_node_a = [0.3],
                   list_fault_resistance = [0.001],
                   list_no_line_circuit  = [1],

                   list_edges_to_have_fault   = [ 8 ],
                   clear_fault_selection_list = [1],

                   basekV          = 1.0,    
                   use_pu_in_PQ    = true,
                   line_data_in_pu = true,

                   use_init_u0 = false,
                   use_nlsolve = false,

                   pf_alg  = NewtonRaphson(),
                   ode_alg = Rodas4(),
                   dae_alg = IDA(),

                   abstol = 1e-12,
                   reltol = 1e-12,

                   dt = 0.0001,
                   Δt = 1.0 / 2^(4),
                list_network_status,
                nt_system_dynamics_fun_type) ,
                (:system_sol,
                 :model_syms,
                 :gens_nodes_names,
                 :SM_gens_nodes_names,
                 :non_gens_nodes_names,
                 :sim_timespan) )... )


generic_δ_a_plot = getproperty(
    generic_plot_line_outage, :δ_a_plot)


generic_ω_a_plot  = getproperty(
    generic_plot_line_outage, :ω_a_plot)


# propertynames( plot_line_outage )

# (:δ_a_plot, :ω_a_plot, :ed_dash_a_plot,
#  :eq_dash_a_plot, :vr1_a_plot, :vr2_a_plot,
#  :vf_tilade_a_plot, :xg1_a_plot, :xg2_a_plot,
#  :plot_gens_vh, :plot_gens_θh,
#  :plot_non_gens_vh, :plot_non_gens_θh)


#---------------------------------------------------
#---------------------------------------------------

list_network_status = 
    [:pre_fault_state,
     :fault_state,
     :post_fault_state ]

 
plot_line_outage =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            get_line_loss_outage_wt_or_no_ref_adjs(
                :line_outage, # outage_type,    

                   on_fault_time,
                   clear_fault_time,

                   line_outage_time,
                   generation_adjustment_time,

                   net_data_by_components_file;

                   timespan,    

                   with_faults,   

                   components_libs_dir,

                   data_dir, 

                   list_fault_point_from_node_a = [0.3],
                   list_fault_resistance = [0.001],
                   list_no_line_circuit  = [1],

                   list_edges_to_have_fault   = [ 8 ],
                   clear_fault_selection_list = [1],

                   basekV          = 1.0,    
                   use_pu_in_PQ    = true,
                   line_data_in_pu = true,

                   use_init_u0 = false,
                   use_nlsolve = false,

                   pf_alg  = NewtonRaphson(),
                   ode_alg = Rodas4(),
                   dae_alg = IDA(),

                   abstol = 1e-12,
                   reltol = 1e-12,

                   dt = 0.0001,
                   Δt = 1.0 / 2^(4),
                   list_network_status ) ,
                (:system_sol,
                 :model_syms,
                 :gens_nodes_names,
                 :SM_gens_nodes_names,
                 :non_gens_nodes_names,
                 :sim_timespan) )... )


δ_a_plot = getproperty(
    plot_line_outage, :δ_a_plot)


ω_a_plot  = getproperty(
    plot_line_outage, :ω_a_plot)


# propertynames( plot_line_outage )

# (:δ_a_plot, :ω_a_plot, :ed_dash_a_plot,
#  :eq_dash_a_plot, :vr1_a_plot, :vr2_a_plot,
#  :vf_tilade_a_plot, :xg1_a_plot, :xg2_a_plot,
#  :plot_gens_vh, :plot_gens_θh,
#  :plot_non_gens_vh, :plot_non_gens_θh)

#---------------------------------------------------


plot_line_outage_wt_pref_adjs =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            get_line_loss_outage_wt_or_no_ref_adjs(
                :line_outage_wt_pref_adjs, # outage_type,    

                   on_fault_time,
                   clear_fault_time,

                   line_outage_time,
                   generation_adjustment_time,

                   net_data_by_components_file;

                   timespan,    

                   with_faults,   

                   components_libs_dir,

                   data_dir, 

                   list_fault_point_from_node_a = [0.3],
                   list_fault_resistance = [0.001],
                   list_no_line_circuit  = [1],

                   list_edges_to_have_fault   = [ 8 ],
                   clear_fault_selection_list = [1],

                   basekV          = 1.0,    
                   use_pu_in_PQ    = true,
                   line_data_in_pu = true,

                   use_init_u0 = false,
                   use_nlsolve = false,

                   pf_alg  = NewtonRaphson(),
                   ode_alg = Rodas4(),
                   dae_alg = IDA(),

                   abstol = 1e-12,
                   reltol = 1e-12,

                   dt = 0.0001,
                   Δt = 1.0 / 2^(4),
                   list_network_status ) ,
                (:system_sol,
                 :model_syms,
                 :gens_nodes_names,
                 :SM_gens_nodes_names,
                 :non_gens_nodes_names,
                 :sim_timespan) )... )


wt_pref_adjs_δ_a_plot =
    getproperty(
        plot_line_outage_wt_pref_adjs,
        :δ_a_plot)


wt_pref_adjs_ω_a_plot  =
    getproperty(
        plot_line_outage_wt_pref_adjs,
        :ω_a_plot)


#---------------------------------------------------


plot_line_outage_wt_vpref_adjs =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            get_line_loss_outage_wt_or_no_ref_adjs(
                :line_outage_wt_vpref_adjs, # outage_type,    

                   on_fault_time,
                   clear_fault_time,

                   line_outage_time,
                   generation_adjustment_time,

                   net_data_by_components_file;

                   timespan,    

                   with_faults,   

                   components_libs_dir,

                   data_dir, 

                   list_fault_point_from_node_a = [0.3],
                   list_fault_resistance = [0.001],
                   list_no_line_circuit  = [1],

                   list_edges_to_have_fault   = [ 8 ],
                   clear_fault_selection_list = [1],

                   basekV          = 1.0,    
                   use_pu_in_PQ    = true,
                   line_data_in_pu = true,

                   use_init_u0 = false,
                   use_nlsolve = false,

                   pf_alg  = NewtonRaphson(),
                   ode_alg = Rodas4(),
                   dae_alg = IDA(),

                   abstol = 1e-12,
                   reltol = 1e-12,

                   dt = 0.0001,
                   Δt = 1.0 / 2^(4),
                   list_network_status ) ,
                (:system_sol,
                 :model_syms,
                 :gens_nodes_names,
                 :SM_gens_nodes_names,
                 :non_gens_nodes_names,
                 :sim_timespan) )... )



wt_vpref_adjs_δ_a_plot =
    getproperty(
        plot_line_outage_wt_vpref_adjs,
        :δ_a_plot)


wt_vpref_adjs_ω_a_plot  =
    getproperty(
        plot_line_outage_wt_vpref_adjs,
        :ω_a_plot)


#---------------------------------------------------
#---------------------------------------------------


list_network_status = 
    [:pre_fault_state,
     :fault_state,
     :post_fault_state ]


ntuple_status_steady_state_data =
    get_ntuple_status_steady_state_data(
        ;with_faults =
            with_faults,
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
            line_data_in_pu,
        list_network_status =
            list_network_status)

#---------------------------------------------------

line_loss_outage = 
    get_line_loss_outage_wt_or_no_ref_adjs(
        :line_outage, # outage_type,       
        ntuple_status_steady_state_data,
        line_outage_time,
        generation_adjustment_time;
        sim_timespan =
            sim_timespan,    
        dae_alg =
            dae_alg,
        abstol  =
            abstol,
        reltol  =
            reltol )

plot_line_outage =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            line_loss_outage,
            (:system_sol,
             :model_syms,
             :gens_nodes_names,
             :SM_gens_nodes_names,
             :non_gens_nodes_names,
             :sim_timespan) )... )

δ_a_plot = getproperty(
    plot_line_outage,
    :δ_a_plot)


ω_a_plot = getproperty(
    plot_line_outage,
    :ω_a_plot)


# (:δ_a_plot, :ω_a_plot, :ed_dash_a_plot,
#  :eq_dash_a_plot, :vr1_a_plot, :vr2_a_plot,
#  :vf_tilade_a_plot, :xg1_a_plot, :xg2_a_plot,
#  :plot_gens_vh, :plot_gens_θh,
#  :plot_non_gens_vh, :plot_non_gens_θh)

# list_outage_type =
#     [:line_outage,
#      :line_outage_wt_pref_adjs,
#      :line_outage_wt_vpref_adjs]



#---------------------------------------------------


plot_line_outage_wt_pref_adjs =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            get_line_loss_outage_wt_or_no_ref_adjs(
                :line_outage_wt_pref_adjs, # outage_type,       
                ntuple_status_steady_state_data,
                line_outage_time,
                generation_adjustment_time;
                sim_timespan =
                    sim_timespan,    
                dae_alg =
                    dae_alg,
                abstol  =
                    abstol,
                reltol  =
                    reltol ) ,
                (:system_sol,
                 :model_syms,
                 :gens_nodes_names,
                 :SM_gens_nodes_names,
                 :non_gens_nodes_names,
                 :sim_timespan) )... )


plot_line_outage_wt_pref_adjs =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            plot_line_outage_wt_pref_adjs,
            (:system_sol,
             :model_syms,
             :gens_nodes_names,
             :SM_gens_nodes_names,
             :non_gens_nodes_names,
             :sim_timespan) )... )

δ_a_plot = getproperty(
    plot_line_outage_wt_pref_adjs,
    :δ_a_plot)


ω_a_plot = getproperty(
    plot_line_outage_wt_pref_adjs,
    :ω_a_plot)


# (:δ_a_plot, :ω_a_plot, :ed_dash_a_plot,
#  :eq_dash_a_plot, :vr1_a_plot, :vr2_a_plot,
#  :vf_tilade_a_plot, :xg1_a_plot, :xg2_a_plot,
#  :plot_gens_vh, :plot_gens_θh,
#  :plot_non_gens_vh, :plot_non_gens_θh)

#---------------------------------------------------


plot_line_outage_wt_vref_adjs =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            get_line_loss_outage_wt_or_no_ref_adjs(
        :line_outage_wt_vpref_adjs, # outage_type,       
            ntuple_status_steady_state_data,
            line_outage_time,
            generation_adjustment_time;
            sim_timespan =
                sim_timespan,    
            dae_alg =
                dae_alg,
            abstol  =
                abstol,
            reltol  =
                reltol ) ,
            (:system_sol,
             :model_syms,
             :gens_nodes_names,
             :SM_gens_nodes_names,
             :non_gens_nodes_names,
             :sim_timespan) )... )


plot_line_outage_wt_vpref_adjs =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            plot_line_outage_wt_vref_adjs ,
            (:system_sol,
             :model_syms,
             :gens_nodes_names,
             :SM_gens_nodes_names,
             :non_gens_nodes_names,
             :sim_timespan) )... )

δ_a_plot = getproperty(
    plot_line_outage_wt_vpref_adjs,
    :δ_a_plot)


ω_a_plot = getproperty(
    plot_line_outage_wt_vpref_adjs,
    :ω_a_plot)


# (:δ_a_plot, :ω_a_plot, :ed_dash_a_plot,
#  :eq_dash_a_plot, :vr1_a_plot, :vr2_a_plot,
#  :vf_tilade_a_plot, :xg1_a_plot, :xg2_a_plot,
#  :plot_gens_vh, :plot_gens_θh,
#  :plot_non_gens_vh, :plot_non_gens_θh)

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

#---------------------------------------------------
#---------------------------------------------------

(;sta_pf_red_sol,
 dyn_pf_fun_kwd_net_idxs) =
    NamedTupleTools.select(
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :pre_fault_state),
            :static_prefault_paras),
        (:sta_pf_red_sol,
         :dyn_pf_fun_kwd_net_idxs))

#---------------------------------------------------

(;pf_P_gens,
  pf_Q_gens,
  vh,
  θh,
  θh_deg) =
    NamedTupleTools.select(
        sta_pf_red_sol,
        (:pf_P_gens,
         :pf_Q_gens,
         :vh,
         :θh,
         :θh_deg))

(gens_nodes_idx,
 non_gens_nodes_idx,
 gens_nodes_with_loc_loads_idx,) =
    NamedTupleTools.select(
        dyn_pf_fun_kwd_net_idxs,
        (:gens_nodes_idx,
         :non_gens_nodes_idx,
         :gens_nodes_with_loc_loads_idx))

#---------------------------------------------------

gens_vh =
    round.( vh[gens_nodes_idx]; digits=4)


gens_θh =
    round.( θh[gens_nodes_idx]; digits=4)


gens_θh_deg =
    round.( θh_deg[gens_nodes_idx]; digits=4)

#---------------------------------------------------

non_gens_vh =
    round.( vh[non_gens_nodes_idx]; digits=4)

non_gens_θh =
    round.( θh[non_gens_nodes_idx]; digits=4)

non_gens_θh_deg =
    round.( θh_deg[non_gens_nodes_idx]; digits=4)

#---------------------------------------------------

t_pf_P_gens = round.( pf_P_gens; digits=4)

t_pf_Q_gens = round.( pf_Q_gens; digits=4)

#---------------------------------------------------

# Gens

pf_gens_results =
    hcat(gens_vh,
         gens_θh_deg,
         t_pf_P_gens,
         t_pf_Q_gens)

#---------------------------------------------------

# Non gens

pf_non_gens_results =
    hcat(non_gens_vh,
         non_gens_θh_deg )

#---------------------------------------------------


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

#---------------------------------------------------

# generic_network_fault_pertubation_plot =
#     get_guick_group_vars_plots_dae_or_ode_sol(
#         ;include_v_θ_plot =
#             false,
#          sim_timespan =
#              (0, timespan),
#          get_generic_network_fault_pertubation(
#              ;case_name = "case9",
#              timespan   = timespan,

#              on_fault_time = 5.0,
#              clear_fault_time = 7.0,

#              list_fault_point_from_node_a = [0.3],
#              list_fault_resistance = [0.001],
#              list_no_line_circuit =  [1],

#              list_edges_to_have_fault = [ 8 ],
#              clear_fault_selection_list = [1],

#              basekV = 1.0,    
#              use_pu_in_PQ = true,
#              line_data_in_pu = true,

#              with_faults =
#                  false,
#              use_state_in_on_clear_fault =
#                  false,
#              return_extended_results =
#                  false,

#              json_net_data_by_components_file =
#                  json_net_data_by_components_file,
#              components_libs_dir =
#                  components_libs_dir,
#              data_dir =
#                  data_dir )...)


# save_network_pertubation_sim_plot(
#     "case9";
#     sim_timespan = (0, timespan),
#     figure_dir,
#     sim_type = "network-pertubation",
#     line_in_fault_name = "line-8",
#     get_generic_network_fault_pertubation_by_cb(
#              ;case_name = "case9",
#              timespan   = timespan,

#              on_fault_time = 5.0,
#              clear_fault_time = 7.0,

#              list_fault_point_from_node_a = [0.3],
#              list_fault_resistance = [0.001],
#              list_no_line_circuit =  [1],

#              list_edges_to_have_fault = [ 8 ],
#              clear_fault_selection_list = [1],

#              basekV = 1.0,    
#              use_pu_in_PQ = true,
#              line_data_in_pu = true,

#              with_faults =
#                  false,
#              use_state_in_on_clear_fault =
#                  false,
#              return_extended_results =
#                  false,

#              json_net_data_by_components_file =
#                  json_net_data_by_components_file,
#              components_libs_dir =
#                  components_libs_dir,
#              data_dir =
#                  data_dir )...)




#---------------------------------------------------
#---------------------------------------------------
