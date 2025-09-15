#=
# [Example of Ordinary Power Flow with IEEE 5 Bus Test System](@id ieee-5-bus-ordinary-power-flow)
=#

#---------------------------------------------------
#---------------------------------------------------

using Revise

# using Pkg

using ePowerSim

#---------------------------------------------------
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



case_name = "case5"


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


sim_type  = "sim-dynmaic-simulation"

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


abstol      = 1e-12

reltol      = 1e-12


#---------------------------------------------------
#---------------------------------------------------

system_simulation_parameters_nothing =    
    get_system_simulation_parameters(
        net_data_by_components_file,
        nothing;
        components_libs_dir =
            components_libs_dir,
        basekV = 1.0,    
        use_pu_in_PQ = true,
        line_data_in_pu = true,

        use_init_u0 = false,
        use_nlsolve = false,
        pf_alg = NewtonRaphson(),
        abstol = abstol,
        reltol = reltol ) 

#---------

system_simulation_parameters =
    get_system_simulation_parameters(
        net_data_by_components_file;
        components_libs_dir =
            components_libs_dir,
        basekV          = 1.0,    
        use_pu_in_PQ    = true,
        line_data_in_pu = true,

        use_init_u0 = false,    
        use_nlsolve = false,    
        pf_alg      = NewtonRaphson(),    
        abstol      = abstol,
         reltol     = reltol  )

#---------------------------------------------------

(;u0_model_states_init,
 model_syms,
 model_mass_matrix,
 model_bool_dae_vars,

 ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
 plants_cb_paras_switches,

 generic_system_dynamics_kwd_para,

 gens_nodes_names,
 SM_gens_nodes_names,
 non_gens_nodes_names,

 cb_states,

 plants_states_syms,

 gens_nodes_idx,

 state_labels,
 algebraic_vars_labels,
 network_vars_labels,

 dyn_pf_fun_kwd_n2s_idxs,
 dyn_pf_fun_kwd_net_idxs,

 # Ybr_cal_and_edges_orientation,
 Ynet_wt_nodes_idx_wt_adjacent_nodes) =
     NamedTupleTools.select(
         get_system_simulation_parameters(
             net_data_by_components_file;
             components_libs_dir =
                 components_libs_dir,
             basekV = basekV,    
             use_pu_in_PQ = use_pu_in_PQ,
             line_data_in_pu = line_data_in_pu,

             use_init_u0 = use_init_u0,
             use_nlsolve = use_nlsolve,

             pf_alg = pf_alg,

             abstol = abstol,
             reltol = reltol),

         (:u0_model_states_init,
          :model_syms,
          :model_mass_matrix,
          :model_bool_dae_vars,

          :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
          :plants_cb_paras_switches,

          :generic_system_dynamics_kwd_para,

          :gens_nodes_names,
          :SM_gens_nodes_names,
          :non_gens_nodes_names,

          :cb_states,

          :plants_states_syms,
          :gens_nodes_idx,

          :state_labels,
          :algebraic_vars_labels,
          :network_vars_labels,

          :dyn_pf_fun_kwd_n2s_idxs,
          :dyn_pf_fun_kwd_net_idxs,

          # :Ybr_cal_and_edges_orientation,
          :Ynet_wt_nodes_idx_wt_adjacent_nodes))

#---------------------------------------------------

(;loc_load_exist,
 dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
 dyn_pf_fun_kwd_n2s_idxs,
 dyn_pf_fun_kwd_net_idxs,

 gens_state_vars_idx_in_state,
 state_vars_and_i_dq_Idx_in_state,

 state_labels,
 algebraic_vars_labels) =
     NamedTupleTools.select(
         generic_system_dynamics_kwd_para ,
         (:loc_load_exist,
          :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
          :dyn_pf_fun_kwd_n2s_idxs,
          :dyn_pf_fun_kwd_net_idxs,

          :gens_state_vars_idx_in_state,
          :state_vars_and_i_dq_Idx_in_state,

          :state_labels,
          :algebraic_vars_labels))

#---------------------------------------------------
# ode_generic_system_dynamics_by_ode_pf_funcs!
#---------------------------------------------------

model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

model_dynamics_kwd_para =
    ode_generic_model_dynamics_kwd_para

ODE_system_dynamics_fun! =
    ode_generic_system_dynamics_by_ode_pf_funcs!



system_sol =
    DifferentialEquations.solve(
        ODEProblem(
    ODEFunction(
    (dx,x,p,t) ->
        ODE_system_dynamics_fun!(
            dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                model_dynamics_kwd_para);
    syms =
        model_syms,
mass_matrix = model_mass_matrix ),
    u0_model_states_init,
    sim_timespan,
    model_dynamics_para ),
        ode_alg,
        abstol = abstol,
        reltol = reltol )

nt_ode_results =
    (;system_sol,
     model_syms,
     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,
     sim_timespan)

#---------------------------------------------------
# ode_generic_system_model_by_funcs_dynamics!
#---------------------------------------------------

model_dynamics_para =
    (;generic_model_dynamics_para =
    ωref0_vref0_porder0_id_iq_vh,
      plants_cb_paras_switches)

model_dynamics_kwd_para =
    ode_generic_model_dynamics_kwd_para

ODE_system_dynamics_fun! =
    ode_generic_system_model_by_funcs_dynamics!

system_sol =
    DifferentialEquations.solve(
        ODEProblem(
    ODEFunction(
    (dx,x,p,t) ->
        ODE_system_dynamics_fun!(
            dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                model_dynamics_kwd_para);
    syms =
        model_syms,
mass_matrix = model_mass_matrix ),
    u0_model_states_init,
    sim_timespan,
    model_dynamics_para ),
        ode_alg,
        abstol = abstol,
        reltol = reltol )

nt_ode_results =
    (;system_sol,
     model_syms,
     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,
     sim_timespan)

#---------------------------------------------------
# ode_generic_system_model_by_funcs_dynamics!
#---------------------------------------------------

model_dynamics_para =
    (;generic_model_dynamics_para =
    ωref0_vref0_porder0_id_iq_vh,
      plants_cb_paras_switches)

model_dynamics_kwd_para =
    ode_generic_model_dynamics_kwd_para

ODE_system_dynamics_fun! =
    ode_generic_system_model_by_funcs_dynamics!



system_sol =
    DifferentialEquations.solve(
        ODEProblem(
    ODEFunction(
    (dx,x,p,t) ->
        ODE_system_dynamics_fun!(
            dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                model_dynamics_kwd_para);
    syms =
        model_syms,
mass_matrix = model_mass_matrix ),
    u0_model_states_init,
    sim_timespan,
    model_dynamics_para ),
        ode_alg,
        abstol = abstol,
        reltol = reltol )

nt_ode_results =
    (;system_sol,
     model_syms,
     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,
     sim_timespan)

#---------------------------------------------------
#---------------------------------------------------
# DAE system dyamanics simulation
#---------------------------------------------------
#---------------------------------------------------

#---------------------------------------------------
# dae_generic_system_dynamics_by_dae_pf_funcs!
#---------------------------------------------------

model_dynamics_para =
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll


dae_model_dynamics_kwd_para =
    dae_generic_model_kwd_para


DAE_system_dynamics_fun! =
    dae_generic_system_dynamics_by_dae_pf_funcs!


system_sol =
    DifferentialEquations.solve(
        DAEProblem(
    DAEFunction(
    (res, dx,x,p,t) ->
        DAE_system_dynamics_fun!(
            res, dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                dae_model_dynamics_kwd_para);
    syms =
        model_syms),
    du0_model_states_init,
    u0_model_states_init,
    sim_timespan,
    model_dynamics_para,
    differential_vars =
        model_bool_dae_vars ),
        dae_alg,
        callback = cb,
        abstol   = abstol,
        reltol   = reltol )


nt_dae_results =  (;system_sol,
        model_syms,
        gens_nodes_names,
        SM_gens_nodes_names,
        non_gens_nodes_names,
        sim_timespan)


#---------------------------------------------------
# dae_generic_system_model_by_funcs_dynamics!
#---------------------------------------------------

model_dynamics_para =
    (;generic_model_dynamics_para =
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
      plants_cb_paras_switches)

dae_model_dynamics_kwd_para =
    dae_generic_model_kwd_para


DAE_system_dynamics_fun! =
    dae_generic_system_model_by_funcs_dynamics!


system_sol =
    DifferentialEquations.solve(
        DAEProblem(
    DAEFunction(
    (res, dx,x,p,t) ->
        DAE_system_dynamics_fun!(
            res, dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                dae_model_dynamics_kwd_para);
    syms =
        model_syms),
    du0_model_states_init,
    u0_model_states_init,
    sim_timespan,
    model_dynamics_para,
    differential_vars =
        model_bool_dae_vars ),
        dae_alg,
        callback = cb,
        abstol   = abstol,
        reltol   = reltol )


nt_dae_results =  (;system_sol,
        model_syms,
        gens_nodes_names,
        SM_gens_nodes_names,
        non_gens_nodes_names,
        sim_timespan)

#---------------------------------------------------
# dae_generic_system_model_dynamics!
#---------------------------------------------------

model_dynamics_para =
    (;generic_model_dynamics_para =
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
      plants_cb_paras_switches)

dae_model_dynamics_kwd_para =
    dae_generic_model_dynamics_kwd_para


DAE_system_dynamics_fun! =
    dae_generic_system_model_dynamics!


system_sol =
    DifferentialEquations.solve(
        DAEProblem(
    DAEFunction(
    (res, dx,x,p,t) ->
        DAE_system_dynamics_fun!(
            res, dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                dae_model_dynamics_kwd_para);
    syms =
        model_syms),
    du0_model_states_init,
    u0_model_states_init,
    sim_timespan,
    model_dynamics_para,
    differential_vars =
        model_bool_dae_vars ),
        dae_alg,
        callback = cb,
        abstol   = abstol,
        reltol   = reltol )


nt_dae_results =  (;system_sol,
        model_syms,
        gens_nodes_names,
        SM_gens_nodes_names,
        non_gens_nodes_names,
        sim_timespan)

#---------------------------------------------------
#---------------------------------------------------
# Mass matrix ODE system dyamanics simulation
#---------------------------------------------------
#---------------------------------------------------

#---------------------------------------------------
# mm_ode_generic_system_dynamics_by_ode_pf_funcs!
#---------------------------------------------------

model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

model_dynamics_kwd_para =
    mm_generic_model_kwd_para

mm_ODE_system_dynamics_fun! =
    mm_ode_generic_system_dynamics_by_ode_pf_funcs!


system_sol =
    DifferentialEquations.solve(
        ODEProblem(
    ODEFunction(
    (dx,x,p,t) ->
        mm_ODE_system_dynamics_fun!(
            dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                model_dynamics_kwd_para);
        syms =
            model_syms,
        mass_matrix =
            model_mass_matrix ),
            u0_model_states_init,
            sim_timespan,
            model_dynamics_para ),
        ode_alg,
        abstol = abstol,
        reltol = reltol )

nt_mm_ode_results =
    (;system_sol,
     model_syms,
     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,
     sim_timespan)

#---------------------------------------------------
# ode_generic_system_model_by_funcs_dynamics!
#---------------------------------------------------

model_dynamics_para =
    (;generic_model_dynamics_para =
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
      plants_cb_paras_switches)

model_dynamics_kwd_para =
    mm_generic_model_kwd_para

##
mm_ODE_system_dynamics_fun! =
    mm_ode_generic_system_model_by_funcs_dynamics!


system_sol =
    DifferentialEquations.solve(
        ODEProblem(
    ODEFunction(
    (dx,x,p,t) ->
        mm_ODE_system_dynamics_fun!(
            dx, x,
            model_dynamics_para,
            t;
            kwd_para =
                model_dynamics_kwd_para);
        syms =
            model_syms,
        mass_matrix =
            model_mass_matrix ),
            u0_model_states_init,
            sim_timespan,
            model_dynamics_para ),
        ode_alg,
        abstol = abstol,
        reltol = reltol )

nt_mm_ode_results =
    (;system_sol,
     model_syms,
     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,
     sim_timespan)

#---------------------------------------------------
#---------------------------------------------------
