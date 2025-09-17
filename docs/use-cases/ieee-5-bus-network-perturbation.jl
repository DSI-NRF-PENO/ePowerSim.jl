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
# Dynamic modeling and simulations
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


# json_net_data_by_components_file =
#     "net-static-data-avr-sauer-gov-sauer.json"


json_net_data_by_components_file =
    dynamic_net_data_by_components_file

#---------------------------------------------------
#---------------------------------------------------

case_name = "case9"


# case_name = "case14"


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

sim_type  = "sim-network-perturbation"

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

"""

# Possible system_status

system_status = :pre_fault_state,
system_status = :fault_state
system_status = :post_fault_state

"""


line_outage_pertubation_by_mm_ode = 
sim_line_outage_pertubation_by_mm_ode(
    :line_outage_wt_pref_adjs, # outage_type

    on_fault_time,
    clear_fault_time,

    line_outage_time,
    generation_adjustment_time;

    net_data_by_components_file,

    timespan = 50,

    components_libs_dir,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit  = [1],

    list_edges_to_have_fault  = [ 8 ],
    clear_fault_selection_list = [1] )


#----------------------------------------
#----------------------------------------


list_network_status = 
            [:pre_fault_state,
             :post_fault_state]

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
            list_network_status )


pre_fault_steady_state_data =
    get_a_status_steady_state_data(
        :pre_fault_state;
        with_faults =
            false,

        # :pre_fault_state
        # :fault_state 
        # :post_fault_state ,
        # system_status = :pre_fault_state,

        net_data_by_components_file =
            net_data_by_components_file,
        components_libs_dir =
            components_libs_dir,

        timespan         = 10.0,
        on_fault_time    = 9.0,
        clear_fault_time = 9.001,

        list_fault_point_from_node_a = [0.3],
        list_fault_resistance        = [0.001],
        list_no_line_circuit         =  [4],

        list_edges_to_have_fault   = [ 2 ],
        clear_fault_selection_list = [1],

        basekV          = 1.0,    
        use_pu_in_PQ    = true,
        line_data_in_pu = true,

        #--------------------------------------    

        use_init_u0 = false,

        use_nlsolve = false,

        pf_alg        = NewtonRaphson(),

        #--------------------------------------    

        ode_alg       = Rodas4(),
        # ode_alg     = FBDF()
        # ode_alg     = QNDF()
        # ode_alg     = radau()
        # ode_alg     = RadauIIA5()
        # ode_alg     = DFBDF()

        dae_alg       = IDA(),

        abstol        = 1e-12,
        reltol        = 1e-12,

        dt            = 0.01)

#--------------

fault_steady_state_data =
    get_a_status_steady_state_data(
        :fault_state;
        with_faults =
            false,
        
        # :pre_fault_state
        # :fault_state 
        # :post_fault_state ,
        # system_status = :pre_fault_state,

        net_data_by_components_file =
            net_data_by_components_file,
        components_libs_dir =
            components_libs_dir,

        timespan         = 10.0,
        on_fault_time    = 9.0,
        clear_fault_time = 9.001,

        list_fault_point_from_node_a = [0.3],
        list_fault_resistance        = [0.001],
        list_no_line_circuit         =  [4],

        list_edges_to_have_fault   = [ 2 ],
        clear_fault_selection_list = [1],

        basekV          = 1.0,    
        use_pu_in_PQ    = true,
        line_data_in_pu = true,

        #--------------------------------------    

        use_init_u0 = false,

        use_nlsolve = false,

        pf_alg        = NewtonRaphson(),

        #--------------------------------------    

        ode_alg       = Rodas4(),
        # ode_alg     = FBDF()
        # ode_alg     = QNDF()
        # ode_alg     = radau()
        # ode_alg     = RadauIIA5()
        # ode_alg     = DFBDF()

        dae_alg       = IDA(),

        abstol        = 1e-12,
        reltol        = 1e-12,

        dt            = 0.01)

#--------------

post_fault_steady_state_data =
    get_a_status_steady_state_data(
        :post_fault_state;
        with_faults =
            false,

        # :pre_fault_state
        # :fault_state 
        # :post_fault_state ,
        # system_status = :pre_fault_state,

        net_data_by_components_file =
            net_data_by_components_file,
        components_libs_dir =
            components_libs_dir,

        timespan         = 10.0,
        on_fault_time    = 9.0,
        clear_fault_time = 9.001,

        list_fault_point_from_node_a = [0.3],
        list_fault_resistance        = [0.001],
        list_no_line_circuit         =  [4],

        list_edges_to_have_fault   = [ 2 ],
        clear_fault_selection_list = [1],

        basekV          = 1.0,    
        use_pu_in_PQ    = true,
        line_data_in_pu = true,

        #--------------------------------------    

        use_init_u0 = false,

        use_nlsolve = false,

        pf_alg        = NewtonRaphson(),

        #--------------------------------------    

        ode_alg       = Rodas4(),
        # ode_alg     = FBDF()
        # ode_alg     = QNDF()
        # ode_alg     = radau()
        # ode_alg     = RadauIIA5()
        # ode_alg     = DFBDF()

        dae_alg       = IDA(),

        abstol        = 1e-12,
        reltol        = 1e-12,

        dt            = 0.01)

# status_steady_state_parameters =
#     get_status_steady_state_parameters(
#         net_data_by_components_file;
#         components_libs_dir =
#             components_libs_dir,
        
#         basekV = 1.0,
        
#         use_pu_in_PQ = true,
        
#         line_data_in_pu = true,

#         use_init_u0 = false,    
#         use_nlsolve = false,

#         pf_alg = NewtonRaphson(),

#         abstol =
#             abstol,
        
#         reltol =
#             reltol,

#         on_fault_time    =
#             on_fault_time,
        
#         clear_fault_time =
#             clear_fault_time,

#         list_fault_point_from_node_a =
#             list_fault_point_from_node_a,
        
#         list_fault_resistance =
#             list_fault_resistance,
        
#         list_no_line_circuit =
#             list_no_line_circuit,

#         list_edges_to_have_fault   =
#             list_edges_to_have_fault,
        
#         clear_fault_selection_list =
#             clear_fault_selection_list,

#         with_faults = false )

# #----------------------------------------

# #----------------------------------------

# ntuple_status_steady_state_data =
#     get_ntuple_status_steady_state_data(
#     ;with_faults = with_faults,
#     net_data_by_components_file =
#         net_data_by_components_file,
#     components_libs_dir =
#         components_libs_dir,
    
#         timespan =
#             timespan,
#         on_fault_time =
#             on_fault_time,
#         clear_fault_time =
#             clear_fault_time,
    
#         list_fault_point_from_node_a =
#             list_fault_point_from_node_a,
#         list_fault_resistance =
#             list_fault_resistance,
#         list_no_line_circuit =
#             list_no_line_circuit,

#         list_edges_to_have_fault =
#             list_edges_to_have_fault,
#         clear_fault_selection_list =
#             clear_fault_selection_list,
    
#         basekV =
#             1.0,    
#         use_pu_in_PQ =
#             use_pu_in_PQ,
#         line_data_in_pu =
#             line_data_in_pu)


# (;sta_pf_red_sol,
#  dyn_pf_fun_kwd_net_idxs) =
#     NamedTupleTools.select(
#         getproperty(
#             getproperty(
#                 ntuple_status_steady_state_data,
#                 :pre_fault_state),
#             :static_prefault_paras),
#         (:sta_pf_red_sol,
#          :dyn_pf_fun_kwd_net_idxs))

#---------------------------------------------------
#---------------------------------------------------
