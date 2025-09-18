#=
# [Example of Aggregate Indices and Parameters with IEEE 5 Bus Test System](@id ieee-5-bus-aggregate-indices-and-parameters)
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


sim_type  = "sim-ordinary-powerflow"

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


abstol      = 1e-12

reltol      = 1e-12


#---------------------------------------------------
#---------------------------------------------------

system_net_static_data =
    get_system_net_static_data(
        case_name ;
        script_dir="",
        data_dir = "",
        json_net_data_by_components_file =
            json_net_data_by_components_file,
        components_libs_dir = "",
        basekV              = 1.0,    
        use_pu_in_PQ        = true,
        line_data_in_pu     = true,
        pf_alg              =
            NewtonRaphson(),
        no_lines_fault = 1)


#----------------------------------------

generic_system_simulation_parameters = 
    get_generic_system_simulation_parameters(
        net_data_by_components_file;
        components_libs_dir =
            components_libs_dir,
        basekV = 1.0,    
        use_pu_in_PQ      = true,
        opf_use_pu_in_PQ  = true,
        line_data_in_pu   = true,
        with_faults       = false,
        pf_alg            = NewtonRaphson(),
        ode_alg           = Rodas4(),
        dae_alg           = IDA() )

#----------------------------------------

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

#----------------------------------------

status_steady_state_parameters =
    get_status_steady_state_parameters(
        net_data_by_components_file;
        components_libs_dir =
            components_libs_dir,
        
        basekV = 1.0,
        
        use_pu_in_PQ = true,
        
        line_data_in_pu = true,

        use_init_u0 = false,    
        use_nlsolve = false,

        pf_alg = NewtonRaphson(),

        abstol =
            abstol,
        
        reltol =
            reltol,

        on_fault_time    =
            on_fault_time,
        
        clear_fault_time =
            clear_fault_time,

        list_fault_point_from_node_a =
            list_fault_point_from_node_a,
        
        list_fault_resistance =
            list_fault_resistance,
        
        list_no_line_circuit =
            list_no_line_circuit,

        list_edges_to_have_fault   =
            list_edges_to_have_fault,
        
        clear_fault_selection_list =
            clear_fault_selection_list,

        with_faults = false )

#----------------------------------------

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

#----------------------------------------

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
            1.0,    
        use_pu_in_PQ =
            use_pu_in_PQ,
        line_data_in_pu =
            line_data_in_pu)


#---------------------------------------------------
#---------------------------------------------------
