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


sim_type  = "sim-network-graph"

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


#-------------------------------

edges_orientation =
    getproperty(
        system_net_static_data,
        :edges_orientation)

#-------------------------------

#-------------------------------
# GraphRecipes
#-------------------------------
# SimpleDiGraph
# SimpleGraph

# https://docs.juliaplots.org/stable/generated/graph_attributes/#graph_attributes
# https://docs.juliaplots.org/stable/GraphRecipes/examples/#graph_examples

using GraphRecipes

 
G_net_graph = SimpleDiGraph(
    Edge.(edges_orientation))

dict_edges_label =
    Dict( a_key => a_value
          for (a_value, a_key) in
              enumerate(edges_orientation))


# plt_G_net = graphplot(G_net_graph)

plt_G_net2 = plot(G_net_graph;
                  names=collect(1:nv(G_net_graph)),
                  edgelabel = dict_edges_label,
                  edgelabel_offset=0.0) 

#---------------------------------------------------
#---------------------------------------------------
