#=
# [Example of Distributed Slack Power Flow with IEEE 5 Bus Test System] (@id ieee-5-bus-distributed-slack-power-flow)
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

sauer_net_wt_avr_string    = "net-static-data-avr-sauer-"

rtds_net_wt_avr_string     = "net-static-data-avr-rtds-"

sauer_gov_string           = "gov-sauer"

ieee_tgov2sauer_gov_string = "gov-ieee-tgov2sauer"

ieee_tgov1_gov_string      = "gov-ieee-tgov1"

#---------------------------------------------------

"""
net_wt_avr_string = sauer_net_wt_avr_string
gov_string        = sauer_gov_string

net_wt_avr_string = rtds_net_wt_avr_string
gov_string        = ieee_tgov1_gov_string

dynamic_net_data_by_components_file =
    "$(net_wt_avr_string)"*
    "$(gov_string)" *
    ".json"

json_net_data_by_components_file =
    dynamic_net_data_by_components_file

"""


json_net_data_by_components_file =
    "dynamic-net-data-gen-sauer-avr-t1-cb-sauer-gov-t1-cb-sauer.json"

#---------------------------------------------------

case_name = "case5"

# case_name = "case9"

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

# sim_type  = "$(case_name)-"*"network-pertubation-"*
#     "-$(net_wt_avr_string)-$(gov_string)"

sim_type  = "$(case_name)-"*"distributed"


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

sim_red_model_distributed_slack_pf(
    net_data_by_components_file;
    components_libs_dir = "",
    basekV            = 1.0,    
    use_pu_in_PQ      = true,
    line_data_in_pu   = true,
    with_faults       = false,
    wt_network_loss_by_sta_pf_PQ_para_bool = true,
    fractional_digits = 6,
    pf_alg            = NewtonRaphson(),
    abstol      = 1e-12,

    reltol      = 1e-12)

#---------------------------------------------------
#---------------------------------------------------


(ds_slack_value,
 df_gens,
 df_non_gens,
 df_all_nodes_result) =
     NamedTupleTools.select(
         red_model_distributed_slack_pf,
         (:ds_slack_value,
          :df_gens,
          :df_non_gens,
          :df_all_nodes_result))


#---------------------------------------------------
#---------------------------------------------------
