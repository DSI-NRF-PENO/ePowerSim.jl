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

(;plant_generators_data_from_json,
 plant_loads_data_from_json,
 plant_transmission_data_from_json,
 edge_data_from_json,
 shunt_data_from_json,
 baseMVA,
 all_nodes_idx,
 n2s_all_nodes_idx) =
     NamedTupleTools.select(
         system_net_static_data,
         (:plant_generators_data_from_json,
          :plant_loads_data_from_json,
          :plant_transmission_data_from_json,
          :edge_data_from_json,
          :shunt_data_from_json,
          :baseMVA,
          :all_nodes_idx,
          :n2s_all_nodes_idx )  )

#-------------------------------

sh_Ynet_wt_nodes_idx_wt_adjacent_nodes =
    get_Ynet_sp_sh(
        edge_data_from_json,
        shunt_data_from_json;
        baseMVA = baseMVA,
        basekV = 1.0,
        baseShunt = baseMVA,
        line_data_in_pu = true)


Ynet_wt_nodes_idx_wt_adjacent_nodes =
    get_Ynet(
        edge_data_from_json,
        shunt_data_from_json;
        baseMVA = baseMVA,
        basekV = 1.0,
        baseShunt = baseMVA,
        line_data_in_pu = true)

r_Ynet =
    round_up_Ynet(getproperty(
        Ynet_wt_nodes_idx_wt_adjacent_nodes,
        :Ynet))

# julia> r_Ynet
# 5-element Vector{Vector{ComplexF64}}:
#  [22.2507 - 222.4844im, -3.5235 + 35.2348im, -3.2569 + 32.569im, -15.4703 + 154.703im]
#  [12.6911 - 126.8979im, -3.5235 + 35.2348im, -9.1676 + 91.6758im]
#  [12.5013 - 124.9999im, -9.1676 + 91.6758im, -3.3337 + 33.3367im]
#  [9.9242 - 99.2324im, -3.2569 + 32.569im, -3.3337 + 33.3367im, -3.3337 + 33.3367im]
#  [18.804 - 188.0206im, -15.4703 + 154.703im, -3.3337 + 33.3367im]

r_nodes_idx_with_adjacent_nodes_idx =
    getproperty(
        Ynet_wt_nodes_idx_wt_adjacent_nodes,
        :nodes_idx_with_adjacent_nodes_idx)

# julia> nodes_idx_with_adjacent_nodes_idx
# 5-element Vector{Vector{Int64}}:
#  [1, 2, 4, 5]
#  [2, 1, 3]
#  [3, 2, 4]
#  [4, 1, 3, 5]
#  [5, 1, 4]

Ybus = getproperty(get_Ybus(
    edge_data_from_json,
    shunt_data_from_json;
    basekV = 1.0,
    baseMVA = baseMVA,
    line_data_in_pu = true ), :Ybus)

# julia> Ybus
# 5×5 SparseMatrixCSC{ComplexF64, Int64} with 17 stored entries:
#   2225.07-22248.4im  -352.348+3523.48im           ⋅           -325.69+3256.9im   -1547.03+15470.3im
#  -352.348+3523.48im   1269.11-12689.8im  -916.758+9167.58im           ⋅                   ⋅    
#           ⋅          -916.758+9167.58im   1250.13-12500.0im  -333.367+3333.67im           ⋅    
#   -325.69+3256.9im            ⋅          -333.367+3333.67im   992.424-9923.24im  -333.367+3333.67im
#  -1547.03+15470.3im           ⋅                   ⋅          -333.367+3333.67im    1880.4-18802.1im


#-------------------------------

(sta_red_vh_θh_0,
 pf_PQ_param,
 pf_kw_para) =
     NamedTupleTools.select(
         system_net_static_data,
         (:sta_red_vh_θh_0,
          :pf_PQ_param,
          :pf_kw_para))

generic_red_sol_kwd_para =
    getproperty(system_net_static_data,
                :generic_red_sol_kwd_para)

pf_kw_net_para =
    getproperty(
        getproperty(
            generic_red_sol_kwd_para,
            :pf_kw_para),
        :pf_kw_net_para)


#-------------------------------
# Ordinary Power flow Ynet
#-------------------------------

Ynet_pf_fun_mismatch =
    # get_ΔI_mismatch_by_Ynet
    get_ΔPQ_mismatch_by_Ynet


 Ynet_pf_sol =  NonlinearSolve.solve(
    NonlinearProblem(
        NonlinearFunction( ( g, x, p ) ->
            Ynet_pf_fun_mismatch(
                g, x, p;
                pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para ),
                Ynet_wt_nodes_idx_wt_adjacent_nodes =
                    Ynet_wt_nodes_idx_wt_adjacent_nodes)),
        sta_red_vh_θh_0,
        pf_PQ_param ),
    pf_alg ;
    abstol =
        abstol,
    reltol =
        reltol )

r_Ynet_pf_sol =
    round.(Ynet_pf_sol;digits=4)

# julia> r_Ynet_pf_sol
# 5-element Vector{Float64}:
#   0.9893
#   0.0571
#  -0.0133
#  -0.0086
#   0.0718

(;pf_P_gens, pf_Q_gens,
 pf_P_g_gens, pf_Q_g_gens,
 vh, θh, θh_deg ) = NamedTupleTools.select(
         get_results_static_pf_red_sol_u(
             Ynet_pf_sol;
             generic_red_sol_kwd_para =
                 generic_red_sol_kwd_para,
             baseMVA = baseMVA,
             basekV = 1.0) ,
             (:pf_P_gens, :pf_Q_gens,             
              :pf_P_g_gens, :pf_Q_g_gens,
             
              :vh, :θh, :θh_deg ) )


r_pf_P_gens =
    round.(pf_P_gens ;digits=4)

# julia> r_pf_P_gens
# 4-element Vector{Float64}:
#  2.1
#  3.2349
#  0.0503
#  4.6651

r_pf_Q_gens =
    round.(pf_Q_gens ;digits=4)

# julia> r_pf_Q_gens =
#            round.(pf_Q_gens ;digits=4)
# 4-element Vector{Float64}:
#   0.3073
#   1.9465
#   1.8412
#  -0.3821

r_pf_P_g_gens =
    round.(pf_P_g_gens ;digits=4)

# julia> r_pf_P_g_gens
# 4-element Vector{Float64}:
#  2.1
#  3.2349
#  0.0503
#  4.6651

r_pf_Q_g_gens =
    round.(pf_Q_g_gens ;digits=4)

# julia> r_pf_Q_g_gens
# 4-element Vector{Float64}:
#   0.3073
#   1.9465
#   1.8412
#  -0.3821

r_vh =
    round.(vh ;digits=4)

# julia> r_vh 
# 5-element Vector{Float64}:
#  1.0
#  0.9893
#  1.0
#  1.0
#  1.0


r_θh =
    round.(θh ;digits=4)

# julia> r_θh 
# 5-element Vector{Float64}:
#   0.0571
#  -0.0133
#  -0.0086
#   0.0
#   0.0718

r_θh_deg =
    round.(θh_deg ;digits=4)

# julia> r_θh_deg
# 5-element Vector{Float64}:
#   3.2734
#  -0.7593
#  -0.4923
#   0.0
#   4.112


#-------------------------------
# Ynet  NonlinearSolve with Jac
#-------------------------------

red_ΔPQ_x = similar(sta_red_vh_θh_0)


Jac_row_size =
    Jac_col_size = length( sta_red_vh_θh_0 )


Jac_vh_θh =
    spzeros( Jac_row_size, Jac_col_size )


Ynet_pf_fun_mismatch =
    get_generic_sta_pf_ΔPQ_mismatch
    # get_a_model_integrated_pf_sta_ΔPQ_mismatch


 jac_Ynet_pf_sol =  NonlinearSolve.solve(
    NonlinearProblem(
        NonlinearFunction( ( g, x, p ) ->
            Ynet_pf_fun_mismatch(
                g, x, p;
                pf_kw_para =
                    pf_kw_para),
                           jac = (Jac_vh_θh, x, p) ->
                               sta_pf_Jac!_df_dx_by_Ynet_or_Yπ_net!(
                                   Jac_vh_θh, x, p;
                                   pf_kw_para =
                                       pf_kw_para,
                                   func =
                                       Ynet_pf_fun_mismatch,
                                   net_addmitance_tuple =
                                       Ynet_wt_nodes_idx_wt_adjacent_nodes,
                                   by_Ynet_or_Yπ_net =
                                       :Ynet ) ),
        sta_red_vh_θh_0,
        pf_PQ_param ),
    NewtonRaphson() ;
    abstol =
        abstol,
    reltol =
        reltol )

#-------------------------------

Ynet_pf_fun_mismatch =
    # get_ΔI_mismatch_by_Ynet
    get_ΔPQ_mismatch_by_Ynet


 jac_Ynet_pf_sol =  NonlinearSolve.solve(
    NonlinearProblem(
        NonlinearFunction( ( g, x, p ) ->
            Ynet_pf_fun_mismatch(
                g, x, p;
                pf_kw_para =  disaggregate_sta_pf_keywords_parameter(
                                           pf_kw_para),
                Ynet_wt_nodes_idx_wt_adjacent_nodes =
                    Ynet_wt_nodes_idx_wt_adjacent_nodes),
                           jac = (Jac_vh_θh, x, p) ->
                               sta_pf_Jac!_df_dx_by_Ynet_or_Yπ_net!(
                                   Jac_vh_θh, x, p;
                                   pf_kw_para = pf_kw_para,
                                   func =
                                       Ynet_pf_fun_mismatch,
                                   net_addmitance_tuple =
                                       Ynet_wt_nodes_idx_wt_adjacent_nodes,
                                   by_Ynet_or_Yπ_net =
                                       :Ynet ) ),
        sta_red_vh_θh_0,
        pf_PQ_param ),
    NewtonRaphson() ;
    abstol =
        abstol,
    reltol =
        reltol )

#-------------------------------
# Ybus  NonlinearSolve with Jac
#-------------------------------


red_ΔPQ_x = similar(sta_red_vh_θh_0)


Jac_row_size =
    Jac_col_size = length( sta_red_vh_θh_0 )


Jac_vh_θh =
    spzeros( Jac_row_size, Jac_col_size )


Ybus_pf_fun_mismatch =
    get_ΔPQ_mismatch_by_sparse_Ybus
    # get_ΔPQ_mismatch_by_Ybus
    # get_ΔI_mismatch_by_Ybus


jac_nonlinearsolve_pf_sol = NonlinearSolve.solve(
    NonlinearProblem(
        NonlinearFunction(
            ( g, x, p ) ->
                Ybus_pf_fun_mismatch(
                g, x, p;
                pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
                    Ybus = Ybus,
                    use_autodiff =
                        false ),
        jac = (Jac_vh_θh, x, p) -> sta_pf_Jac!(
            Jac_vh_θh, x, p;
            Ybus =
                Ybus,
            pf_kw_para =
                disaggregate_sta_pf_keywords_parameter(
                    pf_kw_para ) )  ) ,
        sta_red_vh_θh_0,
        pf_PQ_param),
    NewtonRaphson();
        abstol = 1e-10,
    reltol = 1e-10)


#---------------------------------------------------
#---------------------------------------------------
