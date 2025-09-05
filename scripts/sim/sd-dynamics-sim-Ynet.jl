# ####################################################
# using Pkg
# PowerSystemsCoSim_folder = joinpath(@__DIR__,"../..")
# cd(PowerSystemsCoSim_folder)
# Pkg.activate(PowerSystemsCoSim_folder)

#---------------------------------------------------
# global settings
#---------------------------------------------------

freq = 60

Ωb = 2 * pi * freq

ωs = Ωb 

ω_ref0 = ωs

basekV = 1.0

#---------------------------------------------------


using Pkg

using ePowerSim


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

# ode_alg   = ImplicitMidpoint()

dt          = 0.0001

Δt          = 1.0 / 2^(4)

dae_alg     = IDA()

abstol      = 1e-12

reltol      = 1e-12

#---------------------------------------------------
## fault time and data
#---------------------------------------------------

on_fault_time = 10.617 # 10.0

Δt_clear_time =
    0.5

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

pertubation_factor = 1.10

restoration_factor = 1.0

pertubation_time   = on_fault_time

restoration_time   = clear_fault_time

Δt1 = 1.5

Δt2 = 1.5

#---------------------------------------------------

list_fault_point_from_node_a = [0.3]

list_fault_resistance        = [0.001]

list_no_line_circuit         = [1]

list_edges_to_have_fault     = [ 3 ] # [ 4 ]

clear_fault_selection_list   = [1]

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


sim_type  = "sim-Ynet"


#---------------------------------------------------

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

#---------------------------------------------------


tex_filename =
    joinpath(results_dir,
                 "$(sim_type)-" *
                 "Ynet.tex")

#---------------------------------------------------

benchmark_filename =
    joinpath(results_dir,
             "$(sim_type).json")

#---------------------------------------------------


# https://github.com/JuliaTeX/TikzPictures.jl

#---------------------------------------------------

sd_dynamics_sim_csv_filename =
    joinpath(
        results_dir,
            "$(sim_type)-" *
            "Ynet.csv")

benchmark_time_ΔI_mismatch_csv_filename =
    joinpath(
        results_dir,
            "$(sim_type)-" *
            "benchmark-pf-mismatch-by-I-time.csv")

benchmark_time_ΔPQ_mismatch_csv_filename =
    joinpath(
        results_dir,
            "$(sim_type)-" *
            "benchmark-pf-mismatch-by-PQ-time.csv")

benchmark_ΔI_mismatch_csv_filename =
    joinpath(
        results_dir,
            "$(sim_type)-" *
            "benchmark-pf-mismatch-by-I-time-log10-ns.csv")

benchmark_ΔPQ_mismatch_csv_filename =
    joinpath(
        results_dir,
            "$(sim_type)-" *
            "benchmark-pf-mismatch-by-PQ-time-log10-ns.csv")

#----------------------------------------    
#----------------------------------------    

static_net_data_by_components_file =
    "opf-pf-net-default-static-data.json"

#----------------------------------------    

net_type = "dynamic-net-data"

gen_type = "gen-sauer"

gov_type = "gov-t1-cb-sauer"

avr_type = "avr-t1-cb-sauer"

data_ext = "json"

dynamic_net_data_by_components_file =
    "$(net_type)-"*
    "$(gen_type)-"*
    "$(avr_type)-"*
    "$(gov_type)"*
    ".$(data_ext)"

#---------------------------------------- 
#----------------------------------------    
# Benchmark
#----------------------------------------
#---------------------------------------- 

# using BenchmarkTools, BenchmarkPlots

# using StatsPlots

# using DataFrames

# using Statistics 

#---------------------------------------- 

SUITE = BenchmarkGroup()

SUITE[:ΔI_mismatch] =
    BenchmarkGroup(["ΔI mismatch" ])


SUITE[:ΔPQ_mismatch] =
    BenchmarkGroup(["ΔPQ mismatch" ])


SUITE[:network_admittance] =
    BenchmarkGroup(["network admittances" ])

#---------------------------------------- 

network_cases =
    ["case4gs",
     "case5",
     "case9",
     "case14",
     "case30",
     "case39",
     "case57",
     "case118",
     "case300",
     "case1354pegase",
     "case2869pegase",
     # "case9241pegase",
     # "case13659pegase"
     
     ]

#---------------------------------------- 

nt_admittance_vector_funcs =
    (;get_Ynet,
     get_Ynet_sp_sh,
     get_Yπ_net,
     get_Ybus ) 

#---------------------------------------- 

mismatch_funcs_types =
    (;get_ΔI_mismatch_by_Ynet,
     # get_ΔI_mismatch_by_Yπ_net,
     get_ΔI_mismatch_by_Ybus,

     get_ΔPQ_mismatch_by_Ynet,
     # get_ΔPQ_mismatch_by_Yπ_net,
     get_ΔPQ_mismatch_by_Ybus,
     get_ΔPQ_mismatch_by_sparse_Ybus)

#---------------------------------------- 

create_benchmarkgroups_admittance_vector_funcs(
    SUITE, network_cases,
    nt_admittance_vector_funcs  )

create_benchmarkgroups_pf_mismatch(
    SUITE, network_cases,
    mismatch_funcs_types  )

# clear_empty!(SUITE) 

results = run(SUITE, verbose = true, seconds = 450 )

BenchmarkTools.save(benchmark_filename, results)

# t_plt = plot( results["network_admittance"] )

#---------------------------------------- 
#---------------------------------------- 

results[:ΔI_mismatch]

results[:ΔPQ_mismatch]

results[:network_admittance]

#---------------------------------------- 
#---------------------------------------- 

plot_benchmark_results(
    results;
    bm_group_sym =
        :ΔI_mismatch,
    
    fun_string_Ybus =
        "get_ΔI_mismatch_by_Ybus",
    fun_string_Ynet =
        "get_ΔI_mismatch_by_Ynet",
    
    network_cases =
        ["case4gs", "case5", "case9",
         "case14",  "case30", "case39",
         "case57", "case118", "case300",
         "case1354pegase", "case2869pegase"],
    
    net_cases_abr =
        ["c4", "c5", "c9", "c14", "c30",
         "c39", "c57", "c118", "c300", "c1354", "c2869"],
    x_coords =
        [4, 5, 9, 14, 30, 39, 57, 118, 300, 1354, 2869],
    
    ctg_string_by_Ybus =
        "ΔI mismatch by Ybus",
    ctg_string_by_Ynet =
        "ΔI mismatch by Ynet",

    benchmark_Δ_mismatch_csv_filename =
        benchmark_ΔI_mismatch_csv_filename,
        # benchmark_ΔPQ_mismatch_csv_filename

    benchmark_time_Δ_mismatch_csv_filename =
        benchmark_time_ΔI_mismatch_csv_filename,
    
    groupedbar_filename_string =
        "ΔI-mismatch-Ynet-Ybus-groupedbar.pdf",
    scatter_filename_string =
        "ΔI-mismatch-Ynet-Ybus-scatter.pdf",
    
    lable_scatter_mismatch_by_Ybus =
        "ΔI mismatch by Ybus",
    lable_scatter_mismatch_by_Ynet =
        "ΔI mismatch by Ynet" )



plot_benchmark_results(
    results;
    bm_group_sym =
        :ΔPQ_mismatch,
    
    fun_string_Ybus =
        "get_ΔPQ_mismatch_by_Ybus",
    fun_string_Ynet =
        "get_ΔPQ_mismatch_by_Ynet",
    
    network_cases =
        ["case4gs", "case5", "case9",
         "case14",  "case30", "case39",
         "case57", "case118", "case300",
         "case1354pegase", "case2869pegase"],
    
    net_cases_abr =
        ["c4", "c5", "c9", "c14", "c30",
         "c39", "c57", "c118", "c300", "c1354", "c2869"],
    x_coords =
        [4, 5, 9, 14, 30, 39, 57, 118, 300, 1354, 2869],
    
    ctg_string_by_Ybus =
        "ΔPQ mismatch by Ybus",
    ctg_string_by_Ynet =
        "ΔPQ mismatch by Ynet",

    benchmark_Δ_mismatch_csv_filename =
        benchmark_ΔPQ_mismatch_csv_filename,

    benchmark_time_Δ_mismatch_csv_filename =
        benchmark_time_ΔPQ_mismatch_csv_filename,
    
    groupedbar_filename_string =
        "ΔPQ-mismatch-Ynet-Ybus-groupedbar.pdf",
    scatter_filename_string =
        "ΔPQ-mismatch-Ynet-Ybus-scatter.pdf",
    
    lable_scatter_mismatch_by_Ybus =
        "ΔPQ mismatch by Ybus",
    lable_scatter_mismatch_by_Ynet =
        "ΔPQ mismatch by Ynet" )



# """

# df_ΔI_mismatch = DataFrame( results[:ΔI_mismatch] )

# tran_df_ΔI_mismatch =
#     transform(
#         df_ΔI_mismatch,
#         :first => ByRow(identity) => [:ΔI_mismatch, :case],
#         :second => (t -> log10.(time.(t))) => :Time )

# st_df_ΔI_mismatch =
#     DataFrames.select(
#         tran_df_ΔI_mismatch,
#         :ΔI_mismatch,
#         :case,
#         :Time)

# # Select rows where column :ΔI_mismatch == "get_ΔI_mismatch_by_Ybus"

# s_Ybus_ΔI_mismatch_by_ΔI_mismatch_cases =
#     subset(st_df_ΔI_mismatch,
#            :ΔI_mismatch => ByRow(x -> x ==
#                "get_ΔI_mismatch_by_Ybus" ))


# s_Ynet_ΔI_mismatch_by_ΔI_mismatch_cases =
#     subset(st_df_ΔI_mismatch,
#            :ΔI_mismatch => ByRow(x -> x ==
#                "get_ΔI_mismatch_by_Ynet" ))

# #--------------------

# cases_sort_order =
#     network_cases

# Ybus_sort_df =
#     s_Ybus_ΔI_mismatch_by_ΔI_mismatch_cases[
#         indexin(cases_sort_order,
#                 s_Ybus_ΔI_mismatch_by_ΔI_mismatch_cases.case),:]

# Ynet_sort_df =
#     s_Ynet_ΔI_mismatch_by_ΔI_mismatch_cases[
#         indexin(cases_sort_order,
#                 s_Ynet_ΔI_mismatch_by_ΔI_mismatch_cases.case),:]


# net_cases =
#     ["c4", "c5", "c9", "c14", "c30",
#      "c39", "c57", "c118", "c300", "c1354", "c2869"]

# x_coords = [4, 5, 9, 14, 30, 39, 57, 118, 300, 1354, 2869]

# Ybus_t_values = Ybus_sort_df.Time

# Ynet_t_values = Ynet_sort_df.Time

# labels = net_cases

# # https://github.com/JuliaPlots/StatsPlots.jl#grouped-bar-plots

# ctg = repeat(["ΔI mismatch by Ybus", "ΔI mismatch by Ynet"], inner = 11)

# nam = repeat(net_cases, outer = 2)


# figure_groupedbar_filename =
#     joinpath(figure_dir,"ΔI-mismatch-Ynet-Ybus-groupedbar.pdf")

# figure_scatter_filename =
#     joinpath(figure_dir,"ΔI-mismatch-Ynet-Ybus-scatter.pdf")

# ΔI_mismatch_groupedbar =
#     StatsPlots.groupedbar(
#         nam, [Ybus_t_values Ynet_t_values],
#         group = ctg,
#         xlabel = "Cases",
#         ylabel = "Time log10([ns])",
#         title = "ΔI mismatch by admittances",
#         bar_width = 0.67,
#         lw = 0,
#         framestyle = :box)

# plt_scatter_ΔI_mismatch =
#     scatter(x_coords, Ybus_t_values, 
#         # series_annotations=text.(labels, :bottom),
#         # Annotate points with string labels
#         xlabel="cases [number of buses]", 
#         ylabel="Time log10([ns])", 
#             label="ΔI mismatch by Ybus")

# scatter!(plt_scatter_ΔI_mismatch,
#          x_coords, Ynet_t_values,
#          label="ΔI mismatch by Ynet")


# # plt_bar_ΔI_mismatch_by_Ybus =
# #     bar(net_cases, Ybus_t_values, 
# #         xlabel="Cases [number of buses]", 
# #         ylabel="Time log10([ns])",
# #         title="ΔI mismatch by Ybus")

# #--------------------

# # Select rows based on multiple columns

# m_ΔI_mismatch_by_ΔI_mismatch_cases =
#     subset(
#         st_df_ΔI_mismatch,
#         AsTable([:ΔI_mismatch, :case]) =>
#             ByRow(x -> x.ΔI_mismatch ==
#             "get_ΔI_mismatch_by_Ybus" && x.Time > 0)) 
# """



#---------------------------------------- 
#----------------------------------------    
# Static power flow
#----------------------------------------
#---------------------------------------- 

case_name = "case14"

# case_name = "case300" 

# case_name = "case1354pegase"

# case_name = "case13659pegase"

json_case_dir =
    joinpath(
        data_dir,
        case_name,
        "json" )

# components_libs_dir

json_net_data_by_components_file =
    static_net_data_by_components_file


if json_net_data_by_components_file == ""

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

# #-------------------------------
# #-------------------------------

# (;
# plant_generators_data_from_json,
# plant_loads_data_from_json,
# plant_transmission_data_from_json,
 
# edge_data_from_json,
# shunt_data_from_json,
# baseMVA_data_from_json,
# gencost_data_from_json) =
#    NamedTupleTools.select(
#        get_net_data_by_components_from_json_file(
#            net_data_by_components_file;
#            in_components_type_sym =
#                false ),
#        (:plant_generators_data_from_json,
#         :plant_loads_data_from_json,
#         :plant_transmission_data_from_json,
        
#         :edge_data_from_json,
#         :shunt_data_from_json,
#         :baseMVA_data_from_json,
#         :gencost_data_from_json ))

# baseMVA =
#     baseMVA_data_from_json

#-------------------------------
#-------------------------------

system_net_static_data =
    get_system_net_static_data(
        case_name ;
        script_dir=
            script_dir,
        data_dir = "",
        json_net_data_by_components_file =
            static_net_data_by_components_file,
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
#-------------------------------


sh_Ynet_wt_nodes_idx_wt_adjacent_nodes =
    get_Ynet_by_trimed_Shunt(
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


Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
    get_Yπ_net(
        edge_data_from_json,
        shunt_data_from_json;
        baseMVA = baseMVA,
        basekV = 1.0,
        baseShunt = baseMVA,
        line_data_in_pu = true,
        orientated_bool = true )


Ybus = getproperty(get_Ybus(
    edge_data_from_json,
    shunt_data_from_json;
    basekV = 1.0,
    baseMVA = baseMVA,
    line_data_in_pu = true ), :Ybus)


Ynet_wt_nodes_idx_wt_adjacent_nodes_diff =
    get_Ynet_wt_nodes_idx_wt_adjacent_nodes_diff(
        sh_Ynet_wt_nodes_idx_wt_adjacent_nodes,
        Ynet_wt_nodes_idx_wt_adjacent_nodes)

#-------------------------------
# size of net admittance vector
#-------------------------------

# using HypothesisTests


# (data_memory_size_Ybus,
#  idx_memory_size_Ybus,
#  data_to_idx_ratio_Ybus) =
#     get_size_Ybus(Ybus)


# (data_memory_size_Ynet,
#  idx_memory_size_Ynet,
#  data_to_idx_ratio_Ynet) =
#     get_size_Ynet_wt_nodes_idx_wt_adjacent_nodes(
#         Ynet_wt_nodes_idx_wt_adjacent_nodes)


# (data_memory_size_Yπ_net,
#  idx_memory_size_Yπ_net,
#  data_to_idx_ratio_Yπ_net) =
#     get_size_Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes(
#         Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes)


#-------------------------------
# Benchark for net admittance vectors
#-------------------------------

bm_Ybus = @benchmarkable get_Ybus(
    $edge_data_from_json,
    $shunt_data_from_json;
    basekV = 1.0,
    baseMVA = $baseMVA,
    line_data_in_pu = true )

bm_sh_Ynet = @benchmarkable  get_Ynet_sp_sh(
        $edge_data_from_json,
        $shunt_data_from_json;
        baseMVA = $baseMVA,
        basekV = 1.0,
        baseShunt = $baseMVA,
        line_data_in_pu = true)


bm_Ynet = @benchmarkable  get_Ynet(
        $edge_data_from_json,
        $shunt_data_from_json;
        baseMVA = $baseMVA,
        basekV = 1.0,
        baseShunt = $baseMVA,
        line_data_in_pu = true)

tune!(bm_Ybus)

tune!(bm_sh_Ynet)

tune!(bm_Ynet)

r_Ybus    = median(run(bm_Ybus ))

r_sh_Ynet = median(run(bm_sh_Ynet))

r_Ynet    = median(run(bm_Ynet))

ratio(r_Ybus, r_Ynet)

judge(r_Ybus, r_Ynet)

ratio(r_sh_Ynet, r_Ynet)

judge(r_sh_Ynet, r_Ynet)


#-------------------------------
#-------------------------------
# Power flow 
#-------------------------------
#-------------------------------


# ΔI_mismatch_funcs_types =
#     (;get_ΔI_mismatch_by_Ynet,
#      get_ΔI_mismatch_by_Yπ_net,
#      get_ΔI_mismatch_by_Ybus )


# ΔPQ_mismatch_funcs_types =
#     (;get_ΔPQ_mismatch_by_Ynet,
#      get_ΔPQ_mismatch_by_Yπ_net,
#      get_ΔPQ_mismatch_by_Ybus,
#      get_ΔPQ_mismatch_by_sparse_Ybus )


#-------------------------------


(sta_red_vh_θh_0,
 pf_PQ_param,
 pf_kw_para) =
     NamedTupleTools.select(
         system_net_static_data,
         (:sta_red_vh_θh_0,
          :pf_PQ_param,
          :pf_kw_para))


#-------------------------------
# Power flow Ynet
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
                    Ynet_wt_nodes_idx_wt_adjacent_nodes ) ),
        sta_red_vh_θh_0,
        pf_PQ_param ),
    pf_alg ;
    abstol =
        abstol,
    reltol =
        reltol )

#-------------------------------
#Power flow Yπ_net
#-------------------------------

Yπ_net_pf_fun_mismatch =
    get_ΔPQ_mismatch_by_S_inj_Yπ_net

    # get_ΔI_mismatch_by_Yπ_net
    # get_ΔPQ_mismatch_by_Yπ_net

 Yπ_net_pf_sol =  NonlinearSolve.solve(
    NonlinearProblem(
        NonlinearFunction( ( g, x, p ) ->
            Yπ_net_pf_fun_mismatch(
                g, x, p;
                pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
            Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
        Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes ) ),
        sta_red_vh_θh_0,
        pf_PQ_param ),
    pf_alg ;
    abstol =
        abstol,
    reltol =
        reltol )

#-------------------------------
# Power flow Ybus
#-------------------------------

Ybus_pf_fun_mismatch =
    get_ΔPQ_mismatch_by_Ybus
    # get_ΔPQ_mismatch_by_sparse_Ybus
    # get_ΔI_mismatch_by_Ybus


 Ybus_pf_sol =  NonlinearSolve.solve(
    NonlinearProblem(
        NonlinearFunction( ( g, x, p ) ->
            Ybus_pf_fun_mismatch(
                g, x, p;
                pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
                Ybus = Ybus,
            use_autodiff = true ) ),
        sta_red_vh_θh_0,
        pf_PQ_param ),
    pf_alg ;
    abstol =
        abstol,
    reltol =
        reltol )


#-------------------------------
# NonlinearSolve autodiff
#-------------------------------

pf_fun_mismatch =
    get_ΔPQ_mismatch_by_sparse_Ybus


autodiff_nonlinearsolve_pf_sol = NonlinearSolve.solve(
    NonlinearProblem(
        NonlinearFunction(
            ( g, x, p ) ->
                pf_fun_mismatch(
                g, x, p;
                pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
                Ybus = Ybus) ) ,
        sta_red_vh_θh_0,
        pf_PQ_param),
    NewtonRaphson();
        abstol = 1e-10,
        reltol = 1e-10)

#-------------------------------
# nlsolve autodiff
#-------------------------------

pf_fun_mismatch =
    get_ΔPQ_mismatch_by_sparse_Ybus


autodiff_nlsolve_pf_sol = nlsolve((g, x) ->
    pf_fun_mismatch(
        g, x, pf_PQ_param;
        Ybus =
            Ybus,
        pf_kw_para =
            disaggregate_sta_pf_keywords_parameter(
                pf_kw_para )),
        sta_red_vh_θh_0,
        method=:trust_region,
        autodiff=:forward )


#-------------------------------
#-------------------------------
# External Jac through NonlinearFunction
#-------------------------------
#-------------------------------

# https://docs.sciml.ai/NonlinearSolve/stable/basics/nonlinear_functions/

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
# Yπ_net  NonlinearSolve with Jac
#-------------------------------

Yπ_net_pf_fun_mismatch =
    get_ΔPQ_mismatch_by_Yπ_net
    # get_ΔI_mismatch_by_Yπ_net
    # get_ΔPQ_mismatch_by_S_inj_Yπ_net

 jac_Ynet_pf_sol =  NonlinearSolve.solve(
    NonlinearProblem(
        NonlinearFunction( ( g, x, p ) ->
            Yπ_net_pf_fun_mismatch(
                g, x, p;
                pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
                Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes  =
                    Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes ),
                           jac = (Jac_vh_θh, x, p) ->
                               sta_pf_Jac!_df_dx_by_Ynet_or_Yπ_net!(
                                   Jac_vh_θh, x, p;
                                   pf_kw_para = pf_kw_para,
                                   func =
                                       Yπ_net_pf_fun_mismatch,
                                net_addmitance_tuple =
                                   Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
                                   by_Ynet_or_Yπ_net =
                                       :Yπ_net ) ),
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
 
#-------------------------------
# Nlsolve with Jac
#-------------------------------

red_ΔPQ_x = similar(sta_red_vh_θh_0)


Jac_row_size =
    Jac_col_size = length( sta_red_vh_θh_0 )


Jac_vh_θh =
    spzeros( Jac_row_size, Jac_col_size )


Ybus_pf_fun_mismatch =
    get_ΔPQ_mismatch_by_sparse_Ybus


Ybus_jac_nlsolve_pf_sol = nlsolve(
    (red_ΔPQ_x, red) ->
        Ybus_pf_fun_mismatch(
            red_ΔPQ_x,
            red,
            pf_PQ_param;
            Ybus =
                Ybus,
            pf_kw_para =
                disaggregate_sta_pf_keywords_parameter(
                    pf_kw_para),
            use_autodiff =
                false),
    (Jac_vh_θh, red )->
        sta_pf_Jac!(
            Jac_vh_θh,
            red,
            pf_PQ_param;            
            Ybus=
                Ybus,
            pf_kw_para =
                disaggregate_sta_pf_keywords_parameter(
                    pf_kw_para) ) ,
    sta_red_vh_θh_0;
    ftol = 1e-10 )


#-------------------------------
#-------------------------------
# Jacobian-Free Newton-Krylov
#-------------------------------
#-------------------------------


#-------------------------------
# Ybus Using Jacobian-Free Newton-Krylov
#-------------------------------

using SparseConnectivityTracer

Ybus_pf_fun_mismatch =
    get_ΔPQ_mismatch_by_sparse_Ybus

 KrylovJL_GMRES_pf_sol =  NonlinearSolve.solve(
    NonlinearProblem(
        NonlinearFunction( ( g, x, p ) ->
            Ybus_pf_fun_mismatch(
                g, x, p;
                pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
                Ybus =
                    Ybus) ;sparsity = TracerSparsityDetector()  ),
        sta_red_vh_θh_0,
        pf_PQ_param ),
    NewtonRaphson(linsolve = KrylovJL_GMRES() ) ;
    abstol =
        abstol,
    reltol =
        reltol )

#-------------------------------

nlsolve_sol_properties =
    (:method, :initial_x, :zero, :residual_norm,
     :iterations, :x_converged, :xtol, :f_converged,
     :ftol, :trace, :f_calls, :g_calls )

nonlinearsolve_sol_properties =
    (:u, :resid, :prob, :alg, :retcode,
     :original, :left, :right, :stats, :trace )

jac_nlsolve_pf_sol.method

autodiff_nlsolve_pf_sol.method

jac_nonlinearsolve_pf_sol.stats

autodiff_nonlinearsolve_pf_sol.stats




# """


# (;
# plant_generators_data_from_json,
# plant_loads_data_from_json,
# plant_transmission_data_from_json,
 
# edge_data_from_json,
# shunt_data_from_json,
# baseMVA_data_from_json,
# gencost_data_from_json) =
#    NamedTupleTools.select(
#        get_net_data_by_components_from_json_file(
#            net_data_by_components_file;
#            in_components_type_sym =
#                false ),
#        (:plant_generators_data_from_json,
#         :plant_loads_data_from_json,
#         :plant_transmission_data_from_json,
        
#         :edge_data_from_json,
#         :shunt_data_from_json,
#         :baseMVA_data_from_json,
#         :gencost_data_from_json ))

# baseMVA =
#     baseMVA_data_from_json


# #-------------------------------
# #-------------------------------

# # (;edges_fbus, edges_tbus) =
# #     get_edges_nt_fbus_tbus_by_json(
# #         edge_data_from_json)

# # edges_orientation =
# #     get_edges_orientation_by_generic(
# #     edges_fbus, edges_tbus)


# t_edges_orientations =
#     get_edges_orientation_by_generic(
#     values(get_edges_nt_fbus_tbus_by_json(
#         edge_data_from_json))... )


# t_Cnb = get_Cnb_by_orientations(
#         t_edges_orientations)


# t_nodes_incident_edges =
#     get_nodes_incident_edges_by_orientations(
#         t_edges_orientations )


# #-------------------------------


# get_Ybus(
#     edge_data_from_json,
#     shunt_data_from_json;
#     basekV = 1.0,
#     baseMVA = baseMVA,
#     line_data_in_pu = true )


#  get_Ynet(
#     edge_data_from_json,
#     shunt_data_from_json;
#     baseMVA = baseMVA,
#     basekV = 1.0,
#     baseShunt = baseMVA,
#     line_data_in_pu = true)


# get_Yπ_net(
#     edge_data_from_json,
#     shunt_data_from_json;
#     baseMVA = baseMVA,
#     basekV = 1.0,
#     baseShunt = baseMVA,
#     line_data_in_pu = true,
#     orientated_bool = false )


# net_nodes_type_idxs =
#    get_net_nodes_type_idxs_by_json(
#        plant_generators_data_from_json,
#        plant_loads_data_from_json,
#        plant_transmission_data_from_json )


# # dyn_pf_fun_kwd_net_idxs =
# #    NamedTupleTools.select(
# #        net_nodes_type_idxs,
# #        (:slack_gens_nodes_idx,
# #         :non_slack_gens_nodes_idx,
# #         :gens_nodes_idx,
# #         :non_gens_nodes_idx,
# #         :gens_with_loc_load_idx,
# #         :gens_nodes_with_loc_loads_idx,
# #         :all_nodes_idx))


# # dyn_pf_fun_kwd_n2s_idxs =
# #    NamedTupleTools.select(
# #        get_dict_net_streamlined_idx_by_nodes_type_idxs(
# #            net_nodes_type_idxs ),
# #        (:n2s_slack_gens_idx,
# #         :n2s_non_slack_gens_idx,
# #         :n2s_gens_idx,
# #         :n2s_non_gens_idx,
# #         :n2s_gens_with_loc_load_idxs,
# #         :n2s_all_nodes_idx))


# all_nodes_idx =
#     getproperty(
#         net_nodes_type_idxs,
#         :all_nodes_idx)

# n2s_all_nodes_idx =
#     get_n2s_any(all_nodes_idx) 


# #-------------------------------

# gens_para_sequence_order =
#    (:components_data,
#     :gen)

# gens_generic_selections =
#    (
#     :Sn,:vh,
#     :P, :Q,
#     :Pmin, :Pmax,
#     :Qmin, :Qmax,
#     :vmin, :vmax )

# #----------------------------------------

# generic_each_gen_para =
#     get_components_properties_by_json(
#         plant_generators_data_from_json;
#         sequence_order =
#              gens_para_sequence_order,
#          selections =
#               gens_generic_selections )

# "To make sure Vector{NamedTuple} is returned
#  instead of Vector{Any}"
# generic_each_gen_para =
#     NamedTuple[
#         item for item in
#             generic_each_gen_para]

# #-------------------------------

# get_selected_vec_nt_to_vec_vec(
#    generic_each_gen_para,
#    nothing;
#     selections =
#         gens_generic_selections,
#     vec_datatype = Float64 )


# ode_gens_generic_para =
#      get_selected_comps_ode_para_by_json(
#          plant_generators_data_from_json;
#          sequence_order =
#              gens_para_sequence_order,
#          selections =
#              gens_generic_selections )

# #-------------------------------

# generic_gens_para =
#     get_ode_gens_generic_para(
#         plant_generators_data_from_json;
#         sequence_order =
#             gens_para_sequence_order,
#         selections =
#             gens_generic_selections)


# system_net_static_data =
#     get_system_net_static_data(
#         case_name;
#         script_dir=
#             script_dir,
#         data_dir = "",
#         json_net_data_by_components_file =
#             static_net_data_by_components_file,
#         components_libs_dir = "",
#         basekV              = 1.0,    
#         use_pu_in_PQ        = true,
#         line_data_in_pu     = true,
#         pf_alg              =
#             NewtonRaphson(),
#         no_lines_fault = 1)

# #-------------------------------

# (;edge_data_from_json,
#  shunt_data_from_json,
#  baseMVA,
#  all_nodes_idx,
#  n2s_all_nodes_idx) =
#      NamedTupleTools.select(
#          system_net_static_data,
#          (:edge_data_from_json,
#           :shunt_data_from_json,
#           :baseMVA,
#           :all_nodes_idx,
#           :n2s_all_nodes_idx )  )


# #-------------------------------
# #-------------------------------

# # Ynet_wt_nodes_idx_wt_adjacent_nodes =
# #     getproperty(
# #         system_net_static_data,
# #         :Ynet_wt_nodes_idx_wt_adjacent_nodes )


# edges_orientation =
#     getproperty(
#         system_net_static_data,
#         :edges_orientation)

# # (Ynet,
# #  nodes_idx_with_adjacent_nodes_idx ) =
# #      NamedTupleTools.select(
# #          getproperty(
# #              system_net_static_data,
# #              :Ynet_wt_nodes_idx_wt_adjacent_nodes ),
# #          (:Ynet,
# #           :nodes_idx_with_adjacent_nodes_idx))


# gens_nodes_idx =
#     getproperty(getproperty(
#         system_net_static_data,
#         :dyn_pf_fun_kwd_net_idxs ),
#                 :gens_nodes_idx)

# all_nodes_idx =
#     getproperty(getproperty(
#         system_net_static_data,
#         :dyn_pf_fun_kwd_net_idxs ),
#                 :all_nodes_idx)

# n2s_all_nodes_idx =
#     getproperty(getproperty(
#         system_net_static_data,
#         :dyn_pf_fun_kwd_n2s_idxs ),
#                 :n2s_all_nodes_idx)

# #-------------------------------
# #-------------------------------



# list_faulted_line_a_b_orientation =
#    edges_orientation[
#        list_edges_to_have_fault ] 


# #----------------------------------------
# #----------------------------------------


# # on_fault_net_para =
# #     make_lines_faults_data_set(
# #         Ynet,
# #         nodes_idx_with_adjacent_nodes_idx ,
# #         all_nodes_idx,
# #         n2s_all_nodes_idx,

# #         list_faulted_line_a_b_orientation ,
# #         list_fault_point_from_node_a,
# #         list_fault_resistance,
# #         list_no_line_circuit )


# #  I am splatting the values of namedtuple here
# on_fault_net_para =
#     make_lines_faults_data_set(
#         values(NamedTupleTools.select(
#          getproperty(
#              system_net_static_data,
#              :Ynet_wt_nodes_idx_wt_adjacent_nodes ),
#          (:Ynet,
#           :nodes_idx_with_adjacent_nodes_idx)))... ,

#         all_nodes_idx,
#         n2s_all_nodes_idx,

#         list_faulted_line_a_b_orientation ,
#         list_fault_point_from_node_a,
#         list_fault_resistance,
#         list_no_line_circuit )

# (;faulty_Ynet,
#  faulty_nodes_idx_with_adjacent_nodes_idx,) =
#      NamedTupleTools.select(
#          on_fault_net_para,
#          (:faulty_Ynet,
#           :faulty_nodes_idx_with_adjacent_nodes_idx))



# # (:faulty_Ynet,
# #  :faulty_nodes_idx_with_adjacent_nodes_idx,
# #  :faulty_all_nodes_idx,
# #  :n2s_faulty_all_nodes_idx,
# #  :fault_nodes_idx,
# #  :n2s_fault_nodes_idx,
# #  :list_Ya_nkf,
# #  :list_Ynkf_b,
# #  :list_faulty_line_Yl,
# #  :list_healthy_lines_Yl,
# #  :list_node_b_idx_in_a_node_row,
# #  :list_node_a_idx_in_b_node_row,
# #  :list_faulted_line_a_b_orientation,
# #  :list_fault_point_from_node_a,
# #  :list_fault_resistance,
# #  :list_no_line_circuit)

# #----------------------------------------
# # write to tex
# #----------------------------------------


# rounded_Ynet = round_up_Ynet(
#     Ynet;
#     fractional_digits=4 )

# rounded_faulty_Ynet = round_up_Ynet(
#     faulty_Ynet;
#     fractional_digits=4 )

# object_for_tex =
#     (;rounded_Ynet,
#      nodes_idx_with_adjacent_nodes_idx,
#      rounded_faulty_Ynet,
#      faulty_nodes_idx_with_adjacent_nodes_idx)

# write_vector_or_matrix_to_tex(
#     object_for_tex,
#     tex_filename)


# #----------------------------------------
# # clear fault
# #----------------------------------------


# cleared_selected_lines_faults_net_para =
#    get_cleared_selected_lines_faults_data_set(
#        clear_fault_selection_list;
#        deepcopy(on_fault_net_para)...)

# (;post_clear_fault_Ynet,
#  post_clear_fault_nodes_idx_with_adjacent_nodes_idx) =
#      NamedTupleTools.select(
#          cleared_selected_lines_faults_net_para ,
#          (:post_clear_fault_Ynet,
#           :post_clear_fault_nodes_idx_with_adjacent_nodes_idx))


# #----------------------------------------
# # write to tex
# #----------------------------------------

# rounded_post_clear_fault_Ynet =
#     round_up_Ynet(post_clear_fault_Ynet;
#                   fractional_digits=4 )

# object_for_tex =
#     (; rounded_post_clear_fault_Ynet,
#      post_clear_fault_nodes_idx_with_adjacent_nodes_idx)



# write_vector_or_matrix_to_tex(
#     object_for_tex,
#     tex_filename)

# #----------------------------------------

# # (:pre_clear_fault_Ynet,
# #  :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
# #  :pre_clear_fault_all_nodes_idx,
# #  :n2s_pre_clear_fault_all_nodes_idx,
# #  :pre_clear_fault_nodes_idx, :n2s_pre_clear_fault_nodes_idx,
# #  :pre_clear_list_Ya_nkf, :pre_clear_list_Ynkf_b,
# #  :pre_clear_list_faulty_line_Yl,
# #  :pre_clear_list_healthy_lines_Yl,
# #  :pre_clear_list_node_b_idx_in_a_node_row,
# #  :pre_clear_list_node_a_idx_in_b_node_row,
# #  :pre_clear_list_faulted_line_a_b_orientation,
# #  :pre_clear_list_fault_point_from_node_a,
# #  :pre_clear_list_fault_resistance,
# #  :pre_clear_list_no_line_circuit,
# #  :post_clear_fault_Ynet,
# #  :post_clear_fault_nodes_idx_with_adjacent_nodes_idx,
# #  :post_clear_fault_all_nodes_idx,
# #  :n2s_post_clear_fault_all_nodes_idx,
# #  :post_clear_fault_nodes_idx,
# #  :n2s_post_clear_fault_nodes_idx,
# #  :faulty_Ynet, :faulty_nodes_idx_with_adjacent_nodes_idx,
# #  :faulty_all_nodes_idx, :n2s_faulty_all_nodes_idx,
# #  :fault_nodes_idx, :n2s_fault_nodes_idx,
# #  :list_faulted_line_a_b_orientation,
# #  :list_fault_point_from_node_a,
# #  :list_fault_resistance,
# #  :list_no_line_circuit,
# #  :list_Ya_nkf, :list_Ynkf_b,
# #  :list_faulty_line_Yl,
# #  :list_healthy_lines_Yl,
# #  :list_node_b_idx_in_a_node_row,
# #  :list_node_a_idx_in_b_node_row)


# #-------------------------------
# # GraphRecipes
# #-------------------------------
# # SimpleDiGraph
# # SimpleGraph

# # https://docs.juliaplots.org/stable/generated/graph_attributes/#graph_attributes
# # https://docs.juliaplots.org/stable/GraphRecipes/examples/#graph_examples

# using GraphRecipes

 
# G_net_graph = SimpleDiGraph(
#     Edge.(edges_orientation))

# plt_G_net = graphplot(G_net_graph)

# dict_edges_label =
#     Dict( a_key => a_value
#           for (a_value, a_key) in
#               enumerate(edges_orientation))

# plt_G_net2 = plot(G_net_graph;
#                   names=collect(1:nv(G_net_graph)),
#                   edgelabel = dict_edges_label,
#                   edgelabel_offset=0.0) 

# #-------------------------------
# # StatsPlots
# #-------------------------------

# # https://docs.juliaplots.org/stable/generated/statsplots/

# using StatsPlots

# using DataFrames, IndexedTables


# #-------------------------------
# #-------------------------------


# (sta_red_vh_θh_0,
#  pf_PQ_param,
#  pf_kw_para) =
#      NamedTupleTools.select(
#          system_net_static_data,
#          (:sta_red_vh_θh_0,
#           :pf_PQ_param,
#           :pf_kw_para))

# #-------------------------------

# pf_fun_mismatch =
#     # get_a_model_integrated_pf_sta_ΔI_mismatch_generic
#     # get_generic_sta_pf_ΔI_mismatch
#     # get_a_model_integrated_pf_sta_ΔPQ_mismatch_generic
#     # get_generic_sta_pf_ΔPQ_mismatch
#     get_a_model_integrated_pf_sta_ΔPQ_mismatch


#  pf_sol =  NonlinearSolve.solve(
#     NonlinearProblem(
#         NonlinearFunction( ( g, x, p ) ->
#             pf_fun_mismatch(
#                 g, x, p;
#                 pf_kw_para = pf_kw_para )),
#         sta_red_vh_θh_0,
#         pf_PQ_param ),
#     NewtonRaphson() ;
#     abstol =
#         abstol,
#     reltol =
#         reltol )

# #-------------------------------
# # Results
# #-------------------------------

# generic_red_sol_kwd_para =
#     getproperty(system_net_static_data,
#                 :generic_red_sol_kwd_para)

# # getproperty(getproperty(generic_red_sol_kwd_para, :pf_kw_para), :pf_kw_net_para)

# (;pf_P_gens, pf_Q_gens,
#  pf_P_g_gens, pf_Q_g_gens,
#  vh, θh, θh_deg
#  # gens_nodes_idx,
#  # transformed_slack_gens_nodes_idx,
#  # transformed_gens_nodes_idx,
#  # transformed_non_gens_nodes_idx,
#  # transformed_all_nodes_idx
#   ) = NamedTupleTools.select(
#          get_results_static_pf_red_sol_u(
#              pf_sol;
#              generic_red_sol_kwd_para =
#                  generic_red_sol_kwd_para,
#              baseMVA = baseMVA,
#              basekV = 1.0) ,
#              (:pf_P_gens, :pf_Q_gens,             
#               :pf_P_g_gens, :pf_Q_g_gens,
             
#               :vh, :θh, :θh_deg
              
#              # :gens_nodes_idx,
#              # :transformed_slack_gens_nodes_idx,
#              # :transformed_gens_nodes_idx,
#              # :transformed_non_gens_nodes_idx,
#              # :transformed_all_nodes_idx
#               ) )

# #-------------------------------
# # Ynet, Yπ_net and Ybus
# #-------------------------------

# (;edge_data_from_json,
#  shunt_data_from_json,
#  baseMVA,
#  all_nodes_idx,
#  n2s_all_nodes_idx) =
#      NamedTupleTools.select(
#          system_net_static_data,
#          (:edge_data_from_json,
#           :shunt_data_from_json,
#           :baseMVA,
#           :all_nodes_idx,
#           :n2s_all_nodes_idx )  )

# #-------------------------------

# Ybus_and_related_matrices =
#     get_Ybus_and_related_matrices(
#         edge_data_from_json,
#         shunt_data_from_json;
#         baseMVA = baseMVA,
#         basekV = 1.0,
#         line_data_in_pu = true)


# Ynet_and_related_vectors =
#     get_Ynet_and_related_vectors(
#         edge_data_from_json,
#         # shunt_data_from_json;
#         # all_nodes_idx,
#         n2s_all_nodes_idx,
#         baseMVA = baseMVA,
#         baseShunt = baseMVA,
#         basekV = 1.0,
#         line_data_in_pu = true)


# Yπ_net_and_related_vectors =
#     get_Yπ_net_and_related_vectors(
#         edge_data_from_json,
#         shunt_data_from_json;
#         # all_nodes_idx,
#         # n2s_all_nodes_idx,
#         baseMVA = baseMVA,
#         baseShunt = baseMVA,
#         basekV = 1.0,
#         line_data_in_pu = true,
#         orientated_bool = false)

# #-------------------------------


# (;
#  Ybus,) =
#      NamedTupleTools.select(
#          Ybus_and_related_matrices,
#          (
#           :Ybus,))



# #-------------------------------

# (;
#  Ynet_wt_nodes_idx_wt_adjacent_nodes, ) =
#      NamedTupleTools.select(
#          Ynet_and_related_vectors,
#          (
#           :Ynet_wt_nodes_idx_wt_adjacent_nodes, ))


# (;Ynet,
#  nodes_idx_with_adjacent_nodes_idx) =
#     NamedTupleTools.select(
#         Ynet_wt_nodes_idx_wt_adjacent_nodes,
#         (:Ynet,
#          :nodes_idx_with_adjacent_nodes_idx))


# #-------------------------------



# (;
#  Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes, ) =
#      NamedTupleTools.select(
#          Yπ_net_and_related_vectors,
#          (
#           :Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes, ) )


# (;Yπ_net,
#  Yshunt,
#  nodes_idx_with_adjacent_nodes_idx ) =
#     NamedTupleTools.select(
#         Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
#         (:Yπ_net,
#          :Yshunt,
#          :nodes_idx_with_adjacent_nodes_idx ))




# """


# net_admittance_and_related_matrices_type =
#     (;
#      Ynet_and_related_vectors  =
#          get_Ynet_and_related_vectors(
#              edge_data_from_json,
#              shunt_data_from_json;
#              # all_nodes_idx,
#              # n2s_all_nodes_idx,
#              baseMVA = baseMVA,
#              baseShunt = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true ),

#      Yπ_net_and_related_vectors =
#          get_Yπ_net_and_related_vectors(
#              edge_data_from_json,
#              shunt_data_from_json;
#              # all_nodes_idx,
#              # n2s_all_nodes_idx,
#              baseMVA = baseMVA,
#              baseShunt = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true,
#              orientated_bool = false ),

#      Ybus_and_related_matrices =
#          get_Ybus_and_related_matrices(
#              edge_data_from_json,
#              shunt_data_from_json;
#              baseMVA = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true) )

#-------------------------------

# net_admittance_vector_type =
#     (;
#      Ynet_wt_nodes_idx_wt_adjacent_nodes  =
#          getproperty(get_Ynet_and_related_vectors(
#              edge_data_from_json,
#              shunt_data_from_json;
#              # all_nodes_idx,
#              # n2s_all_nodes_idx,
#              baseMVA = baseMVA,
#              baseShunt = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true ),
#                  :Ynet_wt_nodes_idx_wt_adjacent_nodes),

#      Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
#          getproperty(get_Yπ_net_and_related_vectors(
#              edge_data_from_json,
#              shunt_data_from_json;
#              # all_nodes_idx,
#              # n2s_all_nodes_idx,
#              baseMVA = baseMVA,
#              baseShunt = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true,
#              orientated_bool = false ),
#         :Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes),
#      Ybus =
#          getproperty(get_Ybus_and_related_matrices(
#              edge_data_from_json,
#              shunt_data_from_json;
#              baseMVA = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true ),
#                      :Ybus))

#-------------------------------

# net_admittance_vector_type =            
#    (; Ynet_wt_nodes_idx_wt_adjacent_nodes =
#         get_Ynet( edge_data_from_json,
#               shunt_data_from_json;
#               # all_nodes_idx,
#               # n2s_all_nodes_idx,
#               baseMVA = baseMVA,
#               baseShunt = baseMVA,
#               basekV = 1.0,
#               line_data_in_pu = true ),

#     Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
#         get_Yπ_net(edge_data_from_json,
#                  shunt_data_from_json;
#                  # all_nodes_idx,
#                  # n2s_all_nodes_idx,
#                  baseMVA = baseMVA,
#                  baseShunt = baseMVA,
#                  basekV = 1.0,
#                  line_data_in_pu = true,
#                  orientated_bool = false),

#     Ybus = getproperty(get_Ybus(edge_data_from_json,
#               shunt_data_from_json;
#               baseMVA = baseMVA,
#               basekV = 1.0,
#               line_data_in_pu = true ), :Ybus))


#-------------------------------
#-------------------------------


# net_admittance_and_related_matrices_type =
#     (;
#      Ynet_and_related_vectors  =
#          get_Ynet_and_related_vectors(
#              edge_data_from_json,
#              shunt_data_from_json;
#              # all_nodes_idx,
#              # n2s_all_nodes_idx,
#              baseMVA = baseMVA,
#              baseShunt = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true ),
     
#      Yπ_net_and_related_vectors =
#          get_Yπ_net_and_related_vectors(
#              edge_data_from_json,
#              # shunt_data_from_json;
#              # all_nodes_idx,
#              n2s_all_nodes_idx,
#              baseMVA = baseMVA,
#              baseShunt = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true,
#              orientated_bool = false ),
     
#      Ybus_and_related_matrices =
#          get_Ybus_and_related_matrices(
#              edge_data_from_json,
#              shunt_data_from_json;
#              baseMVA = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true) )

#-------------------------------

# net_admittance_vector_type =
#     (;
#      Ynet_wt_nodes_idx_wt_adjacent_nodes  =
#          getproperty(get_Ynet_and_related_vectors(
#              edge_data_from_json,
#              shunt_data_from_json;
#              baseMVA = baseMVA,
#              baseShunt = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true),
#                      :Ynet_wt_nodes_idx_wt_adjacent_nodes),
     
#      Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
#          getproperty(get_Yπ_net_and_related_vectors(
#              edge_data_from_json,
#              shunt_data_from_json;
#              baseMVA = baseMVA,
#              baseShunt = baseMVA,
#              basekV = 1.0,
#              line_data_in_pu = true,
#              orientated_bool = false),
#                      :Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes) ,
#      Ybus =
#     getproperty(get_Ybus_and_related_matrices(
#         edge_data_from_json,
#         shunt_data_from_json;
#         baseMVA = baseMVA,
#         basekV = 1.0,
#         line_data_in_pu = true), :Ybus))


#-------------------------------
#-------------------------------


# system_net_data_wt_static_parameters =
#     get_system_net_data_wt_static_parameters(
#         case_name;
#         script_dir,
#         data_dir =
#             data_dir,
#         json_net_data_by_components_file =
#             json_net_data_by_components_file,
#         components_libs_dir =
#             components_libs_dir )

# (;plant_generators_data_from_json,
#  plant_loads_data_from_json,
#  plant_transmission_data_from_json,
#  edge_data_from_json,
#  shunt_data_from_json,
 
#  dyn_pf_fun_kwd_net_idxs,
#  dyn_pf_fun_kwd_n2s_idxs,

#  all_nodes_idx,
#  n2s_all_nodes_idx,

#  pf_sta_ΔPQ_mismatch_parameters,
#  kwd_sta_sta_ΔPQ_sol_by_json,
#  generic_red_sol_kwd_para,
#  generic_dyn_sol_kwd_para,
#  states_Idx_syms_wt_funcs,
#  baseMVA) =
#      NamedTupleTools.select(
#          system_net_data_wt_static_parameters,
#          (:plant_generators_data_from_json,
#           :plant_loads_data_from_json,
#           :plant_transmission_data_from_json,
#           :edge_data_from_json,
#           :shunt_data_from_json,
          

#           :dyn_pf_fun_kwd_net_idxs,
#           :dyn_pf_fun_kwd_n2s_idxs,

#           :all_nodes_idx,
#           :n2s_all_nodes_idx,

#           :pf_sta_ΔPQ_mismatch_parameters,
#           :kwd_sta_sta_ΔPQ_sol_by_json,
#           :generic_red_sol_kwd_para,
#           :generic_dyn_sol_kwd_para,
#           :states_Idx_syms_wt_funcs,
#           :baseMVA))


# (pf_kw_para,
# red_types_Idxs_etc,
# pf_PQ_param) =
#     NamedTupleTools.select(
#         pf_sta_ΔPQ_mismatch_parameters,
#         (:pf_kw_para,
#          :red_types_Idxs_etc,
#          :pf_PQ_param) )

# #----------------------------------------

# pretty_summarysize(x) = Base.format_bytes(
#     Base.summarysize(x))

# getBytes(x::DataType) = sizeof(x)

# function getBytes(x)
#    total = 0;
#    fieldNames = fieldnames(typeof(x));
#    if fieldNames == []
#       return sizeof(x);
#    else
#      for fieldName in fieldNames
#         total += getBytes(getfield(x,fieldName));
#      end
#      return total;
#    end
# end

# #----------------------------------------
# #----------------------------------------


# Ybus_and_related_matrices =
#     get_Ybus_and_related_matrices(
#         edge_data_from_json,
#         shunt_data_from_json;
#         basekV = 1.0,
#         baseMVA = baseMVA,
#          line_data_in_pu = true)


# Ynet_and_related_vectors =
#     get_Ynet_and_related_vectors(
#         edge_data_from_json,
#         shunt_data_from_json;
#         # all_nodes_idx,
#         # n2s_all_nodes_idx,
#         baseMVA = baseMVA,
#         basekV = 1.0,
#         baseShunt = baseMVA,
#         line_data_in_pu = true)


# Yπ_net_and_related_vectors =
#     get_Yπ_net_and_related_vectors(
#         edge_data_from_json,
#         shunt_data_from_json;
#         # all_nodes_idx,
#         # n2s_all_nodes_idx,
#         baseMVA = baseMVA,
#         basekV = 1.0,
#         baseShunt = baseMVA,
#         line_data_in_pu = true,
#         orientated_bool = false)

# #---------------------------------------------------


# (;Yff,
#  Ytf,
#  Yft,
#  Ytt,
#  Y_sh,
#  Cf,
#  Ct,
#  Yf,
#  Yt,
#  Ybus) =
#      NamedTupleTools.select(
#          Ybus_and_related_matrices,
#          (:Yff,
#           :Ytf,
#           :Yft,
#           :Ytt,
#           :Y_sh,
#           :Cf,
#           :Ct,
#           :Yf,
#           :Yt,
#           :Ybus))

# #---------------------------------------------------


# (;nodes_idx_and_Yshunt,
#  edges_orientation,
#  nodes_incident_edges,
#  edges_Ybr_cal,
#  edges_orientation_and_edges_Ybr_cal,
#  nodes_incident_edges_and_orientation,
#  nodes_idx_with_adjacent_nodes_idx,
#  nodes_incident_edges_orientation_and_Ybr_cal,
#  Ynet_no_shunt,
#  Ynet_wt_nodes_idx_wt_adjacent_nodes ) =
#      NamedTupleTools.select(
#          Ynet_and_related_vectors,
#          (:nodes_idx_and_Yshunt ,
#           :edges_orientation,
#           :nodes_incident_edges,
#           :edges_Ybr_cal,
#           :edges_orientation_and_edges_Ybr_cal,
#           :nodes_incident_edges_and_orientation,
#           :nodes_idx_with_adjacent_nodes_idx,
#           :nodes_incident_edges_orientation_and_Ybr_cal,
#           :Ynet_no_shunt,
#           :Ynet_wt_nodes_idx_wt_adjacent_nodes ))


# (Ynet, nodes_idx_with_adjacent_nodes_idx) =
#     NamedTupleTools.select(
#         Ynet_wt_nodes_idx_wt_adjacent_nodes,
#         (:Ynet,
#          :nodes_idx_with_adjacent_nodes_idx))

# #---------------------------------------------------


# (;nodes_idx_and_Yshunt ,
#  edges_orientation,
#  nodes_incident_edges,
#  edges_Ybr_cal,
#  edges_orientation_and_edges_Ybr_cal,
#  nodes_incident_edges_and_orientation,
#  nodes_idx_with_adjacent_nodes_idx,
#  nodes_incident_edges_orientation_and_Ybr_cal,
#  Yπ_net,
#  Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes ) =
#      NamedTupleTools.select(
#          Yπ_net_and_related_vectors,
#          (:nodes_idx_and_Yshunt ,
#           :edges_orientation,
#           :nodes_incident_edges,
#           :edges_Ybr_cal,
#           :edges_orientation_and_edges_Ybr_cal,
#           :nodes_incident_edges_and_orientation,
#           :nodes_idx_with_adjacent_nodes_idx,
#           :nodes_incident_edges_orientation_and_Ybr_cal,
#           :Yπ_net,
#           :Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes ) )


# (;Yπ_net,
#  Yshunt,
#  nodes_idx_with_adjacent_nodes_idx ) =
#     NamedTupleTools.select(
#         Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
#         (:Yπ_net,
#          :Yshunt,
#          :nodes_idx_with_adjacent_nodes_idx ))

# pretty_summarysize(Ybus )


# pretty_summarysize(
#     Ynet_wt_nodes_idx_wt_adjacent_nodes )


# pretty_summarysize(Ynet )


# pretty_summarysize(
#     Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes )


# pretty_summarysize(Yπ_net )


# pretty_summarysize(
#     nodes_idx_with_adjacent_nodes_idx )


# #-----------------------------------------------    
# # Powerflow func and prob
# #-----------------------------------------------

# pf_sol =
#     get_pf_sta_ΔPQ_mismatch_sol_by_generic(
#         pf_PQ_param;
#         kwd_para =
#             kwd_sta_sta_ΔPQ_sol_by_json )

# #-------------------------------
# #-------------------------------


# (;pf_alg,
#  pf_kw_para,
#  red_vh_Idxs,
#  red_θh_Idxs,
#  sta_red_vh_θh_0) =
#      NamedTupleTools.select(
#          kwd_sta_sta_ΔPQ_sol_by_json,
#          (:pf_alg,
#           :pf_kw_para,
#           :red_vh_Idxs,
#           :red_θh_Idxs,
#           :sta_red_vh_θh_0))

# #-------------------------------
# #-------------------------------


# pf_fun_mismatch =
#     # get_a_model_integrated_pf_sta_ΔI_mismatch_generic
#     # get_generic_sta_pf_ΔI_mismatch
#     # get_a_model_integrated_pf_sta_ΔPQ_mismatch_generic
#     # get_generic_sta_pf_ΔPQ_mismatch
#     get_a_model_integrated_pf_sta_ΔPQ_mismatch


#  pf_sol =  NonlinearSolve.solve(
#     NonlinearProblem(
#         NonlinearFunction( ( g, x, p ) ->
#             pf_fun_mismatch(
#                 g, x, p;
#                 pf_kw_para = pf_kw_para )),
#         sta_red_vh_θh_0,
#         pf_PQ_param ),
#     pf_alg ;
#     abstol =
#         abstol,
#     reltol =
#         reltol )

# #-------------------------------
# #-------------------------------


# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_S_inj_Yπ_net

#     # get_ΔI_mismatch_by_Yπ_net
#     # get_ΔPQ_mismatch_by_Yπ_net
#     # get_ΔPQ_mismatch_by_ph_qh_Yπ_net

#  Yπ_net_pf_sol =  NonlinearSolve.solve(
#     NonlinearProblem(
#         NonlinearFunction( ( g, x, p ) ->
#             pf_fun_mismatch(
#                 g, x, p;
#                 pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#             Yπ_net_wt_Yshunt_wt_nodes_idx =
#         Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
#     use_Yπ_net_bool = true ) ),
#         sta_red_vh_θh_0,
#         pf_PQ_param ),
#     pf_alg ;
#     abstol =
#         abstol,
#     reltol =
#         reltol )

# #-------------------------------
# #-------------------------------

# Ybus_pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_Ybus
#     # get_ΔI_mismatch_by_Ybus


#  pf_sol =  NonlinearSolve.solve(
#     NonlinearProblem(
#         NonlinearFunction( ( g, x, p ) ->
#             pf_fun_mismatch(
#                 g, x, p;
#                 pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#             Ybus = Ybus,
#     use_Ybus_bool = true ) ),
#         sta_red_vh_θh_0,
#         pf_PQ_param ),
#     pf_alg ;
#     abstol =
#         abstol,
#     reltol =
#         reltol )

# #-------------------------------
# #-------------------------------

# pf_fun_mismatch =
#     get_ΔI_mismatch_by_Ynet
#     # get_ΔPQ_mismatch_by_Ynet


#  Ynet_pf_sol =  NonlinearSolve.solve(
#     NonlinearProblem(
#         NonlinearFunction( ( g, x, p ) ->
#             pf_fun_mismatch(
#                 g, x, p;
#                 pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para ) ) ),
#         sta_red_vh_θh_0,
#         pf_PQ_param ),
#     pf_alg ;
#     abstol =
#         abstol,
#     reltol =
#         reltol )

# #-------------------------------
# # NonlinearSolve autodiff
# #-------------------------------

# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_Ybus


# autodiff_nonlinearsolve_pf_sol = NonlinearSolve.solve(
#     NonlinearProblem(
#         NonlinearFunction(
#             ( g, x, p ) ->
#                 pf_fun_mismatch(
#                 g, x, p;
#                 pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#                 Ybus = Ybus) ) ,
#         sta_red_vh_θh_0,
#         pf_PQ_param),
#     NewtonRaphson();
#         abstol = 1e-10,
#         reltol = 1e-10)

# #-------------------------------
# # nlsolve autodiff
# #-------------------------------

# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_Ybus


# autodiff_nlsolve_pf_sol = nlsolve((g, x) ->
#     pf_fun_mismatch(
#         g, x, pf_PQ_param;
#         Ybus =
#             Ybus,
#         pf_kw_para =
#             disaggregate_sta_pf_keywords_parameter(
#                 pf_kw_para )),
#         sta_red_vh_θh_0,
#         method=:trust_region,
#         autodiff=:forward )
              
# #-------------------------------
# # NonlinearSolve with Jac
# #-------------------------------

# red_ΔPQ_x = similar(sta_red_vh_θh_0)


# Jac_row_size =
#     Jac_col_size = length( sta_red_vh_θh_0 )


# Jac_vh_θh =
#     spzeros( Jac_row_size, Jac_col_size )


# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_Ybus


# jac_nonlinearsolve_pf_sol = NonlinearSolve.solve(
#     NonlinearProblem(
#         NonlinearFunction(
#             ( g, x, p ) ->
#                 pf_fun_mismatch(
#                 g, x, p;
#                 pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#                     Ybus = Ybus,
#                     use_autodiff =
#                         true),
#         jac = (Jac_vh_θh, x, p) -> sta_pf_Jac!(
#             Jac_vh_θh, x, p ;
#             Ybus =
#                 Ybus,
#             pf_kw_para =
#                 disaggregate_sta_pf_keywords_parameter(
#                     pf_kw_para) ) ) ,
#         sta_red_vh_θh_0,
#         pf_PQ_param),
#     NewtonRaphson();
#         abstol = 1e-10,
#     reltol = 1e-10)

# #-------------------------------
# # nlsolve with Jac
# #-------------------------------

# red_ΔPQ_x = similar(sta_red_vh_θh_0)


# Jac_row_size =
#     Jac_col_size = length( sta_red_vh_θh_0 )


# Jac_vh_θh =
#     spzeros( Jac_row_size, Jac_col_size )


# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_Ybus


# jac_nlsolve_pf_sol = nlsolve(
#     (red_ΔPQ_x, red) ->
#         pf_fun_mismatch(
#             red_ΔPQ_x,
#             red,
#             pf_PQ_param;
#             Ybus =
#                 Ybus,
#             pf_kw_para =
#                 disaggregate_sta_pf_keywords_parameter(
#                     pf_kw_para),
#             use_autodiff =
#                 false),
#     (Jac_vh_θh, red )->
#         sta_pf_Jac!(
#             Jac_vh_θh,
#             red,
#             pf_PQ_param;            
#             Ybus=
#                 Ybus,
#             pf_kw_para =
#                 disaggregate_sta_pf_keywords_parameter(
#                     pf_kw_para) ) ,
#     sta_red_vh_θh_0;
#     ftol = 1e-10 )

# #-------------------------------
# #-------------------------------

# nlsolve_sol_properties =
#     (:method, :initial_x, :zero, :residual_norm,
#      :iterations, :x_converged, :xtol, :f_converged,
#      :ftol, :trace, :f_calls, :g_calls )

# nonlinearsolve_sol_properties =
#     (:u, :resid, :prob, :alg, :retcode,
#      :original, :left, :right, :stats, :trace )

# jac_nlsolve_pf_sol.method

# autodiff_nlsolve_pf_sol.method

# jac_nonlinearsolve_pf_sol.stats

# autodiff_nonlinearsolve_pf_sol.stats

# #-------------------------------
# # Using Jacobian-Free Newton-Krylov
# #-------------------------------

# using SparseConnectivityTracer

# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_Ybus

#  KrylovJL_GMRES_pf_sol =  NonlinearSolve.solve(
#     NonlinearProblem(
#         NonlinearFunction( ( g, x, p ) ->
#             pf_fun_mismatch(
#                 g, x, p;
#                 pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#                 Ybus =
#                     Ybus) ;sparsity = TracerSparsityDetector()  ),
#         sta_red_vh_θh_0,
#         pf_PQ_param ),
#     NewtonRaphson(linsolve = KrylovJL_GMRES() ) ;
#     abstol =
#         abstol,
#     reltol =
#         reltol )

# #-------------------------------
# #-------------------------------
# # experimental 
# #-------------------------------
# #-------------------------------

# using SparseConnectivityTracer, ADTypes

# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_Ybus


# f! = (du, u) -> pf_fun_mismatch(
#                     du, u, pf_PQ_param;
#                     pf_kw_para =
#                         disaggregate_sta_pf_keywords_parameter(
#                             pf_kw_para),
#                     Ybus =
#                         Ybus)

# du0 = similar(sta_red_vh_θh_0)

# jac_sparsity = ADTypes.jacobian_sparsity(
#     f!, du0, sta_red_vh_θh_0, TracerSparsityDetector())

# KLUFactorization_pf_sol = NonlinearSolve.solve(
#     NonlinearProblem(
#         NonlinearFunction(
#             ( g, x, p ) ->
#                 pf_fun_mismatch(
#                     g, x, p;
#                     pf_kw_para =
#                         disaggregate_sta_pf_keywords_parameter(
#                             pf_kw_para),
#                     Ybus =
#                         Ybus) ;sparsity = TracerSparsityDetector() ) ,
#         sta_red_vh_θh_0,
#         pf_PQ_param),
#     NewtonRaphson(linsolve = KLUFactorization());
#         abstol = 1e-10,
#         reltol = 1e-10 )

# #-------------------------------

# # https://docs.sciml.ai/NonlinearSolve/stable/tutorials/large_systems/

# using SparseConnectivityTracer , ADTypes


# f! = (du, u) -> get_ΔPQ_mismatch_by_Ybus(
#     du, u, (pf_PQ_param, Ybus);
#     pf_kw_para =
#         disaggregate_sta_pf_keywords_parameter(
#             pf_kw_para),
#     use_Ybus_bool =
#         true )

# du0 = similar(sta_red_vh_θh_0)

# jac_sparsity = ADTypes.jacobian_sparsity(
#     f!, du0, sta_red_vh_θh_0,
#     TracerSparsityDetector() )

# #---------------------------------------------------

# (sp_ybus_I, sp_ybus_J, sp_ybus_nzv) = findnz(Ybus)

# sp_ybus_nzv_real = real.(sp_ybus_nzv)

# sp_ybus_nzv_imag = imag.(sp_ybus_nzv)

# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_decomp_Ybus

# pf_prob_Ybus_sparse =
#     NonlinearProblem(
#         NonlinearFunction(
#     (g, x, p ) ->
#             pf_fun_mismatch(
#                 g, x, p;

#                 pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#                 sp_ybus_I =
#                     sp_ybus_I,
#                 sp_ybus_J=
#                     sp_ybus_J);
#     jac_prototype = SparseMatrixCSC) ,
#         sta_red_vh_θh_0,
#         (pf_PQ_param,
#          sp_ybus_nzv_real,
#          sp_ybus_nzv_imag) ;
#         abstol = 1e-10,
#         reltol = 1e-10)

# #-------------------------------

# using Symbolics


# red_ΔPQ_x = similar(sta_red_vh_θh_0)

# red_vh_θh_DiffCache =
#     DiffCache(sta_red_vh_θh_0)


# ( sp_ybus_I , sp_ybus_J , sp_ybus_nzv) = findnz(Ybus)

# sp_ybus_nzv_real = real.(sp_ybus_nzv)

# sp_ybus_nzv_imag = imag.(sp_ybus_nzv)

# I_bus = zeros(size(Ybus)[1] )

# I_bus_cache = DiffCache(I_bus, length(I_bus) )

# sparse_nzvalues_DiffCache =
#     DiffCache(sp_ybus_nzv)

# sp_ybus_nzv_real_cache =
#     DiffCache(sp_ybus_nzv_real)

# sp_ybus_nzv_imag_cache =
#     DiffCache(sp_ybus_nzv_imag)


# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_decomp_Ybus


# sd = SymbolicsSparsityDetection()

# adtype = AutoSparse(AutoFiniteDiff())

# f = (y,x) -> pf_fun_mismatch(
#                 y, x, pf_PQ_param;
#                 pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para,
#                     ),
#                     sp_ybus_nzv_real =
#                         sp_ybus_nzv_real,
#                     sp_ybus_nzv_imag =
#                         sp_ybus_nzv_imag,
#                     sp_ybus_I =
#                         sp_ybus_I,
#                     sp_ybus_J =
#                         sp_ybus_J,
#                     sp_ybus_nzv_real_cache =
#                         sp_ybus_nzv_real_cache,
#                     sp_ybus_nzv_imag_cache =
#                         sp_ybus_nzv_imag_cache,
#                     sparse_nzvalues_DiffCache =
#                         sparse_nzvalues_DiffCache,
#                     I_bus_cache =
#                         I_bus_cache,
#                     I_bus = I_bus)

# cache = sparse_jacobian_cache(adtype, sd, f, red_ΔPQ_x, sta_red_vh_θh_0 )

# J = sparse_jacobian(adtype, cache, f, red_ΔPQ_x, sta_red_vh_θh_0)


# #-----------------------------------------------


# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_decomp_Ybus

    
# full_df_dx = ForwardDiff.jacobian(
#     ( x ) ->
#         pf_fun_mismatch(
#             red_ΔPQ_x, x, pf_PQ_param;
#         pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#                     Ybus =
#                         Ybus,
#                     sp_ybus_nzv_real =
#                         sp_ybus_nzv_real,
#                     sp_ybus_nzv_imag =
#                         sp_ybus_nzv_imag,
#                     sp_ybus_I =
#                         sp_ybus_I,
#                     sp_ybus_J =
#                         sp_ybus_J,
#                     red_vh_θh_DiffCache =
#                         red_vh_θh_DiffCache ,
#                     sp_ybus_nzv_real_cache =
#                         sp_ybus_nzv_real_cache,
#                     sp_ybus_nzv_imag_cache =
#                         sp_ybus_nzv_imag_cache,
#                     sparse_nzvalues_DiffCache =
#                         sparse_nzvalues_DiffCache,
#                     I_bus_cache =
#                         I_bus_cache,
#                     I_bus = I_bus),
#     sta_red_vh_θh_0 )

    
# full_df_dx = ForwardDiff.jacobian(
#     ( dx, x ) ->
#         pf_fun_mismatch(
#             dx, x, pf_PQ_param;
#         pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#                     Ybus =
#                         Ybus,
#                     sp_ybus_nzv_real =
#                         sp_ybus_nzv_real,
#                     sp_ybus_nzv_imag =
#                         sp_ybus_nzv_imag,
#                     sp_ybus_I =
#                         sp_ybus_I,
#                     sp_ybus_J =
#                         sp_ybus_J,
#                     red_vh_θh_DiffCache =
#                         red_vh_θh_DiffCache ,
#                     sp_ybus_nzv_real_cache =
#                         sp_ybus_nzv_real_cache,
#                     sp_ybus_nzv_imag_cache =
#                         sp_ybus_nzv_imag_cache,
#                     sparse_nzvalues_DiffCache =
#                         sparse_nzvalues_DiffCache,
#                     I_bus_cache =
#                         I_bus_cache,
#                     I_bus = I_bus),
#     red_ΔPQ_x, sta_red_vh_θh_0 )

    
# full_df_dp = ForwardDiff.jacobian(
#     ( dx, p ) ->
#         pf_fun_mismatch(
#             dx, sta_red_vh_θh_0, p;
#         pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#                     Ybus =
#                         Ybus,
#                     sp_ybus_nzv_real =
#                         sp_ybus_nzv_real,
#                     sp_ybus_nzv_imag =
#                         sp_ybus_nzv_imag,
#                     sp_ybus_I =
#                         sp_ybus_I,
#                     sp_ybus_J =
#                         sp_ybus_J,
#                     red_vh_θh_DiffCache =
#                         red_vh_θh_DiffCache ,
#                     sp_ybus_nzv_real_cache =
#                         sp_ybus_nzv_real_cache,
#                     sp_ybus_nzv_imag_cache =
#                         sp_ybus_nzv_imag_cache,
#                     sparse_nzvalues_DiffCache =
#                         sparse_nzvalues_DiffCache,
#                     I_bus_cache =
#                         I_bus_cache,
#                     I_bus = I_bus),
#     red_ΔPQ_x, pf_PQ_param )




# ff = NonlinearFunction(
#             ( g, x, p ) ->
#                 pf_fun_mismatch(
#                 g, x, p;
#                 pf_kw_para =
#                     disaggregate_sta_pf_keywords_parameter(
#                         pf_kw_para),
#                     Ybus =
#                         Ybus,
#                     sp_ybus_nzv_real =
#                         sp_ybus_nzv_real,
#                     sp_ybus_nzv_imag =
#                         sp_ybus_nzv_imag,
#                     sp_ybus_I =
#                         sp_ybus_I,
#                     sp_ybus_J =
#                         sp_ybus_J,
#                     red_vh_θh_DiffCache =
#                         red_vh_θh_DiffCache ,
#                     sp_ybus_nzv_real_cache =
#                         sp_ybus_nzv_real_cache,
#                     sp_ybus_nzv_imag_cache =
#                         sp_ybus_nzv_imag_cache,
#                     sparse_nzvalues_DiffCache =
#                         sparse_nzvalues_DiffCache,
#                     I_bus_cache =
#                         I_bus_cache,
#                     I_bus = I_bus) )


# pf_sol = NonlinearSolve.solve(
#     NonlinearProblem(ff,        
#         sta_red_vh_θh_0,
#         pf_PQ_param ),
#     NewtonRaphson( );
#         abstol = 1e-10,
#     reltol = 1e-10)


# pf_sol = NonlinearSolve.solve(
#     NonlinearProblem(ff,        
#         sta_red_vh_θh_0,
#         pf_PQ_param ),
#     NewtonRaphson(linsolve = KrylovJL_GMRES() );
#         abstol = 1e-10,
#     reltol = 1e-10)


# pf_sol = NonlinearSolve.solve(
#     NonlinearProblem(ff,        
#         sta_red_vh_θh_0,
#         pf_PQ_param ),
#     NewtonRaphson(linsolve = KLUFactorization());
#         abstol = 1e-10,
#     reltol = 1e-10)

# #-------------------------------

# pf_fun_mismatch =
#     get_ΔPQ_mismatch_by_sparse_Ybus



# #----------------------------------------
# #----------------------------------------
# # Result
# #----------------------------------------
# #----------------------------------------    

# get_results_static_pf_red_sol_u(
#     pf_sol;
#     generic_red_sol_kwd_para =
#         generic_red_sol_kwd_para,
#     baseMVA = 1.0,
#     basekV = 1.0)


# generic_results_pf_sta_red_sol =
#     get_generic_results_pf_sta_red_sol_u(
#         pf_sol;
#         generic_red_sol_kwd_para =
#             generic_red_sol_kwd_para,
#         baseMVA =
#             baseMVA,
#         basekV =
#             1.0 )
    
# (pf_P_gens,
#  pf_Q_gens,
#  vh,
#  θh,
#  gens_vh,
#  gens_θh,
#  θh_deg) =
#     NamedTupleTools.select(
#         generic_results_pf_sta_red_sol,
#         (:pf_P_gens,
#          :pf_Q_gens,
#          :vh,
#          :θh,
#          :gens_vh,
#          :gens_θh,
#          :θh_deg) )


# #----------------------------------------    
# #----------------------------------------    

# (;
#  slack_gens_nodes_idx,
#  non_slack_gens_nodes_idx,
#  gens_nodes_idx,
#  non_gens_nodes_idx,
#  gens_with_loc_load_idx,
#  all_nodes_idx ) =
#      NamedTupleTools.select(
#          dyn_pf_fun_kwd_net_idxs,
#          (
#           :slack_gens_nodes_idx,
#           :non_slack_gens_nodes_idx,
#           :gens_nodes_idx,
#           :non_gens_nodes_idx,
#           :gens_with_loc_load_idx,
#           :all_nodes_idx ))

# (;n2s_slack_gens_idx,
#  n2s_non_slack_gens_idx,
#  n2s_gens_idx,
#  n2s_non_gens_idx,
#  n2s_gens_with_loc_load_idxs,
#  n2s_all_nodes_idx ) =
#      NamedTupleTools.select(
#          dyn_pf_fun_kwd_n2s_idxs,
#          (:n2s_slack_gens_idx,
#           :n2s_non_slack_gens_idx,
#           :n2s_gens_idx,
#           :n2s_non_gens_idx,
#           :n2s_gens_with_loc_load_idxs,
#           :n2s_all_nodes_idx ) )

# loc_load_exist =
#     getproperty(
#         pf_kw_para,
#         :loc_load_exist)

# (;P_gens_sta_para_Idxs,
#  Q_gens_sta_para_Idxs,
#  P_non_gens_sta_para_Idxs,
#  Q_non_gens_sta_para_Idxs,
#  P_g_loc_load_sta_para_Idxs,
#  Q_g_loc_load_sta_para_Idxs ) =
#      NamedTupleTools.select(
#          getproperty(
#              pf_kw_para,
#              :pf_kw_PQ_para_idxs),
#          (:P_gens_sta_para_Idxs,
#           :Q_gens_sta_para_Idxs,
#           :P_non_gens_sta_para_Idxs,
#           :Q_non_gens_sta_para_Idxs,
#           :P_g_loc_load_sta_para_Idxs,
#           :Q_g_loc_load_sta_para_Idxs ) )



# P_gens =
#     pf_PQ_param[
#         P_gens_sta_para_Idxs ]

# Q_gens =
#     pf_PQ_param[
#         Q_gens_sta_para_Idxs ]

# P_non_gens  =
#     pf_PQ_param[
#         P_non_gens_sta_para_Idxs ]

# Q_non_gens = 
#     pf_PQ_param[
#         Q_non_gens_sta_para_Idxs ]

# if loc_load_exist == true

#     P_g_loc_load =
#         pf_PQ_param[
#             P_g_loc_load_sta_para_Idxs ]

#     Q_g_loc_load =
#         pf_PQ_param[
#             Q_g_loc_load_sta_para_Idxs ]

# else

#     P_g_loc_load = [0.0]

#     Q_g_loc_load = [0.0]
# end


# lines_ph_pk_qh_qk =
#     get_ph_pk_qh_qk_by_Ybr_cal(
#         vh,
#         θh,
#         edges_Ybr_cal;
#         edges_orientation,
#         n2s_all_nodes_idx)


# # map(()-> (), edges_orientation, Yπ_net)


# ∑_S_injection_by_Yπ_net =
#     round.(get_nodes_∑_Sh_injection_by_Yπ_net(
#         vh,
#         θh;
#         Yπ_net,
#         Yshunt,
#         nodes_idx_with_adjacent_nodes_idx,
#         n2s_all_nodes_idx ); digits=4)


# round.(real.(∑_S_injection_by_Yπ_net );digits=4)

# round.(imag.(∑_S_injection_by_Yπ_net );digits=4)


# P_nodes_injections_lines_by_ynet =
#     [ nth_idx ∈ non_gens_nodes_idx ?
#      (
#      vh[ n2s_all_nodes_idx[nth_idx]] *
#      sum([ vh[ n2s_all_nodes_idx[ idx]] *
#      abs(ynj) *
#      cos( θh[ n2s_all_nodes_idx[nth_idx]] -
#      θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
#                for (ynj, idx) in
#                    zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                         nodes_idx_with_adjacent_nodes_idx[
#                             n2s_all_nodes_idx[nth_idx]])])) :
#                             nth_idx ∈ gens_with_loc_load_idx ?
#      (
#       vh[ n2s_all_nodes_idx[nth_idx]] *
#       sum([ vh[ n2s_all_nodes_idx[idx]] *
#       abs(ynj) *
#       cos(θh[ n2s_all_nodes_idx[ nth_idx]] -
#       θh[ n2s_all_nodes_idx[ idx]] -
#       angle(ynj) )
#             for (ynj, idx) in
#                 zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
#                      nodes_idx_with_adjacent_nodes_idx[
#                          n2s_all_nodes_idx[ nth_idx]])])) :
#     (
#     vh[ n2s_all_nodes_idx[nth_idx]] *
#     sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#     cos(θh[ n2s_all_nodes_idx[nth_idx]] -
#     θh[ n2s_all_nodes_idx[idx]] -
#     angle(ynj) )
#           for (ynj, idx) in
#               zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
#                    nodes_idx_with_adjacent_nodes_idx[
#                        n2s_all_nodes_idx[ nth_idx] ]) ]))
#                     for nth_idx in all_nodes_idx ] 

# red_P_mismatch =
#     P_mismatch[
#         setdiff(transformed_all_nodes_idx,
#                 transformed_slack_gens_nodes_idx)]


# P_nodes_injections_lines_by_yπ_net =
#     round.([ nth_idx ∈ non_gens_nodes_idx ?
#      (
#      ( vh[ n2s_all_nodes_idx[nth_idx]] )^2 * real(
#          Yshunt[ n2s_all_nodes_idx[nth_idx] ] ) +
#      sum([ get_a_node_ph_injection_by_Yπ_net(
#                vh[ n2s_all_nodes_idx[nth_idx]],
#                θh[ n2s_all_nodes_idx[nth_idx]],
#                vh[ n2s_all_nodes_idx[ idx]],
#                θh[ n2s_all_nodes_idx[ idx]]; y_π)

#                for (y_π, idx) in
#                    zip( Yπ_net[ n2s_all_nodes_idx[nth_idx] ],
#                         nodes_idx_with_adjacent_nodes_idx[
#                             n2s_all_nodes_idx[nth_idx]]) ]) ) :
#                             nth_idx ∈ gens_with_loc_load_idx ?
#      (
#       (vh[ n2s_all_nodes_idx[nth_idx]] )^2 * real(
#          Yshunt[ n2s_all_nodes_idx[nth_idx] ] ) -
#       sum([ get_a_node_ph_injection_by_Yπ_net(
#                vh[ n2s_all_nodes_idx[nth_idx]],
#                θh[ n2s_all_nodes_idx[nth_idx]],
#                vh[ n2s_all_nodes_idx[ idx]],
#                θh[ n2s_all_nodes_idx[ idx]]; y_π)
#             for (y_π, idx) in
#                 zip( Yπ_net[ n2s_all_nodes_idx[nth_idx]],
#                      nodes_idx_with_adjacent_nodes_idx[
#                          n2s_all_nodes_idx[ nth_idx]])])) :
#     (
#     (vh[ n2s_all_nodes_idx[nth_idx]] )^2 * real(
#          Yshunt[ n2s_all_nodes_idx[nth_idx] ] ) -
#       sum([ get_a_node_ph_injection_by_Yπ_net(
#                vh[ n2s_all_nodes_idx[nth_idx]],
#                θh[ n2s_all_nodes_idx[nth_idx]],
#                vh[ n2s_all_nodes_idx[ idx]],
#                θh[ n2s_all_nodes_idx[ idx]]; y_π)
#           for (y_π, idx) in
#               zip( Yπ_net[ n2s_all_nodes_idx[ nth_idx] ],
#                    nodes_idx_with_adjacent_nodes_idx[
#                        n2s_all_nodes_idx[ nth_idx] ]) ]))
#                     for nth_idx in all_nodes_idx ];digits=4) 


# P_nodal_components = [
#     nth_idx ∈ non_gens_nodes_idx ?
#         P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] :
#         nth_idx ∈ gens_with_loc_load_idx ?
#         P_gens[ n2s_gens_idx[nth_idx]]  -
#         P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]]  : 
#         P_gens[ n2s_gens_idx[nth_idx]]        
#     for nth_idx in all_nodes_idx ]


# Q_nodal_components = [
#     nth_idx ∈ non_gens_nodes_idx ?
#         Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] :
#         nth_idx ∈ gens_with_loc_load_idx ?
#         Q_gens[ n2s_gens_idx[ nth_idx ]] -
#         Q_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] :
#         Q_gens[ n2s_gens_idx[nth_idx]]    
#     for nth_idx in all_nodes_idx ]



# P_nodes_injections_lines = round.([
#     nth_idx ∈ non_gens_nodes_idx ?
#         (
#         real(∑_S_injection_by_Yπ_net[  n2s_non_gens_idx[ nth_idx ] ] ) ) :
#         nth_idx ∈ gens_with_loc_load_idx ?
#         (
#             real(∑_S_injection_by_Yπ_net[
#                 n2s_gens_with_loc_load_idxs[ nth_idx ]] )) : 
#         (
#         real(∑_S_injection_by_Yπ_net[ n2s_gens_idx[ nth_idx ] ])  )        
#     for nth_idx in all_nodes_idx ];digits=4)




# Q_nodes_injections_lines = [ nth_idx ∈ non_gens_nodes_idx ?
#     (
#     (vh[ n2s_all_nodes_idx[nth_idx]] )^2 * imag(
#          Yshunt[ n2s_all_nodes_idx[nth_idx] ] ) +
#     sum([ get_an_edge_qh_injection(
#                vh[ n2s_all_nodes_idx[nth_idx]],
#                θh[ n2s_all_nodes_idx[nth_idx]],
#                vh[ n2s_all_nodes_idx[ idx]],
#                θh[ n2s_all_nodes_idx[ idx]]; y_π)
#                for (y_π, idx) in
#                    zip( Yπ_net[ n2s_all_nodes_idx[nth_idx] ],
#                         nodes_idx_with_adjacent_nodes_idx[
#                             n2s_all_nodes_idx[nth_idx]])])) :
#                             nth_idx ∈ gens_with_loc_load_idx ?
#     ( (vh[ n2s_all_nodes_idx[nth_idx]] )^2 * imag(
#          Yshunt[ n2s_all_nodes_idx[nth_idx] ] ) -
#     sum([ get_an_edge_qh_injection(
#                vh[ n2s_all_nodes_idx[nth_idx]],
#                θh[ n2s_all_nodes_idx[nth_idx]],
#                vh[ n2s_all_nodes_idx[ idx]],
#                θh[ n2s_all_nodes_idx[ idx]]; y_π)
#                for (y_π, idx) in
#                    zip(Yπ_net[ n2s_all_nodes_idx[nth_idx]],
#                        nodes_idx_with_adjacent_nodes_idx[
#                            n2s_all_nodes_idx[nth_idx]])])) :
#    (
#     (vh[ n2s_all_nodes_idx[nth_idx]] )^2 * imag(
#          Yshunt[ n2s_all_nodes_idx[nth_idx] ] ) +
#     sum([ get_an_edge_qh_injection(
#                vh[ n2s_all_nodes_idx[nth_idx]],
#                θh[ n2s_all_nodes_idx[nth_idx]],
#                vh[ n2s_all_nodes_idx[ idx]],
#                θh[ n2s_all_nodes_idx[ idx]]; y_π)
#               for (y_π, idx) in
#                   zip(Yπ_net[ n2s_all_nodes_idx[nth_idx]],
#                       nodes_idx_with_adjacent_nodes_idx[
#                           n2s_all_nodes_idx[nth_idx]])] ) ) 
#                for nth_idx in all_nodes_idx ]



# Q_nodes_without_lines = round.([ nth_idx ∈ non_gens_nodes_idx ?
#     Q_non_gens[ n2s_non_gens_idx[ nth_idx ]]  :
#     nth_idx ∈ gens_with_loc_load_idx ?
#     (pf_Q_gens[ n2s_gens_idx[ nth_idx ]] - P_g_loc_load[
#         n2s_gens_with_loc_load_idxs[nth_idx]]) :
#     pf_Q_gens[ n2s_gens_idx[nth_idx]]    
#                for nth_idx in all_nodes_idx ];digits=4)


# Q_nodes_without_lines = round.([
#     Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
#         imag(∑_S_injection_by_Yπ_net[
#             n2s_non_gens_idx[ nth_idx ]] ) 
#                for nth_idx in all_nodes_idx
#                    if nth_idx ∈ non_gens_nodes_idx  ];digits=4)

# Q_mismatch = [
#     Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
#         imag(∑_S_injection_by_Yπ_net[
#             n2s_non_gens_idx[ nth_idx ]] ) 
#                for nth_idx in all_nodes_idx
#                    if nth_idx ∈ non_gens_nodes_idx  ]

# # pf_P_gens
# # pf_Q_gens

# #----------------------------------------
# #----------------------------------------


# #---------------------------------------------------

# uh = vh .* exp.(im * θh)

# Ia = get_nodes_∑_ynj_x_vj(
#     vh,
#     θh,
#     Ynet,
#     nodes_idx_with_adjacent_nodes_idx,
#     n2s_all_nodes_idx )

# Ib = get_nodes_∑_ynj_x_vj(
#     uh,
#     Ynet,
#     nodes_idx_with_adjacent_nodes_idx )

# Ic = get_nodes_∑_ynj_x_vj_by_Yπ_net(
#     vh,
#     θh,
#     Yπ_net,
#     Yshunt,
#     nodes_idx_with_adjacent_nodes_idx,
#     n2s_all_nodes_idx )

# Id = get_nodes_∑_ynj_x_vj_by_Yπ_net(
#     uh,
#     Yπ_net,
#     Yshunt,
#     nodes_idx_with_adjacent_nodes_idx )

# #---------------------------------------------------

# get_nodes_∑_S_injection_by_Yπ_net(
#     vh,
#     θh;
#     Yπ_net,
#     Yshunt,
#     nodes_idx_with_adjacent_nodes_idx,
#     n2s_all_nodes_idx )

# #---------------------------------------------------
# #---------------------------------------------------
# # Sensitivity of net loss function to
# # vh, θh  and transmission line parameters
# #---------------------------------------------------
# #---------------------------------------------------

# status_steady_state_parameters =
#     get_status_steady_state_parameters(
#         net_data_by_components_file;
#         components_libs_dir =
#             components_libs_dir ,    
#         pf_alg = NewtonRaphson(),

#         abstol =
#             abstol,    
#         reltol =
#             reltol,
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
#         basekV = 1.0,    
#         use_pu_in_PQ = true,
#         line_data_in_pu = true,

#         use_init_u0 = false,    
#         use_nlsolve = false,

#         with_faults = false )



# (Ynet,
#  nodes_idx_with_adjacent_nodes_idx) =
#     NamedTupleTools.select(
#         getproperty(status_steady_state_parameters,
#                     :Ynet_wt_nodes_idx_wt_adjacent_nodes),
#         (:Ynet,
#          :nodes_idx_with_adjacent_nodes_idx) )


# (faulty_Ynet,
#  faulty_nodes_idx_with_adjacent_nodes_idx) =
#     NamedTupleTools.select(
#         getproperty(status_steady_state_parameters,
#                     :on_fault_net_para),
#         (:faulty_Ynet,
#          :faulty_nodes_idx_with_adjacent_nodes_idx) )


# (pre_clear_fault_Ynet,
#  pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
#  post_clear_fault_Ynet,
#  post_clear_fault_nodes_idx_with_adjacent_nodes_idx) =
#      NamedTupleTools.select(
#          getproperty(status_steady_state_parameters,
#                      :cleared_selected_lines_faults_net_para),
#          (:pre_clear_fault_Ynet,
#           :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
#           :post_clear_fault_Ynet,
#           :post_clear_fault_nodes_idx_with_adjacent_nodes_idx))


# r_Ynet =
#     map((x) -> round.(x; digits=4), Ynet)

# r_faulty_Ynet =
#     map((x) -> round.(x; digits=4), faulty_Ynet)

# r_pre_clear_fault_Ynet =
#     map((x) -> round.(x; digits=4), pre_clear_fault_Ynet)

# r_post_clear_fault_Ynet =
#     map((x) -> round.(x; digits=4), post_clear_fault_Ynet)




# list_network_status = 
#         [:pre_fault_state,
#          :post_fault_state]






# ntuple_status_steady_state_data =
#     get_ntuple_status_steady_state_data(
#         ;with_faults =
#             with_faults,
#         net_data_by_components_file =
#             net_data_by_components_file,
#         components_libs_dir =
#             components_libs_dir,

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
#             basekV,    
#         use_pu_in_PQ =
#             use_pu_in_PQ,
#         line_data_in_pu =
#             line_data_in_pu,
#         list_network_status =
#             list_network_status )


# a_status_steady_state_data =
#     get_a_status_steady_state_data(
#         :pre_fault_state;
#         with_faults = false,
#         net_data_by_components_file =
#             net_data_by_components_file,
#         components_libs_dir =
#             components_libs_dir,        
        
#         timespan   = timespan,

#         on_fault_time = on_fault_time,
#         clear_fault_time = clear_fault_time,

#             list_fault_point_from_node_a =
#                 list_fault_point_from_node_a,
#             list_fault_resistance =
#                 list_fault_resistance,
#             list_no_line_circuit =
#                 list_no_line_circuit,

#             list_edges_to_have_fault =
#                 list_edges_to_have_fault,
#             clear_fault_selection_list =
#                 clear_fault_selection_list,

#         basekV = 1.0,    
#         use_pu_in_PQ = true,
#     line_data_in_pu = true)



# (baseMVA,
#  basekV,

#  vh,
#  θh,
#  u0_model_states_init,
#  model_syms,
#  model_mass_matrix,

#  ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
#  plants_cb_paras_switches,

#  generic_system_dynamics_kwd_para,

#  gens_nodes_names,
#  SM_gens_nodes_names,
#  non_gens_nodes_names,

#  cb_states,
 
#  dyn_pf_fun_kwd_n2s_idxs,
#  dyn_pf_fun_kwd_net_idxs,
 
#  Pg_Qg_Png_Qng_Pll_Qll_Idx,
#  Png_Qng_Pll_Qll_Idx,
#  Pg_Png_Qng_Idx,

#  pf_vh_θh_idx_and_idx2Idx,
 
#  dyn_pf_flat_vh_flat_θh_Idx,

#  edges_r,
#  edges_x,
#  edges_b,
#  edges_ratio,
#  edges_angle,
#  edges_type,
#  Gs,
#  Bs,
 
#  Ynet_wt_nodes_idx_wt_adjacent_nodes,
#  Ybr_cal_and_edges_orientation,

#  Pg_Qg_Png_Qng_Pll_Qll,
#  loc_load_exist,
#  slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

#  edges_r_x_b_ratio_angle_idx,
#  Ynet_rows_Idxs_in_flattend,
#  Ynet_real_imag_Idxs_in_flattend,
 
#  scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
#  scale_Pg_Png_Qng_Idx,
#  dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

#  sta_pf_PQ_para,
#  generic_gens_para,
#  ode_gens_para ) =
#      NamedTupleTools.select(
#          getproperty(
#          a_status_steady_state_data,
#         :static_prefault_paras),             
#          (:baseMVA,
#           :basekV,

#           :vh,
#           :θh,
#           :u0_model_states_init,
#           :model_syms,
#           :model_mass_matrix,

#           :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
#           :plants_cb_paras_switches,

#           :generic_system_dynamics_wt_fault_kwd_para,

#           :gens_nodes_names,
#           :SM_gens_nodes_names,
#           :non_gens_nodes_names,

#           :cb_states,
          
#           :dyn_pf_fun_kwd_n2s_idxs,
#           :dyn_pf_fun_kwd_net_idxs,
          
#           :Pg_Qg_Png_Qng_Pll_Qll_Idx,
#           :Png_Qng_Pll_Qll_Idx,
#           :Pg_Png_Qng_Idx,
          
#           :pf_vh_θh_idx_and_idx2Idx,
          
#           :dyn_pf_flat_vh_flat_θh_Idx,

#           :edges_r,
#           :edges_x,
#           :edges_b,
#           :edges_ratio,
#           :edges_angle,
#           :edges_type,
#           :Gs,
#           :Bs,
          
#           :Ynet_wt_nodes_idx_wt_adjacent_nodes,
#           :Ybr_cal_and_edges_orientation,
          
#           :Pg_Qg_Png_Qng_Pll_Qll,
#           :loc_load_exist,
          
#           :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

#           :edges_r_x_b_ratio_angle_idx,
#           :Ynet_rows_Idxs_in_flattend,
#           :Ynet_real_imag_Idxs_in_flattend,

#           :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
#           :scale_Pg_Png_Qng_Idx,
#           :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

#           :sta_pf_PQ_para,
#           :generic_gens_para,
#           :ode_gens_para))



# (non_gens_vh_idx,
#  non_slack_gens_θh_idx,
#  non_gens_θh_idx,

#  red_vh_Idxs,
#  red_non_slack_gens_θh_idx2Idx,
#  red_non_gens_θh_idx2Idx,

#  red_non_gens_vh_Idxs,
#  red_non_slack_gens_θh_Idxs,
#  red_non_gens_θh_Idxs,

#  red_slack_value_Idxs,

#  gens_nodes_idx,

#  non_gens_nodes_idx,
#  non_slack_gens_and_non_gens_idx,
#  gens_nodes_with_loc_loads_idx,
#  all_nodes_idx,
 
#  n2s_gens_idx,
#  n2s_gens_with_loc_load_idxs,
#  n2s_non_gens_idx,
#  n2s_all_nodes_idx ) =
#      NamedTupleTools.select(
#          pf_vh_θh_idx_and_idx2Idx,
#          (:non_gens_vh_idx,
#           :non_slack_gens_θh_idx,
#           :non_gens_θh_idx,

#           :red_vh_Idxs,
#           :red_non_slack_gens_θh_idx2Idx,
#           :red_non_gens_θh_idx2Idx,

#           :red_non_gens_vh_Idxs,
#           :red_non_slack_gens_θh_Idxs,
#           :red_non_gens_θh_Idxs,

#           :red_slack_value_Idxs,

#           :gens_nodes_idx,
          
#           :non_gens_nodes_idx,
#           :non_slack_gens_and_non_gens_idx,
#           :gens_nodes_with_loc_loads_idx,
#           :all_nodes_idx,

#           :n2s_gens_idx,
#           :n2s_gens_with_loc_load_idxs,
#           :n2s_non_gens_idx,
#           :n2s_all_nodes_idx))



# #--------------------------------------------
# #--------------------------------------------


# (;Ynet,
#  nodes_idx_with_adjacent_nodes_idx) =
#      NamedTupleTools.select(
#          Ynet_wt_nodes_idx_wt_adjacent_nodes,
#          (:Ynet,
#           :nodes_idx_with_adjacent_nodes_idx))


# (;edges_Ybr_cal,
#  edges_orientation) =
#      NamedTupleTools.select(
#          Ybr_cal_and_edges_orientation,
#          (:edges_Ybr_cal,
#           :edges_orientation ))

# #--------------------------------------------
# #--------------------------------------------

# (;dyn_P_gens_Idxs,
#  dyn_Q_gens_Idxs,
#  dyn_P_non_gens_Idxs,
#  dyn_Q_non_gens_Idxs,
#  dyn_P_gens_loc_load_Idxs,
#  dyn_Q_gens_loc_load_Idxs) =
#      NamedTupleTools.select(
#          Pg_Qg_Png_Qng_Pll_Qll_Idx,
#          (:dyn_P_gens_Idxs,
#           :dyn_Q_gens_Idxs,
#           :dyn_P_non_gens_Idxs,
#           :dyn_Q_non_gens_Idxs,
#           :dyn_P_gens_loc_load_Idxs,
#           :dyn_Q_gens_loc_load_Idxs))


# #--------------------------------------------

# red_vh_θh_idx =
#     getproperty(
#         pf_vh_θh_idx_and_idx2Idx,
#         :red_vh_θh_idx)

# #--------------------------------------------

# (;Ynet_rows_Idxs_in_flattend,
#  Ynet_real_imag_Idxs_in_flattend ) =
#      NamedTupleTools.select(
#          get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
#              Ynet),
#          (:Ynet_rows_Idxs_in_flattend,
#           :Ynet_real_imag_Idxs_in_flattend))

# (;Ynet_rows_Idxs_in_flattend,
# Ynet_real_imag_Idxs_in_flattend ) =
#     NamedTupleTools.select(
#         get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
#             getproperty(
#                 Ynet_wt_nodes_idx_wt_adjacent_nodes,
#                 :Ynet)),
#         (:Ynet_rows_Idxs_in_flattend,
#          :Ynet_real_imag_Idxs_in_flattend) )

# (;Ynet_real_Idxs,
#  Ynet_imag_Idxs) =
#     NamedTupleTools.select(
#         Ynet_real_imag_Idxs_in_flattend,
#         (:Ynet_real_Idxs,
#          :Ynet_imag_Idxs))

# #--------------------------------------------


# (;r_Idxs, x_Idxs, b_Idxs,
#  ratio_Idxs, angle_Idxs) =
#      NamedTupleTools.select(
#          edges_r_x_b_ratio_angle_idx,
#          (:r_Idxs, :x_Idxs, :b_Idxs,
#           :ratio_Idxs, :angle_Idxs))


# (vh_Idxs, θh_Idxs) =
#     NamedTupleTools.select(
#         dyn_pf_flat_vh_flat_θh_Idx,
#         (:dyn_pf_vh_Idxs,
#          :dyn_pf_θh_Idxs))


# #--------------------------------------------
# #--------------------------------------------

# Pg_inj_Qg_inj_Png_Qng_Pll_Qll =
#     get_Pg_inj_Qg_inj_Png_Qng_Pll_Qll(
#         Pg_Qg_Png_Qng_Pll_Qll;
#         loc_load_exist,
#         Pg_Qg_Png_Qng_Pll_Qll_Idx,    
#         gens_nodes_idx,
#         n2s_gens_idx,
#         n2s_gens_with_loc_load_idxs)


# Pg_inj_Qg_inj_Png_Qng =
#     getproperty(
#         Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
#         :Pg_inj_Qg_inj_Png_Qng)


# (;gens_Pg_inj,
#  gens_Qg_inj,
#  P_non_gens,
#  Q_non_gens,
#  Pg_inj_Qg_inj_Png_Qng) =
#      NamedTupleTools.select(
#          Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
#          (:gens_Pg_inj,
#           :gens_Qg_inj,
#           :P_non_gens,
#           :Q_non_gens,
#           :Pg_inj_Qg_inj_Png_Qng))


# #--------------------------------------------
# #--------------------------------------------


# Pg =
#     Pg_Qg_Png_Qng_Pll_Qll[dyn_P_gens_Idxs]

# Qg =
#     Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_gens_Idxs]

# Png =
#     Pg_Qg_Png_Qng_Pll_Qll[dyn_P_non_gens_Idxs]
 
# Qng =
#     Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_non_gens_Idxs]

# Pll =
#     Pg_Qg_Png_Qng_Pll_Qll[dyn_P_gens_loc_load_Idxs]

# Qll =
#     Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_gens_loc_load_Idxs]


# #--------------------------------------------
# # Loss participation factors
# #--------------------------------------------

# gens_Sn =
#     getproperty(
#         ode_gens_para,
#         :Sn)


# gens_loss_participation =
#     get_loss_particpation_by_gens_rating(
#         gens_Sn, Pg)

# gens_loss_participation_by_loading =
#     get_loss_particpation_by_gens_loading( Pg)

# #--------------------------------------------
# # power disturbance resolution participation
# #--------------------------------------------

# gens_active_power_particpation_by_rating =
#     get_gens_active_power_particpation_by_rating(
#         gens_Sn, Pg)

# gens_active_power_particpation_by_loading =
#     get_gens_active_power_particpation_by_loading( Pg )

# gens_reactive_power_particpation_by_rating =
#     get_gens_reactive_power_particpation_by_rating(
#         gens_Sn, Qg)

# gens_reactive_power_particpation_by_loading =
#     get_gens_reactive_power_particpation_by_loading( Qg)


# active_power_disturbance_resolution_participation =
#     gens_active_power_particpation_by_rating


# reactive_power_disturbance_resolution_participation =
#     gens_reactive_power_particpation_by_loading


# #--------------------------------------------
# # kwd parameters
# #--------------------------------------------

# pf_model_kwd_para =
#     (;loc_load_exist,
#      slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     
#      Ynet_wt_nodes_idx_wt_adjacent_nodes,

#      gens_loss_participation,
#      active_power_disturbance_resolution_participation,
#      reactive_power_disturbance_resolution_participation,

#      sta_pf_PQ_para,     

#      dyn_pf_fun_kwd_n2s_idxs,
#      dyn_pf_fun_kwd_net_idxs,
     
#      Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      Pg_Png_Qng_Idx,
#      pf_vh_θh_idx_and_idx2Idx,
     
#      scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      scale_Pg_Png_Qng_Idx,
#      dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

#      edges_type,
#      edges_orientation,
     
#      edges_r_x_b_ratio_angle_idx,
#      Ynet_rows_Idxs_in_flattend,
#      Ynet_real_imag_Idxs_in_flattend
#      )


# generic_sol_kwd_para =
#     (;
#      Ybr_cal_and_edges_orientation,
#      ode_gens_para,
#      sta_pf_PQ_para,
#      pf_model_kwd_para )


# #--------------------------------------------
# #--------------------------------------------
# # Initial values and parameters
# #--------------------------------------------
# #--------------------------------------------


# vh_θh =
#     [vh;
#      θh]



# vh_θh_slack_value =
#     [vh_θh;
#      [0.0]]

# red_vh_θh_slack_value =
#     [vh_θh[non_gens_vh_idx];
#      vh_θh[non_slack_gens_θh_idx];
#      vh_θh[non_gens_θh_idx];
#      [0.0]]


# non_gens_vh_all_θh_slack_value =
#     [vh[non_gens_nodes_idx];
#      θh;
#      [0.0]]    

# #--------------------------------------------

# Pg_inj_Qg_inj_Png_Qng =
#     [gens_Pg_inj;
#      gens_Qg_inj;
#      P_non_gens;
#      Q_non_gens ]

# Pg_inj_Png_Qng =
#     [gens_Pg_inj;
#      P_non_gens;
#      Q_non_gens ]

# slack_gens_nodes_idx =
#     getproperty(
#         dyn_pf_fun_kwd_net_idxs,
#         :slack_gens_nodes_idx)[1]

# # Subtract net loss from slack node Pg,
# # loss will subsequently be distributed
 
# ds_gens_Pg_inj = @set gens_Pg_inj[slack_gens_nodes_idx] =
#     gens_Pg_inj[slack_gens_nodes_idx] -
#     get_total_P_network_loss(
#         vh, θh, Ynet; nodes_idx_with_adjacent_nodes_idx,
#         n2s_all_nodes_idx,
#         all_nodes_idx )


# ds_Pg_inj_Qg_inj_Png_Qng =
#     [ds_gens_Pg_inj;
#      gens_Qg_inj;
#      P_non_gens;
#      Q_non_gens ]

# ds_Pg_inj_Png_Qng =
#     [ds_gens_Pg_inj;
#      P_non_gens;
#      Q_non_gens]

# ds_scale_Pg_inj_Qn_inj_Png_Qng =
#     [[1.0];
#      ds_Pg_inj_Qg_inj_Png_Qng]

# ds_scale_Pg_inj_Png_Qng =
#     [[1.0];
#      ds_Pg_inj_Png_Qng ]

# #--------------------------------------------
# #  net loss functions 
# #--------------------------------------------

# sum_line_losses =
#     get_sum_line_losses(
#         vh_θh,
#         edges_r_x_b_ratio_angle;    
#         edges_type,
#         edges_orientation,
#         n2s_all_nodes_idx,
#         dyn_pf_flat_vh_flat_θh_Idx,
#         edges_r_x_b_ratio_angle_idx)

# #--------------------------------------------


# total_P_network_loss_by_flattened_Ynet =
#     get_total_P_network_loss_by_flattened_Ynet(
#         vh_θh, Ynet_real_imag_flattend;
#         dyn_pf_flat_vh_flat_θh_Idx,
#         nodes_idx_with_adjacent_nodes_idx,
#         Ynet_rows_Idxs_in_flattend,
#         Ynet_real_imag_Idxs_in_flattend,
#         n2s_all_nodes_idx,
#         all_nodes_idx )

# total_Q_network_loss_by_flattened_Ynet =
#     get_total_Q_network_loss_by_flattened_Ynet(
#         vh_θh, Ynet_real_imag_flattend;
#         dyn_pf_flat_vh_flat_θh_Idx,
#         nodes_idx_with_adjacent_nodes_idx,
#         Ynet_rows_Idxs_in_flattend,
#         Ynet_real_imag_Idxs_in_flattend,
#         n2s_all_nodes_idx,
#         all_nodes_idx )


# #--------------------------------------------
# #  vh θh by per line losses sensitivity 
# #--------------------------------------------

# lines_total_losses =
#     get_lines_total_losses(
#         vh,
#         θh;
#         Ybr_cal_and_edges_orientation,
#          n2s_all_nodes_idx)


# #--------------------------------------------

# #--------------------------------------------
# # Loss participation factor based on gens ratings
# #--------------------------------------------

# # Synchrnous condensers will not participate

# gens_loss_participation_factor =
#     getproperty(
#         gens_loss_participation,
#         :gens_loss_participation_factor)

# #--------------------------------------------


# net_P_loss = loc_load_exist == true ?
#     sum(Pg) - sum(Pll) - sum(Png) : sum(Pg) - sum(Png)


# net_Q_loss = loc_load_exist == true ?
#     sum(Qg) - sum(Qll) - sum(Qng) : sum(Qg) - sum(Qng)

# #--------------------------------------------
# #--------------------------------------------


# t_df_dx = ForwardDiff.gradient(
#     ( vh_θh_x ) ->
#         get_total_P_network_loss_by_flattened_Ynet(
#             vh_θh_x, Ynet_real_imag_flattend ;
#             dyn_pf_flat_vh_flat_θh_Idx,
#             nodes_idx_with_adjacent_nodes_idx,
#             Ynet_rows_Idxs_in_flattend,
#             Ynet_real_imag_Idxs_in_flattend,
#             n2s_all_nodes_idx,
#             all_nodes_idx ) ,
#     vh_θh )


# t_df_dp = ForwardDiff.gradient(
#     ( p ) ->
#         get_total_P_network_loss_by_flattened_Ynet(
#             vh_θh, p ;
#             dyn_pf_flat_vh_flat_θh_Idx,
#             nodes_idx_with_adjacent_nodes_idx,
#             Ynet_rows_Idxs_in_flattend,
#             Ynet_real_imag_Idxs_in_flattend,
#             n2s_all_nodes_idx,
#             all_nodes_idx ) ,
#     Ynet_real_imag_flattend )


# # t_df_dx_ 

# # dx_dp  = -(svd( t_df_dx )) \ t_df_dp


# #--------------------------------------------
# #  vh θh by per node losses sensitivity 
# #--------------------------------------------

# n_df_dx = ForwardDiff.jacobian(
#     ( vh_θh_x ) ->
#         get_losses_per_node_by_flattened_Ynet(
#             vh_θh_x, Ynet_real_imag_flattend ;
#             dyn_pf_flat_vh_flat_θh_Idx,
#             nodes_idx_with_adjacent_nodes_idx,
#             Ynet_rows_Idxs_in_flattend,
#             Ynet_real_imag_Idxs_in_flattend,
#             n2s_all_nodes_idx,
#             all_nodes_idx ),
#     vh_θh )

# n_df_dp = ForwardDiff.jacobian(
#     ( p ) ->
#         get_losses_per_node_by_flattened_Ynet(
#             vh_θh, p ;
#             dyn_pf_flat_vh_flat_θh_Idx,
#             nodes_idx_with_adjacent_nodes_idx,
#             Ynet_rows_Idxs_in_flattend,
#             Ynet_real_imag_Idxs_in_flattend,
#             n2s_all_nodes_idx,
#             all_nodes_idx ) ,
#     Ynet_real_imag_flattend )

# n_dx_dp  = -( svd( n_df_dx ) \ n_df_dp )



# df_dx = ForwardDiff.jacobian(
#     ( vh_θh_x ) ->
#         get_losses_per_line(
#             vh_θh_x,
#             edges_r_x_b_ratio_angle;    
#             edges_type,
#             edges_orientation,
#             n2s_all_nodes_idx,
#             dyn_pf_flat_vh_flat_θh_Idx,
#             edges_r_x_b_ratio_angle_idx),
#             vh_θh )

# df_dp = ForwardDiff.jacobian(
#     ( p ) ->
#         get_losses_per_line(
#             vh_θh,
#             p;    
#             edges_type,
#             edges_orientation,
#             n2s_all_nodes_idx,
#             dyn_pf_flat_vh_flat_θh_Idx,
#             edges_r_x_b_ratio_angle_idx) ,
#     edges_r_x_b_ratio_angle )

# dx_dp  = -( svd( df_dx ) \ df_dp )

# #--------------------------------------------


# network_cases_dynamic =
#     ["case4gs",
#      "case5",
#      # "case5mul",
#      "case9",
#      "case14",
#      "case30",
#      "case39",
#      "case57",
#      "case118x",
#      "case300",
#      "case1354pegase",
#      "case2869pegase",
#      "case9241pegase",
#      "case13659pegase"
#      ]


# #-------------------------------
# #-------------------------------

# # (;plant_generators_data_from_json,
# #  plant_loads_data_from_json,
# #  plant_transmission_data_from_json,
# #  edge_data_from_json,
# #  shunt_data_from_json,
# #  baseMVA,
# #  gencost_data_from_json,

# #  edges_fbus,
# #  edges_tbus,
# #  edges_r,
# #  edges_x,
# #  edges_b,
# #  edges_ratio,
# #  edges_angle,
# #  edges_type,

# #  shunt_idx,
# #  Gs,
# #  Bs,

# #  edges_orientation,
# #  edges_Ybr_cal,
# #  Ybr_cal_and_edges_orientation,
# #  Ynet_wt_nodes_idx_wt_adjacent_nodes,
# #  Ynet_rows_Idxs_in_flattend,
# #  Ynet_real_imag_Idxs_in_flattend ,

# #  P_gens,
# #  Q_gens,
# #  P_non_gens,
# #  Q_non_gens,
# #  P_g_loc_load,
# #  Q_g_loc_load,
# #  loc_load_exist,

# #  generic_gens_para,
# #  generic_each_gen_para,

# # # pf_generic_gens_para,

# # net_nodes_type_idxs,

# # dyn_pf_fun_kwd_net_idxs,
# # dyn_pf_fun_kwd_n2s_idxs,

# # nodes_with_demands_idx,

# # all_nodes_idx,
# # n2s_all_nodes_idx,

# # sta_pf_PQ_para,
# # gens_vh_slack_θh_para,

# # sta_pf_vars_and_paras_idx,
# # pf_sta_ΔPQ_mismatch_parameters,

# # sta_red_vh_θh_0,

# # pf_PQ_param,
# # pf_kw_para,

# # kwd_sta_sta_ΔPQ_sol_by_json,

# # generic_red_sol_kwd_para,
# # # generic_dyn_sol_kwd_para,
           
# #  static_Idx_and_syms ) =
# #      NamedTupleTools.select(
# #          system_net_static_data,
# #          (:plant_generators_data_from_json,
# #           :plant_loads_data_from_json,
# #           :plant_transmission_data_from_json,
# #            :edge_data_from_json,
# #            :shunt_data_from_json,
# #            :baseMVA,
# #            :gencost_data_from_json,

# #            :edges_fbus,
# #            :edges_tbus,
# #            :edges_r,
# #            :edges_x,
# #            :edges_b,
# #            :edges_ratio,
# #            :edges_angle,
# #            :edges_type,

# #            :shunt_idx,
# #            :Gs,
# #            :Bs,
           
# #            :edges_orientation,
# #            :edges_Ybr_cal,
# #            :Ybr_cal_and_edges_orientation,
# #            :Ynet_wt_nodes_idx_wt_adjacent_nodes,
# #            :Ynet_rows_Idxs_in_flattend,
# #            :Ynet_real_imag_Idxs_in_flattend ,

# #            :P_gens,
# #            :Q_gens,
# #            :P_non_gens,
# #            :Q_non_gens,
# #            :P_g_loc_load,
# #            :Q_g_loc_load,
# #            :loc_load_exist,

# #            :generic_gens_para,
# #            :generic_each_gen_para,

# #            # :pf_generic_gens_para,

# #            :net_nodes_type_idxs,
           
# #            :dyn_pf_fun_kwd_net_idxs,
# #            :dyn_pf_fun_kwd_n2s_idxs,

# #            :nodes_with_demands_idx,

# #            :all_nodes_idx,
# #            :n2s_all_nodes_idx,

# #            :sta_pf_PQ_para,
# #            :gens_vh_slack_θh_para,

# #            :sta_pf_vars_and_paras_idx,
# #            :pf_sta_ΔPQ_mismatch_parameters,

# #            :sta_red_vh_θh_0,

# #            :pf_PQ_param,
# #            :pf_kw_para,
           
# #            :kwd_sta_sta_ΔPQ_sol_by_json,
           
# #            :generic_red_sol_kwd_para,
# #            # :generic_dyn_sol_kwd_para,
           
# #            :static_Idx_and_syms) )


# #--------------------------------------------



# """

# #---------------------------------------------------
# # save julia object to latex file
# #---------------------------------------------------



# tuple_julia_object =
#     (;Yint,
#       Yred,
#       δg,
#       Eg,
#       Cinm,
#       Dinm,
#       Aω_matrix,
#       A_matrix)

# names_julia_object =
#     propertynames(tuple_julia_object)

# open(tex_filename, "a") do file_handle
    
#     for (name_object, a_julia_object) in
#         zip(names_julia_object,
#             tuple_julia_object)
        
#         write(file_handle,"\n $(String(name_object)) = ")
        
#         write(file_handle, latexify(
#             a_julia_object;
#             fmt=FancyNumberFormatter()))
#     end
    
# end




# pf_tuple_julia_object =
#     (;pf_gens_results,
#       pf_non_gens_results)

# pf_names_julia_object =
#     propertynames(pf_tuple_julia_object)

# open(tex_filename, "a") do file_handle
    
#     for (name_object, a_julia_object) in
#         zip(pf_names_julia_object,
#             pf_tuple_julia_object)
        
#         write(file_handle,"\n $(String(name_object)) = ")
        
#         write( file_handle, latexify(
#             a_julia_object;
#             fmt=FancyNumberFormatter()))
#     end
    
# end

# """
