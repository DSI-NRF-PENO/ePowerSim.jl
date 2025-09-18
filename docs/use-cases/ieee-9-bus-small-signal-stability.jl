#=
# [Example of Small-signal Stability Analysis Flow with IEEE 9 Bus Test System](@id ieee-9-bus-small-signal-stability)
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

using DifferentialEquations

using OrdinaryDiffEq, Sundials, ODEInterfaceDiffEq

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
# global settings
#---------------------------------------------------

freq = 60

Ωb = 2 * pi * freq

ωs = Ωb 

ω_ref0 = ωs

basekV = 1.0

#---------------------------------------------------
#---------------------------------------------------
# Reading network data
#---------------------------------------------------
#---------------------------------------------------

"""
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

"""


net_wt_avr_string = "net-static-data-avr-sauer"

gov_string        = "gov-sauer"

dynamic_net_data_by_components_file =
    "$(net_wt_avr_string)-$(gov_string).json"

json_net_data_by_components_file =
    dynamic_net_data_by_components_file

#---------------------------------------------------

case_name = "case9"

# This ensures a lowercase name
case_name = lowercase(case_name) 

#---------------------------------------------------

sim_type  = "$(case_name)-"*"small-signal-stability-"*
    "-$(net_wt_avr_string)-$(gov_string)"

#---------------------------------------------------

case_data_dir =
   joinpath( data_dir,
             case_name)

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
# Simulation Period
#---------------------------------------------------

timespan      = 20.0
    

time_start    = 0.0

time_final    = timespan

tspan         = (0.0, timespan)

sim_timespan  = (0.0, timespan)

plot_timespan = (0.0, timespan)

#---------------------------------------------------
## variables, parameters, indices, intialisation
#---------------------------------------------------

reduced_and_linearised_model_parameters =
    get_generic_reduced_and_linearised_model_parameters(
        net_data_by_components_file,
        # baseMVA,
        basekV,
        use_pu_in_PQ,
        line_data_in_pu,
        pf_alg;
        tt = 0.0,
        in_components_type_sym =
            false)


(;system_init,
 ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
 linearised_dynamic_kwd_para,
 generic_electro_mechanical_oscillation_indicies) =
    NamedTupleTools.select(
        reduced_and_linearised_model_parameters,
        (:system_init,
         :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
         :linearised_dynamic_kwd_para,
         :generic_electro_mechanical_oscillation_indicies))


(Yint,
 Yred,
 δg,
 Eg,
 Cinm,
 Dinm,
 Aω_matrix,
 A_matrix) =
     NamedTupleTools.select(
         generic_electro_mechanical_oscillation_indicies,
         (:Yint,
          :Yred,
          :δg,
          :Eg,
          :Cinm,
          :Dinm,
          :Aω_matrix,
          :A_matrix ))


tt = 0.0

(; Asys, ϕ, Q,
 JAE, JLF, Asys_dash,
 sys_eigvalues,
 printed_sys_eigvalues,
 sys_eigvalues_dash,
 printed_sys_eigvalues_dash,
     ) =
      NamedTupleTools.select(
         get_generic_Asys_linearised_dynamic_model(
             system_init,
             tt,
             ω_ref_v_ref_p_order_Png_Qng_Pll_Qll;
             use_init_condition =
                 true,
             kwd_para =
                 linearised_dynamic_kwd_para ),
          (:Asys, :ϕ, :Q,
           :JAE, :JLF, :Asys_dash,
           :sys_eigvalues,
           :printed_sys_eigvalues,
           :sys_eigvalues_dash,
           :printed_sys_eigvalues_dash ))


(PF_Asys,
 printed_PF_Asys,
 M_diag,
 printed_M_diag,
 eig_values,
 printed_eig_values,
 eigvecs_left,
 eigvecs_right,
 inv_eigvecs_right) =
     NamedTupleTools.select(
         get_generic_small_signal_stability_indices(
        Asys),
         (:PF_Asys,
          :printed_PF_Asys,
          :M_diag,
          :printed_M_diag,
          :eig_values,
          :printed_eig_values,
          :eigvecs_left,
          :eigvecs_right,
          :inv_eigvecs_right))


generic_small_signal_stability_indices_dash =
    get_generic_small_signal_stability_indices(
        Asys_dash)

(PF_Asys_dash,
 printed_PF_Asys_dash,
 M_diag_dash,
 printed_M_diag_dash,
 eig_values_dash,
 printed_eig_values_dash,
 eigvecs_left_dash,
 eigvecs_right_dash,
 inv_eigvecs_righ_dasht) =
     NamedTupleTools.select(
         get_generic_small_signal_stability_indices(
        Asys_dash),
         (:PF_Asys,
          :printed_PF_Asys,
          :M_diag,
          :printed_M_diag,
          :eig_values,
          :printed_eig_values,
          :eigvecs_left,
          :eigvecs_right,
          :inv_eigvecs_right))

#---------------------------------------------------
## plot eigen values
#---------------------------------------------------


plt_eigenvalues =
    plot(real(printed_eig_values),
         imag(printed_eig_values),
         seriestype = :scatter,
         label = "eigenvalues",
         # title = "Eigenvalues of IEEE14",
         xlabel = "Real Part",
         ylabel = "Imaginary Part")

savefig(plt_eigenvalues,
        joinpath(figure_dir,
             "$(case_name)-" *
                 "$(sim_type)-" *
                 "eig-values.pdf"))

###

x_part = real.(printed_eig_values)
y_part = imag.(printed_eig_values)
z_part = norm.(eigvecs_left)


s_plt_eigenvalues =
    scatter(
        x_part,
        y_part,
        marker_z=z_part,
        markercolors=:blues,
        label = "eigenvalues",
        xlabel = "eigenvalues real part",
        ylabel = "eigenvalues imaginary part")

savefig(s_plt_eigenvalues,
        joinpath(figure_dir,
             "$(case_name)-" *
                 "$(sim_type)-" *
                 "eig-values-heatmap.pdf"))

#---------------------------------------------------
# save julia object to latex file
#---------------------------------------------------


tuple_julia_object =
    (;Yint,
      Yred,
      δg,
      Eg,
      Cinm,
      Dinm,
      Aω_matrix,
      A_matrix)

names_julia_object =
    propertynames(tuple_julia_object)

open(tex_filename, "a") do file_handle
    
    for (name_object, a_julia_object) in
        zip(names_julia_object,
            tuple_julia_object)
        
        write(file_handle,"\n $(String(name_object)) = ")
        
        write(file_handle, latexify(
            a_julia_object;
            fmt=FancyNumberFormatter()))
    end
    
end


#---------------------------------------------------
#-------------------------------------------------------

# (t2_printed_sys_eigvalues,
#  ) =
#       NamedTupleTools.select(
#          get_generic_linearised_dynamic_model(
#              system_init,
#              ω_ref_v_ref_p_order_Png_Qng_Pll_Qll;
#              kwd_para =
#                  linearised_dynamic_kwd_para ),
#          (:printed_sys_eigvalues, ))


    
# (t3_printed_sys_eigvalues,
#  ) =
#       NamedTupleTools.select(
#          get_generic_linearised_dynamic_model(
#              system_init,
#              tt,
#              ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
#              nothing;
#              use_init_condition =
#                  true,
#              kwd_para =
#                  linearised_dynamic_kwd_para ),
#           (:printed_sys_eigvalues, ))

#---------------------------------------------------


# (t1_printed_sys_eigvalues,
#  t1_printed_sys_eigvalues_dash
#      ) =
#       NamedTupleTools.select(
#           get_generic_linearised_dynamic_model(
#               system_init,
#               tt,
#               ω_ref_v_ref_p_order_Png_Qng_Pll_Qll;
#               use_init_condition =
#                   true,
#               kwd_para =
#                   linearised_dynamic_kwd_para ),
#           (:printed_sys_eigvalues,
#            :printed_sys_eigvalues_dash))


# str1 = latexify(
#      δg;
#     fmt=FancyNumberFormatter())

# str2 = latexify(
#     Yint;
#     fmt=FancyNumberFormatter())

# str3 = latexify(
#     Yint;
#     fmt=SiunitxNumberFormatter(
#         format_options=
#             "round-mode=places,round-precision=1",
#         version=3))

# open(tex_filename, "a") do file_handle
    
#     for a_julia_object in list_julia_object
        
#         write(file_handle,"\n")
        
#         write(file_handle, latexify(
#             a_julia_object;
#             fmt=x->round(x, sigdigits=4)))
#     end
    
# end


# tex_file_hanle = open(tex_filename, "a")

# write(tex_file_hanle, " Hello World , Julia welcomes you")

# close(tex_file_hanle)


#-------------------------------------------------------
#-------------------------------------------------------
