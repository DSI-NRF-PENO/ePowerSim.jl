
#---------------------------------------------------
# global settings
#---------------------------------------------------

freq = 60

Ωb = 2 * pi * freq

ωs = Ωb 

ω_ref0 = ωs

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


# using PDMats

#---------------------------------------------------

Random.seed!(123)

#---------------------------------------------------
# base setting and some booleans 
#---------------------------------------------------

use_pu_in_PQ    = true

line_data_in_pu = true

with_faults     = false

#---------------------------------------------------
# Solvers
#---------------------------------------------------

pf_alg            = NewtonRaphson()

ode_alg           = Rodas4()

# ode_alg          = ImplicitMidpoint()
  
dae_alg           = IDA()


dt                = 0.01

Δt                = 1.0 / 2^(4)


abstol            = 1e-12

reltol            = 1e-12

fractional_digits = 6

#---------------------------------------------------
#---------------------------------------------------
## simulation cases
#---------------------------------------------------
#---------------------------------------------------

sauer_net_wt_avr_string    = "net-static-data-avr-sauer-"

rtds_net_wt_avr_string     = "net-static-data-avr-rtds-"

sauer_gov_string           = "gov-sauer"

ieee_tgov2sauer_gov_string = "gov-ieee-tgov2sauer"

ieee_tgov1_gov_string      = "gov-ieee-tgov1"

#---------------------------------------------------

net_wt_avr_string = rtds_net_wt_avr_string

gov_string        = ieee_tgov1_gov_string

json_net_data_by_components_file =
    "$(net_wt_avr_string)"*
    "$(gov_string)" *
    ".json"

#---------------------------------------------------
## simulation case
#---------------------------------------------------

# case_name = "case9"

case_name = "case14"

#---------------------------------------------------

json_case_dir =
    joinpath(
        data_dir,
        case_name,
        "json")


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

case_file = joinpath(data
    "matpower-data",
    "$(case_name).m" )

#---------------------------------------------------

sim_type = "optimial-powerflow"

#---------------------------------------------------

cd( script_dir )

results_dir =
    joinpath(
        script_dir, 
        "results", sim_type)

if !(isdir(results_dir))
    
    mkpath(results_dir)
    
end

figure_dir =
    joinpath(
        script_dir,
        "figure", sim_type)

if !(isdir(figure_dir))
    
    mkpath(figure_dir)
    
end

#---------------------------------------------------

tex_filename =
    joinpath(results_dir,
             "$(case_name)-" *
                 "$(sim_type)-" *
                 "validation.tex")

#---------------------------------------------------
# generic system simulation parameters 
#---------------------------------------------------

generic_system_simulation_parameters =
    get_generic_system_simulation_parameters(
        net_data_by_components_file;
        components_libs_dir,
        basekV          = 1.0,    
        use_pu_in_PQ    = true,
        line_data_in_pu = true,
        with_faults     = false )

#---------------------------------------------------
# Optimal powerflow generic system simulation parameters
#---------------------------------------------------

(df_opf,
 opf_objval_solution,
 opf_status) =
     NamedTupleTools.select(
         get_opf_wt_generic_system_simulation_parameters(
             net_data_by_components_file;
             components_libs_dir ),
         (:df_opf,
          :objval_solution,
          :status))


opf_sim_type = "opf"

opf_sim_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "$(sim_type)-" *
            "$(opf_sim_type).csv")



CSV.write(opf_sim_csv_filename,
          df_opf )



open( tex_filename , "a") do file_handle

    write(file_handle,"\n" )

    write(file_handle," Optimial power flow type: $(opf_sim_type)" )

    write(file_handle,"\n" )
    
    write(file_handle,"Objective value $(opf_objval_solution)" )
    
    write(file_handle,"\n" )
    
    write(file_handle,"Solution status $(opf_status)" )

    write(file_handle,"\n" )
        
    write(file_handle,
          get_df2tex(
              df_opf) )
    
end

#---------------------------------------------------
# Optimal powerflow relaxation generic system simulation parameters
#---------------------------------------------------

(df_opf_by_relaxation,
 sdp_relaxation_lower_bound,
 pf_by_relaxation_status ) =
     NamedTupleTools.select(
         get_opf_by_relaxation_wt_generic_system_simulation_parameters(
             net_data_by_components_file;
             components_libs_dir ),
         (:df_opf_by_relaxation,
          :sdp_relaxation_lower_bound,
          :status ))


opf_relax_sim_type = "opf-relaxation"


opf_relax_sim_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "$(sim_type)-" *
            "$(opf_relax_sim_type).csv")



CSV.write(opf_relax_sim_csv_filename ,
          df_opf_by_relaxation )


open( tex_filename , "a") do file_handle

    write(file_handle,"\n" )

    write(file_handle," Optimial power flow type: $(opf_relax_sim_type)" )

    write(file_handle,"\n" )
    
    write(file_handle,"Objective value $(sdp_relaxation_lower_bound)" )
    
    write(file_handle,"\n" )
    
    write(file_handle,"Solution status $(pf_by_relaxation_status)" )

    write(file_handle,"\n" )
        
    write(file_handle,
          get_df2tex(
              df_opf_by_relaxation) )
    
end


#----------------------------------------
#----------------------------------------
# old work
#----------------------------------------
#----------------------------------------    

#############################################################################

#----------------------------------------
# opf
#----------------------------------------    

# Bus = all_nodes_idx

df_opf = opf( case_file )

df_opf[!, :] =
    round.(
        df_opf[:, :],
        digits=fractional_digits)

opf_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "sim-opf.csv")

CSV.write(opf_csv_filename,
          df_opf )

#----------------------------------------
# opf by relaxation
#----------------------------------------    

df_opf_by_relaxation =
        opf_by_relaxation( case_file )

df_opf_by_relaxation[!, :] =
    round.(
        df_opf_by_relaxation[:, :],
        digits=fractional_digits)

opf_by_relaxation_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "sim-opf-by-relaxation.csv")

CSV.write(opf_by_relaxation_csv_filename,
          df_opf_by_relaxation )



#----------------------------------------
# opf by line loss
#----------------------------------------    

(;slack_gens_nodes_idx,
 non_slack_gens_nodes_idx,
 gens_nodes_idx,
 non_gens_nodes_idx,
 gens_nodes_with_loc_loads_idx,
 all_nodes_idx,

 nodes_with_demands_idx,

 dyn_pf_fun_kwd_net_idxs,

 no_nodes,
 no_gens,
 no_edges,

 gens_cost_coeff_ascen,
 gens_cost_coeff_decen,

 gens_installed_capacity,

 gens_Pmax,
 gens_Pmin,
 gens_Qmax,
 gens_Qmin,

 sch_Pg,
 nodes_Pd,
 nodes_Qd,

 P_Gen_lb,
 P_Gen_ub,
 Q_Gen_lb,
 Q_Gen_ub,

 P_Demand,
 Q_Demand,
 S_Demand,

 branch_fbus,
 branch_tbus,
 branch_r,
 branch_x,
 branch_b,

 edges_orientation,

 branch_data,

 Ybus) =
     NamedTupleTools.select(
         get_opf_net_optimisation_parameters(
             case_file ),
         (:slack_gens_nodes_idx,
          :non_slack_gens_nodes_idx,
          :gens_nodes_idx,
          :non_gens_nodes_idx,
          :gens_nodes_with_loc_loads_idx,
          :all_nodes_idx,

          :nodes_with_demands_idx,

          :dyn_pf_fun_kwd_net_idxs,

          :no_nodes,
          :no_gens,
          :no_edges,

          :gens_cost_coeff_ascen,
          :gens_cost_coeff_decen,

          :gens_installed_capacity,

          :gens_Pmax,
          :gens_Pmin,
          :gens_Qmax,
          :gens_Qmin,

          :sch_Pg,
          :nodes_Pd,
          :nodes_Qd,

          :P_Gen_lb,
          :P_Gen_ub,
          :Q_Gen_lb,
          :Q_Gen_ub,

          :P_Demand,
          :Q_Demand,
          :S_Demand,

          :branch_fbus,
          :branch_tbus,
          :branch_r,
          :branch_x,
          :branch_b,

          :edges_orientation,

          :branch_data,

          :Ybus))


#----------------------------------------
# opf
#----------------------------------------    

# Bus = all_nodes_idx

df_opf = opf( case_file )

df_opf[!, :] =
    round.(
        df_opf[:, :],
        digits=fractional_digits)

opf_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "sim-opf.csv")

CSV.write(opf_csv_filename,
          df_opf )

#----------------------------------------
# opf by relaxation
#----------------------------------------    

df_opf_by_relaxation =
        opf_by_relaxation( case_file )

df_opf_by_relaxation[!, :] =
    round.(
        df_opf_by_relaxation[:, :],
        digits=fractional_digits)

opf_by_relaxation_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "sim-opf-by-relaxation.csv")

CSV.write(opf_by_relaxation_csv_filename,
          df_opf_by_relaxation )


#----------------------------------------
# opf by line loss
#----------------------------------------    

result_opf_by_line_loss =
    opf_by_line_loss( case_file )

#---------------------------------------------------

# If = Yf * V

# It = Yt * V

# Ibranches = If + It

# # I^2 * R

# Line_loses =
#     (abs.(Ibranches)).^2 .* branch_r

# Line_loses =
#     (abs.(Ibranches)).^2 .* real.(Zbranch)

#---------------------------------------------------

opf_net_optimisation_parameters =
    get_opf_net_optimisation_parameters(
        case_file )


(;
 # gens_nodes_idx,
 # all_nodes_idx,
 # dyn_pf_fun_kwd_net_idxs,

 # gens_installed_capacity,
 
 # gens_Pmax,
 # gens_Pmin,
 # gens_Qmax,
 # gens_Qmin,

 # sch_Pg,
 nodes_Pd,
 nodes_Qd,
 
 gens_cost_coeff_ascen) =
    NamedTupleTools.select(
        opf_net_optimisation_parameters,
        (# :gens_nodes_idx,
         # :all_nodes_idx,
         # :dyn_pf_fun_kwd_net_idxs,

         # :gens_installed_capacity,

         # :gens_Pmax,
         # :gens_Pmin,
         # :gens_Qmax,
         # :gens_Qmin,

         # :sch_Pg,
         :nodes_Pd,
         :nodes_Qd,

         :gens_cost_coeff_ascen))

#---------------------------------------------------

opf_by_para =
    opf(
        nodes_Pd,
        nodes_Qd;
        opf_net_optimisation_parameters =
            opf_net_optimisation_parameters,
         round_digits = 4 )

opf_by_relaxation_by_para =
    opf_by_relaxation(
        nodes_Pd,
        nodes_Qd;
        opf_net_optimisation_parameters =
            opf_net_optimisation_parameters,
        round_digits = 4)

#---------------------------------------------------

(opf_by_s_df_header,
 opf_by_scenario_df) =
    get_opf_by_scenario(
    case_file )

opf_by_scenario_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "opf-by-scenario.csv")

CSV.write(opf_by_scenario_csv_filename,
          opf_by_scenario_df )

opf_by_s_header_idx_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "opf-by-sc-header-idx.csv")

CSV.write(opf_by_s_header_idx_csv_filename,
          DataFrame(OrderedDict(
              :idx =>
                  1:length(opf_by_s_df_header),
              :parameters =>
                  opf_by_s_df_header)) )

#---------------------------------------------------

(relx_opf_df_header,
 relx_opf_by_scenario_df ) =
    get_opf_by_relaxation_by_scenario(
    case_file )

opf_by_relx_scenario_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "opf-by-relx-scenario.csv")

CSV.write(opf_by_relx_scenario_csv_filename,
          relx_opf_by_scenario_df )

opf_by_relx_s_header_idx_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "opf-by-relx-sc-header-idx.csv")

CSV.write(opf_by_relx_s_header_idx_csv_filename,
          DataFrame(OrderedDict(
              :idx =>
                  1:length(relx_opf_df_header),
              :parameters =>
                  relx_opf_df_header)) )

#---------------------------------------------------
# Save results to files
#---------------------------------------------------

wind_gens_cost_scale  = 0.98    
solar_gens_cost_scale = 0.90
    
wind_gens_capacity_scale  = 0.001    
solar_gens_capacity_scale = 0.05

active_power_demand_deviation_scale   = 0.5
reactive_power_demand_deviation_scale = 0.5

(economic_dispatch_header,
 df_economic_dispatch_by_scenario) =
    get_economic_dispatch_by_scenario(
        case_file,

        wind_gens_cost_scale,    
        solar_gens_cost_scale,

        wind_gens_capacity_scale,    
        solar_gens_capacity_scale,

        active_power_demand_deviation_scale,
        reactive_power_demand_deviation_scale;
        no_scenarios = 100 )

ed_by_scenario_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "ed-by-scenario.csv")

CSV.write(ed_by_scenario_csv_filename,
          df_economic_dispatch_by_scenario )

economic_dispatch_header_idx_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "ed-header-sym-idx.csv")

CSV.write(economic_dispatch_header_idx_csv_filename,
          DataFrame(OrderedDict(
              :idx =>
                  1:length(economic_dispatch_header),
              :parameters =>
                  economic_dispatch_header)) )

#---------------------------------------------------

wind_gens_cost_scale  = 0.98    
solar_gens_cost_scale = 0.90
    
wind_gens_capacity_scale  = 0.001
solar_gens_capacity_scale = 0.05

active_power_demand_deviation_scale   = 0.5
reactive_power_demand_deviation_scale = 0.5

(unit_commitment_header,
 df_unit_commitment_by_scenario) =
    get_unit_commitment_by_scenario(
        case_file,

        wind_gens_cost_scale,    
        solar_gens_cost_scale,

        wind_gens_capacity_scale,    
        solar_gens_capacity_scale,

        active_power_demand_deviation_scale,
        reactive_power_demand_deviation_scale;
        no_scenarios = 100 )


uc_by_scenario_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "uc-by-scenario.csv")

CSV.write(uc_by_scenario_csv_filename,
          df_unit_commitment_by_scenario )


unit_commitment_header_idx_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "uc-header-sym-idx.csv")

CSV.write(unit_commitment_header_idx_csv_filename,
          DataFrame(OrderedDict(
              :idx =>
                  1:length(unit_commitment_header),
              :parameters =>
                  unit_commitment_header)) )

#---------------------------------------------------
#---------------------------------------------------

# Annualized capital cost
annualized_capital_cost_I = 90_000

# Operation cost per MWh
operation_cost_per_MWh_C  = 60

# Hours per year
no_hours_in_a_year_τ = 8_760

# # Number of Scenario
number_of_scenarios  = 5

# cal Scenario probabilities

# scenario_probabilities_θ =
#     [0.2, 0.2, 0.2, 0.2, 0.2]

# Utility function coefficients

utility_func_coefficients_A =
    [300, 350, 400, 450, 500]

# Utility function coefficients

utility_func_coefficients_B = 1

no_of_cases = 100

(mixed_complementarity_header,
 df_mixed_complementarity)  =
    get_mixed_complementarity(
        annualized_capital_cost_I,
        operation_cost_per_MWh_C,
        no_hours_in_a_year_τ,
        number_of_scenarios,
        # scenario_probabilities_θ,
        utility_func_coefficients_A,
        utility_func_coefficients_B;
        round_digits = 4,
        no_of_cases = 100)

mixed_complementarity_csv_filename =
    joinpath(
        results_dir,
        "mixed-complementarity.csv")

CSV.write(mixed_complementarity_csv_filename,
          df_mixed_complementarity )

mixed_complementarity_header_idx_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "mc-header-sym-idx.csv")

CSV.write(mixed_complementarity_header_idx_csv_filename,
          DataFrame(OrderedDict(
              :idx =>
                  1:length(mixed_complementarity_header),
              :parameters =>
                  mixed_complementarity_header)) )

#---------------------------------------------------
####################################################
#---------------------------------------------------
