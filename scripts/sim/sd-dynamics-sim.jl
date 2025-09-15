
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

list_fault_point_from_node_a = [0.01]

list_fault_resistance        = [0.001]

list_no_line_circuit         = [1]

list_edges_to_have_fault     = [ 2 ] # [ 4 ]

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
## simulation cases
#---------------------------------------------------
#---------------------------------------------------

# json_net_data_by_components_file =
#     "net-static-data-avr-sauer-gov-sauer.json"

# json_net_data_by_components_file =
#     "net-static-data-avr-rtds-gov-ieee-tgov2sauer.json"

# json_net_data_by_components_file =
#     "net-static-data-avr-rtds-gov-ieee-tgov1.json"

# json_net_data_by_components_file =
#     "net-static-data-avr-rtds-gov-sauer.json"

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


sim_type  = "sudden-load-change-validation"*
    "-$(net_wt_avr_string)-$(gov_string)"

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
             "$(case_name)-" *
                 "$(sim_type)-" *
                 "validation.tex")

#---------------------------------------------------

sd_dynamics_sim_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "$(sim_type)-" *
            "sim-sd-dynamics.csv")

#---------------------------------------------------
# sudden load change
#---------------------------------------------------

result_guick_sol =
    sim_sudden_load_change(
        ;case_name,
        script_dir,
        json_net_data_by_components_file,        
        sim_type,
        timespan=35,

        target_bus_name="bus5",

        pertubation_factor,
        restoration_factor,
        
        pertubation_time,
        restoration_time,
        
        Δt1,
        Δt2,

        fractional_digits=6,
        
        dt = 0.001,
        use_saveat = true,
        ts=0.001,

        ode_alg = RadauIIA5()
        # radau5()

        # RadauIIA5()
        # Trapezoid()
        # ImplicitEuler()
        # ImplicitMidpoint() # Rodas4()
    )


# Rodas4()
# ImplicitMidpoint(autodiff=false)
# ImplicitMidpoint()

 (;result_sol_plot,
  target_parameter_plot) =
      NamedTupleTools.select(
          result_guick_sol,
          (:result_sol_plot,
           :target_parameter_plot))

(;δ_a_plot,
 ω_a_plot,
 plot_gens_vh,
 plot_gens_θh,
 plot_non_gens_vh,
 plot_non_gens_θh) =
     NamedTupleTools.select(
         result_sol_plot ,
         (:δ_a_plot,
          :ω_a_plot,
          :plot_gens_vh,
          :plot_gens_θh,
          :plot_non_gens_vh,
          :plot_non_gens_θh ) )

#---------------------------------------------------
# Optimal powerflow generic system simulation parameters
#---------------------------------------------------

opf_wt_wt_generic_sys_sim_para =
    get_opf_wt_generic_system_simulation_parameters(
        net_data_by_components_file;
        components_libs_dir )


opf_by_relaxa_wt_generic_sys_sim_para =
    get_opf_by_relaxation_wt_generic_system_simulation_parameters(
        net_data_by_components_file;
        components_libs_dir )

#---------------------------------------------------
# Optimisal powerflow case file
#---------------------------------------------------

case_file = joinpath(data,
    "matpower-data",
    "$(case_name).m")

opf_net_optimisation_parameters =
    get_opf_net_optimisation_parameters(
        case_file )

#---------------------------------------------------
# opf
#---------------------------------------------------

(;df_opf,
 objval_solution,
 status) = NamedTupleTools.select(
     opf( case_file ),
     (:df_opf,
      :objval_solution,
      :status))

df_opf[!, :] =
    round.(
        df_opf[:, :],
        digits= fractional_digits )


open( tex_filename , "a") do file_handle

    write(file_handle,"\n" )

    write(file_handle,"objective value $(objval_solution)" )

    write(file_handle,"\n" )
        
    write(file_handle,
          get_df2tex(
              df_opf) )
    
    # write(file_handle,
    #       get_csv2tex(
    #           opf_csv_filename) )
    
end


opf_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "sim-opf.csv")

CSV.write(opf_csv_filename,
          df_opf )

#---------------------------------------------------
# opf by relaxation
#---------------------------------------------------

(;df_opf_by_relaxation,
 sdp_relaxation_lower_bound
 status) =
    NamedTupleTools.select(
        opf_by_relaxation( case_file ),
        (:df_opf_by_relaxation,
         :sdp_relaxation_lower_bound
         :status))

df_opf_by_relaxation[!, :] =
    round.(
        df_opf_by_relaxation[:, :],
        digits=fractional_digits)

open(tex_filename, "a") do file_handle

    write(file_handle,"\n" )

    write(file_handle,
          "objective value $(sdp_relaxation_lower_bound)" )

    write(file_handle,"\n" )
    
    write(file_handle,
          get_df2tex(
              df_opf_by_relaxation) )
    
end


opf_by_relaxation_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "sim-opf-by-relaxation.csv")

CSV.write(opf_by_relaxation_csv_filename,
          df_opf_by_relaxation )

#---------------------------------------------------
#---------------------------------------------------

# case 9

# Objective value (feasible solution)       : 5296.69

# Objective value (W & V relax. lower bound): 5298.15

#---------------------------------------------------

# case 14

# Objective value (feasible solution)       : 8137.07

# Objective value (W & V relax. lower bound): 8139.15

#---------------------------------------------------
#---------------------------------------------------

red_model_distributed_slack_pf =
    sim_red_model_distributed_slack_pf(
        ;system_status=
            :pre_fault_state,
        case_name=
            case_name,

        with_faults = false,

        json_net_data_by_components_file =
            json_net_data_by_components_file,

        data_dir=
            data_dir,

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

        basekV          = 1.0,    
        use_pu_in_PQ    = true,
        line_data_in_pu = true )


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


open(tex_filename, "a") do file_handle

    write(file_handle,"\n" )

    write(file_handle,"ds-slack-value = $(ds_slack_value)" )

    write(file_handle,"\n" )
    
    write(file_handle,
          get_df2tex(
              df_gens) )

    write(file_handle,"\n" )
    
    write(file_handle,
          get_df2tex(
              df_non_gens) )

    write(file_handle,"\n" )
    
    write(file_handle,
          get_df2tex(
              df_all_nodes_result) )        
    
end

distributed_slack_pf_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "distributed_slack_pf.csv")


CSV.write(distributed_slack_pf_csv_filename,
          df_all_nodes_result )

#---------------------------------------------------
# Sudden load change line outage pertubation
#---------------------------------------------------

sim_result_sudden_load_change_line_outage_pertubation =
    get_sim_mm_sudden_load_or_line_outage_pertubation(

    pertubation_factor,
    restoration_factor,

    pertubation_time,
    restoration_time;
        
    Δt1 = 1.5,
    Δt2 = 1.5,
    
    net_data_by_components_file,

    target_bus_name = "bus4",
    # pertubation_factor = 1.10,
    # restoration_factor = 1.0,
    
    sim_type,

    outage_type =
        :line_outage,
    
    on_fault_time,
    clear_fault_time,
    
    line_outage_time,
    generation_adjustment_time,
    
    components_libs_dir,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit  = [1],

    list_edges_to_have_fault   = [ 8 ],
    clear_fault_selection_list = [1],
        
    list_network_status = 
            [:pre_fault_state,
             :post_fault_state],
    
    with_faults =
        false,    
    
    # # base setting and some booleans
    
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    
    timespan      = 20,    
    dt            = 0.01,
    Δt            = 1.0 / 2^(4),

    # # digits rounding
    
    fractional_digits = 6 )

#---------------------------------------------------

sim_result_sudden_load_change_line_outage_pertubation =
    get_sim_mm_sudden_load_or_line_outage_pertubation(
        ;case_name,
        json_net_data_by_components_file =
            json_net_data_by_components_file,

        pertubation_factor,
        restoration_factor,
        pertubation_time,
        restoration_time,

        sim_type,
        script_dir,

        target_bus_name = "bus4",

        

        outage_type =
            :line_outage,

        Δt1 = 1.5,
        Δt2 = 1.5,

        on_fault_time,
        clear_fault_time,

        line_outage_time,
        generation_adjustment_time,

        components_libs_dir,

        list_fault_point_from_node_a = [0.3],
        list_fault_resistance = [0.001],
        list_no_line_circuit  = [1],

        list_edges_to_have_fault   = [ 8 ],
        clear_fault_selection_list = [1],

        list_network_status = 
                [:pre_fault_state,
                 :post_fault_state],

        with_faults =
            false,    

        # # base setting and some booleans

        basekV = 1.0,
        use_pu_in_PQ = true,
        line_data_in_pu = true,

        use_init_u0 = false,
        use_nlsolve = false,

        # # solvers and settings

        pf_alg        = NewtonRaphson(),    
        ode_alg       = Rodas4(),
        dae_alg       = IDA(),

        abstol        = 1e-12,
        reltol        = 1e-12,

        # # Simulation timespan

        timespan      = 20,    
        dt            = 0.01,
        Δt            = 1.0 / 2^(4),

        # # digits rounding

        fractional_digits = 6 )

#---------------------------------------------------
# Sudden load change
#---------------------------------------------------

sim_results_sudden_load_change =
    sim_mm_sudden_load_change(
        1.10,  # pertubation_factor
        1.0,   # restoration_factor
        10.0,  # pertubation_time,
        10.5;  # restoration_time
    
        Δt1 = 1.5,
        Δt2 = 1.5,
        
        net_data_by_components_file,

        # target_bus_name = "bus5",
        target_bus_name = "bus5",
        # pertubation_factor = 1.10,
        # restoration_factor = 1.0,

        sim_type,

        components_libs_dir )


rtds_plot_sudden_load_change =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            sim_results_sudden_load_change,
            (:system_sol,
             :model_syms,
             :gens_nodes_names,
             :SM_gens_nodes_names,
             :non_gens_nodes_names,
             :sim_timespan ) )... )

rtds_plot_sudden_target_parameter_plot =
    getproperty(
        sim_results_sudden_load_change,
        :target_parameter_plot)


δ_a_plot =
    getproperty(
        rtds_plot_sudden_load_change,
        :δ_a_plot)


ω_a_plot =
    getproperty(
        rtds_plot_sudden_load_change,
        :ω_a_plot)


vf_tilade_a_plot =
    getproperty(
        rtds_plot_sudden_load_change,
        :vf_tilade_a_plot)

xg2_a_plot =
    getproperty(
        rtds_plot_sudden_load_change,
        :xg2_a_plot)

plot_gens_vh =
    getproperty(
        rtds_plot_sudden_load_change,
        :plot_gens_vh)

plot_gens_θh =
    getproperty(
        rtds_plot_sudden_load_change,
        :plot_gens_θh)

plot_non_gens_vh =
    getproperty(
        rtds_plot_sudden_load_change,
        :plot_non_gens_vh)

plot_non_gens_θh =
    getproperty(
        rtds_plot_sudden_load_change,
        :plot_non_gens_θh)

#---------------------------------------------------

sim_sudden_load_change(
    net_data_by_components_file;

    # target_bus_name = "bus5",
    target_bus_name = "bus4",
    pertubation_factor = 1.10,
    restoration_factor = 1.0,

    pertubation_time = 10.0,
    restoration_time = 10.5,
    
    Δt1 = 1.5,
    Δt2 = 1.5,
    
    
    sim_type,

    case_name,
    components_libs_dir,
    data_dir,
    json_case_dir,
    
    results_dir,
    figure_dir,
    sd_dynamics_sim_csv_filename,

    # # base setting and some booleans
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    timespan      = 20,    
    dt            = 0.001,
    Δt            = 1.0 / 2^(4),

    # # digits rounding
    fractional_digits = 6 )


#---------------------------------------------------
## system fault states
#---------------------------------------------------

"""

# Possible system_status

system_status = :pre_fault_state,
system_status = :fault_state
system_status = :post_fault_state

"""

#---------------------------------------------------
## Test cases ode mass matrix and dae
#---------------------------------------------------

# case_name = "case9"

case_name = "case9"

# json_net_data_by_components_file =
#     "net-static-data-avr-sauer-gov-sauer.json"

json_net_data_by_components_file =
    "net-static-data-avr-rtds-gov-ieee-tgov2sauer.json"
    # "net-static-data-avr-rtds-gov-ieee-tgov1.json"

# "net-static-data-avr-rscad-gov_ieee_tgov1.json"


json_case_dir =
    joinpath(data_dir
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
#---------------------------------------------------

test_status_steady_state_data =
    get_a_status_steady_state_data(
        :pre_fault_state;
        with_faults = false,
        net_data_by_components_file,
        components_libs_dir,        

        timespan,

        on_fault_time,
        clear_fault_time,

        list_fault_point_from_node_a,
        list_fault_resistance,
        list_no_line_circuit,

        list_edges_to_have_fault,
        clear_fault_selection_list,

        basekV,    
        use_pu_in_PQ,
        line_data_in_pu )

#---------------------------------------------------

# sauer_test_Pg_Qg_Png_Qng_Pll_Qll_Idx

(t1_dyn_P_gens_Idxs,
 t1_dyn_Q_gens_Idxs,
 t1_dyn_P_non_gens_Idxs,
 t1_dyn_Q_non_gens_Idxs,
 t1_dyn_P_gens_loc_load_Idxs,
 t1_dyn_Q_gens_loc_load_Idxs) =
     NamedTupleTools.select(
        getproperty(
            getproperty(
                test_status_steady_state_data,
            :static_prefault_paras),
        :Pg_Qg_Png_Qng_Pll_Qll_Idx),
     (:dyn_P_gens_Idxs, :dyn_Q_gens_Idxs,
      :dyn_P_non_gens_Idxs,
      :dyn_Q_non_gens_Idxs,
      :dyn_P_gens_loc_load_Idxs,
      :dyn_Q_gens_loc_load_Idxs) )


# sauer_test_dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx =

(t2_dyn_ω_ref_Idx,
 t2_dyn_v_ref_Idx,
 t2_dyn_p_order_Idx,
 t2_dyn_Png_Idx,
 t2_dyn_Qng_Idx,
 t2_dyn_Pll_Idx,
 t2_dyn_Qll_Idx) =
     NamedTupleTools.select(
         getproperty(
             getproperty(
                 test_status_steady_state_data,
                 :static_prefault_paras),
             :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx),
         (:dyn_ω_ref_Idx,
          :dyn_v_ref_Idx,
          :dyn_p_order_Idx,
          :dyn_Png_Idx,
          :dyn_Qng_Idx,
          :dyn_Pll_Idx,
          :dyn_Qll_Idx) )


# sauer_test_dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx

(t3_dyn_δ_Idx,
 t3_dyn_ed_dash_Idx,
 t3_dyn_eq_dash_Idx,
 t3_dyn_Png_Idx,
 t3_dyn_Qng_Idx,
 t3_dyn_Pll_Idx, dyn_Qll_Idx) =
     NamedTupleTools.select(
         getproperty(
             getproperty(
                 test_status_steady_state_data,
                 :static_prefault_paras),
             :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx),
         (:dyn_δ_Idx,
          :dyn_ed_dash_Idx,
          :dyn_eq_dash_Idx,
          :dyn_Png_Idx,
          :dyn_Qng_Idx,
          :dyn_Pll_Idx,
          :dyn_Qll_Idx) )

#---------------------------------------------------

t1_ode_gens_para = getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :ode_gens_para)


(t1_H,
 t1_D,
 t1_X_d,
 t1_X_q,
 t1_X_d_dash,
 t1_X_q_dash,
 t1_T_d_dash,
 t1_T_q_dash,
 t1_Sn) =
     NamedTupleTools.select(
         t1_ode_gens_para,
         (:H,
          :D,
          :X_d,
          :X_q,
          :X_d_dash,
          :X_q_dash,
          :T_d_dash,
          :T_q_dash,
          :Sn) )

t1_sauer_test_pf_P_gens =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :pf_P_gens)


t1_sauer_test_pf_Q_gens =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :pf_Q_gens)


t1_p_order =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :p_order)

t1_gens_vh =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_vh)


t1_gens_θh =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_θh)


t1_gens_id =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_id)


t1_gens_iq =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_iq)


t1_gens_δ =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_δ)


t1_gens_ed_dash =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_ed_dash)


t1_gens_eq_dash =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_eq_dash)


t1_sauer_test_vh  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :vh)

t1_sauer_test_θh  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :θh)


t1_sauer_test_Pg_Qg_Png_Qng_Pll_Qll  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
         :Pg_Qg_Png_Qng_Pll_Qll)


t1_sauer_test_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
         :δ_ed_dash_eq_dash_Png_Qng_Pll_Qll)


t1_sauer_test_ωref0_vref0_porder0_id_iq_vh  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
         :ωref0_vref0_porder0_id_iq_vh)

t1_sauer_test_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
         :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll)

get_dynamic_model_Pq_Qg(
    t1_gens_δ,
    t1_gens_ed_dash,
    t1_gens_eq_dash,
    t1_gens_vh,
    t1_gens_θh,
    zeros(length(t1_gens_θh)),
    t1_X_d_dash,
    t1_X_q_dash,
    t1_gens_id,
    t1_gens_iq )

#---------------------------------------------------

t2_sauer_test_pf_P_gens =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :pf_P_gens)

t2_sauer_test_pf_Q_gens =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :pf_Q_gens)


t2_p_order =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :p_order)


t2_gens_vh =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_vh)


t2_gens_θh =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_θh)


t2_gens_id =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_id)


t2_gens_iq =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_iq)


t2_gens_δ =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_δ)


t2_gens_ed_dash =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_ed_dash)


t2_gens_eq_dash =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :gens_eq_dash)



t2_sauer_test_vh  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :vh)

t2_sauer_test_θh  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
        :θh)

t2_sauer_test_Pg_Qg_Png_Qng_Pll_Qll  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
         :Pg_Qg_Png_Qng_Pll_Qll)

t2_sauer_test_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
         :δ_ed_dash_eq_dash_Png_Qng_Pll_Qll)

t2_sauer_test_ωref0_vref0_porder0_id_iq_vh  =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
         :ωref0_vref0_porder0_id_iq_vh)

t2_sauer_test_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll =
    getproperty(
        getproperty(
            test_status_steady_state_data,
            :static_prefault_paras),
         :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll)

get_dynamic_model_Pq_Qg(
    t2_gens_δ,
    t2_gens_ed_dash,
    t2_gens_eq_dash,
    t2_gens_vh,
    t2_gens_θh,
    zeros(length(t2_gens_θh)),
    t1_X_d_dash,
    t1_X_q_dash,
    t2_gens_id,
    t2_gens_iq )


system_simulation_parameters =
    get_system_simulation_parameters(
        net_data_by_components_file #, nothing
        ;
        components_libs_dir =
            components_libs_dir,
        basekV = 1.0,    
        use_pu_in_PQ = true,
        line_data_in_pu = true,

        use_init_u0 = false,

        use_nlsolve = false,

        pf_alg        = NewtonRaphson(),

        abstol        = 1e-12,

        reltol        = 1e-12 )


(t_pf_P_gens,
 t_pf_Q_gens) =
     NamedTupleTools.select(
         system_simulation_parameters,
         (:pf_P_gens,
            :pf_Q_gens))

#---------------------------------------------------
#---------------------------------------------------

# "gov": {
#     "p_max": 2.2,
#     "R": 0.02,
#     "p_min": 0.1,
#     "Ts": 0.2,
#     "Tc": 0.4
# },


"""
sim_model_by_mass_matrix_ode_by_funcs_dynamics!
sim_model_by_mass_matrix_by_ode_pf_funcs!

not working
sim_model_by_mass_matrix_ode_by_model_dynamics!

"""
rtds_plot_sim_model_mass_matrix_ode =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;sim_model_by_mass_matrix_ode_by_funcs_dynamics!(
            :pre_fault_state ;
            with_faults = false,
            
            net_data_by_components_file,
            components_libs_dir,
            
            timespan=100,
            on_fault_time,
            clear_fault_time,

            list_fault_point_from_node_a,
            list_fault_resistance,
            list_no_line_circuit,

            list_edges_to_have_fault,
            clear_fault_selection_list,


            ode_alg = Rodas4()
        )...)

mm_δ_a_plot =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :δ_a_plot)


mm_ω_a_plot =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :ω_a_plot)

mm_vf_tilade_a_plot =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :vf_tilade_a_plot)

mm_xg2_a_plot =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :xg2_a_plot)

mm_plot_gens_vh =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :plot_gens_vh)

mm_plot_gens_θh =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :plot_gens_θh)

mm_plot_non_gens_vh =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :plot_non_gens_vh)

mm_plot_non_gens_θh =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :plot_non_gens_θh)

#---------------------------------------------------
#---------------------------------------------------

# sim_model_by_dae_pf_funcs!

# sim_model_by_dae_funcs_dynamics!

# sim_model_by_dae_model_dynamics!

rtds_plot_sim_model_dae =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;sim_model_by_dae_funcs_dynamics!(
            :pre_fault_state ;
            with_faults = false,
            
            net_data_by_components_file,
            components_libs_dir,

            timespan=10,
            on_fault_time,
            clear_fault_time,

            list_fault_point_from_node_a,
            list_fault_resistance,
            list_no_line_circuit,

            list_edges_to_have_fault,
            clear_fault_selection_list )...)

δ_a_plot =
    getproperty(
        rtds_plot_sim_model_dae,
        :δ_a_plot)


ω_a_plot =
    getproperty(
        rtds_plot_sim_model_dae,
        :ω_a_plot)


vf_tilade_a_plot =
    getproperty(
        rtds_plot_sim_model_dae,
        :vf_tilade_a_plot)

xg2_a_plot =
    getproperty(
        rtds_plot_sim_model_dae,
        :xg2_a_plot)

plot_gens_vh =
    getproperty(
        rtds_plot_sim_model_dae,
        :plot_gens_vh)

plot_gens_θh =
    getproperty(
        rtds_plot_sim_model_dae,
        :plot_gens_θh)

plot_non_gens_vh =
    getproperty(
        rtds_plot_sim_model_dae,
        :plot_non_gens_vh)

plot_non_gens_θh =
    getproperty(
        rtds_plot_sim_model_dae,
        :plot_non_gens_θh)

#---------------------------------------------------
# line outage pertubation
#---------------------------------------------------

# outage_type = :line_outage,    
# :line_outage_wt_pref_adjs
# :line_outage_wt_vpref_adjs

rtds_plot_sim_line_outage_pertubation =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            sim_line_outage_pertubation_by_dae(
                :line_outage, # outage_type

                on_fault_time,
                clear_fault_time,

                line_outage_time,
                generation_adjustment_time;

                net_data_by_components_file,

                timespan = 20,

                components_libs_dir,

                list_fault_point_from_node_a = [0.3],
                list_fault_resistance = [0.001],
                list_no_line_circuit  = [1],

                list_edges_to_have_fault  = [ 8 ],
                clear_fault_selection_list = [1] ),

            (:system_sol,
            :model_syms,
            :gens_nodes_names,
            :SM_gens_nodes_names,
            :non_gens_nodes_names,:sim_timespan))...)


δ_a_plot =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :δ_a_plot)


ω_a_plot =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :ω_a_plot)


vf_tilade_a_plot =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :vf_tilade_a_plot)

xg2_a_plot =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :xg2_a_plot)

plot_gens_vh =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :plot_gens_vh)

plot_gens_θh =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :plot_gens_θh)

plot_non_gens_vh =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :plot_non_gens_vh)

plot_non_gens_θh =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :plot_non_gens_θh)

#---------------------------------------------------

rtds_plot_sim_line_outage_pertubation =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
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
                clear_fault_selection_list = [1] ),

            (:system_sol,
            :model_syms,
            :gens_nodes_names,
            :SM_gens_nodes_names,
            :non_gens_nodes_names,:sim_timespan))...)


δ_a_plot =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :δ_a_plot)


ω_a_plot =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :ω_a_plot)


vf_tilade_a_plot =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :vf_tilade_a_plot)

xg2_a_plot =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :xg2_a_plot)

plot_gens_vh =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :plot_gens_vh)

plot_gens_θh =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :plot_gens_θh)

plot_non_gens_vh =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :plot_non_gens_vh)

plot_non_gens_θh =
    getproperty(
        rtds_plot_sim_line_outage_pertubation,
        :plot_non_gens_θh)


#---------------------------------------------------
#---------------------------------------------------


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


δ_a_plot =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :δ_a_plot)


ω_a_plot =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :ω_a_plot)

plot_gens_vh =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :plot_gens_vh)

plot_gens_θh =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :plot_gens_θh)

plot_non_gens_vh =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :plot_non_gens_vh)

plot_non_gens_θh =
    getproperty(
        rtds_plot_sim_model_mass_matrix_ode,
        :plot_non_gens_θh)

#---------------------------------------------------
## ntuple_status_steady_state_data
#---------------------------------------------------

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
            [:pre_fault_state,
             :post_fault_state] )


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

# propertynames(sta_pf_red_sol)

# (:gens_current_injection, :gens_loc_load_current, :gens_nodes_network_current, :non_gens_nodes_network_current, :nodes_network_current, :S_gens, :pf_P_gens, :pf_Q_gens, :gens_S, :pf_P_g_gens, :pf_Q_g_gens, :Igen, :Inet, :Iinj, :pu_Igen, :Sbus_n, :GenSinj, :If, :It, :Ibranches, :vh, :θh, :θh_deg, :Vbus, :gens_E, :gens_mag_E, :gens_ang_E, :gens_vh, :gens_θh, :gens_nodes_idx)


# propertynames(dyn_pf_fun_kwd_net_idxs)

# (:slack_gens_nodes_idx, :non_slack_gens_nodes_idx, :gens_nodes_idx, :non_gens_nodes_idx, :gens_nodes_with_loc_loads_idx, :all_nodes_idx)

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

# IEEE 9

# 3×4 Matrix{Float64}:

#  vh     θh      P_gens   Q_gens

#  1.04   0.0     0.7164   0.2705
#  1.025  9.28    1.63     0.0665
#  1.025  4.6648  0.85    -0.1086


# IEEE 14

# 5×4 Matrix{Float64}:

#  vh       θh       P_gens   Q_gens

#  1.06     0.0      2.3239  -0.1655
#  1.045   -4.9826   0.4      0.4356
#  1.01   -12.7251   0.0      0.2508
#  1.07   -14.2209  -0.0      0.1273
#  1.09   -13.3596  -0.0      0.1762

#---------------------------------------------------

# Non gens

pf_non_gens_results =
    hcat(non_gens_vh,
         non_gens_θh_deg )

# IEEE 9

# 6×2 Matrix{Float64}:

#   vh     θh_deg

#  1.0258  -2.2168
#  1.0127  -3.6874
#  1.0324   1.9667
#  1.0159   0.7275
#  1.0258   3.7197
#  0.9956  -3.9888

# IEEE 14

# 9×2 Matrix{Float64}:

#   vh     θh_deg

#  1.0177  -10.3129
#  1.0195   -8.7739
#  1.0615  -13.3596
#  1.0559  -14.9385
#  1.051   -15.0973
#  1.0569  -14.7906
#  1.0552  -15.0756
#  1.0504  -15.1563
#  1.0355  -16.0336

#---------------------------------------------------
#---------------------------------------------------
# Sudden change in load
#---------------------------------------------------
#---------------------------------------------------

#---------------------------------------------------
## simulation case
#---------------------------------------------------

# sim_type  = "simulation-validation"

sim_type = "validation-sudden-load-pertubation"

case_name = "case14"

# case_name = "case9"

json_net_data_by_components_file =
    "net-static-data-avr-rtds-gov-ieee-tgov1.json"

# "net-static-data-avr-sauer-gov-sauer.json"
# "net-static-data-avr-rscad-gov_ieee_tgov1.json"

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
             "$(case_name)-" *
                 "$(sim_type)-" *
                 "validation.tex")

#---------------------------------------------------

sd_dynamics_sim_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "$(sim_type)-" *
            "sim-sd-dynamics.csv")

# sauer_plot_sudden_load_change

#---------------------------------------------------

sim_results_sudden_load_change =
    sim_sudden_load_change(
        net_data_by_components_file;

        target_bus_name = "bus4",
        pertubation_factor = 1.10,
        restoration_factor = 1.0,

        sim_type,

        case_name,
        components_libs_dir,
        data_dir,
        json_case_dir,

        results_dir,
        figure_dir,
        sd_dynamics_sim_csv_filename )


rtds_plot_sudden_load_change =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            sim_results_sudden_load_change,
            (:system_sol,
             :model_syms,
             :gens_nodes_names,
             :SM_gens_nodes_names,
             :non_gens_nodes_names,
             :sim_timespan ) )... )

rtds_plot_sudden_target_parameter_plot =
    getproperty(
        sim_results_sudden_load_change,
        :target_parameter_plot)

# rtds_plot_sudden_load_change =
#     sim_sudden_load_change(
#         ;case_name,
#         json_net_data_by_components_file,
#         sim_type,
#         timespan)

# rtds_plot_sudden_load_change.δ_a_plot


# (:δ_a_plot, :ω_a_plot,
#  :ed_dash_a_plot, :eq_dash_a_plot,
#  :vr1_a_plot, :vr2_a_plot,
#  :vf_tilade_a_plot, :xg1_a_plot,
#  :xg2_a_plot, :plot_gens_vh,
#  :plot_gens_θh, :plot_non_gens_vh,
#  :plot_non_gens_θh )    


#---------------------------------------------------
#---------------------------------------------------
# Network line loss pertubation
#---------------------------------------------------
#---------------------------------------------------

#---------------------------------------------------
## simulation case
#---------------------------------------------------

sim_type = "validation-line-loss-pertubation"

case_name = "case14"

# case_name = "case9"

json_net_data_by_components_file =
    "net-static-data-avr-rtds-gov-ieee-tgov1.json"

# "net-static-data-avr-sauer-gov-sauer.json"
# "net-static-data-avr-rscad-gov_ieee_tgov1.json"

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
             "$(case_name)-" *
                 "$(sim_type)-" *
                 "validation.tex")

#---------------------------------------------------

sd_dynamics_sim_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "$(sim_type)-" *
            "sim-sd-dynamics.csv")


#---------------------------------------------------
# Simulation Period
#---------------------------------------------------

timespan      = 20.0

time_start    = 0.0

time_final    = timespan

dt            = 0.0001

Δt            = 1.0 / 2^(4)

tspan         = (0.0, timespan)

sim_timespan  = (0.0, timespan)

plot_timespan = (0.0, timespan)

#---------------------------------------------------
# base setting and some booleans 
#---------------------------------------------------

basekV          = 1.0

use_pu_in_PQ    = true

line_data_in_pu = true

#---------------------------------------------------
## solvers and settings
#---------------------------------------------------

use_init_u0 = false

use_nlsolve = false

pf_alg      = NewtonRaphson()

#---------------------------------------------------

ode_alg       = Rodas4()

# ode_alg       = ImplicitMidpoint()

dae_alg       = IDA()

abstol        = 1e-12

reltol        = 1e-12

#---------------------------------------------------
## fault time and data
#---------------------------------------------------

with_faults                   = false

Δt_generation_adjustment_time = 0.2

Δt_clear_time                 = 0.04

on_fault_time                 = 10.0

clear_fault_time =
    on_fault_time +
    Δt_clear_time

line_outage_time =
    on_fault_time

generation_adjustment_time =
    line_outage_time +
    Δt_generation_adjustment_time

#---------------------------------------------------

list_fault_point_from_node_a = [0.01]

list_fault_resistance        = [0.001]

list_no_line_circuit         =  [1]

list_edges_to_have_fault     = [ 4 ]

clear_fault_selection_list   = [1]

#---------------------------------------------------
# pertubation line loss outage wt_or_no_ref_adjs
#---------------------------------------------------

results_sim_line_loss_pertubation =
    sim_line_loss_pertubation(
        net_data_by_components_file;
        components_libs_dir,
        with_faults,
        sim_type,

        case_name,
        components_libs_dir,
        data_dir,

        # # fault or pertubation
        on_fault_time,
        clear_fault_time,
        line_outage_time,
        generation_adjustment_time,

        list_fault_point_from_node_a = [0.01],
        list_fault_resistance      = [0.001],
        list_no_line_circuit       =  [1],
        list_edges_to_have_fault   = [ 4 ],
        clear_fault_selection_list = [1],

        # # solvers and settings
        pf_alg        = NewtonRaphson(),    
        ode_alg       = Rodas4(),
        dae_alg       = IDA(),
        abstol        = 1e-12,
        reltol        = 1e-12,

        # # Simulation timespan
        timespan      = 20,    
        dt            = 0.01,
        Δt            = 1.0 / 2^(4),

        # # base setting and some booleans
        basekV          = 1.0,
        use_pu_in_PQ    = true,
        line_data_in_pu = true,

        use_init_u0 = false,
        use_nlsolve = false,

        # # digits rounding
        fractional_digits   = 6,
        list_network_status = 
                    [:pre_fault_state,
                     :post_fault_state],

        # results folders
        results_dir = nothing,
        figure_dir  = nothing,
        sd_dynamics_sim_csv_filename = nothing )

#---------------------------------------------------
#---------------------------------------------------

list_network_status = 
            [:pre_fault_state,
             :post_fault_state] 

outage_type = :line_outage

t_line_loss_outage =
    get_line_loss_outage_wt_or_no_ref_adjs(
        outage_type,

        on_fault_time,
        clear_fault_time,

        line_outage_time,
        generation_adjustment_time,

        net_data_by_components_file;

        timespan,

        components_libs_dir,
        data_dir,

        list_fault_point_from_node_a,
        list_fault_resistance,
        list_no_line_circuit,

        list_edges_to_have_fault,
        clear_fault_selection_list,

        list_network_status =
            list_network_status )

plot_outage_type =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            get_line_loss_outage_wt_or_no_ref_adjs(
                outage_type,

                on_fault_time,
                clear_fault_time,
                
                line_outage_time,
                generation_adjustment_time,
                
                net_data_by_components_file;

                timespan,
                
                components_libs_dir,
                data_dir,

                list_fault_point_from_node_a,
                list_fault_resistance,
                list_no_line_circuit,

                list_edges_to_have_fault,
                clear_fault_selection_list,

                list_network_status =
                    list_network_status ) ,
            (:system_sol,
             :model_syms,
             :gens_nodes_names,
             :SM_gens_nodes_names,
             :non_gens_nodes_names,
             :sim_timespan) )... )


plot_outage_type =
    get_guick_single_vars_plots_dae_or_ode_sol(
        ;NamedTupleTools.select(
            get_line_loss_outage_wt_or_no_ref_adjs(
                outage_type,
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
             :sim_timespan))... )

#---------------------------------------------------
#---------------------------------------------------

list_outage_type =
    [:line_outage,
     :line_outage_wt_pref_adjs,
     :line_outage_wt_vpref_adjs]

dict_outage_type_sol =
    Dict()

dict_outage_type_plot =
    Dict{Symbol, NTuple}()

for an_outage_type in list_outage_type

    dict_outage_type_plot[an_outage_type] =
        get_guick_single_vars_plots_dae_or_ode_sol(
            ;NamedTupleTools.select(
                get_line_loss_outage_wt_or_no_ref_adjs(
                    an_outage_type,
                    line_outage_time,
                    generation_adjustment_time;
                    case_name = case_name,
                    json_net_data_by_components_file =
                         json_net_data_by_components_file,
                    components_libs_dir =
                    components_libs_dir,
                    data_dir =
                        data_dir,    

                    with_faults =
                        with_faults,    

                    timespan   =
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
                    [:pre_fault_state,
                     :post_fault_state] ) ,
                (:system_sol,
                 :model_syms,
                 :gens_nodes_names,
                 :SM_gens_nodes_names,
                 :non_gens_nodes_names,
                 :sim_timespan) )... )
        
end

#---------------------------------------------------
#---------------------------------------------------

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
            [:pre_fault_state,
             :post_fault_state] )

result_Ynet_pertubation =
    sim_Ynet_pertubation(
        on_fault_time,
        clear_fault_time,    
        ntuple_status_steady_state_data)

#---------------------------------------------------
# save julia object to latex file
#---------------------------------------------------

pf_tuple_julia_object =
    (;pf_gens_results,
      pf_non_gens_results)

pf_names_julia_object =
    propertynames(pf_tuple_julia_object)

open(tex_filename, "a") do file_handle
    
    for (name_object, a_julia_object) in
        zip(pf_names_julia_object,
            pf_tuple_julia_object)
        
        write(file_handle,"\n $(String(name_object)) = ")
        
        write( file_handle, latexify(
            a_julia_object;
            fmt=FancyNumberFormatter()))
    end
    
end

#---------------------------------------------------
#---------------------------------------------------

#---------------------------------------------------
# case model 
#---------------------------------------------------

# fault_status =
#     (no_fault = 0,
#      on_fault = 1,
#      clear_fault = 2,
#      partial_clear_fault = 3)

# system_fault_status =
#     [ fault_status.no_fault]

#---------------------------------------------------

# list_network_status = list_system_status =
#     [:pre_fault_state,
#      :fault_state,
#      :post_fault_state]

#     list_network_status = 
#         [:pre_fault_state,
#          :fault_state,
#          :post_fault_state ]


# dict_status_steady_state_data =
#     Dict{Symbol, NamedTuple }()

# for a_system_status in list_network_status

#     dict_status_steady_state_data[a_system_status] =
#         get_a_status_steady_state_data(
#             a_system_status;
#             with_faults = true,
#             # case_name = "case9",
#             net_data_by_components_file =
#                 net_data_by_components_file,
#             components_libs_dir =
#                 components_libs_dir,

#             timespan      = timespan,
#             on_fault_time = on_fault_time,
#             clear_fault_time = clear_fault_time,    

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

#             basekV = basekV,    
#             use_pu_in_PQ =
#                 use_pu_in_PQ,
#             line_data_in_pu =
#                 line_data_in_pu)

# end

# status_steady_state_data =
#     NamedTupleTools.namedtuple(
#         dict_status_steady_state_data)



# timespan   = 50


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
