module ePowerSim

include("headers/sd-dynamics-imports.jl")

include("headers/sd-dynamics-header.jl")

# Write your package code here.

# ------------------------------------------------------
#  System Dynamics Data Type
# ------------------------------------------------------

export AbstractPowerSystemComponent

export SdGen, SdNonGen, SdGov, SdAvr, SdPss

export SdNonGenPlant, SdGenPlant, SdBranchElement



# ------------------------------------------------------
#  Init function
# ------------------------------------------------------

export plants_generic_model_init_func
export get_generic_system_dynamics_init_and_refs
export get_init_a_gen_full_model

export get_init_gens_full_model
export get_state_init_SC_generic_model

export a_SM_plant_generic_model_init_func
export a_SC_plant_generic_model_init_func

# ------------------------------------------------------
#  Algebraic functions 
# ------------------------------------------------------

export algebraic_generic_pf_ΔI_mismatch!
export algebraic_generic_pf_ΔI_mismatch_sol(
export pertubation_algebraic_generic_pf_ΔPQ_mismatch!

export algebraic_generic_pf_ΔPQ_mismatch!
export algebraic_generic_pf_ΔPQ_mismatch_sol


# ------------------------------------------------------
#  ODE functions 
# ------------------------------------------------------

# Partitioned-explicit (PE) method.

export ode_a_SM_plant_generic_model_func!
export ode_a_SC_plant_generic_model_func!

export ode_gens_plants_generic_model_func!

export ode_generic_system_model_by_funcs_dynamics!
export ode_generic_system_model_dynamics!
export ode_generic_system_dynamics_by_ode_pf_funcs!

export ode_gens_plants_generic_model_func!

export ode_a_gen_generic_model_func!


export ode_a_gen_generic_model_by_ext_idq_func!
export ode_gens_generic_model_by_ext_idq_func!

export ode_generic_model_func!
export ode_generic_model_dynamics!

export ode_gens_plants_generic_model_wt_cb_state_func!

export ode_test_a_SM_plant_generic_model_func!
export ode_cb_test_gens_plants_generic_model_func!

# ------------------------------------------------------
#  DAE functions 
# ------------------------------------------------------

# Simultaneous implicit DAE (SI) method

export dae_a_SM_plant_generic_model_func!
export dae_a_SC_plant_generic_model_func!
export dae_a_gen_generic_model_func!

export dae_a_gen_generic_model_by_ext_idq_func!
export dae_gens_generic_model_by_ext_idq_func!

export dae_generic_model_func!
export dae_generic_model_dynamics!
export dae_gens_plants_generic_model_func!

# ------------------------------------------------------
#  Systems ODE functions by mass matrix
# ------------------------------------------------------

export ode_a_SC_plant_generic_model_func!
export ode_a_SM_plant_generic_model_func!
export ode_gens_plants_generic_model_func!

export mm_ode_generic_system_model_by_funcs_dynamics!
export mm_ode_generic_system_dynamics_by_ode_pf_funcs!

# ------------------------------------------------------
#  Systems DAE functions 
# ------------------------------------------------------

export dae_a_SC_plant_generic_model_func!
export dae_a_SM_plant_generic_model_func!
export dae_gens_plants_generic_model_func!

export dae_generic_system_dynamics_by_dae_pf_funcs!
export dae_generic_system_model_by_funcs_dynamics!
export dae_generic_system_model_dynamics!

export get_generic_namedtuple_per_plant_para_wt_kwd_para

# ------------------------------------------------------
#  Systems Dynamics functions 
# ------------------------------------------------------

export mm_ode_generic_system_dynamics_by_ode_pf_funcs!
export mm_ode_generic_system_model_by_funcs_dynamics!
export mm_ode_generic_system_model_dynamics!

export mm_line_outage_generic_dynamics_wt_pre_post_fault_by_ode_pf_funcs!

export line_outage_generic_dynamics_wt_pre_post_fault_by_dae_pf_funcs!

export dae_generic_system_dynamics_by_dae_pf_funcs!
export dae_generic_system_model_by_funcs_dynamics!
export dae_generic_system_model_dynamics!

export Ynet_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!

export get_generic_network_fault_pertubation
export get_net_status_Yint_and_Yred_matrices

#-------------------------------------------------------
# callbacks
#-------------------------------------------------------

export on_fault_condition
export clear_fault_condition

export on_line_outage_condition
export on_line_outage_affect!

export on_generation_adjustment_condition
export on_generation_adjustment_affect!

export on_fault_Ynet_affect!
export clear_fault_Ynet_affect!

#-------------------------------------------------------
# simulations
#-------------------------------------------------------

export plot_benchmark_results
export create_benchmarkgroups_admittance_vector_funcs
export create_benchmarkgroups_pf_mismatch

export sim_full_model_distributed_slack_pf
export sim_red_model_distributed_slack_pf

export sim_model_by_mass_matrix_ode_by_model_dynamics!
export sim_model_by_mass_matrix_by_ode_pf_funcs!

export sim_model_by_dae_funcs_dynamics!
export sim_model_by_dae_model_dynamics!

export sim_sudden_load_change
export sim_mm_sudden_load_change

export get_sim_mm_sudden_load_or_line_outage_pertubation

# export sim_network_Ynet_pertubation

export sim_line_loss_pertubation

export sim_line_outage_pertubation_by_dae
export sim_line_outage_pertubation_by_mm_ode

export sim_line_outage_pertubation

export sim_Ynet_pertubation


#-------------------------------------------------------
# data
#-------------------------------------------------------

export get_system_net_static_data
export get_system_net_data_wt_static_parameters

export get_a_status_steady_state_data
export get_system_simulation_parameters
export get_status_steady_state_parameters
export get_ntuple_status_steady_state_data
export get_generic_system_simulation_parameters
export get_system_simulation_parameters_wt_faults

export get_net_data_by_components_from_json_file
export get_edges_nt_fbus_tbus_by_json
export get_edges_orientation_by_generic
export get_Cnb_by_orientations
export get_nodes_incident_edges_by_orientations

export get_net_nodes_type_idxs_by_json
export get_dict_net_streamlined_idx_by_nodes_type_idxs

export get_components_properties_by_json
export get_selected_vec_nt_to_vec_vec
export get_selected_comps_ode_para_by_json
export get_ode_gens_generic_para

export make_lines_faults_data_set
export get_cleared_selected_lines_faults_data_set

export get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend
export get_Pg_inj_Qg_inj_Png_Qng_Pll_Qll
export get_dynamic_model_Pq_Qg

export get_loss_particpation_by_gens_rating
export get_loss_particpation_by_gens_loading

export get_gens_active_power_particpation_by_rating
export get_gens_reactive_power_particpation_by_rating

export get_gens_active_power_particpation_by_loading
export get_gens_reactive_power_particpation_by_loading

export get_sum_line_losses
export get_losses_per_line
export get_lines_total_losses

export get_total_P_network_loss_by_flattened_Ynet
export get_total_Q_network_loss_by_flattened_Ynet

export get_losses_per_node_by_flattened_Ynet

export pertubation_by_itegrator

export get_sol_auxilliary_results
export make_plot_gens_streamedlined_auxilliary_results
export save_sol_auxilliary_results_plot

export round_up_Ynet
export write_vector_or_matrix_to_tex

export get_make_df_header_generic_model_dynamics_para

export get_P_or_Q_idx_in_generic_model_dynamics_para
export get_a_node_state_algb_vars_indices_in_system

export get_line_loss_outage_wt_or_no_ref_adjs


#-------------------------------------------------------
# opf
#-------------------------------------------------------

export opf
export opf_by_relaxation
export get_opf_by_relaxation_by_scenario

export get_opf_wt_generic_system_simulation_parameters
export get_opf_by_relaxation_wt_generic_system_simulation_parameters
export get_opf_net_optimisation_parameters

export opf_by_line_loss

export get_economic_dispatch_by_scenario
export get_unit_commitment_by_scenario
export get_mixed_complementarity


#-------------------------------------------------------
# benchmarks
#-------------------------------------------------------

export create_benchmarkgroups_admittance_vector_funcs
export create_benchmarkgroups_pf_mismatch
export plot_benchmark_results

#-------------------------------------------------------
# results
#-------------------------------------------------------

export get_guick_single_vars_plots_dae_or_ode_sol
export get_guick_group_vars_plots_dae_or_ode_sol
export save_network_pertubation_sim_plot

export get_generic_results_ds_pf_sol_u
export get_generic_results_conti_or_ds_pf_sol_u

export get_generic_results_ds_pf_red_sol_u
export get_generic_results_conti_or_ds_pf_red_sol_u

#-------------------------------------------------------
# static mismatch
#-------------------------------------------------------

export disaggregate_sta_pf_keywords_parameter

export get_ΔPQ_mismatch_by_Ybus
export get_ΔI_mismatch_by_Ybus
export get_ΔPQ_mismatch_by_sparse_Ybus

export get_ΔPQ_mismatch_by_Ynet
export get_ΔI_mismatch_by_Ynet

export get_ΔPQ_mismatch_by_S_inj_Yπ_net

export get_ΔI_mismatch_by_Yπ_net
export get_ΔPQ_mismatch_by_Yπ_net

export get_generic_sta_pf_ΔPQ_mismatch
export get_generic_sta_pf_ΔI_mismatch

export get_a_model_integrated_pf_sta_ΔI_mismatch_generic
export get_a_model_integrated_pf_sta_ΔPQ_mismatch_generic
export get_a_model_integrated_pf_sta_ΔPQ_mismatch

export get_results_static_pf_red_sol_u
export get_generic_results_pf_sta_red_sol_u

export get_pf_sta_ΔPQ_mismatch_sol_by_generic

export get_model_distributed_slack_pf_ΔPQ_mismatch!
export get_red_model_distributed_slack_pf_ΔPQ_mismatch!

#-------------------------------------------------------
# Jacobian
#-------------------------------------------------------

export sta_pf_Jac!
export sta_pf_Jac!_df_dx_by_Ynet_or_Yπ_net!
export sta_pf_Jac!_df_dp_by_Ynet_or_Yπ_net!

# get_ph_pk_qh_qk_by_Ybr_cal
# get_nodes_∑_ynj_x_vj
# get_nodes_∑_ynj_x_vj_by_Yπ_net

#-------------------------------------------------------
# others
#-------------------------------------------------------

export get_Ynet_by_trimed_Shunt
export get_Ynet
export get_Yπ_net
export get_Ybus
export get_Ynet_sp_sh

export get_size_Ybus

export get_Ynet_wt_nodes_idx_wt_adjacent_nodes_diff
export get_size_Ynet_wt_nodes_idx_wt_adjacent_nodes
export get_size_Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes

#-------------------------------------------------------
# Experimental magnetic equations
#-------------------------------------------------------

# export common_equations
# export stator_equations_by_flux_dynamics

# export stator_equations_by_zero_flux_dynamics
# export stator_equations_by_small_ω_approximation

# export magnetic_equations_by_sauer_pai_model
# export algebraic_magnetic_equations_by_sauer_pai_model

# export magnetic_equations_by_marconato_model
# export algebraic_magnetic_equations_by_marconato_model

# export magnetic_equations_by_anderson_fouad_model
# export algebraic_magnetic_equations_by_anderson_fouad_model


end
