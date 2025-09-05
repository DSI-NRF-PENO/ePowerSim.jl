


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
## simulation case
#---------------------------------------------------

case_name = "case14"

# case_name = "case9"

timespan = 50.0

sim_type = "line-loss-pertubation"


fractional_digits = 4

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

json_net_data_by_components_file =
    "net-static-data-avr-sauer-gov-sauer.json"


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


sd_dynamics_sim_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "$(sim_type)-" *
            "sim-sd-dynamics.csv")

#---------------------------------------------------
## system fault states
#---------------------------------------------------

# """

# # Possible system_status

# system_status = :pre_fault_state,
# system_status = :fault_state
# system_status = :post_fault_state

# """

#---------------------------------------------------
## fault data
#---------------------------------------------------

Δt_clear_time = 0.02

Δt_generation_adjustment_time = 0.2

#---------------------------------------------------

on_fault_time = 10.0

clear_fault_time =
    on_fault_time + Δt_clear_time

#---------------------------------------------------

line_outage_time = 10.0

generation_adjustment_time =
    line_outage_time + Δt_generation_adjustment_time

#---------------------------------------------------

with_faults = false

list_fault_point_from_node_a = [0.01]
list_fault_resistance = [0.001]
list_no_line_circuit =  [1]

list_edges_to_have_fault = [ 8 ]
clear_fault_selection_list = [1]

#---------------------------------------------------
# base setting and some booleans 
#---------------------------------------------------

basekV = 1.0

use_pu_in_PQ = true

line_data_in_pu = true

#---------------------------------------------------
# Simulation Period
#---------------------------------------------------

time_start = 0.0

time_final = timespan

dt = 0.0001

Δt = 1.0 / 2^(4)

tspan =
    (0.0, timespan)

sim_timespan =
    (0.0, timespan)

plot_timespan =
    (0.0, timespan)

#---------------------------------------------------
## solvers and settings
#---------------------------------------------------

use_init_u0 = false

use_nlsolve = false

pf_alg = NewtonRaphson()

#---------------------------------------------------

ode_alg       = Rodas4()

# ode_alg      = ImplicitMidpoint()

dae_alg       = IDA()

abstol        = 1e-12

reltol        = 1e-12

#---------------------------------------------------
## functions
#---------------------------------------------------

#---------------------------------------------------

list_outage_type =
    [:line_outage,
     :line_outage_wt_pref_adjs,
     :line_outage_wt_vpref_adjs]


dict_outage_type_sol = Dict() # Dict{Symbol, NTuple}()

for an_outage_type in list_outage_type

    an_outage_type_sol =
        get_line_loss_outage_wt_or_no_ref_adjs(
            an_outage_type,
            line_outage_time,
            generation_adjustment_time;
            case_name = "case9",
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
                line_data_in_pu )

    dict_outage_type_sol[an_outage_type] =
        an_outage_type_sol

    save_network_pertubation_sim_plot(
        case_name;
        figure_dir,
        sim_type =
            "$(sim_type)",
        line_in_fault_name =
            "line-8",
        include_v_θ_plot =
            true,
        sub_folder =
            "$(an_outage_type)",
        NamedTupleTools.select(
            an_outage_type_sol,
            (:system_sol,
            :model_syms,
            :gens_nodes_names,
            :SM_gens_nodes_names,
            :non_gens_nodes_names,
            :sim_timespan))...)

    save_sol_auxilliary_results_plot(
        case_name;
         NamedTupleTools.select(
            an_outage_type_sol,
             (:system_sol,
              :state_labels,
              :algebraic_vars_labels,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs))...,
        figure_dir,
        sim_type =
            "$(sim_type)",
        line_in_fault_name =
            "line-8",
        sub_folder =
            "$(an_outage_type)" )

    sd_dynamics_sim_df =
        DataFrame(an_outage_type_sol.system_sol)

    sd_dynamics_sim_df[!, :] =
        round.(
            sd_dynamics_sim_df[:, :],
            digits=fractional_digits)

    sd_dynamics_sim_csv_filename =
            joinpath(results_dir,
                     "$(case_name)-" *
                         "$(sim_type)-" *
                         "$(an_outage_type).csv")

    CSV.write(sd_dynamics_sim_csv_filename,
              sd_dynamics_sim_df )
        
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
            line_data_in_pu )

#---------------------------------------------------

post_sta_PQ =
    getproperty(
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :pre_fault_state),
            :static_prefault_paras),
        :post_sta_PQ)

post_dyn_PQ =
    getproperty(
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :post_fault_state),
            :dynamic_status_paras),
        :post_dyn_PQ)

(;pf_P_gens,
 P_non_gens,
 Q_non_gens,
 P_g_loc_load,
 Q_g_loc_load) =
     NamedTupleTools.select(
         post_sta_PQ,
         (:pf_P_gens,
          :P_non_gens,
          :Q_non_gens,
          :P_g_loc_load,
          :Q_g_loc_load))

(;dyn_pf_P_gens,
 ) =
     NamedTupleTools.select(
         post_dyn_PQ,
         (:dyn_pf_P_gens,         
          ))

#---------------------------------------------------

Total_demands = sum(P_non_gens) + (
    length(P_g_loc_load) == 0 ? 0 : sum(
        P_g_loc_load))

Δ_genP =
    sum(dyn_pf_P_gens) - sum(pf_P_gens)

pre_fault_loss =
    sum(pf_P_gens) - Total_demands

post_fault_loss =
    sum(dyn_pf_P_gens) - Total_demands

Δ_loss =
    post_fault_loss - pre_fault_loss


"""

IEEE 9 bus

Total_demands = 3.15

Δ_genP = 0.018878627948037074

pre_fault_loss = 0.04641021474482132

post_fault_loss = 0.0652888426928584

Δ_loss = 0.018878627948037074

"""

#---------------------------------------------------

"""

line_loss_outage_wt_or_no_ref_adjs =
    get_line_loss_outage_wt_or_no_ref_adjs(
        :line_outage,
        ntuple_status_steady_state_data,
        line_outage_time,
        generation_adjustment_time; 
        sim_timespan = sim_timespan,    
        dae_alg = dae_alg,
        abstol  = abstol,
        reltol  = reltol )


auxilliary_results =
    get_sol_auxilliary_results(
        getproperty(
            line_loss_outage_wt_or_no_ref_adjs,
            :system_sol);
             NamedTupleTools.select(
                 line_loss_outage_wt_or_no_ref_adjs ,
                 (:state_labels,
                  :algebraic_vars_labels,

                  :dyn_pf_fun_kwd_n2s_idxs,
                  :dyn_pf_fun_kwd_net_idxs,

                  :Ybr_cal_and_edges_orientation))... )

branches_current_type_plot =
    make_plot_branches_current_type(
        getproperty(auxilliary_results, :t),
        getproperty(auxilliary_results, :t_i_b),
        getproperty(auxilliary_results, :t_α);
        ib_type = :i_b,
        α_type = :α)


streamedlined_auxilliary_results_plots =
    make_plot_streamedlined_auxilliary_results(
        ; NamedTupleTools.select(
            get_sol_auxilliary_results(
                getproperty(
                    line_loss_outage_wt_or_no_ref_adjs,
                    :system_sol);
                NamedTupleTools.select(
                    line_loss_outage_wt_or_no_ref_adjs,
                    (:state_labels,
                     :algebraic_vars_labels,

                     :dyn_pf_fun_kwd_n2s_idxs,
                     :dyn_pf_fun_kwd_net_idxs,

                     :Ybr_cal_and_edges_orientation))...),
            (:t,
             :gens_nodes_idx,
             :vd,
             :vq,
             :ph,
             :qh,
             :i_b,
             :α,))...)

#---------------------------------------------------

(edges_r,
 edges_x,
 Ybr_cal_and_edges_orientation) =
    NamedTupleTools.select(
        line_loss_outage_wt_or_no_ref_adjs,
        (:edges_r,
         :edges_x,
         :Ybr_cal_and_edges_orientation))

(edges_orientation,) =
    NamedTupleTools.select(
        Ybr_cal_and_edges_orientation,
        (:edges_orientation,))

(time, vh, θh, from_idxs, to_idxs, i_b, t_i_b) =
    NamedTupleTools.select(
        auxilliary_results,
        (:t, :vh, :θh, :from_idxs, :to_idxs, :i_b, :t_i_b))

#---------------------------------------------------

vec_t_i_b_t0 = t_i_b[:,1]

sum(vec_t_i_b_t0.^2 .* edges_r)

vec_i_b_t0 = i_b[:,1]

sum(vec_i_b_t0.^2 .* edges_r)


vec_vh = vh[:, 1]

vec_θh = θh[:, 1]

vec_uh = vec_vh .* exp.(im * vec_θh)

edges_I =
    (edges_r .+ im * edges_x) .* (
        vec_uh[from_idxs] - vec_uh[to_idxs])


edges_I_abs = abs.(edges_I)

sum(edges_I_abs.^2 .* edges_r)

#---------------------------------------------------
"""

# propertynames( line_loss_outage_wt_or_no_ref_adjs )

# (:system_sol, :model_syms, :gens_nodes_names, :SM_gens_nodes_names, :non_gens_nodes_names, :sim_timespan, :state_labels, :algebraic_vars_labels, :dyn_pf_fun_kwd_n2s_idxs, :dyn_pf_fun_kwd_net_idxs, :edges_r, :edges_x, :edges_b, :edges_ratio, :edges_angle, :edges_type, :Gs, :Bs, :Ybr_cal_and_edges_orientation, :Ynet_wt_nodes_idx_wt_adjacent_nodes)

# propertynames( auxilliary_results )

# (:t, :gens_nodes_idx, :vd, :vq, :ph, :qh, :vh, :θh, :from_idxs, :to_idxs, :i_b, :α, :I_from_I_to, :f_i_b, :f_α, :t_i_b, :t_α)


"""

(t_I_from_I_to,) =
    NamedTupleTools.select(
        auxilliary_results,
        (:I_from_I_to,))

last(split(string(:f_i_b),"_"))

plt_f_t_i_b = plot(
    getproperty(auxilliary_results, :t),
    getproperty(auxilliary_results, :t_i_b)', lw=1 )


plt_f_t_α = plot(
    getproperty(auxilliary_results, :t),
    getproperty(auxilliary_results, :t_α)', lw=1 )


auxilliary_results =
    get_sol_auxilliary_results(
        getproperty(
            line_loss_outage_wt_or_no_ref_adjs,
            :system_sol);
             NamedTupleTools.select(
                 line_loss_outage_wt_or_no_ref_adjs ,
                 (:state_labels,
                  :algebraic_vars_labels,

                  :dyn_pf_fun_kwd_n2s_idxs,
                  :dyn_pf_fun_kwd_net_idxs,

                  :Ybr_cal_and_edges_orientation))... )


(t_state_labels,                 
t_algebraic_vars_labels,

t_dyn_pf_fun_kwd_n2s_idxs,
t_dyn_pf_fun_kwd_net_idxs,
 t_Ybr_cal_and_edges_orientation,) =
     NamedTupleTools.select(
          line_loss_outage_wt_or_no_ref_adjs ,
         ( :state_labels,
           :algebraic_vars_labels,

           :dyn_pf_fun_kwd_n2s_idxs,
           :dyn_pf_fun_kwd_net_idxs,

           :Ybr_cal_and_edges_orientation,))


n2s_all_nodes_idx =
    NamedTupleTools.select(
        t_dyn_pf_fun_kwd_n2s_idxs,
        :n2s_all_nodes_idx)

all_nodes_idx =
    NamedTupleTools.select(
        t_dyn_pf_fun_kwd_net_idxs,
        :all_nodes_idx)

t_vh = ones(length(all_nodes_idx), 3)

t_θh = zeros(length(all_nodes_idx), 3)

( t_edges_Ybr_cal,
  t_edges_orientation) =
    NamedTupleTools.select(
        t_Ybr_cal_and_edges_orientation,
        (:edges_Ybr_cal,
         :edges_orientation))

t_I_from_I_to_by_vh_θh =
    get_vec_I_from_I_to_by_vh_θh(
        t_vh, t_θh ;
        n2s_all_nodes_idx,
        NamedTupleTools.select(
        t_Ybr_cal_and_edges_orientation,
        (:edges_Ybr_cal,
         :edges_orientation))... )



_, t_ncol_If_It = size(t_I_from_I_to_by_vh_θh)

t_I_from_I_to_by_vh_θh[:,2]

t_branches_currents =
    [[map( (a_if, a_it) -> a_if + a_it, a_if_vec, a_it_vec )
      for (a_if_vec, a_it_vec) in
          zip(first.(first.(t_I_from_I_to_by_vh_θh[:,idx])),
              second.(first.(t_I_from_I_to_by_vh_θh[:,idx])))] 
     for idx in 1:t_ncol_If_It ]



t_branches_currents =
    [map( (a_if, a_it) -> a_if .+ a_it, a_if_vec, a_it_vec )
      for (a_if_vec, a_it_vec) in
          zip(first.(first.(t_I_from_I_to_by_vh_θh)),
              second.(first.(t_I_from_I_to_by_vh_θh)))]


first.(first.(t_I_from_I_to_by_vh_θh[:,1]))
second.(first.(t_I_from_I_to_by_vh_θh[:,1]))

t_branches_currents_by_vh_θh =
    get_vec_branches_currents_by_vh_θh(
        t_vh, t_θh ;
        n2s_all_nodes_idx,
        NamedTupleTools.select(
        t_Ybr_cal_and_edges_orientation,
        (:edges_Ybr_cal,
         :edges_orientation))... )

"""
#---------------------------------------------------
#---------------------------------------------------


"""

sol_auxilliary_results =
    get_sol_auxilliary_results(
        an_outage_type_sol.system_sol;
        NamedTupleTools.select(
            an_outage_type_sol,
            (:state_labels,
             :algebraic_vars_labels,
             :dyn_pf_fun_kwd_n2s_idxs,
             :dyn_pf_fun_kwd_net_idxs))...)

"""

#---------------------------------------------------
#---------------------------------------------------

"""
# :line_outage

line_outage_plot =
    get_guick_group_vars_plots_dae_or_ode_sol(
        ;include_v_θ_plot =
            true,
        get_line_loss_outage_wt_or_no_ref_adjs(
            :line_outage,
        line_outage_time,
        generation_adjustment_time;
        case_name = "case9",
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
            line_data_in_pu )...)


# :line_outage_wt_pref_adjs


line_outage_wt_pref_adjs_plot =
    get_guick_group_vars_plots_dae_or_ode_sol(
        ;include_v_θ_plot =
            true,
        get_line_loss_outage_wt_or_no_ref_adjs(
            :line_outage_wt_pref_adjs,
        line_outage_time,
        generation_adjustment_time;
        case_name = "case9",
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
            line_data_in_pu )...)


# :line_outage_wt_vpref_adjs


line_outage_wt_vpref_adjs_plot =
    get_guick_group_vars_plots_dae_or_ode_sol(
        ;include_v_θ_plot =
            true,
        get_line_loss_outage_wt_or_no_ref_adjs(
            :line_outage_wt_vpref_adjs,
        line_outage_time,
        generation_adjustment_time;
        case_name = "case9",
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
            line_data_in_pu )...)

"""

#----------------------------------------
#----------------------------------------

# # :line_outage_wt_pref_adjs
# # :line_outage_wt_vpref_adjs
# # :line_outage

# # outage_type = :line_outage

# line_outage_plot =
#     get_guick_group_vars_plots_dae_or_ode_sol(
#         ;include_v_θ_plot =
#             true,
#          get_line_loss_outage_wt_or_no_ref_adjs(
#              :line_outage,
#              ntuple_status_steady_state_data,
#              line_outage_time,
#              generation_adjustment_time; 
#              sim_timespan = sim_timespan,    
#              dae_alg = dae_alg,
#              abstol  = abstol,
#              reltol  = reltol )...)


# line_outage_wt_pref_adjs_plot =
#     get_guick_group_vars_plots_dae_or_ode_sol(
#         ;include_v_θ_plot =
#             true,
#          get_line_loss_outage_wt_or_no_ref_adjs(
#              :line_outage_wt_pref_adjs,
#              ntuple_status_steady_state_data,
#              line_outage_time,
#              generation_adjustment_time; 
#              sim_timespan = sim_timespan,    
#              dae_alg = dae_alg,
#              abstol  = abstol,
#              reltol  = reltol )...)


# line_outage_wt_vpref_adjs_plot =
#     get_guick_group_vars_plots_dae_or_ode_sol(
#         ;include_v_θ_plot =
#             true,
#          get_line_loss_outage_wt_or_no_ref_adjs(
#              :line_outage_wt_vpref_adjs,
#              ntuple_status_steady_state_data,
#              line_outage_time,
#              generation_adjustment_time; 
#              sim_timespan = sim_timespan,    
#              dae_alg = dae_alg,
#              abstol  = abstol,
#              reltol  = reltol )...)


# line_loss_outage_wt_or_no_ref_adjs =
#     get_line_loss_outage_wt_or_no_ref_adjs(
#         :line_outage,
#         ntuple_status_steady_state_data,
#         line_outage_time,
#         generation_adjustment_time; 
#         sim_timespan = sim_timespan,    
#         dae_alg = dae_alg,
#         abstol  = abstol,
#         reltol  = reltol )


# #----------------------------------------
# #----------------------------------------

# (;system_fault_status,
#  generic_system_dynamics_wt_fault_kwd_para,
#  Ynet_wt_nodes_idx_wt_adjacent_nodes,
#  on_fault_net_para,
#  cleared_selected_lines_faults_net_para,

#  ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
 
#  model_bool_dae_vars_wt_fault,
#  model_syms_wt_fault,         
#  u0_model_states_init_wt_fault,

#  cb_states,
#  plants_cb_paras_switches,

#  nodes_names,
#  gens_nodes_names,
#  non_gens_nodes_names,
#  SM_gens_nodes_names,
#  SC_gens_nodes_names,

#  ωref0_vref0_porder0_id_iq_vh_Idx,
#  dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

#  dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx) =
#     NamedTupleTools.select(
#         getproperty(
#             getproperty(
#                 ntuple_status_steady_state_data,
#                 :pre_fault_state),
#             :static_prefault_paras),        
#         (:system_fault_status,
#          :generic_system_dynamics_wt_fault_kwd_para,
#          :Ynet_wt_nodes_idx_wt_adjacent_nodes,
#          :on_fault_net_para,
#          :cleared_selected_lines_faults_net_para,

#          :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

#          :model_bool_dae_vars_wt_fault,
#          :model_syms_wt_fault,         
#          :u0_model_states_init_wt_fault,

#          :cb_states,
#          :plants_cb_paras_switches,

#          :nodes_names,
#          :gens_nodes_names,
#          :non_gens_nodes_names,
#          :SM_gens_nodes_names,
#          :SC_gens_nodes_names,

#          :ωref0_vref0_porder0_id_iq_vh_Idx,
#          :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
         
#          :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx))

# #----------------------------------------

# (post_outage_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,) =
#     NamedTupleTools.select(
#         getproperty(
#             getproperty(
#                 ntuple_status_steady_state_data,
#                 :post_fault_state),
#             :dynamic_status_paras),
#         (:ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,))

# #----------------------------------------
# #----------------------------------------

# (Ynet, ) =
#      NamedTupleTools.select(
#     Ynet_wt_nodes_idx_wt_adjacent_nodes,
#          (:Ynet, ) )

# #----------------------------------------

# # (faulty_Ynet,
# #  faulty_nodes_idx_with_adjacent_nodes_idx) =
# #     NamedTupleTools.select(
# #     on_fault_net_para,
# #         (:faulty_Ynet,
# #          :faulty_nodes_idx_with_adjacent_nodes_idx))

# (fault_Ynet,
#  post_fault_Ynet) =
#     NamedTupleTools.select(
#     cleared_selected_lines_faults_net_para,
#         (:pre_clear_fault_Ynet,
#          :post_clear_fault_Ynet,))

# #----------------------------------------

# (;dyn_ω_ref_Idx,
#  dyn_v_ref_Idx,
#  dyn_p_order_Idx,
#  dyn_Png_Idx,
#  dyn_Qng_Idx,
#  dyn_Pll_Idx,
#  dyn_Qll_Idx ) =
#      NamedTupleTools.select(
#          dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
#          (:dyn_ω_ref_Idx,
#           :dyn_v_ref_Idx,
#           :dyn_p_order_Idx,
#           :dyn_Png_Idx,
#           :dyn_Qng_Idx,
#           :dyn_Pll_Idx,
#           :dyn_Qll_Idx))

# #---------------------------------------
# #---------------------------------------
# # line loss with only porder_adj
# #---------------------------------------
# #---------------------------------------

# gens_porder_adj =
#     post_outage_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
#         dyn_p_order_Idx]

# #---------------------------------------

# system_fault_status[1] = 0

# cb_line_outage = DiscreteCallback(
#     (u, t, integrator) ->
#         on_line_outage_condition(
#             u, t, integrator,
#             line_outage_time),

#    on_line_outage_affect!;
#     save_positions=(true, true),
#     initializealg =
#         ShampineCollocationInit() )


# cb_gens_porder_adjustment = DiscreteCallback(
#     (u, t, integrator) ->
#         on_generation_adjustment_condition(
#             u, t, integrator,
#             generation_adjustment_time),

#     (integrator) -> on_generation_adjustment_affect!(
#         integrator,
#         gens_porder_adj,
#         dyn_p_order_Idx );
#     save_positions=(true, true),
#     initializealg =
#         ShampineCollocationInit() )


# cb_line_outage_and_cb_gens_porder_adjustment =
#     CallbackSet(cb_line_outage,
#     cb_gens_porder_adjustment)

# tstop_line_outage_and_gens_porder_adjustment =
#     [line_outage_time,
#      generation_adjustment_time]

# #---------------------------------------------------

# generic_model_dynamics_para =
#     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

# # system_fault_status[1] = 0

# model_dynamics_para =
#     (;generic_model_dynamics_para,
#      Ynet,
#      fault_Ynet,
#      post_fault_Ynet,
#      system_fault_status,
#      plants_cb_paras_switches)

# model_dynamics_kwd_para =
#     generic_system_dynamics_wt_fault_kwd_para

# #----------------------------------------

# model_bool_dae_vars =
#     model_bool_dae_vars_wt_fault

# model_syms =
#     model_syms_wt_fault

# u0_model_states_init =
#     u0_model_states_init_wt_fault

# #----------------------------------------

# du0_model_states_init =
#     zeros(length(u0_model_states_init))

# res = similar(u0_model_states_init)

# #---------------------------------------------------

# system_dynamics_fun! =
#     line_loss_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!

# #---------------------------------------------------

# system_sol =
#     DifferentialEquations.solve(
#         DAEProblem(
#     DAEFunction(
#     (res, dx, x, p, t) ->
#         system_dynamics_fun!(
#             res, dx, x,
#             model_dynamics_para,
#             t;
#             kwd_para =
#                 model_dynamics_kwd_para);
#     syms =
#         model_syms),
#     du0_model_states_init,
#     u0_model_states_init,
#     sim_timespan,
#     model_dynamics_para,
#     differential_vars =
#         model_bool_dae_vars,
#         callback = cb_states),
#         dae_alg,
#         callback =
#             cb_line_outage_and_cb_gens_porder_adjustment,
#         tstops =
#             tstop_line_outage_and_gens_porder_adjustment,
#         abstol = abstol,
#         reltol = reltol )

# line_outage_wt_gens_porder_adjs_plot =
#     get_guick_group_vars_plots_dae_or_ode_sol(
#         ;include_v_θ_plot =
#             true,
#          sim_timespan =
#              (0, timespan),
#         system_sol,
#         model_syms,
#         gens_nodes_names,
#         SM_gens_nodes_names,
#         non_gens_nodes_names)



# #----------------------------------------
# #----------------------------------------
# # line loss with vref_and_porder_adj
# #---------------------------------------
# #----------------------------------------


# (;
#  Ynet_wt_nodes_idx_wt_adjacent_nodes,
#  on_fault_net_para,
#  cleared_selected_lines_faults_net_para,

#  ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
 
#  model_bool_dae_vars_wt_fault,
#  u0_model_states_init_wt_fault,

#  cb_states,
#  plants_cb_paras_switches) =
#     NamedTupleTools.select(
#         getproperty(
#             getproperty(
#                 ntuple_status_steady_state_data,
#                 :pre_fault_state),
#             :static_prefault_paras),        
#         (
#          :Ynet_wt_nodes_idx_wt_adjacent_nodes,
#          :on_fault_net_para,
#          :cleared_selected_lines_faults_net_para,

#          :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

#          :model_bool_dae_vars_wt_fault,
#          :u0_model_states_init_wt_fault,

#          :cb_states,
#          :plants_cb_paras_switches ))

# (post_outage_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,) =
#     NamedTupleTools.select(
#         getproperty(
#             getproperty(
#                 ntuple_status_steady_state_data,
#                 :post_fault_state),
#             :dynamic_status_paras),
#         (:ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,))

# v_ref_adj =
#     post_outage_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
#         dyn_v_ref_Idx]

# #---------------------------------------

# system_fault_status[1] = 0

# #---------------------------------------

# generic_model_dynamics_para =
#     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

# # system_fault_status[1] = 0

# model_dynamics_para =
#     (;generic_model_dynamics_para,
#      Ynet,
#      fault_Ynet,
#      post_fault_Ynet,
#      system_fault_status,
#      plants_cb_paras_switches)

# model_dynamics_kwd_para =
#     generic_system_dynamics_wt_fault_kwd_para

# #----------------------------------------

# model_bool_dae_vars =
#     model_bool_dae_vars_wt_fault

# model_syms =
#     model_syms_wt_fault

# u0_model_states_init =
#     u0_model_states_init_wt_fault

# #----------------------------------------

# du0_model_states_init =
#     zeros(length(u0_model_states_init))

# res = similar(u0_model_states_init)

# #---------------------------------------------------

# system_dynamics_fun! =
#     line_loss_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!

# #---------------------------------------------------

# system_sol =
#     DifferentialEquations.solve(
#         DAEProblem(
#     DAEFunction(
#     (res, dx, x, p, t) ->
#         system_dynamics_fun!(
#             res, dx, x,
#             model_dynamics_para,
#             t;
#             kwd_para =
#                 model_dynamics_kwd_para);
#     syms =
#         model_syms),
#     du0_model_states_init,
#     u0_model_states_init,
#     sim_timespan,
#     model_dynamics_para,
#     differential_vars =
#         model_bool_dae_vars,
#         callback = cb_states),
#         dae_alg,
#         callback =
#             cb_line_outage_and_cb_vref_and_porder_adj,
#         tstops =
#             tstop_line_outage_and_vref_and_porder_adj,
#         abstol = abstol,
#         reltol = reltol )

# line_outage_and_vref_and_porder_adj_plot =
#     get_guick_group_vars_plots_dae_or_ode_sol(
#         ;include_v_θ_plot =
#             true,
#          sim_timespan =
#              (0, timespan),
#         system_sol,
#         model_syms,
#         gens_nodes_names,
#         SM_gens_nodes_names,
#         non_gens_nodes_names)

# #---------------------------------------------------
# # line outage no restoration 
# #---------------------------------------------------

# system_fault_status[1] = 0

# cb_line_outage = DiscreteCallback(
#     (u, t, integrator) ->
#         on_line_outage_condition(
#             u, t, integrator,
#             line_outage_time),

#    on_line_outage_affect!;
#     save_positions=(true, true),
#     initializealg =
#         ShampineCollocationInit() )


# cb_set = CallbackSet(
#     cb_line_outage,
#     cb_states)


# faults_times =
#     [ line_outage_time ]

# model_dynamics_para =
#     (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
#      Ynet,
#      fault_Ynet,
#      post_fault_Ynet,
#      system_fault_status )          

# model_dynamics_kwd_para =
#     generic_system_dynamics_wt_fault_kwd_para

# system_dynamics_fun! =
#     Ynet_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!


# #----------------------------------------

# model_bool_dae_vars =
#     model_bool_dae_vars_wt_fault

# model_syms =
#     model_syms_wt_fault

# u0_model_states_init =
#     u0_model_states_init_wt_fault

# #----------------------------------------

# du0_model_states_init =
#     zeros(length(u0_model_states_init))

# res = similar(u0_model_states_init)

# #----------------------------------------
# #----------------------------------------

# system_sol =
#     DifferentialEquations.solve(
#         DAEProblem(
#     DAEFunction(
#     (res, dx, x, p, t) ->
#         system_dynamics_fun!(
#             res, dx, x,
#             model_dynamics_para,
#             t;
#             kwd_para =
#                 model_dynamics_kwd_para);
#     syms =
#         model_syms),
#     du0_model_states_init,
#     u0_model_states_init,
#     sim_timespan,
#     model_dynamics_para,
#     differential_vars =
#         model_bool_dae_vars,
#         callback = cb_states),
#         dae_alg,
#         callback = cb_line_outage,
#         tstops = faults_times,
#         abstol = abstol,
#         reltol = reltol )


# line_outage_no_adjs_plot =
#     get_guick_group_vars_plots_dae_or_ode_sol(
#         ;include_v_θ_plot =
#             true,
#          sim_timespan =
#              (0, timespan),
#         system_sol,
#         model_syms,
#         gens_nodes_names,
#         SM_gens_nodes_names,
#         non_gens_nodes_names)


#---------------------------------------------------
#---------------------------------------------------

# system_fault_status[1] = 0



# static_prefault_paras =
#     getproperty(
#         getproperty(
#             ntuple_status_steady_state_data,
#             :pre_fault_state),
#             :static_prefault_paras)

#----------------------------------------

# fault_paras =
#     getproperty(
#         getproperty(
#             ntuple_status_steady_state_data,
#             :fault_state),
#         :dynamic_status_paras)

# #----------------------------------------

# post_fault_paras =
#     getproperty(
#         getproperty(
#             ntuple_status_steady_state_data,
#             :post_fault_state),
#         :dynamic_status_paras)



# (s_pf_P_gens,
#  s_pf_Q_gens,
#  s_vh,
#  s_θh,
#  s_gens_vh,
#  s_gens_θh,
#  s_gens_id,
#  s_gens_iq,
 
#  s_gens_mag_E,
#  s_gens_ang_E,
#  s_post_sta_PQ,
#  s_Yred,
#  s_Yint,

#  s_ω_ref, s_v_ref, s_p_order,
#  s_gens_i_d, s_gens_i_q,
 
#  s_flat_vh_flat_θh_id_iq_u0,
#  s_flat_vh_flat_θh_id_iq_vfh_θfh,
#  s_gens_δ,
#  s_gens_ed_dash,
#  s_gens_eq_dash) =
#      NamedTupleTools.select(
#         getproperty(
#             getproperty(
#                 ntuple_status_steady_state_data,
#                 :pre_fault_state),
#             :static_prefault_paras),
#         (:pf_P_gens,
#          :pf_Q_gens,
#          :vh,
#          :θh,
#          :gens_vh,
#          :gens_θh,
#          :gens_id,
#          :gens_iq,
         
#          :gens_mag_E,
#          :gens_ang_E,
#          :post_sta_PQ,
#          :Yred,
#          :Yint,
         
#          :ω_ref, :v_ref, :p_order,
#          :gens_i_d, :gens_i_q,
         
#          :flat_vh_flat_θh_id_iq_u0,
#          :flat_vh_flat_θh_id_iq_vfh_θfh,
#          :gens_δ,
#          :gens_ed_dash,
#          :gens_eq_dash ) )


# (f_dyn_pf_P_gens,
#  f_dyn_pf_Q_gens,
#  f_dyn_vh,
#  f_dyn_θh,
#  f_dyn_gens_vh,
#  f_dyn_gens_θh,
#  f_dyn_gens_id,
#  f_dyn_gens_iq,
#  f_dyn_gens_mag_E,
#  f_dyn_gens_ang_E,
#  f_post_dyn_PQ,
#  f_dyn_Yred,
#  f_dyn_Yint,
#  f_flat_vh_flat_θh_id_iq_vfh_θfh,
#  f_dyn_gens_δ,
#  f_dyn_gens_ed_dash,
#  f_dyn_gens_eq_dash) =
#      NamedTupleTools.select(
#         getproperty(
#             getproperty(
#                 ntuple_status_steady_state_data,
#                 :fault_state),
#             :dynamic_status_paras),
#         (:dyn_pf_P_gens,
#          :dyn_pf_Q_gens,
#          :dyn_vh,
#          :dyn_θh,
#          :dyn_gens_vh,
#          :dyn_gens_θh,
#          :dyn_gens_id,
#          :dyn_gens_iq,
#          :dyn_gens_mag_E,
#          :dyn_gens_ang_E,
#          :post_dyn_PQ,
#          :dyn_Yred, :dyn_Yint,
#          :flat_vh_flat_θh_id_iq_vfh_θfh,
#          :dyn_gens_δ,
#          :dyn_gens_ed_dash,
#          :dyn_gens_eq_dash ) )

# #----------------------------------------

# (p_dyn_pf_P_gens,
#  p_dyn_pf_Q_gens,
#  p_dyn_vh,
#  p_dyn_θh,
#  p_dyn_gens_vh,
#  p_dyn_gens_θh,
#  p_dyn_gens_id,
#  p_dyn_gens_iq,
#  p_dyn_gens_mag_E,
#  p_dyn_gens_ang_E,
#  p_post_dyn_PQ,
#  p_dyn_Yred,
#  p_dyn_Yint,
#  p_flat_vh_flat_θh_id_iq_vfh_θfh,
#  p_dyn_gens_δ,
#  p_dyn_gens_ed_dash,
#  p_dyn_gens_eq_dash) =
#      NamedTupleTools.select(
#         getproperty(
#             getproperty(
#                 ntuple_status_steady_state_data,
#                 :post_fault_state),
#             :dynamic_status_paras),
#         (:dyn_pf_P_gens,
#          :dyn_pf_Q_gens,
#          :dyn_vh,
#          :dyn_θh,
#          :dyn_gens_vh,
#          :dyn_gens_θh,
#          :dyn_gens_id,
#          :dyn_gens_iq,
#          :dyn_gens_mag_E,
#          :dyn_gens_ang_E,
#          :post_dyn_PQ,
#          :dyn_Yred, :dyn_Yint,
#          :flat_vh_flat_θh_id_iq_vfh_θfh,
#          :dyn_gens_δ,
#          :dyn_gens_ed_dash,
#          :dyn_gens_eq_dash) )


# mismatch_faulty_fault_Ynet =
#     map(x -> round.(x; digits=4),
#         [t1_array - t2_array
#          for (t1_array, t2_array) in
#              zip(faulty_Ynet,
#                  fault_Ynet)])

# t_Ynet =
#     map(x -> round.(x; digits=4),
#         Ynet)

# t_faulty_Ynet =
#     map(x -> round.(x; digits=4),
#         faulty_Ynet)

# t_post_fault_Ynet =
#     map(x -> round.(x; digits=4),
#         post_fault_Ynet)

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

