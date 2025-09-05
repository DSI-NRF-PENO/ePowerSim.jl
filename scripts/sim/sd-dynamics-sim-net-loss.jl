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

#---------------------------------------------------
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

sim_type = "net-loss-and-sensitivity"

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
## solvers and settings
#---------------------------------------------------

use_init_u0 = false

use_nlsolve = false

pf_alg        = NewtonRaphson()

#---------------------------------------------------

ode_alg       = Rodas4()

# ode_alg       = ImplicitMidpoint()

dae_alg       = IDA()

abstol        = 1e-12

reltol        = 1e-12

#---------------------------------------------------
# base setting and some booleans 
#---------------------------------------------------

# """

# basekV = 1.0

# use_pu_in_PQ = true

# line_data_in_pu = true

# #---------------------------------------------------
# ## fault data
# #---------------------------------------------------

# on_fault_time = 9.0

# clear_fault_time = 9.02

# with_faults = false

# list_fault_point_from_node_a = [0.01]

# list_fault_resistance = [0.001]

# list_no_line_circuit =  [1]

# list_edges_to_have_fault = [ 2 ]

# clear_fault_selection_list = [1]

# """

#---------------------------------------------------
# web resources
#---------------------------------------------------

# https://www.wias-berlin.de/people/fuhrmann/AdSciComp-WS2324/assets/nb09-ode-dae.html


(vh,
 θh,x
 u0_model_states_init,
 model_syms,
 model_mass_matrix,

 ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
 plants_cb_paras_switches,

 generic_system_dynamics_kwd_para,

 gens_nodes_names,
 SM_gens_nodes_names,
 non_gens_nodes_names,

 cb_states,
 
 dyn_pf_fun_kwd_n2s_idxs,
 dyn_pf_fun_kwd_net_idxs,
 
 Pg_Qg_Png_Qng_Pll_Qll_Idx,
 Png_Qng_Pll_Qll_Idx,
 Pg_Png_Qng_Idx,

 pf_vh_θh_idx_and_idx2Idx,
 
 dyn_pf_flat_vh_flat_θh_Idx,

 edges_r,
 edges_x,
 edges_b,
 edges_ratio,
 edges_angle,
 edges_type,
 Gs,
 Bs,
 
 Ynet_wt_nodes_idx_wt_adjacent_nodes,
 Ybr_cal_and_edges_orientation,

 Pg_Qg_Png_Qng_Pll_Qll,
 loc_load_exist,
 slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

 edges_r_x_b_ratio_angle_idx,
 Ynet_rows_Idxs_in_flattend,
 Ynet_real_imag_Idxs_in_flattend,
 
 scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
 scale_Pg_Png_Qng_Idx,
 dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

 sta_pf_PQ_para,
 generic_gens_para,
 ode_gens_para ) =
     NamedTupleTools.select(
         getproperty(
        get_a_status_steady_state_data(
        :pre_fault_state;
        with_faults = false,
        net_data_by_components_file =
            net_data_by_components_file,
        components_libs_dir =
            components_libs_dir,        
        
        timespan   = timespan,

        on_fault_time = 5.0,
        clear_fault_time = 7.0,

        list_fault_point_from_node_a = [0.3],
        list_fault_resistance = [0.001],
        list_no_line_circuit =  [1],

        list_edges_to_have_fault = [ 8 ],
        clear_fault_selection_list = [1],

        basekV = 1.0,    
        use_pu_in_PQ = true,
        line_data_in_pu = true),
        :static_prefault_paras),             
         (:vh,
          :θh,
          :u0_model_states_init,
          :model_syms,
          :model_mass_matrix,

          :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
          :plants_cb_paras_switches,

          :generic_system_dynamics_wt_fault_kwd_para,

          :gens_nodes_names,
          :SM_gens_nodes_names,
          :non_gens_nodes_names,

          :cb_states,
          
          :dyn_pf_fun_kwd_n2s_idxs,
          :dyn_pf_fun_kwd_net_idxs,
          
          :Pg_Qg_Png_Qng_Pll_Qll_Idx,
          :Png_Qng_Pll_Qll_Idx,
          :Pg_Png_Qng_Idx,
          
          :pf_vh_θh_idx_and_idx2Idx,
          
          :dyn_pf_flat_vh_flat_θh_Idx,

          :edges_r,
          :edges_x,
          :edges_b,
          :edges_ratio,
          :edges_angle,
          :edges_type,
          :Gs,
          :Bs,
          
          :Ynet_wt_nodes_idx_wt_adjacent_nodes,
          :Ybr_cal_and_edges_orientation,
          
          :Pg_Qg_Png_Qng_Pll_Qll,
          :loc_load_exist,
          
          :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

          :edges_r_x_b_ratio_angle_idx,
          :Ynet_rows_Idxs_in_flattend,
          :Ynet_real_imag_Idxs_in_flattend,

          :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
          :scale_Pg_Png_Qng_Idx,
          :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

          :sta_pf_PQ_para,
          :generic_gens_para,
          :ode_gens_para))



(non_gens_vh_idx,
 non_slack_gens_θh_idx,
 non_gens_θh_idx,

 red_vh_Idxs,
 red_non_slack_gens_θh_idx2Idx,
 red_non_gens_θh_idx2Idx,

 red_non_gens_vh_Idxs,
 red_non_slack_gens_θh_Idxs,
 red_non_gens_θh_Idxs,

 red_slack_value_Idxs,

 gens_nodes_idx,

 non_gens_nodes_idx,
 non_slack_gens_and_non_gens_idx,
 gens_nodes_with_loc_loads_idx,
 all_nodes_idx,
 
 n2s_gens_idx,
 n2s_gens_with_loc_load_idxs,
 n2s_non_gens_idx,
 n2s_all_nodes_idx ) =
     NamedTupleTools.select(
         pf_vh_θh_idx_and_idx2Idx,
         (:non_gens_vh_idx,
          :non_slack_gens_θh_idx,
          :non_gens_θh_idx,

          :red_vh_Idxs,
          :red_non_slack_gens_θh_idx2Idx,
          :red_non_gens_θh_idx2Idx,

          :red_non_gens_vh_Idxs,
          :red_non_slack_gens_θh_Idxs,
          :red_non_gens_θh_Idxs,

          :red_slack_value_Idxs,

          :gens_nodes_idx,
          
          :non_gens_nodes_idx,
          :non_slack_gens_and_non_gens_idx,
          :gens_nodes_with_loc_loads_idx,
          :all_nodes_idx,

          :n2s_gens_idx,
          :n2s_gens_with_loc_load_idxs,
          :n2s_non_gens_idx,
          :n2s_all_nodes_idx))


(;Ynet,
 nodes_idx_with_adjacent_nodes_idx) =
     NamedTupleTools.select(
         Ynet_wt_nodes_idx_wt_adjacent_nodes,
         (:Ynet,
          :nodes_idx_with_adjacent_nodes_idx))


(;edges_Ybr_cal,
 edges_orientation) =
     NamedTupleTools.select(
         Ybr_cal_and_edges_orientation,
         (:edges_Ybr_cal,
          :edges_orientation ))


#--------------------------------------------
#--------------------------------------------


(;dyn_P_gens_Idxs,
 dyn_Q_gens_Idxs,
 dyn_P_non_gens_Idxs,
 dyn_Q_non_gens_Idxs,
 dyn_P_gens_loc_load_Idxs,
 dyn_Q_gens_loc_load_Idxs) =
     NamedTupleTools.select(
         Pg_Qg_Png_Qng_Pll_Qll_Idx,
         (:dyn_P_gens_Idxs,
          :dyn_Q_gens_Idxs,
          :dyn_P_non_gens_Idxs,
          :dyn_Q_non_gens_Idxs,
          :dyn_P_gens_loc_load_Idxs,
          :dyn_Q_gens_loc_load_Idxs))


#--------------------------------------------


red_vh_θh_idx =
    getproperty(
        pf_vh_θh_idx_and_idx2Idx,
        :red_vh_θh_idx)


#--------------------------------------------


(;Ynet_rows_Idxs_in_flattend,
 Ynet_real_imag_Idxs_in_flattend ) =
     NamedTupleTools.select(
         get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
             Ynet),
         (:Ynet_rows_Idxs_in_flattend,
          :Ynet_real_imag_Idxs_in_flattend))


(;Ynet_rows_Idxs_in_flattend,
Ynet_real_imag_Idxs_in_flattend ) =
    NamedTupleTools.select(
        get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
            getproperty(
                Ynet_wt_nodes_idx_wt_adjacent_nodes,
                :Ynet)),
        (:Ynet_rows_Idxs_in_flattend,
         :Ynet_real_imag_Idxs_in_flattend) )


(;Ynet_real_Idxs,
 Ynet_imag_Idxs) =
    NamedTupleTools.select(
        Ynet_real_imag_Idxs_in_flattend,
        (:Ynet_real_Idxs,
         :Ynet_imag_Idxs))


#--------------------------------------------


(;r_Idxs, x_Idxs, b_Idxs,
 ratio_Idxs, angle_Idxs) =
     NamedTupleTools.select(
         edges_r_x_b_ratio_angle_idx,
         (:r_Idxs, :x_Idxs, :b_Idxs,
          :ratio_Idxs, :angle_Idxs))


(vh_Idxs, θh_Idxs) =
    NamedTupleTools.select(
        dyn_pf_flat_vh_flat_θh_Idx,
        (:dyn_pf_vh_Idxs,
         :dyn_pf_θh_Idxs))


#--------------------------------------------
#--------------------------------------------

Pg =
    Pg_Qg_Png_Qng_Pll_Qll[dyn_P_gens_Idxs]

Qg =
    Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_gens_Idxs]

Png =
    Pg_Qg_Png_Qng_Pll_Qll[dyn_P_non_gens_Idxs]
 
Qng =
    Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_non_gens_Idxs]

Pll =
    Pg_Qg_Png_Qng_Pll_Qll[dyn_P_gens_loc_load_Idxs]

Qll =
    Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_gens_loc_load_Idxs]


#--------------------------------------------

gens_Sn = getproperty(ode_gens_para, :Sn)


gens_loss_participation =
    get_loss_particpation_by_gens_rating(
        gens_Sn, Pg)

gens_loss_participation_by_loading =
    get_loss_particpation_by_gens_loading( Pg)

gens_active_power_particpation_by_rating =
    get_gens_active_power_particpation_by_rating(
        gens_Sn, Pg)

gens_active_power_particpation_by_loading =
    get_gens_active_power_particpation_by_loading( Pg )

gens_reactive_power_particpation_by_rating =
    get_gens_reactive_power_particpation_by_rating(
        gens_Sn, Qg)

gens_reactive_power_particpation_by_loading =
    get_gens_reactive_power_particpation_by_loading( Qg)


active_power_disturbance_resolution_participation =
    gens_active_power_particpation_by_rating

reactive_power_disturbance_resolution_participation =
    gens_reactive_power_particpation_by_loading

pf_model_kwd_para =
    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     
     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     gens_loss_participation,
     active_power_disturbance_resolution_participation,
     reactive_power_disturbance_resolution_participation,

     sta_pf_PQ_para,     

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,
     pf_vh_θh_idx_and_idx2Idx,
     
     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     scale_Pg_Png_Qng_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx )


generic_red_sol_kwd_para =
    (;
     Ybr_cal_and_edges_orientation,
     ode_gens_para,
     sta_pf_PQ_para,
     pf_model_kwd_para )


#--------------------------------------------

Pg_inj_Qg_inj_Png_Qng_Pll_Qll =
    get_Pg_inj_Qg_inj_Png_Qng_Pll_Qll(
        Pg_Qg_Png_Qng_Pll_Qll;
        loc_load_exist,
        Pg_Qg_Png_Qng_Pll_Qll_Idx,    
        gens_nodes_idx,
        n2s_gens_idx,
        n2s_gens_with_loc_load_idxs)


Pg_inj_Qg_inj_Png_Qng =
    getproperty(
        Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
        :Pg_inj_Qg_inj_Png_Qng)

(;gens_Pg_inj,
 gens_Qg_inj,
 P_non_gens,
 Q_non_gens,
 Pg_inj_Qg_inj_Png_Qng) =
     NamedTupleTools.select(
         Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
         (:gens_Pg_inj,
          :gens_Qg_inj,
          :P_non_gens,
          :Q_non_gens,
          :Pg_inj_Qg_inj_Png_Qng))

#--------------------------------------------
#--------------------------------------------
# initial values and parameters
#--------------------------------------------
#--------------------------------------------


vh_θh =
    [vh;
     θh]

#--------------------------------------------


Pg_inj_Qg_inj_Png_Qng =
    [gens_Pg_inj;
     gens_Qg_inj;
     P_non_gens;
     Q_non_gens ]


Pg_inj_Png_Qng =
    [gens_Pg_inj;
     P_non_gens;
     Q_non_gens ]

slack_gens_nodes_idx =
    getproperty(
        dyn_pf_fun_kwd_net_idxs,
        :slack_gens_nodes_idx)[1]

# Subtract net loss from slack node Pg,
# loss will subsequently be distributed
 
ds_gens_Pg_inj = @set gens_Pg_inj[slack_gens_nodes_idx] =
    gens_Pg_inj[slack_gens_nodes_idx] - get_total_P_network_loss(
        vh,
        θh,
        Ynet;
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx,
        all_nodes_idx )


ds_Pg_inj_Qg_inj_Png_Qng =
    [ds_gens_Pg_inj;
     gens_Qg_inj;
     P_non_gens;
     Q_non_gens ]

ds_Pg_inj_Png_Qng =
    [ds_gens_Pg_inj;
     P_non_gens;
     Q_non_gens]

ds_scale_Pg_inj_Qn_inj_Png_Qng =
    [[1.0];
     ds_Pg_inj_Qg_inj_Png_Qng]

ds_scale_Pg_inj_Png_Qng =
    [[1.0];
     ds_Pg_inj_Png_Qng ]

vh_θh_slack_value =
    [vh_θh;
     [0.0]]

red_vh_θh_slack_value =
    [vh_θh[non_gens_vh_idx];
     vh_θh[non_slack_gens_θh_idx];
     vh_θh[non_gens_θh_idx];
     [0.0]]


non_gens_vh_all_θh_slack_value =
    [vh[non_gens_nodes_idx];
     θh;
     [0.0]]    

#--------------------------------------------
# Loss participation factor based on gens ratings
#--------------------------------------------

# Synchrnous condensers will not participate

gens_loss_participation_factor =
    getproperty(
        gens_loss_participation,
        :gens_loss_participation_factor)



#--------------------------------------------


net_P_loss = loc_load_exist == true ?
    sum(Pg) - sum(Pll) - sum(Png) : sum(Pg) - sum(Png)


net_Q_loss = loc_load_exist == true ?
    sum(Qg) - sum(Qll) - sum(Qng) : sum(Qg) - sum(Qng)

#--------------------------------------------

total_P_network_loss =
    get_total_P_network_loss(
        vh_θh,
        Ynet;
        dyn_pf_flat_vh_flat_θh_Idx,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx,
        all_nodes_idx )


total_Q_network_loss =
    get_total_Q_network_loss(
        vh_θh,
        Ynet;
        dyn_pf_flat_vh_flat_θh_Idx,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx,
        all_nodes_idx )


#--------------------------------------------


total_P_network_loss =
    get_total_P_network_loss(
        vh,
        θh,
        Ynet;
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx,
        all_nodes_idx )


total_Q_network_loss =
    get_total_Q_network_loss(
        vh,
        θh,
        Ynet;
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx,
        all_nodes_idx )


#--------------------------------------------


ph_pk_qh_qk =
    get_ph_pk_qh_qk(
        vh,
        θh;
        Ybr_cal_and_edges_orientation,
        n2s_all_nodes_idx)


p_h = round.(
    getproperty(
        ph_pk_qh_qk,
        :p_h); digits=4)


q_h = round.(
    getproperty(
        ph_pk_qh_qk,
        :q_h); digits=4)


p_k = round.(
    getproperty(
        ph_pk_qh_qk,
        :p_k); digits=4)


q_k = round.(
    getproperty(
        ph_pk_qh_qk,
        :q_k); digits=4)


#--------------------------------------------

total_line_loss =
    getproperty(
        get_lines_losses_cal(
            vh, θh;
            Ybr_cal_and_edges_orientation,
            n2s_all_nodes_idx),
        :total_line_loss)

#--------------------------------------------

total_line_loss =
    getproperty(
        get_lines_losses_cal(
            vh,
            θh,
            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle;    
            edges_type,
            edges_orientation,
            n2s_all_nodes_idx),
        :total_line_loss)

#--------------------------------------------

Ynet_real_imag_Idxs_in_flattend

Ynet_real = real.(Ynet)

Ynet_imag = imag.(Ynet)

Ynet_real_imag_flattend =
    [[real.(Ynet)...;];
     [imag.(Ynet)...;]]

#--------------------------------------------

t_df_dx = ForwardDiff.gradient(
    ( vh_θh_x ) ->
        get_total_P_network_loss_by_flattened_Ynet(
            vh_θh_x, Ynet_real_imag_flattend ;
            dyn_pf_flat_vh_flat_θh_Idx,
            nodes_idx_with_adjacent_nodes_idx,
            Ynet_rows_Idxs_in_flattend,
            Ynet_real_imag_Idxs_in_flattend,
            n2s_all_nodes_idx,
            all_nodes_idx ) ,
    vh_θh )


t_df_dp = ForwardDiff.gradient(
    ( p ) ->
        get_total_P_network_loss_by_flattened_Ynet(
            vh_θh, p ;
            dyn_pf_flat_vh_flat_θh_Idx,
            nodes_idx_with_adjacent_nodes_idx,
            Ynet_rows_Idxs_in_flattend,
            Ynet_real_imag_Idxs_in_flattend,
            n2s_all_nodes_idx,
            all_nodes_idx ) ,
    Ynet_real_imag_flattend )


# t_df_dx_ 

# dx_dp  = -(svd( t_df_dx )) \ t_df_dp

#--------------------------------------------
#  vh θh by per node losses sensitivity 
#--------------------------------------------

n_df_dx = ForwardDiff.jacobian(
    ( vh_θh_x ) ->
        get_get_losses_per_node_by_flattened_Ynet(
            vh_θh_x, Ynet_real_imag_flattend ;
            dyn_pf_flat_vh_flat_θh_Idx,
            nodes_idx_with_adjacent_nodes_idx,
            Ynet_rows_Idxs_in_flattend,
            Ynet_real_imag_Idxs_in_flattend,
            n2s_all_nodes_idx,
            all_nodes_idx ),
    vh_θh )

n_df_dp = ForwardDiff.jacobian(
    ( p ) ->
        get_get_losses_per_node_by_flattened_Ynet(
            vh_θh, p ;
            dyn_pf_flat_vh_flat_θh_Idx,
            nodes_idx_with_adjacent_nodes_idx,
            Ynet_rows_Idxs_in_flattend,
            Ynet_real_imag_Idxs_in_flattend,
            n2s_all_nodes_idx,
            all_nodes_idx ) ,
    Ynet_real_imag_flattend )

n_dx_dp  = -( svd( n_df_dx ) \ n_df_dp )

#--------------------------------------------
#  vh θh by per line losses sensitivity 
#--------------------------------------------

edges_r_x_b_ratio_angle =
    [edges_r;
     edges_x;
     edges_b;
     edges_ratio;
     edges_angle]

total_line_losses =
    sum( get_losses_per_line(
        vh_θh,
        edges_r_x_b_ratio_angle;    
        edges_type,
        edges_orientation,
        n2s_all_nodes_idx,
        dyn_pf_flat_vh_flat_θh_Idx,
        edges_r_x_b_ratio_angle_idx) )

df_dx = ForwardDiff.jacobian(
    ( vh_θh_x ) ->
        get_losses_per_line(
            vh_θh_x,
            edges_r_x_b_ratio_angle;    
            edges_type,
            edges_orientation,
            n2s_all_nodes_idx,
            dyn_pf_flat_vh_flat_θh_Idx,
            edges_r_x_b_ratio_angle_idx),
            vh_θh )

df_dp = ForwardDiff.jacobian(
    ( p ) ->
        get_losses_per_line(
            vh_θh,
            p;    
            edges_type,
            edges_orientation,
            n2s_all_nodes_idx,
            dyn_pf_flat_vh_flat_θh_Idx,
            edges_r_x_b_ratio_angle_idx) ,
    edges_r_x_b_ratio_angle )

dx_dp  = -( svd( df_dx ) \ df_dp )

#--------------------------------------------

# sensitivy to lines resistances

dx_dp_vhg_by_r  =
    round.( dx_dp[ vh_Idxs[gens_nodes_idx],
                  r_Idxs]; digits=4)

dx_dp_θhg_by_r  =
    round.( dx_dp[ θh_Idxs[gens_nodes_idx],
                  r_Idxs]; digits=4)

dx_dp_vhng_by_r  =
    round.( dx_dp[ vh_Idxs[non_gens_nodes_idx],
                  r_Idxs]; digits=4)

dx_dp_θhng_by_r  =
    round.( dx_dp[ θh_Idxs[non_gens_nodes_idx],
                  r_Idxs]; digits=4)

# sensitivy to lines reactance

dx_dp_vhg_by_x  =
    round.( dx_dp[ vh_Idxs[gens_nodes_idx],
                  x_Idxs]; digits=4)

dx_dp_θhg_by_x  =
    round.( dx_dp[ θh_Idxs[gens_nodes_idx],
                  x_Idxs]; digits=4)

dx_dp_vhng_by_x  =
    round.( dx_dp[ vh_Idxs[non_gens_nodes_idx],
                  x_Idxs]; digits=4)

dx_dp_θhng_by_x  =
    round.( dx_dp[ θh_Idxs[non_gens_nodes_idx],
                  x_Idxs]; digits=4)

#--------------------------------------------

dx_dp_x  =
    round.( dx_dp[:, x_Idxs]; digits=4)

dx_dp_b  =
    round.( dx_dp[:, b_Idxs]; digits=4)

dx_dp_ratio  =
    round.( dx_dp[:, ratio_Idxs]; digits=4)

dx_dp_angle  =
    round.( dx_dp[:, angle_Idxs]; digits=4)

#--------------------------------------------

# """

# edges_size = length(edges_r)

# edges_r_x_b_ratio_angle_idx =
#     get_edges_r_x_b_ratio_angle_idx(
#         edges_size)

# """

# #--------------------------------------------

# """

#  Check:

# t_Ynet_real =
#     Ynet_real_imag_flattend[ Ynet_real_Idxs ]

# t_Ynet_imag =
#     Ynet_real_imag_flattend[ Ynet_imag_Idxs ]


# tt_Ynet = [ t_Ynet_real[idx] + im * t_Ynet_imag[idx]  
#           for idx in Ynet_rows_Idxs_in_flattend ]

# [a .- b for (a,b) in zip(tt_Ynet, Ynet)]

# """

#--------------------------------------------

# df_dx = ForwardDiff.jacobian(
#     ( vh_θh_x ) ->
#         get_total_P_network_loss(
#             vh_θh_x, Ynet;
#             dyn_pf_flat_vh_flat_θh_Idx,
#             nodes_idx_with_adjacent_nodes_idx,
#             n2s_all_nodes_idx,
#             all_nodes_idx ),
#     vh_θh )


# df_dp = ForwardDiff.jacobian(
#     ( p ) ->
#         get_total_P_network_loss(
#             vh_θh, p;
#             dyn_pf_flat_vh_flat_θh_Idx,
#             nodes_idx_with_adjacent_nodes_idx,
#             n2s_all_nodes_idx,
#             all_nodes_idx ),
#     Ynet )


# intg_dx_dp  = -(svd( intg_df_dx )) \ intg_df_dp


#--------------------------------------------
#--------------------------------------------
