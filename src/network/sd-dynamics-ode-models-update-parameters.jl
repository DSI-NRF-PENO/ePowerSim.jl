# (c) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.4:  pg 135 - 136, 6.121 - 6.123


#-----------------------------------------------------
#-----------------------------------------------------
#  update functions
#-----------------------------------------------------
#-----------------------------------------------------

function update_im_dynamic_id_iq_pg_vh!(
    gens_dynamic_id_iq_pg_vh_by_vhθh,
    stateDiffCache,
    im_para_aux_inputs )

    (;im_nodes_voltage_Idx,
     im_vars_Idx_in_state,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     im_vars_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view) =
         im_para_aux_inputs

    (;gens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges)  =
        im_nodes_voltage_Idx

    im_vars_view_in_state[:] .=
        get_gen_nodes_im_vars_from_state(
            stateDiffCache,
            im_vars_Idx_in_state )

    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        im_vars_view_in_state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    # gen_nodes_δ_ω_ed_dash_eq_dash_views .= industrial_model_pure_states_view_in_state[ nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state ] 
    
    gens_vh_θh_view[:] .= get_gens_vh_θh( nodes_pf_U_view, gens_idx )
    
    gens_dynamic_id_iq_pg_vh_by_vhθh .=
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )
    
    return nothing
    
end



function update_dynamic_id_iq_pg_vh!(
    gens_dynamic_id_iq_pg_vh_by_vhθh,
    stateDiffCache,
    industrial_model_para_aux_inputs )

    (;industrial_model_nodes_voltage_Idx,
     industrial_model_pure_states_Idx,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     industrial_model_pure_states_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view) =
         industrial_model_para_aux_inputs

    (;gens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges)  =
        industrial_model_nodes_voltage_Idx
        
    industrial_model_pure_states_view_in_state[:] .=
        get_industrial_gen_nodes_pure_states_from_state(
            stateDiffCache,
            industrial_model_pure_states_Idx )

    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        industrial_model_pure_states_view_in_state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    # gen_nodes_δ_ω_ed_dash_eq_dash_views .= industrial_model_pure_states_view_in_state[ nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state ] 
    
    gens_vh_θh_view[:] .=
        get_gens_vh_θh( nodes_pf_U_view, gens_idx )
    
    gens_dynamic_id_iq_pg_vh_by_vhθh .=
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )
    
    return nothing
    
end

function update_dynamic_τm_vf!(
    vec_τm_vf, stateDiffCache,
    industrial_model_para_aux_inputs )

    (;industrial_model_nodes_voltage_Idx,
     industrial_model_pure_states_Idx,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     industrial_model_pure_states_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view) =
         industrial_model_para_aux_inputs

    (;lgens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges)  =
         industrial_model_nodes_voltage_Idx
    
    # gen_nodes_state_views .= get_gen_nodes_state_views(stateDiffCache, gens_nodes_collection )
        
    industrial_model_pure_states_view_in_state[:] .=
        get_industrial_gen_nodes_pure_states_from_state(
            stateDiffCache,
            industrial_model_pure_states_Idx )

    # update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!( gen_nodes_δ_ω_ed_dash_eq_dash_views, industrial_model_pure_states_view_in_state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    gen_nodes_δ_ω_ed_dash_eq_dash_views[:] .=
        industrial_model_pure_states_view_in_state[
            nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state ] 
    
    gens_vh_θh_view[:] .=
        get_gens_vh_θh( nodes_pf_U_view, gens_idx )

    # nodes_u_Idx_in_ranges
    
    vec_τm_vf .=
        get_gens_dynamic_τm_vf(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view  )
    
    return nothing
    
end

                         

#-----------------------------------------------------
#-----------------------------------------------------
#  Intitialisations
#-----------------------------------------------------
#-----------------------------------------------------


"""


id_iq =
    invZ_dq(ra, X_d_dash, X_q_dash) *
    [ed_dash - vh * sin(δ - θh),
     eq_dash - vh * cos(δ - θh)]

pg =
    ed_dash * id_iq[1] + eq_dash * id_iq[2] +
    ( X_q_dash - X_d_dash ) *  *(id_iq...)

id_iq =
    get_dynamic_idq_vhθh(
        vh, θh, δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )

pg =
    get_dynamic_pg_from_id_iq(
        id_iq..., δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)


"""



#---------------------------------------------------
# Industrial
#---------------------------------------------------


function get_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init(
    netd
    ; maxiter=40,
    ftol=1000*eps(),
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false )
    
    # -----------------------------------------------
    # -----------------------------------------------

    dict_sys_to_industry =
        get_net_to_industrial_model_indices_dict(
            netd )

    (; pure_states_Idx_in_system,
     ur_ui_Idx_in_system,
     industrial_Idx,
     industrial_model_pure_states_Idx,
     industrial_model_ur_ui_Idx,
     pure_states_and_ur_ui_Idx_in_system,
     industrial_model_state_Idx,
     net_to_industrial_idx_conversion_dict) =
        get_industrial_model_indices_and_conversion_dict(
            netd  )

    # ---------------------------------------------------
        
    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_δ_ω_ed_dash_eq_dash_Idx(
            netd )

    nodes_ω_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_ω_ed_dash_eq_dash_Idx(
            netd )

    nodes_δ_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_δ_ed_dash_eq_dash_Idx(
            netd )

    #------------------------------------------------- 
    
    industrial_model_each_gen_nodes_pure_states_idx_in_state =
        get_industrial_gens_pure_states_indices_in_state(
            netd )
    
    industrial_model_each_gen_nodes_stab_states_idx_in_state =
        get_industrial_gens_stab_states_indices_in_state(
            netd )
    
    
    # ---------------------------------------------------

    nodes_cb_sw =
        get_nodes_cb_sw( netd.nodes )

    gens_nodes =
        get_gens_nodes( netd.nodes )

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #---------------------------------------------------

    gen_nodes_ra_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list =
                [ :ra, :X_d_dash, :X_q_dash ] )

    gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list =
                [ :ra, :X_d, :X_q, :X_d_dash, :X_q_dash ] )

    gen_nodes_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes; param_list =
                [ :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q,
                  :X_d_dash, :X_q_dash,
                  :T_d_dash, :T_q_dash ],
            gens_view_only = true )
    
    gen_nodes_sub_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes; param_list = [
                :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q,
                :X_d_dash,
                :X_q_dash, :X_d_2dash,
                :X_q_2dash,
                :T_d_dash, :T_q_dash,
                :T_d_2dash, :T_q_2dash ],
            gens_view_only = true )
   
    # -----------------------------------------------------
    # -----------------------------------------------------

    # pf_net_param = get_powerflow_net_parameters( netd )

    pf_net_param =
        get_industrial_model_powerflow_net_parameters(
            netd )
    
    (; pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) = pf_net_param

    (;Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net

    (; slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (; ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
        pf_param_views

    load_trans_nodes_Idx_and_vlimits =
        pf_limits

    (;slack_bus_idx,
     nodes_u_Idx,
     gens_idx,
     ur_IDX,
     ui_IDX,
     vh_IDX,
     θh_IDX,
     red_vh_θh_idx,
     ur_idx,
     ui_idx,
     ur_ui_idx) =
         pf_Idx


    #--------------------------------------------------    
    #--------------------------------------------------

    state = zeros(
        length(
            generate_industrial_model_sym(
                ; nodes = netd.nodes ) ) )

    #-------------------------------------------------   
    
    gen_nodes_δ_ω_ed_dash_eq_dash_views =
        get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )
    
    industrial_model_pure_states_view_in_state =
        get_industrial_pure_states_view_in_state(
            state, industrial_model_pure_states_Idx )
    
    #-------------------------------------------------

    state_view = view(state, 1:length(state))

    # ----------------------------------------------
    
    # for pf flat start, ur = vh = 1,  ui = θh = 0

    # remove init_pf = true
        
    
    if init_pf == true 

        state_view[ur_idx] .=
            ones(  length( ur_IDX ))
        state_view[ui_idx] .=
            zeros( length( ui_IDX ))

    end
    
    # ----------------------------------------------

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]
        
    #-----------------------------------------------
    
    nodes_u_view  =
        [ view(state, nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ]

    # ---------------------------------------------
    
    nodes_pf_U_view  =
        [ view(pf_state , nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ] 
    
    # ---------------------------------------------
    
    uh_state =
        state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]
    
    x0_ur_ui =
        [state_view[ur_ui_idx][ ur_IDX ]...;
         state_view[ur_ui_idx][ ui_IDX ]...]
        
    x0_vh_θh =
        [abs.(uh_state)...; angle.(uh_state)...]

    working_vh_θh_view =
        view(x0_vh_θh, 1:length(x0_vh_θh))

    x0_vh_view       = @view x0_vh_θh[vh_IDX]

    x0_θh_view       = @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0      = [ red_vh_θh_0_view; ]

    mismatch         = similar( red_vh_θh_0 )
    
    # ----------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))
    
    # ---------------------------------------------
    
    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet  =
        zeros(ComplexF64, length( uh_state ))
     
    Inet_view  =
        view( Inet, 1:length( Inet ) )

    Iinj  =
        zeros(ComplexF64, length( uh_state ))

    Iinj_view =
        view(Iinj, 1:length( Iinj ))

    idq_wt_pad =
        zeros(ComplexF64, length( uh_state ))

    idq_wt_pad_view =
        view(idq_wt_pad, 1:length( uh_state ) )
    
    #----------------------------------------------

    global_pf_views =
        (; working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view )

    sd_pf_views =
        (; working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          gen_nodes_δ_ω_ed_dash_eq_dash_views )

    # ----------------------------------------------
    # global_pf param
    # ----------------------------------------------
    
    global_pf_param =
        (; pf_net_param,
          sd_pf_views,
          mismatch )
    
    branches_name  = collect(keys( netd.edges ))
    
    nodes_name = collect(keys( netd.nodes ))    

    #-----------------------------------------------
    #-----------------------------------------------

    # maxiter=40
    # ftol=1000*eps()
    # xtol=1000*eps()
    # init_pf = true
    # with_δ_ed_eq = false
    
    named_tup_pf_result =
        power_balance_powerflow(
            x0_vh_θh,
            mismatch,
            sd_pf_views,
            (nodes_name, branches_name) ,
            pf_net_param
            ; maxiter=maxiter,
            ftol=ftol ,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )

    #-----------------------------------------------

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init =
        named_tup_pf_result.branch_dict_init
    
    #---------------------------------------------    
    #---------------------------------------------

    state .=
        industrial_model_init_operationpoint(
            netd, bus_dict_init; pure = :pure )

    #--------------------------------------------- 

    # gen_nodes_δ_ω_ed_dash_eq_dash_views .= get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash( state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    
    #--------------------------------------------
    
    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)
    
    #--------------------------------------------

    # nodes_f_t  = external_get_nodes_or_edges_f_t(
    #     netd.nodes, bus_dict_init )

    #--------------------------------------------
    #-------------------------------------------
 
    # gens_vh_θh_post_pf =
    #     get_gens_vh_θh_post_pf(
    #         gens_nodes_collection ,
    #         bus_dict_init )
        
    # gens_ur_ui_post_pf =
    #     get_gens_ur_ui_post_pf(
    #         gens_nodes_collection,
    #         bus_dict_init )
    
    gens_vh_θh =
        get_gens_vh_θh(
            nodes_pf_U_view, gens_idx )

    gens_vh_θh_view =
        @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    """
    gens_ur_ui = get_gens_ur_ui(nodes_pf_U_view, gens_idx )

    gens_ur_ui_view = @view gens_ur_ui[ 1:length(gens_ur_ui ) ]

    """

    #---------------------------------------------
    #---------------------------------------------
    
    idq_wt_pad_view[gens_idx] .=
        [ industrial_model_get_dynamic_idq_θ_π_vhθh(
            vh_θh, δ_ω_ed_eq,
            ra_Xd_dash_Xq_dash )
          for (vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash) in
              zip( gens_vh_θh_view ,
                   gen_nodes_δ_ω_ed_dash_eq_dash_views,
                   ra_Xd_dash_Xq_dash_view[gens_idx] ) ]
  

    #---------------------------------------------

    # gens_nodes_ωs_τm_v_ref =
    #     get_gens_nodes_ωs_τm_v_ref(
    #         gens_nodes_collection,
    #         bus_dict_init )

    # gens_nodes_ωs_τm_v_ref_view =
    #     view( gens_nodes_ωs_τm_v_ref,
    #           1:length( gens_nodes_ωs_τm_v_ref ) )

    #----------------------------------------------

    gens_nodes_ωs_τm_vref_porder =
        get_gens_nodes_ωs_τm_vref_porder(
            gens_nodes_collection, bus_dict_init )

    gens_nodes_ωs_τm_vref_porder_view =
        view(
            gens_nodes_ωs_τm_vref_porder,
            1:length(gens_nodes_ωs_τm_vref_porder) )
    
    #-----------------------------------------------

    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )

    
    #-----------------------------------------------
    #-----------------------------------------------

    gens_dynamic_id_iq_pg_vh_by_vhθh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )

    #-----------------------------------------------  

    nodes_u_Idx_in_ranges =
        get_nodes_u_Idx_in_ranges(
            nodes_u_Idx )    
    #-----------------------------------------------
    
    non_gens_idx =
        get_load_trans_nodes_Idx(
            netd.nodes )
    
    #------------------------------------------------

    industrial_model_nodes_voltage_Idx =
        (; gens_idx,
         non_gens_idx,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges )
    
    industrial_model_misc_Idx =
        (; nodes_u_Idx,
         nodes_u_Idx_in_ranges,
         industrial_model_pure_states_Idx,
         industrial_model_each_gen_nodes_pure_states_idx_in_state,
         gens_nodes_collection )

    para_update_gen_Ax_aux =
        (; industrial_model_pure_states_view_in_state,
         industrial_model_each_gen_nodes_pure_states_idx_in_state,
         industrial_model_pure_states_Idx )

    industrial_model_para_aux_inputs =
        (; industrial_model_nodes_voltage_Idx,
         industrial_model_pure_states_Idx,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         gen_nodes_ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
         industrial_model_pure_states_view_in_state,
         gen_nodes_δ_ω_ed_dash_eq_dash_views,
         gens_vh_θh_view, nodes_pf_U_view )
    

    industrial_model_pf_para =
        (; gens_dynamic_id_iq_pg_vh_by_vhθh,
         gens_nodes_ωs_τm_vref_porder_view )

    
    """ need by their views """

    industrial_model_ωs_τm_vref_vhθh_idq =
        (; gens_nodes_ωs_τm_vref_porder,
         gens_vh_θh, idq_wt_pad )
    

    industrial_model_dyn_pf_up_para =
        (; nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_δ_ed_dash_eq_dash_Idxs_in_state )
    
    #----------------------------------------------

    industrial_model_idq_pf_cal =
        (;  idq_wt_pad_view,
           gens_idx )
    
    #---------------------------------------------

    return (; nodes_cb_sw, state, global_pf_param,
            named_tup_pf_result,
            industrial_model_misc_Idx,
            para_update_gen_Ax_aux,
            industrial_model_para_aux_inputs,
            industrial_model_pf_para,
            industrial_model_ωs_τm_vref_vhθh_idq ,
            industrial_model_dyn_pf_up_para,
            industrial_model_idq_pf_cal,
            gens_nodes_τm_vf  )
    
end


function get_industrial_model_pf_param_views_and_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen  =false )


    dict_sys_to_industry =
        get_net_to_industrial_model_indices_dict(
            netd; no_control_device = only_gen  )

    (;pure_states_Idx_in_system, ur_ui_Idx_in_system,
     industrial_Idx, industrial_model_pure_states_Idx,
     industrial_model_ur_ui_Idx,
     pure_states_and_ur_ui_Idx_in_system,
     industrial_model_state_Idx,
     net_to_industrial_idx_conversion_dict) =
        get_industrial_model_indices_and_conversion_dict(
            netd; no_control_device = only_gen  )

    # -----------------------------------------------
    # -----------------------------------------------  

    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_δ_ω_ed_dash_eq_dash_Idx(
            netd; no_control_device = only_gen )

    nodes_ω_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_ω_ed_dash_eq_dash_Idx(
            netd; no_control_device = only_gen )

    nodes_δ_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_δ_ed_dash_eq_dash_Idx(
            netd; no_control_device = only_gen )

    #-----------------------------------------------

    industrial_model_each_gen_nodes_pure_states_idx_in_state =
        get_industrial_gens_pure_states_indices_in_state(
            netd; no_control_device = only_gen )

    industrial_model_each_gen_nodes_stab_states_idx_in_state =
        get_industrial_gens_stab_states_indices_in_state(
            netd; no_control_device = only_gen )

    #----------------------------------------------

    gen_nodes_ra_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d_dash, :X_q_dash ] )

    gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list =
                [ :ra, :X_d, :X_q,
                  :X_d_dash, :X_q_dash ] )

    gen_nodes_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list =
                [ :D, :H, :ωs, :ra, :xℓ, :X_d,
                  :X_q, :X_d_dash, :X_q_dash,
                  :T_d_dash, :T_q_dash ],
            gens_view_only = only_gen )

    gen_nodes_sub_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list =
                [ :D, :H, :ωs, :ra, :xℓ, :X_d,
                  :X_q, :X_d_dash, :X_q_dash,
                  :X_d_2dash, :X_q_2dash, :T_d_dash,
                  :T_q_dash, :T_d_2dash, :T_q_2dash ],
            gens_view_only = only_gen )

    # -----------------------------------------------
    # -----------------------------------------------

    # pf_net_param = get_powerflow_net_parameters( netd )

    pf_net_param = get_industrial_model_powerflow_net_parameters(
        netd; no_control_device = only_gen )

    (; pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits, pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) = pf_net_param

    (; Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net

    (;slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
        pf_idx_and_state

    (;ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits =
        pf_limits

    (;slack_bus_idx,
     nodes_u_Idx,
     gens_idx,
     ur_IDX,
     ui_IDX,
     vh_IDX,
     θh_IDX,
     red_vh_θh_idx,
     ur_idx,
     ui_idx,
     ur_ui_idx) =
         pf_Idx


    #---------------------------------------------    
    #---------------------------------------------

    state = zeros(
        length( generate_industrial_model_sym(
            ; nodes = netd.nodes,
            no_control_device = only_gen ) ) )

    #---------------------------------------------   

    gen_nodes_δ_ω_ed_dash_eq_dash_views =
        get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            state,
            nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    industrial_model_pure_states_view_in_state =
        get_industrial_pure_states_view_in_state(
            state, industrial_model_pure_states_Idx )

    #---------------------------------------------------- 
    #----------------------------------------------------  

    state_view = view(state, 1:length(state))

    # ----------------------------------------------------

    # for pf flat start, ur = vh = 1,  ui = θh = 0

    # remove init_pf = true


    if init_pf == true 

        state_view[ur_idx] .= ones(  length( ur_IDX ))
        state_view[ui_idx] .= zeros( length( ui_IDX ))

    end

    # ----------------------------------------------------

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]

    #----------------------------------------------------

    nodes_u_view  = [
        view(state, nodes_u_Idx[Ind])
        for Ind in
            collect(1:length(nodes_u_Idx)) ]

    # ---------------------------------------------------- 

    nodes_pf_U_view  = [
        view(pf_state , nodes_u_Idx[Ind])
        for Ind in
            collect(1:length(nodes_u_Idx)) ] 

    # ----------------------------------------------------

    uh_state = state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]

    x0_ur_ui = [state_view[ur_ui_idx][ ur_IDX ]...;
                state_view[ur_ui_idx][ ui_IDX ]...]

    x0_vh_θh =
        [abs.(uh_state)...;
         angle.(uh_state)...]

    working_vh_θh_view =
        view(x0_vh_θh, 1:length(x0_vh_θh))

    x0_vh_view  = @view x0_vh_θh[vh_IDX]

    x0_θh_view  = @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view =
        @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0  = [ red_vh_θh_0_view; ]

    mismatch  = similar( red_vh_θh_0 )

    # ----------------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))

    # ----------------------------------------------------

    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet =
        zeros(ComplexF64, length( uh_state ))

    Inet_view  =
        view( Inet, 1:length( Inet ) )

    Iinj =
        zeros(ComplexF64, length( uh_state ))

    Iinj_view  =
        view(Iinj, 1:length( Iinj ))

    idq_wt_pad =
        zeros(ComplexF64, length( uh_state ))

    idq_wt_pad_view =
        view(idq_wt_pad, 1:length( uh_state ) )

    #----------------------------------------------------

    global_pf_views =
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view, Iinj_view )

    sd_pf_views =
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          gen_nodes_δ_ω_ed_dash_eq_dash_views )

    # ---------------------------------------------------
    # global_pf param
    # ---------------------------------------------------

    global_pf_param =
        ( pf_net_param,
          sd_pf_views, mismatch )

    branches_name =
        collect(keys( netd.edges ))

    nodes_name =
        collect(keys( netd.nodes ))    

    #----------------------------------------------------
    #----------------------------------------------------

    named_tup_pf_result =
        power_balance_powerflow(
            x0_vh_θh, mismatch, sd_pf_views,
            (nodes_name, branches_name) ,
            pf_net_param ;
            maxiter=maxiter,
            ftol=ftol ,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )

    #----------------------------------------------------

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init =
        named_tup_pf_result.branch_dict_init

    #---------------------------------------------
    #---------------------------------------------

    state .= industrial_model_init_operationpoint(
        netd, bus_dict_init; pure = :pure, no_control_device = only_gen )

    return (; dict_sys_to_industry, pure_states_Idx_in_system, ur_ui_Idx_in_system, industrial_Idx, industrial_model_pure_states_Idx, industrial_model_ur_ui_Idx, pure_states_and_ur_ui_Idx_in_system, industrial_model_state_Idx, net_to_industrial_idx_conversion_dict, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, nodes_ω_ed_dash_eq_dash_Idxs_in_state, nodes_δ_ed_dash_eq_dash_Idxs_in_state, industrial_model_each_gen_nodes_pure_states_idx_in_state, industrial_model_each_gen_nodes_stab_states_idx_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view, gen_nodes_dyn_param_view, gen_nodes_sub_dyn_param_view, pf_net_param, ra_Xd_dash_Xq_dash_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, industrial_model_pure_states_view_in_state, state_view, pf_state, nodes_u_view, nodes_pf_U_view, x0_vh_θh, working_vh_θh_view, red_vh_θh_0_view, mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, branches_name, nodes_name, named_tup_pf_result, bus_dict_init, branch_dict_init, state, gens_idx, nodes_u_Idx )
    
    
end



function get_industrial_model_pf_param_views_and_init_with_or_no_controller(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen =false  )

    # ----------------------------------------------
    # ----------------------------------------------

    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)

    gens_nodes = get_gens_nodes( netd.nodes  )

    non_gens_nodes = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #----------------------------------------------
    #---------------------------------------------- 

    (;dict_sys_to_industry,
     pure_states_Idx_in_system,
     ur_ui_Idx_in_system,
     industrial_Idx,
     industrial_model_pure_states_Idx,
     industrial_model_ur_ui_Idx,
     pure_states_and_ur_ui_Idx_in_system,
     industrial_model_state_Idx,
     net_to_industrial_idx_conversion_dict,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_δ_ed_dash_eq_dash_Idxs_in_state,
     industrial_model_each_gen_nodes_pure_states_idx_in_state,
     industrial_model_each_gen_nodes_stab_states_idx_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     gen_nodes_dyn_param_view,
     gen_nodes_sub_dyn_param_view,
     pf_net_param,
     ra_Xd_dash_Xq_dash_view,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     industrial_model_pure_states_view_in_state,
     state_view, pf_state, nodes_u_view,
     nodes_pf_U_view,
     x0_vh_θh,
     working_vh_θh_view,
     red_vh_θh_0_view,
     mismatch,
     Jac_vh_θh,
     Inet,
     Inet_view,
     Iinj,
     Iinj_view,
     idq_wt_pad,
     idq_wt_pad_view,
     global_pf_views,
     sd_pf_views,
     global_pf_param,
     branches_name,
     nodes_name,
     named_tup_pf_result,
     bus_dict_init,
     branch_dict_init,
     state,
     gens_idx,
     nodes_u_Idx) =
        get_industrial_model_pf_param_views_and_init(
            netd;
            maxiter=40,
            ftol=1000*eps(),
            xtol=1000*eps(),
            init_pf = true,
            with_δ_ed_eq = false,
            only_gen  =false)
    
    #---------------------------------------------------- 
    #---------------------------------------------------- 
    
    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )


    #----------------------------------------------------

    nodes_cb_sw =
        get_nodes_cb_sw(netd.nodes)

    #----------------------------------------------------

    gens_vh_θh = get_gens_vh_θh(
        nodes_pf_U_view, gens_idx )

    gens_vh_θh_view =
        @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    #---------------------------------------------------
    #---------------------------------------------------

    # gen_uh  = (named_tup_pf_result.Vbus)[ gen_idx ]

    idq_wt_pad_view[gens_idx] .=  [
        industrial_model_get_dynamic_idq_θ_π_vhθh(
            vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
        for ( vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in
            zip( gens_vh_θh_view ,
                 gen_nodes_δ_ω_ed_dash_eq_dash_views,
                 ra_Xd_dash_Xq_dash_view[gens_idx] ) ]

    #---------------------------------------------------

    gens_nodes_ωs_τm_vref_porder =
        get_gens_nodes_ωs_τm_vref_porder(
            gens_nodes_collection,
            bus_dict_init )

    gens_nodes_ωs_τm_vref_porder_view =
        view( gens_nodes_ωs_τm_vref_porder,
              1:length(gens_nodes_ωs_τm_vref_porder) )

    #---------------------------------------------------

    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )

    vec_τm_vf_views =
        view( gens_nodes_τm_vf,
              1:length(gens_nodes_τm_vf) )
    
    #---------------------------------------------------   
    #---------------------------------------------------

    gens_dynamic_id_iq_pg_vh_by_vhθh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )

    gens_dynamic_id_iq_pg_vh_by_vhθh_view =
        view( gens_dynamic_id_iq_pg_vh_by_vhθh,
              1:length(gens_dynamic_id_iq_pg_vh_by_vhθh) )

    """
    gens_dynamic_id_iq_pg_vh_by_ur_ui =
        get_gens_dynamic_id_iq_pg_vh_by_ur_ui(
            gens_ur_ui_post_pf,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )
    
    """
    
    #---------------------------------------------------  

    nodes_u_Idx_in_ranges =
        get_nodes_u_Idx_in_ranges(
            nodes_u_Idx )
    
    #---------------------------------------------------

    non_gens_idx =
        get_load_trans_nodes_Idx(
            netd.nodes )

    #---------------------------------------------------

    industrial_model_nodes_voltage_Idx =
        (; gens_idx,
         non_gens_idx,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges )

    industrial_model_misc_Idx =
        (; nodes_u_Idx,
         nodes_u_Idx_in_ranges,
         industrial_model_pure_states_Idx,
         industrial_model_each_gen_nodes_pure_states_idx_in_state,
         gens_nodes_collection )

    para_update_gen_Ax_aux =
        (; industrial_model_pure_states_view_in_state,
         industrial_model_each_gen_nodes_pure_states_idx_in_state,
         industrial_model_pure_states_Idx )

    industrial_model_para_aux_inputs =
        (; industrial_model_nodes_voltage_Idx,
         industrial_model_pure_states_Idx,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         gen_nodes_ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
         industrial_model_pure_states_view_in_state,
         gen_nodes_δ_ω_ed_dash_eq_dash_views,
         gens_vh_θh_view,
         nodes_pf_U_view )

    if only_gen == false
        
        industrial_model_pf_para =
            (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
             gens_nodes_ωs_τm_vref_porder_view )
        
    else
        industrial_model_pf_para =
            (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
             gens_nodes_ωs_τm_vref_porder_view,
             vec_τm_vf_views )
        
    end


    """ need by their views """

    industrial_model_ωs_τm_vref_vhθh_idq =
        (; gens_dynamic_id_iq_pg_vh_by_vhθh,
         gens_nodes_ωs_τm_vref_porder,
         gens_nodes_τm_vf,
         gens_vh_θh,
         idq_wt_pad )


    industrial_model_dyn_pf_up_para =
        (; nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_δ_ed_dash_eq_dash_Idxs_in_state )

    #---------------------------------------------------

    industrial_model_idq_pf_cal =
        (;  idq_wt_pad_view, gens_idx )

    #---------------------------------------------------

    return (; nodes_cb_sw,
            state,
            global_pf_param,
            named_tup_pf_result,
            industrial_model_misc_Idx,
            para_update_gen_Ax_aux,
            industrial_model_para_aux_inputs,
            industrial_model_pf_para,
            industrial_model_ωs_τm_vref_vhθh_idq,
            industrial_model_dyn_pf_up_para,
            industrial_model_idq_pf_cal )

    #---------------------------------------------------
    #---------------------------------------------------

end


#---------------------------------------------------
# v2
#---------------------------------------------------


function v2_get_im_model_pf_param_views_and_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen  =false)

    #-------------------------------
    
    # dynamics_case = case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    # dynamics_case =
    #     case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    
    # netd     = NetworkData( dynamics_case()... )
    # only_gen = true
    # maxiter  = 40
    # ftol     = 1000*eps()
    # xtol     = 1000*eps()
    # init_pf  = true
    # with_δ_ed_eq = false
    
    # #-------------------------------
    

    dict_sys_to_im = get_net_to_im_indices_dict( netd  )

    (; im_vars_indices_in_system,
     pure_states_Idx_in_system,
     im_algebraic_vars_Idx_in_system,
     ur_ui_Idx_in_system,
     im_vars_and_ur_ui_Idx_in_system,
     im_vars_Idx_in_state,
     nodes_ur_ui_Idx_in_state,
     im_state_Idx,
     each_gens_im_vars_Idx_in_state,
     net_to_im_idx_conversion_dict)  =
         get_im_indices_and_conversion_dict(
             netd  )
    
    # ---------------------------------------------------
    # ---------------------------------------------------    

    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
        get_im_δ_ω_ed_dash_eq_dash_Idx(netd )

    nodes_ω_ed_dash_eq_dash_Idxs_in_state =
        get_im_ω_ed_dash_eq_dash_Idx(netd )

    nodes_δ_ed_dash_eq_dash_Idxs_in_state =
        get_im_δ_ed_dash_eq_dash_Idx(netd )

    #---------------------------------------------------
    #---------------------------------------------------

    gen_nodes_ra_Xd_dash_Xq_dash_view =
        get_im_model_gens_params_view_in_param_values(
        netd.nodes_param,
        netd.nodes; param_list = [
            :ra, :X_d_dash, :X_q_dash ] )

    gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_im_model_gens_params_view_in_param_values(
        netd.nodes_param,
        netd.nodes; param_list = [
            :ra, :X_d, :X_q, :X_d_dash, :X_q_dash ] )

    gen_nodes_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes
            ; param_list = [
                :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q,
                :X_d_dash, :X_q_dash,
                :T_d_dash, :T_q_dash ])

    gen_nodes_sub_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes
            ; param_list = [
                :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q,
                :X_d_dash, :X_q_dash,
                :X_d_2dash, :X_q_2dash,
                :T_d_dash, :T_q_dash,
                :T_d_2dash, :T_q_2dash ] )

    # -----------------------------------------------------
    # -----------------------------------------------------

    pf_net_param =
        get_im_model_powerflow_net_parameters( netd  )

    (;pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits, pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) =
         pf_net_param

    (;Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
        pf_net

    (;slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
        pf_idx_and_state

    (; ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
        pf_param_views

    load_trans_nodes_Idx_and_vlimits =
        pf_limits

    (;slack_bus_idx,
     nodes_u_Idx,
     gens_idx,
     ur_IDX,
     ui_IDX,
     vh_IDX,
     θh_IDX,
     red_vh_θh_idx,
     ur_idx,
     ui_idx,
     ur_ui_idx) = pf_Idx


    #----------------------------------------------------    
    #----------------------------------------------------

    state = zeros(length(
        generate_im_sym(
            ; nodes = netd.nodes )  ) )

    #----------------------------------------------------   

    gen_nodes_δ_ω_ed_dash_eq_dash_views =
        get_im_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            state,
            nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    # im_model_pure_states_view_in_state =
    #     get_im_pure_states_view_in_state(
    #         state,
    #         pure_states_Idx )

    im_vars_view_in_state =
        get_im_vars_view_in_state( state, im_vars_Idx_in_state )

    #----------------------------------------------------  

    state_view = view(state, 1:length(state))

    # ----------------------------------------------------

    # for pf flat start, ur = vh = 1,  ui = θh = 0

    if init_pf == true 

        state_view[ur_idx] .= ones(  length( ur_IDX ))
        state_view[ui_idx] .= zeros( length( ui_IDX ))

    end

    # ----------------------------------------------------

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]

    #----------------------------------------------------

    nodes_u_view  = [
        view(state, nodes_u_Idx[Ind])
        for Ind in collect(1:length(nodes_u_Idx)) ]

    # ---------------------------------------------------- 

    nodes_pf_U_view  = [
        view(pf_state , nodes_u_Idx[Ind])
        for Ind in collect(1:length(nodes_u_Idx)) ] 

    # ----------------------------------------------------

    uh_state = state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]

    x0_ur_ui = [
        state_view[ur_ui_idx][ ur_IDX ]...;
        state_view[ur_ui_idx][ ui_IDX ]...]

    x0_vh_θh      = [abs.(uh_state)...; angle.(uh_state)...]

    working_vh_θh_view = view(x0_vh_θh, 1:length(x0_vh_θh))

    x0_vh_view       = @view x0_vh_θh[vh_IDX]

    x0_θh_view       = @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0      = [ red_vh_θh_0_view; ]

    mismatch         = similar( red_vh_θh_0 )

    # ----------------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))

    # ----------------------------------------------------

    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet            = zeros(ComplexF64, length( uh_state ))

    Inet_view       =  view( Inet, 1:length( Inet ) )

    Iinj            = zeros(ComplexF64, length( uh_state ))

    Iinj_view       =  view(Iinj, 1:length( Iinj ))

    idq_wt_pad      = zeros(ComplexF64, length( uh_state ))

    idq_wt_pad_view = view(idq_wt_pad,1:length(uh_state) )


    # idq_wt_pad =
    #     [[0.0, 0.0]
    #      for idx in 1:length(uh_state )]

    # idq_wt_pad_view =
    #     view(idq_wt_pad, 1:length( uh_state ) )
    
    #----------------------------------------------------

    global_pf_views = (
        working_vh_θh_view,
        nodes_pf_U_view,
        Inet_view,
        Iinj_view )

    sd_pf_views = (
        working_vh_θh_view,
        nodes_pf_U_view,
        Inet_view,
        Iinj_view,
        gen_nodes_δ_ω_ed_dash_eq_dash_views )

    # ---------------------------------------------------
    # global_pf param
    # ---------------------------------------------------

    global_pf_param = (
        pf_net_param,
        sd_pf_views,
        mismatch )

    branches_name  = collect(keys( netd.edges ))

    nodes_name     = collect(keys( netd.nodes ))    

    #----------------------------------------------------

    named_tup_pf_result = power_balance_powerflow(
        x0_vh_θh,
        mismatch,
        sd_pf_views,
        (nodes_name, branches_name) ,
        pf_net_param;
        maxiter=maxiter,
        ftol=ftol,
        xtol=xtol,
        with_δ_ed_eq = with_δ_ed_eq )

    #----------------------------------------------------

    bus_dict_init = named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init

    #---------------------------------------------------

    state .= im_model_init_operationpoint(
        netd, bus_dict_init  )

    return (
        ; im_vars_indices_in_system,
        pure_states_Idx_in_system,
        im_algebraic_vars_Idx_in_system,
        ur_ui_Idx_in_system,
        im_vars_and_ur_ui_Idx_in_system,
        im_vars_Idx_in_state,
        nodes_ur_ui_Idx_in_state,
        im_state_Idx,
        each_gens_im_vars_Idx_in_state,
        net_to_im_idx_conversion_dict,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_δ_ed_dash_eq_dash_Idxs_in_state,
        gen_nodes_ra_Xd_dash_Xq_dash_view,
        gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
        gen_nodes_dyn_param_view,
        gen_nodes_sub_dyn_param_view,
        pf_net_param,
        ra_Xd_dash_Xq_dash_view,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        im_vars_view_in_state,
        state_view,
        pf_state,
        nodes_u_view,
        nodes_pf_U_view,
        x0_vh_θh,
        working_vh_θh_view,
        red_vh_θh_0_view,
        mismatch,
        Jac_vh_θh,
        Inet,
        Inet_view,
        Iinj,
        Iinj_view,
        idq_wt_pad,
        idq_wt_pad_view,
        global_pf_views,
        sd_pf_views,
        global_pf_param,
        branches_name,
        nodes_name,
        named_tup_pf_result,
        bus_dict_init,
        branch_dict_init,
        state,
        gens_idx,
        nodes_u_Idx )
    
    
end



function v2_get_im_sys_pf_param_views_and_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen =false  )

    # ---------------------------------------------------

    nodes_cb_sw           = get_nodes_cb_sw(netd.nodes)

    gens_nodes            = get_gens_nodes( netd.nodes  )

    non_gens_nodes        = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #---------------------------------------------------

    (; im_vars_indices_in_system,
     pure_states_Idx_in_system,
     im_algebraic_vars_Idx_in_system,
     ur_ui_Idx_in_system,
     im_vars_and_ur_ui_Idx_in_system,
     im_vars_Idx_in_state,
     nodes_ur_ui_Idx_in_state,
     im_state_Idx,
     each_gens_im_vars_Idx_in_state,
     net_to_im_idx_conversion_dict,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_δ_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     gen_nodes_dyn_param_view,
     gen_nodes_sub_dyn_param_view,
     pf_net_param,
     ra_Xd_dash_Xq_dash_view,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     im_vars_view_in_state,
     state_view,
     pf_state,
     nodes_u_view,
     nodes_pf_U_view,
     x0_vh_θh,
     working_vh_θh_view,
     red_vh_θh_0_view,
     mismatch,
     Jac_vh_θh,
     Inet,
     Inet_view,
     Iinj, Iinj_view,
     idq_wt_pad,
     idq_wt_pad_view,
     global_pf_views,
     sd_pf_views,
     global_pf_param,
     branches_name,
     nodes_name,
     named_tup_pf_result,
     bus_dict_init,
     branch_dict_init,
     state,
     gens_idx,
     nodes_u_Idx)  =
        v2_get_im_model_pf_param_views_and_init(
            netd;
            maxiter=40,
            ftol=1000*eps() ,
            xtol=1000*eps(),
            init_pf = true,
            with_δ_ed_eq = false,
            only_gen = only_gen )
    
    #---------------------------------------------------- 
    
    update_im_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )


    #----------------------------------------------------

    nodes_cb_sw =
        get_nodes_cb_sw(netd.nodes)

    #----------------------------------------------------

    gens_vh_θh =
        get_gens_vh_θh( nodes_pf_U_view, gens_idx )

    gens_vh_θh_view =
        @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    #---------------------------------------------------
    
    
    dyn_idq = [ get_pf_dynamic_idq_θ_π_vhθh(
            vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
                for ( vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in
                    zip(
                        gens_vh_θh_view ,
                        gen_nodes_δ_ω_ed_dash_eq_dash_views,
                        ra_Xd_dash_Xq_dash_view[gens_idx] ) ]

    idq_wt_pad_view[gens_idx] .= dyn_idq
    
    # for idx in gens_idx
        
    #     idq_wt_pad_view[idx] = dyn_idq[idx]
        
    # end
    
    #---------------------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        get_gens_nodes_ωs_ωref0_vref0_porder0(
        gens_nodes_collection,
        bus_dict_init )

    gens_nodes_ωs_ωref0_vref0_porder0_view = view(
        gens_nodes_ωs_ωref0_vref0_porder0,
        1:length( gens_nodes_ωs_ωref0_vref0_porder0 ) )

    #---------------------------------------------------

    gens_nodes_τm_vf =
        get_gens_τm_vf( gens_nodes_collection,
                        bus_dict_init )

    vec_τm_vf_views =
        view( gens_nodes_τm_vf,
              1:length(gens_nodes_τm_vf) )
    
    #---------------------------------------------------   

    gens_dynamic_id_iq_pg_vh_by_vhθh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
        gens_vh_θh_view,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        gen_nodes_ra_Xd_dash_Xq_dash_view )

    gens_dynamic_id_iq_pg_vh_by_vhθh_view =
        view( gens_dynamic_id_iq_pg_vh_by_vhθh,
            1:length(gens_dynamic_id_iq_pg_vh_by_vhθh) )
    
    #---------------------------------------------------  

    nodes_u_Idx_in_ranges = get_nodes_u_Idx_in_ranges(
        nodes_u_Idx )
    
    #---------------------------------------------------

    non_gens_idx = get_load_trans_nodes_Idx( netd.nodes )

    #---------------------------------------------------

    
    im_nodes_voltage_Idx = (
        ; gens_idx,
        non_gens_idx,
        nodes_u_Idx,
        nodes_u_Idx_in_ranges )

    
    im_misc_Idx = (
        ; nodes_u_Idx,
        nodes_u_Idx_in_ranges,
        im_vars_Idx_in_state, 
        each_gens_im_vars_Idx_in_state, 
        gens_nodes_collection )

    
    para_update_gen_Ax_aux = (
        ; im_vars_view_in_state, 
        each_gens_im_vars_Idx_in_state, 
        im_vars_Idx_in_state ) 

    
    im_para_aux_inputs = (
        ; im_nodes_voltage_Idx,
        im_vars_Idx_in_state, 
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
        gen_nodes_ra_Xd_dash_Xq_dash_view,
        gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
        im_vars_view_in_state,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        gens_vh_θh_view,
        nodes_pf_U_view )

    
    im_pf_para = (
        ; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
        gens_nodes_ωs_ωref0_vref0_porder0_view  )


    """ need by their views """

    im_ωs_ωref_vref_vhθh_idq_etc =  (
        ; gens_dynamic_id_iq_pg_vh_by_vhθh,
        gens_nodes_ωs_ωref0_vref0_porder0,
        gens_nodes_τm_vf,
        gens_vh_θh,
        idq_wt_pad )

    im_dyn_pf_up_para = (
        ; nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_δ_ed_dash_eq_dash_Idxs_in_state )

    #---------------------------------------------------

    im_idq_pf_cal = (
        ; idq_wt_pad_view,
        gens_idx )

    #---------------------------------------------------

    return (
        ;nodes_cb_sw,
        state,
        global_pf_param,
        named_tup_pf_result,
        im_misc_Idx, 
        para_update_gen_Ax_aux, 
        im_para_aux_inputs, 
        im_pf_para, 
        im_ωs_ωref_vref_vhθh_idq_etc, 
        im_dyn_pf_up_para, 
        im_idq_pf_cal ) 
    
end

