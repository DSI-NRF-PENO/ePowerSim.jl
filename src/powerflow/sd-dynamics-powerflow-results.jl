# (c) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.4:  pg 135 - 136, 6.121 - 6.123


########################################################
# ------------------------------------------------------
#  Powerflow results
# ------------------------------------------------------
########################################################


function extract_sol_results( sol, netd )

    # extract_sol_results( netd, sol )

    # -------------------------------------------------
    # -------------------------------------------------
        
    # dynamics_case = case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix
    
    # netd = NetworkData( dynamics_case()... )
    
    # -------------------------------------------------
    # -------------------------------------------------

    """
    Idx_converstion_fun:

       `net_to_im_model_indices(idx, dict_conv)`

        `net_to_industrial_model_indices(idx, dict_conv)`

    dict_conv:

        dict_sys_to_im =
            get_net_to_im_indices_dict( netd  )

        dict_sys_to_industry =
            get_net_to_industrial_model_indices_dict(
                 netd; no_control_device = false   )

    get_a_model_streamlined_powerflow_net_parameters(
        netd;
        dict_sys_to_model_Idx = dict_sys_to_model_Idx,
        Idx_converstion_fun = Idx_converstion_fun )
    

    """
    
    # -------------------------------------------------

    pf_net_param =
        get_streamlined_powerflow_net_parameters(
            netd )
    
    #----------------------------------------------------
    
    (pf_net, pf_idx_and_state, pf_param_views,
     pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx ,
     pf_net_misc) = pf_net_param

    (Ybus, Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal, edges_orientation) = pf_net

    (slack_vh, gens_vh, gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) = pf_idx_and_state

    (ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view,
     P_Q_gens_view, P_Q_gens_loc_load_view,
     P_Q_non_gens_view) = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx, nodes_u_Idx, gens_idx,
     ur_IDX, ui_IDX, vh_IDX, θh_IDX,
     red_vh_θh_idx, ur_idx,
     ui_idx, ur_ui_idx) = pf_Idx

    # ----------------------------------------------------

    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist

    # ----------------------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx

    slack_gens_nodes_idx =
        net_comp_type_idx.slack_gens_nodes_idx
    
    non_slack_gens_nodes_idx =
        net_comp_type_idx.non_slack_gens_nodes_idx
    
    gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    gens_with_loc_load_idx =
        loc_load_exist == false ?
        [] : 
        net_comp_type_idx.gens_with_loc_load_idx

    all_nodes_idx =
        net_comp_type_idx.all_nodes_idx
    
    # ---------------------------------------------------
    
    idx_and_Idx =
        pf_and_dyn_idx_and_Idx.hybrid_pf_etc_idx_and_Idx

    red_vh_Idxs =
        idx_and_Idx.red_vh_Idxs

    red_θh_Idxs =
        idx_and_Idx.red_θh_Idxs

    non_slack_gens_θh_idx2Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx
    
    non_slack_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx_in_Idx
    
    non_gens_θh_idx2Idx =
        idx_and_Idx.non_gens_θh_idx2Idx
    
    non_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_gens_θh_idx2Idx_in_Idx
        
    
    # -------------------------------------------------- 

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx
    
    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash
    
    #----------------------------------------------------
    
    models_networks_labels =
        get_models_networks_labels(netd )
    
    #----------------------------------------------------

    network_vars_labels_hybrid =
        models_networks_labels.network_vars_labels_hybrid
    
    network_vars_labels_industrial =
        models_networks_labels.network_vars_labels_industrial
    
    network_vars_labels_im =
        models_networks_labels.network_vars_labels_im

   #----------------------------------------------------

    syms_hybrid =
        network_vars_labels_hybrid

    syms_industrial =
        network_vars_labels_industrial

    syms_im =
        network_vars_labels_im
            
    #----------------------------------------------------

    state =
        zeros(length( syms_hybrid ))

    #----------------------------------------------------
    #----------------------------------------------------
    
    """
    get_some_states_Idxs_of_a_node_in_net_state(
        ; network_vars_labels = network_vars_labels,
        bus_name = bus_name,
        some_states_syms =
            [:δ, :ed_dash, :eq_dash]  )
    
    get_gens_nodes_some_state_vars_Idxs_in_net_state(
        netd.nodes;
        network_vars_labels =
            network_vars_labels_hybrid,
        some_state_vars =
            [:δ, :ω, :ed_dash, :eq_dash],
        no_control_device = false )

    """
    
    #----------------------------------------------------

    models_gens_nodes_some_vars_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )       

    nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid =
        models_gens_nodes_some_vars_Idxs.Idxs_hybrid

    """
    nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial =
        models_gens_nodes_some_vars_Idxs.Idxs_industrial
    
    nodes_δ_ω_ed_dash_eq_dash_Idxs_im =
        models_gens_nodes_some_vars_Idxs.Idxs_im


    """
    
    #----------------------------------------------------

    """

    industrial_δ_ω_ed_dash_eq_dash_view =
        get_gen_nodes_δ_ω_ed_dash_eq_dash_views(
           state, nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial)

    industrial_δ_ω_ed_dash_eq_dash =
        get_gen_nodes_δ_ω_ed_dash_eq_dash(
           state, nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial)


    im_δ_ω_ed_dash_eq_dash_view =
        get_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            state, nodes_δ_ω_ed_dash_eq_dash_Idxs_im )

    im_δ_ω_ed_dash_eq_dash =
        get_gen_nodes_δ_ω_ed_dash_eq_dash(
            state, nodes_δ_ω_ed_dash_eq_dash_Idxs_im )
    
    """

    hybrid_δ_ω_ed_dash_eq_dash_view =
        get_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            state, nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid )

    hybrid_δ_ω_ed_dash_eq_dash =
        get_gen_nodes_δ_ω_ed_dash_eq_dash(
            state, nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid )
    
    # ---------------------------------------------------
    # ---------------------------------------------------
    
    δ_ω_ed_dash_eq_dash_view =
        hybrid_δ_ω_ed_dash_eq_dash_view
    
    # ---------------------------------------------------

    state_view =
        view(state, 1:length(state))
    
    # ----------------------------------------------------
    # x  = similar( state_x0 )
    #
    # Decoupling pf_state from state_x0
    # ----------------------------------------------------  

    pf_state = [ state; ]
        
    #----------------------------------------------------
    
    nodes_u_view  =
        [ view(state, nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ]

    # ----------------------------------------------------  
    
    nodes_pf_U_view  =
        [ view(pf_state , nodes_u_Idx[Ind])
          for Ind in collect(1:length(nodes_u_Idx)) ] 
    
    # ----------------------------------------------------
    
    uh_state =
        state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]
    
    x0_ur_ui =
        [state_view[ur_ui_idx][ ur_IDX ]...;
         state_view[ur_ui_idx][ ui_IDX ]...]
        
    x0_vh_θh =
        [abs.(uh_state)...;
         angle.(uh_state)...]
    
    
    x0_vh_view =
        @view x0_vh_θh[vh_IDX]

    x0_θh_view  =
        @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view =
        @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0  =
        [ red_vh_θh_0_view; ]
     
    # ----------------------------------------------------
    
    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet =
        zeros(ComplexF64,
              length( uh_state ))
    
    Inet_view =
        view( Inet,
              1:length( Inet ) )

    Iinj =
        zeros(ComplexF64,
              length( uh_state ))

    Iinj_view =
        view(Iinj,
             1:length( Iinj ))
            
    # -------------------------------------------------      
    # Results
    # -------------------------------------------------    
    
    red_vh_θh = sol.u

    # -------------------------------------------------    

    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.(non_slack_gens_Idx_and_vh)
    
    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]
    
    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )
    
    # --------------------------------------------------- 

    uh = [uh_slack...;
          uh_non_slack...;
          uh_non_gens... ]

    #------------------------------------------

    Vm   = abs.(uh)
    Vθ   = angle.(uh)
    Vbus = uh
        
    # --------------------------------------------------- 

    I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
        uh,
        Ynet,
        nodes_node_idx_and_incident_edges_other_node_idx
    )


    if loc_load_exist == true

        Igen =
            I_sum_ynj_vj +
            get_gens_local_load_current(
                uh,
                P_Q_gens_loc_load_view,
                gens_with_loc_load_idx;
                loc_load_exist = loc_load_exist )
    else

        Igen = I_sum_ynj_vj
    end

    S_gens = x_from_xr_xi.( P_Q_gens_view )
    
    gens_S  =
        uh[gens_idx] .* conj.(Igen[gens_idx])

    P_Q_nodes_view[slack_bus_idx] .=
        [real(gens_S[slack_bus_idx]),
         imag(gens_S[slack_bus_idx])]

    for (idx, gen_S) in enumerate( gens_S)

        # P_Q_gens_view[idx][2] = imag(gen_S)
        
        P_Q_gens_view[idx] .= [ real( S_gens[idx]), imag(gen_S) ]
        
    end

    # ------------------------------------------------
    # update Inet_view  and Iinj_view
    # ------------------------------------------------

    Inet_view .= I_sum_ynj_vj

    # Iinj_view .=
    #     get_Iinj( uh,
    #               P_Q_non_gens_view,
    #               P_Q_gens_loc_load_view,
    #               Inet_view )

    Iinj_view .= get_Iinj( uh,
                           P_Q_non_gens_view,
                           P_Q_gens_loc_load_view,
                           ( gens_nodes_idx,
                             non_gens_nodes_idx,
                             gens_with_loc_load_idx ),
                           Inet_view;
                           loc_load_exist = loc_load_exist  )

    
    # ------------------------------------------------------ 
    # update nodes_pf_U_view
    # ------------------------------------------------------ 

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .=
            [ real( uk ), imag( uk ) ]

    end
    
    # ------------------------------------------------    
    # ------------------------------------------------    
    
    t_Sbus  = uh .* conj.( Inet_view )


    # Sbus_n  = uh .* conj.( Inet_inj )

    # GenSinj = Sbus_n + S_non_gens + S_gens_loc_load


    GenSinj =
            zeros(ComplexF64, length(t_Sbus))

    GenSinj .= t_Sbus 
    
    if loc_load_exist == true        

        GenSinj[non_gens_nodes_idx] .=
            t_Sbus[non_gens_nodes_idx] +
            x_from_xr_xi.( P_Q_non_gens_view )

        GenSinj[gens_with_loc_load_idx] .=
            t_Sbus[gens_with_loc_load_idx] +
            x_from_xr_xi.( P_Q_gens_loc_load_view )
        
    else
                
        GenSinj[non_gens_nodes_idx] .=
            t_Sbus[non_gens_nodes_idx] -
            x_from_xr_xi.( P_Q_non_gens_view )
        
    end
    
        
    Ifrom_Ito =
        get_Ifrom_Ito(
            nodes_pf_U_view,
            edges_Ybr_cal,
            edges_orientation)

    If = first.(Ifrom_Ito)
    
    It = last.(Ifrom_Ito)

    Ibranches = If + It


    # -------------------------------------------------
    
    branches_name =
        collect(keys( netd.edges ))
    
    nodes_name =
        collect(keys( netd.nodes ))    

    # -------------------------------------------------

    
    bus_dict_Iinj =
        OrderedDict(
            name => [ih_r, ih_i]
            for (name, ih_r, ih_i) in
                zip(nodes_name, real.( Iinj_view ),
                    imag.( Iinj_view )))

    bus_dict_Inet_inj =
        OrderedDict(
            name => [ih_r, ih_i]
            for (name, ih_r, ih_i) in
                zip(nodes_name,
                    real.(  Inet_view ),
                    imag.(  Inet_view )))     

    branch_dict_init =
        OrderedDict(
            name => (
                real(i_f),
                imag(i_f),
                real(i_t),
                imag(i_t) )
            for (name, i_f, i_t) in
                zip(branches_name, If, It))

    bus_dict_init =
        OrderedDict(
            name => (vh, θh, ph, qh, ih_r,
                     ih_i, pg, qg, ig_r, ig_i )
            for (name, vh, θh, ph, qh, ih_r,
                 ih_i, pg, qg, ig_r, ig_i ) in
                zip(
                    nodes_name,
                    abs.(uh),
                    angle.(uh),
                    real.( t_Sbus),
                    imag.( t_Sbus),
                    real.( Inet_view ),
                    imag.( Inet_view ),
                    real.( GenSinj),
                    imag.( GenSinj),
                    real.( Iinj_view ),
                    imag.( Iinj_view )) )

    #------------------------------------------


    named_tup_pf_result = (
        ; P_Q_gens_view,
        P_Q_non_gens_view,
        P_Q_gens_loc_load_view,
        Vm,
        Vθ,
        Vbus,
        Ibranches,
        Ybus ,
        t_Sbus,
        GenSinj,
        bus_dict_init,
        branch_dict_init,
        bus_dict_Iinj,
        bus_dict_Inet_inj,
        If,
        It ) 
    #----------------------------------------------------
    #----------------------------------------------------

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init =
        named_tup_pf_result.branch_dict_init

    pf_init_dict = Dict(
        "bus_dict_init" => bus_dict_init,
        "branch_dict_init" => branch_dict_init )
    
    #---------------------------------------------------
    #---------------------------------------------------

    """
    states init:

    init_operationpoint(
         netd, pf_init_dict)

    external_init_operationpoint(
         nd::NetworkData, bus_dict_init,
         branch_dict_init )

    industrial_model_init_operationpoint(
         netd,
         bus_dict_init
         ;pure = :pure,
         no_control_device = false )

    im_model_init_operationpoint(
         netd,
         bus_dict_init )
    
    """
    
    state .=
        external_init_operationpoint(
            netd,
            bus_dict_init,
            branch_dict_init )

    #---------------------------------------------------

    idq_θ_π_vhθh =
        [ get_pf_dynamic_idq_θ_π_vhθh(
            vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
          for ( vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in
              zip( abs.(uh[gens_nodes_idx]),
                   angle.(uh[gens_nodes_idx]),
                   δ_ω_ed_dash_eq_dash_view ,
                   ra_Xd_dash_Xq_dash_view ) ]

    # --------------------------------------------------- 

    # Igen = idq_θ_π_vhθh 
   
    #---------------------------------------------------

    nt_results = (
        ; Vm,
        Vθ,
        Vbus,
        idq_θ_π_vhθh,
        P_Q_gens_view,
        Inet_view,
        Iinj_view,
        state,
        named_tup_pf_result
         )
        
    return nt_results
    
    #---------------------------------------------------

    
end


# -------------------------------------------
# Get results
# -------------------------------------------

function get_pf_gens_uh_from_pf_sol(
    pf_sol,
    gens_uh_Q_from_red_sol_para)

    red_vh_θh = pf_sol.u

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) =
         gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx) =
         n2s_idxs

    
    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    pf_uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    pf_gens_uh = [ pf_uh[  idx  ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]
    
    return pf_gens_uh

end


function get_pf_gens_vh_θh_from_pf_sol(
    pf_sol,
    gens_uh_Q_from_red_sol_para)

    pf_gens_uh = get_pf_gens_uh_from_pf_sol(
        pf_sol,
        gens_uh_Q_from_red_sol_para)

    vec_pf_gens_vh_θh =  [
        [ abs( a_uh ), angle( a_uh ) ]
        for a_uh in
            pf_gens_uh ]
    
    return vec_pf_gens_vh_θh
    
end


function get_pf_gens_flat_vh_θh_from_pf_sol(
    pf_sol,
    gens_uh_Q_from_red_sol_para)

    pf_gens_uh = get_pf_gens_uh_from_pf_sol(
        pf_sol,
        gens_uh_Q_from_red_sol_para)

    flat_pf_gens_vh_θh = [
        abs.(pf_gens_uh);
        angle.(pf_gens_uh) ]
    
    return flat_pf_gens_vh_θh
end


function get_pf_gens_uh_from_pf_sol_u(
    pf_sol_u,
    gens_uh_Q_from_red_sol_para)

    red_vh_θh = pf_sol_u

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) =
         gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx) =
         n2s_idxs
    
    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    pf_uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    pf_gens_uh = [ pf_uh[  idx  ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]
    
    return pf_gens_uh
end



function get_pf_gens_vh_θh_from_pf_sol_u(
    pf_sol_u,
    gens_uh_Q_from_red_sol_para)

    pf_gens_uh =
        get_pf_gens_uh_from_pf_sol_u(
            pf_sol_u,
            gens_uh_Q_from_red_sol_para )

    pf_gens_vh_θh =  [
        [abs( a_uh ), angle( a_uh )]
          for a_uh in pf_gens_uh ]
    
    return pf_gens_vh_θh
end



function get_pf_gens_flat_vh_θh_from_pf_sol_u(
    pf_sol_u,
    gens_uh_Q_from_red_sol_para)

    pf_gens_uh =
        get_pf_gens_uh_from_pf_sol_u(
            pf_sol_u,
            gens_uh_Q_from_red_sol_para )

    pf_gens_flat_vh_θh =
        [abs.( pf_gens_uh );
         angle.( pf_gens_uh )]
    
    # flat_pf_gens_vh_θh_alt = [
    #     [ [abs( a_uh ), angle( a_uh )]
    #       for a_uh in pf_gens_uh ] ...; ]
    
    return pf_gens_flat_vh_θh
end




function get_red_vh_θh_idq_from_pf_sol(
    pf_sol,
    δ_ω_ed_dash_eq_dash_view,
    gens_uh_Q_from_red_sol_para )

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) = gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx) =
         n2s_idxs


    # ---------------------------------------------------
    
    red_vh_θh_idq = pf_sol.u

    red_vh_θh =
        red_vh_θh_idq[flat_red_vh_θh_0_Idx]
    
    idq_0 =
        red_vh_θh_idq[flat_idq_0_Idx]


    # ---------------------------------------------------
    
    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    pf_uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    pf_gens_uh = [ pf_uh[  idx  ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]

    # ---------------------------------------------------
    
    pf_gens_idq =
        [  get_pf_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( pf_gens_uh ),
                   angle.( pf_gens_uh ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    # ---------------------------------------------------
    
    flat_pf_gens_idq =
        [ first.(pf_gens_idq);
          second.(pf_gens_idq) ]

    flat_pf_gens_idq_alt = [pf_gens_idq...;]

    # ---------------------------------------------------
    
    red_vh_θh_idq =
        [ red_vh_θh;
          flat_pf_gens_idq ]

    # red_vh_θh_idq_alt =
    #     [ red_vh_θh;
    #       flat_pf_gens_idq_alt  ]

    return red_vh_θh_idq
    
end



function get_red_vh_θh_idq_from_pf_sol_u(
    pf_sol_u,
    δ_ω_ed_dash_eq_dash_view,
    gens_uh_Q_from_red_sol_para )

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) = gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx) =
         n2s_idxs

    
    # ---------------------------------------------------
    
    red_vh_θh_idq = pf_sol_u

    red_vh_θh =
        red_vh_θh_idq[flat_red_vh_θh_0_Idx]
    
    idq_0 =
        red_vh_θh_idq[flat_idq_0_Idx]


    # ---------------------------------------------------
    
    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    pf_uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    pf_gens_uh = [ pf_uh[  idx  ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]

    # ---------------------------------------------------
    
    pf_gens_idq =
        [  get_pf_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( pf_gens_uh ),
                   angle.( pf_gens_uh ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    # ---------------------------------------------------
    
    flat_pf_gens_idq = [ first.(pf_gens_idq); second.(pf_gens_idq) ]

    flat_pf_gens_idq_alt = [pf_gens_idq...;]

    # ---------------------------------------------------
    
    red_vh_θh_idq =
        [ red_vh_θh;
          flat_pf_gens_idq ]

    # red_vh_θh_idq_alt =
    #     [ red_vh_θh;
    #       flat_pf_gens_idq_alt  ]

    return red_vh_θh_idq
    
end



function get_pf_gens_uh_and_red_vh_θh_idq_from_pf_sol_u(
    pf_sol_u,
    δ_ω_ed_dash_eq_dash_view,
    gens_uh_Q_from_red_sol_para )

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) = gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         n2s_idxs

    # ---------------------------------------------------
    
    red_vh_θh_idq = pf_sol_u

    red_vh_θh =
        red_vh_θh_idq[flat_red_vh_θh_0_Idx]
    
    idq_0 =
        red_vh_θh_idq[flat_idq_0_Idx]


    # ---------------------------------------------------
    
    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    pf_uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    pf_gens_uh = [ pf_uh[  idx  ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]

    # ---------------------------------------------------
    
    pf_gens_idq =
        [  get_pf_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( pf_gens_uh ),
                   angle.( pf_gens_uh ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    # ---------------------------------------------------
    
    flat_pf_gens_idq =
        [ first.(pf_gens_idq);
          second.(pf_gens_idq) ]

    flat_pf_gens_idq_alt = [pf_gens_idq...;]

    # ---------------------------------------------------
    
    pf_red_vh_θh_idq =
        [ red_vh_θh;
          flat_pf_gens_idq ]

    return (; pf_red_vh_θh_idq, pf_gens_uh)
    
end


function get_pf_gens_flat_vh_θh_and_red_vh_θh_idq_from_pf_sol_u(
    pf_sol_u,
    δ_ω_ed_dash_eq_dash_view,
    gens_uh_Q_from_red_sol_para )

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) = gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         n2s_idxs

    # ---------------------------------------------------
    
    red_vh_θh_idq = pf_sol_u

    red_vh_θh =
        red_vh_θh_idq[flat_red_vh_θh_0_Idx]
    
    idq_0 =
        red_vh_θh_idq[flat_idq_0_Idx]


    # ---------------------------------------------------
    
    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    pf_uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    pf_gens_uh = [ pf_uh[  idx  ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]

    
    # pf_gens_vh_θh = [ [abs(a_uh), angle(a_uh)]
    #                for a_uh  in
    #                    pf_gens_uh ]

    # ---------------------------------------------------
    
    pf_gens_idq =
        [  get_pf_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( pf_gens_uh ),
                   angle.( pf_gens_uh ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    # ---------------------------------------------------
    
    flat_pf_gens_idq = [ first.(pf_gens_idq);
                         second.(pf_gens_idq) ]

    # flat_pf_gens_idq_alt = [pf_gens_idq...;]

    pf_gens_flat_vh_θh =
        [abs.( pf_gens_uh ) ;
         angle.( pf_gens_uh )]
        

    # ---------------------------------------------------
    
    pf_red_vh_θh_idq =
        [ red_vh_θh;
          flat_pf_gens_idq ]

    return (; pf_red_vh_θh_idq, pf_gens_flat_vh_θh  )
    
end




function get_uh_vh_θh_idq_from_pf_sol(
    pf_sol,
    δ_ω_ed_dash_eq_dash_view,
    gens_uh_Q_from_red_sol_para )

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) = gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         n2s_idxs

    # ---------------------------------------------------
    
    red_vh_θh_idq = pf_sol.u

    red_vh_θh =
        red_vh_θh_idq[flat_red_vh_θh_0_Idx]
    
    idq_0 =
        red_vh_θh_idq[flat_idq_0_Idx]


    # ---------------------------------------------------
    
    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    pf_uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    pf_gens_uh = [ pf_uh[  idx  ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]
    
    # ---------------------------------------------------
    
    pf_gens_idq =
        [  get_pf_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( pf_gens_uh ),
                   angle.( pf_gens_uh ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    return (; pf_gens_uh, pf_gens_idq)
    
end



function get_nodes_uh_from_pf_sol_u(
    pf_sol_u,    
    gens_uh_Q_from_red_sol_para )

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) = gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         n2s_idxs

    # ---------------------------------------------------
    
    red_vh_θh_idq = pf_sol_u

    red_vh_θh =
        red_vh_θh_idq[flat_red_vh_θh_0_Idx]

    # ---------------------------------------------------
    
    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    pf_uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    return pf_uh
    
end



function get_nodes_ur_ui_from_pf_sol_u(
    pf_sol_u,
    gens_uh_Q_from_red_sol_para )

    pf_uh = get_nodes_uh_from_pf_sol_u(
        pf_sol_u,    
        gens_uh_Q_from_red_sol_para )

    nodes_ur_ui =
        [[real(a_uh), imag(a_uh)]
         for a_uh in pf_uh ]
    
    return nodes_ur_ui
    
end



function get_nodes_vh_θh_from_pf_sol_u(
    pf_sol_u,
    gens_uh_Q_from_red_sol_para )

    pf_uh = get_nodes_uh_from_pf_sol_u(
        pf_sol_u,    
        gens_uh_Q_from_red_sol_para )
    
    nodes_vh_θh =
        [ [ abs(a_uh), angle(a_uh) ]
         for a_uh in pf_uh ]
    
    return nodes_vh_θh
    
end


function get_results_from_vh_θh_pf_sol(
    pf_sol,
    δ_ω_ed_dash_eq_dash,
    gens_uh_Q_from_red_sol_para,
    full_idxs_and_full_vh_θh )
    
    #-------------------------------

    return get_results_from_vh_θh_pf_sol_u(
        pf_sol.u, δ_ω_ed_dash_eq_dash,
        gens_uh_Q_from_red_sol_para,
        full_idxs_and_full_vh_θh )

end
 

function get_results_from_vh_θh_pf_sol_u(
    pf_sol_u,
    δ_ω_ed_dash_eq_dash,
    gens_uh_Q_from_red_sol_para,
    full_idxs_and_full_vh_θh )
    
    #-------------------------------

    full_vh_θh = pf_sol_u
    
    #-------------------------------

    ra_Xd_dash_Xq_dash =
        gens_uh_Q_from_red_sol_para.ra_Xd_dash_Xq_dash

    #-------------------------------

    full_gens_vh_Idxs =
        full_idxs_and_full_vh_θh.full_gens_vh_Idxs

    full_gens_θh_Idxs =
        full_idxs_and_full_vh_θh.full_gens_θh_Idxs

    full_non_gens_nodes_vh_Idxs =
        full_idxs_and_full_vh_θh.full_non_gens_nodes_vh_Idxs

    full_non_gens_nodes_θh_Idxs =
        full_idxs_and_full_vh_θh.full_non_gens_nodes_θh_Idxs
    
    #-------------------------------
    
    full_gens_vh =
        full_vh_θh[ full_gens_vh_Idxs ]

    full_gens_θh =
        full_vh_θh[ full_gens_θh_Idxs ]

    non_gens_vh = full_non_gens_nodes_vh =
        full_vh_θh[ full_non_gens_nodes_vh_Idxs ]

    non_gens_θh  = full_non_gens_nodes_θh =
        full_vh_θh[ full_non_gens_nodes_θh_Idxs ]

    vh = [full_gens_vh; full_non_gens_nodes_vh]

    θh = [full_gens_θh; full_non_gens_nodes_θh ]


    #--------------------------------------

    gens_δ = first.( δ_ω_ed_dash_eq_dash )

    gens_ed_dash = third.( δ_ω_ed_dash_eq_dash )

    gens_eq_dash = fourth.( δ_ω_ed_dash_eq_dash )

    # -------------------------------------

    gens_ra = first.(ra_Xd_dash_Xq_dash)

    gens_Xd_dash = second.(ra_Xd_dash_Xq_dash)

    gens_Xq_dash = third.(ra_Xd_dash_Xq_dash)        

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
        for (a_vh, a_θh,
             a_δ, a_ed_dash, a_eq_dash,
             a_ra, a_X_d_dash, a_X_q_dash ) in
            zip( full_gens_vh, full_gens_θh,
                 gens_δ, gens_ed_dash, gens_eq_dash,
                 gens_ra, gens_Xd_dash, gens_Xq_dash )  ]

    gens_i_d = first.(gens_id_iq )

    gens_i_q = last.( gens_id_iq  )


    
    P_gens_cal = [
        get_a_gen_active_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh,
            a_δ,
            a_id, a_iq )
        for ( a_vh, a_θh, a_δ,  a_id, a_iq ) in
            zip( full_gens_vh, full_gens_θh,
                 gens_δ,
                 gens_i_d, gens_i_q ) ]

    Q_gens_cal = [
        get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh,
            a_δ,  a_id, a_iq )
        for ( a_vh, a_θh, a_δ,  a_id, a_iq ) in
            zip( full_gens_vh, full_gens_θh,
                 gens_δ,
                 gens_i_d, gens_i_q ) ]

    # -------------------------------------

    return (;vh,
            θh,
            
            gens_i_d,
            gens_i_q,
            
            P_gens_cal,
            Q_gens_cal )
end
 

# -------------------------------------


function get_results_from_vh_θh_id_iq_pf_sol(
    pf_sol,
    δ_ω_ed_dash_eq_dash,
    gens_uh_Q_from_red_sol_para,
    intg_idxs_and_full_vh_θh )
    
    #-------------------------------

    return get_results_from_vh_θh_id_iq_pf_sol_u(
        pf_sol.u, δ_ω_ed_dash_eq_dash,
        gens_uh_Q_from_red_sol_para,
        intg_idxs_and_full_vh_θh )
end


function get_results_from_vh_θh_id_iq_pf_sol_u(
    pf_sol_u,
    δ_ω_ed_dash_eq_dash,
    gens_uh_Q_from_red_sol_para,
    intg_idxs_and_full_vh_θh )
    
    #-------------------------------

    intg_vh_θh_id_iq = pf_sol_u
    
    #-------------------------------

    ra_Xd_dash_Xq_dash =
        gens_uh_Q_from_red_sol_para.ra_Xd_dash_Xq_dash

    #-------------------------------

    intg_gens_vh_Idxs =
        intg_idxs_and_full_vh_θh.intg_gens_vh_Idxs

    intg_gens_θh_Idxs =
        intg_idxs_and_full_vh_θh.intg_gens_θh_Idxs

    intg_non_gens_nodes_vh_Idxs =
        intg_idxs_and_full_vh_θh.intg_non_gens_nodes_vh_Idxs

    intg_non_gens_nodes_θh_Idxs =
        intg_idxs_and_full_vh_θh.intg_non_gens_nodes_θh_Idxs

    intg_gen_id_Idxs =
        intg_idxs_and_full_vh_θh.intg_gen_id_Idxs

    intg_gen_iq_Idxs =
        intg_idxs_and_full_vh_θh.intg_gen_iq_Idxs
    
    #-------------------------------
    
    intg_gens_vh =
        intg_vh_θh_id_iq[ intg_gens_vh_Idxs ]

    intg_gens_θh =
        intg_vh_θh_id_iq[ intg_gens_θh_Idxs ]

    non_gens_vh = intg_non_gens_nodes_vh =
        intg_vh_θh_id_iq[ intg_non_gens_nodes_vh_Idxs ]

    non_gens_θh  = intg_non_gens_nodes_θh =
        intg_vh_θh_id_iq[ intg_non_gens_nodes_θh_Idxs ]

    vh = [intg_gens_vh; intg_non_gens_nodes_vh]

    θh = [intg_gens_θh; intg_non_gens_nodes_θh ]

    intg_gens_id =
        intg_vh_θh_id_iq[intg_gen_id_Idxs]

    intg_gens_iq =
        intg_vh_θh_id_iq[ intg_gen_iq_Idxs ]
    
    #--------------------------------------

    gens_δ = first.( δ_ω_ed_dash_eq_dash )

    gens_ed_dash = third.( δ_ω_ed_dash_eq_dash )

    gens_eq_dash = fourth.( δ_ω_ed_dash_eq_dash )

    # -------------------------------------

    gens_ra = first.(ra_Xd_dash_Xq_dash)

    gens_Xd_dash = second.(ra_Xd_dash_Xq_dash)

    gens_Xq_dash = third.(ra_Xd_dash_Xq_dash)        

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
        for (a_vh, a_θh,
             a_δ, a_ed_dash, a_eq_dash,
             a_ra, a_X_d_dash, a_X_q_dash ) in
            zip( intg_gens_vh, intg_gens_θh,
                 gens_δ, gens_ed_dash, gens_eq_dash,
                 gens_ra, gens_Xd_dash, gens_Xq_dash )  ]

    gens_i_d = first.(gens_id_iq )

    gens_i_q = last.( gens_id_iq  )

    
    P_gens_cal = [
        get_a_gen_active_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh,
            a_δ,
            a_id, a_iq )
        for ( a_vh, a_θh, a_δ,  a_id, a_iq ) in
            zip( intg_gens_vh, intg_gens_θh,
                 gens_δ,
                 gens_i_d, gens_i_q ) ]

    Q_gens_cal = [
        get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh,
            a_δ,  a_id, a_iq )
        for ( a_vh, a_θh, a_δ,  a_id, a_iq ) in
            zip( intg_gens_vh, intg_gens_θh,
                 gens_δ,
                 gens_i_d, gens_i_q ) ]

    # -------------------------------------
    
    P_gens_intg = [
        get_a_gen_active_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh,
            a_δ,
            a_id, a_iq )
        for ( a_vh, a_θh, a_δ,  a_id, a_iq ) in
            zip( intg_gens_vh, intg_gens_θh,
                 gens_δ,
                 intg_gens_id, intg_gens_iq ) ]

    Q_gens_intg = [
        get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh,
            a_δ,  a_id, a_iq )
        for ( a_vh, a_θh, a_δ,  a_id, a_iq ) in
            zip( intg_gens_vh, intg_gens_θh,
                 gens_δ,
                 intg_gens_id, intg_gens_iq ) ]
    
    # -------------------------------------
        
    return (;vh,
            θh,
            
            gens_i_d,
            gens_i_q,
            
            intg_gens_id,
            intg_gens_iq,
            
            P_gens_cal,
            Q_gens_cal,
            
            P_gens_intg,
            Q_gens_intg
            )
end
 


function get_results_from_red_vh_θh_iq_iq_pf_sol(
    pf_sol,
    δ_ω_ed_dash_eq_dash,
    gens_uh_Q_from_red_sol_para )


    return get_results_from_red_vh_θh_iq_iq_pf_sol_u(
        pf_sol.u, δ_ω_ed_dash_eq_dash,
        gens_uh_Q_from_red_sol_para )
    
end


function get_results_from_red_vh_θh_iq_iq_pf_sol_u(
    pf_sol,
    δ_ω_ed_dash_eq_dash,
    gens_uh_Q_from_red_sol_para )

    red_vh_θh_iq_iq = pf_sol.u

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) =
         gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         n2s_idxs

    # ---------------------------------------------------

    non_slack_gens_vh = last.(
            non_slack_gens_Idx_and_vh)

    # ---------------------------------------------------
    
    idq_flat =
            red_vh_θh_iq_iq[
                flat_idq_0_Idx ]
    
    gens_id = idq_flat[ gens_id_Idx ]

    gens_iq = idq_flat[ gens_iq_Idx ]

    # ---------------------------------------------------
    
    red_vh_θh =
            red_vh_θh_iq_iq[
                flat_red_vh_θh_0_Idx ]
    
    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]
    
    non_slack_gens_θh =
        red_vh_θh[
            non_slack_gens_θh_idx2Idx ]

    non_gens_θh =
        red_vh_θh[
            non_gens_θh_idx2Idx ]        


    vh = [
        idx ∈ slack_gens_nodes_idx ?
            slack_vh :
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_vh[
                n2s_non_slack_gens_idx[ idx]] : 
            non_gens_vh[ n2s_non_gens_idx[ idx ] ]
               for idx in 1:nodes_size ]

    θh = [
        idx ∈ slack_gens_nodes_idx ?
            0.0 :
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx]] : 
            non_gens_θh[ n2s_non_gens_idx[ idx ] ]
               for idx in 1:nodes_size ]

    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]

    
    uh = vh .* exp.(im * θh )

    gens_uh = uh[gens_nodes_idx]

        
    # ---------------------------------------------------

    gens_δ = first.( δ_ω_ed_dash_eq_dash )

    # -------------------------------------
    
    P_gens_red = [
        get_a_gen_active_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh,
            a_δ,
            a_id, a_iq )
        for ( a_vh, a_θh, a_δ,  a_id, a_iq ) in
            zip( gens_vh, gens_θh,
                 gens_δ,
                 gens_id, gens_iq ) ]

    Q_gens_red = [
        get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh,
            a_δ,  a_id, a_iq )
        for ( a_vh, a_θh, a_δ,  a_id, a_iq ) in
            zip( gens_vh, gens_θh,
                 gens_δ,
                 gens_id, gens_iq ) ]

    # ---------------------------------------------------
    # ---------------------------------------------------
    
    return (; vh,
            θh,
            gens_id,
            gens_iq,
            P_gens_red,
            Q_gens_red )
end
 


# ---------------------------------------------------


function get_results_S_gens_from_pf_sta_red_sol_u(
    pf_sol,
    red_sol_para )

    #-------------------------------

    red_vh_θh = pf_sol.u

    #-------------------------------
    
    (; edges_Ybr_cal,
     edges_orientation,
     nodes_name,
     branches_name,
     all_nodes_idx,
     pf_kw_para) =
          red_sol_para

    #-------------------------------

    (;loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs) =
         pf_kw_para
        
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh) =
         pf_kw_gens_vh_slack_θh_para

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          pf_kw_net_para

     (; red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          pf_kw_var_idxs

     (; P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs ) =
          pf_kw_PQ_para_idxs

     (;slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx #
      ) =
          pf_kw_nodes_types_idxs    

     (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx ) =
         pf_kw_n2s_idxs
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )

    # all_nodes_idx =
    #     sort([gens_nodes_idx; non_gens_nodes_idx])
    
    #-------------------------------

    non_slack_gens_θh =
        red_vh_θh[
            red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[
            red_non_gens_θh_idx2Idx ]

    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ] ]
        for idx in all_nodes_idx ]

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]
    
    #-------------------------------

    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,
            
            Ynet,            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx )

    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx )
        

    #-------------------------------

    gens_nodes_network_current =
        [ get_a_node_∑_ynj_x_vj(
            vh, θh, Ynet[ gen_idx ],
            nodes_idx_with_adjacent_nodes_idx[ gen_idx ],
            n2s_all_nodes_idx )
          for gen_idx in
              gens_nodes_idx ]

    non_gens_nodes_network_current =
        [ get_a_node_∑_ynj_x_vj(
            vh, θh, Ynet[ non_gen_idx ],
            nodes_idx_with_adjacent_nodes_idx[non_gen_idx],
            n2s_all_nodes_idx )
          for non_gen_idx in
              non_gens_nodes_idx ]
    
    nodes_network_current =
        [ get_a_node_∑_ynj_x_vj(
            vh, θh, Ynet[ idx ],
            nodes_idx_with_adjacent_nodes_idx[ idx ],
            n2s_all_nodes_idx )
          for idx in
              all_nodes_idx ]
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    
    
    #-------------------------------

    return (; vh, θh, gens_uh,
            gens_current_injection,
            gens_loc_load_current,
            gens_nodes_network_current,
            non_gens_nodes_network_current,
            nodes_network_current,
            pf_P_gens,
            pf_Q_gens )
end



function get_results_from_pf_sta_red_sol_u(
    pf_sol,
    red_sol_para )

    #-------------------------------    
    
    (; edges_Ybr_cal,
     edges_orientation,
     nodes_name,
     branches_name,
     pf_kw_para,
     P_g_loc_load,
     Q_g_loc_load,
     P_non_gens,
     Q_non_gens) =
          red_sol_para

    red_vh_θh = pf_sol.u

    #-------------------------------

    (;loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs) =
         pf_kw_para
        
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh) =
         pf_kw_gens_vh_slack_θh_para

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          pf_kw_net_para

     (; red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          pf_kw_var_idxs

     (; P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs) =
          pf_kw_PQ_para_idxs

     (;slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx) =
          pf_kw_nodes_types_idxs    

     (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx ) =
         pf_kw_n2s_idxs
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------
    
    # non_slack_gens_vh = last.(
    #     non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[
            red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[
            red_non_gens_θh_idx2Idx ]

    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ] ]
        for idx in all_nodes_idx ]
    

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]
    
    #-------------------------------

    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx )

    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx )

    #-------------------------------

    gens_nodes_network_current =
        [ get_a_node_∑_ynj_x_vj(
            vh, θh, Ynet[ gen_idx ],
            nodes_idx_with_adjacent_nodes_idx[ gen_idx ],
            n2s_all_nodes_idx )
          for gen_idx in
              gens_nodes_idx ]

    non_gens_nodes_network_current =
        [ get_a_node_∑_ynj_x_vj(
            vh, θh, Ynet[ non_gen_idx ],
            nodes_idx_with_adjacent_nodes_idx[non_gen_idx],
            n2s_all_nodes_idx )
          for non_gen_idx in
              non_gens_nodes_idx ]
    
    nodes_network_current =
        [ get_a_node_∑_ynj_x_vj(
            vh, θh, Ynet[ idx ],
            nodes_idx_with_adjacent_nodes_idx[idx],
            n2s_all_nodes_idx )
          for idx in
              all_nodes_idx ]

    #-------------------------------

    S_gens =
        gens_uh .* conj.(gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    
    
    #-------------------------------

    Inet = get_Inet_inj(
        uh, Ynet,
        nodes_idx_with_adjacent_nodes_idx )

    Iinj = get_Iinj(
        uh, P_non_gens,
        Q_non_gens,
        P_g_loc_load, Q_g_loc_load,
        pf_kw_nodes_types_idxs, pf_kw_n2s_idxs,    
        loc_load_exist, Inet)
    
    Sbus_n  = uh .* conj.( nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n, P_non_gens, Q_non_gens,
        P_g_loc_load, Q_g_loc_load,
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs,    
        loc_load_exist )

    Ifrom_Ito = get_I_from_I_to(
        uh, edges_Ybr_cal, edges_orientation)

    If = first.(Ifrom_Ito)
    It = last.(Ifrom_Ito)

    Ibranches = If + It
    
    #-------------------------------

    bus_dict_Iinj = OrderedDict( name => [ih_r, ih_i] for (name, ih_r, ih_i) in zip(nodes_name, real.( Iinj ), imag.( Iinj )))

    bus_dict_Inet_inj = OrderedDict( name => [ih_r, ih_i] for (name, ih_r, ih_i) in zip(nodes_name, real.(  Inet ), imag.(  Inet )))     

    branch_dict_init = OrderedDict(name => ( real(i_f), imag(i_f), real(i_t), imag(i_t) ) for (name, i_f, i_t) in zip(branches_name, If, It))

    bus_dict_init = OrderedDict( name => (vh, θh, ph, qh, ih_r, ih_i, pg, qg, ig_r, ig_i ) for (name, vh, θh, ph, qh, ih_r, ih_i, pg, qg, ig_r, ig_i ) in zip(nodes_name,  abs.(uh),  angle.(uh), real.(Sbus_n), imag.(Sbus_n), real.( Inet ), imag.( Inet ), real.(GenSinj), imag.(GenSinj), real.( Iinj), imag.( Iinj )) )


    #-------------------------------

    pf_init_dict = Dict(
        "bus_dict_init" => bus_dict_init,
        "branch_dict_init" => branch_dict_init )

    #-------------------------------
    
    Vm   = vh
    Vθ   = θh
    Vbus = uh

    #-------------------------------
    #-------------------------------
    
    # return (; vh,
    #         θh,
    #         P_gens
    #         Q_gens )
    
    return  (; gens_current_injection,
             gens_loc_load_current,
             gens_nodes_network_current,
             non_gens_nodes_network_current,
             nodes_network_current,
             
             S_gens,
             pf_P_gens,
             pf_Q_gens,
             
             Inet,
             Iinj,
             
             Sbus_n,
             GenSinj,
             
             If,
             It,
             Ibranches,
             
             bus_dict_Iinj,
             bus_dict_Inet_inj,
             branch_dict_init,
             bus_dict_init,
             
             pf_init_dict,
             
             Vm, Vθ, Vbus )
    
end



function get_generic_results_pf_spcm_red_sol_u(
    pf_sol;
    generic_red_sol_kwd_para =
        generic_red_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0 )

    #-------------------------------
    # pu

    Ibase =  baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    red_vh_θh = pf_sol.u

    #-------------------------------    
    
    (; Ybr_cal_and_edges_orientation,
     sta_pf_PQ_para,
     ode_gens_generic_para,
     pf_kw_para
      ) =
          generic_red_sol_kwd_para

    #-------------------------------

    (; edges_Ybr_cal,
     edges_orientation) =
         Ybr_cal_and_edges_orientation    

    #-------------------------------

    # (;P_gens,
    #  Q_gens,
    #  P_non_gens,
    #  Q_non_gens,
    #  P_g_loc_load,
    #  Q_g_loc_load,
    #  loc_load_exist) =
    #      sta_pf_PQ_para

    P_non_gens = sta_pf_PQ_para.P_non_gens
    Q_non_gens = sta_pf_PQ_para.Q_non_gens
    P_g_loc_load = sta_pf_PQ_para.P_g_loc_load
    Q_g_loc_load = sta_pf_PQ_para.Q_g_loc_load
    
    #-------------------------------


    X_d_dash = ode_gens_generic_para.X_d_dash

    X_q_dash = ode_gens_generic_para.X_q_dash

    X_d = ode_gens_generic_para.X_d

    X_q = ode_gens_generic_para.X_q
    
    #-------------------------------

    (;loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs) =
         pf_kw_para
        
    #-------------------------------

    (;slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh) =
         pf_kw_gens_vh_slack_θh_para

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          pf_kw_net_para

     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          pf_kw_var_idxs

     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs) =
          pf_kw_PQ_para_idxs

     (;slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx) =
          pf_kw_nodes_types_idxs    

     (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx ) =
         pf_kw_n2s_idxs

    
    #-------------------------------

    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]


    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------
    
    non_slack_gens_θh =
        red_vh_θh[
            red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[
            red_non_gens_θh_idx2Idx ]

    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ] ]
        for idx in all_nodes_idx ]

    θh_deg = (180/pi) * θh  

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]
    
    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_Bnj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx )

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx )
    
    gens_current_injection =
        get_gens_current_injection_by_Bnj(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx )

    Igen = gens_nodes_network_current +
        gens_loc_load_current
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    

    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens = real.( gens_S )

    pf_Q_g_gens = imag.( gens_S )

    #-------------------------------

    Inet = get_Inet_inj_by_Bnj(
        uh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    Iinj = get_Iinj(
        uh,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs,    
        loc_load_exist,
        
        Inet)
    
    Sbus_n  = uh .* conj.( nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs,    
        loc_load_exist )

    Ifrom_Ito = get_I_from_I_to_by_Bnj(
        uh, edges_Ybr_cal,
        edges_orientation,
        n2s_all_nodes_idx)

    If = first.(Ifrom_Ito)
    It = last.(Ifrom_Ito)

    Ibranches = If + It
    
    #-------------------------------
    
    # Vm   = vh
    # Vθ   = θh
    
    Vbus = uh

    #-------------------------------

    # pu_X_d_dash = X_d_dash ./ Zbase
    
    pu_Igen = Igen ./ Ibase
    
    E =
        Vbus[transformed_gens_nodes_idx] .+
        im * X_d_dash 
    
    #-------------------------------

    gens_vh = vh[gens_nodes_idx]
    
    gens_θh = θh[gens_nodes_idx]
    
    #-------------------------------
    
    return  (; gens_current_injection,
             gens_loc_load_current,
             gens_nodes_network_current,
             non_gens_nodes_network_current,
             nodes_network_current,
             
             S_gens,
             pf_P_gens, pf_Q_gens,
             gens_S,
             
             pf_P_g_gens,
             pf_Q_g_gens,

             Igen, Inet, Iinj, pu_Igen,
             
             Sbus_n,
             GenSinj,
             
             If, It, Ibranches,
             
             vh, θh, θh_deg, Vbus, E,

             gens_vh, gens_θh,
             
             gens_nodes_idx )
    
end



function get_results_static_pf_red_sol_u(
    pf_sol;
    generic_red_sol_kwd_para =
        generic_red_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0)

    #-------------------------------
    # pu

    Ibase =  baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    red_vh_θh = pf_sol.u

    #-------------------------------    
    
    (;Ybr_cal_and_edges_orientation,
     sta_pf_PQ_para,
     pf_kw_para,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
          NamedTupleTools.select(
              generic_red_sol_kwd_para,
              (:Ybr_cal_and_edges_orientation,
               :sta_pf_PQ_para,
               :pf_kw_para,
               :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #-------------------------------

    (;edges_Ybr_cal,
     edges_orientation) =
         NamedTupleTools.select(
             Ybr_cal_and_edges_orientation,
             (:edges_Ybr_cal,
              :edges_orientation))
                                

    #-------------------------------

    # (;P_gens,
    #  Q_gens,
    #  P_non_gens,
    #  Q_non_gens,
    #  P_g_loc_load,
    #  Q_g_loc_load,
    #  loc_load_exist) =
    #      sta_pf_PQ_para

    # P_non_gens = sta_pf_PQ_para.P_non_gens
    # Q_non_gens = sta_pf_PQ_para.Q_non_gens
    # P_g_loc_load = sta_pf_PQ_para.P_g_loc_load
    # Q_g_loc_load = sta_pf_PQ_para.Q_g_loc_load



    (;P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             (:P_non_gens,
              :Q_non_gens,
              :P_g_loc_load,
              :Q_g_loc_load))
    
    #-------------------------------

    (;loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     # pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs) =
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :pf_kw_gens_vh_slack_θh_para,
              # :pf_kw_net_para,
              :pf_kw_var_idxs,
              :pf_kw_PQ_para_idxs,
              :pf_kw_nodes_types_idxs,
              :pf_kw_n2s_idxs))
        
    #-------------------------------

    (;slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh) =
         NamedTupleTools.select(
             pf_kw_gens_vh_slack_θh_para,
             (:slack_gens_vh,
              :slack_gens_θh,

              :gens_vh,
              :non_slack_gens_vh))

     # (;Ynet,
     #  nodes_idx_with_adjacent_nodes_idx) =
     #      NamedTupleTools.select(
     #          Pf_kw_net_para,
     #          (:Ynet,
     #           :nodes_idx_with_adjacent_nodes_idx))

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              Ynet_wt_nodes_idx_wt_adjacent_nodes,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))

    
     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          NamedTupleTools.select(
              pf_kw_var_idxs,
              (:red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx ) )

     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs) =
          NamedTupleTools.select(
              pf_kw_PQ_para_idxs,
              (:P_gens_sta_para_Idxs,
               :Q_gens_sta_para_Idxs,
               :P_non_gens_sta_para_Idxs,
               :Q_non_gens_sta_para_Idxs,
               :P_g_loc_load_sta_para_Idxs,
               :Q_g_loc_load_sta_para_Idxs))

     (;slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx) =
          NamedTupleTools.select(
              pf_kw_nodes_types_idxs,
              (:slack_gens_nodes_idx,
               :non_slack_gens_nodes_idx,
               :gens_nodes_idx,
               :non_gens_nodes_idx,
               :gens_with_loc_load_idx,
               :all_nodes_idx))    

     (;n2s_slack_gens_idx,
      n2s_non_slack_gens_idx,
      n2s_gens_idx,
      n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx ) =
          NamedTupleTools.select(
              pf_kw_n2s_idxs,
              (:n2s_slack_gens_idx,
               :n2s_non_slack_gens_idx,
               :n2s_gens_idx,
               :n2s_non_gens_idx,
               :n2s_gens_with_loc_load_idxs,
               :n2s_all_nodes_idx ))

    
    #-------------------------------

    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]


    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]

    #-------------------------------
    
    non_slack_gens_θh =
        red_vh_θh[
            red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[
            red_non_gens_θh_idx2Idx ]

    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ] ]
        for idx in all_nodes_idx ]

    θh_deg = (180/pi) * θh  

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]
    
    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_ynj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx )

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx )
    
    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx )

    Igen = gens_nodes_network_current +
        gens_loc_load_current
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(
            gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    

    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens = real.( gens_S )

    pf_Q_g_gens = imag.( gens_S )

    #-------------------------------

    Inet = get_Inet_inj(
        uh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    Iinj = get_Iinj(
        uh,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs,    
        loc_load_exist,
        
        Inet)
    
    Sbus_n  = uh .* conj.(
        nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs,    
        loc_load_exist )

    Ifrom_Ito = get_I_from_I_to(
        uh, edges_Ybr_cal,
        edges_orientation,
        n2s_all_nodes_idx)

    If = first.(Ifrom_Ito)
    
    It = last.(Ifrom_Ito)

    Ibranches = If + It
    
    #-------------------------------
    
    # Vm   = vh
    # Vθ   = θh
    
    Vbus = uh

    #-------------------------------
    #-------------------------------

    gens_vh = vh[gens_nodes_idx]
    
    gens_θh = θh[gens_nodes_idx]
    
    #-------------------------------
    
    return  (;gens_current_injection,
             gens_loc_load_current,
             gens_nodes_network_current,
             non_gens_nodes_network_current,
             nodes_network_current,
             
             S_gens,
             
             pf_P_gens, pf_Q_gens,
             gens_S,
             
             pf_P_g_gens,
             pf_Q_g_gens,

             Igen, Inet, Iinj,
             
             Sbus_n,
             GenSinj,
             
             If, It, Ibranches,
             
             vh, θh, θh_deg,
             
             Vbus,
             
             gens_vh, gens_θh,

             gens_nodes_idx,

             transformed_slack_gens_nodes_idx,
             transformed_gens_nodes_idx,
             transformed_non_gens_nodes_idx,
             transformed_all_nodes_idx)
    
end


function get_generic_results_pf_sta_red_sol_u(
    pf_sol;
    generic_red_sol_kwd_para =
        generic_red_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0)

    #-------------------------------
    # pu

    Ibase =  baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    red_vh_θh = pf_sol.u

    #-------------------------------    
    
    (;Ybr_cal_and_edges_orientation,
     sta_pf_PQ_para,
     ode_gens_generic_para,
     pf_kw_para ) =
          NamedTupleTools.select(
              generic_red_sol_kwd_para,
              (:Ybr_cal_and_edges_orientation,
               :sta_pf_PQ_para,
               :ode_gens_generic_para,
               :pf_kw_para))

    #-------------------------------

    (; edges_Ybr_cal,
     edges_orientation) =
         Ybr_cal_and_edges_orientation    

    #-------------------------------

    # (;P_gens,
    #  Q_gens,
    #  P_non_gens,
    #  Q_non_gens,
    #  P_g_loc_load,
    #  Q_g_loc_load,
    #  loc_load_exist) =
    #      sta_pf_PQ_para

    P_non_gens = sta_pf_PQ_para.P_non_gens
    Q_non_gens = sta_pf_PQ_para.Q_non_gens
    P_g_loc_load = sta_pf_PQ_para.P_g_loc_load
    Q_g_loc_load = sta_pf_PQ_para.Q_g_loc_load
    
    #-------------------------------


    X_d_dash = ode_gens_generic_para.X_d_dash

    X_q_dash = ode_gens_generic_para.X_q_dash

    X_d = ode_gens_generic_para.X_d

    X_q = ode_gens_generic_para.X_q
    
    #-------------------------------

    (;loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs) =
         pf_kw_para
        
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh) =
         pf_kw_gens_vh_slack_θh_para

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          pf_kw_net_para

     (; red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          pf_kw_var_idxs

     (; P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs) =
          pf_kw_PQ_para_idxs

     (;slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx) =
          pf_kw_nodes_types_idxs    

     (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx ) =
         pf_kw_n2s_idxs

    
    #-------------------------------

    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]


    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------
    
    non_slack_gens_θh =
        red_vh_θh[
            red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[
            red_non_gens_θh_idx2Idx ]

    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ] ]
        for idx in all_nodes_idx ]

    θh_deg = (180/pi) * θh  

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]
    
    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_ynj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx )

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx )
    
    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx )

    Igen = gens_nodes_network_current +
        gens_loc_load_current
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(
            gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    

    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens = real.( gens_S )

    pf_Q_g_gens = imag.( gens_S )

    #-------------------------------

    Inet = get_Inet_inj(
        uh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    Iinj = get_Iinj(
        uh,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs,    
        loc_load_exist,
        
        Inet)
    
    Sbus_n  = uh .* conj.(
        nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs,    
        loc_load_exist )

    Ifrom_Ito = get_I_from_I_to(
        uh, edges_Ybr_cal,
        edges_orientation,
        n2s_all_nodes_idx)

    If = first.(Ifrom_Ito)
    
    It = last.(Ifrom_Ito)

    Ibranches = If + It
    
    #-------------------------------
    
    # Vm   = vh
    # Vθ   = θh
    
    Vbus = uh

    #-------------------------------

    # pu_X_d_dash = X_d_dash ./ Zbase
    
    pu_Igen = Igen ./ Ibase
    
    gens_E =
        Vbus[transformed_gens_nodes_idx] .+
        im * X_d_dash .* pu_Igen

    gens_mag_E = abs.(gens_E )
    gens_ang_E = angle.(gens_E)

    #-------------------------------

    gens_vh = vh[gens_nodes_idx]
    
    gens_θh = θh[gens_nodes_idx]
    
    #-------------------------------
    
    return  (;gens_current_injection,
             gens_loc_load_current,
             gens_nodes_network_current,
             non_gens_nodes_network_current,
             nodes_network_current,
             
             S_gens,
             pf_P_gens, pf_Q_gens,
             gens_S,
             
             pf_P_g_gens,
             pf_Q_g_gens,

             Igen, Inet, Iinj, pu_Igen,
             
             Sbus_n,
             GenSinj,
             
             If, It, Ibranches,
             
             vh, θh, θh_deg, Vbus, gens_E,
             gens_mag_E, gens_ang_E,
             

             gens_vh, gens_θh,
             
             gens_nodes_idx )
    
end



function get_generic_results_pf_sta_red_sol_u(
    pf_sol,
    Ynet_wt_nodes_idx_wt_adjacent_nodes;
    generic_red_sol_kwd_para =
        generic_red_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0,
    wt_branch_current = false )

    #-------------------------------
    # pu

    Ibase =  baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    red_vh_θh = pf_sol.u

    #-------------------------------
    
    if wt_branch_current == true
    
        (;Ybr_cal_and_edges_orientation,
         sta_pf_PQ_para,
         ode_gens_generic_para,
         pf_kw_para
          ) =
              NamedTupleTools.select(
                  generic_red_sol_kwd_para,
                  (:Ybr_cal_and_edges_orientation,
                   :sta_pf_PQ_para,
                   :ode_gens_generic_para,
                   :pf_kw_para))

        (;edges_Ybr_cal,
         edges_orientation) =
             Ybr_cal_and_edges_orientation            
    else

        (;
         # Ybr_cal_and_edges_orientation,
         sta_pf_PQ_para,
         ode_gens_generic_para,
         pf_kw_para
          ) =
              NamedTupleTools.select(
                  generic_red_sol_kwd_para,
                  (# :Ybr_cal_and_edges_orientation,
                   :sta_pf_PQ_para,
                   :ode_gens_generic_para,
                   :pf_kw_para))
            
    end

    P_non_gens =
        sta_pf_PQ_para.P_non_gens
    
    Q_non_gens =
        sta_pf_PQ_para.Q_non_gens
    
    P_g_loc_load =
        sta_pf_PQ_para.P_g_loc_load
    
    Q_g_loc_load =
        sta_pf_PQ_para.Q_g_loc_load
    
    #-------------------------------

    X_d_dash =
        ode_gens_generic_para.X_d_dash

    X_q_dash =
        ode_gens_generic_para.X_q_dash

    X_d =
        ode_gens_generic_para.X_d

    X_q =
        ode_gens_generic_para.X_q
    
    #-------------------------------

    (;loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     # pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs) =
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :pf_kw_gens_vh_slack_θh_para,
              # :pf_kw_net_para,
              :pf_kw_var_idxs,
              :pf_kw_PQ_para_idxs,
              :pf_kw_nodes_types_idxs,
              :pf_kw_n2s_idxs ))
        
    #-------------------------------

    (;slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh) =
         pf_kw_gens_vh_slack_θh_para

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          Ynet_wt_nodes_idx_wt_adjacent_nodes

     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx) =
          pf_kw_var_idxs

     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs) =
          pf_kw_PQ_para_idxs

     (;slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx) =
          pf_kw_nodes_types_idxs    

     (;n2s_slack_gens_idx,
      n2s_non_slack_gens_idx,
      n2s_gens_idx,
      n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx ) =
         pf_kw_n2s_idxs
    
    #-------------------------------

    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------
    
    non_slack_gens_θh =
        red_vh_θh[
            red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[
            red_non_gens_θh_idx2Idx ]

    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ] ]
        for idx in all_nodes_idx ]

    θh_deg = (180/pi) * θh  

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]
    
    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_ynj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx)

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx)
    
    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx )

    Igen = gens_nodes_network_current +
        gens_loc_load_current
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(
            gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    
    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens =
        real.( gens_S )

    pf_Q_g_gens =
        imag.( gens_S )

    #-------------------------------

    Inet = get_Inet_inj(
        uh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    Iinj = get_Iinj(
        uh,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs,    
        loc_load_exist,
        
        Inet)
    
    Sbus_n  = uh .* conj.(
        nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs,    
        loc_load_exist )

    if wt_branch_current == true
    
        Ifrom_Ito = get_I_from_I_to(
            uh, edges_Ybr_cal,
            edges_orientation,
            n2s_all_nodes_idx)

        If = first.( Ifrom_Ito )

        It = last.( Ifrom_Ito )

        Ibranches = If + It

    end
    
    #-------------------------------
    
    Vbus = uh

    #-------------------------------

    # pu_X_d_dash = X_d_dash ./ Zbase
    
    pu_Igen = Igen ./ Ibase
    
    gens_E =
        Vbus[transformed_gens_nodes_idx] .+
        im * X_d_dash .* pu_Igen

    
    gens_mag_E = abs.(gens_E )
    
    gens_ang_E = angle.(gens_E)
     
    #-------------------------------

    gens_vh =
        vh[gens_nodes_idx]
    
    gens_θh =
        θh[gens_nodes_idx]

    #-------------------------------

    if  wt_branch_current == true

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 pf_P_gens, pf_Q_gens,
                 gens_S,

                 pf_P_g_gens,
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 GenSinj,

                 If, It, Ibranches,

                 vh, θh, θh_deg, Vbus, gens_E,
                 
                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 gens_nodes_idx )
        
        
    else

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 
                 pf_P_gens, pf_Q_gens,
                 
                 gens_S,

                 pf_P_g_gens,
                 
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 
                 GenSinj,

                 vh, θh, θh_deg, Vbus, gens_E,

                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 gens_nodes_idx )
        
   end
    
    
end


function get_generic_results_dyn_pf_sol_u(
    pf_sol,
    Ynet_wt_nodes_idx_wt_adjacent_nodes;
    generic_dyn_sol_kwd_para =
        generic_dyn_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0,
    system_status =
        :nothing,
    on_fault_net_para =
        :nothing,
    dyn_pf_vh_θh_id_iq_vhf_θhf_Idx =
        :nothing )

    #-------------------------------
    # pu

    Ibase =  baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    pf_sol_u = pf_sol.u

    #-------------------------------
        
    (;loc_load_exist,
    sta_pf_PQ_para,
    ode_gens_generic_para,
    dyn_pf_flat_vh_flat_θh_id_iq_Idx,
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs) =
         NamedTupleTools.select(
             generic_dyn_sol_kwd_para,
             (:loc_load_exist,
              :sta_pf_PQ_para,
              :ode_gens_generic_para,
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs))

    #-------------------------------

    dyn_pf_fun_kwd_net_idxs =
        deepcopy(dyn_pf_fun_kwd_net_idxs)

    dyn_pf_fun_kwd_n2s_idxs =
        deepcopy(dyn_pf_fun_kwd_n2s_idxs)
    
    #-------------------------------
    
    (P_non_gens,
     Q_non_gens, 
     P_g_loc_load,
     Q_g_loc_load) =
        NamedTupleTools.select(
            sta_pf_PQ_para,
            (:P_non_gens,
             :Q_non_gens,
             :P_g_loc_load,
             :Q_g_loc_load ) )
    
    #-------------------------------

    (; X_d_dash,
     X_q_dash,
     X_d,
     X_q ) =
        NamedTupleTools.select(
            ode_gens_generic_para,
            (:X_d_dash,
             :X_q_dash,
             :X_d,
             :X_q))
    
    #-------------------------------
    
   (;Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        Ynet_wt_nodes_idx_wt_adjacent_nodes
    
    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs ,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))
    
   (;n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))

    pf_kw_nodes_types_idxs =
        (;gens_nodes_idx,
         non_gens_nodes_idx,
         gens_with_loc_load_idx,
         all_nodes_idx)
    
    
    pf_kw_n2s_idxs =
        (;n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx)
    

    (dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
        (:dyn_pf_vh_Idxs,
         :dyn_pf_θh_Idxs,
         :dyn_pf_id_Idxs,
         :dyn_pf_iq_Idxs))
    
    #-------------------------------
    
    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]
    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    
    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
        
    #-------------------------------

    vh      = pf_sol_u[dyn_pf_vh_Idxs]
    θh      = pf_sol_u[dyn_pf_θh_Idxs]
    gens_id = pf_sol_u[dyn_pf_id_Idxs]
    gens_iq = pf_sol_u[dyn_pf_iq_Idxs]

    θh_deg = (180/pi) * θh  

    #-------------------------------

    if system_status == :fault_state
        
        (;dyn_pf_vh_Idxs,
         dyn_pf_θh_Idxs,
         dyn_pf_id_Idxs,
         dyn_pf_iq_Idxs,
         dyn_pf_vhf_Idxs,
         dyn_pf_θhf_Idxs) = 
             NamedTupleTools.select(
                 dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs,
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θhf_Idxs))

        dyn_pf_vh  = pf_sol_u[dyn_pf_vh_Idxs]
        dyn_pf_θh  = pf_sol_u[dyn_pf_θh_Idxs]
        gens_id    = pf_sol_u[dyn_pf_id_Idxs]
        gens_iq    = pf_sol_u[dyn_pf_iq_Idxs]
        dyn_pf_vhf = pf_sol_u[dyn_pf_vhf_Idxs]
        dyn_pf_θhf = pf_sol_u[dyn_pf_θhf_Idxs]

        vh = [dyn_pf_vh; dyn_pf_vhf]
        θh = [dyn_pf_θh; dyn_pf_θhf]

        θh_deg = (180/pi) * θh  

        #-------------------------------

        (Ynet,
         nodes_idx_with_adjacent_nodes_idx,
            
         all_nodes_idx,
         n2s_all_nodes_idx,

         fault_nodes_idx,
         n2s_fault_nodes_idx) =
            NamedTupleTools.select(
                on_fault_net_para,
                (:faulty_Ynet,
                 :faulty_nodes_idx_with_adjacent_nodes_idx,
            
                 :faulty_all_nodes_idx,
                 :n2s_faulty_all_nodes_idx,
                 :fault_nodes_idx,
                 :n2s_fault_nodes_idx))


        non_gens_nodes_idx =
            [non_gens_nodes_idx;
             fault_nodes_idx]

        n2s_non_gens_idx =
            OrderedDict(
                node_idx => idx
                for (idx, node_idx) in
                    enumerate(
                        non_gens_nodes_idx) )
        
        pf_kw_nodes_types_idxs =
            (;gens_nodes_idx,
             non_gens_nodes_idx,
             gens_with_loc_load_idx,
             all_nodes_idx,
             
             fault_nodes_idx)

        pf_kw_n2s_idxs =
            (;n2s_gens_idx,
             n2s_non_gens_idx,
             n2s_gens_with_loc_load_idxs,
             n2s_all_nodes_idx,
             
             n2s_fault_nodes_idx)

        #-------------------------------

        transformed_gens_nodes_idx = [
            n2s_all_nodes_idx[idx]
            for idx in gens_nodes_idx ]

        transformed_non_gens_nodes_idx = [
            n2s_all_nodes_idx[idx]
            for idx in non_gens_nodes_idx ]

        transformed_all_nodes_idx = [
            n2s_all_nodes_idx[idx]
            for idx in all_nodes_idx ]

        #-------------------------------
                
    end        
    
    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]
    
    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_ynj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx)

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx)
    
    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx )

    Igen = gens_nodes_network_current +
        gens_loc_load_current
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(
            gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    
    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens =
        real.( gens_S )

    pf_Q_g_gens =
        imag.( gens_S )

    #-------------------------------

    if system_status == :fault_state
        Inet = get_Inet_inj(
            uh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx )

        Iinj = get_Iinj(
            uh,

            P_non_gens,
            Q_non_gens,

            P_g_loc_load,
            Q_g_loc_load,

            pf_kw_nodes_types_idxs,
            pf_kw_n2s_idxs,    
            loc_load_exist,

            Inet;
            system_status =
                :fault_state)

        Sbus_n  = uh .* conj.(
            nodes_network_current)

        GenSinj = get_GenSinj(
            Sbus_n,

            P_non_gens,
            Q_non_gens,

            P_g_loc_load,
            Q_g_loc_load,

            pf_kw_nodes_types_idxs,
            pf_kw_n2s_idxs,    
            loc_load_exist;
            system_status =
                :fault_state )
        
    else
    
        Inet = get_Inet_inj(
            uh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx )

        Iinj = get_Iinj(
            uh,

            P_non_gens,
            Q_non_gens,

            P_g_loc_load,
            Q_g_loc_load,

            pf_kw_nodes_types_idxs,
            pf_kw_n2s_idxs,    
            loc_load_exist,

            Inet)

        Sbus_n  = uh .* conj.(
            nodes_network_current)

        GenSinj = get_GenSinj(
            Sbus_n,

            P_non_gens,
            Q_non_gens,

            P_g_loc_load,
            Q_g_loc_load,

            pf_kw_nodes_types_idxs,
            pf_kw_n2s_idxs,    
            loc_load_exist )
        
    end

    #-------------------------------
    
    Vbus = uh

    #-------------------------------

    # pu_X_d_dash = X_d_dash ./ Zbase
    
    pu_Igen = Igen ./ Ibase
    
    gens_E =
        Vbus[transformed_gens_nodes_idx] .+
        im * X_d_dash .* pu_Igen

    gens_mag_E = abs.(gens_E)

    gens_ang_E = angle.(gens_E)

    #-------------------------------

    gens_vh =
        vh[gens_nodes_idx]
    
    gens_θh =
        θh[gens_nodes_idx]

    #-------------------------------

    return  (;gens_current_injection,
             gens_loc_load_current,
             gens_nodes_network_current,
             non_gens_nodes_network_current,
             nodes_network_current,

             vh, θh, θh_deg, Vbus,

             Igen, Inet, Iinj, Sbus_n,
             
             pu_Igen, GenSinj,

             S_gens, pf_P_gens, pf_Q_gens,

             gens_S, pf_P_g_gens, pf_Q_g_gens,

             gens_E, gens_mag_E, gens_ang_E,

             gens_vh, gens_θh,
             
             gens_id, gens_iq,

             gens_nodes_idx )    
    
end



# ---------------------------------------------------

function get_uh_Q_vh_θh_and_from_sol(
    sol,
    δ_ω_ed_dash_eq_dash,
    gens_uh_Q_from_red_sol_para )

    red_vh_θh = pf_sol.u

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) = gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         n2s_idxs

    nodes_size =
        sum( length.( [gens_nodes_idx,
                     non_gens_nodes_idx ] ) )
    
    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    # uh = [uh_slack...;
    #   uh_non_slack...;
    #   uh_non_gens... ]

    

    uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    # gens_uh = [ uh[  idx  ]
    #             for idx  in 1:nodes_size
    #                 if idx ∈ gens_nodes_idx]

    gens_uh = uh[gens_nodes_idx]
    
    # ---------------------------------------------------
    
    gens_idq =
        [  get_pf_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            a_ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                a_ra_Xd_dash_Xq_dash ) in
              zip( abs.( gens_uh ),
                   angle.( gens_uh ),
                   δ_ω_ed_dash_eq_dash,
                   ra_Xd_dash_Xq_dash ) ]

    gens_i_d = first.( gens_idq )
        
    gens_i_q = second.( gens_idq )

    gens_idq_flat = [ gens_i_d; gens_i_q ]

    # ---------------------------------------------------

    vec_red_vh_θh_idq =
        [ red_vh_θh, gens_idq_flat ] 

    dim_vec_red_vh_θh_idq =
        length.( vec_red_vh_θh_idq)
    
    _,_, flat_red_vh_θh_idq_Idx =
        create_size_offset_Idx(
            dim_vec_red_vh_θh_idq )
    
    (flat_red_vh_θh_Idxs, flat_idq_Idxs) =
        flat_red_vh_θh_idq_Idx  

    flat_red_vh_θh_idq_Idxs =
        (; flat_red_vh_θh_Idxs,
         flat_idq_Idxs )
    
    red_vh_θh_idq = flat_red_vh_θh_idq =
        [ vec_red_vh_θh_idq...; ]    

    red_vh_θh_idq_and_Idxs =
        (; red_vh_θh_idq,
         flat_red_vh_θh_idq,
         flat_red_vh_θh_idq_Idxs )
    
    # ---------------------------------------------------

    vec_vh_θh_idq_sauer = [ gens_idq_flat,
                            gens_vh,
                            angle.(uh_slack),
                            red_vh_θh  ]

    dim_vec_vh_θh_idq_sauer =
        length.( vec_vh_θh_idq_sauer )
    
    _,_,flat_vh_θh_idq_sauer_Idxs =
        create_size_offset_Idx(
            dim_vec_vh_θh_idq_sauer )
    
    flat_gens_idq_sauer_Idx, flat_gens_vh_sauer_Idx, flat_θh_slack_sauer_Idx, flat_red_vh_θh_sauer_Idx =
        flat_vh_θh_idq_sauer_Idxs
    
    flat_vh_θh_idq_sauer_idxs =
        (; flat_gens_idq_sauer_Idx,
         flat_gens_vh_sauer_Idx,
         flat_θh_slack_sauer_Idx,
         flat_red_vh_θh_sauer_Idx )

    flat_vh_θh_idq_sauer = [ vec_vh_θh_idq_sauer...;]

    sauer_vh_θh_idq_and_idxs =
        (; flat_vh_θh_idq_sauer,
         flat_vh_θh_idq_sauer_idxs )

    # ---------------------------------------------------

    vec_vh_θh_alt =
        [ [abs(a_uh), angle(a_uh)]
          for a_uh in uh  ]

    vh = abs.(uh)

    θh = angle.(uh)
    
    vh_θh = [ vh ;  θh  ]

    vh_θh_alt = [ vec_vh_θh_alt...; ]

    vec_gens_vh_θh = vec_gens_vh_θh_alt =
        [ [abs.( a_uh ),
          angle.( a_uh )]
          for a_uh in
              gens_uh ]
    
    flat_gens_vh_θh =
        [ abs.( gens_uh );
          angle.( gens_uh ) ]

    flat_gens_vh_θh_alt = [ vec_gens_vh_θh_alt...; ]    
      
    # ---------------------------------------------------
    
    # # Sauer alt
    
    # vec_vh_θh_alt_sauer = [
    #     [ abs(a_uh), angle(a_uh) for a_uh in uh_slack ],
    #     [ abs(a_uh), angle(a_uh) for a_uh in uh_non_slack ],
    #     [ abs(a_uh), angle(a_uh) for a_uh in uh_non_gens ] ]
    
    # dim_vec_vh_θh_alt_sauer =
    #     length.( vec_vh_θh_alt_sauer )
    
    # _,_,vec_vh_θh_alt_sauer_Idxs =
    #     create_size_offset_Idx(
    #         dim_vec_vh_θh_alt_sauer )
    
    # flat_vh_θh_alt_sauer_slack_gens_Idx, flat_vh_θh_alt_sauer_non_slack_gens_Idx, flat_vh_θh_alt_sauer_non_gens_Idx =
    #     vec_vh_θh_alt_sauer_Idxs

    # flat_vh_θh_alt_sauer_Idxs =
    #     (; flat_vh_θh_alt_sauer_slack_gens_Idx,
    #      flat_vh_θh_alt_sauer_non_slack_gens_Idx,
    #      flat_vh_θh_alt_sauer_non_gens_Idx )

    # flat_vh_θh_alt_sauer = [ vec_vh_θh_alt_sauer...; ]

    # vh_θh_alt_sauer_and_Idxs =
    #     (; flat_vh_θh_alt_sauer,
    #      flat_vh_θh_alt_sauer_Idxs)

    # ---------------------------------------------------

    idq_net =
        [ get_pf_dyn_idq_θ_π_vhθh(
            vh, θh, δ_ω_ed_eq...,
            a_ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                a_ra_Xd_dash_Xq_dash ) in
              zip( abs.( gens_uh ),
                   angle.( gens_uh ),
                   δ_ω_ed_dash_eq_dash,
                   ra_Xd_dash_Xq_dash ) ]
    
    gens_S  =
        uh[ gens_nodes_idx ] .*
        conj.( x_from_xr_xi.(
            idq_net ) )

    gens_Q = imag.(gens_S)
    
    # ---------------------------------------------------

    gens_δ = first.( δ_ω_ed_dash_eq_dash )

    idq_net_cal =  [ get_pf_dyn_idq_net( id, iq, δ )
                 for (id, iq, δ) in
                     zip( gens_i_d ,
                          gens_i_q, gens_δ ) ]

    
    gens_S_cal  =
        uh[ gens_nodes_idx ] .*
        conj.( x_from_xr_xi.(
            idq_net_cal ) )

    gens_Q_cal = imag.(gens_S_cal)
    
    # ---------------------------------------------------

    #------------------------------------------        
    #------------------------------------------

    return (; uh,
            gens_uh,
            gens_idq,            
            idq_net,
            gens_S, 
            gens_Q,
            vh,
            θh,
            idq_net_cal,
            gens_S_cal,
            gens_Q_cal  )
end
 
# 

function get_integrated_uh_Q_red_vh_θh_and_Idxs_from_sol(
    sol,
    δ_ω_ed_dash_eq_dash,
    gens_uh_Q_from_red_sol_para )

    red_vh_θh = sol.u

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) = gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         n2s_idxs

    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )    

    # uh = [uh_slack...;
    #   uh_non_slack...;
    #   uh_non_gens... ]

    uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    gens_uh = [ uh[  idx  ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]
    
    # ---------------------------------------------------
    
    gens_idq =
        [  get_pf_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( gens_uh ),
                   angle.( gens_uh ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    gens_i_d = first.( gens_idq )
        
    gens_i_q = second.( gens_idq )

    gens_idq_flat = [ gens_i_d; gens_i_q ]

    # ---------------------------------------------------

    vec_red_vh_θh_idq =
        [ red_vh_θh, gens_idq_flat ] 

    dim_vec_red_vh_θh_idq
    
    _,_, flat_red_vh_θh_idq_Idx =
        create_size_offset_Idx(
            dim_vec_red_vh_θh_idq )
    
    flat_red_vh_θh_Idxs, flat_idq_Idxs = flat_red_vh_θh_idq_Idx  

    flat_red_vh_θh_idq_Idxs =
        (; flat_red_vh_θh_Idxs,
         flat_idq_Idxs )
    
    red_vh_θh_idq = flat_red_vh_θh_idq =
        [ vec_red_vh_θh_idq...; ]    

    red_vh_θh_idq_and_Idxs =
        (; red_vh_θh_idq,
         flat_red_vh_θh_idq,
         flat_red_vh_θh_idq_Idxs )
    
    # ---------------------------------------------------

    vec_vh_θh_idq_sauer = [ gens_idq_flat,
                            gens_vh,
                            angle.(uh_slack),
                            red_vh_θh  ]

    dim_vec_vh_θh_idq_sauer =
        length.( vec_vh_θh_idq_sauer )
    
    _,_,flat_vh_θh_idq_sauer_Idxs =
        create_size_offset_Idx(
            dim_vec_vh_θh_idq_sauer )
    
    flat_gens_idq_sauer_Idx, flat_gens_vh_sauer_Idx, flat_θh_slack_sauer_Idx, flat_red_vh_θh_sauer_Idx =
        flat_vh_θh_idq_sauer_Idxs
    
    flat_vh_θh_idq_sauer_idxs =
        (; flat_gens_idq_sauer_Idx,
         flat_gens_vh_sauer_Idx,
         flat_θh_slack_sauer_Idx,
         flat_red_vh_θh_sauer_Idx )

    flat_vh_θh_idq_sauer = [ vec_vh_θh_idq_sauer...;]

    sauer_vh_θh_idq_and_idxs =
        (; flat_vh_θh_idq_sauer,
         flat_vh_θh_idq_sauer_idxs )

    # ---------------------------------------------------

    vec_vh_θh_alt =
        [ [abs(a_uh), angle(a_uh)]
          for a_uh in
              uh  ]

    vh_θh = [ abs.(uh); angle.(uh) ]

    vh_θh_alt = [ vec_vh_θh_alt...; ]

    vec_gens_vh_θh = vec_gens_vh_θh_alt =
        [ [abs.( a_uh ),
          angle.( a_uh )]
          for a_uh in
              gens_uh ]
    
    flat_gens_vh_θh =
        [ abs.( gens_uh );
          angle.( gens_uh ) ]

    flat_gens_vh_θh_alt = [ vec_gens_vh_θh_alt...; ]    
      
    # ---------------------------------------------------
    
    # Sauer alt
    
    vec_vh_θh_alt_sauer = [
        [ abs(a_uh), angle(a_uh) for a_uh in uh_slack ],
        [ abs(a_uh), angle(a_uh) for a_uh in uh_non_slack ],
        [ abs(a_uh), angle(a_uh) for a_uh in uh_non_gens ] ]
    
    dim_vec_vh_θh_alt_sauer =
        length.( vec_vh_θh_alt_sauer )
    
    _,_,vec_vh_θh_alt_sauer_Idxs =
        create_size_offset_Idx(
            dim_vec_vh_θh_alt_sauer )
    
    flat_vh_θh_alt_sauer_slack_gens_Idx, flat_vh_θh_alt_sauer_non_slack_gens_Idx, flat_vh_θh_alt_sauer_non_gens_Idx =
        vec_vh_θh_alt_sauer_Idxs

    flat_vh_θh_alt_sauer_Idxs =
        (; flat_vh_θh_alt_sauer_slack_gens_Idx,
         flat_vh_θh_alt_sauer_non_slack_gens_Idx,
         flat_vh_θh_alt_sauer_non_gens_Idx )

    flat_vh_θh_alt_sauer = [ vec_vh_θh_alt_sauer...; ]

    vh_θh_alt_sauer_and_Idxs =
        (; flat_vh_θh_alt_sauer,
         flat_vh_θh_alt_sauer_Idxs)

    # ---------------------------------------------------

    # idq_net =
    #     [ get_pf_dyn_idq_θ_π_vhθh(
    #         vh, θh, δ_ω_ed_eq...,
    #         ra_Xd_dash_Xq_dash...)
    #       for ( vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in
    #           zip( abs.( gens_uh ),
    #                angle.( gens_uh ),
    #                δ_ω_ed_dash_eq_dash_view,
    #                ra_Xd_dash_Xq_dash_view ) ]


    idq_net =
        [ get_pf_dyn_idq_θ_π_vhθh(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash...)
          for ( vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in
              zip( abs.( gens_uh ),
                   angle.( gens_uh ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    
    gens_S  =
        uh[ gens_nodes_idx ] .*
        conj.( x_from_xr_xi.(
            idq_net ) )

    gens_Q = imag.(gens_S)

    #------------------------------------------

    pf_uh, pf_gens_idq, pf_idq_net = uh, gens_idq, idq_net

    pf_gens_uh, flat_pf_gens_vh_θh = gens_uh, flat_gens_vh_θh

    
    
    #------------------------------------------

    return (; pf_uh,
            uh_slack,
            uh_non_slack,
            uh_non_gens,
            pf_gens_idq,
            pf_idq_net,
            red_vh_θh_idq_and_Idxs,
            sauer_vh_θh_idq_and_idxs,
            vh_θh_alt_sauer_and_Idxs,
            pf_vh_θh,
            vec_vh_θh_alt,
            vh_θh_alt,            
            gens_Q,
            pf_gens_uh,
            flat_pf_gens_vh_θh,
            vec_gens_vh_θh,
            flat_gens_vh_θh_alt )
end



function get_uh_and_gens_Q_from_red_ΔPQ_Δidq_pf_sol(
    sol, gens_uh_Q_from_red_sol_para
     )

    # gens_Q_from_red_pf_para
    
    #------------------------------------------

    # (; slack_vh,
    #  non_slack_gens_Idx_and_vh,
    #  red_vh_Idxs,
    #  red_vh_θh_0_Idx,
    #  idq_0_Idx,
    #  non_gens_θh_idx2Idx,
    #  gens_idx,
    #  n2s_gens_idx,
    #  non_slack_gens_θh_idx2Idx,
    #  δ_ω_ed_dash_eq_dash,
    #  ra_Xd_dash_Xq_dash_view,
    #  gens_id_Idx,
    #  gens_iq_Idx) = gens_Q_from_red_pf_para


    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) =
         NamedTupleTools.select(
             gens_uh_Q_from_red_sol_para,
             (:slack_vh,
              :gens_vh,
              :non_slack_gens_Idx_and_vh,
              :ra_Xd_dash_Xq_dash_view,         
              :red_vh_Idxs,
              :red_θh_Idxs,     
              :red_vh_θh_idx,
              :flat_red_vh_θh_0_Idx,
              :flat_idq_0_Idx,
              :gens_id_Idx,
              :gens_iq_Idx,         
              :non_gens_θh_idx2Idx,
              :non_slack_gens_θh_idx2Idx,         
              :net_idxs,
              :n2s_idxs  ))

    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         NamedTupleTools.select(
             net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :load_nodes_idx,
              :transmission_nodes_idx,
              :non_gens_nodes_idx,
              :load_trans_nodes_idx,
              :gens_with_loc_load_idx))
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             n2s_idxs,
             (:n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_load_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_transmission_idxs,
              :n2s_all_nodes_idx ))


    red_vh_θh_0_idq = sol.u


    red_vh_θh_0 =
        red_vh_θh_0_idq[ flat_red_vh_θh_0_Idx ]

    idq_flat =
        red_vh_θh_0_idq[ flat_idq_0_Idx ]

    #-------------------------------


    uh_slack = [slack_vh * exp(im * 0)]

    non_slack_gens_vh =
        last.( non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh_0[ non_slack_gens_θh_idx2Idx ]

    uh_non_slack =
        non_slack_gens_vh .* exp.(im * non_slack_gens_θh)

    non_gens_vh =
        red_vh_θh_0[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh_0[ non_gens_θh_idx2Idx ]

    uh_non_gens = non_gens_vh .* exp.(im * non_gens_θh )

    # -------------------------------------

    # uh = [uh_slack...;
    #   uh_non_slack...;
    #   uh_non_gens... ]


    uh = [ idx ∈ slack_gens_nodes_idx ?
        uh_slack[n2s_slack_gens_idx[ idx] ] :
        idx ∈ non_slack_gens_nodes_idx ?
        uh_non_slack[ n2s_non_slack_gens_idx[ idx ]] : 
        uh_non_gens[ n2s_non_gens_idx[ idx ] ]
           for idx in 1:nodes_size ]

    gens_uh = [ uh[ idx ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]
    
    # -------------------------------------
    
    gens_i_d = idq_flat[ gens_id_Idx ]

    gens_i_q = idq_flat[ gens_iq_Idx ]


    gens_idq = [ [  i_d, i_q  ]
        for ( i_d, i_q ) in
            zip( gens_i_d, gens_i_q ) ]

    # -------------------------------------

    gens_δ = first.( δ_ω_ed_dash_eq_dash )

    idq_net =  [ get_pf_dyn_idq_net( idq, δ )
                 for (idq, δ) in
                     zip( gens_idq, gens_δ ) ]

    gens_S  =
        uh[gens_idx] .*
        conj.( x_from_xr_xi.(
            idq_net ) )

    # -------------------------------------

    vh_θh_0_idq = [ abs.(uh);
                    angle.(uh);
                    [gens_idq...;] ]

    gens_Q = imag.(gens_S)
    
        return (; uh, vh_θh_0_idq, gens_Q)


end


#---------------------------------------------------
#---------------------------------------------------

function get_vh_θh_from_red_pf_sol_u(
    pf_sol;
    kw_para =
        red_pf_kw_para )

    #-------------------------------   

    red_vh_θh = pf_sol.u

    #-------------------------------

    (;
     gens_vh_slack_θh_para,
     var_idxs,
     nodes_types_idxs,
     n2s_idxs ) =
         kw_para
        
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh) =
         gens_vh_slack_θh_para

     (; red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          var_idxs

     (;slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx) =
          nodes_types_idxs    

     (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
          n2s_idxs
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------
    
    # non_slack_gens_vh = last.(
    #     non_slack_gens_Idx_and_vh )

    non_slack_gens_θh =
        red_vh_θh[
            red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  =
        red_vh_θh[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh[
            red_non_gens_θh_idx2Idx ]

    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in 1:nodes_size ]

    θh = [
        idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ] ]
        for idx in 1:nodes_size ]

    
    #-------------------------------
    
    return (; vh, θh)
    
    #-------------------------------
    
end


function get_generic_results_ds_pf_red_sol_u(
    pf_sol;
    generic_sol_kwd_para =
        generic_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0,
    wt_branch_current = false,
    with_ref_θh  = true )

    #-------------------------------
    # pu
    #-------------------------------

    Ibase =  baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    red_vh_θh_slack_value = pf_sol.u

    #-------------------------------
    
    if wt_branch_current == true
    
        (;Ybr_cal_and_edges_orientation,
         ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para
          ) =
              NamedTupleTools.select(
                  generic_sol_kwd_para,
                  (:Ybr_cal_and_edges_orientation,
                   :ode_gens_para,
                   :sta_pf_PQ_para,
                   :pf_model_kwd_para))

        (;edges_Ybr_cal,
         edges_orientation) =
             Ybr_cal_and_edges_orientation            
    else

        (;ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para
          ) =
              NamedTupleTools.select(
                  generic_sol_kwd_para,
                  (:ode_gens_para,
                   :sta_pf_PQ_para,
                   :pf_model_kwd_para
                   ))
            
    end
    
    #-------------------------------
    #-------------------------------


    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     
     pf_vh_θh_idx_and_idx2Idx,

     gens_loss_participation) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              
              :pf_vh_θh_idx_and_idx2Idx,

              :gens_loss_participation))


    (;participating_gens,
     gens_loss_participation_factor) =
        NamedTupleTools.select(
            gens_loss_participation,
            (:participating_gens,
            :gens_loss_participation_factor))
    
    
     (dyn_pf_flat_vh_flat_θh_Idx,
      non_gens_vh_idx,
      non_slack_gens_θh_idx,
      non_gens_θh_idx,
      
      red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx,

      red_non_gens_vh_Idxs,
      red_non_slack_gens_θh_Idxs,
      red_non_gens_θh_Idxs,

      red_slack_value_Idxs,

      slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      
      gens_nodes_with_loc_loads_idx,
      gens_with_loc_load_idx,
      
      all_nodes_idx,

      n2s_slack_gens_idx,
      n2s_non_slack_gens_idx,
      n2s_gens_idx,
      n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx,

      dyn_pf_fun_kwd_net_idxs,
      dyn_pf_fun_kwd_n2s_idxs,

      non_gens_vh_gens_θh_non_gens_θh_slack_value_Idxs) =
          NamedTupleTools.select(
              pf_vh_θh_idx_and_idx2Idx,
              (:dyn_pf_flat_vh_flat_θh_Idx,
               :non_gens_vh_idx,
               :non_slack_gens_θh_idx,
               :non_gens_θh_idx,
               
               :red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx,

               :red_non_gens_vh_Idxs,
               :red_non_slack_gens_θh_Idxs,
               :red_non_gens_θh_Idxs,

               :red_slack_value_Idxs,

               :slack_gens_nodes_idx,
               :non_slack_gens_nodes_idx,
               :gens_nodes_idx,
               :non_gens_nodes_idx,
               
               :gens_nodes_with_loc_loads_idx,
               :gens_with_loc_load_idx,
               
               :all_nodes_idx,

               :n2s_slack_gens_idx,
               :n2s_non_slack_gens_idx,
               :n2s_gens_idx,
               :n2s_non_gens_idx,
               :n2s_gens_with_loc_load_idxs,
               :n2s_all_nodes_idx,

               :dyn_pf_fun_kwd_net_idxs,
               :dyn_pf_fun_kwd_n2s_idxs,

               :non_gens_vh_gens_θh_non_gens_θh_slack_value_Idxs))


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

    #-------------------------------
    #-------------------------------

    (;P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             (:P_non_gens,
              :Q_non_gens,
              :P_g_loc_load,
              :Q_g_loc_load))
    
    #-------------------------------

    (;X_d_dash,
     X_q_dash,
     X_d,
     X_q ) =
         NamedTupleTools.select(
             ode_gens_para,
             (:X_d_dash,
              :X_q_dash,
              :X_d,
              :X_q))
    
    #-------------------------------
    
     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          Ynet_wt_nodes_idx_wt_adjacent_nodes

    #-------------------------------

    (;slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
             (:slack_gens_vh,
              :slack_gens_θh,

              :gens_vh,
              :non_slack_gens_vh ))
    
    #-------------------------------
    #-------------------------------
    
    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------

    if with_ref_θh == true
    
        non_gens_vh  =
            red_vh_θh_slack_value[ red_non_gens_vh_Idxs ]

        non_slack_gens_θh =
            red_vh_θh_slack_value[
                red_non_slack_gens_θh_Idxs ]

        non_gens_θh =
            red_vh_θh_slack_value[
                red_non_gens_θh_Idxs ]

        slack_value =
            red_vh_θh_slack_value[
                red_slack_value_Idxs]

        #-------------------------------

        vh = [
            idx ∈ gens_nodes_idx ?
                gens_vh[
                    n2s_gens_idx[ idx ] ]  :
                        non_gens_vh[
                            n2s_non_gens_idx[ idx ] ]
                   for idx in all_nodes_idx ]

        θh = [ idx ∈ slack_gens_nodes_idx ?
                slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
                idx ∈ non_slack_gens_nodes_idx ?
                non_slack_gens_θh[
                    n2s_non_slack_gens_idx[ idx] ] :
                        non_gens_θh[
                            n2s_non_gens_idx[ idx ] ]
            for idx in all_nodes_idx ]
        
    else

        (non_gens_vh_Idxs,
         gens_θh_Idxs,       
         non_gens_θh_Idxs,
         slack_value_Idxs) =
             NamedTupleTools.select(
                 non_gens_vh_gens_θh_non_gens_θh_slack_value_Idxs,
                 (:non_gens_vh_Idxs,
                 :gens_θh_Idxs,       
                 :non_gens_θh_Idxs,
                 :slack_value_Idxs)) 

        non_gens_vh =
            non_gens_vh_all_θh_slack_value[
                non_gens_vh_Idxs ]

        gens_θh =
            non_gens_vh_all_θh_slack_value[
                gens_θh_Idxs ]

        non_gens_θh =
            non_gens_vh_all_θh_slack_value[
                non_gens_θh_Idxs ]

        slack_value =
            non_gens_vh_all_θh_slack_value[
                slack_value_Idxs]

        #-------------------------------

        vh = [ idx ∈ gens_nodes_idx ?
                gens_vh[
                    n2s_gens_idx[ idx ]] :
                        non_gens_vh[
                            n2s_non_gens_idx[ idx ]]
                   for idx in all_nodes_idx ]

        θh = [ idx ∈  gens_nodes_idx  ?
                gens_θh[
                    n2s_gens_idx[ idx ]] :
                        non_gens_θh[
                            n2s_non_gens_idx[ idx ]]
            for idx in all_nodes_idx ]

        
    end

    #-------------------------------

    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    θh_deg = (180/pi) * θh  

    #-------------------------------

    gens_vh =
        vh[gens_nodes_idx]
    
    gens_θh =
        θh[gens_nodes_idx]

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]

    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_ynj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx)

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loss_current_contribution =
        gens_loss_participation ./ (
            gens_vh .* cos.(gens_θh ) )
    
    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx)
    
    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx ) +
                gens_loss_current_contribution

    Igen = gens_nodes_network_current +
        gens_loc_load_current +
        gens_loss_current_contribution
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(
            gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    
    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens =
        real.( gens_S )

    pf_Q_g_gens =
        imag.( gens_S )

    #-------------------------------

    Inet = get_Inet_inj(
        uh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    Iinj = get_Iinj(
        uh,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist,
        
        Inet)
    
    Sbus_n  = uh .* conj.(
        nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist )

    if wt_branch_current == true
    
        Ifrom_Ito = get_I_from_I_to(
            uh, edges_Ybr_cal,
            edges_orientation,
            n2s_all_nodes_idx)

        If = first.( Ifrom_Ito )

        It = last.( Ifrom_Ito )

        Ibranches = If + It

    end
    
    #-------------------------------
    
    Vbus = uh

    #-------------------------------

    # pu_X_d_dash = X_d_dash ./ Zbase
    
    pu_Igen = Igen ./ Ibase
    
    gens_E =
        Vbus[transformed_gens_nodes_idx] .+
        im * X_d_dash .* pu_Igen

    
    gens_mag_E = abs.(gens_E )
    
    gens_ang_E = angle.(gens_E)
     
    #-------------------------------

    if  wt_branch_current == true

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 pf_P_gens,
                 pf_Q_gens,
                 gens_S,

                 pf_P_g_gens,
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 GenSinj,

                 If, It, Ibranches,

                 vh, θh, θh_deg, Vbus, gens_E,
                 
                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 slack_value,

                 gens_nodes_idx )
        
        
    else

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 
                 pf_P_gens,
                 pf_Q_gens,
                 
                 gens_S,

                 pf_P_g_gens,
                 
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 
                 GenSinj,

                 vh, θh, θh_deg, Vbus, gens_E,

                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 slack_value,

                 gens_nodes_idx )
        
   end
    
    
end



function get_generic_results_conti_or_ds_pf_red_sol_u(
    pf_sol,
    ds_wt_scale_or_Pg_inj_Png_Qng;
    generic_sol_kwd_para =
        generic_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0,
    wt_branch_current = false,
    with_ref_θh  = true,

    normalise_θh_to_zero_ref = true,

    wt_scaling = false,
    scaling = :nothing, # :multiplicative, :additive
    scaling_type = :nothing # :P , :Q, :P_and_Q
    )

    #-------------------------------
    # pu
    #-------------------------------

    Ibase =  baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    red_vh_θh_slack_value = pf_sol.u

    #-------------------------------
    
    if wt_branch_current == true
    
        (;Ybr_cal_and_edges_orientation,
         ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para
          ) =
              NamedTupleTools.select(
                  generic_sol_kwd_para,
                  (:Ybr_cal_and_edges_orientation,
                   :ode_gens_para,
                   :sta_pf_PQ_para,
                   :pf_model_kwd_para))

        (;edges_Ybr_cal,
         edges_orientation) =
             Ybr_cal_and_edges_orientation            
    else

        (;ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para
          ) =
              NamedTupleTools.select(
                  generic_sol_kwd_para,
                  (:ode_gens_para,
                   :sta_pf_PQ_para,
                   :pf_model_kwd_para
                   ))
            
    end
    
    #-------------------------------
    #-------------------------------


    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     
     Pg_Qg_Png_Qng_Pll_Qll_Idx,

     scale_Pg_Png_Qng_Idx,     
     Pg_Png_Qng_Idx,
     
     pf_vh_θh_idx_and_idx2Idx,

     gens_loss_participation,

     active_power_disturbance_resolution_participation,
     sta_pf_PQ_para ) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,

              :scale_Pg_Png_Qng_Idx,
              :Pg_Png_Qng_Idx,
              
              :pf_vh_θh_idx_and_idx2Idx,

              :gens_loss_participation,

              :active_power_disturbance_resolution_participation,     
              :sta_pf_PQ_para ))

    #-------------------------------
    #-------------------------------
    
     (dyn_pf_flat_vh_flat_θh_Idx,
      non_gens_vh_idx,
      non_slack_gens_θh_idx,
      non_gens_θh_idx,
      
      red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx,

      red_non_gens_vh_Idxs,
      red_non_slack_gens_θh_Idxs,
      red_non_gens_θh_Idxs,

      red_slack_value_Idxs,

      slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      
      gens_nodes_with_loc_loads_idx,
      gens_with_loc_load_idx,
      
      all_nodes_idx,

      n2s_slack_gens_idx,
      n2s_non_slack_gens_idx,
      n2s_gens_idx,
      n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx,
      
      non_gens_vh_gens_θh_non_gens_θh_slack_value_Idxs) =
          NamedTupleTools.select(
              pf_vh_θh_idx_and_idx2Idx,
              (:dyn_pf_flat_vh_flat_θh_Idx,
               :non_gens_vh_idx,
               :non_slack_gens_θh_idx,
               :non_gens_θh_idx,
               
               :red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx,

               :red_non_gens_vh_Idxs,
               :red_non_slack_gens_θh_Idxs,
               :red_non_gens_θh_Idxs,

               :red_slack_value_Idxs,

               :slack_gens_nodes_idx,
               :non_slack_gens_nodes_idx,
               :gens_nodes_idx,
               :non_gens_nodes_idx,
               
               :gens_nodes_with_loc_loads_idx,
               :gens_with_loc_load_idx,
               
               :all_nodes_idx,

               :n2s_slack_gens_idx,
               :n2s_non_slack_gens_idx,
               :n2s_gens_idx,
               :n2s_non_gens_idx,
               :n2s_gens_with_loc_load_idxs,
               :n2s_all_nodes_idx,

               :non_gens_vh_gens_θh_non_gens_θh_slack_value_Idxs))

    #-------------------------------

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

    #-------------------------------        
    #-------------------------------

    (;X_d_dash,
     X_q_dash,
     X_d,
     X_q ) =
         NamedTupleTools.select(
             ode_gens_para,
             (:X_d_dash,
              :X_q_dash,
              :X_d,
              :X_q))
        
    #-------------------------------
    #-------------------------------
        
     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          Ynet_wt_nodes_idx_wt_adjacent_nodes

    #-------------------------------

    (;slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
             (:slack_gens_vh,
              :slack_gens_θh,

              :gens_vh,
              :non_slack_gens_vh ))

    #-------------------------------

    (;participating_gens,
     gens_loss_participation_factor) =
        NamedTupleTools.select(
            gens_loss_participation,
            (:participating_gens,
            :gens_loss_participation_factor))

    (;gens_active_power_particpation_factor,) =
        NamedTupleTools.select(
           active_power_disturbance_resolution_participation,
            (:gens_active_power_particpation_factor,))
    
    #-------------------------------
    #-------------------------------

    
    if wt_scaling == false

        (# Pg_Idxs,
         Png_Idxs,
         Qng_Idxs) =
             NamedTupleTools.select(
                 Pg_Png_Qng_Idx,
                 (# :dyn_P_gens_Idxs,
                  :dyn_P_non_gens_Idxs,
                  :dyn_Q_non_gens_Idxs))

        P_non_gens =
            ds_wt_scale_or_Pg_inj_Png_Qng[
                Png_Idxs]

        Q_non_gens =
            ds_wt_scale_or_Pg_inj_Png_Qng[
                Qng_Idxs]                 
        
    else

        (scale_Idxs,
         # Pg_Idxs,
         Png_Idxs,
         Qng_Idxs) =
              NamedTupleTools.select(
                  scale_Pg_Png_Qng_Idx,
                  (:dyn_scale_Idxs,
                   # :dyn_P_gens_Idxs,
                   :dyn_P_non_gens_Idxs,
                   :dyn_Q_non_gens_Idxs))

        if scaling == :additive

            P_non_gens = scaling_type == :P ?
                (1 .+ scale) .*
                ds_wt_scale_or_Pg_inj_Png_Qng[
                    Png_Idxs] :
                        ds_wt_scale_or_Pg_inj_Png_Qng[
                    Png_Idxs]

            Q_non_gens = scaling_type == :Q ?
                (1 .+ scale) .*
                ds_wt_scale_or_Pg_inj_Png_Qng[
                    Qng_Idxs] :
                        ds_wt_scale_or_Pg_inj_Png_Qng[
                    Qng_Idxs]

            if scaling_type == :P_and_Q

                P_non_gens = (1 .+ scale) .*
                    ds_wt_scale_or_Pg_inj_Png_Qng[
                        Png_Idxs]

                Q_non_gens = (1 .+ scale) .*
                    ds_wt_scale_or_Pg_inj_Png_Qng[
                        Qng_Idxs] 
            end

        elseif scaling == :multiplicative

            P_non_gens = scaling_type == :P ?
                scale .*
                ds_wt_scale_or_Pg_inj_Png_Qng[
                    Png_Idxs] :
                        ds_wt_scale_or_Pg_inj_Png_Qng[
                    Png_Idxs]

            Q_non_gens = scaling_type == :Q ?
                scale .*
                ds_wt_scale_or_Pg_inj_Png_Qng[
                    Qng_Idxs] :
                        ds_wt_scale_or_Pg_inj_Png_Qng[
                        Qng_Idxs]

            if scaling_type == :P_and_Q

                P_non_gens = scale .*
                    ds_wt_scale_or_Pg_inj_Png_Qng[
                    Png_Idxs]

                Q_non_gens = scale .*
                    ds_wt_scale_or_Pg_inj_Png_Qng[
                    Qng_Idxs] 
            end

        else

            P_non_gens =
                ds_wt_scale_or_Pg_inj_Png_Qng[
                Png_Idxs]

            Q_non_gens =
                ds_wt_scale_or_Pg_inj_Png_Qng[
                Qng_Idxs]         
        end
        
        
    end

    
    #-------------------------------
    #-------------------------------

    (Png_base,
     Qng_base,
     P_g_loc_load,
     Q_g_loc_load) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             (:P_non_gens,
              :Q_non_gens,
              :P_g_loc_load,
              :Q_g_loc_load))

    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)

    scaled_P_loading = sum(P_non_gens)

    scaled_Q_loading = sum(Q_non_gens)

    
    ΔP_loading = scaled_P_loading - base_P_loading 
    
    ΔQ_loading = scaled_Q_loading - base_Q_loading 

    gens_active_power_disturbance_particpation =
        ΔP_loading * gens_active_power_particpation_factor

    #-------------------------------
    #-------------------------------
    
    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in
            gens_nodes_idx ]

    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in
            non_gens_nodes_idx ]
    

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in
            slack_gens_nodes_idx ]

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in
            all_nodes_idx ]
    
    transformed_gens_with_loc_load_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_with_loc_load_idx ]

    #-------------------------------

    if with_ref_θh == true
    
        non_gens_vh  =
            red_vh_θh_slack_value[
                red_non_gens_vh_Idxs ]

        non_slack_gens_θh =
            red_vh_θh_slack_value[
                red_non_slack_gens_θh_Idxs ]

        non_gens_θh =
            red_vh_θh_slack_value[
                red_non_gens_θh_Idxs ]

        slack_value =
            red_vh_θh_slack_value[
                red_slack_value_Idxs]

        #-------------------------------

        vh = [
            idx ∈ gens_nodes_idx ?
                gens_vh[
                    n2s_gens_idx[ idx ] ]  :
                        non_gens_vh[
                            n2s_non_gens_idx[ idx ] ]
                   for idx in all_nodes_idx ]

        θh_unzero_ref = θh = [ idx ∈ slack_gens_nodes_idx ?
                slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
                idx ∈ non_slack_gens_nodes_idx ?
                non_slack_gens_θh[
                    n2s_non_slack_gens_idx[ idx] ] :
                        non_gens_θh[
                            n2s_non_gens_idx[ idx ] ]
            for idx in all_nodes_idx ]
        
    else

        (non_gens_vh_Idxs,
         gens_θh_Idxs,       
         non_gens_θh_Idxs,
         slack_value_Idxs) =
             NamedTupleTools.select(
                 non_gens_vh_gens_θh_non_gens_θh_slack_value_Idxs,
                 (:non_gens_vh_Idxs,
                 :gens_θh_Idxs,       
                 :non_gens_θh_Idxs,
                 :slack_value_Idxs)) 

        non_gens_vh =
            non_gens_vh_all_θh_slack_value[
                non_gens_vh_Idxs ]

        gens_θh =
            non_gens_vh_all_θh_slack_value[
                gens_θh_Idxs ]

        non_gens_θh =
            non_gens_vh_all_θh_slack_value[
                non_gens_θh_Idxs ]

        slack_value =
            non_gens_vh_all_θh_slack_value[
                slack_value_Idxs]

        #-------------------------------

        vh = [ idx ∈ gens_nodes_idx ?
                gens_vh[
                    n2s_gens_idx[ idx ]] :
                        non_gens_vh[
                            n2s_non_gens_idx[ idx ]]
                   for idx in all_nodes_idx ]

        θh_unzero_ref = [ idx ∈  gens_nodes_idx  ?
                gens_θh[
                    n2s_gens_idx[ idx ]] :
                        non_gens_θh[
                            n2s_non_gens_idx[ idx ]]
            for idx in all_nodes_idx ]
        
        θh = normalise_θh_to_zero_ref == true ?
            normalise_angle_θh(
                θh_unzero_ref,
                transformed_slack_gens_nodes_idx ) :
                    θh_unzero_ref
        
    end


    #-------------------------------

    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    θh_deg = (180/pi) * θh  

    #-------------------------------

    gens_vh =
        vh[gens_nodes_idx]
    
    gens_θh =
        θh[gens_nodes_idx]

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]

    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_ynj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx)

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loss_current_contribution =
        gens_loss_participation ./ (
            gens_vh .* cos.(gens_θh ) )


    gens_active_current_disturbance_particpation =
        gens_active_power_disturbance_particpation ./ (
            gens_vh .* cos.(gens_θh ) )
    
    #-------------------------------

    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx)
    
    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx ) +
                gens_loss_current_contribution +
                gens_active_current_disturbance_particpation

    Igen = gens_nodes_network_current +
        gens_loc_load_current +
        gens_loss_current_contribution +
        gens_active_current_disturbance_particpation
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(
            gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    
    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens =
        real.( gens_S )

    pf_Q_g_gens =
        imag.( gens_S )

    #-------------------------------

    Inet = get_Inet_inj(
        uh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    Iinj = get_Iinj(
        uh,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist,
        
        Inet)
    
    Sbus_n  = uh .* conj.(
        nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist )

    if wt_branch_current == true
    
        Ifrom_Ito = get_I_from_I_to(
            uh, edges_Ybr_cal,
            edges_orientation,
            n2s_all_nodes_idx)

        If = first.( Ifrom_Ito )

        It = last.( Ifrom_Ito )

        Ibranches = If + It

    end
    
    #-------------------------------
    
    Vbus = uh

    #-------------------------------

    # pu_X_d_dash = X_d_dash ./ Zbase
    
    pu_Igen = Igen ./ Ibase
    
    gens_E =
        Vbus[transformed_gens_nodes_idx] .+
        im * X_d_dash .* pu_Igen

    
    gens_mag_E = abs.(gens_E )
    
    gens_ang_E = angle.(gens_E)
     
    #-------------------------------

    if  wt_branch_current == true

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 pf_P_gens,
                 pf_Q_gens,
                 gens_S,

                 pf_P_g_gens,
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 GenSinj,

                 If, It, Ibranches,

                 vh, θh, θh_deg, Vbus, gens_E,
                 
                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 slack_value,

                 gens_nodes_idx,
    
                 with_ref_θh,
                 wt_scaling,
                 scaling,
                 scaling_type,
                 normalise_θh_to_zero_ref,

                 base_Q_loading,
                 base_P_loading,
                 scaled_P_loading,
                 scaled_Q_loading,

                 θh_unzero_ref,
                 gens_loss_participation,
                 gens_loss_current_contribution,     
                 gens_active_power_disturbance_particpation,
                 gens_active_current_disturbance_particpation)
        
        
    else

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 
                 pf_P_gens,
                 pf_Q_gens,
                 
                 gens_S,

                 pf_P_g_gens,
                 
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 
                 GenSinj,

                 vh, θh, θh_deg, Vbus, gens_E,

                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 slack_value,

                 gens_nodes_idx,
    
                 with_ref_θh,
                 wt_scaling,
                 scaling,
                 scaling_type,
                 normalise_θh_to_zero_ref,

                 base_Q_loading,
                 base_P_loading,
                 scaled_P_loading,
                 scaled_Q_loading,

                 θh_unzero_ref,
                 gens_loss_participation,
                 gens_loss_current_contribution,     
                 gens_active_power_disturbance_particpation,
                 gens_active_current_disturbance_particpation)
                 
        
   end
    
    
end


#--------------------------------------------------
#--------------------------------------------------


function get_generic_results_ds_pf_sol_u(
    pf_sol;
    generic_sol_kwd_para =
        generic_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0,
    wt_branch_current = false,
    normalise_θh_to_zero_ref = true )

    #-------------------------------
    # pu
    #-------------------------------

    Ibase =  baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    vh_θh_slack_value = pf_sol.u

    #-------------------------------
    #-------------------------------
    
    if wt_branch_current == true
    
        (;Ybr_cal_and_edges_orientation,
         ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para) =
              NamedTupleTools.select(
                  generic_sol_kwd_para,
                  (:Ybr_cal_and_edges_orientation,
                   :ode_gens_para,
                   :sta_pf_PQ_para,
                   :pf_model_kwd_para))

        (;edges_Ybr_cal,
         edges_orientation) =
             Ybr_cal_and_edges_orientation            
    else

        (;ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para) =
              NamedTupleTools.select(
                  generic_sol_kwd_para,
                  (:ode_gens_para,
                   :sta_pf_PQ_para,
                   :pf_model_kwd_para ))
            
    end
    
    #-------------------------------
    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
     
     pf_vh_θh_idx_and_idx2Idx,

     gens_loss_participation,
     active_power_disturbance_resolution_participation,
     reactive_power_disturbance_resolution_participation,

     sta_pf_PQ_para) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,
              
              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
              
              :pf_vh_θh_idx_and_idx2Idx,

              :gens_loss_participation,
              :active_power_disturbance_resolution_participation,
              :reactive_power_disturbance_resolution_participation,

              :sta_pf_PQ_para ))

    #-------------------------------
    #-------------------------------

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              Ynet_wt_nodes_idx_wt_adjacent_nodes,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))
    
    #-------------------------------

    (;participating_gens,
     gens_loss_participation_factor) =
        NamedTupleTools.select(
            gens_loss_participation,
            (:participating_gens,
            :gens_loss_participation_factor))

    #-------------------------------
    #-------------------------------
    
   (;n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #-------------------------------

   (;slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_with_loc_load_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))
    
    #-------------------------------

    (vh_Idxs,
     θh_Idxs,
     slack_value_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_slack_value_Idxs))
    
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]

    # normalise_θh_to_zero_ref
    
    θh  = normalise_θh_to_zero_ref == true ?
        normalise_angle_θh(
            vh_θh_slack_value[θh_Idxs],
            slack_gens_nodes_idx ) :
                vh_θh_slack_value[θh_Idxs]
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]
        
    #-------------------------------
    #-------------------------------

    (;P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             (:P_non_gens,
              :Q_non_gens,
              :P_g_loc_load,
              :Q_g_loc_load))
    
    #-------------------------------

    (;X_d_dash,
     X_q_dash,
     X_d,
     X_q ) =
         NamedTupleTools.select(
             ode_gens_para,
             (:X_d_dash,
              :X_q_dash,
              :X_d,
              :X_q))
    
    #-------------------------------
    
    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------
    #-------------------------------

    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    θh_deg = (180/pi) * θh  

    #-------------------------------

    gens_vh =
        vh[gens_nodes_idx]
    
    gens_θh =
        θh[gens_nodes_idx]

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]

    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_ynj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx)

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loss_current_contribution =
        gens_loss_participation ./ (gens_vh .* cos.(
            gens_θh ) )
    
    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx)
    
    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx ) +
                gens_loss_current_contribution

    Igen = gens_nodes_network_current +
        gens_loc_load_current +
        gens_loss_current_contribution
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(
            gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    
    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens =
        real.( gens_S )

    pf_Q_g_gens =
        imag.( gens_S )

    #-------------------------------

    Inet = get_Inet_inj(
        uh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    Iinj = get_Iinj(
        uh,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist,
        
        Inet)
    
    Sbus_n  = uh .* conj.(
        nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist )

    if wt_branch_current == true
    
        Ifrom_Ito = get_I_from_I_to(
            uh, edges_Ybr_cal,
            edges_orientation,
            n2s_all_nodes_idx)

        If = first.( Ifrom_Ito )

        It = last.( Ifrom_Ito )

        Ibranches = If + It

    end
    
    #-------------------------------
    
    Vbus = uh

    #-------------------------------

    # pu_X_d_dash = X_d_dash ./ Zbase
    
    pu_Igen = Igen ./ Ibase
    
    gens_E =
        Vbus[transformed_gens_nodes_idx] .+
        im * X_d_dash .* pu_Igen

    
    gens_mag_E = abs.(gens_E )
    
    gens_ang_E = angle.(gens_E)
     
    #-------------------------------

    if  wt_branch_current == true

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 pf_P_gens,
                 pf_Q_gens,
                 gens_S,

                 pf_P_g_gens,
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 GenSinj,

                 If, It, Ibranches,

                 vh, θh, θh_deg, Vbus, gens_E,
                 
                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 slack_value,

                 gens_nodes_idx )
        
        
    else

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 
                 pf_P_gens,
                 pf_Q_gens,
                 
                 gens_S,

                 pf_P_g_gens,
                 
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 
                 GenSinj,

                 vh, θh, θh_deg, Vbus, gens_E,

                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 slack_value,

                 gens_nodes_idx )
        
   end
    
    
end


function get_generic_results_conti_or_ds_pf_sol_u(
    pf_sol,
    ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng;
    generic_sol_kwd_para =
        generic_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0,
    wt_branch_current = false,
    normalise_θh_to_zero_ref = true,
    
    wt_scaling = false,
    scaling = :nothing, # :multiplicative, :additive
    scaling_type = :nothing # :P , :Q, :P_and_Q
    )

    #-------------------------------
    # pu
    #-------------------------------

    Ibase = baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    vh_θh_slack_value = pf_sol.u

    #-------------------------------
    #-------------------------------
    
    if wt_branch_current == true
    
        (;Ybr_cal_and_edges_orientation,
         ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para) =
              NamedTupleTools.select(
                  generic_sol_kwd_para,
                  (:Ybr_cal_and_edges_orientation,
                   :ode_gens_para,
                   :sta_pf_PQ_para,
                   :pf_model_kwd_para))

        (;edges_Ybr_cal,
         edges_orientation) =
             Ybr_cal_and_edges_orientation            
    else

        (;ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para) =
              NamedTupleTools.select(
                  generic_sol_kwd_para,
                  (:ode_gens_para,
                   :sta_pf_PQ_para,
                   :pf_model_kwd_para ))
            
    end
    
    #-------------------------------
    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
     
     pf_vh_θh_idx_and_idx2Idx,

     gens_loss_participation,
     active_power_disturbance_resolution_participation,
     reactive_power_disturbance_resolution_participation,

     sta_pf_PQ_para) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,
              
              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
              
              :pf_vh_θh_idx_and_idx2Idx,

              :gens_loss_participation,
              :active_power_disturbance_resolution_participation,
              :reactive_power_disturbance_resolution_participation,

              :sta_pf_PQ_para ))


    #-------------------------------
    
   (;n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))

    #-------------------------------

   (;slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_with_loc_load_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))

    #-------------------------------
    #-------------------------------

    (;X_d_dash,
     X_q_dash,
     X_d,
     X_q ) =
         NamedTupleTools.select(
             ode_gens_para,
             (:X_d_dash,
              :X_q_dash,
              :X_d,
              :X_q))
    
    #-------------------------------    
    #-------------------------------

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              Ynet_wt_nodes_idx_wt_adjacent_nodes,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))
    
    #-------------------------------

    (;participating_gens,
     gens_loss_participation_factor) =
        NamedTupleTools.select(
            gens_loss_participation,
            (:participating_gens,
            :gens_loss_participation_factor))


    (;gens_active_power_particpation_factor,) =
        NamedTupleTools.select(
            active_power_disturbance_resolution_participation,
            (:gens_active_power_particpation_factor,))


    (;gens_reactive_power_particpation_factor,) =
        NamedTupleTools.select(
            reactive_power_disturbance_resolution_participation,
            (:gens_reactive_power_particpation_factor,))    
    
    #-------------------------------
    #-------------------------------
    
    if wt_scaling == false

        (# Pg_Idxs,
         # Qg_Idxs,
         Png_Idxs,
         Qng_Idxs) =
             NamedTupleTools.select(
                 Pg_Qg_Png_Qng_Pll_Qll_Idx,
                 (# :dyn_P_gens_Idxs,
                  # :dyn_Q_gens_Idxs,
                  :dyn_P_non_gens_Idxs,
                  :dyn_Q_non_gens_Idxs))

        P_non_gens =
            ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                Png_Idxs]

        Q_non_gens =
            ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                Qng_Idxs]                 
        
    else

        (scale_Idxs,
         # Pg_Idxs,
         # Qg_Idxs,
         Png_Idxs,
         Qng_Idxs) =
              NamedTupleTools.select(
                  scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
                  (:dyn_scale_Idxs,
                   # :dyn_P_gens_Idxs,
                   # :dyn_Q_gens_Idxs,
                   :dyn_P_non_gens_Idxs,
                   :dyn_Q_non_gens_Idxs))

        if scaling == :additive

            P_non_gens = scaling_type == :P ?
                (1 .+ scale) .*
                ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                    Png_Idxs] :
                        ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                    Png_Idxs]

            Q_non_gens = scaling_type == :Q ?
                (1 .+ scale) .*
                ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                    Qng_Idxs] :
                        ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                    Qng_Idxs]

            if scaling_type == :P_and_Q

                P_non_gens = (1 .+ scale) .*
                    ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                        Png_Idxs]

                Q_non_gens = (1 .+ scale) .*
                    ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                        Qng_Idxs] 
            end

        elseif scaling == :multiplicative

            P_non_gens = scaling_type == :P ?
                scale .*
                ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                    Png_Idxs] :
                        ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                    Png_Idxs]

            Q_non_gens = scaling_type == :Q ?
                scale .*
                ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                    Qng_Idxs] :
                        ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                        Qng_Idxs]

            if scaling_type == :P_and_Q

                P_non_gens = scale .*
                    ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                    Png_Idxs]

                Q_non_gens = scale .*
                    ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                    Qng_Idxs] 
            end

        else

            P_non_gens =
                ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                Png_Idxs]

            Q_non_gens =
                ds_wt_scale_or_Pg_inj_Qn_inj_Png_Qng[
                Qng_Idxs]         
        end
        
        
    end
    
    #-------------------------------

    (Png_base,
     Qng_base,
     P_g_loc_load,
     Q_g_loc_load) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             (:P_non_gens,
              :Q_non_gens,
              :P_g_loc_load,
              :Q_g_loc_load))


    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)

    scaled_P_loading = sum(P_non_gens)

    scaled_Q_loading = sum(Q_non_gens)

    
    ΔP_loading = scaled_P_loading - base_P_loading 
    
    ΔQ_loading = scaled_Q_loading - base_Q_loading 

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
        ΔQ_loading * gens_reactive_power_particpation_factor

    gens_power_disturbance_particpation =
        gens_active_power_particpation +
        im * gens_reactive_power_particpation
    
    #-------------------------------

    (vh_Idxs,
     θh_Idxs,
     slack_value_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_slack_value_Idxs))
    
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]

    # normalise_θh_to_zero_ref

    θh_unzero_ref  =
        vh_θh_slack_value[θh_Idxs]
    
    θh  = normalise_θh_to_zero_ref == true ?
        normalise_angle_θh(
            θh_unzero_ref,
            slack_gens_nodes_idx ) :
                θh_unzero_ref 
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]

    #-------------------------------

    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    
    θh_deg = (180/pi) * θh  
    
    #-------------------------------
    
    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]

    
    transformed_gens_with_loc_load_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_with_loc_load_idx ]
    
    #-------------------------------
    #-------------------------------

    gens_vh =
        vh[transformed_gens_nodes_idx]
    
    gens_θh =
        θh[transformed_gens_nodes_idx]

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]

    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_ynj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx)

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loss_current_contribution =
        gens_loss_participation ./ (gens_vh .* cos.(
            gens_θh ) )

    gens_current_disturbance_particpation =
        conj.(gens_power_disturbance_particpation) ./conj.(
            gens_uh)
    
    #-------------------------------
    
    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx)
    
    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_with_loc_load_idx,
            gens_nodes_idx ) +
                gens_loss_current_contribution +
                gens_current_disturbance_particpation

    Igen = gens_nodes_network_current +
        gens_loc_load_current +
        gens_loss_current_contribution +
        gens_current_disturbance_particpation
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(
            gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    
    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens =
        real.( gens_S )

    pf_Q_g_gens =
        imag.( gens_S )

    #-------------------------------

    Inet = get_Inet_inj(
        uh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    Iinj = get_Iinj(
        uh,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist,
        
        Inet)
    
    Sbus_n  = uh .* conj.(
        nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist )

    if wt_branch_current == true
    
        Ifrom_Ito = get_I_from_I_to(
            uh, edges_Ybr_cal,
            edges_orientation,
            n2s_all_nodes_idx)

        If = first.( Ifrom_Ito )

        It = last.( Ifrom_Ito )

        Ibranches = If + It

    end
    
    #-------------------------------
    
    Vbus = uh

    #-------------------------------

    # pu_X_d_dash = X_d_dash ./ Zbase
    
    pu_Igen = Igen ./ Ibase
    
    gens_E =
        Vbus[transformed_gens_nodes_idx] .+
        im * X_d_dash .* pu_Igen

    
    gens_mag_E = abs.(gens_E )
    
    gens_ang_E = angle.(gens_E)
     
    #-------------------------------

    if  wt_branch_current == true

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 pf_P_gens,
                 pf_Q_gens,
                 gens_S,

                 pf_P_g_gens,
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 GenSinj,

                 If, It, Ibranches,

                 vh, θh, θh_deg, Vbus, gens_E,
                 
                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 slack_value,

                 gens_nodes_idx,
                 
                 wt_scaling,
                 scaling,
                 scaling_type,
                 normalise_θh_to_zero_ref,

                 base_P_loading,
                 base_Q_loading,
                 scaled_P_loading,
                 scaled_Q_loading,
                  
                 θh_unzero_ref,
                 gens_loss_participation,
                 gens_loss_current_contribution,
                 gens_power_disturbance_particpation,
                 gens_current_disturbance_particpation)
        
        
    else

        return  (;gens_current_injection,
                 gens_loc_load_current,
                 gens_nodes_network_current,
                 non_gens_nodes_network_current,
                 nodes_network_current,

                 S_gens,
                 
                 pf_P_gens,
                 pf_Q_gens,
                 
                 gens_S,

                 pf_P_g_gens,
                 
                 pf_Q_g_gens,

                 Igen, Inet, Iinj, pu_Igen,

                 Sbus_n,
                 
                 GenSinj,

                 vh, θh, θh_deg, Vbus, gens_E,

                 gens_mag_E, gens_ang_E,

                 gens_vh, gens_θh,

                 slack_value,

                 gens_nodes_idx,
                 
                 wt_scaling,
                 scaling,
                 scaling_type,
                 normalise_θh_to_zero_ref,

                 base_P_loading,
                 base_Q_loading,
                 scaled_P_loading,
                 scaled_Q_loading,
                 
                 θh_unzero_ref,
                 gens_loss_participation,
                 gens_loss_current_contribution,
                 gens_power_disturbance_particpation,
                 gens_current_disturbance_particpation)
                 
        
   end
        
end


function get_generic_results_ds_comparray_pf_sol_u(
    pf_sol,
    ds_pf_setpoint_model_para_wt_Ynet_flattend_comparray;
    generic_sol_kwd_para =
        generic_sol_kwd_para,
    baseMVA = 1.0,
    basekV = 1.0,
    Pg_inj_Qg_inj_bool = false,
    normalise_θh_to_zero_ref = true,
    
    wt_scaling = false,
    scaling = :nothing, # :multiplicative, :additive
    scaling_type = :nothing # :P , :Q, :P_and_Q
    )
    
    # -------------------------------------

    if Pg_inj_Qg_inj_bool == true
        
        @unpack Pg_net_inj, Qg_net_inj, P_non_gens, Q_non_gens, Ynet_real_imag_flattend =
            ds_pf_setpoint_model_para_wt_Ynet_flattend_comparray

        P_gens = Pg_net_inj

        Q_gens = Qg_net_inj
    else
        
        @unpack P_gens, Q_gens, P_non_gens, Q_non_gens, P_g_loc_load, Q_g_loc_load, Ynet_real_imag_flattend =
            ds_pf_setpoint_model_para_wt_Ynet_flattend_comparray
        
    end
    
    
    
    #-------------------------------
    # pu
    #-------------------------------

    Ibase = baseMVA / basekV

    Zbase = basekV / Ibase

    Ybase = Ibase / basekV
    
    #-------------------------------    
        
    vh_θh_slack_value = pf_sol.u

    #-------------------------------
    #-------------------------------
    
    # if wt_branch_current == true

    #     (;edges_Ybr_cal,
    #      edges_orientation) =
    #        get_edges_Ybr_cal_and_edges_orientation_by_generic(
    #              edges_fbus, edges_tbus,
    #              edges_r, edges_x, edges_b,
    #              edges_ratio, edges_angle, Gs, Bs;
    #              baseMVA = baseMVA, basekV = basekV )
            
    # end


   (;ode_gens_para,
    sta_pf_PQ_para,
    pf_model_kwd_para) =
         NamedTupleTools.select(
             generic_sol_kwd_para,
             (:ode_gens_para,
              :sta_pf_PQ_para,
              :pf_model_kwd_para ))
    
    #-------------------------------
    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
     
     pf_vh_θh_idx_and_idx2Idx,

     Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend,

     gens_loss_participation,
     active_power_disturbance_resolution_participation,
     reactive_power_disturbance_resolution_participation,

     sta_pf_PQ_para) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,
              
              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
              
              :pf_vh_θh_idx_and_idx2Idx,

              :Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend,

              :gens_loss_participation,
          :active_power_disturbance_resolution_participation,
          :reactive_power_disturbance_resolution_participation,

              :sta_pf_PQ_para ))


    #-------------------------------
    
   (;n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))

    #-------------------------------

   (;slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_with_loc_load_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))

    #-------------------------------
    #-------------------------------

    (;X_d_dash,
     X_q_dash,
     X_d,
     X_q ) =
         NamedTupleTools.select(
             ode_gens_para,
             (:X_d_dash,
              :X_q_dash,
              :X_d,
              :X_q))
    
    #-------------------------------    
    #-------------------------------

     (; # Ynet,
      nodes_idx_with_adjacent_nodes_idx,) =
          NamedTupleTools.select(
              Ynet_wt_nodes_idx_wt_adjacent_nodes,
              (# :Ynet,
               :nodes_idx_with_adjacent_nodes_idx, ))

    (Ynet_real_Idxs,
     Ynet_imag_Idxs ) =
         Ynet_real_imag_Idxs_in_flattend

    Ynet_real =
        Ynet_real_imag_flattend[ Ynet_real_Idxs ]

    Ynet_imag =
        Ynet_real_imag_flattend[ Ynet_imag_Idxs ]


    Ynet = [Ynet_real[idx] + im * Ynet_imag[idx]  
              for idx in Ynet_rows_Idxs_in_flattend ]
    
    #-------------------------------

    (;participating_gens,
     gens_loss_participation_factor) =
        NamedTupleTools.select(
            gens_loss_participation,
            (:participating_gens,
            :gens_loss_participation_factor))


    (;gens_active_power_particpation_factor,) =
        NamedTupleTools.select(
            active_power_disturbance_resolution_participation,
            (:gens_active_power_particpation_factor,))


    (;gens_reactive_power_particpation_factor,) =
        NamedTupleTools.select(
            reactive_power_disturbance_resolution_participation,
            (:gens_reactive_power_particpation_factor,))    
    
    #-------------------------------
    #-------------------------------

    if scaling == :additive

        P_non_gens = scaling_type == :P ?
            (1 .+ scale) .* P_non_gens : P_non_gens

        Q_non_gens = scaling_type == :Q ?
            (1 .+ scale) .* Q_non_gens : Q_non_gens

        if scaling_type == :P_and_Q

            P_non_gens = (1 .+ scale) .* P_non_gens

            Q_non_gens = (1 .+ scale) .* Q_non_gens

        end

    elseif scaling == :multiplicative

        P_non_gens = scaling_type == :P ?
            scale .* P_non_gens : P_non_gens

        Q_non_gens = scaling_type == :Q ?
            scale .* Q_non_gens : Q_non_gens

        if scaling_type == :P_and_Q

            P_non_gens = scale .* P_non_gens

            Q_non_gens = scale .* Q_non_gens
        end

    else

        nothing
    end
    
    #-------------------------------

    (Png_base,
     Qng_base,
     P_g_loc_load,
     Q_g_loc_load) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             (:P_non_gens,
              :Q_non_gens,
              :P_g_loc_load,
              :Q_g_loc_load))


    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)

    scaled_P_loading = sum(P_non_gens)

    scaled_Q_loading = sum(Q_non_gens)

    
    ΔP_loading = scaled_P_loading - base_P_loading 
    
    ΔQ_loading = scaled_Q_loading - base_Q_loading 

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
        ΔQ_loading * gens_reactive_power_particpation_factor

    gens_power_disturbance_particpation =
        gens_active_power_particpation +
        im * gens_reactive_power_particpation
    
    #-------------------------------

    (vh_Idxs,
     θh_Idxs,
     slack_value_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_slack_value_Idxs))
    
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]

    # normalise_θh_to_zero_ref

    θh_unzero_ref  =
        vh_θh_slack_value[θh_Idxs]
    
    θh  = normalise_θh_to_zero_ref == true ?
        normalise_angle_θh(
            θh_unzero_ref,
            slack_gens_nodes_idx ) :
                θh_unzero_ref 
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]

    #-------------------------------

    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    
    θh_deg = (180/pi) * θh  
    
    #-------------------------------
    
    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]

    
    transformed_gens_with_loc_load_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_with_loc_load_idx ]
    
    #-------------------------------
    #-------------------------------

    gens_vh =
        vh[transformed_gens_nodes_idx]
    
    gens_θh =
        θh[transformed_gens_nodes_idx]

    #-------------------------------
    
    uh = vh .* exp.(im * θh)

    gens_uh = uh[ gens_nodes_idx ]

    #-------------------------------
    
    nodes_network_current =
        get_nodes_∑_ynj_x_vj(
            vh,
            θh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx)

    gens_nodes_network_current =
        nodes_network_current[
            transformed_gens_nodes_idx ]

    non_gens_nodes_network_current =
        nodes_network_current[
            transformed_non_gens_nodes_idx ]

    #-------------------------------

    gens_loss_current_contribution =
        gens_loss_participation ./ (gens_vh .* cos.(
            gens_θh ) )

    gens_current_disturbance_particpation =
        conj.(gens_power_disturbance_particpation) ./conj.(
            gens_uh)
    
    #-------------------------------
    
    gens_loc_load_current =
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,
                        
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,
            
            gens_nodes_idx,
            gens_with_loc_load_idx)
    
    gens_current_injection =
        get_gens_current_injection(
            vh,
            θh,           
            Ynet,
            
            nodes_idx_with_adjacent_nodes_idx,

            P_g_loc_load,
            Q_g_loc_load,
            
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, 

            gens_with_loc_load_idx,
            gens_nodes_idx ) +
                gens_loss_current_contribution +
                gens_current_disturbance_particpation

    Igen = gens_nodes_network_current +
        gens_loc_load_current +
        gens_loss_current_contribution +
        gens_current_disturbance_particpation
    
    #-------------------------------

    S_gens =
        gens_uh .* conj.(
            gens_current_injection)

    pf_P_gens = real.( S_gens )

    pf_Q_gens = imag.( S_gens )    
    
    #-------------------------------

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.( Igen )

    
    pf_P_g_gens =
        real.( gens_S )

    pf_Q_g_gens =
        imag.( gens_S )

    #-------------------------------

    Inet = get_Inet_inj(
        uh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    Iinj = get_Iinj(
        uh,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist,
        
        Inet)
    
    Sbus_n  = uh .* conj.(
        nodes_network_current)

    GenSinj = get_GenSinj(
        Sbus_n,
        
        P_non_gens,
        Q_non_gens,
        
        P_g_loc_load,
        Q_g_loc_load,
        
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs,
        
        loc_load_exist )

    
    #-------------------------------
    
    Vbus = uh

    #-------------------------------

    # pu_X_d_dash = X_d_dash ./ Zbase
    
    pu_Igen = Igen ./ Ibase
    
    gens_E =
        Vbus[transformed_gens_nodes_idx] .+
        im * X_d_dash .* pu_Igen

    
    gens_mag_E = abs.(gens_E )
    
    gens_ang_E = angle.(gens_E)
     
    #-------------------------------

    return  (;gens_current_injection,
             gens_loc_load_current,
             gens_nodes_network_current,
             non_gens_nodes_network_current,
             nodes_network_current,

             S_gens,

             pf_P_gens,
             pf_Q_gens,

             gens_S,

             pf_P_g_gens,

             pf_Q_g_gens,

             Igen, Inet, Iinj, pu_Igen,

             Sbus_n,

             GenSinj,

             vh, θh, θh_deg, Vbus, gens_E,

             gens_mag_E, gens_ang_E,

             gens_vh, gens_θh,

             slack_value,

             gens_nodes_idx,

             wt_scaling,
             scaling,
             scaling_type,
             normalise_θh_to_zero_ref,

             base_P_loading,
             base_Q_loading,
             scaled_P_loading,
             scaled_Q_loading,

             θh_unzero_ref,
             gens_loss_participation,
             gens_loss_current_contribution,
             gens_power_disturbance_particpation,
             gens_current_disturbance_particpation)
        
end

#--------------------------------------------
#--------------------------------------------
