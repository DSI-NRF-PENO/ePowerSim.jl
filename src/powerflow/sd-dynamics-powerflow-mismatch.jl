# (c) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.4:  pg 135 - 136, 6.121 - 6.123


########################################################
# ------------------------------------------------------
#  Jacobian
# ------------------------------------------------------
########################################################

function sta_pf_Jac!_df_dx_by_Ynet_or_Yπ_net!(
    Jac_vh_θh,
    red_vh_θh_x,
    pf_PQ_param ;
    pf_kw_para,
    
    # This could be any mismatech func for Ynet_or_Yπ_net
    # e.g get_generic_sta_pf_ΔPQ_mismatch

    # get_ΔPQ_mismatch_by_Ynet    
    # get_ΔI_mismatch_by_Yπ_net
    # get_ΔPQ_mismatch_by_Ynet
    # get_ΔI_mismatch_by_Ynet    
    func = get_ΔPQ_mismatch_by_Ynet,
    
    # Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
    # Ynet_wt_nodes_idx_wt_adjacent_nodes,    
    net_addmitance_tuple =
        Ynet_wt_nodes_idx_wt_adjacent_nodes,
        
    # :generic_Ynet, :Ynet or :Yπ_net
    by_Ynet_or_Yπ_net =
        :Ynet )
    
    #-------------------------------

    ΔPQ = similar(red_vh_θh_x)
    
    if by_Ynet_or_Yπ_net == :generic_Ynet
        
        # return
        Jac_vh_θh .= ForwardDiff.jacobian(
        ( g, x ) ->
            func(
                g, x, pf_PQ_param
                ;pf_kw_para =
                    pf_kw_para ),
            ΔPQ,
            red_vh_θh_x )

    elseif by_Ynet_or_Yπ_net == :Ynet 
        Ynet_wt_nodes_idx_wt_adjacent_nodes =
            net_addmitance_tuple
        # return
        Jac_vh_θh .= ForwardDiff.jacobian(
        ( g, x ) ->
            func(
                g, x, pf_PQ_param
                ;pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
                Ynet_wt_nodes_idx_wt_adjacent_nodes =
                    Ynet_wt_nodes_idx_wt_adjacent_nodes),
            ΔPQ,
            red_vh_θh_x )
        
    elseif by_Ynet_or_Yπ_net == :Yπ_net

        Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
            net_addmitance_tuple
        #return
         Jac_vh_θh .= ForwardDiff.jacobian(
        ( g, x ) ->
            func(
                g, x, pf_PQ_param
                ;pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
                Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
                    Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes ),
            ΔPQ,
            red_vh_θh_x )

    else

        nothing
        
    end
    

    return nothing

end

# ------------------------------------------------

function sta_pf_Jac!_df_dp_by_Ynet_or_Yπ_net!(
    Jac_vh_θh,
    red_vh_θh_x,
    pf_PQ_param ;
    pf_kw_para,
    
    # This could be any mismatech func for Ynet_or_Yπ_net
    # e.g get_generic_sta_pf_ΔPQ_mismatch

    # get_ΔPQ_mismatch_by_Ynet    
    # get_ΔI_mismatch_by_Yπ_net
    # get_ΔPQ_mismatch_by_Ynet
    # get_ΔI_mismatch_by_Ynet    
    func = get_ΔPQ_mismatch_by_Ynet,
    
    # Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
    # Ynet_wt_nodes_idx_wt_adjacent_nodes,    
    net_addmitance_tuple =
        Ynet_wt_nodes_idx_wt_adjacent_nodes,
        
    # :generic_Ynet, :Ynet or :Yπ_net
    by_Ynet_or_Yπ_net =
        :Ynet )
    
    #-------------------------------

    ΔPQ = similar(red_vh_θh_x)

    if by_Ynet_or_Yπ_net == :generic_Ynet
        
        Jac_vh_θh .= ForwardDiff.jacobian(
        ( g , p ) ->
            func(
                g, red_vh_θh_x, p
                ;pf_kw_para = pf_kw_para ),
             ΔPQ ,
            pf_PQ_param )

    elseif by_Ynet_or_Yπ_net == :Ynet

        Ynet_wt_nodes_idx_wt_adjacent_nodes =
            net_addmitance_tuple
        
        Jac_vh_θh .= ForwardDiff.jacobian(
        ( g , p ) ->
            func(
                g, red_vh_θh_x, p
                ;pf_kw_para = disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
                Ynet_wt_nodes_idx_wt_adjacent_nodes =
                    Ynet_wt_nodes_idx_wt_adjacent_nodes),
             ΔPQ ,
            pf_PQ_param )
        
    elseif by_Ynet_or_Yπ_net == :Yπ_net
        
        Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
            net_addmitance_tuple
        
        Jac_vh_θh .= ForwardDiff.jacobian(
        ( ΔPQ, p ) ->
            func(
                f, red_vh_θh_x, p
                ;pf_kw_para =
                    disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para),
                Yπ_net_wt_Yshunt_wt_nodes_idx =
                    Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes  ),
            ΔPQ,
            pf_PQ_param )

    else

        nothing
        
    end
    return nothing
    
end

# ------------------------------------------------

function sta_pf_Jac!(
    Jac,
    red_vh_θh_x,
    pf_PQ_param ;
    Ybus = Ybus,
    pf_kw_para )

    # ------------------------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     transformed_all_nodes_idx) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              :transformed_all_nodes_idx) )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]

    # -------------------------------------------------

    uh     = vh .* exp.(im * θh)
    
    ν      = abs.( uh )
    
    Λ      = inv.(uh)

    Λ_conj = inv.( conj.(uh) )
        
    Θ      = angle.(uh)

    E      = Diagonal( inv.(ν) ) * uh     

    Ibus   = Ybus * uh
    
    nodes_size = length( Ibus )
    
    # ------------------------------------------------    
    # Jacobian: See Matpower Technical Note 2, page 8,
    # eq 29 and 32
    
    Gc_vh = Diagonal(uh) * ( Diagonal(conj.(Ibus)) +
        conj(Ybus) * Diagonal( conj.(uh) ) ) *
        Diagonal(inv.(ν))
    
    Gc_θh = im * Diagonal(uh) * (Diagonal(conj.(Ibus)) -
        conj(Ybus) * Diagonal( conj.(uh) ) )
    
    P_Gc_vh = real.(
        Gc_vh[setdiff(1:end,
                      transformed_slack_gens_nodes_idx),
              setdiff(1:end,
                      transformed_gens_nodes_idx)] )
    
    P_Gc_θh = real.(
        Gc_θh[setdiff(1:end,
                      transformed_slack_gens_nodes_idx),
              setdiff(1:end,
                      transformed_slack_gens_nodes_idx)] )
    
    Q_Gc_vh = imag.(
        Gc_vh[setdiff(1:end,
                      transformed_gens_nodes_idx),
              setdiff(1:end,
                      transformed_gens_nodes_idx) ] )
    
    Q_Gc_θh = imag.(
        Gc_θh[setdiff(1:end,
                      transformed_gens_nodes_idx),
              setdiff(1:end,
                      transformed_slack_gens_nodes_idx)] )
    
    Jac  .= vcat((hcat(P_Gc_vh, P_Gc_θh )),
                 (hcat(Q_Gc_vh, Q_Gc_θh)))

    return nothing
end



########################################################
# ------------------------------------------------------
#  Powerflow mismatch
# ------------------------------------------------------
########################################################

#-----------------------------------------------------
# new
#-----------------------------------------------------

function dynamic_power_balance_mismatch_with_vh_θh!(
    red_ΔPQ,
    red_vh_θh_0_view,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view ),
    pf_net_param;
    with_δ_ed_eq = true )
    
    # ---------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal, edges_orientation) =
         pf_net
    
    (slack_vh, gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
        pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
     nodes_u_Idx,
     gens_idx, ur_IDX,
     ui_IDX, vh_IDX,
     θh_IDX, red_vh_θh_idx,
     ur_idx, ui_idx, ur_ui_idx) =
         pf_Idx

    
    # ----------------------------------------------

    working_vh_θh_view[ red_vh_θh_idx ] .=
        red_vh_θh_0_view

    uh =  working_vh_θh_view[ vh_IDX ] .*
        exp.(im * working_vh_θh_view[ θh_IDX ])
        
    # ---------------------------------------------- 

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    # S_gens_loc_load = x_from_xr_xi.( P_Q_gens_loc_load_view )
    # S_non_gens = x_from_xr_xi.( P_Q_non_gens_view )

    # --------------------------------------------- 

    if with_δ_ed_eq == true
        
        idq_θ_π_vhθh = [
            get_dynamic_idq_θ_π_vhθh(
                vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
            for (vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash) in
                zip( abs.(uh), angle.(uh),
                     δ_ω_ed_dash_eq_dash_view ,
                     ra_Xd_dash_Xq_dash_view)]

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh, Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx)

        current_mismatch = get_nodes_current_mismatch(
            uh, (P_Q_non_gens_view,
                 P_Q_gens_loc_load_view),
            idq_θ_π_vhθh, I_sum_ynj_vj )
        
    else

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh, Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )
    
        current_mismatch = get_nodes_current_mismatch(
            uh, (P_Q_gens_view, P_Q_non_gens_view,
                 P_Q_gens_loc_load_view), I_sum_ynj_vj )

    end
        
    # ----------------------------------------------
    
    power_mismatch = uh .* conj.(current_mismatch)

    red_P_mismatch = real.(
        power_mismatch[setdiff(1:end, slack_bus_idx)])
    
    red_Q_mismatch = imag.(
        power_mismatch[setdiff(1:end,gens_idx)])
        
    red_ΔPQ .= vcat(red_P_mismatch, red_Q_mismatch)

   # ------------------------------------------------
    # update Inet_view and Iinj_view
    # -----------------------------------------------

    Inet_view .= I_sum_ynj_vj

    Iinj_view .= get_Iinj(
        uh, P_Q_non_gens_view,
        P_Q_gens_loc_load_view,
        Inet_view )
    
    # ------------------------------------------------   
    # update Q 
    # ------------------------------------------------

    if with_δ_ed_eq == true
        
        Igen = idq_θ_π_vhθh
        
        gens_S  = uh[gens_idx] .* conj.(Igen[gens_idx])
        
        P_Q_nodes_view[slack_bus_idx] .=
            [real(gens_S[slack_bus_idx]),
             imag(gens_S[slack_bus_idx])]

        for (idx, gen_S) in zip(gens_idx, gens_S)

            P_Q_gens_view[idx] .=
                [real( S_gens[idx]), imag(gen_S)]
       end
    else
  
        Igen = I_sum_ynj_vj +
            get_gens_local_load_current(
                uh, P_Q_gens_loc_load_view) 

        gens_S  = uh[gens_idx] .* conj.(Igen[gens_idx])

        P_Q_nodes_view[slack_bus_idx] .= [
            real(gens_S[slack_bus_idx]),
            imag(gens_S[slack_bus_idx])]

        for (idx, gen_S) in zip(gens_idx, gens_S)

            P_Q_gens_view[idx] .= [
                real( S_gens[idx]), imag(gen_S)]
        end
    end
    
    # ------------------------------------------------ 
    # update nodes_pf_U_view
    # ------------------------------------------------ 

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .=
            [ real( uk ), imag( uk ) ]

    end

    return nothing
    
end



function power_balance_mismatch_with_vh_θh!(
    red_ΔPQ,
    red_vh_θh_0_view,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view),
      pf_net_param;
      with_δ_ed_eq = true )

    # ------------------------------------------------
    # ------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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

    
    # ------------------------------------------------

    working_vh_θh_view[ red_vh_θh_idx ] .=
        red_vh_θh_0_view

    uh =  working_vh_θh_view[ vh_IDX ] .* exp.(
        im * working_vh_θh_view[ θh_IDX ])
        
    # ------------------------------------------------ 

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    # ------------------------------------------------- 

    if with_δ_ed_eq == true
        
        idq_θ_π_vhθh = [
            get_dynamic_idq_θ_π_vhθh(
                vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
            for (vh,θh,δ_ω_ed_eq,ra_Xd_dash_Xq_dash) in
                zip( abs.(uh),
                     angle.(uh),
                     δ_ω_ed_dash_eq_dash_view ,
                     ra_Xd_dash_Xq_dash_view ) ]

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh,
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )

        current_mismatch =
            get_nodes_current_mismatch(
                uh,
                (P_Q_non_gens_view,
                 P_Q_gens_loc_load_view),
                idq_θ_π_vhθh,
                I_sum_ynj_vj )
        
    else 

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh,
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )
    
        current_mismatch =
            get_nodes_current_mismatch(
                uh,
                (P_Q_gens_view,
                 P_Q_non_gens_view,
                 P_Q_gens_loc_load_view),
                I_sum_ynj_vj )

    end
        
    # -----------------------------------------------
    
    power_mismatch =
        uh .* conj.(current_mismatch)

    red_P_mismatch =
        real.(power_mismatch[
            setdiff(1:end, slack_bus_idx)])
    
    red_Q_mismatch =
        imag.(power_mismatch[
            setdiff(1:end,gens_idx)])
        
    red_ΔPQ .=
        vcat(red_P_mismatch,
             red_Q_mismatch)

   # ------------------------------------------------
    # update Inet_view  and Iinj_view
    # -----------------------------------------------

    Inet_view .= I_sum_ynj_vj

    Iinj_view .= get_Iinj(
        uh,
        P_Q_non_gens_view,
        P_Q_gens_loc_load_view,
        Inet_view )
    
    # ------------------------------------------------   
    # update Q 
    # ------------------------------------------------

    if with_δ_ed_eq == true
        
        Igen = idq_θ_π_vhθh 
        
        gens_S  = uh[gens_idx] .*
            conj.(Igen[gens_idx])
        
        P_Q_nodes_view[slack_bus_idx] .= [
            real(gens_S[slack_bus_idx]),
            imag(gens_S[slack_bus_idx])]

       for (idx, gen_S) in zip(gens_idx, gens_S)

           P_Q_gens_view[idx] .= [
               real( S_gens[idx]), imag(gen_S)]
       end
        
    else
  
        Igen = I_sum_ynj_vj +
            get_gens_local_load_current(
                uh, P_Q_gens_loc_load_view) 

        gens_S = uh[gens_idx] .* conj.(Igen[gens_idx])

        P_Q_nodes_view[slack_bus_idx] .=
            [real(gens_S[slack_bus_idx]),
             imag(gens_S[slack_bus_idx])]

        for (idx, gen_S) in zip(gens_idx, gens_S)

            P_Q_gens_view[idx] .= [
                real( S_gens[idx]), imag(gen_S) ]
        end
    end
    
    # ---------------------------------------------- 
    # update nodes_pf_U_view
    # ---------------------------------------------- 

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .= [
            real( uk ), imag( uk ) ]

    end

    return nothing
    
end



#-----------------------------------------------------
# new
#-----------------------------------------------------

function dynamic_power_balance_Jac_with_vh_θh!(
    Jac,
    red_vh_θh_0_view,
    working_vh_θh_view,
    pf_net_param )

    # ------------------------------------------------

    pf_net, pf_idx_and_state, _, _, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    # ------------------------------------------------
    
    Ybus, _, _, _, _ = pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (slack_bus_idx,
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
    
   # -------------------------------------------------- 

    working_vh_θh_view[ red_vh_θh_idx ] .=
        red_vh_θh_0_view

    uh = working_vh_θh_view[ vh_IDX ] .*
        exp.(im * working_vh_θh_view[ θh_IDX ])
    
    # -------------------------------------------------

    ν      = abs.( uh )
    
    Λ      = inv.(uh)

    Λ_conj = inv.( conj.(uh) )
        
    Θ      = angle.(uh)

    E      = Diagonal( inv.(ν) ) * uh     

    Ibus   = Ybus * uh
    
    nodes_size = length( Ibus )
    
    # ------------------------------------------------    
    # Jacobian: See Matpower Technical Note 2, page 8,
    # eq 29 and 32
    
    Gc_vh = Diagonal(uh) * ( Diagonal(conj.(Ibus)) +
        conj(Ybus) * Diagonal( conj.(uh) ) ) *
        Diagonal(inv.(ν))
    
    Gc_θh = im * Diagonal(uh) * (Diagonal(conj.(Ibus)) -
        conj(Ybus) * Diagonal( conj.(uh) ) )
    
    P_Gc_vh = real.(
        Gc_vh[setdiff(1:end, slack_bus_idx),
              setdiff(1:end, gens_idx)] )
    
    P_Gc_θh = real.(
        Gc_θh[setdiff(1:end, slack_bus_idx),
              setdiff(1:end, slack_bus_idx)] )
    
    Q_Gc_vh = imag.(
        Gc_vh[setdiff(1:end, gens_idx),
              setdiff(1:end, gens_idx) ] )
    
    Q_Gc_θh = imag.(
        Gc_θh[setdiff(1:end, gens_idx),
              setdiff(1:end, slack_bus_idx)] )
    
    Jac  .= vcat((hcat(P_Gc_vh, P_Gc_θh )),
                 (hcat(Q_Gc_vh, Q_Gc_θh)))

    return nothing
end


function power_balance_Jac_with_vh_θh!(
    Jac,
    red_vh_θh_0_view,
    working_vh_θh_view,
    pf_net_param )

    # --------------------------------------------------

    pf_net, pf_idx_and_state, _, _, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    # -------------------------------------------------
    
    Ybus, _, _, _, _ = pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state)=
         pf_idx_and_state

    (slack_bus_idx,
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
    
   # -------------------------------------------------- 

    working_vh_θh_view[ red_vh_θh_idx ] .=
        red_vh_θh_0_view

    uh = working_vh_θh_view[ vh_IDX ] .*
        exp.(im * working_vh_θh_view[ θh_IDX ])
    
    # -------------------------------------------------

    ν      = abs.( uh )
    
    Λ      = inv.(uh)

    Λ_conj = inv.( conj.(uh) )
        
    Θ      = angle.(uh)

    E      = Diagonal( inv.(ν) ) * uh     

    Ibus   = Ybus * uh
    
    nodes_size = length( Ibus )
    
    # ------------------------------------------------    
    # Jacobian: See Matpower Technical Note 2, page 8, eq 29 and 32
    
    Gc_vh = Diagonal(uh) * ( Diagonal(conj.(Ibus)) +
        conj(Ybus) * Diagonal( conj.(uh) ) ) *
        Diagonal(inv.(ν))
    
    Gc_θh = im * Diagonal(uh) * (
        Diagonal(conj.(Ibus)) - conj(Ybus) *
            Diagonal( conj.(uh) ) )
    
    P_Gc_vh = real.(
        Gc_vh[setdiff(1:end, slack_bus_idx),
              setdiff(1:end, gens_idx)] )
    
    P_Gc_θh = real.(
        Gc_θh[setdiff(1:end, slack_bus_idx),
              setdiff(1:end, slack_bus_idx)] )
    
    Q_Gc_vh = imag.(
        Gc_vh[setdiff(1:end, gens_idx),
              setdiff(1:end, gens_idx) ] )
    
    Q_Gc_θh = imag.(
        Gc_θh[setdiff(1:end, gens_idx),
              setdiff(1:end, slack_bus_idx)] )
    
    Jac  .= vcat((hcat(P_Gc_vh, P_Gc_θh )),
                 (hcat(Q_Gc_vh, Q_Gc_θh)))

    return nothing
end


#---------------------
# Dynamic powerflow function for nodes_edges_dynamics!
# --------------------

# -------
# Begin im model pf
# -------


function im_model_dynamic_power_balance_mismatch_with_vh_θh!(
    red_ΔPQ,
    red_vh_θh_0_view,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view ),
    pf_net_param,
    im_model_idq_pf_cal;
    with_δ_ed_eq = true )

    # ------------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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

    # ------------------------------------------------------

    idq_wt_pad_view, _ = im_model_idq_pf_cal
    
    # ------------------------------------------------------
    
    working_vh_θh_view[ red_vh_θh_idx ] .=
        red_vh_θh_0_view

    uh =  working_vh_θh_view[ vh_IDX ] .*
        exp.(im * working_vh_θh_view[ θh_IDX ])
        
    # ------------------------------------------------------ 

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    # S_gens_loc_load = x_from_xr_xi.(
    #     P_Q_gens_loc_load_view )
    
    # S_non_gens = x_from_xr_xi.( P_Q_non_gens_view )

    # ------------------------------------------------------ 

    if with_δ_ed_eq == true

        idq_wt_pad_view[gens_idx] .= [
            get_dynamic_idq_θ_π_vhθh(
                vh, θh, δ_ω_ed_eq...,
                ra_Xd_dash_Xq_dash... )
            for (vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash) in
                zip( abs.(uh[gens_idx]),
                     angle.(uh[gens_idx]),
                     δ_ω_ed_dash_eq_dash_view,
                     filter((x)->length(x) != 1,
                            ra_Xd_dash_Xq_dash_view) )  ]

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh,
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )

        # get_nodes_current_mismatch
        current_mismatch = get_nodes_current_mismatch_idq_θπ(
            uh,
            (P_Q_non_gens_view,
             P_Q_gens_loc_load_view),
            idq_wt_pad_view,
            I_sum_ynj_vj )
        
    else

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh,
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )
    
        current_mismatch = get_nodes_current_mismatch(
            uh,
            (P_Q_gens_view,
             P_Q_non_gens_view,
             P_Q_gens_loc_load_view),
            I_sum_ynj_vj )

    end
        
    # ------------------------------------------------------
    
    power_mismatch = uh .* conj.(current_mismatch)

    red_P_mismatch = real.(power_mismatch[
        setdiff(1:end, slack_bus_idx)])
    
    red_Q_mismatch = imag.(power_mismatch[
        setdiff(1:end,gens_idx)])
        
    red_ΔPQ .= vcat(red_P_mismatch, red_Q_mismatch)

   # ------------------------------------------------------
    # update Inet_view and Iinj_view
    # -----------------------------------------------------

    Inet_view .= I_sum_ynj_vj

    Iinj_view .= get_Iinj(
        uh, P_Q_non_gens_view,
        P_Q_gens_loc_load_view,
        Inet_view )
    
    # -----------------------------------------------------   
    # update Q 
    # -----------------------------------------------------


    if with_δ_ed_eq == true
        
        Igen = idq_wt_pad_view
        
        gens_S  = uh[gens_idx] .* conj.(
            Igen[gens_idx])
        
        P_Q_nodes_view[slack_bus_idx] .= [
            real(gens_S[slack_bus_idx]),
            imag(gens_S[slack_bus_idx])]

       for (idx, gen_S) in zip(gens_idx, gens_S)

           P_Q_gens_view[idx] .= [
               real( S_gens[idx]), imag(gen_S)]
       end
        
    else
  
        Igen = I_sum_ynj_vj +
            get_gens_local_load_current(
                uh, P_Q_gens_loc_load_view) 

        gens_S  = uh[gens_idx] .* conj.(Igen[gens_idx])

        P_Q_nodes_view[slack_bus_idx] .= [
            real(gens_S[slack_bus_idx]),
            imag(gens_S[slack_bus_idx])]

        for (idx, gen_S) in zip(gens_idx, gens_S)

            P_Q_gens_view[idx] .= [
                real( S_gens[idx]),
                imag(gen_S)]
        end
    end
    
    # ------------------------------------------------------ 
    # update nodes_pf_U_view
    # ------------------------------------------------------ 

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .= [
            real( uk ), imag( uk ) ]

    end

    return nothing
    
end


function im_model_power_balance_mismatch_with_vh_θh!(
    red_ΔPQ,
    red_vh_θh_0_view,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view ),
    pf_net_param,
    im_model_idq_pf_cal;
    with_δ_ed_eq = true )

    # ----------------------------------------------------
    # ----------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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
    
    # ----------------------------------------------------

    idq_wt_pad_view, _ =
        im_model_idq_pf_cal
    
    # ----------------------------------------------------

    working_vh_θh_view[ red_vh_θh_idx ] .=
        red_vh_θh_0_view

    uh =
        working_vh_θh_view[ vh_IDX ] .*
        exp.(im * working_vh_θh_view[ θh_IDX ])
        
    # --------------------------------------------------

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    # --------------------------------------------------

    if with_δ_ed_eq == true

        idq_wt_pad_view[gens_idx] .=
            [ get_dynamic_idq_θ_π_vhθh(
                vh, θh, δ_ω_ed_eq...,
                ra_Xd_dash_Xq_dash... )
              for (vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash) in
                  zip( abs.(uh[gens_idx]),
                       angle.(uh[gens_idx]),
                       δ_ω_ed_dash_eq_dash_view,
                       filter((x)->length(x) != 1,
                              ra_Xd_dash_Xq_dash_view) ) ]

        I_sum_ynj_vj =
            get_nodes_∑_ynj_x_vj(
                uh, Ynet,
                nodes_node_idx_and_incident_edges_other_node_idx )

        # get_nodes_current_mismatch
        current_mismatch =
             get_nodes_current_mismatch_idq_θπ(
                uh,
                (P_Q_non_gens_view,
                 P_Q_gens_loc_load_view),
                 idq_wt_pad_view ,
                 I_sum_ynj_vj )
        
    else

        I_sum_ynj_vj =
            get_nodes_∑_ynj_x_vj(
                uh, Ynet,
                nodes_node_idx_and_incident_edges_other_node_idx )
    
        current_mismatch =
            get_nodes_current_mismatch(
                uh,
                (P_Q_gens_view,
                 P_Q_non_gens_view,
                 P_Q_gens_loc_load_view),
                I_sum_ynj_vj )

    end
        
    # ------------------------------------------------------
    
    power_mismatch =
        uh .* conj.(current_mismatch)

    red_P_mismatch =
        real.(
            power_mismatch[
                setdiff(1:end, slack_bus_idx)])
    
    red_Q_mismatch =
        imag.(power_mismatch[setdiff(1:end,gens_idx)])
        
    red_ΔPQ .=
        vcat(red_P_mismatch, red_Q_mismatch)

   # ----------------------------------------------------
    # update Inet_view  and Iinj_view
    # ---------------------------------------------------

    Inet_view .=
        I_sum_ynj_vj

    Iinj_view .=
        get_Iinj(
            uh,
            P_Q_non_gens_view,
            P_Q_gens_loc_load_view,
            Inet_view )
    
    # --------------------------------------------------
    # update Q 
    # --------------------------------------------------


    if with_δ_ed_eq == true
        
        # Igen = idq_θ_π_vhθh

        Igen =
            idq_wt_pad_view 
        
        gens_S  =
            uh[gens_idx] .*
            conj.(Igen[gens_idx])
        
        P_Q_nodes_view[slack_bus_idx] .=
            [real(gens_S[slack_bus_idx]),
             imag(gens_S[slack_bus_idx])]

       for (idx, gen_S) in zip(gens_idx, gens_S)

           P_Q_gens_view[idx] .=
               [real( S_gens[idx]), imag(gen_S)]
       end
    else
  
        Igen =
            I_sum_ynj_vj +
            get_gens_local_load_current(
                uh, P_Q_gens_loc_load_view) 

        gens_S  =
            uh[gens_idx] .* conj.(Igen[gens_idx])

        P_Q_nodes_view[slack_bus_idx] .=
            [real(gens_S[slack_bus_idx]),
             imag(gens_S[slack_bus_idx])]

        for (idx, gen_S) in zip(gens_idx, gens_S)

            P_Q_gens_view[idx] .=
                [real( S_gens[idx]), imag(gen_S)]
        end
    end
    
    # ---------------------------------------------------
    # update nodes_pf_U_view
    # ---------------------------------------------------

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .=
            [ real( uk ), imag( uk ) ]

    end

    return nothing
    
end


# -------
# End im model pf
# -------


#-------------------------------------------------------

# -------
# Begin industrial model pf
# -------


function industrial_model_dynamic_power_balance_mismatch_with_vh_θh!( red_ΔPQ, red_vh_θh_0_view, ( working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view ), pf_net_param, industrial_model_idq_pf_cal; with_δ_ed_eq = true )

    # ------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) = pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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

    # ------------------------------------------------

    idq_wt_pad_view, _ = industrial_model_idq_pf_cal
    
    # ------------------------------------------------
    
    working_vh_θh_view[ red_vh_θh_idx ] .=
        red_vh_θh_0_view

    uh =  working_vh_θh_view[ vh_IDX ] .*
        exp.(im * working_vh_θh_view[ θh_IDX ])
        
    # ---------------------------------------------- 

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    # S_gens_loc_load = x_from_xr_xi.(
    #     P_Q_gens_loc_load_view )
    
    # S_non_gens = x_from_xr_xi.(
    #     P_Q_non_gens_view )

    # ----------------------------------------------- 

    if with_δ_ed_eq == true

        idq_wt_pad_view[gens_idx] .= [
            industrial_model_get_dynamic_idq_θ_π_vhθh(
                vh, θh, δ_ω_ed_eq...,
                ra_Xd_dash_Xq_dash... )
            for (vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash) in
                zip( abs.(uh[gens_idx]),
                     angle.(uh[gens_idx]),
                     δ_ω_ed_dash_eq_dash_view,
                     filter((x)->length(x) != 1,
                            ra_Xd_dash_Xq_dash_view)) ]

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh,
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )

        # get_nodes_current_mismatch
        # current_mismatch =
        #      get_nodes_current_mismatch(
        #          uh,
        #          (P_Q_non_gens_view,
        #           P_Q_gens_loc_load_view),
        #          idq_wt_pad_view,
        #          I_sum_ynj_vj )

        current_mismatch =
             get_current_mismatch_by_idq_with_loc_load(
                 uh,
                 (P_Q_non_gens_view,
                  P_Q_gens_loc_load_view,
                 idq_wt_pad_view,
                 I_sum_ynj_vj ))
        
    else

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh,
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )
    
        current_mismatch =
            get_nodes_current_mismatch(
                uh,
                (P_Q_gens_view,
                 P_Q_non_gens_view,
                 P_Q_gens_loc_load_view),
                I_sum_ynj_vj )

    end
        
    # -------------------------------------------------
    
    power_mismatch = uh .* conj.(current_mismatch)

    red_P_mismatch = real.(power_mismatch[
        setdiff(1:end, slack_bus_idx)])
    
    red_Q_mismatch = imag.(power_mismatch[
        setdiff(1:end,gens_idx)])
        
    red_ΔPQ .= vcat(red_P_mismatch, red_Q_mismatch)

    # -------------------------------------------------
    # update Inet_view and Iinj_view
    # -------------------------------------------------

    Inet_view .= I_sum_ynj_vj

    Iinj_view .= get_Iinj(
        uh,
        P_Q_non_gens_view,
        P_Q_gens_loc_load_view,
        Inet_view )
    
    # -----------------------------------------------   
    # update Q 
    # -----------------------------------------------


    if with_δ_ed_eq == true
        
        Igen = idq_wt_pad_view
        
        gens_S  = uh[gens_idx] .* conj.(
            Igen[gens_idx])
        
        P_Q_nodes_view[slack_bus_idx] .= [
            real(gens_S[slack_bus_idx]),
            imag(gens_S[slack_bus_idx])]

       for (idx, gen_S) in zip(gens_idx, gens_S)

           P_Q_gens_view[idx] .= [real(
               S_gens[idx]), imag(gen_S)]
       end
    else
  
        Igen = I_sum_ynj_vj +
            get_gens_local_load_current(
                uh, P_Q_gens_loc_load_view) 

        gens_S  = uh[gens_idx] .* conj.(
            Igen[gens_idx])

        P_Q_nodes_view[slack_bus_idx] .= [
            real(gens_S[slack_bus_idx]),
            imag(gens_S[slack_bus_idx])]

        for (idx, gen_S) in zip(gens_idx, gens_S)

            P_Q_gens_view[idx] .= [
                real( S_gens[idx]), imag(gen_S)]
        end
    end
    
    # ------------------------------------------------ 
    # update nodes_pf_U_view
    # ------------------------------------------------ 

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .= [
            real( uk ), imag( uk ) ]

    end

    return nothing
    
end


function industrial_model_power_balance_mismatch_with_vh_θh!(
    red_ΔPQ,
    red_vh_θh_0_view,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view),
    pf_net_param,
    industrial_model_idq_pf_cal;
    with_δ_ed_eq = true )

    # -----------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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
    
    # ------------------------------------------------------

    idq_wt_pad_view, _ = industrial_model_idq_pf_cal
    
    # ------------------------------------------------------

    working_vh_θh_view[ red_vh_θh_idx ] .=
        red_vh_θh_0_view

    uh =  working_vh_θh_view[ vh_IDX ] .*
        exp.(im * working_vh_θh_view[
            θh_IDX ])
        
    # ------------------------------------------------------ 

    S_gens = x_from_xr_xi.(
        P_Q_gens_view )

    # ------------------------------------------------------ 

    if with_δ_ed_eq == true

        idq_wt_pad_view[gens_idx] .= [
            industrial_model_get_dynamic_idq_θ_π_vhθh(
                vh, θh, δ_ω_ed_eq...,
                ra_Xd_dash_Xq_dash...)
            for (vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash) in
                zip( abs.(uh[gens_idx]),
                     angle.(uh[gens_idx]),
                     δ_ω_ed_dash_eq_dash_view,
                     filter((x)->length(x) != 1,
                            ra_Xd_dash_Xq_dash_view))]

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh,
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )

        # get_nodes_current_mismatch
        
        current_mismatch =  get_nodes_current_mismatch(
            uh,
            (P_Q_non_gens_view,
             P_Q_gens_loc_load_view),
            idq_wt_pad_view ,
            I_sum_ynj_vj )
        
    else

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh,
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx)
    
        current_mismatch = get_nodes_current_mismatch(
            uh,
            (P_Q_gens_view,
             P_Q_non_gens_view,
             P_Q_gens_loc_load_view),
            I_sum_ynj_vj )

    end
        
    # ------------------------------------------------------
    
    power_mismatch = uh .*
        conj.(current_mismatch)

    red_P_mismatch = real.(power_mismatch[
        setdiff(1:end, slack_bus_idx)])
    
    red_Q_mismatch = imag.(power_mismatch[
        setdiff(1:end,gens_idx)])
        
    red_ΔPQ .= vcat(red_P_mismatch,
                    red_Q_mismatch)

   # ------------------------------------------------------
    # update Inet_view  and Iinj_view
    # -----------------------------------------------------

    Inet_view .= I_sum_ynj_vj

    Iinj_view .= get_Iinj(
        uh, P_Q_non_gens_view,
        P_Q_gens_loc_load_view,
        Inet_view )
    
    # -----------------------------------------------------   
    # update Q 
    # -----------------------------------------------------

    if with_δ_ed_eq == true
        
        # Igen = idq_θ_π_vhθh

        Igen = idq_wt_pad_view 
        
        gens_S  = uh[gens_idx] .* conj.(
            Igen[gens_idx])
        
        P_Q_nodes_view[slack_bus_idx] .= [
            real(gens_S[slack_bus_idx]),
            imag(gens_S[slack_bus_idx])]

       for (idx, gen_S) in zip(gens_idx, gens_S)

           P_Q_gens_view[idx] .= [
               real( S_gens[idx]),
               imag(gen_S)]
       end
    else
  
        Igen = I_sum_ynj_vj +
            get_gens_local_load_current(
                uh, P_Q_gens_loc_load_view) 

        gens_S  = uh[gens_idx] .* conj.(
            Igen[gens_idx])

        P_Q_nodes_view[slack_bus_idx] .= [
            real(gens_S[slack_bus_idx]),
            imag(gens_S[slack_bus_idx])]

        for (idx, gen_S) in zip(gens_idx, gens_S)

            P_Q_gens_view[idx] .= [real(
                S_gens[idx]), imag(gen_S)]
        end
    end
    
    # ------------------------------------------------------ 
    # update nodes_pf_U_view
    # ------------------------------------------------------ 

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .= [
            real( uk ), imag( uk ) ]

    end

    return nothing
    
end

# -------
# End industrial  pf
# -------


########################################################
# ------------------------------------------------------
#  Important Powerflow Mismatch
# ------------------------------------------------------
########################################################

function get_ΔPQ_Δidq_mismatch(
    uh,
    Ynet,
    PQ_non_gens_pad,    
    δ_ω_ed_dash_eq_dash,
    ra_Xd_dash_Xq_dash_view,
    idq,
    nodes_node_idx_and_incident_edges_other_node_idx,
    slack_bus_idx,
    gens_nodes_idx )

    nodes_size = length(
        nodes_node_idx_and_incident_edges_other_node_idx )
                
   idq_wt_pad =
       [ idx ∈ gens_nodes_idx ? 
       idq[ n2s_gens_idx[idx] ] :
             [0.0, 0.0]
         for idx in 1:nodes_size ]
    
    I_sum_ynj_vj =
        [ sum( [ ynj * uh[node_idx]
                for (ynj, node_idx) in
                    zip( Y_bus_vec,
                         nth_node_idx_and_adj_nodes_idx )] )
          for (Y_bus_vec, nth_node_idx_and_adj_nodes_idx ) in
              zip(Ynet, nodes_node_idx_and_incident_edges_other_node_idx ) ]
    
    current_mismatch =
        I_sum_ynj_vj  +
        ( conj.( x_from_xr_xi.( PQ_non_gens_pad ) )) ./ ( conj.( uh )) - x_from_xr_xi.( idq_wt_pad )

    power_mismatch =
        uh .* conj.(current_mismatch)

    red_P_mismatch =
        real.(
            power_mismatch[setdiff(
                1:end, slack_bus_idx)])

    red_Q_mismatch =
        imag.(
            power_mismatch[setdiff(
                1:end,gens_idx)])
    
    gens_stator_equations =
        [ get_a_gen_stator_equations(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash..., id_iq )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash, id_iq ) in
              zip( abs.(uh[gens_idx]),
                   angle.(uh[gens_idx]),
                   δ_ω_ed_dash_eq_dash,
                   ra_Xd_dash_Xq_dash_view, idq ) ]

    flat_gens_stator_equations =
        [gens_stator_equations...;]

    ΔPQ_Δidq_mismatch =
        vcat(red_P_mismatch,
             red_Q_mismatch,
             flat_gens_stator_equations )    

end



function get_a_model_integrated_nll_dyn_intg_ΔPQ_Δidq_mismatch(
    intg_ΔPQ_id_iq, intg_vh_θh_id_iq, integ_param;
    intg_dyn_pf_fun_kwd_para =
        intg_dyn_pf_fun_kwd_para)


    #-------------------------------

    (;
     intg_vars_Idxs,
     dyn_pf_fun_kwd_para) =
         intg_dyn_pf_fun_kwd_para

    #-------------------------------
    
   (;intg_gens_vh_Idxs,
    intg_gens_θh_Idxs,
    intg_non_gens_nodes_vh_Idxs,
    intg_non_gens_nodes_θh_Idxs,
    intg_gen_id_Idxs,
    intg_gen_iq_Idxs) =
        intg_vars_Idxs

    #-------------------------------    
    
    (;loc_load_exist,
     dyn_pf_fun_kwd_nll_para_vars_Idxs,
     dyn_pf_fun_kwd_wll_para_vars_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         dyn_pf_fun_kwd_para
    
    #-------------------------------

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #-------------------------------

   # (;
   #  n2s_slack_gens_idx,
   #  n2s_non_slack_gens_idx,
   #  n2s_gens_idx,
   #  n2s_non_gens_idx,
   #  n2s_gens_with_loc_load_idxs,
   #  n2s_all_nodes_idx ) =
   #      dyn_pf_fun_kwd_net_idxs

    
    #-------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------    
    #-------------------------------

   (;P_gens_dyn_para_Idxs,
     Q_gens_dyn_para_Idxs,
     P_non_gens_dyn_para_Idxs,
     Q_non_gens_dyn_para_Idxs,
     δ_ed_eq_pf_dyn_para_Idxs,
     P_g_loc_load_dyn_para_Idxs,
     Q_g_loc_load_dyn_para_Idxs                   
    ) =
        dyn_pf_fun_kwd_wll_para_vars_Idxs 
    
    #-------------------------------
    #-------------------------------
    
    P_gens =
        integ_param[ P_gens_dyn_para_Idxs ]

    Q_gens =
        integ_param[ Q_gens_dyn_para_Idxs ]

    P_non_gens  =
        integ_param[ P_non_gens_dyn_para_Idxs ]

    Q_non_gens = 
        integ_param[ Q_non_gens_dyn_para_Idxs ]

    flat_δ_ω_ed_dash_eq_dash =
        integ_param[ δ_ed_eq_pf_dyn_para_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            integ_param[ P_g_loc_load_dyn_para_Idxs ]

        Q_g_loc_load =
            integ_param[ Q_g_loc_load_dyn_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------
    
    #-------------------------------
    #-------------------------------

    if length(flat_δ_ω_ed_dash_eq_dash) !=
        length(δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed )

        δ_ω_ed_dash_eq_dash =
            [ flat_δ_ω_ed_dash_eq_dash[idx]
              for idx in
                  δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed ]

    else

        δ_ω_ed_dash_eq_dash =
            flat_δ_ω_ed_dash_eq_dash

    end

    #-------------------------------

    intg_gens_vh =
        intg_vh_θh_id_iq[ intg_gens_vh_Idxs ]

    intg_gens_θh =
        intg_vh_θh_id_iq[ intg_gens_θh_Idxs ]

    intg_non_gens_nodes_vh =
        intg_vh_θh_id_iq[ intg_non_gens_nodes_vh_Idxs ]

    intg_non_gens_nodes_θh =
        intg_vh_θh_id_iq[ intg_non_gens_nodes_θh_Idxs ]

    gens_i_d =
        intg_vh_θh_id_iq[ intg_gen_id_Idxs ]

    gens_i_q =
        intg_vh_θh_id_iq[ intg_gen_iq_Idxs ]


    vh = [
        idx ∈ gens_nodes_idx ?
            intg_gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    intg_non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in 1:nodes_size ]

    θh = [
        idx ∈ gens_nodes_idx ?
            intg_gens_θh[
                n2s_gens_idx[ idx] ] :
                    intg_non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in 1:nodes_size ]

    # -------------------------------------

    gens_δ = first.( δ_ω_ed_dash_eq_dash )

    gens_ed_dash = third.( δ_ω_ed_dash_eq_dash )

    gens_eq_dash = fourth.( δ_ω_ed_dash_eq_dash )

    # -------------------------------------

    gens_ra = first.( gens_nodes_ra_Xd_dash_Xq_dash  )

    gens_Xd_dash = second.( gens_nodes_ra_Xd_dash_Xq_dash  )

    gens_Xq_dash = third.( gens_nodes_ra_Xd_dash_Xq_dash  )        

    # -----------------------------------
    # gens real part stator equation mismatch
    # -----------------------------------

    gens_stator_equations_mismatch_real =
        [ get_a_gen_real_stator_equations_mismatch(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash,
            a_id, a_iq )
          for ( a_vh, a_θh,
                a_δ, a_ed_dash, a_eq_dash,
                a_ra, a_X_d_dash, a_X_q_dash,
                a_id, a_iq ) in
              zip( intg_gens_vh, intg_gens_θh,
                   gens_δ, gens_ed_dash, gens_eq_dash,
                   gens_ra, gens_Xd_dash, gens_Xq_dash,
                   gens_i_d, gens_i_q ) ]

    # -----------------------------------
    # gens imag part stator equation mismatch
    # -----------------------------------

    gens_stator_equations_mismatch_imag =
        [ get_a_gen_imag_stator_equations_mismatch(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash,
            a_id, a_iq )
          for ( a_vh, a_θh,
                a_δ, a_ed_dash, a_eq_dash,
                a_ra, a_X_d_dash, a_X_q_dash,
                a_id, a_iq ) in
              zip( intg_gens_vh, intg_gens_θh,
                   gens_δ, gens_ed_dash, gens_eq_dash,
                   gens_ra, gens_Xd_dash, gens_Xq_dash,
                   gens_i_d, gens_i_q ) ]

    # -----------------------------------
    # gens active power mismatch
    # -----------------------------------

    gens_P_mismatch_no_loc_load =
        [get_nth_gen_P_pf_mismatch_no_loc_load(    
            vh,
            θh,    
            gens_δ[n2s_gens_idx[nth_idx]],
            gens_i_d[n2s_gens_idx[nth_idx]],
            gens_i_q[n2s_gens_idx[nth_idx]],  
            Ynet[ nth_idx ] ,
            nodes_node_idx_and_incident_edges_other_node_idx[
                nth_idx ], nth_idx )

         for nth_idx in
                gens_nodes_idx ]


    # -----------------------------------
    # gens reactive power mismatch
    # -----------------------------------

    gens_Q_mismatch_no_loc_load =
        [get_nth_gen_Q_pf_mismatch_no_loc_load(    
            vh,
            θh,    
            gens_δ[n2s_gens_idx[nth_idx]],
            gens_i_d[n2s_gens_idx[nth_idx]],
            gens_i_q[n2s_gens_idx[nth_idx]],  
            Ynet[ nth_idx ] ,
            nodes_node_idx_and_incident_edges_other_node_idx[
                nth_idx ], nth_idx )

         for nth_idx in
                gens_nodes_idx ]

    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------

    non_gens_P_mismatch =
        [ get_nth_non_gen_P_pf_mismatch( 
            vh,
            θh,           
            Ynet[ nth_idx ],
            nodes_node_idx_and_incident_edges_other_node_idx[
                nth_idx ], P_non_gens[ n2s_non_gens_idx[ nth_idx ] ],
            nth_idx )

          for nth_idx in
                non_gens_nodes_idx ]

    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------

    non_gens_Q_mismatch =
        [ get_nth_non_gen_Q_pf_mismatch( 
            vh,
            θh,           
            Ynet[ nth_idx ],
            nodes_node_idx_and_incident_edges_other_node_idx[
                nth_idx ],  Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ],
            nth_idx )

          for nth_idx in
                non_gens_nodes_idx ]

    # ------------------------------------

    intg_ΔPQ_id_iq .=
        vcat(gens_P_mismatch_no_loc_load,
             non_gens_P_mismatch,
             gens_Q_mismatch_no_loc_load,                 
             non_gens_Q_mismatch,
             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag)

    return nothing

end


function get_a_model_integrated_nll_dyn_full_ΔPQ_mismatch(
    full_ΔPQ, full_vh_θh, integ_param;
    full_dyn_pf_fun_kwd_para  =
        full_dyn_pf_fun_kwd_para  )

    #-------------------------------

    (;
     full_kwd_para,
     dyn_pf_fun_kwd_para) =
         full_dyn_pf_fun_kwd_para

    (; full_vars_Idxs,
     full_gens_id_iq ) =
         full_kwd_para

    #-------------------------------

   (;full_gens_vh_Idxs,
    full_gens_θh_Idxs,
    full_non_gens_nodes_vh_Idxs,
    full_non_gens_nodes_θh_Idxs) =
        full_vars_Idxs

    #-------------------------------
    
   (;gens_i_d_0,
    gens_i_q_0 ) =
        full_gens_id_iq 
    
    #-------------------------------
    
    (;loc_load_exist,
     dyn_pf_fun_kwd_nll_para_vars_Idxs,
     dyn_pf_fun_kwd_wll_para_vars_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         dyn_pf_fun_kwd_para
    
    #-------------------------------

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #-------------------------------

   # (;
   #  n2s_slack_gens_idx,
   #  n2s_non_slack_gens_idx,
   #  n2s_gens_idx,
   #  n2s_non_gens_idx,
   #  n2s_gens_with_loc_load_idxs,
   #  n2s_all_nodes_idx ) =
   #      dyn_pf_fun_kwd_net_idxs

    
    #-------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------
    
    #-------------------------------

   (;P_gens_dyn_para_Idxs,
     Q_gens_dyn_para_Idxs,
     P_non_gens_dyn_para_Idxs,
     Q_non_gens_dyn_para_Idxs,
     δ_ed_eq_pf_dyn_para_Idxs,
     P_g_loc_load_dyn_para_Idxs,
     Q_g_loc_load_dyn_para_Idxs                   
    ) =
        dyn_pf_fun_kwd_wll_para_vars_Idxs 
    
    #-------------------------------
    #-------------------------------
    
    P_gens =
        integ_param[ P_gens_dyn_para_Idxs ]

    Q_gens =
        integ_param[ Q_gens_dyn_para_Idxs ]

    P_non_gens  =
        integ_param[ P_non_gens_dyn_para_Idxs ]

    Q_non_gens = 
        integ_param[ Q_non_gens_dyn_para_Idxs ]

    flat_δ_ω_ed_dash_eq_dash =
        integ_param[ δ_ed_eq_pf_dyn_para_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            integ_param[ P_g_loc_load_dyn_para_Idxs ]

        Q_g_loc_load =
            integ_param[ Q_g_loc_load_dyn_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------
    
    #-------------------------------
    #-------------------------------

    if length(flat_δ_ω_ed_dash_eq_dash) !=
        length(δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed )

        δ_ω_ed_dash_eq_dash =
            [ flat_δ_ω_ed_dash_eq_dash[idx]
              for idx in
                  δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed ]

    else

        δ_ω_ed_dash_eq_dash =
            flat_δ_ω_ed_dash_eq_dash

    end

    #-------------------------------

    full_gens_vh =
        full_vh_θh[ full_gens_vh_Idxs ]

    full_gens_θh =
        full_vh_θh[ full_gens_θh_Idxs ]

    non_gens_vh = full_non_gens_nodes_vh =
        full_vh_θh[ full_non_gens_nodes_vh_Idxs ]

    non_gens_θh  = full_non_gens_nodes_θh =
        full_vh_θh[ full_non_gens_nodes_θh_Idxs ]

    vh = [
        idx ∈ gens_nodes_idx ?
            full_gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    full_non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in 1:nodes_size ]

    θh = [
        idx ∈ gens_nodes_idx ?
            full_gens_θh[
                n2s_gens_idx[ idx] ] :
                    full_non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in 1:nodes_size ]

    # -------------------------------------

    gens_δ = first.( δ_ω_ed_dash_eq_dash )

    gens_ed_dash = third.( δ_ω_ed_dash_eq_dash )

    gens_eq_dash = fourth.( δ_ω_ed_dash_eq_dash )

    # -------------------------------------

    gens_ra = first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash = second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash = third.( gens_nodes_ra_Xd_dash_Xq_dash )

    # -------------------------------------

    # # from closure

    gens_i_d = gens_i_d_0

    gens_i_q = gens_i_q_0

    # -----------------------------------
    # gens active power mismatch
    # -----------------------------------

    gens_P_mismatch_no_loc_load =
        [get_nth_gen_P_pf_mismatch_no_loc_load(    
            vh,
            θh,    
            gens_δ[ n2s_gens_idx[ nth_idx ] ],
            gens_i_d[ n2s_gens_idx[ nth_idx ] ],
            gens_i_q[ n2s_gens_idx[ nth_idx ] ],  
            Ynet[ nth_idx ] ,
            nodes_node_idx_and_incident_edges_other_node_idx[
                nth_idx ], nth_idx )

         for nth_idx in
                gens_nodes_idx ]


    # -----------------------------------
    # gens reactive power mismatch
    # -----------------------------------

    gens_Q_mismatch_no_loc_load =
        [get_nth_gen_Q_pf_mismatch_no_loc_load(    
            vh,
            θh,    
            gens_δ[ n2s_gens_idx[ nth_idx ] ],
            gens_i_d[ n2s_gens_idx[ nth_idx ] ],
            gens_i_q[ n2s_gens_idx[ nth_idx ] ],  
            Ynet[ nth_idx ] ,
            nodes_node_idx_and_incident_edges_other_node_idx[
                nth_idx ], nth_idx )

         for nth_idx in
                gens_nodes_idx ]

    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------

    non_gens_P_mismatch =
        [ get_nth_non_gen_P_pf_mismatch( 
            vh,
            θh,           
            Ynet[ nth_idx ],
            nodes_node_idx_and_incident_edges_other_node_idx[
                nth_idx ], P_non_gens[ n2s_non_gens_idx[ nth_idx ] ],
            nth_idx )

          for nth_idx in
                non_gens_nodes_idx ]

    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------

    non_gens_Q_mismatch =
        [ get_nth_non_gen_Q_pf_mismatch( 
            vh,
            θh,           
            Ynet[ nth_idx ],
            nodes_node_idx_and_incident_edges_other_node_idx[
                nth_idx ], Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ],
            nth_idx )

          for nth_idx in
                non_gens_nodes_idx ]

    # ------------------------------------

    full_ΔPQ .=
        vcat(gens_P_mismatch_no_loc_load,
             non_gens_P_mismatch,
             gens_Q_mismatch_no_loc_load,
             non_gens_Q_mismatch)

    return nothing

end


function get_a_model_integrated_pf_spcm_ΔPQ_mismatch_generic(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para )

    (; loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs ) =
         pf_kw_para
        
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh ) =
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


    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]


    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    
    # nodes_size = length( gens_nodes_idx ) +
    #     length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]


    # -----------------------------------
    # alternative
    # -----------------------------------

    # uh = vh .* exp.(im * θh)

    # -----------------------------------

    I_sum_Bnj_vj = get_nodes_∑_Bnj_x_vj(
        vh,
        θh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    # ----------------------------------------------   
    # update Q 
    # ------------------------------------------------

    Igen = I_sum_Bnj_vj[ transformed_gens_nodes_idx ] +
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,

            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_nodes_idx,
            gens_with_loc_load_idx )
    
    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .*
        conj.(
            Igen )
    
    # gens_S  = uh[gens_nodes_idx] .* conj.(Igen )

    P_gens_update = [ idx ∈ slack_gens_nodes_idx ?
        real( gens_S[ n2s_gens_idx[idx]]) : P_gens[
            n2s_gens_idx[idx]]
                      for idx in gens_nodes_idx]

    Q_gens_update = [ imag( gens_S[ n2s_gens_idx[idx]]) 
                      for idx in gens_nodes_idx]

    # ----------------------------------------------- 

    current_mismatch =
        get_nodes_current_mismatch(
            vh,
            θh,
            P_gens_update,
            Q_gens_update,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,    
            ( gens_nodes_idx,
              non_gens_nodes_idx,
              gens_with_loc_load_idx,
               all_nodes_idx, # transformed_all_nodes_idx #
              ),
            ( n2s_gens_idx,
             n2s_non_gens_idx,
              n2s_gens_with_loc_load_idxs,
              n2s_all_nodes_idx #
              ),
            I_sum_Bnj_vj;
            loc_load_exist =
                loc_load_exist)
    
    # -----------------------------------------------

    power_mismatch =
        vh .* exp.(im * θh) .* conj.(current_mismatch)

    red_P_mismatch =
        real.(power_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)])

    red_Q_mismatch =
        imag.(power_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)])

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)

    return nothing


end

#---------------------------------------------------
#---------------------------------------------------
# Mismatch ΔI
#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
# mismatch by Ynet
#---------------------------------------------------

 
function get_ΔI_mismatch_by_Ynet(
    red_ΔI_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para,
    Ynet_wt_nodes_idx_wt_adjacent_nodes =
        Ynet_wt_nodes_idx_wt_adjacent_nodes)

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     # Ynet,
     # nodes_idx_with_adjacent_nodes_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx ) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
                     :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              :Ynet,
              :nodes_idx_with_adjacent_nodes_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx ) )
    
    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]

    # -----------------------------------


    (;Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx))

    #-------------------------------
    
    I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
        vh,
        θh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    # ------------------------------------------------   

    Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
        get_gens_loc_load_current(
            vh,
            θh,           
            P_g_loc_load,
            Q_g_loc_load,

            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_nodes_idx,
            gens_with_loc_load_idx )

    

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]
             ) .* conj.(
            Igen )
    
    # gens_S  = uh[gens_nodes_idx] .* conj.(Igen )

    P_gens_update = [ idx ∈ slack_gens_nodes_idx ?
        real( gens_S[ n2s_gens_idx[idx] ] ) : P_gens[
            n2s_gens_idx[idx]]
                      for idx in gens_nodes_idx]

    Q_gens_update = [ imag( gens_S[ n2s_gens_idx[idx]]) 
                      for idx in gens_nodes_idx]

    # ----------------------------------------------- 

    current_mismatch =
        get_nodes_current_mismatch(
            vh,
            θh,
            P_gens_update,
            Q_gens_update,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,    
            ( gens_nodes_idx,
              non_gens_nodes_idx,
              gens_with_loc_load_idx,
              # transformed_all_nodes_idx #
               all_nodes_idx, 
              ),
            ( n2s_gens_idx,
             n2s_non_gens_idx,
              n2s_gens_with_loc_load_idxs,
              n2s_all_nodes_idx #
              ),
            I_sum_ynj_vj;
            loc_load_exist =
                loc_load_exist)
    
    # -----------------------------------------------

    red_real_mismatch =
        real.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)])

    red_imag_mismatch =
        imag.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)])

    red_ΔI_x .=
        vcat(red_real_mismatch,
             red_imag_mismatch)

    return nothing


end


#---------------------------------------------------
# mismatch by Yπ_net
#---------------------------------------------------

function get_ΔI_mismatch_by_Yπ_net(
    red_ΔI_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para,
    Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes  =
        Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes )

    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     # Ynet,
     # nodes_idx_with_adjacent_nodes_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx ) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
                     :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              # :Ynet,
              # :nodes_idx_with_adjacent_nodes_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx) )

    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]

    # -----------------------------------
    
   (;Yπ_net,
    Yshunt,
    nodes_idx_with_adjacent_nodes_idx ) =
        NamedTupleTools.select(
            Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes ,
            (:Yπ_net,
             :Yshunt,
             :nodes_idx_with_adjacent_nodes_idx ))

    # -----------------------------------
    
    I_sum_ynj_vj =
        get_nodes_∑_ynj_x_vj_by_Yπ_net(
            vh,
            θh,
            Yπ_net,
            Yshunt,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx )
    
    # ------------------------------------------------   
    # update Q 
    # ------------------------------------------------

    Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
        get_gens_loc_load_current(
            vh,
            θh,           
            P_g_loc_load,
            Q_g_loc_load,

            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_nodes_idx,
            gens_with_loc_load_idx )

    
    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .* conj.(
            Igen )
    
    # gens_S  = uh[gens_nodes_idx] .* conj.(Igen )

    P_gens_update = [ idx ∈ slack_gens_nodes_idx ?
        real( gens_S[ n2s_gens_idx[idx] ] ) : P_gens[
            n2s_gens_idx[idx]]
                      for idx in gens_nodes_idx]

    Q_gens_update = [ imag( gens_S[ n2s_gens_idx[idx]]) 
                      for idx in gens_nodes_idx]

    # ----------------------------------------------- 

    current_mismatch =
        get_nodes_current_mismatch(
            vh,
            θh,
            P_gens_update,
            Q_gens_update,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,    
            ( gens_nodes_idx,
              non_gens_nodes_idx,
              gens_with_loc_load_idx,
               all_nodes_idx, # transformed_all_nodes_idx #
              ),
            ( n2s_gens_idx,
             n2s_non_gens_idx,
              n2s_gens_with_loc_load_idxs,
              n2s_all_nodes_idx #
              ),
            I_sum_ynj_vj;
            loc_load_exist =
                loc_load_exist)
    
    # -----------------------------------------------

    # power_mismatch =
    #     vh .* exp.(im * θh) .* conj.(current_mismatch)

    red_real_mismatch =
        real.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)])

    red_Imag_mismatch =
        imag.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)])

    red_ΔI_x .=
        vcat(red_real_mismatch,
             red_Imag_mismatch)

    return nothing


end


#---------------------------------------------------
# mismatch by Ybus
#---------------------------------------------------

function get_ΔI_mismatch_by_Ybus(
    red_ΔI_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para,
    Ybus =
        Ybus,
    use_autodiff =
        true)

    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     # Ynet,
     # nodes_idx_with_adjacent_nodes_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
                     :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              # :Ynet,
              # :nodes_idx_with_adjacent_nodes_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx) )

    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]


    # -----------------------------------

    I_sum_ynj_vj =
        get_nodes_∑_ynj_x_vj_by_Ybus(
            vh,
            θh;
            Ybus)
    
    # ------------------------------------------------   
    # update Q 
    # ------------------------------------------------

    Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
        get_gens_loc_load_current(
            vh,
            θh,           
            P_g_loc_load,
            Q_g_loc_load,

            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_nodes_idx,
            gens_with_loc_load_idx )

    
    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .* conj.(
            Igen )
    
    # gens_S  = uh[gens_nodes_idx] .* conj.(Igen )

    P_gens_update = [ idx ∈ slack_gens_nodes_idx ?
        real( gens_S[ n2s_gens_idx[idx] ] ) : P_gens[
            n2s_gens_idx[idx]]
                      for idx in gens_nodes_idx]

    Q_gens_update = [ imag( gens_S[ n2s_gens_idx[idx]]) 
                      for idx in gens_nodes_idx]

    # ----------------------------------------------- 

    current_mismatch =
        get_nodes_current_mismatch(
            vh,
            θh,
            P_gens_update,
            Q_gens_update,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,    
            ( gens_nodes_idx,
              non_gens_nodes_idx,
              gens_with_loc_load_idx,
               all_nodes_idx, # transformed_all_nodes_idx #
              ),
            ( n2s_gens_idx,
             n2s_non_gens_idx,
              n2s_gens_with_loc_load_idxs,
              n2s_all_nodes_idx #
              ),
            I_sum_ynj_vj;
            loc_load_exist =
                loc_load_exist)
        
    # -----------------------------------------------

    red_real_mismatch =
        real.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)])

    red_imag_mismatch =
        imag.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)])

    red_ΔI_x .=
        vcat(red_real_mismatch,
             red_imag_mismatch)

    return nothing


end

#---------------------------------------------------
#---------------------------------------------------
# Mismatch ΔPQ
#---------------------------------------------------
#---------------------------------------------------

#---------------------------------------------------
# mismatch ΔPQ by Ynet
#---------------------------------------------------

function get_ΔPQ_mismatch_by_Ynet(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param
    ;pf_kw_para =
        pf_kw_para,
    Ynet_wt_nodes_idx_wt_adjacent_nodes =
        Ynet_wt_nodes_idx_wt_adjacent_nodes)

    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     # Ynet,
     # nodes_idx_with_adjacent_nodes_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx ) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,

              :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              # :Ynet,
              # :nodes_idx_with_adjacent_nodes_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx) )

    #-------------------------------
    #-------------------------------

    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]

    #-------------------------------
        
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    # Note components_injection is assigned negative
    # value or positive value based on wether it is
    # a gen or non-gen
    
    P_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (P_non_gens[ n2s_non_gens_idx[
                nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]])  : 
            -(P_gens[ n2s_gens_idx[nth_idx]])        
        for nth_idx in all_nodes_idx ]


    Q_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (Q_non_gens[ n2s_non_gens_idx[ nth_idx]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(Q_gens[ n2s_gens_idx[ nth_idx ]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]]) :
            -(Q_gens[ n2s_gens_idx[nth_idx]])    
        for nth_idx in all_nodes_idx ]
    
    #-------------------------------    

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              Ynet_wt_nodes_idx_wt_adjacent_nodes,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))

    #-------------------------------    

    P_line_injection =
        get_nodes_real_uh_x_∑_ynj_x_vj(
            vh, θh;
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            all_nodes_idx,
            n2s_all_nodes_idx)

    Q_line_injection =
        get_nodes_imag_uh_x_∑_ynj_x_vj(
            vh, θh;
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            all_nodes_idx,
            n2s_all_nodes_idx)

    #-------------------------------    
    
    P_mismatch = P_line_injection .+
        P_nodal_components_injection
    
    Q_mismatch = Q_line_injection .+
        Q_nodal_components_injection
    
    # ----------------------------------------------- 

    red_P_mismatch =
        P_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)]

    red_Q_mismatch =
        Q_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)]

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)
    
    return nothing


end

#---------------------------------------------------
# mismatch ΔPQ by Yπ_net
#---------------------------------------------------

function get_ΔPQ_mismatch_by_Yπ_net(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param
    ;pf_kw_para =
        pf_kw_para,
    Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
        Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes )

    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     # Ynet,
     # nodes_idx_with_adjacent_nodes_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx ) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
                     :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              # :Ynet,
              # :nodes_idx_with_adjacent_nodes_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx) )

    #-------------------------------
    #-------------------------------

    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]

    #-------------------------------
        
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    # Note components_injection is assigned negative
    # value or positive value based on wether it is
    # a gen or non-gen
    
    P_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (P_non_gens[ n2s_non_gens_idx[
                nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]])  : 
            -(P_gens[ n2s_gens_idx[nth_idx]])        
        for nth_idx in all_nodes_idx ]


    Q_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (Q_non_gens[ n2s_non_gens_idx[ nth_idx]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(Q_gens[ n2s_gens_idx[ nth_idx ]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]]) :
            -(Q_gens[ n2s_gens_idx[nth_idx]])    
        for nth_idx in all_nodes_idx ]
    
    #-------------------------------


    (;Yπ_net,
     Yshunt,
     nodes_idx_with_adjacent_nodes_idx ) =
         NamedTupleTools.select(
             Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
             (:Yπ_net,
              :Yshunt,
              :nodes_idx_with_adjacent_nodes_idx ))

    #-------------------------------

     ∑_S_injection_by_Yπ_net =
        get_nodes_∑_Sh_injection_by_Yπ_net(
            vh, θh; Yπ_net, Yshunt,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx )

    P_line_injection =
        real.(∑_S_injection_by_Yπ_net)

    Q_line_injection =
        imag.(∑_S_injection_by_Yπ_net)

    P_mismatch = P_line_injection .+
        P_nodal_components_injection
    
    Q_mismatch = Q_line_injection .+
        Q_nodal_components_injection
    
    # ----------------------------------------------- 

    red_P_mismatch =
        P_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)]

    red_Q_mismatch =
        Q_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)]

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)
    

    return nothing


end


function get_ΔPQ_mismatch_by_S_inj_Yπ_net(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param;
    pf_kw_para =
        pf_kw_para,
    Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
        Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes )

    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx ) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
                     :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx) )
    
    #-------------------------------
    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
           for idx in all_nodes_idx ]
    
    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    # -----------------------------------

    # Note components_injection is assigned negative
    # value or positive value based on wether it is
    # a gen or non-gen
    
    P_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (P_non_gens[ n2s_non_gens_idx[ nth_idx ] ]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]])  : 
            -(P_gens[ n2s_gens_idx[nth_idx]])        
        for nth_idx in all_nodes_idx ]


    Q_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(Q_gens[ n2s_gens_idx[ nth_idx ]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]]) :
            -(Q_gens[ n2s_gens_idx[nth_idx]])    
        for nth_idx in all_nodes_idx ]

    # -----------------------------------

    (;Yπ_net,
     Yshunt,
     nodes_idx_with_adjacent_nodes_idx ) =
         NamedTupleTools.select(
             Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
             (:Yπ_net,
              :Yshunt,
              :nodes_idx_with_adjacent_nodes_idx ))

    # -----------------------------------

     ∑_S_injection_by_Yπ_net =
        get_nodes_∑_Sh_injection_by_Yπ_net(
            vh, θh; Yπ_net, Yshunt,
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx )

    P_line_injection =
        real.(∑_S_injection_by_Yπ_net)

    Q_line_injection =
        imag.(∑_S_injection_by_Yπ_net)

    P_mismatch = P_line_injection .+
        P_nodal_components_injection
    
    Q_mismatch = Q_line_injection .+
        Q_nodal_components_injection
    
    # ----------------------------------------------- 

    red_P_mismatch =
        P_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)]

    red_Q_mismatch =
        Q_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)]

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)
    

    return nothing

end

#---------------------------------------------------
# mismatch ΔPQ by Ybus
#---------------------------------------------------


function get_ΔPQ_mismatch_by_Ybus(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param
    ;pf_kw_para =
        pf_kw_para,
    Ybus =
        Ybus,
    use_autodiff = true )
    
    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     # Ynet,
     # nodes_idx_with_adjacent_nodes_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
                     :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              # :Ynet,
              # :nodes_idx_with_adjacent_nodes_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx) )

    #-------------------------------
    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]


    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
        
    # -----------------------------------

    # Note components_injection is assigned negative
    # value or positive value based on wether it is
    # a gen or non-gen
    
    P_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (P_non_gens[ n2s_non_gens_idx[ nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]])  : 
            -(P_gens[ n2s_gens_idx[nth_idx]])        
        for nth_idx in all_nodes_idx ]


    Q_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(Q_gens[ n2s_gens_idx[ nth_idx ]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]]) :
            -(Q_gens[ n2s_gens_idx[nth_idx]])    
        for nth_idx in all_nodes_idx ]
    
    # -----------------------------------

    U_bus = vh .* exp.(im * θh)

    if use_autodiff == true

        I_bus = Matrix(Ybus) * U_bus
        
    else

        I_bus = Ybus * U_bus

        
    end

    S_bus = U_bus .* conj.(I_bus)

    P_line_injection = real.( S_bus )

    Q_line_injection = imag.( S_bus )
    
    # -----------------------------------
    
    P_mismatch = P_line_injection .+
        P_nodal_components_injection
    
    Q_mismatch = Q_line_injection .+
        Q_nodal_components_injection
    
    # ----------------------------------------------- 

    red_P_mismatch =
        P_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)]

    red_Q_mismatch =
        Q_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)]

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)

    return nothing


end



function get_ΔPQ_mismatch_by_sparse_Ybus(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param;
    Ybus =
        Ybus,
    pf_kw_para =
        pf_kw_para,
    use_autodiff = true )
    
    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
      gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx,
     
     transformed_red_P_mismatch_idx,
     transformed_red_Q_mismatch_idx) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx,

              :transformed_red_P_mismatch_idx,
              :transformed_red_Q_mismatch_idx) )

    #-------------------------------
    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]


    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
        
    # -----------------------------------

    # Note components_injection is assigned negative
    # value or positive value based on wether it is
    # a gen or non-gen
    
    P_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (P_non_gens[ n2s_non_gens_idx[ nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]])  : 
            -(P_gens[ n2s_gens_idx[nth_idx]])        
        for nth_idx in all_nodes_idx ]

    Q_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(Q_gens[ n2s_gens_idx[ nth_idx ]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]]) :
            -(Q_gens[ n2s_gens_idx[nth_idx]])    
        for nth_idx in all_nodes_idx ]
    
    # -----------------------------------

    if use_autodiff == true

        S_bus = (vh .* exp.(im * θh)) .* conj.(
            Matrix(Ybus) * ( vh .* exp.(im * θh) ) )
        
    else

        S_bus = (vh .* exp.(im * θh)) .* conj.(
            Ybus * ( vh .* exp.(im * θh) ) )
        
    end
    
    P_line_injection = real.( S_bus )

    Q_line_injection = imag.( S_bus )
        
    # -----------------------------------
    
    P_mismatch = P_line_injection .+
        P_nodal_components_injection
    
    Q_mismatch = Q_line_injection .+
        Q_nodal_components_injection

    # -----------------------------------
    
    red_P_mismatch =
        P_mismatch[ transformed_red_P_mismatch_idx ]

    red_Q_mismatch =
        Q_mismatch[ transformed_red_Q_mismatch_idx ]

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)

    return nothing


end

#---------------------------------------------------
# experimental  mismatch ΔPQ 
#---------------------------------------------------


#---------------------------------------------------

function get_ΔPQ_mismatch_by_experimental_ph_qh_Yπ_net(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para,
    Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
        Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
    use_Yπ_net_bool = true )

    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     # Ynet,
     # nodes_idx_with_adjacent_nodes_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx ) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
                     :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              # :Ynet,
              # :nodes_idx_with_adjacent_nodes_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx) )
        
    #-------------------------------
    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]

    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    # -----------------------------------

    # Note components_injection is assigned negative
    # value or positive value based on wether it is
    # a gen or non-gen
    
    P_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (P_non_gens[ n2s_non_gens_idx[
                nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]])  : 
            -(P_gens[ n2s_gens_idx[nth_idx]])        
        for nth_idx in all_nodes_idx ]


    Q_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (Q_non_gens[ n2s_non_gens_idx[ nth_idx]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(Q_gens[ n2s_gens_idx[ nth_idx ]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]]) :
            -(Q_gens[ n2s_gens_idx[nth_idx]])    
        for nth_idx in all_nodes_idx ]
    
    # -----------------------------------

    if use_Yπ_net_bool == true

        (;Yπ_net,
         Yshunt,
         nodes_idx_with_adjacent_nodes_idx ) =
             NamedTupleTools.select(
                 Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
                 (:Yπ_net,
                  :Yshunt,
                  :nodes_idx_with_adjacent_nodes_idx ))

         ∑_S_injection_by_Yπ_net =
            get_nodes_∑_Sh_injection_by_Yπ_net(
                vh, θh; Yπ_net, Yshunt,
                nodes_idx_with_adjacent_nodes_idx,
                n2s_all_nodes_idx )

        P_line_injection =
            real.(∑_S_injection_by_Yπ_net)

        Q_line_injection =
            imag.(∑_S_injection_by_Yπ_net)

    else

         (;Ynet,
          nodes_idx_with_adjacent_nodes_idx) =
              NamedTupleTools.select(
                  pf_kw_para,
                  (:Ynet,
                   :nodes_idx_with_adjacent_nodes_idx))

        P_line_injection =
            get_nodes_real_uh_x_∑_ynj_x_vj(
                vh, θh; Ynet,
                nodes_idx_with_adjacent_nodes_idx,
                all_nodes_idx,
                n2s_all_nodes_idx)

        Q_line_injection =
            get_nodes_imag_uh_x_∑_ynj_x_vj(
                vh, θh; Ynet,
                nodes_idx_with_adjacent_nodes_idx,
                all_nodes_idx,
                n2s_all_nodes_idx)

    end

    P_mismatch = P_line_injection .+
        P_nodal_components_injection
    
    Q_mismatch = Q_line_injection .+
        Q_nodal_components_injection
    
    # ----------------------------------------------- 

    red_P_mismatch =
        P_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)]

    red_Q_mismatch =
        Q_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)]

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)
    

    return nothing

end


function get_ΔPQ_mismatch_by_sparse_experimental_Ybus(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param;
    pf_kw_para =
        pf_kw_para,
    Ybus = Ybus,
    sp_ybus_nzv_real =
        sp_ybus_nzv_real,
    sp_ybus_nzv_imag =
        sp_ybus_nzv_imag,
    sp_ybus_I=
        sp_ybus_I,
    sp_ybus_J=
        sp_ybus_J,
    red_vh_θh_DiffCache =
        red_vh_θh_DiffCache ,
    sp_ybus_nzv_real_cache =
        sp_ybus_nzv_real_cache,
    sp_ybus_nzv_imag_cache =
        sp_ybus_nzv_imag_cache,
    sparse_nzvalues_DiffCache =
        sparse_nzvalues_DiffCache,
    I_bus_cache = I_bus_cache,
    I_bus = I_bus)
    
    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh,

     # Ynet,
     # nodes_idx_with_adjacent_nodes_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     P_gens_sta_para_Idxs,
     Q_gens_sta_para_Idxs,
     P_non_gens_sta_para_Idxs,
     Q_non_gens_sta_para_Idxs,
     P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs,

     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx,

     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     transformed_slack_gens_nodes_idx,
     # transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     # transformed_non_gens_nodes_idx,
     # transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx) = 
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
                     :slack_gens_vh,
              :slack_gens_θh,
              :gens_vh,
              :non_slack_gens_vh,

              # :Ynet,
              # :nodes_idx_with_adjacent_nodes_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :P_gens_sta_para_Idxs,
              :Q_gens_sta_para_Idxs,
              :P_non_gens_sta_para_Idxs,
              :Q_non_gens_sta_para_Idxs,
              :P_g_loc_load_sta_para_Idxs,
              :Q_g_loc_load_sta_para_Idxs,

              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx,

              :n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx,

              :transformed_slack_gens_nodes_idx,
              # :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              # :transformed_non_gens_nodes_idx,
              # :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx) )

    #-------------------------------
    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    red_vh_θh_DiffCache =
        get_tmp(red_vh_θh_DiffCache,
                red_vh_θh_x)

    red_vh_θh_DiffCache .= red_vh_θh_x
    
    non_slack_gens_θh = red_vh_θh_DiffCache[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_DiffCache[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_DiffCache[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]


    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
        
    # -----------------------------------

    # Note components_injection is assigned negative
    # value or positive value based on wether it is
    # a gen or non-gen
    
    P_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (P_non_gens[ n2s_non_gens_idx[ nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]])  : 
            -(P_gens[ n2s_gens_idx[nth_idx]])        
        for nth_idx in all_nodes_idx ]


    Q_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(Q_gens[ n2s_gens_idx[ nth_idx ]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]]) :
            -(Q_gens[ n2s_gens_idx[nth_idx]])    
        for nth_idx in all_nodes_idx ]
    
    # -----------------------------------
    
    # U_bus = vh .* exp.(im * θh)

    # I_bus = Matrix(Ybus) * U_bus
    
    # sparse_nzvalues_DiffCache =
    #         get_tmp(sparse_nzvalues_DiffCache,
    #                 sp_ybus_nzv )

    # sparse_nzvalues_DiffCache .= sp_ybus_nzv

    # sp_ybus_nzv_real_cache =
    #     get_tmp(sp_ybus_nzv_real_cache,
    #             sp_ybus_nzv_real )

    # sp_ybus_nzv_real_cache .=
    #     sp_ybus_nzv_real

    # sp_ybus_nzv_imag_cache =
    #     get_tmp(sp_ybus_nzv_imag_cache,
    #             sp_ybus_nzv_imag )
    
    # sp_ybus_nzv_imag_cache .=
    #     sp_ybus_nzv_imag

    # sp_ybus_nzv =
    #     sp_ybus_nzv_real_cache .+
    #     im * sp_ybus_nzv_imag_cache

    # # sp_ybus_nzv =
    # #     sp_ybus_nzv_real .+ im * sp_ybus_nzv_imag
    
    # Ybus = SparseArrays.sparse(
    #     sp_ybus_I,
    #     sp_ybus_J,
    #     sp_ybus_nzv)

    # I_bus_cache = get_tmp(I_bus_cache,I_bus)
    
    # mul!(I_bus_cache, Ybus, U_bus)
    
    # S_bus = U_bus .* conj.(I_bus_cache)
    
    # U_bus  =
    #     vh .* exp.(im * θh)

    # I_bus_g =
    #     Ybus_G * vh
    
    # I_bus_b =
    #     Ybus_B * vh

    # I_bus_t =
    #     (I_bus_g .+ im * I_bus_b
    #      ) .* exp.(im * θh)
    
    # I_bus_t = Ybus * U_bus

    # S_bus = U_bus .* conj.(I_bus_t)

    # S_bus = vh .* exp.(im * θh) .* (I_bus_g .-
    #     im * I_bus_b) .* exp.(im * θh)

    # S_bus = vh .* (I_bus_g .- im * I_bus_b)

    # P_line_injection = real.( S_bus )

    # Q_line_injection = imag.( S_bus )

    Ybus_G = real.(Ybus)
    
    Ybus_B = imag.(Ybus)

    # sq_vh = spdiagm((vh).^2)

    # P_line_injection =
    #     Array(diag(sq_vh * Ybus_G))

    # Q_line_injection =
    #     Array(diag((-sq_vh * Ybus_B)))

    # P_line_injection =
    #     Array(diag(spdiagm(vh) * (Ybus_G * vh) ))

    # Q_line_injection =
    #     Array(diag(-spdiagm(vh) * (Ybus_B * vh )))

    
    # Ybus_G = real.(Ybus)
    
    # Ybus_B = imag.(Ybus)


    # P_line_injection =
    #     Array(vh .* (Ybus_G * (vh ) )
    #     # Array(diag(spdiagm(vh) * (Ybus_G * vh) ))

    # Q_line_injection =
    #     Array(-vh  .* (Ybus_B * vh ) )
    #     # Array(diag(-spdiagm(vh) * (Ybus_B * vh )))
    

    # α = cos.(θh)
    
    # β = sin.(θh)
    
    # ( (vh .* α) .+ im * (vh .*  β) ) .* (Ybus_G - im * Ybus_B) * (
    #      ( (vh .* α) .- im * (vh .*  β) ) ) 

    # Ybus_G * (vh .* α) - im * Ybus_B * (vh .* α)

    # .- (Ybus_G * im * (vh .*  β) -  Ybus_B * (vh .*  β) )
    

    # A =       (Ybus_G * (vh .* α) .+ Ybus_B * (vh .*  β))
    
    # B = -im * (Ybus_B * (vh .* α) .+ Ybus_G * (vh .*  β))
    
    # C =       (Ybus_B * (vh .* α) .+ Ybus_G * (vh .*  β))

    # ( (vh .* α) .+ im * (vh .*  β) ) .* (A .- im * C)

    # ((vh .* α) .* A) .+ (im * (vh .* β) .* A ) .- im * ((vh .* α) .* C)
    # ((vh .*  β) .* C)

    # ((vh .* α) .* A) .+ ((vh .*  β) .* C) .- im * (
    #     ((vh .* α) .* C) - ((vh .* β) .* A ) )
    
    # Ybus_G = real.(Ybus)
    
    # Ybus_B = imag.(Ybus)

    # α = cos.(θh)
    
    # β = sin.(θh)

    # A = (Ybus_G * (vh .* α) .+ Ybus_B * (vh .*  β))
        
    # C = (Ybus_B * (vh .* α) .+ Ybus_G * (vh .*  β))

    # P_line_injection = ((vh .* α) .* A) .+ ((vh .*  β) .* C)

    # Q_line_injection = ((vh .* β) .* A ) .- ((vh .* α) .* C)   

    
    # -----------------------------------
    
    P_mismatch = P_line_injection .+
        P_nodal_components_injection
    
    Q_mismatch = Q_line_injection .+
        Q_nodal_components_injection

    # -----------------------------------
    
    red_P_mismatch =
        P_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)]

    red_Q_mismatch =
        Q_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)]

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)

    return nothing

end

#---------------------------------------------------
#---------------------------------------------------

function get_a_model_integrated_pf_sta_ΔI_mismatch_generic(
    red_ΔI_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para )

    (; loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs ) =
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :pf_kw_gens_vh_slack_θh_para,
              :pf_kw_net_para,
              :pf_kw_var_idxs,
              :pf_kw_PQ_para_idxs,
              :pf_kw_nodes_types_idxs,
              :pf_kw_n2s_idxs))

    
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             pf_kw_gens_vh_slack_θh_para,
             (:slack_gens_vh,
              :slack_gens_θh,

              :gens_vh,
              :non_slack_gens_vh ))
    
     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              pf_kw_net_para,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))

     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          NamedTupleTools.select(
              pf_kw_var_idxs,
              (:red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx ))

    
     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs ) =
          NamedTupleTools.select(
              pf_kw_PQ_para_idxs,
              (:P_gens_sta_para_Idxs,
               :Q_gens_sta_para_Idxs,
               :P_non_gens_sta_para_Idxs,
               :Q_non_gens_sta_para_Idxs,
               :P_g_loc_load_sta_para_Idxs,
               :Q_g_loc_load_sta_para_Idxs ))

     (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx ) =
         NamedTupleTools.select(
             pf_kw_nodes_types_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx ))    

     (; n2s_slack_gens_idx,
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


    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]


    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    
    # nodes_size = length( gens_nodes_idx ) +
    #     length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]


    # -----------------------------------

    I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
        vh,
        θh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    # ------------------------------------------------   
    # update Q 
    # ------------------------------------------------

    Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,

            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_nodes_idx,
            gens_with_loc_load_idx )

    

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]
             ) .* conj.(
            Igen )
    
    # gens_S  = uh[gens_nodes_idx] .* conj.(Igen )

    P_gens_update = [ idx ∈ slack_gens_nodes_idx ?
        real( gens_S[ n2s_gens_idx[idx]]) : P_gens[
            n2s_gens_idx[idx]]
                      for idx in gens_nodes_idx]

    Q_gens_update = [ imag( gens_S[ n2s_gens_idx[idx]]) 
                      for idx in gens_nodes_idx]

    # ----------------------------------------------- 

    current_mismatch =
        get_nodes_current_mismatch(
            vh,
            θh,
            P_gens_update,
            Q_gens_update,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,    
            ( gens_nodes_idx,
              non_gens_nodes_idx,
              gens_with_loc_load_idx,
               all_nodes_idx, # transformed_all_nodes_idx #
              ),
            ( n2s_gens_idx,
             n2s_non_gens_idx,
              n2s_gens_with_loc_load_idxs,
              n2s_all_nodes_idx #
              ),
            I_sum_ynj_vj;
            loc_load_exist =
                loc_load_exist)
    
    # -----------------------------------------------

    # power_mismatch =
    #     vh .* exp.(im * θh) .* conj.(current_mismatch)

    # red_P_mismatch =
    #     real.(power_mismatch[
    #         setdiff(transformed_all_nodes_idx,
    #                 transformed_slack_gens_nodes_idx)])

    # red_Q_mismatch =
    #     imag.(power_mismatch[
    #         setdiff(transformed_all_nodes_idx,
    #                 transformed_gens_nodes_idx)])

    # red_ΔPQ_x .=
    #     vcat(red_P_mismatch,
    #          red_Q_mismatch)

    red_real_current_mismatch =
        real.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)])

    red_imag_current_mismatch =
        imag.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)])

    red_ΔI_x .=
        vcat(red_real_current_mismatch,
             red_imag_current_mismatch)

    return nothing


end


function get_generic_sta_pf_ΔI_mismatch(
    red_ΔI_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para )

    (; loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs ) =
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :pf_kw_gens_vh_slack_θh_para,
              :pf_kw_net_para,
              :pf_kw_var_idxs,
              :pf_kw_PQ_para_idxs,
              :pf_kw_nodes_types_idxs,
              :pf_kw_n2s_idxs))

    
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             pf_kw_gens_vh_slack_θh_para,
             (:slack_gens_vh,
              :slack_gens_θh,

              :gens_vh,
              :non_slack_gens_vh ))
    
     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              pf_kw_net_para,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))

     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          NamedTupleTools.select(
              pf_kw_var_idxs,
              (:red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx ))

    
     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs ) =
          NamedTupleTools.select(
              pf_kw_PQ_para_idxs,
              (:P_gens_sta_para_Idxs,
               :Q_gens_sta_para_Idxs,
               :P_non_gens_sta_para_Idxs,
               :Q_non_gens_sta_para_Idxs,
               :P_g_loc_load_sta_para_Idxs,
               :Q_g_loc_load_sta_para_Idxs ))

     (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx ) =
         NamedTupleTools.select(
             pf_kw_nodes_types_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx ))    

     (; n2s_slack_gens_idx,
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


    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]


    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------
        
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    
    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]

    # -----------------------------------

    I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
        vh,
        θh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    # ------------------------------------------------   
    # update Q 
    # ------------------------------------------------

    Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
        get_gens_loc_load_current(
            vh,
            θh,           
            P_g_loc_load,
            Q_g_loc_load,

            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx, #

            gens_nodes_idx,
            gens_with_loc_load_idx )

    

    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]) .* conj.(
            Igen )
    
    # gens_S  = uh[gens_nodes_idx] .* conj.(Igen )

    P_gens_update = [ idx ∈ slack_gens_nodes_idx ?
        real( gens_S[ n2s_gens_idx[idx] ] ) : P_gens[
            n2s_gens_idx[idx]]
                      for idx in gens_nodes_idx]

    Q_gens_update = [ imag( gens_S[ n2s_gens_idx[idx]]) 
                      for idx in gens_nodes_idx]

    # ----------------------------------------------- 

    current_mismatch =
        get_nodes_current_mismatch(
            vh,
            θh,
            P_gens_update,
            Q_gens_update,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,    
            ( gens_nodes_idx,
              non_gens_nodes_idx,
              gens_with_loc_load_idx,
               all_nodes_idx, # transformed_all_nodes_idx #
              ),
            ( n2s_gens_idx,
             n2s_non_gens_idx,
              n2s_gens_with_loc_load_idxs,
              n2s_all_nodes_idx #
              ),
            I_sum_ynj_vj;
            loc_load_exist =
                loc_load_exist)
    
    # -----------------------------------------------

    # power_mismatch =
    #     vh .* exp.(im * θh) .* conj.(current_mismatch)

    # red_P_mismatch =
    #     real.(power_mismatch[
    #         setdiff(transformed_all_nodes_idx,
    #                 transformed_slack_gens_nodes_idx)])

    # red_Q_mismatch =
    #     imag.(power_mismatch[
    #         setdiff(transformed_all_nodes_idx,
    #                 transformed_gens_nodes_idx)])

    # red_ΔPQ_x .=
    #     vcat(red_P_mismatch,
    #          red_Q_mismatch)

    red_real_current_mismatch =
        real.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)])

    red_imag_current_mismatch =
        imag.(current_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)])

    red_ΔI_x .=
        vcat(red_real_current_mismatch,
             red_imag_current_mismatch)

    return nothing


end

# -----------------------------------------------

function get_a_model_integrated_pf_sta_ΔPQ_mismatch_generic(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para )

    (; loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs ) =
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :pf_kw_gens_vh_slack_θh_para,
              :pf_kw_net_para,
              :pf_kw_var_idxs,
              :pf_kw_PQ_para_idxs,
              :pf_kw_nodes_types_idxs,
              :pf_kw_n2s_idxs))

    
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             pf_kw_gens_vh_slack_θh_para,
             (:slack_gens_vh,
              :slack_gens_θh,

              :gens_vh,
              :non_slack_gens_vh ))
    
     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              pf_kw_net_para,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))

     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          NamedTupleTools.select(
              pf_kw_var_idxs,
              (:red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx ))

    
     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs ) =
          NamedTupleTools.select(
              pf_kw_PQ_para_idxs,
              (:P_gens_sta_para_Idxs,
               :Q_gens_sta_para_Idxs,
               :P_non_gens_sta_para_Idxs,
               :Q_non_gens_sta_para_Idxs,
               :P_g_loc_load_sta_para_Idxs,
               :Q_g_loc_load_sta_para_Idxs ))

     (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx ) =
         NamedTupleTools.select(
             pf_kw_nodes_types_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx ))    

     (; n2s_slack_gens_idx,
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


    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]


    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------

    
    # nodes_size = length( gens_nodes_idx ) +
    #     length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]

    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    # Note components_injection is assigned negative
    # value or positive value based on wether it is
    # a gen or non-gen
    
    P_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (P_non_gens[ n2s_non_gens_idx[
                nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]])  : 
            -(P_gens[ n2s_gens_idx[nth_idx]])        
        for nth_idx in all_nodes_idx ]


    Q_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (Q_non_gens[ n2s_non_gens_idx[ nth_idx]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(Q_gens[ n2s_gens_idx[ nth_idx ]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]]) :
            -(Q_gens[ n2s_gens_idx[nth_idx]])    
        for nth_idx in all_nodes_idx ]
    
    #-------------------------------    

    P_line_injection =
        get_nodes_real_uh_x_∑_ynj_x_vj(
            vh, θh; Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            all_nodes_idx,
            n2s_all_nodes_idx)

    Q_line_injection =
        get_nodes_imag_uh_x_∑_ynj_x_vj(
            vh, θh; Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            all_nodes_idx,
            n2s_all_nodes_idx)

    #-------------------------------    
    
    P_mismatch = P_line_injection .+
        P_nodal_components_injection
    
    Q_mismatch = Q_line_injection .+
        Q_nodal_components_injection
    
    # ----------------------------------------------- 

    red_P_mismatch =
        P_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)]

    red_Q_mismatch =
        Q_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)]

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)
    
    return nothing


end


function get_generic_sta_pf_ΔPQ_mismatch(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para )

    (; loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs ) =
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :pf_kw_gens_vh_slack_θh_para,
              :pf_kw_net_para,
              :pf_kw_var_idxs,
              :pf_kw_PQ_para_idxs,
              :pf_kw_nodes_types_idxs,
              :pf_kw_n2s_idxs))

    
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             pf_kw_gens_vh_slack_θh_para,
             (:slack_gens_vh,
              :slack_gens_θh,

              :gens_vh,
              :non_slack_gens_vh ))
    
     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              pf_kw_net_para,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))

     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          NamedTupleTools.select(
              pf_kw_var_idxs,
              (:red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx ))

    
     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs ) =
          NamedTupleTools.select(
              pf_kw_PQ_para_idxs,
              (:P_gens_sta_para_Idxs,
               :Q_gens_sta_para_Idxs,
               :P_non_gens_sta_para_Idxs,
               :Q_non_gens_sta_para_Idxs,
               :P_g_loc_load_sta_para_Idxs,
               :Q_g_loc_load_sta_para_Idxs ))

     (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx ) =
         NamedTupleTools.select(
             pf_kw_nodes_types_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx ))    

     (; n2s_slack_gens_idx,
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


    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]


    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------

    
    # nodes_size = length( gens_nodes_idx ) +
    #     length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  = red_vh_θh_x[
        red_vh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_idx2Idx ]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]

    #-------------------------------
    
    P_gens =
        pf_PQ_param[
            P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[
            Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[
            P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[
            Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    # Note components_injection is assigned negative
    # value or positive value based on wether it is
    # a gen or non-gen
    
    P_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (P_non_gens[ n2s_non_gens_idx[
                nth_idx ]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]])  : 
            -(P_gens[ n2s_gens_idx[nth_idx]])        
        for nth_idx in all_nodes_idx ]


    Q_nodal_components_injection = [
        nth_idx ∈ non_gens_nodes_idx ?
            (Q_non_gens[ n2s_non_gens_idx[ nth_idx]]) :
            nth_idx ∈ gens_with_loc_load_idx ?
            -(Q_gens[ n2s_gens_idx[ nth_idx ]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx]]) :
            -(Q_gens[ n2s_gens_idx[nth_idx]])    
        for nth_idx in all_nodes_idx ]
    
    #-------------------------------    

    P_line_injection =
        get_nodes_real_uh_x_∑_ynj_x_vj(
            vh, θh; Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            all_nodes_idx,
            n2s_all_nodes_idx)

    Q_line_injection =
        get_nodes_imag_uh_x_∑_ynj_x_vj(
            vh, θh; Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            all_nodes_idx,
            n2s_all_nodes_idx)

    #-------------------------------    
    
    P_mismatch = P_line_injection .+
        P_nodal_components_injection
    
    Q_mismatch = Q_line_injection .+
        Q_nodal_components_injection
    
    # ----------------------------------------------- 

    red_P_mismatch =
        P_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)]

    red_Q_mismatch =
        Q_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)]

    red_ΔPQ_x .=
        vcat(red_P_mismatch,
             red_Q_mismatch)
    
    return nothing


end

# -----------------------------------------------


function get_a_model_integrated_pf_sta_ΔPQ_mismatch(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_PQ_param
    ; pf_kw_para =
        pf_kw_para )

    #-------------------------------

    (; loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs ) =
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :pf_kw_gens_vh_slack_θh_para,
              :pf_kw_net_para,
              :pf_kw_var_idxs,
              :pf_kw_PQ_para_idxs,
              :pf_kw_nodes_types_idxs,
              :pf_kw_n2s_idxs))

    
    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             pf_kw_gens_vh_slack_θh_para,
             (:slack_gens_vh,
              :slack_gens_θh,

              :gens_vh,
              :non_slack_gens_vh ))
    
     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              pf_kw_net_para,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))

     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          NamedTupleTools.select(
              pf_kw_var_idxs,
              (:red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx ))

    
     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs ) =
          NamedTupleTools.select(
              pf_kw_PQ_para_idxs,
              (:P_gens_sta_para_Idxs,
               :Q_gens_sta_para_Idxs,
               :P_non_gens_sta_para_Idxs,
               :Q_non_gens_sta_para_Idxs,
               :P_g_loc_load_sta_para_Idxs,
               :Q_g_loc_load_sta_para_Idxs ))

     (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx ) =
         NamedTupleTools.select(
             pf_kw_nodes_types_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx ))    

     (; n2s_slack_gens_idx,
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


    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]


    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------
    
    P_gens =
        pf_PQ_param[ P_gens_sta_para_Idxs ]

    Q_gens =
        pf_PQ_param[ Q_gens_sta_para_Idxs ]

    P_non_gens  =
        pf_PQ_param[ P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[ Q_non_gens_sta_para_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[ P_g_loc_load_sta_para_Idxs ]

        Q_g_loc_load =
            pf_PQ_param[ Q_g_loc_load_sta_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
    
    #-------------------------------

    non_slack_gens_θh =
        red_vh_θh_x[
            red_non_slack_gens_θh_idx2Idx ]

    non_gens_vh  =
        red_vh_θh_x[ red_vh_Idxs ]

    non_gens_θh =
        red_vh_θh_x[
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
    
    # -----------------------------------
    # non slack gens active power mismatch
    # -----------------------------------
    
    if loc_load_exist == true

        non_slack_gens_P_mismatch = [
            nth_idx ∈ gens_with_loc_load_idx ? 
            (P_gens[ n2s_gens_idx[nth_idx] ] -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[
                nth_idx ] ] -
                vh[nth_idx] * sum( [
                    vh[idx] *
                        abs(ynj) *
                        cos(θh[nth_idx] - θh[idx] - angle(ynj))
                    for (ynj, idx) in
                        zip( Ynet[ nth_idx ],
                             nodes_idx_with_adjacent_nodes_idx[
                                 nth_idx])])) :
                                 (P_gens[ n2s_gens_idx[nth_idx] ] -
                vh[nth_idx] * sum( [
                    vh[idx] *
                        abs(ynj) *
                        cos(θh[nth_idx] - θh[idx] - angle(ynj))
                    for (ynj, idx) in
                        zip( Ynet[ nth_idx ],
                             nodes_idx_with_adjacent_nodes_idx[
                                 nth_idx])]))
            for nth_idx in
                non_slack_gens_nodes_idx ]
        
    else
    
    non_slack_gens_P_mismatch = [
        P_gens[ n2s_gens_idx[nth_idx] ] -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_idx_with_adjacent_nodes_idx[
                             nth_idx ])])
        for nth_idx in
            non_slack_gens_nodes_idx ]

    end

    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------
            
    non_gens_P_mismatch = [
        P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_idx_with_adjacent_nodes_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]

    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------
        
    non_gens_Q_mismatch = [
        Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_idx_with_adjacent_nodes_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]

    # ------------------------------------

    red_ΔPQ_x .=
        vcat(non_slack_gens_P_mismatch,
             non_gens_P_mismatch,
             non_gens_Q_mismatch)

    return nothing

end


# -----------------------------------------------
# -----------------------------------------------

# function get_a_model_integrated_pf_sta_ΔPQ_mismatch_generic2(
#     red_ΔPQ_x,
#     red_vh_θh_x,
#     pf_PQ_param
#     ; pf_kw_para =
#         pf_kw_para
#     )

#     (; loc_load_exist,
#      pf_kw_gens_vh_slack_θh_para,
#      pf_kw_net_para,
#      pf_kw_var_idxs,
#      pf_kw_PQ_para_idxs,
#      pf_kw_nodes_types_idxs,
#      pf_kw_n2s_idxs ) =
#          pf_kw_para
        
#     #-------------------------------

#     (; slack_gens_vh,
#      slack_gens_θh,

#      gens_vh,
#      non_slack_gens_vh ) =
#          pf_kw_gens_vh_slack_θh_para

#      (;Ynet,
#       nodes_idx_with_adjacent_nodes_idx) =
#           pf_kw_net_para

#      (; red_vh_Idxs,
#       red_non_slack_gens_θh_idx2Idx,
#       red_non_gens_θh_idx2Idx ) =
#           pf_kw_var_idxs
    
#      (; P_gens_sta_para_Idxs,
#       Q_gens_sta_para_Idxs,
#       P_non_gens_sta_para_Idxs,
#       Q_non_gens_sta_para_Idxs,
#       P_g_loc_load_sta_para_Idxs,
#       Q_g_loc_load_sta_para_Idxs ) =
#           pf_kw_PQ_para_idxs

#      (;slack_gens_nodes_idx,
#      non_slack_gens_nodes_idx,
#      gens_nodes_idx,
#      non_gens_nodes_idx,
#      gens_with_loc_load_idx,
#      all_nodes_idx) =
#           pf_kw_nodes_types_idxs    

#      (; n2s_slack_gens_idx,
#      n2s_non_slack_gens_idx,
#      n2s_gens_idx,
#      n2s_non_gens_idx,
#      n2s_gens_with_loc_load_idxs,
#      n2s_all_nodes_idx ) =
#          pf_kw_n2s_idxs
    
#     #-------------------------------
    
#     P_gens =
#         pf_PQ_param[ P_gens_sta_para_Idxs ]

#     Q_gens =
#         pf_PQ_param[ Q_gens_sta_para_Idxs ]

#     P_non_gens  =
#         pf_PQ_param[ P_non_gens_sta_para_Idxs ]

#     Q_non_gens = 
#         pf_PQ_param[ Q_non_gens_sta_para_Idxs ]
    
#     if loc_load_exist == true

#         P_g_loc_load =
#             pf_PQ_param[ P_g_loc_load_sta_para_Idxs ]

#         Q_g_loc_load =
#             pf_PQ_param[ Q_g_loc_load_sta_para_Idxs ]
        
#     else

#         P_g_loc_load = [0.0]
            
#         Q_g_loc_load = [0.0]
#     end
    
#     #-------------------------------

#     nodes_size = length( gens_nodes_idx ) +
#         length( non_gens_nodes_idx )
    
#     #-------------------------------

#     non_slack_gens_θh = red_vh_θh_x[
#         red_non_slack_gens_θh_idx2Idx ]

#     non_gens_vh  = red_vh_θh_x[
#         red_vh_Idxs ]

#     non_gens_θh = red_vh_θh_x[
#         red_non_gens_θh_idx2Idx ]

#     vh = [ idx ∈ gens_nodes_idx ?
#             gens_vh[
#                 n2s_gens_idx[ idx ] ]  :
#                     non_gens_vh[
#                         n2s_non_gens_idx[ idx ]]
#                for idx in all_nodes_idx ]

#     θh = [ idx ∈ slack_gens_nodes_idx ?
#             slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
#             idx ∈ non_slack_gens_nodes_idx ?
#             non_slack_gens_θh[
#                 n2s_non_slack_gens_idx[ idx] ] :
#                     non_gens_θh[
#                         n2s_non_gens_idx[ idx ]]
#         for idx in all_nodes_idx ]


#     # -----------------------------------
#     # alternative
#     # -----------------------------------

#     # uh = vh .* exp.(im * θh)

#     # -----------------------------------

#     I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
#         vh,
#         θh,
#         Ynet,
#         nodes_idx_with_adjacent_nodes_idx )

#     # ------------------------------------------------ 
#     # update Q 
#     # ------------------------------------------------

#     Igen = I_sum_ynj_vj[gens_nodes_idx] +
#         get_gens_loc_load_current(
#             vh,
#             θh,           
#             P_g_loc_load,
#             Q_g_loc_load,

#             n2s_gens_with_loc_load_idxs,

#             gens_nodes_idx,
#             gens_with_loc_load_idx )

    

#     gens_S  = vh[gens_nodes_idx] .*
#         exp.(im * θh[gens_nodes_idx]) .* conj.(Igen )
    
#     # gens_S  = uh[gens_nodes_idx] .* conj.(Igen )

#     P_gens_update = [ idx ∈ slack_gens_nodes_idx ?
#         real( gens_S[ n2s_gens_idx[idx] ] ) : P_gens[
#             n2s_gens_idx[idx]]
#                       for idx in gens_nodes_idx]

#     Q_gens_update = [ imag( gens_S[ n2s_gens_idx[idx]) 
#                       for idx in gens_nodes_idx]

#     # ------------------------------------------------ 

#     current_mismatch =
#         get_nodes_current_mismatch(
#             vh,
#             θh,
#             P_gens_update,
#             Q_gens_update,
#             P_non_gens,
#             Q_non_gens,
#             P_g_loc_load,
#             Q_g_loc_load,    
#             ( gens_nodes_idx,
#               non_gens_nodes_idx,
#               gens_with_loc_load_idx,
#               all_nodes_idx),
#             ( n2s_gens_idx,
#              n2s_non_gens_idx,
#              n2s_gens_with_loc_load_idxs ),
#             I_sum_ynj_vj;
#             loc_load_exist =
#                 loc_load_exist)
    
#     # -----------------------------------------------

#     power_mismatch =
#         vh .* exp.(im * θh) .* conj.(current_mismatch)

#     red_P_mismatch =
#         real.(power_mismatch[
#             setdiff(all_nodes_idx,
#                     slack_gens_nodes_idx)])

#     red_Q_mismatch =
#         imag.(power_mismatch[
#             setdiff(all_nodes_idx,
#                     gens_nodes_idx)])

#     red_ΔPQ_x .=
#         vcat(red_P_mismatch,
#              red_Q_mismatch)
#     return nothing
# end


# -----------------------------------------------
# Others
# -----------------------------------------------

function get_a_model_integrated_dyn_pf_full_ΔPQ_mismatch(
    full_ΔPQ,
    full_vh_θh,
    integ_param
    ;full_dyn_pf_fun_kwd_para  =
        full_dyn_pf_fun_kwd_para   )

    #-------------------------------

    (;
     full_kwd_para,
     dyn_pf_fun_kwd_para ) =
         full_dyn_pf_fun_kwd_para

    (;full_vars_Idxs,
     full_gens_id_iq ) =
         full_kwd_para

    #-------------------------------

   (;full_gens_vh_Idxs,
    full_gens_θh_Idxs,
    full_non_gens_nodes_vh_Idxs,
    full_non_gens_nodes_θh_Idxs) =
        full_vars_Idxs

    #-------------------------------
    
   (;gens_i_d_0,
     gens_i_q_0 ) =
        full_gens_id_iq 
    
    #-------------------------------
    
    (;loc_load_exist,
     dyn_pf_fun_kwd_nll_para_vars_Idxs,
     dyn_pf_fun_kwd_wll_para_vars_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
     gens_nodes_ra_Xd_dash_Xq_dash) =
         dyn_pf_fun_kwd_para
    
    #-------------------------------
    #-------------------------------

   (;P_gens_dyn_para_Idxs,
     Q_gens_dyn_para_Idxs,
     P_non_gens_dyn_para_Idxs,
     Q_non_gens_dyn_para_Idxs,
     δ_ed_eq_pf_dyn_para_Idxs,
     P_g_loc_load_dyn_para_Idxs,
     Q_g_loc_load_dyn_para_Idxs                   
    ) =
        dyn_pf_fun_kwd_wll_para_vars_Idxs 


   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 

    #-------------------------------
    
   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------
    #-------------------------------
    #-------------------------------
    
    P_gens =
        integ_param[ P_gens_dyn_para_Idxs ]

    Q_gens =
        integ_param[ Q_gens_dyn_para_Idxs ]

    P_non_gens  =
        integ_param[ P_non_gens_dyn_para_Idxs ]

    Q_non_gens = 
        integ_param[ Q_non_gens_dyn_para_Idxs ]

    flat_δ_ω_ed_dash_eq_dash =
        integ_param[ δ_ed_eq_pf_dyn_para_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            integ_param[ P_g_loc_load_dyn_para_Idxs ]

        Q_g_loc_load =
            integ_param[ Q_g_loc_load_dyn_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    #-------------------------------

    if length(flat_δ_ω_ed_dash_eq_dash) !=
        length(δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed )

        δ_ω_ed_dash_eq_dash =
            [ flat_δ_ω_ed_dash_eq_dash[idx]
              for idx in
                  δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed ]

    else

        δ_ω_ed_dash_eq_dash =
            flat_δ_ω_ed_dash_eq_dash

    end
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
        
    #-------------------------------
    #-------------------------------

    full_gens_vh =
        full_vh_θh[ full_gens_vh_Idxs ]

    full_gens_θh =
        full_vh_θh[ full_gens_θh_Idxs ]

    non_gens_vh = full_non_gens_nodes_vh =
        full_vh_θh[ full_non_gens_nodes_vh_Idxs ]

    non_gens_θh  = full_non_gens_nodes_θh =
        full_vh_θh[ full_non_gens_nodes_θh_Idxs ]

    vh = [
        idx ∈ gens_nodes_idx ?
            full_gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    full_non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ gens_nodes_idx ?
            full_gens_θh[
                n2s_gens_idx[ idx] ] :
                    full_non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    # -------------------------------------

    gens_δ = first.( δ_ω_ed_dash_eq_dash )

    gens_ed_dash = third.( δ_ω_ed_dash_eq_dash )

    gens_eq_dash = fourth.( δ_ω_ed_dash_eq_dash )
    
    # -------------------------------------

    gens_ra = first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash = second.( gens_nodes_ra_Xd_dash_Xq_dash  )

    gens_Xq_dash = third.( gens_nodes_ra_Xd_dash_Xq_dash  ) 

    
    # -------------------------------------

    gens_i_d = gens_i_d_0

    gens_i_q = gens_i_q_0
    
    # -----------------------------------
    # gens active power mismatch
    # -----------------------------------
    
    if loc_load_exist == true

        gens_P_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ] -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) ) :
                                 ( (gens_i_d[ n2s_gens_idx[nth_idx] ] *
                                 vh[nth_idx] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) )
                                 for nth_idx in
                gens_nodes_idx ]  

    else

        gens_P_mismatch = [ (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) for nth_idx in
                gens_nodes_idx ]  
        
        
    end

    # -----------------------------------
    # gens reactive power mismatch
    # -----------------------------------
    
    if loc_load_exist == true
        
        gens_Q_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ nth_idx ] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx] ] - 
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj) )
                for (ynj, idx) in
                    zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])]) ) :
                                 ((gens_i_d[ n2s_gens_idx[nth_idx]] * vh[nth_idx] *
                                 cos(gens_δ[n2s_gens_idx[nth_idx]]-θh[nth_idx])) -
                                 (gens_i_q[ n2s_gens_idx[nth_idx]] * vh[ nth_idx] *
                                 sin(gens_δ[n2s_gens_idx[nth_idx]]-θh[nth_idx])) -
                                         vh[nth_idx] *
                                         sum([vh[idx] * abs(ynj) *
                                         sin(θh[nth_idx]-θh[idx]-angle(ynj))
                                              for (ynj, idx) in
                                                  zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])]) ) 
                            for nth_idx in
                                gens_nodes_idx]        
        
    else

        gens_Q_mismatch = [ (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ nth_idx ] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])])
                            for nth_idx in
                                gens_nodes_idx ]
        
    end

    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------
        
    non_gens_P_mismatch = [
        P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------
        
    non_gens_Q_mismatch = [
        Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------

    full_ΔPQ .=
        vcat(gens_P_mismatch,
             non_gens_P_mismatch,
             gens_Q_mismatch,                 
             non_gens_Q_mismatch )
    
    # # -----------------------------------
    # # gens active power mismatch
    # # -----------------------------------

    # gens_P_mismatch =
    #     [  nth_idx ∈ gens_nodes_with_loc_loads_idx ?
    #     get_nth_gen_P_pf_mismatch_wt_loc_load(
    #         vh,
    #         θh,    
    #         gens_δ[n2s_gens_idx[nth_idx]],
    #         gens_i_d[n2s_gens_idx[nth_idx]],
    #         gens_i_q[n2s_gens_idx[nth_idx]],  
    #         Ynet[ nth_idx ] ,
    #         nodes_node_idx_and_incident_edges_other_node_idx[
    #             nth_idx ],
    #         P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ],
    #         nth_idx ) :
    #             get_nth_gen_P_pf_mismatch_no_loc_load(
    #                 vh,
    #                 θh,    
    #                 gens_δ[n2s_gens_idx[nth_idx]],
    #                 gens_i_d[n2s_gens_idx[nth_idx]],
    #                 gens_i_q[n2s_gens_idx[nth_idx]],  
    #                 Ynet[ nth_idx ] ,
    #                 nodes_node_idx_and_incident_edges_other_node_idx[
    #                     nth_idx ],
    #                 nth_idx )

    #      for nth_idx in
    #             gens_nodes_idx ]

    # # -----------------------------------
    # # gens reactive power mismatch
    # # -----------------------------------

    # gens_Q_mismatch =
    #     [ nth_idx ∈ gens_nodes_with_loc_loads_idx ?
    #     get_nth_gen_Q_pf_mismatch_wt_loc_load(
    #         vh,
    #         θh,    
    #         gens_δ[n2s_gens_idx[nth_idx]],
    #         gens_i_d[n2s_gens_idx[nth_idx]],
    #         gens_i_q[n2s_gens_idx[nth_idx]],  
    #         Ynet[ nth_idx ] ,
    #         nodes_node_idx_and_incident_edges_other_node_idx[
    #             nth_idx ],
    #         Q_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ],
    #         nth_idx ) :
    #             get_nth_gen_Q_pf_mismatch_no_loc_load(    
    #                 vh,
    #                 θh,    
    #                 gens_δ[n2s_gens_idx[nth_idx]],
    #                 gens_i_d[n2s_gens_idx[nth_idx]],
    #                 gens_i_q[n2s_gens_idx[nth_idx]],  
    #                 Ynet[ nth_idx ] ,
    #                 nodes_node_idx_and_incident_edges_other_node_idx[
    #                     nth_idx ],
    #                 nth_idx )

    #      for nth_idx in
    #             gens_nodes_idx ]

    # # -----------------------------------
    # # non_gens active power mismatch        
    # # -----------------------------------

    # non_gens_P_mismatch =
    #     [ get_nth_non_gen_P_pf_mismatch( 
    #         vh,
    #         θh,           
    #         Ynet[ nth_idx ],
    #         nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ],
    #         P_non_gens[ n2s_non_gens_idx[ nth_idx ] ],
    #         nth_idx )

    #       for nth_idx in
    #             non_gens_nodes_idx ]

    # # ------------------------------------
    # # non_gens reactive power mismatch                
    # # ------------------------------------

    # non_gens_Q_mismatch =
    #     [ get_nth_non_gen_Q_pf_mismatch( 
    #         vh,
    #         θh,           
    #         Ynet[ nth_idx ],
    #         nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ],
    #         Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ],
    #         nth_idx )

    #       for nth_idx in
    #             non_gens_nodes_idx ]

    # # ------------------------------------

    # full_ΔPQ .=
    #     vcat(gens_P_mismatch,
    #          non_gens_P_mismatch,
    #          gens_Q_mismatch,                 
    #          non_gens_Q_mismatch )

    return nothing

end



function get_a_model_integrated_dyn_pf_intg_ΔPQ_Δidq_mismatch(
    intg_ΔPQ_id_iq, intg_vh_θh_id_iq, integ_param;
    intg_dyn_pf_fun_kwd_para =
        intg_dyn_pf_fun_kwd_para  )

    #-------------------------------

    (;
     intg_vars_Idxs,
     dyn_pf_fun_kwd_para ) =
         intg_dyn_pf_fun_kwd_para

    #-------------------------------
    
   (;intg_gens_vh_Idxs,
    intg_gens_θh_Idxs,
    intg_non_gens_nodes_vh_Idxs,
    intg_non_gens_nodes_θh_Idxs,
    intg_gen_id_Idxs,
    intg_gen_iq_Idxs) =
        intg_vars_Idxs

    #-------------------------------    
    
    (;loc_load_exist,
     dyn_pf_fun_kwd_nll_para_vars_Idxs,
     dyn_pf_fun_kwd_wll_para_vars_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         dyn_pf_fun_kwd_para
    
    #-------------------------------

   (;P_gens_dyn_para_Idxs,
     Q_gens_dyn_para_Idxs,
     P_non_gens_dyn_para_Idxs,
     Q_non_gens_dyn_para_Idxs,
     δ_ed_eq_pf_dyn_para_Idxs,
     P_g_loc_load_dyn_para_Idxs,
     Q_g_loc_load_dyn_para_Idxs                   
    ) =
        dyn_pf_fun_kwd_wll_para_vars_Idxs 
    
    #-------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #-------------------------------

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------

   # (;
   #  n2s_slack_gens_idx,
   #  n2s_non_slack_gens_idx,
   #  n2s_gens_idx,
   #  n2s_non_gens_idx,
   #  n2s_gens_with_loc_load_idxs,
   #  n2s_all_nodes_idx ) =
   #      dyn_pf_fun_kwd_net_idxs

    
    #-------------------------------
    #-------------------------------
    
    P_gens =
        integ_param[ P_gens_dyn_para_Idxs ]

    Q_gens =
        integ_param[ Q_gens_dyn_para_Idxs ]

    P_non_gens  =
        integ_param[ P_non_gens_dyn_para_Idxs ]

    Q_non_gens = 
        integ_param[ Q_non_gens_dyn_para_Idxs ]

    flat_δ_ω_ed_dash_eq_dash =
        integ_param[ δ_ed_eq_pf_dyn_para_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            integ_param[ P_g_loc_load_dyn_para_Idxs ]

        Q_g_loc_load =
            integ_param[ Q_g_loc_load_dyn_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #-------------------------------

    if length(flat_δ_ω_ed_dash_eq_dash) !=
        length(δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed )

        δ_ω_ed_dash_eq_dash =
            [ flat_δ_ω_ed_dash_eq_dash[idx]
              for idx in
                  δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed ]

    else

        δ_ω_ed_dash_eq_dash =
            flat_δ_ω_ed_dash_eq_dash

    end
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
        
    #-------------------------------

    intg_gens_vh =
        intg_vh_θh_id_iq[ intg_gens_vh_Idxs ]

    intg_gens_θh =
        intg_vh_θh_id_iq[ intg_gens_θh_Idxs ]

    intg_non_gens_nodes_vh =
        intg_vh_θh_id_iq[ intg_non_gens_nodes_vh_Idxs ]

    intg_non_gens_nodes_θh =
        intg_vh_θh_id_iq[ intg_non_gens_nodes_θh_Idxs ]

    gens_i_d =
        intg_vh_θh_id_iq[ intg_gen_id_Idxs ]

    gens_i_q =
        intg_vh_θh_id_iq[ intg_gen_iq_Idxs ]


    vh = [
        idx ∈ gens_nodes_idx ?
            intg_gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    intg_non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ gens_nodes_idx ?
            intg_gens_θh[
                n2s_gens_idx[ idx] ] :
                    intg_non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    # -------------------------------------

    gens_δ = first.( δ_ω_ed_dash_eq_dash )

    gens_ed_dash = third.( δ_ω_ed_dash_eq_dash )

    gens_eq_dash = fourth.( δ_ω_ed_dash_eq_dash )

    # -------------------------------------

    gens_ra = first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash = second.( gens_nodes_ra_Xd_dash_Xq_dash  )

    gens_Xq_dash = third.( gens_nodes_ra_Xd_dash_Xq_dash  )        

    # -----------------------------------
    # gens real part stator equation mismatch
    # -----------------------------------

    # gens_stator_equations_mismatch_real =
    #     [ get_a_gen_real_stator_equations_mismatch(
    #         a_vh, a_θh,
    #         a_δ, a_ed_dash, a_eq_dash,
    #         a_ra, a_X_d_dash, a_X_q_dash,
    #         a_id, a_iq )
    #       for ( a_vh, a_θh,
    #             a_δ, a_ed_dash, a_eq_dash,
    #             a_ra, a_X_d_dash, a_X_q_dash,
    #             a_id, a_iq ) in
    #           zip( intg_gens_vh, intg_gens_θh,
    #                gens_δ, gens_ed_dash, gens_eq_dash,
    #                gens_ra, gens_Xd_dash, gens_Xq_dash,
    #                gens_i_d, gens_i_q ) ]

    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            vh[idx] * sin(gens_δ[n2s_gens_idx[idx]] - θh[idx]) -
            gens_ra[n2s_gens_idx[idx]] * gens_i_d[n2s_gens_idx[idx]] +
            gens_Xq_dash[ n2s_gens_idx[idx]] * gens_i_q[ n2s_gens_idx[idx]]
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    # gens imag part stator equation mismatch
    # -----------------------------------

    # gens_stator_equations_mismatch_imag =
    #     [ get_a_gen_imag_stator_equations_mismatch(
    #         a_vh, a_θh,
    #         a_δ, a_ed_dash, a_eq_dash,
    #         a_ra, a_X_d_dash, a_X_q_dash,
    #         a_id, a_iq )
    #       for ( a_vh, a_θh,
    #             a_δ, a_ed_dash, a_eq_dash,
    #             a_ra, a_X_d_dash, a_X_q_dash,
    #             a_id, a_iq ) in
    #           zip( intg_gens_vh, intg_gens_θh,
    #                gens_δ, gens_ed_dash, gens_eq_dash,
    #                gens_ra, gens_Xd_dash, gens_Xq_dash,
    #                gens_i_d, gens_i_q ) ]


    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] -
            vh[ idx] * cos(gens_δ[ n2s_gens_idx[idx]] - θh[idx]) -
            gens_ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[idx]] -
            gens_Xd_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]]
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    # gens active power mismatch
    # -----------------------------------


    # gens_P_mismatch =
    #     [  nth_idx ∈ gens_nodes_with_loc_loads_idx ?
    #     get_nth_gen_P_pf_mismatch_wt_loc_load(
    #         vh,
    #         θh,    
    #         gens_δ[ n2s_gens_idx[nth_idx] ],
    #         gens_i_d[n2s_gens_idx[nth_idx]],
    #         gens_i_q[n2s_gens_idx[nth_idx]],  
    #         Ynet[ nth_idx ] ,
    #         nodes_node_idx_and_incident_edges_other_node_idx[
    #             nth_idx ],
    #         P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ],
    #         nth_idx ) :
    #             get_nth_gen_P_pf_mismatch_no_loc_load(
    #                 vh,
    #                 θh,    
    #                 gens_δ[n2s_gens_idx[nth_idx]],
    #                 gens_i_d[n2s_gens_idx[nth_idx]],
    #                 gens_i_q[n2s_gens_idx[nth_idx]],  
    #                 Ynet[ nth_idx ] ,
    #                 nodes_node_idx_and_incident_edges_other_node_idx[
    #                     nth_idx ],
    #                 nth_idx )

    #      for nth_idx in
    #             gens_nodes_idx ]
    
    if loc_load_exist == true

        gens_P_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ] -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) ) :
                                 ( (gens_i_d[ n2s_gens_idx[ nth_idx] ] *
                                 vh[nth_idx] *
                                 sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
                                 (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
                                 cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
                                 vh[nth_idx] * sum( [
                                     vh[idx] *
                                         abs(ynj) *
                                         cos(θh[nth_idx] - θh[idx] - angle(ynj))
                                     for (ynj, idx) in
                                         zip( Ynet[ nth_idx ],
                                              nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ])]) )

            for nth_idx in gens_nodes_idx ]  

    else

        gens_P_mismatch = [ (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) for nth_idx in
                gens_nodes_idx ]  
        
        
    end

    # -----------------------------------
    # gens reactive power mismatch
    # -----------------------------------


    # gens_Q_mismatch =
    #     [ nth_idx ∈ gens_nodes_with_loc_loads_idx ?
    #     get_nth_gen_Q_pf_mismatch_wt_loc_load(
    #         vh,
    #         θh,    
    #         gens_δ[n2s_gens_idx[nth_idx]],
    #         gens_i_d[n2s_gens_idx[nth_idx]],
    #         gens_i_q[n2s_gens_idx[nth_idx]],  
    #         Ynet[ nth_idx ] ,
    #         nodes_node_idx_and_incident_edges_other_node_idx[
    #             nth_idx ],
    #         Q_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ],
    #         nth_idx ) :
    #             get_nth_gen_Q_pf_mismatch_no_loc_load(    
    #                 vh,
    #                 θh,    
    #                 gens_δ[n2s_gens_idx[nth_idx]],
    #                 gens_i_d[n2s_gens_idx[nth_idx]],
    #                 gens_i_q[n2s_gens_idx[nth_idx]],  
    #                 Ynet[ nth_idx ] ,
    #                 nodes_node_idx_and_incident_edges_other_node_idx[
    #                     nth_idx ],
    #                 nth_idx )

    #      for nth_idx in
    #             gens_nodes_idx ]
    #
    
    if loc_load_exist == true
        
        gens_Q_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ nth_idx ] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx] ] - 
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj) )
                for (ynj, idx) in
                    zip( Ynet[ nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])]) ) :
                                 ((gens_i_d[ n2s_gens_idx[ nth_idx]] *
                                 vh[nth_idx] *
                                 cos(gens_δ[ n2s_gens_idx[ nth_idx]] - θh[nth_idx])) -
                                 (gens_i_q[ n2s_gens_idx[ nth_idx]] * vh[ nth_idx] *
                                 sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[nth_idx])) -
                                         vh[nth_idx] *
                                         sum([vh[idx] * abs(ynj) *
                                         sin(θh[nth_idx]-θh[idx]-angle(ynj))
                                              for (ynj, idx) in
                                                  zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])]) ) 
                            for nth_idx in
                                gens_nodes_idx]        
        
    else

        gens_Q_mismatch = [ (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ nth_idx ] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])])
                            for nth_idx in
                                gens_nodes_idx ]
        
    end
    
    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------

    # non_gens_P_mismatch =
    #     [ get_nth_non_gen_P_pf_mismatch( 
    #         vh,
    #         θh,           
    #         Ynet[ nth_idx ],
    #         nodes_node_idx_and_incident_edges_other_node_idx[
    #             nth_idx ], P_non_gens[ n2s_non_gens_idx[ nth_idx ] ],
    #         nth_idx )

    #       for nth_idx in
    #             non_gens_nodes_idx ]
        
    non_gens_P_mismatch = [
        P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------

    # non_gens_Q_mismatch =
    #     [ get_nth_non_gen_Q_pf_mismatch( 
    #         vh,
    #         θh,           
    #         Ynet[ nth_idx ],
    #         nodes_node_idx_and_incident_edges_other_node_idx[
    #             nth_idx ],  Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ],
    #         nth_idx )

    #       for nth_idx in
    #             non_gens_nodes_idx ]

        
    non_gens_Q_mismatch = [
        Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------

    intg_ΔPQ_id_iq .=
        vcat(gens_P_mismatch,
             non_gens_P_mismatch,
             gens_Q_mismatch,                 
             non_gens_Q_mismatch,
             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag)

    return nothing

end


#########################################################
#########################################################


"""

Ig_i_loc_load_real_B = ( id_i * sin(δ_i) + iq_i * cos(δ_i) )

Ig_i_loc_load_imag_B = ( iq_i * sin(δ_i) - id_i * cos(δ_i) )


Ig_i_loc_load_real_A =
    ( Ploc_i * cos(θ_i) + Qloc_i * sin(θ_i) ) / vh_i 


Ig_i_loc_load_imag_A =
    ( Ploc_i * sin(θ_i) - Qloc_i * cos(θ_i)) / vh_i

Ig_i_loc_load_real_C = sum( vh_k * Y_ik * cos( θ_k + α_ik  ) )

Ig_i_loc_load_real_C = -sum( vh_k * Y_ik * sin( θ_k + α_ik  ) )

"""
function get_a_model_integrated_intra_dyn_current_balance_pf_ΔI_idq_mismatch(
    ΔI_id_iq, vh_θh_id_iq, intra_dyn_pf_mismatch_flat_para;
    intra_dyn_pf_mismatch_kwd_para =
        intra_dyn_pf_mismatch_kwd_para  )
    
    #-------------------------------

    (;
     loc_load_exist,
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     non_pre_ordered_dyn_pf_vars_Idxs) =
         intra_dyn_pf_mismatch_kwd_para 

    #-------------------------------

  (;
   gens_vh_idxs,
   gens_θh_idxs,

   non_gens_nodes_vh_idxs,
   non_gens_nodes_θh_idxs,

   gens_id_idxs,
    gens_iq_idxs,
    
    flat_vh_flat_θh_flat_id_iq_Idx) =
        non_pre_ordered_dyn_pf_vars_Idxs 

    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        dyn_pf_vars_Idxs

    #-------------------------------    
    
    (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
    
    #-------------------------------

   (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) =
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs 
    
    #-------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #-------------------------------

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------
    
    P_gens =
        intra_dyn_pf_mismatch_flat_para[
            dyn_P_gens_Idxs ]

    Q_gens =
        intra_dyn_pf_mismatch_flat_para[
            dyn_Q_gens_Idxs ]

    P_non_gens  =
        intra_dyn_pf_mismatch_flat_para[
            dyn_P_non_gens_Idxs ]

    Q_non_gens = 
        intra_dyn_pf_mismatch_flat_para[
            dyn_Q_non_gens_Idxs ]

    flat_δ_ed_dash_eq_dash =
        intra_dyn_pf_mismatch_flat_para[
            dyn_δ_ed_eq_pf_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            intra_dyn_pf_mismatch_flat_para[
                dyn_P_g_loc_load_Idxs ]

        Q_g_loc_load =
            intra_dyn_pf_mismatch_flat_para[
                dyn_Q_g_loc_load_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #-------------------------------

    if length(flat_δ_ed_dash_eq_dash) !=
        length( δ_ed_dash_eq_dash_Idxs_in_flattend )

        δ_ed_dash_eq_dash =
            [ flat_δ_ed_dash_eq_dash[idx]
              for idx in
                  δ_ed_dash_eq_dash_Idxs_in_flattend ]

    else

        δ_ed_dash_eq_dash =
            flat_δ_ed_dash_eq_dash

    end
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
        
    #-------------------------------

    # gens_vh =
    #     vh_θh_id_iq[ gens_vh_Idxs ]

    # gens_θh =
    #     vh_θh_id_iq[ gens_θh_Idxs ]

    # non_gens_nodes_vh =
    #     vh_θh_id_iq[ non_gens_nodes_vh_Idxs ]

    # non_gens_nodes_θh =
    #     vh_θh_id_iq[ non_gens_nodes_θh_Idxs ]

    # gens_i_d =
    #     vh_θh_id_iq[ gen_id_Idxs ]

    # gens_i_q =
    #     vh_θh_id_iq[ gen_iq_Idxs ]
    
    #-------------------------------


    # gens_vh =
    #     vh_θh_id_iq[ gens_vh_idxs ]

    # gens_θh =
    #     vh_θh_id_iq[ gens_θh_idxs ]

    # non_gens_nodes_vh =
    #     vh_θh_id_iq[ non_gens_nodes_vh_idxs ]

    # non_gens_nodes_θh =
    #     vh_θh_id_iq[ non_gens_nodes_θh_idxs ]

    # gens_i_d =
    #     vh_θh_id_iq[ gens_id_idxs ]

    # gens_i_q =
    #     vh_θh_id_iq[ gens_iq_idxs ]

    # #-------------------------------
    
    flat_vh_Idx, flat_θh_Idx, flat_id_Idx, flat_iq_Idx =
        flat_vh_flat_θh_flat_id_iq_Idx
    

    flat_vh =
        vh_θh_id_iq[ flat_vh_Idx ]

    flat_θh =
        vh_θh_id_iq[ flat_θh_Idx ]

    gens_i_d =
        vh_θh_id_iq[ flat_id_Idx ]

    gens_i_q =
        vh_θh_id_iq[ flat_iq_Idx ]
    

    gens_vh =
        flat_vh[ gens_nodes_idx ]

    gens_θh =
        flat_θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[ non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[ non_gens_nodes_idx ]
    
    #-------------------------------

    
    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ gens_nodes_idx ?
            gens_θh[
                n2s_gens_idx[ idx] ] :
                    non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    # -------------------------------------

    gens_δ = first.( δ_ed_dash_eq_dash )

    gens_ed_dash = second.( δ_ed_dash_eq_dash )

    gens_eq_dash = third.( δ_ed_dash_eq_dash )

    # -------------------------------------

    gens_ra = first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash = second.( gens_nodes_ra_Xd_dash_Xq_dash  )

    gens_Xq_dash = third.( gens_nodes_ra_Xd_dash_Xq_dash  )        

    # -----------------------------------
    # gens real part stator equation mismatch
    # -----------------------------------

    """ ed_dash - vh * sin(δ - θ) - ra * id  + Xq_dash * iq """
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            vh[idx] * sin(gens_δ[n2s_gens_idx[idx]] - θh[idx]) -
            gens_ra[n2s_gens_idx[idx]] * gens_i_d[n2s_gens_idx[idx]] +
            gens_Xq_dash[ n2s_gens_idx[idx]] * gens_i_q[ n2s_gens_idx[idx]]
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    # gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θ) - ra * iq  - Xd_dash * id """
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] -
            vh[ idx] * cos(gens_δ[ n2s_gens_idx[idx]] - θh[idx]) -
            gens_ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[idx]] -
            gens_Xd_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]]
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    # gens nodes  current real part mismatch
    # -----------------------------------
    
    """
    id_i  * sin(δ_i) + iq_i * cos(δ_i) + ( Ploc_i * cos(θ_i) + Qloc_i * sin(θ_i) ) / vh_i =
        ∑( vh_k * Y_ik * cos(θ_k + β_ik ) )

    """
    if loc_load_exist == true

        gens_I_real_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * sin( gens_δ[ n2s_gens_idx[nth_idx] ] )) +
                (gens_i_q[ n2s_gens_idx[nth_idx] ] * cos(gens_δ[ n2s_gens_idx[nth_idx] ] )) )  -
                ( P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ] * cos( θh[nth_idx] ) +
                Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] * sin(θh[nth_idx]))/vh[nth_idx] -
                sum([ vh[idx] * abs(ynj) * cos(θh[idx] + angle(ynj))
                      for (ynj, idx) in
                          zip( Ynet[ nth_idx ],
                               nodes_node_idx_and_incident_edges_other_node_idx[nth_idx])])) :
                ( ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * sin( gens_δ[n2s_gens_idx[nth_idx]])) +
                (gens_i_q[ n2s_gens_idx[nth_idx] ] * cos(gens_δ[ n2s_gens_idx[nth_idx] ] )) ) -
                sum([ vh[idx] * abs(ynj) * cos(θh[idx] + angle(ynj))
                       for (ynj, idx) in
                           zip( Ynet[ nth_idx ],
                                nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ])]) )
            for nth_idx in gens_nodes_idx ]  

    else

        gens_I_real_mismatch  = [
            ( ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * sin( gens_δ[ n2s_gens_idx[nth_idx] ] )) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * cos(gens_δ[ n2s_gens_idx[nth_idx] ] )) ) -
            sum([ vh[idx] * abs(ynj) * cos(θh[idx] + angle(ynj))
                   for (ynj, idx) in
                       zip( Ynet[ nth_idx ],
                            nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ])]) )
            for nth_idx in gens_nodes_idx ]  
    end

    # -----------------------------------
    # gens nodes  current imag part mismatch
    # -----------------------------------

    """
    -(id_i * cos(δ_i) - iq_i * sin(δ_i) )  + ( Ploc_i * sin(θ_i) - Qloc_i * cos(θ_i) ) / vh_i =
        -∑( vh_k * Y_ik * sin(θ_k + β_ik ) )

    """
    if loc_load_exist == true 
        
        gens_I_imag_mismatch  = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( -( (gens_i_d[ n2s_gens_idx[nth_idx] ] * cos(gens_δ[ n2s_gens_idx[nth_idx] ] ) ) -
                (gens_i_q[ n2s_gens_idx[nth_idx] ] *  sin(gens_δ[ n2s_gens_idx[nth_idx] ] ))) -
                ( P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ] * sin( θh[nth_idx] ) -
                Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] * cos(θh[nth_idx]))/vh[nth_idx] - 
                sum( [ vh[idx] * abs(ynj) * sin( θh[idx] + angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ nth_idx],
                                nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx])]) ) :
                ( -( (gens_i_d[ n2s_gens_idx[nth_idx] ] * cos(gens_δ[ n2s_gens_idx[nth_idx] ] )) -
                (gens_i_q[ n2s_gens_idx[nth_idx] ] * sin(gens_δ[ n2s_gens_idx[nth_idx] ] )) ) +
                sum( [ vh[idx] * abs(ynj) * sin( θh[idx] + angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ nth_idx],
                                nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx])]) )
            for nth_idx in gens_nodes_idx]        
        
    else

        gens_I_imag_mismatch =
            [ ( -( (gens_i_d[ n2s_gens_idx[nth_idx] ] * cos(gens_δ[ n2s_gens_idx[nth_idx] ] )) -
            (gens_i_q[ n2s_gens_idx[nth_idx] ] *  sin(gens_δ[ n2s_gens_idx[nth_idx] ] )) ) -
            sum( [ vh[idx] * abs(ynj) * sin( θh[idx] + angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ nth_idx],
                            nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx])]) )
              for nth_idx in gens_nodes_idx ]
        
    end
    
    # -----------------------------------
    # non gens nodes  current real part mismatch        
    # -----------------------------------

    """

     ( Pl_i * cos(θ_i) + Ql_i * sin(θ_i) ) / vh_i =
        ∑( vh_k * Y_ik * cos(θ_k + β_ik ) )


    """

    non_gens_I_real_mismatch = [
            ( ( P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] * cos( θh[nth_idx] ) +
                Q_non_gens[ n2s_non_gens_idx[ nth_idx ]]  * sin(θh[nth_idx]))/vh[nth_idx] +
            sum([ vh[idx] * abs(ynj) * cos(θh[idx] + angle(ynj))
                   for (ynj, idx) in
                       zip( Ynet[ nth_idx ],
                            nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ])]) )
        for nth_idx in non_gens_nodes_idx ]
    
    # ------------------------------------
    # non gens nodes  current imag part mismatch                
    # ------------------------------------

    """

    ( Pl_i * sin(θ_i) - Ql_i * cos(θ_i) ) / vh_i =
        -∑( vh_k * Y_ik * sin(θ_k + β_ik ) )

    """

    non_gens_I_imag_mismatch = [  ( P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] * sin( θh[ nth_idx ] ) -
                Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ] * cos( θh[ nth_idx ]))/vh[ nth_idx ] +
                sum( [ vh[idx] * abs(ynj) * sin( θh[idx] + angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ nth_idx],
                                nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx])]) 
            for nth_idx in non_gens_nodes_idx ]  
    
    # ------------------------------------
        
    ΔI_id_iq .=
        vcat(gens_I_real_mismatch,
             non_gens_I_real_mismatch,
             gens_I_imag_mismatch,
             non_gens_I_imag_mismatch,
             
             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag)

    return nothing

end




function get_a_model_integrated_init_dyn_pf_ΔPQ_mismatch(
    ΔPQ, vh_θh,
    init_dyn_pf_flat_para;
    init_dyn_pf_mismatch_kwd_para =
        init_dyn_pf_mismatch_kwd_para  )

    
    #-------------------------------

    (;
     loc_load_exist,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,     
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     non_pre_ordered_dyn_pf_vars_Idxs) =
         init_dyn_pf_mismatch_kwd_para

    (;
     gens_vh_idxs,
     gens_θh_idxs,

     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,

     gens_id_idxs,
     gens_iq_idxs,
     flat_vh_flat_θh_flat_id_iq_Idx ) =
         non_pre_ordered_dyn_pf_vars_Idxs
    
    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #-------------------------------
    
    # (;
    #  loc_load_exist,
    #  dyn_pf_vars_Idxs,
    #  dyn_pf_Idxs_kwd_para,
    #  gens_nodes_ra_Xd_dash_Xq_dash) =
    #      intra_dyn_pf_mismatch_kwd_para 

    #-------------------------------
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        dyn_pf_vars_Idxs

    #-------------------------------    
    
    (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
    
    #-------------------------------

   (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) =
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs 
    
    #-------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #-------------------------------

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------
    
    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
        
    #-------------------------------

    # gens_vh =
    #     vh_θh[ gens_vh_Idxs ]

    # gens_θh =
    #     vh_θh[ gens_θh_Idxs ]

    # non_gens_nodes_vh =
    #     vh_θh[ non_gens_nodes_vh_Idxs ]

    # non_gens_nodes_θh =
    #     vh_θh[ non_gens_nodes_θh_Idxs ]


    # gens_vh =
    #     vh_θh[ gens_vh_idxs ]

    # gens_θh =
    #     vh_θh[ gens_θh_idxs ]

    # non_gens_nodes_vh =
    #     vh_θh[ non_gens_nodes_vh_idxs ]

    # non_gens_nodes_θh =
    #     vh_θh[ non_gens_nodes_θh_idxs ]

    # #-------------------------------
    
    flat_vh_Idx, flat_θh_Idx, flat_id_Idx, flat_iq_Idx =
        flat_vh_flat_θh_flat_id_iq_Idx
    
    flat_vh =
        vh_θh[ flat_vh_Idx ]

    flat_θh =
        vh_θh[ flat_θh_Idx ]    

    gens_vh =
        flat_vh[ gens_nodes_idx ]

    gens_θh =
        flat_θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[ non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[ non_gens_nodes_idx ]
    
    #-------------------------------
    
    
    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ gens_nodes_idx ?
            gens_θh[
                n2s_gens_idx[ idx] ] :
                    non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    # -------------------------------------
    
    # -----------------------------------
    # gens active power mismatch
    # -----------------------------------
    
    """
    P_gens + Ploc_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true

        gens_P_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( P_gens[n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ] -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) ) :
                                 ( 
                                 P_gens[n2s_gens_idx[nth_idx]]  -
                                 vh[nth_idx] * sum( [
                                     vh[idx] *
                                         abs(ynj) *
                                         cos(θh[nth_idx] - θh[idx] - angle(ynj))
                                     for (ynj, idx) in
                                         zip( Ynet[ nth_idx ],
                                              nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ])]) )

            for nth_idx in gens_nodes_idx ]  

    else

        gens_P_mismatch = [ P_gens[n2s_gens_idx[nth_idx]]  -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) for nth_idx in
                gens_nodes_idx ]  
        
        
    end

    # -----------------------------------
    # gens reactive power mismatch
    # -----------------------------------

    """
    Q_gens + Qloc_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true
        
        gens_Q_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( Q_gens[n2s_gens_idx[nth_idx]] -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx] ] - 
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj) )
                for (ynj, idx) in
                    zip( Ynet[ nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])]) ) :
                                 ( Q_gens[n2s_gens_idx[nth_idx]] -
                                         vh[nth_idx] *
                                         sum([vh[idx] * abs(ynj) *
                                         sin(θh[nth_idx]-θh[idx]-angle(ynj))
                                              for (ynj, idx) in
                                                  zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])]) ) 
                            for nth_idx in
                                gens_nodes_idx]        
        
    else

        gens_Q_mismatch = [ Q_gens[n2s_gens_idx[nth_idx]] -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])])
                            for nth_idx in
                                gens_nodes_idx ]
        
    end
    
    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------

    """
    Pl_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )
    """
        
    non_gens_P_mismatch = [
        P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------

    """
    Ql_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    non_gens_Q_mismatch = [
        Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------

    ΔPQ .=
        vcat(gens_P_mismatch,
             non_gens_P_mismatch,
             gens_Q_mismatch,                 
             non_gens_Q_mismatch)

    return nothing

end


function get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch(
    ΔPQ_id_iq, vh_θh_id_iq, intra_dyn_pf_mismatch_flat_para;
    intra_dyn_pf_mismatch_kwd_para = intra_dyn_pf_mismatch_kwd_para  )
    
    #-------------------------------

    (;
     loc_load_exist,
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     non_pre_ordered_dyn_pf_vars_Idxs) =
         intra_dyn_pf_mismatch_kwd_para 

    #-------------------------------


    (;
     gens_vh_idxs,
     gens_θh_idxs,

     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,

     gens_id_idxs,
     gens_iq_idxs,
     
     flat_vh_flat_θh_flat_id_iq_Idx) =
         non_pre_ordered_dyn_pf_vars_Idxs
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        dyn_pf_vars_Idxs

    #-------------------------------    
    
    (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
    
    #-------------------------------

   (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) =
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs 
    
    #-------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #-------------------------------

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------
    
    P_gens =
        intra_dyn_pf_mismatch_flat_para[
            dyn_P_gens_Idxs ]

    Q_gens =
        intra_dyn_pf_mismatch_flat_para[
            dyn_Q_gens_Idxs ]

    P_non_gens  =
        intra_dyn_pf_mismatch_flat_para[
            dyn_P_non_gens_Idxs ]

    Q_non_gens = 
        intra_dyn_pf_mismatch_flat_para[
            dyn_Q_non_gens_Idxs ]

    flat_δ_ed_dash_eq_dash =
        intra_dyn_pf_mismatch_flat_para[
            dyn_δ_ed_eq_pf_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            intra_dyn_pf_mismatch_flat_para[
                dyn_P_g_loc_load_Idxs ]

        Q_g_loc_load =
            intra_dyn_pf_mismatch_flat_para[
                dyn_Q_g_loc_load_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #-------------------------------

    if length(flat_δ_ed_dash_eq_dash) !=
        length( δ_ed_dash_eq_dash_Idxs_in_flattend )

        δ_ed_dash_eq_dash =
            [ flat_δ_ed_dash_eq_dash[idx]
              for idx in
                  δ_ed_dash_eq_dash_Idxs_in_flattend ]

    else

        δ_ed_dash_eq_dash =
            flat_δ_ed_dash_eq_dash

    end
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
        
    #-------------------------------

    # gens_vh =
    #     vh_θh_id_iq[ gens_vh_Idxs ]

    # gens_θh =
    #     vh_θh_id_iq[ gens_θh_Idxs ]

    # non_gens_nodes_vh =
    #     vh_θh_id_iq[ non_gens_nodes_vh_Idxs ]

    # non_gens_nodes_θh =
    #     vh_θh_id_iq[ non_gens_nodes_θh_Idxs ]

    # gens_i_d =
    #     vh_θh_id_iq[ gen_id_Idxs ]

    # gens_i_q =
    #     vh_θh_id_iq[ gen_iq_Idxs ]


    # gens_vh =
    #     vh_θh_id_iq[ gens_vh_idxs ]

    # gens_θh =
    #     vh_θh_id_iq[ gens_θh_idxs ]

    # non_gens_nodes_vh =
    #     vh_θh_id_iq[ non_gens_nodes_vh_idxs ]

    # non_gens_nodes_θh =
    #     vh_θh_id_iq[ non_gens_nodes_θh_idxs ]

    # gens_i_d =
    #     vh_θh_id_iq[ gens_id_idxs ]

    # gens_i_q =
    #     vh_θh_id_iq[ gens_iq_idxs ]

    # #-------------------------------
    
    flat_vh_Idx, flat_θh_Idx, flat_id_Idx, flat_iq_Idx =
        flat_vh_flat_θh_flat_id_iq_Idx
    

    flat_vh =
        vh_θh_id_iq[ flat_vh_Idx ]

    flat_θh =
        vh_θh_id_iq[ flat_θh_Idx ]

    gens_i_d =
        vh_θh_id_iq[ flat_id_Idx ]

    gens_i_q =
        vh_θh_id_iq[ flat_iq_Idx ]
    

    gens_vh =
        flat_vh[ gens_nodes_idx ]

    gens_θh =
        flat_θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[ non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[ non_gens_nodes_idx ]
        
    #-------------------------------
    
    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ gens_nodes_idx ?
            gens_θh[
                n2s_gens_idx[ idx] ] :
                    non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    # -------------------------------------

    gens_δ = first.( δ_ed_dash_eq_dash )

    gens_ed_dash = second.( δ_ed_dash_eq_dash )

    gens_eq_dash = third.( δ_ed_dash_eq_dash )

    # -------------------------------------

    gens_ra = first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash = second.( gens_nodes_ra_Xd_dash_Xq_dash  )

    gens_Xq_dash = third.( gens_nodes_ra_Xd_dash_Xq_dash  )        

    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------

    """ ed_dash - vh * sin(δ - θ) - ra * id + Xq_dash * iq """
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            vh[idx] * sin(gens_δ[n2s_gens_idx[idx]] - θh[idx]) -
            gens_ra[n2s_gens_idx[idx]] * gens_i_d[n2s_gens_idx[idx]] +
            gens_Xq_dash[ n2s_gens_idx[idx]] * gens_i_q[ n2s_gens_idx[idx]]
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θ) - ra * iq - Xd_dash * id """
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] -
            vh[ idx] * cos(gens_δ[ n2s_gens_idx[idx]] - θh[idx]) -
            gens_ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[idx]] -
            gens_Xd_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]]
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    # gens active power mismatch
    # -----------------------------------
    
    """
    id_i * vh_i * sin(δ_i - θ_i) + iq_i * vh_i * cos(δ_i - θ_i) + Ploc_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true

        gens_P_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ] -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) ) :
                                 ( (gens_i_d[ n2s_gens_idx[ nth_idx] ] *
                                 vh[nth_idx] *
                                 sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
                                 (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
                                 cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
                                 vh[nth_idx] * sum( [
                                     vh[idx] *
                                         abs(ynj) *
                                         cos(θh[nth_idx] - θh[idx] - angle(ynj))
                                     for (ynj, idx) in
                                         zip( Ynet[ nth_idx ],
                                              nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ])]) )

            for nth_idx in gens_nodes_idx ]  

    else

        gens_P_mismatch = [ (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) for nth_idx in
                gens_nodes_idx ]  
        
        
    end

    # -----------------------------------
    # gens reactive power mismatch
    # -----------------------------------

    """
    id_i * vh_i * cos(δ_i - θ_i) - iq_i * vh_i * sin(δ_i - θ_i) + Qloc_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true
        
        gens_Q_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ nth_idx ] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx] ] - 
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj) )
                for (ynj, idx) in
                    zip( Ynet[ nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])]) ) :
                                 ((gens_i_d[ n2s_gens_idx[ nth_idx]] *
                                 vh[nth_idx] *
                                 cos(gens_δ[ n2s_gens_idx[ nth_idx]] - θh[nth_idx])) -
                                 (gens_i_q[ n2s_gens_idx[ nth_idx]] * vh[ nth_idx] *
                                 sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[nth_idx])) -
                                         vh[nth_idx] *
                                         sum([vh[idx] * abs(ynj) *
                                         sin(θh[nth_idx]-θh[idx]-angle(ynj))
                                              for (ynj, idx) in
                                                  zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])]) ) 
                            for nth_idx in
                                gens_nodes_idx]        
        
    else

        gens_Q_mismatch = [ (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ nth_idx ] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])])
                            for nth_idx in
                                gens_nodes_idx ]
        
    end
    
    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------

    """
    Pl_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )
    """
        
    non_gens_P_mismatch = [
        P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------

    """
    Ql_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    non_gens_Q_mismatch = [
        Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------

    ΔPQ_id_iq .=
        vcat(gens_P_mismatch,
             non_gens_P_mismatch,
             gens_Q_mismatch,                 
             non_gens_Q_mismatch,
             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag)

    return nothing

end




function get_a_model_integrated_intra_dyn_pf_ΔPQ_mismatch(
    ΔPQ, flat_vh_flat_θh, intra_dyn_pf_mismatch_flat_para;
    intra_dyn_pf_mismatch_kwd_para =
        intra_dyn_pf_mismatch_kwd_para  )
    
    #-------------------------------

    (;
     loc_load_exist,
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     non_pre_ordered_dyn_pf_vars_Idxs) =
         intra_dyn_pf_mismatch_kwd_para 

    #-------------------------------


    (;
     gens_vh_idxs,
     gens_θh_idxs,

     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,

     gens_id_idxs,
     gens_iq_idxs,
     
     flat_vh_flat_θh_flat_id_iq_Idx) =
         non_pre_ordered_dyn_pf_vars_Idxs
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        dyn_pf_vars_Idxs

    #-------------------------------    
    
    (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
    
    #-------------------------------

   (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) =
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs 
    
    #-------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #-------------------------------

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------
    
    P_gens =
        intra_dyn_pf_mismatch_flat_para[
            dyn_P_gens_Idxs ]

    Q_gens =
        intra_dyn_pf_mismatch_flat_para[
            dyn_Q_gens_Idxs ]

    P_non_gens  =
        intra_dyn_pf_mismatch_flat_para[
            dyn_P_non_gens_Idxs ]

    Q_non_gens = 
        intra_dyn_pf_mismatch_flat_para[
            dyn_Q_non_gens_Idxs ]

    flat_δ_ed_dash_eq_dash =
        intra_dyn_pf_mismatch_flat_para[
            dyn_δ_ed_eq_pf_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            intra_dyn_pf_mismatch_flat_para[
                dyn_P_g_loc_load_Idxs ]

        Q_g_loc_load =
            intra_dyn_pf_mismatch_flat_para[
                dyn_Q_g_loc_load_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #-------------------------------

    if length(flat_δ_ed_dash_eq_dash) !=
        length( δ_ed_dash_eq_dash_Idxs_in_flattend )

        δ_ed_dash_eq_dash =
            [ flat_δ_ed_dash_eq_dash[idx]
              for idx in
                  δ_ed_dash_eq_dash_Idxs_in_flattend ]

    else

        δ_ed_dash_eq_dash =
            flat_δ_ed_dash_eq_dash

    end
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
        
    #-------------------------------
    #-------------------------------
    
    flat_vh_Idx, flat_θh_Idx, flat_id_Idx, flat_iq_Idx =
        flat_vh_flat_θh_flat_id_iq_Idx
    
    flat_vh =
        flat_vh_flat_θh[ flat_vh_Idx ]

    flat_θh =
        flat_vh_flat_θh[ flat_θh_Idx ]

    #-------------------------------    

    gens_vh =
        flat_vh[ gens_nodes_idx ]

    gens_θh =
        flat_θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[ non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[ non_gens_nodes_idx ]
        
    #-------------------------------
    
    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ gens_nodes_idx ?
            gens_θh[
                n2s_gens_idx[ idx] ] :
                    non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    # -------------------------------------

    gens_δ = first.( δ_ed_dash_eq_dash )

    gens_ed_dash = second.( δ_ed_dash_eq_dash )

    gens_eq_dash = third.( δ_ed_dash_eq_dash )

    # -------------------------------------

    gens_ra = first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash = second.(
        gens_nodes_ra_Xd_dash_Xq_dash  )

    gens_Xq_dash = third.(
        gens_nodes_ra_Xd_dash_Xq_dash )

    # -------------------------------------
    
    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
        for (a_vh, a_θh,
             a_δ, a_ed_dash, a_eq_dash,
             a_ra, a_X_d_dash, a_X_q_dash ) in
            zip( gens_vh, gens_θh,
                 gens_δ, gens_ed_dash, gens_eq_dash,
                 gens_ra, gens_Xd_dash,  gens_Xq_dash ) ]

    gens_i_d =
        first.( gens_id_iq )
    
    gens_i_q =
        last.( gens_id_iq )

    # # -----------------------------------
    # # gens real part stator equation mismatch
    # # -----------------------------------

    # """ ed_dash - vh * sin(δ - θ) - ra * id  + Xq_dash * iq """
    
    # gens_stator_equations_mismatch_real = [
    #     gens_ed_dash[ n2s_gens_idx[ idx ]] -
    #         vh[idx] * sin(gens_δ[n2s_gens_idx[idx]] - θh[idx]) -
    #         gens_ra[n2s_gens_idx[idx]] * gens_i_d[n2s_gens_idx[idx]] +
    #         gens_Xq_dash[ n2s_gens_idx[idx]] * gens_i_q[ n2s_gens_idx[idx]]
    #     for idx in gens_nodes_idx ]  
    
    # # -----------------------------------
    # # gens imag part stator equation mismatch
    # # -----------------------------------

    # """ eq_dash - vh * cos(δ - θ) - ra * iq  - Xd_dash * id """
    
    # gens_stator_equations_mismatch_imag = [
    #     gens_eq_dash[ n2s_gens_idx[ idx ]] -
    #         vh[ idx] * cos(gens_δ[ n2s_gens_idx[idx]] - θh[idx]) -
    #         gens_ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[idx]] -
    #         gens_Xd_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]]
    #     for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    # gens active power mismatch
    # -----------------------------------
    
    """
    id_i * vh_i * sin(δ_i - θ_i) + iq_i * vh_i * cos(δ_i - θ_i) + Ploc_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true

        gens_P_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx] ] -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) ) :
                                 ( (gens_i_d[ n2s_gens_idx[ nth_idx] ] *
                                 vh[nth_idx] *
                                 sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
                                 (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
                                 cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
                                 vh[nth_idx] * sum( [
                                     vh[idx] *
                                         abs(ynj) *
                                         cos(θh[nth_idx] - θh[idx] - angle(ynj))
                                     for (ynj, idx) in
                                         zip( Ynet[ nth_idx ],
                                              nodes_node_idx_and_incident_edges_other_node_idx[ nth_idx ])]) )

            for nth_idx in gens_nodes_idx ]  

    else

        gens_P_mismatch = [ (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[nth_idx])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[nth_idx]))  -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])]) for nth_idx in
                gens_nodes_idx ]  
        
        
    end

    # -----------------------------------
    # gens reactive power mismatch
    # -----------------------------------

    """
    id_i * vh_i * cos(δ_i - θ_i) - iq_i * vh_i * sin(δ_i - θ_i) + Qloc_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true
        
        gens_Q_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ nth_idx ] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx] ] - 
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj) )
                for (ynj, idx) in
                    zip(Ynet[nth_idx],
                                      nodes_node_idx_and_incident_edges_other_node_idx[nth_idx])])) : ((gens_i_d[ n2s_gens_idx[ nth_idx]] *
                                 vh[nth_idx] *
                                 cos(gens_δ[ n2s_gens_idx[ nth_idx]] - θh[nth_idx])) -
                                 (gens_i_q[ n2s_gens_idx[ nth_idx]] * vh[ nth_idx] *
                                 sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[nth_idx])) -
                                         vh[nth_idx] *
                                         sum([vh[idx] * abs(ynj) *
                                         sin(θh[nth_idx]-θh[idx]-angle(ynj))
                                              for (ynj, idx) in
                                                  zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])]) ) 
                            for nth_idx in
                                gens_nodes_idx]        
        
    else

        gens_Q_mismatch = [ (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[nth_idx] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ nth_idx ] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ nth_idx ])) -
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin( θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[nth_idx],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx])])
                            for nth_idx in
                                gens_nodes_idx ]
        
    end
    
    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------

    """
    Pl_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )
    """
        
    non_gens_P_mismatch = [
        P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    cos(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------

    """
    Ql_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    non_gens_Q_mismatch = [
        Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[nth_idx] * sum( [
                vh[idx] *
                    abs(ynj) *
                    sin(θh[nth_idx] - θh[idx] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ nth_idx ],
                         nodes_node_idx_and_incident_edges_other_node_idx[
                             nth_idx ])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------

    ΔPQ .=
        vcat(gens_P_mismatch,
             non_gens_P_mismatch,
             gens_Q_mismatch,                 
             non_gens_Q_mismatch
             )

    return nothing

end

#---------------------------------------------------------
##########################################################
#---------------------------------------------------------

function get_power_balance_mismatch_generic_with_vh_θh!(
    red_ΔPQ,
    red_vh_θh_0_view,
    (working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view),
    pf_net_param;
    with_δ_ed_eq =
        true )

    # ----------------------------------------------
    # ----------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_idx_with_adjacent_nodes_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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
    
    # ----------------------------------------------
    # ----------------------------------------------

    slack_gens_θh =
        zeros(length( slack_bus_idx))
    
    # ----------------------------------------------

    all_nodes_idx = first.(
        nodes_idx_with_adjacent_nodes_idx)
    
    non_gen_nodes_idx =
        setdiff(all_nodes_idx, gens_idx)
        
    non_slack_gen_nodes_idx =
        setdiff(gens_idx, slack_bus_idx)

    non_slack_gen_and_non_gen_nodes_idx =
        sort([non_slack_gen_nodes_idx;
              non_gen_nodes_idx ])

    _,_, θh_idx_non_slack_gen_and_non_gen_nodes =
        create_size_offset_Idx(
            length.(non_slack_gen_nodes_idx,
                    non_gen_nodes_idx))

    θh_idx_non_slack_gen =
        θh_idx_non_slack_gen_and_non_gen_nodes[1]
    
    θh_idx_non_gen_nodes =
        θh_idx_non_slack_gen_and_non_gen_nodes[2]
    
    _, _, vh_θh_idx_red_vh_θh =
        create_size_offset_Idx(
            length.(non_gen_nodes_idx,
                    non_slack_gen_and_non_gen_nodes_idx))

    non_gen_vh_idx =
        vh_θh_idx_red_vh_θh[1]
    
    non_slack_and_non_gen_θh_idx =
        vh_θh_idx_red_vh_θh[2]
    
    # ----------------------------------------------

    n2s_slack_gens_idx =
        get_a_n2s_dict(slack_bus_idx )
            
    n2s_gens_idx =
        get_a_n2s_dict(gens_idx )
           
    n2s_non_gen_nodes_idx =
        get_a_n2s_dict( non_gen_nodes_idx )
    
    n2s_all_nodes_idx =
        get_a_n2s_dict( all_nodes_idx )
    
    # ----------------------------------------------

    slack_gens_nodes_idx = slack_bus_idx
    
    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_bus_idx ]
    
    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_idx ]

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    transformed_non_gen_nodes_idx = [
        n2s_all_nodes_idx[idx]
        
        for idx in non_gen_nodes_idx ]
    # ----------------------------------------------

    non_gens_vh  =
        red_vh_θh_0_view[ non_gen_vh_idx]

    non_slack_and_non_gen_θh =
        red_vh_θh_0_view[non_slack_and_non_gen_θh_idx]
    
    non_slack_gens_θh =
        non_slack_and_non_gen_θh[θh_idx_non_slack_gen]
    
    non_gens_θh =
        non_slack_and_non_gen_θh[θh_idx_non_gen_nodes]


    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[idx]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]
    
    
    # ----------------------------------------------

    working_vh_θh_view[ red_vh_θh_idx ] .=
        red_vh_θh_0_view

    uh =  working_vh_θh_view[ vh_IDX ] .* exp.(
        im * working_vh_θh_view[ θh_IDX ])
        
    # ----------------------------------------------

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    # ----------------------------------------------

    if with_δ_ed_eq == true
        
        idq_θ_π_vhθh = [
            get_dynamic_idq_θ_π_vhθh(
                vh, θh,
                δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
            for (vh,θh,δ_ω_ed_eq,ra_Xd_dash_Xq_dash) in
                zip(abs.(uh),
                    angle.(uh),
                    δ_ω_ed_dash_eq_dash_view,
                    ra_Xd_dash_Xq_dash_view)]

        I_sum_ynj_vj =
            get_nodes_∑_ynj_x_vj(
                uh,
                Ynet,
                nodes_idx_with_adjacent_nodes_idx )

        current_mismatch =
            get_nodes_current_mismatch(
                uh,
                (P_Q_non_gens_view,
                 P_Q_gens_loc_load_view),
                idq_θ_π_vhθh,
                I_sum_ynj_vj )
        
    else

        I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
            uh,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
        current_mismatch =
            get_nodes_current_mismatch(
                uh,
                (P_Q_gens_view,
                 P_Q_non_gens_view,
                 P_Q_gens_loc_load_view),
                I_sum_ynj_vj )

    end
        
    # ----------------------------------------------
    
    power_mismatch =
        uh .* conj.(current_mismatch)

    red_P_mismatch =
        real.(power_mismatch[
            setdiff(1:end, slack_bus_idx)])
    
    red_Q_mismatch =
        imag.(power_mismatch[
            setdiff(1:end,gens_idx)])
        
    red_ΔPQ .=
        vcat(red_P_mismatch,
             red_Q_mismatch)

    # ------------------------------------------------------
    # update Inet_view  and Iinj_view
    # -----------------------------------------------------

    Inet_view .= I_sum_ynj_vj

    Iinj_view .= get_Iinj(
        uh,
        P_Q_non_gens_view,
        P_Q_gens_loc_load_view,
        Inet_view )
    
    # ----------------------------------------------------   
    # update Q 
    # ----------------------------------------------------


    if with_δ_ed_eq == true
        
        Igen = idq_θ_π_vhθh 
        
        gens_S  = uh[gens_idx] .*
            conj.(Igen[gens_idx])
        
        P_Q_nodes_view[slack_bus_idx] .= [
            real(gens_S[slack_bus_idx]),
            imag(gens_S[slack_bus_idx])]

       for (idx, gen_S) in zip(gens_idx, gens_S)

           P_Q_gens_view[idx] .= [
               real( S_gens[idx]), imag(gen_S)]
       end
    else
  
        Igen = I_sum_ynj_vj +
            get_gens_local_load_current(
                uh, P_Q_gens_loc_load_view) 

        gens_S  = uh[gens_idx] .* conj.(Igen[gens_idx])

        P_Q_nodes_view[slack_bus_idx] .=
            [real(gens_S[slack_bus_idx]),
             imag(gens_S[slack_bus_idx])]

        for (idx, gen_S) in zip(gens_idx, gens_S)

            P_Q_gens_view[idx] .= [
                real( S_gens[idx]), imag(gen_S)]
        end
    end
    
    # ------------------------------------------------------ 
    # update nodes_pf_U_view
    # ------------------------------------------------------ 

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .= [
            real( uk ), imag( uk ) ]

    end

    return nothing
    
end



#-----------------------------------------------------
# mismatch models with idx that could be a union of
# string, symbol or integer
#-----------------------------------------------------


function get_a_init_dyn_pf_model_ΔPQ_mismatch(
    ΔPQ, vh_θh,
    init_dyn_pf_flat_para;
    init_dyn_pf_mismatch_kwd_para =
        init_dyn_pf_mismatch_kwd_para  )
    
    #-------------------------------

    (;
     loc_load_exist,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,     
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     non_pre_ordered_dyn_pf_vars_Idxs) =
         init_dyn_pf_mismatch_kwd_para

    (;
     gens_vh_idxs,
     gens_θh_idxs,

     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,

     gens_id_idxs,
     gens_iq_idxs,
     flat_vh_flat_θh_flat_id_iq_Idx ) =
         non_pre_ordered_dyn_pf_vars_Idxs
    
    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #-------------------------------
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        dyn_pf_vars_Idxs

    #-------------------------------    
    
    (;dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
    
    #-------------------------------

   (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) =
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs 
    
    #-------------------------------

   (;n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------

   (;slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #-------------------------------

   (Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------
    
    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
        
    #-------------------------------
    
    (flat_vh_Idx,
     flat_θh_Idx,
     flat_id_Idx,
     flat_iq_Idx) =
        flat_vh_flat_θh_flat_id_iq_Idx
    
    flat_vh =
        vh_θh[ flat_vh_Idx ]

    flat_θh =
        vh_θh[ flat_θh_Idx ]    

    gens_vh =
        flat_vh[ gens_nodes_idx ]

    gens_θh =
        flat_θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[ non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[ non_gens_nodes_idx ]
    
    #-------------------------------
        
    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ gens_nodes_idx ?
            gens_θh[
                n2s_gens_idx[ idx] ] :
                    non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]
    
    # -----------------------------------
    # gens active power mismatch
    # -----------------------------------
    
    """
    P_gens + Ploc_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true

        gens_P_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                (P_gens[ n2s_gens_idx[nth_idx]]  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [
                vh[ n2s_all_nodes_idx[idx]] *
                    abs(ynj) *
                    cos(θh[ n2s_all_nodes_idx[ nth_idx]] -
                    θh[ n2s_all_nodes_idx[idx]] -
                    angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                         nodes_idx_with_adjacent_nodes_idx[
                             n2s_all_nodes_idx[nth_idx]])])) :
                                 (P_gens[n2s_gens_idx[nth_idx]] -
                                 vh[n2s_all_nodes_idx[nth_idx]] * sum([
                                     vh[ n2s_all_nodes_idx[idx]] *
                                         abs(ynj) *
                                         cos(θh[ n2s_all_nodes_idx[nth_idx]] -
                                         θh[ n2s_all_nodes_idx[idx]] -
                                         angle(ynj))
                                     for (ynj, idx) in
                                         zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                             nodes_idx_with_adjacent_nodes_idx[
                                                 nth_idx ])]) )

            for nth_idx in gens_nodes_idx ]  

    else

        gens_P_mismatch = [
            P_gens[n2s_gens_idx[nth_idx]]  -
                vh[ n2s_all_nodes_idx[nth_idx]] * sum([
                    vh[ n2s_all_nodes_idx[idx]] *
                        abs(ynj) *
                        cos(θh[ n2s_all_nodes_idx[nth_idx]] -
                        θh[ n2s_all_nodes_idx[idx]] -
                        angle(ynj))
                    for (ynj, idx) in
                        zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])])
            for nth_idx in gens_nodes_idx ]  
        
        
    end

    # -----------------------------------
    # gens reactive power mismatch
    # -----------------------------------

    """
    Q_gens + Qloc_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true
        
        gens_Q_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                (Q_gens[n2s_gens_idx[nth_idx]] -
                Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx] ] - 
                vh[ n2s_all_nodes_idx[nth_idx]] * sum( [
                    vh[n2s_all_nodes_idx[idx]] *
                       abs(ynj) *
                       sin( θh[ n2s_all_nodes_idx[nth_idx]] -
                       θh[ n2s_all_nodes_idx[idx]] -
                       angle(ynj) )
                    for (ynj, idx) in
                        zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])])) :
                                    ( Q_gens[n2s_gens_idx[nth_idx]] -
                                    vh[ n2s_all_nodes_idx[nth_idx]] *
                                    sum([vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                                    sin(θh[ n2s_all_nodes_idx[nth_idx]] -
                                    θh[ n2s_all_nodes_idx[idx]] -
                                    angle(ynj))
                                         for (ynj, idx) in
                                             zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                                 nodes_idx_with_adjacent_nodes_idx[
                                                      n2s_all_nodes_idx[nth_idx]])]) ) 
            for nth_idx in gens_nodes_idx]        
        
    else

        gens_Q_mismatch = [
            Q_gens[n2s_gens_idx[nth_idx]] -
                vh[ n2s_all_nodes_idx[nth_idx]] * sum( [
                    vh[ n2s_all_nodes_idx[idx]] *
                        abs(ynj) *
                        sin( θh[ n2s_all_nodes_idx[nth_idx]] -
                        θh[ n2s_all_nodes_idx[idx]] -
                        angle(ynj))
                    for (ynj, idx) in
                        zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])])

            for nth_idx in gens_nodes_idx ]
        
    end
    
    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------

    """
    Pl_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )
    """
        
    non_gens_P_mismatch = [
        P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [
                vh[ n2s_all_nodes_idx[idx]] *
                    abs(ynj) *
                    cos(θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[idx]] -
                    angle(ynj))
                for (ynj, idx) in
                    zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                        nodes_idx_with_adjacent_nodes_idx[
                             n2s_all_nodes_idx[nth_idx]])])
        for nth_idx in non_gens_nodes_idx ]
    
    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------

    """
    Ql_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    non_gens_Q_mismatch = [
        Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [
                vh[n2s_all_nodes_idx[idx]] *
                   abs(ynj) *
                   sin(θh[ n2s_all_nodes_idx[nth_idx]] -
                   θh[ n2s_all_nodes_idx[idx]] -
                   angle(ynj))
                for (ynj, idx) in
                    zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                        nodes_idx_with_adjacent_nodes_idx[
                             n2s_all_nodes_idx[nth_idx]])])
        for nth_idx in non_gens_nodes_idx ]
    
    # ------------------------------------

    ΔPQ .=
        vcat(gens_P_mismatch,
             non_gens_P_mismatch,
             gens_Q_mismatch,                 
             non_gens_Q_mismatch)

    return nothing

end




function get_a_intra_dyn_pf_model_ΔPQ_Δidq_mismatch(
    ΔPQ_id_iq, vh_θh_id_iq, intra_dyn_pf_mismatch_flat_para;
    intra_dyn_pf_mismatch_kwd_para = intra_dyn_pf_mismatch_kwd_para  )
    
    #-------------------------------

    (;
     loc_load_exist,
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     non_pre_ordered_dyn_pf_vars_Idxs) =
         intra_dyn_pf_mismatch_kwd_para 

    #-------------------------------


    (;
     gens_vh_idxs,
     gens_θh_idxs,

     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,

     gens_id_idxs,
     gens_iq_idxs,
     
     flat_vh_flat_θh_flat_id_iq_Idx) =
         non_pre_ordered_dyn_pf_vars_Idxs
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        dyn_pf_vars_Idxs

    #-------------------------------    
    
    (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
    
    #-------------------------------

   (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) =
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs 
    
    #-------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #-------------------------------

   ( Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        dyn_pf_fun_kwd_net_para
    
    #-------------------------------
    
    P_gens =
        intra_dyn_pf_mismatch_flat_para[
            dyn_P_gens_Idxs ]

    Q_gens =
        intra_dyn_pf_mismatch_flat_para[
            dyn_Q_gens_Idxs ]

    P_non_gens  =
        intra_dyn_pf_mismatch_flat_para[
            dyn_P_non_gens_Idxs ]

    Q_non_gens = 
        intra_dyn_pf_mismatch_flat_para[
            dyn_Q_non_gens_Idxs ]

    flat_δ_ed_dash_eq_dash =
        intra_dyn_pf_mismatch_flat_para[
            dyn_δ_ed_eq_pf_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            intra_dyn_pf_mismatch_flat_para[
                dyn_P_g_loc_load_Idxs ]

        Q_g_loc_load =
            intra_dyn_pf_mismatch_flat_para[
                dyn_Q_g_loc_load_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #-------------------------------

    if length(flat_δ_ed_dash_eq_dash) !=
        length( δ_ed_dash_eq_dash_Idxs_in_flattend )

        δ_ed_dash_eq_dash =
            [ flat_δ_ed_dash_eq_dash[idx]
              for idx in
                  δ_ed_dash_eq_dash_Idxs_in_flattend ]

    else

        δ_ed_dash_eq_dash =
            flat_δ_ed_dash_eq_dash

    end
    
    #-------------------------------

    nodes_size = length( gens_nodes_idx ) +
        length( non_gens_nodes_idx )
        
    #-------------------------------
    
    flat_vh_Idx, flat_θh_Idx, flat_id_Idx, flat_iq_Idx =
        flat_vh_flat_θh_flat_id_iq_Idx
    

    flat_vh =
        vh_θh_id_iq[ flat_vh_Idx ]

    flat_θh =
        vh_θh_id_iq[ flat_θh_Idx ]

    gens_i_d =
        vh_θh_id_iq[ flat_id_Idx ]

    gens_i_q =
        vh_θh_id_iq[ flat_iq_Idx ]
    

    gens_vh =
        flat_vh[ gens_nodes_idx ]

    gens_θh =
        flat_θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[ non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[ non_gens_nodes_idx ]
        
    #-------------------------------
    
    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ gens_nodes_idx ?
            gens_θh[
                n2s_gens_idx[ idx] ] :
                    non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    # -------------------------------------

    gens_δ = first.( δ_ed_dash_eq_dash )

    gens_ed_dash = second.( δ_ed_dash_eq_dash )

    gens_eq_dash = third.( δ_ed_dash_eq_dash )

    # -------------------------------------

    gens_ra = first.(gens_nodes_ra_Xd_dash_Xq_dash)

    gens_Xd_dash = second.(gens_nodes_ra_Xd_dash_Xq_dash)

    gens_Xq_dash = third.(gens_nodes_ra_Xd_dash_Xq_dash)   

    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------

    """ ed_dash - vh * sin(δ - θ) - ra * id + Xq_dash * iq """
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            vh[ n2s_all_nodes_idx[idx]] *
            sin(gens_δ[n2s_gens_idx[idx]] -
            θh[ n2s_all_nodes_idx[idx]]) -
            gens_ra[n2s_gens_idx[idx]] *
            gens_i_d[n2s_gens_idx[idx]] +
            gens_Xq_dash[ n2s_gens_idx[idx]] *
            gens_i_q[ n2s_gens_idx[idx]]
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θ) - ra * iq - Xd_dash * id """
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] -
            vh[ n2s_all_nodes_idx[idx]] *
            cos(gens_δ[ n2s_gens_idx[idx]] -
            θh[ n2s_all_nodes_idx[idx]]) -
            gens_ra[ n2s_gens_idx[ idx]] *
            gens_i_q[ n2s_gens_idx[idx]] -
            gens_Xd_dash[ n2s_gens_idx[idx]] *
            gens_i_d[ n2s_gens_idx[ idx]]
        for idx in gens_nodes_idx ]  

    
    # -----------------------------------
    # gens active power mismatch
    # -----------------------------------
    
    """
    id_i * vh_i * sin(δ_i - θ_i) + iq_i * vh_i * cos(δ_i - θ_i) + Ploc_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true

        gens_P_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ( (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[ n2s_all_nodes_idx[nth_idx]] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ n2s_all_nodes_idx[nth_idx]])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ n2s_all_nodes_idx[nth_idx]] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[ n2s_all_nodes_idx[nth_idx]]))  -
            P_g_loc_load[ n2s_gens_with_loc_load_idxs[ n2s_all_nodes_idx[nth_idx]] ] -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [
                vh[ n2s_all_nodes_idx[idx]] *
                    abs(ynj) *
                    cos(θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                         nodes_idx_with_adjacent_nodes_idx[
                             n2s_all_nodes_idx[nth_idx] ])]) ) :
                                 ( (gens_i_d[ n2s_gens_idx[ nth_idx] ] *
                                 vh[ n2s_all_nodes_idx[nth_idx]] *
                                 sin(gens_δ[ n2s_gens_idx[nth_idx] ] -
                                 θh[ n2s_all_nodes_idx[nth_idx]])) +
                                 (gens_i_q[ n2s_gens_idx[nth_idx] ] *
                                 vh[ n2s_all_nodes_idx[nth_idx]] *
                                 cos(gens_δ[ n2s_gens_idx[nth_idx] ]-
                                 θh[ n2s_all_nodes_idx[nth_idx]]))  -
                                 vh[nth_idx] * sum( [
                                     vh[idx] *
                                         abs(ynj) *
                                         cos(θh[ n2s_all_nodes_idx[nth_idx]] -
                                         θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
                                     for (ynj, idx) in
                                         zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                              nodes_idx_with_adjacent_nodes_idx[
                                                  nth_idx ])]) )

            for nth_idx in gens_nodes_idx ]  

    else

        gens_P_mismatch = [
            (gens_i_d[ n2s_gens_idx[nth_idx] ] * vh[ n2s_all_nodes_idx[nth_idx]] *
            sin(gens_δ[ n2s_gens_idx[nth_idx] ] - θh[ n2s_all_nodes_idx[nth_idx]])) +
            (gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ n2s_all_nodes_idx[nth_idx]] *
            cos(gens_δ[ n2s_gens_idx[nth_idx] ]- θh[ n2s_all_nodes_idx[nth_idx]]))  -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [
                vh[ n2s_all_nodes_idx[idx]] *
                    abs(ynj) *
                    cos(θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[idx]] -
                    angle(ynj))
                for (ynj, idx) in
                    zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                         nodes_idx_with_adjacent_nodes_idx[
                             n2s_all_nodes_idx[nth_idx] ])])
                            for nth_idx in gens_nodes_idx ]  
        
        
    end

    # -----------------------------------
    # gens reactive power mismatch
    # -----------------------------------

    """
    id_i * vh_i * cos(δ_i - θ_i) - iq_i * vh_i * sin(δ_i - θ_i) + Qloc_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    if loc_load_exist == true
        
        gens_Q_mismatch = [
            nth_idx ∈ gens_nodes_with_loc_loads_idx ?
                ((gens_i_d[ n2s_gens_idx[nth_idx]] *
                vh[ n2s_all_nodes_idx[nth_idx]] *
                cos(gens_δ[ n2s_gens_idx[nth_idx]] -
                θh[ n2s_all_nodes_idx[nth_idx]])) -
                (gens_i_q[ n2s_gens_idx[nth_idx]] *
                vh[ n2s_all_nodes_idx[nth_idx]] *
                sin(gens_δ[ n2s_gens_idx[nth_idx]] -
                θh[ n2s_all_nodes_idx[nth_idx]])) -
                Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx]] - 
                vh[ n2s_all_nodes_idx[nth_idx]] * sum( [
                    vh[n2s_all_nodes_idx[idx]] *
                       abs(ynj) *
                       sin( θh[ n2s_all_nodes_idx[nth_idx]] -
                       θh[ n2s_all_nodes_idx[idx]] -
                       angle(ynj) )
                    for (ynj, idx) in
                        zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])]) ) :
                                    ((gens_i_d[ n2s_gens_idx[ nth_idx]] *
                                    vh[ n2s_all_nodes_idx[nth_idx]] *
                                    cos(gens_δ[ n2s_gens_idx[ nth_idx]] -
                                    θh[ n2s_all_nodes_idx[nth_idx]])) -
                                    (gens_i_q[ n2s_gens_idx[ nth_idx]] *
                                    vh[ n2s_all_nodes_idx[nth_idx]] *
                                    sin(gens_δ[ n2s_gens_idx[nth_idx]] -
                                    θh[ n2s_all_nodes_idx[nth_idx]])) -
                                       vh[ n2s_all_nodes_idx[nth_idx]] *
                                       sum([vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                                       sin(θh[ n2s_all_nodes_idx[nth_idx]] -
                                       θh[n2s_all_nodes_idx[idx]] - angle(ynj))
                                            for (ynj, idx) in
                                                zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                                    nodes_idx_with_adjacent_nodes_idx[
                                                        n2s_all_nodes_idx[nth_idx]])])) 
            for nth_idx in gens_nodes_idx]        
        
    else

        gens_Q_mismatch = [
            (gens_i_d[ n2s_gens_idx[nth_idx]] *
                vh[ n2s_all_nodes_idx[nth_idx]] *
                cos(gens_δ[ n2s_gens_idx[nth_idx]] -
                θh[ n2s_all_nodes_idx[nth_idx]])) -
                (gens_i_q[ n2s_gens_idx[nth_idx]] *
                vh[ n2s_all_nodes_idx[nth_idx]] *
                sin(gens_δ[ n2s_gens_idx[nth_idx]] -
                θh[ n2s_all_nodes_idx[nth_idx]])) -
                vh[ n2s_all_nodes_idx[nth_idx]] * sum( [
                    vh[ n2s_all_nodes_idx[idx]] *
                        abs(ynj) *
                        sin( θh[ n2s_all_nodes_idx[
                            nth_idx]] -
                        θh[ n2s_all_nodes_idx[idx]] -
                        angle(ynj))
                    for (ynj, idx) in
                        zip(Ynet[ n2s_all_nodes_idx[
                            nth_idx]],
                       nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[
                                    nth_idx]])])
            for nth_idx in gens_nodes_idx ]
        
    end
    
    # -----------------------------------
    # non_gens active power mismatch        
    # -----------------------------------

    """
    Pl_i =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )
    """
        
    non_gens_P_mismatch = [
        P_non_gens[n2s_non_gens_idx[nth_idx]] +
            vh[n2s_all_nodes_idx[nth_idx]] * sum([
                vh[ n2s_all_nodes_idx[idx]] *
                    abs(ynj) *
                    cos(θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[idx]] -
                    angle(ynj))
                for (ynj, idx) in
                    zip(Ynet[n2s_all_nodes_idx[nth_idx]],
                        nodes_idx_with_adjacent_nodes_idx[
                            n2s_all_nodes_idx[nth_idx]])])
        for nth_idx in
                non_gens_nodes_idx ]
    
    # ------------------------------------
    # non_gens reactive power mismatch                
    # ------------------------------------

    """
    Ql_i =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    non_gens_Q_mismatch = [
        Q_non_gens[n2s_non_gens_idx[nth_idx]] +
            vh[n2s_all_nodes_idx[nth_idx]] * sum([
                vh[ n2s_all_nodes_idx[idx]] *
                    abs(ynj) *
                    sin(θh[n2s_all_nodes_idx[nth_idx]] -
                    θh[n2s_all_nodes_idx[idx]] -
                    angle(ynj))
                for (ynj, idx) in
                    zip(Ynet[n2s_all_nodes_idx[nth_idx]],
                        nodes_idx_with_adjacent_nodes_idx[
                            n2s_all_nodes_idx[nth_idx]])])
        for nth_idx in non_gens_nodes_idx ]
    
    # ------------------------------------

    ΔPQ_id_iq .=
        vcat(gens_P_mismatch,
             non_gens_P_mismatch,
             gens_Q_mismatch,                 
             non_gens_Q_mismatch,
             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag)

    return nothing
    
end


