# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247

####################################################


#-----------------------------------------------------
#-----------------------------------------------------
# generic model ode : Sauer
#-----------------------------------------------------
#-----------------------------------------------------
         
function dae_a_gen_generic_model_by_ext_idq_func!(
    res, dx, x, gen_vh_id_iq_V_ref_Tm_para, t;
    kwd_para =
        gen_ode_kwd_para )
    
    vh, i_d, i_q, V_ref, Tm =
        gen_vh_id_iq_V_ref_Tm_para

    
    (gen_para,
     avr_para,
     ωs) =
         kwd_para

    
    (H,
     D,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
        gen_para
    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
        avr_para
    
    δ, ω, ed_dash, eq_dash, E_fd, R_f, V_R = x

    res[1] = ω - ωs - dx[1]
    
    # res[2] = ( ωs/(2*H)) * (Tm - ed_dash * i_d - eq_dash *
    #     i_q - (X_q_dash - X_d_dash) * i_d * i_q -  D * (ω - ωs)) 
    #     - dx[2]
    
    res[2] = ( ωs/(2*H)) * (Tm - ed_dash * i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q) 
    
    res[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash) - dx[3]

    res[4] = (1/T_d_dash) *
        (E_fd - (X_d - X_d_dash) * i_d - eq_dash) - dx[4]

    res[5] = (1/Te) * (V_R - (Ke + Sevf(Ae, Be, E_fd)) *
        E_fd  ) - dx[5]
    
    res[6] = (1/Tf) * (Kf * E_fd /Tf -  R_f ) - dx[6]

    
    res[7] = (1/Ta) * ( Ka * R_f - Ka * Kf * E_fd /Tf +
        Ka * (V_ref - vh) - V_R) - dx[7]
    
    return nothing    

end



function dae_generic_model_by_ext_idq_func!(
    res, dx, x, gens_vh_id_iq_V_ref_Tm_para, t;
    kwd_para =
        kwd_para  )

    (;gens_para,
     avrs_para,
     ωs,
     state_vars_idx,
     ode_vh_id_iq_V_ref_Tm_Idx ) =
         kwd_para

    (;ode_vh_Idx,
     ode_id_Idx,
     ode_iq_Idx,
     ode_V_ref_Idx,
     ode_Tm_Idx ) =
         ode_vh_id_iq_V_ref_Tm_Idx

    #----------------------------------------
    
    gens_vh =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_vh_Idx]
    
    gens_i_d =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_id_Idx]
    
    gens_i_q =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_iq_Idx]
    
    V_ref =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_V_ref_Idx]
    
    Tm =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_Tm_Idx]

    #----------------------------------------
    
    # (H,
    #  D,

    #  ra,
     
    #  X_d,
    #  X_q,
     
    #  X_d_dash,
    #  X_q_dash,
     
    #  T_d_dash,
    #  T_q_dash ) =
    #      gens_para

    

    (;H,
     D,

     # ra,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gens_para,
             (:H,
              :D,
              # :ra,
              :X_d,
              :X_q,     
              :X_d_dash,
              :X_q_dash,     
              :T_d_dash,
              :T_q_dash))
       
    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
        avrs_para
    

    for (gen_vh,
         gen_i_d, gen_i_q,
         gen_V_ref, gen_Tm,
         gen_H, gen_D,
         gen_X_d, gen_X_q,
         gen_X_d_dash, gen_X_q_dash,
         gen_T_d_dash, gen_T_q_dash,
         avr_Ka, avr_Ta,
         avr_Ke, avr_Te,
         avr_Kf, avr_Tf,
         avr_Ae, avr_Be,
         avr_Tr,
         state_var_idx) in
        zip( gens_vh,
             gens_i_d, gens_i_q,
             V_ref, Tm,
             H, D,
             X_d, X_q,
             X_d_dash, X_q_dash,
             T_d_dash, T_q_dash,
             Ka, Ta,
             Ke, Te,
             Kf, Tf,
             Ae, Be,
             Tr,
             state_vars_idx )

        gen_vh_id_iq_V_ref_Tm_para = [gen_vh;
               gen_i_d;
               gen_i_q;
               gen_V_ref;
               gen_Tm] 
    
        gen_para = (gen_H, gen_D,
                    gen_X_d, gen_X_q,
                    gen_X_d_dash, gen_X_q_dash,
                    gen_T_d_dash, gen_T_q_dash, )
        
        avr_para = (avr_Ka, avr_Ta,
                    avr_Ke, avr_Te,
                    avr_Kf, avr_Tf,
                    avr_Ae, avr_Be,
                    avr_Tr)
        
        a_kwd_para = (;gen_para,
                      avr_para,
                      ωs)

        dae_a_gen_generic_model_by_ext_idq_func!(
            res[state_var_idx],
            dx[state_var_idx],
            x[state_var_idx],
            gen_vh_id_iq_V_ref_Tm_para,
            t;
            kwd_para =
                a_kwd_para  )

    end

    return nothing    

end


#-----------------------------------------------------
# generic model system dynamics
#-----------------------------------------------------


function dae_SC_generic_model_dynamics!(
    res,
    dx,
    x,
    dae_SC_generic_model_dynamics_para,
    t;
    kwd_para =
        dae_SC_generic_model_dynamics_kwd_para )

    #----------------------------------------

    (;ωs,
     loc_load_exist,
     state_vars_idx,
     SC_generic_model_states_idx_in_state_Idx,

     dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,

     id_iq_pg_vh_Idx,
     ωs_ωref0_vref0_porder0_Idx,

     SC_generic_model_vars_wt_i_dq_Idx_in_state,
     SC_generic_model_states_comp_idxs_in_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs,

     gens_para,
     avrs_para,
     Ynet_wt_nodes_idx_wt_adjacent_nodes ) =
         NamedTupleTools.select(
             kwd_para,
             (:ωs,
              :loc_load_exist,
              :state_vars_idx,
              :SC_generic_model_states_idx_in_state_Idx,
              
              :dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
              
              :id_iq_pg_vh_Idx,
              :ωs_ωref0_vref0_porder0_Idx,
              
              :SC_generic_model_vars_wt_i_dq_Idx_in_state,
              :SC_generic_model_states_comp_idxs_in_Idx,
              
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              :non_pre_ordered_pf_vars_Idxs,
              
              :gens_para,
              :avrs_para,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes) )

    #----------------------------------------
    
    # needed
    (;dyn_V_ref_Idx,
     dyn_Tm_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
             (:dyn_V_ref_Idx,
              :dyn_Tm_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             SC_generic_model_vars_wt_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    #----------------------------------------
    
    # (;δ_idx_in_state,
    #  ed_dash_idx_in_state,
    #  eq_dash_idx_in_state ) =
    #      NamedTupleTools.select(
    #          flux_decay_model_states_comp_idxs_in_Idx,
    #          (:δ_idx_in_state,
    #           :ed_dash_idx_in_state,
    #           :eq_dash_idx_in_state ))

    
    (δ_idx_in_state,
     ω_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state,
     E_fd_idx_in_state,
     R_f_idx_in_state,
     V_R_idx_in_state ) =
         NamedTupleTools.select(
             SC_generic_model_states_comp_idxs_in_Idx,
             (:δ_idx_in_state,
              :ω_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state,
              :E_fd_idx_in_state,
              :R_f_idx_in_state,
              :V_R_idx_in_state))
    
    # -------------------------------------

   (;
    n2s_slack_gens_idx,
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
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))  

    #----------------------------------------
    
    (; gens_vh_idxs,
     gens_θh_idxs,

     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,

     flat_vh_flat_θh_flat_id_iq_Idx) =     
         NamedTupleTools.select(
             non_pre_ordered_pf_vars_Idxs,
             (:gens_vh_idxs,
              :gens_θh_idxs,
              
              :non_gens_nodes_vh_idxs,
              :non_gens_nodes_θh_idxs,

              :flat_vh_flat_θh_flat_id_iq_Idx))
        
    (;flat_vh_Idx,
     flat_θh_Idx,
     flat_id_Idx,
     flat_iq_Idx) =
         NamedTupleTools.select(
             flat_vh_flat_θh_flat_id_iq_Idx,
             (:flat_vh_Idx,
              :flat_θh_Idx,
              :flat_id_Idx,
              :flat_iq_Idx))
    
    #----------------------------------------    
    
    (H,
     D,
     ra,
     X_d,
     X_q,     
     X_d_dash,
     X_q_dash,     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gens_para,
             (:H,
              :D,
              :ra,
              :X_d,
              :X_q,     
              :X_d_dash,
              :X_q_dash,     
              :T_d_dash,
              :T_q_dash))
    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
         NamedTupleTools.select(
             avrs_para,
             (:Ka, :Ta,
              :Ke, :Te,
              :Kf, :Tf,
              :Ae, :Be,
              :Tr))
    
    #----------------------------------------

   (Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx)) 
    
    #----------------------------------------

    # vh_θh_Idx_in_state =
    #     [vh_Idx_in_state;
    #      θh_Idx_in_state]

    #----------------------------------------

    state_vars = x[state_var_Idx_in_state]
    
    flat_vh  = x[vh_Idx_in_state]
    
    flat_θh  = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    
    gens_i_q = x[iq_Idx_in_state]
    
    #----------------------------------------
    
    gens_δ = x[δ_idx_in_state]
        
    gens_ed_dash = x[ed_dash_idx_in_state]
    
    gens_eq_dash = x[eq_dash_idx_in_state]
        
    #----------------------------------------
    
    gens_vh =
        flat_vh[
            gens_nodes_idx ]

    gens_θh =
        flat_θh[
            gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[
            non_gens_nodes_idx ]
            
    #----------------------------------------
    
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
    
    #----------------------------------------
        
    V_ref = dae_SC_generic_model_dynamics_para[
        dyn_V_ref_Idx]
    
    Tm = dae_SC_generic_model_dynamics_para[
        dyn_Tm_Idx]
    
    P_non_gens = dae_SC_generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = dae_SC_generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            dae_SC_generic_model_dynamics_para[
                dyn_Pll_Idx] 
        
        Q_g_loc_load =
            dae_SC_generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    
    ode_res_views = [view(res, idx)
                     for idx in state_vars_idx ]
    
    ode_dx_views = [view(dx, idx)
                    for idx in state_vars_idx ]
    
    ode_x_views  = [view( x, idx)
                    for idx in state_vars_idx ]
        
    #----------------------------------------
    # ode parameters and kwd para
    #----------------------------------------
    
    gens_ode_kwd_para = [
        (gen_para =
            (gen_H, gen_D,
             gen_X_d, gen_X_q,
             gen_X_d_dash, gen_X_q_dash,
             gen_T_d_dash, gen_T_q_dash ),
         avr_para =
             (avr_Ka, avr_Ta,
              avr_Ke, avr_Te,
              avr_Kf, avr_Tf,
              avr_Ae, avr_Be,
              avr_Tr),
         ωs = ωs )
        for (gen_H, gen_D,
             gen_X_d, gen_X_q,
             gen_X_d_dash, gen_X_q_dash,
             gen_T_d_dash, gen_T_q_dash,
             avr_Ka, avr_Ta,
             avr_Ke, avr_Te,
             avr_Kf, avr_Tf,
             avr_Ae, avr_Be,
             avr_Tr) in
            zip(H, D,
             X_d, X_q,
             X_d_dash, X_q_dash,
             T_d_dash, T_q_dash,
             Ka, Ta,
             Ke, Te,
             Kf, Tf,
             Ae, Be,
             Tr )]    
    
    gens_vh_id_iq_V_ref_Tm_para = [
        [a_vh, a_id, a_iq, a_vref, a_tm]
        for (a_vh, a_id, a_iq, a_vref, a_tm) in
            zip(gens_vh,
            gens_i_d,
            gens_i_q,
            V_ref,
                Tm) ]
    
    #----------------------------------------        
    #----------------------------------------
    # ode 
    #----------------------------------------
    #----------------------------------------
    
    for (ode_res, ode_dx, ode_x,
         gen_vh_id_iq_V_ref_Tm_para,
         gen_ode_kwd_para) in
        zip(ode_res_views,
            ode_dx_views,
            ode_x_views,
            gens_vh_id_iq_V_ref_Tm_para,
            gens_ode_kwd_para)
        
        dae_a_gen_generic_model_by_ext_idq_func!(
            ode_res,
            ode_dx,
            ode_x,
            gen_vh_id_iq_V_ref_Tm_para,
            t;
            kwd_para =
                gen_ode_kwd_para )
        
    end
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    # -----------------------------------
    # active power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """
    
    res[vh_Idx_in_state] .= [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * (gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -          
          sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]))  -
                                   P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]]) :
        (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]) ))
                        for nth_idx in all_nodes_idx ] 

    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    res[θh_Idx_in_state]  .= [ nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          ( vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -           
           sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]]) :
           (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])]))) 
                              for nth_idx in all_nodes_idx ]

    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """
    
    res[id_Idx_in_state] .= [
     gens_ed_dash[ n2s_gens_idx[ idx ]] -
     gens_vh[ n2s_gens_idx[idx]] * sin( gens_δ[n2s_gens_idx[idx]] - gens_θh[ n2s_gens_idx[idx]]) -
     ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
     X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """
    
    res[iq_Idx_in_state] .= [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
        gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[ n2s_gens_idx[idx]] - gens_θh[ n2s_gens_idx[idx]]) -
        ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
        X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  
    
    # ------------------------------------
    
    return nothing

end


#-----------------------------------------------------
# flux decay  model ode
#-----------------------------------------------------


function dae_a_gen_flux_decay!(
    res, dx, x, p, t;
    kwd_para =
        kwd_para )

    # (vh, θh, i_d, i_q, V_ref) = p
    
    (vh, θh, V_ref, Tm) = p

    (gens_para,
     avrs_para,
     ωs) =
         kwd_para
    
    (H,
     X_d, X_q, X_d_dash,
     T_d_dash ) =
        gens_para

    (Ka, Ta) =
        avrs_para
    
    δ, ω, ed_dash, eq_dash, E_fd, R_f, V_R = x

    # Ig = (eq_dash - vh * exp(im * θh)) / (im * X_d_dash)
    # I_gen_dq =  Ig * exp(im * ( -δ + π/2) )
    # i_d = real(I_gen_dq )
    # i_q = imag(I_gen_dq )


    i_d = (eq_dash - vh * cos(δ - θh)) / X_d_dash
    
    i_q =  (vh * sin(δ - θh)) / X_q   

    res[1] = ω - ωs - dx[1]
    
    res[2] =( ωs/(2*H)) * (Tm -  eq_dash * i_q -
        (X_q - X_d_dash) * i_d * i_q) - dx[2]

    res[3] = (1/T_d_dash) *
        (E_fd - eq_dash - (X_d - X_d_dash) * i_d ) - dx[3]
    
    res[4] = (1/Ta) * (  Ka * ( V_ref - vh) - E_fd ) - dx[4]

    return nothing    

end



function dae_a_gen_flux_decay!(
    res, dx, x, p, t;
    kwd_para =
        kwd_para )

    # (vh, θh, i_d, i_q, V_ref) = p
    
    (vh, θh, V_ref, Tm) = p

    (gens_para,
     avrs_para,
     ωs) =
         kwd_para
    
    (H,
     X_d, X_q, X_d_dash,
     T_d_dash ) =
        gens_para

    (Ka, Ta) =
        avrs_para
    
    δ, ω, eq_dash, E_fd = x

    # Ig = (eq_dash - vh * exp(im * θh)) / (im * X_d_dash)
    # I_gen_dq =  Ig * exp(im * ( -δ + π/2) )
    # i_d = real(I_gen_dq )
    # i_q = imag(I_gen_dq )


    i_d = (eq_dash - vh * cos(δ - θh)) / X_d_dash
    
    i_q =  (vh * sin(δ - θh)) / X_q   

    res[1] = ω - ωs - dx[1]
    
    res[2] =( ωs/(2*H)) * (Tm -  eq_dash * i_q -
        (X_q - X_d_dash) * i_d * i_q) - dx[2]

    res[3] = (1/T_d_dash) *
        (E_fd - eq_dash - (X_d - X_d_dash) * i_d ) - dx[3]
    
    res[4] = (1/Ta) * (  Ka * ( V_ref - vh) - E_fd ) - dx[4]

    return nothing    

end


function dae_flux_decay!(
    res, dx, x, p, t;
    kwd_para =
        kwd_para  )

    # (vh, θh, i_d, i_q, V_ref) = p
    
    (gens_vh, gens_θh,
     V_ref, Tm) = p

    (gens_para,
     avrs_para,
     state_vars_idx,
     ωs) =
         kwd_para
    
    (H,
     X_d, X_q, X_d_dash,
     T_d_dash ) =
        gens_para

    ( Ka, Ta ) =
        avrs_para
    

    for (gen_vh, gen_θh,
         gen_V_ref, gen_Tm,
         gen_H,  gen_X_d, gen_X_q,
         gen_X_d_dash, gen_T_d_dash,
         avr_Ka, avr_Ta,
         state_var_idx) in
        zip(gens_vh, gens_θh,
            V_ref, Tm,
            H,  X_d, X_q,
            X_d_dash, T_d_dash,
            Ka, Ta,
            state_vars_idx )

        a_p =
            (gen_vh,
             gen_θh,
             gen_V_ref,
             gen_Tm) 
    
        gen_para =
            (gen_H,
             gen_X_d,
             gen_X_q,
             gen_X_d_dash,
             gen_T_d_dash )
        
        avr_para = (avr_Ka,
                    avr_Ta)
        
        a_kwd_para = (gen_para,
                      avr_para,
                      ωs)

        dae_a_gen_flux_decay!(
            res[state_var_idx],
            dx[state_var_idx],
            x[state_var_idx],
            a_p, t;
            kwd_para =
                a_kwd_para  )

    end

    return nothing    

end


function dae_a_gen_flux_decay_by_ext_idq_func!(
    res,dx, x, gen_vh_id_iq_V_ref_Tm_para, t;
    kwd_para =
        gen_ode_kwd_para )
    
    vh, i_d, i_q, V_ref, Tm =
        gen_vh_id_iq_V_ref_Tm_para

    (gen_para,
     avr_para,
     ωs) =
         kwd_para
    
    (H,
     X_d,
     X_q,
     X_d_dash,
     T_d_dash ) =
        gen_para

    (Ka, Ta) =
        avr_para
    
    δ, ω, eq_dash, E_fd = x


    res[1] = ω - ωs - dx[1]
    
    res[2] =( ωs/(2*H)) * (Tm -  eq_dash * i_q -
        (X_q - X_d_dash) * i_d * i_q) - dx[2]

    res[3] = (1/T_d_dash) *
        (E_fd - eq_dash - (X_d - X_d_dash) * i_d ) - dx[3]
    
    res[4] = (1/Ta) * (  Ka * ( V_ref - vh) - E_fd ) - dx[4]
    
    return nothing    

end


function dae_flux_decay_by_ext_idq_func!(
    res, dx, x, gens_vh_id_iq_V_ref_Tm_para, t;
    kwd_para =
        kwd_para  )

    (;gens_para,
     avrs_para,
     ωs,
     state_vars_idx,
     ode_vh_id_iq_V_ref_Tm_Idx ) =
         kwd_para

    (;ode_vh_Idx,
     ode_id_Idx,
     ode_iq_Idx,
     ode_V_ref_Idx,
     ode_Tm_Idx ) =
         ode_vh_id_iq_V_ref_Tm_Idx

    #----------------------------------------
    
    gens_vh =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_vh_Idx]
    
    gens_i_d =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_id_Idx]
    
    gens_i_q =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_iq_Idx]
    
    V_ref =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_V_ref_Idx]
    
    Tm =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_Tm_Idx]

    #----------------------------------------
    
    (H,
     X_d,
     X_q,
     X_d_dash, T_d_dash ) =
        gens_para

    ( Ka, Ta ) =
        avrs_para
    

    for (gen_vh, gen_i_d, gen_i_q,
         gen_V_ref, gen_Tm,
         gen_H,  gen_X_d, gen_X_q,
         gen_X_d_dash, gen_T_d_dash,
         avr_Ka, avr_Ta,
         state_var_idx) in
        zip(gens_vh, gens_i_d, gens_i_q,
            V_ref, Tm,
            H,  X_d, X_q,
            X_d_dash, T_d_dash,
            Ka, Ta,
            state_vars_idx )

        gen_vh_id_iq_V_ref_Tm_para = [gen_vh;
               gen_i_d;
               gen_i_q;
               gen_V_ref;
               gen_Tm] 
    
        gen_para = (gen_H,
                    gen_X_d,
                    gen_X_q,
                    gen_X_d_dash,
                    gen_T_d_dash )
        
        avr_para = (avr_Ka,
                    avr_Ta)
        
        a_kwd_para = (;gen_para,
                      avr_para,
                      ωs)

        dae_a_gen_flux_decay_by_ext_idq_func!(
            res[state_var_idx],
            dx[state_var_idx],
            x[state_var_idx],
            gen_vh_id_iq_V_ref_Tm_para,
            t;
            kwd_para =
                a_kwd_para  )

    end

    return nothing    

end


#-----------------------------------------------------
# flux decay system dynamicss
#-----------------------------------------------------


function dae_flux_decay_model_dynamics!(
    res,
    dx,
    x,
    dae_flux_decay_model_dynamics_para,
    t;
    kwd_para =
        dae_flux_decay_model_dynamics_kwd_para  )

    #----------------------------------------

    (;ωs,
     loc_load_exist,
     state_vars_idx,
     dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
     
     flux_decay_model_vars_wt_i_dq_Idx_in_state,
     flux_decay_model_states_comp_idxs_in_Idx,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs,

     gens_para,
     avrs_para,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
          NamedTupleTools.select(kwd_para,
             (:ωs,
     :loc_load_exist,
     :state_vars_idx,
     :dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
     
     :flux_decay_model_vars_wt_i_dq_Idx_in_state,
     :flux_decay_model_states_comp_idxs_in_Idx,
     :dyn_pf_fun_kwd_n2s_idxs,
     :dyn_pf_fun_kwd_net_idxs,
     :non_pre_ordered_pf_vars_Idxs,

     :gens_para,
     :avrs_para,
     :Ynet_wt_nodes_idx_wt_adjacent_nodes) )

    #----------------------------------------    
    #----------------------------------------
    
    # needed
    (;dyn_V_ref_Idx,
     dyn_Tm_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
             (:dyn_V_ref_Idx,
              :dyn_Tm_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             flux_decay_model_vars_wt_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    #----------------------------------------
    
    (;δ_idx_in_state,
     ω_idx_in_state,
     eq_dash_idx_in_state,
     E_fd_idx_in_state ) =
         NamedTupleTools.select(
             flux_decay_model_states_comp_idxs_in_Idx,
             (:δ_idx_in_state,
              :ω_idx_in_state,
              :eq_dash_idx_in_state,
              :E_fd_idx_in_state ))

    # -------------------------------------

   (;
    n2s_slack_gens_idx,
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
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))  

    #----------------------------------------
    
    (;
     gens_vh_idxs,
     gens_θh_idxs,

     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,

     flat_vh_flat_θh_flat_id_iq_Idx) =     
         NamedTupleTools.select(
             non_pre_ordered_pf_vars_Idxs,
             (:gens_vh_idxs,
              :gens_θh_idxs,
              
              :non_gens_nodes_vh_idxs,
              :non_gens_nodes_θh_idxs,

              :flat_vh_flat_θh_flat_id_iq_Idx))
        
    (;flat_vh_Idx,
     flat_θh_Idx,
     flat_id_Idx,
     flat_iq_Idx) =
         NamedTupleTools.select(
             flat_vh_flat_θh_flat_id_iq_Idx,
             (:flat_vh_Idx,
              :flat_θh_Idx,
              :flat_id_Idx,
              :flat_iq_Idx))
    
    #----------------------------------------    
    #----------------------------------------
    
    (;H,
     X_d,
     X_q,
     X_d_dash,
     T_d_dash) =
         NamedTupleTools.select(
             gens_para,
             (:H,
              :X_d,
              :X_q,
              :X_d_dash,
              :T_d_dash))

    
    (;Ka,
     Ta) =
         NamedTupleTools.select(
             avrs_para,
             (:Ka,
              :Ta))
    
    #----------------------------------------

   ( Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx)) 

    #----------------------------------------
    
    # system_ode_para_kwd_para =
    #     (gens_para,
    #      avrs_para,
    #      state_vars_idx,
    #      ωs)

    #----------------------------------------
    
    (;H,
     X_d,
     X_q,
     X_d_dash,
     T_d_dash) =
         NamedTupleTools.select(
             gens_para,
             (:H,
              :X_d,
              :X_q,
              :X_d_dash,
              :T_d_dash))

    
    (;Ka,
     Ta) =
         NamedTupleTools.select(
             avrs_para,
             (:Ka,
              :Ta))
    
    #----------------------------------------

    vh_θh_Idx_in_state =
        [vh_Idx_in_state;
         θh_Idx_in_state]

    #----------------------------------------
    #----------------------------------------

    state_vars = x[state_var_Idx_in_state]
    
    flat_vh = x[vh_Idx_in_state]
    flat_θh = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    gens_i_q = x[iq_Idx_in_state]
    
    #----------------------------------------
    #----------------------------------------
    
    gens_δ = x[δ_idx_in_state]
    
    # gens_ω = x[ω_idx_in_state]
    
    gens_eq_dash = x[eq_dash_idx_in_state]
    
    # gens_E_fd = x[E_fd_idx_in_state]
    
    #----------------------------------------
    #----------------------------------------
    
    gens_vh =
        flat_vh[
            gens_nodes_idx ]

    gens_θh =
        flat_θh[
            gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[
            non_gens_nodes_idx ]
            
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
    
    #----------------------------------------
    #----------------------------------------
        
    V_ref = dae_flux_decay_model_dynamics_para[
        dyn_V_ref_Idx]
    
    Tm = dae_flux_decay_model_dynamics_para[
        dyn_Tm_Idx]
    
    P_non_gens = dae_flux_decay_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = dae_flux_decay_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            dae_flux_decay_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            dae_flux_decay_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    #----------------------------------------
    
    ode_res_views = [view(res, idx)
                     for idx in state_vars_idx ]
    
    ode_dx_views = [view(dx, idx)
                    for idx in state_vars_idx ]
    
    ode_x_views  = [view( x, idx)
                    for idx in state_vars_idx ]
    
    # #----------------------------------------
    # # PreallocationTools.jl
    # # https://github.com/SciML/PreallocationTools.jl
    # #----------------------------------------

    # stateDiffCache =
    #     get_tmp(stateDiffCache, x)

    # stateDiffCache .= x 
        
    #----------------------------------------
    # ode parameters and kwd para
    #----------------------------------------
    
    gens_ode_kwd_para = [
        (gen_para =
            (gen_H, gen_X_d,
             gen_X_q, gen_X_d_dash,
             gen_T_d_dash ),
         avr_para =
             (avr_Ka,
              avr_Ta),
         ωs = ωs )
        for (gen_H,  gen_X_d, gen_X_q,
             gen_X_d_dash, gen_T_d_dash,
             avr_Ka, avr_Ta) in
            zip( H,  X_d, X_q,
                 X_d_dash, T_d_dash,
                 Ka, Ta)]    
    
    gens_vh_id_iq_V_ref_Tm_para = [
        [a_vh, a_id, a_iq, a_vref, a_tm]
        for (a_vh, a_id, a_iq, a_vref, a_tm) in
            zip(gens_vh,
            gens_i_d,
            gens_i_q,
            V_ref,
                Tm) ]
    
    #----------------------------------------        
    #----------------------------------------
    # ode 
    #----------------------------------------
    #----------------------------------------
    
    for (ode_res, ode_dx, ode_x,
         gen_vh_id_iq_V_ref_Tm_para,
         gen_ode_kwd_para) in
        zip(ode_res_views,
            ode_dx_views,
            ode_x_views,
            gens_vh_id_iq_V_ref_Tm_para,
            gens_ode_kwd_para)
        
        dae_a_gen_flux_decay_by_ext_idq_func!(
            ode_res,
            ode_dx,
            ode_x,
            gen_vh_id_iq_V_ref_Tm_para,
            t;
            kwd_para =
                gen_ode_kwd_para )
    end

    # #--------------------------------------

    # # alternative
    
    # for (state_var_idx,
    #      gen_vh_id_iq_V_ref_Tm_para,
    #      gen_ode_kwd_para) in
    #     zip(state_vars_idx,
    #         gens_vh_id_iq_V_ref_Tm_para,
    #         gens_ode_kwd_para)
        
    #     dae_a_gen_flux_decay_by_ext_idq_func!(
    #         res[state_var_idx],
    #         dx[state_var_idx],
    #         x[state_var_idx],
    #         gen_vh_id_iq_V_ref_Tm_para,
    #         t;
    #         kwd_para =
    #             gen_ode_kwd_para )
    # end
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    # -----------------------------------
    # active power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
    #      (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
    #       vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
    #       cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
    #                for (ynj, idx) in
    #                    zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
    #                         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
    #                             nth_idx ∈ gens_nodes_with_loc_loads_idx ?
    #      (gens_i_d[ n2s_gens_idx[ nth_idx]] * vh[ n2s_all_nodes_idx[ nth_idx]] *
    #       sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
    #       gens_i_q[ n2s_gens_idx[ nth_idx]] * vh[ n2s_all_nodes_idx[ nth_idx]] * 
    #       cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
    #       P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
    #       vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
    #             cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
    #                   for (ynj, idx) in
    #                       zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
    #                          nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
    #     (gens_i_d[ n2s_gens_idx[nth_idx]] * vh[ n2s_all_nodes_idx[nth_idx]] * 
    #      sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
    #      gens_i_q[ n2s_gens_idx[nth_idx]] * vh[ n2s_all_nodes_idx[nth_idx]] * 
    #      cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
    #      vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
    #      cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
    #                    for (ynj, idx) in
    #                        zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
    #                          nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ]))
    #                     for nth_idx in all_nodes_idx ] 

    
    res[vh_Idx_in_state] .= [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * (gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -          
          sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]))  -
                                   P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]]) :
        (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]) ))
                        for nth_idx in all_nodes_idx ] 


    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------

    
    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # Q_mismatch  = [ nth_idx ∈ non_gens_nodes_idx ?
    #     ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
    #         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
    #         sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
    #                for (ynj, idx) in
    #                    zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
    #                         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
    #                             nth_idx ∈ gens_nodes_with_loc_loads_idx ?
    #       (gens_i_d[ n2s_gens_idx[nth_idx]] * vh[ n2s_all_nodes_idx[nth_idx]] *
    #        cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
    #        gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ n2s_all_nodes_idx[nth_idx]] *
    #        sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
    #        Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
    #        vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
    #        sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
    #                for (ynj, idx) in
    #                        zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
    #                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
    #        (gens_i_d[ n2s_gens_idx[nth_idx]] * vh[ n2s_all_nodes_idx[nth_idx]] *
    #         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
    #         gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ n2s_all_nodes_idx[nth_idx]] *
    #         sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
    #         vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
    #         sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
    #                    for (ynj, idx) in
    #                        zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
    #                           nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) 
    #                           for nth_idx in all_nodes_idx ]

    res[θh_Idx_in_state]  .= [ nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          ( vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -           
           sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]]) :
           (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])]))) 
                              for nth_idx in all_nodes_idx ]


    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - Xd_dash * id """
    
    res[id_Idx_in_state] .= [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - (
            gens_vh[ n2s_gens_idx[idx]] * cos(
                gens_δ[ n2s_gens_idx[idx]] -
                    gens_θh[ n2s_gens_idx[idx]])) - (
                        X_d_dash[ n2s_gens_idx[idx]] *
                            gens_i_d[ n2s_gens_idx[ idx]])
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ Xq * iq - vh * sin(δ - θh)  """
    
    res[iq_Idx_in_state] .= [
        X_q[ n2s_gens_idx[idx]] * gens_i_q[
            n2s_gens_idx[idx]] -
                gens_vh[ n2s_gens_idx[idx]] * sin(
                    gens_δ[n2s_gens_idx[idx]] -
                        gens_θh[ n2s_gens_idx[idx]])
                    
        for idx in gens_nodes_idx ]  

    # ------------------------------------

    # res[vh_Idx_in_state] .= P_mismatch
    
    # res[θh_Idx_in_state] .= Q_mismatch
    
    # res[id_Idx_in_state] .=
    #     gens_stator_equations_mismatch_real
    
    # res[iq_Idx_in_state] .=
    #     gens_stator_equations_mismatch_imag
    
    return nothing

end


# ---------------------------------------------------
# Comments
# ---------------------------------------------------


