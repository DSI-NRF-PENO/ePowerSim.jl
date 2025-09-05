# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247


####################################################


state_vars_syms_flux_decay = [:δ, :ω, :eq_dash, :E_fd ]

state_vars_syms_spcm = [:δ, :ω]

state_vars_syms_internal_mode = [:δ, :ω]


# ---------------------------------------------------
####################################################
# ---------------------------------------------------


function get_init_a_gen_classical_model(
    gen_vh,
    gen_θh,
    P_g,
    Q_g,
    X_d_dash )

    # E_q
    
    Ig = (P_g - im * Q_g) / (
        gen_vh * exp(-im * gen_θh))

    E_gen = gen_vh * exp(im * gen_θh) + im * X_d_dash * Ig

    E = abs( E_gen )

    δ = angle( E_gen )

    i_dq =  abs(Ig) * exp(im * (angle(Ig) - δ + π/2) )

    v_dq =  gen_vh * exp(im * (gen_θh - δ + π/2) )

    i_d = real(i_dq )

    i_q = imag(i_dq )

    v_d = real(v_dq )

    v_q = imag(v_dq )

    Tm = E * gen_vh * sin(δ - gen_θh) / X_d_dash

    return (; δ, E, i_d, i_q, v_d, v_q, Tm )

end


function get_init_gens_classical_model(
    gens_vh,
    gens_θh,
    pf_P_g_gens,
    pf_Q_g_gens,
    gens_Xd_dash)

    return [
        get_init_a_gen_classical_model(
            gen_vh, gen_θh, P_g, Q_g, Xd_dash )
        for (gen_vh, gen_θh, P_g, Q_g, 
             Xd_dash ) in
            zip(gens_vh, gens_θh,
                pf_P_g_gens, pf_Q_g_gens,
                gens_Xd_dash ) ]

end


function get_state_init_classical_model(
    init_δ,
    ωs)

    (δ, ) =
        init_δ

    return [[ [δ_i, ωs]
             for δ_i in δ]...;]
    
end


#-----------------------------------------------------
# model flux
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
    
    δ, ω, eq_dash, E_fd = x

    # Ig = (eq_dash - vh * exp(im * θh)) / (im * X_d_dash)
    # I_gen_dq =  Ig * exp(im * ( -δ + π/2) )
    # i_d = real(I_gen_dq )
    # i_q = imag(I_gen_dq )


    i_d = (eq_dash - vh * cos(δ - θh)) / X_d_dash
    
    i_q =  (vh * sin(δ - θh)) / X_q   

    res[1] = ω - ωs - dx[1]
    
    res[2] = (ωs/(2*H)) * (Tm -  eq_dash * i_q -
        (X_q - X_d_dash) * i_d * i_q) - dx[2]

    res[3] = (1/T_d_dash) *
        (E_fd - eq_dash - (X_d - X_d_dash) * i_d ) - dx[3]
    
    res[4] = (1/Ta) * (Ka * (V_ref - vh) - E_fd ) - dx[4]

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
    

    for (gen_vh, gen_θh, gen_V_ref, gen_Tm, gen_H,  gen_X_d, gen_X_q, gen_X_d_dash, gen_T_d_dash, avr_Ka, avr_Ta, state_var_idx) in zip(gens_vh, gens_θh, V_ref, Tm, H,  X_d, X_q, X_d_dash, T_d_dash, Ka, Ta, state_vars_idx )

        a_p =
            (gen_vh, gen_θh,
             gen_V_ref, gen_Tm) 
    
        gen_para =
            (gen_H,
             gen_X_d, gen_X_q, gen_X_d_dash,
             gen_T_d_dash )
        
        avr_para = (avr_Ka, avr_Ta)
        
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
    
    res[2] =(ωs/(2*H)) * (Tm -  eq_dash * i_q -
        (X_q - X_d_dash) * i_d * i_q) - dx[2]

    res[3] = (1/T_d_dash) *
        (E_fd - eq_dash - (X_d - X_d_dash) * i_d ) - dx[3]
    
    res[4] = (1/Ta) * (Ka * (V_ref - vh) - E_fd ) - dx[4]
    
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
    

    for (gen_vh, gen_i_d, gen_i_q, gen_V_ref, gen_Tm, gen_H,  gen_X_d, gen_X_q, gen_X_d_dash, gen_T_d_dash, avr_Ka, avr_Ta, state_var_idx) in zip(gens_vh, gens_i_d, gens_i_q, V_ref, Tm, H,  X_d, X_q, X_d_dash, T_d_dash, Ka, Ta, state_vars_idx )

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
    
    Tm = flux_decay_model_dynamics_para[
        dyn_Tm_Idx]
    
    P_non_gens = flux_decay_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = flux_decay_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            flux_decay_model_dynamics_para[
                dyn_Qll_Idx]
        
        Q_g_loc_load =
            flux_decay_model_dynamics_para[
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
    #                          nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ]))f
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

# ------------------------------------------------------
# ------------------------------------------------------
# Partitioned explicit
# ------------------------------------------------------
# ------------------------------------------------------

function get_flux_decay_model_id_iq(
    gens_vh,
    gens_θh,
    gens_δ,
    gens_eq_dash,
    gens_Xq,
    gens_Xd_dash )
    
    #-------------------------------
    
    gens_i_q = gens_vh * sin(gens_δ - gens_θh) /
        gens_Xq

    gens_i_d = ( gens_eq_dash - gens_vh * cos(
        gens_δ - gens_θh ) ) / gens_Xd_dash

    return (; gens_i_d, gens_i_q)
    
end


function get_flux_decay_model_id_iq_Pg_Qg(
    gens_vh, gens_θh, gens_δ, gens_eq_dash;
    kwd_para = nt_X_d_dash_X_q  )

    #-------------------------------
    
    (gens_Xd_dash,
     gens_Xq ) =
         NamedTupleTools.select(
             kwd_para,
             (:X_d_dash,
              :X_q ) )
    
    #-------------------------------
    
    gens_iq = gens_vh .* sin.(gens_δ - gens_θh) ./
        gens_Xq

    gens_id = ( gens_eq_dash - gens_vh .* cos.(
        gens_δ - gens_θh ) ) ./ gens_Xd_dash
    
    #----------------------------------------
    
    P_gens = get_gens_Pg_from_δ_i_dq_vh_θh(
        gens_vh, gens_θh, gens_δ,
        gens_id, gens_iq )
    
     Q_gens = get_gens_Qg_from_δ_i_dq_vh_θh(
         gens_vh, gens_θh, gens_δ,
         gens_id, gens_iq )

    # -------------------------------------

    return (; gens_iq,
            gens_id,
            P_gens,
            Q_gens )

end


#-----------------------------------------------------
# Reduced order models
#-----------------------------------------------------

# ---------------------------------------------------
#####################################################
# ---------------------------------------------------
# flux decay model 
# ---------------------------------------------------
#####################################################
# ---------------------------------------------------


# ---------------------------------------------------
# ---------------------------------------------------
# flux decay model dynamics
# ---------------------------------------------------
# ---------------------------------------------------

function get_flux_decay_model_Pq_Qg_id_iq_by_parts(
    init_flat_vh_flat_θh,
    flux_decay_dyn_pf_δ_eq_dash_flat_para,
    flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para;
    by_part_dyn_pf_mismatch_kwd_para =
        by_part_dyn_pf_mismatch_kwd_para )

    #----------------------------------------

    (; pf_solver,
     use_nlsolve,
     dyn_pf_mismatch_kwd_para,
     sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx )  =
         by_part_dyn_pf_mismatch_kwd_para
    
    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver

    (;loc_load_exist,
     dyn_pf_mismatch_vars_kwd_para,
     dyn_pf_mismatch_idx_kwd_para ) =
         dyn_pf_mismatch_kwd_para
    
    #-------------------------------
    
    ode_gens_para =
        getproperty(
            dyn_pf_mismatch_vars_kwd_para,
            :ode_gens_para )

    
    (gens_Xd_dash,
     gens_Xq ) =
         NamedTupleTools.select(
             ode_gens_para,
             (:X_d_dash,
              :X_q ) )
    
    #-------------------------------
    
    gens_nodes_idx =
        getproperty(
            getproperty(
                dyn_pf_mismatch_idx_kwd_para,
                :dyn_pf_fun_kwd_net_idxs),
            :gens_nodes_idx)
    
    non_pre_ordered_pf_vars_Idxs =
        getproperty(
            dyn_pf_mismatch_idx_kwd_para,
            :non_pre_ordered_pf_vars_Idxs )

    #-------------------------------
    
    (gens_vh_idxs,
     gens_θh_idxs,
     flat_vh_flat_θh_flat_id_iq_Idx) =
         NamedTupleTools.select(
             non_pre_ordered_pf_vars_Idxs,
             (:gens_vh_idxs,
              :gens_θh_idxs,
              :flat_vh_flat_θh_flat_id_iq_Idx ) )
    
    
    (flat_vh_Idx,
     flat_θh_Idx,
     flat_id_Idx,
     flat_iq_Idx) =
        flat_vh_flat_θh_flat_id_iq_Idx
    
    #-------------------------------
    
    # dyn_pf_δ_eq_dash_0_P_Q_idx
    # dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx
    
    dyn_pf_δ_eq_dash_Idx =
        getproperty(
            sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx,
            :dyn_pf_δ_eq_dash_Idx )
    
   ( dyn_pf_δ_Idxs,
     dyn_pf_eq_dash_Idxs ) =
         dyn_pf_δ_eq_dash_Idx
    
    #-------------------------------
    #-------------------------------
    
    gens_vh = init_flat_vh_flat_θh[
        flat_vh_Idx][gens_nodes_idx]

    gens_θh = init_flat_vh_flat_θh[
        flat_θh_Idx][gens_nodes_idx]
    
    #-------------------------------

    # gens_δ, gens_eq_dash

    gens_δ =
         flux_decay_dyn_pf_δ_eq_dash_flat_para[
            dyn_pf_δ_Idxs ]
    
    gens_eq_dash =
         flux_decay_dyn_pf_δ_eq_dash_flat_para[
            dyn_pf_eq_dash_Idxs ]
    
    #-------------------------------
    
    gens_iq = gens_vh .* sin.(gens_δ - gens_θh) ./
        gens_Xq

    gens_id = ( gens_eq_dash - gens_vh .* cos.(
        gens_δ - gens_θh ) ) ./ gens_Xd_dash
    
    #----------------------------------------
    
    P_gens = get_gens_Pg_from_δ_i_dq_vh_θh(
        gens_vh, gens_θh, gens_δ,
        gens_id, gens_iq )
    
     Q_gens = get_gens_Qg_from_δ_i_dq_vh_θh(
         gens_vh, gens_θh, gens_δ,
         gens_id, gens_iq )

    return (;P_gens, Q_gens, gens_id, gens_iq)
end

# ---------------------------------------------------
# flux decay model state init
# ---------------------------------------------------

function get_init_a_gen_flux_decay_model(
    gen_vh,
    gen_θh,
    P_g,
    Q_g,
    X_d,
    X_q,
    X_d_dash,
    X_q_dash,
    Ka )

    # E_q
    
    Ig = (P_g - im * Q_g) / (
        gen_vh * exp(-im * gen_θh))

    E_gen = gen_vh * exp(im * gen_θh) + im * X_q * Ig

    δ = angle( E_gen )

    i_dq =  abs(Ig) * exp(im * (angle(Ig) - δ + π/2) )

    v_dq =  gen_vh * exp(im * (gen_θh - δ + π/2) )

    i_d = real(i_dq )

    i_q = imag(i_dq )

    v_d = real(v_dq )

    v_q = imag(v_dq )

    eq_dash = gen_vh * cos( δ - gen_θh ) +
        X_d_dash * i_d

    E_fd = eq_dash + (X_d - X_d_dash) * i_d

    V_ref = ((E_fd / Ka) + gen_vh)

    Tm = eq_dash * i_q + (X_q - X_d_dash) * i_d * i_q

    return (;δ,
            eq_dash,
            E_fd,
            V_ref, Tm,
            i_d, i_q,
            v_d, v_q )

end


function get_init_flux_decay_model(
    gens_vh,
    gens_θh,
    pf_P_g_gens,
    pf_Q_g_gens,
    gens_Xd,
    gens_Xq,
    gens_Xd_dash,
    gens_Xq_dash,
    avrs_Ka)

    return [
        get_init_a_gen_flux_decay_model(
            gen_vh, gen_θh, P_g, Q_g,  Xd, Xq,
            Xd_dash, Xq_dash, Ka )
        for (gen_vh, gen_θh, P_g, Q_g, 
             Xd, Xq, Xd_dash, Xq_dash, Ka) in
            zip(gens_vh, gens_θh,
                pf_P_g_gens, pf_Q_g_gens,
                gens_Xd, gens_Xq,
                gens_Xd_dash, gens_Xq_dash, avrs_Ka ) ]

end


function get_state_init_flux_decay(
    init_flux_decay_δ_eq_E_fd,
    ωs)

    (δ, eq_dash, E_fd) =
        init_flux_decay_δ_eq_E_fd

    return [[ [δ_i, ωs, eq_dash_i, E_fd_i]
             for (δ_i, eq_dash_i, E_fd_i) in
                 zip(δ, eq_dash, E_fd) ]...;]
    
end


# ---------------------------------------------------
# flux decay voltage ternminal functions
# ---------------------------------------------------


function a_gen_flux_decay_model_vtf_func!(
    dx,
    x,
    gen_vh_θh,
    t )

    #----------------------------------------    

    #vh, θh = gen_vh_θh

    dx .= gen_vh_θh - x
    
    return nothing
    

end


function a_non_gen_flux_decay_model_vtf_func!(
    dx,
    x,
    non_gen_vh_θh,
    t )

    dx .= non_gen_vh_θh - x

    return nothing


end


function gen_flux_decay_model_vtf_func!(
    dx,
    x,
    post_pf_vec_gens_vh_θh,
    t;
    idx_in_state =
        vec_gens_vh_θh_Idx_in_state )

    for (idx, gen_vh_θh) in
        zip(vec_gens_vh_θh_Idx_in_state,
            post_pf_vec_gens_vh_θh)

        a_gen_flux_decay_model_vtf_func!(
            dx[idx],
            x[idx],
            gen_vh_θh,
            t)

    end
        
    return nothing
    
end

function non_gen_flux_decay_model_vtf_func!(
    dx,
    x,
    post_pf_vec_non_gens_vh_θh,
    t;
    idx_in_state =
        vec_non_gens_vh_θh_Idx_in_state )

    for (idx, non_gen_vh_θh) in
        zip(vec_non_gens_vh_θh_Idx_in_state,
            post_pf_vec_non_gens_vh_θh)

        a_non_gen_flux_decay_model_vtf_func!(
            dx[idx],
            x[idx],
            non_gen_vh_θh,
            t)

    end
        
    return nothing
    
end



function vtf_a_non_gen_flux_decay_model_func!(
    dx,
    x,
    non_gen_vh_θh,
    t
    )

    dx .= non_gen_vh_θh - x
    
    return nothing

end


function vtf_a_gen_flux_decay_model_func!(
    dx,
    x,
    gen_vh_θh_δ_eq_dash,
    t
    )

    #----------------------------------------    

    gen_vh, gen_θh, gen_δ, gen_eq_dash =
        gen_vh_θh_δ_eq_dash
    
    #----------------------------------------    

    α = [ sin(gen_δ) cos(gen_δ);
          -cos(gen_δ) sin(gen_δ)]

    gen_vd = gen_vh * sin(gen_δ - gen_θh)
    
    gen_vq = gen_vh * cos(gen_δ - gen_θh)

    # vd_vq = [gen_vd, gen_vq]

    # β =  α * vd_vq
    
    # γ =  cartesian_to_polar( β )
       
    # # v_dq =  gen_vh * exp(im * (gen_θh - gen_δ + π/2) )

    v_dq = gen_vd + im * gen_vq

    uh = v_dq * exp(im * (gen_δ - π/2) )
    
    γ =  cartesian_to_polar([real(uh), imag(uh)])
   
    dx .= γ - x
    
    return nothing
    
    # (gen_Xq,
    #  gens_Xd_dash) =
    #      NamedTupleTools.select(
    #          vtf_kwd_para,

    #          (:gen_X_q,
    #           :gen_X_d_dash))
        
    # gen_iq = gen_vh * sin(gen_δ - gen_θh) /
    #     gen_Xq

    # gen_id = ( gens_eq_dash - gens_vh * cos(
    #     gens_δ - gens_θh ) ) / gens_Xd_dash

    # id_iq = [gen_id, gen_iq]
    
    # zdq  = 
    
    # x_vh, x_θh = x[idx_in_state]

    # vh, θh = gen_vh_θh

    """
    
    #----------------------------------------

    β = [vh * cos(θh) - x_vh * cos(x_θh),
         vh * sin(θh) - x_vh * sin(x_θh)]

    dx[idx_in_state] .=  cartesian_to_polar( β )
    
    #----------------------------------------

    """
    

end



function all_nodes_flux_decay_model_vtf_func!(
    dx,
    x,
    post_pf_flat_vh_θh,
    t;
    vtf_all_nodes_kwd_para =
        vtf_all_nodes_kwd_para )

    (vh_Idx_in_state,
     θh_Idx_in_state,
     flat_vh_Idx,
     flat_θh_Idx) = vtf_all_nodes_kwd_para

    dx[vh_Idx_in_state] .=
        post_pf_flat_vh_θh[flat_vh_Idx] -
        x[vh_Idx_in_state]

    dx[θh_Idx_in_state] .=
        post_pf_flat_vh_θh[flat_θh_Idx] -
        x[θh_Idx_in_state]
    
    return nothing
    

end


#-----------------------------------------------------
# flux decay model ode
#-----------------------------------------------------


function ode_a_gen_flux_decay!(
    dx, x, p, t;
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

    dx[1] = ω - ωs
    
    dx[2] =( ωs/(2*H)) * (Tm -  eq_dash * i_q -
        (X_q - X_d_dash) * i_d * i_q)

    dx[3] = (1/T_d_dash) *
        (E_fd - eq_dash - (X_d - X_d_dash) * i_d )
    
    dx[4] = (1/Ta) * (Ka * ( V_ref - vh) - E_fd )

    return nothing    

end


function ode_flux_decay!(
    dx, x, p, t;
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
    

    for (gen_vh, gen_θh, gen_V_ref, gen_Tm, gen_H,  gen_X_d, gen_X_q, gen_X_d_dash, gen_T_d_dash, avr_Ka, avr_Ta, state_var_idx) in zip(gens_vh, gens_θh, V_ref, Tm, H,  X_d, X_q, X_d_dash, T_d_dash, Ka, Ta, state_vars_idx )

        a_p =
            (gen_vh, gen_θh,
             gen_V_ref, gen_Tm) 
    
        gen_para =
            (gen_H,
             gen_X_d, gen_X_q, gen_X_d_dash,
             gen_T_d_dash )
        
        avr_para = (avr_Ka, avr_Ta)
        
        a_kwd_para = (gen_para,
                      avr_para,
                      ωs)

        ode_a_gen_flux_decay!(
            dx[state_var_idx],
            x[state_var_idx],
            a_p, t;
            kwd_para =
                a_kwd_para  )

    end

    return nothing    

end


function ode_a_gen_flux_decay_by_ext_idq_func!(
    dx, x, gen_vh_id_iq_V_ref_Tm_para, t;
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

    dx[1] = ω - ωs
    
    dx[2] = ( ωs/(2 * H)) * (Tm - eq_dash * i_q -
        (X_q - X_d_dash) * i_d * i_q)

    dx[3] = (1/T_d_dash) *
        ( E_fd - (X_d - X_d_dash) * i_d  - eq_dash)
    
    dx[4] = (1/Ta) * (Ka * (V_ref - vh) - E_fd  )

    return nothing    

end


function ode_flux_decay_by_ext_idq_func!(
    dx, x, gens_vh_id_iq_V_ref_Tm_para, t;
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
    

    for (gen_vh, gen_i_d, gen_i_q, gen_V_ref, gen_Tm, gen_H,  gen_X_d, gen_X_q, gen_X_d_dash, gen_T_d_dash, avr_Ka, avr_Ta, state_var_idx) in zip(gens_vh, gens_i_d, gens_i_q, V_ref, Tm, H,  X_d, X_q, X_d_dash, T_d_dash, Ka, Ta, state_vars_idx )

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

        ode_a_gen_flux_decay_by_ext_idq_func!(
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

"""

 use `sd_dynamics_model_by_per_gen_by_vh_θh_func!` or 
`sd_dynamics_model_by_per_gen_by_ur_ui_func!`
 as an example

"""


function flux_decay_model_dynamics!(
    dx,
    x,
    flux_decay_model_dynamics_para,
    t;
    kwd_para =
        flux_decay_model_dynamics_kwd_para  )

    #----------------------------------------

    (;dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
      stateDiffCache,
      ode_flux_decay_kwd_para,
      by_part_dyn_pf_mismatch_kwd_para,
      sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx,
      flux_decay_model_states_Idx_wt_idxs_in_Idx) =
          NamedTupleTools.select(kwd_para,
             (:dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
              :stateDiffCache,
              :ode_flux_decay_kwd_para,
              :by_part_dyn_pf_mismatch_kwd_para,
              :sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx,
              :flux_decay_model_states_Idx_wt_idxs_in_Idx) )

    #----------------------------------------

    loc_load_exist =
        getproperty(
            getproperty(
                by_part_dyn_pf_mismatch_kwd_para,
                :dyn_pf_mismatch_kwd_para),
            :loc_load_exist)

    #----------------------------------------

    (slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx) =
         NamedTupleTools.select(
             getproperty(
                 getproperty(
                     getproperty(
                         by_part_dyn_pf_mismatch_kwd_para,
                         :dyn_pf_mismatch_kwd_para),
                     :dyn_pf_mismatch_idx_kwd_para),
                 :dyn_pf_fun_kwd_net_idxs),
             (:slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx)) 
    
    #----------------------------------------

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

    V_ref = flux_decay_model_dynamics_para[
        dyn_V_ref_Idx]
    
    Tm = flux_decay_model_dynamics_para[
        dyn_Tm_Idx]
    
    Png = flux_decay_model_dynamics_para[
        dyn_Png_Idx]
    
    Qng = flux_decay_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true    
        Pll = flux_decay_model_dynamics_para[
            dyn_Qll_Idx]
        
        Qll = flux_decay_model_dynamics_para[
            dyn_Qll_Idx]        
    else
    
    Pll = []
    Qll = []
        
    end
    
    #----------------------------------------
    
    (gens_para,
     avrs_para,
     state_vars_idx) =
         NamedTupleTools.select(
                 ode_flux_decay_kwd_para,
             (:ode_flux_decay_gens_para,
              :ode_flux_decay_avrs_para,
              :state_vars_idx) )

    system_ode_para_kwd_para =
        (gens_para,
         avrs_para,
         state_vars_idx,
         ωs)

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

    (dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
     non_pre_ordered_pf_vars_Idxs) =
         NamedTupleTools.select(        
             getproperty(
                 getproperty(
                     by_part_dyn_pf_mismatch_kwd_para,
                     :dyn_pf_mismatch_kwd_para),
                 :dyn_pf_mismatch_idx_kwd_para),
             (:dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
              :non_pre_ordered_pf_vars_Idxs) )

    (;
     gens_vh_idxs,
     gens_θh_idxs,
     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,
     gens_id_idxs,
     gens_iq_idxs,     
     flat_vh_flat_θh_flat_id_iq_Idx) =
         NamedTupleTools.select(
             non_pre_ordered_pf_vars_Idxs,
             (:gens_vh_idxs,
              :gens_θh_idxs,
              :non_gens_nodes_vh_idxs,
              :non_gens_nodes_θh_idxs,
              :gens_id_idxs,
              :gens_iq_idxs,     
              :flat_vh_flat_θh_flat_id_iq_Idx))
          
    
    (flat_vh_Idx,
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
    
   (;flux_decay_model_states_Idx,
    flux_decay_model_states_comp_idxs_in_Idx,
    flux_decay_model_vars_Idx_in_state ) =
        NamedTupleTools.select(
            flux_decay_model_states_Idx_wt_idxs_in_Idx,
            (:flux_decay_model_states_Idx,
             :flux_decay_model_states_comp_idxs_in_Idx,
             :flux_decay_model_vars_Idx_in_state))
    
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
    
    (;state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state ) =
         NamedTupleTools.select(
             flux_decay_model_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :vh_Idx_in_state,
              :θh_Idx_in_state))


    vh_θh_Idx_in_state =
        [vh_Idx_in_state;
         θh_Idx_in_state]

    #----------------------------------------

    vec_vh_θh_Idx_in_state =
        [[a_vh_Idx, a_θh_Idx]
         for (a_vh_Idx, a_θh_Idx) in
             zip(vh_Idx_in_state,
                 θh_Idx_in_state) ]

    vec_gens_vh_θh_Idx_in_state =
        vec_vh_θh_Idx_in_state[
            gens_nodes_idx ]

    vec_non_gens_vh_θh_Idx_in_state =
        vec_vh_θh_Idx_in_state[
            non_gens_nodes_idx]

    #----------------------------------------

    ode_dx_views = [view(dx, idx) for idx in state_vars_idx ]
    ode_x_views  = [view( x, idx) for idx in state_vars_idx ]
    
    #----------------------------------------
   
    vtf_dx_gens_vh_θh_views = [
        view( dx, a_vh_θh_Idx )
        for a_vh_θh_Idx in
            vec_gens_vh_θh_Idx_in_state ]
    
    vtf_x_gens_vh_θh_views = [
        view( x, a_vh_θh_Idx )
        for  a_vh_θh_Idx in
           vec_gens_vh_θh_Idx_in_state ]

    #----------------------------------------

    vtf_dx_non_gens_vh_θh_views = [
        view( dx, a_vh_θh_Idx )
        for a_vh_θh_Idx in
            vec_non_gens_vh_θh_Idx_in_state ]
    
    vtf_x_non_gens_vh_θh_views = [
        view( x, a_vh_θh_Idx )
        for  a_vh_θh_Idx in
            vec_non_gens_vh_θh_Idx_in_state ]

    #----------------------------------------
    # PreallocationTools.jl
    # https://github.com/SciML/PreallocationTools.jl
    #----------------------------------------

    stateDiffCache =
        get_tmp(stateDiffCache, x)

    stateDiffCache .= x 
    
    #----------------------------------------
    
    δ = stateDiffCache[
        δ_idx_in_state]
    
    eq_dash = stateDiffCache[
            eq_dash_idx_in_state]
    
    #----------------------------------------
    
    flat_vh = stateDiffCache[
        vh_Idx_in_state]
    
    flat_θh = stateDiffCache[
        θh_Idx_in_state ]

    #----------------------------------------
        
    gens_vh = flat_vh[ gens_nodes_idx ]
    
    gens_θh = flat_θh[ gens_nodes_idx ]
    
    #----------------------------------------

    vh_θh = [flat_vh;
             flat_θh ]

    #----------------------------------------
    
    gens_i_d = (eq_dash - gens_vh .* cos.(
        δ - gens_θh)) ./ X_d_dash

    gens_i_q =
        (gens_vh .* sin.(δ - gens_θh)) ./ X_q

    #----------------------------------------
    
    # gens_Pg_inj =
    #     get_gens_nodes_Pg_net_inj(
    #         δ,
    #         gens_vh,
    #         gens_θh,
    #         gens_i_d,
    #         gens_i_q,
    #         P_g_loc_load,
    #         Q_g_loc_load,

    #         n2s_gens_idx,
    #         gens_nodes_with_loc_loads_idx,

    #         gens_nodes_idx)
    
    # gens_Qg_inj =
    #     get_gens_nodes_Qg_net_inj(
    #         δ,
    #         gens_vh,
    #         gens_θh,
    #         gens_i_d,
    #         gens_i_q,
    #         P_g_loc_load,
    #         Q_g_loc_load,

    #         n2s_gens_idx,
    #         gens_nodes_with_loc_loads_idx,

    #         gens_nodes_idx)

    
    Pg =
        get_gens_Pg_from_δ_i_dq_vh_θh(
            gens_vh, gens_θh,
            δ,
            gens_i_d, gens_i_q )


    Qg =
        get_gens_Qg_from_δ_i_dq_vh_θh(
            gens_vh, gens_θh,
            δ,
            gens_i_d, gens_i_q )

    
    gens_Pg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        Pg[ n2s_gens_idx[ idx]] -
        P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx] ] :
        Pg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]

    gens_Qg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        Qg[ n2s_gens_idx[ idx]] -
        Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx] ] :
        Qg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]
    
    
    dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll =
        [gens_Pg_inj;
         gens_Qg_inj;
         Png;
         Qng;
         Pll;
         Qll]
    
    #----------------------------------------
    
    pf_sol = get_flux_decay_model_pf_ΔPQ_func!(
        vh_θh,
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll,
        by_part_dyn_pf_mismatch_kwd_para =
            by_part_dyn_pf_mismatch_kwd_para )

    #----------------------------------------


    # (pf_solver,
    #  dyn_pf_mismatch_kwd_para) =
    #     NamedTupleTools.select(
    #         by_part_dyn_pf_mismatch_kwd_para,
    #         (:pf_solver,
    #          :dyn_pf_mismatch_kwd_para ) )

    # (pf_alg,
    #  abstol,
    #  reltol) =
    #     NamedTupleTools.select(
    #         pf_solver,
    #         (:pf_alg,
    #          :abstol,
    #          :reltol ) )
    
    # #----------------------------------------

    # pf_sol = NonlinearSolve.solve(
    #     NonlinearProblem(
    #         NonlinearFunction( ( g, x, p ) ->
    #             get_a_flux_decay_model_pf_ΔPQ_mismatch!(
    #                 g, x, p;
    #                 flux_decay_dyn_pf_mismatch_kwd_para =
    #                     dyn_pf_mismatch_kwd_para ) ),
    #         stateDiffCache[vh_θh_Idx_in_state],
    #         dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll ),
    #     pf_alg )
        
    #----------------------------------------

    post_pf_vh = pf_sol[ flat_vh_Idx ]
    
    # post_pf_θh = pf_sol[ flat_θh_Idx ]

    post_pf_θh =
        normalise_angle_θh(
            pf_sol[flat_θh_Idx ],
            slack_gens_nodes_idx )

    post_pf_vec_vh_θh =
        [ [vh, θh]
          for (vh, θh) in
              zip( post_pf_vh,
                   post_pf_θh)]

    #----------------------------------------
    
    post_pf_vec_gens_vh_θh =
        post_pf_vec_vh_θh[
            gens_nodes_idx ]

    post_pf_vec_non_gens_vh_θh =
        post_pf_vec_vh_θh[
            non_gens_nodes_idx ]

    #----------------------------------------
    
    post_pf_gens_vh = post_pf_vh[
        gens_nodes_idx ]
    
    post_pf_gens_θh = post_pf_θh[
        gens_nodes_idx ]

    post_pf_non_gens_vh = post_pf_vh[
        non_gens_nodes_idx ]
    
    post_pf_non_gens_θh = post_pf_θh[
        non_gens_nodes_idx ]

    #----------------------------------------
    
    post_pf_id = ( eq_dash - post_pf_gens_vh .* cos.(
        δ - post_pf_gens_θh )) ./ X_d_dash

    post_pf_iq =
        ( post_pf_gens_vh .* sin.(
            δ - post_pf_gens_θh )) ./ X_q
    
    #----------------------------------------
    # ode
    #----------------------------------------

    for (ode_dx_view, ode_x_view, gen_vh, gen_i_d, gen_i_q,
         gen_V_ref, gen_Tm,
         gen_H,  gen_X_d,
         gen_X_q, gen_X_d_dash,
         gen_T_d_dash,
         avr_Ka, avr_Ta,
         state_var_idx) in
        zip(ode_dx_views, ode_x_views,
            gens_vh, post_pf_id, post_pf_iq,
            V_ref, Tm, H,
            X_d, X_q, X_d_dash, T_d_dash,
            Ka, Ta,
            state_vars_idx )

        a_δ, a_ω, a_eq_dash, a_E_fd = x[state_var_idx]
        
        an_ode_para =
            (gen_vh,
             gen_i_d, gen_i_q,
             gen_V_ref, gen_Tm)
        
        a_gen_kwd_para =
            (gen_H,
             gen_X_d, gen_X_q,
             gen_X_d_dash, gen_T_d_dash )
        
        an_avr_kwd_para =
            (avr_Ka,
             avr_Ta)
        
        a_kwd_para =
            (a_gen_kwd_para,
             an_avr_kwd_para,
             ωs)

        ode_a_gen_flux_decay!(
            ode_dx_view,
            ode_x_view,
            an_ode_para,
            t;
            kwd_para =
                a_kwd_para )
        
    end
    
    #----------------------------------------
    # gens vtf
    #----------------------------------------

    for (vtf_dx_gen_vh_θh_view,
         vtf_x_gen_vh_θh_view,
         a_gen_vh_θh) in
        zip(vtf_dx_gens_vh_θh_views,
            vtf_x_gens_vh_θh_views,
            post_pf_vec_gens_vh_θh )

        a_gen_flux_decay_model_vtf_func!(
            vtf_dx_gen_vh_θh_view,
            vtf_x_gen_vh_θh_view,
            a_gen_vh_θh,
            t )        
    end

    #----------------------------------------
    # non-gens vtf
    #----------------------------------------

    for (vtf_dx_non_gen_vh_θh_view,
         vtf_x_non_gen_vh_θh_view,
         a_vh_θh_Idx,
         a_non_gen_vh_θh) in
        zip(vtf_dx_non_gens_vh_θh_views,
            vtf_x_non_gens_vh_θh_views,
            vec_non_gens_vh_θh_Idx_in_state,
             post_pf_vec_non_gens_vh_θh )

        a_non_gen_flux_decay_model_vtf_func!(
            vtf_dx_non_gen_vh_θh_view,
            vtf_x_non_gen_vh_θh_view,
            a_non_gen_vh_θh,
            t )        
    end

    return nothing

end


function flux_decay_model_ode_wt_pf_vtf_dynamics!(
    dx,
    x,
    flux_decay_model_dynamics_para,
    t;
    kwd_para =
        model_dynamics_kwd_para )

    
    (; # flux_decay_model_states_init,
       # counter_array,
     
     loc_load_exist,
     
     dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
     stateDiffCache,
     by_part_dyn_pf_mismatch_kwd_para,
    
     slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx,

     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     dyn_V_ref_Idx,
     dyn_Tm_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx,
     
     gens_para,
     avrs_para,
     state_vars_idx,
     
     system_ode_para_kwd_para,
     
     H,
     X_d,
     X_q,
     X_d_dash,
     T_d_dash,
     
     Ka,
     Ta,
     
     dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
     non_pre_ordered_pf_vars_Idxs,
     
     gens_vh_idxs,
     gens_θh_idxs,
     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,
     
     gens_id_idxs,
     gens_iq_idxs,     
     flat_vh_flat_θh_flat_id_iq_Idx,
     
     flat_vh_Idx,
     flat_θh_Idx,
     flat_id_Idx,
     flat_iq_Idx,
     
     flux_decay_model_states_Idx,
     flux_decay_model_states_comp_idxs_in_Idx,
     flux_decay_model_vars_Idx_in_state,
     
     δ_idx_in_state,
     ω_idx_in_state,
     eq_dash_idx_in_state,
     E_fd_idx_in_state,
     
     state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,

     ode_vh_id_iq_V_ref_Tm_Idx,
     flux_decay_ode_para_kwd_para) =
        model_dynamics_kwd_para
    
    #----------------------------------------
    
    V_ref = flux_decay_model_dynamics_para[
        dyn_V_ref_Idx]
    
    Tm = flux_decay_model_dynamics_para[
        dyn_Tm_Idx]
    
    Png = flux_decay_model_dynamics_para[
        dyn_Png_Idx]
    
    Qng = flux_decay_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true    
        Pll = flux_decay_model_dynamics_para[
            dyn_Qll_Idx]
        
        Qll = flux_decay_model_dynamics_para[
            dyn_Qll_Idx]        
    else
    
    Pll = []
    Qll = []
        
    end
    
    #----------------------------------------

    vh_θh_Idx_in_state =
        [vh_Idx_in_state;
         θh_Idx_in_state]

    #----------------------------------------

    vec_vh_θh_Idx_in_state =
        [[a_vh_Idx, a_θh_Idx]
         for (a_vh_Idx, a_θh_Idx) in
             zip(vh_Idx_in_state,
                 θh_Idx_in_state) ]

    vec_gens_vh_θh_Idx_in_state =
        vec_vh_θh_Idx_in_state[
            gens_nodes_idx ]

    vec_non_gens_vh_θh_Idx_in_state =
        vec_vh_θh_Idx_in_state[
            non_gens_nodes_idx]

    #----------------------------------------

    ode_dx_views = [view(dx, idx)
                    for idx in
                        state_vars_idx ]

    ode_x_views  = [view( x, idx)
                    for idx in
                        state_vars_idx ]
    
    #----------------------------------------
   
    vtf_dx_gens_vh_θh_views = [
        view( dx, a_vh_θh_Idx )
        for a_vh_θh_Idx in
            vec_gens_vh_θh_Idx_in_state ]
    
    vtf_x_gens_vh_θh_views = [
        view( x, a_vh_θh_Idx )
        for  a_vh_θh_Idx in
           vec_gens_vh_θh_Idx_in_state ]

    #----------------------------------------

    vtf_dx_non_gens_vh_θh_views = [
        view( dx, a_vh_θh_Idx )
        for a_vh_θh_Idx in
            vec_non_gens_vh_θh_Idx_in_state ]
    
    vtf_x_non_gens_vh_θh_views = [
        view( x, a_vh_θh_Idx )
        for  a_vh_θh_Idx in
            vec_non_gens_vh_θh_Idx_in_state ]

    #----------------------------------------
    
    vtf_dx_vh_views = @view dx[ vh_Idx_in_state] 

    vtf_x_vh_views  = @view x[ vh_Idx_in_state] 

    vtf_dx_θh_views = @view dx[ θh_Idx_in_state]
    
    vtf_x_θh_views  = @view x[ θh_Idx_in_state]
    
    #----------------------------------------
    # PreallocationTools.jl
    # https://github.com/SciML/PreallocationTools.jl
    #----------------------------------------
        
    # counter = counter_array[1]

    #--------------------------------------------

    # stateDiffCache = get_tmp(stateDiffCache, x)

    # if counter == 1
    #     stateDiffCache .= flux_decay_model_states_init 
    #     flat_vh = stateDiffCache[ vh_Idx_in_state]
    #     flat_θh = stateDiffCache[ θh_Idx_in_state ]
    #     counter_array .+= 1
    # else
    #     stateDiffCache .= x
    #     flat_vh = stateDiffCache[ vh_Idx_in_state]
    #     flat_θh = stateDiffCache[ θh_Idx_in_state ]
    # end

    #----------------------------------------
    
    δ = x[ δ_idx_in_state]
    
    eq_dash = x[ eq_dash_idx_in_state]

    @show  δ
    
    @show eq_dash
    
    #----------------------------------------
            
    flat_vh = x[ vh_Idx_in_state]
    
    flat_θh = x[ θh_Idx_in_state ]

    vh_θh   = [flat_vh; flat_θh ]
    
    #----------------------------------------
        
    gens_vh = flat_vh[ gens_nodes_idx ]
    
    gens_θh = flat_θh[ gens_nodes_idx ]
    
    #----------------------------------------
     
    gens_i_d = (eq_dash - gens_vh .* cos.(
        δ - gens_θh)) ./ X_d_dash

    gens_i_q =
        (gens_vh .* sin.(δ - gens_θh)) ./ X_q

    #----------------------------------------
    
    # gens_Pg_inj =
    #     get_gens_nodes_Pg_net_inj(
    #         δ,
    #         gens_vh,
    #         gens_θh,
    #         gens_i_d,
    #         gens_i_q,
    #         P_g_loc_load,
    #         Q_g_loc_load,

    #         n2s_gens_idx,
    #         gens_nodes_with_loc_loads_idx,

    #         gens_nodes_idx)
    
    # gens_Qg_inj =
    #     get_gens_nodes_Qg_net_inj(
    #         δ,
    #         gens_vh,
    #         gens_θh,
    #         gens_i_d,
    #         gens_i_q,
    #         P_g_loc_load,
    #         Q_g_loc_load,

    #         n2s_gens_idx,
    #         gens_nodes_with_loc_loads_idx,

    #         gens_nodes_idx)

    #--------------------------------------    
    
    # Pg =
    #     get_gens_Pg_from_δ_i_dq_vh_θh(
    #         gens_vh, gens_θh,
    #         δ,
    #         gens_i_d, gens_i_q )

    # Qg =
    #     get_gens_Qg_from_δ_i_dq_vh_θh(
    #         gens_vh, gens_θh,
    #         δ,
    #         gens_i_d, gens_i_q )
    
    #--------------------------------------    
    
    Pg = gens_i_d .* gens_vh .* sin.(δ - gens_θh) +
         gens_i_q .* gens_vh .* cos.(δ - gens_θh)
    
    Qg = gens_i_d .* gens_vh .* cos.(δ - gens_θh) -
        gens_i_q .* gens_vh .* sin.(δ - gens_θh)
    

    #--------------------------------------    
    
    gens_Pg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        Pg[ n2s_gens_idx[ idx]] -
        P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] :
        Pg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]

    gens_Qg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        Qg[ n2s_gens_idx[ idx]] -
        Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] :
        Qg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]
    
    
    dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll =
        [gens_Pg_inj;
         gens_Qg_inj;
         Png;
         Qng;
         Pll;
         Qll]
    
    #----------------------------------------
    
    pf_sol = get_flux_decay_model_pf_ΔPQ_func!(
        vh_θh,
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll,
        by_part_dyn_pf_mismatch_kwd_para =
            by_part_dyn_pf_mismatch_kwd_para )

    # pf_sol = 1.0 .* vh_θh
    
    #----------------------------------------
    #----------------------------------------

    post_pf_vh = pf_sol[ flat_vh_Idx ]

    # post_pf_θh = pf_sol[ flat_θh_Idx ]
    
    pre_norm_post_pf_θh = pf_sol[ flat_θh_Idx ]

    post_pf_θh =
        normalise_angle_θh(
             pre_norm_post_pf_θh,
             slack_gens_nodes_idx )

    @show post_pf_vh

    @show post_pf_θh
    
    post_pf_flat_vh_θh =
        [post_pf_vh; post_pf_θh]
    
    post_pf_vec_vh_θh =
        [ [vh, θh]
          for (vh, θh) in
              zip( post_pf_vh,
                   post_pf_θh)]
     
    #----------------------------------------
    
    post_pf_vec_gens_vh_θh =
        post_pf_vec_vh_θh[
            gens_nodes_idx ]

    post_pf_vec_non_gens_vh_θh =
        post_pf_vec_vh_θh[
            non_gens_nodes_idx ]

    #----------------------------------------
    
    post_pf_gens_vh = post_pf_vh[
        gens_nodes_idx ]
    
    post_pf_gens_θh = post_pf_θh[
        gens_nodes_idx ]

    post_pf_non_gens_vh = post_pf_vh[
        non_gens_nodes_idx ]
    
    post_pf_non_gens_θh = post_pf_θh[
        non_gens_nodes_idx ]

    #----------------------------------------
    
    post_pf_id = ( eq_dash - post_pf_gens_vh .* cos.(
        δ - post_pf_gens_θh )) ./ X_d_dash

    post_pf_iq =
        ( post_pf_gens_vh .* sin.(
            δ - post_pf_gens_θh )) ./ X_q
    
    #----------------------------------------    
    #----------------------------------------
    # case 0
    #----------------------------------------
    #----------------------------------------
    
    (;gens_para,
     avrs_para,
     ωs) =
         NamedTupleTools.select(
             flux_decay_ode_para_kwd_para,
             (:gens_para,
              :avrs_para,
              :ωs))
    
    (H,
     X_d,
     X_q,
     X_d_dash,
     T_d_dash ) =
         NamedTupleTools.select(
             gens_para,
             (:H,
              :X_d,
              :X_q,
              :X_d_dash,
              :T_d_dash))

    ( Ka,
      Ta ) =
          NamedTupleTools.select(
              avrs_para,
              (:Ka,
               :Ta))
    
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
            zip(post_pf_gens_vh,
            post_pf_id,
            post_pf_iq,
            V_ref,
            Tm) ]

    for (ode_dx, ode_x,
         gen_vh_id_iq_V_ref_Tm_para,
         gen_ode_kwd_para) in
        zip(ode_dx_views,
            ode_x_views,
            gens_vh_id_iq_V_ref_Tm_para,
            gens_ode_kwd_para)
        
        ode_a_gen_flux_decay_by_ext_idq_func!(
            ode_dx,
            ode_x,
            gen_vh_id_iq_V_ref_Tm_para,
            t;
            kwd_para =
                gen_ode_kwd_para )
    end
    
    #----------------------------------------
    # all nodes vtf
    #----------------------------------------

    vtf_dx_θh_views .=
        post_pf_flat_vh_θh[flat_θh_Idx] -
        vtf_x_θh_views

    vtf_dx_vh_views .=
        post_pf_flat_vh_θh[flat_vh_Idx] -
        vtf_x_vh_views


    return nothing

end



function v4_flux_decay_model_ode_wt_vtf_no_pf_dynamics!(
    dx,
    x,
    flux_decay_model_dynamics_para,
    t;
    kwd_para =
        model_dynamics_kwd_para )
    
    (;flux_decay_model_states_init,
     counter_array,
     loc_load_exist,
     
     dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
     stateDiffCache,
     by_part_dyn_pf_mismatch_kwd_para,
    
     slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx,

     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx,

     dyn_V_ref_Idx,
     dyn_Tm_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx,
     
     gens_para,
     avrs_para,
     state_vars_idx,
     
     system_ode_para_kwd_para,
     
     H,
     X_d,
     X_q,
     X_d_dash,
     T_d_dash,
     
     Ka,
     Ta,
     
     dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
     non_pre_ordered_pf_vars_Idxs,
     
     gens_vh_idxs,
     gens_θh_idxs,
     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,
     
     gens_id_idxs,
     gens_iq_idxs,     
     flat_vh_flat_θh_flat_id_iq_Idx,
     
     flat_vh_Idx,
     flat_θh_Idx,
     flat_id_Idx,
     flat_iq_Idx,
     
     flux_decay_model_states_Idx,
     flux_decay_model_states_comp_idxs_in_Idx,
     flux_decay_model_vars_Idx_in_state,
     
     δ_idx_in_state,
     ω_idx_in_state,
     eq_dash_idx_in_state,
     E_fd_idx_in_state,
     
     state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,

     ode_vh_id_iq_V_ref_Tm_Idx,
     flux_decay_ode_para_kwd_para) =
        model_dynamics_kwd_para
    
    #----------------------------------------
    
    V_ref = flux_decay_model_dynamics_para[
        dyn_V_ref_Idx]
    
    Tm = flux_decay_model_dynamics_para[
        dyn_Tm_Idx]
    
    Png = flux_decay_model_dynamics_para[
        dyn_Png_Idx]
    
    Qng = flux_decay_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true    
        Pll = flux_decay_model_dynamics_para[
            dyn_Qll_Idx]
        
        Qll = flux_decay_model_dynamics_para[
            dyn_Qll_Idx]        
    else
    
    Pll = []
    Qll = []
        
    end
    
    #----------------------------------------

    vh_θh_Idx_in_state =
        [vh_Idx_in_state;
         θh_Idx_in_state]

    #----------------------------------------

    vec_vh_θh_Idx_in_state =
        [[a_vh_Idx, a_θh_Idx]
         for (a_vh_Idx, a_θh_Idx) in
             zip(vh_Idx_in_state,
                 θh_Idx_in_state) ]

    vec_gens_vh_θh_Idx_in_state =
        vec_vh_θh_Idx_in_state[
            gens_nodes_idx ]

    vec_non_gens_vh_θh_Idx_in_state =
        vec_vh_θh_Idx_in_state[
            non_gens_nodes_idx]

    #----------------------------------------

    ode_dx_views = [view(dx, idx)
                    for idx in state_vars_idx ]

    ode_x_views  = [view( x, idx)
                    for idx in state_vars_idx ]
    
    #----------------------------------------
   
    vtf_dx_gens_vh_θh_views = [
        view( dx, a_vh_θh_Idx )
        for a_vh_θh_Idx in
            vec_gens_vh_θh_Idx_in_state ]
    
    vtf_x_gens_vh_θh_views = [
        view( x, a_vh_θh_Idx )
        for  a_vh_θh_Idx in
           vec_gens_vh_θh_Idx_in_state ]

    #----------------------------------------

    vtf_dx_non_gens_vh_θh_views = [
        view( dx, a_vh_θh_Idx )
        for a_vh_θh_Idx in
            vec_non_gens_vh_θh_Idx_in_state ]
    
    vtf_x_non_gens_vh_θh_views = [
        view( x, a_vh_θh_Idx )
        for  a_vh_θh_Idx in
            vec_non_gens_vh_θh_Idx_in_state ]

    #----------------------------------------
    
    vtf_dx_vh_views = @view dx[vh_Idx_in_state] 

    vtf_x_vh_views  = @view x[vh_Idx_in_state] 

    vtf_dx_θh_views = @view dx[θh_Idx_in_state]
    
    vtf_x_θh_views  = @view x[θh_Idx_in_state]
    
    #----------------------------------------
    # PreallocationTools.jl
    # https://github.com/SciML/PreallocationTools.jl
    #----------------------------------------
        
    counter = counter_array[1]

    #--------------------------------------------

    # stateDiffCache = get_tmp(stateDiffCache, x)

    # if counter == 1
        
    #     stateDiffCache .= flux_decay_model_states_init 

    #     flat_vh = stateDiffCache[ vh_Idx_in_state]
    
    #     flat_θh = stateDiffCache[ θh_Idx_in_state ]

    #     counter_array .+= 1
    # else

    #     stateDiffCache .= x
            
    #     flat_vh = stateDiffCache[ vh_Idx_in_state]
    
    #     flat_θh = stateDiffCache[ θh_Idx_in_state ]

        
    # end

    #----------------------------------------
    
    # pf_δ = stateDiffCache[ δ_idx_in_state]
    
    # pf_eq_dash = stateDiffCache[ eq_dash_idx_in_state]
    
    δ = x[ δ_idx_in_state]
    
    eq_dash = x[ eq_dash_idx_in_state]

    @show  δ
    
    @show eq_dash
    
    #----------------------------------------
            
    flat_vh = x[ vh_Idx_in_state]
    
    flat_θh = x[ θh_Idx_in_state ]

    vh_θh = [flat_vh; flat_θh ]
    
    #----------------------------------------
        
    gens_vh = flat_vh[ gens_nodes_idx ]
    
    gens_θh = flat_θh[ gens_nodes_idx ]
    
    #----------------------------------------

    #----------------------------------------
     
    gens_i_d = (eq_dash - gens_vh .* cos.(
        δ - gens_θh)) ./ X_d_dash

    gens_i_q =
        (gens_vh .* sin.(δ - gens_θh)) ./ X_q

    #----------------------------------------
    
    # gens_Pg_inj =
    #     get_gens_nodes_Pg_net_inj(
    #         δ,
    #         gens_vh,
    #         gens_θh,
    #         gens_i_d,
    #         gens_i_q,
    #         P_g_loc_load,
    #         Q_g_loc_load,

    #         n2s_gens_idx,
    #         gens_nodes_with_loc_loads_idx,

    #         gens_nodes_idx)
    
    # gens_Qg_inj =
    #     get_gens_nodes_Qg_net_inj(
    #         δ,
    #         gens_vh,
    #         gens_θh,
    #         gens_i_d,
    #         gens_i_q,
    #         P_g_loc_load,
    #         Q_g_loc_load,

    #         n2s_gens_idx,
    #         gens_nodes_with_loc_loads_idx,

    #         gens_nodes_idx)

    #--------------------------------------    
    
    # Pg =
    #     get_gens_Pg_from_δ_i_dq_vh_θh(
    #         gens_vh, gens_θh,
    #         δ,
    #         gens_i_d, gens_i_q )

    # Qg =
    #     get_gens_Qg_from_δ_i_dq_vh_θh(
    #         gens_vh, gens_θh,
    #         δ,
    #         gens_i_d, gens_i_q )
    
    #--------------------------------------    
    
    Pg = gens_i_d .* gens_vh .* sin.(δ - gens_θh) +
         gens_i_q .* gens_vh .* cos.(δ - gens_θh)
    
    Qg = gens_i_d .* gens_vh .* cos.(δ - gens_θh) -
        gens_i_q .* gens_vh .* sin.(δ - gens_θh)
    

    #--------------------------------------    
    
    gens_Pg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        Pg[ n2s_gens_idx[ idx]] -
        P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] :
        Pg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]

    gens_Qg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        Qg[ n2s_gens_idx[ idx]] -
        Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] :
        Qg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]
    
    
    dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll =
        [gens_Pg_inj;
         gens_Qg_inj;
         Png;
         Qng;
         Pll;
         Qll]
    
    #----------------------------------------
    
    pf_sol = get_flux_decay_model_pf_ΔPQ_func!(
        vh_θh,
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll,
        by_part_dyn_pf_mismatch_kwd_para =
            by_part_dyn_pf_mismatch_kwd_para )

    # pf_sol = 1.0 .* vh_θh
    
    #----------------------------------------
    #----------------------------------------

    post_pf_vh = pf_sol[ flat_vh_Idx ]

    # post_pf_θh = pf_sol[ flat_θh_Idx ]
    
    pre_norm_post_pf_θh = pf_sol[ flat_θh_Idx ]

    post_pf_θh =
        normalise_angle_θh(
             pre_norm_post_pf_θh,
             slack_gens_nodes_idx )

    @show post_pf_vh

    @show post_pf_θh
    
    post_pf_flat_vh_θh =
        [post_pf_vh; post_pf_θh]
    
    post_pf_vec_vh_θh =
        [ [vh, θh]
          for (vh, θh) in
              zip( post_pf_vh,
                   post_pf_θh)]
     
    #----------------------------------------
    
    post_pf_vec_gens_vh_θh =
        post_pf_vec_vh_θh[
            gens_nodes_idx ]

    post_pf_vec_non_gens_vh_θh =
        post_pf_vec_vh_θh[
            non_gens_nodes_idx ]

    #----------------------------------------
    
    post_pf_gens_vh = post_pf_vh[
        gens_nodes_idx ]
    
    post_pf_gens_θh = post_pf_θh[
        gens_nodes_idx ]

    post_pf_non_gens_vh = post_pf_vh[
        non_gens_nodes_idx ]
    
    post_pf_non_gens_θh = post_pf_θh[
        non_gens_nodes_idx ]

    #----------------------------------------
    
    post_pf_id = ( eq_dash - post_pf_gens_vh .* cos.(
        δ - post_pf_gens_θh )) ./ X_d_dash

    post_pf_iq =
        ( post_pf_gens_vh .* sin.(
            δ - post_pf_gens_θh )) ./ X_q

    #----------------------------------------
    # ode 
    #----------------------------------------
    
    (;gens_para,
     avrs_para,
     ωs) =
         NamedTupleTools.select(
             flux_decay_ode_para_kwd_para,
             (:gens_para,
              :avrs_para,
              :ωs))
    
    (H,
     X_d,
     X_q,
     X_d_dash,
     T_d_dash ) =
         NamedTupleTools.select(
             gens_para,
             (:H,
              :X_d,
              :X_q,
              :X_d_dash,
              :T_d_dash))

    ( Ka,
      Ta ) =
          NamedTupleTools.select(
              avrs_para,
              (:Ka,
               :Ta))
    
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
            zip(post_pf_gens_vh,
            post_pf_id,
            post_pf_iq,
            V_ref,
            Tm) ]

    for (ode_dx, ode_x,
         gen_vh_id_iq_V_ref_Tm_para,
         gen_ode_kwd_para) in
        zip(ode_dx_views,
            ode_x_views,
            gens_vh_id_iq_V_ref_Tm_para,
            gens_ode_kwd_para)
        
        ode_a_gen_flux_decay_by_ext_idq_func!(
            ode_dx,
            ode_x,
            gen_vh_id_iq_V_ref_Tm_para,
            t;
            kwd_para =
                gen_ode_kwd_para )
    end

    # vtf_dx_θh_views
    
    #----------------------------------------
    # all nodes vtf
    #----------------------------------------

    vtf_dx_θh_views .=
        post_pf_flat_vh_θh[flat_θh_Idx] -
        vtf_x_θh_views

    vtf_dx_vh_views .=
        post_pf_flat_vh_θh[flat_vh_Idx] -
        vtf_x_vh_views


    return nothing

end

 
#-----------------------------------------------------
# flux decay model powwer flow 
#-----------------------------------------------------

#-----------------------------------------------------
# begin working dynamic powerflow
#-----------------------------------------------------


function get_flux_decay_model_pf_ΔPQ_func!(
    init_flat_vh_flat_θh,
    dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll;
    by_part_dyn_pf_mismatch_kwd_para =
        by_part_dyn_pf_mismatch_kwd_para )

    #----------------------------------------

    (; pf_solver,
     use_nlsolve,
     dyn_pf_mismatch_kwd_para,
     sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx ) =
         by_part_dyn_pf_mismatch_kwd_para
    
    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------
    
    pf_fun_mismatch =
        get_a_flux_decay_model_pf_ΔPQ_mismatch!
    
    #----------------------------------------
    
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll;
                flux_decay_dyn_pf_mismatch_kwd_para =
                    dyn_pf_mismatch_kwd_para),
                      init_flat_vh_flat_θh,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        return sol
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    flux_decay_dyn_pf_mismatch_kwd_para =
                        dyn_pf_mismatch_kwd_para)),
            init_flat_vh_flat_θh,
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll),
                                    pf_alg )

        return sol
        
    end
        
end


function get_flux_decay_model_pf_ΔPQ_func!(
    init_flat_vh_flat_θh,
    flux_decay_dyn_pf_δ_eq_dash_flat_para,
    flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para;
    by_part_dyn_pf_mismatch_kwd_para =
        by_part_dyn_pf_mismatch_kwd_para )

    #----------------------------------------

    (; pf_solver,
     use_nlsolve,
     dyn_pf_mismatch_kwd_para,
     sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx ) =
         by_part_dyn_pf_mismatch_kwd_para
    
    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------

    (gens_Xd_dash,
     gens_Xq ) = NamedTupleTools.select(
        getproperty(
            getproperty(          
              getproperty(
                  by_part_dyn_pf_mismatch_kwd_para,
                  :dyn_pf_mismatch_kwd_para),
              :dyn_pf_mismatch_vars_kwd_para),
            :ode_gens_para),
        (:X_d_dash,
         :X_q))
    
    #----------------------------------------

    (;
     gens_vh_idxs,
     gens_θh_idxs
     # ,
     # non_gens_nodes_vh_idxs,
     # non_gens_nodes_θh_idxs,
     # flat_vh_flat_θh_flat_id_iq_Idx
     ) =  NamedTupleTools.select(
             getproperty(
                 getproperty(
                     getproperty(
                         by_part_dyn_pf_mismatch_kwd_para,
                         :dyn_pf_mismatch_kwd_para),
                     :dyn_pf_mismatch_idx_kwd_para),
             :non_pre_ordered_pf_vars_Idxs) ,
             (:gens_vh_idxs,
              :gens_θh_idxs
              # ,              
              # :non_gens_nodes_vh_idxs,
              # :non_gens_nodes_θh_idxs,
              # :flat_vh_flat_θh_flat_id_iq_Idx
              ))
    
    #-------------------------------
    
    # (;flat_vh_Idx,
    #   flat_θh_Idx) =
    #      NamedTupleTools.select(
    #          flat_vh_flat_θh_flat_id_iq_Idx,
    #          (:flat_vh_Idx,
    #           :flat_θh_Idx))

    #----------------------------------------
    
    (;dyn_pf_δ_Idxs,
     dyn_pf_eq_dash_Idxs) =
        NamedTupleTools.select(
            getproperty(
                sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx,
                :dyn_pf_δ_eq_dash_Idx),
            (:dyn_pf_δ_Idxs,
             :dyn_pf_eq_dash_Idxs))

    #----------------------------------------
    

    # (; dyn_pf_δ_Idxs,
    #  dyn_pf_eq_dash_Idxs,

    #  dyn_pf_P_non_gens_Idxs,
    #  dyn_pf_Q_non_gens_Idxs,

    #  dyn_pf_P_gens_loc_load_Idxs,
    #  dyn_pf_Q_gens_loc_load_Idxs ) =
    #      NamedTupleTools.select(
    #                       getproperty(
    #              getproperty(
    #                  getproperty(
    #                      by_part_dyn_pf_mismatch_kwd_para,
    #                      :dyn_pf_mismatch_kwd_para),
    #                  :dyn_pf_mismatch_idx_kwd_para),
    #              :dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx),
    
    #          (:dyn_pf_δ_Idxs,
    #           :dyn_pf_eq_dash_Idxs,

    #           :dyn_pf_P_non_gens_Idxs,
    #           :dyn_pf_Q_non_gens_Idxs,

    #           :dyn_pf_P_gens_loc_load_Idxs,
    #           :dyn_pf_Q_gens_loc_load_Idxs))
    
    #----------------------------------------

    (slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx ) =
         NamedTupleTools.select(
             getproperty(
                 getproperty(
                     getproperty(
                         by_part_dyn_pf_mismatch_kwd_para,
                         :dyn_pf_mismatch_kwd_para),
                     :dyn_pf_mismatch_idx_kwd_para),
                 :dyn_pf_fun_kwd_net_idxs),
             (:slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx )) 

    #----------------------------------------    


    (n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             getproperty(
                 getproperty(
                     getproperty(
                         by_part_dyn_pf_mismatch_kwd_para,
                         :dyn_pf_mismatch_kwd_para),
                     :dyn_pf_mismatch_idx_kwd_para),
                 :dyn_pf_fun_kwd_n2s_idxs),
             (:n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx)) 

    #----------------------------------------
    #----------------------------------------
    
    gens_vh = init_flat_vh_flat_θh[
        gens_vh_idxs]

    gens_θh = init_flat_vh_flat_θh[
        gens_θh_idxs]

    #----------------------------------------

    gens_δ =
        flux_decay_dyn_pf_δ_eq_dash_flat_para[
            dyn_pf_δ_Idxs]
    
    gens_eq_dash =
        flux_decay_dyn_pf_δ_eq_dash_flat_para[
            dyn_pf_eq_dash_Idxs]

    #----------------------------------------
    
    gens_i_d = (gens_eq_dash - gens_vh .* cos.(
        gens_δ - gens_θh )) ./ gens_Xd_dash

    gens_i_q = gens_vh .* sin.(gens_δ - gens_θh) ./
        gens_Xq

    #----------------------------------------

    gens_Pg = get_gens_Pg_from_δ_i_dq_vh_θh(
        gens_vh, gens_θh,
        gens_δ,
        gens_i_d, gens_i_q )

    gens_Qg = get_gens_Qg_from_δ_i_dq_vh_θh(
        gens_vh, gens_θh,
        gens_δ,
        gens_i_d, gens_i_q )

    #----------------------------------------    

    Pg_net_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        gens_Pg[ n2s_gens_idx[ idx]] -
        P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx] ] :
        gens_Pg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]

    Qg_net_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        gens_Qg[ n2s_gens_idx[ idx]] -
        Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx] ] :
        gens_Qg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]
    
    dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll =
        [Pg_net_inj;
         Qg_net_inj;
         flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para ]
    
    #----------------------------------------    
    #----------------------------------------
    
    pf_fun_mismatch =
        get_a_flux_decay_model_pf_ΔPQ_mismatch!
    
    #----------------------------------------
    
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll;
                flux_decay_dyn_pf_mismatch_kwd_para =
                    dyn_pf_mismatch_kwd_para),
                      init_flat_vh_flat_θh,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        return sol
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    flux_decay_dyn_pf_mismatch_kwd_para =
                        dyn_pf_mismatch_kwd_para)),
            init_flat_vh_flat_θh,
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll),
                                    pf_alg )

        return sol
        
    end
        
end


function get_flux_decay_model_pf_ΔI_func!(
    init_flat_vh_flat_θh,
    flux_decay_dyn_pf_δ_eq_dash_flat_para,
    flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para;
    by_part_dyn_pf_mismatch_kwd_para =
        by_part_dyn_pf_mismatch_kwd_para )

    #----------------------------------------

    (; pf_solver,
     use_nlsolve,
     dyn_pf_mismatch_kwd_para,
     sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx ) =
         by_part_dyn_pf_mismatch_kwd_para
    
    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------

    (gens_Xd_dash,
     gens_Xq ) = NamedTupleTools.select(
        getproperty(
            getproperty(          
              getproperty(
                  by_part_dyn_pf_mismatch_kwd_para,
                  :dyn_pf_mismatch_kwd_para),
              :dyn_pf_mismatch_vars_kwd_para),
            :ode_gens_para),
        (:X_d_dash,
         :X_q))
    
    #----------------------------------------

    (;
     gens_vh_idxs,
     gens_θh_idxs ) =
         NamedTupleTools.select(
             getproperty(
                 getproperty(
                     getproperty(
                         by_part_dyn_pf_mismatch_kwd_para,
                         :dyn_pf_mismatch_kwd_para),
                     :dyn_pf_mismatch_idx_kwd_para),
             :non_pre_ordered_pf_vars_Idxs) ,
             (:gens_vh_idxs,
              :gens_θh_idxs ))
    
    #----------------------------------------
    
    (;dyn_pf_δ_Idxs,
     dyn_pf_eq_dash_Idxs) =
        NamedTupleTools.select(
            getproperty(
                sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx,
                :dyn_pf_δ_eq_dash_Idx),
            (:dyn_pf_δ_Idxs,
             :dyn_pf_eq_dash_Idxs))

    #----------------------------------------
    #----------------------------------------

    (slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx ) =
         NamedTupleTools.select(
             getproperty(
                 getproperty(
                     getproperty(
                         by_part_dyn_pf_mismatch_kwd_para,
                         :dyn_pf_mismatch_kwd_para),
                     :dyn_pf_mismatch_idx_kwd_para),
                 :dyn_pf_fun_kwd_net_idxs),
             (:slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx )) 

    #----------------------------------------    


    (n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             getproperty(
                 getproperty(
                     getproperty(
                         by_part_dyn_pf_mismatch_kwd_para,
                         :dyn_pf_mismatch_kwd_para),
                     :dyn_pf_mismatch_idx_kwd_para),
                 :dyn_pf_fun_kwd_n2s_idxs),
             (:n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx)) 

    #----------------------------------------
    #----------------------------------------
    
    gens_vh = init_flat_vh_flat_θh[
        gens_vh_idxs]

    gens_θh = init_flat_vh_flat_θh[
        gens_θh_idxs]

    #----------------------------------------

    gens_δ =
        flux_decay_dyn_pf_δ_eq_dash_flat_para[
            dyn_pf_δ_Idxs]
    
    gens_eq_dash =
        flux_decay_dyn_pf_δ_eq_dash_flat_para[
            dyn_pf_eq_dash_Idxs]

    #----------------------------------------
    
    gens_i_d = (gens_eq_dash - gens_vh .* cos.(
        gens_δ - gens_θh ) ) ./ gens_Xd_dash

    gens_i_q = gens_vh .* sin.(gens_δ - gens_θh) ./
        gens_Xq

    #----------------------------------------
    #----------------------------------------    
       
    gens_I_D_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        (gens_i_q[ n2s_gens_idx[idx]] * cos( dyn_δ[ n2s_gens_idx[idx]]) +
         gens_i_d[ n2s_gens_idx[idx]] * sin( dyn_δ[ n2s_gens_idx[idx]]) -
          (P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] * cos(gens_θh[ n2s_gens_idx[ idx]]) +
           Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] * sin(gens_θh[ n2s_gens_idx[ idx]]))/
           gens_vh[ n2s_gens_idx[idx]]) :
          (gens_i_q[ n2s_gens_idx[idx]] * cos( dyn_δ[ n2s_gens_idx[idx]]) +
           gens_i_d[ n2s_gens_idx[idx]] * sin( dyn_δ[ n2s_gens_idx[idx]])) 
                  for idx in gens_nodes_idx ]

    gens_I_Q_inj = [ idx ∈ gens_nodes_with_loc_loads_idx ?
        (gens_i_q[ n2s_gens_idx[idx]] * sin( dyn_δ[ n2s_gens_idx[idx]]) -
         gens_i_d[ n2s_gens_idx[idx]] * cos( dyn_δ[ n2s_gens_idx[idx]]) -
         (P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] *
         sin(gens_θh[ n2s_gens_idx[ idx]]) -
       Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] * cos(gens_θh[ n2s_gens_idx[ idx]]))/
       gens_vh[ n2s_gens_idx[idx]]) :
       (gens_i_q[ n2s_gens_idx[idx]] * sin( dyn_δ[ n2s_gens_idx[idx]]) -
        gens_i_d[ n2s_gens_idx[idx]] * cos( dyn_δ[ n2s_gens_idx[idx]])) 
                  for idx in gens_nodes_idx ]
    
    dyn_pf_flat_I_D_inj_I_Q_inj_Png_Qng_Pgll_Qgll =
        [gens_I_D_inj;
         gens_I_Q_inj;
         flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para ]
    
    #----------------------------------------    
    #----------------------------------------
    
    pf_fun_mismatch =
        get_a_flux_decay_model_pf_ΔI_mismatch!
    
    #----------------------------------------
    
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            dyn_pf_flat_I_D_inj_I_Q_inj_Png_Qng_Pgll_Qgll;
                flux_decay_dyn_pf_mismatch_kwd_para =
                    dyn_pf_mismatch_kwd_para),
                      init_flat_vh_flat_θh,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        return sol
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    flux_decay_dyn_pf_mismatch_kwd_para =
                        dyn_pf_mismatch_kwd_para)),
            init_flat_vh_flat_θh,
            dyn_pf_flat_I_D_inj_I_Q_inj_Png_Qng_Pgll_Qgll),
                                    pf_alg )

        return sol
        
    end
        
end


function get_a_flux_decay_model_pf_ΔPQ_mismatch!(
    ΔPQ,
    vh_θh,
    dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll;
    flux_decay_dyn_pf_mismatch_kwd_para =
        dyn_pf_mismatch_kwd_para  )
    
    # -------------------------------------

    (; loc_load_exist,
     dyn_pf_mismatch_vars_kwd_para,
     dyn_pf_mismatch_idx_kwd_para ) =
         NamedTupleTools.select(
             flux_decay_dyn_pf_mismatch_kwd_para,
             (:loc_load_exist,
              :dyn_pf_mismatch_vars_kwd_para,
              :dyn_pf_mismatch_idx_kwd_para)) 
    
    # -------------------------------------

    (Ynet_wt_nodes_idx_wt_adjacent_nodes,
      ode_gens_para) =
          NamedTupleTools.select(
              dyn_pf_mismatch_vars_kwd_para,
              (:Ynet_wt_nodes_idx_wt_adjacent_nodes,
               :ode_gens_para))

    # -------------------------------------

   ( Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx)) 
    
    # -------------------------------------

    (;
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs ) =
         NamedTupleTools.select(
             dyn_pf_mismatch_idx_kwd_para,
             (
              :Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              :non_pre_ordered_pf_vars_Idxs ) )
    
    #-----------------

   (;Pg_Idxs,
    Qg_Idxs,
    Png_Idxs,
    Qng_Idxs,
    Pgll_Idxs,
    Qgll_Idxs ) =
        NamedTupleTools.select(
            Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
            (:Pg_Idxs,
             :Qg_Idxs,
             :Png_Idxs,
             :Qng_Idxs,
             :Pgll_Idxs,
             :Qgll_Idxs))

    #-----------------

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
    
    #-----------------

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
    
    #-----------------

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
    
    #-------------------------------
    
    (;flat_vh_Idx,
      flat_θh_Idx) =
         NamedTupleTools.select(
             flat_vh_flat_θh_flat_id_iq_Idx,
             (:flat_vh_Idx,
              :flat_θh_Idx))
    
    #-------------------------------
    #-------------------------------

    Pg_net_inj =
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
            Pg_Idxs]

    Qg_net_inj =
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
            Qg_Idxs]
    
    P_non_gens =
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
            Png_Idxs]

    Q_non_gens =
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
            Qng_Idxs]

    
    if loc_load_exist == true

        P_g_loc_load =
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
                Pgll_Idxs]

        Q_g_loc_load =
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
                Pgll_Idxs]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #-------------------------------
    

    # Pg_net_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
    #     P_gens[ n2s_gens_idx[ idx]] -
    #     P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx] ] :
    #     P_gens[ n2s_gens_idx[ idx]]
    #               for idx in gens_nodes_idx ]


    # Qg_net_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
    #     Q_gens[ n2s_gens_idx[ idx]] -
    #     Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx] ] :
    #     Q_gens[ n2s_gens_idx[ idx] ]
    #               for idx in gens_nodes_idx ]

    
    #-------------------------------
    
    flat_vh =
        vh_θh[
            flat_vh_Idx ]

    flat_θh =
        vh_θh[
            flat_θh_Idx ]

    #-------------------------------
    
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

    
    # -----------------------------------
    # active power mismatch
    # -----------------------------------
    
    """
    Gens nodes
    Pg_i_inj =
        vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    Load nodes
    Pl_i =
        -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    
    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
    (P_non_gens[ n2s_non_gens_idx[ nth_idx ]] + vh[ n2s_all_nodes_idx[nth_idx]] *
            sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                    cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
                for (ynj, idx) in
                    zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
                        nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx] ])])) :     
     (Pg_net_inj[ n2s_gens_idx[ nth_idx] ] -
     vh[n2s_all_nodes_idx[nth_idx]] *
     sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
              for (ynj, idx) in zip(
                  Ynet[ n2s_all_nodes_idx[nth_idx] ],
                  nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) )
    for nth_idx in all_nodes_idx ]

    
    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------
    
    """
    Gens nodes
    Qg_i_inj =
        vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    Load nodes
    Ql_i =
        -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """
    
    Q_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
      (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
      vh[ n2s_all_nodes_idx[nth_idx]] *
      sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
      sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
                for (ynj, idx) in
                    zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
                        nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :    
     (Qg_net_inj[ n2s_gens_idx[ nth_idx] ] -
      vh[ n2s_all_nodes_idx[ nth_idx]] *
      sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
      sin( θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj)) 
                    for (ynj, idx) in
                        zip(Ynet[ n2s_all_nodes_idx[ nth_idx]],
                           nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])] ) )
    for nth_idx in all_nodes_idx ]
    
    # ------------------------------------

    ΔPQ .=
        vcat(P_mismatch,
             Q_mismatch)

    return nothing

end

#-----------------------------------------------------

function get_a_flux_decay_model_pf_ΔI_mismatch!(
    ΔI, vh_θh,
    flux_decay_dyn_pf_δ_ed_dash_non_gen_PQ;
    flux_decay_dyn_pf_mismatch_kwd_para =
        dyn_pf_mismatch_kwd_para  )
        
    # -------------------------------------

    (;loc_load_exist,
     dyn_pf_mismatch_vars_kwd_para,
     dyn_pf_mismatch_idx_kwd_para ) =
         NamedTupleTools.select(
             flux_decay_dyn_pf_mismatch_kwd_para,
             ( :loc_load_exist,
               :dyn_pf_mismatch_vars_kwd_para,
               :dyn_pf_mismatch_idx_kwd_para)) 
    
    # -------------------------------------

    (;Ynet_wt_nodes_idx_wt_adjacent_nodes,
      ode_gens_para) =
          NamedTupleTools.select(
              dyn_pf_mismatch_vars_kwd_para,
              (:Ynet_wt_nodes_idx_wt_adjacent_nodes,
               :ode_gens_para))

    #-----------------

   ( Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx))

    #-----------------
    
    (gens_Xd_dash,
     gens_Xq ) =
         NamedTupleTools.select(
             ode_gens_para,
             (:X_d_dash,
              :X_q ) )
    
    # -------------------------------------

    (; dyn_pf_δ_eq_dash_0_P_Q_idx,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs,
     dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,
     dyn_pf_idq_Idx,
     dyn_pf_idq_δ_ed_dash_non_gen_PQ_Idx     
     ) =
         NamedTupleTools.select(
             dyn_pf_mismatch_idx_kwd_para,
             (:dyn_pf_δ_eq_dash_0_P_Q_idx,
              :Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              :non_pre_ordered_pf_vars_Idxs,
              :dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,
              :dyn_pf_idq_Idx,
              :dyn_pf_idq_δ_ed_dash_non_gen_PQ_Idx ) )

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
    

    (;
     gens_vh_idxs,
     gens_θh_idxs,

     non_gens_nodes_vh_idxs,
     non_gens_nodes_θh_idxs,

     gens_id_idxs,
     gens_iq_idxs,
     
     flat_vh_flat_θh_flat_id_iq_Idx) =
         NamedTupleTools.select(
             non_pre_ordered_pf_vars_Idxs,
             (:gens_vh_idxs,
              :gens_θh_idxs,

              :non_gens_nodes_vh_idxs,
              :non_gens_nodes_θh_idxs,

              :gens_id_idxs,
              :gens_iq_idxs,

              :flat_vh_flat_θh_flat_id_iq_Idx)) 
    
    #-------------------------------
    
    (flat_vh_Idx,
     flat_θh_Idx) =
         NamedTupleTools.select(
             flat_vh_flat_θh_flat_id_iq_Idx,
             (:flat_vh_Idx,
              :flat_θh_Idx)) 

    #-------------------------------

    (; dyn_pf_δ_Idxs,
     dyn_pf_eq_dash_Idxs,

     dyn_pf_P_non_gens_Idxs,
     dyn_pf_Q_non_gens_Idxs,

     dyn_pf_P_gens_loc_load_Idxs,
     dyn_pf_Q_gens_loc_load_Idxs ) =
         NamedTupleTools.select(
             dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
             (:dyn_pf_δ_Idxs,
              :dyn_pf_eq_dash_Idxs,

              :dyn_pf_P_non_gens_Idxs,
              :dyn_pf_Q_non_gens_Idxs,

              :dyn_pf_P_gens_loc_load_Idxs,
              :dyn_pf_Q_gens_loc_load_Idxs))
    
    #-------------------------------
    
    gens_I_D_inj =
        flux_decay_dyn_pf_δ_ed_dash_non_gen_PQ[
            dyn_pf_δ_Idxs]
    
    gens_I_Q_inj =
        flux_decay_dyn_pf_δ_ed_dash_non_gen_PQ[
            dyn_pf_eq_dash_Idxs]
    
    P_non_gens =
        flux_decay_dyn_pf_δ_ed_dash_non_gen_PQ[
            dyn_pf_P_non_gens_Idxs]
    
    Q_non_gens =
        flux_decay_dyn_pf_δ_ed_dash_non_gen_PQ[
            dyn_pf_Q_non_gens_Idxs]

    
    if loc_load_exist == true
    
        P_g_loc_load =
            flux_decay_dyn_pf_δ_ed_dash_non_gen_PQ[
                dyn_pf_P_gens_loc_load_Idxs]
    
        Q_g_loc_load =
            flux_decay_dyn_pf_δ_ed_dash_non_gen_PQ[
                dyn_pf_Q_gens_loc_load_Idxs]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #-------------------------------
    
    flat_vh =
        vh_θh[
            flat_vh_Idx ]

    flat_θh =
        vh_θh[
            flat_θh_Idx ]
    
    #-------------------------------
    
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
    
    # -----------------------------------
    # -----------------------------------

    """

    (id_i  + im * iq_i) * exp( im * ( δ_i - π/2 )) +
    (Ploc_i - Qloc_i) / ( vh_i exp( -im * θ_i )) =
     ∑( vh_k * Y_ik * exp( im * ( θ_k + β_ik )) )

    for non gen nodes, id_i = iq_i = 0
    
    """

    # -----------------------------------
    # nodes  current real part mismatch
    # -----------------------------------
    
    """

    iq_i * cos(δ_i) + id_i * sin(δ_i) + 
    (Ploc_i * cos(θ_i) + Qloc_i * sin(θ_i))/vh_i =
        ∑( vh_k * Y_ik * cos(θ_k + β_ik) )

    # -----------------------------------
    
   """

    I_real_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
        ( (P_non_gens[ n2s_non_gens_idx[ nth_idx]] * cos(θh[ n2s_all_nodes_idx[ nth_idx]]) +
           Q_non_gens[ n2s_non_gens_idx[ nth_idx]] * sin(θh[ n2s_all_nodes_idx[ nth_idx]]))/
        vh[ n2s_all_nodes_idx[ nth_idx]] +
            sum( [vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * cos(θh[n2s_all_nodes_idx[idx]] + angle(ynj))
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
         ( gens_I_D_inj[ n2s_gens_idx[nth_idx]] -
           sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * cos(θh[ n2s_all_nodes_idx[idx]] + angle(ynj))
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx] ]) ] ) )
                        for nth_idx in all_nodes_idx ]

    
    # -----------------------------------
    # nodes  current imag part mismatch
    # -----------------------------------

    """

     iq_i * sin(δ_i) - id_i * cos(δ_i) +
     (Ploc_i * sin(θ_i) - Qloc_i * cos(θ_i)) / vh_i =
        ∑( vh_k * Y_ik * sin( θ_k + β_ik ) )

    # -----------------------------------

    """
    I_imag_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
      ( (P_non_gens[ n2s_non_gens_idx[ nth_idx ]] * sin( θh[ n2s_all_nodes_idx[nth_idx]]) -
         Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] * cos( θh[ n2s_all_nodes_idx[nth_idx]])) /
        vh[ n2s_all_nodes_idx[nth_idx]] +
        sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * sin(θh[ n2s_all_nodes_idx[idx]] + angle(ynj))
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
        (gens_I_Q_inj[ n2s_gens_idx[nth_idx]] -
         sum( [vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * sin( θh[ n2s_all_nodes_idx[idx]] + angle(ynj))
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) 
                              for nth_idx in all_nodes_idx ]

    # ------------------------------------
        
    ΔI .=
        vcat(I_real_mismatch,
             I_imag_mismatch )

    return nothing

end


# ---------------------------------------------------
# end working dynamic powerflow
# ---------------------------------------------------


# ---------------------------------------------------
#####################################################
# ---------------------------------------------------
# structure preserving classical model
# ---------------------------------------------------
#####################################################
# ---------------------------------------------------


function get_a_gen_spcm_Pg(
    δ, E,
    vh, θh,
    X_d_dash )

    Pg = E * vh * sin(
         δ - θh) / X_d_dash
    
    return Pg
end



function get_a_gen_spcm_Qg(
    δ, E,
    vh, θh,
    Xd_dash )
    
    Qg = -vh^2 / Xd_dash + E * vh * cos(
        δ - θh ) / Xd_dash

    return Qg
end

#-----------------------------------------------------

function spcm_Pg(
    gens_δ, gens_E,
    gens_vh, gens_θh,
    gens_Xd_dash )

    Pg = gens_E .* gens_vh .*
        sin.(gens_δ - gens_θh ) ./ gens_Xd_dash
    
    return Pg
end


function spcm_Qg(
    gens_δ, gens_E,
    gens_vh, gens_θh,
    gens_Xd_dash )
    
    Qg = -gens_vh.^2 ./gens_Xd_dash .+
        gens_E .* gens_vh .*
        cos.( gens_δ - gens_θh) ./ gens_Xd_dash

    return Qg
end


#-----------------------------------------------------
# structure preserving classical model state init
#-----------------------------------------------------


# get_init_gens_spcm =
#     get_init_gens_classical_model

# get_state_init_spcm =
#     get_state_init_classical_model



#-----------------------------------------------------
# structure preserving classical model ode
#-----------------------------------------------------

function a_gen_spcm_swing!(
    dx, x, p, t;
    kwd_para = kwd_para )

    vh, θh, E, Tm = p

    (H, Xd_dash, ωs) = kwd_para

    δ, ω = x

    dx[1] = ω - ωs
    
    dx[2] = ( ωs / ( 2 * H )) *
        (Tm - E * vh * sin( δ - θh ) / Xd_dash)

    return nothing
    
end


function ode_spcm_swing!(
    dx, x, p, t;
    kwd_para = kwd_para )

    ( gens_vh,
      gens_θh,
      gens_E,
      gens_Tm ) = p

    (gens_para,
     ωs,
     state_vars_idx) =
          kwd_para
    
    ( gens_H,
      gens_X_d_dash ) =
          gens_para

    dx_views =
        [view(dx, idx)
                for idx in state_vars_idx]

    x_views =
        [view(x, idx)
                for idx in state_vars_idx]
    
    for (a_dx_view, a_x_view, vh, θh,
         E, Tm,
         H, Xd_dash,
         δ_idx_i, ω_idx_i ) in
        zip(dx_views, x_views, gens_vh, gens_θh,
            gens_E, gens_Tm,
            gens_H, gens_Xd_dash,
            δ_idx, ω_idx )

        a_gen_spcm_para = (vh, θh, E, Tm)
        
        a_gen_spcm_kwd_para = (H, Xd_dash, ωs)
        
        a_gen_spcm_swing!(
            # dx[state_var_idx],
            # x[state_var_idx],
            a_dx_view,
            a_x_view,
            a_gen_spcm_para,
            t;
            kwd_para =
                a_gen_spcm_kwd_para )

    end

    return nothing
    
end


#-----------------------------------------------------
# structure preserving classical model powerflow
#-----------------------------------------------------


function get_spcm_pf_ΔPQ_func!(
    init_flat_vh_flat_θh,
    dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll;
    by_part_dyn_pf_mismatch_kwd_para =
        by_part_dyn_pf_mismatch_kwd_para )

    #----------------------------------------

    (; pf_solver,
     use_nlsolve,
     dyn_pf_mismatch_kwd_para,
     sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx ) =
         by_part_dyn_pf_mismatch_kwd_para
    
    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------
    
    pf_fun_mismatch =
        get_spcm_pf_ΔPQ_mismatch
    
    #----------------------------------------
    
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll;
                mismatch_kwd_para =
                    dyn_pf_mismatch_kwd_para),
                      init_flat_vh_flat_θh,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        return sol
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    mismatch_kwd_para =
                        dyn_pf_mismatch_kwd_para)),
            init_flat_vh_flat_θh,
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll),
                                    pf_alg )

        return sol
        
    end
        
end



function get_spcm_pf_ΔPQ_mismatch(
    ΔPQ,
    vh_θh,
    dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll;
    mismatch_kwd_para =
        dyn_pf_mismatch_kwd_para)
    
    
    # -------------------------------------

    (;loc_load_exist,
     dyn_pf_mismatch_vars_kwd_para,
     dyn_pf_mismatch_idx_kwd_para ) =
         NamedTupleTools.select(
             mismatch_kwd_para,
             (:loc_load_exist,
              :dyn_pf_mismatch_vars_kwd_para,
              :dyn_pf_mismatch_idx_kwd_para)) 
    
    # -------------------------------------

    (Ynet_wt_nodes_idx_wt_adjacent_nodes,
      ode_gens_para) =
          NamedTupleTools.select(
              dyn_pf_mismatch_vars_kwd_para,
              (:Ynet_wt_nodes_idx_wt_adjacent_nodes,
               :ode_gens_para))

    # -------------------------------------

   (Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx)) 
    
    # -------------------------------------

    (;
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs ) =
         NamedTupleTools.select(
             dyn_pf_mismatch_idx_kwd_para,
             (
              :Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              :non_pre_ordered_pf_vars_Idxs ) )
    
    #-----------------

   (;Pg_Idxs,
    Qg_Idxs,
    Png_Idxs,
    Qng_Idxs,
    Pgll_Idxs,
    Qgll_Idxs ) =
        NamedTupleTools.select(
            Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
            (:Pg_Idxs,
             :Qg_Idxs,
             :Png_Idxs,
             :Qng_Idxs,
             :Pgll_Idxs,
             :Qgll_Idxs))

    #-----------------

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
    
    #-----------------

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
    
    #-----------------

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
    
    #-------------------------------
    
    (;flat_vh_Idx,
      flat_θh_Idx) =
         NamedTupleTools.select(
             flat_vh_flat_θh_flat_id_iq_Idx,
             (:flat_vh_Idx,
              :flat_θh_Idx))
    
    #-------------------------------
    #-------------------------------

    Pg_net_inj =
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
            Pg_Idxs]

    Qg_net_inj =
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
            Qg_Idxs]
    
    P_non_gens =
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
            Png_Idxs]

    Q_non_gens =
        dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
            Qng_Idxs]

    
    if loc_load_exist == true

        P_g_loc_load =
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
                Pgll_Idxs]

        Q_g_loc_load =
            dyn_pf_flat_Pg_inj_Qg_inj_Png_Qng_Pgll_Qgll[
                Pgll_Idxs]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #-------------------------------
    
    flat_vh =
        vh_θh[
            flat_vh_Idx ]

    flat_θh =
        vh_θh[
            flat_θh_Idx ]

    #-------------------------------
    
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

    
    # -----------------------------------
    # active power mismatch
    # -----------------------------------
    
    """
    Gens nodes
    Pg_i_inj =
        vh_i * ∑(vh_k * B_ik * sin(θ_i - θ_k ) )

    Load nodes
    Pl_i =
        -vh_i * ∑(vh_k * B_ik * sin(θ_i - θ_k ) )

    """
    
    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
        (P_non_gens[ n2s_non_gens_idx[ nth_idx ]] + vh[
            n2s_all_nodes_idx[nth_idx]] *
            sum([ vh[ n2s_all_nodes_idx[idx]] * imag(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[
                n2s_all_nodes_idx[idx]] )
                for (ynj, idx) in
                    zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
                        nodes_idx_with_adjacent_nodes_idx[
                            n2s_all_nodes_idx[nth_idx] ])])) :     
     (Pg_net_inj[ n2s_gens_idx[ nth_idx] ] -
     vh[n2s_all_nodes_idx[nth_idx]] *
     sum([ vh[ n2s_all_nodes_idx[idx]] * imag(ynj) *
     sin( θh[ n2s_all_nodes_idx[nth_idx]] -
     θh[ n2s_all_nodes_idx[idx]] )
              for (ynj, idx) in zip(
                  Ynet[ n2s_all_nodes_idx[nth_idx] ],
                  nodes_idx_with_adjacent_nodes_idx[
                      n2s_all_nodes_idx[nth_idx]])] ) )
    for nth_idx in all_nodes_idx ]

    
    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------
    
    """
    Gens nodes
    Qg_i_inj =
        -vh_i * ∑(vh_k * B_ik * cos(θ_i - θ_k) )

    Load nodes
    Ql_i =
        vh_i * ∑(vh_k * B_ik * cos(θ_i - θ_k ) )

    """
    
    Q_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
      (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] -
      vh[ n2s_all_nodes_idx[nth_idx]] *
      sum([ vh[ n2s_all_nodes_idx[idx]] * imag(ynj) *
      cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[
          n2s_all_nodes_idx[idx]] )
                for (ynj, idx) in
                    zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
                        nodes_idx_with_adjacent_nodes_idx[
                            n2s_all_nodes_idx[nth_idx]])] ) ) :    
     (Qg_net_inj[ n2s_gens_idx[ nth_idx] ] +
      vh[ n2s_all_nodes_idx[ nth_idx]] *
      sum([ vh[ n2s_all_nodes_idx[ idx]] * imag(ynj) *
      cos( θh[ n2s_all_nodes_idx[ nth_idx]] - θh[
          n2s_all_nodes_idx[idx]] ) 
                    for (ynj, idx) in
                        zip(Ynet[ n2s_all_nodes_idx[ nth_idx]],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[ nth_idx]])] ) )
    for nth_idx in all_nodes_idx ]
    
    # ------------------------------------

    ΔPQ .=
        vcat(P_mismatch,
             Q_mismatch)

    return nothing

end


# ---------------------------------------------------
#####################################################
# ---------------------------------------------------
# internal node model 
# ---------------------------------------------------
#####################################################
# ---------------------------------------------------

# ---------------------------------------------------
# internal node model state init
# ---------------------------------------------------

# get_init_gens_gen_inm =
#     get_init_gens_classical_model

# get_state_init_inm =
#     get_state_init_classical_model

#-----------------------------------------------------


#-----------------------------------------------------
# internal node model
#-----------------------------------------------------

# https://discourse.julialang.org/t/vector-valued-function-in-multidimensional-list-comprehension-differentiable-code/68772

# https://discourse.julialang.org/t/vector-valued-function-in-multidimensional-list-comprehension-differentiable-code/68772


function get_a_gen_Ig(
    gen_vh, gen_θh, Pg, Qg  )

    return ( Pg - im * Qg ) / (
        gen_vh * exp(-im * gen_θh) )

end


function get_a_gen_Eg_wt_δg(
    gen_vh, gen_θh,
    Pg, Qg,
    X_d_dash )

    E_g = gen_vh * exp(im * gen_θh) +
        im * X_d_dash * get_a_gen_Ig(
            gen_vh, gen_θh, Pg, Qg  )

    return (Eg = abs(E_g), δg = angle(E_g) )

end


function get_Cinm_Dinm(
    Eg, Yint)

    rows_size, cols_size = size(Yint)

    Cinm = reshape( [ i==j ? 0.0 :
        Eg[i] * Eg[j] * imag(Yint[i, j])
              for i in 1:cols_size, j in 1:cols_size  ],
            (rows_size, cols_size) ) 


    Dinm = reshape([ Eg[i] * Eg[j] * real(Yint[i, j])
              for i in 1:rows_size, j in 1:cols_size ],
                   (rows_size, cols_size) )

    return (;Cinm, Dinm)
    

end


function get_minus_∂Pei_∂δj(
    i, j, δg, Cinm, Dinm)

    rows_size, cols_size = size(Cinm)

    return  i != j ?
        C[i,j] * cos(δg[i] - δg[j]) - D[i,j] *
        sin(δg[i] - δg[j]) :
        -sum([C[i,m] * cos(δg[i] - δg[m]) -
        D[i,m] * sin(δg[i] - δg[m])
          for m in 1:cols_size if m != i])
end



function get_minus_∂Pei_∂δj(
    i, j, δg, Eg, Yint, nothing)

    rows_size, cols_size = size(Yint)

    return  i != j ? Eg[i] * Eg[j] * abs(Yint[i,j]) *
        sin(angle(Yint[i,j]) - (δg[i] - δg[j])) :
        -sum( [Eg[i] * Eg[m] * abs(Yint[i,m]) *
        sin( angle(Yint[i,m]) - ( δg[i] - δg[m]))
              for m in 1:cols_size if m != i ] )
end

            
function get_minus_∂Pei_∂δj_matrix(
    δg,
    Cinm,
    Dinm)

    rows_size, cols_size = size(Cinm)

    return reshape(
        [ get_minus_∂Pei_∂δj(i, j, δg, Cinm, Dinm)
          for i in 1:rows_size, j in
              1:cols_size ],
                   (rows_size, cols_size))

end

            

function get_minus_∂Pei_∂δj_matrix(
    δg,
    Eg,
    Yint,
    nothing)

    rows_size, cols_size = size(Yint)

    return reshape(
        [ get_minus_∂Pei_∂δj(i, j, δg, Eg, Yint, nothing)
          for i in 1:rows_size, j in
              1:cols_size ],
                   (rows_size, cols_size))
end


function get_Aω_matrix(
    δg,
    H,    
    Cinm,
    Dinm,
    ωs)

    minus_∂Pei_∂δj_matrix =
        get_minus_∂Pei_∂δj_matrix(
            δg, Cinm, Dinm)

    M_inv = (ωs/2) .* Diagonal(inv.(H))

    Aω_matrix = M_inv * minus_∂Pei_∂δj_matrix
    

    return Aω_matrix

end


function get_Aω_matrix(
    δg,
    Eg,
    H,    
    Yint,
    ωs,
    nothing)

    minus_∂Pei_∂δj_matrix =
        get_minus_∂Pei_∂δj_matrix(
    δg, Eg, Yint, nothing )

    M_inv = (ωs/2) .* Diagonal(inv.(H))

    Aω_matrix = M_inv * minus_∂Pei_∂δj_matrix
    
    return Aω_matrix

end
            


function get_electro_mechanical_mode_matrices(
    δg,
    Eg,
    H,    
    Yint,
    ωs )

    rows_size, cols_size = size(Yint)

    Aω_matrix = get_Aω_matrix(
        δg, Eg, H, Yint, ωs, nothing)

    # A = vcat( hcat( zeros( rows_size,rows_size),
    #                 Matrix(1.0I, rows_size, rows_size) ),
    #           hcat( Aω_matrix,
    #                 zeros(rows_size, rows_size )) )


    A_matrix = vcat( hcat( zeros( rows_size, rows_size),
                    1.0 * LinearAlgebra.I( rows_size) ),
              hcat( Aω_matrix,
                    zeros(rows_size, rows_size )) )
    
    return (; Aω_matrix, A_matrix)

end
            

#-----------------------------------------------------
# internal node model ode
#-----------------------------------------------------


function get_Pe_inm(
    gens_δ,
    gens_E,
    Yint)
    
    E_gens = gens_E .* exp.(im * gens_δ)
    
    I_gens = Yint * E_gens
    
    return real.( E_gens .* conj.(I_gens))
    

end


function get_Pe_inm(
    gens_δ,
    gens_E,
    Gint,
    Bint)

    # Gint = real.(Yint)
    
    # Bint = imag.(Yint)

    no_of_gens = length(gens_δ)

    return  [
        sum([gens_E[idx_i] * gens_E[idx_j] * (
            Gint[idx_i,idx_j] *
                cos(gens_δ[idx_i] - gens_δ[idx_j]) +
            Bint[idx_i,idx_j] *
                sin(gens_δ[idx_i] -  gens_δ[idx_j]) )
             for idx_j in 1:no_of_gens])
        for idx_i in 1:no_of_gens ]

end


function ode_a_gen_inm_swing!(
    dx,
    x,
    a_gen_inm_para, t;
    kwd_para =
         a_gen_inm_kwd_para )

    ( a_Δω, a_Tm,
      a_Pe_i ) =
        a_gen_inm_para

    (a_H,
     a_D,
      ωs) =
          kwd_para    

    δ, ω = x

    dx[1] = ω + a_Δω - ωs 

    dx[2 ] = ( ωs/(2 * a_H) ) *
        ( (a_Tm  - a_Pe_i)  -
         a_D * ( ω - ωs + a_Δω )

          )


    return nothing
    
end


function ode_gens_inm_swing!(
    dx, x, p, t;
    kwd_para = kwd_para )

    (; Δω,
     gens_Tm,
      Pe ) = p

    ( gens_H,
      gens_D,
      ωs,
      state_vars_idx) =
          kwd_para

    dx_views = [view(dx, idx)
                for idx in state_vars_idx ]
    
    x_views =  [view(x, idx)
                for idx in state_vars_idx ]
    
    for (a_dx_view, a_x_view, a_Δω,a_Tm,a_Pe_i,a_H,a_D ) in
        zip(dx_views, x_views, Δω, gens_Tm, Pe,
            gens_H, gens_D )

        a_gen_inm_para =
            (a_Δω, a_Tm, a_Pe_i)
        
        a_gen_inm_kwd_para =
            (a_H, a_D, ωs)
        
        ode_a_gen_inm_swing!(
            # dx[state_var_idx],
            # x[state_var_idx],
            a_dx_view,
            a_x_view,
            a_gen_inm_para,
            t;
            kwd_para =
                a_gen_inm_kwd_para )

    end

    return nothing
    
end



function ode_gens_inm_by_Yint_swing!(
    dx, x, p, t;
    kwd_para = kwd_para )

    (; gens_δ,
      gens_E,
      gens_Tm,
      Yint ) = p

    # Yint = fault_state_Yint[1]
    
    ( ωs,
      gens_para,
      state_vars_idx,
      gens_nodes_idx) =
          kwd_para
    
    ( gens_H,
      gens_X_d_dash ) =
          gens_para

    # see eq 7.209, Sauer

    # note gens_E = abs(E_gens)
    
    E_gens = gens_E .* exp.(im * gens_δ)
    
    I_gens = Yint * E_gens
    
    Pe = real.( E_gens .* conj.(I_gens))
    
    dx_views = [view(dx, idx)
                for idx in state_vars_idx ]
    
    x_views =  [view(x, idx)
                for idx in state_vars_idx ]

    for (a_dx_view, a_x_view, Tm, Pe_i,  H, ) in
        zip(dx_views, x_views, gens_Tm, Pe, gens_H )

        a_gen_inm_para =
            (Tm, Pe_i)
        
        a_gen_inm_kwd_para =
            (H, ωs)
        
        ode_a_gen_inm_swing!(
            # dx[state_var_idx],
            # x[state_var_idx],
            a_dx_view,
            a_x_view,
            a_gen_inm_para,
            t;
            kwd_para =
                a_gen_inm_kwd_para )

    end

    return nothing
    
end


#-----------------------------------------------------
# internal node model powerflow
#-----------------------------------------------------

# N/A


#-----------------------------------------------------
#-----------------------------------------------------


function get_spcm_parameters_and_idx(
    net_nodes_type_idxs )
    
    #--------------------------------------
    
    flat_vh_flat_θh_Idx =
        flat_vh_flat_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_vh,
              dim_θh ];
            dims_given = true )

    (nodes_flat_vh_Idx,
     nodes_flat_θh_Idx) =
         flat_vh_flat_θh_Idx

    gens_vh_idxs =
        nodes_flat_vh_Idx[ gens_nodes_idx ]
    
    gens_θh_idxs =
        nodes_flat_θh_Idx[gens_nodes_idx]

    non_gens_nodes_vh_idxs =
        nodes_flat_vh_Idx[ non_gens_nodes_idx ]
    
    non_gens_nodes_θh_idxs =
        nodes_flat_vh_Idx[ non_gens_nodes_idx ]

    non_pre_ordered_pf_vars_Idxs =
        (;gens_vh_idxs,
         gens_θh_idxs,
         non_gens_nodes_vh_idxs,
         non_gens_nodes_θh_idxs,
         flat_vh_flat_θh_Idx)
    
    #--------------------------------------
        
    gens_vh_idx_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ no_of_gens,
              no_of_gens ] ;
            dims_given = true )

    
    flat_Pg_flat_Qg_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_Pg,
              dim_Qg ] ;
            dims_given = true )

    dim_Pg_Qg_Png_Qng_Pll_Qll  =
            [ dim_P_gens,
              dim_Q_gens,
              dim_P_non_gens,
              dim_Q_non_gens,
              dim_P_g_loc_load,
              dim_Q_g_loc_load]
    
    pf_Pg_Qg_Png_Qng_Pll_Qll_Idxs =
        get_vars_or_paras_Idxs_in_flattend(
            dim_Pg_Qg_Png_Qng_Pll_Qll;
            dims_given = true )

    (
     P_gens_Idxs,
     Q_gens_Idxs,
     P_non_gens_Idxs,
     Q_non_gens_Idxs,
     P_g_loc_load_Idxs,
        Q_g_loc_load_Idxs ) =
            pf_Pg_Qg_Png_Qng_Pll_Qll_Idxs
    

end



# ------------------------------------------------------
# ------------------------------------------------------


function get_a_case_generic_para_and_Idxs_by_json(
    ;case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true)

    case_name = "case9"

    #--------------------------------------

    # case_name = "case9"
        
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")
    
    if (json_net_data_by_components_file == "") || (
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
    
    #-----------------------------------------------

    net_generic_parameters_and_idx =
        get_net_generic_parameters_and_idx(
            net_data_by_components_file;
                    
            basekV = basekV,    
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu,
            
            in_components_type_sym =
                false )
    
    (;data_by_components,            
     net_nodes_type_idxs,
     n2s_idxs,
     net_generic_idx,
     net_generic_parameters ) =
         net_generic_parameters_and_idx
end


function get_selected_generic_Idxs_wt_nodes_Idx_by_json(
    ; case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    
    selected_Idxs_syms =
        (:loc_load_exist,
         :dyn_pf_δ_eq_dash_0_P_Q_idx,
         :Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         :dyn_pf_fun_kwd_n2s_idxs,
         :dyn_pf_fun_kwd_net_idxs,
         :non_pre_ordered_pf_vars_Idxs,

         :dyn_pf_δ_eq_dash_Idx,
         :dyn_pf_Png_Qng_Pll_Qll_Idx,

         :dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,

         :flux_decay_model_states_Idx,
         :flux_decay_model_states_comp_idxs_in_Idx,

         :flux_decay_model_vars_Idx_in_state))

    #--------------------------------------    
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
    
    #----------------------------------------
    
   selected_Idxs =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                    net_data_by_components_file;
                    in_components_type_sym =
                        false),
                 :net_generic_idx),
             selected_Idxs_syms )

    #----------------------------------------

    if :dyn_pf_fun_kwd_net_idxs ∈ selected_Idxs_syms
    
    (gens_nodes_idx,
     non_gens_nodes_idx,
        all_nodes_idx) =
            NamedTupleTools.select(
                dyn_pf_fun_kwd_net_idxs,
                (:gens_nodes_idx,
                 :non_gens_nodes_idx,
                 :all_nodes_idx))         
    else

       dyn_pf_fun_kwd_net_idxs =
             NamedTupleTools.select(
                 getproperty(
                     get_net_generic_parameters_and_idx(
                        net_data_by_components_file;
                        in_components_type_sym =
                            false),
                     :net_generic_idx),
                 (:dyn_pf_fun_kwd_net_idxs,) )
    
    (gens_nodes_idx,
     non_gens_nodes_idx,
        all_nodes_idx) =
            NamedTupleTools.select(
                dyn_pf_fun_kwd_net_idxs,
                (:gens_nodes_idx,
                 :non_gens_nodes_idx,
                 :all_nodes_idx))                 
        
    end

    return (; gens_nodes_idx,
            non_gens_nodes_idx,
            all_nodes_idx,
            selected_Idxs )

end


function get_selected_generic_parameters_by_json(
    ; case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    selected_para_syms =
        (:Ynet_wt_nodes_idx_wt_adjacent_nodes,
         :ode_gens_para,

         :pf_sta_ΔPQ_mismatch_parameters,

         :Ybr_cal_and_edges_orientation,
         :sta_pf_PQ_para,
         :ode_gens_generic_para,
         :baseMVA,

         :generic_govs_para,
         :generic_avrs_para,
         :ode_flux_decay_avrs_para,
         :dyn_pf_mismatch_vars_kwd_para) )

    #--------------------------------------    
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")

    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
    
    
    #----------------------------------------
    
         return NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                     net_data_by_components_file;
            
                    basekV = basekV,    
                    use_pu_in_PQ = use_pu_in_PQ,
                    line_data_in_pu = line_data_in_pu,
                     
                    in_components_type_sym =
                        false),
                 :net_generic_parameters),
             selected_para_syms )

end


function get_selected_para_and_Idxs_wt_nodes_Idx_by_json(
    ; case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    selected_para_syms =
        (:Ynet_wt_nodes_idx_wt_adjacent_nodes,
         :ode_gens_para,

         :pf_sta_ΔPQ_mismatch_parameters,

         :Ybr_cal_and_edges_orientation,
         :sta_pf_PQ_para,
         :ode_gens_generic_para,
         :baseMVA,

         :generic_govs_para,
         :generic_avrs_para,
         :ode_flux_decay_avrs_para ) ,
    
    selected_Idxs_syms =
        (:loc_load_exist,
         :dyn_pf_δ_eq_dash_0_P_Q_idx,
         :Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         :dyn_pf_fun_kwd_n2s_idxs,
         :dyn_pf_fun_kwd_net_idxs,
         :non_pre_ordered_pf_vars_Idxs,

         :dyn_pf_δ_eq_dash_Idx,
         :dyn_pf_Png_Qng_Pll_Qll_Idx,

         :dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,

         :flux_decay_model_states_Idx,
         :flux_decay_model_states_comp_idxs_in_Idx,

         :flux_decay_model_vars_Idx_in_state))

    #--------------------------------------    
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")

    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
    
        
    #----------------------------------------

    selected_para =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                     net_data_by_components_file;
            
                    basekV = basekV,    
                    use_pu_in_PQ = use_pu_in_PQ,
                    line_data_in_pu = line_data_in_pu,
                     
                    in_components_type_sym =
                        false),
                 :net_generic_parameters),
             selected_para_syms )
    
    dyn_pf_mismatch_vars_kwd_para =
        (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para )
    
    #----------------------------------------    

   selected_Idxs =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                    net_data_by_components_file;
                    in_components_type_sym =
                        false),
                 :net_generic_idx),
             selected_Idxs_syms )

    #----------------------------------------

    if :dyn_pf_fun_kwd_net_idxs ∈ selected_Idxs_syms
    
    (gens_nodes_idx,
     non_gens_nodes_idx,
        all_nodes_idx) =
            NamedTupleTools.select(
                dyn_pf_fun_kwd_net_idxs,
                (:gens_nodes_idx,
                 :non_gens_nodes_idx,
                 :all_nodes_idx))         
    else

       dyn_pf_fun_kwd_net_idxs =
             NamedTupleTools.select(
                 getproperty(
                     get_net_generic_parameters_and_idx(
                        net_data_by_components_file;
                        in_components_type_sym =
                            false),
                     :net_generic_idx),
                 (:dyn_pf_fun_kwd_net_idxs,) )
    
    (gens_nodes_idx,
     non_gens_nodes_idx,
        all_nodes_idx) =
            NamedTupleTools.select(
                dyn_pf_fun_kwd_net_idxs,
                (:gens_nodes_idx,
                 :non_gens_nodes_idx,
                 :all_nodes_idx))                 
        
    end

    return (; gens_nodes_idx,
            non_gens_nodes_idx,
            all_nodes_idx,
            selected_para,
            selected_Idxs )

end


#---------------------------------------------------
#---------------------------------------------------

function get_flux_decay_dynamic_pf_sol_para_by_json(
    ; case_name     = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",
    basekV          = 1.0,    
    use_pu_in_PQ    = true,
    line_data_in_pu = true,
        
    pf_alg        = NewtonRaphson(),
    ode_alg       = Rodas4(),
    
    abstol        = 1e-12,
    reltol        = 1e-12,
    
    dt            = 0.01,
    timespan      = 10.0,
    
    tspan         = (0.0, timespan),
    sim_timespan  = (0.0, timespan),
    plot_timespan = (0.0, timespan))
    
    #--------------------------------------    
    
    pf_solver =
        (; pf_alg,
         abstol,
         reltol)        

    use_nlsolve = false
    
    #--------------------------------------    
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
    
    
    
    #----------------------------------------
    # parameters
    #----------------------------------------

    (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
     ode_gens_para,

     pf_sta_ΔPQ_mismatch_parameters,

     Ybr_cal_and_edges_orientation,
     sta_pf_PQ_para,
     ode_gens_generic_para,
     baseMVA,

     generic_govs_para,
     generic_avrs_para,
     ode_flux_decay_avrs_para) =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                     net_data_by_components_file;
                     
                     basekV = basekV,    
                     use_pu_in_PQ = use_pu_in_PQ,
                     line_data_in_pu = line_data_in_pu,
                     
                    in_components_type_sym =
                        false),
                 :net_generic_parameters),
             (:Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :ode_gens_para,

              :pf_sta_ΔPQ_mismatch_parameters,
              
              :Ybr_cal_and_edges_orientation,
              :sta_pf_PQ_para,
              :ode_gens_generic_para,
              :baseMVA,

              :generic_govs_para,
              :generic_avrs_para,
              :ode_flux_decay_avrs_para) )
    
    #----------------------------------------    
    # Indices
    #----------------------------------------    

    (loc_load_exist,
     dyn_pf_δ_eq_dash_0_P_Q_idx,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs,
     # non_pre_ordered_dyn_pf_vars_Idxs,

     dyn_pf_δ_eq_dash_Idx,
     dyn_pf_Png_Qng_Pll_Qll_Idx,

     dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,

     dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,

     flux_decay_model_states_Idx,
     flux_decay_model_states_comp_idxs_in_Idx,

     flux_decay_model_vars_Idx_in_state) =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                     net_data_by_components_file;
                     
                     basekV = basekV,    
                     use_pu_in_PQ = use_pu_in_PQ,
                     line_data_in_pu = line_data_in_pu,
                     
                    in_components_type_sym =
                        false),
                 :net_generic_idx),
             (:loc_load_exist,
              :dyn_pf_δ_eq_dash_0_P_Q_idx,
              :Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              :non_pre_ordered_pf_vars_Idxs,
              # :non_pre_ordered_dyn_pf_vars_Idxs,

              :dyn_pf_δ_eq_dash_Idx,
              :dyn_pf_Png_Qng_Pll_Qll_Idx,

              :dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,

              :dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,

              :flux_decay_model_states_Idx,
              :flux_decay_model_states_comp_idxs_in_Idx,

              :flux_decay_model_vars_Idx_in_state) )

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx) =  NamedTupleTools.select(
        dyn_pf_fun_kwd_net_idxs,
         (:gens_nodes_idx,
          :non_gens_nodes_idx,
          :all_nodes_idx)) 

    #----------------------------------------

    state_vars_syms_flux_decay =
        [:δ, :ω, :eq_dash, :E_fd ]
    
    ode_flux_decay_state_sym =
        generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            state_vars_syms_flux_decay;
            label_prefix = "bus" )

    pf_flux_decay_state_sym =
        [generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:vh];
            label_prefix = "bus");
         generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:θh];
            label_prefix = "bus")]

    flux_decay_model_syms =
        [ode_flux_decay_state_sym;
         pf_flux_decay_state_sym]
    
    #----------------------------------------
    # states in system_dynamics_flux_decay_model
    #----------------------------------------
    
    flux_decay_model_mass_matrix =
        DAE_MassMatrix(
            length(ode_flux_decay_state_sym),
            length(pf_flux_decay_state_sym) )
    
    #----------------------------------------
    
    dyn_pf_mismatch_vars_kwd_para =
        (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para )

    #----------------------------------------
    
    dyn_pf_mismatch_idx_kwd_para =
        (;dyn_pf_δ_eq_dash_0_P_Q_idx,
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         non_pre_ordered_pf_vars_Idxs,
         dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
         dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx)
    
    #----------------------------------------

    dyn_pf_mismatch_kwd_para =
        (; loc_load_exist,
         dyn_pf_mismatch_vars_kwd_para,
         dyn_pf_mismatch_idx_kwd_para )

    #----------------------------------------

    sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx =
        (;dyn_pf_δ_eq_dash_Idx,
         dyn_pf_Png_Qng_Pll_Qll_Idx)

    flux_decay_model_states_Idx_wt_idxs_in_Idx =
        (;flux_decay_model_states_Idx,
         flux_decay_model_states_comp_idxs_in_Idx,
         flux_decay_model_vars_Idx_in_state )

    #----------------------------------------    
    #----------------------------------------
    
    by_part_dyn_pf_mismatch_kwd_para =
        (;pf_solver,
         use_nlsolve,
         dyn_pf_mismatch_kwd_para,
         sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx )

    #----------------------------------------

    (pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
         NamedTupleTools.select(
             pf_sta_ΔPQ_mismatch_parameters,
             (:pf_kw_para,
              :red_types_Idxs_etc,
              :pf_PQ_param) )

    #----------------------------------------
    
    (red_vh_Idxs,
     red_θh_Idxs) =
         NamedTupleTools.select(
             red_types_Idxs_etc,
             (:red_vh_Idxs,
              :red_θh_Idxs) )
    
    #----------------------------------------
    
    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs))]

    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
          
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------
    
    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (; Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_results_pf_sta_red_sol =
        get_generic_results_pf_sta_red_sol_u(
            pf_sol;
            generic_red_sol_kwd_para =
                generic_red_sol_kwd_para,
            baseMVA =
                baseMVA,
            basekV =
                1.0 )

    #----------------------------------------    
    #----------------------------------------    
    
    (pf_P_gens, pf_Q_gens,
     pf_P_g_gens, pf_Q_g_gens,
     Igen,
     vh, θh, θh_deg,
     gens_vh, gens_θh,
     gens_nodes_idx) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens, :pf_Q_gens,
             :pf_P_g_gens, :pf_Q_g_gens,
             :Igen,
             :vh, :θh, :θh_deg,
             :gens_vh, :gens_θh,
             :gens_nodes_idx ) )

    #----------------------------------------
    
    # P_g = pf_P_gens #./ baseMVA
    
    # Q_g = pf_Q_gens # ./ baseMVA
        
    #----------------------------------------
    # Dynamic simulation
    #----------------------------------------

    (X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            ode_gens_generic_para,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash ) )
    
    #----------------------------------------
    
    (Ka,
     Ta) = ode_flux_decay_avrs_para

    #----------------------------------------
    
    nt_init_flux_decay =
        get_init_flux_decay_model(
            gens_vh,
            gens_θh,            
            pf_P_g_gens,
            pf_Q_g_gens,            
            X_d,
            X_q,            
            X_d_dash,
            X_q_dash,
            Ka)
    
    #----------------------------------------
    
    init_flux_decay_vars =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections =
            (:δ, :eq_dash, :E_fd,
             :V_ref, :Tm ) )
    
    init_flux_decay_V_ref_Tm =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm ) )

    
    init_flux_decay_δ_eq_E_fd =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:δ, :eq_dash, :E_fd ) )

   init_flux_decay_i_d_i_q =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:i_d, :i_q ) )

    #----------------------------------------

    (V_ref,
     Tm,
     gens_i_d,
     gens_i_q) =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm, :i_d, :i_q ))
    
    #----------------------------------------
    
    ode_flux_decay_gens_para =
        NamedTupleTools.select(
            ode_gens_generic_para,
                (:H,  :X_d, :X_q, :X_d_dash,
                 :T_d_dash ) )
        
    state_vars_idx =
        get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
            state_vars_syms_flux_decay,
            gens_nodes_idx )

    #----------------------------------------

    ode_flux_decay_para =
        ( gens_vh,
          gens_θh,
          V_ref,
          Tm)

    ode_flux_decay_kwd_para = (;
        ode_flux_decay_gens_para,
        ode_flux_decay_avrs_para,
        state_vars_idx,
        ωs)
     
    #----------------------------------------    
    # ode states init
    #----------------------------------------    
    
    state_init_flux_decay =
        get_state_init_flux_decay(
            init_flux_decay_δ_eq_E_fd,
            ωs)
    
    #----------------------------------------    
    
    (dyn_δ,
     dyn_eq_dash_0) =
        NamedTupleTools.select(
            init_flux_decay_δ_eq_E_fd,
            (:δ,
             :eq_dash ) )

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

    if loc_load_exist == true

        flux_decay_dyn_pf_flat_para =
            Float64[dyn_δ;
             dyn_eq_dash_0;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load ]
        
    else

        flux_decay_dyn_pf_flat_para =
            Float64[dyn_δ;
             dyn_eq_dash_0;
             P_non_gens;
             Q_non_gens ]
        
    end

    # -----------------------------------

    # flux_decay_dyn_pf_flat_para =
    #     convert(Vector{Float64},
    #             flux_decay_dyn_pf_flat_para )


    flux_decay_dyn_pf_δ_eq_dash_flat_para =
        Float64[dyn_δ;
         dyn_eq_dash_0]
    
    flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para =
        Float64[P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load ]

    flux_decay_dyn_pf_Pg_Qg_Png_Qng_Pgll_Qgll_flat_para =
        Float64[pf_P_g_gens;
         pf_Q_g_gens;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load ]

    #----------------------------------------    
    
    (flat_vh_Idx,
     flat_θh_Idx) =
         NamedTupleTools.select(
             getproperty(
                 getproperty(
                     dyn_pf_mismatch_idx_kwd_para,
                     :non_pre_ordered_pf_vars_Idxs ),
                 :flat_vh_flat_θh_flat_id_iq_Idx),
             (:flat_vh_Idx,
              :flat_θh_Idx))
    
    #----------------------------------------    
    
    init_flat_vh_flat_θh = [vh; θh]

    # init_flat_vh_flat_θh_id_iq =
    #     Float64[vh;
    #      θh;
    #      init_flux_decay_i_d_i_q.i_d;
    #      init_flux_decay_i_d_i_q.i_q ]

    #----------------------------------------
    
    flux_decay_model_Pq_Qg_id_iq =
        get_flux_decay_model_Pq_Qg_id_iq_by_parts(
            init_flat_vh_flat_θh,
            flux_decay_dyn_pf_δ_eq_dash_flat_para,
            flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para;
            by_part_dyn_pf_mismatch_kwd_para =
                by_part_dyn_pf_mismatch_kwd_para )
    
    #----------------------------------------    
    
    pf_sol_post_state_init =
        get_flux_decay_model_pf_wt_ext_idq_by_parts_func!(
            init_flat_vh_flat_θh,
            flux_decay_dyn_pf_δ_eq_dash_flat_para,
            flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para;
            by_part_dyn_pf_mismatch_kwd_para =
                by_part_dyn_pf_mismatch_kwd_para )

    #----------------------------------------    

    pf_sol_post_state_init_vh =
        pf_sol_post_state_init[
            flat_vh_Idx]
    
    pf_sol_post_state_init_θh =
        pf_sol_post_state_init[
            flat_θh_Idx]
    
    pf_state_init_flux_decay = [
        pf_sol_post_state_init_vh;
        pf_sol_post_state_init_θh]
    
    #----------------------------------------    

    flux_decay_model_states_init =
        [state_init_flux_decay;
         pf_state_init_flux_decay]

    #----------------------------------------    
    flux_decay_model_dynamics_para =
        ComponentVector{Float64}(
            V_ref = V_ref,
            Tm = Tm,
            Png = P_non_gens,
            Qng = Q_non_gens, 
            Pll = P_g_loc_load,
            Qll = Q_g_loc_load)
    
    #----------------------------------------
    # PreallocationTools.jl
    # https://github.com/SciML/PreallocationTools.jl
    #----------------------------------------
    
    stateDiffCache = DiffCache(
        flux_decay_model_states_init)
    
    #----------------------------------------
    # flux_decay_model_dynamics_kwd_para
    #----------------------------------------
    
    flux_decay_model_dynamics_kwd_para =
        (;stateDiffCache,
         ode_flux_decay_kwd_para,
         by_part_dyn_pf_mismatch_kwd_para,
         sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx,         
         flux_decay_model_states_Idx_wt_idxs_in_Idx )
    
    #----------------------------------------    

    return (;

            sta_pf_PQ_para,
            pf_PQ_param,
            
            ode_gens_generic_para, #
            
            ode_flux_decay_avrs_para,
            
            generic_red_sol_kwd_para,
            
            generic_results_pf_sta_red_sol,
                        
            state_vars_syms_flux_decay,
            
            init_flux_decay_i_d_i_q,
            
            state_init_flux_decay,
            
            nt_init_flux_decay,

            init_flux_decay_vars,
            init_flux_decay_V_ref_Tm,

            ode_flux_decay_gens_para,
            state_vars_idx,
            
            ode_flux_decay_para,
            ode_flux_decay_kwd_para,

            init_flux_decay_δ_eq_E_fd,
            
            # init_flat_vh_flat_θh_id_iq,
            
            init_flat_vh_flat_θh,

            flux_decay_dyn_pf_flat_para,
            flux_decay_dyn_pf_δ_eq_dash_flat_para,
            flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para,
            flux_decay_dyn_pf_Pg_Qg_Png_Qng_Pgll_Qgll_flat_para,
            by_part_dyn_pf_mismatch_kwd_para,
            non_pre_ordered_pf_vars_Idxs,
            
            flux_decay_model_Pq_Qg_id_iq,
            
            flux_decay_model_states_init,
            flux_decay_model_dynamics_para,
            flux_decay_model_dynamics_kwd_para)
end


function get_a_case_Y_aug_matrices_by_json(
    case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    pf_alg  = NewtonRaphson(),
    abstol  = 1e-12,
    reltol  = 1e-12)

    #--------------------------------------    
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")
    
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
    
    
    #----------------------------------------

    (;pf_sta_ΔPQ_mismatch_parameters,
     Ybr_cal_and_edges_orientation,
     sta_pf_PQ_para,
     ode_gens_generic_para,
     baseMVA) =
         NamedTupleTools.select(
            getproperty(
                get_net_generic_parameters_and_idx(
                    net_data_by_components_file;
                    
                     basekV = basekV,    
                     use_pu_in_PQ = use_pu_in_PQ,
                     line_data_in_pu = line_data_in_pu,
                    
                    in_components_type_sym =
                        false),
                :net_generic_parameters),
         (:pf_sta_ΔPQ_mismatch_parameters,
          :Ybr_cal_and_edges_orientation,
          :sta_pf_PQ_para,
          :ode_gens_generic_para,
          :baseMVA) )
    
    #----------------------------------------

    (;pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
        NamedTupleTools.select(
            pf_sta_ΔPQ_mismatch_parameters,
            (:pf_kw_para,
             :red_types_Idxs_etc,
             :pf_PQ_param))

    #----------------------------------------

    (;red_vh_Idxs,
     red_θh_Idxs) =
         NamedTupleTools.select(
             red_types_Idxs_etc,
             (:red_vh_Idxs,
              :red_θh_Idxs))
    
    #----------------------------------------

    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs))]

    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
    
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------

    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (; Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_results_pf_sta_red_sol =
        get_generic_results_pf_sta_red_sol_u(
            pf_sol;
            generic_red_sol_kwd_para =
                generic_red_sol_kwd_para,
            baseMVA =
                baseMVA,
            basekV =
                1.0 )
    
    # (; gens_current_injection,
    # gens_loc_load_current,
    # gens_nodes_network_current,
    # non_gens_nodes_network_current,
    # nodes_network_current,
    # S_gens, pf_P_gens, pf_Q_gens, gens_S,
    # Igen, Inet, Iinj, pu_Igen,
    # Sbus_n, GenSinj,
    # If, It, Ibranches,
    # vh, θh, θh_deg, Vbus, E,
    # gens_vh, gens_θh,
    # gens_nodes_idx ) =
    #     generic_results_pf_sta_red_sol

    #----------------------------------------    

    (pf_P_gens, pf_Q_gens,
     vh, θh,
     gens_vh, gens_θh,
     gens_nodes_idx ) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens, :pf_Q_gens,
             :vh, :θh,
             :gens_vh, :gens_θh,
             :gens_nodes_idx ) )
    
    #----------------------------------------    

    # (;     
    #  P_non_gens,
    #  Q_non_gens,
    #  P_g_loc_load,
    #  Q_g_loc_load,
    #  loc_load_exist ) =
    #      NamedTupleTools.select(
    #          sta_pf_PQ_para,
    #          (:P_non_gens,
    #           :Q_non_gens,
    #           :P_g_loc_load,
    #           :Q_g_loc_load,
    #           :loc_load_exist) )

    # if loc_load_exist == true
    #     post_pf_PQ_param =
    #         Float64[pf_P_gens;
    #                 pf_Q_gens;
    #                 P_non_gens;
    #                 Q_non_gens;
    #                 P_g_loc_load;
    #                 Q_g_loc_load]
    # else
    #     post_pf_PQ_param =
    #         Float64[pf_P_gens;
    #                 pf_Q_gens;
    #                 P_non_gens;
    #                 Q_non_gens]
    # end
    
    #----------------------------------------    
    #----------------------------------------    

    ( X_d_dash, ) =
         NamedTupleTools.select(
            ode_gens_generic_para,
                ( :X_d_dash, ) )
            
    #----------------------------------------
    #########################################
    #----------------------------------------

    y_aug_kw_para =
        (;X_d_dash,
         pf_kw_para )
    
    Y_aug_matrices =
        get_Y_aug_matrices(
            pf_PQ_param,
            vh;
            y_aug_kw_para =
                y_aug_kw_para )

    (;Yred,
     Y_internal_nodes) =
        NamedTupleTools.select(
            Y_aug_matrices,
            (:Yred,
             :Y_internal_nodes))

    return (;Y_aug_matrices,
            generic_results_pf_sta_red_sol)
end


# ------------------------------------------------------

function get_static_pf_sol_by_json(
    ; case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",    
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    pf_alg  = NewtonRaphson())

    #----------------------------------------
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
    
    
    #----------------------------------------
    
    (pf_sta_ΔPQ_mismatch_parameters,
     sta_pf_PQ_para,
     ode_gens_generic_para,
     baseMVA,
     Ybr_cal_and_edges_orientation)  =
        NamedTupleTools.select(
            getproperty(
                get_net_generic_parameters_and_idx(
                    net_data_by_components_file;

                    basekV = basekV,    
                    use_pu_in_PQ = use_pu_in_PQ,
                    line_data_in_pu = line_data_in_pu,
                    
                    in_components_type_sym =
                        false),
                :net_generic_parameters),
            (:pf_sta_ΔPQ_mismatch_parameters,
             :sta_pf_PQ_para,
             :ode_gens_generic_para,             
             :baseMVA,
             :Ybr_cal_and_edges_orientation))

    
    #----------------------------------------

    (;pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
         NamedTupleTools.select(
             pf_sta_ΔPQ_mismatch_parameters,
             (:pf_kw_para,
              :red_types_Idxs_etc,
              :pf_PQ_param))

    loc_load_exist = getproperty(
        pf_kw_para,
        :loc_load_exist)

    #----------------------------------------

    (;red_vh_Idxs,
     red_θh_Idxs) =
        NamedTupleTools.select(
            red_types_Idxs_etc,
            (:red_vh_Idxs,
             :red_θh_Idxs))
    
    #----------------------------------------

    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs))]

    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
          
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------

    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (; Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_results_pf_sta_red_sol =
        get_generic_results_pf_sta_red_sol_u(
            pf_sol;
            generic_red_sol_kwd_para =
                generic_red_sol_kwd_para,
            baseMVA =
                baseMVA,
            basekV =
                1.0 )
    
    # (; gens_current_injection,
    # gens_loc_load_current,
    # gens_nodes_network_current,
    # non_gens_nodes_network_current,
    # nodes_network_current,
    # S_gens, pf_P_gens, pf_Q_gens, gens_S,
    # Igen, Inet, Iinj, pu_Igen,
    # Sbus_n, GenSinj,
    # If, It, Ibranches,
    # vh, θh, θh_deg, Vbus, E,
    # gens_vh, gens_θh,
    # gens_nodes_idx ) =
    #     generic_results_pf_sta_red_sol

    #----------------------------------------    
    
    (pf_P_gens, pf_Q_gens,
     vh, θh, θh_deg,
     gens_vh, gens_θh,
     gens_nodes_idx) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens, :pf_Q_gens,
             :vh, :θh, :θh_deg,
             :gens_vh, :gens_θh,
             :gens_nodes_idx ) )
    
    #----------------------------------------    
    #----------------------------------------    

    return (vh,
            θh,
            θh_deg,
            pf_P_gens,
            pf_Q_gens,
            gens_nodes_idx)
    
end


function get_flux_decay_dynamic_pf_sol_by_json(
    ;case_name      = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",    
    basekV          = 1.0,    
    use_pu_in_PQ    = true,
    line_data_in_pu = true,
        
    pf_alg          = NewtonRaphson(),
    ode_alg         = Rodas4(),
    
    abstol          = 1e-12,
    reltol          = 1e-12,
    
    dt              = 0.01,
    timespan        = 10.0,
    
    tspan           = (0.0, timespan),
    sim_timespan    = (0.0, timespan),
    plot_timespan   = (0.0, timespan),
    
    wt_id_iq        = false )

    #--------------------------------------
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end


    #--------------------------------------

    flux_decay_dynamic_pf_sol_para =
        get_flux_decay_dynamic_pf_sol_para_by_json(
            ;case_name      = case_name,
            components_libs_dir=components_libs_dir,
            json_net_data_by_components_file =
                json_net_data_by_components_file,
            data_dir = data_dir,            
            basekV          = 1.0,    
            use_pu_in_PQ    = true,
            line_data_in_pu = true,

            pf_alg          = NewtonRaphson(),
            ode_alg         = Rodas4(),

            abstol          = 1e-12,
            reltol          = 1e-12,

            dt              = 0.01,
            timespan        = 10.0,

            tspan           = (0.0, timespan),
            sim_timespan    = (0.0, timespan),
            plot_timespan   = (0.0, timespan))
    

    #--------------------------------------    
    
    (;generic_results_pf_sta_red_sol,
     init_flat_vh_flat_θh,
     # init_flat_vh_flat_θh_id_iq,
     flux_decay_dyn_pf_flat_para,
     flux_decay_dyn_pf_δ_eq_dash_flat_para,
     flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para,
     flux_decay_dyn_pf_Pg_Qg_Png_Qng_Pgll_Qgll_flat_para,
     by_part_dyn_pf_mismatch_kwd_para,
     non_pre_ordered_pf_vars_Idxs,
     
     init_flux_decay_i_d_i_q,
     ode_gens_generic_para,
     init_flux_decay_δ_eq_E_fd) =
         NamedTupleTools.select(
             flux_decay_dynamic_pf_sol_para,
             (:generic_results_pf_sta_red_sol,
              :init_flat_vh_flat_θh,
              # :init_flat_vh_flat_θh_id_iq,
              :flux_decay_dyn_pf_flat_para,
              :flux_decay_dyn_pf_δ_eq_dash_flat_para,
              :flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para,
              :flux_decay_dyn_pf_Pg_Qg_Png_Qng_Pgll_Qgll_flat_para,
              :by_part_dyn_pf_mismatch_kwd_para,
              :non_pre_ordered_pf_vars_Idxs,
              
              :init_flux_decay_i_d_i_q,
              :ode_gens_generic_para,
              :init_flux_decay_δ_eq_E_fd
              ))

    #----------------------------------------    
    
    (;flat_vh_Idx,
     flat_θh_Idx,
     flat_id_Idx,
     flat_iq_Idx) =
        NamedTupleTools.select(
            getproperty(
                non_pre_ordered_pf_vars_Idxs,
                :flat_vh_flat_θh_flat_id_iq_Idx) ,
               (:flat_vh_Idx,
                :flat_θh_Idx,
                :flat_id_Idx,
                :flat_iq_Idx))
    
    #----------------------------------------    

    if wt_id_iq == true

        pf_sol_post_state_init =
           get_flux_decay_model_pf_wt_ext_idq_by_parts_func!(
                init_flat_vh_flat_θh,
                flux_decay_dyn_pf_δ_eq_dash_flat_para,
                flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para;
                by_part_dyn_pf_mismatch_kwd_para =
                    by_part_dyn_pf_mismatch_kwd_para )

        pf_sol_post_state_init_vh =
            pf_sol_post_state_init[
                flat_vh_Idx]

        pf_sol_post_state_init_θh =
            pf_sol_post_state_init[
                flat_θh_Idx]

        pf_sol_post_state_init_id =
            pf_sol_post_state_init[
                flat_id_Idx]

        pf_sol_post_state_init_iq =
            pf_sol_post_state_init[
                flat_iq_Idx]

        return (; pf_sol_post_state_init_vh,
                pf_sol_post_state_init_θh,
                pf_sol_post_state_init_id,
                pf_sol_post_state_init_iq )
        

    else
        
        pf_sol_post_state_init =
            get_flux_decay_model_pf_ΔPQ_func!(
                init_flat_vh_flat_θh,
         flux_decay_dyn_pf_Pg_Qg_Png_Qng_Pgll_Qgll_flat_para;
            by_part_dyn_pf_mismatch_kwd_para =
                by_part_dyn_pf_mismatch_kwd_para)

        pf_sol_post_state_init_vh =
            pf_sol_post_state_init[
                flat_vh_Idx]

        pf_sol_post_state_init_θh =
            pf_sol_post_state_init[
                flat_θh_Idx]

        return (; pf_sol_post_state_init_vh,
                pf_sol_post_state_init_θh )                

    end

end

# ------------------------------------------------------
# ------------------------------------------------------

function get_ode_flux_decay_dynamics_by_json(
    case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",    
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    pf_alg        = NewtonRaphson(),
    ode_alg       = Rodas4(),
    
    abstol        = 1e-12,
    reltol        = 1e-12,
    
    dt            = 0.01,
    timespan      = 10.0 )

    #--------------------------------------    
    
    tspan         = (0.0, timespan)
    
    sim_timespan  = (0.0, timespan)
    
    plot_timespan = (0.0, timespan)
    
    #--------------------------------------    
    
    pf_solver =
        (; pf_alg,
         abstol,
         reltol)        

    use_nlsolve = false

    #--------------------------------------

    # ImplicitMidpoint(),
    # ImplicitMidpoint(autodiff=false)
    # ode_alg         = Rodas4()
    
    #--------------------------------------
    #--------------------------------------    
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
    
    
    #----------------------------------------

    (;pf_sta_ΔPQ_mismatch_parameters,
     Ybr_cal_and_edges_orientation,
     sta_pf_PQ_para,
     ode_gens_generic_para,
     baseMVA,
     generic_govs_para,
     generic_avrs_para,
     ode_flux_decay_avrs_para) =
         NamedTupleTools.select(
            getproperty(
                get_net_generic_parameters_and_idx(
                    net_data_by_components_file;
                    
                     basekV = basekV,    
                     use_pu_in_PQ = use_pu_in_PQ,
                     line_data_in_pu = line_data_in_pu,
                    
                    in_components_type_sym =
                        false),
                :net_generic_parameters),
         (:pf_sta_ΔPQ_mismatch_parameters,
          :Ybr_cal_and_edges_orientation,
          :sta_pf_PQ_para,
          :ode_gens_generic_para,
          :baseMVA,
          :generic_govs_para,
          :generic_avrs_para,
          :ode_flux_decay_avrs_para) )
    
    #----------------------------------------

    (;pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
        NamedTupleTools.select(
            pf_sta_ΔPQ_mismatch_parameters,
            (:pf_kw_para,
             :red_types_Idxs_etc,
             :pf_PQ_param))

    #----------------------------------------
    
    (red_vh_Idxs,
     red_θh_Idxs) =
         NamedTupleTools.select(
             red_types_Idxs_etc,
             (:red_vh_Idxs,
              :red_θh_Idxs) )
    
    #----------------------------------------

    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs))]

    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
    
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------

    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (; Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_results_pf_sta_red_sol =
        get_generic_results_pf_sta_red_sol_u(
            pf_sol;
            generic_red_sol_kwd_para =
                generic_red_sol_kwd_para,
            baseMVA =
                baseMVA,
            basekV =
                1.0 )

    #----------------------------------------    

    (pf_P_gens, pf_Q_gens,
     vh, θh,
     gens_vh, gens_θh,
     gens_nodes_idx ) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens, :pf_Q_gens,
             :vh, :θh,
             :gens_vh, :gens_θh,
             :gens_nodes_idx ) )
        
    #----------------------------------------
    # Dynamic simulation
    #----------------------------------------

    (X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            ode_gens_generic_para,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash ) )
    
    #----------------------------------------
    # Flux decay
    #----------------------------------------
    
    (Ka,
     Ta) = ode_flux_decay_avrs_para

    #----------------------------------------
    
    nt_init_flux_decay =
        get_init_flux_decay_model(
            gens_vh,
            gens_θh,            
            pf_P_gens,
            pf_Q_gens,            
            X_d,
            X_q,            
            X_d_dash,
            X_q_dash,
            Ka)
    
    #----------------------------------------
    
    init_flux_decay_vars =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections =
            (:δ, :eq_dash, :E_fd,
             :V_ref, :Tm ) )
    
    init_flux_decay_V_ref_Tm =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm ) )
    
    init_flux_decay_δ_eq_E_fd =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:δ, :eq_dash, :E_fd ) )

    #----------------------------------------

    (V_ref,
     Tm) =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm
                          ))

    #----------------------------------------
    
    ode_flux_decay_gens_para =
        NamedTupleTools.select(
            ode_gens_generic_para,
                (:H,  :X_d, :X_q, :X_d_dash,
                 :T_d_dash ) )

    state_vars_syms_flux_decay =
        [:δ, :ω, :eq_dash, :E_fd ]
        
    state_vars_idx =
        get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
            state_vars_syms_flux_decay,
            gens_nodes_idx )

    #----------------------------------------

    ode_flux_decay_para =
        ( ; gens_vh,
          gens_θh,
          V_ref,
          Tm)

    ode_flux_decay_kwd_para = (
        ode_flux_decay_gens_para,
        ode_flux_decay_avrs_para,
        state_vars_idx,
        ωs)
    
    ode_flux_decay_state_sym =
        generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            state_vars_syms_flux_decay;
            label_prefix = "bus" )

    state_init_flux_decay =
        get_state_init_flux_decay(
            init_flux_decay_δ_eq_E_fd,
            ωs)
        
    #----------------------------------------
    #----------------------------------------
    
    ode_fun_flux_decay = ODEFunction(
        (dx,x,p,t) ->
            ode_flux_decay!(
                dx, x,
                ode_flux_decay_para, t;
                kwd_para =
                    ode_flux_decay_kwd_para);
        syms =
            ode_flux_decay_state_sym)
    
    ode_prob_flux_decay = ODEProblem(
        ode_fun_flux_decay,
        state_init_flux_decay,
        sim_timespan,
        ode_flux_decay_para )
    
    ode_sol_flux_decay =
        DifferentialEquations.solve(
            ode_prob_flux_decay,
            ode_alg,
            abstol = abstol,
            reltol = reltol )

    #----------------------------------------

    plot_flux_decay =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :eq_dash, :E_fd ],
            tspan = sim_timespan,
            fmt = :png)

    
    δ_plot_flux_decay =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ ],
            tspan = sim_timespan,
            fmt = :png)


    ω_plot_flux_decay =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:ω ],
            tspan = sim_timespan,
            fmt = :png)

    
    eq_dash_plot_flux_decay =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    
    E_fd_plot_flux_decay =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:E_fd ],
            tspan = sim_timespan,
            fmt = :png)

    vars_plots_flux_decay = [
        δ_plot_flux_decay,
        ω_plot_flux_decay,
        eq_dash_plot_flux_decay,
        E_fd_plot_flux_decay]

    plt_layout = (2,2)
    
    plt = plot(
        vars_plots_flux_decay...;
        layout = plt_layout,
        size = (1000, 500),
        lw = 3,
        xlabel = "t[s]")
end


function get_ode_flux_decay_wt_ext_id_iq_dynamics_by_json(
    case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",    
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    pf_alg        = NewtonRaphson(),
    ode_alg       = Rodas4(),
    
    abstol        = 1e-12,
    reltol        = 1e-12,
    
    dt            = 0.01,
    timespan      = 10.0
    )

    #--------------------------------------    

    # case_name = "case9"                        
    # basekV = 1.0    
    # use_pu_in_PQ = true
    # line_data_in_pu = true
    # pf_alg        = NewtonRaphson()
    # ode_alg       = Rodas4()
    
    # abstol        = 1e-12
    # reltol        = 1e-12
    
    # dt            = 0.01
    # timespan      = 10.0
    
    #--------------------------------------    
    
    tspan         = (0.0, timespan)
    
    sim_timespan  = (0.0, timespan)
    
    plot_timespan = (0.0, timespan)
    
    #--------------------------------------    
    
    pf_solver =
        (; pf_alg,
         abstol,
         reltol)        

    use_nlsolve = false

    #--------------------------------------

    # ImplicitMidpoint(),
    # ImplicitMidpoint(autodiff=false)
    # ode_alg         = Rodas4()
    
    #--------------------------------------
    #--------------------------------------    
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
            
    #----------------------------------------

    (;pf_sta_ΔPQ_mismatch_parameters,
     Ybr_cal_and_edges_orientation,
     sta_pf_PQ_para,
     ode_gens_generic_para,
     baseMVA,
     generic_govs_para,
     generic_avrs_para,
     ode_flux_decay_avrs_para) =
         NamedTupleTools.select(
            getproperty(
                get_net_generic_parameters_and_idx(
                    net_data_by_components_file;
                    
                     basekV = basekV,    
                     use_pu_in_PQ = use_pu_in_PQ,
                     line_data_in_pu = line_data_in_pu,
                    
                    in_components_type_sym =
                        false),
                :net_generic_parameters),
         (:pf_sta_ΔPQ_mismatch_parameters,
          :Ybr_cal_and_edges_orientation,
          :sta_pf_PQ_para,
          :ode_gens_generic_para,
          :baseMVA,
          :generic_govs_para,
          :generic_avrs_para,
          :ode_flux_decay_avrs_para) )
    
    #----------------------------------------

    (;pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
        NamedTupleTools.select(
            pf_sta_ΔPQ_mismatch_parameters,
            (:pf_kw_para,
             :red_types_Idxs_etc,
             :pf_PQ_param))

    #----------------------------------------
    
    (red_vh_Idxs,
     red_θh_Idxs) =
         NamedTupleTools.select(
             red_types_Idxs_etc,
             (:red_vh_Idxs,
              :red_θh_Idxs) )
    
    #----------------------------------------

    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs))]

    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
    
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------

    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (; Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_results_pf_sta_red_sol =
        get_generic_results_pf_sta_red_sol_u(
            pf_sol;
            generic_red_sol_kwd_para =
                generic_red_sol_kwd_para,
            baseMVA =
                baseMVA,
            basekV =
                1.0 )

    #----------------------------------------    

    (pf_P_gens, pf_Q_gens,
     vh, θh,
     gens_vh, gens_θh,
     gens_nodes_idx ) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens, :pf_Q_gens,
             :vh, :θh,
             :gens_vh, :gens_θh,
             :gens_nodes_idx ) )
        
    #----------------------------------------
    # Dynamic simulation
    #----------------------------------------

    (X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            ode_gens_generic_para,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash ) )
    
    #----------------------------------------
    # Flux decay
    #----------------------------------------
    
    (Ka,
     Ta) = ode_flux_decay_avrs_para

    #----------------------------------------
    
    nt_init_flux_decay =
        get_init_flux_decay_model(
            gens_vh,
            gens_θh,            
            pf_P_gens,
            pf_Q_gens,            
            X_d,
            X_q,            
            X_d_dash,
            X_q_dash,
            Ka)
    
    #----------------------------------------
    
    init_flux_decay_vars =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections =
            (:δ, :eq_dash, :E_fd,
             :V_ref, :Tm ) )
    
    init_flux_decay_V_ref_Tm =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm ) )

    
    init_flux_decay_δ_eq_E_fd =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:δ, :eq_dash, :E_fd ) )

   init_flux_decay_i_d_i_q =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:i_d, :i_q ) )

    #----------------------------------------

    (V_ref,
     Tm,
     gens_i_d,
     gens_i_q) =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm, :i_d, :i_q ))

    #----------------------------------------
    
    ode_flux_decay_gens_para =
        NamedTupleTools.select(
            ode_gens_generic_para,
                (:H,  :X_d, :X_q, :X_d_dash,
                 :T_d_dash ) )

    state_vars_syms_flux_decay =
        [:δ, :ω, :eq_dash, :E_fd ]
        
    state_vars_idx =
        get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
            state_vars_syms_flux_decay,
            gens_nodes_idx )

    #----------------------------------------

    # ode_flux_decay_para =
    #     ( ; gens_vh,
    #       gens_θh,
    #       V_ref,
    #       Tm)
    
    ode_flux_decay_by_ext_idq_para =
        (; gens_vh,
         gens_i_d,
         gens_i_q,
         V_ref,
         Tm )

    ode_flux_decay_kwd_para = (
        ode_flux_decay_gens_para,
        ode_flux_decay_avrs_para,
        state_vars_idx,
        ωs)
    
    ode_flux_decay_state_sym =
        generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            state_vars_syms_flux_decay;
            label_prefix = "bus" )

    state_init_flux_decay =
        get_state_init_flux_decay(
            init_flux_decay_δ_eq_E_fd,
            ωs)
        
    #----------------------------------------
    #----------------------------------------
    
    ode_fun_flux_decay_by_ext_idq = ODEFunction(
        (dx,x,p,t) ->
            ode_flux_decay_by_ext_idq_func!(
                dx, x,
                ode_flux_decay_by_ext_idq_para, t;
                kwd_para =
                    ode_flux_decay_kwd_para);
        syms =
            ode_flux_decay_state_sym)
    
    ode_prob_flux_decay_by_ext_idq = ODEProblem(
        ode_fun_flux_decay_by_ext_idq,
        state_init_flux_decay,
        sim_timespan,
        ode_flux_decay_by_ext_idq_para )
    
    ode_sol_flux_decay_by_ext_idq =
        DifferentialEquations.solve(
            ode_prob_flux_decay_by_ext_idq,
            ode_alg,
            abstol = abstol,
            reltol = reltol )

    #----------------------------------------

    plot_flux_decay_by_ext_idq =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :eq_dash, :E_fd ],
            tspan = sim_timespan,
            fmt = :png)

    plt_layout = (2,2)

    vars_plots_flux_decay_by_ext_idq = [
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ;sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ ],
            tspan = sim_timespan,
            fmt = :png),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:ω ],
            tspan = sim_timespan,
            fmt = :png),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:eq_dash ],
            tspan = sim_timespan,
            fmt = :png),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:E_fd ],
            tspan = sim_timespan,
            fmt = :png)]
    
    plt_by_ext_idq = plot(
        vars_plots_flux_decay_by_ext_idq...;
        layout = plt_layout,
        size = (1000, 500),
        lw = 3,
        xlabel = "t[s]")

    return (; ode_sol_flux_decay_by_ext_idq,
            plt_by_ext_idq)

    """
    δ_plot_flux_decay_by_ext_idq =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ ],
            tspan = sim_timespan,
            fmt = :png)


    ω_plot_flux_decay_by_ext_idq =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:ω ],
            tspan = sim_timespan,
            fmt = :png)

    
    eq_dash_plot_flux_decay_by_ext_idq =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    
    E_fd_plot_flux_decay_by_ext_idq =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:E_fd ],
            tspan = sim_timespan,
            fmt = :png)

    """
    
end


function get_flux_decay_dynamic_sol_by_json(
    ; case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
        
    pf_alg        = NewtonRaphson(),
    ode_alg       = Rodas4(),
    
    abstol        = 1e-12,
    reltol        = 1e-12,
    
    dt            = 0.01,
    timespan      = 10.0 )

    #--------------------------------------    
    
    tspan         = (0.0, timespan)
    
    sim_timespan  = (0.0, timespan)
    
    plot_timespan = (0.0, timespan)
    
    #--------------------------------------    
    
    pf_solver =
        (; pf_alg,
         abstol,
         reltol)        

    use_nlsolve = false
    
    #--------------------------------------
    #--------------------------------------

    # ImplicitMidpoint(),
    # ImplicitMidpoint(autodiff=false)
    # ode_alg         = Rodas4()

    #--------------------------------------
    #--------------------------------------

    # flux_decay_dynamic_pf_sol_para =
    #     get_flux_decay_dynamic_pf_sol_para_by_json(
    #         ;case_name      = case_name,
    #         basekV          = 1.0,    
    #         use_pu_in_PQ    = true,
    #         line_data_in_pu = true,

    #         pf_alg          = NewtonRaphson(),
    #         ode_alg         = Rodas4(),

    #         abstol          = 1e-12,
    #         reltol          = 1e-12,

    #         dt              = 0.01,
    #         timespan        = 10.0,

    #         tspan           = (0.0, timespan),
    #         sim_timespan    = (0.0, timespan),
    #         plot_timespan   = (0.0, timespan))
    
    # #--------------------------------------    
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
    
    
    #----------------------------------------
    # parameters
    #----------------------------------------

    (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
     ode_gens_para,

     pf_sta_ΔPQ_mismatch_parameters,

     Ybr_cal_and_edges_orientation,
     sta_pf_PQ_para,
     ode_gens_generic_para,
     baseMVA,

     generic_govs_para,
     generic_avrs_para,
     ode_flux_decay_avrs_para) =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                     net_data_by_components_file;
                     
                     basekV = basekV,    
                     use_pu_in_PQ = use_pu_in_PQ,
                     line_data_in_pu = line_data_in_pu,
                     
                    in_components_type_sym =
                        false),
                 :net_generic_parameters),
             (:Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :ode_gens_para,

              :pf_sta_ΔPQ_mismatch_parameters,
              
              :Ybr_cal_and_edges_orientation,
              :sta_pf_PQ_para,
              :ode_gens_generic_para,
              :baseMVA,

              :generic_govs_para,
              :generic_avrs_para,
              :ode_flux_decay_avrs_para) )
    
    #----------------------------------------    
    # Indices
    #----------------------------------------    

    (loc_load_exist,
     dyn_pf_δ_eq_dash_0_P_Q_idx,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs,

     dyn_pf_δ_eq_dash_Idx,
     dyn_pf_Png_Qng_Pll_Qll_Idx,

     dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,

     dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,

     flux_decay_model_states_Idx,
     flux_decay_model_states_comp_idxs_in_Idx,

     flux_decay_model_vars_Idx_in_state,
     dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx) =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                     net_data_by_components_file;
                     
                     basekV = basekV,    
                     use_pu_in_PQ = use_pu_in_PQ,
                     line_data_in_pu = line_data_in_pu,
                     
                    in_components_type_sym =
                        false),
                 :net_generic_idx),
             (:loc_load_exist,
              :dyn_pf_δ_eq_dash_0_P_Q_idx,
              :Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              :non_pre_ordered_pf_vars_Idxs,

              :dyn_pf_δ_eq_dash_Idx,
              :dyn_pf_Png_Qng_Pll_Qll_Idx,

              :dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,

              :dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,

              :flux_decay_model_states_Idx,
              :flux_decay_model_states_comp_idxs_in_Idx,

              :flux_decay_model_vars_Idx_in_state,
              :dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx) )

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx) =  NamedTupleTools.select(
        dyn_pf_fun_kwd_net_idxs,
         (:gens_nodes_idx,
          :non_gens_nodes_idx,
          :all_nodes_idx)) 

    #----------------------------------------

    state_vars_syms_flux_decay =
        [:δ, :ω, :eq_dash, :E_fd ]
    
    ode_flux_decay_state_sym =
        generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            state_vars_syms_flux_decay;
            label_prefix = "bus" )

    pf_flux_decay_state_sym =
        [generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:vh];
            label_prefix = "bus");
         generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:θh];
            label_prefix = "bus")]

    flux_decay_model_syms =
        [ode_flux_decay_state_sym;
         pf_flux_decay_state_sym]
    
    #----------------------------------------
    # states in system_dynamics_flux_decay_model
    #----------------------------------------
    
    flux_decay_model_mass_matrix =
        DAE_MassMatrix(
            length(ode_flux_decay_state_sym),
            length(pf_flux_decay_state_sym) )
    
    #----------------------------------------
    
    dyn_pf_mismatch_vars_kwd_para =
        (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para )

    #----------------------------------------
    
    dyn_pf_mismatch_idx_kwd_para =
        (;dyn_pf_δ_eq_dash_0_P_Q_idx,
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         non_pre_ordered_pf_vars_Idxs,
         dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
         dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx)
    
    #----------------------------------------

    dyn_pf_mismatch_kwd_para =
        (; loc_load_exist,
         dyn_pf_mismatch_vars_kwd_para,
         dyn_pf_mismatch_idx_kwd_para )

    #----------------------------------------

    sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx =
        (;dyn_pf_δ_eq_dash_Idx,
         dyn_pf_Png_Qng_Pll_Qll_Idx)

    flux_decay_model_states_Idx_wt_idxs_in_Idx =
        (;flux_decay_model_states_Idx,
         flux_decay_model_states_comp_idxs_in_Idx,
         flux_decay_model_vars_Idx_in_state )

    #----------------------------------------    
    #----------------------------------------
    
    by_part_dyn_pf_mismatch_kwd_para =
        (;pf_solver,
         use_nlsolve,
         dyn_pf_mismatch_kwd_para,
         sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx )

    #----------------------------------------

    (pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
         NamedTupleTools.select(
             pf_sta_ΔPQ_mismatch_parameters,
             (:pf_kw_para,
              :red_types_Idxs_etc,
              :pf_PQ_param) )

    #----------------------------------------
    
    (red_vh_Idxs,
     red_θh_Idxs) =
         NamedTupleTools.select(
             red_types_Idxs_etc,
             (:red_vh_Idxs,
              :red_θh_Idxs) )
    
    #----------------------------------------
    
    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs))]

    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
          
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------
    
    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (; Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_results_pf_sta_red_sol =
        get_generic_results_pf_sta_red_sol_u(
            pf_sol;
            generic_red_sol_kwd_para =
                generic_red_sol_kwd_para,
            baseMVA =
                baseMVA,
            basekV =
                1.0 )

    #----------------------------------------    
    #----------------------------------------    
    
    (pf_P_gens, pf_Q_gens,
     pf_P_g_gens, pf_Q_g_gens,
     Igen,
     vh, θh, θh_deg,
     gens_vh, gens_θh,
     gens_nodes_idx) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens, :pf_Q_gens,
             :pf_P_g_gens, :pf_Q_g_gens,
             :Igen,
             :vh, :θh, :θh_deg,
             :gens_vh, :gens_θh,
             :gens_nodes_idx ) )
        
    #----------------------------------------
    # Dynamic simulation
    #----------------------------------------

    (X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            ode_gens_generic_para,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash ) )
    
    #----------------------------------------
    
    (Ka,
     Ta) = ode_flux_decay_avrs_para

    #----------------------------------------
    
    nt_init_flux_decay =
        get_init_flux_decay_model(
            gens_vh,
            gens_θh,            
            pf_P_g_gens,
            pf_Q_g_gens,            
            X_d,
            X_q,            
            X_d_dash,
            X_q_dash,
            Ka)
    
    #----------------------------------------
    
    init_flux_decay_vars =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:δ,
                          :eq_dash, :E_fd,
                          :V_ref, :Tm ) )
    
    init_flux_decay_V_ref_Tm =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm ) )
    
    init_flux_decay_δ_eq_E_fd =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:δ, :eq_dash, :E_fd ) )

   init_flux_decay_i_d_i_q =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:i_d, :i_q ) )

    #----------------------------------------

    (V_ref,
     Tm,
     gens_i_d,
     gens_i_q) =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm, :i_d, :i_q ))
    
    #----------------------------------------
    
    ode_flux_decay_gens_para =
        NamedTupleTools.select(
            ode_gens_generic_para,
                (:H,  :X_d, :X_q, :X_d_dash,
                 :T_d_dash ) )
        
    state_vars_idx =
     get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
            state_vars_syms_flux_decay,
            gens_nodes_idx )

    #----------------------------------------

    ode_flux_decay_para =
        ( ;gens_vh,
          gens_θh,
          V_ref,
          Tm)

    ode_flux_decay_kwd_para = (;
        ode_flux_decay_gens_para,
        ode_flux_decay_avrs_para,
        state_vars_idx,
        ωs)
     
    #----------------------------------------    
    # ode states init
    #----------------------------------------    
    
    state_init_flux_decay =
        get_state_init_flux_decay(
            init_flux_decay_δ_eq_E_fd,
            ωs)
    
    #----------------------------------------    
    
    (dyn_δ,
     dyn_eq_dash) =
        NamedTupleTools.select(
            init_flux_decay_δ_eq_E_fd,
            (:δ,
             :eq_dash ) )

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

    if loc_load_exist == true

        flux_decay_dyn_pf_flat_para =
            Float64[dyn_δ;
             dyn_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load ]
        
    else

        flux_decay_dyn_pf_flat_para =
            Float64[dyn_δ;
             dyn_eq_dash;
             P_non_gens;
             Q_non_gens ]
        
    end

    # -----------------------------------

    flux_decay_dyn_pf_δ_eq_dash_flat_para =
        Float64[dyn_δ;
         dyn_eq_dash]
    
    flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para =
        Float64[P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load ]

    flux_decay_dyn_pf_Pg_Qg_Png_Qng_Pgll_Qgll_flat_para =
        Float64[pf_P_g_gens;
         pf_Q_g_gens;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load ]

    #----------------------------------------    
        
    (;flat_vh_Idx,
     flat_θh_Idx) =
        NamedTupleTools.select(
            getproperty(
                non_pre_ordered_pf_vars_Idxs,
                :flat_vh_flat_θh_flat_id_iq_Idx) ,
               (:flat_vh_Idx,
                :flat_θh_Idx))

    #----------------------------------------    

    init_flat_vh_flat_θh = [vh; θh]
    
    #----------------------------------------    
    
    # pf_sol_post_state_init =
    #     v2_get_flux_decay_model_pf_wt_ext_idq_by_parts_func!(
    #         init_flat_vh_flat_θh,
    #         flux_decay_dyn_pf_δ_eq_dash_flat_para,
    #         flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para;
    #         by_part_dyn_pf_mismatch_kwd_para =
    #             by_part_dyn_pf_mismatch_kwd_para )

    # pf_sol_post_state_init_vh =
    #     pf_sol_post_state_init[
    #         flat_vh_Idx]
    
    # pf_sol_post_state_init_θh =
    #     pf_sol_post_state_init[
    #         flat_θh_Idx]
    
    # pf_state_init_flux_decay = [
    #     pf_sol_post_state_init_vh;
    #     pf_sol_post_state_init_θh ]
    
    #----------------------------------------    
    # pf states init
    #----------------------------------------    


    pf_state_init_flux_decay = [
        vh;
        θh ]
    
    #----------------------------------------
    # states in system_dynamics_flux_decay_model
    #----------------------------------------
    
    """

    flux_decay_model_states_init =
        [state_init_flux_decay;
         init_flat_vh_flat_θh ]

    """

    flux_decay_model_states_init =
        [state_init_flux_decay;
         pf_state_init_flux_decay ]

    #----------------------------------------
    # flux_decay_model_dynamics_para
    #----------------------------------------

   # https://docs.sciml.ai/ComponentArrays/stable/quickstart/
    
    ComponentVector_flux_decay_model_dynamics_para =
        ComponentVector{Float64}(
            V_ref = V_ref,
            Tm = Tm,
            Png = P_non_gens,
            Qng = Q_non_gens, 
            Pll = P_g_loc_load,
            Qll = Q_g_loc_load)

    flux_decay_model_dynamics_para =
        Float64[V_ref;
                Tm;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]
    
    #----------------------------------------
    # PreallocationTools.jl
    # https://github.com/SciML/PreallocationTools.jl
    #----------------------------------------
    
    stateDiffCache = DiffCache(
        flux_decay_model_states_init)
    
    #----------------------------------------
    # flux_decay_model_dynamics_kwd_para
    #----------------------------------------
    
    flux_decay_model_dynamics_kwd_para =
        (;dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
         stateDiffCache,
         ode_flux_decay_kwd_para,
         by_part_dyn_pf_mismatch_kwd_para,
         sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx,
         flux_decay_model_states_Idx_wt_idxs_in_Idx )
    
    #----------------------------------------
    # system dyamanicl flux decay simulation
    #----------------------------------------    

    # flux_decay_model_dynamics!
    
    DAE_dynamics_fun! = v2_flux_decay_model_dynamics!
    
    model_dynamics_fun = ODEFunction(
        (dx,x,p,t) ->
            DAE_dynamics_fun!(
                dx, x,
                flux_decay_model_dynamics_para, t;
                kwd_para =
                    flux_decay_model_dynamics_kwd_para);
        syms =
            flux_decay_model_syms,
        mass_matrix =
            flux_decay_model_mass_matrix)
    
    model_prob = ODEProblem(
        model_dynamics_fun,
        flux_decay_model_states_init,
        sim_timespan,
        flux_decay_model_dynamics_para )
    
    model_sol =
        DifferentialEquations.solve(
            model_prob,
            ode_alg,
            abstol = abstol,
            reltol = reltol )

    
    # model_sol =
    #     DifferentialEquations.solve(
    #         model_prob,
    #         ImplicitMidpoint(),
    #         dt = 0.01,
    #         abstol = abstol,
    #         reltol = reltol )
    
    
end


function get_flux_decay_ode_and_system_dynamic_sol_by_json(
    ; case_name = "case9",
    components_libs_dir="",
    json_net_data_by_components_file ="",
    data_dir = "",    
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,

        
    pf_alg        = NewtonRaphson(),
    ode_alg       = Rodas4(),
    
    abstol        = 1e-12,
    reltol        = 1e-12,
    
    dt            = 0.01,
    timespan      = 10.0 )

    #--------------------------------------    
    
    tspan         = (0.0, timespan)
    
    sim_timespan  = (0.0, timespan)
    
    plot_timespan = (0.0, timespan)
    
    #--------------------------------------    
    
    pf_solver =
        (; pf_alg,
         abstol,
         reltol)        

    use_nlsolve = false

    #--------------------------------------
    #--------------------------------------

    # ImplicitMidpoint(),
    # ImplicitMidpoint(autodiff=false)
    # ode_alg         = Rodas4()

    #--------------------------------------

    # flux_decay_dynamic_pf_sol_para =
    #     get_flux_decay_dynamic_pf_sol_para_by_json(
    #         ;case_name      = case_name,
    #         basekV          = 1.0,    
    #         use_pu_in_PQ    = true,
    #         line_data_in_pu = true,

    #         pf_alg          = NewtonRaphson(),
    #         ode_alg         = Rodas4(),

    #         abstol          = 1e-12,
    #         reltol          = 1e-12,

    #         dt              = 0.01,
    #         timespan        = 10.0,

    #         tspan           = (0.0, timespan),
    #         sim_timespan    = (0.0, timespan),
    #         plot_timespan   = (0.0, timespan))
    

    
    #--------------------------------------    
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")
    
    if (json_net_data_by_components_file == "") || (
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
            
    #----------------------------------------
    # parameters
    #----------------------------------------

    (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
     ode_gens_para,

     pf_sta_ΔPQ_mismatch_parameters,

     Ybr_cal_and_edges_orientation,
     sta_pf_PQ_para,
     ode_gens_generic_para,
     baseMVA,

     generic_govs_para,
     generic_avrs_para,
     ode_flux_decay_avrs_para) =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                     net_data_by_components_file;
                     
                     basekV = basekV,    
                     use_pu_in_PQ = use_pu_in_PQ,
                     line_data_in_pu = line_data_in_pu,
                     
                    in_components_type_sym =
                        false),
                 :net_generic_parameters),
             (:Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :ode_gens_para,

              :pf_sta_ΔPQ_mismatch_parameters,
              
              :Ybr_cal_and_edges_orientation,
              :sta_pf_PQ_para,
              :ode_gens_generic_para,
              :baseMVA,

              :generic_govs_para,
              :generic_avrs_para,
              :ode_flux_decay_avrs_para) )
    
    #----------------------------------------    
    # Indices
    #----------------------------------------    

    (loc_load_exist,
     dyn_pf_δ_eq_dash_0_P_Q_idx,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs,

     dyn_pf_δ_eq_dash_Idx,
     dyn_pf_Png_Qng_Pll_Qll_Idx,

     dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,

     dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,

     flux_decay_model_states_Idx,
     flux_decay_model_states_comp_idxs_in_Idx,

     flux_decay_model_vars_Idx_in_state,
     dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx) =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                     net_data_by_components_file;
                     
                     basekV = basekV,    
                     use_pu_in_PQ = use_pu_in_PQ,
                     line_data_in_pu = line_data_in_pu,
                     
                    in_components_type_sym =
                        false),
                 :net_generic_idx),
             (:loc_load_exist,
              :dyn_pf_δ_eq_dash_0_P_Q_idx,
              :Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              :non_pre_ordered_pf_vars_Idxs,

              :dyn_pf_δ_eq_dash_Idx,
              :dyn_pf_Png_Qng_Pll_Qll_Idx,

              :dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,

              :dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,

              :flux_decay_model_states_Idx,
              :flux_decay_model_states_comp_idxs_in_Idx,

              :flux_decay_model_vars_Idx_in_state,
              :dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx) )

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx) =  NamedTupleTools.select(
        dyn_pf_fun_kwd_net_idxs,
         (:gens_nodes_idx,
          :non_gens_nodes_idx,
          :all_nodes_idx)) 

    #----------------------------------------

    state_vars_syms_flux_decay =
        [:δ, :ω, :eq_dash, :E_fd ]
    
    ode_flux_decay_state_sym =
        generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            state_vars_syms_flux_decay;
            label_prefix = "bus" )

    pf_flux_decay_state_sym =
        [generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:vh];
            label_prefix = "bus");
         generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:θh];
            label_prefix = "bus")]

    flux_decay_model_syms =
        [ode_flux_decay_state_sym;
         pf_flux_decay_state_sym]
    
    #----------------------------------------
    # states in system_dynamics_flux_decay_model
    #----------------------------------------
    
    flux_decay_model_mass_matrix =
        DAE_MassMatrix(
            length(ode_flux_decay_state_sym),
            length(pf_flux_decay_state_sym) )
    
    #----------------------------------------
    
    dyn_pf_mismatch_vars_kwd_para =
        (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para )

    #----------------------------------------
    
    dyn_pf_mismatch_idx_kwd_para =
        (;dyn_pf_δ_eq_dash_0_P_Q_idx,
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         non_pre_ordered_pf_vars_Idxs,
         dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
         dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx)
    
    #----------------------------------------

    dyn_pf_mismatch_kwd_para =
        (; loc_load_exist,
         dyn_pf_mismatch_vars_kwd_para,
         dyn_pf_mismatch_idx_kwd_para )

    #----------------------------------------

    sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx =
        (;dyn_pf_δ_eq_dash_Idx,
         dyn_pf_Png_Qng_Pll_Qll_Idx)

    flux_decay_model_states_Idx_wt_idxs_in_Idx =
        (;flux_decay_model_states_Idx,
         flux_decay_model_states_comp_idxs_in_Idx,
         flux_decay_model_vars_Idx_in_state )

    #----------------------------------------    
    #----------------------------------------
    
    by_part_dyn_pf_mismatch_kwd_para =
        (;pf_solver,
         use_nlsolve,
         dyn_pf_mismatch_kwd_para,
         sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx )

    #----------------------------------------

    (pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
         NamedTupleTools.select(
             pf_sta_ΔPQ_mismatch_parameters,
             (:pf_kw_para,
              :red_types_Idxs_etc,
              :pf_PQ_param) )

    #----------------------------------------
    
    (red_vh_Idxs,
     red_θh_Idxs) =
         NamedTupleTools.select(
             red_types_Idxs_etc,
             (:red_vh_Idxs,
              :red_θh_Idxs) )
    
    #----------------------------------------
    
    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs))]

    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
          
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------
    
    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (; Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_results_pf_sta_red_sol =
        get_generic_results_pf_sta_red_sol_u(
            pf_sol;
            generic_red_sol_kwd_para =
                generic_red_sol_kwd_para,
            baseMVA =
                baseMVA,
            basekV =
                1.0 )

    #----------------------------------------    
    #----------------------------------------    
    
    (pf_P_gens, pf_Q_gens,
     pf_P_g_gens, pf_Q_g_gens,
     Igen,
     vh, θh, θh_deg,
     gens_vh, gens_θh,
     gens_nodes_idx) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens, :pf_Q_gens,
             :pf_P_g_gens, :pf_Q_g_gens,
             :Igen,
             :vh, :θh, :θh_deg,
             :gens_vh, :gens_θh,
             :gens_nodes_idx ) )
        
    #----------------------------------------
    # Dynamic simulation
    #----------------------------------------

    (X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            ode_gens_generic_para,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash ) )
    
    #----------------------------------------
    
    (Ka,
     Ta) = ode_flux_decay_avrs_para

    #----------------------------------------
    
    nt_init_flux_decay =
        get_init_flux_decay_model(
            gens_vh,
            gens_θh,            
            pf_P_g_gens,
            pf_Q_g_gens,            
            X_d,
            X_q,            
            X_d_dash,
            X_q_dash,
            Ka)
    
    #----------------------------------------
    
    init_flux_decay_vars =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:δ,
                          :eq_dash, :E_fd,
                          :V_ref, :Tm ) )
    
    init_flux_decay_V_ref_Tm =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm ) )
    
    init_flux_decay_δ_eq_E_fd =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:δ, :eq_dash, :E_fd ) )

   init_flux_decay_i_d_i_q =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:i_d, :i_q ) )

    #----------------------------------------

    (V_ref,
     Tm,
     gens_i_d,
     gens_i_q) =
        get_selected_vec_nt_to_vec_vec(
            nt_init_flux_decay, nothing;
            selections = (:V_ref, :Tm, :i_d, :i_q ))
    
    #----------------------------------------
    
    ode_flux_decay_gens_para =
        NamedTupleTools.select(
            ode_gens_generic_para,
                (:H,  :X_d, :X_q, :X_d_dash,
                 :T_d_dash ) )
        
    state_vars_idx =
     get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
            state_vars_syms_flux_decay,
            gens_nodes_idx )
    
    #----------------------------------------
    
    ode_flux_decay_para =
        ( ;gens_vh,
          gens_θh,
          V_ref,
          Tm)
    
    ode_flux_decay_by_ext_idq_para =
        (; gens_vh,
         gens_i_d,
         gens_i_q,
         V_ref,
         Tm )

    ode_flux_decay_kwd_para = (;
        ode_flux_decay_gens_para,
        ode_flux_decay_avrs_para,
        state_vars_idx,
        ωs)
     
    #----------------------------------------    
    # ode states init
    #----------------------------------------    
    
    state_init_flux_decay =
        get_state_init_flux_decay(
            init_flux_decay_δ_eq_E_fd,
            ωs)
    
    #----------------------------------------    
    #----------------------------------------    
    
    (dyn_δ,
     dyn_eq_dash) =
        NamedTupleTools.select(
            init_flux_decay_δ_eq_E_fd,
            (:δ,
             :eq_dash ) )

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

    if loc_load_exist == true

        flux_decay_dyn_pf_flat_para =
            Float64[dyn_δ;
             dyn_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load ]
        
    else

        flux_decay_dyn_pf_flat_para =
            Float64[dyn_δ;
             dyn_eq_dash;
             P_non_gens;
             Q_non_gens ]
        
    end

    # -----------------------------------

    flux_decay_dyn_pf_δ_eq_dash_flat_para =
        Float64[dyn_δ;
         dyn_eq_dash]
    
    flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para =
        Float64[P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load ]

    flux_decay_dyn_pf_Pg_Qg_Png_Qng_Pgll_Qgll_flat_para =
        Float64[pf_P_g_gens;
         pf_Q_g_gens;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load ]

    #----------------------------------------    

    init_flat_vh_flat_θh = [vh; θh]

    #----------------------------------------    
    
    pf_sol_post_state_init =
        get_flux_decay_model_pf_wt_ext_idq_by_parts_func!(
            init_flat_vh_flat_θh,
            flux_decay_dyn_pf_δ_eq_dash_flat_para,
            flux_decay_dyn_pf_Png_Qng_Pll_Qll_flat_para;
            by_part_dyn_pf_mismatch_kwd_para =
                by_part_dyn_pf_mismatch_kwd_para )
    
    #----------------------------------------    
        
    (;flat_vh_Idx,
     flat_θh_Idx) =
        NamedTupleTools.select(
            getproperty(
                non_pre_ordered_pf_vars_Idxs,
                :flat_vh_flat_θh_flat_id_iq_Idx) ,
               (:flat_vh_Idx,
                :flat_θh_Idx))
    
    #----------------------------------------    
    # pf states init
    #----------------------------------------    

    pf_sol_post_state_init_vh =
        pf_sol_post_state_init[
            flat_vh_Idx]
    
    pf_sol_post_state_init_θh =
        pf_sol_post_state_init[
            flat_θh_Idx]

    #----------------------------------------
    
    # pf_state_init_flux_decay = [
    #     pf_sol_post_state_init_vh;
    #     pf_sol_post_state_init_θh ]


    pf_state_init_flux_decay = [
        vh;
        θh ]
    
    #----------------------------------------
    # states in system_dynamics_flux_decay_model
    #----------------------------------------
    
    """

    flux_decay_model_states_init =
        [state_init_flux_decay;
         init_flat_vh_flat_θh ]

    """

    flux_decay_model_states_init =
        [state_init_flux_decay;
         pf_state_init_flux_decay ]

    #----------------------------------------
    # flux_decay_model_dynamics_para
    #----------------------------------------

   # https://docs.sciml.ai/ComponentArrays/stable/quickstart/
    
    ComponentVector_flux_decay_model_dynamics_para =
        ComponentVector{Float64}(
            V_ref = V_ref,
            Tm = Tm,
            Png = P_non_gens,
            Qng = Q_non_gens, 
            Pll = P_g_loc_load,
            Qll = Q_g_loc_load)

    flux_decay_model_dynamics_para =
        Float64[V_ref;
                Tm;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]
    
    #----------------------------------------
    # PreallocationTools.jl
    # https://github.com/SciML/PreallocationTools.jl
    #----------------------------------------
    
    stateDiffCache = DiffCache(
        flux_decay_model_states_init)
    
    #----------------------------------------
    # flux_decay_model_dynamics_kwd_para
    #----------------------------------------
    
    flux_decay_model_dynamics_kwd_para =
        (;dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
         stateDiffCache,
         ode_flux_decay_kwd_para,
         by_part_dyn_pf_mismatch_kwd_para,
         sep_dyn_pf_δ_eq_dash_wt_PQng_PQll_Idx,
         flux_decay_model_states_Idx_wt_idxs_in_Idx )
    
    #----------------------------------------
    # ode flux decay model simulation
    #----------------------------------------    
    
    ode_fun_flux_decay = ODEFunction(
        (dx,x,p,t) ->
            ode_flux_decay!(
                dx, x,
                ode_flux_decay_para, t;
                kwd_para =
                    ode_flux_decay_kwd_para);
        syms =
            ode_flux_decay_state_sym)
    
    ode_prob_flux_decay = ODEProblem(
        ode_fun_flux_decay,
        state_init_flux_decay,
        sim_timespan,
        ode_flux_decay_para )
    
    ode_sol_flux_decay =
        DifferentialEquations.solve(
            ode_prob_flux_decay,
            ode_alg,
            abstol = abstol,
            reltol = reltol )

    #----------------------------------------

    plot_flux_decay =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :eq_dash, :E_fd ],
            tspan = sim_timespan,
            fmt = :png)

    plt_layout = (2,2)

    vars_plots_flux_decay = [
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ;sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ ],
            tspan = sim_timespan,
            fmt = :png),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:ω ],
            tspan = sim_timespan,
            fmt = :png),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:eq_dash ],
            tspan = sim_timespan,
            fmt = :png),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:E_fd ],
            tspan = sim_timespan,
            fmt = :png)]
    
    plt = plot(
        vars_plots_flux_decay...;
        layout = plt_layout,
        size = (1000, 500),
        lw = 3,
        xlabel = "t[s]")

    
    # δ_plot_flux_decay =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ; sol =
    #             ode_sol_flux_decay,
    #         network_vars_labels =
    #             ode_flux_decay_state_sym,
    #         nodes_name = ["bus1", "bus2", "bus3"],
    #         vars = [:δ ],
    #         tspan = sim_timespan,
    #         fmt = :png)


    # ω_plot_flux_decay =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ; sol =
    #             ode_sol_flux_decay,
    #         network_vars_labels =
    #             ode_flux_decay_state_sym,
    #         nodes_name = ["bus1", "bus2", "bus3"],
    #         vars = [:ω ],
    #         tspan = sim_timespan,
    #         fmt = :png)

    
    # eq_dash_plot_flux_decay =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ; sol =
    #             ode_sol_flux_decay,
    #         network_vars_labels =
    #             ode_flux_decay_state_sym,
    #         nodes_name = ["bus1", "bus2", "bus3"],
    #         vars = [:eq_dash ],
    #         tspan = sim_timespan,
    #         fmt = :png)

    
    # E_fd_plot_flux_decay =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ; sol =
    #             ode_sol_flux_decay,
    #         network_vars_labels =
    #             ode_flux_decay_state_sym,
    #         nodes_name = ["bus1", "bus2", "bus3"],
    #         vars = [:E_fd ],
    #         tspan = sim_timespan,
    #         fmt = :png)
    
    #----------------------------------------
    # ode flux decay model with ext id_iq simulation
    #----------------------------------------    
    
    ode_fun_flux_decay_by_ext_idq = ODEFunction(
        (dx,x,p,t) ->
            ode_flux_decay_by_ext_idq_func!(
                dx, x,
                ode_flux_decay_by_ext_idq_para, t;
                kwd_para =
                    ode_flux_decay_kwd_para);
        syms =
            ode_flux_decay_state_sym)
    
    ode_prob_flux_decay_by_ext_idq = ODEProblem(
        ode_fun_flux_decay_by_ext_idq,
        state_init_flux_decay,
        sim_timespan,
        ode_flux_decay_by_ext_idq_para )
    
    ode_sol_flux_decay_by_ext_idq =
        DifferentialEquations.solve(
            ode_prob_flux_decay_by_ext_idq,
            ode_alg,
            abstol = abstol,
            reltol = reltol )

    #----------------------------------------

    plot_flux_decay_by_ext_idq =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :eq_dash, :E_fd ],
            tspan = sim_timespan,
            fmt = :png)

    plt_layout = (2,2)

    vars_plots_flux_decay_by_ext_idq = [
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ;sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ ],
            tspan = sim_timespan,
            fmt = :png),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:ω ],
            tspan = sim_timespan,
            fmt = :png),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:eq_dash ],
            tspan = sim_timespan,
            fmt = :png),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                ode_sol_flux_decay_by_ext_idq,
            network_vars_labels =
                ode_flux_decay_state_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:E_fd ],
            tspan = sim_timespan,
            fmt = :png)]
    
    plt_by_ext_idq = plot(
        vars_plots_flux_decay_by_ext_idq...;
        layout = plt_layout,
        size = (1000, 500),
        lw = 3,
        xlabel = "t[s]")

    # δ_plot_flux_decay_by_ext_idq =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ; sol =
    #             ode_sol_flux_decay_by_ext_idq,
    #         network_vars_labels =
    #             ode_flux_decay_state_sym,
    #         nodes_name = ["bus1", "bus2", "bus3"],
    #         vars = [:δ ],
    #         tspan = sim_timespan,
    #         fmt = :png)


    # ω_plot_flux_decay_by_ext_idq =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ; sol =
    #             ode_sol_flux_decay_by_ext_idq,
    #         network_vars_labels =
    #             ode_flux_decay_state_sym,
    #         nodes_name = ["bus1", "bus2", "bus3"],
    #         vars = [:ω ],
    #         tspan = sim_timespan,
    #         fmt = :png)

    
    # eq_dash_plot_flux_decay_by_ext_idq =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ; sol =
    #             ode_sol_flux_decay_by_ext_idq,
    #         network_vars_labels =
    #             ode_flux_decay_state_sym,
    #         nodes_name = ["bus1", "bus2", "bus3"],
    #         vars = [:eq_dash ],
    #         tspan = sim_timespan,
    #         fmt = :png)

    
    # E_fd_plot_flux_decay_by_ext_idq =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ; sol =
    #             ode_sol_flux_decay_by_ext_idq,
    #         network_vars_labels =
    #             ode_flux_decay_state_sym,
    #         nodes_name = ["bus1", "bus2", "bus3"],
    #         vars = [:E_fd ],
    #         tspan = sim_timespan,
    #         fmt = :png)

    #----------------------------------------
    # system dyamanicl flux decay simulation
    #----------------------------------------    

    # flux_decay_model_dynamics!
    
    DAE_dynamics_fun! = v2_flux_decay_model_dynamics!
    
    model_dynamics_fun = ODEFunction(
        (dx,x,p,t) ->
            DAE_dynamics_fun!(
                dx, x,
                flux_decay_model_dynamics_para, t;
                kwd_para =
                    flux_decay_model_dynamics_kwd_para);
        syms =
            flux_decay_model_syms,
        mass_matrix =
            flux_decay_model_mass_matrix)
    
    model_prob = ODEProblem(
        model_dynamics_fun,
        flux_decay_model_states_init,
        sim_timespan,
        flux_decay_model_dynamics_para )
    
    model_sol =
        DifferentialEquations.solve(
            model_prob,
            ode_alg,
            abstol = abstol,
            reltol = reltol )
    
    # model_sol =
    #     DifferentialEquations.solve(
    #         model_prob,
    #         ImplicitMidpoint(),
    #         dt = 0.01,
    #         abstol = abstol,
    #         reltol = reltol )
    
    
end

