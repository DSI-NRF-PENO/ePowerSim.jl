# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247


#####################################################
# ---------------------------------------------------
#  Plants init
# ---------------------------------------------------
#####################################################


function plants_generic_model_init_func(
    gens_vh,
    gens_θh,
    pf_P_gens,
    pf_Q_gens,
    ωs;
    kwd_para =
        plants_init_kwd_para )

    (gens_init_para,
     avrs_init_para,
     govs_init_para,
     comps_init_funs,
     ) =
         kwd_para
             
    gens_init_fun =
        [an_item.gen_init_fun for an_item in
             comps_init_funs]

    avrs_init_fun = 
        [an_item.avr_init_fun for an_item in
             comps_init_funs]

    govs_init_fun = 
        [an_item.gov_init_fun for an_item in
             comps_init_funs]

    states_init_wt_ref =
        [  Symbol(split(String(
            nameof(gov_init_fun)),
                        "__")[2]) == :nothing  ?
            a_SC_plant_generic_model_init_func(
                    gen_vh, gen_θh,
                    pf_P_gen, pf_Q_gen,
                    ωs;
                    kwd_para =
                        (gen_init_para,
                         avr_init_para,
                         (gen_init_fun,
                          avr_init_fun ) ) ) :
                        a_SM_plant_generic_model_init_func(
                            gen_vh, gen_θh,
                            pf_P_gen, pf_Q_gen,
                            ωs;
                            kwd_para =
                                (gen_init_para,
                                 avr_init_para,
                                 gov_init_para,
                                 (gen_init_fun,
                                 avr_init_fun,
                                  gov_init_fun)  ))
          for (gen_vh, gen_θh,
             pf_P_gen, pf_Q_gen,
               gen_init_para, avr_init_para,
               gov_init_para,
               gen_init_fun, avr_init_fun,
               gov_init_fun ) in
              zip( gens_vh, gens_θh,
                   
                   pf_P_gens, pf_Q_gens,
                   
                   gens_init_para,
                   avrs_init_para,
                   govs_init_para,
                   
                   gens_init_fun,
                   avrs_init_fun,
                   govs_init_fun) ]

    #----------------------------------------
    #----------------------------------------
    
    # init_states = first.( states_init_wt_ref )

    # plants_states_init = Float64[init_states...;]
    
    # plants_refs = ( second.( states_init_wt_ref ))
    
    plants_states_init =
        Float64[ [a_plant_init_states.plant_states_init
                for a_plant_init_states in
                    states_init_wt_ref]...; ]
    
    plants_refs = Tuple( [a_plant_init_states.plant_ref
                for a_plant_init_states in
                    states_init_wt_ref] )
    
    plants_states_init_per_plant =
        [a_plant_init_states.plant_states_init
                for a_plant_init_states in
                    states_init_wt_ref]
    
    plants_refs_per_plant =
        [a_plant_init_states.plant_ref
                for a_plant_init_states in
                    states_init_wt_ref] 
    
    #----------------------------------------

    return (;plants_states_init,
            plants_refs,
            plants_states_init_per_plant,
            plants_refs_per_plant)
    
end

# ---------------------------------------------------
# plants init
# ---------------------------------------------------

function get_generic_system_dynamics_init_and_refs(
    vh,
    θh,
    pf_P_gens,
    pf_Q_gens;
    kwd_para =
        generic_network_init_kwd_para)

    (;gens_nodes_idx,
     gens_state_vars_idx_in_state,
     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,
     comps_init_funs,
     ωs) =
         kwd_para

    #----------------------------------------

    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))

    #----------------------------------------
    
    gens_vh = vh[gens_nodes_idx]
    
    gens_θh = θh[gens_nodes_idx]

    #----------------------------------------
    # Init
    #----------------------------------------
    
    plants_init_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para ,
         comps_init_funs )
    
    plants_generic_states_init_wt_ref =
        plants_generic_model_init_func(
            gens_vh,
            gens_θh,
            pf_P_gens,
            pf_Q_gens,
            ωs;
            kwd_para =
                plants_init_kwd_para )

    (plants_states_init,
     plants_states_init_per_plant,
     plants_refs ) =
         NamedTupleTools.select(
             plants_generic_states_init_wt_ref,
             (:plants_states_init,
              :plants_states_init_per_plant,
              :plants_refs))


    ( nt_vec_per_paras,
      vec_vec_per_paras ) =
          get_nt_vec_wt_vec_vec_per_paras(
        plants_refs ;
        nt_syms =
            (:ωs,
             :ω_ref,
             :v_ref,
             :p_order,
             :i_d,
             :i_q,
             :gen_vh ) )

    (ω_ref,
     v_ref,
     p_order,
     gens_i_d,
     gens_i_q ) =
         NamedTupleTools.select(
             nt_vec_per_paras, (
                 :ω_ref,
                 :v_ref,
                 :p_order,
                 :i_d,
                 :i_q))

    #----------------------------------------
    # System model init
    #----------------------------------------
    
    system_init =
        Float64[plants_states_init;
                vh;
                θh;
                gens_i_d;
                gens_i_q]

    #----------------------------------------
    
    δ = system_init[
        δ_idx_in_state]
        
    ed_dash = system_init[
        ed_dash_idx_in_state]
    
    eq_dash = system_init[
        eq_dash_idx_in_state]
    
    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh]

    #----------------------------------------
    
    ω_ref_v_ref_p_order =
        Float64[ω_ref;
                v_ref;
                p_order]

    #----------------------------------------
    
    return (;system_init,
            plants_states_init_per_plant,
            
            δ,
            ed_dash,
            eq_dash,
            
            ω_ref,
            v_ref,
            p_order,
            
            gens_i_d,
            gens_i_q,
            gens_vh )
 
end

# ---------------------------------------------------
# ---------------------------------------------------

function get_init_a_gen_full_model(
    gen_vh,
    gen_θh,
    P_g,
    Q_g,
    gen_para,
    avr_para,
    # govs_para,    
    ωs
    )

    (ra, X_d, X_q, X_d_dash, X_q_dash) =
        gen_para

    (Ka, Ke, Kf, Ta, Te, Tf, Ae, Be ) =
        avr_para
    
    # Sauer Page 187,    
    
    Ig = (P_g - im * Q_g) / (
        gen_vh * exp(-im * gen_θh))

    E_gen = gen_vh * exp(im * gen_θh) +
        (ra + im * X_q) * Ig

    δ = angle( E_gen )

    E = abs(E_gen)

    i_dq =  abs(Ig) * exp(im * (angle(Ig) - δ + π/2) )

    v_dq =  gen_vh * exp(im * (gen_θh - δ + π/2) )

    i_d = real(i_dq )

    i_q = imag(i_dq )

    v_d = real(v_dq )

    v_q = imag(v_dq )

    #   ed_dash = (X_q - X_q_dash) * i_q
    
    ed_dash = v_d  + ra * i_d - X_q_dash * i_q
    
    eq_dash = v_q  + ra * i_q + X_d_dash * i_d

    E_fd = eq_dash + (X_d - X_d_dash) * i_d

    
    V_R = (Ke + Sevf(Ae, Be,  E_fd)) * E_fd

    R_f = Kf/Tf * E_fd

    V_ref =  gen_vh + (V_R / Ka)
    
    
    Tm = ed_dash * i_d + eq_dash * i_q +
        (X_q_dash - X_d_dash) * i_d * i_q

    return (; δ, ed_dash, eq_dash,
            E_fd, V_R, R_f, V_ref,
            Tm,
            v_d, v_q,
            i_d, i_q)

end


function get_init_gens_full_model(
    gens_vh,
    gens_θh,
    gens_P_g,
    gens_Q_g,
    gens_para,
    # govs_para,
    avrs_para,
    ωs )
    
    (;ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash) =  NamedTupleTools.select(
         gens_para,
         (:ra,
          :X_d,
          :X_q,
          :X_d_dash,
          :X_q_dash))

    (;Ka , Ke, Kf,
     Ta, Te, Tf, Ae, Be ) =
        NamedTupleTools.select(
            avrs_para,
            (:Ka,:Ke, :Kf,
             :Ta, :Te, :Tf,
             :Ae, :Be))
    
    return [ get_init_a_gen_full_model(
        a_vh, a_θh, a_P_g, a_Q_g,
        (a_ra, a_X_d, a_X_q, a_X_d_dash, a_X_q_dash),
        (a_Ka, a_Ke, a_Kf,
         a_Ta, a_Te, a_Tf,
         a_Ae, a_Be ), ωs)
             for (a_vh, a_θh, a_P_g, a_Q_g,
                  a_ra, a_X_d, a_X_q,
                  a_X_d_dash, a_X_q_dash,
                  a_Ka, a_Ke, a_Kf,
                  a_Ta, a_Te, a_Tf,
                  a_Ae, a_Be ) in zip(
                      gens_vh, gens_θh,
                      gens_P_g, gens_Q_g,
                      ra, X_d, X_q,
                      X_d_dash, X_q_dash,
                      Ka ,Ke, Kf,
                      Ta, Te, Tf, Ae, Be ) ]

    """
    init_vec = []

    for (a_vh, a_θh, a_P_g, a_Q_g,
         a_ra, a_X_d, a_X_q,
         a_X_d_dash, a_X_q_dash,
         a_Ka, a_Ke, a_Kf,
         a_Ta, a_Te, a_Tf,
         a_Ae, a_Be ) in zip(
             gens_vh, gens_θh,
             gens_P_g, gens_Q_g,
             ra X_d, X_q,
             X_d_dash, X_q_dash,
             Ka ,Ke, Kf,
             Ta, Te, Tf, Ae, Be )

        a_gen_para =
            (a_ra, a_X_d, a_X_q,
             a_X_d_dash, a_X_q_dash)
        
        a_avr_para = (a_Ka, a_Ke, a_Kf,
                      a_Ta, a_Te, a_Tf,
                      a_Ae, a_Be )

        push!(init_vec,
              get_init_a_gen_full_model(
                  a_vh, a_θh,
                  a_P_g, a_Q_g,
                  a_gen_para,
                  a_avr_para,
                  ωs))
    end

    """
end


function get_state_init_SC_generic_model(
    init_SC_generic_model,
    ωs)

    
    (δ, ed_dash, eq_dash, E_fd, R_f, V_R) =
        init_SC_generic_model

    return [
        [ [δ_i, ωs, ed_dash_i, eq_dash_i,
           E_fd_i, R_f_i, V_R_i]
             for (δ_i, ed_dash_i, eq_dash_i,
           E_fd_i, R_f_i, V_R_i) in
                 zip(δ, ed_dash, eq_dash,
                     E_fd, R_f, V_R) ]...;]


end


#-----------------------------------------------------
# comment
#-----------------------------------------------------
