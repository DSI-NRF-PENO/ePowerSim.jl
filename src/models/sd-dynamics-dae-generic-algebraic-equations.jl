# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247


#####################################################
# ---------------------------------------------------
# DAE algebraic equations functions
# ---------------------------------------------------
#####################################################


#-----------------------------------------------------
# powerflow by current mismatch
#-----------------------------------------------------


function algebraic_generic_pf_ΔI_mismatch!(
    dx,
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_kwd_para  )

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes
     ) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes ))

    #----------------------------------------

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
    
    #----------------------------------------    

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
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------
    
    (Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         Ynet_wt_nodes_idx_wt_adjacent_nodes

    #----------------------------------------
    #----------------------------------------
    
    # du_vh  =  dx[dyn_pf_vh_Idxs]
    
    # du_θh  =  dx[dyn_pf_θh_Idxs]
    
    # du_gens_i_d  =  dx[dyn_pf_id_Idxs]
    
    # du_gens_i_q  =  dx[dyn_pf_iq_Idxs]
    

    vh = x[dyn_pf_vh_Idxs]
    
    θh = x[dyn_pf_θh_Idxs]
    
    gens_i_d = x[dyn_pf_id_Idxs]
    
    gens_i_q = x[dyn_pf_iq_Idxs]
    
    #----------------------------------------

    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]

    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------
    
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

    # # dx[dyn_pf_vh_Idxs]
    
    I_real_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
        ( (P_non_gens[ n2s_non_gens_idx[ nth_idx]] * cos(θh[ n2s_all_nodes_idx[ nth_idx]]) +
           Q_non_gens[ n2s_non_gens_idx[ nth_idx]] * sin(θh[ n2s_all_nodes_idx[ nth_idx]]))/
        vh[ n2s_all_nodes_idx[ nth_idx]] +
            sum( [vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * cos(θh[n2s_all_nodes_idx[idx]] + angle(ynj))
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) : nth_idx ∈ gens_nodes_with_loc_loads_idx ? 
         ( (gens_i_q[ n2s_gens_idx[ nth_idx]] * cos( gens_δ[ n2s_gens_idx[ nth_idx]]) +
         gens_i_d[ n2s_gens_idx[ nth_idx]] * sin( gens_δ[ n2s_gens_idx[ nth_idx]]) -
          (P_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx]] * cos(gens_θh[ n2s_gens_idx[ nth_idx]]) +
           Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx]] * sin(gens_θh[ n2s_gens_idx[ nth_idx]]))/
           gens_vh[ n2s_gens_idx[ nth_idx]])  -
           sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * cos(θh[ n2s_all_nodes_idx[idx]] + angle(ynj))
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx] ]) ] ) ) : ( (gens_i_q[ n2s_gens_idx[ nth_idx]] * cos( gens_δ[ n2s_gens_idx[ nth_idx]]) +
         gens_i_d[ n2s_gens_idx[ nth_idx]] * sin( gens_δ[ n2s_gens_idx[ nth_idx]]) ) -  sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * cos(θh[ n2s_all_nodes_idx[idx]] + angle(ynj))
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

    # dx[dyn_pf_θh_Idxs]
        
    I_imag_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
      ( (P_non_gens[ n2s_non_gens_idx[ nth_idx ]] * sin( θh[ n2s_all_nodes_idx[nth_idx]]) -
         Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] * cos( θh[ n2s_all_nodes_idx[nth_idx]])) /
        vh[ n2s_all_nodes_idx[nth_idx]] +
        sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * sin(θh[ n2s_all_nodes_idx[idx]] + angle(ynj))
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) : nth_idx ∈ gens_nodes_with_loc_loads_idx ?
        ( (gens_i_q[ n2s_gens_idx[nth_idx]] * sin( gens_δ[ n2s_gens_idx[nth_idx]]) -
         gens_i_d[ n2s_gens_idx[nth_idx]] * cos( gens_δ[ n2s_gens_idx[nth_idx]]) -
         (P_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx ]] *
         sin(gens_θh[ n2s_gens_idx[ nth_idx]]) -
       Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ nth_idx]] * cos(gens_θh[ n2s_gens_idx[ nth_idx]]))/
       gens_vh[ n2s_gens_idx[ nth_idx]]) -
         sum( [vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * sin( θh[ n2s_all_nodes_idx[idx]] + angle(ynj))
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) : ( (gens_i_q[ n2s_gens_idx[ nth_idx]] * sin( gens_δ[ n2s_gens_idx[ nth_idx]]) -
         gens_i_d[ n2s_gens_idx[ nth_idx]] * cos( gens_δ[ n2s_gens_idx[ nth_idx]]) ) -
         sum( [vh[ n2s_all_nodes_idx[idx]] * abs(ynj) * sin( θh[ n2s_all_nodes_idx[idx]] + angle(ynj))
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) 
                              for nth_idx in all_nodes_idx ]
    
    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """
               
    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin( gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[ n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  
    
    # ------------------------------------
        
    dx .=
        vcat(I_real_mismatch,
             I_imag_mismatch,
             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag)

    return nothing

end


function algebraic_generic_pf_ΔI_mismatch_sol(
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_sol_kwd_para )
    
    
    (dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:dyn_pf_flat_vh_flat_θh_id_iq_Idx,
              :flat_vh_flat_θh_id_iq_u0,
              :dyn_pf_solver,
              :algebraic_generic_model_kwd_para ))

    #----------------------------------------    
    
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol)  =
         NamedTupleTools.select(
             dyn_pf_solver,
             (:use_init_u0,
              :use_nlsolve,
              :pf_alg,
              :abstol,
              :reltol))

    #----------------------------------------        

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
   
    idx_range =
        first(dyn_pf_vh_Idxs):last(dyn_pf_iq_Idxs)
            
    #----------------------------------------    

    if use_init_u0 == true
        
        flat_vh_flat_θh_id_iq =
            flat_vh_flat_θh_id_iq_u0
    else

        flat_vh_flat_θh_id_iq = x
        
    end

    pf_fun_mismatch =
        algebraic_generic_pf_ΔI_mismatch!
        # ode_algebraic_generic_model_func!

    #----------------------------------------    
    
    
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_kwd_para),
                      flat_vh_flat_θh_id_iq,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    kwd_para =
                        algebraic_generic_model_kwd_para)),
            flat_vh_flat_θh_id_iq,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para),
                                    pf_alg )
        
    end
    

end

#-----------------------------------------------------
# powerflow by power mismatch
#-----------------------------------------------------



function pertubation_algebraic_generic_pf_ΔPQ_mismatch!(
    dx,
    x,
    (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para, Ynet);
    kwd_para =
        algebraic_generic_model_kwd_para,
    nodes_idx_with_adjacent_nodes_idx =
        nodes_idx_with_adjacent_nodes_idx)

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para
     ) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para ))

    #----------------------------------------
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
    
    #----------------------------------------    

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
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------
    
    # (Ynet,
    #  nodes_idx_with_adjacent_nodes_idx) =
    #      Ynet_wt_nodes_idx_wt_adjacent_nodes

    #----------------------------------------
    #----------------------------------------

     
    vh  = x[dyn_pf_vh_Idxs]
    
    θh  = x[dyn_pf_θh_Idxs]
    
    gens_i_d  = x[dyn_pf_id_Idxs]
    
    gens_i_q  = x[dyn_pf_iq_Idxs]

    #----------------------------------------

    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]

    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
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

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]
    
    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ]))
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

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) 
                              for nth_idx in all_nodes_idx ]

    
    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin( gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[ n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    
    # ------------------------------------
        
    dx .=
        vcat(P_mismatch,
             Q_mismatch,
             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag)
    
    # ------------------------------------

    return nothing

end


function algebraic_generic_pf_ΔPQ_mismatch!(
    dx,
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_kwd_para  )

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes
     ) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes ))

    #----------------------------------------
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
    
    #----------------------------------------    

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
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------
    
    (Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         Ynet_wt_nodes_idx_wt_adjacent_nodes

    #----------------------------------------
    #----------------------------------------

     
    vh  = x[dyn_pf_vh_Idxs]
    
    θh  = x[dyn_pf_θh_Idxs]
    
    gens_i_d  = x[dyn_pf_id_Idxs]
    
    gens_i_q  = x[dyn_pf_iq_Idxs]

    #----------------------------------------

    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]

    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
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

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]
    
    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ]))
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

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) 
                              for nth_idx in all_nodes_idx ]

    
    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin( gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[ n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    
    # ------------------------------------
        
    dx .=
        vcat(P_mismatch,
             Q_mismatch,
             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag)
    
    # ------------------------------------

    return nothing

end

# ode_pf_algebraic_generic_model_sol_ext

function algebraic_generic_pf_ΔPQ_mismatch_sol(
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_sol_kwd_para )

    
    (dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:dyn_pf_flat_vh_flat_θh_id_iq_Idx,
              :flat_vh_flat_θh_id_iq_u0,
              :dyn_pf_solver,
              :algebraic_generic_model_kwd_para ))

    #----------------------------------------    
    
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol)  =
         NamedTupleTools.select(
             dyn_pf_solver,
             (:use_init_u0,
              :use_nlsolve,
              :pf_alg,
              :abstol,
              :reltol))

    #----------------------------------------        

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
   
    idx_range =
        first(dyn_pf_vh_Idxs):last(dyn_pf_iq_Idxs)
            
    #----------------------------------------    

    if use_init_u0 == true
        
        flat_vh_flat_θh_id_iq =
            flat_vh_flat_θh_id_iq_u0
    else

        flat_vh_flat_θh_id_iq = x
        
    end

    pf_fun_mismatch =
        algebraic_generic_pf_ΔPQ_mismatch!
        # ode_algebraic_generic_model_func!

    #----------------------------------------    
        
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_kwd_para),
                      flat_vh_flat_θh_id_iq,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    kwd_para =
                        algebraic_generic_model_kwd_para)),
            flat_vh_flat_θh_id_iq,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para),
                                    pf_alg )

        
    end
    

end


#-----------------------------------------------------
#-----------------------------------------------------

#-----------------------------------------------------
# comment
#-----------------------------------------------------




