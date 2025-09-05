# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

####################################################

# using Pkg

# QuickDynamicTest =
#     joinpath(@__DIR__,"..","..")

# cd(QuickDynamicTest)

# Pkg.activate( QuickDynamicTest )

# #---------------------------------------------------
# #---------------------------------------------------

#---------------------------------------------------
#---------------------------------------------------


"""
https://discourse.julialang.org/t/bifurcationkit-jl-automatic-bifurcation-diagrams-in-julia/42192

https://discourse.julialang.org/t/producing-fold-bifurcation-diagram-using-birfurcationkit/96543

https://discourse.julialang.org/t/bifurcationkit-jl-automatic-bifurcation-diagrams-in-julia/42192/2

https://github.com/JuliaDynamics/Attractors.jl

https://discourse.julialang.org/t/too-many-bifurcation-points-marked/123504


"""

#---------------------------------------------------


# get_nodes_incident_edges

# get_nodes_incident_edges_by_orientations

# get_Cnb_by_orientations(edges_orientations)

# get_Cbn_by_orientations(edges_orientations)

# get_nodes_incident_edges_by_orientations( edges_orientations )


#---------------------------------------------------
#---------------------------------------------------
# Distributed slack bus
#---------------------------------------------------
#---------------------------------------------------


function get_red_model_distributed_slack_pf_ΔPQ_mismatch!(
    red_vh_θh_slack_value,
    ds_Pg_inj_Png_Qng;
    pf_model_kwd_para =
        pf_model_kwd_para )
    
    # -------------------------------------

    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,
     pf_vh_θh_idx_and_idx2Idx,

     gens_loss_participation) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,
              :pf_vh_θh_idx_and_idx2Idx,
              :gens_loss_participation))

    (;participating_gens,
     gens_loss_participation_factor) =
        NamedTupleTools.select(
            gens_loss_participation,
            (:participating_gens,
            :gens_loss_participation_factor))
    
    #-------------------------------
    
     (Pg_Idxs,
      Png_Idxs,
      Qng_Idxs) =
          NamedTupleTools.select(
              Pg_Png_Qng_Idx,
              (:dyn_P_gens_Idxs,
               :dyn_P_non_gens_Idxs,
               :dyn_Q_non_gens_Idxs))
    
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
      
      all_nodes_idx,

      n2s_slack_gens_idx,
      n2s_non_slack_gens_idx,
      n2s_gens_idx,
      n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx ) =
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
               
               :all_nodes_idx,

               :n2s_slack_gens_idx,
               :n2s_non_slack_gens_idx,
               :n2s_gens_idx,
               :n2s_non_gens_idx,
               :n2s_gens_with_loc_load_idxs,
               :n2s_all_nodes_idx))

    #-------------------------------

    (vh_Idxs,
     θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    
    #-------------------------------

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              Ynet_wt_nodes_idx_wt_adjacent_nodes,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))
        
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
    
    non_gens_vh  = red_vh_θh_slack_value[
        red_non_gens_vh_Idxs ]

    non_slack_gens_θh = red_vh_θh_slack_value[
        red_non_slack_gens_θh_Idxs ]

    non_gens_θh = red_vh_θh_slack_value[
        red_non_gens_θh_Idxs ]
    
    slack_value =
        red_vh_θh_slack_value[
            red_slack_value_Idxs]

    vh = [ idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ]] :
                    non_gens_vh[
                        n2s_non_gens_idx[ idx ]]
               for idx in all_nodes_idx ]

    θh = [ idx ∈ slack_gens_nodes_idx ?
            slack_gens_θh[ n2s_slack_gens_idx[ idx ]] : 
            idx ∈ non_slack_gens_nodes_idx ?
            non_slack_gens_θh[
                n2s_non_slack_gens_idx[ idx] ] :
                    non_gens_θh[
                        n2s_non_gens_idx[ idx ]]
        for idx in all_nodes_idx ]


    # gens_θh =
    #     θh[ gens_nodes_idx ]

    #-------------------------------

    Pg_net_inj =
        ds_Pg_inj_Png_Qng[
            Pg_Idxs]
    
    P_non_gens =
        ds_Pg_inj_Png_Qng[
            Png_Idxs]

    Q_non_gens =
        ds_Pg_inj_Png_Qng[
            Qng_Idxs]

    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor
    
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
    
    P_mismatch = [
        nth_idx ∈ non_gens_nodes_idx ?
            ( P_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[ nth_idx ]] *
            sum([ vh[ n2s_all_nodes_idx[ idx ]] *
            abs( ynj ) *
            cos(θh[ n2s_all_nodes_idx[ nth_idx ]] -
            θh[n2s_all_nodes_idx[ idx ]] - angle( ynj ))
                  for ( ynj, idx ) in
                    zip(Ynet[ n2s_all_nodes_idx[ nth_idx ]],
                        nodes_idx_with_adjacent_nodes_idx[
                            n2s_all_nodes_idx[nth_idx]])])) :
                                ( Pg_net_inj[ n2s_gens_idx[
                                    nth_idx] ] +
                                        gens_loss_participation[
                                            n2s_gens_idx[ nth_idx] ] -
     vh[n2s_all_nodes_idx[ nth_idx ]] *
     sum([ vh[ n2s_all_nodes_idx[ idx ]] *
     abs( ynj ) *
     cos( θh[ n2s_all_nodes_idx[ nth_idx ]] -
     θh[ n2s_all_nodes_idx[ idx ]] -
     angle( ynj ) )
           for (ynj, idx) in
               zip(
                  Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                   nodes_idx_with_adjacent_nodes_idx[
                       n2s_all_nodes_idx[nth_idx]])] ) )
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
    
    Q_mismatch = [
        (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
      vh[ n2s_all_nodes_idx[ nth_idx ]] *
      sum([ vh[ n2s_all_nodes_idx[ idx ]] *
      abs( ynj ) *
      sin(θh[ n2s_all_nodes_idx[ nth_idx ]] -
      θh[ n2s_all_nodes_idx[ idx ]] - angle(ynj))
                for (ynj, idx) in
                    zip(Ynet[ n2s_all_nodes_idx[ nth_idx ] ],
                        nodes_idx_with_adjacent_nodes_idx[
                            n2s_all_nodes_idx[nth_idx]])] ))
                   for nth_idx in all_nodes_idx
                       if nth_idx ∈ non_gens_nodes_idx  ]
    
    # ------------------------------------


    """
    Slack value
    
    """
    Δslack_value = slack_value .-
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )
    
    # ------------------------------------
    
    return vcat(P_mismatch,
                Q_mismatch,
                Δslack_value)

end


function get_model_distributed_slack_pf_ΔPQ_mismatch!(    
    vh_θh_slack_value,
    ds_Pg_inj_Qn_inj_Png_Qng;
    pf_model_kwd_para =
        pf_model_kwd_para )
    
    # -------------------------------------

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

    (Png_base, Qng_base,
     Pll_base, Qll_base) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             ( :P_non_gens, :Q_non_gens,
               :P_g_loc_load, :Q_g_loc_load))

    # base_P_loading = loc_load_exist == true ?
    #     sum([Png_base;Pll_base]) : sum(Png_base)
    
    # base_Q_loading = loc_load_exist == true ?
    #     sum([Qng_base;Qll_base]) : sum(Qng_base)


    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)
    
    #-------------------------------
    
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
    
    #-------------------------------

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
        
    #-------------------------------
    
     (Pg_Idxs,
      Qg_Idxs,
      Png_Idxs,
      Qng_Idxs) =
          NamedTupleTools.select(
              Pg_Qg_Png_Qng_Pll_Qll_Idx,
              (:dyn_P_gens_Idxs,
               :dyn_Q_gens_Idxs,
               :dyn_P_non_gens_Idxs,
               :dyn_Q_non_gens_Idxs))
    
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
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]
    
    #-------------------------------

    Pg_net_inj =
        ds_Pg_inj_Qn_inj_Png_Qng[
            Pg_Idxs]

    Qg_net_inj =
        ds_Pg_inj_Qn_inj_Png_Qng[
            Qg_Idxs]
    
    P_non_gens =
        ds_Pg_inj_Qn_inj_Png_Qng[
            Png_Idxs]

    Q_non_gens =
        ds_Pg_inj_Qn_inj_Png_Qng[
            Qng_Idxs]

    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    ΔP_loading = base_P_loading - sum(P_non_gens)
    
    ΔQ_loading = base_Q_loading - sum(Q_non_gens)

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
        ΔQ_loading * gens_reactive_power_particpation_factor
    
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
    
    P_mismatch = [
        nth_idx ∈ non_gens_nodes_idx ?
            ( P_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[ nth_idx ]] *
            sum([ vh[ n2s_all_nodes_idx[ idx ]] *
            abs( ynj ) *
            cos(θh[ n2s_all_nodes_idx[ nth_idx ]] -
            θh[n2s_all_nodes_idx[ idx ]] - angle( ynj ))
                  for ( ynj, idx ) in
                    zip(Ynet[ n2s_all_nodes_idx[ nth_idx ]],
                        nodes_idx_with_adjacent_nodes_idx[
                            n2s_all_nodes_idx[nth_idx]])])) :
                                ( Pg_net_inj[ n2s_gens_idx[
                                    nth_idx] ] +
                                        gens_loss_participation[
                                            n2s_gens_idx[ nth_idx] ] +
                                                gens_active_power_particpation[
                                                    n2s_gens_idx[ nth_idx] ] -
     vh[n2s_all_nodes_idx[ nth_idx ]] *
     sum([ vh[ n2s_all_nodes_idx[ idx ]] *
     abs( ynj ) *
     cos( θh[ n2s_all_nodes_idx[ nth_idx ]] -
     θh[ n2s_all_nodes_idx[ idx ]] -
     angle( ynj ) )
           for (ynj, idx) in
               zip(
                  Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                   nodes_idx_with_adjacent_nodes_idx[
                       n2s_all_nodes_idx[nth_idx]])] ) )
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

    Q_mismatch = [
        nth_idx ∈ non_gens_nodes_idx ?
        (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum([ vh[ n2s_all_nodes_idx[idx]] *
        abs(ynj) *
        sin(θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
              for (ynj, idx) in
                  zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                       nodes_idx_with_adjacent_nodes_idx[
                           n2s_all_nodes_idx[nth_idx]])])) : 
                               (Qg_net_inj[n2s_gens_idx[ nth_idx]] +
                               gens_reactive_power_particpation[
                                   n2s_gens_idx[ nth_idx] ]  -
                               vh[ n2s_all_nodes_idx[nth_idx]] *
                               sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                               sin(θh[ n2s_all_nodes_idx[nth_idx]] -
                               θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                                      for (ynj, idx) in
                                          zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                              nodes_idx_with_adjacent_nodes_idx[
                                                  n2s_all_nodes_idx[nth_idx]])] ) ) 
        for nth_idx in all_nodes_idx ]
    
    # ------------------------------------


    """
    Slack value
    
    """
    Δslack_value = slack_value .-
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )
    
    # ------------------------------------
    
    return vcat(P_mismatch,
                Q_mismatch,
                Δslack_value)

end

# ------------------------------------
# By Ynet flattend
# ------------------------------------


function get_model_distributed_slack_pf_ΔPQ_mismatch_Ynet_flattend(    
    vh_θh_slack_value,
    ds_pf_setpoint_model_para_wt_Ynet_flattend_comparray;
    pf_model_kwd_para =
        pf_model_kwd_para )
    
    # -------------------------------------
    
    @unpack P_gens, Q_gens, P_non_gens, Q_non_gens, P_g_loc_load, Q_g_loc_load, Ynet_real_imag_flattend = ds_pf_setpoint_model_para_wt_Ynet_flattend_comparray
    
    # -------------------------------------

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


    (;gens_reactive_power_particpation_factor,) =
        NamedTupleTools.select(
         reactive_power_disturbance_resolution_participation,
            (:gens_reactive_power_particpation_factor,))    

    #-------------------------------
    
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
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
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
        Ynet_real_imag_flattend[
            Ynet_real_Idxs ]

    Ynet_imag =
        Ynet_real_imag_flattend[
            Ynet_imag_Idxs ]
    
    Ynet = [Ynet_real[idx] + im * Ynet_imag[idx]  
            for idx in
                Ynet_rows_Idxs_in_flattend ]
    
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
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]
    
    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------


    (Png_base, Qng_base,
     Pll_base, Qll_base) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             ( :P_non_gens, :Q_non_gens,
               :P_g_loc_load, :Q_g_loc_load))

    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)
    
    #-------------------------------

    if loc_load_exist == true
        
        Δ_Pll_loading = P_g_loc_load - Pll_base 

        Δ_Qll_loading = Q_g_loc_load - Qll_base

    end
    
    
    ΔP_loading = base_P_loading - sum(P_non_gens)
    
    ΔQ_loading = base_Q_loading - sum(Q_non_gens)

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
        ΔQ_loading * gens_reactive_power_particpation_factor
    
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
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
         vh[ n2s_all_nodes_idx[nth_idx]] *
         sum([ vh[ n2s_all_nodes_idx[ idx]] *
         abs(ynj) *
         cos( θh[ n2s_all_nodes_idx[nth_idx]] -
         θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (P_gens[ n2s_gens_idx[nth_idx]] +
         Δ_Pll_loading[ n2s_gens_with_loc_load_idxs[nth_idx]] +
         gens_loss_participation[ n2s_gens_idx[ nth_idx] ] +
         gens_active_power_particpation[ n2s_gens_idx[ nth_idx]]  -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]]  -
          vh[ n2s_all_nodes_idx[nth_idx]] *
          sum([ vh[ n2s_all_nodes_idx[idx]] *
          abs(ynj) *
          cos(θh[ n2s_all_nodes_idx[ nth_idx]] -
          θh[ n2s_all_nodes_idx[ idx]] -
          angle(ynj) )
                for (ynj, idx) in
                    zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                         nodes_idx_with_adjacent_nodes_idx[
                             n2s_all_nodes_idx[ nth_idx]])])) :
        (P_gens[ n2s_gens_idx[nth_idx]] + gens_loss_participation[
                                            n2s_gens_idx[ nth_idx] ] +
                                                gens_active_power_particpation[
                                                    n2s_gens_idx[ nth_idx] ]  -
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
        cos(θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
              for (ynj, idx) in
                  zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                       nodes_idx_with_adjacent_nodes_idx[
                           n2s_all_nodes_idx[ nth_idx] ]) ]))
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
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
        sin(θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
        (Q_gens[ n2s_gens_idx[ nth_idx ]] +
        Δ_Qll_loading[ n2s_gens_with_loc_load_idxs[ nth_idx ]] +        
        gens_reactive_power_particpation[ n2s_gens_idx[ nth_idx ] ] -
        Q_g_loc_load[n2s_gens_with_loc_load_idxs[ nth_idx ]] -
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum( [ vh[ n2s_all_nodes_idx[idx]] *
        abs(ynj) *
        sin( θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
                   for (ynj, idx) in
                       zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                           nodes_idx_with_adjacent_nodes_idx[
                               n2s_all_nodes_idx[nth_idx]])])) :
       (Q_gens[ n2s_gens_idx[nth_idx]] +
        gens_reactive_power_particpation[ n2s_gens_idx[ nth_idx] ] -
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum( [ vh[ n2s_all_nodes_idx[idx]] *
        abs(ynj) *
        sin( θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
                  for (ynj, idx) in
                      zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                          nodes_idx_with_adjacent_nodes_idx[
                              n2s_all_nodes_idx[nth_idx]])] ) ) 
                   for nth_idx in all_nodes_idx ]
    
    # ------------------------------------


    """
    Slack value
    
    """
    Δslack_value = slack_value .-
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )
    
    # ------------------------------------
    
    return vcat(P_mismatch,
                Q_mismatch,
                Δslack_value)

end




function get_inj_model_distributed_slack_pf_ΔPQ_mismatch_Ynet_flattend(    
    vh_θh_slack_value,
    ds_pf_setpoint_inj_model_para_wt_Ynet_flattend_comparray;
    pf_model_kwd_para =
        pf_model_kwd_para )
    
    # -------------------------------------
    
    @unpack Pg_net_inj, Qg_net_inj, P_non_gens, Q_non_gens, Ynet_real_imag_flattend = ds_pf_setpoint_inj_model_para_wt_Ynet_flattend_comparray
    
    # -------------------------------------

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


    (;gens_reactive_power_particpation_factor,) =
        NamedTupleTools.select(
         reactive_power_disturbance_resolution_participation,
            (:gens_reactive_power_particpation_factor,))    

    #-------------------------------

    (Png_base, Qng_base,
     Pll_base, Qll_base) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             ( :P_non_gens, :Q_non_gens,
               :P_g_loc_load, :Q_g_loc_load))

    # base_P_loading = loc_load_exist == true ?
    #     sum([Png_base;Pll_base]) : sum(Png_base)
    
    # base_Q_loading = loc_load_exist == true ?
    #     sum([Qng_base;Qll_base]) : sum(Qng_base)


    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)
    
    #-------------------------------
    
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
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
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

    (vh_Idxs,
     θh_Idxs,
     slack_value_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_slack_value_Idxs))
    
    #-------------------------------
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]
    
    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    ΔP_loading = base_P_loading - sum(P_non_gens)
    
    ΔQ_loading = base_Q_loading - sum(Q_non_gens)

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
        ΔQ_loading * gens_reactive_power_particpation_factor
    
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
    
    P_mismatch = [
        nth_idx ∈ non_gens_nodes_idx ?
            ( P_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[ nth_idx ]] *
            sum([ vh[ n2s_all_nodes_idx[ idx ]] *
            abs( ynj ) *
            cos(θh[ n2s_all_nodes_idx[ nth_idx ]] -
            θh[n2s_all_nodes_idx[ idx ]] - angle( ynj ))
                  for ( ynj, idx ) in
                    zip(Ynet[ n2s_all_nodes_idx[ nth_idx ]],
                        nodes_idx_with_adjacent_nodes_idx[
                            n2s_all_nodes_idx[nth_idx]])])) :
                                ( Pg_net_inj[ n2s_gens_idx[
                                    nth_idx] ] +
                                        gens_loss_participation[
                                            n2s_gens_idx[ nth_idx] ] +
                                                gens_active_power_particpation[
                                                    n2s_gens_idx[ nth_idx] ] -
     vh[n2s_all_nodes_idx[ nth_idx ]] *
     sum([ vh[ n2s_all_nodes_idx[ idx ]] *
     abs( ynj ) *
     cos( θh[ n2s_all_nodes_idx[ nth_idx ]] -
     θh[ n2s_all_nodes_idx[ idx ]] -
     angle( ynj ) )
           for (ynj, idx) in
               zip(
                  Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                   nodes_idx_with_adjacent_nodes_idx[
                       n2s_all_nodes_idx[nth_idx]])] ) )
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

    Q_mismatch = [
        nth_idx ∈ non_gens_nodes_idx ?
        (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum([ vh[ n2s_all_nodes_idx[idx]] *
        abs(ynj) *
        sin(θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
              for (ynj, idx) in
                  zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                       nodes_idx_with_adjacent_nodes_idx[
                           n2s_all_nodes_idx[nth_idx]])])) : 
                               (Qg_net_inj[n2s_gens_idx[ nth_idx]] +
                               gens_reactive_power_particpation[
                                   n2s_gens_idx[ nth_idx] ]  -
                               vh[ n2s_all_nodes_idx[nth_idx]] *
                               sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                               sin(θh[ n2s_all_nodes_idx[nth_idx]] -
                               θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                                      for (ynj, idx) in
                                          zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                              nodes_idx_with_adjacent_nodes_idx[
                                                  n2s_all_nodes_idx[nth_idx]])] ) ) 
        for nth_idx in all_nodes_idx ]

    
    # ------------------------------------


    """
    Slack value
    
    """
    Δslack_value = slack_value .-
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )
    
    # ------------------------------------
    
    return vcat(P_mismatch,
                Q_mismatch,
                Δslack_value)

end



function get_model_distributed_slack_pf_ΔI_mismatch_Ynet_flattend(
    vh_θh_slack_value,
    pf_setpoint_model_para_wt_Ynet_flattend_comparray;
    pf_model_kwd_para =
        pf_model_kwd_para )
    
    # -------------------------------------
    
    @unpack P_gens, Q_gens, P_non_gens, Q_non_gens, P_g_loc_load, Q_g_loc_load, Ynet_real_imag_flattend  = pf_setpoint_model_para_wt_Ynet_flattend_comparray
    
    # -------------------------------------

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
    #-------------------------------

    (;slack_gens_vh,
     slack_gens_θh,

     # gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
             (:slack_gens_vh,
              :slack_gens_θh,

              # :gens_vh,
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

    (;gens_reactive_power_particpation_factor,) =
        NamedTupleTools.select(
         reactive_power_disturbance_resolution_participation,
            (:gens_reactive_power_particpation_factor,))    
    
    #-------------------------------
    
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
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
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

    (vh_Idxs,
     θh_Idxs,
     slack_value_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_slack_value_Idxs))
    
    #-------------------------------
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]

    gens_vh = vh[transformed_gens_nodes_idx]

    gens_θh = θh[transformed_gens_nodes_idx]

    gens_uh = gens_vh .* exp.(im * gens_θh)
    
    #-------------------------------

    (Png_base, Qng_base,
     Pll_base, Qll_base) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             ( :P_non_gens, :Q_non_gens,
               :P_g_loc_load, :Q_g_loc_load))

    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)
    
    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    ΔP_loading = base_P_loading - sum(P_non_gens)
    
    ΔQ_loading = base_Q_loading - sum(Q_non_gens)

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
       ΔQ_loading * gens_reactive_power_particpation_factor


    gens_power_disturbance_particpation =
        gens_active_power_particpation +
        im * gens_reactive_power_particpation
    
    # -----------------------------------
    # active current mismatch
    # -----------------------------------

    I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
        vh,
        θh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    # ------------------------------------------------   
    # update 
    # ------------------------------------------------

    gens_loss_current_contribution =
        gens_loss_participation ./ ( gens_vh .* cos.(
            gens_θh ) )

    gens_current_disturbance_particpation =
      conj.(gens_power_disturbance_particpation) ./ conj.(
            gens_uh )

    Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,

            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,

            gens_nodes_idx,
            gens_with_loc_load_idx ) +
                gens_loss_current_contribution  +
                 gens_current_disturbance_particpation
    
    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]
             ) .* conj.(
            Igen )

    # P_gens_update = [ real(
    #     gens_S[ n2s_gens_idx[idx]]) 
    #                   for idx in gens_nodes_idx]

    # Q_gens_update = [ imag(
    #     gens_S[ n2s_gens_idx[idx]]) 
    #                   for idx in gens_nodes_idx]

    # P_gens_update = real.(gens_S)

    # Q_gens_update = imag.(gens_S)

    P_gens_update = P_gens +
        gens_loss_participation +
        gens_active_power_particpation


    Q_gens_update = Q_gens +
        gens_reactive_power_particpation
    
    # -----------------------------------------------
    
    P_ΔP = P_non_gens .+ ΔP_loading
    
    Q_ΔP = Q_non_gens .+ ΔQ_loading
    
    current_mismatch =
        get_nodes_current_mismatch(
            vh,
            θh,
            P_gens_update,
            Q_gens_update,
            P_ΔP,
            Q_ΔP,
            # P_non_gens,
            # Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,    
            ( gens_nodes_idx,
              non_gens_nodes_idx,
              gens_with_loc_load_idx,
               all_nodes_idx,
              ),
            ( n2s_gens_idx,
             n2s_non_gens_idx,
              n2s_gens_with_loc_load_idxs,
              n2s_all_nodes_idx
              ),
            I_sum_ynj_vj;
            loc_load_exist =
                loc_load_exist)

    """
    Slack value
    
    """
    Δslack_value = slack_value .-
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )
    
    # ------------------------------------
    
    return vcat(real.(current_mismatch),
                imag.(current_mismatch),
                Δslack_value)
    
    
end


# ------------------------------------
# By edges parameters
# ------------------------------------


function get_model_distributed_slack_pf_ΔPQ_mismatch_edges_para(
    vh_θh_slack_value,
    ds_pf_setpoint_model_para_wt_edge_paras_comparray ;
    pf_model_kwd_para =
        pf_model_kwd_para,
    line_data_in_pu = true )
    
    # -------------------------------------

    @unpack P_gens, Q_gens, P_non_gens, Q_non_gens, P_g_loc_load, Q_g_loc_load, edges_r, edges_x, edges_b, edges_ratio, edges_angle, Gs, Bs = ds_pf_setpoint_model_para_wt_edge_paras_comparray 
    
    # -------------------------------------

    (;baseMVA,
     basekV,
     
     loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
     
     pf_vh_θh_idx_and_idx2Idx,

     shunt_idx,
     edges_fbus,
     edges_tbus,
     
     edges_type,
     edges_orientation,

     gens_loss_participation,
     active_power_disturbance_resolution_participation,
     reactive_power_disturbance_resolution_participation,

     sta_pf_PQ_para) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:baseMVA,
              :basekV,
              
              :loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,
              
              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
              
              :pf_vh_θh_idx_and_idx2Idx,

              :shunt_idx,
              :edges_fbus,
              :edges_tbus,
              
              :edges_type,
              :edges_orientation,

              :gens_loss_participation,
       :active_power_disturbance_resolution_participation,
       :reactive_power_disturbance_resolution_participation,

              :sta_pf_PQ_para ))

    #-------------------------------
    #-------------------------------

    (;slack_gens_vh,
     slack_gens_θh,

     # gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
             (:slack_gens_vh,
              :slack_gens_θh,

              # :gens_vh,
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

    (;gens_reactive_power_particpation_factor,) =
        NamedTupleTools.select(
         reactive_power_disturbance_resolution_participation,
            (:gens_reactive_power_particpation_factor,))    
    
    #-------------------------------
    
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
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
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

    if line_data_in_pu == true

        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx ) =
             get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
                 edges_fbus,
                 edges_tbus,
                 edges_r,
                 edges_x,
                 edges_b,
                 edges_ratio,
                 edges_angle,
                 Gs,
                 Bs;
                 edges_type,
                 all_nodes_idx,
                 n2s_all_nodes_idx,
                 
                 baseMVA = 1.0,
                 basekV = 1.0,
                 baseShunt = baseMVA)        
    else

        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx ) =
             get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
                 edges_fbus,
                 edges_tbus,
                 edges_r,
                 edges_x,
                 edges_b,
                 edges_ratio,
                 edges_angle,
                 Gs,
                 Bs;
                 edges_type,
                 all_nodes_idx,
                 n2s_all_nodes_idx,
                 baseMVA = baseMVA,
                 basekV = basekV,
                 baseShunt = baseMVA )        
    end
    
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
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]
    
    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------


    (Png_base, Qng_base,
     Pll_base, Qll_base) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             ( :P_non_gens, :Q_non_gens,
               :P_g_loc_load, :Q_g_loc_load))

    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)
    
    #-------------------------------

    if loc_load_exist == true
        
        Δ_Pll_loading = P_g_loc_load - Pll_base 

        Δ_Qll_loading = Q_g_loc_load - Qll_base

    end
    
    
    ΔP_loading = base_P_loading - sum(P_non_gens)
    
    ΔQ_loading = base_Q_loading - sum(Q_non_gens)

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
        ΔQ_loading * gens_reactive_power_particpation_factor
    
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
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
         vh[ n2s_all_nodes_idx[nth_idx]] *
         sum([ vh[ n2s_all_nodes_idx[ idx]] *
         abs(ynj) *
         cos( θh[ n2s_all_nodes_idx[nth_idx]] -
         θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (P_gens[ n2s_gens_idx[nth_idx]] +
         Δ_Pll_loading[ n2s_gens_with_loc_load_idxs[nth_idx]] +
         gens_loss_participation[ n2s_gens_idx[ nth_idx] ] +
         gens_active_power_particpation[ n2s_gens_idx[ nth_idx]]  -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]]  -
          vh[ n2s_all_nodes_idx[nth_idx]] *
          sum([ vh[ n2s_all_nodes_idx[idx]] *
          abs(ynj) *
          cos(θh[ n2s_all_nodes_idx[ nth_idx]] -
          θh[ n2s_all_nodes_idx[ idx]] -
          angle(ynj) )
                for (ynj, idx) in
                    zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                         nodes_idx_with_adjacent_nodes_idx[
                             n2s_all_nodes_idx[ nth_idx]])])) :
        (P_gens[ n2s_gens_idx[nth_idx]] + gens_loss_participation[
                                            n2s_gens_idx[ nth_idx] ] +
                                                gens_active_power_particpation[
                                                    n2s_gens_idx[ nth_idx] ]  -
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
        cos(θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
              for (ynj, idx) in
                  zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                       nodes_idx_with_adjacent_nodes_idx[
                           n2s_all_nodes_idx[ nth_idx] ]) ]))
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
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
        sin(θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
        (Q_gens[ n2s_gens_idx[ nth_idx ]] +
        Δ_Qll_loading[ n2s_gens_with_loc_load_idxs[ nth_idx ]] +        
        gens_reactive_power_particpation[ n2s_gens_idx[ nth_idx ] ] -
        Q_g_loc_load[n2s_gens_with_loc_load_idxs[ nth_idx ]] -
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum( [ vh[ n2s_all_nodes_idx[idx]] *
        abs(ynj) *
        sin( θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
                   for (ynj, idx) in
                       zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                           nodes_idx_with_adjacent_nodes_idx[
                               n2s_all_nodes_idx[nth_idx]])])) :
       (Q_gens[ n2s_gens_idx[nth_idx]] +
        gens_reactive_power_particpation[ n2s_gens_idx[ nth_idx] ] -
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum( [ vh[ n2s_all_nodes_idx[idx]] *
        abs(ynj) *
        sin( θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
                  for (ynj, idx) in
                      zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                          nodes_idx_with_adjacent_nodes_idx[
                              n2s_all_nodes_idx[nth_idx]])] ) ) 
                   for nth_idx in all_nodes_idx ]
    
    # ------------------------------------


    """
    Slack value
    
    """
    Δslack_value = slack_value .-
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )
    
    # ------------------------------------
    
    return vcat(P_mismatch,
                Q_mismatch,
                Δslack_value)

end




function get_model_distributed_slack_pf_ΔI_mismatch_by_edges_para(
    vh_θh_slack_value,
    ds_pf_setpoint_model_para_wt_edge_paras_comparray;
    pf_model_kwd_para =
        pf_model_kwd_para,
    line_data_in_pu = true )
    
    # -------------------------------------

    @unpack P_gens, Q_gens, P_non_gens, Q_non_gens, P_g_loc_load, Q_g_loc_load, edges_r, edges_x, edges_b, edges_ratio, edges_angle, Gs, Bs = ds_pf_setpoint_model_para_wt_edge_paras_comparray 
    
    # -------------------------------------

    (;baseMVA,
     basekV,
     
     loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
     
     pf_vh_θh_idx_and_idx2Idx,

     shunt_idx,
     edges_fbus,
     edges_tbus,
     
     edges_type,
     edges_orientation,

     gens_loss_participation,
     active_power_disturbance_resolution_participation,
     reactive_power_disturbance_resolution_participation,

     sta_pf_PQ_para) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:baseMVA,
              :basekV,
              
              :loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,
              
              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
              
              :pf_vh_θh_idx_and_idx2Idx,

              :shunt_idx,
              :edges_fbus,
              :edges_tbus,
              
              :edges_type,
              :edges_orientation,

              :gens_loss_participation,
       :active_power_disturbance_resolution_participation,
       :reactive_power_disturbance_resolution_participation,

              :sta_pf_PQ_para ))

    #-------------------------------
    #-------------------------------

    (;slack_gens_vh,
     slack_gens_θh,

     # gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
             (:slack_gens_vh,
              :slack_gens_θh,

              # :gens_vh,
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

    (;gens_reactive_power_particpation_factor,) =
        NamedTupleTools.select(
         reactive_power_disturbance_resolution_participation,
            (:gens_reactive_power_particpation_factor,))    
    
    #-------------------------------
    
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
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
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

    # (;Ynet,
    #  nodes_idx_with_adjacent_nodes_idx ) =
    #      NamedTupleTools.select(
    #          Ynet_wt_nodes_idx_wt_adjacent_nodes,
    #          (:Ynet,
    #           :nodes_idx_with_adjacent_nodes_idx))

    if line_data_in_pu == true

        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx ) =
             get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
                 edges_fbus,
                 edges_tbus,
                 edges_r,
                 edges_x,
                 edges_b,
                 edges_ratio,
                 edges_angle,
                 Gs,
                 Bs;
                 edges_type,
                 all_nodes_idx,
                 n2s_all_nodes_idx,
                 
                 baseMVA = 1.0,
                 basekV = 1.0,
                 baseShunt = baseMVA )        
    else

        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx ) =
             get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
                 edges_fbus,
                 edges_tbus,
                 edges_r,
                 edges_x,
                 edges_b,
                 edges_ratio,
                 edges_angle,
                 Gs,
                 Bs;
                 edges_type,
                 all_nodes_idx,
                 n2s_all_nodes_idx,
                 
                 baseMVA = baseMVA,
                 basekV = basekV,
                 baseShunt = baseMVA )        
    end
    
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
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]

    gens_vh = vh[transformed_gens_nodes_idx]

    gens_θh = θh[transformed_gens_nodes_idx]

    gens_uh = gens_vh .* exp.(im * gens_θh)
    
    #-------------------------------

    (Png_base, Qng_base,
     Pll_base, Qll_base) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             ( :P_non_gens, :Q_non_gens,
               :P_g_loc_load, :Q_g_loc_load))

    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)
    
    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    ΔP_loading = base_P_loading - sum(P_non_gens)
    
    ΔQ_loading = base_Q_loading - sum(Q_non_gens)

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
       ΔQ_loading * gens_reactive_power_particpation_factor


    gens_power_disturbance_particpation =
        gens_active_power_particpation +
        im * gens_reactive_power_particpation
    
    # -----------------------------------
    # active current mismatch
    # -----------------------------------

    I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
        vh,
        θh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    # ------------------------------------------------   
    # update 
    # ------------------------------------------------

    gens_loss_current_contribution =
        gens_loss_participation ./ ( gens_vh .* cos.(
            gens_θh ) )

    gens_current_disturbance_particpation =
      conj.(gens_power_disturbance_particpation) ./ conj.(
            gens_uh )

    Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,

            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,

            gens_nodes_idx,
            gens_with_loc_load_idx ) +
                gens_loss_current_contribution  +
                 gens_current_disturbance_particpation
    
    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]
             ) .* conj.(
            Igen )

    # P_gens_update = [ real(
    #     gens_S[ n2s_gens_idx[idx]]) 
    #                   for idx in gens_nodes_idx]

    # Q_gens_update = [ imag(
    #     gens_S[ n2s_gens_idx[idx]]) 
    #                   for idx in gens_nodes_idx]

    # P_gens_update = real.(gens_S)

    # Q_gens_update = imag.(gens_S)

    P_gens_update = P_gens +
        gens_loss_participation +
        gens_active_power_particpation


    Q_gens_update = Q_gens +
        gens_reactive_power_particpation
    
    # -----------------------------------------------
    
    P_ΔP = P_non_gens .+ ΔP_loading
    
    Q_ΔP = Q_non_gens .+ ΔQ_loading
    
    current_mismatch =
        get_nodes_current_mismatch(
            vh,
            θh,
            P_gens_update,
            Q_gens_update,
            P_ΔP,
            Q_ΔP,
            # P_non_gens,
            # Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,    
            ( gens_nodes_idx,
              non_gens_nodes_idx,
              gens_with_loc_load_idx,
               all_nodes_idx,
              ),
            ( n2s_gens_idx,
             n2s_non_gens_idx,
              n2s_gens_with_loc_load_idxs,
              n2s_all_nodes_idx
              ),
            I_sum_ynj_vj;
            loc_load_exist =
                loc_load_exist)

    """
    Slack value
    
    """
    Δslack_value = slack_value .-
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )
    
    # ------------------------------------
    
    return vcat(real.(current_mismatch),
                imag.(current_mismatch),
                Δslack_value)
    
    
end

# ------------------------------------
# ------------------------------------

function get_distributed_slack_pf_ΔPQ_mismatch_by_edges_para(    
    vh_θh_slack_value,
    ds_pf_PQ_setpoint_para_comparray;
    pf_model_kwd_para =
        pf_model_kwd_para,
    Ynet_wt_nodes_idx_wt_adjacent_nodes =
        Ynet_wt_nodes_idx_wt_adjacent_nodes )
    
    # -------------------------------------

    @unpack P_gens, Q_gens, P_non_gens, Q_non_gens, P_g_loc_load, Q_g_loc_load = ds_pf_PQ_setpoint_para_comparray
    
    # -------------------------------------

    (;baseMVA,
     basekV,

     loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

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
             (:baseMVA,
              :basekV,

              :loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

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


    (;gens_reactive_power_particpation_factor,) =
        NamedTupleTools.select(
         reactive_power_disturbance_resolution_participation,
            (:gens_reactive_power_particpation_factor,))    

    #-------------------------------
    
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
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
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

     (; Ynet,
      nodes_idx_with_adjacent_nodes_idx,) =
          NamedTupleTools.select(
              Ynet_wt_nodes_idx_wt_adjacent_nodes,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx, ))

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
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]
    
    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------


    (Png_base, Qng_base,
     Pll_base, Qll_base) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             ( :P_non_gens, :Q_non_gens,
               :P_g_loc_load, :Q_g_loc_load))

    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)
    
    #-------------------------------

    if loc_load_exist == true
        
        Δ_Pll_loading = P_g_loc_load - Pll_base 

        Δ_Qll_loading = Q_g_loc_load - Qll_base

    end
    
    
    ΔP_loading = base_P_loading - sum(P_non_gens)
    
    ΔQ_loading = base_Q_loading - sum(Q_non_gens)

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
        ΔQ_loading * gens_reactive_power_particpation_factor
    
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
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
         vh[ n2s_all_nodes_idx[nth_idx]] *
         sum([ vh[ n2s_all_nodes_idx[ idx]] *
         abs(ynj) *
         cos( θh[ n2s_all_nodes_idx[nth_idx]] -
         θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (P_gens[ n2s_gens_idx[nth_idx]] +
         Δ_Pll_loading[ n2s_gens_with_loc_load_idxs[nth_idx]] +
         gens_loss_participation[ n2s_gens_idx[ nth_idx] ] +
         gens_active_power_particpation[ n2s_gens_idx[ nth_idx]]  -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]]  -
          vh[ n2s_all_nodes_idx[nth_idx]] *
          sum([ vh[ n2s_all_nodes_idx[idx]] *
          abs(ynj) *
          cos(θh[ n2s_all_nodes_idx[ nth_idx]] -
          θh[ n2s_all_nodes_idx[ idx]] -
          angle(ynj) )
                for (ynj, idx) in
                    zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                         nodes_idx_with_adjacent_nodes_idx[
                             n2s_all_nodes_idx[ nth_idx]])])) :
        (P_gens[ n2s_gens_idx[nth_idx]] + gens_loss_participation[
                                            n2s_gens_idx[ nth_idx] ] +
                                                gens_active_power_particpation[
                                                    n2s_gens_idx[ nth_idx] ]  -
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
        cos(θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
              for (ynj, idx) in
                  zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                       nodes_idx_with_adjacent_nodes_idx[
                           n2s_all_nodes_idx[ nth_idx] ]) ]))
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
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
        sin(θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[
                                n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
        (Q_gens[ n2s_gens_idx[ nth_idx ]] +
        Δ_Qll_loading[ n2s_gens_with_loc_load_idxs[ nth_idx ]] +        
        gens_reactive_power_particpation[ n2s_gens_idx[ nth_idx ] ] -
        Q_g_loc_load[n2s_gens_with_loc_load_idxs[ nth_idx ]] -
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum( [ vh[ n2s_all_nodes_idx[idx]] *
        abs(ynj) *
        sin( θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
                   for (ynj, idx) in
                       zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                           nodes_idx_with_adjacent_nodes_idx[
                               n2s_all_nodes_idx[nth_idx]])])) :
       (Q_gens[ n2s_gens_idx[nth_idx]] +
        gens_reactive_power_particpation[ n2s_gens_idx[ nth_idx] ] -
        vh[ n2s_all_nodes_idx[nth_idx]] *
        sum( [ vh[ n2s_all_nodes_idx[idx]] *
        abs(ynj) *
        sin( θh[ n2s_all_nodes_idx[nth_idx]] -
        θh[ n2s_all_nodes_idx[idx]] -
        angle(ynj) )
                  for (ynj, idx) in
                      zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                          nodes_idx_with_adjacent_nodes_idx[
                              n2s_all_nodes_idx[nth_idx]])] ) ) 
                   for nth_idx in all_nodes_idx ]
    
    # ------------------------------------


    """
    Slack value
    
    """
    Δslack_value = slack_value .-
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )
    
    # ------------------------------------
    
    return vcat(P_mismatch,
                Q_mismatch,
                Δslack_value)

end



function get_distributed_slack_pf_ΔI_mismatch_by_edges_para(
    vh_θh_slack_value,
    ds_pf_PQ_setpoint_para_comparray;
    pf_model_kwd_para =
        pf_model_kwd_para,
    Ynet_wt_nodes_idx_wt_adjacent_nodes =
        Ynet_wt_nodes_idx_wt_adjacent_nodes )
    
    # -------------------------------------

    @unpack P_gens, Q_gens, P_non_gens, Q_non_gens, P_g_loc_load, Q_g_loc_load = ds_pf_PQ_setpoint_para_comparray
    
    # -------------------------------------

    (;baseMVA,
     basekV,
     
     loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
     
     pf_vh_θh_idx_and_idx2Idx,

     edges_type,
     edges_orientation,

     gens_loss_participation,
     active_power_disturbance_resolution_participation,
     reactive_power_disturbance_resolution_participation,

     sta_pf_PQ_para) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:baseMVA,
              :basekV,
              
              :loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,
              
              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
              
              :pf_vh_θh_idx_and_idx2Idx,

              :edges_type,
              :edges_orientation,

              :gens_loss_participation,
       :active_power_disturbance_resolution_participation,
       :reactive_power_disturbance_resolution_participation,

              :sta_pf_PQ_para ))

    #-------------------------------
    #-------------------------------

    (;slack_gens_vh,
     slack_gens_θh,

     # gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
             (:slack_gens_vh,
              :slack_gens_θh,

              # :gens_vh,
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

    (;gens_reactive_power_particpation_factor,) =
        NamedTupleTools.select(
         reactive_power_disturbance_resolution_participation,
            (:gens_reactive_power_particpation_factor,))    
    
    #-------------------------------
    
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
    
    #-------------------------------

   (;
    slack_gens_nodes_idx,
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

    (;Ynet,
     nodes_idx_with_adjacent_nodes_idx ) =
         NamedTupleTools.select(
             Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx))
    
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
    #-------------------------------
    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]

    gens_vh = vh[transformed_gens_nodes_idx]

    gens_θh = θh[transformed_gens_nodes_idx]

    gens_uh = gens_vh .* exp.(im * gens_θh)
    
    #-------------------------------

    (Png_base, Qng_base,
     Pll_base, Qll_base) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             ( :P_non_gens, :Q_non_gens,
               :P_g_loc_load, :Q_g_loc_load))

    base_P_loading =  sum(Png_base)
    
    base_Q_loading =  sum(Qng_base)
    
    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    ΔP_loading = base_P_loading - sum(P_non_gens)
    
    ΔQ_loading = base_Q_loading - sum(Q_non_gens)

    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
    gens_reactive_power_particpation =
       ΔQ_loading * gens_reactive_power_particpation_factor


    gens_power_disturbance_particpation =
        gens_active_power_particpation +
        im * gens_reactive_power_particpation
    
    # -----------------------------------
    # active current mismatch
    # -----------------------------------

    I_sum_ynj_vj = get_nodes_∑_ynj_x_vj(
        vh,
        θh,
        Ynet,
        nodes_idx_with_adjacent_nodes_idx,
        n2s_all_nodes_idx )

    # ------------------------------------------------   
    # update 
    # ------------------------------------------------

    gens_loss_current_contribution =
        gens_loss_participation ./ ( gens_vh .* cos.(
            gens_θh ) )

    gens_current_disturbance_particpation =
      conj.(gens_power_disturbance_particpation) ./ conj.(
            gens_uh )

    Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
        get_gens_loc_load_current(
            vh,
            θh,
            
            P_g_loc_load,
            Q_g_loc_load,

            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,

            gens_nodes_idx,
            gens_with_loc_load_idx ) +
                gens_loss_current_contribution  +
                 gens_current_disturbance_particpation
    
    gens_S  = vh[transformed_gens_nodes_idx] .*
        exp.(im * θh[transformed_gens_nodes_idx]
             ) .* conj.(
            Igen )

    # P_gens_update = [ real(
    #     gens_S[ n2s_gens_idx[idx]]) 
    #                   for idx in gens_nodes_idx]

    # Q_gens_update = [ imag(
    #     gens_S[ n2s_gens_idx[idx]]) 
    #                   for idx in gens_nodes_idx]

    # P_gens_update = real.(gens_S)

    # Q_gens_update = imag.(gens_S)

    P_gens_update = P_gens +
        gens_loss_participation +
        gens_active_power_particpation


    Q_gens_update = Q_gens +
        gens_reactive_power_particpation
    
    # -----------------------------------------------
    
    P_ΔP = P_non_gens .+ ΔP_loading
    
    Q_ΔP = Q_non_gens .+ ΔQ_loading
    
    current_mismatch =
        get_nodes_current_mismatch(
            vh,
            θh,
            P_gens_update,
            Q_gens_update,
            P_ΔP,
            Q_ΔP,
            # P_non_gens,
            # Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,    
            ( gens_nodes_idx,
              non_gens_nodes_idx,
              gens_with_loc_load_idx,
               all_nodes_idx,
              ),
            ( n2s_gens_idx,
             n2s_non_gens_idx,
              n2s_gens_with_loc_load_idxs,
              n2s_all_nodes_idx
              ),
            I_sum_ynj_vj;
            loc_load_exist =
                loc_load_exist)

    """
    Slack value
    
    """
    Δslack_value = slack_value .-
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )
    
    # ------------------------------------
    
    return vcat(real.(current_mismatch),
                imag.(current_mismatch),
                Δslack_value)
        
end

#----------------------------------------    
#----------------------------------------    

function get_distributed_slack_pf_mismatch_by_edges_para_sol(
    vh_θh_slack_value,
    ds_pf_PQ_setpoint_para_comparray,
    ds_pf_edge_para_comparray;
    pf_model_kwd_para =
        pf_model_kwd_para,
    line_data_in_pu = true,
    external_Ynet = nothing,    
    external_Ynet_bool = false )

    # -------------------------------------

    # @unpack edges_r, edges_x, edges_b, edges_ratio, edges_angle, Gs, Bs = ds_pf_edge_para_comparray
    
    # -------------------------------------
    
    (;baseMVA,
     basekV,
     
     edges_fbus,
     edges_tbus,
     edges_type,
     
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:baseMVA,
              :basekV,

              :edges_fbus,
              :edges_tbus,
              :edges_type,
              
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs))

    
   (; n2s_all_nodes_idx, ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            ( :n2s_all_nodes_idx, ))
    
    #-------------------------------

   (; all_nodes_idx, ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            ( :all_nodes_idx, ))

    #-------------------------------
    
    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------

    if  external_Ynet_bool == true
        
        Ynet_wt_nodes_idx_wt_adjacent_nodes =
            external_Ynet
    else


    # -------------------------------------
        
        Ynet_wt_nodes_idx_wt_adjacent_nodes =
            get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data(
                ;NamedTuple(ds_pf_edge_para_comparray)...,
                edges_fbus,
                edges_tbus,
                edges_type,
                all_nodes_idx,
                n2s_all_nodes_idx,
                baseMVA=baseMVA,
                basekV=1.0,
                line_data_in_pu =
                    true )        
    end
    

    #----------------------------------------    
    #----------------------------------------    

    pf_fun_mismatch =
        # get_distributed_slack_pf_ΔI_mismatch_by_edges_para
        get_distributed_slack_pf_ΔPQ_mismatch_by_edges_para

    #----------------------------------------    


    return NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( (x, p) ->
                pf_fun_mismatch(
                    x, p;
                    pf_model_kwd_para =
                        pf_model_kwd_para,
                    Ynet_wt_nodes_idx_wt_adjacent_nodes =
                        Ynet_wt_nodes_idx_wt_adjacent_nodes )),
            vh_θh_slack_value,
            ds_pf_PQ_setpoint_para_comparray ),
        pf_alg )

    

end



function get_distributed_slack_pf_ΔI_or_ΔPQ_mismatch_by_edges_para_sol(
    vh_θh_slack_value,
    ds_pf_PQ_setpoint_para_comparray,
    edges_paras;
    pf_model_kwd_para =
        pf_model_kwd_para,
    line_data_in_pu = true,
    external_Ynet = nothing,    
    external_Ynet_bool = false,
    pf_fun_mismatch =
        get_distributed_slack_pf_ΔI_mismatch_by_edges_para
       # get_distributed_slack_pf_ΔPQ_mismatch_by_edges_para
    )

    # -------------------------------------

    (;baseMVA,
     basekV,     
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:baseMVA,
              :basekV,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs))

    
   (; n2s_all_nodes_idx, ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            ( :n2s_all_nodes_idx, ))
    
    #-------------------------------

   (; all_nodes_idx, ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            ( :all_nodes_idx, ))

    #-------------------------------
    
    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    #-------------------------------

    if  external_Ynet_bool == true
        
        Ynet_wt_nodes_idx_wt_adjacent_nodes =
            external_Ynet
    else
        
        Ynet_wt_nodes_idx_wt_adjacent_nodes =
            get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data(
                ;edges_paras...,
                all_nodes_idx,
                n2s_all_nodes_idx,
                baseMVA=baseMVA,
                basekV=1.0,
                line_data_in_pu =
                    true )        
        
    end
    

    #----------------------------------------    
    #----------------------------------------    

    # pf_fun_mismatch =
    #     get_distributed_slack_pf_ΔI_mismatch_by_edges_para
    #    # get_distributed_slack_pf_ΔPQ_mismatch_by_edges_para

    #----------------------------------------    


    ds_pf_sol = NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( (x, p) ->
                pf_fun_mismatch(
                    x, p;
                    pf_model_kwd_para =
                        pf_model_kwd_para,
                    Ynet_wt_nodes_idx_wt_adjacent_nodes =
                        Ynet_wt_nodes_idx_wt_adjacent_nodes )),
            vh_θh_slack_value,
            ds_pf_PQ_setpoint_para_comparray ),
        pf_alg )

    return (; ds_pf_sol, Ynet_wt_nodes_idx_wt_adjacent_nodes )

end



# ------------------------------------
# ------------------------------------

function get_red_oop_pf_model_ΔPQ_mismatch_by_Ynet_flattend(
    red_vh_θh,
    pf_setpoint_model_para_wt_Ynet_flattend_comparray;
    pf_model_kwd_para =
        pf_model_kwd_para )

    #-------------------------------

    @unpack P_gens, Q_gens, P_non_gens, Q_non_gens, P_g_loc_load, Q_g_loc_load, Ynet_real_imag_flattend = pf_setpoint_model_para_wt_Ynet_flattend_comparray

    #-------------------------------

    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     pf_vh_θh_idx_and_idx2Idx,

     Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :pf_vh_θh_idx_and_idx2Idx,

              :Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend ))
    
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
              :non_slack_gens_vh))

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
    #-------------------------------
    

     (;dyn_pf_flat_vh_flat_θh_Idx,
      non_gens_vh_idx,
      non_slack_gens_θh_idx,
      non_gens_θh_idx,

      red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx,

      red_non_gens_vh_Idxs,
      red_non_slack_gens_θh_Idxs,
      red_non_gens_θh_Idxs,

      slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      
      gens_with_loc_load_idx,
      gens_nodes_with_loc_loads_idx,
      
      all_nodes_idx,

      n2s_slack_gens_idx,
      n2s_non_slack_gens_idx,
      n2s_gens_idx,
      n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx ) =
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

               :slack_gens_nodes_idx,
               :non_slack_gens_nodes_idx,
               :gens_nodes_idx,
               :non_gens_nodes_idx,
               
               :gens_with_loc_load_idx,
               :gens_nodes_with_loc_loads_idx,
               
               :all_nodes_idx,

               :n2s_slack_gens_idx,
               :n2s_non_slack_gens_idx,
               :n2s_gens_idx,
               :n2s_non_gens_idx,
               :n2s_gens_with_loc_load_idxs,
               :n2s_all_nodes_idx))
    
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
        
    non_gens_vh  = red_vh_θh[
        red_non_gens_vh_Idxs ]

    non_slack_gens_θh = red_vh_θh[
        red_non_slack_gens_θh_Idxs ]

    non_gens_θh = red_vh_θh[
        red_non_gens_θh_Idxs ]

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

    Q_gens_update = [ imag(
        gens_S[ n2s_gens_idx[idx]]) 
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

    power_mismatch =
        vh .* exp.(im * θh) .* conj.(
            current_mismatch)

    red_P_mismatch =
        real.(power_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)])

    red_Q_mismatch =
        imag.(power_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)])

    red_ΔPQ_x =
        vcat(red_P_mismatch,
             red_Q_mismatch)

    return red_ΔPQ_x


end



function get_red_oop_pf_model_ΔPQ_mismatch_by_edges_para(
    # red_ΔPQ_x,
    red_vh_θh,
    pf_setpoint_model_para_wt_edge_paras_comparray;
    pf_model_kwd_para =
        pf_model_kwd_para,
    line_data_in_pu = true )

    @unpack P_gens, Q_gens, P_non_gens, Q_non_gens, P_g_loc_load, Q_g_loc_load, edges_fbus, edges_tbus, edges_r, edges_x, edges_b, edges_ratio, edges_angle, Gs, Bs = pf_setpoint_model_para_wt_edge_paras_comparray

    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     pf_vh_θh_idx_and_idx2Idx,

     edges_type ) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :pf_vh_θh_idx_and_idx2Idx,

              :edges_type ))
    
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
              :non_slack_gens_vh))

    #-------------------------------

     (;dyn_pf_flat_vh_flat_θh_Idx,
      non_gens_vh_idx,
      non_slack_gens_θh_idx,
      non_gens_θh_idx,

      red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx,

      red_non_gens_vh_Idxs,
      red_non_slack_gens_θh_Idxs,
      red_non_gens_θh_Idxs,

      slack_gens_nodes_idx,
      non_slack_gens_nodes_idx,
      gens_nodes_idx,
      non_gens_nodes_idx,
      
      gens_with_loc_load_idx,
      gens_nodes_with_loc_loads_idx,
      
      all_nodes_idx,

      n2s_slack_gens_idx,
      n2s_non_slack_gens_idx,
      n2s_gens_idx,
      n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx ) =
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

               :slack_gens_nodes_idx,
               :non_slack_gens_nodes_idx,
               :gens_nodes_idx,
               :non_gens_nodes_idx,
               
               :gens_with_loc_load_idx,
               :gens_nodes_with_loc_loads_idx,
               
               :all_nodes_idx,

               :n2s_slack_gens_idx,
               :n2s_non_slack_gens_idx,
               :n2s_gens_idx,
               :n2s_non_gens_idx,
               :n2s_gens_with_loc_load_idxs,
               :n2s_all_nodes_idx))
    
    #-------------------------------    
    #-------------------------------

    if line_data_in_pu == true

        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx ) =
             get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
                 edges_fbus,
                 edges_tbus,
                 edges_r,
                 edges_x,
                 edges_b,
                 edges_ratio,
                 edges_angle,
                 Gs,
                 Bs;
                 edges_type,
                 all_nodes_idx,
                 n2s_all_nodes_idx,
                 baseMVA = 1.0,
                 basekV = 1.0,
                 baseShunt = baseMVA )
        
    else

        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx ) =
             get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
                 edges_fbus,
                 edges_tbus,
                 edges_r,
                 edges_x,
                 edges_b,
                 edges_ratio,
                 edges_angle,
                 Gs,
                 Bs;
                 edges_type,
                 all_nodes_idx,
                 n2s_all_nodes_idx,
                 baseMVA = baseMVA,
                 basekV = basekV,
                 baseShunt = baseMVA )
        
    end

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
        
    non_gens_vh  = red_vh_θh[
        red_non_gens_vh_Idxs ]

    non_slack_gens_θh = red_vh_θh[
        red_non_slack_gens_θh_Idxs ]

    non_gens_θh = red_vh_θh[
        red_non_gens_θh_Idxs ]

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

    Q_gens_update = [ imag(
        gens_S[ n2s_gens_idx[idx]]) 
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

    power_mismatch =
        vh .* exp.(im * θh) .* conj.(
            current_mismatch)

    red_P_mismatch =
        real.(power_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_slack_gens_nodes_idx)])

    red_Q_mismatch =
        imag.(power_mismatch[
            setdiff(transformed_all_nodes_idx,
                    transformed_gens_nodes_idx)])

    red_ΔPQ_x =
        vcat(red_P_mismatch,
             red_Q_mismatch)

    return red_ΔPQ_x


end



function get_red_iip_pf_model_ΔPQ_mismatch_generic(
    red_ΔPQ_x,
    red_vh_θh_x,
    pf_model_para;
    pf_model_kwd_para =
        pf_model_kwd_para )

    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     pf_vh_θh_idx_and_idx2Idx) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :pf_vh_θh_idx_and_idx2Idx))

    
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
              :non_slack_gens_vh))

     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              Ynet_wt_nodes_idx_wt_adjacent_nodes,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))
    
     (;dyn_P_gens_Idxs,
      dyn_Q_gens_Idxs,
      dyn_P_non_gens_Idxs,
      dyn_Q_non_gens_Idxs,
      dyn_P_gens_loc_load_Idxs,
      dyn_Q_gens_loc_load_Idxs ) =
          NamedTupleTools.select(
              Pg_Qg_Png_Qng_Pll_Qll_Idx,
              (:dyn_P_gens_Idxs,
               :dyn_Q_gens_Idxs,
               :dyn_P_non_gens_Idxs,
               :dyn_Q_non_gens_Idxs,
               :dyn_P_gens_loc_load_Idxs,
               :dyn_Q_gens_loc_load_Idxs))


     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx,

      red_non_gens_vh_Idxs,
      red_non_slack_gens_θh_Idxs,
      red_non_gens_θh_Idxs,

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
      n2s_all_nodes_idx ) =
          NamedTupleTools.select(
              pf_vh_θh_idx_and_idx2Idx,
              (:red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx,

               :red_non_gens_vh_Idxs,
               :red_non_slack_gens_θh_Idxs,
               :red_non_gens_θh_Idxs,

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
               :n2s_all_nodes_idx))
    
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
        pf_model_para[
            dyn_P_gens_Idxs ]

    Q_gens =
        pf_model_para[
            dyn_Q_gens_Idxs ]

    P_non_gens  =
        pf_model_para[
            dyn_P_non_gens_Idxs ]

    Q_non_gens = 
        pf_model_para[
            dyn_Q_non_gens_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_model_para[
                dyn_P_gens_loc_load_Idxs ]

        Q_g_loc_load =
            pf_model_para[
                dyn_Q_gens_loc_load_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #-------------------------------    
    #-------------------------------
    
    non_gens_vh  = red_vh_θh_x[
        red_non_gens_vh_Idxs ]

    non_slack_gens_θh = red_vh_θh_x[
        red_non_slack_gens_θh_Idxs ]

    non_gens_θh = red_vh_θh_x[
        red_non_gens_θh_Idxs ]

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

    Q_gens_update = [ imag(
        gens_S[ n2s_gens_idx[idx]]) 
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
              all_nodes_idx,
              # transformed_all_nodes_idx #
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

    power_mismatch =
        vh .* exp.(im * θh) .* conj.(
            current_mismatch)

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
# Comments
#---------------------------------------------------
