# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

####################################################

# using Pkg

# QuickDynamicTest =
#     joinpath(@__DIR__,"..","..")

# cd(QuickDynamicTest)

# Pkg.activate( QuickDynamicTest )

#---------------------------------------------------
#---------------------------------------------------

# using Plots
# using BifurcationKit

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
# continuation powerflow 
#---------------------------------------------------
#---------------------------------------------------

"""

Worked for powerflow but did not work for
continuation powerflow.

"""
function get_scaling_red_continuation_distributed_slack_pf_ΔPQ_mismatch!(
    
    red_vh_θh_slack_value,
    ds_scale_Pg_inj_Png_Qng;
    pf_model_kwd_para =
        pf_model_kwd_para )
    
    # -------------------------------------

    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     
     scale_Pg_Png_Qng_Idx,
     
     Pg_Png_Qng_Idx,
     pf_vh_θh_idx_and_idx2Idx,

     gens_loss_participation) =
         NamedTupleTools.select(
             pf_model_kwd_para,
             (:loc_load_exist,
              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              
              :scale_Pg_Png_Qng_Idx,
              
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
    
    (scale_Idxs,
     Pg_Idxs,
     Png_Idxs,
     Qng_Idxs) =
          NamedTupleTools.select(
              scale_Pg_Png_Qng_Idx,
              (:dyn_scale_Idxs,
               :dyn_P_gens_Idxs,
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
    
    non_gens_vh = red_vh_θh_slack_value[
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


    gens_θh =
        θh[ gens_nodes_idx ]

    #-------------------------------

    scale =
        ds_scale_Pg_inj_Png_Qng[
            scale_Idxs]

    Pg_net_inj =
        scale .* ds_scale_Pg_inj_Png_Qng[
            Pg_Idxs]
    
    P_non_gens =
        scale .* ds_scale_Pg_inj_Png_Qng[
            Png_Idxs]

    Q_non_gens =
        scale .* ds_scale_Pg_inj_Png_Qng[
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

"""

Worked for powerflow and continuation powerflow.

"""
function get_load_scaling_red_continuation_distributed_slack_pf_ΔPQ_mismatch!(    
    non_gens_vh_all_θh_slack_value,
    ds_scale_Pg_inj_Png_Qng;
    pf_model_kwd_para =
        pf_model_kwd_para,
    scaling_type = :P, # :Q, :P_and_Q
    scaling = :additive # :additive, :multiplicative
    )
    
    # -------------------------------------

    (;loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
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
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              
              :scale_Pg_Png_Qng_Idx,              
              :Pg_Png_Qng_Idx,
              
              :pf_vh_θh_idx_and_idx2Idx,
              
              :gens_loss_participation,

              :active_power_disturbance_resolution_participation,     
              :sta_pf_PQ_para ))
    
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
               
               :all_nodes_idx,

               :n2s_slack_gens_idx,
               :n2s_non_slack_gens_idx,
               :n2s_gens_idx,
               :n2s_non_gens_idx,
               :n2s_gens_with_loc_load_idxs,
               :n2s_all_nodes_idx,

              :non_gens_vh_gens_θh_non_gens_θh_slack_value_Idxs))

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

    
    #-------------------------------

    (Png_base,) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             (:P_non_gens, ))

    base_P_loading =  sum(Png_base)
    
    #-------------------------------
    
    (scale_Idxs,
     Pg_Idxs,
     Png_Idxs,
     Qng_Idxs) =
          NamedTupleTools.select(
              scale_Pg_Png_Qng_Idx,
              (:dyn_scale_Idxs,
               :dyn_P_gens_Idxs,
               :dyn_P_non_gens_Idxs,
               :dyn_Q_non_gens_Idxs))
    
    #-------------------------------

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
            slack_value_Idxs ]

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

    #-------------------------------

    scale =
        ds_scale_Pg_inj_Png_Qng[
            scale_Idxs]

    Pg_net_inj =
         ds_scale_Pg_inj_Png_Qng[
            Pg_Idxs]

    if scaling == :additive

        P_non_gens = scaling_type == :P ?
            (1 .+ scale) .* ds_scale_Pg_inj_Png_Qng[
                Png_Idxs] : ds_scale_Pg_inj_Png_Qng[
                Png_Idxs]

        Q_non_gens = scaling_type == :Q ?
            (1 .+ scale) .* ds_scale_Pg_inj_Png_Qng[
                Qng_Idxs] : ds_scale_Pg_inj_Png_Qng[
                    Qng_Idxs]

        if scaling_type == :P_and_Q

            P_non_gens = ( 1 .+ scale ) .*
                ds_scale_Pg_inj_Png_Qng[
                Png_Idxs ]

            Q_non_gens = ( 1 .+ scale ) .*
                ds_scale_Pg_inj_Png_Qng[
                Qng_Idxs ] 
        end
        
    elseif scaling == :multiplicative
    
        P_non_gens = scaling_type == :P ?
            scale .* ds_scale_Pg_inj_Png_Qng[
                Png_Idxs] : ds_scale_Pg_inj_Png_Qng[
                Png_Idxs]

        Q_non_gens = scaling_type == :Q ?
            scale .* ds_scale_Pg_inj_Png_Qng[
                Qng_Idxs] : ds_scale_Pg_inj_Png_Qng[
                    Qng_Idxs]

        if scaling_type == :P_and_Q

            P_non_gens = scale .*
                ds_scale_Pg_inj_Png_Qng[
                Png_Idxs ]

            Q_non_gens =  scale .*
                ds_scale_Pg_inj_Png_Qng[
                Qng_Idxs ] 
        end
        
    else
        
        P_non_gens = 
                ds_scale_Pg_inj_Png_Qng[
                    Png_Idxs ]
        Q_non_gens = 
                ds_scale_Pg_inj_Png_Qng[
                Qng_Idxs ] 
        
    end

    # -----------------------------------
            
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor
    
    # -----------------------------------
    
    ΔP_loading = sum(P_non_gens) - base_P_loading 
    
    gens_active_power_particpation =
        ΔP_loading * gens_active_power_particpation_factor
    
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
                                            n2s_gens_idx[ nth_idx] ] + gens_active_power_particpation[ n2s_gens_idx[ nth_idx] ] -
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


#---------------------------------------------------
#---------------------------------------------------


function get_load_scaling_continuation_distributed_slack_pf_ΔPQ_mismatch!(    
    vh_θh_slack_value,
    ds_scale_Pg_inj_Qn_inj_Png_Qng;
    pf_model_kwd_para =
        pf_model_kwd_para,
    scaling_type = :P,  # :Q, :P_and_Q
    scaling = :additive, # :additive, :multiplicative
    )
    
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

    (Png_base, Qng_base) =
         NamedTupleTools.select(
             sta_pf_PQ_para,
             ( :P_non_gens, :Q_non_gens))

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
    
    (scale_Idxs,
     Pg_Idxs,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs) =
          NamedTupleTools.select(
              scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              (:dyn_scale_Idxs,
               :dyn_P_gens_Idxs,
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
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]

    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    #-------------------------------
    #-------------------------------

    scale =
        ds_scale_Pg_inj_Qn_inj_Png_Qng[
            scale_Idxs]
    
    Pg_net_inj =
        ds_scale_Pg_inj_Qn_inj_Png_Qng[
            Pg_Idxs]

    Qg_net_inj =
        ds_scale_Pg_inj_Qn_inj_Png_Qng[
            Qg_Idxs]

    if scaling == :additive

        P_non_gens = scaling_type == :P ?
            (1 .+ scale) .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Png_Idxs] : ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Png_Idxs]

        Q_non_gens = scaling_type == :Q ?
            (1 .+ scale) .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Qng_Idxs] : ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Qng_Idxs]

        if scaling_type == :P_and_Q

            P_non_gens = (1 .+ scale) .*
                ds_scale_Pg_inj_Qn_inj_Png_Qng[Png_Idxs]

            Q_non_gens = (1 .+ scale) .*
                ds_scale_Pg_inj_Qn_inj_Png_Qng[Qng_Idxs] 
        end
        
    elseif scaling == :multiplicative
    
        P_non_gens = scaling_type == :P ?
            scale .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Png_Idxs] : ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Png_Idxs]

        Q_non_gens = scaling_type == :Q ?
            scale .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Qng_Idxs] : ds_scale_Pg_inj_Qn_inj_Png_Qng[
                    Qng_Idxs]

        if scaling_type == :P_and_Q

            P_non_gens = scale .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Png_Idxs]

            Q_non_gens = scale .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Qng_Idxs] 
        end
        
    else

            P_non_gens = ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Png_Idxs]

            Q_non_gens = ds_scale_Pg_inj_Qn_inj_Png_Qng[
                Qng_Idxs]         
    end
    
    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    ΔP_loading = sum(P_non_gens) - base_P_loading 
    
    ΔQ_loading = sum(Q_non_gens) - base_Q_loading 

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


function get_scaling_continuation_distributed_slack_pf_ΔPQ_mismatch!(    
    vh_θh_slack_value,
    ds_scale_Pg_inj_Qn_inj_Png_Qng;
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
    
    (scale_Idxs,
     Pg_Idxs,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs) =
          NamedTupleTools.select(
              scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              (:dyn_scale_Idxs,
               :dyn_P_gens_Idxs,
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
    
    slack_value =
        vh_θh_slack_value[
            slack_value_Idxs]

    
    vh  = vh_θh_slack_value[
        vh_Idxs ]
    
    θh  = vh_θh_slack_value[
        θh_Idxs ]
    
    #-------------------------------

    #-------------------------------

    scale =
        ds_scale_Pg_inj_Qn_inj_Png_Qng[
            scale_Idxs]
    
    Pg_net_inj =
        scale .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
            Pg_Idxs]

    Qg_net_inj =
        scale .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
            Qg_Idxs]
    
    P_non_gens =
        scale .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
            Png_Idxs]

    Q_non_gens =
        scale .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
            Qng_Idxs]

    #-------------------------------
        
    gens_loss_participation =
        slack_value .* gens_loss_participation_factor

    #-------------------------------
    
    ΔP_loading = sum(P_non_gens) - base_P_loading 
    
    ΔQ_loading = sum(Q_non_gens) - base_Q_loading 

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

#---------------------------------------------------
# Comments
#---------------------------------------------------

# function get_continuation_case3_distributed_slack_pf_ΔPQ_mismatch!(    
#     vh_θh_slack_value,
#     ds_scale_Pg_inj_Qn_inj_Png_Qng;
#     pf_model_kwd_para =
#         pf_model_kwd_para,
#     scaling_type =
#         :P # :Q, :P_and_Q
#     )
    
#     # -------------------------------------

#     (;loc_load_exist,
#      slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#      Ynet_wt_nodes_idx_wt_adjacent_nodes,

#      dyn_pf_fun_kwd_n2s_idxs,
#      dyn_pf_fun_kwd_net_idxs,
     
#      Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      Pg_Png_Qng_Idx,

#      scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
     
#      pf_vh_θh_idx_and_idx2Idx,

#      gens_loss_participation,
#      active_power_disturbance_resolution_participation,
#      reactive_power_disturbance_resolution_participation,

#      sta_pf_PQ_para) =
#          NamedTupleTools.select(
#              pf_model_kwd_para,
#              (:loc_load_exist,
#               :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#               :Ynet_wt_nodes_idx_wt_adjacent_nodes,

#               :dyn_pf_fun_kwd_n2s_idxs,
#               :dyn_pf_fun_kwd_net_idxs,
              
#               :Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               :Pg_Png_Qng_Idx,
              
#               :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
              
#               :pf_vh_θh_idx_and_idx2Idx,

#               :gens_loss_participation,
#               :active_power_disturbance_resolution_participation,
#               :reactive_power_disturbance_resolution_participation,

#               :sta_pf_PQ_para ))

#     #-------------------------------
#     #-------------------------------

#     (;slack_gens_vh,
#      slack_gens_θh,

#      gens_vh,
#      non_slack_gens_vh ) =
#          NamedTupleTools.select(
#              slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#              (:slack_gens_vh,
#               :slack_gens_θh,

#               :gens_vh,
#               :non_slack_gens_vh ))
        
#     #-------------------------------

#      (;Ynet,
#       nodes_idx_with_adjacent_nodes_idx) =
#           NamedTupleTools.select(
#               Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               (:Ynet,
#                :nodes_idx_with_adjacent_nodes_idx))
    
#     #-------------------------------

#     (;participating_gens,
#      gens_loss_participation_factor) =
#         NamedTupleTools.select(
#             gens_loss_participation,
#             (:participating_gens,
#             :gens_loss_participation_factor))


#     (;gens_active_power_particpation_factor,) =
#         NamedTupleTools.select(
#             active_power_disturbance_resolution_participation,
#             (:gens_active_power_particpation_factor,))


#     (;gens_reactive_power_particpation_factor,) =
#         NamedTupleTools.select(
#             reactive_power_disturbance_resolution_participation,
#             (:gens_reactive_power_particpation_factor,))    

#     #-------------------------------

#     (Png_base, Qng_base,
#      Pll_base, Qll_base) =
#          NamedTupleTools.select(
#              sta_pf_PQ_para,
#              ( :P_non_gens, :Q_non_gens,
#                :P_g_loc_load, :Q_g_loc_load))

#     # base_P_loading = loc_load_exist == true ?
#     #     sum([Png_base;Pll_base]) : sum(Png_base)
    
#     # base_Q_loading = loc_load_exist == true ?
#     #     sum([Qng_base;Qll_base]) : sum(Qng_base)


#     base_P_loading =  sum(Png_base)
    
#     base_Q_loading =  sum(Qng_base)
    
#     #-------------------------------
    
#    (;
#     n2s_slack_gens_idx,
#     n2s_non_slack_gens_idx,
#     n2s_gens_idx,
#     n2s_non_gens_idx,
#     n2s_gens_with_loc_load_idxs,
#     n2s_all_nodes_idx ) =
#         NamedTupleTools.select(
#             dyn_pf_fun_kwd_n2s_idxs,
#             (:n2s_slack_gens_idx,
#              :n2s_non_slack_gens_idx,
#              :n2s_gens_idx,
#              :n2s_non_gens_idx,
#              :n2s_gens_with_loc_load_idxs,
#              :n2s_all_nodes_idx))
    
#     #-------------------------------

#    (;
#     slack_gens_nodes_idx,
#     non_slack_gens_nodes_idx,
#     gens_nodes_idx,
#     non_gens_nodes_idx,
#     gens_nodes_with_loc_loads_idx,
#     all_nodes_idx ) =
#         NamedTupleTools.select(
#             dyn_pf_fun_kwd_net_idxs,
#             (:slack_gens_nodes_idx,
#              :non_slack_gens_nodes_idx,
#              :gens_nodes_idx,
#              :non_gens_nodes_idx,
#              :gens_nodes_with_loc_loads_idx,
#              :all_nodes_idx))
        
#     #-------------------------------
    
#     (scale_Idxs,
#      Pg_Idxs,
#      Qg_Idxs,
#      Png_Idxs,
#      Qng_Idxs) =
#           NamedTupleTools.select(
#               scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               (:dyn_scale_Idxs,
#                :dyn_P_gens_Idxs,
#                :dyn_Q_gens_Idxs,
#                :dyn_P_non_gens_Idxs,
#                :dyn_Q_non_gens_Idxs))
    
#     #-------------------------------

#     (vh_Idxs,
#      θh_Idxs,
#      slack_value_Idxs) =
#          NamedTupleTools.select(
#              dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
#              (:dyn_pf_vh_Idxs,
#               :dyn_pf_θh_Idxs,
#               :dyn_slack_value_Idxs))
    
#     #-------------------------------
    
#     slack_value =
#         vh_θh_slack_value[
#             slack_value_Idxs]

    
#     vh  = vh_θh_slack_value[
#         vh_Idxs ]
    
#     θh  = vh_θh_slack_value[
#         θh_Idxs ]
    
#     #-------------------------------
    
#     # cal_vh  = vh_θh_slack_value[
#     #     vh_Idxs ]
    
#     # non_gens_vh  = cal_vh[
#     #     non_gens_nodes_idx ]
    
#     # cal_θh  = vh_θh_slack_value[
#     #     θh_Idxs ]
    
#     # gens_θh = cal_θh[gens_nodes_idx]

#     # non_slack_gens_θh = cal_θh[
#     #     non_slack_gens_nodes_idx ]

#     # non_gens_θh = cal_θh[
#     #     non_gens_nodes_idx ]

#     # vh = [ idx ∈ gens_nodes_idx ?
#     #         gens_vh[
#     #             n2s_gens_idx[ idx ]] :
#     #                 non_gens_vh[
#     #                     n2s_non_gens_idx[ idx ]]
#     #            for idx in all_nodes_idx ]

#     # θh = [ idx ∈ slack_gens_nodes_idx ?
#     #         slack_gens_θh[ n2s_slack_gens_idx[ idx ]] : 
#     #         idx ∈ non_slack_gens_nodes_idx ?
#     #         non_slack_gens_θh[
#     #             n2s_non_slack_gens_idx[ idx] ] :
#     #                 non_gens_θh[
#     #                     n2s_non_gens_idx[ idx ]]
#     #     for idx in all_nodes_idx ]


#     #-------------------------------

#     scale =
#         ds_scale_Pg_inj_Qn_inj_Png_Qng[
#             scale_Idxs]
    
#     Pg_net_inj =
#          ds_scale_Pg_inj_Qn_inj_Png_Qng[
#             Pg_Idxs]

#     Qg_net_inj =
#         ds_scale_Pg_inj_Qn_inj_Png_Qng[
#             Qg_Idxs]
    
#     P_non_gens = scaling_type == :P ?
#         (1 .+ scale) .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
#             Png_Idxs] : ds_scale_Pg_inj_Qn_inj_Png_Qng[
#             Png_Idxs]

#     Q_non_gens = scaling_type == :Q ?
#         (1 .+ scale) .* ds_scale_Pg_inj_Qn_inj_Png_Qng[
#             Qng_Idxs] : ds_scale_Pg_inj_Qn_inj_Png_Qng[
#             Qng_Idxs]
    
#     if scaling_type == :P_and_Q
        
#         P_non_gens = (1 .+ scale) .*
#             ds_scale_Pg_inj_Qn_inj_Png_Qng[Png_Idxs]
        
#         Q_non_gens = (1 .+ scale) .*
#             ds_scale_Pg_inj_Qn_inj_Png_Qng[Qng_Idxs] 
#     end

#     #-------------------------------
        
#     gens_loss_participation =
#         slack_value .* gens_loss_participation_factor

#     #-------------------------------
    
#     ΔP_loading = sum(P_non_gens) - base_P_loading
    
#     ΔQ_loading = sum(Q_non_gens) - base_Q_loading 

#     gens_active_power_particpation =
#         ΔP_loading * gens_active_power_particpation_factor
    
#     gens_reactive_power_particpation =
#         ΔQ_loading * gens_reactive_power_particpation_factor
    
#     # -----------------------------------
#     # active power mismatch
#     # -----------------------------------
    
#     """
#     Gens nodes
#     Pg_i_inj =
#         vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

#     Load nodes
#     Pl_i =
#         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

#     """
    
#     P_mismatch = [
#         nth_idx ∈ non_gens_nodes_idx ?
#             ( P_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
#             vh[ n2s_all_nodes_idx[ nth_idx ]] *
#             sum([ vh[ n2s_all_nodes_idx[ idx ]] *
#             abs( ynj ) *
#             cos(θh[ n2s_all_nodes_idx[ nth_idx ]] -
#             θh[n2s_all_nodes_idx[ idx ]] - angle( ynj ))
#                   for ( ynj, idx ) in
#                     zip(Ynet[ n2s_all_nodes_idx[ nth_idx ]],
#                         nodes_idx_with_adjacent_nodes_idx[
#                             n2s_all_nodes_idx[nth_idx]])])) :
#                                 ( Pg_net_inj[ n2s_gens_idx[
#                                     nth_idx] ] +
#                                         gens_loss_participation[
#                                             n2s_gens_idx[ nth_idx] ] +
#                                                 gens_active_power_particpation[
#                                                     n2s_gens_idx[ nth_idx] ] -
#      vh[n2s_all_nodes_idx[ nth_idx ]] *
#      sum([ vh[ n2s_all_nodes_idx[ idx ]] *
#      abs( ynj ) *
#      cos( θh[ n2s_all_nodes_idx[ nth_idx ]] -
#      θh[ n2s_all_nodes_idx[ idx ]] -
#      angle( ynj ) )
#            for (ynj, idx) in
#                zip(
#                   Ynet[ n2s_all_nodes_idx[ nth_idx] ],
#                    nodes_idx_with_adjacent_nodes_idx[
#                        n2s_all_nodes_idx[nth_idx]])] ) )
#         for nth_idx in all_nodes_idx ]
    
#     # -----------------------------------
#     # reactive power mismatch
#     # -----------------------------------
    
#     """
#     Gens nodes
#     Qg_i_inj =
#         vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     Load nodes
#     Ql_i =
#         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     """

#     Q_mismatch = [
#         nth_idx ∈ non_gens_nodes_idx ?
#         (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
#         vh[ n2s_all_nodes_idx[nth_idx]] *
#         sum([ vh[ n2s_all_nodes_idx[idx]] *
#         abs(ynj) *
#         sin(θh[ n2s_all_nodes_idx[nth_idx]] -
#         θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
#               for (ynj, idx) in
#                   zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                        nodes_idx_with_adjacent_nodes_idx[
#                            n2s_all_nodes_idx[nth_idx]])])) : 
#                                (Qg_net_inj[n2s_gens_idx[ nth_idx]] +
#                                gens_reactive_power_particpation[
#                                    n2s_gens_idx[ nth_idx] ]  -
#                                vh[ n2s_all_nodes_idx[nth_idx]] *
#                                sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#                                sin(θh[ n2s_all_nodes_idx[nth_idx]] -
#                                θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
#                                       for (ynj, idx) in
#                                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
#                                               nodes_idx_with_adjacent_nodes_idx[
#                                                   n2s_all_nodes_idx[nth_idx]])] ) ) 
#         for nth_idx in all_nodes_idx ]
    
#     # ------------------------------------


#     """
#     Slack value
    
#     """
#     Δslack_value = slack_value .-
#         get_total_P_network_loss(
#             vh,
#             θh,
#             Ynet;
#             nodes_idx_with_adjacent_nodes_idx,
#             n2s_all_nodes_idx,
#             all_nodes_idx )
    
#     # ------------------------------------
    
#     return vcat(P_mismatch,
#                 Q_mismatch,
#                 Δslack_value)

# end

# #---------------------------------------------------
# # Continuation powerflow
# #---------------------------------------------------

# function get_continuation_pf_sta_ΔPQ_mismatch(    
#     red_vh_θh_x,
#     pf_PQ_param
#     ;pf_kw_para =
#         pf_kw_para )

#     #-------------------------------
    
#     red_ΔPQ_x = similar( red_vh_θh_x )

#     #-------------------------------
        
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
#       gens_with_loc_load_idx,
#       all_nodes_idx) =
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
#             pf_PQ_param[
#                 P_g_loc_load_sta_para_Idxs ]

#         Q_g_loc_load =
#             pf_PQ_param[
#                 Q_g_loc_load_sta_para_Idxs ]
        
#     else

#         P_g_loc_load = [0.0]
            
#         Q_g_loc_load = [0.0]
#     end
    
#     #-------------------------------

#     nodes_size = length( gens_nodes_idx ) +
#         length( non_gens_nodes_idx )
    
#     #-------------------------------

#     non_slack_gens_θh =
#         red_vh_θh_x[
#             red_non_slack_gens_θh_idx2Idx ]

#     non_gens_vh  =
#         red_vh_θh_x[ red_vh_Idxs ]

#     non_gens_θh =
#         red_vh_θh_x[
#             red_non_gens_θh_idx2Idx ]

#     vh = [
#         idx ∈ gens_nodes_idx ?
#             gens_vh[
#                 n2s_gens_idx[ idx ] ]  :
#                     non_gens_vh[
#                         n2s_non_gens_idx[ idx ] ]
#                for idx in 1:nodes_size ]

#     θh = [
#         idx ∈ slack_gens_nodes_idx ?
#             slack_gens_θh[ n2s_slack_gens_idx[idx] ] : 
#             idx ∈ non_slack_gens_nodes_idx ?
#             non_slack_gens_θh[
#                 n2s_non_slack_gens_idx[ idx] ] :
#                     non_gens_θh[
#                         n2s_non_gens_idx[ idx ] ]
#         for idx in 1:nodes_size ]
    
#     # -----------------------------------
#     # non slack gens active power mismatch
#     # -----------------------------------
    
#     if loc_load_exist == true

#         non_slack_gens_P_mismatch = [
#             nth_idx ∈ gens_with_loc_load_idx ? 
#             (P_gens[ n2s_gens_idx[nth_idx] ] -
#             P_g_loc_load[
#                 n2s_gens_with_loc_load_idxs[nth_idx ]] -
#                 vh[nth_idx] * sum( [
#                     vh[idx] *
#                         abs(ynj) *
#                         cos(θh[nth_idx] - θh[idx] -
#                         angle(ynj))
#                     for (ynj, idx) in
#                         zip( Ynet[ nth_idx ],
#                              nodes_idx_with_adjacent_nodes_idx[
#                                  nth_idx])])) :
#                                  (P_gens[ n2s_gens_idx[nth_idx]] -
#                 vh[nth_idx] * sum( [
#                     vh[idx] *
#                         abs(ynj) *
#                         cos(θh[nth_idx] - θh[idx] -
#                         angle(ynj))
#                     for (ynj, idx) in
#                         zip( Ynet[ nth_idx ],
#                              nodes_idx_with_adjacent_nodes_idx[nth_idx])]))
#             for nth_idx in
#                 non_slack_gens_nodes_idx ]
        
#     else
    
#     non_slack_gens_P_mismatch = [
#         P_gens[ n2s_gens_idx[nth_idx] ] -
#             vh[nth_idx] * sum( [
#                 vh[idx] *
#                     abs(ynj) *
#                     cos(θh[nth_idx] - θh[idx] -
#                     angle(ynj))
#                 for (ynj, idx) in
#                     zip( Ynet[ nth_idx ],
#                          nodes_idx_with_adjacent_nodes_idx[
#                              nth_idx ])])
#         for nth_idx in
#             non_slack_gens_nodes_idx ]

#     end

#     # -----------------------------------
#     # non_gens active power mismatch        
#     # -----------------------------------
    
        
#     non_gens_P_mismatch = [
#         P_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
#             vh[nth_idx] * sum( [
#                 vh[idx] *
#                     abs(ynj) *
#                     cos(θh[nth_idx] - θh[idx] - angle(ynj))
#                 for (ynj, idx) in
#                     zip( Ynet[ nth_idx ],
#                          nodes_idx_with_adjacent_nodes_idx[
#                              nth_idx ])])
#         for nth_idx in
#                 non_gens_nodes_idx ]

#     # ------------------------------------
#     # non_gens reactive power mismatch                
#     # ------------------------------------
        
#     non_gens_Q_mismatch = [
#         Q_non_gens[ n2s_non_gens_idx[ nth_idx ] ] +
#             vh[nth_idx] * sum( [
#                 vh[idx] *
#                     abs(ynj) *
#                     sin(θh[nth_idx] - θh[idx] -
#                     angle(ynj))
#                 for (ynj, idx) in
#                     zip( Ynet[ nth_idx ],
#                          nodes_idx_with_adjacent_nodes_idx[
#                              nth_idx ])])
#         for nth_idx in
#                 non_gens_nodes_idx ]

#     # ------------------------------------

#     red_ΔPQ_x .=
#         vcat(non_slack_gens_P_mismatch,
#              non_gens_P_mismatch,
#              non_gens_Q_mismatch)

#     return red_ΔPQ_x

# end


# #---------------------------------------------------
# #---------------------------------------------------

# function get_iip_continuation_pf_ΔPQ_mismatch(
#     ΔPQ,
#     vh_θh,
#     pf_model_para;
#     pf_model_kwd_para =
#         pf_model_kwd_para  )
    
#     #-------------------------------

#     (;loc_load_exist,
#      Ynet_wt_nodes_idx_wt_adjacent_nodes,
#      Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      pf_vh_θh_idx_and_idx2Idx) =
#           NamedTupleTools.select(
#              pf_model_kwd_para,
#              (:loc_load_exist,
#               :Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               :Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               :pf_vh_θh_idx_and_idx2Idx))

#     #-------------------------------
    
#      (Pg_Idx,
#      Qg_Idxs,
#      Png_Idxs,
#      Qng_Idxs,
#      Pgll_Idxs,
#      Qgll_Idxs ) =
#           NamedTupleTools.select(
#               Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               (:dyn_P_gens_Idxs,
#                :dyn_Q_gens_Idxs,
#                :dyn_P_non_gens_Idxs,
#                :dyn_Q_non_gens_Idxs,
#                :dyn_P_gens_loc_load_Idxs,
#                :dyn_Q_gens_loc_load_Idxs))
    
#     #-------------------------------

#      (dyn_pf_flat_vh_flat_θh_Idx,
#       non_gens_vh_idx,
#       non_slack_gens_θh_idx,
#       non_gens_θh_idx,
      
#       red_vh_Idxs,
#       red_non_slack_gens_θh_idx2Idx,
#       red_non_gens_θh_idx2Idx,

#       red_non_gens_vh_Idxs,
#       red_non_slack_gens_θh_Idxs,
#       red_non_gens_θh_Idxs,

#       slack_gens_nodes_idx,
#       non_slack_gens_nodes_idx,
#       gens_nodes_idx,
#       non_gens_nodes_idx,
      
#       gens_nodes_with_loc_loads_idx,
      
#       all_nodes_idx,

#       n2s_slack_gens_idx,
#       n2s_non_slack_gens_idx,
#       n2s_gens_idx,
#       n2s_non_gens_idx,
#       n2s_gens_with_loc_load_idxs,
#       n2s_all_nodes_idx ) =
#           NamedTupleTools.select(
#               pf_vh_θh_idx_and_idx2Idx,
#               (:dyn_pf_flat_vh_flat_θh_Idx,
#                :non_gens_vh_idx,
#                :non_slack_gens_θh_idx,
#                :non_gens_θh_idx,
               
#                :red_vh_Idxs,
#                :red_non_slack_gens_θh_idx2Idx,
#                :red_non_gens_θh_idx2Idx,

#                :red_non_gens_vh_Idxs,
#                :red_non_slack_gens_θh_Idxs,
#                :red_non_gens_θh_Idxs,

#                :slack_gens_nodes_idx,
#                :non_slack_gens_nodes_idx,
#                :gens_nodes_idx,
#                :non_gens_nodes_idx,
               
#                :gens_nodes_with_loc_loads_idx,
               
#                :all_nodes_idx,

#                :n2s_slack_gens_idx,
#                :n2s_non_slack_gens_idx,
#                :n2s_gens_idx,
#                :n2s_non_gens_idx,
#                :n2s_gens_with_loc_load_idxs,
#                :n2s_all_nodes_idx))

#     #-------------------------------

#     (;dyn_pf_vh_Idxs,
#      dyn_pf_θh_Idxs) =
#         NamedTupleTools.select(
#             dyn_pf_flat_vh_flat_θh_Idx,
#             (:dyn_pf_vh_Idxs,
#              :dyn_pf_θh_Idxs))
    
#     #-------------------------------

#      (;Ynet,
#       nodes_idx_with_adjacent_nodes_idx) =
#           NamedTupleTools.select(
#               Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               (:Ynet,
#                :nodes_idx_with_adjacent_nodes_idx))
    
#     #-------------------------------
    
#     P_gens =
#         pf_model_para[ Pg_Idx ]

#     Q_gens =
#         pf_model_para[ Qg_Idxs ]

#     P_non_gens  =
#         pf_model_para[ Png_Idxs ]

#     Q_non_gens = 
#         pf_model_para[ Qng_Idxs ]
    
#     if loc_load_exist == true

#         P_g_loc_load =
#             pf_model_para[ Pgll_Idxs ]

#         Q_g_loc_load =
#             pf_model_para[ Qgll_Idxs ]
        
#     else

#         P_g_loc_load = [0.0]
            
#         Q_g_loc_load = [0.0]
#     end
    
#     #-------------------------------   
    
#     vh =
#         vh_θh[ dyn_pf_vh_Idxs ]

#     θh =
#         vh_θh[ dyn_pf_θh_Idxs ]    

#     gens_vh =
#         vh[ gens_nodes_idx ]

#     gens_θh =
#         θh[ gens_nodes_idx ]

#     non_gens_nodes_vh =
#         vh[ non_gens_nodes_idx ]

#     non_gens_nodes_θh =
#         θh[ non_gens_nodes_idx ]

#     # -----------------------------------
#     # gens active power mismatch
#     # -----------------------------------
    
#     """
#     P_gens + Ploc_i =
#         vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

#     """
#     # # P_mismatch
#     # dx[dyn_pf_vh_Idxs]
    
#     P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
#          (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
#           vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
#           cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
#                    for (ynj, idx) in
#                        zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) : nth_idx ∈ gens_nodes_with_loc_loads_idx ?
#          ( P_gens[ n2s_gens_idx[ nth_idx]] -
#           P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
#           vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#                 cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
#                       for (ynj, idx) in
#                           zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
#                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
#         ( P_gens[ n2s_gens_idx[ nth_idx]] -
#          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#          cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
#                        for (ynj, idx) in
#                            zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
#                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ]))
#                         for nth_idx in all_nodes_idx ] 


#     # -----------------------------------
#     # reactive power mismatch
#     # -----------------------------------

    
#     """
#     gens:

#     id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
#             vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     loads:

#     Ql_i =
#          -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     """

#     # # Q_mismatch

#     Q_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
#         ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
#             vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#             sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
#                    for (ynj, idx) in
#                        zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
#                                 nth_idx ∈ gens_nodes_with_loc_loads_idx ?
#           ( Q_gens[ n2s_gens_idx[ nth_idx]] -
#            Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
#            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#            sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
#                    for (ynj, idx) in
#                            zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
#                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) : ( Q_gens[ n2s_gens_idx[ nth_idx]] - vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#             sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
#                        for (ynj, idx) in
#                            zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
#                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) 
#                               for nth_idx in all_nodes_idx ]


#     # ------------------------------------

#     ΔPQ .=
#         vcat(P_mismatch,
#              Q_mismatch)

#     return nothing

# end


# #---------------------------------------------------

# function get_a_oop_pf_model_sta_ΔPQ_mismatch_generic(
#     # red_ΔPQ_x,
#     red_vh_θh_x,
#     pf_model_para;
#     pf_model_kwd_para =
#         pf_model_kwd_para )

#     (;loc_load_exist,
#      slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#      Ynet_wt_nodes_idx_wt_adjacent_nodes,
#      Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      pf_vh_θh_idx_and_idx2Idx) =
#          NamedTupleTools.select(
#              pf_model_kwd_para,
#              (:loc_load_exist,
#               :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#               :Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               :Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               :pf_vh_θh_idx_and_idx2Idx))

    
#     #-------------------------------

#     (;slack_gens_vh,
#      slack_gens_θh,

#      gens_vh,
#      non_slack_gens_vh ) =
#          NamedTupleTools.select(
#              slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#              (:slack_gens_vh,
#               :slack_gens_θh,

#               :gens_vh,
#               :non_slack_gens_vh))

#      (;Ynet,
#       nodes_idx_with_adjacent_nodes_idx) =
#           NamedTupleTools.select(
#               Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               (:Ynet,
#                :nodes_idx_with_adjacent_nodes_idx))
    
#      (dyn_P_gens_Idxs,
#       dyn_Q_gens_Idxs,
#       dyn_P_non_gens_Idxs,
#       dyn_Q_non_gens_Idxs,
#       dyn_P_gens_loc_load_Idxs,
#       dyn_Q_gens_loc_load_Idxs ) =
#           NamedTupleTools.select(
#               Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               (:dyn_P_gens_Idxs,
#                :dyn_Q_gens_Idxs,
#                :dyn_P_non_gens_Idxs,
#                :dyn_Q_non_gens_Idxs,
#                :dyn_P_gens_loc_load_Idxs,
#                :dyn_Q_gens_loc_load_Idxs))


#      (;red_vh_Idxs,
#       red_non_slack_gens_θh_idx2Idx,
#       red_non_gens_θh_idx2Idx,

#       red_non_gens_vh_Idxs,
#       red_non_slack_gens_θh_Idxs,
#       red_non_gens_θh_Idxs,

#       slack_gens_nodes_idx,
#       non_slack_gens_nodes_idx,
#       gens_nodes_idx,
#       non_gens_nodes_idx,
#       gens_with_loc_load_idx,
#       all_nodes_idx,

#       n2s_slack_gens_idx,
#       n2s_non_slack_gens_idx,
#       n2s_gens_idx,
#       n2s_non_gens_idx,
#       n2s_gens_with_loc_load_idxs,
#       n2s_all_nodes_idx ) =
#           NamedTupleTools.select(
#               pf_vh_θh_idx_and_idx2Idx,
#               (:red_vh_Idxs,
#                :red_non_slack_gens_θh_idx2Idx,
#                :red_non_gens_θh_idx2Idx,

#                :red_non_gens_vh_Idxs,
#                :red_non_slack_gens_θh_Idxs,
#                :red_non_gens_θh_Idxs,

#                :slack_gens_nodes_idx,
#                :non_slack_gens_nodes_idx,
#                :gens_nodes_idx,
#                :non_gens_nodes_idx,
#                :gens_with_loc_load_idx,
#                :all_nodes_idx,

#                :n2s_slack_gens_idx,
#                :n2s_non_slack_gens_idx,
#                :n2s_gens_idx,
#                :n2s_non_gens_idx,
#                :n2s_gens_with_loc_load_idxs,
#                :n2s_all_nodes_idx))
    
#     #-------------------------------

#     transformed_gens_nodes_idx = [
#         n2s_all_nodes_idx[idx]
#         for idx in gens_nodes_idx ]


#     transformed_slack_gens_nodes_idx = [
#         n2s_all_nodes_idx[idx]
#         for idx in slack_gens_nodes_idx ]


#     transformed_all_nodes_idx = [
#         n2s_all_nodes_idx[idx]
#         for idx in all_nodes_idx ]
    
#     #-------------------------------
    
#     P_gens =
#         pf_model_para[
#             dyn_P_gens_Idxs ]

#     Q_gens =
#         pf_model_para[
#             dyn_Q_gens_Idxs ]

#     P_non_gens  =
#         pf_model_para[
#             dyn_P_non_gens_Idxs ]

#     Q_non_gens = 
#         pf_model_para[
#             dyn_Q_non_gens_Idxs ]
    
#     if loc_load_exist == true

#         P_g_loc_load =
#             pf_model_para[
#                 dyn_P_gens_loc_load_Idxs ]

#         Q_g_loc_load =
#             pf_model_para[
#                 dyn_Q_gens_loc_load_Idxs ]
        
#     else

#         P_g_loc_load = [0.0]
            
#         Q_g_loc_load = [0.0]
#     end
    
#     #-------------------------------

    
#     # nodes_size = length( gens_nodes_idx ) +
#     #     length( non_gens_nodes_idx )
    
#     #-------------------------------
    
#     non_gens_vh  = red_vh_θh_x[
#         red_non_gens_vh_Idxs ]

#     non_slack_gens_θh = red_vh_θh_x[
#         red_non_slack_gens_θh_Idxs ]

#     non_gens_θh = red_vh_θh_x[
#         red_non_gens_θh_Idxs ]

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
#         nodes_idx_with_adjacent_nodes_idx,
#         n2s_all_nodes_idx )

#     # ------------------------------------------------   
#     # update Q 
#     # ------------------------------------------------

#     Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
#         get_gens_loc_load_current(
#             vh,
#             θh,
            
#             P_g_loc_load,
#             Q_g_loc_load,

#             n2s_gens_with_loc_load_idxs,
#             n2s_all_nodes_idx, #

#             gens_nodes_idx,
#             gens_with_loc_load_idx )

    
#     gens_S  = vh[transformed_gens_nodes_idx] .*
#         exp.(im * θh[transformed_gens_nodes_idx]
#              ) .* conj.(
#             Igen )
    
#     # gens_S  = uh[gens_nodes_idx] .* conj.(Igen )

#     P_gens_update = [ idx ∈ slack_gens_nodes_idx ?
#         real( gens_S[ n2s_gens_idx[idx]]) : P_gens[
#             n2s_gens_idx[idx]]
#                       for idx in gens_nodes_idx]

#     Q_gens_update = [ imag(
#         gens_S[ n2s_gens_idx[idx]]) 
#                       for idx in gens_nodes_idx]

#     # ----------------------------------------------- 

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
#                all_nodes_idx, # transformed_all_nodes_idx #
#               ),
#             ( n2s_gens_idx,
#              n2s_non_gens_idx,
#               n2s_gens_with_loc_load_idxs,
#               n2s_all_nodes_idx #
#               ),
#             I_sum_ynj_vj;
#             loc_load_exist =
#                 loc_load_exist)
    
#     # -----------------------------------------------

#     power_mismatch =
#         vh .* exp.(im * θh) .* conj.(
#             current_mismatch)

#     red_P_mismatch =
#         real.(power_mismatch[
#             setdiff(transformed_all_nodes_idx,
#                     transformed_slack_gens_nodes_idx)])

#     red_Q_mismatch =
#         imag.(power_mismatch[
#             setdiff(transformed_all_nodes_idx,
#                     transformed_gens_nodes_idx)])

#     red_ΔPQ_x =
#         vcat(red_P_mismatch,
#              red_Q_mismatch)

#     return red_ΔPQ_x


# end


# function get_a_iip_pf_model_sta_ΔPQ_mismatch_generic(
#     red_ΔPQ_x,
#     red_vh_θh_x,
#     pf_model_para;
#     pf_model_kwd_para =
#         pf_model_kwd_para )

#     (;loc_load_exist,
#      slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#      Ynet_wt_nodes_idx_wt_adjacent_nodes,
#      Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      pf_vh_θh_idx_and_idx2Idx) =
#          NamedTupleTools.select(
#              pf_model_kwd_para,
#              (:loc_load_exist,
#               :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#               :Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               :Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               :pf_vh_θh_idx_and_idx2Idx))

    
#     #-------------------------------

#     (;slack_gens_vh,
#      slack_gens_θh,

#      gens_vh,
#      non_slack_gens_vh ) =
#          NamedTupleTools.select(
#              slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#              (:slack_gens_vh,
#               :slack_gens_θh,

#               :gens_vh,
#               :non_slack_gens_vh))

#      (;Ynet,
#       nodes_idx_with_adjacent_nodes_idx) =
#           NamedTupleTools.select(
#               Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               (:Ynet,
#                :nodes_idx_with_adjacent_nodes_idx))
    
#      (;dyn_P_gens_Idxs,
#       dyn_Q_gens_Idxs,
#       dyn_P_non_gens_Idxs,
#       dyn_Q_non_gens_Idxs,
#       dyn_P_gens_loc_load_Idxs,
#       dyn_Q_gens_loc_load_Idxs ) =
#           NamedTupleTools.select(
#               Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               (:dyn_P_gens_Idxs,
#                :dyn_Q_gens_Idxs,
#                :dyn_P_non_gens_Idxs,
#                :dyn_Q_non_gens_Idxs,
#                :dyn_P_gens_loc_load_Idxs,
#                :dyn_Q_gens_loc_load_Idxs))


#      (;red_vh_Idxs,
#       red_non_slack_gens_θh_idx2Idx,
#       red_non_gens_θh_idx2Idx,

#       red_non_gens_vh_Idxs,
#       red_non_slack_gens_θh_Idxs,
#       red_non_gens_θh_Idxs,

#       slack_gens_nodes_idx,
#       non_slack_gens_nodes_idx,
#       gens_nodes_idx,
#       non_gens_nodes_idx,
#       gens_with_loc_load_idx,
#       all_nodes_idx,

#       n2s_slack_gens_idx,
#       n2s_non_slack_gens_idx,
#       n2s_gens_idx,
#       n2s_non_gens_idx,
#       n2s_gens_with_loc_load_idxs,
#       n2s_all_nodes_idx ) =
#           NamedTupleTools.select(
#               pf_vh_θh_idx_and_idx2Idx,
#               (:red_vh_Idxs,
#                :red_non_slack_gens_θh_idx2Idx,
#                :red_non_gens_θh_idx2Idx,

#                :red_non_gens_vh_Idxs,
#                :red_non_slack_gens_θh_Idxs,
#                :red_non_gens_θh_Idxs,

#                :slack_gens_nodes_idx,
#                :non_slack_gens_nodes_idx,
#                :gens_nodes_idx,
#                :non_gens_nodes_idx,
#                :gens_with_loc_load_idx,
#                :all_nodes_idx,

#                :n2s_slack_gens_idx,
#                :n2s_non_slack_gens_idx,
#                :n2s_gens_idx,
#                :n2s_non_gens_idx,
#                :n2s_gens_with_loc_load_idxs,
#                :n2s_all_nodes_idx))
    
#     #-------------------------------

#     transformed_gens_nodes_idx = [
#         n2s_all_nodes_idx[idx]
#         for idx in gens_nodes_idx ]


#     transformed_slack_gens_nodes_idx = [
#         n2s_all_nodes_idx[idx]
#         for idx in slack_gens_nodes_idx ]


#     transformed_all_nodes_idx = [
#         n2s_all_nodes_idx[idx]
#         for idx in all_nodes_idx ]
    
#     #-------------------------------
    
#     P_gens =
#         pf_model_para[
#             dyn_P_gens_Idxs ]

#     Q_gens =
#         pf_model_para[
#             dyn_Q_gens_Idxs ]

#     P_non_gens  =
#         pf_model_para[
#             dyn_P_non_gens_Idxs ]

#     Q_non_gens = 
#         pf_model_para[
#             dyn_Q_non_gens_Idxs ]
    
#     if loc_load_exist == true

#         P_g_loc_load =
#             pf_model_para[
#                 dyn_P_gens_loc_load_Idxs ]

#         Q_g_loc_load =
#             pf_model_para[
#                 dyn_Q_gens_loc_load_Idxs ]
        
#     else

#         P_g_loc_load = [0.0]
            
#         Q_g_loc_load = [0.0]
#     end
    
#     #-------------------------------    
#     #-------------------------------
    
#     non_gens_vh  = red_vh_θh_x[
#         red_non_gens_vh_Idxs ]

#     non_slack_gens_θh = red_vh_θh_x[
#         red_non_slack_gens_θh_Idxs ]

#     non_gens_θh = red_vh_θh_x[
#         red_non_gens_θh_Idxs ]

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
#         nodes_idx_with_adjacent_nodes_idx,
#         n2s_all_nodes_idx )

#     # ------------------------------------------------   
#     # update Q 
#     # ------------------------------------------------

#     Igen = I_sum_ynj_vj[ transformed_gens_nodes_idx ] +
#         get_gens_loc_load_current(
#             vh,
#             θh,
            
#             P_g_loc_load,
#             Q_g_loc_load,

#             n2s_gens_with_loc_load_idxs,
#             n2s_all_nodes_idx, #

#             gens_nodes_idx,
#             gens_with_loc_load_idx )

    

#     gens_S  = vh[transformed_gens_nodes_idx] .*
#         exp.(im * θh[transformed_gens_nodes_idx]
#              ) .* conj.(
#             Igen )
    
#     # gens_S  = uh[gens_nodes_idx] .* conj.(Igen )

#     P_gens_update = [ idx ∈ slack_gens_nodes_idx ?
#         real( gens_S[ n2s_gens_idx[idx]]) : P_gens[
#             n2s_gens_idx[idx]]
#                       for idx in gens_nodes_idx]

#     Q_gens_update = [ imag(
#         gens_S[ n2s_gens_idx[idx]]) 
#                       for idx in gens_nodes_idx]

#     # ----------------------------------------------- 

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
#               all_nodes_idx,
#               # transformed_all_nodes_idx #
#               ),
#             ( n2s_gens_idx,
#              n2s_non_gens_idx,
#               n2s_gens_with_loc_load_idxs,
#               n2s_all_nodes_idx #
#               ),
#             I_sum_ynj_vj;
#             loc_load_exist =
#                 loc_load_exist)
    
#     # -----------------------------------------------

#     power_mismatch =
#         vh .* exp.(im * θh) .* conj.(
#             current_mismatch)

#     red_P_mismatch =
#         real.(power_mismatch[
#             setdiff(transformed_all_nodes_idx,
#                     transformed_slack_gens_nodes_idx)])

#     red_Q_mismatch =
#         imag.(power_mismatch[
#             setdiff(transformed_all_nodes_idx,
#                     transformed_gens_nodes_idx)])

#     red_ΔPQ_x .=
#         vcat(red_P_mismatch,
#              red_Q_mismatch)

#     return nothing


# end


# #---------------------------------------------------


# function get_continuation_model_pf_ΔPQ_mismatch!(
#     ΔPQ,
#     vh_θh,
#     Pg_inj_Qg_inj_Png_Qng;
#     pf_model_kwd_para =
#         pf_model_kwd_para  )
    
#     # -------------------------------------

#     (;loc_load_exist,
#      Ynet_wt_nodes_idx_wt_adjacent_nodes,
#      Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      pf_vh_θh_idx_and_idx2Idx) =
#          NamedTupleTools.select(
#              pf_model_kwd_para,
#              (:loc_load_exist,
#               :Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               :Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               :pf_vh_θh_idx_and_idx2Idx))

#     #-------------------------------
    
#      (Pg_Idxs,
#      Qg_Idxs,
#      Png_Idxs,
#      Qng_Idxs,
#      Pgll_Idxs,
#      Qgll_Idxs ) =
#           NamedTupleTools.select(
#               Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               (:dyn_P_gens_Idxs,
#                :dyn_Q_gens_Idxs,
#                :dyn_P_non_gens_Idxs,
#                :dyn_Q_non_gens_Idxs,
#                :dyn_P_gens_loc_load_Idxs,
#                :dyn_Q_gens_loc_load_Idxs))
    
#     #-------------------------------
    

#      (dyn_pf_flat_vh_flat_θh_Idx,
#       non_gens_vh_idx,
#       non_slack_gens_θh_idx,
#       non_gens_θh_idx,
      
#       red_vh_Idxs,
#       red_non_slack_gens_θh_idx2Idx,
#       red_non_gens_θh_idx2Idx,

#       red_non_gens_vh_Idxs,
#       red_non_slack_gens_θh_Idxs,
#       red_non_gens_θh_Idxs,

#       slack_gens_nodes_idx,
#       non_slack_gens_nodes_idx,
#       gens_nodes_idx,
#       non_gens_nodes_idx,
      
#       gens_nodes_with_loc_loads_idx,
      
#       all_nodes_idx,

#       n2s_slack_gens_idx,
#       n2s_non_slack_gens_idx,
#       n2s_gens_idx,
#       n2s_non_gens_idx,
#       n2s_gens_with_loc_load_idxs,
#       n2s_all_nodes_idx ) =
#           NamedTupleTools.select(
#               pf_vh_θh_idx_and_idx2Idx,
#               (:dyn_pf_flat_vh_flat_θh_Idx,
#                :non_gens_vh_idx,
#                :non_slack_gens_θh_idx,
#                :non_gens_θh_idx,
               
#                :red_vh_Idxs,
#                :red_non_slack_gens_θh_idx2Idx,
#                :red_non_gens_θh_idx2Idx,

#                :red_non_gens_vh_Idxs,
#                :red_non_slack_gens_θh_Idxs,
#                :red_non_gens_θh_Idxs,

#                :slack_gens_nodes_idx,
#                :non_slack_gens_nodes_idx,
#                :gens_nodes_idx,
#                :non_gens_nodes_idx,
               
#                :gens_nodes_with_loc_loads_idx,
               
#                :all_nodes_idx,

#                :n2s_slack_gens_idx,
#                :n2s_non_slack_gens_idx,
#                :n2s_gens_idx,
#                :n2s_non_gens_idx,
#                :n2s_gens_with_loc_load_idxs,
#                :n2s_all_nodes_idx))

#     #-------------------------------

#     (vh_Idxs,
#      θh_Idxs) =
#         NamedTupleTools.select(
#             dyn_pf_flat_vh_flat_θh_Idx,
#             (:dyn_pf_vh_Idxs,
#              :dyn_pf_θh_Idxs))
    
#     #-------------------------------


#      (;Ynet,
#       nodes_idx_with_adjacent_nodes_idx) =
#           NamedTupleTools.select(
#               Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               (:Ynet,
#                :nodes_idx_with_adjacent_nodes_idx))
    
    
#     #-------------------------------
#     #-------------------------------

#     Pg_net_inj =
#         Pg_inj_Qg_inj_Png_Qng[
#             Pg_Idxs]

#     Qg_net_inj =
#         Pg_inj_Qg_inj_Png_Qng[
#             Qg_Idxs]
    
#     P_non_gens =
#         Pg_inj_Qg_inj_Png_Qng[
#             Png_Idxs]

#     Q_non_gens =
#         Pg_inj_Qg_inj_Png_Qng[
#             Qng_Idxs]

#     #-------------------------------    
#     #-------------------------------
     
#     flat_vh =
#         vh_θh[
#             vh_Idxs ]

#     flat_θh =
#         vh_θh[
#             θh_Idxs ]

#     #-------------------------------
    
#     gens_vh =
#         flat_vh[
#             gens_nodes_idx ]

#     gens_θh =
#         flat_θh[
#             gens_nodes_idx ]

#     non_gens_nodes_vh =
#         flat_vh[
#             non_gens_nodes_idx ]

#     non_gens_nodes_θh =
#         flat_θh[
#             non_gens_nodes_idx ]
        
#     #-------------------------------
    
#     vh = [
#         idx ∈ gens_nodes_idx ?
#             gens_vh[
#                 n2s_gens_idx[ idx ] ]  :
#                     non_gens_nodes_vh[
#                         n2s_non_gens_idx[ idx ] ]
#                for idx in all_nodes_idx ]

#     θh = [
#         idx ∈ gens_nodes_idx ?
#             gens_θh[
#                 n2s_gens_idx[ idx] ] :
#                     non_gens_nodes_θh[
#                         n2s_non_gens_idx[ idx ] ]
#         for idx in all_nodes_idx ]

    
#     # -----------------------------------
#     # active power mismatch
#     # -----------------------------------
    
#     """
#     Gens nodes
#     Pg_i_inj =
#         vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

#     Load nodes
#     Pl_i =
#         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

#     """

    
#     P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
#     (P_non_gens[ n2s_non_gens_idx[ nth_idx ]] + vh[ n2s_all_nodes_idx[nth_idx]] *
#             sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#                     cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
#                 for (ynj, idx) in
#                     zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx] ])])) :     
#      (Pg_net_inj[ n2s_gens_idx[ nth_idx] ] -
#      vh[n2s_all_nodes_idx[nth_idx]] *
#      sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#                     cos( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
#               for (ynj, idx) in zip(
#                   Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                   nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) )
#     for nth_idx in all_nodes_idx ]

    
#     # -----------------------------------
#     # reactive power mismatch
#     # -----------------------------------
    
#     """
#     Gens nodes
#     Qg_i_inj =
#         vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     Load nodes
#     Ql_i =
#         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     """
    
#     Q_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
#       (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
#       vh[ n2s_all_nodes_idx[nth_idx]] *
#       sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#       sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
#                 for (ynj, idx) in
#                     zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :    
#      (Qg_net_inj[ n2s_gens_idx[ nth_idx] ] -
#       vh[ n2s_all_nodes_idx[ nth_idx]] *
#       sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
#       sin( θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj)) 
#                     for (ynj, idx) in
#                         zip(Ynet[ n2s_all_nodes_idx[ nth_idx]],
#                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])] ) )
#     for nth_idx in all_nodes_idx ]
    
#     # ------------------------------------

#     ΔPQ .=
#         vcat(P_mismatch,
#              Q_mismatch)

#     return nothing

# end


# function get_continuation_model_pf_ΔPQ_mismatch!(    
#     vh_θh,
#     Pg_inj_Qg_inj_Png_Qng;
#     pf_model_kwd_para =
#         pf_model_kwd_para  )
    
#     # -------------------------------------

#     (;loc_load_exist,
#      Ynet_wt_nodes_idx_wt_adjacent_nodes,
#      Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      pf_vh_θh_idx_and_idx2Idx) =
#          NamedTupleTools.select(
#              pf_model_kwd_para,
#              (:loc_load_exist,
#               :Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               :Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               :pf_vh_θh_idx_and_idx2Idx))

#     #-------------------------------
    
#      (Pg_Idxs,
#      Qg_Idxs,
#      Png_Idxs,
#      Qng_Idxs,
#      Pgll_Idxs,
#      Qgll_Idxs ) =
#           NamedTupleTools.select(
#               Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               (:dyn_P_gens_Idxs,
#                :dyn_Q_gens_Idxs,
#                :dyn_P_non_gens_Idxs,
#                :dyn_Q_non_gens_Idxs,
#                :dyn_P_gens_loc_load_Idxs,
#                :dyn_Q_gens_loc_load_Idxs))
    
#     #-------------------------------
    

#      (dyn_pf_flat_vh_flat_θh_Idx,
#       non_gens_vh_idx,
#       non_slack_gens_θh_idx,
#       non_gens_θh_idx,
      
#       red_vh_Idxs,
#       red_non_slack_gens_θh_idx2Idx,
#       red_non_gens_θh_idx2Idx,

#       red_non_gens_vh_Idxs,
#       red_non_slack_gens_θh_Idxs,
#       red_non_gens_θh_Idxs,

#       slack_gens_nodes_idx,
#       non_slack_gens_nodes_idx,
#       gens_nodes_idx,
#       non_gens_nodes_idx,
      
#       gens_nodes_with_loc_loads_idx,
      
#       all_nodes_idx,

#       n2s_slack_gens_idx,
#       n2s_non_slack_gens_idx,
#       n2s_gens_idx,
#       n2s_non_gens_idx,
#       n2s_gens_with_loc_load_idxs,
#       n2s_all_nodes_idx ) =
#           NamedTupleTools.select(
#               pf_vh_θh_idx_and_idx2Idx,
#               (:dyn_pf_flat_vh_flat_θh_Idx,
#                :non_gens_vh_idx,
#                :non_slack_gens_θh_idx,
#                :non_gens_θh_idx,
               
#                :red_vh_Idxs,
#                :red_non_slack_gens_θh_idx2Idx,
#                :red_non_gens_θh_idx2Idx,

#                :red_non_gens_vh_Idxs,
#                :red_non_slack_gens_θh_Idxs,
#                :red_non_gens_θh_Idxs,

#                :slack_gens_nodes_idx,
#                :non_slack_gens_nodes_idx,
#                :gens_nodes_idx,
#                :non_gens_nodes_idx,
               
#                :gens_nodes_with_loc_loads_idx,
               
#                :all_nodes_idx,

#                :n2s_slack_gens_idx,
#                :n2s_non_slack_gens_idx,
#                :n2s_gens_idx,
#                :n2s_non_gens_idx,
#                :n2s_gens_with_loc_load_idxs,
#                :n2s_all_nodes_idx))

#     #-------------------------------

#     (vh_Idxs,
#      θh_Idxs) =
#         NamedTupleTools.select(
#             dyn_pf_flat_vh_flat_θh_Idx,
#             (:dyn_pf_vh_Idxs,
#              :dyn_pf_θh_Idxs))
    
#     #-------------------------------


#      (;Ynet,
#       nodes_idx_with_adjacent_nodes_idx) =
#           NamedTupleTools.select(
#               Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               (:Ynet,
#                :nodes_idx_with_adjacent_nodes_idx))
    
    
#     #-------------------------------
#     #-------------------------------

#     Pg_net_inj =
#         Pg_inj_Qg_inj_Png_Qng[
#             Pg_Idxs]

#     Qg_net_inj =
#         Pg_inj_Qg_inj_Png_Qng[
#             Qg_Idxs]
    
#     P_non_gens =
#         Pg_inj_Qg_inj_Png_Qng[
#             Png_Idxs]

#     Q_non_gens =
#         Pg_inj_Qg_inj_Png_Qng[
#             Qng_Idxs]


#     #-------------------------------
    

#     # Pg_net_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
#     #     P_gens[ n2s_gens_idx[ idx]] -
#     #     P_g_loc_load[ n2s_gens_with_loc_load_idxs[idx]] :
#     #     P_gens[ n2s_gens_idx[ idx]]
#     #               for idx in gens_nodes_idx ]


#     # Qg_net_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
#     #     Q_gens[ n2s_gens_idx[ idx]] -
#     #     Q_g_loc_load[ n2s_gens_with_loc_load_idxs[idx]] :
#     #     Q_gens[ n2s_gens_idx[ idx] ]
#     #               for idx in gens_nodes_idx ]

    
#     #-------------------------------
     
#     flat_vh =
#         vh_θh[
#             vh_Idxs ]

#     flat_θh =
#         vh_θh[
#             θh_Idxs ]

#     #-------------------------------
    
#     gens_vh =
#         flat_vh[
#             gens_nodes_idx ]

#     gens_θh =
#         flat_θh[
#             gens_nodes_idx ]

#     non_gens_nodes_vh =
#         flat_vh[
#             non_gens_nodes_idx ]

#     non_gens_nodes_θh =
#         flat_θh[
#             non_gens_nodes_idx ]
        
#     #-------------------------------
    
#     vh = [
#         idx ∈ gens_nodes_idx ?
#             gens_vh[
#                 n2s_gens_idx[ idx ] ]  :
#                     non_gens_nodes_vh[
#                         n2s_non_gens_idx[ idx ] ]
#                for idx in all_nodes_idx ]

#     θh = [
#         idx ∈ gens_nodes_idx ?
#             gens_θh[
#                 n2s_gens_idx[ idx] ] :
#                     non_gens_nodes_θh[
#                         n2s_non_gens_idx[ idx ] ]
#         for idx in all_nodes_idx ]

    
#     # -----------------------------------
#     # active power mismatch
#     # -----------------------------------
    
#     """
#     Gens nodes
#     Pg_i_inj =
#         vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

#     Load nodes
#     Pl_i =
#         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

#     """

    
#     P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
#     (P_non_gens[ n2s_non_gens_idx[ nth_idx ]] + vh[ n2s_all_nodes_idx[nth_idx]] *
#             sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#                     cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
#                 for (ynj, idx) in
#                     zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx] ])])) :     
#      (Pg_net_inj[ n2s_gens_idx[ nth_idx] ] -
#      vh[n2s_all_nodes_idx[nth_idx]] *
#      sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#                     cos( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
#               for (ynj, idx) in zip(
#                   Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                   nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) )
#     for nth_idx in all_nodes_idx ]

    
#     # -----------------------------------
#     # reactive power mismatch
#     # -----------------------------------
    
#     """
#     Gens nodes
#     Qg_i_inj =
#         vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     Load nodes
#     Ql_i =
#         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     """
    
#     Q_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
#       (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
#       vh[ n2s_all_nodes_idx[nth_idx]] *
#       sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
#       sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj))
#                 for (ynj, idx) in
#                     zip(Ynet[ n2s_all_nodes_idx[nth_idx] ],
#                         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :    
#      (Qg_net_inj[ n2s_gens_idx[ nth_idx] ] -
#       vh[ n2s_all_nodes_idx[ nth_idx]] *
#       sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
#       sin( θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj)) 
#                     for (ynj, idx) in
#                         zip(Ynet[ n2s_all_nodes_idx[ nth_idx]],
#                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])] ) )
#     for nth_idx in all_nodes_idx ]
    
#     # ------------------------------------

#     return vcat(
#         P_mismatch, Q_mismatch)

# end




# function get_red_continuation_model_pf_ΔPQ_mismatch!(
    
#     red_vh_θh,
#     Pg_inj_Png_Qng;
#     pf_model_kwd_para =
#         pf_model_kwd_para  )
    
#     # -------------------------------------

#     (;loc_load_exist,
#      slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#      Ynet_wt_nodes_idx_wt_adjacent_nodes,
#      Pg_Qg_Png_Qng_Pll_Qll_Idx,
#      Pg_Png_Qng_Idx,
#      pf_vh_θh_idx_and_idx2Idx) =
#          NamedTupleTools.select(
#              pf_model_kwd_para,
#              (:loc_load_exist,
#               :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#               :Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               :Pg_Qg_Png_Qng_Pll_Qll_Idx,
#               :Pg_Png_Qng_Idx,
#               :pf_vh_θh_idx_and_idx2Idx ))

#     #-------------------------------
    
#      (Pg_Idxs,
#      Png_Idxs,
#      Qng_Idxs) =
#           NamedTupleTools.select(
#               Pg_Png_Qng_Idx,
#               (:dyn_P_gens_Idxs,
#                :dyn_P_non_gens_Idxs,
#                :dyn_Q_non_gens_Idxs))
    
#     #-------------------------------
    
#      (dyn_pf_flat_vh_flat_θh_Idx,
#       non_gens_vh_idx,
#       non_slack_gens_θh_idx,
#       non_gens_θh_idx,
      
#       red_vh_Idxs,
#       red_non_slack_gens_θh_idx2Idx,
#       red_non_gens_θh_idx2Idx,

#       red_non_gens_vh_Idxs,
#       red_non_slack_gens_θh_Idxs,
#       red_non_gens_θh_Idxs,

#       slack_gens_nodes_idx,
#       non_slack_gens_nodes_idx,
#       gens_nodes_idx,
#       non_gens_nodes_idx,
      
#       gens_nodes_with_loc_loads_idx,
      
#       all_nodes_idx,

#       n2s_slack_gens_idx,
#       n2s_non_slack_gens_idx,
#       n2s_gens_idx,
#       n2s_non_gens_idx,
#       n2s_gens_with_loc_load_idxs,
#       n2s_all_nodes_idx ) =
#           NamedTupleTools.select(
#               pf_vh_θh_idx_and_idx2Idx,
#               (:dyn_pf_flat_vh_flat_θh_Idx,
#                :non_gens_vh_idx,
#                :non_slack_gens_θh_idx,
#                :non_gens_θh_idx,
               
#                :red_vh_Idxs,
#                :red_non_slack_gens_θh_idx2Idx,
#                :red_non_gens_θh_idx2Idx,

#                :red_non_gens_vh_Idxs,
#                :red_non_slack_gens_θh_Idxs,
#                :red_non_gens_θh_Idxs,

#                :slack_gens_nodes_idx,
#                :non_slack_gens_nodes_idx,
#                :gens_nodes_idx,
#                :non_gens_nodes_idx,
               
#                :gens_nodes_with_loc_loads_idx,
               
#                :all_nodes_idx,

#                :n2s_slack_gens_idx,
#                :n2s_non_slack_gens_idx,
#                :n2s_gens_idx,
#                :n2s_non_gens_idx,
#                :n2s_gens_with_loc_load_idxs,
#                :n2s_all_nodes_idx))

#     #-------------------------------

#     (vh_Idxs,
#      θh_Idxs) =
#         NamedTupleTools.select(
#             dyn_pf_flat_vh_flat_θh_Idx,
#             (:dyn_pf_vh_Idxs,
#              :dyn_pf_θh_Idxs))
    
#     #-------------------------------

#      (;Ynet,
#       nodes_idx_with_adjacent_nodes_idx) =
#           NamedTupleTools.select(
#               Ynet_wt_nodes_idx_wt_adjacent_nodes,
#               (:Ynet,
#                :nodes_idx_with_adjacent_nodes_idx))
        
#     #-------------------------------
#     #-------------------------------

#     Pg_net_inj =
#         Pg_inj_Png_Qng[
#             Pg_Idxs]
    
#     P_non_gens =
#         Pg_inj_Png_Qng[
#             Png_Idxs]

#     Q_non_gens =
#         Pg_inj_Png_Qng[
#             Qng_Idxs]

#     #-------------------------------

#     (;slack_gens_vh,
#      slack_gens_θh,

#      gens_vh,
#      non_slack_gens_vh ) =
#          NamedTupleTools.select(
#              slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
#              (:slack_gens_vh,
#               :slack_gens_θh,

#               :gens_vh,
#               :non_slack_gens_vh ))
    
#     #-------------------------------
    
#     non_gens_vh  = red_vh_θh[
#         red_non_gens_vh_Idxs ]

#     non_slack_gens_θh = red_vh_θh[
#         red_non_slack_gens_θh_Idxs ]

#     non_gens_θh = red_vh_θh[
#         red_non_gens_θh_Idxs ]

#     vh = [ idx ∈ gens_nodes_idx ?
#             gens_vh[
#                 n2s_gens_idx[ idx ]] :
#                     non_gens_vh[
#                         n2s_non_gens_idx[ idx ]]
#                for idx in all_nodes_idx ]

#     θh = [ idx ∈ slack_gens_nodes_idx ?
#             slack_gens_θh[ n2s_slack_gens_idx[ idx ]] : 
#             idx ∈ non_slack_gens_nodes_idx ?
#             non_slack_gens_θh[
#                 n2s_non_slack_gens_idx[ idx] ] :
#                     non_gens_θh[
#                         n2s_non_gens_idx[ idx ]]
#         for idx in all_nodes_idx ]


#     gens_θh =
#         θh[ gens_nodes_idx ]
    
#     # -----------------------------------
#     # active power mismatch
#     # -----------------------------------
    
#     """
#     Gens nodes
#     Pg_i_inj =
#         vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

#     Load nodes
#     Pl_i =
#         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

#     """
    
#     P_mismatch = [
#         nth_idx ∈ non_gens_nodes_idx ?
#             ( P_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
#             vh[ n2s_all_nodes_idx[ nth_idx ]] *
#             sum([ vh[ n2s_all_nodes_idx[ idx ]] *
#             abs( ynj ) *
#             cos(θh[ n2s_all_nodes_idx[ nth_idx ]] -
#             θh[n2s_all_nodes_idx[ idx ]] - angle( ynj ))
#                   for ( ynj, idx ) in
#                     zip(Ynet[ n2s_all_nodes_idx[ nth_idx ]],
#                         nodes_idx_with_adjacent_nodes_idx[
#                             n2s_all_nodes_idx[nth_idx]])])) : ( Pg_net_inj[ n2s_gens_idx[ nth_idx] ] -
#      vh[n2s_all_nodes_idx[ nth_idx ]] *
#      sum([ vh[ n2s_all_nodes_idx[ idx ]] *
#      abs( ynj ) *
#      cos( θh[ n2s_all_nodes_idx[ nth_idx ]] -
#      θh[ n2s_all_nodes_idx[ idx ]] -
#      angle( ynj ) )
#            for (ynj, idx) in
#                zip(
#                   Ynet[ n2s_all_nodes_idx[ nth_idx] ],
#                    nodes_idx_with_adjacent_nodes_idx[
#                        n2s_all_nodes_idx[nth_idx]])] ) )
#         for nth_idx in all_nodes_idx ]
    
#     # -----------------------------------
#     # reactive power mismatch
#     # -----------------------------------
    
#     """
#     Gens nodes
#     Qg_i_inj =
#         vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     Load nodes
#     Ql_i =
#         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

#     """
    
#     Q_mismatch = [
#         (Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
#       vh[ n2s_all_nodes_idx[ nth_idx ]] *
#       sum([ vh[ n2s_all_nodes_idx[ idx ]] *
#       abs( ynj ) *
#       sin(θh[ n2s_all_nodes_idx[ nth_idx ]] -
#       θh[ n2s_all_nodes_idx[ idx ]] - angle(ynj))
#                 for (ynj, idx) in
#                     zip(Ynet[ n2s_all_nodes_idx[ nth_idx ] ],
#                         nodes_idx_with_adjacent_nodes_idx[
#                             n2s_all_nodes_idx[nth_idx]])] ))
#                    for nth_idx in all_nodes_idx
#                        if nth_idx ∈ non_gens_nodes_idx  ]
    
#     # ------------------------------------

#     return vcat(P_mismatch,
#              Q_mismatch)

# end

# #---------------------------------------------------

# function driver_continuation_pf_sta_ΔPQ_mismatch_sol_by_mpc(
#     case_file )


#     case_file = case_9_file
    
#     #-----------------------------------------------
    
                    
#     pf_alg  = NewtonRaphson()

#     #-----------------------------------------------

#     mpc_baseMVA =
#         get_matpower_scalar_as_iobuffer_by_case_file(
#             case_file;
#             type_key_string = "mpc_baseMVA" )[1]

#     mpc_bus =
#         get_matpower_mpc_type_iobuffer_by_case_file(
#             case_file; type_key= "mpc.bus" )


#     # mpc_gencost =
#     #     get_matpower_mpc_type_iobuffer_by_case_file(
#     #         case_file; type_key= "mpc.gencost" )


#     mpc_gen =
#         get_matpower_mpc_type_iobuffer_by_case_file(
#             case_file; type_key= "mpc.gen" )


#     mpc_branch =
#         get_matpower_mpc_type_iobuffer_by_case_file(
#             case_file; type_key= "mpc.branch" )

#     #-----------------------------------------------
#     #-----------------------------------------------

#     # (; pf_kw_para,
#     #  pf_PQ_param,
#     #  red_types_Idxs_etc,
#     #  net_para ) =
#     #     get_pf_sta_ΔPQ_mismatch_parameters_by_mpc(
#     #         mpc_gen, mpc_bus,
#     #         mpc_branch, mpc_baseMVA )


#     pf_sta_ΔPQ_mismatch_parameters =
#         get_pf_sta_ΔPQ_mismatch_parameters_by_mpc(
#             mpc_gen, mpc_bus,
#             mpc_branch, mpc_baseMVA )
    
#     pf_kw_para =
#         pf_sta_ΔPQ_mismatch_parameters.pf_kw_para

#     pf_PQ_param =
#         pf_sta_ΔPQ_mismatch_parameters.pf_PQ_param

#     red_types_Idxs_etc =
#         pf_sta_ΔPQ_mismatch_parameters.red_types_Idxs_etc

#     net_para =
#         pf_sta_ΔPQ_mismatch_parameters.net_para

#     #-----------------------------------------------
#     #-----------------------------------------------

#     red_vh_Idxs = red_types_Idxs_etc.red_vh_Idxs

#     red_θh_Idxs = red_types_Idxs_etc.red_θh_Idxs

#     #-----------------------------------------------
    
#     loc_load_exist =
#         pf_kw_para.loc_load_exist
    
#     gens_vh_slack_θh_para =
#         pf_kw_para.pf_kw_gens_vh_slack_θh_para
    
#     var_idxs =
#         pf_kw_para.pf_kw_var_idxs
    
#     nodes_types_idxs =
#         pf_kw_para.pf_kw_nodes_types_idxs
    
#     n2s_idxs =
#         pf_kw_para.pf_kw_n2s_idxs

#     #----------------------------------------

#     sta_red_vh_θh_0 =
#         [ ones(length(red_vh_Idxs));
#           zeros(length(red_θh_Idxs)) ]
    
#     #----------------------------------------

#     kwd_sta_sta_ΔPQ_sol_by_mpc =
#         (;
#          pf_alg,
#          pf_kw_para,
#          red_vh_Idxs,
#          red_θh_Idxs,
#          sta_red_vh_θh_0) 
    
#     red_pf_kw_para =
#         (;
#          gens_vh_slack_θh_para,
#          var_idxs,
#          nodes_types_idxs,
#          n2s_idxs )
    
#     #----------------------------------------
#     # Powerflow func and prob
#     #----------------------------------------

#     pf_sol =
#         get_pf_sta_ΔPQ_mismatch_sol_by_mpc(
#             pf_PQ_param;
#             kwd_para =
#                 kwd_sta_sta_ΔPQ_sol_by_mpc )

#     #----------------------------------------
#     # Results    
#     #----------------------------------------

#     tup_vh_θh =
#         get_vh_θh_from_red_pf_sol_u(
#             pf_sol;
#             kw_para =
#                 red_pf_kw_para )

    
#     #----------------------------------------
#     # Powerflow func and prob
#     #----------------------------------------

#     oop_pf_sol =
#         get_oop_pf_sta_ΔPQ_mismatch_sol_by_mpc(
#             pf_PQ_param;
#             kwd_para =
#                 kwd_sta_sta_ΔPQ_sol_by_mpc )

#     #----------------------------------------
#     # Results    
#     #----------------------------------------

#     tup_vh_θh =
#         get_vh_θh_from_red_pf_sol_u(
#             oop_pf_sol;
#             kw_para =
#                 red_pf_kw_para )

#     #----------------------------------------
#     # Continuation
#     #----------------------------------------

#     #Insert Explanation here
    
#     recordFromSolution(x, p; k...) =
#         (u1 = x[1], u2 = x[2],u3 = x[3],
#          u4 = x[4], u5 = x[5], u6 = x[6]) 

    
#     Fbp = get_continuation_pf_sta_ΔPQ_mismatch
    
#     # bifurcation problem
    
#     prob = BifurcationProblem(
#         (x,p)-> Fbp(
#             x, p;
#             pf_kw_para =
#                 pf_kw_para),
#         oop_pf_sol.u,
#         pf_PQ_param,
#         # specify the continuation parameter or its index
#         8,
#         # record_from_solution = (x, p; k...) -> x[1]

#         record_from_solution = recordFromSolution
#     )    

#     # options for continuation
    
#     opts_br = ContinuationPar(

#         # parameter interval
#         p_max = 350.0, p_min = 80.0,
        
#         # detect bifurcations with bisection method
#         # we increase the precision of the bisection
#         n_inversion = 4)

#     # automatic bifurcation diagram computation
    
#     diagram = bifurcationdiagram(prob, PALC(),

#     # very important parameter. This specifies the
#     # maximum amount of recursion

#     # when computing the bifurcation diagram. It means
#     # we allow computing branches of branches

#     # at most in the present case.
#             2, opts_br, )
    
#     plot(diagram)
    
# end



 
# function continuation_powerflow()

#     Fbp = get_continuation_pf_sta_ΔPQ_mismatch
    
#     # bifurcation problem
    
#     prob = BifurcationProblem(
#         (x,p)-> Fbp( x, p; pf_kw_para = pf_kw_para),
#         red_vh_θh_x,
#         pf_PQ_param,

#         # specify the continuation parameter
#         # or its index
        
#         1,
#         record_from_solution =
#             (x, p; k...) -> x[1])    

#     # options for continuation
    
#     opts_br = ContinuationPar(

#         # parameter interval
#         p_max = 2.5, p_min = 0.9,
        
#         # detect bifurcations with bisection method
#         # we increase the precision of the bisection
#         n_inversion = 4)

#     # automatic bifurcation diagram computation
    
#     diagram = bifurcationdiagram(prob, PALC(),

#     # very important parameter. This specifies the
#     # maximum amount of recursion

#     # when computing the bifurcation diagram. It means
#     # we allow computing branches of branches

#     # at most in the present case.
#             2, opts_br, )
    

#     plot(diagram)
    
    
# end

