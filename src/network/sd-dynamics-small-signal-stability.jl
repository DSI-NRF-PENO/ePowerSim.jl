# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.4:  pg 135 - 136, 6.121 - 6.123

####################################################


function get_sensitivity_and_stability_para_Idxs_by_per_gen_func(
    netd;
    abstol = 1e-14,
    reltol = 1e-14,    
    pf_alg  = NewtonRaphson(),
    ode_alg_2 = Rodas4(),
    # algr_name = "rodas4" 
    ode_alg = ImplicitMidpoint(),
    algr_name = "ImplicitMidpoint",
    dt = 0.01,
    sim_timespan    = (0.0, 10.0),
    only_gen         = false,
    no_control_device = false,
    all_wt_size_and_dims = false,
    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial
    Idxs_type        = :Idxs_im, 
    diagnostics_data = false,
    dynamics_by_gens = false,
    dynamics_by_per_gen = true,
    streamlined =  true)
    
    #-----------------------------------------------
    #-----------------------------------------------


    paras_and_Idxs_set =
        get_dynamics_paras_and_Idxs_set(
            netd; Idxs_type =
                Idxs_type,
            no_control_device =
                no_control_device,
            only_gen  =
                only_gen,
            all_wt_size_and_dims =
                all_wt_size_and_dims,
            pf_alg  =
                pf_alg,
            streamlined =
                streamlined )
    

    (; loc_load_exist,
     gens_Ax_update_parameters,
     im_vars_Idx_in_state,
     nodes_ur_ui_Idx_in_state,
     nodes_u_Idx,
     ur_ui_Idx_in_state,

     gens_nodes_u_Idx_in_ranges,
     non_gens_nodes_u_Idx_in_ranges,

     flat_ur_flat_ui_Idx,

     state,

     vec_Ax_views,

     each_gens_im_vars_Idx_in_state,
     nodes_state_Idx,
     nodes_δ_ed_dash_eq_dash_Idxs,

     gens_nodes_collection,

     f_dyn_ode_pf_para_Idx,

     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,

     ur_ui_idx_in_Idx,
     vh_θh_idx_in_Idx,
     
     nodes_u_Idx_in_ranges,

     intra_flat_ur_ui_Idx,
     intra_flat_vh_θh_Idx, 
     intra_dyn_pf_flat_para_Idx,

     im_mass_matrix,
     im_ode_mass_matrix,
     im_model_pf_mass_matrix,
     labels_and_symbols,
     gens_nodes_im_vars_labels, 
     net_bus_volts_labels,
     im_net_states_and_var_labels,
     im_sym,
     im_ode_sym,

     #
     post_pf_idxs,
     intra_pf_kwd_para,
     intra_dyn_pf_mismatch_kwd_para,

     dynamics_by_gens_kwd_para,
     dynamics_by_per_gen_kwd_para,
     states_and_matrices,
     δ_ed_dash_eq_dash_Idxs_in_flattend,


     flat_δ_ed_dash_eq_dash_Idxs_in_state,
     #

     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,

     gens_nodes_idx,
     non_gens_nodes_idx,

     vtf_para_and_idxs,
     vtf_gens_fun_kwd_para,
     ode_per_gen_models_func_kwd_paras,

     init_dyn_pf_flat_para,
     flat_ωs_ωref0_vref0_porder0,

     #
     gens_nodes_τm_vf,
     gens_dynamic_id_iq_pg_vh,
     gens_i_d_0, gens_i_q_0,     
     gens_vd, gens_vq,
     gens_ph, gens_qh,
     #
     
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,

     #
     vec_im_pure_states_Idxs,
     im_pure_states_idxs,
     vec_im_τm_tilade_vf_states_Idxs,
     im_τm_tilade_vf_states_idxs,
     
     vec_im_sauer_states_Idxs,
     flat_im_sauer_states_Idxs,
     vec_im_sauer_δ_eq_dash_ed_dash_Idxs,
     flat_im_sauer_δ_eq_dash_ed_dash_Idxs,
     #
     

     Idx_sets
     )  =
         paras_and_Idxs_set
    
    #----------------------------------------
    #----------------------------------------

    dim_Pg = Idx_sets.sizes_and_dims.dim_Pg
    dim_Qg = Idx_sets.sizes_and_dims.dim_Qg

    dim_Png = Idx_sets.sizes_and_dims.dim_Png
    dim_Qng = Idx_sets.sizes_and_dims.dim_Qng

    dim_gens_id = Idx_sets.sizes_and_dims.dim_gens_id
    dim_gens_iq = Idx_sets.sizes_and_dims.dim_gens_iq
    
    dim_vh = Idx_sets.sizes_and_dims.dim_vh

    dim_θh = Idx_sets.sizes_and_dims.dim_θh

    dim_gens_id = Idx_sets.sizes_and_dims.dim_gens_id
    
    dim_gens_iq = Idx_sets.sizes_and_dims.dim_gens_iq
    
    dim_id_iq_pg_vh = Idx_sets.sizes_and_dims.dim_id_iq_pg_vh
    
    dim_flat_gens_ωs_ωref0_vref0_porder0 =
        Idx_sets.sizes_and_dims.dim_flat_gens_ωs_ωref0_vref0_porder0

    dim_flat_gens_id_iq_pg_vh =
        Idx_sets.sizes_and_dims.dim_flat_gens_id_iq_pg_vh
    
    dim_flat_net_vh_θh =
        Idx_sets.sizes_and_dims.dim_flat_net_vh_θh
    
    dim_init_dyn_pf_flat_para =
        Idx_sets.sizes_and_dims.dim_init_dyn_pf_flat_para
    
    dim_intra_dyn_pf_flat_para =
        Idx_sets.sizes_and_dims.dim_intra_dyn_pf_flat_para

    dim_flat_gens_δ_ed_dash_eq_dash =
        Idx_sets.sizes_and_dims.dim_flat_gens_δ_ed_dash_eq_dash
    
    dim_pure_states = length(
        im_net_states_and_var_labels.gens_nodes_pure_states_labels)
    
    dim_im_states = length(
        im_net_states_and_var_labels.gens_nodes_im_vars_labels)
    
    dim_P_Q =
        2 * Idx_sets.sizes_and_dims.node_size

    
    no_of_gens =
        Idx_sets.sizes_and_dims.no_of_gens

    
    no_of_non_gens =
        Idx_sets.sizes_and_dims.no_of_non_gens
    
    #----------------------------------------    
    
    gens_vh_idx_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ no_of_gens,
              no_of_gens ] ;
            dims_given = true )

    
    gens_id_iq_pg_vh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_id_iq_pg_vh
              for a_gen in
                  gens_nodes_idx  ]
            ; dims_given = true )
    
    
    flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_flat_gens_ωs_ωref0_vref0_porder0,
              dim_flat_gens_id_iq_pg_vh ] ;
            dims_given = true )

    
    flat_vh_flat_θh_flat_id_iq_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_vh,
              dim_θh,
              dim_gens_id,
              dim_gens_iq ] ;
            dims_given = true )
    
    
    flat_vh_flat_θh_Idx =
        flat_vh_flat_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_vh,
              dim_θh ] ;
            dims_given = true )

    
    flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx = get_vars_or_paras_Idxs_in_flattend(
        [dim_flat_net_vh_θh,
         dim_flat_gens_ωs_ωref0_vref0_porder0,
         dim_init_dyn_pf_flat_para];
        dims_given = true )

    
    flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_flat_net_vh_θh,
             dim_flat_gens_ωs_ωref0_vref0_porder0,
             dim_init_dyn_pf_flat_para];
            dims_given = true )

    
    flat_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_flat_gens_ωs_ωref0_vref0_porder0,
             dim_init_dyn_pf_flat_para];
            dims_given = true )

    
    flat_vh_θh_flat_ωs_ωref0_vref0_porder0_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_flat_net_vh_θh,
             dim_flat_gens_ωs_ωref0_vref0_porder0];
            dims_given = true )

    
    flat_vh_θh_flat_init_dyn_pf_para_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_flat_net_vh_θh,
              dim_init_dyn_pf_flat_para];
            dims_given = true )

    
    flat_vh_θh_intra_dyn_pf_para_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_flat_net_vh_θh,
             dim_intra_dyn_pf_flat_para];
            dims_given = true )

    
    flat_vh_θh_flat_δ_ed_dash_eq_dash_init_dyn_pf_para_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_flat_net_vh_θh,
             dim_flat_gens_δ_ed_dash_eq_dash,
             dim_init_dyn_pf_flat_para];
            dims_given = true )

    
    pure_state_net_vh_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_pure_states,
             dim_flat_net_vh_θh];
            dims_given = true )
        
    im_state_net_vh_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_im_states,
             dim_flat_net_vh_θh];
            dims_given = true )
    
    flat_P_Q_flat_δ_ed_dash_eq_dash_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_P_Q,
             dim_flat_gens_δ_ed_dash_eq_dash] ;
            dims_given = true )

    
    flat_δ_ed_dash_eq_dash_intra_dyn_pf_flat_para_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_flat_gens_δ_ed_dash_eq_dash,
             dim_intra_dyn_pf_flat_para];
            dims_given = true )

    
    flat_Pg_Png_Qg_Qng_id_iq_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_Pg,
             dim_Png,
             dim_Qg,
             dim_Qng,
             dim_gens_id,
             dim_gens_iq];
            dims_given = true )

    
    flat_Pg_flat_Qg_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_Pg,
             dim_Qg];
            dims_given = true )
    
    #----------------------------------------------- 
    #-----------------------------------------------
        
    # vh_θh_idx_in_Idx = ur_ui_idx_in_Idx
    
    #---------------------------------------------
    #---------------------------------------------

    nodes_vh_θh_Idx_in_state =
        nodes_ur_ui_Idx_in_state
    
    nodes_vh_Idx_in_state =
        first.(nodes_u_Idx_in_ranges)

    nodes_θh_Idx_in_state =
        last.(nodes_u_Idx_in_ranges)

    flat_vh_idx_in_flat_vh_θh =
        flat_vh_idx_in_Idx =
        first.(vh_θh_idx_in_Idx)

    flat_θh_idx_in_flat_vh_θh =
        flat_θh_idx_in_Idx =
        last.(vh_θh_idx_in_Idx)

    flat_vh_idx_flat_θh_idx_in_flat_vh_θh =
        [flat_vh_idx_in_flat_vh_θh,
         flat_θh_idx_in_flat_vh_θh]

    flat_vh_flat_θh_Idx_in_state =
        [nodes_vh_Idx_in_state;
         nodes_θh_Idx_in_state]

    #----------------------------------------------- 
    # idxs needed for small signal stability analysis
    #----------------------------------------------- 

    gens_vh_θh_idx_in_flat_vh_θh =
        [vh_θh_idx_in_Idx[gens_nodes_idx]...;]
    
    non_gens_vh_θh_idx_in_flat_vh_θh =
        [vh_θh_idx_in_Idx[non_gens_nodes_idx]...;]
    
    consecutive_vh_θh_idxs, vec_a_vh_a_θh_idxs =
        convert_to_consecutive_idxs(
            flat_vh_flat_θh_Idx )

    gens_consecutive_vh_θh_idxs =
        [vec_a_vh_a_θh_idxs[gens_nodes_idx]...;]

    non_gens_consecutive_vh_θh_idxs =
        [vec_a_vh_a_θh_idxs[non_gens_nodes_idx]...;]

    im_state_idx_in_∂f∂x, net_vh_θh_idx_in_∂f∂x =
        im_state_net_vh_θh_idx_in_Idx
         
   vh_θh_idx_in_∂f∂p, ωs_ωref0_vref0_porder0_idx_in_∂f∂p =
       flat_vh_θh_flat_ωs_ωref0_vref0_porder0_Idx

   P_Q_idx_in_∂g∂p, δ_ed_dash_eq_dash_idx_in_∂g∂p =
       flat_P_Q_flat_δ_ed_dash_eq_dash_idx_in_Idx
    
    #----------------------------------------------
    # for non pre ordered nodes
    #----------------------------------------------

    (nodes_flat_vh_Idx,
     nodes_flat_θh_Idx,
     flat_gens_id_Idx,
     flat_gens_iq_Idx)  =
        flat_vh_flat_θh_flat_id_iq_Idx

    gens_vh_idxs =
        nodes_flat_vh_Idx[gens_nodes_idx]
    
    gens_θh_idxs =
        nodes_flat_θh_Idx[gens_nodes_idx]

    non_gens_nodes_vh_idxs =
        nodes_flat_vh_Idx[non_gens_nodes_idx]
    
    non_gens_nodes_θh_idxs =
        nodes_flat_vh_Idx[non_gens_nodes_idx]

    gens_id_idxs = flat_gens_id_Idx
    
    gens_iq_idxs = flat_gens_iq_Idx    
    
    non_pre_ordered_dyn_pf_vars_Idxs =
        (; gens_vh_idxs,
         gens_θh_idxs,

         non_gens_nodes_vh_idxs,
         non_gens_nodes_θh_idxs,

         gens_id_idxs,
         gens_iq_idxs,
         flat_vh_flat_θh_flat_id_iq_Idx)

    #----------------------------------------
    # indices
    #----------------------------------------

   (flat_ωs_ωref0_vref0_porder0_idx_in_∂f∂p,
    flat_id_iq_pg_vh_idx_in_∂f∂p ) =
          flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx

    gens_id_iq_pg_vh_idx_in_∂f∂p =
        collect.([ flat_id_iq_pg_vh_idx_in_∂f∂p[idx]
          for idx in
                   gens_id_iq_pg_vh_idx_in_Idx ])

    gens_id_iq_idx_in_∂f∂p =
        [[[first(idx), second(idx) ]
         for idx in
             gens_id_iq_pg_vh_idx_in_∂f∂p ]...;]

    gens_ph_idx_in_∂f∂p =
        [third(idx)
         for idx in
             gens_id_iq_pg_vh_idx_in_∂f∂p ]

    gens_vh_idx_in_∂f∂p =
        [last(idx)
         for idx in
             gens_id_iq_pg_vh_idx_in_∂f∂p ]

   (flat_vh_idx_in_∂g∂x,
    flat_θh_idx_in_∂g∂x,
    flat_id_idx_in_∂g∂x,
    flat_iq_idx_in_∂g∂x) =
        flat_vh_flat_θh_flat_id_iq_Idx

    flat_vh_flat_θh_idx_in_∂g∂x =
        [ flat_vh_idx_in_∂g∂x,
          flat_θh_idx_in_∂g∂x ]

    flat_id_flat_iq_idx_in_∂g∂x =
        [ flat_id_idx_in_∂g∂x,
          flat_iq_idx_in_∂g∂x ]

    # consecutive_flat_vh_flat_θh_idx_in_∂g∂x

    # consecutive_flat_vh_flat_θh_idx_in_∂g∂x
    # flat_vh_flat_θh_idx_in_∂g∂x
    
    ( consecutive_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
      vec_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0 ) =
            convert_to_consecutive_idxs(
                flat_Pg_flat_Qg_Idx )

    
    ( consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
      vec_flat_vh_flat_θh_idx_in_∂g∂x ) =
            convert_to_consecutive_idxs(
                flat_vh_flat_θh_idx_in_∂g∂x )

    
    gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x =
        [vec_flat_vh_flat_θh_idx_in_∂g∂x[
            gens_nodes_idx ]...;]

    
    non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x =
        [vec_flat_vh_flat_θh_idx_in_∂g∂x[
            non_gens_nodes_idx ]...;]

    
    ( consecutive_id_iq_idx_in_∂g∂x,
      vec_id_iq_idx_in_∂g∂x ) =
            convert_to_consecutive_idxs(
                flat_id_flat_iq_idx_in_∂g∂x )

    
    ( consecutive_gens_vh_idx_θh_idx_in_Idx,
      vec_gens_vh_idx_θh_idx_in_Idx ) =
            convert_to_consecutive_idxs(
                gens_vh_idx_θh_idx_in_Idx )


    ( consecutive_flat_vh_θh_idx_in_∂g∂x,
      vec_flat_vh_θh_idx_in_∂g∂x ) =
            convert_to_consecutive_idxs(
                [ first.(vh_θh_idx_in_Idx),
                  last.(vh_θh_idx_in_Idx) ] )

    
    gens_consecutive_flat_vh_θh_idx_in_∂g∂x =
        [vec_flat_vh_θh_idx_in_∂g∂x[
            gens_nodes_idx ]...;]

    
    non_gens_consecutive_flat_vh_θh_idx_in_∂g∂x =
        [vec_flat_vh_θh_idx_in_∂g∂x[
            non_gens_nodes_idx ]...;]

    
    (Pg_idx_in_mismatch,
     Png_idx_in_mismatch,
     Qg_idx_in_mismatch,
     Qng_idx_in_mismatch,
     id_idx_in_mismatch,
     iq_idx_in_mismatch ) =
         flat_Pg_Png_Qg_Qng_id_iq_Idx

    
    Pg_idx_Qg_idx_in_mismatch =
        [ Pg_idx_in_mismatch,
          Qg_idx_in_mismatch ]

    
    ( consecutive_Pg_idx_Qg_idx_in_mismatch,
      vec_Pg_idx_Qg_idx_in_mismatch ) =
            convert_to_consecutive_idxs(
                Pg_idx_Qg_idx_in_mismatch )

    
    Png_idx_Qng_idx_in_mismatch =
        [ Png_idx_in_mismatch,
          Qng_idx_in_mismatch ]

    
    ( consecutive_Png_idx_Qng_idx_in_mismatch,
      vec_Png_idx_Qng_idx_in_mismatch ) =
            convert_to_consecutive_idxs(
                Png_idx_Qng_idx_in_mismatch )

    
    id_idx_iq_idx_in_mismatch =
        [ id_idx_in_mismatch,
          iq_idx_in_mismatch ]

    
    ( consecutive_id_idx_iq_idx_in_mismatch,
      vec_id_idx_iq_idx_in_mismatch ) =
            convert_to_consecutive_idxs(
                id_idx_iq_idx_in_mismatch )

    #----------------------------------------

    (; Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs

    #----------------------------------------

    (; dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
     ) = dyn_pf_P_Q_δ_etc_kwd_para_Idxs 


    intra_flat_δ_ed_dash_eq_dash_idx_in_∂g∂p =
        dyn_δ_ed_eq_pf_Idxs


    gens_PQ_idx_in_Idxs =
        [ dyn_P_gens_Idxs,
          dyn_Q_gens_Idxs ]

    ( consecutive_dyn_gens_PQ_idx_in_Idxs,
      vec_dyn_gens_PQ_idx_in_Idxs ) =
            convert_to_consecutive_idxs(
                gens_PQ_idx_in_Idxs )


    non_gens_PQ_idx_in_Idxs =
        [ dyn_P_non_gens_Idxs,
          dyn_Q_non_gens_Idxs ]

    ( consecutive_dyn_non_gens_PQ_idx_in_Idxs,
      vec_dyn_non_gens_PQ_idx_in_Idxs ) =
            convert_to_consecutive_idxs(
                non_gens_PQ_idx_in_Idxs )

    if loc_load_exist == true
           
        gens_loc_load_PQ_idx_in_Idxs =
            [ dyn_P_g_loc_load_Idxs,
              dyn_Q_g_loc_load_Idxs ]

        ( consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs,
          vec_dyn_gens_loc_load_PQ_idx_in_Idxs ) =
                convert_to_consecutive_idxs(
                    gens_loc_load_PQ_idx_in_Idxs )
    else
        
        consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs =
            nothing
    end
    
    #----------------------------------------------
    # param
    #----------------------------------------------

    sim_state_x0 = copy(state)
        
    flat_ur_ui_x0 =
        sim_state_x0[ nodes_ur_ui_Idx_in_state ]

    vec_nodes_ur_ui_x0 =
        [ flat_ur_ui_x0[idx]
          for idx in
              ur_ui_idx_in_Idx ]
    
    vec_nodes_vh_θh_x0 =
        cartesian_to_polar.(vec_nodes_ur_ui_x0 )

    flat_vh_θh_x0 = [vec_nodes_vh_θh_x0...;]
    
    
    flat_ur_x0 = first.(vec_nodes_ur_ui_x0)

    flat_ui_x0 = last.(vec_nodes_ur_ui_x0)
    
    flat_ur_flat_ui_x0 = [flat_ur_x0; flat_ui_x0]
    

    flat_vh_x0 = first.(vec_nodes_vh_θh_x0)    

    flat_θh_x0 = last.(vec_nodes_vh_θh_x0)

    flat_vh_flat_θh_x0 = [flat_vh_x0; flat_θh_x0]

        
    sim_state_x0[nodes_ur_ui_Idx_in_state] .=
        flat_vh_θh_x0

    #-----------------------------------------------
    
    stateDiffCache =
        DiffCache(
            similar( sim_state_x0 ))
    
    flat_ur_ui_DiffCache =
        DiffCache(
            similar( flat_ur_ui_x0 ))

    
    flat_ur_flat_ui_DiffCache =
        DiffCache(
            similar( flat_ur_flat_ui_x0 ))
    
    flat_vh_θh_DiffCache =
        DiffCache(
            similar( flat_vh_θh_x0 ))
    
    flat_vh_flat_θh_DiffCache  =
        DiffCache(
            similar(  flat_vh_flat_θh_x0 ))

    _,_, Ax_sparse_nzvalues =
        findnz( sparse( Matrix(vec_Ax_views[ 1 ])))

    Ax_sparse_nzvalues_DiffCache =
        DiffCache( similar( Ax_sparse_nzvalues ) )

    flat_δ_ed_dash_eq_dash_x0 =
        sim_state_x0[
            flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    flat_δ_ed_dash_eq_dash_DiffCache =
        DiffCache( similar(
            flat_δ_ed_dash_eq_dash_x0  ) )
        
    #-----------------------------------------------

    flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh =
        [ flat_ωs_ωref0_vref0_porder0...;
          gens_dynamic_id_iq_pg_vh... ]
    
    flat_ur_ui_wt_init_dyn_pf_flat_para  =
        [flat_ur_ui_x0;
         init_dyn_pf_flat_para ]
    
    flat_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para =
        [flat_ωs_ωref0_vref0_porder0;
         init_dyn_pf_flat_para ]
    
    flat_vh_θh_flat_ωs_ωref0_vref0_porder0 =
        [flat_vh_θh_x0;
         flat_ωs_ωref0_vref0_porder0 ]

    flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0 =
        [flat_vh_flat_θh_x0;
         flat_ωs_ωref0_vref0_porder0 ]

    flat_vh_θh_init_dyn_pf_flat_para =
        [flat_vh_θh_x0 ;
         init_dyn_pf_flat_para ]

    flat_vh_flat_θh_init_dyn_pf_flat_para =
        [flat_vh_flat_θh_x0 ;
         init_dyn_pf_flat_para ]
    
    flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para =
        [flat_vh_θh_x0;
         flat_ωs_ωref0_vref0_porder0;
         init_dyn_pf_flat_para ]    

    flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para =
        [flat_vh_flat_θh_x0;
         flat_ωs_ωref0_vref0_porder0;
         init_dyn_pf_flat_para ]

    flat_vh_flat_θh_flat_δ_ed_dash_eq_dash_init_dyn_pf_flat_para =
        [flat_vh_flat_θh_x0;
         flat_δ_ed_dash_eq_dash_x0;
         init_dyn_pf_flat_para ]

    #----------------------------------------------
    #----------------------------------------------
    
    # vh_θh_idx_order = :consecutive
    
    vh_θh_idx_order = :flat_vh_first_flat_θh_last
    
    Pg_Qg_external_control =
        false

    use_nlsolve =
        false

    #----------------------------------------------
    
       pf_solver =
        (; pf_alg,
         abstol,
         reltol )
    
    #----------------------------------------------

    counter_array = [1]
    
    #----------------------------------------------

    init_dyn_pf_mismatch_kwd_para =
        (; loc_load_exist,
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,     
         dyn_pf_vars_Idxs,
         dyn_pf_Idxs_kwd_para,
         gens_nodes_ra_Xd_dash_Xq_dash,
         non_pre_ordered_dyn_pf_vars_Idxs
         )
    
    init_pf_by_vh_θh_by_parts_kwd =
        (;
         pf_solver,
         use_nlsolve,
         init_dyn_pf_mismatch_kwd_para
         )

    init_pf_by_vh_θh_kwd =
        (;
         flat_vh_θh_flat_init_dyn_pf_para_Idx,
         pf_solver,
         use_nlsolve,
         init_dyn_pf_mismatch_kwd_para
         )

    #----------------------------------------------

    intra_dyn_pf_mismatch_kwd_para =
        (;
         loc_load_exist,
         dyn_pf_vars_Idxs,
         dyn_pf_Idxs_kwd_para,
         gens_nodes_ra_Xd_dash_Xq_dash,
         non_pre_ordered_dyn_pf_vars_Idxs
         )
    # flat_vh_flat_θh_idx_in_Idx
    intra_pf_by_vh_θh_by_parts_kwd =
        (;         
         Pg_Qg_external_control,
         loc_load_exist,

         gens_nodes_ra_Xd_dash_Xq_dash,

         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,

         vh_θh_idx_order,
         vh_θh_idx_in_Idx,
         flat_vh_flat_θh_Idx,

         gens_nodes_idx,
         
         δ_ed_dash_eq_dash_Idxs_in_flattend,
         
         intra_dyn_pf_mismatch_kwd_para,
         pf_solver,
         use_nlsolve
         )
    
    intra_pf_by_vh_θh_kwd =
        (;
         flat_vh_θh_flat_δ_ed_dash_eq_dash_init_dyn_pf_para_Idx,
         
         Pg_Qg_external_control,
         loc_load_exist,

         gens_nodes_ra_Xd_dash_Xq_dash,

         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,

         vh_θh_idx_order,
         vh_θh_idx_in_Idx,
         flat_vh_flat_θh_Idx,

         gens_nodes_idx,
         
         δ_ed_dash_eq_dash_Idxs_in_flattend,
         
         intra_dyn_pf_mismatch_kwd_para,
         pf_solver,
         use_nlsolve
         )

    intra_dyn_pf_kwd_para =
        (;                  
         Pg_Qg_external_control,
         loc_load_exist,

         gens_nodes_ra_Xd_dash_Xq_dash,

         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,

         vh_θh_idx_order,
         vh_θh_idx_in_Idx,

         flat_vh_flat_θh_Idx,

         gens_nodes_idx,

         δ_ed_dash_eq_dash_Idxs_in_flattend,

         intra_dyn_pf_mismatch_kwd_para,
         use_nlsolve,
         pf_solver

         )

    intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd =
        (;
         flat_vh_flat_θh_flat_id_iq_Idx,

         intra_dyn_pf_kwd_para
         )
    
    #-----------------------------------------------
    
    ode_vtf_by_vh_θh_para_kwd =
        (;
         gens_nodes_ra_Xd_dash_Xq_dash,
         gens_Ax_update_parameters,
         ode_per_gen_models_func_kwd_paras,
         vtf_gens_fun_kwd_para
          )
    
    ode_vtf_by_vh_θh_Idxs_kwd =
        (;
         flat_vh_θh_flat_ωs_ωref0_vref0_porder0_Idx,

         flat_vh_idx_in_Idx,
         flat_θh_idx_in_Idx,
         
         flat_vh_flat_θh_Idx,

         gens_ωs_ωref0_vref0_porder0_idx_in_Idx,

         gens_nodes_idx,
         non_gens_nodes_idx,
         
         nodes_state_Idx,

         gens_nodes_u_Idx_in_ranges,
         non_gens_nodes_u_Idx_in_ranges,

         nodes_vh_θh_Idx_in_state,
         im_vars_Idx_in_state,
         each_gens_im_vars_Idx_in_state,

         nodes_δ_ed_dash_eq_dash_Idxs,
         flat_δ_ed_dash_eq_dash_Idxs_in_state )

    ode_vtf_by_vh_θh_kwd_para =
        (;
         ode_vtf_by_vh_θh_para_kwd,
         ode_vtf_by_vh_θh_Idxs_kwd )

    ode_only_per_gen_id_iq_pg_vh_kwd =
        (;
         flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,
         gens_id_iq_pg_vh_idx_in_Idx,
         ode_vtf_by_vh_θh_kwd_para )

    #-----------------------------------------------
    #-----------------------------------------------

    sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para =
        (;
         post_pf_idxs,
         counter_array,

         sim_state_x0,
         flat_vh_flat_θh_x0,
         flat_δ_ed_dash_eq_dash_x0,

         stateDiffCache,
         flat_vh_flat_θh_DiffCache,     
         Ax_sparse_nzvalues_DiffCache,
         flat_δ_ed_dash_eq_dash_DiffCache,

         flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx,
         init_pf_by_vh_θh_by_parts_kwd,

         ode_vtf_by_vh_θh_kwd_para )
    
    #----------------------------------------------
    #----------------------------------------------
    
    sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para =
        (;
         flat_vh_flat_θh_flat_id_iq_Idx,
         
         post_pf_idxs,
         counter_array,

         sim_state_x0,
         flat_vh_flat_θh_x0,
         flat_δ_ed_dash_eq_dash_x0,

         stateDiffCache,
         flat_vh_flat_θh_DiffCache,     
         Ax_sparse_nzvalues_DiffCache,
         flat_δ_ed_dash_eq_dash_DiffCache,

         flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx,
         init_pf_by_vh_θh_by_parts_kwd,

         intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd,

         ode_only_per_gen_id_iq_pg_vh_kwd
         )

    
    sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd =
        (;
         
         sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para,
     
         # flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,

         sim_state_x0,
         im_mass_matrix,
         im_sym,
         sim_timespan,
         ode_alg,
         dt )

    sd_per_gen_ode_vtf_pf_by_flat_vh_θh_kwd =
        (;
         sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para,
         # sd_per_gen_ode_vtf_pf_kwd_para,

         # flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,

         sim_state_x0,
         im_mass_matrix,
         im_sym,     
         ode_alg,
         dt,
         sim_timespan
         )


    ode_vtf_by_vh_θh_kwd =
        (;
         ode_vtf_by_vh_θh_kwd_para,
         # flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0,

         im_mass_matrix,
         im_sym,
         sim_state_x0,
         sim_timespan,
         ode_alg,
         dt
         ) 

    ode_vtf_by_flat_vh_θh_kwd =
        (;
         ode_vtf_by_vh_θh_kwd_para,
         # flat_vh_θh_flat_ωs_ωref0_vref0_porder0,

         im_mass_matrix,
         im_sym,
         sim_state_x0,
         sim_timespan,
         ode_alg,
         dt
         )

    ode_vtf_by_flat_vh_flat_θh_kwd =
        (;
         ode_vtf_by_vh_θh_kwd_para,
         # flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0,

         im_mass_matrix,
         im_sym,
         sim_state_x0,
         ode_alg,
         dt
         )

    intra_dyn_pf_by_vh_θh_δ_ed_eq_by_parts_sense_kwd =
        (;
         sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd,
         # flat_vh_flat_θh_x0,
         
         flat_δ_ed_dash_eq_dash_Idxs_in_state,
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         
         # init_dyn_pf_flat_para,
         intra_pf_by_vh_θh_by_parts_kwd,
         intra_dyn_pf_mismatch_kwd_para
         )

    init_dyn_pf_by_vh_θh_by_parts_sense_kwd =
        (;
         # flat_vh_flat_θh_x0,
         # init_dyn_pf_flat_para,
         init_pf_by_vh_θh_by_parts_kwd,
         init_dyn_pf_mismatch_kwd_para     
         )

    intra_dyn_pf_by_vh_θh_by_parts_sense_kwd =
        (;
         sim_timespan,
         sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd ,

         flat_vh_flat_θh_x0,
         flat_δ_ed_dash_eq_dash_Idxs_in_state, 
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,

         init_dyn_pf_flat_para,
         intra_pf_by_vh_θh_by_parts_kwd,
         intra_dyn_pf_mismatch_kwd_para
         )


    intra_dyn_current_balance_pf_by_vh_θh_sense_kwd =
        (;
         sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd,
         flat_vh_flat_θh_x0,
         flat_δ_ed_dash_eq_dash_Idxs_in_state,
         
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,

         init_dyn_pf_flat_para,
         intra_pf_by_vh_θh_by_parts_kwd,
         intra_dyn_pf_mismatch_kwd_para
         )

    stability_sd_dynamics_wt_intra_dyn_pf_kwd =
        (;
         loc_load_exist,
         sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para,

         # flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,
         # flat_ωs_ωref0_vref0_porder0,

         sim_state_x0,
         # flat_vh_flat_θh_x0,
         # init_dyn_pf_flat_para,
         gens_nodes_ra_Xd_dash_Xq_dash,

         im_mass_matrix,
         im_sym, 
         ode_alg,
         dt,

         flat_vh_flat_θh_Idx,

         gens_nodes_idx,
         non_gens_nodes_idx,

         flat_δ_ed_dash_eq_dash_Idxs_in_state,
         δ_ed_dash_eq_dash_Idxs_in_flattend,

         intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd,
         intra_dyn_pf_mismatch_kwd_para,

         im_state_idx_in_∂f∂x,
         net_vh_θh_idx_in_∂f∂x,

         consecutive_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
         vec_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
         

         #
         vec_flat_vh_flat_θh_idx_in_∂g∂x,
         vec_flat_vh_θh_idx_in_∂g∂x,
         vec_id_iq_idx_in_∂g∂x,

         consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
         consecutive_flat_vh_θh_idx_in_∂g∂x,
         #
         
         gens_consecutive_flat_vh_θh_idx_in_∂g∂x,
         non_gens_consecutive_flat_vh_θh_idx_in_∂g∂x,

         gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
         non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
         consecutive_id_iq_idx_in_∂g∂x,

         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         dyn_pf_P_Q_δ_etc_kwd_para_Idxs,

         consecutive_Pg_idx_Qg_idx_in_mismatch,
         consecutive_Png_idx_Qng_idx_in_mismatch,      
         consecutive_id_idx_iq_idx_in_mismatch,

         consecutive_dyn_gens_PQ_idx_in_Idxs,
         consecutive_dyn_non_gens_PQ_idx_in_Idxs,
         consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs,
         dyn_δ_ed_eq_pf_Idxs,   

         # ∂f∂x_xp_by_vh_θh_∂f∂x,

         im_pure_states_idxs,
         im_τm_tilade_vf_states_idxs,
         )


    stability_ode_per_gen_wt_intra_dyn_pf_kwd =
        (; loc_load_exist,
         
         sim_state_x0,
         #flat_vh_flat_θh_x0,

         im_vars_Idx_in_state,

         ode_only_per_gen_id_iq_pg_vh_kwd,

         # flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh,

         # flat_ωs_ωref0_vref0_porder0,

         im_ode_mass_matrix,
         im_ode_sym,

         gens_nodes_ra_Xd_dash_Xq_dash,


         flat_vh_flat_θh_Idx,

         gens_nodes_idx,
         non_gens_nodes_idx,

         flat_δ_ed_dash_eq_dash_Idxs_in_state,
         δ_ed_dash_eq_dash_Idxs_in_flattend,

         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         dyn_pf_P_Q_δ_etc_kwd_para_Idxs,

         # init_dyn_pf_flat_para,
         intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd,
         intra_dyn_pf_mismatch_kwd_para,

         im_state_idx_in_∂f∂x,

         flat_id_iq_pg_vh_idx_in_∂f∂p,
         flat_ωs_ωref0_vref0_porder0_idx_in_∂f∂p,
         gens_id_iq_idx_in_∂f∂p,
         gens_vh_idx_in_∂f∂p,
         
         gens_ph_idx_in_∂f∂p,
         
         consecutive_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
         vec_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,         

         #
         vec_flat_vh_flat_θh_idx_in_∂g∂x,
         vec_flat_vh_θh_idx_in_∂g∂x,
         vec_id_iq_idx_in_∂g∂x,

         consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
         consecutive_flat_vh_θh_idx_in_∂g∂x,
         consecutive_id_iq_idx_in_∂g∂x,
         #         
         
         consecutive_gens_vh_idx_θh_idx_in_Idx,

         gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
         
         non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,

         intra_flat_δ_ed_dash_eq_dash_idx_in_∂g∂p,

         consecutive_Pg_idx_Qg_idx_in_mismatch,
         consecutive_Png_idx_Qng_idx_in_mismatch,
         consecutive_id_idx_iq_idx_in_mismatch,

         im_pure_states_idxs,
         im_τm_tilade_vf_states_idxs
         )

    vars =
        (;
         state,
         sim_state_x0,
         flat_ur_ui_x0,
         flat_vh_θh_x0,
         flat_ur_flat_ui_x0,
         flat_vh_flat_θh_x0,
         flat_δ_ed_dash_eq_dash_x0,
         init_dyn_pf_flat_para,
         flat_ωs_ωref0_vref0_porder0,
         vec_nodes_ur_ui_x0,
         vec_nodes_vh_θh_x0,
         gens_dynamic_id_iq_pg_vh)

    return (;
            vars,
            
            init_dyn_pf_mismatch_kwd_para,

            init_pf_by_vh_θh_by_parts_kwd,

            init_pf_by_vh_θh_kwd,

            intra_dyn_pf_mismatch_kwd_para,

            intra_pf_by_vh_θh_by_parts_kwd,

            intra_pf_by_vh_θh_kwd,

            intra_dyn_pf_kwd_para,

            intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd,

            ode_vtf_by_vh_θh_para_kwd,

            ode_vtf_by_vh_θh_Idxs_kwd,

            ode_vtf_by_vh_θh_kwd_para,

            ode_only_per_gen_id_iq_pg_vh_kwd,

            sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para,

            sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para,

            sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd,

             sd_per_gen_ode_vtf_pf_by_flat_vh_θh_kwd,

            ode_vtf_by_vh_θh_kwd,

            ode_vtf_by_flat_vh_θh_kwd,

            ode_vtf_by_flat_vh_flat_θh_kwd,

            intra_dyn_pf_by_vh_θh_δ_ed_eq_by_parts_sense_kwd,

            init_dyn_pf_by_vh_θh_by_parts_sense_kwd,

            intra_dyn_pf_by_vh_θh_by_parts_sense_kwd,

            intra_dyn_current_balance_pf_by_vh_θh_sense_kwd,

            stability_sd_dynamics_wt_intra_dyn_pf_kwd,

            stability_ode_per_gen_wt_intra_dyn_pf_kwd
            )    
    

end



function simulate_sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_flat_θh(flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para; kwd = sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd  )


     (; sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para,

      # flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,

      sim_state_x0,
      im_mass_matrix,
      im_sym,
      sim_timespan,
      ode_alg,
      dt
      ) = kwd
    
    #------------------------------------------------
    #------------------------------------------------
    # sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_func!
    #------------------------------------------------
    #------------------------------------------------

    sd_per_gen_ode_vtf_pf_kwd_para =
        sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para
    
    sim_func  =
        sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_func!
    
    #------------------------------------------------

    ode_vtf_dynamics_sim_sol = DifferentialEquations.solve(
        ODEProblem(
            ODEFunction{true}(
                (dx, x, p, t) ->
                    sim_func(dx, x, p, t;
                             kwd_para =
                                 sd_per_gen_ode_vtf_pf_kwd_para
                         );
                    mass_matrix =
                        im_mass_matrix ,
                    syms =
                        im_sym ),
                sim_state_x0 ,
                sim_timespan,
                sim_fun_para ),
            ode_alg,
            dt=dt)
        

    
end




function simulate_sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh(flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para ; kwd = sd_per_gen_ode_vtf_pf_by_flat_vh_θh_kwd )

     (;
      sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para,
      # sd_per_gen_ode_vtf_pf_kwd_para,

      sim_state_x0,
      im_mass_matrix,
      im_sym,     
      ode_alg,
      dt,
      sim_timespan
      ) = kwd
    
    #----------------------------------------------
    #----------------------------------------------
    # sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    #----------------------------------------------
    #----------------------------------------------
    
    sd_per_gen_ode_vtf_pf_kwd_para =
        sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para 
    
    sim_func  =
        sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    
    #---------------------------------------------

     sd_dynamics_per_gen_ode_vtf_pf_sim_sol =
         DifferentialEquations.solve(
             ODEProblem(
                 ODEFunction{true}(
                     (dx, x, p, t) ->
                         sim_func(dx, x, p, t;
                                  kwd_para =
                                      sd_per_gen_ode_vtf_pf_kwd_para
                              );
                         mass_matrix =
                             im_mass_matrix ,
                         syms =
                             im_sym ),
                     sim_state_x0 ,
                     sim_timespan,
                     sim_fun_para ),
                 ode_alg,
                 dt=dt )
    
    #------------------------------------------------
    # plot
    #------------------------------------------------

    plot_idxs_a =
        get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = sd_dynamics_per_gen_ode_vtf_pf_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = sd_dynamics_per_gen_ode_vtf_pf_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = sd_dynamics_per_gen_ode_vtf_pf_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            # vars = [:δ, :ω, :ed_dash, :eq_dash ],
            vars = [:δ, :ω ],
            tspan = sim_timespan,
            fmt = :png)
    
    #-----------------------------------------------
    # sensitivity
    #-----------------------------------------------

    sd_dynamics_per_gen_ode_vtf_pf_prob =
        ODEForwardSensitivityProblem(
            (dx, x, p, t) ->
                sim_func(
                    dx, x, p, t;
                    kwd_para =
                        sd_per_gen_ode_vtf_pf_kwd_para ),
            sim_state_x0,
            sim_timespan,
            sim_fun_para;
            sensealg =
                ForwardDiffSensitivity() )

    sd_dynamics_per_gen_ode_vtf_pf_probsim_sol =
        DifferentialEquations.solve(
            sd_dynamics_per_gen_ode_vtf_pf_prob,
            ode_alg, dt = dt )
    
    #---------------------------------------------
    # sensitivity
    #---------------------------------------------

    tt   = sim_timespan[1]
    sd_x = sd_dynamics_per_gen_ode_vtf_pf_sim_sol( tt )
   sd_dx = similar(sd_x)
    
    sd_per_gen_ode_vtf_pf_kwd_para =
        sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para 
    
    sim_func  =
        sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    
    #---------------------------------------------
        
    sd_per_gen_ode_vtf_pf_∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> sim_func(
            dx, x, sim_fun_para, tt;
            kwd_para = sd_per_gen_ode_vtf_pf_kwd_para ),
        sd_dx, sd_x )
    
    sd_per_gen_ode_vtf_pf_∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> sim_func(
            dx, sd_x, p, tt;
            kwd_para = sd_per_gen_ode_vtf_pf_kwd_para ),
        sd_dx, sim_fun_para )

    sd_per_gen_ode_vtf_pf_dxdp =
        -(svd(sd_per_gen_ode_vtf_pf_∂f∂x)\
        sd_per_gen_ode_vtf_pf_∂f∂p)
            
end




function simulate_ode_vtf_by_vh_θh(
    flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0;
    kwd = ode_vtf_by_vh_θh_kwd )

   (;
    ode_vtf_by_vh_θh_kwd_para,

    im_mass_matrix,
    im_sym,
    sim_state_x0,
    sim_timespan,
    ode_alg,
    dt
    ) = kwd
    
    #------------------------------------------------
    # ode_vtf_by_vh_θh_func!
    #------------------------------------------------
    
    sim_fun_ode_vtf_by_vh_θh_kwd_para =
        ode_vtf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0 
    
    sim_func  =
        ode_vtf_by_vh_θh_func!
    
    #-----------------------------------------------

    ode_vtf_dynamics_sim_sol = DifferentialEquations.solve(
        ODEProblem(
            ODEFunction{true}(
                (dx, x, p, t) ->
                    sim_func(dx, x, p, t;
                             kwd_para =
                                 sim_fun_ode_vtf_by_vh_θh_kwd_para 
                         );
                    mass_matrix =
                        im_mass_matrix ,
                    syms =
                        im_sym),
                sim_state_x0 ,
                sim_timespan,
                sim_fun_para ),
            ode_alg,
            dt=dt)
    
    #---------------------------------------------
    # sensitivity
    #---------------------------------------------

    tt = sim_timespan[1]
    sd_x = ode_vtf_dynamics_sim_sol( tt )
   sd_dx = similar(sd_x)
    
    
    sim_fun_ode_vtf_by_vh_θh_kwd_para =
        ode_vtf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0 

    
    sim_func  =
        ode_vtf_by_flat_vh_flat_θh_func!
    
    #---------------------------------------------

    ode_vtf_by_flat_vh_flat_θh_∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> sim_func(
            dx, x, sim_fun_para, tt;
            kwd_para = sim_fun_ode_vtf_by_vh_θh_kwd_para ),
        sd_dx, sd_x )

    
    ode_vtf_by_flat_vh_flat_θh_∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> sim_func(
            dx, sd_x, p, tt;
            kwd_para = sim_fun_ode_vtf_by_vh_θh_kwd_para ),
        sd_dx, sim_fun_para )

    
end




function simulate_ode_vtf_by_flat_vh_θh(
    flat_vh_θh_flat_ωs_ωref0_vref0_porder0;
    kwd = ode_vtf_by_flat_vh_θh_kwd )

    (;
     ode_vtf_by_vh_θh_kwd_para,
     

     im_mass_matrix,
     im_sym,
     sim_state_x0,
     sim_timespan,
     ode_alg,
     dt,

     ) = kwd
    
    #---------------------------------------------
    # ode_vtf_by_flat_vh_θh_func!
    #----------------------------------------------
    
    sim_fun_ode_vtf_by_vh_θh_kwd_para =
        ode_vtf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_θh_flat_ωs_ωref0_vref0_porder0 
    
    sim_func  =
        ode_vtf_by_flat_vh_θh_func!

    #----------------------------------------------

    ode_vtf_by_flat_vh_θh_sim_sol = DifferentialEquations.solve(
        ODEProblem(
            ODEFunction{true}(
                (dx, x, p, t) ->
                    sim_func(dx, x, p, t;
                             kwd_para =
                                 sim_fun_ode_vtf_by_vh_θh_kwd_para 
                         );
                    mass_matrix =
                        im_mass_matrix ,
                    syms =
                        im_sym ),
                sim_state_x0 ,
                sim_timespan,
                sim_fun_para ),
            ode_alg,
            dt=dt)

    #---------------------------------------------
    # sensitivity
    #---------------------------------------------

    tt = sim_timespan[1]
    sd_x = ode_vtf_by_flat_vh_θh_sim_sol( tt )
   sd_dx = similar(sd_x)
    
    sim_fun_ode_vtf_by_vh_θh_kwd_para =
        ode_vtf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_θh_flat_ωs_ωref0_vref0_porder0 
    
    sim_func  =
        ode_vtf_by_flat_vh_θh_func!
    
    #---------------------------------------------

    ode_vtf_by_flat_vh_θh_∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> sim_func(
            dx, x, sim_fun_para, tt;
            kwd_para = sim_fun_ode_vtf_by_vh_θh_kwd_para ),
        sd_dx,
        sd_x )
    
    ode_vtf_by_flat_vh_θh_∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> sim_func(
            dx, sd_x, p, tt;
            kwd_para = sim_fun_ode_vtf_by_vh_θh_kwd_para ),
        sd_dx,
        sim_fun_para )

    ode_vtf_by_flat_vh_θh_dxdp =
        -(svd(ode_vtf_by_flat_vh_θh_∂f∂x)\
        ode_vtf_by_flat_vh_θh_∂f∂p)
        
end



function simulate_ode_vtf_by_flat_vh_flat_θh(
    flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0;
    kwd = ode_vtf_by_flat_vh_flat_θh_kwd )

    (;
     ode_vtf_by_vh_θh_kwd_para,
     

     im_mass_matrix,
     im_sym,
     sim_state_x0,
     ode_alg,
     dt
     ) = kwd
    
    #-----------------------------------------------
    # ode_vtf_by_flat_vh_flat_θh_func!
    #-----------------------------------------------
    
    sim_fun_ode_vtf_by_vh_θh_kwd_para =
        ode_vtf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0 
    
    sim_func  =
        ode_vtf_by_flat_vh_flat_θh_func!

    ode_vtf_by_flat_vh_flat_θh_sim_sol =
        DifferentialEquations.solve(
            ODEProblem(
                ODEFunction{true}(
                    (dx, x, p, t) ->
                        sim_func(dx, x, p, t;
                                 kwd_para =
                                     sim_fun_ode_vtf_by_vh_θh_kwd_para 
                             );
                    mass_matrix =
                        im_mass_matrix ,
                    syms =
                        im_sym ),
                sim_state_x0 ,
                sim_timespan,
                sim_fun_para ),
            ode_alg,
            dt=dt)

    #---------------------------------------------
    # sensitivity
    #---------------------------------------------

    tt   = sim_timespan[1]
    sd_x = ode_vtf_by_flat_vh_flat_θh_sim_sol( tt )
   sd_dx = similar(sd_x)
        
    sim_fun_ode_vtf_by_vh_θh_kwd_para =
        ode_vtf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0 
    
    #---------------------------------------------
        
    ode_vtf_by_flat_vh_flat_θh_∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> sim_func(
            dx, x, sim_fun_para, tt;
            kwd_para =
                sim_fun_ode_vtf_by_vh_θh_kwd_para ),
        sd_dx, sd_x )
    
    ode_vtf_by_flat_vh_flat_θh_∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> sim_func(
            dx, sd_x, p, tt;
            kwd_para =
                sim_fun_ode_vtf_by_vh_θh_kwd_para ),
        sd_dx, sim_fun_para )

    ode_vtf_by_flat_vh_flat_θh_dxdp =
        -(svd(ode_vtf_by_flat_vh_flat_θh_∂f∂x)\
        ode_vtf_by_flat_vh_flat_θh_∂f∂p)
    
    
end



function get_intra_dyn_pf_by_vh_θh_δ_ed_eq_by_parts_sensitivity(init_flat_vh_flat_θh, flat_δ_ed_dash_eq_dash, init_dyn_pf_flat_para ; kwd = intra_dyn_pf_by_vh_θh_δ_ed_eq_by_parts_sense_kwd )

    (;
     sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd,
     
     flat_δ_ed_dash_eq_dash_Idxs_in_state,     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     
     intra_pf_by_vh_θh_by_parts_kwd,
     intra_dyn_pf_mismatch_kwd_para
     ) = kwd
    
    #---------------------------------------------
    # intra by_vh_θh_δ_ed_eq sensitivity
    #---------------------------------------------

    (;
        Pg_Idx,
        Qg_Idxs,
        Png_Idxs,
        Qng_Idxs,
        Pgll_Idxs,
        Qgll_Idxs ) =
            Pg_Qg_Png_Qng_Pgll_Qgll_Idxs

    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    if loc_load_exist == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash;
                P_g_loc_load;
                Q_g_loc_load ]                  

    else
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]   

    end

    intra_by_vh_θh_δ_ed_eq_pf_sol =
        intra_dyn_pf_by_vh_θh_δ_ed_eq_by_parts_func!(
            init_flat_vh_flat_θh,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_pf_by_vh_θh_by_parts_kwd )    

    intra_by_vh_θh_δ_ed_eq_pf_x  =
        intra_by_vh_θh_δ_ed_eq_pf_sol.u

    intra_by_vh_θh_δ_ed_eq_pf_ΔPQ =
        similar( intra_by_vh_θh_δ_ed_eq_pf_x )


    intra_by_vh_θh_δ_ed_eq_pf_fun =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_mismatch


    intra_by_vh_θh_δ_ed_eq_dyn_sd_∂g∂x = ForwardDiff.jacobian(
        ( pf_ΔPQ, x ) ->
            intra_by_vh_θh_δ_ed_eq_pf_fun(
            pf_ΔPQ, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_by_vh_θh_δ_ed_eq_pf_ΔPQ,
        intra_by_vh_θh_δ_ed_eq_pf_x )


    intra_by_vh_θh_δ_ed_eq_dyn_sd_∂g∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ,  p ) -> intra_by_vh_θh_δ_ed_eq_pf_fun(
            pf_ΔPQ, intra_by_vh_θh_δ_ed_eq_pf_x, p;
            intra_dyn_pf_mismatch_kwd_para =
                intra_dyn_pf_mismatch_kwd_para ),
        intra_by_vh_θh_δ_ed_eq_pf_ΔPQ,
        intra_dyn_pf_mismatch_flat_para )


       intra_by_vh_θh_δ_ed_eq_dyn_sd_dxdp =
           -(svd( intra_by_vh_θh_δ_ed_eq_dyn_sd_∂g∂x )\
           intra_by_vh_θh_δ_ed_eq_dyn_sd_∂g∂p)

end




function get_intra_dyn_pf_by_vh_θh_δ_ed_eq_sensitivity(
    flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para ; kwd = intra_dyn_pf_by_vh_θh_δ_ed_eq_by_parts_sense_kwd )

    (;
     sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd,
     
     flat_δ_ed_dash_eq_dash_Idxs_in_state,     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     
     intra_pf_by_vh_θh_by_parts_kwd,
     intra_dyn_pf_mismatch_kwd_para
     ) = kwd

    (;

     sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para,

     # flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,

     sim_state_x0,
     im_mass_matrix,
     im_sym,
     sim_timespan,
     ode_alg,
     dt ) =
         sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd
    
    #-----------------------------------------------

    flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx =
        sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para.flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx

    flat_vh_θh_idx, flat_ωs_ωref0_vref0_porder0_idx, flat_init_dyn_pf_para_idx = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx

   
    init_flat_vh_flat_θh =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[flat_vh_θh_idx]
    
    flat_ωs_ωref0_vref0_porder0 =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_init_dyn_pf_para_idx ]

    
    #----------------------------------------------
    #----------------------------------------------
    # sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    #----------------------------------------------
    #----------------------------------------------
    
    sd_per_gen_ode_vtf_pf_kwd_para =
        sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para 
    
    sim_func  =
        sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    
    #----------------------------------------------

     sd_dynamics_per_gen_ode_vtf_pf_sim_sol =
         DifferentialEquations.solve(
             ODEProblem(
                 ODEFunction{true}(
                     (dx, x, p, t) ->
                         sim_func(dx, x, p, t;
                                  kwd_para =
                                      sd_per_gen_ode_vtf_pf_kwd_para
                              );
                         mass_matrix =
                             im_mass_matrix ,
                         syms =
                             im_sym ),
                     sim_state_x0 ,
                     sim_timespan,
                     sim_fun_para ),
                 ode_alg,
                 dt=dt )
        
    #---------------------------------------------
    # intra sensitivity
    #---------------------------------------------

    tt = sim_timespan[1]
    sd_x = sd_dynamics_per_gen_ode_vtf_pf_sim_sol( tt )
    sd_dx = similar(sd_x)
    
    flat_δ_ed_dash_eq_dash =
        sd_x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    flat_δ_ed_dash_eq_dash_x0 =
        sim_state_x0[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    #---------------------------------------------
    # intra by_vh_θh_δ_ed_eq sensitivity
    #---------------------------------------------

    (;
        Pg_Idx,
        Qg_Idxs,
        Png_Idxs,
        Qng_Idxs,
        Pgll_Idxs,
        Qgll_Idxs ) =
            Pg_Qg_Png_Qng_Pgll_Qgll_Idxs

    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    if loc_load_exist == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash;
                P_g_loc_load;
                Q_g_loc_load ]                  

    else
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]   

    end

    intra_by_vh_θh_δ_ed_eq_pf_sol =
        intra_dyn_pf_by_vh_θh_δ_ed_eq_by_parts_func!(
            init_flat_vh_flat_θh,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_pf_by_vh_θh_by_parts_kwd )    

    intra_by_vh_θh_δ_ed_eq_pf_x  =
        intra_by_vh_θh_δ_ed_eq_pf_sol.u

    intra_by_vh_θh_δ_ed_eq_pf_ΔPQ =
        similar( intra_by_vh_θh_δ_ed_eq_pf_x )


    intra_by_vh_θh_δ_ed_eq_pf_fun =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_mismatch


    intra_by_vh_θh_δ_ed_eq_dyn_sd_∂g∂x = ForwardDiff.jacobian(
        ( pf_ΔPQ, x ) ->
            intra_by_vh_θh_δ_ed_eq_pf_fun(
            pf_ΔPQ, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_by_vh_θh_δ_ed_eq_pf_ΔPQ,
        intra_by_vh_θh_δ_ed_eq_pf_x )


    intra_by_vh_θh_δ_ed_eq_dyn_sd_∂g∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ,  p ) -> intra_by_vh_θh_δ_ed_eq_pf_fun(
            pf_ΔPQ, intra_by_vh_θh_δ_ed_eq_pf_x, p;
            intra_dyn_pf_mismatch_kwd_para =
                intra_dyn_pf_mismatch_kwd_para ),
        intra_by_vh_θh_δ_ed_eq_pf_ΔPQ,
        intra_dyn_pf_mismatch_flat_para )


       intra_by_vh_θh_δ_ed_eq_dyn_sd_dxdp =
           -(svd( intra_by_vh_θh_δ_ed_eq_dyn_sd_∂g∂x )\
           intra_by_vh_θh_δ_ed_eq_dyn_sd_∂g∂p)

end




function get_init_dyn_pf_by_vh_θh_by_parts_sensitivity(
    flat_vh_flat_θh_x0,
    init_dyn_pf_flat_para;
    kwd = init_dyn_pf_by_vh_θh_by_parts_sense_kwd )

    (;
     init_pf_by_vh_θh_by_parts_kwd,
     init_dyn_pf_mismatch_kwd_para     
     ) = kwd

    #---------------------------------------------
    # init sensitivity
    #---------------------------------------------

    init_flat_vh_flat_θh =
        flat_vh_flat_θh_x0

    init_pf_sol =
        init_dyn_pf_by_vh_θh_by_parts_func!(
            init_flat_vh_flat_θh,
            init_dyn_pf_flat_para;
            init_pf_by_vh_θh_by_parts_kwd =
                init_pf_by_vh_θh_by_parts_kwd )

    init_pf_x   = init_pf_sol.u
    init_pf_ΔPQ = similar( init_pf_x )

    init_pf_fun  =
        get_a_model_integrated_init_dyn_pf_ΔPQ_mismatch
    
    init_dyn_sd_∂f∂x = ForwardDiff.jacobian(
        (pf_ΔPQ, x ) ->
            init_pf_fun(
            pf_ΔPQ, x, init_dyn_pf_flat_para;
            init_dyn_pf_mismatch_kwd_para =
                init_dyn_pf_mismatch_kwd_para ),
        init_pf_ΔPQ, init_pf_x )

    init_dyn_sd_∂f∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ, p ) ->
            init_pf_fun(
            pf_ΔPQ, init_pf_x, p;
            init_dyn_pf_mismatch_kwd_para =
                init_dyn_pf_mismatch_kwd_para ),
        init_pf_ΔPQ, init_dyn_pf_flat_para )

     init_dyn_sh_dxdp =
        -(svd(init_dyn_sd_∂f∂x)\init_dyn_sd_∂f∂p)

    
end



function get_intra_dyn_pf_by_vh_θh_by_parts_sensitivity(
    init_flat_vh_flat_θh,
    # flat_ωref0_vref0_porder0,
    flat_ωs_ωref0_vref0_porder0,
    init_dyn_pf_flat_para ;
    kwd = intra_dyn_pf_by_vh_θh_δ_ed_eq_by_parts_sense_kwd )

    (;
     # sim_timespan,
     sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd,
     
     flat_δ_ed_dash_eq_dash_Idxs_in_state, 
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     
     intra_pf_by_vh_θh_by_parts_kwd,
     intra_dyn_pf_mismatch_kwd_para
     ) = kwd


     (;
      sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para,
      # sd_per_gen_ode_vtf_pf_kwd_para,

      sim_state_x0,
      im_mass_matrix,
      im_sym,     
      ode_alg,
      dt,
      sim_timespan
      ) = sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd

    # flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para = [init_flat_vh_flat_θh; flat_ωref0_vref0_porder0; init_dyn_pf_flat_para ]


    flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para = [init_flat_vh_flat_θh; flat_ωs_ωref0_vref0_porder0; init_dyn_pf_flat_para ]
    
    # 
    #--------------------------------------------------
    #--------------------------------------------------
    # sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    #--------------------------------------------------
    #--------------------------------------------------
    
    sd_per_gen_ode_vtf_pf_kwd_para =
        sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para 
    
    sim_func  =
        sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    
    #-----------------------------------------------------

     sd_dynamics_per_gen_ode_vtf_pf_sim_sol =
         DifferentialEquations.solve(
             ODEProblem(
                 ODEFunction{true}(
                     (dx, x, p, t) ->
                         sim_func(dx, x, p, t;
                                  kwd_para =
                                      sd_per_gen_ode_vtf_pf_kwd_para
                              );
                         mass_matrix =
                             im_mass_matrix ,
                         syms =
                             im_sym ),
                     sim_state_x0 ,
                     sim_timespan,
                     sim_fun_para ),
                 ode_alg,
                 dt=dt )
        
    #---------------------------------------------
    # intra sensitivity
    #---------------------------------------------

    tt = sim_timespan[1]
    sd_x = sd_dynamics_per_gen_ode_vtf_pf_sim_sol( tt )
    sd_dx = similar(sd_x)
    
    flat_δ_ed_dash_eq_dash =
        sd_x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    flat_δ_ed_dash_eq_dash_x0 =
        sim_state_x0[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    (;
        Pg_Idx,
        Qg_Idxs,
        Png_Idxs,
        Qng_Idxs,
        Pgll_Idxs,
        Qgll_Idxs ) =
            Pg_Qg_Png_Qng_Pgll_Qgll_Idxs

    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    if loc_load_exist == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash;
                P_g_loc_load;
                Q_g_loc_load ]                  

    else
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]   

    end

    intra_pf_sol =
        intra_dyn_pf_by_vh_θh_by_parts_func!(
            init_flat_vh_flat_θh,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_pf_by_vh_θh_by_parts_kwd )    

    intra_pf_x = intra_pf_sol.u
   intra_pf_ΔPQ_id_iq = similar( intra_pf_x )


    intra_pf_fun =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch


    intra_dyn_sd_∂f∂x = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq, x ) ->
            intra_pf_fun(
            pf_ΔPQ_id_iq, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_pf_ΔPQ_id_iq,
        intra_pf_x )

    intra_dyn_sd_∂f∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq,  p ) -> intra_pf_fun(
            pf_ΔPQ_id_iq, intra_pf_x, p;
            intra_dyn_pf_mismatch_kwd_para =
                intra_dyn_pf_mismatch_kwd_para ),
        intra_pf_ΔPQ_id_iq,
        intra_dyn_pf_mismatch_flat_para )

       intra_dyn_sd_dxdp =
        -(svd( intra_dyn_sd_∂f∂x)\ intra_dyn_sd_∂f∂p)
    
end
 

function get_intra_dyn_current_balance_pf_by_vh_θh_sensitivity(flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para ;
    kwd = intra_dyn_current_balance_pf_by_vh_θh_sense_kwd )

    (;
     sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd,
     
     flat_δ_ed_dash_eq_dash_Idxs_in_state,
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     
     intra_pf_by_vh_θh_by_parts_kwd,
     intra_dyn_pf_mismatch_kwd_para
     ) = kwd


     (;
      sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para,
      # sd_per_gen_ode_vtf_pf_kwd_para,

      sim_state_x0,
      im_mass_matrix,
      im_sym,     
      ode_alg,
      dt,
      sim_timespan
      ) = sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd

    #--------------------------------------------------

    flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx =
        sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para.flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx

    flat_vh_θh_idx, flat_ωs_ωref0_vref0_porder0_idx, flat_init_dyn_pf_para_idx = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx

   
    init_flat_vh_flat_θh =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[flat_vh_θh_idx]
    
    flat_ωs_ωref0_vref0_porder0 =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_init_dyn_pf_para_idx ]

    
    #--------------------------------------------------
    #--------------------------------------------------
    # sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    #--------------------------------------------------
    #--------------------------------------------------
    
    sd_per_gen_ode_vtf_pf_kwd_para =
        sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para 
    
    sim_func  =
        sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    
    #-----------------------------------------------------

     sd_dynamics_per_gen_ode_vtf_pf_sim_sol =
         DifferentialEquations.solve(
             ODEProblem(
                 ODEFunction{true}(
                     (dx, x, p, t) ->
                         sim_func(dx, x, p, t;
                                  kwd_para =
                                      sd_per_gen_ode_vtf_pf_kwd_para
                              );
                         mass_matrix =
                             im_mass_matrix ,
                         syms =
                             im_sym ),
                     sim_state_x0 ,
                     sim_timespan,
                     sim_fun_para ),
                 ode_alg,
                 dt=dt )
        
    #---------------------------------------------
    # current mismatch sensitivity
    #---------------------------------------------

    tt    = sim_timespan[1]
    sd_x  = sd_dynamics_per_gen_ode_vtf_pf_sim_sol( tt )
    sd_dx = similar(sd_x)
    
    flat_δ_ed_dash_eq_dash =
        sd_x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    flat_δ_ed_dash_eq_dash_x0 =
        sim_state_x0[flat_δ_ed_dash_eq_dash_Idxs_in_state]


    (;
        Pg_Idx,
        Qg_Idxs,
        Png_Idxs,
        Qng_Idxs,
        Pgll_Idxs,
        Qgll_Idxs ) =
            Pg_Qg_Png_Qng_Pgll_Qgll_Idxs

    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    if loc_load_exist == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash;
                P_g_loc_load;
                Q_g_loc_load ]                  

    else
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]   

    end
    
    current_balance_dyn_pf_sol =
        intra_dyn_current_balance_pf_by_vh_θh_by_parts_func!(
            init_flat_vh_flat_θh,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_pf_by_vh_θh_by_parts_kwd )
    

    current_balance_pf_x = current_balance_dyn_pf_sol.u

    current_balance_pf_ΔPQ_id_iq = similar( current_balance_pf_x )
    
    # current_balance_pf_ΔPQ_id_iq = similar( intra_pf_x )
    

    current_balance_pf_fun =
        get_a_model_integrated_intra_dyn_current_balance_pf_ΔI_idq_mismatch


    current_balance_dyn_sd_∂f∂x = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq, x ) ->
            current_balance_pf_fun(
                pf_ΔPQ_id_iq,
                x,
                intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        current_balance_pf_ΔPQ_id_iq,
        current_balance_pf_x )

    
    current_balance_dyn_sd_∂f∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq,  p ) -> current_balance_pf_fun(
            pf_ΔPQ_id_iq,
            current_balance_pf_x,
            p;
            intra_dyn_pf_mismatch_kwd_para =
                intra_dyn_pf_mismatch_kwd_para ),
        current_balance_pf_ΔPQ_id_iq,
        intra_dyn_pf_mismatch_flat_para )

       current_balance_dyn_sd_dxdp =
           -(svd(current_balance_dyn_sd_∂f∂x)\
           current_balance_dyn_sd_∂f∂p)
    
end



function get_intra_dyn_current_balance_pf_by_vh_θh_by_parts_sensitivity(
    init_flat_vh_flat_θh,
    flat_ωs_ωref0_vref0_porder0,
    init_dyn_pf_flat_para ;
    kwd = intra_dyn_current_balance_pf_by_vh_θh_sense_kwd )

    (;
     sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd,
     
     flat_δ_ed_dash_eq_dash_Idxs_in_state,
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     
     intra_pf_by_vh_θh_by_parts_kwd,
     intra_dyn_pf_mismatch_kwd_para
     ) = kwd


     (;
      sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para,
      # sd_per_gen_ode_vtf_pf_kwd_para,

      sim_state_x0,
      im_mass_matrix,
      im_sym,     
      ode_alg,
      dt,
      sim_timespan
      ) = sd_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_kwd

    flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para = [ init_flat_vh_flat_θh; flat_ωs_ωref0_vref0_porder0; init_dyn_pf_flat_para ]
    
    #--------------------------------------------------
    #--------------------------------------------------
    # sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    #--------------------------------------------------
    #--------------------------------------------------
    
    sd_per_gen_ode_vtf_pf_kwd_para =
        sd_per_gen_ode_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para 
    
    sim_func  =
        sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!
    
    #-----------------------------------------------------

     sd_dynamics_per_gen_ode_vtf_pf_sim_sol =
         DifferentialEquations.solve(
             ODEProblem(
                 ODEFunction{true}(
                     (dx, x, p, t) ->
                         sim_func(dx, x, p, t;
                                  kwd_para =
                                      sd_per_gen_ode_vtf_pf_kwd_para
                              );
                         mass_matrix =
                             im_mass_matrix ,
                         syms =
                             im_sym ),
                     sim_state_x0 ,
                     sim_timespan,
                     sim_fun_para ),
                 ode_alg,
                 dt=dt )
        
    #---------------------------------------------
    # current mismatch sensitivity
    #---------------------------------------------

    tt    = sim_timespan[1]
    sd_x  = sd_dynamics_per_gen_ode_vtf_pf_sim_sol( tt )
    sd_dx = similar(sd_x)
    
    flat_δ_ed_dash_eq_dash =
        sd_x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    flat_δ_ed_dash_eq_dash_x0 =
        sim_state_x0[flat_δ_ed_dash_eq_dash_Idxs_in_state]


    (;
        Pg_Idx,
        Qg_Idxs,
        Png_Idxs,
        Qng_Idxs,
        Pgll_Idxs,
        Qgll_Idxs ) =
            Pg_Qg_Png_Qng_Pgll_Qgll_Idxs

    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    if loc_load_exist == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash;
                P_g_loc_load;
                Q_g_loc_load ]                  

    else
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]   

    end
    
    current_balance_dyn_pf_sol =
        intra_dyn_current_balance_pf_by_vh_θh_by_parts_func!(
            init_flat_vh_flat_θh,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_pf_by_vh_θh_by_parts_kwd )
    

    current_balance_pf_x = current_balance_dyn_pf_sol.u
    
    current_balance_pf_ΔPQ_id_iq = similar( current_balance_pf_x )
    
   # current_balance_pf_ΔPQ_id_iq = similar( intra_pf_x )
    

    current_balance_pf_fun =
        get_a_model_integrated_intra_dyn_current_balance_pf_ΔI_idq_mismatch


    current_balance_dyn_sd_∂f∂x = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq, x ) ->
            current_balance_pf_fun(
                pf_ΔPQ_id_iq,
                x,
                intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        current_balance_pf_ΔPQ_id_iq,
        current_balance_pf_x )

    
    current_balance_dyn_sd_∂f∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq,  p ) -> current_balance_pf_fun(
            pf_ΔPQ_id_iq,
            current_balance_pf_x,
            p;
            intra_dyn_pf_mismatch_kwd_para =
                intra_dyn_pf_mismatch_kwd_para ),
        current_balance_pf_ΔPQ_id_iq,
        intra_dyn_pf_mismatch_flat_para )

       current_balance_dyn_sd_dxdp =
           -(svd(current_balance_dyn_sd_∂f∂x)\
           current_balance_dyn_sd_∂f∂p)
    
end



function get_small_signal_stability_sd_dynamics_wt_intra_dyn_pf(flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para ;kwd = stability_sd_dynamics_wt_intra_dyn_pf_kwd )

    (;
     loc_load_exist,
     
     sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para,

     # flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,

     # flat_ωs_ωref0_vref0_porder0,
     
     sim_state_x0,
     # flat_vh_flat_θh_x0,
     # init_dyn_pf_flat_para,
     gens_nodes_ra_Xd_dash_Xq_dash,

     im_mass_matrix,
     im_sym, 
     ode_alg,
     dt,

     flat_vh_flat_θh_Idx,

     gens_nodes_idx,
     non_gens_nodes_idx,

     flat_δ_ed_dash_eq_dash_Idxs_in_state,
     δ_ed_dash_eq_dash_Idxs_in_flattend,

     intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd,
     intra_dyn_pf_mismatch_kwd_para,

     im_state_idx_in_∂f∂x,
     net_vh_θh_idx_in_∂f∂x,

     consecutive_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
     vec_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,     

     #
     vec_flat_vh_flat_θh_idx_in_∂g∂x,
     vec_flat_vh_θh_idx_in_∂g∂x,
     vec_id_iq_idx_in_∂g∂x,
     
     consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
     consecutive_flat_vh_θh_idx_in_∂g∂x,
     #
     
     gens_consecutive_flat_vh_θh_idx_in_∂g∂x,
     non_gens_consecutive_flat_vh_θh_idx_in_∂g∂x,

     gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
     non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
     consecutive_id_iq_idx_in_∂g∂x,

     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,

     consecutive_Pg_idx_Qg_idx_in_mismatch,
     consecutive_Png_idx_Qng_idx_in_mismatch,      
     consecutive_id_idx_iq_idx_in_mismatch,

     consecutive_dyn_gens_PQ_idx_in_Idxs,
     consecutive_dyn_non_gens_PQ_idx_in_Idxs,
     consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs,
     dyn_δ_ed_eq_pf_Idxs,   

     # ∂f∂x_xp_by_vh_θh_∂f∂x,

     im_pure_states_idxs,
     im_τm_tilade_vf_states_idxs,
     ) = kwd 

    #--------------------------------------------------

    flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx = sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para.flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx

    flat_vh_θh_idx, flat_ωs_ωref0_vref0_porder0_idx, flat_init_dyn_pf_para_idx = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx

   
    init_flat_vh_flat_θh =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[flat_vh_θh_idx]
    
    flat_ωs_ωref0_vref0_porder0 =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_init_dyn_pf_para_idx ]

    #--------------------------------------------------

    
    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    
    #--------------------------------------------------
    #--------------------------------------------------
    # sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!
    #--------------------------------------------------
    #--------------------------------------------------
    
    sim_fun_kwd_para =
        sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para

    # flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para

    sim_func  =
        sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!
    
    #-----------------------------------------------------

    sd_per_gen_ode_sep_vtf_pf_sim_sol =
        DifferentialEquations.solve(
            ODEProblem(
                ODEFunction{true}(
                    (dx, x, p, t) ->
                        sim_func(dx, x, p, t;
                                 kwd_para =
                                     sim_fun_kwd_para
                         );
                    mass_matrix =
                        im_mass_matrix ,
                    syms =
                        im_sym ),
                sim_state_x0 ,
                sim_timespan,
                sim_fun_para ),
            ode_alg,
            dt=dt )        
    

    #--------------------------------------------------
    #--------------------------------------------------
    # Method 1
    #--------------------------------------------------
    #--------------------------------------------------

    #---------------------------------------------
    #---------------------------------------------
    # sensitivity
    # sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!
    #---------------------------------------------
    #---------------------------------------------

    tt   = sim_timespan[1]
    sd_x = sd_per_gen_ode_sep_vtf_pf_sim_sol( tt )
           
   sd_dx = similar(sd_x)
    
    sim_fun_kwd_para =
        sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para

    sim_func  =
        sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!

    #---------------------------------------------
        
    sd_per_gen_ode_sep_vtf_pf_∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> sim_func(
            dx, x, sim_fun_para, tt;
            kwd_para = sim_fun_kwd_para ),
        sd_dx, sd_x )
    
    sd_per_gen_ode_sep_vtf_pf_∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> sim_func(
            dx, sd_x, p, tt;
            kwd_para = sim_fun_kwd_para ),
        sd_dx, sim_fun_para )

    sd_per_gen_ode_sep_vtf_pf_dxdp =
        -( svd( sd_per_gen_ode_sep_vtf_pf_∂f∂x )\
        sd_per_gen_ode_sep_vtf_pf_∂f∂p )


    #---------------------------------------------
    #---------------------------------------------
    # sensitivity
    # intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq
    #---------------------------------------------
    #---------------------------------------------


    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------

    #----------------------------------------    

    (flat_vh_idx_in_flat_vh_flat_θh,
     flat_θh_idx_in_flat_vh_flat_θh) =
         flat_vh_flat_θh_Idx
    
    flat_vh = init_flat_vh_flat_θh[
        flat_vh_idx_in_flat_vh_flat_θh ]

    flat_θh = init_flat_vh_flat_θh[
        flat_θh_idx_in_flat_vh_flat_θh ]
    
    gens_vh = flat_vh[ gens_nodes_idx ]

    gens_θh = flat_θh[ gens_nodes_idx ]        

    #----------------------------------------    

    flat_δ_ed_dash_eq_dash =
        sd_x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    flat_δ_ed_dash_eq_dash_x0 =
        sim_state_x0[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    #----------------------------------------    
    #----------------------------------------    

    gens_δ_ed_dash_eq_dash = [
        flat_δ_ed_dash_eq_dash[idx]
        for idx in
            δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    
    gens_δ =
        first.( gens_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_δ_ed_dash_eq_dash )

    #----------------------------------------

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh, a_θh, a_δ,
            a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
        for (a_vh, a_θh,
             a_δ, a_ed_dash, a_eq_dash,
             a_ra, a_X_d_dash, a_X_q_dash ) in
            zip( gens_vh, gens_θh,
                 gens_δ, gens_ed_dash, gens_eq_dash,
                 gens_ra, gens_Xd_dash, gens_Xq_dash ) ]

    gens_i_d =
        first.( gens_id_iq )
    
    gens_i_q =
        last.( gens_id_iq )

    #----------------------------------------    
    
    init_flat_vh_flat_θh_cal_id_iq =
        [flat_vh; flat_θh;
         gens_i_d; gens_i_q  ]

    #----------------------------------------   

    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    if loc_load_exist == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens; Q_gens; P_non_gens;
                Q_non_gens; flat_δ_ed_dash_eq_dash;
                P_g_loc_load; Q_g_loc_load ]
    else
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens; Q_gens; P_non_gens;
                Q_non_gens; flat_δ_ed_dash_eq_dash ]   

    end
    
    #----------------------------------------
    #----------------------------------------

    intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_sol =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_func!(
            init_flat_vh_flat_θh_cal_id_iq,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            kwd_para =
                intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd )

    #----------------------------------------

    intra_dyn_id_iq_δ_pf_x =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_sol.u

    intra_dyn_id_iq_δ_pf_ΔPQ_id_iq =
        similar( intra_dyn_id_iq_δ_pf_x )
    

    intra_dyn_id_iq_δ_pf_fun =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch


    intra_dyn_id_iq_δ_pf_∂g∂x = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq, x ) ->
            intra_dyn_id_iq_δ_pf_fun( pf_ΔPQ_id_iq, x,
                intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_dyn_id_iq_δ_pf_ΔPQ_id_iq,
        intra_dyn_id_iq_δ_pf_x  )

    
    intra_dyn_id_iq_δ_pf_∂g∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq,  p ) ->
            intra_dyn_id_iq_δ_pf_fun(
                pf_ΔPQ_id_iq,
                intra_dyn_id_iq_δ_pf_x,
                p;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_dyn_id_iq_δ_pf_ΔPQ_id_iq,
        intra_dyn_pf_mismatch_flat_para )

       intra_dyn_id_iq_δ_pf_∂x∂p =
           -(svd( intra_dyn_id_iq_δ_pf_∂g∂x )\
            intra_dyn_id_iq_δ_pf_∂g∂p)

    #---------------------------------------------
    # sd_per_gen_ode_sep_vtf_pf_∂f∂x
    #---------------------------------------------

    # xp : pure states

    # xa : algebraic states

    ∂f∂x_xp_by_xp_∂f∂x = Ax_states_by_states_∂f∂x =
        sd_per_gen_ode_sep_vtf_pf_∂f∂x[
            im_state_idx_in_∂f∂x,
            im_state_idx_in_∂f∂x ]

    ∂f∂x_xp_by_vh_θh_∂f∂x = Ax_states_by_vh_θh_∂f∂x = 
        sd_per_gen_ode_sep_vtf_pf_∂f∂x[
            im_state_idx_in_∂f∂x,
            net_vh_θh_idx_in_∂f∂x ]

    ∂f∂x_vh_θh_by_xp_∂f∂x = Ax_vh_θh_by_states_∂f∂x =
        sd_per_gen_ode_sep_vtf_pf_∂f∂x[
            net_vh_θh_idx_in_∂f∂x,
            im_state_idx_in_∂f∂x]

    ∂f∂x_vh_θh_by_vh_θh_∂f∂x = Ax_vh_θh_by_vh_θh_∂f∂x = 
        sd_per_gen_ode_sep_vtf_pf_∂f∂x[
            net_vh_θh_idx_in_∂f∂x,
            net_vh_θh_idx_in_∂f∂x]

    #---------------------------------------------
    # sd_per_gen_ode_sep_vtf_pf_∂f∂p
    #---------------------------------------------

    ∂f∂p_xp_by_ωs_ωref0_vref0_porder0_∂f∂p =
        Ax_states_by_ωs_ωref0_vref0_porder0_∂f∂p =
        sd_per_gen_ode_sep_vtf_pf_∂f∂p[
            im_state_idx_in_∂f∂x,
            :]

    ∂f∂p_xa_by_ωs_ωref0_vref0_porder0_∂f∂p =
        Ax_vh_θh_by_ωs_ωref0_vref0_porder0_∂f∂p =
            sd_per_gen_ode_sep_vtf_pf_∂f∂p[
                net_vh_θh_idx_in_∂f∂x,
                :]

    #---------------------------------------------
    # intra_dyn_id_iq_δ_pf_∂g∂x
    #---------------------------------------------

    ∂g∂x_id_iq_by_id_iq =
        Mx_gens_id_iq_δ_by_gens_id_iq_δ_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            consecutive_id_iq_idx_in_∂g∂x ]

    ∂g∂x_id_iq_by_vθ_g =
        Mx_gens_id_iq_δ_by_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x]

    ∂g∂x_id_iq_by_vθ_ng =
        Mx_gens_id_iq_δ_by_non_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    #---------------------------------------------

    ∂g∂x_g_by_id_iq =
        Mx_gens_PQ_by_gens_id_iq_δ_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch ,
            consecutive_id_iq_idx_in_∂g∂x ]

    ∂g∂x_g_by_vθ_g =
        Mx_gens_PQ_by_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch ,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]


    ∂g∂x_g_by_vθ_ng =
        Mx_gens_PQ_by_non_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch ,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]


    #---------------------------------------------

    ∂g∂x_ng_by_id_iq =
        Mx_non_gens_PQ_by_gens_id_iq_δ_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch ,
            consecutive_id_iq_idx_in_∂g∂x ]
    
    ∂g∂x_ng_by_vθ_g =
        Mx_non_gens_PQ_by_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch ,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    ∂g∂x_ng_by_vθ_ng =
        Mx_non_gens_PQ_by_non_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch ,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    #---------------------------------------------
    # intra_dyn_id_iq_δ_pf_∂g∂p
    #---------------------------------------------

    ∂g∂p_id_iq_by_gens_PQ =
        Mx_gens_id_iq_δ_by_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_id_idx_iq_idx_in_mismatch ,
             consecutive_dyn_gens_PQ_idx_in_Idxs ]

    ∂g∂p_id_iq_by_non_gens_PQ =
        Mx_gens_id_iq_δ_by_non_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_id_idx_iq_idx_in_mismatch,
            consecutive_dyn_non_gens_PQ_idx_in_Idxs ]


    ∂g∂p_id_iq_by_δ_ed_eq =
        Mx_gens_id_iq_δ_by_dyn_δ_ed_eq_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_id_idx_iq_idx_in_mismatch,
            dyn_δ_ed_eq_pf_Idxs ]


    if loc_load_exist == true
        
        ∂g∂p_id_iq_by_PQll_g =
            Mx_gens_id_iq_δ_by_gens_loc_load_PQ_∂g∂p =
            intra_dyn_id_iq_δ_pf_∂g∂p[
                consecutive_id_idx_iq_idx_in_mismatch,
                consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs ]

    end

    #---------------------------------------------

    ∂g∂p_g_by_PQ_g = Mx_gens_PQ_by_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            consecutive_dyn_gens_PQ_idx_in_Idxs ]

    ∂g∂p_g_by_PQ_ng = Mx_gens_PQ_by_non_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            consecutive_dyn_non_gens_PQ_idx_in_Idxs ]

    ∂g∂p_g_by_δ_ed_eq = Mx_gens_PQ_by_dyn_δ_ed_eq_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            dyn_δ_ed_eq_pf_Idxs ]


    if loc_load_exist == true
        
        ∂g∂p_g_by_PQll_g =
            Mx_gens_PQ_by_gens_loc_load_PQ_∂g∂p =
            intra_dyn_id_iq_δ_pf_∂g∂p[
                consecutive_Pg_idx_Qg_idx_in_mismatch,
                consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs ]

    end

    #---------------------------------------------

    ∂g∂p_ng_by_PQ_g = Mx_non_gens_PQ_by_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            consecutive_dyn_gens_PQ_idx_in_Idxs ]

    ∂g∂p_ng_by_PQ_ng = Mx_non_gens_PQ_by_non_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            consecutive_dyn_non_gens_PQ_idx_in_Idxs ]

    ∂g∂p_ng_by_δ_ed_eq = Mx_non_gens_PQ_by_dyn_δ_ed_eq_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            dyn_δ_ed_eq_pf_Idxs ]

    if loc_load_exist == true
        
         ∂g∂p_ng_by_PQll_g =
             Mx_non_gens_PQ_by_gens_loc_load_PQ_∂g∂p =
            intra_dyn_id_iq_δ_pf_∂g∂p[
                consecutive_Png_idx_Qng_idx_in_mismatch,
                consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs ]

    end

    #---------------------------------------------

    # Define Q_vθ_g as:
    # ∂f∂x_g_by_vθ_g
    Q_vθ_g = ( ∂g∂x_g_by_vθ_g -
        ∂g∂x_g_by_vθ_ng *
        ( ∂g∂x_ng_by_vθ_ng\∂g∂x_ng_by_vθ_g) -
            ∂g∂x_g_by_id_iq *
            (∂g∂x_id_iq_by_id_iq\∂g∂x_id_iq_by_vθ_g)) 

    # Define Q_δ_ed_eq as:
    
    Q_δ_ed_eq =
        (∂g∂p_g_by_δ_ed_eq -
        ∂g∂x_g_by_id_iq *
        (∂g∂x_id_iq_by_id_iq\
        ∂g∂p_id_iq_by_δ_ed_eq ))

    Asys =  copy(∂f∂x_xp_by_xp_∂f∂x)

    Asys[:, flat_δ_ed_dash_eq_dash_Idxs_in_state ] .=
        Asys[:, flat_δ_ed_dash_eq_dash_Idxs_in_state] -
        ∂f∂x_xp_by_vh_θh_∂f∂x[
            :, gens_consecutive_flat_vh_θh_idx_in_∂g∂x] *
                ( Q_vθ_g \ Q_δ_ed_eq ) 

    #---------------------------------------------

    Axp =
        Asys[im_pure_states_idxs,
             im_pure_states_idxs ]

    Axa =
        Asys[im_pure_states_idxs,
             im_τm_tilade_vf_states_idxs ]

    Aap =
        Asys[im_τm_tilade_vf_states_idxs,
             im_pure_states_idxs ]

    Aaa =
        Asys[im_τm_tilade_vf_states_idxs,
             im_τm_tilade_vf_states_idxs ]

    
    pure_Asys = Axp - Axa * (Aaa \ Aap )

    #---------------------------------------------

    eigen_object = eigen( pure_Asys )

    eig_values = eigen_object.values

    printed_eig_values =
        round.(eig_values; digits=4)

    eigvecs_right = eigen_object.vectors

    # inv_eigvecs_right = inv(eigvecs_right)
    
    inv_eigvecs_right =
        eigvecs_right \ I( size( eigvecs_right )[1] )

    #--------------------------------------------------
    # states associated with eig_values
    #--------------------------------------------------

    M_diag = inv_eigvecs_right * pure_Asys * eigvecs_right

    printed_M_diag =
        round.( M_diag; digits = 4 )

    # tup_state_var_eig_value =
    #     [(a_state_sym, an_eig_value)
    #      for ( a_state_sym, an_eig_value )  in
    #          zip( im_ode_sym, printed_eig_values )]

    #--------------------------------------------------
    # participation factor
    #--------------------------------------------------

    PF_pure_Asys  = get_participation_factors( pure_Asys )

    printed_PF_Asys =
        round.( PF_pure_Asys; digits = 4 )

    #--------------------------------------------------

    eig_values, eigvecs_left, eigvecs_right =
        get_eigens(pure_Asys)


    #--------------------------------------------------

    im_pure_states_syms =
        im_sym[im_pure_states_idxs]
    
    #--------------------------------------------------
    #--------------------------------------------------
    return (; pure_Asys, eig_values,
            eigvecs_left, eigvecs_right,
            inv_eigvecs_right,
            M_diag, PF_pure_Asys,
            im_pure_states_syms )
    
end



function get_small_signal_stability_analysis_ode_per_gen_wt_intra_dyn_pf(init_flat_vh_flat_θh, init_dyn_pf_flat_para, flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh ; kwd = stability_ode_per_gen_wt_intra_dyn_pf_kwd )

    (;
     loc_load_exist,
     
     sim_state_x0,
     # flat_vh_flat_θh_x0,

     im_vars_Idx_in_state,

     ode_only_per_gen_id_iq_pg_vh_kwd,

     # flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh,

     # flat_ωs_ωref0_vref0_porder0,

     im_ode_mass_matrix,
     im_ode_sym,

     gens_nodes_ra_Xd_dash_Xq_dash,

     flat_vh_flat_θh_Idx,

     gens_nodes_idx,
     non_gens_nodes_idx,

     flat_δ_ed_dash_eq_dash_Idxs_in_state,
     δ_ed_dash_eq_dash_Idxs_in_flattend,

     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,

     # init_dyn_pf_flat_para,
     intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd,
     intra_dyn_pf_mismatch_kwd_para,

     im_state_idx_in_∂f∂x,

     flat_id_iq_pg_vh_idx_in_∂f∂p,
     flat_ωs_ωref0_vref0_porder0_idx_in_∂f∂p,
     gens_id_iq_idx_in_∂f∂p,
     gens_vh_idx_in_∂f∂p,

     gens_ph_idx_in_∂f∂p,

     consecutive_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
     vec_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
     

     #
     vec_flat_vh_flat_θh_idx_in_∂g∂x,
     vec_flat_vh_θh_idx_in_∂g∂x,
     vec_id_iq_idx_in_∂g∂x,

     consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
     consecutive_flat_vh_θh_idx_in_∂g∂x,
     consecutive_id_iq_idx_in_∂g∂x,
     #         
     
     consecutive_gens_vh_idx_θh_idx_in_Idx,     

     gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,     
     non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,

     
     intra_flat_δ_ed_dash_eq_dash_idx_in_∂g∂p,

     consecutive_Pg_idx_Qg_idx_in_mismatch,
     consecutive_Png_idx_Qng_idx_in_mismatch,
     consecutive_id_idx_iq_idx_in_mismatch,

     im_pure_states_idxs,
     im_τm_tilade_vf_states_idxs
     ) = kwd

    #-----------------------------------------------------

    flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_only_per_gen_id_iq_pg_vh_kwd.flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx


    flat_ωs_ωref0_vref0_porder0_Idx, flat_id_iq_pg_vh_Idx =
        flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx

        
     flat_ωs_ωref0_vref0_porder0 =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh[
            flat_ωs_ωref0_vref0_porder0_Idx ]
    
    flat_id_iq_pg_vh =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh[
            flat_id_iq_pg_vh_Idx ]
    
    
    #-----------------------------------------------------

    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs

    (;
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
     ) = dyn_pf_P_Q_δ_etc_kwd_para_Idxs 

    
    gens_PQ_idx_in_Idxs =
        [ dyn_P_gens_Idxs,
          dyn_Q_gens_Idxs ]

    
    ( consecutive_dyn_gens_PQ_idx_in_Idxs,
      vec_dyn_gens_PQ_idx_in_Idxs ) =
            convert_to_consecutive_idxs(
                gens_PQ_idx_in_Idxs )

    consecutive_Pg_Qg_idx_in_∂g∂p =
        consecutive_dyn_gens_PQ_idx_in_Idxs
    
    non_gens_PQ_idx_in_Idxs =
        [ dyn_P_non_gens_Idxs,
          dyn_Q_non_gens_Idxs ]

    ( consecutive_dyn_non_gens_PQ_idx_in_Idxs,
      vec_dyn_non_gens_PQ_idx_in_Idxs ) =
            convert_to_consecutive_idxs(
                non_gens_PQ_idx_in_Idxs )

    consecutive_Png_Qng_idx_in_∂g∂p =
        consecutive_dyn_non_gens_PQ_idx_in_Idxs
    
    if loc_load_exist == true
           
        gens_loc_load_PQ_idx_in_Idxs =
            [ dyn_P_g_loc_load_Idxs,
              dyn_Q_g_loc_load_Idxs ]

        ( consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs,
          vec_dyn_gens_loc_load_PQ_idx_in_Idxs ) =
                convert_to_consecutive_idxs(
                    gens_loc_load_PQ_idx_in_Idxs )
        
        consecutive_Pgll_Qgll_idx_in_∂g∂p =
            consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs

    end
    
    #--------------------------------------------------
    # ode_only_per_gen_model_by_id_iq_pg_vh_func!
    #--------------------------------------------------

    ode_only_sim_state_x0 =
        sim_state_x0[ im_vars_Idx_in_state ]

    sim_fun_kwd_para =
        ode_only_per_gen_id_iq_pg_vh_kwd
    
    sim_fun_para =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh 
    
    sim_func  =
        ode_only_per_gen_model_by_id_iq_pg_vh_func!
    
    #-----------------------------------------------------

    ode_only_per_gen_by_id_iq_pg_vh_sim_sol =
        DifferentialEquations.solve(
        ODEProblem(
            ODEFunction{true}(
                (dx, x, p, t) ->
                    sim_func(dx, x, p, t;
                             kwd_para =
                                 sim_fun_kwd_para 
                         );
                mass_matrix =
                    # im_mass_matrix,
                    im_ode_mass_matrix,
                syms =
                    # im_sym
                    im_ode_sym ),

            # sim_state_x0,            
            # ode_sim_state_x0,
            ode_only_sim_state_x0,
            sim_timespan,
            sim_fun_para ),
            ode_alg,
            dt = dt)
    
    #---------------------------------------------
    # sensitivity
    #---------------------------------------------

    tt   = sim_timespan[1]
    sd_x = ode_only_per_gen_by_id_iq_pg_vh_sim_sol( tt )
 
    sd_dx = similar(sd_x)
    
    sim_fun_kwd_para =
        ode_only_per_gen_id_iq_pg_vh_kwd
            
    sim_fun_para =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh 
        
    sim_func  =
        ode_only_per_gen_model_by_id_iq_pg_vh_func!
        
    #---------------------------------------------

    ode_only_per_gen_by_id_iq_pg_vh_∂f∂x =
        ForwardDiff.jacobian(
        (dx, x ) -> sim_func(
            dx, x, sim_fun_para, tt;
            kwd_para = sim_fun_kwd_para ),
        sd_dx, sd_x )


    ode_only_per_gen_by_id_iq_pg_vh_∂f∂p =
        ForwardDiff.jacobian(
        (dx, p ) -> sim_func(
            dx, sd_x, p, tt;
            kwd_para = sim_fun_kwd_para ),
        sd_dx, sim_fun_para )

    #---------------------------------------------
    # intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq sensitivity
    #---------------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------
    #----------------------------------------    

    (flat_vh_idx_in_flat_vh_flat_θh,
     flat_θh_idx_in_flat_vh_flat_θh) =
         flat_vh_flat_θh_Idx
    
    flat_vh = init_flat_vh_flat_θh[
        flat_vh_idx_in_flat_vh_flat_θh ]

    flat_θh = init_flat_vh_flat_θh[
        flat_θh_idx_in_flat_vh_flat_θh ]
    
    gens_vh = flat_vh[ gens_nodes_idx ]

    gens_θh = flat_θh[ gens_nodes_idx ]        

    #----------------------------------------    

    flat_δ_ed_dash_eq_dash =
        sd_x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    flat_δ_ed_dash_eq_dash_x0 =
        sim_state_x0[flat_δ_ed_dash_eq_dash_Idxs_in_state]

    #----------------------------------------    
    #----------------------------------------    

    gens_δ_ed_dash_eq_dash = [
        flat_δ_ed_dash_eq_dash[idx]
        for idx in
            δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    
    gens_δ =
        first.( gens_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_δ_ed_dash_eq_dash )

    #----------------------------------------

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh, a_θh, a_δ,
            a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
        for (a_vh, a_θh,
             a_δ, a_ed_dash, a_eq_dash,
             a_ra, a_X_d_dash, a_X_q_dash ) in
            zip( gens_vh, gens_θh,
                 gens_δ, gens_ed_dash, gens_eq_dash,
                 gens_ra, gens_Xd_dash, gens_Xq_dash ) ]

    gens_i_d =
        first.( gens_id_iq )
    
    gens_i_q =
        last.( gens_id_iq )

    #----------------------------------------    
    
    init_flat_vh_flat_θh_cal_id_iq =
        [flat_vh; flat_θh;
         gens_i_d; gens_i_q  ]

    #----------------------------------------   

    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    if loc_load_exist == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens; Q_gens; P_non_gens;
                Q_non_gens; flat_δ_ed_dash_eq_dash;
                P_g_loc_load; Q_g_loc_load ]
        

    else
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens; Q_gens; P_non_gens;
                Q_non_gens; flat_δ_ed_dash_eq_dash ]   

    end
    
    #----------------------------------------
    #----------------------------------------

    intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_sol =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_func!(
            init_flat_vh_flat_θh_cal_id_iq,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            kwd_para =
                intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd)

    #----------------------------------------

    intra_dyn_id_iq_δ_pf_x =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_sol.u

    intra_dyn_id_iq_δ_pf_ΔPQ_id_iq =
        similar( intra_dyn_id_iq_δ_pf_x )
    

    intra_dyn_id_iq_δ_pf_fun =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch


    intra_dyn_id_iq_δ_pf_∂g∂x = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq, x ) ->
            intra_dyn_id_iq_δ_pf_fun( pf_ΔPQ_id_iq, x,
                intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_dyn_id_iq_δ_pf_ΔPQ_id_iq,
        intra_dyn_id_iq_δ_pf_x  )

    
    intra_dyn_id_iq_δ_pf_∂g∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq,  p ) ->
            intra_dyn_id_iq_δ_pf_fun(
                pf_ΔPQ_id_iq,
                intra_dyn_id_iq_δ_pf_x,
                p;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_dyn_id_iq_δ_pf_ΔPQ_id_iq,
        intra_dyn_pf_mismatch_flat_para )

       intra_dyn_id_iq_δ_pf_∂x∂p =
           -(svd( intra_dyn_id_iq_δ_pf_∂g∂x )\
            intra_dyn_id_iq_δ_pf_∂g∂p)

    #----------------------------------------
    # ode f
    #----------------------------------------

    # A1
    ∂f∂x_im_state_by_im_state =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂x[
            im_state_idx_in_∂f∂x,
            im_state_idx_in_∂f∂x ]
    

    # B1_gens_id_iq
    ∂f∂p_im_state_by_gens_id_iq =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            gens_id_iq_idx_in_∂f∂p ]

    # F
    ∂f∂p_pg = ∂f∂p_im_state_by_gens_ph =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            gens_ph_idx_in_∂f∂p ]

    ∂f∂p_pg_extended_wt_zero_qh =
        hcat( ∂f∂p_pg, zeros( size(∂f∂p_pg )) )

    ∂f∂p_pg_qh0 =
        ∂f∂p_pg_extended_wt_zero_qh[
            :,
            consecutive_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0 ]
    
    #---------------------------------------------

    # B2_gens_vh 
    ∂f∂p_im_state_by_gens_vh =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            gens_vh_idx_in_∂f∂p ]

    # B2_gens_vh_extended_wt_zero_θh
    ∂f∂p_im_state_by_gens_vh_extended_wt_zero_θh =
        hcat(∂f∂p_im_state_by_gens_vh,
             zeros(size(∂f∂p_im_state_by_gens_vh)) )

    # B2
    # B2_gens_vh_wt_zero_θh
    ∂f∂p_im_state_by_gens_vh_wt_zero_θh =
        ∂f∂p_im_state_by_gens_vh_extended_wt_zero_θh[
            :, consecutive_gens_vh_idx_θh_idx_in_Idx ]

    #---------------------------------------------

    # E_∂f∂p_ωs_ωref0_vref0_porder0
    ∂f∂p_im_state_by_ωs_ωref0_vref0_porder0 =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            flat_ωs_ωref0_vref0_porder0_idx_in_∂f∂p ]

    ∂f∂p_ωs_ωref0_vref0_porder0 =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            flat_ωs_ωref0_vref0_porder0_idx_in_∂f∂p ]

    #----------------------------------------
    # pf ∂g∂x
    #----------------------------------------

    # D4
    ∂g∂x_Pg_Qg_mismatch_by_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    # D5
    ∂g∂x_Pg_Qg_mismatch_by_non_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    # D3
    ∂g∂x_Pg_Qg_mismatch_by_id_iq =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            consecutive_id_iq_idx_in_∂g∂x ]

    #---------------------------------------------

    # D6
    ∂g∂x_Png_Qng_mismatch_by_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    # D7
    ∂g∂x_Png_Qng_mismatch_by_non_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    ∂g∂x_Png_Qng_mismatch_by_id_iq =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            consecutive_id_iq_idx_in_∂g∂x ]

    #---------------------------------------------

    # D2
    ∂g∂x_id_iq_mismatch_by_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    ∂g∂x_id_iq_mismatch_by_non_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    # D1
    ∂g∂x_id_iq_mismatch_by_id_iq =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            consecutive_id_iq_idx_in_∂g∂x ]

    #----------------------------------------
    # pf ∂g∂p
    #----------------------------------------

    # C2
    ∂g∂p_Pg_Qg_mismatch_by_δ_ed_dash_eq_dash =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            intra_flat_δ_ed_dash_eq_dash_idx_in_∂g∂p ]

    # C1
    ∂g∂p_id_iq_mismatch_by_δ_ed_dash_eq_dash =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_id_idx_iq_idx_in_mismatch,
            intra_flat_δ_ed_dash_eq_dash_idx_in_∂g∂p ]

    ∂g∂p_Png_Qng_mismatch_by_Png_Qng =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Png_idx_Qng_idx_in_mismatch,
             consecutive_Png_Qng_idx_in_∂g∂p ]


    #---------------------------------------------

    #
    # id_iq_mis_kg_id_iq_by_vθ_g =
    #     -(∂g∂x_id_iq \ ∂g∂x_vθ_g)

    id_iq_mis_kg_id_iq_by_vθ_g =
        -(∂g∂x_id_iq_mismatch_by_id_iq \
        ∂g∂x_id_iq_mismatch_by_gens_vh_θh)


    # id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash =
    #     -(∂g∂x_id_iq \
    #         ∂g∂p_δ_ed_dash_eq_dash)


    id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash =
        -(∂g∂x_id_iq_mismatch_by_id_iq \
            ∂g∂p_id_iq_mismatch_by_δ_ed_dash_eq_dash)

    #
    # Png_Qng_mis_k_vθ_ng_by_vθ_g =
    #     -(∂g∂x_vθ_ng \∂g∂x_vθ_g)

    Png_Qng_mis_k_vθ_ng_by_vθ_g =
        -(∂g∂x_Png_Qng_mismatch_by_non_gens_vh_θh \
        ∂g∂x_Png_Qng_mismatch_by_gens_vh_θh)


    # Png_Qng_mis_k_vθ_ng_by_Png_Qng =
    #     -(∂g∂x_vθ_ng \ ∂g∂p_Png_Qng)

    Png_Qng_mis_k_vθ_ng_by_Png_Qng =
        -(∂g∂x_Png_Qng_mismatch_by_non_gens_vh_θh \
        ∂g∂p_Png_Qng_mismatch_by_Png_Qng )

    #
    # Pg_Qg_mis_k_vθ_g =
    #     (∂g∂x_vθ_g + ∂g∂x_vθ_ng *
    #     Png_Qng_mis_k_vθ_ng_by_vθ_g + ∂g∂x_id_iq *
    #     id_iq_mis_kg_id_iq_by_vθ_g)

    Pg_Qg_mis_k_vθ_g =
        (∂g∂x_Pg_Qg_mismatch_by_gens_vh_θh +
        ∂g∂x_Pg_Qg_mismatch_by_non_gens_vh_θh *
        Png_Qng_mis_k_vθ_ng_by_vθ_g +
        ∂g∂x_Pg_Qg_mismatch_by_id_iq *
        id_iq_mis_kg_id_iq_by_vθ_g )

    # Pg_Qg_mis_k_δ_ed_dash_eq_dash =
    #     (∂g∂p_δ_ed_dash_eq_dash + ∂g∂x_id_iq *
    #     id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash)

    Pg_Qg_mis_k_δ_ed_dash_eq_dash =
        (∂g∂p_Pg_Qg_mismatch_by_δ_ed_dash_eq_dash +
        ∂g∂x_Pg_Qg_mismatch_by_id_iq *
        id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash)

    # Pg_Qg_mis_k_Png_Qng =
    #     ∂g∂x_vθ_ng * Png_Qng_mis_k_vθ_ng_by_Png_Qng

    Pg_Qg_mis_k_Png_Qng =
        ∂g∂x_Pg_Qg_mismatch_by_non_gens_vh_θh *
        Png_Qng_mis_k_vθ_ng_by_Png_Qng


    #
    Pg_Qg_mis_vθ_g_by_δ_ed_dash_eq_dash =
        -(Pg_Qg_mis_k_vθ_g\Pg_Qg_mis_k_δ_ed_dash_eq_dash)

    Pg_Qg_mis_vθ_g_by_Png_Qng =
        -(Pg_Qg_mis_k_vθ_g\Pg_Qg_mis_k_Png_Qng)

    if loc_load_exist == true

        Pg_Qg_mis_vθ_g_by_Pgll_Qgll =
            -(Pg_Qg_mis_k_vθ_g\∂g∂p_Pgll_Qgll)
    end

    #
    kg_id_iq_by_δ_ed_dash_eq_dash =
        (id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash +
        id_iq_mis_kg_id_iq_by_vθ_g *
        Pg_Qg_mis_vθ_g_by_δ_ed_dash_eq_dash)

    kg_id_iq_by_Png_Qng =
        id_iq_mis_kg_id_iq_by_vθ_g *
        Pg_Qg_mis_vθ_g_by_Png_Qng

    if loc_load_exist == true

         kg_id_iq_by_Pgll_Qgll =
             id_iq_mis_kg_id_iq_by_vθ_g *
             Pg_Qg_mis_vθ_g_by_Pgll_Qgll
    end

    # 
    # kf_id_iq_by_δ_ed_dash_eq_dash =
    #     ∂f∂p_id_iq *
    #     id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash 

    kf_id_iq_by_δ_ed_dash_eq_dash =
        ∂f∂p_im_state_by_gens_id_iq *
        id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash 

    # kf_id_iq_by_vθ_g =
    #     ∂f∂p_id_iq * id_iq_mis_kg_id_iq_by_vθ_g

    kf_id_iq_by_vθ_g =
        ∂f∂p_im_state_by_gens_id_iq *
        id_iq_mis_kg_id_iq_by_vθ_g

    #
    # kf_δ_ed_dash_eq_dash =
    #     ( kf_id_iq_by_δ_ed_dash_eq_dash +
    #     (kf_id_iq_by_vθ_g + ∂f∂p_vh_g_wt_θh0) *
    #     Pg_Qg_mis_vθ_g_by_δ_ed_dash_eq_dash ) 

    kf_δ_ed_dash_eq_dash =
        ( kf_id_iq_by_δ_ed_dash_eq_dash +
        (kf_id_iq_by_vθ_g +
        ∂f∂p_im_state_by_gens_vh_wt_zero_θh) *
        Pg_Qg_mis_vθ_g_by_δ_ed_dash_eq_dash ) 

    # kf_Png_Qng =
    #     (kf_id_iq_by_vθ_g + ∂f∂p_vh_g_wt_θh0 ) *
    #     Pg_Qg_mis_vθ_g_by_Png_Qng 

    kf_Png_Qng =
        (kf_id_iq_by_vθ_g +
        ∂f∂p_im_state_by_gens_vh_wt_zero_θh ) *
        Pg_Qg_mis_vθ_g_by_Png_Qng 

    if loc_load_exist == true

         kf_Pgll_Qgll =
             (kf_id_iq_by_vθ_g + ∂f∂p_vh_g_wt_θh0) *
             Pg_Qg_mis_vθ_g_by_Pgll_Qgll

    end

    #---------------------------------------------

    Asys = copy( ∂f∂x_im_state_by_im_state )

    Asys[:, flat_δ_ed_dash_eq_dash_Idxs_in_state] .=
        Asys[:, flat_δ_ed_dash_eq_dash_Idxs_in_state ] +
        kf_δ_ed_dash_eq_dash


    sys_eigvalues = eigvals(Asys)


    #---------------------------------------------

    Axp =
        Asys[im_pure_states_idxs,
             im_pure_states_idxs ]

    Axa =
        Asys[im_pure_states_idxs,
             im_τm_tilade_vf_states_idxs ]

    Aap =
        Asys[im_τm_tilade_vf_states_idxs,
             im_pure_states_idxs ]

    Aaa =
        Asys[im_τm_tilade_vf_states_idxs,
             im_τm_tilade_vf_states_idxs ]

    pure_Asys = Axp - Axa * (Aaa \ Aap )

    #---------------------------------------------

    pure_Asys_eigen_object = eigen( pure_Asys )

    pure_Asys_eig_values = pure_Asys_eigen_object.values

    printed_pure_Asys_eig_values =
        round.(pure_Asys_eig_values; digits=4)
    

    pure_Asys_eigvecs_right =
        pure_Asys_eigen_object.vectors

    # inv_pure_Asys_eigvecs_right =
    #     inv(pure_Asys_eigvecs_right)

    inv_pure_Asys_eigvecs_right =
        pure_Asys_eigvecs_right \
        LinearAlgebra.I(size(pure_Asys_eigvecs_right)[1] )
    
    #--------------------------------------------------
    # states associated with eig_values
    #--------------------------------------------------

    M_diag = inv_pure_Asys_eigvecs_right *
        pure_Asys * pure_Asys_eigvecs_right

    printed_M_diag =
        round.( M_diag; digits = 4 )

    #--------------------------------------------------
    # participation factor
    #--------------------------------------------------

    PF_pure_Asys  = get_participation_factors(pure_Asys)

    printed_PF_pure_Asys =
        round.( PF_pure_Asys; digits = 4 )

    #--------------------------------------------------

    eig_values, eigvecs_left, eigvecs_right =
        get_eigens( pure_Asys )

    #--------------------------------------------------

    im_pure_states_syms =
        im_ode_sym[im_pure_states_idxs]

    #--------------------------------------------------
    #--------------------------------------------------


    ∂f∂p_pg_qh0_1 = ∂f∂p_pg_qh0[
        im_pure_states_idxs, :]

    ∂f∂p_pg_qh0_2 = ∂f∂p_pg_qh0[
        im_τm_tilade_vf_states_idxs, :]

    ∂f∂p_pg_qh0_f =
        (∂f∂p_pg_qh0_1 - Axa *(Aaa\∂f∂p_pg_qh0_2) )


    kf_Png_Qng_1 = kf_Png_Qng[
        im_pure_states_idxs, :]

    kf_Png_Qng_2 = kf_Png_Qng[
        im_τm_tilade_vf_states_idxs, :]

    kf_Png_Qng_f = (kf_Png_Qng_1 - Axa *(Aaa\kf_Png_Qng_2) )

    ∂f∂p_ωs_ωref0_vref0_porder0_1 =
        ∂f∂p_ωs_ωref0_vref0_porder0[
            im_pure_states_idxs, :]

    ∂f∂p_ωs_ωref0_vref0_porder0_2 =
        ∂f∂p_ωs_ωref0_vref0_porder0[
            im_τm_tilade_vf_states_idxs, :]

    ∂f∂p_ωs_ωref0_vref0_porder0_f =
        (∂f∂p_ωs_ωref0_vref0_porder0_1 -
        Axa *(Aaa\∂f∂p_ωs_ωref0_vref0_porder0_2 ) )

    if loc_load_exist == true

        kf_Pgll_Qgll_1 = kf_Pgll_Qgll[
            im_pure_states_idxs, :]

        kf_Pgll_Qgll_2 = kf_Pgll_Qgll[
            im_τm_tilade_vf_states_idxs, :]

        kf_Pgll_Qgll_f =
            (kf_Pgll_Qgll_1 - Axa *(Aaa\kf_Pgll_Qgll_2) )
        
    end

    #--------------------------------------------------

    pure_states =
        sim_state_x0[ im_pure_states_idxs ]

    PQ_g = [ [ [ a_P, a_Q ] for ( a_P, a_Q ) in
                   zip(P_gens,  Q_gens) ]...;]

    PQ_ng = [ [ [ a_P, a_Q ] for ( a_P, a_Q ) in
                   zip(P_non_gens, Q_non_gens) ]...;]

    if loc_load_exist == true

        PQ_g_ll = [ [ [ a_P, a_Q ] for ( a_P, a_Q ) in
                       zip(P_g_loc_load, Q_g_loc_load) ]...;]
    end

    #--------------------------------------------------

    Δpure_states =
        similar( pure_states )

    Δx     = similar( pure_states )

    ΔPQ_g  = similar( PQ_g )

    ΔPQ_ng = similar( PQ_ng )

    Δωs_ωref0_vref0_porder0  =
        similar( flat_ωs_ωref0_vref0_porder0 )

    if loc_load_exist == true

        ΔPQll_g = similar( PQ_g_ll )
        
    end

    #--------------------------------------------------
    #--------------------------------------------------

    if loc_load_exist == true

        Δx = pure_Asys * Δpure_states +
            ∂f∂p_pg_qh0_f   * ΔPQ_g   +
            kf_Png_Qng_f  * ΔPQ_ng  +
            kf_Pgll_Qgll_f * ΔPQll_g +
            ∂f∂p_ωs_ωref0_vref0_porder0_f *
            Δωs_ωref0_vref0_porder0
        
    else
        
        Δx = pure_Asys    * Δpure_states +
            ∂f∂p_pg_qh0_f * ΔPQ_g   +
            kf_Png_Qng_f  * ΔPQ_ng  +
            ∂f∂p_ωs_ωref0_vref0_porder0_f *
            Δωs_ωref0_vref0_porder0
        
    end


    if loc_load_exist == true

        return (; pure_Asys, M_diag,
                PF_pure_Asys, eig_values,
                eigvecs_left, eigvecs_right,
                inv_pure_Asys_eigvecs_right,
                ∂f∂p_pg_qh0_f, kf_Png_Qng_f,
                ∂f∂p_ωs_ωref0_vref0_porder0_f,
                kf_Pgll_Qgll_f, im_pure_states_syms  )
        
    else
        
        return (; pure_Asys, M_diag,
                PF_pure_Asys, eig_values,
                eigvecs_left, eigvecs_right,
                inv_pure_Asys_eigvecs_right,
                ∂f∂p_pg_qh0_f, kf_Png_Qng_f,
                ∂f∂p_ωs_ωref0_vref0_porder0_f,
                im_pure_states_syms )
        
    end
    
end




function get_quasi_static_analysis_sd_dynamics_wt_intra_dyn_pf(flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para ;kwd = stability_sd_dynamics_wt_intra_dyn_pf_kwd )

    (;
     loc_load_exist,
     
     sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para,

     # flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,

     # flat_ωs_ωref0_vref0_porder0,
     
     sim_state_x0,
     # flat_vh_flat_θh_x0,
     # init_dyn_pf_flat_para,
     gens_nodes_ra_Xd_dash_Xq_dash,

     im_mass_matrix,
     im_sym, 
     ode_alg,
     dt,

     flat_vh_flat_θh_Idx,

     gens_nodes_idx,
     non_gens_nodes_idx,

     flat_δ_ed_dash_eq_dash_Idxs_in_state,
     δ_ed_dash_eq_dash_Idxs_in_flattend,

     intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd,
     intra_dyn_pf_mismatch_kwd_para,

     im_state_idx_in_∂f∂x,
     net_vh_θh_idx_in_∂f∂x,

     consecutive_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
     vec_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,     

     #
     vec_flat_vh_flat_θh_idx_in_∂g∂x,
     vec_flat_vh_θh_idx_in_∂g∂x,
     vec_id_iq_idx_in_∂g∂x,
     
     consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
     consecutive_flat_vh_θh_idx_in_∂g∂x,
     #
     
     gens_consecutive_flat_vh_θh_idx_in_∂g∂x,
     non_gens_consecutive_flat_vh_θh_idx_in_∂g∂x,

     gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
     non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
     consecutive_id_iq_idx_in_∂g∂x,

     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,

     consecutive_Pg_idx_Qg_idx_in_mismatch,
     consecutive_Png_idx_Qng_idx_in_mismatch,      
     consecutive_id_idx_iq_idx_in_mismatch,

     consecutive_dyn_gens_PQ_idx_in_Idxs,
     consecutive_dyn_non_gens_PQ_idx_in_Idxs,
     consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs,
     dyn_δ_ed_eq_pf_Idxs,   

     # ∂f∂x_xp_by_vh_θh_∂f∂x,

     im_pure_states_idxs,
     im_τm_tilade_vf_states_idxs,
     ) = kwd 

    #--------------------------------------------------

    flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx = sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para.flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx

    flat_vh_θh_idx, flat_ωs_ωref0_vref0_porder0_idx, flat_init_dyn_pf_para_idx = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx

   
    init_flat_vh_flat_θh =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[flat_vh_θh_idx]
    
    flat_ωs_ωref0_vref0_porder0 =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_init_dyn_pf_para_idx ]

    #--------------------------------------------------

    
    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    
    #--------------------------------------------------
    #--------------------------------------------------
    # sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!
    #--------------------------------------------------
    #--------------------------------------------------
    
    sim_fun_kwd_para =
        sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para

    # flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para

    sim_func  =
        sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!
    
    #-----------------------------------------------------

    sd_per_gen_ode_sep_vtf_pf_sim_sol =
        DifferentialEquations.solve(
            ODEProblem(
                ODEFunction{true}(
                    (dx, x, p, t) ->
                        sim_func(dx, x, p, t;
                                 kwd_para =
                                     sim_fun_kwd_para
                         );
                    mass_matrix =
                        im_mass_matrix ,
                    syms =
                        im_sym ),
                sim_state_x0 ,
                sim_timespan,
                sim_fun_para ),
            ode_alg,
            dt=dt )        
    

    #--------------------------------------------------
    #--------------------------------------------------
    # Method 1
    #--------------------------------------------------
    #--------------------------------------------------

    #---------------------------------------------
    #---------------------------------------------
    # sensitivity
    # sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!
    #---------------------------------------------
    #---------------------------------------------

    tt   = sim_timespan[1]
    sd_x = sd_per_gen_ode_sep_vtf_pf_sim_sol( tt )
           
   sd_dx = similar(sd_x)
    
    sim_fun_kwd_para =
        sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para

    sim_func  =
        sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!

    #---------------------------------------------
        
    sd_per_gen_ode_sep_vtf_pf_∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> sim_func(
            dx, x, sim_fun_para, tt;
            kwd_para = sim_fun_kwd_para ),
        sd_dx, sd_x )
    
    sd_per_gen_ode_sep_vtf_pf_∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> sim_func(
            dx, sd_x, p, tt;
            kwd_para = sim_fun_kwd_para ),
        sd_dx, sim_fun_para )

    sd_per_gen_ode_sep_vtf_pf_dxdp =
        -( svd( sd_per_gen_ode_sep_vtf_pf_∂f∂x )\
        sd_per_gen_ode_sep_vtf_pf_∂f∂p )


    #---------------------------------------------
    #---------------------------------------------
    # sensitivity
    # intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq
    #---------------------------------------------
    #---------------------------------------------


    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------

    #----------------------------------------    

    (flat_vh_idx_in_flat_vh_flat_θh,
     flat_θh_idx_in_flat_vh_flat_θh) =
         flat_vh_flat_θh_Idx
    
    flat_vh = init_flat_vh_flat_θh[
        flat_vh_idx_in_flat_vh_flat_θh ]

    flat_θh = init_flat_vh_flat_θh[
        flat_θh_idx_in_flat_vh_flat_θh ]
    
    gens_vh = flat_vh[ gens_nodes_idx ]

    gens_θh = flat_θh[ gens_nodes_idx ]        

    #----------------------------------------    

    flat_δ_ed_dash_eq_dash =
        sd_x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    flat_δ_ed_dash_eq_dash_x0 =
        sim_state_x0[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    #----------------------------------------    
    #----------------------------------------    

    gens_δ_ed_dash_eq_dash = [
        flat_δ_ed_dash_eq_dash[idx]
        for idx in
            δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    
    gens_δ =
        first.( gens_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_δ_ed_dash_eq_dash )

    #----------------------------------------

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh, a_θh, a_δ,
            a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
        for (a_vh, a_θh,
             a_δ, a_ed_dash, a_eq_dash,
             a_ra, a_X_d_dash, a_X_q_dash ) in
            zip( gens_vh, gens_θh,
                 gens_δ, gens_ed_dash, gens_eq_dash,
                 gens_ra, gens_Xd_dash, gens_Xq_dash ) ]

    gens_i_d =
        first.( gens_id_iq )
    
    gens_i_q =
        last.( gens_id_iq )

    #----------------------------------------    
    
    init_flat_vh_flat_θh_cal_id_iq =
        [flat_vh; flat_θh;
         gens_i_d; gens_i_q  ]

    #----------------------------------------   

    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    if loc_load_exist == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens; Q_gens; P_non_gens;
                Q_non_gens; flat_δ_ed_dash_eq_dash;
                P_g_loc_load; Q_g_loc_load ]                  

    else
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens; Q_gens; P_non_gens;
                Q_non_gens; flat_δ_ed_dash_eq_dash ]   

    end
    
    #----------------------------------------
    #----------------------------------------

    intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_sol =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_func!(
            init_flat_vh_flat_θh_cal_id_iq,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            kwd_para =
                intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd )

    #----------------------------------------

    intra_dyn_id_iq_δ_pf_x =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_sol.u

    intra_dyn_id_iq_δ_pf_ΔPQ_id_iq =
        similar( intra_dyn_id_iq_δ_pf_x )
    

    intra_dyn_id_iq_δ_pf_fun =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch


    intra_dyn_id_iq_δ_pf_∂g∂x = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq, x ) ->
            intra_dyn_id_iq_δ_pf_fun( pf_ΔPQ_id_iq, x,
                intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_dyn_id_iq_δ_pf_ΔPQ_id_iq,
        intra_dyn_id_iq_δ_pf_x  )

    
    intra_dyn_id_iq_δ_pf_∂g∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq,  p ) ->
            intra_dyn_id_iq_δ_pf_fun(
                pf_ΔPQ_id_iq,
                intra_dyn_id_iq_δ_pf_x,
                p;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_dyn_id_iq_δ_pf_ΔPQ_id_iq,
        intra_dyn_pf_mismatch_flat_para )

       intra_dyn_id_iq_δ_pf_∂x∂p =
           -(svd( intra_dyn_id_iq_δ_pf_∂g∂x )\
            intra_dyn_id_iq_δ_pf_∂g∂p)

    #---------------------------------------------
    # sd_per_gen_ode_sep_vtf_pf_∂f∂x
    #---------------------------------------------

    # xp : pure states

    # xa : algebraic states

    ∂f∂x_xp_by_xp_∂f∂x = Ax_states_by_states_∂f∂x =
        sd_per_gen_ode_sep_vtf_pf_∂f∂x[
            im_state_idx_in_∂f∂x,
            im_state_idx_in_∂f∂x ]

    ∂f∂x_xp_by_vh_θh_∂f∂x = Ax_states_by_vh_θh_∂f∂x = 
        sd_per_gen_ode_sep_vtf_pf_∂f∂x[
            im_state_idx_in_∂f∂x,
            net_vh_θh_idx_in_∂f∂x ]

    ∂f∂x_vh_θh_by_xp_∂f∂x = Ax_vh_θh_by_states_∂f∂x =
        sd_per_gen_ode_sep_vtf_pf_∂f∂x[
            net_vh_θh_idx_in_∂f∂x,
            im_state_idx_in_∂f∂x]

    ∂f∂x_vh_θh_by_vh_θh_∂f∂x = Ax_vh_θh_by_vh_θh_∂f∂x = 
        sd_per_gen_ode_sep_vtf_pf_∂f∂x[
            net_vh_θh_idx_in_∂f∂x,
            net_vh_θh_idx_in_∂f∂x]

    #---------------------------------------------
    # sd_per_gen_ode_sep_vtf_pf_∂f∂p
    #---------------------------------------------

    ∂f∂p_xp_by_ωs_ωref0_vref0_porder0_∂f∂p =
        Ax_states_by_ωs_ωref0_vref0_porder0_∂f∂p =
        sd_per_gen_ode_sep_vtf_pf_∂f∂p[
            im_state_idx_in_∂f∂x,
            :]

    ∂f∂p_xa_by_ωs_ωref0_vref0_porder0_∂f∂p =
        Ax_vh_θh_by_ωs_ωref0_vref0_porder0_∂f∂p =
            sd_per_gen_ode_sep_vtf_pf_∂f∂p[
                net_vh_θh_idx_in_∂f∂x,
                :]

    #---------------------------------------------
    # intra_dyn_id_iq_δ_pf_∂g∂x
    #---------------------------------------------

    ∂g∂x_id_iq_by_id_iq =
        Mx_gens_id_iq_δ_by_gens_id_iq_δ_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            consecutive_id_iq_idx_in_∂g∂x ]

    ∂g∂x_id_iq_by_vθ_g =
        Mx_gens_id_iq_δ_by_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x]

    ∂g∂x_id_iq_by_vθ_ng =
        Mx_gens_id_iq_δ_by_non_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    #---------------------------------------------

    ∂g∂x_g_by_id_iq =
        Mx_gens_PQ_by_gens_id_iq_δ_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch ,
            consecutive_id_iq_idx_in_∂g∂x ]

    ∂g∂x_g_by_vθ_g =
        Mx_gens_PQ_by_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch ,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]


    ∂g∂x_g_by_vθ_ng =
        Mx_gens_PQ_by_non_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch ,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]


    #---------------------------------------------

    ∂g∂x_ng_by_id_iq =
        Mx_non_gens_PQ_by_gens_id_iq_δ_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch ,
            consecutive_id_iq_idx_in_∂g∂x ]

    ∂g∂x_ng_by_vθ_g =
        Mx_non_gens_PQ_by_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch ,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]


    ∂g∂x_ng_by_vθ_ng =
        Mx_non_gens_PQ_by_non_gens_vh_θh_∂g∂x =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch ,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    #---------------------------------------------
    # intra_dyn_id_iq_δ_pf_∂g∂p
    #---------------------------------------------

    ∂g∂p_id_iq_by_gens_PQ =
        Mx_gens_id_iq_δ_by_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_id_idx_iq_idx_in_mismatch ,
             consecutive_dyn_gens_PQ_idx_in_Idxs ]

    ∂g∂p_id_iq_by_non_gens_PQ =
        Mx_gens_id_iq_δ_by_non_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_id_idx_iq_idx_in_mismatch,
            consecutive_dyn_non_gens_PQ_idx_in_Idxs ]


    ∂g∂p_id_iq_by_δ_ed_eq =
        Mx_gens_id_iq_δ_by_dyn_δ_ed_eq_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_id_idx_iq_idx_in_mismatch,
            dyn_δ_ed_eq_pf_Idxs ]


    if loc_load_exist == true
        
        ∂g∂p_id_iq_by_PQll_g =
            Mx_gens_id_iq_δ_by_gens_loc_load_PQ_∂g∂p =
            intra_dyn_id_iq_δ_pf_∂g∂p[
                consecutive_id_idx_iq_idx_in_mismatch,
                consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs ]

    end

    #---------------------------------------------

    ∂g∂p_g_by_PQ_g = Mx_gens_PQ_by_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            consecutive_dyn_gens_PQ_idx_in_Idxs ]

    ∂g∂p_g_by_PQ_ng = Mx_gens_PQ_by_non_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            consecutive_dyn_non_gens_PQ_idx_in_Idxs ]

    ∂g∂p_g_by_δ_ed_eq = Mx_gens_PQ_by_dyn_δ_ed_eq_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            dyn_δ_ed_eq_pf_Idxs ]


    if loc_load_exist == true
        
        ∂g∂p_g_by_PQll_g =
            Mx_gens_PQ_by_gens_loc_load_PQ_∂g∂p =
            intra_dyn_id_iq_δ_pf_∂g∂p[
                consecutive_Pg_idx_Qg_idx_in_mismatch,
                consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs ]

    end

    #---------------------------------------------

    ∂g∂p_ng_by_PQ_g = Mx_non_gens_PQ_by_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            consecutive_dyn_gens_PQ_idx_in_Idxs ]

    ∂g∂p_ng_by_PQ_ng = Mx_non_gens_PQ_by_non_gens_PQ_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            consecutive_dyn_non_gens_PQ_idx_in_Idxs ]

    ∂g∂p_ng_by_δ_ed_eq = Mx_non_gens_PQ_by_dyn_δ_ed_eq_∂g∂p =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            dyn_δ_ed_eq_pf_Idxs ]

    if loc_load_exist == true
        
         ∂g∂p_ng_by_PQll_g =
             Mx_non_gens_PQ_by_gens_loc_load_PQ_∂g∂p =
            intra_dyn_id_iq_δ_pf_∂g∂p[
                consecutive_Png_idx_Qng_idx_in_mismatch,
                consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs ]

    end

    #---------------------------------------------

    # Define Q_vθ_g as:
    
    Q_vθ_g = ( ∂g∂x_g_by_vθ_g -
        ∂g∂x_g_by_vθ_ng *
        ( ∂g∂x_ng_by_vθ_ng\∂g∂x_ng_by_vθ_g) -
            ∂g∂x_g_by_id_iq *
            (∂g∂x_id_iq_by_id_iq\∂g∂x_id_iq_by_vθ_g)) 

    # Define Q_δ_ed_eq as:
    
    Q_δ_ed_eq =
        (∂g∂p_g_by_δ_ed_eq -
        ∂g∂x_g_by_id_iq *
        (∂g∂x_id_iq_by_id_iq\
        ∂g∂p_id_iq_by_δ_ed_eq ))

    Asys =  copy(∂f∂x_xp_by_xp_∂f∂x)

    Asys[:, flat_δ_ed_dash_eq_dash_Idxs_in_state ] .=
        Asys[:, flat_δ_ed_dash_eq_dash_Idxs_in_state] -
        ∂f∂x_xp_by_vh_θh_∂f∂x[
            :, gens_consecutive_flat_vh_θh_idx_in_∂g∂x] *
                ( Q_vθ_g \ Q_δ_ed_eq ) 

    #---------------------------------------------

    Axp =
        Asys[im_pure_states_idxs,
             im_pure_states_idxs ]

    Axa =
        Asys[im_pure_states_idxs,
             im_τm_tilade_vf_states_idxs ]

    Aap =
        Asys[im_τm_tilade_vf_states_idxs,
             im_pure_states_idxs ]

    Aaa =
        Asys[im_τm_tilade_vf_states_idxs,
             im_τm_tilade_vf_states_idxs ]

    
    pure_Asys = Axp - Axa * (Aaa \ Aap )

    #---------------------------------------------

    eigen_object = eigen( pure_Asys )

    eig_values = eigen_object.values

    printed_eig_values =
        round.(eig_values; digits=4)


    eigvecs_right = eigen_object.vectors

    # inv_eigvecs_right = inv(eigvecs_right)
    
    inv_eigvecs_right =
        eigvecs_right \ LinearAlgebra.I( size(eigvecs_right)[1] )

    #--------------------------------------------------
    # states associated with eig_values
    #--------------------------------------------------

    M_diag = inv_eigvecs_right * pure_Asys * eigvecs_right

    printed_M_diag =
        round.( M_diag; digits = 4 )

    # tup_state_var_eig_value =
    #     [(a_state_sym, an_eig_value)
    #      for ( a_state_sym, an_eig_value )  in
    #          zip( im_ode_sym, printed_eig_values )]

    #--------------------------------------------------
    # participation factor
    #--------------------------------------------------

    PF_Asys  = get_participation_factors( pure_Asys )

    printed_PF_Asys =
        round.( PF_Asys; digits = 4 )

    #--------------------------------------------------

    eig_values, eigvecs_left, eigvecs_right =
        get_eigens(pure_Asys)

    #--------------------------------------------------
    #--------------------------------------------------


    Δxp = similar( sim_state_x0[ im_state_idx_in_∂f∂x ])

    Δvθ = similar( flat_vh_flat_θh_x0[
        consecutive_flat_vh_flat_θh_idx_in_∂g∂x ] )

    Δid_iq = similar(
        gens_id_iq[ consecutive_id_iq_idx_in_∂g∂x ] )

    Δωs_ωref_vref_porder =
        similar( flat_ωs_ωref0_vref0_porder0 )

    Δvθ_g = similar( flat_vh_flat_θh_x0[
        gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ])
    
    Δvθ_ng = similar(flat_vh_flat_θh_x0[
        non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ])

    ΔPQ_g = [[[a_P, a_Q]
            for (a_P, a_Q) in
                zip(P_gens, Q_gens)]...;]

    ΔPQ_ng = [[[a_P, a_Q]
            for (a_P, a_Q) in
                zip(P_non_gens, Q_non_gens)]...;]

    Δδ_ed_eq = similar( flat_δ_ed_dash_eq_dash )

    if loc_load_exist == true
        
        ΔPQll_g = [[[a_P, a_Q]
                for (a_P, a_Q) in
                    zip(P_g_loc_load, Q_g_loc_load)]...;]

    else

        ΔPQll_g = similar( zero( length(ΔPQ_g) ) )
        
    end

   
    #--------------------------------------------------
    
    Δf_xp =
        ∂f∂x_xp_by_xp_∂f∂x * Δxp +
        ∂f∂x_xp_by_vh_θh_∂f∂x * Δvθ +
        (∂f∂p_xp_by_ωs_ωref0_vref0_porder0_∂f∂p *
        Δωs_ωref_vref_porder )

    # 0 = Δf_vθ
    Δf_vθ = 
        ∂f∂x_vh_θh_by_xp_∂f∂x * Δxp +
        ∂f∂x_vh_θh_by_vh_θh_∂f∂x * Δvθ +
        (∂f∂p_xa_by_ωs_ωref0_vref0_porder0_∂f∂p *
        Δωs_ωref_vref_porder )


    Δvθ  =
        -(∂f∂x_vh_θh_by_vh_θh_∂f∂x\
        ∂f∂x_vh_θh_by_xp_∂f∂x) * Δxp -
        (∂f∂x_vh_θh_by_vh_θh_∂f∂x\
        ∂f∂p_xa_by_ωs_ωref0_vref0_porder0_∂f∂p) *
        Δωs_ωref_vref_porder


    #---------------------------------------------
    # intra_dyn_id_iq_δ_pf_∂g∂x
    #---------------------------------------------


    """
    0 = ΔPQ_id_iq =
        ∂g∂x * Δvh_θh_id_ig +
        ∂g∂p * ΔPg_Qg_Png_Qng_δ_ed_eq_Pgll_Qgll
    """
    #---------------------------------------------
    
    ΔPQ_id_iq_rows_g =
        ∂g∂x_g_by_vθ_g  * Δvθ_g  +
        ∂g∂x_g_by_vθ_ng * Δvθ_ng +
        ∂g∂x_g_by_id_iq * Δid_iq +

        ∂g∂p_g_by_PQ_g    * ΔPQ_g +
        ∂g∂p_g_by_δ_ed_eq * Δδ_ed_eq +
        ∂g∂p_g_by_PQll_g  * ΔPQll_g
        
    ΔPQ_id_iq_rows_ng =
        ∂g∂x_ng_by_vθ_g  * Δvθ_g +
        ∂g∂x_ng_by_vθ_ng * Δvθ_ng +

        ∂g∂p_ng_by_PQ_ng * ΔPQ_ng
        
    ΔPQ_id_iq_rows_id_iq =
        ∂g∂x_id_iq_by_vθ_g   * Δvθ_g +
        ∂g∂x_id_iq_by_id_iq  * Δid_iq +
        
        ∂g∂p_id_iq_by_δ_ed_eq * Δδ_ed_eq

    #---------------------------------------------

    Δid_iq =
        -(∂g∂x_id_iq_by_id_iq\∂g∂x_id_iq_by_vθ_g ) * Δvθ_g -
        (∂g∂x_id_iq_by_id_iq\∂g∂p_id_iq_by_by_δ_ed_eq ) *
        Δδ_ed_eq

    vθ_ng =
        -(∂g∂x_ng_by_vθ_ng \ ∂g∂x_ng_by_vθ_g) * Δvθ_g -
        (∂g∂x_ng_by_vθ_ng \ ∂g∂p_ng_by_PQ_ng) * ΔPQ_ng

    ΔPQ_g =
        ( ∂f∂x_g_by_vθ_g -
        ∂g∂x_g_by_vθ_ng * ( ∂g∂x_ng_by_vθ_ng\∂g∂x_ng_by_vθ_g) -
        ∂g∂x_g_by_id_iq *
        (∂g∂x_id_iq_by_id_iq\∂g∂x_id_iq_by_vθ_g)) * Δvθ_g +

        (∂g∂p_g_by_δ_ed_eq - ∂g∂x_g_by_id_iq *
        (∂g∂x_id_iq_by_id_iq\∂g∂p_id_iq_by_δ_ed_eq ))  *
        Δδ_ed_eq +
        ∂g∂p_g_by_PQ_g * ΔPQ_g +
        ∂g∂p_g_by_PQll_g  * ΔPQll_g -
        (∂g∂x_ng_by_vθ_ng \ ∂g∂p_ng_by_PQ_ng) * ΔPQ_ng

   # Define  Q_vθ_g as:
    
    Q_vθ_g = ( ∂f∂x_g_by_vθ_g -
        ∂g∂x_g_by_vθ_ng *
        ( ∂g∂x_ng_by_vθ_ng\∂g∂x_ng_by_vθ_g) -
            ∂g∂x_g_by_id_iq *
            (∂g∂x_id_iq_by_id_iq\∂g∂x_id_iq_by_vθ_g)) 

   # Define  Q_δ_ed_eq as:
    
    Q_δ_ed_eq =
        (∂g∂p_g_by_δ_ed_eq -
        ∂g∂x_g_by_id_iq *
        (∂g∂x_id_iq_by_id_iq\
        ∂g∂p_id_iq_by_δ_ed_eq ))

   # 0 =  ΔPQ_g

   ΔPQ_g = Q_vθ_g * Δvθ_g  + Q_δ_ed_eq * Δδ_ed_eq +
        
            ∂g∂p_g_by_PQ_g * ΔPQ_g +
            ∂g∂p_g_by_PQll_g  * ΔPQll_g -
            (∂g∂x_ng_by_vθ_ng \ ∂g∂p_ng_by_PQ_ng) * ΔPQ_ng

   Δvθ_g = -(Q_vθ_g\Q_δ_ed_eq) * Δδ_ed_eq -         
            (Q_vθ_g\∂g∂p_g_by_PQ_g) * ΔPQ_g -
            (Q_vθ_g\∂g∂p_g_by_PQll_g)  * ΔPQll_g +
            (Q_vθ_g\(∂g∂x_ng_by_vθ_ng \ ∂g∂p_ng_by_PQ_ng)) *
            ΔPQ_ng

    #---------------------------------------------

    Δf_xp =
        ∂f∂x_xp_by_xp_∂f∂x * Δxp +
        ∂f∂x_xp_by_vh_θh_∂f∂x * Δvθ +
        (∂f∂p_xp_by_ωs_ωref0_vref0_porder0_∂f∂p *
        Δωs_ωref_vref_porder )


    Δf_xp =
        ∂f∂x_xp_by_xp_∂f∂x * Δxp -
        ∂f∂x_xp_by_vh_θh_∂f∂x * (
            Q_vθ_g\Q_δ_ed_eq) * Δδ_ed_eq + 

        (∂f∂p_xp_by_ωs_ωref0_vref0_porder0_∂f∂p *
        Δωs_ωref_vref_porder ) -
        
        ∂f∂x_xp_by_vh_θh_∂f∂x * (
            Q_vθ_g\∂g∂p_g_by_PQ_g) * ΔPQ_g - 

        ∂f∂x_xp_by_vh_θh_∂f∂x * (
            Q_vθ_g\∂g∂p_g_by_PQll_g) * ΔPQll_g +

        ∂f∂x_xp_by_vh_θh_∂f∂x * (
            Q_vθ_g\( ∂g∂x_ng_by_vθ_ng \
                ∂g∂p_ng_by_PQ_ng)) * ΔPQ_ng 
    
    #--------------------------------------------------       #--------------------------------------------------
    
end



function get_quasi_static_analysis_ode_per_gen_wt_intra_dyn_pf(init_flat_vh_flat_θh, init_dyn_pf_flat_para, flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh; kwd = stability_ode_per_gen_wt_intra_dyn_pf_kwd )

    (;
     loc_load_exist,
     
     sim_state_x0,
     # flat_vh_flat_θh_x0,

     im_vars_Idx_in_state,

     ode_only_per_gen_id_iq_pg_vh_kwd,

     # flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh,

     # flat_ωs_ωref0_vref0_porder0,

     im_ode_mass_matrix,
     im_ode_sym,

     gens_nodes_ra_Xd_dash_Xq_dash,

     flat_vh_flat_θh_Idx,

     gens_nodes_idx,
     non_gens_nodes_idx,

     flat_δ_ed_dash_eq_dash_Idxs_in_state,
     δ_ed_dash_eq_dash_Idxs_in_flattend,

     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,

     # init_dyn_pf_flat_para,
     intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd,
     intra_dyn_pf_mismatch_kwd_para,

     im_state_idx_in_∂f∂x,

     flat_id_iq_pg_vh_idx_in_∂f∂p,
     flat_ωs_ωref0_vref0_porder0_idx_in_∂f∂p,
     gens_id_iq_idx_in_∂f∂p,
     gens_vh_idx_in_∂f∂p,

     gens_ph_idx_in_∂f∂p,

     consecutive_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
     vec_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0,
     

     #
     vec_flat_vh_flat_θh_idx_in_∂g∂x,
     vec_flat_vh_θh_idx_in_∂g∂x,
     vec_id_iq_idx_in_∂g∂x,

     consecutive_flat_vh_flat_θh_idx_in_∂g∂x,
     consecutive_flat_vh_θh_idx_in_∂g∂x,
     consecutive_id_iq_idx_in_∂g∂x,
     #         
     
     consecutive_gens_vh_idx_θh_idx_in_Idx,     

     gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,     
     non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x,

     
     intra_flat_δ_ed_dash_eq_dash_idx_in_∂g∂p,

     consecutive_Pg_idx_Qg_idx_in_mismatch,
     consecutive_Png_idx_Qng_idx_in_mismatch,
     consecutive_id_idx_iq_idx_in_mismatch,

     im_pure_states_idxs,
     im_τm_tilade_vf_states_idxs
     ) = kwd

    #-----------------------------------------------------

    flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_only_per_gen_id_iq_pg_vh_kwd.flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx


    flat_ωs_ωref0_vref0_porder0_Idx, flat_id_iq_pg_vh_Idx =
        flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx

        
     flat_ωs_ωref0_vref0_porder0 =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh[
            flat_ωs_ωref0_vref0_porder0_Idx ]
    
    flat_id_iq_pg_vh =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh[
            flat_id_iq_pg_vh_Idx ]
    
    
    #-----------------------------------------------------

    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs

    (;
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
     ) = dyn_pf_P_Q_δ_etc_kwd_para_Idxs 

    
    gens_PQ_idx_in_Idxs =
        [ dyn_P_gens_Idxs,
          dyn_Q_gens_Idxs ]

    
    ( consecutive_dyn_gens_PQ_idx_in_Idxs,
      vec_dyn_gens_PQ_idx_in_Idxs ) =
            convert_to_consecutive_idxs(
                gens_PQ_idx_in_Idxs )

    consecutive_Pg_Qg_idx_in_∂g∂p =
        consecutive_dyn_gens_PQ_idx_in_Idxs
    
    non_gens_PQ_idx_in_Idxs =
        [ dyn_P_non_gens_Idxs,
          dyn_Q_non_gens_Idxs ]

    ( consecutive_dyn_non_gens_PQ_idx_in_Idxs,
      vec_dyn_non_gens_PQ_idx_in_Idxs ) =
            convert_to_consecutive_idxs(
                non_gens_PQ_idx_in_Idxs )

    consecutive_Png_Qng_idx_in_∂g∂p =
        consecutive_dyn_non_gens_PQ_idx_in_Idxs
    
    if loc_load_exist == true
           
        gens_loc_load_PQ_idx_in_Idxs =
            [ dyn_P_g_loc_load_Idxs,
              dyn_Q_g_loc_load_Idxs ]

        ( consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs,
          vec_dyn_gens_loc_load_PQ_idx_in_Idxs ) =
                convert_to_consecutive_idxs(
                    gens_loc_load_PQ_idx_in_Idxs )
        
        consecutive_Pgll_Qgll_idx_in_∂g∂p =
            consecutive_dyn_gens_loc_load_PQ_idx_in_Idxs

    end
    
    #--------------------------------------------------
    # ode_only_per_gen_model_by_id_iq_pg_vh_func!
    #--------------------------------------------------

    ode_only_sim_state_x0 =
        sim_state_x0[ im_vars_Idx_in_state ]

    sim_fun_kwd_para =
        ode_only_per_gen_id_iq_pg_vh_kwd
    
    sim_fun_para =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh 
    
    sim_func  =
        ode_only_per_gen_model_by_id_iq_pg_vh_func!
    
    #-----------------------------------------------------

    ode_only_per_gen_by_id_iq_pg_vh_sim_sol =
        DifferentialEquations.solve(
        ODEProblem(
            ODEFunction{true}(
                (dx, x, p, t) ->
                    sim_func(dx, x, p, t;
                             kwd_para =
                                 sim_fun_kwd_para 
                         );
                mass_matrix =
                    # im_mass_matrix,
                    im_ode_mass_matrix,
                syms =
                    # im_sym
                    im_ode_sym ),

            # sim_state_x0,            
            # ode_sim_state_x0,
            ode_only_sim_state_x0,
            sim_timespan,
            sim_fun_para ),
            ode_alg,
            dt = dt)
    
    #---------------------------------------------
    # sensitivity
    #---------------------------------------------

    tt   = sim_timespan[1]
    sd_x = ode_only_per_gen_by_id_iq_pg_vh_sim_sol( tt )
 
    sd_dx = similar(sd_x)
    
    sim_fun_kwd_para =
        ode_only_per_gen_id_iq_pg_vh_kwd
            
    sim_fun_para =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh 
        
    sim_func  =
        ode_only_per_gen_model_by_id_iq_pg_vh_func!
        
    #---------------------------------------------

    ode_only_per_gen_by_id_iq_pg_vh_∂f∂x =
        ForwardDiff.jacobian(
        (dx, x ) -> sim_func(
            dx, x, sim_fun_para, tt;
            kwd_para = sim_fun_kwd_para ),
        sd_dx, sd_x )


    ode_only_per_gen_by_id_iq_pg_vh_∂f∂p =
        ForwardDiff.jacobian(
        (dx, p ) -> sim_func(
            dx, sd_x, p, tt;
            kwd_para = sim_fun_kwd_para ),
        sd_dx, sim_fun_para )

    #---------------------------------------------
    # intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq sensitivity
    #---------------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------
    #----------------------------------------    

    (flat_vh_idx_in_flat_vh_flat_θh,
     flat_θh_idx_in_flat_vh_flat_θh) =
         flat_vh_flat_θh_Idx
    
    flat_vh = init_flat_vh_flat_θh[
        flat_vh_idx_in_flat_vh_flat_θh ]

    flat_θh = init_flat_vh_flat_θh[
        flat_θh_idx_in_flat_vh_flat_θh ]
    
    gens_vh = flat_vh[ gens_nodes_idx ]

    gens_θh = flat_θh[ gens_nodes_idx ]        

    #----------------------------------------    

    flat_δ_ed_dash_eq_dash =
        sd_x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]

    flat_δ_ed_dash_eq_dash_x0 =
        sim_state_x0[flat_δ_ed_dash_eq_dash_Idxs_in_state]

    #----------------------------------------    
    #----------------------------------------    

    gens_δ_ed_dash_eq_dash = [
        flat_δ_ed_dash_eq_dash[idx]
        for idx in
            δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    
    gens_δ =
        first.( gens_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_δ_ed_dash_eq_dash )

    #----------------------------------------

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh, a_θh, a_δ,
            a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
        for (a_vh, a_θh,
             a_δ, a_ed_dash, a_eq_dash,
             a_ra, a_X_d_dash, a_X_q_dash ) in
            zip( gens_vh, gens_θh,
                 gens_δ, gens_ed_dash, gens_eq_dash,
                 gens_ra, gens_Xd_dash, gens_Xq_dash ) ]

    gens_i_d =
        first.( gens_id_iq )
    
    gens_i_q =
        last.( gens_id_iq )

    #----------------------------------------    
    
    init_flat_vh_flat_θh_cal_id_iq =
        [flat_vh; flat_θh;
         gens_i_d; gens_i_q  ]

    #----------------------------------------   

    P_gens =
        init_dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        init_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        init_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        init_dyn_pf_flat_para[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            init_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            init_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end


    if loc_load_exist == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens; Q_gens; P_non_gens;
                Q_non_gens; flat_δ_ed_dash_eq_dash;
                P_g_loc_load; Q_g_loc_load ]
        

    else
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens; Q_gens; P_non_gens;
                Q_non_gens; flat_δ_ed_dash_eq_dash ]   

    end
    
    #----------------------------------------
    #----------------------------------------

    intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_sol =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_func!(
            init_flat_vh_flat_θh_cal_id_iq,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            kwd_para =
                intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd)

    #----------------------------------------

    intra_dyn_id_iq_δ_pf_x =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_sol.u

    intra_dyn_id_iq_δ_pf_ΔPQ_id_iq =
        similar( intra_dyn_id_iq_δ_pf_x )
    

    intra_dyn_id_iq_δ_pf_fun =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch


    intra_dyn_id_iq_δ_pf_∂g∂x = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq, x ) ->
            intra_dyn_id_iq_δ_pf_fun( pf_ΔPQ_id_iq, x,
                intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_dyn_id_iq_δ_pf_ΔPQ_id_iq,
        intra_dyn_id_iq_δ_pf_x  )

    
    intra_dyn_id_iq_δ_pf_∂g∂p = ForwardDiff.jacobian(
        ( pf_ΔPQ_id_iq,  p ) ->
            intra_dyn_id_iq_δ_pf_fun(
                pf_ΔPQ_id_iq,
                intra_dyn_id_iq_δ_pf_x,
                p;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
        intra_dyn_id_iq_δ_pf_ΔPQ_id_iq,
        intra_dyn_pf_mismatch_flat_para )

       intra_dyn_id_iq_δ_pf_∂x∂p =
           -(svd( intra_dyn_id_iq_δ_pf_∂g∂x )\
            intra_dyn_id_iq_δ_pf_∂g∂p)

    #----------------------------------------
    # ode f
    #----------------------------------------

    # A1
    ∂f∂x_im_state_by_im_state =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂x[
            im_state_idx_in_∂f∂x,
            im_state_idx_in_∂f∂x ]
    

    # B1_gens_id_iq
    ∂f∂p_im_state_by_gens_id_iq =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            gens_id_iq_idx_in_∂f∂p ]

    # F
    ∂f∂p_pg = ∂f∂p_im_state_by_gens_ph =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            gens_ph_idx_in_∂f∂p ]

    ∂f∂p_pg_extended_wt_zero_qh =
        hcat( ∂f∂p_pg, zeros( size(∂f∂p_pg )) )

    ∂f∂p_pg_qh0 =
        ∂f∂p_pg_extended_wt_zero_qh[
            :, consecutive_flat_Pg_flat_Qg_idx_in_∂f∂p_pg_qh0 ]
    #---------------------------------------------

    # B2_gens_vh 
    ∂f∂p_im_state_by_gens_vh =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            gens_vh_idx_in_∂f∂p ]

    # B2_gens_vh_extended_wt_zero_θh
    ∂f∂p_im_state_by_gens_vh_extended_wt_zero_θh =
        hcat(∂f∂p_im_state_by_gens_vh,
             zeros(size(∂f∂p_im_state_by_gens_vh)) )

    # B2
    # B2_gens_vh_wt_zero_θh
    ∂f∂p_im_state_by_gens_vh_wt_zero_θh =
        ∂f∂p_im_state_by_gens_vh_extended_wt_zero_θh[
            :, consecutive_gens_vh_idx_θh_idx_in_Idx ]

    #---------------------------------------------

    # E_∂f∂p_ωs_ωref0_vref0_porder0
    ∂f∂p_im_state_by_ωs_ωref0_vref0_porder0 =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            flat_ωs_ωref0_vref0_porder0_idx_in_∂f∂p ]

    ∂f∂p_ωs_ωref0_vref0_porder0 =
        ode_only_per_gen_by_id_iq_pg_vh_∂f∂p[
            im_state_idx_in_∂f∂x,
            flat_ωs_ωref0_vref0_porder0_idx_in_∂f∂p ]

    #----------------------------------------
    # pf ∂g∂x
    #----------------------------------------

    # D4
    ∂g∂x_Pg_Qg_mismatch_by_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    # D5
    ∂g∂x_Pg_Qg_mismatch_by_non_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    # D3
    ∂g∂x_Pg_Qg_mismatch_by_id_iq =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            consecutive_id_iq_idx_in_∂g∂x ]

    #---------------------------------------------

    # D6
    ∂g∂x_Png_Qng_mismatch_by_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    # D7
    ∂g∂x_Png_Qng_mismatch_by_non_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    ∂g∂x_Png_Qng_mismatch_by_id_iq =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_Png_idx_Qng_idx_in_mismatch,
            consecutive_id_iq_idx_in_∂g∂x ]

    #---------------------------------------------

    # D2
    ∂g∂x_id_iq_mismatch_by_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]

    ∂g∂x_id_iq_mismatch_by_non_gens_vh_θh =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            non_gens_consecutive_flat_vh_flat_θh_idx_in_∂g∂x ]       
    # D1
    ∂g∂x_id_iq_mismatch_by_id_iq =
        intra_dyn_id_iq_δ_pf_∂g∂x[
            consecutive_id_idx_iq_idx_in_mismatch,
            consecutive_id_iq_idx_in_∂g∂x ]

    #----------------------------------------
    # pf ∂g∂p
    #----------------------------------------

    # C2
    ∂g∂p_Pg_Qg_mismatch_by_δ_ed_dash_eq_dash =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Pg_idx_Qg_idx_in_mismatch,
            intra_flat_δ_ed_dash_eq_dash_idx_in_∂g∂p ]

    # C1
    ∂g∂p_id_iq_mismatch_by_δ_ed_dash_eq_dash =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_id_idx_iq_idx_in_mismatch,
            intra_flat_δ_ed_dash_eq_dash_idx_in_∂g∂p ]

    ∂g∂p_Png_Qng_mismatch_by_Png_Qng =
        intra_dyn_id_iq_δ_pf_∂g∂p[
            consecutive_Png_idx_Qng_idx_in_mismatch,
             consecutive_Png_Qng_idx_in_∂g∂p ]


    #---------------------------------------------


    #
    # id_iq_mis_kg_id_iq_by_vθ_g =
    #     -(∂g∂x_id_iq \ ∂g∂x_vθ_g)

    id_iq_mis_kg_id_iq_by_vθ_g =
        -(∂g∂x_id_iq_mismatch_by_id_iq \
        ∂g∂x_id_iq_mismatch_by_gens_vh_θh)


    # id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash =
    #     -(∂g∂x_id_iq \
    #         ∂g∂p_δ_ed_dash_eq_dash)


    id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash =
        -(∂g∂x_id_iq_mismatch_by_id_iq \
            ∂g∂p_id_iq_mismatch_by_δ_ed_dash_eq_dash)

    #
    # Png_Qng_mis_k_vθ_ng_by_vθ_g =
    #     -(∂g∂x_vθ_ng \∂g∂x_vθ_g)

    Png_Qng_mis_k_vθ_ng_by_vθ_g =
        -(∂g∂x_Png_Qng_mismatch_by_non_gens_vh_θh \
        ∂g∂x_Png_Qng_mismatch_by_gens_vh_θh)


    # Png_Qng_mis_k_vθ_ng_by_Png_Qng =
    #     -(∂g∂x_vθ_ng \ ∂g∂p_Png_Qng)

    Png_Qng_mis_k_vθ_ng_by_Png_Qng =
        -(∂g∂x_Png_Qng_mismatch_by_non_gens_vh_θh \
        ∂g∂p_Png_Qng_mismatch_by_Png_Qng )

    #
    # Pg_Qg_mis_k_vθ_g =
    #     (∂g∂x_vθ_g + ∂g∂x_vθ_ng *
    #     Png_Qng_mis_k_vθ_ng_by_vθ_g + ∂g∂x_id_iq *
    #     id_iq_mis_kg_id_iq_by_vθ_g)

    Pg_Qg_mis_k_vθ_g =
        (∂g∂x_Pg_Qg_mismatch_by_gens_vh_θh +
        ∂g∂x_Pg_Qg_mismatch_by_non_gens_vh_θh *
        Png_Qng_mis_k_vθ_ng_by_vθ_g +
        ∂g∂x_Pg_Qg_mismatch_by_id_iq *
        id_iq_mis_kg_id_iq_by_vθ_g )

    # Pg_Qg_mis_k_δ_ed_dash_eq_dash =
    #     (∂g∂p_δ_ed_dash_eq_dash + ∂g∂x_id_iq *
    #     id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash)

    Pg_Qg_mis_k_δ_ed_dash_eq_dash =
        (∂g∂p_Pg_Qg_mismatch_by_δ_ed_dash_eq_dash +
        ∂g∂x_Pg_Qg_mismatch_by_id_iq *
        id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash)

    # Pg_Qg_mis_k_Png_Qng =
    #     ∂g∂x_vθ_ng * Png_Qng_mis_k_vθ_ng_by_Png_Qng

    Pg_Qg_mis_k_Png_Qng =
        ∂g∂x_Pg_Qg_mismatch_by_non_gens_vh_θh *
        Png_Qng_mis_k_vθ_ng_by_Png_Qng


    #
    Pg_Qg_mis_vθ_g_by_δ_ed_dash_eq_dash =
        -(Pg_Qg_mis_k_vθ_g\Pg_Qg_mis_k_δ_ed_dash_eq_dash)

    Pg_Qg_mis_vθ_g_by_Png_Qng =
        -(Pg_Qg_mis_k_vθ_g\Pg_Qg_mis_k_Png_Qng)

    if loc_load_exist == true

        Pg_Qg_mis_vθ_g_by_Pgll_Qgll =
            -(Pg_Qg_mis_k_vθ_g\∂g∂p_Pgll_Qgll)
    end

    #
    kg_id_iq_by_δ_ed_dash_eq_dash =
        (id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash +
        id_iq_mis_kg_id_iq_by_vθ_g *
        Pg_Qg_mis_vθ_g_by_δ_ed_dash_eq_dash)

    kg_id_iq_by_Png_Qng =
        id_iq_mis_kg_id_iq_by_vθ_g *
        Pg_Qg_mis_vθ_g_by_Png_Qng

    if loc_load_exist == true

         kg_id_iq_by_Pgll_Qgll =
             id_iq_mis_kg_id_iq_by_vθ_g *
             Pg_Qg_mis_vθ_g_by_Pgll_Qgll
    end

    # 
    # kf_id_iq_by_δ_ed_dash_eq_dash =
    #     ∂f∂p_id_iq *
    #     id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash 

    kf_id_iq_by_δ_ed_dash_eq_dash =
        ∂f∂p_im_state_by_gens_id_iq *
        id_iq_mis_kg_id_iq_by_δ_ed_dash_eq_dash 

    # kf_id_iq_by_vθ_g =
    #     ∂f∂p_id_iq * id_iq_mis_kg_id_iq_by_vθ_g

    kf_id_iq_by_vθ_g =
        ∂f∂p_im_state_by_gens_id_iq *
        id_iq_mis_kg_id_iq_by_vθ_g

    #
    # kf_δ_ed_dash_eq_dash =
    #     ( kf_id_iq_by_δ_ed_dash_eq_dash +
    #     (kf_id_iq_by_vθ_g + ∂f∂p_vh_g_wt_θh0) *
    #     Pg_Qg_mis_vθ_g_by_δ_ed_dash_eq_dash ) 

    kf_δ_ed_dash_eq_dash =
        ( kf_id_iq_by_δ_ed_dash_eq_dash +
        (kf_id_iq_by_vθ_g +
        ∂f∂p_im_state_by_gens_vh_wt_zero_θh) *
        Pg_Qg_mis_vθ_g_by_δ_ed_dash_eq_dash ) 

    # kf_Png_Qng =
    #     (kf_id_iq_by_vθ_g + ∂f∂p_vh_g_wt_θh0 ) *
    #     Pg_Qg_mis_vθ_g_by_Png_Qng 

    kf_Png_Qng =
        (kf_id_iq_by_vθ_g +
        ∂f∂p_im_state_by_gens_vh_wt_zero_θh ) *
        Pg_Qg_mis_vθ_g_by_Png_Qng 

    if loc_load_exist == true

         kf_Pgll_Qgll =
             (kf_id_iq_by_vθ_g + ∂f∂p_vh_g_wt_θh0) *
             Pg_Qg_mis_vθ_g_by_Pgll_Qgll

    end

    #---------------------------------------------

    Asys = copy( ∂f∂x_im_state_by_im_state )

    Asys[:, flat_δ_ed_dash_eq_dash_Idxs_in_state] .=
        Asys[:, flat_δ_ed_dash_eq_dash_Idxs_in_state ] +
        kf_δ_ed_dash_eq_dash


    sys_eigvalues = eigvals(Asys)

    #---------------------------------------------

    Axp =
        Asys[im_pure_states_idxs,
             im_pure_states_idxs ]

    Axa =
        Asys[im_pure_states_idxs,
             im_τm_tilade_vf_states_idxs ]

    Aap =
        Asys[im_τm_tilade_vf_states_idxs,
             im_pure_states_idxs ]

    Aaa =
        Asys[im_τm_tilade_vf_states_idxs,
             im_τm_tilade_vf_states_idxs ]

    pure_Asys = Axp - Axa * (Aaa \ Aap )

    #---------------------------------------------

    pure_Asys_eigen_object = eigen( pure_Asys )

    pure_Asys_eig_values = pure_Asys_eigen_object.values

    printed_pure_Asys_eig_values =
        round.(pure_Asys_eig_values; digits=4)
    

    pure_Asys_eigvecs_right =
        pure_Asys_eigen_object.vectors

    # inv_pure_Asys_eigvecs_right =
    #     inv(pure_Asys_eigvecs_right)

    inv_pure_Asys_eigvecs_right =
        pure_Asys_eigvecs_right \
        LinearAlgebra.I(size(pure_Asys_eigvecs_right)[1] )
    
    #--------------------------------------------------
    # states associated with eig_values
    #--------------------------------------------------

    M_diag = inv_pure_Asys_eigvecs_right *
        pure_Asys * pure_Asys_eigvecs_right

    printed_M_diag =
        round.( M_diag; digits = 4 )

    # tup_state_var_eig_value =
    #     [(a_state_sym, an_eig_value)
    #      for ( a_state_sym, an_eig_value )  in
    #          zip( im_ode_sym, printed_eig_values )]

    #--------------------------------------------------
    # participation factor
    #--------------------------------------------------

    PF_pure_Asys  = get_participation_factors( pure_Asys )

    printed_PF_pure_Asys =
        round.( PF_pure_Asys; digits = 4 )

    #--------------------------------------------------

    eig_values, eigvecs_left, eigvecs_right =
        get_eigens( pure_Asys )

    #--------------------------------------------------

    im_pure_states_syms =
        im_ode_sym[im_pure_states_idxs]

    #--------------------------------------------------
    #--------------------------------------------------


    ∂f∂p_pg_qh0_1 = ∂f∂p_pg_qh0[
        im_pure_states_idxs, :]

    ∂f∂p_pg_qh0_2 = ∂f∂p_pg_qh0[
        im_τm_tilade_vf_states_idxs, :]

    ∂f∂p_pg_qh0_f =
        (∂f∂p_pg_qh0_1 - Axa *(Aaa\∂f∂p_pg_qh0_2) )


    kf_Png_Qng_1 = kf_Png_Qng[
        im_pure_states_idxs, :]

    kf_Png_Qng_2 = kf_Png_Qng[
        im_τm_tilade_vf_states_idxs, :]

    kf_Png_Qng_f = (kf_Png_Qng_1 -
        Axa *(Aaa\kf_Png_Qng_2) )

    ∂f∂p_ωs_ωref0_vref0_porder0_1 =
        ∂f∂p_ωs_ωref0_vref0_porder0[
            im_pure_states_idxs, :]

    ∂f∂p_ωs_ωref0_vref0_porder0_2 =
        ∂f∂p_ωs_ωref0_vref0_porder0[
            im_τm_tilade_vf_states_idxs, :]

    ∂f∂p_ωs_ωref0_vref0_porder0_f =
        (∂f∂p_ωs_ωref0_vref0_porder0_1 -
        Axa *(Aaa\∂f∂p_ωs_ωref0_vref0_porder0_2 ) )

    if loc_load_exist == true

        kf_Pgll_Qgll_1 = kf_Pgll_Qgll[
            im_pure_states_idxs, :]

        kf_Pgll_Qgll_2 = kf_Pgll_Qgll[
            im_τm_tilade_vf_states_idxs, :]

        kf_Pgll_Qgll_f =
            (kf_Pgll_Qgll_1 - Axa *(Aaa\kf_Pgll_Qgll_2) )
        
    end

    #--------------------------------------------------

    pure_states =
        sim_state_x0[ im_pure_states_idxs ]

    PQ_g = [ [ [ a_P, a_Q ] for ( a_P, a_Q ) in
                   zip(P_gens,  Q_gens) ]...;]

    PQ_ng = [ [ [ a_P, a_Q ] for ( a_P, a_Q ) in
                   zip(P_non_gens, Q_non_gens) ]...;]

    if loc_load_exist == true

        PQ_g_ll = [ [ [ a_P, a_Q ] for ( a_P, a_Q ) in
                       zip(P_g_loc_load, Q_g_loc_load) ]...;]
    end

    #--------------------------------------------------

    Δpure_states =
        similar( pure_states )

    Δx     = similar( pure_states )

    ΔPQ_g  = similar( PQ_g )

    ΔPQ_ng = similar( PQ_ng )

    Δωs_ωref0_vref0_porder0  =
        similar( flat_ωs_ωref0_vref0_porder0 )

    if loc_load_exist == true

        ΔPQll_g = similar( PQ_g_ll )
        
    end

    #--------------------------------------------------
    #--------------------------------------------------

    if loc_load_exist == true

        Δx = pure_Asys * Δpure_states +
            ∂f∂p_pg_qh0_f   * ΔPQ_g   +
            kf_Png_Qng_f  * ΔPQ_ng  +
            kf_Pgll_Qgll_f * ΔPQll_g +
            ∂f∂p_ωs_ωref0_vref0_porder0_f *
            Δωs_ωref0_vref0_porder0
        
    else
        
        Δx = pure_Asys    * Δpure_states +
            ∂f∂p_pg_qh0_f * ΔPQ_g   +
            kf_Png_Qng_f  * ΔPQ_ng  +
            ∂f∂p_ωs_ωref0_vref0_porder0_f *
            Δωs_ωref0_vref0_porder0
        
    end



    if loc_load_exist == true

        return (; pure_Asys, M_diag,
                PF_pure_Asys, eig_values,
                eigvecs_left, eigvecs_right,
                inv_pure_Asys_eigvecs_right,
                ∂f∂p_pg_qh0_f, kf_Png_Qng_f,
                ∂f∂p_ωs_ωref0_vref0_porder0_f,
                kf_Pgll_Qgll_f, im_pure_states_syms  )
        
    else
        
        return (; pure_Asys, M_diag,
                PF_pure_Asys, eig_values,
                eigvecs_left, eigvecs_right,
                inv_pure_Asys_eigvecs_right,
                ∂f∂p_pg_qh0_f, kf_Png_Qng_f,
                ∂f∂p_ωs_ωref0_vref0_porder0_f,
                im_pure_states_syms )
        
    end
    
end



function simulate_sd_dynamics_model_by_per_gen_func(
    flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para;
    kwd = sd_dynamics_model_by_per_gen_kwd  )

    
    (;
     sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para,
     sim_state_x0,
     sim_timespan,
     im_sym,
     im_mass_matrix,
     ode_alg,
     dt) = kwd

    #--------------------------------------------------
    #--------------------------------------------------
    # sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!
    #--------------------------------------------------
    #--------------------------------------------------
    
    sim_fun_kwd_para =
        sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para
    
    sim_fun_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para

    # flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para

    sim_func  =
        sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!

    #-----------------------------------------------------

    sd_per_gen_ode_sep_vtf_pf_sim_sol =
        DifferentialEquations.solve(
            ODEProblem(
                ODEFunction{true}(
                    (dx, x, p, t) ->
                        sim_func(dx, x, p, t;
                                 kwd_para =
                                     sim_fun_kwd_para
                         );
                    mass_matrix =
                        im_mass_matrix ,
                    syms =
                        im_sym ),
                sim_state_x0 ,
                sim_timespan,
                sim_fun_para ),
            ode_alg,
            dt=dt)        
    
    #-----------------------------------------------------
    # plot
    #-----------------------------------------------------

    (; gens_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    getproperty(
                        getproperty(
                            sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para,
                            :init_pf_by_vh_θh_by_parts_kwd),
                        :init_dyn_pf_mismatch_kwd_para ),
                    :dyn_pf_Idxs_kwd_para),
                :dyn_pf_fun_kwd_net_idxs ),
            (:gens_nodes_idx,
             :non_gens_nodes_idx,
             :all_nodes_idx))
       
    
    (;n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_all_nodes_idx) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    getproperty(
                        getproperty(
                            sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para,
                            :init_pf_by_vh_θh_by_parts_kwd),
                        :init_dyn_pf_mismatch_kwd_para ),
                    :dyn_pf_Idxs_kwd_para),
                :dyn_pf_fun_kwd_n2s_idxs ),
            (:n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_all_nodes_idx))
    
    gens_nodes_names =
        ["bus$(n2s_gens_idx[idx])"
         for idx in gens_nodes_idx ]
    
    plot_single_machine =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = sd_per_gen_ode_sep_vtf_pf_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_multiple_machine =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = sd_per_gen_ode_sep_vtf_pf_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = gens_nodes_names,
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            # vars = [:δ, :ω ],
            tspan = sim_timespan,
            fmt = :png)

    plot_indices_of_vars =
        get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = sd_per_gen_ode_sep_vtf_pf_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ])

    return(; sd_per_gen_ode_sep_vtf_pf_sim_sol,
           im_sym, gens_nodes_idx,
           non_gens_nodes_idx, all_nodes_idx,
           n2s_gens_idx, n2s_all_nodes_idx,
           n2s_non_gens_idx, plot_multiple_machine )

end


function get_sensitivity_and_stability_analysis_by_per_gen_func(
    init_flat_vh_flat_θh,
    init_dyn_pf_flat_para,
    flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh ,
    flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para;
    kwd = sensitivity_and_stability_analysis_kwd  )


    (; stability_ode_per_gen_wt_intra_dyn_pf_kwd,
     stability_sd_dynamics_wt_intra_dyn_pf_kwd)  = kwd

    #----------------------------------------------------
    #----------------------------------------------------
    
    stability_sd_dynamics_wt_intra_dyn_pf =
        get_small_signal_stability_sd_dynamics_wt_intra_dyn_pf(flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para; kwd = stability_sd_dynamics_wt_intra_dyn_pf_kwd )

    stability_ode_per_gen_wt_intra_dyn_pf =
        get_small_signal_stability_analysis_ode_per_gen_wt_intra_dyn_pf( init_flat_vh_flat_θh, init_dyn_pf_flat_para, flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh ; kwd = stability_ode_per_gen_wt_intra_dyn_pf_kwd )

    #-----------------------------------------------------

    sd_eig_values = round.(
        stability_sd_dynamics_wt_intra_dyn_pf.eig_values;digits=4)

    sd_PF_pure_Asys = round.(
        stability_sd_dynamics_wt_intra_dyn_pf.PF_pure_Asys;digits=4)

    sd_pure_Asys = round.(
        stability_sd_dynamics_wt_intra_dyn_pf.pure_Asys;digits=4)

    sd_im_pure_states_syms =
        stability_sd_dynamics_wt_intra_dyn_pf.im_pure_states_syms


    sd_eig_values_in_states_participation =
        get_eig_values_in_states_participation(
            sd_eig_values, sd_PF_pure_Asys,
            sd_im_pure_states_syms;
            participation_threshold = 0.2 )

    #-----------------------------------------------------


    ode_eig_values =
        round.(stability_ode_per_gen_wt_intra_dyn_pf.eig_values;digits=4)

    ode_PF_pure_Asys =
        round.(stability_ode_per_gen_wt_intra_dyn_pf.PF_pure_Asys;digits=4)

    ode_pure_Asys =
        round.(stability_ode_per_gen_wt_intra_dyn_pf.pure_Asys;digits=4)

    ode_im_pure_states_syms =
        stability_ode_per_gen_wt_intra_dyn_pf.im_pure_states_syms


    ode_eig_values_in_states_participation =
        get_eig_values_in_states_participation(
            ode_eig_values, ode_PF_pure_Asys,
            ode_im_pure_states_syms;
            participation_threshold = 0.2 )

    #-----------------------------------------------------

    return (; stability_sd_dynamics_wt_intra_dyn_pf,
            stability_ode_per_gen_wt_intra_dyn_pf,
            sd_eig_values_in_states_participation,
            ode_eig_values_in_states_participation )
    

end



function get_net_sensitivity_and_stability_analysis_by_per_gen_func(
    netd;
    abstol = 1e-14,
    reltol = 1e-14,    
    pf_alg  = NewtonRaphson(),
    ode_alg_2 = Rodas4(),
    # algr_name = "rodas4" 
    ode_alg = ImplicitMidpoint(),
    algr_name = "ImplicitMidpoint",
    dt = 0.01,
    sim_timespan    = (0.0, 10.0),
    only_gen         = false,
    no_control_device = false,
    all_wt_size_and_dims = false,
    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial
    Idxs_type        = :Idxs_im, 
    diagnostics_data = false,
    dynamics_by_gens = false,
    dynamics_by_per_gen = true,
    streamlined =  true)
    
    #-----------------------------------------------
    #-----------------------------------------------

    sensitivity_and_stability_para_Idxs_by_per_gen =
        get_sensitivity_and_stability_para_Idxs_by_per_gen_func(
        netd;
        abstol = abstol,
        reltol = reltol,    
        pf_alg  = pf_alg,
        ode_alg_2 = ode_alg_2,
        # algr_name = "rodas4" 
        ode_alg = ode_alg,
        algr_name = algr_name,
        dt = dt,
        sim_timespan = sim_timespan,
        only_gen = only_gen,
        no_control_device = no_control_device,
        all_wt_size_and_dims = all_wt_size_and_dims,
        # :Idxs_hybrid, :Idxs_im, :Idxs_industrial
        Idxs_type        = Idxs_type, 
        diagnostics_data = diagnostics_data,
        dynamics_by_gens = dynamics_by_gens,
        dynamics_by_per_gen = dynamics_by_per_gen,
        streamlined =  streamlined)

    #-----------------------------------------------------
    #-----------------------------------------------------

    stability_sd_dynamics_wt_intra_dyn_pf_kwd = sensitivity_and_stability_para_Idxs_by_per_gen.stability_sd_dynamics_wt_intra_dyn_pf_kwd

    
    #-----------------------------------------------------

    stability_ode_per_gen_wt_intra_dyn_pf_kwd = sensitivity_and_stability_para_Idxs_by_per_gen.stability_ode_per_gen_wt_intra_dyn_pf_kwd
    
    #-----------------------------------------------------
  
    vars =
        sensitivity_and_stability_para_Idxs_by_per_gen.vars

    init_flat_vh_flat_θh = vars.flat_vh_flat_θh_x0

    flat_ωs_ωref0_vref0_porder0 =
        vars.flat_ωs_ωref0_vref0_porder0

    init_dyn_pf_flat_para = vars.init_dyn_pf_flat_para

    gens_dynamic_id_iq_pg_vh = vars.gens_dynamic_id_iq_pg_vh
    
    flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para = [ init_flat_vh_flat_θh; flat_ωs_ωref0_vref0_porder0; init_dyn_pf_flat_para]
        
    flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh =
        [flat_ωs_ωref0_vref0_porder0;
         [gens_dynamic_id_iq_pg_vh...;] ]

    
    #-----------------------------------------------------
    
    stability_sd_dynamics_wt_intra_dyn_pf =
        get_small_signal_stability_sd_dynamics_wt_intra_dyn_pf(flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para ;kwd = stability_sd_dynamics_wt_intra_dyn_pf_kwd )

    stability_ode_per_gen_wt_intra_dyn_pf =
        get_small_signal_stability_analysis_ode_per_gen_wt_intra_dyn_pf( init_flat_vh_flat_θh, init_dyn_pf_flat_para, flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh ; kwd = stability_ode_per_gen_wt_intra_dyn_pf_kwd )

    #-----------------------------------------------------

    sd_eig_values = round.(
        stability_sd_dynamics_wt_intra_dyn_pf.eig_values;digits=4)

    sd_PF_pure_Asys = round.(
        stability_sd_dynamics_wt_intra_dyn_pf.PF_pure_Asys;digits=4)

    sd_pure_Asys = round.(
        stability_sd_dynamics_wt_intra_dyn_pf.pure_Asys;digits=4)

    sd_im_pure_states_syms =
        stability_sd_dynamics_wt_intra_dyn_pf.im_pure_states_syms


    sd_eig_values_in_states_participation =
        get_eig_values_in_states_participation(
            sd_eig_values, sd_PF_pure_Asys,
            sd_im_pure_states_syms;
            participation_threshold = 0.2 )

    #-----------------------------------------------------


    ode_eig_values =
        round.(stability_ode_per_gen_wt_intra_dyn_pf.eig_values;digits=4)

    ode_PF_pure_Asys =
        round.(stability_ode_per_gen_wt_intra_dyn_pf.PF_pure_Asys;digits=4)

    ode_pure_Asys =
        round.(stability_ode_per_gen_wt_intra_dyn_pf.pure_Asys;digits=4)

    ode_im_pure_states_syms =
        stability_ode_per_gen_wt_intra_dyn_pf.im_pure_states_syms


    ode_eig_values_in_states_participation =
        get_eig_values_in_states_participation(
            ode_eig_values, ode_PF_pure_Asys,
            ode_im_pure_states_syms;
            participation_threshold = 0.2 )

    #-----------------------------------------------------

    return (; stability_sd_dynamics_wt_intra_dyn_pf,
            stability_ode_per_gen_wt_intra_dyn_pf,
            sd_eig_values_in_states_participation,
            ode_eig_values_in_states_participation )
    

end

