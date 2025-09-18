# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247


function get_eig_values_in_states_participation(
    eig_values,
    PF_pure_Asys,
    im_pure_states_syms;
    participation_threshold = 0.2 )

    rows_size_PF_pure_Asys,cols_size_PF_pure_Asys =
        size(PF_pure_Asys)

    cols_PF_pure_Asys =
        [PF_pure_Asys[:,a_col_idx]
         for a_col_idx in
             1:cols_size_PF_pure_Asys]


    rows_PF_pure_Asys =
        [PF_pure_Asys[a_row_idx,:]
         for a_row_idx in
             1:rows_size_PF_pure_Asys]


    indices_of_PF_cols_greter_than =
        [ findall(>( participation_threshold  ),
                  PF_pure_Asys[:, a_col] )
          for a_col  in
              1:cols_size_PF_pure_Asys ]


    indices_of_PF_rows_greter_than =
        [ findall(>( participation_threshold ),
                  PF_pure_Asys[a_row, :] )
          for a_row  in
              1:rows_size_PF_pure_Asys ]


    tup_eig_value_state_vars_PF_cols_greter_than =
        [(eig_values[idx_eig_value],
          im_pure_states_syms[idx_row_greater_than],
          
         a_col_PF_pure_Asys[idx_row_greater_than] )
         for (idx_eig_value, idx_row_greater_than,
              a_col_PF_pure_Asys) in
             zip(1:cols_size_PF_pure_Asys,
                 indices_of_PF_cols_greter_than,
                 cols_PF_pure_Asys)]


    tup_state_vars_eig_value_PF_rows_greter_than =
        [(im_pure_states_syms[state_idx],
          eig_values[idx_col_greater_than],
         a_row_PF_pure_Asys[idx_col_greater_than] )
         for (state_idx, idx_col_greater_than,
              a_row_PF_pure_Asys) in
             zip(1:rows_size_PF_pure_Asys,
                 indices_of_PF_rows_greter_than,
                 rows_PF_pure_Asys)]

    return (
        ;tup_eig_value_state_vars_PF_cols_greter_than,
        tup_state_vars_eig_value_PF_rows_greter_than )
end

# ---------------------------------------------------


"""

https://arblib.org/acb_mat.html#acb-mat-eigenvalues

https://github.com/kalmarek/Arblib.jl

matrix_rows, matrix_cols = size(system_matrix)

eig_val = Arblib.AcbVector(zeros(matrix_rows))

Acb_system_matrix = Arblib.AcbMatrix(system_matrix)

eigvecs_left  = similar(Acb_system_matrix)

eigvecs_right = similar(Acb_system_matrix)


return_code = Arblib.approx_eig_qr!(
    eig_val,
    eigvecs_left,
    eigvecs_right,
    Acb_system_matrix,
    Mag(),
    0,
    Arblib._precision(Acb_system_matrix) )

"""

function get_eigens_via_arblib(
    system_matrix;
    prec = nothing  )


    matrix_rows, matrix_cols = size(system_matrix)

    eig_values = Arblib.AcbVector(zeros(matrix_rows))

    Acb_system_matrix = Arblib.AcbMatrix(system_matrix)

    eigvecs_left  = similar(Acb_system_matrix)

    eigvecs_right = similar(Acb_system_matrix)

    if prec ==  nothing
        prec = Arblib._precision(Acb_system_matrix)
    else
        prec =  prec
    end
    
        
    return_code =
        Arblib.approx_eig_qr!(
            eig_values,
            eigvecs_left,
            eigvecs_right,
            Acb_system_matrix,
            Mag(),
            0,
            prec )

    return (; eig_values, eigvecs_left,
            eigvecs_right, return_code )
    
end



"""
    get_eigens(system_matrix)

This function returns eigen values, left and right
eigen vectors.

"""
function get_eigens(system_matrix)

    # https://ralphas.github.io/GenericSchur.jl/stable/
    
    system_matrix = system_matrix .+ 0im
    
    schur_object = schur( system_matrix )

    eig_values = schur_object.values

    eigvecs_right = eigvecs( schur_object )

    eigvecs_left = eigvecs(schur_object,left=true)

    return (; eig_values, eigvecs_left, eigvecs_right )
    

end


"""
https://tobydriscoll.net/fnc-julia/linsys/norms.html

https://discourse.julialang.org/t/computing-left-eigenvectors-with-high-precision/72949

https://ralphas.github.io/GenericSchur.jl/stable/

`get_participation_factors`
This function returns the participation factor matrix

Sauer: see page 232, equation 8.84
    """

"""

    get_participation_factors(
        system_matrix )

Returns the participation factor matrix.

"""
function get_participation_factors(
    system_matrix )

    (; eig_values, eigvecs_left, eigvecs_right ) =
        get_eigens(system_matrix)


    participation_factor =
        abs.(eigvecs_left) .* abs.(eigvecs_right)

    for i in eachindex(eig_values)

        λ_i_pf = participation_factor[:,i] ./
            (abs.(eigvecs_left[:,i]') *
            abs.(eigvecs_right[:,i] ))
        
        max_λ_i_pf = max(abs.(λ_i_pf)... )
        
        participation_factor[:,i] .= λ_i_pf ./ max_λ_i_pf
            
        # max_pf = max(participation_factor[:,i]... )

        # participation_factor[:,i] .=
        #     participation_factor[:,i] ./  max_pf 

    end

    return participation_factor
    
end



#####################################################
# ---------------------------------------------------
# utility linearisation stability
# ---------------------------------------------------
#####################################################


"""
    get_generic_stability_static_powerflow(
        pf_PQ_param ;
        kwd_para =
            stability_static_powerflow_kwd_para)


Returns a tuple of `(vh, θh, pf_P_gens, pf_Q_gens)` based on powerflow.

"""
function get_generic_stability_static_powerflow(
    pf_PQ_param ;
    kwd_para =
        stability_static_powerflow_kwd_para)

    (;
     pf_alg,
     loc_load_exist,
     
     ode_gens_generic_para,
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     pf_sta_ΔPQ_mismatch_parameters,    
     Ybr_cal_and_edges_orientation,
     baseMVA,
     basekV) =
         kwd_para
    
    #----------------------------------------    
    # Static power flow
    #----------------------------------------

    (dyn_P_gens_Idxs,
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

        
    (pf_kw_para,
     red_types_Idxs_etc) =
         NamedTupleTools.select(
             pf_sta_ΔPQ_mismatch_parameters,
             (:pf_kw_para,
              :red_types_Idxs_etc ) )

    (red_vh_Idxs,
     red_θh_Idxs) =
         NamedTupleTools.select(
             red_types_Idxs_etc,
             (:red_vh_Idxs,
              :red_θh_Idxs) )

    #----------------------------------------

    P_gens =
        pf_PQ_param[
            dyn_P_gens_Idxs]
    
    Q_gens =
        pf_PQ_param[
            dyn_Q_gens_Idxs]
    
    P_non_gens =
        pf_PQ_param[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        pf_PQ_param[
            dyn_Q_non_gens_Idxs]
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            pf_PQ_param[
                dyn_P_gens_loc_load_Idxs]
        
    else
        
        P_g_loc_load = []
        
        Q_g_loc_load = []
        
    end
    

    sta_pf_PQ_para =
        (;P_gens,
         Q_gens,
         P_non_gens,
         Q_non_gens,
         P_g_loc_load,
         Q_g_loc_load,
         loc_load_exist)
    
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
        (;Ybr_cal_and_edges_orientation,
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
    
    (vh,
     θh,
     pf_P_gens,
     pf_Q_gens) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:vh,
             :θh,
             :pf_P_gens,
             :pf_Q_gens))
    
    return (vh, θh, pf_P_gens, pf_Q_gens)

end



"""
    get_generic_electro_mechanical_oscillation_indicies(
       pf_P_gens, pf_Q_gens,
       vh, gens_vh, gens_θh, ωs;
       kwd_para =
           inm_electro_mech_oscill_kwd_para )


Returns a nampedruple of `(;Yint, Yred, δg, Eg, Cinm, Dinm, Aω_matrix, A_matrix)`.
"""
function get_generic_electro_mechanical_oscillation_indicies(
    pf_P_gens, pf_Q_gens,
    vh,
    # θh,
    gens_vh, gens_θh, ωs;
    kwd_para =
        inm_electro_mech_oscill_kwd_para )

    (; generic_gens_para,
     pf_sta_ΔPQ_mismatch_parameters ) =
         kwd_para

    #--------------------------------------

    (pf_kw_para,
     pf_PQ_param) =
         NamedTupleTools.select(
             pf_sta_ΔPQ_mismatch_parameters,
             (:pf_kw_para,
              :pf_PQ_param) )
    
    #--------------------------------------
    
    (;H,
     X_d_dash,
     X_q_dash) =
        NamedTupleTools.select(
            get_selected_vec_nt_to_vec_vec(
                generic_gens_para,
                nothing;
                selections =
                    (:H,
                     :X_d_dash,
                     :X_q_dash ),
                vec_datatype = Float64 ),
            (:H,
             :X_d_dash,
             :X_q_dash ) )

    #----------------------------------------

    y_aug_kw_para =
        (; X_d_dash, pf_kw_para )
    
    (Yred,
     Yint) =
        NamedTupleTools.select(
            get_Y_aug_matrices(
                pf_PQ_param,
                vh;
                y_aug_kw_para =
                    y_aug_kw_para ),
            (:Yred,
             :Y_internal_nodes ))

    #----------------------------------------

    
    (;Eg,
     δg) =
       NamedTupleTools.select(
            get_selected_vec_nt_to_vec_vec(
                get_a_gen_Eg_wt_δg.(
                    gens_vh, gens_θh,
                    pf_P_gens, pf_Q_gens,
                    X_d_dash  ),
                nothing;
                selections =
                    (:Eg,
                     :δg ),
                vec_datatype = Float64 ),
            (:Eg,
             :δg ) )

        
    #----------------------------------------
    #----------------------------------------
    
    Cinm, Dinm = get_Cinm_Dinm( Eg, Yint)

    #----------------------------------------

    Aω_matrix = get_Aω_matrix(
        δg, Eg, H, Yint, ωs, nothing)
    

    (; Aω_matrix, A_matrix) =
        get_electro_mechanical_mode_matrices(
            δg, Eg, H, Yint, ωs )

    return (;Yint, Yred,
            δg, Eg,
            Cinm, Dinm,
            Aω_matrix, A_matrix)
    
end


function get_generic_reduced_and_linearised_model_parameters(
    net_data_by_components_file,    
    # baseMVA,
    basekV,
    use_pu_in_PQ,
    line_data_in_pu,
    pf_alg;
    components_libs_dir = "",
    tt = 0.0,
    in_components_type_sym =
        false)
    
    #----------------------------------------
    
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

    #--------------------------------------

    (;
     pf_sta_ΔPQ_mismatch_parameters,
     sta_pf_PQ_para,
     baseMVA,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     Ybr_cal_and_edges_orientation,
     
     ode_gens_generic_para ,
     ode_gens_para,

     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,

     dyn_pf_mismatch_vars_kwd_para ) =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                     net_data_by_components_file;
            
                    basekV = basekV,    
                    use_pu_in_PQ = use_pu_in_PQ,
                    line_data_in_pu = line_data_in_pu,
                     
                    in_components_type_sym =
                        in_components_type_sym),
                 :net_generic_parameters),
             (:pf_sta_ΔPQ_mismatch_parameters,
              :sta_pf_PQ_para,
              :baseMVA,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :Ybr_cal_and_edges_orientation,

              :ode_gens_generic_para ,
              :ode_gens_para,

              :generic_gens_para,
              :generic_avrs_para,
              :generic_govs_para,

              :dyn_pf_mismatch_vars_kwd_para) )
    
    #----------------------------------------    

    (loc_load_exist,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs,

     dyn_pf_δ_eq_dash_Idx,
     dyn_pf_Png_Qng_Pll_Qll_Idx,

     dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,

     dyn_pf_idq_Idx,
     dyn_pf_idq_δ_ed_dash_non_gen_PQ_Idx,
     
     dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
          
     id_iq_pg_vh_Idx,
     ωs_ωref0_vref0_porder0_Idx,
     
     ωref0_vref0_porder0_Idx,
     ωref0_vref0_porder0_id_iq_vh_Idx,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx) =
         NamedTupleTools.select(
             getproperty(
                 get_net_generic_parameters_and_idx(
                    net_data_by_components_file;
                    in_components_type_sym =
                        false),
                 :net_generic_idx),
             (:loc_load_exist,
              :Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              :non_pre_ordered_pf_vars_Idxs,

              :dyn_pf_δ_eq_dash_Idx,
              :dyn_pf_Png_Qng_Pll_Qll_Idx,

              :dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,

              :dyn_pf_idq_Idx,
              :dyn_pf_idq_δ_ed_dash_non_gen_PQ_Idx,
              
              :dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
              
              :id_iq_pg_vh_Idx,
              :ωs_ωref0_vref0_porder0_Idx,

              :ωref0_vref0_porder0_Idx,
              :ωref0_vref0_porder0_id_iq_vh_Idx,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx))

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
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
    

    #----------------------------------------
    #----------------------------------------
    
    (;
     state_vars_idx,
     vec_comp_states_Idx,

     plants_states_syms ,
     generic_state_sym, 

     state_labels,
     algebraic_vars_labels,
     network_vars_labels,
     
     generic_model_states_comp_idxs_in_Idx ,
     generic_model_vars_wt_i_dq_Idx_in_state ,

     comps_callback_paras_funs,     
     comps_init_funs,
     comps_output_funs,
     comps_dyn_funs,
     ode_comps_dyn_funs,
     dae_comps_dyn_funs,

     algebraic_state_sym,
     model_syms,

     nodes_names,
     gens_nodes_names,
     non_gens_nodes_names,

     SM_gens_nodes_names,
     SC_gens_nodes_names,
            
     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,
     
     model_mass_matrix,     
     model_bool_dae_vars,

     ode_gens_mass_matrix,
     ode_gens_bool_dae_vars,

     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Png_Qng_Pll_Qll_Idx,

     dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     dyn_vh_id_iq_V_ref_Tm_Idx,
     dyn_V_ref_Tm_id_iq_vh_Idx,

     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
     dyn_vh_id_iq_ωref0_vref0_porder0_Idx,

     dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_3_gens_type_paras_Idx,
     dyn_5_gens_type_paras_Idx,

     dyn_2_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
     dyn_3_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
                 
     system_states_idx_kwd_para,
     system_paras_idx_kwd_para,
     plants_dyn_fun_idx_kwd_para,
     plants_algeb_fun_idx_kwd_para ) =
         NamedTupleTools.select(
             get_states_Idx_syms_wt_functions(
                 net_data_by_components_file,
                 dyn_pf_fun_kwd_net_idxs,
                 dyn_pf_fun_kwd_n2s_idxs;
                 components_libs_dir =
                     components_libs_dir),
             (:state_vars_idx,
              :vec_comp_states_Idx,

              :plants_states_syms,
              :generic_state_sym, 

              :state_labels,
              :algebraic_vars_labels,
              :network_vars_labels,

              :generic_model_states_comp_idxs_in_Idx ,
              :generic_model_vars_wt_i_dq_Idx_in_state ,

              :comps_callback_paras_funs,     
              :comps_init_funs,
              :comps_output_funs,
              :comps_dyn_funs,
              :ode_comps_dyn_funs,
              :dae_comps_dyn_funs,

              :algebraic_state_sym,
              :model_syms,

              :nodes_names,
              :gens_nodes_names,
              :non_gens_nodes_names,

              :SM_gens_nodes_names,
              :SC_gens_nodes_names,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :model_mass_matrix,     
              :model_bool_dae_vars,

              :ode_gens_mass_matrix,
              :ode_gens_bool_dae_vars,

              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Png_Qng_Pll_Qll_Idx,

              :dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :dyn_vh_id_iq_V_ref_Tm_Idx,
              :dyn_V_ref_Tm_id_iq_vh_Idx,

              :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
              :dyn_vh_id_iq_ωref0_vref0_porder0_Idx,

              :dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_3_gens_type_paras_Idx,
              :dyn_5_gens_type_paras_Idx,

              :dyn_2_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
              :dyn_3_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,

              :system_states_idx_kwd_para,
              :system_paras_idx_kwd_para,
              :plants_dyn_fun_idx_kwd_para,
              :plants_algeb_fun_idx_kwd_para))
    
        
    #----------------------------------------    
    # Static power flow
    #----------------------------------------

    (pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
         NamedTupleTools.select(
             pf_sta_ΔPQ_mismatch_parameters,
             (:pf_kw_para,
              :red_types_Idxs_etc,
              :pf_PQ_param) )

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
        (;Ybr_cal_and_edges_orientation,
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
    
    (pf_P_gens,
     pf_Q_gens,
     vh,
     θh,
     gens_vh,
     gens_θh) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens,
             :pf_Q_gens,
             :vh,
             :θh,
             :gens_vh,
             :gens_θh ) )

    
    #----------------------------------------
    # System load para
    #----------------------------------------

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

    #----------------------------------------

     pf_generic_gens_para =
        NamedTupleTools.select(
            get_selected_vec_nt_to_vec_vec(
                generic_gens_para,
                nothing;
                selections =
                    (:ra,
                     :X_d,
                     :X_q,     
                     :X_d_dash,
                     :X_q_dash ),
                vec_datatype = Float64 ),
            (:ra,
              :X_d,
              :X_q,     
              :X_d_dash,
             :X_q_dash ) )
    
    #----------------------------------------
    # callbacks
    #----------------------------------------

    plants_cb_paras_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para,
         comps_callback_paras_funs )

    #----------------------------------------
    
    """
    generic_model_callback_paras =
        plants_generic_model_callback_paras_func(
            state_vars_idx,
            plants_states_syms;
            kwd_para =
                plants_cb_paras_kwd_para )

    plants_cb_paras_switches =
        getproperty(
            generic_model_callback_paras,
            :plants_cb_paras_switches)
 
    
    avrs_govs_cb_sw =
        getproperty(
            generic_model_callback_paras,
            :plants_avr_gov_cb_para_sw_in_plant)

    avrs_govs_cb_sw_Idx =
        getproperty(
            generic_model_callback_paras,
            :plants_avr_gov_cb_para_sw_idx_in_plant )

    cb = cb_fun_make_state_callbacks(
        generic_model_callback_paras)
    """
    
    #----------------------------------------

    (plants_cb_paras_switches,
     list_selected_plants_state_event_cb_paras,
     list_selected_plants_state_affect_cb_paras,

     avrs_govs_cb_sw,
     avrs_govs_cb_sw_Idx ) =
         NamedTupleTools.select(
             plants_generic_model_callback_paras_func(
                 state_vars_idx,
                 plants_states_syms;
                 kwd_para =
                     plants_cb_paras_kwd_para ) ,
             (:plants_cb_paras_switches,
              :list_selected_plants_state_event_cb_paras,
              :list_selected_plants_state_affect_cb_paras,
              
              :plants_avr_gov_cb_para_sw_in_plant,
              :plants_avr_gov_cb_para_sw_idx_in_plant))
    
    cb = cb_fun_make_state_callbacks(
        list_selected_plants_state_event_cb_paras,
        list_selected_plants_state_affect_cb_paras )


    #----------------------------------------    
    # algebraic mismatch equation kwd_para
    #----------------------------------------    

    algebraic_generic_model_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes )    
    
    #----------------------------------------
    # generic model_dynamics_kwd_para
    #----------------------------------------
    
    ode_plants_kwd_para =
        (;state_vars_idx,

         ωref0_vref0_porder0_id_iq_vh_Idx,
         
         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         vec_comp_states_Idx,
         
         avrs_govs_cb_sw,
         avrs_govs_cb_sw_Idx,
         
         ode_comps_dyn_funs,
         comps_output_funs,
         ωs)

    
    #----------------------------------------
    #----------------------------------------
    # Dynamic simulation functions tests
    #----------------------------------------

    #----------------------------------------
    #########################################
    #----------------------------------------


    stability_static_powerflow_kwd_para =
        (;
         pf_alg,
         loc_load_exist,
     
         ode_gens_generic_para,
         Pg_Qg_Png_Qng_Pll_Qll_Idx,
         pf_sta_ΔPQ_mismatch_parameters,    
         Ybr_cal_and_edges_orientation,
         baseMVA,
         basekV)
    
    generic_stability_static_powerflow =
        get_generic_stability_static_powerflow(
            pf_PQ_param ;
            kwd_para =
                stability_static_powerflow_kwd_para)


    #----------------------------------------
    
    generic_network_init_kwd_para =
        (;gens_nodes_idx,
         gens_state_vars_idx_in_state,
         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,
         comps_init_funs,
         ωs)

    generic_system_dynamics_init_and_refs =
        get_generic_system_dynamics_init_and_refs(
        vh,
        θh,
        pf_P_gens,
        pf_Q_gens;
        kwd_para =
            generic_network_init_kwd_para)

    #----------------------------------------

    inm_electro_mech_oscill_kwd_para =
        (;generic_gens_para,
         pf_sta_ΔPQ_mismatch_parameters )

    
    generic_electro_mechanical_oscillation_indicies =
        get_generic_electro_mechanical_oscillation_indicies(
            pf_P_gens, pf_Q_gens,
            vh,
            gens_vh, gens_θh, ωs;
            kwd_para =
                inm_electro_mech_oscill_kwd_para )

    (Yint,
     Yred,
     δg,
     Eg,
     Cinm,
     Dinm,
     Aω_matrix,
     A_matrix) =
         NamedTupleTools.select(
             generic_electro_mechanical_oscillation_indicies,
             (:Yint,
              :Yred,
              :δg,
              :Eg,
              :Cinm,
              :Dinm,
              :Aω_matrix,
              :A_matrix ))
    
    #----------------------------------------

    (system_init,
     gens_δ,
     ed_dash,
     eq_dash,

     ω_ref,
     v_ref,
     p_order,

     gens_i_d,
     gens_i_q,
     gens_vh) =
         NamedTupleTools.select(
             generic_system_dynamics_init_and_refs,
             (:system_init,
              :δ,
              :ed_dash,
              :eq_dash,

              :ω_ref,
              :v_ref,
              :p_order,

              :gens_i_d,
              :gens_i_q,
              :gens_vh))

    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll =
        [ω_ref;v_ref;p_order;
         P_non_gens; Q_non_gens; 
         P_g_loc_load; Q_g_loc_load]


    # tt = 0.0
    
    linearised_dynamic_kwd_para =
        (;loc_load_exist,
      
         pf_generic_gens_para,         
         Ynet_wt_nodes_idx_wt_adjacent_nodes,
         
         state_vars_and_i_dq_Idx_in_state,
         state_algebraic_vars_Idx_in_state,
         gens_state_vars_idx_in_state,

         dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
         
         dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
         
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_n2s_idxs,

         ode_plants_kwd_para)


    return (;system_init,
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
            linearised_dynamic_kwd_para,
            generic_electro_mechanical_oscillation_indicies)
end



function get_generic_Asys_linearised_dynamic_model(
    sol,
    tt,
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll;
    use_init_condition =
        false,
    kwd_para =
        linearised_dynamic_kwd_para )

    if use_init_condition == true

        sol_t = sol
    else
        sol_t = sol(tt)
    end
    

    (;loc_load_exist,

     pf_generic_gens_para,         
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,
     gens_state_vars_idx_in_state,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
     
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_n2s_idxs,

     ode_plants_kwd_para) =
         kwd_para

    #----------------------------------------
    
    algebraic_generic_model_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes )    

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------    
    
    (f_ωref0_Idx,
     f_vref0_Idx,
     f_porder0_Idx,
     f_id_Idx,
     f_iq_Idx,
     f_vg_Idx ) =
         NamedTupleTools.select(
             dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
             (:dyn_ωref0_Idx,
              :dyn_vref0_Idx,
              :dyn_porder0_Idx,
              :dyn_id_Idx,
              :dyn_iq_Idx,
              :dyn_vh_Idx))
    
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

    (slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))

    (;n2s_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    

    #----------------------------------------

    N_g = length(gens_nodes_idx)

    #----------------------------------------    

    f_u_idx = first(f_ωref0_Idx):last(f_porder0_Idx)
    
    f_i_dq_idx = first(f_id_Idx):last(f_iq_Idx)

    # f_vg_Idx
    
    #----------------------------------------    

    g_vθ_idx =
        first(dyn_pf_vh_Idxs):last(dyn_pf_θh_Idxs)

    g_i_dq_idx =
        first(dyn_pf_id_Idxs):last(dyn_pf_iq_Idxs)

    #----------------------------------------    

    g_x_δ_edash_idx =
        first(dyn_δ_Idxs):last(dyn_eq_dash_Idxs)

    g_PQl_idx = loc_load_exist == true ? 
        (first(dyn_P_non_gens_Idxs):last(
            dyn_Q_gens_loc_load_Idxs)) :
        (first(dyn_P_non_gens_Idxs):last(
            dyn_Q_non_gens_Idxs))
    
    #----------------------------------------    

    δ_ed_dash_eq_dash_idx_in_state =
        [δ_idx_in_state;
         ed_dash_idx_in_state;
         eq_dash_idx_in_state ]

    #----------------------------------------    

    dyn_pf_gens_vh_Idxs =
        dyn_pf_vh_Idxs[gens_nodes_idx]

    dyn_pf_non_gens_vh_Idxs =
        dyn_pf_vh_Idxs[non_gens_nodes_idx]

    dyn_pf_slack_θh_Idxs =
        dyn_pf_θh_Idxs[slack_gens_nodes_idx]

    dyn_pf_non_slack_θh_Idxs =
        setdiff(collect(dyn_pf_θh_Idxs),
                dyn_pf_θh_Idxs[slack_gens_nodes_idx] )
    
    JLF_idx = [dyn_pf_non_gens_vh_Idxs ;
               dyn_pf_non_slack_θh_Idxs]
    
    algeb_by_state_link_cols_idx =
        [dyn_pf_gens_vh_Idxs;
         # dyn_pf_slack_θh_Idxs;
         dyn_pf_id_Idxs;
         dyn_pf_iq_Idxs]

    # δ_ed_dash_eq_dash_idx_in_state
    state_by_algeb_link_cols_idx =
        [δ_idx_in_state;
         ed_dash_idx_in_state;
         eq_dash_idx_in_state]
    
    #----------------------------------------    
    #----------------------------------------

    vh  =  sol_t[ vh_Idx_in_state]
    
    θh  =  sol_t[ θh_Idx_in_state]
    
    gens_i_d =  sol_t[id_Idx_in_state]
    
    gens_i_q =  sol_t[iq_Idx_in_state]


    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]
    
    #----------------------------------------
        
    gens_δ = sol_t[
        δ_idx_in_state]
        
    gens_ed_dash = sol_t[
        ed_dash_idx_in_state]
    
    gens_eq_dash = sol_t[
        eq_dash_idx_in_state]    

    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_v_ref_Idx]
    
    p_order =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_p_order_Idx]
    
    P_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Png_Idx]
    
    Q_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    #----------------------------------------
    
    ωref_vref_porder0_id_iq_vh =
        Float64[ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]


    #----------------------------------------    

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        Float64[gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------    

    f_x  = sol_t[
        state_var_Idx_in_state]
    
    f_dx = similar(f_x)

    g_y = sol_t[
        algebraic_var_Idx_in_state]

    g_dy = similar(g_y)

    #----------------------------------------
    
    ∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> ode_gens_plants_generic_model_func!(
            dx, x,
            ωref_vref_porder0_id_iq_vh,
            tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, f_x )

    
    ∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> ode_gens_plants_generic_model_func!(
            dx, f_x, p, tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, ωref_vref_porder0_id_iq_vh )
    
    #----------------------------------------
    
    ∂g∂y = ForwardDiff.jacobian(
        (dy, y ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, y,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                    algebraic_generic_model_kwd_para ),
        g_dy, g_y )

    
    ∂g∂p = ForwardDiff.jacobian(
        (dy, p ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, g_y, p;
            kwd_para = algebraic_generic_model_kwd_para
                     ),
        g_dy, δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para )

    #----------------------------------------
    
    Asys = copy(∂f∂x)
    
    Q = ∂f∂p[state_var_Idx_in_state, f_u_idx ]
    M = ∂f∂p[state_var_Idx_in_state, f_i_dq_idx ]
    N = ∂f∂p[state_var_Idx_in_state, f_vg_Idx ]
    
    #----------------------------------------

    A = ∂g∂y[g_vθ_idx, g_vθ_idx ]
    B = ∂g∂y[g_vθ_idx, g_i_dq_idx ]
    C = ∂g∂y[g_i_dq_idx, g_vθ_idx ]
    D = ∂g∂y[g_i_dq_idx, g_i_dq_idx ]

    E = ∂g∂p[g_vθ_idx, g_x_δ_edash_idx]
    F = ∂g∂p[g_vθ_idx, g_PQl_idx]
    G = ∂g∂p[g_i_dq_idx, g_x_δ_edash_idx]
    H = ∂g∂p[g_i_dq_idx, g_PQl_idx]

    #----------------------------------------

    I = A - B * (D\C)
    J = E - B * (D\G)
    K = F - B * (D\H)

    α = I\J
    β = I\K

    #----------------------------------------

    R = M * (D \ C) * α - M * (D \ G)
    S = M * (D \ C) * β - M * (D \ H)
    
    γ = α[gens_nodes_idx, :]
    ρ = β[gens_nodes_idx, :]


    ψ = R - N * γ
    ϕ = S - N * ρ


    Asys[:, δ_ed_dash_eq_dash_idx_in_state ] .=
        Asys[:, δ_ed_dash_eq_dash_idx_in_state ] + ψ

    # Δ⨰ = Asys Δx +  ϕ ΔPQng_PQll + Q Δu
    
    # where Δu = Δωref0_vref0_porder0
    
    #----------------------------------------
    #----------------------------------------

    size_states = length(state_var_Idx_in_state)
    size_algebr = length(algebraic_var_Idx_in_state)

    JAE = ∂g∂y

    JLF = JAE[JLF_idx, JLF_idx ]
    
    
    Γ = zeros(size_states, size_algebr )
    Λ = zeros(size_algebr, size_states )

    Γ[state_var_Idx_in_state,
           algeb_by_state_link_cols_idx] .=
               hcat(∂f∂p[state_var_Idx_in_state, f_vg_Idx ],
                    ∂f∂p[state_var_Idx_in_state,
                         f_i_dq_idx ] )
    
    # Λ[algebraic_var_Idx_in_state,
    #        state_by_algeb_link_cols_idx] .=
    #            ∂g∂p[:, g_x_δ_edash_idx]

    
    Λ[:, state_by_algeb_link_cols_idx] .=
               ∂g∂p[:, g_x_δ_edash_idx]
    
    # Asys = A - B * ( JAE \ C)
    
    Asys_dash = ∂f∂x - Γ * ( ∂g∂y \ Λ )
    
    #----------------------------------------
    #----------------------------------------
    
    sys_eigvalues = eigvals(Asys)

    printed_sys_eigvalues =
        round.( sys_eigvalues; digits = 4 )

    
    sys_eigvalues_dash  = eigvals(Asys_dash )

    printed_sys_eigvalues_dash  =
        round.( sys_eigvalues_dash ; digits = 4 )
    
    return (; Asys, ϕ, Q, JAE, JLF, Asys_dash,
            sys_eigvalues,
            printed_sys_eigvalues,
            sys_eigvalues_dash,
            printed_sys_eigvalues_dash )
    
end



"""
Seems not to be providing correct answer
"""
function test_get_generic_linearised_dynamic_model(
    sol,
    tt,
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
    nothing;
    use_init_condition =
        false,
    kwd_para =
        linearised_dynamic_kwd_para )

    if use_init_condition == true

        sol_t = sol
    else
        sol_t = sol(tt)
    end

    (;loc_load_exist,

     pf_generic_gens_para,         
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,
     gens_state_vars_idx_in_state,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
     
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_n2s_idxs,

     ode_plants_kwd_para) =
         kwd_para

    #----------------------------------------
    
    algebraic_generic_model_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes )    

    #----------------------------------------
    
    # (;loc_load_exist,
     
    #  state_vars_and_i_dq_Idx_in_state,
    #  state_algebraic_vars_Idx_in_state,
    #  gens_state_vars_idx_in_state,
     
    #  dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
    #  dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
    #  dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
    #  dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
    #  dyn_pf_fun_kwd_net_idxs,
    #  dyn_pf_fun_kwd_n2s_idxs,

    #  ode_plants_kwd_para) =
    #      kwd_para

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------    
    
    (f_ωref0_Idx,
     f_vref0_Idx,
     f_porder0_Idx,
     f_id_Idx,
     f_iq_Idx,
     f_vg_Idx ) =
         NamedTupleTools.select(
             dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
             (:dyn_ωref0_Idx,
              :dyn_vref0_Idx,
              :dyn_porder0_Idx,
              :dyn_id_Idx,
              :dyn_iq_Idx,
              :dyn_vh_Idx))
    
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

    (slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))

    (;n2s_slack_gens_idx,
     n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------        
    #----------------------------------------    

    f_u_idx = first(f_ωref0_Idx):last(f_porder0_Idx)
    
    f_i_dq_idx = first(f_id_Idx):last(f_iq_Idx)

    
    #----------------------------------------    

    g_vθ_idx =
        first(dyn_pf_vh_Idxs):last(dyn_pf_θh_Idxs)

    g_i_dq_idx =
        first(dyn_pf_id_Idxs):last(dyn_pf_iq_Idxs)
    
    vθ_g_idx  =
        [dyn_pf_vh_Idxs[gens_nodes_idx];
         dyn_pf_θh_Idxs[gens_nodes_idx]]
    
    vθ_ng_idx =
        [dyn_pf_vh_Idxs[non_gens_nodes_idx];
         dyn_pf_θh_Idxs[non_gens_nodes_idx]]

    ###
    
    g_rows_P_idx = dyn_pf_vh_Idxs
    g_rows_Q_idx = dyn_pf_θh_Idxs
    
    g_rows_P_g_idx = g_rows_P_idx[gens_nodes_idx]
    g_rows_Q_g_idx = g_rows_Q_idx[gens_nodes_idx]
    
    g_rows_P_ng_idx = g_rows_P_idx[non_gens_nodes_idx]
    g_rows_Q_ng_idx = g_rows_Q_idx[non_gens_nodes_idx]

    g_rows_idq_idx =
        first(dyn_pf_id_Idxs):last(dyn_pf_iq_Idxs)
    
    #----------------------------------------    

    g_x_δ_edash_idx =
        first(dyn_δ_Idxs):last(dyn_eq_dash_Idxs)
    
    g_PQ_ng_idx = 
        first(dyn_P_non_gens_Idxs):last(
            dyn_Q_non_gens_Idxs)

    g_PQll_idx = loc_load_exist == true ? 
        (first(dyn_P_gens_loc_load_Idxs):last(
            dyn_Q_gens_loc_load_Idxs)) : []

    #----------------------------------------    

    δ_ed_dash_eq_dash_idx_in_state =
        [δ_idx_in_state;
         ed_dash_idx_in_state;
         eq_dash_idx_in_state ]

    #----------------------------------------    

    # dyn_pf_vh_Idxs
    # dyn_pf_θh_Idxs
    # dyn_pf_id_Idxs
    # dyn_pf_iq_Idxs

    # slack_gens_nodes_idx
    
    # gens_nodes_idx
    
    # non_gens_nodes_idx

    dyn_pf_gens_vh_Idxs =
        dyn_pf_vh_Idxs[gens_nodes_idx]

    dyn_pf_non_gens_vh_Idxs =
        dyn_pf_vh_Idxs[non_gens_nodes_idx]

    dyn_pf_slack_θh_Idxs =
        dyn_pf_θh_Idxs[slack_gens_nodes_idx]

    dyn_pf_non_slack_θh_Idxs =
        setdiff(collect(dyn_pf_θh_Idxs),
                dyn_pf_θh_Idxs[slack_gens_nodes_idx] )
    
    JLF_idx = [dyn_pf_non_gens_vh_Idxs ;
               dyn_pf_non_slack_θh_Idxs]
    
    algeb_by_state_link_cols_idx =
        [dyn_pf_gens_vh_Idxs;
         # dyn_pf_slack_θh_Idxs;
         dyn_pf_id_Idxs;
         dyn_pf_iq_Idxs]

    # δ_ed_dash_eq_dash_idx_in_state
    state_by_algeb_link_cols_idx =
        [δ_idx_in_state;
         ed_dash_idx_in_state;
         eq_dash_idx_in_state]
    
    #----------------------------------------    
    #----------------------------------------

    vh  =  sol_t[ vh_Idx_in_state]
    
    θh  =  sol_t[ θh_Idx_in_state]
    
    gens_i_d =  sol_t[id_Idx_in_state]
    
    gens_i_q =  sol_t[iq_Idx_in_state]


    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]
    
    #----------------------------------------
        
    gens_δ = sol_t[
        δ_idx_in_state]
        
    gens_ed_dash = sol_t[
        ed_dash_idx_in_state]
    
    gens_eq_dash = sol_t[
        eq_dash_idx_in_state]    

    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_v_ref_Idx]
    
    p_order =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_p_order_Idx]
    
    P_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Png_Idx]
    
    Q_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    
    N_g  = length(gens_nodes_idx)
    N_ng = length(non_gens_nodes_idx)
    
    #----------------------------------------
    
    ωref_vref_porder0_id_iq_vh =
        Float64[ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]


    #----------------------------------------    

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        Float64[gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------    

    f_x  = sol_t[
        state_var_Idx_in_state]
    
    f_dx = similar(f_x)

    g_y = sol_t[
        algebraic_var_Idx_in_state]

    g_dy = similar(g_y)

    #----------------------------------------
    
    ∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> ode_gens_plants_generic_model_func!(
            dx, x,
            ωref_vref_porder0_id_iq_vh,
            tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, f_x )

    
    ∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> ode_gens_plants_generic_model_func!(
            dx, f_x, p, tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, ωref_vref_porder0_id_iq_vh )
    
    #----------------------------------------
    
    ∂g∂y = ForwardDiff.jacobian(
        (dy, y ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, y,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                    algebraic_generic_model_kwd_para ),
        g_dy, g_y )

    
    ∂g∂p = ForwardDiff.jacobian(
        (dy, p ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, g_y, p;
            kwd_para = algebraic_generic_model_kwd_para
                     ),
        g_dy, δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para )

    #----------------------------------------

    a =  ∂g∂p[g_rows_P_ng_idx, g_PQ_ng_idx]\
            ∂g∂y[g_rows_P_ng_idx, vθ_g_idx ] 
    
    b = ∂g∂p[g_rows_P_ng_idx, g_PQ_ng_idx]\
            ∂g∂y[g_rows_P_ng_idx, vθ_ng_idx ] 
    
    c = ∂g∂y[g_rows_Q_ng_idx, vθ_g_idx] -
        ∂g∂p[g_rows_Q_ng_idx, g_PQ_ng_idx ] * a
    
    d = ∂g∂y[g_rows_Q_ng_idx, vθ_ng_idx] -
        ∂g∂p[g_rows_Q_ng_idx, g_PQ_ng_idx ] * b


    e =  ∂g∂y[g_rows_idq_idx, g_i_dq_idx]\
            ∂g∂y[g_rows_idq_idx, vθ_g_idx ] 


    f =  ∂g∂y[g_rows_idq_idx, g_i_dq_idx]\
            ∂g∂p[g_rows_idq_idx, g_x_δ_edash_idx ] 

    h = ∂g∂y[g_rows_P_g_idx, vθ_g_idx ] -
        ∂g∂y[g_rows_P_g_idx, vθ_ng_idx ] * (d\c) -
        ∂g∂y[g_rows_P_g_idx, g_i_dq_idx ] * e

    i = ∂g∂p[g_rows_P_g_idx, g_x_δ_edash_idx ] -
        ∂g∂y[g_rows_P_g_idx, g_i_dq_idx ] * f

    j = ∂g∂y[g_rows_Q_g_idx, vθ_g_idx] -
        ∂g∂y[g_rows_Q_g_idx, vθ_ng_idx] * (d\c) -
        ∂g∂y[g_rows_Q_g_idx, g_i_dq_idx] * e

    k = ∂g∂p[g_rows_Q_g_idx, g_x_δ_edash_idx ] -
        ∂g∂y[g_rows_Q_g_idx, g_i_dq_idx] * f

    l = loc_load_exist == true ? (
        h - ∂g∂p[g_rows_P_g_idx, g_PQll_idx ] * (
            ∂g∂p[g_rows_Q_g_idx, g_PQll_idx ] \ j) ) :
                j # j + h
    
    m = loc_load_exist == true ? (
        i - ∂g∂p[g_rows_P_g_idx, g_PQll_idx ] * (
            ∂g∂p[g_rows_Q_g_idx, g_PQll_idx ] \ k) ) :
                k # k + i

    N = (l\m)

    Γ = ∂f∂p[:, f_i_dq_idx] * (e * N - f) -
        ∂f∂p[:, f_vg_Idx] * N[ 1:N_g, :]

    # Asys = ∂f∂x[state_var_Idx_in_state,
    #             state_var_Idx_in_state ]

    Asys = copy(∂f∂x)
    
    #
    Asys[:,  δ_ed_dash_eq_dash_idx_in_state] .=
        Asys[:, δ_ed_dash_eq_dash_idx_in_state] + Γ

    # Asys[:, g_x_δ_edash_idx] =
    #     Asys[:, g_x_δ_edash_idx] + Γ
    
    Λ = ∂f∂u =  ∂f∂p[:, f_u_idx]
    
    #----------------------------------------

    # Δ⨰ = Asys Δx + Λ Δu
    
    # where Δu = Δωref0_vref0_porder0
    
    #----------------------------------------
    #----------------------------------------

    sys_eigvalues = eigvals(Asys)

    printed_sys_eigvalues =
        round.( sys_eigvalues; digits = 4 )
    
    return (; Asys,  Λ,
            sys_eigvalues,
            printed_sys_eigvalues )
    
end


function get_generic_linearised_dynamic_model(
    sol,
    tt,
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
    nothing;
    use_init_condition =
        false,
    kwd_para =
        linearised_dynamic_kwd_para )

    if use_init_condition == true

        sol_t = sol
    else
        sol_t = sol(tt)
    end
    

    (;loc_load_exist,

     pf_generic_gens_para,         
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,
     gens_state_vars_idx_in_state,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
     
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_n2s_idxs,

     ode_plants_kwd_para) =
         kwd_para

    #----------------------------------------
    
    algebraic_generic_model_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes )    

    #----------------------------------------
    
    # (;loc_load_exist,
     
    #  state_vars_and_i_dq_Idx_in_state,
    #  state_algebraic_vars_Idx_in_state,
    #  gens_state_vars_idx_in_state,
     
    #  dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
    #  dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
    #  dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
    #  dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
    #  dyn_pf_fun_kwd_net_idxs,
    #  dyn_pf_fun_kwd_n2s_idxs,

    #  ode_plants_kwd_para) =
    #      kwd_para

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------    
    
    (f_ωref0_Idx,
     f_vref0_Idx,
     f_porder0_Idx,
     f_id_Idx,
     f_iq_Idx,
     f_vg_Idx ) =
         NamedTupleTools.select(
             dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
             (:dyn_ωref0_Idx,
              :dyn_vref0_Idx,
              :dyn_porder0_Idx,
              :dyn_id_Idx,
              :dyn_iq_Idx,
              :dyn_vh_Idx))
    
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

    (slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))

    (;n2s_slack_gens_idx,
     n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------        
    #----------------------------------------    

    f_u_idx = first(f_ωref0_Idx):last(f_porder0_Idx)
    
    f_i_dq_idx = first(f_id_Idx):last(f_iq_Idx)

    
    #----------------------------------------    

    g_vθ_idx =
        first(dyn_pf_vh_Idxs):last(dyn_pf_θh_Idxs)
    
    vθ_g_idx  =
        [dyn_pf_vh_Idxs[gens_nodes_idx];
         dyn_pf_θh_Idxs[gens_nodes_idx]]
    
    vθ_ng_idx =
        [dyn_pf_vh_Idxs[non_gens_nodes_idx];
         dyn_pf_θh_Idxs[non_gens_nodes_idx]]

    g_i_dq_idx =
        first(dyn_pf_id_Idxs):last(dyn_pf_iq_Idxs)

    ###
    
    g_rows_P_idx = dyn_pf_vh_Idxs
    g_rows_Q_idx = dyn_pf_θh_Idxs
    
    g_rows_P_g_idx = g_rows_P_idx[gens_nodes_idx]
    g_rows_Q_g_idx = g_rows_Q_idx[gens_nodes_idx]
    
    g_rows_P_ng_idx = g_rows_P_idx[non_gens_nodes_idx]
    g_rows_Q_ng_idx = g_rows_Q_idx[non_gens_nodes_idx]

    g_rows_idq_idx =
        first(dyn_pf_id_Idxs):last(dyn_pf_iq_Idxs)
    
    #----------------------------------------    

    g_x_δ_edash_idx =
        first(dyn_δ_Idxs):last(dyn_eq_dash_Idxs)
    
    g_PQ_ng_idx = 
        first(dyn_P_non_gens_Idxs):last(
            dyn_Q_gens_loc_load_Idxs)

    g_PQll_idx = loc_load_exist == true ? 
        (first(dyn_P_gens_loc_load_Idxs):last(
            dyn_Q_gens_loc_load_Idxs)) : []

    #----------------------------------------    

    δ_ed_dash_eq_dash_idx_in_state =
        [δ_idx_in_state;
         ed_dash_idx_in_state;
         eq_dash_idx_in_state ]

    #----------------------------------------    

    # dyn_pf_vh_Idxs
    # dyn_pf_θh_Idxs
    # dyn_pf_id_Idxs
    # dyn_pf_iq_Idxs

    # slack_gens_nodes_idx
    
    # gens_nodes_idx
    
    # non_gens_nodes_idx

    dyn_pf_gens_vh_Idxs =
        dyn_pf_vh_Idxs[gens_nodes_idx]

    dyn_pf_non_gens_vh_Idxs =
        dyn_pf_vh_Idxs[non_gens_nodes_idx]

    dyn_pf_slack_θh_Idxs =
        dyn_pf_θh_Idxs[slack_gens_nodes_idx]

    dyn_pf_non_slack_θh_Idxs =
        setdiff(collect(dyn_pf_θh_Idxs),
                dyn_pf_θh_Idxs[slack_gens_nodes_idx] )
    
    JLF_idx = [dyn_pf_non_gens_vh_Idxs ;
               dyn_pf_non_slack_θh_Idxs]
    
    algeb_by_state_link_cols_idx =
        [dyn_pf_gens_vh_Idxs;
         # dyn_pf_slack_θh_Idxs;
         dyn_pf_id_Idxs;
         dyn_pf_iq_Idxs]

    # δ_ed_dash_eq_dash_idx_in_state
    state_by_algeb_link_cols_idx =
        [δ_idx_in_state;
         ed_dash_idx_in_state;
         eq_dash_idx_in_state]
    
    #----------------------------------------    
    #----------------------------------------

    vh  =  sol_t[ vh_Idx_in_state]
    
    θh  =  sol_t[ θh_Idx_in_state]
    
    gens_i_d =  sol_t[id_Idx_in_state]
    
    gens_i_q =  sol_t[iq_Idx_in_state]


    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]
    
    #----------------------------------------
        
    gens_δ = sol_t[
        δ_idx_in_state]
        
    gens_ed_dash = sol_t[
        ed_dash_idx_in_state]
    
    gens_eq_dash = sol_t[
        eq_dash_idx_in_state]    

    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_v_ref_Idx]
    
    p_order =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_p_order_Idx]
    
    P_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Png_Idx]
    
    Q_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    #----------------------------------------
    
    ωref_vref_porder0_id_iq_vh =
        Float64[ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]


    #----------------------------------------    

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        Float64[gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------    

    f_x  = sol_t[
        state_var_Idx_in_state]
    
    f_dx = similar(f_x)

    g_y = sol_t[
        algebraic_var_Idx_in_state]

    g_dy = similar(g_y)

    #----------------------------------------
    
    ∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> ode_gens_plants_generic_model_func!(
            dx, x,
            ωref_vref_porder0_id_iq_vh,
            tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, f_x )

    
    ∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> ode_gens_plants_generic_model_func!(
            dx, f_x, p, tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, ωref_vref_porder0_id_iq_vh )
    
    #----------------------------------------
    
    ∂g∂y = ForwardDiff.jacobian(
        (dy, y ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, y,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                    algebraic_generic_model_kwd_para ),
        g_dy, g_y )

    
    ∂g∂p = ForwardDiff.jacobian(
        (dy, p ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, g_y, p;
            kwd_para = algebraic_generic_model_kwd_para
                     ),
        g_dy, δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para )

    #----------------------------------------

    λ = loc_load_exist == true ? (
        ∂g∂p[g_rows_Q_g_idx, g_PQll_idx]\
            ∂g∂y[g_rows_Q_g_idx, vθ_g_idx ]) : []
    
    ζ = loc_load_exist == true ? (
        ∂g∂p[g_rows_Q_g_idx, g_PQll_idx]\
            ∂g∂y[g_rows_Q_g_idx, vθ_ng_idx ]) : []
    
    ν = loc_load_exist == true ? (
        ∂g∂p[g_rows_Q_g_idx, g_PQll_idx]\
            ∂g∂y[g_rows_Q_g_idx, g_i_dq_idx ]) : []
    
    η = loc_load_exist == true ? (
        ∂g∂p[g_rows_Q_g_idx, g_PQll_idx]\
            ∂g∂p[g_rows_Q_g_idx, g_x_δ_edash_idx ]) : []

    χ = loc_load_exist == true ? (
        ∂g∂y[g_rows_P_g_idx, vθ_g_idx] -
            ∂g∂p[g_rows_P_g_idx, g_PQll_idx] * (
                ∂g∂p[g_rows_Q_g_idx, g_PQll_idx]\
                    ∂g∂y[g_rows_Q_g_idx, vθ_g_idx ]) ) :
                    ∂g∂y[g_rows_P_g_idx, vθ_g_idx]
    
    κ = loc_load_exist == true ? (
        ∂g∂y[g_rows_P_g_idx, vθ_ng_idx] -
            ∂g∂p[g_rows_P_g_idx, g_PQll_idx] * (
                ∂g∂p[g_rows_Q_g_idx, g_PQll_idx]\
                    ∂g∂y[g_rows_Q_g_idx, vθ_ng_idx ]) ) :
                    ∂g∂y[g_rows_P_g_idx, vθ_ng_idx]
    
    ψ = loc_load_exist == true ? (
        ∂g∂y[g_rows_P_g_idx, g_i_dq_idx] -
            ∂g∂p[g_rows_P_g_idx, g_PQll_idx] * (
                ∂g∂p[g_rows_Q_g_idx, g_PQll_idx] \
                    ∂g∂p[g_rows_Q_g_idx, g_i_dq_idx])) :
                    ∂g∂y[g_rows_P_g_idx, g_i_dq_idx]
    
    ξ = loc_load_exist == true ? (
        ∂g∂p[g_rows_P_g_idx, g_x_δ_edash_idx] -
            ∂g∂p[g_rows_P_g_idx, g_PQll_idx] * (
                ∂g∂p[g_rows_Q_g_idx, g_PQll_idx] \
                    ∂g∂p[g_rows_Q_g_idx,
                         g_x_δ_edash_idx])) :
                    ∂g∂p[g_rows_P_g_idx, g_x_δ_edash_idx]

    α = ∂g∂p[g_rows_Q_ng_idx, g_PQ_ng_idx] \
        ∂g∂y[g_rows_Q_ng_idx, vθ_g_idx ]
    
    β = ∂g∂p[g_rows_Q_ng_idx, g_PQ_ng_idx] \
        ∂g∂y[g_rows_Q_ng_idx, vθ_ng_idx ]

    σ = ∂g∂y[g_rows_P_ng_idx, vθ_g_idx ] -
        ∂g∂p[g_rows_P_ng_idx, g_PQ_ng_idx] * α
    
    ρ = ∂g∂y[g_rows_P_ng_idx, vθ_ng_idx ] -
        ∂g∂p[g_rows_P_ng_idx, g_PQ_ng_idx] * β

    A = loc_load_exist == true ? (
        ∂g∂y[g_rows_P_g_idx, vθ_g_idx] -
            ∂g∂y[g_rows_P_g_idx, vθ_ng_idx] * (ρ\σ) -
            ∂g∂p[g_rows_P_g_idx, g_PQll_idx] * (λ - ζ *
            (ρ\σ))) :
            (∂g∂y[g_rows_P_g_idx, vθ_g_idx] -
            ∂g∂y[g_rows_P_g_idx, vθ_ng_idx] * (ρ\σ) )
    
    B = loc_load_exist == true ? (
        ∂g∂y[g_rows_P_g_idx, g_i_dq_idx] -
            ∂g∂p[g_rows_P_g_idx, g_PQll_idx] * ν) :
            ∂g∂y[g_rows_P_g_idx, g_i_dq_idx]
    
    C = loc_load_exist == true ? (
        ∂g∂p[g_rows_P_g_idx, g_x_δ_edash_idx] -
            ∂g∂p[g_rows_P_g_idx, g_PQll_idx] * η) :
            ∂g∂p[g_rows_P_g_idx, g_x_δ_edash_idx]

    υ = ∂g∂y[g_rows_idq_idx, g_i_dq_idx] \
        ∂g∂y[g_rows_idq_idx, vθ_g_idx]
    
    ϕ = ∂g∂y[g_rows_idq_idx, g_i_dq_idx] \
        ∂g∂p[g_rows_idq_idx, g_x_δ_edash_idx]

    D = A - B * υ
    E = C - B * ϕ

    F = -(D\E)
    
    Γ = ∂f∂x[:, f_i_dq_idx] * (υ * F - ϕ) -
        ∂f∂p[:, f_vg_Idx] * F[gens_nodes_idx, :]

    Asys = ∂f∂x[state_var_Idx_in_state,
                state_var_Idx_in_state ]

    # Asys[:, δ_ed_dash_eq_dash_idx_in_state ] =
    #     Asys[:, δ_ed_dash_eq_dash_idx_in_state ] +  Γ

    Asys[:, δ_ed_dash_eq_dash_idx_in_state ] +=  Γ

    Λ = ∂f∂u =  ∂f∂p[:, f_u_idx]
    
    #----------------------------------------

    # Δ⨰ = Asys Δx + Λ Δu
    
    # where Δu = Δωref0_vref0_porder0
    
    #----------------------------------------
    #----------------------------------------

    
    sys_eigvalues = eigvals(Asys)

    printed_sys_eigvalues =
        round.( sys_eigvalues; digits = 4 )

    return (; Asys,  Λ, sys_eigvalues,
            printed_sys_eigvalues )
    
end


function get_generic_linearised_dynamic_model(
    sol,
    tt,
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll;
    use_init_condition =
        false,
    kwd_para =
        linearised_dynamic_kwd_para )

    if use_init_condition == true

        sol_t = sol
    else
        sol_t = sol(tt)
    end
    

    (;loc_load_exist,

     pf_generic_gens_para,         
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,
     gens_state_vars_idx_in_state,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
     
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_n2s_idxs,

     ode_plants_kwd_para) =
         kwd_para

    #----------------------------------------
    
    algebraic_generic_model_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes )    

    #----------------------------------------
    
    # (;loc_load_exist,
     
    #  state_vars_and_i_dq_Idx_in_state,
    #  state_algebraic_vars_Idx_in_state,
    #  gens_state_vars_idx_in_state,
     
    #  dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
    #  dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
    #  dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
    #  dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
    #  dyn_pf_fun_kwd_net_idxs,
    #  dyn_pf_fun_kwd_n2s_idxs,

    #  ode_plants_kwd_para) =
    #      kwd_para

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------    
    
    (f_ωref0_Idx,
     f_vref0_Idx,
     f_porder0_Idx,
     f_id_Idx,
     f_iq_Idx,
     f_vg_Idx ) =
         NamedTupleTools.select(
             dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
             (:dyn_ωref0_Idx,
              :dyn_vref0_Idx,
              :dyn_porder0_Idx,
              :dyn_id_Idx,
              :dyn_iq_Idx,
              :dyn_vh_Idx))
    
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

    (slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))

    (;n2s_slack_gens_idx,
     n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    

    #----------------------------------------        
    #----------------------------------------    

    f_u_idx = first(f_ωref0_Idx):last(f_porder0_Idx)
    
    f_i_dq_idx = first(f_id_Idx):last(f_iq_Idx)

    # f_vg_Idx
    
    #----------------------------------------    

    g_vθ_idx =
        first(dyn_pf_vh_Idxs):last(dyn_pf_θh_Idxs)

    g_i_dq_idx =
        first(dyn_pf_id_Idxs):last(dyn_pf_iq_Idxs)

    #----------------------------------------    

    g_x_δ_edash_idx =
        first(dyn_δ_Idxs):last(dyn_eq_dash_Idxs)

    g_PQl_idx = loc_load_exist == true ? 
        (first(dyn_P_non_gens_Idxs):last(
            dyn_Q_gens_loc_load_Idxs)) :
        (first(dyn_P_non_gens_Idxs):last(
            dyn_Q_non_gens_Idxs))
    
    #----------------------------------------    

    δ_ed_dash_eq_dash_idx_in_state =
        [δ_idx_in_state;
         ed_dash_idx_in_state;
         eq_dash_idx_in_state ]

    #----------------------------------------    

    # dyn_pf_vh_Idxs
    # dyn_pf_θh_Idxs
    # dyn_pf_id_Idxs
    # dyn_pf_iq_Idxs

    # slack_gens_nodes_idx
    
    # gens_nodes_idx
    
    # non_gens_nodes_idx

    dyn_pf_gens_vh_Idxs =
        dyn_pf_vh_Idxs[gens_nodes_idx]

    dyn_pf_non_gens_vh_Idxs =
        dyn_pf_vh_Idxs[non_gens_nodes_idx]

    dyn_pf_slack_θh_Idxs =
        dyn_pf_θh_Idxs[slack_gens_nodes_idx]

    dyn_pf_non_slack_θh_Idxs =
        setdiff(collect(dyn_pf_θh_Idxs),
                dyn_pf_θh_Idxs[slack_gens_nodes_idx] )
    
    JLF_idx = [dyn_pf_non_gens_vh_Idxs ;
               dyn_pf_non_slack_θh_Idxs]
    
    algeb_by_state_link_cols_idx =
        [dyn_pf_gens_vh_Idxs;
         # dyn_pf_slack_θh_Idxs;
         dyn_pf_id_Idxs;
         dyn_pf_iq_Idxs]

    # δ_ed_dash_eq_dash_idx_in_state
    state_by_algeb_link_cols_idx =
        [δ_idx_in_state;
         ed_dash_idx_in_state;
         eq_dash_idx_in_state]
    
    #----------------------------------------    
    #----------------------------------------

    vh  =  sol_t[ vh_Idx_in_state]
    
    θh  =  sol_t[ θh_Idx_in_state]
    
    gens_i_d =  sol_t[id_Idx_in_state]
    
    gens_i_q =  sol_t[iq_Idx_in_state]


    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]
    
    #----------------------------------------
        
    gens_δ = sol_t[
        δ_idx_in_state]
        
    gens_ed_dash = sol_t[
        ed_dash_idx_in_state]
    
    gens_eq_dash = sol_t[
        eq_dash_idx_in_state]    

    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_v_ref_Idx]
    
    p_order =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_p_order_Idx]
    
    P_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Png_Idx]
    
    Q_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    #----------------------------------------
    
    ωref_vref_porder0_id_iq_vh =
        Float64[ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]


    #----------------------------------------    

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        Float64[gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------    

    f_x  = sol_t[
        state_var_Idx_in_state]
    
    f_dx = similar(f_x)

    g_y = sol_t[
        algebraic_var_Idx_in_state]

    g_dy = similar(g_y)

    #----------------------------------------
    
    ∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> ode_gens_plants_generic_model_func!(
            dx, x,
            ωref_vref_porder0_id_iq_vh,
            tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, f_x )

    
    ∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> ode_gens_plants_generic_model_func!(
            dx, f_x, p, tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, ωref_vref_porder0_id_iq_vh )
    
    #----------------------------------------
    
    ∂g∂y = ForwardDiff.jacobian(
        (dy, y ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, y,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                    algebraic_generic_model_kwd_para ),
        g_dy, g_y )

    
    ∂g∂p = ForwardDiff.jacobian(
        (dy, p ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, g_y, p;
            kwd_para = algebraic_generic_model_kwd_para
                     ),
        g_dy, δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para )

    #----------------------------------------

    # L = ∂f∂x[state_var_Idx_in_state,
    #          state_var_Idx_in_state ]
    
    # Asys = ∂f∂x[state_var_Idx_in_state,
    #             state_var_Idx_in_state ]

    
    Asys = copy(∂f∂x)
    
    Q = ∂f∂p[state_var_Idx_in_state, f_u_idx ]
    M = ∂f∂p[state_var_Idx_in_state, f_i_dq_idx ]
    N = ∂f∂p[state_var_Idx_in_state, f_vg_Idx ]
    
    #----------------------------------------

    A = ∂g∂y[g_vθ_idx, g_vθ_idx ]
    B = ∂g∂y[g_vθ_idx, g_i_dq_idx ]
    C = ∂g∂y[g_i_dq_idx, g_vθ_idx ]
    D = ∂g∂y[g_i_dq_idx, g_i_dq_idx ]

    E = ∂g∂p[g_vθ_idx, g_x_δ_edash_idx]
    F = ∂g∂p[g_vθ_idx, g_PQl_idx]
    G = ∂g∂p[g_i_dq_idx, g_x_δ_edash_idx]
    H = ∂g∂p[g_i_dq_idx, g_PQl_idx]

    
    #----------------------------------------

    I = A - B * (D\C)
    J = E - B * (D\G)
    K = F - B * (D\H)

    α = I\J
    β = I\K

    #----------------------------------------

    R = M * (D \ C) * α - M * (D \ G)
    S = M * (D \ C) * β - M * (D \ H)

    γ = α[gens_nodes_idx, :]
    ρ = β[gens_nodes_idx, :]


    ψ = R - N * γ
    ϕ = S - N * ρ


    Asys[:, δ_ed_dash_eq_dash_idx_in_state ] .=
        Asys[:, δ_ed_dash_eq_dash_idx_in_state ] + ψ

    # Δ⨰ = Asys Δx +  ϕ ΔPQng_PQll + Q Δu
    
    # where Δu = Δωref0_vref0_porder0
    
    #----------------------------------------
    #----------------------------------------

    size_states = length(state_var_Idx_in_state)
    size_algebr = length(algebraic_var_Idx_in_state)

    JAE = ∂g∂y

    JLF = JAE[JLF_idx, JLF_idx ]
    
    
    Γ = zeros(size_states, size_algebr )
    Λ = zeros(size_algebr, size_states )

    Γ[state_var_Idx_in_state,
           algeb_by_state_link_cols_idx] .=
               hcat(∂f∂p[state_var_Idx_in_state, f_vg_Idx ],
                    ∂f∂p[state_var_Idx_in_state,
                         f_i_dq_idx ] )
    
    # Λ[algebraic_var_Idx_in_state,
    #        state_by_algeb_link_cols_idx] .=
    #            ∂g∂p[:, g_x_δ_edash_idx]

    
    Λ[:, state_by_algeb_link_cols_idx] .=
               ∂g∂p[:, g_x_δ_edash_idx]
    
    # Asys = A - B * ( JAE \ C)
    
    Asys_dash = ∂f∂x - Γ * ( ∂g∂y \ Λ )
    
    #----------------------------------------
    #----------------------------------------
    
    sys_eigvalues = eigvals(Asys)

    printed_sys_eigvalues =
        round.( sys_eigvalues; digits = 4 )

    
    sys_eigvalues_dash  = eigvals(Asys_dash )

    printed_sys_eigvalues_dash  =
        round.( sys_eigvalues_dash ; digits = 4 )
    
    return (; Asys, ϕ, Q, JAE, JLF, Asys_dash,
            sys_eigvalues,
            printed_sys_eigvalues,
            sys_eigvalues_dash,
            printed_sys_eigvalues_dash)
    
end


function get_generic_linearised_dynamic_model(
    sol_tt_or_system_init,
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll;
    kwd_para =
        linearised_dynamic_kwd_para )

    sol_t = sol_tt_or_system_init

    (;loc_load_exist,

     pf_generic_gens_para,         
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,
     gens_state_vars_idx_in_state,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
     
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_n2s_idxs,

     ode_plants_kwd_para) =
         kwd_para

    #----------------------------------------
    
    algebraic_generic_model_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes )    
    
    #----------------------------------------

        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
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

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
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
    
    
    #----------------------------------------    

    g_vθ_idx =
        first(dyn_pf_vh_Idxs):last(dyn_pf_θh_Idxs)

    g_i_dq_idx =
        first(dyn_pf_id_Idxs):last(dyn_pf_iq_Idxs)

    #----------------------------------------    

    g_x_δ_edash_idx =
        first(dyn_δ_Idxs):last(dyn_eq_dash_Idxs)

    g_PQl_idx = loc_load_exist == true ? 
        (first(dyn_P_non_gens_Idxs):last(
            dyn_Q_gens_loc_load_Idxs)) :
                (first(dyn_P_non_gens_Idxs):last(
                    dyn_Q_non_gens_Idxs))

    #----------------------------------------    
    
    (f_ωref0_Idx,
     f_vref0_Idx,
     f_porder0_Idx,
     f_id_Idx,
     f_iq_Idx,
     f_vg_Idx ) =
         NamedTupleTools.select(
             dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
             (:dyn_ωref0_Idx,
              :dyn_vref0_Idx,
              :dyn_porder0_Idx,
              :dyn_id_Idx,
              :dyn_iq_Idx,
              :dyn_vh_Idx))

    # f_ωref0_Idx =
    #     ωref0_vref0_porder0_id_iq_vh_Idx.ωref0_Idx
    
    # f_vref0_Idx =
    #     ωref0_vref0_porder0_id_iq_vh_Idx.vref0_Idx
    
    # f_porder0_Idx =
    #     ωref0_vref0_porder0_id_iq_vh_Idx.porder0_Idx
    
    # f_id_Idx =
    #     ωref0_vref0_porder0_id_iq_vh_Idx.id_Idx
    
    # f_iq_Idx =
    #     ωref0_vref0_porder0_id_iq_vh_Idx.iq_Idx
    
    # f_vg_Idx =
    #     ωref0_vref0_porder0_id_iq_vh_Idx.vh_Idx

    f_u_idx = first(f_ωref0_Idx):last(f_porder0_Idx)
    f_i_dq_idx = first(f_id_Idx):last(f_iq_Idx)


    #----------------------------------------

    δ_ed_dash_eq_dash_idx_in_state =
        [δ_idx_in_state;
         ed_dash_idx_in_state;
         eq_dash_idx_in_state ]

    
    #----------------------------------------    
    #----------------------------------------

    vh  =  sol_t[ vh_Idx_in_state]
    
    θh  =  sol_t[ θh_Idx_in_state]
    
    gens_i_d =  sol_t[id_Idx_in_state]
    
    gens_i_q =  sol_t[iq_Idx_in_state]


    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]
    
    #----------------------------------------
        
    gens_δ = sol_t[
        δ_idx_in_state]
        
    gens_ed_dash = sol_t[
        ed_dash_idx_in_state]
    
    gens_eq_dash = sol_t[
        eq_dash_idx_in_state]    

    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_v_ref_Idx]
    
    p_order =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_p_order_Idx]
    
    P_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Png_Idx]
    
    Q_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    #----------------------------------------
    
    ωref_vref_porder_id_iq_vh =
        Float64[ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]


    #----------------------------------------    

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        Float64[gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------    

    f_x  = sol_t[
        state_var_Idx_in_state]
    
    f_dx = similar(f_x)

    g_y = sol_t[
        algebraic_var_Idx_in_state]

    g_dy = similar(g_y)

    #----------------------------------------
    
    ∂f∂x = ForwardDiff.jacobian(
        (dx, x ) -> ode_gens_plants_generic_model_func!(
            dx, x,
            ωref_vref_porder_id_iq_vh,
            tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, f_x )

    
    ∂f∂p = ForwardDiff.jacobian(
        (dx, p ) -> ode_gens_plants_generic_model_func!(
            dx, f_x, p, tt;
            kwd_para =
                    ode_plants_kwd_para ),
        f_dx, ωref_vref_porder_id_iq_vh )
    
    #----------------------------------------
    
    ∂g∂y = ForwardDiff.jacobian(
        (dy, y ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, y,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                    algebraic_generic_model_kwd_para ),
        g_dy, g_y )

    
    ∂g∂p = ForwardDiff.jacobian(
        (dy, p ) -> algebraic_generic_pf_ΔPQ_mismatch!(
            dy, g_y, p;
            kwd_para = algebraic_generic_model_kwd_para
                     ),
        g_dy, δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para )

    #----------------------------------------

    A = ∂g∂y[g_vθ_idx, g_vθ_idx ]
    B = ∂g∂y[g_vθ_idx, g_i_dq_idx ]
    C = ∂g∂y[g_i_dq_idx, g_vθ_idx ]
    D = ∂g∂y[g_i_dq_idx, g_i_dq_idx ]

    E = ∂g∂p[g_vθ_idx, g_x_δ_edash_idx]
    F = ∂g∂p[g_vθ_idx, g_PQl_idx]
    G = ∂g∂p[g_i_dq_idx, g_x_δ_edash_idx]
    H = ∂g∂p[g_i_dq_idx, g_PQl_idx]

    I = A - B * (D\C)
    J = E - B * (D\G)
    K = F - B * (D\H)

    α = I\J
    β = I\K

    Asys = ∂f∂x[state_var_Idx_in_state,
                state_var_Idx_in_state ]
    Q = ∂f∂p[state_var_Idx_in_state, f_u_idx ]
    M = ∂f∂p[state_var_Idx_in_state, f_i_dq_idx ]
    N = ∂f∂p[state_var_Idx_in_state, f_vg_Idx ]

    R = M * (D \ C) * α - M * (D \ G)
    S = M * (D \ C) * β - M * (D \ H)

    γ = α[gens_nodes_idx, :]
    ρ = β[gens_nodes_idx, :]


    ψ = R - N * γ
    ϕ = S - N * ρ


    Asys[:, δ_ed_dash_eq_dash_idx_in_state ] .=
        Asys[:, δ_ed_dash_eq_dash_idx_in_state ] + ψ

    # Δ⨰ = Asys Δx +  ϕ ΔPQng_PQll + Q Δu
    
    # where Δu = Δωref0_vref0_porder0

    #----------------------------------------
    
    sys_eigvalues = eigvals(Asys)

    printed_sys_eigvalues =
        round.( sys_eigvalues; digits = 4 )

    #----------------------------------------
    
    return (; Asys, ϕ, Q,
            sys_eigvalues,
            printed_sys_eigvalues )
    
end



function get_generic_small_signal_stability_indices(
    Asys)

    eigen_object = eigen( Asys )

    eig_values = eigen_object.values

    printed_eig_values =
        round.(eig_values; digits=4)

    eigvecs_right = eigen_object.vectors

    # inv_eigvecs_right = inv(eigvecs_right)
    
    inv_eigvecs_right =
        eigvecs_right \ LinearAlgebra.I(
            size(eigvecs_right)[1] )


    #--------------------------------------------------
    # states associated with eig_values
    #--------------------------------------------------

    M_diag = inv_eigvecs_right * Asys * eigvecs_right

    printed_M_diag =
        round.( M_diag; digits = 4 )

    #--------------------------------------------------
    # participation factor
    #--------------------------------------------------

    PF_Asys  = get_participation_factors( Asys )

    printed_PF_Asys =
        round.( PF_Asys; digits = 4 )

    #--------------------------------------------------

    eig_values, eigvecs_left, eigvecs_right =
        get_eigens(Asys)


    #--------------------------------------------------

    return (; PF_Asys, printed_PF_Asys,
            M_diag, printed_M_diag,
            eig_values, printed_eig_values,
            eigvecs_left, eigvecs_right,
            inv_eigvecs_right )
    

end
     


#-------------------------------------------------------
#-------------------------------------------------------
