# (c) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.4:  pg 135 - 136, 6.121 - 6.123



########################################################
# ------------------------------------------------------
#  ODE models
# ------------------------------------------------------
########################################################

#-------------------------------------------------------
# one model im
#-------------------------------------------------------

function ode_one_im_model_with_only_δ_func!(
    dx, x,
    (stateDiffCache, ode_para), t)

    (; Ax_Bx_Cx_views,
     gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view,
     node_i) = ode_para

    (; vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views)  = Ax_Bx_Cx_views

    Ax_row_size, Ax_col_size =
        size(vec_Ax_views[ node_i ])

    Bx_row_size, Bx_col_size =
        size(vec_Bx_views[ node_i ])

    Cx_row_size, Cx_col_size =
        size(vec_Cx_views[ node_i ])
    
    
    id_iq_pg =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view[
            node_i ]

    ωs_ωref0_vref0_porder0 =
        gens_nodes_ωs_ωref0_vref0_porder0_view[
            node_i ]

    Ax_row_size, Ax_col_size
    
    dx_gen_δ = dx[ 1 ]
    
    x_gen_δ  = x[ 1 ]
        
    Ax_m = reshape(
        Matrix( vec_Ax_views[ node_i ])[
            1 , 1:Ax_col_size ],
        1,  Ax_col_size)
    
    Bx_m = reshape(
        Matrix( vec_Bx_views[ node_i ])[
            1, 1:Bx_col_size ],
        1, Bx_col_size) 
    
    Cx_m = reshape(
        Matrix( vec_Cx_views[ node_i ])[
            1, 1:Cx_col_size ],
        1, Cx_col_size) 

    dx[1:1] .= Ax_m * stateDiffCache +
        Bx_m * id_iq_pg  +
        Cx_m * ωs_ωref0_vref0_porder0

    return nothing
    
end


function ode_one_im_model_without_δ_func!(
    dx,
    x,
    ode_para,    
    t)

    (; Ax_Bx_Cx_views,
     gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view,
     node_i ) = ode_para

    (; vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views ) = Ax_Bx_Cx_views

    Ax_row_size, Ax_col_size =
        size(vec_Ax_views[ node_i ])

    Bx_row_size, Bx_col_size =
        size(vec_Bx_views[ node_i ])

    Cx_row_size, Cx_col_size =
        size(vec_Cx_views[ node_i ])

        
    Ax_m = Matrix(
        vec_Ax_views[ node_i ] )[
            2:Ax_row_size, 2:Ax_col_size ]
    
    Bx_m = Matrix(
        vec_Bx_views[ node_i ] )[
            2:Bx_row_size, 1:Bx_col_size ]
    
    Cx_m = Matrix(
        vec_Cx_views[ node_i ])[
            2:Cx_row_size, 1:Cx_col_size ]   
    
    id_iq_pg =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view[
            node_i ]

    ωs_ωref0_vref0_porder0 =
        gens_nodes_ωs_ωref0_vref0_porder0_view[
            node_i ]

    Ax_row_size, Ax_col_size
    
    dx_gen_no_δ = @view dx[ 2:Ax_row_size ]
    
    x_gen_no_δ  = @view x[ 2:Ax_row_size ]
    
    dx_gen_no_δ .= Ax_m * x_gen_no_δ +
        Bx_m * id_iq_pg  +
        Cx_m * ωs_ωref0_vref0_porder0
    
    return nothing
    
end


function ode_one_im_model_func!(
    dx,
    x,
    ode_para,
    t)


    (Ax_Bx_Cx_views,
     gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view,
     node_i) =
         ode_para

    (vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views)  =
         Ax_Bx_Cx_views
        
    Ax_m = vec_Ax_views[ node_i ]
    
    Bx_m = vec_Bx_views[ node_i ]
    
    Cx_m = vec_Cx_views[ node_i ]
    
    
    id_iq_pg =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view[ node_i ]

    ωs_ωref0_vref0_porder0 =
        gens_nodes_ωs_ωref0_vref0_porder0_view[ node_i ]
        
    dx .=
        Ax_m * x +
        Bx_m * id_iq_pg  +
        Cx_m * ωs_ωref0_vref0_porder0

    return nothing
    
end


function dynamics_one_im_model!(
    dx, x, para, t )

    (; netd,
     nodes_cb_sw,
     ode_para,
     plant_i,
     node_i,
     idxs,
     stateDiffCache) = para

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs) = idxs

    (; Ax_Bx_Cx_views,
     gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view,
     node_i ) =  ode_para 

    (; vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views)  = Ax_Bx_Cx_views

    
    vec_Ax_i =
        vec_Ax_views[ node_i ]
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------

    stateDiffCache =
        get_tmp(stateDiffCache, x)
    
    update_a_im_plant_system_matrices!(
        vec_Ax_i,
        stateDiffCache,
        plant_i )    
    
    #--------------------------------------------
    # ode gens pure state
    #--------------------------------------------

    ode_one_im_model_func!(
        dx,
        x,
        ode_para,
        t)

    # #--------------------------------------------
    # # Another method for all nodes votages
    # #--------------------------------------------

    # # im_all_nodes_voltage_terminal_func!(
    # #     dx,
    # #     x,
    # #     (nodes_pf_U_view, nodes_u_Idx_in_ranges),
    # #     t)
            
    # #--------------------------------------------
    # # non gens terminal voltages
    # #--------------------------------------------

    # non_gens_im_voltage_terminal_func!(
    #     dx, x, im_para_aux_inputs, t)
    
    # #--------------------------------------------
    # # gens terminal voltages
    # #--------------------------------------------

    # gens_im_voltage_terminal_func!(
    #     dx, x, im_para_aux_inputs, t)
    
    # #--------------------------------------------
    # #--------------------------------------------    
    
    return nothing
    
end



function t_ode_one_im_model_func!(
    dx, x, ode_para, t;
    sim_fun_kwd_para = sim_fun_kwd_para  )

    #--------------------------------------------
    
    (; Ax_m,
     Bx_m,
     Cx_m,
     stateDiffCache,
     sim_state_x0,
     plant_i,
     sim_fun_para_Idxs,
     ra_Xd_dash_Xq_dash_i ) = 
         sim_fun_kwd_para

    #--------------------------------------------

    vh_θh_i_Idx,
     ωs_ωref0_vref0_porder0_i_Idx  =
         sim_fun_para_Idxs

    #--------------------------------------------

    stateDiffCache = get_tmp(stateDiffCache, x)

    stateDiffCache .= sim_state_x0
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    update_a_im_plant_system_matrices!(
        Ax_m,
        stateDiffCache,
        plant_i )    

    #--------------------------------------------

    ωs_ωref0_vref0_porder0  =
        ode_para[ ωs_ωref0_vref0_porder0_i_Idx ] 

    vh_θh_i =
        ode_para[ vh_θh_i_Idx ]

    #--------------------------------------------

    vh, θh, = vh_θh_i

    ra, X_d_dash, X_q_dash = ra_Xd_dash_Xq_dash_i
    
    #--------------------------------------------
    
    id_iq =
        invZ_dq(ra, X_d_dash, X_q_dash) *
        [ x[3] - vh * sin( x[1] - θh),
          x[4] - vh * cos( x[1] - θh )]
    
    pg =  x[3] * id_iq[1] + x[4] * id_iq[2] +
        ( X_q_dash - X_d_dash ) *  *( id_iq... )

    id_iq_pg_vh = [ [id_iq, pg, vh ]...; ]
        
    dx .=
        Ax_m * x +
        Bx_m * id_iq_pg_vh  +
        Cx_m * ωs_ωref0_vref0_porder0

    return nothing
    
end


function s_ode_one_im_model_func!(
    dx, x , poi , t;
    poi_idx = poi_idx,
    sim_fun_kwd_para =
        s_sim_fun_kwd_para  )

    #--------------------------------------------
    
    (; Ax_m,
     Bx_m,
     Cx_m,
     stateDiffCache,
     poi_DiffCache,
     sim_state_x0,
     plant_i,
     ra_Xd_dash_Xq_dash_i,
     pois_dyn,
     poi_dyn_Idxs ) =
         sim_fun_kwd_para
    
    #--------------------------------------------

    stateDiffCache =
        get_tmp(stateDiffCache, x)

    stateDiffCache .=
        sim_state_x0

    # poi_DiffCache = get_tmp(poi_DiffCache, poi )

    # poi_DiffCache .= poi
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    update_a_im_plant_system_matrices!(
        Ax_m,
        stateDiffCache,
        plant_i )    

    #--------------------------------------------

    (; vh_θh_i_Idx,
     ωs_ωref0_vref0_porder0_i_Idx ) =
         poi_dyn_Idxs

    #--------------------------------------------

    (; gens_vh_θh_i,
     ωs_ωref0_vref0_porder0_i ) =
         pois_dyn

    #--------------------------------------------


    param_Idxs = collect(poi_dyn_Idxs)
        

    param_values = collect(pois_dyn)


    updated_para =
        [
            param_idx == poi_idx ?
                poi : param_value
                # poi_DiffCache : param_value
            for ( param_idx, param_value ) in
                zip(
                    param_Idxs,
                    param_values ) ]

    gens_vh_θh_i, ωs_ωref0_vref0_porder0_i = updated_para
        
    #--------------------------------------------
    #--------------------------------------------

    vh, θh, =  gens_vh_θh_i

    ra, X_d_dash, X_q_dash = ra_Xd_dash_Xq_dash_i
    
    #--------------------------------------------
    
    id_iq =
        invZ_dq(ra, X_d_dash, X_q_dash) *
        [ x[3] - vh * sin( x[1] - θh),
          x[4] - vh * cos( x[1] - θh )]
    
    pg =  x[3] * id_iq[1] + x[4] * id_iq[2] +
        ( X_q_dash - X_d_dash ) *  *( id_iq... )

    id_iq_pg_vh = [ [id_iq, pg, vh ]...; ]
        
    dx .=
        Ax_m * x +
        Bx_m * id_iq_pg_vh  +
        Cx_m * ωs_ωref0_vref0_porder0_i

    return nothing
    
end


function s_oop_ode_one_im_model_func!(
    x , poi , t;
    poi_idx = poi_idx,
    sim_fun_kwd_para = s_sim_fun_kwd_para  ) 

    #--------------------------------------------
    
    (; Ax_m,
     Bx_m,
     Cx_m,
     stateDiffCache,
     poi_DiffCache,
     sim_state_x0,
     plant_i,
     ra_Xd_dash_Xq_dash_i,
     pois_dyn,
     poi_dyn_Idxs ) =
         sim_fun_kwd_para
    
    #--------------------------------------------

    stateDiffCache =
        get_tmp(stateDiffCache, x)

    stateDiffCache .=
        sim_state_x0

    # poi_DiffCache = get_tmp(poi_DiffCache, poi )

    # poi_DiffCache .= poi
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    update_a_im_plant_system_matrices!(
        Ax_m,
        stateDiffCache,
        plant_i )    

    #--------------------------------------------

    (; vh_θh_i_Idx,
     ωs_ωref0_vref0_porder0_i_Idx ) =
         poi_dyn_Idxs

    #--------------------------------------------

    (; gens_vh_θh_i,
     ωs_ωref0_vref0_porder0_i ) =
         pois_dyn

    #--------------------------------------------


    param_Idxs = collect(poi_dyn_Idxs)
        

    param_values = collect(pois_dyn)


    updated_para =
        [
            param_idx == poi_idx ?
                poi : param_value
                # poi_DiffCache : param_value
            for ( param_idx, param_value ) in
                zip(
                    param_Idxs,
                    param_values ) ]

    gens_vh_θh_i, ωs_ωref0_vref0_porder0_i = updated_para
        
    #--------------------------------------------
    #--------------------------------------------

    vh, θh, =  gens_vh_θh_i

    ra, X_d_dash, X_q_dash = ra_Xd_dash_Xq_dash_i
    
    #--------------------------------------------
    
    id_iq =
        invZ_dq(ra, X_d_dash, X_q_dash) *
        [ x[3] - vh * sin( x[1] - θh),
          x[4] - vh * cos( x[1] - θh )]
    
    pg =  x[3] * id_iq[1] + x[4] * id_iq[2] +
        ( X_q_dash - X_d_dash ) *  *( id_iq... )

    id_iq_pg_vh = [ [id_iq, pg, vh ]...; ]
        
    dx =
        Ax_m * x +
        Bx_m * id_iq_pg_vh  +
        Cx_m * ωs_ωref0_vref0_porder0_i
    
end

#-------------------------------------------------------
# system model industrial
#-------------------------------------------------------

function ode_industrial_model_func!(
    dx, x, (industrial_model_pure_states_Idx,
            ode_fun_para), t)

   
    dx_gen  = @view dx[
        industrial_model_pure_states_Idx ]
    
    x_gen   = @view x[
        industrial_model_pure_states_Idx ]
    

    (;Ax_matrix,
     Bx_matrix,
     Cx_matrix,
     industrial_model_pf_para) =
         ode_fun_para

    (; gens_dynamic_id_iq_pg_vh,
     gens_nodes_ωs_τm_vref_porder_view) =
         industrial_model_pf_para
        
    dx_gen .=
        Ax_matrix * x_gen + Bx_matrix *
        [ gens_dynamic_id_iq_pg_vh...; ] +
        Cx_matrix * [ gens_nodes_ωs_τm_vref_porder_view...;]

    return nothing
    
end

function ode_industrial_model_func_no_controllers!(
    dx, x,
    (ode_fun_para,
     industrial_model_pure_states_Idx,
     only_gen),
    t)


    dx_gen  = @view dx[
        industrial_model_pure_states_Idx ]

    x_gen   = @view x[
        industrial_model_pure_states_Idx ]
    
    if only_gen == false
        
        (Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         industrial_model_pf_para) =
             ode_fun_para

        (id_iq_pg_vh,
         ωs_τm_vref_porder) =
             industrial_model_pf_para

        dx_gen .=
            Ax_matrix * x_gen +
            Bx_matrix * [ id_iq_pg_vh...; ] +
            Cx_matrix * [ ωs_τm_vref_porder...;]
        
    else
        
        (Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix,
         industrial_model_pf_para) = ode_fun_para

        (id_iq_pg_vh,
         ωs_τm_vref_porder,
         τm_vf) =
             industrial_model_pf_para

        dx_gen .= Ax_matrix * x_gen +
            Bx_matrix * [ id_iq_pg_vh...; ] +
            Cx_matrix * [ ωs_τm_vref_porder...; ] +
            Ax_τm_vf_matrix * [ τm_vf...; ]

    end
    

    return nothing
    
end

#-------------------------------------------------------
# system  model im
#-------------------------------------------------------

function system_flat_agg_ode_im_model_func!(
    dx,
    x,
    sim_fun_system_ode_flat_agg_para, t;
    sim_fun_kwd_para =
        sim_fun_system_kwd_flat_agg_para )

     (; vec_Ax_views,
      vec_Bx_views,
      vec_Cx_views,
      
      gen_nodes_ra_Xd_dash_Xq_dash_view,      
      stateDiffCache_gens,
      sim_state_x0_gens,
      gens_nodes_collection ,
      gens_sim_fun_gen_para_Idxs,
      
      sim_fun_system_ode_flat_agg_para_Idxs,         
      sys_states_Idxs_and_mat_Idxs,         
      gens_nodes_id_iq_pg_vh_idx_in_Idx,
      gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
      gens_nodes_vh_θh_indx_in_Idx ) =
          sim_fun_kwd_para

    (; vh_θh_Idx,     
     ωs_ωref0_vref0_porder0_Idx ) =
         sim_fun_system_ode_flat_agg_para_Idxs
    
    f_gens_vh_θh =
        sim_fun_system_ode_flat_agg_para[
            vh_θh_Idx]
    
    f_ωs_ωref0_vref0_porder0 =
        sim_fun_system_ode_flat_agg_para[
            ωs_ωref0_vref0_porder0_Idx ] 
    
    gens_vh_θh =
        [ f_gens_vh_θh[idx]
          for idx in
              gens_nodes_vh_θh_indx_in_Idx ]

    gens_ωs_ωref0_vref0_porder0 =
        [ f_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx ]

    gens_vh_θh_ωs_ωref0_vref0_porder0 =
        [ [a_vh_θh...; a_ωs_ωref0_vref0_porder0...]
          for (a_vh_θh, a_ωs_ωref0_vref0_porder0) in
              zip( gens_vh_θh,
                   gens_ωs_ωref0_vref0_porder0) ]

    
    (;nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs ) =
         sys_states_Idxs_and_mat_Idxs
       
    dx_gens  = [ view(dx, idx)
                 for idx in  nodes_state_Idx ]
    
    x_gens   = [ view(x, idx)
                 for idx in nodes_state_Idx ]
    
    vec_sim_fun_kwd_para = [
        (; Ax_m, Bx_m, Cx_m,
         stateDiffCache,
         sim_state_x0,
         plant_i,
         sim_fun_para_Idxs,
         ra_Xd_dash_Xq_dash_i )
        for ( Ax_m, Bx_m, Cx_m,
              stateDiffCache,
              sim_state_x0,
              plant_i,
              sim_fun_para_Idxs,
              ra_Xd_dash_Xq_dash_i ) in
            zip(vec_Ax_views,vec_Bx_views,vec_Cx_views,
                 stateDiffCache_gens,
                 sim_state_x0_gens,
                 gens_nodes_collection,
                 gens_sim_fun_gen_para_Idxs,
                 gen_nodes_ra_Xd_dash_Xq_dash_view  ) ]

    #---------------------------------------------------
    
    for (a_dx, a_x, a_ode_para, a_sim_fun_kwd_para ) in
        zip( dx_gens, x_gens,
             gens_vh_θh_ωs_ωref0_vref0_porder0,
             vec_sim_fun_kwd_para )

        t_ode_one_im_model_func!(
            a_dx, a_x, a_ode_para, t;
            sim_fun_kwd_para =
                a_sim_fun_kwd_para )
    end

    return nothing
end



function system_ode_im_model_func!(
    dx, x,
    system_ode_para, t;
    sim_fun_kwd_para =
        sim_fun_system_kwd_para )

    (; vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     stateDiffCache_gens,
     sim_state_x0_gens,
     gens_nodes_collection,
     gens_sim_fun_gens_para_Idxs,
     sim_fun_system_para_Idxs,         
     sys_states_Idxs_and_mat_Idxs,
     para_update_gen_Ax_aux) =
         sim_fun_kwd_para

    #--------------------------------------------------    
    
    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs) =
         sys_states_Idxs_and_mat_Idxs

    #--------------------------------------------------

    # dx = similar( sim_state_x0 )
    
    dx_gens  = [ view(dx, idx)
                 for idx in nodes_state_Idx ]
    
    #-----------------------------------------------------

    # x = similar( sim_state_x0 )
    
    x_gens   = [ view(x, idx)
                 for idx in nodes_state_Idx ]

    #-----------------------------------------------------

    gens_ode_para = [
        system_ode_para[idx]
        for idx in sim_fun_system_para_Idxs[1] ]
    
    #-----------------------------------------------------
    
    vec_sim_fun_kwd_para = [
        (; Ax_m, Bx_m, Cx_m,
         stateDiffCache,
         sim_state_x0,
         plant_i,
         sim_fun_para_Idxs,
         ra_Xd_dash_Xq_dash_i )
        for ( Ax_m, Bx_m, Cx_m,
              stateDiffCache,
              sim_state_x0,
              plant_i,
              sim_fun_para_Idxs,
              ra_Xd_dash_Xq_dash_i ) in
            zip( vec_Ax_views,
                 vec_Bx_views,
                 vec_Cx_views,
                 stateDiffCache_gens,
                 sim_state_x0_gens,
                 gens_nodes_collection,
                 gens_sim_fun_gens_para_Idxs,
                 gen_nodes_ra_Xd_dash_Xq_dash_view  ) ]

    #-----------------------------------------------------
    
    for (a_dx, a_x, a_ode_para, a_sim_fun_kwd_para ) in
        zip( dx_gens,
             x_gens,
             gens_ode_para,
             vec_sim_fun_kwd_para )

        t_ode_one_im_model_func!(
            a_dx,
            a_x,
            a_ode_para, t;
            sim_fun_kwd_para =
                a_sim_fun_kwd_para )
    end

    return nothing
end


## 
function system_ode_poi_im_model_func!(
    dx,
    x,
    system_ode_para, t;
    poi_idx = poi_idx, 
    sim_fun_system_kwd_para =
        sim_fun_system_kwd_para_wt_poi )

    (; vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views,
     gen_nodes_ra_Xd_dash_Xq_dash_view,     
     stateDiffCache_gens,
     sim_state_x0_gens,
     gens_nodes_collection,
     gens_sim_fun_gens_para_Idxs,
     sim_fun_system_para_Idxs,         
     states_and_mat_Idxs ) =
         sim_fun_system_kwd_para

    (; vh_θh_Idx,     
     ωs_ωref0_vref0_porder0_Idx ) =
         sim_fun_system_para_Idxs

    f_gens_vh_θh =
        sim_fun_system_para[vh_θh_Idx]
    
    f_ωs_ωref0_vref0_porder0 =
        sim_fun_system_para[ωs_ωref0_vref0_porder0_Idx ] 
    

    gens_vh_θh =
        [ f_gens_vh_θh[idx]
          for idx in
              gens_nodes_vh_θh_indx_in_Idx]

    ωs_ωref0_vref0_porder0 =
        [ f_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx ]

    
    (;nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs) =
         states_and_mat_Idxs
       
    dx_gens  = [ view(dx, idx)
                 for idx in
                     nodes_state_Idx ]
    
    x_gens   = [ view(x, idx)
                 for idx in
                     nodes_state_Idx ]

    stateDiffCache_gens =
        [stateDiffCache_gens[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]
    
    vec_sim_fun_kwd_para = [
        (; Ax_m, Bx_m, Cx_m,
         stateDiffCache,
         sim_state_x0,
         plant_i,
         sim_fun_para_Idxs,
         ra_Xd_dash_Xq_dash_i )
        for ( Ax_m, Bx_m, Cx_m,
              stateDiffCache,
              sim_state_x0,
              plant_i,
              sim_fun_para_Idxs,
              ra_Xd_dash_Xq_dash_i ) in
            zip( vec_Ax_views,
                 vec_Bx_views,
                 vec_Cx_views,
                 stateDiffCache_gens,
                 sim_state_x0_gens,
                 gens_nodes_collection,
                 gen_nodes_ra_Xd_dash_Xq_dash_view  ) ]

    return nothing
end


function ode_im_model_func!(
    dx,
    x,
    ode_fun_para,
    t;
    im_vars_Idx_in_state =
        im_vars_Idx_in_state )
   
    dx_gen  = @view dx[ im_vars_Idx_in_state ]
    
    x_gen   = @view x[ im_vars_Idx_in_state ]
    
    (Ax_matrix,
     Bx_matrix,
     Cx_matrix,
     im_pf_para) = ode_fun_para

    id_iq_pg_vh, ωs_ω_ref_vref_porder = im_pf_para
        
    dx_gen .=
        Ax_matrix * x_gen +
        Bx_matrix * [ id_iq_pg_vh...; ] +
        Cx_matrix * [ ωs_ω_ref_vref_porder...;]

    return nothing
    
end

#-----------------------------------------------------
#  ode_model_func!,
# ode_per_gen_model_func, ode_gens_model_func!
#-----------------------------------------------------


# function sd_ode_one_im_model_func!(
#     dx, x, ode_para, t;
#     ode_per_gen_model_func_kwd_para =
#         ode_per_gen_model_func_kwd_para  )

#     (;
#      Ax_view,
#      Bx_view,
#      Cx_view,
     
#      per_gen_ωs_ωref0_vref0_porder0_Idx,
#      per_gen_id_iq_pg_vh_Idx ) =
#          ode_per_gen_model_func_kwd_para
    
#     ωs_ωref0_vref0_porder0 =
#         per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh[
#             per_gen_ωs_ωref0_vref0_porder0_Idx ]

#     id_iq_pg_vh =
#         per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh[
#             per_gen_id_iq_pg_vh_Idx ]
        
#     dx_gen .=
#         Ax_view * x_gen +
#         Bx_view * id_iq_pg_vh +
#         Cx_view * ωs_ωref0_vref0_porder0

    
#     #--------------------------------------------
    
#     (; Ax_m,
#      Bx_m,
#      Cx_m,
#      stateDiffCache,
#      sim_state_x0,
#      plant_i,
#      sim_fun_para_Idxs,
#      ra_Xd_dash_Xq_dash_i ) = 
#          sim_fun_kwd_para

#     #--------------------------------------------

#     vh_θh_i_Idx,
#      ωs_ωref0_vref0_porder0_i_Idx  =
#          sim_fun_para_Idxs

#     #--------------------------------------------

#     stateDiffCache = get_tmp(stateDiffCache, x)

#     stateDiffCache .= sim_state_x0
    
#     #--------------------------------------------
#     # update Ax, 
#     #--------------------------------------------
    
#     update_a_im_plant_system_matrices!(
#         Ax_m,
#         stateDiffCache,
#         plant_i )    

#     #--------------------------------------------

#     ωs_ωref0_vref0_porder0  =
#         ode_para[ ωs_ωref0_vref0_porder0_i_Idx ] 

#     vh_θh_i =
#         ode_para[ vh_θh_i_Idx ]

#     #--------------------------------------------

#     vh, θh, = vh_θh_i

#     ra, X_d_dash, X_q_dash = ra_Xd_dash_Xq_dash_i
    
#     #--------------------------------------------
    
#     id_iq =
#         invZ_dq(ra, X_d_dash, X_q_dash) *
#         [ x[3] - vh * sin( x[1] - θh),
#           x[4] - vh * cos( x[1] - θh )]
    
#     pg =  x[3] * id_iq[1] + x[4] * id_iq[2] +
#         ( X_q_dash - X_d_dash ) *  *( id_iq... )

#     id_iq_pg_vh = [ [id_iq, pg, vh ]...; ]
        
#     dx .=
#         Ax_m * x +
#         Bx_m * id_iq_pg_vh  +
#         Cx_m * ωs_ωref0_vref0_porder0

#     return nothing
    
# end



function ode_vtf_by_vh_θh_func!(
    dx,
    x,
    flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0,
    t;
    kwd_para  =
        ode_vtf_by_vh_θh_kwd_para )
    
    #----------------------------------------    

    (;
     ode_vtf_by_vh_θh_para_kwd,
     ode_vtf_by_vh_θh_Idxs_kwd ) =
         kwd_para
        
    (;
     gens_nodes_ra_Xd_dash_Xq_dash,
     gens_Ax_update_parameters,
     ode_per_gen_models_func_kwd_paras,
     vtf_gens_fun_kwd_para
      ) = ode_vtf_by_vh_θh_para_kwd
    
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
     flat_δ_ed_dash_eq_dash_Idxs_in_state) =
        ode_vtf_by_vh_θh_Idxs_kwd

    #----------------------------------------

    flat_vh_θh_idx, flat_ωs_ωref0_vref0_porder0_idx =
        flat_vh_θh_flat_ωs_ωref0_vref0_porder0_Idx
    
    #----------------------------------------    

    flat_vh_flat_θh =
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0[
            flat_vh_θh_idx]
    
     flat_ωs_ωref0_vref0_porder0 =
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0[
            flat_ωs_ωref0_vref0_porder0_idx]
        
    #----------------------------------------    

    flat_vh_idx_in_flat_Idx, flat_θh_idx_in_flat_Idx =
        flat_vh_flat_θh_Idx
    
    flat_vh = flat_vh_flat_θh[ flat_vh_idx_in_flat_Idx]

    flat_θh = flat_vh_flat_θh[ flat_θh_idx_in_flat_Idx ]
    
    gens_ωs_ωref0_vref0_porder0 =
        [ flat_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
        
    #----------------------------------------    

    vec_vh_θh = [[ a_vh, a_θh]
                 for ( a_vh, a_θh ) in
                     zip( flat_vh, flat_θh )]
    
    #----------------------------------------
        
    gens_vh = flat_vh[gens_nodes_idx]

    gens_θh = flat_θh[gens_nodes_idx]

    gens_vh_θh = vec_vh_θh[ gens_nodes_idx ]

    #----------------------------------------   
    
    non_gens_vh = flat_vh[ non_gens_nodes_idx ]

    non_gens_θh = flat_θh[non_gens_nodes_idx]

    non_gens_vh_θh = vec_vh_θh[
        non_gens_nodes_idx ]

    #----------------------------------------    
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              each_gens_im_vars_Idx_in_state ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              each_gens_im_vars_Idx_in_state ]
    
    #----------------------------------------

    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
           
    # #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]

    #----------------------------------------
    #----------------------------------------
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state )
    
    im_vars_in_state = x[ im_vars_Idx_in_state ]

    #----------------------------------------

    gens_δ_ed_dash_eq_dash =
        [ x[idx]
         for idx in
             nodes_δ_ed_dash_eq_dash_Idxs ]

    flat_δ_ed_dash_eq_dash =
        x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    #----------------------------------------
    #----------------------------------------

    gens_δ =
        first.( gens_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
         for (a_vh, a_θh,
              a_δ, a_ed_dash, a_eq_dash,
              a_ra, a_X_d_dash, a_X_q_dash ) in
                zip( gens_vh, gens_θh,
                     gens_δ, gens_ed_dash, gens_eq_dash,
                     gens_ra, gens_Xd_dash, gens_Xq_dash )]
    
    #----------------------------------------
    
    for ( dx_gen,
          x_gen,
          a_gen_im_vars_Idx_in_state,
          a_id_iq_pg_vh,
          a_ωs_ωref0_vref0_porder0,
          ode_per_gen_model_func_kwd_para,
          a_Ax_update_para ) in

        zip(list_dx_gens_states_views,
            list_x_gens_states_views,
            each_gens_im_vars_Idx_in_state,
            id_iq_pg_vh,
            gens_ωs_ωref0_vref0_porder0,             
            ode_per_gen_models_func_kwd_paras,
            gens_Ax_update_parameters )

        (; Ae, Be,
         Ke, Te,
         vf_tilade_idx_in_Idx ) =
            a_Ax_update_para

        vf_tilade =
            x[a_gen_im_vars_Idx_in_state][
                vf_tilade_idx_in_Idx]
        
        # a_Ax_view, a_Bx_view, a_Cx_view, _, _ =
        #     ode_per_gen_model_func_kwd_para
        
        a_Ax_view, _, _, _, _ =
            ode_per_gen_model_func_kwd_para
        
        updated_Ax = get_a_im_updated_Ax(
            a_Ax_view, vf_tilade, a_Ax_update_para )

        # a_Ax_view[:,:] .= updated_Ax

        ode_per_gen_model_by_vh_θh_by_part_func!(
            dx_gen, x_gen,
            (a_ωs_ωref0_vref0_porder0,
             a_id_iq_pg_vh,
             updated_Ax
             # a_Ax_view
             ),
            t;
            ode_per_gen_model_func_kwd_para =
                ode_per_gen_model_func_kwd_para )
        
    end
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    vtf_gens_vh_θh_δ_ed_dash_eq_dash =
        [[a_gen_vh_θh; a_gen_δ_ed_dash_eq_dash]
         for (a_gen_vh_θh, a_gen_δ_ed_dash_eq_dash) in
             zip(gens_vh_θh, gens_δ_ed_dash_eq_dash) ]
    
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ed_dash_eq_dash,
            t ;
            vtf_kwd_para =
                vtf_kwd_para )
        
    end
    
    #----------------------------------------    
    # vtf non gen nodes
    #----------------------------------------    

    for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh) in
        zip(
            vtf_dx_non_gens_u_views,
            vtf_x_non_gens_u_views,
            non_gens_vh_θh)

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end

    return nothing
    
    #--------------------------------------------
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    #--------------------------------------------   
    
    # for ( a_gen_im_vars_Idx_in_state,
    #       ode_per_gen_model_func_kwd_para,
    #       a_Ax_update_para ) in

    #     zip( each_gens_im_vars_Idx_in_state, 
    #          ode_per_gen_models_func_kwd_paras,
    #          gens_Ax_update_parameters )
        
    #     a_Ax, a_Bx, a_Cx, _, _ =
    #         ode_per_gen_model_func_kwd_para

    #     (sparse_row_idxs,
    #      sparse_col_idxs,
    #      Ax_sparse_nzvalues ) =
    #          findnz( sparse( Matrix( a_Ax )) )

    #     (; Ae, Be,
    #      Ke, Te,
    #      vf_tilade_idx_in_Idx ) =
    #         a_Ax_update_para
        
    #     γ_idx_in_sp_nzv =
    #         find_V_idx_in_sparse_matrix_IJV(
    #             vf_tilade_idx_in_Idx,
    #             vf_tilade_idx_in_Idx,
    #             sparse_row_idxs,
    #             sparse_col_idxs )
    #     vf_tilade =
    #         x[a_gen_im_vars_Idx_in_state][
    #             vf_tilade_idx_in_Idx]

    #     update_γ  = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

    #     Ax_sparse_nzvalues =
    #         [Ax_sparse_nzvalues[1:γ_idx_in_sp_nzv-1];
    #          [update_γ];
    #          Ax_sparse_nzvalues[γ_idx_in_sp_nzv+1:end] ]
        
    #     updated_Ax =
    #         sparse( sparse_row_idxs,
    #                 sparse_col_idxs,
    #                 Ax_sparse_nzvalues)

    #     a_Ax_view[:,:] .= updated_Ax
    # end

    
end


function ode_vtf_by_flat_vh_flat_θh_func!(
    dx,
    x,
    flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0,
    t;
    kwd_para  =
        ode_vtf_by_vh_θh_kwd_para )
    
    #----------------------------------------    

    (;
     ode_vtf_by_vh_θh_para_kwd,
     ode_vtf_by_vh_θh_Idxs_kwd ) =
         kwd_para
        
    (;
     gens_nodes_ra_Xd_dash_Xq_dash,
     gens_Ax_update_parameters,
     ode_per_gen_models_func_kwd_paras,
     vtf_gens_fun_kwd_para
      ) = ode_vtf_by_vh_θh_para_kwd
    
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
     flat_δ_ed_dash_eq_dash_Idxs_in_state) =
        ode_vtf_by_vh_θh_Idxs_kwd

    #----------------------------------------

    flat_vh_θh_idx, flat_ωs_ωref0_vref0_porder0_idx =
        flat_vh_θh_flat_ωs_ωref0_vref0_porder0_Idx
    
    #----------------------------------------    

    flat_vh_flat_θh =
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0[
            flat_vh_θh_idx]
    
     flat_ωs_ωref0_vref0_porder0 =
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0[
            flat_ωs_ωref0_vref0_porder0_idx]
        
    #----------------------------------------    

    flat_vh_idx_in_flat_Idx, flat_θh_idx_in_flat_Idx =
        flat_vh_flat_θh_Idx
    
    flat_vh =
        flat_vh_flat_θh[
            flat_vh_idx_in_flat_Idx]

    flat_θh =
        flat_vh_flat_θh[
            flat_θh_idx_in_flat_Idx ]
    
    gens_ωs_ωref0_vref0_porder0 =
        [ flat_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
        
    #----------------------------------------    

    vec_vh_θh = [[ a_vh, a_θh]
                 for ( a_vh, a_θh ) in
                     zip( flat_vh, flat_θh )]
    
    #----------------------------------------
        
    gens_vh = flat_vh[gens_nodes_idx]

    gens_θh = flat_θh[gens_nodes_idx]

    gens_vh_θh = vec_vh_θh[ gens_nodes_idx ]

    #----------------------------------------   
    
    non_gens_vh = flat_vh[ non_gens_nodes_idx ]

    non_gens_θh = flat_θh[non_gens_nodes_idx]

    non_gens_vh_θh = vec_vh_θh[
        non_gens_nodes_idx ]

    #----------------------------------------    
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              each_gens_im_vars_Idx_in_state ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              each_gens_im_vars_Idx_in_state ]
    
    #----------------------------------------

    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
           
    # #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]

    #----------------------------------------
    #----------------------------------------
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state )
    
    im_vars_in_state = x[ im_vars_Idx_in_state ]

    #----------------------------------------

    gens_δ_ed_dash_eq_dash =
        [ x[idx]
         for idx in
             nodes_δ_ed_dash_eq_dash_Idxs ]

    flat_δ_ed_dash_eq_dash =
        x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    #----------------------------------------
    #----------------------------------------

    gens_δ =
        first.( gens_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
         for (a_vh, a_θh,
              a_δ, a_ed_dash, a_eq_dash,
              a_ra, a_X_d_dash, a_X_q_dash ) in
                zip( gens_vh, gens_θh,
                     gens_δ, gens_ed_dash, gens_eq_dash,
                     gens_ra, gens_Xd_dash, gens_Xq_dash )]
    
    #----------------------------------------
    
    for ( dx_gen,
          x_gen,
          a_gen_im_vars_Idx_in_state,
          a_id_iq_pg_vh,
          a_ωs_ωref0_vref0_porder0,
          ode_per_gen_model_func_kwd_para,
          a_Ax_update_para ) in

        zip(list_dx_gens_states_views,
            list_x_gens_states_views,
            each_gens_im_vars_Idx_in_state,
            id_iq_pg_vh,
            gens_ωs_ωref0_vref0_porder0,             
            ode_per_gen_models_func_kwd_paras,
            gens_Ax_update_parameters )

        (; Ae, Be,
         Ke, Te,
         vf_tilade_idx_in_Idx ) =
            a_Ax_update_para

        vf_tilade =
            x[a_gen_im_vars_Idx_in_state][
                vf_tilade_idx_in_Idx]
        
        # a_Ax_view, a_Bx_view, a_Cx_view, _, _ =
        #     ode_per_gen_model_func_kwd_para
        
        a_Ax_view, _, _, _, _ =
            ode_per_gen_model_func_kwd_para
        
        updated_Ax = get_a_im_updated_Ax(
            a_Ax_view, vf_tilade, a_Ax_update_para )

        # a_Ax_view[:,:] .= updated_Ax

        ode_per_gen_model_by_vh_θh_by_part_func!(
            dx_gen, x_gen,
            (a_ωs_ωref0_vref0_porder0,
             a_id_iq_pg_vh,
             updated_Ax
             # a_Ax_view
             ),
            t;
            ode_per_gen_model_func_kwd_para =
                ode_per_gen_model_func_kwd_para )
        
    end
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    vtf_gens_vh_θh_δ_ed_dash_eq_dash =
        [[a_gen_vh_θh; a_gen_δ_ed_dash_eq_dash]
         for (a_gen_vh_θh, a_gen_δ_ed_dash_eq_dash) in
             zip(gens_vh_θh, gens_δ_ed_dash_eq_dash) ]
    
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ed_dash_eq_dash,
            t ;
            vtf_kwd_para =
                vtf_kwd_para )
        
    end
    
    #----------------------------------------    
    # vtf non gen nodes
    #----------------------------------------    

    for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh) in
        zip(
            vtf_dx_non_gens_u_views,
            vtf_x_non_gens_u_views,
            non_gens_vh_θh)

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end

    return nothing
    
end



function ode_vtf_by_flat_vh_θh_func!(
    dx,
    x,
    flat_vh_θh_flat_ωs_ωref0_vref0_porder0,
    t;
    kwd_para  =
        ode_vtf_by_vh_θh_kwd_para )
    
    #----------------------------------------    

    (;
     ode_vtf_by_vh_θh_para_kwd,
     ode_vtf_by_vh_θh_Idxs_kwd ) =
         kwd_para
        
    (;
     gens_nodes_ra_Xd_dash_Xq_dash,
     gens_Ax_update_parameters,
     ode_per_gen_models_func_kwd_paras,
     vtf_gens_fun_kwd_para
      ) = ode_vtf_by_vh_θh_para_kwd
    
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
     flat_δ_ed_dash_eq_dash_Idxs_in_state) =
        ode_vtf_by_vh_θh_Idxs_kwd

    #----------------------------------------

    flat_vh_θh_idx, flat_ωs_ωref0_vref0_porder0_idx =
        flat_vh_θh_flat_ωs_ωref0_vref0_porder0_Idx
    
    #----------------------------------------    

    flat_vh_flat_θh =
        flat_vh_θh_flat_ωs_ωref0_vref0_porder0[
            flat_vh_θh_idx]
    
     flat_ωs_ωref0_vref0_porder0 =
        flat_vh_θh_flat_ωs_ωref0_vref0_porder0[
            flat_ωs_ωref0_vref0_porder0_idx]
        
    #----------------------------------------    
   
    flat_vh = flat_vh_flat_θh[ flat_vh_idx_in_Idx ]

    flat_θh = flat_vh_flat_θh[ flat_θh_idx_in_Idx ]
    
    gens_ωs_ωref0_vref0_porder0 =
        [ flat_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
        
    #----------------------------------------    

    vec_vh_θh = [[ a_vh, a_θh]
                 for ( a_vh, a_θh ) in
                     zip( flat_vh, flat_θh )]
    
    #----------------------------------------
        
    gens_vh = flat_vh[gens_nodes_idx]

    gens_θh = flat_θh[gens_nodes_idx]

    gens_vh_θh = vec_vh_θh[ gens_nodes_idx ]

    #----------------------------------------   
    
    non_gens_vh = flat_vh[ non_gens_nodes_idx ]

    non_gens_θh = flat_θh[non_gens_nodes_idx]

    non_gens_vh_θh = vec_vh_θh[
        non_gens_nodes_idx ]

    #----------------------------------------    
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              each_gens_im_vars_Idx_in_state ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              each_gens_im_vars_Idx_in_state ]
    
    #----------------------------------------

    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
           
    # #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]

    #----------------------------------------
    #----------------------------------------
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state )
    
    im_vars_in_state = x[ im_vars_Idx_in_state ]

    #----------------------------------------

    gens_δ_ed_dash_eq_dash =
        [ x[idx]
         for idx in
             nodes_δ_ed_dash_eq_dash_Idxs ]

    flat_δ_ed_dash_eq_dash =
        x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    #----------------------------------------
    #----------------------------------------

    gens_δ =
        first.( gens_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
         for (a_vh, a_θh,
              a_δ, a_ed_dash, a_eq_dash,
              a_ra, a_X_d_dash, a_X_q_dash ) in
                zip( gens_vh, gens_θh,
                     gens_δ, gens_ed_dash, gens_eq_dash,
                     gens_ra, gens_Xd_dash, gens_Xq_dash )]
    
    #----------------------------------------
    
    for ( dx_gen,
          x_gen,
          a_gen_im_vars_Idx_in_state,
          a_id_iq_pg_vh,
          a_ωs_ωref0_vref0_porder0,
          ode_per_gen_model_func_kwd_para,
          a_Ax_update_para ) in

        zip(list_dx_gens_states_views,
            list_x_gens_states_views,
            each_gens_im_vars_Idx_in_state,
            id_iq_pg_vh,
            gens_ωs_ωref0_vref0_porder0,             
            ode_per_gen_models_func_kwd_paras,
            gens_Ax_update_parameters )

        (; Ae, Be,
         Ke, Te,
         vf_tilade_idx_in_Idx ) =
            a_Ax_update_para

        vf_tilade =
            x[a_gen_im_vars_Idx_in_state][
                vf_tilade_idx_in_Idx]
        
        # a_Ax_view, a_Bx_view, a_Cx_view, _, _ =
        #     ode_per_gen_model_func_kwd_para
        
        a_Ax_view, _, _, _, _ =
            ode_per_gen_model_func_kwd_para
        
        updated_Ax = get_a_im_updated_Ax(
            a_Ax_view, vf_tilade, a_Ax_update_para )

        # a_Ax_view[:,:] .= updated_Ax

        ode_per_gen_model_by_vh_θh_by_part_func!(
            dx_gen, x_gen,
            (a_ωs_ωref0_vref0_porder0,
             a_id_iq_pg_vh,
             updated_Ax
             # a_Ax_view
             ),
            t;
            ode_per_gen_model_func_kwd_para =
                ode_per_gen_model_func_kwd_para )
        
    end
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    vtf_gens_vh_θh_δ_ed_dash_eq_dash =
        [[a_gen_vh_θh; a_gen_δ_ed_dash_eq_dash]
         for (a_gen_vh_θh, a_gen_δ_ed_dash_eq_dash) in
             zip(gens_vh_θh, gens_δ_ed_dash_eq_dash) ]
    
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ed_dash_eq_dash,
            t ;
            vtf_kwd_para =
                vtf_kwd_para )
        
    end
    
    #----------------------------------------    
    # vtf non gen nodes
    #----------------------------------------    

    for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh) in
        zip(
            vtf_dx_non_gens_u_views,
            vtf_x_non_gens_u_views,
            non_gens_vh_θh)

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end

    return nothing
    
end



function ode_model_func!(
    dx,
    x,
    ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    t ;
    ode_model_func_kwd_para =
        ode_model_func_kwd_para )

    (Ax,
     Bx,
     Cx,
     ωs_ωref0_vref0_porder0_Idx,
     id_iq_pg_vh_Idx) =
        ode_model_func_kwd_para
    
    ωs_ωref0_vref0_porder0 =
        ωs_ωref0_vref0_porder0_id_iq_pg_vh[
            ωs_ωref0_vref0_porder0_Idx ]

    id_iq_pg_vh =
        ωs_ωref0_vref0_porder0_id_iq_pg_vh[
            id_iq_pg_vh_Idx ]
        
    dx .=
        Ax * x +
        Bx * id_iq_pg_vh +
        Cx * ωs_ωref0_vref0_porder0

    return nothing
    
end


function ode_per_gen_model_func!(
    dx_gen,
    x_gen,
    per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    t ;
    ode_per_gen_model_func_kwd_para =
        ode_per_gen_model_func_kwd_para )

    (;
     Ax_view,
     Bx_view,
     Cx_view,
     
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     per_gen_id_iq_pg_vh_Idx ) =
         ode_per_gen_model_func_kwd_para
    
    ωs_ωref0_vref0_porder0 =
        per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh[
            per_gen_ωs_ωref0_vref0_porder0_Idx ]

    id_iq_pg_vh =
        per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh[
            per_gen_id_iq_pg_vh_Idx ]
        
    dx_gen .=
        Ax_view * x_gen +
        Bx_view * id_iq_pg_vh +
        Cx_view * ωs_ωref0_vref0_porder0

    return nothing
    
end


function ode_per_gen_model_by_vh_θh_by_part_func!(
    dx_gen,
    x_gen,
    (ωs_ωref0_vref0_porder0,
    id_iq_pg_vh,
    updated_Ax),
    t;
    ode_per_gen_model_func_kwd_para =
        ode_per_gen_model_func_kwd_para )

         Ax_view, Bx_view, Cx_view, _, _ =
         ode_per_gen_model_func_kwd_para
        
    dx_gen .=
        updated_Ax * x_gen +
        Bx_view * id_iq_pg_vh +
        Cx_view * ωs_ωref0_vref0_porder0

    return nothing
    
end


function ode_gens_model_func!(
    dx_gens,
    x_gens,
    per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    t;
    ode_per_para_model_func_kwd_para =
        ode_per_para_model_func_kwd_para )

    (;
     Ax_matrix,
     Bx_matrix,
     Cx_matrix,
     per_para_ωs_ωref0_vref0_porder0_Idx,
     per_para_id_iq_pg_vh_Idx )  =
         ode_per_para_model_func_kwd_para
    
    ωs_ωref0_vref0_porder0 =
        per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh[
            per_para_ωs_ωref0_vref0_porder0_Idx ]

    id_iq_pg_vh =
        per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh[
            per_para_id_iq_pg_vh_Idx ]
        
    dx_gens .=
        Ax_matrix * x_gens +
        Bx_matrix * id_iq_pg_vh +
        Cx_matrix * ωs_ωref0_vref0_porder0

    return nothing
    
end



#-----------------------------------------------------
######################################################
#-----------------------------------------------------
#  Algebraic models
#-----------------------------------------------------
######################################################
#-----------------------------------------------------

#---------------------------------------------------
# im
#---------------------------------------------------

function non_gens_im_voltage_terminal_func!(
    dx, x, im_para_aux_inputs, t)

    im_nodes_voltage_Idx, _, _, _, _,_, _, _, nodes_pf_U_view =
        im_para_aux_inputs

    _, non_gens_idx, _, nodes_u_Idx_in_ranges  =
        im_nodes_voltage_Idx         
    
    dx_non_gen_ur_ui_views = [
        view( dx, idx ) for idx in
            nodes_u_Idx_in_ranges[ non_gens_idx ] ]
    
    x_non_gen_ur_ui_views   = [
        view( x,  idx ) for idx in
            nodes_u_Idx_in_ranges[ non_gens_idx ] ]
    
    for ( du, u, u_pf ) in zip(
        dx_non_gen_ur_ui_views ,
        x_non_gen_ur_ui_views,
        nodes_pf_U_view[ non_gens_idx ]  )

        a_non_gen_im_voltage_terminal_func!(
            du, u, u_pf, t)
    end
    
end


function gens_im_voltage_terminal_func!(
    dx, x, im_para_aux_inputs, t)

    im_nodes_voltage_Idx, _, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, _,_, gen_nodes_δ_ω_ed_dash_eq_dash_views, gens_vh_θh_view, nodes_pf_U_view = im_para_aux_inputs

    (gens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges) =
         im_nodes_voltage_Idx         
    
    dx_gen_ur_ui_views =
        [ view( dx, idx )
          for idx in
              nodes_u_Idx_in_ranges[ gens_idx] ]
    
    x_gen_ur_ui_views =
        [ view( x,  idx )
          for idx in
              nodes_u_Idx_in_ranges[ gens_idx] ]
    
    for (du, u, gens_vh_θh,
         gen_nodes_δ_ω_ed_dash_eq_dash,
         gen_nodes_ra_Xd_dash_Xq_dash) in
        zip(
            dx_gen_ur_ui_views,
            x_gen_ur_ui_views,
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )

        a_gen_im_voltage_terminal_func!(
            du, u, (gens_vh_θh,
                    gen_nodes_δ_ω_ed_dash_eq_dash,
                    gen_nodes_ra_Xd_dash_Xq_dash), t) 

    end
    
end


function a_gen_im_voltage_terminal_func!(
    dx, x, (gen_vh_θh,
            gen_node_δ_ω_ed_dash_eq_dash,
            gen_node_ra_Xd_dash_Xq_dash), t) 
            
    (δ,
     ω,
     ed_dash,
     eq_dash) =
        gen_node_δ_ω_ed_dash_eq_dash

    id_iq = get_dynamic_idq_vhθh(
        gen_vh_θh...,
        gen_node_δ_ω_ed_dash_eq_dash...,
        gen_node_ra_Xd_dash_Xq_dash... )

    zdq  = Z_dq( gen_node_ra_Xd_dash_Xq_dash... )

    dx .= [ sin(δ) cos(δ); -cos(δ) sin(δ)] *
        ( [ed_dash, eq_dash] - zdq * id_iq ) - x 
        
    return nothing

end


function a_non_gen_im_voltage_terminal_func!(
    dx, x, u_pf, t)

    u_r, u_i = u_pf   

    dx[1] = u_r  - x[1]
    dx[2] = u_i  - x[2]
    
    return nothing

end

function im_all_nodes_voltage_terminal_func!(
    dx, x, (nodes_pf_U_view,
            nodes_u_Idx_in_ranges), t)
        
    dx_non_gen_ur_ui_views = [
        view(dx, idx)
        for idx in
            nodes_u_Idx_in_ranges ]
    
    x_non_gen_ur_ui_views  = [
        view(x, idx)
        for idx in
            nodes_u_Idx_in_ranges ]
   
    for  ( dx_pf, x_pf, u_pf ) in zip(
        dx_non_gen_ur_ui_views,
        x_non_gen_ur_ui_views,
        nodes_pf_U_view )

        algebraic_industrial_model_func!(
            dx_pf, x_pf, u_pf, t)
        
    end
    
end

#---------------------------------------------------

function non_gens_industrial_model_voltage_terminal_func!(
    dx, x, industrial_model_para_aux_inputs, t)

    # industrial_model_nodes_voltage_Idx, _, _, _, _,_, _, _, nodes_pf_U_view = industrial_model_para_aux_inputs


    (;industrial_model_nodes_voltage_Idx,
     industrial_model_pure_states_Idx,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     industrial_model_pure_states_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view) =
         industrial_model_para_aux_inputs
    
     _, non_gens_idx, _, nodes_u_Idx_in_ranges  = industrial_model_nodes_voltage_Idx         
    
    dx_non_gen_ur_ui_views = [
        view( dx, idx )
        for idx in  nodes_u_Idx_in_ranges[ non_gens_idx ] ]
    
    x_non_gen_ur_ui_views   = [
        view( x,  idx )
        for idx in nodes_u_Idx_in_ranges[ non_gens_idx ] ]
    
    for ( du, u, u_pf ) in zip(
        dx_non_gen_ur_ui_views ,
        x_non_gen_ur_ui_views,
        nodes_pf_U_view[ non_gens_idx ]  )

        a_non_gen_industrial_model_voltage_terminal_func!(
            du, u, u_pf, t)
    end
    
end


function gens_industrial_model_voltage_terminal_func!(
    dx, x, industrial_model_para_aux_inputs, t)

    # industrial_model_nodes_voltage_Idx, _, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, _,_, gen_nodes_δ_ω_ed_dash_eq_dash_views, gens_vh_θh_view, nodes_pf_U_view = industrial_model_para_aux_inputs

    (;industrial_model_nodes_voltage_Idx,
     industrial_model_pure_states_Idx,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     industrial_model_pure_states_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view) =
         industrial_model_para_aux_inputs
    
    (gens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges)  =
         industrial_model_nodes_voltage_Idx         
    
    dx_gen_ur_ui_views =
        [ view( dx, idx )
          for idx in
              nodes_u_Idx_in_ranges[gens_idx] ]
    
    x_gen_ur_ui_views  =
        [ view( x,  idx )
          for idx in
              nodes_u_Idx_in_ranges[gens_idx] ]
    
    for (du, u,
         gens_vh_θh,
         gen_nodes_δ_ω_ed_dash_eq_dash,
         gen_nodes_ra_Xd_dash_Xq_dash) in
        zip( dx_gen_ur_ui_views,
             x_gen_ur_ui_views,
             gens_vh_θh_view,
             gen_nodes_δ_ω_ed_dash_eq_dash_views,
             gen_nodes_ra_Xd_dash_Xq_dash_view )

        a_gen_industrial_model_voltage_terminal_func!(
            du, u,
            (gens_vh_θh,
             gen_nodes_δ_ω_ed_dash_eq_dash,
             gen_nodes_ra_Xd_dash_Xq_dash), t) 

    end
    
end


function a_gen_industrial_model_voltage_terminal_func!(
    dx,
    x,
    (gen_vh_θh,
     gen_node_δ_ω_ed_dash_eq_dash,
     gen_node_ra_Xd_dash_Xq_dash),
    t) 
            
    (δ,
     ω,
     ed_dash,
     eq_dash) =
         gen_node_δ_ω_ed_dash_eq_dash

    id_iq = get_dynamic_idq_vhθh(
        gen_vh_θh...,
        gen_node_δ_ω_ed_dash_eq_dash...,
        gen_node_ra_Xd_dash_Xq_dash... )

    zdq  =
        Z_dq( gen_node_ra_Xd_dash_Xq_dash... )

    dx .= [ sin(δ) cos(δ); -cos(δ) sin(δ)] *
        ( [ed_dash, eq_dash] - zdq * id_iq ) - x 
        
    return nothing

end


function a_non_gen_industrial_model_voltage_terminal_func!(
    dx, x, u_pf, t)

    u_r, u_i = u_pf   

    dx[1] = u_r  - x[1]
    dx[2] = u_i  - x[2]
    
    return nothing

end


function algebraic_industrial_model_func!(
    dx, x, algb_fun_para, t)

    u_r, u_i = algb_fun_para    

    dx[1] = u_r  - x[1]
    dx[2] = u_i  - x[2]
    
    return nothing

end


function industrial_model_all_nodes_voltage_terminal_func!(
    dx, x, (nodes_pf_U_view,
            nodes_u_Idx_in_ranges), t)
        
    dx_non_gen_ur_ui_views = [
        view(dx, idx)
        for idx in
            nodes_u_Idx_in_ranges ]
    
    x_non_gen_ur_ui_views  = [
        view(x, idx)
        for idx in
            nodes_u_Idx_in_ranges ]
   
    for  ( dx_pf, x_pf, u_pf ) in
        zip( dx_non_gen_ur_ui_views,
             x_non_gen_ur_ui_views,
             nodes_pf_U_view )

        algebraic_industrial_model_func!(
            dx_pf, x_pf, u_pf, t)
        
    end
    
end



function a_voltage_terminal_by_vh_θh_func!(
    alg_dx,
    alg_x,
    a_vh_θh,
    t )

    a_vh, a_θh = a_vh_θh

    @. alg_dx =
        [ a_vh * cos(a_θh), a_vh * sin(a_θh)] - alg_x
    
    return nothing

end


function a_voltage_terminal_by_ur_ui_func!(
    alg_dx,
    alg_x,
    a_ur_ui, t)
    
    @. alg_dx = a_ur_ui - alg_x
    
    return nothing

end


function a_non_gen_voltage_terminal_ur_ui_by_vh_θh_func!(
    alg_dx_non_gen,
    alg_x_non_gen,
    non_gen_vh_θh,
    t)

    a_vh, a_θh = non_gen_vh_θh

   @. alg_dx_non_gen =
         [ a_vh * cos(a_θh), a_vh * sin(a_θh) ] -
         alg_x_non_gen

    
    return nothing

end


function a_non_gen_voltage_terminal_by_vh_θh_func!(
    alg_dx_non_gen,
    alg_x_non_gen,
    non_gen_vh_θh,
    t)

    a_vh, a_θh = non_gen_vh_θh

   # @. alg_dx_non_gen =
   #       [ a_vh * cos(a_θh), a_vh * sin(a_θh) ] -
   #       alg_x_non_gen


   @. alg_dx_non_gen = non_gen_vh_θh - alg_x_non_gen
    
    return nothing

end


function a_non_gen_voltage_terminal_by_ur_ui_func!(
    alg_dx_non_gen,
    alg_x_non_gen,
    non_gen_ur_ui, t)

    @. alg_dx_non_gen = non_gen_ur_ui - alg_x_non_gen
    
    return nothing

end


function a_gen_voltage_terminal_by_vh_θh_by_part_func!(
    alg_dx_gen,
    alg_x_gen,
    gen_vh_θh,
    gen_node_δ_ed_dash_eq_dash,
    t ;
    vtf_kwd_para =
        vtf_kwd_para )

    (;
     gen_ra_Xd_dash_Xq_dash,
     a_gen_vtf_vh_θh_Idx,
     a_gen_vtf_δ_ed_dash_eq_dash_Idx ) =
         vtf_kwd_para

    #----------------------------------------    

    ra, Xd_dash, Xq_dash = gen_ra_Xd_dash_Xq_dash

    vh, θh = gen_vh_θh

    δ, ed_dash, eq_dash = gen_node_δ_ed_dash_eq_dash
    
    #----------------------------------------
    
    id_iq = get_a_gen_dyn_idq(
            vh, θh,
            δ, ed_dash, eq_dash, ra,
            Xd_dash, Xq_dash )
    
    zdq  = Z_dq( gen_ra_Xd_dash_Xq_dash... )

    α = [ sin(δ) cos(δ); -cos(δ) sin(δ)]
    
    alg_dx_gen .=
        α * ([ed_dash, eq_dash] - zdq * id_iq) - alg_x_gen

    
    # alg_dx_gen .=
    #     gen_vh_θh - alg_x_gen
    
    #----------------------------------------
    
    """

    #----------------------------------------    

    ra, Xd_dash, Xq_dash = gen_ra_Xd_dash_Xq_dash

    vh, θh = gen_vh_θh

    δ, ed_dash, eq_dash = gen_node_δ_ed_dash_eq_dash
    
    #----------------------------------------
    
    id_iq = get_a_gen_dyn_idq(
            vh, θh,
            δ, ed_dash, eq_dash, ra,
            Xd_dash, Xq_dash )
    
    zdq  = Z_dq( gen_ra_Xd_dash_Xq_dash... )

    α = [ sin(δ) cos(δ); -cos(δ) sin(δ)]
    
    alg_dx_gen .=
        α * ([ed_dash, eq_dash] - zdq * id_iq) - alg_x_gen
        


    # See page 160, Sauer


   # id_iq = get_dynamic_idq_vhθh(
   #      gen_vh_θh...,
   #      gen_node_δ_ω_ed_dash_eq_dash...,
   #      gen_ra_Xd_dash_Xq_dash... )

   #  zdq  = Z_dq( gen_ra_Xd_dash_Xq_dash... )

   #  dx .= [ sin(δ) cos(δ); -cos(δ) sin(δ)] *
   #      ( [ed_dash, eq_dash] - zdq * id_iq ) - x 
 
    
    vd_vq = [ed_dash, eq_dash] -  zdq * id_iq

    @. alg_dx_gen = 
        [ sin(δ) cos(δ); -cos(δ) sin(δ)] * vd_vq -
         alg_x_gen
    
    #----------------------------------------    

    a_vh, a_θh = gen_vh_θh

    δ, ed_dash, eq_dash = gen_node_δ_ed_dash_eq_dash

    ra, Xd_dash, Xq_dash = gen_ra_Xd_dash_Xq_dash
    
    #----------------------------------------    

    δ = first(gen_node_δ_ed_dash_eq_dash)
    
    ed_dash = second( gen_node_δ_ed_dash_eq_dash )
    
    eq_dash = third( gen_node_δ_ed_dash_eq_dash )
    
    #----------------------------------------    

    ra =
        first( gen_ra_Xd_dash_Xq_dash )

    Xd_dash =
        second(gen_ra_Xd_dash_Xq_dash )

    Xq_dash =
        third( gen_ra_Xd_dash_Xq_dash )        

    #----------------------------------------    

     @. alg_dx_gen =
         [a_vh * cos(a_θh),a_vh * sin(a_θh)] - alg_x_gen

    """
     
    return nothing

end




function a_gen_voltage_terminal_by_vh_θh_func!(
    alg_dx_gen,
    alg_x_gen,
    gen_vh_θh_δ_ed_dash_eq_dash,
    t ;
    vtf_kwd_para =
        vtf_kwd_para )

    (;
     gen_ra_Xd_dash_Xq_dash,
     a_gen_vtf_vh_θh_Idx,
     a_gen_vtf_δ_ed_dash_eq_dash_Idx ) =
         vtf_kwd_para

    #----------------------------------------    

    ra, Xd_dash, Xq_dash = gen_ra_Xd_dash_Xq_dash
    
    gen_vh_θh = gen_vh_θh_δ_ed_dash_eq_dash[
        a_gen_vtf_vh_θh_Idx]

    vh, θh = gen_vh_θh

    gen_node_δ_ed_dash_eq_dash =
        gen_vh_θh_δ_ed_dash_eq_dash[
            a_gen_vtf_δ_ed_dash_eq_dash_Idx]

    δ, ed_dash, eq_dash = gen_node_δ_ed_dash_eq_dash
    
    #----------------------------------------

    # See page 160, Sauer


   # id_iq = get_dynamic_idq_vhθh(
   #      gen_vh_θh...,
   #      gen_node_δ_ω_ed_dash_eq_dash...,
   #      gen_ra_Xd_dash_Xq_dash... )

   #  zdq  = Z_dq( gen_ra_Xd_dash_Xq_dash... )

   #  dx .= [ sin(δ) cos(δ); -cos(δ) sin(δ)] *
   #      ( [ed_dash, eq_dash] - zdq * id_iq ) - x 
 
    
    id_iq = get_a_gen_dyn_idq(
            vh, θh,
            δ, ed_dash, eq_dash, ra,
            Xd_dash, Xq_dash )
    
    zdq  = Z_dq( gen_ra_Xd_dash_Xq_dash... )

    α = [ sin(δ) cos(δ); -cos(δ) sin(δ)]
    
    alg_dx_gen .=
        α * ([ed_dash, eq_dash] - zdq * id_iq) - alg_x_gen
        

    #----------------------------------------
    
    """
    
    vd_vq = [ed_dash, eq_dash] -  zdq * id_iq

    @. alg_dx_gen = 
        [ sin(δ) cos(δ); -cos(δ) sin(δ)] * vd_vq -
         alg_x_gen
    
    #----------------------------------------    

    a_vh, a_θh = gen_vh_θh

    δ, ed_dash, eq_dash = gen_node_δ_ed_dash_eq_dash

    ra, Xd_dash, Xq_dash = gen_ra_Xd_dash_Xq_dash
    
    #----------------------------------------    

    δ = first(gen_node_δ_ed_dash_eq_dash)
    
    ed_dash = second( gen_node_δ_ed_dash_eq_dash )
    
    eq_dash = third( gen_node_δ_ed_dash_eq_dash )
    
    #----------------------------------------    

    ra =
        first( gen_ra_Xd_dash_Xq_dash )

    Xd_dash =
        second(gen_ra_Xd_dash_Xq_dash )

    Xq_dash =
        third( gen_ra_Xd_dash_Xq_dash )        

    #----------------------------------------    

     @. alg_dx_gen =
         [a_vh * cos(a_θh),a_vh * sin(a_θh)] - alg_x_gen

    """
     
    return nothing

end


function a_gen_voltage_terminal_by_ur_ui_func!(
    alg_dx_gen,
    alg_x_gen,
    gen_ur_ui_δ_ed_dash_eq_dash,
    t
    ; vtf_kwd_para =
        vtf_kwd_para )    

    (;
     gen_ra_Xd_dash_Xq_dash,
     a_gen_vtf_vh_θh_Idx,
     a_gen_vtf_δ_ed_dash_eq_dash_Idx ) =
         vtf_kwd_para

    gen_ur_ui = gen_ur_ui_δ_ed_dash_eq_dash[
        a_gen_vtf_vh_θh_Idx]

    δ_ed_dash_eq_dash =
        gen_ur_ui_δ_ed_dash_eq_dash[
            a_gen_vtf_δ_ed_dash_eq_dash_Idx]

    gen_vh_θh =
        cartesian_to_polar(gen_ur_ui)

    vh, θh = gen_vh_θh
    
    δ, ed_dash, eq_dash

    ra, Xd_dash, Xq_dash = gen_ra_Xd_dash_Xq_dash
    
    id_iq = get_a_gen_dyn_idq(
        vh, θh,
        δ, ed_dash, eq_dash,
        ra, Xd_dash, Xq_dash )
    
    zdq  = Z_dq( gen_ra_Xd_dash_Xq_dash... )

    α = [ sin(δ) cos(δ); -cos(δ) sin(δ)]
    
    alg_dx_gen .=
        α *( [ ed_dash, eq_dash ] - zdq * id_iq) - alg_x_gen
    
     # @. alg_dx_gen = gen_ur_ui - alg_x_gen
        
    return nothing

end



#-----------------------------------------------------
#-----------------------------------------------------
#  update functions
#-----------------------------------------------------
#-----------------------------------------------------

function update_im_dynamic_id_iq_pg_vh!(
    gens_dynamic_id_iq_pg_vh_by_vhθh,
    stateDiffCache,
    im_para_aux_inputs )

    (;im_nodes_voltage_Idx,
     im_vars_Idx_in_state,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     im_vars_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view) =
         im_para_aux_inputs

    (;gens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges)  =
        im_nodes_voltage_Idx

    im_vars_view_in_state[:] .=
        get_gen_nodes_im_vars_from_state(
            stateDiffCache,
            im_vars_Idx_in_state )

    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        im_vars_view_in_state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    # gen_nodes_δ_ω_ed_dash_eq_dash_views .= industrial_model_pure_states_view_in_state[ nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state ] 
    
    gens_vh_θh_view[:] .= get_gens_vh_θh( nodes_pf_U_view, gens_idx )
    
    gens_dynamic_id_iq_pg_vh_by_vhθh .=
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )
    
    return nothing
    
end



function update_dynamic_id_iq_pg_vh!(
    gens_dynamic_id_iq_pg_vh_by_vhθh,
    stateDiffCache,
    industrial_model_para_aux_inputs )

    (;industrial_model_nodes_voltage_Idx,
     industrial_model_pure_states_Idx,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     industrial_model_pure_states_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view) =
         industrial_model_para_aux_inputs

    (;gens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges)  =
        industrial_model_nodes_voltage_Idx
        
    industrial_model_pure_states_view_in_state[:] .=
        get_industrial_gen_nodes_pure_states_from_state(
            stateDiffCache,
            industrial_model_pure_states_Idx )

    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        industrial_model_pure_states_view_in_state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    # gen_nodes_δ_ω_ed_dash_eq_dash_views .= industrial_model_pure_states_view_in_state[ nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state ] 
    
    gens_vh_θh_view[:] .=
        get_gens_vh_θh( nodes_pf_U_view, gens_idx )
    
    gens_dynamic_id_iq_pg_vh_by_vhθh .=
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )
    
    return nothing
    
end

function update_dynamic_τm_vf!(
    vec_τm_vf, stateDiffCache,
    industrial_model_para_aux_inputs )

    (;industrial_model_nodes_voltage_Idx,
     industrial_model_pure_states_Idx,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     industrial_model_pure_states_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view) =
         industrial_model_para_aux_inputs

    (;lgens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges)  =
         industrial_model_nodes_voltage_Idx
    
    # gen_nodes_state_views .= get_gen_nodes_state_views(stateDiffCache, gens_nodes_collection )
        
    industrial_model_pure_states_view_in_state[:] .=
        get_industrial_gen_nodes_pure_states_from_state(
            stateDiffCache,
            industrial_model_pure_states_Idx )

    # update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!( gen_nodes_δ_ω_ed_dash_eq_dash_views, industrial_model_pure_states_view_in_state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    gen_nodes_δ_ω_ed_dash_eq_dash_views[:] .=
        industrial_model_pure_states_view_in_state[
            nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state ] 
    
    gens_vh_θh_view[:] .=
        get_gens_vh_θh( nodes_pf_U_view, gens_idx )

    # nodes_u_Idx_in_ranges
    
    vec_τm_vf .=
        get_gens_dynamic_τm_vf(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view  )
    
    return nothing
    
end


#-----------------------------------------------------
#-----------------------------------------------------
# System model functions
#-----------------------------------------------------
#-----------------------------------------------------


function dynamics_im_model!(dx, x, para, t)

    
    (; netd,
     nodes_cb_sw,
     global_pf_param,
     counter_array,
     state,
     stateDiffCache,
     dyn_global_pf_options,
     para_update_gen_Ax_aux,
     para_model_fun,
     im_misc_Idx,
     im_para_aux_inputs,
     im_dyn_pf_up_para,
     im_idq_pf_cal) = para
    
    #--------------------------------------------

    (;nodes_u_Idx,
     nodes_u_Idx_in_ranges,
     im_vars_Idx_in_state,
     each_gens_im_vars_Idx_in_state,
     gens_nodes_collection) =
         im_misc_Idx

    #--------------------------------------------
    
    vec_Ax_views, ode_fun_para = para_model_fun

    (;Ax_matrix,
     Bx_matrix,
     Cx_matrix,
     im_pf_para) = ode_fun_para
    
    (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para
        
    #--------------------------------------------

    pf_net_param, sd_pf_views, mismatch = global_pf_param

    (; working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view) =
         sd_pf_views
    
    #--------------------------------------------    

    im_vars_view_in_state, _, _ = para_update_gen_Ax_aux
    
    #--------------------------------------------    
    #--------------------------------------------    
    # Powerflow
    #--------------------------------------------
    #--------------------------------------------
        
    counter = counter_array[1]

    #--------------------------------------------

    stateDiffCache = get_tmp(stateDiffCache, x)

    if counter == 1
        
        stateDiffCache .= state

        counter_array .+= 1

    end


    #--------------------------------------------

    industrial_dyn_powerflow(
    nd,
    stateDiffCache,
    global_pf_param,
    im_dyn_pf_up_para,
    im_idq_pf_cal;
    dyn_global_pf_options... )

    #--------------------------------------------
    # update  W
    #--------------------------------------------

    # update_dynamic_id_iq_pg_vh!(
    #     gens_dynamic_id_iq_pg_vh_by_vhθh_view,
    #     stateDiffCache,
    #     im_para_aux_inputs )

    update_im_dynamic_id_iq_pg_vh!(
        gens_dynamic_id_iq_pg_vh_by_vhθh_view,
        stateDiffCache,
        im_para_aux_inputs )
    
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    im_vars_view_in_state .=
        get_gen_nodes_im_vars_from_state(
            stateDiffCache,
            im_vars_Idx_in_state )

    
    update_gens_nodes_im_Ax_system_matrices!(
        vec_Ax_views,
        im_vars_view_in_state,
        each_gens_im_vars_Idx_in_state,
        gens_nodes_collection )
    
    #--------------------------------------------
    # ode gens pure state
    #--------------------------------------------

    # ode_im_model_func!(
    #     dx,
    #     x,
    #     ( im_vars_Idx_in_state,
    #       ode_fun_para),
    #     t)
    
    ode_im_model_func!(
        dx, x, ode_fun_para, t;
        im_vars_Idx_in_state =
            im_vars_Idx_in_state )
    
    #--------------------------------------------
    # Another method for all nodes votages
    #--------------------------------------------

    # im_all_nodes_voltage_terminal_func!(
    #     dx,
    #     x,
    #     (nodes_pf_U_view, nodes_u_Idx_in_ranges),
    #     t)
            
    #--------------------------------------------
    # non gens terminal voltages
    #--------------------------------------------

    non_gens_im_voltage_terminal_func!(
        dx, x, im_para_aux_inputs, t)
    
    #--------------------------------------------
    # gens terminal voltages
    #--------------------------------------------

    gens_im_voltage_terminal_func!(
        dx, x, im_para_aux_inputs, t)
    
    #--------------------------------------------
    #--------------------------------------------    
    
    return nothing
    
end


function dynamics_industrial_model!(dx, x, para, t)
    
    (;netd,
     nodes_cb_sw,
     global_pf_param,
     counter_array,
     state,
     stateDiffCache,
     dyn_global_pf_options,
     para_update_gen_Ax_aux,
     para_model_fun,
     industrial_model_misc_Idx,
     industrial_model_para_aux_inputs,
     industrial_model_dyn_pf_up_para,
     industrial_model_idq_pf_cal) =
         para
    
    #--------------------------------------------

    (;nodes_u_Idx,
     nodes_u_Idx_in_ranges,
     industrial_model_pure_states_Idx,
     industrial_model_each_gen_nodes_pure_states_idx_in_state,
     gens_nodes_collection) =
         industrial_model_misc_Idx

    #--------------------------------------------
    
    (;vec_Ax_views,
     ode_fun_para) = para_model_fun

    (; Ax_matrix,
     Bx_matrix,
     Cx_matrix,
     industrial_model_pf_para) =
         ode_fun_para

    (; gens_dynamic_id_iq_pg_vh_by_vhθh,
     gens_nodes_ωs_τm_vref_porder_view) =
         industrial_model_pf_para
        
    #--------------------------------------------

    (;pf_net_param,
     sd_pf_views,
     mismatch) =
         global_pf_param

    (; working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view) =
         sd_pf_views
    
    #--------------------------------------------    

    industrial_model_pure_states_view_in_state, _, _ =
        para_update_gen_Ax_aux
    
    #--------------------------------------------    
    #--------------------------------------------    
    # Powerflow
    #--------------------------------------------
    #--------------------------------------------
        
    counter = counter_array[1]

    #--------------------------------------------

    stateDiffCache =
        get_tmp(stateDiffCache, x)

    if counter == 1
        
        stateDiffCache .= state
    end

    counter_array .+= 1

    #--------------------------------------------

    industrial_dyn_powerflow(
    netd,
    stateDiffCache,
    global_pf_param,
    industrial_model_dyn_pf_up_para,
    industrial_model_idq_pf_cal;
    dyn_global_pf_options... )

    #--------------------------------------------
    # update  W
    #--------------------------------------------

    update_dynamic_id_iq_pg_vh!(
        gens_dynamic_id_iq_pg_vh_by_vhθh,
        stateDiffCache,
        industrial_model_para_aux_inputs )
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    industrial_model_pure_states_view_in_state .=
        get_industrial_gen_nodes_pure_states_from_state(
            stateDiffCache,
            industrial_model_pure_states_Idx )
    
    update_gens_nodes_Ax_system_matrices!(
        vec_Ax_views,
        industrial_model_pure_states_view_in_state,
        industrial_model_each_gen_nodes_pure_states_idx_in_state,
        gens_nodes_collection )

    #--------------------------------------------
    # ode gens pure state
    #--------------------------------------------

    ode_industrial_model_func!(
        dx,
        x,
        (industrial_model_pure_states_Idx,
         ode_fun_para),
        t)
    
    #--------------------------------------------
    # Another method for all nodes votages
    #--------------------------------------------

    # industrial_model_all_nodes_voltage_terminal_func!(
    #     dx, x,
    #     (nodes_pf_U_view, nodes_u_Idx_in_ranges),
    #     t)
            
    #--------------------------------------------
    # non gens terminal voltages
    #--------------------------------------------

    non_gens_industrial_model_voltage_terminal_func!(
        dx, x, industrial_model_para_aux_inputs, t)
    
    #--------------------------------------------
    # gens terminal voltages
    #--------------------------------------------

    gens_industrial_model_voltage_terminal_func!(
        dx, x, industrial_model_para_aux_inputs, t)
    
    return nothing
    
end

#-------------------------------------------------------


function dynamics_industrial_model_with_or_no_controller!(
    dx,
    x,
    para,
    t)
    
    (;netd,
     nodes_cb_sw,
     global_pf_param,
     counter_array,
     state,
     stateDiffCache,
     dyn_global_pf_options,
     para_update_gen_Ax_aux,
     para_model_fun,
     industrial_model_misc_Idx,
     industrial_model_para_aux_inputs,
     industrial_model_dyn_pf_up_para,
     industrial_model_idq_pf_cal,
     only_gen) = para
    
    (;nodes_u_Idx,
     nodes_u_Idx_in_ranges,
     industrial_model_pure_states_Idx,
     industrial_model_each_gen_nodes_pure_states_idx_in_state,
     gens_nodes_collection) =
         industrial_model_misc_Idx

    (;pf_net_param,
     sd_pf_views,
     mismatch) =
         global_pf_param

    (; working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view) =
         sd_pf_views
        
    industrial_model_pure_states_view_in_state, _, _ =
        para_update_gen_Ax_aux


    #--------------------------------------------    
    #--------------------------------------------    
    # Powerflow
    #--------------------------------------------
    #--------------------------------------------

    counter = counter_array[1]

    #--------------------------------------------

    stateDiffCache = get_tmp(stateDiffCache, x)

    if counter == 1

        stateDiffCache .= state
    end

    counter_array .+= 1

    #--------------------------------------------
    # Power flow calculation
    #--------------------------------------------

    industrial_dyn_powerflow(
    netd,
    stateDiffCache,
    global_pf_param,
    industrial_model_dyn_pf_up_para,
    industrial_model_idq_pf_cal;
    dyn_global_pf_options... )

    #--------------------------------------------
    # Filtering of pure states from all states 
    #--------------------------------------------    
    
    industrial_model_pure_states_view_in_state .=
        get_industrial_gen_nodes_pure_states_from_state(
            stateDiffCache,
            industrial_model_pure_states_Idx )
    
    if only_gen == false
            
        #--------------------------------------------

        vec_Ax_views, ode_fun_para =
            para_model_fun

        _, _, _, industrial_model_pf_para = ode_fun_para

        #--------------------------------------------
        #--------------------------------------------

        (;gens_dynamic_id_iq_pg_vh_by_vhθh_view,
         gens_nodes_ωs_τm_vref_porder_view) =
             industrial_model_pf_para

         #--------------------------------------------
        # update  W
        #--------------------------------------------

        update_dynamic_id_iq_pg_vh!(
            gens_dynamic_id_iq_pg_vh_by_vhθh_view,
            stateDiffCache,
            industrial_model_para_aux_inputs )

        #--------------------------------------------
        # update Ax, 
        #--------------------------------------------

        update_gens_nodes_Ax_system_matrices!(
            vec_Ax_views,
            industrial_model_pure_states_view_in_state,
            industrial_model_each_gen_nodes_pure_states_idx_in_state,
            gens_nodes_collection; only_gen = only_gen )

        #--------------------------------------------
        # gens pure state
        #--------------------------------------------

        # dx_gen  = @view dx[
        #     industrial_model_pure_states_Idx ]

        # x_gen   = @view x[
        #     industrial_model_pure_states_Idx ]
        
        # ode_industrial_model_func_no_controllers!(
        #     dx_gen,
        #     x_gen,
        #     (ode_fun_para, only_gen),
        #     t)

        ode_industrial_model_func_no_controllers!(
            dx, x,
            (ode_fun_para,
             industrial_model_pure_states_Idx,
             only_gen),
            t)
        
        #--------------------------------------------
        # Another method for all nodes votages
        #--------------------------------------------

        # industrial_model_all_nodes_voltage_terminal_func!(
        #     dx,
        #     x,
        #     (nodes_pf_U_view, nodes_u_Idx_in_ranges)
        #     , t)

        #--------------------------------------------
        # non gens terminal voltages
        #--------------------------------------------


        non_gens_industrial_model_voltage_terminal_func!(
        dx, x, industrial_model_para_aux_inputs, t)

        #--------------------------------------------
        # gens terminal voltages
        #--------------------------------------------

        gens_industrial_model_voltage_terminal_func!(
            dx,
            x,
            industrial_model_para_aux_inputs,
            t)
        
    else
 
        #--------------------------------------------

        vec_Ax_views, ode_fun_para = para_model_fun

        (;Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix,
         industrial_model_pf_para) =
             ode_fun_para

        #--------------------------------------------

        (;gens_dynamic_id_iq_pg_vh_by_vhθh_view,
         gens_nodes_ωs_τm_vref_porder_view,
         vec_τm_vf_views) =
             industrial_model_pf_para

        #--------------------------------------------
        # update  W
        #--------------------------------------------

        update_dynamic_id_iq_pg_vh!(
            gens_dynamic_id_iq_pg_vh_by_vhθh_view,
            stateDiffCache,
            industrial_model_para_aux_inputs )

        #--------------------------------------------
        # update Ax, 
        #--------------------------------------------

        update_gens_nodes_Ax_system_matrices!(
            vec_Ax_views,
            industrial_model_pure_states_view_in_state,
            industrial_model_each_gen_nodes_pure_states_idx_in_state, gens_nodes_collection; only_gen = only_gen  )

        #--------------------------------------------
        # gens pure state
        #--------------------------------------------
        
        # dx_gen  =
        #     @view dx[ industrial_model_pure_states_Idx ]

        # x_gen   =
        #     @view x[ industrial_model_pure_states_Idx ]

        # ode_industrial_model_func_no_controllers!(
        #     dx_gen,
        #     x_gen,
        #     (ode_fun_para, only_gen),
        #     t)


        ode_industrial_model_func_no_controllers!(
            dx, x,
            (ode_fun_para,
             industrial_model_pure_states_Idx,
             only_gen),
            t)
        
        #--------------------------------------------
        # Another method for all nodes votages
        #--------------------------------------------

        # industrial_model_all_nodes_voltage_terminal_func!(
        #     dx,
        #     x,
        #     (nodes_pf_U_view, nodes_u_Idx_in_ranges)
        #     , t)

        #--------------------------------------------
        # non gens terminal voltages
        #--------------------------------------------
        
        non_gens_industrial_model_voltage_terminal_func!(
        dx, x, industrial_model_para_aux_inputs, t)

        #--------------------------------------------
        # gens terminal voltages
        #--------------------------------------------

        gens_industrial_model_voltage_terminal_func!(
            dx,
            x,
            industrial_model_para_aux_inputs,
            t)
        
    end
    
    return nothing
    
end

#-----------------------------------------------------

function dynamics_model!(
    dx,
    x,
    f_gens_vh_θh_ωs_ωref0_vref0_porder0_idq_pg_vh_dyn_pf,
    t
    ; dynamics_model_kwd_para =
        dynamics_model_kwd_para)

    (;
     f_gens_vh_θh_Idx,
     f_gens_ωs_ωref0_vref0_porder0_idq_pg_vh_Idx,
     f_dyn_pf_fun_flat_para_Idx,

     Ax_matrix,
     Bx_matrix,
     Cx_matrix,
     
     sub_ωs_ωref0_vref0_porder0_Idx,
     sub_id_iq_pg_vh_Idx,
     
     gens_nodes_ra_Xd_dash_Xq_dash,     
     a_gen_vtf_vh_θh_Idx,
     a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx,     
     
     vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views,
     
     nodes_u_Idx,
     nodes_u_Idx_in_ranges,
     im_vars_Idx_in_state,
     each_gens_im_vars_Idx_in_state,
     gens_nodes_collection,          

     gens_nodes_u_Idx_in_ranges,
     non_gens_nodes_u_Idx_in_ranges,
     
     stateDiffCache,
     state,
     counter_array,

     # get_net_nodes_type_idxs(netd)
     net_nodes_type_idxs,
     
     # get_dict_n2s_streamlined_idx(netd)    
     n2s_streamlined_idx, 

     states_Idxs_and_mat_Idxs,


     dyn_ωs_ωref0_vref0_porder0_Idx,
     dyn_id_iq_pg_vh_Idx,
     
     each_gens_im_vars_Idx_in_state,
     im_vars_Idx_in_state,
     
     nodes_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges,

     nodes_δ_ω_ed_dash_eq_dash_Idxs,

     para_update_gen_Ax_aux,
     
     a_gen_vtf_vh_θh_Idx,
     a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx

     )  = dynamics_model_kwd_para
    
    #--------------------------------------------

    f_gens_vh_θh =
        f_gens_vh_θh_ωs_ωref0_vref0_porder0_idq_pg_vh_dyn_pf[
            f_gens_vh_θh_Idx]

    f_ωs_ωref0_vref0_porder0_idq_pg_vh =
        f_gens_vh_θh_ωs_ωref0_vref0_porder0_idq_pg_vh_dyn_pf[
            f_gens_ωs_ωref0_vref0_porder0_idq_pg_vh_Idx]

    f_dyn_pf_fun_flat_para =
        f_gens_vh_θh_ωs_ωref0_vref0_porder0_idq_pg_vh_dyn_pf[
            f_dyn_pf_fun_flat_para_Idx]
    
    #--------------------------------------------

    ode_model_func_kwd_para =
        (;
         Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         sub_ωs_ωref0_vref0_porder0_Idx,
         sub_id_iq_pg_vh_Idx ) 

    # f_ωs_ωref0_vref0_porder0_id_iq_pg_vh
    # non_gen_vh_θh
    # gen_vh_θh_δ_ω_ed_dash_eq_dash

    # intg_vh_θh_id_iq
    # integ_param = dyn_pf_fun_flat_para
    # intg_dyn_pf_fun_kwd_para 
             
    #--------------------------------------------

    ode_model_fun_kwd_para =
        (;
         Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         sub_ωs_ωref0_vref0_porder0_Idx,
         sub_id_iq_pg_vh_Idx )
    
    #--------------------------------------------
    
    gens_voltage_terminal_fun_kwd_para = [
        (;
         gen_ra_Xd_dash_Xq_dash,
         a_gen_vtf_vh_θh_Idx,
         a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx )

        for gen_ra_Xd_dash_Xq_dash in
            gens_nodes_ra_Xd_dash_Xq_dash ]
    
    #--------------------------------------------


    vtf_kwd_para =
        (;
         gen_ra_Xd_dash_Xq_dash,
         a_gen_vtf_vh_θh_Idx,
         a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx )
    
    #--------------------------------------------

    (;
     nodes_u_Idx,
     nodes_u_Idx_in_ranges,
     im_vars_Idx_in_state,
     each_gens_im_vars_Idx_in_state,
     gens_nodes_collection ) =
         im_misc_Idx

    #--------------------------------------------
    
    vec_Ax_views, ode_fun_para = para_model_fun

    (;
     Ax_matrix,
     Bx_matrix,
     Cx_matrix,
     im_pf_para) = ode_fun_para

    # gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_ωref0_vref0_porder0_view
    
    (;
     gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para
        
    #--------------------------------------------

    pf_net_param, sd_pf_views, mismatch = global_pf_param

    (;
     working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view) =
         sd_pf_views
    
    #--------------------------------------------    

    im_vars_view_in_state, _, _ = para_update_gen_Ax_aux
    
    #--------------------------------------------    
    #--------------------------------------------    
    # Powerflow
    #--------------------------------------------
    #--------------------------------------------
        
    counter = counter_array[1]

    #--------------------------------------------

    stateDiffCache = get_tmp(stateDiffCache, x)

    if counter == 1
        
        stateDiffCache .= state
    end

    counter_array .+= 1

    #--------------------------------------------

    industrial_dyn_powerflow(
    nd,
    stateDiffCache,
    global_pf_param,
    im_dyn_pf_up_para,
    im_idq_pf_cal;
    dyn_global_pf_options... )

    #--------------------------------------------
    # update  W
    #--------------------------------------------

    # update_dynamic_id_iq_pg_vh!(
    #     gens_dynamic_id_iq_pg_vh_by_vhθh_view,
    #     stateDiffCache,
    #     im_para_aux_inputs )

    update_im_dynamic_id_iq_pg_vh!(
    gens_dynamic_id_iq_pg_vh_by_vhθh_view,
    stateDiffCache,
    im_para_aux_inputs )
    
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------

    # im_vars_view_in_state,
    # im_vars_Idx_in_state,
    # each_gens_im_vars_Idx_in_state,
    
    im_vars_view_in_state .=
        get_gen_nodes_im_vars_from_state(
            stateDiffCache,
            im_vars_Idx_in_state )
    
    # update_gens_nodes_Ax_system_matrices!(
    #     vec_Ax_views,
    #     im_vars_view_in_state,
    #     each_gens_im_vars_Idx_in_state,
    #     gens_nodes_collection )


    update_gens_nodes_im_Ax_system_matrices!(
        vec_Ax_views,
        im_vars_view_in_state,
        each_gens_im_vars_Idx_in_state,
        gens_nodes_collection )
    
    #--------------------------------------------
    # ode gens pure state
    #--------------------------------------------

    # ode_im_model_func!(
    #     dx,
    #     x,
    #     ( im_vars_Idx_in_state,
    #       ode_fun_para),
    #     t)
     
     ode_im_model_func!(
         dx, x, ode_fun_para, t;
         im_vars_Idx_in_state =
             im_vars_Idx_in_state )
     
    #--------------------------------------------
    # Another method for all nodes votages
    #--------------------------------------------

    # im_all_nodes_voltage_terminal_func!(
    #     dx,
    #     x,
    #     (nodes_pf_U_view, nodes_u_Idx_in_ranges),
    #     t)
            
    #--------------------------------------------
    # non gens terminal voltages
    #--------------------------------------------

    non_gens_im_voltage_terminal_func!(
        dx, x, im_para_aux_inputs, t)
    
    #--------------------------------------------
    # gens terminal voltages
    #--------------------------------------------

    gens_im_voltage_terminal_func!(
        dx, x, im_para_aux_inputs, t)
    
    #--------------------------------------------
    #--------------------------------------------    
    
    return nothing
    
end

#-----------------------------------------------------

function dynamics_by_per_gen_ode_model_and_powerflow_func!(
    dx,
    x,
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
    t ;
    kwd_para  =
        kwd_para )


    (;sim_state_x0,
     
     stateDiffCache,
     
     pf_solver,
     
     Pg_Qg_external_control,
     
     counter_array,
     
     dynamics_by_per_gen_kwd_para ) =
        kwd_para

    (;
     
    nodes_state_Idx,

    im_vars_Idx_in_state,

    loc_load_exist,

    gens_nodes_idx,

    non_gens_nodes_idx,

    gens_nodes_ra_Xd_dash_Xq_dash,

    nodes_δ_ω_ed_dash_eq_dash_Idxs,

    disaggretation_idxs,

     pre_pf_idx_and_para ,

    intra_ode_pf_fun_kwd_para,

    post_pf_idxs,

    vtf_para_and_idxs,

    ode_per_gen_models_func_kwd_paras,

     gens_nodes_collection,

     vec_Ax_views
     

     ) =
         dynamics_by_per_gen_kwd_para 

    #----------------------------------------    
    #----------------------------------------    

    (;
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     a_gen_vtf_vh_θh_Idx,
     a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx,
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     per_gen_id_iq_pg_vh_Idx,
     per_para_ωs_ωref0_vref0_porder0_Idx,
     per_para_id_iq_pg_vh_Idx,
     f_ωs_ωref0_vref0_porder0_Idx,
     f_dyn_pf_para_Idx,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs) =
         disaggretation_idxs

    (;
     nodes_u_Idx_in_ranges,         
     dyn_pf_fun_kwd_wll_para_vars_Idxs ) =
         pre_pf_idx_and_para 


    (; pf_alg,
     abstol,
     reltol) =
         pf_solver

    (;
     pf_fun_mismatch,
     intg_dyn_pf_fun_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     
     loc_load_exist,
     
     dyn_pf_fun_flat_para,
     dyn_pf_fun_kwd_wll_para_vars_Idxs,
     nodes_δ_ω_ed_dash_eq_dash_Idxs,
     nodes_u_Idx_in_ranges

     ) = intra_ode_pf_fun_kwd_para 
    

    (;
     intg_vh_Idx,
     intg_θh_Idx,
     intg_id_Idx,
     intg_iq_Idx ) =  post_pf_idxs

    (;
     vtf_gens_fun_kwd_para,
     gens_nodes_u_Idx_in_ranges,
     non_gens_nodes_u_Idx_in_ranges,
     ) =
         vtf_para_and_idxs
    
    #----------------------------------------    
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]
    
    #----------------------------------------    
    #----------------------------------------

    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
    
    #----------------------------------------    
    #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           non_gens_nodes_u_Idx_in_ranges ]
    
    #----------------------------------------    
    #----------------------------------------    

    f_gens_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_ωs_ωref0_vref0_porder0_Idx]
    
    dyn_pf_flat_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_dyn_pf_para_Idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ f_gens_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------
    #----------------------------------------    
    # dyn_pf_fun_kwd_para
    
    pf_fun_kwd_para =
        intg_dyn_pf_fun_kwd_para.pf_fun_kwd_para

    loc_load_exist =
        pf_fun_kwd_para.loc_load_exist    
    
    dyn_pf_fun_kwd_wll_para_vars_Idxs =
        pf_fun_kwd_para.dyn_pf_fun_kwd_wll_para_vars_Idxs

   (;P_gens_dyn_para_Idxs,
     Q_gens_dyn_para_Idxs,
     P_non_gens_dyn_para_Idxs,
     Q_non_gens_dyn_para_Idxs,
     δ_ed_eq_pf_dyn_para_Idxs,
     P_g_loc_load_dyn_para_Idxs,
     Q_g_loc_load_dyn_para_Idxs                   
    ) =
        dyn_pf_fun_kwd_wll_para_vars_Idxs


    #----------------------------------------    

    gens_nodes_δ_ω_ed_dash_eq_dash =
        [ x[idx]
         for idx in
             nodes_δ_ω_ed_dash_eq_dash_Idxs ]
    
     flat_δ_ω_ed_dash_eq_dash =
        [gens_nodes_δ_ω_ed_dash_eq_dash...;]
    
    #----------------------------------------    
    #----------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_ed_dash =
        third.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_eq_dash =
        fourth.( gens_nodes_δ_ω_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.(gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    
    #----------------------------------------    
    
     (Pg_Idx,
      Qg_Idxs,
      Png_Idxs,
      Qng_Idxs,
      Pgll_Idxs, Qgll_Idxs) =
        Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------    
    
    P_gens =
        dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        dyn_pf_flat_para[ Qng_Idxs ]    
    
    if loc_load_exist == true

        P_g_loc_load =
            dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            dyn_pf_flat_para[ Qgll_Idxs ]

        dyn_pf_flat_para = [ P_gens...;
            Q_gens...;            
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ω_ed_dash_eq_dash...;
            P_g_loc_load...;
            Q_g_loc_load...] 
    else

        dyn_pf_flat_para = [
            P_gens...;
            Q_gens...;
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ω_ed_dash_eq_dash... ]     

    end
    
    #----------------------------------------    
     
    pf_sol_integ =
        dyn_intg_intra_ode_powerflow(
            pf_fun_mismatch;
            intg_dyn_pf_fun_kwd_para =
                intg_dyn_pf_fun_kwd_para,
            intg_vh_θh_id_iq =
                intg_vh_θh_id_iq,
            dyn_pf_fun_flat_para =
                dyn_pf_fun_flat_para,
            pf_alg = pf_alg )
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    vh_post_pf =
        pf_sol_integ[intg_vh_Idx]
    
    θh_post_pf =
        pf_sol_integ[intg_θh_Idx]
    
    gens_id_post_pf =
        pf_sol_integ[intg_id_Idx]
    
    gens_iq_post_pf =
        pf_sol_integ[intg_iq_Idx]
    
    #----------------------------------------

    gens_vh_post_pf =
        vh_post_pf[gens_nodes_idx]

    gens_θh_post_pf =
        θh_post_pf[gens_nodes_idx]

    gens_vh_θh_post_pf =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(gens_vh_post_pf ,
                 gens_θh_post_pf)]
    
    #----------------------------------------

    non_gens_vh_post_pf =
        vh_post_pf[non_gens_nodes_idx]

    non_gens_θh_post_pf =
        θh_post_pf[non_gens_nodes_idx]

    non_gens_vh_θh_post_pf =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(non_gens_vh_post_pf ,
                 non_gens_θh_post_pf)]

    #----------------------------------------

    gens_id_iq_post_pf =
        [[a_gen_id, a_gen_iq]
         for (a_gen_id, a_gen_iq) in
             zip(gens_id_post_pf,
                 gens_iq_post_pf)]
    
    id_iq_pg_vh_post_pf =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_post_pf,
            gens_nodes_δ_ω_ed_dash_eq_dash,
            gens_nodes_ra_Xd_dash_Xq_dash )

    #----------------------------------------
    #----------------------------------------

    per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
        get_per_vars_or_paras_to_per_node(
            [ gens_nodes_ωs_ωref0_vref0_porder0,
             id_iq_pg_vh_post_pf ] )
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    # im_pure_states_view_in_state = stateDiffCache[im_vars_Idx_in_state]
    
    # update_gens_nodes_im_Ax_system_matrices!(
    #     vec_Ax_views,
    #     im_vars_Idx_in_state,
    #     nodes_state_Idx,
    #     gens_nodes_collection )
    
    #----------------------------------------
    #----------------------------------------

    for ( dx_gen, x_gen,
            per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
            ode_per_gen_model_func_kwd_para ) in
        zip(list_dx_gens_states_views,
            list_x_gens_states_views,
            per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
            ode_per_gen_models_func_kwd_paras)

        ode_per_gen_model_func!(
            dx_gen,
            x_gen,
            per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
            t
            ;ode_per_gen_model_func_kwd_para =
                ode_per_gen_model_func_kwd_para )

    end
    
    #----------------------------------------
    #----------------------------------------

    vtf_gens_vh_θh_δ_ω_ed_dash_eq_dash =
        get_per_vars_or_paras_to_per_node(
            [ gens_vh_θh_post_pf,
              gens_nodes_δ_ω_ed_dash_eq_dash ] )
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ω_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ω_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ω_ed_dash_eq_dash,
            t
            ; vtf_kwd_para =
                vtf_kwd_para )
    end

    #----------------------------------------    
    # vtf non gen nodes
    #----------------------------------------

    for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh) in
        zip(
            vtf_dx_non_gens_u_views,
            vtf_x_non_gens_u_views,
            non_gens_vh_θh)

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end

    
    #----------------------------------------
    #----------------------------------------
    

    return nothing
                                                   
    
    
end



function dynamics_by_ode_gens_model_and_powerflow_func!(
    dx,
    x,
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
    t ;
    kwd_para  =
        kwd_para  )

    (;
     sim_state_x0,     
     stateDiffCache,
     pf_solver,
     Pg_Qg_external_control,
     counter_array,
     dynamics_by_gens_kwd_para) =
        kwd_para

    
    (;   
     sim_state_x0,

     stateDiffCache,

     nodes_state_Idx,

    im_vars_Idx_in_state,

    loc_load_exist,

    gens_nodes_idx,

    non_gens_nodes_idx,

    gens_nodes_ra_Xd_dash_Xq_dash,

    nodes_δ_ω_ed_dash_eq_dash_Idxs,

    disaggretation_idxs,

    pre_pf_idx_and_para ,

     intra_ode_pf_fun_kwd_para,
     
    post_pf_idxs,

    vtf_para_and_idxs,

     ode_per_para_model_func_kwd_para,

     gens_nodes_collection,

     vec_Ax_views
     
     
     ) =
         dynamics_by_gens_kwd_para 

    #----------------------------------------    
    #----------------------------------------    
    
    (;
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     a_gen_vtf_vh_θh_Idx,
     a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx,
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     per_gen_id_iq_pg_vh_Idx,
     per_para_ωs_ωref0_vref0_porder0_Idx,
     per_para_id_iq_pg_vh_Idx,
     f_ωs_ωref0_vref0_porder0_Idx,
     f_dyn_pf_para_Idx,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs) =
         disaggretation_idxs


    (;
     nodes_u_Idx_in_ranges,         
     dyn_pf_fun_kwd_wll_para_vars_Idxs ) =
         pre_pf_idx_and_para 

    
    (; pf_alg,
     abstol,
     reltol) =
         pf_solver

    (;
     pf_fun_mismatch,
     intg_dyn_pf_fun_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     
     loc_load_exist,
     
     dyn_pf_fun_flat_para,
     dyn_pf_fun_kwd_wll_para_vars_Idxs,
     nodes_δ_ω_ed_dash_eq_dash_Idxs,
     nodes_u_Idx_in_ranges

     ) = intra_ode_pf_fun_kwd_para 
    

    (;
     intg_vh_Idx,
     intg_θh_Idx,
     intg_id_Idx,
     intg_iq_Idx ) =  post_pf_idxs

    (;
     vtf_gens_fun_kwd_para,
     gens_nodes_u_Idx_in_ranges,
     non_gens_nodes_u_Idx_in_ranges,
     ) =
         vtf_para_and_idxs

    (;
     Ax_matrix,
     Bx_matrix,
     Cx_matrix,
     per_para_ωs_ωref0_vref0_porder0_Idx,
     per_para_id_iq_pg_vh_Idx ) =
         ode_per_para_model_func_kwd_para
    
    #----------------------------------------
    #----------------------------------------    
    
    dx_gens_states_view =
        view( dx, im_vars_Idx_in_state ) 

    x_gens_states_view =
        view( x, im_vars_Idx_in_state ) 

    #----------------------------------------    
    #----------------------------------------

    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
    
    #----------------------------------------    
    #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           non_gens_nodes_u_Idx_in_ranges ]
    
    #----------------------------------------    
    #----------------------------------------    

    f_gens_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_ωs_ωref0_vref0_porder0_Idx]
    
    dyn_pf_flat_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_dyn_pf_para_Idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ f_gens_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------
    #----------------------------------------    
    # dyn_pf_fun_kwd_para
    
    pf_fun_kwd_para =
        intg_dyn_pf_fun_kwd_para.pf_fun_kwd_para

    loc_load_exist =
        pf_fun_kwd_para.loc_load_exist    
    
    dyn_pf_fun_kwd_wll_para_vars_Idxs =
        pf_fun_kwd_para.dyn_pf_fun_kwd_wll_para_vars_Idxs

   (;P_gens_dyn_para_Idxs,
     Q_gens_dyn_para_Idxs,
     P_non_gens_dyn_para_Idxs,
     Q_non_gens_dyn_para_Idxs,
     δ_ed_eq_pf_dyn_para_Idxs,
     P_g_loc_load_dyn_para_Idxs,
     Q_g_loc_load_dyn_para_Idxs                   
    ) =
        dyn_pf_fun_kwd_wll_para_vars_Idxs 
    
    #----------------------------------------    
    
    P_gens =
        dyn_pf_flat_para[ P_gens_dyn_para_Idxs ]

    Q_gens =
        dyn_pf_flat_para[ Q_gens_dyn_para_Idxs ]

    P_non_gens  =
        dyn_pf_flat_para[ P_non_gens_dyn_para_Idxs ]

    Q_non_gens = 
        dyn_pf_flat_para[ Q_non_gens_dyn_para_Idxs ]

    flat_δ_ω_ed_dash_eq_dash =
        dyn_pf_flat_para[ δ_ed_eq_pf_dyn_para_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            dyn_pf_flat_para[ P_g_loc_load_dyn_para_Idxs ]

        Q_g_loc_load =
            dyn_pf_flat_para[ Q_g_loc_load_dyn_para_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #----------------------------------------

    stateDiffCache = get_tmp(stateDiffCache, x)

    counter = counter_array[1]
    
    if counter == 1
        
        stateDiffCache .= sim_state_x0

        vh =  [ abs( x_from_xr_xi(
            stateDiffCache[ idx ] ) )
                for idx in
                    nodes_u_Idx_in_ranges ]

        θh =  [ angle( x_from_xr_xi(
            stateDiffCache[ idx ] ) )
                for idx in
                    nodes_u_Idx_in_ranges ]

        gens_vh = vh[gens_nodes_idx]

        gens_θh = θh[gens_nodes_idx]

        counter_array .+= 1

    else
        
        stateDiffCache .= x

        vh =  [ abs( x_from_xr_xi(
            stateDiffCache[ idx ] ) )
                for idx in
                    nodes_u_Idx_in_ranges ]

        θh =  [ angle( x_from_xr_xi(
            stateDiffCache[ idx ] ) )
                for idx in
                    nodes_u_Idx_in_ranges ]

        gens_vh = vh[gens_nodes_idx]

        gens_θh = θh[gens_nodes_idx]
        
        
    end

    
    #----------------------------------------    
    #----------------------------------------
    
    # vh =  [ abs( x_from_xr_xi( x[ idx ] ) )
    #         for idx in
    #             nodes_u_Idx_in_ranges ]

    # θh =  [ angle( x_from_xr_xi( x[ idx ] ) )
    #         for idx in
    #             nodes_u_Idx_in_ranges ]
    
    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]
    
    #----------------------------------------

    gens_nodes_δ_ω_ed_dash_eq_dash =
        [ x[idx]
         for idx in
             nodes_δ_ω_ed_dash_eq_dash_Idxs ]
    
    #----------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_ed_dash =
        third.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_eq_dash =
        fourth.( gens_nodes_δ_ω_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )        

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

    gens_dyn_ph = get_gens_dynamic_ph_by_vhθh_δ_idq(
        gens_vh,
        gens_θh,
        gens_δ,
        gens_i_d,
        gens_i_q )
    
    gens_dyn_qh = get_gens_dynamic_qh_by_vhθh_δ_idq(
        gens_vh,
        gens_θh,
        gens_δ,
        gens_i_d,
        gens_i_q )
    
    #----------------------------------------
    #----------------------------------------
    

    if loc_load_exist == true

        dyn_pf_fun_flat_para = [
            # gens_dyn_ph...;
            # gens_dyn_qh...;
            P_gens...;
            Q_gens...;            
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ω_ed_dash_eq_dash...;
            P_g_loc_load...;
            Q_g_loc_load...]                  
    else

        dyn_pf_fun_flat_para = [
            # gens_dyn_ph...;
            # gens_dyn_qh...;
            P_gens...;
            Q_gens...;            
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ω_ed_dash_eq_dash... ] 
    end

    
    #----------------------------------------    
    #----------------------------------------
    
    vh_θh_id_iq =
        [ vh, θh, gens_i_d, gens_i_q ]
    
    #----------------------------------------        
    #----------------------------------------
    
    pf_sol_integ = NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( ( g, x, p ) ->
                pf_fun_mismatch(
                    g, x, p;
                    intg_dyn_pf_fun_kwd_para =
                        intg_dyn_pf_fun_kwd_para )
                               ),
            vh_θh_id_iq,
            dyn_pf_fun_flat_para ),
        pf_alg )

    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------
    
     vh_post_pf =
         pf_sol_integ.u[ vh_Idx ]
    
     θh_post_pf =
         pf_sol_integ.u[ θh_Idx ]
    
     gens_id_post_pf =
         pf_sol_integ.u[ id_Idx ]
    
     gens_iq_post_pf =
         pf_sol_integ.u[ iq_Idx ]
    
    #----------------------------------------

     gens_vh_post_pf =
         vh_post_pf[gens_nodes_idx]

     gens_θh_post_pf =
         θh_post_pf[ gens_nodes_idx ]

    gens_vh_θh_post_pf =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(gens_vh_post_pf , gens_θh_post_pf)]
    
    #----------------------------------------

     non_gens_vh_post_pf =
         vh_post_pf[ non_gens_nodes_idx ]

     non_gens_θh_post_pf =
         θh_post_pf[ non_gens_nodes_idx ]

    non_gens_vh_θh_post_pf =
        [[a_vh, a_θh]
         for (a_vh, a_θh) in
             zip(non_gens_vh_post_pf ,
                 non_gens_θh_post_pf)]

    #----------------------------------------

    gens_id_iq_post_pf =
        [[a_gen_id, a_gen_iq]
         for (a_gen_id, a_gen_iq) in
             zip(gens_id_post_pf , gens_iq_post_pf)]
    
    id_iq_pg_vh_post_pf =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_post_pf,
            gens_nodes_δ_ω_ed_dash_eq_dash,
            gens_nodes_ra_Xd_dash_Xq_dash )

    #----------------------------------------
    #----------------------------------------
    
    per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
        get_a_flattened_by_per_vars_or_paras(
            [ gens_nodes_ωs_ωref0_vref0_porder0,
              id_iq_pg_vh_post_pf ])

    #----------------------------------------
    
    # or
    
    # per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     [ gens_nodes_ωs_ωref0_vref0_porder0...;
    #       id_iq_pg_vh_post_pf... ]
    
    #----------------------------------------
    #----------------------------------------
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    # im_pure_states_view_in_state = stateDiffCache[im_vars_Idx_in_state]
    
    # update_gens_nodes_im_Ax_system_matrices!(
    #     vec_Ax_views,
    #     im_vars_Idx_in_state,
    #     nodes_state_Idx,
    #     gens_nodes_collection )

        
    #----------------------------------------
    #----------------------------------------
    
    ode_gens_model_func!(
        dx_gens_states_view,
        x_gens_states_view,
        per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
        t ;
        ode_per_para_model_func_kwd_para =
            ode_per_para_model_func_kwd_para )
    
    #----------------------------------------
    #----------------------------------------

    vtf_gens_vh_θh_δ_ω_ed_dash_eq_dash =
        get_per_vars_or_paras_to_per_node(
            [ gens_vh_θh_post_pf,
              gens_nodes_δ_ω_ed_dash_eq_dash ] )
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ω_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ω_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ω_ed_dash_eq_dash,
            t
            ; vtf_kwd_para = vtf_kwd_para )
    end

    #----------------------------------------    
    # vtf non gen nodes
    #----------------------------------------

     for (vtf_non_gen_du,
          vtf_non_gen_u,
          non_gen_vh_θh) in
         zip(vtf_dx_non_gens_u_views,
             vtf_x_non_gens_u_views,
             non_gens_vh_θh )

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end

    
    #----------------------------------------
    #----------------------------------------
    

    return nothing
                                                   
    
    
end


function dynamics_by_per_gen_ode_model_only_func!(
    dx,
    x,
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
    t ;
    kwd_para  =
        kwd_para )

    (;
     vh,
     
     θh,
               
     pf_solver,
     
     Pg_Qg_external_control,
     
     counter_array,
     
     dynamics_by_per_gen_kwd_para ) =
        kwd_para

    
    #----------------------------------------    

    
    (;

     loc_load_exist,
     
     gens_nodes_idx,

     non_gens_nodes_idx,

     each_gens_im_vars_Idx_in_state,
     
     nodes_state_Idx,

     im_vars_Idx_in_state,

     nodes_ur_ui_Idx_in_state,

     nodes_u_Idx_in_ranges,

     gens_nodes_ra_Xd_dash_Xq_dash,

     nodes_δ_ed_dash_eq_dash_Idxs,
     
     post_pf_idxs,

     disaggretation_idxs,

     intra_pf_kwd_para,

     vtf_para_and_idxs,

     gens_nodes_collection,

     vec_Ax_views,

     ode_per_gen_models_func_kwd_paras     

     ) =
         dynamics_by_per_gen_kwd_para 

    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol ) =
         pf_solver
    
    #----------------------------------------    

    # (;
    #  nodes_u_Idx_in_ranges,
     
    #  dyn_pf_fun_kwd_wll_para_vars_Idxs
    #  ) =
    #      pre_pf_idx_and_para 
    
    #----------------------------------------    
    
    (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (;
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     
     a_gen_vtf_vh_θh_Idx,
     
     a_gen_vtf_δ_ed_dash_eq_dash_Idx,
     
          
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     
     per_gen_id_iq_pg_vh_Idx,
     
     
     per_para_ωs_ωref0_vref0_porder0_Idx,
     
     per_para_id_iq_pg_vh_Idx,
     
     
     f_ωs_ωref0_vref0_porder0_Idx,
     
     f_dyn_pf_para_Idx,
     f_dyn_ode_pf_para_Idx,
     
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
     
     ) =
         disaggretation_idxs
    
    #----------------------------------------    

    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 

    #----------------------------------------       
    
    (;
     vtf_gens_fun_kwd_para,
     
     gens_nodes_u_Idx_in_ranges,
     
     non_gens_nodes_u_Idx_in_ranges,
     ) =
         vtf_para_and_idxs
    
    #----------------------------------------        
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]
    
    # #----------------------------------------

    # vtf_dx_gens_u_views = [
    #     view( dx, node_u_Idx )
    #     for node_u_Idx in
    #         gens_nodes_u_Idx_in_ranges ]
    
    # vtf_x_gens_u_views = [
    #     view( x, node_u_Idx )
    #     for  node_u_Idx in
    #        gens_nodes_u_Idx_in_ranges ]
    
    # #----------------------------------------

    # vtf_dx_non_gens_u_views = [
    #     view( dx, node_u_Idx )
    #     for node_u_Idx in
    #         non_gens_nodes_u_Idx_in_ranges ]
    
    # vtf_x_non_gens_u_views = [
    #     view( x, node_u_Idx )
    #     for  node_u_Idx in
    #        non_gens_nodes_u_Idx_in_ranges ]

    #----------------------------------------

    # nodes_ur_ui_view =
    #     view( x, nodes_ur_ui_Idx_in_state )

    # nodes_ur_ui =
    #     x[ nodes_ur_ui_Idx_in_state ]

    #----------------------------------------    
    
    # vh =  [ abs( x_from_xr_xi(
    #     nodes_ur_ui[ idx ] ) )
    #         for idx in
    #             ur_ui_idx_in_Idx ]

    # θh =  [ angle( x_from_xr_xi(
    #     nodes_ur_ui[ idx ] ) )
    #         for idx in
    #             ur_ui_idx_in_Idx ]
    
    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]

    # gens_vh_θh =
    #     [[a_vh, a_θh]
    #      for (a_vh, a_θh) in
    #          zip(gens_vh, gens_θh) ]
        
    # non_gens_vh = vh[non_gens_nodes_idx]

    # non_gens_θh = θh[non_gens_nodes_idx]

    # non_gens_vh_θh =
    #     [[a_vh, a_θh]
    #      for (a_vh, a_θh) in
    #          zip(non_gens_vh, non_gens_θh) ]
    
    #----------------------------------------
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]
        
    non_gens_vh = vh[non_gens_nodes_idx]

    non_gens_θh = θh[non_gens_nodes_idx]    

    
    #----------------------------------------    
    #----------------------------------------    

    f_gens_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_ωs_ωref0_vref0_porder0_Idx]
    
    dyn_pf_flat_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_dyn_ode_pf_para_Idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ f_gens_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    
   
    gens_nodes_δ_ed_dash_eq_dash =
        [ x[idx]
         for idx in
             nodes_δ_ed_dash_eq_dash_Idxs ]
    
     flat_δ_ed_dash_eq_dash =
        [gens_nodes_δ_ed_dash_eq_dash...;]
    
    #----------------------------------------    
    
    Pg_Idx, Qg_Idxs, Png_Idxs, Qng_Idxs, Pgll_Idxs, Qgll_Idxs =
        Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------    
    
    P_gens =
        dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        dyn_pf_flat_para[ Qng_Idxs ]    
    
    if loc_load_exist == true

        P_g_loc_load =
            dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            dyn_pf_flat_para[ Qgll_Idxs ]

        dyn_pf_flat_para = [ P_gens...;
            Q_gens...;            
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ed_dash_eq_dash...;
            P_g_loc_load...;
            Q_g_loc_load...] 
    else

        dyn_pf_flat_para = [
            P_gens...;
            Q_gens...;
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ed_dash_eq_dash... ]                

    end

    #--------------------------------------------
    #--------------------------------------------
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    # im_pure_states_view_in_state =
    #     stateDiffCache[
    #         im_vars_Idx_in_state]
    
    # update_gens_nodes_im_Ax_system_matrices!(
    #     vec_Ax_views,
    #     im_vars_Idx_in_state,
    #     nodes_state_Idx,
    #     gens_nodes_collection )
    
    #----------------------------------------
    #----------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_nodes_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_nodes_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
         for (a_vh,
              a_θh,
              a_δ,
              a_ed_dash,
              a_eq_dash,
              a_ra,
              a_X_d_dash,
              a_X_q_dash ) in
                zip( gens_vh,
                     gens_θh,
                     gens_δ,
                     gens_ed_dash,
                     gens_eq_dash,
                     gens_ra,
                     gens_Xd_dash,
                     gens_Xq_dash ) ]
    
    #----------------------------------------

    
    per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
        [[a_gen_ωs_ωref0_vref0_porder0;
          a_gen_id_iq_pg_vh_post_pf]
         for (a_gen_ωs_ωref0_vref0_porder0,
              a_gen_id_iq_pg_vh_post_pf) in
             zip( gens_nodes_ωs_ωref0_vref0_porder0,
                  id_iq_pg_vh ) ] 

    
    #----------------------------------------
    #----------------------------------------
    
    for ( dx_gen, x_gen,
            per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
            ode_per_gen_model_func_kwd_para ) in
        zip(list_dx_gens_states_views,
            list_x_gens_states_views,
            per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
            ode_per_gen_models_func_kwd_paras)

        ode_per_gen_model_func!(
            dx_gen,
            x_gen,
            per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
            t
            ;ode_per_gen_model_func_kwd_para =
                ode_per_gen_model_func_kwd_para )

    end
    
    #----------------------------------------
    #----------------------------------------

    # vtf_gens_vh_θh_δ_ed_dash_eq_dash =
    #     get_per_vars_or_paras_to_per_node(
    #         [ gens_vh_θh,
    #           gens_nodes_δ_ed_dash_eq_dash ] )
    
    # #----------------------------------------    
    # # vtf gen nodes
    # #----------------------------------------

    # for (vtf_gen_dx, vtf_gen_x,
    #      gen_vh_θh_δ_ed_dash_eq_dash,
    #      vtf_kwd_para) in zip(
    #          vtf_dx_gens_u_views,
    #          vtf_x_gens_u_views,
    #          vtf_gens_vh_θh_δ_ed_dash_eq_dash,
    #          vtf_gens_fun_kwd_para )

    #     a_gen_voltage_terminal_by_vh_θh_func!(
    #         vtf_gen_dx,
    #         vtf_gen_x,
    #         gen_vh_θh_δ_ed_dash_eq_dash,
    #         t
    #         ; vtf_kwd_para =
    #             vtf_kwd_para )
        
    # end

    # #----------------------------------------    
    # # vtf non gen nodes
    # #----------------------------------------

    # for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh) in
    #     zip(
    #         vtf_dx_non_gens_u_views,
    #         vtf_x_non_gens_u_views,
    #         non_gens_vh_θh)

    #      a_non_gen_voltage_terminal_by_vh_θh_func!(
    #          vtf_non_gen_du,
    #          vtf_non_gen_u,
    #          non_gen_vh_θh,
    #          t )
    # end

    return nothing

    
end


function dynamics_by_per_gen_ode_model_and_vtf_only_func!(
    dx,
    x,
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
    t ;
    kwd_para  =
        kwd_para )

    (;
     vh,

     θh,

     pf_solver,
     
     Pg_Qg_external_control,
     
     counter_array,
     
     dynamics_by_per_gen_kwd_para ) =
        kwd_para

    
    #----------------------------------------    

    
    (;

     loc_load_exist,
     
     gens_nodes_idx,

     non_gens_nodes_idx,

     each_gens_im_vars_Idx_in_state,
     
     nodes_state_Idx,

     im_vars_Idx_in_state,

     nodes_ur_ui_Idx_in_state,

     nodes_u_Idx_in_ranges,

     gens_nodes_ra_Xd_dash_Xq_dash,

     nodes_δ_ed_dash_eq_dash_Idxs,
     
     post_pf_idxs,

     disaggretation_idxs,

     intra_pf_kwd_para,

     vtf_para_and_idxs,

     gens_nodes_collection,

     vec_Ax_views,

     ode_per_gen_models_func_kwd_paras     

     ) =
         dynamics_by_per_gen_kwd_para 

    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol ) =
         pf_solver
    
    #----------------------------------------    
    
    (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (;
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     
     a_gen_vtf_vh_θh_Idx,
     
     a_gen_vtf_δ_ed_dash_eq_dash_Idx,
     
          
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     
     per_gen_id_iq_pg_vh_Idx,
     
     
     per_para_ωs_ωref0_vref0_porder0_Idx,
     
     per_para_id_iq_pg_vh_Idx,
     
     
     f_ωs_ωref0_vref0_porder0_Idx,
     
     f_dyn_pf_para_Idx,
     f_dyn_ode_pf_para_Idx,
     
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
     
     ) =
         disaggretation_idxs
    
    #----------------------------------------    

    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 

    #----------------------------------------       
    
    (;
     vtf_gens_fun_kwd_para,
     
     gens_nodes_u_Idx_in_ranges,
     
     non_gens_nodes_u_Idx_in_ranges,
     ) =
         vtf_para_and_idxs
    
    #----------------------------------------        
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]
    
    # #----------------------------------------

    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
    
    # #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]

    #----------------------------------------

    # nodes_ur_ui_view =
    #     view( x, nodes_ur_ui_Idx_in_state )

    # nodes_ur_ui =
    #     x[ nodes_ur_ui_Idx_in_state ]

    #----------------------------------------        
    #----------------------------------------
    
    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]
        
    # non_gens_vh = vh[non_gens_nodes_idx]

    # non_gens_θh = θh[non_gens_nodes_idx]    

    
    # vh =  [ abs( x_from_xr_xi(
    #     nodes_ur_ui[ idx ] ) )
    #         for idx in
    #             ur_ui_idx_in_Idx ]

    # θh =  [ angle( x_from_xr_xi(
    #     nodes_ur_ui[ idx ] ) )
    #         for idx in
    #             ur_ui_idx_in_Idx ]


    # vh =  [ abs( x_from_xr_xi(
    #     x[ idx ] ) )
    #         for idx in
    #             nodes_u_Idx_in_ranges ]

    # θh =  [ angle( x_from_xr_xi(
    #     x[ idx ] ) )
    #         for idx in
    #             nodes_u_Idx_in_ranges ]
    
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]

    gens_vh_θh =
        [[a_vh, a_θh]
         for (a_vh, a_θh) in
             zip(gens_vh, gens_θh) ]
        
    non_gens_vh = vh[non_gens_nodes_idx]

    non_gens_θh = θh[non_gens_nodes_idx]

    non_gens_vh_θh =
        [[a_vh, a_θh]
         for (a_vh, a_θh) in
             zip(non_gens_vh, non_gens_θh) ]
    
    #----------------------------------------    
    #----------------------------------------    

    f_gens_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_ωs_ωref0_vref0_porder0_Idx]
    
    dyn_pf_flat_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_dyn_ode_pf_para_Idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ f_gens_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    
   
    gens_nodes_δ_ed_dash_eq_dash =
        [ x[idx]
         for idx in
             nodes_δ_ed_dash_eq_dash_Idxs ]
    
     flat_δ_ed_dash_eq_dash =
        [gens_nodes_δ_ed_dash_eq_dash...;]
    
    #----------------------------------------    
    
    Pg_Idx, Qg_Idxs, Png_Idxs, Qng_Idxs, Pgll_Idxs, Qgll_Idxs =
        Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------    
    
    P_gens =
        dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        dyn_pf_flat_para[ Qng_Idxs ]    
    
    if loc_load_exist == true

        P_g_loc_load =
            dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            dyn_pf_flat_para[ Qgll_Idxs ]

        dyn_pf_flat_para = [ P_gens...;
            Q_gens...;            
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ed_dash_eq_dash...;
            P_g_loc_load...;
            Q_g_loc_load...] 
    else

        dyn_pf_flat_para = [
            P_gens...;
            Q_gens...;
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ed_dash_eq_dash... ]               

    end

    #--------------------------------------------
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    #--------------------------------------------

    
    # im_pure_states_view_in_state =
    #     stateDiffCache[
    #         im_vars_Idx_in_state]
    
    # update_gens_nodes_im_Ax_system_matrices!(
    #     vec_Ax_views,
    #     im_vars_Idx_in_state,
    #     nodes_state_Idx,
    #     gens_nodes_collection )
    
    #----------------------------------------
    #----------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_nodes_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_nodes_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
         for (a_vh,
              a_θh,
              a_δ,
              a_ed_dash,
              a_eq_dash,
              a_ra,
              a_X_d_dash,
              a_X_q_dash ) in
                zip( gens_vh,
                     gens_θh,
                     gens_δ,
                     gens_ed_dash,
                     gens_eq_dash,
                     gens_ra,
                     gens_Xd_dash,
                     gens_Xq_dash ) ]
    
    #----------------------------------------

    
    per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
        [[a_gen_ωs_ωref0_vref0_porder0;
          a_gen_id_iq_pg_vh_post_pf]
         for (a_gen_ωs_ωref0_vref0_porder0,
              a_gen_id_iq_pg_vh_post_pf) in
             zip( gens_nodes_ωs_ωref0_vref0_porder0,
                  id_iq_pg_vh ) ] 
    
    #----------------------------------------
    #----------------------------------------
    
    for ( dx_gen, x_gen,
            per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
            ode_per_gen_model_func_kwd_para ) in
        zip(list_dx_gens_states_views,
            list_x_gens_states_views,
            per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
            ode_per_gen_models_func_kwd_paras )

        ode_per_gen_model_func!(
            dx_gen,
            x_gen,
            per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
            t
            ;ode_per_gen_model_func_kwd_para =
                ode_per_gen_model_func_kwd_para )

    end
    
    #----------------------------------------
    #----------------------------------------

    vtf_gens_vh_θh_δ_ed_dash_eq_dash =
        get_per_vars_or_paras_to_per_node(
            [ gens_vh_θh,
              gens_nodes_δ_ed_dash_eq_dash ] )
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ed_dash_eq_dash,
            t
            ; vtf_kwd_para =
                vtf_kwd_para )
        
    end
 
    #----------------------------------------    
    # vtf non gen nodes

    for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh) in
        zip(
            vtf_dx_non_gens_u_views,
            vtf_x_non_gens_u_views,
            non_gens_vh_θh)

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end

    return nothing

    
end

 
  
function dynamics_by_per_gen_powerflow_only_func!(
    dx,
    x,
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
    t;
    kwd_para =
        kwd_para )

    (;
     stand_alone_pf,
     
     gens_nodes_δ_ed_dash_eq_dash,
     
     sim_state_x0,
     
     ur_ui_DiffCache,
     
     pf_solver,

     use_nlsolve,
     
     Pg_Qg_external_control,
     
     counter_array,
     
     dynamics_by_per_gen_kwd_para ) =
        kwd_para
    
    #----------------------------------------    
    
    (;

     loc_load_exist,
     
     gens_nodes_idx,

     non_gens_nodes_idx,

     each_gens_im_vars_Idx_in_state,
     
     nodes_state_Idx,

     im_vars_Idx_in_state,

     nodes_ur_ui_Idx_in_state,

     nodes_u_Idx_in_ranges,

     gens_nodes_ra_Xd_dash_Xq_dash,

     nodes_δ_ed_dash_eq_dash_Idxs,
     
     post_pf_idxs,

     disaggretation_idxs,

     intra_pf_kwd_para,

     vtf_para_and_idxs,

     gens_nodes_collection,

     vec_Ax_views,

     ode_per_gen_models_func_kwd_paras     

     ) =
         dynamics_by_per_gen_kwd_para 

    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol ) =
         pf_solver
    
    #----------------------------------------    

    # (;
    #  nodes_u_Idx_in_ranges,
     
    #  dyn_pf_fun_kwd_wll_para_vars_Idxs
    #  ) =
    #      pre_pf_idx_and_para 
    
    #----------------------------------------    
    
    (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (;
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     
     a_gen_vtf_vh_θh_Idx,
     
     a_gen_vtf_δ_ed_dash_eq_dash_Idx,
     
          
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     
     per_gen_id_iq_pg_vh_Idx,
     
     
     per_para_ωs_ωref0_vref0_porder0_Idx,
     
     per_para_id_iq_pg_vh_Idx,
     
     
     f_ωs_ωref0_vref0_porder0_Idx,
     
     f_dyn_pf_para_Idx,
     f_dyn_ode_pf_para_Idx,
     
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
     
     ) =
         disaggretation_idxs
    
    #----------------------------------------    

    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 

    #----------------------------------------       
    
    (;
     vtf_gens_fun_kwd_para,
     
     gens_nodes_u_Idx_in_ranges,
     
     non_gens_nodes_u_Idx_in_ranges,
     ) =
         vtf_para_and_idxs
    
    #----------------------------------------        
    #----------------------------------------

    if stand_alone_pf == true

        vtf_dx_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            ur_ui_idx_in_Idx ]
    

        vtf_x_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            ur_ui_idx_in_Idx ]
         
    #----------------------------------------
   
        first_idx = first(first(ur_ui_idx_in_Idx))

        last_idx = last(last(ur_ui_idx_in_Idx))

        u_r_i = first_idx:last_idx
        
        nodes_ur_ui = x[ u_r_i  ]
            
    else
         
        vtf_dx_gens_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                gens_nodes_u_Idx_in_ranges ]

        vtf_x_gens_u_views = [
            view( x, node_u_Idx )
            for  node_u_Idx in
               gens_nodes_u_Idx_in_ranges ]

        # #----------------------------------------

        vtf_dx_non_gens_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                non_gens_nodes_u_Idx_in_ranges ]

        vtf_x_non_gens_u_views = [
            view( x, node_u_Idx )
            for  node_u_Idx in
                non_gens_nodes_u_Idx_in_ranges ]

        #----------------------------------------        

        vtf_dx_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                nodes_u_Idx_in_ranges ]

        vtf_x_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                nodes_u_Idx_in_ranges ]

        #----------------------------------------        
                
        nodes_ur_ui = x[nodes_ur_ui_Idx_in_state]
        
    end
    
    #----------------------------------------
    #----------------------------------------

    f_gens_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_ωs_ωref0_vref0_porder0_Idx]
    
    dyn_pf_flat_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_dyn_ode_pf_para_Idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ f_gens_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    
    
     flat_δ_ed_dash_eq_dash =
        [gens_nodes_δ_ed_dash_eq_dash...;]
    
    #----------------------------------------    
    
    (;Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs) =
        Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------    
    
    P_gens =
        dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        dyn_pf_flat_para[ Qng_Idxs ]    
    
    if loc_load_exist == true

        P_g_loc_load =
            dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            dyn_pf_flat_para[ Qgll_Idxs ]

        intra_dyn_pf_flat_para = [ P_gens...;
            Q_gens...;            
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ed_dash_eq_dash...;
            P_g_loc_load...;
            Q_g_loc_load...] 
    else

        intra_dyn_pf_flat_para = [
            P_gens...;
            Q_gens...;
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ed_dash_eq_dash... ]   

    end

    #--------------------------------------------
    #--------------------------------------------

    intra_dyn_ur_ui = [nodes_ur_ui...;]

    pf_sol_integ =
        dyn_pf_func!(
            intra_dyn_ur_ui,
            intra_dyn_pf_flat_para;
            ur_ui_DiffCache  =
                ur_ui_DiffCache,
            Pg_Qg_external_control =
                Pg_Qg_external_control,
            intra_pf_kwd_para =
                intra_pf_kwd_para,
            pf_solver =
                pf_solver,
            use_nlsolve =
                use_nlsolve)
    
    """
    
    vh =  [ abs( x_from_xr_xi(
        nodes_ur_ui[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]

    θh =  [ angle( x_from_xr_xi(
        nodes_ur_ui[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]

    gens_vh_θh =
        [[a_vh, a_θh]
         for (a_vh, a_θh) in
             zip(gens_vh, gens_θh) ]
        
    non_gens_vh = vh[non_gens_nodes_idx]

    non_gens_θh = θh[non_gens_nodes_idx]

    non_gens_vh_θh =
        [[a_vh, a_θh]
         for (a_vh, a_θh) in
             zip(non_gens_vh, non_gens_θh) ]
    
    #----------------------------------------
    #----------------------------------------

    ur_ui_DiffCache =
        get_tmp(ur_ui_DiffCache, x)
        
    #----------------------------------------    

    intra_ur_ui_wt_pf_flat_para =
        [ nodes_ur_ui...; dyn_pf_flat_para... ]

    
    #----------------------------------------    
    
    pf_sol_integ =
        intra_dyn_pf_func!(
            intra_ur_ui_wt_pf_flat_para;
            
            ur_ui_DiffCache =
                ur_ui_DiffCache,
            
            Pg_Qg_external_control =
                Pg_Qg_external_control,
            
            intra_pf_kwd_para =
                intra_pf_kwd_para,
            
            pf_solver =
                pf_solver,
            
            use_nlsolve =
                use_nlsolve
        )

    """
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    vh_post_pf =
        pf_sol_integ[vh_Idx]
    
    θh_post_pf =
        pf_sol_integ[θh_Idx]
    
    gens_id_post_pf =
        pf_sol_integ[id_Idx]
    
    gens_iq_post_pf =
        pf_sol_integ[iq_Idx]
    
    #----------------------------------------

    nodes_vh_θh_post_pf =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(vh_post_pf ,
                 θh_post_pf)]
    
    #----------------------------------------

    gens_vh_post_pf =
        vh_post_pf[gens_nodes_idx]

    gens_θh_post_pf =
        θh_post_pf[gens_nodes_idx]

    gens_vh_θh_post_pf =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(gens_vh_post_pf ,
                 gens_θh_post_pf)]
    
    #----------------------------------------

    non_gens_vh_post_pf =
        vh_post_pf[non_gens_nodes_idx]

    non_gens_θh_post_pf =
        θh_post_pf[non_gens_nodes_idx]

    non_gens_vh_θh_post_pf =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(non_gens_vh_post_pf ,
                 non_gens_θh_post_pf)]

    #----------------------------------------

    gens_id_iq_post_pf =
        [[a_gen_id, a_gen_iq]
         for (a_gen_id, a_gen_iq) in
             zip(gens_id_post_pf,
                 gens_iq_post_pf)]

    #----------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_nodes_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_nodes_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh_post_pf =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
         for (a_vh,
              a_θh,
              a_δ,
              a_ed_dash,
              a_eq_dash,
              a_ra,
              a_X_d_dash,
              a_X_q_dash ) in
                zip( gens_vh_post_pf,
                     gens_θh_post_pf,
                     gens_δ,
                     gens_ed_dash,
                     gens_eq_dash,
                     gens_ra,
                     gens_Xd_dash,
                     gens_Xq_dash ) ]
    
    #----------------------------------------
    #----------------------------------------

    per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
        [[a_gen_ωs_ωref0_vref0_porder0;
          a_gen_id_iq_pg_vh_post_pf]
         for (a_gen_ωs_ωref0_vref0_porder0,
              a_gen_id_iq_pg_vh_post_pf) in
             zip( gens_nodes_ωs_ωref0_vref0_porder0,
                  id_iq_pg_vh_post_pf) ] 

    #----------------------------------------    
    # vtf  nodes
    #----------------------------------------

    for ( a_vtf_du, a_vtf_u, a_vh_θh ) in
        zip(
            vtf_dx_u_views,
            vtf_x_u_views,
            nodes_vh_θh_post_pf )

         a_voltage_terminal_by_vh_θh_func!(
              a_vtf_du, a_vtf_u, a_vh_θh, t )
    end

    #----------------------------------------
    #----------------------------------------

    return nothing
    
end
 
  
function dynamics_by_per_gen_steady_state_powerflow_only_func!(
    dx,
    x,
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
    t;
    kwd_para =
        kwd_para )

    (;
     stand_alone_pf,
     
     gens_nodes_δ_ed_dash_eq_dash,
     
     sim_state_x0,

     u_stateDiffCache,
     
     ur_ui_DiffCache,
         
     vh_θh_DiffCache,
     
     pf_solver,

     use_nlsolve,
     
     Pg_Qg_external_control,
     
     counter_array,
     
     dynamics_by_per_gen_kwd_para ) =
        kwd_para
    
    #----------------------------------------    
    
    (;

     loc_load_exist,
     
     gens_nodes_idx,

     non_gens_nodes_idx,

     each_gens_im_vars_Idx_in_state,
     
     nodes_state_Idx,

     im_vars_Idx_in_state,

     nodes_ur_ui_Idx_in_state,

     nodes_u_Idx_in_ranges,

     gens_nodes_ra_Xd_dash_Xq_dash,

     nodes_δ_ed_dash_eq_dash_Idxs,
     
     post_pf_idxs,

     disaggretation_idxs,

     intra_pf_kwd_para,

     vtf_para_and_idxs,

     gens_nodes_collection,

     vec_Ax_views,

     ode_per_gen_models_func_kwd_paras     

     ) =
         dynamics_by_per_gen_kwd_para 

    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol ) =
         pf_solver
    
    #----------------------------------------    
    
    (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (;
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     
     a_gen_vtf_vh_θh_Idx,
     
     a_gen_vtf_δ_ed_dash_eq_dash_Idx,
     
          
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     
     per_gen_id_iq_pg_vh_Idx,
     
     
     per_para_ωs_ωref0_vref0_porder0_Idx,
     
     per_para_id_iq_pg_vh_Idx,
     
     
     f_ωs_ωref0_vref0_porder0_Idx,
     
     f_dyn_pf_para_Idx,
     f_dyn_ode_pf_para_Idx,
     
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
     
     ) =
         disaggretation_idxs
    
    #----------------------------------------    

    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 

    #----------------------------------------       
    
    (;
     vtf_gens_fun_kwd_para,
     
     gens_nodes_u_Idx_in_ranges,
     
     non_gens_nodes_u_Idx_in_ranges,
     ) =
         vtf_para_and_idxs
    
    #----------------------------------------
    #----------------------------------------

    f_gens_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_ωs_ωref0_vref0_porder0_Idx]
    
    dyn_pf_flat_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_dyn_ode_pf_para_Idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ f_gens_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    
    
     flat_δ_ed_dash_eq_dash =
        [ gens_nodes_δ_ed_dash_eq_dash...;]
    
    #----------------------------------------    
    
    Pg_Idx, Qg_Idxs, Png_Idxs, Qng_Idxs, Pgll_Idxs, Qgll_Idxs =
        Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------    
    
    if loc_load_exist == true
    
        P_gens =
            dyn_pf_flat_para[ Pg_Idx ]

        Q_gens =
            dyn_pf_flat_para[ Qg_Idxs ]

        P_non_gens  =
            dyn_pf_flat_para[ Png_Idxs ]

        Q_non_gens = 
            dyn_pf_flat_para[ Qng_Idxs ]    

        P_g_loc_load =
            dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            dyn_pf_flat_para[ Qgll_Idxs ]

        intra_dyn_pf_flat_para = [ P_gens...;
            Q_gens...;            
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ed_dash_eq_dash...;
            P_g_loc_load...;
            Q_g_loc_load...] 
    else
        P_gens =
            dyn_pf_flat_para[ Pg_Idx ]

        Q_gens =
            dyn_pf_flat_para[ Qg_Idxs ]

        P_non_gens  =
            dyn_pf_flat_para[ Png_Idxs ]

        Q_non_gens = 
            dyn_pf_flat_para[ Qng_Idxs ]    

        intra_dyn_pf_flat_para = [
            P_gens...;
            Q_gens...;
            P_non_gens...;
            Q_non_gens...;
            gens_nodes_δ_ed_dash_eq_dash... ]   

    end

    #--------------------------------------------
    #--------------------------------------------
    
    if stand_alone_pf == true

        vtf_dx_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                ur_ui_idx_in_Idx ]
    

        vtf_x_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                ur_ui_idx_in_Idx ]
         
    #----------------------------------------
   
        first_idx = first(first(ur_ui_idx_in_Idx))

        last_idx = last(last(ur_ui_idx_in_Idx))

        state_idx_range_ur_ui = first_idx:last_idx
        
        nodes_ur_ui = x[ state_idx_range_ur_ui  ]
            
    else
         
        vtf_dx_gens_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                gens_nodes_u_Idx_in_ranges ]

        vtf_x_gens_u_views = [
            view( x, node_u_Idx )
            for  node_u_Idx in
               gens_nodes_u_Idx_in_ranges ]

        # #----------------------------------------

        vtf_dx_non_gens_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                non_gens_nodes_u_Idx_in_ranges ]

        vtf_x_non_gens_u_views = [
            view( x, node_u_Idx )
            for  node_u_Idx in
                non_gens_nodes_u_Idx_in_ranges ]

        #----------------------------------------        

        vtf_dx_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                nodes_u_Idx_in_ranges ]

        vtf_x_u_views = [
            view( dx, node_u_Idx )
            for node_u_Idx in
                nodes_u_Idx_in_ranges ]

        #----------------------------------------        
                
        nodes_ur_ui = x[ nodes_ur_ui_Idx_in_state ]
        
    end
    
    #--------------------------------------------
    #--------------------------------------------
    
    u_stateDiffCache = get_tmp(u_stateDiffCache, x)

    counter = counter_array[1]

    if counter == 1
        
        if stand_alone_pf == true

            u_stateDiffCache .= sim_state_x0

            vh =  [ abs( x_from_xr_xi(
                u_stateDiffCache[ idx ] ) )
                    for idx in
                        ur_ui_idx_in_Idx ]

            θh =  [ angle( x_from_xr_xi(
                u_stateDiffCache[ idx ] ) )
                    for idx in
                        ur_ui_idx_in_Idx ]
            
        else

            u_stateDiffCache .=
                sim_state_x0[ nodes_ur_ui_Idx_in_state ]

            vh =  [ abs( x_from_xr_xi(
                u_stateDiffCache[ idx ] ) )
                    for idx in
                        ur_ui_idx_in_Idx ]

            θh =  [ angle( x_from_xr_xi(
                u_stateDiffCache[ idx ] ) )
                    for idx in
                        ur_ui_idx_in_Idx ]            
            
        end

        intra_dyn_ur_ui = [u_stateDiffCache...;]

        intra_vh_θh = [ vh; θh ]

        counter_array .+= 1
        
        
    else
        
        if stand_alone_pf == true

            u_stateDiffCache .= x

            vh =  [ abs( x_from_xr_xi(
                u_stateDiffCache[ idx ] ) )
                    for idx in
                        ur_ui_idx_in_Idx ]

            θh =  [ angle( x_from_xr_xi(
                u_stateDiffCache[ idx ] ) )
                    for idx in
                        ur_ui_idx_in_Idx ]
            
        else


            u_stateDiffCache .= x[ nodes_ur_ui_Idx_in_state ]

            vh =  [ abs( x_from_xr_xi(
                u_stateDiffCache[ idx ] ) )
                    for idx in
                        ur_ui_idx_in_Idx ]

            θh =  [ angle( x_from_xr_xi(
                u_stateDiffCache[ idx ] ) )
                    for idx in
                        ur_ui_idx_in_Idx ]            
            
        end

        
        intra_dyn_ur_ui = [u_stateDiffCache...;]

        intra_vh_θh = [ vh; θh ]
        
    end

    #--------------------------------------------    
    #--------------------------------------------
    
    pf_sol_integ =
        dyn_pf_func!(
            intra_dyn_ur_ui, intra_dyn_pf_flat_para;
            
            ur_ui_DiffCache  =
                ur_ui_DiffCache,
            
            Pg_Qg_external_control =
                Pg_Qg_external_control,
            
            intra_pf_kwd_para =
                intra_pf_kwd_para,
            
            pf_solver =
                pf_solver,
            
            use_nlsolve =
                use_nlsolve )

    
    # pf_sol_integ =
    #     dyn_pf_by_vh_θh_func!(
    #         intra_vh_θh, intra_dyn_pf_flat_para;
    
    #         vh_θh_DiffCache =
    #             vh_θh_DiffCache,
    
    # Pg_Qg_external_control =
    #         Pg_Qg_external_control,

    #     intra_pf_kwd_para =
    #         intra_pf_kwd_para,

    #     pf_solver =
    #         pf_solver,

    #     use_nlsolve =
    #         use_nlsolve,

    #     post_pf_idxs =
    #         post_pf_idxs )


    
    """
    
    vh =  [ abs( x_from_xr_xi(
        nodes_ur_ui[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]

    θh =  [ angle( x_from_xr_xi(
        nodes_ur_ui[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]

    gens_vh_θh =
        [[a_vh, a_θh]
         for (a_vh, a_θh) in
             zip(gens_vh, gens_θh) ]
        
    non_gens_vh = vh[non_gens_nodes_idx]

    non_gens_θh = θh[non_gens_nodes_idx]

    non_gens_vh_θh =
        [[a_vh, a_θh]
         for (a_vh, a_θh) in
             zip(non_gens_vh, non_gens_θh) ]
    
    #----------------------------------------
    #----------------------------------------

    ur_ui_DiffCache =
        get_tmp(ur_ui_DiffCache, x)
        
    #----------------------------------------    

    intra_ur_ui_wt_pf_flat_para =
        [ nodes_ur_ui...; dyn_pf_flat_para... ]

    
    #----------------------------------------    
    
    pf_sol_integ =
        intra_dyn_pf_func!(
            intra_ur_ui_wt_pf_flat_para;
            
            ur_ui_DiffCache =
                ur_ui_DiffCache,
            
            Pg_Qg_external_control =
                Pg_Qg_external_control,
            
            intra_pf_kwd_para =
                intra_pf_kwd_para,
            
            pf_solver =
                pf_solver,
            
            use_nlsolve =
                use_nlsolve
        )

    """
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    vh_post_pf =
        pf_sol_integ[vh_Idx]
    
    θh_post_pf =
        pf_sol_integ[θh_Idx]
    
    gens_id_post_pf =
        pf_sol_integ[id_Idx]
    
    gens_iq_post_pf =
        pf_sol_integ[iq_Idx]
    
    #----------------------------------------

    nodes_vh_θh_post_pf =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(vh_post_pf ,
                 θh_post_pf)]
    
    #----------------------------------------

    vec_nodes_ur_ui_post_pf =
        polar_to_cartesian.( nodes_vh_θh_post_pf )
    
    #----------------------------------------

    gens_vh_post_pf =
        vh_post_pf[ gens_nodes_idx ]
 
    gens_θh_post_pf =
        θh_post_pf[ gens_nodes_idx ]

    gens_vh_θh_post_pf = [ [a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(gens_vh_post_pf ,
                 gens_θh_post_pf)]
    
    #----------------------------------------

    non_gens_vh_post_pf =
        vh_post_pf[ non_gens_nodes_idx ]

    non_gens_θh_post_pf =
        θh_post_pf[ non_gens_nodes_idx ]

    non_gens_vh_θh_post_pf = [ [a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(non_gens_vh_post_pf ,
                 non_gens_θh_post_pf)]

    #----------------------------------------

    gens_id_iq_post_pf =
        [ [a_gen_id, a_gen_iq]
         for (a_gen_id, a_gen_iq) in
             zip(gens_id_post_pf,
                 gens_iq_post_pf)]

    #----------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_nodes_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_nodes_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh_post_pf =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
         for (a_vh,
              a_θh,
              a_δ,
              a_ed_dash,
              a_eq_dash,
              a_ra,
              a_X_d_dash,
              a_X_q_dash ) in
                zip( gens_vh_post_pf,
                     gens_θh_post_pf,
                     gens_δ,
                     gens_ed_dash,
                     gens_eq_dash,
                     gens_ra,
                     gens_Xd_dash,
                     gens_Xq_dash ) ]
    
    #----------------------------------------
    #----------------------------------------

    per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
        [ [a_gen_ωs_ωref0_vref0_porder0;
          a_gen_id_iq_pg_vh_post_pf]
         for (a_gen_ωs_ωref0_vref0_porder0,
              a_gen_id_iq_pg_vh_post_pf) in
             zip( gens_nodes_ωs_ωref0_vref0_porder0,
                  id_iq_pg_vh_post_pf) ] 

    #----------------------------------------    
    # vtf  nodes
    #----------------------------------------

    """
    
    for ( a_vtf_du, a_vtf_u, a_vh_θh ) in
        zip(
            vtf_dx_u_views,
            vtf_x_u_views,
            nodes_vh_θh_post_pf )

         a_voltage_terminal_by_vh_θh_func!(
              a_vtf_du, a_vtf_u, a_vh_θh, t )
    end


    """
    
    for ( a_vtf_du, a_vtf_u, a_ur_ui ) in
        zip(
            vtf_dx_u_views,
            vtf_x_u_views,
            vec_nodes_ur_ui_post_pf )

         a_voltage_terminal_by_ur_ui_func!(
              a_vtf_du, a_vtf_u, a_ur_ui, t )
    end

    return nothing
    
end


function sd_dynamics_per_gen_ode_vtf_pf_by_vh_θh_func!(
    dx,
    x,
    flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,
    t;
    kwd_para = kwd_para )

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
     
     ode_vtf_by_vh_θh_kwd_para ) =
        kwd_para
    
    #----------------------------------------    

   (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
     
    (;
     pf_solver,
     use_nlsolve,
     init_dyn_pf_mismatch_kwd_para)  =
         init_pf_by_vh_θh_by_parts_kwd

    (;
     loc_load_exist,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,     
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash
     ) =
         init_dyn_pf_mismatch_kwd_para

    #----------------------------------------    
    #----------------------------------------    

     (;
      ode_vtf_by_vh_θh_para_kwd,
      ode_vtf_by_vh_θh_Idxs_kwd ) =
          ode_vtf_by_vh_θh_kwd_para 
    
    #----------------------------------------    
     
    (;
     gens_nodes_ra_Xd_dash_Xq_dash,
     gens_Ax_update_parameters,
     ode_per_gen_models_func_kwd_paras,
     vtf_gens_fun_kwd_para
     ) =
         ode_vtf_by_vh_θh_para_kwd
    
    (;
     flat_vh_θh_flat_ωs_ωref0_vref0_porder0_Idx,

     flat_vh_idx_in_Idx,
     flat_θh_idx_in_Idx,

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
     flat_δ_ed_dash_eq_dash_Idxs_in_state ) =
         ode_vtf_by_vh_θh_Idxs_kwd

    #----------------------------------------    
    
    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------    
    #----------------------------------------    

    flat_vh_flat_θh_idx, flat_ωs_ωref0_vref0_porder0_idx, flat_init_dyn_pf_para_idx = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx
    
    #----------------------------------------    

    flat_vh_flat_θh = flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_vh_flat_θh_idx ]
    
    flat_ωs_ωref0_vref0_porder0 =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[
            flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[
            flat_init_dyn_pf_para_idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ flat_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    

    flat_δ_ed_dash_eq_dash =
        x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    gens_nodes_δ_ed_dash_eq_dash =
        [ flat_δ_ed_dash_eq_dash[idx]
          for idx in
              δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    #----------------------------------------        

    """
    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    #----------------------------------------

    
    list_vtf_dx_nodes_u_views = [
        view( dx, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    
    list_vtf_x_nodes_u_views = [
        view( x, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    
    #----------------------------------------


    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
    
    #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           non_gens_nodes_u_Idx_in_ranges ]

    """
    #----------------------------------------        
    #----------------------------------------    
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state)
    
    im_vars_in_state = x[im_vars_Idx_in_state]
        
    #----------------------------------------    
    #----------------------------------------    


    """

    flat_δ_ed_dash_eq_dash_DiffCache =
        get_tmp(flat_δ_ed_dash_eq_dash_DiffCache,
                x[flat_δ_ed_dash_eq_dash_Idxs_in_state])
    
    flat_δ_ed_dash_eq_dash_DiffCache .=
        flat_δ_ed_dash_eq_dash_x0
    
    pf_sol =
        intra_dyn_pf_by_vh_θh_by_parts_func!(
            flat_vh_flat_θh,
            flat_δ_ed_dash_eq_dash_DiffCache,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_dyn_pf_kwd_para )

    """

    pf_sol = init_dyn_pf_by_vh_θh_by_parts_func!(
            flat_vh_flat_θh,
            init_dyn_pf_flat_para;
            init_pf_by_vh_θh_by_parts_kwd =
                init_pf_by_vh_θh_by_parts_kwd )
    
    
    if use_nlsolve == true
        
        intra_pf_sol = pf_sol.zero
        
    else
        
        intra_pf_sol = pf_sol.u
        
    end
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    flat_vh_post_pf =
        intra_pf_sol[vh_Idx]
    
    flat_θh_post_pf =
        intra_pf_sol[θh_Idx]

    flat_vh_flat_θh_post_pf =
        [ flat_vh_post_pf; flat_θh_post_pf]

    #----------------------------------------

    vec_vh_θh_post_pf =
        [ [ a_vh, a_θh ]
          for (a_vh, a_θh) in
              zip( flat_vh_post_pf, flat_θh_post_pf ) ]

    #----------------------------------------

    flat_vh_flat_θh_flat_gens_ωs_ωref0_vref0_porder0 =
        [flat_vh_flat_θh_post_pf;
         flat_ωs_ωref0_vref0_porder0 ]
    
    #----------------------------------------
    #----------------------------------------

    ode_vtf_by_vh_θh_func!(
        dx, x,
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0,
        t;
        kwd_para  =
            ode_vtf_by_vh_θh_kwd_para )
    
    
    return nothing


    # nodes_vh_θh_view =
    #     view( x, nodes_vh_θh_Idx_in_state )

    # nodes_vh_θh =
    #     x[ nodes_vh_θh_Idx_in_state  ]

    
    # gens_id_post_pf =
    #     intra_pf_sol[id_Idx]
    
    # gens_iq_post_pf =
    #     intra_pf_sol[iq_Idx]

    # (;
    #  counter_array,
     
    #  sim_state_x0,
    #  flat_vh_flat_θh_x0,
     
    #  stateDiffCache,
    #  flat_vh_flat_θh_DiffCache,     
    #  Ax_sparse_nzvalues_DiffCache,
     
    #  flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,
    #  init_pf_by_vh_θh_by_parts_kwd,
     
    #  ode_vtf_by_vh_θh_kwd_para,
     
    #  dynamics_by_per_gen_kwd_para ) =
    #     kwd_para

    
   # (;
   #  flat_net_vh_θh_flat_ωs_ωref0_vref0_porder0_Idx,
   #  gens_Ax_update_parameters,
   #  dynamics_by_per_gen_kwd_para) =
   #      ode_vtf_by_vh_θh_kwd_para 

    # (;
    #  flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state,

    #  nodes_vh_Idx_in_state,
    #  nodes_θh_Idx_in_state,
    #  nodes_vh_θh_Idx_in_state,

    #  vh_θh_idx_order,
    #  flat_vh_flat_θh_Idx,

    #  Pg_Qg_external_control,

    #  intra_pf_kwd_para,
    #  pf_solver,
    #  use_nlsolve,
    #  post_pf_idxs,
    #   Pg_Qg_Png_Qng_Pgll_Qgll_Idxs) =
    #       intra_dyn_pf_kwd_para
    
    #----------------------------------------    
    #----------------------------------------    
    
    # (;

    #  loc_load_exist,
     
    #  gens_nodes_idx,

    #  non_gens_nodes_idx,

    #  each_gens_im_vars_Idx_in_state,
     
    #  nodes_state_Idx,

    #  im_vars_Idx_in_state,

    #  nodes_ur_ui_Idx_in_state,

    #  nodes_u_Idx_in_ranges,

    #  gens_nodes_ra_Xd_dash_Xq_dash,

    #  nodes_δ_ed_dash_eq_dash_Idxs,
     
    #  post_pf_idxs,

    #  disaggretation_idxs,

    #  intra_pf_kwd_para,

    #  vtf_para_and_idxs,

    #  gens_nodes_collection,

    #  vec_Ax_views,

    #  ode_per_gen_models_func_kwd_paras     

    #  ) =
    #      dynamics_by_per_gen_kwd_para 

    # nodes_vh_θh_Idx_in_state =
    #     nodes_ur_ui_Idx_in_state
    
    # (;
    #  vh_Idx,
    #  θh_Idx,
    #  id_Idx,
    #  iq_Idx ) =
    #      post_pf_idxs
    
    # #----------------------------------------    
    
    # (;
    #  gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     
    #  a_gen_vtf_vh_θh_Idx,
     
    #  a_gen_vtf_δ_ed_dash_eq_dash_Idx,
     
          
    #  per_gen_ωs_ωref0_vref0_porder0_Idx,
     
    #  per_gen_id_iq_pg_vh_Idx,
     
     
    #  per_para_ωs_ωref0_vref0_porder0_Idx,
     
    #  per_para_id_iq_pg_vh_Idx,
     
     
    #  f_ωs_ωref0_vref0_porder0_Idx,
     
    #  f_dyn_pf_para_Idx,
    #  f_dyn_ode_pf_para_Idx,
     
     
    #  Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
     
    #  ) =
    #      disaggretation_idxs

    # #----------------------------------------    

    # (;          
    #  intra_flat_ur_ui_Idx,

    #  intra_flat_vh_θh_Idx,
     
    #  intra_dyn_pf_flat_para_Idx,

    #  intra_dyn_pf_mismatch_kwd_para,

    #  ur_ui_idx_in_Idx,
          
    #  nodes_u_Idx_in_ranges,
     
    #  nodes_δ_ed_dash_eq_dash_Idxs
     
    #  ) = intra_pf_kwd_para 

    # #----------------------------------------       
    
    # vh_θh_idx_in_Idx =
    #     ur_ui_idx_in_Idx
    
    
    # #----------------------------------------       
    
    # (;
    #  vtf_gens_fun_kwd_para,
     
    #  gens_nodes_u_Idx_in_ranges,
     
    #  non_gens_nodes_u_Idx_in_ranges,
    #  ) =
    #      vtf_para_and_idxs

    
    # #----------------------------------------    
    
    # (Pg_Idx,
    #  Qg_Idxs,
    #  Png_Idxs,
    #  Qng_Idxs,
    #  Pgll_Idxs,
    #  Qgll_Idxs) =
    #     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
        
    
end



function sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_flat_θh_func!(
    dx,
    x,
    flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,
    t;
    kwd_para = kwd_para )

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
     
     ode_vtf_by_vh_θh_kwd_para ) =
        kwd_para
    
    #----------------------------------------    

   (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
     
    (;
     pf_solver,
     use_nlsolve,
     init_dyn_pf_mismatch_kwd_para)  =
         init_pf_by_vh_θh_by_parts_kwd

    (;
     loc_load_exist,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,     
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     non_pre_ordered_dyn_pf_vars_Idxs
     ) =
         init_dyn_pf_mismatch_kwd_para

    #----------------------------------------    
    #----------------------------------------    

     (;
      ode_vtf_by_vh_θh_para_kwd,
      ode_vtf_by_vh_θh_Idxs_kwd ) =
          ode_vtf_by_vh_θh_kwd_para 
    
    #----------------------------------------    
     
    (;
     gens_nodes_ra_Xd_dash_Xq_dash,
     gens_Ax_update_parameters,
     ode_per_gen_models_func_kwd_paras,
     vtf_gens_fun_kwd_para
     ) =
         ode_vtf_by_vh_θh_para_kwd
    
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
     flat_δ_ed_dash_eq_dash_Idxs_in_state ) =
         ode_vtf_by_vh_θh_Idxs_kwd

    #----------------------------------------    
    
    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------    
    #----------------------------------------    

    flat_vh_flat_θh_idx, flat_ωs_ωref0_vref0_porder0_idx, flat_init_dyn_pf_para_idx = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx
    
    #----------------------------------------    

    flat_vh_flat_θh = flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_vh_flat_θh_idx ]
    
    flat_ωs_ωref0_vref0_porder0 =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[
            flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para =
        flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[
            flat_init_dyn_pf_para_idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ flat_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    

    flat_δ_ed_dash_eq_dash =
        x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    gens_nodes_δ_ed_dash_eq_dash =
        [ flat_δ_ed_dash_eq_dash[idx]
          for idx in
              δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    #----------------------------------------        

    """
    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    #----------------------------------------

    
    list_vtf_dx_nodes_u_views = [
        view( dx, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    
    list_vtf_x_nodes_u_views = [
        view( x, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    
    #----------------------------------------


    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
    
    #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           non_gens_nodes_u_Idx_in_ranges ]

    """
    #----------------------------------------        
    #----------------------------------------    
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state)
    
    im_vars_in_state = x[im_vars_Idx_in_state]
        
    #----------------------------------------    
    #----------------------------------------    


    """

    flat_δ_ed_dash_eq_dash_DiffCache =
        get_tmp(flat_δ_ed_dash_eq_dash_DiffCache,
                x[flat_δ_ed_dash_eq_dash_Idxs_in_state])
    
    flat_δ_ed_dash_eq_dash_DiffCache .=
        flat_δ_ed_dash_eq_dash_x0
    
    pf_sol =
        intra_dyn_pf_by_vh_θh_by_parts_func!(
            flat_vh_flat_θh,
            flat_δ_ed_dash_eq_dash_DiffCache,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_dyn_pf_kwd_para )

    """

    pf_sol = init_dyn_pf_by_vh_θh_by_parts_func!(
            flat_vh_flat_θh,
            init_dyn_pf_flat_para;
            init_pf_by_vh_θh_by_parts_kwd =
                init_pf_by_vh_θh_by_parts_kwd )
    
    
    if use_nlsolve == true
        
        intra_pf_sol = pf_sol.zero
        
    else
        
        intra_pf_sol = pf_sol.u
        
    end
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    flat_vh_post_pf =
        intra_pf_sol[vh_Idx]
    
    flat_θh_post_pf =
        intra_pf_sol[θh_Idx]

    flat_vh_flat_θh_post_pf =
        [ flat_vh_post_pf; flat_θh_post_pf]

    #----------------------------------------

    vec_vh_θh_post_pf =
        [ [ a_vh, a_θh ]
          for (a_vh, a_θh) in
              zip( flat_vh_post_pf, flat_θh_post_pf ) ]

    #----------------------------------------

    flat_vh_flat_θh_flat_gens_ωs_ωref0_vref0_porder0 =
        [flat_vh_flat_θh_post_pf;
         flat_ωs_ωref0_vref0_porder0 ]
    
    #----------------------------------------
    #----------------------------------------

    ode_vtf_by_vh_θh_func!(
        dx, x,
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0,
        t;
        kwd_para  =
            ode_vtf_by_vh_θh_kwd_para )
    
    
    return nothing
    
end


function sd_dynamics_per_gen_ode_vtf_pf_by_flat_vh_θh_func!(
    dx,
    x,
    flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para,
    t;
    kwd_para = kwd_para )

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
     
     ode_vtf_by_vh_θh_kwd_para ) =
        kwd_para
    
    #----------------------------------------    

   (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
     
    (;
     pf_solver,
     use_nlsolve,
     init_dyn_pf_mismatch_kwd_para)  =
         init_pf_by_vh_θh_by_parts_kwd

    (;
     loc_load_exist,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,     
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     non_pre_ordered_dyn_pf_vars_Idxs
     ) =
         init_dyn_pf_mismatch_kwd_para

    #----------------------------------------    
    #----------------------------------------    

     (;
      ode_vtf_by_vh_θh_para_kwd,
      ode_vtf_by_vh_θh_Idxs_kwd ) =
          ode_vtf_by_vh_θh_kwd_para 
    
    #----------------------------------------    
     
    (;
     gens_nodes_ra_Xd_dash_Xq_dash,
     gens_Ax_update_parameters,
     ode_per_gen_models_func_kwd_paras,
     vtf_gens_fun_kwd_para
     ) =
         ode_vtf_by_vh_θh_para_kwd
    
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
     flat_δ_ed_dash_eq_dash_Idxs_in_state ) =
         ode_vtf_by_vh_θh_Idxs_kwd

    #----------------------------------------        
    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------    
    #----------------------------------------    

    flat_vh_θh_idx, flat_ωs_ωref0_vref0_porder0_idx, flat_init_dyn_pf_para_idx = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx
    
    #----------------------------------------    

    flat_vh_θh = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_vh_θh_idx ]
    
    flat_ωs_ωref0_vref0_porder0 = flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para = flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_init_dyn_pf_para_idx ]

    #----------------------------------------    

    # flat_vh_idx_in_flat_Idx, flat_θh_idx_in_flat_Idx =
    #     flat_vh_flat_θh_Idx
    
    flat_vh = flat_vh_θh[ flat_vh_idx_in_Idx ]

    flat_θh = flat_vh_θh[ flat_θh_idx_in_Idx ]
    
    
    flat_vh_flat_θh = [flat_vh; flat_θh ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ flat_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    

    flat_δ_ed_dash_eq_dash =
        x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    gens_nodes_δ_ed_dash_eq_dash =
        [ flat_δ_ed_dash_eq_dash[idx]
          for idx in
              δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    #----------------------------------------        
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state)
    
    im_vars_in_state = x[im_vars_Idx_in_state]
        
    #----------------------------------------    
    #----------------------------------------    


    """

    flat_δ_ed_dash_eq_dash_DiffCache =
        get_tmp(flat_δ_ed_dash_eq_dash_DiffCache,
                x[flat_δ_ed_dash_eq_dash_Idxs_in_state])
    
    flat_δ_ed_dash_eq_dash_DiffCache .=
        flat_δ_ed_dash_eq_dash_x0
    
    pf_sol =
        intra_dyn_pf_by_vh_θh_by_parts_func!(
            flat_vh_flat_θh,
            flat_δ_ed_dash_eq_dash_DiffCache,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_dyn_pf_kwd_para )

    """

    pf_sol = init_dyn_pf_by_vh_θh_by_parts_func!(
            flat_vh_flat_θh,
            init_dyn_pf_flat_para;
            init_pf_by_vh_θh_by_parts_kwd =
                init_pf_by_vh_θh_by_parts_kwd )
        
    if use_nlsolve == true
        
        intra_pf_sol = pf_sol.zero
        
    else
        
        intra_pf_sol = pf_sol.u
        
    end
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    flat_vh_post_pf =
        intra_pf_sol[vh_Idx]
    
    flat_θh_post_pf =
        intra_pf_sol[θh_Idx]

    flat_vh_flat_θh_post_pf =
        [ flat_vh_post_pf; flat_θh_post_pf]

    #----------------------------------------

    vec_vh_θh_post_pf =
        [ [ a_vh, a_θh ]
          for (a_vh, a_θh) in
              zip( flat_vh_post_pf, flat_θh_post_pf ) ]

    #----------------------------------------

    flat_vh_flat_θh_flat_gens_ωs_ωref0_vref0_porder0 =
        [ flat_vh_flat_θh_post_pf;
         flat_ωs_ωref0_vref0_porder0 ]
    
    #----------------------------------------
    #----------------------------------------

    ode_vtf_by_vh_θh_func!(
        dx, x,
        flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0,
        t;
        kwd_para  =
            ode_vtf_by_vh_θh_kwd_para )
        
    return nothing
    
end


function ode_only_per_gen_model_by_id_iq_pg_vh_func!(
    dx, x,
    flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh,
    t;
    kwd_para  =
        ode_only_per_gen_id_iq_pg_vh_kwd )

    #----------------------------------------    

    (;
     flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,
     gens_id_iq_pg_vh_idx_in_Idx,
     ode_vtf_by_vh_θh_kwd_para ) =
         kwd_para
    
    #----------------------------------------    
    
    (;
     ode_vtf_by_vh_θh_para_kwd,
     ode_vtf_by_vh_θh_Idxs_kwd ) =
         ode_vtf_by_vh_θh_kwd_para
        
    (;
     gens_nodes_ra_Xd_dash_Xq_dash,
     gens_Ax_update_parameters,
     ode_per_gen_models_func_kwd_paras,
     vtf_gens_fun_kwd_para
      ) = ode_vtf_by_vh_θh_para_kwd
    
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
     flat_δ_ed_dash_eq_dash_Idxs_in_state) =
         ode_vtf_by_vh_θh_Idxs_kwd
    
    #----------------------------------------    

    flat_ωs_ωref0_vref0_porder0_Idx,flat_id_iq_pg_vh_Idx =
        flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx

        
     flat_ωs_ωref0_vref0_porder0 =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh[
            flat_ωs_ωref0_vref0_porder0_Idx ]
    
    flat_id_iq_pg_vh =
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh[
            flat_id_iq_pg_vh_Idx ]

    #----------------------------------------    
    
    gens_ωs_ωref0_vref0_porder0 =
        [ flat_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]

    
    gens_id_iq_pg_vh =
        [ flat_id_iq_pg_vh[idx]
          for idx in
              gens_id_iq_pg_vh_idx_in_Idx ]
        
    #----------------------------------------    
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              each_gens_im_vars_Idx_in_state ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              each_gens_im_vars_Idx_in_state ]
    
    #----------------------------------------    
    
    for ( dx_gen,
          x_gen,
          a_gen_im_vars_Idx_in_state,
          a_id_iq_pg_vh,
          a_ωs_ωref0_vref0_porder0,
          ode_per_gen_model_func_kwd_para,
          a_Ax_update_para ) in

        zip(list_dx_gens_states_views,
            list_x_gens_states_views,
            each_gens_im_vars_Idx_in_state,
            gens_id_iq_pg_vh,
            gens_ωs_ωref0_vref0_porder0,             
            ode_per_gen_models_func_kwd_paras,
            gens_Ax_update_parameters )

        (; Ae, Be,
         Ke, Te,
         vf_tilade_idx_in_Idx ) =
            a_Ax_update_para

        vf_tilade =
            x[a_gen_im_vars_Idx_in_state][
                vf_tilade_idx_in_Idx]
        
        # a_Ax_view, a_Bx_view, a_Cx_view, _, _ =
        #     ode_per_gen_model_func_kwd_para
        
        a_Ax_view, _, _, _, _ =
            ode_per_gen_model_func_kwd_para
        
        updated_Ax = get_a_im_updated_Ax(
            a_Ax_view, vf_tilade, a_Ax_update_para )

        # a_Ax_view[:,:] .= updated_Ax

        ode_per_gen_model_by_vh_θh_by_part_func!(
            dx_gen, x_gen,
            (a_ωs_ωref0_vref0_porder0,
             a_id_iq_pg_vh,
             updated_Ax
             # a_Ax_view
             ),
            t;
            ode_per_gen_model_func_kwd_para =
                ode_per_gen_model_func_kwd_para )
    
    end
    
end



function sd_dynamics_per_gen_ode_sep_vtf_pf_by_flat_vh_θh_func!( dx, x, flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para, t; kwd_para = sd_per_gen_ode_sep_vtf_pf_by_vh_θh_kwd_para )

    # flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para

    # flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para
    #----------------------------------------    
    
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
     ) =
        kwd_para
    #----------------------------------------    

   (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
     
    (;
     pf_solver,
     use_nlsolve,
     init_dyn_pf_mismatch_kwd_para)  =
         init_pf_by_vh_θh_by_parts_kwd

    (;
     loc_load_exist,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,     
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash,
     non_pre_ordered_dyn_pf_vars_Idxs
     ) =
         init_dyn_pf_mismatch_kwd_para

    #----------------------------------------    

    δ_ed_dash_eq_dash_Idxs_in_flattend = intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd.intra_dyn_pf_kwd_para.δ_ed_dash_eq_dash_Idxs_in_flattend

    #----------------------------------------    

   (;
    flat_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,
    gens_id_iq_pg_vh_idx_in_Idx,
    ode_vtf_by_vh_θh_kwd_para ) =
        ode_only_per_gen_id_iq_pg_vh_kwd

    (;
     ode_vtf_by_vh_θh_para_kwd,
     ode_vtf_by_vh_θh_Idxs_kwd ) =
         ode_vtf_by_vh_θh_kwd_para
        
    (;
     gens_nodes_ra_Xd_dash_Xq_dash,
     gens_Ax_update_parameters,
     ode_per_gen_models_func_kwd_paras,
     vtf_gens_fun_kwd_para
      ) = ode_vtf_by_vh_θh_para_kwd
    
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
     flat_δ_ed_dash_eq_dash_Idxs_in_state) =
         ode_vtf_by_vh_θh_Idxs_kwd

    #----------------------------------------
    # vtf
    #----------------------------------------

    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
           
    #----------------------------------------    

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           non_gens_nodes_u_Idx_in_ranges ]
            
    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    
    #----------------------------------------    

    flat_vh_θh_idx, flat_ωs_ωref0_vref0_porder0_idx, flat_init_dyn_pf_para_idx = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx
    
    #----------------------------------------    

    """
    flat_vh_θh = flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_vh_θh_idx ]
    
    flat_ωs_ωref0_vref0_porder0 = flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para = flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_init_dyn_pf_para_idx ]

    """

    flat_vh_flat_θh = flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_vh_θh_idx ]
    
    flat_ωs_ωref0_vref0_porder0 = flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para = flat_vh_flat_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[ flat_init_dyn_pf_para_idx ]
            
    #----------------------------------------    
    
    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ flat_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------        
    #----------------------------------------    
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state)
    
    im_vars_in_state = x[im_vars_Idx_in_state]
        
    #----------------------------------------    
    
    flat_δ_ed_dash_eq_dash =
        x[ flat_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    gens_nodes_δ_ed_dash_eq_dash =
        [ flat_δ_ed_dash_eq_dash[idx]
          for idx in
              δ_ed_dash_eq_dash_Idxs_in_flattend ]

    #----------------------------------------    

    flat_δ_ed_dash_eq_dash_DiffCache =
        get_tmp(flat_δ_ed_dash_eq_dash_DiffCache,
                flat_δ_ed_dash_eq_dash )
    
    flat_δ_ed_dash_eq_dash_DiffCache .=
        flat_δ_ed_dash_eq_dash
    
    #----------------------------------------    
    
    gens_δ =
        first.( gens_nodes_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_nodes_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_nodes_δ_ed_dash_eq_dash )    
    
    #----------------------------------------    
    #----------------------------------------    

    flat_vh_flat_θh_DiffCache =
        get_tmp(flat_vh_flat_θh_DiffCache,
                x[nodes_vh_θh_Idx_in_state] )

    flat_vh_flat_θh_DiffCache .=
        x[nodes_vh_θh_Idx_in_state]


    flat_vh_from_state =
        flat_vh_flat_θh_DiffCache[
            flat_vh_idx_in_Idx]
    
    flat_θh_from_state =
        flat_vh_flat_θh_DiffCache[
            flat_θh_idx_in_Idx]
   
    gens_vh_from_state =
        flat_vh_from_state[
            gens_nodes_idx ]

    gens_θh_from_state =
        flat_θh_from_state[
            gens_nodes_idx ]
    
    #----------------------------------------

    (flat_vh_idx_in_flat_vh_flat_θh,
     flat_θh_idx_in_flat_vh_flat_θh) =
         flat_vh_flat_θh_Idx

    
    flat_vh = flat_vh_flat_θh[
        flat_vh_idx_in_flat_vh_flat_θh ]

    flat_θh = flat_vh_flat_θh[
        flat_θh_idx_in_flat_vh_flat_θh ]

    
    gens_vh = flat_vh_flat_θh[ gens_nodes_idx ]

    gens_θh = flat_vh_flat_θh[ gens_nodes_idx ]        
    
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

    pf_sol =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_func!(
            init_flat_vh_flat_θh_cal_id_iq,
            flat_δ_ed_dash_eq_dash_DiffCache,
            init_dyn_pf_flat_para;
            kwd_para =
                intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd)
    
        
    if use_nlsolve == true
        
        intra_pf_sol = pf_sol.zero
        
    else
        
        intra_pf_sol = pf_sol.u
        
    end
    
    #----------------------------------------    
    #----------------------------------------    
    #----------------------------------------    
    
    """
    #----------------------------------------     

    gens_id_iq_from_state = [
        get_a_gen_dyn_idq(
            a_vh, a_θh, a_δ,
            a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
        for (a_vh, a_θh,
             a_δ, a_ed_dash, a_eq_dash,
             a_ra, a_X_d_dash, a_X_q_dash ) in
            zip( gens_vh_from_state, gens_θh_from_state,
                 gens_δ, gens_ed_dash, gens_eq_dash,
                 gens_ra, gens_Xd_dash, gens_Xq_dash ) ]

    gens_i_d_from_state =
        first.( gens_id_iq_from_state )
    
    gens_i_q_from_state =
        last.( gens_id_iq_from_state )
        
    init_flat_vh_flat_θh_cal_id_iq =
        [flat_vh_from_state,
         flat_θh_from_state,
         gens_i_d_from_state,
         gens_i_q_from_state ]


    pf_sol =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_func!(
            init_flat_vh_flat_θh_cal_id_iq,
            flat_δ_ed_dash_eq_dash_DiffCache,
            init_dyn_pf_flat_para;
            kwd_para =
                intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd)
            
    if use_nlsolve == true
        
        intra_pf_sol = pf_sol.zero
        
    else
        
        intra_pf_sol = pf_sol.u
        
    end

    """
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    flat_vh_post_pf =
        intra_pf_sol[vh_Idx]
    
    flat_θh_post_pf =
        intra_pf_sol[θh_Idx]
    
    gens_i_d_post_pf = intra_pf_sol[id_Idx]
    
    gens_i_q_post_pf = intra_pf_sol[iq_Idx]
    
    #----------------------------------------

    vec_vh_θh_post_pf =
        [ [ a_vh, a_θh ]
          for (a_vh, a_θh) in
              zip( flat_vh_post_pf, flat_θh_post_pf ) ]
    
    #----------------------------------------
    
    gens_vh_post_pf =
        flat_vh_post_pf[ gens_nodes_idx ]

    gens_θh_post_pf =
        flat_θh_post_pf[ gens_nodes_idx ]

    non_gens_vh_post_pf =
        flat_vh_post_pf[ non_gens_nodes_idx]

    non_gens_θh_post_pf =
        flat_θh_post_pf[ non_gens_nodes_idx]

    #----------------------------------------

    gens_vh_θh_post_pf =
        vec_vh_θh_post_pf[ gens_nodes_idx ]

    non_gens_vh_θh_post_pf =
        vec_vh_θh_post_pf[ non_gens_nodes_idx ]
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    vtf_gens_vh_θh_δ_ed_dash_eq_dash =
        [[a_gen_vh_θh; a_gen_δ_ed_dash_eq_dash]
         for (a_gen_vh_θh, a_gen_δ_ed_dash_eq_dash) in
             zip( gens_vh_θh_post_pf,
                  gens_nodes_δ_ed_dash_eq_dash ) ]
    
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ed_dash_eq_dash,
            t ;
            vtf_kwd_para =
                vtf_kwd_para )
        
    end
    
    #----------------------------------------    
    # vtf non gen nodes
    #----------------------------------------    

    for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh) in
        zip(
            vtf_dx_non_gens_u_views,
            vtf_x_non_gens_u_views,
            non_gens_vh_θh_post_pf )

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end

    #----------------------------------------
    #----------------------------------------
    
    """
    
    vtf_vh =
        flat_vh_flat_θh_DiffCache[flat_vh_idx_in_Idx]

    vtf_θh =
        flat_vh_flat_θh_DiffCache[flat_θh_idx_in_Idx]
    
    gens_vh_vtf = vtf_vh[ gens_nodes_idx ]

    gens_θh_vtf = vtf_θh[ gens_nodes_idx ]

    #----------------------------------------
    
    gens_id_iq_pg_vh_post_pf =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
         for (a_vh,  a_θh,
              a_δ, a_ed_dash, a_eq_dash,
              a_ra, a_X_d_dash, a_X_q_dash ) in
                zip( gens_vh_vtf, gens_θh_vtf,
                     gens_δ, gens_ed_dash, gens_eq_dash,
                     gens_ra, gens_Xd_dash, gens_Xq_dash )]

    #----------------------------------------
    #----------------------------------------

    gens_dyn_ph = get_gens_dynamic_ph_by_vhθh_δ_idq(
        gens_vh_post_pf,
        gens_θh_post_pf,
        gens_δ,
        gens_i_d_post_pf,
        gens_i_q_post_pf
    )
    
    gens_dyn_qh = get_gens_dynamic_qh_by_vhθh_δ_idq(
        gens_vh_post_pf,
        gens_θh_post_pf,
        gens_δ,
        gens_i_d_post_pf,
        gens_i_q_post_pf
    )
    
    gens_id_iq_pg_vh_post_pf =
        [[ a_id,  a_iq, a_pg, a_vh]
         for ( a_id,  a_iq, a_pg, a_vh ) in
                zip( gens_i_d_post_pf, gens_i_q_post_pf,
                     gens_dyn_ph, gens_vh_from_state) ]

    
    #----------------------------------------    
    
    gens_id_iq_pg_vh_post_pf =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh, a_θh,
            a_δ, a_ed_dash, a_eq_dash,
            a_ra, a_X_d_dash, a_X_q_dash )
         for (a_vh,  a_θh,
              a_δ, a_ed_dash, a_eq_dash,
              a_ra, a_X_d_dash, a_X_q_dash ) in
                zip( gens_vh_post_pf, gens_θh_post_pf,
                     gens_δ, gens_ed_dash, gens_eq_dash,
                     gens_ra, gens_Xd_dash, gens_Xq_dash )]

    """
    #----------------------------------------

    gens_dyn_ph = get_gens_dynamic_ph_by_vhθh_δ_idq(
        gens_vh_post_pf,
        gens_θh_post_pf,
        gens_δ,
        gens_i_d_post_pf,
        gens_i_q_post_pf
    )
        
    gens_id_iq_pg_vh_post_pf =
        [[ a_id,  a_iq, a_pg, a_vh]
         for ( a_id,  a_iq, a_pg, a_vh ) in
                zip( gens_i_d_post_pf, gens_i_q_post_pf,
                     gens_dyn_ph, gens_vh_from_state) ]

    #----------------------------------------
    
    per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
        [[a_gen_ωs_ωref0_vref0_porder0;
          a_gen_id_iq_pg_vh_post_pf]
         for (a_gen_ωs_ωref0_vref0_porder0,
              a_gen_id_iq_pg_vh_post_pf) in
             zip( gens_nodes_ωs_ωref0_vref0_porder0,
                  gens_id_iq_pg_vh_post_pf ) ] 

    #----------------------------------------
    # ode
    #----------------------------------------

    flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh =
        [ gens_nodes_ωs_ωref0_vref0_porder0...;
           gens_id_iq_pg_vh_post_pf... ]
    
    #----------------------------------------
    #----------------------------------------

    ode_only_per_gen_model_by_id_iq_pg_vh_func!(
        dx, x,
        flat_gens_ωs_ωref0_vref0_porder0_flat_id_iq_pg_vh,
        t;
        kwd_para  =
            ode_only_per_gen_id_iq_pg_vh_kwd )
    
    #----------------------------------------
        
    return nothing

    #----------------------------------------
    
end



function sd_dynamics_model_by_per_gen_func!(
    dx,
    x,
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
    t;
    kwd_para = kwd_para )

    (;Ax_sparse_nzvalues_DiffCache,
     intra_dyn_pf_kwd_para,     
     counter_array,     
     dynamics_by_per_gen_kwd_para ) =
        kwd_para

    #----------------------------------------    
    
    (;
     sim_state_x0,
     nodes_ur_ui_x0,

     stateDiffCache,
     ur_ui_DiffCache,

     flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state,

     nodes_ur_Idx_in_state,
     nodes_ui_Idx_in_state,

     non_consecutive_ur_ui_idxs,
          
     ur_ui_idx_order,
     flat_ur_flat_ui_Idx,
     Pg_Qg_external_control,
     intra_pf_kwd_para,
     pf_solver,
     use_nlsolve,
     post_pf_idxs,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs ) =
         intra_dyn_pf_kwd_para

    #----------------------------------------    
    
    (;

     loc_load_exist,
     
     gens_nodes_idx,

     non_gens_nodes_idx,

     each_gens_im_vars_Idx_in_state,
     
     nodes_state_Idx,

     im_vars_Idx_in_state,

     nodes_ur_ui_Idx_in_state,

     nodes_u_Idx_in_ranges,

     gens_nodes_ra_Xd_dash_Xq_dash,

     nodes_δ_ed_dash_eq_dash_Idxs,
     
     post_pf_idxs,

     disaggretation_idxs,

     intra_pf_kwd_para,

     vtf_para_and_idxs,

     gens_nodes_collection,

     vec_Ax_views,

     ode_per_gen_models_func_kwd_paras     

     ) =
         dynamics_by_per_gen_kwd_para 

     
    
    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------    
    
    (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (;
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     
     a_gen_vtf_vh_θh_Idx,
     
     a_gen_vtf_δ_ed_dash_eq_dash_Idx,
     
          
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     
     per_gen_id_iq_pg_vh_Idx,
     
     
     per_para_ωs_ωref0_vref0_porder0_Idx,
     
     per_para_id_iq_pg_vh_Idx,
     
     
     f_ωs_ωref0_vref0_porder0_Idx,
     
     f_dyn_pf_para_Idx,
     f_dyn_ode_pf_para_Idx,
     
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
     
     ) =
         disaggretation_idxs
    
    #----------------------------------------    

    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 

    #----------------------------------------       
    
    (;
     vtf_gens_fun_kwd_para,
     
     gens_nodes_u_Idx_in_ranges,
     
     non_gens_nodes_u_Idx_in_ranges,
     ) =
         vtf_para_and_idxs
    
    #----------------------------------------    

    f_gens_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_ωs_ωref0_vref0_porder0_Idx ]
    
    init_dyn_pf_flat_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_dyn_ode_pf_para_Idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ f_gens_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    

    
    # gens_nodes_δ_ed_dash_eq_dash =
    #     [ x[idx]
    #      for idx in
    #          nodes_δ_ed_dash_eq_dash_Idxs ]
    
    #  flat_δ_ed_dash_eq_dash =
    #     [gens_nodes_δ_ed_dash_eq_dash...;]


    flat_δ_ed_dash_eq_dash =
        x[ flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    gens_nodes_δ_ed_dash_eq_dash =
        [ flat_δ_ed_dash_eq_dash[idx]
          for idx in
              δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    
    (Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs) =
        Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
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

        intra_dyn_pf_flat_para = [ P_gens;
            Q_gens;            
            P_non_gens;
            Q_non_gens;
            flat_δ_ed_dash_eq_dash;
            P_g_loc_load;
            Q_g_loc_load] 
    else

        intra_dyn_pf_flat_para = [
            P_gens;
            Q_gens;
            P_non_gens;
            Q_non_gens;
            flat_δ_ed_dash_eq_dash ]   

    end
    
    #----------------------------------------        
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    #----------------------------------------

    list_vtf_dx_nodes_u_views = [
        view( dx, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    
    list_vtf_x_nodes_u_views = [
        view( x, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    
    #----------------------------------------

    """
    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
    
    #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           non_gens_nodes_u_Idx_in_ranges ]
    """
    
    #----------------------------------------        
    #----------------------------------------    
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state)
    
    im_vars_in_state = x[im_vars_Idx_in_state]
    
    #----------------------------------------

    nodes_ur_ui_view =
        view( x, nodes_ur_ui_Idx_in_state )

    nodes_ur_ui =
        x[ nodes_ur_ui_Idx_in_state ]

    stateDiffCache =
        get_tmp(stateDiffCache, x)
    
    ur_ui_DiffCache =
        get_tmp( ur_ui_DiffCache, nodes_ur_ui )

    counter = counter_array[1]
    
    if counter == 1
     
        stateDiffCache .= sim_state_x0

        ur_ui_DiffCache .= nodes_ur_ui_x0

        counter_array .+= 1
        
    else
        
        stateDiffCache .= x

        ur_ui_DiffCache .= nodes_ur_ui_view
        
    end
    
    #----------------------------------------    
    #----------------------------------------    

    """
    intra_pf_sol =
        intra_dyn_pf_by_ur_ui_by_parts_func!(
            ur_ui_DiffCache,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_dyn_pf_kwd_para )
    """

    intra_pf_sol =
        intra_dyn_pf_by_ur_ui_by_stateDiff_func!(
            stateDiffCache,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_dyn_pf_kwd_para )    
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    vh_post_pf =
        intra_pf_sol[vh_Idx]
    
    θh_post_pf =
        intra_pf_sol[θh_Idx]
    
    gens_id_post_pf =
        intra_pf_sol[id_Idx]
    
    gens_iq_post_pf =
        intra_pf_sol[iq_Idx]

    #----------------------------------------

    nodes_vh_θh_post_pf =
        [ [ a_vh, a_θh ]
          for (a_vh, a_θh) in
              zip( vh_post_pf, θh_post_pf ) ]

    #----------------------------------------

    nodes_ur_ui_post_pf =
        polar_to_cartesian.( nodes_vh_θh_post_pf )
    
    #----------------------------------------

    gens_vh_post_pf =
        vh_post_pf[ gens_nodes_idx ]

    gens_θh_post_pf =
        θh_post_pf[ gens_nodes_idx ]

    gens_vh_θh_post_pf =
        [ [ a_gen_vh, a_gen_θh ]
         for ( a_gen_vh, a_gen_θh ) in
             zip( gens_vh_post_pf ,
                  gens_θh_post_pf )]
    
    #----------------------------------------

    non_gens_vh_post_pf =
        vh_post_pf[ non_gens_nodes_idx ]

    non_gens_θh_post_pf =
        θh_post_pf[ non_gens_nodes_idx ]

    non_gens_vh_θh_post_pf =
        [[ a_gen_vh, a_gen_θh ]
         for ( a_gen_vh, a_gen_θh ) in
             zip( non_gens_vh_post_pf ,
                  non_gens_θh_post_pf )]

    #----------------------------------------

    gens_id_iq_post_pf =
        [ [ a_gen_id, a_gen_iq ]
         for ( a_gen_id, a_gen_iq ) in
             zip( gens_id_post_pf,
                  gens_iq_post_pf )]

    #----------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_nodes_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_nodes_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh_post_pf =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
         for (a_vh,
              a_θh,
              a_δ,
              a_ed_dash,
              a_eq_dash,
              a_ra,
              a_X_d_dash,
              a_X_q_dash ) in
                zip( gens_vh_post_pf,
                     gens_θh_post_pf,
                     gens_δ,
                     gens_ed_dash,
                     gens_eq_dash,
                     gens_ra,
                     gens_Xd_dash,
                     gens_Xq_dash ) ]
    
    #----------------------------------------
    #----------------------------------------
    
    # stateDiffCache_gen =
    #     [stateDiffCache[idx]
    #      for idx in
    #          each_gens_im_vars_Idx_in_state ]
    
    # per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     get_per_vars_or_paras_to_per_node(
    #         [ gens_nodes_ωs_ωref0_vref0_porder0,
    #          id_iq_pg_vh_post_pf ] )

    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    # update_gens_nodes_im_Ax_system_matrices!(
    #     vec_Ax_views,
    #     im_vars_in_state,
    #     each_gens_im_vars_Idx_in_state,
    #     gens_nodes_collection )


    # for (a_Ax, a_gen_state_Idx, a_gen_node ) in
    #     zip( vec_Ax_views,
    #          each_gens_im_vars_Idx_in_state,
    #          gens_nodes_collection )

    #     update_a_im_plant_Ax_system_matrices!(
    #         a_Ax, x[a_gen_state_Idx], a_gen_node )
        
    # end
        
    #----------------------------------------
    #----------------------------------------

    # per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     [[a_gen_ωs_ωref0_vref0_porder0;
    #       a_gen_id_iq_pg_vh_post_pf]
    #      for (a_gen_ωs_ωref0_vref0_porder0,
    #           a_gen_id_iq_pg_vh_post_pf) in
    #          zip( gens_nodes_ωs_ωref0_vref0_porder0,
    #               id_iq_pg_vh_post_pf) ] 
        
    # for ( dx_gen,
    #       x_gen,
    #       per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    #       ode_per_gen_model_func_kwd_para ) in

    #     zip( list_dx_gens_states_views,
    #          list_x_gens_states_views, 
    #          per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    #          ode_per_gen_models_func_kwd_paras)
        
    #     ode_per_gen_model_func!(
    #         dx_gen,
    #         x_gen,
    #         per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    #         t
    #         ;ode_per_gen_model_func_kwd_para =
    #             ode_per_gen_model_func_kwd_para )
        
    # end
        
    #----------------------------------------
    #----------------------------------------
    
    for ( a_gen_im_vars_Idx_in_state,
          a_id_iq_pg_vh,
          a_ωs_ωref0_vref0_porder0,
          ode_per_gen_model_func_kwd_para,
          a_gen_node ) in

        zip( each_gens_im_vars_Idx_in_state,
             id_iq_pg_vh_post_pf,
             gens_nodes_ωs_ωref0_vref0_porder0,             
             ode_per_gen_models_func_kwd_paras,
             gens_nodes_collection)
        
        a_Ax, a_Bx, a_Cx, _, _ =
            ode_per_gen_model_func_kwd_para

        # update_a_im_plant_Ax_system_matrices!(
        #     a_Ax,
        #     x[a_gen_im_vars_Idx_in_state],
        #     a_gen_node )

        (sparse_row_idxs,
         sparse_col_idxs,
         Ax_sparse_nzvalues ) =
             findnz( sparse( Matrix( a_Ax )) )
        
        Ax_sparse_nzvalues_DiffCache =
            get_tmp(Ax_sparse_nzvalues_DiffCache,
                    Ax_sparse_nzvalues )
        
        Ax_sparse_nzvalues_DiffCache .=
            Ax_sparse_nzvalues

        a_stateDiffCache =
            stateDiffCache[a_gen_im_vars_Idx_in_state]
        
        (; update_γ, γ_idx_in_sp_nzv ) =
            get_a_im_Ax_update_value_and_Idx(
                sparse_row_idxs, sparse_col_idxs,
                a_stateDiffCache, a_gen_node )

        Ax_sparse_nzvalues =
            [Ax_sparse_nzvalues[1:γ_idx_in_sp_nzv-1];
             [update_γ];
             Ax_sparse_nzvalues[γ_idx_in_sp_nzv+1:end] ]
        
        updated_Ax =
            sparse( sparse_row_idxs,
                    sparse_col_idxs,
                    Ax_sparse_nzvalues)
        

            dx[a_gen_im_vars_Idx_in_state] .=
            updated_Ax * a_stateDiffCache +
            a_Bx * a_id_iq_pg_vh +
            a_Cx * a_ωs_ωref0_vref0_porder0
        
    end

    #----------------------------------------
    # all nodes
    #----------------------------------------

    for ( a_vtf_dx, a_vtf_x, a_node_ur_ui_post_pf ) in
        zip( list_vtf_dx_nodes_u_views,
             list_vtf_x_nodes_u_views,
             nodes_ur_ui_post_pf )

        a_voltage_terminal_by_ur_ui_func!(
            a_vtf_dx, a_vtf_x, a_node_ur_ui_post_pf, t)

    end
            
    #----------------------------------------
    #----------------------------------------
    
    """
    vtf_gens_vh_θh_δ_ed_dash_eq_dash =
        get_per_vars_or_paras_to_per_node(
            [ gens_vh_θh_post_pf,
              gens_nodes_δ_ed_dash_eq_dash ] )
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        
        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ed_dash_eq_dash,
            t
            ; vtf_kwd_para =
                vtf_kwd_para )
    end

    #----------------------------------------    
    # vtf non gen nodes
    #----------------------------------------

    for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh) in
        zip(
            vtf_dx_non_gens_u_views,
            vtf_x_non_gens_u_views,
            non_gens_vh_θh_post_pf )

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end

    """

    
    return nothing
    
end


function sd_dynamics_model_by_per_gen_by_vh_θh_func!(
    dx,
    x,
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
    t;
    kwd_para = kwd_para )

    (;Ax_sparse_nzvalues_DiffCache,
     intra_dyn_pf_kwd_para,     
     counter_array,     
     dynamics_by_per_gen_kwd_para ) =
        kwd_para

    #----------------------------------------    
    
    (;
     sim_state_x0_by_θh_x0,
     nodes_vh_θh_x0,

     stateDiffCache,     
     vh_θh_DiffCache,

     flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state,

     nodes_vh_Idx_in_state,     
     nodes_θh_Idx_in_state,

     non_consecutive_vh_θh_idxs,
          
     vh_θh_idx_order,
     flat_vh_flat_θh_Idx,
     Pg_Qg_external_control,
     intra_pf_kwd_para,
     pf_solver,
     use_nlsolve,
     post_pf_idxs,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs ) =
         intra_dyn_pf_kwd_para

    #----------------------------------------    
    
    (;

     loc_load_exist,
     
     gens_nodes_idx,

     non_gens_nodes_idx,

     each_gens_im_vars_Idx_in_state,
     
     nodes_state_Idx,

     im_vars_Idx_in_state,

     nodes_ur_ui_Idx_in_state,

     nodes_u_Idx_in_ranges,

     gens_nodes_ra_Xd_dash_Xq_dash,

     nodes_δ_ed_dash_eq_dash_Idxs,
     
     post_pf_idxs,

     disaggretation_idxs,

     intra_pf_kwd_para,

     vtf_para_and_idxs,

     gens_nodes_collection,

     vec_Ax_views,

     ode_per_gen_models_func_kwd_paras     

     ) =
         dynamics_by_per_gen_kwd_para 

    nodes_vh_θh_Idx_in_state =
        nodes_ur_ui_Idx_in_state
    
    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------    
    
    (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (;
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     
     a_gen_vtf_vh_θh_Idx,
     
     a_gen_vtf_δ_ed_dash_eq_dash_Idx,
     
          
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     
     per_gen_id_iq_pg_vh_Idx,
     
     
     per_para_ωs_ωref0_vref0_porder0_Idx,
     
     per_para_id_iq_pg_vh_Idx,
     
     
     f_ωs_ωref0_vref0_porder0_Idx,
     
     f_dyn_pf_para_Idx,
     f_dyn_ode_pf_para_Idx,
     
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
     
     ) =
         disaggretation_idxs

    #----------------------------------------    

    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 

    #----------------------------------------       
    
    vh_θh_idx_in_Idx =
        ur_ui_idx_in_Idx
    
    
    #----------------------------------------       
    
    (;
     vtf_gens_fun_kwd_para,
     
     gens_nodes_u_Idx_in_ranges,
     
     non_gens_nodes_u_Idx_in_ranges,
     ) =
         vtf_para_and_idxs
    
    #----------------------------------------    

    f_gens_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_ωs_ωref0_vref0_porder0_Idx ]
    
    init_dyn_pf_flat_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_dyn_ode_pf_para_Idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ f_gens_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    

    
    # gens_nodes_δ_ed_dash_eq_dash =
    #     [ x[idx]
    #      for idx in
    #          nodes_δ_ed_dash_eq_dash_Idxs ]
    
    #  flat_δ_ed_dash_eq_dash =
    #     [gens_nodes_δ_ed_dash_eq_dash...;]


    flat_δ_ed_dash_eq_dash =
        x[ flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    gens_nodes_δ_ed_dash_eq_dash =
        [ flat_δ_ed_dash_eq_dash[idx]
          for idx in
              δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    
    (Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs) =
        Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
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

        intra_dyn_pf_flat_para = [ P_gens;
            Q_gens;            
            P_non_gens;
            Q_non_gens;
            flat_δ_ed_dash_eq_dash;
            P_g_loc_load;
            Q_g_loc_load] 
    else

        intra_dyn_pf_flat_para = [
            P_gens;
            Q_gens;
            P_non_gens;
            Q_non_gens;
            flat_δ_ed_dash_eq_dash ]   

    end
    
    #----------------------------------------        
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    #----------------------------------------

    """
    list_vtf_dx_nodes_u_views = [
        view( dx, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    
    list_vtf_x_nodes_u_views = [
        view( x, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    """
    #----------------------------------------


    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
    
    #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           non_gens_nodes_u_Idx_in_ranges ]

    
    #----------------------------------------        
    #----------------------------------------    
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state)
    
    im_vars_in_state = x[im_vars_Idx_in_state]
    
    #----------------------------------------

    nodes_vh_θh_view =
        view( x, nodes_vh_θh_Idx_in_state )

    nodes_vh_θh =
        x[ nodes_vh_θh_Idx_in_state  ]

    stateDiffCache =
        get_tmp(stateDiffCache, x)
    
    vh_θh_DiffCache =
        get_tmp( vh_θh_DiffCache, nodes_vh_θh )

    counter = counter_array[1]
    
    if counter == 1
     
        stateDiffCache .= sim_state_x0

        vh_θh_DiffCache .= nodes_vh_θh_x0

        counter_array .+= 1
        
    else
        
        stateDiffCache .= x

        vh_θh_DiffCache .= x[ nodes_vh_θh_Idx_in_state ]
        
    end
    
    #----------------------------------------    
    #----------------------------------------    

    """
    intra_pf_sol =
        intra_dyn_pf_by_vh_θh_by_parts_func!(
            vh_θh_DiffCache,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_dyn_pf_kwd_para )
    """

    intra_pf_sol =
        intra_dyn_pf_by_vh_θh_by_stateDiff_func!(
            stateDiffCache,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_dyn_pf_kwd_para )    

    if use_nlsolve == true
        
        intra_pf_sol = pf_sol.zero
        
    else
        
        intra_pf_sol = pf_sol.u
        
    end
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    vh_post_pf =
        intra_pf_sol[vh_Idx]
    
    θh_post_pf =
        intra_pf_sol[θh_Idx]
    
    gens_id_post_pf =
        intra_pf_sol[id_Idx]
    
    gens_iq_post_pf =
        intra_pf_sol[iq_Idx]

    #----------------------------------------

    nodes_vh_θh_post_pf =
        [ [ a_vh, a_θh ]
          for (a_vh, a_θh) in
              zip( vh_post_pf, θh_post_pf ) ]

    #----------------------------------------

    # nodes_ur_ui_post_pf =
    #     polar_to_cartesian.( nodes_vh_θh_post_pf )
    
    #----------------------------------------

    gens_vh_post_pf =
        vh_post_pf[ gens_nodes_idx ]

    gens_θh_post_pf =
        θh_post_pf[ gens_nodes_idx ]

    gens_vh_θh_post_pf =
        [ [ a_gen_vh, a_gen_θh ]
         for ( a_gen_vh, a_gen_θh ) in
             zip( gens_vh_post_pf ,
                  gens_θh_post_pf )]
    
    #----------------------------------------

    non_gens_vh_post_pf =
        vh_post_pf[ non_gens_nodes_idx ]

    non_gens_θh_post_pf =
        θh_post_pf[ non_gens_nodes_idx ]

    non_gens_vh_θh_post_pf =
        [[ a_gen_vh, a_gen_θh ]
         for ( a_gen_vh, a_gen_θh ) in
             zip( non_gens_vh_post_pf ,
                  non_gens_θh_post_pf )]

    #----------------------------------------

    gens_id_iq_post_pf =
        [ [ a_gen_id, a_gen_iq ]
         for ( a_gen_id, a_gen_iq ) in
             zip( gens_id_post_pf,
                  gens_iq_post_pf )]

    #----------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_nodes_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_nodes_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh_post_pf =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
         for (a_vh,
              a_θh,
              a_δ,
              a_ed_dash,
              a_eq_dash,
              a_ra,
              a_X_d_dash,
              a_X_q_dash ) in
                zip( gens_vh_post_pf,
                     gens_θh_post_pf,
                     gens_δ,
                     gens_ed_dash,
                     gens_eq_dash,
                     gens_ra,
                     gens_Xd_dash,
                     gens_Xq_dash ) ]
    
    #----------------------------------------
    #----------------------------------------
    
    # stateDiffCache_gen =
    #     [stateDiffCache[idx]
    #      for idx in
    #          each_gens_im_vars_Idx_in_state ]
    
    # per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     get_per_vars_or_paras_to_per_node(
    #         [ gens_nodes_ωs_ωref0_vref0_porder0,
    #          id_iq_pg_vh_post_pf ] )

    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    # update_gens_nodes_im_Ax_system_matrices!(
    #     vec_Ax_views,
    #     im_vars_in_state,
    #     each_gens_im_vars_Idx_in_state,
    #     gens_nodes_collection )


    # for (a_Ax, a_gen_state_Idx, a_gen_node ) in
    #     zip( vec_Ax_views,
    #          each_gens_im_vars_Idx_in_state,
    #          gens_nodes_collection )

    #     update_a_im_plant_Ax_system_matrices!(
    #         a_Ax, x[a_gen_state_Idx], a_gen_node )
        
    # end
        
    #----------------------------------------
    #----------------------------------------

    # per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     [[a_gen_ωs_ωref0_vref0_porder0;
    #       a_gen_id_iq_pg_vh_post_pf]
    #      for (a_gen_ωs_ωref0_vref0_porder0,
    #           a_gen_id_iq_pg_vh_post_pf) in
    #          zip( gens_nodes_ωs_ωref0_vref0_porder0,
    #               id_iq_pg_vh_post_pf) ] 
        
    # for ( dx_gen,
    #       x_gen,
    #       per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    #       ode_per_gen_model_func_kwd_para ) in

    #     zip( list_dx_gens_states_views,
    #          list_x_gens_states_views, 
    #          per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    #          ode_per_gen_models_func_kwd_paras)
        
    #     ode_per_gen_model_func!(
    #         dx_gen,
    #         x_gen,
    #         per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    #         t
    #         ;ode_per_gen_model_func_kwd_para =
    #             ode_per_gen_model_func_kwd_para )
        
    # end
        
    #----------------------------------------
    #----------------------------------------
    
    for ( a_gen_im_vars_Idx_in_state,
          a_id_iq_pg_vh,
          a_ωs_ωref0_vref0_porder0,
          ode_per_gen_model_func_kwd_para,
          a_gen_node ) in

        zip( each_gens_im_vars_Idx_in_state,
             id_iq_pg_vh_post_pf,
             gens_nodes_ωs_ωref0_vref0_porder0,             
             ode_per_gen_models_func_kwd_paras,
             gens_nodes_collection)
        
        a_Ax, a_Bx, a_Cx, _, _ =
            ode_per_gen_model_func_kwd_para

        # update_a_im_plant_Ax_system_matrices!(
        #     a_Ax,
        #     x[a_gen_im_vars_Idx_in_state],
        #     a_gen_node )

        (sparse_row_idxs,
         sparse_col_idxs,
         Ax_sparse_nzvalues ) =
             findnz( sparse( Matrix( a_Ax )) )
        
        Ax_sparse_nzvalues_DiffCache =
            get_tmp(Ax_sparse_nzvalues_DiffCache,
                    Ax_sparse_nzvalues )
        
        Ax_sparse_nzvalues_DiffCache .=
            Ax_sparse_nzvalues

        a_stateDiffCache =
            stateDiffCache[a_gen_im_vars_Idx_in_state]
        
        (; update_γ, γ_idx_in_sp_nzv ) =
            get_a_im_Ax_update_value_and_Idx(
                sparse_row_idxs, sparse_col_idxs,
                a_stateDiffCache, a_gen_node )

        Ax_sparse_nzvalues =
            [Ax_sparse_nzvalues[1:γ_idx_in_sp_nzv-1];
             [update_γ];
             Ax_sparse_nzvalues[γ_idx_in_sp_nzv+1:end] ]
        
        updated_Ax =
            sparse( sparse_row_idxs,
                    sparse_col_idxs,
                    Ax_sparse_nzvalues)
        

            dx[a_gen_im_vars_Idx_in_state] .=
            updated_Ax * a_stateDiffCache +
            a_Bx * a_id_iq_pg_vh +
            a_Cx * a_ωs_ωref0_vref0_porder0
        
    end

    #----------------------------------------
    # all nodes
    #----------------------------------------

    """
    for ( a_vtf_dx, a_vtf_x, a_node_ur_ui_post_pf ) in
        zip( list_vtf_dx_nodes_u_views,
             list_vtf_x_nodes_u_views,
             nodes_vh_θh_post_pf )

        a_voltage_terminal_by_vh_θh_func!(
            a_vtf_dx, a_vtf_x, a_node_vh_θh_post_pf, t)

    end

    """
    
    #----------------------------------------
    #----------------------------------------
    

    vtf_gens_vh_θh_δ_ed_dash_eq_dash =
        get_per_vars_or_paras_to_per_node(
            [ gens_vh_θh_post_pf,
              gens_nodes_δ_ed_dash_eq_dash ] )
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        
        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ed_dash_eq_dash,
            t
            ; vtf_kwd_para =
                vtf_kwd_para )
    end

    #----------------------------------------    
    # vtf non gen nodes
    #----------------------------------------

    for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh) in
        zip(
            vtf_dx_non_gens_u_views,
            vtf_x_non_gens_u_views,
            non_gens_vh_θh_post_pf )

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end
    
    return nothing
    
end



function sd_dynamics_model_by_per_gen_by_ur_ui_func!(
    dx,
    x,
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
    t;
    kwd_para = kwd_para )

    (;Ax_sparse_nzvalues_DiffCache,
     intra_dyn_pf_kwd_para,     
     counter_array,     
     dynamics_by_per_gen_kwd_para ) =
        kwd_para

    #----------------------------------------    
    
    (;
     sim_state_x0,
     nodes_ur_ui_x0,

     stateDiffCache,
     ur_ui_DiffCache,

     flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state,

     nodes_ur_Idx_in_state,
     nodes_ui_Idx_in_state,

     non_consecutive_ur_ui_idxs,
          
     ur_ui_idx_order,
     flat_ur_flat_ui_Idx,
     Pg_Qg_external_control,
     intra_pf_kwd_para,
     pf_solver,
     use_nlsolve,
     post_pf_idxs,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs ) =
         intra_dyn_pf_kwd_para

    #----------------------------------------    
    
    (;

     loc_load_exist,
     
     gens_nodes_idx,

     non_gens_nodes_idx,

     each_gens_im_vars_Idx_in_state,
     
     nodes_state_Idx,

     im_vars_Idx_in_state,

     nodes_ur_ui_Idx_in_state,

     nodes_u_Idx_in_ranges,

     gens_nodes_ra_Xd_dash_Xq_dash,

     nodes_δ_ed_dash_eq_dash_Idxs,
     
     post_pf_idxs,

     disaggretation_idxs,

     intra_pf_kwd_para,

     vtf_para_and_idxs,

     gens_nodes_collection,

     vec_Ax_views,

     ode_per_gen_models_func_kwd_paras     

     ) =
         dynamics_by_per_gen_kwd_para 

     
    
    #----------------------------------------
    # checked
    #----------------------------------------       

    (; pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------    
    
    (;
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (;
     gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
     
     a_gen_vtf_vh_θh_Idx,
     
     a_gen_vtf_δ_ed_dash_eq_dash_Idx,
     
          
     per_gen_ωs_ωref0_vref0_porder0_Idx,
     
     per_gen_id_iq_pg_vh_Idx,
     
     
     per_para_ωs_ωref0_vref0_porder0_Idx,
     
     per_para_id_iq_pg_vh_Idx,
     
     
     f_ωs_ωref0_vref0_porder0_Idx,
     
     f_dyn_pf_para_Idx,
     f_dyn_ode_pf_para_Idx,
     
     
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
     
     ) =
         disaggretation_idxs
    
    #----------------------------------------    

    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 

    #----------------------------------------       
    
    (;
     vtf_gens_fun_kwd_para,
     
     gens_nodes_u_Idx_in_ranges,
     
     non_gens_nodes_u_Idx_in_ranges,
     ) =
         vtf_para_and_idxs
    
    #----------------------------------------    

    f_gens_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_ωs_ωref0_vref0_porder0_Idx ]
    
    init_dyn_pf_flat_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para[
            f_dyn_ode_pf_para_Idx ]
    
    #----------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        [ f_gens_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_ωs_ωref0_vref0_porder0_idx_in_Idx ]
    
    #----------------------------------------    

    
    # gens_nodes_δ_ed_dash_eq_dash =
    #     [ x[idx]
    #      for idx in
    #          nodes_δ_ed_dash_eq_dash_Idxs ]
    
    #  flat_δ_ed_dash_eq_dash =
    #     [gens_nodes_δ_ed_dash_eq_dash...;]


    flat_δ_ed_dash_eq_dash =
        x[ flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    gens_nodes_δ_ed_dash_eq_dash =
        [ flat_δ_ed_dash_eq_dash[idx]
          for idx in
              δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------    
    
    (Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs) =
        Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
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

        intra_dyn_pf_flat_para = [ P_gens;
            Q_gens;            
            P_non_gens;
            Q_non_gens;
            flat_δ_ed_dash_eq_dash;
            P_g_loc_load;
            Q_g_loc_load] 
    else

        intra_dyn_pf_flat_para = [
            P_gens;
            Q_gens;
            P_non_gens;
            Q_non_gens;
            flat_δ_ed_dash_eq_dash ]   

    end
    
    #----------------------------------------        
    #----------------------------------------    

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    #----------------------------------------

    list_vtf_dx_nodes_u_views = [
        view( dx, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    
    list_vtf_x_nodes_u_views = [
        view( x, node_u_Idx)
        for node_u_Idx in
            nodes_u_Idx_in_ranges]
    
    #----------------------------------------

   
    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]
    
    #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           non_gens_nodes_u_Idx_in_ranges ]
    
    
    #----------------------------------------        
    #----------------------------------------    
    
    im_vars_view_in_state =
        view(x, im_vars_Idx_in_state)
    
    im_vars_in_state = x[im_vars_Idx_in_state]
    
    #----------------------------------------

    nodes_ur_ui_view =
        view( x, nodes_ur_ui_Idx_in_state )

    nodes_ur_ui =
        x[ nodes_ur_ui_Idx_in_state ]

    stateDiffCache =
        get_tmp(stateDiffCache, x)
    
    ur_ui_DiffCache =
        get_tmp( ur_ui_DiffCache, nodes_ur_ui )

    counter = counter_array[1]
    
    if counter == 1
     
        stateDiffCache .= sim_state_x0

        ur_ui_DiffCache .= nodes_ur_ui_x0

        counter_array .+= 1
        
    else
        
        stateDiffCache .= x

        ur_ui_DiffCache .= nodes_ur_ui_view

    end
    
    #----------------------------------------    
    #----------------------------------------    

    """
    intra_pf_sol =
        intra_dyn_pf_by_ur_ui_by_parts_func!(
            ur_ui_DiffCache,
            flat_δ_ed_dash_eq_dash,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_dyn_pf_kwd_para )
    """

    intra_pf_sol =
        intra_dyn_pf_by_ur_ui_by_stateDiff_func!(
            stateDiffCache,
            init_dyn_pf_flat_para;
            intra_dyn_pf_kwd_para =
                intra_dyn_pf_kwd_para )    
    
    #----------------------------------------
    # pf solution extractions and calc
    #----------------------------------------

    vh_post_pf =
        intra_pf_sol[vh_Idx]
    
    θh_post_pf =
        intra_pf_sol[θh_Idx]
    
    gens_id_post_pf =
        intra_pf_sol[id_Idx]
    
    gens_iq_post_pf =
        intra_pf_sol[iq_Idx]

    #----------------------------------------

    nodes_vh_θh_post_pf =
        [ [ a_vh, a_θh ]
          for (a_vh, a_θh) in
              zip( vh_post_pf, θh_post_pf ) ]

    #----------------------------------------

    nodes_ur_ui_post_pf =
        polar_to_cartesian.( nodes_vh_θh_post_pf )
    
    #----------------------------------------

    gens_vh_post_pf =
        vh_post_pf[ gens_nodes_idx ]

    gens_θh_post_pf =
        θh_post_pf[ gens_nodes_idx ]

    gens_vh_θh_post_pf =
        [ [ a_gen_vh, a_gen_θh ]
         for ( a_gen_vh, a_gen_θh ) in
             zip( gens_vh_post_pf ,
                  gens_θh_post_pf )]
    
    #----------------------------------------

    non_gens_vh_post_pf =
        vh_post_pf[ non_gens_nodes_idx ]

    non_gens_θh_post_pf =
        θh_post_pf[ non_gens_nodes_idx ]

    non_gens_vh_θh_post_pf =
        [[ a_gen_vh, a_gen_θh ]
         for ( a_gen_vh, a_gen_θh ) in
             zip( non_gens_vh_post_pf ,
                  non_gens_θh_post_pf )]

    #----------------------------------------

    gens_id_iq_post_pf =
        [ [ a_gen_id, a_gen_iq ]
         for ( a_gen_id, a_gen_iq ) in
             zip( gens_id_post_pf,
                  gens_iq_post_pf )]

    #----------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_nodes_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_nodes_δ_ed_dash_eq_dash )    
    
    #----------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------    

    id_iq_pg_vh_post_pf =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
         for (a_vh,
              a_θh,
              a_δ,
              a_ed_dash,
              a_eq_dash,
              a_ra,
              a_X_d_dash,
              a_X_q_dash ) in
                zip( gens_vh_post_pf,
                     gens_θh_post_pf,
                     gens_δ,
                     gens_ed_dash,
                     gens_eq_dash,
                     gens_ra,
                     gens_Xd_dash,
                     gens_Xq_dash ) ]
    
    #----------------------------------------
    #----------------------------------------
    
    # stateDiffCache_gen =
    #     [stateDiffCache[idx]
    #      for idx in
    #          each_gens_im_vars_Idx_in_state ]
    
    # per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     get_per_vars_or_paras_to_per_node(
    #         [ gens_nodes_ωs_ωref0_vref0_porder0,
    #          id_iq_pg_vh_post_pf ] )

    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    # update_gens_nodes_im_Ax_system_matrices!(
    #     vec_Ax_views,
    #     im_vars_in_state,
    #     each_gens_im_vars_Idx_in_state,
    #     gens_nodes_collection )


    # for (a_Ax, a_gen_state_Idx, a_gen_node ) in
    #     zip( vec_Ax_views,
    #          each_gens_im_vars_Idx_in_state,
    #          gens_nodes_collection )

    #     update_a_im_plant_Ax_system_matrices!(
    #         a_Ax, x[a_gen_state_Idx], a_gen_node )
        
    # end
        
    #----------------------------------------
    #----------------------------------------

    # per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     [[a_gen_ωs_ωref0_vref0_porder0;
    #       a_gen_id_iq_pg_vh_post_pf]
    #      for (a_gen_ωs_ωref0_vref0_porder0,
    #           a_gen_id_iq_pg_vh_post_pf) in
    #          zip( gens_nodes_ωs_ωref0_vref0_porder0,
    #               id_iq_pg_vh_post_pf) ] 
        
    # for ( dx_gen,
    #       x_gen,
    #       per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    #       ode_per_gen_model_func_kwd_para ) in

    #     zip( list_dx_gens_states_views,
    #          list_x_gens_states_views, 
    #          per_node_gens_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    #          ode_per_gen_models_func_kwd_paras)
        
    #     ode_per_gen_model_func!(
    #         dx_gen,
    #         x_gen,
    #         per_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh,
    #         t
    #         ;ode_per_gen_model_func_kwd_para =
    #             ode_per_gen_model_func_kwd_para )
        
    # end
        
    #----------------------------------------
    #----------------------------------------
    
    for ( a_gen_im_vars_Idx_in_state,
          a_id_iq_pg_vh,
          a_ωs_ωref0_vref0_porder0,
          ode_per_gen_model_func_kwd_para,
          a_gen_node ) in

        zip( each_gens_im_vars_Idx_in_state,
             id_iq_pg_vh_post_pf,
             gens_nodes_ωs_ωref0_vref0_porder0,         
             ode_per_gen_models_func_kwd_paras,
             gens_nodes_collection )
        
        a_Ax, a_Bx, a_Cx, _, _ =
            ode_per_gen_model_func_kwd_para

        # update_a_im_plant_Ax_system_matrices!(
        #     a_Ax,
        #     x[a_gen_im_vars_Idx_in_state],
        #     a_gen_node )

        (sparse_row_idxs,
         sparse_col_idxs,
         Ax_sparse_nzvalues ) =
             findnz( sparse( Matrix( a_Ax )) )
        
        Ax_sparse_nzvalues_DiffCache =
            get_tmp(Ax_sparse_nzvalues_DiffCache,
                    Ax_sparse_nzvalues )
        
        Ax_sparse_nzvalues_DiffCache .=
            Ax_sparse_nzvalues

        a_stateDiffCache =
            stateDiffCache[a_gen_im_vars_Idx_in_state]
        
        (; update_γ, γ_idx_in_sp_nzv ) =
            get_a_im_Ax_update_value_and_Idx(
                sparse_row_idxs, sparse_col_idxs,
                a_stateDiffCache, a_gen_node )

        Ax_sparse_nzvalues =
            [Ax_sparse_nzvalues[1:γ_idx_in_sp_nzv-1];
             [update_γ];
             Ax_sparse_nzvalues[γ_idx_in_sp_nzv+1:end] ]
        
        updated_Ax =
            sparse( sparse_row_idxs,
                    sparse_col_idxs,
                    Ax_sparse_nzvalues)
        

            dx[a_gen_im_vars_Idx_in_state] .=
            updated_Ax * x[a_gen_im_vars_Idx_in_state] +
            a_Bx * a_id_iq_pg_vh +
            a_Cx * a_ωs_ωref0_vref0_porder0
        
    end

    """
    #----------------------------------------
    # all nodes
    #----------------------------------------

    for ( a_vtf_dx, a_vtf_x, a_node_ur_ui_post_pf ) in
        zip( list_vtf_dx_nodes_u_views,
             list_vtf_x_nodes_u_views,
             nodes_ur_ui_post_pf )

        a_voltage_terminal_by_ur_ui_func!(
            a_vtf_dx, a_vtf_x, a_node_ur_ui_post_pf, t)

    end

    """
    #----------------------------------------
    #----------------------------------------
    

    vtf_gens_vh_θh_δ_ed_dash_eq_dash =
        get_per_vars_or_paras_to_per_node(
            [ gens_vh_θh_post_pf,
              gens_nodes_δ_ed_dash_eq_dash ] )
    
    #----------------------------------------    
    # vtf gen nodes
    #----------------------------------------

    for (vtf_gen_dx, vtf_gen_x,
         gen_vh_θh_δ_ed_dash_eq_dash,
         vtf_kwd_para) in zip(
             vtf_dx_gens_u_views,
             vtf_x_gens_u_views,
             vtf_gens_vh_θh_δ_ed_dash_eq_dash,
             vtf_gens_fun_kwd_para )

        
        a_gen_voltage_terminal_by_vh_θh_func!(
            vtf_gen_dx,
            vtf_gen_x,
            gen_vh_θh_δ_ed_dash_eq_dash,
            t
            ; vtf_kwd_para =
                vtf_kwd_para )
    end

    #----------------------------------------    
    # vtf non gen nodes
    #----------------------------------------

    for (vtf_non_gen_du, vtf_non_gen_u, non_gen_vh_θh ) in
        zip(
            vtf_dx_non_gens_u_views,
            vtf_x_non_gens_u_views,
            non_gens_vh_θh_post_pf )

         a_non_gen_voltage_terminal_by_vh_θh_func!(
             vtf_non_gen_du,
             vtf_non_gen_u,
             non_gen_vh_θh,
             t )
    end

    
    return nothing
    
end



#-----------------------------------------------------
#-----------------------------------------------------
# Simulation functions
#-----------------------------------------------------
#-----------------------------------------------------


function simulate_a_dynamics_industrial_model(
    ; dynamics_case = dynamics_case,
    only_gen = false,
    sim_timespan  = sim_timespan,
    ode_alg       = ode_alg,
    dyn_global_pf_options = dyn_global_pf_options,
    sta_global_pf_options = sta_global_pf_options )

    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    #-----------------------------------------------------
    #----------------------------------------------------- 
    
    # Dyn_Nodes, Dyn_Branches = dynamics_case()

    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------
    #-----------------------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #-----------------------------------------------------

    # (network_bus_names,
    #  non_gens_bus_names,
    #  gens_bus_names,
    #  Loads_bus_names,
    #  Trans_bus_names) =
    #      make_case_buses_names(
    #          ; case_fun = dynamics_case )

    # (net_bus_volts_labels,
    #  gens_nodes_pure_states_labels,
    #  gens_nodes_stab_states_label,
    #  gens_nodes_algebraic_and_states_labels) =
    #     generate_industrial_model_labels(
    #         network_bus_names,
    #         gens_nodes_collection )

    # (industrial_model_sym,
    #  industrial_model_mass_matrix) =
    #     generate_industrial_model_sym_and_mass_matrix(
    #         gens_nodes_pure_states_labels,
    #         net_bus_volts_labels)
    
    #-----------------------------------------------------
    
    (network_bus_names,
     non_gens_bus_names,
     gens_bus_names,
     Loads_bus_names,
     Trans_bus_names) =
         make_case_buses_names( netd.nodes )

    net_class_names =
        (;
         network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-----------------------------------------------------

    (net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_stab_states_label,
     gens_nodes_algebraic_and_states_labels) =
        generate_industrial_model_labels(
            ; nodes =  netd.nodes )

    net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_stab_states_label,
         gens_nodes_algebraic_and_states_labels )
    
    #-----------------------------------------------------

    (industrial_model_sym,
     industrial_model_mass_matrix) =
        generate_industrial_model_sym_and_mass_matrix(
            ; nodes = netd.nodes )

    #-----------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         net_states_and_var_labels,
         industrial_model_sym )
    
    #-----------------------------------------------------
    
    industrial_model_pf_sys_param_sys_views_sys_industrial_model_init =
        get_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init(
            netd; sta_global_pf_options... )

    #-----------------------------------------------------
    #-----------------------------------------------------

    (nodes_cb_sw, state, global_pf_param,
     named_tup_pf_result,
     industrial_model_misc_Idx,
     para_update_gen_Ax_aux,
     industrial_model_para_aux_inputs,
     industrial_model_pf_para,
     industrial_model_ωs_τm_vref_vhθh_idq,
     industrial_model_dyn_pf_up_para,
     industrial_model_idq_pf_cal,
     gens_nodes_τm_vf)  =
         industrial_model_pf_sys_param_sys_views_sys_industrial_model_init

    
    # -----------------------------------------------------

    pf_net_param, sd_pf_views, _ = global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ = pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views = sd_pf_views

    
    #-----------------------------------------------------
    # Stability
    #-----------------------------------------------------  
    #-----------------------------------------------------


    if only_gen == false
        Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs =
            create_gens_nodes_aggregate_system_matrices_Idx_and_views(  gens_nodes_collection; only_gen = only_gen )

        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views)  = Ax_Bx_Cx_views
        
        (Ax_matrix,
         Bx_matrix,
         Cx_matrix) = Ax_Bx_Cx_matrix

        (nodes_state_Idx,
         Bx_idxs,
         Cx_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen =
                only_gen,
            vec_Ax_τm_vf_views =
                nothing )
        
    else
        (Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
            create_gens_nodes_aggregate_system_matrices_Idx_and_views(  gens_nodes_collection; only_gen = only_gen )
        
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         vec_Ax_τm_vf_views) =
             Ax_Bx_Cx_views
        
        (Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix) =
            Ax_Bx_Cx_matrix
        
        (nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         τm_vf_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen =
                only_gen,
            vec_Ax_τm_vf_views =
                vec_Ax_τm_vf_views )
        
    end
    
    #-----------------------------------------------------
    #-----------------------------------------------------

    (gens_dynamic_id_iq_pg_vh_by_vhθh,
     gens_nodes_ωs_τm_vref_porder_view) =
         industrial_model_pf_para

    
    #-----------------------------------------------------

    ode_fun_para =
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         industrial_model_pf_para  )
    
    para_model_fun = (vec_Ax_views, ode_fun_para )

    #-----------------------------------------------------

    sim_state_x0  = state
    
    sim_timespan  = sim_timespan

    #-----------------------------------------------------
    
    dyn_global_pf_options = dyn_global_pf_options
    
    counter_array  = [1]
    
    stateDiffCache = similar( sim_state_x0 )
    
    sim_fun_para   =
        (; netd,
         nodes_cb_sw,
         global_pf_param,
         counter_array,
         state,
         stateDiffCache,
         dyn_global_pf_options,
         para_update_gen_Ax_aux,
         para_model_fun,
         industrial_model_misc_Idx,
         industrial_model_para_aux_inputs,
         industrial_model_dyn_pf_up_para,
         industrial_model_idq_pf_cal  )
    
    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------
    #-----------------------------------------------------

    sim_func  =
        dynamics_industrial_model!
    
    #---------------------------------------------------    
        
    sim_ode_func! = ODEFunction{true}(
        sim_func;
        mass_matrix =
            industrial_model_mass_matrix,
        syms =
            industrial_model_sym )
    
    sim_prob = ODEProblem(
        sim_ode_func!,
        sim_state_x0,
        sim_timespan,
        sim_fun_para )

    sim_sol =
        DifferentialEquations.solve(sim_prob, ode_alg )

    return (; sim_sol, para_net_names_labels_syms) 
      
end


function simulate_a_dynamics_industrial_model_with_or_no_controllers(
    ; dynamics_case = nothing,
    only_gen = nothing,
    sim_timespan  = nothing,
    ode_alg = ode_alg,
    dyn_global_pf_options = nothing ,
    sta_global_pf_options = nothing,
    with_cb = false)

    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    netd  = NetworkData( dynamics_case()... )
    
    #---------------------------------------------------    

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #---------------------------------------------------
    
    (network_bus_names, non_gens_bus_names,
     gens_bus_names, Loads_bus_names,
     Trans_bus_names) =
        make_case_buses_names(  netd.nodes )
    
    net_class_names =
        (; network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    (net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_stab_states_label,
     gens_nodes_algebraic_and_states_labels) =
        generate_industrial_model_labels(
            ; nodes =
                netd.nodes,
            no_control_device =
                only_gen )

    net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_stab_states_label,
         gens_nodes_algebraic_and_states_labels )

    #-----------------------------------------------------

    (industrial_model_sym,
     industrial_model_mass_matrix) =
        generate_industrial_model_sym_and_mass_matrix(
            ; nodes =
                netd.nodes,
            no_control_device =
                only_gen )

    #-----------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         net_states_and_var_labels,
         industrial_model_sym )

    #-----------------------------------------------------

    industrial_model_pf_sys_param_sys_views_sys_industrial_model_init =
        get_industrial_model_pf_param_views_and_init_with_or_no_controller(
            netd; sta_global_pf_options...,
            only_gen=only_gen )

    #-----------------------------------------------------

    (nodes_cb_sw, state,
     global_pf_param,
     named_tup_pf_result,
     industrial_model_misc_Idx,
     para_update_gen_Ax_aux,
     industrial_model_para_aux_inputs,
     industrial_model_pf_para,
     industrial_model_ωs_τm_vref_vhθh_idq ,
     industrial_model_dyn_pf_up_para,
     industrial_model_idq_pf_cal)  =
         industrial_model_pf_sys_param_sys_views_sys_industrial_model_init


    # -----------------------------------------------------

    pf_net_param, sd_pf_views, _ = global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ = pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views = sd_pf_views

    #-----------------------------------------------------
    #-----------------------------------------------------


    if only_gen == false
        
        #-------------------------------------------------
        # Stability
        #------------------------------------------------- 

        (;Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
            create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        vec_Ax_views, vec_Bx_views, vec_Cx_views  =
            Ax_Bx_Cx_views
        
        Ax_matrix, Bx_matrix, Cx_matrix =
            Ax_Bx_Cx_matrix

        nodes_state_Idx, Bx_idxs, Cx_idxs =
            idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views, vec_Bx_views, vec_Cx_views),
            gens_nodes_collection; only_gen =
                only_gen, vec_Ax_τm_vf_views = nothing )

        #-------------------------------------------------
        #-------------------------------------------------

        ode_fun_para =
            (; Ax_matrix, Bx_matrix,
             Cx_matrix,
             industrial_model_pf_para  )

        para_model_fun =
            (; vec_Ax_views, ode_fun_para )

    else
        
        #------------------------------------------------
        # Stability
        #-------------------------------------------------  
        #-------------------------------------------------
        
        (;Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
            create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        (;vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         vec_Ax_τm_vf_views) =
             Ax_Bx_Cx_views
        
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix) =
             Ax_Bx_Cx_matrix
        
        (;nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         τm_vf_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen =
                only_gen,
            vec_Ax_τm_vf_views =
                vec_Ax_τm_vf_views )

        #-------------------------------------------------

        # update_gens_τm_vf!(
        #     vec_τm_vf_views, gens_nodes_τm_vf )

        #-------------------------------------------------

        ode_fun_para =
            (; Ax_matrix, Bx_matrix,
             Cx_matrix, Ax_τm_vf_matrix,
             industrial_model_pf_para  )

        para_model_fun =
            (; vec_Ax_views,
             ode_fun_para )

    end
    
    #-----------------------------------------------------

    sim_state_x0  = state

    sim_timespan  = sim_timespan

    #-----------------------------------------------------

    dyn_global_pf_options = dyn_global_pf_options

    counter_array  = [1]

    stateDiffCache = similar( sim_state_x0 )

    sim_fun_para = (
        ; netd,
        nodes_cb_sw,
        global_pf_param,
        counter_array,
        state,
        stateDiffCache,
        dyn_global_pf_options,
        para_update_gen_Ax_aux,
        para_model_fun,
        industrial_model_misc_Idx,
        industrial_model_para_aux_inputs,        
        industrial_model_dyn_pf_up_para,
        industrial_model_idq_pf_cal,                
        only_gen )

    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------

    sim_func =
        dynamics_industrial_model_with_or_no_controller!

    # dynamics_industrial_model!

    #---------------------------------------------------    

    sim_ode_func! = ODEFunction{true}(
        sim_func;
        mass_matrix =
            industrial_model_mass_matrix,
        syms =
            industrial_model_sym )

    sim_prob = ODEProblem(sim_ode_func!,
                          sim_state_x0,
                          sim_timespan,
                          sim_fun_para )

    if with_cb == true
        
        cb = industrial_model_fun_make_state_callbacks(
            collect( values( netd.nodes)) )

        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg, callback = cb )
        
    else
        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg  )

    end
    
    return (;sim_sol, para_net_names_labels_syms) 
      
end



function v2_simulate_a_dynamics_im_model(
    ; dynamics_case = dynamics_case,
    only_gen = false,
    sim_timespan  = sim_timespan,
    ode_alg       = ode_alg,
    dyn_global_pf_options = dyn_global_pf_options,
    sta_global_pf_options = sta_global_pf_options,
    with_cb = false
        )

    # ------------------------------------------------
    # ------------------------------------------------
    # Diagnosis
    # ------------------------------------------------ 
    # -------------------------------------------------

     # dynamics_case =
     #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    #--------------------------------------------------
    #-------------------------------------------------- 
    
    # Dyn_Nodes, Dyn_Branches = dynamics_case()

    netd  = NetworkData( dynamics_case()... )
    
    #-------------------------------------------------
    #-------------------------------------------------

    non_gens_nodes = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #-------------------------------------------------
    #-------------------------------------------------
    
    (;network_bus_names,
     non_gens_bus_names,
     gens_bus_names,
     Loads_bus_names,
     Trans_bus_names) =
         make_case_buses_names(
             netd.nodes )
    
    net_class_names =
        (; network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-------------------------------------------------

    (;net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
        generate_im_model_labels(
            ; nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #--------------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes = netd.nodes )

    #-------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    
    #-------------------------------------------------
    
    im_sys_pf_param_views_and_init =
        v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )

    #-------------------------------------------------
    #-------------------------------------------------

    (; nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para,
     im_idq_pf_cal)  =
        im_sys_pf_param_views_and_init
    
    # -------------------------------------------------

    pf_net_param, sd_pf_views, _ =
        global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx =
        pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
        pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
        sd_pf_views
    
    #--------------------------------------------------
    # Stability
    #--------------------------------------------------  

    (Ax_Bx_Cx_views,
     Ax_Bx_Cx_matrix,
     idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix =
        Ax_Bx_Cx_matrix

    (nodes_state_Idx,
     Bx_idxs, Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs) =
        idxs

    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )
    
    #-------------------------------------------------
    #-------------------------------------------------

    ode_fun_para = (
        ; Ax_matrix,
        Bx_matrix,
        Cx_matrix,
        im_pf_para )
    
    para_model_fun = (
        ; vec_Ax_views,
        ode_fun_para )

    #-------------------------------------------------

    sim_state_x0  = state
    
    sim_timespan  = sim_timespan

    #-------------------------------------------------
    
    dyn_global_pf_options = dyn_global_pf_options
    
    counter_array  = [1]
    
    stateDiffCache = similar( sim_state_x0 )
    
    sim_fun_para = (
        ; netd,
        nodes_cb_sw,
        global_pf_param,
        counter_array,
        state,
        stateDiffCache,
        dyn_global_pf_options,
        para_update_gen_Ax_aux,
        para_model_fun,
        im_misc_Idx,
        im_para_aux_inputs,
        im_dyn_pf_up_para,
        im_idq_pf_cal )
    
    #-----------------------------------------------
    # dynamics simulation
    #-----------------------------------------------
    #-----------------------------------------------

    sim_func = dynamics_im_model!
    
    #-----------------------------------------------
    
    #-----------------------------------------------
        
    sim_ode_func! = ODEFunction{true}(
        sim_func; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
    sim_prob = ODEProblem(sim_ode_func!,
                          sim_state_x0,
                          sim_timespan,
                          sim_fun_para )


    if with_cb == true
        
        cb = industrial_model_fun_make_state_callbacks(
            collect( values( netd.nodes)) )

        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg, callback = cb )
        
    else
        
        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg )
    end
    
    #-------------------------------------------------
    #-------------------------------------------------
    
    return sim_sol, para_net_names_labels_syms 
      
end



function t_v2_simulate_one_plant_dynamics_im_model(
    ; dynamics_case =
        dynamics_case,
    only_gen =
        false,
    sim_timespan =
        sim_timespan,
    ode_alg =
        ode_alg,
    dyn_global_pf_options =
        dyn_global_pf_options,
    sta_global_pf_options =
        sta_global_pf_options,
    with_cb =
        false
        )

    # ------------------------------------------------
    # Diagnosis
    # ------------------------------------------------ 

    # ImplicitMidpoint(), ImplicitMidpoint(autodiff=false),
    
    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)

    sim_model_type        = "im-model"
    
    # ode_alg             = Rodas4()
    # algr_name           = "rodas4" 
    
    ode_alg               = ImplicitMidpoint()    
    algr_name             = "ImplicitMidpoint"

    dt                    = 0.01

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = false

    with_cb  = false

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
        
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #-----------------------------------------------------
    
    ; (network_bus_names,
       non_gens_bus_names,
       gens_bus_names,
       Loads_bus_names,
       Trans_bus_names) =
           make_case_buses_names(
               netd.nodes )
    
    net_class_names =
        (;network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-------------------------------------------------

    (; net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #-------------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes = netd.nodes )

    #-------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    
    #-------------------------------------------------

    (; nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para, im_idq_pf_cal)  =
         v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )
    
    #-------------------------------------------------
    
    (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para

    gens_vh_θh_view =
        im_para_aux_inputs.gens_vh_θh_view
    
    # --------------------------------------------------

    pf_net_param, sd_pf_views, _ =
        global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx =
        pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
        pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
        sd_pf_views

    
    #-------------------------------------------------
    # Stability
    #-------------------------------------------------  
    
    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, states_and_mat_Idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    (; Ax_matrix,
     Bx_matrix,
     Cx_matrix) = Ax_Bx_Cx_matrix

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs)  =
         states_and_mat_Idxs
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )

    #-------------------------------------------------
    
    node_i = 1
    
    gen_node_i_idx =
        nodes_state_Idx[ node_i ]
    
    #--------------------------------------------------

    Ax_m = vec_Ax_views[ node_i ]

    Bx_m = vec_Bx_views[ node_i ]

    Cx_m = vec_Cx_views[ node_i ]
    
    #--------------------------------------------------
    
    plant_i =
        gens_nodes_collection[ node_i ]      

    #--------------------------------------------------

    plant_i_eigvalues =
        eigvals( Matrix( Ax_m ))

    #--------------------------------------------------

    ra_Xd_dash_Xq_dash_i =
        ra_Xd_dash_Xq_dash_view[ node_i ]

    #---------------------------------------------------
    
    id_iq_pg_vh_i =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view[ node_i ]

    ωs_ωref0_vref0_porder0_i =
        gens_nodes_ωs_ωref0_vref0_porder0_view[ node_i ]

    gens_vh_θh_i = [ gens_vh_θh_view[ node_i ]...;]

    #-------------------------------------------------

    dims_gens_vh_θh_ωs_ωref0_vref0_porder0_i =
        length.([ gens_vh_θh_i, ωs_ωref0_vref0_porder0_i ])

    _,_, vh_θh_ωs_ωref0_vref0_porder0_i_Idx =
        create_size_offset_Idx(
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0_i,;
            counter = 0)

    vh_θh_i_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_i_Idx[1]

    ωs_ωref0_vref0_porder0_i_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_i_Idx[2]
        
    kwd_dyn_vh_θh_ωs_ωref0_vref0_etc_i =
        (; gens_vh_θh_i,
         ωs_ωref0_vref0_porder0_i )

    #--------------------------------------------------
    
    poi_dyn_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )

    pois_dyn =
        kwd_dyn_vh_θh_ωs_ωref0_vref0_etc_i    

    #-------------------------------------------------

    sim_fun_para_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )
        
    sim_fun_para  =
        vcat( gens_vh_θh_i,
             ωs_ωref0_vref0_porder0_i )
    
    #--------------------------------------------------

    sim_state_x0 =
        state[ gen_node_i_idx ]

    chunk_size = length(sim_state_x0 )
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0),
    #         chunk_size )

    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0) )

    stateDiffCache = similar(sim_state_x0)

    #--------------------------------------------------

    """
    poi_dyn_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )

    pois_dyn =
        (; gens_vh_θh_i,
         ωs_ωref0_vref0_porder0_i )
           
    """

    s_poi =  s_oop_poi = pois_dyn.ωs_ωref0_vref0_porder0_i
    
    poi_idx = poi_dyn_Idxs.ωs_ωref0_vref0_porder0_i_Idx

    poi_DiffCache = DiffCache(similar( s_poi ))

    #-------------------------------------------------
    
    sim_timespan =
        sim_timespan

    #---------------------------------------------------

    s_sim_fun_kwd_para =
        (; Ax_m,
         Bx_m,
         Cx_m,
         stateDiffCache,
         poi_DiffCache,
         sim_state_x0,
         plant_i,
         ra_Xd_dash_Xq_dash_i,
         pois_dyn,
         poi_dyn_Idxs )
    
    t_sim_fun_kwd_para =
        (; Ax_m,
         Bx_m,
         Cx_m,
         stateDiffCache,
         sim_state_x0,
         plant_i,
         sim_fun_para_Idxs,
         ra_Xd_dash_Xq_dash_i )
    
    #------------------------------------------------
    
    #--------------------------------------------------
    # dynamics simulation
    #--------------------------------------------------
            
    im_sym_i =
        im_sym[ gen_node_i_idx ]

    im_mass_matrix_i =
        im_mass_matrix[gen_node_i_idx,
                       gen_node_i_idx]

    #---------------------------------------------------
    # case 1
    #---------------------------------------------------

    t_sim_func! = t_ode_one_im_model_func!
    
    #-------------------------------------------------   
        
    t_sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> t_sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para = t_sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix_i,
        syms = im_sym_i )
    
    t_sim_prob = ODEProblem(
        t_sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )
    
    #--------------------------------------------------
    # simulate 
    #--------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    t_sim_sol = DifferentialEquations.solve(
        t_sim_prob, ode_alg, dt=dt  )

    #-------------------------------------------------
    # plot
    #-------------------------------------------------
    
    bus_i_name = network_bus_names[ node_i ]

    t_plot_δ_ω_ed_dash_eq_dash =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = t_sim_sol,
        node_syms_labels = im_sym_i,
        bus_name = bus_i_name,
        vars = [ :δ, :ω, :ed_dash, :eq_dash ],
        tspan = (0.0, 10.0),
        fmt = :png )

    p_δ_ω_ed_dash_eq_dash_idxs =
        last.(get_a_node_state_algb_vars_indices_in_syms(
            ; node_syms_labels = im_sym_i,
            bus_name = bus_i_name,
            vars = [ :δ, :ω, :ed_dash, :eq_dash ] ))


    #-------------------------------------------------
    # sensitivity
    #-------------------------------------------------

    # case 1
    
    t_sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            t_sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para = t_sim_fun_kwd_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg = ForwardDiffSensitivity() )


    t_sa_sim_sol = DifferentialEquations.solve(
        t_sens_prob , ode_alg, dt=dt  )


    #--------------------------------------------------
    # case 2
    #--------------------------------------------------

    s_sim_func! = s_ode_one_im_model_func!
    
    #-------------------------------------------------   
        
    s_sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> s_sim_func!(
            dx, x, p, t;
            poi_idx = poi_idx,
            sim_fun_kwd_para = s_sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix_i,
        syms = im_sym_i )
    
    s_sim_prob = ODEProblem(
        s_sim_ode_func,
        sim_state_x0,
        sim_timespan,
        s_poi  )
    
    #-------------------------------------------------
    # simulate 
    #--------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    s_sim_sol = DifferentialEquations.solve(
        s_sim_prob, ode_alg, dt=dt  )

    #--------------------------------------------------
    # plot
    #--------------------------------------------------
    
    bus_i_name = network_bus_names[ node_i ]

    s_plot_δ_ω_ed_dash_eq_dash =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = s_sim_sol,
        node_syms_labels = im_sym_i,
        bus_name = bus_i_name,
        vars = [ :δ, :ω, :ed_dash, :eq_dash ],
        tspan = (0.0, 10.0),
        fmt = :png )

    p_δ_ω_ed_dash_eq_dash_idxs =
        last.(get_a_node_state_algb_vars_indices_in_syms(
            ; node_syms_labels = im_sym_i,
            bus_name = bus_i_name,
            vars = [ :δ, :ω, :ed_dash, :eq_dash ] ))


    #-------------------------------------------------
    # sensitivity
    #-------------------------------------------------

    # case 2
    
    s_sa_sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            s_sim_func!(
                dx, x, p, t;
                poi_idx = poi_idx,
                sim_fun_kwd_para = s_sim_fun_kwd_para ),
        sim_state_x0,
        sim_timespan,
        s_poi; sensealg = ForwardDiffSensitivity() )


    s_sa_sim_sol = DifferentialEquations.solve(
        s_sa_sens_prob , ode_alg, dt=dt  )

    #--------------------------------------------------
    #--------------------------------------------------

    """

    https://stackoverflow.com/questions/74653454/why-am-i-getting-a-mutating-arrays-is-not-supported-error-here

    https://discourse.julialang.org/t/sensitivities-with-respect-to-initial-conditions-in-differentialequations-jl/25555/12
    https://github.com/FluxML/Tracker.jl

    https://discourse.julialang.org/t/discrete-adjoint-sensitivity-analysis-for-odes-in-differentialequations-jl/100007/3
    https://docs.sciml.ai/Overview/dev/highlevels/array_libraries/

    """
    
    return nothing

      
end
                         


function t_v2_simulate_sysyem_dynamic_im_model(
    ; dynamics_case =
        dynamics_case,
    only_gen =
        false,
    sim_timespan =
        sim_timespan,
    ode_alg =
        ode_alg,
    dyn_global_pf_options =
        dyn_global_pf_options,
    sta_global_pf_options =
        sta_global_pf_options,
    with_cb =
        false
        )

    # ------------------------------------------------
    # Diagnosis
    # ------------------------------------------------ 

    # ImplicitMidpoint(), ImplicitMidpoint(autodiff=false),
    
    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)

    sim_model_type        = "im-model"
    
    # ode_alg             = Rodas4()
    # algr_name           = "rodas4" 
    
    ode_alg               = ImplicitMidpoint()    
    algr_name             = "ImplicitMidpoint"

    dt                    = 0.01

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = false

    with_cb  = false

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
        
    #--------------------------------------------------
    #--------------------------------------------------
    
    netd  = NetworkData( dynamics_case()... )
    
    #-------------------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #--------------------------------------------------
    
    (; network_bus_names,
       non_gens_bus_names,
       gens_bus_names,
       Loads_bus_names,
       Trans_bus_names) =
           make_case_buses_names(
               netd.nodes )
    
    net_class_names =
        (;network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #--------------------------------------------------

    (; net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #-------------------------------------------------
    
    # im_sym, im_mass_matrix =
    #     generate_im_sym_and_mass_matrix(
    #         ; nodes = netd.nodes )

    (; im_model_ode_mass_matrix,
     im_model_pf_mass_matrix) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes = netd.nodes )

    (; gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes = netd.nodes  )

    im_sym =
        gens_nodes_im_vars_labels

    im_mass_matrix =
        im_model_ode_mass_matrix
        
    #-------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    
    #--------------------------------------------------

    (; nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para, im_idq_pf_cal)  =
         v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )
    
    #-------------------------------------------------
    
    (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para

    gens_vh_θh_view =
        im_para_aux_inputs.gens_vh_θh_view
    
    # ------------------------------------------------

    pf_net_param, sd_pf_views, _ =
        global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx =
        pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
        pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
        sd_pf_views


    (; each_gens_im_vars_Idx_in_state,
     im_vars_Idx_in_state,
     im_vars_view_in_state ) =
         para_update_gen_Ax_aux
    
    #-------------------------------------------------
    # Stability
    #-------------------------------------------------  
    
     (Ax_Bx_Cx_views,
      Ax_Bx_Cx_matrix,
      states_and_mat_Idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

     (nodes_state_Idx, Bx_idxs, Cx_idxs,
      id_iq_ph_vh_idxs, ω_ref_ωref_v_ref_idxs)  =
          states_and_mat_Idxs

    sys_states_Idxs_and_mat_Idxs =
        (; nodes_state_Idx,
         Bx_idxs, Cx_idxs,
         id_iq_ph_vh_idxs,
         ω_ref_ωref_v_ref_idxs )
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )

    #------------------------------------------------

    Ax_m = vec_Ax_views

    Bx_m = vec_Bx_views

    Cx_m = vec_Cx_views
    
    #-----------------------------------------------    

    # plant_i_eigvalues =
    #     eigvals( Matrix( vec_Ax_views[i] ))

    #-----------------------------------------------
    #-----------------------------------------------

    """

    #------------------------------------------------
    # case 0
    #------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i]) for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ]

    gens_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gens_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #--------------------------------------------------
    
   dims_vh_θh_ωs_ωref0_vref0_porder0 =  length.( [ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ] )

    system_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), [dims_vh_θh_ωs_ωref0_vref0_porder0]  )

    vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(system_size_offset_Idx)
    
    sim_fun_system_para_Idxs =
        vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #--------------------------------------------------

    system_ode_para = [[ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ]...;]
    
    #--------------------------------------------------

    sim_state_x0 =
        state[im_vars_Idx_in_state]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0),
    #         chunk_size )
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0) )

    #-------------------------------------------------

    stateDiffCache_gens =
        [ stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]

    #--------------------------------------------------
    
    sim_fun_system_kwd_para =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection,
         gens_sim_fun_gens_para_Idxs,
         sim_fun_system_para_Idxs,         
         ra_Xd_dash_Xq_dash_view,
         sys_states_Idxs_and_mat_Idxs,
         para_update_gen_Ax_aux)

    #-------------------------------------------------
    
    sim_func! = system_ode_im_model_func!

    sim_fun_para = system_ode_para
    
    #------------------------------------------------- 
        
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para = sim_fun_system_kwd_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )
    
    #-------------------------------------------------
    # simulate 
    #-------------------------------------------------

    # system_sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, ode_alg, dt=dt  )

    #-------------------------------------------------
    # sensitivity
    #-------------------------------------------------
    
    sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para =
                    sim_fun_system_kwd_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg = ForwardDiffSensitivity() )


    sa_sim_sol = DifferentialEquations.solve(
        sens_prob , ode_alg, dt=dt  )

    #-------------------------------------------------
    # plot
    #--------------------------------------------------

    plot_idxs_a = get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = system_sim_sol,
        node_syms_labels = im_sym,
        bus_name = "bus1",
        vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    # plot_a = make_plot_of_buses_volt_mag(
    #     ; sol = system_sim_sol,
    #     network_vars_labels = im_sym,
    #     nodes_name = ["bus1", "bus2"],
    #     vars = [:u_r, :u_i],
    #     tspan = sim_timespan,
    #     fmt = :png)
    

    # plot_b = make_plot_of_buses_volt_angle(
    #     ; sol = system_sim_sol,
    #     network_vars_labels = im_sym,
    #     nodes_name = ["bus1", "bus2"],
    #     vars = [:u_r, :u_i],
    #     tspan = sim_timespan,
    #     fmt = :png)

    """

    #-------------------------------------------------
    # case 1
    #-------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i]) for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ]

    gens_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gen_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-------------------------------------------------
    
    
    id_iq_pg_vh =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view

    ωs_ωref0_vref0_porder0 =
        gens_nodes_ωs_ωref0_vref0_porder0_view
        
    gens_vh_θh = gens_vh_θh_view
    
    #--------------------------------------------------

    dims_nodes_id_iq_pg_vh =
        length.( id_iq_pg_vh )

    _,_, nodes_id_iq_pg_vh_idx_in_Idx =
        create_size_offset_Idx(
            dims_nodes_id_iq_pg_vh;
        counter = 0)

    #------------------------------------------------

    dims_nodes_ωs_ωref0_vref0_porder0 =
        length.( ωs_ωref0_vref0_porder0 )

    _,_, nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        create_size_offset_Idx(
            dims_nodes_ωs_ωref0_vref0_porder0;
        counter = 0)

    #-------------------------------------------------

    dims_nodes_vh_θh =
        length.( gens_vh_θh )

    _,_, nodes_vh_θh_indx_in_Idx =
        create_size_offset_Idx(
            dims_nodes_vh_θh;
        counter = 0)
    
    #--------------------------------------------------

    f_gens_vh_θh =
        [ gens_vh_θh...;]
    
    f_ωs_ωref0_vref0_porder0 =
        [ωs_ωref0_vref0_porder0...;]
    
    dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =
        length.(
            [f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0])

    _,_, vh_θh_ωs_ωref0_vref0_porder0_Idx =
        create_size_offset_Idx(
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0;
            counter = 0)

    vh_θh_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[1]

    ωs_ωref0_vref0_porder0_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[2]
        
    kwd_dyn_vh_θh_ωs_ωref0_vref0_etc =
        (; f_gens_vh_θh,
         f_ωs_ωref0_vref0_porder0 )

    #------------------------------------------------

    sim_fun_system_ode_flat_agg_para_Idxs =
        (; vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx )
        
    sim_fun_system_ode_flat_agg_para  =
        vcat( f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0 )
    

    #--------------------------------------------------

    sim_state_x0 =
        state[im_vars_Idx_in_state]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)

    #-------------------------------------------------

    stateDiffCache_gens =
        [stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]
    
    #-------------------------------------------------
    
    # sim_fun_system_kwd_flat_agg_para =
    #     (;
    #      vec_Ax_views,
    #      vec_Bx_views,
    #      vec_Cx_views,
    #      stateDiffCache_gens,
    #      sim_state_x0_gens,
    #      gens_nodes_collection ,
    #      gens_sim_fun_gen_para_Idxs,
    #      sim_fun_system_ode_flat_agg_para_Idxs,         
    #      ra_Xd_dash_Xq_dash_view,
         
    #      sys_states_Idxs_and_mat_Idxs,         
    #      nodes_id_iq_pg_vh_idx_in_Idx,
         
    #      nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
    #      nodes_vh_θh_indx_in_Idx )

    
    gen_nodes_ra_Xd_dash_Xq_dash_view =
        ra_Xd_dash_Xq_dash_view,

    gens_nodes_id_iq_pg_vh_idx_in_Idx =
        nodes_id_iq_pg_vh_idx_in_Idx
    
    gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        nodes_ωs_ωref0_vref0_porder0_idx_in_Idx
    
    gens_nodes_vh_θh_indx_in_Idx  =
        nodes_vh_θh_indx_in_Idx   
    
    sim_fun_system_kwd_flat_agg_para =
        (;
         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,

         gen_nodes_ra_Xd_dash_Xq_dash_view,
         
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection ,
         gens_sim_fun_gen_para_Idxs,
         
         sim_fun_system_ode_flat_agg_para_Idxs, 
         sys_states_Idxs_and_mat_Idxs,         
         gens_nodes_id_iq_pg_vh_idx_in_Idx,
         gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         gens_nodes_vh_θh_indx_in_Idx )
    
    #-------------------------------------------------
    
    sim_func! =
        system_flat_agg_ode_im_model_func!

    sim_fun_para =
        sim_fun_system_ode_flat_agg_para
    
    #--------------------------------------------------
        
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para =
                sim_fun_system_kwd_flat_agg_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )

    #-------------------------------------------------
    # simulate 
    #-------------------------------------------------

    # system_sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg )

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, ode_alg, dt=dt  )

    #--------------------------------------------------
    # sensitivity
    #--------------------------------------------------
    
    sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para =
                    sim_fun_system_kwd_flat_agg_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg =
            ForwardDiffSensitivity() )

    sim_sol = DifferentialEquations.solve(
        sens_prob, ode_alg, dt=dt  )
    
    #---------------------------------------------------
    # plot
    #--------------------------------------------------

    plot_idxs_a = get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = system_sim_sol,
        node_syms_labels = im_sym,
        bus_name = "bus1",
        vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    #-------------------------------------------------
    # case 2
    #-------------------------------------------------
    
    poi_dyn_Idxs =
        (; vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx )

    pois_dyn =
        (; f_gens_vh_θh,
         f_ωs_ωref0_vref0_porder0 )   

    #-------------------------------------------------
    
    sim_fun_system_kwd_para_wt_poi =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views, 
         stateDiffCache,
         poi_DiffCache,
         sim_state_x0,
         gens_nodes_collection ,
         ra_Xd_dash_Xq_dash,
         pois_dyn,
         poi_dyn_Idxs,
         states_and_mat_Idxs,         
         nodes_id_iq_pg_vh_idx_in_Idx,
         nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         nodes_vh_θh_indx_in_Idx
)
    
    #------------------------------------------------

    """
    poi_dyn_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )

    pois_dyn =
        (; gens_vh_θh_i,
         ωs_ωref0_vref0_porder0_i )
           
    """

    s_poi =  s_oop_poi = pois_dyn.ωs_ωref0_vref0_porder0
    
    poi_idx = poi_dyn_Idxs.ωs_ωref0_vref0_porder0_Idx

    poi_DiffCache = DiffCache(similar( s_poi ))

    #--------------------------------------------------
    
    sim_timespan =
        sim_timespan
    #----------------------------------------------
    # case 3
    #------------------------------------------------
    #------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------

    #------------------------------------------------
    # case 1
    #-------------------------------------------------

    t_sim_func! = t_ode_im_model_func!
    
    #------------------------------------------------   
        
    t_sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> t_sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para = t_sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
    t_sim_prob = ODEProblem(
        t_sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )
    
    #-------------------------------------------------
    # simulate 
    #-------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    t_sim_sol = DifferentialEquations.solve(
        t_sim_prob, ode_alg, dt=dt  )

    #------------------------------------------------
    # plot
    #------------------------------------------------
    
    bus_i_name = network_bus_names[ node_i ]

    t_plot_δ_ω_ed_dash_eq_dash =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = t_sim_sol,
        node_syms_labels = im_sym_i,
        bus_name = bus_i_name,
        vars = [ :δ, :ω, :ed_dash, :eq_dash ],
        tspan = (0.0, 10.0),
        fmt = :png )

    p_δ_ω_ed_dash_eq_dash_idxs =
        last.(get_a_node_state_algb_vars_indices_in_syms(
            ; node_syms_labels = im_sym_i,
            bus_name = bus_i_name,
            vars = [ :δ, :ω, :ed_dash, :eq_dash ] ))


    #-----------------------------------------------------
    # sensitivity
    #-----------------------------------------------------

    # case 1
    
    t_sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            t_sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para = t_sim_fun_kwd_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg = ForwardDiffSensitivity() )


    t_sa_sim_sol = DifferentialEquations.solve(
        t_sens_prob , ode_alg, dt=dt  )


    #-----------------------------------------------------
    # case 2
    #-----------------------------------------------------

    s_sim_func! = s_ode_one_im_model_func!
    
    #----------------------------------------------------   
        
    s_sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> s_sim_func!(
            dx, x, p, t;
            poi_idx = poi_idx,
            sim_fun_kwd_para = s_sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix_i,
        syms = im_sym_i )
    
    s_sim_prob = ODEProblem(
        s_sim_ode_func,
        sim_state_x0,
        sim_timespan,
        s_poi  )
    
    #-----------------------------------------------------
    # simulate 
    #-----------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    s_sim_sol = DifferentialEquations.solve(
        s_sim_prob, ode_alg, dt=dt  )

    #-----------------------------------------------------
    # plot
    #-----------------------------------------------------
    
    bus_i_name = network_bus_names[ node_i ]

    s_plot_δ_ω_ed_dash_eq_dash =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = s_sim_sol,
        node_syms_labels = im_sym_i,
        bus_name = bus_i_name,
        vars = [ :δ, :ω, :ed_dash, :eq_dash ],
        tspan = (0.0, 10.0),
        fmt = :png )

    p_δ_ω_ed_dash_eq_dash_idxs =
        last.(get_a_node_state_algb_vars_indices_in_syms(
            ; node_syms_labels = im_sym_i,
            bus_name = bus_i_name,
            vars = [ :δ, :ω, :ed_dash, :eq_dash ] ))


    #-----------------------------------------------------
    # sensitivity
    #-----------------------------------------------------

    # case 2
    
    s_sa_sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            s_sim_func!(
                dx, x, p, t;
                poi_idx = poi_idx,
                sim_fun_kwd_para = s_sim_fun_kwd_para ),
        sim_state_x0,
        sim_timespan,
        s_poi; sensealg = ForwardDiffSensitivity() )


    s_sa_sim_sol = DifferentialEquations.solve(
        s_sa_sens_prob , ode_alg, dt=dt  )

    #-----------------------------------------------------
    #-----------------------------------------------------

    """

    https://stackoverflow.com/questions/74653454/why-am-i-getting-a-mutating-arrays-is-not-supported-error-here

    https://discourse.julialang.org/t/sensitivities-with-respect-to-initial-conditions-in-differentialequations-jl/25555/12
    https://github.com/FluxML/Tracker.jl

    https://discourse.julialang.org/t/discrete-adjoint-sensitivity-analysis-for-odes-in-differentialequations-jl/100007/3
    https://docs.sciml.ai/Overview/dev/highlevels/array_libraries/

    """
    
    return nothing

      
end
                         

#-----------------------------------------------------
#-----------------------------------------------------
#  Intitialisations
#-----------------------------------------------------
#-----------------------------------------------------


"""


id_iq =
    invZ_dq(ra, X_d_dash, X_q_dash) *
    [ed_dash - vh * sin(δ - θh),
     eq_dash - vh * cos(δ - θh)]

pg =
    ed_dash * id_iq[1] + eq_dash * id_iq[2] +
    ( X_q_dash - X_d_dash ) *  *(id_iq...)

id_iq =
    get_dynamic_idq_vhθh(
        vh, θh, δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )

pg =
    get_dynamic_pg_from_id_iq(
        id_iq..., δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)


"""



#---------------------------------------------------
# Industrial
#---------------------------------------------------


function get_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init(
    netd
    ; maxiter=40,
    ftol=1000*eps(),
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false )
    
    # -----------------------------------------------
    # -----------------------------------------------

    dict_sys_to_industry =
        get_net_to_industrial_model_indices_dict(
            netd )

    (; pure_states_Idx_in_system,
     ur_ui_Idx_in_system,
     industrial_Idx,
     industrial_model_pure_states_Idx,
     industrial_model_ur_ui_Idx,
     pure_states_and_ur_ui_Idx_in_system,
     industrial_model_state_Idx,
     net_to_industrial_idx_conversion_dict) =
        get_industrial_model_indices_and_conversion_dict(
            netd  )

    # ---------------------------------------------------
        
    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_δ_ω_ed_dash_eq_dash_Idx(
            netd )

    nodes_ω_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_ω_ed_dash_eq_dash_Idx(
            netd )

    nodes_δ_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_δ_ed_dash_eq_dash_Idx(
            netd )

    #------------------------------------------------- 
    
    industrial_model_each_gen_nodes_pure_states_idx_in_state =
        get_industrial_gens_pure_states_indices_in_state(
            netd )
    
    industrial_model_each_gen_nodes_stab_states_idx_in_state =
        get_industrial_gens_stab_states_indices_in_state(
            netd )
    
    
    # ---------------------------------------------------

    nodes_cb_sw =
        get_nodes_cb_sw( netd.nodes )

    gens_nodes =
        get_gens_nodes( netd.nodes )

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #---------------------------------------------------

    gen_nodes_ra_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list =
                [ :ra, :X_d_dash, :X_q_dash ] )

    gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list =
                [ :ra, :X_d, :X_q, :X_d_dash, :X_q_dash ] )

    gen_nodes_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes; param_list =
                [ :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q,
                  :X_d_dash, :X_q_dash,
                  :T_d_dash, :T_q_dash ],
            gens_view_only = true )
    
    gen_nodes_sub_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes; param_list = [
                :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q,
                :X_d_dash,
                :X_q_dash, :X_d_2dash,
                :X_q_2dash,
                :T_d_dash, :T_q_dash,
                :T_d_2dash, :T_q_2dash ],
            gens_view_only = true )
   
    # -----------------------------------------------------
    # -----------------------------------------------------

    # pf_net_param = get_powerflow_net_parameters( netd )

    pf_net_param =
        get_industrial_model_powerflow_net_parameters(
            netd )
    
    (; pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) = pf_net_param

    (;Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net

    (; slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (; ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
        pf_param_views

    load_trans_nodes_Idx_and_vlimits =
        pf_limits

    (;slack_bus_idx,
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


    #--------------------------------------------------    
    #--------------------------------------------------

    state = zeros(
        length(
            generate_industrial_model_sym(
                ; nodes = netd.nodes ) ) )

    #-------------------------------------------------   
    
    gen_nodes_δ_ω_ed_dash_eq_dash_views =
        get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )
    
    industrial_model_pure_states_view_in_state =
        get_industrial_pure_states_view_in_state(
            state, industrial_model_pure_states_Idx )
    
    #-------------------------------------------------

    state_view = view(state, 1:length(state))

    # ----------------------------------------------
    
    # for pf flat start, ur = vh = 1,  ui = θh = 0

    # remove init_pf = true
        
    
    if init_pf == true 

        state_view[ur_idx] .=
            ones(  length( ur_IDX ))
        state_view[ui_idx] .=
            zeros( length( ui_IDX ))

    end
    
    # ----------------------------------------------

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]
        
    #-----------------------------------------------
    
    nodes_u_view  =
        [ view(state, nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ]

    # ---------------------------------------------
    
    nodes_pf_U_view  =
        [ view(pf_state , nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ] 
    
    # ---------------------------------------------
    
    uh_state =
        state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]
    
    x0_ur_ui =
        [state_view[ur_ui_idx][ ur_IDX ]...;
         state_view[ur_ui_idx][ ui_IDX ]...]
        
    x0_vh_θh =
        [abs.(uh_state)...; angle.(uh_state)...]

    working_vh_θh_view =
        view(x0_vh_θh, 1:length(x0_vh_θh))

    x0_vh_view       = @view x0_vh_θh[vh_IDX]

    x0_θh_view       = @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0      = [ red_vh_θh_0_view; ]

    mismatch         = similar( red_vh_θh_0 )
    
    # ----------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))
    
    # ---------------------------------------------
    
    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet  =
        zeros(ComplexF64, length( uh_state ))
     
    Inet_view  =
        view( Inet, 1:length( Inet ) )

    Iinj  =
        zeros(ComplexF64, length( uh_state ))

    Iinj_view =
        view(Iinj, 1:length( Iinj ))

    idq_wt_pad =
        zeros(ComplexF64, length( uh_state ))

    idq_wt_pad_view =
        view(idq_wt_pad, 1:length( uh_state ) )
    
    #----------------------------------------------

    global_pf_views =
        (; working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view )

    sd_pf_views =
        (; working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          gen_nodes_δ_ω_ed_dash_eq_dash_views )

    # ----------------------------------------------
    # global_pf param
    # ----------------------------------------------
    
    global_pf_param =
        (; pf_net_param,
          sd_pf_views,
          mismatch )
    
    branches_name  = collect(keys( netd.edges ))
    
    nodes_name = collect(keys( netd.nodes ))    

    #-----------------------------------------------
    #-----------------------------------------------

    # maxiter=40
    # ftol=1000*eps()
    # xtol=1000*eps()
    # init_pf = true
    # with_δ_ed_eq = false
    
    named_tup_pf_result =
        power_balance_powerflow(
            x0_vh_θh,
            mismatch,
            sd_pf_views,
            (nodes_name, branches_name) ,
            pf_net_param
            ; maxiter=maxiter,
            ftol=ftol ,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )

    #-----------------------------------------------

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init =
        named_tup_pf_result.branch_dict_init
    
    #---------------------------------------------    
    #---------------------------------------------

    state .=
        industrial_model_init_operationpoint(
            netd, bus_dict_init; pure = :pure )

    #--------------------------------------------- 

    # gen_nodes_δ_ω_ed_dash_eq_dash_views .= get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash( state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    
    #--------------------------------------------
    
    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)
    
    #--------------------------------------------

    # nodes_f_t  = external_get_nodes_or_edges_f_t(
    #     netd.nodes, bus_dict_init )

    #--------------------------------------------
    #-------------------------------------------
 
    # gens_vh_θh_post_pf =
    #     get_gens_vh_θh_post_pf(
    #         gens_nodes_collection ,
    #         bus_dict_init )
        
    # gens_ur_ui_post_pf =
    #     get_gens_ur_ui_post_pf(
    #         gens_nodes_collection,
    #         bus_dict_init )
    
    gens_vh_θh =
        get_gens_vh_θh(
            nodes_pf_U_view, gens_idx )

    gens_vh_θh_view =
        @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    """
    gens_ur_ui = get_gens_ur_ui(nodes_pf_U_view, gens_idx )

    gens_ur_ui_view = @view gens_ur_ui[ 1:length(gens_ur_ui ) ]

    """

    #---------------------------------------------
    #---------------------------------------------
    
    idq_wt_pad_view[gens_idx] .=
        [ industrial_model_get_dynamic_idq_θ_π_vhθh(
            vh_θh, δ_ω_ed_eq,
            ra_Xd_dash_Xq_dash )
          for (vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash) in
              zip( gens_vh_θh_view ,
                   gen_nodes_δ_ω_ed_dash_eq_dash_views,
                   ra_Xd_dash_Xq_dash_view[gens_idx] ) ]
  

    #---------------------------------------------

    # gens_nodes_ωs_τm_v_ref =
    #     get_gens_nodes_ωs_τm_v_ref(
    #         gens_nodes_collection,
    #         bus_dict_init )

    # gens_nodes_ωs_τm_v_ref_view =
    #     view( gens_nodes_ωs_τm_v_ref,
    #           1:length( gens_nodes_ωs_τm_v_ref ) )

    #----------------------------------------------

    gens_nodes_ωs_τm_vref_porder =
        get_gens_nodes_ωs_τm_vref_porder(
            gens_nodes_collection, bus_dict_init )

    gens_nodes_ωs_τm_vref_porder_view =
        view(
            gens_nodes_ωs_τm_vref_porder,
            1:length(gens_nodes_ωs_τm_vref_porder) )
    
    #-----------------------------------------------

    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )

    
    #-----------------------------------------------
    #-----------------------------------------------

    gens_dynamic_id_iq_pg_vh_by_vhθh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )

    #-----------------------------------------------  

    nodes_u_Idx_in_ranges =
        get_nodes_u_Idx_in_ranges(
            nodes_u_Idx )    
    #-----------------------------------------------
    
    non_gens_idx =
        get_load_trans_nodes_Idx(
            netd.nodes )
    
    #------------------------------------------------

    industrial_model_nodes_voltage_Idx =
        (; gens_idx,
         non_gens_idx,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges )
    
    industrial_model_misc_Idx =
        (; nodes_u_Idx,
         nodes_u_Idx_in_ranges,
         industrial_model_pure_states_Idx,
         industrial_model_each_gen_nodes_pure_states_idx_in_state,
         gens_nodes_collection )

    para_update_gen_Ax_aux =
        (; industrial_model_pure_states_view_in_state,
         industrial_model_each_gen_nodes_pure_states_idx_in_state,
         industrial_model_pure_states_Idx )

    industrial_model_para_aux_inputs =
        (; industrial_model_nodes_voltage_Idx,
         industrial_model_pure_states_Idx,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         gen_nodes_ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
         industrial_model_pure_states_view_in_state,
         gen_nodes_δ_ω_ed_dash_eq_dash_views,
         gens_vh_θh_view, nodes_pf_U_view )
    

    industrial_model_pf_para =
        (; gens_dynamic_id_iq_pg_vh_by_vhθh,
         gens_nodes_ωs_τm_vref_porder_view )

    
    """ need by their views """

    industrial_model_ωs_τm_vref_vhθh_idq =
        (; gens_nodes_ωs_τm_vref_porder,
         gens_vh_θh, idq_wt_pad )
    

    industrial_model_dyn_pf_up_para =
        (; nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_δ_ed_dash_eq_dash_Idxs_in_state )
    
    #----------------------------------------------

    industrial_model_idq_pf_cal =
        (;  idq_wt_pad_view,
           gens_idx )
    
    #---------------------------------------------

    return (; nodes_cb_sw, state, global_pf_param,
            named_tup_pf_result,
            industrial_model_misc_Idx,
            para_update_gen_Ax_aux,
            industrial_model_para_aux_inputs,
            industrial_model_pf_para,
            industrial_model_ωs_τm_vref_vhθh_idq ,
            industrial_model_dyn_pf_up_para,
            industrial_model_idq_pf_cal,
            gens_nodes_τm_vf  )
    
end


function get_industrial_model_pf_param_views_and_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen  =false )


    dict_sys_to_industry =
        get_net_to_industrial_model_indices_dict(
            netd; no_control_device = only_gen  )

    (;pure_states_Idx_in_system, ur_ui_Idx_in_system,
     industrial_Idx, industrial_model_pure_states_Idx,
     industrial_model_ur_ui_Idx,
     pure_states_and_ur_ui_Idx_in_system,
     industrial_model_state_Idx,
     net_to_industrial_idx_conversion_dict) =
        get_industrial_model_indices_and_conversion_dict(
            netd; no_control_device = only_gen  )

    # -----------------------------------------------
    # -----------------------------------------------  

    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_δ_ω_ed_dash_eq_dash_Idx(
            netd; no_control_device = only_gen )

    nodes_ω_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_ω_ed_dash_eq_dash_Idx(
            netd; no_control_device = only_gen )

    nodes_δ_ed_dash_eq_dash_Idxs_in_state =
        get_industrial_δ_ed_dash_eq_dash_Idx(
            netd; no_control_device = only_gen )

    #-----------------------------------------------

    industrial_model_each_gen_nodes_pure_states_idx_in_state =
        get_industrial_gens_pure_states_indices_in_state(
            netd; no_control_device = only_gen )

    industrial_model_each_gen_nodes_stab_states_idx_in_state =
        get_industrial_gens_stab_states_indices_in_state(
            netd; no_control_device = only_gen )

    #----------------------------------------------

    gen_nodes_ra_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d_dash, :X_q_dash ] )

    gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list =
                [ :ra, :X_d, :X_q,
                  :X_d_dash, :X_q_dash ] )

    gen_nodes_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list =
                [ :D, :H, :ωs, :ra, :xℓ, :X_d,
                  :X_q, :X_d_dash, :X_q_dash,
                  :T_d_dash, :T_q_dash ],
            gens_view_only = only_gen )

    gen_nodes_sub_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list =
                [ :D, :H, :ωs, :ra, :xℓ, :X_d,
                  :X_q, :X_d_dash, :X_q_dash,
                  :X_d_2dash, :X_q_2dash, :T_d_dash,
                  :T_q_dash, :T_d_2dash, :T_q_2dash ],
            gens_view_only = only_gen )

    # -----------------------------------------------
    # -----------------------------------------------

    # pf_net_param = get_powerflow_net_parameters( netd )

    pf_net_param = get_industrial_model_powerflow_net_parameters(
        netd; no_control_device = only_gen )

    (; pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits, pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) = pf_net_param

    (; Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net

    (;slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
        pf_idx_and_state

    (;ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits =
        pf_limits

    (;slack_bus_idx,
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


    #---------------------------------------------    
    #---------------------------------------------

    state = zeros(
        length( generate_industrial_model_sym(
            ; nodes = netd.nodes,
            no_control_device = only_gen ) ) )

    #---------------------------------------------   

    gen_nodes_δ_ω_ed_dash_eq_dash_views =
        get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            state,
            nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    industrial_model_pure_states_view_in_state =
        get_industrial_pure_states_view_in_state(
            state, industrial_model_pure_states_Idx )

    #---------------------------------------------------- 
    #----------------------------------------------------  

    state_view = view(state, 1:length(state))

    # ----------------------------------------------------

    # for pf flat start, ur = vh = 1,  ui = θh = 0

    # remove init_pf = true


    if init_pf == true 

        state_view[ur_idx] .= ones(  length( ur_IDX ))
        state_view[ui_idx] .= zeros( length( ui_IDX ))

    end

    # ----------------------------------------------------

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]

    #----------------------------------------------------

    nodes_u_view  = [
        view(state, nodes_u_Idx[Ind])
        for Ind in
            collect(1:length(nodes_u_Idx)) ]

    # ---------------------------------------------------- 

    nodes_pf_U_view  = [
        view(pf_state , nodes_u_Idx[Ind])
        for Ind in
            collect(1:length(nodes_u_Idx)) ] 

    # ----------------------------------------------------

    uh_state = state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]

    x0_ur_ui = [state_view[ur_ui_idx][ ur_IDX ]...;
                state_view[ur_ui_idx][ ui_IDX ]...]

    x0_vh_θh =
        [abs.(uh_state)...;
         angle.(uh_state)...]

    working_vh_θh_view =
        view(x0_vh_θh, 1:length(x0_vh_θh))

    x0_vh_view  = @view x0_vh_θh[vh_IDX]

    x0_θh_view  = @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view =
        @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0  = [ red_vh_θh_0_view; ]

    mismatch  = similar( red_vh_θh_0 )

    # ----------------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))

    # ----------------------------------------------------

    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet =
        zeros(ComplexF64, length( uh_state ))

    Inet_view  =
        view( Inet, 1:length( Inet ) )

    Iinj =
        zeros(ComplexF64, length( uh_state ))

    Iinj_view  =
        view(Iinj, 1:length( Iinj ))

    idq_wt_pad =
        zeros(ComplexF64, length( uh_state ))

    idq_wt_pad_view =
        view(idq_wt_pad, 1:length( uh_state ) )

    #----------------------------------------------------

    global_pf_views =
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view, Iinj_view )

    sd_pf_views =
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          gen_nodes_δ_ω_ed_dash_eq_dash_views )

    # ---------------------------------------------------
    # global_pf param
    # ---------------------------------------------------

    global_pf_param =
        ( pf_net_param,
          sd_pf_views, mismatch )

    branches_name =
        collect(keys( netd.edges ))

    nodes_name =
        collect(keys( netd.nodes ))    

    #----------------------------------------------------
    #----------------------------------------------------

    named_tup_pf_result =
        power_balance_powerflow(
            x0_vh_θh, mismatch, sd_pf_views,
            (nodes_name, branches_name) ,
            pf_net_param ;
            maxiter=maxiter,
            ftol=ftol ,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )

    #----------------------------------------------------

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init =
        named_tup_pf_result.branch_dict_init

    #---------------------------------------------
    #---------------------------------------------

    state .= industrial_model_init_operationpoint(
        netd, bus_dict_init; pure = :pure, no_control_device = only_gen )

    return (; dict_sys_to_industry, pure_states_Idx_in_system, ur_ui_Idx_in_system, industrial_Idx, industrial_model_pure_states_Idx, industrial_model_ur_ui_Idx, pure_states_and_ur_ui_Idx_in_system, industrial_model_state_Idx, net_to_industrial_idx_conversion_dict, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, nodes_ω_ed_dash_eq_dash_Idxs_in_state, nodes_δ_ed_dash_eq_dash_Idxs_in_state, industrial_model_each_gen_nodes_pure_states_idx_in_state, industrial_model_each_gen_nodes_stab_states_idx_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view, gen_nodes_dyn_param_view, gen_nodes_sub_dyn_param_view, pf_net_param, ra_Xd_dash_Xq_dash_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, industrial_model_pure_states_view_in_state, state_view, pf_state, nodes_u_view, nodes_pf_U_view, x0_vh_θh, working_vh_θh_view, red_vh_θh_0_view, mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, branches_name, nodes_name, named_tup_pf_result, bus_dict_init, branch_dict_init, state, gens_idx, nodes_u_Idx )
    
    
end



function get_industrial_model_pf_param_views_and_init_with_or_no_controller(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen =false  )

    # ----------------------------------------------
    # ----------------------------------------------

    nodes_cb_sw           = get_nodes_cb_sw(netd.nodes)

    gens_nodes            = get_gens_nodes( netd.nodes  )

    non_gens_nodes        = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #----------------------------------------------
    #---------------------------------------------- 

    (;dict_sys_to_industry,
     pure_states_Idx_in_system,
     ur_ui_Idx_in_system,
     industrial_Idx,
     industrial_model_pure_states_Idx,
     industrial_model_ur_ui_Idx,
     pure_states_and_ur_ui_Idx_in_system,
     industrial_model_state_Idx,
     net_to_industrial_idx_conversion_dict,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_δ_ed_dash_eq_dash_Idxs_in_state,
     industrial_model_each_gen_nodes_pure_states_idx_in_state,
     industrial_model_each_gen_nodes_stab_states_idx_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     gen_nodes_dyn_param_view,
     gen_nodes_sub_dyn_param_view,
     pf_net_param,
     ra_Xd_dash_Xq_dash_view,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     industrial_model_pure_states_view_in_state,
     state_view, pf_state, nodes_u_view,
     nodes_pf_U_view,
     x0_vh_θh,
     working_vh_θh_view,
     red_vh_θh_0_view,
     mismatch,
     Jac_vh_θh,
     Inet,
     Inet_view,
     Iinj,
     Iinj_view,
     idq_wt_pad,
     idq_wt_pad_view,
     global_pf_views,
     sd_pf_views,
     global_pf_param,
     branches_name,
     nodes_name,
     named_tup_pf_result,
     bus_dict_init,
     branch_dict_init,
     state,
     gens_idx,
     nodes_u_Idx) =
        get_industrial_model_pf_param_views_and_init(
            netd;
            maxiter=40,
            ftol=1000*eps(),
            xtol=1000*eps(),
            init_pf = true,
            with_δ_ed_eq = false,
            only_gen  =false)
    
    #---------------------------------------------------- 
    #---------------------------------------------------- 
    
    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )


    #----------------------------------------------------

    nodes_cb_sw =
        get_nodes_cb_sw(netd.nodes)

    #----------------------------------------------------

    gens_vh_θh = get_gens_vh_θh(
        nodes_pf_U_view, gens_idx )

    gens_vh_θh_view =
        @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    #---------------------------------------------------
    #---------------------------------------------------

    # gen_uh  = (named_tup_pf_result.Vbus)[ gen_idx ]

    idq_wt_pad_view[gens_idx] .=  [
        industrial_model_get_dynamic_idq_θ_π_vhθh(
            vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
        for ( vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in
            zip( gens_vh_θh_view ,
                 gen_nodes_δ_ω_ed_dash_eq_dash_views,
                 ra_Xd_dash_Xq_dash_view[gens_idx] ) ]

    #---------------------------------------------------

    gens_nodes_ωs_τm_vref_porder =
        get_gens_nodes_ωs_τm_vref_porder(
            gens_nodes_collection,
            bus_dict_init )

    gens_nodes_ωs_τm_vref_porder_view =
        view( gens_nodes_ωs_τm_vref_porder,
              1:length(gens_nodes_ωs_τm_vref_porder) )

    #---------------------------------------------------

    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )

    vec_τm_vf_views =
        view( gens_nodes_τm_vf,
              1:length(gens_nodes_τm_vf) )
    
    #---------------------------------------------------   
    #---------------------------------------------------

    gens_dynamic_id_iq_pg_vh_by_vhθh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )

    gens_dynamic_id_iq_pg_vh_by_vhθh_view =
        view( gens_dynamic_id_iq_pg_vh_by_vhθh,
              1:length(gens_dynamic_id_iq_pg_vh_by_vhθh) )

    """
    gens_dynamic_id_iq_pg_vh_by_ur_ui = get_gens_dynamic_id_iq_pg_vh_by_ur_ui( gens_ur_ui_post_pf, gen_nodes_δ_ω_ed_dash_eq_dash_views, gen_nodes_ra_Xd_dash_Xq_dash_view )

    """
    #---------------------------------------------------  

    nodes_u_Idx_in_ranges =
        get_nodes_u_Idx_in_ranges(
            nodes_u_Idx )
    
    #---------------------------------------------------

    non_gens_idx =
        get_load_trans_nodes_Idx(
            netd.nodes )

    #---------------------------------------------------

    industrial_model_nodes_voltage_Idx =
        (; gens_idx,
         non_gens_idx,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges )

    industrial_model_misc_Idx =
        (; nodes_u_Idx,
         nodes_u_Idx_in_ranges,
         industrial_model_pure_states_Idx,
         industrial_model_each_gen_nodes_pure_states_idx_in_state,
         gens_nodes_collection )

    para_update_gen_Ax_aux =
        (; industrial_model_pure_states_view_in_state,
         industrial_model_each_gen_nodes_pure_states_idx_in_state,
         industrial_model_pure_states_Idx )

    industrial_model_para_aux_inputs =
        (; industrial_model_nodes_voltage_Idx,
         industrial_model_pure_states_Idx,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         gen_nodes_ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
         industrial_model_pure_states_view_in_state,
         gen_nodes_δ_ω_ed_dash_eq_dash_views,
         gens_vh_θh_view,
         nodes_pf_U_view )

    if only_gen == false
        
        industrial_model_pf_para =
            (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
             gens_nodes_ωs_τm_vref_porder_view )
        
    else
        industrial_model_pf_para =
            (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
             gens_nodes_ωs_τm_vref_porder_view,
             vec_τm_vf_views )
        
    end


    """ need by their views """

    industrial_model_ωs_τm_vref_vhθh_idq =
        (; gens_dynamic_id_iq_pg_vh_by_vhθh,
         gens_nodes_ωs_τm_vref_porder,
         gens_nodes_τm_vf,
         gens_vh_θh,
         idq_wt_pad )


    industrial_model_dyn_pf_up_para =
        (; nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_δ_ed_dash_eq_dash_Idxs_in_state )

    #---------------------------------------------------

    industrial_model_idq_pf_cal =
        (;  idq_wt_pad_view, gens_idx )

    #---------------------------------------------------

    return (; nodes_cb_sw,
            state,
            global_pf_param,
            named_tup_pf_result,
            industrial_model_misc_Idx,
            para_update_gen_Ax_aux,
            industrial_model_para_aux_inputs,
            industrial_model_pf_para,
            industrial_model_ωs_τm_vref_vhθh_idq,
            industrial_model_dyn_pf_up_para,
            industrial_model_idq_pf_cal )

    #---------------------------------------------------
    #---------------------------------------------------

end


#---------------------------------------------------
# v2
#---------------------------------------------------


function v2_get_im_model_pf_param_views_and_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen  =false)

    #-------------------------------
    
    # dynamics_case = case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    # dynamics_case =
    #     case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    
    # netd     = NetworkData( dynamics_case()... )
    # only_gen = true
    # maxiter  = 40
    # ftol     = 1000*eps()
    # xtol     = 1000*eps()
    # init_pf  = true
    # with_δ_ed_eq = false
    
    # #-------------------------------
    

    dict_sys_to_im = get_net_to_im_indices_dict( netd  )

    (; im_vars_indices_in_system,
     pure_states_Idx_in_system,
     im_algebraic_vars_Idx_in_system,
     ur_ui_Idx_in_system,
     im_vars_and_ur_ui_Idx_in_system,
     im_vars_Idx_in_state,
     nodes_ur_ui_Idx_in_state,
     im_state_Idx,
     each_gens_im_vars_Idx_in_state,
     net_to_im_idx_conversion_dict)  =
         get_im_indices_and_conversion_dict(
             netd  )
    
    # ---------------------------------------------------
    # ---------------------------------------------------    

    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
        get_im_δ_ω_ed_dash_eq_dash_Idx(netd )

    nodes_ω_ed_dash_eq_dash_Idxs_in_state =
        get_im_ω_ed_dash_eq_dash_Idx(netd )

    nodes_δ_ed_dash_eq_dash_Idxs_in_state =
        get_im_δ_ed_dash_eq_dash_Idx(netd )

    #---------------------------------------------------
    #---------------------------------------------------

    gen_nodes_ra_Xd_dash_Xq_dash_view =
        get_im_model_gens_params_view_in_param_values(
        netd.nodes_param,
        netd.nodes; param_list = [
            :ra, :X_d_dash, :X_q_dash ] )

    gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_im_model_gens_params_view_in_param_values(
        netd.nodes_param,
        netd.nodes; param_list = [
            :ra, :X_d, :X_q, :X_d_dash, :X_q_dash ] )

    gen_nodes_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes
            ; param_list = [
                :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q,
                :X_d_dash, :X_q_dash,
                :T_d_dash, :T_q_dash ])

    gen_nodes_sub_dyn_param_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes
            ; param_list = [
                :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q,
                :X_d_dash, :X_q_dash,
                :X_d_2dash, :X_q_2dash,
                :T_d_dash, :T_q_dash,
                :T_d_2dash, :T_q_2dash ] )

    # -----------------------------------------------------
    # -----------------------------------------------------

    pf_net_param =
        get_im_model_powerflow_net_parameters( netd  )

    (;pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits, pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) =
         pf_net_param

    (;Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
        pf_net

    (;slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
        pf_idx_and_state

    (; ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
        pf_param_views

    load_trans_nodes_Idx_and_vlimits =
        pf_limits

    (;slack_bus_idx,
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


    #----------------------------------------------------    
    #----------------------------------------------------

    state = zeros(length(
        generate_im_sym(
            ; nodes = netd.nodes )  ) )

    #----------------------------------------------------   

    gen_nodes_δ_ω_ed_dash_eq_dash_views =
        get_im_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            state,
            nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    # im_model_pure_states_view_in_state =
    #     get_im_pure_states_view_in_state(
    #         state,
    #         pure_states_Idx )

    im_vars_view_in_state =
        get_im_vars_view_in_state( state, im_vars_Idx_in_state )

    #----------------------------------------------------  

    state_view = view(state, 1:length(state))

    # ----------------------------------------------------

    # for pf flat start, ur = vh = 1,  ui = θh = 0

    if init_pf == true 

        state_view[ur_idx] .= ones(  length( ur_IDX ))
        state_view[ui_idx] .= zeros( length( ui_IDX ))

    end

    # ----------------------------------------------------

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]

    #----------------------------------------------------

    nodes_u_view  = [
        view(state, nodes_u_Idx[Ind])
        for Ind in collect(1:length(nodes_u_Idx)) ]

    # ---------------------------------------------------- 

    nodes_pf_U_view  = [
        view(pf_state , nodes_u_Idx[Ind])
        for Ind in collect(1:length(nodes_u_Idx)) ] 

    # ----------------------------------------------------

    uh_state = state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]

    x0_ur_ui = [
        state_view[ur_ui_idx][ ur_IDX ]...;
        state_view[ur_ui_idx][ ui_IDX ]...]

    x0_vh_θh      = [abs.(uh_state)...; angle.(uh_state)...]

    working_vh_θh_view = view(x0_vh_θh, 1:length(x0_vh_θh))

    x0_vh_view       = @view x0_vh_θh[vh_IDX]

    x0_θh_view       = @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0      = [ red_vh_θh_0_view; ]

    mismatch         = similar( red_vh_θh_0 )

    # ----------------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))

    # ----------------------------------------------------

    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet            = zeros(ComplexF64, length( uh_state ))

    Inet_view       =  view( Inet, 1:length( Inet ) )

    Iinj            = zeros(ComplexF64, length( uh_state ))

    Iinj_view       =  view(Iinj, 1:length( Iinj ))

    idq_wt_pad      = zeros(ComplexF64, length( uh_state ))

    idq_wt_pad_view = view(idq_wt_pad,1:length(uh_state) )


    # idq_wt_pad =
    #     [[0.0, 0.0]
    #      for idx in 1:length(uh_state )]

    # idq_wt_pad_view =
    #     view(idq_wt_pad, 1:length( uh_state ) )
    
    #----------------------------------------------------

    global_pf_views = (
        working_vh_θh_view,
        nodes_pf_U_view,
        Inet_view,
        Iinj_view )

    sd_pf_views = (
        working_vh_θh_view,
        nodes_pf_U_view,
        Inet_view,
        Iinj_view,
        gen_nodes_δ_ω_ed_dash_eq_dash_views )

    # ---------------------------------------------------
    # global_pf param
    # ---------------------------------------------------

    global_pf_param = (
        pf_net_param,
        sd_pf_views,
        mismatch )

    branches_name  = collect(keys( netd.edges ))

    nodes_name     = collect(keys( netd.nodes ))    

    #----------------------------------------------------

    named_tup_pf_result = power_balance_powerflow(
        x0_vh_θh,
        mismatch,
        sd_pf_views,
        (nodes_name, branches_name) ,
        pf_net_param;
        maxiter=maxiter,
        ftol=ftol,
        xtol=xtol,
        with_δ_ed_eq = with_δ_ed_eq )

    #----------------------------------------------------

    bus_dict_init = named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init

    #---------------------------------------------------

    state .= im_model_init_operationpoint(
        netd, bus_dict_init  )

    return (
        ; im_vars_indices_in_system,
        pure_states_Idx_in_system,
        im_algebraic_vars_Idx_in_system,
        ur_ui_Idx_in_system,
        im_vars_and_ur_ui_Idx_in_system,
        im_vars_Idx_in_state,
        nodes_ur_ui_Idx_in_state,
        im_state_Idx,
        each_gens_im_vars_Idx_in_state,
        net_to_im_idx_conversion_dict,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_δ_ed_dash_eq_dash_Idxs_in_state,
        gen_nodes_ra_Xd_dash_Xq_dash_view,
        gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
        gen_nodes_dyn_param_view,
        gen_nodes_sub_dyn_param_view,
        pf_net_param,
        ra_Xd_dash_Xq_dash_view,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        im_vars_view_in_state,
        state_view,
        pf_state,
        nodes_u_view,
        nodes_pf_U_view,
        x0_vh_θh,
        working_vh_θh_view,
        red_vh_θh_0_view,
        mismatch,
        Jac_vh_θh,
        Inet,
        Inet_view,
        Iinj,
        Iinj_view,
        idq_wt_pad,
        idq_wt_pad_view,
        global_pf_views,
        sd_pf_views,
        global_pf_param,
        branches_name,
        nodes_name,
        named_tup_pf_result,
        bus_dict_init,
        branch_dict_init,
        state,
        gens_idx,
        nodes_u_Idx )
    
    
end



function v2_get_im_sys_pf_param_views_and_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen =false  )

    # ---------------------------------------------------

    nodes_cb_sw           = get_nodes_cb_sw(netd.nodes)

    gens_nodes            = get_gens_nodes( netd.nodes  )

    non_gens_nodes        = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #---------------------------------------------------

    (; im_vars_indices_in_system,
     pure_states_Idx_in_system,
     im_algebraic_vars_Idx_in_system,
     ur_ui_Idx_in_system,
     im_vars_and_ur_ui_Idx_in_system,
     im_vars_Idx_in_state,
     nodes_ur_ui_Idx_in_state,
     im_state_Idx,
     each_gens_im_vars_Idx_in_state,
     net_to_im_idx_conversion_dict,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_δ_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     gen_nodes_dyn_param_view,
     gen_nodes_sub_dyn_param_view,
     pf_net_param,
     ra_Xd_dash_Xq_dash_view,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     im_vars_view_in_state,
     state_view,
     pf_state,
     nodes_u_view,
     nodes_pf_U_view,
     x0_vh_θh,
     working_vh_θh_view,
     red_vh_θh_0_view,
     mismatch,
     Jac_vh_θh,
     Inet,
     Inet_view,
     Iinj, Iinj_view,
     idq_wt_pad,
     idq_wt_pad_view,
     global_pf_views,
     sd_pf_views,
     global_pf_param,
     branches_name,
     nodes_name,
     named_tup_pf_result,
     bus_dict_init,
     branch_dict_init,
     state,
     gens_idx,
     nodes_u_Idx)  =
        v2_get_im_model_pf_param_views_and_init(
            netd;
            maxiter=40,
            ftol=1000*eps() ,
            xtol=1000*eps(),
            init_pf = true,
            with_δ_ed_eq = false,
            only_gen = only_gen )
    
    #---------------------------------------------------- 
    
    update_im_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )


    #----------------------------------------------------

    nodes_cb_sw =
        get_nodes_cb_sw(netd.nodes)

    #----------------------------------------------------

    gens_vh_θh =
        get_gens_vh_θh( nodes_pf_U_view, gens_idx )

    gens_vh_θh_view =
        @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    #---------------------------------------------------
    
    
    dyn_idq = [ get_pf_dynamic_idq_θ_π_vhθh(
            vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
                for ( vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in
                    zip(
                        gens_vh_θh_view ,
                        gen_nodes_δ_ω_ed_dash_eq_dash_views,
                        ra_Xd_dash_Xq_dash_view[gens_idx] ) ]

    idq_wt_pad_view[gens_idx] .= dyn_idq
    
    # for idx in gens_idx
        
    #     idq_wt_pad_view[idx] = dyn_idq[idx]
        
    # end
    
    #---------------------------------------------------

    gens_nodes_ωs_ωref0_vref0_porder0 =
        get_gens_nodes_ωs_ωref0_vref0_porder0(
        gens_nodes_collection,
        bus_dict_init )

    gens_nodes_ωs_ωref0_vref0_porder0_view = view(
        gens_nodes_ωs_ωref0_vref0_porder0,
        1:length( gens_nodes_ωs_ωref0_vref0_porder0 ) )

    #---------------------------------------------------

    gens_nodes_τm_vf =
        get_gens_τm_vf( gens_nodes_collection,
                        bus_dict_init )

    vec_τm_vf_views =
        view( gens_nodes_τm_vf,
              1:length(gens_nodes_τm_vf) )
    
    #---------------------------------------------------   

    gens_dynamic_id_iq_pg_vh_by_vhθh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
        gens_vh_θh_view,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        gen_nodes_ra_Xd_dash_Xq_dash_view )

    gens_dynamic_id_iq_pg_vh_by_vhθh_view =
        view( gens_dynamic_id_iq_pg_vh_by_vhθh,
            1:length(gens_dynamic_id_iq_pg_vh_by_vhθh) )
    
    #---------------------------------------------------  

    nodes_u_Idx_in_ranges = get_nodes_u_Idx_in_ranges(
        nodes_u_Idx )
    
    #---------------------------------------------------

    non_gens_idx = get_load_trans_nodes_Idx( netd.nodes )

    #---------------------------------------------------

    
    im_nodes_voltage_Idx = (
        ; gens_idx,
        non_gens_idx,
        nodes_u_Idx,
        nodes_u_Idx_in_ranges )

    
    im_misc_Idx = (
        ; nodes_u_Idx,
        nodes_u_Idx_in_ranges,
        im_vars_Idx_in_state, 
        each_gens_im_vars_Idx_in_state, 
        gens_nodes_collection )

    
    para_update_gen_Ax_aux = (
        ; im_vars_view_in_state, 
        each_gens_im_vars_Idx_in_state, 
        im_vars_Idx_in_state ) 

    
    im_para_aux_inputs = (
        ; im_nodes_voltage_Idx,
        im_vars_Idx_in_state, 
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
        gen_nodes_ra_Xd_dash_Xq_dash_view,
        gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
        im_vars_view_in_state,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        gens_vh_θh_view,
        nodes_pf_U_view )

    
    im_pf_para = (
        ; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
        gens_nodes_ωs_ωref0_vref0_porder0_view  )


    """ need by their views """

    im_ωs_ωref_vref_vhθh_idq_etc =  (
        ; gens_dynamic_id_iq_pg_vh_by_vhθh,
        gens_nodes_ωs_ωref0_vref0_porder0,
        gens_nodes_τm_vf,
        gens_vh_θh,
        idq_wt_pad )

    im_dyn_pf_up_para = (
        ; nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_δ_ed_dash_eq_dash_Idxs_in_state )

    #---------------------------------------------------

    im_idq_pf_cal = (
        ; idq_wt_pad_view,
        gens_idx )

    #---------------------------------------------------

    return (
        ;nodes_cb_sw,
        state,
        global_pf_param,
        named_tup_pf_result,
        im_misc_Idx, 
        para_update_gen_Ax_aux, 
        im_para_aux_inputs, 
        im_pf_para, 
        im_ωs_ωref_vref_vhθh_idq_etc, 
        im_dyn_pf_up_para, 
        im_idq_pf_cal ) 
    
end








#-----------------------------------------------------
#-----------------------------------------------------
# Simulations
#-----------------------------------------------------
#-----------------------------------------------------


#-----------------------------------------------------
#-----------------------------------------------------
# Simulation functions
#-----------------------------------------------------
#-----------------------------------------------------


function simulate_a_dynamics_industrial_model(
    ; dynamics_case = dynamics_case,
    only_gen = false,
    sim_timespan  = sim_timespan,
    ode_alg       = ode_alg,
    dyn_global_pf_options = dyn_global_pf_options,
    sta_global_pf_options = sta_global_pf_options )

    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    #-----------------------------------------------------
    #----------------------------------------------------- 
    
    # Dyn_Nodes, Dyn_Branches = dynamics_case()

    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------
    #-----------------------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #-----------------------------------------------------

    # (network_bus_names,
    #  non_gens_bus_names,
    #  gens_bus_names,
    #  Loads_bus_names,
    #  Trans_bus_names) =
    #      make_case_buses_names(
    #          ; case_fun = dynamics_case )

    # (net_bus_volts_labels,
    #  gens_nodes_pure_states_labels,
    #  gens_nodes_stab_states_label,
    #  gens_nodes_algebraic_and_states_labels) =
    #     generate_industrial_model_labels(
    #         network_bus_names,
    #         gens_nodes_collection )

    # (industrial_model_sym,
    #  industrial_model_mass_matrix) =
    #     generate_industrial_model_sym_and_mass_matrix(
    #         gens_nodes_pure_states_labels,
    #         net_bus_volts_labels)
    
    #-----------------------------------------------------
    
    (network_bus_names,
     non_gens_bus_names,
     gens_bus_names,
     Loads_bus_names,
     Trans_bus_names) =
         make_case_buses_names( netd.nodes )

    net_class_names =
        (;
         network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-----------------------------------------------------

    (net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_stab_states_label,
     gens_nodes_algebraic_and_states_labels) =
        generate_industrial_model_labels(
            ; nodes =  netd.nodes )

    net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_stab_states_label,
         gens_nodes_algebraic_and_states_labels )
    
    #-----------------------------------------------------

    (industrial_model_sym,
     industrial_model_mass_matrix) =
        generate_industrial_model_sym_and_mass_matrix(
            ; nodes = netd.nodes )

    #-----------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         net_states_and_var_labels,
         industrial_model_sym )
    
    #-----------------------------------------------------
    
    industrial_model_pf_sys_param_sys_views_sys_industrial_model_init =
        get_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init(
            netd; sta_global_pf_options... )

    #-----------------------------------------------------
    #-----------------------------------------------------

    (nodes_cb_sw, state, global_pf_param,
     named_tup_pf_result,
     industrial_model_misc_Idx,
     para_update_gen_Ax_aux,
     industrial_model_para_aux_inputs,
     industrial_model_pf_para,
     industrial_model_ωs_τm_vref_vhθh_idq,
     industrial_model_dyn_pf_up_para,
     industrial_model_idq_pf_cal,
     gens_nodes_τm_vf)  =
         industrial_model_pf_sys_param_sys_views_sys_industrial_model_init

    
    # -----------------------------------------------------

    pf_net_param, sd_pf_views, _ = global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ = pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views = sd_pf_views

    
    #-----------------------------------------------------
    # Stability
    #-----------------------------------------------------  
    #-----------------------------------------------------


    if only_gen == false
        Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs =
            create_gens_nodes_aggregate_system_matrices_Idx_and_views(  gens_nodes_collection; only_gen = only_gen )

        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views)  = Ax_Bx_Cx_views
        
        (Ax_matrix,
         Bx_matrix,
         Cx_matrix) = Ax_Bx_Cx_matrix

        (nodes_state_Idx,
         Bx_idxs,
         Cx_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen =
                only_gen,
            vec_Ax_τm_vf_views =
                nothing )
        
    else
        (Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
            create_gens_nodes_aggregate_system_matrices_Idx_and_views(  gens_nodes_collection; only_gen = only_gen )
        
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         vec_Ax_τm_vf_views) =
             Ax_Bx_Cx_views
        
        (Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix) =
            Ax_Bx_Cx_matrix
        
        (nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         τm_vf_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen =
                only_gen,
            vec_Ax_τm_vf_views =
                vec_Ax_τm_vf_views )
        
    end
    
    #-----------------------------------------------------
    #-----------------------------------------------------

    (gens_dynamic_id_iq_pg_vh_by_vhθh,
     gens_nodes_ωs_τm_vref_porder_view) =
         industrial_model_pf_para

    
    #-----------------------------------------------------

    ode_fun_para =
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         industrial_model_pf_para  )
    
    para_model_fun = (vec_Ax_views, ode_fun_para )

    #-----------------------------------------------------

    sim_state_x0  = state
    
    sim_timespan  = sim_timespan

    #-----------------------------------------------------
    
    dyn_global_pf_options = dyn_global_pf_options
    
    counter_array  = [1]
    
    stateDiffCache = similar( sim_state_x0 )
    
    sim_fun_para   =
        (; netd,
         nodes_cb_sw,
         global_pf_param,
         counter_array,
         state,
         stateDiffCache,
         dyn_global_pf_options,
         para_update_gen_Ax_aux,
         para_model_fun,
         industrial_model_misc_Idx,
         industrial_model_para_aux_inputs,
         industrial_model_dyn_pf_up_para,
         industrial_model_idq_pf_cal  )
    
    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------
    #-----------------------------------------------------

    sim_func  =
        dynamics_industrial_model!
    
    #---------------------------------------------------    
        
    sim_ode_func! = ODEFunction{true}(
        sim_func;
        mass_matrix =
            industrial_model_mass_matrix,
        syms =
            industrial_model_sym )
    
    sim_prob = ODEProblem(
        sim_ode_func!,
        sim_state_x0,
        sim_timespan,
        sim_fun_para )

    sim_sol =
        DifferentialEquations.solve(sim_prob, ode_alg )

    return (; sim_sol, para_net_names_labels_syms) 
      
end



function simulate_a_dynamics_industrial_model_with_or_no_controllers(
    ; dynamics_case = nothing,
    only_gen = nothing,
    sim_timespan  = nothing,
    ode_alg = ode_alg,
    dyn_global_pf_options = nothing ,
    sta_global_pf_options = nothing,
    with_cb = false)

    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    netd  = NetworkData( dynamics_case()... )
    
    #---------------------------------------------------    

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #---------------------------------------------------
    
    (network_bus_names, non_gens_bus_names,
     gens_bus_names, Loads_bus_names,
     Trans_bus_names) =
        make_case_buses_names(  netd.nodes )
    
    net_class_names =
        (; network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    (net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_stab_states_label,
     gens_nodes_algebraic_and_states_labels) =
        generate_industrial_model_labels(
            ; nodes =
                netd.nodes,
            no_control_device =
                only_gen )

    net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_stab_states_label,
         gens_nodes_algebraic_and_states_labels )

    #-----------------------------------------------------

    (industrial_model_sym,
     industrial_model_mass_matrix) =
        generate_industrial_model_sym_and_mass_matrix(
            ; nodes =
                netd.nodes,
            no_control_device =
                only_gen )

    #-----------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         net_states_and_var_labels,
         industrial_model_sym )

    #-----------------------------------------------------

    industrial_model_pf_sys_param_sys_views_sys_industrial_model_init =
        get_industrial_model_pf_param_views_and_init_with_or_no_controller(
            netd; sta_global_pf_options...,
            only_gen=only_gen )

    #-----------------------------------------------------

    (nodes_cb_sw, state,
     global_pf_param,
     named_tup_pf_result,
     industrial_model_misc_Idx,
     para_update_gen_Ax_aux,
     industrial_model_para_aux_inputs,
     industrial_model_pf_para,
     industrial_model_ωs_τm_vref_vhθh_idq ,
     industrial_model_dyn_pf_up_para,
     industrial_model_idq_pf_cal)  =
         industrial_model_pf_sys_param_sys_views_sys_industrial_model_init


    # -----------------------------------------------------

    pf_net_param, sd_pf_views, _ = global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ = pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views = sd_pf_views

    #-----------------------------------------------------
    #-----------------------------------------------------


    if only_gen == false
        
        #-------------------------------------------------
        # Stability
        #------------------------------------------------- 

        (;Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
            create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        vec_Ax_views, vec_Bx_views, vec_Cx_views  =
            Ax_Bx_Cx_views
        
        Ax_matrix, Bx_matrix, Cx_matrix =
            Ax_Bx_Cx_matrix

        nodes_state_Idx, Bx_idxs, Cx_idxs =
            idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views, vec_Bx_views, vec_Cx_views),
            gens_nodes_collection; only_gen =
                only_gen, vec_Ax_τm_vf_views = nothing )

        #-------------------------------------------------
        #-------------------------------------------------

        ode_fun_para =
            (; Ax_matrix, Bx_matrix,
             Cx_matrix,
             industrial_model_pf_para  )

        para_model_fun =
            (; vec_Ax_views, ode_fun_para )

    else
        
        #------------------------------------------------
        # Stability
        #-------------------------------------------------  
        #-------------------------------------------------
        
        (;Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
            create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        (;vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         vec_Ax_τm_vf_views) =
             Ax_Bx_Cx_views
        
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix) =
             Ax_Bx_Cx_matrix
        
        (;nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         τm_vf_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen =
                only_gen,
            vec_Ax_τm_vf_views =
                vec_Ax_τm_vf_views )

        #-------------------------------------------------

        # update_gens_τm_vf!(
        #     vec_τm_vf_views, gens_nodes_τm_vf )

        #-------------------------------------------------

        ode_fun_para =
            (; Ax_matrix, Bx_matrix,
             Cx_matrix, Ax_τm_vf_matrix,
             industrial_model_pf_para  )

        para_model_fun =
            (; vec_Ax_views,
             ode_fun_para )

    end
    
    #-----------------------------------------------------

    sim_state_x0  = state

    sim_timespan  = sim_timespan

    #-----------------------------------------------------

    dyn_global_pf_options = dyn_global_pf_options

    counter_array  = [1]

    stateDiffCache = similar( sim_state_x0 )

    sim_fun_para = (
        ; netd,
        nodes_cb_sw,
        global_pf_param,
        counter_array,
        state,
        stateDiffCache,
        dyn_global_pf_options,
        para_update_gen_Ax_aux,
        para_model_fun,
        industrial_model_misc_Idx,
        industrial_model_para_aux_inputs,        
        industrial_model_dyn_pf_up_para,
        industrial_model_idq_pf_cal,                
        only_gen )

    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------

    sim_func =
        dynamics_industrial_model_with_or_no_controller!

    # dynamics_industrial_model!

    #---------------------------------------------------    

    sim_ode_func! = ODEFunction{true}(
        sim_func;
        mass_matrix =
            industrial_model_mass_matrix,
        syms =
            industrial_model_sym )

    sim_prob = ODEProblem(sim_ode_func!,
                          sim_state_x0,
                          sim_timespan,
                          sim_fun_para )

    if with_cb == true
        
        cb = industrial_model_fun_make_state_callbacks(
            collect( values( netd.nodes)) )

        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg, callback = cb )
        
    else
        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg  )

    end
    
    return (;sim_sol, para_net_names_labels_syms) 
      
end



function v2_simulate_a_dynamics_im_model(
    ; dynamics_case = dynamics_case,
    only_gen = false,
    sim_timespan  = sim_timespan,
    ode_alg       = ode_alg,
    dyn_global_pf_options = dyn_global_pf_options,
    sta_global_pf_options = sta_global_pf_options,
    with_cb = false
        )

    # ------------------------------------------------
    # ------------------------------------------------
    # Diagnosis
    # ------------------------------------------------ 
    # -------------------------------------------------

     # dynamics_case =
     #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    #--------------------------------------------------
    #-------------------------------------------------- 
    
    # Dyn_Nodes, Dyn_Branches = dynamics_case()

    netd  = NetworkData( dynamics_case()... )
    
    #-------------------------------------------------
    #-------------------------------------------------

    non_gens_nodes = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #-------------------------------------------------
    #-------------------------------------------------
    
    (;network_bus_names,
     non_gens_bus_names,
     gens_bus_names,
     Loads_bus_names,
     Trans_bus_names) =
         make_case_buses_names(
             netd.nodes )
    
    net_class_names =
        (; network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-------------------------------------------------

    (;net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
        generate_im_model_labels(
            ; nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #--------------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes = netd.nodes )

    #-------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    
    #-------------------------------------------------
    
    im_sys_pf_param_views_and_init =
        v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )

    #-------------------------------------------------
    #-------------------------------------------------

    (; nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para,
     im_idq_pf_cal)  =
        im_sys_pf_param_views_and_init
    
    # -------------------------------------------------

    pf_net_param, sd_pf_views, _ =
        global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx =
        pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
        pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
        sd_pf_views
    
    #--------------------------------------------------
    # Stability
    #--------------------------------------------------  

    (Ax_Bx_Cx_views,
     Ax_Bx_Cx_matrix,
     idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix =
        Ax_Bx_Cx_matrix

    (nodes_state_Idx,
     Bx_idxs, Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs) =
        idxs

    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )
    
    #-------------------------------------------------
    #-------------------------------------------------

    ode_fun_para = (
        ; Ax_matrix,
        Bx_matrix,
        Cx_matrix,
        im_pf_para )
    
    para_model_fun = (
        ; vec_Ax_views,
        ode_fun_para )

    #-------------------------------------------------

    sim_state_x0  = state
    
    sim_timespan  = sim_timespan

    #-------------------------------------------------
    
    dyn_global_pf_options = dyn_global_pf_options
    
    counter_array  = [1]
    
    stateDiffCache = similar( sim_state_x0 )
    
    sim_fun_para = (
        ; netd,
        nodes_cb_sw,
        global_pf_param,
        counter_array,
        state,
        stateDiffCache,
        dyn_global_pf_options,
        para_update_gen_Ax_aux,
        para_model_fun,
        im_misc_Idx,
        im_para_aux_inputs,
        im_dyn_pf_up_para,
        im_idq_pf_cal )
    
    #-----------------------------------------------
    # dynamics simulation
    #-----------------------------------------------
    #-----------------------------------------------

    sim_func = dynamics_im_model!
    
    #-----------------------------------------------
    
    #-----------------------------------------------
        
    sim_ode_func! = ODEFunction{true}(
        sim_func; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
    sim_prob = ODEProblem(sim_ode_func!,
                          sim_state_x0,
                          sim_timespan,
                          sim_fun_para )


    if with_cb == true
        
        cb = industrial_model_fun_make_state_callbacks(
            collect( values( netd.nodes)) )

        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg, callback = cb )
        
    else
        
        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg )
    end
    
    #-------------------------------------------------
    #-------------------------------------------------
    
    return sim_sol, para_net_names_labels_syms 
      
end



function t_v2_simulate_one_plant_dynamics_im_model(
    ; dynamics_case =
        dynamics_case,
    only_gen =
        false,
    sim_timespan =
        sim_timespan,
    ode_alg =
        ode_alg,
    dyn_global_pf_options =
        dyn_global_pf_options,
    sta_global_pf_options =
        sta_global_pf_options,
    with_cb =
        false
        )

    # ------------------------------------------------
    # Diagnosis
    # ------------------------------------------------ 

    # ImplicitMidpoint(), ImplicitMidpoint(autodiff=false),
    
    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)

    sim_model_type        = "im-model"
    
    # ode_alg             = Rodas4()
    # algr_name           = "rodas4" 
    
    ode_alg               = ImplicitMidpoint()    
    algr_name             = "ImplicitMidpoint"

    dt                    = 0.01

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = false

    with_cb  = false

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
        
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #-----------------------------------------------------
    
    ; (network_bus_names,
       non_gens_bus_names,
       gens_bus_names,
       Loads_bus_names,
       Trans_bus_names) =
           make_case_buses_names(
               netd.nodes )
    
    net_class_names =
        (;network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-------------------------------------------------

    (; net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #-------------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes = netd.nodes )

    #-------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    
    #-------------------------------------------------

    (; nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para, im_idq_pf_cal)  =
         v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )
    
    #-------------------------------------------------
    
    (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para

    gens_vh_θh_view =
        im_para_aux_inputs.gens_vh_θh_view
    
    # --------------------------------------------------

    pf_net_param, sd_pf_views, _ =
        global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx =
        pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
        pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
        sd_pf_views

    
    #-------------------------------------------------
    # Stability
    #-------------------------------------------------  
    
    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, states_and_mat_Idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    (; Ax_matrix,
     Bx_matrix,
     Cx_matrix) = Ax_Bx_Cx_matrix

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs)  =
         states_and_mat_Idxs
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )

    #-------------------------------------------------
    
    node_i = 1
    
    gen_node_i_idx =
        nodes_state_Idx[ node_i ]
    
    #--------------------------------------------------

    Ax_m = vec_Ax_views[ node_i ]

    Bx_m = vec_Bx_views[ node_i ]

    Cx_m = vec_Cx_views[ node_i ]
    
    #--------------------------------------------------
    
    plant_i =
        gens_nodes_collection[ node_i ]      

    #--------------------------------------------------

    plant_i_eigvalues =
        eigvals( Matrix( Ax_m ))

    #--------------------------------------------------

    ra_Xd_dash_Xq_dash_i =
        ra_Xd_dash_Xq_dash_view[ node_i ]

    #---------------------------------------------------
    
    id_iq_pg_vh_i =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view[ node_i ]

    ωs_ωref0_vref0_porder0_i =
        gens_nodes_ωs_ωref0_vref0_porder0_view[ node_i ]

    gens_vh_θh_i = [ gens_vh_θh_view[ node_i ]...;]

    #-------------------------------------------------

    dims_gens_vh_θh_ωs_ωref0_vref0_porder0_i =
        length.([ gens_vh_θh_i, ωs_ωref0_vref0_porder0_i ])

    _,_, vh_θh_ωs_ωref0_vref0_porder0_i_Idx =
        create_size_offset_Idx(
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0_i,;
            counter = 0)

    vh_θh_i_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_i_Idx[1]

    ωs_ωref0_vref0_porder0_i_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_i_Idx[2]
        
    kwd_dyn_vh_θh_ωs_ωref0_vref0_etc_i =
        (; gens_vh_θh_i,
         ωs_ωref0_vref0_porder0_i )

    #--------------------------------------------------
    
    poi_dyn_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )

    pois_dyn =
        kwd_dyn_vh_θh_ωs_ωref0_vref0_etc_i    

    #-------------------------------------------------

    sim_fun_para_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )
        
    sim_fun_para  =
        vcat( gens_vh_θh_i,
             ωs_ωref0_vref0_porder0_i )
    
    #--------------------------------------------------

    sim_state_x0 =
        state[ gen_node_i_idx ]

    chunk_size = length(sim_state_x0 )
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0),
    #         chunk_size )

    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0) )

    stateDiffCache = similar(sim_state_x0)

    #--------------------------------------------------

    """
    poi_dyn_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )

    pois_dyn =
        (; gens_vh_θh_i,
         ωs_ωref0_vref0_porder0_i )
           
    """

    s_poi =  s_oop_poi = pois_dyn.ωs_ωref0_vref0_porder0_i
    
    poi_idx = poi_dyn_Idxs.ωs_ωref0_vref0_porder0_i_Idx

    poi_DiffCache = DiffCache(similar( s_poi ))

    #-------------------------------------------------
    
    sim_timespan =
        sim_timespan

    #---------------------------------------------------

    s_sim_fun_kwd_para =
        (; Ax_m,
         Bx_m,
         Cx_m,
         stateDiffCache,
         poi_DiffCache,
         sim_state_x0,
         plant_i,
         ra_Xd_dash_Xq_dash_i,
         pois_dyn,
         poi_dyn_Idxs )
    
    t_sim_fun_kwd_para =
        (; Ax_m,
         Bx_m,
         Cx_m,
         stateDiffCache,
         sim_state_x0,
         plant_i,
         sim_fun_para_Idxs,
         ra_Xd_dash_Xq_dash_i )
    
    #------------------------------------------------
    
    #--------------------------------------------------
    # dynamics simulation
    #--------------------------------------------------
            
    im_sym_i =
        im_sym[ gen_node_i_idx ]

    im_mass_matrix_i =
        im_mass_matrix[gen_node_i_idx,
                       gen_node_i_idx]

    #---------------------------------------------------
    # case 1
    #---------------------------------------------------

    t_sim_func! = t_ode_one_im_model_func!
    
    #-------------------------------------------------   
        
    t_sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> t_sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para = t_sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix_i,
        syms = im_sym_i )
    
    t_sim_prob = ODEProblem(
        t_sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )
    
    #--------------------------------------------------
    # simulate 
    #--------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    t_sim_sol = DifferentialEquations.solve(
        t_sim_prob, ode_alg, dt=dt  )

    #-------------------------------------------------
    # plot
    #-------------------------------------------------
    
    bus_i_name = network_bus_names[ node_i ]

    t_plot_δ_ω_ed_dash_eq_dash =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = t_sim_sol,
        node_syms_labels = im_sym_i,
        bus_name = bus_i_name,
        vars = [ :δ, :ω, :ed_dash, :eq_dash ],
        tspan = (0.0, 10.0),
        fmt = :png )

    p_δ_ω_ed_dash_eq_dash_idxs =
        last.(get_a_node_state_algb_vars_indices_in_syms(
            ; node_syms_labels = im_sym_i,
            bus_name = bus_i_name,
            vars = [ :δ, :ω, :ed_dash, :eq_dash ] ))


    #-------------------------------------------------
    # sensitivity
    #-------------------------------------------------

    # case 1
    
    t_sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            t_sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para = t_sim_fun_kwd_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg = ForwardDiffSensitivity() )


    t_sa_sim_sol = DifferentialEquations.solve(
        t_sens_prob , ode_alg, dt=dt  )


    #--------------------------------------------------
    # case 2
    #--------------------------------------------------

    s_sim_func! = s_ode_one_im_model_func!
    
    #-------------------------------------------------   
        
    s_sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> s_sim_func!(
            dx, x, p, t;
            poi_idx = poi_idx,
            sim_fun_kwd_para = s_sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix_i,
        syms = im_sym_i )
    
    s_sim_prob = ODEProblem(
        s_sim_ode_func,
        sim_state_x0,
        sim_timespan,
        s_poi  )
    
    #-------------------------------------------------
    # simulate 
    #--------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    s_sim_sol = DifferentialEquations.solve(
        s_sim_prob, ode_alg, dt=dt  )

    #--------------------------------------------------
    # plot
    #--------------------------------------------------
    
    bus_i_name = network_bus_names[ node_i ]

    s_plot_δ_ω_ed_dash_eq_dash =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = s_sim_sol,
        node_syms_labels = im_sym_i,
        bus_name = bus_i_name,
        vars = [ :δ, :ω, :ed_dash, :eq_dash ],
        tspan = (0.0, 10.0),
        fmt = :png )

    p_δ_ω_ed_dash_eq_dash_idxs =
        last.(get_a_node_state_algb_vars_indices_in_syms(
            ; node_syms_labels = im_sym_i,
            bus_name = bus_i_name,
            vars = [ :δ, :ω, :ed_dash, :eq_dash ] ))


    #-------------------------------------------------
    # sensitivity
    #-------------------------------------------------

    # case 2
    
    s_sa_sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            s_sim_func!(
                dx, x, p, t;
                poi_idx = poi_idx,
                sim_fun_kwd_para = s_sim_fun_kwd_para ),
        sim_state_x0,
        sim_timespan,
        s_poi; sensealg = ForwardDiffSensitivity() )


    s_sa_sim_sol = DifferentialEquations.solve(
        s_sa_sens_prob , ode_alg, dt=dt  )

    #--------------------------------------------------
    #--------------------------------------------------

    """

    https://stackoverflow.com/questions/74653454/why-am-i-getting-a-mutating-arrays-is-not-supported-error-here

    https://discourse.julialang.org/t/sensitivities-with-respect-to-initial-conditions-in-differentialequations-jl/25555/12
    https://github.com/FluxML/Tracker.jl

    https://discourse.julialang.org/t/discrete-adjoint-sensitivity-analysis-for-odes-in-differentialequations-jl/100007/3
    https://docs.sciml.ai/Overview/dev/highlevels/array_libraries/

    """
    
    return nothing

      
end



function t_v2_simulate_system_dynamic_im_model(
    ; dynamics_case =
        dynamics_case,
    only_gen =
        false,
    sim_timespan =
        sim_timespan,
    ode_alg =
        ode_alg,
    dyn_global_pf_options =
        dyn_global_pf_options,
    sta_global_pf_options =
        sta_global_pf_options,
    with_cb =
        false
        )

    # ------------------------------------------------
    # Diagnosis
    # ------------------------------------------------ 

    # ImplicitMidpoint(), ImplicitMidpoint(autodiff=false),
    
    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)

    sim_model_type        = "im-model"
    
    # ode_alg             = Rodas4()
    # algr_name           = "rodas4" 
    
    ode_alg               = ImplicitMidpoint()    
    algr_name             = "ImplicitMidpoint"

    dt                    = 0.01

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = false

    with_cb  = false

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
        
    #--------------------------------------------------
    #--------------------------------------------------
    
    netd  = NetworkData( dynamics_case()... )
    
    #-------------------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #--------------------------------------------------
    
    (; network_bus_names,
       non_gens_bus_names,
       gens_bus_names,
       Loads_bus_names,
       Trans_bus_names) =
           make_case_buses_names(
               netd.nodes )
    
    net_class_names =
        (;network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #--------------------------------------------------

    (; net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #-------------------------------------------------
    
    # im_sym, im_mass_matrix =
    #     generate_im_sym_and_mass_matrix(
    #         ; nodes = netd.nodes )

    (; im_model_ode_mass_matrix,
     im_model_pf_mass_matrix) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes = netd.nodes )

    (; gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes = netd.nodes  )

    im_sym =
        gens_nodes_im_vars_labels

    im_mass_matrix =
        im_model_ode_mass_matrix
        
    #-------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    
    #--------------------------------------------------

    (; nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para, im_idq_pf_cal)  =
         v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )
    
    #-------------------------------------------------
    
    (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para

    gens_vh_θh_view =
        im_para_aux_inputs.gens_vh_θh_view
    
    # ------------------------------------------------

    pf_net_param, sd_pf_views, _ =
        global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx =
        pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
        pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
        sd_pf_views


    (; each_gens_im_vars_Idx_in_state,
     im_vars_Idx_in_state,
     im_vars_view_in_state ) =
         para_update_gen_Ax_aux
    
    #-------------------------------------------------
    # Stability
    #-------------------------------------------------  
    
     (Ax_Bx_Cx_views,
      Ax_Bx_Cx_matrix,
      states_and_mat_Idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

     (nodes_state_Idx, Bx_idxs, Cx_idxs,
      id_iq_ph_vh_idxs, ω_ref_ωref_v_ref_idxs)  =
          states_and_mat_Idxs

    sys_states_Idxs_and_mat_Idxs =
        (; nodes_state_Idx,
         Bx_idxs, Cx_idxs,
         id_iq_ph_vh_idxs,
         ω_ref_ωref_v_ref_idxs )
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )

    #------------------------------------------------

    Ax_m = vec_Ax_views

    Bx_m = vec_Bx_views

    Cx_m = vec_Cx_views
    
    #-----------------------------------------------    

    # plant_i_eigvalues =
    #     eigvals( Matrix( vec_Ax_views[i] ))

    #-----------------------------------------------
    #-----------------------------------------------

    """

    #------------------------------------------------
    # case 0
    #------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i]) for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ]

    gens_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gens_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #--------------------------------------------------
    
   dims_vh_θh_ωs_ωref0_vref0_porder0 =  length.( [ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ] )

    system_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), [dims_vh_θh_ωs_ωref0_vref0_porder0]  )

    vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(system_size_offset_Idx)
    
    sim_fun_system_para_Idxs =
        vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #--------------------------------------------------

    system_ode_para = [[ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ]...;]
    
    #--------------------------------------------------

    sim_state_x0 =
        state[im_vars_Idx_in_state]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0),
    #         chunk_size )
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0) )

    #-------------------------------------------------

    stateDiffCache_gens =
        [ stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]

    #--------------------------------------------------
    
    sim_fun_system_kwd_para =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection,
         gens_sim_fun_gens_para_Idxs,
         sim_fun_system_para_Idxs,         
         ra_Xd_dash_Xq_dash_view,
         sys_states_Idxs_and_mat_Idxs,
         para_update_gen_Ax_aux)

    #-------------------------------------------------
    
    sim_func! = system_ode_im_model_func!

    sim_fun_para = system_ode_para
    
    #------------------------------------------------- 
        
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para = sim_fun_system_kwd_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )
    
    #-------------------------------------------------
    # simulate 
    #-------------------------------------------------

    # system_sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, ode_alg, dt=dt  )

    #-------------------------------------------------
    # sensitivity
    #-------------------------------------------------
    
    sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para =
                    sim_fun_system_kwd_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg = ForwardDiffSensitivity() )


    sa_sim_sol = DifferentialEquations.solve(
        sens_prob , ode_alg, dt=dt  )

    #-------------------------------------------------
    # plot
    #--------------------------------------------------

    plot_idxs_a = get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = system_sim_sol,
        node_syms_labels = im_sym,
        bus_name = "bus1",
        vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    # plot_a = make_plot_of_buses_volt_mag(
    #     ; sol = system_sim_sol,
    #     network_vars_labels = im_sym,
    #     nodes_name = ["bus1", "bus2"],
    #     vars = [:u_r, :u_i],
    #     tspan = sim_timespan,
    #     fmt = :png)
    

    # plot_b = make_plot_of_buses_volt_angle(
    #     ; sol = system_sim_sol,
    #     network_vars_labels = im_sym,
    #     nodes_name = ["bus1", "bus2"],
    #     vars = [:u_r, :u_i],
    #     tspan = sim_timespan,
    #     fmt = :png)

    """

    #-------------------------------------------------
    # case 1
    #-------------------------------------------------
    
    dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [
        length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i])
        for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in
            zip( gens_vh_θh_view,
                 gens_nodes_ωs_ωref0_vref0_porder0_view) ]

    gens_size_offset_Idx =  map(
        (dim) -> create_size_offset_Idx(
            dim; counter = 0 ),
        dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gen_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-------------------------------------------------
    
    
    id_iq_pg_vh =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view

    ωs_ωref0_vref0_porder0 =
        gens_nodes_ωs_ωref0_vref0_porder0_view
        
    gens_vh_θh = gens_vh_θh_view
    
    #--------------------------------------------------

    dims_nodes_id_iq_pg_vh =
        length.( id_iq_pg_vh )

    _,_, nodes_id_iq_pg_vh_idx_in_Idx =
        create_size_offset_Idx(
            dims_nodes_id_iq_pg_vh;
        counter = 0)

    #------------------------------------------------

    dims_nodes_ωs_ωref0_vref0_porder0 =
        length.( ωs_ωref0_vref0_porder0 )

    _,_, nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        create_size_offset_Idx(
            dims_nodes_ωs_ωref0_vref0_porder0;
        counter = 0)

    #-------------------------------------------------

    dims_nodes_vh_θh =
        length.( gens_vh_θh )

    _,_, nodes_vh_θh_indx_in_Idx =
        create_size_offset_Idx(
            dims_nodes_vh_θh;
        counter = 0)
    
    #--------------------------------------------------

    f_gens_vh_θh =
        [ gens_vh_θh...;]
    
    f_ωs_ωref0_vref0_porder0 =
        [ωs_ωref0_vref0_porder0...;]
    
    dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =
        length.(
            [f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0])

    _,_, vh_θh_ωs_ωref0_vref0_porder0_Idx =
        create_size_offset_Idx(
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0;
            counter = 0)

    vh_θh_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[1]

    ωs_ωref0_vref0_porder0_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[2]
        
    kwd_dyn_vh_θh_ωs_ωref0_vref0_etc =
        (; f_gens_vh_θh,
         f_ωs_ωref0_vref0_porder0 )

    #------------------------------------------------

    sim_fun_system_ode_flat_agg_para_Idxs =
        (; vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx )
        
    sim_fun_system_ode_flat_agg_para  =
        vcat( f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0 )
    

    #--------------------------------------------------

    sim_state_x0 =
        state[im_vars_Idx_in_state]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)

    #-------------------------------------------------

    stateDiffCache_gens =
        [stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]
    
    #-------------------------------------------------
    
    # sim_fun_system_kwd_flat_agg_para =
    #     (;
    #      vec_Ax_views,
    #      vec_Bx_views,
    #      vec_Cx_views,
    #      stateDiffCache_gens,
    #      sim_state_x0_gens,
    #      gens_nodes_collection ,
    #      gens_sim_fun_gen_para_Idxs,
    #      sim_fun_system_ode_flat_agg_para_Idxs,         
    #      ra_Xd_dash_Xq_dash_view,
         
    #      sys_states_Idxs_and_mat_Idxs,         
    #      nodes_id_iq_pg_vh_idx_in_Idx,
         
    #      nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
    #      nodes_vh_θh_indx_in_Idx )

    
    gen_nodes_ra_Xd_dash_Xq_dash_view =
        ra_Xd_dash_Xq_dash_view,

    gens_nodes_id_iq_pg_vh_idx_in_Idx =
        nodes_id_iq_pg_vh_idx_in_Idx
    
    gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        nodes_ωs_ωref0_vref0_porder0_idx_in_Idx
    
    gens_nodes_vh_θh_indx_in_Idx  =
        nodes_vh_θh_indx_in_Idx   
    
    sim_fun_system_kwd_flat_agg_para =
        (;
         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,

         gen_nodes_ra_Xd_dash_Xq_dash_view,
         
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection ,
         gens_sim_fun_gen_para_Idxs,
         
         sim_fun_system_ode_flat_agg_para_Idxs, 
         sys_states_Idxs_and_mat_Idxs,         
         gens_nodes_id_iq_pg_vh_idx_in_Idx,
         gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         gens_nodes_vh_θh_indx_in_Idx )
    
    #-------------------------------------------------
    
    sim_func! =
        system_flat_agg_ode_im_model_func!

    sim_fun_para =
        sim_fun_system_ode_flat_agg_para
    
    #--------------------------------------------------
        
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para =
                sim_fun_system_kwd_flat_agg_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )

    #-------------------------------------------------
    # simulate 
    #-------------------------------------------------

    # system_sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg )

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, ode_alg, dt=dt  )

    #--------------------------------------------------
    # sensitivity
    #--------------------------------------------------
    
    sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para =
                    sim_fun_system_kwd_flat_agg_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg =
            ForwardDiffSensitivity() )

    sim_sol = DifferentialEquations.solve(
        sens_prob, ode_alg, dt=dt  )
    
    #---------------------------------------------------
    # plot
    #--------------------------------------------------

    plot_idxs_a = get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = system_sim_sol,
        node_syms_labels = im_sym,
        bus_name = "bus1",
        vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    #-------------------------------------------------
    # case 2
    #-------------------------------------------------
    
    poi_dyn_Idxs =
        (; vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx )

    pois_dyn =
        (; f_gens_vh_θh,
         f_ωs_ωref0_vref0_porder0 )   

    #-------------------------------------------------
    
    sim_fun_system_kwd_para_wt_poi =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views, 
         stateDiffCache,
         poi_DiffCache,
         sim_state_x0,
         gens_nodes_collection ,
         ra_Xd_dash_Xq_dash,
         pois_dyn,
         poi_dyn_Idxs,
         states_and_mat_Idxs,         
         nodes_id_iq_pg_vh_idx_in_Idx,
         nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         nodes_vh_θh_indx_in_Idx
)
    
    #------------------------------------------------

    """
    poi_dyn_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )

    pois_dyn =
        (; gens_vh_θh_i,
         ωs_ωref0_vref0_porder0_i )
           
    """

    s_poi =  s_oop_poi = pois_dyn.ωs_ωref0_vref0_porder0
    
    poi_idx = poi_dyn_Idxs.ωs_ωref0_vref0_porder0_Idx

    poi_DiffCache = DiffCache(similar( s_poi ))

    #--------------------------------------------------
    
    sim_timespan =
        sim_timespan
    #----------------------------------------------
    # case 3
    #------------------------------------------------
    #------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------

    #------------------------------------------------
    # case 1
    #-------------------------------------------------

    t_sim_func! = t_ode_im_model_func!
    
    #------------------------------------------------   
        
    t_sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> t_sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para = t_sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
    t_sim_prob = ODEProblem(
        t_sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )
    
    #-------------------------------------------------
    # simulate 
    #-------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    t_sim_sol = DifferentialEquations.solve(
        t_sim_prob, ode_alg, dt=dt  )

    #------------------------------------------------
    # plot
    #------------------------------------------------
    
    bus_i_name = network_bus_names[ node_i ]

    t_plot_δ_ω_ed_dash_eq_dash =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = t_sim_sol,
        node_syms_labels = im_sym_i,
        bus_name = bus_i_name,
        vars = [ :δ, :ω, :ed_dash, :eq_dash ],
        tspan = (0.0, 10.0),
        fmt = :png )

    p_δ_ω_ed_dash_eq_dash_idxs =
        last.(get_a_node_state_algb_vars_indices_in_syms(
            ; node_syms_labels = im_sym_i,
            bus_name = bus_i_name,
            vars = [ :δ, :ω, :ed_dash, :eq_dash ] ))


    #-----------------------------------------------------
    # sensitivity
    #-----------------------------------------------------

    # case 1
    
    t_sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            t_sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para = t_sim_fun_kwd_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg = ForwardDiffSensitivity() )


    t_sa_sim_sol = DifferentialEquations.solve(
        t_sens_prob , ode_alg, dt=dt  )


    #-----------------------------------------------------
    # case 2
    #-----------------------------------------------------

    s_sim_func! = s_ode_one_im_model_func!
    
    #----------------------------------------------------   
        
    s_sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> s_sim_func!(
            dx, x, p, t;
            poi_idx = poi_idx,
            sim_fun_kwd_para = s_sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix_i,
        syms = im_sym_i )
    
    s_sim_prob = ODEProblem(
        s_sim_ode_func,
        sim_state_x0,
        sim_timespan,
        s_poi  )
    
    #-----------------------------------------------------
    # simulate 
    #-----------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    s_sim_sol = DifferentialEquations.solve(
        s_sim_prob, ode_alg, dt=dt  )

    #-----------------------------------------------------
    # plot
    #-----------------------------------------------------
    
    bus_i_name = network_bus_names[ node_i ]

    s_plot_δ_ω_ed_dash_eq_dash =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = s_sim_sol,
        node_syms_labels = im_sym_i,
        bus_name = bus_i_name,
        vars = [ :δ, :ω, :ed_dash, :eq_dash ],
        tspan = (0.0, 10.0),
        fmt = :png )

    p_δ_ω_ed_dash_eq_dash_idxs =
        last.(get_a_node_state_algb_vars_indices_in_syms(
            ; node_syms_labels = im_sym_i,
            bus_name = bus_i_name,
            vars = [ :δ, :ω, :ed_dash, :eq_dash ] ))


    #-----------------------------------------------------
    # sensitivity
    #-----------------------------------------------------

    # case 2
    
    s_sa_sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            s_sim_func!(
                dx, x, p, t;
                poi_idx = poi_idx,
                sim_fun_kwd_para = s_sim_fun_kwd_para ),
        sim_state_x0,
        sim_timespan,
        s_poi; sensealg = ForwardDiffSensitivity() )


    s_sa_sim_sol = DifferentialEquations.solve(
        s_sa_sens_prob , ode_alg, dt=dt  )

    #-----------------------------------------------------
    #-----------------------------------------------------

    """

    https://stackoverflow.com/questions/74653454/why-am-i-getting-a-mutating-arrays-is-not-supported-error-here

    https://discourse.julialang.org/t/sensitivities-with-respect-to-initial-conditions-in-differentialequations-jl/25555/12
    https://github.com/FluxML/Tracker.jl

    https://discourse.julialang.org/t/discrete-adjoint-sensitivity-analysis-for-odes-in-differentialequations-jl/100007/3
    https://docs.sciml.ai/Overview/dev/highlevels/array_libraries/

    """
    
    return nothing

      
end
                         

#-----------------------------------------------------
#-----------------------------------------------------
#  Drivers
#-----------------------------------------------------
#-----------------------------------------------------

function driver_create_gens_nodes_aggregate_system_matrices_Idx_and_views()

    only_gen = false
    
    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    #-----------------------------------------------------

    netd  =
        NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )
        
    #-----------------------------------------------------
    # Stability
    #----------------------------------------------------- 

    if only_gen == false
        (; Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) = create_gens_nodes_aggregate_system_matrices_Idx_and_views(
             gens_nodes_collection;
             only_gen = only_gen )
        

        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views )  = Ax_Bx_Cx_views
        
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix ) = Ax_Bx_Cx_matrix

        (; nodes_state_Idx,
         Bx_idxs,
         Cx_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen = only_gen,
            vec_Ax_τm_vf_views = nothing )
        
    else
        (; Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
             create_gens_nodes_aggregate_system_matrices_Idx_and_views(
                 gens_nodes_collection;
                 only_gen = only_gen )
        
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         vec_Ax_τm_vf_views) =
             Ax_Bx_Cx_views
        
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix) =
             Ax_Bx_Cx_matrix
        
        (; nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         τm_vf_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen = only_gen,
            vec_Ax_τm_vf_views = vec_Ax_τm_vf_views )
        
    end

    #-----------------------------------------------------
    #-----------------------------------------------------
    
    return vec_Ax_views, Ax_matrix
      
end


"""

vec_Ax_views, Ax_matrix = driver_create_gens_nodes_aggregate_system_matrices_Idx_and_views()

"""


#-------------------------------------------------------
#  Industrial
#-------------------------------------------------------



function driver_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init(  )

     dynamics_case = case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix

    

   netd  = NetworkData( dynamics_case()... )

    (;
     nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     industrial_model_misc_Idx,
     para_update_gen_Ax_aux,
     industrial_model_para_aux_inputs,
     industrial_model_pf_para,
     industrial_model_ωs_τm_vref_vhθh_idq,
     industrial_model_dyn_pf_up_para,
     industrial_model_idq_pf_cal,
     gens_nodes_τm_vf) =
        get_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init(
            netd; maxiter=40, ftol=1000*eps() ,
            xtol=1000*eps(), init_pf = true,
            with_δ_ed_eq = false )

    return nothing
    
end


#---------------------------------------------------


function driver_get_industrial_model_pf_param_views_and_init( )

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd =
        NetworkData( dynamics_case()... )

    only_gen = true

    (dict_sys_to_industry,
     pure_states_Idx_in_system,
     ur_ui_Idx_in_system,
     industrial_Idx,
     industrial_model_pure_states_Idx,
     industrial_model_ur_ui_Idx,
     pure_states_and_ur_ui_Idx_in_system,
     industrial_model_state_Idx,
     net_to_industrial_idx_conversion_dict,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_δ_ed_dash_eq_dash_Idxs_in_state,
     industrial_model_each_gen_nodes_pure_states_idx_in_state,
     industrial_model_each_gen_nodes_stab_states_idx_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     gen_nodes_dyn_param_view,
     gen_nodes_sub_dyn_param_view,
     pf_net_param,
     ra_Xd_dash_Xq_dash_view,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     industrial_model_pure_states_view_in_state,
     state_view,
     pf_state,
     nodes_u_view,
     nodes_pf_U_view,
     x0_vh_θh,
     working_vh_θh_view,
     red_vh_θh_0_view,
     mismatch,
     Jac_vh_θh,
     Inet,
     Inet_view,
     Iinj,
     Iinj_view,
     idq_wt_pad,
     idq_wt_pad_view,
     global_pf_views,
     sd_pf_views,
     global_pf_param,
     branches_name,
     nodes_name,
     named_tup_pf_result,
     bus_dict_init,
     branch_dict_init,
     state,
     gens_idx,
     nodes_u_Idx) =
        get_industrial_model_pf_param_views_and_init(
            netd; maxiter=40,
            ftol=1000*eps() ,
            xtol=1000*eps(), init_pf = true,
            with_δ_ed_eq = false,
            only_gen = only_gen )

    return nothing
    
end


function driver_get_industrial_model_pf_param_views_and_init_with_or_no_controller( )

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd  =
        NetworkData( dynamics_case()... )

    only_gen = true

    (;nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     industrial_model_misc_Idx,
     para_update_gen_Ax_aux,
     industrial_model_para_aux_inputs,
     industrial_model_pf_para,
     industrial_model_ωs_τm_vref_vhθh_idq ,
     industrial_model_dyn_pf_up_para,
     industrial_model_idq_pf_cal) =
        get_industrial_model_pf_param_views_and_init_with_or_no_controller(
            netd; maxiter=40, ftol=1000*eps() ,
            xtol=1000*eps(), init_pf = true,
            with_δ_ed_eq = false,
            only_gen = only_gen )

               
    return nothing
    
end


"""
Driver for"
       `dynamics_industrial_model`
Test:
     case_IEEE_9_Bus_sauer_dynamic_plant_SM_v6_P_t2_system_matrix,

     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
"""

#-----------------------------------------

function driver_simulation_dynamics_industrial_model()


    base_dir = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan = (0.0, 10.0)
    sim_timespan = (0.0, 10.0)
    ode_alg = Rodas4()
    sim_model_type = "industrial-model"

    only_gen = false

    algr_name = "rodas4"

    dyn_global_pf_options = (
        ; maxiter = 40,
        ftol=1000*eps() ,
        xtol=1000*eps(),
        with_δ_ed_eq = true )

    sta_global_pf_options = (
        ; maxiter=40, ftol=1000*eps(),
        xtol=1000*eps(),
        init_pf = true,
        with_δ_ed_eq = false )
    

    # list_dynamics_case = [        
    #     case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
    #     case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,
    #     case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix,
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad,
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P,
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_t2_millano,
    #    ]

    list_dynamics_case = [
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix,
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P,
        case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer

  ]

    list_sim_plots = [ ] 
    
    # Threads.@threads
    for a_dynamics_case in list_dynamics_case

        (sim_sol,
         net_names_labels_syms) =
             simulate_a_dynamics_industrial_model(
                 ; dynamics_case =
                     a_dynamics_case,
                 only_gen =
                     only_gen,
                 sim_timespan  =
                     sim_timespan,
                 ode_alg =
                     ode_alg,
                 dyn_global_pf_options =
                     dyn_global_pf_options,
                 sta_global_pf_options =
                     sta_global_pf_options )
        
        
        a_list_sim_plots =
            make_plots_for_industrial_model(
                ; case_fun =
                    a_dynamics_case,
                sim_model_type =
                    sim_model_type,
                sim_sol =
                    sim_sol,
                para_net_names_labels_syms =
                    net_names_labels_syms,
                tspan =
                    plot_tspan,
                base_dir =
                    base_dir,
                algr_name =
                    algr_name )
        
        push!( list_sim_plots, a_list_sim_plots )
    end

    return nothing
    
    # return list_sim_plots
        
end


#----------------------------------------------------
 

"""

Test driver for:

`simulate_a_dynamics_industrial_model_with_or_no_controllers`

    list_dynamics_case = [        
        case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix,
        case_IEEE_14_Bus_dynamic_plants_v6_P_rscad,
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P,
        case_IEEE_14_Bus_dynamic_plants_v6_P_t2_millano ]

    list_dynamics_case = [
        case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix
  ]


"""
function driver_simulation_dynamics_industrial_model_with_or_no_controllers()

    base_dir = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan = (0.0, 10.0)
    sim_timespan = (0.0, 10.0)
    ode_alg = Rodas4()
    sim_model_type = "industrial-model"
    algr_name = "rodas4"

    dyn_global_pf_options = (
        ; maxiter =
            40,
        ftol =
            1000*eps(),
        xtol =
            1000*eps(),
        with_δ_ed_eq =
            true )

    sta_global_pf_options = (
        ; maxiter = 40,
        ftol = 1000*eps(),
        xtol = 1000*eps(),
        init_pf = true,
        with_δ_ed_eq = false )

    only_gen = true
    
    with_cb  = false
    
    list_dynamics_case = [
        case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix, case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P ]

    list_sim_plots = [ ]
    
    for a_dynamics_case in list_dynamics_case
        
        sim_sol, net_names_labels_syms =
            simulate_a_dynamics_industrial_model_with_or_no_controllers( ; dynamics_case         = a_dynamics_case,
              only_gen              = only_gen,
              sim_timespan          = sim_timespan,
              ode_alg               = ode_alg,
              dyn_global_pf_options = dyn_global_pf_options,
            sta_global_pf_options = sta_global_pf_options,
        with_cb  = with_cb )
        
        a_list_sim_plots =
            make_plots_for_industrial_model(
                ; case_fun =
                    a_dynamics_case,
                sim_model_type =
                    sim_model_type,
                sim_sol =
                    sim_sol,
                para_net_names_labels_syms =
                    net_names_labels_syms,
                tspan =
                    plot_tspan,
                base_dir =
                    base_dir,
                algr_name =
                    algr_name)
        
        push!( list_sim_plots, a_list_sim_plots )
    end
    
    return list_sim_plots
    
    # return nothing

        
end


#-------------------------------------------------------
#  v2 im
#-------------------------------------------------------


function driver_v2_get_im_model_pf_param_views_and_init( )

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd  = NetworkData( dynamics_case()... )

    only_gen = false

    (;im_vars_indices_in_system,
     pure_states_Idx_in_system,
     im_algebraic_vars_Idx_in_system,
     ur_ui_Idx_in_system,
     im_vars_and_ur_ui_Idx_in_system,
     im_vars_Idx_in_state,
     nodes_ur_ui_Idx_in_state,
     im_state_Idx,
     each_gens_im_vars_Idx_in_state,
     net_to_im_idx_conversion_dict,
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_δ_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     gen_nodes_dyn_param_view,
     gen_nodes_sub_dyn_param_view,
     pf_net_param,
     ra_Xd_dash_Xq_dash_view,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     im_vars_view_in_state,
     state_view,
     pf_state,
     nodes_u_view,
     nodes_pf_U_view,
     x0_vh_θh,
     working_vh_θh_view,
     red_vh_θh_0_view,
     mismatch,
     Jac_vh_θh,
     Inet,
     Inet_view,
     Iinj,
     Iinj_view,
     idq_wt_pad,
     idq_wt_pad_view,
     global_pf_views,
     sd_pf_views,
     global_pf_param,
     branches_name,
     nodes_name,
     named_tup_pf_result,
     bus_dict_init,
     branch_dict_init,
     state,
     gens_idx,
     nodes_u_Idx) =
        v2_get_im_model_pf_param_views_and_init(
            netd;
            maxiter=40,
            ftol=1000*eps() ,
            xtol=1000*eps(),
            init_pf = true,
            with_δ_ed_eq = false,
            only_gen = only_gen )

    return nothing
    
end


function driver_v2_get_im_sys_pf_param_views_and_init( )

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd  = NetworkData( dynamics_case()... )

    only_gen = false

    (;nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para,
     im_idq_pf_cal) =
        v2_get_im_sys_pf_param_views_and_init(
            netd; maxiter=40,
            ftol=1000*eps() ,
            xtol=1000*eps(),
            init_pf = true,
            with_δ_ed_eq = false,
            only_gen = only_gen )

               
    return nothing
    
end



function driver_v2_simulate_a_dynamics_im_model()

    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)
    ode_alg               = Rodas4()
    sim_model_type        = "im-model"
    algr_name             = "rodas4"

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = false

    with_cb  = false
    

    list_sim_plots = []
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer

    sim_sol, para_net_names_labels_syms =
        v2_simulate_a_dynamics_im_model(
            ; dynamics_case =
                dynamics_case,
            only_gen =
                false,
            sim_timespan  =
                sim_timespan,
            ode_alg =
                ode_alg,
            dyn_global_pf_options =
                dyn_global_pf_options,
            sta_global_pf_options =
                sta_global_pf_options,
            with_cb =
                with_cb)

    a_list_sim_plots =
        make_plots_for_industrial_model(
            ; case_fun = dynamics_case,
            sim_model_type =
                sim_model_type,
            sim_sol =
                sim_sol,
            para_net_names_labels_syms =
                para_net_names_labels_syms,
            tspan =
                plot_tspan,
            base_dir =
                base_dir,
            algr_name =
                algr_name )

    push!( list_sim_plots,
           a_list_sim_plots )    

    # return nothing
    return list_sim_plots
end


function driver_v2_simulate_dynamics_im_model_cases()

    base_dir = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan = (0.0, 10.0)
    sim_timespan = (0.0, 10.0)
    ode_alg = Rodas4()
    sim_model_type = "im-model"
    algr_name = "rodas4"

    dyn_global_pf_options = (
        ; maxiter = 40,
        ftol = 1000*eps(),
        xtol = 1000*eps(),
        with_δ_ed_eq = true )

    sta_global_pf_options = (
        ; maxiter  = 40,
        ftol = 1000*eps(),
        xtol = 1000*eps(),
        init_pf = true,
        with_δ_ed_eq = false )

    only_gen = false
    
    with_cb  = false
    
    list_dynamics_case = [
        case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix,
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
  ]

    list_sim_plots = [ ]
    
    for a_dynamics_case in list_dynamics_case

        sim_sol, net_names_labels_syms =
            v2_simulate_a_dynamics_im_model(
                ; dynamics_case = a_dynamics_case,
                only_gen =
                    only_gen,
                sim_timespan  =
                    sim_timespan,
                ode_alg =
                    ode_alg,
                dyn_global_pf_options =
                    dyn_global_pf_options,
                sta_global_pf_options =
                    sta_global_pf_options,
                with_cb =
                    with_cb)
        

        a_list_sim_plots =
            make_plots_for_industrial_model(
                ; case_fun =
                    a_dynamics_case,
                  sim_model_type =
                      sim_model_type,
                  sim_sol =
                      sim_sol,
                  para_net_names_labels_syms =
                      net_names_labels_syms,
                  tspan =
                      plot_tspan,
                  base_dir =
                      base_dir,
                  algr_name =
                      algr_name )
        
        push!( list_sim_plots, a_list_sim_plots )
    end
    
    # return list_sim_plots
    
    # return nothing

        
end



function t_driver_create_gens_nodes_aggregate_system_matrices_Idx_and_views()

    only_gen = false
    
    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    #-----------------------------------------------------

    netd  =
        NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )
        
    #-----------------------------------------------------
    # Stability
    #----------------------------------------------------- 

    if only_gen == false
        (; Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
             create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )
        

        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views )  =
             Ax_Bx_Cx_views
        
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix ) =
             Ax_Bx_Cx_matrix

        (; nodes_state_Idx,
         Bx_idxs,
         Cx_idxs) =
             idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen = only_gen,
            vec_Ax_τm_vf_views = nothing )
        
    else
        (; Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
             create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )
        
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         vec_Ax_τm_vf_views) =
             Ax_Bx_Cx_views
        
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix) =
             Ax_Bx_Cx_matrix
        
        (; nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         τm_vf_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen = only_gen,
            vec_Ax_τm_vf_views =
                vec_Ax_τm_vf_views )
        
    end

    #-----------------------------------------------------
    #-----------------------------------------------------
    
    return vec_Ax_views, Ax_matrix
      
end


#-----------------------------------------------------
#-----------------------------------------------------
#  Simulation
#-----------------------------------------------------
#-----------------------------------------------------


function v2_simulate_one_plant_dynamics_im_model(
    ; dynamics_case =
        dynamics_case,
    only_gen =
        false,
    sim_timespan =
        sim_timespan,
    ode_alg =
        ode_alg,
    dyn_global_pf_options =
        dyn_global_pf_options,
    sta_global_pf_options =
        sta_global_pf_options,
    with_cb =
        false
        )

    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 

    # ImplicitMidpoint(), ImplicitMidpoint(autodiff=false),
    
    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)

    sim_model_type        = "im-model"
    
    # ode_alg              = Rodas4()
    # algr_name           = "rodas4" 
    
    ode_alg               = ImplicitMidpoint()    
    algr_name             = "ImplicitMidpoint"

    dt                    = 0.01

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = false

    with_cb  = false

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
        
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #-----------------------------------------------------
    
    ; (network_bus_names,
       non_gens_bus_names,
       gens_bus_names,
       Loads_bus_names,
       Trans_bus_names) =
           make_case_buses_names(
               netd.nodes )
    
    net_class_names =
        (;network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-----------------------------------------------------

    (; net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #-----------------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes = netd.nodes )

    #-----------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    
    #-----------------------------------------------------

    (; nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para, im_idq_pf_cal)  =
         v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )
    
    # -----------------------------------------------------

    pf_net_param, sd_pf_views, _ =
        global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx =
        pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
        pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
        sd_pf_views

    
    #-----------------------------------------------------
    # Stability
    #-----------------------------------------------------  
    
    (Ax_Bx_Cx_views,
     Ax_Bx_Cx_matrix,
     idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs)  = idxs
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )

    #-----------------------------------------------------
    
    node_i = 1
    
    gen_node_i_idx =
        nodes_state_Idx[ node_i ]
    
    #-----------------------------------------------------

    plant_i =
        gens_nodes_collection[ node_i ]      

    #-----------------------------------------------------

    plant_i_eigvalues =
        eigvals(Matrix(vec_Ax_views[ node_i ]))
    
    #-----------------------------------------------------
            
    im_sym_i =
        im_sym[ gen_node_i_idx ]

    im_mass_matrix_i =
        im_mass_matrix[gen_node_i_idx,
                       gen_node_i_idx]
    
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para

    ode_para =
        (; Ax_Bx_Cx_views,
        gens_dynamic_id_iq_pg_vh_by_vhθh_view,
        gens_nodes_ωs_ωref0_vref0_porder0_view,
        node_i )
    
    #-----------------------------------------------------

    ode_fun_para =
        (; Ax_matrix,
        Bx_matrix,
        Cx_matrix,
        gens_dynamic_id_iq_pg_vh_by_vhθh_view,
        gens_nodes_ωs_ωref0_vref0_porder0_view )

    #-----------------------------------------------------

    sim_state_x0 =
        state[ gen_node_i_idx ]
    
    sim_timespan =
        sim_timespan

    #-----------------------------------------------------

    stateDiffCache =
        similar(sim_state_x0)
    
    netd_nothing = nothing

    sim_fun_para = (
        ; netd_nothing,
        nodes_cb_sw,
        ode_para,
        plant_i,
        node_i,
        idxs,
        stateDiffCache
       )
    
    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------

    sim_func = dynamics_one_im_model!
    
    #----------------------------------------------------- 
    #-----------------------------------------------------
        
    sim_ode_func! = ODEFunction{true}(
        sim_func; mass_matrix = im_mass_matrix_i,
        syms = im_sym_i )
    
    sim_prob = ODEProblem(sim_ode_func!,
                          sim_state_x0,
                          sim_timespan,
                          sim_fun_para )
    
    # #-----------------------------------------------------

    # # sim_sol = DifferentialEquations.solve(
    # #    sim_prob, ode_alg  )

    # sim_sol = DifferentialEquations.solve(
    #     sim_prob, ode_alg, dt=dt  )

    # plot_idx = [1,2,3,4]
    
    # t_plot = plot( sim_sol, idxs = plot_idx)

    # #-----------------------------------------------------

    (; netd_nothing,
     nodes_cb_sw,
     ode_para,
     plant_i,
     node_i,
     idxs,
     stateDiffCache ) = sim_fun_para

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs) = idxs

    (Ax_Bx_Cx_views,
     gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view,
     node_i) =  ode_para 

    (; vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views)  =
        Ax_Bx_Cx_views

    
    vec_Ax_i =
        vec_Ax_views[ node_i ]

    #--------------------------------------------

    x  = sim_state_x0
    
    dx = similar(x)
    
    #--------------------------------------------

    stateDiffCache =
        get_tmp(stateDiffCache, x)

    stateDiffCache .= sim_state_x0
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    update_a_im_plant_system_matrices!(
        vec_Ax_i,
        stateDiffCache,
        plant_i )    

    #--------------------------------------------

    (t_Ax_Bx_Cx_views,
     t_gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     t_gens_nodes_ωs_ωref0_vref0_porder0_view, node_i) =
         ode_para

    (t_vec_Ax_views,
     t_vec_Bx_views,
     t_vec_Cx_views) =
         t_Ax_Bx_Cx_views
        
    t_Ax_m = t_vec_Ax_views[ node_i ]
    
    t_Bx_m = t_vec_Bx_views[ node_i ]
    
    t_Cx_m = t_vec_Cx_views[ node_i ]
    
    
    t_id_iq_pg_vh =
        t_gens_dynamic_id_iq_pg_vh_by_vhθh_view[ node_i ]

    t_ωs_ωref0_vref0_porder0 =
        t_gens_nodes_ωs_ωref0_vref0_porder0_view[ node_i ]

    #--------------------------------------------

    # ode_one_im_model_func!(
    #     dx,
    #     x,
    #     ode_one_fun_para ,
    #     t)

    #--------------------------------------------
    
    # #--------------------------------------------
    # # ode gens pure state
    # #--------------------------------------------

    # tt = 0.0
    
    # ode_one_im_model_func!(
    #     dx,
    #     x,
    #     ode_para ,
    #     tt)


    # ode_one_im_model_without_δ_func!(
    # dx,
    # x,
    # ode_para,    
    # t)
    
    # ode_one_im_model_with_only_δ_func!(
    # dx,
    # x,    
    # (stateDiffCache, ode_para),
    # t)

    # #--------------------------------------------
    # # Another method for all nodes votages
    # #--------------------------------------------

    # # im_all_nodes_voltage_terminal_func!(
    # #     dx,
    # #     x,
    # #     (nodes_pf_U_view, nodes_u_Idx_in_ranges),
    # #     t)
            
    # #--------------------------------------------
    # # non gens terminal voltages
    # #--------------------------------------------

    # non_gens_im_voltage_terminal_func!(
    #     dx, x, im_para_aux_inputs, t)
    
    # #--------------------------------------------
    # # gens terminal voltages
    # #--------------------------------------------

    # gens_im_voltage_terminal_func!(
    #     dx, x, im_para_aux_inputs, t)
    
    # #--------------------------------------------
    # #--------------------------------------------    
    
    return nothing

      
end
                         

#-----------------------------------------------------
#-----------------------------------------------------
# combines dynamic and powerflow
#-----------------------------------------------------
#-----------------------------------------------------


function combined_dyn_and_pf_simulation()

    
    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)

    sim_model_type        = "im-model"
    
    # ode_alg             = Rodas4()
    # algr_name           = "rodas4" 
    
    ode_alg               = ImplicitMidpoint()    
    algr_name             = "ImplicitMidpoint"   
    dt                    = 0.01

    pf_alg                = NewtonRaphson()

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = false

    with_cb  = false

    Idxs_type = :Idxs_hybrid  # :Idxs_im

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
        
    #----------------------------------------------
    #----------------------------------------------
    
    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )


    #-----------------------------------------------
    # pf init
    #----------------------------------------------- 

    """ other options are

            Idxs_type = :Idxs_industrial

            Idxs_type = :Idxs_im
    """
    
    (; nll_dyn_iip_pf_param, ) =
        get_integrated_nll_dyn_iip_pf_param(
            netd;
            Idxs_type = Idxs_type ) # :Idxs_im

    (; red_ΔPQ_Δidq,
     red_vh_θh_0_idq,         
     gens_uh_Q_from_red_sol_para,
     kwd_net_param ) =
         nll_dyn_iip_pf_param

    #------------------------------------------------ 
    #------------------------------------------------

    pf_idx_and_state =
        kwd_net_param.pf_net_param.pf_idx_and_state

    slack_vh =
        pf_idx_and_state.slack_vh
    
    gens_vh =
        pf_idx_and_state.gens_vh
    
    gens_Idx_and_vh =
        pf_idx_and_state.gens_Idx_and_vh
    
    non_slack_gens_Idx_and_vh =
        pf_idx_and_state.non_slack_gens_Idx_and_vh

    #-------------------------------------------------
    
    # pf_Idx = kwd_net_param.pf_net_param.pf_Idx

    # red_vh_θh_idx = pf_Idx.red_vh_θh_idx
    
    #------------------------------------------------

    pf_and_dyn_idx_and_Idx =
        kwd_net_param.pf_net_param.pf_and_dyn_idx_and_Idx
    
    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist

    #------------------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx

    slack_gens_nodes_idx =
        net_comp_type_idx.slack_gens_nodes_idx

    slack_bus_idx =
        slack_gens_nodes_idx

    non_slack_gens_nodes_idx =
        net_comp_type_idx.non_slack_gens_nodes_idx

    gens_idx = gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    # gens_with_loc_load_idx =
    #     net_comp_type_idx.gens_with_loc_load_idx

    gens_with_loc_load_idx =
        loc_load_exist == false ?
        [] : 
        net_comp_type_idx.gens_with_loc_load_idx

    all_nodes_idx =
        net_comp_type_idx.all_nodes_idx

    #------------------------------------------------

    nodes_size = sum(
        [length(gens_nodes_idx),
         length(non_gens_nodes_idx) ])

    #------------------------------------------------

    dict_n2s =
        net_comp_type_idx.dicts_net_to_streamlined_idx
    
    n2s_slack_gens_idx =
        dict_n2s.dict_n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        dict_n2s.dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx =
        dict_n2s.dict_n2s_gens_idx

    n2s_non_gens_idx =
        dict_n2s.dict_n2s_non_gens_idx

    n2s_gens_with_loc_load_idxs =
        dict_n2s.dict_n2s_gens_with_loc_load_idxs
    
    n2s_load_idx =
        dict_n2s.dict_n2s_load_idx
    
    n2s_transmission_idxs =
        dict_n2s.dict_n2s_transmission_idxs
    
    #-----------------------------------------------

    if Idxs_type == :Idxs_hybrid
        
        idx_and_Idx =
            pf_and_dyn_idx_and_Idx.hybrid_pf_etc_idx_and_Idx
    elseif Idxs_type == :Idxs_industrial

        idx_and_Idx =
            pf_and_dyn_idx_and_Idx.industrial_model_pf_idx_and_Idx
    elseif Idxs_type == :Idxs_im

        idx_and_Idx =
            pf_and_dyn_idx_and_Idx.im_model_pf_idx_and_Idx
    else
         nothing
    end

    red_vh_θh_idx =
        idx_and_Idx.red_vh_θh_idx
    
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

    #-----------------------------------------------

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx

    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash

    #----------------------------------------------
    # dyn init
    #----------------------------------------------
    
    (; network_bus_names,
       non_gens_bus_names,
       gens_bus_names,
       Loads_bus_names,
       Trans_bus_names) =
           make_case_buses_names(
               netd.nodes )
    
    net_class_names =
        (;network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #-----------------------------------------------

    (; net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ;nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #-----------------------------------------------
    
    # im_sym, im_mass_matrix =
    #     generate_im_sym_and_mass_matrix(
    #         ; nodes = netd.nodes )

    (; im_model_ode_mass_matrix,
     im_model_pf_mass_matrix) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes = netd.nodes )

    (; gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes = netd.nodes  )

    im_sym =
        gens_nodes_im_vars_labels

    im_mass_matrix =
        im_model_ode_mass_matrix
        
    #------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    
    #------------------------------------------------
    #------------------------------------------------
    
    (; nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para,
     im_idq_pf_cal )  =
         v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )
    
    #------------------------------------------------
    
    (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para

    # ------------------------------------------------

    (; im_nodes_voltage_Idx,
     im_vars_Idx_in_state, 
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     im_vars_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view ) =
         im_para_aux_inputs

    (; gens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges ) =
         im_nodes_voltage_Idx
    
    # -----------------------------------------------

    δ_ω_ed_dash_eq_dash =
        gen_nodes_δ_ω_ed_dash_eq_dash_views
    
    id_iq_pg_vh =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view

    ωs_ωref0_vref0_porder0 =
        gens_nodes_ωs_ωref0_vref0_porder0_view
        
    gens_vh_θh = gens_vh_θh_view
    
    # -----------------------------------------------

    # pf_net_param, _, _ = global_pf_param

    # _, _, _, _, _, pf_and_dyn_idx_and_Idx = pf_net_param

    # _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
    #     pf_param_views

    # _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
    #     sd_pf_views


    # (; each_gens_im_vars_Idx_in_state,
    #  im_vars_Idx_in_state,
    #  im_vars_view_in_state ) =
    #      para_update_gen_Ax_aux
    
    #-----------------------------------------------
    # Stability
    #----------------------------------------------- 
    
    (Ax_Bx_Cx_views,
     Ax_Bx_Cx_matrix,
     states_and_mat_Idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    (vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views)  =
        Ax_Bx_Cx_views
    
    (Ax_matrix,
     Bx_matrix,
     Cx_matrix) =
         Ax_Bx_Cx_matrix

    (;nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs)  =
         states_and_mat_Idxs

    sys_states_Idxs_and_mat_Idxs =
        (; nodes_state_Idx,
         Bx_idxs, Cx_idxs,
         id_iq_ph_vh_idxs,
         ω_ref_ωref_v_ref_idxs )
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )
    
    #------------------------------------------------
    #------------------------------------------------

    ## alternative
    
    # nodes_ur_ui = [state[idx]
    #                for idx in
    #                    nodes_u_Idx_in_ranges]
    
    # dyn_nodes_vh_θh = [
    #     abs.(x_from_xr_xi.( nodes_ur_ui ))...;
    #     angle.(x_from_xr_xi.( nodes_ur_ui ))... ]

    # dyn_nodes_vh_θh = [
    #     abs.(x_from_xr_xi.(nodes_pf_U_view))...;
    #     angle.(x_from_xr_xi.(nodes_pf_U_view))... ]

    # dyn_red_vh_θh =
    #     dyn_nodes_vh_θh[ red_vh_θh_idx ]

    """

    #-----------------------------------------------
    # case 0
    #-----------------------------------------------

    sim_state_x0 =
        state[im_vars_Idx_in_state]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0),
    #         chunk_size )
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0) )

    #-----------------------------------------------------

    stateDiffCache_gens =
        [ stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]

    #-----------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i]) for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh,  ωs_ωref0_vref0_porder0 ) ]

    gens_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gens_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =  length.( [ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh,  ωs_ωref0_vref0_porder0 ) ] )

    system_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), [dims_gens_vh_θh_ωs_ωref0_vref0_porder0]  )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(system_size_offset_Idx)
    
    sim_fun_system_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------------

    system_ode_para = [[ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh,  ωs_ωref0_vref0_porder0) ]...;]
    
    #-----------------------------------------------------
    
    sim_fun_system_kwd_para =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         # ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_dash_Xq_dash_view,
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection,
         gens_sim_fun_gens_para_Idxs,
         sim_fun_system_para_Idxs,         
         sys_states_Idxs_and_mat_Idxs,
         para_update_gen_Ax_aux )

    #-----------------------------------------------------
    
    sim_func! = system_ode_im_model_func!

    sim_fun_para = system_ode_para
    
    #----------------------------------------------------- 
        
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para = sim_fun_system_kwd_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )
    
    #-----------------------------------------------------
    # simulate 
    #-----------------------------------------------------

    # system_sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, ode_alg, dt=dt  )

    #-----------------------------------------------------
    # sensitivity
    #-----------------------------------------------------
    
    sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para =
                    sim_fun_system_kwd_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg = ForwardDiffSensitivity() )


    sa_sim_sol = DifferentialEquations.solve(
        sens_prob , ode_alg, dt=dt  )

    #-----------------------------------------------------
    # plot
    #-----------------------------------------------------

    plot_idxs_a = get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = system_sim_sol,
        node_syms_labels = im_sym,
        bus_name = "bus1",
        vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    # plot_a = make_plot_of_buses_volt_mag(
    #     ; sol = system_sim_sol,
    #     network_vars_labels = im_sym,
    #     nodes_name = ["bus1", "bus2"],
    #     vars = [:u_r, :u_i],
    #     tspan = sim_timespan,
    #     fmt = :png)
    

    # plot_b = make_plot_of_buses_volt_angle(
    #     ; sol = system_sim_sol,
    #     network_vars_labels = im_sym,
    #     nodes_name = ["bus1", "bus2"],
    #     vars = [:u_r, :u_i],
    #     tspan = sim_timespan,
    #     fmt = :png)

    """
    
    #-----------------------------------------------------
    # case 1
    #-----------------------------------------------------
    
    """
    sim_state_x0 =
        state[ im_vars_Idx_in_state ]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)

    #-----------------------------------------------------

    stateDiffCache_gens =
        [stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]

    #-----------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i]) for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh,  ωs_ωref0_vref0_porder0) ]

    gens_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gen_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------------

    dims_gens_nodes_id_iq_pg_vh =
        length.( id_iq_pg_vh )

    _,_, gens_nodes_id_iq_pg_vh_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_id_iq_pg_vh;
        counter = 0)

    #-----------------------------------------------------

    dims_gens_nodes_ωs_ωref0_vref0_porder0 =
        length.( ωs_ωref0_vref0_porder0 )

    _,_, gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_ωs_ωref0_vref0_porder0;
        counter = 0)

    #-----------------------------------------------------

    dims_gens_nodes_vh_θh =
        length.( gens_vh_θh )

    _,_, gens_nodes_vh_θh_indx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_vh_θh;
        counter = 0)
    
    #-----------------------------------------------------

    f_gens_vh_θh =
        [ gens_vh_θh...;]
    
    f_ωs_ωref0_vref0_porder0 =
        [ωs_ωref0_vref0_porder0...;]
    
    dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =
        length.(
            [f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0])

    _,_, vh_θh_ωs_ωref0_vref0_porder0_Idx =
        create_size_offset_Idx(
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0;
            counter = 0)

    vh_θh_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[1]

    ωs_ωref0_vref0_porder0_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[2]
        
    kwd_dyn_vh_θh_ωs_ωref0_vref0_etc =
        (; f_gens_vh_θh,
         f_ωs_ωref0_vref0_porder0 )

    #-----------------------------------------------------

    sim_fun_system_ode_flat_agg_para_Idxs =
        (; vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx )

    #-----------------------------------------------------
    
    sim_fun_system_kwd_flat_agg_para =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         # ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_dash_Xq_dash_view,         
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection ,
         gens_sim_fun_gen_para_Idxs,
         sim_fun_system_ode_flat_agg_para_Idxs,         
         sys_states_Idxs_and_mat_Idxs,         
         gens_nodes_id_iq_pg_vh_idx_in_Idx,
         gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         gens_nodes_vh_θh_indx_in_Idx )

    #-----------------------------------------------------
    
    sim_func! =
        system_flat_agg_ode_im_model_func!

    #-----------------------------------------------------
    
    sim_fun_para = sim_fun_system_ode_flat_agg_para  =
        vcat( f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0 )
    
    # sim_fun_para =
    #     sim_fun_system_ode_flat_agg_para

    #-----------------------------------------------------
    # dynamic ode func and prob
    #-----------------------------------------------------
            
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para =
                sim_fun_system_kwd_flat_agg_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
         sim_fun_para  )

    #-----------------------------------------------------
    # powerflow func 
    #-----------------------------------------------------
    
    pf_fun_mismatch =
        integrated_nll_dyn_pf_iip_ΔPQ_Δidq_mismatch

    pf_func = NonlinearFunction( ( g, x, p ) ->
        pf_fun_mismatch(
            g, x, p;
            kwd_net_param =
                kwd_net_param,
            Idxs_type =
                Idxs_type  ) )

    #-----------------------------------------------------

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, ode_alg, dt=dt )

    #-----------------------------------------------------
    
    # system_sim_sol = DifferentialEquations.solve(
    #     ODEProblem( ODEFunction{true}(
    #         (dx, x, p, t) ->
    #             sim_func!(
    #                 dx, x, p, t;
    #                 sim_fun_kwd_para =
    #                     sim_fun_system_kwd_flat_agg_para )
    #         ; mass_matrix = im_mass_matrix,
    #         syms = im_sym ),
    #                 sim_state_x0,
    #                 sim_timespan,
    #                 sim_fun_para ),
    #     ode_alg, dt=dt )
      
    #-----------------------------------------------------

   # Coupling betwing ode and pf, via
   # δ_ω_ed_dash_eq_dash system_sim_sol.u

    δ_ω_ed_dash_eq_dash =
        get_im_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            system_sim_sol.u[1,end],
            nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    #------------------------------------------------------

    pf_prob = NonlinearProblem(
        pf_func,
        red_vh_θh_0_idq,
        δ_ω_ed_dash_eq_dash ) 

    pf_sol = NonlinearSolve.solve(
        pf_prob,
        pf_alg )

    #------------------------------------------------------

    # pf_sol = NonlinearSolve.solve(
    #     NonlinearProblem(
    #         NonlinearFunction(
    #             ( g, x, p ) ->
    #                 pf_fun_mismatch(
    #                     g, x, p;
    #                     kwd_net_param =
    #                         kwd_net_param,
    #                     Idxs_type =
    #                         Idxs_type  )),
    #         red_vh_θh_0_idq,
    #         δ_ω_ed_dash_eq_dash ),
    #     pf_alg )

    #-----------------------------------------------

    red_vh_θh_idq =
        get_red_vh_θh_idq_from_pf_sol(
            pf_sol,
            δ_ω_ed_dash_eq_dash,
            gens_uh_Q_from_red_sol_para )
    
    #----------------------------------------------
    #----------------------------------------------
    # Integrator inferface
    #----------------------------------------------
    #----------------------------------------------

    time_final  = 10

    simtime  = time_final # +  dt

    grantedtime = 0

    # simulation_time_to = 0

    reinit!(integ)

    integ =
        init( sim_prob,
             ode_alg,
              dt = dt,
              abstol=1e-12,
              reltol=1e-12,
              tstops = [time_final],
              advance_to_tstop = true)

    dyn_ode_para = integ.p

    # vh_θh_Idx,     
    # ωs_ωref0_vref0_porder0_Idx
    
    #------------------------------------------

    result = []

    while grantedtime <= simtime

        global grantedtime, dyn_pf_δ_ω_ed_dash_eq_dash, pf_sol, red_vh_θh_idq

        global simulation_time_to, updated_pf_prob, pf_gens_vh_θh

        simulation_time_to = grantedtime + dt


        if simulation_time_to < time_final 

            add_tstop!(integ, simulation_time_to)

        else

            add_tstop!(integ, time_final)

        end

        step!(integ, dt, true)

        # Coupling of dynamic simmulation of ode to pf
        # through dyn_pf_δ_ω_ed_dash_eq_dash

        dyn_pf_δ_ω_ed_dash_eq_dash =
            [ integ.u[idx]
              for idx in
                  nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state  ]
        
        #----------------------------------------------
        
         updated_pf_prob = remake( pf_prob,
                       u0 = red_vh_θh_idq,
                       p = dyn_pf_δ_ω_ed_dash_eq_dash )

        pf_sol = NonlinearSolve.solve(
            updated_pf_prob, pf_alg )
        
        #----------------------------------------------

        # Coupling of powerflow simmulation to
        # dynamic simulation through vec_pf_gens_vh_θh.

        # vec_pf_gens_vh_θh

        pf_gens_vh_θh =
            get_pf_gens_uh_from_pf_sol(
                pf_sol, 
                gens_uh_Q_from_red_sol_para )

        # ode para update
        
        dyn_ode_para[vh_θh_Idx] .=
            [pf_gens_vh_θh...;]

        # next iteration initial value for pf
        
        red_vh_θh_idq =
            get_red_vh_θh_idq_from_pf_sol(
                pf_sol,
                dyn_pf_δ_ω_ed_dash_eq_dash,
                gens_uh_Q_from_red_sol_para )
                
        if simulation_time_to >= time_final
            push!(result, integrator.sol)
        end

         grantedtime = simulation_time_to

    end

    system_sim_sol = result[1]

    
    #-----------------------------------------------
    # plot
    #-----------------------------------------------

    plot_idxs_a = get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = system_sim_sol,
        node_syms_labels = im_sym,
        bus_name = "bus1",
        vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    """
    #-----------------------------------------------
    # case 2
    #-----------------------------------------------

    #-----------------------------------------------
    #-----------------------------------------------
    # simulation via IterativeCallback 
    #----------------------------------------------
    #----------------------------------------------


    sim_state_x0 =
        state[ im_vars_Idx_in_state ]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)

    #-----------------------------------------------

    stateDiffCache_gens =
        [stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]

    #-----------------------------------------------

    dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =
        [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i])
          for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in
              zip( gens_vh_θh,  ωs_ωref0_vref0_porder0) ]

    gens_size_offset_Idx =
        map((dim) ->
        create_size_offset_Idx(dim; counter = 0 ),
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gen_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------

    dims_gens_nodes_id_iq_pg_vh =
        length.( id_iq_pg_vh )

    _,_, gens_nodes_id_iq_pg_vh_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_id_iq_pg_vh;
        counter = 0)

    #-----------------------------------------------

    dims_gens_nodes_ωs_ωref0_vref0_porder0 =
        length.( ωs_ωref0_vref0_porder0 )

    _,_, gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_ωs_ωref0_vref0_porder0;
        counter = 0)

    #-----------------------------------------------

    dims_gens_nodes_vh_θh =
        length.( gens_vh_θh )

    _,_, gens_nodes_vh_θh_indx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_vh_θh;
        counter = 0)
    
    #-----------------------------------------------

    f_gens_vh_θh =
        [ gens_vh_θh...;]
    
    f_ωs_ωref0_vref0_porder0 =
        [ωs_ωref0_vref0_porder0...;]
    
    dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =
        length.(
            [f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0])

    _,_, vh_θh_ωs_ωref0_vref0_porder0_Idx =
        create_size_offset_Idx(
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0;
            counter = 0)

    vh_θh_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[1]

    ωs_ωref0_vref0_porder0_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[2]
        
    kwd_dyn_vh_θh_ωs_ωref0_vref0_etc =
        (; f_gens_vh_θh,
         f_ωs_ωref0_vref0_porder0 )

    #-----------------------------------------------

    sim_fun_system_ode_flat_agg_para_Idxs =
        (; vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx )

    #-----------------------------------------------
    
    sim_fun_system_kwd_flat_agg_para =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         # ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_dash_Xq_dash_view,         
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection ,
         gens_sim_fun_gen_para_Idxs,
         sim_fun_system_ode_flat_agg_para_Idxs,         
         sys_states_Idxs_and_mat_Idxs,         
         gens_nodes_id_iq_pg_vh_idx_in_Idx,
         gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         gens_nodes_vh_θh_indx_in_Idx )

    #-----------------------------------------------
    
    sim_func! =
        system_flat_agg_ode_im_model_func!

    #------------------------------------------------
    
    sim_fun_para = sim_fun_system_ode_flat_agg_para  =
        vcat(f_gens_vh_θh, f_ωs_ωref0_vref0_porder0 )
    
    # sim_fun_para =
    #     sim_fun_system_ode_flat_agg_para


    #------------------------------------------------
    # powerflow func 
    #------------------------------------------------
    
    pf_fun_mismatch =
        integrated_nll_dyn_pf_iip_ΔPQ_Δidq_mismatch

    pf_func = NonlinearFunction( ( g, x, p ) ->
        pf_fun_mismatch(
            g, x, p;
            kwd_net_param =
                kwd_net_param,
            Idxs_type =
                Idxs_type  ) )

    pf_prob = NonlinearProblem(
        pf_func,
        red_vh_θh_0_idq,
        gen_nodes_δ_ω_ed_dash_eq_dash_views ) 

    # cache = init( pf_prob )


    t_pf_alg = NewtonRaphson(
        autodiff =
            AutoForwardDiff(
                ; chunksize =
                    NonlinearSolve.pickchunksize(
                        pf_prob.u0),
      ),
    )

     cache = init(pf_prob, t_pf_alg  )

     step!(cache, t_pf_alg )

    # NonlinearSolve.solve!(cache, t_pf_alg )

    pf_sol_u = copy( cache.u )


    #-----------------------------------------------
    # dynamic ode func and prob
    #-----------------------------------------------
            
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para =
                sim_fun_system_kwd_flat_agg_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
         sim_fun_para  )

    #------------------------------------------------

    pf_solve_affect_time( integrator ) =
        integrator.t + pf_iterative_Δt

    function fun_ode_pf_solve_affect!(
        integrator;
        pf_param_affect =
            pf_param_affect )
        
        (; cache,
         gens_uh_Q_from_red_sol_para,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx,
         nodes_state_Idx,
         im_nodes_voltage_Idx,
         im_vars_Idx_in_state,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges,
         pf_alg) =
             pf_param_affect

        δ_ω_ed_dash_eq_dash =
            [ integrator.u[idx]
              for idx in
                  nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state]

       # cache.u0 =  pf_red_vh_θh_idq

       cache.p =  δ_ω_ed_dash_eq_dash
       
       NonlinearSolve.solve!( cache, pf_alg )

       pf_sol_u = copy( cache.u )

       (; pf_red_vh_θh_idq,
        pf_gens_flat_vh_θh ) =
            get_pf_gens_flat_vh_θh_and_red_vh_θh_idq_from_pf_sol_u( pf_sol_u, δ_ω_ed_dash_eq_dash, gens_uh_Q_from_red_sol_para )

       integrator.p[ vh_θh_Idx ] .=
           pf_gens_flat_vh_θh

       cache.u0 .=
           pf_red_vh_θh_idq

        """
        nodes_ur_ui =
            get_nodes_ur_ui_from_pf_sol_u(
                pf_sol_u,
                gens_uh_Q_from_red_sol_para )        

        for idx in nodes_u_Idx_in_ranges
            
            integrator.u[idx] .= nodes_ur_ui[idx]
            
        end

        """

        return nothing
                         
    end

    pf_iterative_Δt = 0.001

    pf_param_affect =
        (; cache,
         gens_uh_Q_from_red_sol_para,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx,
         nodes_state_Idx,
         im_nodes_voltage_Idx,
         im_vars_Idx_in_state,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges,
         pf_alg)

    # https://docs.sciml.ai/DiffEqCallbacks/stable/timed_callbacks/

    cb_periodic =  PeriodicCallback(
        (integrator) ->
            fun_ode_pf_solve_affect!(
                integrator;
                pf_param_affect ),
        pf_iterative_Δt )


    cb_iterative =
        IterativeCallback(
            (integrator) ->
                pf_solve_affect_time(
                    integrator.t,
                    pf_iterative_Δt),
            (integrator) ->
                fun_ode_pf_solve_affect!(
                    integrator;
                    pf_param_affect ) )


    #-----------------------------------------------------
    #-----------------------------------------------------
    # simulation
    #-----------------------------------------------------
    #-----------------------------------------------------

    # system_sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, ode_alg, dt=dt, callback = cb_periodic )

    
    #-----------------------------------------------------
    # plot
    #-----------------------------------------------------

plot_idxs_a =
    get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = system_sim_sol,
        node_syms_labels = im_sym,
        bus_name = "bus1",
        vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    #-----------------------------------------------------
    # sensitivity
    #-----------------------------------------------------
    
    sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para =
                    sim_fun_system_kwd_flat_agg_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg = ForwardDiffSensitivity() )

    sim_sol = DifferentialEquations.solve(
        sens_prob , ode_alg, dt=dt  )

    #-----------------------------------------------------
    

    # # Coupling of dynamic simmulation of ode to pf through
    # # dyn_pf_δ_ω_ed_dash_eq_dash

    
    # dyn_pf_δ_ω_ed_dash_eq_dash =
    #     get_im_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
    #         system_sim_sol.u[1,end],
    #         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    # # dyn_pf_δ_ω_ed_dash_eq_dash =
    # #     [ system_sim_sol.u[idx]
    # #       for idx in
    # #           nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state  ]

    # #-----------------------------------------------------

    # pf_sol = NonlinearSolve.solve(pf_prob, pf_alg )

    # #-----------------------------------------------------
    

    # # Coupling of powerflow simmulation to
    # # dynamic simulation through vec_pf_gens_vh_θh.


    # # vec_pf_gens_vh_θh

    # pf_gens_vh_θh =
    #     get_pf_gens_uh_from_pf_sol(
    #         pf_sol, 
    #             gens_uh_Q_from_red_sol_para )

    # red_vh_θh_idq =
    #     get_red_vh_θh_idq_from_pf_sol(
    #         pf_sol,
    #         dyn_pf_δ_ω_ed_dash_eq_dash,
    #         gens_uh_Q_from_red_sol_para )

    
end


function combined_integrated_dyn_and_pf_simulation()

    
    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)

    sim_model_type        = "im-model"
    
    # ode_alg             = Rodas4()
    # algr_name           = "rodas4" 
    
    ode_alg               = ImplicitMidpoint()    
    algr_name             = "ImplicitMidpoint"   
    dt                    = 0.01

    pf_alg                = NewtonRaphson()

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = false

    with_cb  = false

    Idxs_type = :Idxs_hybrid  # :Idxs_im

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
        
    #--------------------------------------
    #--------------------------------------
    
    netd  = NetworkData( dynamics_case()... )
    
    #--------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )


    #----------------------------------------
    # pf init
    #---------------------------------------- 

    """ other options are

            Idxs_type = :Idxs_industrial

            Idxs_type = :Idxs_im
    """
    
    (; nll_dyn_iip_pf_param, ) =
        get_integrated_nll_dyn_iip_pf_param(
            netd;
            Idxs_type = Idxs_type ) # :Idxs_im

    (; red_ΔPQ_Δidq,
     red_vh_θh_0_idq,         
     gens_uh_Q_from_red_sol_para,
     kwd_net_param ) =
         nll_dyn_iip_pf_param

    #------------------------------------ 
    #-------------------------------------

    pf_idx_and_state =
        kwd_net_param.pf_net_param.pf_idx_and_state

    slack_vh =
        pf_idx_and_state.slack_vh
    
    gens_vh =
        pf_idx_and_state.gens_vh
    
    gens_Idx_and_vh =
        pf_idx_and_state.gens_Idx_and_vh
    
    non_slack_gens_Idx_and_vh =
        pf_idx_and_state.non_slack_gens_Idx_and_vh

    #---------------------------------------
    
    # pf_Idx = kwd_net_param.pf_net_param.pf_Idx

    # red_vh_θh_idx = pf_Idx.red_vh_θh_idx
    
    #-------------------------------------

    pf_and_dyn_idx_and_Idx =
        kwd_net_param.pf_net_param.pf_and_dyn_idx_and_Idx
    
    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist

    #-------------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx

    slack_gens_nodes_idx =
        net_comp_type_idx.slack_gens_nodes_idx

    slack_bus_idx =
        slack_gens_nodes_idx

    non_slack_gens_nodes_idx =
        net_comp_type_idx.non_slack_gens_nodes_idx

    gens_idx = gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    # gens_with_loc_load_idx =
    #     net_comp_type_idx.gens_with_loc_load_idx

    gens_with_loc_load_idx =
        loc_load_exist == false ?
        [] : 
        net_comp_type_idx.gens_with_loc_load_idx

    all_nodes_idx =
        net_comp_type_idx.all_nodes_idx

    #--------------------------------------------

    nodes_size = sum(
        [length(gens_nodes_idx),
         length(non_gens_nodes_idx) ])

    #-------------------------------------------

    dict_n2s =
        net_comp_type_idx.dicts_net_to_streamlined_idx
    
    n2s_slack_gens_idx =
        dict_n2s.dict_n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        dict_n2s.dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx =
        dict_n2s.dict_n2s_gens_idx

    n2s_non_gens_idx =
        dict_n2s.dict_n2s_non_gens_idx

    n2s_gens_with_loc_load_idxs =
        dict_n2s.dict_n2s_gens_with_loc_load_idxs
    
    n2s_load_idx =
        dict_n2s.dict_n2s_load_idx
    
    n2s_transmission_idxs =
        dict_n2s.dict_n2s_transmission_idxs
    
    #-------------------------------------------

    if Idxs_type == :Idxs_hybrid
        
        idx_and_Idx =
            pf_and_dyn_idx_and_Idx.hybrid_pf_etc_idx_and_Idx
    elseif Idxs_type == :Idxs_industrial

        idx_and_Idx =
            pf_and_dyn_idx_and_Idx.industrial_model_pf_idx_and_Idx
    elseif Idxs_type == :Idxs_im

        idx_and_Idx =
            pf_and_dyn_idx_and_Idx.im_model_pf_idx_and_Idx
    else
         nothing
    end

    red_vh_θh_idx =
        idx_and_Idx.red_vh_θh_idx
    
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

    #-------------------------------------------

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx

    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash

    #------------------------------------------
    # dyn init
    #------------------------------------------
    
    (; network_bus_names,
       non_gens_bus_names,
       gens_bus_names,
       Loads_bus_names,
       Trans_bus_names) =
           make_case_buses_names(
               netd.nodes )
    
    net_class_names =
        (;network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #------------------------------------------

    (; net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #------------------------------------------
    
    # im_sym, im_mass_matrix =
    #     generate_im_sym_and_mass_matrix(
    #         ; nodes = netd.nodes )

    (; im_model_ode_mass_matrix,
     im_model_pf_mass_matrix) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes = netd.nodes )

    (; gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes = netd.nodes  )

    im_sym =
        gens_nodes_im_vars_labels

    im_mass_matrix =
        im_model_ode_mass_matrix
        
    #-------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    
    #------------------------------------------
    #------------------------------------------
    
    (;
     nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para,
     im_idq_pf_cal )  =
         v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )
    
    #------------------------------------------------
    
    (;
     gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para

    # ------------------------------------------------

    (;
     im_nodes_voltage_Idx,
     im_vars_Idx_in_state, 
     nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     gen_nodes_ra_Xd_dash_Xq_dash_view,
     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
     im_vars_view_in_state,
     gen_nodes_δ_ω_ed_dash_eq_dash_views,
     gens_vh_θh_view,
     nodes_pf_U_view ) =
         im_para_aux_inputs

    (;
     gens_idx,
     non_gens_idx,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges ) =
         im_nodes_voltage_Idx
    
    # -----------------------------------------------

    δ_ω_ed_dash_eq_dash =
        gen_nodes_δ_ω_ed_dash_eq_dash_views
    
    id_iq_pg_vh =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view

    ωs_ωref0_vref0_porder0 =
        gens_nodes_ωs_ωref0_vref0_porder0_view
        
    gens_vh_θh = gens_vh_θh_view
    
    # -----------------------------------------------

    # pf_net_param, _, _ = global_pf_param

    # _, _, _, _, _, pf_and_dyn_idx_and_Idx = pf_net_param

    # _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
    #     pf_param_views

    # _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
    #     sd_pf_views


    # (; each_gens_im_vars_Idx_in_state,
    #  im_vars_Idx_in_state,
    #  im_vars_view_in_state ) =
    #      para_update_gen_Ax_aux
    
    #-----------------------------------------------
    # Stability
    #----------------------------------------------- 
    
    (Ax_Bx_Cx_views,
     Ax_Bx_Cx_matrix,
     states_and_mat_Idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    (vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views)  =
        Ax_Bx_Cx_views
    
    (Ax_matrix,
     Bx_matrix,
     Cx_matrix) =
         Ax_Bx_Cx_matrix

    (nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs)  =
         states_and_mat_Idxs

    sys_states_Idxs_and_mat_Idxs =
        (; nodes_state_Idx,
         Bx_idxs, Cx_idxs,
         id_iq_ph_vh_idxs,
         ω_ref_ωref_v_ref_idxs )
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )
    
    #-------------------------------------------
    #-------------------------------------------

    ## alternative
    
    # nodes_ur_ui = [state[idx]
    #                for idx in
    #                    nodes_u_Idx_in_ranges]
    
    # dyn_nodes_vh_θh = [
    #     abs.(x_from_xr_xi.( nodes_ur_ui ))...;
    #     angle.(x_from_xr_xi.( nodes_ur_ui ))... ]

    # dyn_nodes_vh_θh = [
    #     abs.(x_from_xr_xi.(nodes_pf_U_view))...;
    #     angle.(x_from_xr_xi.(nodes_pf_U_view))... ]

    # dyn_red_vh_θh =
    #     dyn_nodes_vh_θh[ red_vh_θh_idx ]
    
    #----------------------------------------------- 
    #----------------------------------------------- 
    #-----------------------------------------------
    #-----------------------------------------------

    """

    #-----------------------------------------------
    # case 0
    #-----------------------------------------------

    sim_state_x0 =
        state[im_vars_Idx_in_state]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0),
    #         chunk_size )
    
    # stateDiffCache =
    #     DiffCache(
    #         similar(sim_state_x0) )

    #-----------------------------------------------------

    stateDiffCache_gens =
        [ stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]

    #-----------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i]) for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh,  ωs_ωref0_vref0_porder0 ) ]

    gens_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gens_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =  length.( [ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh,  ωs_ωref0_vref0_porder0 ) ] )

    system_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), [dims_gens_vh_θh_ωs_ωref0_vref0_porder0]  )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(system_size_offset_Idx)
    
    sim_fun_system_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------------

    system_ode_para = [[ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh,  ωs_ωref0_vref0_porder0) ]...;]
    
    #-----------------------------------------------------
    
    sim_fun_system_kwd_para =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         # ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_dash_Xq_dash_view,
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection,
         gens_sim_fun_gens_para_Idxs,
         sim_fun_system_para_Idxs,         
         sys_states_Idxs_and_mat_Idxs,
         para_update_gen_Ax_aux)

    #-----------------------------------------------------
    
    sim_func! = system_ode_im_model_func!

    sim_fun_para = system_ode_para
    
    #----------------------------------------------------- 
        
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para = sim_fun_system_kwd_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
        sim_fun_para  )
    
    #-----------------------------------------------------
    # simulate 
    #-----------------------------------------------------

    # system_sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, ode_alg, dt=dt  )

    #-----------------------------------------------------
    # sensitivity
    #-----------------------------------------------------
    
    sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para =
                    sim_fun_system_kwd_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para; sensealg = ForwardDiffSensitivity() )


    sa_sim_sol = DifferentialEquations.solve(
        sens_prob , ode_alg, dt=dt  )

    #-----------------------------------------------------
    # plot
    #-----------------------------------------------------

    plot_idxs_a = get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = system_sim_sol,
        node_syms_labels = im_sym,
        bus_name = "bus1",
        vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    # plot_a = make_plot_of_buses_volt_mag(
    #     ; sol = system_sim_sol,
    #     network_vars_labels = im_sym,
    #     nodes_name = ["bus1", "bus2"],
    #     vars = [:u_r, :u_i],
    #     tspan = sim_timespan,
    #     fmt = :png)
    

    # plot_b = make_plot_of_buses_volt_angle(
    #     ; sol = system_sim_sol,
    #     network_vars_labels = im_sym,
    #     nodes_name = ["bus1", "bus2"],
    #     vars = [:u_r, :u_i],
    #     tspan = sim_timespan,
    #     fmt = :png)

    """
    
    #-----------------------------------------------------
    # case 1
    #-----------------------------------------------------
    
    """
    sim_state_x0 =
        state[ im_vars_Idx_in_state ]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)

    #-----------------------------------------------------

    stateDiffCache_gens =
        [stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]

    #-----------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i]) for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh,  ωs_ωref0_vref0_porder0) ]

    gens_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gen_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------------

    dims_gens_nodes_id_iq_pg_vh =
        length.( id_iq_pg_vh )

    _,_, gens_nodes_id_iq_pg_vh_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_id_iq_pg_vh;
        counter = 0)

    #-----------------------------------------------------

    dims_gens_nodes_ωs_ωref0_vref0_porder0 =
        length.( ωs_ωref0_vref0_porder0 )

    _,_, gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_ωs_ωref0_vref0_porder0;
        counter = 0)

    #-----------------------------------------------------

    dims_gens_nodes_vh_θh =
        length.( gens_vh_θh )

    _,_, gens_nodes_vh_θh_indx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_vh_θh;
        counter = 0)
    
    #-----------------------------------------------------

    f_gens_vh_θh =
        [ gens_vh_θh...;]
    
    f_ωs_ωref0_vref0_porder0 =
        [ωs_ωref0_vref0_porder0...;]
    
    dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =
        length.(
            [f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0])

    _,_, vh_θh_ωs_ωref0_vref0_porder0_Idx =
        create_size_offset_Idx(
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0;
            counter = 0)

    vh_θh_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[1]

    ωs_ωref0_vref0_porder0_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[2]
        
    kwd_dyn_vh_θh_ωs_ωref0_vref0_etc =
        (; f_gens_vh_θh,
         f_ωs_ωref0_vref0_porder0 )

    #-----------------------------------------------------

    sim_fun_system_ode_flat_agg_para_Idxs =
        (; vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx )

    #-----------------------------------------------------
    
    sim_fun_system_kwd_flat_agg_para =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         # ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_dash_Xq_dash_view,         
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection ,
         gens_sim_fun_gen_para_Idxs,
         sim_fun_system_ode_flat_agg_para_Idxs,         
         sys_states_Idxs_and_mat_Idxs,         
         gens_nodes_id_iq_pg_vh_idx_in_Idx,
         gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         gens_nodes_vh_θh_indx_in_Idx )

    #-----------------------------------------------------
    
    sim_func! =
        system_flat_agg_ode_im_model_func!

    #-----------------------------------------------------
    
    sim_fun_para = sim_fun_system_ode_flat_agg_para  =
        vcat( f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0 )
    
    # sim_fun_para =
    #     sim_fun_system_ode_flat_agg_para

    #-----------------------------------------------------
    # dynamic ode func and prob
    #-----------------------------------------------------
            
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para =
                sim_fun_system_kwd_flat_agg_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
         sim_fun_para  )

    #-----------------------------------------------------
    # powerflow func 
    #-----------------------------------------------------
    
    pf_fun_mismatch =
        integrated_nll_dyn_pf_iip_ΔPQ_Δidq_mismatch

    pf_func = NonlinearFunction( ( g, x, p ) ->
        pf_fun_mismatch(
            g, x, p;
            kwd_net_param =
                kwd_net_param,
            Idxs_type =
                Idxs_type  ) )

    #-----------------------------------------------------

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, ode_alg, dt=dt )

    #-----------------------------------------------------
    
    # system_sim_sol = DifferentialEquations.solve(
    #     ODEProblem( ODEFunction{true}(
    #         (dx, x, p, t) ->
    #             sim_func!(
    #                 dx, x, p, t;
    #                 sim_fun_kwd_para =
    #                     sim_fun_system_kwd_flat_agg_para )
    #         ; mass_matrix = im_mass_matrix,
    #         syms = im_sym ),
    #                 sim_state_x0,
    #                 sim_timespan,
    #                 sim_fun_para ),
    #     ode_alg, dt=dt )
      
    #-----------------------------------------------------

   # Coupling betwing ode and pf, via
   # δ_ω_ed_dash_eq_dash system_sim_sol.u

    δ_ω_ed_dash_eq_dash =
        get_im_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            system_sim_sol.u[1,end],
            nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    #------------------------------------------------------

    pf_prob = NonlinearProblem(
        pf_func,
        red_vh_θh_0_idq,
        δ_ω_ed_dash_eq_dash ) 

    pf_sol = NonlinearSolve.solve(
        pf_prob,
        pf_alg )

    #------------------------------------------------------

    # pf_sol = NonlinearSolve.solve(
    #     NonlinearProblem(
    #         NonlinearFunction(
    #             ( g, x, p ) ->
    #                 pf_fun_mismatch(
    #                     g, x, p;
    #                     kwd_net_param =
    #                         kwd_net_param,
    #                     Idxs_type =
    #                         Idxs_type  )),
    #         red_vh_θh_0_idq,
    #         δ_ω_ed_dash_eq_dash ),
    #     pf_alg )

    #-----------------------------------------------

    red_vh_θh_idq =
        get_red_vh_θh_idq_from_pf_sol(
            pf_sol,
            δ_ω_ed_dash_eq_dash,
            gens_uh_Q_from_red_sol_para )
    
    #----------------------------------------------
    #----------------------------------------------
    # Integrator inferface
    #----------------------------------------------
    #----------------------------------------------

    time_final  = 10

    simtime  = time_final # +  dt

    grantedtime = 0

    # simulation_time_to = 0

    reinit!(integ)

    integ =
        init( sim_prob,
             ode_alg,
              dt = dt,
              abstol=1e-12,
              reltol=1e-12,
              tstops = [time_final],
              advance_to_tstop = true)

    dyn_ode_para = integ.p

    # vh_θh_Idx,     
    # ωs_ωref0_vref0_porder0_Idx
    
    #------------------------------------------

    result = []

    while grantedtime <= simtime

        global grantedtime, dyn_pf_δ_ω_ed_dash_eq_dash, pf_sol, red_vh_θh_idq

        global simulation_time_to, updated_pf_prob, pf_gens_vh_θh

        simulation_time_to = grantedtime + dt


        if simulation_time_to < time_final 

            add_tstop!(integ, simulation_time_to)

        else

            add_tstop!(integ, time_final)

        end

        step!(integ, dt, true)

        # Coupling of dynamic simmulation of ode to pf
        # through dyn_pf_δ_ω_ed_dash_eq_dash

        dyn_pf_δ_ω_ed_dash_eq_dash =
            [ integ.u[idx]
              for idx in
                  nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state  ]
        
        #----------------------------------------------
        
         updated_pf_prob = remake( pf_prob,
                       u0 = red_vh_θh_idq,
                       p = dyn_pf_δ_ω_ed_dash_eq_dash )

        pf_sol = NonlinearSolve.solve(
            updated_pf_prob, pf_alg )
        
        #----------------------------------------------

        # Coupling of powerflow simmulation to
        # dynamic simulation through vec_pf_gens_vh_θh.

        # vec_pf_gens_vh_θh

        pf_gens_vh_θh =
            get_pf_gens_uh_from_pf_sol(
                pf_sol, 
                gens_uh_Q_from_red_sol_para )

        # ode para update
        
        dyn_ode_para[vh_θh_Idx] .=
            [pf_gens_vh_θh...;]

        # next iteration initial value for pf
        
        red_vh_θh_idq =
            get_red_vh_θh_idq_from_pf_sol(
                pf_sol,
                dyn_pf_δ_ω_ed_dash_eq_dash,
                gens_uh_Q_from_red_sol_para )
                
        if simulation_time_to >= time_final
            push!(result, integrator.sol)
        end

         grantedtime = simulation_time_to

    end

    system_sim_sol = result[1]

    
    #-----------------------------------------------
    # plot
    #-----------------------------------------------

    plot_idxs_a = get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
        ; sol = system_sim_sol,
        node_syms_labels = im_sym,
        bus_name = "bus1",
        vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    """

    #-----------------------------------------------
    # case 2
    #-----------------------------------------------

    #-----------------------------------------------
    #-----------------------------------------------
    # simulation via IterativeCallback 
    #----------------------------------------------
    #----------------------------------------------


    sim_state_x0 =
        state[ im_vars_Idx_in_state ]

    chunk_size = length(sim_state_x0 )

    stateDiffCache = similar(sim_state_x0)

    #-----------------------------------------------

    stateDiffCache_gens =
        [stateDiffCache[idx]
         for idx in
             nodes_state_Idx ]

    sim_state_x0_gens =
        [ sim_state_x0[idx]
         for idx in
             nodes_state_Idx ]

    #-----------------------------------------------

    dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =
        [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i])
          for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in
              zip( gens_vh_θh,  ωs_ωref0_vref0_porder0) ]

    gens_size_offset_Idx =
        map((dim) ->
        create_size_offset_Idx(dim; counter = 0 ),
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gen_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------

    dims_gens_nodes_id_iq_pg_vh =
        length.( id_iq_pg_vh )

    _,_, gens_nodes_id_iq_pg_vh_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_id_iq_pg_vh;
        counter = 0)

    #-----------------------------------------------

    dims_gens_nodes_ωs_ωref0_vref0_porder0 =
        length.( ωs_ωref0_vref0_porder0 )

    _,_, gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_ωs_ωref0_vref0_porder0;
        counter = 0)

    #-----------------------------------------------

    dims_gens_nodes_vh_θh =
        length.( gens_vh_θh )

    _,_, gens_nodes_vh_θh_indx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_vh_θh;
        counter = 0)
    
    #-----------------------------------------------

    f_gens_vh_θh =
        [ gens_vh_θh...;]
    
    f_ωs_ωref0_vref0_porder0 =
        [ωs_ωref0_vref0_porder0...;]
    
    dims_gens_vh_θh_ωs_ωref0_vref0_porder0 =
        length.(
            [f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0])

    _,_, vh_θh_ωs_ωref0_vref0_porder0_Idx =
        create_size_offset_Idx(
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0;
            counter = 0)

    vh_θh_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[1]

    ωs_ωref0_vref0_porder0_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_Idx[2]
        
    kwd_dyn_vh_θh_ωs_ωref0_vref0_etc =
        (; f_gens_vh_θh,
         f_ωs_ωref0_vref0_porder0 )

    #-----------------------------------------------

    sim_fun_system_ode_flat_agg_para_Idxs =
        (; vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx )

    #-----------------------------------------------
    
    sim_fun_system_kwd_flat_agg_para =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         # ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_dash_Xq_dash_view,         
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection ,
         gens_sim_fun_gen_para_Idxs,
         sim_fun_system_ode_flat_agg_para_Idxs,         
         sys_states_Idxs_and_mat_Idxs,         
         gens_nodes_id_iq_pg_vh_idx_in_Idx,
         gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         gens_nodes_vh_θh_indx_in_Idx )

    #-----------------------------------------------
    
    sim_func! =
        system_flat_agg_ode_im_model_func!

    #------------------------------------------------
    
    sim_fun_para = sim_fun_system_ode_flat_agg_para  =
        vcat( f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0 )
    
    # sim_fun_para =
    #     sim_fun_system_ode_flat_agg_para


    #------------------------------------------------
    # powerflow func 
    #------------------------------------------------
    
    pf_fun_mismatch =
        integrated_nll_dyn_pf_iip_ΔPQ_Δidq_mismatch

    pf_func = NonlinearFunction( ( g, x, p ) ->
        pf_fun_mismatch(
            g, x, p;
            kwd_net_param =
                kwd_net_param,
            Idxs_type =
                Idxs_type  ) )

    pf_prob = NonlinearProblem(
        pf_func,
        red_vh_θh_0_idq,
        gen_nodes_δ_ω_ed_dash_eq_dash_views ) 

    # cache = init( pf_prob )


    t_pf_alg = NewtonRaphson(
        autodiff =
            AutoForwardDiff(
                ; chunksize =
                    NonlinearSolve.pickchunksize(
                        pf_prob.u0), ),)

     cache = init(pf_prob, t_pf_alg  )

     step!(cache, t_pf_alg )

    # NonlinearSolve.solve!(cache, t_pf_alg )

    pf_sol_u = copy( cache.u )


    #-----------------------------------------------
    # dynamic ode func and prob
    #-----------------------------------------------
            
    sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_func!(
            dx, x, p, t;
            sim_fun_kwd_para =
                sim_fun_system_kwd_flat_agg_para
             )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
     sim_prob = ODEProblem(
        sim_ode_func,
        sim_state_x0,
        sim_timespan,
         sim_fun_para  )

    #------------------------------------------------

    pf_solve_affect_time( integrator ) =
        integrator.t + pf_iterative_Δt

    function fun_ode_pf_solve_affect!(
        integrator;
        pf_param_affect =
            pf_param_affect )
        
        (; cache,
         gens_uh_Q_from_red_sol_para,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx,
         nodes_state_Idx,
         im_nodes_voltage_Idx,
         im_vars_Idx_in_state,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges,
         pf_alg) =
             pf_param_affect

        δ_ω_ed_dash_eq_dash =
            [ integrator.u[idx]
              for idx in
                  nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state]

       # cache.u0 =  pf_red_vh_θh_idq

       cache.p =  δ_ω_ed_dash_eq_dash
       
       NonlinearSolve.solve!( cache, pf_alg )

       pf_sol_u = copy( cache.u )

       (; pf_red_vh_θh_idq,
        pf_gens_flat_vh_θh ) =
            get_pf_gens_flat_vh_θh_and_red_vh_θh_idq_from_pf_sol_u(
                pf_sol_u,
                δ_ω_ed_dash_eq_dash,
                gens_uh_Q_from_red_sol_para )

       integrator.p[ vh_θh_Idx ] .=
           pf_gens_flat_vh_θh

       cache.u0 .=
           pf_red_vh_θh_idq

        """
        nodes_ur_ui =
            get_nodes_ur_ui_from_pf_sol_u(
                pf_sol_u,
                gens_uh_Q_from_red_sol_para )        

        for idx in nodes_u_Idx_in_ranges
            
            integrator.u[idx] .= nodes_ur_ui[idx]
            
        end

        """

        return nothing
                         
    end

    pf_iterative_Δt = 0.001

    pf_param_affect =
        (; cache,
         gens_uh_Q_from_red_sol_para,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx,
         nodes_state_Idx,
         im_nodes_voltage_Idx,
         im_vars_Idx_in_state,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges,
         pf_alg)

    # https://docs.sciml.ai/DiffEqCallbacks/stable/timed_callbacks/

    cb_periodic =  PeriodicCallback(
        (integrator) ->
            fun_ode_pf_solve_affect!(
                integrator;
                pf_param_affect ),
        pf_iterative_Δt )


    cb_iterative =
        IterativeCallback(
            (integrator) ->
                pf_solve_affect_time(
                    integrator.t,
                    pf_iterative_Δt),
            (integrator) ->
                fun_ode_pf_solve_affect!(
                    integrator;
                    pf_param_affect ) )


    #-----------------------------------------------------
    #-----------------------------------------------------
    # simulation
    #-----------------------------------------------------
    #-----------------------------------------------------

    # system_sim_sol = DifferentialEquations.solve(
    #    sim_prob, ode_alg  )

    system_sim_sol =
        DifferentialEquations.solve(
            sim_prob,
            ode_alg,
            dt=dt,
            callback = cb_periodic )

    
    #-----------------------------------------------------
    # plot
    #-----------------------------------------------------

    plot_idxs_a =
        get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(
            ;sol =
                system_sim_sol,
            node_syms_labels =
                im_sym,
            bus_name =
                "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ])

    plot_c =
        make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            node_syms_labels = im_sym,
            bus_name = "bus1",
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)

    plot_d =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol = system_sim_sol,
            network_vars_labels = im_sym,
            nodes_name = ["bus1", "bus2", "bus3"],
            vars = [:δ, :ω, :ed_dash, :eq_dash ],
            tspan = sim_timespan,
            fmt = :png)


    #-----------------------------------------------------
    # sensitivity
    #-----------------------------------------------------
    
    sens_prob = ODEForwardSensitivityProblem(
        (dx, x, p, t) ->
            sim_func!(
                dx, x, p, t;
                sim_fun_kwd_para =
                    sim_fun_system_kwd_flat_agg_para ),
        sim_state_x0,
        sim_timespan,
        sim_fun_para;
        sensealg =
            ForwardDiffSensitivity())

    sim_sol = DifferentialEquations.solve(
        sens_prob,
        ode_alg,
        dt = dt  )

    #-----------------------------------------------------
    

    # # Coupling of dynamic simmulation of ode to pf through
    # # dyn_pf_δ_ω_ed_dash_eq_dash

    
    # dyn_pf_δ_ω_ed_dash_eq_dash =
    #     get_im_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
    #         system_sim_sol.u[1,end],
    #         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    # # dyn_pf_δ_ω_ed_dash_eq_dash =
    # #     [ system_sim_sol.u[idx]
    # #       for idx in
    # #           nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state  ]

    # #-----------------------------------------------------

    # pf_sol = NonlinearSolve.solve(pf_prob, pf_alg )

    # #-----------------------------------------------------
    

    # # Coupling of powerflow simmulation to
    # # dynamic simulation through vec_pf_gens_vh_θh.


    # # vec_pf_gens_vh_θh

    # pf_gens_vh_θh =
    #     get_pf_gens_uh_from_pf_sol(
    #         pf_sol, 
    #             gens_uh_Q_from_red_sol_para )

    # red_vh_θh_idq =
    #     get_red_vh_θh_idq_from_pf_sol(
    #         pf_sol,
    #         dyn_pf_δ_ω_ed_dash_eq_dash,
    #         gens_uh_Q_from_red_sol_para )

    
end


#-----------------------------------------------------
#-----------------------------------------------------
# Diagnosis
#-----------------------------------------------------
#-----------------------------------------------------

function dignosis_v2_simulate_a_dynamics_im_model()


    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)
    ode_alg               = Rodas4()
    sim_model_type        = "industrial-model"
    algr_name             = "rodas4"

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = false

    with_cb  = false
    

    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    

    netd  = NetworkData(
        dynamics_case()... )
    
    #-----------------------------------------------
    #-------------------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #-------------------------------------------------
    #-------------------------------------------------
    
    (; network_bus_names,
     non_gens_bus_names,
     gens_bus_names,
     Loads_bus_names,
     Trans_bus_names) =
         make_case_buses_names(
             netd.nodes )

    net_class_names =
        (; network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #------------------------------------------------

    (;net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes =  netd.nodes )

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )

    #------------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes = netd.nodes )

    #------------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels, im_sym )
    
    #------------------------------------------------
    
    
    im_sys_pf_param_views_and_init =
        v2_get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options... )

    #-------------------------------------------------
    #------------------------------------------------

    (; nodes_cb_sw,
     state,
     global_pf_param,
     named_tup_pf_result,
     im_misc_Idx,
     para_update_gen_Ax_aux,
     im_para_aux_inputs,
     im_pf_para,
     im_ωs_ωref_vref_vhθh_idq_etc,
     im_dyn_pf_up_para,
     im_idq_pf_cal)  =
        im_sys_pf_param_views_and_init
    
    # -----------------------------------------------

    pf_net_param, sd_pf_views, _ =
        global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx =
        pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ =
        pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views =
        sd_pf_views

    
    #-------------------------------------------------
    # Stability
    #-------------------------------------------------  
    #-------------------------------------------------

    (;Ax_Bx_Cx_views,
     Ax_Bx_Cx_matrix,
     idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    (;vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views)  =
        Ax_Bx_Cx_views
    
    (;Ax_matrix,
     Bx_matrix,
     Cx_matrix) =
        Ax_Bx_Cx_matrix

    (;nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs) =
        idxs

    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
         gens_nodes_collection )
    
    #------------------------------------------------
    #------------------------------------------------

    ode_fun_para =
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         im_pf_para )
    
    para_model_fun =
        (; vec_Ax_views,
         ode_fun_para )

    #-------------------------------------------------

    sim_state_x0 =
        state
    
    sim_timespan =
        sim_timespan

    #-------------------------------------------------
    
    dyn_global_pf_options =
        dyn_global_pf_options
    
    counter_array = [1]
    
    stateDiffCache =
        similar( sim_state_x0 )
    
    sim_fun_para =
        (; netd,
        nodes_cb_sw,
        global_pf_param,
        counter_array,
        state,
        stateDiffCache,
        dyn_global_pf_options,
        para_update_gen_Ax_aux,
        para_model_fun,
        im_misc_Idx,
        im_para_aux_inputs,
        im_dyn_pf_up_para,
         im_idq_pf_cal )
     
    """
    #-------------------------------------------------
    # dynamics simulation
    #-------------------------------------------------
    #-------------------------------------------------

    sim_func = dynamics_im_model!
    #-------------------------------------------------
    
    #-------------------------------------------------
        
    sim_ode_func! = ODEFunction{true}(
        sim_func; mass_matrix = im_mass_matrix,
        syms = im_sym )
    
    sim_prob = ODEProblem(sim_ode_func!,
                          sim_state_x0,
                          sim_timespan,
                          sim_fun_para )

    if with_cb == true
        
        cb = industrial_model_fun_make_state_callbacks(
            collect( values( netd.nodes)) )

        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg, callback = cb )
        
    else
        sim_sol = DifferentialEquations.solve(
            sim_prob, ode_alg  )

    end
    
    #--------------------------------------------------
    
    return sim_sol, para_net_names_labels_syms 

     """
    
    (;netd,
     nodes_cb_sw,
     global_pf_param,
     counter_array,
     state,
     stateDiffCache,
     dyn_global_pf_options,
     para_update_gen_Ax_aux,
     para_model_fun,
     im_misc_Idx,
     im_para_aux_inputs,
     im_dyn_pf_up_para,
     im_idq_pf_cal) =
         sim_fun_para
    
    #--------------------------------------------

    (;nodes_u_Idx,
     nodes_u_Idx_in_ranges,
     im_vars_Idx_in_state,
     each_gens_im_vars_Idx_in_state,
     gens_nodes_collection) =
         im_misc_Idx

    #--------------------------------------------
    
    (;vec_Ax_views,
     ode_fun_para) =
         para_model_fun

    (;Ax_matrix,
     Bx_matrix,
     Cx_matrix,
     im_pf_para) =
        ode_fun_para
    
    (;gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para
        
    #--------------------------------------------

    (;pf_net_param,
     sd_pf_views,
     mismatch) =
        global_pf_param

    (;working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view) =
        sd_pf_views
    
    #--------------------------------------------    

    im_vars_view_in_state, _, _ =
        para_update_gen_Ax_aux
    
    #--------------------------------------------    
    #--------------------------------------------    
    # Powerflow
    #--------------------------------------------
    #--------------------------------------------
        
    counter = counter_array[1]

    #--------------------------------------------

    x  = state
    
    dx = similar(x)
    
    #--------------------------------------------

    stateDiffCache = get_tmp(
        stateDiffCache, x)

    if counter == 1
        
        stateDiffCache .= state
    end

    counter_array .+= 1

    #--------------------------------------------

    industrial_dyn_powerflow(
    nd,
    stateDiffCache,
    global_pf_param,
    im_dyn_pf_up_para,
    im_idq_pf_cal;
    dyn_global_pf_options... )

    #--------------------------------------------
    # update  W
    #--------------------------------------------

    update_im_dynamic_id_iq_pg_vh!(
        gens_dynamic_id_iq_pg_vh_by_vhθh_view,
        stateDiffCache,
        im_para_aux_inputs )
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    im_vars_view_in_state .=
        get_gen_nodes_im_vars_from_state(
            stateDiffCache,
            im_vars_Idx_in_state )

    update_gens_nodes_im_Ax_system_matrices!(
        vec_Ax_views,
        im_vars_view_in_state,
        each_gens_im_vars_Idx_in_state,
        gens_nodes_collection )

    #--------------------------------------------
    # ode gens pure state
    #--------------------------------------------

    tt = 0.0

     ode_im_model_func!(
         dx, x, ode_fun_para, tt;
         im_vars_Idx_in_state =
             im_vars_Idx_in_state )
     
    #--------------------------------------------
    # non gens terminal voltages
    #--------------------------------------------

    non_gens_im_voltage_terminal_func!(
        dx,
        x,
        im_para_aux_inputs,
        tt)
    
    #--------------------------------------------
    # gens terminal voltages
    #--------------------------------------------

    gens_im_voltage_terminal_func!(
        dx,
        x,
        im_para_aux_inputs,
        tt)

    #--------------------------------------------
    #--------------------------------------------
    
    x0 = sim_fun_para.state
    
    x0_ps_idx =
        sim_fun_para.im_misc_Idx.im_vars_Idx_in_state
    
    x0_gen_ps_idx =
        sim_fun_para.im_misc_Idx.each_gens_im_vars_Idx_in_state

    # t_vec_Ax, t_vec_Bx, t_vec_Cx = Ax_Bx_Cx

    t_vec_Ax =
        sim_fun_para.para_model_fun[1]

     # .ode_fun_para.Ax_matrix
     
    t_Ax_m =
        sim_fun_para.para_model_fun[2][1]

     # .ode_fun_para.Bx_matrix
     
    t_Bx_m =
        sim_fun_para.para_model_fun[2][2]

     # .ode_fun_para.Cx_matrix
     
    t_Cx_m =
        sim_fun_para.para_model_fun[2][3]

     # .ode_fun_para.im_pf_para
     
    im_pf_para =
        sim_fun_para.para_model_fun[2][4]

    id_iq_pg_vh =
        im_pf_para.gens_dynamic_id_iq_pg_vh_by_vhθh_view

    ωs_ω_ref_vref_porder =
        im_pf_para.gens_nodes_ωs_ωref0_vref0_porder0_view

    t_dx = t_Ax_m * x0[ x0_ps_idx ] +
        t_Bx_m * [id_iq_pg_vh...;] +
        t_Cx_m * [ ωs_ω_ref_vref_porder...; ]

    """
    #--------------------------------------------
    # Ode gens pure state
    #--------------------------------------------

    tt = 0.0
    
    # ode_im_model_func!(
    #     dx,
    #     x,
    #     ( im_vars_Idx_in_state, ode_fun_para),
    #     tt)


     ode_im_model_func!(
         dx, x, ode_fun_para, tt;
         im_vars_Idx_in_state =
             im_vars_Idx_in_state )
     
    #--------------------------------------------
    # Another method for all nodes votages
    #--------------------------------------------

    # im_all_nodes_voltage_terminal_func!(
    #     dx,
    #     x,
    #     (nodes_pf_U_view, nodes_u_Idx_in_ranges),
    #     t)
            
    #--------------------------------------------
    # non gens terminal voltages
    #--------------------------------------------

    non_gens_im_voltage_terminal_func!(
        dx, x, im_para_aux_inputs, tt)
    
    #--------------------------------------------
    # gens terminal voltages
    #--------------------------------------------

    gens_im_voltage_terminal_func!(
        dx, x, im_para_aux_inputs, tt)

    """
    return nothing

   
end


function diagnosis_simulate_a_dynamics_industrial_model_with_or_no_controllers(
    ; dynamics_case = nothing,
    only_gen = nothing,
    sim_timespan  = nothing,
    ode_alg = ode_alg,
    dyn_global_pf_options = nothing ,
    sta_global_pf_options = nothing,
    with_cb  = false)

    # --------------------------------------------------
    # --------------------------------------------------
    # Diagnosis
    # -------------------------------------------------- 
    # --------------------------------------------------

    netd  = NetworkData(
        dynamics_case()... )
    
    #-------------------------------------------------    

    non_gens_nodes =
        get_non_gens_nodes(
            netd.nodes )

    gens_nodes_collection =
        get_gens_nodes(
            netd.nodes )

    #--------------------------------------------------
    
    (;network_bus_names,
     non_gens_bus_names,
     gens_bus_names,
     Loads_bus_names,
     Trans_bus_names) =
         make_case_buses_names(
             netd.nodes )
    
    net_class_names =
        (;network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
    
    #--------------------------------------------------
    #--------------------------------------------------

    if only_gen == false

        (;net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_stab_states_label,
         gens_nodes_algebraic_and_states_labels) =
            generate_industrial_model_labels(
                ; nodes =  netd.nodes )

        net_states_and_var_labels =
            (; net_bus_volts_labels,
             gens_nodes_pure_states_labels,
             gens_nodes_stab_states_label,
             gens_nodes_algebraic_and_states_labels )

        #-------------------------------------------------

        (;industrial_model_sym,
         industrial_model_mass_matrix) =
            generate_industrial_model_sym_and_mass_matrix(
                ; nodes = netd.nodes )

        #-------------------------------------------------

        para_net_names_labels_syms =
            (; net_class_names,
             net_states_and_var_labels,
             industrial_model_sym )

        #-------------------------------------------------

        industrial_model_pf_param_views_and_init_with_or_no_controller =
            get_industrial_model_pf_param_views_and_init_with_or_no_controller( netd; sta_global_pf_options..., only_gen=only_gen )

        #-------------------------------------------------

        (;nodes_cb_sw,
         state,
         global_pf_param,
         named_tup_pf_result,
         industrial_model_misc_Idx,
         para_update_gen_Ax_aux,
         industrial_model_para_aux_inputs,
         industrial_model_pf_para,
         industrial_model_ωs_τm_vref_vhθh_idq,
         industrial_model_dyn_pf_up_para,
         industrial_model_idq_pf_cal)  =
             industrial_model_pf_param_views_and_init_with_or_no_controller

        # -------------------------------------------------

        pf_net_param, sd_pf_views, _ =
            global_pf_param

        pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

        _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ = pf_param_views

        _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views = sd_pf_views
        
        #-------------------------------------------
        # Stability
        #------------------------------------------- 

        (;Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
            create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        (;vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views)  =
             Ax_Bx_Cx_views
        
        (;Ax_matrix,
         Bx_matrix,
         Cx_matrix) =
            Ax_Bx_Cx_matrix

        (;nodes_state_Idx,
         Bx_idxs,
         Cx_idxs) = idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views, vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen =
                only_gen,
            vec_Ax_τm_vf_views =
                nothing )

        #--------------------------------------------
        #--------------------------------------------

        ode_fun_para =
            (; Ax_matrix,
             Bx_matrix,
             Cx_matrix,
             industrial_model_pf_para  )

        para_model_fun =
            (; vec_Ax_views,
             ode_fun_para )

    else

        (;net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_stab_states_label,
         gens_nodes_algebraic_and_states_labels) =
            generate_industrial_model_labels(
                ; nodes =  netd.nodes,
                no_control_device =true )

        net_states_and_var_labels =
            (; net_bus_volts_labels,
             gens_nodes_pure_states_labels,
             gens_nodes_stab_states_label,
             gens_nodes_algebraic_and_states_labels )

        #--------------------------------------------

        (;industrial_model_sym,
         industrial_model_mass_matrix) =
            generate_industrial_model_sym_and_mass_matrix(
                ; nodes = netd.nodes,
                no_control_device =true )

        #--------------------------------------------

        para_net_names_labels_syms =
            (; net_class_names,
             net_states_and_var_labels,
             industrial_model_sym )

        #-------------------------------------------

        industrial_model_pf_param_views_and_init_with_or_no_controller =
            get_industrial_model_pf_param_views_and_init_with_or_no_controller( netd; sta_global_pf_options..., only_gen=only_gen )

        #----------------------------------------------

        (;nodes_cb_sw,
         state,
         global_pf_param,
         named_tup_pf_result,
         industrial_model_misc_Idx,
         para_update_gen_Ax_aux,
         industrial_model_para_aux_inputs,
         industrial_model_pf_para,
         industrial_model_ωs_τm_vref_vhθh_idq ,
         industrial_model_dyn_pf_up_para,
         industrial_model_idq_pf_cal)  =
             industrial_model_pf_param_views_and_init_with_or_no_controller

        # ---------------------------------------------

        pf_net_param, sd_pf_views, _ = global_pf_param

        pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

        _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ = pf_param_views

        _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views = sd_pf_views
        
        #---------------------------------------------
        # Stability
        #--------------------------------------------  
        #---------------------------------------------
        
        (;Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         idxs) =
             create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        (;vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         vec_Ax_τm_vf_views) =
             Ax_Bx_Cx_views
        
        (; Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix) =
             Ax_Bx_Cx_matrix
        
        (; nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         τm_vf_idxs) =
             idxs

        init_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection;
            only_gen =
                only_gen,
            vec_Ax_τm_vf_views =
                vec_Ax_τm_vf_views )

        #---------------------------------------------

        # update_gens_τm_vf!(
        #     vec_τm_vf_views, gens_nodes_τm_vf  )

        #--------------------------------------------

        ode_fun_para =
            (; Ax_matrix,
             Bx_matrix,
             Cx_matrix,
             Ax_τm_vf_matrix,
             industrial_model_pf_para  )

        para_model_fun =
            (; vec_Ax_views,
             ode_fun_para )

    end
    
    #------------------------------------------------

    sim_state_x0  = state

    sim_timespan  = sim_timespan

    #-----------------------------------------------

    dyn_global_pf_options = dyn_global_pf_options

    counter_array  = [1]

    stateDiffCache = similar( sim_state_x0 )

    sim_fun_para = (
        ;
        nodes_cb_sw,
        global_pf_param,
        counter_array,
        state,
        stateDiffCache,
        dyn_global_pf_options,
        para_update_gen_Ax_aux,
        para_model_fun,
        industrial_model_misc_Idx,
        industrial_model_para_aux_inputs,        
        industrial_model_dyn_pf_up_para,
        industrial_model_idq_pf_cal,                
        only_gen )

    return Ax_Bx_Cx_views, sim_fun_para 
      
end


function diagnosis_simulate_an_industrial_model_with_or_no_controllers()

    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)
    ode_alg               = Rodas4()
    sim_model_type        = "industrial-model"
    algr_name             = "rodas4"

    dyn_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        with_δ_ed_eq      = true )

    sta_global_pf_options = (
        ; maxiter         = 40,
        ftol              = 1000*eps(),
        xtol              = 1000*eps(),
        init_pf           = true,
        with_δ_ed_eq      = false )

    only_gen = true

    with_cb  = false
    

    dynamics_case = case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer

    Ax_Bx_Cx,  sim_fun_para  = diagnosis_simulate_a_dynamics_industrial_model_with_or_no_controllers(
        ;
        dynamics_case         = dynamics_case,
        only_gen              = only_gen,
        sim_timespan          = sim_timespan,
        ode_alg               = ode_alg,
        dyn_global_pf_options = dyn_global_pf_options,
        sta_global_pf_options = sta_global_pf_options,
        with_cb  = with_cb )


    x0 = sim_fun_para.state
    x0_ps_idx = sim_fun_para.industrial_model_misc_Idx.industrial_model_pure_states_Idx
    x0_gen_ps_idx = sim_fun_para.industrial_model_misc_Idx.industrial_model_each_gen_nodes_pure_states_idx_in_state

    (t_vec_Ax,
     t_vec_Bx,
     t_vec_Cx) = Ax_Bx_Cx

    vec_Ax =
        sim_fun_para.para_model_fun.vec_Ax_views

    Ax_m =
        sim_fun_para.para_model_fun.ode_fun_para.Ax_matrix
    Bx_m =
        sim_fun_para.para_model_fun.ode_fun_para.Bx_matrix
    Cx_m =
        sim_fun_para.para_model_fun.ode_fun_para.Cx_matrix

    id_iq_pg_vh = sim_fun_para.para_model_fun.ode_fun_para.industrial_model_pf_para.gens_dynamic_id_iq_pg_vh_by_vhθh_view
    
    ωs_τm_vref_porder = sim_fun_para.para_model_fun.ode_fun_para.industrial_model_pf_para.gens_nodes_ωs_τm_vref_porder_view

    t_dx = Ax_m * x0[ x0_ps_idx ] + Bx_m *
        [id_iq_pg_vh...;] + Cx_m * [
            ωs_τm_vref_porder...; ]

    if only_gen == true
        
        t_vec_Ax, t_vec_Bx, t_vec_Cx, t_vec_Ax_τm_vf = Ax_Bx_Cx
        
        Ax_τm_vf_m = sim_fun_para.para_model_fun.ode_fun_para.Ax_τm_vf_matrix
        
        vec_τm_vf = sim_fun_para.para_model_fun.ode_fun_para.industrial_model_pf_para.vec_τm_vf_views

        t2_dx = Ax_m * x0[ x0_ps_idx ] + Bx_m *
            [id_iq_pg_vh...;] + Ax_τm_vf_m * [
                vec_τm_vf...;]
    end
    
    
    return nothing
    
    # return sim_fun_para
    
        
end

     
#-----------------------------------------------------
#-----------------------------------------------------
# Sensitivities
#-----------------------------------------------------
#-----------------------------------------------------

function t_ode_one_im_model_func_dxdp(
    ; sim_state_x0 = sim_state_x0,
     sim_fun! = t_ode_one_im_model_func! ,
     sim_fun_para = sim_fun_para,
     mass_matrix = im_mass_matrix_i,
     syms = im_sym_i,    
     sim_timespan = sim_timespan,
     sim_fun_kwd_para = t_sim_fun_kwd_para,
     ode_alg = ode_alg)

    #--------------------------------------------

     sim_state_x0 = sim_state_x0
     sim_fun! = t_ode_one_im_model_func!
     sim_fun_para = sim_fun_para
     mass_matrix = im_mass_matrix_i
     syms = im_sym_i    
     sim_timespan = sim_timespan
     sim_fun_kwd_para = t_sim_fun_kwd_para
     ode_alg = ode_alg

    #--------------------------------------------
    
     sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_fun!(
            dx, x, p, t;
            sim_fun_kwd_para = sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix_i,
         syms = im_sym_i )

     sim_prob = ODEProblem(
         sim_ode_func,
         sim_state_x0,
         sim_timespan,
         sim_fun_para  )

     sim_sol = DifferentialEquations.solve(
         sim_prob, alg, dt=dt  )


     function f(x)

         prob = remake(sim_prob,
                       u0 = x,
                       p = sim_fun_para )

         DifferentialEquations.solve(prob, alg, dt=dt,
               reltol = 1e-6, abstol = 1e-6,
               saveat = 1)[1, :]

     end


     function g(p)

         prob = remake(sim_prob,
                       u0 = xx,
                       p = p )

         DifferentialEquations.solve(prob, alg, dt=dt,
               reltol = 1e-6, abstol = 1e-6,
               saveat = 1)[1, :]

     end

     t_df_dx = ForwardDiff.jacobian(f, xx)
     
    t_df_dp = t_dg_dp = ForwardDiff.jacobian(
        g, sim_fun_para)

     t_dx_dp  = -(svd( t_df_dx )) \ t_df_dp
    
    #--------------------------------------------

    # https://discourse.julialang.org/t/how-to-autodifferentiate-the-results-of-nlsolve/65680/19
    
     # eltype(p).(u0)
    
     tt = 5.0
     
     xx = sim_sol( tt  )
    
     function t_dfdp( p ) 
         # dx = similar( xx  )
         
         dx = zero(xx)         
         
         sim_fun( dx, xx, p, tt;
                  sim_fun_kwd_para =
                      sim_fun_kwd_para )
         dx
     end

     
    function t_dfdx( x )
        
         dx = similar( x  )
         # dx = zero(x)
         sim_fun( dx, x, sim_fun_para, tt;
                  sim_fun_kwd_para =
                      sim_fun_kwd_para )
         dx
     end

     
    tf_df_dx =
        ForwardDiff.jacobian(
            t_dfdx, xx)

    tf_df_dp =
        ForwardDiff.jacobian(
            t_dfdp, [sim_fun_para...;] )

    tt_dx_dp =
        -(svd(tt_df_dx) \ tt_df_dp)

    return (; t_df_dx,
            t_df_dp,
            t_dx_dp,
            tt_df_dx,
            tt_df_dp,
            tt_dx_dp)

end



function s_ode_one_im_model_func_dxdp(
    ; sim_state_x0 = sim_state_x0,
     sim_fun! = s_ode_one_im_model_func! ,
     poi = [poi...] ,
     poi_idx = poi_idx,
     mass_matrix = im_mass_matrix_i,
     syms = im_sym_i,    
     sim_timespan = sim_timespan,
     sim_fun_kwd_para = s_sim_fun_kwd_para ,
     alg = alg)

    # #--------------------------------------------

    #  sim_state_x0 = sim_state_x0 
    #  sim_fun! = s_ode_one_im_model_func! 
    #  poi = s_poi 
    #  poi_idx = poi_idx 
    #  mass_matrix = im_mass_matrix_i 
    #  syms = im_sym_i 
    #  sim_timespan = sim_timespan 
    #  sim_fun_kwd_para = s_sim_fun_kwd_para 
    #  ode_alg =  ode_alg
     
    # #--------------------------------------------
     
     sto = similar( sim_state_x0  )

     sim_ode_func = ODEFunction{true}(
        (dx, x, p, t) -> sim_fun!(
            dx, x, p, t; poi_idx = poi_idx,
            sim_fun_kwd_para = sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix_i,
         syms = im_sym_i )

     sim_prob = ODEProblem(
         sim_ode_func,
         sim_state_x0,
         sim_timespan,
         poi  )

     sim_sol = DifferentialEquations.solve(
         sim_prob, ode_alg, dt=dt  )



     function s_f(x)

         prob = remake(sim_prob,
                       u0 = x,
                       p = poi )

         DifferentialEquations.solve(prob, ode_alg, dt=dt,
               reltol = 1e-6, abstol = 1e-6,
               saveat = 1)[1, :]

     end


     function s_g(p)

         prob = remake(sim_prob,
                       u0 = xx,
                       p = p )

         DifferentialEquations.solve(prob, ode_alg, dt=dt,
               reltol = 1e-6, abstol = 1e-6,
               saveat = 1)[1, :]

     end

     s_df_dx = ForwardDiff.jacobian(s_f, xx)
     
     s_df_dp = dgdp = ForwardDiff.jacobian(s_g, poi)

     s_dx_dp = -(svd(s_df_dx) \ s_df_dp)
     
     #--------------------------------------------
     
     tt = 5.0
     
     xx = sim_sol( tt  )
     
     # p_ = eltype(xx).(poi)
     
     # function s_dfdp( p::Vector{T} )  where T <:AbstractFloat
     
     function s_dfdp( p )
         
         dx = similar( xx  )
         
         # dx = zero(xx)         
         
         sim_fun!( dx, xx, p, tt;
                  poi_idx = poi_idx,
                  sim_fun_kwd_para =
                      sim_fun_kwd_para )
         dx
     end

     
     function s_dfdx( x )
         
         dx = similar( x  )
         
         # dx = zero(x)
         
         sim_fun!( dx, x, s_poi, tt;
                  poi_idx = poi_idx,
                  sim_fun_kwd_para =
                      sim_fun_kwd_para )
         dx
     end

     
     ss_df_dx = ForwardDiff.jacobian(s_dfdx, xx)

     ss_df_dp = ForwardDiff.jacobian(s_dfdp, [poi...;] )


    return (; s_df_dx,
            s_df_dp,
            s_dx_dp,
            ss_df_dx,
            ss_df_dp,
            ss_dx_dp)
     
end



function s_oop_ode_one_im_model_func_dxdp(
    ; sim_state_x0 = sim_state_x0,
     sim_fun! = s_oop_ode_one_im_model_func!,
     poi = [poi...],
     poi_idx = poi_idx,
     mass_matrix = im_mass_matrix_i,
     syms = im_sym_i,    
     sim_timespan = sim_timespan,
     sim_fun_kwd_para = s_sim_fun_kwd_para,
     ode_alg = ode_alg)
     
    # #--------------------------------------------

    #  sim_state_x0 = sim_state_x0 
    #  sim_fun! = s_oop_ode_one_im_model_func! 
    #  poi = s_oop_poi 
    #  poi_idx = poi_idx 
    #  mass_matrix = im_mass_matrix_i 
    #  syms = im_sym_i 
    #  sim_timespan = sim_timespan 
    #  sim_fun_kwd_para = s_sim_fun_kwd_para 
    #  ode_alg =  ode_alg
     
    # #--------------------------------------------
     
     sim_ode_func = ODEFunction{false}(
        (x, p, t) -> sim_fun!(
            x, p, t; poi_idx = poi_idx,
            sim_fun_kwd_para = sim_fun_kwd_para
             )
        ; mass_matrix = im_mass_matrix_i,
         syms = im_sym_i )

     sim_prob = ODEProblem(
         sim_ode_func,
         sim_state_x0,
         sim_timespan,
         poi )

     sim_sol = DifferentialEquations.solve(
         sim_prob, ode_alg, dt=dt  )

     function s_oop_f(x)

         prob = remake(sim_prob,
                       u0 = x,
                       p = poi )

         DifferentialEquations.solve(prob, ode_alg, dt=dt,
               reltol = 1e-6, abstol = 1e-6,
               saveat = 1)[1, :]

     end


     function s_oop_g(p)

         prob = remake(sim_prob,
                       u0 = xx,
                       p = p )

         DifferentialEquations.solve(prob, ode_alg, dt=dt,
               reltol = 1e-6, abstol = 1e-6,
               saveat = 1)[1, :]

     end

     so_df_dx = ForwardDiff.jacobian(s_oop_f, xx)
     
     so_df_dp = dgdp = ForwardDiff.jacobian(s_oop_g, poi)

     so_dx_dp = -(svd(so_df_dx) \ so_df_dp)
     
    # #--------------------------------------------

     tt = 5.0
     
     xx = sim_sol( tt  )
     
     function s_oop_dfdp( p ) 
         
         sim_fun!(  xx, p, tt;
                  poi_idx = poi_idx,
                  sim_fun_kwd_para =
                      sim_fun_kwd_para )
     end

     
     function s_oop_dfdx( x )
         
         sim_fun!( x, poi, tt;
                  poi_idx = poi_idx,
                  sim_fun_kwd_para =
                      sim_fun_kwd_para )
     end

     
     soo_df_dx = ForwardDiff.jacobian(s_oop_dfdx, xx)

     soo_df_dp = ForwardDiff.jacobian(s_oop_dfdp, [poi...;] )

     soo_dx_dp =  -(svd(soo_df_dx) \ soo_df_dp)

    return  (; so_df_dx,
             so_df_dp,
             so_dx_dp,
             soo_df_dx,
             soo_df_dp,
             soo_dx_dp )
     
end

#-----------------------------------------------------

function get_dynamics_ode_model_and_powerflow_sim_parameters(
    netd;
    pf_alg = pf_alg,
    algr_name = algr_name,
    abstol = 1e-14,
    reltol = 1e-14,
    only_gen  = false,
    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial
    Idxs_type = :Idxs_im ,
    
    diagnostics_data = false,
    dynamics_by_gens = false,
    dynamics_by_per_gen = false
    )

  
    #-----------------------------------------------

    # pf_alg  = NewtonRaphson()
    
    # # ode_alg       = Rodas4()
    
    # # algr_name     = "rodas4" 
    
    # ode_alg         = ImplicitMidpoint()
    
    # algr_name       = "ImplicitMidpoint"

    # dt              = 0.01

    # sim_timespan    = (0.0, 10.0)
    
    # abstol          = 1e-14

    # reltol          = 1e-14
        
    # only_gen  = false

    # Idxs_type = :Idxs_im # :Idxs_hybrid, :Idxs_im, :Idxs_industrial

    #-----------------------------------------------
    #-----------------------------------------------
    
    (; Ynet,
     nodes_idx_with_adjacent_nodes_idx ) =
         get_Ynet_and_nodes_idx_wt_adjacent_nodes_idx(
             netd )

    nodes_node_idx_and_incident_edges_other_node_idx =
        nodes_idx_with_adjacent_nodes_idx

    #-----------------------------------------------
    
    nodes_types_idxs = net_nodes_types_idxs =
        get_net_nodes_type_idxs(
            netd)
    
    slack_gens_nodes_idx =
        nodes_types_idxs.slack_gens_nodes_idx
    
    non_slack_gens_nodes_idx =
        nodes_types_idxs.non_slack_gens_nodes_idx
    
    gens_nodes_idx =
        nodes_types_idxs.gens_nodes_idx
    
    non_gens_nodes_idx =
        nodes_types_idxs.non_gens_nodes_idx
    
    gens_with_loc_load_idx =
        nodes_types_idxs.gens_nodes_with_loc_loads_idx
    
    all_nodes_idx = nodes_types_idxs.all_nodes_idx


    gens_nodes_with_loc_loads_idx =
        gens_with_loc_load_idx

    #-----------------------------------------------

    integrated_pf_vars_and_para_idx =
        get_a_model_integrated_pf_vars_and_para_idx(
            netd;
            Idxs_type =
                Idxs_type,
            no_control_device =
                false )

    # nodes_types_idxs =
    #     integrated_pf_vars_and_para_idx.nodes_types_idxs

    # n2s_idxs =
    #     integrated_pf_vars_and_para_idx.n2s_idxs

    full_types_Idxs_etc =
        integrated_pf_vars_and_para_idx.full_types_Idxs_etc


    full_gens_vh_Idxs =
        full_types_Idxs_etc.full_gens_vh_Idxs
    
    full_gens_θh_Idxs  =
        full_types_Idxs_etc.full_gens_θh_Idxs
    
    full_non_gens_nodes_vh_Idxs  =
        full_types_Idxs_etc.full_non_gens_nodes_vh_Idxs
    
    full_non_gens_nodes_θh_Idxs  =
        full_types_Idxs_etc.full_non_gens_nodes_θh_Idxs

    # ----
   
    full_vars_Idxs = 
        (;full_gens_vh_Idxs,
         full_gens_θh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_non_gens_nodes_θh_Idxs)
    
    #----------------------------------------    
    
    intg_types_Idxs_etc =
        integrated_pf_vars_and_para_idx.intg_types_Idxs_etc


    intg_gens_vh_Idxs =
        intg_types_Idxs_etc.intg_gens_vh_Idxs
    
    intg_gens_θh_Idxs =
        intg_types_Idxs_etc.intg_gens_θh_Idxs
    
    intg_non_gens_nodes_vh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_vh_Idxs
    
    intg_non_gens_nodes_θh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_θh_Idxs

    intg_gen_id_Idxs =
        intg_types_Idxs_etc.intg_gen_id_Idxs
    
    intg_gen_iq_Idxs =
        intg_types_Idxs_etc.intg_gen_iq_Idxs
    
    # ----
   
    intg_vars_Idxs = 
        (;
         intg_gens_vh_Idxs,
         intg_gens_θh_Idxs,
         intg_non_gens_nodes_vh_Idxs,
         intg_non_gens_nodes_θh_Idxs,
         intg_gen_id_Idxs,
         intg_gen_iq_Idxs )

    intg_kwd_para =
        (; intg_vars_Idxs, )
    
    #----------------------------------------    
    
    u_ur_ui_Idx_in_state =
        integrated_pf_vars_and_para_idx.u_ur_ui_Idx_in_state

   (; slack_ur_ui_Idx_in_state, 
    non_slack_ur_ui_Idx_in_state, 
    ur_ui_Idx_in_state, 
    nodes_u_Idx,
    nodes_u_Idx_in_ranges ) =
        u_ur_ui_Idx_in_state

    
    #-----------------------------------------------

    (;
     gens_nodes_u_Idx_in_ranges,
     non_gens_nodes_u_Idx_in_ranges ) =
         gens_and_non_gens_u_Idx_in_ranges(
             all_nodes_idx,
             gens_nodes_idx,
             non_gens_nodes_idx,
                 nodes_u_Idx_in_ranges)
        
    #-----------------------------------------------

    n2s_idxs = n2s_streamlined_idx =
        get_dict_n2s_streamlined_idx(netd)
    
    n2s_slack_gens_idx =
        n2s_idxs.n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        n2s_idxs.n2s_non_slack_gens_idx
    
    n2s_gens_idx = n2s_idxs.n2s_gens_idx
    
    n2s_non_gens_idx =
        n2s_idxs.n2s_non_gens_idx
    
    n2s_gens_with_loc_load_idxs =
        n2s_idxs.n2s_gens_with_loc_load_idxs
    
    n2s_all_nodes_idx =
        n2s_idxs.n2s_all_nodes_idx

    #-----------------------------------------------

    dyn_PQ_δ_ω_ed_dash_eq_Idxs  =
         get_a_model_integrated_pf_PQ_δ_ω_ed_dash_eq_dash_Idxs(
             netd)

    dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs =
        dyn_PQ_δ_ω_ed_dash_eq_Idxs.dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs
    
    (;
     P_gens_dyn_pf_no_ll_para_Idxs,
     Q_gens_dyn_pf_wt_ll_para_Idxs,
     P_non_gens_dyn_pf_wt_ll_para_Idxs,
     Q_non_gens_dyn_pf_wt_ll_para_Idxs,
     δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs,
     P_g_loc_load_dyn_pf_wt_ll_para_Idxs,
     Q_g_loc_load_dyn_pf_wt_ll_para_Idxs) =
         dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs

    P_gens_dyn_para_Idxs =
        P_gens_dyn_pf_no_ll_para_Idxs
    
    Q_gens_dyn_para_Idxs =
        Q_gens_dyn_pf_wt_ll_para_Idxs
    
    P_non_gens_dyn_para_Idxs =
        P_non_gens_dyn_pf_wt_ll_para_Idxs
    
    Q_non_gens_dyn_para_Idxs =
        Q_non_gens_dyn_pf_wt_ll_para_Idxs
    
    δ_ed_eq_pf_dyn_para_Idxs =
        δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs
    
    P_g_loc_load_dyn_para_Idxs =
        P_g_loc_load_dyn_pf_wt_ll_para_Idxs
    
    Q_g_loc_load_dyn_para_Idxs =
        P_g_loc_load_dyn_pf_wt_ll_para_Idxs
    
    #-----------------------------------------------
        
    loc_load_exist =
        loc_load_exist_bool( netd )
    
    #--------------------------------------------        
    #  label syms and mass_matrix
    #-----------------------------------------------
    
    net_class_names =
           make_case_buses_names(
               netd.nodes )

    #-----------------------------------------------

    im_net_states_and_var_labels =
         generate_im_model_labels(
             ;nodes =
                 netd.nodes )

    #-----------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes =
                netd.nodes )

    (; im_model_ode_mass_matrix,
     im_model_pf_mass_matrix ) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes =
                 netd.nodes )

    (; gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes =
                 netd.nodes  )

    im_ode_sym =
        gens_nodes_im_vars_labels

    im_ode_mass_matrix =
        im_model_ode_mass_matrix
    
    #-----------------------------------------------

    para_net_names_labels_syms =
        (;
         net_class_names,
         im_net_states_and_var_labels,
         im_ode_sym )

    #-----------------------------------------------    
    #-----------------------------------------------    
    
    integrated_pf_PQ_param =
        get_a_model_integrated_pf_PQ_param(
            netd )
    (;
     P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         integrated_pf_PQ_param

    #-----------------------------------------------

    gens_nodes_ra_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            netd.nodes_param,
            netd.nodes;
            param_list =
                [ :ra, :X_d_dash, :X_q_dash ],
            gens_view_only = true )
    
    #-----------------------------------------------
    
    industrial_indices_and_conversion_dict =
        get_industrial_model_indices_and_conversion_dict(
        netd; no_control_device = false )
    
   im_indices_and_conversion_dict =
         get_im_indices_and_conversion_dict(
             netd  )

    indices_and_conversion_dict =
        im_indices_and_conversion_dict

    im_vars_Idx_in_state =
        indices_and_conversion_dict.im_vars_Idx_in_state
    
    #-----------------------------------------------
    
    models_types_gens_δ_ω_ed_dash_eq_dash_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )

    nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_hybrid

    nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_industrial

    nodes_δ_ω_ed_dash_eq_dash_Idxs_im =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_im
    
    nodes_δ_ω_ed_dash_eq_dash_Idxs =
        nodes_δ_ω_ed_dash_eq_dash_Idxs_im
    
    #-----------------------------------------------

    
    models_types_gens_δ_ed_dash_eq_dash_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ed_dash, :eq_dash] )

    nodes_δ_ed_dash_eq_dash_Idxs_hybrid =
        models_types_gens_δ_ed_dash_eq_dash_Idxs.Idxs_hybrid

    nodes_δ_ed_dash_eq_dash_Idxs_industrial =
        models_types_gens_δ_ed_dash_eq_dash_Idxs.Idxs_industrial

    nodes_δ_ed_dash_eq_dash_Idxs_im =
        models_types_gens_δ_ed_dash_eq_dash_Idxs.Idxs_im
    
    nodes_δ_ed_dash_eq_dash_Idxs =
        nodes_δ_ed_dash_eq_dash_Idxs_im
    
    
    #-----------------------------------------------
    
    gens_nodes_collection =
        get_gens_nodes(
            netd.nodes )

    #-----------------------------------------------
    
    δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        get_gens_nodes_some_state_vars_Idxs_in_flattend(
            gens_nodes_collection
            ; some_state_vars =
                [ :δ, :ω, :ed_dash, :eq_dash ] )

    
    δ_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        get_gens_nodes_some_state_vars_Idxs_in_flattend(
            gens_nodes_collection
            ; some_state_vars =
                [ :δ, :ed_dash, :eq_dash ] )
    
    #-----------------------------------------------
    
    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, states_and_mat_Idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(
            gens_nodes_collection )

    (;vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views )  =
        Ax_Bx_Cx_views
    
    (;Ax_matrix,
     Bx_matrix,
     Cx_matrix ) =
         Ax_Bx_Cx_matrix

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ωs_ω_ref_v_ref_p_order_idxs )  =
         states_and_mat_Idxs
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views ),
        gens_nodes_collection )

    im_plants_system_matrices =
        (;
         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         
         Ax_matrix,
         Bx_matrix,
         Cx_matrix )
    
    #-----------------------------------------------
    #-----------------------------------------------
    
    integrated_pf_and_init =
         get_a_model_integrated_sta_powerflow_and_init(
             netd;
             pf_alg  =
                 pf_alg)
        
    #-----------------------------------------------
    
    nodes_name =
        integrated_pf_and_init.nodes_name

    branches_name =
        integrated_pf_and_init.branches_name

    gens_nodes_idx =
        integrated_pf_and_init.gens_nodes_idx

    #-----------------------------------------------

    named_tup_pf_result =
        integrated_pf_and_init.named_tup_pf_result

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init

    pf_init_dict = named_tup_pf_result.pf_init_dict
    
    #-----------------------------------------------
    
    im_state =
        im_model_init_operationpoint(
            netd, bus_dict_init  )
    
    industrial_state =
        industrial_model_init_operationpoint(
            netd, bus_dict_init; pure = :pure,
            no_control_device = only_gen )

    ext_state =
        external_init_operationpoint(
            netd, bus_dict_init,
            branch_dict_init )

    hybrid_state =
        init_operationpoint(netd, pf_init_dict)

    #-----------------------------------------------
    
    state = im_state
    
    #-----------------------------------------------
    
    im_vars_view_in_state =
        get_im_vars_view_in_state(
            state,
            im_vars_Idx_in_state )

    
    im_vars_in_state =
        state[im_vars_Idx_in_state ]
    
    #-----------------------------------------------

    sim_state_x0 = state
    
    # state[ im_vars_Idx_in_state ]

    
    gens_sim_state_x0 =
        [ sim_state_x0[ idx ]
         for idx in
             nodes_state_Idx ]

    #-----------------------------------------------

    
    # stateDiffCache = DiffCache(similar(sim_state_x0))

    stateDiffCache = similar(sim_state_x0) 
    
    #-----------------------------------------------
    
    vh = named_tup_pf_result.Vm
    
    θh = named_tup_pf_result.Vθ

    #-----------------------------------------------

    full_vh_θh = [vh; θh]

    #-----------------------------------------------

    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]


    gens_vh_θh =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(gens_vh , gens_θh)]
    
    #-----------------------------------------------

    non_gens_vh = vh[non_gens_nodes_idx]

    non_gens_θh = θh[non_gens_nodes_idx]

    non_gens_vh_θh =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(non_gens_vh , non_gens_θh)]

    #-----------------------------------------------
    
    gens_nodes_ωs_ωref0_vref0_porder0 =
        get_gens_nodes_ωs_ωref0_vref0_porder0(
        gens_nodes_collection,
        bus_dict_init )


    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )
    
    #-----------------------------------------------

    gens_nodes_δ_ω_ed_dash_eq_dash =
        get_gen_im_nodes_ω_ed_dash_eq_dash(
            state,
            nodes_δ_ω_ed_dash_eq_dash_Idxs )

    #-----------------------------------------------
    
    gens_dynamic_id_iq_pg_vh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
        gens_vh_θh,
        gens_nodes_δ_ω_ed_dash_eq_dash,
            gens_nodes_ra_Xd_dash_Xq_dash )
    
    #-----------------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_ed_dash =
        third.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_eq_dash =
        fourth.( gens_nodes_δ_ω_ed_dash_eq_dash )    
    
    #-----------------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )        

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
        for (a_vh,
             a_θh,
             a_δ,
             a_ed_dash,
             a_eq_dash,
             a_ra,
             a_X_d_dash,
             a_X_q_dash ) in
            zip( gens_vh,
                 gens_θh,
                 gens_δ,
                 gens_ed_dash,
                 gens_eq_dash,
                 gens_ra,
                 gens_Xd_dash,
                 gens_Xq_dash ) ]

    gens_i_d_0 =
        first.( gens_id_iq )
    
    gens_i_q_0 =
        last.( gens_id_iq )

    #-----------------------------------------------

    _, intg_vh_Idx, intg_θh_Idx, intg_id_Idx, intg_iq_Idx  =
        get_flat_intg_vh_θh_id_iq_and_idxs(
            vh, θh, gens_i_d_0, gens_i_q_0 )
    
    #-----------------------------------------------

    gens_vd = [
        a_ed_dash +
            a_Xq_dash * a_iq -
            a_ra * a_id
        for ( a_ed_dash,
              a_Xq_dash,
              a_ra,
              a_id,
              a_iq ) in
            zip(gens_ed_dash,
                gens_Xq_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]

    gens_vq = [
        a_eq_dash - a_Xd_dash * a_id - a_ra * a_id
        for ( a_eq_dash, a_Xd_dash, a_ra, a_id, a_iq ) in
            zip(gens_eq_dash,
                gens_Xd_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]    
    
    #-----------------------------------------------

    gens_ph = [ a_vd * a_id + a_vq * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    gens_qh = [ a_vq * a_id - a_vd * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]

    #-----------------------------------------------

    # a_id * a_vh * sin(a_δ - a_θh) + a_iq * a_vh * cos(a_δ - a_θh)
    
    gens_ph_by_vh_θh_δ_id_iq = [
        get_a_gen_active_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    # a_id * a_vh * cos(a_δ - a_θh) - a_iq * a_vh * sin(a_δ - a_θh)
    
    gens_qh_by_vh_θh_δ_id_iq = [
        get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]
    
    #-----------------------------------------------
    
    S_gens = [
        a_vh * exp(im * (a_θh - a_δ + π/2)) * ( a_id - im * a_iq )
        for ( a_vh, a_θh, a_δ, a_id, a_iq) in
            zip(gens_vh,
                gens_θh,
                gens_δ,
                gens_i_d_0,
                gens_i_q_0) ]
    
    #-----------------------------------------------
    
    gens_Pei = [
        ed_dash * i_d_0 + eq_dash * i_q_0 +
            (Xq_dash - Xd_dash ) * i_d_0 *  i_q_0
        for (ed_dash,eq_dash,Xd_dash,Xq_dash,i_d_0,i_q_0 ) in
            zip(gens_ed_dash, gens_eq_dash,
                gens_Xd_dash,
                gens_Xq_dash,
                gens_i_d_0, gens_i_q_0 ) ]

    #----------------------------------------    
    #----------------------------------------
    # dyn_pf_para
    #----------------------------------------
    #----------------------------------------    

    dyn_pf_wt_ll_para =
        [gens_ph...;
         gens_qh...;
         P_non_gens...;
         Q_non_gens...;
         # gens_nodes_δ_ω_ed_dash_eq_dash...;
         P_g_loc_load...;
         Q_g_loc_load...]


    dyn_pf_no_ll_para =
        [gens_ph...;
         gens_qh...;
         P_non_gens...;
         Q_non_gens...
         # gens_nodes_δ_ω_ed_dash_eq_dash...
             ]

    if loc_load_exist == true

        dyn_pf_fun_flat_para =
            dyn_pf_wt_ll_para
        
    else

        dyn_pf_fun_flat_para =
            dyn_pf_no_ll_para
        
    end
    
    #-----------------------------------------------

    dyn_pf_fun_kwd_net_para = 
        (; 
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx )

    #-----------------------------------------------    
    
    dyn_pf_fun_kwd_net_idxs  =  
        (;
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx
         )

    #----------------------------------------

    dyn_pf_fun_kwd_n2s_idxs =
        (;
         n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx )
    
    #----------------------------------------
    
    dyn_pf_fun_kwd_nll_para_vars_Idxs = 
        (;
         P_gens_dyn_para_Idxs,
         Q_gens_dyn_para_Idxs,
         P_non_gens_dyn_para_Idxs,
         Q_non_gens_dyn_para_Idxs,
         δ_ed_eq_pf_dyn_para_Idxs
         ) 

    dyn_pf_fun_kwd_wll_para_vars_Idxs = 
        (;
         P_gens_dyn_para_Idxs,
         Q_gens_dyn_para_Idxs,
         P_non_gens_dyn_para_Idxs,
         Q_non_gens_dyn_para_Idxs,
         δ_ed_eq_pf_dyn_para_Idxs,
         P_g_loc_load_dyn_para_Idxs,
         Q_g_loc_load_dyn_para_Idxs                   
         ) 
    
    #----------------------------------------

    pf_fun_kwd_para =
        (;
         loc_load_exist,
         dyn_pf_fun_kwd_nll_para_vars_Idxs,
         dyn_pf_fun_kwd_wll_para_vars_Idxs,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_net_para,
         δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
         gens_nodes_ra_Xd_dash_Xq_dash )

    #-----------------------------------------------    
    
    full_gens_id_iq =
        (;
         gens_i_d_0,
         gens_i_q_0 )

    full_kwd_para =
        (;
         full_vars_Idxs,
         full_gens_id_iq )

    #----------------------------------------
    
    full_dyn_pf_fun_kwd_para =
        (;
         full_kwd_para,
         pf_fun_kwd_para)
    
    #----------------------------------------

    intg_dyn_pf_fun_kwd_para =
        (;
         intg_vars_Idxs,
         pf_fun_kwd_para )
            
    #-----------------------------------------------
    
    intg_vh_θh_id_iq =
        [full_vh_θh; gens_i_d_0; gens_i_q_0]

    #-----------------------------------------------

    dyn_pf_fun_var_values =
        (;
         full_vh_θh,
         intg_vh_θh_id_iq) 

    dyn_pf_fun_para_values =
        (;
         dyn_pf_wt_ll_para,
         dyn_pf_no_ll_para,
         dyn_pf_fun_flat_para )

    dyn_pf_fun_kwd_para =
        (;
         full_dyn_pf_fun_kwd_para,
         intg_dyn_pf_fun_kwd_para )

    #-----------------------------------------------
    # aux dyn_pf_fun_kwd_wll_para_vars_Idxs
    #-----------------------------------------------
 
    Pg_Qg_Png_Qng_Pgll_Qgll_Idxs =
        get_vars_or_paras_Idxs_in_flattend(
            [ P_gens, Q_gens,
              P_non_gens, Q_non_gens,
              P_g_loc_load, Q_g_loc_load]
            ; dims_given = false )

    #-----------------------------------------------

    _, gens_vh_θh_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx(
            gens_vh_θh )

    _, non_gens_vh_θh_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx(
            non_gens_vh_θh )

    ( f_gens_ωs_ωref0_vref0_porder0,
      gens_ωs_ωref0_vref0_porder0_idx_in_Idx ) =
        get_a_flattend_vars_or_paras_and_Idx(
            gens_nodes_ωs_ωref0_vref0_porder0 )

    _, gens_dyn_id_iq_pg_vh_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx(
            gens_dynamic_id_iq_pg_vh )

    _, gens_δ_ω_ed_dash_eq_dash_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx(
            gens_nodes_δ_ω_ed_dash_eq_dash )

    #----------------------------------------    

    """ length( gens_vh_θh[1] )  and
    length( gens_nodes_δ_ω_ed_dash_eq_dash[1] ) are dimensions
    of a gen_vh_θh and  gen_δ_ω_ed_dash_eq_dash respectively (2, 4)
    in flattend gen_vh_θh_δ_ω_ed_dash_eq_dash for
    a_gen_voltage_terminal_func """

    a_gen_vtf_vh_θh_δ_ω_ed_dash_eq_dash_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [length( gens_vh_θh[1] ),
             length( gens_nodes_δ_ω_ed_dash_eq_dash[1] ) ]
            ; dims_given = true )

    a_gen_vtf_vh_θh_Idx, a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx =
        a_gen_vtf_vh_θh_δ_ω_ed_dash_eq_dash_Idx

    #----------------------------------------    


    """ length( gens_nodes_ωs_ωref0_vref0_porder0[1] )  and
    length( gens_dynamic_id_iq_pg_vh[1] ) are dimensions
    of a gen_nodes_ωs_ωref0_vref0_porder0 and
    gen_dynamic_id_iq_pg_vh respectively (4, 4)
    in flattend gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh """

    a_gen_per_var_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [length( gens_dynamic_id_iq_pg_vh[1] ),
             length( gens_dynamic_id_iq_pg_vh[1] )
              ]
            ; dims_given = true )

    per_gen_ωs_ωref0_vref0_porder0_Idx, per_gen_id_iq_pg_vh_Idx =
        a_gen_per_var_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx

    #----------------------------------------    

    # per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_paras_and_Idx
    
    ode_per_para_per_node_and_idxs =
        get_per_vars_or_paras_per_node_flat_para_and_idxs(
             [ gens_nodes_ωs_ωref0_vref0_porder0,
              gens_dynamic_id_iq_pg_vh] )

    # 

    # f_per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     ode_per_para_per_node_and_idxs.flattend_per_vars_or_paras

    # 

    per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_per_para_per_node_and_idxs.per_vars_or_paras_Idx

    per_para_ωs_ωref0_vref0_porder0_Idx, per_para_id_iq_pg_vh_Idx =
        per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx
    # 
    # wrong results
    per_para_per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_per_para_per_node_and_idxs.per_vars_or_paras_per_node_Idx

    #----------------------------------------    

    # f_per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh
    
    ode_per_node_per_vars_and_idxs =
        get_per_node_per_vars_or_paras_flat_para_and_idxs(
             [gens_nodes_ωs_ωref0_vref0_porder0,
              gens_dynamic_id_iq_pg_vh ] )

    # 

    # f_per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     ode_per_node_per_vars_and_idxs.flattend_per_node_vars_or_paras

    # 

    per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_per_node_per_vars_and_idxs.per_node_vars_or_paras_Idx

    # 

    per_node_per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_per_node_per_vars_and_idxs.per_node_per_vars_or_paras_Idxs
    
    #----------------------------------------    
    #----------------------------------------    
   
    # f_gens_ωs_ωref0_vref0_porder0_wt_f_dyn_pf_para_and_Idx

    ode_para_and_dyn_pf_para_and_Idx =
        get_per_vars_or_paras_flat_para_and_idxs(
            [f_gens_ωs_ωref0_vref0_porder0,
             dyn_pf_fun_flat_para ] )

    # 

    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para =
        ode_para_and_dyn_pf_para_and_Idx.flattend_per_vars_or_paras

    # 

    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para_Idx =
        ode_para_and_dyn_pf_para_and_Idx.per_vars_or_paras_Idx


    f_ωs_ωref0_vref0_porder0_Idx, f_dyn_pf_para_Idx =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para_Idx

    #----------------------------------------    
    #----------------------------------------

    # f_gens_ωs_ωref0_vref0_porder0
    # f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para

    #----------------------------------------    
    #----------------------------------------    

    pf_fun_mismatch =
        get_a_model_integrated_dyn_pf_intg_ΔPQ_Δidq_mismatch

    #----------------------------------------    
    #----------------------------------------    

    intra_ode_pf_fun_kwd_para =
        (;
         pf_fun_mismatch,
         intg_dyn_pf_fun_kwd_para,
         gens_nodes_ra_Xd_dash_Xq_dash,

         loc_load_exist,

         dyn_pf_fun_flat_para,
         dyn_pf_fun_kwd_wll_para_vars_Idxs,
         nodes_δ_ω_ed_dash_eq_dash_Idxs,
         nodes_u_Idx_in_ranges )  
    

    #----------------------------------------    
    #----------------------------------------    


    ode_per_gen_models_func_kwd_paras =  [
        (;
         Ax_view,
         Bx_view,
         Cx_view,
         per_gen_ωs_ωref0_vref0_porder0_Idx,
         per_gen_id_iq_pg_vh_Idx )
        for (Ax_view, Bx_view, Cx_view) in
            zip(
              vec_Ax_views,
              vec_Bx_views,
              vec_Cx_views, ) ]
        


    ode_per_para_model_func_kwd_para =
        (;
         Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         per_para_ωs_ωref0_vref0_porder0_Idx,
         per_para_id_iq_pg_vh_Idx )
         
    
    vtf_gens_fun_kwd_para = [
        (;
         gen_ra_Xd_dash_Xq_dash,
         a_gen_vtf_vh_θh_Idx,
         a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx )

        for gen_ra_Xd_dash_Xq_dash in
            gens_nodes_ra_Xd_dash_Xq_dash ]


    vtf_para_and_idxs =
        (;
         vtf_gens_fun_kwd_para,
         gens_nodes_u_Idx_in_ranges,
         non_gens_nodes_u_Idx_in_ranges,
         )

    #----------------------------------------
    #----------------------------------------


    disaggretation_idxs =
        (; 
         gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
         a_gen_vtf_vh_θh_Idx,
         a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx,
         per_gen_ωs_ωref0_vref0_porder0_Idx,
         per_gen_id_iq_pg_vh_Idx,
         per_para_ωs_ωref0_vref0_porder0_Idx,
         per_para_id_iq_pg_vh_Idx,
         f_ωs_ωref0_vref0_porder0_Idx,
         f_dyn_pf_para_Idx,
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs)


    pre_pf_idx_and_para =
        (;
         nodes_u_Idx_in_ranges,
         
         dyn_pf_fun_kwd_wll_para_vars_Idxs,
         )

    # pf_solv_para =
    #     (; pf_fun_mismatch,
    #      intg_dyn_pf_fun_kwd_para,
    #      pf_alg,
    #      abstol,
    #      reltol)


    post_pf_idxs =
        (;
         intg_vh_Idx,
         intg_θh_Idx,
         intg_id_Idx,
         intg_iq_Idx )

    #-----------------------------------------------
    #-----------------------------------------------


    states_and_matrices =
        (;
         gens_sim_state_x0,
         hybrid_state,
         ext_state,
         industrial_state ,
         im_state,

         state,

         im_mass_matrix,
         im_ode_mass_matrix,
         im_model_ode_mass_matrix,
         im_model_pf_mass_matrix  )


    labels_and_symbols =
        (;
         net_class_names,
         im_net_states_and_var_labels,
         gens_nodes_im_vars_labels,
         net_bus_volts_labels,
         im_sym,
         im_ode_sym )

    #-----------------------------------------------
    #-----------------------------------------------

    if dynamics_by_gens == true
        
        dynamics_by_gens_kwd_para =
            (;
             sim_state_x0,
             stateDiffCache,
             nodes_state_Idx,
             im_vars_Idx_in_state,
             loc_load_exist,
             gens_nodes_idx,
             non_gens_nodes_idx,
             gens_nodes_ra_Xd_dash_Xq_dash,
             nodes_δ_ω_ed_dash_eq_dash_Idxs,
             disaggretation_idxs,
             pre_pf_idx_and_para ,
             # pf_solv_para ,
             intra_ode_pf_fun_kwd_para,
             post_pf_idxs,
             vtf_para_and_idxs,
             ode_per_para_model_func_kwd_para,
             gens_nodes_collection,
             vec_Ax_views )
        

        return (;
                dynamics_by_gens_kwd_para,
                f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
                states_and_matrices,
                labels_and_symbols )
        
    elseif dynamics_by_per_gen == true

        
        dynamics_by_per_gen_kwd_para  =
            (;
             sim_state_x0,
             stateDiffCache,
             nodes_state_Idx,
             im_vars_Idx_in_state,
             loc_load_exist,
             gens_nodes_idx,
             non_gens_nodes_idx,
             gens_nodes_ra_Xd_dash_Xq_dash,
             nodes_δ_ω_ed_dash_eq_dash_Idxs,
             disaggretation_idxs,
             pre_pf_idx_and_para ,
             # pf_solv_para ,
             intra_ode_pf_fun_kwd_para,
             post_pf_idxs,
             vtf_para_and_idxs,
             ode_per_gen_models_func_kwd_paras,
             gens_nodes_collection,
             vec_Ax_views )


        return (; dynamics_by_per_gen_kwd_para,
                f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
                states_and_matrices,
                labels_and_symbols )

        
    elseif diagnostics_data == true

        recomposition_Idxs =
            (;
             gens_vh_θh_idx_in_Idx,
             non_gens_vh_θh_idx_in_Idx,
             gens_ωs_ωref0_vref0_porder0_idx_in_Idx, 
             gens_dyn_id_iq_pg_vh_idx_in_Idx,
             gens_δ_ω_ed_dash_eq_dash_idx_in_Idx,
             a_gen_vtf_vh_θh_Idx, 
             a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx, 
             per_gen_ωs_ωref0_vref0_porder0_Idx, 
             per_gen_id_iq_pg_vh_Idx, 
             per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx, 

             per_para_ωs_ωref0_vref0_porder0_Idx, 
             per_para_id_iq_pg_vh_Idx, 
             per_para_per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,
             per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,
             per_node_per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,
             f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para_Idx,
             f_ωs_ωref0_vref0_porder0_Idx, 
             f_dyn_pf_para_Idx  )

            pre_init_data =
                (;
                 Ynet,
                 nodes_idx_with_adjacent_nodes_idx,
                 nodes_node_idx_and_incident_edges_other_node_idx,
                 net_nodes_types_idxs,
                 nodes_types_idxs,
                 gens_nodes_idx,
                 non_gens_nodes_idx,
                 gens_nodes_with_loc_loads_idx,
                 integrated_pf_vars_and_para_idx,
                 full_vars_Idxs,
                 intg_vars_Idxs,
                 u_ur_ui_Idx_in_state,
                 slack_ur_ui_Idx_in_state,
                 non_slack_ur_ui_Idx_in_state,
                 ur_ui_Idx_in_state,
                 nodes_u_Idx,
                 nodes_u_Idx_in_ranges,
                 gens_nodes_u_Idx_in_ranges,
                 non_gens_nodes_u_Idx_in_ranges,
                 n2s_idxs,
                 n2s_streamlined_idx,
                 dyn_PQ_δ_ω_ed_dash_eq_Idxs,
                 loc_load_exist,
                 net_class_names,
                 im_net_states_and_var_labels,
                 im_sym,
                 im_mass_matrix,
                 im_model_ode_mass_matrix,
                 im_model_pf_mass_matrix,
                 gens_nodes_im_vars_labels,
                 net_bus_volts_labels,
                 im_ode_sym,
                 im_ode_mass_matrix,
                 para_net_names_labels_syms,
                 integrated_pf_PQ_param,
                 gens_nodes_ra_Xd_dash_Xq_dash,
                 industrial_indices_and_conversion_dict,
                 im_indices_and_conversion_dict,
                 indices_and_conversion_dict,
                 im_vars_Idx_in_state,
                 models_types_gens_δ_ω_ed_dash_eq_dash_Idxs,
                 nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid,
                 nodes_δ_ω_ed_dash_eq_dash_Idxs_im,
                 nodes_δ_ω_ed_dash_eq_dash_Idxs,
                 gens_nodes_collection,
                 δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
                 nodes_state_Idx,
                 Bx_idxs,
                 Cx_idxs,
                 id_iq_ph_vh_idxs,
                 ωs_ω_ref_v_ref_p_order_idxs,
                 states_and_mat_Idxs,
                 vec_Ax_views,
                 vec_Bx_views,
                 vec_Cx_views,
                 Ax_matrix,
                 Bx_matrix,
                 Cx_matrix )

            #-------------------------------------------
            #-------------------------------------------

            init_data =
                (;
                 integrated_pf_and_init,
                 im_state,
                 industrial_state,
                 ext_state,
                 hybrid_state,
                 state,
                 im_vars_view_in_state,
                 im_vars_in_state,
                 sim_state_x0,
                 gens_sim_state_x0) 

            #--------------------------------------------
            #--------------------------------------------

        post_init_data =
            (;
             vh,
             θh,
             full_vh_θh,
             gens_vh,
             gens_θh,
             gens_vh_θh,
             non_gens_vh,
             non_gens_θh,
             non_gens_vh_θh,
             gens_nodes_ωs_ωref0_vref0_porder0,
             gens_nodes_τm_vf,
             gens_nodes_δ_ω_ed_dash_eq_dash,
             gens_dynamic_id_iq_pg_vh,
             gens_δ, gens_ed_dash,
             gens_eq_dash,
             gens_ra,
             gens_Xd_dash,
             gens_Xq_dash,
             gens_id_iq,
             gens_i_d_0,
             gens_i_q_0,
             gens_vd,
             gens_vq,
             gens_ph,
             gens_qh,
             gens_ph_by_vh_θh_δ_id_iq,
             gens_qh_by_vh_θh_δ_id_iq,
             S_gens,
             gens_Pei,
             pf_fun_kwd_para,
             full_dyn_pf_fun_kwd_para,
             intg_dyn_pf_fun_kwd_para,
             intg_vh_θh_id_iq,
             dyn_pf_fun_var_values,
             dyn_pf_fun_para_values,
             dyn_pf_fun_kwd_para )

            #-------------------------------------------
            #-------------------------------------------

            diagnostics_data =
                (;
                 recomposition_Idxs,
                 pre_init_data ,
                 init_data ,
                 post_init_data  )

            return diagnostics_data
        
    else
                
        nothing 
    end

end



function driver_get_dynamics_ode_model_and_powerflow_sim_parameters( )
    
    #-----------------------------------------------
    #-----------------------------------------------

    pf_alg  = NewtonRaphson()
    
    # ode_alg       = Rodas4()    
    # algr_name     = "rodas4" 
    
    ode_alg         = ImplicitMidpoint()
    
    algr_name       = "ImplicitMidpoint"

    dt              = 0.01

    sim_timespan    = (0.0, 10.0)
    
    abstol          = 1e-14

    reltol          = 1e-14
        
    only_gen  = false

    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial

    Idxs_type = :Idxs_im 

    diagnostics_data = false

    dynamics_by_gens = false

    dynamics_by_per_gen = true    

    #-----------------------------------------------
    #-----------------------------------------------
    
    """
    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plants_v6_P_rscad

    """
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    
    #-----------------------------------------------
    #-----------------------------------------------
    
    netd =
        NetworkData( dynamics_case()... )
    
    #-----------------------------------------------
    #-----------------------------------------------
    
    dyn_ode_model_and_pf_sim_para =
        get_dynamics_ode_model_and_powerflow_sim_parameters(
            netd;            
            pf_alg = pf_alg,
            algr_name = algr_name,            
            abstol = abstol,
            reltol = reltol,
            only_gen  = only_gen,
            Idxs_type = :Idxs_im ,
            diagnostics_data = diagnostics_data,            
            dynamics_by_gens = dynamics_by_gens,            
            dynamics_by_per_gen = dynamics_by_per_gen )

    #-----------------------------------------------------


    """
    states_and_matrices =
        (;
         gens_sim_state_x0,
         hybrid_state,
         ext_state,
         industrial_state ,
         im_state,

         state,

         im_mass_matrix,
         im_ode_mass_matrix,
         im_model_ode_mass_matrix,
         im_model_pf_mass_matrix  )


    labels_and_symbols =
        (;
         net_class_names,
         im_net_states_and_var_labels,
         gens_nodes_im_vars_labels,
         net_bus_volts_labels,
         im_sym,
         im_ode_sym )
    """

    #-----------------------------------------------------
    
     dynamics_by_per_gen_kwd_para =
        dyn_ode_model_and_pf_sim_para.dynamics_by_per_gen_kwd_para

    
    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para =
        dyn_ode_model_and_pf_sim_para.f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para
    states_and_matrices =
        dyn_ode_model_and_pf_sim_para.states_and_matrices

    state =
        states_and_matrices.state

    im_mass_matrix =
        states_and_matrices.im_mass_matrix 

    labels_and_symbols =
        dyn_ode_model_and_pf_sim_para.labels_and_symbols

    im_sym = labels_and_symbols.im_sym
    
    #-----------------------------------------------------
    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------
    #-----------------------------------------------------

    sim_state_x0 = state

    sim_timespan = (0.0 , 1.0 )

    sim_fun_para = f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para

    #-----------------------------------------------------    


    pf_solver =
        (; pf_alg,
         abstol,
         reltol)
         
    
    Pg_Qg_external_control = false
    
    counter_array  = [1]

    sim_fun_kwd_para =
        (;
         sim_state_x0,     
         stateDiffCache,
         pf_solver,
         Pg_Qg_external_control,
         counter_array,
         dynamics_by_per_gen_kwd_para )

    #-------------------------------------------------- 

    # dynamics_by_ode_gens_model_and_powerflow_func!
    # dynamics_by_per_gen_ode_model_and_powerflow_func!

    sim_func  =
        dynamics_by_per_gen_ode_model_and_powerflow_func!

    #-----------------------------------------------

    sim_sol =
        DifferentialEquations.solve(ODEProblem(
            ODEFunction{true}(
                (dx, x, p, t) -> sim_func(
                    dx, x, p, t; kwd_para =
                        sim_fun_kwd_para );
                mass_matrix = im_mass_matrix,
                syms = im_sym ),
            sim_state_x0,
            sim_timespan,
            sim_fun_para ), ode_alg, dt=dt )

    #------------------------------------------------    

    
    # sim_ode_func! = ODEFunction{true}( (dx, x, p, t) ->
    #     sim_func(dx, x, p,  t;  kwd_para =
    #     sim_fun_kwd_para  );
    #     mass_matrix = im_mass_matrix,
            
    #     syms = im_sym
    #          )
    
    # sim_prob = ODEProblem(
    #     sim_ode_func!,
    #     sim_state_x0,
    #     sim_timespan,
    #     sim_fun_para )

    # sim_sol =
    #     DifferentialEquations.solve(sim_prob, ode_alg, dt=dt )

    #----------------------------------------
    #########################################
    #----------------------------------------    


    sim_func  =
        dynamics_by_ode_gens_model_and_powerflow_func!

    #------------------------------------------------    

    sim_state_x0 = state

    sim_timespan = (0.0 , 1.0 )

     sim_fun_para =
         f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para

    #-------------------------------------------------    

    counter_array  = [1]

    sim_fun_kwd_para =
        (;
         counter_array,

         dynamics_by_gens_kwd_para )

    #-------------------------------------------------- 

    sim_ode_func! = ODEFunction{true}(
        ( dx, x, p, t) ->
        sim_func(
            dx, x, p, t ;  kwd_para  =
                sim_fun_kwd_para )
        ; mass_matrix = im_mass_matrix,
        syms = im_sym  )

    
    sim_prob = ODEProblem(
        sim_ode_func!,
        sim_state_x0,
        sim_timespan,
        sim_fun_para )

    sim_sol =
        DifferentialEquations.solve(
            sim_prob, ode_alg, dt=dt  )
    
    
end


function simulate_dynamics_ode_model_and_powerflow_func!( )

    pf_alg  = NewtonRaphson()
    
    # ode_alg       = Rodas4()
    
    # algr_name     = "rodas4" 
    
    ode_alg         = ImplicitMidpoint()
    
    algr_name       = "ImplicitMidpoint"

    dt              = 0.01

    sim_timespan    = (0.0, 10.0)
    
    abstol          = 1e-14

    reltol          = 1e-14
        
    only_gen  = false

    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial

    Idxs_type = :Idxs_im 


    diagnostics_data = false
    
    dynamics_by_gens = false
    
    dynamics_by_per_gen = true
    
    #-----------------------------------------------
    
    """
    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plants_v6_P_rscad

    """
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer

    #-----------------------------------------------
    
    netd =
        NetworkData( dynamics_case()... )
    
    #-----------------------------------------------
    #-----------------------------------------------

    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial
    
    dynamics_ode_model_and_powerflow_sim_parameters =
        get_dynamics_ode_model_and_powerflow_sim_parameters(
            netd;
            pf_alg = pf_alg, algr_name = algr_name,
            abstol = abstol, reltol = reltol,
            only_gen  = only_gen, Idxs_type = :Idxs_im ,
            diagnostics_data = diagnostics_data,
            dynamics_by_gens = dynamics_by_gens ,
            dynamics_by_per_gen = dynamics_by_per_gen
        )
    
    #-------------------------------------------
    #-----------------------------------------------
    # dynamics simulation
    #------------------------------------------------
    #---------------------------------------------

    # dynamics_by_ode_gens_model_and_powerflow_func!
    # dynamics_by_per_gen_ode_model_and_powerflow_func!

    sim_func  =
        dynamics_by_per_gen_ode_model_and_powerflow_func!

    #-------------------------------------------------

    sim_state_x0 = state

    sim_timespan = (0.0 , 1.0 )

     sim_fun_para =
         f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para

    #-------------------------------------------------    

    counter_array  = [1]


    sim_fun_kwd_para =
        (;
         counter_array,

         dynamics_by_per_gen_kwd_para )

    #-------------------------------------------------- 
    
     sim_ode_func! = ODEFunction{true}(
         (dx, x, p, t) ->
        sim_func(dx, x, p,  t;  kwd_para =
        sim_fun_kwd_para  );
        mass_matrix =
            im_sym,
        syms =
            im_mass_matrix )
    
    sim_prob = ODEProblem(
        sim_ode_func!,
        sim_state_x0,
        sim_timespan,
        sim_fun_para )

    sim_sol =
        DifferentialEquations.solve(
            sim_prob, ode_alg, dt=dt )

    #----------------------------------------
    #########################################
    #----------------------------------------    


    sim_func  =
        dynamics_by_ode_gens_model_and_powerflow_func!

    #-------------------------------------------------

    sim_state_x0 = state

    sim_timespan = (0.0 , 1.0 )

     sim_fun_para =
         f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para

    #--------------------------------------------------

    counter_array  = [1]

    sim_fun_kwd_para =
        (;
         counter_array,

         dynamics_by_gens_kwd_para )

    #-------------------------------------------------

    sim_ode_func! = ODEFunction{true}( ( dx, x, p, t) ->
        sim_func( dx, x, p,  t;
                 kwd_para  = sim_fun_kwd_para );
        mass_matrix =
            im_sym,
        syms =
            im_mass_matrix )
    
    sim_prob = ODEProblem(
        sim_ode_func!,
        sim_state_x0,
        sim_timespan,
        sim_fun_para )

    sim_sol =
        DifferentialEquations.solve(
            sim_prob, ode_alg )
    
end


function driver_and_diagnostics_for_dynamics_ode_model_and_powerflow_func!( )
  
    #-----------------------------------------------

    pf_alg  = NewtonRaphson()
    
    # ode_alg       = Rodas4()
    
    # algr_name     = "rodas4" 
    
    ode_alg         = ImplicitMidpoint()
    
    algr_name       = "ImplicitMidpoint"

    dt              = 0.01

    sim_timespan    = (0.0, 10.0)
    
    abstol          = 1e-14

    reltol          = 1e-14
        
    only_gen  = false

    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial
    
    Idxs_type = :Idxs_im 

    #-----------------------------------------------
    #-----------------------------------------------

    diagnostics_data = false
    
    dynamics_by_gens = false
    
    dynamics_by_per_gen = true
    
    #-----------------------------------------------
    #-----------------------------------------------
    
    """
    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plants_v6_P_rscad

    """
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer

    #-----------------------------------------------
    
    netd =
        NetworkData( dynamics_case()... )

    
    #-----------------------------------------------
    #-----------------------------------------------
    
    (; Ynet,
     nodes_idx_with_adjacent_nodes_idx ) =
         get_Ynet_and_nodes_idx_wt_adjacent_nodes_idx(
             netd )

    nodes_node_idx_and_incident_edges_other_node_idx =
        nodes_idx_with_adjacent_nodes_idx

    #-----------------------------------------------
    
    nodes_types_idxs = net_nodes_types_idxs =
        get_net_nodes_type_idxs(
            netd)
    
    slack_gens_nodes_idx =
        nodes_types_idxs.slack_gens_nodes_idx
    
    non_slack_gens_nodes_idx =
        nodes_types_idxs.non_slack_gens_nodes_idx
    
    gens_nodes_idx =
        nodes_types_idxs.gens_nodes_idx
    
    non_gens_nodes_idx =
        nodes_types_idxs.non_gens_nodes_idx
    
    gens_with_loc_load_idx =
        nodes_types_idxs.gens_nodes_with_loc_loads_idx
    
    all_nodes_idx = nodes_types_idxs.all_nodes_idx


    gens_nodes_with_loc_loads_idx =
        gens_with_loc_load_idx

    #-----------------------------------------------

    integrated_pf_vars_and_para_idx =
        get_a_model_integrated_pf_vars_and_para_idx(
            netd;
            Idxs_type =
                Idxs_type,
            no_control_device =
                false )

    # nodes_types_idxs =
    #     integrated_pf_vars_and_para_idx.nodes_types_idxs

    # n2s_idxs =
    #     integrated_pf_vars_and_para_idx.n2s_idxs

    full_types_Idxs_etc =
        integrated_pf_vars_and_para_idx.full_types_Idxs_etc


    full_gens_vh_Idxs =
        full_types_Idxs_etc.full_gens_vh_Idxs
    
    full_gens_θh_Idxs  =
        full_types_Idxs_etc.full_gens_θh_Idxs
    
    full_non_gens_nodes_vh_Idxs  =
        full_types_Idxs_etc.full_non_gens_nodes_vh_Idxs
    
    full_non_gens_nodes_θh_Idxs  =
        full_types_Idxs_etc.full_non_gens_nodes_θh_Idxs

    # ----
   
    full_vars_Idxs = 
        (;full_gens_vh_Idxs,
         full_gens_θh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_non_gens_nodes_θh_Idxs)
    
    #----------------------------------------    
    
    intg_types_Idxs_etc =
        integrated_pf_vars_and_para_idx.intg_types_Idxs_etc


    intg_gens_vh_Idxs =
        intg_types_Idxs_etc.intg_gens_vh_Idxs
    
    intg_gens_θh_Idxs =
        intg_types_Idxs_etc.intg_gens_θh_Idxs
    
    intg_non_gens_nodes_vh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_vh_Idxs
    
    intg_non_gens_nodes_θh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_θh_Idxs

    intg_gen_id_Idxs =
        intg_types_Idxs_etc.intg_gen_id_Idxs
    
    intg_gen_iq_Idxs =
        intg_types_Idxs_etc.intg_gen_iq_Idxs
    
    # ----
   
    intg_vars_Idxs = 
        (;
         intg_gens_vh_Idxs,
         intg_gens_θh_Idxs,
         intg_non_gens_nodes_vh_Idxs,
         intg_non_gens_nodes_θh_Idxs,
         intg_gen_id_Idxs,
         intg_gen_iq_Idxs )

    intg_kwd_para =
        (; intg_vars_Idxs, )
    
    #----------------------------------------    
    
    u_ur_ui_Idx_in_state =
        integrated_pf_vars_and_para_idx.u_ur_ui_Idx_in_state

   (; slack_ur_ui_Idx_in_state, 
    non_slack_ur_ui_Idx_in_state, 
    ur_ui_Idx_in_state, 
    nodes_u_Idx,
    nodes_u_Idx_in_ranges ) =
        u_ur_ui_Idx_in_state

    
    #-----------------------------------------------

    (;
     gens_nodes_u_Idx_in_ranges,
     non_gens_nodes_u_Idx_in_ranges ) =
         gens_and_non_gens_u_Idx_in_ranges(
             all_nodes_idx,
             gens_nodes_idx,
             non_gens_nodes_idx,
                 nodes_u_Idx_in_ranges)
        
    #-----------------------------------------------

    n2s_idxs = n2s_streamlined_idx =
        get_dict_n2s_streamlined_idx(netd)
    
    n2s_slack_gens_idx =
        n2s_idxs.n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        n2s_idxs.n2s_non_slack_gens_idx
    
    n2s_gens_idx = n2s_idxs.n2s_gens_idx
    
    n2s_non_gens_idx =
        n2s_idxs.n2s_non_gens_idx
    
    n2s_gens_with_loc_load_idxs =
        n2s_idxs.n2s_gens_with_loc_load_idxs
    
    n2s_all_nodes_idx =
        n2s_idxs.n2s_all_nodes_idx

    #-----------------------------------------------

    dyn_PQ_δ_ω_ed_dash_eq_Idxs  =
         get_a_model_integrated_pf_PQ_δ_ω_ed_dash_eq_dash_Idxs(
             netd)

    dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs =
        dyn_PQ_δ_ω_ed_dash_eq_Idxs.dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs
    
    (;
     P_gens_dyn_pf_no_ll_para_Idxs,
     Q_gens_dyn_pf_wt_ll_para_Idxs,
     P_non_gens_dyn_pf_wt_ll_para_Idxs,
     Q_non_gens_dyn_pf_wt_ll_para_Idxs,
     δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs,
     P_g_loc_load_dyn_pf_wt_ll_para_Idxs,
     Q_g_loc_load_dyn_pf_wt_ll_para_Idxs) =
         dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs

    P_gens_dyn_para_Idxs =
        P_gens_dyn_pf_no_ll_para_Idxs
    
    Q_gens_dyn_para_Idxs =
        Q_gens_dyn_pf_wt_ll_para_Idxs
    
    P_non_gens_dyn_para_Idxs =
        P_non_gens_dyn_pf_wt_ll_para_Idxs
    
    Q_non_gens_dyn_para_Idxs =
        Q_non_gens_dyn_pf_wt_ll_para_Idxs
    
    δ_ed_eq_pf_dyn_para_Idxs =
        δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs
    
    P_g_loc_load_dyn_para_Idxs =
        P_g_loc_load_dyn_pf_wt_ll_para_Idxs
    
    Q_g_loc_load_dyn_para_Idxs =
        P_g_loc_load_dyn_pf_wt_ll_para_Idxs
    
    #-----------------------------------------------
        
    loc_load_exist =
        loc_load_exist_bool( netd )
    
    #----------------------------------------------- 
    #  label syms and mass_matrix
    #-----------------------------------------------
    
    net_class_names =
           make_case_buses_names(
               netd.nodes )

    #-----------------------------------------------

    im_net_states_and_var_labels =
         generate_im_model_labels(
             ;nodes =
                 netd.nodes )

    #-----------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes =
                netd.nodes )

    (; im_model_ode_mass_matrix,
     im_model_pf_mass_matrix ) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes =
                 netd.nodes )

    (; gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes =
                 netd.nodes  )

    im_ode_sym =
        gens_nodes_im_vars_labels

    im_ode_mass_matrix =
        im_model_ode_mass_matrix
        
    #-----------------------------------------------

    para_net_names_labels_syms =
        (;
         net_class_names,
         im_net_states_and_var_labels,
         im_ode_sym )

    #-----------------------------------------------    
    #-----------------------------------------------    
    
    integrated_pf_PQ_param =
        get_a_model_integrated_pf_PQ_param(
            netd )
    (;
     P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         integrated_pf_PQ_param

    #-----------------------------------------------

    gens_nodes_ra_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            netd.nodes_param,
            netd.nodes;
            param_list =
                [ :ra, :X_d_dash, :X_q_dash ],
            gens_view_only = true )
    
    #-----------------------------------------------

    industrial_indices_and_conversion_dict =
        get_industrial_model_indices_and_conversion_dict(
        netd; no_control_device = false )
    
   im_indices_and_conversion_dict =
         get_im_indices_and_conversion_dict(
             netd  )

    indices_and_conversion_dict =
        im_indices_and_conversion_dict

    im_vars_Idx_in_state =
        indices_and_conversion_dict.im_vars_Idx_in_state
    
    #-----------------------------------------------
    
    models_types_gens_δ_ω_ed_dash_eq_dash_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )

    nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_hybrid

    nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_industrial

    nodes_δ_ω_ed_dash_eq_dash_Idxs_im =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_im
    
    nodes_δ_ω_ed_dash_eq_dash_Idxs =
        nodes_δ_ω_ed_dash_eq_dash_Idxs_im
    
    #-----------------------------------------------
    
    gens_nodes_collection =
        get_gens_nodes(
            netd.nodes )

    #-----------------------------------------------
    
    δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        get_gens_nodes_some_state_vars_Idxs_in_flattend(
            gens_nodes_collection
            ; some_state_vars =
                [ :δ, :ω, :ed_dash, :eq_dash ] )
    
    #-----------------------------------------------
    
    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, states_and_mat_Idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(
            gens_nodes_collection )

    (;vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views )  =
        Ax_Bx_Cx_views
    
    (;Ax_matrix,
     Bx_matrix,
     Cx_matrix ) =
         Ax_Bx_Cx_matrix

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ωs_ω_ref_v_ref_p_order_idxs )  =
         states_and_mat_Idxs
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views ),
        gens_nodes_collection )

    im_plants_system_matrices =
        (;
         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         
         Ax_matrix,
         Bx_matrix,
         Cx_matrix )
    
    #-----------------------------------------------
    #-----------------------------------------------
    
    integrated_pf_and_init =
         get_a_model_integrated_sta_powerflow_and_init(
             netd;
             pf_alg  =
                 pf_alg)
        
    #-----------------------------------------------
    
    nodes_name =
        integrated_pf_and_init.nodes_name

    branches_name =
        integrated_pf_and_init.branches_name

    gens_nodes_idx =
        integrated_pf_and_init.gens_nodes_idx

    #-----------------------------------------------

    named_tup_pf_result =
        integrated_pf_and_init.named_tup_pf_result

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init

    pf_init_dict = named_tup_pf_result.pf_init_dict
    
    #-----------------------------------------------
    
    im_state =
        im_model_init_operationpoint(
            netd, bus_dict_init  )
    
    industrial_state =
        industrial_model_init_operationpoint(
            netd, bus_dict_init; pure = :pure,
            no_control_device = only_gen )

    ext_state =
        external_init_operationpoint(
            netd, bus_dict_init,
            branch_dict_init )

    hybrid_state =
        init_operationpoint(netd, pf_init_dict)
    
    #-----------------------------------------------
    
    state = im_state
    
    #-----------------------------------------------
    
    im_vars_view_in_state =
        get_im_vars_view_in_state(
            state,
            im_vars_Idx_in_state )

    
    im_vars_in_state =
        state[im_vars_Idx_in_state ]
    
    #-----------------------------------------------

    sim_state_x0 = state
    
    # state[ im_vars_Idx_in_state ]

    
    gens_sim_state_x0 =
        [ sim_state_x0[ idx ]
         for idx in
             nodes_state_Idx ]

    #-----------------------------------------------

    
    stateDiffCache = DiffCache(similar(sim_state_x0))
    
    #-----------------------------------------------
    
    vh = named_tup_pf_result.Vm
    
    θh = named_tup_pf_result.Vθ

    #-----------------------------------------------

    full_vh_θh = [vh; θh]

    #-----------------------------------------------

    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]


    gens_vh_θh =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(gens_vh , gens_θh)]
    
    #-----------------------------------------------

    non_gens_vh = vh[non_gens_nodes_idx]

    non_gens_θh = θh[non_gens_nodes_idx]

    non_gens_vh_θh =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(non_gens_vh , non_gens_θh)]

    #-----------------------------------------------
    
    gens_nodes_ωs_ωref0_vref0_porder0 =
        get_gens_nodes_ωs_ωref0_vref0_porder0(
        gens_nodes_collection,
        bus_dict_init )


    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )
    
    #-----------------------------------------------

    gens_nodes_δ_ω_ed_dash_eq_dash =
        get_gen_im_nodes_ω_ed_dash_eq_dash(
            state,
            nodes_δ_ω_ed_dash_eq_dash_Idxs )

    #-----------------------------------------------
    
    gens_dynamic_id_iq_pg_vh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
        gens_vh_θh,
        gens_nodes_δ_ω_ed_dash_eq_dash,
            gens_nodes_ra_Xd_dash_Xq_dash )
    
    #-----------------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_ed_dash =
        third.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_eq_dash =
        fourth.( gens_nodes_δ_ω_ed_dash_eq_dash )    
    
    #-----------------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )        

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
        for (a_vh,
             a_θh,
             a_δ,
             a_ed_dash,
             a_eq_dash,
             a_ra,
             a_X_d_dash,
             a_X_q_dash ) in
            zip( gens_vh,
                 gens_θh,
                 gens_δ,
                 gens_ed_dash,
                 gens_eq_dash,
                 gens_ra,
                 gens_Xd_dash,
                 gens_Xq_dash ) ]

    gens_i_d_0 =
        first.( gens_id_iq )
    
    gens_i_q_0 =
        last.( gens_id_iq )

    #-----------------------------------------------

    _, intg_vh_Idx, intg_θh_Idx, intg_id_Idx, intg_iq_Idx  =
        get_flat_intg_vh_θh_id_iq_and_idxs(
            vh, θh, gens_i_d_0, gens_i_q_0 )
    
    #-----------------------------------------------

    gens_vd = [
        a_ed_dash +
            a_Xq_dash * a_iq -
            a_ra * a_id
        for ( a_ed_dash,
              a_Xq_dash,
              a_ra,
              a_id,
              a_iq ) in
            zip(gens_ed_dash,
                gens_Xq_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]

    gens_vq = [
        a_eq_dash - a_Xd_dash * a_id - a_ra * a_id
        for ( a_eq_dash, a_Xd_dash, a_ra, a_id, a_iq ) in
            zip(gens_eq_dash,
                gens_Xd_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]    
    
    #-----------------------------------------------

    gens_ph = [ a_vd * a_id + a_vq * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    gens_qh = [ a_vq * a_id - a_vd * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]

    #-----------------------------------------------

    # a_id * a_vh * sin(a_δ - a_θh) + a_iq * a_vh * cos(a_δ - a_θh)
    
    gens_ph_by_vh_θh_δ_id_iq = [
        get_a_gen_active_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    # a_id * a_vh * cos(a_δ - a_θh) - a_iq * a_vh * sin(a_δ - a_θh)
    
    gens_qh_by_vh_θh_δ_id_iq = [
        get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]
    
    #-----------------------------------------------
    
    S_gens = [
        a_vh * exp(im * (a_θh - a_δ + π/2)) * ( a_id - im * a_iq )
        for ( a_vh, a_θh, a_δ, a_id, a_iq) in
            zip(gens_vh,
                gens_θh,
                gens_δ,
                gens_i_d_0,
                gens_i_q_0) ]
    
    #-----------------------------------------------
    
    gens_Pei = [
        ed_dash * i_d_0 + eq_dash * i_q_0 +
            (Xq_dash - Xd_dash ) * i_d_0 *  i_q_0
        for (ed_dash,eq_dash,Xd_dash,Xq_dash,i_d_0,i_q_0 ) in
            zip(gens_ed_dash, gens_eq_dash,
                gens_Xd_dash,
                gens_Xq_dash,
                gens_i_d_0, gens_i_q_0 ) ]

    #----------------------------------------    
    #----------------------------------------
    # dyn_pf_para
    #----------------------------------------
    #----------------------------------------    

    dyn_pf_wt_ll_para =
        [gens_ph...;
         gens_qh...;
         P_non_gens...;
         Q_non_gens...;
         gens_nodes_δ_ω_ed_dash_eq_dash...;
         P_g_loc_load...;
         Q_g_loc_load...]


    dyn_pf_no_ll_para =
        [gens_ph...;
         gens_qh...;
         P_non_gens...;
         Q_non_gens...;
         gens_nodes_δ_ω_ed_dash_eq_dash... ]

    if loc_load_exist == true

        dyn_pf_fun_flat_para =
            dyn_pf_wt_ll_para
        
    else

        dyn_pf_fun_flat_para =
            dyn_pf_no_ll_para
        
    end
    
    #-----------------------------------------------

    dyn_pf_fun_kwd_net_para = 
        (; 
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx )

    #-----------------------------------------------    
    
    dyn_pf_fun_kwd_net_idxs  =  
        (;
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx
         )

    #----------------------------------------

    dyn_pf_fun_kwd_n2s_idxs =
        (;
         n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx )
    
    #----------------------------------------
    
    dyn_pf_fun_kwd_nll_para_vars_Idxs = 
        (;
         P_gens_dyn_para_Idxs,
         Q_gens_dyn_para_Idxs,
         P_non_gens_dyn_para_Idxs,
         Q_non_gens_dyn_para_Idxs,
         δ_ed_eq_pf_dyn_para_Idxs
         ) 

    dyn_pf_fun_kwd_wll_para_vars_Idxs = 
        (;
         P_gens_dyn_para_Idxs,
         Q_gens_dyn_para_Idxs,
         P_non_gens_dyn_para_Idxs,
         Q_non_gens_dyn_para_Idxs,
         δ_ed_eq_pf_dyn_para_Idxs,
         P_g_loc_load_dyn_para_Idxs,
         Q_g_loc_load_dyn_para_Idxs                   
         ) 
    
    #----------------------------------------

    pf_fun_kwd_para =
        (;
         loc_load_exist,
         dyn_pf_fun_kwd_nll_para_vars_Idxs,
         dyn_pf_fun_kwd_wll_para_vars_Idxs,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_net_para,
         δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
         gens_nodes_ra_Xd_dash_Xq_dash )

    #-----------------------------------------------    
    
    full_gens_id_iq =
        (;
         gens_i_d_0,
         gens_i_q_0 )

    full_kwd_para =
        (;
         full_vars_Idxs,
         full_gens_id_iq )

    #----------------------------------------
    
    full_dyn_pf_fun_kwd_para =
        (;
         full_kwd_para,
         pf_fun_kwd_para)
    
    #----------------------------------------

    intg_dyn_pf_fun_kwd_para =
        (;
         intg_vars_Idxs,
         pf_fun_kwd_para )

            
    #-----------------------------------------------
    
    intg_vh_θh_id_iq =
        [full_vh_θh; gens_i_d_0; gens_i_q_0]

    #-----------------------------------------------

    dyn_pf_fun_var_values =
        (;
         full_vh_θh,
         intg_vh_θh_id_iq) 

    dyn_pf_fun_para_values =
        (;
         dyn_pf_wt_ll_para,
         dyn_pf_no_ll_para,
         dyn_pf_fun_flat_para )

    dyn_pf_fun_kwd_para =
        (;
         full_dyn_pf_fun_kwd_para,
         intg_dyn_pf_fun_kwd_para )

    #-----------------------------------------------
    #-----------------------------------------------

    _, gens_vh_θh_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx(
            gens_vh_θh )

    _, non_gens_vh_θh_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx(
            non_gens_vh_θh )

    ( f_gens_ωs_ωref0_vref0_porder0,
      gens_ωs_ωref0_vref0_porder0_idx_in_Idx ) =
        get_a_flattend_vars_or_paras_and_Idx(
            gens_nodes_ωs_ωref0_vref0_porder0 )

    _, gens_dyn_id_iq_pg_vh_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx(
            gens_dynamic_id_iq_pg_vh )

    _, gens_δ_ω_ed_dash_eq_dash_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx(
            gens_nodes_δ_ω_ed_dash_eq_dash )

    #----------------------------------------    

    """ length( gens_vh_θh[1] )  and
    length( gens_nodes_δ_ω_ed_dash_eq_dash[1] ) are dimensions
    of a gen_vh_θh and  gen_δ_ω_ed_dash_eq_dash respectively (2, 4)
    in flattend gen_vh_θh_δ_ω_ed_dash_eq_dash for
    a_gen_voltage_terminal_func """

    a_gen_vtf_vh_θh_δ_ω_ed_dash_eq_dash_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [length( gens_vh_θh[1] ),
             length( gens_nodes_δ_ω_ed_dash_eq_dash[1] ) ]
            ; dims_given = true )

    a_gen_vtf_vh_θh_Idx, a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx =
        a_gen_vtf_vh_θh_δ_ω_ed_dash_eq_dash_Idx

    #----------------------------------------    


    """ length( gens_nodes_ωs_ωref0_vref0_porder0[1] )  and
    length( gens_dynamic_id_iq_pg_vh[1] ) are dimensions
    of a gen_nodes_ωs_ωref0_vref0_porder0 and
    gen_dynamic_id_iq_pg_vh respectively (4, 4)
    in flattend gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh """

    a_gen_per_var_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [length( gens_dynamic_id_iq_pg_vh[1] ),
             length( gens_dynamic_id_iq_pg_vh[1] )
              ]
            ; dims_given = true )

    per_gen_ωs_ωref0_vref0_porder0_Idx, per_gen_id_iq_pg_vh_Idx =
        a_gen_per_var_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx

    #----------------------------------------    

    # per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_paras_and_Idx
    
    ode_per_para_per_node_and_idxs =
        get_per_vars_or_paras_per_node_flat_para_and_idxs(
             [ gens_nodes_ωs_ωref0_vref0_porder0,
              gens_dynamic_id_iq_pg_vh] )

    # 

    # f_per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     ode_per_para_per_node_and_idxs.flattend_per_vars_or_paras

    # 

    per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_per_para_per_node_and_idxs.per_vars_or_paras_Idx

    per_para_ωs_ωref0_vref0_porder0_Idx, per_para_id_iq_pg_vh_Idx =
        per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx
    # 
    # wrong results
    per_para_per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_per_para_per_node_and_idxs.per_vars_or_paras_per_node_Idx

    #----------------------------------------    

    # f_per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh
    
    ode_per_node_per_vars_and_idxs =
        get_per_node_per_vars_or_paras_flat_para_and_idxs(
             [gens_nodes_ωs_ωref0_vref0_porder0,
              gens_dynamic_id_iq_pg_vh ] )

    # 

    # f_per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
    #     ode_per_node_per_vars_and_idxs.flattend_per_node_vars_or_paras

    # 

    per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_per_node_per_vars_and_idxs.per_node_vars_or_paras_Idx

    # 

    per_node_per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        ode_per_node_per_vars_and_idxs.per_node_per_vars_or_paras_Idxs
    
    #----------------------------------------    
    #----------------------------------------    
   
    # f_gens_ωs_ωref0_vref0_porder0_wt_f_dyn_pf_para_and_Idx

    ode_para_and_dyn_pf_para_and_Idx =
        get_per_vars_or_paras_flat_para_and_idxs(
            [f_gens_ωs_ωref0_vref0_porder0,
             dyn_pf_fun_flat_para ] )

    # 

    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para =
        ode_para_and_dyn_pf_para_and_Idx.flattend_per_vars_or_paras

    # 

    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para_Idx =
        ode_para_and_dyn_pf_para_and_Idx.per_vars_or_paras_Idx


    f_ωs_ωref0_vref0_porder0_Idx, f_dyn_pf_para_Idx =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para_Idx

    #----------------------------------------    
    #----------------------------------------

    # f_gens_ωs_ωref0_vref0_porder0

    # f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para

    #----------------------------------------    
    #----------------------------------------    

    pf_fun_mismatch =
        get_a_model_integrated_dyn_pf_intg_ΔPQ_Δidq_mismatch

    #----------------------------------------    
    #----------------------------------------    


    ode_per_gen_models_func_kwd_paras =  [
        (;
         Ax_view,
         Bx_view,
         Cx_view,
         per_gen_ωs_ωref0_vref0_porder0_Idx,
         per_gen_id_iq_pg_vh_Idx )
        for (Ax_view, Bx_view, Cx_view) in
            zip(
              vec_Ax_views,
              vec_Bx_views,
              vec_Cx_views, ) ]
        


    ode_per_para_model_func_kwd_para =
        (;
         Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         per_para_ωs_ωref0_vref0_porder0_Idx,
         per_para_id_iq_pg_vh_Idx )
         

    
    vtf_gens_fun_kwd_para = [
        (;
         gen_ra_Xd_dash_Xq_dash,
         a_gen_vtf_vh_θh_Idx,
         a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx )

        for gen_ra_Xd_dash_Xq_dash in
            gens_nodes_ra_Xd_dash_Xq_dash ]


    vtf_para_and_idxs =
        (;
         vtf_gens_fun_kwd_para,
         gens_nodes_u_Idx_in_ranges,
         non_gens_nodes_u_Idx_in_ranges,
         )

    #----------------------------------------
    #----------------------------------------


    disaggretation_idxs =
        (; 
         gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
         a_gen_vtf_vh_θh_Idx,
         a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx,
         per_gen_ωs_ωref0_vref0_porder0_Idx,
         per_gen_id_iq_pg_vh_Idx,
         per_para_ωs_ωref0_vref0_porder0_Idx,
         per_para_id_iq_pg_vh_Idx,
         f_ωs_ωref0_vref0_porder0_Idx,
         f_dyn_pf_para_Idx )


    pre_pf_idx_and_para =
        (;
         nodes_u_Idx_in_ranges,
         
         dyn_pf_fun_kwd_wll_para_vars_Idxs,
         )

    pf_solv_para =
        (; pf_fun_mismatch,
         intg_dyn_pf_fun_kwd_para,
         pf_alg,
         abstol,
         reltol)


    post_pf_idxs =
        (; intg_vh_Idx,

         intg_θh_Idx,

         intg_id_Idx,

         intg_iq_Idx )


    #----------------------------------------    
    #----------------------------------------

    dynamics_by_per_gen_kwd_para  =
        (;
                 
         sim_state_x0,
         stateDiffCache,
         nodes_state_Idx,
         im_vars_Idx_in_state,
         loc_load_exist,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_ra_Xd_dash_Xq_dash,
         nodes_δ_ω_ed_dash_eq_dash_Idxs,
         disaggretation_idxs,
         pre_pf_idx_and_para ,
         pf_solv_para ,
         post_pf_idxs,
         vtf_para_and_idxs,
         ode_per_gen_models_func_kwd_paras,
         gens_nodes_collection,
         vec_Ax_views )


    dynamics_by_gens_kwd_para =
        (;
                  
         sim_state_x0,
         stateDiffCache,
         nodes_state_Idx,
         im_vars_Idx_in_state,
         loc_load_exist,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_ra_Xd_dash_Xq_dash,
         nodes_δ_ω_ed_dash_eq_dash_Idxs,
         disaggretation_idxs,
         pre_pf_idx_and_para ,
         pf_solv_para ,
         post_pf_idxs,
         vtf_para_and_idxs,
         ode_per_para_model_func_kwd_para,
         gens_nodes_collection,
         vec_Ax_views )

    #-------------------------------------------------
    #-------------------------------------------------
    # dynamics simulation
    #-------------------------------------------------
    #-------------------------------------------------

    # dynamics_by_ode_gens_model_and_powerflow_func!
    # dynamics_by_per_gen_ode_model_and_powerflow_func!

    sim_func  =
        dynamics_by_per_gen_ode_model_and_powerflow_func!

    #-------------------------------------------------

    sim_state_x0 = state

    sim_timespan = (0.0 , 1.0 )

     sim_fun_para =
         f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para

    #------------------------------------------------    

    counter_array  = [1]


    sim_fun_kwd_para =
        (;
         counter_array,
         dynamics_by_per_gen_kwd_para )

    #------------------------------------------------    
    
     sim_ode_func! = ODEFunction{true}(
         (dx, x, p, t) ->
        sim_func(dx, x, p, t;  kwd_para =
        sim_fun_kwd_para );
        mass_matrix =
            im_sym,
        syms =
            im_mass_matrix )
    
    sim_prob = ODEProblem(
        sim_ode_func!,
        sim_state_x0,
        sim_timespan,
        sim_fun_para )

    sim_sol =
        DifferentialEquations.solve(
            sim_prob, ode_alg, dt=dt )

    #----------------------------------------
    #########################################
    #----------------------------------------    


    sim_func  =
        dynamics_by_ode_gens_model_and_powerflow_func!

    #------------------------------------------------    

    sim_state_x0 = state

    sim_timespan = (0.0 , 1.0 )

     sim_fun_para =
         f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para

    #------------------------------------------------    


    counter_array  = [1]


    sim_fun_kwd_para =
        (;
         counter_array,

         dynamics_by_gens_kwd_para )

    #-------------------------------------------------    

    sim_ode_func! = ODEFunction{true}( ( dx, x, p, t) ->
        sim_func( dx, x, p,  t;
                 kwd_para  = sim_fun_kwd_para );
        mass_matrix =
            im_sym,
        syms =
            im_mass_matrix )
    
    sim_prob = ODEProblem(
        sim_ode_func!,
        sim_state_x0,
        sim_timespan,
        sim_fun_para )

    sim_sol =
        DifferentialEquations.solve(sim_prob, ode_alg )

    #----------------------------------------
    #########################################
    #----------------------------------------    

    recomposition_Idxs =
        (; gens_vh_θh_idx_in_Idx,
         
         non_gens_vh_θh_idx_in_Idx,
         
         gens_ωs_ωref0_vref0_porder0_idx_in_Idx, 
         
         gens_dyn_id_iq_pg_vh_idx_in_Idx,
         
         gens_δ_ω_ed_dash_eq_dash_idx_in_Idx,
         
         a_gen_vtf_vh_θh_Idx, 
         
         a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx, 

         per_gen_ωs_ωref0_vref0_porder0_Idx, 

         per_gen_id_iq_pg_vh_Idx, 

         per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx, 

         per_para_ωs_ωref0_vref0_porder0_Idx, 

         per_para_id_iq_pg_vh_Idx, 

         per_para_per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,

         per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,

         per_node_per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,

         f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para_Idx,

         f_ωs_ωref0_vref0_porder0_Idx, 

         f_dyn_pf_para_Idx 
         )


    #----------------------------------------    
    #----------------------------------------

    pre_init_data =
        (;
         Ynet,
         nodes_idx_with_adjacent_nodes_idx,
         nodes_node_idx_and_incident_edges_other_node_idx,
         net_nodes_types_idxs,
         nodes_types_idxs,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         integrated_pf_vars_and_para_idx,
         full_vars_Idxs,
         intg_vars_Idxs,
         u_ur_ui_Idx_in_state,
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges,
         gens_nodes_u_Idx_in_ranges,
         non_gens_nodes_u_Idx_in_ranges,
         n2s_idxs,
         n2s_streamlined_idx,
         dyn_PQ_δ_ω_ed_dash_eq_Idxs,
         loc_load_exist,
         net_class_names,
         im_net_states_and_var_labels,
         im_sym,
         im_mass_matrix,
         im_model_ode_mass_matrix,
         im_model_pf_mass_matrix,
         gens_nodes_im_vars_labels,
         net_bus_volts_labels,
         im_ode_sym,
         im_ode_mass_matrix,
         para_net_names_labels_syms,
         integrated_pf_PQ_param,
         gens_nodes_ra_Xd_dash_Xq_dash,
         industrial_indices_and_conversion_dict,
         im_indices_and_conversion_dict,
         indices_and_conversion_dict,
         im_vars_Idx_in_state,
         models_types_gens_δ_ω_ed_dash_eq_dash_Idxs,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_im,
         nodes_δ_ω_ed_dash_eq_dash_Idxs,
         gens_nodes_collection,
         δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
         nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         id_iq_ph_vh_idxs,
         ωs_ω_ref_v_ref_p_order_idxs,
         states_and_mat_Idxs,
         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         Ax_matrix,
         Bx_matrix,
         Cx_matrix)

    #----------------------------------------    
    #----------------------------------------    

    init_data =
        (;
         integrated_pf_and_init,
         im_state,
         industrial_state,
         ext_state,
         hybrid_state,
         state,
         im_vars_view_in_state,
         im_vars_in_state,
         sim_state_x0,
         gens_sim_state_x0) 


    #----------------------------------------    
    #----------------------------------------    

     post_init_data =
         (;
          vh,
          θh,
          full_vh_θh,
          gens_vh,
          gens_θh,
          gens_vh_θh,
          non_gens_vh,
          non_gens_θh,
          non_gens_vh_θh,
          gens_nodes_ωs_ωref0_vref0_porder0,
          gens_nodes_τm_vf,
          gens_nodes_δ_ω_ed_dash_eq_dash,
          gens_dynamic_id_iq_pg_vh,
          gens_δ,
          gens_ed_dash,
          gens_eq_dash,
          gens_ra,
          gens_Xd_dash,
          gens_Xq_dash,
          gens_id_iq,
          gens_i_d_0,
          gens_i_q_0,
          gens_vd,
          gens_vq,
          gens_ph,
          gens_qh,
          gens_ph_by_vh_θh_δ_id_iq,
          gens_qh_by_vh_θh_δ_id_iq,
          S_gens,
          gens_Pei,
          pf_fun_kwd_para,
          full_dyn_pf_fun_kwd_para,
          intg_dyn_pf_fun_kwd_para,
          intg_vh_θh_id_iq,
          dyn_pf_fun_var_values,
          dyn_pf_fun_para_values,
          dyn_pf_fun_kwd_para )


    #----------------------------------------    
    #----------------------------------------

    diagnostics_data =
        (;
         recomposition_Idxs,
         pre_init_data ,
         init_data ,
         post_init_data  )

    #----------------------------------------    
    #----------------------------------------

end




#####################################################

function driver_and_diagnostics_for_ode_and_vtf_func!( )

  
    """
    #-----------------------------------------------

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    
    only_gen  = false

    pf_alg  = NewtonRaphson()

    Idxs_type = :Idxs_im # :Idxs_hybrid, :Idxs_im, :Idxs_industrial

    #-----------------------------------------------
    
    netd =
        NetworkData( dynamics_case()... )

    """
    
    #-----------------------------------------------
    #-----------------------------------------------
    
    (; Ynet,
     nodes_idx_with_adjacent_nodes_idx ) =
         get_Ynet_and_nodes_idx_wt_adjacent_nodes_idx(
             netd )

    nodes_node_idx_and_incident_edges_other_node_idx =
        nodes_idx_with_adjacent_nodes_idx

    #-----------------------------------------------
    
    nodes_types_idxs = net_nodes_types_idxs =
        get_net_nodes_type_idxs(
            netd)
    
    slack_gens_nodes_idx =
        nodes_types_idxs.slack_gens_nodes_idx
    
    non_slack_gens_nodes_idx =
        nodes_types_idxs.non_slack_gens_nodes_idx
    
    gens_nodes_idx =
        nodes_types_idxs.gens_nodes_idx
    
    non_gens_nodes_idx =
        nodes_types_idxs.non_gens_nodes_idx
    
    gens_with_loc_load_idx =
        nodes_types_idxs.gens_nodes_with_loc_loads_idx
    
    all_nodes_idx = nodes_types_idxs.all_nodes_idx


    gens_nodes_with_loc_loads_idx =
        gens_with_loc_load_idx

    #-----------------------------------------------

    integrated_pf_vars_and_para_idx =
        get_a_model_integrated_pf_vars_and_para_idx(
            netd;
            Idxs_type =
                Idxs_type,
            no_control_device =
                false )

    # nodes_types_idxs =
    #     integrated_pf_vars_and_para_idx.nodes_types_idxs

    # n2s_idxs =
    #     integrated_pf_vars_and_para_idx.n2s_idxs

    full_types_Idxs_etc =
        integrated_pf_vars_and_para_idx.full_types_Idxs_etc


    full_gens_vh_Idxs =
        full_types_Idxs_etc.full_gens_vh_Idxs
    
    full_gens_θh_Idxs  =
        full_types_Idxs_etc.full_gens_θh_Idxs
    
    full_non_gens_nodes_vh_Idxs  =
        full_types_Idxs_etc.full_non_gens_nodes_vh_Idxs
    
    full_non_gens_nodes_θh_Idxs  =
        full_types_Idxs_etc.full_non_gens_nodes_θh_Idxs

    # ----
   
    full_vars_Idxs = 
        (;full_gens_vh_Idxs,
         full_gens_θh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_non_gens_nodes_θh_Idxs)
    
    #----------------------------------------    
    
    intg_types_Idxs_etc =
        integrated_pf_vars_and_para_idx.intg_types_Idxs_etc


    intg_gens_vh_Idxs =
        intg_types_Idxs_etc.intg_gens_vh_Idxs
    
    intg_gens_θh_Idxs =
        intg_types_Idxs_etc.intg_gens_θh_Idxs
    
    intg_non_gens_nodes_vh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_vh_Idxs
    
    intg_non_gens_nodes_θh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_θh_Idxs

    intg_gen_id_Idxs =
        intg_types_Idxs_etc.intg_gen_id_Idxs
    
    intg_gen_iq_Idxs =
        intg_types_Idxs_etc.intg_gen_iq_Idxs
    
    # ----
   
    intg_vars_Idxs = 
        (;
         intg_gens_vh_Idxs,
         intg_gens_θh_Idxs,
         intg_non_gens_nodes_vh_Idxs,
         intg_non_gens_nodes_θh_Idxs,
         intg_gen_id_Idxs,
         intg_gen_iq_Idxs )

    intg_kwd_para =
        (; intg_vars_Idxs, )
    
    #----------------------------------------    
    
    u_ur_ui_Idx_in_state =
        integrated_pf_vars_and_para_idx.u_ur_ui_Idx_in_state

   (; slack_ur_ui_Idx_in_state, 
    non_slack_ur_ui_Idx_in_state, 
    ur_ui_Idx_in_state, 
    nodes_u_Idx,
    nodes_u_Idx_in_ranges ) =
        u_ur_ui_Idx_in_state

    
    #-----------------------------------------------

    (;
     gens_nodes_u_Idx_in_ranges,
     non_gens_nodes_u_Idx_in_ranges ) =
         gens_and_non_gens_u_Idx_in_ranges(
             all_nodes_idx,
             gens_nodes_idx,
             non_gens_nodes_idx,
                 nodes_u_Idx_in_ranges)
        
    #-----------------------------------------------

    n2s_idxs = n2s_streamlined_idx =
        get_dict_n2s_streamlined_idx(netd)
    
    n2s_slack_gens_idx =
        n2s_idxs.n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        n2s_idxs.n2s_non_slack_gens_idx
    
    n2s_gens_idx = n2s_idxs.n2s_gens_idx
    
    n2s_non_gens_idx =
        n2s_idxs.n2s_non_gens_idx
    
    n2s_gens_with_loc_load_idxs =
        n2s_idxs.n2s_gens_with_loc_load_idxs
    
    n2s_all_nodes_idx =
        n2s_idxs.n2s_all_nodes_idx

    #-----------------------------------------------

    dyn_PQ_δ_ω_ed_dash_eq_Idxs  =
         get_a_model_integrated_pf_PQ_δ_ω_ed_dash_eq_dash_Idxs(
             netd)

    dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs =
        dyn_PQ_δ_ω_ed_dash_eq_Idxs.dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs
    
    (;
     P_gens_dyn_pf_no_ll_para_Idxs,
     Q_gens_dyn_pf_wt_ll_para_Idxs,
     P_non_gens_dyn_pf_wt_ll_para_Idxs,
     Q_non_gens_dyn_pf_wt_ll_para_Idxs,
     δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs,
     P_g_loc_load_dyn_pf_wt_ll_para_Idxs,
     Q_g_loc_load_dyn_pf_wt_ll_para_Idxs) =
         dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs

    P_gens_dyn_para_Idxs =
        P_gens_dyn_pf_no_ll_para_Idxs
    
    Q_gens_dyn_para_Idxs =
        Q_gens_dyn_pf_wt_ll_para_Idxs
    
    P_non_gens_dyn_para_Idxs =
        P_non_gens_dyn_pf_wt_ll_para_Idxs
    
    Q_non_gens_dyn_para_Idxs =
        Q_non_gens_dyn_pf_wt_ll_para_Idxs
    
    δ_ed_eq_pf_dyn_para_Idxs =
        δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs
    
    P_g_loc_load_dyn_para_Idxs =
        P_g_loc_load_dyn_pf_wt_ll_para_Idxs
    
    Q_g_loc_load_dyn_para_Idxs =
        P_g_loc_load_dyn_pf_wt_ll_para_Idxs
    
    #-----------------------------------------------
        
    loc_load_exist =
        loc_load_exist_bool( netd )
    
    #-----------------------------------------------        
    #  label syms and mass_matrix
    #-----------------------------------------------
    
    net_class_names =
           make_case_buses_names(
               netd.nodes )

    #-----------------------------------------------

    im_net_states_and_var_labels =
         generate_im_model_labels(
             ;nodes =
                 netd.nodes )

    #-----------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes =
                netd.nodes )

    (; im_model_ode_mass_matrix,
     im_model_pf_mass_matrix ) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes =
                 netd.nodes )

    (; gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes =
                 netd.nodes  )

    im_ode_sym =
        gens_nodes_im_vars_labels

    im_ode_mass_matrix =
        im_model_ode_mass_matrix
        
    #-----------------------------------------------

    para_net_names_labels_syms =
        (;
         net_class_names,
         im_net_states_and_var_labels,
         im_ode_sym )

    #-----------------------------------------------    
    #-----------------------------------------------    
    
    integrated_pf_PQ_param =
        get_a_model_integrated_pf_PQ_param(
            netd )
    (;
     P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         integrated_pf_PQ_param

    #-----------------------------------------------

    gens_nodes_ra_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            netd.nodes_param,
            netd.nodes;
            param_list =
                [ :ra, :X_d_dash, :X_q_dash ],
            gens_view_only = true )
    
    #-----------------------------------------------

    industrial_indices_and_conversion_dict =
        get_industrial_model_indices_and_conversion_dict(
        netd; no_control_device = false )
    
   im_indices_and_conversion_dict =
         get_im_indices_and_conversion_dict(
             netd  )

    indices_and_conversion_dict =
        im_indices_and_conversion_dict

    im_vars_Idx_in_state =
        indices_and_conversion_dict.im_vars_Idx_in_state
    
    #-----------------------------------------------
    
    models_types_gens_δ_ω_ed_dash_eq_dash_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )

    nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_hybrid

    nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_industrial

    nodes_δ_ω_ed_dash_eq_dash_Idxs_im =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_im
    
    nodes_δ_ω_ed_dash_eq_dash_Idxs =
        nodes_δ_ω_ed_dash_eq_dash_Idxs_im
    
    #-----------------------------------------------
    
    gens_nodes_collection =
        get_gens_nodes(
            netd.nodes )

    #-----------------------------------------------
    
    δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        get_gens_nodes_some_state_vars_Idxs_in_flattend(
            gens_nodes_collection
            ; some_state_vars =
                [ :δ, :ω, :ed_dash, :eq_dash ] )
    
    #-----------------------------------------------
    
    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, states_and_mat_Idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(
            gens_nodes_collection )

    (;vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views )  =
        Ax_Bx_Cx_views
    
    (;Ax_matrix,
     Bx_matrix,
     Cx_matrix ) =
         Ax_Bx_Cx_matrix

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ωs_ω_ref_v_ref_p_order_idxs )  =
         states_and_mat_Idxs
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views ),
        gens_nodes_collection )

    im_plants_system_matrices =
        (;
         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         
         Ax_matrix,
         Bx_matrix,
         Cx_matrix )


    
    #-----------------------------------------------
    #-----------------------------------------------
    
    integrated_pf_and_init =
         get_a_model_integrated_sta_powerflow_and_init(
             netd;
             pf_alg  =
                 pf_alg)
        
    #-----------------------------------------------
    
    nodes_name =
        integrated_pf_and_init.nodes_name

    branches_name =
        integrated_pf_and_init.branches_name

    gens_nodes_idx =
        integrated_pf_and_init.gens_nodes_idx

    #-----------------------------------------------

    named_tup_pf_result =
        integrated_pf_and_init.named_tup_pf_result

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init

    pf_init_dict = named_tup_pf_result.pf_init_dict
    
    #-----------------------------------------------
    
    im_state =
        im_model_init_operationpoint(
            netd, bus_dict_init  )
    
    industrial_state =
        industrial_model_init_operationpoint(
            netd, bus_dict_init; pure = :pure,
            no_control_device = only_gen )

    ext_state =
        external_init_operationpoint(
            netd, bus_dict_init,
            branch_dict_init )

    hybrid_state =
        init_operationpoint(netd, pf_init_dict)
    
    #-----------------------------------------------
    
    state = im_state
    
    #-----------------------------------------------
    
    im_vars_view_in_state =
        get_im_vars_view_in_state(
            state,
            im_vars_Idx_in_state )

    
    im_vars_in_state =
        state[im_vars_Idx_in_state ]
    
    #-----------------------------------------------

    sim_state_x0 =
        state[ im_vars_Idx_in_state ]

    
    gens_sim_state_x0 =
        [ sim_state_x0[ idx ]
         for idx in
             nodes_state_Idx ]
    
    #-----------------------------------------------
    
    vh = named_tup_pf_result.Vm
    
    θh = named_tup_pf_result.Vθ

    #-----------------------------------------------

    full_vh_θh = [vh; θh]

    #-----------------------------------------------

    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]


    gens_vh_θh =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(gens_vh , gens_θh)]
    
    #-----------------------------------------------

    non_gens_vh = vh[non_gens_nodes_idx]

    non_gens_θh = θh[non_gens_nodes_idx]

    non_gens_vh_θh =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(non_gens_vh , non_gens_θh)]

    #-----------------------------------------------
    
    gens_nodes_ωs_ωref0_vref0_porder0 =
        get_gens_nodes_ωs_ωref0_vref0_porder0(
        gens_nodes_collection,
        bus_dict_init )


    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )
    
    #-----------------------------------------------

    gens_nodes_δ_ω_ed_dash_eq_dash =
        get_gen_im_nodes_ω_ed_dash_eq_dash(
            state,
            nodes_δ_ω_ed_dash_eq_dash_Idxs )

    #-----------------------------------------------
    
    gens_dynamic_id_iq_pg_vh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
        gens_vh_θh,
        gens_nodes_δ_ω_ed_dash_eq_dash,
            gens_nodes_ra_Xd_dash_Xq_dash )
    
    #-----------------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_ed_dash =
        third.( gens_nodes_δ_ω_ed_dash_eq_dash )

    gens_eq_dash =
        fourth.( gens_nodes_δ_ω_ed_dash_eq_dash )    
    
    #-----------------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )        

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
        for (a_vh,
             a_θh,
             a_δ,
             a_ed_dash,
             a_eq_dash,
             a_ra,
             a_X_d_dash,
             a_X_q_dash ) in
            zip( gens_vh,
                 gens_θh,
                 gens_δ,
                 gens_ed_dash,
                 gens_eq_dash,
                 gens_ra,
                 gens_Xd_dash,
                 gens_Xq_dash ) ]

    gens_i_d_0 =
        first.( gens_id_iq )
    
    gens_i_q_0 =
        last.( gens_id_iq )
    
    #-----------------------------------------------

    gens_vd = [
        a_ed_dash +
            a_Xq_dash * a_iq -
            a_ra * a_id
        for ( a_ed_dash,
              a_Xq_dash,
              a_ra,
              a_id,
              a_iq ) in
            zip(gens_ed_dash,
                gens_Xq_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]

    gens_vq = [
        a_eq_dash - a_Xd_dash * a_id - a_ra * a_id
        for ( a_eq_dash, a_Xd_dash, a_ra, a_id, a_iq ) in
            zip(gens_eq_dash,
                gens_Xd_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]    
    
    #-----------------------------------------------

    gens_ph = [ a_vd * a_id + a_vq * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    gens_qh = [ a_vq * a_id - a_vd * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]

    #-----------------------------------------------

    # a_id * a_vh * sin(a_δ - a_θh) + a_iq * a_vh * cos(a_δ - a_θh)
    
    gens_ph_by_vh_θh_δ_id_iq = [
        get_a_gen_active_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    # a_id * a_vh * cos(a_δ - a_θh) - a_iq * a_vh * sin(a_δ - a_θh)
    
    gens_qh_by_vh_θh_δ_id_iq = [
        get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]
    
    #-----------------------------------------------
    
    S_gens = [
        a_vh * exp(im * (a_θh - a_δ + π/2)) * ( a_id - im * a_iq )
        for ( a_vh, a_θh, a_δ, a_id, a_iq) in
            zip(gens_vh,
                gens_θh,
                gens_δ,
                gens_i_d_0,
                gens_i_q_0) ]
    
    #-----------------------------------------------
    
    gens_Pei = [
        ed_dash * i_d_0 + eq_dash * i_q_0 +
            (Xq_dash - Xd_dash ) * i_d_0 *  i_q_0
        for (ed_dash,eq_dash,Xd_dash,Xq_dash,i_d_0,i_q_0 ) in
            zip(gens_ed_dash, gens_eq_dash,
                gens_Xd_dash,
                gens_Xq_dash,
                gens_i_d_0, gens_i_q_0 ) ]

    #----------------------------------------    
    #----------------------------------------
    # dyn_pf_para
    #----------------------------------------
    #----------------------------------------    

    dyn_pf_wt_ll_para =
        [gens_ph...;
         gens_qh...;
         P_non_gens...;
         Q_non_gens...;
         gens_nodes_δ_ω_ed_dash_eq_dash...;
         P_g_loc_load...;
         Q_g_loc_load...]


    dyn_pf_no_ll_para =
        [gens_ph...;
         gens_qh...;
         P_non_gens...;
         Q_non_gens...;
         gens_nodes_δ_ω_ed_dash_eq_dash... ]

    if loc_load_exist == true

        dyn_pf_fun_flat_para =
            dyn_pf_wt_ll_para
        
    else

        dyn_pf_fun_flat_para =
            dyn_pf_no_ll_para
        
    end
    
    #-----------------------------------------------

    dyn_pf_fun_kwd_net_para = 
        (; 
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx )

    #-----------------------------------------------    
    
    dyn_pf_fun_kwd_net_idxs  =  
        (;
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx
         )

    #----------------------------------------

    dyn_pf_fun_kwd_n2s_idxs =
        (;
         n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx )
    
    #----------------------------------------
    
    dyn_pf_fun_kwd_nll_para_vars_Idxs = 
        (;
         P_gens_dyn_para_Idxs,
         Q_gens_dyn_para_Idxs,
         P_non_gens_dyn_para_Idxs,
         Q_non_gens_dyn_para_Idxs,
         δ_ed_eq_pf_dyn_para_Idxs
         ) 

    dyn_pf_fun_kwd_wll_para_vars_Idxs = 
        (;
         P_gens_dyn_para_Idxs,
         Q_gens_dyn_para_Idxs,
         P_non_gens_dyn_para_Idxs,
         Q_non_gens_dyn_para_Idxs,
         δ_ed_eq_pf_dyn_para_Idxs,
         P_g_loc_load_dyn_para_Idxs,
         Q_g_loc_load_dyn_para_Idxs                   
         ) 
    
    #----------------------------------------

    pf_fun_kwd_para =
        (;
         loc_load_exist,
         dyn_pf_fun_kwd_nll_para_vars_Idxs,
         dyn_pf_fun_kwd_wll_para_vars_Idxs,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_net_para,
         δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
         gens_nodes_ra_Xd_dash_Xq_dash )

    #-----------------------------------------------    
    
    full_gens_id_iq =
        (;
         gens_i_d_0,
         gens_i_q_0 )

    full_kwd_para =
        (;
         full_vars_Idxs,
         full_gens_id_iq )

    #----------------------------------------
    
    full_dyn_pf_fun_kwd_para =
        (;
         full_kwd_para,
         pf_fun_kwd_para)
    
    #----------------------------------------

    intg_dyn_pf_fun_kwd_para =
        (;
         intg_vars_Idxs,
         pf_fun_kwd_para )

            
    #-----------------------------------------------
    
    intg_vh_θh_id_iq =
        [full_vh_θh; gens_i_d_0; gens_i_q_0]

    #-----------------------------------------------

    dyn_pf_fun_var_values =
        (;
         full_vh_θh,
         intg_vh_θh_id_iq) 

    dyn_pf_fun_para_values =
        (;
         dyn_pf_wt_ll_para,
         dyn_pf_no_ll_para,
         dyn_pf_fun_flat_para )

    dyn_pf_fun_kwd_para =
        (;
         full_dyn_pf_fun_kwd_para,
         intg_dyn_pf_fun_kwd_para )

    #-----------------------------------------------
    #-----------------------------------------------

  
  # (;
  #    gens_nodes_vh_θh_idx_in_Idx,
  #    gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
  #    gens_nodes_id_iq_pg_vh_idx_in_Idx ) =
  #        get_ode_flat_para_Idxs_in_Idxs(
  #            gens_vh_θh,
  #            gens_nodes_ωs_ωref0_vref0_porder0,
  #            gens_dynamic_id_iq_pg_vh )

    f_gens_nodes_vh_θh, gens_nodes_vh_θh_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx( gens_vh_θh )

    f_non_gens_vh_θh, non_gens_vh_θh_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx( non_gens_vh_θh )

    f_gens_nodes_ωs_ωref0_vref0_porder0, gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx( gens_nodes_ωs_ωref0_vref0_porder0 )

    f_gens_dynamic_id_iq_pg_vh, gens_dynamic_id_iq_pg_vh_idx_in_Idx =
        get_a_flattend_vars_or_paras_and_Idx( gens_dynamic_id_iq_pg_vh )

    f_gens_nodes_δ_ω_ed_dash_eq_dash, gens_nodes_δ_ω_ed_dash_eq_dash_idx_in_Idx = get_a_flattend_vars_or_paras_and_Idx(gens_nodes_δ_ω_ed_dash_eq_dash)

    #----------------------------------------    

    """ length(gens_vh_θh[1])  and
    length(gens_nodes_δ_ω_ed_dash_eq_dash[1]) are dimensions
    of a gen_vh_θh and  gen_δ_ω_ed_dash_eq_dash respectively (2, 4)
    in flattend gen_vh_θh_δ_ω_ed_dash_eq_dash for
    a_gen_voltage_terminal_func """

    a_gen_vtf_vh_θh_δ_ω_ed_dash_eq_dash_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [length(gens_vh_θh[1]),
             length(gens_nodes_δ_ω_ed_dash_eq_dash[1]) ]
            ; dims_given = true )

    a_gen_vtf_vh_θh_Idx,  a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx =
        a_gen_vtf_vh_θh_δ_ω_ed_dash_eq_dash_Idx

    #----------------------------------------    

    gens_per_node_vtf_para_and_Idx =
        get_per_node_flat_para_and_idxs(
            [gens_vh_θh,
             gens_nodes_δ_ω_ed_dash_eq_dash ] )


    gens_per_node_vtf_flat_paras =
        gens_per_node_vtf_para_and_Idx.flattend_per_node_vars_or_paras

    gens_per_node_vtf_paras_flat_Idx =
        gens_per_node_vtf_para_and_Idx.per_node_vars_or_paras_Idx


    #----------------------------------------    
    
    # flat_per_vars_vh_θh_ωs_ωref0_vref0_porder0_id_iq_pg_vh
    
    ode_per_vars_per_node_and_idxs =
        get_per_vars_or_paras_per_node_flat_para_and_idxs(
             [gens_vh_θh,
              gens_nodes_ωs_ωref0_vref0_porder0,
              gens_dynamic_id_iq_pg_vh] )

    ode_per_vars_flattend =
        ode_per_vars_per_node_and_idxs.flattend_per_vars_or_paras

    ode_per_vars_or_paras_Idx =
        ode_per_vars_per_node_and_idxs.per_vars_or_paras_Idx
    
    ode_per_vars_or_paras_per_node_Idx =
        ode_per_vars_per_node_and_idxs.per_vars_or_paras_per_node_Idx

    #----------------------------------------    

    # flat_per_node_vh_θh_ωs_ωref0_vref0_porder0_id_iq_pg_vh
    
    ode_per_node_per_vars_and_idxs =
        get_per_node_per_vars_or_paras_flat_para_and_idxs(
             [gens_vh_θh,
             gens_nodes_ωs_ωref0_vref0_porder0,
             gens_dynamic_id_iq_pg_vh ] )


    ode_per_node_flattend =
        ode_per_node_per_vars_and_idxs.flattend_per_node_vars_or_paras

    ode_per_node_paras_Idx =
        ode_per_node_per_vars_and_idxs.per_node_vars_or_paras_Idx

    ode_per_node_per_vars_or_paras_Idx =
        ode_per_node_per_vars_and_idxs.per_node_per_vars_or_paras_Idxs
    

    #----------------------------------------    

    ode_per_vars_wt_dyn_pf_para_and_Idx =
        get_per_vars_or_paras_flat_para_and_idxs(
            [ode_per_vars_flattend,
             dyn_pf_fun_flat_para ] )

    ode_per_vars_wt_dyn_pf_para =
        ode_per_vars_wt_dyn_pf_para_and_Idx.flattend_per_vars_or_paras

    ode_per_vars_wt_dyn_pf_Idx =
        ode_per_vars_wt_dyn_pf_para_and_Idx.per_vars_or_paras_Idx

    #----------------------------------------    

    ode_per_node_wt_dyn_pf_para_and_Idx =
        get_per_node_flat_para_and_idxs(
            [ode_per_node_flattend,
             dyn_pf_fun_flat_para ] )

    ode_per_node_wt_dyn_pf_flat_para =
        ode_per_node_wt_dyn_pf_para_and_Idx.flattend_per_node_vars_or_paras

    ode_per_node_wt_dyn_pf_flat_Idx =
        ode_per_node_wt_dyn_pf_para_and_Idx.per_node_vars_or_paras_Idx

    #----------------------------------------    
    # f_non_gens_vh_θh
    #----------------------------------------    

    ode_per_node_wt_dyn_pf_wt_non_gens_vh_θh_para_and_Idx =
        get_per_node_flat_para_and_idxs(
            [ode_per_node_wt_dyn_pf_flat_para,
             f_non_gens_vh_θh ] )

    ode_per_node_wt_dyn_pf_wt_non_gens_vh_θh_para =
        ode_per_node_wt_dyn_pf_wt_non_gens_vh_θh_para_and_Idx.flattend_per_node_vars_or_paras

    ode_per_node_wt_dyn_pf_wt_non_gens_vh_θh_Idx =
        ode_per_node_wt_dyn_pf_wt_non_gens_vh_θh_para_and_Idx.per_node_vars_or_paras_Idx

    #----------------------------------------    

   
   # gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx

    ode_para_and_dyn_pf_para_and_Idx =
        get_per_vars_or_paras_flat_para_and_idxs(
            [f_gens_nodes_ωs_ωref0_vref0_porder0,
             dyn_pf_fun_flat_para ] )

    ode_para_and_dyn_pf_para =
        ode_para_and_dyn_pf_para_and_Idx.flattend_per_vars_or_paras

    ode_para_and_dyn_pf_Idx =
        ode_para_and_dyn_pf_para_and_Idx.per_vars_or_paras_Idx

     #----------------------------------------
     #########################################
     #----------------------------------------    

    ode_per_gens_wt_dyn_pf_Idx, non_gens_vh_θh_para_Idx =
        ode_per_node_wt_dyn_pf_wt_non_gens_vh_θh_Idx

    ode_per_gens_wt_dyn_pf_para =
        ode_per_node_wt_dyn_pf_wt_non_gens_vh_θh_para[
            ode_per_gens_wt_dyn_pf_Idx]

    non_gens_vh_θh_para =
        ode_per_node_wt_dyn_pf_wt_non_gens_vh_θh_para[
            non_gens_vh_θh_para_Idx]

    #----------------------------------------    
    #----------------------------------------

    ode_per_gens_para_Idx, dyn_pf_flat_para_Idx =
        ode_per_node_wt_dyn_pf_flat_Idx

    ode_per_gens_para =
        ode_per_node_wt_dyn_pf_para[
            ode_per_gens_para_Idx]

    dyn_pf_flat_para =
        ode_per_node_wt_dyn_pf_para[
            dyn_pf_flat_para_Idx ]

    #----------------------------------------    
    #----------------------------------------

    list_ode_per_gens_para =
        [ ode_per_gens_para[idx]
          for idx in
              ode_per_node_paras_Idx]

    #----------------------------------------    
    #----------------------------------------

    list_vtf_non_gens_para =
        [non_gens_vh_θh_para[idx]
         for idx in
             non_gens_vh_θh_idx_in_Idx ]

    #----------------------------------------    
    #----------------------------------------

    dx_gens_states_view =
        view( dx, im_vars_Idx_in_state ) 

    x_gens_states_view =
        view( x, im_vars_Idx_in_state ) 

    #----------------------------------------    
    #----------------------------------------

    list_dx_gens_states_views =
        [ view( dx, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]

    list_x_gens_states_views =
        [ view( x, a_gen_state_Idx )
          for a_gen_state_Idx in
              nodes_state_Idx ]
    
    #----------------------------------------    
    #----------------------------------------

    vtf_dx_non_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            non_gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_non_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           non_gens_nodes_u_Idx_in_ranges ]

    #----------------------------------------    
    #----------------------------------------

    vtf_dx_gens_u_views = [
        view( dx, node_u_Idx )
        for node_u_Idx in
            gens_nodes_u_Idx_in_ranges ]
    
    vtf_x_gens_u_views = [
        view( x, node_u_Idx )
        for  node_u_Idx in
           gens_nodes_u_Idx_in_ranges ]

    #----------------------------------------    
    #----------------------------------------

    vtf_gens_vh_θh_δ_ω_ed_dash_eq_dash =
        [ gens_per_node_vtf_flat_paras[
            per_node_vtf_paras_flat_Idx]
         for per_node_vtf_paras_flat_Idx in
             gens_per_node_vtf_paras_flat_Idx ]

    #----------------------------------------    
    #----------------------------------------
    
    vtf_gens_fun_kwd_para = [
        (;
         gen_ra_Xd_dash_Xq_dash,
         a_gen_vtf_vh_θh_Idx,
         a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx )

        for gen_ra_Xd_dash_Xq_dash in
            gens_nodes_ra_Xd_dash_Xq_dash ]

    #----------------------------------------    
    #----------------------------------------

    for (vtf_dx, vtf_x, vtf_para, vtf_kwd_para) in zip(
        vtf_dx_gens_u_views,
        vtf_x_gens_u_views,
        vtf_gens_vh_θh_δ_ω_ed_dash_eq_dash,
        vtf_gens_fun_kwd_para )

        a_gen_voltage_terminal_func!(
            vtf_dx, vtf_x, vtf_para
            ; vtf_kwd_para = vtf_kwd_para )
    end


end


