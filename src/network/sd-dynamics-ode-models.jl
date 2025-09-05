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
        [ gens_dynamic_id_iq_pg_vh...;] +
        Cx_matrix * [gens_nodes_ωs_τm_vref_porder_view...;]

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
            Bx_matrix * [id_iq_pg_vh...; ] +
            Cx_matrix * [ωs_τm_vref_porder...;]
        
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
            Bx_matrix * [id_iq_pg_vh...;] +
            Cx_matrix * [ωs_τm_vref_porder...;] +
            Ax_τm_vf_matrix * [τm_vf...;]

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
        sim_fun_system_kwd_flat_agg_para)

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

    (;vh_θh_Idx,     
     ωs_ωref0_vref0_porder0_Idx) =
         sim_fun_system_ode_flat_agg_para_Idxs
    
    f_gens_vh_θh =
        sim_fun_system_ode_flat_agg_para[
            vh_θh_Idx]
    
    f_ωs_ωref0_vref0_porder0 =
        sim_fun_system_ode_flat_agg_para[
            ωs_ωref0_vref0_porder0_Idx] 
    
    gens_vh_θh =
        [ f_gens_vh_θh[idx]
          for idx in
              gens_nodes_vh_θh_indx_in_Idx]

    gens_ωs_ωref0_vref0_porder0 =
        [ f_ωs_ωref0_vref0_porder0[idx]
          for idx in
              gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx ]

    gens_vh_θh_ωs_ωref0_vref0_porder0 =
        [[a_vh_θh...;
          a_ωs_ωref0_vref0_porder0...]
          for (a_vh_θh, a_ωs_ωref0_vref0_porder0) in
              zip(gens_vh_θh,
                  gens_ωs_ωref0_vref0_porder0)]

    
    (;nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs ) =
         sys_states_Idxs_and_mat_Idxs
       
    dx_gens  = [view(dx, idx)
                for idx in  nodes_state_Idx]
    
    x_gens   = [view(x, idx)
                for idx in nodes_state_Idx]
    
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
                 gen_nodes_ra_Xd_dash_Xq_dash_view) ]

    #---------------------------------------------------
    
    for (a_dx, a_x, a_ode_para, a_sim_fun_kwd_para) in
        zip( dx_gens, x_gens,
             gens_vh_θh_ωs_ωref0_vref0_porder0,
             vec_sim_fun_kwd_para )

        t_ode_one_im_model_func!(
            a_dx, a_x, a_ode_para, t;
            sim_fun_kwd_para =
                a_sim_fun_kwd_para)
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

    (working_vh_θh_view,
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

    im_dyn_powerflow(
        netd,
        stateDiffCache,
        global_pf_param,
        im_dyn_pf_up_para,
        im_idq_pf_cal;
        dyn_global_pf_options... )

    # industrial_dyn_powerflow(
    # netd,
    # stateDiffCache,
    # global_pf_param,
    # im_dyn_pf_up_para,
    # im_idq_pf_cal;
    # dyn_global_pf_options... )

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
    
    (Pg_Idx, Qg_Idxs,
     Png_Idxs, Qng_Idxs,
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
    
    (Pg_Idx, Qg_Idxs,
     Png_Idxs, Qng_Idxs,
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
    
    (Pg_Idx, Qg_Idxs,
     Png_Idxs, Qng_Idxs,
     Pgll_Idxs, Qgll_Idxs) =
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
    
    # gens_nodes_δ_ed_dash_eq_dash =
    #     [ flat_δ_ed_dash_eq_dash[idx]
    #       for idx in
    #           δ_ed_dash_eq_dash_Idxs_in_flattend ]

    
    # gens_nodes_δ_ed_dash_eq_dash =
    #     [ flat_δ_ed_dash_eq_dash[idx]
    #       for idx in
    #           nodes_δ_ed_dash_eq_dash_Idxs ]
    
    
    
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

    flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0 =
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

    (flat_vh_θh_idx,
     flat_ωs_ωref0_vref0_porder0_idx,
     flat_init_dyn_pf_para_idx) =
         flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx
    
    #----------------------------------------    

    flat_vh_θh =
        flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[
            flat_vh_θh_idx ]
    
    flat_ωs_ωref0_vref0_porder0 =
        flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[
            flat_ωs_ωref0_vref0_porder0_idx ]
    
    init_dyn_pf_flat_para =
        flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_flat_para[
            flat_init_dyn_pf_para_idx ]

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
    
    # gens_nodes_δ_ed_dash_eq_dash =
    #     [ flat_δ_ed_dash_eq_dash[idx]
    #       for idx in
    #           δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
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

    flat_vh_flat_θh_flat_ωs_ωref0_vref0_porder0 =
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

    (flat_vh_θh_idx,
     flat_ωs_ωref0_vref0_porder0_idx,
     flat_init_dyn_pf_para_idx) =
         flat_vh_θh_ωs_ωref0_vref0_porder0_init_dyn_pf_para_Idx
    
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
             zip( gens_vh_post_pf,
                  gens_θh_post_pf ) ]
    
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
            [gens_vh_θh_post_pf,
             gens_nodes_δ_ed_dash_eq_dash] )
    
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

