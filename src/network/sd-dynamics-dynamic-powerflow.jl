# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.4:  pg 135 - 136, 6.121 - 6.123

########################################################
#-------------------------------------------------------
# Begining of functions defiitions
#-------------------------------------------------------
########################################################


function get_static_powerflow_sol_by_json(
    pf_PQ_param;
    kwd_para =
        kwd_sta_sta_ΔPQ_sol_by_pf_type,
    maxiters = 1e5,
    abstol = 1000*eps(),
    reltol = 1000*eps()
    )
    
    #-----------------------------------------------

    (;
     pf_alg,
     pf_kw_para,
     red_vh_Idxs,
     red_θh_Idxs,
     sta_red_vh_θh_0) =
         kwd_para

    #-----------------------------------------------
    
    #-----------------------------------------------    
    # Powerflow func and prob
    #-----------------------------------------------
    
    pf_fun_mismatch =
        get_a_model_integrated_pf_sta_ΔPQ_mismatch_generic
    
    return  NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( ( g, x, p ) ->
                pf_fun_mismatch(
                    g, x, p;
                    pf_kw_para =
                        pf_kw_para ) ),
            sta_red_vh_θh_0,
            pf_PQ_param ),
        pf_alg;
        maxiters =
            maxiters )

    
end



function get_pf_spcm_ΔPQ_mismatch_sol_by_generic(
    pf_PQ_param;
    kwd_para =
        kwd_sta_sta_ΔPQ_sol_by_pf_type,
    # maxiters = 1e5,
    abstol = 1000*eps(),
    reltol = 1000*eps()
    )
    
    #-----------------------------------------------

    (;
     pf_alg,
     pf_kw_para,
     red_vh_Idxs,
     red_θh_Idxs,
     sta_red_vh_θh_0) =
         kwd_para
    
    #-----------------------------------------------    
    # Powerflow func and prob
    #-----------------------------------------------
    
    pf_fun_mismatch =
        get_a_model_integrated_pf_spcm_ΔPQ_mismatch_generic

    # pf_sol
    
    return  NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( ( g, x, p ) ->
                pf_fun_mismatch(
                    g, x, p;
                    pf_kw_para =
                        pf_kw_para ) ),
            sta_red_vh_θh_0,
            pf_PQ_param ),
        pf_alg ;
        abstol =
            abstol,
        reltol =
            reltol )
    
end



function get_pf_sta_ΔPQ_mismatch_sol_by_generic(
    pf_PQ_param;
    kwd_para =
        kwd_sta_sta_ΔPQ_sol_by_pf_type,
    # maxiters = 1e5,
    abstol = 1000*eps(),
    reltol = 1000*eps()
    )
    
    #-----------------------------------------------

    (;
     pf_alg,
     pf_kw_para,
     red_vh_Idxs,
     red_θh_Idxs,
     sta_red_vh_θh_0) =
         kwd_para
    
    #-----------------------------------------------    
    # Powerflow func and prob
    #-----------------------------------------------
    
    # pf_fun_mismatch =
    #     get_a_model_integrated_pf_sta_ΔPQ_mismatch
    
    pf_fun_mismatch =
        get_a_model_integrated_pf_sta_ΔPQ_mismatch_generic

    # pf_sol
    
    return  NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( ( g, x, p ) ->
                pf_fun_mismatch(
                    g, x, p;
                    pf_kw_para =
                        pf_kw_para ) ),
            sta_red_vh_θh_0,
            pf_PQ_param ),
        pf_alg ;
        abstol =
            abstol,
        reltol =
            reltol )
    
end


function get_oop_pf_sta_ΔPQ_mismatch_sol_by_generic(
    pf_PQ_param;
    kwd_para =
        kwd_sta_sta_ΔPQ_sol_by_mpc )
    
    #-----------------------------------------------

    (;
     pf_alg,
     pf_kw_para,
     red_vh_Idxs,
     red_θh_Idxs,
     sta_red_vh_θh_0) =
         kwd_para
        
    #-----------------------------------------------    
    # Powerflow func and prob
    #-----------------------------------------------
    
    pf_fun_mismatch =
        get_a_model_integrated_pf_sta_ΔPQ_mismatch_generic

    # pf_sol
    
    return  NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( ( x, p ) ->
                pf_fun_mismatch(
                    x, p;
                    pf_kw_para =
                        pf_kw_para)),
            sta_red_vh_θh_0,
            pf_PQ_param ),
        pf_alg )
    
end

#-----------------------------------------------------
# new 
#-----------------------------------------------------

function dynamic_power_balance_powerflow_with_vh_θh!(
    x0_vh_θh,
    mismatch,
    (working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view),
    pf_net_param ;
    maxiter= 1e5,
    ftol=1000*eps(),
    xtol=1000*eps(),
    with_δ_ed_eq = true )

    # ---------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
     nodes_u_Idx,
     gens_idx,
     ur_IDX,
     ui_IDX,
     vh_IDX,
     θh_IDX,
     red_vh_θh_idx,
     ur_idx,
     ui_idx,
     ur_ui_idx)=
         pf_Idx
    
    #------------------------------------------    
    
    x0_vh_view    = @view x0_vh_θh[vh_IDX]
    x0_θh_view    = @view x0_vh_θh[θh_IDX]

    # Get a view of a reduced vh_θh
    
    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]
    
    #------------------------------------------    
    
    gens_uh = x0_vh_view[ gens_idx ] .*
        exp.( im * x0_θh_view[ gens_idx ] )
    
    x0_vh_view[gens_idx] .= gens_vh
    x0_θh_view[gens_idx] .= angle.( gens_uh )

    red_vh_θh_0 = [ red_vh_θh_0_view; ]

    #------------------------------------------    
    
    Jac_row_size =
        Jac_col_size = length( red_vh_θh_idx )
    
    Jac_vh_θh =
        zeros( Jac_row_size, Jac_col_size )

    #------------------------------------------    

    dynamic_power_balance_Jac_with_vh_θh!(
        Jac_vh_θh,
        red_vh_θh_0,
        working_vh_θh_view,
        pf_net_param )

    dynamic_power_balance_mismatch_with_vh_θh!(
        mismatch,
        red_vh_θh_0,
        (working_vh_θh_view,
         nodes_pf_U_view,
         Inet_view,
         Iinj_view,
         δ_ω_ed_dash_eq_dash_view),
        pf_net_param;
        with_δ_ed_eq = with_δ_ed_eq )
    
    #------------------------------------------    

    sol = nlsolve(
        (mismatch, red) ->
            dynamic_power_balance_mismatch_with_vh_θh!(
                mismatch,
                red,
                (working_vh_θh_view,
                 nodes_pf_U_view,
                 Inet_view,
                 Iinj_view,
                 δ_ω_ed_dash_eq_dash_view),
                pf_net_param;
                with_δ_ed_eq = with_δ_ed_eq),
        (Jac_vh_θh, red )->
            dynamic_power_balance_Jac_with_vh_θh!(
                Jac_vh_θh,
                red,
                working_vh_θh_view,
                pf_net_param),
        red_vh_θh_0;
        ftol = ftol )

    return sol
    
end

#-----------------------------------------------------

function power_balance_powerflow_with_vh_θh!(
    x0_vh_θh,
    mismatch,
    (working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view),
    pf_net_param;
    maxiter=1e5,
    ftol=1000*eps(),
    xtol=1000*eps(),
    with_δ_ed_eq = false)

    # -----------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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
    
    # ------------------------------------------------   
    
    x0_vh_view    = @view x0_vh_θh[vh_IDX]
    x0_θh_view    = @view x0_vh_θh[θh_IDX]
 
    # Get a view of a reduced vh_θh
    
    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx ]

    # -----------------------------------------------   

    gens_uh = x0_vh_view[gens_idx] .*
        exp.(im * x0_θh_view[gens_idx])
    
    x0_vh_view[gens_idx] .= gens_vh
    x0_θh_view[gens_idx] .= angle.( gens_uh )

    red_vh_θh_0 = [ red_vh_θh_0_view;]

    #------------------------------------------    
    #------------------------------------------    
    
    Jac_row_size = Jac_col_size =
        length( red_vh_θh_idx )
    
    Jac_vh_θh = zeros(
        Jac_row_size,
        Jac_col_size)

    power_balance_Jac_with_vh_θh!(
        Jac_vh_θh,
        red_vh_θh_0,
        working_vh_θh_view,
        pf_net_param ) 

    power_balance_mismatch_with_vh_θh!(
        mismatch,
        red_vh_θh_0,
        (working_vh_θh_view,
         nodes_pf_U_view,
         Inet_view,
         Iinj_view,
         δ_ω_ed_dash_eq_dash_view),
        pf_net_param;
        with_δ_ed_eq = with_δ_ed_eq )

    #------------------------------------------    

    sol = nlsolve(
        ( mismatch,
          red ) ->  power_balance_mismatch_with_vh_θh!(
              mismatch,
              red,
              (working_vh_θh_view,
               nodes_pf_U_view,
               Inet_view,
               Iinj_view,
               δ_ω_ed_dash_eq_dash_view ),
              pf_net_param;
              with_δ_ed_eq = with_δ_ed_eq),
        (Jac_vh_θh, red)->
            power_balance_Jac_with_vh_θh!(
                Jac_vh_θh,
                red,
                working_vh_θh_view,
                pf_net_param),
        red_vh_θh_0;
        ftol =  ftol )

    return sol # sol.zero
    
end

# ------------------------------------------------------
# ------------------------------------------------------

function dynamic_power_balance_powerflow_from_state_with_vh_θh!( x,  mismatch, (working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view), pf_param ; maxiter=1e5, ftol=1000*eps(),  xtol=1000*eps(), with_δ_ed_eq = true, alg = SimpleNewtonRaphson()  )

     _, _, _, _, pf_Idx = pf_param
     
    (slack_bus_idx,
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

    #------------------------------------------------
    
    uh_state_x0 = x[ur_ui_idx][ ur_IDX ] +
        im * x[ur_ui_idx][ ui_IDX ]
            
    x0_vh_θh = [abs.(uh_state_x0)...;
                angle.(uh_state_x0)...]
    
    #------------------------------------------------ 

    dynamic_power_balance_powerflow_with_vh_θh!(
        x0_vh_θh,
        mismatch,
        (working_vh_θh_view,
         nodes_pf_U_view,
         Inet_view,
         Iinj_view,
         δ_ω_ed_dash_eq_dash_view),
        pf_net_param;
        maxiter = 1e5,
        ftol = ftol,
        xtol = xtol,
        with_δ_ed_eq = with_δ_ed_eq, alg = alg  )

    

    # return nodes_pf_U_view, Inet_view

    return nothing
    
end

# ------------------------------------------------------
# New
# ------------------------------------------------------

function dynamic_power_balance_powerflow(
    x0_vh_θh,
    mismatch,
    (working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view),
    pf_net_param;
    maxiter=1e5,
    ftol=1000*eps(),
    xtol=1000*eps(),
    with_δ_ed_eq = true)

    # ---------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
     nodes_u_Idx,
     gens_idx,
     ur_IDX,
     ui_IDX,
     vh_IDX,
     θh_IDX,
     red_vh_θh_idx,
     ur_idx,ui_idx,
     ur_ui_idx) =
         pf_Idx
    
    # ---------------------------------------------------
    # solve powerflow 
    # ---------------------------------------------------

    sol_with_vh_θh =
        dynamic_power_balance_powerflow_with_vh_θh!(
            x0_vh_θh,
            mismatch,
            ( working_vh_θh_view,
              nodes_pf_U_view,
              Inet_view,
              Iinj_view,
              δ_ω_ed_dash_eq_dash_view ),
            pf_net_param ;
            maxiter=maxiter,
            ftol=ftol,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )

    
    # ---------------------------------------------------
    # Results 
    # ---------------------------------------------------
    
    red_vh_θh = sol_with_vh_θh.zero

    working_vh_θh_view[ red_vh_θh_idx ] .= red_vh_θh

    uh =  working_vh_θh_view[ vh_IDX ] .*
        exp.(im *  working_vh_θh_view[ θh_IDX ])

    # --------------------------------------------------

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    #------------------------------------------

    Inet_view .= get_Inet_inj(
        uh,
        Ynet,
        nodes_node_idx_and_incident_edges_other_node_idx)
    
    Iinj_view .=  get_Iinj(
        uh,
        P_Q_non_gens_view,
        P_Q_gens_loc_load_view,
        Inet_view  )

    
    # -------------------------------------------------- 
    # update nodes_pf_U_view
    # -------------------------------------------------- 

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .=
            [ real( uk ), imag( uk ) ]

    end

    # ------------------------------------------------- 
    
    return  nothing

end


function power_balance_powerflow(
    x0_vh_θh,
    mismatch,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view ),
    (nodes_name,
     branches_name),
    pf_net_param;
    maxiter=1e5,
    ftol=1000*eps(),
    xtol=1000*eps(),
    with_δ_ed_eq = false )

    # ---------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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
    
    # ---------------------------------------------------
    # solve powerflow 
    # ---------------------------------------------------

    sol_with_vh_θh =
        power_balance_powerflow_with_vh_θh!(
            x0_vh_θh,
            mismatch,
            (working_vh_θh_view,
             nodes_pf_U_view,
             Inet_view,
             Iinj_view,
             δ_ω_ed_dash_eq_dash_view),
            pf_net_param;
            with_δ_ed_eq =
                with_δ_ed_eq )
    
    # ---------------------------------------------------
    # Results 
    # ---------------------------------------------------
    
    red_vh_θh = sol_with_vh_θh.zero

    working_vh_θh_view[ red_vh_θh_idx ] .= red_vh_θh

    uh =
        working_vh_θh_view[ vh_IDX ] .*
        exp.( im *  working_vh_θh_view[ θh_IDX ])

    # --------------------------------------------------

    S_gens = x_from_xr_xi.(P_Q_gens_view)

    S_gens_loc_load =
        x_from_xr_xi.(
            P_Q_gens_loc_load_view)
    
    S_non_gens = x_from_xr_xi.(
        P_Q_non_gens_view)

    #------------------------------------------

    Inet_view .= get_Inet_inj(
        uh,
        Ynet,
        nodes_node_idx_and_incident_edges_other_node_idx)
    
    Iinj_view .= get_Iinj(
        uh,
        P_Q_non_gens_view,
        P_Q_gens_loc_load_view,
        Inet_view )

    Sbus = uh .* conj.( Inet_view )

    GenSinj = Sbus + x_from_xr_xi.(P_Q_non_gens_view) +
        x_from_xr_xi.( P_Q_gens_loc_load_view )
    
    # Inet_inj = get_Inet_inj(
    #     uh, Ynet,
    #   nodes_node_idx_and_incident_edges_other_node_idx)
    
    # Iinj = get_Iinj(
    #     uh,
    #     S_non_gens,
    #     S_gens_loc_load,
    #     Inet_inj)
    
    # Sbus_n  = uh .* conj.( Inet_inj )
    # GenSinj = Sbus + S_non_gens + S_gens_loc_load
        
    Ifrom_Ito = get_Ifrom_Ito(
        nodes_pf_U_view,
        edges_Ybr_cal,
        edges_orientation)

    If = first.(Ifrom_Ito)
    
    It = last.(Ifrom_Ito)

    Ibranches = If + It
    
    #------------------------------------------

    bus_dict_Iinj = OrderedDict(
        name => [ih_r, ih_i]
        for (name, ih_r, ih_i) in
            zip(nodes_name,
                real.( Iinj_view ),
                imag.( Iinj_view )))

    bus_dict_Inet_inj = OrderedDict(
        name => [ih_r, ih_i]
        for (name, ih_r, ih_i) in
            zip(nodes_name,
                real.(Inet_view ),
                imag.(  Inet_view )))     

    branch_dict_init = OrderedDict(
        name => (real(i_f),
                 imag(i_f),
                 real(i_t),
                 imag(i_t))
        for (name, i_f, i_t) in
            zip(branches_name, If, It))

    bus_dict_init = OrderedDict(
        name => (vh, θh, ph, qh, ih_r,
                 ih_i, pg, qg, ig_r, ig_i)
        for (name, vh, θh, ph, qh, ih_r,
             ih_i, pg, qg, ig_r, ig_i ) in
            zip(nodes_name,
                abs.( uh),
                angle.( uh),
                real.( Sbus),
                imag.( Sbus),
                real.( Inet_view ),
                imag.( Inet_view ),
                real.( GenSinj),
                imag.( GenSinj),
                real.( Iinj_view ),
                imag.( Iinj_view )) )

    #------------------------------------------

    Vm   = abs.(uh)
    Vθ   = angle.(uh)
    Vbus = uh
    
    #------------------------------------------

    return  (; S_gens,
             S_non_gens,
             S_gens_loc_load,
             Vm,
             Vθ,
             Vbus,
             Ibranches,
             Ybus ,
             Sbus,
             GenSinj,
             bus_dict_init,
             branch_dict_init,
             bus_dict_Iinj,
             bus_dict_Inet_inj,
             If,
             It) 

end

#---------------------
# Dynamic powerflow function for nodes_edges_dynamics!
# --------------------

# ---------------------------------------------------
# Begin im model pf
# ---------------------------------------------------

function im_model_dynamic_power_balance_powerflow_with_vh_θh!(
    x0_vh_θh,
    mismatch,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view ),
    pf_net_param,
    im_model_idq_pf_cal;
    maxiter=1e5,
    ftol=1000*eps(),
    xtol=1000*eps(),
    with_δ_ed_eq = true )

    # ---------------------------------------------------

     pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx, pf_net_misc = pf_net_param
    
    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net
    
    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx
    
    #------------------------------------------    
    
    x0_vh_view =
        @view x0_vh_θh[vh_IDX]
    x0_θh_view =
        @view x0_vh_θh[θh_IDX]

    # Get a view of a reduced vh_θh
    
    red_vh_θh_0_view =
        @view x0_vh_θh[ red_vh_θh_idx  ]
    
    #------------------------------------------    
    
    gens_uh =
        x0_vh_view[ gens_idx ] .*
        exp.( im * x0_θh_view[ gens_idx ] )
    
    x0_vh_view[gens_idx] .=
        gens_vh
    
    x0_θh_view[gens_idx] .=
        angle.( gens_uh )

    red_vh_θh_0 =
        [ red_vh_θh_0_view; ]

    #------------------------------------------    
    
    Jac_row_size = Jac_col_size = length( red_vh_θh_idx )
    
    Jac_vh_θh =
        zeros( Jac_row_size, Jac_col_size )

    #------------------------------------------    

    dynamic_power_balance_Jac_with_vh_θh!(
        Jac_vh_θh, red_vh_θh_0,
        working_vh_θh_view, pf_net_param )

    im_model_dynamic_power_balance_mismatch_with_vh_θh!(
        mismatch, red_vh_θh_0,
        ( working_vh_θh_view, nodes_pf_U_view,
          Inet_view, Iinj_view,
          δ_ω_ed_dash_eq_dash_view ),
        pf_net_param, im_model_idq_pf_cal;
        with_δ_ed_eq = with_δ_ed_eq )
    
    #------------------------------------------    

    sol =
        nlsolve(
            ( mismatch, red ) -> im_model_dynamic_power_balance_mismatch_with_vh_θh!(
                mismatch, red,
                ( working_vh_θh_view,
                  nodes_pf_U_view,
                  Inet_view,
                  Iinj_view,
                  δ_ω_ed_dash_eq_dash_view ),
                pf_net_param, im_model_idq_pf_cal;
                with_δ_ed_eq = with_δ_ed_eq ),
            ( Jac_vh_θh, red )-> dynamic_power_balance_Jac_with_vh_θh!(
                Jac_vh_θh,
                red,
                working_vh_θh_view,
                pf_net_param ),
            red_vh_θh_0;
            ftol =  ftol )

    return sol
    
end



function im_model_power_balance_powerflow_with_vh_θh!(
    x0_vh_θh, mismatch,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view  ),
    pf_net_param,
    im_model_idq_pf_cal;
    maxiter=1e5,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    with_δ_ed_eq = false  )

    # ------------------------------------------------------

     pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param
    
    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net
    
    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx
    
    # ------------------------------------------------------   
    
    x0_vh_view    = @view x0_vh_θh[vh_IDX]
    x0_θh_view    = @view x0_vh_θh[θh_IDX]
 
    # Get a view of a reduced vh_θh
    
    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]

    # ------------------------------------------------------   

    gens_uh =
        x0_vh_view[gens_idx] .*
        exp.(im * x0_θh_view[gens_idx])
    
    x0_vh_view[gens_idx] .= gens_vh
    x0_θh_view[gens_idx] .= angle.( gens_uh )

    red_vh_θh_0 = [ red_vh_θh_0_view;]

    #------------------------------------------    
    #------------------------------------------    
    
    Jac_row_size = Jac_col_size = length( red_vh_θh_idx )
    
    Jac_vh_θh = zeros(Jac_row_size, Jac_col_size)

    power_balance_Jac_with_vh_θh!(
        Jac_vh_θh,
        red_vh_θh_0,
        working_vh_θh_view,
        pf_net_param ) 

    im_model_power_balance_mismatch_with_vh_θh!(
        mismatch,
        red_vh_θh_0,
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          δ_ω_ed_dash_eq_dash_view ),
        pf_net_param,
        im_model_idq_pf_cal;
        with_δ_ed_eq = with_δ_ed_eq )

    #------------------------------------------    

    sol = nlsolve(
        ( mismatch, red ) ->  im_model_power_balance_mismatch_with_vh_θh!(
            mismatch,
            red,
            ( working_vh_θh_view,
              nodes_pf_U_view,
              Inet_view,
              Iinj_view,
              δ_ω_ed_dash_eq_dash_view ),
            pf_net_param,
            im_model_idq_pf_cal;
            with_δ_ed_eq = with_δ_ed_eq ) ,
        (Jac_vh_θh, red)-> power_balance_Jac_with_vh_θh!(
            Jac_vh_θh,
            red,
            working_vh_θh_view,
            pf_net_param ),
        red_vh_θh_0;
        ftol =  ftol )

    return sol # sol.zero
    
end


function im_model_power_balance_powerflow(
    x0_vh_θh,
    mismatch,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view ),
    pf_net_param,
    im_model_idq_pf_cal;
    maxiter=1e5,
    ftol=1000*eps(),
    xtol=1000*eps(),
    with_δ_ed_eq = true )

    # ---------------------------------------------------

     pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param
    
    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net
    
    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx
    
    # ---------------------------------------------------
    # solve powerflow 
    # ---------------------------------------------------

    sol_with_vh_θh =
        im_model_power_balance_powerflow_with_vh_θh!(
            x0_vh_θh,
            mismatch,
            ( working_vh_θh_view,
              nodes_pf_U_view,
              Inet_view,
              Iinj_view,
              δ_ω_ed_dash_eq_dash_view ),
            pf_net_param,
            im_model_idq_pf_cal;
            with_δ_ed_eq = with_δ_ed_eq )
    
    # ---------------------------------------------------
    # Results 
    # ---------------------------------------------------
    
    red_vh_θh =
        sol_with_vh_θh.zero

    working_vh_θh_view[ red_vh_θh_idx ] .= red_vh_θh

    uh =
        working_vh_θh_view[ vh_IDX ] .*
        exp.(im *  working_vh_θh_view[ θh_IDX ])

    # --------------------------------------------------

    S_gens = x_from_xr_xi.(  P_Q_gens_view )

    # --------------------------------------------------
    
    return nothing
end


function im_model_dynamic_power_balance_powerflow(
    x0_vh_θh,
    mismatch,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view ),
    pf_net_param,
    im_model_idq_pf_cal;
    maxiter=1e5,
    ftol=1000*eps(),
    xtol=1000*eps() ,
    with_δ_ed_eq = true )

    # ---------------------------------------------------

     pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param
    
    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net
    
    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx
    
    # ---------------------------------------------------
    # solve powerflow 
    # ---------------------------------------------------

    sol_with_vh_θh =
        im_model_dynamic_power_balance_powerflow_with_vh_θh!(
            x0_vh_θh,
            mismatch,
            ( working_vh_θh_view,
              nodes_pf_U_view,
              Inet_view,
              Iinj_view,
              δ_ω_ed_dash_eq_dash_view ),
            pf_net_param,
            im_model_idq_pf_cal;
            maxiter=maxiter,
            ftol=ftol,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )

    
    # ---------------------------------------------------
    # Results 
    # ---------------------------------------------------
    
    red_vh_θh =
        sol_with_vh_θh.zero

    working_vh_θh_view[ red_vh_θh_idx ] .= red_vh_θh

    uh =
        working_vh_θh_view[ vh_IDX ] .*
        exp.(im *  working_vh_θh_view[ θh_IDX ])

    # --------------------------------------------------

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    #------------------------------------------

    Inet_view .=
        get_Inet_inj(
            uh,
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )
    
    Iinj_view .=
        get_Iinj(
            uh,
            P_Q_non_gens_view,
            P_Q_gens_loc_load_view,
            Inet_view  )
    
    # ---------------------------------------------------
    # update nodes_pf_U_view
    # ---------------------------------------------------

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .= [ real( uk ), imag( uk ) ]

    end

    return  nothing

end


# -------
# End im model pf
# -------

#-------------------------------------------------------

# -------
# Begin industrial model pf
# -------

function industrial_model_dynamic_power_balance_powerflow_with_vh_θh!( x0_vh_θh, mismatch, ( working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view ), pf_net_param, industrial_model_idq_pf_cal ; maxiter=1e5, ftol=1000*eps(), xtol=1000*eps(), with_δ_ed_eq = true )

    # ------------------------------------------------------

     pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param
    
    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net
    
    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx
    
    #------------------------------------------    
    
    x0_vh_view    = @view x0_vh_θh[vh_IDX]
    x0_θh_view    = @view x0_vh_θh[θh_IDX]

    # Get a view of a reduced vh_θh
    
    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]
    
    #------------------------------------------    
    
    gens_uh = x0_vh_view[ gens_idx ] .*
        exp.( im * x0_θh_view[ gens_idx ] )
    
    x0_vh_view[gens_idx] .= gens_vh
    x0_θh_view[gens_idx] .= angle.( gens_uh )

    red_vh_θh_0 = [ red_vh_θh_0_view; ]

    #------------------------------------------    
    
    Jac_row_size = Jac_col_size = length( red_vh_θh_idx )
    
    Jac_vh_θh = zeros( Jac_row_size, Jac_col_size )

    #------------------------------------------    

    dynamic_power_balance_Jac_with_vh_θh!(
        Jac_vh_θh,
        red_vh_θh_0,
        working_vh_θh_view,
        pf_net_param )

    industrial_model_dynamic_power_balance_mismatch_with_vh_θh!(
        mismatch,
        red_vh_θh_0,
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          δ_ω_ed_dash_eq_dash_view ),
        pf_net_param,
        industrial_model_idq_pf_cal ;
        with_δ_ed_eq = with_δ_ed_eq )
    
    #------------------------------------------    

    sol = nlsolve(
        ( mismatch, red ) -> industrial_model_dynamic_power_balance_mismatch_with_vh_θh!(
            mismatch, red,
            ( working_vh_θh_view,
              nodes_pf_U_view,
              Inet_view,
              Iinj_view,
              δ_ω_ed_dash_eq_dash_view ),
            pf_net_param,
            industrial_model_idq_pf_cal;
            with_δ_ed_eq = with_δ_ed_eq ),
        ( Jac_vh_θh, red ) ->
            dynamic_power_balance_Jac_with_vh_θh!(
                Jac_vh_θh, red,
                working_vh_θh_view,
                pf_net_param ),
        red_vh_θh_0;
        ftol =  ftol )

    return sol
    
end



function industrial_model_power_balance_powerflow_with_vh_θh!(
    x0_vh_θh, mismatch,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view, Iinj_view,
      δ_ω_ed_dash_eq_dash_view),
    pf_net_param,
    industrial_model_idq_pf_cal;
    maxiter=1e5,
    ftol=1000*eps(),
    xtol=1000*eps(),
    with_δ_ed_eq = false  )

    # --------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) = pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net
    
    (slack_vh, gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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
    
    # ------------------------------------------------------   
    
    x0_vh_view = @view x0_vh_θh[vh_IDX]
    x0_θh_view = @view x0_vh_θh[θh_IDX]
 
    # Get a view of a reduced vh_θh
    
    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]

    # ------------------------------------------------------   

    gens_uh = x0_vh_view[gens_idx] .*
        exp.(im * x0_θh_view[gens_idx])
    
    x0_vh_view[gens_idx] .= gens_vh
    x0_θh_view[gens_idx] .= angle.( gens_uh )

    red_vh_θh_0 = [ red_vh_θh_0_view;]

    #------------------------------------------    
    #------------------------------------------    
    
    Jac_row_size = Jac_col_size = length( red_vh_θh_idx )
    
    Jac_vh_θh = zeros(Jac_row_size, Jac_col_size)

    power_balance_Jac_with_vh_θh!(
        Jac_vh_θh,
        red_vh_θh_0,
        working_vh_θh_view,
        pf_net_param ) 

    industrial_model_power_balance_mismatch_with_vh_θh!(
        mismatch,
        red_vh_θh_0,
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          δ_ω_ed_dash_eq_dash_view ),
        pf_net_param,
        industrial_model_idq_pf_cal;
        with_δ_ed_eq = with_δ_ed_eq )

    #------------------------------------------    

    sol = nlsolve(
        ( mismatch, red ) ->  industrial_model_power_balance_mismatch_with_vh_θh!(
            mismatch, red,
            ( working_vh_θh_view,
              nodes_pf_U_view,
              Inet_view,
              Iinj_view,
              δ_ω_ed_dash_eq_dash_view ),
            pf_net_param,
            industrial_model_idq_pf_cal;
            with_δ_ed_eq = with_δ_ed_eq ) ,
        (Jac_vh_θh, red)->
            power_balance_Jac_with_vh_θh!(
                Jac_vh_θh, red,
                working_vh_θh_view,
                pf_net_param ),
        red_vh_θh_0;
        ftol =  ftol )

    return sol # sol.zero
    
end


function industrial_model_power_balance_powerflow(
    x0_vh_θh,
    mismatch,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view,
      Iinj_view,
      δ_ω_ed_dash_eq_dash_view ),
    pf_net_param,
    industrial_model_idq_pf_cal;
    maxiter=1e5,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    with_δ_ed_eq = true )

    # ---------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) = pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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
    
    # ---------------------------------------------------
    # solve powerflow 
    # ---------------------------------------------------

    sol_with_vh_θh =
        industrial_model_power_balance_powerflow_with_vh_θh!(x0_vh_θh, mismatch, ( working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view ), pf_net_param, industrial_model_idq_pf_cal ; with_δ_ed_eq = with_δ_ed_eq )
    
    # ---------------------------------------------------
    # Results 
    # ---------------------------------------------------
    
    red_vh_θh = sol_with_vh_θh.zero

    working_vh_θh_view[ red_vh_θh_idx ] .= red_vh_θh

    uh =  working_vh_θh_view[ vh_IDX ] .*
        exp.(im *  working_vh_θh_view[ θh_IDX ])

    # --------------------------------------------------

    S_gens = x_from_xr_xi.(  P_Q_gens_view )

    # --------------------------------------------------
    
    return nothing
end


function industrial_model_dynamic_power_balance_powerflow(
    x0_vh_θh, mismatch,
    ( working_vh_θh_view,
      nodes_pf_U_view,
      Inet_view, Iinj_view,
      δ_ω_ed_dash_eq_dash_view ),
    pf_net_param,
    industrial_model_idq_pf_cal ;
    maxiter=1e5, ftol=1000*eps() ,
    xtol=1000*eps() ,
    with_δ_ed_eq = true )

    # ---------------------------------------------------

    (pf_net, pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) = pf_net_param
    
    (Ybus, Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net
    
    (slack_vh, gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
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
    
    # ---------------------------------------------------
    # solve powerflow 
    # ---------------------------------------------------

    sol_with_vh_θh =
        industrial_model_dynamic_power_balance_powerflow_with_vh_θh!(
            x0_vh_θh, mismatch,
            ( working_vh_θh_view,
              nodes_pf_U_view,
              Inet_view,
              Iinj_view,
              δ_ω_ed_dash_eq_dash_view ),
            pf_net_param,
            industrial_model_idq_pf_cal;
            maxiter=maxiter,
            ftol=ftol,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )

    
    # -------------------------------------------------
    # Results 
    # -------------------------------------------------
    
    red_vh_θh = sol_with_vh_θh.zero

    working_vh_θh_view[ red_vh_θh_idx ] .= red_vh_θh

    uh =
        working_vh_θh_view[ vh_IDX ] .*
        exp.(im *  working_vh_θh_view[ θh_IDX ])

    # -----------------------------------------

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    #------------------------------------------

    Inet_view .=
        get_Inet_inj(
            uh, Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx )
    
    Iinj_view .=
        get_Iinj(
            uh,
            P_Q_non_gens_view,
            P_Q_gens_loc_load_view,
            Inet_view  )
    
    # ---------------------------------------------- 
    # update nodes_pf_U_view
    # ---------------------------------------------- 

    for (k, uk) in enumerate(uh)

        nodes_pf_U_view[ k ] .=
            [ real( uk ), imag( uk ) ]

    end

    # ------------------------------------------------ 
    

    return  nothing

end

# -------
# End industrial model pf
# -------

# -------------------------------------------------------
# ------------------------------------------------------

function get_power_balance_powerflow_results(
    x0_vh_θh,
    pf_param,
    nodes_name,
    branches_name
    ;with_δ_ed_eq =
        false  )

    # ---------------------------------------------------

    pf_net, pf_idx_and_state, pf_views, pf_limits, pf_Idx = pf_param
     
    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net
    
     slack_vh,  gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state,  non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    working_vh_θh_view, nodes_u_view, nodes_pf_U_view,  nodes_idx_and_δ_ω_ed_dash_eq_dash_view, ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view, Inet_view, Iinj_view =  pf_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx,  ur_idx,  ui_idx, ur_ui_idx  = pf_Idx
    
    # -----------------------------------------------
    # -----------------------------------------------
    # Results 
    # -----------------------------------------------
    # -----------------------------------------------  
    
    sol_with_vh_θh =
        power_balance_powerflow_with_vh_θh!(
            x0_vh_θh, pf_param
            ; with_δ_ed_eq =
                with_δ_ed_eq )
    
    red_vh_θh =
        sol_with_vh_θh.zero

    working_vh_θh_view[ red_vh_θh_idx ] .= red_vh_θh

    uh =
        working_vh_θh_view[ vh_IDX ] .* exp.(
            im *  working_vh_θh_view[ θh_IDX ])


    # ---------------------------------------------
    # ---------------------------------------------

    S_gens =
        x_from_xr_xi.(
            P_Q_gens_view)

    S_gens_loc_load =
        x_from_xr_xi.(
            P_Q_gens_loc_load_view)

    S_non_gens =
        x_from_xr_xi.(
            P_Q_non_gens_view)

    # ----------------------------------------------    
    
    Inet_inj = [
        sum([ ynj * vj
              for (ynj, vj) in
                  zip( Y_bus_vec,
                       uh[nth_node_idx_and_adj_nodes_idx] )] )
        for ( Y_bus_vec,
              nth_node_idx_and_adj_nodes_idx ) in
            
            zip( Ynet,
                 nodes_node_idx_and_incident_edges_other_node_idx) ]
    
    Iinj = Inet_inj +
        (conj.( S_non_gens ))./ ( conj.( uh )) +
        (conj.( S_gens_loc_load ))./ (conj.( uh )) 

    Sbus_n  =
        uh .* conj.( Inet_inj )

    GenSinj =
        Sbus_n +
        S_non_gens +
        S_gens_loc_load

    Ifrom_Ito = [
        y_π * [ x_from_xr_xi(
            nodes_pf_U_view[ orient[1] ]),
                x_from_xr_xi(
                    nodes_pf_U_view[ orient[2] ])]
        for (y_π, orient ) in
            zip(edges_Ybr_cal,
                edges_orientation)]

    If = first.(Ifrom_Ito)
    It = last.(Ifrom_Ito)

    Ibranches = If + It
    
    # -------------------------------------------

    bus_dict_Iinj = OrderedDict(
        name => [ih_r, ih_i]
        for ( name, ih_r, ih_i) in
            zip( nodes_name,
                 real.(Iinj),
                 imag.(Iinj) ))

    bus_dict_Inet_inj =
        OrderedDict(
            name => [ih_r, ih_i]
            for (name, ih_r, ih_i) in
                zip( nodes_name,
                     real.(Inet_inj),
                     imag.(Inet_inj))) 

    branch_dict_init = OrderedDict(
        name => (
            real(i_f),
            imag(i_f),
            real(i_t),
            imag(i_t) )
        for (name, i_f, i_t) in
            zip( branches_name,
                 If,
                 It ) )

    bus_dict_init = OrderedDict(
        name => (vh,
                 θh,
                 ph,
                 qh,
                 ih_r,
                 ih_i,
                 pg,
                 qg,
                 ig_r,
                 ig_i )
        for (name,
             vh,
             θh,
             ph,
             qh,
             ih_r,
             ih_i,
             pg,
             qg,
             ig_r,
             ig_i ) in
            zip(nodes_name,  abs.(uh),
                angle.(uh),
                real.(Sbus_n),
                imag.(Sbus_n),
                real.(Inet_inj),
                imag.(Inet_inj),
                real.(GenSinj),
                imag.(GenSinj),
                real.(Iinj),
                imag.(Iinj)) )

    return Dict(
        "Iinj" => Iinj,
        "Inet_inj" => Inet_inj,
        "Sg" => S_gens,
        "Sd" => S_non_gens,
        "Sloc" => S_gens_loc_load,
        "Vm" => abs.(uh),
        "Vθ" => angle.(uh),
        "Vbus" => uh,
        "Ibranches" => Ibranches,
        "Ybus" => Ybus ,
        "Sbus" => Sbus_n,
        "GenSinj" => GenSinj,
        "dict_init" =>Dict(
            "bus_dict_init" => bus_dict_init,
            "branch_dict_init" =>  branch_dict_init,
            "bus_dict_Iinj" => bus_dict_Iinj,
            "bus_dict_Inet_inj" => bus_dict_Inet_inj ),
        "Ifrom" => If,
        "Ito" => It,
        "Ynet"=> Ynet )

end


function get_power_balance_powerflow_results(
    x0_vh_θh,
    (pf_param,
     nodes_name,
     branches_name )
    ; with_δ_ed_eq =
        false  )

    # ---------------------------------------------------

    pf_net, pf_idx_and_state, pf_views, pf_limits, pf_Idx = pf_param
     
    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net
    
     slack_vh,  gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state,  non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    working_vh_θh_view, nodes_u_view, nodes_pf_U_view, δ_ω_ed_dash_eq_dash_view, ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view, Inet_view, Iinj_view =  pf_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx,  ur_idx,  ui_idx, ur_ui_idx  = pf_Idx
    
    
    # ---------------------------------------------------
    # ---------------------------------------------------
    # Results 
    # ---------------------------------------------------
    # ---------------------------------------------------
    
    sol_with_vh_θh =
        power_balance_powerflow_with_vh_θh!(
            x0_vh_θh, pf_param
            ; with_δ_ed_eq =
                with_δ_ed_eq )
    
    red_vh_θh =
        sol_with_vh_θh.zero

    working_vh_θh_view[ red_vh_θh_idx ] .= red_vh_θh

    uh =
        working_vh_θh_view[ vh_IDX ] .*
        exp.(im * working_vh_θh_view[ θh_IDX ])

    # -----------------------------------------------

    S_gens =
        x_from_xr_xi.(
            P_Q_gens_view)

    S_gens_loc_load =
        x_from_xr_xi.(
            P_Q_gens_loc_load_view)

    S_non_gens =
        x_from_xr_xi.(
            P_Q_non_gens_view)

    # -----------------------------------------------    
    
    Inet_inj = [
        sum([ ynj * vj
              for (ynj, vj) in
                  zip(Y_bus_vec,
                      uh[nth_node_idx_and_adj_nodes_idx]) ])
        for (Y_bus_vec,
             nth_node_idx_and_adj_nodes_idx ) in
            zip(Ynet,
                nodes_node_idx_and_incident_edges_other_node_idx ) ]
    
    Iinj =
        Inet_inj + (conj.( S_non_gens ))./ ( conj.( uh )) +
        ( conj.( S_gens_loc_load ))./ (conj.( uh )) 

    Sbus_n  =
        uh .* conj.( Inet_inj )

    GenSinj =
        Sbus_n +
        S_non_gens +
        S_gens_loc_load

    Ifrom_Ito = [
        y_π * [ x_from_xr_xi( nodes_pf_U_view[ orient[1] ]),
                x_from_xr_xi( nodes_pf_U_view[ orient[2] ])]
        for (y_π, orient ) in
            zip(edges_Ybr_cal,
                edges_orientation)]

    If = first.( Ifrom_Ito)
    It = last.( Ifrom_Ito)

    Ibranches = If + It
    
    # ----------------------------------------------------

    bus_dict_Iinj = OrderedDict(
        name => [ih_r, ih_i]
        for (name, ih_r, ih_i) in
            zip( nodes_name,
                 real.(Iinj),
                 imag.(Iinj)))

    bus_dict_Inet_inj = OrderedDict(
        name => [ih_r, ih_i]
        for (name, ih_r, ih_i) in
            zip( nodes_name,
                 real.(Inet_inj),
                 imag.(Inet_inj))) 

    branch_dict_init = OrderedDict(
        name => ( real(i_f), imag(i_f), real(i_t), imag(i_t) )
        for (name, i_f, i_t) in
            zip(branches_name, If, It))

    bus_dict_init = OrderedDict(
        name => (vh,
                 θh,
                 ph,
                 qh,
                 ih_r,
                 ih_i,
                 pg,
                 qg,
                 ig_r,
                 ig_i )
        for (name,
             vh,
             θh,
             ph,
             qh,
             ih_r,
             ih_i,
             pg,
             qg,
             ig_r,
             ig_i ) in
            zip(nodes_name,
                abs.(uh),
                angle.(uh),
                real.(Sbus_n),
                imag.(Sbus_n),
                real.(Inet_inj),
                imag.(Inet_inj),
                real.(GenSinj),
                imag.(GenSinj),
                real.(Iinj),
                imag.(Iinj)) )

    return Dict(
        "Iinj" => Iinj,
        "Inet_inj" => Inet_inj,
        "Sg" => S_gens,
        "Sd" => S_non_gens,
        "Sloc" => S_gens_loc_load,
        "Vm" => abs.(uh),
        "Vθ" => angle.(uh),
        "Vbus" => uh,
        "Ibranches" => Ibranches,
        "Ybus" => Ybus ,
        "Sbus" => Sbus_n,
        "GenSinj" => GenSinj,
        "dict_init" => Dict(
            "bus_dict_init" => bus_dict_init,
            "branch_dict_init" => branch_dict_init,
            "bus_dict_Iinj" => bus_dict_Iinj,
            "bus_dict_Inet_inj" => bus_dict_Inet_inj),
        "Ifrom" => If,
        "Ito" => It,
        "Ynet"=> Ynet )

end


# ------------------------------------------------------
# ------------------------------------------------------


function intra_dyn_pf_by_vh_θh_by_parts_func!(
    intra_dyn_vh_θh, flat_δ_ed_dash_eq_dash,
    dyn_pf_flat_para ;
    intra_dyn_pf_by_vh_θh_by_parts_kwd =
        intra_dyn_pf_by_vh_θh_by_parts_kwd)

    (;
     vh_θh_DiffCache,
     Pg_Qg_external_control,
     intra_pf_kwd_para,
     pf_solver,
     use_nlsolve,
     post_pf_idxs,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs ) =
         intra_dyn_pf_by_vh_θh_by_parts_kwd 
    
    #----------------------------------------    

    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver

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
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (
     loc_load_exist,
     vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         intra_dyn_pf_mismatch_kwd_para

    #----------------------------------------
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        vars_Idxs
    
    #----------------------------------------    
    
  (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
  
    #----------------------------------------    

    (;     
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) = 
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs

    #----------------------------------------    

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 

    #----------------------------------------    
    
    (; Pg_Idx,
      Qg_Idxs,
      Png_Idxs,
      Qng_Idxs,
      Pgll_Idxs,
      Qgll_Idxs ) = Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------    

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para

    #----------------------------------------     
    #----------------------------------------    
    
    P_gens =
        dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        intra_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        intra_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        intra_dyn_pf_flat_para[ Qng_Idxs ]

    # flat_δ_ed_dash_eq_dash =
    #     intra_dyn_pf_flat_para[ dyn_δ_ed_eq_pf_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            intra_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            intra_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #----------------------------------------      
    #----------------------------------------

    # if use_nlsolve == true

    #     vh = intra_dyn_vh_θh[ vh_Idx ]
    #     θh = intra_dyn_vh_θh[ θh_Idx ]
        
    # else
    
    #     vh_θh_DiffCache =
    #         get_tmp( vh_θh_DiffCache,
    #                 intra_dyn_vh_θh )
    #     vh_θh_DiffCache .= intra_dyn_vh_θh
    #     vh = vh_θh_DiffCache[ vh_Idx ]
    #     θh = vh_θh_DiffCache[ θh_Idx ]
        
    # end

    #----------------------------------------
    #----------------------------------------
    
    vh_θh_DiffCache =
        get_tmp( vh_θh_DiffCache, intra_dyn_vh_θh )

    vh_θh_DiffCache .= intra_dyn_vh_θh
    
    vh = vh_θh_DiffCache[vh_Idx]

    θh = vh_θh_DiffCache[θh_Idx]
        
    #----------------------------------------    
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]
    
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

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
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

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...]                  
            
        end
    end
    
    #----------------------------------------
    
    vh_θh_id_iq =
        [ vh; θh; gens_i_d; gens_i_q ]

    #--------------------------------------------

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch    
    
    #--------------------------------------------

    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para),
                      vh_θh_id_iq,
                      method=:trust_region,
                      autodiff=:forward )
        return sol.zero
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para  ) ),
                vh_θh_id_iq,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        return sol.u
        
    end
        
end



#-------------------------------------------------------

function intra_dyn_pf_by_vh_θh_func!(
    intra_dyn_vh_θh, intra_dyn_pf_flat_para ;
    
    vh_θh_DiffCache =
        vh_θh_DiffCache,
    
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver = pf_solver,
    
    use_nlsolve =
        false,
    
    post_pf_idxs =
        post_pf_idxs

    )

    #----------------------------------------    

    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver

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
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (
     loc_load_exist,
     vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         intra_dyn_pf_mismatch_kwd_para

    #----------------------------------------
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        vars_Idxs
    
    #----------------------------------------    
    
  (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
  
    #----------------------------------------    

    (;     
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) = 
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs

    #----------------------------------------    

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #----------------------------------------    

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para

    #----------------------------------------     
    #----------------------------------------    
    
    P_gens =
        intra_dyn_pf_flat_para[ dyn_P_gens_Idxs ]

    Q_gens =
        intra_dyn_pf_flat_para[ dyn_Q_gens_Idxs ]

    P_non_gens  =
        intra_dyn_pf_flat_para[ dyn_P_non_gens_Idxs ]

    Q_non_gens = 
        intra_dyn_pf_flat_para[ dyn_Q_non_gens_Idxs ]

    flat_δ_ed_dash_eq_dash =
        intra_dyn_pf_flat_para[ dyn_δ_ed_eq_pf_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            intra_dyn_pf_flat_para[ dyn_P_g_loc_load_Idxs ]

        Q_g_loc_load =
            intra_dyn_pf_flat_para[ dyn_Q_g_loc_load_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #----------------------------------------      
    #----------------------------------------

    # if use_nlsolve == true

    #     vh = intra_dyn_vh_θh[ vh_Idx ]
    #     θh = intra_dyn_vh_θh[ θh_Idx ]
        
    # else
    
    #     vh_θh_DiffCache =
    #         get_tmp( vh_θh_DiffCache,
    #                 intra_dyn_vh_θh )
    #     vh_θh_DiffCache .= intra_dyn_vh_θh
    #     vh = vh_θh_DiffCache[ vh_Idx ]
    #     θh = vh_θh_DiffCache[ θh_Idx ]
        
    # end

    #----------------------------------------
    #----------------------------------------
    
    vh_θh_DiffCache =
        get_tmp( vh_θh_DiffCache, intra_dyn_vh_θh )

    vh_θh_DiffCache .= intra_dyn_vh_θh
    
    vh = vh_θh_DiffCache[vh_Idx]

    θh = vh_θh_DiffCache[θh_Idx]
        
    #----------------------------------------    
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]
    
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

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
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

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...]                  
            
        end
    end
    
    #----------------------------------------
    
    vh_θh_id_iq =
        [ vh; θh; gens_i_d; gens_i_q ]

    #--------------------------------------------

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch    
    
    #--------------------------------------------

    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para),
                      vh_θh_id_iq,
                      method=:trust_region,
                      autodiff=:forward )
        return sol.zero
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para  ) ),
                vh_θh_id_iq,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        return sol.u
        
    end
        
end



# ------------------------------------------------------


function intra_dyn_pf_by_vh_θh_func!(
    intra_vh_θh_wt_intra_dyn_pf_flat_para;
    
    vh_θh_DiffCache =
        vh_θh_DiffCache,
    
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver = pf_solver,
    
    use_nlsolve =
        false,
    
    post_pf_idxs =
        post_pf_idxs

    )

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
    
    intra_dyn_vh_θh =
        intra_vh_θh_wt_intra_dyn_pf_flat_para[
            intra_flat_vh_θh_Idx]
    
    intra_dyn_pf_flat_para =
        intra_vh_θh_wt_intra_dyn_pf_flat_para[
            intra_dyn_pf_flat_para_Idx ]
    
    #----------------------------------------

    return intra_dyn_pf_by_vh_θh_func!(
        intra_dyn_vh_θh, intra_dyn_pf_flat_para ;
        vh_θh_DiffCache =
            vh_θh_DiffCache,
        Pg_Qg_external_control =
            Pg_Qg_external_control,
        intra_pf_kwd_para =
            intra_pf_kwd_para,
        pf_solver =
            pf_solver,    
    use_nlsolve =
        use_nlsolve,    
    post_pf_idxs =
        post_pf_idxs )


    
end


# ------------------------------------------------------


function intra_dyn_pf_by_ur_ui_func!(
    intra_ur_ui, intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache =
        ur_ui_DiffCache,
    
    vh_θh_DiffCache =
            vh_θh_DiffCache,
    
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver = pf_solver,
    
    use_nlsolve =
        false,
    post_pf_idxs =
        post_pf_idxs

    )

    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 


    ur_ui_DiffCache =
        get_tmp(ur_ui_DiffCache,
                intra_ur_ui )

    ur_ui_DiffCache .= intra_ur_ui 
    
    vh =  [ abs( x_from_xr_xi(
        ur_ui_DiffCache[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]

    θh =  [ angle( x_from_xr_xi(
        ur_ui_DiffCache[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]

    
    intra_dyn_vh_θh = [vh; θh ]

    return  intra_dyn_pf_by_vh_θh_func!(
        intra_dyn_vh_θh, intra_dyn_pf_flat_para ;
        vh_θh_DiffCache =
            vh_θh_DiffCache,
        Pg_Qg_external_control =
            Pg_Qg_external_control,
        intra_pf_kwd_para =
            intra_pf_kwd_para,
        pf_solver =
            pf_solver,    
        use_nlsolve =
            use_nlsolve,    
        post_pf_idxs =
            post_pf_idxs )
    
end

# ------------------------------------------------------

function intra_dyn_pf_by_ur_ui_func!(
    intra_ur_ui_wt_intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache =
        ur_ui_DiffCache,
    
    vh_θh_DiffCache =
            vh_θh_DiffCache,
    
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver = pf_solver,
    
    use_nlsolve =
        false,
    post_pf_idxs =
        post_pf_idxs

    )

    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 

    
    intra_dyn_ur_ui =
        intra_ur_ui_wt_intra_dyn_pf_flat_para[
            intra_flat_ur_ui_Idx]
    
    intra_dyn_pf_flat_para =
        intra_ur_ui_wt_intra_dyn_pf_flat_para[
            intra_dyn_pf_flat_para_Idx ]

    return intra_dyn_pf_by_ur_ui_func!(
        intra_ur_ui, intra_dyn_pf_flat_para;
    
        ur_ui_DiffCache =
            ur_ui_DiffCache,

        vh_θh_DiffCache =
                vh_θh_DiffCache,

        Pg_Qg_external_control =
            Pg_Qg_external_control,

        intra_pf_kwd_para =
            intra_pf_kwd_para,

        pf_solver = pf_solver,

        use_nlsolve =
            use_nlsolve,
        
        post_pf_idxs =
            post_pf_idxs )
    
end

# ------------------------------------------------------
# ------------------------------------------------------


function dyn_pf_by_vh_θh_func!(
    intra_vh_θh, intra_dyn_pf_flat_para;
    
    vh_θh_DiffCache =
        vh_θh_DiffCache,
    
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver =
        pf_solver,
    
    use_nlsolve =
        false,

    post_pf_idxs =
        post_pf_idxs

    )
        
    return  intra_dyn_pf_by_vh_θh_func!(
        intra_dyn_vh_θh, intra_dyn_pf_flat_para ;
        
        vh_θh_DiffCache =
            vh_θh_DiffCache,
        
        Pg_Qg_external_control =
            Pg_Qg_external_control,
        
        intra_pf_kwd_para =
            intra_pf_kwd_para,
        
        pf_solver =
            pf_solver,
        
        use_nlsolve =
            use_nlsolve,
        
        post_pf_idxs =
            post_pf_idxs )
    
end


# ------------------------------------------------------

function dyn_pf_by_vh_θh_func!(
    intra_vh_θh_wt_intra_dyn_pf_flat_para;
    
    vh_θh_DiffCache =
        vh_θh_DiffCache,
    
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver =
        pf_solver,
    
    use_nlsolve =
        false,

    post_pf_idxs =
        post_pf_idxs

    )
    
    return intra_dyn_pf_by_vh_θh_func!(
        intra_vh_θh_wt_intra_dyn_pf_flat_para;
        
        vh_θh_DiffCache =
            vh_θh_DiffCache,
        
        Pg_Qg_external_control =
            Pg_Qg_external_control,
        
        intra_pf_kwd_para =
            intra_pf_kwd_para,
        
        pf_solver =
            pf_solver,
        
        use_nlsolve =
            use_nlsolve,
        
        post_pf_idxs =
            post_pf_idxs)
end

# ------------------------------------------------------

function dyn_pf_by_ur_ui_func!(
    intra_ur_ui, intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache =
        ur_ui_DiffCache,

    vh_θh_DiffCache =

        vh_θh_DiffCache,
        
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver =
        pf_solver,
    
    use_nlsolve =
        false,

    post_pf_idxs =
        post_pf_idxs

    )


    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 


    ur_ui_DiffCache =
        get_tmp(ur_ui_DiffCache,
                intra_ur_ui )

    ur_ui_DiffCache .= intra_ur_ui
    
    vh =  [ abs( x_from_xr_xi(
        ur_ui_DiffCache[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]

    θh =  [ angle( x_from_xr_xi(
        ur_ui_DiffCache[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]
    
    intra_dyn_vh_θh = [ vh; θh ]

    return intra_dyn_pf_by_vh_θh_func!(
        intra_dyn_vh_θh, intra_dyn_pf_flat_para ;
   
        vh_θh_DiffCache =

            vh_θh_DiffCache,

        Pg_Qg_external_control =
            Pg_Qg_external_control,

        intra_pf_kwd_para =
            intra_pf_kwd_para,

        pf_solver = pf_solver,

        use_nlsolve =
            use_nlsolve,

        post_pf_idxs =
            post_pf_idxs

        )

end

# ------------------------------------------------------

function dyn_pf_by_ur_ui_func!(
    intra_ur_ui_wt_intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache =
        ur_ui_DiffCache,

    vh_θh_DiffCache =

        vh_θh_DiffCache,
        
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver =
        pf_solver,
    
    use_nlsolve =
        false,

    post_pf_idxs =
        post_pf_idxs

    )


    (;          
     intra_flat_ur_ui_Idx,

     intra_flat_vh_θh_Idx,
     
     intra_dyn_pf_flat_para_Idx,

     intra_dyn_pf_mismatch_kwd_para,

     ur_ui_idx_in_Idx,
          
     nodes_u_Idx_in_ranges,
     
     nodes_δ_ed_dash_eq_dash_Idxs
     
     ) = intra_pf_kwd_para 

    
    intra_ur_ui =
        intra_ur_ui_wt_intra_dyn_pf_flat_para[
            intra_flat_ur_ui_Idx]
    
    intra_dyn_pf_flat_para =
        intra_ur_ui_wt_intra_dyn_pf_flat_para[
            intra_dyn_pf_flat_para_Idx ]
    
    return dyn_pf_by_ur_ui_func!(
        intra_ur_ui, intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache =
        ur_ui_DiffCache,

    vh_θh_DiffCache =

        vh_θh_DiffCache,
        
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver =
        pf_solver,
    
    use_nlsolve =
        use_nlsolve,

    post_pf_idxs =
        post_pf_idxs

    )

end

# ------------------------------------------------------

function intra_dyn_pf_func!(
    intra_dyn_ur_ui, intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache =
        ur_ui_DiffCache,
    
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver = pf_solver,
    
    use_nlsolve = false

    )

    #----------------------------------------    

    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver

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
    
    (
     loc_load_exist,
     vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         intra_dyn_pf_mismatch_kwd_para

    #----------------------------------------
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        vars_Idxs
    
    #----------------------------------------    
    
  (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
  
    #----------------------------------------    

    (;     
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) = 
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs

    #----------------------------------------    

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 
    
    #----------------------------------------    

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para

    #----------------------------------------     
    #----------------------------------------    
    
    P_gens =
        intra_dyn_pf_flat_para[ dyn_P_gens_Idxs ]

    Q_gens =
        intra_dyn_pf_flat_para[ dyn_Q_gens_Idxs ]

    P_non_gens  =
        intra_dyn_pf_flat_para[ dyn_P_non_gens_Idxs ]

    Q_non_gens = 
        intra_dyn_pf_flat_para[ dyn_Q_non_gens_Idxs ]

    flat_δ_ed_dash_eq_dash =
        intra_dyn_pf_flat_para[ dyn_δ_ed_eq_pf_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            intra_dyn_pf_flat_para[ dyn_P_g_loc_load_Idxs ]

        Q_g_loc_load =
            intra_dyn_pf_flat_para[ dyn_Q_g_loc_load_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #----------------------------------------      
    #----------------------------------------

    # if use_nlsolve == true

    #     vh =  [ abs( x_from_xr_xi(
    #         intra_dyn_ur_ui[ idx ] ) )
    #             for idx in
    #                 ur_ui_idx_in_Idx ]

    #     θh =  [ angle( x_from_xr_xi(
    #         intra_dyn_ur_ui[ idx ] ) )
    #             for idx in
    #                 ur_ui_idx_in_Idx ]
        
    # else
    
    #     ur_ui_DiffCache =
    #         get_tmp(ur_ui_DiffCache, intra_dyn_ur_ui )

    #     ur_ui_DiffCache .= intra_dyn_ur_ui

    #     vh =  [ abs( x_from_xr_xi(
    #         ur_ui_DiffCache[ idx ] ) )
    #             for idx in
    #                 ur_ui_idx_in_Idx ]

    #     θh =  [ angle( x_from_xr_xi(
    #         ur_ui_DiffCache[ idx ] ) )
    #             for idx in
    #                 ur_ui_idx_in_Idx ]
        
    # end


    # vh =  [ abs( x_from_xr_xi(
    #     intra_dyn_ur_ui[ idx ] ) )
    #         for idx in
    #             ur_ui_idx_in_Idx ]

    # θh =  [ angle( x_from_xr_xi(
    #     intra_dyn_ur_ui[ idx ] ) )
    #         for idx in
    #             ur_ui_idx_in_Idx ]
    

    ur_ui_DiffCache =
        get_tmp(ur_ui_DiffCache, intra_dyn_ur_ui )

    ur_ui_DiffCache .= intra_dyn_ur_ui

    vh =  [ abs( x_from_xr_xi(
        ur_ui_DiffCache[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]

    θh =  [ angle( x_from_xr_xi(
        ur_ui_DiffCache[ idx ] ) )
            for idx in
                ur_ui_idx_in_Idx ]
    
    #----------------------------------------
    
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]
    
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

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
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

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...]                  
            
        end
    end
    
    #----------------------------------------
    
    vh_θh_id_iq =
        [ vh; θh; gens_i_d; gens_i_q ]

    #--------------------------------------------

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch    
    
    #--------------------------------------------

    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para),
                      vh_θh_id_iq,
                      method=:trust_region,
                      autodiff=:forward )
        return sol.zero
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para  ) ),
                vh_θh_id_iq,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        return sol.u
        
    end
        
end


# ------------------------------------------------------


function intra_dyn_pf_func!(
    intra_ur_ui_wt_intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache =
        ur_ui_DiffCache,
    
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver = pf_solver,
    
    use_nlsolve = false

    )
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

    intra_dyn_ur_ui =
        intra_ur_ui_wt_intra_dyn_pf_flat_para[
            intra_flat_ur_ui_Idx]
    
    intra_dyn_pf_flat_para =
        intra_ur_ui_wt_intra_dyn_pf_flat_para[
            intra_dyn_pf_flat_para_Idx ]
    
    #----------------------------------------
   

    return intra_dyn_pf_func!(
        intra_dyn_ur_ui, intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache =
        ur_ui_DiffCache,
    
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver = pf_solver,
    
    use_nlsolve = use_nlsolve

    )
        
end


# ------------------------------------------------------

function dyn_pf_func!(
     intra_dyn_ur_ui, intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache  =
        ur_ui_DiffCache,
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver = pf_solver,
    
    use_nlsolve = false

    )

    return intra_dyn_pf_func!(
            intra_dyn_ur_ui, intra_dyn_pf_flat_para;
            ur_ui_DiffCache =
                ur_ui_DiffCache,
            Pg_Qg_external_control =
                Pg_Qg_external_control,
            intra_pf_kwd_para =
                intra_pf_kwd_para,
            pf_solver =
                pf_solver,
            use_nlsolve =
                use_nlsolve)
end



function dyn_pf_func!(
     intra_ur_ui_wt_intra_dyn_pf_flat_para;
    
    ur_ui_DiffCache  =
        ur_ui_DiffCache,
    Pg_Qg_external_control =
        Pg_Qg_external_control,
    
    intra_pf_kwd_para =
        intra_pf_kwd_para,
    
    pf_solver = pf_solver,
    
    use_nlsolve = false

    )

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

    intra_dyn_ur_ui =
        intra_ur_ui_wt_intra_dyn_pf_flat_para[
            intra_flat_ur_ui_Idx]
    
    intra_dyn_pf_flat_para =
        intra_ur_ui_wt_intra_dyn_pf_flat_para[
            intra_dyn_pf_flat_para_Idx ]
    
    #----------------------------------------
   
    
    return intra_dyn_pf_func!(
            intra_dyn_ur_ui, intra_dyn_pf_flat_para;
            ur_ui_DiffCache =
                ur_ui_DiffCache,
            Pg_Qg_external_control =
                Pg_Qg_external_control,
            intra_pf_kwd_para =
                intra_pf_kwd_para,
            pf_solver =
                pf_solver,
            use_nlsolve =
                use_nlsolve)
end


# ------------------------------------------------------
# ------------------------------------------------------

function dyn_pf_dfdp_func!(
    intra_dyn_pf_flat_para;
    ur_ui_DiffCache =
        dyn_pf_sol,
    dyn_pf_sense_kwd =
        dyn_pf_sense_kwd )

    (;
     ur_ui_Cache,
     Pg_Qg_external_control,
     intra_pf_kwd_para ,
     pf_solver,
     use_nlsolve ) =
         dyn_pf_sense_kwd
    
    ForwardDiff.jacobian(
        ( p ) -> dyn_pf_func!( ur_ui_DiffCache, p;
                               ur_ui_Cache  =
                                   ur_ui_Cache,
                               Pg_Qg_external_control =
                                   Pg_Qg_external_control,   
                               intra_pf_kwd_para =
                                   intra_pf_kwd_para,   
                               pf_solver =
                                   pf_solver,   
                               use_nlsolve =
                                   use_nlsolve ),
        intra_dyn_pf_flat_para )

    
end


function dyn_pf_dfdx_func!(
    dyn_pf_sol;
    intra_dyn_pf_flat_para =
        intra_dyn_pf_flat_para,
    dyn_pf_sense_kwd =
        dyn_pf_sense_kwd )

    (;
     ur_ui_Cache,
     Pg_Qg_external_control,
     intra_pf_kwd_para ,
     pf_solver,
     use_nlsolve ) =
         dyn_pf_sense_kwd

    
    
    dyn_pf_∂f∂x = ForwardDiff.jacobian(
        ( x ) -> dyn_pf_func!( dyn_pf_sol,
                               intra_dyn_pf_flat_para;
                               ur_ui_Cache =
                                   ur_ui_Cache,
                               Pg_Qg_external_control =
                                   Pg_Qg_external_control,   
                               intra_pf_kwd_para =
                                   intra_pf_kwd_para,   
                               pf_solver =
                                   pf_solver,   
                               use_nlsolve =
                                   use_nlsolve ),
        dyn_pf_sol )


end

#----------------------------------------------------

function intra_dyn_pf_ΔPQ_Δidq_mismatch_dfdx(
    post_pf_sol;
    intra_dyn_pf_mismatch_flat_para =
        intra_dyn_pf_flat_para,
    intra_dyn_pf_mismatch_kwd_para =
        intra_dyn_pf_mismatch_kwd_para )

    ΔPQ_Δidq = similar( post_pf_sol )

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch
    
    ∂f∂x = ForwardDiff.jacobian(
        ( ΔPQ_Δidq, x ) ->
            pf_fun_mismatch(
                ΔPQ_Δidq , x,
                intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para  ),
         ΔPQ_Δidq, post_pf_sol )


end



function intra_dyn_pf_by_vh_θh_by_parts_dfdp_func!(
    intra_dyn_pf_flat_para ;
    post_pf_sol =
        post_pf_sol,
    intra_dyn_pf_mismatch_kwd_para =
        intra_dyn_pf_mismatch_kwd_para )

    ΔPQ_Δidq = similar( post_pf_sol )

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch
    
    ∂f∂p = ForwardDiff.jacobian(
        ( ΔPQ_Δidq, p ) ->
            pf_fun_mismatch(
                ΔPQ_Δidq , post_pf_sol, p;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para  ),
         ΔPQ_Δidq, intra_dyn_pf_flat_para )


end

#########################################################


function intra_dyn_pf_by_ur_ui_by_parts_func!(
    ur_ui_DiffCache,
    flat_δ_ed_dash_eq_dash,
    dyn_pf_flat_para ;
    intra_dyn_pf_kwd_para =
        intra_dyn_pf_kwd_para )

    #----------------------------------------

    (;
     ur_ui_idx_order,
     flat_ur_flat_ui_Idx,
     Pg_Qg_external_control,
     intra_pf_kwd_para,
     pf_solver,
     use_nlsolve,
     post_pf_idxs,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs )  =
         intra_dyn_pf_kwd_para

    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver

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
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (
     loc_load_exist,
     vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         intra_dyn_pf_mismatch_kwd_para

    #----------------------------------------
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        vars_Idxs
    
    #----------------------------------------    
    
  (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
  
    #----------------------------------------    

    (;     
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) = 
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs

    #----------------------------------------    

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 

    #----------------------------------------    
    
    (; Pg_Idx,
      Qg_Idxs,
      Png_Idxs,
      Qng_Idxs,
      Pgll_Idxs,
      Qgll_Idxs ) = Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------    

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para

    #----------------------------------------     
    #----------------------------------------    
    
    P_gens =
        dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        intra_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        intra_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        intra_dyn_pf_flat_para[ Qng_Idxs ]

    # flat_δ_ed_dash_eq_dash =
    #     intra_dyn_pf_flat_para[ dyn_δ_ed_eq_pf_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            intra_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            intra_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #----------------------------------------

    
    if ur_ui_idx_order == :consecutive
        vec_nodes_ur_ui =
        [ ur_ui_DiffCache[idx]
          for idx in  ur_ui_idx_in_Idx ]
    else
        
        flat_ur_Idx, flat_ui_Idx = flat_ur_flat_ui_Idx

        flat_ur = ur_ui_DiffCache[ flat_ur_Idx ]
        
        flat_ui = ur_ui_DiffCache[ flat_ui_Idx ]

        vec_nodes_ur_ui =
            [[ur, ui] for (ur, ui) in
                 zip( flat_ur, flat_ui )]        
    end

    #----------------------------------------
    
    vec_vh_θh =
        cartesian_to_polar.( vec_nodes_ur_ui  )

    vh = first.( vec_vh_θh )

    θh = last.( vec_vh_θh )
    
    #----------------------------------------
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]
    
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

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
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

    if loc_load_exist == true

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...]                  
            
        end
    end
    
    #----------------------------------------
    
    vh_θh_id_iq =
        [ vh; θh; gens_i_d; gens_i_q ]

    #--------------------------------------------

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch    
    
    #--------------------------------------------

    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para),
                      vh_θh_id_iq,
                      method=:trust_region,
                      autodiff=:forward )
        return sol.zero
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para
                    ) ),
                
                vh_θh_id_iq,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        return sol.u
        
    end
        
end


#-----------------------------------------------


function intra_dyn_pf_by_ur_ui_by_stateDiff_func!(
    stateDiffCache, dyn_pf_flat_para ;
    intra_dyn_pf_kwd_para =
        intra_dyn_pf_kwd_para )

    #----------------------------------------

    (;
     
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
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs )  =
         intra_dyn_pf_kwd_para

    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver

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
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (
     loc_load_exist,
     vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         intra_dyn_pf_mismatch_kwd_para

    #----------------------------------------
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        vars_Idxs
    
    #----------------------------------------    
    
  (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
  
    #----------------------------------------    

    (;     
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) = 
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs

    #----------------------------------------    

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 

    #----------------------------------------    
    
    (; Pg_Idx,
      Qg_Idxs,
      Png_Idxs,
      Qng_Idxs,
      Pgll_Idxs,
      Qgll_Idxs ) = Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------    

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para

    #----------------------------------------     
    #----------------------------------------    
    
    P_gens =
        dyn_pf_flat_para[ Pg_Idx ]

    Q_gens =
        intra_dyn_pf_flat_para[ Qg_Idxs ]

    P_non_gens  =
        intra_dyn_pf_flat_para[ Png_Idxs ]

    Q_non_gens = 
        intra_dyn_pf_flat_para[ Qng_Idxs ]

    # flat_δ_ed_dash_eq_dash =
    #     intra_dyn_pf_flat_para[ dyn_δ_ed_eq_pf_Idxs ]
    
    
    if loc_load_exist == true

        P_g_loc_load =
            intra_dyn_pf_flat_para[ Pgll_Idxs ]

        Q_g_loc_load =
            intra_dyn_pf_flat_para[ Qgll_Idxs ]
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end
    
    #----------------------------------------
    
    # δ_ed_dash_eq_dash_Idxs_in_flattend
    
    # gens_nodes_δ_ed_dash_eq_dash =
    #     [ stateDiffCache[ idx ]
    #       for idx in
    #           nodes_δ_ed_dash_eq_dash_Idxs ]
    
    # flat_Idx_nodes_δ_ed_dash_eq_dash_Idxs =
    #     [ nodes_δ_ed_dash_eq_dash_Idxs...;]

    flat_δ_ed_dash_eq_dash =
        stateDiffCache[
            flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    gens_δ_ed_dash_eq_dash =
        [ flat_δ_ed_dash_eq_dash[idx]
          for idx in
              δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------

    """
    
    nodes_ur_Idx = first.(nodes_u_Idx_in_ranges)
        
    nodes_ui_Idx = last.(nodes_u_Idx_in_ranges)

    flat_nodes_ur_flat_nodes_ui_Idx  =
        [ nodes_ur_Idx; nodes_ui_Idx ]

    flat_nodes_ur_flat_nodes_ui =
        stateDiffCache[ flat_nodes_ur_flat_nodes_ui_Idx ]

    flat_ur_Idx, flat_ui_Idx =  flat_ur_flat_ui_Idx

    flat_nodes_ur =
        flat_nodes_ur_flat_nodes_ui[ flat_ur_Idx ]

    flat_nodes_ui =
        flat_nodes_ur_flat_nodes_ui[ flat_ui_Idx ]

    """
    
    if ur_ui_idx_order == :consecutive

        # vec_nodes_ur_ui =
        # [ stateDiffCache[idx]
        #   for idx in  nodes_u_Idx_in_ranges ]

        # nodes_ur_Idx_in_state =
        #     first.(nodes_u_Idx_in_ranges)
        
        # nodes_ui_Idx_in_state =
        #     last.(nodes_u_Idx_in_ranges)
        
        flat_ur = stateDiffCache[ nodes_ur_Idx_in_state  ]
        
        flat_ui = stateDiffCache[ nodes_ui_Idx_in_state ]

        vec_nodes_ur_ui =
            [[ur, ui] for (ur, ui) in
                 zip( flat_ur, flat_ui )]        
        
    else

        flat_ur_Idx, flat_ui_Idx =  flat_ur_flat_ui_Idx
        
        nodes_ur_Idx =
            non_consecutive_ur_ui_idxs[ flat_ur_Idx ]
        
        nodes_ui_Idx =
            non_consecutive_ur_ui_idxs[ flat_ui_Idx ]
        
        flat_ur = stateDiffCache[ nodes_ur_Idx  ]
        
        flat_ui = stateDiffCache[ nodes_ui_Idx ]

        vec_nodes_ur_ui =
            [[ur, ui] for (ur, ui) in
                 zip( flat_ur, flat_ui )]        
    end

    #----------------------------------------
    
    vec_vh_θh =
        cartesian_to_polar.( vec_nodes_ur_ui  )

    vh = first.( vec_vh_θh )

    θh = last.( vec_vh_θh )
    
    #----------------------------------------
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]
    
    #----------------------------------------

    # gens_δ_ed_dash_eq_dash = [
    #     flat_δ_ed_dash_eq_dash[idx]
    #     for idx in
    #         δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
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

    if loc_load_exist == true

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...] 
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...] 
            
        end
    end
    
    #----------------------------------------
    
    vh_θh_id_iq =
        [ vh; θh; gens_i_d; gens_i_q ]

    #--------------------------------------------

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch    
    
    #--------------------------------------------

    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para),
                      vh_θh_id_iq,
                      method=:trust_region,
                      autodiff=:forward )
        return sol.zero
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para
                    ) ),
                
                vh_θh_id_iq,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        return sol.u
        
    end
        
end


#-----------------------------------------------


function intra_dyn_pf_by_vh_θh_by_stateDiff_func!(
    stateDiffCache,
    init_dyn_pf_flat_para ;
    intra_dyn_pf_kwd_para =
        intra_dyn_pf_kwd_para )

    #----------------------------------------

    (;     
     flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state,

     nodes_vh_Idx_in_state,
     nodes_θh_Idx_in_state,
     nodes_vh_θh_Idx_in_state,

     vh_θh_idx_order,
     flat_vh_flat_θh_Idx,
     
     Pg_Qg_external_control,
     
     intra_pf_kwd_para,
     pf_solver,
     use_nlsolve,
     post_pf_idxs,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs )  =
         intra_dyn_pf_kwd_para

    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver

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
     vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         post_pf_idxs
    
    #----------------------------------------    
    
    (
     loc_load_exist,
     vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         intra_dyn_pf_mismatch_kwd_para

    #----------------------------------------
    
   (gens_vh_Idxs,
    gens_θh_Idxs,
    
    non_gens_nodes_vh_Idxs,
    non_gens_nodes_θh_Idxs,
    
    gen_id_Idxs,
    gen_iq_Idxs) =
        vars_Idxs
    
    #----------------------------------------    
    
  (;
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_net_para,
     δ_ed_dash_eq_dash_Idxs_in_flattend
      ) =
         dyn_pf_Idxs_kwd_para
  
    #----------------------------------------    

    (;     
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) = 
        dyn_pf_P_Q_δ_etc_kwd_para_Idxs
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        dyn_pf_fun_kwd_n2s_idxs

    #----------------------------------------    

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx
    ) = dyn_pf_fun_kwd_net_idxs 

    #----------------------------------------    
    
    (; Pg_Idx,
      Qg_Idxs,
      Png_Idxs,
      Qng_Idxs,
      Pgll_Idxs,
      Qgll_Idxs ) = Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------    

   (; 
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx) =
        dyn_pf_fun_kwd_net_para

    #----------------------------------------     
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
    
    #----------------------------------------

    flat_δ_ed_dash_eq_dash =
        stateDiffCache[
            flat_nodes_δ_ed_dash_eq_dash_Idxs_in_state ]
    
    gens_δ_ed_dash_eq_dash =
        [ flat_δ_ed_dash_eq_dash[idx]
          for idx in
              δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
    #----------------------------------------
    
    if vh_θh_idx_order == :consecutive
        
        flat_vh = stateDiffCache[ nodes_vh_Idx_in_state  ]
        
        flat_θh = stateDiffCache[ nodes_θh_Idx_in_state ]

        vec_nodes_vh_θh =
            [[a_vh, a_θh] for (a_vh, a_θh) in
                 zip( flat_vh, flat_θh )]        
        
    else

        flat_vh_Idx, flat_θh_Idx  = flat_vh_flat_θh_Idx
                
        flat_vh =
            stateDiffCache[
                nodes_vh_θh_Idx_in_state ][
                    flat_vh_Idx ]
        
        flat_θh = stateDiffCache[
            nodes_vh_θh_Idx_in_state ][
                flat_θh_Idx ]

        vec_nodes_vh_θh =
            [[vh, θh] for (vh, θh) in
                 zip( flat_vh, flat_θh )]        
    end

    #----------------------------------------

    vh = first.( vec_nodes_vh_θh )

    θh = last.( vec_nodes_vh_θh )
    
    #----------------------------------------
    
    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]
    
    #----------------------------------------

    # gens_δ_ed_dash_eq_dash = [
    #     flat_δ_ed_dash_eq_dash[idx]
    #     for idx in
    #         δ_ed_dash_eq_dash_Idxs_in_flattend ]
    
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

    if loc_load_exist == true

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...;
                P_g_loc_load...;
                Q_g_loc_load...]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens...;
                Q_gens...;            
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...] 
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph...;
                gens_dyn_qh...;
                P_non_gens...;
                Q_non_gens...;
                gens_δ_ed_dash_eq_dash...] 
            
        end
    end
    
    #----------------------------------------
    
    vh_θh_id_iq =
        [ vh; θh; gens_i_d; gens_i_q ]

    #--------------------------------------------

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch    
    
    #--------------------------------------------

    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para),
                      vh_θh_id_iq,
                      method=:trust_region,
                      autodiff=:forward )
        return sol
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para
                    ) ),
                
                vh_θh_id_iq,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        return sol
        
    end
        
end

#----------------------------------------    
#----------------------------------------    

function intra_dyn_current_balance_pf_by_vh_θh_by_parts_func!(
    init_flat_vh_flat_θh,
    flat_δ_ed_dash_eq_dash,
    init_dyn_pf_flat_para;
    intra_dyn_pf_kwd_para =
        intra_dyn_pf_kwd_para    
    )

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
     pf_solver) =
         intra_dyn_pf_kwd_para
    
    #----------------------------------------
    #----------------------------------------

    (;
     pf_alg,
     abstol,
     reltol ) =
         pf_solver
    
    #----------------------------------------
    
    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
    
    #----------------------------------------     
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

    #----------------------------------------
    
    flat_vh_idx_in_flat_vh_flat_θh, flat_θh_idx_in_flat_vh_flat_θh = flat_vh_flat_θh_Idx
        
    flat_vh = init_flat_vh_flat_θh[
        flat_vh_idx_in_flat_vh_flat_θh ]
        
    flat_θh = init_flat_vh_flat_θh[
        flat_θh_idx_in_flat_vh_flat_θh ]

    #----------------------------------------
    
    gens_vh = flat_vh[ gens_nodes_idx ]

    gens_θh = flat_θh[ gens_nodes_idx ]
        
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

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
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

    if loc_load_exist == true

        if Pg_Qg_external_control == true
        
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
                gens_dyn_ph;
                gens_dyn_qh;
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash;
                P_g_loc_load;
                Q_g_loc_load ]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]   
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph;
                gens_dyn_qh;
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]  
            
        end
    end
    
    #----------------------------------------
    
    vh_θh_id_iq =
        [ init_flat_vh_flat_θh; gens_i_d; gens_i_q ]
    
    #----------------------------------------    

    
    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_current_balance_pf_ΔI_idq_mismatch    
    
    #--------------------------------------------

    if use_nlsolve == true

        return nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para),
                      vh_θh_id_iq,
                      method=:trust_region,
                      autodiff=:forward )
        
    else

        return NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para  ) ),
                vh_θh_id_iq,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        
        
    end
        
end




function init_dyn_pf_by_vh_θh_by_parts_func!(
    init_flat_vh_flat_θh,
    init_dyn_pf_flat_para;
    init_pf_by_vh_θh_by_parts_kwd =
        init_pf_by_vh_θh_by_parts_kwd )

    #----------------------------------------

    (;     
     pf_solver,
     use_nlsolve,
     init_dyn_pf_mismatch_kwd_para )  =
         init_pf_by_vh_θh_by_parts_kwd

    (;
     loc_load_exist,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,     
     dyn_pf_vars_Idxs,
     dyn_pf_Idxs_kwd_para,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         init_dyn_pf_mismatch_kwd_para
    
    (;
     pf_alg,
     abstol,
     reltol) =
         pf_solver
    
    #----------------------------------------
    #----------------------------------------
    
    pf_fun_mismatch =
        get_a_model_integrated_init_dyn_pf_ΔPQ_mismatch
    

    #----------------------------------------
    
    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, init_dyn_pf_flat_para;
                init_dyn_pf_mismatch_kwd_para =
                    init_dyn_pf_mismatch_kwd_para),
                      init_flat_vh_flat_θh,
                      method=:trust_region,
                      autodiff=:forward )
        return sol
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        init_dyn_pf_mismatch_kwd_para  =
                            init_dyn_pf_mismatch_kwd_para
                    ) ),
                
                init_flat_vh_flat_θh,
                init_dyn_pf_flat_para ),
            pf_alg )

        return sol
        
    end
        
end


function intra_dyn_pf_by_vh_θh_by_parts_func!(
    init_flat_vh_flat_θh,
    flat_δ_ed_dash_eq_dash,
    init_dyn_pf_flat_para;
    intra_dyn_pf_kwd_para =
        intra_dyn_pf_kwd_para )

    #----------------------------------------

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
     pf_solver) =
         intra_dyn_pf_kwd_para
    
    #----------------------------------------

    (;
     pf_alg,
     abstol,
     reltol ) =
         pf_solver
    
    #----------------------------------------
    
    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
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
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #----------------------------------------
    
    # if vh_θh_idx_order == :consecutive
    #     vec_nodes_vh_θh =
    #     [ vh_θh_post_init_pf[idx]
    #       for idx in  vh_θh_idx_in_Idx ]
    # else
        
    #     flat_vh_Idx, flat_θh_Idx = flat_vh_flat_θh_Idx

    #     flat_vh = vh_θh_post_init_pf[ flat_vh_Idx ]
        
    #     flat_θh = vh_θh_post_init_pf[ flat_θh_Idx ]

    #     vec_nodes_vh_θh =
    #         [[vh, θh] for (vh, θh) in
    #              zip( flat_vh, flat_θh )]        
    # end

    
    flat_vh_idx_in_flat_vh_flat_θh, flat_θh_idx_in_flat_vh_flat_θh = flat_vh_flat_θh_Idx
        
    flat_vh = init_flat_vh_flat_θh[
        flat_vh_idx_in_flat_vh_flat_θh ]
        
    flat_θh = init_flat_vh_flat_θh[
        flat_θh_idx_in_flat_vh_flat_θh ]

    #----------------------------------------
    
    gens_vh = flat_vh[ gens_nodes_idx ]

    gens_θh = flat_θh[ gens_nodes_idx ]
        
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

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
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

    if loc_load_exist == true

        if Pg_Qg_external_control == true
        
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
                gens_dyn_ph;
                gens_dyn_qh;
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash;
                P_g_loc_load;
                Q_g_loc_load ]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]   
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph;
                gens_dyn_qh;
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]  
            
        end
    end
    
    #----------------------------------------
    
    vh_θh_id_iq =
        [ init_flat_vh_flat_θh; gens_i_d; gens_i_q ]

    #--------------------------------------------

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch    
    
    #--------------------------------------------

    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
                      vh_θh_id_iq,
                      method=:trust_region,
                      autodiff=:forward )
        return sol
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para
                    ) ),
                
                vh_θh_id_iq,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        return sol
        
    end
        
end



function intra_dyn_pf_by_vh_θh_δ_ed_eq_by_parts_func!(
    init_flat_vh_flat_θh,
    flat_δ_ed_dash_eq_dash,
    init_dyn_pf_flat_para;
    intra_dyn_pf_kwd_para =
        intra_dyn_pf_kwd_para )

    #----------------------------------------

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
     pf_solver) =
         intra_dyn_pf_kwd_para
    
    #----------------------------------------

    (;
     pf_alg,
     abstol,
     reltol ) =
         pf_solver
    
    #----------------------------------------
    
    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
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
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #----------------------------------------
    
    flat_vh_idx_in_flat_vh_flat_θh, flat_θh_idx_in_flat_vh_flat_θh = flat_vh_flat_θh_Idx
        
    flat_vh = init_flat_vh_flat_θh[
        flat_vh_idx_in_flat_vh_flat_θh ]
        
    flat_θh = init_flat_vh_flat_θh[
        flat_θh_idx_in_flat_vh_flat_θh ]

    #----------------------------------------
    
    gens_vh = flat_vh[ gens_nodes_idx ]

    gens_θh = flat_θh[ gens_nodes_idx ]
        
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

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
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

    if loc_load_exist == true

        if Pg_Qg_external_control == true
        
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
                gens_dyn_ph;
                gens_dyn_qh;
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash;
                P_g_loc_load;
                Q_g_loc_load ]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]   
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph;
                gens_dyn_qh;
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]  
            
        end
    end

    #--------------------------------------------

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_mismatch   
    
    #--------------------------------------------

    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
                      init_flat_vh_flat_θh,
                      method=:trust_region,
                      autodiff=:forward )
        return sol
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para
                    ) ),
                
                init_flat_vh_flat_θh,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        return sol
        
    end
        
end


function intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_func!(
    init_flat_vh_flat_θh_id_iq,
    flat_δ_ed_dash_eq_dash,
    init_dyn_pf_flat_para;
    kwd_para =
        intra_dyn_pf_by_vh_θh_id_iq_δ_ed_eq_by_parts_kwd )

    #----------------------------------------

    (;
     flat_vh_flat_θh_flat_id_iq_Idx,

     intra_dyn_pf_kwd_para) =
         kwd_para
     
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
     ) =
         intra_dyn_pf_kwd_para
    
    #----------------------------------------

    (;
     pf_alg,
     abstol,
     reltol ) =
         pf_solver
    
    #----------------------------------------
    
    (;
     Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
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
        
    else

        P_g_loc_load = [0.0]
            
        Q_g_loc_load = [0.0]
    end

    #----------------------------------------
    
    # flat_vh_idx_in_flat_vh_flat_θh, flat_θh_idx_in_flat_vh_flat_θh = flat_vh_flat_θh_Idx

    flat_vh_Idx, flat_θh_Idx, flat_id_Idx, flat_iq_Idx =
        flat_vh_flat_θh_flat_id_iq_Idx
        
    flat_vh = init_flat_vh_flat_θh_id_iq[
        flat_vh_Idx ]
        
    flat_θh = init_flat_vh_flat_θh_id_iq[
        flat_θh_Idx ]

    flat_id = init_flat_vh_flat_θh_id_iq[
        flat_id_Idx ]

    flat_iq = init_flat_vh_flat_θh_id_iq[
        flat_iq_Idx ]

    #----------------------------------------
    
    gens_vh = flat_vh[ gens_nodes_idx ]

    gens_θh = flat_θh[ gens_nodes_idx ]
        
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

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )
    
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

    if loc_load_exist == true

        if Pg_Qg_external_control == true
        
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
                gens_dyn_ph;
                gens_dyn_qh;
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash;
                P_g_loc_load;
                Q_g_loc_load ]                  
            
        end
    else

        if Pg_Qg_external_control == true
        
            intra_dyn_pf_mismatch_flat_para = [
                P_gens;
                Q_gens;            
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]   
            
        else
        
            intra_dyn_pf_mismatch_flat_para = [
                gens_dyn_ph;
                gens_dyn_qh;
                P_non_gens;
                Q_non_gens;
                flat_δ_ed_dash_eq_dash ]  
            
        end
    end

    #--------------------------------------------

    pf_fun_mismatch =
        get_a_model_integrated_intra_dyn_pf_ΔPQ_Δidq_mismatch
    
    #--------------------------------------------

    if use_nlsolve == true

        sol = nlsolve((g, x) ->
            pf_fun_mismatch(
                g, x, intra_dyn_pf_mismatch_flat_para;
                intra_dyn_pf_mismatch_kwd_para =
                    intra_dyn_pf_mismatch_kwd_para ),
                      init_flat_vh_flat_θh_id_iq,
                      method=:trust_region,
                      autodiff=:forward )
        return sol
        
    else

        sol = NonlinearSolve.solve(
            NonlinearProblem(
                NonlinearFunction( ( g, x, p ) ->
                    pf_fun_mismatch(
                        g, x, p;
                        intra_dyn_pf_mismatch_kwd_para  =
                            intra_dyn_pf_mismatch_kwd_para
                    ) ),
                
                init_flat_vh_flat_θh_id_iq,
                intra_dyn_pf_mismatch_flat_para ),
            pf_alg )

        return sol
        
    end
        
end


#-------------------------------------------------------
#-------------------------------------------------------
# Dynamic powerflow function for nodes_edges_dynamics!
# ------------------------------------------------------
#-------------------------------------------------------

function dyn_powerflow(
    netd,
    stateDiffCache,
    global_pf_param;
    maxiter=1e5,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    with_δ_ed_eq = true )

    (pf_net_param,
     sd_pf_views,
     mismatch) =
         global_pf_param

    _, _, _, _, pf_Idx, pf_and_dyn_idx_and_Idx, pf_net_misc = pf_net_param

    (slack_bus_idx,
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
    
    (working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view) =
         sd_pf_views
    
    δ_ω_ed_dash_eq_dash_view .=
        get_nodes_δ_ω_ed_dash_eq_dash_view(
            stateDiffCache,
            netd )
    
    uh_state = stateDiffCache[ur_ui_idx][ ur_IDX ] +
        im * stateDiffCache[ur_ui_idx][ui_IDX] 
    
    x0_vh_θh = [abs.(uh_state)...; angle.(uh_state)...]
    
    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]
    red_vh_θh_0      = [ red_vh_θh_0_view; ]

    dynamic_power_balance_powerflow(
        x0_vh_θh,
        mismatch,
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          δ_ω_ed_dash_eq_dash_view ),
        pf_net_param ;
        maxiter=maxiter,
        ftol=ftol ,
        xtol=xtol,
        with_δ_ed_eq = with_δ_ed_eq )

    return nothing
    
end


# -------------------------------------------------------
# -------------------------------------------------------
# Begin im model pf
# -------------------------------------------------------
# -------------------------------------------------------


function im_dyn_powerflow(
    netd,
    stateDiffCache,
    global_pf_param,
    im_model_dyn_pf_up_para,
    im_model_idq_pf_cal;
    maxiter=1e5, ftol=1000*eps() ,
    xtol=1000*eps(),
    with_δ_ed_eq = true )

    (pf_net_param,
     sd_pf_views,
     mismatch) = global_pf_param

    _, _, _, _, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    (slack_bus_idx,
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
    
    (working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view) =
         sd_pf_views

    (nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_δ_ed_dash_eq_dash_Idxs_in_state) =
         im_model_dyn_pf_up_para

   # idq_wt_pad_view, gen_idx = im_model_idq_pf_cal
    
    # δ_ω_ed_dash_eq_dash_view .= get_nodes_δ_ω_ed_dash_eq_dash_view( stateDiffCache, netd )

    update_im_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        δ_ω_ed_dash_eq_dash_view,
        stateDiffCache,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )
    
    uh_state = stateDiffCache[ur_ui_idx][ ur_IDX ] +
        im * stateDiffCache[ur_ui_idx][ui_IDX] 
    
    x0_vh_θh = [abs.(uh_state)...; angle.(uh_state)...]
    
    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]
    red_vh_θh_0      = [ red_vh_θh_0_view; ]

    im_model_dynamic_power_balance_powerflow(
        x0_vh_θh, mismatch,
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          δ_ω_ed_dash_eq_dash_view ),
        pf_net_param,
        im_model_idq_pf_cal ;
        maxiter=maxiter,
        ftol=ftol ,
        xtol=xtol,
        with_δ_ed_eq = with_δ_ed_eq )

    return nothing
    
end


# -------------------------------------------------------
# -------------------------------------------------------
# End im model pf
# -------------------------------------------------------
# -------------------------------------------------------



# -------------------------------------------------------
# -------------------------------------------------------
# Begin industrial model pf
# -------------------------------------------------------
# -------------------------------------------------------


function industrial_dyn_powerflow(
    netd,
    stateDiffCache,
    global_pf_param,
    industrial_model_dyn_pf_up_para,
    industrial_model_idq_pf_cal;
    maxiter=1e5,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    with_δ_ed_eq = true )

    pf_net_param, sd_pf_views, mismatch =
        global_pf_param

    _, _, _, _, pf_Idx, pf_and_dyn_idx_and_Idx, pf_net_misc = pf_net_param

    # pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx, pf_net_misc = pf_net_param

    (slack_bus_idx,
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
    
    (working_vh_θh_view,
     nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     δ_ω_ed_dash_eq_dash_view) =
         sd_pf_views

    (nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_ω_ed_dash_eq_dash_Idxs_in_state,
     nodes_δ_ed_dash_eq_dash_Idxs_in_state) =
         industrial_model_dyn_pf_up_para

   # idq_wt_pad_view, gen_idx = industrial_model_idq_pf_cal
    
    # δ_ω_ed_dash_eq_dash_view .= get_nodes_δ_ω_ed_dash_eq_dash_view( stateDiffCache, netd )

    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        δ_ω_ed_dash_eq_dash_view,
        stateDiffCache,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )
    
    uh_state =
        stateDiffCache[ur_ui_idx][ ur_IDX ] +
        im * stateDiffCache[ur_ui_idx][ui_IDX] 
    
    x0_vh_θh = [abs.(uh_state)...;
                angle.(uh_state)...]
    
    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]
    red_vh_θh_0      = [ red_vh_θh_0_view; ]

    industrial_model_dynamic_power_balance_powerflow(
        x0_vh_θh,
        mismatch,
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view, Iinj_view,
          δ_ω_ed_dash_eq_dash_view ),
        pf_net_param,
        industrial_model_idq_pf_cal;
        maxiter=maxiter,
        ftol=ftol,
        xtol=xtol,
        with_δ_ed_eq = with_δ_ed_eq )

    return nothing
    
end

# -------------------------------------------------------
# -------------------------------------------------------
# End industrial model pf
# -------------------------------------------------------
# -------------------------------------------------------
 
function power_balance_reduced_dict_powerflow(
    netd;
    with_δ_ed_eq =
        false )

    # dynamics_case = new_case_IEEE_9_Bus_dynamic_plant_SM_v6_P
  
    # edges_Ybr_cal, edges_orientation, Ynet, nodes_node_idx_and_incident_edges_other_node_idx = get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )

    (edges_Ybr_cal,
     Ynet,
     edges_orientation,
     nodes_node_idx_and_incident_edges_other_node_idx) =
        get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------

    Ybus =  get_Ybus_from_Ynet(
        Ynet,
        nodes_node_idx_and_incident_edges_other_node_idx)
      
    #------------------------------------------

    slack_bus_idx =
        get_components_slack_Idx(
            netd.nodes )

    slack_vh  =
        get_components_slack_vh(
            netd.nodes )

    gens_Idx_and_vh =
        get_generators_Idx_and_vh(
            netd.nodes )

    non_slack_gens_Idx_and_vh =
        get_non_slack_generators_Idx_and_vh(
            netd.nodes )

    slack_ur_ui_Idx_in_state =
        get_components_slack_ur_ui_Idx_in_state(
            netd.nodes )

    non_slack_ur_ui_Idx_in_state =
        get_components_no_slack_ur_ui_Idx_in_state(
            netd.nodes )
    
    ur_ui_Idx_in_state =
        get_components_ur_ui_Idx_in_state(
            netd.nodes )
    
    nodes_u_Idx =
        get_components_ur_ui_Idx_in_state(
            netd.nodes )

    #------------------------------------------        
    
    gens_idx = first.(gens_Idx_and_vh)
    
    gens_vh  = last.(gens_Idx_and_vh)

    #------------------------------------------        
        
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   = [length(ur_idx), length(ui_idx)]
    
    ur_ui_offset = create_offsets(ur_ui_dims)
    
    ur_ui_IDX    = create_idxs(
        ur_ui_offset, ur_ui_dims)
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    red_vh_θh_idx =
        [setdiff(vh_IDX, gens_idx)...;
         setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ----------------------------------------------------
    
    load_trans_nodes_Idx_and_vlimits =
        get_load_trans_nodes_Idx_and_vlimits(
            netd.nodes )
    
    # ----------------------------------------------------

    syms = get_network_vars_labels(netd)

    # state = zeros(length(syms))

    state_x0 = zeros(length(syms))

    # state_x0 = view(state, 1:length(state))

    state = view(state_x0, 1:length(state_x0))

    # ----------------------------------------------------
    # for pf flat start, ur = vh = 1,  ui = θh = 0

    state[ur_idx] .= ones(  length( vh_IDX ))
    state[ui_idx] .= zeros( length( θh_IDX ))

    # ----------------------------------------------------

    x  = similar(state_x0)
    
    nodes_idx_and_δ_ed_dash_eq_dash_view =
        get_nodes_idx_and_δ_ed_dash_eq_dash_view(
            x, netd )

    nodes_idx_and_δ_ω_ed_dash_eq_dash_view =
        get_nodes_idx_and_δ_ω_ed_dash_eq_dash_view(
            x, netd )
    
    # ----------------------------------------------------
    
    nodes_u_view  =
        [view(x, nodes_u_Idx[Ind])
         for Ind in
             collect(1:length(nodes_u_Idx))]

    pf_state = [ state_x0; ]
    
    nodes_pf_U_view  =
        [ view(pf_state , nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ] 

    
    #------------------------------------------
    
    uh_state_x0 =
        state[ur_ui_idx][ ur_IDX ] +
        im * state[ur_ui_idx][ ui_IDX ]
    
    x0_ur_ui =
        [state[ur_ui_idx][ ur_IDX ]...;
         state[ur_ui_idx][ ui_IDX ]...]
        
    x0_vh_θh =
        [abs.(uh_state_x0)...;
         angle.(uh_state_x0)...]

    working_vh_θh_view =
        view(x0_vh_θh,
             1:length(x0_vh_θh))

    #------------------------------------------
    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii
    #------------------------------------------

    Inet = zeros( ComplexF64, length(uh_state_x0))

    # [ Inet ] # [ view(Inet, 1:length(Inet)) ]
    
    Inet_view =  view( Inet, 1:length(Inet))


    Iinj = zeros( ComplexF64, length(uh_state_x0))

    Iinj_view =  view( Iinj, 1:length(Iinj))
    

    # ----------------------------------------------------

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list =
                [ :ra, :X_d, :X_q,
                  :X_d_dash, :X_q_dash ] )

    ra_Xd_dash_Xq_dash_view  =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list =
                [ :ra, :X_d_dash,
                  :X_q_dash  ] )
        
    #------------------------------------------        

    ra_Xd_Xq_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d, :X_q ] )
        
    P_Q_nodes_view =
        get_components_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )
    
    P_Q_gens_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )

    P_Q_gens_loc_load_view =
        get_gens_loc_load_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :loc_P, :loc_Q ] )

    P_Q_non_gens_view =
        get_non_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )

    #------------------------------------------
    
    pf_net = (
        Ybus,
        Ynet,
        nodes_node_idx_and_incident_edges_other_node_idx,
        edges_Ybr_cal,
        edges_orientation )
    
    pf_idx_and_state =
        ( slack_vh,
          gens_vh,
          gens_Idx_and_vh,
          non_slack_gens_Idx_and_vh,
          slack_ur_ui_Idx_in_state,
          non_slack_ur_ui_Idx_in_state,
          ur_ui_Idx_in_state)    

    pf_views =
        (working_vh_θh_view,
         nodes_u_view,
         nodes_pf_U_view,
         nodes_idx_and_δ_ω_ed_dash_eq_dash_view,
         ra_Xd_Xq_view,
         ra_Xd_dash_Xq_dash_view,
         ra_Xd_Xq_Xd_dash_Xq_dash_view,
         P_Q_nodes_view,
         P_Q_gens_view,
         P_Q_gens_loc_load_view,
         P_Q_non_gens_view,
         Inet_view,
         Iinj_view )

    pf_limits =
        (load_trans_nodes_Idx_and_vlimits, )

    pf_Idx =
        ( slack_bus_idx,
          nodes_u_Idx,
          gens_idx,
          ur_IDX,
          ui_IDX,
          vh_IDX,
          θh_IDX,
          red_vh_θh_idx,
          ur_idx,
          ui_idx,
          ur_ui_idx )
    
    pf_param =
        ( pf_net,
          pf_idx_and_state,
          pf_views,
          pf_limits,
          pf_Idx )

    dyn_pf_param =
        ( nodes_pf_U_view,
          Inet_view,
          Iinj_view )

    # ---------------------------------------------------
    # ---------------------------------------------------
    # Results 
    # ---------------------------------------------------
    # --------------------------------------------------  
    
    sol_with_vh_θh =
        power_balance_powerflow_with_vh_θh!(
            x0_vh_θh,
            pf_param;
            with_δ_ed_eq =
                with_δ_ed_eq )
    
    red_vh_θh =
        sol_with_vh_θh.zero

    working_vh_θh_view[ red_vh_θh_idx ] .= red_vh_θh

    uh = working_vh_θh_view[ vh_IDX ] .*
        exp.( im *  working_vh_θh_view[ θh_IDX ])

    # --------------------------------------------------
      
    branches_name  = collect( keys( netd.edges ))
    
    nodes_name = collect( keys( netd.nodes ))    

    nodes_branches_names = ( nodes_name,  branches_name )
    
    # --------------------------------------------------

    S_gens =
        x_from_xr_xi.(P_Q_gens_view)

    S_gens_loc_load =
        x_from_xr_xi.(P_Q_gens_loc_load_view)

    S_non_gens =
        x_from_xr_xi.(P_Q_non_gens_view)

    # ------------------------------------------------- 
    
    Inet_inj = [ sum(
        [ ynj * vj
          for (ynj, vj) in
              zip(
                  Y_bus_vec,
                  uh[nth_node_idx_and_adj_nodes_idx])])
                 for (Y_bus_vec,
                      nth_node_idx_and_adj_nodes_idx) in zip( Ynet, nodes_node_idx_and_incident_edges_other_node_idx ) ]
    
    Iinj = Inet_inj +
        (conj.( S_non_gens ))./ ( conj.( uh )) +
        (conj.( S_gens_loc_load ))./ (conj.( uh )) 

    Sbus_n  = uh .* conj.( Inet_inj )

    GenSinj = Sbus_n + S_non_gens + S_gens_loc_load
    
    # the matrix in view y_π, is extracted with [1]
    
    Ifrom_Ito =
        [ y_π * [x_from_xr_xi(nodes_pf_U_view[orient[1]]),
                 x_from_xr_xi(nodes_pf_U_view[orient[2]])]
          for (y_π, orient ) in
              zip(edges_Ybr_cal,
                  edges_orientation)]

    If = first.(Ifrom_Ito)
    
    It = last.(Ifrom_Ito)

    Ibranches = If + It
    
    # ---------------------------------------------

    bus_dict_Iinj = OrderedDict(
        name => [ih_r, ih_i]
        for (name, ih_r, ih_i) in
            zip(nodes_name,
                real.(Iinj), imag.(Iinj)))

    bus_dict_Inet_inj = OrderedDict(
        name => [ih_r, ih_i]
        for (name, ih_r, ih_i) in
            zip(nodes_name, real.(Inet_inj),
                imag.(Inet_inj))) 

    branch_dict_init = OrderedDict(
        name => ( real(i_f), imag(i_f),
                  real(i_t), imag(i_t) )
        for (name, i_f, i_t) in
            zip(branches_name, If, It))

    bus_dict_init = OrderedDict(
        name => (vh, θh, ph, qh,
                 ih_r, ih_i, pg, qg,
                 ig_r, ig_i )
        for (name,vh,θh,ph,qh,ih_r,
             ih_i,pg,qg,ig_r,ig_i) in
            zip(nodes_name,
                abs.(uh),
                angle.(uh),
                real.(Sbus_n),
                imag.(Sbus_n),
                real.(Inet_inj),
                imag.(Inet_inj),
                real.(GenSinj),
                imag.(GenSinj),
                real.(Iinj),
                imag.(Iinj)) )

    #-------------------------------------------------

    return ( dyn_pf_param,
            pf_param,
            Dict("Iinj" => Iinj,
                 "Inet_inj" => Inet_inj,
                 "Sg" => S_gens,
                 "Sd" => S_non_gens,
                 "Sloc" => S_gens_loc_load,
                 "Vm" => abs.(uh),
                 "Vθ" => angle.(uh),
                 "Vbus" => uh,
                 "Ibranches" => Ibranches,
                 "Ybus" => Ybus ,
                 "Sbus" => Sbus_n,
                 "GenSinj" => GenSinj,
                 "dict_init" =>Dict(
                     "bus_dict_init" =>
                         bus_dict_init,
                     "branch_dict_init" =>
                         branch_dict_init,
                     "bus_dict_Iinj" =>
                         bus_dict_Iinj,
                     "bus_dict_Inet_inj" =>
                         bus_dict_Inet_inj),
                 "Ifrom" => If, "Ito" => It,
                 "Ynet"=>  Ynet,
                 "nodes_branches_names" =>
                     nodes_branches_names))

end

# ------------------------------------------------------

function power_balance_reduced_dict_powerflow(
    ; dynamics_case = file )    

    return  power_balance_reduced_dict_powerflow(
        NetworkData( dynamics_case()... ) )
    
end

# ------------------------------------------------------
