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

#-------------------------------------------------------
#-------------------------------------------------------
