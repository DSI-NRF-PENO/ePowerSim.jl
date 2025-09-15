# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


#---------------------------------------------------
####################################################
#---------------------------------------------------


# """
#     fdiff_derivatives(f::Function) -> Tuple{Function,Function}

# Return a tuple of functions that evaluate the gradient and Hessian of `f` using
# ForwardDiff.jl.
# """
# function fdiff_derivatives(f::Function)
#     function ∇f(g::AbstractVector{T}, x::Vararg{T,N}) where {T,N}
#         ForwardDiff.gradient!(g, y -> f(y...), collect(x))
#         return
#     end
#     function ∇²f(H::AbstractMatrix{T}, x::Vararg{T,N}) where {T,N}
#         h = ForwardDiff.hessian(y -> f(y...), collect(x))
#         for i in 1:N, j in 1:i
#             H[i, j] = h[i, j]
#         end
#         return
#     end
#     return ∇f, ∇²f
# end


# https://github.com/JuliaStats/PDMats.jl/

# https://juliastats.org/Distributions.jl/stable/multivariate/

function get_active_power_demand_scenario(
    base_Pd;
    scale_factor = 1.0,
    wt_cov = false,
    cov_mat = nothing )
    
    μ_vec = copy(base_Pd )

    # Random.seed!(123)

    if wt_cov == false
        
        σ_vec = ones(length(base_Pd))
    
        Σ_cov = PDiagMat(σ_vec)
        
    else
        
        Σ_cov = cov_mat
        
    end
    

    demand_distribution = MvNormal(μ_vec, Σ_cov)

    scaled_demand =
        scale_factor .* rand(demand_distribution)

    return base_Pd + scaled_demand

end


function get_power_demand_scenario(
    base_power_demand;
    scale_factor = 1.0,
    wt_cov = false,
    cov_mat = nothing )
    
    μ_vec = copy(base_power_demand )

    # Random.seed!(123)

    if wt_cov == false
        
        σ_vec = ones(length(base_power_demand))
    
        Σ_cov = PDiagMat(σ_vec)
        
    else
        
        Σ_cov = cov_mat
        
    end
    

    demand_distribution = MvNormal(μ_vec, Σ_cov)

    scaled_demand =
        scale_factor .* rand(demand_distribution)

    return base_power_demand + scaled_demand

end


function get_wind_gens_power_forecast(
    sites_mean_power;
    scale_factor = 1.0,
    wt_cov = false,
    cov_mat = nothing )
    
    μ_vec = copy( sites_mean_power )

    if wt_cov == false
        
        σ_vec = ones(length(sites_mean_power))
    
        Σ_cov = PDiagMat(σ_vec)
        
    else
        
        Σ_cov = cov_mat
        
    end
    
    wind_power_distribution = MvNormal(μ_vec, Σ_cov)

    scaled_power =
        scale_factor .* rand(wind_power_distribution)

    return sites_mean_power + scaled_power

end


function get_solar_gens_power_forecast(
    sites_mean_power;
    scale_factor = 1.0,
    wt_cov = false,
    cov_mat = nothing )
    
    μ_vec = copy( sites_mean_power )

    if wt_cov == false
        
        σ_vec = ones(length(sites_mean_power))
    
        Σ_cov = PDiagMat(σ_vec)
        
    else
        
        Σ_cov = cov_mat
        
    end
    
    solar_power_distribution = MvNormal(μ_vec, Σ_cov)

    scaled_power =
        scale_factor .* rand(solar_power_distribution)

    return sites_mean_power .+ scaled_power

end

#---------------------------------------------------

function thermal_cost_function(pg, c0, c1, c2)
    if pg <= c0
        return pg
    else
        return pg + c1 * (pg - c2)^2
    end
end

function a_gen_cost_fun(pg, c0, c1, c2 )

    return c0 + c1 * pg + c2 * pg^2
end


function a_wind_gen_cost_fun(
    wg, c0, c1, c2)

    return c0 + c1 * wg + c2 * wg^2

end


function a_solar_gen_cost_func(
    sg, c0, c1, c2)

    return c0 + c1 * sg + c2 * sg^2

end


function a_demand_wind_solar_scenario(
    demand::Float64,
    wind::Float64,
    solar::Float64)
    
    return (demand = demand,
            wind = wind,
            solar = solar )
    
end


#---------------------------------------------------

"""
This was taken from JuMP tutoria:

    thermal_cost_function(g)

A user-defined thermal cost function in pure-Julia! You can include nonlinearities, and even things like control flow.

!!! warning
    It's still up to you to make sure that the function has a meaningful
    derivative.
"""
function thermal_cost_function(g)
    if g <= 500
        return g
    else
        return g + 1e-2 * (g - 500)^2
    end
end


"""
This was taken from JuMP tutoria:
"""
function ThermalGenerator(
    min::Float64,
    max::Float64,
    fixed_cost::Float64,
    variable_cost::Float64,)
    
    return (
        min = min,
        max = max,
        fixed_cost = fixed_cost,
        variable_cost = variable_cost,
    )
end

"""
This was taken from JuMP tutoria:
"""
function Scenario(
    demand::Float64,
    wind::Float64)
    
    return (demand = demand,
            wind = wind)
end

"""
This was taken from JuMP tutoria:
"""
function WindGenerator(variable_cost::Float64)

    return (variable_cost = variable_cost,)

end

"""
This was taken from JuMP tutoria:
"""
function scale_generator_cost(
    g,
    scale)

    return ThermalGenerator(
        g.min,
        g.max, g.fixed_cost,
        scale * g.variable_cost)
end

"""
This was taken from JuMP tutoria:
"""
function scale_demand(
    scenario,
    scale)
    return Scenario(
        scale * scenario.demand,
        scenario.wind)
end

#---------------------------------------------------
#---------------------------------------------------

# transpose(t_pg ) * P * t_pg + transpose(Q) * t_pg

function make_obj_quadratic_form(mpc_gencost )

   gens_cost_coeff_ascen =
        [ (coeff_0, coeff_1, coeff_2)
          for (coeff_0, coeff_1, coeff_2) in
              zip(mpc_gencost.c0,
                  mpc_gencost.var"...",
                  mpc_gencost.var"c(n-1)" )]

    P = diagm(last.( gens_cost_coeff_ascen ))

    Q = second.( gens_cost_coeff_ascen )

    function obj_quadratic_form(gens_pg)

        return transpose(gens_pg ) * P * gens_pg +
            transpose(Q) * gens_pg

    end
end

# obj_quadratic_form =
#     make_obj_quadratic_form(mpc_gencost )


function obj_quadratic_form(gens_pg)

    return transpose(gens_pg ) * P * gens_pg +
        transpose(Q) * gens_pg

end

#---------------------------------------------------


function get_opf_wt_generic_system_simulation_parameters(
    net_data_by_components_file;
    components_libs_dir,
    basekV            = 1.0,    
    use_pu_in_PQ      = true,
    opf_use_pu_in_PQ  = false,
    opf_use_pu_in_df  = true,
    line_data_in_pu   = true,
    with_faults       = false,
    fractional_digits = 6,
    pf_alg            = NewtonRaphson(),
    abstol            = 1e-12,

    reltol            = 1e-12)


    
    (; # dynamic_parameters,
      # pf_parameters,
      # states_Idx_syms_wt_functions,
     opf_parameters, ) =
         NamedTupleTools.select(
             get_generic_system_simulation_parameters(
                 net_data_by_components_file;
                 components_libs_dir,
                 basekV           = basekV,    
                 use_pu_in_PQ     = use_pu_in_PQ,
                 opf_use_pu_in_PQ = opf_use_pu_in_PQ,
                 line_data_in_pu  = line_data_in_pu,
                 with_faults      = with_faults,
                 pf_alg           = pf_alg ),
             (# :dynamic_parameters,
              # :pf_parameters,
              # :states_Idx_syms_wt_functions,
              :opf_parameters, ) )
    
            (;baseMVA,
             basekV,
             slack_gens_nodes_idx,
             non_slack_gens_nodes_idx,
             gens_nodes_idx,
             non_gens_nodes_idx,
             gens_nodes_with_loc_loads_idx,
             all_nodes_idx,

             nodes_with_demands_idx,

             dyn_pf_fun_kwd_net_idxs,
             dyn_pf_fun_kwd_n2s_idxs,

             no_nodes,
             no_gens,
             no_edges,

             ds_generic_gens_para,
             gens_vh_slack_θh_para,

             ode_gens_para,
             generic_gens_para,
             opf_generic_gens_para,

             gencost_data,

             gens_cost_coeff_ascen,
             gens_cost_coeff_decen,

             gens_installed_capacity,

             gens_Pmax,
             gens_Pmin,
             gens_Qmax,
             gens_Qmin,

             sch_Pg,
             nodes_Pd,
             nodes_Qd,

             P_Gen_lb,
             P_Gen_ub,
             Q_Gen_lb,
             Q_Gen_ub,

             P_Demand,
             Q_Demand,
             S_Demand,

             branches_fbus,
             branches_tbus,
             branches_r,
             branches_x,
             branches_b,
             branches_ratio,
             branches_angle,
             branches_type,

             edges_orientation,

             incidence_matrix,

             ys,
             y_c,
             y_sh,
             inv_τ,
             θ_shift,
             Y_sh,
             Cf,
             Ct,
             Yf,
             Yt,
             Ybus) =
                 NamedTupleTools.select(
                     opf_parameters,
                     (:baseMVA,
                      :basekV,
                      :slack_gens_nodes_idx,
                      :non_slack_gens_nodes_idx,
                      :gens_nodes_idx,
                      :non_gens_nodes_idx,
                      :gens_nodes_with_loc_loads_idx,
                      :all_nodes_idx,
            
                      :nodes_with_demands_idx,

                      :dyn_pf_fun_kwd_net_idxs,
                      :dyn_pf_fun_kwd_n2s_idxs,

                      :no_nodes,
                      :no_gens,
                      :no_edges,

                      :ds_generic_gens_para,
                      :gens_vh_slack_θh_para,

                      :ode_gens_para,
                      :generic_gens_para,
                      :opf_generic_gens_para,

                      :gencost_data,
                      
                      :gens_cost_coeff_ascen,
                      :gens_cost_coeff_decen,

                      :gens_installed_capacity,

                      :gens_Pmax,
                      :gens_Pmin,
                      :gens_Qmax,
                      :gens_Qmin,

                      :sch_Pg,
                      :nodes_Pd,
                      :nodes_Qd,

                      :P_Gen_lb,
                      :P_Gen_ub,
                      :Q_Gen_lb,
                      :Q_Gen_ub,

                      :P_Demand,
                      :Q_Demand,
                      :S_Demand,

                      :branches_fbus,
                      :branches_tbus,
                      :branches_r,
                      :branches_x,
                      :branches_b,
                      :branches_ratio,
                      :branches_angle,
                      :branches_type,

                      :edges_orientation,

                      :incidence_matrix,

                      :ys,
                      :y_c,
                      :y_sh,
                      :inv_τ,
                      :θ_shift,
                      :Y_sh,
                      :Cf,
                      :Ct,
                      :Yf,
                      :Yt,
                      :Ybus ))

        
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    model = Model(Ipopt.Optimizer)
    
    set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    @operator(model, gen_cost_fun, 4, a_gen_cost_fun)
    
    #---------------------------------------------------

    @variable(
        model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )

    # 
    #---------------------------------------------------

    @variable(model, V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model, [i in 1:no_nodes],
                0.9^2 <= real(V[i])^2 + imag(V[i])^2 <= 1.1^2)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0);

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0);

    @constraint(model, S_G - S_Demand .== V .* conj(Ybus * V))

    #---------------------------------------------------

    P_G = real(S_G)
    
    @objective(
        model,
        Min,
        sum(gen_cost_fun(P_G[i],
                         gens_cost_coeff_ascen[i]...)
            for i in 1:no_gens), );
           
    #---------------------------------------------------

    optimize!(model)

    # assert_is_solved_and_feasible(model)

    solution_summary(model)

    status = termination_status(model)
    
    #---------------------------------------------------

    objval_solution =
        round(objective_value(model);
              digits = fractional_digits)


    println("Status :" * " $(status)")
    
    println("Objective value (feasible solution) :" *
        " $(objval_solution)")
    
    df_opf = opf_use_pu_in_df == true ? DataFrames.DataFrame(
        ; Bus = all_nodes_idx,
        
        # ComplexPowerGen =
        #     round.(value.(S_G);
        #            digits = fractional_digits),
        
        vh =
            round.(abs.(value.(V));
                   digits = fractional_digits),
        
        θh_rad =
            round.(angle.(value.(V));
                   digits = fractional_digits),
        
        θh_deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits = fractional_digits),
    
        Pg =
            real.(round.( (value.(S_G)) ./ baseMVA;
                         digits = fractional_digits)),

        Qg =
            imag.(round.( (value.(S_G))./ baseMVA ;
                          digits = fractional_digits)), ) :
                              DataFrames.DataFrame(
                                  ; Bus = all_nodes_idx,

                                   vh =
                                       round.(abs.(value.(V));
                                              digits = fractional_digits),

                                   θh_rad =
                                       round.(angle.(value.(V));
                                              digits = fractional_digits),

                                   θh_deg =
                                       round.(rad2deg.(angle.(value.(V)));
                                              digits = fractional_digits),

                                   Pg =
                                       real.(round.( (value.(S_G));
                                                    digits = fractional_digits)),

                                   Qg =
                                       imag.(round.( (value.(S_G));
                                                    digits = fractional_digits)), ) 


    return (;df_opf, objval_solution, status)
    
end


#---------------------------------------------------


function get_opf_by_relaxation_wt_generic_system_simulation_parameters(
    net_data_by_components_file;
    components_libs_dir,
    basekV            = 1.0,    
    use_pu_in_PQ      = true,
    opf_use_pu_in_PQ  = false,
    opf_use_pu_in_df  = true,
    line_data_in_pu   = true,
    with_faults       = false,
    fractional_digits = 6,
    pf_alg            = NewtonRaphson(),
    abstol            = 1e-12,

    reltol            = 1e-12,

    tol_gap_rel       = 1e-3,
    tol_feas          = 1e-3,
    tol_ktratio       = 1e-3)


    
    (;# dynamic_parameters,
      # pf_parameters,
      # states_Idx_syms_wt_functions,
     opf_parameters,) =
         NamedTupleTools.select(
             get_generic_system_simulation_parameters(
                 net_data_by_components_file;
                 components_libs_dir,
                 basekV           = basekV ,    
                 use_pu_in_PQ     = use_pu_in_PQ,
                 opf_use_pu_in_PQ = opf_use_pu_in_PQ,
                 line_data_in_pu  = line_data_in_pu,
                 with_faults      = with_faults,
                 pf_alg           =  pf_alg ),
             (# :dynamic_parameters,
              # :pf_parameters,
              # :states_Idx_syms_wt_functions,
              :opf_parameters, ))
            (;baseMVA,
             basekV,
             slack_gens_nodes_idx,
             non_slack_gens_nodes_idx,
             gens_nodes_idx,
             non_gens_nodes_idx,
             gens_nodes_with_loc_loads_idx,
             all_nodes_idx,

             nodes_with_demands_idx,

             dyn_pf_fun_kwd_net_idxs,
             dyn_pf_fun_kwd_n2s_idxs,

             no_nodes,
             no_gens,
             no_edges,

             ds_generic_gens_para,
             gens_vh_slack_θh_para,

             ode_gens_para,
             generic_gens_para,
             opf_generic_gens_para,

             gencost_data,

             gens_cost_coeff_ascen,
             gens_cost_coeff_decen,

             gens_installed_capacity,

             gens_Pmax,
             gens_Pmin,
             gens_Qmax,
             gens_Qmin,

             sch_Pg,
             nodes_Pd,
             nodes_Qd,

             P_Gen_lb,
             P_Gen_ub,
             Q_Gen_lb,
             Q_Gen_ub,

             P_Demand,
             Q_Demand,
             S_Demand,

             branches_fbus,
             branches_tbus,
             branches_r,
             branches_x,
             branches_b,
             branches_ratio,
             branches_angle,
             branches_type,

             edges_orientation,

             incidence_matrix,

             ys,
             y_c,
             y_sh,
             inv_τ,
             θ_shift,
             Y_sh,
             Cf,
             Ct,
             Yf,
             Yt,
             Ybus) =
                 NamedTupleTools.select(
                     opf_parameters,
                     (:baseMVA,
                      :basekV,
                      :slack_gens_nodes_idx,
                      :non_slack_gens_nodes_idx,
                      :gens_nodes_idx,
                      :non_gens_nodes_idx,
                      :gens_nodes_with_loc_loads_idx,
                      :all_nodes_idx,
            
                      :nodes_with_demands_idx,

                      :dyn_pf_fun_kwd_net_idxs,
                      :dyn_pf_fun_kwd_n2s_idxs,

                      :no_nodes,
                      :no_gens,
                      :no_edges,

                      :ds_generic_gens_para,
                      :gens_vh_slack_θh_para,

                      :ode_gens_para,
                      :generic_gens_para,
                      :opf_generic_gens_para,

                      :gencost_data,
                      
                      :gens_cost_coeff_ascen,
                      :gens_cost_coeff_decen,

                      :gens_installed_capacity,

                      :gens_Pmax,
                      :gens_Pmin,
                      :gens_Qmax,
                      :gens_Qmin,

                      :sch_Pg,
                      :nodes_Pd,
                      :nodes_Qd,

                      :P_Gen_lb,
                      :P_Gen_ub,
                      :Q_Gen_lb,
                      :Q_Gen_ub,

                      :P_Demand,
                      :Q_Demand,
                      :S_Demand,

                      :branches_fbus,
                      :branches_tbus,
                      :branches_r,
                      :branches_x,
                      :branches_b,
                      :branches_ratio,
                      :branches_angle,
                      :branches_type,

                      :edges_orientation,

                      :incidence_matrix,

                      :ys,
                      :y_c,
                      :y_sh,
                      :inv_τ,
                      :θ_shift,
                      :Y_sh,
                      :Cf,
                      :Ct,
                      :Yf,
                      :Yt,
                      :Ybus ))
    #---------------------------------------------------


    E(k, n) =
        SparseArrays.sparse([k], [n], 1, no_nodes, no_nodes);
    
    #---------------------------------------------------
    # Relaxations and better objective bounds
    #---------------------------------------------------

    model = Model(Clarabel.Optimizer)

    set_silent(model)
    
    set_attribute(model, "tol_gap_rel", tol_gap_rel)
    set_attribute(model, "tol_feas",    tol_feas)
    set_attribute(model, "tol_ktratio", tol_ktratio)
    
    #---------------------------------------------------

    # coeff of quadratic

    # x' * P * x + Q' * x + C' * ones(length(x))
    
    P =
        sparse(
            gens_nodes_idx,
            gens_nodes_idx,
            last.(gens_cost_coeff_ascen),
            no_nodes,no_nodes )

    Q =
        SparseArrays.sparsevec(
            gens_nodes_idx,
            second.(gens_cost_coeff_ascen),
            no_nodes )
    
    C =
        SparseArrays.sparsevec(
            gens_nodes_idx,
            first.(gens_cost_coeff_ascen),
            no_nodes )
    
    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    # @operator(model, gen_cost_fun, 4, a_gen_cost_fun)

    # @operator(model, oqf_fun, no_gens, obj_quadratic_form)
    
    
    #---------------------------------------------------

    @variable( model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )

    @variable(model,
              W[1:no_nodes, 1:no_nodes] in HermitianPSDCone())

    @variable(model,
              V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model,
                [i in 1:no_nodes],
                0.9^2 <= real(W[i, i]) <= 1.1^2)

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0)

    @constraint(model,
                0.9 <= real(V[slack_gens_nodes_idx[1]]) <= 1.1)

    @constraint(model,
                LinearAlgebra.Hermitian([1 V'; V W]) in HermitianPSDCone())

    # 2 x 2 minor inequalities:
    @constraint(
        model,
        [i in 1:no_nodes],
        [0.5, real(W[i, i]), real(V[i]), imag(V[i])] in RotatedSecondOrderCone()
    )

    @constraint(
        model,
        [i in 1:no_nodes],
        S_G[i] - S_Demand[i] ==
            LinearAlgebra.tr((conj(Ybus) * E(i, i)) * W),
    )

    #---------------------------------------------------

    P_G = real(S_G)

    #---------------------------------------------------

    @objective(
        model,
        Min,

        sum( C[i] + Q[i] * P_G[i] + P[i,i] * P_G[i]^2
            for i in 1:no_nodes )
        ,)

    
    #---------------------------------------------------

    optimize!(model)


    status = termination_status(model)
    
    #---------------------------------------------------

    # assert_is_solved_and_feasible(
    #     model; allow_almost = true)

    sdp_relaxation_lower_bound =
        round(objective_value(model);
              digits = fractional_digits)

    println(
        "Objective value (W & V relax. lower bound): " *
            "$sdp_relaxation_lower_bound",
    )

    #---------------------------------------------------

    W_1 = SparseArrays.sparse(
        round.(value.(W);
               digits = fractional_digits))

    #---------------------------------------------------

    df_opf_by_relaxation = opf_use_pu_in_df  == true ?
        DataFrames.DataFrame(
            ;Bus = all_nodes_idx,
             vh =
                 round.(abs.(value.(V)); digits =
                 fractional_digits),
             θh_rad =
                 round.(angle.(value.(V));
                        digits =
                            fractional_digits),
             θh_deg =
                 round.(rad2deg.(angle.(value.(V)));
                        digits =
                            fractional_digits),
             Pg =
                 real.(round.( (value.(S_G)) ./ baseMVA;
                              digits = fractional_digits)),

             Qg =
                 imag.(round.( (value.(S_G)) ./ baseMVA;
                              digits = fractional_digits)), ) :
                                  DataFrames.DataFrame(
                                      ;Bus = all_nodes_idx,
                                      vh =
                                          round.(abs.(value.(V)); digits =
                                          fractional_digits),
                                      θh_rad =
                                          round.(angle.(value.(V));
                                                 digits =
                                                     fractional_digits),
                                      θh_deg =
                                          round.(rad2deg.(angle.(value.(V)));
                                                 digits =
                                                     fractional_digits),
                                      Pg =
                                          real.(round.(value.(S_G);
                                                       digits = fractional_digits)),

                                      Qg =
                                          imag.(round.(value.(S_G);
                                                       digits = fractional_digits)), )
    
    return (;df_opf_by_relaxation,
            sdp_relaxation_lower_bound,
            status)
    
end


#---------------------------------------------------

function opf(
    case_file;
    fractional_digits=4)

    mpc_baseMVA =
        get_matpower_scalar_as_iobuffer_by_case_file(
            case_file; type_key_string = "mpc_baseMVA" )[1]

    #---------------------------------------------------

    mpc_bus = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.bus" )


    mpc_gencost =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.gencost" )


    mpc_gen =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.gen" )


    mpc_branch =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.branch" )

    #---------------------------------------------------

    slack_gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 3 ]


    non_slack_gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 2 ]

    gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 3 || node_type == 2  ]


    non_gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 1 ]


    gens_nodes_with_loc_loads_idx =
        [a_node for (a_node, node_type, node_Pd) in
             zip( mpc_bus.bus_i, mpc_bus.type, mpc_bus.Pd)
             if (node_type == 3 || node_type == 2) &&
                 (node_Pd != 0.0 || node_Pd != 0) ]


    nodes_with_demands_idx =
        [a_node for (a_node, node_Pd, node_Qd) in
             zip( mpc_bus.bus_i, mpc_bus.Pd, mpc_bus.Qd)
             if  (node_Pd != 0.0 || node_Pd != 0) ||
                 (node_Qd != 0.0 || node_Qd != 0)]

    all_nodes_idx = copy(mpc_bus.bus_i)

    #---------------------------------------------------

    no_nodes = length(mpc_bus.bus_i)
    
    #---------------------------------------------------

    dyn_pf_fun_kwd_net_idxs =
        (;
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx,

         nodes_with_demands_idx ) 

    #---------------------------------------------------

    no_gens = length(gens_nodes_idx )

    #---------------------------------------------------
    
    gens_cost_coeff_ascen =
        [ (coeff_0, coeff_1, coeff_2)
          for (coeff_0, coeff_1, coeff_2) in
              zip(mpc_gencost.c0,
                  mpc_gencost.var"...",
                  mpc_gencost.var"c(n-1)" )]


    gens_cost_coeff_decen =
        [ ( coeff_2, coeff_1,  coeff_0)
          for (coeff_2, coeff_1,  coeff_0) in
              zip( mpc_gencost.var"c(n-1)",
                   mpc_gencost.var"...",
                   mpc_gencost.c0
                   )]
    
    #---------------------------------------------------
    
    # gens_Pmax = mpc_gen.Pmax[ gens_nodes_idx ]

    # gens_Pmin = mpc_gen.Pmin[ gens_nodes_idx ]

    # gens_Qmax = mpc_gen.Qmax[ gens_nodes_idx ]

    # gens_Qmin = mpc_gen.Qmin[ gens_nodes_idx ]

    
    gens_Pmax = mpc_gen.Pmax

    gens_Pmin = mpc_gen.Pmin

    gens_Qmax = mpc_gen.Qmax

    gens_Qmin = mpc_gen.Qmin


    sch_Pg = mpc_gen.Pg
    
    gens_Pmax = [(sch_Pg[idx] == 0.0 || sch_Pg[idx] == 0 ) ?
        sch_Pg[idx] : pg for (idx, pg) in
                 enumerate(gens_Pmax)]
    
    #---------------------------------------------------

    nodes_Pd = mpc_bus.Pd[nodes_with_demands_idx]

    nodes_Qd = mpc_bus.Qd[nodes_with_demands_idx]

    #---------------------------------------------------
    #---------------------------------------------------

    # Real generation:
    # lower (`lb`) and upper (`ub`) bounds

    P_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmin, no_nodes)
    P_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmax, no_nodes)

    # Reactive generation:
    # lower (`lb`) and upper (`ub`) bounds

    Q_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmin, no_nodes)

    Q_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmax, no_nodes)

    # Power demand levels
    # (real, reactive, and complex form)

    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Pd , no_nodes)

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Qd , no_nodes)

    S_Demand = P_Demand + im * Q_Demand

    #---------------------------------------------------

    branch_data  =
        DataFrame(Any[ mpc_branch.fbus,
                    mpc_branch.tbus,
                    mpc_branch.r,
                    mpc_branch.x,
                    mpc_branch.b ],
                  [:F_BUS, :T_BUS,
                   :BR_R, :BR_X,
                   :BR_Bc] )

    #---------------------------------------------------

    no_edges = size(branch_data, 1)

    #---------------------------------------------------

    incidence_matrix  =
        SparseArrays.sparse(
            branch_data.F_BUS,
            1:no_edges, 1, no_nodes, no_edges) +

                SparseArrays.sparse(
                    branch_data.T_BUS,
                    1:no_edges, -1, no_nodes, no_edges)

    #---------------------------------------------------
    # Matpower
    #---------------------------------------------------

    ys =
         mpc_baseMVA ./ (mpc_branch.r + im * mpc_branch.x) 
       

    y_c =
        1 / 2 * (im *  mpc_branch.b) *
        mpc_baseMVA

    y_sh =  (mpc_bus.Gs .+ im *  mpc_bus.Bs)/mpc_baseMVA

    inv_τ =
        [ (a_ratio == 0.0 || a_ratio == 0) ? 1.0 : 1/a_ratio
              for a_ratio in
                  mpc_branch.ratio  ]
    
    θ_shift = mpc_branch.angle
    
    Yff = (ys +  y_c) .* (inv_τ).^2

    Ytf = -ys .* (inv_τ ./ exp.(im * θ_shift))

    Yft = -ys .* (inv_τ ./ exp.(-im * θ_shift))

    Ytt = (ys +  y_c)

    Y_sh = SparseArrays.spdiagm(y_sh)

    Cf = SparseArrays.sparse(
        1:no_edges,
        branch_data.F_BUS,
             1,  no_edges, no_nodes) 

    Ct = SparseArrays.sparse(
        1:no_edges,
        branch_data.T_BUS,
             1,  no_edges, no_nodes) 
    
    Yf = SparseArrays.spdiagm( Yff ) * Cf +
        SparseArrays.spdiagm( Yft ) * Ct

    Yt = SparseArrays.spdiagm( Ytf ) * Cf +
        SparseArrays.spdiagm( Ytt ) * Ct

    Ybus = Cf' * Yf + Ct' * Yt + Y_sh
        
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    model = Model(Ipopt.Optimizer)
    
    set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    @operator(model, gen_cost_fun, 4, a_gen_cost_fun)
    
    #---------------------------------------------------

    @variable(
        model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )

    #---------------------------------------------------

    @variable(model, V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model, [i in 1:no_nodes],
                0.9^2 <= real(V[i])^2 + imag(V[i])^2 <= 1.1^2)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0);

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0);

    @constraint(model, S_G - S_Demand .== V .* conj(Ybus * V))

    #---------------------------------------------------

    P_G = real(S_G)
    
    @objective(
        model,
        Min,
        sum(gen_cost_fun(P_G[i],
                         gens_cost_coeff_ascen[i]...)
            for i in 1:no_gens), );
           
    #---------------------------------------------------

    optimize!(model)

    # assert_is_solved_and_feasible(model)

    solution_summary(model)

    status = termination_status(model)
    
    #---------------------------------------------------

    objval_solution =
        round(objective_value(model);
              digits = fractional_digits)


    println("Status :" * " $(status)")
    
    println("Objective value (feasible solution) :" *
        " $(objval_solution)")
    
    df_opf = DataFrames.DataFrame(
        ; Bus = all_nodes_idx,
        
        # ComplexPowerGen =
        #     round.(value.(S_G);
        #            digits = fractional_digits),

        P =
            real.(round.(value.(S_G);
                         digits = fractional_digits)),

        Q =
            imag.(round.(value.(S_G);
                         digits = fractional_digits)),
        
        VoltageMagnitude =
            round.(abs.(value.(V));
                   digits = fractional_digits),
        
        VoltageAngle_Rad =
            round.(angle.(value.(V));
                   digits = fractional_digits),
        
        VoltageAngle_Deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits = fractional_digits), )


    return (;df_opf, objval_solution, status)
end


function opf_by_relaxation(
    case_file;
    fractional_digits = 4 )

    # case_file = case_9_file
    
    mpc_baseMVA =
        get_matpower_scalar_as_iobuffer_by_case_file(
            case_file; type_key_string = "mpc_baseMVA" )[1]

    #---------------------------------------------------


    mpc_bus = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.bus" )


    mpc_gencost =
        get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.gencost" )


    mpc_gen = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.gen" )


    mpc_branch = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.branch" )

    #---------------------------------------------------

    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_mpc( mpc_bus  )
    
    slack_gens_nodes_idx =
        net_nodes_type_idxs.slack_gens_nodes_idx

    non_slack_gens_nodes_idx =
        net_nodes_type_idxs.non_slack_gens_nodes_idx
    
    gens_nodes_idx =
        net_nodes_type_idxs.gens_nodes_idx

    non_gens_nodes_idx =
        net_nodes_type_idxs.non_gens_nodes_idx

    gens_nodes_with_loc_loads_idx =
        net_nodes_type_idxs.gens_nodes_with_loc_loads_idx

    all_nodes_idx =
        net_nodes_type_idxs.all_nodes_idx

    nodes_with_demands_idx =
        net_nodes_type_idxs.nodes_with_demands_idx

    #---------------------------------------------------
    
    no_nodes = length(mpc_bus.bus_i)

    #---------------------------------------------------

    E(k, n) =
        SparseArrays.sparse([k], [n], 1, no_nodes, no_nodes);
    
    #---------------------------------------------------

    dyn_pf_fun_kwd_net_idxs =
        (;
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx,
         nodes_with_demands_idx ) 

    #---------------------------------------------------

    no_gens = length(gens_nodes_idx )

    #---------------------------------------------------
    
    gens_cost_coeff_ascen =
        [ (coeff_0, coeff_1, coeff_2)
          for (coeff_0, coeff_1, coeff_2) in
              zip(mpc_gencost.c0,
                  mpc_gencost.var"...",
                  mpc_gencost.var"c(n-1)" )]


    gens_cost_coeff_decen =
        [ ( coeff_2, coeff_1,  coeff_0)
          for (coeff_2, coeff_1,  coeff_0) in
              zip( mpc_gencost.var"c(n-1)",
                   mpc_gencost.var"...",
                   mpc_gencost.c0
                   )]
    
    #---------------------------------------------------

    # gens_Pmax = mpc_gen.Pmax[ gens_nodes_idx ]

    # gens_Pmin = mpc_gen.Pmin[ gens_nodes_idx ]

    # gens_Qmax = mpc_gen.Qmax[ gens_nodes_idx ]

    # gens_Qmin = mpc_gen.Qmin[ gens_nodes_idx ]


    gens_Pmax = mpc_gen.Pmax

    gens_Pmin = mpc_gen.Pmin

    gens_Qmax = mpc_gen.Qmax

    gens_Qmin = mpc_gen.Qmin

    sch_Pg = mpc_gen.Pg
    
    gens_Pmax =
        [(sch_Pg[idx] == 0.0 || sch_Pg[idx] == 0 ) ?
        sch_Pg[idx] : pg for (idx, pg) in
                 enumerate(gens_Pmax)]
    #---------------------------------------------------

    nodes_Pd = mpc_bus.Pd[nodes_with_demands_idx]

    nodes_Qd = mpc_bus.Qd[nodes_with_demands_idx]

    #---------------------------------------------------
    #---------------------------------------------------

    # Real generation:
    # lower (`lb`) and upper (`ub`) bounds

    P_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmin, no_nodes)
    P_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmax, no_nodes)

    # Reactive generation:
    # lower (`lb`) and upper (`ub`) bounds

    Q_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmin, no_nodes)

    Q_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmax, no_nodes)

    # Power demand levels
    # (real, reactive, and complex form)

    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Pd , no_nodes)

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Qd , no_nodes)

    S_Demand = P_Demand + im * Q_Demand

    #---------------------------------------------------

    branch_data  =
        DataFrame(Any[ mpc_branch.fbus,
                       mpc_branch.tbus,
                       mpc_branch.r,
                       mpc_branch.x,
                       mpc_branch.b,
                       mpc_branch.ratio,
                       mpc_branch.angle],
                  [:F_BUS, :T_BUS,
                   :BR_R, :BR_X,
                   :BR_Bc, :BR_ratio,
                   :BR_angle ] )

    #---------------------------------------------------

    no_edges = size(branch_data)[1]

    #---------------------------------------------------
    # Matpower
    #---------------------------------------------------

    ys =
         mpc_baseMVA ./ (mpc_branch.r + im * mpc_branch.x) 
       

    y_c =
        1 / 2 * (im *  mpc_branch.b) *
        mpc_baseMVA

    y_sh =  (mpc_bus.Gs .+ im *  mpc_bus.Bs)/mpc_baseMVA

    inv_τ =
        [ (a_ratio == 0.0 || a_ratio == 0) ? 1.0 : 1/a_ratio
              for a_ratio in
                  mpc_branch.ratio  ]
    
    θ_shift = mpc_branch.angle
    
    Yff = (ys +  y_c) .* (inv_τ).^2

    Ytf = -ys .* (inv_τ ./ exp.(im * θ_shift))

    Yft = -ys .* (inv_τ ./ exp.(-im * θ_shift))

    Ytt = (ys +  y_c)

    Y_sh = SparseArrays.spdiagm(y_sh)

    Cf = SparseArrays.sparse(
        1:no_edges,
        branch_data.F_BUS,
             1,  no_edges, no_nodes) 

    Ct = SparseArrays.sparse(
        1:no_edges,
        branch_data.T_BUS,
             1,  no_edges, no_nodes) 
    
    Yf = SparseArrays.spdiagm( Yff ) * Cf +
        SparseArrays.spdiagm( Yft ) * Ct

    Yt = SparseArrays.spdiagm( Ytf ) * Cf +
        SparseArrays.spdiagm( Ytt ) * Ct

    Ybus = Cf' * Yf + Ct' * Yt + Y_sh
    
    #---------------------------------------------------
    # Relaxations and better objective bounds
    #---------------------------------------------------

    model = Model(Clarabel.Optimizer)

    set_attribute(model, "tol_gap_rel", 1e-3)
    set_attribute(model, "tol_feas", 1e-3)
    set_attribute(model, "tol_ktratio", 1e-3)

    #---------------------------------------------------

    # coeff of quadratic

    # x' * P * x + Q' * x + C' * ones(length(x))
    
    P =
        sparse(
            gens_nodes_idx,
            gens_nodes_idx,
            last.(gens_cost_coeff_ascen),
            no_nodes,no_nodes )

    Q =
        SparseArrays.sparsevec(
            gens_nodes_idx,
            second.(gens_cost_coeff_ascen),
            no_nodes )
    
    C =
        SparseArrays.sparsevec(
            gens_nodes_idx,
            first.(gens_cost_coeff_ascen),
            no_nodes )
    
    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    # @operator(model, gen_cost_fun, 4, a_gen_cost_fun)

    # @operator(model, oqf_fun, no_gens, obj_quadratic_form)
    
    
    #---------------------------------------------------

    @variable( model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )

    @variable(model,
              W[1:no_nodes, 1:no_nodes] in HermitianPSDCone())

    @variable(model,
              V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model,
                [i in 1:no_nodes],
                0.9^2 <= real(W[i, i]) <= 1.1^2)

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0)

    @constraint(model,
                0.9 <= real(V[slack_gens_nodes_idx[1]]) <= 1.1)

    @constraint(model,
                LinearAlgebra.Hermitian([1 V'; V W]) in HermitianPSDCone())

    # 2 x 2 minor inequalities:
    @constraint(
        model,
        [i in 1:no_nodes],
        [0.5, real(W[i, i]), real(V[i]), imag(V[i])] in RotatedSecondOrderCone()
    )

    @constraint(
        model,
        [i in 1:no_nodes],
        S_G[i] - S_Demand[i] ==
            LinearAlgebra.tr((conj(Ybus) * E(i, i)) * W),
    )

    #---------------------------------------------------

    P_G = real(S_G)

    #---------------------------------------------------

    @objective(
        model,
        Min,

        sum( C[i] + Q[i] * P_G[i] + P[i,i] * P_G[i]^2
            for i in 1:no_nodes )
        ,)

    
    #---------------------------------------------------

    optimize!(model)


    status = termination_status(model)
    
    #---------------------------------------------------

    # assert_is_solved_and_feasible(
    #     model; allow_almost = true)

    sdp_relaxation_lower_bound =
        round(objective_value(model);
              digits = fractional_digits)

    println(
        "Objective value (W & V relax. lower bound): " *
            "$sdp_relaxation_lower_bound",
    )

    #---------------------------------------------------

    W_1 = SparseArrays.sparse(
        round.(value.(W);
               digits = fractional_digits))

    #---------------------------------------------------

    df_opf_by_relaxation = DataFrames.DataFrame(
        ;Bus = all_nodes_idx,
        Magnitude =
            round.(abs.(value.(V)); digits =
            round_fractional_digits),
        AngleRad =
            round.(angle.(value.(V));
                   digits =
                       fractional_digits),
        AngleDeg =
            round.(rad2deg.(angle.(value.(V)));
                   digits =
                       fractional_digits), )
    
    return (;df_opf_by_relaxation,
            sdp_relaxation_lower_bound,
            status )
        
end



function opf(
    nodes_Pd,
    nodes_Qd;
    opf_net_optimisation_parameters =
        opf_net_optimisation_parameters,
    round_digits = 4 )


    (;no_nodes,
     no_gens,
     
     P_Gen_lb,
     P_Gen_ub,
     Q_Gen_lb,
     Q_Gen_ub,

     gens_nodes_idx,
     nodes_with_demands_idx,
     
     slack_gens_nodes_idx,
     dyn_pf_fun_kwd_net_idxs,

     # S_Demand,
     Ybus,

     gens_cost_coeff_ascen) =
        NamedTupleTools.select(
            opf_net_optimisation_parameters,
            (:no_nodes,
             :no_gens,
             
             :P_Gen_lb,
             :P_Gen_ub,
             
             :Q_Gen_lb,
             :Q_Gen_ub,

             :gens_nodes_idx,
             :nodes_with_demands_idx,

             :slack_gens_nodes_idx,
             :dyn_pf_fun_kwd_net_idxs,

             # :S_Demand,
             :Ybus,

             :gens_cost_coeff_ascen))

    
    #---------------------------------------------------
    
    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            nodes_Pd,
            no_nodes)

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            nodes_Qd,
            no_nodes)

    S_Demand = P_Demand + im * Q_Demand
    
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    model = Model(Ipopt.Optimizer)
    
    set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    @operator(model, gen_cost_fun, 4, a_gen_cost_fun)
    
    #---------------------------------------------------

    @variable(
        model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )

    #---------------------------------------------------

    @variable(model, V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model, [i in 1:no_nodes],
                0.9^2 <= real(V[i])^2 + imag(V[i])^2 <= 1.1^2)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0);

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0);

    @constraint(model, S_G - S_Demand .== V .* conj(Ybus * V))

    #---------------------------------------------------

    P_G = real(S_G)

    # Q_G = imag(S_G)

    # #---------------------------------------------------

    # @constraint(model, [i in 1:no_nodes],
    #             P_Gen_lb[i] <= P_G[i]  <= P_Gen_ub[i])

    # @constraint(model, [i in 1:no_nodes],
    #             Q_Gen_lb[i] <= Q_G[i]  <= Q_Gen_ub[i])

    #---------------------------------------------------
    
    @objective(
        model,
        Min,
        sum(gen_cost_fun(P_G[i],
                         gens_cost_coeff_ascen[i]...)
            for i in 1:no_gens), );
           
    #---------------------------------------------------

    optimize!(model)

    # assert_is_solved_and_feasible(model)

    solution_summary(model)

    status = termination_status(model)
    
    #---------------------------------------------------

    objval_solution =
        round(objective_value(model);
              digits = round_digits)

    println("Objective value (feasible solution) :" *
        " $(objval_solution)")

    return (
        ; status = status,
        
        objval_solution =
            objval_solution,
        
        Gen_S =
            round.(value.(S_G)[gens_nodes_idx];
                   digits = round_digits),
        Gen_P =
            round.(real.(value.(S_G)[gens_nodes_idx]);
                   digits = round_digits),
        Gen_Q =
            round.(imag.(value.(S_G)[gens_nodes_idx]);
                   digits = round_digits),
        
        V_mag =
            round.(abs.(value.(V));
                   digits = round_digits),
        
        V_ang_rad =
            round.(angle.(value.(V));
                   digits =
                       round_digits),

        Gen_V_mag =
            round.(abs.(value.(V)[gens_nodes_idx]);
                   digits = round_digits),

        Gen_ang_rad =
            round.(angle.(value.(V)[gens_nodes_idx]);
                   digits = round_digits),
        
        V_ang_deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits =
                       round_digits), )

end



function opf_by_relaxation(
    nodes_Pd,
    nodes_Qd;
    opf_net_optimisation_parameters =
        opf_net_optimisation_parameters,
    round_digits = 4)


    (;no_nodes,
     no_gens,
     gens_nodes_idx,
     
     P_Gen_lb,
     P_Gen_ub,
     
     Q_Gen_lb,
     Q_Gen_ub,

     nodes_with_demands_idx,
     
     slack_gens_nodes_idx,
     dyn_pf_fun_kwd_net_idxs,

     # S_Demand,
     Ybus,

     gens_cost_coeff_ascen) =
        NamedTupleTools.select(
            opf_net_optimisation_parameters,
            (:no_nodes,
             :no_gens,
             :gens_nodes_idx,
             
             :P_Gen_lb,
             :P_Gen_ub,
             
             :Q_Gen_lb,
             :Q_Gen_ub,

             :nodes_with_demands_idx,

             :slack_gens_nodes_idx,
             :dyn_pf_fun_kwd_net_idxs,

             # :S_Demand,
             :Ybus,

             :gens_cost_coeff_ascen))

    
    #---------------------------------------------------
    
    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            nodes_Pd,
            no_nodes)

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            nodes_Qd,
            no_nodes)

    S_Demand = P_Demand + im * Q_Demand
    
    
    #---------------------------------------------------
    # Relaxations and better objective bounds
    #---------------------------------------------------

    model = Model(Clarabel.Optimizer)

    set_attribute(model, "tol_gap_rel", 1e-3)
    set_attribute(model, "tol_feas", 1e-3)
    set_attribute(model, "tol_ktratio", 1e-3)

    #---------------------------------------------------

    # coeff of quadratic

    # x' * P * x + Q' * x + C' * ones(length(x))
    
    P =
        sparse(
            gens_nodes_idx,
            gens_nodes_idx,
            last.(gens_cost_coeff_ascen),
            no_nodes,no_nodes )

    Q =
        SparseArrays.sparsevec(
            gens_nodes_idx,
            second.(gens_cost_coeff_ascen),
            no_nodes )
    
    C =
        SparseArrays.sparsevec(
            gens_nodes_idx,
            first.(gens_cost_coeff_ascen),
            no_nodes )
    
    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    # @operator(model, gen_cost_fun, 4, a_gen_cost_fun)

    # @operator(model, oqf_fun, no_gens, obj_quadratic_form)
        
    #---------------------------------------------------

    E(k, n) =
        SparseArrays.sparse([k], [n], 1, no_nodes, no_nodes);
    
    #---------------------------------------------------

    @variable( model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )

    @variable(model,
              W[1:no_nodes, 1:no_nodes] in HermitianPSDCone())

    @variable(model,
              V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model,
                [i in 1:no_nodes],
                0.9^2 <= real(W[i, i]) <= 1.1^2)

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0)

    @constraint(model,
                0.9 <= real(V[slack_gens_nodes_idx[1]]) <= 1.1)

    @constraint(model,
                LinearAlgebra.Hermitian([1 V'; V W]) in HermitianPSDCone())

    # 2 x 2 minor inequalities:
    @constraint(
        model,
        [i in 1:no_nodes],
        [0.5, real(W[i, i]), real(V[i]), imag(V[i])] in RotatedSecondOrderCone()
    )

    @constraint(
        model,
        [i in 1:no_nodes],
        S_G[i] - S_Demand[i] ==
            LinearAlgebra.tr((conj(Ybus) * E(i, i)) * W),
    )

    #---------------------------------------------------

    P_G = real(S_G)

    #---------------------------------------------------

    @objective(
        model,
        Min,

        sum( C[i] + Q[i] * P_G[i] + P[i,i] * P_G[i]^2
            for i in 1:no_nodes )
        ,)

    
    #---------------------------------------------------

    optimize!(model)

    status = termination_status(model)
    
    #---------------------------------------------------

    # assert_is_solved_and_feasible(
    #     model; allow_almost = true)

    sdp_relaxation_lower_bound =
        round(objective_value(model);
              digits = round_digits)

    println(
        "Objective value (W & V relax. lower bound): " *
            "$sdp_relaxation_lower_bound",
    )

    #---------------------------------------------------

    W_1 = SparseArrays.sparse(round.(value.(W);
                                     digits = round_digits))

    #---------------------------------------------------

    return(
        ;status = status,

        objval_solution =
            sdp_relaxation_lower_bound,
        
         Gen_S =
            round.(value.(S_G)[gens_nodes_idx];
                   digits = round_digits),

        Gen_P =
            round.(real.(value.(S_G)[gens_nodes_idx]);
                   digits = round_digits),

        Gen_Q =
            round.(imag.(value.(S_G)[gens_nodes_idx]);
                   digits = round_digits),
        
        V_mag =
            round.(abs.(value.(V));
                   digits = round_digits),
        
        V_ang_rad =
            round.(angle.(value.(V));
                   digits =
                       round_digits),

        Gen_V_mag =
            round.(abs.(value.(V)[gens_nodes_idx]);
                   digits = round_digits),

        Gen_ang_rad =
            round.(angle.(value.(V)[gens_nodes_idx]);
                   digits = round_digits),
        
        V_ang_deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits = round_digits), )
        
end


function opf(
    nodes_Pd,
    nodes_Qd,
    Ynet,
    nodes_idx_with_adjacent_nodes_idx;
    opf_net_optimisation_parameters =
        opf_net_optimisation_parameters,
    round_digits = 4)


    (;no_nodes,
     no_gens,
     
     P_Gen_lb,
     P_Gen_ub,
     Q_Gen_lb,
     Q_Gen_ub,
     
     slack_gens_nodes_idx,
     dyn_pf_fun_kwd_net_idxs,

     # S_Demand,
     # Ybus,

     gens_cost_coeff_ascen) =
        NamedTupleTools.select(
            opf_net_optimisation_parameters,
            (:no_nodes,
             :no_gens,
             
             :P_Gen_lb,
             :P_Gen_ub,
             
             :Q_Gen_lb,
             :Q_Gen_ub,

             :slack_gens_nodes_idx,
             :dyn_pf_fun_kwd_net_idxs,

             # :S_Demand,
             # :Ybus,

             :gens_cost_coeff_ascen))
    

    #---------------------------------------------------
    
    YBus =  get_Ybus_from_Ynet(
        Ynet,
        nodes_idx_with_adjacent_nodes_idx)
    
    
    #---------------------------------------------------
    
    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            nodes_Pd,
            no_nodes)

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            nodes_Qd,
            no_nodes)

    S_Demand = P_Demand + im * Q_Demand
    
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    model = Model(Ipopt.Optimizer)
    
    set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    @operator(model, gen_cost_fun, 4, a_gen_cost_fun)
    
    #---------------------------------------------------

    @variable(
        model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )

    #---------------------------------------------------

    @variable(model, V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model, [i in 1:no_nodes],
                0.9^2 <= real(V[i])^2 + imag(V[i])^2 <= 1.1^2)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0);

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0);

    @constraint(model, S_G - S_Demand .== V .* conj(Ybus * V))

    #---------------------------------------------------

    P_G = real(S_G)

    Q_G = imag(S_G)

    #---------------------------------------------------

    @constraint(model, [i in 1:no_nodes],
                P_Gen_lb[i] <= P_G[i]  <= P_Gen_ub[i])

    @constraint(model, [i in 1:no_nodes],
                Q_Gen_lb[i] <= Q_G[i]  <= Q_Gen_ub[i])

    #---------------------------------------------------
    
    @objective(
        model,
        Min,
        sum(gen_cost_fun(P_G[i],
                         gens_cost_coeff_ascen[i]...)
            for i in 1:no_gens), );
           
    #---------------------------------------------------

    optimize!(model)

    # assert_is_solved_and_feasible(model)

    solution_summary(model)

    status = termination_status(model)
    
    #---------------------------------------------------

    objval_solution =
        round(objective_value(model);
              digits = round_digits)

    println("Objective value (feasible solution) :" *
        " $(objval_solution)")

    return (
        ; status = status,
        
        objval_solution =
            objval_solution,
        
        Gen_S =
            round.(value.(S_G)[gens_nodes_idx];
                   digits = round_digits),
        Gen_P =
            round.(real.(value.(S_G)[gens_nodes_idx]);
                   digits = round_digits),
        Gen_Q =
            round.(imag.(value.(S_G)[gens_nodes_idx]);
                   digits = round_digits),
        
        V_mag =
            round.(abs.(value.(V));
                   digits = round_digits),
        
        V_ang_rad =
            round.(angle.(value.(V));
                   digits =
                       round_digits),

        Gen_V_mag =
            round.(abs.(value.(V)[gens_nodes_idx]);
                   digits = round_digits),

        Gen_ang_rad =
            round.(angle.(value.(V)[gens_nodes_idx]);
                   digits = round_digits),
        
        V_ang_deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits =
                       round_digits), )

end



function opf_by_relaxation(
    nodes_Pd,
    nodes_Qd,
    Ynet,
    nodes_idx_with_adjacent_nodes_idx;
    opf_net_optimisation_parameters =
        opf_net_optimisation_parameters,
    round_digits = 4)


    (;no_nodes,
     no_gens,
     gens_nodes_idx,
     
     P_Gen_lb,
     P_Gen_ub,
     
     Q_Gen_lb,
     Q_Gen_ub,
     
     slack_gens_nodes_idx,
     dyn_pf_fun_kwd_net_idxs,

     # S_Demand,
     # Ybus,

     gens_cost_coeff_ascen) =
        NamedTupleTools.select(
            opf_net_optimisation_parameters,
            (:no_nodes,
             :no_gens,
             :gens_nodes_idx,
             
             :P_Gen_lb,
             :P_Gen_ub,
             
             :Q_Gen_lb,
             :Q_Gen_ub,

             :slack_gens_nodes_idx,
             :dyn_pf_fun_kwd_net_idxs,

             # :S_Demand,
             # :Ybus,

             :gens_cost_coeff_ascen))

    

    #---------------------------------------------------
    
    YBus =  get_Ybus_from_Ynet(
        Ynet,
        nodes_idx_with_adjacent_nodes_idx)
    
    
    #---------------------------------------------------
    
    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            nodes_Pd,
            no_nodes)

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            nodes_Qd,
            no_nodes)

    S_Demand = P_Demand + im * Q_Demand
    
    #---------------------------------------------------
    # Relaxations and better objective bounds
    #---------------------------------------------------

    model = Model(Clarabel.Optimizer)

    set_attribute(model, "tol_gap_rel", 1e-3)
    set_attribute(model, "tol_feas", 1e-3)
    set_attribute(model, "tol_ktratio", 1e-3)

    #---------------------------------------------------

    # coeff of quadratic

    # x' * P * x + Q' * x + C' * ones(length(x))
    
    P =
        sparse(
            gens_nodes_idx,
            gens_nodes_idx,
            last.(gens_cost_coeff_ascen),
            no_nodes,no_nodes )

    Q =
        SparseArrays.sparsevec(
            gens_nodes_idx,
            second.(gens_cost_coeff_ascen),
            no_nodes )
    
    C =
        SparseArrays.sparsevec(
            gens_nodes_idx,
            first.(gens_cost_coeff_ascen),
            no_nodes )
    
    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    # @operator(model, gen_cost_fun, 4, a_gen_cost_fun)

    # @operator(model, oqf_fun, no_gens, obj_quadratic_form)
    
    
    #---------------------------------------------------

    E(k, n) =
        SparseArrays.sparse([k], [n], 1, no_nodes, no_nodes);
    
    #---------------------------------------------------


    @variable( model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )

    @variable(model,
              W[1:no_nodes, 1:no_nodes] in HermitianPSDCone())

    @variable(model,
              V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model,
                [i in 1:no_nodes],
                0.9^2 <= real(W[i, i]) <= 1.1^2)

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0)

    @constraint(model,
                0.9 <= real(V[slack_gens_nodes_idx[1]]) <= 1.1)

    @constraint(model,
                LinearAlgebra.Hermitian([1 V'; V W]) in HermitianPSDCone())

    # 2 x 2 minor inequalities:
    @constraint(
        model,
        [i in 1:no_nodes],
        [0.5, real(W[i, i]), real(V[i]), imag(V[i])] in RotatedSecondOrderCone()
    )

    @constraint(
        model,
        [i in 1:no_nodes],
        S_G[i] - S_Demand[i] ==
            LinearAlgebra.tr((conj(Ybus) * E(i, i)) * W),
    )

    #---------------------------------------------------

    P_G = real(S_G)

    #---------------------------------------------------

    @objective(
        model,
        Min,

        sum( C[i] + Q[i] * P_G[i] + P[i,i] * P_G[i]^2
            for i in 1:no_nodes )
        ,)

    
    #---------------------------------------------------

    optimize!(model)


    status = termination_status(model)
    
    #---------------------------------------------------

    # assert_is_solved_and_feasible(
    #     model; allow_almost = true)

    sdp_relaxation_lower_bound =
        round(objective_value(model);
              digits = round_digits)

    println(
        "Objective value (W & V relax. lower bound): " *
            "$sdp_relaxation_lower_bound",
    )

    #---------------------------------------------------

    W_1 = SparseArrays.sparse(round.(value.(W);
                                     digits = round_digits))

    #---------------------------------------------------

    return(
        ;status = status,

        objval_solution =
            sdp_relaxation_lower_bound,
        
         Gen_S =
            round.(value.(S_G)[gens_nodes_idx];
                   digits = round_digits),

        Gen_P =
            round.(real.(value.(S_G)[gens_nodes_idx]);
                   digits = round_digits),

        Gen_Q =
            round.(imag.(value.(S_G)[gens_nodes_idx]);
                   digits = round_digits),
        
        V_mag =
            round.(abs.(value.(V));
                   digits = round_digits),
        
        V_ang_rad =
            round.(angle.(value.(V));
                   digits =
                       round_digits),

        Gen_V_mag =
            round.(abs.(value.(V)[gens_nodes_idx]);
                   digits = round_digits),

        Gen_ang_rad =
            round.(angle.(value.(V)[gens_nodes_idx]);
                   digits = round_digits),
        
        V_ang_deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits = round_digits), )
        
end


function opf_by_line_loss( case_file )

   (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
            
     nodes_with_demands_idx,

     dyn_pf_fun_kwd_net_idxs,

     no_nodes,
     no_gens,
     no_edges,

     gens_cost_coeff_ascen,
     gens_cost_coeff_decen,

     gens_installed_capacity,

     gens_Pmax,
     gens_Pmin,
     gens_Qmax,
     gens_Qmin,
            
     sch_Pg,
     nodes_Pd,
     nodes_Qd,

     P_Gen_lb,
     P_Gen_ub,
     Q_Gen_lb,
     Q_Gen_ub,

     P_Demand,
     Q_Demand,
     S_Demand,

     branch_fbus,
     branch_tbus,
     branch_r,
     branch_x,
     branch_b,

     edges_orientation,

     branch_data,
 
     Ybus) =
         NamedTupleTools.select(
             get_opf_net_optimisation_parameters(
                 case_file ),
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,

              :nodes_with_demands_idx,

              :dyn_pf_fun_kwd_net_idxs,

              :no_nodes,
              :no_gens,
              :no_edges,

              :gens_cost_coeff_ascen,
              :gens_cost_coeff_decen,

              :gens_installed_capacity,

              :gens_Pmax,
              :gens_Pmin,
              :gens_Qmax,
              :gens_Qmin,

              :sch_Pg,
              :nodes_Pd,
              :nodes_Qd,

              :P_Gen_lb,
              :P_Gen_ub,
              :Q_Gen_lb,
              :Q_Gen_ub,

              :P_Demand,
              :Q_Demand,
              :S_Demand,

              :branch_fbus,
              :branch_tbus,
              :branch_r,
              :branch_x,
              :branch_b,

              :edges_orientation,

              :branch_data,

              :Ybus))

    n2s_gens_nodes =
            OrderedDict( a_node => idx
                         for (idx, a_node) in
                             enumerate(gens_nodes_idx))
    
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    model = Model(Ipopt.Optimizer)
    
    # set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
        
    # get_a_line_loss_r_i(
    #     edge_r, edge_x,    
    #     f_uh_real, f_uh_imag,
    #     t_uh_real, t_uh_imag)
    
    @operator(model, loss_fun, 6, get_a_line_loss_by_ur_ui)
    
    #---------------------------------------------------

    @variable(
        model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )

    #---------------------------------------------------

    @variable(model, V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    P_G = real(S_G)
    
    @constraint(model, [i in non_slack_gens_nodes_idx],
              P_G[i]  == sch_Pg[n2s_gens_nodes[i]] )

    
    @constraint(model, [i in 1:no_nodes],
                0.9^2 <= real(V[i])^2 + imag(V[i])^2 <= 1.1^2)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0);

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0);

    @constraint(model, S_G - S_Demand .== V .* conj(Ybus * V))

    #---------------------------------------------------

    V_real = real(V)
    V_imag = imag(V)
    
    @objective(
        model,
        Min,
        sum(loss_fun(
            [branch_r[i],
             branch_x[i],
             V_real[branch_fbus[i]],
             V_imag[branch_fbus[i]],
             V_real[branch_tbus[i]],
             V_imag[branch_tbus[i]] ]...)
            for i in 1:length(branch_r) ), );

    #---------------------------------------------------

    optimize!(model)

    # assert_is_solved_and_feasible(model)

    solution_summary(model)

    status = termination_status(model)
    
    #---------------------------------------------------

    objval_solution =
        round(objective_value(model); digits = 2)


    println("Status :" * " $(status)")
    
    println("Objective value (feasible solution) :" *
        " $(objval_solution)")

    return DataFrames.DataFrame(
        ; Bus = all_nodes_idx,
        
        ComplexPowerGen =
            round.(value.(S_G); digits = 2),
        
        VoltageMagnitude =
            round.(abs.(value.(V)); digits = 2),
        
        VoltageAngle_Rad =
            round.(angle.(value.(V));
                   digits = 2),
        
        VoltageAngle_Deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits = 2), )

end


function get_opf_by_scenario(
    case_file;
    active_power_demand_deviation_scale = 1.1,
    reactive_power_demand_deviation_scale = 1.1,
    round_digits = 4,
    no_scenarios = 100 )


    opf_net_optimisation_parameters =
        get_opf_net_optimisation_parameters(
            case_file )

    
    (;gens_nodes_idx,
     all_nodes_idx,
     dyn_pf_fun_kwd_net_idxs,

     gens_installed_capacity,

     gens_Pmax,
     gens_Pmin,
     gens_Qmax,
     gens_Qmin,

     sch_Pg,
     nodes_Pd,
     nodes_Qd,

     gens_cost_coeff_ascen) =
        NamedTupleTools.select(
            opf_net_optimisation_parameters,
            (:gens_nodes_idx,
             :all_nodes_idx,
             :dyn_pf_fun_kwd_net_idxs,

             :gens_installed_capacity,

             :gens_Pmax,
             :gens_Pmin,
             :gens_Qmax,
             :gens_Qmin,

             :sch_Pg,
             :nodes_Pd,
             :nodes_Qd,

             :gens_cost_coeff_ascen))

    #---------------------------------------------------

    (;
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx)  =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))

    #---------------------------------------------------

    base_Pd = copy(nodes_Pd)
    
    base_Qd = copy(nodes_Qd)
    
    #---------------------------------------------------
        
    n2s_gens_nodes =
        OrderedDict(
            a_node => idx
            for (idx, a_node) in
                enumerate(
                    gens_nodes_idx))
    
    n2s_nodes_with_demands =
        OrderedDict(
            a_node => idx
            for (idx, a_node) in
                enumerate(
                    nodes_with_demands_idx)) 

    #---------------------------------------------------

    df_header_sym = [
        [:status,
         :demand_P,
         :demand_Q,
         :objval_solution]...;
        
     [ Symbol("P_gen$(idx)")
       for idx in  gens_nodes_idx]...;
     [ Symbol("Q_gen$(idx)")
       for idx in  gens_nodes_idx]...;
          
     [ Symbol("V_mag_node$(idx)_P")
       for idx in  all_nodes_idx]...;
     [ Symbol("V_ang_rad_node$(idx)_P")
       for idx in  all_nodes_idx]...;        

     [ Symbol("scenario_Pd_node$(idx)")
       for idx in  nodes_with_demands_idx]...;
       [ Symbol("scenario_Qd_node$(idx)")
       for idx in  nodes_with_demands_idx]...]

    opf_df =
        DataFrame(
            OrderedDict( a_header == :status ?
                a_header => [] : 
                a_header => Float64[]
                for a_header in
                    df_header_sym ))
    
    #---------------------------------------------------

    for a_scenario in 1:no_scenarios

        nodes_Pd = get_power_demand_scenario(
            base_Pd;
            scale_factor =
                active_power_demand_deviation_scale)

        nodes_Qd =  get_power_demand_scenario(
            base_Qd;
            scale_factor =
                reactive_power_demand_deviation_scale)

        opf_by_scenario = opf(
            nodes_Pd,
            nodes_Qd;
            opf_net_optimisation_parameters =
                opf_net_optimisation_parameters,
            round_digits =
                round_digits)

        push!(opf_df,
              tuple( [[opf_by_scenario.status,
                       sum(nodes_Pd),
                       sum(nodes_Qd),
                      opf_by_scenario.objval_solution];
                      
                      opf_by_scenario.Gen_P;
                      opf_by_scenario.Gen_Q;
                      
                      opf_by_scenario.V_mag;
                      opf_by_scenario.V_ang_rad;
                      nodes_Pd;
                      nodes_Qd]... ) )
    end

    return (df_header_sym, opf_df)
    
end


function get_opf_by_relaxation_by_scenario(
    case_file;
    active_power_demand_deviation_scale = 1.1,
    reactive_power_demand_deviation_scale = 1.1,
    round_digits = 4,
    no_scenarios = 100 )


    opf_net_optimisation_parameters =
        get_opf_net_optimisation_parameters(
            case_file )

    
    (;gens_nodes_idx,
     all_nodes_idx,
     dyn_pf_fun_kwd_net_idxs,

     gens_installed_capacity,

     gens_Pmax,
     gens_Pmin,
     gens_Qmax,
     gens_Qmin,

     sch_Pg,
     nodes_Pd,
     nodes_Qd,

     gens_cost_coeff_ascen) =
        NamedTupleTools.select(
            opf_net_optimisation_parameters,
            (:gens_nodes_idx,
             :all_nodes_idx,
             :dyn_pf_fun_kwd_net_idxs,

             :gens_installed_capacity,

             :gens_Pmax,
             :gens_Pmin,
             :gens_Qmax,
             :gens_Qmin,

             :sch_Pg,
             :nodes_Pd,
             :nodes_Qd,

             :gens_cost_coeff_ascen))

    #---------------------------------------------------

    (;
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx)  =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))

    #---------------------------------------------------

    base_Pd = copy(nodes_Pd)
    
    base_Qd = copy(nodes_Qd)
    
    #---------------------------------------------------
        
    n2s_gens_nodes =
        OrderedDict(
            a_node => idx
            for (idx, a_node) in
                enumerate(
                    gens_nodes_idx))
    
    n2s_nodes_with_demands =
        OrderedDict(
            a_node => idx
            for (idx, a_node) in
                enumerate(
                    nodes_with_demands_idx)) 

    #---------------------------------------------------

    df_header_sym = [
        [:status,
         :demand_P,
         :demand_Q,
         :objval_solution]...;
        
     [ Symbol("P_gen$(idx)")
       for idx in  gens_nodes_idx]...;
     [ Symbol("Q_gen$(idx)")
       for idx in  gens_nodes_idx]...;
          
     [ Symbol("V_mag_node$(idx)_P")
       for idx in  all_nodes_idx]...;
     [ Symbol("V_ang_rad_node$(idx)_P")
       for idx in  all_nodes_idx]...;        

     [ Symbol("scenario_Pd_node$(idx)")
       for idx in  nodes_with_demands_idx]...;
       [ Symbol("scenario_Qd_node$(idx)")
       for idx in  nodes_with_demands_idx]...]

    opf_df =
        DataFrame(
            OrderedDict( a_header == :status ?
                a_header => [] : 
                a_header => Float64[]
                for a_header in
                    df_header_sym ))
    
    #---------------------------------------------------

    for a_scenario in 1:no_scenarios

        nodes_Pd = get_power_demand_scenario(
            base_Pd;
            scale_factor =
                active_power_demand_deviation_scale)

        nodes_Qd = base_Qd
        
        # nodes_Qd =  get_power_demand_scenario(
        #     base_Qd;
        #     scale_factor =
        #         reactive_power_demand_deviation_scale)

        opf_by_scenario = opf_by_relaxation(
            nodes_Pd,
            nodes_Qd;
            opf_net_optimisation_parameters =
                opf_net_optimisation_parameters,
            round_digits =
                round_digits)

        push!(opf_df,
              tuple( [[opf_by_scenario.status,
                       sum(nodes_Pd),
                       sum(nodes_Qd),
                      opf_by_scenario.objval_solution];
                      
                      opf_by_scenario.Gen_P;
                      opf_by_scenario.Gen_Q;
                      
                      opf_by_scenario.V_mag;
                      opf_by_scenario.V_ang_rad;
                      nodes_Pd;
                      nodes_Qd]... ) )
    end

    return (df_header_sym, opf_df)
    
end


#---------------------------------------------------
#---------------------------------------------------


function solve_economic_dispatch_wt_net_constraint(
    case_file)

    mpc_baseMVA =
        get_matpower_scalar_as_iobuffer_by_case_file(
            case_file; type_key_string = "mpc_baseMVA" )[1]

    #---------------------------------------------------

    mpc_bus = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.bus" )


    mpc_gencost =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.gencost" )


    mpc_gen =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.gen" )


    mpc_branch =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.branch" )

    #---------------------------------------------------

    slack_gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 3 ]


    non_slack_gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 2 ]

    gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 3 || node_type == 2  ]


    non_gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 1 ]


    gens_nodes_with_loc_loads_idx =
        [a_node for (a_node, node_type, node_Pd) in
             zip( mpc_bus.bus_i, mpc_bus.type, mpc_bus.Pd)
             if (node_type == 3 || node_type == 2) &&
                 (node_Pd != 0.0 || node_Pd != 0) ]


    nodes_with_demands_idx =
        [a_node for (a_node, node_Pd, node_Qd) in
             zip( mpc_bus.bus_i, mpc_bus.Pd, mpc_bus.Qd)
             if  (node_Pd != 0.0 || node_Pd != 0) ||
                 (node_Qd != 0.0 || node_Qd != 0)]

    all_nodes_idx = copy(mpc_bus.bus_i)

    #---------------------------------------------------

    no_nodes = length(mpc_bus.bus_i)
    
    #---------------------------------------------------

    dyn_pf_fun_kwd_net_idxs =
        (;
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx,

         nodes_with_demands_idx ) 

    #---------------------------------------------------

    no_gens = length(gens_nodes_idx )

    #---------------------------------------------------
    
    gens_cost_coeff_ascen =
        [ (coeff_0, coeff_1, coeff_2)
          for (coeff_0, coeff_1, coeff_2) in
              zip(mpc_gencost.c0,
                  mpc_gencost.var"...",
                  mpc_gencost.var"c(n-1)" )]


    gens_cost_coeff_decen =
        [ ( coeff_2, coeff_1,  coeff_0)
          for (coeff_2, coeff_1,  coeff_0) in
              zip( mpc_gencost.var"c(n-1)",
                   mpc_gencost.var"...",
                   mpc_gencost.c0
                   )]

    wind_gens_cost_scale  = 0.5
    
    solar_gens_cost_scale = 0.5
    
    wind_gens_cost_coeff =
        [wind_gens_cost_scale .* a_gen_cost_coeff_ascen
         for a_gen_cost_coeff_ascen in
             gens_cost_coeff_ascen] 
            
    solar_gens_cost_coeff  =
        [solar_gens_cost_scale .* a_gen_cost_coeff_ascen
         for a_gen_cost_coeff_ascen in
             gens_cost_coeff_ascen]
    
    #---------------------------------------------------
    
    # gens_Pmax = mpc_gen.Pmax[ gens_nodes_idx ]

    # gens_Pmin = mpc_gen.Pmin[ gens_nodes_idx ]

    # gens_Qmax = mpc_gen.Qmax[ gens_nodes_idx ]

    # gens_Qmin = mpc_gen.Qmin[ gens_nodes_idx ]

    
    gens_Pmax = mpc_gen.Pmax

    gens_Pmin = mpc_gen.Pmin

    gens_Qmax = mpc_gen.Qmax

    gens_Qmin = mpc_gen.Qmin

    sch_Pg = mpc_gen.Pg

    # check if scheduled power is 0, if it is, set Pmax = 0
    
    gens_Pmax = [(sch_Pg[idx] == 0.0 || sch_Pg[idx] == 0 ) ?
        sch_Pg[idx] : pg for (idx, pg) in
                 enumerate(gens_Pmax)]
    
    #---------------------------------------------------

    wind_gens_capacity_scale  = 0.1
    
    solar_gens_capacity_scale = 0.1
    
    wind_gens_Pmin =
        zeros(length(gens_Pmin))
    
    wind_gens_Pmax =
        get_wind_gens_power_forecast(
            mpc_gen.Pmax;
            scale_factor =
                wind_gens_capacity_scale )

    solar_gens_Pmin =
        zeros(length(gens_Pmin))
    
    solar_gens_Pmax =
        get_solar_gens_power_forecast(
            mpc_gen.Pmax;
            scale_factor =
                solar_gens_capacity_scale )
    
    #---------------------------------------------------

    base_Pd = mpc_bus.Pd[nodes_with_demands_idx]

    nodes_Pd = get_active_power_demand_scenario(
        base_Pd;
        scale_factor = 0.1 )

    nodes_Qd = mpc_bus.Qd[nodes_with_demands_idx]

    #---------------------------------------------------

    # Real generation:
    # lower (`lb`) and upper (`ub`) bounds

    P_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmin, no_nodes)
    P_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmax, no_nodes)

    # Reactive generation:
    # lower (`lb`) and upper (`ub`) bounds

    Q_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmin, no_nodes)

    Q_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmax, no_nodes)

    # Power demand levels
    # (real, reactive, and complex form)

    #---------------------------------------------------

    P_wind_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, wind_gens_Pmin, no_nodes)
    
    P_wind_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, wind_gens_Pmax, no_nodes)
    
    #---------------------------------------------------

    P_solar_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, solar_gens_Pmin, no_nodes)
    
    P_solar_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, solar_gens_Pmax, no_nodes)
    
    #---------------------------------------------------
    
    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Pd , no_nodes)

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Qd , no_nodes)

    S_Demand = P_Demand + im * Q_Demand

    #---------------------------------------------------

    branch_data  =
        DataFrame(Any[ mpc_branch.fbus,
                    mpc_branch.tbus,
                    mpc_branch.r,
                    mpc_branch.x,
                    mpc_branch.b ],
                  [:F_BUS, :T_BUS,
                   :BR_R, :BR_X,
                   :BR_Bc] )

    #---------------------------------------------------

    no_edges = size(branch_data, 1)

    #---------------------------------------------------

    incidence_matrix  =
        SparseArrays.sparse(
            branch_data.F_BUS,
            1:no_edges, 1, no_nodes, no_edges) +

                SparseArrays.sparse(
                    branch_data.T_BUS,
                    1:no_edges, -1, no_nodes, no_edges)

    #---------------------------------------------------
    # Matpower
    #---------------------------------------------------

    ys =
         mpc_baseMVA ./ (mpc_branch.r + im * mpc_branch.x) 
       

    y_c =
        1 / 2 * (im *  mpc_branch.b) *
        mpc_baseMVA

    y_sh =  (mpc_bus.Gs .+ im *  mpc_bus.Bs)/mpc_baseMVA

    inv_τ =
        [ (a_ratio == 0.0 || a_ratio == 0) ? 1.0 : 1/a_ratio
              for a_ratio in
                  mpc_branch.ratio  ]
    
    θ_shift = mpc_branch.angle
    
    Yff = (ys +  y_c) .* (inv_τ).^2

    Ytf = -ys .* (inv_τ ./ exp.(im * θ_shift))

    Yft = -ys .* (inv_τ ./ exp.(-im * θ_shift))

    Ytt = (ys +  y_c)

    Y_sh = SparseArrays.spdiagm(y_sh)

    Cf = SparseArrays.sparse(
        1:no_edges,
        branch_data.F_BUS,
             1,  no_edges, no_nodes) 

    Ct = SparseArrays.sparse(
        1:no_edges,
        branch_data.T_BUS,
             1,  no_edges, no_nodes) 
    
    Yf = SparseArrays.spdiagm( Yff ) * Cf +
        SparseArrays.spdiagm( Yft ) * Ct

    Yt = SparseArrays.spdiagm( Ytf ) * Cf +
        SparseArrays.spdiagm( Ytt ) * Ct

    Ybus = Cf' * Yf + Ct' * Yt + Y_sh
        
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    model = Model(Ipopt.Optimizer)
    
    set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    @operator(model,
              gen_cost_fun,
              4,
              a_gen_cost_fun)
    
    @operator(model,
              wind_gen_cost_fun,
              4,
              a_wind_gen_cost_fun)
    
    @operator(model,
              solar_gen_cost_fun,
              4,
              a_solar_gen_cost_func)
    
    #---------------------------------------------------

    @variable(
        model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )


    @variable(
        model,
        P_wind_Gen_lb[i] <= wind_G[i in 1:no_nodes] <= P_wind_Gen_ub[i],
    )


    @variable(
        model,
        P_solar_Gen_lb[i] <= solar_G[i in 1:no_nodes] <= P_solar_Gen_ub[i],
    )
    
    #---------------------------------------------------

    @variable(model, V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model, [i in 1:no_nodes],
                0.9^2 <= real(V[i])^2 + imag(V[i])^2 <= 1.1^2)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0);

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0);


    # @constraint(model, wind_G + solar_G + real.(S_G) - real.(S_Demand) .== real.(V .* conj(Ybus * V)) )


    # @constraint(model, imag.(S_G) - imag.(S_Demand) .== imag.(V .* conj(Ybus * V)))
    
    @constraint(model, wind_G + solar_G + S_G - S_Demand .== V .* conj(Ybus * V))

    #---------------------------------------------------

    P_G = real(S_G)

    # @objective(
    #     model,
    #     Min,
    #     sum(gen_cost_fun(P_G[i],
    #                      gens_cost_coeff_ascen[i]...)
    #         for i in 1:no_gens) +
    #             sum(wind_gen_cost_fun(wind_G[i],
    #                      wind_gens_cost_coeff[i]...)
    #                 for i in 1:no_gens) +
    #                     sum(solar_gen_cost_fun(solar_G[i],
    #                      solar_gens_cost_coeff[i]...)
    #         for i in 1:no_gens), );



    @objective(
        model,
        Min,
        sum(
            gen_cost_fun(
                P_G[i],
                gens_cost_coeff_ascen[i]...) +
                    wind_gen_cost_fun(
                        wind_G[i],
                        wind_gens_cost_coeff[i]...) +
                            solar_gen_cost_fun(
                                solar_G[i],
                                solar_gens_cost_coeff[i]...)
            for i in 1:no_gens),);
    
    #---------------------------------------------------

    optimize!(model)

    solution_summary(model)

    #---------------------------------------------------
    
    objval_solution =
        round(objective_value(model); digits = 4)

    println("Objective value (feasible solution) :" *
        " $(objval_solution)")
    
    #---------------------------------------------------    

    wind_spill = wind_gens_Pmax -
        value.(wind_G)[gens_nodes_idx]

    solar_spill = solar_gens_Pmax -
        value.(solar_G)[gens_nodes_idx]

    spill_df = DataFrames.DataFrame(
        ; Bus = gens_nodes_idx,
        wind_spill =
            round.(wind_spill;
                   digits = 4),
        
        solar_spill =
            round.(solar_spill;
                   digits = 4), )

    dispatch_scenario_df = DataFrames.DataFrame(
        ; Bus = gens_nodes_idx,
                
        gens_P =
            round.(real.(value.(S_G))[gens_nodes_idx];
                   digits = 4),
        wind_P =
            round.(value.(wind_G)[gens_nodes_idx];
                   digits = 4),
        
        solar_P =
            round.(value.(solar_G)[gens_nodes_idx];
                   digits = 4), )

    ed_opf_df = DataFrames.DataFrame(
        ; Bus = 1:no_nodes,
        ComplexPowerGen =
            round.(value.(S_G);
                   digits = 4),        
        wind_Gen =
            round.(value.(wind_G);
                   digits = 4),        
        solar_Gen =
            round.(value.(solar_G); digits = 4),
        
        VoltageMagnitude =
            round.(abs.(value.(V)); digits = 4),

        VoltageAngle_Rad =
            round.(angle.(value.(V)); digits = 4), 
        
        VoltageAngle_Deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits = 4), )
    
    demands_scenario_df = DataFrames.DataFrame(
        ; Bus = nodes_with_demands_idx,
        base_Pd =
            round.(base_Pd; digits = 4),
        
        nodes_Pd =
            round.(nodes_Pd; digits = 4),
        
        nodes_Qd =
            round.(nodes_Qd; digits = 4), )

    
    return (;
            spill_df,
            dispatch_scenario_df,
            ed_opf_df,
            demands_scenario_df)

end



function solve_economic_dispatch_by_net_parameters(
    net_optimisation_parameters;
    round_digits = 4)

    (;no_nodes,
     no_gens,
     no_edges,
     dyn_pf_fun_kwd_net_idxs,
     
     P_Gen_lb,
     P_Gen_ub,
     Q_Gen_lb,
     Q_Gen_ub,

     wind_gens_Pmax,
     solar_gens_Pmax,

     P_wind_Gen_lb,
     P_wind_Gen_ub,
     P_solar_Gen_lb,
     P_solar_Gen_ub,

     base_Pd,
     base_Qd,
     nodes_Pd,
     nodes_Qd,
     
     S_Demand,
     Ybus,

     gens_cost_coeff_ascen,
     wind_gens_cost_coeff,
     solar_gens_cost_coeff) =
        NamedTupleTools.select(
            net_optimisation_parameters,
            (:no_nodes,
             :no_gens,
             :no_edges,             
             :dyn_pf_fun_kwd_net_idxs,
             
             :P_Gen_lb,
             :P_Gen_ub,
             :Q_Gen_lb,
             :Q_Gen_ub,

             :wind_gens_Pmax,
             :solar_gens_Pmax,

             :P_wind_Gen_lb,
             :P_wind_Gen_ub,
             :P_solar_Gen_lb,
             :P_solar_Gen_ub,

             :base_Pd,
             :base_Qd,
             :nodes_Pd,
             :nodes_Qd,
             
             :S_Demand,
             :Ybus,

             :gens_cost_coeff_ascen,
             :wind_gens_cost_coeff,
             :solar_gens_cost_coeff))

    (;
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx)  =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))
    
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    model = Model(Ipopt.Optimizer)
    
    set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    @operator(model,
              gen_cost_fun,
              4,
              a_gen_cost_fun)
    
    @operator(model,
              wind_gen_cost_fun,
              4,
              a_wind_gen_cost_fun)
    
    @operator(model,
              solar_gen_cost_fun,
              4,
              a_solar_gen_cost_func)
    
    #---------------------------------------------------

    @variable(
        model,
        S_G[i in 1:no_nodes] in ComplexPlane(),
        lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
        upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
    )


    @variable(
        model,
        P_wind_Gen_lb[i] <= wind_G[i in 1:no_nodes] <= P_wind_Gen_ub[i],
    )


    @variable(
        model,
        P_solar_Gen_lb[i] <= solar_G[i in 1:no_nodes] <= P_solar_Gen_ub[i],
    )
    
    #---------------------------------------------------

    @variable(model, V[1:no_nodes] in ComplexPlane(),
              start = 1.0 + 0.0im)

    #---------------------------------------------------

    @constraint(model, [i in 1:no_nodes],
                0.9^2 <= real(V[i])^2 + imag(V[i])^2 <= 1.1^2)

    @constraint(model,
                imag(V[slack_gens_nodes_idx[1]]) == 0);

    @constraint(model,
                real(V[slack_gens_nodes_idx[1]]) >= 0);

    @constraint(model, wind_G + solar_G + real.(S_G) - real.(S_Demand) .== real.(V .* conj(Ybus * V)))


    @constraint(model, imag.(S_G) - imag.(S_Demand) .== imag.(V .* conj(Ybus * V)))

    # @constraint(model, wind_G + solar_G + S_G - S_Demand .== V .* conj(Ybus * V))
    
    #---------------------------------------------------

    P_G = real(S_G)

    @objective(
        model,
        Min,
        sum(
            gen_cost_fun(
                P_G[i],
                gens_cost_coeff_ascen[i]...) +
                    wind_gen_cost_fun(
                        wind_G[i],
                        wind_gens_cost_coeff[i]...) +
                            solar_gen_cost_fun(
                                solar_G[i],
                                solar_gens_cost_coeff[i]...)
            for i in 1:no_gens),);
    
    #---------------------------------------------------

    optimize!(model)

    solution_summary(model)

    #---------------------------------------------------
        
    println("Objective value (feasible solution) :" *
        " $(round(objective_value(model);
              digits = round_digits))")
    
    #---------------------------------------------------    

    return (; objval_solution =
        round(objective_value(model);
              digits = round_digits),
     
     wind_spill = round.(wind_gens_Pmax -
        value.(wind_G)[gens_nodes_idx];
                   digits = round_digits),

    solar_spill = round.(solar_gens_Pmax -
        value.(solar_G)[gens_nodes_idx];
                         digits = round_digits),
    gens_P =
            round.(real.(value.(S_G))[gens_nodes_idx];
                   digits = round_digits),
    wind_P =
            round.(value.(wind_G)[gens_nodes_idx];
                   digits = round_digits),
    solar_P =
            round.(value.(solar_G)[gens_nodes_idx];
                   digits = round_digits),
     
    ComplexPowerGen =
            round.(value.(S_G);
                   digits = round_digits),
    wind_Gen =
            round.(value.(wind_G);
                   digits = round_digits),
    solar_Gen =
        round.(value.(solar_G);
               digits = 4),
        
     VoltageMagnitude =
            round.(abs.(value.(V));
                   digits = round_digits),

    VoltageAngle_Rad =
        round.(angle.(value.(V));
               digits = round_digits), 
        
    VoltageAngle_Deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits = round_digits),    
    base_Pd =
        round.(base_Pd;
               digits = round_digits),
    
    base_Qd =
        round.(base_Qd;
               digits = round_digits),

    nodes_Pd =
        round.(nodes_Pd;
               digits = round_digits),
        
            nodes_Qd =
            round.(nodes_Qd;
                   digits = round_digits))

end


function solve_economic_dispatch_by_parameters(
    net_optimisation_parameters;
    round_digits = 4)

    (;
     P_Demand,
     Q_Demand,
     S_Demand,

     wind_gens_Pmax,
     solar_gens_Pmax,
     
     P_wind_Gen_ub,
     P_solar_Gen_ub,
     
     nodes_Pd,
     nodes_Qd,
     
     P_solar_Gen_lb,
     P_wind_Gen_lb,

     no_nodes,
     no_gens,
     no_edges,
     dyn_pf_fun_kwd_net_idxs,
     
     P_Gen_lb,
     P_Gen_ub,
     Q_Gen_lb,
     Q_Gen_ub,

     base_Pd,
     base_Qd,
     
     Ybus,

     gens_cost_coeff_ascen,
     wind_gens_cost_coeff,
     solar_gens_cost_coeff) =
        NamedTupleTools.select(
            net_optimisation_parameters,
            (:P_Demand,
             :Q_Demand,
             :S_Demand,

             :wind_gens_Pmax,
             :solar_gens_Pmax,

             :P_wind_Gen_ub,
             :P_solar_Gen_ub,

             :nodes_Pd,
             :nodes_Qd,

             :P_solar_Gen_lb,
             :P_wind_Gen_lb,

             :no_nodes,
             :no_gens,
             :no_edges,
             :dyn_pf_fun_kwd_net_idxs,

             :P_Gen_lb,
             :P_Gen_ub,
             :Q_Gen_lb,
             :Q_Gen_ub,

             :base_Pd,
             :base_Qd,

             :Ybus,

             :gens_cost_coeff_ascen,
             :wind_gens_cost_coeff,
             :solar_gens_cost_coeff))

    (;
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx)  =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))
    
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    model = Model(Ipopt.Optimizer)
    # model = Model(HiGHS.Optimizer)
    
    set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    @operator(model,
              gen_cost_fun,
              4,
              a_gen_cost_fun)
    
    @operator(model,
              wind_gen_cost_fun,
              4,
              a_wind_gen_cost_fun)
    
    @operator(model,
              solar_gen_cost_fun,
              4,
              a_solar_gen_cost_func)
    
    #---------------------------------------------------

    @variable(
        model,
        P_Gen_lb[i] <= P_G[i in 1:no_nodes] <= P_Gen_ub[i],
    )

    @variable(
        model,
        P_wind_Gen_lb[i] <= wind_G[i in 1:no_nodes] <= P_wind_Gen_ub[i],
    )


    @variable(
        model,
        P_solar_Gen_lb[i] <= solar_G[i in 1:no_nodes] <= P_solar_Gen_ub[i],
    )

    #---------------------------------------------------

    @constraint(
        model,
        sum(wind_G[i] + solar_G[i] + P_G[i] for i in 1:no_gens) == sum(P_Demand[i] for i in 1:no_gens ))


    # @constraint(
    #     model,
    #     sum(wind_G[i] for i in 1:no_gens) + sum(solar_G[i] for i in 1:no_gens) + sum( P_G[i] for i in 1:no_gens) == sum(P_Demand[i] for i in 1:no_gens ))
        
    #---------------------------------------------------
    
    @objective(
        model,
        Min,
        sum(
            gen_cost_fun(
                P_G[i],
                gens_cost_coeff_ascen[i]...) +
                    wind_gen_cost_fun(
                        wind_G[i],
                        wind_gens_cost_coeff[i]...) +
                            solar_gen_cost_fun(
                                solar_G[i],
                                solar_gens_cost_coeff[i]...)
            for i in 1:no_gens),);
    
    #---------------------------------------------------

    optimize!(model)

    solution_summary(model)

    #---------------------------------------------------
        
    println("Objective value (feasible solution) :" *
        " $(round(objective_value(model);
              digits = round_digits))")
    
    #---------------------------------------------------    

    return (; objval_solution =
        round(objective_value(model);
              digits = round_digits),
     
     wind_spill = round.(wind_gens_Pmax -
        value.(wind_G)[gens_nodes_idx];
                   digits = round_digits),

    solar_spill = round.(solar_gens_Pmax -
        value.(solar_G)[gens_nodes_idx];
                         digits = round_digits),
    gens_P =
            round.(value.(P_G)[gens_nodes_idx];
                   digits = round_digits),
    wind_P =
            round.(value.(wind_G)[gens_nodes_idx];
                   digits = round_digits),
    solar_P =
            round.(value.(solar_G)[gens_nodes_idx];
                   digits = round_digits),

    wind_Gen =
            round.(value.(wind_G);
                   digits = round_digits),
    solar_Gen =
        round.(value.(solar_G);
               digits = 4),
   
    base_Pd =
        round.(base_Pd;
               digits = round_digits),
    
    base_Qd =
        round.(base_Qd;
               digits = round_digits),

    nodes_Pd =
        round.(nodes_Pd;
               digits = round_digits),
        
            nodes_Qd =
            round.(nodes_Qd;
                   digits = round_digits))

end



function solve_economic_dispatch_by_scenario(     
    nodes_Pd,
    nodes_Qd,

    wind_gens_Pmax,
    solar_gens_Pmax;
    net_optimisation_parameters =
        net_optimisation_parameters,
    round_digits = 4,
    
    n2s_gens_nodes =
        n2s_gens_nodes,
    n2s_nodes_with_demands =
        n2s_nodes_with_demands,
    n2s_wind_nodes =
        n2s_wind_nodes,
    n2s_solar_nodes =
        n2s_solar_nodes,
    
    wind_sources_nodes = [],
    solar_sources_nodes = [],
    
    wind_gens_cost_coeff =
        wind_gens_cost_coeff,
    solar_gens_cost_coeff =
        solar_gens_cost_coeff,
    
    by_non_thermal_nodes =
        false )

    (;     
     # nodes_Pd,
     # nodes_Qd,
     
     # P_Demand,
     # Q_Demand,
     # S_Demand,

     # wind_gens_Pmax,
     # solar_gens_Pmax,
     
     # P_wind_Gen_ub,
     # P_solar_Gen_ub,
     
     # P_solar_Gen_lb,
     # P_wind_Gen_lb,

     no_nodes,
     no_gens,
     no_edges,
     dyn_pf_fun_kwd_net_idxs,
     
     P_Gen_lb,
     P_Gen_ub,
     Q_Gen_lb,
     Q_Gen_ub,

     base_Pd,
     base_Qd,
     
     Ybus,

     # wind_gens_cost_coeff,
     # solar_gens_cost_coeff,
     
     gens_cost_coeff_ascen) =
        NamedTupleTools.select(
            net_optimisation_parameters,
            (
             # :nodes_Pd,
             # :nodes_Qd,

             # :P_Demand,
             # :Q_Demand,
             # :S_Demand,

             # :wind_gens_Pmax,
             # :solar_gens_Pmax,

             # :P_wind_Gen_ub,
             # :P_solar_Gen_ub,

             # :P_solar_Gen_lb,
             # :P_wind_Gen_lb,

             :no_nodes,
             :no_gens,
             :no_edges,
             :dyn_pf_fun_kwd_net_idxs,

             :P_Gen_lb,
             :P_Gen_ub,
             :Q_Gen_lb,
             :Q_Gen_ub,

             :base_Pd,
             :base_Qd,

             :Ybus,
                
             # :wind_gens_cost_coeff,
             # :solar_gens_cost_coeff,
             :gens_cost_coeff_ascen ))

    (;
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx)  =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))


    #---------------------------------------------------

    # n2s_gens_nodes =
    #         OrderedDict( a_node => idx
    #                      for (idx, a_node) in
    #                          enumerate(gens_nodes_idx))
    
    # n2s_nodes_with_demands =
    #         OrderedDict( a_node => idx
    #                      for (idx, a_node) in
    #                          enumerate(nodes_with_demands_idx)) 

    # if by_non_thermal_nodes == true

    #     n2s_wind_nodes =
    #         OrderedDict( a_node => idx
    #                      for (idx, a_node) in
    #                          enumerate( wind_sources_nodes )) 
    #     n2s_solar_nodes =
    #         OrderedDict( a_node => idx
    #                      for (idx, a_node) in
    #                          enumerate( solar_sources_nodes ))
    # end
    
    #---------------------------------------------------
    
    if by_non_thermal_nodes == false

        wind_gens_Pmin =
            zeros(length(wind_gens_Pmax))

        solar_gens_Pmin =
            zeros(length(solar_gens_Pmax))

        P_wind_Gen_lb =
            SparseArrays.sparsevec(
                gens_nodes_idx, wind_gens_Pmin, no_nodes)

        P_solar_Gen_lb =
            SparseArrays.sparsevec(
                gens_nodes_idx, solar_gens_Pmin, no_nodes)

        P_wind_Gen_ub =
            SparseArrays.sparsevec(
                gens_nodes_idx, wind_gens_Pmax, no_nodes)

        P_solar_Gen_ub =
            SparseArrays.sparsevec(
                gens_nodes_idx, solar_gens_Pmax, no_nodes)
        
    else
        
        wind_gens_Pmin =
            zeros(length( wind_sources_nodes ))

        solar_gens_Pmin =
            zeros(length( solar_sources_nodes ))

        P_wind_Gen_lb =
            SparseArrays.sparsevec(
                 wind_sources_nodes, wind_gens_Pmin, no_nodes)

        P_solar_Gen_lb =
            SparseArrays.sparsevec(
                solar_sources_nodes, solar_gens_Pmin, no_nodes)

        P_wind_Gen_ub =
            SparseArrays.sparsevec(
                 wind_sources_nodes, wind_gens_Pmax, no_nodes)

        P_solar_Gen_ub =
            SparseArrays.sparsevec(
                solar_sources_nodes, solar_gens_Pmax, no_nodes)
        
    end

    #---------------------------------------------------
    
    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Pd , no_nodes)

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            # nodes_Qd,
            base_Qd,
            no_nodes)

    S_Demand = P_Demand + im * Q_Demand
    
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    model = Model(Ipopt.Optimizer)
    # model = Model(HiGHS.Optimizer)
    
    set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    @operator(model,
              gen_cost_fun,
              4,
              a_gen_cost_fun)
    
    @operator(model,
              wind_gen_cost_fun,
              4,
              a_wind_gen_cost_fun)
    
    @operator(model,
              solar_gen_cost_fun,
              4,
              a_solar_gen_cost_func)
    
    #---------------------------------------------------

    @variable(
        model,
        P_Gen_lb[i] <= P_G[i in 1:no_nodes] <= P_Gen_ub[i],
    )

    @variable(
        model,
        P_wind_Gen_lb[i] <= wind_G[i in 1:no_nodes] <= P_wind_Gen_ub[i],
    )


    @variable(
        model,
        P_solar_Gen_lb[i] <= solar_G[i in 1:no_nodes] <= P_solar_Gen_ub[i],
    )

    #---------------------------------------------------

    if by_non_thermal_nodes == false

        @constraint(
            model,
            sum(wind_G[i] + solar_G[i] + P_G[i]
                for i in gens_nodes_idx) ==
                    sum(P_Demand[i]
                        for i in
                            nodes_with_demands_idx ))

        @objective(
            model,
            Min,
            sum(
                gen_cost_fun(
                    P_G[i],
                    gens_cost_coeff_ascen[
                        n2s_gens_nodes[i]]...) +
                        wind_gen_cost_fun(
                            wind_G[i],
                            wind_gens_cost_coeff[
                                n2s_gens_nodes[i]]...) +
                                solar_gen_cost_fun(
                                    solar_G[i],
                                    solar_gens_cost_coeff[
                                        n2s_gens_nodes[i]]...)
                for i in gens_nodes_idx ),);        
        
    else
        
        @constraint(
            model,
            sum(wind_G[i] for i in wind_sources_nodes) +
                sum(solar_G[i] for i in solar_sources_nodes) +
                sum( P_G[i] for i in gens_nodes_idx) ==
                sum(P_Demand[i]
                    for i in nodes_with_demands_idx ))


        @objective(
            model,
            Min,
            sum(
                gen_cost_fun(
                    P_G[i],
                    gens_cost_coeff_ascen[n2s_gens_nodes[i]]...)
                for i in gens_nodes_idx ) +
                    sum(
                        wind_gen_cost_fun(
                            wind_G[i],
                            wind_gens_cost_coeff[
                                n2s_wind_nodes[i]]...)
                        for i in wind_sources_nodes ) +
                            sum(
                                solar_gen_cost_fun(
                                    solar_G[i],
                                    solar_gens_cost_coeff[
                                        n2s_solar_nodes[i]]...)
                for i in solar_sources_nodes),);
        
    end
    
    #---------------------------------------------------    
    #---------------------------------------------------

    optimize!(model)

    solution_summary(model)

    #---------------------------------------------------
        
    println("Objective value (feasible solution) :" *
        " $(round(objective_value(model);
              digits = round_digits))")
    
    #---------------------------------------------------    
    #---------------------------------------------------    
    
    return (;demands = 
        sum(round.(nodes_Pd; digits = round_digits)),
            
            total_cost =
        round(objective_value(model);
              digits = round_digits),

            thermal_gens_cost =  round.(
                [a_gen_cost_fun(
                    pg, gens_cost_coeff_ascen[idx]...)
                 for (idx, pg) in
                     enumerate(value.(P_G)[gens_nodes_idx])];
                digits = round_digits),

            wind_gens_cost =  round.(
                [a_wind_gen_cost_fun(
                    wg, wind_gens_cost_coeff[idx]...)
                 for (idx, wg) in
                     enumerate(value.(wind_G)[
                         by_non_thermal_nodes == false ?
                             gens_nodes_idx :
                             wind_sources_nodes ])];
                digits = round_digits),

            solar_gens_cost = round.(
                [a_solar_gen_cost_func(
                    sg, solar_gens_cost_coeff[idx]...)
                 for (idx, sg) in
                     enumerate(value.(solar_G)[
                         by_non_thermal_nodes == false ?
                             gens_nodes_idx :
                             solar_sources_nodes ]) ];
                digits = round_digits),
                 
            wind_spill = round.(wind_gens_Pmax -
                value.(wind_G)[ by_non_thermal_nodes == false ?
                gens_nodes_idx : wind_sources_nodes ];
                   digits = round_digits),

            solar_spill = round.(solar_gens_Pmax -
                value.(solar_G)[ by_non_thermal_nodes == false ?
                gens_nodes_idx : solar_sources_nodes ];
                         digits = round_digits),
            thermal_P =
                    round.(value.(P_G)[gens_nodes_idx];
                           digits = round_digits),
            wind_P =
                round.(value.(wind_G)[
                    by_non_thermal_nodes == false ?
                        gens_nodes_idx :
                        wind_sources_nodes ];
                           digits = round_digits),
            solar_P =
                round.(value.(solar_G)[
                    by_non_thermal_nodes == false ?
                        gens_nodes_idx :
                        solar_sources_nodes];
                           digits = round_digits),

            thermal_Gen =
                    round.(value.(P_G);
                           digits = round_digits),
            
            wind_Gen =
                    round.(value.(wind_G);
                           digits = round_digits),
            solar_Gen =
                round.(value.(solar_G);
                       digits = 4),

            base_Pd =
                round.(base_Pd;
                       digits = round_digits),

            nodes_Pd =
                round.(nodes_Pd;
                       digits = round_digits) )

end



function get_economic_dispatch_by_scenario(
    case_file,

    wind_gens_cost_scale,    
    solar_gens_cost_scale,

    wind_gens_capacity_scale,    
    solar_gens_capacity_scale,

    active_power_demand_deviation_scale,
    reactive_power_demand_deviation_scale;
    no_scenarios = 100 )


    net_optimisation_parameters =
        get_renewable_energy_net_optimisation_parameters(
            case_file,

            wind_gens_cost_scale,    
            solar_gens_cost_scale,

            wind_gens_capacity_scale,    
            solar_gens_capacity_scale,

            active_power_demand_deviation_scale,
            reactive_power_demand_deviation_scale )


    (nodes_Pd,
     nodes_Qd,

     wind_gens_Pmax,
     solar_gens_Pmax,

     base_Pd,
     base_Qd,

     wind_gens_cost_coeff,
     solar_gens_cost_coeff,
     gens_cost_coeff_ascen,

     gens_installed_capacity,

     dyn_pf_fun_kwd_net_idxs) =
         NamedTupleTools.select(
             net_optimisation_parameters,
             (:nodes_Pd,
              :nodes_Qd,

              :wind_gens_Pmax,
              :solar_gens_Pmax,

              :base_Pd,
              :base_Qd,

              :wind_gens_cost_coeff,
              :solar_gens_cost_coeff,              
              :gens_cost_coeff_ascen,

              :gens_installed_capacity,

              :dyn_pf_fun_kwd_net_idxs))

    #---------------------------------------------------

    (;
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx)  =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))

    #---------------------------------------------------
    
    by_non_thermal_nodes = false

    #---------------------------------------------------

    
    n2s_gens_nodes =
            OrderedDict( a_node => idx
                         for (idx, a_node) in
                             enumerate(
                                 gens_nodes_idx))
    
    n2s_nodes_with_demands =
            OrderedDict(a_node => idx
                         for (idx, a_node) in
                             enumerate(
                                 nodes_with_demands_idx)) 


    if by_non_thermal_nodes == true

        n2s_wind_nodes =
            OrderedDict( a_node => idx
                         for (idx, a_node) in
                             enumerate(
                                 wind_sources_nodes )) 
        n2s_solar_nodes =
            OrderedDict(a_node => idx
                         for (idx, a_node) in
                             enumerate(
                                 solar_sources_nodes ))
    else

        wind_sources_nodes = gens_nodes_idx

        solar_sources_nodes = gens_nodes_idx
        
        n2s_wind_nodes = n2s_gens_nodes
        
        n2s_solar_nodes = n2s_gens_nodes
        
    end
    
    #--------------------------------------------------- 

    df_header_sym = [[:demand, :total_cost]...;
     [ Symbol("thermal_cost_gen$(a_t_gen)")
       for a_t_gen in  gens_nodes_idx]...;
     [ Symbol("wind_cost_gen$(a_t_gen)")
       for a_t_gen in  wind_sources_nodes]...;
     [ Symbol("solar_cost_gen$(a_t_gen)")
       for a_t_gen in  solar_sources_nodes]...;
     
     [ Symbol("wind_spill$(a_t_gen)")
       for a_t_gen in  wind_sources_nodes]...;
     [ Symbol("solar_spill$(a_t_gen)")
       for a_t_gen in  solar_sources_nodes]...;
     
     [ Symbol("thermal_P_gen$(a_t_gen)")
       for a_t_gen in  gens_nodes_idx]...;
     [ Symbol("wind_P_gen$(a_t_gen)")
       for a_t_gen in  wind_sources_nodes]...;
     [ Symbol("solar_P_gen$(a_t_gen)")
       for a_t_gen in  solar_sources_nodes]...;

     [ Symbol("base_Pd_node$(idx)")
       for idx in  nodes_with_demands_idx]...;
     [ Symbol("scenario_Pd_node$(idx)")
       for idx in  nodes_with_demands_idx]...]

    economic_dispatch_df =
        DataFrame(
            OrderedDict(
                a_header => Float64[]
                for a_header in
                    df_header_sym ))
    
    #---------------------------------------------------
   
    wind_gens_cost_coeff =
        [wind_gens_cost_scale .* a_gen_cost_coeff_ascen
         for a_gen_cost_coeff_ascen in
             gens_cost_coeff_ascen] 
            
    solar_gens_cost_coeff  =
        [solar_gens_cost_scale .* a_gen_cost_coeff_ascen
         for a_gen_cost_coeff_ascen in
             gens_cost_coeff_ascen]

    
    #---------------------------------------------------

    
    
    for a_scenario in 1:no_scenarios

        nodes_Pd = get_power_demand_scenario(
            base_Pd;
            scale_factor =
                active_power_demand_deviation_scale)

        nodes_Qd =  get_power_demand_scenario(
            base_Qd;
            scale_factor =
                reactive_power_demand_deviation_scale)

        wind_gens_Pmax =
            get_wind_gens_power_forecast(
                gens_installed_capacity;
                scale_factor =
                    wind_gens_capacity_scale )

        solar_gens_Pmax =
            get_solar_gens_power_forecast(
                gens_installed_capacity;
                scale_factor =
                    solar_gens_capacity_scale )    

        ed_by_scenario =
            solve_economic_dispatch_by_scenario(     
                nodes_Pd,
                nodes_Qd,

                wind_gens_Pmax,
                solar_gens_Pmax;

                net_optimisation_parameters =
                    net_optimisation_parameters,
                round_digits = 4,

                n2s_gens_nodes =
                    n2s_gens_nodes,
                n2s_nodes_with_demands =
                    n2s_nodes_with_demands,
                n2s_wind_nodes =
                    n2s_wind_nodes,
                n2s_solar_nodes =
                    n2s_solar_nodes,

                wind_sources_nodes =
                    wind_sources_nodes,
                solar_sources_nodes =
                    solar_sources_nodes,

                wind_gens_cost_coeff =
                    wind_gens_cost_coeff,
                solar_gens_cost_coeff =
                    solar_gens_cost_coeff,

                by_non_thermal_nodes =
                    by_non_thermal_nodes )

        push!(economic_dispatch_df,
              tuple( [[sum(ed_by_scenario.nodes_Pd),
                      ed_by_scenario.total_cost];
                      ed_by_scenario.thermal_gens_cost;
                      ed_by_scenario.wind_gens_cost;
                      ed_by_scenario.solar_gens_cost;
                      ed_by_scenario.wind_spill;
                      ed_by_scenario.solar_spill;
                      ed_by_scenario.thermal_P;
                      ed_by_scenario.wind_P;
                      ed_by_scenario.solar_P;
                      ed_by_scenario.base_Pd;
                      ed_by_scenario.nodes_Pd]... ) )
    end

    return (; df_header_sym, economic_dispatch_df)
    
end


function solve_unit_commitment_by_scenario(     
    nodes_Pd,
    nodes_Qd,

    wind_gens_Pmax,
    solar_gens_Pmax;
    net_optimisation_parameters =
        net_optimisation_parameters,
    round_digits = 4,
    
    n2s_gens_nodes =
        n2s_gens_nodes,
    n2s_nodes_with_demands =
        n2s_nodes_with_demands,
    n2s_wind_nodes =
        n2s_wind_nodes,
    n2s_solar_nodes =
        n2s_solar_nodes,
    
    wind_sources_nodes = [],
    solar_sources_nodes = [],
    
    wind_gens_cost_coeff =
        wind_gens_cost_coeff,
    solar_gens_cost_coeff =
        solar_gens_cost_coeff,
    
    by_non_thermal_nodes =
        false )

    (;     
     # nodes_Pd,
     # nodes_Qd,
     
     # P_Demand,
     # Q_Demand,
     # S_Demand,

     # wind_gens_Pmax,
     # solar_gens_Pmax,
     
     # P_wind_Gen_ub,
     # P_solar_Gen_ub,
     
     # P_solar_Gen_lb,
     # P_wind_Gen_lb,

     no_nodes,
     no_gens,
     no_edges,
     dyn_pf_fun_kwd_net_idxs,
     
     P_Gen_lb,
     P_Gen_ub,
     Q_Gen_lb,
     Q_Gen_ub,

     base_Pd,
     base_Qd,
     
     Ybus,

     # wind_gens_cost_coeff,
     # solar_gens_cost_coeff,
     
     gens_cost_coeff_ascen) =
        NamedTupleTools.select(
            net_optimisation_parameters,
            (
             # :nodes_Pd,
             # :nodes_Qd,

             # :P_Demand,
             # :Q_Demand,
             # :S_Demand,

             # :wind_gens_Pmax,
             # :solar_gens_Pmax,

             # :P_wind_Gen_ub,
             # :P_solar_Gen_ub,

             # :P_solar_Gen_lb,
             # :P_wind_Gen_lb,

             :no_nodes,
             :no_gens,
             :no_edges,
             :dyn_pf_fun_kwd_net_idxs,

             :P_Gen_lb,
             :P_Gen_ub,
             :Q_Gen_lb,
             :Q_Gen_ub,

             :base_Pd,
             :base_Qd,

             :Ybus,
                
             # :wind_gens_cost_coeff,
             # :solar_gens_cost_coeff,
             :gens_cost_coeff_ascen ))

    (;
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx)  =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))


    #---------------------------------------------------

    # n2s_gens_nodes =
    #         OrderedDict( a_node => idx
    #                      for (idx, a_node) in
    #                          enumerate(gens_nodes_idx))
    
    # n2s_nodes_with_demands =
    #         OrderedDict( a_node => idx
    #                      for (idx, a_node) in
    #                          enumerate(nodes_with_demands_idx))

    # if by_non_thermal_nodes == true

    #     n2s_wind_nodes =
    #         OrderedDict( a_node => idx
    #                      for (idx, a_node) in
    #                          enumerate( wind_sources_nodes )) 
    #     n2s_solar_nodes =
    #         OrderedDict( a_node => idx
    #                      for (idx, a_node) in
    #                          enumerate( solar_sources_nodes ))
    # end
    
    #---------------------------------------------------
    
    if by_non_thermal_nodes == false

        wind_gens_Pmin =
            zeros(length(wind_gens_Pmax))

        solar_gens_Pmin =
            zeros(length(solar_gens_Pmax))

        P_wind_Gen_lb =
            SparseArrays.sparsevec(
                gens_nodes_idx, wind_gens_Pmin, no_nodes)

        P_solar_Gen_lb =
            SparseArrays.sparsevec(
                gens_nodes_idx, solar_gens_Pmin, no_nodes)

        P_wind_Gen_ub =
            SparseArrays.sparsevec(
                gens_nodes_idx, wind_gens_Pmax, no_nodes)

        P_solar_Gen_ub =
            SparseArrays.sparsevec(
                gens_nodes_idx, solar_gens_Pmax, no_nodes)
        
    else
        
        wind_gens_Pmin =
            zeros(length( wind_sources_nodes ))

        solar_gens_Pmin =
            zeros(length( solar_sources_nodes ))

        P_wind_Gen_lb =
            SparseArrays.sparsevec(
                 wind_sources_nodes, wind_gens_Pmin, no_nodes)

        P_solar_Gen_lb =
            SparseArrays.sparsevec(
                solar_sources_nodes, solar_gens_Pmin, no_nodes)

        P_wind_Gen_ub =
            SparseArrays.sparsevec(
                 wind_sources_nodes, wind_gens_Pmax, no_nodes)

        P_solar_Gen_ub =
            SparseArrays.sparsevec(
                solar_sources_nodes, solar_gens_Pmax, no_nodes)
        
    end

    #---------------------------------------------------
    
    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Pd , no_nodes)

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx,
            # nodes_Qd,
            base_Qd,
            no_nodes)

    S_Demand = P_Demand + im * Q_Demand
    
    #---------------------------------------------------
    # Optimal power flow
    #---------------------------------------------------

    # https://stackoverflow.com/questions/74212037/cant-solve-my-non-linear-optimization-problem-julia

    # https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers
    
    nl_solver =
        optimizer_with_attributes(
            Ipopt.Optimizer, "print_level"=>0)

    minlp_solver =
        optimizer_with_attributes(
            Juniper.Optimizer, "nl_solver"=>nl_solver)
    
    model = Model(minlp_solver)
    
    # model = Model(Ipopt.Optimizer)
    # model = Model(HiGHS.Optimizer)
    
    set_silent(model)

    #---------------------------------------------------

    # a_gen_cost_fun(pg, c0, c1, c2)
    
    @operator(model,
              gen_cost_fun,
              4,
              a_gen_cost_fun)
    
    @operator(model,
              wind_gen_cost_fun,
              4,
              a_wind_gen_cost_fun)
    
    @operator(model,
              solar_gen_cost_fun,
              4,
              a_solar_gen_cost_func)
    
    #---------------------------------------------------

    @variable(
        model,
        P_Gen_lb[i] <= P_G[i in 1:no_nodes] <= P_Gen_ub[i],
    )

    @variable(
        model,
        P_wind_Gen_lb[i] <= wind_G[i in 1:no_nodes] <= P_wind_Gen_ub[i],
    )


    @variable(
        model,
        P_solar_Gen_lb[i] <= solar_G[i in 1:no_nodes] <= P_solar_Gen_ub[i],
    )

    # @variable(model, u[i in gens_nodes_idx], Bin)

    @variable(model, u[i = 1:no_gens ], Bin)
    
    # @constraint(
    #     model,
    #     [i = 1:no_gens ], P_G[i] <= P_Gen_ub[i] * u[i ])

    # @constraint(
    #     model,
    #     [i = 1:no_gens ], P_G[i] >= P_Gen_lb[i] * u[i ])

    @constraint(
        model,
        [i in gens_nodes_idx ], P_G[i] <= P_Gen_ub[i] * u[
            n2s_gens_nodes[i] ])

    @constraint(
        model,
        [i in gens_nodes_idx ], P_G[i] >= P_Gen_lb[i] * u[
            n2s_gens_nodes[i]])
    
    #---------------------------------------------------

    if by_non_thermal_nodes == false

        @constraint(
            model,
            sum(wind_G[i] + solar_G[i] + P_G[i]
                for i in gens_nodes_idx) ==
                    sum(P_Demand[i]
                        for i in
                            nodes_with_demands_idx ))

        @objective(
            model,
            Min,
            sum(
                gen_cost_fun(
                    P_G[i],
                    gens_cost_coeff_ascen[
                        n2s_gens_nodes[i]]...) +
                        wind_gen_cost_fun(
                            wind_G[i],
                            wind_gens_cost_coeff[
                                n2s_gens_nodes[i]]...) +
                                solar_gen_cost_fun(
                                    solar_G[i],
                                    solar_gens_cost_coeff[
                                        n2s_gens_nodes[i]]...)
                for i in gens_nodes_idx ) +
                    # sum(gens_cost_coeff_ascen[
                    #     n2s_gens_nodes[i]][2] * u[
                    #         n2s_gens_nodes[i]]
                    #     for i in gens_nodes_idx )
                    sum((gens_cost_coeff_ascen[i])[1] * u[i ]
                        for i in 1:no_gens )
            ,);        
        
    else
        
        @constraint(
            model,
            sum(wind_G[i] for i in wind_sources_nodes) +
                sum(solar_G[i] for i in solar_sources_nodes) +
                sum( P_G[i] for i in gens_nodes_idx) ==
                sum(P_Demand[i]
                    for i in nodes_with_demands_idx ))


        @objective(
            model,
            Min,
            sum(
                gen_cost_fun(
                    P_G[i],
                    gens_cost_coeff_ascen[
                        n2s_gens_nodes[i]]...)
                for i in gens_nodes_idx ) +
                    sum(
                        wind_gen_cost_fun(
                            wind_G[i],
                            wind_gens_cost_coeff[
                                n2s_wind_nodes[i]]...)
                        for i in wind_sources_nodes ) +
                            sum(
                                solar_gen_cost_fun(
                                    solar_G[i],
                                    solar_gens_cost_coeff[
                                        n2s_solar_nodes[i]]...)
                for i in solar_sources_nodes) +
                    # sum(gens_cost_coeff_ascen[
                    #     n2s_gens_nodes[i]][1] * u[
                    #         n2s_gens_nodes[i]]
                    #     for i in gens_nodes_idx )
                    sum((gens_cost_coeff_ascen[i])[1] * u[i ]
                        for i in 1:no_gens )
            ,);
        
    end
    
    #---------------------------------------------------    
    #---------------------------------------------------

    optimize!(model)

    status = termination_status(model)

    # if status != OPTIMAL
    #     return (status = status,)
    # end

    # @assert primal_status(model) == FEASIBLE_POINT

    solution_summary(model)

    #---------------------------------------------------
        
    println("Objective value (feasible solution) :" *
        " $(round(objective_value(model);
              digits = round_digits))")
    
    #---------------------------------------------------    
    #---------------------------------------------------    
    
    return (; status = status,
            u = value.(u),            
            demands = 
                sum(round.(nodes_Pd; digits = round_digits)),
            
            total_cost =
        round(objective_value(model);
              digits = round_digits),

            thermal_gens_cost =  round.(
                [a_gen_cost_fun(
                    pg, gens_cost_coeff_ascen[idx]...)
                 for (idx, pg) in
                     enumerate(value.(P_G)[gens_nodes_idx])];
                digits = round_digits),

            wind_gens_cost =  round.(
                [a_wind_gen_cost_fun(
                    wg, wind_gens_cost_coeff[idx]...)
                 for (idx, wg) in
                     enumerate(value.(wind_G)[
                         by_non_thermal_nodes == false ?
                             gens_nodes_idx :
                             wind_sources_nodes ])];
                digits = round_digits),

            solar_gens_cost = round.(
                [a_solar_gen_cost_func(
                    sg, solar_gens_cost_coeff[idx]...)
                 for (idx, sg) in
                     enumerate(value.(solar_G)[
                         by_non_thermal_nodes == false ?
                             gens_nodes_idx :
                             solar_sources_nodes ]) ];
                digits = round_digits),
                 
            wind_spill = round.(wind_gens_Pmax -
                value.(wind_G)[ by_non_thermal_nodes == false ?
                gens_nodes_idx : wind_sources_nodes ];
                   digits = round_digits),

            solar_spill = round.(solar_gens_Pmax -
                value.(solar_G)[ by_non_thermal_nodes == false ?
                gens_nodes_idx : solar_sources_nodes ];
                         digits = round_digits),
            thermal_P =
                    round.(value.(P_G)[gens_nodes_idx];
                           digits = round_digits),
            wind_P =
                round.(value.(wind_G)[
                    by_non_thermal_nodes == false ?
                        gens_nodes_idx :
                        wind_sources_nodes ];
                           digits = round_digits),
            solar_P =
                round.(value.(solar_G)[
                    by_non_thermal_nodes == false ?
                        gens_nodes_idx :
                        solar_sources_nodes];
                           digits = round_digits),

            thermal_Gen =
                    round.(value.(P_G);
                           digits = round_digits),
            
            wind_Gen =
                    round.(value.(wind_G);
                           digits = round_digits),
            solar_Gen =
                round.(value.(solar_G);
                       digits = 4),

            base_Pd =
                round.(base_Pd;
                       digits = round_digits),

            nodes_Pd =
                round.(nodes_Pd;
                       digits = round_digits) )

end



function get_unit_commitment_by_scenario(
    case_file,

    wind_gens_cost_scale,    
    solar_gens_cost_scale,

    wind_gens_capacity_scale,    
    solar_gens_capacity_scale,

    active_power_demand_deviation_scale,
    reactive_power_demand_deviation_scale;
    no_scenarios = 100 )


    net_optimisation_parameters =
        get_renewable_energy_net_optimisation_parameters(
            case_file,

            wind_gens_cost_scale,    
            solar_gens_cost_scale,

            wind_gens_capacity_scale,    
            solar_gens_capacity_scale,

            active_power_demand_deviation_scale,
            reactive_power_demand_deviation_scale )


    (nodes_Pd,
     nodes_Qd,

     wind_gens_Pmax,
     solar_gens_Pmax,

     base_Pd,
     base_Qd,

     wind_gens_cost_coeff,
     solar_gens_cost_coeff,
     gens_cost_coeff_ascen,

     gens_installed_capacity,

     dyn_pf_fun_kwd_net_idxs) =
         NamedTupleTools.select(
             net_optimisation_parameters,
             (:nodes_Pd,
              :nodes_Qd,

              :wind_gens_Pmax,
              :solar_gens_Pmax,

              :base_Pd,
              :base_Qd,

              :wind_gens_cost_coeff,
              :solar_gens_cost_coeff,              
              :gens_cost_coeff_ascen,

              :gens_installed_capacity,

              :dyn_pf_fun_kwd_net_idxs))

    #---------------------------------------------------

    (;
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx)  =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))

    #---------------------------------------------------
    
    by_non_thermal_nodes = false

    #---------------------------------------------------

    
    n2s_gens_nodes =
            OrderedDict( a_node => idx
                         for (idx, a_node) in
                             enumerate(
                                 gens_nodes_idx))
    
    n2s_nodes_with_demands =
            OrderedDict(a_node => idx
                         for (idx, a_node) in
                             enumerate(
                                 nodes_with_demands_idx)) 


    if by_non_thermal_nodes == true

        n2s_wind_nodes =
            OrderedDict( a_node => idx
                         for (idx, a_node) in
                             enumerate(
                                 wind_sources_nodes )) 
        n2s_solar_nodes =
            OrderedDict(a_node => idx
                         for (idx, a_node) in
                             enumerate(
                                 solar_sources_nodes ))
    else

        wind_sources_nodes = gens_nodes_idx

        solar_sources_nodes = gens_nodes_idx
        
        n2s_wind_nodes = n2s_gens_nodes
        
        n2s_solar_nodes = n2s_gens_nodes
        
    end
    
    #--------------------------------------------------- 

    df_header_sym = [[:status, :demand, :total_cost]...;
      [ Symbol("u_gen$(a_t_gen)")
       for a_t_gen in  gens_nodes_idx]...;               
     [ Symbol("thermal_cost_gen$(a_t_gen)")
       for a_t_gen in  gens_nodes_idx]...;
     [ Symbol("wind_cost_gen$(a_t_gen)")
       for a_t_gen in  wind_sources_nodes]...;
     [ Symbol("solar_cost_gen$(a_t_gen)")
       for a_t_gen in  solar_sources_nodes]...;
     
     [ Symbol("wind_spill$(a_t_gen)")
       for a_t_gen in  wind_sources_nodes]...;
     [ Symbol("solar_spill$(a_t_gen)")
       for a_t_gen in  solar_sources_nodes]...;
     
     [ Symbol("thermal_P_gen$(a_t_gen)")
       for a_t_gen in  gens_nodes_idx]...;
     [ Symbol("wind_P_gen$(a_t_gen)")
       for a_t_gen in  wind_sources_nodes]...;
     [ Symbol("solar_P_gen$(a_t_gen)")
       for a_t_gen in  solar_sources_nodes]...;

     [ Symbol("base_Pd_node$(idx)")
       for idx in  nodes_with_demands_idx]...;
     [ Symbol("scenario_Pd_node$(idx)")
       for idx in  nodes_with_demands_idx]...]

    unit_commitment_df =
        DataFrame(
            OrderedDict( a_header == :status ?
                a_header => [] : a_header => Float64[]
                for a_header in
                    df_header_sym ))
    
    #---------------------------------------------------
   
    wind_gens_cost_coeff =
        [wind_gens_cost_scale .* a_gen_cost_coeff_ascen
         for a_gen_cost_coeff_ascen in
             gens_cost_coeff_ascen] 
            
    solar_gens_cost_coeff  =
        [solar_gens_cost_scale .* a_gen_cost_coeff_ascen
         for a_gen_cost_coeff_ascen in
             gens_cost_coeff_ascen]

    
    #---------------------------------------------------
       
    for a_scenario in 1:no_scenarios

        nodes_Pd = get_power_demand_scenario(
            base_Pd;
            scale_factor =
                active_power_demand_deviation_scale)

        nodes_Qd =  get_power_demand_scenario(
            base_Qd;
            scale_factor =
                reactive_power_demand_deviation_scale)

        wind_gens_Pmax =
            get_wind_gens_power_forecast(
                gens_installed_capacity;
                scale_factor =
                    wind_gens_capacity_scale )

        solar_gens_Pmax =
            get_solar_gens_power_forecast(
                gens_installed_capacity;
                scale_factor =
                    solar_gens_capacity_scale )    

        uc_by_scenario =
            solve_unit_commitment_by_scenario(     
                nodes_Pd,
                nodes_Qd,

                wind_gens_Pmax,
                solar_gens_Pmax;

                net_optimisation_parameters =
                    net_optimisation_parameters,
                round_digits = 4,

                n2s_gens_nodes =
                    n2s_gens_nodes,
                n2s_nodes_with_demands =
                    n2s_nodes_with_demands,
                n2s_wind_nodes =
                    n2s_wind_nodes,
                n2s_solar_nodes =
                    n2s_solar_nodes,

                wind_sources_nodes =
                    wind_sources_nodes,
                solar_sources_nodes =
                    solar_sources_nodes,

                wind_gens_cost_coeff =
                    wind_gens_cost_coeff,
                solar_gens_cost_coeff =
                    solar_gens_cost_coeff,

                by_non_thermal_nodes =
                    by_non_thermal_nodes )

        # if :status ∉ propertynames(uc_by_scenario)
        # end

        push!(unit_commitment_df,
              tuple( [[uc_by_scenario.status,
                       uc_by_scenario.demands,
                       uc_by_scenario.total_cost];
                      uc_by_scenario.u;
                      uc_by_scenario.thermal_gens_cost;
                      uc_by_scenario.wind_gens_cost;
                      uc_by_scenario.solar_gens_cost;
                      uc_by_scenario.wind_spill;
                      uc_by_scenario.solar_spill;
                      uc_by_scenario.thermal_P;
                      uc_by_scenario.wind_P;
                      uc_by_scenario.solar_P;
                      uc_by_scenario.base_Pd;
                      uc_by_scenario.nodes_Pd]... ) )
        
    end

    return (; df_header_sym, unit_commitment_df)
    
end


#---------------------------------------------------
# mixed complementarity problems
#---------------------------------------------------


function solve_mixed_complementarity(
    annualized_capital_cost_I,
    operation_cost_per_MWh_C,
    no_hours_in_a_year_τ,
    scenario_probabilities_θ,
    utility_func_coefficients_A,
    utility_func_coefficients_B;
    round_digits = 4)

    # # Annualized capital cost
    
    # annualized_capital_cost_I = 90_000
    
    # # Operation cost per MWh
    
    # operation_cost_per_MWh_C = 60

    # # Hours per year
    # no_hours_in_a_year_τ = 8_760

    # # Scenario probabilities
    # scenario_probabilities_θ =
    #     [0.2, 0.2, 0.2, 0.2, 0.2]

    # # Utility function coefficients
    
    # utility_func_coefficients_A =
    #     [300, 350, 400, 450, 500]

    # # Utility function coefficients
    # utility_func_coefficients_B = 1
    

    model = Model(PATHSolver.Optimizer)

    set_silent(model)

    # Installed capacity
    @variable(model, x >= 0, start = 1)
    
    # Consumption
    @variable(model, Q[ω = 1:5] >= 0, start = 1)
    
    # Production
    @variable(model, Y[ω = 1:5] >= 0, start = 1)

    # Electricity price
    @variable(model, P[ω = 1:5], start = 1)

    # Capital scarcity margin
    @variable(model, μ[ω = 1:5] >= 0, start = 1)  

    # Unit investment cost equals annualized scarcity
    # profit or investment is 0
    
    @constraint(
        model,
        annualized_capital_cost_I -
            no_hours_in_a_year_τ *
            scenario_probabilities_θ' * μ ⟂ x)

    # Difference between price and scarcity margin is
    # equal to operation cost
    
    @constraint(
        model,
        [ω = 1:5], operation_cost_per_MWh_C -
            (P[ω] - μ[ω]) ⟂ Y[ω])

    # Price is equal to consumer's marginal utility
    
    @constraint(
        model,
        [ω = 1:5],
        P[ω] - (utility_func_coefficients_A[ω] -
            utility_func_coefficients_B * Q[ω]) ⟂ Q[ω])

    # Production is equal to consumption
    
    @constraint(
        model,
        [ω = 1:5],
        Y[ω] - Q[ω] ⟂ P[ω])

    # Production does not exceed capacity
    
    @constraint(model, [ω = 1:5], x - Y[ω] ⟂ μ[ω])

    optimize!(model)

    status = termination_status(model)
    
    # assert_is_solved_and_feasible(model)
    
    solution_summary(model)

    return (;
            status = status,
            
            installed_capacity =
                round(value(x);
                      digits = round_digits),
            
            consumption_scenarios  =
                round.(value.(Q);
                       digits = round_digits),
            
            production_scenarios  =
                round.(value.(Y);
                       digits = round_digits),
            
            price_scenarios  =
                round.(value.(P);
                       digits = round_digits),
            
            capital_scarcity_margin_scenarios =
                round.(value.(μ);
                       digits = round_digits)
            )
    
end



function get_mixed_complementarity(
    annualized_capital_cost_I,
    operation_cost_per_MWh_C,
    no_hours_in_a_year_τ,
    number_of_scenarios,
    # scenario_probabilities_θ,
    utility_func_coefficients_A,
    utility_func_coefficients_B;
    round_digits = 4,
    no_of_cases = 100)

    # # Annualized capital cost
    
    # annualized_capital_cost_I = 90_000
    
    # # Operation cost per MWh
    
    # operation_cost_per_MWh_C = 60

    # # Hours per year
    
    # no_hours_in_a_year_τ = 8_760
    
    # scenario_probabilities_θ =
    #     [0.2, 0.2, 0.2, 0.2, 0.2]
    
    # # Utility function coefficients
    
    # utility_func_coefficients_A =
    #     [300, 350, 400, 450, 500]

    # # Utility function coefficients
    
    # utility_func_coefficients_B = 1
    
    # # # Number of Scenario
    
    # number_of_scenarios = 5
    
    # cals Scenario probabilities
    
    # alpha =
    #     1/2 .* ones(number_of_scenarios)

    alpha =
        ones(number_of_scenarios)
    
    #--------------------------------------------------- 

    df_header_sym = [
        [:status,
         :installed_capacity]...;
      [ Symbol("consumption_scenario$(idx)")
        for idx in 1:number_of_scenarios  ]...;
        
     [ Symbol("production_scenario$(idx)")
       for idx in 1:number_of_scenarios ]...;
        
     [ Symbol("price_scenario$(idx)")
       for idx in 1:number_of_scenarios ]...;
        
     [ Symbol("capital_scarcity_margin_scenario$(idx)")
       for idx in 1:number_of_scenarios ]...]

    mixed_complementarity_df =
        DataFrame(
            OrderedDict( a_header == :status ?
                a_header => [] : a_header => Float64[]
                for a_header in
                    df_header_sym ))

    # no_of_cases = 100
    
    #---------------------------------------------------
    
    for a_case in 1:no_of_cases

        scenario_probabilities_θ =
            rand(Dirichlet(alpha))

        mc =
            solve_mixed_complementarity(
                annualized_capital_cost_I,
                operation_cost_per_MWh_C,
                no_hours_in_a_year_τ,
                scenario_probabilities_θ,
                utility_func_coefficients_A,
                utility_func_coefficients_B)

        push!(mixed_complementarity_df,
              tuple( [[mc.status,
                       mc.installed_capacity];
                      mc.consumption_scenarios;
                      mc.production_scenarios;
                      mc.price_scenarios;
                      mc.capital_scarcity_margin_scenarios
                      ]...) )
        
    end

    #---------------------------------------------------
    
    return (; df_header_sym, mixed_complementarity_df)
    
end



function solve_unit_commitment(
    generators::Vector,
    wind,
    scenario)

    model = Model(HiGHS.Optimizer)

    set_silent(model)

    N = length(generators)

    @variable(
        model, 0 <= g[i = 1:N] <= generators[i].max)

    @variable(model, 0 <= w <= scenario.wind)

    @constraint(
        model,
        sum(g[i] for i in 1:N) + w == scenario.demand)

    # !!! New: add binary on-off variables
    # for each generator

    @variable(model, u[i = 1:N], Bin)

    @constraint(
        model,
        [i = 1:N], g[i] <= generators[i].max * u[i])

    @constraint(
        model,
        [i = 1:N], g[i] >= generators[i].min * u[i])

    @objective(
        model,
        Min,
        sum(generators[i].variable_cost * g[i]
            for i in 1:N) + wind.variable_cost * w +
        # !!! new
        sum(generators[i].fixed_cost * u[i] for i in 1:N)
    )

    optimize!(model)

    status = termination_status(model)

    if status != OPTIMAL
        return (status = status,)
    end

    @assert primal_status(model) == FEASIBLE_POINT

    return (
        status = status,
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        u = value.(u),
        total_cost = objective_value(model), )
end


function solve_nonlinear_economic_dispatch(
    generators::Vector,
    wind,
    scenario;
    silent::Bool = false, )

    model = Model(Ipopt.Optimizer)

    if silent
        set_silent(model)
    end

    @operator(
        model,
        op_tcf,
        1,
        thermal_cost_function)

    N = length(generators)

    @variable(
        model,
        generators[i].min <= g[i = 1:N] <=
            generators[i].max)

    @variable(
        model,
        0 <= w <= scenario.wind)

    @objective(
        model,
        Min,
        sum(generators[i].variable_cost * op_tcf(g[i])
            for i in 1:N) + wind.variable_cost * w, )

    @constraint(
        model,
        sum(g[i]
            for i in 1:N) + sqrt(w) == scenario.demand)

    optimize!(model)

    # assert_is_solved_and_feasible(model)
    return (
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        total_cost = objective_value(model), )
end


# function solve_economic_dispatch(
#     generators::Vector,
#     wind,
#     scenario)

#     # Define the economic dispatch (ED) model
#     model = Model(HiGHS.Optimizer)
    
#     set_silent(model)

#     # Define decision variables
#     # power output of generators

#     N = length(generators)
#     @variable(
#         model,
#         generators[i].min <= g[i = 1:N] <= generators[i].max)
#     # wind power injection
#     @variable(model, 0 <= w <= scenario.wind)

#     # Define the objective function
#     @objective(
#         model,
#         Min,
#         sum(generators[i].variable_cost * g[i]
#             for i in 1:N) +
#         wind.variable_cost * w, )

#     # Define the power balance constraint
#     @constraint(
#         model,
#         sum(g[i] for i in 1:N) + w == scenario.demand)

#     # Solve statement
#     optimize!(model)

#     # assert_is_solved_and_feasible(model)

#     # return the optimal value of
#     # the objective function and its minimizers
#     return (
#         g = value.(g),
#         w = value(w),
#         wind_spill = scenario.wind - value(w),
#         total_cost = objective_value(model), )
# end


# function solve_economic_dispatch_inplace(
#     generators::Vector,
#     wind,
#     scenario,
#     scale::AbstractVector{Float64}, )
#     obj_out = Float64[]
#     w_out = Float64[]
#     g1_out = Float64[]
#     g2_out = Float64[]

#     # This function only works for two generators
#     @assert length(generators) == 2

#     model = Model(HiGHS.Optimizer)
#     set_silent(model)
#     N = length(generators)
#     @variable(
#         model,
#         generators[i].min <= g[i = 1:N] <=
#             generators[i].max)

#     @variable(model, 0 <= w <= scenario.wind)

#     @objective(
#         model,
#         Min,
#         sum(
#             generators[i].variable_cost * g[i]
#             for i in 1:N) + wind.variable_cost * w, )

#     @constraint(
#         model,
#         sum(g[i]
#             for i in 1:N) + w == scenario.demand)

#     for c_g1_scale in scale
#         @objective(
#             model,
#             Min,
#             c_g1_scale * generators[1].variable_cost *
#                 g[1] +
#             generators[2].variable_cost * g[2] +
#             wind.variable_cost * w, )

#         optimize!(model)
#         # assert_is_solved_and_feasible(model)
#         push!(obj_out, objective_value(model))
#         push!(w_out, value(w))
#         push!(g1_out, value(g[1]))
#         push!(g2_out, value(g[2]))
#     end

#     df = DataFrames.DataFrame(;
#         scale = scale,
#         dispatch_G1 = g1_out,
#         dispatch_G2 = g2_out,
#         dispatch_wind = w_out,
#         spillage_wind = scenario.wind .- w_out,
#         total_cost = obj_out,
#     )
#     return df
# end

#---------------------------------------------------
# mixed complementarity problems
#---------------------------------------------------

# https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/complementarity/#Electricity-consumption


function mixed_complementarity()

    # Annualized capital cost
    
    I = 90_000
    
    # Operation cost per MWh
    
    C = 60

    # Hours per year
    
    τ = 8_760

    # Scenario probabilities
    θ = [0.2, 0.2, 0.2, 0.2, 0.2]

    # Utility function coefficients
    
    A = [300, 350, 400, 450, 500]

    # Utility function coefficients
    B = 1                          

    model = Model(PATHSolver.Optimizer)

    set_silent(model)

    # Installed capacity
    @variable(model, x >= 0, start = 1)
    
    # Consumption
    @variable(model, Q[ω = 1:5] >= 0, start = 1)
    
    # Production
    @variable(model, Y[ω = 1:5] >= 0, start = 1)

    # Electricity price
    @variable(model, P[ω = 1:5], start = 1)

    # Capital scarcity margin
    @variable(model, μ[ω = 1:5] >= 0, start = 1)  

    # Unit investment cost equals annualized scarcity
    # profit or investment is 0
    @constraint(model, I - τ * θ' * μ ⟂ x)

    # Difference between price and scarcity margin is
    # equal to operation cost
    @constraint(model, [ω = 1:5], C - (P[ω] - μ[ω]) ⟂ Y[ω])

    # Price is equal to consumer's marginal utility
    @constraint(
        model,
        [ω = 1:5],
        P[ω] - (A[ω] - B * Q[ω]) ⟂ Q[ω])

    # Production is equal to consumption
    @constraint(model, [ω = 1:5], Y[ω] - Q[ω] ⟂ P[ω])

    # Production does not exceed capacity
    @constraint(model, [ω = 1:5], x - Y[ω] ⟂ μ[ω])

    optimize!(model)

    # assert_is_solved_and_feasible(model)
    solution_summary(model)
    
end

#---------------------------------------------------
#---------------------------------------------------



function case8_driver_examples_from_jump()

    # mixed complementarity problems

    # Annualized capital cost
    I = 90_000
    
    # Operation cost per MWh
    C = 60

    # Hours per year
    τ = 8_760

    # Scenario probabilities
    θ = [0.2, 0.2, 0.2, 0.2, 0.2]

    # Utility function coefficients
    A = [300, 350, 400, 450, 500]

    # Utility function coefficients
    B = 1                          

    model = Model(PATHSolver.Optimizer)

    set_silent(model)

    # Installed capacity
    @variable(model, x >= 0, start = 1)
    
    # Consumption
    @variable(model, Q[ω = 1:5] >= 0, start = 1)
    
    # Production
    @variable(model, Y[ω = 1:5] >= 0, start = 1)

    # Electricity price
    @variable(model, P[ω = 1:5], start = 1)

    # Capital scarcity margin
    @variable(model, μ[ω = 1:5] >= 0, start = 1)  

    # Unit investment cost equals annualized scarcity
    # profit or investment is 0
    @constraint(model, I - τ * θ' * μ ⟂ x)

    # Difference between price and scarcity margin is
    # equal to operation cost
    @constraint(model, [ω = 1:5], C - (P[ω] - μ[ω]) ⟂ Y[ω])

    # Price is equal to consumer's marginal utility
    @constraint(
        model,
        [ω = 1:5],
        P[ω] - (A[ω] - B * Q[ω]) ⟂ Q[ω])

    # Production is equal to consumption
    @constraint(model, [ω = 1:5], Y[ω] - Q[ω] ⟂ P[ω])

    # Production does not exceed capacity
    @constraint(model, [ω = 1:5], x - Y[ω] ⟂ μ[ω])
    
    optimize!(model)
    
    # assert_is_solved_and_feasible(model)
    solution_summary(model)

end

# op_sol = case8_driver_examples_from_jump()

# ed_opf =
#     solve_economic_dispatch_wt_net_constraint(
#     case_file)

# net_optimisation_parameters =
#     get_renewable_energy_net_optimisation_parameters(
#         case_file,

#         wind_gens_cost_scale,    
#         solar_gens_cost_scale,

#         wind_gens_capacity_scale,    
#         solar_gens_capacity_scale,

#         active_power_demand_deviation_scale,
#         reactive_power_demand_deviation_scale )


# economic_dispatch_by_net_parameters =
#     solve_economic_dispatch_by_net_parameters(
#         net_optimisation_parameters;
#         round_digits = 4)

# economic_dispatch_by_parameters_nothing =
#     solve_economic_dispatch_by_parameters(
#         net_optimisation_parameters;
#         round_digits = 4)



# # getproperty(
# #     getproperty(
# #         net_optimisation_parameters,
# #         :dyn_pf_fun_kwd_net_idxs),
# #     :nodes_with_demands_idx)


# # aa = getproperty(
# #         net_optimisation_parameters,
# #         :gens_Pmax)


# # getproperty(
# #     getproperty(
# #         net_optimisation_parameters,
# #         :dyn_pf_fun_kwd_net_idxs),
# #     :nodes_with_demands_idx)

# (nodes_Pd,
#  nodes_Qd,
 
#  wind_gens_Pmax,
#  solar_gens_Pmax,
 
#  wind_gens_cost_coeff,
#  solar_gens_cost_coeff) =
#      NamedTupleTools.select(
#          net_optimisation_parameters,
#          (:nodes_Pd,
#           :nodes_Qd,

#           :wind_gens_Pmax,
#           :solar_gens_Pmax,
          
#           :wind_gens_cost_coeff,
#           :solar_gens_cost_coeff))

# tt_economic_dispatch_by_scenario =
#     solve_economic_dispatch_by_scenario(     
#         nodes_Pd,
#         nodes_Qd,

#         wind_gens_Pmax,
#         solar_gens_Pmax;
#         net_optimisation_parameters =
#             net_optimisation_parameters,
#         round_digits = 4,
#         wind_gens_cost_coeff =
#             wind_gens_cost_coeff,
#         solar_gens_cost_coeff =
#             solar_gens_cost_coeff,
#         by_non_thermal_nodes =
#         false)

# function get_economic_dispatch_by_scenario(
#     case_file,

#     wind_gens_cost_scale,    
#     solar_gens_cost_scale,

#     wind_gens_capacity_scale,    
#     solar_gens_capacity_scale,

#     active_power_demand_deviation_scale,
#     reactive_power_demand_deviation_scale)


#     net_optimisation_parameters =
#         get_renewable_energy_net_optimisation_parameters(
#             case_file,

#             wind_gens_cost_scale,    
#             solar_gens_cost_scale,

#             wind_gens_capacity_scale,    
#             solar_gens_capacity_scale,

#             active_power_demand_deviation_scale,
#             reactive_power_demand_deviation_scale )


#     (nodes_Pd,
#      nodes_Qd,

#      wind_gens_Pmax,
#      solar_gens_Pmax,

#      base_Pd,
#      base_Qd,

#      wind_gens_cost_coeff,
#      solar_gens_cost_coeff,
#      gens_cost_coeff_ascen,

#      gens_installed_capacity,

#      dyn_pf_fun_kwd_net_idxs) =
#          NamedTupleTools.select(
#              net_optimisation_parameters,
#              (:nodes_Pd,
#               :nodes_Qd,

#               :wind_gens_Pmax,
#               :solar_gens_Pmax,

#               :base_Pd,
#               :base_Qd,

#               :wind_gens_cost_coeff,
#               :solar_gens_cost_coeff,              
#               :gens_cost_coeff_ascen,

#               :gens_installed_capacity,

#               :dyn_pf_fun_kwd_net_idxs))

#     #---------------------------------------------------

#     (;
#      slack_gens_nodes_idx,
#      non_slack_gens_nodes_idx,
#      gens_nodes_idx,
#      non_gens_nodes_idx,
#      gens_nodes_with_loc_loads_idx,
#      all_nodes_idx,
#      nodes_with_demands_idx)  =
#          NamedTupleTools.select(
#              dyn_pf_fun_kwd_net_idxs,
#              (:slack_gens_nodes_idx,
#               :non_slack_gens_nodes_idx,
#               :gens_nodes_idx,
#               :non_gens_nodes_idx,
#               :gens_nodes_with_loc_loads_idx,
#               :all_nodes_idx,
#               :nodes_with_demands_idx))

#     #---------------------------------------------------
    
#     by_non_thermal_nodes = false

#     #---------------------------------------------------

    
#     n2s_gens_nodes =
#             OrderedDict( a_node => idx
#                          for (idx, a_node) in
#                              enumerate(
#                                  gens_nodes_idx))
    
#     n2s_nodes_with_demands =
#             OrderedDict(a_node => idx
#                          for (idx, a_node) in
#                              enumerate(
#                                  nodes_with_demands_idx)) 


#     if by_non_thermal_nodes == true

#         n2s_wind_nodes =
#             OrderedDict( a_node => idx
#                          for (idx, a_node) in
#                              enumerate(
#                                  wind_sources_nodes )) 
#         n2s_solar_nodes =
#             OrderedDict(a_node => idx
#                          for (idx, a_node) in
#                              enumerate(
#                                  solar_sources_nodes ))
#     else

#         wind_sources_nodes = gens_nodes_idx

#         solar_sources_nodes = gens_nodes_idx
        
#         n2s_wind_nodes = n2s_gens_nodes
        
#         n2s_solar_nodes = n2s_gens_nodes
        
#     end
    
#     #--------------------------------------------------- 

#     df_header_sym = [[:demand, :total_cost]...;
#      [ Symbol("thermal_gen$(a_t_gen)_cost")
#        for a_t_gen in  gens_nodes_idx]...;
#      [ Symbol("wind_gen$(a_t_gen)_cost")
#        for a_t_gen in  wind_sources_nodes]...;
#      [ Symbol("solar_gen$(a_t_gen)_cost")
#        for a_t_gen in  solar_sources_nodes]...;
     
#      [ Symbol("wind_spill$(a_t_gen)")
#        for a_t_gen in  wind_sources_nodes]...;
#      [ Symbol("solar_spill$(a_t_gen)")
#        for a_t_gen in  solar_sources_nodes]...;
     
#      [ Symbol("thermal_gen$(a_t_gen)_P")
#        for a_t_gen in  gens_nodes_idx]...;
#      [ Symbol("wind_gen$(a_t_gen)_P")
#        for a_t_gen in  wind_sources_nodes]...;
#      [ Symbol("solar_gen$(a_t_gen)_P")
#        for a_t_gen in  solar_sources_nodes]...;

#      [ Symbol("node$(idx)_base_Pd")
#        for idx in  nodes_with_demands_idx]...;
#      [ Symbol("node$(idx)_scenario_Pd")
#        for idx in  nodes_with_demands_idx]...]

#     economic_dispatch_df =
#         DataFrame(
#             OrderedDict(
#                 a_header => Float64[]
#                 for a_header in
#                     df_header_sym ))
    
#     #---------------------------------------------------
   
#     wind_gens_cost_coeff =
#         [wind_gens_cost_scale .* a_gen_cost_coeff_ascen
#          for a_gen_cost_coeff_ascen in
#              gens_cost_coeff_ascen] 
            
#     solar_gens_cost_coeff  =
#         [solar_gens_cost_scale .* a_gen_cost_coeff_ascen
#          for a_gen_cost_coeff_ascen in
#              gens_cost_coeff_ascen]

    
#     #---------------------------------------------------

#     no_scenarios = 100
    
#     for a_scenario in 1:no_scenarios

#         nodes_Pd = get_power_demand_scenario(
#             base_Pd;
#             scale_factor =
#                 active_power_demand_deviation_scale)

#         nodes_Qd =  get_power_demand_scenario(
#             base_Qd;
#             scale_factor =
#                 reactive_power_demand_deviation_scale)

#         wind_gens_Pmax =
#             get_wind_gens_power_forecast(
#                 gens_installed_capacity;
#                 scale_factor =
#                     wind_gens_capacity_scale )

#         solar_gens_Pmax =
#             get_solar_gens_power_forecast(
#                 gens_installed_capacity;
#                 scale_factor =
#                     solar_gens_capacity_scale )    

#         ed_by_scenario =
#             solve_economic_dispatch_by_scenario(     
#                 nodes_Pd,
#                 nodes_Qd,

#                 wind_gens_Pmax,
#                 solar_gens_Pmax;

#                 net_optimisation_parameters =
#                     net_optimisation_parameters,
#                 round_digits = 4,

#                 n2s_gens_nodes =
#                     n2s_gens_nodes,
#                 n2s_nodes_with_demands =
#                     n2s_nodes_with_demands,
#                 n2s_wind_nodes =
#                     n2s_wind_nodes,
#                 n2s_solar_nodes =
#                     n2s_solar_nodes,

#                 wind_sources_nodes =
#                     wind_sources_nodes,
#                 solar_sources_nodes =
#                     solar_sources_nodes,

#                 wind_gens_cost_coeff =
#                     wind_gens_cost_coeff,
#                 solar_gens_cost_coeff =
#                     solar_gens_cost_coeff,

#                 by_non_thermal_nodes =
#                     by_non_thermal_nodes )

#         push!(economic_dispatch_df,
#               tuple( [[sum(ed_by_scenario.nodes_Pd),
#                       ed_by_scenario.total_cost];
#                       ed_by_scenario.thermal_gens_cost;
#                       ed_by_scenario.wind_gens_cost;
#                       ed_by_scenario.solar_gens_cost;
#                       ed_by_scenario.wind_spill;
#                       ed_by_scenario.solar_spill;
#                       ed_by_scenario.thermal_P;
#                       ed_by_scenario.wind_P;
#                       ed_by_scenario.solar_P;
#                       ed_by_scenario.base_Pd;
#                       ed_by_scenario.nodes_Pd]... ) )
#     end
    

#     header_names_economic_dispatch_df =
#         names(economic_dispatch_df)

#     dispatch_plot = StatsPlots.@df(
#         economic_dispatch_df,
#         Plots.plot(
#             :thermal_gen1_cost,
#             [:thermal_gen1_P],
#             labels = ["G1"],
#             title = "Thermal Dispatch",
#             legend = :bottomright,
#             linewidth = 3,
#             xlabel = "Demand",
#             ylabel = "Dispatch [MW]",
#         ),
#     )

#     wind_and_solar_plot = StatsPlots.@df(
#         economic_dispatch_df,
#         Plots.plot(
#             :demand,
#             [:thermal_gen1_cost,
#              :thermal_gen2_cost,
#              :wind_gen1_cost,
#              :wind_gen2_cost,
#              :solar_gen1_cost,
#              :solar_gen2_cost],
#             labels = ["thermal_gen1_cost"  "thermal_gen2_cost" "wind_gen1_cost"  "wind_gen2_cost" "solar_gen1_cost"  "solar_gen2_cost"],
#             title = "Wind",
#             legend = :bottomright,
#             linewidth = 3,
#             xlabel = "Demand [MW]",
#             ylabel = "Cost [R]",
#         ),
#     )

#     plt1 = Plots.plot(dispatch_plot, wind_and_solar_plot)
    
    
# end
