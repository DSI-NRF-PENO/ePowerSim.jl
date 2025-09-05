# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


# using Pkg

# QuickDynamicTest =
#     joinpath(@__DIR__,"..","..")

# cd(QuickDynamicTest)

# Pkg.activate( QuickDynamicTest )


#---------------------------------------------------
####################################################
#---------------------------------------------------

#---------------------------------------------------

# thermal generators:

function ThermalGenerator(
    min::Float64,
    max::Float64,
    fixed_cost::Float64,
    variable_cost::Float64,
)
    return (
        min = min,
        max = max,
        fixed_cost = fixed_cost,
        variable_cost = variable_cost,
    )
end

# A wind generator
WindGenerator(variable_cost::Float64) =
    (variable_cost = variable_cost,)

# a scenario
function Scenario(demand::Float64, wind::Float64)
    return (demand = demand, wind = wind)
end

scenario = Scenario(1500.0, 200.0)

# economic_dispatch,

# economic_dispatch
function solve_economic_dispatch(
    generators::Vector, wind, scenario)

    # Define the economic dispatch (ED) model

    model = Model(HiGHS.Optimizer)

    set_silent(model)

    # Define decision variables
    # power output of generators

    N = length(generators)

    @variable(model,
              generators[i].min <= g[i = 1:N] <=
                  generators[i].max)

    # wind power injection

    @variable(model, 0 <= w <= scenario.wind)

    # Define the objective function

    @objective(
        model,
        Min,
        sum(generators[i].variable_cost * g[i]
            for i in 1:N) +
                wind.variable_cost * w, )

    # Define the power balance constraint
    @constraint(model,
                sum(g[i]
                    for i in 1:N) + w == scenario.demand )

    # Solve statement

    optimize!(model)

    # assert_is_solved_and_feasible(model)

    # return the optimal value of the objective function and its minimizers

    return (
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        total_cost = objective_value(model), )
end

# Economic dispatch with adjustable incremental costs
function scale_generator_cost(g, scale)

    return ThermalGenerator(
        g.min, g.max, g.fixed_cost,
        scale * g.variable_cost)
end

function scale_demand(scenario, scale)

    return Scenario(
        scale * scenario.demand, scenario.wind)
end

# Modifying the JuMP model in-place
function solve_economic_dispatch_inplace(
    generators::Vector,
    wind,
    scenario,
    scale::AbstractVector{Float64},
)
    obj_out = Float64[]
    w_out   = Float64[]
    g1_out  = Float64[]
    g2_out  = Float64[]

    # This function only works for two generators

    @assert length(generators) == 2

    model = Model(HiGHS.Optimizer)

    set_silent(model)

    N = length(generators)

    @variable(model,
              generators[i].min <= g[i = 1:N] <=
                  generators[i].max)

    @variable(model, 0 <= w <= scenario.wind)

    @objective(
        model,
        Min,
        sum(generators[i].variable_cost * g[i]
            for i in 1:N) + wind.variable_cost * w, )

    @constraint(model, sum(g[i]
                           for i in 1:N) + w ==
                               scenario.demand)
    for c_g1_scale in scale
        @objective(
            model,
            Min,
            c_g1_scale * generators[1].variable_cost * g[1] +
            generators[2].variable_cost * g[2] +
            wind.variable_cost * w, )

        optimize!(model)

        # assert_is_solved_and_feasible(model)

        push!(obj_out, objective_value(model))

        push!(w_out, value(w))

        push!(g1_out, value(g[1]))

        push!(g2_out, value(g[2]))

    end

    df = DataFrames.DataFrame(
        ; scale = scale,
        dispatch_G1 = g1_out,
        dispatch_G2 = g2_out,
        dispatch_wind = w_out,
        spillage_wind = scenario.wind .- w_out,
        total_cost = obj_out, )

    return df
end



# Unit commitment

function solve_unit_commitment(
    generators::Vector, wind, scenario)

    model = Model(HiGHS.Optimizer)

    set_silent(model)

    N = length(generators)

    @variable(model, 0 <= g[i = 1:N] <= generators[i].max)

    @variable(model, 0 <= w <= scenario.wind)

    @constraint(model, sum(g[i] for i in 1:N) + w ==
        scenario.demand)

    # !!! New: add binary on-off variables for each generator
    @variable(model, u[i = 1:N], Bin)

    @constraint(model, [i = 1:N], g[i] <=
        generators[i].max * u[i])

    @constraint(model, [i = 1:N], g[i] >=
        generators[i].min * u[i])

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
        total_cost = objective_value(model),
    )
end

# Nonlinear economic dispatch

"""
    thermal_cost_function(g)

A user-defined thermal cost function in pure-Julia! You can include
nonlinearities, and even things like control flow.

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

function solve_nonlinear_economic_dispatch(
    generators::Vector,
    wind,
    scenario;
    silent::Bool = false, )
    model = Model(Ipopt.Optimizer)
    if silent
        set_silent(model)
    end

    @operator(model, op_tcf, 1, thermal_cost_function)
    N = length(generators)

    @variable(model, generators[i].min <= g[i = 1:N] <=
        generators[i].max)
    @variable(model, 0 <= w <= scenario.wind)

    @objective(
        model,
        Min,
        sum(generators[i].variable_cost * op_tcf(g[i])
            for i in 1:N) +
                wind.variable_cost * w, )

    @constraint(model, sum(g[i]
                           for i in 1:N) + sqrt(w) ==
                               scenario.demand)
    optimize!(model)
    # assert_is_solved_and_feasible(model)
    return (
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        total_cost = objective_value(model),
    )
end

#---------------------------------------------------
#---------------------------------------------------



function a_gen_cost_fun(pg, c0, c1, c2)

    return c0 + c1 * pg + c2 * pg^2
end


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
#---------------------------------------------------

function opf( case_file )

    mpc_baseMVA =
        get_matpower_scalar_as_iobuffer_by_case_file(
            case_file; type_key_string = "mpc_baseMVA" )[1]

    #---------------------------------------------------

    mpc_bus = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.bus" )


    mpc_gencost = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.gencost" )


    mpc_gen = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.gen" )


    mpc_branch = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file; type_key= "mpc.branch" )

    #---------------------------------------------------
    #---------------------------------------------------
    
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_mpc( mpc_bus  )
    
    (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx ) =
         NamedTupleTools.select(
             net_nodes_type_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))

    #---------------------------------------------------

    no_nodes = length( all_nodes_idx )
    
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
    
    gens_Pmax = mpc_gen.Pmax[ gens_nodes_idx ]

    gens_Pmin = mpc_gen.Pmin[ gens_nodes_idx ]

    gens_Qmax = mpc_gen.Qmax[ gens_nodes_idx ]

    gens_Qmin = mpc_gen.Qmin[ gens_nodes_idx ]

    #---------------------------------------------------

    nodes_Pd = mpc_bus.Pd[nodes_with_demands_idx]

    nodes_Qd = mpc_bus.Qd[nodes_with_demands_idx]

    #---------------------------------------------------
    #---------------------------------------------------

    # Real generation: lower (`lb`) and upper (`ub`) bounds

    P_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmin, no_nodes)
    P_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmax, no_nodes)

    # Reactive generation: lower (`lb`) and upper (`ub`) bounds

    Q_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmin, no_nodes)

    Q_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmax, no_nodes)

    # Power demand levels (real, reactive, and complex form)

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

    # @objective(
    #     model,
    #     Min,
    #     (0.11 * P_G[1]^2 + 5 * P_G[1] + 150) +
    #     (0.085 * P_G[2]^2 + 1.2 * P_G[2] + 600) +
    #     (0.1225 * P_G[3]^2 + P_G[3] + 335),
    # );


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

    #---------------------------------------------------

    objval_solution =
        round(objective_value(model); digits = 2)

    println("Objective value (feasible solution) :" *
        " $(objval_solution)")

    return DataFrames.DataFrame(
        ; Bus = 1:no_nodes,
        ComplexPowerGen =
            round.(value.(S_G); digits = 2),
        
        VoltageMagnitude =
            round.(abs.(value.(V)); digits = 2),
        
        VoltageAngle_Deg =
            round.(rad2deg.(angle.(value.(V)));
                   digits = 2), )

end


function opf_by_relaxation( case_file )

    # case_file = case_9_file
    
    mpc_baseMVA =
        get_matpower_scalar_as_iobuffer_by_case_file(
            case_file;
            type_key_string = "mpc_baseMVA" )[1]

    #---------------------------------------------------


    mpc_bus = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file;
        type_key= "mpc.bus" )


    mpc_gencost =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.gencost" )


    mpc_gen = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file;
        type_key= "mpc.gen" )


    mpc_branch = get_matpower_mpc_type_iobuffer_by_case_file(
        case_file;
        type_key= "mpc.branch" )

    #---------------------------------------------------
    
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_mpc( mpc_bus  )
    
    (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx ) =
         NamedTupleTools.select(
             net_nodes_type_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))

x    #---------------------------------------------------
    
    no_nodes = length( all_nodes_idx )

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

    no_gens = length( gens_nodes_idx )

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
                   mpc_gencost.c0 )]
    
    #---------------------------------------------------

    gens_Pmax = mpc_gen.Pmax[ gens_nodes_idx ]

    gens_Pmin = mpc_gen.Pmin[ gens_nodes_idx ]

    gens_Qmax = mpc_gen.Qmax[ gens_nodes_idx ]

    gens_Qmin = mpc_gen.Qmin[ gens_nodes_idx ]

    #---------------------------------------------------

    nodes_Pd = mpc_bus.Pd[nodes_with_demands_idx]

    nodes_Qd = mpc_bus.Qd[nodes_with_demands_idx]

    #---------------------------------------------------
    #---------------------------------------------------

    # Real generation: lower (`lb`) and upper (`ub`) bounds

    P_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmin, no_nodes)
    P_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmax, no_nodes)

    # Reactive generation: lower (`lb`) and upper (`ub`) bounds

    Q_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmin, no_nodes)

    Q_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmax, no_nodes)

    # Power demand levels (real, reactive, and complex form)

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
    set_attribute(model, "tol_feas",    1e-3)
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
        [0.5, real(W[i, i]), real(V[i]), imag(V[i])] in RotatedSecondOrderCone() )

    @constraint(
        model,
        [i in 1:no_nodes],
        S_G[i] - S_Demand[i] ==
            LinearAlgebra.tr((conj(Ybus) * E(i, i)) * W), )

    #---------------------------------------------------

    P_G = real(S_G)

    #---------------------------------------------------

    @objective(
        model,
        Min,
        
        # (0.11 * P_G[1]^2 + 5 * P_G[1] + 150) +
        # (0.085 * P_G[2]^2 + 1.2 * P_G[2] + 600) +
        # (0.1225 * P_G[3]^2 + P_G[3] + 335)

        sum( C[i] + Q[i] * P_G[i] + P[i,i] * P_G[i]^2
            for i in 1:no_nodes ) ,)

    
    #---------------------------------------------------

    optimize!(model)

    #---------------------------------------------------

    # assert_is_solved_and_feasible(model; allow_almost = true)

    sdp_relaxation_lower_bound =
        round(objective_value(model); digits = 2)

    println( "Objective value (W & V relax. lower bound): " *
        "$sdp_relaxation_lower_bound", )

    #---------------------------------------------------

    W_1 = SparseArrays.sparse( round.(value.(W); digits = 2) )

    #---------------------------------------------------

    return DataFrames.DataFrame(
        ; Bus = 1:no_nodes,
        Magnitude = round.(abs.(value.(V)); digits = 2),
        AngleDeg = round.(rad2deg.(angle.(value.(V)));
                          digits = 2), )

    #---------------------------------------------------
    
end

#---------------------------------------------------
# Test
#---------------------------------------------------

function driver_opf()
    case_9_file = "/Users/yusufaa/myDocument/Research/DSI-NRF/dev_packages/Quick-Dynamic-Test/src/model3-dynamics/../../src/data-dir/matpower_data/case9.m"

    df_opf = opf( case_9_file )

    df_opf_by_relaxation =
        opf_by_relaxation( case_9_file )
end

#---------------------------------------------------
####################################################
#---------------------------------------------------


# case_9_matpower_dict = get_matpower_dict_by_case_file(
#     case_9_file )

# matpower_dict_keys  =
#     keys( get_matpower_dict_by_case_file(
#         case_9_matpower_dict ) )


function driver_examples_from_jump()
    
    #---------------------------------------------------
    #---------------------------------------------------
    # Economic dispatch
    #---------------------------------------------------
    #---------------------------------------------------

    # https://jump.dev/JuMP.jl/stable/tutorials/applications/power_systems/


    # thermal generators:

    function ThermalGenerator(
        min::Float64,
        max::Float64,
        fixed_cost::Float64,
        variable_cost::Float64,
    )
        return (
            min = min,
            max = max,
            fixed_cost = fixed_cost,
            variable_cost = variable_cost,
        )
    end

    generators = [
        ThermalGenerator(0.0, 1000.0, 1000.0, 50.0),
        ThermalGenerator(300.0, 1000.0, 0.0, 100.0),
    ]

    # A wind generator

    WindGenerator(variable_cost::Float64) =
        (variable_cost = variable_cost,)

    wind_generator = WindGenerator(50.0)

    # a scenario

    function Scenario(demand::Float64, wind::Float64)
        return (demand = demand, wind = wind)
    end

    scenario = Scenario(1500.0, 200.0)

    # economic_dispatch,


    function solve_economic_dispatch(
        generators::Vector, wind, scenario)
        
        # Define the economic dispatch (ED) model
        
        model = Model(HiGHS.Optimizer)
        
        set_silent(model)
        
        # Define decision variables
        # power output of generators
        
        N = length(generators)
        
        @variable(model,
                  generators[i].min <= g[i = 1:N] <=
                      generators[i].max)
        
        # wind power injection
        
        @variable(model, 0 <= w <= scenario.wind)
        
        # Define the objective function
        
        @objective(
            model,
            Min,
            sum(generators[i].variable_cost * g[i]
                for i in 1:N) +
                    wind.variable_cost * w, )
        
        # Define the power balance constraint
        @constraint(model,
                    sum(g[i]
                        for i in 1:N) + w == scenario.demand )

        # Solve statement

        optimize!(model)

        # assert_is_solved_and_feasible(model)

        # return the optimal value of the objective function and its minimizers

        return (
            g = value.(g),
            w = value(w),
            wind_spill = scenario.wind - value(w),
            total_cost = objective_value(model), )
    end

    solution = solve_economic_dispatch(
        generators, wind_generator, scenario);

    println("Dispatch of Generators: ", solution.g, " MW")
    println("Dispatch of Wind: ", solution.w, " MW")
    println("Wind spillage: ", solution.wind_spill, " MW")
    println("Total cost: \$", solution.total_cost)

    #---------------------------------------------------

    # Economic dispatch with adjustable incremental costs

    function scale_generator_cost(g, scale)
        
        return ThermalGenerator(
            g.min, g.max, g.fixed_cost,
            scale * g.variable_cost)
    end

    start = time()
    
    c_g_scale_df = DataFrames.DataFrame(
        ; # Scale factor
        scale = Float64[],
        
        # Dispatch of Generator 1 [MW]
        dispatch_G1 = Float64[],
        
        # Dispatch of Generator 2 [MW]
        dispatch_G2 = Float64[],
        
        # Dispatch of Wind [MW]
        dispatch_wind = Float64[],
        
        # Spillage of Wind [MW]
        spillage_wind = Float64[],
        
        # Total cost [$]
        total_cost = Float64[], )
    
    for c_g1_scale in 0.5:0.1:3.0

        # Update the incremental cost of the first
        # generator at every iteration.
        
        new_generators =
            scale_generator_cost.(
                generators, [c_g1_scale, 1.0])
        
        # Solve the economic-dispatch problem with
        # the updated incremental cost
        
        sol = solve_economic_dispatch(
            new_generators, wind_generator, scenario)
        
        push!( c_g_scale_df,
               (c_g1_scale, sol.g[1], sol.g[2],
                sol.w, sol.wind_spill, sol.total_cost), )
    end
    
    print(string("elapsed time: ", time() - start, " seconds"))

    c_g_scale_df

    #---------------------------------------------------

    # Modifying the JuMP model in-place

    function solve_economic_dispatch_inplace(
        generators::Vector,
        wind,
        scenario,
        scale::AbstractVector{Float64},
    )
        obj_out = Float64[]
        w_out   = Float64[]
        g1_out  = Float64[]
        g2_out  = Float64[]
        
        # This function only works for two generators
        
        @assert length(generators) == 2
        
        model = Model(HiGHS.Optimizer)
        
        set_silent(model)
        
        N = length(generators)
        
        @variable(model,
                  generators[i].min <= g[i = 1:N] <=
                      generators[i].max)
        
        @variable(model, 0 <= w <= scenario.wind)
        
        @objective(
            model,
            Min,
            sum(generators[i].variable_cost * g[i]
                for i in 1:N) + wind.variable_cost * w, )
        
        @constraint(model, sum(g[i]
                               for i in 1:N) + w ==
                                   scenario.demand)
        for c_g1_scale in scale
            @objective(
                model,
                Min,
                c_g1_scale * generators[1].variable_cost * g[1] +
                generators[2].variable_cost * g[2] +
                wind.variable_cost * w, )
            
            optimize!(model)
            
            # assert_is_solved_and_feasible(model)
            
            push!(obj_out, objective_value(model))
            
            push!(w_out, value(w))
            
            push!(g1_out, value(g[1]))
            
            push!(g2_out, value(g[2]))
            
        end
        
        df = DataFrames.DataFrame(
            ; scale = scale,
            dispatch_G1 = g1_out,
            dispatch_G2 = g2_out,
            dispatch_wind = w_out,
            spillage_wind = scenario.wind .- w_out,
            total_cost = obj_out, )
        
        return df
    end

    start = time()
    inplace_df = solve_economic_dispatch_inplace(
        generators,
        wind_generator,
        scenario,
        0.5:0.1:3.0,
    )
    print(string("elapsed time: ", time() - start, " seconds"))

    inplace_df

    #---------------------------------------------------

    # Inefficient usage of wind generators

    demand_scale_df = DataFrames.DataFrame(;
        demand = Float64[],
        dispatch_G1 = Float64[],
        dispatch_G2 = Float64[],
        dispatch_wind = Float64[],
        spillage_wind = Float64[],
        total_cost = Float64[],
    )

    function scale_demand(scenario, scale)
        
        return Scenario(
            scale * scenario.demand, scenario.wind)
    end

    for demand_scale in 0.2:0.1:1.4
        new_scenario = scale_demand(
            scenario, demand_scale)
        
        sol = solve_economic_dispatch(
            generators, wind_generator, new_scenario)
        
        push!(
            demand_scale_df,
            (
                new_scenario.demand,
                sol.g[1],
                sol.g[2],
                sol.w,
                sol.wind_spill,
                sol.total_cost,
            ),
        )
    end

    demand_scale_df

    #---------------------------------------------------

    dispatch_plot = StatsPlots.@df(
        demand_scale_df,
        Plots.plot(
            :demand,
            [:dispatch_G1, :dispatch_G2],
            labels = ["G1" "G2"],
            title = "Thermal Dispatch",
            legend = :bottomright,
            linewidth = 3,
            xlabel = "Demand",
            ylabel = "Dispatch [MW]",
        ),
    )

    wind_plot = StatsPlots.@df(
        demand_scale_df,
        Plots.plot(
            :demand,
            [:dispatch_wind, :spillage_wind],
            labels = ["Dispatch" "Spillage"],
            title = "Wind",
            legend = :bottomright,
            linewidth = 3,
            xlabel = "Demand [MW]",
            ylabel = "Energy [MW]",
        ),
    )

    Plots.plot(dispatch_plot, wind_plot)

    #---------------------------------------------------

    # Unit commitment

    function solve_unit_commitment(
        generators::Vector, wind, scenario)
        
        model = Model(HiGHS.Optimizer)

        set_silent(model)

        N = length(generators)

        @variable(model, 0 <= g[i = 1:N] <= generators[i].max)

        @variable(model, 0 <= w <= scenario.wind)

        @constraint(model, sum(g[i] for i in 1:N) + w ==
            scenario.demand)
        
        # !!! New: add binary on-off variables for each generator
        @variable(model, u[i = 1:N], Bin)
        
        @constraint(model, [i = 1:N], g[i] <=
            generators[i].max * u[i])
        
        @constraint(model, [i = 1:N], g[i] >=
            generators[i].min * u[i])
        
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
            total_cost = objective_value(model),
        )
    end

    solution = solve_unit_commitment(
        generators, wind_generator, scenario)

    println("Dispatch of Generators: ", solution.g, " MW")
    println("Commitments of Generators: ", solution.u)
    println("Dispatch of Wind: ", solution.w, " MW")
    println("Wind spillage: ", solution.wind_spill, " MW")
    println("Total cost: \$", solution.total_cost)

    # Unit commitment as a function of demand

    uc_df = DataFrames.DataFrame(;
        demand = Float64[],
        commitment_G1 = Float64[],
        commitment_G2 = Float64[],
        dispatch_G1 = Float64[],
        dispatch_G2 = Float64[],
        dispatch_wind = Float64[],
        spillage_wind = Float64[],
        total_cost = Float64[],
    )

    for demand_scale in 0.2:0.1:1.4
        new_scenario = scale_demand(scenario, demand_scale)
        sol = solve_unit_commitment(generators, wind_generator, new_scenario)
        if sol.status == OPTIMAL
            push!(
                uc_df,
                (
                    new_scenario.demand,
                    sol.u[1],
                    sol.u[2],
                    sol.g[1],
                    sol.g[2],
                    sol.w,
                    sol.wind_spill,
                    sol.total_cost,
                ),
            )
        end
        println("Status: $(sol.status) for demand_scale = $(demand_scale)")
    end

    uc_df



    #---------------------------------------------------

    commitment_plot = StatsPlots.@df(
        uc_df,
        Plots.plot(
            :demand,
            [:commitment_G1, :commitment_G2],
            labels = ["G1" "G2"],
            title = "Commitment",
            legend = :bottomright,
            linewidth = 3,
            xlabel = "Demand [MW]",
            ylabel = "Commitment decision {0, 1}",
        ),
    )

    dispatch_plot = StatsPlots.@df(
        uc_df,
        Plots.plot(
            :demand,
            [:dispatch_G1, :dispatch_G2, :dispatch_wind],
            labels = ["G1" "G2" "Wind"],
            title = "Dispatch [MW]",
            legend = :bottomright,
            linewidth = 3,
            xlabel = "Demand",
            ylabel = "Dispatch [MW]",
        ),
    )

    Plots.plot(commitment_plot, dispatch_plot)


    #---------------------------------------------------

    # Nonlinear economic dispatch

    """
        thermal_cost_function(g)

    A user-defined thermal cost function in pure-Julia! You can include
    nonlinearities, and even things like control flow.

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

    function solve_nonlinear_economic_dispatch(
        generators::Vector,
        wind,
        scenario;
        silent::Bool = false, )
        model = Model(Ipopt.Optimizer)
        if silent
            set_silent(model)
        end
        
        @operator(model, op_tcf, 1, thermal_cost_function)
        N = length(generators)
        
        @variable(model, generators[i].min <= g[i = 1:N] <=
            generators[i].max)
        @variable(model, 0 <= w <= scenario.wind)
        
        @objective(
            model,
            Min,
            sum(generators[i].variable_cost * op_tcf(g[i])
                for i in 1:N) +
                    wind.variable_cost * w, )
        
        @constraint(model, sum(g[i]
                               for i in 1:N) + sqrt(w) ==
                                   scenario.demand)
        optimize!(model)
        # assert_is_solved_and_feasible(model)
        return (
            g = value.(g),
            w = value(w),
            wind_spill = scenario.wind - value(w),
            total_cost = objective_value(model),
        )
    end

    solution = solve_nonlinear_economic_dispatch(
            generators, wind_generator, scenario)


    wind_cost = 0.0:1:100
    wind_dispatch = Float64[]
    for c in wind_cost
        sol = solve_nonlinear_economic_dispatch(
            generators,
            WindGenerator(c),
            scenario;
            silent = true,
        )
        push!(wind_dispatch, sol.w)
    end

    Plots.plot(
        wind_cost,
        wind_dispatch;
        xlabel = "Cost",
        ylabel = "Dispatch [MW]",
        label = false,
    )



    #---------------------------------------------------

    #---------------------------------------------------
    #---------------------------------------------------
    # mixed complementarity problems
    #---------------------------------------------------
    #---------------------------------------------------

    # https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/complementarity/#Electricity-consumption

    #---------------------------------------------------
    # Electricity consumption
    #---------------------------------------------------

    # https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/complementarity/#Electricity-consumption


    I = 90_000                     # Annualized capital cost
    C = 60                         # Operation cost per MWh
    τ = 8_760                      # Hours per year
    θ = [0.2, 0.2, 0.2, 0.2, 0.2]  # Scenario probabilities

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
    @constraint(model, [ω = 1:5], P[ω] - (A[ω] - B * Q[ω]) ⟂ Q[ω])

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
