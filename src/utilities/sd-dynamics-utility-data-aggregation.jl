# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247


#####################################################
# ---------------------------------------------------
# OPF Data aggregation utility functions 
# ---------------------------------------------------
#####################################################

"""
    get_renewable_energy_net_optimisation_parameters(
        case_file,

        wind_gens_cost_scale,    
        solar_gens_cost_scale,

        wind_gens_capacity_scale,    
        solar_gens_capacity_scale,

        active_power_demand_deviation_scale,
        reactive_power_demand_deviation_scale)

Return renewable energy network optimisation parameters for optimal power flow, unit commitment or economic dispatch


"""
function get_renewable_energy_net_optimisation_parameters(
    case_file,
    
    wind_gens_cost_scale,    
    solar_gens_cost_scale,
    
    wind_gens_capacity_scale,    
    solar_gens_capacity_scale,

    active_power_demand_deviation_scale,
    reactive_power_demand_deviation_scale)

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


    mpc_gen =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.gen" )


    mpc_branch =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.branch" )


    #---------------------------------------------------

    edges_orientation =
        get_edges_orientation_by_generic(
            mpc_branch.fbus,
            mpc_branch.tbus )
    
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


    gens_installed_capacity = copy(mpc_gen.Pmax) 
    
    gens_Pmax = mpc_gen.Pmax

    gens_Pmin = mpc_gen.Pmin

    gens_Qmax = mpc_gen.Qmax

    gens_Qmin = mpc_gen.Qmin

    #---------------------------------------------------

    sch_Pg = mpc_gen.Pg

    # check if scheduled power is 0, if it is, set Pmax = 0
    
    gens_Pmax = [(sch_Pg[idx] == 0.0 || sch_Pg[idx] == 0 ) ?
        sch_Pg[idx] : pg for (idx, pg) in
                 enumerate(gens_Pmax)]
    
    #---------------------------------------------------
    
    wind_gens_Pmin =
        zeros(length(gens_Pmin))

    solar_gens_Pmin =
        zeros(length(gens_Pmin))

    
    # wind_gens_Pmax =
    #     wind_gens_power_forecast
    
    # solar_gens_Pmax =
    #     solar_gens_power_forecast
    
    wind_gens_Pmax =
        get_wind_gens_power_forecast(
            mpc_gen.Pmax;
            scale_factor =
                wind_gens_capacity_scale )    
    solar_gens_Pmax =
        get_solar_gens_power_forecast(
            mpc_gen.Pmax;
            scale_factor =
                solar_gens_capacity_scale )
    
    #---------------------------------------------------

    base_Pd = mpc_bus.Pd[nodes_with_demands_idx]

    base_Qd = mpc_bus.Qd[nodes_with_demands_idx]
    
    nodes_Pd = get_power_demand_scenario(
        base_Pd;
        scale_factor =
            active_power_demand_deviation_scale)

    nodes_Qd =  get_power_demand_scenario(
        base_Qd;
        scale_factor =
            reactive_power_demand_deviation_scale)

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
            nodes_with_demands_idx,
            # nodes_Qd,
            base_Qd,
            no_nodes)

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

    branch_r = mpc_branch.r
    
    branch_x = mpc_branch.x

    branch_b = mpc_branch.b

    bus_Gs = mpc_bus.Gs
    bus_Bs = mpc_bus.Bs
    
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

    return (;            
            nodes_Pd,
            nodes_Qd,

            P_Demand,
            Q_Demand,
            S_Demand,
            
            solar_gens_Pmax,
            wind_gens_Pmax,
            
            P_wind_Gen_ub,
            P_solar_Gen_ub,
            
            mpc_baseMVA,

            no_nodes,
            no_gens,
            no_edges,

            dyn_pf_fun_kwd_net_idxs,
            
            gens_cost_coeff_ascen,
            gens_cost_coeff_decen,
            wind_gens_cost_coeff,
            solar_gens_cost_coeff,
            
            solar_gens_Pmin,            
            wind_gens_Pmin,

            gens_installed_capacity,
            
            
            gens_Pmin,            
            gens_Pmax,
            
            gens_Qmin,            
            gens_Qmax,

            base_Pd,
            base_Qd,            

            P_wind_Gen_lb,
            P_solar_Gen_lb,

            P_Gen_lb,
            P_Gen_ub,
            
            Q_Gen_lb,
            Q_Gen_ub,

            branch_r,
            branch_x,
            branch_b,
            bus_Gs,
            bus_Bs,

            edges_orientation,
            
            inv_τ,
            θ_shift,

            ys,
            y_c,
            y_sh,

            Y_sh,

            Cf,
            Ct,
            
            Yf,
            Yt,
            Ybus)
        
end


"""
    get_opf_net_optimisation_parameters(
        case_file )

Returns network optimisation parameters for optimal power flow, unit commitment or economic dispatch
"""
function get_opf_net_optimisation_parameters(
        case_file )

    mpc_baseMVA =
        get_matpower_scalar_as_iobuffer_by_case_file(
            case_file;
            type_key_string =
                "mpc_baseMVA" )[1]

    #---------------------------------------------------

    mpc_bus =
        get_matpower_mpc_type_iobuffer_by_case_file(
            case_file;
            type_key= "mpc.bus" )


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


    edges_orientation =
        get_edges_orientation_by_generic(
            mpc_branch.fbus,
            mpc_branch.tbus )
    
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
        [a_node for ( a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 3 || node_type == 2  ]


    non_gens_nodes_idx =
        [a_node for ( a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 1 ]


    gens_nodes_with_loc_loads_idx =
        [a_node for ( a_node, node_type, node_Pd) in
             zip( mpc_bus.bus_i, mpc_bus.type,
                  mpc_bus.Pd)
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
        (;slack_gens_nodes_idx,
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
                   mpc_gencost.c0 )]
    
    #---------------------------------------------------
    
    # gens_Pmax = mpc_gen.Pmax[ gens_nodes_idx ]

    # gens_Pmin = mpc_gen.Pmin[ gens_nodes_idx ]

    # gens_Qmax = mpc_gen.Qmax[ gens_nodes_idx ]

    # gens_Qmin = mpc_gen.Qmin[ gens_nodes_idx ]

    
    gens_installed_capacity =
        mpc_gen.Pmax ./ mpc_baseMVA

    gens_Pmax = copy(gens_installed_capacity)

    gens_Pmin = mpc_gen.Pmin ./ mpc_baseMVA

    gens_Qmax = mpc_gen.Qmax ./ mpc_baseMVA

    gens_Qmin = mpc_gen.Qmin ./ mpc_baseMVA


    sch_Pg = mpc_gen.Pg ./ mpc_baseMVA
    
    gens_Pmax =
        [(sch_Pg[idx] == 0.0 || sch_Pg[idx] == 0 ) ?
        sch_Pg[idx] : pg for (idx, pg) in
                 enumerate(gens_Pmax)]
    
    #---------------------------------------------------

    nodes_Pd =
        mpc_bus.Pd[nodes_with_demands_idx] ./
        mpc_baseMVA

    nodes_Qd =
        mpc_bus.Qd[nodes_with_demands_idx] ./
        mpc_baseMVA

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

    branch_fbus = mpc_branch.fbus
    branch_tbus = mpc_branch.tbus
    branch_r = mpc_branch.r
    branch_x = mpc_branch.x
    branch_b = mpc_branch.b
    
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
        mpc_baseMVA ./ (mpc_branch.r +
        im * mpc_branch.x) 
       

    y_c =
        1 / 2 * (im *  mpc_branch.b) *
        mpc_baseMVA

    y_sh =  (mpc_bus.Gs .+ im *  mpc_bus.Bs
             )/mpc_baseMVA

    inv_τ =
        [ (a_ratio == 0.0 || a_ratio == 0) ? 1.0 : 1/a_ratio
              for a_ratio in
                  mpc_branch.ratio  ]
    
    θ_shift = mpc_branch.angle
    
    Yff = (ys +  y_c) .* (inv_τ).^2

    Ytf = -ys .* ( inv_τ ./ exp.(im * θ_shift) )

    Yft = -ys .* ( inv_τ ./ exp.(-im * θ_shift) )

    Ytt = ( ys +  y_c )

    Y_sh = SparseArrays.spdiagm( y_sh )

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

    return (;
            mpc_baseMVA,
            mpc_bus,
            mpc_gencost,
            mpc_gen,
            mpc_branch,

            slack_gens_nodes_idx,
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
            Ybus )
    
end



"""
    get_system_net_data_wt_static_parameters(
        case_name;
        <keyword arguments>)


It is used to simplify static model parameters or data of a system that are needed for static analyses.


# Arguments

-`case_name`: the case name
- `script_dir::String=""`: working folder
- `data_dir::String=""`: data folder where all cases are stored
- `json_net_data_by_components_file`: the network data file
- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,
`pf_alg`: power flow solver

"""
function get_system_net_data_wt_static_parameters(
    case_name
    ;script_dir = "",
    data_dir ="",
    json_net_data_by_components_file = "",
    components_libs_dir = "",
    basekV              = 1.0,    
    use_pu_in_PQ        = true,
    line_data_in_pu     = true,
    pf_alg              = NewtonRaphson())

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
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath( data_dir,
                  "converted-data",
                  case_name,)
    
    json_case_dir =
        joinpath(case_data_dir,
                  "json")

    
    if (json_net_data_by_components_file == "") || (
        json_net_data_by_components_file == nothing)

        net_data_by_components_file =
            joinpath(
                json_case_dir,
                "net_data_by_components_file.json")
    else

        net_data_by_components_file =
            joinpath(
                json_case_dir,
                json_net_data_by_components_file)

    end

    #----------------------------------------    
    
    (;plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json,
    edge_data_from_json,
    shunt_data_from_json,
    baseMVA_data_from_json) =
       NamedTupleTools.select(
           get_net_data_by_components_from_json_file(
               net_data_by_components_file;
               in_components_type_sym =
                   false ),
           (:plant_generators_data_from_json,
            :plant_loads_data_from_json,
            :plant_transmission_data_from_json,
            :edge_data_from_json,
            :shunt_data_from_json,
            :baseMVA_data_from_json))

    #------------------------------------------------

    baseMVA = baseMVA_data_from_json

    #------------------------------------------------
    # edges
    #------------------------------------------------

    # (branches_fbus,
    #  branches_tbus,
    #  branches_r,
    #  branches_x,
    #  branches_b,
    #  branches_ratio,
    #  branches_angle,
    #  branches_type) =
    #       get_edges_ftbus_and_generic_data_by_json(
    #          edge_data_from_json )

    # (branches_fbus,
    #  branches_tbus) =
    #      get_edges_fbus_tbus_by_json(
    #          edge_data_from_json)

    # (edges_r, edges_x, edges_b,
    #  edges_ratio, edges_angle, edges_type) =
    #      get_edges_generic_data_by_json(
    #          edge_data_from_json )

    (edges_fbus,
    edges_tbus,
    edges_r,
    edges_x,
    edges_b,
    edges_ratio,
    edges_angle,
    edges_type) =
        get_edges_ftbus_and_generic_data_by_json(
            edge_data_from_json )

    # (Gs, Bs) =
    #     get_nodes_shunts_Gs_and_Bs_by_json(
    #         shunt_data_from_json)

    (shunt_idx,
    Gs,
    Bs) =
       get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
           shunt_data_from_json)

    #------------------------------------------------

    (;edges_orientation,
    edges_Ybr_cal,
    Ybr_cal_and_edges_orientation,
    Ynet_wt_nodes_idx_wt_adjacent_nodes) =
        NamedTupleTools.select(
            get_transmission_network_parameters_by_json(
                plant_generators_data_from_json,
                plant_loads_data_from_json,
                plant_transmission_data_from_json,
                edge_data_from_json,
                shunt_data_from_json;
                baseMVA =
                    baseMVA,
                basekV =
                    basekV,
                use_pu_in_PQ =
                    use_pu_in_PQ,
                line_data_in_pu =
                    line_data_in_pu ),
            (:edges_orientation,
             :edges_Ybr_cal,
             :Ybr_cal_and_edges_orientation,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #------------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
    Ynet_real_imag_Idxs_in_flattend ) =
        NamedTupleTools.select(
            get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                getproperty(
                    Ynet_wt_nodes_idx_wt_adjacent_nodes,
                    :Ynet)),
            (:Ynet_rows_Idxs_in_flattend,
             :Ynet_real_imag_Idxs_in_flattend) )

    #------------------------------------------------
    # nodes
    #------------------------------------------------

    (;
    P_gens,
    Q_gens,
    P_non_gens,
    Q_non_gens,
    P_g_loc_load,
    Q_g_loc_load,
    loc_load_exist ) =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA =
                baseMVA,
            use_pu_in_PQ =
                use_pu_in_PQ )

    #------------------------------------------------

    ode_gens_para_selections  =
       (:H,
        :D,
        :X_d,
        :X_q,                  
        :X_d_dash,
        :X_q_dash,
        :T_d_dash,
        :T_q_dash, :Sn )

    ode_gens_para_sequence_order =
       (:components_data,
        :gen)

    ode_gens_generic_selections =
       (:H,
        :D,
        :ra,
        :xℓ,
        :X_d,
        :X_q,

        :X_d_dash,
        :X_q_dash,

        :X_d_2dash,
        :X_q_2dash,

        :T_d_dash,
        :T_q_dash,

        :Sn,

        :T_d_2dash,
        :T_q_2dash,

        :vh,
        :P,
        :Q,

        :Pmin,
        :Pmax,
        :Qmin,
        :Qmax,
        :vmin,
        :vmax )

    ode_gens_generic_sequence_order =
       (:components_data, :gen)

    govs_and_avrs_sequence_order =
       ( :components_data,)

    govs_and_avrs_selections =
       ( :gov, :avr )

    #----------------------------------------

    ode_gens_generic_para =
        get_ode_gens_generic_para(
            plant_generators_data_from_json;
            sequence_order =
                ode_gens_generic_sequence_order,
            selections =
                ode_gens_generic_selections)

    #------------------------------------------------
    #------------------------------------------------

    (;generic_gens_para,
    generic_govs_para,
    generic_avrs_para) =
        get_generic_gens_avr_gov_para(
            plant_generators_data_from_json;
            gens_sequence_order =
                ode_gens_generic_sequence_order,
           gens_selections =
               ode_gens_generic_selections,
           govs_and_avrs_sequence_order =
               govs_and_avrs_sequence_order,
           govs_and_avrs_selections =
               govs_and_avrs_selections)

    #------------------------------------------------
    #------------------------------------------------

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
                    :X_q_dash, :Sn ),
               vec_datatype = Float64 ),
           (:ra,
            :X_d,
            :X_q,     
            :X_d_dash,
            :X_q_dash, :Sn ) )

    #------------------------------------------------
    #------------------------------------------------

    opf_generic_gens_para =
       NamedTupleTools.select(
           get_selected_vec_nt_to_vec_vec(
               generic_gens_para,
               nothing;
               selections =
                   (:vh,
                    :P, :Q,
                    :Pmin, :Pmax,
                    :Qmin, :Qmax,
                    :vmin, :vmax,
                    :Sn ),
               vec_datatype = Float64 ),
           (:vh,
            :P, :Q,
            :Pmin, :Pmax,
            :Qmin, :Qmax,
            :vmin, :vmax,
            :Sn ) )

    #----------------------------------------

    ode_gens_para =
       NamedTupleTools.select(
           ode_gens_generic_para,
           (:H,
            :D,
            :X_d,
            :X_q,
            :X_d_dash,
            :X_q_dash,
            :T_d_dash,
            :T_q_dash, :Sn))

    net_nodes_type_idxs =
       get_net_nodes_type_idxs_by_json(
           plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json )

    dyn_pf_fun_kwd_net_idxs =
       NamedTupleTools.select(
           net_nodes_type_idxs,
           (:slack_gens_nodes_idx,
            :non_slack_gens_nodes_idx,
            :gens_nodes_idx,
            :non_gens_nodes_idx,
            :gens_with_loc_load_idx,
            :gens_nodes_with_loc_loads_idx,
            :all_nodes_idx))

    dyn_pf_fun_kwd_n2s_idxs =
       NamedTupleTools.select(
           get_dict_net_streamlined_idx_by_nodes_type_idxs(
               net_nodes_type_idxs ),
           (:n2s_slack_gens_idx,
            :n2s_non_slack_gens_idx,
            :n2s_gens_idx,
            :n2s_non_gens_idx,
            :n2s_gens_with_loc_load_idxs,
            :n2s_all_nodes_idx))

    nodes_with_demands_idx =
       getproperty(net_nodes_type_idxs,
                   :nodes_with_demands_idx )

    #----------------------------------------

    all_nodes_idx =
        getproperty(
            dyn_pf_fun_kwd_net_idxs,
            :all_nodes_idx )

    n2s_all_nodes_idx =
        getproperty(
            dyn_pf_fun_kwd_n2s_idxs,
            :n2s_all_nodes_idx )
    
    #----------------------------------------

    # dyn_pf_mismatch_vars_kwd_para = 
    #     (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
    #       ode_gens_para )

    sta_pf_PQ_para =
       get_pf_PQ_param_by_json(
           plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json;
           baseMVA =
               baseMVA,
           use_pu_in_PQ =
               use_pu_in_PQ)


    gens_vh_slack_θh_para =
       get_gens_vh_slack_θh_para_by_json(
           plant_generators_data_from_json )


    sta_pf_vars_and_paras_idx =
       get_sta_pf_vars_and_paras_idx_by_json(
           plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json )


    pf_sta_ΔPQ_mismatch_parameters =
       get_pf_sta_ΔPQ_mismatch_parameters_by_json(
           plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json,
           edge_data_from_json,
           shunt_data_from_json;
           baseMVA = baseMVA,
           basekV = basekV,
           use_pu_in_PQ = use_pu_in_PQ,
           line_data_in_pu = line_data_in_pu)

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
    # Powerflow func and prob
    #----------------------------------------

    kwd_sta_sta_ΔPQ_sol_by_json =
       (;
        pf_alg,
        pf_kw_para,
        red_vh_Idxs,
        red_θh_Idxs,
        sta_red_vh_θh_0) 

    #----------------------------------------
    # states_Idx, syms and functions
    #----------------------------------------
    states_Idx_syms_wt_funcs =
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

             :state_vars_and_i_dq_wt_fault_Idx_in_state,

             :state_algebraic_vars_Idx_in_state,
             :state_algebraic_vars_wt_fault_Idx_in_state,

             :model_mass_matrix,     
             :model_bool_dae_vars,

             :ode_gens_mass_matrix,
             :ode_gens_bool_dae_vars,

             :Pg_Qg_Png_Qng_Pll_Qll_Idx,
             :Png_Qng_Pll_Qll_Idx,
             :Pg_Png_Qng_Idx,

             :dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

             :dyn_vh_id_iq_V_ref_Tm_Idx,
             :dyn_V_ref_Tm_id_iq_vh_Idx,

             :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
             :dyn_vh_id_iq_ωref0_vref0_porder0_Idx,

             :ωref0_vref0_porder0_id_iq_vh_Idx,

             :dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
             :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
             :pf_vh_θh_idx_and_idx2Idx,

             :dyn_pf_flat_vh_flat_θh_Idx,

             :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

             :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

             :dyn_pf_vh_vhf_Idx,
             :dyn_pf_θh_θhf_Idx,

             :system_states_idx_kwd_para,
             :system_paras_idx_kwd_para,
             :plants_dyn_fun_idx_kwd_para,
             :plants_algeb_fun_idx_kwd_para,

             :id_iq_pg_vh_Idx,

             :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
             :scale_Pg_Png_Qng_Idx,
             :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx))

    dyn_pf_flat_vh_flat_θh_id_iq_Idx =
        getproperty(states_Idx_syms_wt_funcs,
                    :dyn_pf_flat_vh_flat_θh_id_iq_Idx)
    
    #----------------------------------------
    # Results    
    #----------------------------------------

    generic_red_sol_kwd_para =
       (;Ybr_cal_and_edges_orientation,
         sta_pf_PQ_para,
         ode_gens_generic_para,
         pf_kw_para ) 

    generic_dyn_sol_kwd_para =
       (;loc_load_exist,
        sta_pf_PQ_para,
        ode_gens_generic_para,
        dyn_pf_flat_vh_flat_θh_id_iq_Idx,
        dyn_pf_fun_kwd_n2s_idxs,
        dyn_pf_fun_kwd_net_idxs)


    return(; plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json,
           edge_data_from_json,
           shunt_data_from_json,
           baseMVA,

           dyn_pf_fun_kwd_net_idxs,
           dyn_pf_fun_kwd_n2s_idxs,

           all_nodes_idx,
           n2s_all_nodes_idx,

           pf_sta_ΔPQ_mismatch_parameters,
           kwd_sta_sta_ΔPQ_sol_by_json,
           generic_red_sol_kwd_para,
           generic_dyn_sol_kwd_para,
           states_Idx_syms_wt_funcs )
   
end


"""
    get_system_net_static_data(
        case_name;
        <keyword arguments>)


It is used to simplify static model parameters or data of a system that are needed for static analyses.


# Arguments

-`case_name`: the case name
- `script_dir::String=""`: working folder
- `data_dir::String=""`: data folder where all cases are stored
- `json_net_data_by_components_file`: the network data file
- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,

`use_init_u0`: the boolean variable that determines if initial state u0 should be used in a power flow.
`use_nlsolve`: the boolean variable that determines if `nlsolve` should be used in power flow.

`pf_alg`: power flow solver
`no_lines_fault`: the number of faults
"""
function get_system_net_static_data(
    case_name
    ;script_dir = "",
    data_dir = "",
    json_net_data_by_components_file = "",
    components_libs_dir = "",
    basekV              = 1.0,    
    use_pu_in_PQ        = true,
    line_data_in_pu     = true,
    pf_alg              = NewtonRaphson(),

    no_lines_fault = 1 )

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
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath(data_dir,
                 "converted-data",
                 case_name)
    
    json_case_dir =
        joinpath(case_data_dir, "json")

    if (json_net_data_by_components_file == "") || (
        json_net_data_by_components_file == nothing)

        net_data_by_components_file =
            joinpath(
                json_case_dir,
                "net_data_by_components_file.json")
    else

        net_data_by_components_file =
            joinpath(
                json_case_dir,
                json_net_data_by_components_file)

    end

    #----------------------------------------    
    
    (;plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json,
    edge_data_from_json,
    shunt_data_from_json,
     baseMVA_data_from_json,
     gencost_data_from_json) =
       NamedTupleTools.select(
           get_net_data_by_components_from_json_file(
               net_data_by_components_file;
               in_components_type_sym =
                   false ),
           (:plant_generators_data_from_json,
            :plant_loads_data_from_json,
            :plant_transmission_data_from_json,
            :edge_data_from_json,
            :shunt_data_from_json,
            :baseMVA_data_from_json,
            :gencost_data_from_json))

    #------------------------------------------------

    baseMVA = baseMVA_data_from_json

    #------------------------------------------------
    # edges
    #------------------------------------------------

    # (branches_fbus,
    #  branches_tbus,
    #  branches_r,
    #  branches_x,
    #  branches_b,
    #  branches_ratio,
    #  branches_angle,
    #  branches_type) =
    #       get_edges_ftbus_and_generic_data_by_json(
    #          edge_data_from_json )

    # (branches_fbus,
    #  branches_tbus) =
    #      get_edges_fbus_tbus_by_json(
    #          edge_data_from_json)

    # (edges_r, edges_x, edges_b,
    #  edges_ratio, edges_angle, edges_type) =
    #      get_edges_generic_data_by_json(
    #          edge_data_from_json )

    (edges_fbus,
    edges_tbus,
    edges_r,
    edges_x,
    edges_b,
    edges_ratio,
    edges_angle,
    edges_type) =
        get_edges_ftbus_and_generic_data_by_json(
            edge_data_from_json )

    # (Gs, Bs) =
    #     get_nodes_shunts_Gs_and_Bs_by_json(
    #         shunt_data_from_json)

    (shunt_idx,
    Gs,
    Bs) =
       get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
           shunt_data_from_json)

    #------------------------------------------------
    
    (;edges_orientation,
    edges_Ybr_cal,
    Ybr_cal_and_edges_orientation,
    Ynet_wt_nodes_idx_wt_adjacent_nodes) =
        NamedTupleTools.select(
            get_transmission_network_parameters_by_json(
                plant_generators_data_from_json,
                plant_loads_data_from_json,
                plant_transmission_data_from_json,
                edge_data_from_json,
                shunt_data_from_json;
                baseMVA =
                    baseMVA,
                basekV =
                    basekV,
                use_pu_in_PQ =
                    use_pu_in_PQ,
                line_data_in_pu =
                    line_data_in_pu ),
            (:edges_orientation,
             :edges_Ybr_cal,
             :Ybr_cal_and_edges_orientation,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    # (;Ynet,
    #   nodes_idx_with_adjacent_nodes_idx) =
    #       NamedTupleTools.select(
    #           Ynet_wt_nodes_idx_wt_adjacent_nodes,
    #           (:Ynet,
    #            :nodes_idx_with_adjacent_nodes_idx))
    
    #------------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend ) =
        NamedTupleTools.select(
            get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                getproperty(
                    Ynet_wt_nodes_idx_wt_adjacent_nodes,
                    :Ynet)),
            (:Ynet_rows_Idxs_in_flattend,
             :Ynet_real_imag_Idxs_in_flattend) )
    
    #------------------------------------------------
    # nodes
    #------------------------------------------------

    (;
    P_gens,
    Q_gens,
    P_non_gens,
    Q_non_gens,
    P_g_loc_load,
    Q_g_loc_load,
    loc_load_exist ) =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA =
                baseMVA,
            use_pu_in_PQ =
                use_pu_in_PQ )

    
    #------------------------------------------------

    gens_para_sequence_order =
       (:components_data,
        :gen)

    gens_generic_selections =
       (
        :Sn,:vh,
        :P, :Q,
        :Pmin, :Pmax,
        :Qmin, :Qmax,
        :vmin, :vmax )

    #----------------------------------------

    generic_each_gen_para =
        get_components_properties_by_json(
            plant_generators_data_from_json;
            sequence_order =
                 gens_para_sequence_order,
             selections =
                  gens_generic_selections )


    "To make sure Vector{NamedTuple} is returned
     instead of Vector{Any}"
    generic_each_gen_para =
        NamedTuple[
            item for item in
                generic_each_gen_para]
    
    #----------------------------------------

    generic_gens_para =
        get_ode_gens_generic_para(
            plant_generators_data_from_json;
            sequence_order =
                gens_para_sequence_order,
            selections =
                gens_generic_selections)

    # "To make sure Vector{NamedTuple} is returned
    #  instead of Vector{Any}"
    # generic_gens_para =
    #     NamedTuple[
    #         item for item in
    #             generic_gens_para]

    #------------------------------------------------

    # selected_generic_gens_para =
    #     get_selected_comps_ode_para_by_json(
    #          plant_generators_data_from_json;
    #          sequence_order =
    #              gens_para_sequence_order,
    #          selections =
    #              gens_generic_selections)
    
    #------------------------------------------------

    # pf_generic_gens_para =
    #    NamedTupleTools.select(
    #        get_selected_vec_nt_to_vec_vec(
    #            generic_gens_para,
    #            nothing;
    #            selections =
    #                (:Sn, :vh,
    #                 :P, :Q,
    #                 :Pmin, :Pmax,
    #                 :Qmin, :Qmax,
    #                 :vmin, :vmax ),
    #            vec_datatype = Float64 ),
    #                    (:Sn, :vh,
    #                     :P, :Q,
    #                     :Pmin, :Pmax,
    #                     :Qmin, :Qmax,
    #                     :vmin, :vmax ) )

    #------------------------------------------------
    #------------------------------------------------

    net_nodes_type_idxs =
       get_net_nodes_type_idxs_by_json(
           plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json )

    dyn_pf_fun_kwd_net_idxs =
       NamedTupleTools.select(
           net_nodes_type_idxs,
           (:slack_gens_nodes_idx,
            :non_slack_gens_nodes_idx,
            :gens_nodes_idx,
            :non_gens_nodes_idx,
            :gens_with_loc_load_idx,
            :gens_nodes_with_loc_loads_idx,
            :all_nodes_idx))

    dyn_pf_fun_kwd_n2s_idxs =
       NamedTupleTools.select(
           get_dict_net_streamlined_idx_by_nodes_type_idxs(
               net_nodes_type_idxs ),
           (:n2s_slack_gens_idx,
            :n2s_non_slack_gens_idx,
            :n2s_gens_idx,
            :n2s_non_gens_idx,
            :n2s_gens_with_loc_load_idxs,
            :n2s_all_nodes_idx))

    nodes_with_demands_idx =
       getproperty(net_nodes_type_idxs,
                   :nodes_with_demands_idx )


    #----------------------------------------

    all_nodes_idx =
        getproperty(
            dyn_pf_fun_kwd_net_idxs,
            :all_nodes_idx )

    n2s_all_nodes_idx =
        getproperty(
            dyn_pf_fun_kwd_n2s_idxs,
            :n2s_all_nodes_idx )
    
    #----------------------------------------

    # dyn_pf_mismatch_vars_kwd_para = 
    #     (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
    #       ode_gens_para )

    sta_pf_PQ_para =
       get_pf_PQ_param_by_json(
           plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json;
           baseMVA =
               baseMVA,
           use_pu_in_PQ =
               use_pu_in_PQ)

    
    gens_vh_slack_θh_para =
       get_gens_vh_slack_θh_para_by_json(
           plant_generators_data_from_json )


    sta_pf_vars_and_paras_idx =
       get_sta_pf_vars_and_paras_idx_by_json(
           plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json )


    pf_sta_ΔPQ_mismatch_parameters =
       get_pf_sta_ΔPQ_mismatch_parameters_by_json(
           plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json,
           edge_data_from_json,
           shunt_data_from_json;
           baseMVA = baseMVA,
           basekV = basekV,
           use_pu_in_PQ = use_pu_in_PQ,
           line_data_in_pu = line_data_in_pu)

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
    # Powerflow func and prob
    #----------------------------------------

    kwd_sta_sta_ΔPQ_sol_by_json =
       (;
        pf_alg,
        pf_kw_para,
        red_vh_Idxs,
        red_θh_Idxs,
        sta_red_vh_θh_0) 
    
    #----------------------------------------
    # Results    
    #----------------------------------------

    generic_red_sol_kwd_para =
        (;Ybr_cal_and_edges_orientation,
          Ynet_wt_nodes_idx_wt_adjacent_nodes,
          sta_pf_PQ_para,
          pf_kw_para) 

    # generic_dyn_sol_kwd_para =
    #    (;loc_load_exist,
    #     sta_pf_PQ_para,
    #     # ode_gens_generic_para,
    #     # dyn_pf_flat_vh_flat_θh_id_iq_Idx,
    #     dyn_pf_fun_kwd_n2s_idxs,
    #     dyn_pf_fun_kwd_net_idxs)


    #----------------------------------------
    # states_Idx, syms and functions
    #----------------------------------------

    static_Idx_and_syms =
        get_static_Idx_and_syms(
            dyn_pf_fun_kwd_net_idxs,
            dyn_pf_fun_kwd_n2s_idxs;
            no_lines_fault =
                no_lines_fault)

    
    return (;plant_generators_data_from_json,
           plant_loads_data_from_json,
           plant_transmission_data_from_json,
           edge_data_from_json,
           shunt_data_from_json,
           baseMVA,
           gencost_data_from_json,

           edges_fbus,
           edges_tbus,
           edges_r,
           edges_x,
           edges_b,
           edges_ratio,
           edges_angle,
           edges_type,

           shunt_idx,
           Gs,
           Bs,
           
           edges_orientation,
           edges_Ybr_cal,
           Ybr_cal_and_edges_orientation,
           Ynet_wt_nodes_idx_wt_adjacent_nodes,
           Ynet_rows_Idxs_in_flattend,
           Ynet_real_imag_Idxs_in_flattend ,

           P_gens,
           Q_gens,
           P_non_gens,
           Q_non_gens,
           P_g_loc_load,
           Q_g_loc_load,
           loc_load_exist,

            generic_gens_para,
            generic_each_gen_para,

           # pf_generic_gens_para,

           net_nodes_type_idxs,
           
           dyn_pf_fun_kwd_net_idxs,
           dyn_pf_fun_kwd_n2s_idxs,

           nodes_with_demands_idx,

           all_nodes_idx,
           n2s_all_nodes_idx,

           sta_pf_PQ_para,
           gens_vh_slack_θh_para,

           sta_pf_vars_and_paras_idx,
           pf_sta_ΔPQ_mismatch_parameters,

            sta_red_vh_θh_0,

            pf_PQ_param,
            pf_kw_para,
           
           kwd_sta_sta_ΔPQ_sol_by_json,
           
           generic_red_sol_kwd_para,
           # generic_dyn_sol_kwd_para,
           
           static_Idx_and_syms )
    
    # merge((;plant_generators_data_from_json,
    #         plant_loads_data_from_json,
    #         plant_transmission_data_from_json,
    #         edge_data_from_json,
    #         shunt_data_from_json,
    #         baseMVA,
    #         gencost_data_from_json,

    #         generic_gens_para,

    #         selected_generic_gens_para,
           
    #          Ynet_rows_Idxs_in_flattend,
    #          Ynet_real_imag_Idxs_in_flattend,

    #          P_gens, Q_gens,
    #          P_non_gens, Q_non_gens,
    #          P_g_loc_load, Q_g_loc_load,
    #          loc_load_exist,

    #         # pf_PQ_param,
    #         red_vh_Idxs,
    #         red_θh_Idxs),

    #       net_nodes_type_idxs,
    #       dyn_pf_fun_kwd_net_idxs,
          
    #       sta_pf_PQ_para,
    #       gens_vh_slack_θh_para,
    #       sta_pf_vars_and_paras_idx,
    #       pf_sta_ΔPQ_mismatch_parameters,
    #        generic_gens_para,
    #       selected_generic_gens_para,

    #       static_Idx_and_syms )
    
    
end




"""
    get_pf_streamedlined_simulation_parameters(
        net_data_by_components_file;
        <keyword arguments>)

It is used to simplify generation of parameters or data of a system that are needed for  power flow analyses.


# Arguments

- `net_data_by_components_file`: the network data file
- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,

`use_init_u0`: the boolean variable that determines if initial state u0 should be used in a power flow.
`use_nlsolve`: the boolean variable that determines if `nlsolve` should be used in power flow.

`pf_alg`: power flow solver

"""
function get_pf_streamedlined_simulation_parameters(
    net_data_by_components_file;
    components_libs_dir =
        nothing,
    basekV = 1.0,    
    use_pu_in_PQ      = true,
    line_data_in_pu   = true,
    with_faults       = false,
    pf_alg            = NewtonRaphson(),
    no_lines_fault    = 1)
    
    #--------------------------------------    
    
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

    (;plant_generators_data_from_json,
     plant_loads_data_from_json,
     plant_transmission_data_from_json,
     edge_data_from_json,
     shunt_data_from_json,
     baseMVA_data_from_json,
     gencost_data_from_json) =
        NamedTupleTools.select(
            get_net_data_by_components_from_json_file(
                net_data_by_components_file;
                in_components_type_sym =
                    false ),
            (:plant_generators_data_from_json,
             :plant_loads_data_from_json,
             :plant_transmission_data_from_json,
             :edge_data_from_json,
             :shunt_data_from_json,
             :baseMVA_data_from_json,
             :gencost_data_from_json))

    baseMVA = baseMVA_data_from_json

    #------------------------------------------------
    #------------------------------------------------
        
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )
    
    #------------------------------------------------
    
    gencost_data =
        get_gencost_data_by_json(
            gencost_data_from_json)

    gens_cost_coeff_ascen =
        get_gens_cost_coeff_in_ascen(
            gencost_data )
        
    
    gens_cost_coeff_decen =
        get_gens_cost_coeff_in_decen(
            gencost_data)

    #------------------------------------------------


    (;edges_orientation,
     edges_Ybr_cal,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json,
                 edge_data_from_json,
                 shunt_data_from_json;
                 baseMVA =
                     baseMVA,
                 basekV =
                     basekV,
                 use_pu_in_PQ =
                     use_pu_in_PQ,
                 line_data_in_pu =
                     line_data_in_pu ),
             (:edges_orientation,
              :edges_Ybr_cal,
              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))
    
    # Ynet_wt_nodes_idx_wt_adjacent_nodes = 
    #     get_Ynet(
    #         edge_data_from_json,
    #         shunt_data_from_json;
    #         baseMVA = baseMVA,
    #         basekV = basekV,
    #         baseShunt = baseMVA ,
    #         line_data_in_pu = line_data_in_pu)

    #------------------------------------------------
    
    Ybus =
        get_Ybus(
            edge_data_from_json,
            shunt_data_from_json;
            basekV = basekV,
            baseMVA = baseMVA,
            line_data_in_pu =
                line_data_in_pu )

    #------------------------------------------------
    
    (branches_fbus,
     branches_tbus,
     branches_r,
     branches_x,
     branches_b,
     branches_ratio,
     branches_angle,
     branches_type) =
     get_edges_ftbus_and_generic_data_by_json(
         edge_data_from_json )

    
    (edges_fbus, edges_tbus) =
         get_edges_fbus_tbus_by_json(
             edge_data_from_json)    
    
    #------------------------------------------------
    
    (;P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist ) =
         get_pf_PQ_param_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json;
             baseMVA =
                 baseMVA,
             use_pu_in_PQ =
                 use_pu_in_PQ )
        
    #------------------------------------------------
    #------------------------------------------------
    
    (Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         Ynet_wt_nodes_idx_wt_adjacent_nodes

    #------------------------------------------------
        
    no_edges = edges_size = length(edges_fbus)

    edges_r_x_b_ratio_angle_idx =
        get_edges_r_x_b_ratio_angle_idx(
            edges_size)

    #------------------------------------------------
    #------------------------------------------------
    
    ode_gens_para_selections  =
        (:H,
         :D,
         :X_d,
         :X_q,                  
         :X_d_dash,
         :X_q_dash,
         :T_d_dash,
         :T_q_dash, :Sn )

    ode_gens_para_sequence_order =
        (:components_data,
         :gen)
        
    ode_gens_generic_selections =
        (:H,
         :D,
         :ra,
         :xℓ,
         :X_d,
         :X_q,
         
         :X_d_dash,
         :X_q_dash,
         
         :X_d_2dash,
         :X_q_2dash,
         
         :T_d_dash,
         :T_q_dash,
         
         :Sn,
         
         :T_d_2dash,
         :T_q_2dash,
         
         :vh,
         :P,
         :Q,
         
         :Pmin,
         :Pmax,
         :Qmin,
         :Qmax,
         :vmin,
         :vmax )

    ode_gens_generic_sequence_order =
        (:components_data, :gen)
        
    govs_and_avrs_sequence_order =
        ( :components_data,)
    
    govs_and_avrs_selections =
        ( :gov, :avr )
    
    #----------------------------------------
    
    ode_gens_generic_para =
         get_ode_gens_generic_para(
             plant_generators_data_from_json;
             sequence_order =
                 ode_gens_generic_sequence_order,
             selections =
                 ode_gens_generic_selections)

    #------------------------------------------------

    
    # ode_gens_para =
    #     NamedTupleTools.select(
    #         ode_gens_generic_para,
    #         (:H,
    #          :D,
    #          :X_d,
    #          :X_q,
    #          :X_d_dash,
    #          :X_q_dash,
    #          :T_d_dash,
    #          :T_q_dash, :Sn))
    
    #------------------------------------------------

    (;generic_gens_para,) =
        NamedTupleTools.select(
            get_generic_gens_avr_gov_para(
                plant_generators_data_from_json;
                gens_sequence_order =
                    ode_gens_generic_sequence_order,
                gens_selections =
                    ode_gens_generic_selections,
                govs_and_avrs_sequence_order =
                    govs_and_avrs_sequence_order,
                govs_and_avrs_selections =
                    govs_and_avrs_selections),
                (:generic_gens_para,))

    #------------------------------------------------

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
                     :X_q_dash, :Sn ),
                vec_datatype = Float64 ),
            (:ra,
             :X_d,
             :X_q,     
             :X_d_dash,
             :X_q_dash, :Sn ) )

    
    #------------------------------------------------
    #------------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend ) =
         NamedTupleTools.select(
             get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                 getproperty(
                     Ynet_wt_nodes_idx_wt_adjacent_nodes,
                     :Ynet)),
             (:Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend) )

    #------------------------------------------------
    #------------------------------------------------
    
    dyn_pf_fun_kwd_net_idxs =
        NamedTupleTools.select(
            net_nodes_type_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx,
             :nodes_with_demands_idx))
    
    dyn_pf_fun_kwd_n2s_idxs =
        NamedTupleTools.select(
            get_dict_net_streamlined_idx_by_nodes_type_idxs(
                net_nodes_type_idxs ),
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx,
             :n2s_nodes_with_demands_idx))

    #------------------------------------------------
    
    (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx,
             :nodes_with_demands_idx))

    
   (;n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx,
    n2s_nodes_with_demands_idx) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx,
             :n2s_nodes_with_demands_idx))

    #------------------------------------------------
        
    no_nodes = nodes_size = length(all_nodes_idx)

    no_gens  = length(gens_nodes_idx )
    
    #------------------------------------------------    
    
    dyn_pf_mismatch_vars_kwd_para =
        (;Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para = pf_generic_gens_para  )

    sta_pf_PQ_para =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA =
                baseMVA,
            use_pu_in_PQ =
                use_pu_in_PQ)

    
    slack_gens_vh_θh_gens_vh_non_slack_gens_vh = gens_vh_slack_θh_para =
        get_gens_vh_slack_θh_para_by_json(
            plant_generators_data_from_json )

    
    sta_pf_vars_and_paras_idx =
        get_sta_pf_vars_and_paras_idx_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    
    pf_sta_ΔPQ_mismatch_parameters =
        get_pf_sta_ΔPQ_mismatch_parameters_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json;
            baseMVA = baseMVA,
            basekV = basekV,
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu)
    
    #----------------------------------------
    # Yred and Yint
    #----------------------------------------

    # X_d_dash =
    #     getproperty(
    #         ode_gens_para,
    #         :X_d_dash) 

    X_d_dash =
        getproperty(
            pf_generic_gens_para,
            :X_d_dash) 

    y_aug_kw_para =
        (;loc_load_exist,
         X_d_dash,    
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs)

    #----------------------------------------
    # 
    #----------------------------------------

    (gens_vh,
     sch_Pg, sch_Qg,
     gens_Pmin, gens_Pmax,
     gens_Qmin, gens_Qmax) =
         NamedTupleTools.select(
             opf_generic_gens_para,
             (:vh,
              :P, :Q,
              :Pmin, :Pmax,
              :Qmin, :Qmax) )


    
    if use_pu_in_PQ == true

        sch_Pg = sch_Pg ./baseMVA
        sch_Qg = sch_Qg ./baseMVA

        gens_Pmax = gens_Pmax ./baseMVA
        gens_Pmin = gens_Pmin ./baseMVA

        gens_Qmax = gens_Qmax ./baseMVA
        gens_Qmin = gens_Qmin ./baseMVA
        
    end
        
    gens_installed_capacity =
        copy(gens_Pmax)

    
    gens_Pmax =
        [(sch_Pg[idx] == 0.0 || sch_Pg[idx] == 0 ) ?
        sch_Pg[idx] : pg for (idx, pg) in
            enumerate(gens_Pmax)]
    
    #--------------------------------------

    nodes_Pd = loc_load_exist == true ?
        [ idx ∈ gens_nodes_with_loc_loads_idx ?
        P_g_loc_load[n2s_gens_with_loc_load_idxs[idx]] : 
        P_non_gens[n2s_non_gens_idx[idx]] 
        
                 for idx in nodes_with_demands_idx ] : [
        P_non_gens[n2s_non_gens_idx[idx]] 
                 for idx in nodes_with_demands_idx  ] 

    nodes_Qd = loc_load_exist == true ?
        [ idx ∈ gens_nodes_with_loc_loads_idx ?
        Q_g_loc_load[n2s_gens_with_loc_load_idxs[idx]] :
        Q_non_gens[n2s_non_gens_idx[idx]] 
        
                 for idx in nodes_with_demands_idx ] : [
        Q_non_gens[n2s_non_gens_idx[idx]] 
                 for idx in nodes_with_demands_idx  ] 

    if use_pu_in_PQ == false

        nodes_Pd = nodes_Pd .* baseMVA
        nodes_Qd = nodes_Qd .* baseMVA
        
    end
    

    #--------------------------------------
    #----------------------------------------    
    # Indicies
    #----------------------------------------


    Pg_Qg_Png_Qng_Pll_Qll_Idx = 
        get_generic_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)
    

    Png_Qng_Pll_Qll_Idx = 
        get_generic_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)

    
    Pg_Png_Qng_Idx = 
    get_generic_Pg_Png_Qng_Idx(
       dyn_pf_fun_kwd_net_idxs)

    
    pf_vh_θh_idx_and_idx2Idx =
        get_pf_vh_θh_idx_and_idx2Idx(
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs)
    
    
    dyn_pf_flat_vh_flat_θh_Idx =
        get_generic_flat_vh_flat_θh_Idx(
            all_nodes_idx)
    
    
    dyn_pf_flat_vh_flat_θh_id_iq_Idx =
        get_generic_flat_vh_flat_θh_id_iq_Idx(
            gens_nodes_idx,
            all_nodes_idx)

    
    dyn_pf_vh_θh_id_iq_vhf_θhf_Idx =
        get_generic_vh_θh_id_iq_vhf_θhf_Idx(
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault = no_lines_fault)

    
    dyn_pf_vh_vhf_θh_θhf_id_iq_Idx =
         get_generic_vh_vhf_θh_θhf_id_iq_Idx(
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault = no_lines_fault)

    
    dyn_pf_vh_vhf_Idx =
        get_generic_vh_vhf_Idx(
            all_nodes_idx;
            no_lines_fault = no_lines_fault)
    

    dyn_pf_θh_θhf_Idx =
        get_generic_θh_θhf_Idx(
            all_nodes_idx;
            no_lines_fault = no_lines_fault)
    

    id_iq_pg_vh_Idx = 
    get_id_iq_pg_vh_Idx(
        gens_nodes_idx)
    
    scale_Pg_Qg_Png_Qng_Pll_Qll_Idx = 
    get_generic_scale_Pg_Qg_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)

    scale_Pg_Png_Qng_Idx = 
    get_generic_scale_Pg_Png_Qng_Idx(
        dyn_pf_fun_kwd_net_idxs)

    
    dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx =
        get_generic_flat_vh_flat_θh_wt_slack_value_Idx(
            all_nodes_idx)
    
    #--------------------------------------

    (;dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))


    (;dyn_slack_value_Idxs,) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
            (:dyn_slack_value_Idxs, ))

    #----------------------------------------


    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]

    transformed_non_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_slack_gens_nodes_idx ]

    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]

    transformed_gens_with_loc_load_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_with_loc_load_idx ]
        
    # transformed_non_slack_gens_and_non_gens_idx = [
    #     n2s_all_nodes_idx[idx]
    #     for idx in non_slack_gens_and_non_gens_idx]
    
    transformed_nodes_with_demands_idx = [
        n2s_all_nodes_idx[idx]
        for idx in nodes_with_demands_idx ]
    
    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]


    vars_and_paras_Idx =
        (; Pg_Qg_Png_Qng_Pll_Qll_Idx,
         Png_Qng_Pll_Qll_Idx,
         Pg_Png_Qng_Idx,
         pf_vh_θh_idx_and_idx2Idx,

         dyn_pf_flat_vh_flat_θh_Idx,

         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
         dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

         dyn_pf_vh_vhf_Idx,
         dyn_pf_θh_θhf_Idx,

         id_iq_pg_vh_Idx,

         scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
         scale_Pg_Png_Qng_Idx,
         dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

         dyn_pf_vh_Idxs,
         dyn_pf_θh_Idxs,
         dyn_slack_value_Idxs,

         transformed_slack_gens_nodes_idx,
         transformed_non_slack_gens_nodes_idx,
         transformed_gens_nodes_idx,
         transformed_non_gens_nodes_idx,
         transformed_gens_with_loc_load_idx,
         # transformed_non_slack_gens_and_non_gens_idx,
         transformed_nodes_with_demands_idx,
         transformed_all_nodes_idx )    

    
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
    # Powerflow func and prob
    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (;Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_dyn_sol_kwd_para =
        (;loc_load_exist,
         sta_pf_PQ_para,
         ode_gens_generic_para,
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs)


    #--------------------------------------


    return (;baseMVA,
            basekV,
            vars_and_paras_Idx,
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

            generic_gens_para,
            ode_gens_generic_para,
            # ode_gens_para,

            gens_vh,
            gens_vh_slack_θh_para,
            
            slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

            pf_sta_ΔPQ_mismatch_parameters,
            sta_pf_PQ_para,
            pf_PQ_param,
            pf_generic_gens_para,
            
            gens_installed_capacity,

            gens_Pmax,
            gens_Pmin,
            gens_Qmax,
            gens_Qmin,

            sch_Pg,
            sch_Qg,
            nodes_Pd,
            nodes_Qd,

            branches_fbus,
            branches_tbus,
            branches_r,
            branches_x,
            branches_b,
            branches_ratio,
            branches_angle,
            branches_type,

            edges_orientation,

            Ybus,
            Ynet_wt_nodes_idx_wt_adjacent_nodes)

end





"""
    get_opf_streamedlined_simulation_parameters(
        net_data_by_components_file;
        <keyword arguments>)

It is used to simplify generation of parameters or data of a system that are needed for optimal power flow analyses.


# Arguments

- `net_data_by_components_file`: the network data file
- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,

`use_init_u0`: the boolean variable that determines if initial state u0 should be used in a power flow.
`use_nlsolve`: the boolean variable that determines if `nlsolve` should be used in power flow.

`pf_alg`: power flow solver

"""
function get_opf_streamedlined_simulation_parameters(
    net_data_by_components_file;
    components_libs_dir =
        nothing,
    basekV = 1.0,    
    use_pu_in_PQ      = true,
    opf_use_pu_in_PQ  = true,
    line_data_in_pu   = true,
    with_faults       = false )
    
    #--------------------------------------    
    
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

    (;plant_generators_data_from_json,
     plant_loads_data_from_json,
     plant_transmission_data_from_json,
     edge_data_from_json,
     shunt_data_from_json,
     baseMVA_data_from_json,
     gencost_data_from_json) =
        NamedTupleTools.select(
            get_net_data_by_components_from_json_file(
                net_data_by_components_file;
                in_components_type_sym =
                    false ),
            (:plant_generators_data_from_json,
             :plant_loads_data_from_json,
             :plant_transmission_data_from_json,
             :edge_data_from_json,
             :shunt_data_from_json,
             :baseMVA_data_from_json,
             :gencost_data_from_json))

    baseMVA = baseMVA_data_from_json

    #------------------------------------------------
    #------------------------------------------------
        
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )
    
    #------------------------------------------------
    
    gencost_data =
        get_gencost_data_by_json(
            gencost_data_from_json)

    gens_cost_coeff_ascen =
        get_gens_cost_coeff_in_ascen(
            gencost_data )
        
    
    gens_cost_coeff_decen =
        get_gens_cost_coeff_in_decen(
            gencost_data)

    #------------------------------------------------


    # (;edges_orientation,
    #  edges_Ybr_cal,
    #  Ybr_cal_and_edges_orientation,
    #  Ynet_wt_nodes_idx_wt_adjacent_nodes) =
    #      NamedTupleTools.select(
    #          get_transmission_network_parameters_by_json(
    #              plant_generators_data_from_json,
    #              plant_loads_data_from_json,
    #              plant_transmission_data_from_json,
    #              edge_data_from_json,
    #              shunt_data_from_json;
    #              baseMVA =
    #                  baseMVA,
    #              basekV =
    #                  basekV,
    #              use_pu_in_PQ =
    #                  use_pu_in_PQ,
    #              line_data_in_pu =
    #                  line_data_in_pu ),
    #          (:edges_orientation,
    #           :edges_Ybr_cal,
    #           :Ybr_cal_and_edges_orientation,
    #           :Ynet_wt_nodes_idx_wt_adjacent_nodes))
    
    Ynet_wt_nodes_idx_wt_adjacent_nodes = 
        get_Ynet(
            edge_data_from_json,
            shunt_data_from_json;
            baseMVA = baseMVA,
            basekV = basekV,
            baseShunt = baseMVA ,
            line_data_in_pu = line_data_in_pu)

    #------------------------------------------------
    
    Ybus =
        get_Ybus(
            edge_data_from_json,
            shunt_data_from_json;
            basekV = basekV,
            baseMVA = baseMVA,
                line_data_in_pu = line_data_in_pu )

    #------------------------------------------------
    
    (branches_fbus,
     branches_tbus,
     branches_r,
     branches_x,
     branches_b,
     branches_ratio,
     branches_angle,
     branches_type) =
     get_edges_ftbus_and_generic_data_by_json(
         edge_data_from_json )

    
    (edges_fbus, edges_tbus) =
         get_edges_fbus_tbus_by_json(
             edge_data_from_json)    
    
    #------------------------------------------------
    
    (;P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist ) =
         get_pf_PQ_param_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json;
             baseMVA =
                 baseMVA,
             use_pu_in_PQ =
                 use_pu_in_PQ )
        
    #------------------------------------------------
    #------------------------------------------------
    
    (Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         Ynet_wt_nodes_idx_wt_adjacent_nodes

    #------------------------------------------------
        
    no_edges = edges_size = length(edges_fbus)

    edges_r_x_b_ratio_angle_idx =
        get_edges_r_x_b_ratio_angle_idx(
            edges_size)

    #------------------------------------------------

    ode_gens_para_selections  =
        (:vh,
         :P, :Q,
         :Pmin, :Pmax,
         :Qmin, :Qmax,
         :vmin, :vmax,
         :Sn   )

    ode_gens_para_sequence_order =
        (:components_data,
         :gen)
    
    #----------------------------------------
    
    opf_generic_gens_para =
         get_ode_gens_generic_para(
             plant_generators_data_from_json;
             sequence_order =
                 ode_gens_para_sequence_order,
             selections =
                 ode_gens_para_selections )
    
    #------------------------------------------------
    #------------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend ) =
         NamedTupleTools.select(
             get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                 getproperty(
                     Ynet_wt_nodes_idx_wt_adjacent_nodes,
                     :Ynet)),
             (:Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend) )

    
    #----------------------------------------    
    #----------------------------------------    
    
    dyn_pf_fun_kwd_net_idxs =
        NamedTupleTools.select(
            net_nodes_type_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx,
             :nodes_with_demands_idx))
    
    dyn_pf_fun_kwd_n2s_idxs =
        NamedTupleTools.select(
            get_dict_net_streamlined_idx_by_nodes_type_idxs(
                net_nodes_type_idxs ),
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx,
             :n2s_nodes_with_demands_idx))

    #------------------------------------------------

    
    (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx,
             :nodes_with_demands_idx))

    
   (;n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx,
    n2s_nodes_with_demands_idx) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx,
             :n2s_nodes_with_demands_idx))
    
    #------------------------------------------------
    
    no_nodes = nodes_size = length(all_nodes_idx)

    no_gens  = length(gens_nodes_idx )

    #------------------------------------------------


    (gens_vh,
     sch_Pg, sch_Qg,
     gens_Pmin, gens_Pmax,
     gens_Qmin, gens_Qmax) =
         NamedTupleTools.select(
             opf_generic_gens_para,
             (:vh,
              :P, :Q,
              :Pmin, :Pmax,
              :Qmin, :Qmax) )


    
    if opf_use_pu_in_PQ == true

        sch_Pg = sch_Pg ./baseMVA
        sch_Qg = sch_Qg ./baseMVA

        gens_Pmax = gens_Pmax ./baseMVA
        gens_Pmin = gens_Pmin ./baseMVA

        gens_Qmax = gens_Qmax ./baseMVA
        gens_Qmin = gens_Qmin ./baseMVA
        
    end
        
    gens_installed_capacity =
        copy(gens_Pmax)

    
    gens_Pmax =
        [(sch_Pg[idx] == 0.0 || sch_Pg[idx] == 0 ) ?
        sch_Pg[idx] : pg for (idx, pg) in
            enumerate(gens_Pmax)]
    
    #--------------------------------------

    nodes_Pd = loc_load_exist == true ?
        [ idx ∈ gens_nodes_with_loc_loads_idx ?
        P_g_loc_load[n2s_gens_with_loc_load_idxs[idx]] : 
        P_non_gens[n2s_non_gens_idx[idx]] 
        
                 for idx in nodes_with_demands_idx ] : [
        P_non_gens[n2s_non_gens_idx[idx]] 
                 for idx in nodes_with_demands_idx  ] 

    nodes_Qd = loc_load_exist == true ?
        [ idx ∈ gens_nodes_with_loc_loads_idx ?
        Q_g_loc_load[n2s_gens_with_loc_load_idxs[idx]] :
        Q_non_gens[n2s_non_gens_idx[idx]] 
        
                 for idx in nodes_with_demands_idx ] : [
        Q_non_gens[n2s_non_gens_idx[idx]] 
                 for idx in nodes_with_demands_idx  ] 

    if opf_use_pu_in_PQ == false

        nodes_Pd = nodes_Pd .* baseMVA
        nodes_Qd = nodes_Qd .* baseMVA
        
    end
    
   #--------------------------------------


    P_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmin, no_nodes)

    P_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmax, no_nodes)

    #--------------------------------------


    Q_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmin, no_nodes)

    Q_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmax, no_nodes)
    
    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Pd , no_nodes )

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Qd , no_nodes)

    S_Demand = P_Demand + im * Q_Demand

    #--------------------------------------
    #--------------------------------------


    return (;baseMVA,
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

            gens_vh,
            # ds_generic_gens_para,
            gens_vh_slack_θh_para,

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
            sch_Qg,
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

            Ybus,
            Ynet_wt_nodes_idx_wt_adjacent_nodes)


end




"""
    get_generic_system_simulation_parameters(
        net_data_by_components_file;
        <keyword arguments>)

It is used to simplify generation of parameters or data of a system that are needed for simulation.


# Arguments

- `net_data_by_components_file`: the network data file
- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,

`use_init_u0`: the boolean variable that determines if initial state u0 should be used in a power flow.
`use_nlsolve`: the boolean variable that determines if `nlsolve` should be used in power flow.

`pf_alg`: power flow solver
`ode_alg`: ode solver`
`dae_alg`: dae solver

"""
function get_generic_system_simulation_parameters(
    net_data_by_components_file;
    components_libs_dir =
        nothing,
    basekV = 1.0,    
    use_pu_in_PQ      = true,
    opf_use_pu_in_PQ  = true,
    line_data_in_pu   = true,
    with_faults       = false,
    pf_alg            = NewtonRaphson(),
    ode_alg           = Rodas4(),
    dae_alg           = IDA() )
    
    #--------------------------------------    
    
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

    (;plant_generators_data_from_json,
     plant_loads_data_from_json,
     plant_transmission_data_from_json,
     edge_data_from_json,
     shunt_data_from_json,
     baseMVA_data_from_json,
     gencost_data_from_json) =
        NamedTupleTools.select(
            get_net_data_by_components_from_json_file(
                net_data_by_components_file;
                in_components_type_sym =
                    false ),
            (:plant_generators_data_from_json,
             :plant_loads_data_from_json,
             :plant_transmission_data_from_json,
             :edge_data_from_json,
             :shunt_data_from_json,
             :baseMVA_data_from_json,
             :gencost_data_from_json))

    baseMVA = baseMVA_data_from_json
    
    #------------------------------------------------
    
    (;P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist ) =
         get_pf_PQ_param_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json;
             baseMVA =
                 baseMVA,
             use_pu_in_PQ =
                 use_pu_in_PQ )
    
    #----------------------------------------

    ode_gens_para_selections  =
        (:H,
         :D,
         :X_d,
         :X_q,                  
         :X_d_dash,
         :X_q_dash,
         :T_d_dash,
         :T_q_dash, :Sn )

    ode_gens_para_sequence_order =
        (:components_data,
         :gen)
        
    ode_gens_generic_selections =
        (:H,
         :D,
         :ra,
         :xℓ,
         :X_d,
         :X_q,
         
         :X_d_dash,
         :X_q_dash,
         
         :X_d_2dash,
         :X_q_2dash,
         
         :T_d_dash,
         :T_q_dash,
         
         :Sn,
         
         :T_d_2dash,
         :T_q_2dash,
         
         :vh,
         :P,
         :Q,
         
         :Pmin,
         :Pmax,
         :Qmin,
         :Qmax,
         :vmin,
         :vmax )

    ode_gens_generic_sequence_order =
        (:components_data, :gen)
        
    govs_and_avrs_sequence_order =
        ( :components_data,)
    
    govs_and_avrs_selections =
        ( :gov, :avr )
    
    #----------------------------------------
    
    ode_gens_generic_para =
         get_ode_gens_generic_para(
             plant_generators_data_from_json;
             sequence_order =
                 ode_gens_generic_sequence_order,
             selections =
                 ode_gens_generic_selections)
    
    #------------------------------------------------

    (;generic_gens_para,
     generic_govs_para,
     generic_avrs_para) =
         get_generic_gens_avr_gov_para(
             plant_generators_data_from_json;
             gens_sequence_order =
                 ode_gens_generic_sequence_order,
            gens_selections =
                ode_gens_generic_selections,
            govs_and_avrs_sequence_order =
                govs_and_avrs_sequence_order,
            govs_and_avrs_selections =
                govs_and_avrs_selections)

    #------------------------------------------------

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
                     :X_q_dash, :Sn ),
                vec_datatype = Float64 ),
            (:ra,
             :X_d,
             :X_q,     
             :X_d_dash,
             :X_q_dash, :Sn ) )

    #------------------------------------------------

     opf_generic_gens_para =
        NamedTupleTools.select(
            get_selected_vec_nt_to_vec_vec(
                generic_gens_para,
                nothing;
                selections =
                    (:vh,
                     :P, :Q,
                     :Pmin, :Pmax,
                     :Qmin, :Qmax,
                     :vmin, :vmax,
                     :Sn ),
                vec_datatype = Float64 ),
            (:vh,
             :P, :Q,
             :Pmin, :Pmax,
             :Qmin, :Qmax,
             :vmin, :vmax,
             :Sn ) )

    #------------------------------------------------

     ds_generic_gens_para =
        NamedTupleTools.select(
            get_selected_vec_nt_to_vec_vec(
                generic_gens_para,
                nothing;
                selections =
                    (:vh,
                     :P, :Q,
                     :Pmin, :Pmax,
                     :Qmin, :Qmax,
                     :vmin, :vmax,
                     :Sn ),
                vec_datatype = Float64 ),
            (:vh,
             :P, :Q,
             :Pmin, :Pmax,
             :Qmin, :Qmax,
             :vmin, :vmax,
             :Sn ) )

    #------------------------------------------------    
    
    (branches_fbus,
     branches_tbus,
     branches_r,
     branches_x,
     branches_b,
     branches_ratio,
     branches_angle,
     branches_type) =
          get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    (branches_fbus,
     branches_tbus) =
         get_edges_fbus_tbus_by_json(
             edge_data_from_json)

    (edges_fbus, edges_tbus,
     edges_r, edges_x, edges_b,
     edges_ratio, edges_angle, edges_type) =
         get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )
    
    (shunt_idx, Gs, Bs) =
        get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
            shunt_data_from_json)
    
    #------------------------------------------------
    
    gencost_data =
        get_gencost_data_by_json(
            gencost_data_from_json)

    gens_cost_coeff_ascen =
        get_gens_cost_coeff_in_ascen(
            gencost_data )
        
    
    gens_cost_coeff_decen =
        get_gens_cost_coeff_in_decen(
            gencost_data)

    #------------------------------------------------

    (;edges_orientation,
     edges_Ybr_cal,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json,
                 edge_data_from_json,
                 shunt_data_from_json;
                 baseMVA =
                     baseMVA,
                 basekV =
                     basekV,
                 use_pu_in_PQ =
                     use_pu_in_PQ,
                 line_data_in_pu =
                     line_data_in_pu ),
             (:edges_orientation,
              :edges_Ybr_cal,
              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))
    
    #------------------------------------------------
    #------------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend ) =
         NamedTupleTools.select(
             get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                 getproperty(
                     Ynet_wt_nodes_idx_wt_adjacent_nodes,
                     :Ynet)),
             (:Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend) )

    #------------------------------------------------
    #------------------------------------------------
    
    ode_gens_para =
        NamedTupleTools.select(
            ode_gens_generic_para,
            (:H,
             :D,
             :X_d,
             :X_q,
             :X_d_dash,
             :X_q_dash,
             :T_d_dash,
             :T_q_dash, :Sn))
    
    #------------------------------------------------
    #------------------------------------------------
        
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )
    
    dyn_pf_fun_kwd_net_idxs =
        NamedTupleTools.select(
            net_nodes_type_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx,
             :nodes_with_demands_idx))
    
    dyn_pf_fun_kwd_n2s_idxs =
        NamedTupleTools.select(
            get_dict_net_streamlined_idx_by_nodes_type_idxs(
                net_nodes_type_idxs ),
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx,
             :n2s_nodes_with_demands_idx))

    # @show tnodes_with_demands_idx
    
    #------------------------------------------------
    #------------------------------------------------
    
    dyn_pf_mismatch_vars_kwd_para =
        (;Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para )

    sta_pf_PQ_para =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA =
                baseMVA,
            use_pu_in_PQ =
                use_pu_in_PQ)

    
    slack_gens_vh_θh_gens_vh_non_slack_gens_vh = gens_vh_slack_θh_para =
        get_gens_vh_slack_θh_para_by_json(
            plant_generators_data_from_json )

    
    sta_pf_vars_and_paras_idx =
        get_sta_pf_vars_and_paras_idx_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    
    pf_sta_ΔPQ_mismatch_parameters =
        get_pf_sta_ΔPQ_mismatch_parameters_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json;
            baseMVA = baseMVA,
            basekV = basekV,
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu)

    #----------------------------------------
    #----------------------------------------

    (slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,
     nodes_with_demands_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,
              :nodes_with_demands_idx))


   (;n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx,
    n2s_nodes_with_demands_idx) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx,
             :n2s_nodes_with_demands_idx))
    
    #----------------------------------------
    # Yred and Yint
    #----------------------------------------

    X_d_dash =
        getproperty(
            ode_gens_para,
            :X_d_dash) 

    y_aug_kw_para =
        (;loc_load_exist,
         X_d_dash,    
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs)

    #----------------------------------------
    # 
    #----------------------------------------
    
    (Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         Ynet_wt_nodes_idx_wt_adjacent_nodes
    
    edges_orientation =
        getproperty(
            Ybr_cal_and_edges_orientation,
            :edges_orientation)

    #--------------------------------------    
        
    no_edges = edges_size = length(edges_r)

    edges_r_x_b_ratio_angle_idx =
        get_edges_r_x_b_ratio_angle_idx(
            edges_size)
    
    #--------------------------------------
    
    no_nodes = nodes_size = length(all_nodes_idx)

    no_gens  = length(gens_nodes_idx )


    (gens_vh,
     sch_Pg, sch_Qg,
     gens_Pmin, gens_Pmax,
     gens_Qmin, gens_Qmax) =
         NamedTupleTools.select(
             opf_generic_gens_para,
             (:vh,
              :P, :Q,
              :Pmin, :Pmax,
              :Qmin, :Qmax) )


    
    if opf_use_pu_in_PQ == true

        sch_Pg = sch_Pg ./baseMVA
        sch_Qg = sch_Qg ./baseMVA

        gens_Pmax = gens_Pmax ./baseMVA
        gens_Pmin = gens_Pmin ./baseMVA

        gens_Qmax = gens_Qmax ./baseMVA
        gens_Qmin = gens_Qmin ./baseMVA
        
    end
        
    gens_installed_capacity =
        copy(gens_Pmax)

    
    gens_Pmax =
        [(sch_Pg[idx] == 0.0 || sch_Pg[idx] == 0 ) ?
        sch_Pg[idx] : pg for (idx, pg) in
            enumerate(gens_Pmax)]
    
    #--------------------------------------

    nodes_Pd = loc_load_exist == true ?
        [ idx ∈ gens_nodes_with_loc_loads_idx ?
        P_g_loc_load[n2s_gens_with_loc_load_idxs[idx]] : 
        P_non_gens[n2s_non_gens_idx[idx]] 
        
                 for idx in nodes_with_demands_idx ] : [
        P_non_gens[n2s_non_gens_idx[idx]] 
                 for idx in nodes_with_demands_idx  ] 

    nodes_Qd = loc_load_exist == true ?
        [ idx ∈ gens_nodes_with_loc_loads_idx ?
        Q_g_loc_load[n2s_gens_with_loc_load_idxs[idx]] :
        Q_non_gens[n2s_non_gens_idx[idx]] 
        
                 for idx in nodes_with_demands_idx ] : [
        Q_non_gens[n2s_non_gens_idx[idx]] 
                 for idx in nodes_with_demands_idx  ] 

    if opf_use_pu_in_PQ == false

        nodes_Pd = nodes_Pd .* baseMVA
        nodes_Qd = nodes_Qd .* baseMVA
        
    end
    
   #--------------------------------------


    P_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmin, no_nodes)

    P_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Pmax, no_nodes)

    #--------------------------------------


    Q_Gen_lb =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmin, no_nodes)

    Q_Gen_ub =
        SparseArrays.sparsevec(
            gens_nodes_idx, gens_Qmax, no_nodes)

    # Power demand levels
    # (real, reactive, and complex form)

    # @show length(nodes_with_demands_idx), length(nodes_Pd)
    
    P_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Pd , no_nodes )

    Q_Demand =
        SparseArrays.sparsevec(
            nodes_with_demands_idx, nodes_Qd , no_nodes)

    S_Demand = P_Demand + im * Q_Demand

    #--------------------------------------

    incidence_matrix  =
        SparseArrays.sparse(
            edges_fbus,
            1:no_edges, 1, no_nodes, no_edges) +

                SparseArrays.sparse(
                    edges_tbus,
                    1:no_edges, -1, no_nodes, no_edges)

    #------------------------------------------------
    
    ys =
        1.0 ./ (edges_r + im * edges_x) 
       

    y_c =
        1 / 2 * (im *  edges_b)

    y_sh =  (Gs .+ im *  Bs )

    if line_data_in_pu == true

        baseZ = basekV^2/baseMVA

        baseY = 1/baseZ


        Gs = Gs ./  baseY
        Bs = Bs ./  baseY        

        ys   = baseY ./ (edges_r + im * edges_x) 
        
        y_c  = 1 / 2 * (im *  edges_b) * baseY 

        y_sh =  (Gs .+ im *  Bs )

        
    end

    inv_τ =
        [ (a_ratio == 0.0 || a_ratio == 0) ? 1.0 : 1/a_ratio
              for a_ratio in
                  edges_ratio  ]
    
    θ_shift = edges_angle
    
    Yff = (ys +  y_c) .* (inv_τ).^2

    Ytf = -ys .* (inv_τ ./ exp.(im * θ_shift))

    Yft = -ys .* (inv_τ ./ exp.(-im * θ_shift))

    Ytt = (ys +  y_c)

        
    Y_sh = SparseArrays.spdiagm(y_sh)

    
    Cf = SparseArrays.sparse(
        1:no_edges,
        edges_fbus,
             1,  no_edges, no_nodes) 

    Ct = SparseArrays.sparse(
        1:no_edges,
        edges_tbus,
             1,  no_edges, no_nodes) 
    
    Yf = SparseArrays.spdiagm( Yff ) * Cf +
        SparseArrays.spdiagm( Yft )  * Ct

    Yt = SparseArrays.spdiagm( Ytf ) * Cf +
        SparseArrays.spdiagm( Ytt )  * Ct

    Ybus = Cf' * Yf + Ct' * Yt + Y_sh
    
    #----------------------------------------    
    
    states_Idx_syms_wt_functions =
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

              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :model_mass_matrix,     
              :model_bool_dae_vars,

              :ode_gens_mass_matrix,
              :ode_gens_bool_dae_vars,

              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,

              :dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :dyn_vh_id_iq_V_ref_Tm_Idx,
              :dyn_V_ref_Tm_id_iq_vh_Idx,

              :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
              :dyn_vh_id_iq_ωref0_vref0_porder0_Idx,

              :ωref0_vref0_porder0_id_iq_vh_Idx,

              :dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              :pf_vh_θh_idx_and_idx2Idx,
              
              :dyn_pf_flat_vh_flat_θh_Idx,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :system_states_idx_kwd_para,
              :system_paras_idx_kwd_para,
              :plants_dyn_fun_idx_kwd_para,
              :plants_algeb_fun_idx_kwd_para,

              :id_iq_pg_vh_Idx,

              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :scale_Pg_Png_Qng_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx))
    

    (;state_vars_idx,
     plants_states_syms,
     comps_callback_paras_funs,
     comps_init_funs,
     dyn_pf_flat_vh_flat_θh_id_iq_Idx) =
         NamedTupleTools.select(
             states_Idx_syms_wt_functions,
             (:state_vars_idx,
              :plants_states_syms,
              :comps_callback_paras_funs,
              :comps_init_funs,
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx))

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
    # Powerflow func and prob
    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (;Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_dyn_sol_kwd_para =
        (;loc_load_exist,
         sta_pf_PQ_para,
         ode_gens_generic_para,
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs)

    #----------------------------------------
    # callbacks
    #----------------------------------------

    plants_cb_paras_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para,
         comps_callback_paras_funs )

    
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
    
    cb_states = cb_fun_make_state_callbacks(
        list_selected_plants_state_event_cb_paras,
        list_selected_plants_state_affect_cb_paras )

    #----------------------------------------
    # Init
    #----------------------------------------
    
    plants_init_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para ,
         comps_init_funs )


    #--------------------------------------
    

    dynamic_parameters =
        (;states_Idx_syms_wt_functions,
         
         plants_cb_paras_switches,
         cb_states,
         
         baseMVA,
         basekV,

         gens_vh,
         
         ds_generic_gens_para,
         gens_vh_slack_θh_para,

         loc_load_exist,
         # state_vars_idx,

         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx,            
         nodes_with_demands_idx,

         # state_labels,
         # algebraic_vars_labels,
         # network_vars_labels,

         # plants_states_syms,
         # generic_state_sym, 
         # algebraic_state_sym,            
         # model_syms,

         # nodes_names,
         # gens_nodes_names,
         # non_gens_nodes_names,
         # SM_gens_nodes_names,
         # SC_gens_nodes_names,

         # model_mass_matrix,     
         # model_bool_dae_vars,

         # ode_gens_mass_matrix,
         # ode_gens_bool_dae_vars,

         # gens_state_vars_idx_in_state,
         # state_vars_and_i_dq_Idx_in_state,

         # state_algebraic_vars_Idx_in_state,

         # pf_vh_θh_idx_and_idx2Idx,

         # dyn_pf_flat_vh_flat_θh_Idx,

         # dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         # dyn_pf_θh_θhf_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         # dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
         # dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

         # ωref0_vref0_porder0_id_iq_vh_Idx,
         # dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

         # vec_comp_states_Idx,

         # comps_callback_paras_funs,     
         # comps_output_funs,
         # comps_dyn_funs,

         # ode_comps_dyn_funs,
         # dae_comps_dyn_funs,

         # comps_init_funs,

         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         ode_gens_para,

         sta_pf_PQ_para,
         pf_PQ_param,
         pf_generic_gens_para,

         edges_fbus,
         edges_tbus,
         edges_r,
         edges_x,
         edges_b,
         edges_ratio,
         edges_angle,
         edges_type,
         shunt_idx,
         Gs,
         Bs,

         Ynet_wt_nodes_idx_wt_adjacent_nodes,
         Ybr_cal_and_edges_orientation,

         kwd_sta_sta_ΔPQ_sol_by_json,
         generic_red_sol_kwd_para,
         generic_dyn_sol_kwd_para,

         plants_init_kwd_para,
         plants_cb_paras_kwd_para,

         y_aug_kw_para,

         with_faults,

         # Pg_Qg_Png_Qng_Pll_Qll_Idx,
         # Png_Qng_Pll_Qll_Idx,
         # Pg_Png_Qng_Idx,

         edges_r_x_b_ratio_angle_idx,

         Ynet_rows_Idxs_in_flattend,
         Ynet_real_imag_Idxs_in_flattend,

         # id_iq_pg_vh_Idx,

         # scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
         # scale_Pg_Png_Qng_Idx,
         # dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx
         )

    pf_parameters =
        (;baseMVA,
         basekV,
         loc_load_exist,
         
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx,
            
         nodes_with_demands_idx,
            
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_n2s_idxs,

         generic_gens_para,
         ode_gens_para,

         gens_vh,
         
         ds_generic_gens_para,
         gens_vh_slack_θh_para,
         slack_gens_vh_θh_gens_vh_non_slack_gens_vh,
         
         pf_sta_ΔPQ_mismatch_parameters,
         sta_pf_PQ_para,
         pf_PQ_param,
         pf_generic_gens_para,
         opf_generic_gens_para,

         gens_installed_capacity,

         gens_Pmax,
         gens_Pmin,
         gens_Qmax,
         gens_Qmin,

         sch_Pg,
         sch_Qg,
         nodes_Pd,
         nodes_Qd,

         ) 


    opf_parameters =
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

         gens_vh,
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
         sch_Qg,
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
         Ybus )

    return (;dynamic_parameters,
            pf_parameters,
            opf_parameters,
            states_Idx_syms_wt_functions)

end

#####################################################
# ---------------------------------------------------
# General Data aggregation utility functions
# ---------------------------------------------------
#####################################################


"""
    get_system_simulation_parameters(
        net_data_by_components_file,
        nothing;
        <keyword arguments>)

It is used to simplify generation of parameters or data of a system that are needed for simulation.


# Arguments

- `net_data_by_components_file`: the network data file
- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,

`use_init_u0`: the boolean variable that determines if initial state u0 should be used in a power flow.
`use_nlsolve`: the boolean variable that determines if `nlsolve` should be used in power flow.

`pf_alg`: power flow solver

`abstol`: the absolute error tolerance
`reltol`: the relative error tolerance

"""
function get_system_simulation_parameters(
    net_data_by_components_file,
    nothing;
    components_libs_dir = nothing,
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    
    use_nlsolve = false,
    
    pf_alg      = NewtonRaphson(),
    
    abstol      = 1e-12,
    
    reltol      = 1e-12 ) 

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
                        false),
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

              :dyn_pf_mismatch_vars_kwd_para ) )
    
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
    
    (;state_vars_idx,
     vec_comp_states_Idx,

     plants_states_syms,
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
     Pg_Png_Qng_Idx,

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
              :Pg_Png_Qng_Idx,

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
    #----------------------------------------
    # Dynamic simulation
    #----------------------------------------
    #----------------------------------------

    #----------------------------------------
    # Init
    #----------------------------------------
    
    plants_init_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para ,
         comps_init_funs )
    
    plants_generic_states_init_wt_ref =
        plants_generic_model_init_func(
            gens_vh,
            gens_θh,
            pf_P_gens,
            pf_Q_gens,
            ωs;
            kwd_para =
                plants_init_kwd_para )

    (plants_states_init,
     plants_refs ) =
         NamedTupleTools.select(
             plants_generic_states_init_wt_ref,
             (:plants_states_init,
              :plants_refs))


    ( nt_vec_per_paras,
      vec_vec_per_paras ) =
          get_nt_vec_wt_vec_vec_per_paras(
        plants_refs ;
        nt_syms =
            (:ωs,
             :ω_ref,
             :v_ref,
             :p_order,
             :i_d,
             :i_q,
             :gen_vh ) )

    (ω_ref,
     v_ref,
     p_order,
     gens_i_d,
     gens_i_q ) =
         NamedTupleTools.select(
             nt_vec_per_paras, (
                 :ω_ref,
                 :v_ref,
                 :p_order,
                 :i_d,
                 :i_q))

    #----------------------------------------
    # System model init
    #----------------------------------------

    u0_model_states_init =
        Float64[plants_states_init;
                vh;
                θh;
                gens_i_d;
                gens_i_q]

    
    #----------------------------------------
    # callbacks
    #----------------------------------------

    plants_cb_paras_kwd_para =
        (generic_gens_para,
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
    
    cb_states = cb_fun_make_state_callbacks(
        list_selected_plants_state_event_cb_paras,
        list_selected_plants_state_affect_cb_paras )
    
    #----------------------------------------
    # System model para
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
    # model para with callback
    #----------------------------------------

    """

    generic_model_dynamics_para =
        Float64[ω_ref;
                v_ref;
                p_order;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]

    model_dynamics_para =
        (;generic_model_dynamics_para,
          plants_cb_paras_switches)
    """
    
    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh]
    
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll = 
        Float64[ω_ref;
                v_ref;
                p_order;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]

    #----------------------------------------    
    # algebraic equation kwd_para
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
                     :X_q_dash, :Sn ),
                vec_datatype = Float64 ),
            (:ra,
             :X_d,
             :X_q,     
             :X_d_dash,
             :X_q_dash, :Sn ) )

    flat_vh_flat_θh_id_iq_u0 =
        [vh;
         θh;
         gens_i_d;
         gens_i_q ]
    
    dyn_pf_solver =
        (;use_init_u0,
         use_nlsolve,
         pf_alg,
         abstol,
         reltol )

    algebraic_generic_model_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes )
    
    algebraic_generic_model_sol_kwd_para =
        (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         flat_vh_flat_θh_id_iq_u0,
         dyn_pf_solver,
         algebraic_generic_model_kwd_para )

    #----------------------------------------
    #----------------------------------------
    # ODE system dyamanics simulation
    #----------------------------------------    
    #----------------------------------------
    
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

    
    dae_plants_kwd_para =
        (;state_vars_idx,

         ωref0_vref0_porder0_id_iq_vh_Idx,
         
         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         vec_comp_states_Idx,
         
         avrs_govs_cb_sw,
         avrs_govs_cb_sw_Idx,
         
         dae_comps_dyn_funs,
         comps_output_funs,
         ωs)


    plants_kwd_para =
        (;state_vars_idx,

         ωref0_vref0_porder0_id_iq_vh_Idx,

         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         vec_comp_states_Idx,

         avrs_govs_cb_sw,
         avrs_govs_cb_sw_Idx,

         comps_dyn_funs,
         comps_output_funs,

        ωs) 
    
    generic_system_dynamics_kwd_para =
        (;
         gens_nodes_idx,
         ωs,
         loc_load_exist,
         state_vars_idx,

         id_iq_pg_vh_Idx,

         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

         gens_state_vars_idx_in_state,
         state_vars_and_i_dq_Idx_in_state,
         state_algebraic_vars_Idx_in_state,

         dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         dae_plants_kwd_para,
         ode_plants_kwd_para,

         plants_kwd_para,

         algebraic_generic_model_kwd_para,

         algebraic_generic_model_sol_kwd_para,

         state_labels,
         algebraic_vars_labels,
         network_vars_labels )
    

    #----------------------------------------
    #----------------------------------------

    return (;u0_model_states_init,
            model_syms,
            model_mass_matrix,
            model_bool_dae_vars,
            
            algebraic_state_sym,
            ode_gens_mass_matrix,
            ode_gens_bool_dae_vars,

            nodes_names,
            gens_nodes_names,
            non_gens_nodes_names,

            SM_gens_nodes_names,
            SC_gens_nodes_names,

            plants_states_init,

            plants_cb_paras_switches,
            avrs_govs_cb_sw,
            avrs_govs_cb_sw_Idx,
            cb_states,

            ωref0_vref0_porder0_id_iq_vh,
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

            algebraic_generic_model_sol_kwd_para,
            generic_system_dynamics_kwd_para,
            
            pf_P_gens,
            pf_Q_gens,
            vh,
            θh,
            gens_vh,
            gens_θh,
            
            ω_ref,
            v_ref,
            p_order,
            gens_i_d,
            gens_i_q,

            pf_generic_gens_para,
            flat_vh_flat_θh_id_iq_u0,
            
            generic_gens_para,
            generic_avrs_para,
            generic_govs_para,

            # Diag
            
            algebraic_generic_model_kwd_para,
            dae_plants_kwd_para,
            ode_plants_kwd_para,
            plants_kwd_para,

            plants_states_syms,
            gens_nodes_idx,
            
            state_labels,
            algebraic_vars_labels,
            network_vars_labels,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs,
            
            Ybr_cal_and_edges_orientation,
            Ynet_wt_nodes_idx_wt_adjacent_nodes)
    
end


"""
    get_system_simulation_parameters(
        net_data_by_components_file;
        <keyword arguments>)

It is used to simplify generation of parameters or data of a system that are needed for simulation.


# Arguments

- `net_data_by_components_file`: the network data file
- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,

`use_init_u0`: the boolean variable that determines if initial state u0 should be used in a power flow.
`use_nlsolve`: the boolean variable that determines if `nlsolve` should be used in power flow.

`pf_alg`: power flow solver

`abstol`: the absolute error tolerance
`reltol`: the relative error tolerance

"""
function get_system_simulation_parameters(
    net_data_by_components_file;
    components_libs_dir = nothing,
    basekV          = 1.0,    
    use_pu_in_PQ    = true,
    line_data_in_pu = true,

    use_init_u0 = false,    
    use_nlsolve = false,    
    pf_alg      = NewtonRaphson(),    
    abstol      = 1e-12,
    reltol      = 1e-12 ) 

    #--------------------------------------
    
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

    (;plant_generators_data_from_json,
     plant_loads_data_from_json,
     plant_transmission_data_from_json,
     edge_data_from_json,
     shunt_data_from_json,
     baseMVA_data_from_json) =
        NamedTupleTools.select(
            get_net_data_by_components_from_json_file(
                net_data_by_components_file;
                in_components_type_sym =
                    false ),
            (:plant_generators_data_from_json,
             :plant_loads_data_from_json,
             :plant_transmission_data_from_json,
             :edge_data_from_json,
             :shunt_data_from_json,
             :baseMVA_data_from_json))

    #----------------------------------------

    baseMVA = baseMVA_data_from_json
    
    #------------------------------------------------
    
    (;
     P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist ) =
         get_pf_PQ_param_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json;
             baseMVA = baseMVA,
             use_pu_in_PQ = use_pu_in_PQ )
    
    #----------------------------------------
    #----------------------------------------

    ode_gens_para_selections  =
        (:H, :D,
         :X_d, :X_q,                  
         :X_d_dash, :X_q_dash,
         :T_d_dash, :T_q_dash, :Sn )

    ode_gens_para_sequence_order =
        (:components_data, :gen)
        
    ode_gens_generic_selections =
        (:H, :D,
         :ra, :xℓ,
         :X_d, :X_q,
         :X_d_dash,  :X_q_dash,
         :X_d_2dash, :X_q_2dash,
         :T_d_dash,  :T_q_dash, :Sn )

    ode_gens_generic_sequence_order =
        (:components_data, :gen)
        
    govs_and_avrs_sequence_order =
        ( :components_data,)
    
    govs_and_avrs_selections =
        ( :gov, :avr )
    
    #----------------------------------------
    
    ode_gens_generic_para =
         get_ode_gens_generic_para(
             plant_generators_data_from_json;
             sequence_order =
                 ode_gens_generic_sequence_order,
             selections =
                 ode_gens_generic_selections)
    
    #------------------------------------------------
    #------------------------------------------------

    (;generic_gens_para,
     generic_govs_para,
     generic_avrs_para) =
         get_generic_gens_avr_gov_para(
             plant_generators_data_from_json;
             gens_sequence_order =
                 ode_gens_generic_sequence_order,
            gens_selections =
                ode_gens_generic_selections,
            govs_and_avrs_sequence_order =
                govs_and_avrs_sequence_order,
            govs_and_avrs_selections =
                govs_and_avrs_selections)
    
    #------------------------------------------------
    #------------------------------------------------
        
    (branches_fbus, branches_tbus,
     branches_r, branches_x, branches_b,
     branches_ratio, branches_angle, branches_type) =
          get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    (branches_fbus,
     branches_tbus) =
         get_edges_fbus_tbus_by_json(
             edge_data_from_json)
    
    (edges_r, edges_x, edges_b,
     edges_ratio, edges_angle, edges_type) =
         get_edges_generic_data_by_json(
             edge_data_from_json )
    
    (Gs, Bs) =
        get_nodes_shunts_Gs_and_Bs_by_json(
            shunt_data_from_json)

    edges_size = length(edges_r)
    
    edges_r_x_b_ratio_angle_idx =
        get_edges_r_x_b_ratio_angle_idx(
            edges_size)
    
    #------------------------------------------------
    #------------------------------------------------

    (;edges_orientation,
     edges_Ybr_cal,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json,
                 edge_data_from_json,
                 shunt_data_from_json;
                 baseMVA =
                     baseMVA,
                 basekV =
                     basekV,
                 use_pu_in_PQ =
                     use_pu_in_PQ,
                 line_data_in_pu =
                     line_data_in_pu ),
             (:edges_orientation,
              :edges_Ybr_cal,
              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #------------------------------------------------
    #------------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend ) =
         NamedTupleTools.select(
             get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                 getproperty(
                     Ynet_wt_nodes_idx_wt_adjacent_nodes,
                     :Ynet)),
             (:Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend) )
    
    #------------------------------------------------
    #------------------------------------------------
    
    ode_gens_para =
        NamedTupleTools.select(
            ode_gens_generic_para,
            (:H,
             :D,
             :X_d,
             :X_q,
             :X_d_dash,
             :X_q_dash,
             :T_d_dash,
             :T_q_dash,
             :Sn))
    
    #------------------------------------------------
    #------------------------------------------------
        
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )
    
    dyn_pf_fun_kwd_net_idxs =
        NamedTupleTools.select(
            net_nodes_type_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))

    dyn_pf_fun_kwd_n2s_idxs =
        NamedTupleTools.select(
            get_dict_net_streamlined_idx_by_nodes_type_idxs(
                net_nodes_type_idxs ),
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))

    nodes_with_demands_idx =
        getproperty(net_nodes_type_idxs,
                    :nodes_with_demands_idx )
    
    #------------------------------------------------
    #------------------------------------------------
    
    dyn_pf_mismatch_vars_kwd_para =
        (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para )

    sta_pf_PQ_para =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA =
                baseMVA,
            use_pu_in_PQ =
                use_pu_in_PQ)

    
    gens_vh_slack_θh_para =
        get_gens_vh_slack_θh_para_by_json(
            plant_generators_data_from_json )

    
    sta_pf_vars_and_paras_idx =
        get_sta_pf_vars_and_paras_idx_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    
    pf_sta_ΔPQ_mismatch_parameters =
        get_pf_sta_ΔPQ_mismatch_parameters_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json;
            baseMVA = baseMVA,
            basekV = basekV,
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu)

    #------------------------------------------------
    #------------------------------------------------
    
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
     generic_model_vars_wt_i_dq_Idx_in_state,

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
     Pg_Png_Qng_Idx,

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
     plants_algeb_fun_idx_kwd_para,
     
     id_iq_pg_vh_Idx,
     ωref0_vref0_porder0_id_iq_vh_Idx) =
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
              :Pg_Png_Qng_Idx,

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
              :plants_algeb_fun_idx_kwd_para,

              :id_iq_pg_vh_Idx,
              :ωref0_vref0_porder0_id_iq_vh_Idx ))
            
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
    #----------------------------------------
    # Dynamic simulation
    #----------------------------------------
    #----------------------------------------

    #----------------------------------------
    # Init
    #----------------------------------------
    
    plants_init_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para ,
         comps_init_funs )
    
    plants_generic_states_init_wt_ref =
        plants_generic_model_init_func(
            gens_vh,
            gens_θh,
            pf_P_gens,
            pf_Q_gens,
            ωs;
            kwd_para =
                plants_init_kwd_para )

    (plants_states_init,
     plants_refs ) =
         NamedTupleTools.select(
             plants_generic_states_init_wt_ref,
             (:plants_states_init,
              :plants_refs))


    ( nt_vec_per_paras,
      vec_vec_per_paras ) =
          get_nt_vec_wt_vec_vec_per_paras(
        plants_refs ;
        nt_syms =
            (:ωs,
             :ω_ref,
             :v_ref,
             :p_order,
             :i_d,
             :i_q,
             :gen_vh ) )

    (ω_ref,
     v_ref,
     p_order,
     gens_i_d,
     gens_i_q ) =
         NamedTupleTools.select(
             nt_vec_per_paras, (
                 :ω_ref,
                 :v_ref,
                 :p_order,
                 :i_d,
                 :i_q ))

    #----------------------------------------
    # System model init
    #----------------------------------------

    u0_model_states_init =
        Float64[plants_states_init;
                vh;
                θh;
                gens_i_d;
                gens_i_q]
    
    #----------------------------------------
    # callbacks
    #----------------------------------------

    plants_cb_paras_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para,
         comps_callback_paras_funs )

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
    
    cb_states = cb_fun_make_state_callbacks(
        list_selected_plants_state_event_cb_paras,
        list_selected_plants_state_affect_cb_paras )

    #----------------------------------------
    # System model para
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
    # model para with callback
    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh]
    
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll = 
        Float64[ω_ref;
                v_ref;
                p_order;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    #----------------------------------------    
    # algebraic equation kwd_para
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
                     :X_q_dash, :Sn ),
                vec_datatype = Float64 ),
            (:ra,
              :X_d,
              :X_q,     
              :X_d_dash,
             :X_q_dash, :Sn ) )

    flat_vh_flat_θh_id_iq_u0 =
        [vh;
         θh;
         gens_i_d;
         gens_i_q ]
    
    dyn_pf_solver =
        (;use_init_u0,
         use_nlsolve,
         pf_alg,
         abstol,
         reltol )

         
     # dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     # dyn_pf_vh_vhf_Idx,
     # dyn_pf_θh_θhf_Idx,
    
    #--------------------------------------
    
    (Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         Ynet_wt_nodes_idx_wt_adjacent_nodes
    
    edges_orientation =
        getproperty(Ybr_cal_and_edges_orientation,
                    :edges_orientation)

    #--------------------------------------
    #--------------------------------------

    algebraic_generic_model_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes )
    
    algebraic_generic_model_sol_kwd_para =
        (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         flat_vh_flat_θh_id_iq_u0,
         dyn_pf_solver,
         algebraic_generic_model_kwd_para
         )

    #----------------------------------------
    #----------------------------------------
    # ODE system dyamanics simulation
    #----------------------------------------    
    #----------------------------------------
    
    #----------------------------------------
    # generic model_dynamics_kwd_para
    #----------------------------------------

    # ωref0_vref0_porder0_id_iq_vh_Idx =
    #     dyn_ωref0_vref0_porder0_id_iq_vh_Idx

    dae_plants_kwd_para =
        (;state_vars_idx,

         ωref0_vref0_porder0_id_iq_vh_Idx,
         
         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         vec_comp_states_Idx,
         
         avrs_govs_cb_sw,
         avrs_govs_cb_sw_Idx,

         dae_comps_dyn_funs,
         comps_output_funs,
         ωs)

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

    plants_kwd_para =
        (;state_vars_idx,

         ωref0_vref0_porder0_id_iq_vh_Idx,

         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         vec_comp_states_Idx,

         avrs_govs_cb_sw,
         avrs_govs_cb_sw_Idx,


         comps_dyn_funs,
         comps_output_funs,

        ωs) 

    generic_system_dynamics_kwd_para =
        (;
         gens_nodes_idx,
         ωs,
         loc_load_exist,
         state_vars_idx,

         id_iq_pg_vh_Idx,

         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

         gens_state_vars_idx_in_state,
         state_vars_and_i_dq_Idx_in_state,                  
         state_algebraic_vars_Idx_in_state,

         dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         dae_plants_kwd_para,
         ode_plants_kwd_para,

         plants_kwd_para,

         algebraic_generic_model_kwd_para,

         algebraic_generic_model_sol_kwd_para,

         state_labels,
         algebraic_vars_labels,
         network_vars_labels)

    return (;u0_model_states_init,
            model_syms,
            model_mass_matrix,
            model_bool_dae_vars,
            
            algebraic_state_sym,
            generic_state_sym,            
            
            ode_gens_mass_matrix,
            ode_gens_bool_dae_vars,

            nodes_names,
            gens_nodes_names,
            non_gens_nodes_names,

            SM_gens_nodes_names,
            SC_gens_nodes_names,

            plants_states_init,

            plants_cb_paras_switches,
            avrs_govs_cb_sw,
            avrs_govs_cb_sw_Idx,
            cb_states,

            ωref0_vref0_porder0_id_iq_vh,
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

            algebraic_generic_model_sol_kwd_para,
            generic_system_dynamics_kwd_para,
            
            pf_P_gens,
            pf_Q_gens,
            vh,
            θh,
            gens_vh,
            gens_θh,
            
            ω_ref,
            v_ref,
            p_order,
            gens_i_d,
            gens_i_q,

            pf_generic_gens_para,
            flat_vh_flat_θh_id_iq_u0,
            
            generic_gens_para,
            generic_avrs_para,
            generic_govs_para,
            
            # Diag
            algebraic_generic_model_kwd_para,
            dae_plants_kwd_para,
            ode_plants_kwd_para,
            plants_kwd_para,

            plants_states_syms,
            gens_nodes_idx,
            
            state_labels,
            algebraic_vars_labels,
            network_vars_labels,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs,
            
            Ybr_cal_and_edges_orientation,
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            
            edges_r_x_b_ratio_angle_idx,
            
            Ynet_rows_Idxs_in_flattend,
            Ynet_real_imag_Idxs_in_flattend)
                
end


# """
#     get_status_steady_state_parameters(
#         net_data_by_components_file;
#         <keyword arguments>)

# These functions are used to simplify  generation of
# parameters or data of a system that are needed for
# simulation.


# The parameters are pre-fault
# state steady state parameters, fault
# state steady state parameters, post-fault
# state steady state parameters,

# `get_Yint_and_Yred_matrices`

# `get_net_status_Yint_and_Yred_matrices`

# `get_a_status_Yint_Yred`

# `get_status_steady_state_parameters`,

# `get_a_status_steady_state_data`

# """

"""
    get_status_steady_state_parameters(
        net_data_by_components_file;
        <keyword arguments>)

It is used to simplify generation of parameters or data of a system that are needed for simulation.


# Arguments

- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,

`use_init_u0`: the boolean variable that determines if initial state u0 should be used in a power flow.
`use_nlsolve`: the boolean variable that determines if `nlsolve` should be used in power flow.

`pf_alg`: power flow solver

`abstol`: the absolute error tolerance
`reltol`: the relative error tolerance

`on_fault_time`: the on fault time
`clear_fault_time`: the clear fault time

`list_fault_point_from_node_a::Vector{Float64}=[0.3]`: the list containing a ratio of the fault point from the source (from) orientation of lines
`list_fault_resistance::Vector{Float64} = [0.001]`: the list containing fault resistances of each fault in the network.
`list_no_line_circuit::Vector{Float64} = [4]`: the list containing the number of circuits per faulted lines.

`list_edges_to_have_fault::Vector{Int64} = [2]`: the list containing indices of lines that should have a fault.  
`clear_fault_selection_list::Vector{Int64} = [1]`: the list containing indices of faulted lines to be cleared in `list_edges_to_have_fault`.

`with_faults::Bool=false`: a legacy boolean variable.

"""
function get_status_steady_state_parameters(
    net_data_by_components_file;
    components_libs_dir = nothing,
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,    
    use_nlsolve = false,
    
    pf_alg = NewtonRaphson(),
    
    abstol = 1e-12,    
    reltol = 1e-12,
    
    on_fault_time = 9.0,
    clear_fault_time = 9.02,
       
    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [4],

    list_edges_to_have_fault = [ 2 ], 
    clear_fault_selection_list = [ 1 ],

    with_faults = false )
    
    #--------------------------------------    
    
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

    (;plant_generators_data_from_json,
     plant_loads_data_from_json,
     plant_transmission_data_from_json,
     edge_data_from_json,
     shunt_data_from_json,
     baseMVA_data_from_json) =
        NamedTupleTools.select(
            get_net_data_by_components_from_json_file(
                net_data_by_components_file;
                in_components_type_sym =
                    false ),
            (:plant_generators_data_from_json,
             :plant_loads_data_from_json,
             :plant_transmission_data_from_json,
             :edge_data_from_json,
             :shunt_data_from_json,
             :baseMVA_data_from_json))

    #----------------------------------------

    baseMVA = baseMVA_data_from_json
    
    #------------------------------------------------
    
    (;
     P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist ) =
         get_pf_PQ_param_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json;
             baseMVA =
                 baseMVA,
             use_pu_in_PQ =
                 use_pu_in_PQ )

    
    #----------------------------------------
    #----------------------------------------

    ode_gens_para_selections  =
        (:H,
         :D,
         :X_d,
         :X_q,                  
         :X_d_dash,
         :X_q_dash,
         :T_d_dash,
         :T_q_dash, :Sn )

    ode_gens_para_sequence_order =
        (:components_data,
         :gen)
        
    ode_gens_generic_selections =
        (:H,
         :D,
         :ra,
         :xℓ,
         :X_d,
         :X_q,
         
         :X_d_dash,
         :X_q_dash,
         
         :X_d_2dash,
         :X_q_2dash,
         
         :T_d_dash,
         :T_q_dash,
         
         :Sn,
         
         :T_d_2dash,
         :T_q_2dash,
         
         :vh,
         :P,
         :Q,
         
         :Pmin,
         :Pmax,
         :Qmin,
         :Qmax,
         :vmin,
         :vmax )

    ode_gens_generic_sequence_order =
        (:components_data, :gen)
        
    govs_and_avrs_sequence_order =
        ( :components_data,)
    
    govs_and_avrs_selections =
        ( :gov, :avr )
    
    #----------------------------------------
    
    ode_gens_generic_para =
         get_ode_gens_generic_para(
             plant_generators_data_from_json;
             sequence_order =
                 ode_gens_generic_sequence_order,
             selections =
                 ode_gens_generic_selections)
    
    #------------------------------------------------
    #------------------------------------------------

    (;generic_gens_para,
     generic_govs_para,
     generic_avrs_para) =
         get_generic_gens_avr_gov_para(
             plant_generators_data_from_json;
             gens_sequence_order =
                 ode_gens_generic_sequence_order,
            gens_selections =
                ode_gens_generic_selections,
            govs_and_avrs_sequence_order =
                govs_and_avrs_sequence_order,
            govs_and_avrs_selections =
                govs_and_avrs_selections)

    #------------------------------------------------
    #------------------------------------------------

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
                     :X_q_dash, :Sn ),
                vec_datatype = Float64 ),
            (:ra,
             :X_d,
             :X_q,     
             :X_d_dash,
             :X_q_dash, :Sn ) )

    #------------------------------------------------
    #------------------------------------------------

     opf_generic_gens_para =
        NamedTupleTools.select(
            get_selected_vec_nt_to_vec_vec(
                generic_gens_para,
                nothing;
                selections =
                    (:vh,
                     :P, :Q,
                     :Pmin, :Pmax,
                     :Qmin, :Qmax,
                     :vmin, :vmax,
                     :Sn ),
                vec_datatype = Float64 ),
            (:vh,
             :P, :Q,
             :Pmin, :Pmax,
             :Qmin, :Qmax,
             :vmin, :vmax,
             :Sn ) )

    #------------------------------------------------
    #------------------------------------------------

    
    (branches_fbus,
     branches_tbus,
     branches_r,
     branches_x,
     branches_b,
     branches_ratio,
     branches_angle,
     branches_type) =
          get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    (branches_fbus,
     branches_tbus) =
         get_edges_fbus_tbus_by_json(
             edge_data_from_json)
    
    # (edges_r, edges_x, edges_b,
    #  edges_ratio, edges_angle, edges_type) =
    #      get_edges_generic_data_by_json(
    #          edge_data_from_json )

    (edges_fbus, edges_tbus,
     edges_r, edges_x, edges_b,
     edges_ratio, edges_angle, edges_type) =
         get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )
    
    # (Gs, Bs) =
    #     get_nodes_shunts_Gs_and_Bs_by_json(
    #         shunt_data_from_json)
    
    (shunt_idx, Gs, Bs) =
        get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
            shunt_data_from_json)
    
    #------------------------------------------------
    #------------------------------------------------

    (;edges_orientation,
     edges_Ybr_cal,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json,
                 edge_data_from_json,
                 shunt_data_from_json;
                 baseMVA =
                     baseMVA,
                 basekV =
                     basekV,
                 use_pu_in_PQ =
                     use_pu_in_PQ,
                 line_data_in_pu =
                     line_data_in_pu ),
             (:edges_orientation,
              :edges_Ybr_cal,
              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))
    
    #------------------------------------------------
    #------------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend ) =
         NamedTupleTools.select(
             get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                 getproperty(
                     Ynet_wt_nodes_idx_wt_adjacent_nodes,
                     :Ynet)),
             (:Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend) )

    #------------------------------------------------
    #------------------------------------------------
    
    ode_gens_para =
        NamedTupleTools.select(
            ode_gens_generic_para,
            (:H,
             :D,
             :X_d,
             :X_q,
             :X_d_dash,
             :X_q_dash,
             :T_d_dash,
             :T_q_dash, :Sn))
    
    #------------------------------------------------
    #------------------------------------------------
        
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )
    
    dyn_pf_fun_kwd_net_idxs =
        NamedTupleTools.select(
            net_nodes_type_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))
    
    dyn_pf_fun_kwd_n2s_idxs =
        NamedTupleTools.select(
            get_dict_net_streamlined_idx_by_nodes_type_idxs(
                net_nodes_type_idxs ),
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))

    nodes_with_demands_idx =
        getproperty(net_nodes_type_idxs,
                    :nodes_with_demands_idx )
    
    #------------------------------------------------
    #------------------------------------------------
    
    dyn_pf_mismatch_vars_kwd_para =
        (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para )

    sta_pf_PQ_para =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA =
                baseMVA,
            use_pu_in_PQ =
                use_pu_in_PQ)

    
    gens_vh_slack_θh_para =
        get_gens_vh_slack_θh_para_by_json(
            plant_generators_data_from_json )

    
    sta_pf_vars_and_paras_idx =
        get_sta_pf_vars_and_paras_idx_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    
    pf_sta_ΔPQ_mismatch_parameters =
        get_pf_sta_ΔPQ_mismatch_parameters_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json;
            baseMVA = baseMVA,
            basekV = basekV,
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu)

    #------------------------------------------------
    #------------------------------------------------

    list_faulted_line_a_b_orientation =
        edges_orientation[
            list_edges_to_have_fault  ] 
    
    #--------------------------------------

    no_lines_fault =
        length(list_faulted_line_a_b_orientation)

    #--------------------------------------

    no_current_lines_fault =
        no_lines_fault - length(
            clear_fault_selection_list)
    
    #----------------------------------------

    no_cleared_lines_fault =
        length(clear_fault_selection_list)


    #----------------------------------------
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
     generic_model_vars_wt_i_dq_Idx_in_state,

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

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,     

     model_mass_matrix,     
     model_bool_dae_vars,

     ode_gens_mass_matrix,
     ode_gens_bool_dae_vars,
    
     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     dyn_vh_id_iq_V_ref_Tm_Idx,
     dyn_V_ref_Tm_id_iq_vh_Idx,

     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
     dyn_vh_id_iq_ωref0_vref0_porder0_Idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,

     dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

     pf_vh_θh_idx_and_idx2Idx,
     
     dyn_pf_flat_vh_flat_θh_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,
                 
     system_states_idx_kwd_para,
     system_paras_idx_kwd_para,
     plants_dyn_fun_idx_kwd_para,
     plants_algeb_fun_idx_kwd_para,

     id_iq_pg_vh_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     scale_Pg_Png_Qng_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx) =
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

              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :model_mass_matrix,     
              :model_bool_dae_vars,

              :ode_gens_mass_matrix,
              :ode_gens_bool_dae_vars,

              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,

              :dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :dyn_vh_id_iq_V_ref_Tm_Idx,
              :dyn_V_ref_Tm_id_iq_vh_Idx,

              :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
              :dyn_vh_id_iq_ωref0_vref0_porder0_Idx,

              :ωref0_vref0_porder0_id_iq_vh_Idx,

              :dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              :pf_vh_θh_idx_and_idx2Idx,
              
              :dyn_pf_flat_vh_flat_θh_Idx,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :system_states_idx_kwd_para,
              :system_paras_idx_kwd_para,
              :plants_dyn_fun_idx_kwd_para,
              :plants_algeb_fun_idx_kwd_para,

              :id_iq_pg_vh_Idx,

              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :scale_Pg_Png_Qng_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx))
    
        
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
    # Powerflow func and prob
    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (;Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_dyn_sol_kwd_para =
        (;loc_load_exist,
         sta_pf_PQ_para,
         ode_gens_generic_para,
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs)
    
    #----------------------------------------
    # callbacks
    #----------------------------------------

    plants_cb_paras_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para,
         comps_callback_paras_funs )

    #----------------------------------------
    # Init
    #----------------------------------------
    
    plants_init_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para ,
         comps_init_funs )

    
    #----------------------------------------
    # Yred and Yint
    #----------------------------------------

    X_d_dash =
        getproperty(
            ode_gens_para,
            :X_d_dash) 

    y_aug_kw_para =
        (;loc_load_exist,
         X_d_dash,    
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs)
    
    #----------------------------------------
    # dyn solver
    #----------------------------------------
    
    dyn_pf_solver =
        (;use_init_u0,
         use_nlsolve,
         pf_alg,
         abstol,
         reltol )

    #----------------------------------------
    # faults
    #----------------------------------------
    
    (Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         Ynet_wt_nodes_idx_wt_adjacent_nodes
    
    edges_orientation =
        getproperty(
            Ybr_cal_and_edges_orientation,
            :edges_orientation)

    #--------------------------------------
    
    on_fault_net_para =
         make_lines_faults_data_set(
             Ynet,
             nodes_idx_with_adjacent_nodes_idx,

             all_nodes_idx,
             n2s_all_nodes_idx,

             list_faulted_line_a_b_orientation ,
             list_fault_point_from_node_a,
             list_fault_resistance,
             list_no_line_circuit )

    (;fault_nodes_idx,) =
        NamedTupleTools.select(
            on_fault_net_para,
            (:fault_nodes_idx,))
    
    #--------------------------------------

    cleared_selected_lines_faults_net_para =
        get_cleared_selected_lines_faults_data_set(
            clear_fault_selection_list;
            deepcopy(on_fault_net_para)...)

    #--------------------------------------
    
    dyn_pf_vh_θh_id_iq_vhf_θhf_Idx =
        get_generic_vh_θh_id_iq_vhf_θhf_Idx(
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault =
                no_lines_fault)


    dyn_pf_vh_vhf_θh_θhf_id_iq_Idx =
         get_generic_vh_vhf_θh_θhf_id_iq_Idx(
            gens_nodes_idx,
            all_nodes_idx;
             no_lines_fault =
                 no_lines_fault)


    state_algebraic_vars_wt_fault_Idx_in_state =
        get_state_algebraic_vars_wt_fault_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault =
                no_lines_fault)
    
    #--------------------------------------
    
    number_of_faults =
        length(list_edges_to_have_fault)

    twice_number_of_faults =
        2 * number_of_faults

    fault_vh_sym =
        generate_labels_by_nodes_idxs_and_vars(
            fault_nodes_idx,
            [:vh];
            label_prefix = "bus")

    fault_θh_sym =
        generate_labels_by_nodes_idxs_and_vars(
            fault_nodes_idx,
            [:θh];
            label_prefix = "bus")

    fault_nodes_sym =
        [fault_vh_sym;
         fault_θh_sym ]  

    fault_algeb_init =[
        ones(number_of_faults);
        zeros(number_of_faults) ]
    
    #--------------------------------------
        
    edges_size = length(edges_r)

    edges_r_x_b_ratio_angle_idx =
        get_edges_r_x_b_ratio_angle_idx(
            edges_size)
    
    #--------------------------------------
    
    return (;baseMVA,
            basekV,
            
            gens_nodes_idx,
            loc_load_exist,
            state_vars_idx,
            
            state_labels,
            algebraic_vars_labels,
            network_vars_labels,

            plants_states_syms,
            generic_state_sym, 
            algebraic_state_sym,            
            model_syms,

            nodes_names,
            gens_nodes_names,
            non_gens_nodes_names,
            SM_gens_nodes_names,
            SC_gens_nodes_names,

            model_mass_matrix,     
            model_bool_dae_vars,

            ode_gens_mass_matrix,
            ode_gens_bool_dae_vars,

            gens_state_vars_idx_in_state,
            state_vars_and_i_dq_Idx_in_state,

            state_vars_and_i_dq_wt_fault_Idx_in_state,

            state_algebraic_vars_Idx_in_state,
            state_algebraic_vars_wt_fault_Idx_in_state,

            pf_vh_θh_idx_and_idx2Idx,
            
            dyn_pf_flat_vh_flat_θh_Idx,
            
            dyn_pf_flat_vh_flat_θh_id_iq_Idx,
            
            dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
            dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

            dyn_pf_vh_vhf_Idx,
            dyn_pf_θh_θhf_Idx,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs,
            
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

            ωref0_vref0_porder0_id_iq_vh_Idx,
            dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
            
            vec_comp_states_Idx,

            comps_callback_paras_funs,     
            comps_output_funs,
            comps_dyn_funs,

            ode_comps_dyn_funs,
            dae_comps_dyn_funs,
            
            comps_init_funs,            
            generic_gens_para,
            generic_avrs_para,
            generic_govs_para,

            ode_gens_para,
            
            sta_pf_PQ_para,
            pf_PQ_param,
            pf_generic_gens_para,

            edges_fbus,
            edges_tbus,
            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle,
            edges_type,
            shunt_idx,
            Gs,
            Bs,
            
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            Ybr_cal_and_edges_orientation,
            
                        
            dyn_pf_solver,

            kwd_sta_sta_ΔPQ_sol_by_json,
            generic_red_sol_kwd_para,
            generic_dyn_sol_kwd_para,
            
            plants_init_kwd_para,
            plants_cb_paras_kwd_para,

            y_aug_kw_para,

            with_faults,

            clear_fault_selection_list,
            no_lines_fault,
            no_cleared_lines_fault,
            no_current_lines_fault,            
            
            list_fault_resistance,

            on_fault_net_para,
            cleared_selected_lines_faults_net_para,

            fault_nodes_sym,
            fault_algeb_init,

            Pg_Qg_Png_Qng_Pll_Qll_Idx,
            Png_Qng_Pll_Qll_Idx,
            Pg_Png_Qng_Idx,

            edges_r_x_b_ratio_angle_idx,

            Ynet_rows_Idxs_in_flattend,
            Ynet_real_imag_Idxs_in_flattend,

            id_iq_pg_vh_Idx,
            
            scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
            scale_Pg_Png_Qng_Idx,
            dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx)

end

# ------------------------------------------------------

"""
    get_a_status_steady_state_data(
        system_status;
        <keyword arguments>)

It is used to simplify generation of parameters or data of a system that are needed for simulation.

# Arguments

- `net_data_by_components_file`: the network data file
- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,

`use_init_u0`: the boolean variable that determines if initial state u0 should be used in a power flow.
`use_nlsolve`: the boolean variable that determines if `nlsolve` should be used in power flow.
`pf_alg`: power flow solver
`abstol`: the absolute error tolerance
`reltol`: the relative error tolerance
`on_fault_time`: the on fault time
`clear_fault_time`: the clear fault time
`list_fault_point_from_node_a::Vector{Float64}=[0.3]`: the list containing a ratio of the fault point from the source (from) orientation of lines
`list_fault_resistance::Vector{Float64} = [0.001]`: the list containing fault resistances of each fault in the network.
`list_no_line_circuit::Vector{Float64} = [4]`: the list containing the number of circuits per faulted lines
`list_edges_to_have_fault::Vector{Int64} = [2]`: the list containing indices of lines that should have a fault.  
`clear_fault_selection_list::Vector{Int64} = [1]`: the list containing indices of faulted lines to be cleared in `list_edges_to_have_fault`.
`ode_alg`: the ode solver.
dae_alg: the dae solver.
dt: the solve interval.
`with_faults::Bool=false`: a legacy boolean variable.


- `system_status`: can take the form of the following
   - `:pre_fault_state`
   - `:fault_state` 
   - `:post_fault_state`

"""
function get_a_status_steady_state_data(
    system_status;
    with_faults =
        false,
    
    # pre_fault_state
    # fault_state 
    # post_fault_state ,
    # system_status = :pre_fault_state,

    net_data_by_components_file =
        nothing,
    components_libs_dir =
        nothing,
    
    timespan         = 10.0,
    on_fault_time    = 9.0,
    clear_fault_time = 9.001,
    
    list_fault_point_from_node_a = [0.3],
    list_fault_resistance        = [0.001],
    list_no_line_circuit         =  [4],

    list_edges_to_have_fault   = [ 2 ],
    clear_fault_selection_list = [1],
    
    basekV          = 1.0,    
    use_pu_in_PQ    = true,
    line_data_in_pu = true,
    
    #--------------------------------------    

    use_init_u0 = false,
    
    use_nlsolve = false,
    
    pf_alg        = NewtonRaphson(),

    #--------------------------------------    
    
    ode_alg       = Rodas4(),
    # ode_alg     = FBDF()
    # ode_alg     = QNDF()
    # ode_alg     = radau()
    # ode_alg     = RadauIIA5()
    # ode_alg     = DFBDF()

    dae_alg       = IDA(),
    
    abstol        = 1e-12,
    reltol        = 1e-12,
    
    dt            = 0.01 )

    #--------------------------------------
    
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

    
    tspan         = (0.0, timespan)
    
    sim_timespan  = (0.0, timespan)
    
    plot_timespan = (0.0, timespan)
    
    #----------------------------------------        
    
    fault_status =
        (no_fault = 0,
         on_fault = 1,
         clear_fault = 2,
         partial_clear_fault = 3)

    system_fault_status =
        [fault_status.no_fault]
    
    #----------------------------------------

    (;baseMVA,
     basekV,

     gens_nodes_idx,
     loc_load_exist,
     state_vars_idx,

     state_labels,
     algebraic_vars_labels,
     network_vars_labels,

     plants_states_syms,
     generic_state_sym, 
     algebraic_state_sym,            
     model_syms,

     nodes_names,
     gens_nodes_names,
     non_gens_nodes_names,
     SM_gens_nodes_names,
     SC_gens_nodes_names,

     model_mass_matrix,     
     model_bool_dae_vars,

     ode_gens_mass_matrix,
     ode_gens_bool_dae_vars,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,

     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     pf_vh_θh_idx_and_idx2Idx,
     
     dyn_pf_flat_vh_flat_θh_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

     vec_comp_states_Idx,

     comps_callback_paras_funs,     
     comps_output_funs,
     comps_dyn_funs,

     ode_comps_dyn_funs,
     dae_comps_dyn_funs,

     comps_init_funs,            
     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,

     ode_gens_para,

     sta_pf_PQ_para,
     pf_PQ_param,
     pf_generic_gens_para,

     edges_fbus,
     edges_tbus,
     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type,
     shunt_idx,
     Gs,
     Bs,
     
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     Ybr_cal_and_edges_orientation,

     dyn_pf_solver,

     kwd_sta_sta_ΔPQ_sol_by_json,
     generic_red_sol_kwd_para,
     
     generic_dyn_sol_kwd_para,

     plants_init_kwd_para,
     plants_cb_paras_kwd_para,

     y_aug_kw_para,

     with_faults,

     clear_fault_selection_list,
     no_lines_fault,
     no_cleared_lines_fault,
     no_current_lines_fault,
     
     list_fault_resistance,

     on_fault_net_para,
     cleared_selected_lines_faults_net_para,

     fault_nodes_sym,
     fault_algeb_init,

     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     edges_r_x_b_ratio_angle_idx,

     Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend,
     id_iq_pg_vh_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     scale_Pg_Png_Qng_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx) =
         NamedTupleTools.select(
             get_status_steady_state_parameters(
                 net_data_by_components_file;
                 components_libs_dir =
                     components_libs_dir,
                 basekV =
                     basekV,    
                 use_pu_in_PQ =
                     use_pu_in_PQ,
                 line_data_in_pu =
                     line_data_in_pu,

                 use_init_u0 =
                     use_init_u0,    
                 use_nlsolve =
                     use_nlsolve,

                 pf_alg =
                     pf_alg,

                 abstol =
                     abstol,    
                 reltol =
                     reltol,

                 on_fault_time =
                     on_fault_time,
                 clear_fault_time =
                     clear_fault_time,

                 list_fault_point_from_node_a =
                     list_fault_point_from_node_a,
                 list_fault_resistance =
                     list_fault_resistance,
                 list_no_line_circuit =
                     list_no_line_circuit,

                 list_edges_to_have_fault =
                     list_edges_to_have_fault, 
                 clear_fault_selection_list =
                     clear_fault_selection_list,

                 with_faults =
                     with_faults),
             (:baseMVA,
              :basekV,

              :gens_nodes_idx,
              :loc_load_exist,
              :state_vars_idx,

              :state_labels,
              :algebraic_vars_labels,
              :network_vars_labels,

              :plants_states_syms,
              :generic_state_sym, 
              :algebraic_state_sym,            
              :model_syms,

              :nodes_names,
              :gens_nodes_names,
              :non_gens_nodes_names,
              :SM_gens_nodes_names,
              :SC_gens_nodes_names,

              :model_mass_matrix,     
              :model_bool_dae_vars,

              :ode_gens_mass_matrix,
              :ode_gens_bool_dae_vars,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,

              :state_vars_and_i_dq_wt_fault_Idx_in_state,

              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,
              :pf_vh_θh_idx_and_idx2Idx,
              
              :dyn_pf_flat_vh_flat_θh_Idx,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

            :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

              :ωref0_vref0_porder0_id_iq_vh_Idx,
              :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

              :vec_comp_states_Idx,

              :comps_callback_paras_funs,     
              :comps_output_funs,
              :comps_dyn_funs,

              :ode_comps_dyn_funs,
              :dae_comps_dyn_funs,

              :comps_init_funs,            
              :generic_gens_para,
              :generic_avrs_para,
              :generic_govs_para,

              :ode_gens_para,

              :sta_pf_PQ_para,
              :pf_PQ_param,
              :pf_generic_gens_para,

              :edges_fbus,
              :edges_tbus,
              :edges_r,
              :edges_x,
              :edges_b,
              :edges_ratio,
              :edges_angle,
              :edges_type,
              :shunt_idx,
              :Gs,
              :Bs,
              
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :Ybr_cal_and_edges_orientation,

              :dyn_pf_solver,

              :kwd_sta_sta_ΔPQ_sol_by_json,
              :generic_red_sol_kwd_para,
              
              :generic_dyn_sol_kwd_para,

              :plants_init_kwd_para,
              :plants_cb_paras_kwd_para,

              :y_aug_kw_para,

              :with_faults,

              :clear_fault_selection_list,
              :no_lines_fault,
              :no_cleared_lines_fault,
              :no_current_lines_fault,
              
              :list_fault_resistance,

              :on_fault_net_para,
              :cleared_selected_lines_faults_net_para,

              :fault_nodes_sym,
              :fault_algeb_init,
              
              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,

              :edges_r_x_b_ratio_angle_idx,

              :Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend,
              :id_iq_pg_vh_Idx,

              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :scale_Pg_Png_Qng_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,) )

    #----------------------------------------
    #----------------------------------------
    
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol)  =
         NamedTupleTools.select(
             dyn_pf_solver,
             (:use_init_u0,
              :use_nlsolve,
              :pf_alg,
              :abstol,
              :reltol))
    
    selected_gens_paras_by_vec =
        get_selected_vec_nt_to_vec_vec(
            generic_gens_para,
            nothing;
            selections =
                (:D, :H, :xℓ, :X_d_dash),
            vec_datatype =
                Float64 )

    (gens_H,
     gens_D,
     X_d_dash) =
        NamedTupleTools.select(
            selected_gens_paras_by_vec,
            (:H,
             :D,
             :X_d_dash))

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
    # * static power flow 
    #----------------------------------------    
    
    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------

    # generic_red_sol_kwd_para =
    #     (;Ybr_cal_and_edges_orientation,
    #       sta_pf_PQ_para,
    #       ode_gens_generic_para,
    #       pf_kw_para ) 
    
    generic_results_pf_sta_red_sol = sta_pf_red_sol =
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
     gens_θh,
     
     gens_ang_E,
     gens_mag_E) =
        NamedTupleTools.select(
            sta_pf_red_sol,
            (:pf_P_gens,
             :pf_Q_gens,
             :vh,
             :θh,
             
             :gens_vh,
             :gens_θh,
             
             :gens_ang_E,
             :gens_mag_E) )
        
    post_sta_PQ =
        (;pf_P_gens,
         pf_Q_gens,
         P_non_gens,
         Q_non_gens,
         P_g_loc_load,
         Q_g_loc_load )


    Pg_Qg_Png_Qng_Pll_Qll =
        [pf_P_gens;
         pf_Q_gens;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------
    # Yred and Yint
    #----------------------------------------
    
    (Yred, Yint ) =
        get_Yint_and_Yred_matrices(
            vh,
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load;
            y_aug_kw_para =
                y_aug_kw_para )    
    
    #----------------------------------------
    # prefault Init
    #----------------------------------------
    
    # plants_init_kwd_para =
    #     (generic_gens_para ,
    #      generic_avrs_para,
    #      generic_govs_para ,
    #      comps_init_funs )
    
    plants_generic_states_init_wt_ref =
        plants_generic_model_init_func(
            gens_vh,
            gens_θh,
            pf_P_gens,
            pf_Q_gens,
            ωs;
            kwd_para =
                plants_init_kwd_para )

    (plants_states_init,
     plants_refs ) =
         NamedTupleTools.select(
             plants_generic_states_init_wt_ref,
             (:plants_states_init,
              :plants_refs))


    ( nt_vec_per_paras,
      vec_vec_per_paras ) =
          get_nt_vec_wt_vec_vec_per_paras(
        plants_refs ;
        nt_syms =
            (:ωs,
             :ω_ref,
             :v_ref,
             :p_order,
             :i_d,
             :i_q,
             :gen_vh ) )

    (ω_ref,
     v_ref,
     p_order,
     gens_i_d,
     gens_i_q ) =
         NamedTupleTools.select(
             nt_vec_per_paras, (
                 :ω_ref,
                 :v_ref,
                 :p_order,
                 :i_d,
                 :i_q))

    gens_id = gens_i_d
    
    gens_iq = gens_i_q
    
    #----------------------------------------
    # prefault  System model init
    #----------------------------------------

    flat_vh_flat_θh_id_iq_u0 =
        [vh;
         θh;
         gens_i_d;
         gens_i_q]
    
    flat_vh_flat_θh_id_iq_vfh_θfh =
        [flat_vh_flat_θh_id_iq_u0;
         fault_algeb_init]
    
    u0_model_states_init =
        Float64[plants_states_init;
                flat_vh_flat_θh_id_iq_u0]

    u0_model_states_init_wt_fault_algeb_init =
        [plants_states_init;
         flat_vh_flat_θh_id_iq_u0;
         fault_algeb_init]

    #----------------------------------------
    #----------------------------------------

    fault_model_bool_dae_vars =
        DAE_BoolVector(
            0,
            length(fault_nodes_sym) )
    
    model_bool_dae_vars_wt_fault =
        [model_bool_dae_vars;
         fault_model_bool_dae_vars]

    model_syms_wt_fault =
        [model_syms;
         fault_nodes_sym]

    u0_model_states_init_wt_fault =
        [u0_model_states_init;
         fault_algeb_init ]
    

    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))

    #----------------------------------------

    gens_δ = u0_model_states_init[
        δ_idx_in_state]

    gens_ed_dash = u0_model_states_init[
        ed_dash_idx_in_state]

    gens_eq_dash = u0_model_states_init[
        eq_dash_idx_in_state]

    #----------------------------------------    
    
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    #----------------------------------------
    # model para 
    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh]
    
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll = 
        Float64[ω_ref;
                v_ref;
                p_order;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]

    slack_gens_vh_θh_gens_vh_non_slack_gens_vh =
        get_slack_gens_vh_θh_gens_vh_non_slack_gens_vh(
            vh,
            θh,
            dyn_pf_fun_kwd_net_idxs,
            dyn_pf_fun_kwd_n2s_idxs)
    
    #----------------------------------------    
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
                     plants_cb_paras_kwd_para),
             (:plants_cb_paras_switches,
              :list_selected_plants_state_event_cb_paras,
              :list_selected_plants_state_affect_cb_paras,

              :plants_avr_gov_cb_para_sw_in_plant,
              :plants_avr_gov_cb_para_sw_idx_in_plant))

    cb_states = cb_fun_make_state_callbacks(
        list_selected_plants_state_event_cb_paras,
        list_selected_plants_state_affect_cb_paras )

    #----------------------------------------    

    
    (;
     nodes_idx_with_adjacent_nodes_idx,) =
         NamedTupleTools.select(
             Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:nodes_idx_with_adjacent_nodes_idx,))
    
    
    (;
     pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
     post_clear_fault_nodes_idx_with_adjacent_nodes_idx) =
         NamedTupleTools.select(
             cleared_selected_lines_faults_net_para,
             (
              :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
              :post_clear_fault_nodes_idx_with_adjacent_nodes_idx))
    
    #----------------------------------------
    # 
    #----------------------------------------

    algebraic_generic_model_kwd_para =
        (; loc_load_exist,

         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,

         Ynet_wt_nodes_idx_wt_adjacent_nodes
         )


    algebraic_generic_model_sol_kwd_para =
        (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         flat_vh_flat_θh_id_iq_u0,
         dyn_pf_solver,
         algebraic_generic_model_kwd_para)


    algebraic_generic_model_wt_fault_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
         dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
         
         on_fault_net_para,

         clear_fault_selection_list,

         no_lines_fault,
         no_cleared_lines_fault,

         list_fault_resistance,

         with_faults,

         cleared_selected_lines_faults_net_para )
    
    algebraic_generic_model_wt_fault_sol_kwd_para =
        (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         flat_vh_flat_θh_id_iq_u0,
         dyn_pf_solver,
         algebraic_generic_model_wt_fault_kwd_para,
         cleared_selected_lines_faults_net_para
         )

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

    dae_plants_kwd_para =
        (;state_vars_idx,

         ωref0_vref0_porder0_id_iq_vh_Idx,
         
         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         vec_comp_states_Idx,
         
         avrs_govs_cb_sw,
         avrs_govs_cb_sw_Idx,
         
         dae_comps_dyn_funs,
         comps_output_funs,
         ωs)

    generic_system_dynamics_wt_fault_kwd_para =
        (;
         gens_nodes_idx,
         ωs,
         loc_load_exist,
         state_vars_idx,

         # id_iq_pg_vh_Idx,

         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

         gens_state_vars_idx_in_state,
         state_vars_and_i_dq_Idx_in_state,
         
         state_vars_and_i_dq_wt_fault_Idx_in_state,
         
         state_algebraic_vars_Idx_in_state,
         state_algebraic_vars_wt_fault_Idx_in_state,

         dyn_pf_flat_vh_flat_θh_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
         dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

         dyn_pf_vh_vhf_Idx,
         dyn_pf_θh_θhf_Idx,
         
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         ode_plants_kwd_para,
         dae_plants_kwd_para,

         algebraic_generic_model_wt_fault_kwd_para,
         algebraic_generic_model_wt_fault_sol_kwd_para,

         algebraic_generic_model_kwd_para,
         algebraic_generic_model_sol_kwd_para,

         no_lines_fault,
         no_current_lines_fault,
         
         cleared_selected_lines_faults_net_para,

         with_faults,
         generic_results_pf_sta_red_sol,

         Pg_Qg_Png_Qng_Pll_Qll_Idx,
         Png_Qng_Pll_Qll_Idx,
         Pg_Png_Qng_Idx,
         
         id_iq_pg_vh_Idx,

         nodes_idx_with_adjacent_nodes_idx,
         pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
         post_clear_fault_nodes_idx_with_adjacent_nodes_idx)
    
    #----------------------------------------
    # static prefault parameters
    #----------------------------------------

    static_prefault_paras =
        (;baseMVA,
         basekV,

         system_fault_status,

         loc_load_exist,
         sta_pf_PQ_para,
         
         sta_pf_red_sol,
         
         pf_P_gens,
         pf_Q_gens,
         
         vh, θh,
         gens_vh, gens_θh,
         gens_id, gens_iq,
         
         gens_ang_E, gens_mag_E,
         post_sta_PQ,
         Yred, Yint,
         
         plants_generic_states_init_wt_ref,
         plants_states_init, plants_refs,
         
         nt_vec_per_paras, vec_vec_per_paras,
         
         ω_ref, v_ref, p_order, gens_i_d, gens_i_q,
         
         flat_vh_flat_θh_id_iq_u0,
         flat_vh_flat_θh_id_iq_vfh_θfh,
         
         u0_model_states_init,

         gens_δ,
         gens_ed_dash,
         gens_eq_dash,

         δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,
         ωref0_vref0_porder0_id_iq_vh,
         ω_ref_v_ref_p_order_Png_Qng_Pll_Qll, 
         
         algebraic_generic_model_sol_kwd_para,
         algebraic_generic_model_wt_fault_sol_kwd_para,

         ode_plants_kwd_para,
         dae_plants_kwd_para,
         generic_system_dynamics_wt_fault_kwd_para,

         edges_fbus,
         edges_tbus,
         edges_r,
         edges_x,
         edges_b,
         edges_ratio,
         edges_angle,
         edges_type,
         shunt_idx,
         Gs,
         Bs,
         Ybr_cal_and_edges_orientation,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes,
         on_fault_net_para,
         cleared_selected_lines_faults_net_para,

         fault_nodes_sym,
         fault_algeb_init,

         avrs_govs_cb_sw,
         avrs_govs_cb_sw_Idx,
         
         cb_states,
         
         plants_cb_paras_switches,

         gens_H,
         gens_D,
         X_d_dash,
         
         state_labels,
         algebraic_vars_labels,
         network_vars_labels,

         plants_states_syms,
         generic_state_sym, 
         algebraic_state_sym,
         model_syms,

         nodes_names,
         gens_nodes_names,
         non_gens_nodes_names,
         SM_gens_nodes_names,
         SC_gens_nodes_names,

         model_bool_dae_vars_wt_fault,
         model_syms_wt_fault,         
         u0_model_states_init_wt_fault,

         model_mass_matrix,     
         model_bool_dae_vars,

         ode_gens_mass_matrix,
         ode_gens_bool_dae_vars,

         ωref0_vref0_porder0_id_iq_vh_Idx,
         dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         state_vars_idx,
         vec_comp_states_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         Pg_Qg_Png_Qng_Pll_Qll_Idx,
         Png_Qng_Pll_Qll_Idx,
         Pg_Png_Qng_Idx,

         edges_r_x_b_ratio_angle_idx,

         Ynet_rows_Idxs_in_flattend,
         Ynet_real_imag_Idxs_in_flattend,

         pf_vh_θh_idx_and_idx2Idx,
         
         dyn_pf_flat_vh_flat_θh_Idx,

         Pg_Qg_Png_Qng_Pll_Qll,

         slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

         scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
         scale_Pg_Png_Qng_Idx,
         dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

         generic_gens_para,
         ode_gens_para )
    
    #----------------------------------------        
    #----------------------------------------
    
    if system_status == :pre_fault_state

        pf_fun_mismatch =
           pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!

    elseif system_status ==  :fault_state 

        pf_fun_mismatch =
            fault_state_algebraic_generic_pf_ΔPQ_mismatch!

    elseif system_status ==  :post_fault_state

        pf_fun_mismatch =
          post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!
    else
        throw(" system status not known ")
    end

    #----------------------------------------    
    #----------------------------------------    
        
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll;
                kwd_para =
                  algebraic_generic_model_wt_fault_kwd_para),
                      flat_vh_flat_θh_id_iq_vfh_θfh,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    kwd_para =
              algebraic_generic_model_wt_fault_kwd_para)),
            flat_vh_flat_θh_id_iq_vfh_θfh,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll ),
                                    pf_alg )

        
    end

    #----------------------------------------    

    if system_status ==  :fault_state

        dyn_pf_sol_u =
            get_generic_results_dyn_pf_sol_u(
                sol,
                Ynet_wt_nodes_idx_wt_adjacent_nodes;
                generic_dyn_sol_kwd_para =
                    generic_dyn_sol_kwd_para,
                baseMVA = baseMVA,
                basekV = basekV,
                system_status =
                    system_status,
                on_fault_net_para =
                    on_fault_net_para,
                dyn_pf_vh_θh_id_iq_vhf_θhf_Idx =
                    dyn_pf_vh_θh_id_iq_vhf_θhf_Idx )

        
        (;
         dyn_pf_vhf_Idxs,
         dyn_pf_θhf_Idxs) = 
             NamedTupleTools.select(
                 dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                 (
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θhf_Idxs))
        
        fault_algeb_init =
            [sol[dyn_pf_vhf_Idxs];
             sol[dyn_pf_θhf_Idxs] ]
        
    elseif system_status ==  :post_fault_state
            
        dyn_pf_sol_u =
            get_generic_results_dyn_pf_sol_u(
                sol,
                Ynet_wt_nodes_idx_wt_adjacent_nodes;
                generic_dyn_sol_kwd_para =
                    generic_dyn_sol_kwd_para,
                baseMVA =
                    baseMVA,
                basekV =
                    basekV )
    else

        
        return  (;static_prefault_paras, )         
                
    end
        
    (dyn_pf_P_gens, dyn_pf_Q_gens,
     dyn_vh, dyn_θh,
     dyn_gens_vh, dyn_gens_θh,
     dyn_gens_id, dyn_gens_iq,
     dyn_gens_mag_E, dyn_gens_ang_E) =
        NamedTupleTools.select(
            dyn_pf_sol_u,
            (:pf_P_gens, :pf_Q_gens,
             :vh, :θh,
             :gens_vh, :gens_θh,
             :gens_id, :gens_iq,
             :gens_mag_E, :gens_ang_E) )    

    post_dyn_PQ =
        (;dyn_pf_P_gens,
         dyn_pf_Q_gens,
         P_non_gens,
         Q_non_gens,
         P_g_loc_load,
         Q_g_loc_load)
    
    (dyn_Yred, dyn_Yint ) =
        get_Yint_and_Yred_matrices(
            dyn_vh,
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load;
            y_aug_kw_para =
                y_aug_kw_para,
            fault_status =
                system_status ,
            cleared_selected_lines_faults_net_para =
                cleared_selected_lines_faults_net_para)    
    
    #----------------------------------------
    # Init
    #----------------------------------------        
    
    dyn_plants_generic_states_init_wt_ref =
        plants_generic_model_init_func(
            dyn_gens_vh,
            dyn_gens_θh,
            dyn_pf_P_gens,
            dyn_pf_Q_gens,
            ωs;
            kwd_para =
                plants_init_kwd_para )

    (dyn_plants_states_init,
     dyn_plants_refs ) =
         NamedTupleTools.select(
             dyn_plants_generic_states_init_wt_ref,
             (:plants_states_init,
              :plants_refs))

    ( dyn_nt_vec_per_paras,
      dyn_vec_vec_per_paras ) =
          get_nt_vec_wt_vec_vec_per_paras(
        dyn_plants_refs ;
        nt_syms =
            (:ωs,
             :ω_ref,
             :v_ref,
             :p_order,
             :i_d,
             :i_q,
             :gen_vh ) )


    (dyn_ω_ref,
     dyn_v_ref,
     dyn_p_order,
     dyn_gens_i_d,
     dyn_gens_i_q ) =
         NamedTupleTools.select(
             dyn_nt_vec_per_paras,
             (:ω_ref,
              :v_ref,
              :p_order,
              :i_d,
              :i_q ) )
    
    #----------------------------------------
    # System model init
    #----------------------------------------

    flat_vh_flat_θh_id_iq_u0 =
        [dyn_vh;
         dyn_θh;
         dyn_gens_i_d;
         dyn_gens_i_q]

    
    flat_vh_flat_θh_id_iq_vfh_θfh =
        [flat_vh_flat_θh_id_iq_u0;
         fault_algeb_init]

    
    u0_model_states_init =
        Float64[dyn_plants_states_init;
                flat_vh_flat_θh_id_iq_u0]

    
    u0_model_states_init_wt_fault_algeb_init =
        [dyn_plants_states_init;
         flat_vh_flat_θh_id_iq_u0;
         fault_algeb_init]
    
    #----------------------------------------
    # System model para
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))

    #----------------------------------------

    dyn_gens_δ = u0_model_states_init[
        δ_idx_in_state]

    dyn_gens_ed_dash = u0_model_states_init[
        ed_dash_idx_in_state]

    dyn_gens_eq_dash = u0_model_states_init[
        eq_dash_idx_in_state]

    #----------------------------------------    
    
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll =
        [dyn_gens_δ;
         dyn_gens_ed_dash;
         dyn_gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    #----------------------------------------
    # model para 
    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [dyn_ω_ref;
         dyn_v_ref;
         dyn_p_order;
         dyn_gens_i_d;
         dyn_gens_i_q;
         dyn_gens_vh]
    
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll = 
        Float64[dyn_ω_ref;
                dyn_v_ref;
                dyn_p_order;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]

    #----------------------------------------

    algebraic_generic_model_kwd_para =
        (;loc_load_exist,
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         pf_generic_gens_para,
         Ynet_wt_nodes_idx_wt_adjacent_nodes )


    algebraic_generic_model_sol_kwd_para =
        (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         flat_vh_flat_θh_id_iq_u0,
         dyn_pf_solver,
         algebraic_generic_model_kwd_para)


    algebraic_generic_model_wt_fault_kwd_para =
        (;loc_load_exist,                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
         dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
         
         on_fault_net_para,

         clear_fault_selection_list,
         no_lines_fault,
         no_cleared_lines_fault,
         list_fault_resistance,
         with_faults,
         cleared_selected_lines_faults_net_para )
    
    algebraic_generic_model_wt_fault_sol_kwd_para =
        (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         flat_vh_flat_θh_id_iq_u0,
         dyn_pf_solver,
         algebraic_generic_model_wt_fault_kwd_para,
         cleared_selected_lines_faults_net_para
         )


    dynamic_status_paras =
        (;system_status,
         dyn_pf_sol_u,
         dyn_pf_P_gens,
         dyn_pf_Q_gens,
         dyn_vh, dyn_θh,
         dyn_gens_vh,
         dyn_gens_θh,
         dyn_gens_id,
         dyn_gens_iq,
         dyn_gens_mag_E,
         dyn_gens_ang_E,
         post_dyn_PQ,
         dyn_Yred, dyn_Yint,
         
         dyn_plants_generic_states_init_wt_ref,
         dyn_plants_states_init,
         dyn_plants_refs,
         
         dyn_nt_vec_per_paras,
         dyn_vec_vec_per_paras,
         
         dyn_ω_ref, dyn_v_ref, dyn_p_order,
         dyn_gens_i_d, dyn_gens_i_q,
         
         flat_vh_flat_θh_id_iq_u0,
         flat_vh_flat_θh_id_iq_vfh_θfh,
         u0_model_states_init,
         u0_model_states_init_wt_fault_algeb_init,

         dyn_gens_δ,
         dyn_gens_ed_dash,
         dyn_gens_eq_dash,

         δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,
         ωref0_vref0_porder0_id_iq_vh,
         ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,          
         
         algebraic_generic_model_sol_kwd_para,
         algebraic_generic_model_wt_fault_sol_kwd_para,
         Pg_Qg_Png_Qng_Pll_Qll_Idx,
         Png_Qng_Pll_Qll_Idx,
         Pg_Png_Qng_Idx,

         edges_r_x_b_ratio_angle_idx,
         dyn_pf_flat_vh_flat_θh_Idx )
 
    # (; static_prefault_paras, dynamic_status_paras )

    # return system_status == :pre_fault_state ?
    #     (;static_prefault_paras, ) : (
    #         ;dynamic_status_paras, )


    return  (; dynamic_status_paras, )

end

# ---------------------------------------------------


"""
    get_ntuple_status_steady_state_data(
        ;<keyword arguments>)

Returns namedtuples of data for network status in `list_network_status`. The network status in the list are `:pre_fault_state`, `:fault_state`, `:post_fault_state`.


# Arguments

- `net_data_by_components_file`: the network data file
- `components_libs_dir`: the components library folder
- `basekV`: the base voltage
- `use_pu_in_PQ`: the boolean variable that determines if PQ should be in pu.
- `line_data_in_pu`: the boolean variable that informs if line data are in pu,
`pf_alg`: power flow solver
`abstol`: the absolute error tolerance
`reltol`: the relative error tolerance
`on_fault_time`: the on fault time
`clear_fault_time`: the clear fault time
`list_fault_point_from_node_a::Vector{Float64}=[0.3]`: the list containing a ratio of the fault point from the source (from) orientation of lines
`list_fault_resistance::Vector{Float64} = [0.001]`: the list containing fault resistances of each fault in the network.
`list_no_line_circuit::Vector{Float64} = [4]`: the list containing the number of circuits per faulted lines
`list_edges_to_have_fault::Vector{Int64} = [2]`: the list containing indices of lines that should have a fault.  
`clear_fault_selection_list::Vector{Int64} = [1]`: the list containing indices of faulted lines to be cleared in `list_edges_to_have_fault`.
`with_faults::Bool=false`: a legacy boolean variable.
`timespan`: the simulation time.
`list_network_status`: the list of network status.


"""
function get_ntuple_status_steady_state_data(
    ;with_faults =
        false,

    net_data_by_components_file =
        nothing,
    components_libs_dir =
        nothing,
    
    timespan         = 10.0,
    on_fault_time    = 9.0,
    clear_fault_time = 9.001,
    
    list_fault_point_from_node_a = [0.3],
    list_fault_resistance        = [0.001],
    list_no_line_circuit         =  [4],

    list_edges_to_have_fault     = [ 2 ],
    clear_fault_selection_list   = [1],
    
    basekV = 1.0,    
    use_pu_in_PQ =
        true,
    line_data_in_pu =
        true,
    list_network_status =
        nothing )

    #--------------------------------------
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing )

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end


    #--------------------------------------

    if list_network_status == nothing 

        list_network_status = 
            [:pre_fault_state,
             :fault_state,
             :post_fault_state]
    end

    #-------------------------------
    #-------------------------------

    dict_status_steady_state_data =
        Dict{Symbol, NamedTuple }()

    for a_system_status in list_network_status

        dict_status_steady_state_data[a_system_status] =
            get_a_status_steady_state_data(
                a_system_status;
                with_faults =
                    true,
                net_data_by_components_file =
                    net_data_by_components_file,
                components_libs_dir =
                    components_libs_dir,

                timespan         = timespan,
                on_fault_time    = on_fault_time,
                clear_fault_time = clear_fault_time,    

                list_fault_point_from_node_a =
                    list_fault_point_from_node_a,
                list_fault_resistance =
                    list_fault_resistance,
                list_no_line_circuit =
                    list_no_line_circuit,

                list_edges_to_have_fault =
                    list_edges_to_have_fault,
                clear_fault_selection_list =
                    clear_fault_selection_list,

                basekV = basekV,    
                use_pu_in_PQ =
                    use_pu_in_PQ,
                line_data_in_pu =
                    line_data_in_pu)

    end
    
    return NamedTupleTools.namedtuple(
            dict_status_steady_state_data)
        
end

# ------------------------------------------------------
# ------------------------------------------------------
