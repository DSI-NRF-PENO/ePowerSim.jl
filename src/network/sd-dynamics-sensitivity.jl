# (c) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.4:  pg 135 - 136, 6.121 - 6.123

# ####################################################

# using Pkg

# QuickDynamicTest=joinpath(@__DIR__,"..","..")
# cd(QuickDynamicTest)
# Pkg.activate(QuickDynamicTest)

#-------------------------------------------------------
#-------------------------------------------------------
# simulation flags and aux functions
#-------------------------------------------------------
#-------------------------------------------------------
 

"""

:hybrid_pf          = :nodes_fun_idx => 4

:node_pf            = :nodes_fun_idx => 4 

:global_pf          = :nodes_fun_idx => 3 

:network_current_pf = :nodes_fun_idx => 2 

:initial_pf         = :nodes_fun_idx => 1 

"""
function sd_dynamics_simulation(
    ; dynamics_case  = case_fun,
    nodes_fun_idx    = sim_nodes_fun_idx,
    sim_system_dynamics  = sim_system_dynamics,
    sim_nodes_edges_dynamics! = sim_nodes_edges_dynamics!,
    
    with_hybrid_pf = false,
    with_node_pf   = false,
    with_global_pf = false,
    with_network_current_pf = false,
    with_initial_pf = false,
    with_cb         = false,
    fixed_timestep  = false,
    timespan        = 10.0,
    dt              = 0.01,
    alg             = Rodas4(),
    alg_name        = alg_name,
    abstol          = 1e-12,
    reltol          = 1e-12 )

    #----------------------------------------------------    
    #----------------------------------------------------
    # Print flags and simulation functions
    #----------------------------------------------------
    #----------------------------------------------------

    name_case_fun = String(
        nameof( dynamics_case ) )
    
    name_sim_system_dynamics_fun = String(
        nameof( sim_system_dynamics ) )
    
    name_sim_nodes_edges_dynamics! = String(
        nameof( sim_nodes_edges_dynamics! ) )

    println("Simulating case  $(name_case_fun) ")

    println("nodes_fun_idx is $(nodes_fun_idx) ")
    
    println("system_dynamics_fun is $(name_sim_system_dynamics_fun) ")

    println("nodes_edges_dynamics! is $( name_sim_nodes_edges_dynamics! ) ")

    #----------------------------------------------------    
    
    tup_flags = (
        ; fixed_timestep,
        with_cb,
        with_initial_pf,
        with_network_current_pf,
        with_node_pf,
        with_global_pf,
        with_hybrid_pf)

    print_flags( tup_flags )
    
    #----------------------------------------------------
    # Net instatitiation
    #----------------------------------------------------

    # dynamics_case = case_IEEE_9_Bus_sauer_dynamic_plant_SM_idq_I_t2_system_matrix
    
    netd  = NetworkData( dynamics_case()... )

    #----------------------------------------------------

    pf_sys_param_sys_views_sys_init = get_pf_sys_param_sys_views_sys_init(
        netd;
        maxiter      = 40,
        ftol         = 1000*eps(),
        xtol         = 1000*eps(),
        init_pf      = true,
        with_δ_ed_eq = false )

    # ---------------------------------------------------

    # nodes_f_t, nodes_cb_sw, state, global_pf_param, named_tup_pf_result = pf_sys_param_sys_views_sys_init

    nodes_f_t       = pf_sys_param_sys_views_sys_init.nodes_f_t
    nodes_cb_sw     = pf_sys_param_sys_views_sys_init.nodes_cb_sw
    state           = pf_sys_param_sys_views_sys_init.state
    global_pf_param = pf_sys_param_sys_views_sys_init.global_pf_param

    # ---------------------------------------------------

    pf_net_param, _, _ = global_pf_param

    pf_net, _, _, _, _, pf_and_dyn_idx_and_Idx  , pf_net_misc = pf_net_param

    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net
    
    # ---------------------------------------------------
    
    nodes_incident_edges = get_nodes_incident_edges_by_orientations( edges_orientation )

    nodes_incident_edges_Ybr_cal = [ edges_Ybr_cal[a_node_incident_edges] for a_node_incident_edges in nodes_incident_edges]
    
    nodes_incident_edges_orientation = [ edges_orientation[a_node_incident_edges]  for a_node_incident_edges in nodes_incident_edges]
    
     nodes_pf_dyn_param = ( nodes_node_idx_and_incident_edges_other_node_idx, nodes_incident_edges_Ybr_cal, nodes_incident_edges_orientation )
    
    # ---------------------------------------------------
    # Testing
    # ---------------------------------------------------
    
    nodes_fun_params = (nodes_cb_sw,
                        nodes_f_t,
                        netd.nodes_param)

    edges_fun_params = netd.edges_param 

    #--------------------------------------------
    
    nodes_fun_idx = nodes_fun_idx   

    edges_fun_idx = 1    
    
    nodes_and_edges_p = create_params_views_and_Idx(
        netd, nodes_fun_params, edges_fun_params;
        nodes_fun_idx = nodes_fun_idx,
        edges_fun_idx = edges_fun_idx )

    #--------------------------------------------

    selected_node_fun = String(
        nameof(
            nodes_and_edges_p[1][1][1]))
    
    println("selected node function is $(selected_node_fun)")

    println("The solver being used is: $(alg_name)")
    #--------------------------------------------

    if with_hybrid_pf == true


        global_pf_options = (; maxiter = 40, ftol=1000*eps() , xtol=1000*eps(), with_δ_ed_eq = true )
        
        counter_array  = [1]
        
        stateDiffCache = similar( state )

        hybrid_pf_dyn_param = (
            global_pf_param,
            stateDiffCache,
            state,
            counter_array,
            global_pf_options,
            nodes_pf_dyn_param  )

        para  = (
            netd,
            nodes_and_edges_p,
            hybrid_pf_dyn_param )
        
    elseif with_node_pf == true 
        
        para  = (
            netd,
            nodes_and_edges_p,
            nodes_pf_dyn_param )

    elseif with_global_pf == true 

        global_pf_options = (
            ; maxiter = 40,
            ftol=1000*eps() ,
            xtol=1000*eps(),
            with_δ_ed_eq = true )
        
        counter_array  = [1]
        
        stateDiffCache = similar( state )

        para  = (
            netd,
            nodes_and_edges_p,
            global_pf_param,
            stateDiffCache,
            state,
            counter_array,
            global_pf_options )

    elseif with_network_current_pf == true
        
        para  = ( netd,
                  nodes_and_edges_p )
        
    elseif with_initial_pf == true
        
        para  = ( netd,
                  nodes_and_edges_p )
    else
        nothing
        
    end

    #--------------------------------------------    
    
    timespan  = ( 0.0, timespan )

    sd = sim_system_dynamics(
        sim_nodes_edges_dynamics!;
        mass_matrix = get_mass_matrix( netd ),
        syms = get_network_vars_labels( netd ) )
    
    prob  = ODEProblem( sd, state, timespan, para )

    #--------------------------------------------
    
    cb  = fun_make_state_callbacks(
        collect(values(netd.nodes)))

    #--------------------------------------------

    """
    
    QNDF2(), Rodas4(), Rodas5(), TRBDF2(),

    KenCarp4(), QNDF(),
    
    ImplicitMidpoint(), ImplicitMidpoint(autodiff=false),

    """

    if fixed_timestep == true

        if with_cb == true
        
            return DifferentialEquations.solve(
                prob,
                alg,
                callback = cb,
                abstol = abstol,
                reltol = reltol,
                dt = dt )
            
        else
            
            return DifferentialEquations.solve(
                prob,
                alg,
                abstol = abstol,
                reltol = reltol,
                dt = dt )
            
        end
        
    else
        
        if with_cb == true
            
            return DifferentialEquations.solve(
                prob,
                alg,
                callback = cb,
                abstol = abstol,
                reltol = reltol)
            
        else
            return DifferentialEquations.solve(
                prob,
                alg,
                abstol = abstol,
                reltol = reltol)
            
        end
        
    end
    
end


#-------------------------------------------------------
# Performance check
#-------------------------------------------------------

function performance_test_sd_dynamics_simulation(
    with_hybrid_pf = false,
    with_node_pf   = false,
    with_global_pf = false,
    with_network_current_pf = false,
    with_initial_pf = false,
    with_cb         = false,
    fixed_timestep  = false,
    timespan        = 10.0,
    dt              = 0.01,
    alg             = Rodas4(),
    alg_name        = alg_name,
    abstol          = 1e-12,
    reltol          = 1e-12 )

    #----------------------------------------------------    
    #----------------------------------------------------
    # Print flags and simulation functions
    #----------------------------------------------------
    #----------------------------------------------------
    
    #----------------------------------------------------
    # Net instatitiation
    #----------------------------------------------------

    dynamics_case = case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix
    
    netd  = NetworkData( dynamics_case()... )

    #----------------------------------------------------

    pf_sys_param_sys_views_sys_init = get_pf_sys_param_sys_views_sys_init(
        netd;
        maxiter      = 40,
        ftol         = 1000*eps(),
        xtol         = 1000*eps(),
        init_pf      = true,
        with_δ_ed_eq = false )

    # ---------------------------------------------------

    # nodes_f_t, nodes_cb_sw, state, global_pf_param, named_tup_pf_result = pf_sys_param_sys_views_sys_init

    nodes_f_t       = pf_sys_param_sys_views_sys_init.nodes_f_t
    nodes_cb_sw     = pf_sys_param_sys_views_sys_init.nodes_cb_sw
    state           = pf_sys_param_sys_views_sys_init.state
    global_pf_param = pf_sys_param_sys_views_sys_init.global_pf_param

    # named_tup_pf_result = pf_sys_param_sys_views_sys_init.named_tup_pf_result

    # ---------------------------------------------------

    pf_net_param, _, _ = global_pf_param

    pf_net, _, _, _, _, pf_and_dyn_idx_and_Idx  , pf_net_misc = pf_net_param

    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net
    
    # ---------------------------------------------------
    
    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )

    nodes_incident_edges_Ybr_cal = [
        edges_Ybr_cal[a_node_incident_edges]
        for a_node_incident_edges in
            nodes_incident_edges]
    
    nodes_incident_edges_orientation = [
        edges_orientation[a_node_incident_edges]
        for a_node_incident_edges in
            nodes_incident_edges]
    
    nodes_pf_dyn_param = (
        nodes_node_idx_and_incident_edges_other_node_idx,
        nodes_incident_edges_Ybr_cal,
        nodes_incident_edges_orientation )
    
    # ---------------------------------------------------
    # Testing
    # ---------------------------------------------------
    
    nodes_fun_params = (nodes_cb_sw,
                        nodes_f_t,
                        netd.nodes_param)

    edges_fun_params = netd.edges_param 

    #--------------------------------------------

    u0 = state
    du = similar(u0)
    t  = 5
    
    #--------------------------------------------
    
    in_pf, net_cur_pf, glo_pf,nod_pf, hyb_pf= (1,2,3,4,5)
    
    #--------------------------------------------
    # hybrid_pf options
    #--------------------------------------------
    
    nodes_fun_idx = hyb_pf

    edges_fun_idx = 1
    
    #--------------------------------------------
    
    nodes_and_edges_p = create_params_views_and_Idx(
        netd, nodes_fun_params, edges_fun_params;
        nodes_fun_idx = nodes_fun_idx,
        edges_fun_idx = edges_fun_idx )

    #--------------------------------------------
    
    global_pf_options =
        (; maxiter = 40, ftol=1000*eps() ,
         xtol=1000*eps(), with_δ_ed_eq = true )

    counter_array  = [1]

    stateDiffCache = similar( state )

    hybrid_pf_dyn_param = (
        global_pf_param,
        stateDiffCache,
        state,
        counter_array,
        global_pf_options,
        nodes_pf_dyn_param  )

    hybrid_para  = (
        netd,
        nodes_and_edges_p,
        hybrid_pf_dyn_param )


    @time hybrid_pf_nodes_edges_dynamics!(
        du, u0, hybrid_para, t)
    
    @code_warntype  hybrid_pf_nodes_edges_dynamics!(
        du, u0, hybrid_para, t)

    #--------------------------------------------
    # node_pf options
    #--------------------------------------------

    
    nodes_fun_idx = nod_pf

    edges_fun_idx = 1
    
    #--------------------------------------------
    
    nodes_and_edges_p = create_params_views_and_Idx(
        netd, nodes_fun_params, edges_fun_params;
        nodes_fun_idx = nodes_fun_idx,
        edges_fun_idx = edges_fun_idx )

    #--------------------------------------------
    
    node_pf_para  = (
            netd,
            nodes_and_edges_p,
            nodes_pf_dyn_param )

    @time node_pf_nodes_edges_dynamics!(
        du, u0, node_pf_para, t)
    
    @code_warntype  node_pf_nodes_edges_dynamics!(
        du, u0, node_pf_para, t)
    
    #--------------------------------------------
    # global_pf options
    #--------------------------------------------
    
    nodes_fun_idx = glo_pf

    edges_fun_idx = 1
    
    #--------------------------------------------
    
    nodes_and_edges_p = create_params_views_and_Idx(
        netd, nodes_fun_params, edges_fun_params;
        nodes_fun_idx = nodes_fun_idx,
        edges_fun_idx = edges_fun_idx )

    #--------------------------------------------

    global_pf_options = (
        ; maxiter = 40,
        ftol=1000*eps() ,
        xtol=1000*eps(),
        with_δ_ed_eq = true )

    counter_array  = [1]

    stateDiffCache = similar( state )

    global_pf_para  = (
        netd,
        nodes_and_edges_p,
        global_pf_param,
        stateDiffCache,
        state,
        counter_array,
        global_pf_options )

    @time global_pf_nodes_edges_dynamics!(
        du, u0, global_pf_para, t)
    
    @code_warntype  global_pf_nodes_edges_dynamics!(
        du, u0, global_pf_para, t)
    
    #--------------------------------------------
    # network_current_pf options
    #--------------------------------------------
    
    nodes_fun_idx = net_cur_pf

    edges_fun_idx = 1
    
    #--------------------------------------------
    
    nodes_and_edges_p = create_params_views_and_Idx(
        netd, nodes_fun_params, edges_fun_params;
        nodes_fun_idx = nodes_fun_idx,
        edges_fun_idx = edges_fun_idx )

    #--------------------------------------------
        
    network_current_pf_para  = (
        netd, nodes_and_edges_p )

    @time network_current_nodes_edges_dynamics!(
        du, u0, network_current_pf_para, t)
    
    
    @code_warntype  network_current_nodes_edges_dynamics!(
        du, u0, network_current_pf_para, t)
    
    #--------------------------------------------
    # initial_pf options
    #--------------------------------------------
    
    nodes_fun_idx = in_pf

    edges_fun_idx = 1
    
    #--------------------------------------------
    
    nodes_and_edges_p = create_params_views_and_Idx(
        netd, nodes_fun_params, edges_fun_params;
        nodes_fun_idx = nodes_fun_idx,
        edges_fun_idx = edges_fun_idx )

    #--------------------------------------------
        
    initial_pf_para  = ( netd, nodes_and_edges_p)

    @time initial_pf_nodes_edges_dynamics!(
        du, u0, initial_pf_para, t)
    
    @code_warntype  initial_pf_nodes_edges_dynamics!(
        du, u0, initial_pf_para, t)
    
    #--------------------------------------------    
    #--------------------------------------------    
    
end

#-------------------------------------------------------
#-------------------------------------------------------


function create_namedtup_sim_model_types()

    list_nodes_fun_idx = [1, 2, 3, 4, 5]

    list_system_dynamics_func = [
        initial_pf_system_dynamics,
        network_current_system_dynamics,
        global_pf_system_dynamics,
        node_pf_system_dynamics,
        hybrid_pf_system_dynamics ] 

    list_nodes_edges_dynamics_func! = [
        initial_pf_nodes_edges_dynamics!,
        network_current_nodes_edges_dynamics!,
        global_pf_nodes_edges_dynamics!,
        node_pf_nodes_edges_dynamics!,
        hybrid_pf_nodes_edges_dynamics! ]

    list_model_type  = [
        :initial_pf,
        :network_current_pf,
        :global_pf,
        :node_pf,
        :hybrid_pf ]

    return namedtuple(
        [ a_type => namedtuple(
            [:nodes_fun_idx => a_node_fun_idx,
             :system_dynamics=> a_sys_dyn_fun,
             :nodes_edges_dynamics! =>
                 a_nodes_edges_dyn_fun ])
          for(
              a_type,
              a_node_fun_idx,
              a_sys_dyn_fun,
              a_nodes_edges_dyn_fun)  in  zip(
                  list_model_type,
                  list_nodes_fun_idx,
                  list_system_dynamics_func,
                  list_nodes_edges_dynamics_func!) ])
  
end


function create_simulation_model_types_and_flags( )

    #; ntup_sim_model_types = ntup_sim_model_types
    
    list_model_type_flags = [
        (; with_initial_pf = true,),
        (; with_network_current_pf = true,),
        (; with_global_pf = true,),
        (; with_node_pf = true,),
        (; with_hybrid_pf = true, )
        ]

    list_model_type  = [
        :initial_pf,
        :network_current_pf,
        :global_pf,
        :node_pf,
        :hybrid_pf ]

    list_sym_model_type_flags  = [
        :initial_pf_model_flags,
        :network_current_pf_model_flags,
        :global_pf_model_flags,
        :node_pf_model_flags,
        :hybrid_pf_model_flags ]


    ntup_sim_model_types =
        create_namedtup_sim_model_types()
    
    sim_model_types_and_flags = [
        (
            ; nodes_fun_idx  = getproperty(
                ntup_sim_model_types,
                model_type).nodes_fun_idx ,
            sim_system_dynamics = getproperty(
                ntup_sim_model_types,
                model_type).system_dynamics, 
            sim_nodes_edges_dynamics! = getproperty(
                ntup_sim_model_types,
                model_type).nodes_edges_dynamics!,
            model_type_flags...  )
        for ( model_type, model_type_flags ) in zip(
            list_model_type, list_model_type_flags ) ]
    
    return (sim_model_types_and_flags,
            list_sym_model_type_flags)

end


function create_simulation_namedtup_for_model_types_and_flags( )

    # (; sim_model_types_and_flags = sim_model_types_and_flags , list_sym_model_type_flags = list_sym_model_type_flags)

    list_model_type  = [
        :initial_pf,
        :network_current_pf,
        :global_pf,
        :node_pf,
        :hybrid_pf ]

    
    (sim_model_types_and_flags,
     list_sym_model_type_flags) =
         create_simulation_model_types_and_flags( )
    
    # sim_namedtup_for_model_types
    
    sim_namedtup_for_model_types =
        namedtuple([
            model_type => ( sim_model_type_and_flags ,
                            sym_model_type_flags  )
            for ( model_type,
                  sim_model_type_and_flags,
                  sym_model_type_flags ) in
                zip(list_model_type,
                    sim_model_types_and_flags,
                    list_sym_model_type_flags )])

    return  sim_namedtup_for_model_types


end

function create_a_simulation_model_type_and_flags(
    model_type, sim_namedtup_for_model_types )
    
    return getproperty(
        sim_namedtup_for_model_types,
        model_type )
    
end


"""
Test:

# t_sim_model_types_and_flags, t_list_sym_model_type_flags = create_simulation_model_types_and_flags( )


# t_sim_namedtup_for_model_types = create_simulation_namedtup_for_model_types_and_flags( )

initial_pf_sim_model_type_and_flags , initial_pf_sym_model_type_flags = create_a_simulation_model_type_and_flags( :initial_pf, t_sim_namedtup_for_model_types )

network_current_pf_sim_model_type_and_flags , network_current_pf_sym_model_type_flags = create_a_simulation_model_type_and_flags( :network_current_pf, t_sim_namedtup_for_model_types )

global_pf_sim_model_type_and_flags , global_pf_sym_model_type_flags = create_a_simulation_model_type_and_flags( :global_pf, t_sim_namedtup_for_model_types )

node_pf_sim_model_type_and_flags , node_pf_sym_model_type_flags = create_a_simulation_model_type_and_flags( :node_pf, t_sim_namedtup_for_model_types )

hybrid_pf_sim_model_type_and_flags , hybrid_pf_sym_model_type_flags = create_a_simulation_model_type_and_flags( :hybrid_pf, t_sim_namedtup_for_model_types )

"""

#-------------------------------------------------------
# a sinlge model simulation
#-------------------------------------------------------

function simulate_a_dynamics_simulation_model_type(
    case_fun;
    sym_model_type =  sim_sym_model_type,
    sim_namedtup_for_model_types =
        sim_namedtup_for_model_types,
    dt = 0.01,
    alg = ImplicitMidpoint(),
    alg_name = alg_name,
    with_cb = false,    
    fixed_timestep = true,
    timespan = 10.0,    
    plot_tspan = (0.0, 10),
    base_dir = Resultsbase_dir)

    #-----------------------------------------------
    
    case_name = String( nameof( case_fun ) )

    net_syms  = get_network_vars_labels(
        NetworkData( case_fun()... ) )

    (network_bus_names,
     non_gens_bus_names,
     gens_bus_names,
     Loads_bus_names,Trans_bus_names) =
         make_case_buses_names(
        ;case_fun = case_fun )

    (res_folder,
     csv_folder,
     fig_folder,
     sol_folder) =
         make_results_dir(
             ; base_dir=base_dir,
             res_dir = case_name )

    #-----------------------------------------------
        
    (namedtup_sim_model_type_and_flags,
     sym_model_type_flags ) =
         create_a_simulation_model_type_and_flags(
             sym_model_type,
             sim_namedtup_for_model_types )

    #-----------------------------------------------
    
    sol = sd_dynamics_simulation(
        ;dynamics_case = case_fun,
        namedtup_sim_model_type_and_flags...,
        with_cb        = with_cb,
        fixed_timestep = fixed_timestep,
        timespan       = timespan,
        dt             = dt,
        alg            = alg,
        alg_name       = alg_name)

    #-----------------------------------------------
    
    plt_ω_δ_ed_dash_eq_dash =
        make_a_plot_for_syms(
            sol ,
            [:ω, :δ, :ed_dash, :eq_dash],
            net_syms;
            tspan = plot_tspan )

    plt_gov = make_a_plot_for_syms(
        sol,
        [:xg1, :xg2, :τm_tilade],
        net_syms;
        tspan = plot_tspan )

    plt_avr = make_a_plot_for_syms(
        sol,
        [:vm, :vr1, :vr2, :vf_tilade],
        net_syms;
        tspan = plot_tspan )

    pt_gen =
        plot(
            plt_ω_δ_ed_dash_eq_dash;
            size = (1000, 500),
            lw = 2, xlabel = "t[s]" )

    pt_gov =
        plot( plt_gov;  size = (1000, 500),
              lw = 2, xlabel = "t[s]" )

    pt_avr =
        plot( plt_avr;  size = (1000, 500),
              lw = 2, xlabel = "t[s]" )
    
    plt_gens_u_mag =
        make_plot_of_buses_volt_mag(
            ; sol = sol,
            network_vars_labels = net_syms,
            nodes_name = gens_bus_names,
            vars = [:u_r, :u_i],
            tspan = plot_tspan,
            fmt = :png )

    plt_loads_u_mag =
        make_plot_of_buses_volt_mag(
        ; sol = sol,
        network_vars_labels = net_syms,
        nodes_name = Loads_bus_names,
        vars = [:u_r, :u_i],
        tspan = plot_tspan,
        fmt = :png )

    plt_trans_u_mag =
        make_plot_of_buses_volt_mag(
        ; sol = sol,
        network_vars_labels = net_syms,
        nodes_name = Trans_bus_names,
        vars = [:u_r, :u_i],
        tspan = plot_tspan,
        fmt = :png )

    plt_gens_u_ang  =
        make_plot_of_buses_volt_angle(
        ; sol = sol,
        network_vars_labels = net_syms,
        nodes_name = gens_bus_names,
        vars = [:u_r, :u_i],
        tspan = plot_tspan,
        fmt = :png )

    plt_loads_u_ang =
        make_plot_of_buses_volt_angle(
        ; sol = sol,
        network_vars_labels = net_syms,
        nodes_name = Loads_bus_names,
        vars = [:u_r, :u_i],
        tspan = plot_tspan,
        fmt = :png )

    plt_trans_u_ang =
        make_plot_of_buses_volt_angle(
        ; sol = sol,
        network_vars_labels = net_syms,
        nodes_name = Trans_bus_names,
        vars = [:u_r, :u_i],
        tspan = plot_tspan,
        fmt = :png )

    plt_list_volt_and_angle = [
        plt_gens_u_mag, plt_loads_u_mag,
        plt_trans_u_mag, plt_gens_u_ang,
        plt_loads_u_ang, plt_trans_u_ang ]

    plt_vh_θh = plot(
        plt_list_volt_and_angle... ;
        layout = (3,2), size = (1000, 500),
        lw = 3, xlabel = "t[s]" )

    
   list_sim_plots = [pt_gen, pt_gov,
                      pt_avr, plt_vh_θh ]
    
    plots_name    = ["pt_gen", "pt_gov",
                      "pt_avr", "plt_vh_θh" ]    
    
    suffix = :png
    
    for (a_sim_plot, a_plot_name) in
        zip(list_sim_plots, plots_name )
        
        savefig(
            a_sim_plot,
            joinpath(
                fig_folder,
                "$(alg_name)-$(String(sym_model_type))-$(a_plot_name).$(suffix)"))
        
    end

    println("Done with solver: $(alg_name), model type: $(String(sym_model_type)),  case: $(case_name) ")

    return nothing    
    
end


#-------------------------------------------------------


"""

Test driver for:

`simulate_a_dynamics_simulation_model_type`

    Test:    

    solvers:

        QNDF2(),Rodas4(), Rodas5(), TRBDF2(),
        KenCarp4(), QNDF(), ImplicitMidpoint(),
        ImplicitMidpoint(autodiff=false),

    model types are:

        :initial_pf, :network_current_pf,
         :global_pf , :node_pf, hybrid-pf
    
    list_cases_fun = [
        case_IEEE_5_Bus_dynamic_plant_SM_cb_v6_P,
        case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix,
        case_IEEE_14_Bus_dynamic_plants_v6_P_rscad,
        case_IEEE_14_Bus_dynamic_plants_v6_P_t2_millano,
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P]
     
    list_sym_model_type = [ :initial_pf,
                            :network_current_pf,
                            :node_pf,
                            :global_pf,
                            :hybrid_pf  ]
    
    list_alg = [ Rodas4(), Rodas5(),
        ImplicitMidpoint(),
        ImplicitMidpoint(autodiff=false)
        TRBDF2(),  KenCarp4()

    bool_fixed_timestep = [
        false, false,
        true,
        true,
        false,false ]

"""
function driver_simulate_a_dynamics_simulation_model_type( )
    
    dt              = 0.01    
    with_cb         = true    
    timespan        = 10.0    
    plot_tspan      = (0.0, timespan )

    
    Resultsbase_dir =
        joinpath(@__DIR__,"..","..","..",
                 "Simulation-Results",
                 "Quick-Dynamic-Test",
                 "model3-dynamics")

    sim_namedtup_for_model_types =
        create_simulation_namedtup_for_model_types_and_flags( )

   
    # list_cases_fun = [
    #     case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix,
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P ]
         
    # list_sym_model_type = [
    #     :initial_pf,
    #     :network_current_pf,
    #     :node_pf,
    #     :global_pf,
    #     :hybrid_pf  ]    
    
    # list_alg = [
    #     ImplicitMidpoint(),
    #     ImplicitMidpoint(autodiff=false),
    #     Rodas4(),
    #     Rodas5() ]

    # list_alg_names = [
    #     "ImplicitMidpoint",
    #     "ImplicitMidpoint_autodiff_false",
    #     "Rodas4",
    #     "Rodas5" ]

    # bool_fixed_timestep = [        
    #     true,
    #     true,
    #     false,
    #     false ]

    # case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix
    
    # case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    # case_IEEE_5_Bus_dynamic_plant_SM_cb_v6_P,
    
    # case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
    
    # case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,    

    list_cases_fun =
        [ case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix, case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P  ]
    
    
    list_sym_model_type = [ :initial_pf  ]    
    
    list_alg = [ Rodas4() ]

    # ImplicitMidpoint(),
    
    list_alg_names = [ "Rodas4" ]

    bool_fixed_timestep = [ false ]
    
    
    # Threads.@threads 
     for a_case_fun in list_cases_fun
        for a_model_type in list_sym_model_type
            for (alg, alg_name, fixed_timestep ) in
                zip(list_alg, list_alg_names,
                    bool_fixed_timestep )
                try
                    simulate_a_dynamics_simulation_model_type(
                        a_case_fun ;
                        sym_model_type = a_model_type,
                        sim_namedtup_for_model_types =
                            sim_namedtup_for_model_types,
                        dt = dt,
                        alg = alg,
                        alg_name = alg_name,
                        with_cb = with_cb,
                        fixed_timestep = fixed_timestep,
                        timespan = timespan,
                        plot_tspan = plot_tspan,
                        base_dir = Resultsbase_dir)
                catch sim_error
                    println("An error occurred: ",
                            sim_error) 
                end
                
            end
        end
    end
    
    
    println(  " I am done with simulation of dynamic_model_types")

    return nothing
end

"""

driver_simulate_a_dynamics_simulation_model_type( )

"""
#-----------------------------------------------------
######################################################
#-----------------------------------------------------

function make_plots_for_industrial_model(
    ; case_fun = dynamics_case,
    sim_model_type = sim_model_type,
    sim_sol = sim_sol,
    para_net_names_labels_syms =
        para_net_names_labels_syms,
    tspan = plot_tspan,
    base_dir = nothing,
    algr_name = algr_name )

    (net_class_names,
     net_states_and_var_labels,
     industrial_model_sym) =
         para_net_names_labels_syms

    (network_bus_names,
     non_gens_bus_names,
     gens_bus_names,
     Loads_bus_names,
     Trans_bus_names) =
         net_class_names

    (net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_stab_states_label,
     gens_nodes_algebraic_and_states_labels) =
         net_states_and_var_labels


    sim_plt_ω_δ_ed_dash_eq_dash =
        make_a_plot_for_syms(
        sim_sol ,
        [:ω, :δ, :ed_dash, :eq_dash],
        industrial_model_sym;
        tspan = tspan )

    
    sim_plt_gov = make_a_plot_for_syms(
        sim_sol,
        [:xg1, :xg2, :τm_tilade],
        industrial_model_sym;
        tspan =  tspan )

    
    sim_plt_avr = make_a_plot_for_syms(
        sim_sol,
        [:vm, :vr1, :vr2, :vf_tilade],
        industrial_model_sym;
        tspan =  tspan )
    
      
    sim_gens_volt_mag = make_plot_of_buses_volt_mag(
        ; sol = sim_sol,
        network_vars_labels = industrial_model_sym,
        nodes_name = gens_bus_names,
        vars = [:u_r, :u_i],
        tspan = tspan,
        fmt = :png )

    
    sim_gens_volt_angle =
        make_plot_of_buses_volt_angle(
        ; sol = sim_sol,
            network_vars_labels =
                industrial_model_sym,
        nodes_name = gens_bus_names,
        vars = [:u_r, :u_i],
        tspan = tspan,
        fmt = :png )

    
    sim_non_gens_volt_mag =
        make_plot_of_buses_volt_mag(
        ; sol = sim_sol,
            network_vars_labels =
                industrial_model_sym,
        nodes_name = non_gens_bus_names,
        vars = [:u_r, :u_i],
        tspan = tspan,
        fmt = :png )
    
    sim_non_gens_volt_angle =
        make_plot_of_buses_volt_angle(
        ; sol = sim_sol,
            network_vars_labels =
                industrial_model_sym,
        nodes_name = non_gens_bus_names,
        vars = [:u_r, :u_i],
        tspan = tspan,
        fmt = :png )

    list_plts_volt_mag_ang = [
        sim_gens_volt_mag,
        sim_gens_volt_angle,
        sim_non_gens_volt_mag,
        sim_non_gens_volt_angle ]

    list_plts_gov_avr = [
        sim_plt_gov,
        sim_plt_avr]

    #-----------------------------------------------

        plts_volt_mag_ang = plot(
            list_plts_volt_mag_ang... ;
            layout = (2,2),
            size = (1000, 500),
            lw = 3,
            xlabel = "t[s]" )

        plts_gov_avr = plot(
            list_plts_gov_avr... ;
            layout = (2,1),
            size = (1000, 500),
            lw = 3,
            xlabel = "t[s]" )
    
    #-----------------------------------------------

    # folder

    sim_model_type =
        sim_model_type 
    
    case_name =
        String(nameof(case_fun))
    
    results_dir =
        joinpath(sim_model_type, case_name )

    if base_dir==nothing
        
        (res_folder,
         csv_folder,
         fig_folder,
         sol_folder) =
            make_results_dir(
                ; base_dir=@__DIR__,
                res_dir = results_dir )
        
    else
        
        (res_folder,
         csv_folder,
         fig_folder,
         sol_folder) =
            make_results_dir(
                ; base_dir=base_dir,
                res_dir = results_dir )
        
    end
    
    #-----------------------------------------------
    
    suffix = :png
        
    list_plots_sym = [:plt_ω_δ_ed_dash_eq_dash,
                      :volt_mag_ang,
                      :gov_avr]
    
    list_sim_plots = [
    sim_plt_ω_δ_ed_dash_eq_dash,
    plts_volt_mag_ang,
        plts_gov_avr]

    for ( a_sim_plot, a_plot_sym ) in
        zip( list_sim_plots, list_plots_sym )
        savefig(
            a_sim_plot,
            joinpath(
                fig_folder,
                "$(algr_name)-$(sim_model_type)-$(a_plot_sym).$(suffix)") )
    end
            
    return list_sim_plots
    
    
end

#-----------------------------------------------------
#-----------------------------------------------------

function get_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init(
    netd
    ; maxiter=40,
    ftol=1000*eps(),
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false )
    
    # -----------------------------------------------------
    # -----------------------------------------------------

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

    """
    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants =
        get_industrial_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(  gens_nodes_collection )

    nodes_ω_ed_dash_eq_dash_Idxs_in_plants =
        get_industrial_gens_nodes_ω_ed_dash_eq_dash_Idxs_in_plants( gens_nodes_collection  )

    """
    
    # -------------------------------------------------
    
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
    
    industrial_model_each_gen_nodes_pure_states_idx_in_state = get_industrial_gens_pure_states_indices_in_state(netd )
    
    industrial_model_each_gen_nodes_stab_states_idx_in_state = get_industrial_gens_stab_states_indices_in_state(netd )
    
    
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
    #---------------------------------------------------

    gen_nodes_ra_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list = [ :ra, :X_d_dash, :X_q_dash ] )

    gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_industrial_model_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list = [ :ra, :X_d, :X_q, :X_d_dash, :X_q_dash ] )
    
    # gen_nodes_ra_Xd_dash_Xq_dash_view = get_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list = [ :ra, :X_d_dash, :X_q_dash ], gens_view_only = true )

    gen_nodes_dyn_param_view = get_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list = [ :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q, :X_d_dash, :X_q_dash, :T_d_dash, :T_q_dash ], gens_view_only = true )
    
     gen_nodes_sub_dyn_param_view = get_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list = [ :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q, :X_d_dash, :X_q_dash, :X_d_2dash, :X_q_2dash, :T_d_dash, :T_q_dash, :T_d_2dash, :T_q_2dash ], gens_view_only = true )
   
    # -----------------------------------------------------
    # -----------------------------------------------------

    # pf_net_param = get_powerflow_net_parameters( netd )

    pf_net_param =
        get_industrial_model_powerflow_net_parameters(
            netd )
    
    pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx


    #----------------------------------------------------    
    #----------------------------------------------------

    state =
        zeros(
            length(
                generate_industrial_model_sym(
                    ; nodes = netd.nodes )  ) )

    #----------------------------------------------------   
    
    gen_nodes_δ_ω_ed_dash_eq_dash_views =
        get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash_views( state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )
    
    industrial_model_pure_states_view_in_state =
        get_industrial_pure_states_view_in_state(
            state, industrial_model_pure_states_Idx )
    
    #------------------------------------------------- 
    #-------------------------------------------------

    # δ_ed_dash_eq_dash_view = get_industrial_nodes_δ_ed_dash_eq_dash_view(state, netd )

    # # get_nodes_δ_ed_dash_eq_dash_view( state, netd )

    # δ_ω_ed_dash_eq_dash_view = get_industrial_nodes_δ_ω_ed_dash_eq_dash_view(state, netd )

    # # get_nodes_δ_ω_ed_dash_eq_dash_view( state, netd )

    # ω_ed_dash_eq_dash_view = get_industrial_nodes_ω_ed_dash_eq_dash_view(state, netd )
    
    #---------------------------------------------------
    # --------------------------------------------------

    state_view = view(state, 1:length(state))

    # ----------------------------------------------------
    
    # for pf flat start, ur = vh = 1,  ui = θh = 0

    # remove init_pf = true
        
    
    if init_pf == true 

        state_view[ur_idx] .=
            ones(  length( ur_IDX ))
        state_view[ui_idx] .=
            zeros( length( ui_IDX ))

    end
    
    # ----------------------------------------------------

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]
        
    #----------------------------------------------------
    
    nodes_u_view  =
        [ view(state, nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ]

    # ---------------------------------------------------
    
    nodes_pf_U_view  =
        [ view(pf_state , nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ] 
    
    # ----------------------------------------------------
    
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
    
    # ----------------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))
    
    # ----------------------------------------------------
    
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
    
    #----------------------------------------------------

    global_pf_views =
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view )

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
          sd_pf_views,
          mismatch )
    
    branches_name  = collect(keys( netd.edges ))
    
    nodes_name = collect(keys( netd.nodes ))    

    #----------------------------------------------------
    #----------------------------------------------------

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

    #----------------------------------------------------

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init =
        named_tup_pf_result.branch_dict_init
    
    #---------------------------------------------------
    
    #---------------------------------------------------

    state .=
        industrial_model_init_operationpoint(
            netd, bus_dict_init; pure = :pure )

    #---------------------------------------------------- 

    # gen_nodes_δ_ω_ed_dash_eq_dash_views .= get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash( state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!( gen_nodes_δ_ω_ed_dash_eq_dash_views, state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    
    #----------------------------------------------------
    
    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)
    
    #----------------------------------------------------

    # nodes_f_t  = external_get_nodes_or_edges_f_t(
    #     netd.nodes, bus_dict_init )

    #----------------------------------------------------
    #---------------------------------------------------
    
 
    # gens_vh_θh_post_pf = get_gens_vh_θh_post_pf( gens_nodes_collection , bus_dict_init )
        
    # gens_ur_ui_post_pf = get_gens_ur_ui_post_pf( gens_nodes_collection, bus_dict_init )

    gens_vh_θh =
        get_gens_vh_θh(
            nodes_pf_U_view, gens_idx )

    gens_vh_θh_view =
        @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    """
    gens_ur_ui = get_gens_ur_ui(nodes_pf_U_view, gens_idx )

    gens_ur_ui_view = @view gens_ur_ui[ 1:length(gens_ur_ui ) ]

    """
   #---------------------------------------------------
   #----------------------------------------------------

    
    # gen_uh  = (named_tup_pf_result.Vbus)[ gen_idx ]

    idq_wt_pad_view[gens_idx] .=
        [ industrial_model_get_dynamic_idq_θ_π_vhθh(
            vh_θh, δ_ω_ed_eq,
            ra_Xd_dash_Xq_dash )
          for (vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash) in
              zip( gens_vh_θh_view ,
                   gen_nodes_δ_ω_ed_dash_eq_dash_views,
                   ra_Xd_dash_Xq_dash_view[gens_idx] ) ]
  

    #---------------------------------------------------

    # gens_nodes_ωs_τm_v_ref = get_gens_nodes_ωs_τm_v_ref( gens_nodes_collection, bus_dict_init )

    # gens_nodes_ωs_τm_v_ref_view = view( gens_nodes_ωs_τm_v_ref, 1:length( gens_nodes_ωs_τm_v_ref ) )

    #---------------------------------------------------

    gens_nodes_ωs_τm_vref_porder =
        get_gens_nodes_ωs_τm_vref_porder(
            gens_nodes_collection, bus_dict_init )

    gens_nodes_ωs_τm_vref_porder_view =
        view(
            gens_nodes_ωs_τm_vref_porder,
            1:length(gens_nodes_ωs_τm_vref_porder) )
    
    #---------------------------------------------------

    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )

    #------------------------------------------------   
    
    """
     industrial_model_pure_states_in_state = get_industrial_gen_nodes_pure_states_from_state(state, industrial_model_each_gen_nodes_pure_states_idx_in_state )

    gen_nodes_stab_states_from_state  = get_industrial_gen_nodes_stab_states_from_state(state, industrial_model_each_gen_nodes_stab_states_idx_in_state )

    #------------------------------------------------
   
    
    each_gen_nodes_pure_state_views = get_industrial_gen_nodes_pure_state_views( state, industrial_model_each_gen_nodes_pure_states_idx_in_state )

    each_gen_nodes_stab_state_views = get_industrial_gen_nodes_stab_state_views( state, industrial_model_each_gen_nodes_stab_states_idx_in_state )

    #------------------------------------------------
  

    each_gen_nodes_pure_states_from_state = get_industrial_each_gen_nodes_pure_states_from_state( state, industrial_model_each_gen_nodes_pure_states_idx_in_state )
    
    each_gen_nodes_stab_states_from_state = get_industrial_each_gen_nodes_stab_states_from_state( state, industrial_model_each_gen_nodes_stab_states_idx_in_state )

    """
    
    #------------------------------------------------
    #------------------------------------------------

    gens_dynamic_id_iq_pg_vh_by_vhθh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )

    """
    gens_dynamic_id_iq_pg_vh_by_ur_ui = get_gens_dynamic_id_iq_pg_vh_by_ur_ui( gens_ur_ui_post_pf, gen_nodes_δ_ω_ed_dash_eq_dash_views, gen_nodes_ra_Xd_dash_Xq_dash_view )

    """
    #------------------------------------------------  

    nodes_u_Idx_in_ranges =
        get_nodes_u_Idx_in_ranges(
            nodes_u_Idx )
    
    #------------------------------------------------
    
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
         industrial_model_each_gen_nodes_pure_states_idx_in_state, gens_nodes_collection )

    para_update_gen_Ax_aux =
        (; industrial_model_pure_states_view_in_state,
         industrial_model_each_gen_nodes_pure_states_idx_in_state, industrial_model_pure_states_Idx )

    industrial_model_para_aux_inputs =
        (; industrial_model_nodes_voltage_Idx,
         industrial_model_pure_states_Idx,
         nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         gen_nodes_ra_Xd_dash_Xq_dash_view,
         gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
         industrial_model_pure_states_view_in_state,
         gen_nodes_δ_ω_ed_dash_eq_dash_views,
         gens_vh_θh_view, nodes_pf_U_view )
    
    
    # industrial_model_pf_para = (; gens_dynamic_id_iq_pg_vh_by_vhθh, gens_nodes_ωs_τm_v_ref_view )

    industrial_model_pf_para =
        (; gens_dynamic_id_iq_pg_vh_by_vhθh,
         gens_nodes_ωs_τm_vref_porder_view )

    

    """ need by their views """
    
    # industrial_model_ωs_τm_vref_vhθh_idq =  (; gens_nodes_ωs_τm_v_ref, gens_vh_θh, idq_wt_pad  )

    industrial_model_ωs_τm_vref_vhθh_idq =
        (; gens_nodes_ωs_τm_vref_porder,
         gens_vh_θh, idq_wt_pad )
    

    industrial_model_dyn_pf_up_para =
        (; nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_ω_ed_dash_eq_dash_Idxs_in_state,
         nodes_δ_ed_dash_eq_dash_Idxs_in_state )
    
    #------------------------------------------------

    industrial_model_idq_pf_cal =
        (  idq_wt_pad_view,
           gens_idx )
    
    #------------------------------------------------

    return (; nodes_cb_sw,
            state,
            global_pf_param,
            named_tup_pf_result,
            industrial_model_misc_Idx,
            para_update_gen_Ax_aux,
            industrial_model_para_aux_inputs,
            industrial_model_pf_para,
            industrial_model_ωs_τm_vref_vhθh_idq ,
            industrial_model_dyn_pf_up_para,
            industrial_model_idq_pf_cal,
            gens_nodes_τm_vf )

    # return (; nodes_cb_sw, state, global_pf_param, named_tup_pf_result, industrial_model_misc_Idx, para_update_gen_Ax_aux, industrial_model_para_aux_inputs, industrial_model_pf_para, industrial_model_ωs_τm_vref_vhθh_idq , industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal, industrial_model_nodes_voltage_Idx, gens_nodes_τm_vf  )

    #---------------------------------------------------
    #---------------------------------------------------
    
end

#---------------------------------------------------

"""

Test driver for:

`driver_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init`


"""
function driver_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init(  )

     dynamics_case = case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix

    

   netd  = NetworkData( dynamics_case()... )

    nodes_cb_sw, state, global_pf_param, named_tup_pf_result, industrial_model_misc_Idx, para_update_gen_Ax_aux, industrial_model_para_aux_inputs, industrial_model_pf_para, industrial_model_ωs_τm_vref_vhθh_idq , industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal,  gens_nodes_τm_vf = get_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init( netd; maxiter=40, ftol=1000*eps() , xtol=1000*eps(), init_pf = true, with_δ_ed_eq = false )

    return nothing
    
end

"""

driver_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init( )

"""

#---------------------------------------------------
#---------------------------------------------------

function get_industrial_model_pf_param_views_and_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen  =false)


    dict_sys_to_industry = get_net_to_industrial_model_indices_dict( netd; no_control_device = only_gen  )

    pure_states_Idx_in_system, ur_ui_Idx_in_system, industrial_Idx, industrial_model_pure_states_Idx, industrial_model_ur_ui_Idx, pure_states_and_ur_ui_Idx_in_system, industrial_model_state_Idx, net_to_industrial_idx_conversion_dict = get_industrial_model_indices_and_conversion_dict( netd; no_control_device = only_gen  )

    # ---------------------------------------------------
    # ---------------------------------------------------    

    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state = get_industrial_δ_ω_ed_dash_eq_dash_Idx(netd; no_control_device = only_gen )

    nodes_ω_ed_dash_eq_dash_Idxs_in_state   = get_industrial_ω_ed_dash_eq_dash_Idx(netd; no_control_device = only_gen )

    nodes_δ_ed_dash_eq_dash_Idxs_in_state   =  get_industrial_δ_ed_dash_eq_dash_Idx(netd; no_control_device = only_gen )

    #---------------------------------------------------    

    industrial_model_each_gen_nodes_pure_states_idx_in_state = get_industrial_gens_pure_states_indices_in_state(netd; no_control_device = only_gen )

    industrial_model_each_gen_nodes_stab_states_idx_in_state = get_industrial_gens_stab_states_indices_in_state(netd; no_control_device = only_gen )

    #---------------------------------------------------

    gen_nodes_ra_Xd_dash_Xq_dash_view = get_industrial_model_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list = [ :ra, :X_d_dash, :X_q_dash ] )

     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view = get_industrial_model_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list = [ :ra, :X_d, :X_q, :X_d_dash, :X_q_dash ] )

    gen_nodes_dyn_param_view = get_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list = [ :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q, :X_d_dash, :X_q_dash, :T_d_dash, :T_q_dash ], gens_view_only = only_gen )

     gen_nodes_sub_dyn_param_view = get_gens_params_view_in_param_values( netd.nodes_param, netd.nodes; param_list = [ :D, :H, :ωs, :ra, :xℓ, :X_d, :X_q, :X_d_dash, :X_q_dash, :X_d_2dash, :X_q_2dash, :T_d_dash, :T_q_dash, :T_d_2dash, :T_q_2dash ], gens_view_only = only_gen )

    # -----------------------------------------------------
    # -----------------------------------------------------

    # pf_net_param = get_powerflow_net_parameters( netd )

    pf_net_param = get_industrial_model_powerflow_net_parameters( netd; no_control_device = only_gen )

    pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx


    #------------------------------------------------
    #------------------------------------------------

    state = zeros(length( generate_industrial_model_sym(; nodes = netd.nodes, no_control_device = only_gen )  ) )

    #------------------------------------------------   

    gen_nodes_δ_ω_ed_dash_eq_dash_views = get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash_views( state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    industrial_model_pure_states_view_in_state = get_industrial_pure_states_view_in_state(state, industrial_model_pure_states_Idx )

    #------------------------------------------------ 
    #------------------------------------------------  

    state_view = view(state, 1:length(state))

    # -----------------------------------------------

    # for pf flat start, ur = vh = 1,  ui = θh = 0

    # remove init_pf = true


    if init_pf == true 

        state_view[ur_idx] .= ones(  length( ur_IDX ))
        state_view[ui_idx] .= zeros( length( ui_IDX ))

    end

    # -----------------------------------------------

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]

    #------------------------------------------------

    nodes_u_view  = [ view(state, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx)) ]

    # ------------------------------------------------ 

    nodes_pf_U_view  = [ view(pf_state , nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx)) ] 

    # --------------------x-----------------------------

    uh_state = state_view[ur_ui_idx][ ur_IDX ] + im * state_view[ur_ui_idx][ ui_IDX ]

    x0_ur_ui = [state_view[ur_ui_idx][ ur_IDX ]...; state_view[ur_ui_idx][ ui_IDX ]...]

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

    idq_wt_pad_view = view(idq_wt_pad, 1:length( uh_state ) )

    #----------------------------------------------------

    global_pf_views = ( working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view )

    sd_pf_views = ( working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, gen_nodes_δ_ω_ed_dash_eq_dash_views )

    # ---------------------------------------------------
    # global_pf param
    # ---------------------------------------------------

    global_pf_param = ( pf_net_param, sd_pf_views, mismatch )

    branches_name  = collect(keys( netd.edges ))

    nodes_name     = collect(keys( netd.nodes ))    

    #----------------------------------------------------
    #----------------------------------------------------

    named_tup_pf_result = power_balance_powerflow( x0_vh_θh, mismatch, sd_pf_views, (nodes_name, branches_name) , pf_net_param ; maxiter=maxiter, ftol=ftol , xtol=xtol, with_δ_ed_eq = with_δ_ed_eq )

    #----------------------------------------------------

    bus_dict_init    = named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init

    #---------------------------------------------------
    #---------------------------------------------------

    state .= industrial_model_init_operationpoint(
        netd, bus_dict_init; pure = :pure, no_control_device = only_gen )

    return (; dict_sys_to_industry, pure_states_Idx_in_system, ur_ui_Idx_in_system, industrial_Idx, industrial_model_pure_states_Idx, industrial_model_ur_ui_Idx, pure_states_and_ur_ui_Idx_in_system, industrial_model_state_Idx, net_to_industrial_idx_conversion_dict, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, nodes_ω_ed_dash_eq_dash_Idxs_in_state, nodes_δ_ed_dash_eq_dash_Idxs_in_state, industrial_model_each_gen_nodes_pure_states_idx_in_state, industrial_model_each_gen_nodes_stab_states_idx_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view, gen_nodes_dyn_param_view, gen_nodes_sub_dyn_param_view, pf_net_param, ra_Xd_dash_Xq_dash_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, industrial_model_pure_states_view_in_state, state_view, pf_state, nodes_u_view, nodes_pf_U_view, x0_vh_θh, working_vh_θh_view, red_vh_θh_0_view, mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, branches_name, nodes_name, named_tup_pf_result, bus_dict_init, branch_dict_init, state, gens_idx, nodes_u_Idx )
    
    
end


"""

Test driver for:

`get_industrial_model_pf_param_views_and_init`

"""
function driver_get_industrial_model_pf_param_views_and_init( )

    dynamics_case = case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd  = NetworkData( dynamics_case()... )

    only_gen = true

    dict_sys_to_industry, pure_states_Idx_in_system, ur_ui_Idx_in_system, industrial_Idx, industrial_model_pure_states_Idx, industrial_model_ur_ui_Idx, pure_states_and_ur_ui_Idx_in_system, industrial_model_state_Idx, net_to_industrial_idx_conversion_dict, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, nodes_ω_ed_dash_eq_dash_Idxs_in_state, nodes_δ_ed_dash_eq_dash_Idxs_in_state, industrial_model_each_gen_nodes_pure_states_idx_in_state, industrial_model_each_gen_nodes_stab_states_idx_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view, gen_nodes_dyn_param_view, gen_nodes_sub_dyn_param_view, pf_net_param, ra_Xd_dash_Xq_dash_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, industrial_model_pure_states_view_in_state, state_view, pf_state, nodes_u_view, nodes_pf_U_view, x0_vh_θh, working_vh_θh_view, red_vh_θh_0_view, mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, branches_name, nodes_name, named_tup_pf_result, bus_dict_init, branch_dict_init, state, gens_idx, nodes_u_Idx = get_industrial_model_pf_param_views_and_init( netd; maxiter=40, ftol=1000*eps() , xtol=1000*eps(), init_pf = true, with_δ_ed_eq = false,  only_gen = only_gen )

    return nothing
    
end


# driver_get_industrial_model_pf_param_views_and_init( )

function get_industrial_model_pf_param_views_and_init_with_or_no_controller(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen =false  )

    # ---------------------------------------------------
    # ---------------------------------------------------

    nodes_cb_sw           = get_nodes_cb_sw(netd.nodes)

    gens_nodes            = get_gens_nodes( netd.nodes  )

    non_gens_nodes        = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #---------------------------------------------------
    #---------------------------------------------------- 

    dict_sys_to_industry, pure_states_Idx_in_system, ur_ui_Idx_in_system, industrial_Idx, industrial_model_pure_states_Idx, industrial_model_ur_ui_Idx, pure_states_and_ur_ui_Idx_in_system, industrial_model_state_Idx, net_to_industrial_idx_conversion_dict, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, nodes_ω_ed_dash_eq_dash_Idxs_in_state, nodes_δ_ed_dash_eq_dash_Idxs_in_state, industrial_model_each_gen_nodes_pure_states_idx_in_state, industrial_model_each_gen_nodes_stab_states_idx_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view, gen_nodes_dyn_param_view, gen_nodes_sub_dyn_param_view, pf_net_param, ra_Xd_dash_Xq_dash_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, industrial_model_pure_states_view_in_state, state_view, pf_state, nodes_u_view, nodes_pf_U_view, x0_vh_θh, working_vh_θh_view, red_vh_θh_0_view, mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, branches_name, nodes_name, named_tup_pf_result, bus_dict_init, branch_dict_init, state, gens_idx, nodes_u_Idx = get_industrial_model_pf_param_views_and_init(
        netd;
        maxiter=40,
        ftol=1000*eps(),
        xtol=1000*eps(),
        init_pf = true,
        with_δ_ed_eq = false,
        only_gen  =false)
    
    #---------------------------------------------------- 
    #---------------------------------------------------- 
    
    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!( gen_nodes_δ_ω_ed_dash_eq_dash_views, state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )


    #----------------------------------------------------

    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)

    #----------------------------------------------------
    #----------------------------------------------------
    #---------------------------------------------------

    gens_vh_θh = get_gens_vh_θh( nodes_pf_U_view, gens_idx )

    gens_vh_θh_view = @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    #---------------------------------------------------
    #---------------------------------------------------

    # gen_uh  = (named_tup_pf_result.Vbus)[ gen_idx ]

    idq_wt_pad_view[gens_idx] .=  [ industrial_model_get_dynamic_idq_θ_π_vhθh( vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )  for ( vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in zip( gens_vh_θh_view , gen_nodes_δ_ω_ed_dash_eq_dash_views, ra_Xd_dash_Xq_dash_view[gens_idx] ) ]

    #---------------------------------------------------
    #---------------------------------------------------

    gens_nodes_ωs_τm_vref_porder = get_gens_nodes_ωs_τm_vref_porder( gens_nodes_collection, bus_dict_init )

    gens_nodes_ωs_τm_vref_porder_view = view( gens_nodes_ωs_τm_vref_porder, 1:length(gens_nodes_ωs_τm_vref_porder) )

    #---------------------------------------------------

    gens_nodes_τm_vf = get_gens_τm_vf( gens_nodes_collection, bus_dict_init )

    vec_τm_vf_views = view( gens_nodes_τm_vf, 1:length(gens_nodes_τm_vf) )
    
    #---------------------------------------------------   
    #---------------------------------------------------

    gens_dynamic_id_iq_pg_vh_by_vhθh = get_gens_dynamic_id_iq_pg_vh_by_vhθh( gens_vh_θh_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, gen_nodes_ra_Xd_dash_Xq_dash_view )

    gens_dynamic_id_iq_pg_vh_by_vhθh_view = view( gens_dynamic_id_iq_pg_vh_by_vhθh, 1:length(gens_dynamic_id_iq_pg_vh_by_vhθh) )

    """
    gens_dynamic_id_iq_pg_vh_by_ur_ui = get_gens_dynamic_id_iq_pg_vh_by_ur_ui( gens_ur_ui_post_pf, gen_nodes_δ_ω_ed_dash_eq_dash_views, gen_nodes_ra_Xd_dash_Xq_dash_view )

    """
    #---------------------------------------------------  

    nodes_u_Idx_in_ranges = get_nodes_u_Idx_in_ranges( nodes_u_Idx )    
    #---------------------------------------------------

    non_gens_idx = get_load_trans_nodes_Idx( netd.nodes )

    #---------------------------------------------------

    industrial_model_nodes_voltage_Idx = (; gens_idx, non_gens_idx, nodes_u_Idx, nodes_u_Idx_in_ranges )

    industrial_model_misc_Idx = (; nodes_u_Idx, nodes_u_Idx_in_ranges, industrial_model_pure_states_Idx, industrial_model_each_gen_nodes_pure_states_idx_in_state, gens_nodes_collection )

    para_update_gen_Ax_aux = (; industrial_model_pure_states_view_in_state, industrial_model_each_gen_nodes_pure_states_idx_in_state, industrial_model_pure_states_Idx )

    industrial_model_para_aux_inputs = (; industrial_model_nodes_voltage_Idx, industrial_model_pure_states_Idx, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,  industrial_model_pure_states_view_in_state, gen_nodes_δ_ω_ed_dash_eq_dash_views, gens_vh_θh_view, nodes_pf_U_view )

    if only_gen == false
        
        industrial_model_pf_para = (; gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_τm_vref_porder_view )
        
    else
        industrial_model_pf_para = (; gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_τm_vref_porder_view, vec_τm_vf_views )
        
    end


    """ need by their views """

    industrial_model_ωs_τm_vref_vhθh_idq =  (; gens_dynamic_id_iq_pg_vh_by_vhθh, gens_nodes_ωs_τm_vref_porder, gens_nodes_τm_vf, gens_vh_θh, idq_wt_pad )


    industrial_model_dyn_pf_up_para = (; nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, nodes_ω_ed_dash_eq_dash_Idxs_in_state, nodes_δ_ed_dash_eq_dash_Idxs_in_state )

    #---------------------------------------------------

     industrial_model_idq_pf_cal = (;  idq_wt_pad_view, gens_idx )

    #---------------------------------------------------

return (; nodes_cb_sw, state, global_pf_param, named_tup_pf_result, industrial_model_misc_Idx, para_update_gen_Ax_aux, industrial_model_para_aux_inputs, industrial_model_pf_para, industrial_model_ωs_τm_vref_vhθh_idq, industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal )

    #---------------------------------------------------
    #---------------------------------------------------

end


"""

# 

Test driver for:

` get_industrial_model_pf_param_views_and_init_with_or_no_controller`

"""
function driver_get_industrial_model_pf_param_views_and_init_with_or_no_controller( )

    dynamics_case = case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd  = NetworkData( dynamics_case()... )

    only_gen = true

        nodes_cb_sw, state, global_pf_param, named_tup_pf_result, industrial_model_misc_Idx, para_update_gen_Ax_aux, industrial_model_para_aux_inputs, industrial_model_pf_para, industrial_model_ωs_τm_vref_vhθh_idq , industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal = get_industrial_model_pf_param_views_and_init_with_or_no_controller( netd; maxiter=40, ftol=1000*eps() , xtol=1000*eps(), init_pf = true, with_δ_ed_eq = false,  only_gen = only_gen )

               
    return nothing
    
end


"""

driver_get_industrial_model_pf_param_views_and_init_with_or_no_controller( )

"""

#-----------------------------------------------------
#-----------------------------------------------------


function get_im_model_pf_param_views_and_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen  =false)

    #-------------------------------
    
    # dynamics_case = case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    # netd     = NetworkData( dynamics_case()... )
    # only_gen = true
    # maxiter  = 40
    # ftol     = 1000*eps()
    # xtol     = 1000*eps()
    # init_pf  = true
    # with_δ_ed_eq = false
    
    # #-------------------------------
    

    dict_sys_to_im = get_net_to_im_indices_dict( netd  )

    pure_states_Idx_in_system, im_algebraic_vars_Idx_in_system, ur_ui_Idx_in_system, im_Idx, pure_states_Idx, im_algebraic_vars_Idx, nodes_ur_ui_Idx, pure_states_im_algebraic_vars_ur_ui_Idx_in_system, im_state_Idx, net_to_im_idx_conversion_dict = get_im_indices_and_conversion_dict( netd  )

    
    # ---------------------------------------------------
    # ---------------------------------------------------    

    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
        get_im_δ_ω_ed_dash_eq_dash_Idx(netd )

    nodes_ω_ed_dash_eq_dash_Idxs_in_state =
        get_im_ω_ed_dash_eq_dash_Idx(netd )

    nodes_δ_ed_dash_eq_dash_Idxs_in_state =
        get_im_δ_ed_dash_eq_dash_Idx(netd )

    #---------------------------------------------------    

    im_model_each_gen_nodes_pure_states_idx_in_state =
        get_gens_im_pure_states_indices_in_state(netd )

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

    pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx


    #----------------------------------------------------    
    #----------------------------------------------------

    state = zeros(length(
        generate_im_model_sym(
            ; nodes = netd.nodes )  ) )

    #----------------------------------------------------   

    gen_nodes_δ_ω_ed_dash_eq_dash_views =
        get_im_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
            state,
            nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    im_model_pure_states_view_in_state =
        get_im_pure_states_view_in_state(
            state,
            pure_states_Idx )

    #---------------------------------------------------- 
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

    idq_wt_pad_view = view(idq_wt_pad, 1:length( uh_state ) )

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

    # ------------------------------------------------
    # global_pf param
    # ------------------------------------------------

    global_pf_param = (
        pf_net_param,
        sd_pf_views,
        mismatch )

    branches_name  = collect(keys( netd.edges ))

    nodes_name     = collect(keys( netd.nodes ))    

    #-------------------------------------------------
    #-------------------------------------------------

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

    #-------------------------------------------------

    bus_dict_init = named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init

    #-------------------------------------------------
    #-------------------------------------------------

    state .= im_model_init_operationpoint(
        netd, bus_dict_init  )

    return (
        ;pure_states_Idx_in_system,
        ur_ui_Idx_in_system,
        im_Idx,
        pure_states_Idx,
        nodes_ur_ui_Idx,
        pure_states_im_algebraic_vars_ur_ui_Idx_in_system,
        im_state_Idx,
        net_to_im_idx_conversion_dict,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_δ_ed_dash_eq_dash_Idxs_in_state,
        im_model_each_gen_nodes_pure_states_idx_in_state,
        gen_nodes_ra_Xd_dash_Xq_dash_view,
        gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
        gen_nodes_dyn_param_view,
        gen_nodes_sub_dyn_param_view,
        pf_net_param,
        ra_Xd_dash_Xq_dash_view,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        im_model_pure_states_view_in_state,
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


"""

Test driver for:

`get_get_im_model_pf_param_views_and_init`

"""
function driver_get_im_model_pf_param_views_and_init( )

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd  = NetworkData( dynamics_case()... )

    only_gen = true

    pure_states_Idx_in_system, ur_ui_Idx_in_system, industrial_Idx, industrial_model_pure_states_Idx, industrial_model_ur_ui_Idx, pure_states_and_ur_ui_Idx_in_system, industrial_model_state_Idx, net_to_industrial_idx_conversion_dict, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, nodes_ω_ed_dash_eq_dash_Idxs_in_state, nodes_δ_ed_dash_eq_dash_Idxs_in_state, industrial_model_each_gen_nodes_pure_states_idx_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view, gen_nodes_dyn_param_view, gen_nodes_sub_dyn_param_view, pf_net_param, ra_Xd_dash_Xq_dash_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, industrial_model_pure_states_view_in_state, state_view, pf_state, nodes_u_view, nodes_pf_U_view, x0_vh_θh, working_vh_θh_view, red_vh_θh_0_view, mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, branches_name, nodes_name, named_tup_pf_result, bus_dict_init, branch_dict_init, state, gens_idx, nodes_u_Idx = get_im_model_pf_param_views_and_init(
        netd;
        maxiter=40,
        ftol=1000*eps() ,
        xtol=1000*eps(),
        init_pf = true,
        with_δ_ed_eq = false,
        only_gen = only_gen )

    return nothing
    
end


"""

driver_get_im_model_pf_param_views_and_init( )

"""


function get_im_sys_pf_param_views_and_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false ,
    only_gen =false  )

    # -----------------------------------------------

    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)

    gens_nodes = get_gens_nodes( netd.nodes  )

    non_gens_nodes = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #------------------------------------------------

    pure_states_Idx_in_system,
    ur_ui_Idx_in_system,
    industrial_Idx,
    im_pure_states_Idx, #
    im_ur_ui_Idx,
    pure_states_and_ur_ui_Idx_in_system,
    im_state_Idx,
    net_to_industrial_idx_conversion_dict,
    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
    nodes_ω_ed_dash_eq_dash_Idxs_in_state,
    nodes_δ_ed_dash_eq_dash_Idxs_in_state,
    im_each_gen_nodes_pure_states_idx_in_state, #
    gen_nodes_ra_Xd_dash_Xq_dash_view,
    gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
    gen_nodes_dyn_param_view,
    gen_nodes_sub_dyn_param_view,
    pf_net_param,
    ra_Xd_dash_Xq_dash_view,
    gen_nodes_δ_ω_ed_dash_eq_dash_views,
    im_pure_states_view_in_state, #
    state_view, pf_state, nodes_u_view, nodes_pf_U_view, x0_vh_θh, working_vh_θh_view, red_vh_θh_0_view, mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, branches_name, nodes_name, named_tup_pf_result, bus_dict_init, branch_dict_init, state, gens_idx, nodes_u_Idx = get_im_model_pf_param_views_and_init(
        netd;
        maxiter=40,
        ftol=1000*eps() ,
        xtol=1000*eps(),
        init_pf = true,
        with_δ_ed_eq = false,
        only_gen = only_gen )
    
    #---------------------------------------------- 
    
    update_im_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )


    #-----------------------------------------------

    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)

    #-----------------------------------------------

    gens_vh_θh = get_gens_vh_θh( nodes_pf_U_view, gens_idx )

    gens_vh_θh_view = @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    #-----------------------------------------------

    # gen_uh  = (named_tup_pf_result.Vbus)[ gen_idx ]

    idq_wt_pad_view[gens_idx] .=  [
        industrial_model_get_dynamic_idq_θ_π_vhθh(
            vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
        for ( vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in zip(
            gens_vh_θh_view ,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            ra_Xd_dash_Xq_dash_view[gens_idx] ) ]

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
        im_pure_states_Idx, 
        im_each_gen_nodes_pure_states_idx_in_state, 
        gens_nodes_collection )

    
    para_update_gen_Ax_aux = (
        ; im_pure_states_view_in_state, 
        im_each_gen_nodes_pure_states_idx_in_state, 
        im_pure_states_Idx ) 

    
    im_para_aux_inputs = (
        ; im_nodes_voltage_Idx,
        im_pure_states_Idx, 
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
        gen_nodes_ra_Xd_dash_Xq_dash_view,
        gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,
        im_pure_states_view_in_state,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        gens_vh_θh_view,
        nodes_pf_U_view )

    
    im_pf_para = (
        ; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
        gens_nodes_ωs_ωref0_vref0_porder0_view  )


    """ need by their views """

    im_ωs_τm_vref_vhθh_idq =  (
        ; gens_dynamic_id_iq_pg_vh_by_vhθh,
        gens_nodes_ωs_ωref0_vref0_porder0,
        gens_nodes_τm_vf,
        gens_vh_θh, idq_wt_pad )

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
        ; nodes_cb_sw,
        state,
        global_pf_param,
        named_tup_pf_result,
        im_misc_Idx, 
        para_update_gen_Ax_aux, 
        im_para_aux_inputs, 
        im_pf_para, 
        im_ωs_τm_vref_vhθh_idq, 
        im_dyn_pf_up_para, 
        im_idq_pf_cal ) 
    
    #---------------------------------------------------

end


"""

Test driver for:

` get_get_im_sys_pf_param_views_and_init`

"""
function driver_get_im_sys_pf_param_views_and_init( )

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd  = NetworkData( dynamics_case()... )

    only_gen = true

    nodes_cb_sw, state, global_pf_param, named_tup_pf_result, im_misc_Idx, para_update_gen_Ax_aux, im_para_aux_inputs, im_pf_para, im_ωs_τm_vref_vhθh_idq, im_dyn_pf_up_para, im_idq_pf_cal = get_im_sys_pf_param_views_and_init( netd; maxiter=40, ftol=1000*eps() , xtol=1000*eps(), init_pf = true, with_δ_ed_eq = false,  only_gen = only_gen )

               
    return nothing
    
end


"""

driver_get_im_sys_pf_param_views_and_init( )

"""


function diagnosis_simulate_a_dynamic_im_model(
    ; dynamics_case = nothing,
    only_gen = false,
    sim_timespan  = nothing,
    algr = algr,
    dyn_global_pf_options = nothing ,
    sta_global_pf_options = nothing,
    with_cb  = false)

    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    # dynamics_case =
    #     case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer

    # sim_timespan          = 1.0
    
    # dyn_global_pf_options = (
    #     ; maxiter         = 40,
    #     ftol              = 1000*eps(),
    #     xtol              = 1000*eps(),
    #     with_δ_ed_eq      = true )

    # sta_global_pf_options = (
    #     ; maxiter         = 40,
    #     ftol              = 1000*eps(),
    #     xtol              = 1000*eps(),
    #     init_pf           = true,
    #     with_δ_ed_eq      = false )

    # only_gen = true

    # with_cb  = false
        
    #-----------------------------------------------------

    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------    

    non_gens_nodes = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #-----------------------------------------------------
    
    network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names,Trans_bus_names = make_case_buses_names(  netd.nodes )
    net_class_names = (;network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names, Trans_bus_names )
    
    #-----------------------------------------------------
    #-----------------------------------------------------


   net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_im_algebraic_vars_labels, gens_nodes_im_vars_labels = generate_im_model_labels(; nodes =  netd.nodes )

    im_net_states_and_var_labels = (; net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_im_algebraic_vars_labels, gens_nodes_im_vars_labels  )

    #-----------------------------------------------------
                                         
    im_model_sym, im_model_mass_matrix = generate_im_model_sym_and_mass_matrix(; nodes = netd.nodes )

    #-----------------------------------------------------

    para_net_names_labels_syms = (; net_class_names, im_net_states_and_var_labels, im_model_sym )

    #-----------------------------------------------------
       
    nodes_cb_sw, state, global_pf_param, named_tup_pf_result, im_misc_Idx, para_update_gen_Ax_aux, im_para_aux_inputs, im_pf_para, im_ωs_τm_vref_vhθh_idq, im_dyn_pf_up_para, im_idq_pf_cal =
        get_im_sys_pf_param_views_and_init(
            netd; sta_global_pf_options...,
            only_gen = false )

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

    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix =
        Ax_Bx_Cx_matrix

    nodes_state_Idx, Bx_idxs, Cx_idxs = idxs
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views, vec_Bx_views, vec_Cx_views),
        gens_nodes_collection )

    #-----------------------------------------------------
    #-----------------------------------------------------

    ode_fun_para = (
        ; Ax_matrix,
          Bx_matrix,
          Cx_matrix,
          im_pf_para  )

    para_model_fun = (
        ; vec_Ax_views,
        ode_fun_para )
    
    #-----------------------------------------------------

    sim_state_x0  = state

    sim_timespan  = sim_timespan

    #-----------------------------------------------------

    dyn_global_pf_options = dyn_global_pf_options

    counter_array  = [1]

    stateDiffCache = similar( sim_state_x0 )

    sim_fun_para = (
        ;netd,
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
        im_idq_pf_cal,                
        only_gen )

    return Ax_Bx_Cx_views, sim_fun_para

    # return nothing
      
end

"""
Test:

`diagnosis_simulate_a_dynamic_im_model`


"""
function driver_diagnosis_simulate_a_dynamic_im_model()

    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)
    algr                  = Rodas4()
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

    Ax_Bx_Cx,  sim_fun_para  =
        diagnosis_simulate_a_dynamic_im_model(
        ;
        dynamics_case         = dynamics_case,
        only_gen              = only_gen,
        sim_timespan          = sim_timespan,
        algr                  = algr,
        dyn_global_pf_options = dyn_global_pf_options,
        sta_global_pf_options = sta_global_pf_options,
        with_cb  = with_cb )

    #--------------------------------------------
    #--------------------------------------------
    
    nd, nodes_cb_sw, global_pf_param, counter_array, state, stateDiffCache, dyn_global_pf_options, para_update_gen_Ax_aux, para_model_fun, im_misc_Idx, im_para_aux_inputs, im_dyn_pf_up_para, im_idq_pf_cal = sim_fun_para
    
    #--------------------------------------------

    nodes_u_Idx, nodes_u_Idx_in_ranges, im_pure_states_Idx, im_each_gen_nodes_pure_states_idx_in_state, gens_nodes_collection = im_misc_Idx

    #--------------------------------------------
    
    vec_Ax_views, ode_fun_para = para_model_fun

    Ax_matrix, Bx_matrix, Cx_matrix, im_pf_para = ode_fun_para

    # gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_ωref0_vref0_porder0_view
    
    gens_dynamic_id_iq_pg_vh_by_vhθh, gens_nodes_ωs_ωref0_vref0_porder0_view = im_pf_para
        
    #--------------------------------------------

    pf_net_param, sd_pf_views, mismatch = global_pf_param

    working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view = sd_pf_views
    #--------------------------------------------    

    im_pure_states_view_in_state, _, _ = para_update_gen_Ax_aux
    
    #--------------------------------------------    
    #--------------------------------------------    
    # Powerflow
    #--------------------------------------------
    #--------------------------------------------
        
    counter = counter_array[1]

    x  = state
    dx = similar(x)
    
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

    update_dynamic_id_iq_pg_vh!( gens_dynamic_id_iq_pg_vh_by_vhθh, stateDiffCache, im_para_aux_inputs )
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    im_pure_states_view_in_state .= get_industrial_gen_nodes_pure_states_from_state(stateDiffCache, im_pure_states_Idx )
    
    update_gens_nodes_im_Ax_system_matrices!(
        vec_Ax_views,
        im_pure_states_view_in_state,
        im_each_gen_nodes_pure_states_idx_in_state,
        gens_nodes_collection )


    #--------------------------------------------
    # ode gens pure state
    #--------------------------------------------

    tt = 0.0
    
    ode_im_model_func!(
    dx,
    x,
    (im_pure_states_Idx, ode_fun_para),
    tt)

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


    #--------------------------------------------
    #--------------------------------------------

    
    x0 = sim_fun_para.state
    x0_ps_idx = sim_fun_para.im_misc_Idx.im_pure_states_Idx
    x0_gen_ps_idx = sim_fun_para.im_misc_Idx.im_each_gen_nodes_pure_states_idx_in_state

    t_vec_Ax, t_vec_Bx, t_vec_Cx = Ax_Bx_Cx

    t_vec_Ax = sim_fun_para.para_model_fun[1]

    # .ode_fun_para.Ax_matrix
    t_Ax_m = sim_fun_para.para_model_fun[2][1]

    # .ode_fun_para.Bx_matrix
    t_Bx_m = sim_fun_para.para_model_fun[2][2]

    # .ode_fun_para.Cx_matrix
    t_Cx_m = sim_fun_para.para_model_fun[2][3]

    # .ode_fun_para.im_pf_para
    im_pf_para = sim_fun_para.para_model_fun[2][4]

    id_iq_pg_vh = im_pf_para.gens_dynamic_id_iq_pg_vh_by_vhθh_view

    ωs_ω_ref_vref_porder = im_pf_para.gens_nodes_ωs_ωref0_vref0_porder0_view

    t_dx = t_Ax_m * x0[ x0_ps_idx ] + t_Bx_m * [id_iq_pg_vh...;] + t_Cx_m * [ ωs_ω_ref_vref_porder...; ]

    
    return nothing
    
    # return sim_fun_para
    
        
end

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
        view( dx, idx )
        for idx in  nodes_u_Idx_in_ranges[ non_gens_idx ] ]
    
    x_non_gen_ur_ui_views   = [
        view( x,  idx )
        for idx in nodes_u_Idx_in_ranges[ non_gens_idx ] ]
    
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

    gens_idx, non_gens_idx, nodes_u_Idx, nodes_u_Idx_in_ranges = im_nodes_voltage_Idx         
    
    dx_gen_ur_ui_views =
        [ view( dx, idx )
          for idx in
              nodes_u_Idx_in_ranges[gens_idx] ]
    
    x_gen_ur_ui_views =
        [ view( x,  idx )
          for idx in
              nodes_u_Idx_in_ranges[gens_idx] ]
    
    for (du, u, gens_vh_θh, gen_nodes_δ_ω_ed_dash_eq_dash, gen_nodes_ra_Xd_dash_Xq_dash) in zip(
        dx_gen_ur_ui_views, x_gen_ur_ui_views,
        gens_vh_θh_view, gen_nodes_δ_ω_ed_dash_eq_dash_views,
        gen_nodes_ra_Xd_dash_Xq_dash_view )

        a_gen_im_voltage_terminal_func!(
            du, u, (gens_vh_θh,
                    gen_nodes_δ_ω_ed_dash_eq_dash,
                    gen_nodes_ra_Xd_dash_Xq_dash), t) 

    end
    
end


function a_gen_im_voltage_terminal_func!(
    dx, x, (gen_vh_θh, gen_node_δ_ω_ed_dash_eq_dash,
            gen_node_ra_Xd_dash_Xq_dash), t) 
            
    δ, _, ed_dash, eq_dash = gen_node_δ_ω_ed_dash_eq_dash

    id_iq = get_dynamic_idq_vhθh(
        gen_vh_θh...,
        gen_node_δ_ω_ed_dash_eq_dash...,
        gen_node_ra_Xd_dash_Xq_dash... )

    zdq  = Z_dq( gen_node_ra_Xd_dash_Xq_dash... )

    dx .= [ sin(δ) cos(δ); -cos(δ) sin(δ)] * ( [ed_dash, eq_dash] - zdq * id_iq ) - x 
        
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
    dx, x, (nodes_pf_U_view, nodes_u_Idx_in_ranges), t)
        
    dx_non_gen_ur_ui_views = [
        view(dx, idx) for idx in nodes_u_Idx_in_ranges ]
    
    x_non_gen_ur_ui_views  = [
        view(x, idx) for idx in nodes_u_Idx_in_ranges ]
   
    for  ( dx_pf, x_pf, u_pf ) in zip(
        dx_non_gen_ur_ui_views,
        x_non_gen_ur_ui_views, nodes_pf_U_view )

        algebraic_industrial_model_func!(
            dx_pf, x_pf, u_pf, t)
        
    end
    
end



#---------------------------------------------------
#---------------------------------------------------


function non_gens_industrial_model_voltage_terminal_func!(
    dx, x, industrial_model_para_aux_inputs, t)

    industrial_model_nodes_voltage_Idx, _, _, _, _,_, _, _, nodes_pf_U_view = industrial_model_para_aux_inputs

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

    industrial_model_nodes_voltage_Idx, _, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, _,_, gen_nodes_δ_ω_ed_dash_eq_dash_views, gens_vh_θh_view, nodes_pf_U_view = industrial_model_para_aux_inputs

     gens_idx, non_gens_idx, nodes_u_Idx, nodes_u_Idx_in_ranges  = industrial_model_nodes_voltage_Idx         
    
    dx_gen_ur_ui_views = [ view( dx, idx ) for idx in  nodes_u_Idx_in_ranges[gens_idx] ]
    
    x_gen_ur_ui_views   = [ view( x,  idx ) for idx in nodes_u_Idx_in_ranges[gens_idx] ]
    
    for (du, u, gens_vh_θh, gen_nodes_δ_ω_ed_dash_eq_dash, gen_nodes_ra_Xd_dash_Xq_dash) in zip( dx_gen_ur_ui_views, x_gen_ur_ui_views, gens_vh_θh_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, gen_nodes_ra_Xd_dash_Xq_dash_view )

        a_gen_industrial_model_voltage_terminal_func!(
            du, u, (gens_vh_θh, gen_nodes_δ_ω_ed_dash_eq_dash, gen_nodes_ra_Xd_dash_Xq_dash), t) 

    end
    
end


function a_gen_industrial_model_voltage_terminal_func!(
    dx, x, (gen_vh_θh, gen_node_δ_ω_ed_dash_eq_dash,
            gen_node_ra_Xd_dash_Xq_dash), t) 
            
    δ, _, ed_dash, eq_dash = gen_node_δ_ω_ed_dash_eq_dash

    id_iq = get_dynamic_idq_vhθh(
        gen_vh_θh...,
        gen_node_δ_ω_ed_dash_eq_dash...,
        gen_node_ra_Xd_dash_Xq_dash... )

    zdq  = Z_dq( gen_node_ra_Xd_dash_Xq_dash... )

    dx .= [ sin(δ) cos(δ); -cos(δ) sin(δ)] * ( [ed_dash, eq_dash] - zdq * id_iq ) - x 
        
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
    dx, x, (nodes_pf_U_view, nodes_u_Idx_in_ranges)  , t)
        
    dx_non_gen_ur_ui_views = [
        view(dx, idx) for idx in nodes_u_Idx_in_ranges ]
    
    x_non_gen_ur_ui_views  = [
        view(x, idx) for idx in nodes_u_Idx_in_ranges ]
   
    for  ( dx_pf, x_pf, u_pf ) in zip(
        dx_non_gen_ur_ui_views,
        x_non_gen_ur_ui_views, nodes_pf_U_view )

        algebraic_industrial_model_func!(
            dx_pf, x_pf, u_pf, t)
        
    end
    
end



function ode_industrial_model_func!(dx, x, (industrial_model_pure_states_Idx, ode_fun_para), t)

   
    dx_gen  = @view dx[ industrial_model_pure_states_Idx ]
    
    x_gen   = @view x[ industrial_model_pure_states_Idx ]
    

    Ax_matrix, Bx_matrix, Cx_matrix, industrial_model_pf_para = ode_fun_para

    gens_dynamic_id_iq_pg_vh, gens_nodes_ωs_τm_vref_porder_view = industrial_model_pf_para
        
    dx_gen .= Ax_matrix * x_gen + Bx_matrix * [ gens_dynamic_id_iq_pg_vh...; ] + Cx_matrix * [ gens_nodes_ωs_τm_vref_porder_view...;]

    return nothing
    
end



function ode_industrial_model_func_no_controllers!(
    dx, x,
    (ode_fun_para,
     industrial_model_pure_states_Idx, only_gen),
    t)


    dx_gen  = @view dx[
        industrial_model_pure_states_Idx ]

    x_gen   = @view x[
        industrial_model_pure_states_Idx ]
    
    if only_gen == false
        
        Ax_matrix, Bx_matrix, Cx_matrix, industrial_model_pf_para = ode_fun_para

        id_iq_pg_vh, ωs_τm_vref_porder = industrial_model_pf_para

        dx_gen .= Ax_matrix * x_gen + Bx_matrix * [ id_iq_pg_vh...; ] + Cx_matrix * [ ωs_τm_vref_porder...;]
        
    else
        Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix,  industrial_model_pf_para = ode_fun_para

        id_iq_pg_vh, ωs_τm_vref_porder, τm_vf = industrial_model_pf_para

        dx_gen .= Ax_matrix * x_gen + Bx_matrix * [ id_iq_pg_vh...; ] + Cx_matrix * [ ωs_τm_vref_porder...; ] + Ax_τm_vf_matrix * [ τm_vf...; ]

    end
    

    return nothing
    
end



# im_vars_view_in_state,
# im_vars_Idx_in_state,
# each_gens_im_vars_Idx_in_state,

function ode_im_model_func!(
    dx,
    x,
    ( im_vars_Idx_in_state, ode_fun_para),
    t)

   
    dx_gen  = @view dx[ im_vars_Idx_in_state ]
    
    x_gen   = @view x[ im_vars_Idx_in_state ]
    

    Ax_matrix, Bx_matrix, Cx_matrix, im_pf_para = ode_fun_para

    id_iq_pg_vh, ωs_ω_ref_vref_porder = im_pf_para
        
    dx_gen .=
        Ax_matrix * x_gen +
        Bx_matrix * [ id_iq_pg_vh...; ] +
        Cx_matrix * [ ωs_ω_ref_vref_porder...;]

    return nothing
    
end


#

function update_im_dynamic_id_iq_pg_vh!(
    gens_dynamic_id_iq_pg_vh_by_vhθh,
    stateDiffCache,
    im_para_aux_inputs )

    im_nodes_voltage_Idx, im_vars_Idx_in_state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view, im_vars_view_in_state, gen_nodes_δ_ω_ed_dash_eq_dash_views, gens_vh_θh_view, nodes_pf_U_view = im_para_aux_inputs

     # im_nodes_voltage_Idx, industrial_model_pure_states_Idx, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,  industrial_model_pure_states_view_in_state, gen_nodes_δ_ω_ed_dash_eq_dash_views, gens_vh_θh_view, nodes_pf_U_view = industrial_model_para_aux_inputs

    gens_idx, non_gens_idx, nodes_u_Idx, nodes_u_Idx_in_ranges  =
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

     industrial_model_nodes_voltage_Idx, industrial_model_pure_states_Idx, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,  industrial_model_pure_states_view_in_state, gen_nodes_δ_ω_ed_dash_eq_dash_views, gens_vh_θh_view, nodes_pf_U_view = industrial_model_para_aux_inputs

    gens_idx, non_gens_idx, nodes_u_Idx, nodes_u_Idx_in_ranges  =
        industrial_model_nodes_voltage_Idx
        
    industrial_model_pure_states_view_in_state[:] .=
        get_industrial_gen_nodes_pure_states_from_state(
            stateDiffCache, industrial_model_pure_states_Idx )

    update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        industrial_model_pure_states_view_in_state,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    # gen_nodes_δ_ω_ed_dash_eq_dash_views .= industrial_model_pure_states_view_in_state[ nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state ] 
    
    gens_vh_θh_view[:] .= get_gens_vh_θh( nodes_pf_U_view, gens_idx )
    
    gens_dynamic_id_iq_pg_vh_by_vhθh .=
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view, gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )
    
    return nothing
    
end



function update_dynamic_τm_vf!(
    vec_τm_vf, stateDiffCache,
    industrial_model_para_aux_inputs )

     industrial_model_nodes_voltage_Idx, industrial_model_pure_states_Idx, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view,  industrial_model_pure_states_view_in_state, gen_nodes_δ_ω_ed_dash_eq_dash_views, gens_vh_θh_view, nodes_pf_U_view = industrial_model_para_aux_inputs

    gens_idx, non_gens_idx, nodes_u_Idx, nodes_u_Idx_in_ranges  =
        industrial_model_nodes_voltage_Idx
    
    # gen_nodes_state_views .= get_gen_nodes_state_views(stateDiffCache, gens_nodes_collection )
        
    industrial_model_pure_states_view_in_state[:] .= get_industrial_gen_nodes_pure_states_from_state( stateDiffCache, industrial_model_pure_states_Idx )

    # update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!( gen_nodes_δ_ω_ed_dash_eq_dash_views, industrial_model_pure_states_view_in_state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    gen_nodes_δ_ω_ed_dash_eq_dash_views[:] .= industrial_model_pure_states_view_in_state[ nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state ] 
    
    gens_vh_θh_view[:] .= get_gens_vh_θh( nodes_pf_U_view, gens_idx )
    # nodes_u_Idx_in_ranges
    
    vec_τm_vf .= get_gens_dynamic_τm_vf( gens_vh_θh_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view  )
    
    return nothing
    
end




function dynamics_im_model!(dx, x, para, t)

    
    nd, nodes_cb_sw, global_pf_param, counter_array, state, stateDiffCache, dyn_global_pf_options, para_update_gen_Ax_aux, para_model_fun, im_misc_Idx, im_para_aux_inputs, im_dyn_pf_up_para, im_idq_pf_cal = para
    
    #--------------------------------------------

    nodes_u_Idx, nodes_u_Idx_in_ranges, im_vars_Idx_in_state, each_gens_im_vars_Idx_in_state, gens_nodes_collection = im_misc_Idx

    #--------------------------------------------
    
    vec_Ax_views, ode_fun_para = para_model_fun

    Ax_matrix, Bx_matrix, Cx_matrix, im_pf_para = ode_fun_para

    # gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_ωref0_vref0_porder0_view
    
    gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_ωref0_vref0_porder0_view = im_pf_para
        
    #--------------------------------------------

    pf_net_param, sd_pf_views, mismatch = global_pf_param

    working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view = sd_pf_views
    
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
    
    im_vars_view_in_state .= get_gen_nodes_im_vars_from_state(
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

    ode_im_model_func!(
        dx,
        x,
        ( im_vars_Idx_in_state, ode_fun_para),
        t)
    
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

    
    nd, nodes_cb_sw, global_pf_param, counter_array, state, stateDiffCache, dyn_global_pf_options, para_update_gen_Ax_aux, para_model_fun, industrial_model_misc_Idx, industrial_model_para_aux_inputs, industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal = para
    
    #--------------------------------------------

    nodes_u_Idx, nodes_u_Idx_in_ranges, industrial_model_pure_states_Idx, industrial_model_each_gen_nodes_pure_states_idx_in_state, gens_nodes_collection = industrial_model_misc_Idx

    #--------------------------------------------
    
    vec_Ax_views, ode_fun_para = para_model_fun

    Ax_matrix, Bx_matrix, Cx_matrix, industrial_model_pf_para = ode_fun_para

    # gens_dynamic_id_iq_pg_vh_by_vhθh, gens_nodes_ωs_τm_v_ref_view = industrial_model_pf_para

    gens_dynamic_id_iq_pg_vh_by_vhθh, gens_nodes_ωs_τm_vref_porder_view = industrial_model_pf_para
        
    #--------------------------------------------

    pf_net_param, sd_pf_views, mismatch = global_pf_param

    working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view = sd_pf_views
    #--------------------------------------------    

    industrial_model_pure_states_view_in_state, _, _ = para_update_gen_Ax_aux
    
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
    industrial_model_dyn_pf_up_para,
    industrial_model_idq_pf_cal;
    dyn_global_pf_options... )

    #--------------------------------------------
    # update  W
    #--------------------------------------------

    update_dynamic_id_iq_pg_vh!( gens_dynamic_id_iq_pg_vh_by_vhθh, stateDiffCache, industrial_model_para_aux_inputs )
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    industrial_model_pure_states_view_in_state .= get_industrial_gen_nodes_pure_states_from_state(stateDiffCache, industrial_model_pure_states_Idx )
    
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
        (industrial_model_pure_states_Idx, ode_fun_para),
        t)
    
    #--------------------------------------------
    # Another method for all nodes votages
    #--------------------------------------------

    # industrial_model_all_nodes_voltage_terminal_func!(dx, x, (nodes_pf_U_view, nodes_u_Idx_in_ranges)  , t)
            
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
    
    #--------------------------------------------
    #--------------------------------------------    
    
    return nothing
    
end

#-------------------------------------------------------


function dynamics_industrial_model_with_or_no_controller!(
    dx,
    x,
    para,
    t)
    
    nd, nodes_cb_sw, global_pf_param, counter_array, state, stateDiffCache, dyn_global_pf_options, para_update_gen_Ax_aux, para_model_fun, industrial_model_misc_Idx, industrial_model_para_aux_inputs, industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal, only_gen = para
    
    nodes_u_Idx, nodes_u_Idx_in_ranges, industrial_model_pure_states_Idx, industrial_model_each_gen_nodes_pure_states_idx_in_state, gens_nodes_collection = industrial_model_misc_Idx

    pf_net_param, sd_pf_views, mismatch = global_pf_param

    working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view = sd_pf_views
        
    industrial_model_pure_states_view_in_state, _, _ = para_update_gen_Ax_aux

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
    nd,
    stateDiffCache,
    global_pf_param,
    industrial_model_dyn_pf_up_para,
    industrial_model_idq_pf_cal;
    dyn_global_pf_options... )

    #--------------------------------------------
    # Filtering of pure states from all states 
    #--------------------------------------------    
    
    industrial_model_pure_states_view_in_state .= get_industrial_gen_nodes_pure_states_from_state(stateDiffCache, industrial_model_pure_states_Idx )
    
    if only_gen == false
            
        #--------------------------------------------

        vec_Ax_views, ode_fun_para = para_model_fun

        _, _, _, industrial_model_pf_para = ode_fun_para

        #--------------------------------------------
        #--------------------------------------------

        gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_τm_vref_porder_view = industrial_model_pf_para

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

        #--------------------------------------------
        #--------------------------------------------
        
    else
 
        #--------------------------------------------

        vec_Ax_views, ode_fun_para = para_model_fun

        Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix, industrial_model_pf_para = ode_fun_para

        #--------------------------------------------

        gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_τm_vref_porder_view, vec_τm_vf_views = industrial_model_pf_para

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
            gens_nodes_collection; only_gen = only_gen  )

        #--------------------------------------------
        # gens pure state
        #--------------------------------------------
        
        # dx_gen  = @view dx[ industrial_model_pure_states_Idx ]

        # x_gen   = @view x[ industrial_model_pure_states_Idx ]

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

        #--------------------------------------------
        #--------------------------------------------
        
    end

    
    #--------------------------------------------    
    
    return nothing
    
end


#-------------------------------------------------------
########################################################
#-------------------------------------------------------

function simulate_a_dynamics_industrial_model(
    ; dynamics_case = dynamics_case,
    only_gen = false,
    sim_timespan  = sim_timespan,
    algr          = algr,
    dyn_global_pf_options = dyn_global_pf_options,
    sta_global_pf_options = sta_global_pf_options )

    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    # dynamics_case = case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    #-----------------------------------------------------
    #----------------------------------------------------- 
    
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    # Dyn_Nodes, Dyn_Branches = dynamics_case()

    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------
    #-----------------------------------------------------

    non_gens_nodes        = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #-----------------------------------------------------

    # network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names,Trans_bus_names = make_case_buses_names(; case_fun = dynamics_case )

    # net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_stab_states_label, gens_nodes_algebraic_and_states_labels = generate_industrial_model_labels( network_bus_names, gens_nodes_collection )

    # industrial_model_sym, industrial_model_mass_matrix = generate_industrial_model_sym_and_mass_matrix( gens_nodes_pure_states_labels, net_bus_volts_labels)
    
    #-----------------------------------------------------
    
    network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names,Trans_bus_names = make_case_buses_names(  netd.nodes )
    net_class_names = (;network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names, Trans_bus_names )
    
    #-----------------------------------------------------

    net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_stab_states_label, gens_nodes_algebraic_and_states_labels = generate_industrial_model_labels(; nodes =  netd.nodes )

    net_states_and_var_labels = (; net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_stab_states_label, gens_nodes_algebraic_and_states_labels )
    
    #-----------------------------------------------------

    industrial_model_sym, industrial_model_mass_matrix = generate_industrial_model_sym_and_mass_matrix(; nodes = netd.nodes )

    #-----------------------------------------------------

    para_net_names_labels_syms = (; net_class_names, net_states_and_var_labels, industrial_model_sym )
    
    #-----------------------------------------------------
    
   industrial_model_pf_sys_param_sys_views_sys_industrial_model_init =  get_industrial_model_pf_sys_param_sys_views_sys_industrial_model_init( netd; sta_global_pf_options... )

    #-----------------------------------------------------
    #-----------------------------------------------------

     nodes_cb_sw, state, global_pf_param, named_tup_pf_result, industrial_model_misc_Idx, para_update_gen_Ax_aux, industrial_model_para_aux_inputs, industrial_model_pf_para, industrial_model_ωs_τm_vref_vhθh_idq, industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal, gens_nodes_τm_vf  = industrial_model_pf_sys_param_sys_views_sys_industrial_model_init

    
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
        Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs = create_gens_nodes_aggregate_system_matrices_Idx_and_views(  gens_nodes_collection; only_gen = only_gen )

        vec_Ax_views, vec_Bx_views, vec_Cx_views  = Ax_Bx_Cx_views
        Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

        nodes_state_Idx, Bx_idxs, Cx_idxs = idxs

        init_plants_system_matrices_views!((vec_Ax_views, vec_Bx_views, vec_Cx_views), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = nothing )
        
    else
        Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs = create_gens_nodes_aggregate_system_matrices_Idx_and_views(  gens_nodes_collection; only_gen = only_gen )
        
        vec_Ax_views, vec_Bx_views, vec_Cx_views, vec_Ax_τm_vf_views = Ax_Bx_Cx_views
        Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix = Ax_Bx_Cx_matrix
        nodes_state_Idx, Bx_idxs, Cx_idxs, τm_vf_idxs = idxs

        init_plants_system_matrices_views!((vec_Ax_views, vec_Bx_views, vec_Cx_views), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = vec_Ax_τm_vf_views )
        
    end


    """    
    if only_gen == false
        system_views, system_matrices, system_idxs =  make_system_matrices_with_link_matrices( gens_nodes_collection ;only_gen = only_gen )
        vec_Ax_views, Ax_matrix, Bx_matrix, Cx_matrix =  make_industrial_model_system_matrices( system_views, system_matrices, system_idxs, state, gens_nodes_collection ; only_gen = only_gen  )
    else
        system_views, system_matrices, system_idxs =  make_system_matrices_with_link_matrices( gens_nodes_collection ;only_gen = only_gen )
        vec_Ax_views, Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix, vec_τm_vf_views = make_industrial_model_system_matrices( system_views, system_matrices, system_idxs, state, gens_nodes_collection ; only_gen = only_gen  )
    end


    system_views, system_matrices, system_idxs =  make_system_matrices_with_link_matrices( gens_nodes_collection ;only_gen = false )

    update_gens_nodes_linked_system_matrices!( system_views, system_matrices, system_idxs, state, gens_nodes_collection ; only_gen = false )
    
    Ax_Bx_Cx_views, stab_Ax_Bx_Cx_views, stab_link_Ax_Bx_Cx_views = system_views

    Ax_Bx_Cx_matrix,  stab_Ax_Bx_Cx_matrix = system_matrices

    idxs, stab_idxs,  filter_stab_idxs = system_idxs

    #---------------------------------------------------

    vec_Ax_views, vec_Bx_views, vec_Cx_views = Ax_Bx_Cx_views

    vec_stab_Ax_views, vec_stab_Bx_views, vec_stab_Cx_views = stab_Ax_Bx_Cx_views

    vec_stab_link_Ax_views, vec_stab_link_Bx_views, vec_stab_link_Cx_views = stab_link_Ax_Bx_Cx_views
    
    #---------------------------------------------------

    Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

    stab_Ax_matrix, stab_Bx_matrix, stab_Cx_matrix = stab_Ax_Bx_Cx_matrix

    nodes_state_Idx, Bx_idxs, Cx_idxs = idxs

    stab_Ax_idxs, stab_Bx_idxs, stab_Cx_idxs = stab_idxs

    filter_stab_nodes_state_Idx, filter_stab_Bx_row_idxs, filter_stab_Cx_row_idxs = filter_stab_idxs
"""
    
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    # gens_dynamic_id_iq_pg_vh_by_vhθh, gens_nodes_ωs_τm_v_ref_view = industrial_model_pf_para

    gens_dynamic_id_iq_pg_vh_by_vhθh, gens_nodes_ωs_τm_vref_porder_view = industrial_model_pf_para

    
    #-----------------------------------------------------
    
    # ode_fun_para = ( Ax_matrix, Bx_matrix, Cx_matrix, gens_dynamic_id_iq_pg_vh_by_vhθh, gens_nodes_ωs_τm_v_ref_view  )

    ode_fun_para = ( Ax_matrix, Bx_matrix, Cx_matrix, industrial_model_pf_para  )
    
    para_model_fun = (vec_Ax_views, ode_fun_para )

    #-----------------------------------------------------

    sim_state_x0  = state
    
    sim_timespan  = sim_timespan

    #-----------------------------------------------------
    
    dyn_global_pf_options = dyn_global_pf_options
    
    counter_array  = [1]
    
    stateDiffCache = similar( sim_state_x0 )
    
    sim_fun_para   = ( netd, nodes_cb_sw, global_pf_param, counter_array, state, stateDiffCache, dyn_global_pf_options, para_update_gen_Ax_aux, para_model_fun, industrial_model_misc_Idx, industrial_model_para_aux_inputs, industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal  )
    
    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------
    #-----------------------------------------------------

    sim_func          = dynamics_industrial_model!
    
    #-----------------------------------------------------    
        
    sim_ode_func! = ODEFunction{true}(
        sim_func; mass_matrix = industrial_model_mass_matrix,
        syms = industrial_model_sym )
    
    sim_prob = ODEProblem(sim_ode_func!, sim_state_x0, sim_timespan, sim_fun_para )
    
    # cb =  industrial_model_fun_make_state_callbacks( collect( values( netd.nodes)) )

    # sim_sol = DifferentialEquations.solve(sim_prob, algr, callback = cb )

    sim_sol = DifferentialEquations.solve(sim_prob, algr )

    return sim_sol, para_net_names_labels_syms 
      
end


#-------------------------------------------------------
########################################################
#-------------------------------------------------------

#-----------------------------------------

"""
Driver for"
       `dynamics_industrial_model`
Test:
     case_IEEE_9_Bus_sauer_dynamic_plant_SM_v6_P_t2_system_matrix,

     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
"""

#-----------------------------------------

function driver_simulation_dynamics_industrial_model()


    base_dir     = joinpath(@__DIR__,"..","..","..","Simulation-Results")

    # joinpath(@__DIR__,"..","Simulation-Results")
 
    plot_tspan     = (0.0, 10.0)
    sim_timespan   = (0.0, 10.0)
    algr           = Rodas4()
    sim_model_type = "industrial-model"

    only_gen = false

    algr_name      = "rodas4"

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

        sim_sol, net_names_labels_syms = simulate_a_dynamics_industrial_model(
            ; dynamics_case = a_dynamics_case,
            only_gen = only_gen,
            sim_timespan  = sim_timespan,
            algr = algr,
            dyn_global_pf_options = dyn_global_pf_options,
            sta_global_pf_options = sta_global_pf_options )
        
        
        a_list_sim_plots = make_plots_for_industrial_model(
            ; case_fun = a_dynamics_case,
            sim_model_type = sim_model_type,
            sim_sol = sim_sol,
            para_net_names_labels_syms = net_names_labels_syms,
            tspan = plot_tspan,
            base_dir = base_dir,
            algr_name = algr_name )
        
        push!( list_sim_plots, a_list_sim_plots )
    end

    return nothing
    
    # return list_sim_plots
        
end


"""
list_sim_plots = driver_simulation_dynamics_industrial_model()

"""

#----------------------------------------------------
#----------------------------------------------------
 
function simulate_a_dynamics_industrial_model_with_or_no_controllers(
    ; dynamics_case = nothing,
    only_gen = nothing,
    sim_timespan  = nothing,
    algr = algr,
    dyn_global_pf_options = nothing ,
    sta_global_pf_options = nothing,
    with_cb = false)

    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------    

    non_gens_nodes        = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #-----------------------------------------------------
    
    network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names,Trans_bus_names = make_case_buses_names(  netd.nodes )
    net_class_names = (;network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names, Trans_bus_names )
    
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_stab_states_label, gens_nodes_algebraic_and_states_labels = generate_industrial_model_labels(; nodes =  netd.nodes, no_control_device = only_gen )

    net_states_and_var_labels = (; net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_stab_states_label, gens_nodes_algebraic_and_states_labels )

    #-----------------------------------------------------

    industrial_model_sym, industrial_model_mass_matrix = generate_industrial_model_sym_and_mass_matrix(; nodes = netd.nodes, no_control_device = only_gen )

    #-----------------------------------------------------

    para_net_names_labels_syms = (; net_class_names, net_states_and_var_labels, industrial_model_sym )

    #-----------------------------------------------------

   industrial_model_pf_sys_param_sys_views_sys_industrial_model_init =  get_industrial_model_pf_param_views_and_init_with_or_no_controller( netd; sta_global_pf_options..., only_gen=only_gen )

    #-----------------------------------------------------

    nodes_cb_sw, state, global_pf_param, named_tup_pf_result, industrial_model_misc_Idx, para_update_gen_Ax_aux, industrial_model_para_aux_inputs, industrial_model_pf_para, industrial_model_ωs_τm_vref_vhθh_idq , industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal  = industrial_model_pf_sys_param_sys_views_sys_industrial_model_init


    # -----------------------------------------------------

    pf_net_param, sd_pf_views, _ = global_pf_param

    pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ = pf_param_views

    _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views = sd_pf_views

    #-----------------------------------------------------
    #-----------------------------------------------------


    if only_gen == false
        
        #-----------------------------------------------------
        # Stability
        #----------------------------------------------------- 

        Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs = create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        vec_Ax_views, vec_Bx_views, vec_Cx_views  = Ax_Bx_Cx_views
        Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

        nodes_state_Idx, Bx_idxs, Cx_idxs = idxs

        init_plants_system_matrices_views!((vec_Ax_views, vec_Bx_views, vec_Cx_views), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = nothing )

        #-----------------------------------------------------
        #-----------------------------------------------------

        ode_fun_para = ( Ax_matrix, Bx_matrix, Cx_matrix, industrial_model_pf_para  )

        para_model_fun = (vec_Ax_views, ode_fun_para )

    else
        
        #-----------------------------------------------------
        # Stability
        #-----------------------------------------------------  
        #-----------------------------------------------------
        
        Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs = create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        vec_Ax_views, vec_Bx_views, vec_Cx_views, vec_Ax_τm_vf_views = Ax_Bx_Cx_views
        Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix = Ax_Bx_Cx_matrix
        nodes_state_Idx, Bx_idxs, Cx_idxs, τm_vf_idxs = idxs

        init_plants_system_matrices_views!((vec_Ax_views, vec_Bx_views, vec_Cx_views), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = vec_Ax_τm_vf_views )

        #-----------------------------------------------------

        # update_gens_τm_vf!( vec_τm_vf_views, gens_nodes_τm_vf  )

        #-----------------------------------------------------

        ode_fun_para = ( Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix, industrial_model_pf_para  )

        para_model_fun = (vec_Ax_views, ode_fun_para )

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
    #-----------------------------------------------------

    sim_func = dynamics_industrial_model_with_or_no_controller!

    # dynamics_industrial_model!

    #-----------------------------------------------------    

    sim_ode_func! = ODEFunction{true}(
        sim_func; mass_matrix = industrial_model_mass_matrix,
        syms = industrial_model_sym )

    sim_prob = ODEProblem(sim_ode_func!,
                          sim_state_x0,
                          sim_timespan,
                          sim_fun_para )

    if with_cb == true
        
        cb = industrial_model_fun_make_state_callbacks(
            collect( values( netd.nodes)) )

        sim_sol = DifferentialEquations.solve(
            sim_prob, algr, callback = cb )
        
    else
        sim_sol = DifferentialEquations.solve(
            sim_prob, algr  )

    end
    
    return sim_sol, para_net_names_labels_syms 
      
end


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

    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)
    algr                  = Rodas4()
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
    
    list_dynamics_case = [
        case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix,
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
  ]

    list_sim_plots = [ ]
    
    for a_dynamics_case in list_dynamics_case
        
        sim_sol, net_names_labels_syms =  simulate_a_dynamics_industrial_model_with_or_no_controllers(
            ; dynamics_case         = a_dynamics_case,
              only_gen              = only_gen,
              sim_timespan          = sim_timespan,
              algr                  = algr,
              dyn_global_pf_options = dyn_global_pf_options,
            sta_global_pf_options = sta_global_pf_options,
        with_cb  = with_cb )
        
        a_list_sim_plots   = make_plots_for_industrial_model(
            ; case_fun     = a_dynamics_case,
            sim_model_type = sim_model_type,
            sim_sol        = sim_sol,
            para_net_names_labels_syms = net_names_labels_syms,
            tspan          = plot_tspan,
            base_dir       = base_dir,
            algr_name      = algr_name )
        
        push!( list_sim_plots, a_list_sim_plots )
    end
    
    return list_sim_plots
    
    # return nothing

        
end



"""

list_sim_plots = driver_simulation_dynamics_industrial_model_with_or_no_controllers()

"""



function diagnosis_simulate_a_dynamics_industrial_model_with_or_no_controllers(
    ; dynamics_case = nothing,
    only_gen = nothing,
    sim_timespan  = nothing,
    algr = algr,
    dyn_global_pf_options = nothing ,
    sta_global_pf_options = nothing,
    with_cb  = false)

    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------    

    non_gens_nodes = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #-----------------------------------------------------
    
    network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names,Trans_bus_names = make_case_buses_names(  netd.nodes )
    net_class_names = (;network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names, Trans_bus_names )
    
    #-----------------------------------------------------
    #-----------------------------------------------------

    if only_gen == false

        net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_stab_states_label, gens_nodes_algebraic_and_states_labels = generate_industrial_model_labels(; nodes =  netd.nodes )

        net_states_and_var_labels = (; net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_stab_states_label, gens_nodes_algebraic_and_states_labels )

        #-----------------------------------------------------

        industrial_model_sym, industrial_model_mass_matrix = generate_industrial_model_sym_and_mass_matrix(; nodes = netd.nodes )

        #-----------------------------------------------------

        para_net_names_labels_syms = (; net_class_names, net_states_and_var_labels, industrial_model_sym )

        #-----------------------------------------------------

       industrial_model_pf_param_views_and_init_with_or_no_controller =  get_industrial_model_pf_param_views_and_init_with_or_no_controller( netd; sta_global_pf_options..., only_gen=only_gen )

        #-----------------------------------------------------

        nodes_cb_sw, state, global_pf_param, named_tup_pf_result, industrial_model_misc_Idx, para_update_gen_Ax_aux, industrial_model_para_aux_inputs, industrial_model_pf_para, industrial_model_ωs_τm_vref_vhθh_idq , industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal  = industrial_model_pf_param_views_and_init_with_or_no_controller

        # -----------------------------------------------------

        pf_net_param, sd_pf_views, _ = global_pf_param

        pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

        _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ = pf_param_views

        _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views = sd_pf_views
        
        #-----------------------------------------------------
        # Stability
        #----------------------------------------------------- 

        Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs = create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        vec_Ax_views, vec_Bx_views, vec_Cx_views  = Ax_Bx_Cx_views
        Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

        nodes_state_Idx, Bx_idxs, Cx_idxs = idxs

        init_plants_system_matrices_views!((vec_Ax_views, vec_Bx_views, vec_Cx_views), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = nothing )

        #-----------------------------------------------------
        #-----------------------------------------------------

        ode_fun_para = ( Ax_matrix, Bx_matrix, Cx_matrix, industrial_model_pf_para  )

        para_model_fun = (vec_Ax_views, ode_fun_para )

    else

        net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_stab_states_label, gens_nodes_algebraic_and_states_labels = generate_industrial_model_labels(; nodes =  netd.nodes, no_control_device =true )

        net_states_and_var_labels = (; net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_stab_states_label, gens_nodes_algebraic_and_states_labels )

        #-----------------------------------------------------

        industrial_model_sym, industrial_model_mass_matrix = generate_industrial_model_sym_and_mass_matrix(; nodes = netd.nodes, no_control_device =true )

        #-----------------------------------------------------

        para_net_names_labels_syms = (; net_class_names, net_states_and_var_labels, industrial_model_sym )

        #-----------------------------------------------------

       industrial_model_pf_param_views_and_init_with_or_no_controller =  get_industrial_model_pf_param_views_and_init_with_or_no_controller( netd; sta_global_pf_options..., only_gen=only_gen )

        #-----------------------------------------------------

        nodes_cb_sw, state, global_pf_param, named_tup_pf_result, industrial_model_misc_Idx, para_update_gen_Ax_aux, industrial_model_para_aux_inputs, industrial_model_pf_para, industrial_model_ωs_τm_vref_vhθh_idq , industrial_model_dyn_pf_up_para, industrial_model_idq_pf_cal  = industrial_model_pf_param_views_and_init_with_or_no_controller

        # -----------------------------------------------------

        pf_net_param, sd_pf_views, _ = global_pf_param

        pf_net, pf_idx_and_state, pf_param_views, _, _, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

        _, ra_Xd_dash_Xq_dash_view, _, _, _, _, _ = pf_param_views

        _, _, _, _, gen_nodes_δ_ω_ed_dash_eq_dash_views = sd_pf_views

        #-----------------------------------------------------
        
        #-----------------------------------------------------
        # Stability
        #-----------------------------------------------------  
        #-----------------------------------------------------
        
        Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs = create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        vec_Ax_views, vec_Bx_views, vec_Cx_views, vec_Ax_τm_vf_views = Ax_Bx_Cx_views
        Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix = Ax_Bx_Cx_matrix
        nodes_state_Idx, Bx_idxs, Cx_idxs, τm_vf_idxs = idxs

        init_plants_system_matrices_views!((vec_Ax_views, vec_Bx_views, vec_Cx_views), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = vec_Ax_τm_vf_views )

        #-----------------------------------------------------

        # update_gens_τm_vf!( vec_τm_vf_views, gens_nodes_τm_vf  )

        #-----------------------------------------------------

        ode_fun_para = (; Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix, industrial_model_pf_para  )

        para_model_fun = (; vec_Ax_views, ode_fun_para )

    end
    
    #-----------------------------------------------------

    sim_state_x0  = state

    sim_timespan  = sim_timespan

    #-----------------------------------------------------

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

"""
Test:

`diagnosis_simulate_a_dynamics_industrial_model_with_or_no_controllers`


"""

function driver_diagnosis_simulate_a_dynamics_industrial_model_with_or_no_controllers()

    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)
    algr                  = Rodas4()
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
        algr                  = algr,
        dyn_global_pf_options = dyn_global_pf_options,
        sta_global_pf_options = sta_global_pf_options,
        with_cb  = with_cb )


    x0 = sim_fun_para.state
    x0_ps_idx = sim_fun_para.industrial_model_misc_Idx.industrial_model_pure_states_Idx
    x0_gen_ps_idx = sim_fun_para.industrial_model_misc_Idx.industrial_model_each_gen_nodes_pure_states_idx_in_state

    t_vec_Ax, t_vec_Bx, t_vec_Cx = Ax_Bx_Cx

    vec_Ax = sim_fun_para.para_model_fun.vec_Ax_views

    Ax_m = sim_fun_para.para_model_fun.ode_fun_para.Ax_matrix
    Bx_m = sim_fun_para.para_model_fun.ode_fun_para.Bx_matrix
    Cx_m = sim_fun_para.para_model_fun.ode_fun_para.Cx_matrix

    id_iq_pg_vh = sim_fun_para.para_model_fun.ode_fun_para.industrial_model_pf_para.gens_dynamic_id_iq_pg_vh_by_vhθh_view
    
    ωs_τm_vref_porder = sim_fun_para.para_model_fun.ode_fun_para.industrial_model_pf_para.gens_nodes_ωs_τm_vref_porder_view

    t_dx = Ax_m * x0[ x0_ps_idx ] + Bx_m * [id_iq_pg_vh...;] + Cx_m * [ ωs_τm_vref_porder...; ]

    if only_gen == true
        
        t_vec_Ax, t_vec_Bx, t_vec_Cx, t_vec_Ax_τm_vf = Ax_Bx_Cx
        
        Ax_τm_vf_m = sim_fun_para.para_model_fun.ode_fun_para.Ax_τm_vf_matrix
        vec_τm_vf = sim_fun_para.para_model_fun.ode_fun_para.industrial_model_pf_para.vec_τm_vf_views

        t2_dx = Ax_m * x0[ x0_ps_idx ] + Bx_m * [id_iq_pg_vh...;] + Ax_τm_vf_m * [ vec_τm_vf...;]
    end
    
    
    return nothing
    
    # return sim_fun_para
    
        
end

"""

t_sim_fun_para = driver_diagnosis_simulate_a_dynamics_industrial_model_with_or_no_controllers()

"""


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

    im_vars_indices_in_system, pure_states_Idx_in_system, im_algebraic_vars_Idx_in_system, ur_ui_Idx_in_system, im_vars_and_ur_ui_Idx_in_system, im_vars_Idx_in_state, nodes_ur_ui_Idx_in_state, im_state_Idx, each_gens_im_vars_Idx_in_state, net_to_im_idx_conversion_dict = get_im_indices_and_conversion_dict( netd  )
    
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

    pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx


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

    idq_wt_pad_view = view(idq_wt_pad, 1:length( uh_state ) )

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


"""

Test driver for:

`v2_get_im_model_pf_param_views_and_init`

"""
function driver_v2_get_im_model_pf_param_views_and_init( )

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd  = NetworkData( dynamics_case()... )

    only_gen = false

    im_vars_indices_in_system, pure_states_Idx_in_system, im_algebraic_vars_Idx_in_system, ur_ui_Idx_in_system, im_vars_and_ur_ui_Idx_in_system, im_vars_Idx_in_state, nodes_ur_ui_Idx_in_state, im_state_Idx, each_gens_im_vars_Idx_in_state, net_to_im_idx_conversion_dict, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, nodes_ω_ed_dash_eq_dash_Idxs_in_state, nodes_δ_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view, gen_nodes_dyn_param_view, gen_nodes_sub_dyn_param_view, pf_net_param, ra_Xd_dash_Xq_dash_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, im_vars_view_in_state, state_view, pf_state, nodes_u_view, nodes_pf_U_view, x0_vh_θh, working_vh_θh_view, red_vh_θh_0_view, mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, branches_name, nodes_name, named_tup_pf_result, bus_dict_init, branch_dict_init, state, gens_idx, nodes_u_Idx = v2_get_im_model_pf_param_views_and_init(
        netd;
        maxiter=40,
        ftol=1000*eps() ,
        xtol=1000*eps(),
        init_pf = true,
        with_δ_ed_eq = false,
        only_gen = only_gen )

    return nothing
    
end

#--------------------------------------------

"""

driver_v2_get_im_model_pf_param_views_and_init()


"""

#--------------------------------------------

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

im_vars_indices_in_system, pure_states_Idx_in_system, im_algebraic_vars_Idx_in_system, ur_ui_Idx_in_system, im_vars_and_ur_ui_Idx_in_system, im_vars_Idx_in_state, nodes_ur_ui_Idx_in_state, im_state_Idx, each_gens_im_vars_Idx_in_state, net_to_im_idx_conversion_dict, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state, nodes_ω_ed_dash_eq_dash_Idxs_in_state, nodes_δ_ed_dash_eq_dash_Idxs_in_state, gen_nodes_ra_Xd_dash_Xq_dash_view, gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view, gen_nodes_dyn_param_view, gen_nodes_sub_dyn_param_view, pf_net_param, ra_Xd_dash_Xq_dash_view, gen_nodes_δ_ω_ed_dash_eq_dash_views, im_vars_view_in_state, state_view, pf_state, nodes_u_view, nodes_pf_U_view, x0_vh_θh, working_vh_θh_view, red_vh_θh_0_view, mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, branches_name, nodes_name, named_tup_pf_result, bus_dict_init, branch_dict_init, state, gens_idx, nodes_u_Idx  = v2_get_im_model_pf_param_views_and_init(
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

    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)

    #----------------------------------------------------

    gens_vh_θh = get_gens_vh_θh( nodes_pf_U_view, gens_idx )

    gens_vh_θh_view = @view gens_vh_θh[ 1:length(gens_vh_θh ) ]

    #---------------------------------------------------

    # gen_uh  = (named_tup_pf_result.Vbus)[ gen_idx ]

    idq_wt_pad_view[gens_idx] .=  [
        industrial_model_get_dynamic_idq_θ_π_vhθh(
            vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
        for ( vh_θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash ) in zip(
            gens_vh_θh_view ,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            ra_Xd_dash_Xq_dash_view[gens_idx] ) ]

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
        gens_vh_θh, idq_wt_pad )

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
    
    #---------------------------------------------------

end

# 


"""

Test driver for:

` v2_get_im_sys_pf_param_views_and_init`

"""
function driver_v2_get_im_sys_pf_param_views_and_init( )

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    netd  = NetworkData( dynamics_case()... )

    only_gen = false

    nodes_cb_sw, state, global_pf_param, named_tup_pf_result, im_misc_Idx, para_update_gen_Ax_aux, im_para_aux_inputs, im_pf_para, im_ωs_ωref_vref_vhθh_idq_etc, im_dyn_pf_up_para, im_idq_pf_cal = v2_get_im_sys_pf_param_views_and_init( netd; maxiter=40, ftol=1000*eps() , xtol=1000*eps(), init_pf = true, with_δ_ed_eq = false,  only_gen = only_gen )

               
    return nothing
    
end

#--------------------------------------------

"""

driver_v2_get_im_sys_pf_param_views_and_init()

"""
#--------------------------------------------

function v2_simulate_a_dynamics_im_model(
    ; dynamics_case = dynamics_case,
    only_gen = false,
    sim_timespan  = sim_timespan,
    algr          = algr,
    dyn_global_pf_options = dyn_global_pf_options,
    sta_global_pf_options = sta_global_pf_options,
    with_cb = false
        )

    # ----------------------------------------------------
    # ----------------------------------------------------
    # Diagnosis
    # ---------------------------------------------------- 
    # ----------------------------------------------------

    # dynamics_case = case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
    
    #-----------------------------------------------------
    #----------------------------------------------------- 
    
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    # Dyn_Nodes, Dyn_Branches = dynamics_case()

    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------
    #-----------------------------------------------------

    non_gens_nodes        = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #-----------------------------------------------------
    #-----------------------------------------------------
    
    network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names,Trans_bus_names = make_case_buses_names(  netd.nodes )
    net_class_names = (;network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names, Trans_bus_names )
    
    #-----------------------------------------------------

    net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_im_algebraic_vars_labels, gens_nodes_im_vars_labels = generate_im_model_labels(; nodes =  netd.nodes )

    im_net_states_and_var_labels = (; net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_im_algebraic_vars_labels, gens_nodes_im_vars_labels )

    #-----------------------------------------------------
    
    im_sym, im_mass_matrix = generate_im_sym_and_mass_matrix(; nodes = netd.nodes )

    #-----------------------------------------------------

    para_net_names_labels_syms = (; net_class_names, im_net_states_and_var_labels, im_sym )
    
    #-----------------------------------------------------
    
   im_sys_pf_param_views_and_init =  v2_get_im_sys_pf_param_views_and_init( netd; sta_global_pf_options... )

    #-----------------------------------------------------
    #-----------------------------------------------------

    nodes_cb_sw, state, global_pf_param, named_tup_pf_result, im_misc_Idx, para_update_gen_Ax_aux, im_para_aux_inputs, im_pf_para, im_ωs_ωref_vref_vhθh_idq_etc, im_dyn_pf_up_para, im_idq_pf_cal  = im_sys_pf_param_views_and_init
    
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
    #-----------------------------------------------------

    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(
            gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix =
        Ax_Bx_Cx_matrix

    nodes_state_Idx, Bx_idxs, Cx_idxs, id_iq_ph_vh_idxs, ω_ref_ωref_v_ref_idxs =
        idxs

    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views),
        gens_nodes_collection )
    
    #-----------------------------------------------------
    #-----------------------------------------------------

    ode_fun_para = (
        ; Ax_matrix,
        Bx_matrix,
        Cx_matrix,
        im_pf_para )
    
    para_model_fun = (
        ; vec_Ax_views,
        ode_fun_para )

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
        im_misc_Idx,
        im_para_aux_inputs,
        im_dyn_pf_up_para,
        im_idq_pf_cal,
        only_gen  )
    
    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------
    #-----------------------------------------------------

    sim_func = dynamics_im_model!
    
    #-----------------------------------------------------
    
    #-----------------------------------------------------
        
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
            sim_prob, algr, callback = cb )
        
    else
        
        sim_sol = DifferentialEquations.solve(
            sim_prob, algr )
    end
    
    #-----------------------------------------------------
    #-----------------------------------------------------
    
    return sim_sol, para_net_names_labels_syms 
      
end


"""
Test:

`v2_simulate_a_dynamics_im_model`


"""

function driver_v2_simulate_a_dynamics_im_model()

    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)
    algr                  = Rodas4()
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
            ; dynamics_case = dynamics_case,
            only_gen = false,
            sim_timespan  = sim_timespan,
            algr = algr,
            dyn_global_pf_options = dyn_global_pf_options,
            sta_global_pf_options = sta_global_pf_options,
        with_cb = with_cb)


    a_list_sim_plots = make_plots_for_industrial_model(
        ; case_fun = dynamics_case,
        sim_model_type = sim_model_type,
        sim_sol = sim_sol,
        para_net_names_labels_syms = para_net_names_labels_syms,
        tspan = plot_tspan,
        base_dir = base_dir,
        algr_name = algr_name )

    push!( list_sim_plots, a_list_sim_plots )    

    # return nothing
    return list_sim_plots
end

#--------------------------------------------
"""

t_list_sim_plots = driver_v2_simulate_a_dynamics_im_model()

"""
#--------------------------------------------

function driver_v2_simulate_dynamics_im_model_cases()

    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)
    algr                  = Rodas4()
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
    
    list_dynamics_case = [
        case_IEEE_9_Bus_dynamic_plants_v6_P_t2_rscad,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer,
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix,
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P
  ]

    list_sim_plots = [ ]
    
    for a_dynamics_case in list_dynamics_case

        sim_sol, net_names_labels_syms = v2_simulate_a_dynamics_im_model(
            ; dynamics_case = a_dynamics_case,
            only_gen = only_gen,
            sim_timespan  = sim_timespan,
            algr = algr,
            dyn_global_pf_options = dyn_global_pf_options,
            sta_global_pf_options = sta_global_pf_options,
            with_cb = with_cb)

        
        a_list_sim_plots   = make_plots_for_industrial_model(
            ; case_fun     = a_dynamics_case,
            sim_model_type = sim_model_type,
            sim_sol        = sim_sol,
            para_net_names_labels_syms = net_names_labels_syms,
            tspan          = plot_tspan,
            base_dir       = base_dir,
            algr_name      = algr_name )
        
        push!( list_sim_plots, a_list_sim_plots )
    end
    
    # return list_sim_plots
    
    # return nothing

        
end

"""

driver_v2_simulate_dynamics_im_model_cases()

"""


#--------------------------------------------

function dignosis_v2_simulate_a_dynamics_im_model()


    base_dir              = joinpath(
        @__DIR__,"..","..","..",
        "Simulation-Results")
 
    plot_tspan            = (0.0, 10.0)
    sim_timespan          = (0.0, 10.0)
    algr                  = Rodas4()
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
    

    netd  = NetworkData( dynamics_case()... )
    
    #-----------------------------------------------------
    #-----------------------------------------------------

    non_gens_nodes        = get_non_gens_nodes( netd.nodes )

    gens_nodes_collection = get_gens_nodes( netd.nodes )

    #-----------------------------------------------------
    #-----------------------------------------------------
    
    network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names,Trans_bus_names = make_case_buses_names( netd.nodes )

    net_class_names = (; network_bus_names, non_gens_bus_names, gens_bus_names, Loads_bus_names, Trans_bus_names )
    
    #-----------------------------------------------------

    net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_im_algebraic_vars_labels, gens_nodes_im_vars_labels = generate_im_model_labels(; nodes =  netd.nodes )

    im_net_states_and_var_labels = (; net_bus_volts_labels, gens_nodes_pure_states_labels, gens_nodes_im_algebraic_vars_labels, gens_nodes_im_vars_labels )

    #-----------------------------------------------------
    
    im_sym, im_mass_matrix = generate_im_sym_and_mass_matrix(; nodes = netd.nodes )

    #-----------------------------------------------------

    para_net_names_labels_syms = (; net_class_names, im_net_states_and_var_labels, im_sym )
    
    #-----------------------------------------------------
    
    
   im_sys_pf_param_views_and_init =  v2_get_im_sys_pf_param_views_and_init( netd; sta_global_pf_options... )

    #-----------------------------------------------------
    #-----------------------------------------------------

    nodes_cb_sw, state, global_pf_param, named_tup_pf_result, im_misc_Idx, para_update_gen_Ax_aux, im_para_aux_inputs, im_pf_para, im_ωs_ωref_vref_vhθh_idq_etc, im_dyn_pf_up_para, im_idq_pf_cal  = im_sys_pf_param_views_and_init
    
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
    #-----------------------------------------------------

    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(  gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix =
        Ax_Bx_Cx_matrix

    nodes_state_Idx, Bx_idxs, Cx_idxs, id_iq_ph_vh_idxs, ω_ref_ωref_v_ref_idxs =
        idxs

    init_im_plants_system_matrices_views!(
        (vec_Ax_views, vec_Bx_views, vec_Cx_views),
        gens_nodes_collection )
    #-----------------------------------------------------
    #-----------------------------------------------------

    ode_fun_para = (
        ; Ax_matrix,
        Bx_matrix,
        Cx_matrix,
        im_pf_para  )
    
    para_model_fun = (
        ; vec_Ax_views,
        ode_fun_para )

    #-----------------------------------------------------

    sim_state_x0  = state
    
    sim_timespan  = sim_timespan

    #-----------------------------------------------------
    
    dyn_global_pf_options = dyn_global_pf_options
    
    counter_array  = [1]
    
    stateDiffCache = similar( sim_state_x0 )
    
    sim_fun_para   = (
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
        im_idq_pf_cal,
        only_gen  )
    
    # #-----------------------------------------------------
    # # dynamics simulation
    # #-----------------------------------------------------
    # #-----------------------------------------------------

    # sim_func = dynamics_im_model!
    # #-----------------------------------------------------
    
    # #-----------------------------------------------------
        
    # sim_ode_func! = ODEFunction{true}(
    #     sim_func; mass_matrix = im_mass_matrix,
    #     syms = im_sym )
    
    # sim_prob = ODEProblem(sim_ode_func!,
    #                       sim_state_x0,
    #                       sim_timespan,
    #                       sim_fun_para )


    # if with_cb == true
        
    #     cb = industrial_model_fun_make_state_callbacks(
    #         collect( values( netd.nodes)) )

    #     sim_sol = DifferentialEquations.solve(
    #         sim_prob, algr, callback = cb )
        
    # else
    #     sim_sol = DifferentialEquations.solve(
    #         sim_prob, algr  )

    # end
    

    # #-----------------------------------------------------
    
    # return sim_sol, para_net_names_labels_syms 

    
    nd, nodes_cb_sw, global_pf_param, counter_array, state, stateDiffCache, dyn_global_pf_options, para_update_gen_Ax_aux, para_model_fun, im_misc_Idx, im_para_aux_inputs, im_dyn_pf_up_para, im_idq_pf_cal = sim_fun_para
    
    #--------------------------------------------

    nodes_u_Idx, nodes_u_Idx_in_ranges, im_vars_Idx_in_state, each_gens_im_vars_Idx_in_state, gens_nodes_collection = im_misc_Idx

    #--------------------------------------------
    
    vec_Ax_views, ode_fun_para = para_model_fun

    Ax_matrix, Bx_matrix, Cx_matrix, im_pf_para = ode_fun_para

    # gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_ωref0_vref0_porder0_view
    
    gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_ωref0_vref0_porder0_view = im_pf_para
        
    #--------------------------------------------

    pf_net_param, sd_pf_views, mismatch = global_pf_param

    working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view = sd_pf_views
    
    #--------------------------------------------    

    im_vars_view_in_state, _, _ = para_update_gen_Ax_aux
    
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

    tt = 0.0
    
    ode_im_model_func!(
        dx,
        x,
        ( im_vars_Idx_in_state,
          ode_fun_para),
        tt)

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

    t_dx = t_Ax_m * x0[ x0_ps_idx ] + t_Bx_m * [id_iq_pg_vh...;] + t_Cx_m * [ ωs_ω_ref_vref_porder...; ]

    
    # #--------------------------------------------
    # # Ode gens pure state
    # #--------------------------------------------

    # tt = 0.0
    
    # ode_im_model_func!(
    #     dx,
    #     x,
    #     ( im_vars_Idx_in_state, ode_fun_para),
    #     tt)
    
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
    #     dx, x, im_para_aux_inputs, tt)
    
    # #--------------------------------------------
    # # gens terminal voltages
    # #--------------------------------------------

    # gens_im_voltage_terminal_func!(
    #     dx, x, im_para_aux_inputs, tt)
    
    # #--------------------------------------------
    # #--------------------------------------------    
    
    # return nothing
    
end



#-------------------------------------------------------
# one model
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


    Ax_Bx_Cx_views, gens_dynamic_id_iq_pg_vh_by_vhθh_view, gens_nodes_ωs_ωref0_vref0_porder0_view, node_i = ode_para

    vec_Ax_views, vec_Bx_views, vec_Cx_views  = Ax_Bx_Cx_views
        
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

function dynamics_one_im_model!(dx, x, para, t)


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

    # ode_one_im_model_func!(
    #     dx,
    #     x,
    #     ode_one_fun_para ,
    #     t)
    
    #--------------------------------------------
    # ode gens pure state
    #--------------------------------------------

    ode_one_im_model_func!(
        dx,
        x,
        ode_para ,
        t)


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


function v2_simulate_one_plant_dynamics_im_model(
    ; dynamics_case =
        dynamics_case,
    only_gen =
        false,
    sim_timespan =
        sim_timespan,
    algr =
        algr,
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
    
    # algr                  = Rodas4()
    # algr_name             = "rodas4" 
    
    algr                  = ImplicitMidpoint()    
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
    
    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

    nodes_state_Idx, Bx_idxs, Cx_idxs, id_iq_ph_vh_idxs, ω_ref_ωref_v_ref_idxs  = idxs
    
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
    # #    sim_prob, algr  )

    # sim_sol = DifferentialEquations.solve(
    #     sim_prob, algr, dt=dt  )

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

    
    vec_Ax_i = vec_Ax_views[ node_i ]

    #--------------------------------------------

    x  = sim_state_x0
    
    dx = similar(x)
    
    #--------------------------------------------

    stateDiffCache = get_tmp(stateDiffCache, x)

    stateDiffCache .= sim_state_x0
    
    #--------------------------------------------
    # update Ax, 
    #--------------------------------------------
    
    update_a_im_plant_system_matrices!(
        vec_Ax_i,
        stateDiffCache,
        plant_i )    

    #--------------------------------------------

    t_Ax_Bx_Cx_views, t_gens_dynamic_id_iq_pg_vh_by_vhθh_view, t_gens_nodes_ωs_ωref0_vref0_porder0_view, node_i = ode_para

    t_vec_Ax_views, t_vec_Bx_views, t_vec_Cx_views  = t_Ax_Bx_Cx_views
        
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
######################################################
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
         idxs) = create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )
        

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
             create_gens_nodes_aggregate_system_matrices_Idx_and_views(  gens_nodes_collection; only_gen = only_gen )
        
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


#---------------------------------------------------
####################################################
#---------------------------------------------------



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


function system_ode_im_model_func!(
    dx, x,
    system_ode_para, t;
    sim_fun_kwd_para = sim_fun_system_kwd_para )

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
     para_update_gen_Ax_aux) =
         sim_fun_kwd_para

    #-----------------------------------------------------    
    
    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ω_ref_ωref_v_ref_idxs) =
         sys_states_Idxs_and_mat_Idxs

    #-----------------------------------------------------

    # dx = similar( sim_state_x0 )
    
    dx_gens  = [ view(dx, idx)
                 for idx in
                     nodes_state_Idx ]
    
    #-----------------------------------------------------

    # x = similar( sim_state_x0 )
    
    x_gens   = [ view(x, idx)
                 for idx in
                     nodes_state_Idx ]

    #-----------------------------------------------------

    gens_ode_para = [ system_ode_para[idx] for idx in sim_fun_system_para_Idxs[1] ]

    
    #-----------------------------------------------------
    
    vec_sim_fun_kwd_para = [ (; Ax_m, Bx_m, Cx_m, stateDiffCache, sim_state_x0, plant_i, sim_fun_para_Idxs, ra_Xd_dash_Xq_dash_i ) for ( Ax_m, Bx_m, Cx_m, stateDiffCache, sim_state_x0, plant_i, sim_fun_para_Idxs, ra_Xd_dash_Xq_dash_i ) in zip( vec_Ax_views, vec_Bx_views, vec_Cx_views, stateDiffCache_gens, sim_state_x0_gens, gens_nodes_collection, gens_sim_fun_gens_para_Idxs, gen_nodes_ra_Xd_dash_Xq_dash_view  ) ]

    #-----------------------------------------------------
    
    for (a_dx, a_x, a_ode_para, a_sim_fun_kwd_para ) in zip( dx_gens, x_gens, gens_ode_para, vec_sim_fun_kwd_para )

        t_ode_one_im_model_func!(a_dx, a_x, a_ode_para, t;
            sim_fun_kwd_para = a_sim_fun_kwd_para )
    end

    return nothing
end


##

function system_flat_agg_ode_im_model_func!(
    dx, x,
    sim_fun_system_ode_flat_agg_para, t;
    sim_fun_kwd_para =
        sim_fun_system_kwd_flat_agg_para )

     (; vec_Ax_views,
      vec_Bx_views,
      vec_Cx_views,
      gen_nodes_ra_Xd_dash_Xq_dash_view,
      # ra_Xd_dash_Xq_dash_view,       
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
                 for idx in
                     nodes_state_Idx ]
    
    x_gens   = [ view(x, idx)
                 for idx in
                     nodes_state_Idx ]
    
    vec_sim_fun_kwd_para = [ (; Ax_m, Bx_m, Cx_m, stateDiffCache, sim_state_x0, plant_i, sim_fun_para_Idxs, ra_Xd_dash_Xq_dash_i ) for ( Ax_m, Bx_m, Cx_m, stateDiffCache, sim_state_x0, plant_i, sim_fun_para_Idxs, ra_Xd_dash_Xq_dash_i ) in zip( vec_Ax_views, vec_Bx_views, vec_Cx_views, stateDiffCache_gens, sim_state_x0_gens, gens_nodes_collection, gens_sim_fun_gen_para_Idxs, gen_nodes_ra_Xd_dash_Xq_dash_view  ) ]

    #-----------------------------------------------------
    
    for (a_dx, a_x, a_ode_para, a_sim_fun_kwd_para ) in zip( dx_gens, x_gens, gens_vh_θh_ωs_ωref0_vref0_porder0, vec_sim_fun_kwd_para )

        t_ode_one_im_model_func!(a_dx, a_x, a_ode_para, t;
            sim_fun_kwd_para = a_sim_fun_kwd_para )
    end

    return nothing
end



## 
function system_ode_poi_im_model_func!(
    dx, x,
    system_ode_para, t;
    poi_idx = poi_idx, 
    sim_fun_system_kwd_para =
        sim_fun_system_kwd_para_wt_poi )

    (; vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views,
     gen_nodes_ra_Xd_dash_Xq_dash_view,      
     # ra_Xd_dash_Xq_dash,
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
              nodes_vh_θh_indx_in_Idx]

    ωs_ωref0_vref0_porder0 =
        [ f_ωs_ωref0_vref0_porder0[idx]
          for idx in
              nodes_ωs_ωref0_vref0_porder0_idx_in_Idx ]

    

    
    nodes_state_Idx, Bx_idxs, Cx_idxs, id_iq_ph_vh_idxs, ω_ref_ωref_v_ref_idxs = states_and_mat_Idxs
       
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
    
    vec_sim_fun_kwd_para = [ (; Ax_m, Bx_m, Cx_m, stateDiffCache, sim_state_x0, plant_i, sim_fun_para_Idxs, ra_Xd_dash_Xq_dash_i ) for ( Ax_m, Bx_m, Cx_m, stateDiffCache, sim_state_x0, plant_i, sim_fun_para_Idxs, ra_Xd_dash_Xq_dash_i ) in zip( vec_Ax_views, vec_Bx_views, vec_Cx_views, stateDiffCache_gens, sim_state_x0_gens, gens_nodes_collection, gen_nodes_ra_Xd_dash_Xq_dash_view  ) ]

    return nothing
end





function t_ode_one_im_model_func_dxdp(
    ;
      sim_state_x0 = sim_state_x0,
     sim_fun! = t_ode_one_im_model_func! ,
     sim_fun_para = sim_fun_para,
     mass_matrix = im_mass_matrix_i,
     syms = im_sym_i,    
     sim_timespan = sim_timespan,
     sim_fun_kwd_para = t_sim_fun_kwd_para,
     alg = alg)

    #--------------------------------------------

     sim_state_x0 = sim_state_x0
     sim_fun! = t_ode_one_im_model_func!
     sim_fun_para = sim_fun_para
     mass_matrix = im_mass_matrix_i
     syms = im_sym_i    
     sim_timespan = sim_timespan
     sim_fun_kwd_para = t_sim_fun_kwd_para
     alg = algr

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
     
     t_df_dp = t_dg_dp = ForwardDiff.jacobian(g, sim_fun_para)

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

        return (; t_df_dx, t_df_dp, t_dx_dp, tt_df_dx, tt_df_dp, tt_dx_dp)

end



function s_ode_one_im_model_func_dxdp(
    ;
      sim_state_x0 = sim_state_x0,
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
    #  alg =  algr
     
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
         sim_prob, algr, dt=dt  )



     function s_f(x)

         prob = remake(sim_prob,
                       u0 = x,
                       p = poi )

         DifferentialEquations.solve(prob, algr, dt=dt,
               reltol = 1e-6, abstol = 1e-6,
               saveat = 1)[1, :]

     end


     function s_g(p)

         prob = remake(sim_prob,
                       u0 = xx,
                       p = p )

         DifferentialEquations.solve(prob, algr, dt=dt,
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


        return (; s_df_dx, s_df_dp, s_dx_dp, ss_df_dx, ss_df_dp, ss_dx_dp)
     
end



function s_oop_ode_one_im_model_func_dxdp(
    ;
      sim_state_x0 = sim_state_x0,
     sim_fun! = s_oop_ode_one_im_model_func!,
     poi = [poi...],
     poi_idx = poi_idx,
     mass_matrix = im_mass_matrix_i,
     syms = im_sym_i,    
     sim_timespan = sim_timespan,
     sim_fun_kwd_para = s_sim_fun_kwd_para,
     alg = algr)
     
    # #--------------------------------------------

    #  sim_state_x0 = sim_state_x0 
    #  sim_fun! = s_oop_ode_one_im_model_func! 
    #  poi = s_oop_poi 
    #  poi_idx = poi_idx 
    #  mass_matrix = im_mass_matrix_i 
    #  syms = im_sym_i 
    #  sim_timespan = sim_timespan 
    #  sim_fun_kwd_para = s_sim_fun_kwd_para 
    #  alg =  algr
     
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
         sim_prob, algr, dt=dt  )

     function s_oop_f(x)

         prob = remake(sim_prob,
                       u0 = x,
                       p = poi )

         DifferentialEquations.solve(prob, algr, dt=dt,
               reltol = 1e-6, abstol = 1e-6,
               saveat = 1)[1, :]

     end


     function s_oop_g(p)

         prob = remake(sim_prob,
                       u0 = xx,
                       p = p )

         DifferentialEquations.solve(prob, algr, dt=dt,
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

        return  (; so_df_dx, so_df_dp, so_dx_dp, soo_df_dx, soo_df_dp, soo_dx_dp )
     
end



function t_v2_simulate_one_plant_dynamics_im_model(
    ; dynamics_case =
        dynamics_case,
    only_gen =
        false,
    sim_timespan =
        sim_timespan,
    algr =
        algr,
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
    
    # algr                  = Rodas4()
    # algr_name             = "rodas4" 
    
    algr                  = ImplicitMidpoint()    
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
    
    #-----------------------------------------------------
    
    (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para

    gens_vh_θh_view =
        im_para_aux_inputs.gens_vh_θh_view
    
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
    
    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, states_and_mat_Idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

    nodes_state_Idx, Bx_idxs, Cx_idxs, id_iq_ph_vh_idxs, ω_ref_ωref_v_ref_idxs  = states_and_mat_Idxs
    
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

    Ax_m = vec_Ax_views[ node_i ]

    Bx_m = vec_Bx_views[ node_i ]

    Cx_m = vec_Cx_views[ node_i ]
    
    #-----------------------------------------------------
    
    plant_i =
        gens_nodes_collection[ node_i ]      

    #-----------------------------------------------------

    plant_i_eigvalues =
        eigvals( Matrix( Ax_m ))

    #-----------------------------------------------------

    ra_Xd_dash_Xq_dash_i =
        ra_Xd_dash_Xq_dash_view[ node_i ]

    #-----------------------------------------------------
    
    id_iq_pg_vh_i =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view[ node_i ]

    ωs_ωref0_vref0_porder0_i =
        gens_nodes_ωs_ωref0_vref0_porder0_view[ node_i ]

    gens_vh_θh_i = [ gens_vh_θh_view[ node_i ]...;]

    #-----------------------------------------------------

    dims_gens_vh_θh_ωs_ωref0_vref0_porder0_i =
        length.([ gens_vh_θh_i, ωs_ωref0_vref0_porder0_i ])

    _,_, vh_θh_ωs_ωref0_vref0_porder0_i_Idx =
        create_size_offset_Idx(
            dims_gens_vh_θh_ωs_ωref0_vref0_porder0_i,;
            counter = 0 )

    vh_θh_i_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_i_Idx[1]

    ωs_ωref0_vref0_porder0_i_Idx =
        vh_θh_ωs_ωref0_vref0_porder0_i_Idx[2]
        
    kwd_dyn_vh_θh_ωs_ωref0_vref0_etc_i =
        (; gens_vh_θh_i,
         ωs_ωref0_vref0_porder0_i )

    #-----------------------------------------------------
    
    poi_dyn_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )

    pois_dyn =
        kwd_dyn_vh_θh_ωs_ωref0_vref0_etc_i    

    #-----------------------------------------------------

    sim_fun_para_Idxs =
        (; vh_θh_i_Idx,
         ωs_ωref0_vref0_porder0_i_Idx )
        
    sim_fun_para  =
        vcat( gens_vh_θh_i,
             ωs_ωref0_vref0_porder0_i )
    
    #-----------------------------------------------------

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

    #-----------------------------------------------------

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

    #-----------------------------------------------------
    
    sim_timespan =
        sim_timespan

    #-----------------------------------------------------

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
    
    #-----------------------------------------------------
    
    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------
            
    im_sym_i =
        im_sym[ gen_node_i_idx ]

    im_mass_matrix_i =
        im_mass_matrix[gen_node_i_idx,
                       gen_node_i_idx]

    #-----------------------------------------------------
    # case 1
    #-----------------------------------------------------

    t_sim_func! = t_ode_one_im_model_func!
    
    #-----------------------------------------------------   
        
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
    
    #-----------------------------------------------------
    # simulate 
    #-----------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, algr  )

    t_sim_sol = DifferentialEquations.solve(
        t_sim_prob, algr, dt=dt  )

    #-----------------------------------------------------
    # plot
    #-----------------------------------------------------
    
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
        t_sens_prob , algr, dt=dt  )


    #-----------------------------------------------------
    # case 2
    #-----------------------------------------------------

    s_sim_func! = s_ode_one_im_model_func!
    
    #-----------------------------------------------------   
        
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
    #    sim_prob, algr  )

    s_sim_sol = DifferentialEquations.solve(
        s_sim_prob, algr, dt=dt  )

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
        s_sa_sens_prob , algr, dt=dt  )

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
                         

#---------------------------------------------------
####################################################
#---------------------------------------------------



function t_v2_simulate_sysyem_dynamic_im_model(
    ; dynamics_case =
        dynamics_case,
    only_gen =
        false,
    sim_timespan =
        sim_timespan,
    algr =
        algr,
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
    
    # algr                  = Rodas4()
    # algr_name             = "rodas4" 
    
    algr                  = ImplicitMidpoint()    
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
    
    #-----------------------------------------------------
    
    (; gens_dynamic_id_iq_pg_vh_by_vhθh_view,
     gens_nodes_ωs_ωref0_vref0_porder0_view) =
         im_pf_para

    gens_vh_θh_view =
        im_para_aux_inputs.gens_vh_θh_view
    
    # -----------------------------------------------------

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
    
    #-----------------------------------------------------
    # Stability
    #-----------------------------------------------------  
    
    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, states_and_mat_Idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

    vec_Ax_views, vec_Bx_views, vec_Cx_views  =
        Ax_Bx_Cx_views
    
    Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

    nodes_state_Idx, Bx_idxs, Cx_idxs, id_iq_ph_vh_idxs, ω_ref_ωref_v_ref_idxs  = states_and_mat_Idxs

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

    #-----------------------------------------------------

    Ax_m = vec_Ax_views

    Bx_m = vec_Bx_views

    Cx_m = vec_Cx_views
    
    #-----------------------------------------------------    

    # plant_i_eigvalues =
    #     eigvals( Matrix( vec_Ax_views[i] ))

    #-----------------------------------------------------
    #-----------------------------------------------------

    """

    #-----------------------------------------------------
    # case 0
    #-----------------------------------------------------
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i]) for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ]

    gens_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gens_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------------
    
   dims_vh_θh_ωs_ωref0_vref0_porder0 =  length.( [ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ] )

    system_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), [dims_vh_θh_ωs_ωref0_vref0_porder0]  )

    vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(system_size_offset_Idx)
    
    sim_fun_system_para_Idxs =
        vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------------

    system_ode_para = [[ [ vh_θh_i...; ωs_ωref0_vref0_porder0_i...]  for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ]...;]
    
    #-----------------------------------------------------

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
    #    sim_prob, algr  )

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, algr, dt=dt  )

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
        sens_prob , algr, dt=dt  )

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
    
   dims_gens_vh_θh_ωs_ωref0_vref0_porder0 = [ length.([ vh_θh_i, ωs_ωref0_vref0_porder0_i]) for (vh_θh_i, ωs_ωref0_vref0_porder0_i ) in zip( gens_vh_θh_view,  gens_nodes_ωs_ωref0_vref0_porder0_view) ]

    gens_size_offset_Idx =  map((dim) -> create_size_offset_Idx(dim; counter = 0 ), dims_gens_vh_θh_ωs_ωref0_vref0_porder0 )

    gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs =
        third.(gens_size_offset_Idx)
    
    gens_sim_fun_gen_para_Idxs =
        gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs
    
    #-----------------------------------------------------
    
    
    id_iq_pg_vh =
        gens_dynamic_id_iq_pg_vh_by_vhθh_view

    ωs_ωref0_vref0_porder0 =
        gens_nodes_ωs_ωref0_vref0_porder0_view
        
    gens_vh_θh = gens_vh_θh_view
    
    #-----------------------------------------------------

    dims_nodes_id_iq_pg_vh =
        length.( id_iq_pg_vh )

    _,_, nodes_id_iq_pg_vh_idx_in_Idx =
        create_size_offset_Idx(
            dims_nodes_id_iq_pg_vh;
        counter = 0)

    #-----------------------------------------------------

    dims_nodes_ωs_ωref0_vref0_porder0 =
        length.( ωs_ωref0_vref0_porder0 )

    _,_, nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        create_size_offset_Idx(
            dims_nodes_ωs_ωref0_vref0_porder0;
        counter = 0)

    #-----------------------------------------------------

    dims_nodes_vh_θh =
        length.( gens_vh_θh )

    _,_, nodes_vh_θh_indx_in_Idx =
        create_size_offset_Idx(
            dims_nodes_vh_θh;
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
        
    sim_fun_system_ode_flat_agg_para  =
        vcat( f_gens_vh_θh,
             f_ωs_ωref0_vref0_porder0 )
    

    #-----------------------------------------------------

    sim_state_x0 =
        state[im_vars_Idx_in_state]

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
    
    sim_fun_system_kwd_flat_agg_para =
        (; vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         stateDiffCache_gens,
         sim_state_x0_gens,
         gens_nodes_collection ,
         gens_sim_fun_gen_para_Idxs,
         sim_fun_system_ode_flat_agg_para_Idxs,         
         ra_Xd_dash_Xq_dash_view,         
         sys_states_Idxs_and_mat_Idxs,         
         nodes_id_iq_pg_vh_idx_in_Idx,
         nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         nodes_vh_θh_indx_in_Idx )

    #-----------------------------------------------------
    
    sim_func! =
        system_flat_agg_ode_im_model_func!

    sim_fun_para =
        sim_fun_system_ode_flat_agg_para
    
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
    # simulate 
    #-----------------------------------------------------

    # system_sim_sol = DifferentialEquations.solve(
    #    sim_prob, algr  )

    system_sim_sol = DifferentialEquations.solve(
        sim_prob, algr, dt=dt  )

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
        sens_prob , algr, dt=dt  )
    
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


    #-----------------------------------------------------
    # case 2
    #-----------------------------------------------------
    
    poi_dyn_Idxs =
        (; vh_θh_Idx,
         ωs_ωref0_vref0_porder0_Idx )

    pois_dyn =
        (; f_gens_vh_θh,
         f_ωs_ωref0_vref0_porder0 )   

    #-----------------------------------------------------
    
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
    
    #-----------------------------------------------------

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

    #-----------------------------------------------------
    
    sim_timespan =
        sim_timespan
    #-----------------------------------------------------
    # case 3
    #-----------------------------------------------------


    #-----------------------------------------------------
    # case 4
    #-----------------------------------------------------
    
    #-----------------------------------------------------
    #----------------------------------------------------- 
    #-----------------------------------------------------
    #-----------------------------------------------------
    #-----------------------------------------------------
    # dynamics simulation
    #-----------------------------------------------------

    #-----------------------------------------------------
    # case 1
    #-----------------------------------------------------

    t_sim_func! = t_ode_im_model_func!
    
    #-----------------------------------------------------   
        
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
    
    #-----------------------------------------------------
    # simulate 
    #-----------------------------------------------------

    # sim_sol = DifferentialEquations.solve(
    #    sim_prob, algr  )

    t_sim_sol = DifferentialEquations.solve(
        t_sim_prob, algr, dt=dt  )

    #-----------------------------------------------------
    # plot
    #-----------------------------------------------------
    
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
        t_sens_prob , algr, dt=dt  )


    #-----------------------------------------------------
    # case 2
    #-----------------------------------------------------

    s_sim_func! = s_ode_one_im_model_func!
    
    #-----------------------------------------------------   
        
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
    #    sim_prob, algr  )

    s_sim_sol = DifferentialEquations.solve(
        s_sim_prob, algr, dt=dt  )

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
        s_sa_sens_prob , algr, dt=dt  )

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
######################################################
#-----------------------------------------------------    


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
         idxs) = create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )
        

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
             create_gens_nodes_aggregate_system_matrices_Idx_and_views(  gens_nodes_collection; only_gen = only_gen )
        
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
