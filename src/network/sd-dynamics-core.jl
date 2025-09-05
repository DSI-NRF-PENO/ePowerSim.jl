# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


########################################################
# ------------------------------------------------------
#  Dynamics core-functions
# ------------------------------------------------------
########################################################

function get_dae_vars_vector(nd)
    return vcat(get_components_dae_var(nd.nodes), get_components_dae_var(nd.edges))
end


function get_mass_matrix(nd)
    
    vec_nodes_mass_matrix =
        get_components_mass_matrix(nd.nodes)
    
    vec_edges_mass_matrix =
        get_components_mass_matrix(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)
    
    # nodes_eq_size = sum(first.(size.(get_components_mass_matrix(nd.nodes))))
    # edges_eq_size = sum(first.(size.(get_components_mass_matrix(nd.edges))))

    system_eq_size =  nodes_dim_size + edges_dim_size
    
    mass_matrix = sparse(1.0I, system_eq_size, system_eq_size)
    
    for (i, node_mass_matrix) in enumerate(vec_nodes_mass_matrix)
        node_state_Idx = nodes_state_Idx[i]
        copyto!(@view(mass_matrix[node_state_Idx, node_state_Idx]), node_mass_matrix)
    end

    for (i, edge_mass_matrix) in enumerate(vec_edges_mass_matrix)
        edge_state_Idx = edges_state_Idx[i]
        copyto!(@view(mass_matrix[edge_state_Idx,edge_state_Idx]),edge_mass_matrix)
    end

    return mass_matrix
end


########################################################
# ------------------------------------------------------
#  Nodes and edges dynamics function call
# ------------------------------------------------------
########################################################


function edge_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, src_u, dst_u, para = p_arg

    fun_para = (src_u, dst_u, para)
    
    dyn_func(dx, x, fun_para, t)
    
    return nothing
end


function edges_dynamic_func!(dx, x, (edges_dyn_func,
                                         edges_dst_node,
                                         edges_src_node,
                                         nodes_u_Idx, 
                                         edges_param,
                                         edges_state_Idx,
                                         nd), t)  
 
    edges_dx_view    = [view(dx, edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    edges_x_view     = [view(x,  edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
    
    nodes_src_u_view = [nodes_u_view[Ind] for Ind in  edges_src_node]

    nodes_dst_u_view = [nodes_u_view[Ind] for Ind in  edges_dst_node]

    tup_edges_func    = Tuple(edges_dyn_func)
    tup_edges_dx      = Tuple(edges_dx_view)
    tup_edges_x       = Tuple(edges_x_view)
    
    tup_nodes_src_u   = Tuple(nodes_src_u_view)
    tup_nodes_dst_u   = Tuple(nodes_dst_u_view)    
    tup_edges_para    = Tuple(edges_param)    

    zip_edges_arg = zip(tup_edges_dx, tup_edges_x, tup_edges_func, tup_nodes_src_u, tup_nodes_dst_u, tup_edges_para)
    
    for (du, u, func, u_s, u_d, para) in zip_edges_arg
        
        p_arg = (func, u_s, u_d, para)
        
        edge_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    
end


function node_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, cb_sw, src_i, dst_i, f_t, para, plant = p_arg

    fun_para = (cb_sw, src_i, dst_i, f_t, para)
    
    dyn_func(dx, x, (fun_para, plant),  t)
    
    return nothing
end


function nodes_dynamic_func!(dx, x, ( nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind]) for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind]) for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_i_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_i_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_i   = Tuple(edges_src_i_view)
    tup_edges_dst_i   = Tuple(edges_dst_i_view)
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    
    tup_plants        = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_i, tup_edges_dst_i, tup_nodes_f_t, tup_nodes_para, tup_plants)
    
    for (du, u, func, cb_sw, edges_src_i, edges_dst_i, f_t, para, plant) in zip_nodes_arg
        p_arg = (func, cb_sw, edges_src_i, edges_dst_i, f_t, para, plant)
        node_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
        
end

# ------------------------------------------------------
# Aggregation
# ------------------------------------------------------

function nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    # nd, dyn_pf_para, nodes_and_edges_p = p

    nd, nodes_and_edges_p = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param 

 
    #-------------------------------------
    # Dynamic sectionm
    #-------------------------------------

    edges_dynamic_func!(
        dx, x, (edges_dyn_func,
                edges_dst_node,
                edges_src_node,
                nodes_u_Idx,
                edges_param,
                edges_state_Idx,
                nd), t)      
    
    nodes_dynamic_func!(
        dx, x, (nodes_dyn_func,
                nodes_cb_sw,
                nodes_src_edges,
                nodes_dst_edges,
                edges_ih_Idx,
                edges_ik_Idx,
                nodes_f_t,
                nodes_param,
                nodes_state_Idx,
                nd), t)

    return nothing
    
end


function system_dynamics(
    nodes_edges_dynamics!;
    mass_matrix = get_mass_matrix(nd),
    syms = get_network_vars_labels(nd)
    # syms = get_network_state_labels(nd)
                         )

    return ODEFunction{true}(
        nodes_edges_dynamics!;
        mass_matrix = mass_matrix,
        syms = syms)

end



# ------------------------------------------------------
# hybrid
# ------------------------------------------------------


function hybrid_pf_edge_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, src_u, dst_u, para = p_arg

    fun_para = (src_u, dst_u, para)
    
    dyn_func(dx, x, fun_para, t)
    
    return nothing
end


function hybrid_pf_edges_dynamic_func!(dx, x, (edges_dyn_func,
                                         edges_dst_node,
                                         edges_src_node,
                                         nodes_u_Idx, 
                                         edges_param,
                                         edges_state_Idx,
                                         nd), t)  
 
    edges_dx_view    = [view(dx, edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    edges_x_view     = [view(x,  edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
    
    nodes_src_u_view = [nodes_u_view[Ind] for Ind in  edges_src_node]

    nodes_dst_u_view = [nodes_u_view[Ind] for Ind in  edges_dst_node]

    tup_edges_func    = Tuple(edges_dyn_func)
    tup_edges_dx      = Tuple(edges_dx_view)
    tup_edges_x       = Tuple(edges_x_view)
    
    tup_nodes_src_u   = Tuple(nodes_src_u_view)
    tup_nodes_dst_u   = Tuple(nodes_dst_u_view)    
    tup_edges_para    = Tuple(edges_param)    

    zip_edges_arg = zip(tup_edges_dx, tup_edges_x, tup_edges_func, tup_nodes_src_u, tup_nodes_dst_u, tup_edges_para)
    
    for (du, u, func, u_s, u_d, para) in zip_edges_arg
        
        p_arg = (func, u_s, u_d, para)
        
        hybrid_pf_edge_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    
end



function hybrid_pf_node_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, cb_sw, src_i, dst_i, f_t, para, pf_U, Inet, node_idx_and_incident_edges_other_node_idx,  node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view,  plant, = p_arg

    dyn_node_pf_param = ( node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view )

    dyn_global_pf_param = ( pf_U, Inet )
    
    fun_para = (cb_sw, src_i, dst_i, f_t, para, dyn_global_pf_param, dyn_node_pf_param )
    
    dyn_func(dx, x, (fun_para, plant),  t)
    
    return nothing
end


function hybrid_pf_nodes_dynamic_func!(
    dx, x, ( nodes_dyn_func,
             nodes_cb_sw,
             nodes_src_edges,
             nodes_dst_edges,
             edges_ih_Idx,
             edges_ik_Idx,
             nodes_f_t,
             nodes_param,
             nodes_u_Idx,
             nodes_pf_U_view,
             Inet_view,
             nodes_node_idx_and_incident_edges_other_node_idx,
             nodes_incident_edges_Ybr_cal,
             nodes_incident_edges_orientation,
             nodes_state_Idx, nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind])
                     for Ind in
                         collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind])
                     for Ind in
                         collect(1:length(nodes_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind])
                     for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind])
                     for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_i_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_i_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_i   = Tuple(edges_src_i_view)
    tup_edges_dst_i   = Tuple(edges_dst_i_view)
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    tup_nodes_pf_U_view = Tuple(nodes_pf_U_view)
    tup_Inet_view       = Tuple(Inet_view)    
    tup_nodes_node_idx_and_incident_edges_other_node_idx = Tuple(nodes_node_idx_and_incident_edges_other_node_idx)
    tup_nodes_incident_edges_Ybr_cal = Tuple(nodes_incident_edges_Ybr_cal)
    tup_nodes_incident_edges_orientation   = Tuple(nodes_incident_edges_orientation)
    
    tup_plants   = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_i, tup_edges_dst_i, tup_nodes_f_t, tup_nodes_para, tup_nodes_pf_U_view, tup_Inet_view, tup_nodes_node_idx_and_incident_edges_other_node_idx, tup_nodes_incident_edges_Ybr_cal, tup_nodes_incident_edges_orientation, tup_plants)

    # I want to each node to have a view of all nodes voltages
    # this maybe useful for control purpose. I will supply
    # to each node nodes_u_view
    
    for (du, u, func, cb_sw, edges_src_i, edges_dst_i, f_t, para, pf_U,  Inet, node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, plant) in zip_nodes_arg
        p_arg = (func, cb_sw, edges_src_i, edges_dst_i, f_t, para, pf_U,  Inet, node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view, plant)
        hybrid_pf_node_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
        
end


function hybrid_pf_nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    # nd, nodes_and_edges_p,  stateDiffCache, state, counter_array, global_pf_options, nodes_pf_dyn_param  = p

    nd, nodes_and_edges_p, hybrid_pf_dyn_param = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param

    global_pf_param, stateDiffCache, state, counter_array, global_pf_options, nodes_pf_dyn_param = hybrid_pf_dyn_param
    
    nodes_node_idx_and_incident_edges_other_node_idx, nodes_incident_edges_Ybr_cal, nodes_incident_edges_orientation = nodes_pf_dyn_param

    #--------------------------------------------
    #--------------------------------------------    
    # Powerflow
    #--------------------------------------------
    #--------------------------------------------

    # pf_net_param, sd_pf_views, mismatch = global_pf_param

    # working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view = global_pf_views


    pf_net_param, sd_pf_views, mismatch = global_pf_param

    working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view = sd_pf_views
    

    #--------------------------------------------
    
    counter = counter_array[1]

    #--------------------------------------------

    stateDiffCache = get_tmp(stateDiffCache, x)

    if counter == 1
        
        stateDiffCache .= state
    end

    counter_array .+= 1

    #--------------------------------------------

    dyn_powerflow(nd, stateDiffCache, global_pf_param; global_pf_options... )

    
    #-------------------------------------
    # Dynamic sectionm
    #-------------------------------------
    
    hybrid_pf_edges_dynamic_func!(
        dx, x, (edges_dyn_func,
                edges_dst_node,
                edges_src_node,
                nodes_u_Idx,
                edges_param,
                edges_state_Idx,
                nd), t)  
    
    hybrid_pf_nodes_dynamic_func!(
        dx, x, (nodes_dyn_func,
                nodes_cb_sw,
                nodes_src_edges,
                nodes_dst_edges,
                edges_ih_Idx,
                edges_ik_Idx,
                nodes_f_t,
                nodes_param,
                nodes_u_Idx,
                nodes_pf_U_view,
                Inet_view,                
                nodes_node_idx_and_incident_edges_other_node_idx,
                nodes_incident_edges_Ybr_cal,
                nodes_incident_edges_orientation,
                nodes_state_Idx,
                nd), t )

    return nothing
    
end



function hybrid_pf_system_dynamics( hybrid_pf_nodes_edges_dynamics!;  mass_matrix = get_mass_matrix(nd), syms = get_network_vars_labels(nd) )

    return ODEFunction{true}( hybrid_pf_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms)

end


#---------------------------------------------------
#---------------------------------------------------
# Hybrid ENd
# ------------------------------------------------------
# ------------------------------------------------------

# ------------------------------------------------------
# network_current
# ------------------------------------------------------

function network_current_edge_components_dyn_func!(
    dx, x, p_arg, t)
    
    dyn_func, src_u, dst_u, para = p_arg

    fun_para = (src_u, dst_u, para)
    
    dyn_func(dx, x, fun_para, t)
    
    return nothing
end


function network_current_edges_dynamic_func!(
    dx, x, (edges_dyn_func,
            edges_dst_node,
            edges_src_node,
            nodes_u_Idx, 
            edges_param,
            edges_state_Idx,
            nd), t )  
 
    edges_dx_view    = [ view(dx, edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx)) ]

    edges_x_view     = [ view(x,  edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx)) ]

    nodes_u_view     = [ view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
    
    nodes_src_u_view = [ nodes_u_view[Ind] for Ind in  edges_src_node ]

    nodes_dst_u_view = [ nodes_u_view[Ind] for Ind in  edges_dst_node ]

    tup_edges_func    = Tuple(edges_dyn_func)
    tup_edges_dx      = Tuple(edges_dx_view)
    tup_edges_x       = Tuple(edges_x_view)
    
    tup_nodes_src_u   = Tuple(nodes_src_u_view)
    tup_nodes_dst_u   = Tuple(nodes_dst_u_view)    
    tup_edges_para    = Tuple(edges_param)    

    zip_edges_arg = zip(tup_edges_dx, tup_edges_x, tup_edges_func, tup_nodes_src_u, tup_nodes_dst_u, tup_edges_para)
    
    for (du, u, func, u_s, u_d, para) in zip_edges_arg
        
        p_arg = (func, u_s, u_d, para)
        
        network_current_edge_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    
end


function network_current_node_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, cb_sw, src_i, dst_i, f_t, para, plant = p_arg

    fun_para = (cb_sw, src_i, dst_i, f_t, para)
    
    dyn_func(dx, x, (fun_para, plant),  t)
    
    return nothing
end


function network_current_nodes_dynamic_func!(dx, x, ( nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind]) for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind]) for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_i_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_i_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_i   = Tuple(edges_src_i_view)
    tup_edges_dst_i   = Tuple(edges_dst_i_view)
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    
    tup_plants        = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_i, tup_edges_dst_i, tup_nodes_f_t, tup_nodes_para, tup_plants)
    
    for (du, u, func, cb_sw, edges_src_i, edges_dst_i, f_t, para, plant) in zip_nodes_arg
        p_arg = (func, cb_sw, edges_src_i, edges_dst_i, f_t, para, plant)
        network_current_node_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    

    
end



function network_current_nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    nd, nodes_and_edges_p = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param 

 
    #-------------------------------------
    # Dynamic sectionm
    #-------------------------------------
    
    network_current_edges_dynamic_func!(dx, x, (edges_dyn_func,
                                    edges_dst_node,
                                    edges_src_node,
                                    nodes_u_Idx,
                                    edges_param,
                                    edges_state_Idx,
                                    nd), t)  
    
    network_current_nodes_dynamic_func!(dx, x, (nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)

    return nothing
    
end


function network_current_system_dynamics( network_current_nodes_edges_dynamics!;  mass_matrix = get_mass_matrix(nd), syms = get_network_vars_labels(nd) )

    return ODEFunction{true}(network_current_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms)

end


# ------------------------------------------------------
# initial_pf_
# ------------------------------------------------------


function initial_pf_nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    nd, nodes_and_edges_p = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param 

 
    #-------------------------------------
    # Dynamic sectionm
    #-------------------------------------
    
    edges_dynamic_func!(dx, x, (edges_dyn_func,
                                    edges_dst_node,
                                    edges_src_node,
                                    nodes_u_Idx,
                                    edges_param,
                                    edges_state_Idx,
                                    nd), t)  
    
    nodes_dynamic_func!(dx, x, (nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)

    return nothing
    
end


function initial_pf_system_dynamics( initial_pf_nodes_edges_dynamics!;  mass_matrix = get_mass_matrix(nd), syms = get_network_vars_labels(nd) )

    return ODEFunction{true}( initial_pf_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms)

end


# ------------------------------------------------------
# f_t
# ------------------------------------------------------


function f_t_nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    nd, nodes_and_edges_p = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param 

 
    #-------------------------------------
    # Dynamic sectionm
    #-------------------------------------
    
    edges_dynamic_func!(dx, x, (edges_dyn_func,
                                    edges_dst_node,
                                    edges_src_node,
                                    nodes_u_Idx,
                                    edges_param,
                                    edges_state_Idx,
                                    nd), t)  
    
    nodes_dynamic_func!(dx, x, (nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)

    return nothing
    
end


function f_t_system_dynamics( f_t_nodes_edges_dynamics!;  mass_matrix = get_mass_matrix(nd), syms = get_network_vars_labels(nd) )

    return ODEFunction{true}(f_t_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms)

end


# ------------------------------------------------------
# external_pf
# ------------------------------------------------------


function external_pf_edge_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, src_u, dst_u, para = p_arg

    fun_para = (src_u, dst_u, para)
    
    dyn_func(dx, x, fun_para, t)
    
    return nothing
end


function external_pf_edges_dynamic_func!(dx, x, (edges_dyn_func,
                                         edges_dst_node,
                                         edges_src_node,
                                         nodes_u_Idx, 
                                         edges_param,
                                         edges_state_Idx,
                                         nd), t)  
 
    edges_dx_view    = [view(dx, edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    edges_x_view     = [view(x,  edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
    
    nodes_src_u_view = [nodes_u_view[Ind] for Ind in  edges_src_node]

    nodes_dst_u_view = [nodes_u_view[Ind] for Ind in  edges_dst_node]

    tup_edges_func    = Tuple(edges_dyn_func)
    tup_edges_dx      = Tuple(edges_dx_view)
    tup_edges_x       = Tuple(edges_x_view)
    
    tup_nodes_src_u   = Tuple(nodes_src_u_view)
    tup_nodes_dst_u   = Tuple(nodes_dst_u_view)    
    tup_edges_para    = Tuple(edges_param)    

    zip_edges_arg = zip(tup_edges_dx, tup_edges_x, tup_edges_func, tup_nodes_src_u, tup_nodes_dst_u, tup_edges_para)
    
    for (du, u, func, u_s, u_d, para) in zip_edges_arg
        
        p_arg = (func, u_s, u_d, para)
        
        external_pf_edge_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    
end

# node_inc_edges_Ybr_orient

# nodes_idx_and_inc_edges

function external_pf_node_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, cb_sw, src_i, dst_i, f_t, para, pf_U, Inet,  plant = p_arg

    fun_para = (cb_sw, src_i, dst_i, f_t, para, pf_U, Inet )
    
    dyn_func(dx, x, (fun_para, plant),  t)
    
    return nothing
end


function external_pf_nodes_dynamic_func!(dx, x, ( nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_pf_U_view,
                                      Inet_view,
                                      nodes_state_Idx, nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind]) for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind]) for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_i_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_i_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_i   = Tuple(edges_src_i_view)
    tup_edges_dst_i   = Tuple(edges_dst_i_view)
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    tup_nodes_pf_U_view = Tuple(nodes_pf_U_view)
    tup_Inet_view       = Tuple(Inet_view)
    
    tup_plants   = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_i, tup_edges_dst_i, tup_nodes_f_t, tup_nodes_para, tup_nodes_pf_U_view, tup_Inet_view, tup_plants)

    # I want to each node to have a view of all nodes voltages
    # this maybe useful for control purpose. I will supply
    # to each node nodes_u_view
    
    for (du, u, func, cb_sw, edges_src_i, edges_dst_i, f_t, para, pf_U, Inet, plant) in zip_nodes_arg
        p_arg = (func, cb_sw, edges_src_i, edges_dst_i, f_t, para, pf_U, Inet,  plant)
        external_pf_node_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    

    
end


function external_pf_nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    nd, nodes_and_edges_p, dyn_pf_param    = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param 

    #-------------------------------------


    nodes_pf_U_view, Inet_view, Iinj_view = dyn_pf_param 
 
    #-------------------------------------
    # Powerflow
    #-------------------------------------
        
    #-------------------------------------
    # Dynamic sectionm
    #-------------------------------------
    
    external_pf_edges_dynamic_func!(dx, x, (edges_dyn_func,
                                    edges_dst_node,
                                    edges_src_node,
                                    nodes_u_Idx,
                                    edges_param,
                                    edges_state_Idx,
                                    nd), t)  
    
    external_pf_nodes_dynamic_func!(dx, x, (nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_pf_U_view,
                                      Inet_view,
                                    nodes_state_Idx, nd), t)

    return nothing
    
end


function external_pf_system_dynamics(external_pf_nodes_edges_dynamics!; mass_matrix = get_mass_matrix(nd),  syms = get_network_vars_labels(nd) )

    return ODEFunction{true}(external_pf_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms)

end


# ------------------------------------------------------
# global_pf
# ------------------------------------------------------


function global_pf_edge_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, src_u, dst_u, para = p_arg

    fun_para = (src_u, dst_u, para)
    
    dyn_func(dx, x, fun_para, t)
    
    return nothing
end


function global_pf_edges_dynamic_func!(dx, x, (edges_dyn_func,
                                         edges_dst_node,
                                         edges_src_node,
                                         nodes_u_Idx, 
                                         edges_param,
                                         edges_state_Idx,
                                         nd), t)  
 
    edges_dx_view    = [view(dx, edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    edges_x_view     = [view(x,  edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
    
    nodes_src_u_view = [nodes_u_view[Ind] for Ind in  edges_src_node]

    nodes_dst_u_view = [nodes_u_view[Ind] for Ind in  edges_dst_node]

    tup_edges_func    = Tuple(edges_dyn_func)
    tup_edges_dx      = Tuple(edges_dx_view)
    tup_edges_x       = Tuple(edges_x_view)
    
    tup_nodes_src_u   = Tuple(nodes_src_u_view)
    tup_nodes_dst_u   = Tuple(nodes_dst_u_view)    
    tup_edges_para    = Tuple(edges_param)    

    zip_edges_arg = zip(tup_edges_dx, tup_edges_x, tup_edges_func, tup_nodes_src_u, tup_nodes_dst_u, tup_edges_para)
    
    for (du, u, func, u_s, u_d, para) in zip_edges_arg
        
        p_arg = (func, u_s, u_d, para)
        
        global_pf_edge_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    
end


#---------------------------------------------------
# New global pf
#---------------------------------------------------


function global_pf_node_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, cb_sw, src_i, dst_i, f_t, para, pf_U, Inet,  plant = p_arg

    global_pf_param = ( pf_U, Inet )
    
    fun_para = (cb_sw, src_i, dst_i, f_t, para, global_pf_param )
    
    dyn_func(dx, x, (fun_para, plant),  t)
    
    return nothing
end


function global_pf_nodes_dynamic_func!(dx, x, ( nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_u_Idx,
                                      nodes_pf_U_view,
                                      Inet_view,
                                      nodes_state_Idx, nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind]) for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind]) for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_i_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_i_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_i   = Tuple(edges_src_i_view)
    tup_edges_dst_i   = Tuple(edges_dst_i_view)
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    tup_nodes_pf_U_view = Tuple(nodes_pf_U_view)
    tup_Inet_view       = Tuple(Inet_view)
    
    tup_plants   = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_i, tup_edges_dst_i, tup_nodes_f_t, tup_nodes_para, tup_nodes_pf_U_view, tup_Inet_view, tup_plants)

    # I want to each node to have a view of all nodes voltages
    # this maybe useful for control purpose. I will supply
    # to each node nodes_u_view
    
    for (du, u, func, cb_sw, edges_src_i, edges_dst_i, f_t, para, pf_U, Inet, plant) in zip_nodes_arg
        p_arg = (func, cb_sw, edges_src_i, edges_dst_i, f_t, para, pf_U, Inet,  plant)
        global_pf_node_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    
end



function global_pf_nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    nd, nodes_and_edges_p, global_pf_param, stateDiffCache, state, counter_array, global_pf_options  = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param 

    #--------------------------------------------
    #--------------------------------------------    
    # Powerflow
    #--------------------------------------------
    #--------------------------------------------

    # pf_net_param, sd_pf_views, mismatch = global_pf_param

    # working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view = global_pf_views


    pf_net_param, sd_pf_views, mismatch = global_pf_param

    working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view = sd_pf_views
    

    #--------------------------------------------
    
    counter = counter_array[1]

    #--------------------------------------------

    stateDiffCache = get_tmp(stateDiffCache, x)

    if counter == 1
        
        stateDiffCache .= state
    end

    counter_array .+= 1

    #--------------------------------------------

    dyn_powerflow(nd,
                  stateDiffCache,
                  global_pf_param;
                  global_pf_options... )

    #--------------------------------------------    
    # Dynamic sectionm
    #--------------------------------------------
    
    global_pf_edges_dynamic_func!(
        dx, x, (edges_dyn_func,
        edges_dst_node,
        edges_src_node,
        nodes_u_Idx,
        edges_param,
        edges_state_Idx,
        nd), t)  
    
    test_global_pf_nodes_dynamic_func!(
        dx, x, (nodes_dyn_func,
        nodes_cb_sw,
        nodes_src_edges,
        nodes_dst_edges,
        edges_ih_Idx,
        edges_ik_Idx,
        nodes_f_t,
        nodes_param,
        nodes_u_Idx,
        nodes_pf_U_view,
        Inet_view,
        nodes_state_Idx, nd), t)

    return nothing
    
end

#---------------------------------------------------
#---------------------------------------------------


function global_pf_system_dynamics(global_pf_nodes_edges_dynamics!; mass_matrix = get_mass_matrix(nd),  syms = get_network_vars_labels(nd) )

    return ODEFunction{true}(global_pf_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms)

end


# ------------------------------------------------------
# node_pf
# ------------------------------------------------------


function node_pf_edge_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, src_u, dst_u, para = p_arg

    fun_para = (src_u, dst_u, para)
    
    dyn_func(dx, x, fun_para, t)
    
    return nothing
end


function node_pf_edges_dynamic_func!(dx, x, (edges_dyn_func,
                                         edges_dst_node,
                                         edges_src_node,
                                         nodes_u_Idx, 
                                         edges_param,
                                         edges_state_Idx,
                                         nd), t)  
 
    edges_dx_view    = [view(dx, edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    edges_x_view     = [view(x,  edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
    
    nodes_src_u_view = [nodes_u_view[Ind] for Ind in  edges_src_node]

    nodes_dst_u_view = [nodes_u_view[Ind] for Ind in  edges_dst_node]

    tup_edges_func    = Tuple(edges_dyn_func)
    tup_edges_dx      = Tuple(edges_dx_view)
    tup_edges_x       = Tuple(edges_x_view)
    
    tup_nodes_src_u   = Tuple(nodes_src_u_view)
    tup_nodes_dst_u   = Tuple(nodes_dst_u_view)    
    tup_edges_para    = Tuple(edges_param)    

    zip_edges_arg = zip(tup_edges_dx, tup_edges_x, tup_edges_func, tup_nodes_src_u, tup_nodes_dst_u, tup_edges_para)
    
    for (du, u, func, u_s, u_d, para) in zip_edges_arg
        
        p_arg = (func, u_s, u_d, para)
        
        node_pf_edge_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    
end


#---------------------------------------------------
# New node pf
#---------------------------------------------------


function node_pf_node_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, cb_sw, src_i, dst_i, f_t, para, node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view, plant = p_arg

    node_pf_param = ( node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view )

    fun_para = (cb_sw, src_i, dst_i, f_t, para, node_pf_param )
    
    dyn_func(dx, x, (fun_para, plant),  t)
    
    return nothing
end


function node_pf_nodes_dynamic_func!(
    dx, x, ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges,
             nodes_dst_edges, edges_ih_Idx, edges_ik_Idx,
             nodes_f_t, nodes_param, nodes_u_Idx,
             nodes_node_idx_and_incident_edges_other_node_idx,
             nodes_incident_edges_Ybr_cal, nodes_incident_edges_orientation,
             nodes_state_Idx, nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind]) for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind]) for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_i_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_i_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_i   = Tuple(edges_src_i_view)
    tup_edges_dst_i   = Tuple(edges_dst_i_view)
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    tup_nodes_node_idx_and_incident_edges_other_node_idx = Tuple(nodes_node_idx_and_incident_edges_other_node_idx)
    tup_nodes_incident_edges_Ybr_cal = Tuple(nodes_incident_edges_Ybr_cal)
    tup_nodes_incident_edges_orientation   = Tuple(nodes_incident_edges_orientation)
    
    tup_plants   = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_i, tup_edges_dst_i, tup_nodes_f_t, tup_nodes_para, tup_nodes_node_idx_and_incident_edges_other_node_idx, tup_nodes_incident_edges_Ybr_cal, tup_nodes_incident_edges_orientation, tup_plants)

    # I want to each node to have a view of all nodes voltages
    # this maybe useful for control purpose. I will supply
    # to each node nodes_u_view
    
    for (du, u, func, cb_sw, edges_src_i, edges_dst_i, f_t, para, node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, plant) in zip_nodes_arg
        p_arg = (func, cb_sw, edges_src_i, edges_dst_i, f_t, para, node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view, plant)
        node_pf_node_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    

    
end


function node_pf_nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    nd, nodes_and_edges_p,  nodes_pf_dyn_param  = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param


    nodes_node_idx_and_incident_edges_other_node_idx, nodes_incident_edges_Ybr_cal, nodes_incident_edges_orientation = nodes_pf_dyn_param

    #-------------------------------------
    # Dynamic sectionm
    #-------------------------------------
    
    node_pf_edges_dynamic_func!(
        dx, x, (edges_dyn_func, edges_dst_node,
                edges_src_node, nodes_u_Idx,
                edges_param, edges_state_Idx,
                nd), t)  
    
    node_pf_nodes_dynamic_func!(
        dx, x, (nodes_dyn_func, nodes_cb_sw,
                nodes_src_edges, nodes_dst_edges,
                edges_ih_Idx, edges_ik_Idx, nodes_f_t,
                nodes_param, nodes_u_Idx,
                nodes_node_idx_and_incident_edges_other_node_idx,
                nodes_incident_edges_Ybr_cal,
                nodes_incident_edges_orientation,
                nodes_state_Idx,
                nd), t )

    return nothing
    
end


#---------------------------------------------------
#---------------------------------------------------

function node_pf_system_dynamics( node_pf_nodes_edges_dynamics!; mass_matrix = get_mass_matrix(nd),  syms = get_network_vars_labels(nd) )

    return ODEFunction{true}( node_pf_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms )

end



#-------------------------------------------------------
# Test
#-------------------------------------------------------


function test_global_pf_nodes_dynamic_func!(dx, x, ( nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_u_Idx,
                                      nodes_pf_U_view,
                                      Inet_view,
                                      nodes_state_Idx, nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind]) for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind]) for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_i_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_i_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_i   = Tuple(edges_src_i_view)
    tup_edges_dst_i   = Tuple(edges_dst_i_view)
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    tup_nodes_pf_U_view = Tuple(nodes_pf_U_view)
    tup_Inet_view       = Tuple(Inet_view)
    
    tup_plants   = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_i, tup_edges_dst_i, tup_nodes_f_t, tup_nodes_para, tup_nodes_pf_U_view, tup_Inet_view, tup_plants)

    # I want to each node to have a view of all nodes voltages
    # this maybe useful for control purpose. I will supply
    # to each node nodes_u_view
    
    for (du, u, func, cb_sw, edges_src_i, edges_dst_i, f_t, para, pf_U, Inet, plant) in zip_nodes_arg
        p_arg = (func, cb_sw, edges_src_i, edges_dst_i, f_t, para, pf_U, Inet,  plant)
        global_pf_node_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    
end


function test_global_pf_nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    nd, nodes_and_edges_p, global_pf_param, stateDiffCache, state, counter_array, global_pf_options  = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param 

    #--------------------------------------------
    #--------------------------------------------    
    # Powerflow
    #--------------------------------------------
    #--------------------------------------------

    pf_net_param, sd_pf_views, mismatch = global_pf_param

    working_vh_θh_view, nodes_pf_U_view, Inet_view, Iinj_view, δ_ω_ed_dash_eq_dash_view = sd_pf_views

    #--------------------------------------------
    
    counter = counter_array[1]

    #--------------------------------------------

    stateDiffCache = get_tmp(stateDiffCache, x)

    if counter == 1
        
        stateDiffCache .= state
    end

    counter_array .+= 1

    #--------------------------------------------

    dyn_powerflow(nd, stateDiffCache, global_pf_param ; global_pf_options... )

    #--------------------------------------------    
    # Dynamic sectionm
    #--------------------------------------------
    
    global_pf_edges_dynamic_func!(
        dx, x, (edges_dyn_func,
        edges_dst_node,
        edges_src_node,
        nodes_u_Idx,
        edges_param,
        edges_state_Idx,
        nd), t)  
    
    test_global_pf_nodes_dynamic_func!(
        dx, x, (nodes_dyn_func,
        nodes_cb_sw,
        nodes_src_edges,
        nodes_dst_edges,
        edges_ih_Idx,
        edges_ik_Idx,
        nodes_f_t,
        nodes_param,
        nodes_u_Idx,
        nodes_pf_U_view,
        Inet_view,
        nodes_state_Idx, nd), t)

    return nothing
    
end

#---------------------------------------------------


function test_global_pf_system_dynamics( test_global_pf_nodes_edges_dynamics!; mass_matrix = get_mass_matrix(nd),  syms = get_network_vars_labels(nd) )

    return ODEFunction{true}( test_global_pf_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms)

end


# ------------------------------------------------------
# ------------------------------------------------------

function steady_state_nodes_edges_dynamics!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    nd, nodes_and_edges_p = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param 

 
    #-------------------------------------
    # Dynamic sectionm
    #-------------------------------------
    
    edges_dynamic_func!(dx, x, (edges_dyn_func,
                                    edges_dst_node,
                                    edges_src_node,
                                    nodes_u_Idx,
                                    edges_param,
                                    edges_state_Idx,
                                    nd), t)  
    
    nodes_dynamic_func!(dx, x, (nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)

    return nothing
    
end


function steady_state_system_dynamics( steady_state_nodes_edges_dynamics!; mass_matrix = get_mass_matrix(nd), syms = get_network_vars_labels(nd) )

    return ODEFunction{true}( steady_state_nodes_edges_dynamics!;  mass_matrix = mass_matrix, syms = syms)

end


# ------------------------------------------------------
# pf : power flow
# ------------------------------------------------------


function pf_system_dynamics( steady_state_nodes_edges_dynamics!; mass_matrix = get_mass_matrix(nd), syms = get_network_vars_labels(nd) )

    return NonlinearFunction{true}( steady_state_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms)

end


# ------------------------------------------------------
# ss : steady state
# ------------------------------------------------------


function ss_system_dynamics( steady_state_nodes_edges_dynamics!; mass_matrix = get_mass_matrix(nd), syms = get_network_vars_labels(nd) )

    return ODEFunction{true}( steady_state_nodes_edges_dynamics!; mass_matrix = mass_matrix, syms = syms)

end


# ------------------------------------------------------
# ------------------------------------------------------

function nodes_edges_dynamics_with_pf!(dx, x, p, t)
    
    # I passed param seperately so that
    # parameters can be controlled in ODE integrator

    nd,  nodes_and_edges_p, nodes_f_t_func, perm_matrix, dyn_pf_para = p
    
    sd_nodes_param, sd_edges_param = nodes_and_edges_p

     nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx, ur_ui_Idx_in_state  = sd_nodes_param

    edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx = sd_edges_param 


    #-------------------------------------
    # Power flow section
    #-------------------------------------

    comp_U = get_components_u_in_state( x, ur_ui_Idx_in_state)


    x_ΘV = u_to_ΘV( perm_matrix' * comp_U )

    tuple_pf_result = x_get_dynamic_powerflow_sol( x_ΘV, dyn_pf_para... )
    
    # tuple_pf_result = get_dynamic_powerflow_sol( x_ΘV, dyn_pf_para... )

    # U_Pf = perm_matrix' * comp_U
    
    # tuple_pf_result = get_dynamic_powerflow_sol( U_Pf, dyn_pf_para... )
    
    nodes_pf_f_t = get_nodes_pf_f_t( nodes_f_t_func, tuple_pf_result ) 

    # use Accessors.jl to change the first element of nodes_f_t with nodes_pf_f_t

    # https://juliaobjects.github.io/Accessors.jl/dev/docstrings/#Accessors.IndexLens-Tuple{Vararg{Integer}}

    IL = IndexLens(1)
    
    nodes_f_t_m = map( (arg) -> set(arg[1], IL, arg[2]), [(obj, value) for (obj, value) in zip(nodes_f_t, nodes_pf_f_t)] )
    
    #-------------------------------------
    # Dynamic sectionm
    #-------------------------------------

    edges_dynamic_func!(dx, x, (edges_dyn_func,
                                    edges_dst_node,
                                    edges_src_node,
                                    nodes_u_Idx,
                                    edges_param,
                                    edges_state_Idx,
                                    nd), t)  
    
    nodes_dynamic_func!(dx, x, (nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t_m,
                                      # nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)
    return nothing
    
end


function system_dynamics_with_pf(nodes_edges_dynamics_with_pf!; mass_matrix = get_mass_matrix(nd), syms = get_network_vars_labels(nd) )

    return ODEFunction{true}(nodes_edges_dynamics_with_pf!;
                             mass_matrix = mass_matrix,
                             syms = syms)

end


# ------------------------------------------------------
# Diagnostics
# ------------------------------------------------------


function diagnosis_node_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, cb_sw, src_ih, src_ik, dst_ih, dst_ik, f_t, para, plant = p_arg

    fun_para = (cb_sw, (src_ih, src_ik), (dst_ih, dst_ik), f_t, para)
    
    dyn_func(dx, x, (fun_para, plant),  t)
    
    return nothing
end



function diagnosis_nodes_dynamic_func!(dx, x, ( nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind]) for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind]) for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_ih_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]
    edges_src_ik_view = [srcs==[] ? [] : [edges_ik_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_ih_view = [dsts==[] ? [] : [edges_ih_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]
    edges_dst_ik_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_ih  = Tuple(edges_src_ih_view)
    tup_edges_src_ik  = Tuple(edges_src_ik_view)    
    tup_edges_dst_ih  = Tuple(edges_dst_ih_view)
    tup_edges_dst_ik  = Tuple(edges_dst_ik_view)    
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    
    tup_plants        = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_ih, tup_edges_src_ik, tup_edges_dst_ih, tup_edges_dst_ik, tup_nodes_f_t, tup_nodes_para, tup_plants)
    
    for (du, u, func, cb_sw, edges_src_ih, edges_src_ik, edges_dst_ih, edges_dst_ik, f_t, para, plant) in zip_nodes_arg
        p_arg = (func, cb_sw, edges_src_ih, edges_src_ik, edges_dst_ih, edges_dst_ik, f_t, para, plant)
        network_current_node_components_dyn_func!(du, u, p_arg, t)
    end
    
    return nothing
    

    
end



function diagnostics_edges_dynamic_func!(dx, x, (edges_dyn_func,
                                         edges_dst_node,
                                         edges_src_node,
                                         nodes_u_Idx, 
                                         edges_param,
                                         edges_state_Idx,
                                         nd), t)  
 
    edges_dx_view    = [view(dx, edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    edges_x_view     = [view(x,  edges_state_Idx[Ind]) for Ind in collect(1:length(edges_state_Idx))]

    nodes_u_view     = [view(x, nodes_u_Idx[Ind]) for Ind in collect(1:length(nodes_u_Idx))]
    
    nodes_src_u_view = [nodes_u_view[Ind] for Ind in  edges_src_node]

    nodes_dst_u_view = [nodes_u_view[Ind] for Ind in  edges_dst_node]

    tup_edges_func    = Tuple(edges_dyn_func)
    tup_edges_dx      = Tuple(edges_dx_view)
    tup_edges_x       = Tuple(edges_x_view)
    
    tup_nodes_src_u   = Tuple(nodes_src_u_view)
    tup_nodes_dst_u   = Tuple(nodes_dst_u_view)    
    tup_edges_para    = Tuple(edges_param)    

    zip_edges_arg = zip(tup_edges_dx, tup_edges_x, tup_edges_func, tup_nodes_src_u, tup_nodes_dst_u, tup_edges_para)

    param_edge_components_dyn_func = []
    
    for (du, u, func, u_s, u_d, para) in zip_edges_arg

        p_arg = (func, u_s, u_d, para)
        
        push!(param_edge_components_dyn_func, p_arg)
        
    end
    
    return param_edge_components_dyn_func
    
end


function diagnostics_nodes_dynamic_func!(dx, x, ( nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind]) for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind]) for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_i_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_i_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_i   = Tuple(edges_src_i_view)
    tup_edges_dst_i   = Tuple(edges_dst_i_view)
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    
    tup_plants        = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_i, tup_edges_dst_i, tup_nodes_f_t, tup_nodes_para, tup_plants)

    param_node_components_dyn_func = []
    
    for (du, u, func, cb_sw, edges_src_i, edges_dst_i, f_t, para, plant) in zip_nodes_arg
        
        p_arg = (func, cb_sw, edges_src_i, edges_dst_i, f_t, para, plant)

        push!(param_node_components_dyn_func, p_arg)
        
    end
    
    return param_node_components_dyn_func
    

    
end




function diagnosis_network_current_node_components_dyn_func!(dx, x, p_arg, t)
    
    dyn_func, cb_sw, src_ih, src_ik, dst_ih, dst_ik, f_t, para, plant = p_arg

    fun_para = (cb_sw, (src_ih, src_ik), (dst_ih, dst_ik), f_t, para)
    
    dyn_func(dx, x, (fun_para, plant),  t)
    
    return nothing
end



function diagnosis_network_current_nodes_dynamic_func!(dx, x, ( nodes_dyn_func,
                                      nodes_cb_sw,
                                      nodes_src_edges,
                                      nodes_dst_edges,
                                      edges_ih_Idx,
                                      edges_ik_Idx,
                                      nodes_f_t,
                                      nodes_param,
                                      nodes_state_Idx,
                                      nd), t)
    
    nodes_dx_view = [view(dx, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]
   
    edges_ih_view = [view(x, edges_ih_Idx[Ind]) for Ind in collect(1:length(edges_ih_Idx))]
    
    edges_ik_view = [view(x, edges_ik_Idx[Ind]) for Ind in collect(1:length(edges_ik_Idx))]
    
    # a src edge or dst edge may be null, there is a need to take care of that

    edges_src_ih_view = [srcs==[] ? [] : [edges_ih_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_src_ik_view = [srcs==[] ? [] : [edges_ik_view[Ind] for Ind in srcs] for srcs in nodes_src_edges]

    edges_dst_ih_view = [dsts==[] ? [] : [edges_ih_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    edges_dst_ik_view = [dsts==[] ? [] : [edges_ik_view[Ind] for Ind in dsts] for dsts in nodes_dst_edges]

    tup_nodes_func    = Tuple(nodes_dyn_func)
    tup_nodes_dx      = Tuple(nodes_dx_view)
    tup_nodes_x       = Tuple(nodes_x_view)
    
    tup_nodes_cb_sw   = Tuple(nodes_cb_sw)
    tup_edges_src_ih  = Tuple(edges_src_ih_view)
    tup_edges_src_ik  = Tuple(edges_src_ik_view)    
    tup_edges_dst_ih  = Tuple(edges_dst_ih_view)
    tup_edges_dst_ik  = Tuple(edges_dst_ik_view)    
    tup_nodes_f_t     = Tuple(nodes_f_t)
    tup_nodes_para    = Tuple(nodes_param)
    
    tup_plants        = Tuple(collect(values(nd.nodes)))
     
    zip_nodes_arg = zip(tup_nodes_dx, tup_nodes_x, tup_nodes_func, tup_nodes_cb_sw, tup_edges_src_ih, tup_edges_src_ik, tup_edges_dst_ih, tup_edges_dst_ik, tup_nodes_f_t, tup_nodes_para, tup_plants)

    param_node_components_dyn_func = []
    
    for (du, u, func, cb_sw, edges_src_ih, edges_src_ik, edges_dst_ih, edges_dst_ik, f_t, para, plant) in zip_nodes_arg
        p_arg = ( func, cb_sw, edges_src_ih, edges_src_ik, edges_dst_ih, edges_dst_ik, f_t, para, plant )
        

        push!(param_node_components_dyn_func, p_arg)
        
    end
    
    return param_node_components_dyn_func
    

    
end

