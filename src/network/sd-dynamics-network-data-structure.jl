
########################################################
# ------------------------------------------------------
#  System Dynamis Network Data structures
# ------------------------------------------------------
########################################################


"""
```Julia
NetworkData(nodes, edges)
```

NetworkData  describes network data that is composed of nodes and edges.

It incorporate the following fields `nodes`, `edges`. `params_nodes` and `params_edges`. `params_nodes` and `params_edges` fields are vector of vector that contains numeric parameters.
"""
Base.@kwdef struct NetworkData
    nodes
    edges
    nodes_param       = get_components_params_value_vectors(nodes)
    edges_param       = get_components_params_value_vectors(edges)
    nodes_control_sig = get_components_control_sig(nodes)
    edges_control_sig = get_components_control_sig(nodes)
    nodes_output_sig  = get_components_output_sig(nodes)
    edges_output_sig  = get_components_output_sig(edges)

    # nodes_src_edges      = get_nodes_src_edges(edges)
    # nodes_dst_edges      = get_nodes_dst_edges(edges)
    # nodes_incident_edges = get_nodes_incident_edges(edges)

    # edges_orientation   = get_edges_orientation(edges)
    # edges_src_node      = first.(edges_orientation)
    # edges_dst_node      = last.(edges_orientation)

    # nodes_dyn_funcs
    # edges_dyn_funcs
end


function NetworkData(nodes::OrderedDict, edges::OrderedDict)
    bus_keys   = collect(keys(nodes))
    edge_array = collect(values(edges))
    edge_keys  = collect(keys(edges))
    
    # assert that keys are consistent
    @assert all([l.from ∈ bus_keys && l.to ∈ bus_keys for l in values(edges)]) "invalid node key given in edges"

    # Assert ordering of from/to according to index in bus_keys to assure
    @assert all([findfirst(bus_keys .== l.from) < findfirst(bus_keys .== l.to) for l in values(edges)]) "the pairs (from, to) need to be ordered according to the index of the corresponding keys in the node dict"
    
    nodes_param      = get_components_params_value_vectors(nodes)
    edges_param      = get_components_params_value_vectors(edges)
    
    nodes_control_sig = get_components_control_sig(nodes)
    edges_control_sig = get_components_control_sig(edges)
    
    nodes_output_sig  = get_components_output_sig(nodes)
    edges_output_sig  = get_components_output_sig(edges)

    NetworkData(nodes, edges, nodes_param, edges_param, nodes_control_sig, edges_control_sig, nodes_output_sig, edges_output_sig)
    
    # nodes_src_edges      = get_nodes_src_edges(edges)
    # nodes_dst_edges      = get_nodes_dst_edges(edges)
    # nodes_incident_edges = get_nodes_incident_edges(edges)

    # edges_orientation   = get_edges_orientation(edges)
    # edges_src_node      = first.(edges_orientation)
    # edges_dst_node      = last.(edges_orientation)

    # nodes_dyn_func     = get_components_dyn_func(nodes)
    # edges_dyn_func     = get_components_dyn_func(edges)

    # NetworkData(nodes, edges, nodes__param, edges_param,
    #             nodes_control_sig, edges_control_sig,
    #             nodes_output_sig, edges_output_sig,
    #             nodes_src_edges, nodes_dst_edges, nodes_incident_edges,
    #             edges_orientation, edges_src_node, edges_dst_node,
    #             nodes_dyn_func, edges_dyn_func  
    #             )
    

    
end


function NetworkData(nodes::Array, edges::Array)
    # assert that keys are consistent
    @assert all([l.from isa Int for l in edges]) "`from` should be of type `Int`"
    @assert all([l.to isa Int for   l in edges]) "`to` should be of type `Int`"
    @assert all([1 <= l.from <= length(nodes) for l in edges]) "numerical index needs to be between 1 and the number of nodes"
    @assert all([1 <= l.to <= length(nodes) for l in edges]) "numerical index needs to be between 1 and the number of nodes"

    # This could otherwise lead to problems for unsymmetric line types.
    @assert all([l.from < l.to for l in lines]) "the pairs (from, to) need to be ordered according to index value"
 
    nodes_param       = get_components_params_value_vectors(nodes)
    edges_param       = get_components_params_value_vectors(edges)
    nodes_control_sig = get_components_control_sig(nodes)
    edges_control_sig = get_components_control_sig(nodes)
    
    nodes_output_sig  = get_components_output_sig(nodes)
    edges_output_sig  = get_components_output_sig(edges)
    
    NetworkData(nodes, edges, nodes_param, edges_param, nodes_control_sig, edges_control_sig, nodes_output_sig, edges_output_sig)
    
    # nodes_src_edges      = get_nodes_src_edges(edges)
    # nodes_dst_edges      = get_nodes_dst_edges(edges)
    # nodes_incident_edges = get_nodes_incident_edges(edges)

    # edges_orientation   = get_edges_orientation(edges)
    # edges_src_node      = first.(edges_orientation)
    # edges_dst_node      = last.(edges_orientation)
    
    # nodes_dyn_func     = get_components_dyn_func(nodes)
    # edges_dyn_func     = get_components_dyn_func(edges)

    # NetworkData(nodes, edges, nodes_param, edges_param,
    #             nodes_control_sig, edges_control_sig,
    #             nodes_output_sig, edges_output_sig,
    #             nodes_src_edges, nodes_dst_edges, nodes_incident_edges,
    #             edges_orientation, edges_src_node, edges_dst_node,
    #             nodes_dyn_func, edges_dyn_func 
    #             )
 
end


