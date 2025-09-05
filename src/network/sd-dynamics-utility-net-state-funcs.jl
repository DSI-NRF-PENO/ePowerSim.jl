# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


######################################################
#-----------------------------------------------------
#  utility functions with plant, netd, etc as arguments 
#-----------------------------------------------------
######################################################


function get_net_to_industrial_model_indices_dict(
    netd;
    no_control_device = false  )

    if no_control_device == false
        
        pure_states_indices_in_system =
            [ get_gens_pure_states_indices_in_system( netd )...; ]

        nodes_ur_ui_indices_in_system =
            [ get_nodes_ur_ui_indices_in_system(netd )...; ]

        dim_industrial_model_state =
            length.( [ pure_states_indices_in_system,
                       nodes_ur_ui_indices_in_system] )

        _, _, industrial_Idx =
            create_size_offset_Idx(
                dim_industrial_model_state;
                counter = 0)

        industrial_model_pure_states_Idx =
            industrial_Idx[1] 

        industrial_model_ur_ui_Idx =
            industrial_Idx[2]

        #------------------------------------------------
        # indices transformation
        #-----------------------------------------------

        indices_in_system =
            [ pure_states_indices_in_system...;
              nodes_ur_ui_indices_in_system... ]

        indices_in_industrial =
            [collect(industrial_model_pure_states_Idx)...;
             collect(industrial_model_ur_ui_Idx)... ]

        return  OrderedDict{Int64, Int64}(
            sys => ind for (sys, ind) in
                zip( indices_in_system,
                     indices_in_industrial  ) )
        
    else
        
        pure_states_indices_in_system =
            [ get_gens_pure_states_indices_in_system(
                netd; no_control_device = true )...; ]

        nodes_ur_ui_indices_in_system =
            [ get_nodes_ur_ui_indices_in_system(netd )...; ]

        dim_industrial_model_state =
            length.(
                [ pure_states_indices_in_system,
                  nodes_ur_ui_indices_in_system] )

        _, _, industrial_Idx =
            create_size_offset_Idx(
                dim_industrial_model_state;
                counter = 0)

        industrial_model_pure_states_Idx =
            industrial_Idx[1] 

        industrial_model_ur_ui_Idx  =
            industrial_Idx[2]

        #----------------------------------------------------
        # indices transformation
        #----------------------------------------------------

        indices_in_system =
            [ pure_states_indices_in_system...;
              nodes_ur_ui_indices_in_system... ]

        indices_in_industrial =
            [collect(industrial_model_pure_states_Idx)...;
             collect(industrial_model_ur_ui_Idx)... ]


        return  OrderedDict{Int64, Int64}(
            sys => ind for (sys, ind) in
                zip( indices_in_system,
                      indices_in_industrial  ) )
        
    end
    
    
end


function get_industrial_model_indices_and_conversion_dict(
    netd; no_control_device = false )


    if no_control_device == false
        pure_states_Idx_in_system =
            [ get_gens_pure_states_indices_in_system(
                netd )...; ]

        ur_ui_Idx_in_system =
            [ get_nodes_ur_ui_indices_in_system(
                netd )...; ]

        dim_industrial_model_state =
            length.( [ pure_states_Idx_in_system,
                       ur_ui_Idx_in_system] )

        _, _, industrial_Idx =
            create_size_offset_Idx(
                dim_industrial_model_state;
                counter = 0)

        industrial_model_pure_states_Idx =
            industrial_Idx[1] 

        industrial_model_ur_ui_Idx =
            industrial_Idx[2]

        #----------------------------------------------------
        # indices transformation
        #----------------------------------------------------

        pure_states_and_ur_ui_Idx_in_system =
            [ pure_states_Idx_in_system...;
              ur_ui_Idx_in_system... ]

        industrial_model_state_Idx =
            [collect(industrial_model_pure_states_Idx)...;
             collect(industrial_model_ur_ui_Idx)... ]


        net_to_industrial_idx_conversion_dict =
            OrderedDict{Int64, Int64}(
                sys => ind for (sys, ind) in
                    zip( pure_states_and_ur_ui_Idx_in_system,
                          industrial_model_state_Idx  ) )

        return  (; pure_states_Idx_in_system,
                 ur_ui_Idx_in_system,
                 industrial_Idx,
                 industrial_model_pure_states_Idx,
                 industrial_model_ur_ui_Idx,
                 pure_states_and_ur_ui_Idx_in_system,
                 industrial_model_state_Idx,
                 net_to_industrial_idx_conversion_dict )

    else
        pure_states_Idx_in_system =
            [ get_gens_pure_states_indices_in_system(
                netd; no_control_device = true )...; ]

        ur_ui_Idx_in_system =
            [ get_nodes_ur_ui_indices_in_system(
                netd )...; ]

        dim_industrial_model_state =
            length.( [ pure_states_Idx_in_system,
                       ur_ui_Idx_in_system] )

        _, _, industrial_Idx =
            create_size_offset_Idx(
                dim_industrial_model_state;
                counter = 0)

        industrial_model_pure_states_Idx =
            industrial_Idx[1] 

        industrial_model_ur_ui_Idx =
            industrial_Idx[2]

        #----------------------------------------------------
        # indices transformation
        #----------------------------------------------------

        pure_states_and_ur_ui_Idx_in_system =
            [ pure_states_Idx_in_system...;
              ur_ui_Idx_in_system... ]

        industrial_model_state_Idx =
            [collect(industrial_model_pure_states_Idx)...;
             collect(industrial_model_ur_ui_Idx)... ]


        net_to_industrial_idx_conversion_dict =
            OrderedDict{Int64, Int64}(
                sys => ind for (sys, ind) in
                    zip( pure_states_and_ur_ui_Idx_in_system,
                          industrial_model_state_Idx  ) )

        return  (; pure_states_Idx_in_system,
                 ur_ui_Idx_in_system,
                 industrial_Idx,
                 industrial_model_pure_states_Idx,
                 industrial_model_ur_ui_Idx,
                 pure_states_and_ur_ui_Idx_in_system,
                 industrial_model_state_Idx,
                 net_to_industrial_idx_conversion_dict )

    end
    
    # (; pure_states_indices_in_system, nodes_ur_ui_indices_in_system, industrial_Idx, industrial_model_pure_states_Idx, industrial_model_ur_ui_Idx, pure_indices_in_system, pure_indices_in_industrial, net_to_industrial_idx_conversion_dict )
    
end

# ------------------------------------------------------
# ------------------------------------------------------


function get_im_indices_and_conversion_dict( netd  )

    each_gens_im_vars_dim =
        length.(get_gens_im_vars_indices_in_system(netd ))

    
    _, _, each_gens_im_vars_Idx_in_state =
        create_size_offset_Idx( each_gens_im_vars_dim )

    
    im_vars_indices_in_system =
        [get_gens_im_vars_indices_in_system(netd )...;]
    
    pure_states_Idx_in_system = [
        get_gens_pure_states_indices_in_system(
            netd )...; ]

    im_algebraic_vars_Idx_in_system = [
        get_gens_im_algebraic_vars_indices_in_system(
            netd )...; ]
    
    ur_ui_Idx_in_system = [
        get_nodes_ur_ui_indices_in_system(netd )...; ]

    dim_im_vars_and_ur_ui_vars = length.(
        [ im_vars_indices_in_system,
          ur_ui_Idx_in_system] )

    _, _, im_Idx =  create_size_offset_Idx(
        dim_im_vars_and_ur_ui_vars )

    im_vars_Idx_in_state  = im_Idx[1] 
    
    nodes_ur_ui_Idx_in_state = im_Idx[2]

    #----------------------------------------------------
    # indices transformation
    #----------------------------------------------------

    im_vars_and_ur_ui_Idx_in_system = [
        im_vars_indices_in_system...;
        ur_ui_Idx_in_system... ]
    
    im_state_Idx = [
        collect( im_vars_Idx_in_state )...;
        collect( nodes_ur_ui_Idx_in_state )... ]

    net_to_im_idx_conversion_dict =
        OrderedDict{Int64, Int64}(
            sys => ind for (sys, ind) in zip(
                im_vars_and_ur_ui_Idx_in_system,
                im_state_Idx  ) )

    return  (; im_vars_indices_in_system,
             pure_states_Idx_in_system,
             im_algebraic_vars_Idx_in_system,
             ur_ui_Idx_in_system,
             im_vars_and_ur_ui_Idx_in_system,
             im_vars_Idx_in_state,
             nodes_ur_ui_Idx_in_state,             
             im_state_Idx,
             each_gens_im_vars_Idx_in_state,
             net_to_im_idx_conversion_dict )
        
end


function get_net_to_im_indices_dict( netd )

    im_vars_indices_in_system =
        [get_gens_im_vars_indices_in_system(
            netd )...;]
    
    ur_ui_Idx_in_system = [
        get_nodes_ur_ui_indices_in_system(
            netd )...; ]

    dim_im_vars_and_ur_ui_vars = length.(
        [ im_vars_indices_in_system,
          ur_ui_Idx_in_system] )

    _, _, im_Idx =  create_size_offset_Idx(
        dim_im_vars_and_ur_ui_vars )

    im_vars_Idx  = im_Idx[1] 
    
    nodes_ur_ui_Idx = im_Idx[2]

    #----------------------------------------------------
    # indices transformation
    #----------------------------------------------------

    im_vars_and_ur_ui_Idx_in_system = [
        im_vars_indices_in_system...;
        ur_ui_Idx_in_system... ]
    
    im_state_Idx = [
        collect( im_vars_Idx )...;
        collect( nodes_ur_ui_Idx)... ]
    
    return OrderedDict{Int64, Int64}(
        sys => ind for (sys, ind) in
            zip( im_vars_and_ur_ui_Idx_in_system,
                 im_state_Idx  ) )
    
end


# ------------------------------------------------------
# Current balance
# ------------------------------------------------------


function get_Pl(Dyn_Nodes)
    return [a_node.Bus_type == :Load ?
        a_node.P :
        a_node.Bus_type == :Generator &&  a_node.with_loc_load == true ? a_node.P :
        0
            for a_node in collect(values(Dyn_Nodes)) ]
end

function get_Ql(Dyn_Nodes)
    return [a_node.Bus_type == :Load ?
        a_node.Q : a_node.Bus_type == :Generator &&  a_node.with_loc_load == true ? a_node.Q :
        0
            for a_node in collect(values(Dyn_Nodes)) ]
end


function get_Pl_Ql(Dyn_Nodes)
    
    return [a_node.Bus_type == :Load ?
        a_node.P + im * a_node.Q :
        a_node.Bus_type == :Generator &&  a_node.with_loc_load == true ? a_node.P + im * a_node.Q :
        0 + im * 0
            for a_node in collect(values(Dyn_Nodes)) ]
end


function get_ra(Dyn_Nodes)
    
    return [a_node.Bus_type == :Generator ?
        a_node.ra : 0
            for a_node in collect(values(Dyn_Nodes)) ]
end


function get_X_d_dash(Dyn_Nodes)
    
    return [a_node.Bus_type == :Generator ?
        a_node.X_d_dash : 0
            for a_node in collect(values(Dyn_Nodes)) ]
end


function get_X_q_dash(Dyn_Nodes)
    
    return [a_node.Bus_type == :Generator ?
        a_node.X_q_dash : 0
            for a_node in collect(values(Dyn_Nodes)) ]
end


function get_ra_X_d_q_dash(Dyn_Nodes)
    
    return Union{Float64, Vector{Float64,Float64,Float64}}[
        a_node.Bus_type == :Generator ?
            ( a_node.ra, a_node.X_d_dash, a_node.X_q_dash ) :
            0.0
        for a_node in collect(values(Dyn_Nodes)) ]
    
end

#------------------------------------------------
#------------------------------------------------



# function get_some_components_vars_in_a_plant(
#     plant ; some_vars =
#         [:δ, :ω, :ed_dash, :eq_dash] )

#     dict_state_vars_syms = plant.dict_state_vars_syms
#     dict_state_vars_syms
#     
    
#     dict_var_syms_Idx =
#         plant.dict_param_syms_Idx
        
#     return [ dict_state_vars_syms[a_sym]
#              for a_sym in some_vars ]    
        
#     return [ dict_param_syms_Idx[a_sym]
#              for a_sym in some_vars ]    

# end



function get_gens_nodes_some_state_vars_idxs_in_plants(
    nodes_collection;
    some_state_vars =
        [:δ, :ω, :ed_dash, :eq_dash] )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_some_state_vars_idxs_in_a_plant(
            comp; some_state_vars = some_state_vars )
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_some_state_vars_idxs_in_a_plant(
            comp; some_state_vars = some_state_vars )
                 for comp in nodes_collection
                     if comp.Bus_type == :Generator ]
    else
        
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ get_some_state_vars_idxs_in_a_plant(
            comp; some_state_vars = some_state_vars )
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    end

end



function get_gens_nodes_some_state_vars_dims(
    nodes_collection
    ; some_state_vars =
        [ :δ, :ω, :ed_dash, :eq_dash ] )

    gens_nodes = get_gens_nodes( nodes_collection  )

    some_state_vars_idxs_in_plants =
        get_gens_nodes_some_state_vars_idxs_in_plants(
            gens_nodes;
            some_state_vars = some_state_vars )

    return length.( some_state_vars_idxs_in_plants  )
end



function get_gens_nodes_some_state_vars_Idxs_in_flattend(
    nodes_collection
    ; some_state_vars =
        [ :δ, :ω, :ed_dash, :eq_dash ] )

    dims_some_state_vars =
        get_gens_nodes_some_state_vars_dims(
        nodes_collection
        ; some_state_vars =
            some_state_vars )
    
    _,_, some_state_vars_IDX = create_size_offset_Idx(
            dims_some_state_vars;
        counter = 0)
    return some_state_vars_IDX
end


function get_some_components_param_idx_in_a_plant(
    plant ; some_param =
        [:ra, :X_d, :X_q, :X_d_dash, :X_q_dash] )

    dict_param_syms_Idx =
        plant.dict_param_syms_Idx
        
    return [ dict_param_syms_Idx[a_sym]
             for a_sym in some_param ]    

end


function get_some_components_param_in_a_plant(
    plant ; some_param =
        [:ra, :X_d, :X_q, :X_d_dash, :X_q_dash] )

    dict_param_syms_Idx =
        plant.dict_param_syms_Idx
        
    params_idx = [
        dict_param_syms_Idx[a_sym]
        for a_sym in some_param ]

    return [plant.param_values...;][ params_idx ]
    

end


function get_nodes_some_param_dims(
    nodes_collection
    ; some_param =
        [:P, :Q] )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        return length.( [
            get_some_components_param_idx_in_a_plant(
                comp; some_param = some_param )
            for comp in comp_collection_values
               ])
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return length.( [
            get_some_components_param_idx_in_a_plant(
                comp; some_param = some_param )
            for comp in nodes_collection
                ])

    else
        
        comp_collection_values =
            collect(values(nodes_collection))
        
        return length.([
            get_some_components_param_idx_in_a_plant(
                comp; some_param = some_param )
            for comp in comp_collection_values
                 ])
        
    end

end


function get_nodes_some_param(
    nodes_collection
    ; some_param =
        [ :P, :Q ] )

    gens_nodes = get_gens_nodes( nodes_collection  )

    return [ get_some_components_param_in_a_plant(
                comp; some_param = some_param )
            for comp in gens_nodes ]

end


function get_gens_nodes_some_param(
    nodes_collection
    ; some_param =
        [:ra, :X_d, :X_q, :X_d_dash, :X_q_dash] )

    gens_nodes = get_gens_nodes( nodes_collection  )

    return [ get_some_components_param_in_a_plant(
        comp;
        some_param = some_param )
             for comp in
                 gens_nodes ]
end


function get_gens_nodes_some_param_dims(
    nodes_collection
    ; some_param =
        [:ra, :X_d, :X_q, :X_d_dash, :X_q_dash] )

    gens_nodes = get_gens_nodes( nodes_collection  )

    return length.( [
            get_some_components_param_idx_in_a_plant(
                comp; some_param = some_param )
            for comp in gens_nodes ])
end



function get_gens_nodes_with_loc_loads_some_param(
    nodes_collection
    ; some_param =
        [ :loc_P, :loc_Q ] )

    gens_nodes_with_loc_loads =
        get_gens_nodes_with_loc_loads( nodes_collection  )

    return [ gens_nodes_with_loc_loads == [] ?
        nothing :
        get_some_components_param_in_a_plant(
            comp; some_param = some_param )
             for comp in gens_nodes_with_loc_loads ]

end

function get_gens_nodes_with_loc_loads_some_param_dims(
    nodes_collection
    ; some_param =
        [ :loc_P, :loc_Q  ] )

    gens_nodes_with_loc_loads =
        get_gens_nodes_with_loc_loads( nodes_collection  )

    return gens_nodes_with_loc_loads == [] ?  length([]) : length.( [  get_some_components_param_idx_in_a_plant( comp; some_param = some_param ) for comp in gens_nodes_with_loc_loads  ])

end


function get_non_gens_nodes_some_param(
    nodes_collection
    ; some_param =
        [ :P, :Q ] )
    
    non_gens_nodes =
        get_non_gens_nodes( nodes_collection  )

    return [
            get_some_components_param_in_a_plant(
                comp; some_param = some_param )
            for comp in non_gens_nodes ]

end


function get_non_gens_nodes_some_param_dims(
    nodes_collection
    ; some_param =
        [:P, :Q ] )

    
    non_gens_nodes =
        get_non_gens_nodes( nodes_collection  )

    return length.([
            get_some_components_param_idx_in_a_plant(
                comp; some_param = some_param )
            for comp in non_gens_nodes ])

end


function get_load_nodes_some_param(
    nodes_collection
    ; some_param =
        [ :P, :Q ] )

    load_nodes =
        get_load_nodes( nodes_collection  )

    return [
            get_some_components_param_in_a_plant(
                comp; some_param = some_param )
            for comp in load_nodes ]

end



function get_load_nodes_some_param_dims(
    nodes_collection
    ; some_param =
        [:P, :Q] )


    load_nodes =
        get_load_nodes( nodes_collection  )

    return length.([
            get_some_components_param_idx_in_a_plant(
                comp; some_param = some_param )
            for comp in load_nodes  ])

end


function get_transmission_nodes_some_param(
    nodes_collection
    ; some_param =
        [ :P, :Q ] )

    transmission_nodes =
        get_transmission_nodes( nodes_collection  )

    return [
            get_some_components_param_in_a_plant(
                comp; some_param = some_param )
            for comp in transmission_nodes  ]
end



function get_transmission_nodes_some_param_dims(
    nodes_collection
    ; some_param =
        [:P, :Q, :Y_n] )


    transmission_nodes =
        get_transmission_nodes( nodes_collection  )

    return length.([
            get_some_components_param_idx_in_a_plant(
                comp; some_param = some_param )
            for comp in transmission_nodes])

end



#------------------------------------------------
#------------------------------------------------

function get_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(plant)

    return get_some_state_vars_idxs_in_a_plant(
        plant;
        some_state_vars =
            [:δ, :ω, :ed_dash, :eq_dash] )

end



function get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
    nodes_collection )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
                 for comp in nodes_collection
                     if comp.Bus_type == :Generator ]
    else
        
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ get_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    end

end


function get_ω_ed_dash_eq_dash_Idxs_in_a_plant(plant)

    return get_some_state_vars_idxs_in_a_plant(
        plant;
        some_state_vars = [:ω, :ed_dash, :eq_dash] )

end


function get_gens_nodes_ω_ed_dash_eq_dash_Idxs_in_plants(
    nodes_collection )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
                 for comp in nodes_collection
                     if comp.Bus_type == :Generator ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ get_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]        
    end

end



function get_nodes_δ_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    return [ get_a_node_δ_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
             for a_node in collect( values(netd.nodes ))
                 if a_node.Bus_type == :Generator ]
    
end 


function get_nodes_δ_ω_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    return [ get_a_node_δ_ω_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
             for a_node in collect( values(netd.nodes ))
                 if a_node.Bus_type == :Generator]
    
end 


function get_nodes_idx_and_δ_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    return [
        (a_node.Bus_num, get_a_node_δ_ed_dash_eq_dash_indices(
            ; network_vars_labels = network_vars_labels,
            bus_name = a_node.name ))
        for a_node in collect( values(netd.nodes ))
            if  a_node.Bus_type == :Generator  ]
    
end 



function get_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    return [ 
        (a_node.Bus_num, get_a_node_δ_ω_ed_dash_eq_dash_indices(
            ; network_vars_labels = network_vars_labels,
            bus_name = a_node.name )) 
             for a_node in collect( values(netd.nodes ))
                  if a_node.Bus_type == :Generator ]
    
end 



function get_nodes_idx_and_δ_ed_dash_eq_dash(
    state, netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    indices = [ get_a_node_δ_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
                for a_node in
                    collect(values(netd.nodes))]
    
    return Tuple{Int64, Vector{Float64}}[       
            (a_node.Bus_num, state[idx...]) 
        for (idx, a_node) in
            zip( indices,
                 collect(values(netd.nodes)) )
            if  a_node.Bus_type == :Generator ]
    
end



function get_nodes_idx_and_δ_ω_ed_dash_eq_dash(
    state, netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    indices = [ get_a_node_δ_ω_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
                for a_node in
                    collect(values(netd.nodes))]
    
    return Tuple{Int64, Vector{Float64}}[        
            (a_node.Bus_num, state[idx...]) 
        for (idx, a_node) in
            zip( indices, collect(values(netd.nodes)) )
            if  a_node.Bus_type == :Generator ]
    
end


function get_nodes_idx_and_δ_ed_dash_eq_dash_view(
    state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ed_dash_eq_dash_Idx(netd )

    
    return [ 
        (first(idx), view(state , last(idx))) 
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator]
    
end


function get_nodes_idx_and_δ_ed_dash_eq_dash(
    state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ed_dash_eq_dash_Idx(netd )

    
    return [ 
        (first(idx), state[last(idx)] ) 
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator]
    
end



function get_nodes_idx_and_δ_ω_ed_dash_eq_dash_view(
    state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    
    return [
        (first(idx), view(state , last(idx))) 
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator ]
    
end



function get_nodes_idx_and_δ_ω_ed_dash_eq_dash(
    state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    
    return [
        (first(idx), state[last(idx)] ) 
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator ]
    
end


function get_nodes_δ_ed_dash_eq_dash(state, netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    indices = [ get_a_node_δ_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
                for a_node in collect(values(netd.nodes))]
    
    return Vector{Float64}[ state[idx...] 
        for (idx, a_node) in
            zip( indices, collect(values(netd.nodes)) )
            if a_node.Bus_type == :Generator ]
    
end



function get_nodes_δ_ω_ed_dash_eq_dash(state, netd )

    network_vars_labels = get_network_vars_labels(netd )

    indices = [ get_a_node_δ_ω_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
                for a_node in
                    collect(values(netd.nodes))]
    
    return Vector{Float64}[ state[idx...] 
        for (idx, a_node) in
            zip( indices,
                 collect(values(netd.nodes)) )
            if a_node.Bus_type == :Generator ]
    
end


function get_nodes_δ_ed_dash_eq_dash_view(state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ed_dash_eq_dash_Idx(netd )

    
    return [
        view(state , last(idx)) 
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator]
    
end



function get_nodes_δ_ω_ed_dash_eq_dash_view(state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    
    return [ 
        view(state , last(idx))
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator ]
    
end


#------------------------------------------------
#------------------------------------------------


function get_gens_nodes_idx_and_δ_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    return [ 
        (a_node.Bus_num, get_a_node_δ_ed_dash_eq_dash_indices(
            ; network_vars_labels = network_vars_labels,
            bus_name = a_node.name ))
        for a_node in collect( values(netd.nodes ))
            if a_node.Bus_type == :Generator ]
    
end 



function get_gens_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    return [ 
        (a_node.Bus_num, get_a_node_δ_ω_ed_dash_eq_dash_indices(
            ; network_vars_labels = network_vars_labels,
            bus_name = a_node.name )) 
        for a_node in collect( values(netd.nodes ))
            if a_node.Bus_type == :Generator ]
    
end 



function get_gens_nodes_idx_and_δ_ed_dash_eq_dash(
    state, netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    indices = [ get_a_node_δ_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
                for a_node in
                    collect(values(netd.nodes))]
    
    return Tuple{Int64, Vector{Float64}}[        
            (a_node.Bus_num, state[idx...]) 
        for (idx, a_node) in
            zip( indices,
                 collect(values(netd.nodes)) )
            if a_node.Bus_type == :Generator]
    
end



function get_gens_nodes_idx_and_δ_ω_ed_dash_eq_dash(
    state, netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    indices = [ get_a_node_δ_ω_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
                for a_node in
                    collect(values(netd.nodes))]
    
    return Tuple{Int64, Vector{Float64}}[        
            (a_node.Bus_num, state[idx...]) 
        for (idx, a_node) in
            zip( indices, collect(values(netd.nodes)) )
            if a_node.Bus_type == :Generator ]
    
end


function get_gens_nodes_idx_and_δ_ed_dash_eq_dash_view(
    state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ed_dash_eq_dash_Idx(netd )

    
    return [
        (first(idx), view(state , last(idx)))
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator ]
    
end



function get_gens_nodes_idx_and_δ_ω_ed_dash_eq_dash_view(
    state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    
    return [
        (first(idx), view(state , last(idx))) 
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator ]
    
end


function get_gens_nodes_δ_ed_dash_eq_dash(state, netd )

    network_vars_labels =
        get_network_vars_labels(netd )

    indices = [ get_a_node_δ_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
                for a_node in collect(values(netd.nodes))]
    
    return Vector{Float64}[ state[idx...] 
        for (idx, a_node) in
            zip( indices, collect(values(netd.nodes)) )
            if a_node.Bus_type == :Generator]
    
end



function get_gens_nodes_δ_ω_ed_dash_eq_dash(state, netd )

    network_vars_labels = get_network_vars_labels(netd )

    indices = [ get_a_node_δ_ω_ed_dash_eq_dash_indices(
        ; network_vars_labels = network_vars_labels,
        bus_name = a_node.name )
                for a_node in
                    collect(values(netd.nodes))]
    
    return Vector{Float64}[ state[idx...]
        for (idx, a_node) in
            zip( indices,
                 collect(values(netd.nodes)) )
            if a_node.Bus_type == :Generator ]
    
end


function get_gens_nodes_δ_ed_dash_eq_dash_view(state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ed_dash_eq_dash_Idx(netd )

    
    return [
        view(state , last(idx))
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator ]
    
end



function get_gens_nodes_δ_ω_ed_dash_eq_dash_view(
    state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    
    return [
        view(state , last(idx)) 
            for (idx, a_node) in
                zip(nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx,
                    comp_collection )
                if a_node.Bus_type == :Generator ]
    
end


#------------------------------------------------


# function get_gens_nodes_idx_and_δ_ω_ed_dash_eq_dash_view(
#     state, netd )

#     comp_collection =
#         collect(values(netd.nodes) )
    
#     nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
#         get_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    
#     return [
#         (first(idx), view(state , last(idx)))
#             for (idx, a_node) in
#                 zip(nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx,
#                     comp_collection )
#                 if a_node.Bus_type == :Generator  ]
    
# end

#------------------------------------------------
#------------------------------------------------

function get_hybrid_pf_etc_idx_and_Idx(netd)

    #------------------------------------------

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
    
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   = [ length(ur_idx), length(ui_idx) ]
    ur_ui_offset = create_offsets( ur_ui_dims )
    ur_ui_IDX    = create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    # ----------------------------------------------------

    slack_bus_idx =
        get_components_slack_Idx( netd.nodes )
    
    gens_idx = first.(
        get_generators_Idx_and_vh(
            netd.nodes ))        

     # ----------------------------------------------------

     red_vh_θh_idx = [ setdiff(vh_IDX, gens_idx)...;
                setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ---------------------------------------------------   
    
     red_vh_idx = setdiff(vh_IDX, gens_idx)

     red_θh_idx =
         setdiff( θh_IDX, θh_IDX[slack_bus_idx] )

    # ----------------------------------------------------

    red_vh_θh_dims =
        [ length( red_vh_idx ), length( red_θh_idx ) ]    

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(red_vh_θh_dims )

    red_vh_Idxs = red_vh_θh_IDX[1]
    
    red_θh_Idxs = red_vh_θh_IDX[2]
    
    # ----------------------------------------------------

     slack_gens_nodes_idx =
         get_slack_gens_nodes_idx(
             netd.nodes )

     non_slack_gens_nodes_idx =
         get_non_slack_gens_nodes_idx(
             netd.nodes )

     gens_nodes_idx =
         get_gens_nodes_idx( netd.nodes )

     gens_nodes_with_loc_loads_idx =
         get_gens_nodes_with_loc_loads_idx(
             netd.nodes )
    
    loc_load_exist =
        gens_nodes_with_loc_loads_idx == [] ?
        false : true
    
     load_nodes_idx =
         get_load_nodes_idx( netd.nodes )

     transmission_nodes_idx =
         get_transmission_nodes_idx( netd.nodes )

     non_gens_nodes_idx =
         get_non_gens_nodes_idx( netd.nodes )

     load_trans_nodes_idx =
         get_load_trans_nodes_Idx( netd.nodes )

     all_nodes_idx =
         get_all_nodes_idx( netd.nodes  )

     non_slack_gens_and_non_gens_idx =
         setdiff(all_nodes_idx,
                 slack_gens_nodes_idx)

     # --------------------------------------------------   

    full_vh_θh_dims =
        [ length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx ),
          length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx )
          ]
    
     _, _, full_vh_θh_IDX =
        create_size_offset_Idx( full_vh_θh_dims )

    full_slack_gens_vh_Idxs =
        full_vh_θh_IDX[1]
    
    full_non_slack_gens_vh_Idxs =
        full_vh_θh_IDX[2]

    full_non_gens_nodes_vh_Idxs =
        full_vh_θh_IDX[3]
    
    full_slack_gens_θh_Idxs =
        full_vh_θh_IDX[4]

    full_non_slack_gens_θh_Idxs =
        full_vh_θh_IDX[5]
    
    full_non_gens_nodes_θh_Idxs =
        full_vh_θh_IDX[6]

    full_nodes_types_vh_Idxs =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs)

    full_nodes_types_θh_Idxs =
        (;full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs)

    full_nodes_types_vh_and_θh_Idxs =
        (; full_nodes_types_vh_Idxs,
         full_nodes_types_θh_Idxs)
    
     # --------------------------------------------------  
    
    nodes_type_idxs =
        [slack_gens_nodes_idx;
         non_slack_gens_nodes_idx;
         non_gens_nodes_idx ]


    full_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      vh_IDX) )
    
    full_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      θh_IDX ) )
    
    vec_types_full_vh_idx2Idx = [
        [ full_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx, full_non_slack_gens_vh_idx2Idx, full_non_gens_vh_idx2Idx = vec_types_full_vh_idx2Idx
        

    vec_types_full_θh_idx2Idx = [
        [ full_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx, full_non_slack_gens_θh_idx2Idx, full_non_gens_θh_idx2Idx = vec_types_full_θh_idx2Idx


    full_nodes_types_vh_idx2Idx =
        (; full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx)
    
    full_nodes_types_θh_idx2Idx =
        (; full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx)

    full_nodes_types_vh_and_θh_idx2Idx =
        (; full_nodes_types_vh_idx2Idx,
         full_nodes_types_θh_idx2Idx)

    full_nodes_types_dict_vh_and_θh_idx2Idx =
        (; full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx )

    
     # --------------------------------------------------   

    
    full_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      1:length(vh_IDX) ) )

    
    full_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      1:length( θh_IDX ) ) )

    
    vec_types_full_vh_idx2Idx_in_Idx = [
        [ full_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx_in_Idx, full_non_slack_gens_vh_idx2Idx_in_Idx, full_non_gens_vh_idx2Idx_in_Idx = vec_types_full_vh_idx2Idx_in_Idx


    vec_types_full_θh_idx2Idx_in_Idx = [
        [ full_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx_in_Idx, full_non_slack_gens_θh_idx2Idx_in_Idx, full_non_gens_θh_idx2Idx_in_Idx = vec_types_full_θh_idx2Idx_in_Idx


    full_nodes_types_vh_idx2Idx_in_Idx =
        (; full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx   )
    
    full_nodes_types_θh_idx2Idx_in_Idx =
        (; full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx )

    full_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; full_nodes_types_vh_idx2Idx_in_Idx,
         full_nodes_types_θh_idx2Idx_in_Idx )
    
    full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )

     # --------------------------------------------------   

    # full_nodes_types_Idxs_idx2Idx_etc =
    #     (;full_nodes_types_vh_and_θh_Idxs,    
    #      full_nodes_types_dict_vh_and_θh_idx2Idx,    
    #      full_nodes_types_vh_and_θh_idx2Idx,    
    #      full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx,    
    #      full_nodes_types_vh_and_θh_idx2Idx_in_Idx )



    full_nodes_types_Idxs_idx2Idx_etc =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs,
         full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx,    
         full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx,
         full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx,
         full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx,
         full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx,    
         full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx
         )
    
    
    # --------------------------------------------------   

     dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      1:length(red_θh_Idxs) ) )

    vec_red_nodes_types_idxs =
        [ non_slack_gens_nodes_idx,
          non_gens_nodes_idx ]

    vec_types_red_idx2Idx = [
        [ dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            vec_red_nodes_types_idxs  ]

    non_slack_gens_θh_idx2Idx, non_gens_θh_idx2Idx =
        vec_types_red_idx2Idx
    
    # non_slack_gens_θh_idx2Idx =
    #      [ dict_θh_idx2Idx[idx]
    #       for idx in
    #           non_slack_gens_nodes_idx ]

    #  non_gens_θh_idx2Idx =
    #      [ dict_θh_idx2Idx[idx]
    #       for idx in
    #           non_gens_nodes_idx ]

    # -----------------------------------


    vec_types_red_idx2Idx_in_Idx = [
        [ dict_θh_idx2Idx_in_Idx[ idx ]
          for idx in types_idxs ]
        for types_idxs in
            vec_red_nodes_types_idxs  ]

    non_slack_gens_θh_idx2Idx_in_Idx, non_gens_θh_idx2Idx_in_Idx = vec_types_red_idx2Idx_in_Idx
    
     # non_slack_gens_θh_idx2Idx_in_Idx =
     #     [ dict_θh_idx2Idx_in_Idx[idx]
     #       for idx in
     #          non_slack_gens_nodes_idx ]
    
     # non_gens_θh_idx2Idx_in_Idx =
     #     [ dict_θh_idx2Idx_in_Idx[idx]
     #       for idx in
     #          non_gens_nodes_idx ]

    # ----------------------------------------------------
    # red_vh_idx, red_θh_idx, red_vh_Idxs, red_θh_Idxs
    # non_slack_gens_θh_idx2Idx,
    # non_slack_gens_θh_idx2Idx_in_Idx,
    # non_gens_θh_idx2Idx, non_gens_θh_idx2Idx_in_Idx
    # ----------------------------------------------------
    
    return (; slack_ur_ui_Idx_in_state,
            non_slack_ur_ui_Idx_in_state,
            ur_ui_Idx_in_state,
            nodes_u_Idx,
            ur_idx,
            ui_idx,
            ur_ui_idx,
            vh_θh_idx,
            ur_ui_IDX,
            ur_IDX,
            ui_IDX,
            ir_IDX,
            ii_IDX,
            vh_IDX,
            θh_IDX,
            slack_bus_idx,
            gens_idx,
            red_vh_θh_idx,
            red_vh_idx,
            red_θh_idx,
            red_vh_Idxs,
            red_θh_Idxs,
            non_slack_gens_θh_idx2Idx,
            non_slack_gens_θh_idx2Idx_in_Idx,
            non_gens_θh_idx2Idx,
            non_gens_θh_idx2Idx_in_Idx,
            loc_load_exist,
            full_nodes_types_Idxs_idx2Idx_etc)
    
    
end


function get_im_model_pf_idx_and_Idx(netd)


    #------------------------------------------        
    
    dict_sys_to_im  =
        get_net_to_im_indices_dict( netd  )

    #------------------------------------------        

    slack_ur_ui_Idx_in_state =
        net_to_im_model_indices(
            get_components_slack_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_im   )

    non_slack_ur_ui_Idx_in_state =
        net_to_im_model_indices(
            get_components_no_slack_ur_ui_Idx_in_state(netd.nodes),
            dict_sys_to_im   )

    ur_ui_Idx_in_state =
        net_to_im_model_indices(
            get_components_ur_ui_Idx_in_state( netd.nodes ),
            dict_sys_to_im   )

    nodes_u_Idx =
        net_to_im_model_indices(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_im   )
    
    #------------------------------------------        
    #------------------------------------------        
        
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   =
        [ length(ur_idx), length(ui_idx) ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    # ----------------------------------------------------

    slack_bus_idx =
        get_components_slack_Idx( netd.nodes )
    
    gens_idx = first.(
        get_generators_Idx_and_vh(
            netd.nodes ))        

    # ----------------------------------------------------

    red_vh_θh_idx =
        [ setdiff(vh_IDX, gens_idx)...;
          setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]

    # ----------------------------------------------------

    red_vh_idx = setdiff(vh_IDX, gens_idx)
    
    red_θh_idx =
        setdiff(θh_IDX, θh_IDX[slack_bus_idx])

    # ----------------------------------------------------

    red_vh_θh_dims =
        [ length( red_vh_idx ), length( red_θh_idx ) ]    

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(red_vh_θh_dims )

    red_vh_Idxs = red_vh_θh_IDX[1]
    
    red_θh_Idxs = red_vh_θh_IDX[2]

    # ----------------------------------------------------
    
     slack_gens_nodes_idx =
         get_slack_gens_nodes_idx(
             netd.nodes )

     non_slack_gens_nodes_idx =
         get_non_slack_gens_nodes_idx(
             netd.nodes )

     gens_nodes_idx =
         get_gens_nodes_idx( netd.nodes )

     gens_nodes_with_loc_loads_idx =
         get_gens_nodes_with_loc_loads_idx(
             netd.nodes )
    
    loc_load_exist =
        gens_nodes_with_loc_loads_idx == [] ?
        false : true

     load_nodes_idx =
         get_load_nodes_idx( netd.nodes )

     transmission_nodes_idx =
         get_transmission_nodes_idx( netd.nodes )

     non_gens_nodes_idx =
         get_non_gens_nodes_idx( netd.nodes )

     load_trans_nodes_idx =
         get_load_trans_nodes_Idx( netd.nodes )

     all_nodes_idx =
         get_all_nodes_idx( netd.nodes  )

     non_slack_gens_and_non_gens_idx =
         setdiff(all_nodes_idx,
                 slack_gens_nodes_idx)

     # --------------------------------------------------   

    full_vh_θh_dims =
        [ length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx ),
          length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx )
          ]
    
     _, _, full_vh_θh_IDX =
        create_size_offset_Idx( full_vh_θh_dims )

    full_slack_gens_vh_Idxs =
        full_vh_θh_IDX[1]
    
    full_non_slack_gens_vh_Idxs =
        full_vh_θh_IDX[2]

    full_non_gens_nodes_vh_Idxs =
        full_vh_θh_IDX[3]
    
    full_slack_gens_θh_Idxs =
        full_vh_θh_IDX[4]

    full_non_slack_gens_θh_Idxs =
        full_vh_θh_IDX[5]
    
    full_non_gens_nodes_θh_Idxs =
        full_vh_θh_IDX[6]

    full_nodes_types_vh_Idxs =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs)

    full_nodes_types_θh_Idxs =
        (;full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs)

    full_nodes_types_vh_and_θh_Idxs =
        (; full_nodes_types_vh_Idxs,
         full_nodes_types_θh_Idxs)
    
     # --------------------------------------------------  
    
    nodes_type_idxs =
        [slack_gens_nodes_idx;
         non_slack_gens_nodes_idx;
         non_gens_nodes_idx ]


    full_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      vh_IDX) )
    
    full_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      θh_IDX ) )
    
    vec_types_full_vh_idx2Idx = [
        [ full_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx, full_non_slack_gens_vh_idx2Idx, full_non_gens_vh_idx2Idx = vec_types_full_vh_idx2Idx
        

    vec_types_full_θh_idx2Idx = [
        [ full_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx, full_non_slack_gens_θh_idx2Idx, full_non_gens_θh_idx2Idx = vec_types_full_θh_idx2Idx


    full_nodes_types_vh_idx2Idx =
        (; full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx)
    
    full_nodes_types_θh_idx2Idx =
        (; full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx)

    full_nodes_types_vh_and_θh_idx2Idx =
        (; full_nodes_types_vh_idx2Idx,
         full_nodes_types_θh_idx2Idx)

    full_nodes_types_dict_vh_and_θh_idx2Idx =
        (; full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx )

    
     # --------------------------------------------------   

    
    full_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      1:length(vh_IDX) ) )

    
    full_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      1:length( θh_IDX ) ) )

    
    vec_types_full_vh_idx2Idx_in_Idx = [
        [ full_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx_in_Idx, full_non_slack_gens_vh_idx2Idx_in_Idx, full_non_gens_vh_idx2Idx_in_Idx = vec_types_full_vh_idx2Idx_in_Idx


    vec_types_full_θh_idx2Idx_in_Idx = [
        [ full_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx_in_Idx, full_non_slack_gens_θh_idx2Idx_in_Idx, full_non_gens_θh_idx2Idx_in_Idx = vec_types_full_θh_idx2Idx_in_Idx


    full_nodes_types_vh_idx2Idx_in_Idx =
        (; full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx   )
    
    full_nodes_types_θh_idx2Idx_in_Idx =
        (; full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx )

    full_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; full_nodes_types_vh_idx2Idx_in_Idx,
         full_nodes_types_θh_idx2Idx_in_Idx )
    
    full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )

     # --------------------------------------------------   

    # full_nodes_types_Idxs_idx2Idx_etc =
    #     (;full_nodes_types_vh_and_θh_Idxs,    
    #      full_nodes_types_dict_vh_and_θh_idx2Idx,    
    #      full_nodes_types_vh_and_θh_idx2Idx,    
    #      full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx,    
    #      full_nodes_types_vh_and_θh_idx2Idx_in_Idx )

    

    full_nodes_types_Idxs_idx2Idx_etc =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs,
         full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx,    
         full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx,
         full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx,
         full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx,
         full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx,    
         full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx
         )
    
    
     # --------------------------------------------------   

     dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      1:length(red_θh_Idxs) ) )

     non_slack_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_slack_gens_nodes_idx ]

     non_slack_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_slack_gens_nodes_idx ]


     non_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_gens_nodes_idx ]

     non_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_gens_nodes_idx ]

    # ----------------------------------------------------

    # red_vh_idx, red_θh_idx, red_vh_Idxs, red_θh_Idxs
    # non_slack_gens_θh_idx2Idx,
    # non_slack_gens_θh_idx2Idx_in_Idx,
    # non_gens_θh_idx2Idx, non_gens_θh_idx2Idx_in_Idx
    # ----------------------------------------------------
    
    return (; slack_ur_ui_Idx_in_state,
            non_slack_ur_ui_Idx_in_state,
            ur_ui_Idx_in_state,
            nodes_u_Idx,
            ur_idx,
            ui_idx,
            ur_ui_idx,
            vh_θh_idx,
            ur_ui_IDX,
            ur_IDX,
            ui_IDX,
            ir_IDX,
            ii_IDX,
            vh_IDX,
            θh_IDX,
            slack_bus_idx,
            gens_idx,
            red_vh_θh_idx,
            red_vh_idx,
            red_θh_idx,
            red_vh_Idxs,
            red_θh_Idxs,
            non_slack_gens_θh_idx2Idx,
            non_slack_gens_θh_idx2Idx_in_Idx,
            non_gens_θh_idx2Idx,
            non_gens_θh_idx2Idx_in_Idx,
            loc_load_exist,
            full_nodes_types_Idxs_idx2Idx_etc)            


end



function get_industrial_model_pf_idx_and_Idx(
    netd; no_control_device = false )

    #------------------------------------------        


    if no_control_device == false

        dict_sys_to_industry =
            get_net_to_industrial_model_indices_dict( netd  )
        
    else
        dict_sys_to_industry =
            get_net_to_industrial_model_indices_dict(
                netd;  no_control_device = true  )
    end
    

    #------------------------------------------        

    slack_ur_ui_Idx_in_state =
        net_to_industrial_model_indices(
            get_components_slack_ur_ui_Idx_in_state(
                netd.nodes ), dict_sys_to_industry  )

    non_slack_ur_ui_Idx_in_state =
        net_to_industrial_model_indices(
            get_components_no_slack_ur_ui_Idx_in_state(netd.nodes),
            dict_sys_to_industry  )

    ur_ui_Idx_in_state =
        net_to_industrial_model_indices(
            get_components_ur_ui_Idx_in_state( netd.nodes ),
            dict_sys_to_industry  )

    nodes_u_Idx =
        net_to_industrial_model_indices(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )
    
    #------------------------------------------        
    #------------------------------------------        
        
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   =
        [ length(ur_idx), length(ui_idx) ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    # ----------------------------------------------------

    slack_bus_idx =
        get_components_slack_Idx( netd.nodes )
    
    gens_idx = first.(
        get_generators_Idx_and_vh(
            netd.nodes ))        

    # ----------------------------------------------------
    
    red_vh_θh_idx =
        [ setdiff(vh_IDX, gens_idx)...;
          setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ----------------------------------------------------

    red_vh_idx = setdiff(vh_IDX, gens_idx)
    
    red_θh_idx =
        setdiff(θh_IDX, θh_IDX[slack_bus_idx])

    # ----------------------------------------------------

    red_vh_θh_dims =
        [ length( red_vh_idx ), length( red_θh_idx ) ]    

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(red_vh_θh_dims )

    red_vh_Idxs = red_vh_θh_IDX[1]
    
    red_θh_Idxs = red_vh_θh_IDX[2]

    # ----------------------------------------------------

    
    # ----------------------------------------------------

     slack_gens_nodes_idx =
         get_slack_gens_nodes_idx(
             netd.nodes )

     non_slack_gens_nodes_idx =
         get_non_slack_gens_nodes_idx(
             netd.nodes )

     gens_nodes_idx =
         get_gens_nodes_idx( netd.nodes )

     gens_nodes_with_loc_loads_idx =
         get_gens_nodes_with_loc_loads_idx(
             netd.nodes )
    
    loc_load_exist =
        gens_nodes_with_loc_loads_idx == [] ?
        false : true

     load_nodes_idx =
         get_load_nodes_idx( netd.nodes )

     transmission_nodes_idx =
         get_transmission_nodes_idx( netd.nodes )

     non_gens_nodes_idx =
         get_non_gens_nodes_idx( netd.nodes )

     load_trans_nodes_idx =
         get_load_trans_nodes_Idx( netd.nodes )

     all_nodes_idx =
         get_all_nodes_idx( netd.nodes  )

     non_slack_gens_and_non_gens_idx =
         setdiff(all_nodes_idx, slack_gens_nodes_idx)

     # --------------------------------------------------   

    full_vh_θh_dims =
        [ length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx ),
          length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx )
          ]
    
     _, _, full_vh_θh_IDX =
        create_size_offset_Idx( full_vh_θh_dims )

    full_slack_gens_vh_Idxs =
        full_vh_θh_IDX[1]
    
    full_non_slack_gens_vh_Idxs =
        full_vh_θh_IDX[2]

    full_non_gens_nodes_vh_Idxs =
        full_vh_θh_IDX[3]
    
    full_slack_gens_θh_Idxs =
        full_vh_θh_IDX[4]

    full_non_slack_gens_θh_Idxs =
        full_vh_θh_IDX[5]
    
    full_non_gens_nodes_θh_Idxs =
        full_vh_θh_IDX[6]

    full_nodes_types_vh_Idxs =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs )

    full_nodes_types_θh_Idxs =
        (;full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs)

    full_nodes_types_vh_and_θh_Idxs =
        (; full_nodes_types_vh_Idxs,
         full_nodes_types_θh_Idxs)
    
     # --------------------------------------------------  
    
    nodes_type_idxs =
        [ slack_gens_nodes_idx;
         non_slack_gens_nodes_idx;
         non_gens_nodes_idx ]


    full_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      vh_IDX) )
    
    full_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      θh_IDX ) )
    
    vec_types_full_vh_idx2Idx = [
        [ full_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx, full_non_slack_gens_vh_idx2Idx, full_non_gens_vh_idx2Idx = vec_types_full_vh_idx2Idx
        

    vec_types_full_θh_idx2Idx = [
        [ full_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx, full_non_slack_gens_θh_idx2Idx, full_non_gens_θh_idx2Idx = vec_types_full_θh_idx2Idx


    full_nodes_types_vh_idx2Idx =
        (; full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx )
    
    full_nodes_types_θh_idx2Idx =
        (; full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx )

    full_nodes_types_vh_and_θh_idx2Idx =
        (; full_nodes_types_vh_idx2Idx,
         full_nodes_types_θh_idx2Idx )

    full_nodes_types_dict_vh_and_θh_idx2Idx =
        (; full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx )

    
     # --------------------------------------------------   

    
    full_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      1:length(vh_IDX) ) )

    
    full_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( nodes_type_idxs,
                      1:length( θh_IDX ) ) )

    
    vec_types_full_vh_idx2Idx_in_Idx = [
        [ full_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx_in_Idx, full_non_slack_gens_vh_idx2Idx_in_Idx, full_non_gens_vh_idx2Idx_in_Idx = vec_types_full_vh_idx2Idx_in_Idx


    vec_types_full_θh_idx2Idx_in_Idx = [
        [ full_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx_in_Idx, full_non_slack_gens_θh_idx2Idx_in_Idx, full_non_gens_θh_idx2Idx_in_Idx = vec_types_full_θh_idx2Idx_in_Idx


    full_nodes_types_vh_idx2Idx_in_Idx =
        (; full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx   )
    
    full_nodes_types_θh_idx2Idx_in_Idx =
        (; full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx )

    full_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; full_nodes_types_vh_idx2Idx_in_Idx,
         full_nodes_types_θh_idx2Idx_in_Idx )
    
    full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )

     # --------------------------------------------------   

    # full_nodes_types_Idxs_idx2Idx_etc =
    #     (;full_nodes_types_vh_and_θh_Idxs,    
    #      full_nodes_types_dict_vh_and_θh_idx2Idx,    
    #      full_nodes_types_vh_and_θh_idx2Idx,    
    #      full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx,    
    #      full_nodes_types_vh_and_θh_idx2Idx_in_Idx )


    full_nodes_types_Idxs_idx2Idx_etc =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs,
         full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx,    
         full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx,
         full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx,
         full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx,
         full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx,    
         full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )
    
     # --------------------------------------------------   
     # --------------------------------------------------   
    
     dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      1:length(red_θh_Idxs) ) )

     non_slack_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_slack_gens_nodes_idx ]

     non_slack_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_slack_gens_nodes_idx ]


     non_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_gens_nodes_idx ]

     non_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_gens_nodes_idx ]

    # ----------------------------------------------------

    # red_vh_idx, red_θh_idx, red_vh_Idxs, red_θh_Idxs
    # non_slack_gens_θh_idx2Idx,
    # non_slack_gens_θh_idx2Idx_in_Idx,
    # non_gens_θh_idx2Idx, non_gens_θh_idx2Idx_in_Idx
    # ----------------------------------------------------
    
    return (; slack_ur_ui_Idx_in_state,
            non_slack_ur_ui_Idx_in_state,
            ur_ui_Idx_in_state,
            nodes_u_Idx,
            ur_idx,
            ui_idx,
            ur_ui_idx,
            vh_θh_idx,
            ur_ui_IDX,
            ur_IDX,
            ui_IDX,
            ir_IDX,
            ii_IDX,
            vh_IDX,
            θh_IDX,
            slack_bus_idx,
            gens_idx,
            red_vh_θh_idx,
            red_vh_idx,
            red_θh_idx,
            red_vh_Idxs,
            red_θh_Idxs,
            non_slack_gens_θh_idx2Idx,
            non_slack_gens_θh_idx2Idx_in_Idx,
            non_gens_θh_idx2Idx,
            non_gens_θh_idx2Idx_in_Idx,
            loc_load_exist,
            full_nodes_types_Idxs_idx2Idx_etc )            
    
end

#-----------------------------------------------

function get_a_hybrid_state_idx_and_Idx(netd)

     
    #-----------------------------------------
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(netd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(netd.edges)
    nodes_incident_edges = get_nodes_incident_edges(netd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(netd.nodes, netd.edges)

    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset = create_offsets(nodes_dims)

    

    nodes_u_Idx = get_components_ur_ui_Idx_in_state(netd.nodes)
    
    # nodes_u_Idx = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(netd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------
    
    """
    These are needed for voltages at src and dst of an edge in edges_src_volt_view,
    and edges_dst_volt_view an edge has two ends `h` and `k`. `from` nodes are attached
    to the `h` end of an edge, while `to` nodes are attached to the `k` end
    """
    edges_ih_offset  = create_offsets(edges_dims; counter = nodes_dim_size)
    edges_ih_Idx     = create_x_src_idxs(edges_ih_offset)

    edges_ik_offset  = create_offsets(edges_dims; counter = nodes_dim_size)
    edges_ik_Idx     = create_x_dst_idxs(edges_ik_offset)

    #--------------------------------------------------

    gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants =
        get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
            netd.nodes )

    gens_nodes_δ_ω_ed_dash_eq_dash_Idx =
        get_nodes_δ_ω_ed_dash_eq_dash_Idx(netd )

    gens_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
        get_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(
            netd )
    
    #--------------------------------------------------
    
    return (;gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants,
            gens_nodes_δ_ω_ed_dash_eq_dash_Idx,
            gens_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx,
            nodes_state_Idx,
            edges_state_Idx,
            nodes_u_Idx,
            edges_ih_Idx,
            edges_ik_Idx,
            nodes_src_edges,
            nodes_dst_edges,
            nodes_incident_edges,
            edges_orientation,
            edges_src_node,
            edges_dst_node)   

end



function get_a_hybrid_para_idx_and_Idx(netd)
     
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(netd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(netd.edges)
    nodes_incident_edges = get_nodes_incident_edges(netd.edges)
    
    edges_orientation   = get_edges_orientation(netd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx =
        create_comps_non_flattened_vectors_Idx(
            nodes_param)

      
    edges_param_Idx =
        create_comps_non_flattened_vectors_Idx(
            edges_param)

    #---------------------------------------------------
    
    gens_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_Idxs =
        get_gens_params_idxs_in_param_values(
            netd.nodes;
            param_list =
                [ :ra, :X_d, :X_q, :X_d_dash, :X_q_dash ])

        
    gens_nodes_ra_Xd_Xq_Xd_Idxs =
        get_gens_params_idxs_in_param_values(
            netd.nodes;
            param_list =
                [ :ra, :X_d, :X_q ])

    
    gens_nodes_ra_Xd_dash_Xq_dash_Idxs =
        get_gens_params_idxs_in_param_values(
            netd.nodes;
            param_list =
                [ :ra, :X_d_dash, :X_q_dash ])
    

    #------------------------------------------------
    
    """
    These are needed for voltages at src and dst of an edge in edges_src_volt_view,
    and edges_dst_volt_view an edge has two ends `h` and `k`. `from` nodes are attached
    to the `h` end of an edge, while `to` nodes are attached to the `k` end
    """
    edges_ih_offset  = create_offsets(edges_dims; counter = nodes_dim_size)
    edges_ih_Idx     = create_x_src_idxs(edges_ih_offset)

    edges_ik_offset  = create_offsets(edges_dims; counter = nodes_dim_size)
    edges_ik_Idx     = create_x_dst_idxs(edges_ik_offset)

    return (; gens_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_Idxs,
            gens_nodes_ra_Xd_Xq_Xd_Idxs,
            gens_nodes_ra_Xd_dash_Xq_dash_Idxs,
            nodes_param_Idx,
            edges_param_Idx,            
            edges_ih_Idx,
            edges_ik_Idx,
            nodes_src_edges,
            nodes_dst_edges,
            nodes_incident_edges,
            edges_orientation,
            edges_src_node,
            edges_dst_node)
   
end


#------------------------------------------------

function get_pf_and_dyn_idx_and_Idx(netd)
    
    #------------------------------------------            
    _, _, edges_orientation,  nodes_node_idx_and_incident_edges_other_node_idx = get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------        

    hybrid_pf_etc_idx_and_Idx =
        get_hybrid_pf_etc_idx_and_Idx(netd)
    
    #------------------------------------------

    im_model_pf_idx_and_Idx =
        get_im_model_pf_idx_and_Idx(netd)

    #------------------------------------------
    
    industrial_model_pf_idx_and_Idx =
        get_industrial_model_pf_idx_and_Idx(
            netd; no_control_device = false )
    
    #------------------------------------------

    only_gens_industrial_model_pf_idx_and_Idx =
        get_industrial_model_pf_idx_and_Idx(
            netd; no_control_device = true )
    
    #------------------------------------------    
    #------------------------------------------    

    slack_gens_nodes_idx =
        get_slack_gens_nodes_idx(
            netd.nodes )

    non_slack_gens_nodes_idx =
        get_non_slack_gens_nodes_idx(
            netd.nodes )
    
    gens_nodes_idx =
        get_gens_nodes_idx( netd.nodes )

    gens_nodes_with_loc_loads_idx =
        get_gens_nodes_with_loc_loads_idx(
            netd.nodes )
    
    loc_load_exist =
        gens_nodes_with_loc_loads_idx == [] ?
        false : true
        
    load_nodes_idx =
        get_load_nodes_idx( netd.nodes )

    transmission_nodes_idx =
        get_transmission_nodes_idx( netd.nodes )
    
    non_gens_nodes_idx =
        get_non_gens_nodes_idx( netd.nodes )
    
    load_trans_nodes_idx =
        get_load_trans_nodes_Idx( netd.nodes )


     all_nodes_idx =
         get_all_nodes_idx( netd.nodes  )


    #------------------------------------------        
    
        gens_with_loc_load_idx = first.(
            get_gens_with_loc_load_Idx_and_loc_load(
                netd.nodes ) )

    # ----------------------------------------------------
    
    vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
                netd.nodes) )

    
    vec_Idx_gens_nodes_ω_ed_dash_eq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_ω_ed_dash_eq_dash_Idxs_in_plants(
                netd.nodes) )
    
    
    # ----------------------------------------------------

    vec_Idx_ra_Xd_Xq_Xd_dash_Xq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param_dims(
                netd.nodes
                ;some_param =
                    [:ra, :X_d, :X_q,
                     :X_d_dash, :X_q_dash]  ))


    vec_Idx_ra_Xd_dash_Xq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param_dims(
                netd.nodes
                ;some_param =
                    [ :ra, :X_d_dash, :X_q_dash ] ))


    vec_Idx_ra_Xd_Xq =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param_dims(
                netd.nodes
                ;some_param =
                    [ :ra, :X_d, :X_q ] ))
    
    # ------------------------------------------------
    
    vec_Idx_P_Q_nodes =
        get_flattened_to_components_vector_var_Idx(
            get_nodes_some_param_dims(
                netd.nodes
                ; some_param = [:P, :Q] ) )

    
    vec_Idx_P_Q_gens =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P, :Q ] ) )

    
    vec_Idx_P_Q_non_gens =
        get_flattened_to_components_vector_var_Idx(
            get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P, :Q ] ) )

    
    #------------------------------------------
    
    gens_nodes_with_loc_loads_some_param_dims =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param =
                [ :loc_P, :loc_Q  ] )
    
    vec_Idx_P_Q_gens_loc_load =
        gens_nodes_with_loc_loads_some_param_dims == 0 ?
        nothing :
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_with_loc_loads_some_param_dims(
                netd.nodes
                ; some_param = [ :loc_P, :loc_Q  ] ) )

    #------------------------------------------
    #------------------------------------------

    dim_P_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :Q ] ) )
    
    dim_P_load =
        length( get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_load =
        length( get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :Q ] ) )       
    
    dim_δ_ed_eq_pf =
        length( [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;] )

    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_P ] ) == 0 ? 0 :
        length( get_gens_nodes_with_loc_loads_some_param_dims( netd.nodes; some_param = [ :loc_P ] )  )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
        length( get_gens_nodes_with_loc_loads_some_param_dims( netd.nodes; some_param = [ :loc_Q  ] )  )
    
    #----------------------------------------
    #----------------------------------------
    
    dim_with_loc_load_with_δ_ed_NL_pf_para =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_load,
          dim_Q_load,
          dim_δ_ed_eq_pf,
          dim_P_g_loc_load,
          dim_Q_g_loc_load  ]

    _,_, with_loc_load_with_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_with_loc_load_with_δ_ed_NL_pf_para )

    with_loc_load_with_δ_ed_P_gens_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[1]

    with_loc_load_with_δ_ed_Q_gens_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[2]

    with_loc_load_with_δ_ed_P_load_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[3]

    with_loc_load_with_δ_ed_Q_load_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[4]

    with_loc_load_δ_ed_eq_pf_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[5]

    with_loc_load_with_δ_ed_P_g_loc_load_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[6]

    with_loc_load_with_δ_ed_Q_g_loc_load_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[7]

    with_loc_load_with_δ_ed_NL_pf_para_Idxs = (
        ; with_loc_load_with_δ_ed_P_gens_NL_para_Idxs,
        with_loc_load_with_δ_ed_Q_gens_NL_para_Idxs,
        with_loc_load_with_δ_ed_P_load_NL_para_Idxs,
        with_loc_load_with_δ_ed_Q_load_NL_para_Idxs,
        with_loc_load_δ_ed_eq_pf_NL_para_Idxs,
        with_loc_load_with_δ_ed_P_g_loc_load_NL_para_Idxs,
        with_loc_load_with_δ_ed_Q_g_loc_load_NL_para_Idxs)

    #----------------------------------------


    dim_with_loc_load_no_δ_ed_NL_pf_para =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_load,
          dim_Q_load, 
          dim_P_g_loc_load,
          dim_Q_g_loc_load  ]

    _,_, with_loc_load_no_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_with_loc_load_no_δ_ed_NL_pf_para )

    with_loc_load_no_δ_ed_P_gens_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[1]

    with_loc_load_no_δ_ed_Q_gens_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[2]

    with_loc_load_no_δ_ed_P_load_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[3]

    with_loc_load_no_δ_ed_Q_load_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[4]

    with_loc_load_no_δ_ed_P_g_loc_load_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[5]

    with_loc_load_no_δ_ed_Q_g_loc_load_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[6]

    with_loc_load_no_δ_ed_NL_pf_para_Idxs = (
        ; with_loc_load_no_δ_ed_P_gens_NL_para_Idxs,
        with_loc_load_no_δ_ed_Q_gens_NL_para_Idxs,
        with_loc_load_no_δ_ed_P_load_NL_para_Idxs,
        with_loc_load_no_δ_ed_Q_load_NL_para_Idxs,
        with_loc_load_no_δ_ed_P_g_loc_load_NL_para_Idxs,
        with_loc_load_no_δ_ed_Q_g_loc_load_NL_para_Idxs)

    #----------------------------------------


    dim_no_loc_load_with_δ_ed_NL_pf_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_load,
          dim_Q_load,
          dim_δ_ed_eq_pf ]

    _,_, no_loc_load_with_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_no_loc_load_with_δ_ed_NL_pf_para  )

    no_loc_load_with_δ_ed_P_gens_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[1]

    no_loc_load_with_δ_ed_Q_gens_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[2]

    no_loc_load_with_δ_ed_P_load_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[3]

    no_loc_load_with_δ_ed_Q_load_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[4]

    no_loc_load_with_δ_ed_δ_ed_eq_pf_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[5]

    no_loc_load_with_δ_ed_NL_pf_para_Idxs = (
        ; no_loc_load_with_δ_ed_P_gens_NL_para_Idxs,
        no_loc_load_with_δ_ed_Q_gens_NL_para_Idxs,
        no_loc_load_with_δ_ed_P_load_NL_para_Idxs,
        no_loc_load_with_δ_ed_Q_load_NL_para_Idxs,
        no_loc_load_with_δ_ed_δ_ed_eq_pf_NL_para_Idxs )

    
    #----------------------------------------


    dim_no_loc_load_no_δ_ed_NL_pf_para =
        [ dim_P_gens, dim_Q_gens, dim_P_load,
          dim_Q_load  ]

    _,_, no_loc_load_no_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_no_loc_load_no_δ_ed_NL_pf_para )

    no_loc_load_no_δ_ed_P_gens_NL_para_Idxs =
        no_loc_load_no_δ_ed_CA_pf_para_Idxs[1]

    no_loc_load_no_δ_ed_Q_gens_NL_para_Idxs =
        no_loc_load_no_δ_ed_CA_pf_para_Idxs[2]

    no_loc_load_no_δ_ed_P_load_NL_para_Idxs =
        no_loc_load_no_δ_ed_CA_pf_para_Idxs[3]

    no_loc_load_no_δ_ed_Q_load_NL_para_Idxs =
        no_loc_load_no_δ_ed_CA_pf_para_Idxs[4]

    no_loc_load_no_δ_ed_NL_pf_para_Idxs = (
        ; no_loc_load_no_δ_ed_P_gens_NL_para_Idxs,
        no_loc_load_no_δ_ed_Q_gens_NL_para_Idxs,
        no_loc_load_no_δ_ed_P_load_NL_para_Idxs,
        no_loc_load_no_δ_ed_Q_load_NL_para_Idxs)
    
    #------------------------------------------

    dims_δ_ω_ed_dash_eq = length.( get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes))

    _,_, δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq =
        create_size_offset_Idx( dims_δ_ω_ed_dash_eq )

    #------------------------------------------
    
    dicts_net_to_streamlined_idx =
        get_dict_net_to_streamlined_idx(netd)
    
    #------------------------------------------

    NL_pf_para_Idxs = (
        ;  with_loc_load_with_δ_ed_NL_pf_para_Idxs,
        with_loc_load_no_δ_ed_NL_pf_para_Idxs,
        no_loc_load_with_δ_ed_NL_pf_para_Idxs,
        no_loc_load_no_δ_ed_NL_pf_para_Idxs,
        δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq )
    
    #------------------------------------------
    
    net_comp_type_idx =
        (; slack_gens_nodes_idx,
           non_slack_gens_nodes_idx,
           gens_nodes_idx,
           load_nodes_idx,
           transmission_nodes_idx,
           non_gens_nodes_idx,
           load_trans_nodes_idx,
           gens_with_loc_load_idx,
         all_nodes_idx,
         loc_load_exist,
         NL_pf_para_Idxs,
         dicts_net_to_streamlined_idx )
    
    #------------------------------------------
    
    vec_Idx =
        (; vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash,
           vec_Idx_gens_nodes_ω_ed_dash_eq_dash,
           vec_Idx_ra_Xd_Xq_Xd_dash_Xq_dash,
           vec_Idx_ra_Xd_dash_Xq_dash,
           vec_Idx_ra_Xd_Xq,
           vec_Idx_P_Q_nodes,
           vec_Idx_P_Q_gens,
           vec_Idx_P_Q_non_gens,
         vec_Idx_P_Q_gens_loc_load,
         loc_load_exist,
         NL_pf_para_Idxs )
    
    #------------------------------------------
    

    return (; hybrid_pf_etc_idx_and_Idx,
            im_model_pf_idx_and_Idx,
            industrial_model_pf_idx_and_Idx,
            only_gens_industrial_model_pf_idx_and_Idx,
            net_comp_type_idx,
            vec_Idx,
            loc_load_exist,
            NL_pf_para_Idxs,
            dicts_net_to_streamlined_idx)        
end

"""

pf_and_dyn_idx_and_Idx = get_pf_and_dyn_idx_and_Idx(netd)

"""

#------------------------------------------------
#------------------------------------------------
# Industrial model δ_ω_ed_dash_eq_das
#------------------------------------------------
#------------------------------------------------


function get_industrial_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(
    plant)

    return get_some_state_vars_idxs_in_a_plant(
        plant; some_state_vars =
            [:δ, :ω, :ed_dash, :eq_dash] )

end


function get_industrial_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( nodes_collection )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [
            get_industrial_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(
                comp)
            for comp in comp_collection_values
                if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [
            get_industrial_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(
                comp)
            for comp in nodes_collection
                if comp.Bus_type == :Generator ]
    else
        
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [
            get_industrial_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(
                comp)
            for comp in comp_collection_values
                if comp.Bus_type == :Generator ]        
    end

end


function get_industrial_ω_ed_dash_eq_dash_Idxs_in_a_plant(
    plant)

    return get_some_state_vars_idxs_in_a_plant(
        plant; some_state_vars =
            [:ω, :ed_dash, :eq_dash] )

end


function get_industrial_gens_nodes_ω_ed_dash_eq_dash_Idxs_in_plants( nodes_collection )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [
            get_industrial_ω_ed_dash_eq_dash_Idxs_in_a_plant(
                comp)
            for comp in comp_collection_values
                if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [
            get_industrial_ω_ed_dash_eq_dash_Idxs_in_a_plant(
                comp)
            for comp in nodes_collection
                if comp.Bus_type == :Generator ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        return [
            get_industrial_ω_ed_dash_eq_dash_Idxs_in_a_plant(
                comp)
            for comp in comp_collection_values
                if comp.Bus_type == :Generator ]
        
    end

end

#------------------------------------------------

function get_industrial_ω_ed_dash_eq_dash_Idx(
    netd; no_control_device = false )
    

    if no_control_device == false
        
        network_vars_labels =
            generate_industrial_model_sym(
                ; nodes = netd.nodes )

        return [
            get_a_node_states_indices_in_system(
                ; network_vars_labels = network_vars_labels,
                bus_name = a_node.name, vars =
                    [:ω, :ed_dash, :eq_dash])
            for a_node in collect( values(netd.nodes ))
                if a_node.Bus_type == :Generator ]
    else
        
        network_vars_labels = generate_industrial_model_sym(
            ; nodes = netd.nodes, no_control_device = true )

        return [
            get_a_node_states_indices_in_system(
                ; network_vars_labels = network_vars_labels,
                bus_name = a_node.name, vars =
                    [:ω, :ed_dash, :eq_dash])
            for a_node in collect( values(netd.nodes ))
                if a_node.Bus_type == :Generator ]        
    end
        
end 

#------------------------------------------------


function get_industrial_δ_ed_dash_eq_dash_Idx(
    netd; no_control_device = false )

    if no_control_device == false
        
        network_vars_labels =
            generate_industrial_model_sym(
                ; nodes = netd.nodes )

        return [
            get_a_node_δ_ed_dash_eq_dash_indices(
                ; network_vars_labels = network_vars_labels,
                bus_name = a_node.name )
            for a_node in collect( values(netd.nodes ))
                if a_node.Bus_type == :Generator ]
    else
        
        network_vars_labels =
            generate_industrial_model_sym(
                ; nodes = netd.nodes,
                no_control_device = true )

        return [
            get_a_node_δ_ed_dash_eq_dash_indices(
                ; network_vars_labels = network_vars_labels,
                bus_name = a_node.name )
            for a_node in collect( values(netd.nodes ))
                if a_node.Bus_type == :Generator ]        
    end
    
    
end 


function get_industrial_δ_ω_ed_dash_eq_dash_Idx(
    netd; no_control_device = false )

    if no_control_device == false
        network_vars_labels =
            generate_industrial_model_sym(
                ; nodes = netd.nodes )

        return [
            get_a_node_δ_ω_ed_dash_eq_dash_indices(
                ; network_vars_labels = network_vars_labels,
                bus_name = a_node.name )
            for a_node in collect( values(netd.nodes ))
                if a_node.Bus_type == :Generator ]
    else
        
        network_vars_labels =
            generate_industrial_model_sym(
                ; nodes = netd.nodes,
                no_control_device = true )

        return [
            get_a_node_δ_ω_ed_dash_eq_dash_indices(
                ; network_vars_labels = network_vars_labels,
                bus_name = a_node.name )
            for a_node in collect( values(netd.nodes ))
                if a_node.Bus_type == :Generator ]

    end
    
    
end 

function get_industrial_nodes_idx_and_δ_ed_dash_eq_dash_Idx(
    netd )

    network_vars_labels =
        generate_industrial_model_sym(
            ; nodes = netd.nodes )

    return [
        (a_node.Bus_num, get_a_node_δ_ed_dash_eq_dash_indices(
            ; network_vars_labels = network_vars_labels,
            bus_name = a_node.name ))
        for a_node in collect( values(netd.nodes ))
            if a_node.Bus_type == :Generator ]
    
end 


function get_industrial_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(
    netd )

    network_vars_labels = generate_industrial_model_sym(
        ; nodes = netd.nodes )

    return [
        (a_node.Bus_num, get_a_node_δ_ω_ed_dash_eq_dash_indices(
            ; network_vars_labels = network_vars_labels,
            bus_name = a_node.name ))
        for a_node in collect( values(netd.nodes ))
            if a_node.Bus_type == :Generator ]
    
end 



function get_industrial_nodes_idx_and_δ_ed_dash_eq_dash(
    state, netd )

    network_vars_labels = generate_industrial_model_sym(
        ; nodes = netd.nodes )

    indices = [
        ( a_node.Bus_num, get_a_node_δ_ed_dash_eq_dash_indices(
            ; network_vars_labels = network_vars_labels,
            bus_name = a_node.name ))
        for a_node in collect(values(netd.nodes))
            if a_node.Bus_type == :Generator ]
    
    return  [
        (node_idx, state[ idx ] )
        for (node_idx, idx) in indices  ] 
    
end


function get_industrial_nodes_idx_and_δ_ω_ed_dash_eq_dash(
    state, netd )

    network_vars_labels =
        generate_industrial_model_sym(
            ; nodes = netd.nodes )

    indices = [
        ( a_node.Bus_num,
          get_a_node_δ_ω_ed_dash_eq_dash_indices(
              ; network_vars_labels = network_vars_labels,
              bus_name = a_node.name )  )
        for a_node in collect(values(netd.nodes))
            if a_node.Bus_type == :Generator ]
    
    return  [
        (node_idx, state[ idx ] )
        for (node_idx, idx) in indices  ] 
    
end


function get_industrial_nodes_idx_and_δ_ed_dash_eq_dash_view(
    state, netd )

    comp_collection = collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ed_dash_eq_dash_Idx =
        get_industrial_nodes_idx_and_δ_ed_dash_eq_dash_Idx(
            netd )
    
    return [
        ( node_idx, view(state , idx ))
        for (node_idx, idx) in
            nodes_idx_and_δ_ed_dash_eq_dash_Idx ]
    
end


function get_industrial_nodes_idx_and_δ_ed_dash_eq_dash(
    state, netd )

    comp_collection = collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ed_dash_eq_dash_Idx =
        get_industrial_nodes_idx_and_δ_ed_dash_eq_dash_Idx(
            netd )
    
    return [
        ( node_idx, state[idx] )
        for (node_idx, idx) in
            nodes_idx_and_δ_ed_dash_eq_dash_Idx ]
    
end



function get_industrial_nodes_idx_and_δ_ω_ed_dash_eq_dash_view(
    state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
        get_industrial_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    
    return [
        ( a_node_idx, view(state , idx ))
        for ( a_node_idx, idx ) in
            nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx ]
    
end


function get_industrial_nodes_idx_and_δ_ω_ed_dash_eq_dash(
    state, netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
        get_industrial_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    
    return [
        ( a_node_idx, state[idx] )
        for ( a_node_idx, idx ) in
            nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx ]
    
end

#------------------------------------------------

#------------------------------------------------
#------------------------------------------------
# im model δ_ω_ed_dash_eq_das
#------------------------------------------------
#------------------------------------------------


function get_im_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(plant)

    return get_some_state_vars_idxs_in_a_plant(plant; some_state_vars = [:δ, :ω, :ed_dash, :eq_dash] )

end


function get_gens_nodes_im_δ_ω_ed_dash_eq_dash_Idxs_in_plants( nodes_collection )
    
    if isa(nodes_collection, OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        return [
            get_im_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
            for comp in comp_collection_values
                if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [
            get_im_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
            for comp in nodes_collection
                if comp.Bus_type == :Generator ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        return [
            get_im_δ_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
            for comp in comp_collection_values
                if comp.Bus_type == :Generator ]
        
    end

end


function get_im_ω_ed_dash_eq_dash_Idxs_in_a_plant(
    plant)

    return get_some_state_vars_idxs_in_a_plant(
        plant;
        some_state_vars =
            [:ω, :ed_dash, :eq_dash] )

end


function get_gens_nodes_im_ω_ed_dash_eq_dash_Idxs_in_plants(
    nodes_collection )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        return [
            get_im_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
            for comp in comp_collection_values
                if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [
            get_im_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
            for comp in nodes_collection
                if comp.Bus_type == :Generator ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [
            get_im_ω_ed_dash_eq_dash_Idxs_in_a_plant(comp)
            for comp in comp_collection_values
                if comp.Bus_type == :Generator ]
        
    end

end


#------------------------------------------------



function get_im_ω_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        generate_im_sym(; nodes = netd.nodes )

    return [
        get_a_node_states_indices_in_system(
            ; network_vars_labels = network_vars_labels,
            bus_name = a_node.name,
            vars = [:ω, :ed_dash, :eq_dash])
        for a_node in collect( values(netd.nodes ))
            if a_node.Bus_type == :Generator ]
    
end 

#------------------------------------------------


function get_im_δ_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        generate_im_sym(;
                        nodes = netd.nodes )

    return [
        get_a_node_δ_ed_dash_eq_dash_indices(
            ; network_vars_labels = network_vars_labels,
            bus_name = a_node.name )
        for a_node in collect( values(netd.nodes ))
            if a_node.Bus_type == :Generator ]
    
end 


function get_im_δ_ω_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        generate_im_sym(; nodes = netd.nodes )

    return [
        get_a_node_δ_ω_ed_dash_eq_dash_indices(
            ; network_vars_labels =
                network_vars_labels,
            bus_name = a_node.name )
        for a_node in collect( values(netd.nodes ))
            if a_node.Bus_type == :Generator ]
    
end 

function get_im_nodes_idx_and_δ_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        generate_im_sym(;
                        nodes = netd.nodes )

    return [
        (a_node.Bus_num,
         get_a_node_δ_ed_dash_eq_dash_indices(
             ; network_vars_labels =
                 network_vars_labels,
             bus_name = a_node.name ))
        for a_node in collect( values(netd.nodes ))
            if a_node.Bus_type == :Generator ]
    
end 


function get_im_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    network_vars_labels =
        generate_im_sym(;
                        nodes = netd.nodes )

    return [
        (a_node.Bus_num,
         get_a_node_δ_ω_ed_dash_eq_dash_indices(
             ; network_vars_labels =
                 network_vars_labels,
             bus_name = a_node.name ))
        for a_node in collect( values(netd.nodes ))
            if a_node.Bus_type == :Generator ]
    
end 



function get_im_nodes_idx_and_δ_ed_dash_eq_dash(
    state,
    netd )

    network_vars_labels =
        generate_im_sym(
            ; nodes = netd.nodes )

    indices = [
        ( a_node.Bus_num,
          get_a_node_δ_ed_dash_eq_dash_indices(
              ; network_vars_labels =
                  network_vars_labels,
              bus_name = a_node.name ))
        for a_node in collect(values(netd.nodes))
            if a_node.Bus_type == :Generator ]
    
    return  [
        (node_idx, state[ idx ] )
        for (node_idx, idx) in indices  ] 
    
end


function get_im_nodes_idx_and_δ_ω_ed_dash_eq_dash(
    state,
    netd )

    network_vars_labels =
        generate_im_sym(;
                        nodes = netd.nodes )

    indices = [
        ( a_node.Bus_num,
          get_a_node_δ_ω_ed_dash_eq_dash_indices(
              ; network_vars_labels =
                  network_vars_labels,
              bus_name = a_node.name )  )
        for a_node in collect(values(netd.nodes))
            if a_node.Bus_type == :Generator ]
    
    return  [
        (node_idx, state[ idx ] )
        for (node_idx, idx) in indices  ] 
    
end


function get_im_nodes_idx_and_δ_ed_dash_eq_dash_view(
    state,
    netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ed_dash_eq_dash_Idx =
        get_im_nodes_idx_and_δ_ed_dash_eq_dash_Idx(netd )

    
    return [
        ( node_idx, view(state , idx ))
        for (node_idx, idx) in
            nodes_idx_and_δ_ed_dash_eq_dash_Idx   ]
    
end



function get_im_nodes_idx_and_δ_ω_ed_dash_eq_dash_view(
    state,
    netd )

    comp_collection =
        collect(values(netd.nodes) )
    
    nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx =
        get_im_nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx(netd )

    
    return [
        ( a_node_idx, view(state , idx ))
        for (a_node_idx, idx) in
            nodes_idx_and_δ_ω_ed_dash_eq_dash_Idx ]
    
end

#------------------------------------------------


#------------------------------------------------
#------------------------------------------------


function create_gen_nodes_state_views(
    state, gens_nodes_collection )

    _, _,  nodes_state_Idx = create_flattened_type_states_dims_offset_Idx_for_gens_nodes( gens_nodes_collection )

    # gens_nodes_state = [ state[a_node_state_idx] for a_node_state_idx in nodes_state_Idx ]
    
    return [ view(state, a_node_state_idx) for a_node_state_idx in nodes_state_Idx ]
    
end



#------------------------------------------------
#------------------------------------------------



function get_slack_gens_nodes(comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        comp_collection_values = collect(values(comp_collection))

        return [  comp for comp in comp_collection_values if comp.Bus_type == :Generator &&  comp.isa_slack == true ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [ comp  for comp in comp_collection if comp.Bus_type == :Generator && comp.isa_slack == true ]
    else
        comp_collection_values = collect(values(comp_collection))
        return  [ comp for comp in comp_collection_values if comp.Bus_type == :Generator && comp.isa_slack == true ]
        
    end
end


function get_gens_nodes( nodes_collection  )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(values(nodes_collection))

        return [ comp for comp in comp_collection_values if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp for comp in nodes_collection if comp.Bus_type == :Generator ]
    else
        comp_collection_values = collect(values(nodes_collection))
        return [ comp for comp in comp_collection_values if comp.Bus_type == :Generator ]
        
    end
    
end

function get_gens_nodes_with_loc_loads( nodes_collection  )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(values(nodes_collection))

        return [ comp
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator &&
                         comp.with_loc_load == true ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp
                 for comp in nodes_collection
                     if comp.Bus_type == :Generator &&
                         comp.with_loc_load == true ]
    else
        
        comp_collection_values = collect(values(nodes_collection))
        
        return [ comp
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator &&
                         comp.with_loc_load == true ]
        
    end
    
end



function get_non_slack_gens_nodes(comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        comp_collection_values = collect(values(comp_collection))

        return [  comp for comp in comp_collection_values if comp.Bus_type == :Generator &&  comp.isa_slack == false ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [ comp  for comp in nodes_collection if comp.Bus_type == :Generator && comp.isa_slack == false ]
    else
        comp_collection_values = collect(values(comp_collection))
        return  [ comp for comp in comp_collection_values if comp.Bus_type == :Generator && comp.isa_slack == false ]
        
    end
end


function get_non_gens_nodes( nodes_collection  )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(values(nodes_collection))

        return [ comp for comp in comp_collection_values
                    if comp.Bus_type != :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp for comp in nodes_collection
                    if comp.Bus_type != :Generator ]
    else
        comp_collection_values = collect(values(nodes_collection))
        return [ comp for comp in comp_collection_values
                    if comp.Bus_type != :Generator ]
        
    end
    
end


function get_load_nodes( nodes_collection  )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(values(nodes_collection))

        return [ comp for comp in comp_collection_values
                    if comp.Bus_type == :Load ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp for comp in nodes_collection
                    if comp.Bus_type == :Load ]
    else
        comp_collection_values = collect(values(nodes_collection))
        return [ comp for comp in comp_collection_values
                    if comp.Bus_type == :Load ]
        
    end
    
end


function get_transmission_nodes( nodes_collection  )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(
            values(nodes_collection))

        return [ comp for comp in comp_collection_values
                    if comp.Bus_type == :Transmission ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp for comp in nodes_collection
                    if comp.Bus_type == :Transmission ]
    else
        comp_collection_values = collect(
            values(nodes_collection))
        return [ comp for comp in comp_collection_values
                    if comp.Bus_type == :Transmission ]
        
    end
    
end



function get_all_nodes( nodes_collection  )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(
            values(nodes_collection))

        return [ comp for comp in comp_collection_values
                     ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp for comp in nodes_collection
                     ]
    else
        comp_collection_values = collect(
            values(nodes_collection))
        return [ comp for comp in comp_collection_values
                     ]
        
    end
    
end

#------------------------------------------------

function get_slack_gens_nodes_idx(
    nodes_collection  )

    slack_gens_nodes = get_slack_gens_nodes(
        nodes_collection  )
        
    return [a_node.Bus_num
            for a_node in slack_gens_nodes ]
end


function get_gens_nodes_idx( nodes_collection  )

    gens_nodes = get_gens_nodes( nodes_collection  )
        
    return [a_node.Bus_num for a_node in gens_nodes ]
end


function get_gens_nodes_with_loc_loads_idx(
    nodes_collection  )

    gens_nodes_with_loc_loads =
        get_gens_nodes_with_loc_loads(
            nodes_collection  )
        
    return [a_node.Bus_num
            for a_node in
                gens_nodes_with_loc_loads ]
end


function get_non_slack_gens_nodes_idx(comp_collection )

    non_slack_gens =
        get_non_slack_gens_nodes(comp_collection )

    return [a_node.Bus_num for a_node in non_slack_gens  ]

end


function get_load_nodes_idx(
    nodes_collection  )

    load_nodes =
        get_load_nodes( nodes_collection  )
        
    return [a_node.Bus_num
            for a_node in load_nodes ]
end




function get_transmission_nodes_idx(
    nodes_collection  )

    transmission_nodes =
        get_transmission_nodes( nodes_collection  )
        
    return [a_node.Bus_num
            for a_node in transmission_nodes ]
end


function get_non_gens_nodes_idx(
    nodes_collection  )

    non_gens_nodes =
        get_non_gens_nodes( nodes_collection  )
    return [a_node.Bus_num
            for a_node in non_gens_nodes ]
end



function get_all_nodes_idx(
    nodes_collection  )

    all_nodes =
        get_all_nodes( nodes_collection  )
    
    return [a_node.Bus_num for a_node in all_nodes ]
end

#------------------------------------------------


function get_a_gen_post_pf_data(
    gen, bus_dict_init )

    power_flow_data =
        bus_dict_init[gen.name]

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
     
    return ( vh, θh, ph, qh, i_r,
             i_i, pg, qg, ig_r, ig_i )
    
end


function get_gen_post_pf_data(
    gen, power_flow_data )

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
     
    return ( vh, θh, ph, qh,
             i_r, i_i, pg, qg, ig_r, ig_i )
    
end


function get_gens_post_pf_data(
    nodes_collection, bus_dict_init )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        return [ comp.Bus_type == :Generator ?
            bus_dict_init[comp.name] : []
                 for comp in comp_collection_values ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp.Bus_type == :Generator ?
            bus_dict_init[comp.name] : []
                 for comp in nodes_collection ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ comp.Bus_type == :Generator ?
            bus_dict_init[comp.name] : []
                 for comp in comp_collection_values  ]
        
    end
    
end


function get_only_gens_post_pf_data(
    nodes_collection, bus_dict_init )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [  bus_dict_init[comp.name]
                  for comp in comp_collection_values
                      if comp.Bus_type == :Generator  ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [  bus_dict_init[comp.name]
                  for comp in nodes_collection
                      if comp.Bus_type == :Generator   ]
    else
        
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [  bus_dict_init[comp.name]
                  for comp in comp_collection_values
                      if comp.Bus_type == :Generator  ]
        
    end
    
end



function get_gens_vh_θh_post_pf(
    nodes_collection, bus_dict_init )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_vh_θh_post_pf(
            bus_dict_init[comp.name])
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_vh_θh_post_pf(
            bus_dict_init[comp.name])
                 for comp in nodes_collection
                     if comp.Bus_type == :Generator ]
    else
        
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ get_vh_θh_post_pf(
            bus_dict_init[comp.name])
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    end
    
end


function get_gens_ur_ui_post_pf(
    nodes_collection, bus_dict_init )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_ur_ui_post_pf(
            bus_dict_init[comp.name])
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_ur_ui_post_pf(
            bus_dict_init[comp.name])
                 for comp in nodes_collection
                     if comp.Bus_type == :Generator ]
        
    else
        
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ get_ur_ui_post_pf(
            bus_dict_init[comp.name])
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    end
    
end



function get_a_gen_vh_θh_Pg_Qg_from_post_pf(
    gen, bus_dict_init  )

    power_flow_data = bus_dict_init[gen.name]

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
     
    return ( vh, θh, pg, qg )
    
end



function get_gen_vh_θh_Pg_Qg_from_post_pf(
    gen, power_flow_data )

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
     
    return ( vh, θh, pg, qg )
    
end



function get_gens_vh_θh_Pg_Qg_from_post_pf_data(
    nodes_collection, bus_dict_init )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(
            values(nodes_collection))

        return [ comp.Bus_type == :Generator ?
            get_gen_vh_θh_Pg_Qg_from_post_pf(
                comp.Gen,  bus_dict_init[comp.name] ) : []
                 for comp in comp_collection_values ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp.Bus_type == :Generator ?
            get_gen_vh_θh_Pg_Qg_from_post_pf(
                comp.Gen,  bus_dict_init[comp.name] )  : []
                 for comp in nodes_collection ]
    else
        comp_collection_values = collect(
            values(nodes_collection))
        return [ comp.Bus_type == :Generator ?
            get_gen_vh_θh_Pg_Qg_from_post_pf(
                comp.Gen,  bus_dict_init[comp.name] )  : []
                 for comp in comp_collection_values  ]
        
    end
    
end


function get_only_gens_vh_θh_Pg_Qg_from_post_pf_data(
    nodes_collection, bus_dict_init )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(
            values(nodes_collection))

        return [  get_gen_vh_θh_Pg_Qg_from_post_pf(
            comp.Gen,  bus_dict_init[comp.name] )
                  for comp in comp_collection_values
                      if comp.Bus_type == :Generator  ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [  get_gen_vh_θh_Pg_Qg_from_post_pf(
            comp.Gen,  bus_dict_init[comp.name] )
                  for comp in nodes_collection
                      if comp.Bus_type == :Generator ]
    else
        comp_collection_values = collect(
            values(nodes_collection))
        return [  get_gen_vh_θh_Pg_Qg_from_post_pf(
            comp.Gen,  bus_dict_init[comp.name] )
                  for comp in comp_collection_values
                      if comp.Bus_type == :Generator  ]
        
    end
    
end

#------------------------------------------------


function get_non_gen_post_pf_data(
    non_gen, power_flow_data )

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
     
    return ( vh, θh, ph, qh, i_r,
             i_i, pg, qg, ig_r, ig_i )
    
end


function get_non_gens_post_pf_data(
    nodes_collection, bus_dict_init )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(
            values(nodes_collection))

        return [ comp.Bus_type != :Generator ?
            bus_dict_init[comp.name] : []
                 for comp in comp_collection_values ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp.Bus_type != :Generator ?
            bus_dict_init[comp.name] : []
                 for comp in nodes_collection ]
    else
        comp_collection_values = collect(
            values(nodes_collection))
        return [ comp.Bus_type != :Generator ?
            bus_dict_init[comp.name] : []
                 for comp in comp_collection_values  ]
        
    end
    
end


function get_only_non_gens_post_pf_data(
    nodes_collection, bus_dict_init )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(
            values(nodes_collection))

        return [ bus_dict_init[comp.name]
                 for comp in comp_collection_values if
                     comp.Bus_type != :Generator  ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [   bus_dict_init[comp.name]
                   for comp in nodes_collection if
                       comp.Bus_type != :Generator  ]
    else
        comp_collection_values = collect(
            values(nodes_collection))
        return [  bus_dict_init[comp.name]
                  for comp in comp_collection_values if
                      comp.Bus_type != :Generator  ]
        
    end
    
end



function get_a_non_gen_vh_θh_Pg_Qg_from_post_pf(
    non_gen, bus_dict_init  )

    power_flow_data = bus_dict_init[non_gen.name]
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
     
    return ( vh, θh, ph, qh )
    
end


function get_non_gen_vh_θh_Pg_Qg_from_post_pf(
    non_gen, power_flow_data )

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
     
    return ( vh, θh, ph, qh )
    
end


function get_non_gens_vh_θh_Pg_Qg_from_post_pf_data(
    nodes_collection, bus_dict_init )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values = collect(
            values(nodes_collection))

        return [ comp.Bus_type == :Load ?
            get_non_gen_vh_θh_Pg_Qg_from_post_pf(
                comp.Load,
                bus_dict_init[comp.name]) :
                    comp.Bus_type == :Transmission ?
                    get_non_gen_vh_θh_Pg_Qg_from_post_pf(
                        comp.Trans,
                        bus_dict_init[comp.name])  : []
                 for comp in comp_collection_values
                     if comp.Bus_type != :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ comp.Bus_type == :Load ?
            get_non_gen_vh_θh_Pg_Qg_from_post_pf(
                comp.Load,
                bus_dict_init[comp.name]) :
                    comp.Bus_type == :Transmission ?
                    get_non_gen_vh_θh_Pg_Qg_from_post_pf(
                        comp.Trans,
                        bus_dict_init[comp.name]) : []
                 for comp in nodes_collection
                     if comp.Bus_type != :Generator ]
        
    else
        comp_collection_values = collect(
            values(nodes_collection))
        return [ comp.Bus_type == :Load ?
            get_non_gen_vh_θh_Pg_Qg_from_post_pf(
                comp.Load,
                bus_dict_init[comp.name] ) :
                    comp.Bus_type == :Transmission ?
                    get_non_gen_vh_θh_Pg_Qg_from_post_pf(
                        comp.Trans,
                        bus_dict_init[comp.name]) : []
                 for comp in comp_collection_values if
                     comp.Bus_type != :Generator ]
        
    end
    
end


#------------------------------------------------

#------------------------------------------------


function get_gens_δ_id_iq_vd_vq_ed_dash_eq_dash_etc_view_from_pf( netd, bus_dict_init )

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash,
                :X_q_dash ] )

    gens_vh_θh_Pg_Qg =
        get_gens_vh_θh_Pg_Qg_from_post_pf_data(
            netd.nodes, bus_dict_init )

    return [vh_θh_Pg_Qg != [] ? get_component_δ_id_iq_vd_vq_ed_dash_eq_dash_from_pf(vh_θh_Pg_Qg..., ra_Xd_Xq_Xd_dash_Xq_dash...) : [] for (vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash) in zip(gens_vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash_view )]

end


function get_gens_δ_id_iq_vd_vq_ed_dash_eq_dash_etc_view_from_pf( nodes_param_values, nodes_collection, bus_dict_init )

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            nodes_param_values,
            nodes_collection; param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash,
                :X_q_dash ] )

    gens_vh_θh_Pg_Qg =
        get_gens_vh_θh_Pg_Qg_from_post_pf_data(
            nodes_collection, bus_dict_init )

    return [vh_θh_Pg_Qg != [] ? get_component_δ_id_iq_vd_vq_ed_dash_eq_dash_from_pf(vh_θh_Pg_Qg..., ra_Xd_Xq_Xd_dash_Xq_dash...) : [] for (vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash) in zip(gens_vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash_view )]

end


function get_gens_δ_ω_ed_dash_eq_dash_and_τm_etc_view_from_pf( netd, bus_dict_init )

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash,
                :X_q_dash ] )

    gens_vh_θh_Pg_Qg =
        get_gens_vh_θh_Pg_Qg_from_post_pf_data(
            netd.nodes,
            bus_dict_init )

    return [  get_component_δ_ω_ed_dash_eq_dash_from_pf(
        vh_θh_Pg_Qg...,
        ra_Xd_Xq_Xd_dash_Xq_dash...)
              for (vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash) in
                  zip(gens_vh_θh_Pg_Qg,
                      ra_Xd_Xq_Xd_dash_Xq_dash_view )
                  if vh_θh_Pg_Qg != [] ]

end


function get_gens_δ_ω_ed_dash_eq_dash_and_τm_etc_view_from_pf( nodes_param_values, nodes_collection, bus_dict_init )

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            nodes_param_values,
            nodes_collection; param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash,
                :X_q_dash ] )

    gens_vh_θh_Pg_Qg =
        get_gens_vh_θh_Pg_Qg_from_post_pf_data(
            nodes_collection, bus_dict_init )

    return [  get_component_δ_ω_ed_dash_eq_dash_from_pf(
        vh_θh_Pg_Qg...,
        ra_Xd_Xq_Xd_dash_Xq_dash...)
              for (vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash) in
                  zip(gens_vh_θh_Pg_Qg,
                      ra_Xd_Xq_Xd_dash_Xq_dash_view )
                  if vh_θh_Pg_Qg != [] ]

end

#------------------------------------------------

function get_gens_δ_id_iq_vd_vq_ed_dash_eq_dash_etc_from_pf(
    netd, bus_dict_init )

    ra_Xd_Xq_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash,
                :X_q_dash ] )

    gens_vh_θh_Pg_Qg =
        get_gens_vh_θh_Pg_Qg_from_post_pf_data(
            netd.nodes, bus_dict_init )

    return [vh_θh_Pg_Qg != [] ? get_component_δ_id_iq_vd_vq_ed_dash_eq_dash_from_pf(vh_θh_Pg_Qg..., ra_Xd_Xq_Xd_dash_Xq_dash...) : [] for (vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash) in zip(gens_vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash )]

end


function get_gens_δ_id_iq_vd_vq_ed_dash_eq_dash_etc_from_pf(
    nodes_param_values, nodes_collection, bus_dict_init )

    ra_Xd_Xq_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            nodes_param_values,
            nodes_collection; param_list = [
                :ra, :X_d,
                :X_q, :X_d_dash,
                :X_q_dash ] )

    gens_vh_θh_Pg_Qg =
        get_gens_vh_θh_Pg_Qg_from_post_pf_data(
            nodes_collection,
            bus_dict_init )

    return [vh_θh_Pg_Qg != [] ? get_component_δ_id_iq_vd_vq_ed_dash_eq_dash_from_pf(vh_θh_Pg_Qg..., ra_Xd_Xq_Xd_dash_Xq_dash...) : [] for (vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash) in zip(gens_vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash )]

end


function get_gens_δ_ω_ed_dash_eq_dash_and_τm_etc_from_pf(
    netd, bus_dict_init )

    ra_Xd_Xq_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            netd.nodes_param,
            netd.nodes; param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash, :X_q_dash ] )

    gens_vh_θh_Pg_Qg =
        get_gens_vh_θh_Pg_Qg_from_post_pf_data(
            netd.nodes,
            bus_dict_init )

    return [
        get_component_δ_ω_ed_dash_eq_dash_from_pf(
            vh_θh_Pg_Qg...,
            ra_Xd_Xq_Xd_dash_Xq_dash...)
        for (vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash) in
            zip(gens_vh_θh_Pg_Qg,
                ra_Xd_Xq_Xd_dash_Xq_dash )
            if vh_θh_Pg_Qg != [] ]

end


function get_gens_δ_ω_ed_dash_eq_dash_and_τm_etc_from_pf(
    nodes_param_values, nodes_collection, bus_dict_init )

    ra_Xd_Xq_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            nodes_param_values,
            nodes_collection; param_list =
                [ :ra, :X_d, :X_q,
                  :X_d_dash,
                  :X_q_dash ] )

    gens_vh_θh_Pg_Qg =
        get_gens_vh_θh_Pg_Qg_from_post_pf_data(
            nodes_collection, bus_dict_init )

    return [
        get_component_δ_ω_ed_dash_eq_dash_from_pf(
            vh_θh_Pg_Qg...,
            ra_Xd_Xq_Xd_dash_Xq_dash...)
        for (vh_θh_Pg_Qg, ra_Xd_Xq_Xd_dash_Xq_dash) in
            zip(gens_vh_θh_Pg_Qg,
                ra_Xd_Xq_Xd_dash_Xq_dash )
            if vh_θh_Pg_Qg != [] ]

end

#------------------------------------------------


function get_components_post_pf_data(
    comp_collection, bus_dict_init )
    
    if isa(comp_collection, OrderedDict)
        comp_collection_values = collect(
            values(comp_collection))

        return [ bus_dict_init[ comp.name ]
                 for comp in comp_collection_values ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [ bus_dict_init[ comp.name ]
                 for comp in comp_collection ]
    else
        comp_collection_values = collect(
            values(comp_collection))
        return [ bus_dict_init[ comp.name ]
                 for comp in comp_collection_values ]
        
    end
    
end



function get_components_Idx_node_type_and_post_pf_data(
    comp_collection, bus_dict_init )
    
    if isa(comp_collection, OrderedDict)
        comp_collection_values = collect(
            values(comp_collection))

        return [
            ( comp.Bus_num, comp.Bus_type,
              bus_dict_init[ comp.name ] )
            for comp in comp_collection_values ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [
            ( comp.Bus_num, comp.Bus_type,
              bus_dict_init[ comp.name ] )
            for comp in comp_collection ]
    else
        comp_collection_values = collect(
            values(comp_collection))
        return [
            ( comp.Bus_num, comp.Bus_type,
              bus_dict_init[ comp.name ] )
            for comp in comp_collection_values ]
        
    end
    
end


function get_nodes_init_data_from_pf(
    netd, bus_dict_init )

    nodes_post_pf_data = get_components_post_pf_data(
        netd.nodes, bus_dict_init )

    return [
        self_init_state(comp, pf_data)
        for (comp, pf_data )  in
            zip( collect(values( netd.nodes )),
                 nodes_post_pf_data ) ]

end


function get_gens_nodes_init_data_from_pf(
    gens_nodes_collections, bus_dict_init )
   
    gens_nodes_post_pf_data = get_gens_post_pf_data(
        gens_nodes_collections, bus_dict_init )

    return [
        self_init_state(comp, pf_data)
        for (comp, pf_data )  in
            zip( gens_nodes_collections ,
                 gens_nodes_post_pf_data )
            if comp.Bus_type == :Generator]


end



# get_gens_post_pf_data( nodes_collection, bus_dict_init )

function get_non_gens_nodes_init_data_from_pf(
    non_gens_nodes_collections,
    bus_dict_init )
   
    non_gens_nodes_post_pf_data =
        get_non_gens_post_pf_data(
            non_gens_nodes_collections,
            bus_dict_init )

    return [
        self_init_state(comp, pf_data)
        for (comp, pf_data )  in
            zip( non_gens_nodes_collections,
                 non_gens_nodes_post_pf_data )
            if comp.Bus_type != :Generator]


end

#------------------------------------------------


function get_gens_nodes_ωs_τm_v_ref(
    nodes_collection, bus_dict_init )

    if isa(nodes_collection,  OrderedDict)

        comp_collection_values = collect(
            values(nodes_collection))

        gens_nodes_collection = [
            comp
            for comp  in comp_collection_values
                if comp.Bus_type == :Generator ]

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection, bus_dict_init )
   
        return [
            [ ωs, a_node_init.aux[2],
              a_node_init.f_t.avr_f_t[1]]
            for a_node_init in
                gens_nodes_init_data_from_pf ]
        
    elseif  isa(nodes_collection, Union{Array, Vector})

        gens_nodes_collection = [
            comp for comp  in nodes_collection
                if comp.Bus_type == :Generator]

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection,
                bus_dict_init )
   
        return [
            [ ωs, a_node_init.aux[2],
              a_node_init.f_t.avr_f_t[1]]
            for a_node_init in
                gens_nodes_init_data_from_pf ]

    else
        
        comp_collection_values = collect(
            values(nodes_collection))

        gens_nodes_collection = [
            comp for comp  in
                comp_collection_values
                if comp.Bus_type == :Generator]

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection, bus_dict_init )
        
        return [
            [ ωs, a_node_init.aux[2],
              a_node_init.f_t.avr_f_t[1]]
            for a_node_init in
                gens_nodes_init_data_from_pf ]
    end
        
end


function get_gens_nodes_ωs_τm_vref_porder(
    nodes_collection, bus_dict_init )

    if isa(nodes_collection,  OrderedDict)

        comp_collection_values =
            collect(values(nodes_collection))

        gens_nodes_collection = [
            comp for comp  in
                comp_collection_values
                if comp.Bus_type == :Generator ]

        gens_isa_condenser_bool = [
            comp.isa_condenser ?
                true : false
            for comp  in
                gens_nodes_collection ]

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection,
                bus_dict_init )
   
        return [ isa_condenser ?
            [ ωs, a_node_init.aux[2],
              a_node_init.f_t.avr_f_t[1], 0.0 ] :
                  [ ωs, a_node_init.aux[2],
                    a_node_init.f_t.avr_f_t[1],
                    a_node_init.f_t.gov_f_t[1] ]
                 for (a_node_init, isa_condenser ) in
                     zip( gens_nodes_init_data_from_pf,
                          gens_isa_condenser_bool ) ]
        
    elseif  isa(nodes_collection, Union{Array, Vector})

        gens_nodes_collection = [
            comp for comp  in
                nodes_collection
                if comp.Bus_type == :Generator]

        gens_isa_condenser_bool = [
            comp.isa_condenser ? true : false
            for comp  in gens_nodes_collection ]

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection, bus_dict_init )
   
        return [ isa_condenser ?
            [ ωs, a_node_init.aux[2],
              a_node_init.f_t.avr_f_t[1], 0.0 ] :
                  [ ωs, a_node_init.aux[2],
                    a_node_init.f_t.avr_f_t[1],
                    a_node_init.f_t.gov_f_t[1] ]
                 for (a_node_init, isa_condenser ) in
                     zip( gens_nodes_init_data_from_pf,
                          gens_isa_condenser_bool ) ]

    else
        
        comp_collection_values =
            collect(values(nodes_collection))        
        
        gens_nodes_collection = [
            comp for comp  in
                comp_collection_values
                if comp.Bus_type == :Generator]

        gens_isa_condenser_bool = [
            comp.isa_condenser ? true : false
            for comp  in gens_nodes_collection ]

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection, bus_dict_init )
        
   
        return [ isa_condenser ?
            [ ωs, a_node_init.aux[2],
              a_node_init.f_t.avr_f_t[1], 0.0 ] :
                  [ ωs, a_node_init.aux[2],
                    a_node_init.f_t.avr_f_t[1],
                    a_node_init.f_t.gov_f_t[1] ]
                 for (a_node_init, isa_condenser ) in
                     zip( gens_nodes_init_data_from_pf,
                          gens_isa_condenser_bool ) ]
        
    end
        
end


"""

f_t.avr_f_t ≡ [v_ref0]

f_t.gov_f_t ≡ [p_order0, ω_ref0]

"""
function get_gens_nodes_ωs_ωref0_vref0_porder0(
    nodes_collection, bus_dict_init )

    
    if isa(nodes_collection,  OrderedDict)

        comp_collection_values = collect(
            values(nodes_collection))

        gens_nodes_collection = [
            comp
            for comp  in comp_collection_values
                if comp.Bus_type == :Generator ]

        gens_isa_condenser_bool = [
            comp.isa_condenser ? true : false
            for comp  in gens_nodes_collection ]

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
            gens_nodes_collection, bus_dict_init )
   
        return [
            isa_condenser ?
                [ ωs,  0.0, a_node_init.f_t.avr_f_t[1],
                  0.0 ] :
                [ ωs, a_node_init.f_t.gov_f_t[2],
                  a_node_init.f_t.avr_f_t[1],
                  a_node_init.f_t.gov_f_t[1] ]
            for (a_node_init, isa_condenser ) in zip(
                gens_nodes_init_data_from_pf,
                gens_isa_condenser_bool ) ]
        
    elseif  isa(nodes_collection, Union{Array, Vector})

        gens_nodes_collection = [
            comp
            for comp  in nodes_collection
                if comp.Bus_type == :Generator]

        gens_isa_condenser_bool = [
            comp.isa_condenser ?
                true :
                false
            for comp  in gens_nodes_collection ]

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
            gens_nodes_collection, bus_dict_init )
   
        return [
            isa_condenser ?
                [ ωs,  0.0, a_node_init.f_t.avr_f_t[1],
                  0.0 ] :
                [ ωs, a_node_init.f_t.gov_f_t[2],
                  a_node_init.f_t.avr_f_t[1],
                  a_node_init.f_t.gov_f_t[1] ]
            for (a_node_init, isa_condenser ) in zip(
                gens_nodes_init_data_from_pf,
                gens_isa_condenser_bool ) ]

    else
        
        comp_collection_values =
            collect(values(nodes_collection))        
        
        gens_nodes_collection = [
            comp
            for comp  in comp_collection_values
                if comp.Bus_type == :Generator]

        gens_isa_condenser_bool = [
            comp.isa_condenser ?
                true :
                false for comp  in gens_nodes_collection ]

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
            gens_nodes_collection,
            bus_dict_init )        
   
        return [
            isa_condenser ?
                [ ωs,  0.0, a_node_init.f_t.avr_f_t[1],
                  0.0 ] :
                      [ ωs, a_node_init.f_t.gov_f_t[2],
                        a_node_init.f_t.avr_f_t[1],
                        a_node_init.f_t.gov_f_t[1] ]
            for (a_node_init, isa_condenser ) in zip(
                gens_nodes_init_data_from_pf,
                gens_isa_condenser_bool ) ]        
    end
        
end


function get_gens_τm_vf(
    nodes_collection, bus_dict_init )

    if isa(nodes_collection,  OrderedDict)

        comp_collection_values =
            collect(values(nodes_collection))

        gens_nodes_collection = [
            comp for comp  in
                comp_collection_values
                if comp.Bus_type == :Generator ]

        gen_plants_init_state = [
            self_init_state(
                a_comp, bus_dict_init[a_comp.name])
            for a_comp in  gens_nodes_collection ]
        
        gen_plants_init_aux = third.(
            gen_plants_init_state ) # [vf0, pm0]
        
        return [ [ gen_init_aux[2], gen_init_aux[1] ]
                 for gen_init_aux in
                     gen_plants_init_aux ]
        
    elseif  isa(nodes_collection, Union{Array, Vector})

        gens_nodes_collection = [
            comp for comp  in nodes_collection
                if comp.Bus_type == :Generator ]
        
        gen_plants_init_state = [
            self_init_state(
                a_comp, bus_dict_init[a_comp.name] )
            for a_comp in  gens_nodes_collection ]
        
        gen_plants_init_aux = third.(
            gen_plants_init_state ) # [vf0, pm0]
   
        return [
            [ gen_init_aux[2], gen_init_aux[1] ]
            for gen_init_aux in gen_plants_init_aux ]

    else
        
        comp_collection_values = collect(
            values(nodes_collection))        
        
        gens_nodes_collection = [
            comp for comp  in
                comp_collection_values
                if comp.Bus_type == :Generator]
        
        gen_plants_init_state = [
            self_init_state(
                a_comp, bus_dict_init[a_comp.name] )
            for a_comp in  gens_nodes_collection ]
        
        gen_plants_init_aux = third.(
            gen_plants_init_state ) # [vf0, pm0]
   
        return [ [ gen_init_aux[2], gen_init_aux[1] ]
                 for gen_init_aux in
                     gen_plants_init_aux ]
        
    end
        
end

# ------------------------------------------------------
# ------------------------------------------------------
# plants states init, f_t, aux
# ------------------------------------------------------
# ------------------------------------------------------

function get_net_nodes_x0_f_t_aux(
    Dyn_Nodes, pf_init_dict;
    key="bus_dict_init")

    x0_f_t_aux = []

    for (plant,  power_flow_data) in
        zip( collect(values(Dyn_Nodes)),
             collect( values(pf_init_dict[key])))
        
        push!(x0_f_t_aux,
              self_init_state(plant, power_flow_data))
    end

    return x0_f_t_aux

    # vf0p0 = last.(map((plant,  power_flow_data) -> self_init_state(plant,  power_flow_data), zip(collect(values(ieee_14_Dyn_Nodes)), collect(values(pf_init_dict[ "bus_dict_init" ])))))
    
end


function get_net_edges_x0_f_t(
    Dyn_Edges, pf_init_dict;
    key="branch_dict_init")

    x0_f_t = []

    for (plant,  power_flow_data) in
        zip( collect(values(Dyn_Edges)),
             collect( values(pf_init_dict[key])))
        
        push!(x0_f_t,
              self_init_state(plant, power_flow_data))
    end

    return x0_f_t

    # vf0p0 = last.(map((plant,  power_flow_data) -> self_init_state(plant,  power_flow_data), zip(collect(values(ieee_14_Dyn_Nodes)), collect(values(pf_init_dict[ "bus_dict_init" ])))))
    
end


# ------------------------------------------------------
# ------------------------------------------------------
# Componets accessors function
# ------------------------------------------------------
# ------------------------------------------------------

function get_gen_Idx_and_vh( gen )
    
    # gen_Bus_num=@optic _.Bus_num
    
    # gen_vh=@optic _.vh

    # Bus_num = getall(gen, gen_Bus_num[1])
    
    # vh = getall(gen, gen_vh[1])
    
    
    return ( gen.Bus_num, gen.vh )
end


function get_gen_para_ra( gen )
    
    gen_para_ra=@optic _.ra
    
    return getall(gen, gen_para_ra[1])
end

function get_gen_para_X_d( gen )
    
    gen_para_X_d=@optic _.X_d
    
    return getall(gen, gen_para_X_d[1])
end


function get_gen_para_X_q( gen )
    
    gen_para_X_q=@optic _.X_q
    
    return getall(gen, gen_para_X_q[1])
end


function get_gen_para_X_d_dash( gen )
    
    gen_para_X_d_dash=@optic _.X_d_dash
    
    return getall(gen, gen_para_X_d_dash[1])
end


function get_gen_para_X_q_dash( gen )
    
    gen_para_X_q_dash=@optic _.X_q_dash
    
    return getall(gen, gen_para_X_q_dash[1])
end


# ------------------------------------------------------
# ------------------------------------------------------


function get_a_node_name( node )
    
    # gen_Bus_num=@optic _.Bus_num
    
    # gen_vh=@optic _.vh

    # Bus_num = getall(gen, gen_Bus_num[1])
    
    # vh = getall(gen, gen_vh[1])
    
    
    return node.name
end


function get_load_trans_node_Idx( node )
    
    # gen_Bus_num=@optic _.Bus_num
    
    # gen_vh=@optic _.vh

    # Bus_num = getall(gen, gen_Bus_num[1])
    
    # vh = getall(gen, gen_vh[1])
    
    
    return node.Bus_num
end


function get_load_trans_node_Idx_and_PQ( node )
        
    return (node.Bus_num, node.P, node.Q)
end


function get_load_trans_node_Idx_and_vlimits( node )
    
    return (node.Bus_num, node.vmin, node.vmax)
end

# ------------------------------------------------------
# ------------------------------------------------------


function get_component_params_sym( component )
    return component.param
end


# function get_component_params_value( component )
    
#     p_data = Union{ComplexF64,Float64}[]
    
#     params_sym = get_component_params_sym( component )
    
#     for a_param in params_sym
#         para = getproperty(component, a_param)
#         push!(p_data, para)
#     end
    
#     return p_data
# end


function get_component_params_value( component )
       
    return component.param_values
end

# ------------------------------------------------------
# ------------------------------------------------------


function get_component_name( component )
    
    lens_comp_name=@optic _.name
    
    return getall(component, lens_comp_name)[1]
end


function get_component_state_vars_syms( component )
    
    lens_comp_state_vars_syms=@optic _.state_vars_syms
    
    return getall(component, lens_comp_state_vars_syms)[1]
end


function get_component_stab_state_vars_syms( component )
    
    lens_comp_stab_state_vars_syms=@optic _.stab_state_vars_syms
    
    return getall(component, lens_comp_stab_state_vars_syms)[1]
end


function get_component_im_vars_syms( component )
    
    lens_comp_im_vars_syms=@optic _.im_vars_syms
    
    return getall(component, lens_comp_im_vars_syms)[1]
end


function get_component_im_algebraic_vars_syms( component )
    
    lens_comp_im_algebraic_vars_syms=@optic _.im_algebraic_vars_syms
    
    return getall(component, lens_comp_im_algebraic_vars_syms)[1]
end


function get_component_algebraic_vars_syms( component )
    
    lens_comp_algebraic_vars_syms=@optic _.algebraic_vars_syms
    
    return getall(component, lens_comp_algebraic_vars_syms)[1]
end


function get_component_syms( component )
    
    lens_comp_syms=@optic _.syms
    
    return getall(component, lens_comp_syms)[1]
end

function get_component_dict_state_syms( component )
    
    lens_comp_dict_state_syms=@optic _.dict_state_syms
    
    return getall(component, lens_comp_dict_state_syms)[1]
end


function get_component_P_Q_idx( component )
    
    lens_comp_P_Q_idx=@optic _.P_Q_idx
    
    return getall(component, lens_comp_P_Q_idx)[1]
end


# function get_component_P_Q_view( component; load = false )

#     # component = wk_nd.nodes["bus5"].Load
    
#     if load == true
#         idx = get_component_P_Q_idx(component)        
#         P_Q_load = -1 * component.param_values[idx]
#         P_Q_view = view(component.param_values, idx )
#         P_Q_view .= P_Q_load
#         return P_Q_view
#     else
#         return view(component.param_values, get_component_P_Q_idx( component ))
#     end
    
# end

function get_component_P_Q_view( component )

    return view(component.param_values, get_component_P_Q_idx( component ))
    
end


function get_component_vh( component )
    
    lens_comp_vh=@optic _.vh
    
    return getall(component, lens_comp_vh)[1]
end


function get_component_ra_Xd_Xq_idx( component )
    
    lens_comp_ra_Xd_Xq_idx=@optic _.ra_Xd_Xq_idx
    
    return getall(component, lens_comp_ra_Xd_Xq_idx)[1]
end


function get_component_ra_Xd_Xq_view( component )
        
    return view(component.param_values, get_component_ra_Xd_Xq_idx( component ))
end


function get_component_ra_Xd_Xq( component )
        
    return component.param_values[
        get_component_ra_Xd_Xq_idx( component )]
end


function get_component_ra_Xd_Xq_Xd_dash_Xq_dash_idx( component )
    
    lens_comp_ra_Xd_Xq_Xd_dash_Xq_dash_idx=@optic _.ra_Xd_Xq_Xd_dash_Xq_dash_idx
    
    return getall(component,
                  lens_comp_ra_Xd_Xq_Xd_dash_Xq_dash_idx)[1]
end


function get_component_ra_Xd_Xq_Xd_dash_Xq_dash_view( component )
        
    return view(component.param_values,
                get_component_ra_Xd_Xq_Xd_dash_Xq_dash_idx(
                    component ))
end


function get_component_ra_Xd_Xq_Xd_dash_Xq_dash( component )
        
    return component.param_values[
        get_component_ra_Xd_Xq_Xd_dash_Xq_dash_idx(
            component )]
end


function get_component_with_loc_load_bool( component )
    
    lens_comp_with_loc_load=@optic _.with_loc_load
    
    return getall(component, lens_comp_with_loc_load)[1]
end



function get_component_ur_ui_idx( component )
    
    lens_comp_ur_ui_idx=@optic _.ur_ui_idx
    
    return getall(component, lens_comp_ur_ui_idx)[1]
end

function get_component_dim( component )
    
    lens_comp_dim=@optic _.dim
    
    return getall(component, lens_comp_dim)[1]
end


function get_component_state_algebraic_vars_dim( component )
    
    lens_comp_dim=@optic _.dim
    
    return getall(component, lens_comp_dim)[1]
end


function get_component_state_vars_dim( component )
        
    return length(component.state_vars_syms)
end

function get_component_stab_state_vars_dim( component )
        
    return length(component.stab_state_vars_syms)
end


function get_component_im_vars_syms_dim( component )
        
    return length(component.im_vars_syms)
end


function get_component_im_algebraic_vars_syms_dim( component )
        
    return length(component.im_algebraic_vars_syms)
end




function get_component_dims(component)
    lens_comp_dims=@optic _.dims
    return getall(component, lens_comp_dims)[1]
end



function get_component_mass_matrix( component )
    
    lens_comp_mass_matrix=@optic _.mass_matrix
    
    return getall(component, lens_comp_mass_matrix)[1]
end



function get_a_gen_im_mass_matrix( component )

    
    state_vars_dim =
        get_component_state_vars_dim(
            component )
    
    im_algebraic_vars_syms_dim =
        get_component_im_algebraic_vars_syms_dim(
            component )     
    
    return DAE_MassMatrix(
        state_vars_dim,
        im_algebraic_vars_syms_dim )
end


function get_component_dae_var( component )
    
    lens_comp_dae_var=@optic _.dae_var
    
    return getall(component, lens_comp_dae_var)[1]
end


function get_component_dyn_func_by_idx(component;idx = 2)
    
    lens_comp_dyn_func=@optic _.func
    
    return getall(component,
                  lens_comp_dyn_func)[1][idx]
end


function get_component_dyn_func( component )
    
    lens_comp_dyn_func=@optic _.func
    
    return getall(component, lens_comp_dyn_func)[1][1]
end


function get_component_dyn_func_non_flat(component)
    
    lens_comp_dyn_func=@optic _.func
    
    return getall(component,
                  lens_comp_dyn_func)[1][2]
end


function get_component_control_sig_syms(component)
    
    lens_comp_control_sig_syms=@optic _.control_sig_syms
    
    return getall(component,
                  lens_comp_control_sig_syms)[1]
end

function get_component_control_sig(component)
    
    lens_comp_control_sig=@optic _.control_sig
    
    return getall(component,
                  lens_comp_control_sig)[1]
end



function get_component_output_sig_syms(component )
    
    lens_comp_output_sig_syms=@optic _.output_sig_syms
    
    return getall(component,
                  lens_comp_output_sig_syms)[1]
end


function get_component_output_sig(component )
    
    lens_comp_output_sig=@optic _.output_sig
    
    return getall(component,
                  lens_comp_output_sig)[1]
end



function get_component_dict_state_syms(component )
   
    return component.dict_state_syms
end



function get_component_dict_stab_state_syms(component )
   
    return component.dict_stab_state_syms
end


function get_component_cb_state_event_func(component)
   
    return component.cb_state_event_func
end


function get_component_cb_state_affect_func(component )
   
    return component.cb_state_affect_func
end


function get_component_cb_state_syms(component )
   
    return component.cb_state_syms
end



function get_component_cb_state_conditions(component )
   
    return component.cb_state_conditions
end



function get_component_cb_state_values(component )
   
    return component.cb_state_values
end


function get_component_cb_state_sym2Idx(component )
   
    return component.cb_state_sym2Idx
end



function get_component_cb_state_dim(component )
   
    return component.cb_state_dim
end


function get_component_cb_dyn_state_event_func(component )
   
    return component.cb_dyn_state_event_func
end


function get_component_cb_dyn_state_affect_func(
    component )
   
    return component.cb_dyn_state_affect_func
end

function get_component_cb_dyn_state_syms(component )
   
    return component.cb_dyn_state_syms
end


function get_component_cb_dict_dyn_state_syms2sw_Idx(
    component )
   
    return component.cb_dict_dyn_state_syms2sw_Idx
end


function get_component_cb_dyn_state_sw_Idx(
    component )
   
    return component.cb_dyn_state_sw_Idx
end

function get_component_cb_dyn_state_conditions(
    component )
   
    return component.cb_dyn_state_conditions
end


function get_component_cb_dyn_state_values(component )
   
    return component.cb_dyn_state_values
end


function get_component_cb_dyn_state_sym2Idx(component )
   
    return component.cb_dyn_state_sym2Idx
end


function get_component_cb_dyn_param_state_sw(component )
   
    return component.cb_dyn_param_state_sw
end


function get_component_cb_dyn_state_dim(component )
   
    return component.cb_dyn_state_dim
end


function get_component_system_matrices( component )
   
    return component.system_matrices
end


function get_component_system_stab_matrices( component )
   
    return component.system_stab_matrices
end


function get_component_func_system_matrices( component )
   
    return component.func_system_matrices( component )
end


function get_component_func_system_stab_matrices(
    component )
   
    return component.func_system_matrices( component )
end


function get_plant_pure_state_vars_idxs_in_state(
    plant )

    pure_state_var_idxs = []

    state_vars_syms =
        plant.state_vars_syms

    dict_state_syms =
        plant.dict_state_syms

    for state_sym in state_vars_syms
        
        push!(pure_state_var_idxs,
              dict_state_syms[state_sym])
    end
   
    return pure_state_var_idxs
end



function get_plant_stab_state_vars_idxs_in_state( plant )

    stab_state_var_idxs = []

    stab_state_vars_syms =
        plant.stab_state_vars_syms

    dict_state_syms =
        plant.dict_state_syms

    for state_sym in stab_state_vars_syms
        push!( stab_state_var_idxs,
               dict_state_syms[state_sym])
    end
   
    return stab_state_var_idxs
end



function get_gens_nodes_pure_state_vars_idxs_in_state(
    gens_nodes_collection )
   
    return map(
        get_plant_pure_state_vars_idxs_in_state,
        gens_nodes_collection)
end


function get_gens_nodes_stab_state_vars_idxs_in_state(
    gens_nodes_collection )
   
    return map(
        get_plant_stab_state_vars_idxs_in_state,
        gens_nodes_collection)
end




function get_gens_nodes_pure_state_init_from_pf(
    nodes_collection, bus_dict_init )

    if isa(nodes_collection,  OrderedDict)

        comp_collection_values =
            collect(values(nodes_collection))

        gens_nodes_collection = [
            comp for comp  in comp_collection_values
                if comp.Bus_type == :Generator ]

        gens_nodes_pure_state_vars_idxs_in_state =
            get_gens_nodes_pure_state_vars_idxs_in_state(
                gens_nodes_collection )

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection,
                bus_dict_init )
   
        return [
            a_node_init.x0[state_idx]
            for (a_node_init, state_idx) in
                zip( gens_nodes_init_data_from_pf,
                     gens_nodes_pure_state_vars_idxs_in_state)]
        
    elseif  isa(nodes_collection, Union{Array, Vector})

        gens_nodes_collection = [
            comp for comp  in nodes_collection
                if comp.Bus_type == :Generator]

        gens_nodes_pure_state_vars_idxs_in_state =
            get_gens_nodes_pure_state_vars_idxs_in_state(
                gens_nodes_collection )

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection, bus_dict_init )
   
        return [
            a_node_init.x0[state_idx]
            for (a_node_init, state_idx) in
                zip( gens_nodes_init_data_from_pf,
                     gens_nodes_pure_state_vars_idxs_in_state)]

    else
        
        comp_collection_values =
            collect(values(nodes_collection))

        gens_nodes_collection = [
            comp for comp  in comp_collection_values
                if comp.Bus_type == :Generator]

        gens_nodes_pure_state_vars_idxs_in_state =
            get_gens_nodes_pure_state_vars_idxs_in_state(
                gens_nodes_collection )

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection, bus_dict_init )
   
        return [
            a_node_init.x0[state_idx]
            for (a_node_init, state_idx) in
                zip( gens_nodes_init_data_from_pf,
                     gens_nodes_pure_state_vars_idxs_in_state)]
    end
        
end


function get_gens_nodes_stab_state_init_from_pf(
    nodes_collection, bus_dict_init )

    if isa(nodes_collection,  OrderedDict)

        comp_collection_values =
            collect(values(nodes_collection))

        gens_nodes_collection = [
            comp for comp  in comp_collection_values
                if comp.Bus_type == :Generator ]

        gens_nodes_stab_state_vars_idxs_in_state =
            get_gens_nodes_stab_state_vars_idxs_in_state(
                gens_nodes_collection )

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection, bus_dict_init )
   
        return [
            a_node_init.x0[state_idx]
            for (a_node_init, state_idx) in
                zip( gens_nodes_init_data_from_pf,
                     gens_nodes_stab_state_vars_idxs_in_state )]
        
    elseif  isa(nodes_collection, Union{Array, Vector})

        gens_nodes_collection = [
            comp for comp  in nodes_collection
                if comp.Bus_type == :Generator]

        gens_nodes_stab_state_vars_idxs_in_state =
            get_gens_nodes_stab_state_vars_idxs_in_state(
                gens_nodes_collection )

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection, bus_dict_init )
   
        return [
            a_node_init.x0[state_idx]
            for (a_node_init, state_idx) in
                zip( gens_nodes_init_data_from_pf,
                     gens_nodes_stab_state_vars_idxs_in_state)]

    else
        
        comp_collection_values =
            collect(values(nodes_collection))

        gens_nodes_collection = [
            comp for comp  in comp_collection_values
                if comp.Bus_type == :Generator]

        gens_nodes_stab_state_vars_idxs_in_state =
            get_gens_nodes_stab_state_vars_idxs_in_state(
                gens_nodes_collection )

        gens_nodes_init_data_from_pf =
            get_gens_nodes_init_data_from_pf(
                gens_nodes_collection, bus_dict_init )
   
        return [
            a_node_init.x0[state_idx]
            for (a_node_init, state_idx) in
                zip( gens_nodes_init_data_from_pf,
                     gens_nodes_stab_state_vars_idxs_in_state)]
    end
        
end


##############################################

function get_components_params_sym(comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(get_component_params_sym,
                   collect(comp_collection) )
    else
        return map(get_component_params_sym,
                   collect(values(comp_collection)))
    end
end


function get_components_state_algebraic_vars_dim(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_state_algebraic_vars_dim,
            collect(comp_collection) )
    else
        
        return map(
            get_component_state_algebraic_vars_dim,
            collect(values(comp_collection)))
    end
end


function get_components_state_vars_dim(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        return map(
            get_component_state_vars_dim,
            collect(comp_collection) )
    else
        
        return map(
            get_component_state_vars_dim,
            collect(values(comp_collection)))
    end
end


function get_components_stab_state_vars_dim(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_stab_state_vars_dim,
            collect(comp_collection) )
    else
        return map(
            get_component_stab_state_vars_dim,
            collect(values(comp_collection)))
    end
end


function get_components_im_vars_syms_dim(
    comp_collection)
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_im_vars_syms_dim,
            collect(comp_collection) )
    else
        return map(
            get_component_im_vars_syms_dim,
            collect(values(comp_collection)))
    end
end


function get_components_im_algebraic_vars_syms_dim(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_im_algebraic_vars_syms_dim,
            collect(comp_collection) )
    else
        return map(
            get_component_im_algebraic_vars_syms_dim,
            collect(values(comp_collection)))
    end
end


function get_components_state_vars(
    comp_collection)
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_state_vars_syms,
            collect(comp_collection) )
    else
        return map(
            get_component_state_vars_syms,
            collect(values(comp_collection)))
    end
end


function get_components_stab_state_vars(comp_collection)
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_stab_state_vars_syms,
            collect(comp_collection) )
    else
        return map(
            get_component_stab_state_vars_syms,
            collect(values(comp_collection)))
    end
end


function get_components_im_vars(comp_collection)
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_im_vars_syms,
            collect(comp_collection) )
    else
        return map(
            get_component_im_vars_syms,
            collect(values(comp_collection)))
    end
end


function get_components_im_vars_syms(
    comp_collection)
    if isa(comp_collection, Union{Array, Vector})
        
        return map(get_component_im_vars_syms,
                   collect(comp_collection) )
    else
        return map(
            get_component_im_vars_syms,
            collect(values(comp_collection)))
    end
end


function get_components_im_algebraic_vars_syms(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_im_algebraic_vars_syms,
            collect(comp_collection) )
    else
        
        return map(
            get_component_im_algebraic_vars_syms,
            collect(values(comp_collection)))
    end
end

function get_components_dict_state_syms(comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        return map(
            get_component_dict_state_syms,
            collect(comp_collection) )
    else
        return map(
            get_component_dict_state_syms,
            collect(values(comp_collection)))
    end
end


function get_components_dict_stab_state_syms(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_dict_stab_state_syms,
            collect(comp_collection) )
    else
        return map(
            get_component_dict_stab_state_syms,
            collect(values(comp_collection)))
    end
end

#-------------------------------------------------------


function get_network_bus_names(
    comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            comp.Bus_type == :Generator ?
                get_a_node_name( comp.Gen ) :
                comp.Bus_type == :Load ?
                get_a_node_name( comp.Load ) :
                get_a_node_name( comp.Trans )
            for comp in comp_collection_values]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            comp.Bus_type == :Generator ?
                get_a_node_name( comp.Gen ) :
                comp.Bus_type == :Load ?
                get_a_node_name( comp.Load ) :
                get_a_node_name( comp.Trans )
            for comp in comp_collection ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [
            comp.Bus_type == :Generator ?
                get_a_node_name( comp.Gen ) :
                comp.Bus_type == :Load ?
                get_a_node_name( comp.Load ) :
                get_a_node_name( comp.Trans )
            for comp in comp_collection_values ]
        
    end
    
end


function get_generators_bus_names(comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            get_a_node_name( comp.Gen )
            for comp in comp_collection_values
                if comp.Bus_type == :Generator]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            get_a_node_name(comp.Gen )
            for comp in comp_collection
                if comp.Bus_type == :Generator]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [
            get_a_node_name(comp.Gen )
            for comp in comp_collection_values
                if comp.Bus_type == :Generator]
        
    end
end


function get_non_generators_bus_names(
    comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            comp.Bus_type == :Load ?
                get_a_node_name( comp.Load ) :
                get_a_node_name( comp.Trans )
            for comp in comp_collection_values
                if comp.Bus_type != :Generator ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [ comp.Bus_type == :Load ?
            get_a_node_name( comp.Load ) :
            get_a_node_name( comp.Trans )
                  for comp in comp_collection
                      if comp.Bus_type != :Generator ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [
            comp.Bus_type == :Load ?
                get_a_node_name( comp.Load ) :
                get_a_node_name( comp.Trans )
            for comp in comp_collection_values
                if comp.Bus_type != :Generator ]
        
    end
end



function get_Loads_bus_names(comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            get_a_node_name( comp.Load )
            for comp in comp_collection_values
                if comp.Bus_type == :Load]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            get_a_node_name(comp.Load )
            for comp in comp_collection
                if comp.Bus_type == :Load]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [
            get_a_node_name(comp.Load )
            for comp in comp_collection_values
                if comp.Bus_type == :Load]
        
    end
end


function get_Trans_bus_names(comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            get_a_node_name( comp.Trans )
            for comp in comp_collection_values
                if comp.Bus_type == :Transmission]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            get_a_node_name(comp.Trans )
            for comp in comp_collection
                if comp.Bus_type == :Transmission]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [
            get_a_node_name(comp.Trans )
            for comp in comp_collection_values
                if comp.Bus_type == :Transmission]
        
    end
end

#-------------------------------------------------------

function get_components_idx_and_Bus_type(
    comp_collection)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        [( comp.Bus_num, comp.Bus_type)
         for comp in
             comp_collection_values]
        
        return [
            ( comp.Bus_num, comp.Bus_type)
            for comp in
                comp_collection_values]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            ( comp.Bus_num, comp.Bus_type)
            for comp in comp_collection]
        
    else
        
        return  [
            ( comp.Bus_num, comp.Bus_type)
            for comp in
                comp_collection_values]
        
    end
end


function get_components_P_Q_idx(
    comp_collection)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return map(
            (comp) -> comp.Bus_type == :Generator ?
                get_component_P_Q_idx(comp.Gen) :
                comp.Bus_type == :Load ?
                get_component_P_Q_idx(comp.Load) :
                get_component_P_Q_idx(comp.Trans),
            comp_collection_values )
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return map(
            (comp) -> comp.Bus_type == :Generator ?
                get_component_P_Q_idx(comp.Gen) :
                comp.Bus_type == :Load ?
                get_component_P_Q_idx(comp.Load) :
                get_component_P_Q_idx(comp.Trans),
            collect(comp_collection) )
    else
        
        return map(
            (comp) -> comp.Bus_type == :Generator ?
                get_component_P_Q_idx(comp.Gen) :
                comp.Bus_type == :Load ?
                get_component_P_Q_idx(comp.Load) :
                get_component_P_Q_idx(comp.Trans),
            collect(values(comp_collection)))
    end
end


function get_components_P_Q_view(
    comp_collection;
    nodes_P_Q_idx = nothing)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return map(
            (comp) -> comp.Bus_type == :Generator ?
                get_component_P_Q_view(comp.Gen) :
                comp.Bus_type == :Load ?
                get_component_P_Q_view(comp.Load) :
                get_component_P_Q_view(comp.Trans),
            comp_collection_values)
        
    elseif isa(comp_collection,  Vector{Vector} )

        # Gen param are in nodes.param_values[1]
        return [
            typeof(a_param_values) == Vector{Float64} ?
                view(a_param_values, P_Q_idx )  :
                view(a_param_values[1], P_Q_idx)
            for (a_param_values, P_Q_idx) in
                zip(comp_collection, nodes_P_Q_idx) ]

    elseif isa(comp_collection, Union{Array, Vector})
        
        return map(
            (comp) -> comp.Bus_type == :Generator ?
                get_component_P_Q_view(comp.Gen) :
                comp.Bus_type == :Load ?
                get_component_P_Q_view(comp.Load) :
                get_component_P_Q_view(comp.Trans),
            comp_collection )

    else
        
        return map(
            (comp) -> comp.Bus_type == :Generator ?
                get_component_P_Q_view(comp.Gen) :
                comp.Bus_type == :Load ?
                get_component_P_Q_view(comp.Load) :
                get_component_P_Q_view(comp.Trans),
            collect(values(comp_collection)))
    end
end


function get_generators_Idx_and_vh(
    comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            get_gen_Idx_and_vh( comp.Gen )
            for comp in comp_collection_values
                if comp.Bus_type == :Generator]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            get_gen_Idx_and_vh(comp.Gen )
            for comp in comp_collection
                if comp.Bus_type == :Generator]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [
            get_gen_Idx_and_vh(comp.Gen )
            for comp in comp_collection_values
                if comp.Bus_type == :Generator]
        
    end
end



function get_non_slack_generators_Idx_and_vh(
    comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            get_gen_Idx_and_vh( comp.Gen )
            for comp in comp_collection_values
                if comp.Bus_type == :Generator &&
                    comp.isa_slack == false ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            get_gen_Idx_and_vh(comp.Gen )
            for comp in comp_collection
                if comp.Bus_type == :Generator &&
                    comp.isa_slack == false ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [
            get_gen_Idx_and_vh(comp.Gen )
            for comp in comp_collection_values
                if comp.Bus_type == :Generator &&
                    comp.isa_slack == false ]
        
    end
end


#############
# with_loc_load


function get_Idx_and_loc_load(
    comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            comp.Bus_type == :Generator &&
                comp.with_loc_load == true  ?
                ( comp.Bus_num, comp.Loc_load.loc_P +
                im *  comp.Loc_load.loc_Q ) :
                (comp.Bus_num, 0.0 + im * 0.0 )
            for comp in comp_collection_values  ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            comp.Bus_type == :Generator &&
                comp.with_loc_load == true  ?
                ( comp.Bus_num, comp.Loc_load.loc_P +
                im *  comp.Loc_load.loc_Q ) :
                (comp.Bus_num, 0.0 + im * 0.0 )
            for comp in comp_collection  ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [
            comp.Bus_type == :Generator &&
                comp.with_loc_load == true  ?
                ( comp.Bus_num, comp.Loc_load.loc_P +
                im *  comp.Loc_load.loc_Q ) :
                (comp.Bus_num, 0.0 + im * 0.0 )
            for comp in comp_collection_values  ]
        
    end
end



function get_gens_with_loc_load_Idx_and_loc_load(
    comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return Tuple{Int64, Vector{Float64}}[
            ( comp.Bus_num,
              [comp.Loc_load.loc_P, comp.Loc_load.loc_Q] )
            for comp in comp_collection_values
                if comp.Bus_type == :Generator &&
                    comp.with_loc_load == true  ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  Tuple{Int64, Vector{Float64}}[
            ( comp.Bus_num,
              [comp.Loc_load.loc_P, comp.Loc_load.loc_Q] )
            for comp in comp_collection
                if comp.Bus_type == :Generator &&
                    comp.with_loc_load == true  ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return   Tuple{Int64, Vector{Float64}}[
            ( comp.Bus_num,
              [comp.Loc_load.loc_P, comp.Loc_load.loc_Q] )
            for comp in comp_collection_values
                if comp.Bus_type == :Generator &&
                    comp.with_loc_load == true  ]
        
    end
end




function get_load_trans_nodes_Idx(
    comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            comp.Bus_num
            for comp in comp_collection_values
                if comp.Bus_type != :Generator  ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            comp.Bus_num
            for comp in comp_collection
                if comp.Bus_type != :Generator ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        return  [
            comp.Bus_num
            for comp in comp_collection_values
                if comp.Bus_type != :Generator  ]        
    end
end



function get_load_trans_nodes_Idx_and_PQ(comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [ comp.Bus_type == :Transmission ?
            get_load_trans_node_Idx_and_PQ( comp.Trans ) :
            get_load_trans_node_Idx_and_PQ( comp.Load )
                 for comp in comp_collection_values
                     if comp.Bus_type != :Generator  ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [
            comp.Bus_type == :Transmission ?
                get_load_trans_node_Idx_and_PQ( comp.Trans ) :
                get_load_trans_node_Idx_and_PQ( comp.Load )
            for comp in comp_collection
                if comp.Bus_type != :Generator  ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [
            comp.Bus_type == :Transmission ?
                get_load_trans_node_Idx_and_PQ( comp.Trans ) :
                get_load_trans_node_Idx_and_PQ( comp.Load )
            for comp in comp_collection_values
                if comp.Bus_type != :Generator  ]
        
    end
end



function get_load_trans_nodes_Idx_and_vlimits(
    comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            comp.Bus_type == :Transmission ?
                get_load_trans_node_Idx_and_vlimits(
                    comp.Trans ) :
                        get_load_trans_node_Idx_and_vlimits(
                            comp.Load )
            for comp in comp_collection_values
                if comp.Bus_type != :Generator  ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [ comp.Bus_type == :Transmission ?
            get_load_trans_node_Idx_and_vlimits( comp.Trans ) :
            get_load_trans_node_Idx_and_vlimits( comp.Load )
                 for comp in comp_collection
                     if comp.Bus_type != :Generator  ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return [ comp.Bus_type == :Transmission ?
            get_load_trans_node_Idx_and_vlimits( comp.Trans ) :
            get_load_trans_node_Idx_and_vlimits( comp.Load )
                 for comp in comp_collection_values
                     if comp.Bus_type != :Generator  ]
        
    end
end


#-----------------------------------------------------
#-----------------------------------------------------

function get_component_param_idxs_in_param_values(
    comp_collection, comp,
    para::Union{String,Symbol})
    
    list_comps_idx   =
        [a_comp.Bus_num
         for a_comp in
             collect(values(comp_collection)) ]
    
    list_comps_names =
        [a_comp.name
         for a_comp in
             collect(values(comp_collection))]
    
    vec_comps_params_syms =
        get_components_params_sym(comp_collection)
    
    comp_idx = findfirst(
        comp_name -> comp_name == string(comp.name),
        list_comps_names )
    
    list_comp_params =
        vec_comps_params_syms[comp_idx]

    if para ∉ [list_comp_params...;]
        
        # I used zero because I want to use it for test
        # since indices do not start from zero
        
        return (comp_idx, Int64[0]) 
        
    else
        comp_para_idx =
            nd_findfirst(
                comp_para -> comp_para == Symbol(para),
                collect(list_comp_params))
        
        return (comp_idx, comp_para_idx)
    end
   
end


function get_components_param_idxs_in_param_values(
    comp_collection,
    para::Union{String,Symbol} )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
            get_component_param_idxs_in_param_values(
                comp_collection, comp, para )
            for comp in  comp_collection_values  ]
       
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
            get_component_param_idxs_in_param_values(
                comp_collection, comp, para )
            for comp in  comp_collection  ]
        
    else
        comp_collection_values = collect(values(comp_collection))
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
            get_component_param_idxs_in_param_values(
                comp_collection, comp, para )
            for comp in
                comp_collection_values  ]

        
    end
   
end



function get_streamlined_components_param_idxs_in_param_values(
    comp_collection,
    para::Union{String,Symbol},
    gens_only = false, non_gens_only = false  )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        if gens_only == true
            
            return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
                get_component_param_idxs_in_param_values(
                    comp_collection, comp, para )
                for comp in  comp_collection_values
                    if comp.Bus_type == :Generator ]
            
        elseif non_gens_only == true
            
            return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
                get_component_param_idxs_in_param_values(
                    comp_collection, comp, para )
                for comp in  comp_collection_values
                    if comp.Bus_type != :Generator  ]
            
        else
            return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
                get_component_param_idxs_in_param_values(
                    comp_collection, comp, para )
                for comp in
                    comp_collection_values  ]
            
        end

    elseif isa(comp_collection, Union{Array, Vector})


        if gens_only == true

            return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
                get_component_param_idxs_in_param_values(
                    comp_collection, comp, para )
                for comp in  comp_collection
                    if comp.Bus_type == :Generator  ]
            
        elseif non_gens_only == true

            return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
                get_component_param_idxs_in_param_values(
                    comp_collection, comp, para )
                for comp in  comp_collection
                    if comp.Bus_type != :Generator ]
            
        else
            
            return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
                get_component_param_idxs_in_param_values(
                    comp_collection, comp, para )
                for comp in  comp_collection  ]
            
        end

    else
        
        comp_collection_values =
            collect(values(comp_collection))

        if gens_only == true
            
            return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
                get_component_param_idxs_in_param_values(
                    comp_collection, comp, para )
                for comp in  comp_collection_values
                    if comp.Bus_type == :Generator  ]
            
        elseif non_gens_only == true

            return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
                get_component_param_idxs_in_param_values(
                    comp_collection, comp, para )
                for comp in  comp_collection_values
                    if comp.Bus_type != :Generator   ]
            
        else

            return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[
                get_component_param_idxs_in_param_values(
                    comp_collection, comp, para )
                for comp in
                    comp_collection_values  ]
            
        end

    end
   
end


#---------------------------------------------------------


function get_streamlined_components_param_view_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values =
        get_components_param_idxs_in_param_values(
            comp_collection, para )

    return [ typeof( idx ) == Tuple{Int64, Int64} ?
        view( a_param_values, last(idx) ) :
        view( a_param_values[last(idx)[1]] , last(idx)[2] )
             for (a_param_values, idx ) in
                 zip( param_values_collection,
                      param_idxs_in_param_values )
                 if !( typeof(idx) ==
                     Tuple{Int64, Vector{Int64}} ) ]

end



function get_streamlined_components_param_view_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values =
        get_components_param_idxs_in_param_values(
            comp_collection, para )

    return [ typeof( idx ) == Tuple{Int64, Int64} ?
        view( a_param_values, last(idx) ) :
        view( a_param_values[last(idx)[1]] , last(idx)[2] )
             for (a_param_values, idx ) in
                 zip( param_values_collection,
                      param_idxs_in_param_values )
                 if !(typeof(idx) ==
                     Tuple{Int64, Vector{Int64}} ) ]

end

#---------------------------------------------------------

function get_components_param_view_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values =
        get_components_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof(idx) == Tuple{Int64, Vector{Int64}} ?
            [0.0] : typeof( idx ) == Tuple{Int64, Int64} ?
            view( a_param_values, last(idx) ) :
            view( a_param_values[last(idx)[1]] , last(idx)[2] )
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values ) ]

end


function get_components_param_view_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values =
        get_components_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof(idx) == Tuple{Int64, Vector{Int64}} ?
            [0.0] : typeof( idx ) == Tuple{Int64, Int64} ?
            view( a_param_values, last(idx) ) :
            view( a_param_values[last(idx)[1]] , last(idx)[2] )
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values ) ]

end


function get_components_idx_and_param_view_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    components_param_view_in_param_values =
        get_components_param_view_in_param_values(
            param_values_collection, comp_collection, para )

    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_view)
            for (comp, a_view) in
                zip( comp_collection_values,
                     components_param_view_in_param_values ) ]
        
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_view)
            for (comp, a_view) in
                zip( comp_collection,
                     components_param_view_in_param_values ) ]
        
    else
        
        comp_collection_values = collect(values(comp_collection))

        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_view)
            for (comp, a_view) in
                zip( comp_collection_values,
                     components_param_view_in_param_values ) ]
       
    end
   

end


function get_components_idx_and_param_view_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)

    components_param_view_in_param_values =
        get_components_param_view_in_param_values(
            param_values_collection,
            comp_collection, para )

    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_view)
            for (comp, a_view) in
                zip( comp_collection_values,
                     components_param_view_in_param_values ) ]
        
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_view)
            for (comp, a_view) in
                zip( comp_collection,
                     components_param_view_in_param_values ) ]
        
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_view)
            for (comp, a_view) in
                zip( comp_collection_values,
                     components_param_view_in_param_values ) ]
       
    end


end

#---------------------------------------------------------
#---------------------------------------------------------


function get_components_param_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values =
        get_components_param_idxs_in_param_values(
            comp_collection, para )

    return [ typeof(idx) == Tuple{Int64, Vector{Int64}} ?
        [0.0] :
        typeof( idx ) == Tuple{Int64, Int64} ?
        a_param_values[ last(idx) ] :
        a_param_values[last(idx)[1]][ last(idx)[2] ]
             for (a_param_values, idx ) in
                 zip( param_values_collection,
                      param_idxs_in_param_values ) ]

end


function get_components_param_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values =
        get_components_param_idxs_in_param_values(
            comp_collection, para )

    return [ typeof(idx) == Tuple{Int64, Vector{Int64}} ?
        [0.0] :
        typeof( idx ) == Tuple{Int64, Int64} ?
        a_param_values[ last(idx) ] :
        a_param_values[last(idx)[1]][ last(idx)[2] ]
             for (a_param_values, idx ) in
                 zip( param_values_collection,
                      param_idxs_in_param_values ) ]

end


function get_components_idx_and_param_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    components_param_in_param_values =
        get_components_param_in_param_values(
            param_values_collection,
            comp_collection, para )

    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_para)
            for (comp, a_para) in
                zip( comp_collection_values,
                     components_param_in_param_values ) ]
        
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_para)
            for (comp, a_para) in
                zip( comp_collection,
                     components_param_in_param_values ) ]
        
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_para)
            for (comp, a_para) in
                zip( comp_collection_values,
                     components_param_in_param_values ) ]
       
    end
   

end


function get_components_idx_and_param_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)


    components_param_in_param_values =
        get_components_param_in_param_values(
            param_values_collection,
            comp_collection, para )

    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_para)
            for (comp, a_para) in
                zip( comp_collection_values,
                     components_param_in_param_values ) ]
        
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_para)
            for (comp, a_para) in
                zip( comp_collection,
                     components_param_in_param_values ) ]
        
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return Union{ Tuple{Int64, Tuple{Int64, Vector{Int64}}}, Tuple{Int64, Tuple{Int64, Int64}}, Tuple{Int64, Tuple{Int64, Tuple{Int64, Int64}}} }[
            (comp.Bus_num, a_para)
            for (comp, a_para) in
                zip( comp_collection_values,
                     components_param_in_param_values ) ]
       
    end
    

end


#---------------------------------------------------------


function get_nodes_idx_and_Yshunt( comp_collection )

    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))
        
        return Tuple{Int64, ComplexF64 }[ comp.Bus_type == :Generator ? (comp.Bus_num, im * comp.Gen.Y_n) : comp.Bus_type == :Load ? (comp.Bus_num, im * comp.Load.Y_n) : comp.Bus_type == :Transmission ? (comp.Bus_num, im * comp.Trans.Y_n) : (0, im * 0.0)  for comp in  comp_collection_values ]
        
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Tuple{Int64, ComplexF64 }[ comp.Bus_type == :Generator ? (comp.Bus_num, im * comp.Gen.Y_n) : comp.Bus_type == :Load ? (comp.Bus_num, im * comp.Load.Y_n) : comp.Bus_type == :Transmission ? (comp.Bus_num, im * comp.Trans.Y_n) : (0, im * 0.0)  for comp in  comp_collection ]
        
    else
        
        comp_collection_values = collect(values(comp_collection))

        return Tuple{Int64, ComplexF64 }[ comp.Bus_type == :Generator ? (comp.Bus_num, im * comp.Gen.Y_n) : comp.Bus_type == :Load ? (comp.Bus_num, im * comp.Load.Y_n) : comp.Bus_type == :Transmission ? (comp.Bus_num, im * comp.Trans.Y_n) : (0, im * 0.0)  for comp in  comp_collection_values ]
       
    end

end


#---------------------------------------------------------
#---------------------------------------------------------

function get_components_P_param_view_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_view_in_param_values(
        param_values_collection,
        comp_collection; para = :P)

end


function get_components_P_param_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_in_param_values(
        param_values_collection,
        comp_collection; para = :P)

end

function get_components_Q_param_view_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_view_in_param_values(
        param_values_collection,
        comp_collection; para = :Q)

end


function get_components_Q_param_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_in_param_values(
        param_values_collection,
        comp_collection; para = :Q)

end


function get_components_ra_param_view_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_view_in_param_values(
        param_values_collection,
        comp_collection; para = :ra)

end



function get_components_ra_param_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_in_param_values(
        param_values_collection,
        comp_collection; para = :ra)

end


function get_components_Xd_param_view_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_view_in_param_values(
        param_values_collection,
        comp_collection; para = :X_d)

end


function get_components_Xd_param_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_in_param_values(
        param_values_collection,
        comp_collection; para = :X_d)

end


function get_components_Xq_param_view_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_view_in_param_values(
        param_values_collection,
        comp_collection; para = :X_q)

end


function get_components_Xq_param_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_in_param_values(
        param_values_collection,
        comp_collection; para = :X_q)

end



function get_components_X_d_dash_param_view_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_view_in_param_values(
        param_values_collection,
        comp_collection; para = :X_d_dash)

end


function get_components_X_d_dash_param_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_in_param_values(
        param_values_collection,
        comp_collection; para = :X_d_dash)

end



function get_components_X_q_dash_param_view_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_view_in_param_values(
        param_values_collection,
        comp_collection; para = :X_q_dash)

end


function get_components_X_q_dash_param_in_param_values(
    param_values_collection,
    comp_collection)

    return get_components_param_in_param_values(
        param_values_collection,
        comp_collection; para = :X_q_dash)

end



#-----------------------------------------------------
# Generators
#-----------------------------------------------------

function get_gens_param_view_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values =
        get_gens_param_idxs_in_param_values(
            comp_collection, para )

    return [ typeof(idx) == Tuple{Int64, Vector{Int64}} ? [0.0] : typeof( idx ) == Tuple{Int64, Int64} ? view( a_param_values, last(idx) ) : view( a_param_values[last(idx)[1]] , last(idx)[2] ) for (a_param_values, idx ) in zip( param_values_collection, param_idxs_in_param_values ) ]

end



function get_gens_param_view_in_param_values(
    param_values_collection,
    comp_collection;
    para::Union{String,Symbol} =
        nothing)

    param_idxs_in_param_values =
        get_gens_param_idxs_in_param_values(
            comp_collection, para )

    return [ typeof(idx) == Tuple{Int64, Vector{Int64}} ?
        [0.0] :
        typeof( idx ) == Tuple{Int64, Int64} ?
        view( a_param_values, last(idx) ) :
        view( a_param_values[last(idx)[1]] , last(idx)[2] )
             for (a_param_values, idx ) in
                 zip( param_values_collection,
                      param_idxs_in_param_values ) ]
    

end


function get_streamlined_gens_param_view_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values = get_gens_param_idxs_in_param_values(comp_collection, para )

    return [ typeof( idx ) == Tuple{Int64, Int64} ? view( a_param_values, last(idx) ) : view( a_param_values[last(idx)[1]] , last(idx)[2] ) for (a_param_values, idx ) in zip( param_values_collection, param_idxs_in_param_values ) if !( typeof(idx) == Tuple{Int64, Vector{Int64}}) ]

end


function get_streamlined_gens_param_view_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values = get_gens_param_idxs_in_param_values(comp_collection, para )

    return [  typeof( idx ) == Tuple{Int64, Int64} ? view( a_param_values, last(idx) ) : view( a_param_values[last(idx)[1]] , last(idx)[2] ) for (a_param_values, idx ) in zip( param_values_collection, param_idxs_in_param_values ) if !(typeof(idx) == Tuple{Int64, Vector{Int64}}) ]
    

end

#-----------------------------------------------------


function get_gens_param_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values =
        get_gens_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof(idx) == Tuple{Int64, Vector{Int64}} ?
            [0.0] :
            typeof( idx ) == Tuple{Int64, Int64} ?
            a_param_values[ last(idx) ] :
             a_param_values[last(idx)[1]][ last(idx)[2] ]
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values ) ]

end



function get_gens_param_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values =
        get_gens_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof(idx) == Tuple{Int64, Vector{Int64}} ?
            [0.0] : typeof( idx ) == Tuple{Int64, Int64} ?
            a_param_values[ last(idx) ] :
            a_param_values[last(idx)[1]][ last(idx)[2] ]
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values ) ]
    
end


function get_streamlined_gens_param_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values =
        get_gens_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof( idx ) == Tuple{Int64, Int64} ?
            a_param_values[ last(idx) ] :
            a_param_values[last(idx)[1]][ last(idx)[2] ]
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values )
            if !( typeof(idx) == Tuple{Int64, Vector{Int64}}) ]

end


function get_streamlined_gens_param_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values =
        get_gens_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof( idx ) == Tuple{Int64, Int64} ?
            a_param_values[ last(idx) ] :
            a_param_values[last(idx)[1]][ last(idx)[2] ]
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values )
            if !(typeof(idx) ==
                Tuple{Int64, Vector{Int64}}) ]
    

end


#-----------------------------------------------------


function get_gens_param_idxs_in_param_values(
    comp_collection,
    para::Union{String,Symbol} )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ comp.Bus_type == :Generator ? get_component_param_idxs_in_param_values(comp_collection, comp, para ) : (comp.Bus_num, Int64[0])  for comp in  comp_collection_values  ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ comp.Bus_type == :Generator ? get_component_param_idxs_in_param_values(comp_collection, comp, para ) : (comp.Bus_num, Int64[0])  for comp in  comp_collection  ]
        
    else
        
        comp_collection_values = collect(values(comp_collection))
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ comp.Bus_type == :Generator ? get_component_param_idxs_in_param_values(comp_collection, comp, para ) : (comp.Bus_num, Int64[0])  for comp in  comp_collection_values  ]
        
    end
   
end



function get_streamlined_gens_param_idxs_in_param_values(
    comp_collection,
    para::Union{String,Symbol} )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[  get_component_param_idxs_in_param_values(comp_collection, comp, para )  for comp in  comp_collection_values if comp.Bus_type == :Generator ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[  get_component_param_idxs_in_param_values(comp_collection, comp, para )  for comp in  comp_collection if comp.Bus_type == :Generator   ]
        
    else
        comp_collection_values = collect(values(comp_collection))
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ get_component_param_idxs_in_param_values(comp_collection, comp, para )  for comp in  comp_collection_values if  comp.Bus_type == :Generator   ]
        
    end
   
end


function get_streamlined_non_gens_param_idxs_in_param_values(
    comp_collection,
    para::Union{String,Symbol} )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[  get_component_param_idxs_in_param_values(comp_collection, comp, para )  for comp in  comp_collection_values if comp.Bus_type != :Generator ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[  get_component_param_idxs_in_param_values(comp_collection, comp, para )  for comp in  comp_collection if comp.Bus_type != :Generator   ]
        
    else
        comp_collection_values = collect(values(comp_collection))
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ get_component_param_idxs_in_param_values(comp_collection, comp, para )  for comp in  comp_collection_values if  comp.Bus_type != :Generator   ]
        
    end
   
end


#-----------------------------------------------------
# Generators local loads
#-----------------------------------------------------


function get_gens_loc_load_param_idxs_in_param_values(
    comp_collection, para::Union{String,Symbol} )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ comp.Bus_type == :Generator && comp.with_loc_load == true  ? get_component_param_idxs_in_param_values(comp_collection, comp, para ) : (comp.Bus_num, Int64[0])  for comp in  comp_collection_values  ]
        # map((comp) -> get_component_param_idxs_in_param_values(comp_collection, comp, para ),  comp_collection_values )
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ comp.Bus_type == :Generator && comp.with_loc_load == true  ? get_component_param_idxs_in_param_values(comp_collection, comp, para ) : (comp.Bus_num, Int64[0])  for comp in  comp_collection  ]
        # map((comp) -> get_component_param_idxs_in_param_values(comp_collection, comp, para ),  comp_collection )
    else
        comp_collection_values = collect(values(comp_collection))
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ comp.Bus_type == :Generator && comp.with_loc_load == true  ? get_component_param_idxs_in_param_values(comp_collection, comp, para ) : (comp.Bus_num, Int64[0])  for comp in  comp_collection_values  ]
        # map((comp) -> get_component_param_idxs_in_param_values(comp_collection, comp, para ),  comp_collection_values )
    end
   
end



function get_gens_loc_load_param_view_in_param_values(
    param_values_collection,
    comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values = get_gens_loc_load_param_idxs_in_param_values(comp_collection, para )

    return [ typeof(idx) == Tuple{Int64, Vector{Int64}} ? [0.0] : typeof( idx ) == Tuple{Int64, Int64} ? view( a_param_values, last(idx) ) : view( a_param_values[last(idx)[1]] , last(idx)[2] ) for (a_param_values, idx ) in zip( param_values_collection, param_idxs_in_param_values ) ]

end



function get_gens_loc_load_param_view_in_param_values(
    param_values_collection,
    comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values = get_gens_loc_load_param_idxs_in_param_values(comp_collection, para )

    return [ typeof(idx) == Tuple{Int64, Vector{Int64}} ? [0.0] : typeof( idx ) == Tuple{Int64, Int64} ? view( a_param_values, last(idx) ) : view( a_param_values[last(idx)[1]] , last(idx)[2] ) for (a_param_values, idx ) in zip( param_values_collection, param_idxs_in_param_values ) ]
    

end

#-----------------------------------------------------



function get_gens_loc_load_param_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values =
        get_gens_loc_load_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof(idx) == Tuple{Int64, Vector{Int64}} ?
            [0.0] :
            typeof( idx ) == Tuple{Int64, Int64} ?
            a_param_values[ last(idx) ] :
            a_param_values[last(idx)[1]][ last(idx)[2] ]
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values ) ]

end



function get_gens_loc_load_param_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values =
        get_gens_loc_load_param_idxs_in_param_values(
            comp_collection, para )

    return [ typeof(idx) == Tuple{Int64, Vector{Int64}} ?
        [0.0] :
        typeof( idx ) == Tuple{Int64, Int64} ?
        a_param_values[ last(idx) ] :
        a_param_values[last(idx)[1]][ last(idx)[2] ]
             for (a_param_values, idx ) in
                 zip( param_values_collection,
                      param_idxs_in_param_values ) ]

end


#-----------------------------------------------------
# Non-Generators
#-----------------------------------------------------


function get_non_gens_param_idxs_in_param_values(
    comp_collection, para::Union{String,Symbol} )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ comp.Bus_type != :Generator ? get_component_param_idxs_in_param_values(comp_collection, comp, para ) : (comp.Bus_num, Int64[0])  for comp in  comp_collection_values  ]
        # map((comp) -> get_component_param_idxs_in_param_values(comp_collection, comp, para ),  comp_collection_values )
    elseif isa(comp_collection, Union{Array, Vector})
        
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ comp.Bus_type != :Generator ? get_component_param_idxs_in_param_values(comp_collection, comp, para ) : (comp.Bus_num, Int64[0])  for comp in  comp_collection  ]
        # map((comp) -> get_component_param_idxs_in_param_values(comp_collection, comp, para ),  comp_collection )
    else
        comp_collection_values = collect(values(comp_collection))
        return Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}[ comp.Bus_type != :Generator ? get_component_param_idxs_in_param_values(comp_collection, comp, para ) : (comp.Bus_num, Int64[0])  for comp in  comp_collection_values  ]
        # map((comp) -> get_component_param_idxs_in_param_values(comp_collection, comp, para ),  comp_collection_values )
    end
   
end

#-----------------------------------------------------


function get_non_gens_param_view_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values =
        get_non_gens_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof(idx) == Tuple{Int64, Vector{Int64}} ?
            [0.0] : typeof( idx ) == Tuple{Int64, Int64} ?
            view( a_param_values, last(idx) ) :
            view( a_param_values[last(idx)[1]] , last(idx)[2] )
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values ) ]

end



function get_non_gens_param_view_in_param_values(
    param_values_collection, comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values =
        get_non_gens_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof(idx) == Tuple{Int64, Vector{Int64}} ?
            [0.0] : typeof( idx ) == Tuple{Int64, Int64} ?
            view( a_param_values, last(idx) ) :
            view( a_param_values[last(idx)[1]] , last(idx)[2] )
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values ) ]
    
end

#-----------------------------------------------------


function get_non_gens_param_in_param_values(
    param_values_collection, comp_collection,
    para::Union{String,Symbol})

    param_idxs_in_param_values =
        get_non_gens_param_idxs_in_param_values(
            comp_collection,
            para )

    return [
        typeof(idx) == Tuple{Int64, Vector{Int64}} ?
            [0.0] :
            typeof( idx ) == Tuple{Int64, Int64} ?
            a_param_values[ last(idx) ] :
             a_param_values[last(idx)[1]][ last(idx)[2] ]
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values ) ]

end


function get_non_gens_param_in_param_values(
    param_values_collection,
    comp_collection;
    para::Union{String,Symbol} = nothing)

    param_idxs_in_param_values =
        get_non_gens_param_idxs_in_param_values(
            comp_collection, para )

    return [
        typeof(idx) == Tuple{Int64, Vector{Int64}} ?
            [0.0] :
            typeof( idx ) == Tuple{Int64, Int64} ?
            a_param_values[ last(idx) ] :
            a_param_values[last(idx)[1]][ last(idx)[2] ]
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 param_idxs_in_param_values ) ]
    
end


#---------------------------------------------------------
#---------------------------------------------------------

function get_components_params_idxs_in_param_values(
    comp_collection;
    param_list::Vector{Symbol} = [nothing])

    # #-----------------------------
    # comp_collection = wk_nd.nodes
    # param_list = [:P, :Q]
    # #-----------------------------
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))

        param_idxs_list = Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_components_param_idxs_in_param_values( comp_collection_values, a_param ))
        
        end

        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] ))   for params_Idx_tup in zip( Tuple(param_idxs_list)... )]

    elseif isa(comp_collection, Union{Array, Vector})

        param_idxs_list =  Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_components_param_idxs_in_param_values(comp_collection, a_param ))
        
        end
        
        
        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] )) for params_Idx_tup in zip( Tuple(param_idxs_list)... )]
        
    else
        comp_collection_values = collect(values(comp_collection))

        param_idxs_list =  Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_components_param_idxs_in_param_values( comp_collection_values, a_param ))
        
        end

        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] )) for params_Idx_tup in zip( Tuple(param_idxs_list)... )]

    end
   
end

#---------------------------------------------------------


function get_components_params_view_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values =
        get_components_params_idxs_in_param_values(
            comp_collection; param_list = param_list)
    
    return [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ? [0.0] : typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) ]

end

#---------------------------------------------------------

function get_components_params_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values =
        get_components_params_idxs_in_param_values(
            comp_collection; param_list = param_list)
    
    return [
        typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ?
            [0.0] :
            typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ?
            a_param_values[ last(idx)[1:end] ] :
            a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 params_idxs_in_param_values ) ]

end


#---------------------------------------------------------
# Generators
#---------------------------------------------------------


function get_gens_params_idxs_in_param_values(
    comp_collection;
    param_list::Vector{Symbol} = [nothing])

    # #-----------------------------
    # comp_collection = wk_nd.nodes
    # param_list = [:P, :Q]
    # #-----------------------------
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(
            values(comp_collection))

        param_idxs_list = Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_gens_param_idxs_in_param_values( comp_collection_values, a_param ))
        
        end

        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] ))   for params_Idx_tup in zip( Tuple(param_idxs_list)... )]

    elseif isa(comp_collection, Union{Array, Vector})

        param_idxs_list =  Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_gens_param_idxs_in_param_values(comp_collection, a_param ))
        
        end
        
        
        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] )) for params_Idx_tup in zip( Tuple(param_idxs_list)... )]
        
    else
        comp_collection_values = collect(values(comp_collection))

        param_idxs_list =  Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_gens_param_idxs_in_param_values( comp_collection_values, a_param ))
        
        end

        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] )) for params_Idx_tup in zip( Tuple(param_idxs_list)... )]

    end
   
end

#---------------------------------------------------------
#---------------------------------------------------------

function get_gens_params_view_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values =
        get_gens_params_idxs_in_param_values(
            comp_collection;
            param_list = param_list)
    
    return [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ? [0.0] : typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) ]

end



function get_gens_params_view_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing],
    gens_view_only = false )

    params_idxs_in_param_values = get_gens_params_idxs_in_param_values(comp_collection; param_list = param_list)
    if gens_view_only == false

        return [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ? [0.0] : typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) ]
        
    else
        
        return [  typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) if  last(idx)[1] != 0  ]
    end
            
end


function get_industrial_model_gens_params_view_in_param_values( param_values_collection, comp_collection; param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values = get_gens_params_idxs_in_param_values(comp_collection; param_list = param_list)
    
    return filter(!isnothing, [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ? nothing : typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) ])

end



function get_im_model_gens_params_view_in_param_values( param_values_collection, comp_collection; param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values = get_gens_params_idxs_in_param_values(comp_collection; param_list = param_list)
    
    return filter(!isnothing, [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ? nothing : typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) ])

end


#---------------------------------------------------------


function get_gens_params_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values = get_gens_params_idxs_in_param_values(comp_collection; param_list = param_list)
    
    return [
        typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ?
            [0.0] :
            typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ?
            a_param_values[ last(idx)[1:end] ] :
            a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
        for (a_param_values, idx ) in
            zip(
                param_values_collection,
                params_idxs_in_param_values ) ]

end


function get_industrial_model_gens_params_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values =
        get_gens_params_idxs_in_param_values(
            comp_collection;
            param_list = param_list)
    
    return filter(
        !isnothing,
        [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ?
            nothing :
            typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ?
            a_param_values[last(idx)[1:end] ] :
            a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
          for (a_param_values, idx ) in
              zip(
                  param_values_collection,
                  params_idxs_in_param_values ) ])

end



function get_im_model_gens_params_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values =
        get_gens_params_idxs_in_param_values(
            comp_collection;
            param_list = param_list)
    
    return filter(
        !isnothing,
        [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ?
            nothing :
            typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ?
            a_param_values[ last(idx)[1:end] ] :
            a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
          for (a_param_values, idx ) in
              zip(
                  param_values_collection,
                  params_idxs_in_param_values ) ])

end



function get_gens_params_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing],
    gens_view_only = false )

    params_idxs_in_param_values =
        get_gens_params_idxs_in_param_values(
            comp_collection;
            param_list = param_list)
    
    if gens_view_only == false

        return [
            typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ?
                [0.0] :
                typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ?
                a_param_values[ last(idx)[1:end] ] :
                a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
            for (a_param_values, idx ) in
                zip(
                    param_values_collection,
                    params_idxs_in_param_values ) ]
        
    else
        
        return [
            typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ?
                a_param_values[ last(idx)[1:end] ] :
                a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
            for (a_param_values, idx ) in
                zip(
                    param_values_collection,
                    params_idxs_in_param_values )
                if  last(idx)[1] != 0  ]
    end
            
end


#---------------------------------------------------------
# Generators local load
#---------------------------------------------------------


function get_gens_loc_load_params_idxs_in_param_values(
    comp_collection;
    param_list::Vector{Symbol} = [nothing])

    # #-----------------------------
    # comp_collection = wk_nd.nodes
    # param_list = [:P, :Q]
    # #-----------------------------
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))

        param_idxs_list = Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_gens_loc_load_param_idxs_in_param_values( comp_collection_values, a_param ))
        
        end

        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] ))   for params_Idx_tup in zip( Tuple(param_idxs_list)... )]

    elseif isa(comp_collection, Union{Array, Vector})

        param_idxs_list =  Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_gens_loc_load_param_idxs_in_param_values(comp_collection, a_param ))
        
        end
        
        
        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] )) for params_Idx_tup in zip( Tuple(param_idxs_list)... )]
        
    else
        comp_collection_values = collect(values(comp_collection))

        param_idxs_list =  Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_gens_loc_load_param_idxs_in_param_values( comp_collection_values, a_param ))
        
        end

        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] )) for params_Idx_tup in zip( Tuple(param_idxs_list)... )]

    end
   
end


function get_gens_loc_load_params_view_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values = get_gens_loc_load_params_idxs_in_param_values(comp_collection; param_list = param_list)
    
    return [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ? [0.0] : typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) ]

end



function get_streamlined_gens_loc_load_params_view_in_param_values( param_values_collection, comp_collection; param_list::Vector{Symbol} = [nothing],  gens_view_only = true  )

    params_idxs_in_param_values = get_gens_loc_load_params_idxs_in_param_values(comp_collection; param_list = param_list)
    
    if  gens_view_only == false
        return [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ? [0.0] : typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) ]
    else
            
        return [  typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) if !(typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0)  ]
    end
    

end

#---------------------------------------------------------


function get_gens_loc_load_params_in_param_values(
    param_values_collection, comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values =
        get_gens_loc_load_params_idxs_in_param_values(
            comp_collection; param_list = param_list)
    
    return [
        typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ?
            [0.0] :
            typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ?
            a_param_values[ last(idx)[1:end] ] :
            a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
        
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 params_idxs_in_param_values ) ]

end



function get_streamlined_gens_loc_load_params_in_param_values(
    param_values_collection, comp_collection;
    param_list::Vector{Symbol} = [nothing],
    gens_view_only = true  )

    params_idxs_in_param_values =
        get_gens_loc_load_params_idxs_in_param_values(
            comp_collection; param_list = param_list)
    
    if  gens_view_only == false
        
        return [
            typeof( last(idx) ) == Vector{Int64} &&
                last(idx)[1] == 0  ?
                [0.0] :
                typeof( last(idx) ) == Vector{Int64} &&
                last(idx)[1] != 0  ?
                
                a_param_values[ last(idx)[1:end] ] :
                a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
            for (a_param_values, idx ) in
                zip( param_values_collection,
                     params_idxs_in_param_values ) ]
    else
            
        return [
            typeof( last(idx) ) == Vector{Int64} &&
                last(idx)[1] != 0  ?
                a_param_values[ last(idx)[1:end] ] :
                a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
            for (a_param_values, idx ) in
                zip( param_values_collection,
                     params_idxs_in_param_values )
                if !(typeof( last(idx) ) == Vector{Int64} &&
                    last(idx)[1] == 0)  ]
    end
    

end


#---------------------------------------------------------
# non-generators
#---------------------------------------------------------


function get_non_gens_params_idxs_in_param_values(
    comp_collection;
    param_list::Vector{Symbol} = [nothing])

    # #-----------------------------
    # comp_collection = wk_nd.nodes
    # param_list = [:P, :Q]
    # #-----------------------------
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))

        param_idxs_list = Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_non_gens_param_idxs_in_param_values( comp_collection_values, a_param ))
        
        end

        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] ))   for params_Idx_tup in zip( Tuple(param_idxs_list)... )]

    elseif isa(comp_collection, Union{Array, Vector})

        param_idxs_list =  Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_non_gens_param_idxs_in_param_values(comp_collection, a_param ))
        
        end
        
        
        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] )) for params_Idx_tup in zip( Tuple(param_idxs_list)... )]
        
    else
        comp_collection_values = collect(values(comp_collection))

        param_idxs_list =  Vector{ Union{Tuple{Int64, Vector{Int64}},Tuple{Int64, Int64}, Tuple{Int64, Tuple{Int64, Int64}}}}[]

        for a_param in param_list
        
            push!(param_idxs_list, get_non_gens_param_idxs_in_param_values( comp_collection_values, a_param ))
        
        end

        return [ typeof(first(params_Idx_tup)) == Tuple{Int64, Vector{Int64}} ? (first(first( params_Idx_tup )), Int64[0] ) : typeof(first(params_Idx_tup)) == Tuple{Int64, Int64} ?  ( first(first( params_Idx_tup )),  [last(a_para_idx) for a_para_idx in params_Idx_tup ] ) : (first(first( params_Idx_tup )), ( first(last(first( params_Idx_tup ) )), [last(last(a_para_idx)) for a_para_idx in params_Idx_tup ] )) for params_Idx_tup in zip( Tuple(param_idxs_list)... )]

    end
   
end


function get_non_gens_params_view_in_param_values(
    param_values_collection,
    comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values =
        get_non_gens_params_idxs_in_param_values(
            comp_collection;
            param_list = param_list)
    
    return [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ? [0.0] : typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) ]

end


function get_streamlined_non_gens_params_view_in_param_values( param_values_collection, comp_collection; param_list::Vector{Symbol} = [nothing], non_gens_view_only = true )

    params_idxs_in_param_values = get_non_gens_params_idxs_in_param_values(comp_collection; param_list = param_list)

    if non_gens_view_only == false

        return [ typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ? [0.0] : typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) ]
        
    else
        return [  typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ? view( a_param_values, last(idx)[1:end] ) : view( a_param_values[last(idx)[1]] , last(idx)[2][1:end] ) for (a_param_values, idx ) in zip( param_values_collection, params_idxs_in_param_values ) if !(typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0 ) ]
    end
    

end


#---------------------------------------------------------


function get_non_gens_params_in_param_values(
    param_values_collection, comp_collection;
    param_list::Vector{Symbol} = [nothing] )

    params_idxs_in_param_values =
        get_non_gens_params_idxs_in_param_values(
            comp_collection; param_list = param_list)
    
    return [
        typeof( last(idx) ) == Vector{Int64} && last(idx)[1] == 0  ?
            [0.0] :
            typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ?
            a_param_values[ last(idx)[1:end] ] :
            a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
        for (a_param_values, idx ) in
            zip( param_values_collection,
                 params_idxs_in_param_values ) ]

end


function get_streamlined_non_gens_params_in_param_values(
    param_values_collection, comp_collection;
    param_list::Vector{Symbol} = [nothing],
    non_gens_view_only = true )

    params_idxs_in_param_values =
        get_non_gens_params_idxs_in_param_values(
            comp_collection; param_list = param_list)

    if non_gens_view_only == false

        return [
            typeof( last(idx) ) == Vector{Int64} &&
                last(idx)[1] == 0  ?
            [0.0] :
            typeof( last(idx) ) == Vector{Int64} && last(idx)[1] != 0  ?
            a_param_values[ last(idx)[1:end] ] :
            a_param_values[last(idx)[1]][ last(idx)[2][1:end] ]
                 for (a_param_values, idx ) in
                     zip( param_values_collection,
                          params_idxs_in_param_values ) ]
        
    else
        
        return [
            typeof( last(idx) ) == Vector{Int64} &&
                last(idx)[1] != 0  ?
                a_param_values[ last(idx)[1:end] ] :
                a_param_values[last(idx)[1]][last(idx)[2][1:end]]
            for (a_param_values, idx ) in
                zip( param_values_collection,
                     params_idxs_in_param_values )
                if !(typeof( last(idx) ) ==
                    Vector{Int64} && last(idx)[1] == 0 ) ]
    end
    

end



#---------------------------------------------------------
#---------------------------------------------------------


function get_components_isa_generator_bool(
    comp_collection )
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [comp.Bus_type == :Generator ?
            true : false
                for comp in comp_collection_values ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return  [comp.Bus_type == :Generator ?
            true : false
                 for comp in comp_collection ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return  [comp.Bus_type == :Generator ?
            true : false
                 for comp in comp_collection_values ]
        
    end
end

#---------------------------------------------------------


function get_components_Idx_and_ra_Xd_Xq_idx(
    comp_collection)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            comp.Bus_type == :Generator ?
                (comp.Bus_num,
                 get_component_ra_Xd_Xq_idx(comp.Gen) ) :
                     (comp.Bus_num, [])
            for (idx, comp) in
                enumerate(comp_collection_values) ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [
            comp.Bus_type == :Generator ?
                (comp.Bus_num,
                 get_component_ra_Xd_Xq_idx(comp.Gen) ) :
                     (comp.Bus_num, [])
            for (idx, comp) in
                enumerate(comp_collection) ]
    else
        comp_collection_values =
            collect(values(comp_collection))
        
        return [comp.Bus_type == :Generator ?
            (comp.Bus_num,
             get_component_ra_Xd_Xq_idx(comp.Gen) ) :
                 (comp.Bus_num, [])
                for (idx, comp) in
                    enumerate(comp_collection_values) ]
        
    end
end

#---------------------------------------------------------

function get_components_Idx_and_ra_Xd_Xq_view(
    comp_collection; nodes_ra_Xd_Xq_idx = nothing,
    isa_generator_bool = nothing)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return [
            comp.Bus_type == :Generator ?
                (comp.Bus_num,
                 get_component_ra_Xd_Xq_view(comp.Gen)) :
                     (comp.Bus_num, [])
            for (idx, comp) in
                enumerate(comp_collection_values) ]

    elseif isa(comp_collection,  Vector{Vector} )
        
        return [
            (typeof(a_param_values) == Vector{Float64} && isa_gen == true) ?
                (first(ra_Xd_Xq_idx),
                 view(a_param_values, last(ra_Xd_Xq_idx) ) ) :
                     isa_gen == true ?
                     (first(ra_Xd_Xq_idx),
                      view(a_param_values[1], last(ra_Xd_Xq_idx) ) ) :
                          (first(ra_Xd_Xq_idx), [])
            for (isa_gen, ra_Xd_Xq_idx, a_param_values) in
                zip(isa_generator_bool,
                    nodes_ra_Xd_Xq_idx, comp_collection) ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [
            comp.Bus_type == :Generator ?
                (comp.Bus_num,
                 get_component_ra_Xd_Xq_view(comp.Gen)) :
                     (comp.Bus_num, [])
            for (idx, comp) in
                enumerate(comp_collection) ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        return [
            comp.Bus_type == :Generator ?
                ( comp.Bus_num,
                  get_component_ra_Xd_Xq_view(comp.Gen)) :
                      ( comp.Bus_num, [])
            for (idx, comp) in
                enumerate(comp_collection_values) ]
    end
end


#------------------------------------------------------


function get_components_Idx_and_ra_Xd_Xq(
    comp_collection; nodes_ra_Xd_Xq_idx = nothing,
    isa_generator_bool = nothing)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return [
            comp.Bus_type == :Generator ?
                (comp.Bus_num,
                 get_component_ra_Xd_Xq(comp.Gen)) :
                     (comp.Bus_num, [])
            for (idx, comp) in
                enumerate(comp_collection_values) ]

    elseif isa(comp_collection,  Vector{Vector} )
        
        return [
            (typeof(a_param_values) ==
                Vector{Float64} && isa_gen == true) ?
                (first(ra_Xd_Xq_idx),
                  a_param_values[ last(ra_Xd_Xq_idx) ] ) :
                     isa_gen == true ?
                     (first(ra_Xd_Xq_idx),
                      a_param_values[1][last(ra_Xd_Xq_idx)] ) :
                          (first(ra_Xd_Xq_idx), [])
            for (isa_gen, ra_Xd_Xq_idx, a_param_values) in
                zip(isa_generator_bool,
                    nodes_ra_Xd_Xq_idx,
                    comp_collection) ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [
            comp.Bus_type == :Generator ?
                (comp.Bus_num,
                 get_component_ra_Xd_Xq(comp.Gen)) :
                     (comp.Bus_num, [])
            for (idx, comp) in
                enumerate(comp_collection) ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return [
            comp.Bus_type == :Generator ?
                ( comp.Bus_num,
                  get_component_ra_Xd_Xq(comp.Gen)) :
                      ( comp.Bus_num, [])
            for (idx, comp) in
                enumerate(comp_collection_values) ]
    end
end



#------------------------------------------------------


function get_components_Idx_and_ra_Xd_Xq_Xd_dash_Xq_dash_idx(comp_collection)
    
    if isa(comp_collection,  OrderedDict)
        comp_collection_values = collect(values(comp_collection))

        return [comp.Bus_type == :Generator ? ( comp.Bus_num, get_component_ra_Xd_Xq_Xd_dash_Xq_dash_idx(comp.Gen) ) : ( comp.Bus_num, [])  for (idx, comp) in enumerate( comp_collection_values ) ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [comp.Bus_type == :Generator ? (comp.Bus_num, get_component_ra_Xd_Xq_Xd_dash_Xq_dash_idx(comp.Gen) ) : (comp.Bus_num, [])  for (idx, comp) in enumerate(comp_collection) ]
    else
        comp_collection_values = collect(values(comp_collection))
        return [comp.Bus_type == :Generator ? (comp.Bus_num, get_component_ra_Xd_Xq_Xd_dash_Xq_dash_idx(comp.Gen) ) : (comp.Bus_num, [])  for (idx, comp) in enumerate(comp_collection_values) ]
        
    end
end


function get_components_Idx_and_ra_Xd_Xq_Xd_dash_Xq_dash_idx_view(comp_collection; nodes_ra_Xd_Xq_Xd_dash_Xq_dash_idx = nothing, isa_generator_bool = nothing)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))
        return [comp.Bus_type == :Generator ? ( comp.Bus_num, get_component_ra_Xd_Xq_Xd_dash_Xq_dash_view(comp.Gen)) : ( comp.Bus_num, [])  for (idx, comp) in enumerate(comp_collection_values) ]

    elseif isa(comp_collection,  Vector{Vector} )
        
        return [ (typeof(a_param_values) == Vector{Float64} && isa_gen == true) ? (first( ra_Xd_Xq_Xd_dash_Xq_dash_idx ), view(a_param_values, last( ra_Xd_Xq_Xd_dash_Xq_dash_idx ) ) ) : isa_gen == true ? (first( ra_Xd_Xq_Xd_dash_Xq_dash_idx ), view(a_param_values[1], last( ra_Xd_Xq_Xd_dash_Xq_dash_idx ) ) ) : (first( ra_Xd_Xq_Xd_dash_Xq_dash_idx ), [])  for (isa_gen, ra_Xd_Xq_Xd_dash_Xq_dash_idx, a_param_values) in  zip(isa_generator_bool, nodes_ra_Xd_Xq_Xd_dash_Xq_dash_idx, comp_collection) ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        return [ comp.Bus_type == :Generator ? ( comp.Bus_num, get_component_ra_Xd_Xq_Xd_dash_Xq_dash_view(comp.Gen)) : ( comp.Bus_num, [])  for (idx, comp) in enumerate(comp_collection) ]
    else
        comp_collection_values = collect(values(comp_collection))
        return [comp.Bus_type == :Generator ? ( comp.Bus_num, get_component_ra_Xd_Xq_Xd_dash_Xq_dash_view(comp.Gen)) : ( comp.Bus_num, [])  for (idx, comp) in enumerate(comp_collection_values) ]
    end
end


#------------------------------------------------------


function get_components_Idx_and_ra_Xd_Xq_Xd_dash_Xq_dash_idx(
    comp_collection;
    nodes_ra_Xd_Xq_Xd_dash_Xq_dash_idx = nothing,
    isa_generator_bool = nothing)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return [comp.Bus_type == :Generator ?
            ( comp.Bus_num,
              get_component_ra_Xd_Xq_Xd_dash_Xq_dash(
                  comp.Gen)) :
                      ( comp.Bus_num, [])
                for (idx, comp) in
                    enumerate(comp_collection_values) ]

    elseif isa(comp_collection,  Vector{Vector} )        
        
        return [
            (typeof(a_param_values) ==
                Vector{Float64} && isa_gen == true) ?
                (first( ra_Xd_Xq_Xd_dash_Xq_dash_idx ),
                 a_param_values[
                     last( ra_Xd_Xq_Xd_dash_Xq_dash_idx ) ] ) :
                     isa_gen == true ?
                     (first( ra_Xd_Xq_Xd_dash_Xq_dash_idx ),
                      a_param_values[1][
                          last(ra_Xd_Xq_Xd_dash_Xq_dash_idx )]) :
                          (first( ra_Xd_Xq_Xd_dash_Xq_dash_idx ), [])
            for (isa_gen,
                 ra_Xd_Xq_Xd_dash_Xq_dash_idx,
                 a_param_values) in
                zip(isa_generator_bool,
                    nodes_ra_Xd_Xq_Xd_dash_Xq_dash_idx,
                    comp_collection) ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [ comp.Bus_type == :Generator ?
            ( comp.Bus_num,
              get_component_ra_Xd_Xq_Xd_dash_Xq_dash(
                  comp.Gen)) : ( comp.Bus_num, [])
                 for (idx, comp) in
                     enumerate(comp_collection) ]
    else
        comp_collection_values =
            collect(values(comp_collection))
        return [
            comp.Bus_type == :Generator ? (
                comp.Bus_num,
                get_component_ra_Xd_Xq_Xd_dash_Xq_dash(
                    comp.Gen)) : ( comp.Bus_num, [])
            for (idx, comp) in
                enumerate(comp_collection_values) ]
    end
end


#------------------------------------------------------


function get_components_with_loc_load_bool(
    comp_collection)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return map(
            get_component_with_loc_load_bool,
            comp_collection_values)
        
    elseif isa(comp_collection, Union{Array, Vector})

        return map(
            get_component_with_loc_load_bool,
            comp_collection)

    else
        return map(
            get_component_with_loc_load_bool,
            collect(values(comp_collection)))
    end
end


function get_components_loc_load_param_idx(
    comp_collection)
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        with_loc_load_bool =
            get_components_with_loc_load_bool(
                comp_collection_values )
                
        return [
            with_loc_load == true ?
                a_plant.loc_load_param_idx : 0
            for (with_loc_load, a_plant) in
                zip(with_loc_load_bool,
                    comp_collection_values) ]
        
    elseif isa(comp_collection, Union{Array, Vector})

        with_loc_load_bool =
            get_components_with_loc_load_bool(
                comp_collection)
        
        return [with_loc_load == true ?
            a_plant.loc_load_param_idx : 0
                for (with_loc_load, a_plant) in
                    zip(with_loc_load_bool,
                        comp_collection) ]

    else
        comp_collection_values =
            collect(values(comp_collection))
        
        with_loc_load_bool =
            get_components_with_loc_load_bool(
                comp_collection_values  )
        
        return [
            with_loc_load == true ?
                a_plant.loc_load_param_idx : 0
            for (with_loc_load, a_plant) in
                zip(with_loc_load_bool,
                    comp_collection_values) ]
    end
end

#------------------------------------------------------

function get_components_loc_Load_P_Q_view(
    comp_collection;
    nodes_with_loc_load_bool = nothing,
    nodes_loc_load_param_idx = nothing)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values = collect(values(comp_collection))
        return map((comp) -> (comp.Bus_type == :Generator  && comp.with_loc_load) ? get_component_P_Q_view(comp.Loc_load) : [0.0, 0.0], comp_collection_values)
        
    elseif isa(comp_collection,  Vector{Vector} )
        
        return [ with_loc_load == true ? view(a_param_values[loc_load_param_idx], 1:2) :  [0.0, 0.0]  for (with_loc_load, loc_load_param_idx, a_param_values) in  zip( nodes_with_loc_load_bool, nodes_loc_load_param_idx, comp_collection) ]        

    elseif isa(comp_collection, Union{Array, Vector})

        return map((comp) -> (comp.Bus_type == :Generator  && comp.with_loc_load) ? get_component_P_Q_view(comp.Loc_load) : [0.0, 0.0], comp_collection)
        
    else
        
        comp_collection_values = collect(values(comp_collection))
        
        return map((comp) -> (comp.Bus_type == :Generator  && comp.with_loc_load) ? get_component_P_Q_view(comp.Loc_load) : [0.0, 0.0], comp_collection_values)
    end
end

#------------------------------------------------------


function get_components_loc_Load_P_Q(
    comp_collection;
    nodes_with_loc_load_bool = nothing,
    nodes_loc_load_param_idx = nothing)
    
    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return map(
            (comp) -> (comp.Bus_type == :Generator  &&
                comp.with_loc_load) ?
                get_component_P_Q(comp.Loc_load) :
                [0.0, 0.0],
            comp_collection_values)
        
    elseif isa(comp_collection,  Vector{Vector} )
        
        return [
            with_loc_load == true ?
                a_param_values[loc_load_param_idx][ 1:2] :
                [0.0, 0.0]
            for (with_loc_load, loc_load_param_idx,
                 a_param_values) in                
                zip( nodes_with_loc_load_bool,
                     nodes_loc_load_param_idx,
                     comp_collection) ]        

    elseif isa(comp_collection, Union{Array, Vector})

        return map(
            (comp) -> (comp.Bus_type == :Generator  &&
                comp.with_loc_load) ?
                get_component_P_Q(comp.Loc_load) :
                [0.0, 0.0],
            comp_collection)
        
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return map(
            (comp) -> (comp.Bus_type == :Generator  &&
                comp.with_loc_load) ?
                get_component_P_Q(comp.Loc_load) :
                [0.0, 0.0],
            comp_collection_values)
    end
end

#------------------------------------------------------

function get_components_ur_ui_idx(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_ur_ui_idx,
            collect(comp_collection) )
    else
        
        return map(
            get_component_ur_ui_idx,
            collect(values(comp_collection)))
    end
end


function get_components_ur_ui_Idx_in_state(
    comp_collection; counter=0 )
    
    if isa(comp_collection, Union{Array, Vector})

        comps =
            collect(comp_collection)
        
        comps_ur_ui_Idx =
            get_components_ur_ui_idx(comps)
        
        dims =
            get_components_state_algebraic_vars_dim(comps)
        
        return get_ur_ui_Idx_in_state(
            dims, comps_ur_ui_Idx;
            counter=counter )
    else
        
        comps =
            collect(values(comp_collection))
        
        comps_ur_ui_Idx =
            get_components_ur_ui_idx(comps)
        
        dims =
            get_components_state_algebraic_vars_dim(
                comps)
        
        return get_ur_ui_Idx_in_state(
            dims, comps_ur_ui_Idx;
            counter=counter )        

    end        
end


function get_components_slack_Idx(
    comp_collection)

    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return map(
            (comp) -> comp.Bus_num,
            [ comp  for comp in
                 comp_collection_values
                 if (comp.Bus_type == :Generator &&
                     comp.isa_slack == true) ])[1]
        
    elseif isa(comp_collection, Union{Array, Vector})
        return map(
            (comp) -> comp.Bus_num,
            [ comp
              for comp in comp_collection
                  if (comp.Bus_type == :Generator &&
                      comp.isa_slack == true) ])[1]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return map(
            (comp) -> comp.Bus_num,
            [ comp
              for comp in comp_collection_values
                  if (comp.Bus_type == :Generator &&
                      comp.isa_slack == true) ])[1]
    end
end



function get_components_slack_vh(
    comp_collection)

    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return map(
            (comp) -> get_component_vh( comp.Gen ),
            [ comp
              for comp in
                  comp_collection_values
                  if (comp.Bus_type == :Generator &&
                      comp.isa_slack == true) ])[1]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return map(
            (comp) -> get_component_vh( comp.Gen ),
            [ comp
              for comp in comp_collection
                  if (comp.Bus_type == :Generator &&
                      comp.isa_slack == true) ])[1]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return map(
            (comp) -> get_component_vh( comp.Gen ),
            [ comp
              for comp in comp_collection_values
                  if (comp.Bus_type == :Generator &&
                      comp.isa_slack == true) ])[1]
    end
end


function get_components_no_slack_ur_ui_idx(
    comp_collection)

    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            get_component_ur_ui_idx(comp)
            for comp in comp_collection_values
                if (comp.Bus_type == :Generator &&
                    comp.isa_slack == false) ||
                    (comp.Bus_type != :Generator)  ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [
            get_component_ur_ui_idx(comp)
            for comp in comp_collection
                if (comp.Bus_type == :Generator &&
                    comp.isa_slack == false) ||
                    (comp.Bus_type != :Generator)  ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return [
            get_component_ur_ui_idx(comp)
            for comp in comp_collection_values
                if (comp.Bus_type == :Generator &&
                    comp.isa_slack == false) ||
                    (comp.Bus_type != :Generator)  ]
    end
end



function get_components_no_slack_ur_ui_Idx_in_state(
    comp_collection; counter=0 )

    if isa(comp_collection,  OrderedDict)

        comp_collection_values =
            collect(values(comp_collection))

        slack_bus_idx =
            get_components_slack_Idx(
                comp_collection_values )

        if slack_bus_idx == 1

            slack_comp_dims =
                sum([ comp.dim
                      for comp in
                          comp_collection_values
                          if (comp.Bus_type == :Generator &&
                              comp.isa_slack == true) ])

            non_slack_comp_collection =
                [
                    comp for comp in
                        comp_collection_values
                        if (comp.Bus_type == :Generator &&
                            comp.isa_slack == false)||
                            (comp.Bus_type != :Generator)  ]

            no_slack_comps_ur_ui_Idx =
                get_components_no_slack_ur_ui_idx(
                    non_slack_comp_collection)

            no_slack_comps_dims =
                get_components_state_algebraic_vars_dim(
                    non_slack_comp_collection)

            new_counter = counter + slack_comp_dims

            return get_non_slack_ur_ui_Idx_in_state(
                no_slack_comps_dims,
                no_slack_comps_ur_ui_Idx;
                counter = new_counter )
        else

            slack_comp_dims =
                sum(
                    [ comp.dim
                      for comp in
                          comp_collection_values
                          if (comp.Bus_type == :Generator &&
                              comp.isa_slack == true) ])

            non_slack_comp_collection =
                [comp
                 for comp in comp_collection_values
                     if (comp.Bus_type == :Generator &&
                         comp.isa_slack == false)||
                         (comp.Bus_type != :Generator)  ]

            comp_collection_before_slack =
                [comp
                 for  comp in
                     non_slack_comp_collection
                     if  comp.Bus_num < slack_bus_idx ]
            
            comp_collection_after_slack =
                [comp
                 for  comp in
                     non_slack_comp_collection
                     if comp.Bus_num > slack_bus_idx ]

            comps_before_slack_ur_ui_Idx =
                get_components_ur_ui_idx(
                    comp_collection_before_slack )

            comps_after_slack_ur_ui_Idx =
                get_components_ur_ui_idx(
                    comp_collection_after_slack )
            
            comps_before_slack_dims =
                get_components_state_algebraic_vars_dim(
                    comp_collection_before_slack )

            sum_slack_and_comp_before_dims =
                slack_comp_dims +
                sum(comps_before_slack_dims)
            
            comps_after_slack_dims =
                get_components_state_algebraic_vars_dim(
                    comp_collection_after_slack )

            comps_before_slack_ur_ui_Idx_in_state =
                get_ur_ui_Idx_in_state(
                    comps_before_slack_dims,
                    comps_before_slack_ur_ui_Idx;
                    counter = counter)
            
            new_counter = counter +
                sum_slack_and_comp_before_dims
            
            comps_sfter_slack_ur_ui_Idx_in_state =
                get_non_slack_ur_ui_Idx_in_state(
                    comps_after_slack_dims,
                    comps_after_slack_ur_ui_Idx;
                    counter = new_counter )
            
            non_slack_ur_ui_Idx_in_state =
                [ comps_before_slack_ur_ui_Idx_in_state ;
                  comps_sfter_slack_ur_ui_Idx_in_state ]
            return non_slack_ur_ui_Idx_in_state
        end

    elseif isa(comp_collection, Union{Array, Vector})

        slack_bus_idx =
            get_components_slack_Idx(
                comp_collection )

        if slack_bus_idx == 1

            slack_comp_dims =
                sum([ comp.dim
                      for comp in
                          comp_collection
                          if (comp.Bus_type == :Generator &&
                              comp.isa_slack == true) ])

            non_slack_comp_collection =
                [comp
                 for comp in
                     comp_collection
                     if (comp.Bus_type == :Generator &&
                         comp.isa_slack == false)||
                         (comp.Bus_type != :Generator)  ]

            no_slack_comps_ur_ui_Idx =
                get_components_no_slack_ur_ui_idx(
                    non_slack_comp_collection)

            no_slack_comps_dims =
                get_components_state_algebraic_vars_dim(
                    non_slack_comp_collection)

            new_counter = counter + slack_comp_dims

            return get_non_slack_ur_ui_Idx_in_state(
                no_slack_comps_dims,
                no_slack_comps_ur_ui_Idx;
                counter = new_counter )
        else

            slack_comp_dims =
                sum([ comp.dim
                      for comp in comp_collection
                          if (comp.Bus_type == :Generator &&
                              comp.isa_slack == true) ])

            non_slack_comp_collection =
                [comp
                 for comp in
                     comp_collection
                     if (comp.Bus_type == :Generator &&
                         comp.isa_slack == false)||
                         (comp.Bus_type != :Generator)  ]

            comp_collection_before_slack =
                [comp
                 for  comp in
                     non_slack_comp_collection
                     if  comp.Bus_num < slack_bus_idx ]

            comp_collection_after_slack =
                [comp
                 for  comp in
                     non_slack_comp_collection
                     if comp.Bus_num > slack_bus_idx ]

            comps_before_slack_ur_ui_Idx =
                get_components_ur_ui_idx(
                    comp_collection_before_slack )

            comps_after_slack_ur_ui_Idx =
                get_components_ur_ui_idx(
                    comp_collection_after_slack )
            
            comps_before_slack_dims =
                get_components_state_algebraic_vars_dim(
                    comp_collection_before_slack )

            sum_slack_and_comp_before_dims =
                slack_comp_dims +
                sum(comps_before_slack_dims)
            
            comps_after_slack_dims =
                get_components_state_algebraic_vars_dim(
                    comp_collection_after_slack )

            comps_before_slack_ur_ui_Idx_in_state =
                get_ur_ui_Idx_in_state(
                    comps_before_slack_dims,
                    comps_before_slack_ur_ui_Idx;
                    counter = counter)
            
            new_counter =
                counter +
                sum_slack_and_comp_before_dims
            
            comps_sfter_slack_ur_ui_Idx_in_state =
                get_non_slack_ur_ui_Idx_in_state(
                    comps_after_slack_dims,
                    comps_after_slack_ur_ui_Idx;
                    counter = new_counter )
            
            non_slack_ur_ui_Idx_in_state =
                [ comps_before_slack_ur_ui_Idx_in_state ;
                  comps_sfter_slack_ur_ui_Idx_in_state ]
            
            return non_slack_ur_ui_Idx_in_state
        end

    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        slack_bus_idx =
            get_components_slack_Idx(
                comp_collection_values )

        if slack_bus_idx == 1

            slack_comp_dims =
                sum([ comp.dim
                      for comp in
                          comp_collection_values
                          if (comp.Bus_type == :Generator &&
                              comp.isa_slack == true) ])

            non_slack_comp_collection =
                [comp
                 for comp in comp_collection_values
                     if (comp.Bus_type == :Generator &&
                         comp.isa_slack == false)||
                         (comp.Bus_type != :Generator)  ]

            no_slack_comps_ur_ui_Idx =
                get_components_no_slack_ur_ui_idx(
                    non_slack_comp_collection)

            no_slack_comps_dims =
                get_components_state_algebraic_vars_dim(
                    non_slack_comp_collection)

            new_counter =
                counter +
                slack_comp_dims

            return get_non_slack_ur_ui_Idx_in_state(
                no_slack_comps_dims,
                no_slack_comps_ur_ui_Idx;
                counter = new_counter )
        else

            slack_comp_dims =
                sum([ comp.dim
                      for comp in
                          comp_collection_values
                          if (comp.Bus_type == :Generator &&
                              comp.isa_slack == true) ])

            non_slack_comp_collection =
                [comp for comp in
                     comp_collection_values
                     if (comp.Bus_type == :Generator &&
                         comp.isa_slack == false)||
                         (comp.Bus_type != :Generator)  ]

            comp_collection_before_slack =
                [comp
                 for  comp in
                     non_slack_comp_collection
                     if  comp.Bus_num < slack_bus_idx ]
            
            comp_collection_after_slack =
                [comp
                 for  comp in
                     non_slack_comp_collection
                     if comp.Bus_num > slack_bus_idx ]

            comps_before_slack_ur_ui_Idx =
                get_components_ur_ui_idx(
                    comp_collection_before_slack )

            comps_after_slack_ur_ui_Idx =
                get_components_ur_ui_idx(
                    comp_collection_after_slack )
            
            comps_before_slack_dims =
                get_components_state_algebraic_vars_dim(
                    comp_collection_before_slack )

            sum_slack_and_comp_before_dims =
                slack_comp_dims +
                sum(comps_before_slack_dims)
            
            comps_after_slack_dims =
                get_components_state_algebraic_vars_dim(
                    comp_collection_after_slack )

            comps_before_slack_ur_ui_Idx_in_state =
                get_ur_ui_Idx_in_state(
                    comps_before_slack_dims,
                    comps_before_slack_ur_ui_Idx;
                    counter = counter)
            
            new_counter =
                counter +
                sum_slack_and_comp_before_dims
            
            comps_sfter_slack_ur_ui_Idx_in_state =
                get_non_slack_ur_ui_Idx_in_state(
                    comps_after_slack_dims,
                    comps_after_slack_ur_ui_Idx;
                    counter = new_counter )
            
            non_slack_ur_ui_Idx_in_state =
                [ comps_before_slack_ur_ui_Idx_in_state ;
                  comps_sfter_slack_ur_ui_Idx_in_state ]
            
            return non_slack_ur_ui_Idx_in_state
        end
 
    end        
end


# ############

function get_components_slack_ur_ui_idx(
    comp_collection)

    if isa(comp_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(comp_collection))

        return [
            get_component_ur_ui_idx(comp)
            for comp in comp_collection_values
                if (comp.Bus_type == :Generator &&
                    comp.isa_slack == true)  ]
        
    elseif isa(comp_collection, Union{Array, Vector})
        
        return [
            get_component_ur_ui_idx(comp)
            for comp in comp_collection
                if (comp.Bus_type == :Generator &&
                    comp.isa_slack == true) ]
    else
        
        comp_collection_values =
            collect(values(comp_collection))
        
        return [
            get_component_ur_ui_idx(comp)
            for comp in
                comp_collection_values
                if (comp.Bus_type == :Generator &&
                    comp.isa_slack == true) ]
    end
end


function get_components_slack_ur_ui_Idx_in_state(
    comp_collection; counter=0 )

    if isa(comp_collection,  OrderedDict)

        comp_collection_values =
            collect(values(comp_collection))

        slack_comp_collection =
            [ comp
              for comp in
                  comp_collection_values
                  if (comp.Bus_type == :Generator &&
                      comp.isa_slack == true) ]
        
        slack_bus_idx =
            get_components_slack_Idx(
                slack_comp_collection )

        comp_collection_before_slack =
            [comp
             for comp in
                 comp_collection_values
                 if ((comp.Bus_type == :Generator &&
                     comp.isa_slack == false)||
                     (comp.Bus_type != :Generator)) &&
                     comp.Bus_num < slack_bus_idx ]
        
        slack_comps_ur_ui_Idx =
            get_components_ur_ui_idx(
                slack_comp_collection)

        slack_comps_dims =
            get_components_state_algebraic_vars_dim(
                slack_comp_collection)

        slack_offset =
            comp_collection_before_slack == [] ? 0 :
            sum([comp.dim
                 for comp in
                     comp_collection_before_slack])
        
        # get_non_slack_ur_ui_Idx_in_state is used instead
        # of get_ur_ui_Idx_in_state because it takes
        # account of offest in the first item in a list
        
        return get_non_slack_ur_ui_Idx_in_state(
            slack_comps_dims,
            slack_comps_ur_ui_Idx;
            counter = slack_offset )
        
    elseif isa(comp_collection, Union{Array, Vector})

        slack_comp_collection =
            [ comp
              for comp in comp_collection
                  if (comp.Bus_type == :Generator &&
                      comp.isa_slack == true) ]
        
        slack_bus_idx =
            get_components_slack_Idx(
                slack_comp_collection )

        comp_collection_before_slack =
            [ comp
              for comp in
                  enumerate(comp_collection)
                  if ((comp.Bus_type == :Generator &&
                      comp.isa_slack == false)||
                      (comp.Bus_type != :Generator)) &&
                      comp.Bus_num < slack_bus_idx ]
        
        slack_comps_ur_ui_Idx =
            get_components_ur_ui_idx(
                slack_comp_collection)

        slack_comps_dims =
            get_components_state_algebraic_vars_dim(
                slack_comp_collection)

        slack_offset = comp_collection_before_slack == [] ?
            0 :
            sum([comp.dim
                 for comp in
                     comp_collection_before_slack])
        
        return get_non_slack_ur_ui_Idx_in_state(
            slack_comps_dims,
            slack_comps_ur_ui_Idx;
            counter = slack_offset )

    else
        
        comp_collection_values =
            collect(values(comp_collection))

        slack_comp_collection =
            [ comp
              for comp in
                  comp_collection_values
                  if (comp.Bus_type == :Generator &&
                      comp.isa_slack == true) ]
        
        slack_bus_idx =
            get_components_slack_Idx(
                slack_comp_collection )

        comp_collection_before_slack =
            [comp
             for comp in
                 comp_collection_values
                 if ((comp.Bus_type == :Generator &&
                     comp.isa_slack == false)||
                     (comp.Bus_type != :Generator)) &&
                     comp.Bus_num < slack_bus_idx ]
        
        slack_comps_ur_ui_Idx =
            get_components_ur_ui_idx(
                slack_comp_collection)

        slack_comps_dims =
            get_components_state_algebraic_vars_dim(
                slack_comp_collection)

        slack_offset =
            comp_collection_before_slack == [] ?
            0 :
            sum([comp.dim
                 for comp in
                     comp_collection_before_slack])
        
        # get_non_slack_ur_ui_Idx_in_state is used instead
        # of get_ur_ui_Idx_in_state because it takes
        # account of offest in the first item in a list
        
        return get_non_slack_ur_ui_Idx_in_state(
            slack_comps_dims,
            slack_comps_ur_ui_Idx;
            counter = slack_offset )    

    end        
end



function get_components_params_value_vectors(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_params_value,
            collect(comp_collection) )
    else
        return map(
            get_component_params_value,
            collect(values(comp_collection)))
    end
end


function get_components_params_dims(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        vec_components_params =
            map(get_component_params_value,
                collect(comp_collection) )
    else
        vec_components_params =
            map(get_component_params_value,
                collect(values(comp_collection)))
    end
    return length.(vec_components_params)
end


function get_components_control_sig(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_control_sig,
            collect(comp_collection) )
    else
        return map(
            get_component_control_sig,
            collect(values(comp_collection)))
    end    
end


function get_components_control_sig_dims(
    comp_collection)

    return length.(
        get_components_control_sig(
            comp_collection))
end


function get_components_control_sig_syms(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        return map(
            get_component_control_sig_syms,
            collect(comp_collection) )
    else
        return map(
            get_component_control_sig_syms,
            collect(values(comp_collection)))
    end    
end



function get_components_output_sig(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_output_sig,
            collect(comp_collection) )
    else
        
        return map(
            get_component_output_sig,
            collect(values(comp_collection)))
    end    
end


function get_components_output_sig_dims(
    comp_collection)

    return length.(
        get_components_output_sig(
            comp_collection) )
end


function get_components_cb_dyn_param_state_sw(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_cb_dyn_param_state_sw,
            collect(comp_collection) )
    else
        return map(
            get_component_cb_dyn_param_state_sw,
            collect(values(comp_collection)))
    end    
end



function get_components_cb_dyn_param_state_sw_dims(
    comp_collection)

    return length.(
        get_components_cb_dyn_param_state_sw(
            comp_collection) )
end



function get_components_dae_var(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return [
            map(
                get_component_dae_var,
                collect(comp_collection))...;]
    else
        
        return [
            map(get_component_dae_var,
                collect(values(comp_collection)))...;]
    end
end


function get_components_mass_matrix(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_mass_matrix,
            collect(comp_collection) )
    else
        
        return map(
            get_component_mass_matrix,
            collect(values(comp_collection)))
    end
end


# get_component_im_vars_syms_dim( comp )


function get_gens_im_vars_syms_dim(
    nodes_collection)
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_component_im_vars_syms_dim( comp )
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_component_im_vars_syms_dim( comp )
                 for comp in nodes_collection if
                     comp.Bus_type == :Generator ]
    else
        comp_collection_values =
            collect(values(
                nodes_collection))
        
        return [ get_component_im_vars_syms_dim( comp )
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]        
    end
    
end



function get_gens_im_mass_matrix(
    nodes_collection)
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_a_gen_im_mass_matrix( comp )
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_a_gen_im_mass_matrix( comp )
                 for comp in nodes_collection
                     if comp.Bus_type == :Generator ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ get_a_gen_im_mass_matrix( comp )
                 for comp in comp_collection_values
                     if comp.Bus_type == :Generator ]        
    end
    
end


function get_components_system_matrices(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        return map(
            get_component_system_matrices,
            collect(comp_collection) )
    else
        return map(
            get_component_system_matrices,
            collect(values(comp_collection)))
    end
end


function get_components_system_stab_matrices(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_system_stab_matrices,
            collect(comp_collection) )
    else
        
        return map(
            get_component_system_stab_matrices,
            collect(values(comp_collection)))
    end
end



function get_components_func_system_matrices(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_func_system_matrices,
            collect(comp_collection) )
    else
        return map(
            get_component_func_system_matrices,
            collect(values(comp_collection)))
    end
end


function get_components_func_system_stab_matrices(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_func_system_stab_matrices,
            collect(comp_collection) )
    else
        
        return map(
            get_component_func_system_stab_matrices,
            collect(values(comp_collection)))
    end
end


function get_components_dyn_func_by_idx(
    comp_collection; idx = 2)
    
    if isa(comp_collection, Union{Array, Vector})
        return map(
            (comp) ->
                get_component_dyn_func_by_idx(
                    comp; idx = idx),
            collect(comp_collection) )
    else
        
        return map(
            (comp) ->
                get_component_dyn_func_by_idx(
                    comp; idx = idx),
            collect(values(comp_collection)))
    end  
end


# function get_component_dyn_func_by_idx(component;idx = 2)
#     lens_comp_dyn_func=@optic _.func
#     return getall(component, lens_comp_dyn_func)[1][1]
# end


function get_components_dyn_func(
    comp_collection)
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_dyn_func,
            collect(comp_collection) )
    else
        return map(
            get_component_dyn_func,
            collect(values(comp_collection)))
    end    
end

function get_components_dyn_func_non_flat(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_dyn_func_non_flat,
            collect(comp_collection) )
    else
        return map(
            get_component_dyn_func_non_flat,
            collect(values(comp_collection)))
    end    
end

function get_components_cb_state_syms(
    comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        
        return map(
            get_component_cb_state_syms,
            collect(comp_collection) )
    else
        return map(
            get_component_cb_state_syms,
            collect(values(comp_collection)))
    end    
end



function get_components_cb_dyn_state_syms(
    comp_collection)
    if isa(comp_collection, Union{Array, Vector})
        return map(
            get_component_cb_dyn_state_syms,
            collect(comp_collection) )
    else
        return map(
            get_component_cb_dyn_state_syms,
            collect(values(comp_collection)))
    end    
end


# ------------------------------------------------------

function symbolsof(component)
    return component.syms
end

# ------------------------------------------------------


function syms_containing(nd, expr)
    
    if typeof(expr) == String
        
        return [
            s for s in
                get_network_vars_labels(nd)
                if occursin(expr, string(s))]
    else
        
        return [
            s for s in
                get_network_vars_labels(nd)
                if occursin(string(expr), string(s))]
    end
    
end


function idx_containing(nd, expr)
    
    if typeof(expr) == String
        
        return [
            i for (i, s) in
                enumerate(get_network_state_labels(nd))
                if occursin(expr, string(s))]
    else
        return [
            i for (i, s) in
                enumerate(get_network_state_labels(nd))
                if occursin(string(expr), string(s))]
    end
    
end


#-----------------------------------------------------



function make_case_buses_names(
    nodes )

    network_bus_names  =
        get_network_bus_names( nodes )
    
    non_gens_bus_names =
        get_non_generators_bus_names( nodes )

    gens_bus_names =
        get_generators_bus_names( nodes )
    
    Loads_bus_names =
        get_Loads_bus_names( nodes )
    
    Trans_bus_names =
        get_Trans_bus_names( nodes )

    return (;
            network_bus_names,
            non_gens_bus_names,
            gens_bus_names,
            Loads_bus_names,
            Trans_bus_names )
    
end

#-----------------------------------------------------

function generate_industrial_model_labels(
    ; nodes =  nodes,
    no_control_device =false )


    if no_control_device == false
        
        gens_nodes_collection =
            get_gens_nodes( nodes )

        network_bus_names, _, _, _,_ =
            make_case_buses_names( nodes ) 


        net_bus_volts_labels =
            generate_net_bus_volts_labels(
                network_bus_names)

        gens_nodes_algebraic_and_states_labels =
            generate_algebraic_and_states_labels(
                gens_nodes_collection )

        gens_nodes_pure_states_labels =
            generate_pure_states_labels(
                gens_nodes_collection )

        gens_nodes_stab_states_label =
            generate_stab_states_labels(
                gens_nodes_collection )

        return (; net_bus_volts_labels,
                gens_nodes_pure_states_labels,
                gens_nodes_stab_states_label,
                gens_nodes_algebraic_and_states_labels )
        
    else
        
        gens_nodes_collection =
            get_gens_nodes( nodes )

        network_bus_names, _, _, _,_ =
            make_case_buses_names( nodes ) 


        net_bus_volts_labels =
            generate_net_bus_volts_labels(
                network_bus_names)

        gens_nodes_algebraic_and_states_labels =
            generate_algebraic_and_states_labels(
                gens_nodes_collection;
                no_control_device = true )

        gens_nodes_pure_states_labels =
            generate_pure_states_labels(
                gens_nodes_collection;
                no_control_device = true )

        gens_nodes_stab_states_label =
            generate_stab_states_labels(
                gens_nodes_collection;
                no_control_device = true )

        return (; net_bus_volts_labels,
                gens_nodes_pure_states_labels,
                gens_nodes_stab_states_label,
                gens_nodes_algebraic_and_states_labels )
        
    end
            
end



function generate_industrial_model_sym_and_mass_matrix(
    ; nodes = nodes,
    no_control_device =false )

    if no_control_device == false
        
        net_bus_volts_labels, gens_nodes_pure_states_labels, _, _ = generate_industrial_model_labels(; nodes = nodes  )

        industrial_model_sym = vcat(
            gens_nodes_pure_states_labels,
            net_bus_volts_labels )


        industrial_model_mass_matrix = DAE_MassMatrix(
            length( gens_nodes_pure_states_labels ),
            length( net_bus_volts_labels  ))

        return (; industrial_model_sym,
                industrial_model_mass_matrix )
    else
        
        net_bus_volts_labels, gens_nodes_pure_states_labels, _, _ = generate_industrial_model_labels(
            ; nodes = nodes,
            no_control_device = true  )

        industrial_model_sym = vcat(
            gens_nodes_pure_states_labels,
            net_bus_volts_labels )

        industrial_model_mass_matrix = DAE_MassMatrix(
            length( gens_nodes_pure_states_labels ),
            length( net_bus_volts_labels  ))

        return (; industrial_model_sym,
                industrial_model_mass_matrix )        
    end
    
end



function generate_industrial_model_sym(
    ; nodes = nodes,
    no_control_device =false )


    if  no_control_device == false
        net_bus_volts_labels, gens_nodes_pure_states_labels, _, _ = generate_industrial_model_labels(; nodes = nodes  )

        industrial_model_sym = vcat(
            gens_nodes_pure_states_labels,
            net_bus_volts_labels )


        return industrial_model_sym
    else
        net_bus_volts_labels, gens_nodes_pure_states_labels, _, _ = generate_industrial_model_labels(
            ; nodes = nodes,
            no_control_device = true  )

        industrial_model_sym = vcat(
            gens_nodes_pure_states_labels,
            net_bus_volts_labels )

        return industrial_model_sym 

    end
    

end


function generate_industrial_model_mass_matrix(
    ; nodes = nodes,
    no_control_device =false )

    if no_control_device == false
        net_bus_volts_labels, gens_nodes_pure_states_labels, _, _ = generate_industrial_model_labels(; nodes = nodes  )

        industrial_model_mass_matrix  = DAE_MassMatrix(
            length( gens_nodes_pure_states_labels ),
            length( net_bus_volts_labels  ))

        return industrial_model_mass_matrix
    else
        net_bus_volts_labels, gens_nodes_pure_states_labels, _, _ = generate_industrial_model_labels(
            ; nodes = nodes, no_control_device = true  )

        industrial_model_mass_matrix  = DAE_MassMatrix(
            length( gens_nodes_pure_states_labels ),
            length( net_bus_volts_labels  ))

        return industrial_model_mass_matrix

    end
            
end


#-----------------------------------------------------
#-----------------------------------------------------


function generate_im_model_labels(
    ; nodes =  nodes )
        
    gens_nodes_collection =
        get_gens_nodes( nodes )

    network_bus_names, _, _, _,_ =
        make_case_buses_names( nodes ) 


    net_bus_volts_labels =
        generate_net_bus_volts_labels(
            network_bus_names)

    gens_nodes_pure_states_labels =
        generate_pure_states_labels(
            gens_nodes_collection )

    gens_nodes_im_vars_labels =
        generate_im_vars_labels(
            gens_nodes_collection )
    
    gens_nodes_im_algebraic_vars_labels =
        generate_im_algebraic_vars_labels(
            gens_nodes_collection )

    return (; net_bus_volts_labels,
            gens_nodes_pure_states_labels,
            gens_nodes_im_algebraic_vars_labels,
            gens_nodes_im_vars_labels )
            
end



function generate_im_sym_and_mass_matrix(
    ; nodes = nodes )

    (net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes = nodes  )

    im_model_sym = vcat(
        gens_nodes_im_vars_labels,
        net_bus_volts_labels )


    # state_size = length(gens_nodes_pure_states_labels)

    # algebraic_size = sum(
    #     length.(
    #         [ gens_nodes_im_algebraic_vars_labels,
    #           net_bus_volts_labels ]))
    
    # im_model_mass_matrix = DAE_MassMatrix(
    #     state_size,
    #     algebraic_size )

    gens_im_vars_dim =
        get_gens_im_vars_syms_dim(nodes)
    
    gens_im_vars_size =
        sum( gens_im_vars_dim  )

    nodes_pf_vars_size =
        length(net_bus_volts_labels )
    
    system_eq_size =
        gens_im_vars_size + nodes_pf_vars_size
    
    vec_gens_im_mass_matrix =
        get_gens_im_mass_matrix(nodes)

    vec_nodes_pf_vars_mass_matrix =
        [DAE_MassMatrix(0, nodes_pf_vars_size)]

    _, _, gens_im_vars_Idx =
        create_size_offset_Idx(
        gens_im_vars_dim )

    # The Idx for pf_var starts after gens_im_vars_Idx
    _,_, nodes_pf_vars_Idx =
        create_size_offset_Idx(
            [nodes_pf_vars_size];
            counter = gens_im_vars_size)

    im_model_mass_matrix = sparse(1.0I, system_eq_size,
                         system_eq_size)
    

    
    for (i, gen_mass_matrix) in
        enumerate( vec_gens_im_mass_matrix  )

        gen_state_Idx =
            gens_im_vars_Idx[i]
        
        copyto!(
            @view( im_model_mass_matrix[
                gen_state_Idx, gen_state_Idx]),
            gen_mass_matrix)
        
    end

    for (i, nodes_pf_vars_mass_matrix ) in
        enumerate( vec_nodes_pf_vars_mass_matrix )

        pf_vars_Idx = nodes_pf_vars_Idx[i]
        
        copyto!(
            @view( im_model_mass_matrix[
                pf_vars_Idx, pf_vars_Idx] ),
            nodes_pf_vars_mass_matrix )
    end
    
    return (; im_model_sym,
            im_model_mass_matrix )
    
end


function generate_im_ode_and_pf_sym(; nodes = nodes  )

    net_bus_volts_labels,_,_,gens_nodes_im_vars_labels =
        generate_im_model_labels(
            ; nodes = nodes  )

    return (;gens_nodes_im_vars_labels,
            net_bus_volts_labels )

end


function generate_im_sym(; nodes = nodes  )

    net_bus_volts_labels,_,_,gens_nodes_im_vars_labels =
        generate_im_model_labels(
            ; nodes = nodes  )

    im_model_sym = vcat(
        gens_nodes_im_vars_labels,
        net_bus_volts_labels )


    return im_model_sym    

end


function generate_im_mass_matrix(; nodes = nodes )

    (net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes = nodes  )


    # state_size = length(gens_nodes_pure_states_labels)

    # algebraic_size = sum(
    #     length.(
    #         [ gens_nodes_im_algebraic_vars_labels,
    #           net_bus_volts_labels ]))
    
    # im_model_mass_matrix = DAE_MassMatrix(
    #     state_size,
    #     algebraic_size )
    
    # return im_model_mass_matrix


    net_bus_volts_labels, _, _, _ =
        generate_im_model_labels(; nodes = nodes )
    
    gens_im_vars_dim =
        get_gens_im_vars_syms_dim(nodes)

    nodes_pf_vars_size =
        length( net_bus_volts_labels )
    
    gens_im_vars_size =
        sum( gens_im_vars_dim  )
    
    system_eq_size =
        gens_im_vars_size + nodes_pf_vars_size

    vec_nodes_pf_vars_mass_matrix =
        [DAE_MassMatrix(0, nodes_pf_vars_size)]
    
    vec_gens_im_mass_matrix =
        get_gens_im_mass_matrix(nodes)

    _, _, gens_im_vars_Idx =
        create_size_offset_Idx(
        gens_im_vars_dim )

    _,_, nodes_pf_vars_Idx =
        create_size_offset_Idx(
            [nodes_pf_vars_size];
            counter = gens_im_vars_size)

    im_model_mass_matrix =
        sparse(1.0I, system_eq_size,
               system_eq_size)
    

    
    for (i, gen_mass_matrix) in
        enumerate( vec_gens_im_mass_matrix  )

        gen_state_Idx = gens_im_vars_Idx[i]
        
        copyto!(
            @view( im_model_mass_matrix[
                gen_state_Idx, gen_state_Idx]),
            gen_mass_matrix)
        
    end

    for (i, nodes_pf_vars_mass_matrix ) in
        enumerate( vec_nodes_pf_vars_mass_matrix )

        pf_vars_Idx = nodes_pf_vars_Idx[i]
        
        copyto!(
            @view( im_model_mass_matrix[
                pf_vars_Idx, pf_vars_Idx] ),
            nodes_pf_vars_mass_matrix )
    end
    
    return  im_model_mass_matrix 
    
end



function generate_im_ode_and_pf_mass_matrix(
    ; nodes = nodes )

    (net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes = nodes  )


    net_bus_volts_labels, _, _, _ =
        generate_im_model_labels(
            ; nodes = nodes )
    
    gens_im_vars_dim =
        get_gens_im_vars_syms_dim(nodes)

    nodes_pf_vars_size =
        length( net_bus_volts_labels )
    
    gens_im_vars_size =
        sum( gens_im_vars_dim  )
    
    system_ode_eq_size =
        gens_im_vars_size 

    
    system_pf_eq_size =
         nodes_pf_vars_size
    
    vec_nodes_pf_vars_mass_matrix =
        [DAE_MassMatrix(0, nodes_pf_vars_size)]
    
    vec_gens_im_mass_matrix =
        get_gens_im_mass_matrix(nodes)

    _, _, gens_im_vars_Idx =
        create_size_offset_Idx(
        gens_im_vars_dim )

    _,_, nodes_pf_vars_Idx =
        create_size_offset_Idx(
            [nodes_pf_vars_size];
            counter = 0)

    im_model_ode_mass_matrix =
        sparse(1.0I, system_ode_eq_size,
               system_ode_eq_size)


    im_model_pf_mass_matrix =
        sparse(1.0I, system_pf_eq_size,
               system_pf_eq_size)
    

    
    for (i, gen_mass_matrix) in
        enumerate( vec_gens_im_mass_matrix  )

        gen_state_Idx = gens_im_vars_Idx[i]
        
        copyto!(
            @view( im_model_ode_mass_matrix[
                gen_state_Idx, gen_state_Idx]),
            gen_mass_matrix)
        
    end

    for (i, nodes_pf_vars_mass_matrix ) in
        enumerate( vec_nodes_pf_vars_mass_matrix )

        pf_vars_Idx = nodes_pf_vars_Idx[i]
        
        copyto!(
            @view( im_model_pf_mass_matrix[
                pf_vars_Idx, pf_vars_Idx] ),
            nodes_pf_vars_mass_matrix )
    end
    
    return  (; im_model_ode_mass_matrix,
             im_model_pf_mass_matrix)
    
end


#-----------------------------------------------------
#----------------------------------------------------


function generate_im_model_labels(
    network_bus_names,
    gens_nodes_collection )
    
    net_bus_volts_labels =
        generate_net_bus_volts_labels(
            network_bus_names)

    gens_nodes_pure_states_labels =
        generate_pure_states_labels(
            gens_nodes_collection )
    
    gens_nodes_im_vars_labels =
        generate_im_vars_labels(
            gens_nodes_collection )

    gens_nodes_im_algebraic_vars_labels  =
        generate_im_algebraic_vars_labels(
            gens_nodes_collection )

    return (; net_bus_volts_labels,
            gens_nodes_pure_states_labels,
            gens_nodes_im_algebraic_vars_labels ,
            gens_nodes_im_vars_labels )
    
end


#-----------------------------------------------------
#-----------------------------------------------------


function generate_industrial_model_labels(
    network_bus_names,
    gens_nodes_collection
    ; no_control_device =false )

    if no_control_device == false
        
        net_bus_volts_labels =
            generate_net_bus_volts_labels(
                network_bus_names)

        gens_nodes_algebraic_and_states_labels =
            generate_algebraic_and_states_labels(
                gens_nodes_collection )

        gens_nodes_pure_states_labels =
            generate_pure_states_labels(
                gens_nodes_collection )

        gens_nodes_stab_states_label  =
            generate_stab_states_labels(
                gens_nodes_collection )

        return (; net_bus_volts_labels,
                gens_nodes_pure_states_labels,
                gens_nodes_stab_states_label,
                gens_nodes_algebraic_and_states_labels )
    else
        
        net_bus_volts_labels =
            generate_net_bus_volts_labels(
                network_bus_names)

        gens_nodes_algebraic_and_states_labels =
            generate_algebraic_and_states_labels(
                gens_nodes_collection
                ; no_control_device =true  )

        gens_nodes_pure_states_labels =
            generate_pure_states_labels(
                gens_nodes_collection;
                no_control_device =true  )

        gens_nodes_stab_states_label =
            generate_stab_states_labels(
                gens_nodes_collection;
                no_control_device =true  )

        return (; net_bus_volts_labels,
                gens_nodes_pure_states_labels,
                gens_nodes_stab_states_label,
                gens_nodes_algebraic_and_states_labels )
    end
    
end



#-----------------------------------------------------
#----------------------------------------------------


function generate_components_algebraic_labels(
    comp_collection;
    no_control_device = false )
    comps_algb_labels = Symbol[]
    
    if no_control_device == false
        
        for (key, a_comp) in  pairs(comp_collection)
            
            comp_algb_vars_syms =
                get_component_algebraic_vars_syms(a_comp)
            append!( comps_algb_labels,
                     [Symbol(key, "_", sym )
                      for sym in comp_algb_vars_syms])
        end
        
        return comps_algb_labels
        
    else
        for (key, a_comp) in  pairs(comp_collection)
            comp_algb_vars_syms =
                a_comp.Bus_type == :Generator ?
                get_component_algebraic_vars_syms(
                    a_comp.Gen) :
                get_component_algebraic_vars_syms(a_comp)
            append!( comps_algb_labels,
                     [Symbol(key, "_", sym )
                      for sym in comp_algb_vars_syms])
        end
        return comps_algb_labels        
    end    
    
end


 
function generate_algebraic_and_states_labels(
    comp_collection;
    no_control_device = false )
    comps_states_labels = Symbol[]
    
    if no_control_device == false
        
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_state_vars_syms =
                get_component_state_vars_syms(a_comp)
            
            append!( comps_states_labels,
                     [Symbol(a_comp.name, "_", sym )
                      for sym in comp_state_vars_syms])
        end
        
        return comps_states_labels
        
    else
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_state_vars_syms =
                a_comp.Bus_type == :Generator ?
                get_component_state_vars_syms(a_comp.Gen) :
                get_component_state_vars_syms(a_comp)
            
            append!( comps_states_labels,
                     [Symbol(a_comp.name, "_", sym )
                      for sym in comp_state_vars_syms])
            
        end
        return comps_states_labels        
    end
    
end


function generate_pure_states_labels(
    comp_collection;
    no_control_device = false )
    comps_states_labels = Symbol[]
    
    if no_control_device == false
        
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_state_vars_syms =
                a_comp.state_vars_syms
            
            append!( comps_states_labels,
                     [Symbol(a_comp.name, "_", sym )
                      for sym in comp_state_vars_syms])
        end
        
        return comps_states_labels
        
    else
        
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_state_vars_syms =
                a_comp.Gen.state_vars_syms
            
            append!( comps_states_labels,
                     [Symbol(a_comp.name, "_", sym )
                      for sym in comp_state_vars_syms])
        end
        
        return comps_states_labels        
    end
    
end



function generate_im_vars_labels(
    comp_collection;
    no_control_device = false )
    
    comps_im_vars_labels = Symbol[]
    
    if no_control_device == false
        
        for (key, a_comp) in  pairs(comp_collection)
            
            comp_im_vars_syms = a_comp.im_vars_syms
            
            append!( comps_im_vars_labels,
                     [Symbol(a_comp.name, "_", sym )
                      for sym in comp_im_vars_syms])
            
        end
        
        return comps_im_vars_labels
    else
        
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_state_vars_syms =
                a_comp.Gen.state_vars_syms
            
            append!( comps_states_labels,
                     [Symbol(a_comp.name, "_", sym )
                      for sym in comp_state_vars_syms])
        end
        
        return comps_states_labels        
    end
    
end


function generate_im_algebraic_vars_labels(
    comp_collection )

    comps_im_algebraic_vars_labels = Symbol[]

    for (key, a_comp) in
        pairs(comp_collection)
        
        comp_im_algebraic_vars_syms =
            a_comp.im_algebraic_vars_syms
        
        append!(
            comps_im_algebraic_vars_labels,
            [Symbol(a_comp.name, "_", sym )
             for sym in
                 comp_im_algebraic_vars_syms])
    end

    return comps_im_algebraic_vars_labels    
    
end



# im_algebraic_vars_syms


function generate_stab_states_labels(
    comp_collection;
    no_control_device = false )
    
    comps_states_labels = Symbol[]
    
    if no_control_device == false
        
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_state_vars_syms =
                a_comp.stab_state_vars_syms
            
            append!(
                comps_states_labels,
                [Symbol(a_comp.name, "_", sym )
                 for sym in comp_state_vars_syms])
            
        end
        
        return comps_states_labels
        
    else
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_state_vars_syms =
                a_comp.Gen.stab_state_vars_syms
            
            append!( comps_states_labels,
                     [Symbol(a_comp.name, "_", sym )
                      for sym in comp_state_vars_syms])
        end
        
        return comps_states_labels        
    end
    
end



function generate_components_states_labels(
    comp_collection;
    no_control_device = false )

    comps_states_labels = Symbol[]
    
    if no_control_device == false
        
        for (key, a_comp) in
            pairs(comp_collection)
            comp_state_vars_syms =
                get_component_state_vars_syms(
                    a_comp)
            
            append!(
                comps_states_labels,
                [Symbol(key, "_", sym )
                 for sym in comp_state_vars_syms])
        end
        
        return comps_states_labels
        
    else

        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_state_vars_syms =
                a_comp.Bus_type == :Generator ?
                get_component_state_vars_syms(a_comp.Gen) :
                get_component_state_vars_syms(a_comp)
            
            append!(
                comps_states_labels,
                [Symbol(key, "_", sym )
                 for sym in comp_state_vars_syms])
        end
        
        return comps_states_labels        
        
    end
    
end


function generate_components_stab_states_labels(
    comp_collection;
    no_control_device = false)
    
    comps_stab_states_labels = Symbol[]
    
    if no_control_device == false
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_stab_state_vars_syms =
                get_component_stab_state_vars_syms(
                    a_comp)
            
            append!(
                comps_stab_states_labels,
                [Symbol(key, "_", sym )
                 for sym in
                     comp_stab_state_vars_syms])
        end
        
        return comps_stab_states_labels
        
    else
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_stab_state_vars_syms =
                get_component_stab_state_vars_syms(
                    a_comp.Gen)
            
            append!(
                comps_stab_states_labels,
                [Symbol(key, "_", sym )
                 for sym in
                     comp_stab_state_vars_syms])
        end
        
        return comps_stab_states_labels        
    end
    
end


function generate_components_im_vars_syms_labels(
    comp_collection;
    no_control_device = false)
    
    comps_im_vars_syms_labels = Symbol[]
    
    if no_control_device == false
        
        for (key, a_comp) in
            pairs(comp_collection)
            comp_im_vars_syms =
                get_component_im_vars_syms(a_comp)
            
            append!(
                comps_im_vars_syms_labels,
                [Symbol(key, "_", sym )
                 for sym in comp_im_vars_syms])
        end
        
        return comps_im_vars_syms_labels
    else
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_im_vars_syms =
                get_component_state_vars_syms(
                    a_comp.Gen)
            
            append!(
                comps_im_vars_syms_labels,
                [Symbol(key, "_", sym )
                 for sym in comp_im_vars_syms])
        end
        
        return comps_im_vars_syms_labels        
    end
    
end



function generate_components_im_algebraic_vars_syms_labels(
    comp_collection )
    
    comps_im_algebraic_vars_syms_labels = Symbol[]
    
    for (key, a_comp) in
        pairs(comp_collection)
        
        comp_im_algebraic_vars_syms =
            get_component_im_algebraic_vars_syms(
                a_comp)
        
        append!(
            comps_im_algebraic_vars_syms_labels,
            [Symbol(key, "_", sym )
             for sym in
                 comp_im_algebraic_vars_syms])

    end

    return comps_im_algebraic_vars_syms_labels 
end




function generate_components_vars_labels(
    comp_collection; no_control_device = false )
    
    comps_vars_labels = Symbol[]
    
    if no_control_device == false
        
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_vars_syms =
                get_component_syms(a_comp)
            
            append!( comps_vars_labels,
                     [Symbol(key, "_", sym )
                      for sym in comp_vars_syms])
        end
        
        return comps_vars_labels
    else
        
        for (key, a_comp) in
            pairs(comp_collection)
            
            comp_vars_syms =
                a_comp.Bus_type == :Generator ?
                get_component_syms(a_comp.Gen) :
                get_component_syms(a_comp)
            
            append!(
                comps_vars_labels,
                [Symbol(key, "_", sym )
                 for sym in comp_vars_syms])
        end
        
        return comps_vars_labels        
    end
    
end



# function generate_components_vars_labels(comp_collection; no_control_device = false )
#     comps_vars_labels = Symbol[]
#     if no_control_device == false
        
#         for (key, a_comp) in  pairs(comp_collection)
#             comp_vars_syms =  get_component_syms(a_comp)
#             append!( comps_vars_labels, [Symbol(key, "_", sym ) for sym in comp_vars_syms])

#         end
        
#         return comps_vars_labels
#     else
        
#         for (key, a_comp) in  pairs(comp_collection)
            
#             comp_vars_syms = a_comp.Bus_type == :Generator ? get_component_syms(a_comp.Gen) : get_component_syms(a_comp)
#             append!( comps_vars_labels, [Symbol(key, "_", sym ) for sym in comp_vars_syms])
#         end
#         return comps_vars_labels

#     end
    
# end


function get_network_algebraic_labels(
    nd ;  no_control_device = false)
    
    if no_control_device == false
        
        nodes_algebraic_labels =
            generate_components_algebraic_labels(nd.nodes)
        
        edges_algebraic_labels =
            generate_components_algebraic_labels(nd.edges)

        return vcat(nodes_algebraic_labels,
                    edges_algebraic_labels)
    else
        return  nodes_algebraic_labels =
            generate_components_algebraic_labels(
                nd.nodes; no_control_device = true)
    end
    
end


function get_network_state_labels(
    nd;  no_control_device = false )
    
    if no_control_device == false

        nodes_states_labels =
            generate_components_states_labels(nd.nodes)
        
        edges_states_labels =
            generate_components_states_labels(nd.edges)
    
        return vcat(nodes_states_labels,
                    edges_states_labels)
    else

        return nodes_states_labels =
            generate_components_states_labels(
                nd.nodes; no_control_device = true)
        
    end
    
end


function get_network_im_vars_labels(
    nd;  no_control_device = false )
    
    if no_control_device == false
            
        return generate_components_im_vars_syms_labels(
            nd.nodes)
    else

        return generate_components_states_labels(
            nd.nodes; no_control_device = true)
        
    end
    
end


function get_network_im_vars_labels(nd  )
    
    return generate_components_im_algebraic_vars_syms(
        nd.nodes)
    
end


function get_network_vars_labels(
    nd; no_control_device = false )

    if no_control_device == false
        
        nodes_vars_labels =
            generate_components_vars_labels(
                nd.nodes)

        edges_vars_labels =
            generate_components_vars_labels(
                nd.edges)
    
        return vcat(nodes_vars_labels,
                    edges_vars_labels)
    else
        
        return nodes_vars_labels =
            generate_components_vars_labels(
                nd.nodes; no_control_device = true )
    end
                    
end



function algebraic_syms_containing(nd, expr)
    
    if typeof(expr) == String
        
        return [
            s for s in
                get_network_algebraic_labels(nd)
                if occursin(expr, string(s))]
    else
        return [
            s for s in
                get_network_algebraic_labels(nd)
                if occursin(string(expr), string(s))]
    end
    
end


function state_syms_containing(nd, expr)
    
    if typeof(expr) == String
        
        return [
            s for s in
                get_network_state_labels(nd)
                if occursin(expr, string(s))]
    else
        return [
            s for s in
                get_network_state_labels(nd)
                if occursin(string(expr), string(s))]
    end
    
end


function vars_syms_containing(nd, expr)
    if typeof(expr) == String
        return [
            s for s in
                get_network_vars_labels(nd)
                if occursin(expr, string(s))]
    else
        return [
            s for s in
                get_network_vars_labels(nd)
                if occursin(string(expr), string(s))]
    end
    
end


function algebraic_idx_containing(nd, expr)
    if typeof(expr) == String
        
        return [
            i for (i, s) in
                enumerate(
                    get_network_algebraic_labels(nd))
                if occursin(expr, string(s))]
    else
        
        return [
            i for (i, s) in
                enumerate(
                    get_network_algebraic_labels(nd))
                if occursin(string(expr), string(s))]
    end
    
end


function state_idx_containing(nd, expr)
    
    if typeof(expr) == String
        
        return [
            i for (i, s) in
                enumerate(get_network_state_labels(nd))
                if occursin(expr, string(s))]
    else
        return [
            i for (i, s) in
                enumerate(get_network_state_labels(nd))
                if occursin(string(expr), string(s))]
    end
    
end


function vars_idx_containing(nd, expr)

    if typeof(expr) == String

        return [
            i for (i, s) in
                enumerate(get_network_vars_labels(nd))
                if occursin(expr, string(s))]
    else

        return [
            i for (i, s) in
                enumerate(get_network_vars_labels(nd))
                if occursin(string(expr), string(s))]
    end
    
end

