# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

########################################################
# ------------------------------------------------------
#  Network Idx functions
# ------------------------------------------------------
########################################################

function get_nodes_u_Idx_in_ranges( nodes_u_Idx )
    return UnitRange[ first(idx):last(idx) for idx in nodes_u_Idx]
end


"""
my_array = Vector{Symbol}[[:a, :bc, :d], [:ef, :k, :dc], [:de, :bc]]

my_vec = [:a, :bc, :d]

dims = length.(my_array)
offset = create_offsets(dims)
Idx = create_idxs(offset, dims)

findall(x -> x == Symbol("d"), my_array )

findall(x -> x == Symbol("d"), [my_array...;])

function findfirst(testf::Function, A)
    for (i, a) in pairs(A)
        testf(a) && return i
    end
    return nothing
end

function nd_findfirst(testf::Function, A)
    if typeof(A) == Vector{Symbol}
        return findfirst(testf, A)
    else
        for (j, vec) in pairs(A)
           for (i, a) in pairs(vec)
              testf(a) && return (j,i)
           end
        end
    end

    return nothing
end

nd_findfirst(x -> x == Symbol(:d), my_array )

nd_findfirst(x -> x == Symbol(:d), my_vec )

t_idx = nd_findfirst(x -> x == Symbol(:d), my_array )

my_array[CartesianIndex(t_idx)]


"""


"""
```Julia
nd_findfirst(testf::Function, A)
```
It return idx of an array or 2 level array of arrays.

"""
function nd_findfirst(testf::Function, A)
    if typeof(A) == Vector{Symbol}
        return findfirst(testf, A)
    else
        for (j, vec) in pairs(A)
           for (i, a) in pairs(vec)
              testf(a) && return (j,i)
           end
        end
    end

    return nothing
end

"""
```Julia
node_params_idxs(sd, node_name, para_symbol)
```
"""
function node_params_idxs(
    sd::NetworkData,
    node::Union{String,Symbol},
    para::Union{String,Symbol})
    
    list_nodes_idx   = keys(collect(sd.nodes))
    list_nodes_names = keys(sd.nodes)
    vec_nodes_params_syms = get_components_params_sym(
        sd.nodes) # get_nodes_param_syms(sd)
    
    node_idx = findfirst(
        node_name -> node_name == string(node),
        collect(list_nodes_names))
    
    list_node_params = vec_nodes_params_syms[node_idx]

    # node_para_idx = findfirst(
    #     node_para -> node_para == Symbol(para),
    #     collect(list_node_params))
    
    node_para_idx = nd_findfirst(
        node_para -> node_para == Symbol(para),
        collect(list_node_params))
    
   return node_idx, node_para_idx
end



"""
```Julia
component_params_idxs(comp_collection, comp_name, para_symbol)
```
Returns the index of a parameter of a comp. This function is important for changing the
parameter of a comp in  network data since the `parameters` tuple passed to `ODEProblem` could be deeply nested.

```
prob = ODEProblem(nd!, initial_state, timespan, (nd.params_nodes, nd.params_edges) )
```
"""
function component_params_idxs(
    comp_collection, comp::Union{String,Symbol},
    para::Union{String,Symbol})
    
    list_comps_idx   = keys(collect(comp_collection))
    
    list_comps_names = keys(comp_collection)
    
    vec_comps_params_syms = get_components_params_sym(
        comp_collection)
    
    comp_idx = findfirst(
        comp_name -> comp_name == string(comp),
        collect(list_comps_names))
    
    list_comp_params = vec_comps_params_syms[comp_idx]

    # comp_para_idx = findfirst(
    #     comp_para -> comp_para == Symbol(para),
    #     collect(list_comp_params))
    
    comp_para_idx = nd_findfirst(
        comp_para -> comp_para == Symbol(para),
        collect(list_comp_params))
    
   return comp_idx, comp_para_idx
end


function component_get_a_param_idxs(
    comp_collection, comp::Union{String,Symbol},
    para::Union{String,Symbol})
    
    list_comps_idx   = keys(collect(comp_collection))
    
    list_comps_names = keys(comp_collection)
    
    vec_comps_params_syms = get_components_params_sym(
        comp_collection)
    
    comp_idx = findfirst(
        comp_name -> comp_name == string(comp),
        collect(list_comps_names))
    
    list_comp_params = vec_comps_params_syms[comp_idx]

    # comp_para_idx = findfirst(
    #     comp_para -> comp_para == Symbol(para),
    #     collect(list_comp_params))

    comp_para_idx = nd_findfirst(
        comp_para -> comp_para == Symbol(para),
        collect(list_comp_params))
    
   return comp_idx, comp_para_idx
end



function component_get_control_sig_idxs(
    comp_collection,
    comp::Union{String,Symbol},
    control::Union{String,Symbol})
    
    list_comps_idx   = keys(collect(comp_collection))
    
    list_comps_names = keys(comp_collection)
    
    vec_comps_control_sig_syms =
        get_components_control_sig_syms(comp_collection)
    
    comp_idx = findfirst(
        comp_name -> comp_name == string(comp),
        collect(list_comps_names))
    
    list_comps_control_sig =
        vec_comps_control_sig_syms[comp_idx]
    
    # comps_control_sig_idx = findfirst(
    #     control_sig -> control_sig == Symbol(control),
    #     collect(list_comp_control_sig))
    
    comps_control_sig_idx = nd_findfirst(
        control_sig -> control_sig == Symbol(control),
        collect(list_comp_control_sig))
    
   return comp_idx, comp_control_sig_idx
end


function is_param_vector_type_1d( comp_para )
    
    return (typeof(comp_para) == Vector{Float64})|| (typeof(comp_para) == Vector{ComplexF64}) || (typeof(comp_para) == Vector{Int64}) || (typeof(comp_para) == Vector{Symbol}) || (typeof(comp_para) == Vector{ Tuple{Float64, Float64, Float64} }) || (typeof(comp_para) ==  Vector{Union{Float64, Tuple{Float64, Float64, Float64}}} ) || ( typeof(comp_para) == Vector{Union{ComplexF64, Float64, Tuple{Float64, Float64, Float64}}})

end



function get_a_param_value_in_a_component(
    comp_para, para_Idx)

    if is_param_vector_type_1d( comp_para )

        return comp_para[para_Idx]
    else
        
        return comp_para[para_Idx[1]][para_Idx[2]]
    end

    return nothing
    
end



function get_component_params_value_in_idxs!(
    components_para, comp_Idx, para_Idx)

    comp_para = components_para[comp_Idx]

    if is_param_vector_type_1d( comp_para )

        return comp_para[para_Idx]
    else
        
        return comp_para[para_Idx[1]][para_Idx[2]]
    end

    return nothing
    
end


function set_component_params_value_in_idxs!(components_para, comp_Idx, para_Idx, para_value)

    comp_para = components_para[comp_Idx]

    if is_param_vector_type_1d( comp_para )

        comp_para[para_Idx] = para_value
    else
        comp_para[para_Idx[1]][para_Idx[2]] = para_value
    end

    return nothing
    
end



"""
This function returns indices of state_algebraic variables of a node in the node syms list.

It is meant to extract the indices in sol.
"""
function get_a_node_state_algb_vars_indices_in_syms(; node_syms_labels = node_syms_labels, bus_name = bus_name, vars = [:ω, :δ])

    vars_indices = []

    for a_var in vars

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   = findfirst(comp_name -> comp_name == bus_var_label,  node_syms_labels )

        push!(vars_indices, (a_var, bus_var_Idx) )
    end

    return Tuple(vars_indices)
end




"""
This function returns indices of state_algebraic variables of a node in the systems.

It is meant to extract the indices in sol.
"""
function get_a_node_state_algb_vars_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = bus_name, vars = [:ω, :δ])

    vars_indices = []

    for a_var in vars

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   = findfirst(comp_name -> comp_name == bus_var_label,  network_vars_labels)

        push!(vars_indices, bus_var_Idx)
    end

    return Tuple(vars_indices)
end


function get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = bus_name, vars = [:δ, :ed_dash, :eq_dash])

    vars_syms = []

    for a_var in vars

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   = findfirst(comp_name -> comp_name == bus_var_label,  network_vars_labels)
        

        push!(vars_syms, network_vars_labels[bus_var_Idx])
    end

    return vars_syms
end



function get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = bus_name, vars = [:δ, :ed_dash, :eq_dash])

    vars_indices = []

    for a_var in vars

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   = findfirst(comp_name -> comp_name == bus_var_label,  network_vars_labels)

        push!(vars_indices, bus_var_Idx)
    end

    return vars_indices
end


function get_a_node_ur_ui_vars_syms_in_system(plant; network_vars_labels = network_vars_labels )

    return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = [:u_r, :u_i] )

end



function get_a_node_ur_ui_indices_in_system(plant; network_vars_labels = network_vars_labels )

    return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = [:u_r, :u_i] )

end



function get_a_gen_states_vars_syms_in_system(plant; network_vars_labels = network_vars_labels, no_control_device = false )

    if  no_control_device == false

        return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.syms )
    else
        return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.Gens.syms )
    end
    
end


function get_a_gen_states_indices_in_system(plant; network_vars_labels = network_vars_labels, no_control_device = false )

    if no_control_device == false

        return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.syms )
    else
        return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.Gen.syms )
    end
            
end


function get_a_gen_pure_states_vars_syms_in_system(plant; network_vars_labels = network_vars_labels, no_control_device = false  )

    if no_control_device == false

        return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.state_vars_syms )
    else
        return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.Gen.state_vars_syms )
    end            

end


function get_a_gen_pure_states_indices_in_system(plant; network_vars_labels = network_vars_labels, no_control_device = false )

    if no_control_device == false
        
        return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.state_vars_syms )
        
    else
        
        return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.Gen.state_vars_syms )
        
    end            
end


function get_a_gen_im_vars_syms_in_system(plant; network_vars_labels = network_vars_labels, no_control_device = false  )

    if no_control_device == false

        return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.im_vars_syms )
    else
        return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.Gen.state_vars_syms )
    end            

end


function get_a_gen_im_vars_indices_in_system(plant; network_vars_labels = network_vars_labels, no_control_device = false )

    if no_control_device == false
        
        return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.im_vars_syms )
        
    else
        
        return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.Gen.state_vars_syms )
        
    end
    
        

end


function get_a_gen_im_algebraic_vars_syms_in_system(plant; network_vars_labels = network_vars_labels )

    return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.im_algebraic_vars_syms )
end

# im_algebraic_vars_syms

function get_a_gen_im_algebraic_vars_indices_in_system(plant; network_vars_labels = network_vars_labels )

    return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.im_algebraic_vars_syms )
        
end



function get_a_gen_stab_states_vars_syms_in_system(plant; network_vars_labels = network_vars_labels, no_control_device = false )

    if no_control_device == false

        return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.stab_state_vars_syms )
    else
        return get_a_node_states_vars_syms_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.Gen.stab_state_vars_syms )

    end
    

end


function get_a_gen_stab_states_indices_in_system(plant; network_vars_labels = network_vars_labels, no_control_device = false  )

    if no_control_device == false

        return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.stab_state_vars_syms )
    else
        return get_a_node_states_indices_in_system(; network_vars_labels = network_vars_labels, bus_name = plant.name, vars = plant.Gen.stab_state_vars_syms )
    end
    

end

# get_a_node_ur_ui_indices_in_system(plant; network_vars_labels = network_vars_labels )

#--------------------------------------------------------
#--------------------------------------------------------


function get_nodes_ur_ui_vars_syms_in_system(netd  )

    network_vars_labels = get_network_vars_labels(netd )
    
    return [ get_a_node_ur_ui_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes))  ]    
        
end


function get_nodes_ur_ui_indices_in_system(netd )

    network_vars_labels = get_network_vars_labels(netd )
    
    return [ get_a_node_ur_ui_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes))  ]
   
end



function get_gens_states_vars_syms_in_system(netd; no_control_device = false )

    if no_control_device == false

        network_vars_labels = get_network_vars_labels(netd )    
        return [ get_a_gen_states_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = get_network_vars_labels(netd )
        [ get_a_gen_states_vars_syms_in_system( a_node.Gen; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator ]
    end
            
end


function get_gens_states_indices_in_system(netd; no_control_device = false  )

    if no_control_device == false
        network_vars_labels = get_network_vars_labels(netd )

        return [ get_a_gen_states_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = get_network_vars_labels(netd )

        return [ get_a_gen_states_indices_in_system( a_node.Gen; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    end
    
end


function get_gens_im_vars_syms_in_system(netd; no_control_device = false )

    if no_control_device == false

        network_vars_labels = get_network_vars_labels(netd )    
        return [ get_a_gen_im_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = get_network_vars_labels(netd )
        [ get_a_gen_states_vars_syms_in_system( a_node.Gen; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator ]
    end
            
end

# im_algebraic_vars

function get_gens_im_algebraic_vars_indices_in_system(netd  )

    network_vars_labels = get_network_vars_labels(netd )

    return [ get_a_gen_im_algebraic_vars_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    
end


function get_gens_im_algebraic_vars_syms_in_system(netd )

    network_vars_labels = get_network_vars_labels(netd )
    
    return [ get_a_gen_im_algebraic_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    
end



function get_gens_im_vars_indices_in_system(netd; no_control_device = false  )

    if no_control_device == false
        network_vars_labels = get_network_vars_labels(netd )

        return [ get_a_gen_im_vars_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = get_network_vars_labels(netd )

        return [ get_a_gen_states_indices_in_system( a_node.Gen; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    end
    
end

#---------------------------------------------------
 
function get_gens_pure_states_indices_in_system(
    netd; no_control_device = false )

    if no_control_device == false

        network_vars_labels = get_network_vars_labels(netd )
    
        return [ get_a_gen_pure_states_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = get_network_vars_labels(netd )
        
        return [ get_a_gen_pure_states_indices_in_system( a_node.Gen; network_vars_labels = network_vars_labels  ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    end
    
end


function get_gens_pure_states_vars_syms_in_system(netd; no_control_device = false )

    if no_control_device == false

        network_vars_labels = get_network_vars_labels(netd )
    
        return [ get_a_gen_pure_states_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else

        network_vars_labels = get_network_vars_labels(netd )

        return [ get_a_gen_pure_states_vars_syms_in_system( a_node.Gen; network_vars_labels = network_vars_labels  ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    end
                    
end

function get_gens_stab_states_indices_in_system(netd; no_control_device = false )

    if no_control_device == false
        network_vars_labels = get_network_vars_labels(netd )

        return [ get_a_gen_stab_states_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = get_network_vars_labels(netd )

        return [ get_a_gen_stab_states_indices_in_system( a_node.Gen; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]        
    end
    
end

# ---------------------------------------------------
#----------------------------------------------------


function get_industrial_nodes_ur_ui_vars_syms_in_state(netd; no_control_device = false )

    if no_control_device == false
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes )

        return [ get_a_node_ur_ui_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes))  ]
    else
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes, no_control_device = true  )

        return [ get_a_node_ur_ui_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes))  ]

    end
    
end

function get_industrial_nodes_ur_ui_indices_in_state(netd; no_control_device = false )

    if no_control_device == false
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes )

        return [ get_a_node_ur_ui_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes))  ]
    else
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes, no_control_device = true  )

        return [ get_a_node_ur_ui_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes))  ]
    end
    
end



function get_industrial_gens_pure_states_vars_syms_in_state(netd; no_control_device = false )

    if no_control_device == false
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes )

        return [ get_a_gen_pure_states_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes, no_control_device = true  )

        return [ get_a_gen_pure_states_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels, no_control_device = true  ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    end
    
end


function get_industrial_gens_pure_states_indices_in_state(netd; no_control_device = false )

    if no_control_device == false
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes )

        return [ get_a_gen_pure_states_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes, no_control_device = true  )

        return [ get_a_gen_pure_states_indices_in_system( a_node; network_vars_labels = network_vars_labels, no_control_device = true  ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]        
    end
    
end


function get_industrial_gens_stab_states_vars_syms_in_state(netd; no_control_device = false )

    if no_control_device == false
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes )

        return [ get_a_gen_stab_states_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes, no_control_device = true  )

        return [ get_a_gen_stab_states_vars_syms_in_system( a_node; network_vars_labels = network_vars_labels, no_control_device = true  ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]        
    end
    
end


function get_industrial_gens_stab_states_indices_in_state(netd; no_control_device = false )

    if no_control_device == false
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes )

        return [ get_a_gen_stab_states_indices_in_system( a_node; network_vars_labels = network_vars_labels ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    else
        network_vars_labels = generate_industrial_model_sym(; nodes = netd.nodes, no_control_device = true )

        return [ get_a_gen_stab_states_indices_in_system( a_node; network_vars_labels = network_vars_labels, no_control_device = true  ) for a_node in  collect(values(netd.nodes)) if a_node.Bus_type == :Generator  ]
    end
    
end



function get_gens_im_pure_states_vars_syms_in_state(netd )

    network_vars_labels = generate_industrial_model_sym(
        ; nodes = netd.nodes )

    return [ get_a_gen_pure_states_vars_syms_in_system(
        a_node; network_vars_labels = network_vars_labels )
             for a_node in  collect(values(netd.nodes))
                 if a_node.Bus_type == :Generator  ]    
end


function get_gens_im_pure_states_indices_in_state(netd )

    network_vars_labels = generate_industrial_model_sym(
        ; nodes = netd.nodes )

    return [ get_a_gen_pure_states_indices_in_system(
        a_node; network_vars_labels = network_vars_labels )
             for a_node in  collect(values(netd.nodes))
                 if a_node.Bus_type == :Generator  ]    
end

#----------------------------------------------------
#----------------------------------------------------




function get_a_node_δ_ed_dash_eq_dash_indices(; network_vars_labels = network_vars_labels, bus_name = bus_name )

    vars_indices = [ ]

    for a_var in [:δ, :ed_dash, :eq_dash]

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   = findfirst(comp_name -> comp_name == bus_var_label,  network_vars_labels)

        push!(vars_indices, bus_var_Idx)
    end

    return vars_indices
end



function get_a_node_δ_ω_ed_dash_eq_dash_indices(; network_vars_labels = network_vars_labels, bus_name = bus_name )

    vars_indices = [ ]

    for a_var in [:δ, :ω, :ed_dash, :eq_dash]

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   = findfirst(comp_name -> comp_name == bus_var_label,  network_vars_labels)

        push!(vars_indices, bus_var_Idx)
    end

    return vars_indices
end


#--------------------------------------------

function get_some_state_vars_idxs_in_a_plant(plant; some_state_vars = [:δ, :ω, :ed_dash, :eq_dash] ,  no_control_device = false )

    if  no_control_device == false
        dict_state_vars_syms = plant.dict_state_vars_syms

        vars_indices = [  ]

        for a_var in some_state_vars
            push!(vars_indices, dict_state_vars_syms[a_var])

        end

        return vars_indices
    else
        dict_state_vars_syms = plant.Gen.dict_state_vars_syms

        vars_indices = [  ]

        for a_var in some_state_vars
            push!(vars_indices, dict_state_vars_syms[a_var])

        end

        return vars_indices
    end
    
end


function get_gens_nodes_some_state_vars_idx_in_plants( nodes_collection; some_state_vars = [:δ, :ω, :ed_dash, :eq_dash],  no_control_device = false )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_some_state_vars_idxs_in_a_plant(comp; some_state_vars = some_state_vars, no_control_device = no_control_device ) for comp in comp_collection_values if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_some_state_vars_idxs_in_a_plant(comp; some_state_vars = some_state_vars, no_control_device = no_control_device ) for comp in nodes_collection if comp.Bus_type == :Generator ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [get_some_state_vars_idxs_in_a_plant(comp; some_state_vars = some_state_vars, no_control_device = no_control_device ) for comp in comp_collection_values if comp.Bus_type == :Generator ]        
    end

end



function get_non_gens_nodes_some_state_vars_idx_in_plants( nodes_collection; some_state_vars = [:P, :Q ],  no_control_device = false )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_some_state_vars_idxs_in_a_plant(comp; some_state_vars = some_state_vars, no_control_device = no_control_device ) for comp in comp_collection_values if comp.Bus_type != :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_some_state_vars_idxs_in_a_plant(comp; some_state_vars = some_state_vars, no_control_device = no_control_device ) for comp in nodes_collection if comp.Bus_type != :Generator ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [get_some_state_vars_idxs_in_a_plant(comp; some_state_vars = some_state_vars, no_control_device = no_control_device ) for comp in comp_collection_values if comp.Bus_type != :Generator ]        
    end

end


function get_nodes_some_state_vars_idx_in_plants( nodes_collection; some_state_vars = [:P, :Q ],  no_control_device = false )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_some_state_vars_idxs_in_a_plant(comp; some_state_vars = some_state_vars, no_control_device = no_control_device ) for comp in comp_collection_values  ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_some_state_vars_idxs_in_a_plant(comp; some_state_vars = some_state_vars, no_control_device = no_control_device ) for comp in nodes_collection  ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [get_some_state_vars_idxs_in_a_plant(comp; some_state_vars = some_state_vars, no_control_device = no_control_device ) for comp in comp_collection_values  ]        
    end

end


#--------------------------------------------


function get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = bus_name, some_states_syms = [:δ, :ed_dash, :eq_dash]  )

    vars_indices = [ ]

    for a_var in some_states_syms

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   = findfirst(comp_name -> comp_name == bus_var_label,  network_vars_labels)

        push!(vars_indices, bus_var_Idx)
    end

    return vars_indices
end


function get_gens_nodes_some_state_vars_Idxs_in_net_state( nodes_collection; network_vars_labels = network_vars_labels, some_state_vars = [:δ, :ω, :ed_dash, :eq_dash],  no_control_device = false )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = comp.name, some_states_syms = some_state_vars  )  for comp in comp_collection_values if comp.Bus_type == :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = comp.name, some_states_syms = some_state_vars  )  for comp in nodes_collection if comp.Bus_type == :Generator ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = comp.name, some_states_syms = some_state_vars  )   for comp in comp_collection_values if comp.Bus_type == :Generator ]        
    end

end


function get_non_gens_nodes_some_state_vars_Idxs_in_net_state( nodes_collection; network_vars_labels = network_vars_labels, some_state_vars = [:P, :Q],  no_control_device = false )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = comp.name, some_states_syms = some_state_vars  )  for comp in comp_collection_values if comp.Bus_type != :Generator ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = comp.name, some_states_syms = some_state_vars  )  for comp in nodes_collection if comp.Bus_type != :Generator ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = comp.name, some_states_syms = some_state_vars  )   for comp in comp_collection_values if comp.Bus_type != :Generator ]        
    end

end


function get_nodes_some_state_vars_Idxs_in_net_state( nodes_collection; network_vars_labels = network_vars_labels, some_state_vars = [:P, :Q],  no_control_device = false )
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))

        return [ get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = comp.name, some_states_syms = some_state_vars  )  for comp in comp_collection_values ]
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        return [ get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = comp.name, some_states_syms = some_state_vars  )  for comp in nodes_collection ]
    else
        comp_collection_values =
            collect(values(nodes_collection))
        
        return [ get_some_states_Idxs_of_a_node_in_net_state(; network_vars_labels = network_vars_labels, bus_name = comp.name, some_states_syms = some_state_vars  )   for comp in comp_collection_values ]        
    end

end

#--------------------------------------------
#--------------------------------------------

function get_models_networks_labels_and_gens_nodes_some_vars_Idxs(
    netd;
    some_state_vars =
        [:δ, :ω, :ed_dash, :eq_dash] )

#--------------------------------------------    
        
    network_vars_labels_hybrid =
        get_network_vars_labels( netd )

    network_vars_labels_industrial =
        generate_industrial_model_sym(
            ; nodes = netd.nodes )
    
    network_vars_labels_im =
        generate_im_sym(
            ; nodes = netd.nodes )
    
    #----------------------------------------------------

    Idxs_hybrid =
        get_gens_nodes_some_state_vars_Idxs_in_net_state(
            netd.nodes;
            network_vars_labels =
                network_vars_labels_hybrid,
            some_state_vars =
                some_state_vars )

    Idxs_industrial =
        get_gens_nodes_some_state_vars_Idxs_in_net_state(
            netd.nodes;
            network_vars_labels =
                network_vars_labels_industrial,
            some_state_vars = some_state_vars )

    Idxs_im =
        get_gens_nodes_some_state_vars_Idxs_in_net_state(
            netd.nodes;
            network_vars_labels =
                network_vars_labels_im,
            some_state_vars = some_state_vars )

    return (; network_vars_labels_hybrid,
            network_vars_labels_industrial,
            network_vars_labels_im,
            Idxs_hybrid,
            Idxs_industrial,
            Idxs_im  )
    
end

function get_models_networks_labels(netd )

#--------------------------------------------    
        
    network_vars_labels_hybrid =
        get_network_vars_labels( netd )

    network_vars_labels_industrial =
        generate_industrial_model_sym(
            ; nodes = netd.nodes )
    
    network_vars_labels_im =
        generate_im_sym(
            ; nodes = netd.nodes )
    
    #----------------------------------------------------

    return (; network_vars_labels_hybrid,
            network_vars_labels_industrial,
            network_vars_labels_im  )
    
end


function get_models_gens_nodes_some_vars_Idxs(
    netd;
    some_state_vars =
        [:δ, :ω, :ed_dash, :eq_dash] )

#--------------------------------------------    
        
    network_vars_labels_hybrid =
        get_network_vars_labels( netd )

    network_vars_labels_industrial =
        generate_industrial_model_sym(
            ; nodes = netd.nodes )
    
    network_vars_labels_im =
        generate_im_sym(
            ; nodes = netd.nodes )
    
    #----------------------------------------------------

    Idxs_hybrid =
        get_gens_nodes_some_state_vars_Idxs_in_net_state(
            netd.nodes;
            network_vars_labels =
                network_vars_labels_hybrid,
            some_state_vars =
                some_state_vars )

    Idxs_industrial =
        get_gens_nodes_some_state_vars_Idxs_in_net_state(
            netd.nodes;
            network_vars_labels =
                network_vars_labels_industrial,
            some_state_vars =
                some_state_vars )

    Idxs_im =
        get_gens_nodes_some_state_vars_Idxs_in_net_state(
            netd.nodes;
            network_vars_labels =
                network_vars_labels_im,
            some_state_vars =
                some_state_vars )

    return (; Idxs_hybrid,
            Idxs_industrial,
            Idxs_im  )
    
end



function get_models_nodes_some_vars_Idxs(
    netd;
    some_state_vars =
        [ :u_r, :u_i ] )

#--------------------------------------------    
        
    network_vars_labels_hybrid =
        get_network_vars_labels( netd )

    network_vars_labels_industrial =
        generate_industrial_model_sym(
            ; nodes = netd.nodes )
    
    network_vars_labels_im =
        generate_im_sym(
            ; nodes = netd.nodes )
    
    #----------------------------------------------------

    Idxs_hybrid =
        get_nodes_some_state_vars_Idxs_in_net_state(
            netd.nodes;
            network_vars_labels =
                network_vars_labels_hybrid,
            some_state_vars =
                some_state_vars )

    Idxs_industrial =
        get_nodes_some_state_vars_Idxs_in_net_state(
            netd.nodes;
            network_vars_labels =
                network_vars_labels_industrial,
            some_state_vars = some_state_vars )

    Idxs_im =
        get_nodes_some_state_vars_Idxs_in_net_state(
            netd.nodes;
            network_vars_labels =
                network_vars_labels_im,
            some_state_vars =
                some_state_vars )

    return (; Idxs_hybrid,
            Idxs_industrial,
            Idxs_im  )
    
end


#--------------------------------------------
#--------------------------------------------


"""
This function returns indices of state_algebraic variables of specified nodes in the systems.

It is meant to extract the indices in sol.
"""
function get_nodes_state_algb_vars_indices_in_system(
    ; network_vars_labels =
        network_vars_labels,
    nodes_name = ["bus1", "bus2"],
    vars = [:ω, :δ])

    vars_indices = []

    for bus_name in nodes_name
        tup_vars_indices =
            get_a_node_state_algb_vars_indices_in_system(
                ; network_vars_labels =
                    network_vars_labels,
                bus_name =
                    bus_name, vars = vars )
        
        push!(vars_indices, tup_vars_indices)
    end

    vars_indices = [an_idx
                    for an_idx in vars_indices
                        if an_idx[1] != nothing]

    return Tuple(vars_indices)
    
    # if filter_nothing == true
        
    #     vars_indices = [an_idx for an_idx in vars_indices if an_idx[1] != nothing]
        
    #     return Tuple(vars_indices)
    # else
    #     return vars_indices
    # end
            
end


function get_all_nodes_state_algb_vars_indices_in_system(
    nodes;
    network_vars_labels =
        network_vars_labels,
    vars = [:ω, :δ])

    nodes_name = collect(keys(nodes))
    
    vars_indices = []

    for bus_name in nodes_name
        tup_vars_indices =
            get_a_node_state_algb_vars_indices_in_system(
                ; network_vars_labels =
                    network_vars_labels,
                bus_name =
                    bus_name,
                vars = vars )
        
        push!(vars_indices, tup_vars_indices)
        
    end

    vars_indices = [an_idx
                    for an_idx in vars_indices
                        if an_idx[1] != nothing]
        
    return Tuple(vars_indices)
end


#########################################################
# Generic functions
#########################################################


function create_comps_non_flattened_vectors_Idx(comps_values_vectors)

    vectors_dims = ones(Int, length(comps_values_vectors))
    vectors_offset = create_offsets(vectors_dims)
    vectors_Idx = create_idxs(vectors_offset, vectors_dims)

    # non_flattened_vectors_Idx
    
    return vectors_Idx
end


function create_comps_flattened_vectors_Idx(comps_values_vectors)

    vectors_dims    = length.(comps_values_vectors)
    vectors_offset  = create_offsets(vectors_dims)
    vectors_Idx     = create_idxs(vectors_offset, vectors_dims)

    # flattened_vectors_Idx
    
    return vectors_Idx
end


function create_comps_non_flattened_vectors_and_Idx(comps_values_vectors)

    comps_values_vectors_Idx = create_comps_non_flattened_vectors_Idx(comps_values_vectors)

    # non_flattened_vectors, non_flattened_vectors_Idx
    
    return comps_values_vectors, comps_values_vectors_Idx
end



function create_comps_flattened_vectors_and_Idx(comps_values_vectors)

    comps_values_vectors_Idx = create_comps_flattened_vectors_Idx(comps_values_vectors)

    # flattened_vectors, flattened_vectors_Idx
    
    return [comps_values_vectors...;], comps_values_vectors_Idx
end


function create_non_flattened_dims_offset_Idx(comps_vector_type::Vector)

    type_dims = ones(Int, length(comps_vector_type))
    type_offset = create_offsets(type_dims)
    type_Idx = create_idxs(type_offset, type_dims)

    # type_dims, type_offset, type_Idx
    
    return type_dims, type_offset, type_Idx
    
end


function create_flattened_dims_offset_Idx(comps_vector_type::Vector)

    type_dims = length.(comps_vector_type)
    type_offset = create_offsets(type_dims)
    type_Idx = create_idxs(type_offset, type_dims)

    # type_dims, type_offset, type_Idx
    
    return type_dims, type_offset, type_Idx
    
end


function create_comps_non_flattened_type_vectors_Idx(type_func, comp_collection)

    comps_vector_type = type_func(comp_collection)
 
    # non_flattened_vectors_Idx
    
    return create_comps_non_flattened_vectors_Idx(comps_vector_type)
end


function create_comps_flattened_type_vectors_Idx(type_func, comp_collection)

    comps_vector_type = type_func(comp_collection)

    # flattened_vectors_Idx
    
    return create_comps_flattened_vectors_Idx( comps_vector_type )
end


function create_comps_non_flattened_type_vectors_and_Idx(type_func, comp_collection::OrderedDict)

    comps_vector_type = type_func(comp_collection)

    # vector_type, type_Idx
    
    return create_comps_non_flattened_vectors_and_Idx(comps_vector_type)
    
end


function create_comps_flattened_type_vectors_and_Idx(type_func, comp_collection::OrderedDict)

    comps_vector_type = type_func(comp_collection)

    # vector_type, type_Idx
    
    return create_comps_flattened_vectors_and_Idx(comps_vector_type)
    
end


function create_non_flattened_type_dims_offset_Idx(type_func, comp_collection::OrderedDict)

    comps_vector_type = type_func(comp_collection)

    # type_dims, type_offset, type_Idx
    
    return create_non_flattened_dims_offset_Idx(comps_vector_type)
    
end



function  create_flattened_type_dims_offset_Idx(type_func, comp_collection::OrderedDict)

    comps_vector_type = type_func(comp_collection)

    # type_dims, type_offset, type_Idx
    
    return  create_flattened_dims_offset_Idx(comps_vector_type)
    
end


function  create_flattened_type_dims_offset_Idx(type_func, comp_collection)

    comps_vector_type = type_func(comp_collection)

    # type_dims, type_offset, type_Idx
    
    return  create_flattened_dims_offset_Idx(comps_vector_type)
    
end


function create_comps_vector_view(vectors, vectors_Idx)

    return [view(vectors, vectors_Idx[Ind]) for Ind in collect(1:length(vectors_Idx))]
    
end


function create_comps_non_flattened_vector_view(comps_values_vectors, comps_values_vectors_Idx)

    return [view(comps_values_vectors, comps_values_vectors_Idx[Ind]) for Ind in collect(1:length(comps_values_vectors_Idx))]
    
end


function create_comps_flattened_vector_view(flattened_vec, comps_values_vectors_Idx)

    return [view(flattened_vec, comps_values_vectors_Idx[Ind]) for Ind in collect(1:length(comps_values_vectors_Idx))]
    
end

########################################################


function create_nodes_and_edges_non_flattened_type_size_offset_dims_Idx(type_func, nodes::OrderedDict, edges::OrderedDict; state=false)

    nodes_vector_type = type_func(nodes)

    edges_vector_type = type_func(edges)

    if state == true
        
        nodes_dims   = length.(nodes_vector_type)
        edges_dims   = length.(edges_vector_type)
    else
        nodes_dims   = ones(Int, length(nodes_vector_type))
        edges_dims   = ones(Int, length(edges_vector_type))
    end

    nodes_dim_size, nodes_offset, nodes_Idx =  create_size_offset_Idx(nodes_dims)
    
    edges_dim_size, edges_offset, edges_Idx =  create_size_offset_Idx(edges_dims; counter = nodes_dim_size)
    

    return nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_Idx, edges_Idx
end


function create_nodes_and_edges_non_flattened_type_size_offset_dims_Idx(type_func, nd::NetworkData; state=false)
        
    # nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_Idx, edges_Idx
    
    return create_nodes_and_edges_non_flattened_type_size_offset_dims_Idx(type_func, nd.nodes, nd.edges; state = state)
end



function create_nodes_and_edges_flattened_type_size_offset_dims_Idx(type_func, nodes::OrderedDict, edges::OrderedDict; state=false)

    nodes_vector_type = type_func(nodes)

    edges_vector_type = type_func(edges)

    if state == true
        nodes_dims   = nodes_vector_type
        edges_dims   = edges_vector_type
    else
        nodes_dims   = length.(nodes_vector_type)
        edges_dims   = length.(edges_vector_type)
    end

    nodes_dim_size, nodes_offset, nodes_Idx =  create_size_offset_Idx(nodes_dims)
    
    edges_dim_size, edges_offset, edges_Idx =  create_size_offset_Idx(edges_dims; counter = nodes_dim_size)
    

    return nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_Idx, edges_Idx
end


function create_nodes_and_edges_flattened_type_size_offset_dims_Idx(
    type_func,  nd::NetworkData; state=false)

    # nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_Idx, edges_Idx
    
    return create_nodes_and_edges_flattened_type_size_offset_dims_Idx(type_func, nd.nodes, nd.edges; state = state )
end


########################################################
########################################################

########################################################
########################################################


function create_nodes_and_edges_states_Idx(nd::NetworkData)

    # nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_Idx, edges_Idx

    return create_nodes_and_edges_flattened_type_size_offset_dims_Idx(get_components_state_algebraic_vars_dim, nd.nodes, nd.edges; state = true )
        
end

function create_nodes_and_edges_states_Idx( nodes::OrderedDict, edges::OrderedDict)
    # nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_Idx, edges_Idx

    return create_nodes_and_edges_flattened_type_size_offset_dims_Idx(get_components_state_algebraic_vars_dim, nodes, edges; state = true )
    
end

########################################################


function create_nodes_and_edges_dyn_func_Idx(get_components_dyn_func, nd::NetworkData)

     # nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_Idx, edges_Idx
    
    return  create_nodes_and_edges_non_flattened_type_size_offset_dims_Idx(get_components_dyn_func, nd)
    
end


function create_nodes_and_edges_dyn_func_Idx(get_components_dyn_func, nodes::OrderedDict, edges::OrderedDict )

    # (; nodes_fun_size, edges_fun_size, nodes_fun_offset, edges_fun_offset, nodes_fun_dims, edges_fun_dims, nodes_fun_Idx, edges_fun_Idx)
    
    return  create_nodes_and_edges_non_flattened_type_size_offset_dims_Idx(get_components_dyn_func, nodes, edges )
    
end

########################################################
########################################################


"""
get_components_params_value_vectors

get_components_control_sig

get_components_output_sig

get_components_cb_dyn_param_state_sw


"""

########################################################


function create_params_views_and_Idx(nd::NetworkData, nodes_fun_params, edges_fun_params; nodes_fun_idx = 2, edges_fun_idx = 1 )

    nodes_cb_sw, nodes_f_t, nodes_param = nodes_fun_params
    
    edges_param = edges_fun_params

    #-----------------------------------------

    nodes_dyn_func = get_components_dyn_func_by_idx(nd.nodes; idx = nodes_fun_idx )

    edges_dyn_func = get_components_dyn_func_by_idx(nd.edges;idx = edges_fun_idx )
    
    #-----------------------------------------
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(nd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(nd.edges)
    nodes_incident_edges = get_nodes_incident_edges(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)

    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset = create_offsets(nodes_dims)

    nodes_u_Idx = get_components_ur_ui_Idx_in_state(nd.nodes)
    
    # nodes_u_Idx = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(nd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx = create_comps_non_flattened_vectors_Idx(nodes_param)

    #---------------------------------------------------
      

    edges_param_Idx = create_comps_non_flattened_vectors_Idx(edges_param)

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

    #-------------------------------------------------
    #--------------------------------------------------
    
    sd_nodes_param   = ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  )
    
    sd_edges_param  = ( edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx )
        
    return  (sd_nodes_param, sd_edges_param)

                                     
end 


########################################################


function external_pf_create_params_views_and_Idx(nd::NetworkData, nodes_fun_params, edges_fun_params )

    nodes_cb_sw, nodes_f_t, nodes_param = nodes_fun_params
    
    edges_param = edges_fun_params

    #-----------------------------------------

    nodes_dyn_func = get_components_dyn_func_by_idx(nd.nodes; idx = 3)

    edges_dyn_func = get_components_dyn_func_by_idx(nd.edges;idx = 1 )
    
    #-----------------------------------------
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(nd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(nd.edges)
    nodes_incident_edges = get_nodes_incident_edges(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)

    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset = create_offsets(nodes_dims)

    nodes_u_Idx = get_components_ur_ui_Idx_in_state(nd.nodes) 
    # nodes_u_Idx = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(nd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx = create_comps_non_flattened_vectors_Idx(nodes_param)

    #---------------------------------------------------
      

    edges_param_Idx = create_comps_non_flattened_vectors_Idx(edges_param)

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

    #-------------------------------------------------
    #--------------------------------------------------
    
    sd_nodes_param   = ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  )
    
    sd_edges_param  = ( edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx )
        
    return  (sd_nodes_param, sd_edges_param)

                                     
end 

########################################################


function hybrid_pf_create_params_views_and_Idx(nd::NetworkData, nodes_fun_params, edges_fun_params )

    nodes_cb_sw, nodes_f_t, nodes_param = nodes_fun_params
    
    edges_param = edges_fun_params

    #-----------------------------------------

    nodes_dyn_func = get_components_dyn_func_by_idx(nd.nodes; idx = 5)

    edges_dyn_func = get_components_dyn_func_by_idx(nd.edges; idx = 1 )
    
    #-----------------------------------------
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(nd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(nd.edges)
    nodes_incident_edges = get_nodes_incident_edges(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)

    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset = create_offsets(nodes_dims)

    nodes_u_Idx = get_components_ur_ui_Idx_in_state(nd.nodes) 
    # nodes_u_Idx = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(nd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx = create_comps_non_flattened_vectors_Idx(nodes_param)

    #---------------------------------------------------
      

    edges_param_Idx = create_comps_non_flattened_vectors_Idx(edges_param)

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

    #-------------------------------------------------
    #--------------------------------------------------
    
    sd_nodes_param   = ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  )
    
    sd_edges_param  = ( edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx )
        
    return  (sd_nodes_param, sd_edges_param)

                                     
end 


########################################################


function node_pf_create_params_views_and_Idx(nd::NetworkData, nodes_fun_params, edges_fun_params )

    nodes_cb_sw, nodes_f_t, nodes_param = nodes_fun_params
    
    edges_param = edges_fun_params

    #-----------------------------------------

    nodes_dyn_func = get_components_dyn_func_by_idx(nd.nodes; idx = 4)

    edges_dyn_func = get_components_dyn_func_by_idx(nd.edges; idx = 1 )
    
    #-----------------------------------------
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(nd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(nd.edges)
    nodes_incident_edges = get_nodes_incident_edges(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)

    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset = create_offsets(nodes_dims)

    nodes_u_Idx = get_components_ur_ui_Idx_in_state(nd.nodes) 
    # nodes_u_Idx = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(nd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx = create_comps_non_flattened_vectors_Idx(nodes_param)

    #---------------------------------------------------
      

    edges_param_Idx = create_comps_non_flattened_vectors_Idx(edges_param)

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

    #-------------------------------------------------
    #--------------------------------------------------
    
    sd_nodes_param   = ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  )
    
    sd_edges_param  = ( edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx )
        
    return  (sd_nodes_param, sd_edges_param)

                                     
end 


########################################################

function global_pf_create_params_views_and_Idx(nd::NetworkData, nodes_fun_params, edges_fun_params )

    nodes_cb_sw, nodes_f_t, nodes_param = nodes_fun_params
    
    edges_param = edges_fun_params

    #-----------------------------------------

    nodes_dyn_func = get_components_dyn_func_by_idx(nd.nodes; idx = 3)

    edges_dyn_func = get_components_dyn_func_by_idx(nd.edges;idx = 1 )
    
    #-----------------------------------------
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(nd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(nd.edges)
    nodes_incident_edges = get_nodes_incident_edges(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)

    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset = create_offsets(nodes_dims)

    nodes_u_Idx = get_components_ur_ui_Idx_in_state(nd.nodes) 
    # nodes_u_Idx = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(nd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx = create_comps_non_flattened_vectors_Idx(nodes_param)

    #---------------------------------------------------
      

    edges_param_Idx = create_comps_non_flattened_vectors_Idx(edges_param)

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

    #-------------------------------------------------
    #--------------------------------------------------
    
    sd_nodes_param   = ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  )
    
    sd_edges_param  = ( edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx )
        
    return  (sd_nodes_param, sd_edges_param)

                                     
end 

########################################################


function create_params_views_and_Idx_external(nd::NetworkData, nodes_fun_params, edges_fun_params; nodes_fun_idx = 2, edges_fun_idx = 1 )

    nodes_cb_sw, nodes_f_t, nodes_param = nodes_fun_params
    
    edges_param = edges_fun_params

    #-----------------------------------------

    nodes_dyn_func = get_components_dyn_func_by_idx(nd.nodes; idx = nodes_fun_idx )

    edges_dyn_func = get_components_dyn_func_by_idx(nd.edges;idx = edges_fun_idx )
    
    #-----------------------------------------
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(nd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(nd.edges)
    nodes_incident_edges = get_nodes_incident_edges(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)

    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset = create_offsets(nodes_dims)

    nodes_u_Idx = get_components_ur_ui_Idx_in_state(nd.nodes) 
    # nodes_u_Idx = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(nd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx = create_comps_non_flattened_vectors_Idx(nodes_param)

    #---------------------------------------------------
      

    edges_param_Idx = create_comps_non_flattened_vectors_Idx(edges_param)

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

    #-------------------------------------------------
    #--------------------------------------------------
    
    sd_nodes_param   = ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  )
    
    sd_edges_param  = ( edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx )
        
    return  (sd_nodes_param, sd_edges_param)

                                     
end 


########################################################

function create_nd_non_flattened_params_views_and_Idx_external_f_t(nd::NetworkData, nodes_fun_params, edges_fun_params )

    nodes_cb_sw, nodes_f_t, nodes_param = nodes_fun_params
    
    edges_param = edges_fun_params

    #-----------------------------------------

    nodes_dyn_func = get_components_dyn_func_non_flat(nd.nodes)

    edges_dyn_func = get_components_dyn_func(nd.edges)
    
    #-----------------------------------------    
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(nd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(nd.edges)
    nodes_incident_edges = get_nodes_incident_edges(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)


    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset        = create_offsets(nodes_dims)

    nodes_u_Idx          = get_components_ur_ui_Idx_in_state(nd.nodes) 
    # nodes_u_Idx         = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(nd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx = create_comps_non_flattened_vectors_Idx(nodes_param)

    #---------------------------------------------------
      

    edges_param_Idx = create_comps_non_flattened_vectors_Idx(edges_param)

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

    #-------------------------------------------------
    #--------------------------------------------------
    
    sd_nodes_param   = ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  )
    
    sd_edges_param  = ( edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx )
    
    
    return  (sd_nodes_param, sd_edges_param)

                                     
end 


function create_nd_non_flattened_params_views_and_Idx_with_pf(nd::NetworkData, nodes_fun_params, edges_fun_params )

    nodes_cb_sw, nodes_f_t, nodes_param = nodes_fun_params
    
    edges_param = edges_fun_params

    #-----------------------------------------

    nodes_dyn_func = get_components_dyn_func_non_flat(nd.nodes)

    edges_dyn_func = get_components_dyn_func(nd.edges)
    
    #-----------------------------------------

    ur_ui_Idx_in_state = get_components_ur_ui_Idx_in_state(nd.nodes)

    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(nd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(nd.edges)
    nodes_incident_edges = get_nodes_incident_edges(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)


    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset        = create_offsets(nodes_dims)
    
    nodes_u_Idx          = get_components_ur_ui_Idx_in_state(nd.nodes)
    # nodes_u_Idx         = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(nd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx = create_comps_non_flattened_vectors_Idx(nodes_param)

    #---------------------------------------------------
      

    edges_param_Idx = create_comps_non_flattened_vectors_Idx(edges_param)

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

    #-------------------------------------------------
    #--------------------------------------------------
    
    sd_nodes_param   = ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx, ur_ui_Idx_in_state )
    
    sd_edges_param  = ( edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx )
    
    
    return  (sd_nodes_param, sd_edges_param)

                                     
end 


function create_nd_non_flattened_params_views_and_Idx(nd::NetworkData, nodes_fun_params, edges_fun_params )

    nodes_cb_sw, nodes_f_t, nodes_param = nodes_fun_params
    
    edges_param = edges_fun_params

    #-----------------------------------------

    nodes_dyn_func = get_components_dyn_func_non_flat(nd.nodes)

    edges_dyn_func = get_components_dyn_func(nd.edges)
    
    #-----------------------------------------
    #-----------------------------------------
    
    nodes_src_edges      = get_nodes_src_edges(nd.edges)
    nodes_dst_edges      = get_nodes_dst_edges(nd.edges)
    nodes_incident_edges = get_nodes_incident_edges(nd.edges)

    nodes_dim_size, edges_dim_size, nodes_offset, edges_offset,  nodes_dims, edges_dims, nodes_state_Idx, edges_state_Idx = create_nodes_and_edges_states_Idx(nd.nodes, nd.edges)


    # These are needed for voltages at src and dst of an edge in
    # edges_src_volt_view, and edges_dst_volt_view
    
    # nodes_offset        = create_offsets(nodes_dims)
    
    nodes_u_Idx          = get_components_ur_ui_Idx_in_state(nd.nodes)
    # nodes_u_Idx         = create_u_idxs(nodes_offset)
    
    edges_orientation   = get_edges_orientation(nd.edges)
    edges_src_node      = first.(edges_orientation)
    edges_dst_node      = last.(edges_orientation)
    
    #--------------------------------------------------

    nodes_param_Idx = create_comps_non_flattened_vectors_Idx(nodes_param)

    #---------------------------------------------------
      

    edges_param_Idx = create_comps_non_flattened_vectors_Idx(edges_param)

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

    #-------------------------------------------------
    #--------------------------------------------------
    
    sd_nodes_param   = ( nodes_dyn_func, nodes_cb_sw, nodes_src_edges, nodes_dst_edges, edges_ih_Idx, edges_ik_Idx, nodes_f_t, nodes_param, nodes_state_Idx  )
    
    sd_edges_param  = ( edges_dyn_func, edges_dst_node, edges_src_node, nodes_u_Idx, edges_param, edges_state_Idx )
    
    
    return  (sd_nodes_param, sd_edges_param)

                                     
end 



function get_num_of_components(comp_collection)
    
    if isa(comp_collection, Union{Array, Vector})
        return length(collect(comp_collection)) 
    else
        return length(collect(values(comp_collection)))
    end    
end

