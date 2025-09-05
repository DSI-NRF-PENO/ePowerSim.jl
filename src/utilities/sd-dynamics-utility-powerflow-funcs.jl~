# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


####################################################
#---------------------------------------------------
# Begining  of power flow utility  
# with netd, nodes, edges collections as argument
#---------------------------------------------------
####################################################


####################################################
#  Nodes Categorisation
####################################################

function map_node_type(node)

    if string(node.Bus_type) == string(:Slack)
        node_type = 1
    elseif string(node.Bus_type) == string(:Generator)
        node_type = 2
    elseif string(node.Bus_type) == string(:EnergyStorage)
        node_type = 3         
    elseif string(node.Bus_type) == string(:Load)
        node_type = 4
    elseif string(node.Bus_type) == string(:Transmission)
        node_type = 5        
    else
        node_type = 6
    end
    
    return node_type
end

####################################################
# Nodes Type Test
####################################################


function test_slack_or_gen_node_type(node)

    return (string(node.Bus_type ) == string(:Slack) ||
        string(node.Bus_type ) ==
        string(:Generator)) ? true : false
end


function test_slack_node_type(node)

    return string(node.Bus_type ) ==
        string(:Slack) ? true : false
end

function test_gen_node_type(node)

    return string(node.Bus_type ) ==
        string(:Generator) ? true : false 
end


function test_load_node_type(node)

    return string(node.Bus_type ) ==
        string(:Load) ? true : false    
end


function test_transmission_node_type(node)

    return string(node.Bus_type ) ==
        string(:Transmission) ? true : false    
end


function test_energystorage_node_type(node)

    return string( node.Bus_type ) ==
        string(:EnergyStorage) ? true : false    
end


####################################################
# Nodes filters
####################################################

function filter_slack_or_gen_nodes(nodes)

    return filter(test_slack_or_gen_node_type,
                  collect(values(nodes)) )
end

function filter_slack_nodes(nodes)

    return filter(test_slack_node_type,
                  collect(values(nodes)) )
end

function filter_gen_nodes(nodes)

    return filter(test_gen_node_type,
                  collect(values(nodes)) )
    
end

function filter_load_nodes(nodes)

    return filter(test_load_node_type,
                  collect(values(nodes)) )
    
end


function filter_transmission_nodes(nodes)

    return filter(test_transmission_node_type,
                  collect(values(nodes)) )
    
end

function filter_energystorage_nodes(nodes)

    return filter(test_energystorage_node_type,
                  collect(values(nodes)) )
    
end

####################################################
# Accessors Functions
####################################################


function get_node_name(node)
 
    lens_bus_name=@optic _.name
 
    return getall(node, lens_bus_name)[1]
end


function get_node_Bus_name(node)
 
    lens_bus_name=@optic _.name
 
    return getall(node, lens_bus_name)[1]
end


function get_edge_name(edge)
 
    lens_edge_name=@optic _.name
 
    return getall(edge, lens_edge_name)[1]
end


function get_branch_name(branch)
 
    lens_branch_name=@optic _.name
 
    return getall(branch, lens_branch_name)[1]
end


function get_node_Bus_num(node)
 
    lens_bus_num=@optic _.Bus_num
 
    return getall(node, lens_bus_num)[1]
end

function get_node_S(node)
 
    lens_S=@optic _.S
 
    return getall(node, lens_S)[1]
end


function get_node_Pmax(node)
 
    lens_Pmax=@optic _.Pmax
 
    return getall(node, lens_Pmax)[1]
end


function get_node_Pmin(node)
 
    lens_Pmin=@optic _.Pmin
 
    return getall(node, lens_Pmin)[1]
end


function get_node_Qmax(node)
 
    lens_Qmax=@optic _.Qmax
 
    return getall(node, lens_Qmax)[1]
end


function get_node_Qmin(node)
 
    lens_Qmin=@optic _.Qmin
 
    return getall(node, lens_Qmin)[1]
end


function get_node_vmax(node)
 
    lens_vmax=@optic _.vmax
 
    return getall(node, lens_vmax)[1]
end


function get_node_vmin(node)
 
    lens_vmin=@optic _.vmin
 
    return getall(node, lens_vmin)[1]
end


function get_branch_orientation(branch)
    # A lens to get orientation field 
    lens=@optic _.orientation

    return getall(branch, lens)[1]
end


function get_branch_y(branch)
    # A lens to get orientation field 
    lens=@optic _.y

    return getall(branch, lens)[1]
end


function get_branch_y_shunt_km(branch)
    # A lens to get orientation field 
    lens=@optic _.y_shunt_km

    return getall(branch, lens)[1]
end


function get_branch_y_shunt_mk(branch)
    # A lens to get orientation field 
    lens=@optic _.y_shunt_mk

    return getall(branch, lens)[1]
end


function calc_branch_Ybr(
    y, y_shunt_km, y_shunt_mk, y_ratio, t_mk )

    return [(y + 1im*y_shunt_km)*1/(abs(y_ratio))^2  -y*1/conj(y_ratio); -y*1/y_ratio  y + 1im*y_shunt_mk ]
end


function get_branch_Ybr(branch)
    # A lens to get Ybr field 
    lens=@optic _.Ybr

    # Get Ybr field
    # Ybr = getall(branch, lens)[1]

    return getall(branch, lens)[1]
end


function get_branch_yff(branch)
    # A lens to get Ybr field 
    lens=@optic _.Ybr

    # Get Ybr field
    Ybr = getall(branch, lens)[1]

    return Ybr[1,1]
end


function get_branch_bff(branch)
    # A lens to get Ybr field 
    lens=@optic _.BbrDC

    # Get Ybr field
    BbrDC = getall(branch, lens)[1]

    return BbrDC[1,1]
end


function get_branch_yft(branch)
    # A lens to get Ybr field 
    lens=@optic _.Ybr

    # Get Ybr field
    Ybr = getall(branch, lens)[1]

    return Ybr[1,2]
end


function get_branch_bft(branch)
    # A lens to get Ybr field 
    lens=@optic _.BbrDC

    # Get Ybr field
    BbrDC = getall(branch, lens)[1]

    return BbrDC[1,2]
end

function get_branch_ytf(branch)
    # A lens to get Ybr field 
    lens=@optic _.Ybr

    # Get Ybr field
    Ybr = getall(branch, lens)[1]

    return Ybr[2,1]
end

function get_branch_btf(branch)
    # A lens to get Ybr field 
    lens=@optic _.BbrDC

    # Get Ybr field
    BbrDC = getall(branch, lens)[1]

    return BbrDC[2,1]
end

function get_branch_ytt(branch)
    # A lens to get Ybr field 
    lens=@optic _.Ybr

    # Get Ybr field
    Ybr = getall(branch, lens)[1]

    return Ybr[2,2]
end

function get_branch_btt(branch)
    # A lens to get Ybr field 
    lens=@optic _.BbrDC

    # Get Ybr field
    BbrDC = getall(branch, lens)[1]

    return BbrDC[2,2]
end

function get_branch_Pf_θshift(branch)
    # A lens to get Ybr field 
    lens_P_θshift = @optic _.Pshift

    P_θshift = getall(branch, lens_P_θshift)[1]

    Pf_θshift = P_θshift[1]
    
    return Pf_θshift
end

function get_branch_Pt_θshift(branch)
    # A lens to get Ybr field 
    lens_P_θshift = @optic _.Pshift

    P_θshift = getall(branch, lens_P_θshift)[1]

    Pt_θshift = P_θshift[2]
    
    return Pt_θshift
end

####################################################
# Indices and Property Tuples
####################################################


function get_node_Idx_and_type_numeric_tuple(node)
    node_Bus_num = get_node_Bus_num(node)
    
    if string( node.Bus_type ) == string(:Slack)
        node_Idx_and_type_numeric = (node_Bus_num, 1) 
    elseif string(node.Bus_type ) == string(:Generator)
        node_Idx_and_type_numeric = (node_Bus_num, 2)
    elseif string(node.Bus_type ) == string(:EnergyStorage)
        node_Idx_and_type_numeric = (node_Bus_num, 3)
    elseif string(node.Bus_type ) == string(:Load)
        node_Idx_and_type_numeric = (node_Bus_num, 4)
    elseif string(node.Bus_type ) == string(:Transmission)
        node_Idx_and_type_numeric = (node_Bus_num, 5)
    else
        node_Idx_and_type_numeric = (node_Bus_num, 6) 
    end
    
    return node_Idx_and_type_numeric
end


function get_node_Idx(node)
    # A lens to get power eld 
    lens_bus_num=@optic _.Bus_num
 
    return getall(node, lens_bus_num)[1] 
end


function get_node_Idx_and_type_tuple(node)
    # A lens to get power eld 
    lens_bus_num=@optic _.Bus_num
    lens_Bus_type=@optic _.Bus_type
 
    return (getall(node, lens_bus_num)[1],
            getall(node, lens_Bus_type)[1])
end

function get_slack_node_Idx_and_Vm_tuple(node)
    # A lens to get power eld 
    lens_bus_num=@optic _.Bus_num
    lens_Vm=@optic _.Vm

    return (getall(node, lens_bus_num)[1],
            getall(node, lens_Vm)[1] )
end

function get_slack_node_Idx_and_Vθ_tuple(node)
    # A lens to get power eld 
    lens_bus_num=@optic _.Bus_num
    lens_Vθ=@optic _.Vθ

    return (getall(node, lens_bus_num)[1],
            getall(node, lens_Vθ)[1] )
end

function get_gen_bus_Idx_and_power_tuple(node)
    # A lens to get power eld 
    lens_bus_num=@optic _.Bus_num
    lens_P=@optic _.P
    lens_Q=@optic _.Q

    return (getall(node, lens_bus_num)[1],
            getall(node, lens_P)[1] + 1im *
                getall(node, lens_Q)[1])
end

function get_gen_bus_Idx_and_Vm_tuple(node)
    # A lens to get power eld 
    lens_bus_num=@optic _.Bus_num
    lens_Vm=@optic _.Vm

    return (getall(node, lens_bus_num)[1],
            getall(node, lens_Vm)[1])
end

function get_gen_bus_Idx_and_Vθ_tuple(node)
    # A lens to get power eld 
    lens_bus_num=@optic _.Bus_num
    lens_Vm=@optic _.Vθ

    return (getall(node, lens_bus_num)[1],
            getall(node, lens_Vθ)[1])
end

function get_load_bus_Idx_and_power_tuple(node)
    # A lens to get power eld 
    lens_bus_num=@optic _.Bus_num
    lens_P=@optic _.P
    lens_Q=@optic _.Q

    return (getall(node,lens_bus_num)[1],
            getall(node,lens_P)[1] + 1im *
                getall(node,lens_Q)[1])
end

function get_shunt_Idx_and_y_shunt_tuple(shunt)
    # A lens to get orientation field 
    lens_bus_num=@optic _.Bus_num
    lens_y_shunt=@optic _.y_shunt

    return (getall(shunt, lens_bus_num)[1],
            getall(shunt, lens_y_shunt)[1])
end

function get_nodes_Idx_and_type_tuple(nodes)

    node_Idx_and_type_tuple =
        map(get_node_Idx_and_type_tuple,
            values(nodes))
 
    return node_Idx_and_type_tuple
end

function get_gens_bus_Idx_and_power_tuple(nodes)

    gen_nodes = filter_gen_nodes(nodes)
    
    gen_bus_Idx_and_power_tuple =
        map(get_gen_bus_Idx_and_power_tuple,
            gen_nodes)

    return gen_bus_Idx_and_power_tuple
end

function get_loads_bus_Idx_and_power_tuple(nodes)

    load_nodes = filter_load_nodes(nodes)
    
    load_bus_power_tuples =
        map(get_load_bus_Idx_and_power_tuple,
            load_nodes)

    return load_bus_power_tuples
end

function get_nodes_Bus_name(nodes)
    return map(get_node_Bus_name,
               collect(values(nodes)))
end


function get_branches_name(branches)
    return map(get_branch_name,
               collect(values(branches)))
end

####################################################
# Network Nodes Indices
####################################################

function get_slack_and_gen_nodes(nodes)
    
    return filter_slack_or_gen_nodes(nodes)
end



function get_slack_and_gen_nodes_Idx(nodes)

    slack_or_gen_nodes = filter_slack_or_gen_nodes(nodes)
        
    return map(get_node_Idx, slack_or_gen_nodes )
end


function get_slack_nodes_Idx(nodes)

    slack_nodes = filter_slack_nodes(nodes)
    
    slack_node_Idx_and_Vm_tuple = map(get_slack_node_Idx_and_Vm_tuple, slack_nodes)
    
    slack_nodes_Idx = first.(slack_node_Idx_and_Vm_tuple)

    return slack_nodes_Idx
end



function get_gens_bus_Idx(nodes)

    gen_nodes = filter_gen_nodes(nodes)
    
    gen_bus_Idx_Vm_tuple = map(get_gen_bus_Idx_and_Vm_tuple, gen_nodes)
    
    gens_bus_Idx = first.(gen_bus_Idx_Vm_tuple)

    return gens_bus_Idx
end

function get_loads_bus_Idx(nodes)

    load_nodes = filter_load_nodes(nodes)
    
    load_bus_power_tuple = map(get_load_bus_Idx_and_power_tuple, load_nodes)
    
    loads_bus_Idx = first.(load_bus_power_tuple)

    return loads_bus_Idx
end


function get_gen_nodes_Pmax(nodes)

    gen_nodes = filter_gen_nodes(nodes)
    
    return map(get_node_Pmax, gen_nodes)
end

function get_gen_nodes_Pmin(nodes)

    gen_nodes = filter_gen_nodes(nodes)
    
    return map(get_node_Pmin, gen_nodes)
end


function get_gen_nodes_Qmax(nodes)

    gen_nodes = filter_gen_nodes(nodes)
    
    return map(get_node_Qmax, gen_nodes)
end

function get_gen_nodes_Qmin(nodes)

    gen_nodes = filter_gen_nodes(nodes)
    
    return map(get_node_Qmin, gen_nodes)
end

function get_nodes_vmax(nodes)

    # load_nodes = filter_load_nodes(nodes)    
    return map(get_node_vmax, collect(values(nodes)) )
end

function get_nodes_vmin(nodes)

    # load_nodes = filter_load_nodes(nodes)    
    return  map(get_node_vmin, collect(values(nodes)) )
end


function get_nodes_S(nodes)

    # load_nodes = filter_load_nodes(nodes)    
    return  map(get_node_S, collect(values(nodes)) )
end


####################################################
# Network Nodal Properties Values
####################################################

function get_slack_nodes_Vm(nodes)

    slack_nodes = filter_slack_nodes(nodes)
    
    slack_node_Idx_and_Vm_tuple =
        map(get_slack_node_Idx_and_Vm_tuple, slack_nodes)
    
    slack_buses_Vm = last.(slack_node_Idx_and_Vm_tuple)

    return slack_buses_Vm
end


function get_gens_bus_num(nodes)

    gen_nodes = filter_gen_nodes(nodes)

    return map(get_node_Bus_num, gen_nodes)
end

function get_slack_nodes_Vθ(nodes)

    slack_nodes = filter_slack_nodes(nodes)
    
    slack_node_Idx_and_Vθ_tuple =
        map(get_slack_node_Idx_and_Vθ_tuple,
            slack_nodes)
    
    slack_buses_Vθ = last.(slack_node_Idx_and_Vθ_tuple)

    return slack_buses_Vθ
end

function get_gens_bus_Vm(nodes)

    gen_nodes = filter_gen_nodes(nodes)
    
    gen_bus_Idx_Vm_tuple =
        map(get_gen_bus_Idx_and_Vm_tuple, gen_nodes)
    
    gens_bus_Vm = last.(gen_bus_Idx_Vm_tuple)

    return gens_bus_Vm
end

function get_gens_bus_Vθ(nodes)

    gen_nodes = filter_gen_nodes(nodes)
    
    gen_bus_Idx_Vθ_tuple =
        map(get_gen_bus_Idx_and_Vθ_tuple, gen_nodes)
    
    gens_bus_Vθ = last.(gen_bus_Idx_Vθ_tuple)

    return gens_bus_Vθ
end

function get_gens_Sg(nodes)

    nodes_Bus_num =
        map(get_node_Bus_num,
            collect(values(nodes)) )

    num_v = maximum(nodes_Bus_num)    

    gens_bus_Idx_and_power_tuple =
        get_gens_bus_Idx_and_power_tuple(nodes)
    
    gens_bus_idx = first.(gens_bus_Idx_and_power_tuple)
    
    gens_power   = last.(gens_bus_Idx_and_power_tuple)

    Sg = sparsevec(gens_bus_idx, gens_power, num_v )

    return Sg
end


function get_gens_Sg(gen_nodes, net_size)

    gens_bus_Idx_and_power_tuple =
        get_gens_bus_Idx_and_power_tuple(gen_nodes)
    
    gens_bus_idx = first.(gens_bus_Idx_and_power_tuple)
    
    gens_power   = last.(gens_bus_Idx_and_power_tuple)

    Sg = sparsevec(gens_bus_idx, gens_power, net_size)

    return Sg
end

function get_loads_Sd(nodes)
    
    nodes_Bus_num =
        map(get_node_Bus_num, collect(values(nodes)) )

    num_v = maximum(nodes_Bus_num)
    
    loads_bus_power_tuple =
        get_loads_bus_Idx_and_power_tuple(nodes)

    loads_bus_idx = first.(loads_bus_power_tuple)
    
    loads_power   = last.(loads_bus_power_tuple)

    Sd = sparsevec(loads_bus_idx, loads_power, num_v )

    return Sd
end

function get_loads_Sd(load_nodes, net_size)

    loads_bus_power_tuple =
        get_loads_bus_Idx_and_power_tuple(load_nodes)

    loads_bus_idx = first.(loads_bus_power_tuple)
    
    loads_power   = last.(loads_bus_power_tuple)

    Sd = sparsevec(loads_bus_idx, loads_power, net_size)

    return Sd
end

# ---------------------------------------------


function get_slack_and_gens_Cg(nodes)

    num_nodes = length( collect(keys(nodes)) )
        
    slack_or_gen_nodes_Idx =
        get_slack_and_gen_nodes_Idx(nodes)

    num_slack_or_gen_nodes = length( slack_or_gen_nodes_Idx )
    
    Cg = sparse(slack_or_gen_nodes_Idx, collect(1:num_slack_or_gen_nodes), ones(Int,length( slack_or_gen_nodes_Idx )), num_nodes, num_slack_or_gen_nodes )

    return Cg
end


function get_gens_Cg(gen_nodes, num_nodes)

    gens_bus_idx  = get_gens_bus_num(gen_nodes)
    
    Cg = sparse(gens_bus_idx, gens_bus_idx, ones(Int,length(gens_bus_idx)), num_nodes, num_nodes)

    return Cg
end


function get_node_S(node)
 
    lens_S=@optic _.S
 
    return getall(node, lens_S)[1]
end


function get_nodes_S(nodes)
 
 
    return map(get_node_S, collect(values(nodes)) )
end


function get_node_Vm(node)
 
    lens_Vm=@optic _.Vm
 
    return getall(node, lens_Vm)[1]
end


function get_nodes_Vm(nodes)
 
    return map(get_node_Vm,
               collect(values(nodes)) )
end


function get_node_Vθ(node)
 
    lens_Vθ=@optic _.Vθ
 
    return getall(node, lens_Vθ)[1]
end


function get_nodes_Vθ(nodes)
 
    return map(get_node_Vθ,
               collect(values(nodes)) )
end


function get_nodes_Qmax(nodes)
 
    return map(get_node_Qmax,
               collect(values(nodes)) )
end


function get_nodes_Qmin(nodes)
 
    return map(get_node_Qmin,
               collect(values(nodes)) )
    
end



function get_slack_and_gen_Sg(nodes)

    slack_and_gen_nodes =
        get_slack_and_gen_nodes(nodes)
        
    return map(get_node_S,
               slack_and_gen_nodes)
end


function get_Sg(gen_nodes, num_nodes)

    gens_bus_idx =
        get_gens_bus_num(gen_nodes)
        
    gens_power =
        map(get_node_S,
            collect(values(gen_nodes)) )

    return sparsevec(
        gens_bus_idx,
        gens_power,
        num_nodes)
end


function node_Idx_from_name(node_name::String)

    return parse(Int,
                 split(lowercase(node_name),"bus")[2] )
end


function node_Idx_from_name(node)

    node_name = get_node_Bus_name(node)

    return parse(Int,
                 split(lowercase(node_name),
                       "bus")[2] )
end


"""
Get indices in slack nodes, gen nodes and load nodes order.
"""
function get_Idx_slack_gens_loads_order(
    slack_nodes, gen_nodes, load_nodes)

    names_nodes = unique(
        [get_nodes_Bus_name(Slack_nodes);
         get_nodes_Bus_name(Gen_nodes);
         get_nodes_Bus_name(Load_nodes)])

    nodes_Idx   =  node_Idx_from_name.(names_nodes)
    
    return nodes_Idx
end


function get_nodes_type_dims(
    slack_nodes,
    gen_nodes,
    load_nodes,
    branches)

    total_buses_type =
        get_network_nodal_size(branches)
    
    dim_slacks_buses =
        length(collect(values(slack_nodes)))
    
    dim_gens_buses =
        length(collect(values(gen_nodes)))
    
    dim_loads_buses =
        total_buses_type -
        (dim_slacks_buses + dim_gens_buses)
   
    return (dim_slacks_buses,
            dim_gens_buses,
            dim_loads_buses)
end


function get_nodes_type_Idx(
    slack_nodes,
    gen_nodes,
    load_nodes,
    branches)

    (dim_slacks_buses,
     dim_gens_buses,
     dim_loads_buses) =
         get_nodes_type_dims(
             slack_nodes,
             gen_nodes,
             load_nodes, branches)

    dims = [dim_slacks_buses,
            dim_gens_buses,
            dim_loads_buses]
    
    offset     = create_offsets(dims; counter=0)
    Idx        = create_idxs(offset, dims)
    slacks_Idx = Idx[1]
    gens_Idx   = Idx[2]
    loads_Idx  = Idx[3]

    gens_to_load_Idx = first(gens_Idx):last(loads_Idx)

    return (slacks_Idx,
            gens_Idx,
            loads_Idx,
            gens_to_load_Idx)
end

function get_Idx_permutation(nodes_idx)

    # sort a list of nodes_idx
    idx_perm = sortperm(nodes_idxs) 

    # get sorted list
    sorted_nodes_idx = nodes_idx[idx_perm] 
    
    perm_idx = Permutation(sorted_nodes_idx)

    perm_matrix = sparse(Matrix(perm_idx))
    
    return perm_matrix, perm_idx
end


####################################################
# Branches orientations
####################################################

function get_branches_orientations(branches)

    return map(get_branch_orientation,
               collect(values(branches)))
end

function get_branches_from(branches)

    return first.(map(get_branch_orientation,
                      collect(values(branches))))
end

function get_branches_to(branches)

    return last.(map(get_branch_orientation,
                     collect(values(branches))))
end



function get_branches_y(branches)

    return (map(get_branch_y,
                collect(values(branches))))
end



function get_branches_y_shunt_km(branches)

    return (map(get_branch_y_shunt_km,
                collect(values(branches))))
end



function get_branches_y_shunt_mk(branches)

    return (map(get_branch_y_shunt_mk,
                collect(values(branches))))
end


####################################################
# Incidence Matrices Formulations
# IncidenceMatricesFormulations.jl
####################################################

# (load_nodes, net_size)

function get_Cg(nodes)

    nodes_Bus_num = map(get_node_Bus_num,
                        collect(values(nodes)) )
    num_v = maximum(nodes_Bus_num)

    gen_nodes = filter_gen_nodes(nodes)
    
    gens_bus_idx = map(get_node_Bus_num, gen_nodes)
    
    # gens_bus_power_tuple = get_gens_bus_Idx_and_power_tuple(nodes)
    # gens_bus_idx = first.(gens_bus_power_tuple)

    Cg = sparse(gens_bus_idx, gens_bus_idx,
                ones(Int,length(gens_bus_idx)), num_v, num_v)

    return Cg
end

function get_Cg(gen_nodes, net_size)

    # gens_bus_power_tuple = get_gens_bus_Idx_and_power_tuple(gen_nodes)
    # gens_bus_idx = first.(gens_bus_power_tuple)
    
    gens_bus_idx = map(get_node_Bus_num,
                       collect(values(gen_nodes)))
    
    Cg = sparse(gens_bus_idx, gens_bus_idx,
                ones(Int,length(gens_bus_idx)),
                net_size, net_size)

    return Cg
end

function get_Cd(nodes)

    nodes_Bus_num = map(get_node_Bus_num,
                        collect(values(nodes)) )

    num_v = maximum(nodes_Bus_num)
    
    loads_bus_power_tuple =
        get_loads_bus_Idx_and_power_tuple(nodes)
    
    loads_bus_idx = first.(loads_bus_power_tuple)

    Cd = sparse(loads_bus_idx,
                loads_bus_idx,
                ones(Int,length(loads_bus_idx)),
                num_v, num_v)

    return Cd
end


function get_edge_orientation(edge)
    # A lens to get orientation field 
    lens=@optic _.orientation

    return getall(edge, lens)[1]
end


function get_edge_orientation(
    edge_number::Int64,
    edges_collection::Union{OrderedDict,Vector})

    return get_edges_orientation(edges_collection)[
        edge_number]
end



function get_node_src_edges(node_number::Int64, Cnb)

    return [find_src(Cnb[node_number,:])
            for a_row in
                node_number:node_number ]
    
end

function get_edges_orientation(edges_collection)
    
    if isa(edges_collection, Union{Array, Vector})
        return map(get_edge_orientation,
                   collect(edges_collection) )
    else
        return map(get_edge_orientation,
                   collect(values(edges_collection)))
    end  
end


function get_edges_orientation(
    edges_collection::Union{OrderedDict,Vector} )
    
    if isa(edges_collection, Union{Array, Vector})
        return map(get_edge_orientation,
                   collect(edges_collection) )
    else
        return map(get_edge_orientation,
                   collect(values(edges_collection)))
    end  
end


function get_edges_from(
    edges_collection::Union{OrderedDict,Vector})

    return first.(get_edges_orientation(edges_collection))
end


function get_edges_to(
    edges_collection::Union{OrderedDict,Vector})

    return last.(get_edges_orientation(
        edges_collection))
end


function get_network_nodal_size(edges_collection)
    
    edges_orientation =
        get_edges_orientation(edges_collection)
    
    src_v = first.(edges_orientation)
    
    dst_v = last.(edges_orientation)

    # net_size = maximum([src_v;dst_v])

    return maximum([src_v;dst_v])
end


function get_network_nodal_size(
    edges_collection::Union{OrderedDict,Vector})
    
    edges_orientation =
        get_edges_orientation(edges_collection)
    
    src_v = first.(edges_orientation)
    
    dst_v = last.(edges_orientation)

    # num_v    = maximum([src_v;dst_v])
    # net_size = maximum([src_v;dst_v])

    return maximum([src_v;dst_v])
end


function get_Cd(load_nodes, net_size)

    loads_bus_power_tuple =
        get_loads_bus_Idx_and_power_tuple(load_nodes)
    
    loads_bus_idx = first.(loads_bus_power_tuple)

    Cd = sparse(loads_bus_idx,
                loads_bus_idx,
                ones(Int,length(loads_bus_idx)),
                net_size, net_size)

    return Cd
end


function get_Cf(nodes, branches)
    
    src_v = get_branches_from(branches)
    
    dst_v = get_branches_to(branches)
    
    num_v = maximum([src_v;dst_v])
    
    from_orient = get_branches_from(branches)
    
    Cf = sparse(collect(1:length(branches)),
                from_orient, ones(Int,length(branches)),
                length(branches), num_v)
    return Cf
end


function get_Cf(edges_collection::Union{OrderedDict,Vector})
    
    edges_orientation =
        get_edges_orientation(edges_collection)
    
    src_v = first.(edges_orientation)
    
    dst_v = last.(edges_orientation)
    
    num_v = maximum([src_v;dst_v])
    
    Cf    = sparse(collect(1:length(src_v)),
                   src_v, ones(Int,length(src_v)),
                   length(src_v), num_v)
    return Cf
end



function get_Ct(nodes, branches)
    
    src_v = get_branches_from(branches)
    
    dst_v = get_branches_to(branches)
    
    num_v = maximum([src_v;dst_v])
    
    to_orient   = get_branches_to(branches)
    
    Ct = sparse(collect(1:length(branches)),
                to_orient, ones(Int,length(branches)),
                length(branches), num_v)
    return Ct
end


function get_Ct(edges_collection::Union{OrderedDict,Vector})
    
    edges_orientation =
        get_edges_orientation(edges_collection)
    
    src_v = first.(edges_orientation)
    
    dst_v = last.(edges_orientation)
    
    num_v = maximum([src_v;dst_v])
    
    Ct    = sparse(collect(1:length(dst_v)),
                   dst_v, ones(Int,length(dst_v)),
                   length(dst_v), num_v)
    return Ct
end



function get_C_branch_node(nodes, branches)
    
    src_v = get_branches_from(branches)
    
    dst_v = get_branches_to(branches)
    
    num_v = maximum([src_v;dst_v])

    from_orient = get_branches_from(branches)
    
    to_orient   = get_branches_to(branches)

    # nl x nb for from
    Cf = sparse(collect(1:length(branches)),
                from_orient, ones(Int,length(branches)),
                length(branches), num_v)

    # nl x nb for to
    Ct = sparse(collect(1:length(branches)),
                to_orient, ones(Int,length(branches)),
                length(branches), num_v )

    # Branch Node incidence matrix
    Cbn = sparse([collect(1:length(branches));
                  collect(1:length(branches))],
                 [from_orient;to_orient],
                 [ones(Int,length(branches));
                  -ones(Int,length(branches))],
                 length(branches), num_v )

    return Cbn
end

function get_Cnb(nodes, branches)
    
    Cf = get_Cf(nodes, branches)
    
    Ct = get_Ct(nodes, branches)
    
    Cbn = Cf - Ct

    return Cbn'
    
end



function get_Cbn(nodes, branches)
    
    Cf = get_Cf(nodes, branches)
    
    Ct = get_Ct(nodes, branches)
   
    return Cf - Ct
    
end


function get_Cnb(
    edges_collection::Union{OrderedDict,Vector})
    
    edges_orientation =
        get_edges_orientation(
            edges_collection)
    
    src_v = first.(edges_orientation)
    
    dst_v = last.(edges_orientation)
    
    num_v = maximum([src_v;dst_v])
    
    Cf    = sparse(collect(1:length(src_v)),
                   src_v, ones(Int,length(src_v)),
                   length(src_v), num_v)
    
    Ct    = sparse(collect(1:length(dst_v)),
                   dst_v, ones(Int,length(dst_v)),
                   length(dst_v), num_v)
    Cbn = Cf - Ct

    return Cbn'
    
end


function get_nodes_incident_edges(
    edges_collection::Union{OrderedDict,Vector})
    
    edges_orientation =
        get_edges_orientation(
            edges_collection)
    
    src_v = first.(edges_orientation)
    
    dst_v = last.(edges_orientation)
    
    num_v = maximum([src_v;dst_v])
        
    # num_v = length(collect(keys(nodes)))
    
    Cnb = get_Cnb(edges_collection)

    return [find_node_incident_edges(Cnb[a_row,:])
            for a_row in 1:num_v ]
    
end


function get_nodes_src_edges(nodes, branches)

    src_v = get_branches_from(branches)
    dst_v = get_branches_to(branches)
    
    num_v = maximum([src_v;dst_v])
    
    # num_v = length(collect(keys(nodes)))
    
    Cnb = get_Cnb(nodes, branches)

    return [find_src(Cnb[a_row,:])
            for a_row in 1:num_v ]
    
end


function get_nodes_src_edges(
    edges_collection::Union{OrderedDict,Vector})
    
    edges_orientation =
        get_edges_orientation(
            edges_collection)
    
    src_v = first.(edges_orientation)
    
    dst_v = last.(edges_orientation)
    
    num_v = maximum([src_v;dst_v])
        
    # num_v = length(collect(keys(nodes)))
    
    Cnb = get_Cnb(edges_collection)

    return [find_src(Cnb[a_row,:])
            for a_row in 1:num_v ]
    
end




function get_nodes_dst_edges(nodes, branches)
    
    src_v = get_branches_from(branches)
    
    dst_v = get_branches_to(branches)
    
    num_v = maximum([src_v;dst_v])
    # num_v = length(collect(keys(nodes)))
    
    Cnb = get_Cnb(nodes, branches)

    return [find_dst(Cnb[a_row,:])
            for a_row in 1:num_v ]
end


function get_nodes_dst_edges(
    edges_collection::Union{OrderedDict,Vector})
    
    edges_orientation =
        get_edges_orientation(edges_collection)
    
    src_v = first.(edges_orientation)
    
    dst_v = last.(edges_orientation)
    
    num_v = maximum([src_v;dst_v])
        
    # num_v = length(collect(keys(nodes)))
    
    Cnb = get_Cnb(edges_collection)

    return [find_dst(Cnb[a_row,:])
            for a_row in 1:num_v ]
end




function get_node_src_edges(node, Cnb)

    node_number = get_node_Bus_num(node) 

    return [find_src(Cnb[node_number,:])
            for a_row in
                node_number:node_number ]
    
end


function get_node_dst_edges(node, Cnb)

    node_number = get_node_Bus_num(node) 

    return [find_dst(Cnb[node_number,:])
            for a_row in
                node_number:node_number ]
    
end


# function get_node_src_edges(node_number::Int64, Cnb)

#     return [find_src(Cnb[node_number,:])
#             for a_row in
#                 node_number:node_number ]
    
# end


# function get_node_dst_edges(node_number::Int64, Cnb)

#     return [find_dst(Cnb[node_number,:])
#             for a_row in
#                 node_number:node_number ]
    
# end


####################################################
# Admittances Formulations
# AdmittancesFormulations.jl
####################################################

function get_branches_Yff(branches)

    Yff = map(get_branch_yff, collect(values(branches)))

    return Yff # Diagonal(Yff)
end


function get_branches_Yff(
    edges_collection::Union{OrderedDict,Vector})

    Yff = map(get_branch_yff,
              collect(values(edges_collection)))

    return Yff # Diagonal(Yff)
end



function get_branches_Bff(branches)

    Bff = map(get_branch_bff,
              collect(values(branches)))

    return Bff # Diagonal(Bff)
end

function get_branches_Yft(branches)

    Yft = map(get_branch_yft,
              collect(values(branches)))

    return Yft # Diagonal(Yft)
end


function get_branches_Yft(
    edges_collection::Union{OrderedDict,Vector})

    Yft = map(get_branch_yft,
              collect(values(edges_collection)))

    return Yft # Diagonal(Yft)
end



function get_branches_Bft(branches)

    Bft = map(get_branch_bft,
              collect(values(branches)))

    return Bft # Diagonal(Bft)
end

function get_branches_Ytf(branches)

    Ytf = map(get_branch_ytf,
              collect(values(branches)))

    return Ytf # Diagonal(Ytf)
end



function get_branches_Ytf(
    edges_collection::Union{OrderedDict,Vector})

    Ytf = map(get_branch_ytf,
              collect(values(edges_collection)))

    return Ytf # Diagonal(Ytf)
end


function get_branches_Btf(branches)

    Btf = map(get_branch_btf,
              collect(values(branches)))

    return Btf # Diagonal(Btf)
end


function get_branches_Ytt(branches)

    Ytt = map(get_branch_ytt,
              collect(values(branches)))

    return Ytt # Diagonal(Ytt)
end


function get_branches_Ytt(
    edges_collection::Union{OrderedDict,Vector})

    Ytt = map(get_branch_ytt,
              collect(values(edges_collection)))

    return Ytt # Diagonal(Ytt)
end



function get_branches_Btt(branches)

    Btt = map(get_branch_btt,
              collect(values(branches)))

    return Btt # Diagonal(Btt)
end


function get_branches_Pf_θshift(branches)

    Pf_θshift = map(get_branch_Pf_θshift,
                    collect(values(branches)))

    return Pf_θshift
end


function get_branches_Pt_θshift(branches)

    Pt_θshift = map(get_branch_Pt_θshift,
                    collect(values(branches)))

    return Pt_θshift
end


function get_Pbus_θshift(nodes, branches)

    Cf = get_Cf(nodes, branches)
    Ct = get_Ct(nodes, branches)
    Pf_θshift = get_branches_Pf_θshift(branches)
    

    return (Cf - Ct)' * Pf_θshift
end

function get_Yf(nodes, branches)
    
    Cf  = get_Cf(nodes, branches)
    
    Ct  = get_Ct(nodes, branches)
    
    Yff = get_branches_Yff(branches)
    
    Yft = get_branches_Yft(branches)

    # Yf = Yff * Cf + Yft * Ct
    Yf = Diagonal(Yff) * Cf +
        Diagonal(Yft) * Ct
    
    return Yf
end


function get_Yf(
    edges_collection::Union{OrderedDict,Vector})
    
    Cf  = get_Cf(edges_collection)
    
    Ct  = get_Ct(edges_collection)
    
    Yff = get_branches_Yff(edges_collection)
    
    Yft = get_branches_Yft(edges_collection)

    # Yf = Yff * Cf + Yft * Ct
    Yf = Diagonal(Yff) * Cf +
        Diagonal(Yft) * Ct
    
    return Yf
end




function get_Bf(nodes, branches)
    
    Cf  = get_Cf(nodes, branches)
    
    Ct  = get_Ct(nodes, branches)
    
    Bff = get_branches_Bff(branches)
    
    Bft = get_branches_Bft(branches)

    # Yf = Yff * Cf + Yft * Ct
    Bf = Diagonal(Bff) * ( Cf -  Ct)
    return Bf
end


function get_Yt(nodes, branches)
    
    Cf  = get_Cf(nodes, branches)
    
    Ct  = get_Ct(nodes, branches)
    
    Ytf = get_branches_Ytf(branches)
    
    Ytt = get_branches_Ytt(branches)

    # Yt = Ytf * Cf + Ytt * Ct
    Yt = Diagonal(Ytf) * Cf +
        Diagonal(Ytt) * Ct
    
    return Yt
end


function get_Yt(
    edges_collection::Union{OrderedDict,Vector})
    
    Cf  = get_Cf(edges_collection)
    
    Ct  = get_Ct(edges_collection)
    
    Ytf = get_branches_Ytf(edges_collection)
    
    Ytt = get_branches_Ytt(edges_collection)

    # Yt = Ytf * Cf + Ytt * Ct
    Yt = Diagonal(Ytf) * Cf +
        Diagonal(Ytt) * Ct
    
    return Yt
end



function get_Bt(nodes, branches)
    
    Cf  = get_Cf(nodes, branches)
    
    Ct  = get_Ct(nodes, branches)
    
    Btf = get_branches_Btf(branches)
    
    Btt = get_branches_Btt(branches)

    # Yt = Ytf * Cf + Ytt * Ct
    Bt = Diagonal(Btf) * Cf +
        Diagonal(Btt) * Ct
    
    return Bt
end


function get_nodal_Ysh(shunts, nodes)

    nodes_Bus_num = map(get_node_Bus_num,
                        collect(values(nodes)) )
    
    num_v = maximum(nodes_Bus_num)
    
    # nodes_num = length(nodes)
    yshs      =
        map(get_shunt_Idx_and_y_shunt_tuple,
            values(shunts))
    
    bus_num   = first.(yshs)
    
    bus_ysh   = last.(yshs)
    
    ysh_sparse = sparsevec(
        bus_num, bus_ysh, num_v)
    
    # Ysh = spdiagm(ysh_sparse)

    return ysh_sparse # Ysh
end


function get_nodal_Bsh(shunts, nodes)
    
    nodes_Bus_num =
        map(get_node_Bus_num,
            collect(values(nodes)) )
    
    num_v = maximum(nodes_Bus_num)
    
    # nodes_num = length(nodes)
    yshs =
        map(get_shunt_Idx_and_y_shunt_tuple,
            values(shunts))
    
    bus_num = first.(yshs)
    
    bus_ysh = last.(yshs)
    
    bus_bsh = imag.(bus_ysh)
    
    bsh_sparse = sparsevec(
        bus_num, bus_bsh, num_v)
    
    # Ysh = spdiagm(ysh_sparse)

    return bsh_sparse # Ysh
end


function get_nodal_Gsh(shunts, nodes)
        
    nodes_Bus_num =
        map(get_node_Bus_num,
            collect(values(nodes)) )
    
    num_v = maximum(nodes_Bus_num)
    
    # nodes_num = length(nodes)
    yshs = map(
        get_shunt_Idx_and_y_shunt_tuple,
        values(shunts))
    
    bus_num = first.(yshs)
    
    bus_ysh = last.(yshs)
    
    bus_gsh = real.(bus_ysh)
    
    gsh_sparse = sparsevec(bus_num, bus_gsh, num_v)
    
    # Ysh = spdiagm(ysh_sparse)

    return gsh_sparse # Ysh
end


function get_Psh(shunts, nodes)
    
    nodes_Bus_num =
        map(get_node_Bus_num,
            collect(values(nodes)) )
    
    num_v = maximum(nodes_Bus_num)
    
    # nodes_num = length(nodes)
    yshs =
        map(get_shunt_Idx_and_y_shunt_tuple,
            values(shunts))
    
    bus_num = first.(yshs)
    
    bus_ysh = last.(yshs)
    
    bus_gsh = real.(bus_ysh)
    
    gsh_sparse = sparsevec(bus_num, bus_gsh, num_v)
    
    # Ysh = spdiagm(ysh_sparse)

    return gsh_sparse # Ysh
end

function get_Ybus(nodes, branches; shunts=nothing)
    
    Cf  = get_Cf(nodes, branches)
    
    Ct  = get_Ct(nodes, branches)
    
    Yf  = get_Yf(nodes, branches)
    
    Yt  = get_Yt(nodes, branches)
    
    if shunts != nothing
        Ysh = get_nodal_Ysh(shunts, nodes)
        
        Ybus = Cf' * Yf + Ct' * Yt + spdiagm(Ysh) # Ysh
    else
        Ybus = Cf' * Yf + Ct' * Yt 
    end

    return Ybus
end

function get_Bbus(nodes, branches)
    
    Cf  = get_Cf(nodes, branches)
    
    Ct  = get_Ct(nodes, branches)
    
    Bf  = get_Bf(nodes, branches)

    Bbus = (Cf - Ct)' * Bf
    return Bbus
end

####################################################
# Permutation matrix
####################################################

function get_sorted_nodes_idx_and_type(nodes)
    
    nodes_idx_and_type =
        map(get_node_Idx_and_type_numeric_tuple,
            values(nodes))

    # sort a list of (bus_index, bus_type) tupels by bus_type
    idx_perm_nodes_idx_and_type =
        sortperm(nodes_idx_and_type, by = last)

    # get sorted list 
    sorted_nodes_idx_and_type =
        nodes_idx_and_type[idx_perm_nodes_idx_and_type]
    
    # get indices of sorted_nodes
    idx_bus_sorted_by_type =
        map(x -> first(x),
            sorted_nodes_idx_and_type)

    perm_idx = Permutation(idx_bus_sorted_by_type)

    perm_matrix = sparse(Matrix(perm_idx))
    
    return perm_matrix, sorted_nodes_idx_and_type
end


function get_nodes_permutation(nodes)

    nodes_idx_and_type =
        map(get_node_Idx_and_type_numeric_tuple,
            values(nodes))

    # sort a list of (bus_index, bus_type) tupels by bus_type
    idx_perm_nodes_idx_and_type =
        sortperm(nodes_idx_and_type, by = last)

    # get sorted list 
    sorted_nodes_idx_and_type =
        nodes_idx_and_type[
            idx_perm_nodes_idx_and_type]

    # get the index of the first bus with type == 1 i.e gen node
    idx_first_source_bus =
        findfirst(x -> last(x) == 2,
                  sorted_nodes_idx_and_type)

    
    # get the index of the first bus with type == 4 i.e demand node
    idx_first_energystorage_bus =
        findfirst(x -> last(x) == 3,
                  sorted_nodes_idx_and_type)
    
    # get the index of the first bus with type == 4 i.e demand node
    idx_first_demand_bus =
        findfirst(x -> last(x) == 4,
                  sorted_nodes_idx_and_type)

    # range of slack nodes indices
    idx_slack_nodes = 1:idx_first_source_bus-1
    
    # range of source nodes indices
    idx_source_nodes =
        idx_first_source_bus:idx_first_demand_bus-1

     # range of demand nodes indices
    idx_demand_nodes =
        idx_first_demand_bus:length(sorted_nodes_idx_and_type)

    # get indices of sorted_nodes
    idx_bus_sorted_by_type =
        map(x -> first(x),
            sorted_nodes_idx_and_type)

    perm_idx = Permutation(idx_bus_sorted_by_type)

    perm_matrix = sparse(Matrix(perm_idx))

    return (perm_matrix,
            idx_perm_nodes_idx_and_type,
            idx_slack_nodes,
            idx_source_nodes,
            idx_demand_nodes)

end

####################################################
#---------------------------------------------------
# End  of power flow utility
#---------------------------------------------------
####################################################
