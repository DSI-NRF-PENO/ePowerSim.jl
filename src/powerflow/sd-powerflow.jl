# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)



####################################################
# Untested functions
####################################################

function V_R2C(V)
    size_Θ_Vm = length(X)
    N = Int(size_Θ_Vm/2)

    Vm = @view V[1:N]
    Θ  = @view V[N+1:end]    
    Vc = Vm .* exp.(im * Θ)

    return Vc
 
end

function V_C2R(V)

    Vc   = @view V[1:end]

    Θ   = angle.(Vc)
    Vm  = abs.(Vc)

    return (Vm, Θ)
 
end

function Sf( V, Cf, Yf)

    Sf = Diagonal(Cf * Vf) * conj.(Yf) * conj.(V)

    return Sf
end

function St(V, Ct, Yt)

    St = Diagonal(Ct * Vt) * conj.(Yt) * conj.(V)

    return St
end

function Sbus(V, Ybus)
    size_Θ_Vm = length(X)
    N = Int(size_Θ_Vm/2)

    Vm   = @view V[1:N]
    Θ    = @view V[N+1:end]
    
    Vc   = Vm .* exp.(im * Θ)
    Sbus = Diagonal(Vc) * conj.(Ybus) * conj.(Vc)

    return Sbus
end

function SbusC(V, Ybus)

    return Diagonal(V) * conj.(Ybus) * conj.(V)
end

function gS(V, p)

    Cg, Sg, Sd, Ybus  = p
    
    return Diagonal(V) * conj.(Ybus) * conj.(V) + Sd

end

########################################################
# Utility Functions
# UtilityFunctions.jl
########################################################


function z2y(; r =1.0,  x = 1.00, G = 0.0,  B = 0.0)

    z = r + im * x

    y = 1/z

    B_2 = B/2.0
    

    return (y = y,  y_shunt_km = B_2, y_shunt_mk = B_2)

end


function DAE_MassMatrix(state_size, algebraic_size)
    return Diagonal([ones(Int, state_size)...; zeros(Int, algebraic_size)...])
end

function threshold_limits(x, x_max, x_min)
    return x > x_max ? x_max : x < x_min  ? x_min : x
end

"""
Create offsets for stacked array of dimensions dims
"""
function create_offsets(dims; counter=0)::Vector{Int}
    offs = [1 for dim in dims]
    for (i, dim) in enumerate(dims)
        offs[i] = counter
        counter += dim
    end
    offs
end

"""
    create_idxs(offs, dims)

Create indexes for stacked array of dimensions dims
using the offsets offs
"""
function create_idxs(offs, dims)::Vector{UnitRange{Int}}
    idxs = [1+off:off+dim
            for (off, dim) in
                zip(offs, dims)]
end


 Nodes Categorisation


function map_node_type(node)

    if string(typeof(node))     == string(:Slack)
        node_type = 1
    elseif string(typeof(node)) == string(:Generator)
        node_type = 2
    elseif string(typeof(node)) == string(:EnergyStorage)
        node_type = 3         
    elseif string(typeof(node)) == string(:Load)
        node_type = 4
    else
        node_type = 5
    end
    
    return node_type
end


function map_node_type(node)

    if string(node.Bus_type)      == string(:Slack)
        node_type = 1
    elseif string(node.Bus_type)  == string(:Generator)
        node_type = 2
    elseif string(node.Bus_type)  == string(:EnergyStorage)
        node_type = 3         
    elseif string(node.Bus_type)  == string(:Load)
        node_type = 4
    else
        node_type = 5
    end
    
    return node_type
end

####################################################
# Nodes Type Test
####################################################

# function test_slack_node_type(node)

#     return string(typeof(node)) == string(:Slack) ? true : false
# end

# function test_gen_node_type(node)

#     return string(typeof(node)) == string(:Generator) ? true : false 
# end


# function test_load_node_type(node)

#     return string(typeof(node)) == string(:Load) ? true : false    
# end


# function test_energystorage_node_type(node)

#     return string(typeof(node)) == string(:EnergyStorage) ? true : false    
# end



# function test_mystring(mystring)

#     return (mystring == string(:Slack) ||  mystring == string(:Generator))  ? true : false
# end

# test_mystring("Slack")
# test_mystring("Generator")


function test_slack_or_gen_node_type(node)

    return (string(node.Bus_type ) == string(:Slack) || string(node.Bus_type ) == string(:Generator))  ? true : false
end


function test_slack_node_type(node)

    return string(node.Bus_type ) == string(:Slack) ? true : false
end

function test_gen_node_type(node)

    return string(node.Bus_type ) == string(:Generator) ? true : false 
end


function test_load_node_type(node)

    return string(node.Bus_type ) == string(:Load) ? true : false    
end


function test_energystorage_node_type(node)

    return string( node.Bus_type ) == string(:EnergyStorage) ? true : false    
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

    return [(y + 1im*y_shunt_km)*1/(abs(y_ratio))^2  -y*1/conj(y_ratio); -y*1/y_ratio y + 1im*y_shunt_mk ]
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

# function get_node_Idx_and_type_numeric_tuple(node)
#     node_Bus_num = get_node_Bus_num(node)
    
#     if string(typeof(node))     == string(:Slack)
#         node_Idx_and_type_numeric = (node_Bus_num, 1) 
#     elseif string(typeof(node)) == string(:Generator)
#         node_Idx_and_type_numeric = (node_Bus_num, 2)
#     elseif string(typeof(node)) == string(:EnergyStorage)
#         node_Idx_and_type_numeric = (node_Bus_num, 3)         
#     elseif string(typeof(node)) == string(:Load)
#         node_Idx_and_type_numeric = (node_Bus_num, 4) 
#     else
#         node_Idx_and_type_numeric = (node_Bus_num, 5) 
#     end
    
#     return node_Idx_and_type_numeric
# end


function get_node_Idx_and_type_numeric_tuple(node)
    
    node_Bus_num = get_node_Bus_num(node)
    
    if string( node.Bus_type )     == string(:Slack)
        node_Idx_and_type_numeric = (node_Bus_num, 1) 
    elseif string(node.Bus_type ) == string(:Generator)
        node_Idx_and_type_numeric = (node_Bus_num, 2)
    elseif string(node.Bus_type ) == string(:EnergyStorage)
        node_Idx_and_type_numeric = (node_Bus_num, 3)        
    elseif string(node.Bus_type ) == string(:Load)
        node_Idx_and_type_numeric = (node_Bus_num, 4) 
    else
        node_Idx_and_type_numeric = (node_Bus_num, 5) 
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
            getall(node, lens_P)[1] +
                1im * getall(node, lens_Q)[1])
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
            getall(node,lens_P)[1] +
                1im * getall(node,lens_Q)[1])
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
 return map(get_node_Bus_name, collect(values(nodes)))
end


function get_branches_name(branches)
 return map(get_branch_name, collect(values(branches)))
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
    
    slack_node_Idx_and_Vm_tuple =
        map(get_slack_node_Idx_and_Vm_tuple, slack_nodes)
    
    slack_nodes_Idx =
        first.(slack_node_Idx_and_Vm_tuple)

    return slack_nodes_Idx
end



function get_gens_bus_Idx(nodes)

    gen_nodes = filter_gen_nodes(nodes)
    
    gen_bus_Idx_Vm_tuple =
        map(get_gen_bus_Idx_and_Vm_tuple,
            gen_nodes)
    
    gens_bus_Idx = first.(gen_bus_Idx_Vm_tuple)

    return gens_bus_Idx
end

function get_loads_bus_Idx(nodes)

    load_nodes = filter_load_nodes(nodes)
    
    load_bus_power_tuple =
        map(get_load_bus_Idx_and_power_tuple,
            load_nodes)
    
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
        map(get_slack_node_Idx_and_Vm_tuple,
            slack_nodes)
    
    slack_buses_Vm =
        last.(slack_node_Idx_and_Vm_tuple)

    return slack_buses_Vm
end


function get_gens_bus_num(nodes)

    gen_nodes = filter_gen_nodes(nodes)

    return return map(get_node_Bus_num, gen_nodes)
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
        map(get_node_Bus_num, collect(values(nodes)) )

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


function get_Cnb_by_orientations(edges_orientations)
    
    src_v = first.(edges_orientations)
    
    dst_v = last.(edges_orientations)
    
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


function get_Cbn_by_orientations(edges_orientations)
    
    src_v = first.(edges_orientations)
    
    dst_v = last.(edges_orientations)
    
    num_v = maximum([src_v;dst_v])
    
    Cf    = sparse(collect(1:length(src_v)),
                   src_v, ones(Int,length(src_v)),
                   length(src_v), num_v)
    
    Ct    = sparse(collect(1:length(dst_v)),
                   dst_v, ones(Int,length(dst_v)),
                   length(dst_v), num_v)
    Cbn = Cf - Ct

    return Cbn
    
end


find_src(row) = findall((x)->x==1,  row)

find_dst(row) = findall((x)->x==-1, row)


find_node_incident_edges(row) =
    findall((x)->(x==1||x==-1), row)


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


function get_nodes_incident_edges_by_orientations(
    edges_orientations )
        
    src_v = first.(edges_orientations)
    
    dst_v = last.(edges_orientations)
    
    num_v = maximum([src_v;dst_v])
        
    # num_v = length(collect(keys(nodes)))
    
    Cnb = get_Cnb_by_orientations(
        edges_orientations)

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


function get_node_src_edges(node_number::Int64, Cnb)

    return [find_src(Cnb[node_number,:])
            for a_row in
                node_number:node_number ]
    
end


function get_node_dst_edges(node_number::Int64, Cnb)

    return [find_dst(Cnb[node_number,:])
            for a_row in
                node_number:node_number ]
    
end


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
# Network Balance Formulations
# NetworkBalanceFormulations.jl
####################################################

function DC_nodal_power_balance(θ, p)

    (Cg,
     Pg, Pd,
     Gsh, Bbus,
     Pbus_θshift,
     PV_Idx,
     PQ_Idx,
     perm_matrix) = p
    
    Idx_PV_PQ  = [PV_Idx;PQ_Idx]
    
    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)

    vars_dim    = [N_PV_θ, N_PQ_θ]
    vars_offset = create_offsets(vars_dim)
    vars_Idx    = create_idxs(vars_offset, vars_dim)
   
    PV_θ       = @view θ[vars_Idx[1]]
    PQ_θ       = @view θ[vars_Idx[2]]

    θdc        =  [[PV_θ]... ; [PQ_θ]...]

    Ac         =  perm_matrix' * Bbus * perm_matrix 
    Ac         =  Ac[Idx_PV_PQ,Idx_PV_PQ]

    bc         = perm_matrix' * (-Pbus_θshift - Pd - Gsh + Cg * Pg)
    bc         = bc[Idx_PV_PQ]

    return Ac, bc
end 

# -----------------------------------------------
# page 250, milano (10.2)

function AC_nodal_power_balance(dV, V, p)

    (nodes_vmax,
     nodes_vmin,
     Cg,
     Sg,
     Sd,
     Ybus,
     gen_Vm,
     slack_Vm,
     slack_Vθ,
     PV_Idx,
     PQ_Idx,
     perm_matrix) = p
    
    Idx_PV_PQ  =
        first(PV_Idx):last(PQ_Idx) # [PV_Idx; PQ_Idx]

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    vars_dim    = [N_PV_θ, N_PQ_θ, N_PQ_Vm]
    vars_offset = create_offsets(vars_dim)
    vars_Idx    = create_idxs(vars_offset, vars_dim)
    
    PV_θ       = @view V[vars_Idx[1]]
    PQ_θ       = @view V[vars_Idx[2]]
    PQ_Vm      = @view V[vars_Idx[3]]

    dgP        = @view dV[1:length(Idx_PV_PQ)]
    dgQ        = @view dV[1+length(Idx_PV_PQ):end]
    
 
    Vnodes = [[slack_Vm .* exp.(im * slack_Vθ)]...;[gen_Vm .* exp.(im*PV_θ)]... ; [PQ_Vm .* exp.(im * PQ_θ)]...]

    Vabs = abs.(Vnodes)
    Vabs = threshold_limits.(Vabs, nodes_vmax, nodes_vmin)
    Vnodes = Vabs .* exp.(im * angle.(Vnodes))     

    # Create Index set to access different portion of Vnodes array
    Vnodes_dim          = [length(slack_Vm), length(PV_Idx), length(PQ_Idx)]
    Vnodes_offset       = create_offsets(Vnodes_dim)
    Vnodes_Idx          = create_idxs(Vnodes_offset, Vnodes_dim)
    
    Slack_Idx_in_Vnodes = Vnodes_Idx[1]
    PV_V_Idx_in_Vnodes  = Vnodes_Idx[2]
    PQ_V_Idx_in_Vnodes  = Vnodes_Idx[3]

    # PV_θ  .= angle.(Vnodes[PV_V_Idx_in_Vnodes])
    # PQ_θ  .= angle.(Vnodes[PQ_V_Idx_in_Vnodes])
    # PQ_Vm .= abs.(Vnodes[PQ_V_Idx_in_Vnodes])

    ###
    # gS =  (Diagonal(Vnodes) * conj.(Ybus) * conj.(Vnodes) + Sd -  Cg * Sg)
    # gS = perm_matrix * Diagonal(Vnodes) * perm_matrix'  * conj.(Ybus) * perm_matrix * conj.(Vnodes) + Sd -  Cg * Sg

    gS =  Diagonal(Vnodes) * (perm_matrix' * conj.(Ybus) * perm_matrix * conj.(Vnodes)) + perm_matrix' * (Sd - Cg * Sg)
    
    gP = real.(gS)
    gQ = imag.(gS)
    
    dgP   .= gP[Idx_PV_PQ]
    dgQ   .= gQ[PQ_Idx]

    return nothing
   
end

# ---------------------------------------------------

function AC_nodal_power_balance(V, p)

    nodes_vmax, nodes_vmin, Cg, Sg, Sd, Ybus, gen_Vm, slack_Vm, slack_Vθ, PV_Idx, PQ_Idx, perm_matrix = p
    
    Idx_PV_PQ  =  first(PV_Idx):last(PQ_Idx) # [PV_Idx; PQ_Idx]

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    vars_dim    = [N_PV_θ, N_PQ_θ, N_PQ_Vm]
    vars_offset = create_offsets(vars_dim)
    vars_Idx    = create_idxs(vars_offset, vars_dim)
    
    PV_θ       = @view V[vars_Idx[1]]
    PQ_θ       = @view V[vars_Idx[2]]
    PQ_Vm      = @view V[vars_Idx[3]]

    Vnodes = [[slack_Vm .* exp.(im * slack_Vθ)]...;
              [gen_Vm .* exp.(im*PV_θ)]... ;
              [PQ_Vm .* exp.(im * PQ_θ)]...]

    Vabs   = abs.(Vnodes)
    
    Vabs = threshold_limits.(Vabs, nodes_vmax, nodes_vmin)
    
    Vnodes = Vabs .* exp.(im * angle.(Vnodes))     

    # Create Index set to access different
    # portion of Vnodes array
    
    Vnodes_dim = [length(slack_Vm),
                  length(PV_Idx),
                  length(PQ_Idx)]
    
    Vnodes_offset = create_offsets(Vnodes_dim)
    
    Vnodes_Idx = create_idxs(Vnodes_offset, Vnodes_dim)
    
    Slack_Idx_in_Vnodes = Vnodes_Idx[1]
    
    PV_V_Idx_in_Vnodes  = Vnodes_Idx[2]
    
    PQ_V_Idx_in_Vnodes  = Vnodes_Idx[3]
    
    PQ_Vm  .= abs.(Vnodes[PQ_V_Idx_in_Vnodes])

    ##
    # gS =  Diagonal(Vnodes) * conj.( Ybus) * conj.(Vnodes) + Sd -  Cg * Sg
    # gS = perm_matrix * Diagonal(Vnodes) * perm_matrix'  * conj.(Ybus) * perm_matrix * conj.(Vnodes) + Sd -  Cg * Sg

    gS = Diagonal(Vnodes) *
        (perm_matrix' * conj.(Ybus) *
        perm_matrix * conj.(Vnodes)) +
        perm_matrix' * (Sd - Cg * Sg)
    
    gP = real.(gS)
    gQ = imag.(gS)
    
    [gP[Idx_PV_PQ]; gQ[PQ_Idx]]    
end


function AC_nodal_current_balance(dV, V, p)

    nodes_vmax, nodes_vmin, Cg, Sg, Sd, Ybus, gen_Vm, slack_Vm, slack_Vθ, PV_Idx, PQ_Idx, perm_matrix = p

    PV_PQ_Idx  = first(PV_Idx):last(PQ_Idx) # [PV_Idx;PQ_Idx]
    
    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PV_Vm    = length(PV_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    vars_dim    = [N_PV_θ, N_PQ_θ, N_PV_Vm, N_PQ_Vm]
    vars_offset = create_offsets(vars_dim)
    vars_Idx    = create_idxs(vars_offset, vars_dim)
    
    PV_θ       = @view V[vars_Idx[1]]
    PQ_θ       = @view V[vars_Idx[2]]
    PV_Vm      = @view V[vars_Idx[3]]
    PQ_Vm      = @view V[vars_Idx[4]]

    PV_Vm      .= gen_Vm 

    dgM        = @view dV[1:length(PV_PQ_Idx)]
    dgN        = @view dV[1+length(PV_PQ_Idx):end]

    Vnodes = [[slack_Vm .* exp.(im * slack_Vθ)]...;[PV_Vm .* exp.(im*PV_θ)]... ; [PQ_Vm .* exp.(im * PQ_θ)]...]

    Vabs   = abs.(Vnodes)
    Vabs   = threshold_limits.(Vabs, nodes_vmax, nodes_vmin)
    Vnodes = Vabs .* exp.(im * angle.(Vnodes))

    # Create Index set to access different portion of Vnodes array
    Vnodes_dim          = [length(slack_Vm), length(PV_Idx), length(PQ_Idx)]
    Vnodes_offset       = create_offsets(Vnodes_dim)
    Vnodes_Idx          = create_idxs(Vnodes_offset, Vnodes_dim)
    
    Slack_Idx_in_Vnodes = Vnodes_Idx[1]
    PV_V_Idx_in_Vnodes  = Vnodes_Idx[2]
    PQ_V_Idx_in_Vnodes  = Vnodes_Idx[3]

    PV_Vm .= gen_Vm  # abs.(Vnodes[PV_V_Idx_in_Vnodes]) 
    PQ_Vm .= abs.(Vnodes[PQ_V_Idx_in_Vnodes])         

    ##
    # gI =  Ybus * Vnodes  .+  inv(Diagonal( conj.(Vnodes) )) *  conj.(Sd - Cg * Sg)
    gI =  perm_matrix' * Ybus * perm_matrix * Vnodes  .+ Diagonal( inv.(conj.(Vnodes)) ) * ( perm_matrix' * conj.(Sd - Cg * Sg))

    gM     = real.(gI)
    gN     = imag.(gI)
   
    dgM   .= gM[PV_PQ_Idx]
    dgN   .= gN[PV_PQ_Idx]

    return nothing
end

# -----------------------------------------------

function AC_nodal_current_balance(V, p)

    nodes_vmax, nodes_vmin, Cg, Sg, Sd, Ybus, gen_Vm, slack_Vm, slack_Vθ, PV_Idx, PQ_Idx, perm_matrix = p
    
    N_PV_θ      = length(PV_Idx)
    N_PQ_θ      = length(PQ_Idx)
    N_PV_Vm     = length(PV_Idx)
    N_PQ_Vm     = length(PQ_Idx)

    PV_PQ_Idx  = first(PV_Idx):last(PQ_Idx)  # [PV_Idx;PQ_Idx]

    vars_dim    = [N_PV_θ, N_PQ_θ, N_PV_Vm, N_PQ_Vm]
    vars_offset = create_offsets(vars_dim)
    vars_Idx    = create_idxs(vars_offset, vars_dim)
    
    PV_θ       = @view V[vars_Idx[1]]
    PQ_θ       = @view V[vars_Idx[2]]
    PV_Vm      = @view V[vars_Idx[3]]
    PQ_Vm      = @view V[vars_Idx[4]]

    PV_Vm     .= gen_Vm 

    Vnodes = [[slack_Vm .* exp.(im * slack_Vθ)]...;[PV_Vm .* exp.(im*PV_θ)]... ; [PQ_Vm .* exp.(im * PQ_θ)]...]

    Vabs   = abs.(Vnodes)
    Vabs   = threshold_limits.(Vabs, nodes_vmax, nodes_vmin)
    Vnodes = Vabs .* exp.(im * angle.(Vnodes))

    # Create Index set to access different portion of Vnodes array
    Vnodes_dim          = [length(slack_Vm), length(PV_Idx), length(PQ_Idx)]
    Vnodes_offset       = create_offsets(Vnodes_dim)
    Vnodes_Idx          = create_idxs(Vnodes_offset, Vnodes_dim)
    
    Slack_Idx_in_Vnodes = Vnodes_Idx[1]
    PV_V_Idx_in_Vnodes  = Vnodes_Idx[2]
    PQ_V_Idx_in_Vnodes  = Vnodes_Idx[3]

    PV_Vm .= gen_Vm  # abs.(Vnodes[PV_V_Idx_in_Vnodes]) 
    PQ_Vm .= abs.(Vnodes[PQ_V_Idx_in_Vnodes])  

    ##
    # gI   = Ybus * Vnodes .+ inv(Diagonal( conj.(Vnodes) )) *  conj.(Sd - Cg * Sg)
    # gI =  Ybus * perm_matrix * Vnodes  .+ inv( perm_matrix * Diagonal(conj.(Vnodes)) * perm_matrix' ) *  conj.(Sd - Cg * Sg)
    # gI =    Diagonal(Ybus * perm_matrix * Vnodes) + ( perm_matrix * Diagonal(conj.(Vnodes)) * perm_matrix' )\(conj.(Sd - Cg * Sg))

    gI =  perm_matrix' * Ybus * perm_matrix * Vnodes  .+ Diagonal( inv.(conj.(Vnodes)) ) * ( perm_matrix' * conj.(Sd - Cg * Sg))
    
    gM   = real.(gI)
    gN   = imag.(gI)
    
    [gM[PV_PQ_Idx]; gN[PV_PQ_Idx]]    
end

#---------------------------------------------------

####################################################
# Power Flow Formulations
# PowerFlowFormulations.jl
####################################################

function DC_nodal_power_balance_power_flow(nodes, branches, shunts)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    Cg            = get_Cg(nodes)
    Sg            = get_gens_Sg(nodes)
    Pg            = real.(Sg)
    Sd            = get_loads_Sd(nodes)
    Pd            = real.(Sd)
    Pbus_θshift   = get_Pbus_θshift(nodes, branches)
    Gsh           = get_nodal_Gsh(shunts, nodes)
    Bbus          = get_Bbus(nodes, branches)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_name    =  get_nodes_Bus_name(nodes)
    
    Slack_Idx     = idx_slack_nodes
    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)

    Idx_PV_PQ     = [PV_Idx;PQ_Idx]
    nodes_Idx     = [Slack_Idx;PV_Idx;PQ_Idx]

    Θ0            = zeros(length(nodes))

    V0            = Θ0[Idx_PV_PQ] 

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (Cg=Cg, Pg=Pg, Pd=Pd, Gsh=Gsh, Bbus=Bbus, Pbus_θshift=Pbus_θshift, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)

    # https://github.com/SciML/LinearSolve.jl
    # https://docs.sciml.ai/LinearSolve/stable/
    
    A, b       = DC_nodal_power_balance(V0, p)
    prob       = LinearProblem(A, b)

    linsolve   = init(prob)
    sol        = LinearSolve.solve(linsolve)

    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)

    x_dim      = [N_PV_θ, N_PQ_θ]
    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]

    θ          = [slack_Vθ...; PV_θ...;PQ_θ...]
    θ_perm     = perm_matrix * θ

    return Dict("θbus" => θ_perm)
end

#---------------------------------------------------------------------------

########################################################
# Power Flow Data Structures
# PowerFlowDataStructures.jl
########################################################

function Slack!(dx, x, p_agg, t)
    x_c, f_t, y_t, p = p_agg
    
    P      = p[1]
    Q      = p[2]
    Vm     = p[3]
    Vθ     = p[4]

    nothing
end

@kwdef struct Slack
    Bus::String
    name::String     = lowercase(Bus)
    Bus_num::Int64   = parse(Int, split(lowercase(Bus),"bus")[2] )
    kV::Float64      = 1.0
    P::Float64       = 0.0 
    Q::Float64       = 0.0
    S::ComplexF64    = P + im * Q
    Vm::Float64      = 1.0
    Vθ::Float64      = 0.0
    vmax::Float64    = 1.06
    vmin::Float64    = 1.0    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0     
    Bus_type::Symbol = :Slack 

    algebraic_vars_syms::Vector{Symbol} = Symbol[]
    algebraic_vars::Vector{Symbol} = Symbol[]
    state_vars::Vector{Symbol} = Symbol[]
    state_vars_syms::Vector{Symbol} = Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} = Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars), length(algebraic_vars))
    param::Vector{Symbol} = Symbol[:kV, :P, :Q, :Vm, :Vθ, :vmax, :vmim, :Pmax, :Pmin,:Qmax, :Qmin, :Bus_type]
    func::Vector{Function} = Function[Slack!]
end

function Generator!(dx, x, p_agg, t)
    x_c, f_t, y_t, p = p_agg

    P      = p[1]
    Q      = p[2]
    Vm     = p[3]
    Vθ     = p[4]    
    
end

@kwdef struct Generator
    Bus::String
    name::String = lowercase(Bus)
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )
    kV::Float64
    P::Float64
    Q::Float64       = 0.0 # Vector{Float64} = Float64[0.0]
    S::ComplexF64    = P + im * Q
    Vm::Float64      = 1.0
    Vθ::Float64      = 0.0 #  Vector{Float64} = Float64[]
    vmax::Float64    = 1.06
    vmin::Float64    = 0.97    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0     
    Bus_type::Symbol = :Generator

    algebraic_vars_syms::Vector{Symbol} = Symbol[]
    algebraic_vars::Vector{Symbol} = Symbol[]
    state_vars::Vector{Symbol} = Symbol[]
    state_vars_syms::Vector{Symbol} = Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} = Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars), length(algebraic_vars))
    param::Vector{Symbol} = Symbol[:kV, :P, :Q, :Vm, :Vθ, :vmax, :vmim, :Pmax, :Pmin,:Qmax, :Qmin, :Bus_type]
    func::Vector{Function} = Function[Generator!]
end

function Load!(dx, x, p_agg, t)
    x_c, f_t, y_t, p = p_agg

    P      = p[1]
    Q      = p[2]
    Vm     = p[3]
    Vθ     = p[4]
    
    nothing
end

@kwdef struct Load
    Bus::String
    name::String = lowercase(Bus)
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )
    kV::Float64
    P::Float64       = 0.0 
    Q::Float64       = 0.0 
    S::ComplexF64    = P + im * Q
    Vm::Float64      = 1.0 # Vector{Float64} = Float64[]
    Vθ::Float64      = 0.0 # Vector{Float64} = Float64[]
    vmax::Float64    = 1.06
    vmin::Float64    = 0.9    
    Bus_type::Symbol = :Load 

    algebraic_vars_syms::Vector{Symbol} = Symbol[]
    algebraic_vars::Vector{Symbol} = Symbol[]
    state_vars::Vector{Symbol} = Symbol[]
    state_vars_syms::Vector{Symbol} = Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} = Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars), length(algebraic_vars))
    param::Vector{Symbol} = Symbol[:kV, :P, :Q, :Vm, :Vθ, :vmax, :vmim, :Bus_type]
    func::Vector{Function} = Function[Generator!]
end

function ShuntElement!(dx, x, p_agg, t)
    x_c, f_t, y_t, p = p_agg
    
    y_shunt  = p[1]

    nothing
end

@kwdef struct ShuntElement
    Bus::String
    name::String = lowercase(Bus)
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2]) 
    y_shunt::ComplexF64 = 0.0 + im * 0.0  
    Vm::Float64      = 1.0
    Vθ::Float64      = 0.0
    Bus_type::Symbol = :ShuntElement
    algebraic_vars_syms::Vector{Symbol} = Symbol[]
    algebraic_vars::Vector{Symbol} = Symbol[]
    state_vars::Vector{Symbol} = Symbol[]
    state_vars_syms::Vector{Symbol} = Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} = Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars), length(algebraic_vars))
    param::Vector{Symbol} = Symbol[:y_shunt]
    func::Vector{Function} = Function[ShuntElement!]
end

function EnergyStorage!(dx, x, p_agg, t)
    x_c, f_t, y_t, p = p_agg
    
    P      = p[1]
    Q      = p[2]
    Vm     = p[3]
    Vθ     = p[4]

    nothing
end

@kwdef struct EnergyStorage
    Bus::String
    name::String = lowercase(Bus)
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )     
    P::Float64        = 0.0 
    Q::Float64        = 0.0 
    S::ComplexF64     = P + im * Q
    Vm::Float64       = 1.0 
    Vθ::Float64       = 0.0 
    Bus_type::Symbol  = :EnergyStorage

    algebraic_vars_syms::Vector{Symbol} = Symbol[]
    algebraic_vars::Vector{Symbol} = Symbol[]
    state_vars::Vector{Symbol} = Symbol[]
    state_vars_syms::Vector{Symbol} = Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} = Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars), length(algebraic_vars))
    param::Vector{Symbol} = Symbol[:P, :Q, :Vm, :Vθ, :Bus_type]
    func::Vector{Function} = Function[EnergyStorage!]
end

function Line!(dx, x, p_agg, t)
    x_c, f_t, y_t, p = p_agg
    
    from        = p[1]
    to          = p[2]
    y           = p[3]
    y_shunt_km  = p[4]
    y_shunt_mk  = p[5]
    y_ratio     = p[6]

    nothing
end

@kwdef struct Line
    from::String
    to::String
    name::String =  lowercase(from) * "-" * lowercase(to)
    orientation::Tuple{Int64, Int64} = (parse(Int, split(lowercase(from),"bus")[2] ), parse(Int, split(lowercase(to),"bus")[2]) ) 
    y::ComplexF64
    y_shunt_km::Float64 = 0.0
    y_shunt_mk::Float64 = 0.0
    y_ratio::Float64    = 1.0   # = Complex(1.0, 0.0)

    Ybr::Matrix{ComplexF64} = [(y + 1im*y_shunt_km)*1/(abs(y_ratio))^2  -y*1/conj(y_ratio); -y*1/y_ratio y + 1im*y_shunt_mk]
    xs::Float64 = -1*imag(y)/(abs(y))^2
    bi::Float64 = 1/(xs * abs(y_ratio))
    YbrDC::Matrix{ComplexF64} = 1/(im * xs) * [1/(abs(y_ratio))^2  -1/conj(y_ratio); -1/y_ratio 1]
    BbrDC::Matrix{Float64} = bi * [1 -1; -1 1]
    Pshift::Vector{Float64} = angle(y_ratio) * bi * [-1, 1]
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[]
    algebraic_vars::Vector{Symbol} = Symbol[]
    state_vars::Vector{Symbol} = Symbol[]
    state_vars_syms::Vector{Symbol} = Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} = Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars), length(algebraic_vars))
    param::Vector{Symbol} = Symbol[:from, :to, :y, :y_shunt_km, :y_shunt_mk, :y_ratio, :Ybr, :YbrDC]
    func::Vector{Function} = Function[Line!]
end

function Transf!(dx, x, p_agg, t)
    x_c, f_t, y_t, p = p_agg

    from        = p[1]
    to          = p[2]
    y           = p[3]
    y_shunt_km  = p[4]
    y_shunt_mk  = p[5]
    y_ratio     = p[6]
    
    nothing
end

@kwdef struct Transf
    from::String
    to::String
    name::String =  lowercase(from) * "-" * lowercase(to)
    orientation::Tuple{Int64, Int64} = (parse(Int, split(lowercase(from),"bus")[2] ), parse(Int, split(lowercase(to),"bus")[2]) ) 
    y::ComplexF64
    y_shunt_km::Float64 = 0.0
    y_shunt_mk::Float64 = 0.0
    y_ratio::Float64    = 1.0    # ::ComplexF64 = Complex(1.0, 0.0)
    Ybr::Matrix{ComplexF64} = [(y + 1im*y_shunt_km)*1/(abs(y_ratio))^2  -y*1/conj(y_ratio); -y*1/y_ratio y + 1im*y_shunt_mk]
    xs::Float64 = -1*imag(y)/(abs(y))^2
    bi::Float64 = 1/(xs * abs(y_ratio))
    YbrDC::Matrix{ComplexF64} = 1/(im * xs) * [1/(abs(y_ratio))^2  -1/conj(y_ratio); -1/y_ratio 1]
    BbrDC::Matrix{Float64} = bi * [1 -1; -1 1]
    Pshift::Vector{Float64} = angle(y_ratio) * bi * [-1, 1]
        
    algebraic_vars_syms::Vector{Symbol} = Symbol[]
    algebraic_vars::Vector{Symbol} = Symbol[]
    state_vars::Vector{Symbol} = Symbol[]
    state_vars_syms::Vector{Symbol} = Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} = Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars), length(algebraic_vars))
    param::Vector{Symbol} = Symbol[:from, :to, :y, :y_shunt_km, :y_shunt_mk, :y_ratio, :Ybr, :YbrDC]
    func::Vector{Function} = Function[Transf!]
end






"""

function dV_AC_nodal_power_balance_power_flow(nodes, branches, shunts)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    Cg            = get_Cg(nodes)
    Cg_p          = perm_matrix * Cg # perm_matrix * Cg * perm_matrix
    Sg            = get_gens_Sg(nodes)
    Sd            = get_loads_Sd(nodes)
    Sd_p          = perm_matrix * Sd
    
    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0; PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax = nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)    
    
    # sol        = NLSolvers.solve(prob)
    # sol        = DifferentialEquations.solve(prob)
    
    prob       = NonlinearProblem((dx,x,p)-> AC_nodal_power_balance(dx,x,p), x0 , p)
    sol        = NonlinearSolve.solve(prob, NLSolveJL(), abstol=1e-10, reltol=1e-10)   
    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim      = [N_PV_θ, N_PQ_θ, N_PQ_Vm]
    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PQ_Vm      = x[x_Idx[3]]

    V_PQ       = PQ_Vm  .* exp.(im * PQ_θ)
    V_PV       = gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]
    
    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)
    gS_pflow   = gS(V_perm, pgS)
    SBUS       = SbusC(V_perm, Ybus)
    GenSinj    = SBUS + Sd   

    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                  ))    
    
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "Sbus" => SBUS, "gS_pflow" => gS_pflow, "GenSinj" => GenSinj, "dict_init" => dict_init)

end

function dV_AC_nodal_power_balance_power_flow(nodes, branches, shunts, gen_nodes, load_nodes)
    
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    net_size      = get_network_nodal_size(branches)

    Cg            = get_Cg(gen_nodes, net_size)
    Cg_p          = perm_matrix * Cg # perm_matrix * Cg * perm_matrix 
    Sg            = get_gens_Sg(gen_nodes, net_size)

    Cd            = get_Cd(load_nodes, net_size)
    Sd            = get_loads_Sd(load_nodes, net_size)
    Sd_p          = perm_matrix * Sd

    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(gen_nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0; PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax = nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)    
    
    # sol        = NLSolvers.solve(prob)
    # sol        = DifferentialEquations.solve(prob)
    
    prob       = NonlinearProblem((dx,x,p) -> AC_nodal_power_balance(dx,x,p), x0 , p)
    sol        = NonlinearSolve.solve(prob, NLSolveJL(), abstol=1e-10,reltol=1e-10)   
    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim      = [N_PV_θ, N_PQ_θ, N_PQ_Vm]
    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PQ_Vm      = x[x_Idx[3]]

    V_PQ       = PQ_Vm  .* exp.(im * PQ_θ)
    V_PV       = gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]
    
    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Sf         = Diagonal(get_Cf(nodes, branches) * V_perm) * conj.(If)
    St         = Diagonal(get_Ct(nodes, branches) * V_perm) * conj.(It)
    
    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)
    
    SBUS       = SbusC(V_perm, Ybus)
   
    GenSinj    = SBUS + Sd 
    
    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                  ))  
       
    
    return Dict("Sf" => Sf, "St" => St,"Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "gS_pflow" => gS_pflow,  "Sbus" => SBUS, "GenSinj" => GenSinj, "dict_init" => dict_init)

end


function dV_nodal_current_balance_power_flow(nodes, branches, shunts)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    Cg            = get_Cg(nodes)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(nodes)
    Sd            = get_loads_Sd(nodes)
    Sd_p          = perm_matrix * Sd
    
    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0;gen_Vm;PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax=nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)


    prob       = NonlinearProblem((dx,x,p) -> AC_nodal_current_balance(dx,x,p), x0 , p)
    sol        = NonlinearSolve.solve(prob, NLSolveJL(), abstol=1e-10,reltol=1e-10) 
    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PV_Vm    = length(PV_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim    = [N_PV_θ, N_PQ_θ,  N_PV_Vm, N_PQ_Vm]
    x_offset = create_offsets(x_dim)
    x_Idx    = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PV_Vm      = x[x_Idx[3]]
    PQ_Vm      = x[x_Idx[4]]
    
    V_PQ       = PQ_Vm   .* exp.(im * PQ_θ)
    V_PV       = PV_Vm   .* exp.(im * PV_θ)   # gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm  .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]
    
    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)

    SBUS       = SbusC(V_perm, Ybus)
   
    GenSinj    = SBUS + Sd 
    
    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                      ))
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "gS_pflow" => gS_pflow,  "Sbus" => SBUS, "GenSinj" => GenSinj, "dict_init" => dict_init)
    
end


function dV_nodal_current_balance_power_flow(nodes, branches, shunts, gen_nodes, load_nodes)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    net_size      = get_network_nodal_size(branches)

    Cg            = get_Cg(gen_nodes, net_size)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(gen_nodes, net_size)

    Cd            = get_Cd(load_nodes, net_size)
    Sd            = get_loads_Sd(load_nodes, net_size)
    Sd_p          = perm_matrix * Sd 

    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix

    gen_Vm        = get_gens_bus_Vm(gen_nodes)
        
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0;gen_Vm;PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax=nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)


    prob       = NonlinearProblem((dx,x,p) -> AC_nodal_current_balance(dx,x,p), x0 , p)
    sol        = NonlinearSolve.solve(prob, NLSolveJL(), abstol=1e-10,reltol=1e-10) 
    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PV_Vm    = length(PV_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim    = [N_PV_θ, N_PQ_θ,  N_PV_Vm, N_PQ_Vm]
    x_offset = create_offsets(x_dim)
    x_Idx    = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PV_Vm      = x[x_Idx[3]]
    PQ_Vm      = x[x_Idx[4]]
    
    V_PQ       = PQ_Vm   .* exp.(im * PQ_θ)
    V_PV       = PV_Vm   .* exp.(im * PV_θ)   # gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm  .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]
    
    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)

    SBUS       = SbusC(V_perm, Ybus)
    GenSinj    = SBUS + Sd     

    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                  ))         
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "gS_pflow" => gS_pflow,  "Sbus" => SBUS, "GenSinj" => GenSinj, "dict_init" => dict_init)

end


function dV_NLSolve_nodal_power_balance_power_flow(nodes, branches, shunts)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)
    
    Cg            = get_Cg(nodes)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(nodes)
    Sd            = get_loads_Sd(nodes)
    Sd_p          = perm_matrix * Sd
    
    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = perm_matrix * Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0; PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax = nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)   
    
    sol        = nlsolve((dx,x) -> AC_nodal_power_balance(dx,x, p), x0)

    x          = sol.zero

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim      = [N_PV_θ, N_PQ_θ, N_PQ_Vm]
    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PQ_Vm      = x[x_Idx[3]]

    V_PQ       = PQ_Vm  .* exp.(im * PQ_θ)
    V_PV       = gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]
    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)

    SBUS       = SbusC(V_perm, Ybus)
    GenSinj    = SBUS + Sd 

    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                  ))     
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "Sbus" => SBUS, "gS_pflow" => gS_pflow, "GenSinj" => GenSinj, "dict_init" => dict_init)    
    
end


function dV_NLSolve_nodal_current_balance_power_flow(nodes, branches, shunts)

    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    Cg            = get_Cg(nodes)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(nodes)
    Sd            = get_loads_Sd(nodes)
    Sd_p          = perm_matrix * Sd
    
    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = perm_matrix * Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0;gen_Vm;PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax=nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)
    
    sol        = nlsolve((dx,x) -> AC_nodal_current_balance(dx,x, p), x0)
    x          = sol.zero

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PV_Vm    = length(PV_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    #
    x_dim    = [N_PV_θ, N_PQ_θ,  N_PV_Vm, N_PQ_Vm]
    x_offset = create_offsets(x_dim)
    x_Idx    = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PV_Vm      = x[x_Idx[3]]
    PQ_Vm      = x[x_Idx[4]]
    
    V_PQ       = PQ_Vm   .* exp.(im * PQ_θ)
    V_PV       = PV_Vm   .* exp.(im * PV_θ) 
    V_slack    = slack_Vm  .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]
    
    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)

    
    SBUS       = SbusC(V_perm, Ybus)
    GenSinj    = SBUS + Sd 

    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                  ))     
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "Sbus" => SBUS, "gS_pflow" => gS_pflow, "GenSinj" => GenSinj, "dict_init" => dict_init)
    
end    


#---------------------------------------------------------------------------

function AC_nodal_power_balance_power_flow(nodes, branches, shunts)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    Cg            = get_Cg(nodes)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(nodes)
    Sd            = get_loads_Sd(nodes)
    Sd_p          = perm_matrix * Sd
    
    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    = get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0; PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax = nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)    
    
    # sol        = NLSolvers.solve(prob)
    # sol        = DifferentialEquations.solve(prob)
    
    prob       = NonlinearProblem((x,p) -> AC_nodal_power_balance(x,p), x0 , p)
    sol        = NonlinearSolve.solve(prob, NLSolveJL(), abstol=1e-10,reltol=1e-10)   
    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim      = [N_PV_θ, N_PQ_θ, N_PQ_Vm]
    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PQ_Vm      = x[x_Idx[3]]

    V_PQ       = PQ_Vm  .* exp.(im * PQ_θ)
    V_PV       = gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]

    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)
    gS_pflow   = gS(V_perm, pgS)
    
    SBUS       = SbusC(V_perm, Ybus)
    GenSinj    = SBUS + Sd
    
    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                  ))     
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "Sbus" => SBUS, "gS_pflow" => gS_pflow, "GenSinj" => GenSinj, "dict_init" => dict_init)

end


function AC_nodal_power_balance_power_flow(nodes, branches, shunts, gen_nodes, load_nodes)
    
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    net_size      = get_network_nodal_size(branches)

    Cg            = get_Cg(gen_nodes, net_size)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(gen_nodes, net_size)

    Cd            = get_Cd(load_nodes, net_size)
    Sd            = get_loads_Sd(load_nodes, net_size)
    Sd_p          = perm_matrix * Sd 

    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(gen_nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    = get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0; PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax = nodes_vmax, nodes_vmin=nodes_vmin, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)    
    
    # sol        = NLSolvers.solve(prob)
    # sol        = DifferentialEquations.solve(prob)
    
    prob       = NonlinearProblem((x,p) -> AC_nodal_power_balance(x,p), x0 , p)
    sol        = NonlinearSolve.solve(prob, NLSolveJL(), abstol=1e-10,reltol=1e-10)   
    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim      = [N_PV_θ, N_PQ_θ, N_PQ_Vm]
    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PQ_Vm      = x[x_Idx[3]]

    V_PQ       = PQ_Vm  .* exp.(im * PQ_θ)
    V_PV       = gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]
    
    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Sf         = Diagonal(get_Cf(nodes, branches) * V_perm) * conj.(If)
    St         = Diagonal(get_Ct(nodes, branches) * V_perm) * conj.(It)    
    
    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)
    
    SBUS       = SbusC(V_perm, Ybus)
   
    GenSinj    = SBUS + Sd 
    
    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                      ))
    
    return Dict("Sf" => Sf, "St" => St, "Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "gS_pflow" => gS_pflow,  "Sbus" => SBUS, "GenSinj" => GenSinj, "dict_init" => dict_init)

end

"""

"""
Thsi allows changing demand or source power inhection at specific nodes without
changing the structure of the network.
"""


"""
function flexible_nodal_power_balance_power_flow(nodes, branches, shunts; Sd_adaptive=nothing, Sg_adaptive=nothing )
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    Cg            = get_Cg(nodes)
    Cg_p          = perm_matrix * Cg
    
    if Sg_adaptive == nothing
        Sg = get_loads_Sg(nodes)
    else
        Sg = Sg_adaptive
    end
    

    if Sd_adaptive == nothing
        Sd = get_loads_Sd(nodes)
    else
        Sd = Sd_adaptive
    end
    
    Sd_p          = perm_matrix * Sd
    
    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Θ0            = zeros(length(nodes))
    Vm0           = ones(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0; PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax=nodes_vmax_p, nodes_vmin=nodes_vmin, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)

    # sol        = DifferentialEquations.solve(prob)

    prob       = NonlinearProblem((x,p) -> AC_nodal_power_balance(x,p), x0 , p)
    sol        = NonlinearSolve.solve(prob, NLSolveJL())
    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim      = [N_PV_θ, N_PQ_θ, N_PQ_Vm]
    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PQ_Vm      = x[x_Idx[3]]

    V_PQ       = PQ_Vm  .* exp.(im * PQ_θ)
    V_PV       = gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]

    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Sf         = Diagonal(get_Cf(nodes, branches) * V_perm) * conj.(If)
    St         = Diagonal(get_Ct(nodes, branches) * V_perm) * conj.(It)

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)

    SBUS       = SbusC(V_perm, Ybus)
   
    GenSinj    = SBUS + Sd 
    
    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                  ))        
    
    return Dict("Sf" => Sf, "St" => St, "Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "gS_pflow" => gS_pflow,  "Sbus" => SBUS, "GenSinj" => GenSinj, "dict_init" => dict_init)
    
end

function nodal_current_balance_power_flow(nodes, branches, shunts)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)


    Cg            = get_Cg(nodes)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(nodes)
    Sd            = get_loads_Sd(nodes)
    Sd_p          = perm_matrix * Sd
    
    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0;gen_Vm;PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax=nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)


    prob       = NonlinearProblem((x,p) -> AC_nodal_current_balance(x,p), x0 , p)
    sol        = NonlinearSolve.solve(prob, NLSolveJL(), abstol=1e-10,reltol=1e-10) 
    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PV_Vm    = length(PV_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim    = [N_PV_θ, N_PQ_θ,  N_PV_Vm, N_PQ_Vm]
    x_offset = create_offsets(x_dim)
    x_Idx    = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PV_Vm      = x[x_Idx[3]]
    PQ_Vm      = x[x_Idx[4]]
    
    V_PQ       = PQ_Vm   .* exp.(im * PQ_θ)
    V_PV       = PV_Vm   .* exp.(im * PV_θ)   # gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm  .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]

    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)

    SBUS       = SbusC(V_perm, Ybus)
   
    GenSinj    = SBUS + Sd 
    
    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                  ))      
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "gS_pflow" => gS_pflow,  "Sbus" => SBUS, "GenSinj" => GenSinj, "dict_init" => dict_init)    
    
end

function nodal_current_balance_power_flow(nodes, branches, shunts, gen_nodes, load_nodes)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)

    net_size      = get_network_nodal_size(branches)

    Cg            = get_Cg(gen_nodes, net_size)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(gen_nodes, net_size)

    Cd            = get_Cd(load_nodes, net_size)
    Sd            = get_loads_Sd(load_nodes, net_size)
    Sd_p          = perm_matrix * Sd 

    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix

    gen_Vm        = get_gens_bus_Vm(gen_nodes)
        
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0;gen_Vm;PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax=nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)


    prob       = NonlinearProblem((x,p) -> AC_nodal_current_balance(x,p), x0 , p)
    sol        = NonlinearSolve.solve(prob, NLSolveJL(), abstol=1e-10,reltol=1e-10) 
    x          = sol.u

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PV_Vm    = length(PV_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim    = [N_PV_θ, N_PQ_θ,  N_PV_Vm, N_PQ_Vm]
    x_offset = create_offsets(x_dim)
    x_Idx    = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PV_Vm      = x[x_Idx[3]]
    PQ_Vm      = x[x_Idx[4]]
    
    V_PQ       = PQ_Vm   .* exp.(im * PQ_θ)
    V_PV       = PV_Vm   .* exp.(im * PV_θ)   # gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm  .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]

    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)

    SBUS       = SbusC(V_perm, Ybus)
    GenSinj    = SBUS + Sd     

    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                  ))       
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "gS_pflow" => gS_pflow,  "Sbus" => SBUS, "GenSinj" => GenSinj, "dict_init" => dict_init)

end

function NLSolve_nodal_power_balance_power_flow(nodes, branches, shunts)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)
    
    Cg            = get_Cg(nodes)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(nodes)
    Sd            = get_loads_Sd(nodes)
    Sd_p          = perm_matrix * Sd
    
    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0; PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax = nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)   
    
    sol        = nlsolve((x) -> AC_nodal_power_balance(x, p), x0)

    x          = sol.zero

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    x_dim      = [N_PV_θ, N_PQ_θ, N_PQ_Vm]
    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PQ_Vm      = x[x_Idx[3]]

    V_PQ       = PQ_Vm  .* exp.(im * PQ_θ)
    V_PV       = gen_Vm .* exp.(im * PV_θ)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]

    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)

    SBUS       = SbusC(V_perm, Ybus)
   
    GenSinj    = SBUS + Sd 
    
    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                      ))
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "gS_pflow" => gS_pflow,  "Sbus" => SBUS, "GenSinj" => GenSinj, "dict_init" => dict_init)     
    
end


function NLSolve_nodal_current_balance_power_flow(nodes, branches, shunts)
    # Get permutation matrix and associated indices
    perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(nodes)
    
    Cg            = get_Cg(nodes)
    Cg_p          = perm_matrix * Cg
    Sg            = get_gens_Sg(nodes)
    Sd            = get_loads_Sd(nodes)
    Sd_p          = perm_matrix * Sd
    
    Ybus          = get_Ybus(nodes, branches, shunts=shunts)
    Ybus_p        = Ybus *  perm_matrix
    
    gen_Vm        = get_gens_bus_Vm(nodes)
    slack_Vm      = get_slack_nodes_Vm(nodes)
    slack_Vθ      = get_slack_nodes_Vθ(nodes)
    
    nodes_vmax    = get_nodes_vmax(nodes)
    nodes_vmax_p  = perm_matrix' * nodes_vmax
    nodes_vmin    = get_nodes_vmin(nodes)
    nodes_vmin_p  = perm_matrix' * nodes_vmin
    
    nodes_name    =  get_nodes_Bus_name(nodes)

    PV_Idx        = idx_source_nodes # get_gens_bus_Idx(nodes)
    PQ_Idx        = idx_demand_nodes # get_loads_bus_Idx(nodes)
    PV_PQ_Idx     = [PV_Idx;PQ_Idx]

    Vm0           = ones(length(nodes))
    Θ0            = zeros(length(nodes))

    PV_θ0_PQ_θ0   = Θ0[PV_PQ_Idx]
    PQ_Vm0        = Vm0[PQ_Idx]

    x0            = [PV_θ0_PQ_θ0;gen_Vm;PQ_Vm0]

    ########################################################
    # ------------------------------------------------------
    ########################################################

    p  = (nodes_vmax=nodes_vmax_p, nodes_vmin=nodes_vmin_p, Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus, gen_Vm=gen_Vm, slack_Vm=slack_Vm, slack_Vθ=slack_Vθ, PV_Idx=PV_Idx, PQ_Idx=PQ_Idx, perm_matrix=perm_matrix)
    
    sol        = nlsolve((x) -> AC_nodal_current_balance(x, p), x0)
    x          = sol.zero

    N_PV_θ     = length(PV_Idx)
    N_PQ_θ     = length(PQ_Idx)
    N_PV_Vm    = length(PV_Idx)
    N_PQ_Vm    = length(PQ_Idx)

    #
    x_dim    = [N_PV_θ, N_PQ_θ,  N_PV_Vm, N_PQ_Vm]
    x_offset = create_offsets(x_dim)
    x_Idx    = create_idxs(x_offset, x_dim)

    PV_θ       = x[x_Idx[1]]
    PQ_θ       = x[x_Idx[2]]
    PV_Vm      = x[x_Idx[3]]
    PQ_Vm      = x[x_Idx[4]]
    
    V_PQ       = PQ_Vm   .* exp.(im * PQ_θ)
    V_PV       = PV_Vm   .* exp.(im * PV_θ) 
    V_slack    = slack_Vm  .* exp.(im * slack_Vθ )
    V          = [V_slack...; V_PV...;V_PQ...]
    
    V_perm     = perm_matrix * V

    If         = get_Yf(nodes, branches) * V_perm
    It         = get_Yt(nodes, branches) * V_perm

    Ibranches  = If + It

    pgS        = (Cg=Cg, Sg=Sg, Sd=Sd, Ybus=Ybus)

    gS_pflow   = gS(V_perm, pgS)


    SBUS       = SbusC(V_perm, Ybus)
   
    GenSinj    = SBUS + Sd 

    dict_init = OrderedDict(name => (vh, θh, ph, qh) for (name,vh, θh, ph, qh) in zip(nodes_name,
                                                                                  round.( abs.(V_perm); digits=4),
                                                                                  round.( angle.(V_perm); digits=4),
                                                                                  round.( real.(SBUS); digits=4),
                                                                                  round.( imag.(SBUS); digits=4)
                                                                                      ))
    
    return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V_perm), "Vθ" => angle.(V_perm), "Vbus" => V_perm, "Ibranches" => Ibranches, "Ybus" => Ybus, "gS_pflow" => gS_pflow,  "Sbus" => SBUS, "GenSinj" => GenSinj, "dict_init" => dict_init)  
    
end    

"""


# ########################################################
# # ------------------------------------------------------
# # Testing Permuatation
# # ------------------------------------------------------
# """
# b1 = (1,2)
# b2 = (1,3)
# b3 = (2,4)
# b4 = (3,4)
# b5 = (2,3)

# bb     = [b1,b2,b3,b4,b5]
# weight = [2,4,1,4,1] 

# src_v = first.(bb)
# dst_v = last.(bb)

# num_v = maximum([src_v;dst_v])
# num_e = length(bb)

# Cf = sparse(collect(1:num_e), src_v, ones(Int,num_e), num_e, num_v)
# Ct = sparse(collect(1:num_e), dst_v, ones(Int,num_e), num_e, num_v)

# Cbn = Cf - Ct
# Cnb = Cbn'

# find_src(row) = findall((x)->x==1, row)
# find_dst(row) = findall((x)->x==-1, row)

# # src and dst edges
# src_e = [find_src(Cnb[a_row,:]) for a_row in 1:num_v ]
# dst_e = [find_dst(Cnb[a_row,:]) for a_row in 1:num_v ]

# Qv = [[1,2,0,0],[2,3,4,5],[0,4,2,1],[0,5,1,3]]
# Qm = [1 2 0 0; 2 3 4 5; 0 4 2 1; 0 5 1 3]

# y  = [2,1,4,6]
# x  = [2,3,4,1]
# b  = [2,3,4,1]

# A  = Diagonal([1,1,1,1])

# # z  = [4,4,8,7]

# b_idx_perm    = sortperm(b)
# b_sorted      = b[b_idx_perm]
# b_perm_idx    = Permutation(b_idx_perm)
# b_perm_matrix = Matrix(b_perm_idx)

# y_p   = b_perm_matrix * y
# p_A   = b_perm_matrix * A
# A_p   = A * b_perm_matrix
# p_A_p = b_perm_matrix * A * b_perm_matrix

# A     * b        + y                          # original
# A     * b_sorted + y                          # wrong
# A     * b_sorted + y_p                        # wrong
# A_p   * b_sorted + y_p                        # wrong
# p_A   * b_sorted + y_p                        # wrong
# p_A_p * b_sorted + y_p                        # wrong
# b_perm_matrix  * (p_A_p * b_sorted + y_p)     # wrong
# inv(b_perm_matrix) * (p_A_p * b_sorted + y_p) # correct

# """



# """
    
#     # gI   = perm_matrix * Ybus * Vnodes  .+  perm_matrix *  inv(Diagonal( conj.(Vnodes) )) * conj.(Sd - Cg * Sg)
#     # gI   = perm_matrix * Ybus * Vnodes  .+  inv(Diagonal( conj.(Vnodes) )) *  perm_matrix * conj.(Sd - Cg * Sg)
    
#     # gS = perm_matrix * (Diagonal(Vnodes) * conj.(Ybus) * conj.(Vnodes) + Sd - Cg * Sg)
#     # gS =  (Diagonal(Vnodes) * conj.(perm_matrix * Ybus) * conj.(Vnodes) + Sd - perm_matrix * Cg * Sg)
#     # gI   = perm_matrix * Ybus * Vnodes  .+  perm_matrix *  inv(Diagonal( conj.(Vnodes) )) * conj.(Sd - Cg * Sg)
#     # gI   = perm_matrix * Ybus * Vnodes  .+  inv(Diagonal( conj.(Vnodes) )) *  perm_matrix * conj.(Sd - Cg * Sg)
#     # gS = perm_matrix * (Diagonal(Vnodes) * conj.(Ybus) * conj.(Vnodes) + Sd - Cg * Sg)
#     # gS =  (Diagonal(Vnodes) * conj.(perm_matrix * Ybus) * conj.(Vnodes) + Sd - perm_matrix * Cg * Sg)    

# """
