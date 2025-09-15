# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

####################################################

#---------------------------------------------------
#---------------------------------------------------
# utilities  functions
#---------------------------------------------------
#---------------------------------------------------

using InteractiveUtils

"""
    anonymus_func(x)


Returns a place holder for an anonymus function that produce nothing.


It is needed in struct that do not have some functions implemented.
"""
anonymus_func = (x) -> nothing



"""

Converts an expression to a string.

# Example
a = 1
b = @name a

b

"a"
    """

macro name(arg)
    x = string(arg)
    quote
        $x
    end
end


# x = "[gov_ieee_tgov1_cb__overby_param, gov_ieee_tgov1_cb__sauer_param, gov_ieee_m, gov_ieee_tgov1_cb__a_param, gov_ieee_tgov1_cb__param, gov_ieee_tgov1_cb__1_pacb__3_param, gov_t1_cb__millano_param, gov_t1_cb__1_param, gov_t1_cb_sauer__1_pa_5_param]"

# test1 =
#     Symbol.(String.(strip.(
#         split(split(split(x,"[")[2],"]")[1], ","))))


macro macro_get_list_items_name_and_value(arg)
    x = string(arg)
    x = Symbol.(String.(strip.(split(split(split(x,"[")[2],"]")[1], ","))))

    return x, eval(arg)
    
end

# ------------------------------------------------------

"""
    second(elements_container)


Returns the second element in a list.

"""
function second(elements_container)
    return elements_container[2]
end


"""
    third(elements_container)


Returns the third element in a list.

"""
function third(elements_container)
    return elements_container[3]
end


"""
    fourth(elements_container)


Returns the fourth element in a list.

"""
function fourth(elements_container)
    return elements_container[4]
end


"""
    fifth(elements_container)


Returns the fifth element in a list.

"""
function fifth(elements_container)
    return elements_container[5]
end


"""
    sixth(elements_container)


Returns the sixth element in a list.

"""
function sixth(elements_container)
    return elements_container[6]
end

"""
    seventh(elements_container)


Returns the seventh element in a list.

"""
function seventh(elements_container)
    return elements_container[7]
end

"""
    eighth(elements_container)


Returns the eighth element in a list.

"""
function eighth(elements_container)
    return elements_container[8]
end


"""
    nineth(elements_container)


Returns the nineth element in a list.

"""
function nineth(elements_container)
    return elements_container[9]
end

"""
    tenth(elements_container)


Returns the tenth element in a list.

"""
function tenth(elements_container)
    return elements_container[10]
end

#----------------------------------------

"""
    get_dict_first_to_tenth_funs(
        no_nth_funcs)

Returns a dictionary of n numbers of functions `first`, `second` ... `tenth`.
"""
function get_dict_first_to_tenth_funs(
    no_nth_funcs)

    @assert (no_nth_funcs > 0 && no_nth_funcs < 11)
    
    list_funcs =
        Function[first, second, third, fourth,
                 fifth, sixth, seventh, eighth,
                 nineth, tenth]
    
    return Dict{Int64, Function}(
        idx => list_funcs[idx] for idx in 1:no_nth_funcs )

end


#----------------------------------------

"""
    get_n2s_any(
        a_net_group_idxs;
        nothing_bool= false)


Returns "indices to ordinal" dictionaries of indices of any
type of nodes in a network.
"""
function get_n2s_any(
    a_net_group_idxs;
    nothing_bool= false)

    if nothing_bool == false

        return OrderedDict{Union{Symbol,String,Int64},Int64}(
                net_idx =>idx
                for (net_idx, idx) in zip(
                    a_net_group_idxs,
                    collect(1:length(
                        a_net_group_idxs )) ) )        
    else
        
        return OrderedDict{Union{Symbol,String,Int64},
                           Int64}( )
    end

    
end



"""
    get_a_n2s_dict(a_type_idxs )


Returns "indices to ordinal" dictionaries of indices of a list.

# Example

a_type_idx = Union{Symbol,String,Int64}[1,2,"dayo", :yusuff]

n2s_a_type =  get_a_n2s_dict(a_type_idx )

"""
function get_a_n2s_dict(a_type_idxs )

    return OrderedDict{Union{Symbol,String,Int64},Int64}(
        idx_key =>idx
        for (idx_key, idx) in
            zip(a_type_idxs,
                collect(1:length(a_type_idxs))))

end

#----------------------------------------
# Network graph related
#----------------------------------------

"""
    find_src(row)

Returns source nodes.
"""
find_src(row) = findall((x)->x==1,  row)


"""
    find_dst(row)

Returns destination nodes.
"""
find_dst(row) = findall((x)->x==-1, row)


"""
    find_node_incident_edges(row)

Returns a node incident edges.
"""
find_node_incident_edges(row) =
    findall((x)->(x==1||x==-1), row)



"""
    get_edges_orientation(edges_collection)

Returns edges orientation `(src, dst)`.
"""
function get_edges_orientation(edges_collection)
    
    if isa(edges_collection, Union{Array, Vector})
        return map(get_edge_orientation,
                   collect(edges_collection) )
    else
        return map(get_edge_orientation,
                   collect(values(edges_collection)))
    end  
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



function get_nodes_src_edges(
    edges_collection::Union{OrderedDict, Vector})
    
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



"""
    get_Cnb_by_orientations(edges_collection)

Returns nodes to branches connectivity matrix `Cnb`.
"""
function get_Cnb_by_orientations(
    edges_orientations)
    
    src_v = first.(edges_orientations)
    
    dst_v = last.(edges_orientations)
    
    num_v = maximum([src_v;dst_v])
    
    Cf    = sparse( collect(1:length(src_v)),
                   src_v, ones(Int,length(src_v)),
                   length(src_v), num_v)
    
    Ct    = sparse(collect(1:length(dst_v)),
                   dst_v, ones(Int,length(dst_v)),
                   length(dst_v), num_v)
    Cbn = Cf - Ct

    return Cbn'
    
end

"""
    get_Cbn_by_orientations(edges_collection)

Returns branches to nodes connectivity matrix `Cbn`.
"""
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

"""
    get_nodes_incident_edges_by_orientations(
        edges_orientations )


Returns nodes incident edges.
"""
function get_nodes_incident_edges_by_orientations(
    edges_orientations )
        
    src_v = first.(edges_orientations)
    
    dst_v = last.(edges_orientations)
    
    num_v = maximum([src_v;dst_v])
        
    # num_v = length(collect(keys(nodes)))
    
    Cnb = get_Cnb_by_orientations(
        edges_orientations)
    
    # return [find_node_incident_edges(Cnb[a_row,:])
    #         for a_row in 1:num_v ]

    all_nodes_idx =
        sort(unique([src_v;dst_v]))
    
    return [find_node_incident_edges(Cnb[a_row,:])
            for a_row in all_nodes_idx ]
    
end

"""
    get_node_src_edges(node_number::Int64, Cnb)


Returns a node source edges, they are edges whose the node is a source node. 

"""
function get_node_src_edges(node_number::Int64, Cnb)

    return [find_src(Cnb[node_number,:])
            for a_row in
                node_number:node_number ]
    
end

"""
    get_node_dst_edges((node_number::Int64, Cnb)


Returns a node destination edges, they are edges whose the node is a destination node. 

"""
function get_node_dst_edges(node_number::Int64, Cnb)

    return [find_dst(Cnb[node_number,:])
            for a_row in
                node_number:node_number ]
    
end

#-------------------------------------------------------


#----------------------------------------
#----------------------------------------


"""
    get_nodes_incident_edges_to_or_fro(
        to_or_fro_edges_current )

Returns sum of injection currents at a node by its source or destination edges.
"""
function get_nodes_incident_edges_to_or_fro(
    to_or_fro_edges_current )
    
    return sum.([ ih_or_ik == [] ?
        x_from_xr_xi([0.0, 0.0]) :
        x_from_xr_xi.([ [current_set...;]
                        for current_set in ih_or_ik ])
                  for ih_or_ik in
                      to_or_fro_edges_current ])
end

"""
    edges_current_partial_sum(
        ih_or_ik )

Returns sum of current from source or destination edges connected to a node.
"""
function edges_current_partial_sum(
    ih_or_ik )

    if ih_or_ik != []

        return sum( [ ih == [] ?
            x_from_xr_xi([0.0, 0.0]) :
            x_from_xr_xi(ih)
                      for ih in ih_or_ik ] )
    else
        return 0.0 + im * 0.0
    end    
    
end


"""
    dynamic_nodal_current_balance(
        x_s, x_d)


Returns dynamic nodal current balance at a node.
"""
function dynamic_nodal_current_balance(
    x_s, x_d)
    
    if x_s != [] &&  x_d != []
        
        is_r = sum(first.(x_s))
        is_i = sum(last.(x_s))
        id_r = sum(first.(x_d))
        id_i = sum(last.(x_d))

        return is_r + id_r  + im * (is_i + id_i)
        
    elseif x_s == [] &&  x_d != []
        
        id_r = sum(first.(x_d))
        id_i = sum(last.(x_d))

        return id_r  + im *  id_i
        
    elseif x_s != [] &&  x_d == []
        is_r = sum(first.(x_s))
        is_i = sum(last.(x_s))

        return is_r  + im * is_i
        
    else
        throw("The node edges source and destination currents are zero ")
    end
end


"""
    pf_dynamic_nodal_current_balance(
        x_s, x_d)


Returns dynamic nodal current balance at a node.
"""
function pf_dynamic_nodal_current_balance(
    x_s, x_d)
    
    if x_s != [] &&  x_d != []
        is  = sum(x_s)
        id = sum(x_d)
        # @show "A"
        return is + id
    elseif x_s == [] &&  x_d != []

        # @show "B"
        return sum(x_d)
    elseif x_s != [] &&  x_d == []
        # @show "C"
        return sum(x_s) 
    else
        throw("The node edges source and destination currents are zero ")
    end
end


# ------------------------------------------------------
# dq transform for Z
# ------------------------------------------------------

"""
    Z_dq(ra, X_d_dash, X_q_dash)


Returns Z_dq matrix transform for a generator.
"""
function Z_dq(ra, X_d_dash, X_q_dash)

    return [ ra         -X_q_dash;
             X_d_dash    ra]

end



# function invZ_dq(ra, X_d_dash, X_q_dash)

#     return (1.0 / ( ra * ra + X_d_dash * X_q_dash )) *
#         [ra X_q_dash; -X_d_dash ra]

# end

"""
    invZ_dq(ra, X_d_dash, X_q_dash)


Returns inverse Z_dq matrix transform for a generator.
"""
function invZ_dq(ra, X_d_dash, X_q_dash)

    return inv( Z_dq(ra, X_d_dash, X_q_dash) )

end


"""
    network_to_machine_ref_frame(x_r, x_i, δ)


Returns a network to machine reference frame. Network reference frame to machine reference frame Sauer page 160

"""
function network_to_machine_ref_frame(x_r, x_i, δ)

    A = [  sin(δ)  -cos(δ);
           cos(δ)   sin(δ)]
    
    x_ri = [x_r, x_i]

    x_dq = A \ x_ri

    return (x_d = x_dq[1], x_q = x_dq[2])

end


"""
    machine_to_network_ref_frame(x_d, x_q, δ)


Returns a  machine reference frame to network reference frame, Sauer page 160

"""
function machine_to_network_ref_frame(x_d, x_q, δ)

    A = [  sin(δ)    cos(δ);
          -cos(δ)   sin(δ)]
    
    x_dq = [x_d, x_q]

    x_ri = A \ x_dq

    return (x_r = x_ri[1], x_i = x_ri[2])

end

#------------------------------------------------

"""
    update_vec_of_vec_from_comp_axis!(
        vec_of_vec, P_comp, Q_comp )


Updates `vec_of_vec`, by `P_comp`, `Q_comp`.

# Example 

a = [[0.0, 0.0], [1.0, 0.0], [0.1,0.2]]
t_p = [1.2, 3.0, 4.0]
t_q = [2.2, 3.1, 1.0]

update_vec_of_vec_from_comp_axis!(
    a , t_p , t_q )

"""
function update_vec_of_vec_from_comp_axis!(
    vec_of_vec, P_comp, Q_comp )

    for (idx,(P, Q)) in enumerate(zip(P_comp, Q_comp ))
        vec_of_vec[idx] .= [P, Q]
    end
    
end


"""
    update_padded_vec_of_vec_from_comp_axis!(
        padded_vec_of_vec,
        P_comp, Q_comp,
        idxs )

Updates `padded_vec_of_vec`, by `P_comp`, `Q_comp` based on `idxs`.

# Example

nodes_size = 9

padded_vec_of_vec = [[0.0, 0.0] for idx in 1:nodes_size ]

P_comp = [1.2, 3.0, 4.0]
Q_comp = [2.2, 3.1, 1.0]

idxs = [2,5,9]

"""
function update_padded_vec_of_vec_from_comp_axis!(
    padded_vec_of_vec,
    P_comp, Q_comp,
    idxs )

    for (idx, P, Q) in zip(idxs, P_comp, Q_comp )
        padded_vec_of_vec[idx] .= [P, Q]
    end
    
end


"""
    update_vec_of_vec_from_vec_vec!(
        vec_of_vec, PQ_vec_vec, idxs )


Updates `vec_of_vec`, by `PQ_vec_vec` based on `idxs`.

# Example

a = [[0.0, 0.0], [0.0, 0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0] ]

t_p = [[1.2, 3.0], [4.0, 6.0]]

idxs = [2, 4, 6]

"""
function update_vec_of_vec_from_vec_vec!(
    vec_of_vec, PQ_vec_vec, idxs )

    for (idx, PQ_vec) in zip(idxs, PQ_vec_vec )
        vec_of_vec[idx] .= PQ_vec
    end
    
end


"""    update_vec_of_vec_from_flattened!(
        vec_of_vec, flattened_vec, vec_vec_Idx )


Updates `vec_of_vec`, by `flattened_vec` based on `vec_vec_Idx`.
"""
function update_vec_of_vec_from_flattened!(
    vec_of_vec, flattened_vec, vec_vec_Idx )

    for idx in 1:length(vec_vec_Idx )
        vec_of_vec[idx] .=
            flattened_vec[ vec_vec_Idx[idx] ]
    end
    

end

#------------------------------------------------

"""
    get_flattened_to_components_vector_var_Idx(
        vec_of_vec_var )

Returns indices in a flattened vector for each vector in a vector of vector.
"""
function get_flattened_to_components_vector_var_Idx(
    vec_of_vec_var )


    vec_of_vec_var_dim =
        length.( vec_of_vec_var )
    
    _, _, vec_of_vec_var_Idx =
        create_size_offset_Idx(
            vec_of_vec_var_dim)

    return vec_of_vec_var_Idx
    
end

#------------------------------------------------


"""
    get_gen_nodes_ω_ed_dash_eq_dash_views(
        state, nodes_ω_ed_dash_eq_dash_Idxs )


Returns a flattend view of gens_ω, gens_ed_dash, gens_eq_dash based on nodes_ω_ed_dash_eq_dash_Idxs.
"""
function get_gen_nodes_ω_ed_dash_eq_dash_views(
    state, nodes_ω_ed_dash_eq_dash_Idxs ) 

    return [ view(state, idx )
             for idx in
                 nodes_ω_ed_dash_eq_dash_Idxs  ]
    
end


"""
    get_gen_nodes_ω_ed_dash_eq_dash(
        state, nodes_ω_ed_dash_eq_dash_Idxs )


Returns a flattened vector of gens_ω, gens_ed_dash, gens_eq_dash based on nodes_ω_ed_dash_eq_dash_Idxs.
"""
function get_gen_nodes_ω_ed_dash_eq_dash(
    state, nodes_ω_ed_dash_eq_dash_Idxs )
    
    return [state[idx]
            for idx in
                nodes_ω_ed_dash_eq_dash_Idxs  ]
    
end


"""
    get_gen_nodes_δ_ω_ed_dash_eq_dash_views(
        state, nodes_δ_ω_ed_dash_eq_dash_Idxs )


Returns a flattened view of gens_δ, gens_ω, gens_ed_dash, gens_eq_dash based on nodes_δ_ω_ed_dash_eq_dash_Idxs.
"""
function get_gen_nodes_δ_ω_ed_dash_eq_dash_views(
    state, nodes_δ_ω_ed_dash_eq_dash_Idxs )
    
    return [ view(state, idx)
             for idx in
                 nodes_δ_ω_ed_dash_eq_dash_Idxs ]
    
end


"""
    get_gen_nodes_δ_ω_ed_dash_eq_dash(
        state, nodes_δ_ω_ed_dash_eq_dash_Idxs )


Returns a flattened vector of gens_δ, gens_ω, gens_ed_dash, gens_eq_dash based on nodes_δ_ω_ed_dash_eq_dash_Idxs.
"""
function get_gen_nodes_δ_ω_ed_dash_eq_dash(
    state, nodes_δ_ω_ed_dash_eq_dash_Idxs )
    
    return [state[idx]
            for idx in
                nodes_δ_ω_ed_dash_eq_dash_Idxs ]
    
end

#------------------------------------------------


"""
    get_dyn_red_vh_θh_idq(
        uh,
        δ_ω_ed_dash_eq_dash_view;
        <keywords arguments> )

# Arguments

- gens_vh
- ra_Xd_dash_Xq_dash_view
- red_vh_θh_idx
- n2s_gens_idx
- n2s_non_gens_idx
- gens_nodes_idx
- nodes_size

Returns `[ red_vh_θh; gens_idq_flat ]`.
"""
function get_dyn_red_vh_θh_idq(
    uh,
    δ_ω_ed_dash_eq_dash_view;
    gens_vh =
        gens_vh,
    ra_Xd_dash_Xq_dash_view =
        ra_Xd_dash_Xq_dash_view,
    red_vh_θh_idx = red_vh_θh_idx,
    n2s_gens_idx =
        n2s_gens_idx,
    n2s_non_gens_idx =
        n2s_non_gens_idx,
    gens_nodes_idx =
        gens_nodes_idx,
    nodes_size = nodes_size )

    
    x_vh = [ idx ∈ gens_nodes_idx ?
        gens_vh[ n2s_gens_idx[ idx ] ] : abs(
            uh[ n2s_non_gens_idx[idx] ] )
              for idx in 1:nodes_size ]
    
    x_θh = [ angle( uh[ n2s_non_gens_idx[idx] ] )
             for idx in 1:nodes_size ]

    uh  = x_vh .* exp.(im * x_θh)

    vh_θh = vcat( x_vh, x_θh )
     
    red_vh_θh = vh_θh[ red_vh_θh_idx  ]

    #--------------------------------------------
    
    gens_idq =
        [  get_pf_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( uh[gens_nodes_idx]  ),
                   angle.( uh[gens_nodes_idx] ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    gens_i_d = first.( gens_idq )
        
    gens_i_q = second.( gens_idq )

    gens_idq_flat = [ gens_i_d; gens_i_q ]
        
    # ---------------------------------------------------

    # red_vh_θh_idq = [ red_vh_θh; gens_idq_flat ]
    
    return [ red_vh_θh; gens_idq_flat ]
end

#------------------------------------------------

"""
    get_dynamic_idq_vhθh(
        vh_θh,
        δ_ω_ed_eq,
        ra_X_d_dash_X_q_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_dynamic_idq_vhθh(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )


    return  get_dynamic_idq_vhθh(
        vh_θh..., δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_dynamic_idq_vhθh(
        vh,
        θh,
        δ_ω_ed_eq,
        ra_Xd_dash_Xq_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_dynamic_idq_vhθh(
    vh,
    θh,
    δ_ω_ed_eq,
    ra_Xd_dash_Xq_dash )
        
    return  get_dynamic_idq_vhθh(
        vh, θh, δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash...)  
end




# function get_dynamic_idq_vhθh(
#     vh,
#     θh,
#     δ_ω_ed_eq,
#     ra_X_d_dash_X_q_dash )

#     if δ_ω_ed_eq == []
        
#         return [0.0, 0.0]
        
#     else
        
#         δ, ω, ed_dash, eq_dash = δ_ω_ed_eq
        
#         ra, X_d_dash, X_q_dash = ra_X_d_dash_X_q_dash

#         #id_iq = invZ_dq(ra, X_d_dash, X_q_dash) * [ed_dash - vh * sin(δ - θh), eq_dash - vh * cos(δ - θh)]
        
#         return invZ_dq(ra, X_d_dash, X_q_dash) * [
#             ed_dash - vh * sin(δ - θh), eq_dash -
#                 vh * cos(δ - θh)]
#     end
    
# end


# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_dynamic_idq_vhθh(
        vh,
        θh,
        δ,
        ω,
        ed_dash,
        eq_dash,
        ra,
        X_d_dash,
        X_q_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_dynamic_idq_vhθh(
    vh,
    θh,
    δ,
    ω,
    ed_dash,
    eq_dash,
    ra,
    X_d_dash,
    X_q_dash )

    if δ == []
        return [0.0, 0.0]
    else
        return invZ_dq(ra, X_d_dash, X_q_dash) * [
            ed_dash - vh * sin(δ - θh), eq_dash -
                vh * cos(δ - θh)]
    end
    
end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_dynamic_idq_vhθh(
        vh,
        θh,
        δ,
        ed_dash,
        eq_dash,
        ra,
        X_d_dash,
        X_q_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_dynamic_idq_vhθh(
    vh,
    θh,
    δ,
    ed_dash,
    eq_dash,
    ra,
    X_d_dash,
    X_q_dash )

    return invZ_dq(ra, X_d_dash, X_q_dash) * [
        ed_dash - vh * sin(δ - θh), eq_dash -
            vh * cos(δ - θh)]

end

#------------------------------------------------

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_dynamic_idq_θ_π_vhθh(
        vh_θh,
        δ_ω_ed_eq,
        ra_X_d_dash_X_q_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_dynamic_idq_θ_π_vhθh(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )


    return get_dynamic_idq_θ_π_vhθh(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end


# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_dynamic_idq_θ_π_vhθh(
        vh,
        θh,
        δ_ω_ed_eq,
        ra_X_d_dash_X_q_dash )

Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_dynamic_idq_θ_π_vhθh(
    vh,
    θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )

    if δ_ω_ed_eq == []
        
        return 0.0 + im * 0.0
        
    else

        return get_dynamic_idq_θ_π_vhθh(
            vh_θh...,
            δ_ω_ed_eq...,
            ra_X_d_dash_X_q_dash... )
    end
    
end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_dynamic_idq_θ_π_vhθh(
        vh, θh,
        δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_dynamic_idq_θ_π_vhθh(
    vh, θh,
    δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )
        
    id_iq = invZ_dq(ra, X_d_dash, X_q_dash) * [
        ed_dash - vh * sin(δ - θh), eq_dash -
            vh * cos(δ - θh)]

    return  (id_iq[1] + im * id_iq[2]) * exp(
        im * (δ - pi/2))    
end


# function get_dynamic_idq_θ_π_vhθh(
#     vh, θh,
#     δ, ω, ed_dash, eq_dash,
#     ra, X_d_dash, X_q_dash )

#     if δ_ω_ed_eq == []
        
#         return 0.0 + im * 0.0
        
#     else
        
#         # δ, ω, ed_dash, eq_dash = δ_ω_ed_eq
        
#         # ra, X_d_dash, X_q_dash = ra_X_d_dash_X_q_dash
        
#         id_iq = invZ_dq(ra, X_d_dash, X_q_dash) * [
#             ed_dash - vh * sin(δ - θh), eq_dash -
#                 vh * cos(δ - θh)]

#         return  (id_iq[1] + im * id_iq[2]) * exp(
#             im * (δ - pi/2))
#     end
    
# end

#------------------------------------------------

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_pf_dynamic_idq_θ_π_vhθh(
        vh_θh, δ_ω_ed_eq, ra_X_d_dash_X_q_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_pf_dynamic_idq_θ_π_vhθh(
    vh_θh, δ_ω_ed_eq, ra_X_d_dash_X_q_dash )


    return  get_pf_dynamic_idq_θ_π_vhθh(
        vh_θh..., δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_pf_dynamic_idq_θ_π_vhθh(
        vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_pf_dynamic_idq_θ_π_vhθh(
    vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
        
    return  get_pf_dynamic_idq_θ_π_vhθh(
        vh, θh, δ_ω_ed_eq..., ra_Xd_dash_Xq_dash... )  
end


# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_pf_dynamic_idq_θ_π_vhθh(
        vh, θh,
        δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_pf_dynamic_idq_θ_π_vhθh(
    vh, θh,
    δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )
        
    id_iq = invZ_dq(ra, X_d_dash, X_q_dash) * [
        ed_dash - vh * sin(δ - θh), eq_dash -
            vh * cos(δ - θh)]

    return  (id_iq[1] + im * id_iq[2]) * exp(
        im * (δ - pi/2))    
end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_pf_dyn_idq_θ_π_vhθh(
        vh_θh,
        δ_ω_ed_eq,
        ra_X_d_dash_X_q_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_pf_dyn_idq_θ_π_vhθh(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )

    return  get_pf_dyn_idq_θ_π_vhθh(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end


# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_pf_dyn_idq_θ_π_vhθh(
        vh, θh,
        δ_ω_ed_eq,
        ra_Xd_dash_Xq_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_pf_dyn_idq_θ_π_vhθh(
    vh, θh,
    δ_ω_ed_eq,
    ra_Xd_dash_Xq_dash )
        
    return  get_pf_dyn_idq_θ_π_vhθh(
        vh, θh,
        δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash... )  
end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_pf_dyn_idq_θ_π_vhθh(
        vh, θh,
        δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )


Returns a generator direct current `id` and quadrature current `iq`.
"""
function get_pf_dyn_idq_θ_π_vhθh(
    vh, θh,
    δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    return  invZ_dq(ra, X_d_dash, X_q_dash) *
        [ed_dash - vh * sin(δ - θh),
         eq_dash - vh * cos(δ - θh)] *
             exp(im * (δ - π/2))

end

# -------------------------------------------

# function get_pf_dyn_idq_net( idq, δ )
    
#     return  idq * exp(im * (δ - π/2))

# end

"""
    get_pf_dyn_idq_net( id, iq, δ )


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_pf_dyn_idq_net( id, iq, δ )

    idq = (id + im * iq)  * exp(im * (δ - π/2))
    return  [real(idq ), imag(idq )]

end

# @doc (@doc get_pf_dyn_idq_net)
"""
    get_pf_dyn_idq_net( idq, δ )


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_pf_dyn_idq_net( idq, δ )
    
    return  get_pf_dyn_idq_net( idq..., δ )

end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_pf_dyn_idq(
        vh_θh,
        δ_ω_ed_eq,
        ra_X_d_dash_X_q_dash )   


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_pf_dyn_idq(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )

    return  get_pf_dyn_idq(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_pf_dyn_idq(
        vh, θh,
        δ_ω_ed_eq,
        ra_Xd_dash_Xq_dash )


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_pf_dyn_idq(
    vh, θh,
    δ_ω_ed_eq,
    ra_Xd_dash_Xq_dash )
        
    return  get_pf_dyn_idq(
        vh, θh,
        δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash... )  
end


# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_pf_dyn_idq(
        vh, θh, δ, ω,
        ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_pf_dyn_idq(
    vh, θh, δ, ω,
    ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    return  invZ_dq(ra, X_d_dash, X_q_dash) *
        [ed_dash - vh * sin(δ - θh),
         eq_dash - vh * cos(δ - θh)]

end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_a_gen_dyn_idq(
        vh_θh,
        δ_ω_ed_eq,
        ra_X_d_dash_X_q_dash )

Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_a_gen_dyn_idq(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )

    return  get_a_gen_dyn_idq(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_a_gen_dyn_idq(
        vh, θh,
        δ_ω_ed_eq,
        ra_Xd_dash_Xq_dash )


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_a_gen_dyn_idq(
    vh, θh,
    δ_ω_ed_eq,
    ra_Xd_dash_Xq_dash )
        
    return  get_a_gen_dyn_idq(
        vh, θh,
        δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash... )  
end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_a_gen_dyn_idq(
        vh, θh,
        δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_a_gen_dyn_idq(
    vh, θh,
    δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    return  invZ_dq(ra, X_d_dash, X_q_dash) *
        [ed_dash - vh * sin(δ - θh),
         eq_dash - vh * cos(δ - θh)]

end

# @doc (@doc get_dynamic_idq_vhθh)
"""
    get_a_gen_dyn_idq(
        vh, θh,
        δ, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_a_gen_dyn_idq(
    vh, θh,
    δ, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    return  invZ_dq(ra, X_d_dash, X_q_dash) *
        [ed_dash - vh * sin(δ - θh),
         eq_dash - vh * cos(δ - θh)]

end

# -------------------------------------------


"""
    get_dynamic_idq_ur_ui(
        u_r, u_i, δ_ω_ed_eq, ra_X_d_dash_X_q_dash)


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_dynamic_idq_ur_ui(
    u_r, u_i, δ_ω_ed_eq, ra_X_d_dash_X_q_dash)

    if δ_ω_ed_eq == []
        return [0.0, 0.0]
    else
        
        δ, ed_dash, eq_dash = δ_ω_ed_eq
        
        ra, X_d_dash, X_q_dash = ra_X_d_dash_X_q_dash
        
        return invZ_dq(ra, X_d_dash, X_q_dash) * [
            ed_dash, eq_dash] - invZ_dq(
                ra, X_d_dash, X_q_dash) * [
                    sin(δ) -cos(δ); cos(δ) sin(δ)] * [
                        u_r, u_i]
    end
    
end

# @doc (@doc get_dynamic_idq_ur_ui)
"""
    get_dynamic_idq_ur_ui(
        u_r, u_i, δ, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_dynamic_idq_ur_ui(
    u_r, u_i, δ, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash)

    return invZ_dq(ra, X_d_dash, X_q_dash) * [
        ed_dash, eq_dash] - invZ_dq(
            ra, X_d_dash, X_q_dash) * [
                sin(δ) -cos(δ); cos(δ) sin(δ)] * [
                    u_r, u_i]
    
end

# @doc (@doc get_dynamic_idq_ur_ui)
"""
    get_dynamic_idq_θ_π_ur_ui(
        u_r, u_i, δ_ω_ed_eq,
        ra_X_d_dash_X_q_dash)


Returns a generator direct current `id` and quadrature current `iq` in network reference frame.
"""
function get_dynamic_idq_θ_π_ur_ui(
    u_r, u_i, δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash)

    if δ_ω_ed_eq == []
        
        return 0.0 + im * 0.0
        
    else
        
        δ, ω, ed_dash, eq_dash = δ_ω_ed_eq
        
        ra, X_d_dash, X_q_dash = ra_X_d_dash_X_q_dash
        
        id_iq = invZ_dq(ra, X_d_dash, X_q_dash) *
            [ed_dash, eq_dash] -
            invZ_dq(ra, X_d_dash, X_q_dash) * [
                sin(δ) -cos(δ); cos(δ) sin(δ)] * [
                    u_r, u_i]
        
        return  (id_iq[1] + im * id_iq[2]) *
            exp(im * (δ - pi/2))
    end
    
end


"""
    get_dynamic_τe_from_id_iq(
       id, iq, δ, ed_dash, eq_dash,  ra,
       X_d_dash, X_q_dash)


Returns a generator torque  for a given `id`, `iq`, `δ`, `ed_dash`, `eq_dash`,  `ra`, `X_d_dash`, `X_q_dash`.
"""
function get_dynamic_τe_from_id_iq(
    id, iq, δ, ed_dash, eq_dash,  ra,
    X_d_dash, X_q_dash)

    return ed_dash * id + eq_dash * iq +
        ( X_q_dash - X_d_dash ) * id * iq
    
end

"""
    get_dynamic_τe_from_id_iq(
        id, iq, δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)


Returns a generator torque  for a given `id`, `iq`, `δ`, `ed_dash`, `eq_dash`,  `ra`, `X_d_dash`, `X_q_dash`.
"""
# @doc (@doc get_dynamic_τe_from_id_iq)
function get_dynamic_τe_from_id_iq(
    id, iq, δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash)

    return ed_dash * id + eq_dash * iq +
        ( X_q_dash - X_d_dash ) * id * iq
    
end


"""
    get_dynamic_pg_from_id_iq(
       id, iq, δ, ed_dash, eq_dash,  ra,
       X_d_dash, X_q_dash)


Returns a generator active power  for a given `id`, `iq`, `δ`, `ed_dash`, `eq_dash`,  `ra`, `X_d_dash`, `X_q_dash`.
"""
function get_dynamic_pg_from_id_iq(
    id, iq, δ, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash)

    return ed_dash * id + eq_dash * iq +
        ( X_q_dash - X_d_dash ) * id * iq
    
end

# @doc (@doc get_dynamic_τe_from_id_iq)
"""
    get_dynamic_pg_from_id_iq(
        id, iq, δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)


Returns a generator active power  for a given `id`, `iq`, `δ`, `ed_dash`, `eq_dash`,  `ra`, `X_d_dash`, `X_q_dash`.
"""
function get_dynamic_pg_from_id_iq(
    id, iq, δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash)

    return ed_dash * id + eq_dash * iq +
        ( X_q_dash - X_d_dash ) * id * iq
    
end


"""
    get_dynamic_id_iq_pg_vh_by_vhθh(
        vh, θh, δ, ed_dash, eq_dash, ra,
        X_d_dash, X_q_dash )


Returns a list of `id`, `iq`, `pg` and `vh` for a given `vh`, `θh`, `δ`, `ed_dash`, `eq_dash`, `ra`, `X_d_dash`, `X_q_dash`.
"""
function get_dynamic_id_iq_pg_vh_by_vhθh(
    vh, θh, δ, ed_dash, eq_dash, ra,
    X_d_dash, X_q_dash )

    id_iq = get_dynamic_idq_vhθh(
        vh, θh, δ, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )

    id = id_iq[1]

    iq = id_iq[2]

    pg = get_dynamic_pg_from_id_iq(
        id, iq, δ, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)

    return [ id, iq, pg, vh ]
    
end

# @doc (@doc get_dynamic_id_iq_pg_vh_by_vhθh)
"""
    get_dynamic_id_iq_pg_vh_by_vhθh(
        vh, θh, δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )


Returns a list of `id`, `iq`, `pg` and `vh` for a given `vh`, `θh`, `δ`, `ed_dash`, `eq_dash`, `ra`, `X_d_dash`, `X_q_dash`.
"""
function get_dynamic_id_iq_pg_vh_by_vhθh(
    vh, θh, δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    id_iq = get_dynamic_idq_vhθh(
        vh, θh, δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )

    id = id_iq[1]

    iq = id_iq[2]

    pg = get_dynamic_pg_from_id_iq(
        id, iq, δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)

    return [ id, iq, pg, vh ]
    
end


# @doc (@doc get_dynamic_id_iq_pg_vh_by_vhθh)
"""
    get_gens_dynamic_id_iq_pg_vh_by_vhθh(
        gens_vh_θh_post_pf,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        gen_nodes_ra_Xd_dash_Xq_dash_view )


Returns a list of `id`, `iq`, `pg` and `vh` for a given `vh`, `θh`, `δ`, `ed_dash`, `eq_dash`, `ra`, `X_d_dash`, `X_q_dash`.
"""
function get_gens_dynamic_id_iq_pg_vh_by_vhθh(
    gens_vh_θh_post_pf,
    gen_nodes_δ_ω_ed_dash_eq_dash_views,
    gen_nodes_ra_Xd_dash_Xq_dash_view )

    return [ get_dynamic_id_iq_pg_vh_by_vhθh(
        vh_θh..., δ_ω_ed_dash_eq_dash...,
        ra_Xd_dash_Xq_dash...)
             for (vh_θh, δ_ω_ed_dash_eq_dash,
                  ra_Xd_dash_Xq_dash) in
                 zip(gens_vh_θh_post_pf,
                     gen_nodes_δ_ω_ed_dash_eq_dash_views,
                     gen_nodes_ra_Xd_dash_Xq_dash_view ) ]

end

# @doc (@doc get_dynamic_id_iq_pg_vh_by_vhθh)
"""
    get_dynamic_id_iq_pg_vh_by_ur_ui(
        u_r, u_i, δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )


Returns a list of `id`, `iq`, `pg` and `vh` for a given `vh`, `θh`, `δ`, `ed_dash`, `eq_dash`, `ra`, `X_d_dash`, `X_q_dash`.
"""
function get_dynamic_id_iq_pg_vh_by_ur_ui(
    u_r, u_i, δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    id_iq = get_dynamic_idq_ur_ui(
        u_r, u_i, δ, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)

    pg = get_dynamic_pg_from_id_iq(
        id_iq..., δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)

    return [ [ id_iq, pg, abs( u_r + im * u_i ) ]...; ]
    
end


# @doc (@doc get_dynamic_id_iq_pg_vh_by_vhθh)
"""
    get_dynamic_id_iq_pg_vh_by_ur_ui(
        u_r, u_i, δ,
        ed_dash,
        eq_dash,
        ra,
        X_d_dash,
        X_q_dash )


Returns a list of `id`, `iq`, `pg` and `vh` for a given `vh`, `θh`, `δ`, `ed_dash`, `eq_dash`, `ra`, `X_d_dash`, `X_q_dash`.
"""
function get_dynamic_id_iq_pg_vh_by_ur_ui(
    u_r, u_i, δ,
    ed_dash,
    eq_dash,
    ra,
    X_d_dash,
    X_q_dash )

    id_iq = get_dynamic_idq_ur_ui(
        u_r, u_i, δ, ed_dash, eq_dash,  ra,
        X_d_dash, X_q_dash)

    pg = get_dynamic_pg_from_id_iq(
        id_iq..., δ, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash)

    return [ [ id_iq, pg, abs( u_r + im * u_i ) ]...; ]
    
end

# @doc (@doc get_dynamic_id_iq_pg_vh_by_vhθh)
"""
    get_gens_dynamic_id_iq_pg_vh_by_ur_ui(
        gens_ur_ui_post_pf,
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        gen_nodes_ra_Xd_dash_Xq_dash_view )


Returns a list of `id`, `iq`, `pg` and `vh` for a given `vh`, `θh`, `δ`, `ed_dash`, `eq_dash`, `ra`, `X_d_dash`, `X_q_dash`.
"""
function get_gens_dynamic_id_iq_pg_vh_by_ur_ui(
    gens_ur_ui_post_pf,
    gen_nodes_δ_ω_ed_dash_eq_dash_views,
    gen_nodes_ra_Xd_dash_Xq_dash_view )

    return [
        get_dynamic_id_iq_pg_vh_by_ur_ui(
            ur_ui...,
            δ_ω_ed_dash_eq_dash...,
            ra_Xd_dash_Xq_dash...)
        for (ur_ui, δ_ω_ed_dash_eq_dash,
             ra_Xd_dash_Xq_dash) in
            zip(gens_ur_ui_post_pf,
                gen_nodes_δ_ω_ed_dash_eq_dash_views,
                gen_nodes_ra_Xd_dash_Xq_dash_view) ]
    
end

#----------------------------------------


"""
    get_dynamic_ph_by_vhθh_δ_idq(
        vh, θh, δ, id, iq )


Returns a generator active power for a given `vh`, `θh`, `δ`, `id`, and `iq`.
"""
function get_dynamic_ph_by_vhθh_δ_idq(
    vh, θh, δ, id, iq )

    return id * vh * sin( δ - θh ) +
        iq * vh * cos( δ - θh )
    
end


"""
    get_dynamic_ph_by_vhθh_δ_idq(
        vh, θh, δ, ω, ed_dash,
        eq_dash, id, iq )


Returns a generator active power for a given `vh`, `θh`, `δ`, `ω`, `ed_dash`, `eq_dash`, `id`, and `iq`.
"""
function get_dynamic_ph_by_vhθh_δ_idq(
    vh, θh, δ, ω, ed_dash,
    eq_dash, id, iq )

    return id * vh * sin( δ - θh )  +
        iq * vh * cos( δ - θh )
    
end


"""
    get_dynamic_qh_by_vhθh_δ_idq(
        vh, θh, δ, id, iq )


Returns a generator reactive power for a given `vh`, `θh`, `δ`, `id`, and `iq`.
"""
function get_dynamic_qh_by_vhθh_δ_idq(
    vh, θh, δ, id, iq )

    return id * vh * cos( δ - θh ) -
        iq * vh * sin( δ - θh )
    
end


"""
    get_dynamic_qh_by_vhθh_δ_idq(
        vh, θh, δ, ω, ed_dash,
        eq_dash, id, iq )


Returns a generator active power for a given `vh`, `θh`, `δ`, `ω`, `ed_dash`, `eq_dash`, `id`, and `iq`.
"""
function get_dynamic_qh_by_vhθh_δ_idq(
    vh, θh, δ, ω,
    ed_dash, eq_dash,
    id, iq )

    return id * vh * cos( δ - θh ) -
        iq * vh * sin( δ - θh )
    
end

"""
    get_gens_dynamic_ph_by_vhθh_δ_idq(
        gens_vh, gens_θh,
        gens_δ,
        gens_id, gens_iq )

Returns generators active power.
"""
function get_gens_dynamic_ph_by_vhθh_δ_idq(
    gens_vh, gens_θh,
    gens_δ,
    gens_id, gens_iq )

    return [
        get_dynamic_ph_by_vhθh_δ_idq(vh, θh, δ, id, iq)
             for (vh, θh, δ, id, iq) in
                 zip(gens_vh,
                     gens_θh,
                     gens_δ,
                     gens_id,
                     gens_iq ) ]
    
end


"""
    get_gens_dynamic_ph_by_vhθh_δ_idq(
        gens_vh_θh,
        gens_nodes_δ_ω_ed_dash_eq_dash,
        gens_id_iq )


Returns generators active power.
"""
function get_gens_dynamic_ph_by_vhθh_δ_idq(
    gens_vh_θh,
    gens_nodes_δ_ω_ed_dash_eq_dash,
    gens_id_iq )

    return [ get_dynamic_ph_by_vhθh_δ_idq(
        vh_θh...,
        δ_ω_ed_dash_eq_dash...,
        id_iq...)
             for (vh_θh, δ_ω_ed_dash_eq_dash, id_iq) in
                 zip(gens_vh_θh,
                     gens_nodes_δ_ω_ed_dash_eq_dash,
                     gens_id_iq) ]
    
end



"""
    get_gens_dynamic_qh_by_vhθh_δ_idq((
        gens_vh, gens_θh,
        gens_δ,
        gens_id, gens_iq )


Returns generators active power.
"""
function get_gens_dynamic_qh_by_vhθh_δ_idq(
    gens_vh, gens_θh,
    gens_δ,
    gens_id,
    gens_iq )

    return [ get_dynamic_qh_by_vhθh_δ_idq(
        vh, θh, δ, id, iq)
             for (vh, θh, δ, id, iq) in
                 zip(gens_vh, gens_θh,
                     gens_δ, gens_id,
                     gens_iq ) ]
    
end

"""
    get_gens_dynamic_qh_by_vhθh_δ_idq(
        gens_vh_θh,
        gens_nodes_δ_ω_ed_dash_eq_dash,
        gens_id_iq )


Returns generators active power.
"""
function get_gens_dynamic_qh_by_vhθh_δ_idq(
    gens_vh_θh,
    gens_nodes_δ_ω_ed_dash_eq_dash,
    gens_id_iq )

    return [ get_dynamic_qh_by_vhθh_δ_idq(
        vh_θh...,
        δ_ω_ed_dash_eq_dash...,
        id_iq...)
             for (vh_θh, δ_ω_ed_dash_eq_dash, id_iq) in
                 zip(gens_vh_θh,
                     gens_nodes_δ_ω_ed_dash_eq_dash,
                     gens_id_iq) ]
    
end



function get_vh_θh_post_pf(
    power_flow_data )

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    return ( vh, θh )
    
end


function get_ur_ui_post_pf(
    power_flow_data )

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    uh = vh * exp( im * θh )
    return ( real(uh), imag(uh) )
    
end

#------------------------------------------------

function get_gens_ur_ui(
    nodes_pf_U_view, gens_idx )


    return nodes_pf_U_view[gens_idx]
    
end


function get_gens_vh_θh(
    nodes_pf_U_view, gens_idx )

    return [[abs(uh), angle(uh)]
            for uh in x_from_xr_xi.(
                nodes_pf_U_view[gens_idx]  ) ]
    
end


#------------------------------------------------

function get_non_gens_ur_ui(
    nodes_pf_U_view, non_gens_idx )


    return nodes_pf_U_view[non_gens_idx]
    
end


function get_non_gens_vh_θh(
    nodes_pf_U_view, non_gens_idx )

    return [[abs(uh), angle(uh)]
            for uh in x_from_xr_xi.(
                nodes_pf_U_view[ non_gens_idx ] )]
    
end

#------------------------------------------------

function get_component_δ_ed_dash_eq_dash_from_pf(
    vh, θh, Pg, Qg, ra, X_d, X_q, X_d_dash, X_q_dash )

    Sg = Pg + im * Qg
    
    uh = vh * exp(im * θh)
    
    Ig = conj( Sg ) / conj( uh )

    E = uh + (ra + im * X_q) * Ig

    δ = angle(E)

    id_iq = Ig * exp( -im * ( δ - pi/2 ))

    id = real(id_iq)

    iq = imag(id_iq)

    vd_vq = uh * exp( -im * ( δ - pi/2 ))

    vd = real( vd_vq )

    vq = imag( vd_vq )

    ed_dash = ( X_q -  X_q_dash ) * iq

    check_ed_dash = vd + ra * id - X_q_dash * iq
    
    eq_dash = vq + ra * iq + X_d_dash * id

    return [δ, ed_dash, eq_dash] 
    

end


# function get_component_δ_ω_ed_dash_eq_dash_from_pf(
#     vh, θh, Pg, Qg, ra, X_d, X_q, X_d_dash, X_q_dash )

#     Sg = Pg + im * Qg
    
#     uh = vh * exp(im * θh)
    
#     Ig = conj( Sg ) / conj( uh )

#     E = uh + (ra + im * X_q) * Ig

#     δ = angle(E)

#     id_iq = Ig * exp( -im * ( δ - pi/2 ))

#     id = real(id_iq)

#     iq = imag(id_iq)

#     vd_vq = uh * exp( -im * ( δ - pi/2 ))

#     vd = real( vd_vq )

#     vq = imag( vd_vq )

#     ed_dash = ( X_q -  X_q_dash ) * iq

#     check_ed_dash = vd + ra * id - X_q_dash * iq
    
#     eq_dash = vq + ra * iq + X_d_dash * id

#     return [δ, ωs, ed_dash, eq_dash] 
    

# end



function get_component_id_iq_from_pf(
    vh, θh, Pg, Qg, ra, X_d, X_q,
    X_d_dash, X_q_dash )

    Sg = Pg + im * Qg
    
    uh = vh * exp(im * θh)
    
    Ig = conj( Sg ) / conj( uh )

    E = uh + (ra + im * X_q) * Ig

    δ = angle(E)

    id_iq = Ig * exp( -im * ( δ - pi/2 ))

    id = real(id_iq)

    iq = imag(id_iq)

    return [ id, iq ] 
    

end

#------------------------------------------------

function get_component_δ_id_iq_vd_vq_ed_dash_eq_dash_from_pf(vh, θh, Pg, Qg, ra, X_d, X_q, X_d_dash, X_q_dash; ω = ωs )

    Sg = Pg + im * Qg
    
    uh = vh * exp(im * θh)
    
    Ig = conj( Sg ) / conj( uh )

    E = uh + (ra + im * X_q) * Ig

    δ = angle(E)

    id_iq = Ig * exp( -im * ( δ - pi/2 ))

    id = real(id_iq)

    iq = imag(id_iq)

    vd_vq = uh * exp( -im * ( δ - pi/2 ))

    vd = real( vd_vq )

    vq = imag( vd_vq )

    ed_dash = ( X_q -  X_q_dash ) * iq

    check_ed_dash = vd + ra * id - X_q_dash * iq
    
    eq_dash = vq + ra * iq + X_d_dash * id

    Vf = eq_dash + (X_d - X_d_dash) * id

    τe = ed_dash * id + eq_dash * iq +
        (X_q_dash - X_d_dash) * id * iq

    # ω = 1.0

    return (δ, ed_dash, eq_dash), (
        τe, Vf, id, iq, vd, vq, check_ed_dash)
    
end


function get_component_δ_ω_ed_dash_eq_dash_from_pf(
    vh, θh, Pg, Qg, ra, X_d, X_q, X_d_dash,
    X_q_dash; ω = ωs )

    Sg = Pg + im * Qg
    
    uh = vh * exp(im * θh)
    
    Ig = conj( Sg ) / conj( uh )

    E = uh + (ra + im * X_q) * Ig

    δ = angle(E)

    id_iq = Ig * exp( -im * ( δ - pi/2 ))

    id = real(id_iq)

    iq = imag(id_iq)

    vd_vq = uh * exp( -im * ( δ - pi/2 ))

    vd = real( vd_vq )

    vq = imag( vd_vq )

    ed_dash = ( X_q -  X_q_dash ) * iq

    check_ed_dash = vd + ra * id - X_q_dash * iq
    
    eq_dash = vq + ra * iq + X_d_dash * id

    Vf = eq_dash + (X_d - X_d_dash) * id

    τe = ed_dash * id + eq_dash * iq +
        (X_q_dash - X_d_dash) * id * iq

    # ω = 1.0

    return (δ, ω, ed_dash, eq_dash), (
        τe, Vf), (id, iq, vd, vq)
    

end

#------------------------------------------------

function get_dynamic_τm_vf(
    vh, θh, δ, ω, ed_dash, eq_dash,
    ra, X_d, X_q, X_d_dash, X_q_dash  )

    id, iq = get_dynamic_idq_vhθh(
        vh, θh, δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )

    τe = ed_dash * id + eq_dash * iq +
        (X_q_dash - X_d_dash) * id * iq
    
    vf = eq_dash + (X_d - X_d_dash) * id

    return [τe, vf]
        
end


function get_dynamic_τm_vf(
    vh, θh, δ, ed_dash, eq_dash, ra, X_d,
    X_q, X_d_dash, X_q_dash  )

    id, iq = get_dynamic_idq_vhθh(
        vh, θh, δ, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )

    τe = ed_dash * id + eq_dash * iq +
        (X_q_dash - X_d_dash) * id * iq
    
    vf = eq_dash + (X_d - X_d_dash) * id

    return [τe, vf]
        
end


function get_gens_dynamic_τm_vf(
    gens_vh_θh_view,
    gen_nodes_δ_ω_ed_dash_eq_dash_views,
    gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view  )

    return [ get_dynamic_τm_vf(
        vh_θh..., δ_ω_ed_dash_eq_dash...,
        ra_Xd_Xq_Xd_dash_Xq_dash... )
             for (vh_θh, δ_ω_ed_dash_eq_dash,
                  ra_Xd_Xq_Xd_dash_Xq_dash) in
                 zip(gens_vh_θh_view,
                     gen_nodes_δ_ω_ed_dash_eq_dash_views,
                     gen_nodes_ra_Xd_Xq_Xd_dash_Xq_dash_view )]
        
end


#------------------------------------------------
#------------------------------------------------


function get_components_ur_ui_in_state(
    comps_flattened_states,
    ur_ui_Idx_in_state  )

    return [
        comps_flattened_states[inx]
        for inx in
            ur_ui_Idx_in_state ]
end

function get_components_u_in_state(
    comps_flattened_states,
    ur_ui_Idx_in_state  )

    return ur_ui_to_u.(
        [ comps_flattened_states[inx]
          for inx in
              ur_ui_Idx_in_state ])    
end


#------------------------------------------------

"""
    syms_containing(
        network_vars_labels,
        expr;
        syms = true)


Returns the list of symbols of variables containing a given `expr` in network variable labels.
"""
function syms_containing(
    network_vars_labels,
    expr;
    syms = true)
    
    if typeof(expr) == String
        return [
            s for s in
                network_vars_labels
                if occursin(expr, string(s))]
    else
        return [
            s for s in
                network_vars_labels
                if occursin(string(expr), string(s))]
    end
    
end


"""
    idx_containing(
        network_vars_labels,
        expr;
        syms = true)


Returns the list of indices of variables containing a given `expr` in network variable labels.
"""
function idx_containing(
    network_vars_labels,
    expr;
    syms = true)
    
    if typeof(expr) == String
        
        return [
            i for (i, s) in
                enumerate(network_vars_labels)
                if occursin(expr, string(s))]
    else
        return [
            i for (i, s) in
                enumerate(network_vars_labels)
                if occursin(string(expr), string(s))]
    end
    
end

#------------------------------------------------

"""
    generate_net_bus_volts_labels(
        network_bus_names)

Returns buses voltage labels for real part and imaginary part od buses voltages.
"""
function generate_net_bus_volts_labels(
    network_bus_names)

    comps_bus_u_labels = Symbol[]
    
    for a_bus in network_bus_names
        volts_syms = [:u_r, :u_i]
        append!( comps_bus_u_labels,
                 [Symbol(a_bus, "_", sym )
                  for sym in volts_syms ])
    end
    
    return comps_bus_u_labels
        
end

#------------------------------------------------

"""
    get_vars_idxs_in_range_Idxs(
        vec_unit_range,
        list_vars_syms )


Returns vectors of idxs as tuples

gens_nodes_idx = [1,2,3]
gens_state_vars_syms = [:δ, :ω, :eq_dash, :E_fd]

gens_states_Idx =
    get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
        gens_state_vars_syms,
        gens_nodes_idx  )

sum( length.( gens_states_Idx ) )

(δ_idx_in_state,
ω_idx_in_state,
eq_dash_idx_in_state,
E_fd_idx_in_state ) =
    get_vars_idxs_in_range_Idxs(
        gens_states_Idx,
        gens_state_vars_syms )

"""
function get_vars_idxs_in_range_Idxs(
    vec_unit_range,
    list_vars_syms )

    # vec_unit_range = collect.(vec_unit_range)

    dims = length.(vec_unit_range)

    @assert allequal(dims)

    @assert length(list_vars_syms) == dims[1]

    idxs_in_Idxs = Vector{Int64}[ [] for idx in
                        1:length(list_vars_syms) ]

    for idx in 1:length(list_vars_syms)
        for a_vec_idx in 1:length(vec_unit_range)

            push!(idxs_in_Idxs[idx],
                  vec_unit_range[a_vec_idx][ idx ])
        end
        
    end

    return namedtuple(
        OrderedDict( a_sym => vec_idx
             for (a_sym, vec_idx) in
                 zip(list_vars_syms, idxs_in_Idxs) ))
    
    
end

#-------------------------------------------------------

"""
    get_non_null_list(list_of_lists)
        return [ a_list for a_list in
                    list_of_lists
                    if a_list != []]


Returns non empty lists.

ll_1 = [[9,2], [4,5,7,3], [3,6,8,4]]

 get_non_null_list(ll_1)

ll_2 = [[9,2], [], [3,6,8,4]]

get_non_null_list(ll_2)


"""
function get_non_null_list(list_of_lists)
    return [ a_list for a_list in
                list_of_lists
                if a_list != []]
end


"""
    get_non_null_list_and_Idx(
        list_of_lists)


Returns non empty lists and their indices.

ll_1 = [[9,2], [4,5,7,3], [3,6,8,4]]

get_non_null_list_and_Idx(ll_1)

ll_2 = [[9,2], [], [3,6,8,4]]

get_non_null_list_and_Idx(ll_2)

"""
function get_non_null_list_and_Idx(list_of_lists)
    
    Idx_and_a_list = [ (Idx, a_list)
                       for (Idx, a_list) in
                           enumerate(list_of_lists)
                           if a_list != []]
    
    return (Indices = first.(Idx_and_a_list),
            lists = last.(Idx_and_a_list))
end



# ------------------------------------------------------
#  Idx creation functions
# ------------------------------------------------------

"""
    create_offsets(
        dims;
        counter=0)::Vector{Int}


Create offsets for stacked array of dimensions dims.
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
function create_idxs(offs, dims)::Vector{UnitRange{Int} }
    
    idxs = [1+off:off+dim
            for (off, dim) in
                zip(offs, dims)]
end


#------------------------------------------

"""
    get_ur_ui_Idx_in_state(
        dims, comps_ur_ui_Idx; counter=0 )


Returns indices of u_r, u_i for a flattend system states.

The input are dims of components state and
indices of u_r and u_i in each component state

Test:
```math
```
comps_states_syms = [[:u_r,:u_i], [:δ, :ω, :u_r, :u_i],
    [:δ,:u_r,:ω,:u_i]]

comps_states = [[9,2], [4,5,7,3],
    [3,6,8,4]]

comps_flattened_states = [comps_states...;]

dims = length.(comps_states)

comps_ur_u_i_Idx = [[1,2], [3,4], [2,4]]

ur_ui_Idx_in_state = get_ur_ui_Idx_in_state(dims,
    comps_ur_u_i_Idx; counter=0 )

ur_ui_in_state = [comps_flattened_states[inx]
    for inx in
        ur_ui_Idx_in_state ]

comps_u_view  = [view(comps_flattened_states,
    ur_ui_Idx_in_state[Ind])
        for Ind in
            collect(1:length( ur_ui_Idx_in_state ))]

"""
function get_ur_ui_Idx_in_state(
    dims, comps_ur_ui_Idx; counter=0 )

    offs = create_offsets(dims; counter=counter)
    
    return [no == 1 ?
        [idx[1], idx[2]] :
        [idx[1]+off, idx[2]+off]
            for (no, off, dim, idx) in
                zip(collect(1:length(offs)),
                    offs, dims, comps_ur_ui_Idx)]

end


function get_non_slack_ur_ui_Idx_in_state(
    dims, comps_ur_ui_Idx; counter=0 )

    offs = create_offsets(dims; counter=counter)
    
    return [no == 1 ?
        [idx[1]+counter, idx[2]+counter] :
        [idx[1]+off, idx[2]+off]
            for (no, off, dim, idx) in
                zip(collect(1:length(offs)),
                    offs, dims, comps_ur_ui_Idx)]

end

"""
    create_size_dims_offset_Idx(comps_dims)

Create indexes  of dimensions dims.
"""
function create_size_offset_Idx(comps_dims; counter = 0)
    
    comps_size    = sum(comps_dims)
    comps_offset  = create_offsets(comps_dims; counter)
    comps_Idx     = create_idxs(comps_offset, comps_dims)

    return comps_size, comps_offset, comps_Idx
end


"""
    create_u_idxs(offs)

Create indexes of u_r and u_i for stacked array
of dimensions dims using the offsets offs
"""
function create_u_idxs(offs)::Vector{UnitRange{Int}}
    # idxs = [1+off:off+2 for off in offs]

    return [1+off:off+2 for off in offs]
end


"""
    create_i_src_idxs(offs)

Create indexes of src-end (h) flow in an edge, xh_r
and xh_i for stacked array of dimensions dims
using the offsets offs
"""
function create_x_src_idxs(offs)::Vector{UnitRange{Int}}
    # idxs = [1+off:off+2 for off in offs]

    return [1+off:off+2 for off in offs]
end


"""
    create_x_dst_idxs(offs)

Create indexes of dst-end (h) flow in an edge,
xk_r and xk_i for stacked array of dimensions
dims using the offsets offs
"""
function create_x_dst_idxs(offs)::Vector{UnitRange{Int}}
    # idxs = [1+off:off+2 for off in offs]

    return [3+off:off+4 for off in offs]
end


#----------------------------------------

"""
    convert_to_consecutive_idxs(
        flat_vh_flat_θh_Idx )


Converts non consecutive idxs to a consecutive idxs.
"""
function convert_to_consecutive_idxs(
    flat_vh_flat_θh_Idx )
    
    (flat_vh_idx_in_flat_vh_flat_θh,
     flat_θh_idx_in_flat_vh_flat_θh) =
         flat_vh_flat_θh_Idx
    
    vec_a_vh_a_θh_idxs =
        [[a_vh, a_θh ]
         for (a_vh, a_θh) in
             zip( flat_vh_idx_in_flat_vh_flat_θh,
                  flat_θh_idx_in_flat_vh_flat_θh )  ]
    
    flat_consecutive_idxs = [vec_a_vh_a_θh_idxs...;]
    
    return (;flat_consecutive_idxs,
            vec_a_vh_a_θh_idxs )
    
end

"""
    convert_to_non_consecutive_idxs(
        flat_vh_idx_flat_θh_idx_in_flat_vh_θh )


Converts consecutive idxs to non consecutive idxs.
"""
function convert_to_non_consecutive_idxs(
    flat_vh_idx_flat_θh_idx_in_flat_vh_θh )
    
    (flat_vh_idx_in_flat_vh_θh,
     flat_θh_idx_in_flat_vh_θh) =
        flat_vh_idx_flat_θh_idx_in_flat_vh_θh
    
    vec_a_vh_a_θh_idxs =
        [[a_vh, a_θh ]
         for (a_vh, a_θh) in
             zip( flat_vh_idx_in_flat_vh_θh,
                  flat_θh_idx_in_flat_vh_θh )  ]
    
    flat_non_consecutive_idxs =
        [flat_vh_idx_in_flat_vh_θh;
         flat_θh_idx_in_flat_vh_θh ]
    
    return (;flat_non_consecutive_idxs,
            vec_a_vh_a_θh_idxs )
    
end


#------------------------------------------
#------------------------------------------

"""
    DAE_MassMatrix(
        state_size, algebraic_size )

Returns the mass matrix for a differential algebraic equations
"""
function DAE_MassMatrix(
    state_size, algebraic_size )
    
    return Diagonal(
        [ones(Int, state_size)...;
         zeros(Int, algebraic_size)...])
end

"""
    DAE_BoolVector(
        state_size, algebraic_size )

Returns the boolean vector for a differential algebraic equations
"""
function DAE_BoolVector(
    state_size, algebraic_size)
    
    return map(
        (x) -> x ==1 ? true : false,
        [ones(Int, state_size)...;
         zeros(Int, algebraic_size)...] )

end

#------------------------------------------


"""
    Sevf(Ae, Be, vf_tilade)

Returns an excitation.
"""
function Sevf(Ae, Be, vf_tilade)
    
    return Ae * exp(Be * abs(vf_tilade))
    
end


"""
    exciter_saturation_function(
        ( efd_1, S_E_1 ), ( efd_2, S_E_2 ) )

Returns an exciter saturation parameters `A_x`, and `B_x`.


Power system modeling, computation and control.

page:243

"""
function exciter_saturation_function(
    ( efd_1, S_E_1 ), ( efd_2, S_E_2 ) )

    Q = [log( efd_1 * S_E_1 ),
         log( efd_2 * S_E_2 )]
    
    P = [1    efd_1;
         1    efd_2]
    
    x = P \ Q
    
    ln_A_x = x[1]
    
    B_x = x[2]
    A_x = exp(ln_A_x)
    
    return A_x, B_x
end


"""
    exciter_saturation_function(
        (K_e, V_R_max), S_E_max, S_E0_75max)


Returns an exciter saturation parameters `A_x`, and `B_x`.
"""
function exciter_saturation_function(
    (K_e, V_R_max), S_E_max, S_E0_75max)
    
    # SE (Efd) = Ax * e ^(Bx * Efd)
    E_fd_max = V_R_max/(K_e + S_E_max)

    Q = [log(E_fd_max * S_E_max),
         log(0.75 * E_fd_max *  S_E0_75max)]
    
    P = [1   E_fd_max; 1   0.75 * E_fd_max]
    
    x = P \ Q
    
    ln_A_x = x[1]
    
    B_x    = x[2]
    A_x    = exp(ln_A_x)
    
    return A_x, B_x
end


"""
    exciter_saturation_function(
        K_e, V_R_max, S_E_max, S_E0_75max)


Returns an exciter saturation parameters `A_x`, and `B_x`.
"""
function exciter_saturation_function(
    K_e, V_R_max, S_E_max, S_E0_75max)
    
    # SE (Efd) = Ax * e ^(Bx * Efd)
    
    E_fd_max = V_R_max/(K_e + S_E_max)
    
    B_x = (log(S_E_max/S_E0_75max)) /
        (E_fd_max - 0.75 * E_fd_max)
    
    A_x = S_E_max * exp(-B_x * E_fd_max)

    return A_x, B_x
end


#-----------------------------------------------------

"""
    dict_reverse_keys_values_pair(
        a_dict)


Returns a reverse pair of a dictionary.
"""
function dict_reverse_keys_values_pair(
    a_dict)

    dict_keys   = keys(a_dict)
    dict_values = values(a_dict)

    return OrderedDict(
        value => key for (value, key) in
            zip(dict_values, dict_keys) )
end



function get_a_flattened_by_per_node(
    list_vars_or_paras)

    per_node_vars_or_paras = [
        [ [ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para]...;]
                  for var_or_para in
                      zip( list_vars_or_paras...) ]

    return [per_node_vars_or_paras...;]
    

end


"""
    get_a_flattened_by_per_vars_or_paras(
        list_vars_or_paras )


Returns for a vector of vectors `list_vars_or_paras` a  flattened vector.
"""
function get_a_flattened_by_per_vars_or_paras(
    list_vars_or_paras)


    per_vars_or_paras = [
        [per_var_or_para...;]
        for per_var_or_para in
            list_vars_or_paras  ]

    return [per_vars_or_paras...;]
    

end


"""
    get_per_vars_or_paras_to_per_node(
        list_vars_or_paras)


a = [1,2,3]
b = [4,5,6]

test_data = [a, b]

test_a = get_a_flattened_by_per_node( test_data )

test_b = get_a_flattened_by_per_vars_or_paras( test_data )


d = [ [1,2], [3,4], [5,6]]

e = [[10,11,12], [13,14,15], [16,17,18]]

test_data_2 = [d, e]

test_m = get_a_flattened_by_per_node( test_data_2  )

test_n = get_a_flattened_by_per_vars_or_paras( test_data_2  )

"""



"""
    get_per_vars_or_paras_to_per_node(
        list_vars_or_paras)


Converts vars or paras given in per vars or paras format to per node format.

"""
function get_per_vars_or_paras_to_per_node(
    list_vars_or_paras)

    return [
        [[ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para]...;]
                  for var_or_para in
                      zip( list_vars_or_paras...) ]
    
end


"""
    get_per_node_para_to_per_vars_or_paras(
        list_per_node_vars_or_paras,
        dims_vars_or_paras_types )


Converts vars or paras given in per node format to per vars or paras format.


The dimension of each of the per vars or paras should be supplied in a list `dims_vars_or_paras_types`

"""
function get_per_node_para_to_per_vars_or_paras(
    list_per_node_vars_or_paras,
    dims_vars_or_paras_types )
     
    _, _, per_node_vars_or_paras_Idx =
        create_size_offset_Idx(
            dims_vars_or_paras_types;
            counter = 0 )

    nodes_vars_or_paras_per_type = [
        [ a_node_vars_or_paras[idx]
              for idx in
                  per_node_vars_or_paras_Idx]
      for a_node_vars_or_paras in
          list_per_node_vars_or_paras ]
    

    return  [[ a_vars_or_paras_type
               for a_vars_or_paras_type in
                   a_node_vars_or_paras_type ]
             for a_node_vars_or_paras_type in
                 zip(nodes_vars_or_paras_per_type...) ]
        
end



"""
    get_a_flattend_vars_or_paras_and_Idx(
        vec_of_vec_var )


Returns indices in a flattened vectors for a vector of vectors.
 

d = [ [1,2], [3,4], [5,6]]

e = [[10,11,12], [13,14,15], [16,17,18]]

test_data_2 = [d, e]

dims_vars_or_paras_types = length.([first(d), first(e) ] )

test_e = get_per_vars_or_paras_to_per_node( test_data_2  )

test_t = get_per_node_para_to_per_vars_or_paras(
     test_e,
    dims_vars_or_paras_types )

"""
function get_a_flattend_vars_or_paras_and_Idx(
    vec_of_vec_var )

    #-----------------------------------------------
    
    vars_or_paras_Idx =
        get_flattened_to_components_vector_var_Idx(
            vec_of_vec_var )
    
    #-----------------------------------------------

    flattened_vars_or_paras = [vec_of_vec_var...;]

    #-----------------------------------------------
    
    return flattened_vars_or_paras, vars_or_paras_Idx
    
end


"""
    get_per_node_flat_idxs(
        list_vars_or_paras )


Returns for a vector of vectors `list_vars_or_paras` a  per sub vector indices  as `per_node_vars_or_paras_Idx`.
"""
function get_per_node_flat_idxs(
    list_vars_or_paras )

    per_node_vars_or_paras = [
        [ [ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para]...;]
                  for var_or_para in
                      zip( list_vars_or_paras...) ]
    
    dims_list_vars_or_paras =
        length.( per_node_vars_or_paras )
    
    _, _, per_node_vars_or_paras_Idx =
        create_size_offset_Idx(
            dims_list_vars_or_paras;
            counter = 0 )

    return per_node_vars_or_paras_Idx
    
end


"""
    get_per_node_flat_para_and_idxs(
        list_vars_or_paras )


Returns for a vector of vectors `list_vars_or_paras` a flattened vector, per sub vector indices  as `flattend_per_node_vars_or_paras`,  and `per_node_vars_or_paras_Idx` respectively.
"""
function get_per_node_flat_para_and_idxs(
    list_vars_or_paras )

    per_node_vars_or_paras = [
        [ [ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para]...;]
                  for var_or_para in
                      zip( list_vars_or_paras...) ]
    
    #----------------------------------------    
    
    per_node_vars_or_paras_Idx =
        get_per_node_flat_idxs( list_vars_or_paras )
    

    flattend_per_node_vars_or_paras =
        [per_node_vars_or_paras...;]
    
    #----------------------------------------    

    return (;
            flattend_per_node_vars_or_paras,
            per_node_vars_or_paras_Idx )
    
end


"""
    get_per_node_per_vars_or_paras_flat_idxs(
        list_vars_or_paras )


Returns for a vector of vectors `list_vars_or_paras` per sub vector indices and sub vector indices in the flattened vector as `per_node_per_vars_or_paras_Idxs` and `per_node_vars_or_paras_Idx` respectively.
"""
function get_per_node_per_vars_or_paras_flat_idxs(
    list_vars_or_paras )
    
    dims_list_vars_or_paras = [
        length.([ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para])
                  for var_or_para in
                      zip( list_vars_or_paras...) ]
    
    per_node_per_vars_or_paras_size_offset_Idx =
        map((dim) ->
        create_size_offset_Idx(dim; counter = 0 ),
            dims_list_vars_or_paras )

    per_node_per_vars_or_paras_Idxs =
        third.(per_node_per_vars_or_paras_size_offset_Idx)

    
    #----------------------------------------        
    #----------------------------------------    

   _, per_node_vars_or_paras_Idx  =
         get_per_node_flat_para_and_idxs(
             list_vars_or_paras )

    per_node_per_vars_or_paras_Idxs =
        per_node_per_vars_or_paras_Idxs[1]
    
    #----------------------------------------    
    
    return (;                        
            per_node_vars_or_paras_Idx,
            per_node_per_vars_or_paras_Idxs )
    
end


"""
    get_per_node_per_vars_or_paras_flat_para_and_idxs(
        list_vars_or_paras )


Returns for a vector of vectors `list_vars_or_paras` a flattened vector, per sub vector indices and sub vector indices in the flattened vector as `flattend_per_node_vars_or_paras`, `per_node_vars_or_paras_Idx`, `per_node_per_vars_or_paras_Idxs` respectively.
"""
function get_per_node_per_vars_or_paras_flat_para_and_idxs(
    list_vars_or_paras )
    
    dims_list_vars_or_paras = [
        length.([ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para])
                  for var_or_para in
                      zip( list_vars_or_paras...) ]
    
    per_node_per_vars_or_paras_size_offset_Idx =
        map((dim) ->
        create_size_offset_Idx(dim; counter = 0 ),
            dims_list_vars_or_paras )

    per_node_per_vars_or_paras_Idxs =
        third.(per_node_per_vars_or_paras_size_offset_Idx)

    
    #----------------------------------------    

    per_node_per_vars_or_paras = [
        [ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para]
                  for var_or_para in
                      zip( list_vars_or_paras...) ]
    
    #----------------------------------------    

    (; flattend_per_node_vars_or_paras,
     per_node_vars_or_paras_Idx ) =
         get_per_node_flat_para_and_idxs(
             list_vars_or_paras )

    per_node_per_vars_or_paras_Idxs =
        per_node_per_vars_or_paras_Idxs[1]
    
    #----------------------------------------    
    
    return (;            
            flattend_per_node_vars_or_paras,
            per_node_vars_or_paras_Idx,
            per_node_per_vars_or_paras_Idxs )
    
end

"""
    get_per_vars_or_paras_flat_idxs(
        list_vars_or_paras )

Returns indices per parameter or variable `per_vars_or_paras_Idx` for a flattened vector of vectors.
"""
function get_per_vars_or_paras_flat_idxs(
    list_vars_or_paras )

    per_vars_or_paras = [
        [per_var_or_para...;]
        for per_var_or_para in
            list_vars_or_paras  ]
    
    dims_per_vars_or_paras = length.(per_vars_or_paras )
    

    _, _, per_vars_or_paras_Idx =
        create_size_offset_Idx(
            dims_per_vars_or_paras;
            counter = 0 )
    
    #----------------------------------------    
    #----------------------------------------    

    return per_vars_or_paras_Idx 
    
end


"""
    get_per_vars_or_paras_flat_para_and_idxs(
        list_vars_or_paras )


Returns for a vector of vectors a flattened vector and and per sub vector indices as `flattend_per_vars_or_paras`, and `per_vars_or_paras_Idx` respectively.
"""
function get_per_vars_or_paras_flat_para_and_idxs(
    list_vars_or_paras )

    per_vars_or_paras = [
        [per_var_or_para...;]
        for per_var_or_para in
            list_vars_or_paras  ]
    
    dims_per_vars_or_paras = length.(per_vars_or_paras )
    

    _, _, per_vars_or_paras_Idx =
        create_size_offset_Idx(
            dims_per_vars_or_paras;
            counter = 0 )
    
    #----------------------------------------    

    flattend_per_vars_or_paras =
        [per_vars_or_paras...;]
    
    #----------------------------------------    

    return (;
            flattend_per_vars_or_paras,
            per_vars_or_paras_Idx )
    
end


"""
    get_per_vars_or_paras_per_node_flat_idxs(
        list_vars_or_paras )


Returns namedtuple `per_vars_or_paras_Idx`, and `per_vars_or_paras_per_node_Idx`
"""
function get_per_vars_or_paras_per_node_flat_idxs(
    list_vars_or_paras )

    # list_vars_or_paras = [gens_vh_θh,
    #  gens_nodes_ωs_ωref0_vref0_porder0  ]
    
    dims_list_vars_or_paras = [
        length.([ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para])
                  for var_or_para in
                      zip( list_vars_or_paras...) ]


    dims_of_each_vars_or_paras =
        dims_list_vars_or_paras[1]
    
    no_vars_or_paras_in_list =
        length(list_vars_or_paras)
    
    per_node_per_vars_or_paras_size_offset_Idx =
        map((dim) ->
        create_size_offset_Idx(dim; counter = 0 ),
            dims_list_vars_or_paras )

    per_vars_or_paras_per_node_Idx = []

    for dim_a_var_or_para in dims_of_each_vars_or_paras
        
        dim_var_or_para_type = []
        
        for  idx in 1:no_vars_or_paras_in_list

            push!(dim_var_or_para_type,
                      dim_a_var_or_para )
        end
        _,_, var_or_para_type_Idx =
            create_size_offset_Idx(
            dim_var_or_para_type;
                counter = 0 )

        push!(per_vars_or_paras_per_node_Idx ,
                  var_or_para_type_Idx)
    end
    
    
    #---------------------------------------- 

     _, per_vars_or_paras_Idx  =
         get_per_vars_or_paras_flat_para_and_idxs(
             list_vars_or_paras )
    
    #----------------------------------------    

    return (;
            per_vars_or_paras_Idx,
            per_vars_or_paras_per_node_Idx
            )
    
end



"""
    get_per_vars_or_paras_per_node_flat_para_and_idxs(
        list_vars_or_paras )

Returns namedtuples of `flattend_per_vars_or_paras`, `per_vars_or_paras_Idx`, and `per_vars_or_paras_per_node_Idx`.

"""
function get_per_vars_or_paras_per_node_flat_para_and_idxs(
    list_vars_or_paras )

    # list_vars_or_paras = [gens_vh_θh,
    #  gens_nodes_ωs_ωref0_vref0_porder0  ]
    
    dims_list_vars_or_paras = [
        length.([ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para])
                  for var_or_para in
                      zip( list_vars_or_paras...) ]


    dims_of_each_vars_or_paras =
        dims_list_vars_or_paras[1]
    
    no_vars_or_paras_in_list =
        length(list_vars_or_paras)
    
    per_node_per_vars_or_paras_size_offset_Idx =
        map((dim) ->
        create_size_offset_Idx(dim; counter = 0 ),
            dims_list_vars_or_paras )

    per_vars_or_paras_per_node_Idx = []

    for dim_a_var_or_para in dims_of_each_vars_or_paras
        
        dim_var_or_para_type = []
        
        for  idx in 1:no_vars_or_paras_in_list

            push!(dim_var_or_para_type,
                      dim_a_var_or_para )
        end
        _,_, var_or_para_type_Idx =
            create_size_offset_Idx(
            dim_var_or_para_type;
                counter = 0 )

        push!(per_vars_or_paras_per_node_Idx ,
                  var_or_para_type_Idx)
    end
    
    
    #---------------------------------------- 

    (;
     flattend_per_vars_or_paras,
     per_vars_or_paras_Idx ) =
         get_per_vars_or_paras_flat_para_and_idxs(
             list_vars_or_paras )
    
    #----------------------------------------    

    return (;
            flattend_per_vars_or_paras,
            per_vars_or_paras_Idx,
            per_vars_or_paras_per_node_Idx
            )
    
end


# gens_nodes_idx, all_nodes_idx


"""
    get_ode_flat_para_Idxs_in_Idxs(
        gens_vh_θh,
        gens_nodes_ωs_ωref0_vref0_porder0,
        gens_dynamic_id_iq_pg_vh )


Returns namedtuples of `gens_nodes_vh_θh_idx_in_Idx`, `gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx`, `gens_nodes_id_iq_pg_vh_idx_in_Idx`.

"""
function get_ode_flat_para_Idxs_in_Idxs(
    gens_vh_θh,
    gens_nodes_ωs_ωref0_vref0_porder0,
    gens_dynamic_id_iq_pg_vh )

    #-----------------------------------------------

    dims_gens_nodes_vh_θh =
        length.( [gens_nodes_idx, gens_nodes_idx] )
        # length.( gens_vh_θh )

    _,_, gens_nodes_vh_θh_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_vh_θh;
        counter = 0)

    #-----------------------------------------------

    dims_gens_nodes_ωs_ωref0_vref0_porder0 =
        length.([ gens_nodes_idx,
                  gens_nodes_idx,
                  gens_nodes_idx,
                  gens_nodes_idx] )
        # length.( gens_nodes_ωs_ωref0_vref0_porder0 )

    _,_, gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_ωs_ωref0_vref0_porder0;
        counter = 0)

    #-----------------------------------------------

    dims_gens_nodes_id_iq_pg_vh =
        length.([ gens_nodes_idx,
                  gens_nodes_idx,
                  gens_nodes_idx,
                  gens_nodes_idx] )        
        # length.( gens_dynamic_id_iq_pg_vh )

    _,_, gens_nodes_id_iq_pg_vh_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_id_iq_pg_vh;
            counter = 0)

    #-----------------------------------------------

    return (;
            gens_nodes_vh_θh_idx_in_Idx,
            gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
            gens_nodes_id_iq_pg_vh_idx_in_Idx )
    
end


"""
    get_flat_full_vh_θh_and_idxs(
        vh, θh)


Returns a flattened concatenated `full_vh_θh` and indices in the flattened vector.  
"""
function get_flat_full_vh_θh_and_idxs(
    vh, θh )

    vec_full_vh_θh =
        [ vh, θh ]
    
    dims_vec_full_vh_θh  =
        length.( vec_full_vh_θh  )

    _,_, full_vh_θh_Idx =
        create_size_offset_Idx(
            dims_vec_full_vh_θh;
            counter = 0)

    full_vh_Idx, full_θh_Idx = full_vh_θh_Idx

    full_vh_θh = [vec_full_vh_θh...;]

    return  (;
             full_vh_θh,
             full_vh_Idx,
             full_θh_Idx )
        
end


"""
    get_flat_intg_vh_θh_id_iq_and_idxs(
        vh, θh, gens_i_d_0, gens_i_q_0 )


Returns a flattened concatenated `vec_intg_vh_θh_id_iq` and indices in the flattened vector.  
"""
function get_flat_intg_vh_θh_id_iq_and_idxs(
    vh, θh, gens_i_d_0, gens_i_q_0 )


    vec_intg_vh_θh_id_iq =
        [ vh, θh, gens_i_d_0, gens_i_q_0 ]
    
    dims_vec_intg_vh_θh_id_iq  =
        length.( vec_intg_vh_θh_id_iq  )

    _,_, intg_vh_θh_id_iq_Idx =
        create_size_offset_Idx(
            dims_vec_intg_vh_θh_id_iq;
            counter = 0)

    intg_vh_Idx, intg_θh_Idx ,intg_id_Idx, intg_iq_Idx = intg_vh_θh_id_iq_Idx

    intg_vh_θh_id_iq = [vec_intg_vh_θh_id_iq...;]

    return  (;
             intg_vh_θh_id_iq,
             intg_vh_Idx,
             intg_θh_Idx ,
             intg_id_Idx,
             intg_iq_Idx )
        
end


"""
    get_vars_or_paras_Idxs_in_flattend(
        vars_or_param_list_or_dims;
        dims_given = true )


Returns indices of variables or parameters in a flattended vector.
"""
function get_vars_or_paras_Idxs_in_flattend(
    vars_or_param_list_or_dims;
    dims_given = true )

    if dims_given == true
    
        _,_, vars_or_param_Idx_in_flattend =
            create_size_offset_Idx(
                vars_or_param_list_or_dims;
                counter = 0)

        return  vars_or_param_Idx_in_flattend
        
        
    else
    
        dims_vars_or_param_list =
            length.( vars_or_param_list_or_dims )

        _,_, vars_or_param_Idx_in_flattend =
            create_size_offset_Idx(
                dims_vars_or_param_list;
                counter = 0)

        return  vars_or_param_Idx_in_flattend
        
    end
    
end


"""
    gens_and_non_gens_u_Idx_in_ranges(
        all_nodes_idx,
        gens_nodes_idx,
        non_gens_nodes_idx,
        nodes_u_Idx_in_ranges )


Returns generators and non-geenrators voltage indices in form of ranges.
"""
function gens_and_non_gens_u_Idx_in_ranges(
    all_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    nodes_u_Idx_in_ranges )

    dict_nodes_types_u_Idx_in_ranges =
        OrderedDict{Int64, UnitRange{Int64}}(
            net_idx => u_idx_range
            for (net_idx, u_idx_range) in zip(
                all_nodes_idx, nodes_u_Idx_in_ranges ))

    gens_nodes_u_Idx_in_ranges =
        [ dict_nodes_types_u_Idx_in_ranges[idx]
          for idx in gens_nodes_idx ]

    non_gens_nodes_u_Idx_in_ranges =
        [ dict_nodes_types_u_Idx_in_ranges[idx]
          for idx in non_gens_nodes_idx ]


    return (; gens_nodes_u_Idx_in_ranges,
            non_gens_nodes_u_Idx_in_ranges )
                

end


#-------------------------------------------------------

"""
    z2y(; r =1.0, x = 1.00, G = 0.0, B = 0.0)


Returns a transmission admittance `y` for a given `r`, `x`, `G` and `B`.
"""
function z2y(; r =1.0, x = 1.00, G = 0.0, B = 0.0)

    z = r + im * x

    y = 1/z

    B_2 = B/2.0
    

    return (y = y,
            y_shunt_km = B_2,
            y_shunt_mk = B_2)

end


# Untested functions

"""
    V_R2C(V)


Returns a vector of complex voltage for a given concatenated flattened vector of voltage magnitudes and voltage angles.
"""
function V_R2C(V)
    
    size_Θ_Vm = length(X)
    N = Int(size_Θ_Vm/2)

    Vm = @view V[1:N]
    Θ  = @view V[N+1:end]    
    Vc = Vm .* exp.(im * Θ)

    return Vc
 
end


"""
    V_C2R(V)


Returns a tuple of a vector of voltage magnitudes and a vector of voltage angles for a given vector of complex voltage.
"""
function V_C2R(V)

    Vc   = @view V[1:end]

    Θ   = angle.(Vc)
    Vm  = abs.(Vc)

    return (Vm, Θ)
 
end


"""
    Sf(V, Ct, Yt)


Returns complex power `from` nodes.
"""
function Sf( V, Cf, Yf)

    Sf = Diagonal(Cf * Vf) * conj.(Yf) * conj.(V)

    return Sf
end


"""
    St(V, Ct, Yt)


Returns complex power `to` nodes.
"""
function St(V, Ct, Yt)

    St = Diagonal(Ct * Vt) * conj.(Yt) * conj.(V)

    return St
end


"""
    Sbus(V, Ybus)


Returns complex power at nodes.
"""
function Sbus(V, Ybus)
    size_Θ_Vm = length(X)
    N = Int(size_Θ_Vm/2)

    Vm   = @view V[1:N]
    Θ    = @view V[N+1:end]
    
    Vc   = Vm .* exp.(im * Θ)
    Sbus = Diagonal(Vc) * conj.(Ybus) * conj.(Vc)

    return Sbus
end


"""
    SbusC(V, Ybus)


Returns complex power at nodes.
"""
function SbusC(V, Ybus)

    return Diagonal(V) * conj.(Ybus) * conj.(V)
end


"""
    gS(V, p)


Returns complex power injections at nodes.
"""
function gS(V, p)

    Cg, Sg, Sd, Ybus  = p
    
    return Diagonal(V) * conj.(Ybus) * conj.(V) + Sd

end

#----------------------------------------

function VrVi_to_V(x, p)
    
    PV_Vr_Idx, PQ_Vr_Idx, PV_Vi_Idx, PQ_Vi_Idx = p
    
    x_Vr_Idx = first(PV_Vr_Idx):last(PQ_Vr_Idx)
    x_Vi_Idx = first(PV_Vi_Idx):last(PQ_Vi_Idx)
    
    x_Vr = x[x_Vr_Idx]
    x_Vi = x[x_Vi_Idx]
    
    return x_Vr .+ (x_Vi * im)
    
end


function V_to_VrVi(x, (PV_Idx, PQ_Idx))
    
    # PV_Idx, PQ_Idx = p
    
    PV_Vr = real.(x[PV_Idx])
    PV_Vi = imag.(x[PV_Idx])

    PQ_Vr = real.(x[PQ_Idx])
    PQ_Vi = imag.(x[PQ_Idx])

    return  [PV_Vr...;PQ_Vr...;PV_Vi...;PQ_Vi...]
    
end


function V_to_VrVi(V)

    return  [[[real(v), imag(v)] for v in V]...]
    
end


function ΘV_to_V(x, p)
    
    PV_Θ_Idx, PQ_Θ_Idx, PV_V_Idx, PQ_V_Idx = p

    PV_Θ = x[PV_Θ_Idx]
    PQ_Θ = x[PQ_Θ_Idx]
    PV_V = x[PV_V_Idx]
    PQ_V = x[PQ_V_Idx]

    return [[PV_V .* exp.(im * PV_Θ)]...;
            [PQ_V .* exp.(im * PQ_Θ)]...]
    
end


function V_to_ΘV(V, (PV_Idx, PQ_Idx) )
    
    # PV_Idx, PQ_Idx = p

    PV_Θ = angle.(V[PV_Idx])
    PV_V = abs.(V[PV_Idx])

    PQ_Θ = angle.(V[PQ_Idx])
    PQ_V = abs.(V[PQ_Idx])

    return  [PV_Θ...;PQ_Θ...;PV_V...;PQ_V...]
    
end


x_to_V(x, x_Θ_Idx, x_v_Idx) =
    x[x_v_Idx] .* exp.(im * x[x_Θ_Idx] )


V_to_x(V) =
    [angle.(V)...;abs.(V)...]


function x_to_V(x)
    
    n_x_Θ = n_x_v = Int64(length(x)/2 )

    dims       = [n_x_Θ, n_x_v]
    offset     = create_offsets(dims; counter=0)
    Idx        = create_idxs(offset, dims)
    x_Θ_Idx    = Idx[1]
    x_v_Idx    = Idx[2]

    return x[x_v_Idx] .* exp.(im * x[x_Θ_Idx] )
    
end


function ΘV_to_VrVi(x, p)

    # (PV_Idx, PQ_Idx), (PV_Θ_Idx, PQ_Θ_Idx, PV_V_Idx, PQ_V_Idx) = p
    
    p2, p1 = p

    V = ΘV_to_V(x, p1)   

    return V_to_VrVi(V, p2)
    
end


function VrVi_to_ΘV(x, p)

    # ((PV_Idx, PQ_Idx)), (PV_Vr_Idx, PQ_Vr_Idx, PV_Vi_Idx, PQ_Vi_Idx) = p

    p2, p1 = p

    V = VrVi_to_V(x, p1)

    return V_to_ΘV(V, p2)
    
end

#----------------------------------------

# function calc_branch_Ybr(y, y_shunt_km, y_shunt_mk, t_km, t_mk )

#     return [(y + 1im*y_shunt_km)*1/(abs(y_ratio))^2  -y*1/conj(y_ratio); -y*1/y_ratio y + 1im*y_shunt_mk]
# end

"""
    PiModel(
        y, y_shunt_km, y_shunt_mk, t_km, t_mk )


Returns a PiModel for a line.
"""
function PiModel(
    y, y_shunt_km, y_shunt_mk, t_km, t_mk )

    m_km = abs(t_km)
    m_mk = abs(t_mk)
    
    Π = zeros(Complex{Float64}, 2, 2)
    Π[1, 1] = m_km^2 * (y + y_shunt_km) 
    Π[1, 2] = -conj(t_km) * t_mk * y 
    Π[2, 1] = -conj(t_mk) * t_km * y
    Π[2, 2] = m_mk^2 * (y + y_shunt_mk)

    return Π
    
end


"""
    Qmax_Qmin_limit_violation(
        genQ, gen_Qmax, gen_Qmin)


Returns lists of `Qmax_limit_violation` and `Qmim_limit_violation`
"""
function Qmax_Qmin_limit_violation(
    genQ, gen_Qmax, gen_Qmin)
    Qmax_limit_violation = genQ .> gen_Qmax
    Qmim_limit_violation = genQ .< gen_Qmin
    return Qmax_limit_violation, Qmim_limit_violation
end


"""
    flat_reals_to_complex(
        P_flat, Q_flat )


Returns a vector of complex values for a list of real part and a list of imaginary part.
"""
function flat_reals_to_complex(
    P_flat, Q_flat )

    return [P + im * Q
            for (P,Q) in
                zip(P_flat, Q_flat)]

end

#-------------------------------------------------------

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247

"""
    get_gens_vd(
        gens_δ, vh, θh;
        gens_nodes_idx)


Returns generators direct axis voltage.
"""
function get_gens_vd(
    gens_δ, vh, θh;
    gens_nodes_idx)

    return vh[gens_nodes_idx] .* sin.(
        gens_δ - θh[gens_nodes_idx])

end


"""
    get_gens_vq(
        gens_δ, vh, θh;
        gens_nodes_idx)


Returns generators quardrature axis voltage.
"""
function get_gens_vq(
    gens_δ, vh, θh;
    gens_nodes_idx)

    return vh[gens_nodes_idx] .* cos.(
        gens_δ - θh[gens_nodes_idx])

end


"""
     get_a_gen_vd(
        gen_δ, gen_vh, gen_θh)


Returns a generator directive axis voltage.
"""
function get_a_gen_vd(
    gen_δ, gen_vh, gen_θh)

    return gen_vh * sin(gen_δ - gen_θh)
end


"""
     get_a_gen_vq(
        gen_δ, gen_vh, gen_θh)


Returns a generator quardrature axis voltage.
"""
function get_a_gen_vq(
    gen_δ, gen_vh, gen_θh )

    return gen_vh * cos(gen_δ - gen_θh)

end


"""
    get_a_gen_ph(
        vd, vq ,id ,iq)


Returns a generator active power.
"""
function get_a_gen_ph(
    vd,vq,id,iq)

    return vd * id + vq * iq

end


"""
    get_a_gen_qh(
        vd, vq ,id ,iq)


Returns a generator reactive power.
"""
function get_a_gen_qh(
    vd, vq ,id ,iq)

    vq * id - vd * iq

end


"""
    center_of_intertia(
        vec_H, vec_ω, vec_δ, ωs )


Returns center of inertial for `δ_coi`, `ω_coi`, `M_T`, `vec_M`.
"""
function center_of_intertia(
    vec_H, vec_ω, vec_δ, ωs )

    vec_M = (2/ ωs) .* vec_H

    M_T = sum(vec_M)

    δ_coi = (1/M_T) .* sum( vec_M .* vec_δ )
    
    ω_coi = (1/M_T) .* sum( vec_M .* vec_ω )

    return (; δ_coi, ω_coi, M_T, vec_M)

    
end

#-------------------------------------------------------

"""
    threshold_limits(x, x_max, x_min)


Returns `x` if it is within its min and max range, else clip it at thresholds..
"""
function threshold_limits(x, x_max, x_min)
    return x > x_max ? x_max : x < x_min  ? x_min : x
end


"""
    no_limit_violation(x, x_max, x_min)


Returns a boolean if `x` is within its min and max range.
"""
function no_limit_violation(x, x_max, x_min)

    return  (x_min < x)  &&  (x < x_max)
end


"""
    limit_violation(x, x_max, x_min)


Returns a boolean if `x` is outside its min and max range.
"""
function limit_violation(x, x_max, x_min)

    return  (x_min > x)  ||  (x > x_max)
end

#----------------------------------------

"""
    polar_to_cartesian( v_θ )


Returns the cartesian complex values for a polar values.
"""
function polar_to_cartesian( v_θ )

    v, θ = v_θ
    
    return [v * cos(θ), v * sin(θ)]
    
end


"""
    cartesian_to_polar(ur_ui)


Returns the polar values for a cartesian complex values.
"""
function cartesian_to_polar(ur_ui)

    ur, ui = ur_ui
    
    return [abs(ur + im * ui),  angle( ur + im * ui )]
    
end

#----------------------------------------

"""
    ur_ui_to_u(ur_ui)


Returns the complex voltage for a list containing real part and imaginary part of a complex voltage.
"""
function ur_ui_to_u(ur_ui)

    return ur_ui[1] + im * ur_ui[2]
    
end


"""
    xr_xi_to_x(ur_ui)


Returns the complex value for a list containing real part and imaginary part of a complex value.
"""
function xr_xi_to_x(ur_ui)

    return ur_ui[1] + im * ur_ui[2]
    
end


"""
    u_from_ur_ui(ur_ui)


Returns the complex voltage for a list containing real part and imaginary part of a complex voltage.
"""
function u_from_ur_ui(ur_ui)

    return ur_ui[1] + im * ur_ui[2]
    
end


"""
    x_from_xr_xi(ur_ui)


Returns the complex value for a list containing real part and imaginary part of a complex value.
"""
function x_from_xr_xi(ur_ui)

    if length(ur_ui) == 1 && typeof(ur_ui[1]) == Float64
        return ur_ui[1] + im * 0.0
        
    elseif length(ur_ui) == 1 && typeof(ur_ui[1]) == ComplexF64
        return  0.0 + im * ur_ui[1]

    else

        return ur_ui[1] + im * ur_ui[2]
    end
    
end


"""
    conj_x_from_xr_xi(ur_ui)


Returns the complex conjugate for a list containing real part and imaginary part of a complex voltage.
"""
function conj_x_from_xr_xi(ur_ui)

    if length(ur_ui) == 1 && typeof(ur_ui[1]) == Float64
        return ur_ui[1] - im * 0.0
        
    elseif length(ur_ui) == 1 && typeof(ur_ui[1]) == ComplexF64
        return  0.0 - im * ur_ui[1]

    else

        return ur_ui[1] - im * ur_ui[2]
    end
    
end

"""
    u_to_ΘV(u)


Returns a flatted concatenated nodes voltage angles and nodes voltage magnitudes for a list of nodes complex voltages.
"""
function u_to_ΘV(u)

    return [ angle.(u)...; abs.(u)... ]
    
end


"""
    u_to_VΘ(u)


Returns the magnitude and angle for a given complex value.
"""
function u_to_VΘ(u)

    return [ abs(u), angle(u) ]
    
end


"""
    VΘ_to_u(VΘ)


Returns a complex value for a given list of angle and magnitude.
"""
function VΘ_to_u(VΘ)

    return VΘ[1] * exp(im * VΘ[2])
    
end


"""
    ur_ui_to_ΘV(ur_ui)


Returns a flatted concatenated nodes voltage angles and nodes voltage magnitudes for a list of list of nodes voltages real and imaginary parts.
"""
function ur_ui_to_ΘV(ur_ui)

    u = ur_ui_to_u.(ur_ui)

    # Θ = angle.(u) 
    # V = abs.(u)

    return [angle.(u)...;abs.(u)...]
    
end

# ------------------------------------------------------

"""
    get_size_Ybus(
        Ybus)


Returns the memory size of `Ybus`.
"""
function get_size_Ybus(Ybus)

    (sp_ybus_I, sp_ybus_J, sp_ybus_nzv ) =
        findnz(Ybus)

    data_memory_size =
        Base.summarysize(sp_ybus_nzv) 

    idx_memory_size =
                ( Base.summarysize( sp_ybus_I ) +
                Base.summarysize( sp_ybus_J ) )

    data_to_idx_ratio =
        data_memory_size / idx_memory_size
    
    return (;data_memory_size ,
            idx_memory_size, data_to_idx_ratio )

end

#----------------------------------------

"""
    get_size_Ynet_wt_nodes_idx_wt_adjacent_nodes(
        Ynet_wt_nodes_idx_wt_adjacent_nodes)


Returns the memory size of `Ynet_wt_nodes_idx_wt_adjacent_nodes`.
"""
function get_size_Ynet_wt_nodes_idx_wt_adjacent_nodes(
    Ynet_wt_nodes_idx_wt_adjacent_nodes)

    (;Ynet,
     nodes_idx_with_adjacent_nodes_idx ) =
         NamedTupleTools.select(
             Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx ))

    data_memory_size =
        Base.summarysize(Ynet)
    
    idx_memory_size =
                Base.summarysize(
                    nodes_idx_with_adjacent_nodes_idx)
    
    data_to_idx_ratio =
        data_memory_size / idx_memory_size
                    
           
    return (;data_memory_size ,
            idx_memory_size, data_to_idx_ratio )
end

#---------------------------------------- 

"""
    get_size_Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes(
        Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes)


Returns the memory size of `Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes`.
"""
function get_size_Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes(
    Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes)

    (;Yπ_net,
     Yshunt,
     nodes_idx_with_adjacent_nodes_idx) =
         NamedTupleTools.select(
            Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
             (:Yπ_net,
              :Yshunt,
              :nodes_idx_with_adjacent_nodes_idx ))

    data_memory_size =
        (Base.summarysize(Yπ_net) +
        Base.summarysize(Yshunt))

    idx_memory_size =
                Base.summarysize(
                    nodes_idx_with_adjacent_nodes_idx )

    data_to_idx_ratio =
        data_memory_size / idx_memory_size
    
    return (; data_memory_size ,
            idx_memory_size, data_to_idx_ratio  )

end


#---------------------------------------- 

"""
    get_namedtuple_of_findnz(
        Ybus)

Returns a namedtuple of rows indices, columns indices and non-zero values of a sparse matrix.
"""
function get_namedtuple_of_findnz(
    Ybus)

    (sp_I, sp_J, sp_nzv) = findnz(Ybus)

    return (; sp_I, sp_J, sp_nzv)

end


"""
     find_V_idx_in_sparse_matrix_IJV(
         row_idx,
         col_idx,
         sparse_row_idxs,
         sparse_col_idxs )


Returns an index in V of a sparse matrix, where `sparse_row_idxs[idx] = row_idx && sparse_col_idxs[idx] = col_idx`. 
"""
function find_V_idx_in_sparse_matrix_IJV(
    row_idx,
    col_idx,
    sparse_row_idxs,
    sparse_col_idxs )

    tup_sparse_row_col =
        [(a_row, a_col)
         for (a_row, a_col) in
             zip( sparse_row_idxs,
                  sparse_col_idxs ) ]

    V_idx_in_sparse_matrix_IJV =
        findfirst((x) ->
            ( x[1] == row_idx && x[2] == col_idx),
            tup_sparse_row_col )

end


"""
    round_up_Ynet(
        Ynet;
        fractional_digits=4)


Returns a rounded version `Ynet` based on `fractional_digits`.
"""
function round_up_Ynet(
    Ynet;
    fractional_digits=4 )

    map((x) -> round.(
        x;
        digits=fractional_digits),
        Ynet)

end


"""
    write_vector_or_matrix_to_tex(
        tuple_julia_object,
        tex_filename)


Write julia collections in `tuple_julia_object` to a tex file.
"""
function write_vector_or_matrix_to_tex(
    tuple_julia_object,
    tex_filename)
    
    names_julia_object =
        propertynames(tuple_julia_object)

    open(tex_filename, "a") do file_handle

        for (name_object, a_julia_object) in
            zip(names_julia_object,
                tuple_julia_object)

            write(file_handle,
                  '\n')

            write(file_handle,
                  "\n $(String(name_object)) = ")

            write(file_handle,
                  '\n')
            
            write(file_handle,
                  latexify(
                      a_julia_object;
                      fmt=FancyNumberFormatter()))
            
            write(file_handle,
                  '\n')
        end

    end

end


"""
     get_df2tex(df)


Converts  a dataframe to a tex table.
"""
function get_df2tex(
    df)

    
    return latexify(df; env = :table, booktabs = true)
    
end


"""
    get_csv2tex(
        csv_file;
        wt_new_header_bool      = false,
        wt_selected_colums_bool = false,
        new_header              = nothing,
        delim              = ',',
        normalizenames     = true,
        selected_colums    = [] )


Converts a csv data to tex table.
"""
function get_csv2tex(
    csv_file;
    wt_new_header_bool      = false,
    wt_selected_colums_bool = false,
    new_header              = nothing,
    delim              = ',',
    normalizenames     = true,
    selected_colums    = [] )

    if wt_new_header_bool == true

        if wt_selected_colums_bool == true
      
            df = CSV.read(
                csv_file,
                DataFrame;
                header = new_header,
                delim  = delim,
                normalizenames = normalizenames,
                select = selected_colums)
            
        else
      
            df = CSV.read(
                csv_file,
                DataFrame;
                header = new_header,
                delim  = delim,
                normalizenames = normalizenames)
            
        end
        
    else

        if wt_selected_colums_bool == true
        
            df = CSV.read(
                csv_file,
                DataFrame;
                delim = delim,
                normalizenames = normalizenames,
                select = selected_colums)
            
        else
            
            df = CSV.read(
                csv_file,
                DataFrame;
                delim = delim,
                normalizenames = normalizenames)            
        end
        
    end    
    
    return latexify(df; env = :table, booktabs = true)
    
end


"""
    get_pf_transformed_idxs(
        ;slack_gens_nodes_idx,
        non_slack_gens_nodes_idx,
        gens_nodes_idx,
        non_gens_nodes_idx,
        gens_with_loc_load_idx,
        all_nodes_idx,

        n2s_slack_gens_idx,
        n2s_non_slack_gens_idx,
        n2s_gens_idx,
        n2s_non_gens_idx,
        n2s_gens_with_loc_load_idxs,
        n2s_all_nodes_idx )


Returns transformed indices of nodes types.
"""
function get_pf_transformed_idxs(
    ;slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_with_loc_load_idx,
    all_nodes_idx,

    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx )

    
    #-------------------------------

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

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]

    transformed_red_P_mismatch_idx =
        setdiff(transformed_all_nodes_idx,
                transformed_slack_gens_nodes_idx)

    transformed_red_Q_mismatch_idx =
        setdiff(transformed_all_nodes_idx,
                transformed_gens_nodes_idx)

    return (;transformed_slack_gens_nodes_idx,
            transformed_non_slack_gens_nodes_idx,
            transformed_gens_nodes_idx,
            transformed_non_gens_nodes_idx,
            transformed_gens_with_loc_load_idx,
            transformed_all_nodes_idx,
            transformed_red_P_mismatch_idx,
            transformed_red_Q_mismatch_idx)
        
end


"""
    disaggregate_sta_pf_keywords_parameter(
        pf_kw_para)


Returns disaggregated constituent of pf_kw_para.
"""
function disaggregate_sta_pf_keywords_parameter(
    pf_kw_para)

    (;loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs ) =
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :pf_kw_gens_vh_slack_θh_para,
              :pf_kw_net_para,
              :pf_kw_var_idxs,
              :pf_kw_PQ_para_idxs,
              :pf_kw_nodes_types_idxs,
              :pf_kw_n2s_idxs ))

    #-------------------------------

    (; slack_gens_vh,
     slack_gens_θh,

     gens_vh,
     non_slack_gens_vh ) =
         NamedTupleTools.select(
             pf_kw_gens_vh_slack_θh_para,
             (:slack_gens_vh,
              :slack_gens_θh,

              :gens_vh,
              :non_slack_gens_vh ))
    
     (;Ynet,
      nodes_idx_with_adjacent_nodes_idx) =
          NamedTupleTools.select(
              pf_kw_net_para,
              (:Ynet,
               :nodes_idx_with_adjacent_nodes_idx))


     (;red_vh_Idxs,
      red_non_slack_gens_θh_idx2Idx,
      red_non_gens_θh_idx2Idx ) =
          NamedTupleTools.select(
              pf_kw_var_idxs,
              (:red_vh_Idxs,
               :red_non_slack_gens_θh_idx2Idx,
               :red_non_gens_θh_idx2Idx ))
    
     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs ) =
          NamedTupleTools.select(
              pf_kw_PQ_para_idxs,
              (:P_gens_sta_para_Idxs,
               :Q_gens_sta_para_Idxs,
               :P_non_gens_sta_para_Idxs,
               :Q_non_gens_sta_para_Idxs,
               :P_g_loc_load_sta_para_Idxs,
               :Q_g_loc_load_sta_para_Idxs ))

     (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx ) =
         NamedTupleTools.select(
             pf_kw_nodes_types_idxs,
             (:slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx ))    

     (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             pf_kw_n2s_idxs,
             (:n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx ))

    (;transformed_slack_gens_nodes_idx,
     transformed_non_slack_gens_nodes_idx,
     transformed_gens_nodes_idx,
     transformed_non_gens_nodes_idx,
     transformed_gens_with_loc_load_idx,
     transformed_all_nodes_idx,
     transformed_red_P_mismatch_idx,
     transformed_red_Q_mismatch_idx) =
         NamedTupleTools.select(
             get_pf_transformed_idxs(
                 ;slack_gens_nodes_idx,
                 non_slack_gens_nodes_idx,
                 gens_nodes_idx,
                 non_gens_nodes_idx,
                 gens_with_loc_load_idx,
                 all_nodes_idx,

                 n2s_slack_gens_idx,
                 n2s_non_slack_gens_idx,
                 n2s_gens_idx,
                 n2s_non_gens_idx,
                 n2s_gens_with_loc_load_idxs,
                 n2s_all_nodes_idx ),
             (:transformed_slack_gens_nodes_idx,
              :transformed_non_slack_gens_nodes_idx,
              :transformed_gens_nodes_idx,
              :transformed_non_gens_nodes_idx,
              :transformed_gens_with_loc_load_idx,
              :transformed_all_nodes_idx,
              :transformed_red_P_mismatch_idx,
              :transformed_red_Q_mismatch_idx))

    return (;loc_load_exist,
            slack_gens_vh,
            slack_gens_θh,
            gens_vh,
            non_slack_gens_vh,

            Ynet,
            nodes_idx_with_adjacent_nodes_idx,

            red_vh_Idxs,
            red_non_slack_gens_θh_idx2Idx,
            red_non_gens_θh_idx2Idx,

            P_gens_sta_para_Idxs,
            Q_gens_sta_para_Idxs,
            P_non_gens_sta_para_Idxs,
            Q_non_gens_sta_para_Idxs,
            P_g_loc_load_sta_para_Idxs,
            Q_g_loc_load_sta_para_Idxs,

            slack_gens_nodes_idx,
            non_slack_gens_nodes_idx,
            gens_nodes_idx,
            non_gens_nodes_idx,
            gens_with_loc_load_idx,
            all_nodes_idx,

            n2s_slack_gens_idx,
            n2s_non_slack_gens_idx,
            n2s_gens_idx,
            n2s_non_gens_idx,
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,

            transformed_slack_gens_nodes_idx,
            transformed_non_slack_gens_nodes_idx,
            transformed_gens_nodes_idx,
            transformed_non_gens_nodes_idx,
            transformed_gens_with_loc_load_idx,
            transformed_all_nodes_idx,

            transformed_red_P_mismatch_idx,
            transformed_red_Q_mismatch_idx)
    

end


"""
    get_only_adjacent_nodes_idx(
        nodes_idx_with_adjacent_nodes_idx)


Returns only adjacent nodes indicies from nodes_idx_with_adjacent_nodes_idx.

It does not consider the first entry in each row of nodes_idx_with_adjacent_nodes_idx.
"""
function get_only_adjacent_nodes_idx(
    nodes_idx_with_adjacent_nodes_idx)

    return [[nth_node_idx_or_adjacent_node_idx
          for (idx, nth_node_idx_or_adjacent_node_idx) in
              enumerate(
                  nth_node_idx_with_adjacent_nodes_idx)
              if idx !=1 ]
         for nth_node_idx_with_adjacent_nodes_idx in
             nodes_idx_with_adjacent_nodes_idx]
end


"""
    get_nt_components_data( components_data)


Swap the diagonal element of a 2 x 2 matrix.


It is used in `get_Yπ_net` to get an appropriate orientation for elementary admittance matrix `yπ`.
"""
function swap_yπ_diagonal_elements( yπ )

    a11,a21,a12,a22 =  yπ
    
    return [a22 a12; a21 a11 ]
    
end

#-------------------------------------------------------
"""
    get_nt_additional_data(
        additional_data)


Returns additional data as namedtuple.
"""
function get_nt_additional_data(
    additional_data)
    
    return namedtuple(additional_data)
end


"""
    get_nt_components_data(
        components_data)


Returns namedtuples of components names and parameters or data.
"""
function get_nt_components_data(
    components_data)
    
    dict_keys = Tuple(collect(
        keys(components_data)))
    
    dict_values =
        Tuple([namedtuple(a_dict_value)
               for a_dict_value in
                   Tuple(collect(
                       values(components_data)))])
    
    return namedtuple( dict_keys)(
        dict_values)

end


"""
    get_nt_components_type(
        components_type)


Returns namedtuple of components property names and symbols for a component type.
"""
function get_nt_components_type(
    components_type)

    nt_components_type =
        namedtuple(
            components_type)
    
    properties = [ ]

    property_names =
        propertynames(
            nt_components_type)
    
    for a_property_name in
        property_names
        
        a_property = getproperty(
            nt_components_type,
            a_property_name)

        push!(properties,
              a_property)
    end
        
    return namedtuple(Tuple(property_names))(
        Tuple(Symbol.(properties)))

end

#---------------------------------------------------


"""
    recursive_dict_to_namedtuple(
        mapping)


Converts recursively a nested dictionary to a nested namedtuple.

# Example

```Julia

dict_values =
    Dict{Symbol, Any}(
        :V_R_max => 5.2,
        :Tf => 0.35,
        :Be => 1.555,
        :Tr => 0.001,
        :Ta => 0.2,
        :Ka => 20,
        :Te => 0.314,
        :Kf => 0.063,
        :V_R_min => -4.16,
        :Ke => 1,
        :gov => Dict{Symbol, Any}(
            :p_max => 2.2,
            :R => 0.02,
            :p_min => 0.1,
            :Ts => 0.2,
            :Tc => 0.4,
            :avr => Dict{Symbol, Any}(
            :Tr => 0.001,
            :Ta => 0.2,
            :Ka => 20,
            :Te => 0.314) ))

recursive_dict_to_namedtuple(dict_values)

```

"""
function recursive_dict_to_namedtuple(
    mapping)
    
    if isa(mapping, Dict) 
        
        for (key, value) in mapping

            mapping[key] =
                recursive_dict_to_namedtuple(
                    value)
            
        end

        return namedtuple(mapping)
        
    else
        return mapping
    end
                
end


"""
    recursive_nested_namedtuple_wt_dict(
        nested_namedtuple)


Converts recursively a nested dictionary to a nested namedtuple.

# Example

```Julia


nt_values =
    (
        V_R_max = 5.2,
        Tf = 0.35,
        Be = 1.555,
        Tr = 0.001,
        Ta = 0.2,
        Ka = 20,
        Te = 0.314,
        Kf = 0.063,
        V_R_min = -4.16,
        Ke = 1,
        gov = Dict{Symbol, Any}(
            :p_max => 2.2,
            :R => 0.02,
            :p_min => 0.1,
            :Ts => 0.2,
            :Tc => 0.4,
            :avr => Dict{Symbol, Any}(
            :Tr => 0.001,
            :Ta => 0.2,
            :Ka => 20,
                :Te => 0.314) ))

recursive_nested_namedtuple_wt_dict(
    nt_values)

```
"""
function recursive_nested_namedtuple_wt_dict(
    nested_namedtuple)

    properties = [ ]
    
    property_names =
        propertynames(nested_namedtuple)

    for a_property_name in property_names
        a_property =
            getproperty(
                nested_namedtuple,
                a_property_name)
        
        r_nt_property =
            recursive_dict_to_namedtuple(
                a_property)

        push!(properties,
              r_nt_property )
    end
    
    return namedtuple(
        Tuple(property_names))(
            Tuple(properties))

end

"""
    namedtuple_nested_selection(
        nested_namedtuple;
        sequence_order = (:nothing,),
        selections = (:nothing,) )


Returns a set of selection at the same level in a nested namedtuple.

# Example

```Julia

    some_plant_generators_data_from_json =
        [(idx = 1,
          plant_type = "plant_wt_loc_load_v6",
          components_type =
              (gen = "SM_2axis_wt_loc_load_cb_v6",
               avr = "avr_t1_cb",
               gov = "gov_ieee_tgov1_cb"),
          components_data = (
              gen = (T_d_dash = 7.4, D = 0, X_d = 0.8979,
                     vh = 1.06, X_q = 0.646,
                     ωs = 376.99111843077515,
                     T_d_2dash = 0.03, Qmin = 0,
                     Sn = 332.5503871596002,
                     X_d_dash = 0.2995, T_q_dash = 0.3,
                     Ωb = 376.99111843077515,
                     Bus = "bus1", X_q_dash = 0.646,
                     ra = 0, X_q_2dash = 0.4,
                     Q = -16.9, Pmax = 332.4,
                     Pmin = 0, xℓ = 0.2396,
                     P = 232.4, vmax = 1.06,
                     H = 5.148, Qmax = 10,
                     vmin = 0.94, T_q_2dash = 0.033,
                     X_d_2dash = 0.23), )),
         (idx = 2,
          plant_type = "plant_wt_loc_load_v6",
          components_type =
              (loc_load = "loc_Load_t1",
               gen = "SM_2axis_wt_loc_load_cb_v6",
               avr = "avr_t1_cb",
               gov = "gov_ieee_tgov1_cb"),
          components_data = (
              loc_load = (Bus = "bus2",
                          loc_P = 21.7,
                          loc_Q = 12.7),
              gen = (T_d_dash = 6.1, D = 0, X_d = 1.05,
                     vh = 1.045, X_q = 0.98,
                     ωs = 376.99111843077515,
                     T_d_2dash = 0.04, Qmin = -40,
                     Sn = 148.66068747318505,
                     X_d_dash = 0.185, T_q_dash = 0.3,
                     Ωb = 376.99111843077515,
                     Bus = "bus2", X_q_dash = 0.36,
                     ra = 0.0031, X_q_2dash = 0.13,
                     Q = 42.4, Pmax = 140, Pmin = 0,
                     xℓ = 0, P = 40, vmax = 1.06,
                     H = 6.54, Qmax = 50, vmin = 0.94,
                     T_q_2dash=0.099,X_d_2dash=0.13),))]


    sequence_order = (:components_data, :gen)
    
    selections = (:P, :Q)
    
    namedtuple_nested_selection(
        some_plant_generators_data_from_json;
        sequence_order = sequence_order,
        selections = selections )

2-element Vector{Any}:
 (P = 232.4, Q = -16.9)
 (P = 40, Q = 42.4)

```
"""
function namedtuple_nested_selection(
    nested_namedtuple;
    sequence_order = (:nothing,),
    selections = (:nothing,) )

    selected_properties = []
    
    sequence_order_depth =
        length(sequence_order)

    k_stage_selection = nested_namedtuple
    
    for (k, a_sequence_order) in
        enumerate(sequence_order)
        
        # global k_stage_selection
        
        k_stage_selection =
            NamedTupleTools.select.(
                k_stage_selection,
                a_sequence_order)        

        if k == sequence_order_depth
            
            for nt_item in k_stage_selection
                
                Set_property_names =
                    Set(collect(propertynames(nt_item)))
            
                issubset(Set(collect(selections)),
                         Set_property_names) ?
                    push!(selected_properties,
                      NamedTupleTools.select(
                      nt_item,
                      selections ) ) : nothing
            end
            
        end

    end
    
    return selected_properties

end

#---------------------------------------------------

"""
    get_selected_vec_nt_to_vec_vec(
        vec_namedtuple;
        selections = (:nothing, ) )


Convert `vec_namedtuple` to nt_vector, vector of namedtuples
to namedtuples of vector.

"""
function get_selected_vec_nt_to_vec_vec(
    vec_namedtuple;
    selections = (:nothing, ) )

    dim_selections =
        length( selections )
    
    vec_vec = Vector{Union{Float64,ComplexF64}}[
        [] for a_para in 1:dim_selections ]
    
    vec_selected_nt =
        [ NamedTupleTools.select(a_nt, selections)
          for a_nt in
              vec_namedtuple]

    for (idx, a_property) in enumerate(selections)
        for a_namedtuple in vec_selected_nt
            
            # push!(vec_vec[idx],
            #       getproperty(a_namedtuple,
            #                   a_property))

            
            if a_property ∈ (:Sn, :vh, :P, :Q, :Pmin, :Pmax,
                             :Qmin, :Qmax, :vmin, :vmax)

                a_property_value_or_nothing =
                    getproperty(
                          a_namedtuple,
                        a_property)

                a_property_value =
                    a_property_value_or_nothing == nothing ?
                    99999 : a_property_value_or_nothing

                push!(vec_vec[idx], a_property_value )
                
            else

                push!(vec_vec[idx],
                      getproperty(
                          a_namedtuple,
                          a_property) )
                
            end

        end
    end
    
    return namedtuple(
        OrderedDict(a_sym => a_value
             for (a_sym, a_value) in
                 zip(selections, vec_vec)))

end


"""
    get_selected_vec_nt_to_vec_vec(
        vec_namedtuple,
        nothing;
        selections = (:nothing, ),
        vec_datatype = Float64 )


Convert
`vec_namedtuple` to nt_vector, vector of namedtuples to namedtuples of vector.

"""
function get_selected_vec_nt_to_vec_vec(
    vec_namedtuple,
    nothing;
    selections = (:nothing, ),
    vec_datatype = Float64 )

    dim_selections =
        length( selections )
    
    vec_vec = Vector{vec_datatype}[
        [] for a_para in 1:dim_selections ]
    
    vec_selected_nt =
        [ NamedTupleTools.select(a_nt, selections)
          for a_nt in
              vec_namedtuple ]

    for (idx, a_property) in enumerate(selections)
        for a_namedtuple in vec_selected_nt
            
            # push!(vec_vec[idx],
            #       getproperty(a_namedtuple,
            #                   a_property))
            
            if a_property ∈ (:Sn, :vh, :P, :Q, :Pmin, :Pmax,
                             :Qmin, :Qmax, :vmin, :vmax)

                a_property_value_or_nothing =
                    getproperty(
                          a_namedtuple,
                        a_property)

                a_property_value =
                    a_property_value_or_nothing == nothing ?
                    99999 : a_property_value_or_nothing

                push!(vec_vec[idx], a_property_value )
                
            else

                push!(vec_vec[idx],
                      getproperty(
                          a_namedtuple,
                          a_property) )
                

            end
            
        end
    end
    
    return namedtuple(
        OrderedDict(a_sym => a_value
             for (a_sym, a_value) in
                 zip(selections, vec_vec)))

end


"""
    get_nt_vec_wt_vec_vec_per_paras(
        vec_vec_per_node;
        nt_syms = (:nothing, ),
        vec_datatype = Float64 )


Convert `vec_namedtuple` to nt_vector, vector of namedtuples to namedtuples of vector per parameters.

"""
function get_nt_vec_wt_vec_vec_per_paras(
    vec_vec_per_node;
    nt_syms = (:nothing, ),
    vec_datatype = Float64 )

    dim_selections =
        length( nt_syms  )

    @assert dim_selections ==
        length( vec_vec_per_node[1])
    
    vec_vec = Vector{vec_datatype}[
        [] for a_para in 1:dim_selections ]
    
    for (idx, a_sym) in enumerate(nt_syms)
        for a_vec in vec_vec_per_node
            push!(vec_vec[idx],
                  a_vec[idx] )
        end
    end

    nt_vec_per_paras = namedtuple(
        OrderedDict(a_sym => a_value
             for (a_sym, a_value) in
                    zip(nt_syms, vec_vec)))
    
    vec_vec_per_paras = vec_vec
    return (; nt_vec_per_paras,
            vec_vec_per_paras)

end


#---------------------------------------------------


"""
    get_dict_nt_params_from_json_lib_file(
        json_nt_params_libs_file )


Returns a dictionary of namedtuples of components parameters from a json library file
"""
function get_dict_nt_params_from_json_lib_file(
    json_nt_params_libs_file )

    #--------------------------------------

    nt_params_from_file =
        JSON3.read( json_nt_params_libs_file )
    
    #--------------------------------------

    return Dict{Symbol, NamedTuple}(
        a_key => namedtuple(a_dict)
        for ( a_key, a_dict ) in
            convert(Dict{Symbol, Dict{Symbol,Float64}},
                    nt_params_from_file ))
    
end

#---------------------------------------------------

"""
    get_nested_nt_from_nt_wt_dict(
        a_nested_nt )


Converts a dict in a namedtuple to a namedtuple.

# Example

```Julia

tt = (plant_type = "plant_cb_v6",
      components_type = Dict{Symbol, Any}(
          :gen => "SM_2axis_cb_v6",
          :avr => "avr_t1_cb_sauer",
          :gov => "gov_t1_cb_sauer"),
      idx = 1)

get_a_plant_data_json_to_nt(
    json_a_plant;
    in_components_type_sym = true)

get_nested_nt_from_nt_wt_dict( tt )

```

"""
function get_nested_nt_from_nt_wt_dict(
    a_nested_nt )
    
    properties = [ ]

    property_names = propertynames( a_nested_nt )
    
    for a_property_name in property_names
        
        a_property = getproperty(
            a_nested_nt, a_property_name)
        (typeof(a_property) ∈
            (Int64,Float64,Float32,String,Symbol) ||
            isa(a_property,Number)) ? push!(
                properties, a_property) :  push!(
                    properties, namedtuple(a_property))
            
        
    end

    return namedtuple(Tuple(property_names))(
        Tuple(properties))

    
end

#---------------------------------------------------

"""
    get_dict_struct_name_type(
        list_struct_sym_types)

Returns a mapping of data structure symbolic names to data structures. 
"""
function get_dict_struct_name_type(
    list_struct_sym_types)

    return  Dict{Symbol, DataType}(
        nameof(a_type) => a_type for a_type in
            list_struct_sym_types )
    
end

#---------------------------------------------------


"""
    get_abstract_type_dict_subtypes(
        absract_type )


Returns a dictionary of concrete subtypes of absract_type.

The set of abstract types currently defined in
the package are:

SdAvr, SdBranchElement, SdGen, SdGenPlant,
SdGov, SdNonGen, SdNonGenPlant, SdPss.

They are subtypes of

`AbstractPowerSystemComponent`

"""
function get_abstract_type_dict_subtypes(
    absract_type )
    
    return get_dict_struct_name_type(
        subtypes(absract_type))    

end


"""
    get_absract_type_dict_subsubtypes(
        absract_type )


Returns a dictionary of symbol to subsubtypes of
an abstract types, e,g. `AbstractPowerSystemComponent`.
This is used to translate a symbol of a type to type

# Example

```Julia
get_absract_type_dict_subsubtypes(
    absract_type )
```

dict_symbol_types = get_absract_type_dict_subsubtypes(
    AbstractPowerSystemComponent )

`dict_symbol_types[:pss_t2_cb]` will return `pss_t2_cb`

"""
function get_absract_type_dict_subsubtypes(
    absract_type )
    
    # subsubtype =
    #     get_subsubtype( absract_type )
    
    # return get_dict_struct_name_type(
    #     subsubtype )    
            
    return get_dict_struct_name_type(
        get_subsubtype( absract_type ) )    

end


#---------------------------------------------------
#---------------------------------------------------
# Dataframe data type conversion
#---------------------------------------------------
#---------------------------------------------------

"""
    convert_dataframe_selected_cols_types(
        df, cols_types, cols_names )


Converts a selected columns to specific types in a dataframe.

The function was created because of the challenges of
getting `is_slack` column as a Bool in dyn_plant.csv
"""
function convert_dataframe_selected_cols_types(
    df, cols_types, cols_names )
    for (a_type, a_name) in zip(cols_types,cols_names)
        if a_type == Symbol
            
            df[!, a_name] = Symbol.(df[!, a_name])
            
        elseif a_type == String
            
            df[!, a_name] = String.(df[!, a_name])
                        
        elseif a_type == Bool

            if eltype(df[!,a_name]) == String
                df[!, a_name] = convert(
                    Vector{Bool},
                    strip.(lowercase.(
                        df[!,a_name])) .== "true")
            end
            
            
        else

            df[!, a_name] = convert.(a_type, df[!, a_name])

            
        end
        
    end
    
    return df
end


#---------------------------------------------------
# dir and files
#---------------------------------------------------


"""
    get_sub_components_libs_dir(
        components_libs_dir,
        sub_components_strings )


Returns folders where sub-componets are stored.
"""
function get_sub_components_libs_dir(
    components_libs_dir,
    sub_components_strings )

    sub_components_dir =
        [ joinpath(components_libs_dir, a_sub_comp)
         for a_sub_comp in
             sub_components_strings ]

    return Tuple(sub_components_dir)

end


"""
    get_sub_components_libs_files(
        components_libs_dir,
        components_files_string;
        ext = "json" )

Returns folders where sub-componets are stored.
"""
function get_sub_components_libs_files(
    components_libs_dir,
    components_files_string;
    ext = "json" )

    sub_components_dir_string =
        [ strip(split(json_comp, "-")[1])
          for json_comp in
              components_files_string ]

    sub_components_libs_dir =
        get_sub_components_libs_dir(
            components_libs_dir,
            sub_components_dir_string )

    sub_components_libs_files =
        [ joinpath(a_sub_comp_dir,
                   "$(a_sub_comp_string).$(ext)")
          for (a_sub_comp_dir, a_sub_comp_string) in
              zip(sub_components_libs_dir,
                  components_files_string ) ]

    return Tuple(sub_components_libs_files)

end


#---------------------------------------------------


"""
    get_net_nodes_type_idxs_by_json(
        plant_generators_data_from_json,
        plant_loads_data_from_json,
        plant_transmission_data_from_json )

Returns list of indices of various type of nodes in a network.

These are:

`slack_bus_idx,
 gens_idx,
 slack_gens_nodes_idx,
 non_slack_gens_nodes_idx,
 gens_nodes_idx,

 gens_with_loc_load_idx,
 gens_nodes_with_loc_loads_idx,

 loc_load_exist,
 load_nodes_idx,
 transmission_nodes_idx,
 non_gens_nodes_idx,
 all_nodes_idx,
 non_slack_gens_and_non_gens_idx,
 nodes_with_demands_idx`

"""
function get_net_nodes_type_idxs_by_json(
    plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json )

    (slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_with_loc_load_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist) =
         get_gens_plants_and_loc_loads_idx_by_json(
             plant_generators_data_from_json )

    (load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx ) =
         get_non_gen_plants_idx_by_json(
             plant_loads_data_from_json,
             plant_transmission_data_from_json )


    non_slack_gens_and_non_gens_idx =
        sort([non_slack_gens_nodes_idx;
              non_gens_nodes_idx])
    
    nodes_with_demands_idx =
        convert(Vector{Int64}, sort([load_nodes_idx;
              gens_with_loc_load_idx]))
    
    all_nodes_idx =
        sort([gens_nodes_idx;
              non_gens_nodes_idx])

    
    return (;slack_bus_idx,
            gens_idx,
            slack_gens_nodes_idx,
            non_slack_gens_nodes_idx,
            gens_nodes_idx,
            
            gens_with_loc_load_idx,
            gens_nodes_with_loc_loads_idx,
            
            loc_load_exist,
            load_nodes_idx,
            transmission_nodes_idx,
            non_gens_nodes_idx,
            all_nodes_idx,
            non_slack_gens_and_non_gens_idx,
            nodes_with_demands_idx)

end


#---------------------------------------------------
#---------------------------------------------------

"""
    get_dict_net_streamlined_idx_by_nodes_type_idxs(
        net_nodes_type_idxs)


Returns "indices to ordinal" dictionaries of indices of various type of nodes in a network.

These are:

`n2s_slack_gens_idx,
 n2s_non_slack_gens_idx,
 n2s_gens_idx,
 n2s_non_gens_idx,
 n2s_load_idx,
 n2s_gens_with_loc_load_idxs,
 n2s_transmission_idxs,
 n2s_all_nodes_idx,
 n2s_nodes_with_demands_idx`

It aligns indices of node types with the indices of parameters
of node types.

Consider a 14 nodes network, where 6 nodes are generator
nodes, 1 node is a transmission node, and 7 nodes are load
nodes.

The indices of generator nodes in the network are
gen_nodes_idx = [1,2,3,6,8]. The size of `gen_nodes_idx` is 5.

Lets define an array for reactive power,
Qg = [-0.169, 0.42, 0.234, 0.122, 0.174]

How do we access the reactive power of the fifth generator,
whose index in the network is 8? It is evident that Qg[8]
will throw an error.

A way out is to use a dictionary

n2s_type_idx = OrderedDict(
    idx => ord
    for (idx, ord) in
        zip(gen_nodes_idx,
            collect(1:length(gen_nodes_idx)) ))

i.e n2s_type_idx = Dict(1=>1, 2=>2, 3=3, 6=>4, 8=5)

The reactive power of the fifth generator can be subsequently
accessed as :

 Qg[ n2s_type_idx[ 8 ] ]

Note n2s_type_idx[ 8 ] will produced 5.

The beauty of this is that, indices of nodes are not resticted to
numbers. Strings, symbols can be used as indices.

If gen_nodes_idx had been,

gen_nodes_idx =
    ["node-gauteng","node-limpopo","node-wc","node-ec", "node-nw"]


n2s_type_idx = OrderedDict(
    idx => ord
    for (idx, ord) in
        zip(gen_nodes_idx,
            collect(1:length(gen_nodes_idx)) ))

Qg[ n2s_type_idx[ "node-nw" ] ] will produce, the reactive power
of the fifth generator,

A generic function could be defined

function get_n2s_any( a_net_group_idxs)

    return OrderedDict{Union{Symbol,String,Int64},Int64}(
            net_idx =>idx
            for (net_idx, idx) in zip(
                a_net_group_idxs,
                collect(1:length(
                    a_net_group_idxs )) ) )
    
end

"""
function get_dict_net_streamlined_idx_by_nodes_type_idxs(
    net_nodes_type_idxs)

    (;slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     nodes_with_demands_idx)  =
         NamedTupleTools.select(
             net_nodes_type_idxs,
             (:slack_bus_idx,
              :gens_idx,
              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :loc_load_exist,
              :load_nodes_idx,
              :transmission_nodes_idx,
              :non_gens_nodes_idx,
              :all_nodes_idx,
              :non_slack_gens_and_non_gens_idx,
              :nodes_with_demands_idx))
             
    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx

    #------------------------------------------

    streamlined_idx = [ slack_gens_nodes_idx,
                        non_slack_gens_nodes_idx,
                        gens_nodes_idx,
                        non_gens_nodes_idx,
                        load_nodes_idx,
                        all_nodes_idx,
                        nodes_with_demands_idx]
    
    # vec_dict_net_to_streamlined_idx =
    #     [ OrderedDict{Union{Symbol,String,Int64},Int64}(
    #         net_idx =>idx
    #         for (net_idx, idx) in
    #             zip( a_net_group_idxs,
    #                  collect(1:length(
    #                      a_net_group_idxs ))))
    #   for a_net_group_idxs in
    #       streamlined_idx ]

    
    vec_dict_net_to_streamlined_idx =
        [ get_n2s_any( a_net_group_idxs)
          for a_net_group_idxs in
              streamlined_idx ]
    
    (n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_all_nodes_idx,
     n2s_nodes_with_demands_idx) =
         vec_dict_net_to_streamlined_idx

    """
    gens_with_loc_load and transmission nodes
    need special attention in cases
    they do not exist

    """        
    if loc_load_exist == true
        n2s_gens_with_loc_load_idxs =
            get_n2s_any( gens_with_loc_load_idx)

        # OrderedDict{Int64, Int64}(
        #     net_idx =>idx
        #     for (net_idx, idx) in zip(
        #         gens_with_loc_load_idx,
        #         collect(1:length(
        #             gens_with_loc_load_idx ))))
        
    else
        
        n2s_gens_with_loc_load_idxs =
            nothing
    end

        
    if transmission_nodes_idx != []
        
        n2s_transmission_idxs =
            get_n2s_any( transmission_nodes_idx)
        
        # OrderedDict{Int64, Int64}(
        #     net_idx =>idx
        #     for (net_idx, idx) in zip(
        #         transmission_nodes_idx,
        #         collect(1:length(
        #             transmission_nodes_idx ))))
        
    else
        
        n2s_transmission_idxs =
            nothing
        
    end

        return (;n2s_slack_gens_idx,
                n2s_non_slack_gens_idx,
                n2s_gens_idx,
                n2s_non_gens_idx,
                n2s_load_idx,
                n2s_gens_with_loc_load_idxs,
                n2s_transmission_idxs,
                n2s_all_nodes_idx,
                n2s_nodes_with_demands_idx)
    
end


#---------------------------------------------------

# function get_all_nodes_idx_by_mpc(
#     mpc_bus)

#     return sort(copy(mpc_bus.bus_i))
# end

"""
    get_a_n2s_net_group(
        a_net_group_idxs;
        loc_load_exist = false,
        transmission_group = false )



Returns  streamlined idxs conversion dictionary

    n2s_gens_idx =
        get_a_n2s_net_group(gens_nodes_idx)

    n2s_transmission_idxs =
        get_a_n2s_net_group(transmission_nodes_idx;
        transmission_group = true)

    n2s_gens_with_loc_load_idxs =
        get_a_n2s_net_group(gens_with_loc_load_idx;
        loc_load_exist = true)

"""
function get_a_n2s_net_group(
    a_net_group_idxs;
    loc_load_exist = false,
    transmission_group = false )

    if (loc_load_exist == false &&
        transmission_group == false)
        
        return get_n2s_any( a_net_group_idxs)
        
    elseif loc_load_exist != false
        
        return get_n2s_any( a_net_group_idxs)
        
    elseif transmission_group != false
        
        return get_n2s_any( a_net_group_idxs)
        
    else
        return nothing
    end
    
end

#---------------------------------------------------

"""
    get_net_nodes_type_idxs_by_mpc( mpc_bus )


Returns list of indices of various type of nodes in a network.

"""
function get_net_nodes_type_idxs_by_mpc(
    mpc_bus )


    slack_bus_idx = slack_gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 3 ]


    non_slack_gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 2 ]

    gens_idx = gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 3 || node_type == 2  ]


    non_gens_nodes_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 1 ]

    non_slack_gens_and_non_gens_idx =
        [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 2 || node_type == 1  ]

    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx =
        [a_node for (a_node, node_type, node_Pd, node_Qd) in
             zip( mpc_bus.bus_i, mpc_bus.type,
                  mpc_bus.Pd, mpc_bus.Qd )
             if (node_type == 3 || node_type == 2) && (
                 (node_Pd != 0.0 || node_Pd != 0) ||
                 (node_Qd != 0.0 || node_Qd != 0) ) ]


    nodes_with_demands_idx =
        [a_node for (a_node, node_Pd, node_Qd) in
             zip( mpc_bus.bus_i, mpc_bus.Pd, mpc_bus.Qd)
             if  (node_Pd != 0.0 || node_Pd != 0) ||
                 (node_Qd != 0.0 || node_Qd != 0)]


    load_nodes_idx =
        [a_node for (a_node, node_type, node_Pd, node_Qd) in
             zip(mpc_bus.bus_i, mpc_bus.type, mpc_bus.Pd,
                 mpc_bus.Qd)
             if (node_type == 1) && (( (node_Pd != 0.0 || node_Pd != 0) || (node_Qd != 0.0 || node_Qd != 0))) ]
    

    transmission_nodes_idx =
        [a_node for (a_node,a_type, node_Pd, node_Qd) in
             zip(mpc_bus.bus_i, mpc_bus.type, mpc_bus.Pd,
                 mpc_bus.Qd)
             if (a_type == 1) &&
                 (( (node_Pd == 0.0 || node_Pd == 0) &&
                 (node_Qd == 0.0 || node_Qd == 0))) ]
    
    all_nodes_idx = copy(mpc_bus.bus_i)

    loc_load_exist =
        length(gens_nodes_with_loc_loads_idx) != 0 ?
        true : false

    return (; slack_bus_idx,
            gens_idx,
            slack_gens_nodes_idx,
            non_slack_gens_nodes_idx,
            gens_nodes_idx,
            gens_with_loc_load_idx,
            gens_nodes_with_loc_loads_idx,
            loc_load_exist,
            load_nodes_idx,
            transmission_nodes_idx,
            non_gens_nodes_idx,
            all_nodes_idx,
            non_slack_gens_and_non_gens_idx,
            nodes_with_demands_idx)
    
end


"""
    get_dict_net_to_streamlined_idx_by_mpc(mpc_bus)


Returns "indices to ordinal" dictionaries of indices of various type of nodes in a network.

"""
function get_dict_net_to_streamlined_idx_by_mpc(
    mpc_bus)

    (;slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_with_loc_load_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     nodes_with_demands_idx )  =
         NamedTupleTools.select(
             get_net_nodes_type_idxs_by_mpc(
                 mpc_bus  ),
             (:slack_bus_idx,
              :gens_idx,
              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :gens_with_loc_load_idx,
              :gens_nodes_with_loc_loads_idx,
              :loc_load_exist,
              :load_nodes_idx,
              :transmission_nodes_idx,
              :non_gens_nodes_idx,
              :all_nodes_idx,
              :non_slack_gens_and_non_gens_idx,
              :nodes_with_demands_idx ))

    
    # gens_with_loc_load_idx =
    #     gens_nodes_with_loc_loads_idx

    #------------------------------------------

    streamlined_idx = [ slack_gens_nodes_idx,
                        non_slack_gens_nodes_idx,
                        gens_nodes_idx,
                        non_gens_nodes_idx,
                        load_nodes_idx,
                        all_nodes_idx,

                        nodes_with_demands_idx]

    vec_dict_net_to_streamlined_idx =
        [ OrderedDict{Int64,Int64}(
            net_idx =>idx
            for (net_idx, idx) in
                zip( a_net_group_idxs,
                     collect(1:length(
                         a_net_group_idxs ))))
      for a_net_group_idxs in
          streamlined_idx ]

    (dict_n2s_slack_gens_idx,
     dict_n2s_non_slack_gens_idx,
     dict_n2s_gens_idx,
     dict_n2s_non_gens_idx,
     dict_n2s_load_idx,
     dict_n2s_all_nodes_idx,
     
     dict_n2s_nodes_with_demands_idx) =
         vec_dict_net_to_streamlined_idx

    """
    gens_with_loc_load and transmission nodes
    need special attention in cases
    they do not exist

    """        
    if loc_load_exist == true
        dict_n2s_gens_with_loc_load_idxs =
        OrderedDict{Int64, Int64}(
            net_idx =>idx
            for (net_idx, idx) in zip(
                gens_with_loc_load_idx,
                collect(1:length(
                    gens_with_loc_load_idx ))))
    else
        dict_n2s_gens_with_loc_load_idxs = nothing
    end

        
    if transmission_nodes_idx != []
        
        dict_n2s_transmission_idxs =
        OrderedDict{Int64, Int64}(
            net_idx =>idx
            for (net_idx, idx) in zip(
                transmission_nodes_idx,
                collect(1:length(
                    transmission_nodes_idx ))))
    else
        dict_n2s_transmission_idxs = nothing
    end

        return (;dict_n2s_slack_gens_idx,
                dict_n2s_non_slack_gens_idx,
                dict_n2s_gens_idx,
                dict_n2s_non_gens_idx,
                dict_n2s_load_idx,
                dict_n2s_gens_with_loc_load_idxs,
                dict_n2s_transmission_idxs,
                dict_n2s_all_nodes_idx,

                dict_n2s_nodes_with_demands_idx)
    
end


"""
    get_dict_n2s_streamlined_idx_by_mpc(mpc_bus)


Returns "indices to ordinal" dictionaries of indices of various type of nodes in a network.

"""
function get_dict_n2s_streamlined_idx_by_mpc(
    mpc_bus)

    (;
     dict_n2s_slack_gens_idx,
     dict_n2s_non_slack_gens_idx,
     dict_n2s_gens_idx,
     dict_n2s_non_gens_idx,
     dict_n2s_load_idx,
     dict_n2s_gens_with_loc_load_idxs,
     dict_n2s_transmission_idxs,
     dict_n2s_all_nodes_idx,
     dict_n2s_nodes_with_demands_idx) =
         NamedTupleTools.select(
             get_dict_net_to_streamlined_idx_by_mpc(
                 mpc_bus),
             (:dict_n2s_slack_gens_idx,
              :dict_n2s_non_slack_gens_idx,
              :dict_n2s_gens_idx,
              :dict_n2s_non_gens_idx,
              :dict_n2s_load_idx,
              :dict_n2s_gens_with_loc_load_idxs,
              :dict_n2s_transmission_idxs,
              :dict_n2s_all_nodes_idx,
              :dict_n2s_nodes_with_demands_idx))


    n2s_slack_gens_idx  =
        dict_n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx  = dict_n2s_gens_idx
    
    n2s_non_gens_idx = dict_n2s_non_gens_idx
    
    n2s_load_idx  = dict_n2s_load_idx
    
    n2s_gens_with_loc_load_idxs  =
        dict_n2s_gens_with_loc_load_idxs
    
    n2s_transmission_idxs =
        dict_n2s_transmission_idxs
    
    n2s_all_nodes_idx =
        dict_n2s_all_nodes_idx

    n2s_nodes_with_demands_idx =
        dict_n2s_nodes_with_demands_idx
    
    
    return (;
            n2s_slack_gens_idx,         
            n2s_non_slack_gens_idx,     
            n2s_gens_idx,               
            n2s_non_gens_idx,           
            n2s_load_idx,               
            n2s_gens_with_loc_load_idxs,
            n2s_transmission_idxs,      
            n2s_all_nodes_idx,
            n2s_nodes_with_demands_idx)          
    
end

#---------------------------------------------------

"""
    get_similar_collection_diff(
        collection_1, collection_2)


Returns the difference between two list, tuple, collections.
"""
function get_similar_collection_diff(
    collection_1, collection_2)

    if typeof(collection_1) == Vector{Tuple{Int64, Int64}}

        return [collect(a) .- collect(b) for (a,b) in
                    zip(collection_1,
                        collection_2)]
        
    else

        return [a .- b for (a,b) in
                    zip(collection_1,
                        collection_2)]
        
    end

end


"""
    get_Ynet_wt_nodes_idx_wt_adjacent_nodes_diff(
        Ynet_wt_nodes_idx_wt_adjacent_nodes_1,
        Ynet_wt_nodes_idx_wt_adjacent_nodes_2)

Returns the difference between two Ynet_wt_nodes_idx_wt_adjacent_nodes data structure.
"""
function get_Ynet_wt_nodes_idx_wt_adjacent_nodes_diff(
    Ynet_wt_nodes_idx_wt_adjacent_nodes_1,
    Ynet_wt_nodes_idx_wt_adjacent_nodes_2)

    (Ynet_1,
     nodes_idx_with_adjacent_nodes_idx_1 ) =
         NamedTupleTools.select(
             Ynet_wt_nodes_idx_wt_adjacent_nodes_1,
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx))

    (Ynet_2,
     nodes_idx_with_adjacent_nodes_idx_2 ) =
         NamedTupleTools.select(
             Ynet_wt_nodes_idx_wt_adjacent_nodes_2,
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx))

    diff_Ynet = [a .- b for (a,b) in
                    zip(Ynet_1,
                        Ynet_2)]

    diff_nodes_idx_with_adjacent_nodes_idx =
        [a .- b for (a,b) in
                    zip(nodes_idx_with_adjacent_nodes_idx_1,
                       nodes_idx_with_adjacent_nodes_idx_2)]

    return (;diff_Ynet,
            diff_nodes_idx_with_adjacent_nodes_idx)
end

#---------------------------------------------------
#---------------------------------------------------


"""
    get_a_node_state_algb_vars_indices_in_syms(
        ; node_syms_labels = node_syms_labels,
        bus_name = bus_name, vars = [:ω, :δ])


Returns indices of state_algebraic variables of a node in the node syms list.

It is meant to extract the indices in sol.
"""
function get_a_node_state_algb_vars_indices_in_syms(
    ; node_syms_labels = node_syms_labels,
    bus_name = bus_name, vars = [:ω, :δ])

    vars_indices = []

    for a_var in vars

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   =
            findfirst(
                comp_name -> comp_name == bus_var_label,
                node_syms_labels )

        push!(vars_indices, (a_var, bus_var_Idx) )
    end

    return Tuple(vars_indices)
end


"""
    get_a_node_state_algb_vars_indices_in_system(
        ; network_vars_labels = network_vars_labels,
        bus_name = bus_name,
        vars = [:ω, :δ])


Returns indices of state_algebraic variables of a node in the systems.

It is meant to extract the indices in sol.
"""
function get_a_node_state_algb_vars_indices_in_system(
    ; network_vars_labels = network_vars_labels,
    bus_name = bus_name,
    vars = [:ω, :δ])

    vars_indices = []

    for a_var in vars

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx =
            findfirst(
                comp_name -> comp_name == bus_var_label,
                network_vars_labels)

        push!(vars_indices, bus_var_Idx)
    end

    return Tuple(vars_indices)
end


"""
    get_a_node_states_vars_syms_in_system(
        ; network_vars_labels = network_vars_labels,
        bus_name = bus_name,
        vars = [:δ, :ed_dash, :eq_dash])

Returns labels of a list of state variables `vars` of a plant in the network state labels.
"""
function get_a_node_states_vars_syms_in_system(
    ; network_vars_labels = network_vars_labels,
    bus_name = bus_name,
    vars = [:δ, :ed_dash, :eq_dash])

    vars_syms = []

    for a_var in vars

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   = findfirst(
            comp_name -> comp_name == bus_var_label,
            network_vars_labels)
        

        push!(vars_syms, network_vars_labels[bus_var_Idx])
    end

    return vars_syms
end


"""
    get_a_node_states_indices_in_system(
        ; network_vars_labels = network_vars_labels,
        bus_name = bus_name,
        vars = [:δ, :ed_dash, :eq_dash])

Returns indices of a list of state variables `vars` of a plant in the network state.
"""
function get_a_node_states_indices_in_system(
    ; network_vars_labels = network_vars_labels,
    bus_name = bus_name,
    vars = [:δ, :ed_dash, :eq_dash])

    vars_indices = []

    for a_var in vars

        if bus_name == ""
            bus_var_label = Symbol("$(a_var)")
        else
            bus_var_label = Symbol("$(bus_name)_$(a_var)")
        end
        
        bus_var_Idx   = findfirst(comp_name ->
            comp_name == bus_var_label,
                                  network_vars_labels)

        push!(vars_indices, bus_var_Idx)
    end

    return vars_indices
end


"""
    get_nodes_state_algb_vars_indices_in_system(
        ; network_vars_labels =
            network_vars_labels,
        nodes_name = ["bus1", "bus2"],
        vars = [:ω, :δ])


Returns indices of state_algebraic variables of specified nodes in the systems.

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


#-------------------------------------------------------
#-------------------------------------------------------



# function get_test_system_Ynet_size(
#     system_orientations )

#     Cnb = get_Cnb_by_orientations(
#         system_orientations)

#     Cbn = get_Cbn_by_orientations(
#         system_orientations)

#     nodes_incident_edges =
#         get_nodes_incident_edges_by_orientations(
#             system_orientations )

#     Ynet_size =
#         length(nodes_incident_edges ) +
#         sum(length.(nodes_incident_edges ) )

#     return Ynet_size

# end


# function get_test_systems_Ynet_sizes(
#     list_test_systems_names,
#     list_test_systems_orientation )

#     systems_Ynet_sizes = map(
#         get_test_system_Ynet_size,
#         list_test_systems_orientation)

#     systems_no_branches = length.(
#         list_test_systems_orientation)

#     result = Any[]

#     # list_header = [["Test system", "Branches", "Ynet size"]]
    
#     # return vcat(list_header, [ [a_system, no_branches, ynet_size] for (a_system, no_branches, ynet_size) in zip(list_test_systems_names, systems_no_branches, systems_Ynet_sizes) ])

#     return [ [a_system, no_branches, ynet_size ]
#              for (a_system, no_branches, ynet_size) in zip(
#                  list_test_systems_names,
#                  systems_no_branches, systems_Ynet_sizes) ] 

# end


# function get_test_system_Cnb_Cbn_nodes_incident_edges_Ynet_size(system_orientations)

#     Cnb = get_Cnb_by_orientations(
#         system_orientations)

#     Cbn = get_Cbn_by_orientations(
#         system_orientations)

#     nodes_incident_edges =
#         get_nodes_incident_edges_by_orientations(
#             system_orientations )

#     Ynet_size =
#         length(nodes_incident_edges ) +
#         sum(length.(nodes_incident_edges ) )

#     return ( Cnb,
#              Cbn,
#              nodes_incident_edges,
#              Ynet_size )

# end


# ------------------------------------------------------


# function dynamic_nodal_current_balance( (src_ih, src_ik ), (dst_ih, dst_ik ) )

#     # dynamic_nodal_current_balance( x_s, x_d )
    
#     # src_ih, src_ik = x_s
    
#     # dst_ih, dst_ik = x_d

#     node_src_ih = edges_current_partial_sum( src_ih )
        
#     node_src_ik = edges_current_partial_sum( src_ik )
    
#     node_dst_ih = edges_current_partial_sum( dst_ih )

#     node_dst_ik = edges_current_partial_sum( dst_ik )
        

#     return (node_src_ih - node_src_ik) .- (node_dst_ih - node_dst_ik)
    

# end


# function x_dynamic_nodal_current_balance(x_s, x_d)

#     node_x_s = edges_current_partial_sum( x_s )
        
#     node_x_d = edges_current_partial_sum( x_d )

#     return  node_x_s .- node_x_d

# end




# function get_eig_values_in_states_participation(
#     eig_values,
#     PF_pure_Asys,
#     im_pure_states_syms;
#     participation_threshold = 0.2 )

#     rows_size_PF_pure_Asys,cols_size_PF_pure_Asys =
#         size(PF_pure_Asys)

#     cols_PF_pure_Asys =
#         [PF_pure_Asys[:,a_col_idx]
#          for a_col_idx in
#              1:cols_size_PF_pure_Asys]


#     rows_PF_pure_Asys =
#         [PF_pure_Asys[a_row_idx,:]
#          for a_row_idx in
#              1:rows_size_PF_pure_Asys]


#     indices_of_PF_cols_greter_than =
#         [ findall(>( participation_threshold  ),
#                   PF_pure_Asys[:, a_col] )
#           for a_col  in
#               1:cols_size_PF_pure_Asys ]


#     indices_of_PF_rows_greter_than =
#         [ findall(>( participation_threshold ),
#                   PF_pure_Asys[a_row, :] )
#           for a_row  in
#               1:rows_size_PF_pure_Asys ]


#     tup_eig_value_state_vars_PF_cols_greter_than =
#         [(eig_values[idx_eig_value],
#           im_pure_states_syms[idx_row_greater_than],
          
#          a_col_PF_pure_Asys[idx_row_greater_than] )
#          for (idx_eig_value, idx_row_greater_than,
#               a_col_PF_pure_Asys) in
#              zip(1:cols_size_PF_pure_Asys,
#                  indices_of_PF_cols_greter_than,
#                  cols_PF_pure_Asys)]


#     tup_state_vars_eig_value_PF_rows_greter_than =
#         [(im_pure_states_syms[state_idx],
#           eig_values[idx_col_greater_than],
#          a_row_PF_pure_Asys[idx_col_greater_than] )
#          for (state_idx, idx_col_greater_than,
#               a_row_PF_pure_Asys) in
#              zip(1:rows_size_PF_pure_Asys,
#                  indices_of_PF_rows_greter_than,
#                  rows_PF_pure_Asys)]

#     return (
#         ;tup_eig_value_state_vars_PF_cols_greter_than,
#         tup_state_vars_eig_value_PF_rows_greter_than )
# end


# """

# https://arblib.org/acb_mat.html#acb-mat-eigenvalues

# https://github.com/kalmarek/Arblib.jl

# matrix_rows, matrix_cols = size(system_matrix)

# eig_val = Arblib.AcbVector(zeros(matrix_rows))

# Acb_system_matrix = Arblib.AcbMatrix(system_matrix)

# eigvecs_left  = similar(Acb_system_matrix)

# eigvecs_right = similar(Acb_system_matrix)


# return_code = Arblib.approx_eig_qr!(
#     eig_val,
#     eigvecs_left,
#     eigvecs_right,
#     Acb_system_matrix,
#     Mag(),
#     0,
#     Arblib._precision(Acb_system_matrix) )

# """
# function get_eigens_via_arblib(
#     system_matrix; prec = nothing  )


#     matrix_rows, matrix_cols = size(system_matrix)

#     eig_values = Arblib.AcbVector(zeros(matrix_rows))

#     Acb_system_matrix = Arblib.AcbMatrix(system_matrix)

#     eigvecs_left  = similar(Acb_system_matrix)

#     eigvecs_right = similar(Acb_system_matrix)

#     if prec ==  nothing
#         prec = Arblib._precision(Acb_system_matrix)
#     else
#         prec =  prec
#     end
    
        
#     return_code =
#         Arblib.approx_eig_qr!(
#             eig_values,
#             eigvecs_left,
#             eigvecs_right,
#             Acb_system_matrix,
#             Mag(),
#             0,
#             prec )

#     return (; eig_values, eigvecs_left,
#             eigvecs_right, return_code )
    
# end



# """

# `get_eigens`
# This function returns eigen values, left and right
# eigen vectors.

# """
# function get_eigens(system_matrix)

#     # https://ralphas.github.io/GenericSchur.jl/stable/
    
#     system_matrix = system_matrix .+ 0im
    
#     schur_object = schur( system_matrix )

#     eig_values = schur_object.values

#     eigvecs_right = eigvecs( schur_object )

#     eigvecs_left = eigvecs(schur_object,left=true)

#     return (; eig_values, eigvecs_left, eigvecs_right )
    

# end


# """
# https://tobydriscoll.net/fnc-julia/linsys/norms.html

# https://discourse.julialang.org/t/computing-left-eigenvectors-with-high-precision/72949

# https://ralphas.github.io/GenericSchur.jl/stable/

# `get_participation_factors`
# This function returns the participation factor matrix

# Sauer: see page 232, equation 8.84


# """
# function get_participation_factors(
#     system_matrix )

#     (; eig_values, eigvecs_left, eigvecs_right ) =
#         get_eigens(system_matrix)


#     participation_factor =
#         abs.(eigvecs_left) .* abs.(eigvecs_right)

#     for i in eachindex(eig_values)

#         λ_i_pf = participation_factor[:,i] ./
#             (abs.(eigvecs_left[:,i]') *
#             abs.(eigvecs_right[:,i] ))
        
#         max_λ_i_pf = max(abs.(λ_i_pf)... )
        
#         participation_factor[:,i] .= λ_i_pf ./ max_λ_i_pf
            
#         # max_pf = max(participation_factor[:,i]... )

#         # participation_factor[:,i] .=
#         #     participation_factor[:,i] ./  max_pf 

#     end

#     return participation_factor
    
# end


#-------------------------------------------------------
