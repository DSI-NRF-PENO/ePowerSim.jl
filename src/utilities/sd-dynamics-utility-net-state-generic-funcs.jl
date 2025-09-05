# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


#-------------------------------------------------------
#-------------------------------------------------------


"""
This function is a place holder for an anonymus
function that produce nothing. It is needed in
struct that do not have some functions implemented
"""
anonymus_func = (x) -> nothing


#-------------------------------------------------------
#-------------------------------------------------------

"""
This macro converts an expression to a string

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

# test1 = Symbol.(String.(strip.(split(split(split(x,"[")[2],"]")[1], ","))))


macro macro_get_list_items_name_and_value(arg)
    x = string(arg)
    x = Symbol.(String.(strip.(split(split(split(x,"[")[2],"]")[1], ","))))

    return x, eval(arg)
    
end

#----------------------------------------

"""
functions returns second,third and
fourth element in a list

"""

function second(elements_container)
    return elements_container[2]
end



function third(elements_container)
    return elements_container[3]
end


function fourth(elements_container)
    return elements_container[4]
end


function fifth(elements_container)
    return elements_container[5]
end


function sixth(elements_container)
    return elements_container[6]
end


function seventh(elements_container)
    return elements_container[7]
end


function eighth(elements_container)
    return elements_container[8]
end


function nineth(elements_container)
    return elements_container[9]
end


function tenth(elements_container)
    return elements_container[10]
end

#----------------------------------------

function threshold_limits(x, x_max, x_min)
    return x > x_max ? x_max : x < x_min  ? x_min : x
end

function no_limit_violation(x, x_max, x_min)

    return  (x_min < x)  &&  (x < x_max)
end


function limit_violation(x, x_max, x_min)

    return  (x_min > x)  ||  (x > x_max)
end

#----------------------------------------

function polar_to_cartesian( v_θ )

    v, θ = v_θ
    
    return [v * cos(θ), v * sin(θ)]
    
end


function cartesian_to_polar(ur_ui)

    ur, ui = ur_ui
    
    return [abs(ur + im * ui),  angle( ur + im * ui )]
    
end


#----------------------------------------


function ur_ui_to_u(ur_ui)

    return ur_ui[1] + im * ur_ui[2]
    
end


function xr_xi_to_x(ur_ui)

    return ur_ui[1] + im * ur_ui[2]
    
end


function u_from_ur_ui(ur_ui)

    return ur_ui[1] + im * ur_ui[2]
    
end


function x_from_xr_xi(ur_ui)

    if length(ur_ui) == 1 && typeof(ur_ui[1]) == Float64
        return ur_ui[1] + im * 0.0
        
    elseif length(ur_ui) == 1 && typeof(ur_ui[1]) == ComplexF64
        return  0.0 + im * ur_ui[1]

    else

        return ur_ui[1] + im * ur_ui[2]
    end
    
end


function conj_x_from_xr_xi(ur_ui)

    if length(ur_ui) == 1 && typeof(ur_ui[1]) == Float64
        return ur_ui[1] - im * 0.0
        
    elseif length(ur_ui) == 1 && typeof(ur_ui[1]) == ComplexF64
        return  0.0 - im * ur_ui[1]

    else

        return ur_ui[1] - im * ur_ui[2]
    end
    
end


function u_to_ΘV(u)

    return [ angle.(u)...; abs.(u)... ]
    
end


function u_to_VΘ(u)

    return [ abs(u), angle(u) ]
    
end

function VΘ_to_u(VΘ)

    return VΘ[1] * exp(im * VΘ[2])
    
end


function ur_ui_to_ΘV(ur_ui)

    u = ur_ui_to_u.(ur_ui)

    # Θ = angle.(u) 
    # V = abs.(u)

    return [angle.(u)...;abs.(u)...]
    
end


######################################################
#-----------------------------------------------------
#  utility functions without plant, netd, etc as arguments 
#-----------------------------------------------------
######################################################


"""
https://stackoverflow.com/questions/28402989/check-size-in-bytes-of-variable-using-julia

https://discourse.julialang.org/t/is-there-a-package-to-list-memory-consumption-of-selected-data-objects/85019

"""
# varinfo()

# varinfo(r"^v$")

# Base.@locals + Base.summarysize 

# Base.summarysize(ones(10000))



"""
`get_vars_idxs_in_range_Idxs` returns vectors of
idxs as tuples

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
# Getting size of vars:
#-------------------------------------------------------

pretty_summarysize(x) = Base.format_bytes(
    Base.summarysize(x) )

getBytes(x::DataType) = sizeof(x)


function getBytes(x)
    
    total = 0;
    
    fieldNames = fieldnames(typeof(x));
    
   if fieldNames == []
      return sizeof(x);
   else
     for fieldName in fieldNames
        total += getBytes(getfield(x, fieldName));
     end
     return total;
   end
end

#----------------------------------------


"""
This function prints a tuple of flags,
representated as a named tuple

tup_flags = (; fixed_timestep, with_cb,
with_node_pf, with_global_pf )

print_flags( tup_flags )

"""
function print_flags( tup_flags )
    
    for (n, v) ∈ pairs(tup_flags)
           println("$n is $v")

    end
end

#----------------------------------------

"""
This function returns non empty lists

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
This function returns non empty lists
and their indices

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
Create offsets for stacked array of
dimensions dims
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
This function is meant to get indices of
u_r, u_i for a flattend system states:

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

Create indexes  of dimensions dims
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
This function converts non consecutive idxs
to a consecutive idxs
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
This function converts consecutive idxs
to non consecutive idxs
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


function DAE_MassMatrix(
    state_size, algebraic_size )
    
    return Diagonal(
        [ones(Int, state_size)...;
         zeros(Int, algebraic_size)...])
end


function DAE_BoolVector(
    state_size, algebraic_size)
    
    return map(
        (x) -> x ==1 ? true : false,
        [ones(Int, state_size)...;
         zeros(Int, algebraic_size)...] )

end

#------------------------------------------


function Sevf(Ae, Be, vf_tilade)
    
    return Ae * exp(Be * abs(vf_tilade))
    
end

"""
Power system modeling, computation and control

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


function get_a_flattened_by_per_vars_or_paras(
    list_vars_or_paras)


    per_vars_or_paras = [
        [per_var_or_para...;]
        for per_var_or_para in
            list_vars_or_paras  ]

    return [per_vars_or_paras...;]
    

end



"""

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
`get_per_vars_or_paras_to_per_node` converts vars or paras given
in per vars or paras format to per node format.

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
`get_per_node_para_to_per_vars_or_paras` converts vars or paras given
in per node format to per vars or paras format.

The dimension of each of the per vars or paras should be supplied in
a list `dims_vars_or_paras_types`


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

    # dims_vars_or_paras =
    #     length.( vec_of_vec_var )

    # _,_, vars_or_paras_Idx =
    #     create_size_offset_Idx(
    #         dims_vars_or_paras ;
    #         counter = 0)
    
    vars_or_paras_Idx =
        get_flattened_to_components_vector_var_Idx(
            vec_of_vec_var )
    
    #-----------------------------------------------

    flattened_vars_or_paras = [vec_of_vec_var...;]

    #-----------------------------------------------
    
    return flattened_vars_or_paras, vars_or_paras_Idx
    
end

# get_a_flat_vars_or_paras_and_Idx( [ non_gens_vh_θh ] )



function get_per_node_flat_idxs( list_vars_or_paras )

    per_node_vars_or_paras = [
        [ [ per_node_var_or_para
            for per_node_var_or_para in
                var_or_para]...;]
                  for var_or_para in
                      zip( list_vars_or_paras...) ]
    
    dims_list_vars_or_paras = length.( per_node_vars_or_paras )
    
    _, _, per_node_vars_or_paras_Idx =
        create_size_offset_Idx(
            dims_list_vars_or_paras;
            counter = 0 )

    return per_node_vars_or_paras_Idx
    
end



function get_per_node_flat_para_and_idxs( list_vars_or_paras )

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

    return (
            per_vars_or_paras_Idx,
            per_vars_or_paras_per_node_Idx
            )
    
end


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


function get_ode_flat_para_Idxs_in_Idxs(
    gens_vh_θh,
    gens_nodes_ωs_ωref0_vref0_porder0,
    gens_dynamic_id_iq_pg_vh )

    #-----------------------------------------------

    dims_gens_nodes_vh_θh =
        length.( gens_vh_θh )

    _,_, gens_nodes_vh_θh_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_vh_θh;
        counter = 0)

    #-----------------------------------------------

    dims_gens_nodes_ωs_ωref0_vref0_porder0 =
        length.( gens_nodes_ωs_ωref0_vref0_porder0 )

    _,_, gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx =
        create_size_offset_Idx(
            dims_gens_nodes_ωs_ωref0_vref0_porder0;
        counter = 0)

    #-----------------------------------------------

    dims_gens_nodes_id_iq_pg_vh =
        length.( gens_dynamic_id_iq_pg_vh )

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


#------------------------------------------
#------------------------------------------

function get_eig_values_in_states_participation(
    eig_values,
    PF_pure_Asys,
    im_pure_states_syms;
    participation_threshold = 0.2 )

    rows_size_PF_pure_Asys,cols_size_PF_pure_Asys =
        size(PF_pure_Asys)

    cols_PF_pure_Asys =
        [PF_pure_Asys[:,a_col_idx]
         for a_col_idx in
             1:cols_size_PF_pure_Asys]


    rows_PF_pure_Asys =
        [PF_pure_Asys[a_row_idx,:]
         for a_row_idx in
             1:rows_size_PF_pure_Asys]


    indices_of_PF_cols_greter_than =
        [ findall(>( participation_threshold  ),
                  PF_pure_Asys[:, a_col] )
          for a_col  in
              1:cols_size_PF_pure_Asys ]


    indices_of_PF_rows_greter_than =
        [ findall(>( participation_threshold ),
                  PF_pure_Asys[a_row, :] )
          for a_row  in
              1:rows_size_PF_pure_Asys ]


    tup_eig_value_state_vars_PF_cols_greter_than =
        [(eig_values[idx_eig_value],
          im_pure_states_syms[idx_row_greater_than],
          
         a_col_PF_pure_Asys[idx_row_greater_than] )
         for (idx_eig_value, idx_row_greater_than,
              a_col_PF_pure_Asys) in
             zip(1:cols_size_PF_pure_Asys,
                 indices_of_PF_cols_greter_than,
                 cols_PF_pure_Asys)]


    tup_state_vars_eig_value_PF_rows_greter_than =
        [(im_pure_states_syms[state_idx],
          eig_values[idx_col_greater_than],
         a_row_PF_pure_Asys[idx_col_greater_than] )
         for (state_idx, idx_col_greater_than,
              a_row_PF_pure_Asys) in
             zip(1:rows_size_PF_pure_Asys,
                 indices_of_PF_rows_greter_than,
                 rows_PF_pure_Asys)]

    return (
        ;tup_eig_value_state_vars_PF_cols_greter_than,
        tup_state_vars_eig_value_PF_rows_greter_than )
end



"""

https://arblib.org/acb_mat.html#acb-mat-eigenvalues

https://github.com/kalmarek/Arblib.jl

matrix_rows, matrix_cols = size(system_matrix)

eig_val = Arblib.AcbVector(zeros(matrix_rows))

Acb_system_matrix = Arblib.AcbMatrix(system_matrix)

eigvecs_left  = similar(Acb_system_matrix)

eigvecs_right = similar(Acb_system_matrix)


return_code = Arblib.approx_eig_qr!(
    eig_val,
    eigvecs_left,
    eigvecs_right,
    Acb_system_matrix,
    Mag(),
    0,
    Arblib._precision(Acb_system_matrix) )

"""
function get_eigens_via_arblib(
    system_matrix; prec = nothing  )


    matrix_rows, matrix_cols = size(system_matrix)

    eig_values = Arblib.AcbVector(zeros(matrix_rows))

    Acb_system_matrix = Arblib.AcbMatrix(system_matrix)

    eigvecs_left  = similar(Acb_system_matrix)

    eigvecs_right = similar(Acb_system_matrix)

    if prec ==  nothing
        prec = Arblib._precision(Acb_system_matrix)
    else
        prec =  prec
    end
    
        
    return_code =
        Arblib.approx_eig_qr!(
            eig_values,
            eigvecs_left,
            eigvecs_right,
            Acb_system_matrix,
            Mag(),
            0,
            prec )

    return (; eig_values, eigvecs_left,
            eigvecs_right, return_code )
    
end



"""

`get_eigens`
This function returns eigen values, left and right
eigen vectors.

"""
function get_eigens(system_matrix)

    # https://ralphas.github.io/GenericSchur.jl/stable/
    
    system_matrix = system_matrix .+ 0im
    
    schur_object = schur( system_matrix )

    eig_values = schur_object.values

    eigvecs_right = eigvecs( schur_object )

    eigvecs_left = eigvecs(schur_object,left=true)

    return (; eig_values, eigvecs_left, eigvecs_right )
    

end


"""
https://tobydriscoll.net/fnc-julia/linsys/norms.html

https://discourse.julialang.org/t/computing-left-eigenvectors-with-high-precision/72949

https://ralphas.github.io/GenericSchur.jl/stable/

`get_participation_factors`
This function returns the participation factor matrix

Sauer: see page 232, equation 8.84


"""
function get_participation_factors(
    system_matrix )

    (; eig_values, eigvecs_left, eigvecs_right ) =
        get_eigens(system_matrix)


    participation_factor =
        abs.(eigvecs_left) .* abs.(eigvecs_right)

    for i in eachindex(eig_values)

        λ_i_pf = participation_factor[:,i] ./
            (abs.(eigvecs_left[:,i]') *
            abs.(eigvecs_right[:,i] ))
        
        max_λ_i_pf = max(abs.(λ_i_pf)... )
        
        participation_factor[:,i] .= λ_i_pf ./ max_λ_i_pf
            
        # max_pf = max(participation_factor[:,i]... )

        # participation_factor[:,i] .=
        #     participation_factor[:,i] ./  max_pf 

    end

    return participation_factor
    
end


function get_a_im_updated_Ax(
    Ax_view, vf_tilade,
    Ax_update_para ) 

    (sparse_row_idxs,
     sparse_col_idxs,
     Ax_sparse_nzvalues ) =
         findnz( sparse( Matrix( Ax_view )) )

    (; Ae, Be,
     Ke, Te,
     vf_tilade_idx_in_Idx ) =
        Ax_update_para

    γ_idx_in_sp_nzv =
        find_V_idx_in_sparse_matrix_IJV(
            vf_tilade_idx_in_Idx,
            vf_tilade_idx_in_Idx,
            sparse_row_idxs,
            sparse_col_idxs )

    update_γ  = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

    Ax_sparse_nzvalues =
        [ Ax_sparse_nzvalues[1:γ_idx_in_sp_nzv-1];
         [update_γ];
         Ax_sparse_nzvalues[
             γ_idx_in_sp_nzv+1:end] ]

    updated_Ax =
        sparse( sparse_row_idxs,
                sparse_col_idxs,
                Ax_sparse_nzvalues)

    # Ax_view[:,:] .= updated_Ax

    return updated_Ax 
    
end



function get_updated_im_Ax(
    vec_Ax_views,
    gens_vf_tilade,
    gens_Ax_update_parameters )

    list_updated_im_Ax = []
    
    for (
        a_Ax_view,
        gens_vf_tilade,
        a_Ax_update_para ) in

        zip(vec_Ax_views,
            gens_vf_tilade,
            gens_Ax_update_parameters )

        (sparse_row_idxs,
         sparse_col_idxs,
         Ax_sparse_nzvalues ) =
             findnz( sparse( Matrix( a_Ax_view )) )

        (; Ae, Be,
         Ke, Te,
         vf_tilade_idx_in_Idx ) =
            a_Ax_update_para
        
        γ_idx_in_sp_nzv =
            find_V_idx_in_sparse_matrix_IJV(
                vf_tilade_idx_in_Idx,
                vf_tilade_idx_in_Idx,
                sparse_row_idxs,
                sparse_col_idxs )

        update_γ  = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

        Ax_sparse_nzvalues =
            [Ax_sparse_nzvalues[1:γ_idx_in_sp_nzv-1];
             [update_γ];
             Ax_sparse_nzvalues[γ_idx_in_sp_nzv+1:end] ]
        
        updated_Ax =
            sparse( sparse_row_idxs,
                    sparse_col_idxs,
                    Ax_sparse_nzvalues)

        push!(list_updated_im_Ax, updated_Ax)

        a_Ax_view[:,:] .= updated_Ax
                
    end

    return list_updated_im_Ax
    
end



function get_vf_tilade_idx_in_a_gen_im_state(plant)

    dict_im_vars_syms =
        plant.dict_im_vars_syms

    vf_tilade_idx_in_state =
        dict_im_vars_syms[ :vf_tilade ]
    
    return vf_tilade_idx_in_state
    
end


function get_a_im_Ax_update_parameters(plant)

    dict_im_vars_syms =
        plant.dict_im_vars_syms
    
    vf_tilade_idx_in_Idx  =
        dict_im_vars_syms[ :vf_tilade ]
    
    (Ta, Te, Tf,
     Tr, Ka,
     Ke, Kf,
     V_R_max, V_R_min,
     Ae, Be) =
         plant.Exc.param_values

    return (; Ae, Be, Ke, Te,  vf_tilade_idx_in_Idx )
    
end



function get_a_im_plant_system_matrices_update_parameters(
    plant)

    if only_gen_avr_in_plant(plant) || gen_gov_avr_in_plant(plant)
        dict_im_vars_syms = plant.dict_im_vars_syms

        vr1_idx       = dict_im_vars_syms[:vr1]
        vf_tilade_idx = dict_im_vars_syms[:vf_tilade]

        (Ta, Te, Tf,
         Tr, Ka,
         Ke, Kf,
         V_R_max, V_R_min,
         Ae, Be) =
             plant.Exc.param_values
        
        return (; Ae, Be, Ke, Te,  vf_tilade_idx_in_state)
        
    else
        nothing
    end
    
    
end

  

function get_a_im_plant_Ax_update_value(
    vf_tilade, plant)


    dict_im_vars_syms =
        plant.dict_im_vars_syms
    
    vr1_idx =
        dict_im_vars_syms[:vr1]
    
    vf_tilade_idx =
        dict_im_vars_syms[:vf_tilade]

    (Ta, Te, Tf,
     Tr, Ka,
     Ke, Kf,
     V_R_max, V_R_min,
     Ae, Be) =
         plant.Exc.param_values
    
    update_γ = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te  

    return (; update_γ,  vf_tilade_idx_in_state )
        
    
end

    
    
function get_a_im_plant_system_matrices_update_value(
    vf_tilade,  plant )

    if only_gen_avr_in_plant(plant) || gen_gov_avr_in_plant(plant)
        dict_im_vars_syms = plant.dict_im_vars_syms

        vr1_idx = dict_im_vars_syms[:vr1]
        vf_tilade_idx = dict_im_vars_syms[:vf_tilade]

        (Ta, Te, Tf,
         Tr, Ka,
         Ke, Kf,
         V_R_max, V_R_min,
         Ae, Be ) =
             plant.Exc.param_values

        update_γ = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te      

        return (; update_γ, vf_tilade_idx )
    else
        
        nothing

    end

            
end


function get_a_im_Ax_update_value_and_Idx(
    sparse_row_idxs,
    sparse_col_idxs,
    a_stateDiffCache,
    plant )
    
    dict_im_vars_syms = plant.dict_im_vars_syms

    vr1_idx       = dict_im_vars_syms[:vr1]
    
    vf_tilade_idx = dict_im_vars_syms[:vf_tilade]


    (Ta, Te, Tf, Tr, Ka,
     Ke, Kf, V_R_max, V_R_min,
     Ae, Be) =
         plant.Exc.param_values

    vf_tilade = a_stateDiffCache[
        dict_im_vars_syms[:vf_tilade] ]
    
    update_γ = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

    γ_idx_in_sp_nzv =
        find_V_idx_in_sparse_matrix_IJV(
            vf_tilade_idx, vf_tilade_idx,
            sparse_row_idxs, sparse_col_idxs )

    return (; update_γ, γ_idx_in_sp_nzv)
    
    
end

#----------------------------------------

"""
This function finds an ind in V of a sparse matrix, where

sparse_row_idxs[idx] = row_idx && sparse_col_idxs[idx] = col_idx 

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

# # ------------------------------------------------------
# #  makes result dir
# # ------------------------------------------------------


# function make_results_dir(
#     ;base_dir=@__DIR__,
#     res_dir="results_dir")
    
#     res_dir = joinpath(base_dir,res_dir)

#     if !(isdir(res_dir))
#         mkpath(res_dir)
#     end

#     csv_dir = joinpath(res_dir,"csv")
    
#         if !(isdir(csv_dir))
#         mkpath(csv_dir)
#         end
    
#     fig_dir = joinpath(res_dir,"fig")
    
#         if !(isdir(fig_dir))
#         mkpath(fig_dir)
#         end

#     sol_dir = joinpath(res_dir,"sol")
    
#         if !(isdir(sol_dir))
#         mkpath(sol_dir)
#         end
    
#     return (; res_dir, csv_dir,
#             fig_dir, sol_dir)
# end


#--------------------------------------
#--------------------------------------


function get_residue_prototype_size(
    nodes_collection )

    """ any common field can be used Bus_type

        We have 2 x no_nodes eqs - 2 + 2 x stator eqs

        The -2 is when the eqs of P and Q balance for the
        slack node is not includes
        
        the residue size: 2 x no_nodes  + 2 x no_gens - 2
    
    """
    
    if isa(nodes_collection,  OrderedDict)
        
        comp_collection_values =
            collect(values(nodes_collection))
        
        no_of_nodes =
            length([ a_node.Bus_type
                     for a_node in
                         comp_collection_values ])
        
        no_of_gens =
            length([ a_node.Bus_type
                     for a_node in
                         comp_collection_values
                         if a_node.Bus_type == :Generator ])

        return 2 * (no_of_nodes + no_of_gens ) -2
        
    elseif isa(nodes_collection, Union{Array, Vector})
        
        no_of_nodes =
            length([ a_node.Bus_type
                     for a_node in
                         nodes_collection ])
        
        no_of_gens =
            length([ a_node.Bus_type
                     for a_node in
                         nodes_collection
                         if a_node.Bus_type == :Generator ])

        return 2 * (no_of_nodes + no_of_gens ) -2

    else
        
        comp_collection_values =
            collect(values(nodes_collection))
        
        no_of_nodes =
            length([ a_node.Bus_type
                     for a_node in
                         comp_collection_values ])
        
        no_of_gens =
            length([ a_node.Bus_type
                     for a_node in
                         comp_collection_values
                         if a_node.Bus_type == :Generator ])

        return 2 * (no_of_nodes + no_of_gens ) -2

    end

end

#--------------------------------------
#--------------------------------------


"""
`net_to_industrial_model_indices`

This function converts indices in a network
model to indices in industrial model.

example:

ur_idx:

net_to_industrial_model_indices( ur_idx ;
    indices_in_system=pure_indices_in_system,
    indices_in_industrial_model=pure_indices_in_industrial )

nodes_u_Idx:

net_to_industrial_model_indices( nodes_u_Idx ;
    indices_in_system=pure_indices_in_system,
    indices_in_industrial_model=pure_indices_in_industrial )



"""
function net_to_industrial_model_indices(
    system_idx ;
    indices_in_system=pure_indices_in_system,
    indices_in_industrial_model=pure_indices_in_industrial )

    dict_sys_to_industry =
        OrderedDict{Int64, Int64}(
            sys => ind for (sys, ind) in
                zip( indices_in_system,
                     indices_in_industrial_model  ) )

    if typeof(system_idx) == Int64 
        
       return dict_sys_to_industry[system_idx ]

    elseif typeof(system_idx) == Vector{Int64} || typeof(system_idx) == Vector{Any}

        return [dict_sys_to_industry[ idx ]
                for idx in system_idx ]
         
    elseif typeof(system_idx) == Vector{Vector{Int64}} || typeof(system_idx) == Vector{Vector{Any}}

        return [[dict_sys_to_industry[ idx ]
                 for idx in vec_idx ]
                for vec_idx in system_idx ]
        
    else
        return nothing
    end
    
end


function net_to_industrial_model_indices(
    system_idx, dict_sys_to_industry  )

    if typeof(system_idx) == Int64
        
       return dict_sys_to_industry[system_idx ]

    elseif typeof(system_idx) == Vector{Int64} || typeof(system_idx) == Vector{Any}

        return [dict_sys_to_industry[ idx ]
                for idx in system_idx ]
         
    elseif typeof(system_idx) == Vector{Vector{Int64}} || typeof(system_idx) == Vector{Vector{Any}}

        return [[dict_sys_to_industry[ idx ]
                 for idx in vec_idx ]
                for vec_idx in system_idx ]
        
    else
        return nothing
    end
    
end



# ------------------------------------------------------
# ------------------------------------------------------


function net_to_im_model_indices(
    system_idx ;
    indices_in_system=
        im_vars_indices_in_system,
    im_vars_indices=
        im_vars_indices )

    dict_sys_to_industry =
        OrderedDict{Int64, Int64}(
            sys => ind for (sys, ind) in
                zip( indices_in_system,
                     im_vars_indices  ) )

    if typeof(system_idx) == Int64 
        
       return dict_sys_to_industry[system_idx ]

    elseif typeof(system_idx) == Vector{Int64} || typeof(system_idx) == Vector{Any}

        return [dict_sys_to_industry[ idx ]
                for idx in system_idx ]
         
    elseif typeof(system_idx) == Vector{Vector{Int64}} || typeof(system_idx) == Vector{Vector{Any}}

        return [[dict_sys_to_industry[ idx ]
                 for idx in vec_idx ]
                for vec_idx in system_idx ]
        
    else
        return nothing
    end
    
end


function net_to_im_model_indices(
    system_idx, dict_sys_to_industry  )

    if typeof(system_idx) == Int64
        
       return dict_sys_to_industry[system_idx ]

    elseif typeof(system_idx) == Vector{Int64} || typeof(system_idx) == Vector{Any}

        return [dict_sys_to_industry[ idx ]
                for idx in system_idx ]
         
    elseif typeof(system_idx) == Vector{Vector{Int64}} || typeof(system_idx) == Vector{Vector{Any}}

        return [[dict_sys_to_industry[ idx ]
                 for idx in vec_idx ]
                for vec_idx in system_idx ]
        
    else
        return nothing
    end
    
end


# ------------------------------------------------------
# dynamic nodal current balance  
# ------------------------------------------------------


function get_nodes_incident_edges_to_or_fro(
    to_or_fro_edges_current )
    
    return sum.([ ih_or_ik == [] ?
        x_from_xr_xi([0.0, 0.0]) :
        x_from_xr_xi.([ [current_set...;]
                        for current_set in ih_or_ik ])
                  for ih_or_ik in
                      to_or_fro_edges_current ])
end


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

function Z_dq(ra, X_d_dash, X_q_dash)

    return [ ra         -X_q_dash;
             X_d_dash    ra]

end



# function invZ_dq(ra, X_d_dash, X_q_dash)

#     return (1.0 / ( ra * ra + X_d_dash * X_q_dash )) *
#         [ra X_q_dash; -X_d_dash ra]

# end



function invZ_dq(ra, X_d_dash, X_q_dash)

    return inv( Z_dq(ra, X_d_dash, X_q_dash) )

end


"""
network reference frame to machine reference frame
Sauer page 160


"""
function network_to_machine_ref_frame(x_r, x_i, δ)

    A = [  sin(δ)  -cos(δ);
           cos(δ)   sin(δ)]
    
    x_ri = [x_r, x_i]

    x_dq = A \ x_ri

    return (x_d = x_dq[1], x_q = x_dq[2])

end


"""
 machine reference frame to network reference frame
Sauer page 160


"""
function machine_to_network_ref_frame(x_d, x_q, δ)

    A = [  sin(δ)    cos(δ);
          -cos(δ)   sin(δ)]
    
    x_dq = [x_d, x_q]

    x_ri = A \ x_dq

    return (x_r = x_ri[1], x_i = x_ri[2])

end


#------------------------------------------------
#------------------------------------------------


"""

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

nodes_size = 9

padded_vec_of_vec = [[0.0, 0.0] for idx in 1:nodes_size ]

P_comp = [1.2, 3.0, 4.0]
Q_comp = [2.2, 3.1, 1.0]

idxs = [2,5,9]

update_padded_vec_of_vec_from_comp_axis!(
    padded_vec_of_vec,
    P_comp, Q_comp,
    idxs )

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
a = [[0.0, 0.0], [0.0, 0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0] ]

t_p = [[1.2, 3.0], [4.0, 6.0]]

idxs = [2, 4, 6]

update_vec_of_vec_from_vec_vec!(
    a, t_p, idxs )


"""
function update_vec_of_vec_from_vec_vec!(
    vec_of_vec, PQ_vec_vec, idxs )

    for (idx, PQ_vec) in zip(idxs, PQ_vec_vec )
        vec_of_vec[idx] .= PQ_vec
    end
    
end

#------------------------------------------------
#------------------------------------------------

function update_vec_of_vec_from_flattened!(
    vec_of_vec, flattened_vec, vec_vec_Idx )

    for idx in 1:length(vec_vec_Idx )
        vec_of_vec[idx] .=
            flattened_vec[ vec_vec_Idx[idx] ]
    end
    

end


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
#------------------------------------------------


function get_gen_nodes_ω_ed_dash_eq_dash_views(
    state, nodes_ω_ed_dash_eq_dash_Idxs ) 

    return [ view(state, idx )
             for idx in
                 nodes_ω_ed_dash_eq_dash_Idxs  ]
    
end



function get_gen_nodes_ω_ed_dash_eq_dash(
    state, nodes_ω_ed_dash_eq_dash_Idxs )
    
    return [state[idx]
            for idx in
                nodes_ω_ed_dash_eq_dash_Idxs  ]
    
end


function get_gen_nodes_δ_ω_ed_dash_eq_dash_views(
    state, nodes_δ_ω_ed_dash_eq_dash_Idxs )
    
    return [ view(state, idx)
             for idx in
                 nodes_δ_ω_ed_dash_eq_dash_Idxs ]
    
end


function get_gen_nodes_δ_ω_ed_dash_eq_dash(
    state, nodes_δ_ω_ed_dash_eq_dash_Idxs )
    
    return [state[idx]
            for idx in
                nodes_δ_ω_ed_dash_eq_dash_Idxs ]
    
end

#------------------------------------------------
#------------------------------------------------



function get_industrial_gen_nodes_ω_ed_dash_eq_dash_views(
    state, nodes_ω_ed_dash_eq_dash_Idxs ) 

    return [
        view(state, idx )
        for idx in  nodes_ω_ed_dash_eq_dash_Idxs  ]
    
end



function get_industrial_gen_nodes_ω_ed_dash_eq_dash(
    state, nodes_ω_ed_dash_eq_dash_Idxs ) 

    return [
        state[ idx]
        for idx in  nodes_ω_ed_dash_eq_dash_Idxs  ]
    
end



# function get_industrial_gen_nodes_ω_ed_dash_eq_dash(
#     state, nodes_ω_ed_dash_eq_dash_Idxs )
    
#     return [
#         state[idx]
#         for idx in nodes_ω_ed_dash_eq_dash_Idxs  ]
    
# end



function get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
    state, nodes_δ_ω_ed_dash_eq_dash_Idxs )
    
    return [
        view(state, idx)
        for idx in  nodes_δ_ω_ed_dash_eq_dash_Idxs  ]
    
end


function get_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash(
    state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )
    
    return [
        state[idx]
        for idx in nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state  ]
    
end


function update_industrial_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
    gen_nodes_δ_ω_ed_dash_eq_dash_views,
    state, nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    for (a_gen_δ_ω_ed_dash_eq_dash_view, a_gen_δ_ω_ed_dash_eq_dash_Idxs_in_state) in zip(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state)
        
        a_gen_δ_ω_ed_dash_eq_dash_view .=
            state[ a_gen_δ_ω_ed_dash_eq_dash_Idxs_in_state ]
    end

    return nothing
    
end


#------------------------------------------------
#------------------------------------------------


function get_gen_im_nodes_ω_ed_dash_eq_dash_views(
    state,
    nodes_ω_ed_dash_eq_dash_Idxs ) 

    return [
        view(state, idx )
        for idx in  nodes_ω_ed_dash_eq_dash_Idxs  ]
    
end



function get_gen_im_nodes_ω_ed_dash_eq_dash(
    state,
    nodes_ω_ed_dash_eq_dash_Idxs )
    
    return [
        state[idx]
        for idx in nodes_ω_ed_dash_eq_dash_Idxs  ]
    
end



function get_im_each_gen_nodes_δ_ω_ed_dash_eq_dash_views(
    state,
    nodes_δ_ω_ed_dash_eq_dash_Idxs )
    
    return [
        view(state, idx)
        for idx in  nodes_δ_ω_ed_dash_eq_dash_Idxs  ]
    
end


function get_im_each_gen_nodes_δ_ω_ed_dash_eq_dash(
    state,
    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )
    
    return [
        state[idx]
        for idx in  nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state ]
    
end


function update_im_each_gen_nodes_δ_ω_ed_dash_eq_dash!(
    gen_nodes_δ_ω_ed_dash_eq_dash_views,
    state,
    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state )

    for (a_gen_δ_ω_ed_dash_eq_dash_view, a_gen_δ_ω_ed_dash_eq_dash_Idxs_in_state) in zip(
        gen_nodes_δ_ω_ed_dash_eq_dash_views,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state)
        
        a_gen_δ_ω_ed_dash_eq_dash_view .= state[ a_gen_δ_ω_ed_dash_eq_dash_Idxs_in_state ]
    end

    return nothing
    
end

#------------------------------------------------
#------------------------------------------------


function get_gen_nodes_states_from_state(
    state, gens_states_indices_in_system )
   
    return state[[gens_states_indices_in_system...;] ]
    
end


function get_gen_nodes_pure_states_from_state(
    state, gens_pure_states_indices_in_system )
   
    return state[[gens_pure_states_indices_in_system...;] ]
    
end


function get_gen_nodes_stab_states_from_state(
    state, gens_stab_states_indices_in_system )
   
    return state[[gens_stab_states_indices_in_system...;] ]
    
end


function get_each_gen_nodes_states_from_state(
    state, gens_states_indices_in_state )
   
    return [
        state[idx] for idx in
            gens_states_indices_in_state ]
    
end


function get_each_gen_nodes_pure_states_from_state(
    state, gens_pure_states_indices_in_state )
   
    return [
        state[idx] for idx in
            gens_pure_states_indices_in_state ]
    
end


function get_each_gen_nodes_stab_states_from_state(
    state, gens_states_indices_in_state )
   
    return [
        state[idx] for idx in
            gens_states_indices_in_state ]
    
end

#------------------------------------------------


function get_gen_nodes_pure_state_views(
    state, gens_pure_states_indices_in_state )

    return [
        view(state, idx) for idx in
            gens_pure_states_indices_in_state  ]
    
end


function get_gen_nodes_stab_state_views(
    state, gens_stab_states_indices_in_state )
    
    return [
        view(state, idx) for idx in
            gens_stab_states_indices_in_state ]
    
end


function get_gen_nodes_state_views(
    state, gens_states_indices_in_state, gens_idx )

    return [
        view(state, idx)  for idx in
            gens_states_indices_in_state  ]
    
end


function get_gen_nodes_pure_state_views(
    state,
    gens_pure_states_indices_in_state,
    gens_idx )

    return [
        view(state, idx)
        for idx in
            gens_pure_states_indices_in_state  ]
    
end


function get_gen_nodes_stab_state_views(
    state,
    gens_stab_states_indices_in_state,
    gens_idx )

    return [
        view(state, idx)  for idx in
            gens_stab_states_indices_in_state  ]
    
end

#------------------------------------------------



function get_gen_nodes_pure_state(
    state,
    gens_pure_states_indices_in_state )

    return [
        state[idx] for idx in
            gens_pure_states_indices_in_state  ]
    
end


function get_gen_nodes_stab_state(
    state,
    gens_stab_states_indices_in_state )
    
    return [
        state[idx] for idx in
            gens_stab_states_indices_in_state ]
    
end


function get_gen_nodes_state(
    state,
    gens_states_indices_in_state,
    gens_idx )

    return [
        state[idx]  for idx in
            gens_states_indices_in_state  ]
    
end


function get_gen_nodes_pure_state(
    state,
    gens_pure_states_indices_in_state,
    gens_idx )

    return [
        state[idx]
        for idx in
            gens_pure_states_indices_in_state  ]
    
end


function get_gen_nodes_stab_state(
    state, gens_stab_states_indices_in_state, gens_idx )

    return [
        state[ idx]
        for idx in
            gens_stab_states_indices_in_state  ]
    
end



             
#------------------------------------------------
#------------------------------------------------


function get_industrial_ur_ui_view_in_state(
    state, industrial_model_ur_ui_Idx )
   
    return view(state, industrial_model_ur_ui_Idx)
    
end


function get_industrial_pure_states_view_in_state(
    state, industrial_model_pure_states_Idx )
    
    return view(state, industrial_model_pure_states_Idx)
    
end


function get_industrial_stab_states_view_in_state(
    state, industrial_model_stab_states_Idx )
   
    return view(state, industrial_model_stab_states_Idx) 
    
end


#------------------------------------------------



function get_industrial_ur_ui_in_state(
    state, industrial_model_ur_ui_Idx )
   
    return [state[idx] for idx in
                industrial_model_ur_ui_Idx]
    
end


function get_industrial_pure_states_in_state(
    state, industrial_model_pure_states_Idx )
    
    return [state[idx] for idx in
                industrial_model_pure_states_Idx]
    
end


function get_industrial_stab_states_in_state(
    state, industrial_model_stab_states_Idx )
   
            return [state[idx] for idx in
                        industrial_model_stab_states_Idx]
    
end




#------------------------------------------------


function get_im_ur_ui_view_in_state(
    state, im_model_ur_ui_Idx )
   
    return view(state, im_model_ur_ui_Idx)
    
end


function get_im_pure_states_view_in_state(
    state, im_pure_states_Idx )
    
    return view(state, im_pure_states_Idx )
    
end



function get_im_vars_view_in_state(
    state, im_vars_Idx_in_state )
    
    return view(state, im_vars_Idx_in_state )
    
end

#------------------------------------------------


function get_im_ur_ui_in_state(
    state, im_model_ur_ui_Idx )
   
    return [state[idx]
            for idx in
                im_model_ur_ui_Idx]
    
end


function get_im_pure_states_in_state(
    state, im_pure_states_Idx )
    
    return [state[idx]
            for idx in
                im_pure_states_Idx ]
    
end


function get_im_vars_in_state(
    state, im_vars_Idx_in_state )
    
    return [state[idx]
            for idx in
                im_vars_Idx_in_state ]
    
end

#------------------------------------------------


function get_industrial_gen_nodes_states_from_state(
    state, states_indices )
    
    return  state[ [states_indices...;] ] 
    
end


function get_industrial_gen_nodes_pure_states_from_state(
    state, gens_pure_states_indices_in_state )
   
    return state[
        [ gens_pure_states_indices_in_state...; ]  ]
    
end


function get_gen_nodes_im_vars_from_state(
    state, im_vars_Idx_in_state )
   
    return state[
        [ im_vars_Idx_in_state...; ]  ]
    
end



function get_industrial_gen_nodes_stab_states_from_state(
    state, gens_stab_states_indices_in_state )
   
    return state[
        [ gens_stab_states_indices_in_state...; ] ] 
    
end

#------------------------------------------------


function get_industrial_each_gen_nodes_states_from_state(
    state, gens_states_indices_in_state )
   
    return [state[idx]
            for idx in
                gens_states_indices_in_state ]
    
end


function get_industrial_each_gen_nodes_pure_states_from_state(
    state, gens_pure_states_indices_in_state )
   
    return [state[idx]
            for idx in
                gens_pure_states_indices_in_state ]
    
end


function get_industrial_each_gen_nodes_stab_states_from_state(
    state, gens_stab_states_indices_in_state )
   
    return [ state[idx]
             for idx in
                 gens_stab_states_indices_in_state ]
    
end

#------------------------------------------------


function get_industrial_gen_nodes_pure_state_views(
    state,
    gens_pure_states_indices_in_state )

    return [
        view(state, idx)
        for idx in
            gens_pure_states_indices_in_state ]
    
end


function get_industrial_gen_nodes_stab_state_views(
    state,
    gens_stab_states_indices_in_state )
    
    return [
        view(state, idx)
        for idx in
            gens_stab_states_indices_in_state ]
    
end


function get_industrial_gen_nodes_state_views(
    state,
    gens_states_indices_in_system,
    gens_idx )

    return [
        view(state, idx)
        for idx in
            gens_states_indices_in_system  ]
    
end


function get_industrial_gen_nodes_pure_state_views(
    state,
    gens_pure_states_indices_in_system,
    gens_idx )

    return [
        view(state, idx)
        for idx in
            gens_pure_states_indices_in_system  ]
    
end


function get_industrial_gen_nodes_stab_state_views(
    state,
    gens_stab_states_indices_in_system,
    gens_idx )

    return [
        view(state, idx)
        for idx in
            gens_stab_states_indices_in_system  ]
    
end

#------------------------------------------------



function get_industrial_gen_nodes_pure_state(
    state, gens_pure_states_indices_in_state )

    return [ state[ idx ]
             for idx in
                 gens_pure_states_indices_in_state ]
    
end


function get_industrial_gen_nodes_stab_state(
    state,
    gens_stab_states_indices_in_state )
    
    return [ state[ idx ]
             for idx in
                 gens_stab_states_indices_in_state ]
    
end


function get_industrial_gen_nodes_state(
    state,
    gens_states_indices_in_system,
    gens_idx )

    return [ state[ idx ]
             for idx in
                 gens_states_indices_in_system  ]
    
end


function get_industrial_gen_nodes_pure_state(
    state,
    gens_pure_states_indices_in_system,
    gens_idx )

    return [ state[ idx ]
             for idx in
                 gens_pure_states_indices_in_system  ]
    
end


function get_industrial_gen_nodes_stab_state(
    state,
    gens_stab_states_indices_in_system,
    gens_idx )

    return [ state[ idx ]
             for idx in
                 gens_stab_states_indices_in_system  ]
    
end


#------------------------------------------------
#------------------------------------------------


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
    gens_nodes_idx =
        gens_nodes_idx,
    nodes_size = nodes_size )

    
    x_vh = [ idx ∈ gens_nodes_idx ?
        gens_vh[ n2s_gens_idx[ idx ] ] : abs( uh[ idx ] )
              for idx in 1:nodes_size ]
    
    x_θh = [ angle( uh[ idx ] ) for idx in 1:nodes_size ]

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
#------------------------------------------------


function get_dynamic_idq_vhθh(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )


    return  get_dynamic_idq_vhθh(
        vh_θh..., δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end


function get_dynamic_idq_vhθh(
    vh,
    θh,
    δ_ω_ed_eq,
    ra_Xd_dash_Xq_dash )
        
    return  get_dynamic_idq_vhθh(
        vh, θh, δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash... )  
end




function get_dynamic_idq_vhθh(
    vh,
    θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )

    if δ_ω_ed_eq == []
        
        return [0.0, 0.0]
        
    else
        
        δ, ω, ed_dash, eq_dash = δ_ω_ed_eq
        
        ra, X_d_dash, X_q_dash = ra_X_d_dash_X_q_dash

        #id_iq = invZ_dq(ra, X_d_dash, X_q_dash) * [ed_dash - vh * sin(δ - θh), eq_dash - vh * cos(δ - θh)]
        
        return invZ_dq(ra, X_d_dash, X_q_dash) * [
            ed_dash - vh * sin(δ - θh), eq_dash -
                vh * cos(δ - θh)]
    end
    
end


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

function get_dynamic_idq_θ_π_vhθh(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )


    return get_dynamic_idq_θ_π_vhθh(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end


function get_dynamic_idq_θ_π_vhθh(
    vh,
    θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )

    if δ_ω_ed_eq == []
        
        return 0.0 + im * 0.0
        
    else
        
        # δ, ω, ed_dash, eq_dash = δ_ω_ed_eq        
        # ra, X_d_dash, X_q_dash = ra_X_d_dash_X_q_dash
        # id_iq = invZ_dq(ra, X_d_dash, X_q_dash) * [
        #     ed_dash - vh * sin(δ - θh), eq_dash -
        #         vh * cos(δ - θh)]
        # return  (id_iq[1] + im * id_iq[2]) * exp(
        #     im * (δ - pi/2))

        return get_dynamic_idq_θ_π_vhθh(
            vh_θh...,
            δ_ω_ed_eq...,
            ra_X_d_dash_X_q_dash... )
    end
    
end


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

# -------------------------------------------


function industrial_model_get_dynamic_idq_θ_π_vhθh(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )


    return industrial_model_get_dynamic_idq_θ_π_vhθh(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end


function industrial_model_get_dynamic_idq_θ_π_vhθh(
    vh, θh,
    δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )
        
    id_iq = invZ_dq(ra, X_d_dash, X_q_dash) * [
        ed_dash - vh * sin(δ - θh), eq_dash -
            vh * cos(δ - θh)]

    return  (id_iq[1] + im * id_iq[2]) * exp(
        im * (δ - pi/2))    
end

# -------------------------------------------

function get_pf_dynamic_idq_θ_π_vhθh(
    vh_θh, δ_ω_ed_eq, ra_X_d_dash_X_q_dash )


    return  get_pf_dynamic_idq_θ_π_vhθh(
        vh_θh..., δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end


function get_pf_dynamic_idq_θ_π_vhθh(
    vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash )
        
    return  get_pf_dynamic_idq_θ_π_vhθh(
        vh, θh, δ_ω_ed_eq..., ra_Xd_dash_Xq_dash... )  
end

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


# -------------------------------------------


function get_pf_dyn_idq_θ_π_vhθh(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )

    return  get_pf_dyn_idq_θ_π_vhθh(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end


function get_pf_dyn_idq_θ_π_vhθh(
    vh, θh,
    δ_ω_ed_eq,
    ra_Xd_dash_Xq_dash )
        
    return  get_pf_dyn_idq_θ_π_vhθh(
        vh, θh,
        δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash... )  
end


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



function get_pf_dyn_idq_net( idq, δ )
    
    return  get_pf_dyn_idq_net( idq..., δ )

end


function get_pf_dyn_idq_net( id, iq, δ )

    idq = (id + im * iq)  * exp(im * (δ - π/2))
    return  [real(idq ), imag(idq )]

end

# -------------------------------------------

function get_pf_dyn_idq(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )

    return  get_pf_dyn_idq(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end



function get_pf_dyn_idq(
    vh, θh,
    δ_ω_ed_eq,
    ra_Xd_dash_Xq_dash )
        
    return  get_pf_dyn_idq(
        vh, θh,
        δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash... )  
end


function get_pf_dyn_idq(
    vh, θh, δ, ω,
    ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    return  invZ_dq(ra, X_d_dash, X_q_dash) *
        [ed_dash - vh * sin(δ - θh),
         eq_dash - vh * cos(δ - θh)]

end

# -------------------------------------------


function get_a_gen_dyn_idq(
    vh_θh,
    δ_ω_ed_eq,
    ra_X_d_dash_X_q_dash )

    return  get_a_gen_dyn_idq(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_X_d_dash_X_q_dash... )
    
end



function get_a_gen_dyn_idq(
    vh, θh,
    δ_ω_ed_eq,
    ra_Xd_dash_Xq_dash )
        
    return  get_a_gen_dyn_idq(
        vh, θh,
        δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash... )  
end


function get_a_gen_dyn_idq(
    vh, θh,
    δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    return  invZ_dq(ra, X_d_dash, X_q_dash) *
        [ed_dash - vh * sin(δ - θh),
         eq_dash - vh * cos(δ - θh)]

end



function get_a_gen_dyn_idq(
    vh, θh,
    δ, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    return  invZ_dq(ra, X_d_dash, X_q_dash) *
        [ed_dash - vh * sin(δ - θh),
         eq_dash - vh * cos(δ - θh)]

end

# -------------------------------------------


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


function get_dynamic_idq_ur_ui(
    u_r, u_i, δ, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash)

    return invZ_dq(ra, X_d_dash, X_q_dash) * [
        ed_dash, eq_dash] - invZ_dq(
            ra, X_d_dash, X_q_dash) * [
                sin(δ) -cos(δ); cos(δ) sin(δ)] * [
                    u_r, u_i]
    
end



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

# #------------------------------------------------

# function Ri_vh_θh( vh, θh, δ, ω, ed_dash, eq_dash, H, ωs, ra, Xd, Xq, Xd_dash, Xq_dash, T_d_dash, T_q_dash )

#     id_iq = get_dynamic_idq_vhθh( vh, θh, δ, ω, ed_dash, eq_dash, ra, X_d_dash, X_q_dash )

#     id, iq = id_iq

#     τe = get_dynamic_pg_from_id_iq( id, iq, δ, ω, ed_dash, eq_dash,  ra, X_d_dash, X_q_dash)

#     sm = [ 0, -τe * ωs /( 2 * H ), ( Xq - Xq_dash ) * iq / T_q_dash, ( Xd - Xd_dash ) * id / T_d_dash ]

#     gov = [-ω/( Ts * R * ωs ), 0  ]

#     avr = [ -Ka * vh/Ta, 0, -vf_tilade * Sevf(Ae, Be, vf_tilade) / Te ]

# end


# function Ri_ur_ui( ur, ui, δ, ω, ed_dash, eq_dash, ra, Xd, Xq, Xd_dash, Xq_dash )

# end


# #------------------------------------------------


function get_dynamic_τe_from_id_iq(
    id, iq, δ, ed_dash, eq_dash,  ra,
    X_d_dash, X_q_dash)

    return ed_dash * id + eq_dash * iq +
        ( X_q_dash - X_d_dash ) * id * iq
    
end


function get_dynamic_τe_from_id_iq(
    id, iq, δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash)

    return ed_dash * id + eq_dash * iq +
        ( X_q_dash - X_d_dash ) * id * iq
    
end



function get_dynamic_pg_from_id_iq(
    id, iq, δ, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash)

    return ed_dash * id + eq_dash * iq +
        ( X_q_dash - X_d_dash ) * id * iq
    
end


function get_dynamic_pg_from_id_iq(
    id, iq, δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash)

    return ed_dash * id + eq_dash * iq +
        ( X_q_dash - X_d_dash ) * id * iq
    
end


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

function get_dynamic_ph_by_vhθh_δ_idq(
    vh, θh, δ, id, iq )

    return id * vh * sin( δ - θh ) +
        iq * vh * cos( δ - θh )
    
end


function get_dynamic_ph_by_vhθh_δ_idq(
    vh, θh, δ, ω, ed_dash,
    eq_dash, id, iq )

    return id * vh * sin( δ - θh )  +
        iq * vh * cos( δ - θh )
    
end


function get_dynamic_qh_by_vhθh_δ_idq(
    vh, θh, δ, id, iq )

    return id * vh * cos( δ - θh ) -
        iq * vh * sin( δ - θh )
    
end


function get_dynamic_qh_by_vhθh_δ_idq(
    vh,
    θh,
    δ,
    ω,
    ed_dash,
    eq_dash,
    id,
    iq )

    return id * vh * cos( δ - θh ) -
        iq * vh * sin( δ - θh )
    
end



function get_gens_dynamic_ph_by_vhθh_δ_idq(
    gens_vh,
    gens_θh,
    gens_δ,
    gens_id,
    gens_iq )

    return [
        get_dynamic_ph_by_vhθh_δ_idq(vh, θh, δ, id, iq)
             for (vh, θh, δ, id, iq) in
                 zip(gens_vh,
                     gens_θh,
                     gens_δ,
                     gens_id,
                     gens_iq ) ]
    
end



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


#------------------------------------------------

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


function get_component_δ_ω_ed_dash_eq_dash_from_pf(
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

    return [δ, ωs, ed_dash, eq_dash] 
    

end



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
    X_q_dash; ω = ωs  )

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

#------------------------------------------
#------------------------------------------

function update_gens_τm_vf!(
    vec_τm_vf_views, gens_τm_vf  )

    for ( a_gen_τm_vf_view, a_gen_τm_vf) in zip(
        vec_τm_vf_views, gens_τm_vf  )

        a_gen_τm_vf_view .= a_gen_τm_vf

    end
    

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
#------------------------------------------------

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
#------------------------------------------------


function make_case_buses_names(
    ; case_fun = case_fun )

    network_bus_names  =
        get_network_bus_names(
            case_fun()[1] )
    
    non_gens_bus_names =
        get_non_generators_bus_names(
            case_fun()[1] )

    gens_bus_names =
        get_generators_bus_names(
            case_fun()[1] )
    
    Loads_bus_names =
        get_Loads_bus_names(
            case_fun()[1] )
    
    Trans_bus_names =
        get_Trans_bus_names(
            case_fun()[1] )

    return (; network_bus_names,
            non_gens_bus_names,
            gens_bus_names,
            Loads_bus_names,
            Trans_bus_names )
    
end


#------------------------------------------------
#------------------------------------------------


function generate_industrial_model_sym_and_mass_matrix(
    gens_nodes_pure_states_labels,
    net_bus_volts_labels )


    industrial_model_sym = vcat(
        gens_nodes_pure_states_labels,
        net_bus_volts_labels )
    
    
    industrial_model_mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(
        length( gens_nodes_pure_states_labels ),
        length( net_bus_volts_labels  ))

    return (; industrial_model_sym,
            industrial_model_mass_matrix )

end

#------------------------------------------------
#------------------------------------------------


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
#------------------------------------------------


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
             :nodes_edges_dynamics! => a_nodes_edges_dyn_fun])
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


    ntup_sim_model_types = create_namedtup_sim_model_types()
    
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

    # (; sim_model_types_and_flags =
    #     sim_model_types_and_flags,
    #  list_sym_model_type_flags =
    #      list_sym_model_type_flags)

    list_model_type  = [
        :initial_pf,
        :network_current_pf,
        :global_pf,
        :node_pf,
        :hybrid_pf ]

    
    sim_model_types_and_flags, list_sym_model_type_flags =
        create_simulation_model_types_and_flags( )
    
    # sim_namedtup_for_model_types
    
    sim_namedtup_for_model_types = namedtuple([
        model_type => ( sim_model_type_and_flags ,
                        sym_model_type_flags )
        for ( model_type,
              sim_model_type_and_flags,
              sym_model_type_flags ) in
            zip( list_model_type,
                 sim_model_types_and_flags,
                 list_sym_model_type_flags )])

    return  sim_namedtup_for_model_types


end


function create_a_simulation_model_type_and_flags(
    model_type, sim_namedtup_for_model_types )
    
    return getproperty(sim_namedtup_for_model_types,
                       model_type )
    
end

#------------------------------------------------
#------------------------------------------------
