# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# ####################################################

# # ------------------------------------------------------
# # ------------------------------------------------------

# function z2y(; r =1.0,  x = 1.00, G = 0.0,  B = 0.0)

#     z = r + im * x

#     y = 1/z

#     B_2 = B/2.0
    

#     return (y = y,  y_shunt_km = B_2, y_shunt_mk = B_2)

# end

# function Qmax_Qmin_limit_violation(genQ, gen_Qmax, gen_Qmin)
#     Qmax_limit_violation = genQ .> gen_Qmax
#     Qmim_limit_violation = genQ .< gen_Qmin
#     return Qmax_limit_violation, Qmim_limit_violation
# end


# function no_limit_violation(x, x_max, x_min)

#     return  (x_min < x)  &&  (x < x_max)
# end


# function limit_violation(x, x_max, x_min)

#     return  (x_min > x)  ||  (x > x_max)
# end


# # function get_gens_bus_num(gen_nodes)

# #     return map(get_node_Bus_num,
# #                collect(values(gen_nodes)))
# # end


# function get_slack_and_gens_Cg(nodes)

#     num_nodes = length( collect(keys(nodes)) )
        
#     slack_or_gen_nodes_Idx =
#         get_slack_and_gen_nodes_Idx(nodes)

#     num_slack_or_gen_nodes = length( slack_or_gen_nodes_Idx )
    
#     Cg = sparse(slack_or_gen_nodes_Idx, collect(1:num_slack_or_gen_nodes), ones(Int,length( slack_or_gen_nodes_Idx )), num_nodes, num_slack_or_gen_nodes )

#     return Cg
# end



# function get_gens_Cg(gen_nodes, num_nodes)

#     gens_bus_idx  = get_gens_bus_num(gen_nodes)
    
#     Cg = sparse(gens_bus_idx, gens_bus_idx, ones(Int,length(gens_bus_idx)), num_nodes, num_nodes)

#     return Cg
# end


# function get_node_S(node)
 
#     lens_S=@optic _.S
 
#     return getall(node, lens_S)[1]
# end


# function get_nodes_S(nodes)
 
 
#     return map(get_node_S, collect(values(nodes)) )
# end


# function get_node_Vm(node)
 
#     lens_Vm=@optic _.Vm
 
#     return getall(node, lens_Vm)[1]
# end


# function get_nodes_Vm(nodes)
 
#     return map(get_node_Vm, collect(values(nodes)) )
# end


# function get_node_Vθ(node)
 
#     lens_Vθ=@optic _.Vθ
 
#     return getall(node, lens_Vθ)[1]
# end


# function get_nodes_Vθ(nodes)
 
#     return map(get_node_Vθ, collect(values(nodes)) )
# end


# function get_nodes_Qmax(nodes)
 
#     return map(get_node_Qmax, collect(values(nodes)) )
# end


# function get_nodes_Qmin(nodes)
 
#     return map(get_node_Qmin, collect(values(nodes)) )
    
# end


# function get_slack_and_gen_Sg(nodes)

#     slack_and_gen_nodes = get_slack_and_gen_nodes(nodes)
        
#     return map(get_node_S, slack_and_gen_nodes)
# end


# function get_Sg(gen_nodes, num_nodes)

#     gens_bus_idx  = get_gens_bus_num(gen_nodes)
        
#     gens_power   = map(get_node_S, collect(values(gen_nodes)) )

#     return sparsevec(gens_bus_idx, gens_power, num_nodes)
# end


# function node_Idx_from_name(node_name::String)

#     return parse(Int, split(lowercase(node_name),"bus")[2] )
# end


# function node_Idx_from_name(node)

#     node_name = get_node_Bus_name(node)

#     return parse(Int, split(lowercase(node_name),"bus")[2] )
# end


# """
# Get indices in slack nodes, gen nodes and load nodes order.
# """
# function get_Idx_slack_gens_loads_order(
#     slack_nodes, gen_nodes, load_nodes)

#     names_nodes =
#         unique([get_nodes_Bus_name(Slack_nodes);
#                 get_nodes_Bus_name(Gen_nodes);
#                 get_nodes_Bus_name(Load_nodes)])

#     nodes_Idx =
#         node_Idx_from_name.(names_nodes)
    
#     return nodes_Idx
# end


# function get_nodes_type_dims(
#     slack_nodes, gen_nodes, load_nodes, branches)

#     total_buses_type =
#         get_network_nodal_size(branches)
    
#     dim_slacks_buses =
#         length(collect(values(slack_nodes)))
    
#     dim_gens_buses =
#         length(collect(values(gen_nodes)))
    
#     dim_loads_buses =
#         total_buses_type -
#         (dim_slacks_buses + dim_gens_buses)
   
#     return (dim_slacks_buses,
#             dim_gens_buses,
#             dim_loads_buses)
# end


# function get_nodes_type_Idx(
#     slack_nodes,
#     gen_nodes,
#     load_nodes,
#     branches)

#     (dim_slacks_buses,
#      dim_gens_buses,
#      dim_loads_buses) =
#          get_nodes_type_dims(
#              slack_nodes,
#              gen_nodes,
#              load_nodes,
#              branches)

#     dims = [dim_slacks_buses,
#             dim_gens_buses,
#             dim_loads_buses]
    
#     offset = create_offsets(dims; counter=0)
    
#     Idx = create_idxs(offset, dims)
#     slacks_Idx = Idx[1]
#     gens_Idx   = Idx[2]
#     loads_Idx  = Idx[3]

#     gens_to_load_Idx =
#         first(gens_Idx):last(loads_Idx)

#     return (slacks_Idx,
#             gens_Idx,
#             loads_Idx,
#             gens_to_load_Idx)
# end

# function get_Idx_permutation(nodes_idx)

#     # sort a list of nodes_idx
#     idx_perm = sortperm(nodes_idxs) 

#     # get sorted list
#     sorted_nodes_idx = nodes_idx[idx_perm] 
    
#     perm_idx = Permutation(sorted_nodes_idx)

#     perm_matrix = sparse(Matrix(perm_idx))
    
#     return perm_matrix, perm_idx
# end

# # ------------------------------------------------------
# # Similar to pypower
# # ------------------------------------------------------

# function ΘV_to_VrVi(x, p)

#     # (PV_Idx, PQ_Idx), (PV_Θ_Idx, PQ_Θ_Idx, PV_V_Idx, PQ_V_Idx) = p
    
#     p2, p1 = p

#     V = ΘV_to_V(x, p1)   

#     return V_to_VrVi(V, p2)
    
# end

# function ur_ui_to_u(ur_ui)

#     return ur_ui[1] + im * ur_ui[2]
    
# end


# function xr_xi_to_x(ur_ui)

#     return ur_ui[1] + im * ur_ui[2]
    
# end



# function u_from_ur_ui(ur_ui)

#     return ur_ui[1] + im * ur_ui[2]
    
# end


# function x_from_xr_xi(ur_ui)

#     if length(ur_ui) == 1 && typeof(ur_ui[1]) == Float64
#         return ur_ui[1] + im * 0.0
        
#     elseif length(ur_ui) == 1 && typeof(ur_ui[1]) == ComplexF64
#         return  0.0 + im * ur_ui[1]

#     else

#         return ur_ui[1] + im * ur_ui[2]
#     end
    
# end


# function conj_x_from_xr_xi(ur_ui)

#     if length(ur_ui) == 1 && typeof(ur_ui[1]) == Float64
#         return ur_ui[1] - im * 0.0
        
#     elseif length(ur_ui) == 1 && typeof(ur_ui[1]) == ComplexF64
#         return  0.0 - im * ur_ui[1]

#     else

#         return ur_ui[1] - im * ur_ui[2]
#     end
    
# end


# function u_to_ΘV(u)

#     return [ angle.(u)...; abs.(u)... ]
    
# end


# function u_to_VΘ(u)

#     return [ abs(u), angle(u) ]
    
# end

# function VΘ_to_u(VΘ)

#     return VΘ[1] * exp(im * VΘ[2])
    
# end


# function ur_ui_to_ΘV(ur_ui)

#     u = ur_ui_to_u.(ur_ui)

#     # Θ = angle.(u) 
#     # V = abs.(u)

#     return [angle.(u)...;abs.(u)...]
    
# end


# function VrVi_to_ΘV(x, p)

#     # ((PV_Idx, PQ_Idx)), (PV_Vr_Idx, PQ_Vr_Idx, PV_Vi_Idx, PQ_Vi_Idx) = p

#     p2, p1 = p

#     V = VrVi_to_V(x, p1)

#     return V_to_ΘV(V, p2)
    
# end


# function VrVi_to_V(x, p)
    
#     PV_Vr_Idx, PQ_Vr_Idx, PV_Vi_Idx, PQ_Vi_Idx = p
    
#     x_Vr_Idx = first(PV_Vr_Idx):last(PQ_Vr_Idx)
#     x_Vi_Idx = first(PV_Vi_Idx):last(PQ_Vi_Idx)
    
#     x_Vr = x[x_Vr_Idx]
#     x_Vi = x[x_Vi_Idx]
    
#     return x_Vr .+ (x_Vi * im)
    
# end


# function V_to_VrVi(x, (PV_Idx, PQ_Idx))
    
#     # PV_Idx, PQ_Idx = p
    
#     PV_Vr = real.(x[PV_Idx])
#     PV_Vi = imag.(x[PV_Idx])

#     PQ_Vr = real.(x[PQ_Idx])
#     PQ_Vi = imag.(x[PQ_Idx])

#     return  [PV_Vr...;PQ_Vr...;PV_Vi...;PQ_Vi...]
    
# end

# function V_to_VrVi(V)

#     return  [[[real(v), imag(v)] for v in V]...]
    
# end



# function ΘV_to_V(x, p)
    
#     PV_Θ_Idx, PQ_Θ_Idx, PV_V_Idx, PQ_V_Idx = p

#     PV_Θ = x[PV_Θ_Idx]
#     PQ_Θ = x[PQ_Θ_Idx]
#     PV_V = x[PV_V_Idx]
#     PQ_V = x[PQ_V_Idx]

#     return [[PV_V .* exp.(im * PV_Θ)]...;
#             [PQ_V .* exp.(im * PQ_Θ)]...]
    
# end


# function V_to_ΘV(V, (PV_Idx, PQ_Idx) )
    
#     # PV_Idx, PQ_Idx = p

#     PV_Θ = angle.(V[PV_Idx])
#     PV_V = abs.(V[PV_Idx])

#     PQ_Θ = angle.(V[PQ_Idx])
#     PQ_V = abs.(V[PQ_Idx])

#     return  [PV_Θ...;PQ_Θ...;PV_V...;PQ_V...]
    
# end


# x_to_V(x, x_Θ_Idx, x_v_Idx) =
#     x[x_v_Idx] .* exp.(im * x[x_Θ_Idx] )

# V_to_x(V) = [angle.(V)...;abs.(V)...]

# function x_to_V(x)
    
#     n_x_Θ = n_x_v = Int64(length(x)/2 )

#     dims       = [n_x_Θ, n_x_v]
#     offset     = create_offsets(dims; counter=0)
#     Idx        = create_idxs(offset, dims)
#     x_Θ_Idx    = Idx[1]
#     x_v_Idx    = Idx[2]

#     return x[x_v_Idx] .* exp.(im * x[x_Θ_Idx] )
    
# end


# function polar_to_cartesian( v_θ )

#     v, θ = v_θ
    
#     return [v * cos(θ), v * sin(θ)]
    
# end


# function cartesian_to_polar(ur_ui)

#     ur, ui = ur_ui
    
#     return [abs(ur + im * ui),  angle( ur + im * ui )]
    
# end


######################################################
######################################################


function update_Jac!(
    Jac_view, gen_Q_view, p_update_Jac)
    
    # ------------------------------------------------------
    # See page 250, Millano, item 2 : the row corresponding to qh and the
    # column corresponding to vh in the Jacobian matrix gy are zero Only the
    # diagonal element has to be set to 1 to avoid singularity.
    # ------------------------------------------------------

    (Jac_Q_PV_row_Idx,
     x_v_PV_Idx,
     gen_Qmax,
     gen_Qmin) =
         p_update_Jac
    
    Q_nolimit_boolean =
        map((x) -> no_limit_violation(x...),
            zip(gen_Q_view, gen_Qmax, gen_Qmin))

    Q_limit_boolean   =
        map((x) -> limit_violation(x...),
            zip(gen_Q_view, gen_Qmax, gen_Qmin))

    Q_nolimit_violation = gen_Q_view[Q_nolimit_boolean]
    
    Q_limit_violation = gen_Q_view[Q_limit_boolean]
    
    if length(Q_limit_violation) == 0

        x_v_PV_restriction_Idx =
            x_v_PV_Idx[Q_nolimit_boolean]
        
        PV_Q_restriction_Idx   =
            Jac_Q_PV_row_Idx[Q_nolimit_boolean]

        for (qh_Idx, vh_Idx) in
            zip(PV_Q_restriction_Idx,
                x_v_PV_restriction_Idx )
            
            Jac_view[qh_Idx, :     ] .= 0.0
            Jac_view[:,      vh_Idx] .= 0.0
            Jac_view[qh_Idx, vh_Idx] = 1.0
            
        end

    end

    return nothing
end


function update_gen_reactive_power!(
    Sg_view, gens_bus_num, gen_Q_view)
    
    for (ind, gens_bus_Idx) in enumerate(gens_bus_num)
        real_part = real(Sg_view[gens_bus_Idx])
    Sg_view[gens_bus_Idx] = real_part + im *  gen_Q_view[ind]
    end
    return nothing
end


function update_gens_Q!(
    Sg_view, gen_Q_view, p_update_gens_Q)
    
    # ------------------------------------------------------
    # See page 250, Millano, item 2 : the row corresponding to qh and the
    # column corresponding to vh in the Jacobian matrix gy are zeroOnly the
    # diagonal element has to be set to 1 to avoid singularity.
    # ------------------------------------------------------

    gens_bus_num, gen_Qmax, gen_Qmin = p_update_gens_Q
    
    Q_nolimit_boolean =
        map((x) -> no_limit_violation(x...),
            zip(gen_Q_view, gen_Qmax, gen_Qmin))
    
    Q_limit_boolean =
        map((x) -> limit_violation(x...),
            zip(gen_Q_view, gen_Qmax, gen_Qmin))

    Q_nolimit_violation  =
        gen_Q_view[Q_nolimit_boolean]

    Q_limit_violation =
        gen_Q_view[Q_limit_boolean]

    if length(Q_limit_violation) != 0
        
        # gen_Q_update = threshold_limits.(gen_Q_view, gen_Qmax, gen_Qmin)
        # update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_update)
        # # or use gen_Q_view for updates
        
        gen_Q_view .=
            threshold_limits.(
                gen_Q_view, gen_Qmax, gen_Qmin)
        
        update_gen_reactive_power!(
            Sg_view, gens_bus_num, gen_Q_view)
    else
        update_gen_reactive_power!(
            Sg_view, gens_bus_num, gen_Q_view)
    end
        
    return nothing
end


# ------------------------------------------------------
# Algorithm from : William Stevenson, McGraw-Hill, 1994.
# page 347 to 355
# ------------------------------------------------------

P_i(i, V, Y) =
    (abs(V[i]))^2 * real(Y[i,i]) +
    sum([abs( V[i] * V[n] * Y[i,n]) *
    cos( angle(Y[i,n]) + angle(V[n]) - angle(V[i]) )
         for n in 1:length(V) if i != n ] )


Q_i(i, V, Y) =
    -(abs(V[i]))^2 * imag(Y[i,i]) -
    sum([abs( V[i] * V[n] * Y[i,n]) *
    sin( angle(Y[i,n]) + angle(V[n]) - angle(V[i]) )
         for n in 1:length(V) if i != n ] )


Sschd_i(i) = (Cg * Sg - Sd)[i]

Pschd_i(i) = real(Sschd_i(i))

Qschd_i(i) = imag(Sschd_i(i))

ΔP_i(i, V, Y)  = Pschd_i(i) - P_i(i, V, Y)

ΔQ_i(i, V, Y)  = Qschd_i(i) - Q_i(i, V, Y)

# J11

# off diagonal

∂P_i_∂δ_j(i,j,V,Y) = -abs( V[i] * V[j] * Y[i,j]) *
    sin( angle(Y[i,j]) + angle(V[j]) - angle(V[i]) )

# diagonal

# ∂P_i_∂δ_i(i,V,Y) = -sum([∂P_i_∂δ_j(i,n,V,Y) for n in 1:length(V) if i != n])

∂P_i_∂δ_i(i,V,Y) =
    -Q_i(i, V, Y) - (abs(V[i]))^2 * imag(Y[i,i])


# J21

# off diagonal

∂Q_i_∂δ_j(i,j,V,Y) = -abs(V[i] * V[j] * Y[i,j]) *
    cos(angle(Y[i,j]) + angle(V[j]) - angle(V[i]) )


# diagonal

# ∂Q_i_∂δ_i(i,V,Y) = -sum([∂Q_i_∂δ_j(i,n,V,Y) for n in 1:length(V) if i != n ])

∂Q_i_∂δ_i(i,V,Y) = P_i(i, V, Y) - (abs(V[i]))^2 * real(Y[i,i])

# J12

# off diagonal

# ∂P_i_∂Vm_j(i,j,V,Y) = abs( V[i] * Y[i,j]) *
#     cos(angle(Y[i,j]) + angle(V[j]) - angle(V[i]) )

# Vm_j_∂P_i_∂Vm_j(i,j,V,Y) = abs( V[j] * V[i] * Y[i,j]) *
#     cos( angle(Y[i,j]) + angle(V[j]) - angle(V[i]) )

Vm_j_∂P_i_∂Vm_j(i,j,V,Y) = -∂Q_i_∂δ_j(i,j,V,Y)

# diagonal

# ∂P_i_∂Vm_i(i,V,Y) = 2*abs(V[i]) * realY[i,i] + sum([ abs(V[n] * Y[i,n]) *
#     cos(angle(Y[i,n]) + angle(V[n]) - angle(V[i]) ) for n in 1:length(V) if i != n])

# Vm_i_∂P_i_∂Vm_i(i,V,Y) = ∂Q_i_∂δ_i(i,V,Y) + 2 * (abs(V[i]))^2 * real(Y[i,i])

Vm_i_∂P_i_∂Vm_i(i,V,Y) =
    P_i(i, V, Y) + (abs(V[i]))^2 * real(Y[i,i])


# J22

# off diagonal 

∂Q_i_Vm_j(i,j,V,Y) = -abs( V[i]  * Y[i,j]) *
    sin(angle(Y[i,j]) + angle(V[j]) - angle(V[i]) )

 # Vm_j_∂Q_i_∂δ_j(i,j,V,Y) = -abs(V[j] * V[i]  * Y[i,j]) *
 #    sin(angle(Y[i,j]) + angle(V[j]) - angle(V[i]) )

Vm_j_∂Q_i_Vm_j(i,j,V,Y) = ∂P_i_∂δ_j(i,j,V,Y) 

# diagonal

# Vm_i_∂Q_i_∂δ_i(i,V,Y) = -∂P_i_∂δ_i(i,V,Y)  - 2*(abs(V[i]))^2 * imag(Y[i,i])

Vm_i_∂Q_i_Vm_i(i,V,Y) = Q_i(i, V, Y) - (abs(V[i]))^2 * imag(Y[i,i])



# ------------------------------------------------------
# Current injections
# ------------------------------------------------------

function nodal_current_injection(V, p_ΔI)

    perm_matrix, Cg, Sd, Y = p_ΔI

    nodal_current  = Y * V + Diagonal(inv.(conj.(V))) *
        ( perm_matrix' * conj.( Sd ) )

    return perm_matrix * nodal_current  

end


function net_current_injection(V, p_ΔI)

    perm_matrix, Cg, Sd, Y = p_ΔI

    net_current  = Y * V 

    return perm_matrix * net_current  

end



function nodal_current_injection(
    x, Sd_view,
    (slack_Vm, slack_Vθ,
     gen_Vm, p_ΔI, p_Idx, p_x_Θ_Idx))

    perm_matrix, Cg, Sd, Y = p_ΔI

    PV_Idx, PQ_Idx = p_Idx

    x_Θ_PV_Idx, x_Θ_PQ_Idx, x_v_PV_Idx, x_v_PQ_Idx = p_x_Θ_Idx
    
    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    x_Θ_PV = x_reduced[x_Θ_PV_Idx]
    x_Θ_PQ = x_reduced[x_Θ_PQ_Idx]
    x_v_PV = x_reduced[x_v_PQ_Idx]
    x_v_PQ = x_reduced[x_v_PQ_Idx]

    x_v_PV .= gen_Vm

    V = [[slack_Vm .* exp.(im * slack_Vθ)]...;[x_v_PV .* exp.(im * x_Θ_PV)]...; [x_v_PQ .* exp.(im * x_Θ_PQ)]...]
    
    # mismatch  = Y * V - Diagonal(inv.(conj.(V))) * ( perm_matrix' * conj.(Sd_view) )

    nodal_current  = Y * V + Diagonal(inv.(conj.(V))) * ( perm_matrix' * conj.( Sd_view ) )

    return perm_matrix * nodal_current  

end



#------------------------------------------------------#
# ------------------------------------------------------
# Reduced
# ------------------------------------------------------
#------------------------------------------------------#

# ------------------------------------------------------
# Jacobian
# ------------------------------------------------------


function reduced_Jac_formula(V, Y, p_Idx )

    PV_Idx, PQ_Idx = p_Idx

    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    m11 = [[i==j ? ∂P_i_∂δ_i(i,V,Y) : ∂P_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

    M11 = reshape( [m11...;], size(Y) )

    n12 = [[i==j ? Vm_i_∂P_i_∂Vm_i(i,V,Y) : Vm_j_∂P_i_∂Vm_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

    N12 = reshape( [n12...;], size(Y) )

    n21 = [[i==j ? ∂Q_i_∂δ_i(i,V,Y) : ∂Q_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

    N21 = reshape( [n21...;], size(Y) )

    m22 = [[i==j ? Vm_i_∂Q_i_Vm_i(i,V,Y) : Vm_j_∂Q_i_Vm_j(i,j,V,Y) for j in 1:length(V)] for i in  1:length(V)] 

    M22 = reshape( [m22...;], size(Y) )  

    return vcat(hcat(M11'[PV_PQ_Idx, PV_PQ_Idx], N12'[PV_PQ_Idx, PQ_Idx]),
                hcat(N21'[PQ_Idx,    PV_PQ_Idx], M22'[PQ_Idx,    PQ_Idx]))
    
end


# ------------------------------------------------------
# Power mismathc
# ------------------------------------------------------


function reduced_ΔP_ΔQ(V, Sg_view, Sd_view, p_ΔP_ΔQ, p_Idx)

    perm_matrix, Cg, _, Y = p_ΔP_ΔQ

    PV_Idx, PQ_Idx = p_Idx
    
    PV_PQ_Idx      = first(PV_Idx):last(PQ_Idx)
    
    mismatch       = Diagonal(V) * conj.(Y * V) - perm_matrix' * (Cg * Sg_view - Sd_view)

    ΔP_ΔQ = Array( [[ real.(mismatch[PV_Idx]), real.(mismatch[PQ_Idx]), imag.(mismatch[PQ_Idx])]...; ] ) 

end

# ------------------------------------------------------
# Power flow solution algorithm
# ------------------------------------------------------


function reduced_pf_sol(x0,
                        (gen_Q_view, Sg_view, Sd_view), p;
                       maxiter=40,
                       ftol=1000*eps(),
                       xtol=1000*eps() )

    p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx , p_slack = p

    nodes_vmax, nodes_vmin, gen_Vm, p_Idx = p_V_limits

    gens_bus_num, gen_Qmax, gen_Qmin  = p_update_gens_Q

    perm_matrix, Cg, Sd, Y = p_ΔP_ΔQ 

    PV_Idx_p, PQ_Idx_p = p_Idx

    x_Θ_Idx, x_v_Idx, x_Θ_PV_PQ_Idx, x_v_PV_PQ_Idx, x_Θ_PV_Idx, x_v_PV_Idx, x_Θ_PQ_Idx, x_v_PQ_Idx, Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx, Δx_Θ_PV_Idx, Δx_v_PV_Idx, Δx_Θ_PQ_Idx, Δx_v_PQ_Idx = p_x_Θ_Δx_Idx

    x_k_1              = deepcopy(x0)

    x_k_1[x_v_PV_Idx] .= gen_Vm

    x_k                = deepcopy(x0)

    Δx                 = zeros(length( x_k ))

    x_k[x_v_PV_Idx]   .= gen_Vm

    V_k_1              = x_to_V(x_k_1, x_Θ_Idx, x_v_Idx)

    V_k                = x_to_V(x_k, x_Θ_Idx, x_v_Idx)

    ΔΘVm_dim           = sum([length(x_Θ_PV_PQ_Idx), length(x_v_PV_PQ_Idx)])
    
    ΔΘVm = Inf

    ΔP_ΔQ_k = reduced_ΔP_ΔQ( V_k_1, Sg_view, Sd_view, p_ΔP_ΔQ, p_Idx )

    local k  = 0
        
    while (norm(ΔΘVm) > xtol) && (norm(ΔP_ΔQ_k) > ftol)

        x_k_1[x_v_PV_Idx] .= gen_Vm

        V_k_1              .=  x_to_V(x_k_1,  x_Θ_Idx, x_v_Idx)

        Jac_k              = reduced_Jac_formula(V_k_1, Y,  p_Idx)

        ΔP_ΔQ_k            = reduced_ΔP_ΔQ(V_k_1, Sg_view, Sd_view, p_ΔP_ΔQ, p_Idx ) 

        ΔΘVm               =  -Jac_k \ ΔP_ΔQ_k

        Δx[Δx_Θ_PV_PQ_Idx] = ΔΘVm[Δx_Θ_PV_PQ_Idx]

        

        # if  reduced == true
            
        if length(ΔΘVm) < ΔΘVm_dim # length(Δx)
            
            Δx[Δx_v_PQ_Idx] .= ΔΘVm[1+length(Δx_Θ_PV_PQ_Idx):end]
            
        else
            
            Δx[Δx_v_PV_Idx]  .= ΔΘVm[Δx_v_PV_Idx]
            Δx[Δx_v_PQ_Idx]  .= ΔΘVm[Δx_v_PQ_Idx]
            
        end

        x_k[ x_Θ_PV_PQ_Idx ] .= x_k_1[ x_Θ_PV_PQ_Idx ] .+ Δx[ Δx_Θ_PV_PQ_Idx ]
        x_k[ x_v_PQ_Idx ]    .= x_k_1[ x_v_PQ_Idx ] .* (1 .+ Δx[ Δx_v_PQ_Idx ])

        # V_k =  x_to_V(x_k, x_Θ_Idx, x_v_Idx)

        V_k  .=  x_to_V(x_k, x_Θ_Idx, x_v_Idx)

        gen_Q = imag.((Diagonal(V_k) * conj.(Y * V_k) + perm_matrix' * Sd_view)[PV_Idx_p])

        update_gens_Q!(Sg_view, gen_Q, p_update_gens_Q)

        x_k_1 .= x_k

        k += 1
        
        if k == maxiter
            @warn "Maximum number of iterations reached."

            break
        end

    end

    return k, V_k, ΔP_ΔQ_k, Δx   
end

# ------------------------------------------------------
# ------------------------------------------------------


function get_dynamic_powerflow_sol(
    x0::Vector{Float64},
    (Sg_view, Sd_view),
    dyn_powerflow_sol_param )

    Ybus_noperm, nodes_name, branches_name, num_nodes, perm_matrix, gen_Q_view, _, _, p = dyn_powerflow_sol_param
    
    _, _, p_ΔP_ΔQ, _, p_x_Θ_Δx_Idx, p_slack = p

    _, Cg, Sd, Y = p_ΔP_ΔQ
    
    # ------------------------------------------------------
    # Powerflow
    # ------------------------------------------------------
        
    _, V_k, _, _ = reduced_pf_sol(x0, (gen_Q_view, Sg_view, Sd_view), p; maxiter=40, ftol=1000*eps(), xtol=1000*eps())
    
    V = perm_matrix * V_k

    Iinj = nodal_current  = Ybus_noperm * V + Diagonal(inv.(conj.(V))) * ( conj.( Sd_view ) )

    Inet_inj  = net_current  = Ybus_noperm * V 

    Sbus     = Diagonal(V) * conj.(Ybus_noperm * V)

    GenSinj  = Sbus + Sd_view

    return [ (vh, θh, ph, qh, ih_r, ih_i) for (vh, θh, ph, qh, ih_r, ih_i) in zip(abs.(V), angle.(V), real.(Sbus), imag.(Sbus), real.(Inet_inj), imag.(Inet_inj)) ]
end


function get_dynamic_powerflow_sol(
    V0::Vector{ComplexF64},
    (Sg_view, Sd_view),
    dyn_powerflow_sol_param)

    Ybus_noperm, nodes_name, branches_name, perm_matrix, gen_Q_view, _, _, p = dyn_powerflow_sol_param
    
    _, _, p_ΔP_ΔQ, _, p_x_Θ_Δx_Idx = p

    _, Cg, Sd, Y = p_ΔP_ΔQ
    
    # ------------------------------------------------------
    # Powerflow
    # ------------------------------------------------------

    x0 = [angle.(V0)...; abs(V0)...]
    
    _, V_k, _, _ = reduced_pf_sol(x0, (gen_Q_view, Sg_view, Sd_view), p; maxiter=40, ftol=1000*eps(), xtol=1000*eps())
    
    V = perm_matrix * V_k

    Iinj = nodal_current  = Ybus_noperm * V + Diagonal(inv.(conj.(V))) * ( conj.( Sd_view ) )

    Inet_inj  = net_current  = Ybus_noperm * V 

    Sbus     = Diagonal(V) * conj.(Ybus_noperm * V)

    GenSinj  = Sbus + Sd_view

    return [(vh, θh, ph, qh, ih_r, ih_i) for (vh, θh, ph, qh, ih_r, ih_i) in
                zip(abs.(V), angle.(V), real.(Sbus),
                     imag.(Sbus), real.(Inet_inj),
                     imag.(Inet_inj) )]

end


function get_dynamic_reduced_powerflow( case_data_func )

    Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

    perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

    num_nodes     = get_network_nodal_size(Branches)

    gens_S        = get_nodes_S(Gen_nodes)

    gen_Q         = imag(gens_S) 

    gens_bus_num  = get_gens_bus_num(Gen_nodes)

    Cg            = get_Cg(Gen_nodes, num_nodes)

    Sg            = get_Sg(Gen_nodes, num_nodes)

    gen_Vm        = get_nodes_Vm(Gen_nodes) 

    num_gen       =  length(gen_Vm)

    gen_Qmax      = get_nodes_Qmax(Gen_nodes)

    gen_Qmin      = get_nodes_Qmin(Gen_nodes)

    Sd            = get_loads_Sd(Load_nodes)


    # branches_name  = get_branches_name(Branches)   
    # nodes_name     = get_nodes_Bus_name(Nodes)

    branches_name  = collect(keys(Branches))
    
    nodes_name     = collect(keys(Nodes))      
    
    Yf             = get_Yf(Nodes, Branches)
    
    Yt             = get_Yt(Nodes, Branches)
    

    Ybus_noperm   = get_Ybus(Nodes, Branches, shunts=Shunts)

    slack_Vm      = get_nodes_Vm(Slack_nodes)

    slack_Vθ      = get_nodes_Vθ(Slack_nodes)

    nodes_vmax    = get_nodes_vmax(Nodes)

    nodes_vmax_p  = perm_matrix' * nodes_vmax

    nodes_vmin    = get_nodes_vmin(Nodes)

    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    = get_nodes_Bus_name(Nodes)

    Ybus_p        = perm_matrix' * (Ybus_noperm) * perm_matrix

    # ------------------------------------------------------
    # Indices
    # ------------------------------------------------------


    SK_Idx_p      = idx_slack_nodes

    PV_Idx_p      = idx_source_nodes   # gen nodes indices from perm 

    PQ_Idx_p      = idx_demand_nodes   # load nodes indices from perm

    PV_PQ_Idx_p   = first(PV_Idx_p):last(PQ_Idx_p)  

    SK_PV_PQ_V_dim     = [length(SK_Idx_p),length(PV_Idx_p),length(PQ_Idx_p)]
    SK_PV_PQ_V_offset  = create_offsets(SK_PV_PQ_V_dim)
    SK_PV_PQ_V_Idx     = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

    SK_V_Idx           = SK_PV_PQ_V_Idx[1]
    PV_V_Idx           = SK_PV_PQ_V_Idx[2]
    PQ_V_Idx           = SK_PV_PQ_V_Idx[3]

    # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

    x_dim      = [length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx),
                  length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx)]

    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    x_Θ_SK_Idx = x_Idx[1]
    x_Θ_PV_Idx = x_Idx[2]
    x_Θ_PQ_Idx = x_Idx[3]
    
    x_v_SK_Idx = x_Idx[4]
    x_v_PV_Idx = x_Idx[5]
    x_v_PQ_Idx = x_Idx[6]

    x_Θ_Idx       = first(x_Θ_SK_Idx):last(x_Θ_PQ_Idx)
    x_v_Idx       = first(x_v_SK_Idx):last(x_v_PQ_Idx)

    x_Θ_PV_PQ_Idx = first(x_Θ_PV_Idx):last(x_Θ_PQ_Idx)
    x_v_PV_PQ_Idx = first(x_v_PV_Idx):last(x_v_PQ_Idx)

    Δx_dim = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PV_Idx),length(x_v_PQ_Idx)]

    Δx_offset   = create_offsets( Δx_dim )
    Δx_Idx      = create_idxs( Δx_offset, Δx_dim )

    Δx_Θ_PV_Idx = Δx_Idx[1]
    Δx_Θ_PQ_Idx = Δx_Idx[2]

    Δx_v_PV_Idx = Δx_Idx[3]
    Δx_v_PQ_Idx = Δx_Idx[4]

    Δx_Θ_PV_PQ_Idx = first(Δx_Θ_PV_Idx):last(Δx_Θ_PQ_Idx)

    Δx_v_PV_PQ_Idx = first(Δx_v_PV_Idx):last(Δx_v_PQ_Idx)

    # ------------------------------------------------------

    # computing arrays, indices and views
    # ------------------------------------------------------

    gen_Q      = zeros(Float64, num_gen)

    gen_Q_view =  view(gen_Q, 1:num_gen)

    Sg_view    = @view Sg[:]

    update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

    Sd_view    = @view Sd[:]
    
    Y          = Ybus_p

    # ------------------------------------------------------
    # Parameteres
    # ------------------------------------------------------
    
    p_ΔP_ΔQ            = ( perm_matrix, Cg, Sd, Y )

    p_update_gens_Q    = ( gens_bus_num, gen_Qmax, gen_Qmin )

    p_x_Θ_Δx_Idx       = ( x_Θ_Idx,        x_v_Idx,
                           x_Θ_PV_PQ_Idx,  x_v_PV_PQ_Idx,
                           x_Θ_PV_Idx,     x_v_PV_Idx,
                           x_Θ_PQ_Idx,     x_v_PQ_Idx,
                           
                           Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx,
                           Δx_Θ_PV_Idx,    Δx_v_PV_Idx,
                           Δx_Θ_PQ_Idx,    Δx_v_PQ_Idx )

    p_Idx = ( PV_Idx_p, PQ_Idx_p )
    
    p_V_limits = ( nodes_vmax, nodes_vmin, gen_Vm, p_Idx )

    p_slack = ( slack_Vm, slack_Vθ, x_Θ_SK_Idx, x_v_SK_Idx) 

    p = ( p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack )

     dyn_powerflow_sol_param = (; Ybus_noperm, nodes_name, branches_name, num_nodes, perm_matrix, gen_Q_view, Sg_view, Sd_view, p)
    
    # ------------------------------------------------------
    # Initial condition
    # ------------------------------------------------------

    x_Θ                = zeros(num_nodes)

    x_v                = ones(num_nodes)

    x0                 = [x_Θ...; x_v...]

    # ------------------------------------------------------
    # Powerflow
    # ------------------------------------------------------
    
    
    # k, V_k, ΔP_ΔQ_k, Δx = reduced_pf_sol(x0,
    #                                      (gen_Q_view, Sg_view), p;
    #                                      maxiter=40,
    #                                      ftol=1000*eps(),
    #                                      xtol=1000*eps())
    
    k, V_k, ΔP_ΔQ_k, Δx = reduced_pf_sol(x0, (gen_Q_view, Sg_view, Sd_view), p; maxiter=40, ftol=1000*eps(), xtol=1000*eps())
    
    # ------------------------------------------------------
    # Results
    # ------------------------------------------------------

    Inet_inj  = net_current_injection(V_k,  p_ΔP_ΔQ)

    Iinj      = nodal_current_injection(V_k, p_ΔP_ΔQ)
    
    V         = perm_matrix * V_k

    Sbus      = Diagonal(V) * conj.(Ybus_noperm * V)

    GenSinj   = Diagonal(V) * conj.(Ybus_noperm * V) + Sd_view
    
    If         = Yf * V
    
    It         = Yt * V

    Ibranches  = If + It

    bus_dict_Iinj = OrderedDict( name => [ih_r, ih_i] for (name, ih_r, ih_i) in zip(nodes_name, round.(real.(Iinj); digits=4), round.(imag.(Iinj); digits=4)) )    

    branch_dict_init = OrderedDict( name => ( real(i_f), imag(i_f), real(i_t), imag(i_t) ) for (name, i_f, i_t) in zip(branches_name, If, It ) )

    # bus_dict_init = OrderedDict(name =>
    #     (vh, θh, ph, qh, ih_r, ih_i ) for (name, vh, θh, ph, qh, ih_r, ih_i) in zip(nodes_name,
    #                                      round.(abs.(V);        digits=4),
    #                                      round.(angle.(V);      digits=4),
    #                                      round.(real.(Sbus);    digits=4),
    #                                      round.(imag.(Sbus);    digits=4),
    #                                      round.(real.(Inet_inj);digits=4),
    #                                      round.(imag.(Inet_inj);digits=4) ))

    bus_dict_init = OrderedDict(name => (vh, θh, ph, qh, ih_r, ih_i ) for (name, vh, θh, ph, qh, ih_r, ih_i) in zip(nodes_name, round.(abs.(V); digits=4), round.(angle.(V); digits=4), round.(real.(Sbus); digits=4), round.(imag.(Sbus); digits=4), round.(real.(Inet_inj);digits=4), round.(imag.(Inet_inj);digits=4) ))     
        
    return perm_matrix, ((Sg_view, Sd_view), dyn_powerflow_sol_param), Dict("Iinj" => Iinj, "Inet_inj" => Inet_inj, "Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V), "Vθ" => angle.(V), "Vbus" => V, "Ibranches" => Ibranches, "Ybus" => Ybus_noperm , "Sbus" => Sbus, "GenSinj" => GenSinj, "dict_init" =>Dict("bus_dict_init" => bus_dict_init,  "branch_dict_init" =>  branch_dict_init, "bus_dict_Iinj" => bus_dict_Iinj) )     
    
end
                
#-------------------------------------------------------

function reduced_powerflow( case_data_func )

    Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

    perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

    num_nodes     = get_network_nodal_size(Branches)

    gens_S        = get_nodes_S(Gen_nodes)

    gen_Q         = imag(gens_S) 

    gens_bus_num  = get_gens_bus_num(Gen_nodes)

    Cg            = get_Cg(Gen_nodes, num_nodes)

    Sg            = get_Sg(Gen_nodes, num_nodes)

    gen_Vm        = get_nodes_Vm(Gen_nodes) 

    num_gen       =  length(gen_Vm)

    gen_Qmax      = get_nodes_Qmax(Gen_nodes)

    gen_Qmin      = get_nodes_Qmin(Gen_nodes)

    Sd            = get_loads_Sd(Load_nodes)


    # branches_name  = get_branches_name(Branches)   
    # nodes_name     = get_nodes_Bus_name(Nodes)

    branches_name  = collect(keys(Branches))
    
    nodes_name     = collect(keys(Nodes))      
    
    Yf             = get_Yf(Nodes, Branches)
    
    Yt             = get_Yt(Nodes, Branches)
    
    Ybus_noperm   = get_Ybus(Nodes, Branches, shunts=Shunts)

    slack_Vm      = get_nodes_Vm(Slack_nodes)

    slack_Vθ      = get_nodes_Vθ(Slack_nodes)

    nodes_vmax    = get_nodes_vmax(Nodes)

    nodes_vmax_p  = perm_matrix' * nodes_vmax

    nodes_vmin    = get_nodes_vmin(Nodes)

    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    = get_nodes_Bus_name(Nodes)

    Ybus_p        = perm_matrix' * (Ybus_noperm) * perm_matrix

    # ------------------------------------------------------
    # Indices
    # ------------------------------------------------------


    SK_Idx_p      = idx_slack_nodes

    PV_Idx_p      = idx_source_nodes   # gen nodes indices from perm 

    PQ_Idx_p      = idx_demand_nodes   # load nodes indices from perm

    PV_PQ_Idx_p   = first(PV_Idx_p):last(PQ_Idx_p)  

    SK_PV_PQ_V_dim     = [length(SK_Idx_p),length(PV_Idx_p),length(PQ_Idx_p)]
    SK_PV_PQ_V_offset  = create_offsets(SK_PV_PQ_V_dim)
    SK_PV_PQ_V_Idx     = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

    SK_V_Idx           = SK_PV_PQ_V_Idx[1]
    PV_V_Idx           = SK_PV_PQ_V_Idx[2]
    PQ_V_Idx           = SK_PV_PQ_V_Idx[3]

    # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

    x_dim      = [length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx),
                  length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx)]

    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    x_Θ_SK_Idx = x_Idx[1]
    x_Θ_PV_Idx = x_Idx[2]
    x_Θ_PQ_Idx = x_Idx[3]
    
    x_v_SK_Idx = x_Idx[4]
    x_v_PV_Idx = x_Idx[5]
    x_v_PQ_Idx = x_Idx[6]

    x_Θ_Idx       = first(x_Θ_SK_Idx):last(x_Θ_PQ_Idx)
    x_v_Idx       = first(x_v_SK_Idx):last(x_v_PQ_Idx)

    x_Θ_PV_PQ_Idx = first(x_Θ_PV_Idx):last(x_Θ_PQ_Idx)
    x_v_PV_PQ_Idx = first(x_v_PV_Idx):last(x_v_PQ_Idx)

    Δx_dim = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PV_Idx),length(x_v_PQ_Idx)]

    Δx_offset   = create_offsets( Δx_dim )
    Δx_Idx      = create_idxs( Δx_offset, Δx_dim )

    Δx_Θ_PV_Idx = Δx_Idx[1]
    Δx_Θ_PQ_Idx = Δx_Idx[2]

    Δx_v_PV_Idx = Δx_Idx[3]
    Δx_v_PQ_Idx = Δx_Idx[4]

    Δx_Θ_PV_PQ_Idx = first(Δx_Θ_PV_Idx):last(Δx_Θ_PQ_Idx)

    Δx_v_PV_PQ_Idx = first(Δx_v_PV_Idx):last(Δx_v_PQ_Idx)

    # ------------------------------------------------------

    # computing arrays, indices and views
    # ------------------------------------------------------

    gen_Q      = zeros(Float64, num_gen)

    gen_Q_view =  view(gen_Q, 1:num_gen)

    Sg_view    = @view Sg[:]

    update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

    Sd_view    = @view Sd[:]
    
    Y          = Ybus_p

    # ------------------------------------------------------
    # Parameteres
    # ------------------------------------------------------
    
    p_ΔP_ΔQ            = ( perm_matrix, Cg, Sd, Y )

    p_update_gens_Q    = ( gens_bus_num, gen_Qmax, gen_Qmin )

    p_x_Θ_Δx_Idx       = ( x_Θ_Idx,        x_v_Idx,
                           x_Θ_PV_PQ_Idx,  x_v_PV_PQ_Idx,
                           x_Θ_PV_Idx,     x_v_PV_Idx,
                           x_Θ_PQ_Idx,     x_v_PQ_Idx,
                           
                           Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx,
                           Δx_Θ_PV_Idx,    Δx_v_PV_Idx,
                           Δx_Θ_PQ_Idx,    Δx_v_PQ_Idx )

    p_Idx              = (PV_Idx_p, PQ_Idx_p)
    
    p_V_limits         = ( nodes_vmax, nodes_vmin, gen_Vm, p_Idx )


    p_slack =  (slack_Vm, slack_Vθ, x_Θ_SK_Idx, x_v_SK_Idx) 
    
    p                  = ( p_V_limits, p_update_gens_Q, p_ΔP_ΔQ,
                           p_Idx, p_x_Θ_Δx_Idx, p_slack )

    # ------------------------------------------------------
    # Initial condition
    # ------------------------------------------------------

    x_Θ                = zeros(num_nodes)

    x_v                = ones(num_nodes)

    x0                 = [x_Θ...; x_v...]

    # ------------------------------------------------------
    # Powerflow
    # ------------------------------------------------------
    
    
    # k, V_k, ΔP_ΔQ_k, Δx = reduced_pf_sol(x0,
    #                                      (gen_Q_view, Sg_view), p;
    #                                      maxiter=40,
    #                                      ftol=1000*eps(),
    #                                      xtol=1000*eps())

    
    k, V_k, ΔP_ΔQ_k, Δx = reduced_pf_sol(x0,
                                         (gen_Q_view, Sg_view, Sd_view), p;
                                         maxiter=40,
                                         ftol=1000*eps(),
                                         xtol=1000*eps())
    
    # ------------------------------------------------------
    # Results
    # ------------------------------------------------------

    Inet_inj  = net_current_injection(V_k,  p_ΔP_ΔQ)

    Iinj      = nodal_current_injection(V_k, p_ΔP_ΔQ)
    
    V         = perm_matrix * V_k

    Sbus = Diagonal(V) * conj.(Ybus_noperm * V)

    GenSinj = Diagonal(V) * conj.(Ybus_noperm * V) + Sd_view
    
    If         = Yf * V
    
    It         = Yt * V

    Ibranches  = If + It

    bus_dict_Iinj = OrderedDict(name =>
        [ih_r, ih_i] for (name, ih_r, ih_i) in zip(nodes_name, round.(real.(Iinj); digits=4), round.(imag.(Iinj); digits=4) ))    

    branch_dict_init = OrderedDict(name =>
        ( real(i_f), imag(i_f), real(i_t), imag(i_t) ) for (name, i_f, i_t) in zip(branches_name, If, It))

    # bus_dict_init = OrderedDict(name =>
    #     (vh, θh, ph, qh, ih_r, ih_i ) for (name, vh, θh, ph, qh, ih_r, ih_i) in zip(nodes_name,
    #                                      round.(abs.(V);        digits=4),
    #                                      round.(angle.(V);      digits=4),
    #                                      round.(real.(Sbus);    digits=4),
    #                                      round.(imag.(Sbus);    digits=4),
    #                                      round.(real.(Inet_inj);digits=4),
    #                                      round.(imag.(Inet_inj);digits=4) ))   
    
    bus_dict_init = OrderedDict(name =>
        (vh, θh, ph, qh, ih_r, ih_i ) for (name, vh, θh, ph, qh, ih_r, ih_i) in zip(nodes_name,
                                         round.(abs.(V);        digits=4),
                                         round.(angle.(V);      digits=4),
                                         round.(real.(Sbus);    digits=4),
                                         round.(imag.(Sbus);    digits=4),
                                         round.(real.(Inet_inj);digits=4),
                                         round.(imag.(Inet_inj);digits=4) ))      

    # return perm_matrix, k, V_k, ΔP_ΔQ_k, Δx
    
    return Dict("Iinj" => Iinj, "Inet_inj" => Inet_inj, "Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V), "Vθ" => angle.(V), "Vbus" => V, "Ibranches" => Ibranches, "Ybus" => Ybus_noperm , "Sbus" => Sbus, "GenSinj" => GenSinj, "dict_init" =>Dict("bus_dict_init" => bus_dict_init,  "branch_dict_init" =>  branch_dict_init, "bus_dict_Iinj" => bus_dict_Iinj) )     
    
end

#-------------------------------------------------------
##########################################
#-------------------------------------------------------

"""
x_θ_v_red_view = @view x[[collect(x_Θ_PV_PQ_Idx); collect(x_v_PQ_Idx)]]

"""
function x_reduced_ΔP_ΔQ( x_θ_v_red_view, (x, Sg_view, Sd_view, p))

    _, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack = p
    
    perm_matrix, Cg, _, Y = p_ΔP_ΔQ

    PV_Idx, PQ_Idx = p_Idx
    
    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    x_Θ_Idx, x_v_Idx, _, _, _, _, _, _, _, _, _, _, _, _ = p_x_Θ_Δx_Idx

    V = x_to_V(x, x_Θ_Idx, x_v_Idx)

    gen_Q = imag.((Diagonal(V) * conj.(Y * V) + perm_matrix' * Sd_view)[PV_Idx_p])

    update_gens_Q!(Sg_view, gen_Q, p_update_gens_Q)    
    
    mismatch  = Diagonal(V) * conj.(Y * V) - perm_matrix' * ( Cg * Sg_view - Sd_view )
    
    mismatch  = Diagonal(V) * conj.(Y * V) - perm_matrix' * ( Cg * Sg_view - Sd_view )

    return Array( [[real.(mismatch[PV_Idx]), real.(mismatch[PQ_Idx]), imag.(mismatch[PQ_Idx])]...; ] )
    

end



function x_reduced_Jac_formula(x_θ_v_red_view, (x, Y, p_Idx, p_x_Θ_Δx_Idx) )

    PV_Idx, PQ_Idx = p_Idx

    x_Θ_Idx, x_v_Idx, _, _, _, _, _, _, _, _, _, _, _, _ = p_x_Θ_Δx_Idx
    
    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    V = x_to_V(x, x_Θ_Idx, x_v_Idx)

    m11 = [[i==j ? ∂P_i_∂δ_i(i,V,Y) : ∂P_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

    M11 = reshape( [m11...;], size(Y) )

    n12 = [[i==j ? Vm_i_∂P_i_∂Vm_i(i,V,Y) : Vm_j_∂P_i_∂Vm_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

    N12 = reshape( [n12...;], size(Y) )

    n21 = [[i==j ? ∂Q_i_∂δ_i(i,V,Y) : ∂Q_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

    N21 = reshape( [n21...;], size(Y) )

    m22 = [[i==j ? Vm_i_∂Q_i_Vm_i(i,V,Y) : Vm_j_∂Q_i_Vm_j(i,j,V,Y) for j in 1:length(V)] for i in  1:length(V)] 

    M22 = reshape( [m22...;], size(Y) )  

    return vcat(hcat(M11'[PV_PQ_Idx, PV_PQ_Idx], N12'[PV_PQ_Idx, PQ_Idx]),
                hcat(N21'[PQ_Idx,    PV_PQ_Idx], M22'[PQ_Idx,    PQ_Idx]))
    
end



function x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param )

    x0_view, Sg_view, Sd_view, slack_Vm, slack_Vθ, Θ_PV_Θ_PQ_v_PQ_Idx, p = param
    
    p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack = p
    _, _, gen_Vm, _ = p_V_limits
    
    perm_matrix, Cg, _, Y = p_ΔP_ΔQ

    PV_Idx, PQ_Idx = p_Idx
    
    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    x_Θ_Idx, x_v_Idx, _, _, x_Θ_PV_Idx, _, x_Θ_PQ_Idx, x_v_PQ_Idx, _, _, _, _, _, _ = p_x_Θ_Δx_Idx

    # dims_view = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PQ_Idx)]

    # _,_,Θ_PV_Θ_PQ_v_PQ_Idx = create_size_offset_Idx(dims_view)

    Θ_PV = x_θ_v_red_view[Θ_PV_Θ_PQ_v_PQ_Idx[1]]
    Θ_PQ = x_θ_v_red_view[Θ_PV_Θ_PQ_v_PQ_Idx[2]]
    v_PQ = x_θ_v_red_view[Θ_PV_Θ_PQ_v_PQ_Idx[3]]

    V_PQ       = v_PQ  .* exp.(im * Θ_PQ)
    V_PV       = gen_Vm .* exp.(im * Θ_PV)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    
    V          = [V_slack...; V_PV...;V_PQ...]    
    

    x0_view    .= u_to_ΘV(V)
    
    # V = x_to_V(x, x_Θ_Idx, x_v_Idx)
    
    mismatch  = Diagonal(V) * conj.(Y * V) - perm_matrix' * ( Cg * Sg_view - Sd_view )
    
    mismatch  = Diagonal(V) * conj.(Y * V) - perm_matrix' * ( Cg * Sg_view - Sd_view )

    ΔP_ΔQ .= Array( [[real.(mismatch[PV_Idx]), real.(mismatch[PQ_Idx]), imag.(mismatch[PQ_Idx])]...; ] )


    gen_Q = imag.((Diagonal(V) * conj.(Y * V) + perm_matrix' * Sd_view)[PV_Idx])

    update_gens_Q!(Sg_view, gen_Q, p_update_gens_Q)    
    
    return nothing

end



function x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param )

    # x, _, _,_,_,_, p = param

     _, _, _, slack_Vm, slack_Vθ, Θ_PV_Θ_PQ_v_PQ_Idx, p = param
    
    # _, _, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx = p

    p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack = p
    _, _, gen_Vm, _ = p_V_limits    
    
    _, _, _, Y = p_ΔP_ΔQ
    
    PV_Idx, PQ_Idx = p_Idx

    x_Θ_Idx, x_v_Idx, _, _, _, _, _, _, _, _, _, _, _, _ = p_x_Θ_Δx_Idx
    
    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    Θ_PV = x_θ_v_red_view[Θ_PV_Θ_PQ_v_PQ_Idx[1]]
    Θ_PQ = x_θ_v_red_view[Θ_PV_Θ_PQ_v_PQ_Idx[2]]
    v_PQ = x_θ_v_red_view[Θ_PV_Θ_PQ_v_PQ_Idx[3]]

    V_PQ       = v_PQ  .* exp.(im * Θ_PQ)
    V_PV       = gen_Vm .* exp.(im * Θ_PV)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    
    V          = [V_slack...; V_PV...;V_PQ...]      

    # V = x_to_V(x, x_Θ_Idx, x_v_Idx)

    m11 = [[i==j ? ∂P_i_∂δ_i(i,V,Y) : ∂P_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

    M11 = reshape( [m11...;], size(Y) )

    n12 = [[i==j ? Vm_i_∂P_i_∂Vm_i(i,V,Y) : Vm_j_∂P_i_∂Vm_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

    N12 = reshape( [n12...;], size(Y) )

    n21 = [[i==j ? ∂Q_i_∂δ_i(i,V,Y) : ∂Q_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

    N21 = reshape( [n21...;], size(Y) )

    m22 = [[i==j ? Vm_i_∂Q_i_Vm_i(i,V,Y) : Vm_j_∂Q_i_Vm_j(i,j,V,Y) for j in 1:length(V)] for i in  1:length(V)] 

    M22 = reshape( [m22...;], size(Y) )  

    Jac_Matrix .= vcat(hcat(M11'[PV_PQ_Idx, PV_PQ_Idx], N12'[PV_PQ_Idx, PQ_Idx]),
                       hcat(N21'[PQ_Idx,    PV_PQ_Idx], M22'[PQ_Idx,    PQ_Idx]))

    return nothing
    
end

function x_fj!( ΔP_ΔQ, Jac_Matrix, x_θ_v_red_view, param )

    # # x, Sg_view, Sd_view, p = param

    # x, Sg_view, Sd_view, slack_Vm, slack_Vθ, Θ_PV_Θ_PQ_v_PQ_Idx, p = param

    # _, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx = p
    
    # perm_matrix, Cg, _, Y = p_ΔP_ΔQ

    # PV_Idx, PQ_Idx = p_Idx
    
    # PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    # x_Θ_Idx, x_v_Idx, _, _, _, _, _, _, _, _, _, _, _, _ = p_x_Θ_Δx_Idx

    # PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    # V = x_to_V(x, x_Θ_Idx, x_v_Idx)


    x0_view, Sg_view, Sd_view, slack_Vm, slack_Vθ, Θ_PV_Θ_PQ_v_PQ_Idx, p = param
    
    p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack = p
    _, _, gen_Vm, _ = p_V_limits
    
    perm_matrix, Cg, _, Y = p_ΔP_ΔQ

    PV_Idx, PQ_Idx = p_Idx
    
    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    x_Θ_Idx, x_v_Idx, _, _, x_Θ_PV_Idx, _, x_Θ_PQ_Idx, x_v_PQ_Idx, _, _, _, _, _, _ = p_x_Θ_Δx_Idx

    Θ_PV = x_θ_v_red_view[Θ_PV_Θ_PQ_v_PQ_Idx[1]]
    Θ_PQ = x_θ_v_red_view[Θ_PV_Θ_PQ_v_PQ_Idx[2]]
    v_PQ = x_θ_v_red_view[Θ_PV_Θ_PQ_v_PQ_Idx[3]]

    V_PQ       = v_PQ  .* exp.(im * Θ_PQ)
    V_PV       = gen_Vm .* exp.(im * Θ_PV)
    V_slack    = slack_Vm .* exp.(im * slack_Vθ )
    
    V          = [V_slack...; V_PV...; V_PQ...]

    x0_view    .= u_to_ΘV(V)

    gen_Q = imag.((Diagonal(V) * conj.(Y * V) + perm_matrix' * Sd_view)[PV_Idx_p])

    update_gens_Q!(Sg_view, gen_Q, p_update_gens_Q)     

    if !(ΔP_ΔQ == nothing)
        
        mismatch  = Diagonal(V) * conj.(Y * V) - perm_matrix' * ( Cg * Sg_view - Sd_view )

        mismatch  = Diagonal(V) * conj.(Y * V) - perm_matrix' * ( Cg * Sg_view - Sd_view )

        ΔP_ΔQ .= Array( [[real.(mismatch[PV_Idx]), real.(mismatch[PQ_Idx]), imag.(mismatch[PQ_Idx])]...; ] )        
    end
    
    if !(Jac_Matrix == nothing)
        #
        m11 = [[i==j ? ∂P_i_∂δ_i(i,V,Y) : ∂P_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

        M11 = reshape( [m11...;], size(Y) )

        n12 = [[i==j ? Vm_i_∂P_i_∂Vm_i(i,V,Y) : Vm_j_∂P_i_∂Vm_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

        N12 = reshape( [n12...;], size(Y) )

        n21 = [[i==j ? ∂Q_i_∂δ_i(i,V,Y) : ∂Q_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

        N21 = reshape( [n21...;], size(Y) )

        m22 = [[i==j ? Vm_i_∂Q_i_Vm_i(i,V,Y) : Vm_j_∂Q_i_Vm_j(i,j,V,Y) for j in 1:length(V)] for i in  1:length(V)] 

        M22 = reshape( [m22...;], size(Y) )  

        Jac_Matrix .= vcat(hcat(M11'[PV_PQ_Idx, PV_PQ_Idx], N12'[PV_PQ_Idx, PQ_Idx]),
                    hcat(N21'[PQ_Idx,    PV_PQ_Idx], M22'[PQ_Idx,    PQ_Idx]))        
    end

    
end



function x_dynamic_reduced_pf_sol( x0, ( gen_Q_view, Sg_view, Sd_view ), p; maxiter=40, ftol=1000*eps(), xtol=1000*eps() )
    
    p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack = p

    nodes_vmax, nodes_vmin, gen_Vm, p_Idx = p_V_limits

    gens_bus_num, gen_Qmax, gen_Qmin  = p_update_gens_Q

    perm_matrix, Cg, Sd, Y = p_ΔP_ΔQ 

    PV_Idx, PQ_Idx = p_Idx

    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    x_Θ_Idx, x_v_Idx, x_Θ_PV_PQ_Idx, x_v_PV_PQ_Idx, x_Θ_PV_Idx, x_v_PV_Idx, x_Θ_PQ_Idx, x_v_PQ_Idx, Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx, Δx_Θ_PV_Idx, Δx_v_PV_Idx, Δx_Θ_PQ_Idx, Δx_v_PQ_Idx = p_x_Θ_Δx_Idx

    slack_Vm, slack_Vθ, x_Θ_SK_Idx, x_v_SK_Idx = p_slack 

    ####

    x0[x_Θ_SK_Idx]   .= slack_Vθ

    x0[ x_v_SK_Idx]  .= slack_Vm
    
    x0[x_v_PV_Idx]    .= gen_Vm

    x0_view          = @view x0[:]
    
    x_θ_v_red_view   = @view x0[[collect(x_Θ_PV_PQ_Idx); collect(x_v_PQ_Idx)]]
    
    dims_view = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PQ_Idx)]

    _,_, Θ_PV_Θ_PQ_v_PQ_Idx = create_size_offset_Idx(dims_view)    

    # zeros is not used for a mismatch
    
    ΔP_ΔQ = ones(sum([length(PV_Idx),
                      length(PQ_Idx),
                      length(PQ_Idx)]))

    Jac_col_size = Jac_row_size = sum([length(PV_PQ_Idx), length(PQ_Idx)])
    
    Jac_Matrix   = zeros(Jac_row_size, Jac_col_size)


    x_0 = [x_θ_v_red_view;]

    ΔP_ΔQ_0 = similar(x_0)

    param = (x0_view, Sg_view, Sd_view, slack_Vm, slack_Vθ, Θ_PV_Θ_PQ_v_PQ_Idx, p ) 
    
    ####
    
    x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param )

    x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param )

    # ------------------------------------------------------

    # alg = SimpleNewtonRaphson()
    
    # nl_func = NonlinearFunction((ΔP_ΔQ, x_θ_v_red_view, param) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param), jac = (Jac_Matrix, x_θ_v_red_view, param) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param) )
    
    # nl_prob = NonlinearProblem(nl_func, x_0, param)
    # sol1 = NonlinearSolve.solve(nl_prob, alg)
    # sol2 = NonlinearSolve.solve(nl_prob, NewtonRaphson(; concrete_jac=true) )

    # ------------------------------------------------------

    """
    Other optional arguments to nlsolve, available for all algorithms, are:

    xtol: norm difference in x between two successive iterates under which convergence is declared. Default: 0.0.

    ftol: infinite norm of residuals under which convergence is declared. Default: 1e-8.

    iterations: maximum number of iterations. Default: 1_000.
    store_trace: should a trace of the optimization algorithm's state be stored? Default: false.

    show_trace: should a trace of the optimization algorithm's state be shown on STDOUT? Default: false.

    extended_trace: should additifonal algorithm internals be added to the state trace? Default: false.

    """

    sol = nlsolve( (ΔP_ΔQ, x_θ_v_red_view) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param), (Jac_Matrix, x_θ_v_red_view) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param ), x_0; ftol =  ftol )

    # ------------------------------------------------------

    V = x_to_V(x0_view, x_Θ_Idx, x_v_Idx)

    return  sol.iterations, V, sol.residual_norm, sol.x_converged
    

end


function x_get_dynamic_powerflow_sol( x0::Vector{Float64}, (Sg_view, Sd_view), dyn_powerflow_sol_param )

    Ybus_noperm, nodes_name, branches_name, num_nodes, perm_matrix, gen_Q_view, _, _, p = dyn_powerflow_sol_param
    
    _, _, p_ΔP_ΔQ, _, p_x_Θ_Δx_Idx, p_slack = p

    _, Cg, Sd, Y = p_ΔP_ΔQ
    
    # ------------------------------------------------------
    # Powerflow
    # ------------------------------------------------------
    
    _, V_k, _, _ = x_dynamic_reduced_pf_sol(x0, (gen_Q_view, Sg_view, Sd_view), p; maxiter=40, ftol=1000*eps(), xtol=1000*eps())
    
    V = perm_matrix * V_k

    Iinj = nodal_current = Ybus_noperm * V + Diagonal(inv.(conj.(V))) * ( conj.( Sd_view ) )

    Inet_inj = net_current = Ybus_noperm * V 

    Sbus_n = Diagonal(V) * conj.(Ybus_noperm * V)

    GenSinj = Sbus_n + Sd_view

    # return [ (vh, θh, ph, qh, ih_r, ih_i) for (vh, θh, ph, qh, ih_r, ih_i) in zip(abs.(V), angle.(V), real.(Sbus_n), imag.(Sbus_n), real.(Inet_inj), imag.(Inet_inj)) ]

   return [ (vh, θh, ph, qh, ih_r, ih_i, pg, qg, ig_r, ig_i) for (vh, θh, ph, qh, ih_r, ih_i, pg, qg, ig_r, ig_i) in zip(abs.(V), angle.(V), real.(Sbus_n), imag.(Sbus_n), real.(Inet_inj), imag.(Inet_inj), real.(GenSinj), imag.(GenSinj), real.(Iinj), imag.(Iinj) ) ]    
end



function x_get_dynamic_reduced_powerflow( case_data_func )

    # Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

    Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

    perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

    num_nodes     = get_network_nodal_size(Branches)

    gens_S        = get_nodes_S(Gen_nodes)

    gen_Q         = imag(gens_S) 

    gens_bus_num  = get_gens_bus_num(Gen_nodes)

    Cg            = get_Cg(Gen_nodes, num_nodes)

    Sg            = get_Sg(Gen_nodes, num_nodes)

    gen_Vm        = get_nodes_Vm(Gen_nodes) 

    num_gen       =  length(gen_Vm)

    gen_Qmax      = get_nodes_Qmax(Gen_nodes)

    gen_Qmin      = get_nodes_Qmin(Gen_nodes)

    Sd            = get_loads_Sd(Load_nodes)


    # branches_name  = get_branches_name(Branches)   
    # nodes_name     = get_nodes_Bus_name(Nodes)

    branches_name  = collect(keys(Branches))
    
    nodes_name     = collect(keys(Nodes))      
    
    Yf             = get_Yf(Nodes, Branches)
    
    Yt             = get_Yt(Nodes, Branches)
    

    Ybus_noperm   = get_Ybus(Nodes, Branches, shunts=Shunts)

    slack_Vm      = get_nodes_Vm(Slack_nodes)

    slack_Vθ      = get_nodes_Vθ(Slack_nodes)

    nodes_vmax    = get_nodes_vmax(Nodes)

    nodes_vmax_p  = perm_matrix' * nodes_vmax

    nodes_vmin    = get_nodes_vmin(Nodes)

    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    = get_nodes_Bus_name(Nodes)

    Ybus_p        = perm_matrix' * (Ybus_noperm) * perm_matrix

    # ------------------------------------------------------
    # Indices
    # ------------------------------------------------------


    SK_Idx_p      = idx_slack_nodes

    PV_Idx_p      = idx_source_nodes   # gen nodes indices from perm 

    PQ_Idx_p      = idx_demand_nodes   # load nodes indices from perm

    PV_PQ_Idx_p   = first(PV_Idx_p):last(PQ_Idx_p)  

    

    SK_PV_PQ_V_dim     = [length(SK_Idx_p),length(PV_Idx_p),length(PQ_Idx_p)]
    SK_PV_PQ_V_offset  = create_offsets(SK_PV_PQ_V_dim)
    SK_PV_PQ_V_Idx     = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

    SK_V_Idx           = SK_PV_PQ_V_Idx[1]
    PV_V_Idx           = SK_PV_PQ_V_Idx[2]
    PQ_V_Idx           = SK_PV_PQ_V_Idx[3]

    # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

    x_dim      = [length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx),
                  length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx)]

    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    x_Θ_SK_Idx = x_Idx[1]
    x_Θ_PV_Idx = x_Idx[2]
    x_Θ_PQ_Idx = x_Idx[3]
    
    x_v_SK_Idx = x_Idx[4]
    x_v_PV_Idx = x_Idx[5]
    x_v_PQ_Idx = x_Idx[6]

    x_Θ_Idx       = first(x_Θ_SK_Idx):last(x_Θ_PQ_Idx)
    x_v_Idx       = first(x_v_SK_Idx):last(x_v_PQ_Idx)

    x_Θ_PV_PQ_Idx = first(x_Θ_PV_Idx):last(x_Θ_PQ_Idx)
    x_v_PV_PQ_Idx = first(x_v_PV_Idx):last(x_v_PQ_Idx)

    Δx_dim = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PV_Idx),length(x_v_PQ_Idx)]

    Δx_offset   = create_offsets( Δx_dim )
    Δx_Idx      = create_idxs( Δx_offset, Δx_dim )

    Δx_Θ_PV_Idx = Δx_Idx[1]
    Δx_Θ_PQ_Idx = Δx_Idx[2]

    Δx_v_PV_Idx = Δx_Idx[3]
    Δx_v_PQ_Idx = Δx_Idx[4]

    Δx_Θ_PV_PQ_Idx = first(Δx_Θ_PV_Idx):last(Δx_Θ_PQ_Idx)

    Δx_v_PV_PQ_Idx = first(Δx_v_PV_Idx):last(Δx_v_PQ_Idx)

    # ------------------------------------------------------

    # computing arrays, indices and views
    # ------------------------------------------------------

    gen_Q      = zeros(Float64, num_gen)

    gen_Q_view =  view(gen_Q, 1:num_gen)

    Sg_view    = @view Sg[:]

    update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

    Sd_view    = @view Sd[:]
    
    Y          = Ybus_p

    # ------------------------------------------------------
    # Parameteres
    # ------------------------------------------------------
    
    p_ΔP_ΔQ            = ( perm_matrix, Cg, Sd, Y )

    p_update_gens_Q    = ( gens_bus_num, gen_Qmax, gen_Qmin )

    p_x_Θ_Δx_Idx       = ( x_Θ_Idx,        x_v_Idx,
                           x_Θ_PV_PQ_Idx,  x_v_PV_PQ_Idx,
                           x_Θ_PV_Idx,     x_v_PV_Idx,
                           x_Θ_PQ_Idx,     x_v_PQ_Idx,
                           
                           Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx,
                           Δx_Θ_PV_Idx,    Δx_v_PV_Idx,
                           Δx_Θ_PQ_Idx,    Δx_v_PQ_Idx )

    p_Idx = (PV_Idx_p, PQ_Idx_p)
    
    p_V_limits = ( nodes_vmax, nodes_vmin, gen_Vm, p_Idx )

    p_slack = (slack_Vm, slack_Vθ,x_Θ_SK_Idx, x_v_SK_Idx )

    p = ( p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack )

     dyn_powerflow_sol_param = (; Ybus_noperm, nodes_name, branches_name, num_nodes, perm_matrix, gen_Q_view, Sg_view, Sd_view, p)

    PV_Idx, PQ_Idx = p_Idx

    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)
    
    # ------------------------------------------------------
    # Initial condition
    # ------------------------------------------------------
    
    x_Θ                = zeros(num_nodes)

    x_v                = ones(num_nodes)

    x0         = [x_Θ...; x_v...]

    # x0_view   = @view x0[:]

    _, V, _, _ = x_dynamic_reduced_pf_sol(x0, (gen_Q_view, Sg_view, Sd_view), p; maxiter=40, ftol=1000*eps(), xtol=1000*eps())
    
 
    # V = x_to_V(x0, x_Θ_Idx, x_v_Idx)

    # V .= perm_matrix * V

    V = perm_matrix * V
    
    # ------------------------------------------------------
    # Results
    # ------------------------------------------------------
    
    Iinj = nodal_current  = Ybus_noperm * V + Diagonal(inv.(conj.(V))) * ( conj.( Sd_view ) )

    Inet_inj = net_current = Ybus_noperm * V 

    Sbus_n  = Diagonal(V) * conj.(Ybus_noperm * V)

    GenSinj = Sbus_n + Sd_view

    If = Yf * V
    
    It = Yt * V

    Ibranches = If + It

    bus_dict_Iinj = OrderedDict( name => [ih_r, ih_i] for (name, ih_r, ih_i) in zip(nodes_name, real.(Iinj), imag.(Iinj)))

    bus_dict_Inet_inj = OrderedDict( name => [ih_r, ih_i] for (name, ih_r, ih_i) in zip(nodes_name,real.(Inet_inj), imag.(Inet_inj))) 

    branch_dict_init = OrderedDict(name => ( real(i_f), imag(i_f), real(i_t), imag(i_t) ) for (name, i_f, i_t) in zip(branches_name, If, It))

   bus_dict_init = OrderedDict( name => (vh, θh, ph, qh, ih_r, ih_i, pg, qg, ig_r, ig_i ) for (name, vh, θh, ph, qh, ih_r, ih_i, pg, qg, ig_r, ig_i ) in zip(nodes_name,  abs.(V),  angle.(V), real.(Sbus_n), imag.(Sbus_n), real.(Inet_inj), imag.(Inet_inj), real.(GenSinj), imag.(GenSinj), real.(Iinj), imag.(Iinj)) )
        
    return perm_matrix, ((Sg_view, Sd_view), dyn_powerflow_sol_param), Dict("Iinj" => Iinj, "Inet_inj" => Inet_inj, "Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V), "Vθ" => angle.(V), "Vbus" => V, "Ibranches" => Ibranches, "Ybus" => Ybus_noperm , "Sbus" => Sbus_n, "GenSinj" => GenSinj, "dict_init" =>Dict("bus_dict_init" => bus_dict_init,  "branch_dict_init" =>  branch_dict_init, "bus_dict_Iinj" => bus_dict_Iinj, "bus_dict_Inet_inj" => bus_dict_Inet_inj ), "Ifrom" => If, "Ito" => It   )
        
end


#########################################


function test_x_dynamic_reduced_pf_sol(x0, (gen_Q_view, Sg_view, Sd_view), p; maxiter=40, ftol=1000*eps(), xtol=1000*eps())
    
    p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack = p

    nodes_vmax, nodes_vmin, gen_Vm, p_Idx = p_V_limits

    gens_bus_num, gen_Qmax, gen_Qmin  = p_update_gens_Q

    perm_matrix, Cg, Sd, Y = p_ΔP_ΔQ 

    PV_Idx, PQ_Idx = p_Idx

    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

    x_Θ_Idx, x_v_Idx, x_Θ_PV_PQ_Idx, x_v_PV_PQ_Idx, x_Θ_PV_Idx, x_v_PV_Idx, x_Θ_PQ_Idx, x_v_PQ_Idx, Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx, Δx_Θ_PV_Idx, Δx_v_PV_Idx, Δx_Θ_PQ_Idx, Δx_v_PQ_Idx = p_x_Θ_Δx_Idx

    slack_Vm, slack_Vθ, x_Θ_SK_Idx, x_v_SK_Idx = p_slack 
    
    ####

    # p_slack = (slack_Vm, slack_Vθ,x_Θ_SK_Idx, x_v_SK_Idx)

    x0[x_Θ_SK_Idx]   .= slack_Vθ

    x0[ x_v_SK_Idx]  .= slack_Vm
    
    x0[x_v_PV_Idx]    .= gen_Vm

    x0_view          = @view x0[:]
    
    x_θ_v_red_view   = @view x0[[collect(x_Θ_PV_PQ_Idx);
                                 collect(x_v_PQ_Idx)]]
    dims_view = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PQ_Idx)]

    _,_, Θ_PV_Θ_PQ_v_PQ_Idx = create_size_offset_Idx(dims_view)    

    # zeros is not used for a mismatch
    
    ΔP_ΔQ = ones(sum([length(PV_Idx),
                      length(PQ_Idx),
                      length(PQ_Idx)]))

    Jac_col_size = sum([length(PV_PQ_Idx), length(PQ_Idx)])
    
    Jac_row_size = Jac_col_size

    Jac_Matrix   = zeros(Jac_row_size, Jac_col_size)

    x_0 = [x_θ_v_red_view;]

    ΔP_ΔQ_0 = similar(x_0)

    param = (x0_view, Sg_view, Sd_view, slack_Vm, slack_Vθ, Θ_PV_Θ_PQ_v_PQ_Idx, p ) 
    
    ####
    
    x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param )

    x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param )

    # ------------------------------------------------------

    # alg = SimpleNewtonRaphson()
    
    # nl_func = NonlinearFunction((ΔP_ΔQ, x_θ_v_red_view, param) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param), jac = (Jac_Matrix, x_θ_v_red_view, param) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param) )
    
    # nl_prob = NonlinearProblem(nl_func, x_0, param)
    # sol1 = NonlinearSolve.solve(nl_prob, alg)
    # sol2 = NonlinearSolve.solve(nl_prob, NewtonRaphson(; concrete_jac=true) )

    # ------------------------------------------------------
    
    # x_θ_v_red_view .= sol.u

    # V = x_to_V(x0, x_Θ_Idx, x_v_Idx)

    # V .= perm_matrix * V
    
    # # ------------------------------------------------------
    """
    Other optional arguments to nlsolve, available for all algorithms, are:

    xtol: norm difference in x between two successive iterates under which convergence is declared. Default: 0.0.

    ftol: infinite norm of residuals under which convergence is declared. Default: 1e-8.

    iterations: maximum number of iterations. Default: 1_000.
    store_trace: should a trace of the optimization algorithm's state be stored? Default: false.

    show_trace: should a trace of the optimization algorithm's state be shown on STDOUT? Default: false.

    extended_trace: should additifonal algorithm internals be added to the state trace? Default: false.

    """

    sol = nlsolve( (ΔP_ΔQ, x_θ_v_red_view) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param), (Jac_Matrix, x_θ_v_red_view) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param ), x_0; ftol =  ftol )

    # ftol = 1e-12

    # V = x_to_V(x0_view, x_Θ_Idx, x_v_Idx)

    # V_norm = perm_matrix * V
    
    # # ------------------------------------------------------
    
    # df = OnceDifferentiable((ΔP_ΔQ, x_θ_v_red_view) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param ), (Jac_Matrix, x_θ_v_red_view) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param ), x_0, ΔP_ΔQ_0)

    # sol = nlsolve(df, x_0 )

    # V = x_to_V(x0_view, x_Θ_Idx, x_v_Idx)

    # V_norm = perm_matrix * V    

    # # ------------------------------------------------------

    # sol = nlsolve( only_fj!((ΔP_ΔQ, Jac_Matrix, x_θ_v_red_view) -> x_fj!(ΔP_ΔQ, Jac_Matrix, x_θ_v_red_view, param)), x_0)

    # # ------------------------------------------------------

    V = x_to_V(x0_view, x_Θ_Idx, x_v_Idx)

    return  sol.iterations, V, sol.residual_norm, sol.x_converged
    

end



function test_x_get_dynamic_reduced_powerflow( case_data_func )

    # Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

    Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

    perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

    num_nodes     = get_network_nodal_size(Branches)

    gens_S        = get_nodes_S(Gen_nodes)

    gen_Q         = imag(gens_S) 

    gens_bus_num  = get_gens_bus_num(Gen_nodes)

    Cg            = get_Cg(Gen_nodes, num_nodes)

    Sg            = get_Sg(Gen_nodes, num_nodes)

    gen_Vm        = get_nodes_Vm(Gen_nodes) 

    num_gen       =  length(gen_Vm)

    gen_Qmax      = get_nodes_Qmax(Gen_nodes)

    gen_Qmin      = get_nodes_Qmin(Gen_nodes)

    Sd            = get_loads_Sd(Load_nodes)


    # branches_name  = get_branches_name(Branches)   
    # nodes_name     = get_nodes_Bus_name(Nodes)

    branches_name  = collect(keys(Branches))
    
    nodes_name     = collect(keys(Nodes))      
    
    Yf             = get_Yf(Nodes, Branches)
    
    Yt             = get_Yt(Nodes, Branches)
    

    Ybus_noperm   = get_Ybus(Nodes, Branches, shunts=Shunts)

    slack_Vm      = get_nodes_Vm(Slack_nodes)

    slack_Vθ      = get_nodes_Vθ(Slack_nodes)

    nodes_vmax    = get_nodes_vmax(Nodes)

    nodes_vmax_p  = perm_matrix' * nodes_vmax

    nodes_vmin    = get_nodes_vmin(Nodes)

    nodes_vmin_p  = perm_matrix' * nodes_vmin

    nodes_name    = get_nodes_Bus_name(Nodes)

    Ybus_p        = perm_matrix' * (Ybus_noperm) * perm_matrix

    # ------------------------------------------------------
    # Indices
    # ------------------------------------------------------


    SK_Idx_p      = idx_slack_nodes

    PV_Idx_p      = idx_source_nodes   # gen nodes indices from perm 

    PQ_Idx_p      = idx_demand_nodes   # load nodes indices from perm

    PV_PQ_Idx_p   = first(PV_Idx_p):last(PQ_Idx_p)  

    

    SK_PV_PQ_V_dim     = [length(SK_Idx_p),length(PV_Idx_p),length(PQ_Idx_p)]
    SK_PV_PQ_V_offset  = create_offsets(SK_PV_PQ_V_dim)
    SK_PV_PQ_V_Idx     = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

    SK_V_Idx           = SK_PV_PQ_V_Idx[1]
    PV_V_Idx           = SK_PV_PQ_V_Idx[2]
    PQ_V_Idx           = SK_PV_PQ_V_Idx[3]

    # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

    x_dim      = [length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx),
                  length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx)]

    x_offset   = create_offsets(x_dim)
    x_Idx      = create_idxs(x_offset, x_dim)

    x_Θ_SK_Idx = x_Idx[1]
    x_Θ_PV_Idx = x_Idx[2]
    x_Θ_PQ_Idx = x_Idx[3]
    
    x_v_SK_Idx = x_Idx[4]
    x_v_PV_Idx = x_Idx[5]
    x_v_PQ_Idx = x_Idx[6]

    x_Θ_Idx       = first(x_Θ_SK_Idx):last(x_Θ_PQ_Idx)
    x_v_Idx       = first(x_v_SK_Idx):last(x_v_PQ_Idx)

    x_Θ_PV_PQ_Idx = first(x_Θ_PV_Idx):last(x_Θ_PQ_Idx)
    x_v_PV_PQ_Idx = first(x_v_PV_Idx):last(x_v_PQ_Idx)

    Δx_dim = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PV_Idx),length(x_v_PQ_Idx)]

    Δx_offset   = create_offsets( Δx_dim )
    Δx_Idx      = create_idxs( Δx_offset, Δx_dim )

    Δx_Θ_PV_Idx = Δx_Idx[1]
    Δx_Θ_PQ_Idx = Δx_Idx[2]

    Δx_v_PV_Idx = Δx_Idx[3]
    Δx_v_PQ_Idx = Δx_Idx[4]

    Δx_Θ_PV_PQ_Idx = first(Δx_Θ_PV_Idx):last(Δx_Θ_PQ_Idx)

    Δx_v_PV_PQ_Idx = first(Δx_v_PV_Idx):last(Δx_v_PQ_Idx)

    # ------------------------------------------------------

    # computing arrays, indices and views
    # ------------------------------------------------------

    gen_Q      = zeros(Float64, num_gen)

    gen_Q_view =  view(gen_Q, 1:num_gen)

    Sg_view    = @view Sg[:]

    update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

    Sd_view    = @view Sd[:]
    
    Y          = Ybus_p

    # ------------------------------------------------------
    # Parameteres
    # ------------------------------------------------------
    
    p_ΔP_ΔQ            = ( perm_matrix, Cg, Sd, Y )

    p_update_gens_Q    = ( gens_bus_num, gen_Qmax, gen_Qmin )

    p_x_Θ_Δx_Idx       = ( x_Θ_Idx,        x_v_Idx,
                           x_Θ_PV_PQ_Idx,  x_v_PV_PQ_Idx,
                           x_Θ_PV_Idx,     x_v_PV_Idx,
                           x_Θ_PQ_Idx,     x_v_PQ_Idx,
                           
                           Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx,
                           Δx_Θ_PV_Idx,    Δx_v_PV_Idx,
                           Δx_Θ_PQ_Idx,    Δx_v_PQ_Idx )

    p_Idx = (PV_Idx_p, PQ_Idx_p)
    
    p_V_limits = ( nodes_vmax, nodes_vmin, gen_Vm, p_Idx )

    p_slack = (slack_Vm, slack_Vθ,x_Θ_SK_Idx, x_v_SK_Idx )

    p = ( p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack )

     dyn_powerflow_sol_param = (; Ybus_noperm, nodes_name, branches_name, num_nodes, perm_matrix, gen_Q_view, Sg_view, Sd_view, p)

    PV_Idx, PQ_Idx = p_Idx

    PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)
    
    # ------------------------------------------------------
    # Initial condition
    # ------------------------------------------------------
    
    x_Θ                = zeros(num_nodes)

    x_v                = ones(num_nodes)

    x0                 = [x_Θ...; x_v...]

    _, V_norm, _, _ = x_dynamic_reduced_pf_sol( x0, ( gen_Q_view, Sg_view, Sd_view ), p; maxiter=40, ftol=1000*eps(), xtol=1000*eps() )
    
    V = perm_matrix * V_norm
    

    # x0[x_Θ_SK_Idx]   .= slack_Vθ

    # x0[ x_v_SK_Idx]  .= slack_Vm
    
    # x0[x_v_PV_Idx]    .= gen_Vm

    # x0_view          = @view x0[:]
    
    # x_θ_v_red_view   = @view x0[[collect(x_Θ_PV_PQ_Idx);
    #                              collect(x_v_PQ_Idx)]]
    # dims_view = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PQ_Idx)]

    # _,_, Θ_PV_Θ_PQ_v_PQ_Idx = create_size_offset_Idx(dims_view)    

    # # zeros is not used for a mismatch
    
    # ΔP_ΔQ = ones(sum([length(PV_Idx),
    #                   length(PQ_Idx),
    #                   length(PQ_Idx)]))

    # Jac_col_size = sum([length(PV_PQ_Idx), length(PQ_Idx)])
    
    # Jac_row_size = Jac_col_size

    # Jac_Matrix   = zeros(Jac_row_size, Jac_col_size)


    # x_0 = [x_θ_v_red_view;]

    # ΔP_ΔQ_0 = similar(x_0)

    # # param = (x0, Sg_view, Sd_view, p )

    # param = (x0_view, Sg_view, Sd_view, slack_Vm, slack_Vθ, Θ_PV_Θ_PQ_v_PQ_Idx, p ) 
    
    # # ------------------------------------------------------
    # # Powerflow
    # # ------------------------------------------------------ 

    
    # x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param )

    # x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param )

    # # ------------------------------------------------------

    # # alg = SimpleNewtonRaphson()
    
    # # nl_func = NonlinearFunction((ΔP_ΔQ, x_θ_v_red_view, param) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param), jac = (Jac_Matrix, x_θ_v_red_view, param) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param) )
    
    # # nl_prob = NonlinearProblem(nl_func, x_0, param)

    # # sol1 = NonlinearSolve.solve(nl_prob, alg)

    # # sol2 = NonlinearSolve.solve(nl_prob, NewtonRaphson(; concrete_jac=true) )

    # # # ------------------------------------------------------
    
    # # x_θ_v_red_view .= sol.u

    # # V = x_to_V(x0, x_Θ_Idx, x_v_Idx)

    # # V .= perm_matrix * V
    
    # # # ------------------------------------------------------
    # """
    # Other optional arguments to nlsolve, available for all algorithms, are:

    # xtol: norm difference in x between two successive iterates under which convergence is declared. Default: 0.0.

    # ftol: infinite norm of residuals under which convergence is declared. Default: 1e-8.

    # iterations: maximum number of iterations. Default: 1_000.
    # store_trace: should a trace of the optimization algorithm's state be stored? Default: false.

    # show_trace: should a trace of the optimization algorithm's state be shown on STDOUT? Default: false.

    # extended_trace: should additifonal algorithm internals be added to the state trace? Default: false.

    # """

    # sol = nlsolve( (ΔP_ΔQ, x_θ_v_red_view) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param), (Jac_Matrix, x_θ_v_red_view) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param ), x_0; ftol =  1e-12 )

    # # V = x_to_V(x0_view, x_Θ_Idx, x_v_Idx)

    
    # # ------------------------------------------------------
    
    # df = OnceDifferentiable((ΔP_ΔQ, x_θ_v_red_view) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param ), (Jac_Matrix, x_θ_v_red_view) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param ), x_0, ΔP_ΔQ_0)

    # sol = nlsolve(df, x_0 )

    # V = x_to_V(x0_view, x_Θ_Idx, x_v_Idx)

    # V_norm = perm_matrix * V    

    # # ------------------------------------------------------

    # sol = nlsolve( only_fj!( (ΔP_ΔQ, Jac_Matrix, x_θ_v_red_view) -> x_fj!(ΔP_ΔQ, Jac_Matrix, x_θ_v_red_view, param )), x_0 )

    # # ------------------------------------------------------

    # # x_θ_v_red_view .= sol.zero
    # # V = x_to_V(x0, x_Θ_Idx, x_v_Idx)
    # # V .= perm_matrix * V

    # V = x_to_V(x0_view, x_Θ_Idx, x_v_Idx)

    # V_norm = perm_matrix * V
    
    # ------------------------------------------------------
    # Results
    # ------------------------------------------------------
    
    Iinj = nodal_current  = Ybus_noperm * V + Diagonal(inv.(conj.(V))) * ( conj.( Sd_view ) )

    Inet_inj  = net_current  = Ybus_noperm * V 

    Sbus_n     = Diagonal(V) * conj.(Ybus_noperm * V)

    GenSinj  = Sbus_n + Sd_view

    return [(vh, θh, ph, qh, ih_r, ih_i) for (vh, θh, ph, qh, ih_r, ih_i) in zip(abs.(V), angle.(V), real.(Sbus_n), imag.(Sbus_n), real.(Inet_inj), imag.(Inet_inj) )]
    
end


########################################


# function get_reduced_powerflow_sol_param( case_data_func::Function )

#     Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()


#     perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

#     num_nodes = get_network_nodal_size(Branches)

#     gens_S = get_nodes_S(Gen_nodes)

#     gen_Q = imag(gens_S) 

#     gens_bus_num = get_gens_bus_num(Gen_nodes)

#     Cg = get_slack_and_gens_Cg(Nodes)

#     # Cg = get_Cg(Gen_nodes, num_nodes)

#     Sg = get_slack_and_gen_Sg(Nodes)

#     # Sg = get_Sg(Gen_nodes, num_nodes)

#     gen_Vm  = get_nodes_Vm(Gen_nodes) 

#     num_gen =  length(gen_Vm)

#     gen_Qmax = get_nodes_Qmax(Gen_nodes)

#     gen_Qmin = get_nodes_Qmin(Gen_nodes)

#     Sd = get_loads_Sd(Load_nodes)

#     # branches_name = get_branches_name(Branches)
#     # nodes_name    = get_nodes_Bus_name(Nodes)

#     branches_name = collect(keys(Branches))
    
#     nodes_name = collect(keys(Nodes))    
    
#     Yf = get_Yf(Nodes, Branches)
    
#     Yt = get_Yt(Nodes, Branches)

#     Ybus_noperm = get_Ybus(Nodes, Branches, shunts=Shunts)

#     slack_Vm = get_nodes_Vm(Slack_nodes)

#     slack_Vθ = get_nodes_Vθ(Slack_nodes)

#     nodes_vmax = get_nodes_vmax(Nodes)

#     nodes_vmax_p = perm_matrix' * nodes_vmax

#     nodes_vmin = get_nodes_vmin(Nodes)

#     nodes_vmin_p = perm_matrix' * nodes_vmin

#     Ybus_p = perm_matrix' * (Ybus_noperm) * perm_matrix
    
#     return powerflow_data = (; perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes, num_nodes, gens_S, gen_Q, gens_bus_num, Cg, Sg, gen_Vm, num_gen, gen_Qmax, gen_Qmin, Sd, branches_name, nodes_name, Yf, Yt, Ybus_noperm, slack_Vm, slack_Vθ, nodes_vmax, nodes_vmax_p, nodes_vmin, nodes_vmin_p,  Ybus_p )

# end

# # ------------------------------------------------------
# # Power balance
# # ------------------------------------------------------

# function get_reduced_powerflow_sol_param(powerflow_data::NamedTuple)

#     perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes, num_nodes, gens_S, gen_Q, gens_bus_num, Cg, Sg, gen_Vm, num_gen, gen_Qmax, gen_Qmin, Sd, branches_name, nodes_name, Yf, Yt, Ybus_noperm, slack_Vm, slack_Vθ, nodes_vmax, nodes_vmax_p, nodes_vmin, nodes_vmin_p, nodes_name, Ybus_p = powerflow_dat

#     # ------------------------------------------------------
#     # Indices
#     # ------------------------------------------------------

#     SK_Idx = idx_slack_nodes

#     PV_Idx = idx_source_nodes  

#     PQ_Idx = idx_demand_nodes 

#     SK_PV_PQ_V_dim    = [length(SK_Idx),length(PV_Idx),length(PQ_Idx)]
#     SK_PV_PQ_V_offset = create_offsets(SK_PV_PQ_V_dim)
#     SK_PV_PQ_V_Idx    = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

#     SK_V_Idx = SK_PV_PQ_V_Idx[1]
#     PV_V_Idx = SK_PV_PQ_V_Idx[2]
#     PQ_V_Idx = SK_PV_PQ_V_Idx[3]

#     # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

#     # x
    
#     x_dim      = [length(PV_V_Idx), length(PQ_V_Idx), length(PQ_V_Idx)]

#     x_offset   = create_offsets(x_dim)
#     x_Idx      = create_idxs(x_offset, x_dim)

#     x_Θ_PV_Idx = x_Idx[1]
#     x_Θ_PQ_Idx = x_Idx[2]    
#     x_v_PQ_Idx = x_Idx[3]

#     # Δx
    
#     Δx_dim = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_Θ_PV_Idx), length(x_v_PQ_Idx)]

#     Δx_offset   = create_offsets(Δx_dim )
#     Δx_Idx      = create_idxs(Δx_offset, Δx_dim)

#     Δx_Θ_PV_Idx = Δx_Idx[1]
#     Δx_Θ_PQ_Idx = Δx_Idx[2]

#     Δx_v_PV_Idx = Δx_Idx[3]
#     Δx_v_PQ_Idx = Δx_Idx[4]

#     # Jac
    
#     Jac_dim    = [length(PV_Idx), length(PQ_Idx), length(PQ_Idx)]
#     Jac_offset = create_offsets( Jac_dim )
#     Jac_Idx    = create_idxs(Jac_offset, Jac_dim)

#     Jac_P_PV_row_Idx  = Jac_Θ_PV_row_Idx  = Jac_Idx[1]
#     Jac_P_PQ_row_Idx  = Jac_Θ_PQ_row_Idx  = Jac_Idx[2]
#     Jac_Q_PQ_row_Idx  = Jac_v_PQ_row_Idx  = Jac_Idx[3]

    
#     # ------------------------------------------------------
#     # computing arrays, indices and views
#     # ------------------------------------------------------

#     # gen_Q      = zeros(Float64, num_gen)

#     gen_Q_view =  view(gen_Q, 1:num_gen)

#     Sg_view    = @view Sg[:]

#     update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

#     Sd_view    = @view Sd[:]
    
#     Y          = Ybus_p

#     # ------------------------------------------------------
#     # Parameteres
#     # ------------------------------------------------------

#         p_jac              = (Jac_P_PV_row_Idx, Jac_P_PQ_row_Idx, Jac_Q_PQ_row_Idx, Jac_Θ_PV_row_Idx, Jac_Θ_PQ_row_Idx, Jac_v_PQ_row_Idx)

#     p_jac              = (Jac_P_PV_row_Idx, Jac_P_PQ_row_Idx, Jac_Q_PQ_row_Idx)
    
#     p_ΔP_ΔQ            = (perm_matrix, Cg, Sd, Y )

#     p_update_gens_Q    = (gens_bus_num, gen_Qmax, gen_Qmin )

#     p_Idx              = (PV_Idx, PQ_Idx)

#     p_V_limits         = (nodes_vmax, nodes_vmin, gen_Vm, p_Idx)
    
#     p_x_Θ_Δx_Idx       = ( x_Θ_PV_Idx,
#                            x_Θ_PQ_Idx,
#                            x_v_PQ_Idx,
                           
#                            Δx_Θ_PV_Idx,
#                            Δx_Θ_PQ_Idx,
#                            Δx_v_PV_Idx,
#                            Δx_v_PQ_Idx )

#     p_slack = ( slack_Vm, slack_Vθ, x_Θ_SK_Idx, x_v_SK_Idx )

#     p                  = (p_V_limits,
#                           p_update_gens_Q,
#                           p_ΔP_ΔQ,
#                           p_Idx,
#                           p_x_Θ_Δx_Idx, p_jac, p_slack )

#     return powerflow_sol_param = (;num_nodes, perm_matrix, gen_Q_view, Sg_view, Sd_view, p)

# end


# function get_powerflow_sol( powerflow_sol_param )

#     num_nodes, perm_matrix, gen_Q_view, Sg_view, Sd_view, p = powerflow_sol_param
    
#     # ------------------------------------------------------
#     # Initial condition
#     # ------------------------------------------------------

#     x_Θ   = zeros(num_nodes)

#     x_v   = ones(num_nodes)

#     x0    = [x_Θ...; x_v...]

#     # ------------------------------------------------------
#     # Powerflow
#     # ------------------------------------------------------
    
#     k, V_k, ΔP_ΔQ_k, Δx = reduced_pf_sol(x0, (gen_Q_view, Sg_view, Sd_view), p; maxiter=40, ftol=1000*eps(), xtol=1000*eps())

#     V = perm_matrix * V_k

#     return V


# end



# function get_dynamic_powerflow_results_dict(V, Ybus_noperm,  Sd_view, nodes_name, branches_name )
    
#     # ------------------------------------------------------
#     # Results
#     # ------------------------------------------------------

#     Iinj = nodal_current  = Ybus_noperm * V + Diagonal(inv.(conj.(V))) * ( conj.( Sd_view ) )

#     Inet_inj  = net_current  = Ybus_noperm * V 

#     Sbus     = Diagonal(V) * conj.(Ybus_noperm * V)

#     GenSinj  = Sbus + Sd_view
    

#     branch_dict_init = OrderedDict(name =>
#         (real(i_f),imag(i_f), real(i_t), imag(i_t) ) for (name, i_f, i_t) in zip(branches_name, If, It))


#     bus_dict_init = OrderedDict(name =>
#         (vh, θh, ph, qh, ih_r, ih_i ) for (name, vh, θh, ph, qh, ih_r, ih_i) in zip(nodes_name,
#                                          round.(abs.(V);        digits=4),
#                                          round.(angle.(V);      digits=4),
#                                          round.(real.(Sbus);    digits=4),
#                                          round.(imag.(Sbus);    digits=4),
#                                          round.(real.(Inet_inj);digits=4),
#                                          round.(imag.(Inet_inj);digits=4) ))            
    
#     # return perm_matrix, k, V_k, ΔP_ΔQ_k, Δx
    
#     return Dict("Iinj"=> Iinj, "Inet_inj"=> Inet_inj, "Sbus" => Sbus, "GenSinj" => GenSinj, "dict_init" =>Dict("bus_dict_init" => bus_dict_init ))     
    
# end



# function get_powerflow_results_dict(V,  Ybus_noperm, Yf, Yt, Sd, nodes_name, branches_name )
    
#     # ------------------------------------------------------
#     # Results
#     # ------------------------------------------------------

#     Sbus      = Diagonal(V) * conj.(Ybus_noperm * V)

#     GenSinj   = Diagonal(V) * conj.(Ybus_noperm * V) + Sd
        
#     If         = Yf * V
    
#     It         = Yt * V

#     Ibranches  = If + It

#     branch_dict_init = OrderedDict(name =>
#         (real(i_f),imag(i_f), real(i_t), imag(i_t) ) for (name, i_f, i_t) in zip(branches_name, If, It))

#     bus_dict_init = OrderedDict(name =>
#         (vh, θh, ph, qh) for (name,
#                               vh,
#                               θh,
#                               ph,
#                               qh) in zip(nodes_name,
#                                          round.(abs.(V);     digits=4),
#                                          round.(angle.(V);   digits=4),
#                                          round.(real.(Sbus); digits=4),
#                                          round.( imag.(Sbus);digits=4) ))      
    
#     # return perm_matrix, k, V_k, ΔP_ΔQ_k, Δx
    
#     return Dict("Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V), "Vθ" => angle.(V), "Vbus" => V, "Ibranches" => Ibranches, "Ybus" => Ybus_noperm , "Sbus" => Sbus, "GenSinj" => GenSinj, "dict_init" =>Dict("bus_dict_init" => bus_dict_init,  "branch_dict_init" =>  branch_dict_init))     
    
# end



# function get_dynamic_reduced_powerflow_sol_param( case_data_func::Function )

#     Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

#     perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

#     num_nodes     = get_network_nodal_size(Branches)

#     gens_S        = get_nodes_S(Gen_nodes)

#     gen_Q         = imag(gens_S) 

#     gens_bus_num  = get_gens_bus_num(Gen_nodes)

#     Cg            = get_Cg(Gen_nodes, num_nodes)

#     Sg            = get_Sg(Gen_nodes, num_nodes)

#     gen_Vm        = get_nodes_Vm(Gen_nodes) 

#     num_gen       = length(gen_Vm)

#     gen_Qmax      = get_nodes_Qmax(Gen_nodes)

#     gen_Qmin      = get_nodes_Qmin(Gen_nodes)

#     Sd            = get_loads_Sd(Load_nodes)


#     # branches_name  = get_branches_name(Branches)   
#     # nodes_name     = get_nodes_Bus_name(Nodes)

#     branches_name  = collect(keys(Branches))
    
#     nodes_name     = collect(keys(Nodes))      
    
#     Yf             = get_Yf(Nodes, Branches)
    
#     Yt             = get_Yt(Nodes, Branches)
    

#     Ybus_noperm   = get_Ybus(Nodes, Branches, shunts=Shunts)

#     slack_Vm      = get_nodes_Vm(Slack_nodes)

#     slack_Vθ      = get_nodes_Vθ(Slack_nodes)

#     nodes_vmax    = get_nodes_vmax(Nodes)

#     nodes_vmax_p  = perm_matrix' * nodes_vmax

#     nodes_vmin    = get_nodes_vmin(Nodes)

#     nodes_vmin_p  = perm_matrix' * nodes_vmin

#     nodes_name    = get_nodes_Bus_name(Nodes)

#     Ybus_p        = perm_matrix' * (Ybus_noperm) * perm_matrix

#     # ------------------------------------------------------
#     # Indices
#     # ------------------------------------------------------


#     SK_Idx_p      = idx_slack_nodes

#     PV_Idx_p      = idx_source_nodes   # gen nodes indices from perm 

#     PQ_Idx_p      = idx_demand_nodes   # load nodes indices from perm

#     PV_PQ_Idx_p   = first(PV_Idx_p):last(PQ_Idx_p)  

    

#     SK_PV_PQ_V_dim     = [length(SK_Idx_p),length(PV_Idx_p),length(PQ_Idx_p)]
#     SK_PV_PQ_V_offset  = create_offsets(SK_PV_PQ_V_dim)
#     SK_PV_PQ_V_Idx     = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

#     SK_V_Idx           = SK_PV_PQ_V_Idx[1]
#     PV_V_Idx           = SK_PV_PQ_V_Idx[2]
#     PQ_V_Idx           = SK_PV_PQ_V_Idx[3]

#     # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

#     x_dim = [length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx), length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx)]

#     x_offset   = create_offsets(x_dim)
#     x_Idx      = create_idxs(x_offset, x_dim)

#     x_Θ_SK_Idx = x_Idx[1]
#     x_Θ_PV_Idx = x_Idx[2]
#     x_Θ_PQ_Idx = x_Idx[3]
    
#     x_v_SK_Idx = x_Idx[4]
#     x_v_PV_Idx = x_Idx[5]
#     x_v_PQ_Idx = x_Idx[6]

#     x_Θ_Idx       = first(x_Θ_SK_Idx):last(x_Θ_PQ_Idx)
#     x_v_Idx       = first(x_v_SK_Idx):last(x_v_PQ_Idx)

#     x_Θ_PV_PQ_Idx = first(x_Θ_PV_Idx):last(x_Θ_PQ_Idx)
#     x_v_PV_PQ_Idx = first(x_v_PV_Idx):last(x_v_PQ_Idx)

#     Δx_dim = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PV_Idx),length(x_v_PQ_Idx)]

#     Δx_offset   = create_offsets( Δx_dim )
#     Δx_Idx      = create_idxs( Δx_offset, Δx_dim )

#     Δx_Θ_PV_Idx = Δx_Idx[1]
#     Δx_Θ_PQ_Idx = Δx_Idx[2]

#     Δx_v_PV_Idx = Δx_Idx[3]
#     Δx_v_PQ_Idx = Δx_Idx[4]

#     Δx_Θ_PV_PQ_Idx = first(Δx_Θ_PV_Idx):last(Δx_Θ_PQ_Idx)

#     Δx_v_PV_PQ_Idx = first(Δx_v_PV_Idx):last(Δx_v_PQ_Idx)

#     # ------------------------------------------------------

#     # computing arrays, indices and views
#     # ------------------------------------------------------

#     gen_Q      = zeros(Float64, num_gen)

#     gen_Q_view =  view(gen_Q, 1:num_gen)

#     Sg_view    = @view Sg[:]

#     update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

#     Sd_view    = @view Sd[:]
    
#     Y          = Ybus_p

#     # ------------------------------------------------------
#     # Parameteres
#     # ------------------------------------------------------
    
#     p_ΔP_ΔQ            = ( perm_matrix, Cg, Sd, Y )

#     p_update_gens_Q    = ( gens_bus_num, gen_Qmax, gen_Qmin )

#     p_x_Θ_Δx_Idx       = ( x_Θ_Idx,        x_v_Idx,
#                            x_Θ_PV_PQ_Idx,  x_v_PV_PQ_Idx,
#                            x_Θ_PV_Idx,     x_v_PV_Idx,
#                            x_Θ_PQ_Idx,     x_v_PQ_Idx,
                           
#                            Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx,
#                            Δx_Θ_PV_Idx,    Δx_v_PV_Idx,
#                            Δx_Θ_PQ_Idx,    Δx_v_PQ_Idx )

#     p_Idx = ( PV_Idx_p, PQ_Idx_p )
    
#     p_V_limits = ( nodes_vmax, nodes_vmin, gen_Vm, p_Idx )

#     p_slack = (slack_Vm, slack_Vθ,x_Θ_SK_Idx, x_v_SK_Idx )
    
#     p = ( p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack )
    
#     return dyn_powerflow_sol_param = (; Ybus_noperm, nodes_name, branches_name, num_nodes, perm_matrix, gen_Q_view, Sg_view, Sd_view, p )

# end



# function get_powerflow_data( case_data_func )

#     Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()


#     perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

#     num_nodes = get_network_nodal_size(Branches)

#     gens_S = get_nodes_S(Gen_nodes)

#     gen_Q = imag(gens_S) 

#     gens_bus_num = get_gens_bus_num(Gen_nodes)

#     Cg = get_slack_and_gens_Cg(Nodes)

#     # Cg = get_Cg(Gen_nodes, num_nodes)

#     Sg = get_slack_and_gen_Sg(Nodes)

#     # Sg = get_Sg(Gen_nodes, num_nodes)

#     gen_Vm = get_nodes_Vm(Gen_nodes) 

#     num_gen = length(gen_Vm)

#     gen_Qmax = get_nodes_Qmax(Gen_nodes)

#     gen_Qmin = get_nodes_Qmin(Gen_nodes)

#     Sd = get_loads_Sd(Load_nodes)

#     # branches_name = get_branches_name(Branches)
#     # nodes_name = get_nodes_Bus_name(Nodes)

#     branches_name = collect(keys(Branches))
    
#     nodes_name = collect(keys(Nodes))    
    
#     Yf = get_Yf(Nodes, Branches)
#     Yt = get_Yt(Nodes, Branches)

#     Ybus_noperm = get_Ybus(Nodes, Branches, shunts=Shunts)

#     slack_Vm = get_nodes_Vm(Slack_nodes)

#     slack_Vθ = get_nodes_Vθ(Slack_nodes)

#     nodes_vmax = get_nodes_vmax(Nodes)

#     nodes_vmax_p = perm_matrix' * nodes_vmax

#     nodes_vmin = get_nodes_vmin(Nodes)

#     nodes_vmin_p = perm_matrix' * nodes_vmin

#     Ybus_p = perm_matrix' * (Ybus_noperm) * perm_matrix
    
#     return powerflow_data = (; perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes, num_nodes, gens_S, gen_Q, gens_bus_num, Cg, Sg, gen_Vm, num_gen, gen_Qmax, gen_Qmin, Sd, branches_name, nodes_name, Yf, Yt, Ybus_noperm, slack_Vm, slack_Vθ, nodes_vmax, nodes_vmax_p, nodes_vmin, nodes_vmin_p,  Ybus_p )

# end



# function reduced_Jac_formula(x, (Y, p_Idx) )

#     V = x_to_V(x)

#     PV_Idx, PQ_Idx = p_Idx

#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

#     m11 = [[i==j ? ∂P_i_∂δ_i(i,V,Y) : ∂P_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

#     M11 = reshape( [m11...;], size(Y) )

#     n12 = [[i==j ? Vm_i_∂P_i_∂Vm_i(i,V,Y) : Vm_j_∂P_i_∂Vm_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

#     N12 = reshape( [n12...;], size(Y) )

#     n21 = [[i==j ? ∂Q_i_∂δ_i(i,V,Y) : ∂Q_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

#     N21 = reshape( [n21...;], size(Y) )

#     m22 = [[i==j ? Vm_i_∂Q_i_Vm_i(i,V,Y) : Vm_j_∂Q_i_Vm_j(i,j,V,Y) for j in 1:length(V)] for i in  1:length(V)] 

#     M22 = reshape( [m22...;], size(Y) )  

#     return vcat(hcat(M11'[PV_PQ_Idx, PV_PQ_Idx], N12'[PV_PQ_Idx, PQ_Idx]),
#                 hcat(N21'[PQ_Idx,    PV_PQ_Idx], M22'[PQ_Idx,    PQ_Idx]))
    
# end


# function Jac!(Jac, V, Y, p_Idx )

#     PV_Idx, PQ_Idx = p_Idx

#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

#     I = Y * V
    
#     Gs_θ = im * Diagonal(V) * (Diagonal(conj.(I)) - conj.(Y) * Diagonal(conj.(V)) )
    
#     Gs_v = Diagonal(V) * (Diagonal(conj.(I)) + conj.(Y) * Diagonal(conj.(V))) * Diagonal( inv.( abs.(V) ) ) 
    
#     Gs_P   = real.(hcat(Gs_θ[PV_PQ_Idx, PV_PQ_Idx], Gs_v[PV_PQ_Idx, PQ_Idx]))
#     Gs_Q   = imag.(hcat(Gs_θ[PQ_Idx,    PV_PQ_Idx], Gs_v[PQ_Idx,    PQ_Idx]))

#     Jac    .= vcat(Gs_P, Gs_Q)
    
#     return nothing

# end


# function Jac(V, Y, p_Idx )

#     PV_Idx, PQ_Idx = p_Idx

#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

#     I = Y * V
    
#     Gs_θ = im * Diagonal(V) * (Diagonal(conj.(I)) - conj.(Y) * Diagonal(conj.(V)) )
    
#     Gs_v = Diagonal(V) * (Diagonal(conj.(I)) + conj.(Y) * Diagonal(conj.(V))) * Diagonal( inv.( abs.(V) ) ) 
    
#     Gs_P   = real.(hcat(Gs_θ[PV_PQ_Idx, PV_PQ_Idx], Gs_v[PV_PQ_Idx, PQ_Idx]))
#     Gs_Q   = imag.(hcat(Gs_θ[PQ_Idx,    PV_PQ_Idx], Gs_v[PQ_Idx,    PQ_Idx]))
    
#     return vcat(Gs_P, Gs_Q)

# end

# function reduced_ΔP_ΔQ(V, Sg_view, p_ΔP_ΔQ, p_Idx)

#     perm_matrix, Cg, Sd, Y = p_ΔP_ΔQ

#     PV_Idx, PQ_Idx = p_Idx
    
#     PV_PQ_Idx      = first(PV_Idx):last(PQ_Idx)
    
#     mismatch       = Diagonal(V) * conj.(Y * V) - perm_matrix' * (Cg * Sg_view - Sd)

#     ΔP_ΔQ          = Array( [[ real.(mismatch[PV_Idx]),
#                                real.(mismatch[PQ_Idx]),
#                                imag.(mismatch[PQ_Idx])]...; ] )     

# end


# function reduced_pf_sol(x0, (gen_Q_view, Sg_view), p;
#                        maxiter=40, ftol=1000*eps(),
#                        xtol=1000*eps() )

#     p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx , p_slack = p

#     nodes_vmax, nodes_vmin, gen_Vm, p_Idx = p_V_limits

#     gens_bus_num, gen_Qmax, gen_Qmin  = p_update_gens_Q

#     perm_matrix, Cg, Sd, Y = p_ΔP_ΔQ 

#     PV_Idx_p, PQ_Idx_p = p_Idx

#     x_Θ_Idx, x_v_Idx, x_Θ_PV_PQ_Idx, x_v_PV_PQ_Idx, x_Θ_PV_Idx, x_v_PV_Idx, x_Θ_PQ_Idx, x_v_PQ_Idx, Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx, Δx_Θ_PV_Idx, Δx_v_PV_Idx, Δx_Θ_PQ_Idx, Δx_v_PQ_Idx = p_x_Θ_Δx_Idx

#     x_k_1              = deepcopy(x0)

#     x_k_1[x_v_PV_Idx] .= gen_Vm

#     x_k                = deepcopy(x0)

#     Δx                 = zeros(length( x_k ))

#     x_k[x_v_PV_Idx]   .= gen_Vm

#     V_k_1              = x_to_V(x_k_1, x_Θ_Idx, x_v_Idx)

#     V_k                = x_to_V(x_k, x_Θ_Idx, x_v_Idx)

#     ΔΘVm_dim = sum([length(x_Θ_PV_PQ_Idx), length(x_v_PV_PQ_Idx)])
    
#     ΔΘVm = Inf

#     ΔP_ΔQ_k = reduced_ΔP_ΔQ( V_k_1, Sg_view, p_ΔP_ΔQ, p_Idx )

#     local k  = 0
        
#     while (norm(ΔΘVm) > xtol) && (norm(ΔP_ΔQ_k) > ftol)

#         x_k_1[x_v_PV_Idx] .= gen_Vm

#         V_k_1              .=  x_to_V(x_k_1,  x_Θ_Idx, x_v_Idx)

#         Jac_k              = reduced_Jac_formula(V_k_1, Y,  p_Idx)

#         ΔP_ΔQ_k            = reduced_ΔP_ΔQ( V_k_1, Sg_view, p_ΔP_ΔQ, p_Idx ) 

#         ΔΘVm               =  -Jac_k \ ΔP_ΔQ_k

#         Δx[Δx_Θ_PV_PQ_Idx] = ΔΘVm[Δx_Θ_PV_PQ_Idx]

        
#         # if  reduced == true
            
#         if length(ΔΘVm) < ΔΘVm_dim # length(Δx)            
#             Δx[Δx_v_PQ_Idx] .= ΔΘVm[1+length(Δx_Θ_PV_PQ_Idx):end]
            
#         else
            
#             Δx[Δx_v_PV_Idx]  .= ΔΘVm[Δx_v_PV_Idx]
#             Δx[Δx_v_PQ_Idx]  .= ΔΘVm[Δx_v_PQ_Idx]
            
#         end

#         x_k[ x_Θ_PV_PQ_Idx ] .= x_k_1[ x_Θ_PV_PQ_Idx ] .+ Δx[ Δx_Θ_PV_PQ_Idx ]
#         x_k[ x_v_PQ_Idx ]    .= x_k_1[ x_v_PQ_Idx ] .* (1 .+ Δx[ Δx_v_PQ_Idx ])

#         # V_k =  x_to_V(x_k, x_Θ_Idx, x_v_Idx)

#         V_k  .=  x_to_V(x_k, x_Θ_Idx, x_v_Idx)

#         gen_Q = imag.((Diagonal(V_k) * conj.(Y * V_k) + perm_matrix' * Sd)[PV_Idx_p])

#         update_gens_Q!(Sg_view, gen_Q, p_update_gens_Q)

#         x_k_1 .= x_k

#         k += 1
        
#         if k == maxiter
#             @warn "Maximum number of iterations reached."

#             break
#         end

#     end

#     return k, V_k, ΔP_ΔQ_k, Δx   
# end



# #------------------------------------------------------#
# # ------------------------------------------------------
# #  Power flow based on current balance
# # ------------------------------------------------------
# #------------------------------------------------------#


# function get_full_current_powerflow_sol_param(powerflow_data)

#     perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes, num_nodes, gens_S, gen_Q, gens_bus_num, Cg, Sg, gen_Vm, num_gen, gen_Qmax, gen_Qmin, Sd, branches_name, nodes_name, Yf, Yt, Ybus_noperm, slack_Vm, slack_Vθ, nodes_vmax, nodes_vmax_p, nodes_vmin, nodes_vmin_p, nodes_name, Ybus_p = powerflow_dat

#     # ------------------------------------------------------
#     # Indices
#     # ------------------------------------------------------

#     SK_Idx = idx_slack_nodes

#     PV_Idx = idx_source_nodes  

#     PQ_Idx = idx_demand_nodes 

#     SK_PV_PQ_V_dim    = [length(SK_Idx),length(PV_Idx),length(PQ_Idx)]
#     SK_PV_PQ_V_offset = create_offsets(SK_PV_PQ_V_dim)
#     SK_PV_PQ_V_Idx    = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

#     SK_V_Idx = SK_PV_PQ_V_Idx[1]
#     PV_V_Idx = SK_PV_PQ_V_Idx[2]
#     PQ_V_Idx = SK_PV_PQ_V_Idx[3]

#     # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

#     x_dim      = [length(PV_V_Idx), length(PQ_V_Idx), length(PV_V_Idx), length(PQ_V_Idx) ]

#     x_offset   = create_offsets(x_dim)
#     x_Idx      = create_idxs(x_offset, x_dim)

#     x_Θ_PV_Idx = x_Idx[1]
#     x_Θ_PQ_Idx = x_Idx[2]
#     x_v_PV_Idx = x_Idx[3]
#     x_v_PQ_Idx = x_Idx[4]

#     Δx_dim = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_Θ_PV_Idx), length(x_v_PQ_Idx)]

#     Δx_offset   = create_offsets(Δx_dim )
#     Δx_Idx      = create_idxs(Δx_offset, Δx_dim)

#     Δx_Θ_PV_Idx = Δx_Idx[1]
#     Δx_Θ_PQ_Idx = Δx_Idx[2]

#     Δx_v_PV_Idx = Δx_Idx[3]
#     Δx_v_PQ_Idx = Δx_Idx[4]

#     # ------------------------------------------------------
#     # computing arrays, indices and views
#     # ------------------------------------------------------

#     # gen_Q      = zeros(Float64, num_gen)

#     gen_Q_view =  view(gen_Q, 1:num_gen)

#     Sg_view    = @view Sg[:]

#     update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

#     Sd_view    = @view Sd[:]
    
#     Y          = Ybus_p

#     # ------------------------------------------------------
#     # Parameteres
#     # ------------------------------------------------------
    
#     p_ΔI               = (perm_matrix, Cg, Sd, Y )

#     p_update_gens_Q    = (gens_bus_num, gen_Qmax, gen_Qmin )

#     p_Idx              = (SK_Idx, PV_Idx, PQ_Idx)

#     p_V_limits         = (nodes_vmax, nodes_vmin, gen_Vm, p_Idx)
    
#     p_x_Θ_Δx_Idx       = ( x_Θ_PV_Idx,
#                            x_Θ_PQ_Idx,
#                            x_v_PV_Idx,
#                            x_v_PQ_Idx,
                           
#                            Δx_Θ_PV_Idx,
#                            Δx_Θ_PQ_Idx,
#                            Δx_v_PV_Idx,
#                            Δx_v_PQ_Idx )

#     p                  = (p_V_limits,
#                           p_update_gens_Q,
#                           p_ΔI,
#                           p_Idx,
#                           p_x_Θ_Δx_Idx)

#     return powerflow_sol_param = (;num_nodes, perm_matrix, gen_Q_view, Sg_view, Sd_view, p)


# end




# function get_reduced_current_powerflow_sol_param(powerflow_data)

#     perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes, num_nodes, gens_S, gen_Q, gens_bus_num, Cg, Sg, gen_Vm, num_gen, gen_Qmax, gen_Qmin, Sd, branches_name, nodes_name, Yf, Yt, Ybus_noperm, slack_Vm, slack_Vθ, nodes_vmax, nodes_vmax_p, nodes_vmin, nodes_vmin_p, nodes_name, Ybus_p = powerflow_dat

#     # ------------------------------------------------------
#     # Indices
#     # ------------------------------------------------------

#     SK_Idx = idx_slack_nodes

#     PV_Idx = idx_source_nodes  

#     PQ_Idx = idx_demand_nodes 

#     SK_PV_PQ_V_dim    = [length(SK_Idx),length(PV_Idx),length(PQ_Idx)]
#     SK_PV_PQ_V_offset = create_offsets(SK_PV_PQ_V_dim)
#     SK_PV_PQ_V_Idx    = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

#     SK_V_Idx = SK_PV_PQ_V_Idx[1]
#     PV_V_Idx = SK_PV_PQ_V_Idx[2]
#     PQ_V_Idx = SK_PV_PQ_V_Idx[3]

#     # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

#     x_dim      = [length(PV_V_Idx), length(PQ_V_Idx), length(PQ_V_Idx)]

#     x_offset   = create_offsets(x_dim)
#     x_Idx      = create_idxs(x_offset, x_dim)

#     x_Θ_PV_Idx = x_Idx[1]
#     x_Θ_PQ_Idx = x_Idx[2]    
#     x_v_PQ_Idx = x_Idx[3]

#     Δx_dim = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_Θ_PV_Idx), length(x_v_PQ_Idx)]

#     Δx_offset   = create_offsets(Δx_dim )
#     Δx_Idx      = create_idxs(Δx_offset, Δx_dim)

#     Δx_Θ_PV_Idx = Δx_Idx[1]
#     Δx_Θ_PQ_Idx = Δx_Idx[2]

#     Δx_v_PV_Idx = Δx_Idx[3]
#     Δx_v_PQ_Idx = Δx_Idx[4]

#     # ------------------------------------------------------
#     # computing arrays, indices and views
#     # ------------------------------------------------------

#     # gen_Q      = zeros(Float64, num_gen)

#     gen_Q_view =  view(gen_Q, 1:num_gen)

#     Sg_view    = @view Sg[:]

#     update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

#     Sd_view    = @view Sd[:]
    
#     Y          = Ybus_p

#     # ------------------------------------------------------
#     # Parameteres
#     # ------------------------------------------------------
    
#     p_ΔI               = (perm_matrix, Cg, Sd, Y )

#     p_update_gens_Q    = (gens_bus_num, gen_Qmax, gen_Qmin )

#     p_Idx              = (PV_Idx, PQ_Idx)

#     p_V_limits         = (nodes_vmax, nodes_vmin, gen_Vm, p_Idx)
    
#     p_x_Θ_Δx_Idx       = ( x_Θ_PV_Idx,
#                            x_Θ_PQ_Idx,
#                            x_v_PQ_Idx,
                           
#                            Δx_Θ_PV_Idx,
#                            Δx_Θ_PQ_Idx,
#                            Δx_v_PV_Idx,
#                            Δx_v_PQ_Idx )

#     p                  = (p_V_limits,
#                           p_update_gens_Q,
#                           p_ΔI,
#                           p_Idx,
#                           p_x_Θ_Δx_Idx)

#     return powerflow_sol_param = (;num_nodes, perm_matrix, gen_Q_view, Sg_view, Sd_view, p)


# end



# function Jac_Ibus(V, (Y, p_Idx) )

#     PV_Idx, PQ_Idx = p_Idx

#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

#     # E =  exp.(im * angle.(V) )
    
#     Ibus_θ = im * Y * Diagonal(V) 
    
#     Ibus_v =  Y * Diagonal( exp.(im * angle.(V) ) ) 
    
#     Gs_M = hcat(real.(Ibus_θ[PV_PQ_Idx, PV_PQ_Idx]), real.(Ibus_v[PV_PQ_Idx, PV_PQ_Idx]))
    
#     Gs_N = hcat(imag.(Ibus_θ[PV_PQ_Idx, PV_PQ_Idx]), imag.(Ibus_v[PV_PQ_Idx, PV_PQ_Idx]))
    
#     return vcat(Gs_M, Gs_N)
    

# end



# function Jac_Ibus!(Jac, V, (Y, p_Idx) )

#     PV_Idx, PQ_Idx = p_Idx

#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

#     # E =  exp.(im * angle.(V) )
    
#     Ibus_θ = im * Y * Diagonal(V) 
    
#     Ibus_v =  Y * Diagonal( exp.(im * angle.(V) ) ) 
    
#     Gs_M = real.(hcat(Ibus_θ[PV_PQ_Idx, PV_PQ_Idx], Ibus_v[PV_PQ_Idx, PV_PQ_Idx]))
#     Gs_N = imag.(hcat(Ibus_θ[PV_PQ_Idx, PV_PQ_Idx], Ibus_v[PV_PQ_Idx, PV_PQ_Idx]))
    
#     Jac    .= vcat(Gs_M, Gs_N)

#     nothing
    

# end



# function Jac_Ibus(x, Y, (slack_Vm, slack_Vθ, gen_Vm, p_Idx, p_x_Θ_Idx) )

#     PV_Idx, PQ_Idx = p_Idx

#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

#     x_Θ_PV_Idx, x_Θ_PQ_Idx, x_v_PV_Idx, x_v_PQ_Idx = p_x_Θ_Idx

#     x_Θ_PV = x[x_Θ_PV_Idx]
#     x_Θ_PQ = x[x_Θ_PQ_Idx]
#     x_v_PV = x[x_v_PQ_Idx]
#     x_v_PQ = x[x_v_PQ_Idx]

#     x_v_PV .= gen_Vm

#     V = [[slack_Vm .* exp.(im * slack_Vθ)]...;[x_v_PV .* exp.(im * x_Θ_PV)]...; [x_v_PQ .* exp.(im * x_Θ_PQ)]...]    

#     # E =  exp.(im * angle.(V) )
    
#     Ibus_θ = im * Y * Diagonal(V) 
    
#     Ibus_v =  Y * Diagonal( exp.(im * angle.(V) ) ) 
    
#     Gs_M = real.(hcat(Ibus_θ[PV_PQ_Idx, PV_PQ_Idx], Ibus_v[PV_PQ_Idx, PV_PQ_Idx]))
#     Gs_N = imag.(hcat(Ibus_θ[PV_PQ_Idx, PV_PQ_Idx], Ibus_v[PV_PQ_Idx, PV_PQ_Idx]))
    
#     return vcat(Gs_M, Gs_N)

# end


# function Jac_Ibus!(Jac, x, Y, (slack_Vm, slack_Vθ, gen_Vm, p_Idx, p_x_Θ_Idx) )

#     PV_Idx, PQ_Idx = p_Idx

#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

#     x_Θ_PV_Idx, x_Θ_PQ_Idx, x_v_PV_Idx, x_v_PQ_Idx = p_x_Θ_Idx

#     x_Θ_PV = x[x_Θ_PV_Idx]
#     x_Θ_PQ = x[x_Θ_PQ_Idx]
#     x_v_PV = x[x_v_PQ_Idx]
#     x_v_PQ = x[x_v_PQ_Idx]

#     x_v_PV .= gen_Vm

#     V = [[slack_Vm .* exp.(im * slack_Vθ)]...;[x_v_PV .* exp.(im * x_Θ_PV)]...; [x_v_PQ .* exp.(im * x_Θ_PQ)]...]    

#     # E =  exp.(im * angle.(V) )
    
#     Ibus_θ = im * Y * Diagonal(V) 
    
#     Ibus_v =  Y * Diagonal( exp.(im * angle.(V) ) ) 
    
#     Gs_M = real.(hcat(Ibus_θ[PV_PQ_Idx, PV_PQ_Idx], Ibus_v[PV_PQ_Idx, PV_PQ_Idx]))
#     Gs_N = imag.(hcat(Ibus_θ[PV_PQ_Idx, PV_PQ_Idx], Ibus_v[PV_PQ_Idx, PV_PQ_Idx]))

#     Jac    .= vcat(Gs_M, Gs_N)
    
#     nothing

# end



# function full_ΔI(V, (Sg_view, Sd_view, p_ΔI, p_Idx))

#     perm_matrix, Cg, _,Y = p_ΔI

#     PV_Idx, PQ_Idx = p_Idx
    
#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)
    
#     mismatch  = Y * V + Diagonal(inv.(conj.(V))) * ( perm_matrix' * conj.( Sd_view - Cg * Sg_view ) )

#     ΔI = [[real.(mismatch[PV_PQ_Idx]), imag.(mismatch[PV_PQ_Idx])]...;] 

# end


# function full_ΔI(x, Sg_view, Sd_view, (slack_Vm, slack_Vθ, gen_Vm, p_ΔI, p_Idx, p_x_Θ_Idx))

#     perm_matrix, Cg, _,Y = p_ΔI

#     PV_Idx, PQ_Idx = p_Idx

#     x_Θ_PV_Idx, x_Θ_PQ_Idx, x_v_PV_Idx, x_v_PQ_Idx = p_x_Θ_Idx
    
#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

#     x_Θ_PV = x[x_Θ_PV_Idx]
#     x_Θ_PQ = x[x_Θ_PQ_Idx]
#     x_v_PV = x[x_v_PQ_Idx]
#     x_v_PQ = x[x_v_PQ_Idx]

#     x_v_PV .= gen_Vm

#     V = [[slack_Vm .* exp.(im * slack_Vθ)]...;[x_v_PV .* exp.(im * x_Θ_PV)]...; [x_v_PQ .* exp.(im * x_Θ_PQ)]...]
    
#     mismatch  = Y * V + Diagonal(inv.(conj.(V))) * ( perm_matrix' * conj.( Sd_view - Cg * Sg_view ) )

#     ΔI = [[real.(mismatch[PV_PQ_Idx]), imag.(mismatch[PV_PQ_Idx])]...;] 

# end



# function full_ΔI!( ΔI, V, (Sg_view, Sd_view, p_ΔI, p_Idx))

#     perm_matrix, Cg, _,Y = p_ΔI

#     PV_Idx, PQ_Idx = p_Idx
    
#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)
    
#     mismatch  = Y * V + Diagonal(inv.(conj.(V))) * ( perm_matrix' * conj.( Sd_view - Cg * Sg_view ) )

#     ΔI = [[real.(mismatch[PV_PQ_Idx]), imag.(mismatch[PV_PQ_Idx])]...;] 

# end


# function full_ΔI!( ΔI, x, Sg_view, Sd_view, (slack_Vm, slack_Vθ, gen_Vm, p_ΔI, p_Idx, p_x_Θ_Idx))

#     perm_matrix, Cg, _,Y = p_ΔI

#     PV_Idx, PQ_Idx = p_Idx

#     x_Θ_PV_Idx, x_Θ_PQ_Idx, x_v_PV_Idx, x_v_PQ_Idx = p_x_Θ_Idx
    
#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)

#     x_Θ_PV = x[x_Θ_PV_Idx]
#     x_Θ_PQ = x[x_Θ_PQ_Idx]
#     x_v_PV = x[x_v_PQ_Idx]
#     x_v_PQ = x[x_v_PQ_Idx]

#     x_v_PV .= gen_Vm

#     V = [[slack_Vm .* exp.(im * slack_Vθ)]...;[x_v_PV .* exp.(im * x_Θ_PV)]...; [x_v_PQ .* exp.(im * x_Θ_PQ)]...]
    
#     mismatch  = Y * V + Diagonal(inv.(conj.(V))) * ( perm_matrix' * conj.( Sd_view - Cg * Sg_view ) )

#     ΔI = [[real.(mismatch[PV_PQ_Idx]), imag.(mismatch[PV_PQ_Idx])]...;] 

# end


# function full_current_pf_sol(x0,
#                        (gen_Q_view, Sg_view, Sd_view), p;
#                        maxiter=40,
#                        ftol=1000*eps(),
#                        xtol=1000*eps() )

#     # p_update_Jac, p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx = p
    
#     p_update_Jac, p_V_limits, p_update_gens_Q, p_ΔI, p_Idx, p_x_Θ_Δx_Idx = p

#     Jac_Q_PV_row_Idx, _ = p_update_Jac

#     nodes_vmax, nodes_vmin, gen_Vm, p_Idx = p_V_limits

#     gens_bus_num, gen_Qmax, gen_Qmin  = p_update_gens_Q

#     perm_matrix, Cg, Sd, Y = p_ΔI 

#     PV_Idx_p, PQ_Idx_p = p_Idx

#     x_Θ_Idx, x_v_Idx, x_Θ_PV_PQ_Idx, x_v_PV_PQ_Idx, x_Θ_PV_Idx, x_v_PV_Idx, x_Θ_PQ_Idx, x_v_PQ_Idx, Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx, Δx_Θ_PV_Idx, Δx_v_PV_Idx, Δx_Θ_PQ_Idx, Δx_v_PQ_Idx = p_x_Θ_Δx_Idx

#     x_k_1              = deepcopy(x0)

#     x_k_1[x_v_PV_Idx] .= gen_Vm

#     x_k                = deepcopy(x0)

#     Δx                 = zeros(length( x_k ))

#     # @show Δx
    
#     x_k[x_v_PV_Idx]   .= gen_Vm

#     V_k_1              = x_to_V(x_k_1, x_Θ_Idx, x_v_Idx)

#     V_k                = x_to_V(x_k,   x_Θ_Idx, x_v_Idx)

#     ΔΘVm_dim           = sum([length(x_Θ_PV_PQ_Idx), length(x_v_PV_PQ_Idx)])
    
#     ΔΘVm               = ones( ΔΘVm_dim ) # Inf

#     ΔI_k               =  full_ΔI(V_k_1, (Sg_view, Sd_view, p_ΔI, p_Idx))
    

#     local k  = 1
        
#     while (norm(ΔΘVm) > xtol) && (norm(ΔI_k) > ftol)

#         x_k_1[x_v_PV_Idx] .= gen_Vm

#         V_k_1             .=  x_to_V(x_k_1,  x_Θ_Idx, x_v_Idx)

#         Jac_k              = Jac_Ibus(V_k_1, (Y, p_Idx) )
        
#         ΔI_k               =  full_ΔI(V_k_1, (Sg_view, Sd_view, p_ΔI, p_Idx))

#         ΔΘVm               .= -Jac_k \ ΔI_k
        
#         Δx[Δx_Θ_PV_PQ_Idx] .= ΔΘVm[Δx_Θ_PV_PQ_Idx]
        
#         #
#         Δx[Δx_v_PV_PQ_Idx] .= ΔΘVm[Δx_v_PV_PQ_Idx]

#         x_k[ x_Θ_PV_PQ_Idx ] .= x_k_1[ x_Θ_PV_PQ_Idx ] .+ Δx[ Δx_Θ_PV_PQ_Idx ]

#         #
#         x_k[ x_v_PV_PQ_Idx ] .= x_k_1[ x_v_PV_PQ_Idx ] .+ Δx[ Δx_v_PV_PQ_Idx ]
        
#         V_k  .=  x_to_V(x_k, x_Θ_Idx, x_v_Idx)

#         # ------------------------------------------------------
#         # @show V_k[1]

#         gen_Q = imag.((Diagonal(V_k) * conj.(Y * V_k) + perm_matrix' * Sd)[PV_Idx_p])

#         # ------------------------------------------------------
        
#         Q_nolimit_boolean      = map((x) -> no_limit_violation(x...),
#                                      zip(gen_Q, gen_Qmax, gen_Qmin))

#         Q_limit_boolean        = map((x) -> limit_violation(x...),
#                                      zip(gen_Q, gen_Qmax, gen_Qmin))

#         Q_nolimit_violation    = gen_Q[Q_nolimit_boolean]

#         Q_limit_violation      = gen_Q[Q_limit_boolean]

#         if length(Q_limit_violation) != 0

#             gen_Q .= threshold_limits.(gen_Q, gen_Qmax, gen_Qmin)
            
#             update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q)

#         else

#             update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q)
            
#             # x_v_PV_restriction_Idx = x_v_PV_Idx[Q_nolimit_boolean]
#             # PV_Q_restriction_Idx   = Jac_Q_PV_row_Idx[Q_nolimit_boolean]

#             # for (qh_Idx, vh_Idx) in zip(PV_Q_restriction_Idx, x_v_PV_restriction_Idx )
#             #     Jac_k[qh_Idx, :     ] .= 0.0
#             #     Jac_k[:,      vh_Idx] .= 0.0
#             #     Jac_k[qh_Idx, vh_Idx] = 1.0
#             # end
            
#             Jac_k    = Jac_Ibus(V_k, (Y, p_Idx) )
        
#             ΔI_k   =  full_ΔI(V_k_1, (Sg_view, Sd_view, p_ΔI, p_Idx))

#             ΔΘVm  .= -Jac_k \ ΔI_k

#             Δx[Δx_Θ_PV_PQ_Idx] .= ΔΘVm[Δx_Θ_PV_PQ_Idx]

#             #
#             # Δx[Δx_v_PV_PQ_Idx] .= ΔΘVm[Δx_v_PV_PQ_Idx]
            
#             Δx[Δx_v_PV_Idx]  .= ΔΘVm[Δx_v_PV_Idx]
#             Δx[Δx_v_PQ_Idx]  .= ΔΘVm[Δx_v_PQ_Idx]

#             x_k[x_Θ_PV_PQ_Idx] .= x_k_1[ x_Θ_PV_PQ_Idx ] .+ Δx[ Δx_Θ_PV_PQ_Idx ]

#             #
            
#             x_k[x_v_PV_PQ_Idx] .= x_k_1[ x_v_PV_PQ_Idx ] .+ Δx[ Δx_v_PV_PQ_Idx ]
            
#             V_k                 .=  x_to_V(x_k, x_Θ_Idx, x_v_Idx)
#         end

#         update_gens_Q!(Sg_view, gen_Q, p_update_gens_Q)

#         # ------------------------------------------------------

#         x_k_1 .= x_k

#         k += 1
        
#         if k == maxiter
#             @warn "Maximum number of iterations reached."

#             break
#         end

#     end 

#     return k, V_k, ΔI_k, Δx
    
#     # return k, V_k, ΔP_ΔQ_k, Δx   
# end


# function full_current_powerflow( case_data_func )

#     Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

#     perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

#     num_nodes     = get_network_nodal_size(Branches)

#     gens_S        = get_nodes_S(Gen_nodes)

#     gen_Q         = imag(gens_S) 

#     gens_bus_num  = get_gens_bus_num(Gen_nodes)

#     Cg            = get_Cg(Gen_nodes, num_nodes)

#     Sg            = get_Sg(Gen_nodes, num_nodes)

#     gen_Vm        = get_nodes_Vm(Gen_nodes) 

#     num_gen       =  length(gen_Vm)

#     gen_Qmax      = get_nodes_Qmax(Gen_nodes)

#     gen_Qmin      = get_nodes_Qmin(Gen_nodes)

#     Sd            = get_loads_Sd(Load_nodes)

#     Ybus_noperm   = get_Ybus(Nodes, Branches, shunts=Shunts)

#     slack_Vm      = get_nodes_Vm(Slack_nodes)

#     slack_Vθ      = get_nodes_Vθ(Slack_nodes)

#     nodes_vmax    = get_nodes_vmax(Nodes)

#     nodes_vmax_p  = perm_matrix' * nodes_vmax

#     nodes_vmin    = get_nodes_vmin(Nodes)

#     nodes_vmin_p  = perm_matrix' * nodes_vmin

#     nodes_name    =  get_nodes_Bus_name(Nodes)

#     Ybus_p        = perm_matrix' * (Ybus_noperm) * perm_matrix

#     # PV_Idx        = idx_source_nodes   # gen nodes indices from perm 
#     # PQ_Idx        = idx_demand_nodes   # load nodes indices from perm

#     SK_Idx_p      = idx_slack_nodes

#     PV_Idx_p      = idx_source_nodes   # gen nodes indices from perm 

#     PQ_Idx_p      = idx_demand_nodes   # load nodes indices from perm

#     PV_PQ_Idx_p   = first(PV_Idx_p):last(PQ_Idx_p)   # [PV_Idx;PQ_Idx]

#     # ------------------------------------------------------
#     # Indices
#     # ------------------------------------------------------

#     # Note slack bus is remove, hence PV indices start at 1.

#     SK_PV_PQ_V_dim     = [length(SK_Idx_p),length(PV_Idx_p),length(PQ_Idx_p)]
#     SK_PV_PQ_V_offset  = create_offsets(SK_PV_PQ_V_dim)
#     SK_PV_PQ_V_Idx     = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

#     SK_V_Idx           = SK_PV_PQ_V_Idx[1]
#     PV_V_Idx           = SK_PV_PQ_V_Idx[2]
#     PQ_V_Idx           = SK_PV_PQ_V_Idx[3]

#     # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

#     x_dim      = [length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx),
#                   length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx)]

#     x_offset   = create_offsets(x_dim)
#     x_Idx      = create_idxs(x_offset, x_dim)

#     x_Θ_SK_Idx = x_Idx[1]
#     x_Θ_PV_Idx = x_Idx[2]
#     x_Θ_PQ_Idx = x_Idx[3]
#     x_v_SK_Idx = x_Idx[4]
#     x_v_PV_Idx = x_Idx[5]
#     x_v_PQ_Idx = x_Idx[6]

#     x_Θ_Idx       = first(x_Θ_SK_Idx):last(x_Θ_PQ_Idx)
#     x_v_Idx       = first(x_v_SK_Idx):last(x_v_PQ_Idx)

#     x_Θ_PV_PQ_Idx = first(x_Θ_PV_Idx):last(x_Θ_PQ_Idx)
#     x_v_PV_PQ_Idx = first(x_v_PV_Idx):last(x_v_PQ_Idx)

#     Δx_dim      = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PV_Idx),length(x_v_PQ_Idx)]

#     Δx_offset   = create_offsets( Δx_dim )
#     Δx_Idx      = create_idxs( Δx_offset, Δx_dim )

#     Δx_Θ_PV_Idx = Δx_Idx[1]
#     Δx_Θ_PQ_Idx = Δx_Idx[2]

#     Δx_v_PV_Idx = Δx_Idx[3]
#     Δx_v_PQ_Idx = Δx_Idx[4]

#     Δx_Θ_PV_PQ_Idx = first(Δx_Θ_PV_Idx):last(Δx_Θ_PQ_Idx)

#     Δx_v_PV_PQ_Idx = first(Δx_v_PV_Idx):last(Δx_v_PQ_Idx)


#     Jac_dim      = [length(PV_Idx_p), length(PQ_Idx_p),length(PV_Idx_p), length(PQ_Idx_p)]
#     Jac_offset   = create_offsets(Jac_dim )
#     Jac_Idx      = create_idxs(Jac_offset, Jac_dim)

#     Jac_P_PV_row_Idx  = Jac_Idx[1]
#     Jac_P_PQ_row_Idx  = Jac_Idx[2]
#     Jac_Q_PV_row_Idx  = Jac_Idx[3]
#     Jac_Q_PQ_row_Idx  = Jac_Idx[4]

    
#     # ------------------------------------------------------

#     # computing arrays, indices and views
#     # ------------------------------------------------------

#     gen_Q      = zeros(Float64, num_gen)

#     gen_Q_view =  view(gen_Q, 1:num_gen)

#     Sg_view    = @view Sg[:]

#     update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

#     Y          = Ybus_p

#     Sd_view    =  @view Sd[:]
    
#     # ------------------------------------------------------
#     # Parameteres
#     # ------------------------------------------------------

#     p_Idx              = ( PV_Idx_p, PQ_Idx_p )

#     p_ΔI               = ( perm_matrix, Cg, Sd, Y )
#     # p_ΔP_ΔQ            = ( perm_matrix, Cg, Sd, Y )

#     p_update_gens_Q    = ( gens_bus_num, gen_Qmax, gen_Qmin )

#     p_x_Θ_Δx_Idx       = ( x_Θ_Idx,        x_v_Idx,
#                            x_Θ_PV_PQ_Idx,  x_v_PV_PQ_Idx,
#                            x_Θ_PV_Idx,     x_v_PV_Idx,
#                            x_Θ_PQ_Idx,     x_v_PQ_Idx,
#                            Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx,
#                            Δx_Θ_PV_Idx,    Δx_v_PV_Idx,
#                            Δx_Θ_PQ_Idx,    Δx_v_PQ_Idx )

#     p_Idx              = (PV_Idx_p, PQ_Idx_p)
    
#     p_V_limits         = ( nodes_vmax, nodes_vmin, gen_Vm, p_Idx )

#     p_update_Jac       = ( Jac_Q_PV_row_Idx, x_v_PV_Idx )


#     # p                  = ( p_update_Jac, p_V_limits,
#     #                       p_update_gens_Q, p_ΔP_ΔQ,
#     #                       p_Idx, p_x_Θ_Δx_Idx )

#     p_slack = (slack_Vm, slack_Vθ,x_Θ_SK_Idx, x_v_SK_Idx)
    
#     p                  = ( p_update_Jac, p_V_limits,
#                           p_update_gens_Q, p_ΔI,
#                           p_Idx, p_x_Θ_Δx_Idx )   

#     # ------------------------------------------------------
#     # Initial condition
#     # ------------------------------------------------------

#     x_Θ                = zeros(num_nodes)

#     x_v                = ones(num_nodes)

#     x0                 = [x_Θ...; x_v...]

#     # x_k_1              = deepcopy(x0)
#     # x_k_1[x_v_PV_Idx] .= gen_Vm
#     # x_k                = deepcopy(x0)
#     # Δx                 = zeros(length( x_k ))
#     # x_k[x_v_PV_Idx]   .= gen_Vm
#     # V_k                = x_to_V(x_k, x_Θ_Idx, x_v_Idx)

#     k, V_k, ΔI_k, Δx = full_current_pf_sol(x0, (gen_Q_view, Sg_view, Sd_view), p;
#                                       maxiter=40,
#                                       ftol=1000*eps(),
#                                       xtol=1000*eps() )


#     V         = perm_matrix * V_k

#     Sbus      = Diagonal(V) * conj.(Ybus_noperm * V)

#     GenSinj   = Diagonal(V) * conj.(Ybus_noperm * V) + Sd

#     If         = get_Yf(Nodes, Branches) * V
    
#     It         = get_Yt(Nodes, Branches) * V

#     Ibranches  = If + It

#     Iinj = nodal_current_injection(V, p_ΔI)

#     # branches_name  = get_branches_name(Branches)  
#     # nodes_name     = get_nodes_Bus_name(Nodes)

#     branches_name  = collect(keys(Branches))
    
#     nodes_name     = collect(keys(Nodes))    

#     branch_dict_init = OrderedDict(name =>
#         ( real(i_f),imag(i_f), real(i_t), imag(i_t) ) for (name, i_f, i_t) in zip(branches_name, If, It))

#     bus_dict_init = OrderedDict(name =>
#         (vh, θh, ph, qh) for (name,
#                               vh,
#                               θh,
#                               ph,
#                               qh) in zip(nodes_name,
#                                          round.( abs.(V); digits=4),
#                                          round.( angle.(V); digits=4),
#                                          round.( real.(Sbus); digits=4),
#                                          round.( imag.(Sbus); digits=4) ))      
    
    
#     return Dict("Iinj" => Iinj, "Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V), "Vθ" => angle.(V), "Vbus" => V, "Ibranches" => Ibranches, "Ybus" => Ybus_noperm , "Sbus" => Sbus, "GenSinj" => GenSinj, "dict_init" =>Dict("bus_dict_init" => bus_dict_init,  "branch_dict_init" =>  branch_dict_init))
    
    
# end




#########################################

# function x_get_dynamic_reduced_powerflow( case_data_func )

#     Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

#     perm_matrix, idx_perm_nodes_idx_and_type, idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

#     num_nodes     = get_network_nodal_size(Branches)

#     gens_S        = get_nodes_S(Gen_nodes)

#     gen_Q         = imag(gens_S) 

#     gens_bus_num  = get_gens_bus_num(Gen_nodes)

#     Cg            = get_Cg(Gen_nodes, num_nodes)

#     Sg            = get_Sg(Gen_nodes, num_nodes)

#     gen_Vm        = get_nodes_Vm(Gen_nodes) 

#     num_gen       =  length(gen_Vm)

#     gen_Qmax      = get_nodes_Qmax(Gen_nodes)

#     gen_Qmin      = get_nodes_Qmin(Gen_nodes)

#     Sd            = get_loads_Sd(Load_nodes)


#     # branches_name  = get_branches_name(Branches)   
#     # nodes_name     = get_nodes_Bus_name(Nodes)

#     branches_name  = collect(keys(Branches))
    
#     nodes_name     = collect(keys(Nodes))      
    
#     Yf             = get_Yf(Nodes, Branches)
    
#     Yt             = get_Yt(Nodes, Branches)
    

#     Ybus_noperm   = get_Ybus(Nodes, Branches, shunts=Shunts)

#     slack_Vm      = get_nodes_Vm(Slack_nodes)

#     slack_Vθ      = get_nodes_Vθ(Slack_nodes)

#     nodes_vmax    = get_nodes_vmax(Nodes)

#     nodes_vmax_p  = perm_matrix' * nodes_vmax

#     nodes_vmin    = get_nodes_vmin(Nodes)

#     nodes_vmin_p  = perm_matrix' * nodes_vmin

#     nodes_name    = get_nodes_Bus_name(Nodes)

#     Ybus_p        = perm_matrix' * (Ybus_noperm) * perm_matrix

#     # ------------------------------------------------------
#     # Indices
#     # ------------------------------------------------------


#     SK_Idx_p      = idx_slack_nodes

#     PV_Idx_p      = idx_source_nodes   # gen nodes indices from perm 

#     PQ_Idx_p      = idx_demand_nodes   # load nodes indices from perm

#     PV_PQ_Idx_p   = first(PV_Idx_p):last(PQ_Idx_p)  

    

#     SK_PV_PQ_V_dim     = [length(SK_Idx_p),length(PV_Idx_p),length(PQ_Idx_p)]
#     SK_PV_PQ_V_offset  = create_offsets(SK_PV_PQ_V_dim)
#     SK_PV_PQ_V_Idx     = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

#     SK_V_Idx           = SK_PV_PQ_V_Idx[1]
#     PV_V_Idx           = SK_PV_PQ_V_Idx[2]
#     PQ_V_Idx           = SK_PV_PQ_V_Idx[3]

#     # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

#     x_dim      = [length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx),
#                   length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx)]

#     x_offset   = create_offsets(x_dim)
#     x_Idx      = create_idxs(x_offset, x_dim)

#     x_Θ_SK_Idx = x_Idx[1]
#     x_Θ_PV_Idx = x_Idx[2]
#     x_Θ_PQ_Idx = x_Idx[3]
    
#     x_v_SK_Idx = x_Idx[4]
#     x_v_PV_Idx = x_Idx[5]
#     x_v_PQ_Idx = x_Idx[6]

#     x_Θ_Idx       = first(x_Θ_SK_Idx):last(x_Θ_PQ_Idx)
#     x_v_Idx       = first(x_v_SK_Idx):last(x_v_PQ_Idx)

#     x_Θ_PV_PQ_Idx = first(x_Θ_PV_Idx):last(x_Θ_PQ_Idx)
#     x_v_PV_PQ_Idx = first(x_v_PV_Idx):last(x_v_PQ_Idx)

#     Δx_dim = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PV_Idx),length(x_v_PQ_Idx)]

#     Δx_offset   = create_offsets( Δx_dim )
#     Δx_Idx      = create_idxs( Δx_offset, Δx_dim )

#     Δx_Θ_PV_Idx = Δx_Idx[1]
#     Δx_Θ_PQ_Idx = Δx_Idx[2]

#     Δx_v_PV_Idx = Δx_Idx[3]
#     Δx_v_PQ_Idx = Δx_Idx[4]

#     Δx_Θ_PV_PQ_Idx = first(Δx_Θ_PV_Idx):last(Δx_Θ_PQ_Idx)

#     Δx_v_PV_PQ_Idx = first(Δx_v_PV_Idx):last(Δx_v_PQ_Idx)

#     # ------------------------------------------------------

#     # computing arrays, indices and views
#     # ------------------------------------------------------

#     gen_Q      = zeros(Float64, num_gen)

#     gen_Q_view =  view(gen_Q, 1:num_gen)

#     Sg_view    = @view Sg[:]

#     update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

#     Sd_view    = @view Sd[:]
    
#     Y          = Ybus_p

#     # ------------------------------------------------------
#     # Parameteres
#     # ------------------------------------------------------
    
#     p_ΔP_ΔQ            = ( perm_matrix, Cg, Sd, Y )

#     p_update_gens_Q    = ( gens_bus_num, gen_Qmax, gen_Qmin )

#     p_x_Θ_Δx_Idx       = ( x_Θ_Idx,        x_v_Idx,
#                            x_Θ_PV_PQ_Idx,  x_v_PV_PQ_Idx,
#                            x_Θ_PV_Idx,     x_v_PV_Idx,
#                            x_Θ_PQ_Idx,     x_v_PQ_Idx,
                           
#                            Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx,
#                            Δx_Θ_PV_Idx,    Δx_v_PV_Idx,
#                            Δx_Θ_PQ_Idx,    Δx_v_PQ_Idx )

#     p_Idx  = (PV_Idx_p, PQ_Idx_p)
    
#     p_V_limits  = ( nodes_vmax, nodes_vmin, gen_Vm, p_Idx )
    
#     p_slack = (slack_Vm, slack_Vθ,x_Θ_SK_Idx, x_v_SK_Idx )

#     p  = ( p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx, p_slack )

#     p_SK_Idx  = (x_Θ_SK_Idx, x_v_SK_Idx)  

#      dyn_powerflow_sol_param = (; Ybus_noperm, nodes_name, branches_name, num_nodes, perm_matrix, gen_Q_view, Sg_view, Sd_view, p)
    
#     # ------------------------------------------------------
#     # Initial condition
#     # ------------------------------------------------------

#     x_Θ                = zeros(num_nodes)

#     x_v                = ones(num_nodes)

#     x0                 = [x_Θ...; x_v...]

#     x0[x_Θ_SK_Idx]   .= slack_Vθ

#     x0[ x_v_SK_Idx]  .= slack_Vm
    
#     x0[x_v_PV_Idx]    .= gen_Vm

#     x0_view          = @view x0[:]

#     x_θ_v_red_view   = @view x0[[collect(x_Θ_PV_PQ_Idx); collect(x_v_PQ_Idx)]]

#     dims_view = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PQ_Idx)]

#     _,_,Θ_PV_Θ_PQ_v_PQ_Idx = create_size_offset_Idx(dims_view)    

#     # zeros is not used for a mismatch
    
#     ΔP_ΔQ = ones(sum([length(PV_Idx), length(PQ_Idx), length(PQ_Idx)]))

#     Jac_col_size = sum([length(PV_PQ_Idx), length(PQ_Idx)])
    
#     Jac_row_size = Jac_col_size

#     Jac_Matrix = zeros(Jac_row_size, Jac_col_size)
    
#     # ------------------------------------------------------
#     # Powerflow
#     # ------------------------------------------------------

#     x_0 = [x_θ_v_red_view;]

#     ΔP_ΔQ_0 = similar(x_0)

#     param = (x0_view, Sg_view, Sd_view, slack_Vm, slack_Vθ, Θ_PV_Θ_PQ_v_PQ_Idx, p)    
    
#     x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param)

#     x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param )
       
#     # ------------------------------------------------------
#     # ------------------------------------------------------

#     sol = nlsolve((ΔP_ΔQ, x_θ_v_red_view) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param) , (Jac_Matrix, x_θ_v_red_view) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param ) , x_0)
    
#     # ------------------------------------------------------

#     df = OnceDifferentiable( (ΔP_ΔQ, x_θ_v_red_view) -> x_reduced_ΔP_ΔQ!(ΔP_ΔQ, x_θ_v_red_view, param ), (Jac_Matrix, x_θ_v_red_view) -> x_reduced_Jac_formula!(Jac_Matrix, x_θ_v_red_view, param), x_θ_v_red_view, ΔP_ΔQ)

#     sol = nlsolve(df, x_θ_v_red_view  )

#     # ------------------------------------------------------

#     sol = nlsolve(only_fj!((ΔP_ΔQ, Jac_Matrix, x_θ_v_red_view) -> x_fj!( ΔP_ΔQ, Jac_Matrix, x_θ_v_red_view, ( x, Sg_view, Sd_view, p ) )) , x_0 )
    
#     # ------------------------------------------------------

#     x_θ_v_red_view .= sol.zero

#     V .= x_to_V(x0, x_Θ_Idx, x_v_Idx)

#     V .= perm_matrix * V
    
#     # ------------------------------------------------------
    
#     # ------------------------------------------------------
#     # Results
#     # ------------------------------------------------------
#     Iinj = nodal_current  = Ybus_noperm * V + Diagonal(inv.(conj.(V))) * ( conj.( Sd_view ) )

#     Inet_inj  = net_current  = Ybus_noperm * V 

#     Sbus     = Diagonal(V) * conj.(Ybus_noperm * V)

#     GenSinj  = Sbus + Sd_view

#     return [(vh, θh, ph, qh, ih_r, ih_i) for (vh, θh, ph, qh, ih_r, ih_i) in
#                 zip(abs.(V), angle.(V), real.(Sbus),
#                      imag.(Sbus), real.(Inet_inj),
#                     imag.(Inet_inj) )]
    
# end
                


# #------------------------------------------------------#
# # ------------------------------------------------------
# # Full
# # ------------------------------------------------------
# #------------------------------------------------------#

# # ------------------------------------------------------
# # Full power flow
# # ------------------------------------------------------


# function Jac_formula(V, Y)

#     m11 = [[i==j ? ∂P_i_∂δ_i(i,V,Y) : ∂P_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

#     M11 = reshape( [m11...;], size(Y) )

#     n12 = [[i==j ? Vm_i_∂P_i_∂Vm_i(i,V,Y) : Vm_j_∂P_i_∂Vm_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

#     N12 = reshape( [n12...;], size(Y) )

#     n21 = [[i==j ? ∂Q_i_∂δ_i(i,V,Y) : ∂Q_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

#     N21 = reshape( [n21...;], size(Y) )

#     m22 = [[i==j ? Vm_i_∂Q_i_Vm_i(i,V,Y) : Vm_j_∂Q_i_Vm_j(i,j,V,Y) for j in 1:length(V)] for i in  1:length(V)] 

#     M22 = reshape( [m22...;], size(Y) )

#     return vcat(hcat(M11', N12'), hcat(N21', M22'))
    
# end


# function full_Jac_formula(V, Y, p_Idx)

#    PV_Idx, PQ_Idx = p_Idx

#     PV_PQ_Idx = first(PV_Idx):last(PQ_Idx)
    
#     m11 = [[i==j ? ∂P_i_∂δ_i(i,V,Y) : ∂P_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

#     M11 = reshape( [m11...;], size(Y) )

#     n12 = [[i==j ? Vm_i_∂P_i_∂Vm_i(i,V,Y) : Vm_j_∂P_i_∂Vm_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

#     N12 = reshape( [n12...;], size(Y) )

#     n21 = [[i==j ? ∂Q_i_∂δ_i(i,V,Y) : ∂Q_i_∂δ_j(i,j,V,Y) for j in 1:length(V)] for i in 1:length(V) ] 

#     N21 = reshape( [n21...;], size(Y) )

#     m22 = [[i==j ? Vm_i_∂Q_i_Vm_i(i,V,Y) : Vm_j_∂Q_i_Vm_j(i,j,V,Y) for j in 1:length(V)] for i in  1:length(V)] 

#     M22 = reshape( [m22...;], size(Y) )

#     # return vcat(hcat(M11,  N12'),
#     #             hcat(N21', M22')
#     #             )

#     N12'[PV_PQ_Idx, PV_Idx] = zeros(length(PV_PQ_Idx), length(PV_Idx) )
#     N21'[PV_Idx, PV_PQ_Idx] = zeros(length(PV_Idx),    length(PV_PQ_Idx) )
#     M22'[PV_Idx, PV_Idx]    = zeros(length(PV_Idx),    length(PV_Idx) )
    
    
#     return vcat(hcat(M11'[PV_PQ_Idx, PV_PQ_Idx], N12'[PV_PQ_Idx, PV_PQ_Idx]),
#                 hcat(N21'[PV_PQ_Idx, PV_PQ_Idx], M22'[PV_PQ_Idx, PV_PQ_Idx]))    
    
# end


# function full_ΔP_ΔQ(V, Sg_view, p_ΔP_ΔQ, p_Idx)

#     perm_matrix, Cg, Sd, Y = p_ΔP_ΔQ

#     PV_Idx, PQ_Idx = p_Idx
    
#     PV_PQ_Idx      = first(PV_Idx):last(PQ_Idx)
    
#     mismatch       = Diagonal(V) * conj.(Y * V) - perm_matrix' * (Cg * Sg_view - Sd)

#     ΔP_ΔQ          = Array( [[ real.(mismatch[PV_Idx]),
#                                real.(mismatch[PQ_Idx]),
#                                # zeros(length(PV_Idx)),
#                                imag.(mismatch[PV_Idx]),
#                                imag.(mismatch[PQ_Idx]) ]...; ] )     
# end



# function full_pf_sol(x0,
#                        (gen_Q_view, Sg_view), p;
#                        maxiter=40,
#                        ftol=1000*eps(),
#                        xtol=1000*eps() )

#     p_update_Jac, p_V_limits, p_update_gens_Q, p_ΔP_ΔQ, p_Idx, p_x_Θ_Δx_Idx = p

#     Jac_Q_PV_row_Idx, _ = p_update_Jac

#     nodes_vmax, nodes_vmin, gen_Vm, p_Idx = p_V_limits

#     gens_bus_num, gen_Qmax, gen_Qmin  = p_update_gens_Q

#     perm_matrix, Cg, Sd, Y = p_ΔP_ΔQ 

#     PV_Idx_p, PQ_Idx_p = p_Idx

#     x_Θ_Idx, x_v_Idx, x_Θ_PV_PQ_Idx, x_v_PV_PQ_Idx, x_Θ_PV_Idx, x_v_PV_Idx, x_Θ_PQ_Idx, x_v_PQ_Idx, Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx, Δx_Θ_PV_Idx, Δx_v_PV_Idx, Δx_Θ_PQ_Idx, Δx_v_PQ_Idx = p_x_Θ_Δx_Idx

#     x_k_1              = deepcopy(x0)

#     x_k_1[x_v_PV_Idx] .= gen_Vm

#     x_k                = deepcopy(x0)

#     Δx                 = zeros(length( x_k ))

#     # @show Δx
    
#     x_k[x_v_PV_Idx]   .= gen_Vm

#     V_k_1              = x_to_V(x_k_1, x_Θ_Idx, x_v_Idx)

#     V_k                = x_to_V(x_k, x_Θ_Idx, x_v_Idx)

#     ΔΘVm_dim           = sum([length(x_Θ_PV_PQ_Idx), length(x_v_PV_PQ_Idx)])
    
#     ΔΘVm               = ones( ΔΘVm_dim ) # Inf

#     ΔP_ΔQ_k            = full_ΔP_ΔQ( V_k_1, Sg_view, p_ΔP_ΔQ, p_Idx )

#     local k  = 1
        
#     while (norm(ΔΘVm) > xtol) && (norm(ΔP_ΔQ_k) > ftol)

#         x_k_1[x_v_PV_Idx]  .= gen_Vm

#         V_k_1              .=  x_to_V(x_k_1,  x_Θ_Idx, x_v_Idx)

#         Jac_k               = full_Jac_formula(V_k_1, Y, p_Idx)

#         ΔP_ΔQ_k             = full_ΔP_ΔQ( V_k_1, Sg_view, p_ΔP_ΔQ, p_Idx ) 

#         ΔΘVm               .=  -Jac_k \ ΔP_ΔQ_k

#         Δx[Δx_Θ_PV_PQ_Idx] .= ΔΘVm[Δx_Θ_PV_PQ_Idx]
        
#         #
#         Δx[Δx_v_PV_PQ_Idx] .= ΔΘVm[Δx_v_PV_PQ_Idx]

#         # # if reduced == true

#         # if length(ΔΘVm) < ΔΘVm_dim
            
#         #     Δx[Δx_v_PQ_Idx] .= ΔΘVm[1+length(Δx_Θ_PV_PQ_Idx):end]
            
#         # else
           
#         #     Δx[Δx_v_PQ_Idx]  .= ΔΘVm[Δx_v_PQ_Idx]
            
#         # end

#         x_k[ x_Θ_PV_PQ_Idx ] .= x_k_1[ x_Θ_PV_PQ_Idx ] .+ Δx[ Δx_Θ_PV_PQ_Idx ]

#         #
#         x_k[ x_v_PV_Idx ]    .= x_k_1[ x_v_PV_Idx ] .* (1 .+ Δx[ Δx_v_PV_Idx ])

#         x_k[ x_v_PQ_Idx ]    .= x_k_1[ x_v_PQ_Idx ] .* (1 .+ Δx[ Δx_v_PQ_Idx ])

#         V_k  .=  x_to_V(x_k, x_Θ_Idx, x_v_Idx)

#         # @show V_k[1]

#         gen_Q = imag.((Diagonal(V_k) * conj.(Y * V_k) + perm_matrix' * Sd)[PV_Idx_p])

#         # ------------------------------------------------------
        
#         Q_nolimit_boolean      = map((x) -> no_limit_violation(x...),
#                                      zip(gen_Q, gen_Qmax, gen_Qmin))

#         Q_limit_boolean        = map((x) -> limit_violation(x...),
#                                      zip(gen_Q, gen_Qmax, gen_Qmin))

#         Q_nolimit_violation    = gen_Q[Q_nolimit_boolean]

#         Q_limit_violation      = gen_Q[Q_limit_boolean]

#         if length(Q_limit_violation) != 0

#             gen_Q .= threshold_limits.(gen_Q, gen_Qmax, gen_Qmin)
            
#             update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q)

#             # for (ind, gens_bus_Idx) in enumerate(gens_bus_num)
#             #     real_part = real(Sg_view[gens_bus_Idx])
#             #     Sg_view[gens_bus_Idx] = real_part + im *  gen_Q[ind]
#             # end

#         else

#             update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q)

#             # for (ind, gens_bus_Idx) in enumerate(gens_bus_num)
#             #     real_part = real(Sg_view[gens_bus_Idx])
#             #     Sg_view[gens_bus_Idx] = real_part + im *  gen_Q[ind]
#             # end

#             x_v_PV_restriction_Idx = x_v_PV_Idx[Q_nolimit_boolean]
#             PV_Q_restriction_Idx   = Jac_Q_PV_row_Idx[Q_nolimit_boolean]

#             for (qh_Idx, vh_Idx) in zip(PV_Q_restriction_Idx, x_v_PV_restriction_Idx )
#                 Jac_k[qh_Idx, :     ] .= 0.0
#                 Jac_k[:,      vh_Idx] .= 0.0
#                 Jac_k[qh_Idx, vh_Idx] = 1.0
#             end

#             ΔΘVm  .=  -Jac_k \ ΔP_ΔQ_k

#             Δx[Δx_Θ_PV_PQ_Idx] = ΔΘVm[Δx_Θ_PV_PQ_Idx]

#             #
#             Δx[Δx_v_PV_Idx]  .= ΔΘVm[Δx_v_PV_Idx]
            
#             Δx[Δx_v_PQ_Idx]  .= ΔΘVm[Δx_v_PQ_Idx]

#             x_k[ x_Θ_PV_PQ_Idx ] .= x_k_1[ x_Θ_PV_PQ_Idx ] .+ Δx[ Δx_Θ_PV_PQ_Idx ]

#             #
#             x_k[ x_v_PV_Idx ] .= x_k_1[ x_v_PV_Idx ] .* (1 .+ Δx[ Δx_v_PV_Idx ])
            
#             x_k[ x_v_PQ_Idx ] .= x_k_1[ x_v_PQ_Idx ] .* (1 .+ Δx[ Δx_v_PQ_Idx ])

#             V_k                 .=  x_to_V(x_k, x_Θ_Idx, x_v_Idx)
#         end

#         # ------------------------------------------------------
        
#         update_gens_Q!(Sg_view, gen_Q, p_update_gens_Q)

#         x_k_1 .= x_k

#         k += 1
        
#         if k == maxiter
#             @warn "Maximum number of iterations reached."

#             break
#         end

#     end

#     return k, V_k, ΔP_ΔQ_k, Δx   
# end


# function full_powerflow( case_data_func )

#     Slack_nodes, Gen_nodes, Load_nodes, Nodes, Shunts, Branches = case_data_func()

#     perm_matrix, idx_perm_nodes_idx_and_type,idx_slack_nodes, idx_source_nodes, idx_demand_nodes = get_nodes_permutation(Nodes)

#     num_nodes     = get_network_nodal_size(Branches)

#     gens_S        = get_nodes_S(Gen_nodes)

#     gen_Q         = imag(gens_S) 

#     gens_bus_num  = get_gens_bus_num(Gen_nodes)

#     Cg            = get_Cg(Gen_nodes, num_nodes)

#     Sg            = get_Sg(Gen_nodes, num_nodes)

#     gen_Vm        = get_nodes_Vm(Gen_nodes) 

#     num_gen       =  length(gen_Vm)

#     gen_Qmax      = get_nodes_Qmax(Gen_nodes)

#     gen_Qmin      = get_nodes_Qmin(Gen_nodes)

#     Sd            = get_loads_Sd(Load_nodes)

#     Ybus_noperm   = get_Ybus(Nodes, Branches, shunts=Shunts)

#     slack_Vm      = get_nodes_Vm(Slack_nodes)

#     slack_Vθ      = get_nodes_Vθ(Slack_nodes)

#     nodes_vmax    = get_nodes_vmax(Nodes)

#     nodes_vmax_p  = perm_matrix' * nodes_vmax

#     nodes_vmin    = get_nodes_vmin(Nodes)

#     nodes_vmin_p  = perm_matrix' * nodes_vmin

#     nodes_name    =  get_nodes_Bus_name(Nodes)

#     Ybus_p        = perm_matrix' * (Ybus_noperm) * perm_matrix

#     # PV_Idx        = idx_source_nodes   # gen nodes indices from perm 
#     # PQ_Idx        = idx_demand_nodes   # load nodes indices from perm

#     SK_Idx_p      = idx_slack_nodes

#     PV_Idx_p      = idx_source_nodes   # gen nodes indices from perm 

#     PQ_Idx_p      = idx_demand_nodes   # load nodes indices from perm

#     PV_PQ_Idx_p   = first(PV_Idx_p):last(PQ_Idx_p)   # [PV_Idx;PQ_Idx]

#     # ------------------------------------------------------
#     # Indices
#     # ------------------------------------------------------

#     # Note slack bus is remove, hence PV indices start at 1.

#     SK_PV_PQ_V_dim     = [length(SK_Idx_p),length(PV_Idx_p),length(PQ_Idx_p)]
#     SK_PV_PQ_V_offset  = create_offsets(SK_PV_PQ_V_dim)
#     SK_PV_PQ_V_Idx     = create_idxs(SK_PV_PQ_V_offset, SK_PV_PQ_V_dim)

#     SK_V_Idx           = SK_PV_PQ_V_Idx[1]
#     PV_V_Idx           = SK_PV_PQ_V_Idx[2]
#     PQ_V_Idx           = SK_PV_PQ_V_Idx[3]

#     # num_PV_PQ_nodes = sum([length(PV_Idx), length(PQ_Idx)])

#     x_dim      = [length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx),
#                   length(SK_V_Idx),length(PV_V_Idx),length(PQ_V_Idx)]

#     x_offset   = create_offsets(x_dim)
#     x_Idx      = create_idxs(x_offset, x_dim)

#     x_Θ_SK_Idx = x_Idx[1]
#     x_Θ_PV_Idx = x_Idx[2]
#     x_Θ_PQ_Idx = x_Idx[3]
#     x_v_SK_Idx = x_Idx[4]
#     x_v_PV_Idx = x_Idx[5]
#     x_v_PQ_Idx = x_Idx[6]

#     x_Θ_Idx       = first(x_Θ_SK_Idx):last(x_Θ_PQ_Idx)
#     x_v_Idx       = first(x_v_SK_Idx):last(x_v_PQ_Idx)

#     x_Θ_PV_PQ_Idx = first(x_Θ_PV_Idx):last(x_Θ_PQ_Idx)
#     x_v_PV_PQ_Idx = first(x_v_PV_Idx):last(x_v_PQ_Idx)

#     Δx_dim      = [length(x_Θ_PV_Idx), length(x_Θ_PQ_Idx), length(x_v_PV_Idx),length(x_v_PQ_Idx)]

#     Δx_offset   = create_offsets( Δx_dim )
#     Δx_Idx      = create_idxs( Δx_offset, Δx_dim )

#     Δx_Θ_PV_Idx = Δx_Idx[1]
#     Δx_Θ_PQ_Idx = Δx_Idx[2]

#     Δx_v_PV_Idx = Δx_Idx[3]
#     Δx_v_PQ_Idx = Δx_Idx[4]

#     Δx_Θ_PV_PQ_Idx = first(Δx_Θ_PV_Idx):last(Δx_Θ_PQ_Idx)

#     Δx_v_PV_PQ_Idx = first(Δx_v_PV_Idx):last(Δx_v_PQ_Idx)


#     Jac_dim      = [length(PV_Idx_p), length(PQ_Idx_p),length(PV_Idx_p), length(PQ_Idx_p)]
#     Jac_offset   = create_offsets(Jac_dim )
#     Jac_Idx      = create_idxs(Jac_offset, Jac_dim)

#     Jac_P_PV_row_Idx  = Jac_Idx[1]
#     Jac_P_PQ_row_Idx  = Jac_Idx[2]
#     Jac_Q_PV_row_Idx  = Jac_Idx[3]
#     Jac_Q_PQ_row_Idx  = Jac_Idx[4]

    
#     # ------------------------------------------------------

#     # computing arrays, indices and views
#     # ------------------------------------------------------

#     gen_Q      = zeros(Float64, num_gen)

#     gen_Q_view =  view(gen_Q, 1:num_gen)

#     Sg_view    = @view Sg[:]

#     update_gen_reactive_power!(Sg_view, gens_bus_num, gen_Q_view)

#     Y          = Ybus_p

#     # ------------------------------------------------------
#     # Parameteres
#     # ------------------------------------------------------

#     p_Idx              = ( PV_Idx_p, PQ_Idx_p )

#     p_ΔP_ΔQ            = ( perm_matrix, Cg, Sd, Y )

#     p_update_gens_Q    = ( gens_bus_num, gen_Qmax, gen_Qmin )

#     p_x_Θ_Δx_Idx       = ( x_Θ_Idx,        x_v_Idx,
#                            x_Θ_PV_PQ_Idx,  x_v_PV_PQ_Idx,
#                            x_Θ_PV_Idx,     x_v_PV_Idx,
#                            x_Θ_PQ_Idx,     x_v_PQ_Idx,
#                            Δx_Θ_PV_PQ_Idx, Δx_v_PV_PQ_Idx,
#                            Δx_Θ_PV_Idx,    Δx_v_PV_Idx,
#                            Δx_Θ_PQ_Idx,    Δx_v_PQ_Idx )

#     p_Idx              = (PV_Idx_p, PQ_Idx_p)
    
#     p_V_limits         = ( nodes_vmax, nodes_vmin, gen_Vm, p_Idx )

#     p_update_Jac       = ( Jac_Q_PV_row_Idx, x_v_PV_Idx )

#     p_slack = (slack_Vm, slack_Vθ,x_Θ_SK_Idx, x_v_SK_Idx)

#     p                  = ( p_update_Jac, p_V_limits,
#                           p_update_gens_Q, p_ΔP_ΔQ,
#                           p_Idx, p_x_Θ_Δx_Idx )


#     # ------------------------------------------------------
#     # Initial condition
#     # ------------------------------------------------------

#     x_Θ                = zeros(num_nodes)

#     x_v                = ones(num_nodes)

#     x0                 = [x_Θ...; x_v...]

#     # x_k_1              = deepcopy(x0)
#     # x_k_1[x_v_PV_Idx] .= gen_Vm
#     # x_k                = deepcopy(x0)
#     # Δx                 = zeros(length( x_k ))
#     # x_k[x_v_PV_Idx]   .= gen_Vm
#     # V_k                = x_to_V(x_k, x_Θ_Idx, x_v_Idx)

#     k, V_k, ΔP_ΔQ_k, Δx = full_pf_sol(x0, (gen_Q_view, Sg_view), p;
#                                       maxiter=40,
#                                       ftol=1000*eps(),
#                                       xtol=1000*eps() )


#     V         = perm_matrix * V_k

#     Sbus      = Diagonal(V) * conj.(Ybus_noperm * V)

#     GenSinj   = Diagonal(V) * conj.(Ybus_noperm * V) + Sd

#     If         = get_Yf(Nodes, Branches) * V
    
#     It         = get_Yt(Nodes, Branches) * V

#     Ibranches  = If + It

#     Iinj = nodal_current_injection(V, p_ΔP_ΔQ)

#     # branches_name  = get_branches_name(Branches)    
#     # nodes_name    =  get_nodes_Bus_name(Nodes)

#     branches_name  = collect(keys(Branches))    
#     nodes_name     = collect(keys(Nodes))   

#     branch_dict_init = OrderedDict(name =>
#         (real(i_f), imag(i_f), real(i_t), imag(i_t)) for (name, i_f, i_t) in zip(branches_name, If, It))

#     bus_dict_init = OrderedDict(name =>
#         (vh, θh, ph, qh) for (name,
#                               vh,
#                               θh,
#                               ph,
#                               qh) in zip(nodes_name,
#                                          round.( abs.(V); digits=4),
#                                          round.( angle.(V); digits=4),
#                                          round.( real.(Sbus); digits=4),
#                                          round.( imag.(Sbus); digits=4) ))      
    
    
#     return Dict("Iinj" => Iinj, "Sg" => Sg, "Sd" => Sd, "Vm" => abs.(V), "Vθ" => angle.(V), "Vbus" => V, "Ibranches" => Ibranches, "Ybus" => Ybus_noperm , "Sbus" => Sbus, "GenSinj" => GenSinj, "dict_init" =>Dict("bus_dict_init" => bus_dict_init,  "branch_dict_init" =>  branch_dict_init))
    
#     # return perm_matrix, k, V_k, ΔP_ΔQ_k, Δx 
# end

