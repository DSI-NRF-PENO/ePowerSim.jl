# # (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# # Licensed under GNU GPL v3 (see LICENSE file)


# # AA Yusuff : yusufaa@unisa.ac.za

#--------------------------------------------------
#------------------------------------------------

function load_K(
    load_dfP, load_dfQ,
    vh, param_kP, param_kQ)

    dfP_dvh = load_dfP(vh, param_kP )
    dfQ_dvh = load_dfQ(vh, param_kQ )

    K_load = (dfP_dvh - im * dfQ_dvh) / (
        (dfP_dvh)^2 + (dfQ_dvh)^2) 
    
    return K_load
    
end


function dyn_load_fP(vh, param_kP )

    PL0,  kP1, kP2 = param_kP

    return PL0 + kP1 * vh + kP2 * (vh)^2
end


function dyn_load_fQ(vh, param_kQ )

    QL0, kQ1, kQ2 = param_kQ

    return QL0 + kQ1 * vh + kQ2 * (vh)^2
end


function dyn_load_dfP(vh, param_kP )

    PL0,  kP1, kP2 = param_kP

    return kP1 + 2 * kP2 * vh
end


function dyn_load_dfQ(vh, param_kQ )

    QL0, kQ1, kQ2 = param_kQ

    return kQ1 + 2 * kQ2 * vh
end


function dyn_load_model(vh, param_kP, param_kQ)

    return dyn_load_fP(vh, param_kP) +
        im * dyn_load_fQ(vh, param_kQ)
    
end


function dyn_load_K(dyn_load_dfP, dyn_load_dfQ,
                    vh, param_kP, param_kQ)
 
    return load_K(dyn_load_dfP, dyn_load_dfQ,
                  vh, param_kP, param_kQ)
    
end


function dyn_load_K(vh, param_kP, param_kQ)

    dfP_dvh = dyn_load_dfP(vh, param_kP )
    dfQ_dvh = dyn_load_dfQ(vh, param_kQ )

    K_load = (dfP_dvh - im * dfQ_dvh) / (
        (dfP_dvh)^2 + (dfQ_dvh)^2) 
    
    return K_load
    
end


# ------- constant power


function const_power_load_fP(vh, param_kP )

    PL0, _,_ = param_kP

    return PL0
end

function const_power_load_fQ(vh, param_kQ )

    QL0 ,_, _ = param_kQ

    return QL0
end


function const_power_load_dfP(vh, param_kP )

    PL0, _, _  = param_kP

    return 0.0
end

function const_power_load_dfQ(vh, param_kQ )

    QL0, _, _  = param_kQ

    return 0.0
end


function const_power_load_model(vh, param_kP, param_kQ)

    return const_power_load_fP(vh, param_kP) +
        im * const_power_load_fQ(vh, param_kQ)
    
end


function const_power_load_K(
    const_power_load_dfP,
    const_power_load_dfQ,
    vh, param_kP, param_kQ)
 
    
    return load_K(const_power_load_dfP,
                  const_power_load_dfQ,
                  vh, param_kP, param_kQ)
    
end


function const_power_load_K(vh, param_kP, param_kQ)

    dfP_dvh = const_power_load_dfP(vh, param_kP )
    dfQ_dvh = const_power_load_dfQ(vh, param_kQ )

    K_load = (dfP_dvh - im * dfQ_dvh) / (
        (dfP_dvh)^2 + (dfQ_dvh)^2) 
    
    return K_load
    
end


# ------- constant current



function const_current_load_fP(vh, param_kP )

    _, kP1, _ = param_kP

    return kP1 * vh
end

function const_current_load_fQ(vh, param_kQ )

    _, kQ1, _  = param_kQ

    return kQ1 * vh 
end


function const_current_load_dfP(vh, param_kP )

    _, kP1, _  = param_kP

    return kP1 
end

function const_current_load_dfQ(vh, param_kQ )

    _, kQ1, _ = param_kQ

    return kQ1
end


function const_current_load_model(vh, param_kP, param_kQ)

    return const_current_load_fP(vh, param_kP) +
        im * const_current_load_fQ(vh, param_kQ)
    
end


function const_current_load_K(
    const_current_load_dfP, const_current_load_dfQ,
    vh, param_kP, param_kQ)
 
    
    return load_K(
        const_current_load_dfP, const_current_load_dfQ,
        vh, param_kP, param_kQ)
    
end



function const_current_load_K(vh, param_kP, param_kQ)

    dfP_dvh = const_current_load_dfP(vh, param_kP )
    dfQ_dvh = const_current_load_dfQ(vh, param_kQ )

    K_load = (dfP_dvh - im * dfQ_dvh) / (
        (dfP_dvh)^2 + (dfQ_dvh)^2) 
    
    return K_load
    
end


# ------- constant impedance


function const_impedance_load_fP(vh, param_kP )

    _, _, kP2 = param_kP

    return kP2 * (vh)^2
end

function const_impedance_load_fQ(vh, param_kQ )

    _, _, kQ2 = param_kQ

    return kQ2 * (vh)^2
end


function const_impedance_load_dfP(vh, param_kP )

    _, _, kP2 = param_kP

    return 2 * kP2 * vh
end

function const_impedance_load_dfQ(vh, param_kQ )

    _, _, kQ2 = param_kQ

    return 2 * kQ2 * vh
end


function const_impedance_load_model(
    vh, param_kP, param_kQ)

    return const_impedance_load_fP(vh, param_kP) +
        im * const_impedance_load_fQ(vh, param_kQ)
    
end


function const_impedance_load_K(
    const_impedance_load_dfP, const_impedance_load_dfQ,
    vh, param_kP, param_kQ)
 
    return load_K(
        const_impedance_load_dfP, const_impedance_load_dfQ,
        vh, param_kP, param_kQ)
    
end


function const_impedance_load_K(vh, param_kP, param_kQ)

    dfP_dvh = const_impedance_load_dfP(vh, param_kP )
    dfQ_dvh = const_impedance_load_dfQ(vh, param_kQ )

    K_load = (dfP_dvh - im * dfQ_dvh) / (
        (dfP_dvh)^2 + (dfQ_dvh)^2) 
    
    return K_load
    
end


########################################################
# ------------------------------------------------------
#  Static Nodes Components Structures
# ------------------------------------------------------
########################################################


"""
    Slack!(dx, x, p_agg, t)


Generic nodal struct for power flow.
"""

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
    Bus_num::Int64 = parse(Int, split(
        lowercase(Bus),"bus")[2] )
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
    state_vars_syms::Vector{Symbol} =
        Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} =
        Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} =
        DAE_MassMatrix(length(state_vars),
                       length(algebraic_vars))
    param::Vector{Symbol} =
        Symbol[:kV, :P, :Q, :Vm, :Vθ, :vmax,
               :vmim, :Bus_type]
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
    Bus_num::Int64 = parse(Int, split(
        lowercase(Bus),"bus")[2]) 
    y_shunt::ComplexF64 = 0.0 + im * 0.0  
    Vm::Float64      = 1.0
    Vθ::Float64      = 0.0
    Bus_type::Symbol = :ShuntElement
    algebraic_vars_syms::Vector{Symbol} = Symbol[]
    algebraic_vars::Vector{Symbol} = Symbol[]
    state_vars::Vector{Symbol} = Symbol[]
    state_vars_syms::Vector{Symbol} =
        Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} =
        Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} =
        DAE_MassMatrix(length(state_vars),
                       length(algebraic_vars))
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
    Bus_num::Int64 = parse(Int, split(
        lowercase(Bus),"bus")[2] )     
    P::Float64        = 0.0 
    Q::Float64        = 0.0 
    S::ComplexF64     = P + im * Q
    Vm::Float64       = 1.0 
    Vθ::Float64       = 0.0 
    Bus_type::Symbol  = :EnergyStorage

    algebraic_vars_syms::Vector{Symbol} = Symbol[]
    algebraic_vars::Vector{Symbol} = Symbol[]
    state_vars::Vector{Symbol} = Symbol[]
    state_vars_syms::Vector{Symbol} =
        Symbol[state_vars...;algebraic_vars...]
    syms::Vector{Symbol} =
        Symbol[algebraic_vars_syms...;state_vars_syms...]
    dim::Int64 = length(state_vars_syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} =
        DAE_MassMatrix(length(state_vars),
                       length(algebraic_vars))
    param::Vector{Symbol} =
        Symbol[:P, :Q, :Vm,
               :Vθ, :Bus_type]
    func::Vector{Function} =
        Function[EnergyStorage!]
end


########################################################
# ------------------------------------------------------
#  Dynamic Nodes Components Structures
# ------------------------------------------------------
########################################################


"""
Generalisation of network and load dynamics

Sauer: page 130

Constant impedance load: Sauer, page 178 :190
eq 7.81

This is a mutating function vertexfunction! for nodes of type PQAlgebraic. The constructor StaticVertex, ODEVertex or DDEVertex are called in order to turn vertexfunction! into a VertexFunction object compatible with network_dynamics.


Milano: page 257 261, see figure 10.1

"""

function hybrid_pf_PQ_Const_P!(dx, x, p_agg, t)
    
    # cb_sw, src_i, dst_i, f_t, q, node_idx_and_inc_edges, node_inc_edges_Ybr_orient, nodes_u_view =  p_agg

    (cb_sw, src_i, dst_i, f_t, q,
     global_pf_param, node_pf_param) =  p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr,
     node_inc_edges_orient, nodes_u_view)  =
         node_pf_param

    p      = q[1]
    
    P      = p[1]
    Q      = p[2]
    kPL    = p[3]
    kQL    = p[4]
    Y_n    = p[5]
    
    bus_no = p[6]

    #----------------------------

    u_r,  u_i = x
    
    u = u_r + u_i * im     
    
    # #----------------------------        
    # # node_pf
    # #----------------------------
    # # Do not delete
    # #----------------------------    
    
    # U_r   = real(Uk)
    
    # U_i   = imag(Uk)

    #----------------------------        
    # Network current
    #----------------------------    
    
    i = dynamic_nodal_current_balance(src_i, dst_i)  +
        im * u * Y_n
        
    i_mag = abs(i)

    S = P + im * Q

    U = S * i / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)

    #----------------------------

    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i
    
    # #---------------------------
    
    return nothing
end


function node_pf_PQ_Const_P!(dx, x, p_agg, t)
    
    # cb_sw, src_i, dst_i, f_t, q, node_idx_and_inc_edges, node_inc_edges_Ybr_orient, nodes_u_view =  p_agg

    (cb_sw, src_i, dst_i,
     f_t, q, node_pf_param) =  p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr, node_inc_edges_orient,
     nodes_u_view)  = node_pf_param

    p      = q[1]
    
    P      = p[1]
    Q      = p[2]
    kPL    = p[3]
    kQL    = p[4]
    Y_n    = p[5]
    
    bus_no = p[6]

    S = P + im * Q

    #----------------------------

    u_r,  u_i = x
    
    u = u_r + u_i * im     
    
    #----------------------------

    
    my_node_idx =
        node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j =
        node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr

    # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length( edges_Ybr )]
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -
        sum([ykj * vj for (ykj, vj) in zip(Ykj, Uj) ]))

    # calculating Uk twice
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -
        sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]))

    # ------------------------------------------------------
    
    U_r   = real(Uk)
    
    U_i   = imag(Uk)

    #----------------------------        
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i
    
    # #---------------------------
    
    return nothing
end

function global_pf_PQ_Const_P!(dx, x, p_agg, t)
    
    # cb_sw, src_i, dst_i, f_t, q, pf_U, Inet =  p_agg

    (cb_sw, src_i, dst_i,
     f_t, q, global_pf_param) =  p_agg

    pf_U, Inet = global_pf_param

    p      = q[1]
    
    P      = p[1]
    Q      = p[2]
    kPL    = p[3]
    kQL    = p[4]
    Y_n    = p[5]
    
    bus_no = p[6]

    S = P + im * Q
    
    #----------------------------

    u_r = x[1]
    
    u_i = x[2]
    
    u = u_r + u_i * im     
    
    #----------------------------

    i = Inet  # + im * u * Y_n
        
    i_mag = abs(i)
    
    # i_ang = angle(i)

    S = P + im * Q

    # U = S / cong(i)

    U = S * i / (i_mag)^2
    
    dx[1] = real(U) - u_r
    
    dx[2] = imag(U) - u_i
    
    return nothing
end



function network_current_PQ_Const_P!(dx, x, p_agg, t)
    
    (cb_sw, src_i, dst_i,
     f_t, q) =  p_agg

    p      = q[1]
    
    P      = p[1]
    Q      = p[2]
    kPL    = p[3]
    kQL    = p[4]
    Y_n    = p[5]
    
    bus_no = p[6]

    #----------------------------

    u_r, u_i = x
    
    u = u_r + u_i * im

    #----------------------------    
    
    i = dynamic_nodal_current_balance(src_i, dst_i) +
        im * u * Y_n
        
    i_mag = abs(i)

    S = P + im * Q

    U = S * i / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)

    #----------------------------
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i

    #---------------------------
    
    return nothing
end


function initial_pf_PQ_Const_P!(dx, x, p_agg, t)
    
    (cb_sw, src_i, dst_i,
     f_t, q) =  p_agg

    p      = q[1]
    
    P      = p[1]
    Q      = p[2]
    kPL    = p[3]
    kQL    = p[4]
    Y_n    = p[5]
    
    bus_no = p[6]

    #----------------------------

    u_r, u_i = x
    
    u = u_r + u_i * im

    #----------------------------

    vh, θh, ph, qh, i_r, i_i = f_t[1]
    
    i = i_r + im * i_i + im * u * Y_n
        
    i_mag = abs(i)
    
    # i_ang = angle(i)

    S = P + im * Q

    # U = S / cong(i)
    # U = (S / cong(i)) *  (i /i)
    # U = S * i / (i_mag)^2

    U = S * i / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)

    #----------------------------
    
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i

    #---------------------------
    
    return nothing
end


#####

"""
This structure is used to organise data for nodes of type PQAlgebraic.
    """

@kwdef struct PQ_Const_P  <: SdNonGen
    Bus::String = "bus0"
    name::String = lowercase(Bus)     
    P::Float64 = 0.0
    Q::Float64 = 0.0
    kPL::Tuple{Float64, Float64, Float64} = (1.0, 0.0, 0.0) 
    kQL::Tuple{Float64, Float64, Float64} = (1.0, 0.0, 0.0)
    
    Y_n::Float64 = 0.0

    vmax::Float64    = 1.06
    vmin::Float64    = 0.97    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :PQ_Const_P

    Bus_num::Int64 = parse(Int, split(
        lowercase(Bus),"bus")[2] )
    
    Bus_type::Symbol  = :Load    
    
    state_vars_syms::Vector{Symbol} = Symbol[  ]
    algebraic_vars_syms::Vector{Symbol} =
        Symbol[ :u_r, :u_i ]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[
        state_vars_syms...;algebraic_vars_syms...]
    
    dict_state_syms::OrderedDict{Symbol, Int64} =
        OrderedDict(sym => ind for (ind, sym) in
                        enumerate(syms))
    
    ur_ui_idx::Vector{Int64} =
        Int64[dict_state_syms[:u_r], dict_state_syms[:u_i]]
    
    dim::Int64 = length(syms)
    
    mass_matrix::Diagonal{Int64, Vector{Int64}} =
        DAE_MassMatrix(
            length(state_vars_syms),
            length(algebraic_vars_syms))
    
    dae_var::Vector{Bool} = DAE_BoolVector(
        length(state_vars_syms),
        length(algebraic_vars_syms) )
    
    param::Vector{Symbol} =
        Symbol[:P, :Q, :kPL, :kQL,
               :Y_n, :Bus, :vmax, :vmin]
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} =
        OrderedDict(sym => ind
                    for (ind, sym) in
                        enumerate( param ) )

    P_Q_idx::Vector{Int64} =
        Int64[ dict_param_syms_Idx[ :P ],
               dict_param_syms_Idx[ :Q ]]  
    
    param_values::Vector{
        Union{
            Int64,Float64,Tuple{Float64,Float64,Float64}}} =
                Union{
                    Int64, Float64,
                    Tuple{Float64, Float64, Float64}}[
                        P, Q, kPL, kQL, Y_n,
                        parse(Int, split(lowercase(Bus),
                                         "bus")[2]) ,
                        vmax, vmin ]

    func::Vector{Function} =
        Function[ initial_pf_PQ_Const_P!,
                  network_current_PQ_Const_P!,
                  global_pf_PQ_Const_P!,
                  node_pf_PQ_Const_P!,
                  hybrid_pf_PQ_Const_P!]

    control_sig_syms::Vector{Symbol} = Symbol[ ]
    
    control_sig::Vector{Float64} =
        ones(length(control_sig_syms))
    
    output_sig_syms::Vector{Symbol} = Symbol[]
    
    output_sig::Vector{Float64} =
        ones(length(output_sig_syms))
    
    cb_state_event_func::Vector{Function} = Function[]
    cb_state_affect_func::Vector{Function} = Function[]    
    cb_state_syms::Vector{Symbol} = Symbol[]           
    cb_state_conditions::Vector{Float64} = Float64[]
    cb_state_values::Vector{Function} = Function[]    
    # cb_state_values::Vector{Float64}           = Float64[]
    cb_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64 = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function} = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol} = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    # cb_dyn_state_values::Vector{Float64}       = Float64[]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)    


end

function hybrid_pf_Const_I!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i, f_t, q,
     global_pf_param,
     node_pf_param) =  p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr, node_inc_edges_orient,
     nodes_u_view)  = node_pf_param
    
    p      = q[1]
    
    P      = p[1]
    Q      = p[2]
    kPL    = p[3]
    kQL    = p[4]
    Y_n    = p[5]
    
    bus_no = p[6]

    S = P + im * Q

     # ---------------------------------------------------- 

    u_r,  u_i = x
    
    u = u_r + u_i * im     
        
    i = dynamic_nodal_current_balance(src_i, dst_i) +
        im * u * Y_n
    
    i_mag = abs(i)

    i_ang = angle(i)
    
    # ---------------------------------------------------- 
    S = P + im * Q
    
    U = S * i / (i_mag)^2
    
    # U = S * conj(i) / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)

    # ----------------------------------------------------
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i
        
    return nothing    

    
end



function node_pf_Const_I!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i,
     f_t, q, node_pf_param) =  p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr,
     node_inc_edges_orient,
     nodes_u_view)  = node_pf_param
    
    p      = q[1]
    
    P      = p[1]
    Q      = p[2]
    kPL    = p[3]
    kQL    = p[4]
    Y_n    = p[5]
    
    bus_no = p[6]

    S = P + im * Q

    # #----------------------------

    u_r,  u_i = x
    
    u = u_r + u_i * im     
    
    #----------------------------
    
    my_node_idx =
        node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j =
        node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -
        sum([ykj * vj for (ykj, vj) in zip(Ykj, Uj) ]))

    # calculating Uk twice
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -
        sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]))

    # ------------------------------------------------------
    
    
    U_r   = real(Uk)
    
    U_i   = imag(Uk)

    #----------------------------        
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i
    
    # #---------------------------
    
    return nothing    

end



function global_pf_Const_I!(dx, x, p_agg, t)
    
    (cb_sw, src_i, dst_i,
     f_t, q, global_pf_param) =  p_agg

    pf_U, Inet = global_pf_param

    p      = q[1]
    
    P      = p[1]
    Q      = p[2]
    kPL    = p[3]
    kQL    = p[4]
    Y_n    = p[5]
    
    bus_no = p[6]

    S = P + im * Q
    
    #----------------------------

    u_r = x[1]
    
    u_i = x[2]
    
    u = u_r + u_i * im     
    
    #----------------------------

    i = Inet 
        
    i_mag = abs(i)
    
    # i_ang = angle(i)

    S = P + im * Q

    # U = S / cong(i)

    U = S * i / (i_mag)^2
    
    dx[1] = real(U) - u_r
    
    dx[2] = imag(U) - u_i

    return nothing
    
end



function network_current_PQ_Const_I!(dx, x, p_agg, t)
    
    (cb_sw, src_i, dst_i,
     f_t, q) =  p_agg

    p = q[1]
    
    P   = p[1]
    Q   = p[2]
    kPL = p[3]
    kQL = p[4]
    Y_n = p[5]

    # ----------------------------------------------------

    u_r, u_i = x
    
    u = u_r + u_i * im    
    
    # ----------------------------------------------------

    i = dynamic_nodal_current_balance(src_i, dst_i) +
        im * u * Y_n
    
    i_mag = abs(i)

    i_ang = angle(i)
    
    # ---------------------------------------------------- 
    S = P + im * Q
    
    U = S * i / (i_mag)^2
    
    # U = S * conj(i) / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)

    # ----------------------------------------------------
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i

    
    return nothing
end



function initial_pf_PQ_Const_I!(dx, x, p_agg, t)
    
    (cb_sw, src_i, dst_i,
     f_t, q) =  p_agg

    p = q[1]
    
    P   = p[1]
    Q   = p[2]
    kPL = p[3]
    kQL = p[4]
    Y_n = p[5]

    # -------------------------------------    

    u_r, u_i = x
    
    u = u_r + u_i * im
    
    # -------------------------------------    
       
    vh, θh, ph, qh, i_r, i_i = f_t[1]
    
    i = i_r + im * i_i +  im * u * Y_n

    i_mag = abs(i)
    
    i_ang = angle(i)

    # -------------------------------------    

    S = P + im * Q
    
    U = S * i / (i_mag)^2
    
    # U = S * conj(i) / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)

    # -------------------------------------    
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i

    #---------------------------
    
    return nothing
end

#####

"""

It is used to organise data for nodes of type PQAlgebraic.
    """

@kwdef struct PQ_Const_I  <: SdNonGen
    Bus::String  = "bus0"
    name::String = lowercase(Bus)     
    P::Float64 = 0.0
    Q::Float64 = 0.0
    kPL::Tuple{Float64, Float64, Float64} = (0.0, 1.0, 0.0) 
    kQL::Tuple{Float64, Float64, Float64} = (0.0, 1.0, 0.0)

    Y_n::Float64 = 0.0
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.97    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :PQ_Const_I

    Bus_num::Int64 =
        parse(Int, split(lowercase(Bus),"bus")[2] )
    
    Bus_type::Symbol  = :Load    
        
    state_vars_syms::Vector{Symbol} = Symbol[  ]
    
    algebraic_vars_syms::Vector{Symbol} =
        Symbol[ :u_r, :u_i ]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} =
        Symbol[ state_vars_syms...;algebraic_vars_syms...]
    
    dict_state_syms::OrderedDict{Symbol, Int64} =
        OrderedDict(sym => ind
                    for (ind, sym) in enumerate(syms))
    
    ur_ui_idx::Vector{Int64} =
        Int64[dict_state_syms[:u_r],
              dict_state_syms[:u_i]]
    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} =
        DAE_MassMatrix(length(state_vars_syms),
                       length(algebraic_vars_syms))
    
    dae_var::Vector{Bool} =
        DAE_BoolVector(length(state_vars_syms),
                       length(algebraic_vars_syms) )
    
    param::Vector{Symbol} =
        Symbol[:P, :Q, :kPL, :kQL,
               :Y_n, :vmax, :vmin ]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} =
        OrderedDict(sym => ind
                    for (ind, sym) in
                        enumerate( param ) )

    P_Q_idx::Vector{Int64} =
        Int64[ dict_param_syms_Idx[ :P ],
               dict_param_syms_Idx[ :Q ]]  
    
    param_values::Vector{Union{Float64, Tuple{Float64, Float64, Float64}}} = Union{Float64, Tuple{Float64, Float64, Float64}}[P, Q, kPL, kQL, Y_n, vmax, vmin ]

    func::Vector{Function} = Function[initial_pf_PQ_Const_I!, network_current_PQ_Const_I!, global_pf_Const_I!, node_pf_Const_I!, hybrid_pf_Const_I! ]

    control_sig_syms::Vector{Symbol} = Symbol[ ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
    
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]
    cb_state_values::Vector{Function}          = Function[]    
    # cb_state_values::Vector{Float64}           = Float64[]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    # cb_dyn_state_values::Vector{Float64}       = Float64[]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

end



function new_PQ_Const_Z!(dx, x, p_agg, t)
    

    (cb_sw, src_i, dst_i,
     f_t, q, node_pf_param) =  p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr,
     node_inc_edges_orient,
     nodes_u_view)  = node_pf_param

    p = q[1]
    
    P   = p[1]
    Q   = p[2]
    kPL = p[3]
    kQL = p[4]

    Y_n = p[5]

    #----------------------------
    
    u_r = x[1]
    
    u_i = x[2]
    
    u = u_r + u_i * im     

    i = dynamic_nodal_current_balance(src_i, dst_i) +
        u * im * Y_n 


    i_mag = abs(i)
    
    i_ang = angle(i)    


    S = P + im * Q

    # U = S * i / (i_mag)^2

    U = S * conj(i) / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i

    #---------------------------    
    
    return nothing
end


function PQ_Const_Z!(dx, x, p_agg, t)
    
    # cb_sw, f_t, p = p_agg

    cb_sw, src_i, dst_i, f_t, q =  p_agg

    p = q[1]
    
    P   = p[1]
    Q   = p[2]
    kPL = p[3]
    kQL = p[4]

    Y_n = p[5]
    
    #----------------------------
    
    u_r = x[1]
    
    u_i = x[2]
    
    u = u_r + u_i * im     

    i = dynamic_nodal_current_balance(src_i, dst_i) +
        u * im * Y_n 


    i_mag = abs(i)
    
    i_ang = angle(i)    


    S = P + im * Q

    # U = S * i / (i_mag)^2

    U = S * conj(i) / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i

    #---------------------------    
    
    return nothing
end

#####

"""
It is used to organise data for nodes of type PQAlgebraic.
    """

@kwdef struct PQ_Const_Z  <: SdNonGen
    Bus::String = "bus0"
    P::Float64 = 0.0
    Q::Float64 = 0.0   
    name::String  = lowercase(Bus) 
    kPL::Tuple{Float64, Float64, Float64} = (0.0, 0.0, 1.0) 
    kQL::Tuple{Float64, Float64, Float64} = (0.0, 0.0, 1.0)

    Y_n::Float64 = 0.0
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.97    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :PQ_Const_Z

    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )     
    Bus_type::Symbol  = :Load    
        
    state_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i  ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ ]
        
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    
    param::Vector{Symbol} = Symbol[:P, :Q, :kPL, :kQL, :Y_n, :vmax, :vmin]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]  
    
    param_values::Vector{Union{Float64, Tuple{Float64, Float64, Float64}}} = Union{Float64, Tuple{Float64, Float64, Float64}}[P, Q, kPL, kQL, Y_n, vmax, vmin]
    func::Vector{Function} = Function[PQ_Const_Z!, PQ_Const_Z!, new_PQ_Const_Z!]
    control_sig_syms::Vector{Symbol} = Symbol[ ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
    
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]
    cb_state_values::Vector{Function}          = Function[]    
    # cb_state_values::Vector{Float64}           = Float64[]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    # cb_dyn_state_values::Vector{Float64}       = Float64[]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

end


function new_PQ_dyn_load!(dx, x, p_agg, t)
    
    # cb_sw, f_t, p = p_agg

    (cb_sw, src_i, dst_i, f_t, p,
     node_idx_and_inc_edges,
     node_inc_edges_Ybr_orient,
     nodes_u_view) =  p_agg
    
    P   = p[1]
    Q   = p[2]

    kPL = p[3]
    kQL = p[4]


    Y_n = p[5]
    #----------------------------
    #----------------------------

    u_r = x[1]
    
    u_i = x[2]
    
    u = u_r + u_i * im    
    
    i = dynamic_nodal_current_balance(src_i, dst_i) +
        u * im * Y_n  


    i_mag = abs(i)
    
    i_ang = angle(i)    


    S = P + im * Q

    # U = S * i / (i_mag)^2
    
    U = S * conj(i) / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i

    #---------------------------      
    
    return nothing
end

function PQ_dyn_load!(dx, x, p_agg, t)
    
    # cb_sw, f_t, p = p_agg

    cb_sw, src_i, dst_i, f_t, p =  p_agg
    
    P   = p[1]
    Q   = p[2]

    kPL = p[3]
    kQL = p[4]


    Y_n = p[5]
    #----------------------------

    u_r = x[1]
    
    u_i = x[2]
    
    u = u_r + u_i * im    
    
    i = dynamic_nodal_current_balance(src_i, dst_i) +
        u * im * Y_n  


    i_mag = abs(i)
    
    i_ang = angle(i)    


    S = P + im * Q

    # U = S * i / (i_mag)^2
    
    U = S * conj(i) / (i_mag)^2
    
    U_r = real(U)
    
    U_i = imag(U)
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i

    #---------------------------      
    
    return nothing
end

#####

"""

It is used to organise data for nodes of type PQAlgebraic.

    """

@kwdef struct PQ_dyn_load  <: SdNonGen

    Bus::String  = "bus0"
    name::String = lowercase(Bus)     
    P::Float64 = 0.0
    Q::Float64 = 0.0
    Y_n::Float64 = 0.0
    kPL::Tuple{Float64, Float64, Float64} = (0.3, 0.3, 0.4) 
    kQL::Tuple{Float64, Float64, Float64} = (0.3, 0.3, 0.4)

    comp_type::Symbol = :PQ_dyn_load
    
    Bus_num::Int64 =
        parse(Int, split(lowercase(Bus),"bus")[2] )     
    Bus_type::Symbol  = :Load    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.97    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0
    
    state_vars_syms::Vector{Symbol} = Symbol[  ]
    algebraic_vars_syms::Vector{Symbol} =
        Symbol[ :u_r, :u_i  ]
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    
    param::Vector{Symbol} = Symbol[:P, :Q, :kPL, :kQL, :Y_n, :vmax, :vmin]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]  
    
    param_values::Vector{Union{Float64, Tuple{Float64, Float64, Float64}}} = Union{Float64, Tuple{Float64, Float64, Float64}}[P, Q, kPL, kQL, Y_n, vmax, vmin ]
    func::Vector{Function} = Function[PQ_dyn_load!, PQ_dyn_load!, new_PQ_dyn_load!]
    control_sig_syms::Vector{Symbol} = Symbol[ ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
    
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]
    cb_state_values::Vector{Function}          = Function[]    
    # cb_state_values::Vector{Float64}           = Float64[]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    # cb_dyn_state_values::Vector{Float64}       = Float64[]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)    

end

#-------------------------------------------------------


function loc_Load_t1!(dx, x, p_agg, t)
    

    cb_sw, u_gen_ur_ui, f_t, p =  p_agg
    
    loc_P   = p[1]
    
    loc_Q   = p[2]

    u_r, u_i = u_gen_ur_ui

    u = u_r + u_i * im

    u_mag = abs(u)
   
    i_r, i_i = x
        
    i = i_r + i_i * im

    S = loc_P + im * loc_Q

    I = conj(S) * u / (u_mag)^2
    
    dx[1] = real(I) - i_r
    
    dx[2] = imag(I) - i_i

    #---------------------------      
    
    return nothing
end

#####

"""

This structure is used to organise data for local loads.
    """

@kwdef struct loc_Load_t1  <: SdNonGen #PQNode Loc_load

    Bus::String  = "bus0"
    name::String = lowercase(Bus)     
    loc_P::Float64 = 0.0
    loc_Q::Float64 = 0.0

    comp_type::Symbol = :loc_Load_t1

    Bus_type::Symbol  = :loc_Load 

    state_vars_syms::Vector{Symbol} = Symbol[  ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :i_r, :i_i  ]
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[ ]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    
    param::Vector{Symbol} = Symbol[:loc_P, :loc_Q ]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :loc_P ], dict_param_syms_Idx[ :loc_Q ]]  
    
    param_values::Vector{Float64} = Float64[loc_P, loc_Q ]

    func::Vector{Function} = Function[ loc_Load_t1!, loc_Load_t1!, loc_Load_t1!, loc_Load_t1! ]

    control_sig_syms::Vector{Symbol} = Symbol[ ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[ :i_r, :i_i ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
    
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]
    cb_state_values::Vector{Function}          = Function[]    
    # cb_state_values::Vector{Float64}           = Float64[]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    # cb_dyn_state_values::Vector{Float64}       = Float64[]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)    

end



# ------------------------------------------------------
# Transmission nodes
# ------------------------------------------------------



function hybrid_pf_Trans_t1_Node!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i, f_t, p,
     global_pf_param,
     node_pf_param) = p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr,
     node_inc_edges_orient,
     nodes_u_view)  = node_pf_param
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]

    S = P + im * Q

   #----------------------------

    u_r,  u_i = x
    
    u = u_r + u_i * im     
    
    #----------------------------

    
    my_node_idx =
        node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j =
        node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )
    
    sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj)])
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -  sum_ykj_vj )

    # calculating Uk twice
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -
        sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]))

    # ------------------------------------------------------
    
    
    dx[1] = real(Uk) - u_r
    
    dx[2] = imag(Uk) - u_i

   
    return nothing
end


function node_pf_Trans_t1_Node!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i,
     f_t, p, node_pf_param) = p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr,
     node_inc_edges_orient,
     nodes_u_view)  = node_pf_param
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]

    S = P + im * Q

   #----------------------------

    u_r,  u_i = x
    
    u = u_r + u_i * im     
    
    #----------------------------

    
    my_node_idx =
        node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j =
        node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )
    
    sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj)])
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -  sum_ykj_vj )

    # calculating Uk twice
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -  sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]))

    # ------------------------------------------------------
    
    
    dx[1] = real(Uk) - u_r
    
    dx[2] = imag(Uk) - u_i

   
    return nothing
end



function global_pf_Trans_t1_Node!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i,
     f_t, p, global_pf_param) = p_agg

    pf_U, Inet = global_pf_param
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]

    y_shunt = im * Y_n    

    u_r, u_i  = x
    
    u = u_r + u_i * im
    
    U_r = pf_U[1]
    
    U_i = pf_U[2]
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i
   
    return nothing
end


function network_current_Trans_t1_Node!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i,
     f_t, p) = p_agg
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]

    y_shunt = im * Y_n    

    # vh = abs(U)
    # θh = angle(U)

    vh, θh, ph, qh, i_r, i_i  = f_t
    
    uh  =  vh * exp(im * θh)

    i  = 1.0 * dynamic_nodal_current_balance(src_i,dst_i) +
        uh * y_shunt  
        
    u_r, u_i  = x
    
    # dx[1] = vh * cos(θh) - u_r
    
    # dx[2] = vh * sin(θh) - u_i


    dx[1]   = real(uh) - u_r
    
    dx[2]   = imag(uh) - u_i     
     
     
    return nothing
end


function initial_pf_Trans_t1_Node!(dx, x, p_agg, t)

    cb_sw, src_i, dst_i, f_t, p = p_agg
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]

    y_shunt = im * Y_n    

    vh, θh, ph, qh, i_r, i_i  = f_t
    
    uh  =  vh * exp(im * θh)
         
    u_r, u_i  = x
    
    # dx[1]   = vh * cos(θh)  - u_r
    
    # dx[2]   = vh * sin(θh)  - u_i

    dx[1]   = real(uh) - u_r
    
    dx[2]   = imag(uh) - u_i    
    

   
    return nothing
end


@kwdef struct Trans_t1_Node   <: SdNonGen

    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64       = 0.0
    Q::Float64       = 0.0
    Y_n              = 0.0

    vmax::Float64    = 1.06
    vmin::Float64    = 0.97

    comp_type::Symbol = :Trans_t1_Node

    Bus_num::Int64 =
        parse(Int, split(lowercase(Bus),"bus")[2] )
    Bus_type::Symbol  =  :Transmission

    state_vars_syms::Vector{Symbol} = Symbol[  ]
    algebraic_vars_syms::Vector{Symbol} =
        Symbol[ :u_r, :u_i ]
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} =
        Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    param::Vector{Symbol} =
        Symbol[:P, :Q, :Y_n, :vmin, :vmax ]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]  
    
    param_values::Vector{Float64} = Float64[P, Q, Y_n, vmin, vmax]
    func::Vector{Function} = Function[initial_pf_Trans_t1_Node!, network_current_Trans_t1_Node!, global_pf_Trans_t1_Node!, node_pf_Trans_t1_Node!, hybrid_pf_Trans_t1_Node! ]
    control_sig_syms::Vector{Symbol} = Symbol[]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[ ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)
    
end

 

function hybrid_pf_Trans_t2_Node!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i,
     f_t, p, dyn_global_pf_param,
     dyn_node_pf_param) = p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr, node_inc_edges_orient,
     nodes_u_view)  = dyn_node_pf_param
    
    pf_U, Inet = dyn_global_pf_param
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]

    S = P + im * Q

    y_shunt = im * Y_n    

   #----------------------------

    u_r,  u_i = x
    
    u = u_r + u_i * im     
    
    #----------------------------

    
    my_node_idx =
        node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j =
        node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.( nodes_u_view[ nodes_j ] )

    Uk  = u_from_ur_ui( nodes_u_view[ my_node_idx ] )
    
    sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj)])
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -  sum_ykj_vj )

    # calculating Uk twice
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -
        sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]))

    # ------------------------------------------------------
        
    dx[1] = real(Uk) - u_r
    
    dx[2] = imag(Uk) - u_i


    return nothing
end



function node_pf_Trans_t2_Node!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i,
     f_t, p, node_pf_param) = p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr,
     node_inc_edges_orient,
     nodes_u_view) = node_pf_param
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]

    S = P + im * Q

    y_shunt = im * Y_n    

   #----------------------------

    u_r,  u_i = x
    
    u = u_r + u_i * im     
    
    #----------------------------

    
    my_node_idx =
        node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j =
        node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )
    
    sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj)])
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -  sum_ykj_vj )

    # calculating Uk twice
    
    Uk  = (1/Ykk) * (conj(S)/conj(Uk) -  sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]))

    # ------------------------------------------------------
        
    dx[1] = real(Uk) - u_r
    
    dx[2] = imag(Uk) - u_i


    return nothing
end



function global_pf_Trans_t2_Node!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i,
     f_t, p, global_pf_param)  = p_agg

    pf_U, Inet = global_pf_param
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]

    y_shunt = im * Y_n    

    u_r, u_i  = x
    
    u = u_r + u_i * im     
   
    U_r = pf_U[1]
    
    U_i = pf_U[2]
    
    dx[1] = U_r - u_r
    
    dx[2] = U_i - u_i
    

    return nothing
end


function network_current_Trans_t2_Node!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i,
     f_t, p) = p_agg
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]

    y_shunt = im * Y_n    

    # vh = abs(U)
    # θh = angle(U)

    vh, θh, ph, qh, i_r, i_i = f_t

    uh  =  vh * exp(im * θh)

    i  = 1.0 * dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt 
        
    u_r, u_i  = x
        
    # dx[1]   = vh * cos(θh)  - u_r
    
    # dx[2]   = vh * sin(θh)  - u_i 

    dx[1]   = real(uh) - u_r
    
    dx[2]   = imag(uh) - u_i     
   
    return nothing
end



function initial_pf_Trans_t2_Node!(dx, x, p_agg, t)

    (cb_sw, src_i, dst_i,
     f_t, p) = p_agg
    
    P     = p[1]
    Q     = p[2]
    Y_n   = p[3]


    # vh = abs(U)
    # θh = angle(U)

    vh, θh, ph, qh, i_r, i_i  = f_t

    uh  =  vh * exp(im * θh)
        
    u_r, u_i  = x
        
    # dx[1]   = vh * cos(θh)  - u_r
    
    # dx[2]   = vh * sin(θh)  - u_i

  
    dx[1]   = real(uh) - u_r
    
    dx[2]   = imag(uh) - u_i 
        
   
    return nothing
end


@kwdef struct Trans_t2_Node   <: SdNonGen

    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64       = 0.0
    Q::Float64       = 0.0
    Y_n              = 0.0

    vmax::Float64    = 1.06
    vmin::Float64    = 0.97

    comp_type::Symbol = :Trans_t2_Node
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )
    Bus_type::Symbol  = :Transmission
    
    state_vars_syms::Vector{Symbol} = Symbol[  ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i  ]
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    param::Vector{Symbol} = Symbol[:P, :Q, :Y_n, :vmin, :vmax ]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]  
    
    param_values::Vector{Float64} = Float64[P, Q, Y_n, vmin, vmax]
    func::Vector{Function} = Function[initial_pf_Trans_t2_Node!, network_current_Trans_t2_Node!, global_pf_Trans_t2_Node!, node_pf_Trans_t2_Node!, hybrid_pf_Trans_t2_Node!]
    control_sig_syms::Vector{Symbol} = Symbol[]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[ ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)
    
end


# ------------------------------------------------------
# Generators nodes
# ------------------------------------------------------


"""
Sauer: page 35

vh * exp(im * θh) = (vd +im * vq) * exp(im * (δ - π/2))
ih * exp(im * ϕh) = (id +im * iq) * exp(im * (δ - π/2))

Sauer: page 51 and Millano: page 223

vd = -ra * id + ed_dash + X_q_dash * iq
vq = -ra * iq + eq_dash - X_d_dash * id

i.e.

vd = ed_dash - (ra * id  - X_q_dash * iq)

vq = eq_dash - (X_d_dash * id  + ra * iq)

τe = (vd + ra * id) * id +  (vq + ra * iq) * iq

τe = (ed_dash - (ra * id  - X_q_dash * iq) + ra * id) * id +
(eq_dash - (X_d_dash * id  + ra * iq) + ra * iq)  * iq

= ed_dash * id - ra * id^2 + X_q_dash * iq *id +ra * id^2 +
eq_dash * iq - X_d_dash * id * iq - ra * iq^2 + ra * iq^2

= ed_dash * id  + eq_dash * iq + (X_q_dash - X_d_dash) * iq * id



# ------------------------------------------------------

vh * exp(im * θh) = (vd + im * vq) * exp(im * (δ - π/2))

vh * cos(θh) + im * vh * sin(θh) =  (vd + im * vq) *  exp(im * - π/2) * exp(im * δ)

vh * cos(θh) + im * vh * sin(θh) = -im * (vd + im * vq) * exp(im * δ)

vh * cos(θh) + im * vh * sin(θh) = -im * (vd + im * vq) * (cos(δ) + im * sin(δ))

vh * cos(θh) + im * vh * sin(θh) = (vd + im * vq) * (-im * cos(δ) -im * im * sin(δ))

vh * cos(θh) + im * vh * sin(θh) = ( vd + im * vq ) * (sin(δ) -im * cos(δ) )

vh * cos(θh) + im * vh * sin(θh) = vd * sin(δ) + vq * cos(δ) + im * (vq * sin(δ) - vd * cos(δ))

# ------------------------------------------------------

uh_r = vh * cos(θh) = vd * sin(δ) + vq * cos(δ)

uh_i = vh * sin(θh) = vq * sin(δ) - vd * cos(δ)

# ------------------------------------------------------

uh_r = vd * sin(δ) + vq * cos(δ) = (ed_dash - (ra * id  - X_q_dash * iq) ) * sin(δ) + ( eq_dash - (X_d_dash * id  + ra * iq ) ) * cos(δ)

uh_r =  ed_dash *  sin(δ) -  ra * id * sin(δ) +  X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ)

uh_i = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ)

in compact form

[uh_r, uh_i] = [sin(δ) cos(δ); -cos(δ) sin(δ)] * [ed_dash - (ra * id  - X_q_dash * iq) , eq_dash - (X_d_dash * id  + ra * iq)]

uh_r and  uh_i are algebraic variables

du .= [sin(δ) cos(δ); -cos(δ) sin(δ)] * [ed_dash - (ra * id  - X_q_dash * iq) , eq_dash - (X_d_dash * id  + ra * iq)] - [uh_r, uh_i]

# ------------------------------------------------------
Sauer: page 160

eq 7.24 and eq 7.25

# ------------------------------------------------------


The equation below is wrong!!!

uh_r = ((X_q_dash * iq - ra * id) + ed_dash) * sin(δ)

uh_i = ((ra * id - X_q_dash * iq) - ed_dash) * cos(δ)

    
"""

#---------------------------------------------------
# SM_2axis_cb_v6
#---------------------------------------------------



function hybrid_pf_SM_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]
    vh        = p[22]

    y_shunt   = im * Y_n

    M         =  ( 2 * H ) / ωs
          
    
    # States
    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
          
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt    

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    

    # ------------------------------------------------------


    dx[1] =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i

    return nothing
end



function node_pf_SM_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q        = p[21]
    vh       = p[22]    

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91  

    # # ------------------------------------------------------


    # my_node_idx = node_idx_and_inc_edges[1]
    
    # inc_edges_orientation = last.(node_inc_edges_Ybr_orient)

    # # Get the node number of the other end of an edge
    # # that is connected to me.
    
    # nodes_j = [my_node_idx == orient[1] ? orient[2] : orient[1]   for orient in inc_edges_orientation ]

    # edges_Ybr = first.(node_inc_edges_Ybr_orient)

    # # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length(edges_Ybr )]
    
    # Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( inc_edges_orientation , edges_Ybr) ] ) : my_node_idx == first( inc_edges_orientation[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    # Ykk = Ybus_node[1]
    
    # Ykj = Ybus_node[2:end]
    
    # Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    # sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj)])

    # # ------------------------------------------------------


   
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
        
    i =  Ykk * uh + sum_ykj_vj + uh * y_shunt 

    i_mag = abs(i)

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)        

    # ph = (ed_dash - (ra * id - X_q_dash * iq)) * id + (eq_dash - (X_d_dash * id  + ra * iq)) * iq
    
    # qh  = (eq_dash - (X_d_dash * id  + ra * iq)) * id - (ed_dash - (ra * id - X_q_dash * iq) ) * iq    
    

    # Vk = (ph + im * qh - vh^2 * Ykk )/ conj(sum_ykj_vj)

    # Vk = vh * exp(im * angle( Vk ))

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    dx[1] = (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash


    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
     return nothing
end


function global_pf_SM_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param = p_agg

    pf_U, Inet = global_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q        = p[21]
    
    vh       = p[22]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91   

    i  = Inet + uh * y_shunt    

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # ------------------------------------------------------
    
    U_r = pf_U[1]

    U_i = pf_U[2]

    # ------------------------------------------------------
    

    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    dx[1] = (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash


    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
     return nothing
end


function network_current_SM_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q        = p[21]
    vh       = p[22]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
          
    
    # States
    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91   

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt    

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)


    dx[1] =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M #  (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i

    return nothing
end



function initial_pf_SM_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91

    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt    

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    dx[1] =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq )  - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end


@kwdef struct SM_2axis_cb_v6 <: SdGen

    Bus::String   = "bus0"
    name::String  = lowercase(Bus)    
    P::Float64 = 0.0          
    D::Float64 = 0.0         
    H::Float64 = 0.0         
    Ωb::Float64 = 0.0        
    ωs::Float64 = 0.0        
    ra::Float64 = 0.0        
    xℓ::Float64 = 0.0        
    X_d::Float64 = 0.0       
    X_q::Float64  = 0.0      
    X_d_dash::Float64 = 0.0  
    X_q_dash::Float64 = 0.0  
    X_d_2dash::Float64 = 0.0 
    X_q_2dash::Float64 = 0.0 
    T_d_dash::Float64 = 0.0  
    T_q_dash::Float64 = 0.0  
    T_d_2dash::Float64 = 0.0 
    T_q_2dash::Float64 = 0.0 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0

    Q::Float64 = 0.0

    Sn::Float64 = 0.0

    vh::Float64  = 1.0

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_cb_v6

    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator
    
    state_vars_syms::Vector{Symbol} = Symbol[  :δ, :ω, :ed_dash, :eq_dash]
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[:u_r, :u_i ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    # state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω,:ed_dash,:eq_dash, :u_r, :u_i ]
    # algebraic_vars_syms::Vector{Symbol} = Symbol[   ]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[
        state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :vh, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64} = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, vh, Sn]
    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ]
    
    func::Vector{Function} = Function[ initial_pf_SM_2axis_cb_v6!, network_current_SM_2axis_cb_v6!, global_pf_SM_2axis_cb_v6!, node_pf_SM_2axis_cb_v6!, hybrid_pf_SM_2axis_cb_v6! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
     output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end


#---------------------------------------------------
# SM_2axis_cb_idq
#---------------------------------------------------



function hybrid_pf_SM_2axis_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters

    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt   = im * Y_n

    M         =  ( 2 * H ) / ωs
    
    # States
    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
       
    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt
   
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    

    # ------------------------------------------------------

    # State equations

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
    
    return nothing
end


function node_pf_SM_2axis_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q        = p[21]
    vh       = p[22]    

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    # u_r  = x[1],  u_i     = x[2], δ       = x[3], 
    # ω    = x[4],  ed_dash = x[5], eq_dash = x[6]

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
   
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
    
    i =  Ykk * Uk + sum_ykj_vj  + Uk * y_shunt

    i_mag = abs(i)

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)        

    # ph = (ed_dash - (ra * id - X_q_dash * iq)) * id + (eq_dash - (X_d_dash * id  + ra * iq)) * iq
    
    # qh  = (eq_dash - (X_d_dash * id  + ra * iq)) * id - (ed_dash - (ra * id - X_q_dash * iq) ) * iq    
    

    # Vk = (ph + im * qh - vh^2 * Ykk )/ conj(sum_ykj_vj)

    # Vk = vh * exp(im * angle( Vk ))

    
    # ------------------------------------------------------
    

    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[3] = dδ =  Ωb * (ω - ωs)
    
    # dx[4] = (τm_gov - τe - D * (ω - ωs)) *  (Ωb / (2*H))

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end



function global_pf_SM_2axis_cb_idq!(dx, x, p_agg, t)

    #u_gov, u_exc, f_t, p = p_agg
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param = p_agg

    pf_U, Inet = global_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q        = p[21]
    vh       = p[22]    

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    # u_r  = x[1],  u_i     = x[2], δ       = x[3], 
    # ω    = x[4],  ed_dash = x[5], eq_dash = x[6]

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91
    
    i = Inet + uh * y_shunt    

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    

    # ------------------------------------------------------

    U_r = pf_U[1]

    U_i = pf_U[2]
    
    # ------------------------------------------------------
    

    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[3] = dδ =  Ωb * (ω - ωs)
    
    # dx[4] = (τm_gov - τe - D * (ω - ωs)) *  (Ωb / (2*H))

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end



function network_current_SM_2axis_cb_idq!(dx, x, p_agg, t)

    #u_gov, u_exc, f_t, p = p_agg
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    # u_r  = x[1],  u_i     = x[2], δ       = x[3], 
    # ω    = x[4],  ed_dash = x[5], eq_dash = x[6]

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91
   
    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt
   
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[3] = dδ =  Ωb * (ω - ωs)
    
    # dx[4] = (τm_gov - τe - D * (ω - ωs)) *  (Ωb / (2*H))

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end


function initial_pf_SM_2axis_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    # u_r  = x[1],  u_i     = x[2], δ       = x[3], 
    # ω    = x[4],  ed_dash = x[5], eq_dash = x[6]

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
     δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91

    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt       

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[3] = dδ =  Ωb * (ω - ωs)
    
    # dx[4] = (τm_gov - τe - D * (ω - ωs)) *  (Ωb / (2*H))

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end



@kwdef struct SM_2axis_cb_idq <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64 = 0.0          
    D::Float64         
    H::Float64         
    Ωb::Float64        
    ωs::Float64        
    ra::Float64        
    xℓ::Float64        
    X_d::Float64       
    X_q::Float64       
    X_d_dash::Float64  
    X_q_dash::Float64  
    X_d_2dash::Float64 
    X_q_2dash::Float64 
    T_d_dash::Float64  
    T_q_dash::Float64  
    T_d_2dash::Float64 
    T_q_2dash::Float64 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0

    Q::Float64 = 0.0

    Sn::Float64 = 0.0

    vh::Float64  = 1.0

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_cb_idq
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator

    state_vars_syms::Vector{Symbol} = Symbol[  :δ, :ω, :ed_dash, :eq_dash]
    algebraic_vars_syms::Vector{Symbol} = Symbol[:u_r, :u_i ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]

    # state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]
    # algebraic_vars_syms::Vector{Symbol} = Symbol[  ]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :vh, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64}  = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, vh, Sn]

    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ]

    func::Vector{Function} = Function[ initial_pf_SM_2axis_cb_idq!, network_current_SM_2axis_cb_idq!, global_pf_SM_2axis_cb_idq!, node_pf_SM_2axis_cb_idq!, hybrid_pf_SM_2axis_cb_idq! ]

    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr    
    
end


#---------------------------------------------------
# SC_2axis_cb_v6
#---------------------------------------------------



function hybrid_pf_SC_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters

    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
        
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    # Algebraics equations
    
    # Network interface
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # State equations
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
    
    return nothing
end


function node_pf_SC_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    # # ------------------------------------------------------

    # my_node_idx = node_idx_and_inc_edges[1]
    
    # inc_edges_orientation = last.(node_inc_edges_Ybr_orient)

    # # Get the node number of the other end of an edge
    # # that is connected to me.
    
    # nodes_j = [my_node_idx == orient[1] ? orient[2] : orient[1]   for orient in inc_edges_orientation ]

    # edges_Ybr = first.(node_inc_edges_Ybr_orient)

    # # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length(edges_Ybr )]

    # Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( inc_edges_orientation , edges_Ybr) ] ) : my_node_idx == first( inc_edges_orientation[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]
    

    # Ykk = Ybus_node[1]
    # Ykj = Ybus_node[2:end]
    # Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    # sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 


   
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------    
    
    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt 

    i_mag = abs(i)

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)        
    
    # qh  = (eq_dash - (X_d_dash * id  + ra * iq)) * id - (ed_dash - (ra * id - X_q_dash * iq) ) * iq    
    

    #  Vk = (im * qh - vh^2 * Ykk )/ sum_ykj_vj 

    # ------------------------------------------------------
    
    dx[1] = (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # dx[5] = real(Vk ) - u_r
    
    # dx[6] = imag(Vk ) - u_i 
    
     
    return nothing
end



function global_pf_SC_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p, global_pf_param = p_agg

    pf_U, Inet = global_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    i  = Inet  + uh * y_shunt

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    U_r = pf_U[1]

    U_i = pf_U[2]
    
    # ------------------------------------------------------
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M # ( 1 /(2*H) )    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # dx[5] = U_r - u_r
    
    # dx[6] = U_i - u_i  
     
    return nothing
end



function network_current_SC_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
        
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # τe = P
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # State equations
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  ( 1 /(2*H) )    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end


function initial_pf_SC_2axis_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]
    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current
    
    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt 


    # Algebraics equations
    # Network interface
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i  
    
     
    return nothing
end



@kwdef struct SC_2axis_cb_v6 <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64  = 0.0         
    D::Float64   = 0.0        
    H::Float64   = 0.0        
    Ωb::Float64   = 0.0       
    ωs::Float64   = 0.0       
    ra::Float64  = 0.0        
    xℓ::Float64  = 0.0        
    X_d::Float64  = 0.0       
    X_q::Float64   = 0.0      
    X_d_dash::Float64  = 0.0  
    X_q_dash::Float64  = 0.0  
    X_d_2dash::Float64  = 0.0 
    X_q_2dash::Float64  = 0.0 
    T_d_dash::Float64   = 0.0 
    T_q_dash::Float64   = 0.0 
    T_d_2dash::Float64  = 0.0 
    T_q_2dash::Float64  = 0.0 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0
    Q::Float64 = 0.0

    Sn::Float64 = 0.0

    vh::Float64  = 1.0

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SC_2axis_cb_v6
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator

    state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i  ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :vh, :Sn]    

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64}  = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, vh, Sn]

    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ]
    
    func::Vector{Function} = Function[ initial_pf_SC_2axis_cb_v6!, network_current_SC_2axis_cb_v6!, global_pf_SC_2axis_cb_v6!, node_pf_SC_2axis_cb_v6!, hybrid_pf_SC_2axis_cb_v6! ]

    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    Ax::Function         = Ax_SC_gen
    Bx::Function         = Bx_SC_gen
    Cx::Function         = Cx_SC_gen
    Ax_τm_vf::Function   = Ax_SC_gen_τm_vf    

    Ax_gen_avr::Function = Ac_SC_gen_avr
    
    Ax_gen_S_avr_S::Function = SC_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SC_Ax_gen_S_avr_A

    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end


#---------------------------------------------------
# SC_2axis_cb_idq
#---------------------------------------------------


function hybrid_pf_SC_2axis_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
        
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    # Algebraics equations
    
    # Network interface
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # State equations
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
    
     
    return nothing
end


function node_pf_SC_2axis_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    # # ------------------------------------------------------

    # my_node_idx = node_idx_and_inc_edges[1]
    
    # inc_edges_orientation = last.(node_inc_edges_Ybr_orient)

    # # Get the node number of the other end of an edge
    # # that is connected to me.
    
    # nodes_j = [my_node_idx == orient[1] ? orient[2] : orient[1]   for orient in inc_edges_orientation ]

    # edges_Ybr = first.(node_inc_edges_Ybr_orient)

    # # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length(edges_Ybr )]

    # Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( inc_edges_orientation , edges_Ybr) ] ) : my_node_idx == first( inc_edges_orientation[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]
    
    # Ykk = Ybus_node[1]
    # Ykj = Ybus_node[2:end]
    # Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    # sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 
       
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
    
    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt 

    i_mag = abs(i)

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)        
    
    # qh  = (eq_dash - (X_d_dash * id  + ra * iq)) * id - (ed_dash - (ra * id - X_q_dash * iq) ) * iq    
    

    #  Vk = (im * qh - vh^2 * Ykk )/ sum_ykj_vj 

    # ------------------------------------------------------
    
    dx[1] = (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # dx[5] = real(Vk ) - u_r
    
    # dx[6] = imag(Vk ) - u_i 
    
     
    return nothing
end


function global_pf_SC_2axis_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p, global_pf_param = p_agg

    pf_U, Inet = global_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    i  = Inet  + uh * y_shunt

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    U_r = pf_U[1]

    U_i = pf_U[2]
    
    # ------------------------------------------------------
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M # ( 1 /(2*H) )    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # dx[5] = U_r - u_r
    
    # dx[6] = U_i - u_i  
     
    return nothing
end



function network_current_SC_2axis_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
        
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # τe = P
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # State equations
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  ( 1 /(2*H) )    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end


function initial_pf_SC_2axis_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]
    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current
    
    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt 


    # Algebraics equations
    # Network interface
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i  
    
     
    return nothing
end


@kwdef struct SC_2axis_cb_idq <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64  = 0.0         
    D::Float64         
    H::Float64         
    Ωb::Float64        
    ωs::Float64        
    ra::Float64        
    xℓ::Float64        
    X_d::Float64       
    X_q::Float64       
    X_d_dash::Float64  
    X_q_dash::Float64  
    X_d_2dash::Float64 
    X_q_2dash::Float64 
    T_d_dash::Float64  
    T_q_dash::Float64  
    T_d_2dash::Float64 
    T_q_2dash::Float64 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0
    Q::Float64 = 0.0

    Sn::Float64 = 0.0

    vh::Float64    = 1.06

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SC_2axis_cb_idq
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator

    # state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i]
    # algebraic_vars_syms::Vector{Symbol} = Symbol[ ]

    state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i  ]

    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]

    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64}  = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, Sn]

    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ]
    
    func::Vector{Function} = Function[ initial_pf_SC_2axis_cb_idq!, network_current_SC_2axis_cb_idq!, global_pf_SC_2axis_cb_idq!, node_pf_SC_2axis_cb_idq!, hybrid_pf_SC_2axis_cb_idq! ]
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax::Function         = Ax_SC_gen
    Bx::Function         = Bx_SC_gen
    Cx::Function         = Cx_SC_gen
    Ax_τm_vf::Function   = Ax_SC_gen_τm_vf    

    Ax_gen_avr::Function = Ac_SC_gen_avr
    
    Ax_gen_S_avr_S::Function = SC_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SC_Ax_gen_S_avr_A

    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end



#---------------------------------------------------
# SM_2axis_wt_loc_load_cb_v6
#---------------------------------------------------



function hybrid_pf_SM_2axis_wt_loc_load_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters

    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt   = im * Y_n

    M         =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    i  = 1.0 * dynamic_nodal_current_balance(src_i, dst_i) + (uh * y_shunt + loc_i_r + loc_i_i * im)


    # Algebraics equations
    
    # Network interface
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    return nothing
end



function node_pf_SM_2axis_wt_loc_load_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    # # ------------------------------------------------------

    # my_node_idx = node_idx_and_inc_edges[1]
    
    # inc_edges_orientation = last.(node_inc_edges_Ybr_orient)

    # # Get the node number of the other end of an edge
    # # that is connected to me.
    
    # nodes_j = [my_node_idx == orient[1] ? orient[2] : orient[1]   for orient in inc_edges_orientation ]

    # edges_Ybr = first.(node_inc_edges_Ybr_orient)

    # # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length(edges_Ybr )]

    # Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( inc_edges_orientation , edges_Ybr) ] ) : my_node_idx == first( inc_edges_orientation[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    # Ykk = Ybus_node[1]
    
    # Ykj = Ybus_node[2:end]
    
    # Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    # sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj)])


   
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
        
    i =  Ykk * uh + sum_ykj_vj + uh * y_shunt + loc_i_r + loc_i_i * im

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)        

    # ph = (ed_dash - (ra * id - X_q_dash * iq)) * id + (eq_dash - (X_d_dash * id  + ra * iq)) * iq
    
    # qh  = (eq_dash - (X_d_dash * id  + ra * iq)) * id - (ed_dash - (ra * id - X_q_dash * iq) ) * iq    
    

    # Vk = (ph + im * qh - vh^2 * Ykk )/ conj(sum_ykj_vj)

    # Vk = vh * exp(im * angle( Vk ))

    # ------------------------------------------------------
    
    # τe = P
    
    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))   
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash


    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
    
    # dx[5] = real(Vk ) - u_r
    
    # dx[6] = imag(Vk ) - u_i   
      
    return nothing
end



function global_pf_SM_2axis_wt_loc_load_cb_v6!(dx,x,p_agg,t)
    
    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, global_pf_param = p_agg

    pf_U, Inet = global_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii
    
    i  =  Inet + uh * y_shunt + loc_i_r + loc_i_i * im 

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    
    U_r = pf_U[1]

    U_i = pf_U[2]

    # ------------------------------------------------------
    
    # τe = P
    
    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # dx[5] = U_r - u_r
    
    # dx[6] = U_i - u_i
    
    return nothing
    
end


function network_current_SM_2axis_wt_loc_load_cb_v6!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    i  = 1.0 * dynamic_nodal_current_balance(src_i, dst_i) + (uh * y_shunt + loc_i_r + loc_i_i * im)


    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # τe = P
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # State equations
    
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs)) * (Ωb / (2*H))
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)    
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M #  (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
 
     
    return nothing
end


function initial_pf_SM_2axis_wt_loc_load_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]
    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii
    
    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt + loc_i_r + loc_i_i * im

    # Algebraics equations
    # Network interface
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs)) * (Ωb / (2*H))
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - (ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq) - D * (ω - ωs)) / M  # * (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end



@kwdef struct SM_2axis_wt_loc_load_cb_v6 <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64  = 0.0         
    D::Float64   = 0.0       
    H::Float64  = 0.0        
    Ωb::Float64 = 0.0        
    ωs::Float64 = 0.0        
    ra::Float64 = 0.0        
    xℓ::Float64  = 0.0       
    X_d::Float64 = 0.0       
    X_q::Float64  = 0.0      
    X_d_dash::Float64  = 0.0 
    X_q_dash::Float64  = 0.0 
    X_d_2dash::Float64 = 0.0 
    X_q_2dash::Float64 = 0.0 
    T_d_dash::Float64  = 0.0 
    T_q_dash::Float64  = 0.0 
    T_d_2dash::Float64 = 0.0 
    T_q_2dash::Float64 = 0.0 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0
    
    Q::Float64       = 0.0

    Sn::Float64 = 0.0

    vh::Float64      = 0.0

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.97    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_wt_loc_load_cb_v6

    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator

    state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))
    
    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :vh, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]
    
    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64}  = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, vh, Sn]

    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ]
    
    func::Vector{Function} = Function[ initial_pf_SM_2axis_wt_loc_load_cb_v6!, network_current_SM_2axis_wt_loc_load_cb_v6!, global_pf_SM_2axis_wt_loc_load_cb_v6!, node_pf_SM_2axis_wt_loc_load_cb_v6!, hybrid_pf_SM_2axis_wt_loc_load_cb_v6! ]
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end


#---------------------------------------------------
# SC_2axis_wt_loc_load_cb_v6
#---------------------------------------------------



function hybrid_pf_SC_2axis_wt_loc_load_cb_v6!(dx, x, p_agg, t)

    
    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt   = im * Y_n

    M         =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    i  = 1.0 * dynamic_nodal_current_balance(src_i, dst_i) + ( uh * y_shunt + loc_i_r + loc_i_i * im)

    # Algebraics equations
    
    # Network interface

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # State equations
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
     
    return nothing
end



function node_pf_SC_2axis_wt_loc_load_cb_v6!(dx, x, p_agg, t)

    
    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    # # --------------------------------------------------

    # my_node_idx = node_idx_and_inc_edges[1]
    
    # inc_edges_orientation = last.(node_inc_edges_Ybr_orient)

    # # Get the node number of the other end of an edge
    # # that is connected to me.
    
    # nodes_j = [my_node_idx == orient[1] ? orient[2] : orient[1]   for orient in inc_edges_orientation ]

    # edges_Ybr = first.(node_inc_edges_Ybr_orient)

    # # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length(edges_Ybr )]
    
    # Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( inc_edges_orientation , edges_Ybr) ] ) : my_node_idx == first( inc_edges_orientation[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    # Ykk = Ybus_node[1]
    
    # Ykj = Ybus_node[2:end]
    
    # Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    # sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
        
    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt + loc_i_r + loc_i_i * im

    i_mag = abs(i)

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)        
    
    # qh  = (eq_dash - (X_d_dash * id  + ra * iq)) * id - (ed_dash - (ra * id - X_q_dash * iq) ) * iq    
    

    #  Vk = (im * qh - vh^2 * Ykk )/ sum_ykj_vj 

    # ------------------------------------------------------
        
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    

    # dx[5] = real(Vk ) - u_r
    
    # dx[6] = imag(Vk ) - u_i 
    
     
    return nothing
end



function global_pf_SC_2axis_wt_loc_load_cb_v6!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, global_pf_param = p_agg

    
    pf_U, Inet = global_pf_param

      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    i  = Inet + ( uh * y_shunt + loc_i_r + loc_i_i * im)

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)

    # ------------------------------------------------------
    
    U_r = pf_U[1]

    U_i = pf_U[2]    

    # ------------------------------------------------------

    # τe = P
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # State equations
    
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs)) * (Ωb / (2*H))
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M # (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    

    # dx[5] = U_r - u_r
    
    # dx[6] = U_i - u_i  
        
    return nothing
end



function network_current_SC_2axis_wt_loc_load_cb_v6!(dx, x, p_agg, t)

    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii
    
    # i_r, i_i = f_t
    
    # i = i_r + im * i_i + uh * y_shunt +
    #     loc_i_r + loc_i_i * im

    i  = 1.0 * dynamic_nodal_current_balance(src_i, dst_i) + ( uh * y_shunt + loc_i_r + loc_i_i * im)

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # τe = P
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # State equations
    
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs)) * (Ωb / (2*H))
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
     
    return nothing
end



function initial_pf_SC_2axis_wt_loc_load_cb_v6!(dx, x, p_agg, t)

    #u_gov, u_exc, f_t, p = p_agg
    
    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]
    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash,  u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii
    
    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt + loc_i_r + loc_i_i * im

    # Algebraics equations
    # Network interface
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    

    dx[1] =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M # (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i 
    
     
    return nothing
end



@kwdef struct SC_2axis_wt_loc_load_cb_v6 <: SdGen
    Bus::String   = "bus0"
    name::String  = lowercase(Bus)     
    P::Float64    = 0.0         
    D::Float64     = 0.0       
    H::Float64   = 0.0        
    Ωb::Float64   = 0.0       
    ωs::Float64   = 0.0       
    ra::Float64   = 0.0       
    xℓ::Float64   = 0.0       
    X_d::Float64   = 0.0      
    X_q::Float64   = 0.0      
    X_d_dash::Float64   = 0.0 
    X_q_dash::Float64   = 0.0 
    X_d_2dash::Float64  = 0.0 
    X_q_2dash::Float64  = 0.0 
    T_d_dash::Float64   = 0.0 
    T_q_dash::Float64   = 0.0 
    T_d_2dash::Float64  = 0.0 
    T_q_2dash::Float64  = 0.0 
    αp::Float64  = 1.0   
    αq::Float64  = 1.0   
    Y_n::Float64 = 0.0
    
    Q::Float64       = 0.0

    Sn::Float64 = 0.0
    
    vh::Float64      = 1.0

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SC_2axis_wt_loc_load_cb_v6
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator

    state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[  :u_r, :u_i  ]

    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))
    
    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :vh, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]

    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64}  = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, vh, Sn]

    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ]

    func::Vector{Function} = Function[ initial_pf_SC_2axis_wt_loc_load_cb_v6!, network_current_SC_2axis_wt_loc_load_cb_v6!, global_pf_SC_2axis_wt_loc_load_cb_v6!, node_pf_SC_2axis_wt_loc_load_cb_v6!, hybrid_pf_SC_2axis_wt_loc_load_cb_v6! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    Ax::Function         = Ax_SC_gen
    Bx::Function         = Bx_SC_gen
    Cx::Function         = Cx_SC_gen
    Ax_τm_vf::Function   = Ax_SC_gen_τm_vf    

    Ax_gen_avr::Function = Ac_SC_gen_avr
    
    Ax_gen_S_avr_S::Function = SC_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SC_Ax_gen_S_avr_A

    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end


#---------------------------------------------------
# Others
#---------------------------------------------------


function hybrid_pf_SM_2axis_idq!(dx, x, p_agg, t)
    
    (cb_sw, τm_gov, vf_exc, src_i, dst_i,
     f_t, p, global_pf_param,
     node_pf_param) = p_agg

    (node_idx_and_incident_edges_other_node_idx,
     node_inc_edges_Ybr,
     node_inc_edges_orient,
     nodes_u_view)  = node_pf_param
          
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt   = im * Y_n

    M =  (2*H) / ωs        
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91

    # i_r, i_i, vh, θh, ph, qh = f_t
    
    i  = dynamic_nodal_current_balance(src_i, dst_i) +
        uh * y_shunt
   
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[3] = dδ =  Ωb * (ω - ωs)
    
    # dx[4] = (τm_gov - τe - D * (ω - ωs)) *  (Ωb / (2*H))

    dx[1] = dδ =  (ω - ωs)
    
    # dx[2] = ( τm_gov - τe - D * ( ω - ωs ) ) *  ( 1 / (2*H) )

    dx[2] = ( τm_gov - (ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H) ) 
    
    dx[3] = ded_dash = ( -ed_dash + (X_q - X_q_dash) * iq ) / T_q_dash
    
    dx[4] = deq_dash =( -eq_dash - (X_d - X_d_dash) * id + vf_exc ) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) - ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ)  - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ)  + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
     
    return nothing
end



function node_pf_SM_2axis_idq!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
          
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt = im * Y_n

    M = (2*H) / ωs 
    
    # States
    
    # u_r  = x[1],  u_i     = x[2], δ       = x[3], 
    # ω    = x[4],  ed_dash = x[5], eq_dash = x[6]

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im

    # # ---------------------------------------------------
   
    # ------------------------------------------------------
    
    my_node_idx =
        node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j =
        node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj
                       for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
        
    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt

    i_mag = abs(i)

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)        

    # ph = (ed_dash - (ra * id - X_q_dash * iq)) * id + (eq_dash - (X_d_dash * id  + ra * iq)) * iq
    
    # qh  = (eq_dash - (X_d_dash * id  + ra * iq)) * id - (ed_dash - (ra * id - X_q_dash * iq) ) * iq    
    

    # Vk = (ph + im * qh - vh^2 * Ykk )/ conj(sum_ykj_vj)

    # Vk = vh * exp(im * angle( Vk ))
    
    # ------------------------------------------------------
    
    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq    

    dx[1] = dδ =  (ω - ωs)
    
    # dx[2] = (τm_gov - τe - D * (ω - ωs)) *  (1 / (2*H))

    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * ( ω - ωs)) / M  # (1 / (2*H)) 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
    
    # dx[5] = real(Vk ) - u_r
    
    # dx[6] = imag(Vk ) - u_i
     
    return nothing
end



function global_pf_SM_2axis_idq!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param = p_agg
    pf_U, Inet = global_pf_param

          
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt = im * Y_n
       
    M =  (2*H) / ωs 
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # i_r, i_i  = f_t 
    
    i  = Inet + uh * y_shunt       

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    

    U_r = pf_U[1]

    U_i = pf_U[2]

    # ------------------------------------------------------
    # ------------------------------------------------------
    
    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq    

    dx[1] = dδ =  (ω - ωs)
    
    # dx[2] = (τm_gov - τe - D * (ω - ωs)) *  (1 / (2*H))

    dx[2] = (τm_gov - (ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M  # (1 / (2*H)) 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
    
    # dx[5] = U_r - u_r
    
    # dx[6] = U_i - u_i
     
    return nothing
end


function network_current_SM_2axis_idq!(dx, x, p_agg, t)

    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg

          
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt   = im * Y_n

    M =  (2*H) / ωs        
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91

    # i_r, i_i, vh, θh, ph, qh = f_t
    
    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt
   
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[3] = dδ =  Ωb * (ω - ωs)
    
    # dx[4] = (τm_gov - τe - D * (ω - ωs)) *  (Ωb / (2*H))

    dx[1] = dδ =  (ω - ωs)
    
    # dx[2] = ( τm_gov - τe - D * ( ω - ωs ) ) *  ( 1 / (2*H) )

    dx[2] = ( τm_gov - (ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H) ) 
    
    dx[3] = ded_dash = ( -ed_dash + (X_q - X_q_dash) * iq ) / T_q_dash
    
    dx[4] = deq_dash =( -eq_dash - (X_d - X_d_dash) * id + vf_exc ) / T_d_dash
    
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)

    dx[5] = ed_dash * sin(δ) - ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ)  - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ)  + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
     
    return nothing
end



function initial_pf_SM_2axis_idq!(dx, x, p_agg, t)

    #u_gov, u_exc, f_t, p = p_agg
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg

          
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt   = im * Y_n

    M =  (2*H) / ωs        
    
    # States
    
    # u_r  = x[1],  u_i     = x[2], δ       = x[3], 
    # ω    = x[4],  ed_dash = x[5], eq_dash = x[6]

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91

    # i_r, i_i, vh, θh, ph, qh = f_t
    
    i_r, i_i  = f_t
    
    i  = i_r + im * i_i + uh * y_shunt    

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[3] = dδ =  Ωb * (ω - ωs)
    
    # dx[4] = (τm_gov - τe - D * (ω - ωs)) *  (Ωb / (2*H))

    dx[1] = dδ =  (ω - ωs)
    
    # dx[2] = ( τm_gov - τe - D * ( ω - ωs ) ) *  ( 1 / (2*H) )

    dx[2] = ( τm_gov - (ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H) ) 
    
    dx[3] = ded_dash = ( -ed_dash + (X_q - X_q_dash) * iq ) / T_q_dash
    
    dx[4] = deq_dash =( -eq_dash - (X_d - X_d_dash) * id + vf_exc ) / T_d_dash
    
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)

    dx[5] = ed_dash * sin(δ) - ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ)  - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ)  + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
     
    return nothing
end


@kwdef struct SM_2axis_idq <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64  = 0.0         
    D::Float64         
    H::Float64         
    Ωb::Float64        
    ωs::Float64        
    ra::Float64        
    xℓ::Float64        
    X_d::Float64       
    X_q::Float64       
    X_d_dash::Float64  
    X_q_dash::Float64  
    X_d_2dash::Float64 
    X_q_2dash::Float64 
    T_d_dash::Float64  
    T_q_dash::Float64  
    T_d_2dash::Float64 
    T_q_2dash::Float64 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0
    Q::Float64 = 0.0

    Sn::Float64 = 0.0

    vh::Float64  = 1.0

    θh::Float64  = 0.0
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_idq
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator

    state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash ]

    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))
    
    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
        
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :vh, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]

    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
        
    param_values::Vector{Float64}  = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, vh, Sn]

    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ] 

    func::Vector{Function} = Function[ initial_pf_SM_2axis_idq!, network_current_SM_2axis_idq!, global_pf_SM_2axis_idq!, node_pf_SM_2axis_idq!, hybrid_pf_SM_2axis_idq! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
        
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf

    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
     
    
end



function hybrid_pf_SM_2axis_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
          
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt   = im * Y_n           

    M         =  ( 2 * H ) / ωs
    
    # States
        
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations

    zdq  = Z_dq(ra, X_d_dash, X_q_dash)

    # voltage
    
    uh  = u_r + u_i * 1im
    
    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    dx[1] = dδ =  ( ω - ωs )
    
    dx[2] = (τm_gov - (ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 /(2*H))  
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ([ sin(δ) cos(δ); -cos(δ) sin(δ)] * ( [ed_dash, eq_dash] - zdq * [id, iq]) - [u_r, u_i] )[1]

    dx[6] = ([ sin(δ) cos(δ); -cos(δ) sin(δ)] * ( [ed_dash, eq_dash] - zdq * [id, iq]) - [u_r, u_i] )[2]
     
    return nothing
end



function node_pf_SM_2axis_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
          
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]    

    y_shunt   = im * Y_n

    M =  (2*H) / ωs

    zdq  = Z_dq(ra, X_d_dash, X_q_dash)
    
    # States
        
    δ, ω, ed_dash, eq_dash, u_r, u_i = x

    uh  = u_r + u_i * 1im
            
    θh = angle(uh)
    
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
    
    i =  Ykk * Uk + sum_ykj_vj  + Uk * y_shunt
    
    i_mag = abs(i)

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)        

    # ph = (ed_dash - (ra * id - X_q_dash * iq)) * id + (eq_dash - (X_d_dash * id  + ra * iq)) * iq
    
    # qh  = (eq_dash - (X_d_dash * id  + ra * iq)) * id - (ed_dash - (ra * id - X_q_dash * iq) ) * iq    
    

    # Vk = (ph + im * qh - vh^2 * Ykk )/ conj(sum_ykj_vj)

    # Vk = vh * exp(im * angle( Vk ))

    # ------------------------------------------------------

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    dx[1] = dδ =  (ω - ωs)
    
    # dx[2] = (τm_gov - τe - D * (ω - ωs)) *  (1 / (2*H))

     dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))  
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i       

    # dx[5] = real(Vk ) - u_r
    
    # dx[6] = imag(Vk ) - u_i   
    
     
    return nothing
end



function global_pf_SM_2axis_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param = p_agg

    pf_U, Inet = global_pf_param
          
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]    

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
        
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    zdq  = Z_dq(ra, X_d_dash, X_q_dash)
        
    uh  = u_r + u_i * 1im
    
    θh = angle(uh)
    
    # ------------------------------------------------------

    U_r = pf_U[1]

    U_i = pf_U[2]
    
    # i  = Inet[1] + uh * y_shunt

    i  = Inet + uh * y_shunt
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    
    
    # ------------------------------------------------------

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    dx[1] = dδ =  (ω - ωs)
    
    # dx[2] = (τm_gov - τe - D * (ω - ωs)) *  (1 / (2*H))

    # dx[2] = (P - τe - D * (ω - ωs)) *  (1 / (2*H))

     dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))  
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i      

    # dx[5] = U_r - u_r
    
    # dx[6] = U_i - u_i 
     
    return nothing
end


function network_current_SM_2axis_v6!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
          
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt   = im * Y_n           

    M =  ( 2 * H ) / ωs
    
    # States
        
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations

    """
    zdq     = Z_dq(ra, X_d_dash, X_q_dash)

    inv_zdq = invZ_dq(ra, X_d_dash, X_q_dash)

    [u_r, u_i] = [ sin(δ) cos(δ); -cos(δ) sin(δ)] * ( [ed_dash, eq_dash] - zdq * [id, iq])


    du_r = dx[5] = ([ sin(δ) cos(δ); -cos(δ) sin(δ)] * ( [ed_dash, eq_dash] - zdq * [id, iq]) - [u_r, u_i] )[1]

    du_i = dx[6] = ([ sin(δ) cos(δ); -cos(δ) sin(δ)] * ( [ed_dash, eq_dash] - zdq * [id, iq]) - [u_r, u_i] )[2]

    """    

    zdq  = Z_dq(ra, X_d_dash, X_q_dash)
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91   

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    # i  = dynamic_nodal_current_balance(src_i, dst_i) - uh * y_shunt   

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # # linear algebraic equ
    # ed_dash - vd = ra * id  - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    dx[1] = dδ =  ( ω - ωs )
    
    # dx[2] = (τm_gov - τe - D * (ω - ωs)) *  (1 / (2*H))

    dx[2] = (τm_gov - (ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 /(2*H))  
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    # dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    # dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i

    dx[5] = ([ sin(δ) cos(δ); -cos(δ) sin(δ)] * ( [ed_dash, eq_dash] - zdq * [id, iq]) - [u_r, u_i] )[1]

    dx[6] = ([ sin(δ) cos(δ); -cos(δ) sin(δ)] * ( [ed_dash, eq_dash] - zdq * [id, iq]) - [u_r, u_i] )[2]

     
    return nothing
end


function initial_pf_SM_2axis_v6!(dx, x, p_agg, t)

    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg

          
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs

    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)   
      
    # Network interface

    # voltage
    
    uh  = u_r + u_i * 1im
    
    # vh  = abs(uh)
    
    # θh  = angle(uh)

    # vd  = vh * sin(δ - θh)
    # vq  = vh * cos(δ - θh)
    # vdq = vd + im * vq

    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
    
    # current

    # # ALternative method  to determine id and iq
    # Sauer, eq: 3.91

    # i_r, i_i, vh, θh, ph, qh = f_t
    
    i_r, i_i  = f_t
    
    i  = i_r + im * i_i + uh * y_shunt    

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    


    # # linear algebraic equ
    # ed_dash - vd = ra * id         - X_q_dash * iq
    # eq_dash - vq = X_d_dash * id  + ra * iq

    # # solve linear algebraic equ to determine id and iq
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash + vd),
    #         (eq_dash + vq)]
    
    # Idq  = A_dq \ b_dq
    # id   = Idq[1]
    # iq   = Idq[2]

    # ------------------------------------------------------

    # τe   = vd * id + vq  * iq + ra * id * id +  ra * iq * iq
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)
    
    # ------------------------------------------------------

    dx[1] =  (ω - ωs)
    
    # dx[2] = (τm_gov - τe - D * (ω - ωs)) *  (1 / (2*H))

     dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))  
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
     
    return nothing
    
end



@kwdef struct SM_2axis_v6 <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64  = 0.0         
    D::Float64  = 0.0        
    H::Float64  = 0.0        
    Ωb::Float64 = 0.0        
    ωs::Float64 = 0.0        
    ra::Float64 = 0.0        
    xℓ::Float64  = 0.0        
    X_d::Float64 = 0.0       
    X_q::Float64 = 0.0       
    X_d_dash::Float64  = 0.0 
    X_q_dash::Float64 = 0.0  
    X_d_2dash::Float64 = 0.0 
    X_q_2dash::Float64 = 0.0 
    T_d_dash::Float64  = 0.0 
    T_q_dash::Float64 = 0.0  
    T_d_2dash::Float64 = 0.0 
    T_q_2dash::Float64 = 0.0 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0
    Q::Float64 = 0.0

    Sn::Float64 = 0.0

    vh::Float64  = 1.0

    θh::Float64  = 0.0

    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_v6

    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator

    state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i ]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))
    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :vh, :Sn]
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    

    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64}  = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, vh, Sn]

    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ]
    
    func::Vector{Function} = Function[ initial_pf_SM_2axis_v6!, network_current_SM_2axis_v6!, global_pf_SM_2axis_v6!, node_pf_SM_2axis_v6!, hybrid_pf_SM_2axis_v6! ]
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
    
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end


#-----------------------------------------------------



function hybrid_pf_SM_2axis_cb!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param

      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt   = im * Y_n

    M         =  ( 2 * H ) / ωs
    
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm = x    
  
    uh  = u_r + u_i * 1im

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt     

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    
    

    # State equations
    
    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # / (2*H)
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash    

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 
     
    return nothing
end


function node_pf_SM_2axis_cb!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param

      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M =  ( 2 * H ) / ωs

    y_shunt = im * Y_n
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm = x    
  
    uh  = u_r + u_i * 1im
   
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ])

    # ------------------------------------------------------
    
    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    
    

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq 

    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M #  / (2*H)
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash    

    # du   = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 

    # voltage
        
    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)  
    
    
    # vh  = abs(uh)
    # θh  = angle(uh)
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # Algebraics equations
    # 0 = gi(xi, xe, yi , ye, v, θ, η)
    
    # vd = vh * sin(δ - θh)  
    # vq = vh * cos(δ - θh)
    
    # vd  = ed_dash - ra * id + X_q_dash * iq     
    # vq  = eq_dash - ra * iq - X_d_dash * id
          
    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ = dδ  = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
     
    return nothing
end


function global_pf_SM_2axis_cb!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param = p_agg

    pf_U, Inet = global_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M =  ( 2 * H ) / ωs

    y_shunt = im * Y_n
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm = x    
  
    uh  = u_r + u_i * 1im

    i  = Inet + uh * y_shunt     

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    
    

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq 

    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M #  / (2*H)
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash    

    # du   = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 

    # voltage
        
    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)  
    
    
    # vh  = abs(uh)
    # θh  = angle(uh)
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # Algebraics equations
    # 0 = gi(xi, xe, yi , ye, v, θ, η)
    
    # vd = vh * sin(δ - θh)  
    # vq = vh * cos(δ - θh)
    
    # vd  = ed_dash - ra * id + X_q_dash * iq     
    # vq  = eq_dash - ra * iq - X_d_dash * id
          
    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ = dδ  = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
     
    return nothing
end


function network_current_SM_2axis_cb!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M =  ( 2 * H ) / ωs

    y_shunt = im * Y_n
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm = x    
  
    uh  = u_r + u_i * 1im

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt     

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    
    

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq 

    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # / (2*H)
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash    

    # du   = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 

    # voltage
        
    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)  
    
    
    # vh  = abs(uh)
    # θh  = angle(uh)
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # Algebraics equations
    # 0 = gi(xi, xe, yi , ye, v, θ, η)
    
    # vd = vh * sin(δ - θh)  
    # vq = vh * cos(δ - θh)
    
    # vd  = ed_dash - ra * id + X_q_dash * iq     
    # vq  = eq_dash - ra * iq - X_d_dash * id
          
    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ = dδ  = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
     
    return nothing
end


function initial_pf_SM_2axis_cb!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M =  ( 2 * H ) / ωs

    y_shunt = im * Y_n
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm = x    
  
    uh  = u_r + u_i * 1im

    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt   

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)    
    

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq 

    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # / (2*H)
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash    

    # du   = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 

    # voltage
        
    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)  
    
    
    # vh  = abs(uh)
    # θh  = angle(uh)
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # Algebraics equations
    # 0 = gi(xi, xe, yi , ye, v, θ, η)
    
    # vd = vh * sin(δ - θh)  
    # vq = vh * cos(δ - θh)
    
    # vd  = ed_dash - ra * id + X_q_dash * iq     
    # vq  = eq_dash - ra * iq - X_d_dash * id
          
    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ = dδ  = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
     
    return nothing
end


@kwdef struct SM_2axis_cb <: SdGen

    Bus::String      = "bus0"
    name::String     = lowercase(Bus)    
    P::Float64 = 0.0          
    D::Float64         
    H::Float64         
    Ωb::Float64        
    ωs::Float64        
    ra::Float64        
    xℓ::Float64        
    X_d::Float64       
    X_q::Float64       
    X_d_dash::Float64  
    X_q_dash::Float64  
    X_d_2dash::Float64 
    X_q_2dash::Float64 
    T_d_dash::Float64  
    T_q_dash::Float64  
    T_d_2dash::Float64 
    T_q_2dash::Float64 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0

    Q::Float64 = 0.0

    Sn::Float64 = 0.0

    vh::Float64 = 1.0

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_cb
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator
    
    state_vars_syms::Vector{Symbol}     = Symbol[  :δ, :ω, :ed_dash, :eq_dash]
    algebraic_vars_syms::Vector{Symbol} = Symbol[:u_r, :u_i, :vf, :τm ]
    
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    param::Vector{Symbol} = Symbol[:P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q, :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash, :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash, :αp, :αq, :Y_n, :Q, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64} =Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, Sn]
    
    func::Vector{Function} = Function[ initial_pf_SM_2axis_cb!, network_current_SM_2axis_cb!, global_pf_SM_2axis_cb!, node_pf_SM_2axis_cb!, hybrid_pf_SM_2axis_cb! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64}     = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol}  = Symbol[:u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr    
        
end



#-----------------------------------------------------


function hybrid_pf_SM_2axis_wt_loc_load_cb_idq!(dx,x,p_agg,t)
    
    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, global_pf_param,  node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    i  = 1.0 * dynamic_nodal_current_balance(src_i, dst_i) + (uh * y_shunt + loc_i_r + loc_i_i * im)

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # State equations
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
     
    return nothing
end


function node_pf_SM_2axis_wt_loc_load_cb_idq!(dx,x,p_agg,t)
    
    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]
    Q         = p[21]

    vh        = p[22] 

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs

    zdq  = Z_dq(ra, X_d_dash, X_q_dash)

    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    θh = angle(uh)
    
    # current    

    
    loc_i_r, loc_i_i = u_loc_load_ir_ii

    # ---------------------------------------------------

    # # ---------------------------------------------------

    # my_node_idx = node_idx_and_inc_edges[1]
    
    # inc_edges_orientation = last.(node_inc_edges_Ybr_orient)

    # # Get the node number of the other end of an edge
    # # that is connected to me.
    
    # nodes_j = [my_node_idx == orient[1] ? orient[2] : orient[1]   for orient in inc_edges_orientation ]

    # edges_Ybr = first.(node_inc_edges_Ybr_orient)

    # # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length(edges_Ybr )]

    # Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( inc_edges_orientation , edges_Ybr) ] ) : my_node_idx == first( inc_edges_orientation[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    # Ykk = Ybus_node[1]
    
    # Ykj = Ybus_node[2:end]
    
    # Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    # sum_ykj_vj = sum([ykj * vj for (ykj,vj) in zip(Ykj,Uj)])

    # # ---------------------------------------------------
  
   
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
      
    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt + loc_i_r + loc_i_i * im
    
    # ---------------------------------------------------

    # Algebraics equations
    # Network interface
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs)) * (Ωb / (2*H))
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - (ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq) - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    # du  = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end


function global_pf_SM_2axis_wt_loc_load_cb_idq!(dx, x, p_agg, t)
   
    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, global_pf_param = p_agg

    pf_U, Inet = global_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii
    
    i  =  Inet + uh * y_shunt + loc_i_r + loc_i_i * im 

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    
    U_r = pf_U[1]

    U_i = pf_U[2]

    # ------------------------------------------------------
    
    # τe = P
    
    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash


    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end


function network_current_SM_2axis_wt_loc_load_cb_idq!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    i  = 1.0 * dynamic_nodal_current_balance(src_i, dst_i) + (uh * y_shunt + loc_i_r + loc_i_i * im)

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # τe = P
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # State equations
    
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs)) * (Ωb / (2*H))
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)    
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M #  (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
 
     
    return nothing
end


function initial_pf_SM_2axis_wt_loc_load_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]
    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii
    
    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt + loc_i_r + loc_i_i * im

    # Algebraics equations
    # Network interface
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    # τe = P 
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs)) * (Ωb / (2*H))
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - (ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq) - D * (ω - ωs)) / M  # * (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    
     
    return nothing
end


@kwdef struct SM_2axis_wt_loc_load_cb_idq <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64 = 0.0          
    D::Float64  = 0.0        
    H::Float64  = 0.0        
    Ωb::Float64 = 0.0        
    ωs::Float64 = 0.0        
    ra::Float64 = 0.0         
    xℓ::Float64 = 0.0        
    X_d::Float64 = 0.0       
    X_q::Float64  = 0.0      
    X_d_dash::Float64 = 0.0  
    X_q_dash::Float64 = 0.0  
    X_d_2dash::Float64 = 0.0 
    X_q_2dash::Float64 = 0.0 
    T_d_dash::Float64  = 0.0 
    T_q_dash::Float64  = 0.0 
    T_d_2dash::Float64 = 0.0 
    T_q_2dash::Float64 = 0.0 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0

    Q::Float64     = 0.0

    Sn::Float64 = 0.0

    vh::Float64  = 1.0

    θh::Float64  = 0.0
    
    
    vmax::Float64  = 1.06
    vmin::Float64  = 0.97    
    Pmax::Float64  = 0.0
    Pmin::Float64  = 0.0
    Qmax::Float64  = 0.0
    Qmin::Float64  = 0.0

    comp_type::Symbol = :SM_2axis_wt_loc_load_cb_idq

    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator

    # state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i]
    # algebraic_vars_syms::Vector{Symbol} = Symbol[ ]

    state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i ]    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :vh, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64}  = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, vh, Sn]

    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ]
    
    func::Vector{Function} = Function[ initial_pf_SM_2axis_wt_loc_load_cb_idq!, network_current_SM_2axis_wt_loc_load_cb_idq!, global_pf_SM_2axis_wt_loc_load_cb_idq!, node_pf_SM_2axis_wt_loc_load_cb_idq!, hybrid_pf_SM_2axis_wt_loc_load_cb_idq! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

        
    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end

#--------------------------------------------------------


function hybrid_pf_SC_2axis_wt_loc_load_cb_idq!(dx, x, p_agg, t)
   
    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt   = im * Y_n

    M         =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    i  = 1.0 * dynamic_nodal_current_balance(src_i, dst_i) + ( uh * y_shunt + loc_i_r + loc_i_i * im)

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
     
    return nothing
end



function node_pf_SC_2axis_wt_loc_load_cb_idq!(dx, x, p_agg, t)
   
    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    vh        = p[22]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    # # --------------------------------------------------

    # my_node_idx = node_idx_and_inc_edges[1]
    
    # inc_edges_orientation = last.(node_inc_edges_Ybr_orient)

    # # Get the node number of the other end of an edge
    # # that is connected to me.
    
    # nodes_j = [my_node_idx == orient[1] ? orient[2] : orient[1]   for orient in inc_edges_orientation ]

    # edges_Ybr = first.(node_inc_edges_Ybr_orient)

    # # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length(edges_Ybr )]
    
    # Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( inc_edges_orientation , edges_Ybr) ] ) : my_node_idx == first( inc_edges_orientation[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    # Ykk = Ybus_node[1]
    
    # Ykj = Ybus_node[2:end]
    
    # Uj  = u_from_ur_ui.(nodes_u_view[nodes_j])

    # sum_ykj_vj = sum([ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # # --------------------------------------------------


   
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
        
    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt + loc_i_r + loc_i_i * im

    i_mag = abs(i)

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)        
    
    # qh  = (eq_dash - (X_d_dash * id  + ra * iq)) * id - (ed_dash - (ra * id - X_q_dash * iq) ) * iq    
    

    #  Vk = (im * qh - vh^2 * Ykk )/ sum_ykj_vj 

    # ------------------------------------------------------
        
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    

    # dx[5] = real(Vk ) - u_r
    
    # dx[6] = imag(Vk ) - u_i 
     
    return nothing
end



function global_pf_SC_2axis_wt_loc_load_cb_idq!(dx, x, p_agg, t)
    
    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p, global_pf_param = p_agg


    pf_U, Inet = global_pf_param

      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii

    i  = Inet + uh * y_shunt + loc_i_r + loc_i_i * im

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)

    # ------------------------------------------------------
    
    U_r = pf_U[1]

    U_i = pf_U[2]    

    # ------------------------------------------------------

    # τe = P
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # State equations
    
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs)) * (Ωb / (2*H))
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M # (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i    

    # dx[5] = U_r - u_r
    
    # dx[6] = U_i - u_i  
    
         
    return nothing
end



function network_current_SC_2axis_wt_loc_load_cb_idq!(dx, x, p_agg, t)

    #u_gov, u_exc, f_t, p = p_agg
    
    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ         = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

    # u_r, u_i, δ, ω, ed_dash, eq_dash = x
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii
    
    # i_r, i_i = f_t
    
    # i = i_r + im * i_i + uh * y_shunt +
    #     loc_i_r + loc_i_i * im

    i  = 1.0 * dynamic_nodal_current_balance(src_i, dst_i) + ( uh * y_shunt + loc_i_r + loc_i_i * im)

    # Algebraics equations
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------

    # τe = P
    
    # τe   = (vd + ra * id) * id + (vq + ra * iq) * iq

    # State equations
    
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs)) * (Ωb / (2*H))
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M #  (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ  # real(du)
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ  # imag(du)

    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
     
    return nothing
end



function initial_pf_SC_2axis_wt_loc_load_cb_idq!(dx, x, p_agg, t)

    #u_gov, u_exc, f_t, p = p_agg
    
    cb_sw, vf_exc, u_loc_load_ir_ii, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]
    Q         = p[21]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States
    
    δ, ω, ed_dash, eq_dash,  u_r, u_i = x
    
    # voltage
    
    uh  = u_r + u_i * 1im
    
    # current

    loc_i_r, loc_i_i = u_loc_load_ir_ii
    
    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt + loc_i_r + loc_i_i * im

    # Algebraics equations
    # Network interface
    # 0 = gi(xi, xe, yi , ye, v, θ, η) 
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq) 

    # ------------------------------------------------------
    

    dx[1] =  (ω - ωs)
    
    dx[2] =  - D * (ω - ωs) / M # (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i 
    
     
    return nothing
end


@kwdef struct SC_2axis_wt_loc_load_cb_idq <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64 = 0.0          
    D::Float64  = 0.0       
    H::Float64  = 0.0        
    Ωb::Float64 = 0.0        
    ωs::Float64 = 0.0        
    ra::Float64 = 0.0        
    xℓ::Float64  = 0.0       
    X_d::Float64  = 0.0      
    X_q::Float64  = 0.0      
    X_d_dash::Float64 = 0.0  
    X_q_dash::Float64  = 0.0 
    X_d_2dash::Float64 = 0.0 
    X_q_2dash::Float64 = 0.0 
    T_d_dash::Float64  = 0.0 
    T_q_dash::Float64  = 0.0 
    T_d_2dash::Float64 = 0.0 
    T_q_2dash::Float64 = 0.0 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0

    Q::Float64     = 0.0

    Sn::Float64 = 0.0
    
    vh::Float64    = 1.0

    θh::Float64  = 0.0
    
    
    vmax::Float64  = 1.06
    vmin::Float64  = 0.97    
    Pmax::Float64  = 0.0
    Pmin::Float64  = 0.0
    Qmax::Float64  = 0.0
    Qmin::Float64  = 0.0

    comp_type::Symbol = :SC_2axis_wt_loc_load_cb_idq

    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator

    # state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i]
    # algebraic_vars_syms::Vector{Symbol} = Symbol[ ]

    state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[  :u_r, :u_i  ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :vh, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]

    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64}  = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, vh, Sn]

    param_matrix_values::Vector{Float64}  = Float64[D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash ]
    
    func::Vector{Function} = Function[ initial_pf_SC_2axis_wt_loc_load_cb_idq! , network_current_SC_2axis_wt_loc_load_cb_idq!, global_pf_SC_2axis_wt_loc_load_cb_idq!, node_pf_SC_2axis_wt_loc_load_cb_idq!, hybrid_pf_SC_2axis_wt_loc_load_cb_idq! ]
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[
        :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    Ax::Function         = Ax_SC_gen
    Bx::Function         = Bx_SC_gen
    Cx::Function         = Cx_SC_gen
    Ax_τm_vf::Function   = Ax_SC_gen_τm_vf    

    Ax_gen_avr::Function = Ac_SC_gen_avr
    
    Ax_gen_S_avr_S::Function = SC_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SC_Ax_gen_S_avr_A
    
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end

#--------------------------------------------------------


function hybrid_pf_SM_2axis_cb_direct!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M         =  ( 2 * H ) / ωs
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm  = x   

    y_shunt = im * Y_n

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt  
            
    idq = i * exp(-im * (δ - π/2))
    id  = real(idq)
    iq  = imag(idq)    

    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M 
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 
     
    return nothing
end


function node_pf_SM_2axis_cb_direct!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M =  ( 2 * H ) / ωs
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm  = x   

    y_shunt = im * Y_n

   # # ------------------------------------------------------
   #  my_node_idx = node_idx_and_inc_edges[1]
    
   #  inc_edges_orientation = last.(node_inc_edges_Ybr_orient)

   #  # Get the node number of the other end of an edge
   #  # that is connected to me.
    
   #  nodes_j = [my_node_idx == orient[1] ? orient[2] : orient[1]   for orient in inc_edges_orientation ]

   #  edges_Ybr = first.(node_inc_edges_Ybr_orient)

   #  # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length(edges_Ybr )]
    
   #  Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( inc_edges_orientation , edges_Ybr) ] ) : my_node_idx == first( inc_edges_orientation[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

   #  Ykk = Ybus_node[1]
   #  Ykj = Ybus_node[2:end]
   #  Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

   #  sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

   #  # ------------------------------------------------------


   
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
        
    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt
     
    idq = i * exp(-im * (δ - π/2))
    id  = real(idq)
    iq  = imag(idq)    

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 

   # Network interface

    # uh   = u_r + im *  u_i
    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
            
    # i_r, i_i = f_t
    
    # i    = i_r + im * i_i
  
    
    # sh   = αp * ph + αq * im * qh
    
    # uh   = u_r + u_i * im
    # vh   = abs(uh)
    # θh   = angle(uh)
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 
          
    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ 
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ
    
    # Algebraics equations
    # 0 = gi(xi, xe, yi , ye, v, θ, η)

    # did = vd + ra * id - ed_dash - X_q_dash * iq
    # diq = vq + ra * iq - eq_dash + X_d_dash * id
    
    # dvd = vh * sin(δ - θh)  - vd
    # dvq = vh * cos(δ - θh)  - vq
    # dτe = (vd + ra * id) * id + (vq + ra * iq) * iq - τe
        
    # auxiliaries page 328, eq 15.7 and 15.8
        
    # dvf =  vf_exc - vf
    # dτm =  τm_gov - τm
     
    return nothing
end



function global_pf_SM_2axis_cb_direct!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param = p_agg

    pf_U, Inet = global_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M =  ( 2 * H ) / ωs
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm  = x   

    y_shunt = im * Y_n

    i  = Inet + uh * y_shunt  
        
    # i   = f_t[1] + 1im * f_t[2]
    
    idq = i * exp(-im * (δ - π/2))
    id  = real(idq)
    iq  = imag(idq)    

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 

   # Network interface

    # uh   = u_r + im *  u_i
    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
            
    # i_r, i_i = f_t
    
    # i    = i_r + im * i_i
  
    
    # sh   = αp * ph + αq * im * qh
    
    # uh   = u_r + u_i * im
    # vh   = abs(uh)
    # θh   = angle(uh)
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 
          
    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ 
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ
    
    # Algebraics equations
    # 0 = gi(xi, xe, yi , ye, v, θ, η)

    # did = vd + ra * id - ed_dash - X_q_dash * iq
    # diq = vq + ra * iq - eq_dash + X_d_dash * id
    
    # dvd = vh * sin(δ - θh)  - vd
    # dvq = vh * cos(δ - θh)  - vq
    # dτe = (vd + ra * id) * id + (vq + ra * iq) * iq - τe
        
    # auxiliaries page 328, eq 15.7 and 15.8
        
    # dvf =  vf_exc - vf
    # dτm =  τm_gov - τm
     
    return nothing
end



function network_current_SM_2axis_cb_direct!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M =  ( 2 * H ) / ωs
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm  = x   

    y_shunt = im * Y_n

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt  
        
    # i   = f_t[1] + 1im * f_t[2]
    
    idq = i * exp(-im * (δ - π/2))
    id  = real(idq)
    iq  = imag(idq)    

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M #  (1 / (2*H))
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 

   # Network interface

    # uh   = u_r + im *  u_i
    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
            
    # i_r, i_i = f_t
    
    # i    = i_r + im * i_i
  
    
    # sh   = αp * ph + αq * im * qh
    
    # uh   = u_r + u_i * im
    # vh   = abs(uh)
    # θh   = angle(uh)
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 
          
    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ 
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ
    
    # Algebraics equations
    # 0 = gi(xi, xe, yi , ye, v, θ, η)

    # did = vd + ra * id - ed_dash - X_q_dash * iq
    # diq = vq + ra * iq - eq_dash + X_d_dash * id
    
    # dvd = vh * sin(δ - θh)  - vd
    # dvq = vh * cos(δ - θh)  - vq
    # dτe = (vd + ra * id) * id + (vq + ra * iq) * iq - τe
        
    # auxiliaries page 328, eq 15.7 and 15.8
        
    # dvf =  vf_exc - vf
    # dτm =  τm_gov - τm
     
    return nothing
end




function initial_pf_SM_2axis_cb_direct!(dx, x, p_agg, t)

    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M =  ( 2 * H ) / ωs
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm  = x   

    y_shunt = im * Y_n

    i_r, i_i = f_t
    
    i = i_r + im * i_i + uh * y_shunt 
        
    # i   = f_t[1] + 1im * f_t[2]
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    iq  = imag(idq)    

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq

    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M #  (1 / (2*H))
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm 

   # Network interface

    # uh   = u_r + im *  u_i
    # vdq = uh * exp(-im * (δ - π/2))
    # vd  = real(vdq)
    # vq  = imag(vdq)
            
    # i_r, i_i = f_t
    
    # i    = i_r + im * i_i
  
    
    # sh   = αp * ph + αq * im * qh
    
    # uh   = u_r + u_i * im
    # vh   = abs(uh)
    # θh   = angle(uh)
    
    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq 
          
    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    # du    = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ 
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ
    
    # Algebraics equations
    # 0 = gi(xi, xe, yi , ye, v, θ, η)

    # did = vd + ra * id - ed_dash - X_q_dash * iq
    # diq = vq + ra * iq - eq_dash + X_d_dash * id
    
    # dvd = vh * sin(δ - θh)  - vd
    # dvq = vh * cos(δ - θh)  - vq
    # dτe = (vd + ra * id) * id + (vq + ra * iq) * iq - τe
        
    # auxiliaries page 328, eq 15.7 and 15.8
        
    # dvf =  vf_exc - vf
    # dτm =  τm_gov - τm
     
    return nothing
end


@kwdef struct SM_2axis_cb_direct <: SdGen

    Bus::String      = "bus0"
    name::String     = lowercase(Bus)    
    P::Float64  = 0.0         
    D::Float64         
    H::Float64         
    Ωb::Float64        
    ωs::Float64        
    ra::Float64        
    xℓ::Float64        
    X_d::Float64       
    X_q::Float64       
    X_d_dash::Float64  
    X_q_dash::Float64  
    X_d_2dash::Float64 
    X_q_2dash::Float64 
    T_d_dash::Float64  
    T_q_dash::Float64  
    T_d_2dash::Float64 
    T_q_2dash::Float64 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0

    Q::Float64 = 0.0

    Sn::Float64 = 0.0

    vh::Float64 = 1.0

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_cb_direct
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator
    
    state_vars_syms::Vector{Symbol}  = Symbol[ :δ, :ω, :ed_dash, :eq_dash ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[:u_r, :u_i, :vf, :τm  ]

    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    param::Vector{Symbol} = Symbol[:P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q, :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash, :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash, :αp, :αq, :Y_n, :Q, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64} = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, Sn]
    func::Vector{Function} = Function[ initial_pf_SM_2axis_cb_direct!, network_current_SM_2axis_cb_direct!, global_pf_SM_2axis_cb_direct!, node_pf_SM_2axis_cb_direct!, hybrid_pf_SM_2axis_cb_direct! ]
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64}     = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol}  = Symbol[:u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A

    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end



function hybrid_pf_SM_2axis_cb_millano!(dx, x, p_agg, t)

    # 0 = vh * exp(im * θh) + (ra + im * X_d_dash) * (id + im * iq) * exp(im * (δ - π/2)) - (ed_dash + (X_q_dash - X_d_dash) * iq + im * eq_dash) * exp(im * (δ - π/2))
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt   = im * Y_n

    M         =  ( 2 * H ) / ωs
    
    # States

     δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm, τe = x    

    
    # Network interface

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)

    # State equations
    
    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq )  - D * (ω - ωs)) / M #  (1 / (2*H))
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # auxiliaries page 328, eq 15.7 and 15.8
        
    # dvf =  vf_exc - vf
    # dτm =  τm_gov - τm

    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm
    
    dx[9] = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq - τe
     
    return nothing
end


function node_pf_SM_2axis_cb_millano!(dx, x, p_agg, t)

    # 0 = vh * exp(im * θh) + (ra + im * X_d_dash) * (id + im * iq) * exp(im * (δ - π/2)) - (ed_dash + (X_q_dash - X_d_dash) * iq + im * eq_dash) * exp(im * (δ - π/2))
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

     δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm, τe = x    

    
    # Network interface

   # # ------------------------------------------------------
   #  my_node_idx = node_idx_and_inc_edges[1]
    
   #  inc_edges_orientation = last.(node_inc_edges_Ybr_orient)

   #  # Get the node number of the other end of an edge
   #  # that is connected to me.
    
   #  nodes_j = [my_node_idx == orient[1] ? orient[2] : orient[1]   for orient in inc_edges_orientation ]

   #  edges_Ybr = first.(node_inc_edges_Ybr_orient)

   #  # Ybus_node = [ k == 0 ? sum([Ybr[1][1] for Ybr in edges_Ybr]) : edges_Ybr[k][1][2] for k in 0:length(edges_Ybr )]
    
   #  Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( inc_edges_orientation , edges_Ybr) ] ) : my_node_idx == first( inc_edges_orientation[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

   #  Ykk = Ybus_node[1]
   #  Ykj = Ybus_node[2:end]
   #  Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

   #  sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

   #  # ------------------------------------------------------


   
    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ]) 

    # ------------------------------------------------------
        
    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt
    
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq    
    
    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq )  - D * (ω - ωs)) / M #  (1 / (2*H))
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ 
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # auxiliaries page 328, eq 15.7 and 15.8
        
    # dvf =  vf_exc - vf
    # dτm =  τm_gov - τm

    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm
    
    dx[9] = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq - τe

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash

    # du = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # uh   = vh * exp(im * θh)
   
    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    
    # vd = vh * sin(δ - θh) = uh_r * sin(δ)  - uh_i * cos(δ)   
    # vq = vh * cos(δ - θh) = uh_r * cos(δ)  + uh_i * sin(δ) 
    
    # vdq = vd + im * vq
    

    # uh = ( (ed_dash + (X_q_dash - X_d_dash) * iq + im * eq_dash ) - (ra + im * X_d_dash) * (id + im * iq) ) * exp(im * (δ - π/2))

    # exp(im * (δ - π/2)) = sin(δ) - im * cos(δ)

    # uh_r = ((X_q_dash * iq - ra * id) + ed_dash) * sin(δ)

    # uh_i = ((ra * id - X_q_dash * iq) - ed_dash) * cos(δ)
    
    # vd = ed_dash - (ra * id  - X_q_dash * iq)
    
    # vq = eq_dash - (X_d_dash * id  + ra * iq)    

    # vdq = uh * exp(-im * (δ - π/2))
    
    # vd  = real(vdq)
    # vq  = imag(vdq)   
    
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash - vd),
    #         (eq_dash - vq)]
    
    # Idq  = A_dq \ b_dq
    
    # id   = Idq[1]
    
    # iq   = Idq[2]

    # Sauer, eq: 3.91
    
     
    return nothing
end



function global_pf_SM_2axis_cb_millano!(dx, x, p_agg, t)

    # 0 = vh * exp(im * θh) + (ra + im * X_d_dash) * (id + im * iq) * exp(im * (δ - π/2)) - (ed_dash + (X_q_dash - X_d_dash) * iq + im * eq_dash) * exp(im * (δ - π/2))
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, global_pf_param  = p_agg

    pf_U, Inet = global_pf_param
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

     δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm, τe = x    

    
    # Network interface

    i  = Inet + uh * y_shunt

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq    
    
    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq )  - D * (ω - ωs)) / M #  (1 / (2*H))
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ 
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # auxiliaries page 328, eq 15.7 and 15.8
        
    # dvf =  vf_exc - vf
    # dτm =  τm_gov - τm

    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm
    
    dx[9] = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq - τe

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash

    # du = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # uh   = vh * exp(im * θh)
   
    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    
    # vd = vh * sin(δ - θh) = uh_r * sin(δ)  - uh_i * cos(δ)   
    # vq = vh * cos(δ - θh) = uh_r * cos(δ)  + uh_i * sin(δ) 
    
    # vdq = vd + im * vq
    

    # uh = ( (ed_dash + (X_q_dash - X_d_dash) * iq + im * eq_dash ) - (ra + im * X_d_dash) * (id + im * iq) ) * exp(im * (δ - π/2))

    # exp(im * (δ - π/2)) = sin(δ) - im * cos(δ)

    # uh_r = ((X_q_dash * iq - ra * id) + ed_dash) * sin(δ)

    # uh_i = ((ra * id - X_q_dash * iq) - ed_dash) * cos(δ)
    
    # vd = ed_dash - (ra * id  - X_q_dash * iq)
    
    # vq = eq_dash - (X_d_dash * id  + ra * iq)    

    # vdq = uh * exp(-im * (δ - π/2))
    
    # vd  = real(vdq)
    # vq  = imag(vdq)   
    
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash - vd),
    #         (eq_dash - vq)]
    
    # Idq  = A_dq \ b_dq
    
    # id   = Idq[1]
    
    # iq   = Idq[2]

    # Sauer, eq: 3.91
    
     
    return nothing
end



function network_current_SM_2axis_cb_millano!(dx, x, p_agg, t)

    # 0 = vh * exp(im * θh) + (ra + im * X_d_dash) * (id + im * iq) * exp(im * (δ - π/2)) - (ed_dash + (X_q_dash - X_d_dash) * iq + im * eq_dash) * exp(im * (δ - π/2))
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

     δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm, τe = x    

    
    # Network interface

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq    
    
    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq )  - D * (ω - ωs)) / M #  (1 / (2*H))
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ 
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # auxiliaries page 328, eq 15.7 and 15.8
        
    # dvf =  vf_exc - vf
    # dτm =  τm_gov - τm

    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm
    
    dx[9] = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq - τe

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash

    # du = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # uh   = vh * exp(im * θh)
   
    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    
    # vd = vh * sin(δ - θh) = uh_r * sin(δ)  - uh_i * cos(δ)   
    # vq = vh * cos(δ - θh) = uh_r * cos(δ)  + uh_i * sin(δ) 
    
    # vdq = vd + im * vq
    

    # uh = ( (ed_dash + (X_q_dash - X_d_dash) * iq + im * eq_dash ) - (ra + im * X_d_dash) * (id + im * iq) ) * exp(im * (δ - π/2))

    # exp(im * (δ - π/2)) = sin(δ) - im * cos(δ)

    # uh_r = ((X_q_dash * iq - ra * id) + ed_dash) * sin(δ)

    # uh_i = ((ra * id - X_q_dash * iq) - ed_dash) * cos(δ)
    
    # vd = ed_dash - (ra * id  - X_q_dash * iq)
    
    # vq = eq_dash - (X_d_dash * id  + ra * iq)    

    # vdq = uh * exp(-im * (δ - π/2))
    
    # vd  = real(vdq)
    # vq  = imag(vdq)   
    
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash - vd),
    #         (eq_dash - vq)]
    
    # Idq  = A_dq \ b_dq
    
    # id   = Idq[1]
    
    # iq   = Idq[2]

    # Sauer, eq: 3.91
    
     
    return nothing
end



function initial_pf_SM_2axis_cb_millano!(dx, x, p_agg, t)

    # 0 = vh * exp(im * θh) + (ra + im * X_d_dash) * (id + im * iq) * exp(im * (δ - π/2)) - (ed_dash + (X_q_dash - X_d_dash) * iq + im * eq_dash) * exp(im * (δ - π/2))
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    # States

     δ, ω, ed_dash, eq_dash, u_r, u_i, vf, τm, τe = x    

    
    # Network interface

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    
    iq  = imag(idq)

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq    
    
    dx[1] = dδ = (ω - ωs)
    
    dx[2] = (τm - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq )  - D * (ω - ωs)) / M #  (1 / (2*H))
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash
    
    dx[4] = deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash
    
    # dx[5] = -cos(δ) * deq_dash - sin(δ) * ded_dash - u_i * dδ
    
    # dx[6] = cos(δ) *  ded_dash - sin(δ) * deq_dash + u_r * dδ 
    
    dx[5] = ed_dash * sin(δ) -  ra * id * sin(δ) + X_q_dash * iq * sin(δ) + eq_dash * cos(δ) - X_d_dash * id * cos(δ) - ra * iq * cos(δ) - u_r

    dx[6] = eq_dash * sin(δ) - X_d_dash * id * sin(δ) - ra * iq * sin(δ) - ed_dash * cos(δ) + ra * id * cos(δ) - X_q_dash * iq * cos(δ) - u_i
    
    # auxiliaries page 328, eq 15.7 and 15.8
        
    # dvf =  vf_exc - vf
    # dτm =  τm_gov - τm

    dx[7] = vf_exc - vf
    
    dx[8] = τm_gov - τm
    
    dx[9] = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq - τe

    # State equations
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash

    # du = -im * exp(im * δ) * (ded_dash + im * deq_dash) + im * uh * dδ

    # Network interface
    
    # ph_i = gp_hi(xi, xe, yi , ye, v, θ, η)
    # qh_i = gq_hi(xi, xe, yi , ye, v, θ, η)

    # uh   = vh * exp(im * θh)
   
    # ph = vd * id + vq * iq 
    # qh = vq * id - vd * iq
    
    # vd = vh * sin(δ - θh) = uh_r * sin(δ)  - uh_i * cos(δ)   
    # vq = vh * cos(δ - θh) = uh_r * cos(δ)  + uh_i * sin(δ) 
    
    # vdq = vd + im * vq
    

    # uh = ( (ed_dash + (X_q_dash - X_d_dash) * iq + im * eq_dash ) - (ra + im * X_d_dash) * (id + im * iq) ) * exp(im * (δ - π/2))

    # exp(im * (δ - π/2)) = sin(δ) - im * cos(δ)

    # uh_r = ((X_q_dash * iq - ra * id) + ed_dash) * sin(δ)

    # uh_i = ((ra * id - X_q_dash * iq) - ed_dash) * cos(δ)
    
    # vd = ed_dash - (ra * id  - X_q_dash * iq)
    
    # vq = eq_dash - (X_d_dash * id  + ra * iq)    

    # vdq = uh * exp(-im * (δ - π/2))
    
    # vd  = real(vdq)
    # vq  = imag(vdq)   
    
    
    # Algebraics equations
    
    # 0 = gi(xi, xe, yi , ye, v, θ, η)
    
    # A_dq = [ra          -X_q_dash;
    #         X_d_dash    ra]
    
    # b_dq = [(ed_dash - vd),
    #         (eq_dash - vq)]
    
    # Idq  = A_dq \ b_dq
    
    # id   = Idq[1]
    
    # iq   = Idq[2]

    # Sauer, eq: 3.91
    
     
    return nothing
end

@kwdef struct SM_2axis_cb_millano <: SdGen

    Bus::String      = "bus0"
    name::String     = lowercase(Bus)    
    P::Float64  = 0.0         
    D::Float64         
    H::Float64         
    Ωb::Float64        
    ωs::Float64        
    ra::Float64        
    xℓ::Float64        
    X_d::Float64       
    X_q::Float64       
    X_d_dash::Float64  
    X_q_dash::Float64  
    X_d_2dash::Float64 
    X_q_2dash::Float64 
    T_d_dash::Float64  
    T_q_dash::Float64  
    T_d_2dash::Float64 
    T_q_2dash::Float64 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0

    Q::Float64 = 0.0

    Sn::Float64 = 0.0

    vh::Float64 = 1.0

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_cb_millano
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator
    
    state_vars_syms::Vector{Symbol}  = Symbol[ :δ, :ω, :ed_dash, :eq_dash]
    algebraic_vars_syms::Vector{Symbol} = Symbol[:u_r, :u_i, :vf, :τm,  :τe]

    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    param::Vector{Symbol} = Symbol[:P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q, :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash, :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash, :αp, :αq, :Y_n, :Q, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Float64} = Float64[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, Sn]
    func::Vector{Function} = Function[ initial_pf_SM_2axis_cb_millano!, network_current_SM_2axis_cb_millano!, global_pf_SM_2axis_cb_millano!, node_pf_SM_2axis_cb_millano!, hybrid_pf_SM_2axis_cb_millano! ]
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64}     = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol}  = Symbol[:u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    

    
end


function new_SM_2axis_cb_inf!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, node_idx_and_inc_edges, node_inc_edges_Ybr_orient, nodes_u_view = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    vh, θh  = f_t
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x    
       
    uh = vh * exp( im * θh )

    ubus_r = vh * cos(θh)
    
    ubus_i = vh * sin(θh)

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt
   
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    iq  = imag(idq)

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    # State equations
    
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M #  (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ubus_r - u_r 

    dx[6] = ubus_i - u_i
     
    return nothing
end



function SM_2axis_cb_inf!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt = im * Y_n

    M =  ( 2 * H ) / ωs
    
    vh, θh  = f_t
    
    # States
    
    δ, ω, ed_dash, eq_dash, u_r, u_i = x    
       
    uh = vh * exp( im * θh )

    ubus_r = vh * cos(θh)
    
    ubus_i = vh * sin(θh)

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt
   
    idq = i * exp(-im * (δ - π/2))
    
    id  = real(idq)
    iq  = imag(idq)

    # τe = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    # State equations
    
    # dxi = fi(xi, xe, yi , ye, v, θ, η)

    # dδ       = Ωb * (ω - ωs)
    # dω       = (τm - τe - D * (ω - ωs))/(2*H)
    # ded_dash = (-ed_dash + (X_q - X_q_dash) * iq )/T_q_dash 
    # deq_dash = (-eq_dash - (X_d - X_d_dash) * id + vf)/T_d_dash

    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1/(2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash
    
    dx[5] = ubus_r - u_r 

    dx[6] = ubus_i - u_i
     
    return nothing
end



@kwdef struct SM_2axis_cb_inf <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64  = 0.0         
    D::Float64         
    H::Float64         
    Ωb::Float64        
    ωs::Float64        
    ra::Float64        
    xℓ::Float64        
    X_d::Float64       
    X_q::Float64       
    X_d_dash::Float64  
    X_q_dash::Float64  
    X_d_2dash::Float64 
    X_q_2dash::Float64 
    T_d_dash::Float64  
    T_q_dash::Float64  
    T_d_2dash::Float64 
    T_q_2dash::Float64 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0

    Q::Float64 = 0.0

    Sn::Float64 = 0.0
    
    U::ComplexF64 = 1.0 + im * 0.0

    vh::Float64    = 1.06

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_cb_inf

    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator
    
    state_vars_syms::Vector{Symbol} = Symbol[ :δ, :ω, :ed_dash, :eq_dash]
    algebraic_vars_syms::Vector{Symbol} = Symbol[:u_r, :u_i ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q,  :U, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Union{Float64, ComplexF64}} = Union{Float64, ComplexF64}[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Q, Y_n, U, Sn]
    func::Vector{Function} = Function[SM_2axis_cb_inf!, SM_2axis_cb_inf!, new_SM_2axis_cb_inf!]
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[  :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end


function new_SM_2axis_cb_inf_bus!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p, node_idx_and_inc_edges, node_inc_edges_Ybr_orient, nodes_u_view = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    y_shunt   = im * Y_n

    M =  ( 2 * H ) / ωs
    
    vh, θh    = f_t
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i = x    

    uh = vh * exp( im * θh )

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt    
 

    vdq = uh * exp(-im * (δ - π/2))
    
    vd  = real(vdq)
    vq  = imag(vdq)    
    
    idq  = i * exp(-im * (δ - π/2))
    
    id   = real(idq)
    
    iq   = imag(idq)

    # τe   = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    # State equations
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))   
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = vh * cos(θh) - u_r
    
    dx[6] = vh * sin(θh) - u_i
     
     
    return nothing
end



function SM_2axis_cb_inf_bus!(dx, x, p_agg, t)
    
    cb_sw, τm_gov, vf_exc, src_i, dst_i, f_t, p = p_agg
      
    # Parameters
    
    P         = p[1]    
    D         = p[2]
    H         = p[3]
    Ωb        = p[4]
    ωs        = p[5]
    
    ra        = p[6]
    xℓ        = p[7]
    X_d       = p[8]
    X_q       = p[9]
    
    X_d_dash  = p[10]
    X_q_dash  = p[11]
    
    X_d_2dash = p[12]
    X_q_2dash = p[13]
    
    T_d_dash  = p[14]
    T_q_dash  = p[15]
    
    T_d_2dash = p[16]
    T_q_2dash = p[17]

    αp        = p[18]
    αq        = p[19]

    Y_n       = p[20]

    M         =  ( 2 * H ) / ωs

    y_shunt   = im * Y_n 
    vh, θh    = f_t
    
    # States

    δ, ω, ed_dash, eq_dash, u_r, u_i = x    

    uh = vh * exp( im * θh )

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt    
 

    vdq = uh * exp(-im * (δ - π/2))
    
    vd  = real(vdq)
    vq  = imag(vdq)    
    
    idq  = i * exp(-im * (δ - π/2))
    
    id   = real(idq)
    
    iq   = imag(idq)

    # τe   = ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq
    
    # State equations
    
    dx[1] = dδ =  (ω - ωs)
    
    dx[2] = (τm_gov - ( ed_dash * id + eq_dash * iq + (X_q_dash - X_d_dash) * id * iq ) - D * (ω - ωs)) / M # (1 / (2*H))    
    
    dx[3] = ded_dash = (-ed_dash + (X_q - X_q_dash) * iq) / T_q_dash
    
    dx[4] = deq_dash =(-eq_dash - (X_d - X_d_dash) * id + vf_exc) / T_d_dash

    dx[5] = vh * cos(θh) - u_r
    
    dx[6] = vh * sin(θh) - u_i
     
     
    return nothing
end


@kwdef struct SM_2axis_cb_inf_bus <: SdGen
    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    P::Float64   = 0.0        
    D::Float64         
    H::Float64         
    Ωb::Float64        
    ωs::Float64        
    ra::Float64        
    xℓ::Float64        
    X_d::Float64       
    X_q::Float64       
    X_d_dash::Float64  
    X_q_dash::Float64  
    X_d_2dash::Float64 
    X_q_2dash::Float64 
    T_d_dash::Float64  
    T_q_dash::Float64  
    T_d_2dash::Float64 
    T_q_2dash::Float64 
    αp::Float64  = 1.0        
    αq::Float64  = 1.0        
    Y_n::Float64 = 0.0

    Q::Float64 = 0.0

    Sn::Float64 = 0.0
    
    U::ComplexF64 = 1.0 + im * 0.0

    vh::Float64    = 1.06

    θh::Float64  = 0.0
    
    
    vmax::Float64    = 1.06
    vmin::Float64    = 0.90    
    Pmax::Float64    = 0.0
    Pmin::Float64    = 0.0
    Qmax::Float64    = 0.0
    Qmin::Float64    = 0.0

    comp_type::Symbol = :SM_2axis_cb_inf_bus
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator
    
    state_vars_syms::Vector{Symbol} = Symbol[  :δ, :ω, :ed_dash, :eq_dash]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    stab_state_vars_syms::Vector{Symbol} = Symbol[ :ω, :ed_dash, :eq_dash ]

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[
        :P, :D, :H, :Ωb, :ωs, :ra, :xℓ, :X_d, :X_q,
        :X_d_dash, :X_q_dash,:X_d_2dash, :X_q_2dash,
        :T_d_dash, :T_q_dash, :T_d_2dash,:T_q_2dash,
        :αp, :αq, :Y_n, :Q, :U, :Sn]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]


    gen_sub_dyn_param_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :D ], dict_param_syms_Idx[ :H ], dict_param_syms_Idx[ :ωs ], dict_param_syms_Idx[ :ra ],  dict_param_syms_Idx[ :xℓ ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ], dict_param_syms_Idx[ :X_d_2dash ], dict_param_syms_Idx[ :X_q_2dash ], dict_param_syms_Idx[ :T_d_dash ], dict_param_syms_Idx[ :T_q_dash ], dict_param_syms_Idx[ :T_d_2dash ], dict_param_syms_Idx[ :T_q_2dash ] ]
    
    ra_Xd_Xq_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ] ]

    ra_Xd_Xq_Xd_dash_Xq_dash_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :ra ], dict_param_syms_Idx[ :X_d ], dict_param_syms_Idx[ :X_q ], dict_param_syms_Idx[ :X_d_dash ], dict_param_syms_Idx[ :X_q_dash ] ]
    
    param_values::Vector{Union{Float64, ComplexF64}} = Union{Float64, ComplexF64}[P, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash, αp, αq, Y_n, Q, U, Sn]
    func::Vector{Function} = Function[SM_2axis_cb_inf_bus!, SM_2axis_cb_inf_bus!,  new_SM_2axis_cb_inf_bus!]
    control_sig_syms::Vector{Symbol} = Symbol[:i_r, :i_i]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[  :u_r, :u_i, :δ, :ω]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax::Function            = Ax_gen
    Bx::Function            = Bx_gen
    Cx::Function            = Cx_gen
    Ax_τm_vf::Function      = Ax_gen_τm_vf

    Ax_gen_gov::Function    = Ac_gen_gov
    Ax_gen_avr::Function    = Ac_gen_avr

    Ax_gen_S_gov_S::Function = SM_Ax_gen_S_gov_S
    Ax_gen_S_gov_A::Function = SM_Ax_gen_S_gov_A
    
    Ax_gen_S_avr_S::Function = SM_Ax_gen_S_avr_S
    Ax_gen_S_avr_A::Function = SM_Ax_gen_S_avr_A
    
    # for stability analysis
    
    stab_Ax::Function            = stab_Ax_gen
    stab_Bx::Function            = stab_Bx_gen
    stab_Cx::Function            = stab_Cx_gen
    stab_Ax_τm_vf::Function      = stab_Ax_gen_τm_vf
    stab_Ax_gen_gov::Function    = stab_Ac_gen_gov
    stab_Ax_gen_avr::Function    = stab_Ac_gen_avr
    
    
end




function new_Infinite_cb_bus!(dx, x, p_agg, t)

    cb_sw, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
    
    U     = p[1]
    P     = p[2]
    Q     = p[3]
    Y_n   = p[4]    

    vh = abs(U)
    
    θh = angle(U)

    # Iinj, pm_vf_x = in_f_t        
    # i    = Iinj[1] + im * Iinj[2]
    # vh, θh, ph, qh  = pm_vf_x

    i_r, i_i, vh, θh = f_t

    # i = i_r + im * i_i

    uh  =  U 

    y_shunt = im * Y_n

    # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ])

    # ------------------------------------------------------

    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt      
    
    # u     = x[1] + x[2] * im
    
    u_r, u_i  = x

    dx[1]   = vh * cos(θh)  - u_r
    
    dx[2]   = vh * sin(θh)  - u_r    

    # nodal_current_balance = dynamic_nodal_current_balance(x_s, x_d) + Y_n * u
    # igen_net = -1.0 * nodal_current_balance
    # s      = u * conj(igen_net)
    
    # y_t[1] = real(s) # ph
    # y_t[2] = imag(s) # qh
    # y_t[3] = real(igen_net) # islack_r
    # y_t[4] = imag(igen_net) # islack_i
    
    # du     = U - u
    
    # dx[1]  = real(du)
    # dx[2]  = imag(du)
   
    return nothing
end


function Infinite_cb_bus!(dx, x, p_agg, t)

    cb_sw, src_i, dst_i, f_t, p = p_agg
    
    U     = p[1]
    P     = p[2]
    Q     = p[3]
    Y_n   = p[4]    

    vh = abs(U)
    
    θh = angle(U)

    # Iinj, pm_vf_x = in_f_t        
    # i    = Iinj[1] + im * Iinj[2]
    # vh, θh, ph, qh  = pm_vf_x

    i_r, i_i, vh, θh = f_t

    # i = i_r + im * i_i

    uh  =  U 

    y_shunt = im * Y_n    

    i  = dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt  
    
    # u     = x[1] + x[2] * im
    
    u_r     = x[1]
    
    u_i     = x[2]

    dx[1]   = vh * cos(θh)  - u_r
    
    dx[2]   = vh * sin(θh)  - u_r    

    # nodal_current_balance = dynamic_nodal_current_balance(x_s, x_d) + Y_n * u
    # igen_net = -1.0 * nodal_current_balance
    # s      = u * conj(igen_net)
    
    # y_t[1] = real(s) # ph
    # y_t[2] = imag(s) # qh
    # y_t[3] = real(igen_net) # islack_r
    # y_t[4] = imag(igen_net) # islack_i
    
    # du     = U - u
    
    # dx[1]  = real(du)
    # dx[2]  = imag(du)
   
    return nothing
end

@kwdef struct Infinite_cb_bus  <: SdGen

    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    U::ComplexF64    = 1.0 + im * 0.0
    P::Float64       = 0.0
    Q::Float64       = 0.0
    Y_n              = 0.0

    comp_type::Symbol = :Infinite_cb_bus

    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator
    
    state_vars_syms::Vector{Symbol} = Symbol[  ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    param::Vector{Symbol} = Symbol[:U, :P, :Q, :Y_n ]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]  
    
    param_values::Vector{Union{Float64, ComplexF64}} = Union{Float64, ComplexF64}[U, P, Q, Y_n]
    func::Vector{Function} = Function[Infinite_cb_bus!, Infinite_cb_bus!, new_Infinite_cb_bus!]
    control_sig_syms::Vector{Symbol} = Symbol[]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[ :ph, :qh, :islack_r, :islack_i ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)
    
end



function new_Infinite_bus!(dx, x, p_agg, t)

    cb_sw, src_i, dst_i, f_t, p, node_pf_param = p_agg

    node_idx_and_incident_edges_other_node_idx, node_inc_edges_Ybr, node_inc_edges_orient, nodes_u_view  = node_pf_param
    
    U     = p[1]
    P     = p[2]
    Q     = p[3]
    Y_n   = p[4]   

    i_r, i_i, vh, θh = f_t

    # U       = vh * exp(im * θh)
    # i = i_r + im * i_i
    # i  = -1.0 * dynamic_nodal_current_balance(src_i, dst_i)

    y_shunt = im * Y_n

    u_r, u_i  = x    

   # ------------------------------------------------------
    
    my_node_idx = node_idx_and_incident_edges_other_node_idx[1]
    
    nodes_j = node_idx_and_incident_edges_other_node_idx[2:end]

    edges_Ybr = node_inc_edges_Ybr
    
    Ybus_node = [ k == 0 ? sum( [ my_node_idx == first( orient ) ? Ybr[1] : Ybr[4] for (orient, Ybr) in zip( node_inc_edges_orient , edges_Ybr) ] ) : my_node_idx == first( node_inc_edges_orient[k] ) ? (edges_Ybr[k])[3] : (edges_Ybr[k])[2] for k in 0:length( edges_Ybr ) ]

    Ykk = Ybus_node[1]
    
    Ykj = Ybus_node[2:end]
    
    Uj  = u_from_ur_ui.(nodes_u_view[nodes_j] )

    Uk  = u_from_ur_ui( nodes_u_view[my_node_idx] )

    sum_ykj_vj = sum([ ykj * vj for (ykj, vj) in zip(Ykj,Uj) ])

    # ------------------------------------------------------

    i =  Ykk * uh + sum_ykj_vj  + uh * y_shunt      

    
    dx[1]  = real(U) - u_r
    dx[2]  = imag(U) - u_i
   
    return nothing
end


function Infinite_bus!(dx, x, p_agg, t)

    cb_sw, src_i, dst_i, f_t, p = p_agg
    
    U     = p[1]
    P     = p[2]
    Q     = p[3]
    Y_n   = p[4]   

    i_r, i_i, vh, θh = f_t

    # U       = vh * exp(im * θh)
    # i = i_r + im * i_i
    # i  = -1.0 * dynamic_nodal_current_balance(src_i, dst_i)

    y_shunt = im * Y_n    

    i  = dynamic_nodal_current_balance(src_i, dst_i) + U * y_shunt     
    
    u_r, u_i  = x    
    
    dx[1]  = real(U) - u_r
    dx[2]  = imag(U) - u_i
   
    return nothing
end

@kwdef struct Infinite_bus <: SdGen

    Bus::String      = "bus0"
    name::String     = lowercase(Bus)     
    U::ComplexF64    = 1.0 +  im * 0.0
    P::Float64       = 0.0
    Q::Float64       = 0.0
    Y_n              = 0.0

    comp_type::Symbol = :Infinite_bus
    
    Bus_num::Int64 = parse(Int, split(lowercase(Bus),"bus")[2] )    
    Bus_type::Symbol  = :Generator
    
    state_vars_syms::Vector{Symbol} = Symbol[ ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :u_r, :u_i ]
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[]
    
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]    
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    param::Vector{Symbol} = Symbol[:U, :P, :Q, :Y_n ]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( param ) )

    P_Q_idx::Vector{Int64} = Int64[ dict_param_syms_Idx[ :P ], dict_param_syms_Idx[ :Q ]]  
    
    param_values::Vector{Union{Float64, ComplexF64}} = Union{Float64, ComplexF64}[U, P, Q, Y_n]
    func::Vector{Function} = Function[Infinite_bus!, Infinite_bus!, new_Infinite_bus!]
    control_sig_syms::Vector{Symbol} = Symbol[]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[ ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
        
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]

    cb_state_values::Vector{Function}          = Function[]    
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2_sw[sym]  for sym in cb_dyn_state_syms ]
            
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]
    cb_dyn_state_values::Vector{Function}      = Function[]    
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)
    
end

#-----------------------------------------------------
#-----------------------------------------------------
