
########################################################
########################################################
# ------------------------------------------------------
# initialisation
# ------------------------------------------------------
########################################################
########################################################


#########################################################
# ------------------------------------------------------
#  self init
# ------------------------------------------------------
#########################################################

function two_axix_synch_machine_init(p, power_flow_data)

    # ---------------------------------------
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

    # ---------------------------------------

    # power_flow_data

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

    
    # ---------------------------------------
  
    v        = vh * exp(im * θh)

    u0_r      = real(v)
    
    u0_i      = imag(v)

    # s_α      = αp * ph + αq * im * qh
    
    s_α       = αp * pg + αq * im * qg
    
    i0       = conj(s_α) / conj(v)

    ih0      = abs(i0)
    
    ϕ0       = angle(i0)    
    
    E0       = v + (ra + im * X_q) * i0
    
    δ0       = angle(E0)
    
    v_dq     = v * exp(-im * (δ0 - pi/2))
    vd0      = real(v_dq)
    vq0      = imag(v_dq)

    i_dq     = i0 * exp(-im * (δ0 - pi/2))
    id0      = real(i_dq)
    iq0      = imag(i_dq)
    
    # vd0      = vh * sin(δ0 - θh)
    # vq0      = vh * cos(δ0 - θh)
    
    # id0      = ih0 * sin(δ0 - ϕ0)
    # iq0      = ih0 * cos(δ0 - ϕ0)


    ed0_dash =  (X_q - X_q_dash) * iq0
    
    # ed0_dash0 = vd0 + ra * id0 - X_q_dash  * iq0
    
    eq0_dash = vq0 + ra * iq0 + X_d_dash  * id0   

    # avr
    
    vf0 = eq0_dash + (X_d - X_d_dash) * id0
    
    # Gov
    
    # τe0 = (vq0 + ra * iq0) * iq0 + (vd0 + ra * id0) * id0

    τe0 = ed0_dash * id0 + eq0_dash * iq0 + (X_q_dash - X_d_dash) * id0 * iq0
    
    ω0  = ωs  

    pe0 = τe0 * ω0
    
    pm0 = pe0

    (; δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i )

end



function two_axix_synch_condenser_init(p, power_flow_data)

    # ---------------------------------------
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

    # ---------------------------------------

    # power_flow_data

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

    
    # ---------------------------------------
  
    v        = vh * exp(im * θh)

    u0_r      = real(v)
    
    u0_i      = imag(v)

    # s_α      = αp * ph + αq * im * qh
    
    s_α      = αp * pg + αq * im * qg
    
    i0       = conj(s_α) / conj(v)

    ih0      = abs(i0)
    
    ϕ0       = angle(i0)    
    
    E0       = v + (ra + im * X_q) * i0
    
    δ0       = angle(E0)
    
    v_dq     = v * exp(-im * (δ0 - pi/2))
    vd0      = real(v_dq)
    vq0      = imag(v_dq)

    i_dq     = i0 * exp(-im * (δ0 - pi/2))
    id0      = real(i_dq)
    iq0      = imag(i_dq)
    
    # vd0      = vh * sin(δ0 - θh)
    # vq0      = vh * cos(δ0 - θh)
    
    # id0      = ih0 * sin(δ0 - ϕ0)
    # iq0      = ih0 * cos(δ0 - ϕ0)


    ed0_dash =  (X_q - X_q_dash) * iq0
    
    # ed0_dash0 = vd0 + ra * id0 - X_q_dash  * iq0
    
    eq0_dash = vq0 + ra * iq0 + X_d_dash  * id0   

    # avr
    
    vf0 = eq0_dash + (X_d - X_d_dash) * id0
    
    # Gov
    
    # τe0 = (vq0 + ra * iq0) * iq0 + (vd0 + ra * id0) * id0

    τe0 = ed0_dash * id0 + eq0_dash * iq0 + (X_q_dash - X_d_dash) * id0 * iq0
    
    ω0  = ωs  

    pe0 = τe0 * ω0
    
    pm0 = pe0

    (; δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i )

end


# ------------------------------------------------------
#  self init dynamic loads
# ------------------------------------------------------


function f_t_self_PQ_Const_P( power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v         = vh * exp(im * θh)
    
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )

    return [ vh, θh, ph, qh, i_r, i_i ]
    
end


function f_t_self(component::PQ_Const_P,  power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v        = vh * exp(im * θh)
    
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )

    return [ vh, θh, ph, qh, i_r, i_i ], f_t_self_PQ_Const_P
    
end



function self_init_state(component::PQ_Const_P,  power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]
    
    # p = get_component_params_value(component)

    v         = vh * exp(im * θh)
    
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )
        
    state_var     = [  ]
    algebraic_var = [ real(v), imag(v) ]

    return  (x0 = [state_var...;algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],  aux = [ph, qh, kPL, kQL] )
    
    # return  (state_var, algebraic_var) , (ph, qh, kPL, kQL)
end


function self_init_state(component_type::Type{PQ_Const_P},  power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]
    
    # p = get_component_params_value(component)

    v         = vh * exp(im * θh)
    
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )
        
    state_var     = [  ]
    algebraic_var = [ real(v), imag(v) ]

    return  (x0 = [state_var...;algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],  aux = [ph, qh, kPL, kQL] )
    
    # return  (state_var, algebraic_var) , (ph, qh, kPL, kQL)
end


function f_t_self_PQ_Const_I( power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v         = vh * exp(im * θh)
    
    kPL       = (0.0, ph/vh, 0.0 )
    kQL       = (0.0, qh/vh, 0.0 )

    state_var     = [ ]
    algebraic_var = [ real(v), imag(v)]

    return [ vh, θh, ph, qh, i_r, i_i ] 
    
end


function f_t_self(component::PQ_Const_I,  power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v        = vh * exp(im * θh)
    
    kPL       = (0.0, ph/vh, 0.0 )
    kQL       = (0.0, qh/vh, 0.0 )

    state_var     = [  ]
    algebraic_var = [ real(v), imag(v) ]

    return [ vh, θh, ph, qh, i_r, i_i ], f_t_self_PQ_Const_I 
    
end


function self_init_state(component::PQ_Const_I, power_flow_data)

    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]    

    # p = get_component_params_value(component)

    v         = vh * exp(im * θh)
    
    kPL       = (0.0, ph/vh, 0.0 )
    kQL       = (0.0, qh/vh, 0.0 )

    state_var     = [ ]
    algebraic_var = [ real(v), imag(v) ]

    return (x0 = [state_var...; algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],   aux = [ph, qh, kPL, kQL] )
end


function self_init_state(component_type::Type{PQ_Const_I}, power_flow_data)

    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]    

    # p = get_component_params_value(component)

    v         = vh * exp(im * θh)
    
    kPL       = (0.0, ph/vh, 0.0 )
    kQL       = (0.0, qh/vh, 0.0 )

    state_var     = [ ]
    algebraic_var = [ real(v), imag(v) ]

    return (x0 = [state_var...; algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],   aux = [ph, qh, kPL, kQL] )
end



function f_t_self_PQ_Const_Z(  power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]
    
    v         = vh * exp(im * θh)
    
    kPL       = (0.0, 0.0, ph/(vh^2) )
    kQL       = (0.0, 0.0, qh/(vh^2) )

    return [ vh, θh, ph, qh, i_r, i_i ] 
    
end


function f_t_self(component::PQ_Const_Z,  power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]
    
    v         = vh * exp(im * θh)
    
    kPL       = (0.0, 0.0, ph/(vh^2) )
    kQL       = (0.0, 0.0, qh/(vh^2) )

    return [ vh, θh, ph, qh, i_r, i_i ], f_t_self_PQ_Const_Z 
    
end


function self_init_state(component::PQ_Const_Z, power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]    
    
    # p = get_component_params_value(component)

    v        = vh * exp(im * θh)
    
    kPL       = (0.0, 0.0, ph/(vh^2) )
    kQL       = (0.0, 0.0, qh/(vh^2) )

    state_var     = [ ]
    algebraic_var = [ real(v), imag(v) ]

    return ( x0 = [state_var...; algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],  aux = [ph, qh, kPL, kQL] )
    
end


function self_init_state(component_type::Type{PQ_Const_Z}, power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]    
    
    # p = get_component_params_value(component)

    v        = vh * exp(im * θh)
    
    kPL       = (0.0, 0.0, ph/(vh^2) )
    kQL       = (0.0, 0.0, qh/(vh^2) )

    state_var     = [ ]
    algebraic_var = [ real(v), imag(v) ]

    return ( x0 = [state_var...; algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],  aux = [ph, qh, kPL, kQL] )
    
end



function f_t_self_PQ_dyn_load( power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v    = vh * exp(im * θh)
    
    kPL  = (0.0, 0.0, ph/(vh^2) )
    
    kQL  = (0.0, 0.0, qh/(vh^2) )

    return [ vh, θh, ph, qh, i_r, i_i ]
    
    
end

function f_t_self(component::PQ_dyn_load,  power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v    = vh * exp(im * θh)
    
    kPL  = (0.0, 0.0, ph/(vh^2) )
    
    kQL  = (0.0, 0.0, qh/(vh^2) )

    return [ vh, θh, ph, qh, i_r, i_i ], f_t_self_PQ_dyn_load
    
    
end


function self_init_state(component::PQ_dyn_load, power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]
    
    # p = get_component_params_value(component)

    v       = vh * exp(im * θh)
    
    kPL     = (0.0, 0.0, ph/(vh^2) )
    kQL     = (0.0, 0.0, qh/(vh^2) )

    state_var     = [  ]
    algebraic_var = [ real(v), imag(v) ]

    return (x0 = [state_var...;algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],  aux = [ph, qh, kPL, kQL]  )
    
end


function self_init_state(component_type::Type{PQ_dyn_load}, power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]
    
    # p = get_component_params_value(component)

    v       = vh * exp(im * θh)
    
    kPL     = (0.0, 0.0, ph/(vh^2) )
    kQL     = (0.0, 0.0, qh/(vh^2) )

    state_var     = [  ]
    algebraic_var = [ real(v), imag(v) ]

    return (x0 = [state_var...;algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],  aux = [ph, qh, kPL, kQL]  )
    
end


# ------------------------------------------------------
# local loads 
# ------------------------------------------------------


function self_init_state(component::loc_Load_t1,  gen_init_data)


    _, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    u_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    u_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    u    = u_r + im * u_i
    
    p = get_component_params_value(component)


    loc_P     = p[1]
    loc_Q     = p[2]
    loc_S     = loc_P + loc_Q * im
    

    i  = conj(loc_S) / conj(u)


    return [ real(i), imag(i) ], [ u_r, u_i ]
end



function self_init_state(component_type::Type{loc_Load_t1},  gen_init_data, p )


    _, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    u_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    u_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    u    = u_r + im * u_i
    
    # p = get_component_params_value(component)


    loc_P     = p[1]
    loc_Q     = p[2]
    loc_S     = loc_P + loc_Q * im
    

    i  = conj(loc_S) / conj(u)


    return [ real(i), imag(i) ], [ u_r, u_i ]
end


# ------------------------------------------------------
# Transmission nodes
# ------------------------------------------------------



function f_t_self_Trans_Node( power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    return [ vh, θh, ph, qh, i_r, i_i ]
    
end


function f_t_self(component::Union{Trans_t1_Node, Trans_t2_Node },  power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    return [ vh, θh, ph, qh, i_r, i_i ], f_t_self_Trans_Node
    
end


function self_init_state(component::Union{ Trans_t1_Node, Trans_t2_Node },  power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]
    
    # p = get_component_params_value(component)


    # P     = p[1]
    # Q     = p[2]
    # Y_n   = p[3]
    

    # vh = abs(U)
    # θh = angle(U)

    v  = vh * exp(im * θh)

    # dim      = component.dim
    
    # dict_state_syms = component.dict_state_syms

    # state_algebraic_var = zeros(dim)

    # list_state_and_algb_val  = [ real(v), imag(v)]

    # list_state_and_algb_sym  = [ :u_r, :u_i ]
    
    # for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
    #     state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    # end    

    # return ( x0 = state_algebraic_var, f_t = [i_r, i_i, vh, θh], aux = [ph, qh])

    # # since we only have one type of var


    state_var     = [ ]
    algebraic_var = [ real(v), imag(v) ]

    return (x0 = [state_var...; algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],   aux = [ph, qh] )
    
    # return ( x0 = [ real(v), imag(v)], f_t = [ vh, θh, ph, qh, i_r, i_i ], aux = [ph, qh])

end



function self_init_state(
    component_type::Union{ Type{Trans_t1_Node},
                           Type{Trans_t2_Node} },
    power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    
    v  = vh * exp(im * θh)

    # # since we only have one type of var

    state_var     = [ ]
    algebraic_var = [ real(v), imag(v) ]

    return (x0 = [state_var...; algebraic_var...], f_t = [ vh, θh, ph, qh, i_r, i_i ],   aux = [ph, qh] )

end



function f_t_self(component::Union{plant_Transmission_t1, plant_Transmission_t2},  power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    return [ vh, θh, ph, qh, i_r, i_i ], f_t_self_Trans_Node
    
end


function self_init_state(node::Union{ plant_Transmission_t1, plant_Transmission_t2 }, power_flow_data)

    trans_dict_state_syms = node.Trans.dict_state_syms
    
    # Load
        
    trans_state_algb_var, trans_f_t, trans_aux = self_init_state( node.Trans, power_flow_data)

    return ( x0 = [trans_state_algb_var...], f_t = (; trans_f_t,), aux = trans_aux )
    
end



function self_init_state(
    node_type::Union{ Type{plant_Transmission_t1},
                      Type{plant_Transmission_t2} },
    Transmission_type,
    power_flow_data )

    # trans_dict_state_syms = node.Trans.dict_state_syms
    
    # Load
        
    trans_state_algb_var, trans_f_t, trans_aux =
        self_init_state( Transmission_type,
                         power_flow_data)

    return ( x0 = [trans_state_algb_var...], f_t = (; trans_f_t,), aux = trans_aux )
    
end


#----------------------------------------------------------
#################################
#-----------------------------------------------------


function f_t_self_pf_SM( power_flow_data)
    
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

    # return [i_r, i_i,  vh, θh, ph, qh ]
    return [ig_r, ig_i]
    
end


function f_t_self(component::Union{pf_SM_2axis_v6},   power_flow_data )

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

    return  [ig_r, ig_i ], f_t_self_pf_SM 
end



function self_init_state(component::Union{pf_SM_2axis_v6},  power_flow_data )


    # ----------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    # ----------------------------------------------------
    
    dim  = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)    

    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, pg, qg ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i, :ph, :qh ]
    

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end


function self_init_state(
    component_type::Type{pf_SM_2axis_v6},    
    two_axix_synch_machine_parameters,
    dim, dict_state_syms, power_flow_data )


    # ----------------------------------------------------
    
    #two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    # ----------------------------------------------------
    
    # dim  = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)    

    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, pg, qg ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i, :ph, :qh ]
    

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end




function f_t_self_SM( power_flow_data)
    
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

    # return [i_r, i_i,  vh, θh, ph, qh ]
    return [ig_r, ig_i]
    
end


function f_t_self(component::Union{SM_2axis_idq, SM_2axis_v6},  power_flow_data)
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

    return  [ig_r, ig_i ], f_t_self_SM 
end


function self_init_state(component::Union{SM_2axis_idq, SM_2axis_v6}, power_flow_data)

 
    # ----------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    # ----------------------------------------------------
    
    dim  = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [  δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [  :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end


function self_init_state(
    component_type::Union{Type{SM_2axis_idq},
                          Type{SM_2axis_v6}},    
    two_axix_synch_machine_parameters,
    dim,
    dict_state_syms,
    power_flow_data )

 
    # ----------------------------------------------------
    
    # two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    # ----------------------------------------------------
    
    # dim  = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [  δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [  :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end

#-------------------------------------------------------
#-------------------------------------------------------


function f_t_self_SC_2axis_cd( power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    u = vh * exp( im * θh)

    return   [ig_r, ig_i] # [ig_r, ig_i]
end



function f_t_self(component::Union{SC_2axis_cb_idq, SC_2axis_cb_v6},  power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    u = vh * exp( im * θh)

    return  [ig_r, ig_i ], f_t_self_SC_2axis_cd
    
    # return  [ig_r, ig_i ], f_t_self_SM_2axis_wt_loc_load
end


function self_init_state(component::Union{ SC_2axis_cb_idq, SC_2axis_cb_v6 }, power_flow_data)

 
    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------
    
    dim   = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [  δ0, ω0, ed0_dash, eq0_dash,u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i ], aux = [vf0, τe0])
end



function self_init_state(
    component_type::Union{ Type{SC_2axis_cb_idq},
                           Type{SC_2axis_cb_v6} },    

    two_axix_synch_machine_parameters,
    dim,
    dict_state_syms,
    power_flow_data)

 
    #-------------------------------------------------
    
    # two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------
    
    # dim   = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [  δ0, ω0, ed0_dash, eq0_dash,u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i ], aux = [vf0, τe0])
end




#####

function f_t_self_SC_2axis_wt_loc_load( power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    u = vh * exp( im * θh)

    return   [i_r, i_i] # [ig_r, ig_i]
end



function f_t_self(component::Union{SC_2axis_wt_loc_load_cb_idq, SC_2axis_wt_loc_load_cb_v6},  power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    u = vh * exp( im * θh)

    return  [i_r, i_i ], f_t_self_SC_2axis_wt_loc_load
    
    # return  [ig_r, ig_i ], f_t_self_SM_2axis_wt_loc_load
end


function self_init_state(component::Union{ SC_2axis_wt_loc_load_cb_idq, SC_2axis_wt_loc_load_cb_v6 }, power_flow_data)

    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------
    
    dim   = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [i_r, i_i ], aux = [vf0, τe0])
end



function self_init_state(
    component_type::Union{ Type{SC_2axis_wt_loc_load_cb_idq},
                           Type{SC_2axis_wt_loc_load_cb_v6} },
    two_axix_synch_machine_parameters,
    dim,
    dict_state_syms,
    power_flow_data )

    #-------------------------------------------------
    
    # two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------
    
    # dim   = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [i_r, i_i ], aux = [vf0, τe0])
end



#----------------------------------------------------------
#----------------------------------------------------------


function f_t_self_SM_2axis_wt_loc_load( power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    u = vh * exp( im * θh)

    return   [i_r, i_i] # [ig_r, ig_i]
end



function f_t_self(component::Union{SM_2axis_wt_loc_load_cb_idq, SM_2axis_wt_loc_load_cb_v6},  power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    u = vh * exp( im * θh)

    return  [i_r, i_i ], f_t_self_SM_2axis_wt_loc_load
    
    # return  [ig_r, ig_i ], f_t_self_SM_2axis_wt_loc_load
end


function self_init_state(component::Union{ SM_2axis_wt_loc_load_cb_idq, SM_2axis_wt_loc_load_cb_v6 }, power_flow_data)

    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------
    

    dim   = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [i_r, i_i ], aux = [vf0, τe0])

end




function self_init_state(
    component_type::Union{ Type{SM_2axis_wt_loc_load_cb_idq},
                           Type{SM_2axis_wt_loc_load_cb_v6} },
    two_axix_synch_machine_parameters,
    dim,
    dict_state_syms,
    power_flow_data)

    #-------------------------------------------------
    
    # two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------
    

    # dim   = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [i_r, i_i ], aux = [vf0, τe0])

end



#----------------------------------------------------------


function f_t_self_SM_2axis( power_flow_data)
    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    return  [ig_r, ig_i]
end


function f_t_self(component::Union{SM_2axis_cb_idq},  power_flow_data)
    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    return  [ig_r, ig_i], f_t_self_SM_2axis
end


function self_init_state(component::Union{SM_2axis_cb_idq}, power_flow_data)

    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    dim            = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end



function self_init_state(
    component_type::Type{SM_2axis_cb_idq},
    
    two_axix_synch_machine_parameters,
    dim,
    dict_state_syms,
    power_flow_data)

    #-------------------------------------------------
    
    # two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    # dim            = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end



function f_t_self(component::Union{SM_2axis_cb_v6}, power_flow_data)
    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    return [ig_r, ig_i], f_t_self_SM_2axis
end


function self_init_state(component::Union{SM_2axis_cb_v6},
                         power_flow_data)

 
    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  
    
    dim            = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end



function self_init_state(
    component_type::Type{SM_2axis_cb_v6},
    
    two_axix_synch_machine_parameters,
    dim,
    dict_state_syms,
    power_flow_data)

 
    #-------------------------------------------------
    
    # two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  
    
    # dim            = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end


function f_t_self(component::SM_2axis_cb_direct,  power_flow_data)
    
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
    
    return [ig_r, ig_i], f_t_self_SM_2axis
end


function self_init_state(component::SM_2axis_cb_direct,
                         power_flow_data)


    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    dim            = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])    
    
end



function self_init_state(
    component_type::Type{SM_2axis_cb_direct},
    
    two_axix_synch_machine_parameters,
    dim,
    dict_state_syms,
    power_flow_data)


    #-------------------------------------------------
    
    # two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    # dim            = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])    
    
end



function f_t_self(component::SM_2axis_cb_millano,  power_flow_data)
    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]    

    return [ig_r, ig_i], f_t_self_SM_2axis
    
end



function self_init_state(component::SM_2axis_cb_millano,  power_flow_data)

 
    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    dim            = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, vf0, τe0, τe0 ]
    
    list_state_and_algb_sym  = [  :δ, :ω, :ed_dash, :eq_dash , :u_r, :u_i, :u_r, :u_i, :vf, :τm,  :τe ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])    
    
end



function self_init_state(
    component_type::Type{SM_2axis_cb_millano},
    
    two_axix_synch_machine_parameters,
    dim,
    dict_state_syms,
    power_flow_data)

 
    #-------------------------------------------------
    
    # two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    # dim            = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, vf0, τe0, τe0 ]
    
    list_state_and_algb_sym  = [  :δ, :ω, :ed_dash, :eq_dash , :u_r, :u_i, :u_r, :u_i, :vf, :τm,  :τe ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])    
    
end



function f_t_self_SM_cb_inf( power_flow_data)
    
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
    
    return [ig_r, ig_i, vh, θh]
    
end




function f_t_self(component::SM_2axis_cb_inf,  power_flow_data)
    
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
    
    return [ig_r, ig_i, vh, θh], f_t_self_SM_cb_inf
    
end



function self_init_state(component::SM_2axis_cb_inf,  power_flow_data)


    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    dim            = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash0, eq0_dash0, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])    

    
end


function f_t_self(component::SM_2axis_cb_inf_bus,  power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    return [ig_r, ig_i, vh, θh], f_t_self_SM_cb_inf

    
end


function self_init_state(component::SM_2axis_cb_inf_bus,  power_flow_data)

 
    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    dim            = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [ δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end


function f_t_self(component::SM_2axis_cb,  power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    return [ig_r, ig_i], f_t_self_SM_2axis

end


function self_init_state(component::SM_2axis_cb,  power_flow_data)

 
    #-------------------------------------------------
    
    two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    dim            = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [  δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, vf0, τm0 ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash,:u_r, :u_i, :vf, :τm  ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end



function self_init_state(
    component_type::Type{SM_2axis_cb},
    
    two_axix_synch_machine_parameters,
    dim,
    dict_state_syms,
    power_flow_data)

 
    #-------------------------------------------------
    
    # two_axix_synch_machine_parameters = get_component_params_value(component)
    
    δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, τe0, vf0, vh, θh, ph, qh, i_r, i_i, pg, qg, ig_r, ig_i = two_axix_synch_machine_init(two_axix_synch_machine_parameters, power_flow_data)

    #-------------------------------------------------  

    # dim            = component.dim
    
    # dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)
    
    list_state_and_algb_val  = [  δ0, ω0, ed0_dash, eq0_dash, u0_r, u0_i, vf0, τm0 ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash,:u_r, :u_i, :vf, :τm  ]

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [ig_r, ig_i], aux = [vf0, τe0])
    
end



function f_t_self(component::Infinite_cb_bus,  power_flow_data)
    
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

    return [ig_r, ig_i, vh, θh], f_t_self_SM_cb_inf

    
end


function self_init_state(component::Infinite_cb_bus,  power_flow_data)

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
    
    p = get_component_params_value(component)

    U     = p[1]
    P     = p[2]
    Q     = p[3]
    Y_n   = p[4]    

    # vh = abs(U)
    # θh = angle(U)

    v        = vh * exp(im * θh)    

    state_var     = [ ]
    algebraic_var = [ real(v), imag(v) ]


    return ( x0 = [state_var...;algebraic_var...], f_t = [ig_r, ig_i, vh, θh], aux = [vh, θh, ph, qh]) 
end



function f_t_self(component::Infinite_bus,  power_flow_data)
    
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
    
    return [ig_r, ig_i, vh, θh], f_t_self_SM_cb_inf
    
end


function self_init_state(component::Infinite_bus,  power_flow_data)

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
    
    p = get_component_params_value(component)

    U     = p[1]
    P     = p[2]
    Q     = p[3]
    Y_n   = p[4]    

    # vh = abs(U)
    # θh = angle(U)

    v        = vh * exp(im * θh)    

    state_var     = [ ]
    algebraic_var = [ real(v), imag(v) ]


    return ( x0 = [state_var...;algebraic_var...], f_t = [ig_r, ig_i, vh, θh], aux = [vh, θh, ph, qh] ) 

end


# ------------------------------------------------------
#  self init controllers
# ------------------------------------------------------


function self_init_state(component::gov_ieee_tgov1_cb,  gen_init_data)
    
    _, pm0, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    # gen_state_var = [gen_state_var...;gen_algebraic_var...]
    
    ω0  = gen_state_algb_var[ gen_dict_state_syms[ :ω ] ]
    
    p   = get_component_params_value(component)

    # Gov 8
    
    T1       = p[1]
    T2       = p[2]
    T3       = p[3]
    Dt       = p[4]
    p_max    = p[5]
    p_min    = p[6]
    R        = p[7]
        
    # Gov

    ω_ref0      =  ω0_ref = ω0
    
    τm0         = pm0  # / ω0

    τm0_tilade  = τm0


    b = [0.0, τm0_tilade]

    A = [(1 - T2/T3)  -1.0;
         (T2/T3)       1.0
         ]
    
    x = A \ b

    x10     = x[1]
    x20     = x[2]

    p_order0 = x10

    state_var     = [x10, x20]
    
    algebraic_var = [τm0_tilade, ω0_ref]

    return [state_var...;algebraic_var...], [p_order0, ω_ref0]

end


function self_init_state(component_type::Type{gov_ieee_tgov1_cb},  gen_init_data, p)
    
    _, pm0, gen_state_algb_var, gen_dict_state_syms = gen_init_data

    
    ω0  = gen_state_algb_var[ gen_dict_state_syms[ :ω ] ]
    

    # Gov 8
    
    T1       = p[1]
    T2       = p[2]
    T3       = p[3]
    Dt       = p[4]
    p_max    = p[5]
    p_min    = p[6]
    R        = p[7]
        
    # Gov

    ω_ref0      =  ω0_ref = ω0
    
    τm0         = pm0  # / ω0

    τm0_tilade  = τm0


    b = [0.0, τm0_tilade]

    A = [(1 - T2/T3)  -1.0;
         (T2/T3)       1.0
         ]
    
    x = A \ b

    x10     = x[1]
    x20     = x[2]

    p_order0 = x10

    state_var     = [x10, x20]
    
    algebraic_var = [τm0_tilade, ω0_ref]

    return [state_var...;algebraic_var...], [p_order0, ω_ref0]

end



function self_init_state(component::gov_t0_cb,  gen_init_data)
    
    _, pm0, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    # gen_state_var = [gen_state_var...;gen_algebraic_var...]
    
    ω0  = gen_state_algb_var[ gen_dict_state_syms[ :ω ] ]
    
    p   = get_component_params_value(component)

    # Gov 8
    
    T3       = p[1]
    T4       = p[2]
    T5       = p[3]
    Tc       = p[4]
    Ts       = p[5]
    p_max    = p[6]
    p_min    = p[7]
    R        = p[8]
        
    # Gov

    ω_ref0      =  ω0_ref = ω0
    
    τm0         = pm0 # / ω0

    τm0_tilade  = τm0

    τm_hat      = τm0_tilade

    b = [0.0, 0.0, τm0_tilade]

    A = [(1 - T3/Tc)          (-1.0)        0.0;
         (1 - T4/T5)*(T3/Tc)  (1 - T4/T5)   -1.0;
         (T4/T5)*(T3/Tc)      (T4/T5)       1.0
         ]
    
    x = A \ b

    xg10     = x[1]
    xg20     = x[2]
    xg30     = x[3]

    p0_in = xg10
    
    phat0_in = p0_in

    p_order0 = phat0_in - 1/R * (ω0_ref/ω0 - 1.0 )

    state_var     = [xg10, xg20, xg30]
    
    algebraic_var = [τm0_tilade, ω0_ref, phat0_in, p0_in]

    return [state_var...;algebraic_var...], [p_order0, ω_ref0]

end


function self_init_state(component_type::Type{gov_t0_cb},  gen_init_data, p)
    
    _, pm0, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    
    ω0  = gen_state_algb_var[ gen_dict_state_syms[ :ω ] ]
    

    # Gov 8
    
    T3       = p[1]
    T4       = p[2]
    T5       = p[3]
    Tc       = p[4]
    Ts       = p[5]
    p_max    = p[6]
    p_min    = p[7]
    R        = p[8]
        
    # Gov

    ω_ref0      =  ω0_ref = ω0
    
    τm0         = pm0 # / ω0

    τm0_tilade  = τm0

    τm_hat      = τm0_tilade

    b = [0.0, 0.0, τm0_tilade]

    A = [(1 - T3/Tc)          (-1.0)        0.0;
         (1 - T4/T5)*(T3/Tc)  (1 - T4/T5)   -1.0;
         (T4/T5)*(T3/Tc)      (T4/T5)       1.0
         ]
    
    x = A \ b

    xg10     = x[1]
    xg20     = x[2]
    xg30     = x[3]

    p0_in = xg10
    
    phat0_in = p0_in

    p_order0 = phat0_in - 1/R * (ω0_ref/ω0 - 1.0 )

    state_var     = [xg10, xg20, xg30]
    
    algebraic_var = [τm0_tilade, ω0_ref, phat0_in, p0_in]

    return [state_var...;algebraic_var...], [p_order0, ω_ref0]

end



function self_init_state(component::gov_t1_cb,  gen_init_data)
    
    _, pm0, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    # gen_state_var = [gen_state_var...;gen_algebraic_var...]
    
    ω0  = gen_state_algb_var[ gen_dict_state_syms[ :ω ] ]
    
    p   = get_component_params_value(component)

    # Gov 8
    
    T3       = p[1]
    T4       = p[2]
    T5       = p[3]
    Tc       = p[4]
    Ts       = p[5]
    p_max    = p[6]
    p_min    = p[7]
    R        = p[8]
        
    # Gov

    ω_ref0      =  ω0_ref = ω0
    
    τm0         = pm0 # / ω0

    τm0_tilade  = τm0

    τm_hat      = τm0_tilade

    b = [0.0, 0.0, τm0_tilade]

    A = [(1 - T3/Tc)          (-1.0)        0.0;
         (1 - T4/T5)*(T3/Tc)  (1 - T4/T5)   -1.0;
         (T4/T5)*(T3/Tc)      (T4/T5)       1.0
         ]
    
    x = A \ b

    xg10     = x[1]
    xg20     = x[2]
    xg30     = x[3]

    phat0_in = xg10

    p_order0 = phat0_in - 1/R * (ω0_ref/ω0 - 1.0 )

    state_var     = [xg10, xg20, xg30]
    
    algebraic_var = [τm0_tilade, ω0_ref, phat0_in]

    return [state_var...;algebraic_var...], [p_order0, ω_ref0]

end



function self_init_state(component_type::Type{gov_t1_cb},  gen_init_data, p )
    
    _, pm0, gen_state_algb_var, gen_dict_state_syms = gen_init_data

    
    ω0  = gen_state_algb_var[ gen_dict_state_syms[ :ω ] ]
    
    # Gov 8
    
    T3       = p[1]
    T4       = p[2]
    T5       = p[3]
    Tc       = p[4]
    Ts       = p[5]
    p_max    = p[6]
    p_min    = p[7]
    R        = p[8]
        
    # Gov

    ω_ref0      =  ω0_ref = ω0
    
    τm0         = pm0 # / ω0

    τm0_tilade  = τm0

    τm_hat      = τm0_tilade

    b = [0.0, 0.0, τm0_tilade]

    A = [(1 - T3/Tc)          (-1.0)        0.0;
         (1 - T4/T5)*(T3/Tc)  (1 - T4/T5)   -1.0;
         (T4/T5)*(T3/Tc)      (T4/T5)       1.0
         ]
    
    x = A \ b

    xg10     = x[1]
    xg20     = x[2]
    xg30     = x[3]

    phat0_in = xg10

    p_order0 = phat0_in - 1/R * (ω0_ref/ω0 - 1.0 )

    state_var     = [xg10, xg20, xg30]
    
    algebraic_var = [τm0_tilade, ω0_ref, phat0_in]

    return [state_var...;algebraic_var...], [p_order0, ω_ref0]

end



function self_init_state(component::gov_t1_cb_sauer,  gen_init_data)
    
    _, pm0, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    # gen_state_var = [gen_state_var...;gen_algebraic_var...]
    
    ω0  = gen_state_algb_var[ gen_dict_state_syms[ :ω ] ]
    
    p   = get_component_params_value(component)

    # Gov 8
    
    Tc       = p[1]
    Ts       = p[2]
    p_max    = p[3]
    p_min    = p[4]
    R        = p[5]
        
    # Gov

    ω_ref0      =  ω0_ref = ω0
    
    τm0_tilade  = pm0 

    xg20        = τm0_tilade 

    xg10        = xg20 

    p_order0 = xg10 

    state_var     = [xg10, xg20]
    
    algebraic_var = [ω0_ref, τm0_tilade ]

    return [state_var...; algebraic_var...], [p_order0, ω_ref0]

end



function self_init_state(component_type::Type{gov_t1_cb_sauer},  gen_init_data, p)
    
    _, pm0, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    
    ω0  = gen_state_algb_var[ gen_dict_state_syms[ :ω ] ]
    

    # Gov 8
    
    Tc       = p[1]
    Ts       = p[2]
    p_max    = p[3]
    p_min    = p[4]
    R        = p[5]
        
    # Gov

    ω_ref0      =  ω0_ref = ω0
    
    τm0_tilade  = pm0 

    xg20        = τm0_tilade 

    xg10        = xg20 

    p_order0 = xg10 

    state_var     = [xg10, xg20]
    
    algebraic_var = [ω0_ref, τm0_tilade ]

    return [state_var...; algebraic_var...], [p_order0, ω_ref0]

end



function self_init_state(component::avr_t0_cb,
                         gen_init_data )

    vf0, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data

    # gen_state_var = [gen_state_var...;gen_algebraic_var...]
    # vf0 = gen_algebraic_var[6]

    v_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    v_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    v    = v_r + im * v_i

    vh   = abs(v)

    p = get_component_params_value(component)

    # Exc 
    Ta         = p[1]
    Te         = p[2]
    Tf         = p[3]
    Tr         = p[4]
    Ka         = p[5]
    Ke         = p[6]
    Kf         = p[7]
    V_R_max    = p[8]
    V_R_min    = p[9]
    Ae         = p[10]
    Be         = p[11]

    vf_tilade0 = vf0
    
    vm0        = vh

    vr10_hat   = vf_tilade0 * (Ke + Sevf(Ae, Be, vf_tilade0))
    
    vr10       = vr10_hat
    
    vr20       = (Kf / Tf) * vf_tilade0
    
    v_ref0     = (vr10 / Ka) + vm0 - vr20 + (Kf / Tf) * vf_tilade0

    state_var  = [vm0, vr10, vr20, vf_tilade0]
    
    algebraic_var = [v_ref0, vr10_hat, vf0]

    return [state_var...;algebraic_var...], [v_ref0]
    
end



function self_init_state(component_type::Type{avr_t0_cb},
                         gen_init_data, p )

    vf0, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    v_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    v_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    v    = v_r + im * v_i

    vh   = abs(v)


    # Exc 
    Ta         = p[1]
    Te         = p[2]
    Tf         = p[3]
    Tr         = p[4]
    Ka         = p[5]
    Ke         = p[6]
    Kf         = p[7]
    V_R_max    = p[8]
    V_R_min    = p[9]
    Ae         = p[10]
    Be         = p[11]

    vf_tilade0 = vf0
    
    vm0        = vh

    vr10_hat   = vf_tilade0 * (Ke + Sevf(Ae, Be, vf_tilade0))
    
    vr10       = vr10_hat
    
    vr20       = (Kf / Tf) * vf_tilade0
    
    v_ref0     = (vr10 / Ka) + vm0 - vr20 + (Kf / Tf) * vf_tilade0

    state_var  = [vm0, vr10, vr20, vf_tilade0]
    
    algebraic_var = [v_ref0, vr10_hat, vf0]

    return [state_var...;algebraic_var...], [v_ref0]
    
end



function self_init_state(component::avr_t1_cb,
                         gen_init_data )

    vf0, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data

    # gen_state_var = [gen_state_var...;gen_algebraic_var...]
    # vf0 = gen_algebraic_var[6]

    v_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    v_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    v    = v_r + im * v_i

    vh   = abs(v)

    p = get_component_params_value(component)

    # Exc 
    Ta         = p[1]
    Te         = p[2]
    Tf         = p[3]
    Tr         = p[4]
    Ka         = p[5]
    Ke         = p[6]
    Kf         = p[7]
    V_R_max    = p[8]
    V_R_min    = p[9]
    Ae         = p[10]
    Be         = p[11]

    vf_tilade0 = vf0
    
    vm0        = vh
    
    vr10       = vf_tilade0 * (Ke + Sevf(Ae, Be, vf_tilade0))

    
    vr20       = (Kf / Tf) * vf_tilade0
    
    v_ref0     = (vr10 / Ka) + vm0 - vr20 + (Kf / Tf) * vf_tilade0

    state_var  = [vm0, vr10, vr20, vf_tilade0]
    
    algebraic_var = [v_ref0, vf0]

    return [state_var...;algebraic_var...], [v_ref0]
    
end



function self_init_state(component_type::Type{avr_t1_cb},
                         gen_init_data, p )

    vf0, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    v_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    v_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    v    = v_r + im * v_i

    vh   = abs(v)


    # Exc 
    Ta         = p[1]
    Te         = p[2]
    Tf         = p[3]
    Tr         = p[4]
    Ka         = p[5]
    Ke         = p[6]
    Kf         = p[7]
    V_R_max    = p[8]
    V_R_min    = p[9]
    Ae         = p[10]
    Be         = p[11]

    vf_tilade0 = vf0
    
    vm0        = vh
    
    vr10       = vf_tilade0 * (Ke + Sevf(Ae, Be, vf_tilade0))

    
    vr20       = (Kf / Tf) * vf_tilade0
    
    v_ref0     = (vr10 / Ka) + vm0 - vr20 + (Kf / Tf) * vf_tilade0

    state_var  = [vm0, vr10, vr20, vf_tilade0]
    
    algebraic_var = [v_ref0, vf0]

    return [state_var...;algebraic_var...], [v_ref0]
    
end



function self_init_state(component::avr_t1_cb_sauer,  gen_init_data)

    vf0, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data

    # gen_state_var = [gen_state_var...;gen_algebraic_var...]
    # vf0 = gen_algebraic_var[6]

    v_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    v_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    v    = v_r + im * v_i

    vh   = abs(v)

    p = get_component_params_value(component)

    # Exc 
    Ta         = p[1]
    Te         = p[2]
    Tf         = p[3]
    Tr         = p[4]
    Ka         = p[5]
    Ke         = p[6]
    Kf         = p[7]
    V_R_max    = p[8]
    V_R_min    = p[9]
    Ae         = p[10]
    Be         = p[11]

    vf_tilade0 = vf0
        
    vr10       = vf_tilade0 * (Ke + Sevf(Ae, Be, vf_tilade0))
    
    vr20       = (Kf / Tf) * vf_tilade0
    
    v_ref0     = vh + (vr10 / Ka)

    state_var  = [ vr10, vr20, vf_tilade0]
    
    algebraic_var = [v_ref0, vf0]

    return [state_var...;algebraic_var...], [v_ref0]
    
end


function self_init_state(component_type::Type{avr_t1_cb_sauer},  gen_init_data, p)

    vf0, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    v_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    v_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    v    = v_r + im * v_i

    vh   = abs(v)

    # Exc 
    Ta         = p[1]
    Te         = p[2]
    Tf         = p[3]
    Tr         = p[4]
    Ka         = p[5]
    Ke         = p[6]
    Kf         = p[7]
    V_R_max    = p[8]
    V_R_min    = p[9]
    Ae         = p[10]
    Be         = p[11]

    vf_tilade0 = vf0
        
    vr10     = vf_tilade0 * (Ke + Sevf(Ae, Be, vf_tilade0))
    
    vr20       = (Kf / Tf) * vf_tilade0
    
    v_ref0     = vh + (vr10 / Ka)

    state_var  = [ vr10, vr20, vf_tilade0]
    
    algebraic_var = [v_ref0, vf0]

    return [state_var...;algebraic_var...], [v_ref0]
    
end




function self_init_state(component::avr_t1_with_pss_cb, gen_init_data)

    vf0, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data

    # vf0 = gen_algebraic_var[6]

    # gen_state_var = [gen_state_var...;gen_algebraic_var...]
    
    v_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    v_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    v    = v_r + im * v_i

    vh   = abs(v)

    p = get_component_params_value(component)

    # Exc 
    Ta         = p[1]
    Te         = p[2]
    Tf         = p[3]
    Tr         = p[4]
    Ka         = p[5]
    Ke         = p[6]
    Kf         = p[7]
    V_R_max    = p[8]
    V_R_min    = p[9]
    Ae         = p[10]
    Be         = p[11]

    vf_tilade0 = vf0
    
    vm0        = vh
    
    vr1_thr    = vf_tilade0 * (Ke + Sevf(Ae, Be, vf_tilade0))

    vr10       = vr1_thr
    
    vr20       = (Kf / Tf) * vf_tilade0
    
    v_ref0     = vr10 / Ka + vm0 - vr20 + Kf / Tf * vf_tilade0

    state_var  = [vm0, vr10, vr20, vf_tilade0]

    # algebraic_var = [ ]
    algebraic_var = [v_ref0, vf0]

    #  return ([v_ref0,  state_var, algebraic_var)
    return [state_var...;algebraic_var...], [v_ref0]
    
end




function self_init_state(component_type::Type{avr_t1_with_pss_cb}, gen_init_data, p )

    vf0, _, gen_state_algb_var, gen_dict_state_syms = gen_init_data


    
    v_r  = gen_state_algb_var[ gen_dict_state_syms[ :u_r ] ]
    
    v_i  = gen_state_algb_var[ gen_dict_state_syms[ :u_i ] ]

    v    = v_r + im * v_i

    vh   = abs(v)


    # Exc 
    Ta         = p[1]
    Te         = p[2]
    Tf         = p[3]
    Tr         = p[4]
    Ka         = p[5]
    Ke         = p[6]
    Kf         = p[7]
    V_R_max    = p[8]
    V_R_min    = p[9]
    Ae         = p[10]
    Be         = p[11]

    vf_tilade0 = vf0
    
    vm0        = vh
    
    vr1_thr    = vf_tilade0 * (Ke + Sevf(Ae, Be, vf_tilade0))

    vr10       = vr1_thr
    
    vr20       = (Kf / Tf) * vf_tilade0
    
    v_ref0     = vr10 / Ka + vm0 - vr20 + Kf / Tf * vf_tilade0

    state_var  = [vm0, vr10, vr20, vf_tilade0]

    # algebraic_var = [ ]
    algebraic_var = [v_ref0, vf0]

    #  return ([v_ref0,  state_var, algebraic_var)
    return [state_var...;algebraic_var...], [v_ref0]
    
end




function self_init_state(component::pss_t2_cb,  gen_init_data)

    vf0, pe0, gen_state_algb_var, gen_dict_state_syms = gen_init_data

    # gen_state_var = [gen_state_var...;gen_algebraic_var...]
    
    ω0   = gen_state_algb_var[ gen_dict_state_syms[ :ω ]]

    p = get_component_params_value(component)

    # PSS
    
    T1     = p[1]
    T2     = p[2]
    T3     = p[3]
    T4     = p[4]
    Tw     = p[5]
    Kw     = p[6]
    vs_max = p[7]
    vs_min = p[8]

    ωs     = 1.0
    
    ω_ref0 =  ω0
    
    v_SI = (ω_ref0 -  ω0)/ ωs # since ω_ref0 =  ω

    v10  = Kw * v_SI 
    
    v20  =  (1 - T1/T2)*(Kw * v_SI + v10)
    
    v30  = (1 - T3/T4)*(v20 + ((T1/T2) * (Kw * v_SI + v10))) 
    
    vs0  = v30 + (T3/T4) * (v20 + (T1/T2) * (Kw * v_SI + v10))

    # vs0  = 0.0 # Since v_ref0 = v_ref
     
    state_var     = [v10, v20, v30]
    algebraic_var = [vs0]

    return [state_var...;algebraic_var...], [ω_ref0]
    
end



function self_init_state(component_type::Type{pss_t2_cb},  gen_init_data, p )

    vf0, pe0, gen_state_algb_var, gen_dict_state_syms = gen_init_data

    
    ω0   = gen_state_algb_var[ gen_dict_state_syms[ :ω ]]


    # PSS
    
    T1     = p[1]
    T2     = p[2]
    T3     = p[3]
    T4     = p[4]
    Tw     = p[5]
    Kw     = p[6]
    vs_max = p[7]
    vs_min = p[8]

    ωs     = 1.0
    
    ω_ref0 =  ω0
    
    v_SI = (ω_ref0 -  ω0)/ ωs # since ω_ref0 =  ω

    v10  = Kw * v_SI 
    
    v20  =  (1 - T1/T2)*(Kw * v_SI + v10)
    
    v30  = (1 - T3/T4)*(v20 + ((T1/T2) * (Kw * v_SI + v10))) 
    
    vs0  = v30 + (T3/T4) * (v20 + (T1/T2) * (Kw * v_SI + v10))

    # vs0  = 0.0 # Since v_ref0 = v_ref
     
    state_var     = [v10, v20, v30]
    algebraic_var = [vs0]

    return [state_var...;algebraic_var...], [ω_ref0]
    
end



# ------------------------------------------------------
# self init plant load 
# ------------------------------------------------------


function f_t_self(component::plant_PQ_Const_P,  power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v        = vh * exp(im * θh)
    
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )

    return [ vh, θh, ph, qh, i_r, i_i ], f_t_self_PQ_Const_P
    
end



function f_t_self(component::plant_PQ_Const_I,  power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v        = vh * exp(im * θh)
    
    kPL       = (0.0, ph/vh, 0.0 )
    kQL       = (0.0, qh/vh, 0.0 )    

    return [ vh, θh, ph, qh, i_r, i_i ], f_t_self_PQ_Const_I
    
end



function f_t_self(component::plant_PQ_Const_Z,  power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v        = vh * exp(im * θh)
    
    kPL       = (0.0, 0.0, ph/(vh^2) )
    kQL       = (0.0, 0.0, qh/(vh^2) )

    return [ vh, θh, ph, qh, i_r, i_i ], f_t_self_PQ_Const_Z
    
end



function f_t_self(component::plant_PQ_dyn_load,  power_flow_data)
    
    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    v        = vh * exp(im * θh)
    
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )

    return [vh, θh, ph, qh, i_r, i_i], f_t_self_PQ_dyn_load
    
end



function self_init_state(node::Union{ plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load }, power_flow_data)

    load_dict_state_syms = node.Load.dict_state_syms
    
    # Load
        
    load_state_algb_var, load_f_t, load_aux = self_init_state( node.Load, power_flow_data)

    return ( x0 = [load_state_algb_var...], f_t = (; load_f_t,), aux = load_aux )
    
end



function self_init_state(
    node_type::Union{ Type{plant_PQ_Const_P},
                      Type{plant_PQ_Const_I},
                      Type{plant_PQ_Const_Z},
                      Type{plant_PQ_dyn_load} },
    node_load_type,
    power_flow_data)
    

    # dispatch based on load type
    load_state_algb_var, load_f_t, load_aux =
        self_init_state( node_load_type, power_flow_data)

    return ( x0 = [load_state_algb_var...], f_t = (; load_f_t,), aux = load_aux )
    
end


# ------------------------------------------------------
# self init plants  SM 
# ------------------------------------------------------

function f_t_self(component::Union{ pf_plant_SM_v6, plant_SM_idq, plant_SM_v6, plant_SM_system_matrix },  power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]

    return [i_r, i_i ], f_t_self_SM
    
end


function self_init_state(node::Union{ pf_plant_SM_v6, plant_SM_idq, plant_SM_v6, plant_SM_system_matrix }, power_flow_data)

    gen_dict_state_syms = node.Gen.dict_state_syms

    plant_dict_state_syms = node.dict_state_syms
    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 = self_init_state( node.Gen, power_flow_data)

    vf0, pm0 = vf0_pm0

    return ( x0 = [gen_state_algb_var...;], f_t = (; gen_f_t, vf0, pm0  ), aux = [vf0, pm0])
    
end


function self_init_state(
    node_type::Union{ Type{pf_plant_SM_v6},
                      Type{plant_SM_idq},
                      Type{plant_SM_v6},
                      Type{plant_SM_system_matrix} },
    
    gen_type,
    gen_para_values,
    gen_dim,
    gen_dict_state_syms,
    power_flow_data )

    # gen_dict_state_syms = node.Gen.dict_state_syms

    # plant_dict_state_syms = node.dict_state_syms
    
    # dispatch based on Gen type
    gen_state_algb_var,  gen_f_t, vf0_pm0 =
        self_init_state( gen_type,
                         
                         gen_para_values,
                         gen_dim,
                         gen_dict_state_syms,
                         power_flow_data)

    vf0, pm0 = vf0_pm0

    return ( x0 = [gen_state_algb_var...;], f_t = (; gen_f_t, vf0, pm0  ), aux = [vf0, pm0])
    
end


# ------------------------------------------------------
# self init plants  
# ------------------------------------------------------



function f_t_self(component::Union{plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_system_matrix },  power_flow_data)

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]    

    return [i_r, i_i], f_t_self_SC_2axis_cd
    
end



function self_init_state(node::Union{plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_system_matrix}, power_flow_data)

    gen_dict_state_syms = node.Gen.dict_state_syms

    plant_dict_state_syms = node.dict_state_syms
    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 = self_init_state( node.Gen, power_flow_data)

    vf0, pm0 = vf0_pm0

    gen_init_data = (vf0, pm0, gen_state_algb_var, gen_dict_state_syms)
    

    # Exc
    exc_state_algb_var, avr_f_t = self_init_state(node.Exc, gen_init_data )
    
    return ( x0 = [gen_state_algb_var...; exc_state_algb_var...], f_t = (; gen_f_t,  avr_f_t), aux = [vf0, pm0])
    
end



function self_init_state(node_type::Union{Type{plant_no_gov_idq}, Type{plant_no_gov_v6}, Type{plant_no_gov_system_matrix} },  gen_type, gen_para_values, gen_dim, gen_dict_state_syms, avr_type, avr_para_values, power_flow_data )

    # plant_dict_state_syms = node.dict_state_syms
    
    # dispatch on Gen type
    gen_state_algb_var,  gen_f_t, vf0_pm0 =
        self_init_state(
            gen_type,            
            gen_para_values,
            gen_dim,
            gen_dict_state_syms,
            power_flow_data)

    
    vf0, pm0 = vf0_pm0

    gen_init_data = (vf0, pm0, gen_state_algb_var, gen_dict_state_syms)
    

    # dispatch on Exc type
    exc_state_algb_var, avr_f_t =
        self_init_state(
            avr_type,
            gen_init_data,
            avr_para_values )
    
    return ( x0 = [gen_state_algb_var...; exc_state_algb_var...], f_t = (; gen_f_t,  avr_f_t), aux = [vf0, pm0])
    
end


#--------------------------------------------------


function f_t_self_plant_no_gov_no_avr_wt_loc_load(power_flow_data)

    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]

    return [i_r, i_i ]
    
end


function f_t_self(component::Union{plant_no_gov_no_avr_wt_loc_load },  power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]    

    u = vh * exp( im * vh )
    
    return [i_r, i_i ], f_t_self_plant_no_gov_no_avr_wt_loc_load
    
    # return [ig_r, ig_i ], f_t_self_SM_2axis_wt_loc_load
    
end


function self_init_state(node::Union{ plant_no_gov_no_avr_wt_loc_load }, power_flow_data)

    gen_dict_state_syms = node.Gen.dict_state_syms

    plant_dict_state_syms = node.dict_state_syms
    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 = self_init_state( node.Gen, power_flow_data)

    vf0, pm0 = vf0_pm0

    gen_init_data = (vf0, pm0, gen_state_algb_var, gen_dict_state_syms)
    
    # loc_load
    loc_load_state_algb_var, loc_load_f_t = self_init_state(node.Loc_load, gen_init_data )    
    
    return ( x0 = [gen_state_algb_var...;  loc_load_state_algb_var... ], f_t = (; gen_f_t, vf0, pm0, loc_load_f_t), aux = [vf0, pm0])
    
end



function self_init_state(
    node_type::Type{ plant_no_gov_no_avr_wt_loc_load },
    gen_type,
    gen_para_values,
    gen_dim,
    gen_dict_state_syms,
    loc_load_type,
    loc_load_para_values, power_flow_data )

    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 =
        self_init_state(
            gen_type,            
            gen_para_values,
            gen_dim,
            gen_dict_state_syms,
            power_flow_data)

    vf0, pm0 = vf0_pm0

    gen_init_data = (vf0, pm0, gen_state_algb_var, gen_dict_state_syms)
    
    # loc_load
    loc_load_state_algb_var, loc_load_f_t =
        self_init_state(
            loc_load_type, gen_init_data,
            loc_load_para_values )    
    
    return ( x0 = [gen_state_algb_var...;  loc_load_state_algb_var... ], f_t = (; gen_f_t, vf0, pm0, loc_load_f_t), aux = [vf0, pm0])
    
end



#--------------------------------------------------

function f_t_self_SC_2axis_wt_loc_load(power_flow_data)

    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]

    return [i_r, i_i ]
    
end


function f_t_self(component::Union{ plant_no_gov_wt_loc_load, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6 },  power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]    

    u = vh * exp( im * vh )
    
    return [i_r, i_i ], f_t_self_SC_2axis_wt_loc_load
    
    # return [ig_r, ig_i ], f_t_self_SM_2axis_wt_loc_load
    
end



function self_init_state(node::Union{plant_no_gov_wt_loc_load, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6}, power_flow_data)

    gen_dict_state_syms = node.Gen.dict_state_syms

    plant_dict_state_syms = node.dict_state_syms
    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 = self_init_state( node.Gen, power_flow_data)

    vf0, pm0 = vf0_pm0

    gen_init_data = (vf0, pm0, gen_state_algb_var, gen_dict_state_syms)
    
    # Exc
    exc_state_algb_var, avr_f_t = self_init_state(node.Exc, gen_init_data )

    # loc_load
    loc_load_state_algb_var, loc_load_f_t = self_init_state(node.Loc_load, gen_init_data )    
    
    return ( x0 = [gen_state_algb_var...; exc_state_algb_var...; loc_load_state_algb_var... ], f_t = (; gen_f_t, avr_f_t, loc_load_f_t), aux = [vf0, pm0])
    
end




function self_init_state(
    node_type::Union{Type{plant_no_gov_wt_loc_load},
                     Type{plant_no_gov_wt_loc_load_idq},
                     Type{plant_no_gov_wt_loc_load_v6} },
    gen_type,
    gen_para_values,
    gen_dim,
    gen_dict_state_syms,
    avr_type,
    avr_para_values,
    loc_load_type,
    loc_load_para_values,
    power_flow_data
    
    )

    gen_dict_state_syms = node.Gen.dict_state_syms

    plant_dict_state_syms = node.dict_state_syms
    
    # Gen

    gen_state_algb_var,  gen_f_t, vf0_pm0 =
        self_init_state(
            gen_type,            
            gen_para_values,
            gen_dim,
            gen_dict_state_syms,
            power_flow_data)
    
    vf0, pm0 = vf0_pm0

    gen_init_data =
        (vf0, pm0,
         gen_state_algb_var,
         gen_dict_state_syms)
    
    # Exc
    exc_state_algb_var, avr_f_t =
        self_init_state(
            avr_type,    
            gen_init_data,
            avr_para_values )

    # loc_load
    loc_load_state_algb_var, loc_load_f_t =
        self_init_state(
            loc_load_type,    
            gen_init_data,
            loc_load_para_values )    
    
    return (
        x0 = [gen_state_algb_var...;
              exc_state_algb_var...;
              loc_load_state_algb_var... ],
        f_t = (; gen_f_t, avr_f_t, loc_load_f_t),
        aux = [vf0, pm0])
    
end



#--------------------------------------------------

function f_t_self(component::Union{ plant_wt_loc_load, plant_wt_loc_load_idq, plant_wt_loc_load_v6 },  power_flow_data)
    
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    u = vh * exp( im * vh )
    
    return [i_r, i_i ], f_t_self_SM_2axis_wt_loc_load
    
    # return [ig_r, ig_i ], f_t_self_SM_2axis_wt_loc_load
    
end



function self_init_state(node::Union{ plant_wt_loc_load, plant_wt_loc_load_idq, plant_wt_loc_load_v6 }, power_flow_data)

    gen_dict_state_syms = node.Gen.dict_state_syms

    plant_dict_state_syms = node.dict_state_syms
    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 = self_init_state( node.Gen, power_flow_data)

    vf0, pm0 = vf0_pm0

    gen_init_data = (vf0, pm0, gen_state_algb_var, gen_dict_state_syms)
    
    # Gov
    gov_state_algb_var, gov_f_t = self_init_state(node.Gov, gen_init_data)

    # Exc
    exc_state_algb_var, avr_f_t = self_init_state(node.Exc, gen_init_data )

    # loc_load
    loc_load_state_algb_var, loc_load_f_t = self_init_state(node.Loc_load, gen_init_data )    
    
    return ( x0 = [gen_state_algb_var...; gov_state_algb_var...; exc_state_algb_var...; loc_load_state_algb_var... ], f_t = (; gen_f_t, gov_f_t, avr_f_t, loc_load_f_t), aux = [vf0, pm0])
    
end



function self_init_state(
    node_type::Union{ Type{plant_wt_loc_load},
                      Type{plant_wt_loc_load_idq},
                      Type{plant_wt_loc_load_v6} },    
    gen_type,
    gen_para_values,
    gen_dim,
    gen_dict_state_syms,
    gov_type,
    gov_para_values,    
    avr_type,
    avr_para_values,
    loc_load_type,
    loc_load_para_values,
    power_flow_data
    )

    # gen_dict_state_syms = node.Gen.dict_state_syms

    # plant_dict_state_syms = node.dict_state_syms
    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 =
        self_init_state(
            gen_type,            
            gen_para_values,
            gen_dim,
            gen_dict_state_syms,
            power_flow_data)

    vf0, pm0 = vf0_pm0

    gen_init_data =
        (vf0, pm0, gen_state_algb_var,
         gen_dict_state_syms)
    
    # Gov
    gov_state_algb_var, gov_f_t =
        self_init_state(
            gov_type,
            gen_init_data,
            gov_para_values )

    # Exc
    exc_state_algb_var, avr_f_t =
        self_init_state(
            avr_type,
            gen_init_data,
            avr_para_values )

    # loc_load
    loc_load_state_algb_var, loc_load_f_t =
        self_init_state(
            loc_load_type,    
            gen_init_data,
            loc_load_para_values )    
    
    return ( x0 = [gen_state_algb_var...; gov_state_algb_var...; exc_state_algb_var...; loc_load_state_algb_var... ], f_t = (; gen_f_t, gov_f_t, avr_f_t, loc_load_f_t), aux = [vf0, pm0])
    
end


#--------------------------------------------------

function f_t_self(component::Union{ plant_millano, plant_rscad_idq, plant_rscad_v6, plant_rscad, plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb, plant_SM_2axis_cb_system_matrix },  power_flow_data)

    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    return [i_r, i_i], f_t_self_SM_2axis
    
end



function self_init_state(node::Union{ plant_millano, plant_rscad_idq, plant_rscad_v6, plant_rscad, plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb, plant_SM_2axis_cb_system_matrix }, power_flow_data)

    gen_dict_state_syms = node.Gen.dict_state_syms

    plant_dict_state_syms = node.dict_state_syms
    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 = self_init_state( node.Gen, power_flow_data)

    vf0, pm0 = vf0_pm0

    gen_init_data = (vf0, pm0, gen_state_algb_var, gen_dict_state_syms)
    
    # Gov
    gov_state_algb_var, gov_f_t = self_init_state(node.Gov, gen_init_data)

    # Exc
    exc_state_algb_var, avr_f_t = self_init_state(node.Exc, gen_init_data )
    
    return ( x0 = [gen_state_algb_var...;gov_state_algb_var...;exc_state_algb_var...], f_t = (; gen_f_t, gov_f_t, avr_f_t), aux = [vf0, pm0])
    
end




function self_init_state(
    node_type::Union{Type{plant_millano},
                     Type{plant_rscad_idq},
                     Type{plant_rscad_v6},
                     Type{plant_rscad},
                     Type{plant_cb_idq},
                     Type{plant_cb_v6},
                     Type{plant_cb_direct},
                     Type{plant_cb_millano},
                     Type{plant_cb},
                     Type{plant_SM_2axis_cb_system_matrix}},
    
    gen_type,
    gen_para_values,
    gen_dim,
    gen_dict_state_syms,
    gov_type,
    gov_para_values,    
    avr_type,
    avr_para_values,
    power_flow_data
    )
    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 =
        self_init_state(
            gen_type,            
            gen_para_values,
            gen_dim,
            gen_dict_state_syms,
            power_flow_data)

    vf0, pm0 = vf0_pm0

    gen_init_data = (vf0, pm0, gen_state_algb_var, gen_dict_state_syms)
    
    # Gov
    gov_state_algb_var, gov_f_t =
        self_init_state(
            gov_type,
            gen_init_data,
            gov_para_values )
    

    # Exc
    exc_state_algb_var, avr_f_t =
        self_init_state(
            avr_type,    
            gen_init_data,
            avr_para_values )
    
    return ( x0 = [gen_state_algb_var...;gov_state_algb_var...;exc_state_algb_var...], f_t = (; gen_f_t, gov_f_t, avr_f_t), aux = [vf0, pm0])
    
end



function f_t_self(component::Union{plant_cb_inf, plant_cb_inf_bus },  power_flow_data)

    vh  = power_flow_data[1]
    θh  = power_flow_data[2]    
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    return [i_r, i_i, vh, θh], f_t_self_SM_cb_inf
    
end


function self_init_state(node::Union{ plant_cb_inf, plant_cb_inf_bus }, power_flow_data)

    gen_dict_state_syms = node.Gen.dict_state_syms
    
    # Gen
    gen_state_var_algebraic_var, gen_f_t, vf0_pm0_vh_θh = self_init_state( node.Gen, power_flow_data)

    vf0, pm0, vh, θh = vf0_pm0_vh_θh 

    gen_init_data = (vf0, pm0, gen_state_var_algebraic_var, gen_dict_state_syms)
    
    # Gov
    gov_state_var_algebraic_var, gov_f_t = self_init_state( node.Gov, gen_init_data)

    # Exc
    exc_state_var_algebraic_var, avr_f_t = self_init_state(  node.Exc, gen_init_data )


    return (x0 = [gen_state_var_algebraic_var...;gov_state_var_algebraic_var...;exc_state_var_algebraic_var...], f_t =  (; gen_f_t, gov_f_t, avr_f_t), aux = [vf0, pm0, vh, θh])
    
end


function f_t_self(component::Union{ plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer, plant_pss_cb },  power_flow_data)

   
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]    

    return [i_r, i_i], f_t_self_SM_2axis
    
end



function self_init_state(node::Union{ plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer, plant_pss_cb }, power_flow_data)


    gen_dict_state_syms = node.Gen.dict_state_syms
    
    # Gen
    
    gen_state_var_algebraic_var,  gen_f_t, vf0_pm0 = self_init_state( node.Gen, power_flow_data)

    vf0, pm0 = vf0_pm0    

     gen_init_data = (vf0, pm0, gen_state_var_algebraic_var, gen_dict_state_syms)    

    # Gov
    gov_state_var_algebraic_var, gov_f_t = self_init_state( node.Gov, gen_init_data)

    # Exc
    
    exc_state_var_algebraic_var, avr_f_t = self_init_state( node.Exc, gen_init_data )

    # Pss
    
    pss_state_var_algebraic_var, pss_f_t = self_init_state( node.Pss, gen_init_data )
    
    return (x0 =  [gen_state_var_algebraic_var...;gov_state_var_algebraic_var...;exc_state_var_algebraic_var...;pss_state_var_algebraic_var...], f_t =  (; gen_f_t, gov_f_t, avr_f_t, pss_f_t ),  aux =[vf0, pm0])
    
    
end



function self_init_state(
    node_type::Union{ Type{plant_pss_cb_idq},
                      Type{plant_pss_cb_v6},
                      Type{plant_pss_cb_direct},
                      Type{plant_pss_cb_sauer},
                      Type{plant_pss_cb} },    
    gen_type,
    gen_para_values,
    gen_dim,
    gen_dict_state_syms,
    gov_type,
    gov_para_values,    
    avr_type,
    avr_para_values,
    pss_type,
    pss_para_values,
    power_flow_data
    )

    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 =
        self_init_state(
            gen_type,            
            gen_para_values,
            gen_dim,
            gen_dict_state_syms,
            power_flow_data)

    vf0, pm0 = vf0_pm0

    gen_init_data = (vf0, pm0, gen_state_algb_var,
                     gen_dict_state_syms)
    
    # Gov
    gov_state_algb_var, gov_f_t =
        self_init_state(
            gov_type,
            gen_init_data,
            gov_para_values )
    

    # Exc
    exc_state_algb_var, avr_f_t =
        self_init_state(
            avr_type,    
            gen_init_data,
            avr_para_values )
    

    # Pss    
    pss_state_var_algebraic_var, pss_f_t =
        self_init_state(
            pss_type,
            gen_init_data,
            pss_para_values )
    
    return (x0 =  [gen_state_var_algebraic_var...;gov_state_var_algebraic_var...;exc_state_var_algebraic_var...;pss_state_var_algebraic_var...], f_t =  (; gen_f_t, gov_f_t, avr_f_t, pss_f_t ),  aux =[vf0, pm0])
    
    
end



function f_t_self(component::Union{plant_pss_cb_inf, plant_pss_cb_inf_bus },  power_flow_data)

   
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    i_r = power_flow_data[5]
    i_i = power_flow_data[6]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    return [i_r, i_i, vh, θh], f_t_self_SM_cb_inf

    
end



function self_init_state(node::Union{ plant_pss_cb_inf, plant_pss_cb_inf_bus }, power_flow_data)

    gen_dict_state_syms = node.Gen.dict_state_syms
    
    # Gen
    
    gen_state_var_algebraic_var, gen_f_t, vf0_pm0_vh_θh  = self_init_state( node.Gen, power_flow_data)


    vf0, pm0, vh, θh = vf0_pm0_vh_θh     

     gen_init_data = (vf0, pm0, gen_state_var_algebraic_var, gen_dict_state_syms)    

    # Gov
    gov_state_var_algebraic_var, gov_f_t = self_init_state( node.Gov, gen_init_data)

    # Exc
    
    exc_state_var_algebraic_var, avr_f_t = self_init_state( node.Exc, gen_init_data )

    # Pss
    
    pss_state_var_algebraic_var, pss_f_t = self_init_state( node.Pss, gen_init_data )

    
    return ( x0 = [gen_state_var_algebraic_var...;gov_state_var_algebraic_var...;exc_state_var_algebraic_var...;pss_state_var_algebraic_var...], f_t = (; gen_f_t, gov_f_t, avr_f_t, pss_f_t), aux =  [vf0, pm0, vh, θh])
    
    
end



function f_t_self(component::Union{ plant_Infinite_bus, plant_Infinite_cb_bus },  power_flow_data)

   
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
    
    return [i_r, i_i, vh, θh], f_t_self_SM_cb_inf

    
end



function self_init_state(node::Union{ plant_Infinite_bus, plant_Infinite_cb_bus }, power_flow_data)

    gen_dict_state_syms = node.Gen.dict_state_syms
    
    # Gen
        
    gen_state_var_algebraic_var, gen_f_t, vh_θh_ph_qh = self_init_state( node.Gen, power_flow_data)


    return ( x0 = [gen_state_var_algebraic_var...],  f_t = (; gen_f_t, ), aux = [vh, θh, ph, qh])
    
end


# ------------------------------------------------------
# Init for Edges
# ------------------------------------------------------


function self_init_state(branch::PiModelLine, power_flow_data)

    if_r = power_flow_data[1]
    if_i = power_flow_data[2]
    it_r = power_flow_data[3]
    it_i = power_flow_data[4]

    
    # p = get_component_params_value(node)
    
    return (x0 = [if_r, if_i, it_r, it_i], f_t =[0.0, 0.0])
    
    
end



function self_init_state(branch::StaticLine,
                         power_flow_data )

    if_r = power_flow_data[1]
    if_i = power_flow_data[2]
    it_r = power_flow_data[3]
    it_i = power_flow_data[4]
    
    # p = get_component_params_value(node)
    
    return (x0 = [if_r, if_i, it_r, it_i], f_t =[0.0, 0.0] )
    
    
end


function self_init_state(branch::RLLine, power_flow_data)

    if_r = power_flow_data[1]
    if_i = power_flow_data[2]
    it_r = power_flow_data[3]
    it_i = power_flow_data[4]
    
    # p = get_component_params_value(node)
    
    return (x0 = [if_r, if_i, it_r, it_i], f_t =[0.0, 0.0])
    
    
end


function self_init_state(branch::Transformer, power_flow_data)
    
    if_r = power_flow_data[1]
    if_i = power_flow_data[2]
    it_r = power_flow_data[3]
    it_i = power_flow_data[4]
       
    return (x0 = [if_r, if_i, it_r, it_i], f_t =[0.0, 0.0])
    
    
end


function self_init_state(
    branch_type::Union{ Type{PiModelLine},
                        Type{StaticLine},
                        Type{RLLine},
                        Type{Transformer} },
    power_flow_data )

    if_r = power_flow_data[1]
    if_i = power_flow_data[2]
    it_r = power_flow_data[3]
    it_i = power_flow_data[4]
    
    # p = get_component_params_value(node)
    
    return (x0 = [if_r, if_i, it_r, it_i], f_t =[0.0, 0.0])
    
    
end


#########################################################

function self_x0(component,  power_flow_data)

    x0_f_t_aux = self_init_state(component, power_flow_data)
       
    return x0_f_t_aux.x0 
    
end

function self_f_t(component,  power_flow_data)

    x0_f_t_aux = self_init_state(component, power_flow_data)
       
    return x0_f_t_aux.f_t 
    
end

function self_aux(component,  power_flow_data)

    x0_f_t_aux = self_init_state(component, power_flow_data)
       
    return x0_f_t_aux.aux 
    
end

function self_x0_f_t_aux(component,  power_flow_data)
       
    return self_init_state(component, power_flow_data) 
    
end

#########################################################
# plant cb_sw
#########################################################



function self_node_cb_sw(node::Union{ plant_no_gov_no_avr_wt_loc_load } )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )

    return node_cb_sw
    
end




function self_node_cb_sw(node::Union{ plant_no_gov_wt_loc_load, plant_no_gov_wt_loc_load_idq,  plant_no_gov_wt_loc_load_v6 } )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )
    
    # cb_sw = get_component_cb_dyn_param_state_sw( node )

    # gen      = node.Gen
    # avr      = node.Exc
    # loc_load = node.Loc_load    

    # if gen.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0
    #     cb_sw_avr      = cb_sw
    #     cb_sw_gen      = Int64[]
    #     cb_sw_loc_load = Int64[]

    # else

    #     dims_cb_sw = [gen.cb_dyn_state_dim,
    #                   avr.cb_dyn_state_dim,
    #                   loc_load.cb_dyn_state_dim]

    #     cb_sw_offset = create_offsets(dims_cb_sw)
    #     cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

    #     cb_sw_gen      = cb_sw[cb_sw_Idx[1]]
    #     cb_sw_avr      = cb_sw[cb_sw_Idx[2]]
    #     cb_sw_loc_load = cb_sw[cb_sw_Idx[3]]
    # end
    
    # node_cb_sw = (cb_sw_gen, cb_sw_avr, cb_sw_loc_load)
    
    return node_cb_sw
    
end



function self_node_cb_sw(node::Union{ plant_wt_loc_load, plant_wt_loc_load_idq, plant_wt_loc_load_v6 } )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )

    # cb_sw = get_component_cb_dyn_param_state_sw(node)

    # gen      = node.Gen
    # gov      = node.Gov
    # avr      = node.Exc
    # loc_load = node.Loc_load

    # if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0

    #     cb_sw_gen      = Int64[]
    #     cb_sw_gov      = Int64[]
    #     cb_sw_loc_load = Int64[]
        
    #     cb_sw_avr      = cb_sw
        

    # elseif gen.cb_dyn_state_dim == 0  && loc_load.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

    #     cb_sw_gen      = Int64[]
    #     cb_sw_loc_load = Int64[]
        
    #     dims_cb_sw = [gov.cb_dyn_state_dim,
    #                   avr.cb_dyn_state_dim]

    #     cb_sw_offset = create_offsets(dims_cb_sw)
    #     cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

    #     cb_sw_gov      = cb_sw[cb_sw_Idx[1]]
    #     cb_sw_avr      = cb_sw[cb_sw_Idx[2]]
        
    # else
        
    #     dims_cb_sw = [gen.cb_dyn_state_dim,
    #                   gov.cb_dyn_state_dim,
    #                   avr.cb_dyn_state_dim,
    #                   loc_load.cb_dyn_state_dim]

    #     cb_sw_offset = create_offsets(dims_cb_sw)
    #     cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

    #     cb_sw_gen      = cb_sw[cb_sw_Idx[1]]
    #     cb_sw_gov      = cb_sw[cb_sw_Idx[2]]
    #     cb_sw_avr      = cb_sw[cb_sw_Idx[3]]
    #     cb_sw_loc_load = cb_sw[cb_sw_Idx[4]]        
    # end
    
    # node_cb_sw = (cb_sw_gen, cb_sw_gov, cb_sw_avr, cb_sw_loc_load)
    
    return node_cb_sw
    
end


function self_node_cb_sw(node::Union{ pf_plant_SM_v6, plant_SM_idq, plant_SM_v6, plant_SM_system_matrix  } )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )
    
    return node_cb_sw
    
end


function self_node_cb_sw(node::Union{ plant_no_gov_system_matrix } )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )
    
    return node_cb_sw
    
end



function self_node_cb_sw(node::Union{ plant_millano, plant_rscad_idq, plant_rscad_v6, plant_rscad, plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb_inf,  plant_cb_inf_bus, plant_cb, plant_SM_2axis_cb_system_matrix } )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )
    
    # cb_sw = get_component_cb_dyn_param_state_sw( node )

    # gen      = node.Gen
    # gov      = node.Gov
    # avr      = node.Exc

    # if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 

    #     cb_sw_gen      = Int64[]
    #     cb_sw_gov      = Int64[]
        
    #     cb_sw_avr      = cb_sw
        

    # elseif gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

    #     cb_sw_gen      = Int64[]
        
    #     dims_cb_sw = [gov.cb_dyn_state_dim,
    #                   avr.cb_dyn_state_dim]

    #     cb_sw_offset = create_offsets(dims_cb_sw)
    #     cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

    #     cb_sw_gov      = cb_sw[cb_sw_Idx[1]]
    #     cb_sw_avr      = cb_sw[cb_sw_Idx[2]]
        
    # else
        
    #     dims_cb_sw = [gen.cb_dyn_state_dim,
    #                   gov.cb_dyn_state_dim,
    #                   avr.cb_dyn_state_dim]

    #     cb_sw_offset = create_offsets(dims_cb_sw)
    #     cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

    #     cb_sw_gen      = cb_sw[cb_sw_Idx[1]]
    #     cb_sw_gov      = cb_sw[cb_sw_Idx[2]]
    #     cb_sw_avr      = cb_sw[cb_sw_Idx[3]]
    # end
    
    # node_cb_sw = (cb_sw_gen, cb_sw_gov, cb_sw_avr)
    
    return node_cb_sw
    
end

function self_node_cb_sw(node::Union{ plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_pss_cb} )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )
    
    return node_cb_sw
    
end

function self_node_cb_sw(node::Union{ plant_Infinite_bus, plant_Infinite_cb_bus } )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )
    
    return node_cb_sw
    
end

function self_node_cb_sw(node::Union{ plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load  } )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )
    
    return node_cb_sw
    
end

function self_node_cb_sw(node::Union{ plant_Transmission_t1, plant_Transmission_t2 } )

    node_cb_sw = get_component_cb_dyn_param_state_sw( node )
    
    return node_cb_sw
    
end

########################################################
# get nodes cb_sw
########################################################

function  get_nodes_cb_sw(comps_collection)
    
    return  map(get_component_cb_dyn_param_state_sw, collect(values(comps_collection)))
      
end

########################################################
# get nodes  pm_vf_x
#########################################################

"""
This is a generic function for making a tuple of components and values:

get_tuple_comps_pf_init(comps_collection, pf_init_dict, key="" )

is equivalent to :

get_tuple_comps_any_dict(comps_collection, pf_init_dict, key="" )

"""
function  get_tuple_comps_any_dict(comps_collection::OrderedDict, any_dict::Dict; key="" )

    dict_values = any_dict[key]
    
    comps_keys   = collect(keys(comps_collection))
    # comps_names   = map(get_component_name, collect(values(comps_collection)))

    tuple_vec_comps_values = [(comps_collection[a_key], dict_values[a_key]) for a_key in comps_keys ]    
    
    return  tuple_vec_comps_values
      
end


function  get_tuple_comps_pf_init(comps_collection, pf_init_dict; key="" )

    dict_pf_init = pf_init_dict[key]
    
    comps_keys   = collect(keys(comps_collection))

    tuple_vec_comps_pf_init = [(comps_collection[a_key], dict_pf_init[a_key]) for a_key in comps_keys ]    
    
    return  tuple_vec_comps_pf_init
      
end



"""
This is a generic input function

get_edges_pm_vf_x(comps_collection, pf_init_dict; key="branch_dict_init" )

is equivalent to:

get_comps_input(self_edge_pm_vf_x, comps_collection,  pf_init_dict; key="branch_dict_init" )

similarly

get_nodes_pm_vf_x(comps_collection, pf_init_dict; key="bus_dict_init")

is equvalent to:

get_comps_input(self_node_pm_vf_x, comps_collection,  pf_init_dict; key="bus_dict_init" )

To get net current injection use:

get_comp_input(self_node_pf_Inet_inj, comps_collection, pf_init_dict; key="bus_dict_init"  )



"""
function  get_comps_input(input_func, comps_collection, any_dict; key="bus_dict_init"  )

    tuple_vec_comps_values = get_tuple_comps_any_dict(comps_collection, any_dict; key=key )
    
    vec_comps_inputs = map((arg) -> input_func(arg[1], arg[2]), tuple_vec_comps_values )
    
    return  vec_comps_inputs 
      
end


#######################################################
# get nodes x0, f_t, aux, cb_sw
########################################################


function get_componets_f_t_self( comps_collection, pf_f_t_dict;  key = "bus_dict_init" )

    tuple_vec_comps_f_t = get_tuple_comps_pf_init(comps_collection,  pf_f_t_dict; key =  key )

    return first.(map((arg) -> f_t_self(arg[1], arg[2]), tuple_vec_comps_f_t ))
    
end


function get_componets_f_t_self_func( comps_collection, pf_f_t_dict;  key = "bus_dict_init" )

    tuple_vec_comps_f_t = get_tuple_comps_pf_init(comps_collection,  pf_f_t_dict; key = key )

    return last.(map((arg) -> f_t_self(arg[1], arg[2]), tuple_vec_comps_f_t ))
    
end


function get_nodes_f_t_self_and_func( comps_collection, pf_f_t_dict; key = "bus_dict_init" )

    tuple_vec_comps_f_t = get_tuple_comps_pf_init(comps_collection,  pf_f_t_dict; key = key )

    return map((arg) -> f_t_self(arg[1], arg[2]), tuple_vec_comps_f_t )
    
end



function get_nodes_pf_f_t(comps_nodes_f_t_func, tuple_pf_result)

    return Vector{Float64}[f_t_self_func( pf_result) for (f_t_self_func, pf_result) in zip(comps_nodes_f_t_func, tuple_pf_result)]    
    
end


function  get_components_x0( comps_collection,
    pf_init_dict; key = "" )

    tuple_vec_comps_pf_init = get_tuple_comps_pf_init(comps_collection, pf_init_dict; key = key )

    vec_comps_x0 = map((arg) -> self_x0(arg[1], arg[2]), tuple_vec_comps_pf_init)
    
    return  vcat(vec_comps_x0...)
      
end


function  get_components_f_t(comps_collection, pf_init_dict; key = "" )

    tuple_vec_comps_pf_init = get_tuple_comps_pf_init(comps_collection, pf_init_dict; key = key )

    vec_comps_f_t = map((arg) -> self_f_t(arg[1], arg[2]), tuple_vec_comps_pf_init)
    
    return  vec_comps_f_t
      
end



function  get_components_aux(comps_collection, pf_init_dict; key = "" )

    tuple_vec_comps_pf_init = get_tuple_comps_pf_init(comps_collection, pf_init_dict; key = key )

    vec_comps_aux = map((arg) -> self_aux(arg[1], arg[2]), tuple_vec_comps_pf_init)
    
    return  vec_comps_aux
      
end


function  get_components_x0_f_t_aux(comps_collection, pf_init_dict; key = "" )

    tuple_vec_comps_pf_init = get_tuple_comps_pf_init(comps_collection, pf_init_dict; key = key )

    vec_comps_x0_f_t_aux = map((arg) -> self_x0_f_t_aux(arg[1], arg[2]), tuple_vec_comps_pf_init)
    
    return  vec_comps_x0_f_t_aux
      
end

function  get_nodes_x0(nodes_collection, pf_init_dict; key = "bus_dict_init" )
    
    return get_components_x0(nodes_collection, pf_init_dict; key = key )
      
end

function  get_nodes_f_t(nodes_collection, pf_init_dict; key = "bus_dict_init" )

    return get_components_f_t(nodes_collection, pf_init_dict; key = key )
    
      
end



function  get_nodes_aux(nodes_collection, pf_init_dict; key = "bus_dict_init" )

    return get_components_aux(nodes_collection, pf_init_dict; key = key )
      
end

function  get_nodes_x0_f_t_aux(nodes_collection, pf_init_dict; key = "bus_dict_init" )

    return get_components_x0_f_t_aux(nodes_collection, pf_init_dict; key = key )
      
end

function  get_edges_x0(edges_collection, pf_init_dict; key = "branch_dict_init" )
    
    return  get_components_x0(edges_collection, pf_init_dict; key = key )
      
end


function  get_edges_f_t(edges_collection, pf_init_dict; key="branch_dict_init" )

    return get_components_f_t(edges_collection, pf_init_dict; key = key )
      
end


# ------------------------------------------------------
# ------------------------------------------------------


function  external_get_tuple_comps_pf_init(comps_collection,  dict_pf_init )

    
    comps_keys   = collect(keys(comps_collection))

    tuple_vec_comps_pf_init = [(comps_collection[a_key], dict_pf_init[a_key]) for a_key in comps_keys ]    
    
    return  tuple_vec_comps_pf_init
      
end



function  external_get_components_f_t(comps_collection, dict_pf_init )

    tuple_vec_comps_pf_init = external_get_tuple_comps_pf_init(comps_collection, dict_pf_init  )

    vec_comps_f_t = map((arg) -> self_f_t(arg[1], arg[2]), tuple_vec_comps_pf_init)
    
    return  vec_comps_f_t
      
end


function  external_get_nodes_or_edges_f_t(collection, dict_init )


    return external_get_components_f_t(collection, dict_init )
    
      
end



function  external_get_components_x0(comps_collection, dict_pf_init )

    tuple_vec_comps_pf_init = external_get_tuple_comps_pf_init(comps_collection, dict_pf_init  )

    vec_comps_x0 = map((arg) -> self_x0(arg[1], arg[2]), tuple_vec_comps_pf_init)
    
    return  vcat(vec_comps_x0...)
      
end


function  external_get_nodes_or_edges_x0(collection, dict_init  )
    
    return external_get_components_x0(collection, dict_init )
      
end


# ------------------------------------------------------
#  external get_state_init
# ------------------------------------------------------


function external_init_operationpoint(nd::NetworkData, bus_dict_init, branch_dict_init )

    nodes_state_init = external_get_nodes_or_edges_x0(nd.nodes, bus_dict_init)

    edges_state_init = external_get_nodes_or_edges_x0(nd.edges, branch_dict_init)

    state_init = vcat(nodes_state_init, edges_state_init)
    
    return state_init 
end


# ------------------------------------------------------
#   industrial model get_state_init
# ------------------------------------------------------


function  industrial_model_get_components_x0(comps_collection, dict_pf_init, netd::NetworkData; pure = :pure,  no_control_device = false )

    if  no_control_device == false
        tuple_vec_comps_pf_init = external_get_tuple_comps_pf_init(comps_collection, dict_pf_init  )

        vec_comps_x0 = map((arg) -> self_x0(arg[1], arg[2]), tuple_vec_comps_pf_init)

        x0 = vcat(vec_comps_x0...)

        if  pure == :pure

            return  vcat(x0[[ get_gens_pure_states_indices_in_system( netd )...; ] ], x0[ [ get_nodes_ur_ui_indices_in_system(netd )...; ]  ])
        elseif pure == :stab
            return vcat(x0[ [ get_gens_stab_states_indices_in_system( netd )...; ] ], x0[[ get_nodes_ur_ui_indices_in_system(netd )...; ] ])
        else
            return nothing
        end
    else
        tuple_vec_comps_pf_init = external_get_tuple_comps_pf_init(comps_collection, dict_pf_init  )

        vec_comps_x0 = map((arg) -> self_x0(arg[1], arg[2]), tuple_vec_comps_pf_init)

        x0 = vcat(vec_comps_x0...)

        if  pure == :pure

            return  vcat(x0[[ get_gens_pure_states_indices_in_system( netd; no_control_device = true  )...; ] ], x0[ [ get_nodes_ur_ui_indices_in_system(netd  )...; ]  ])
        elseif pure == :stab
            return vcat(x0[ [ get_gens_stab_states_indices_in_system( netd; no_control_device = true  )...; ] ], x0[[ get_nodes_ur_ui_indices_in_system(netd )...; ] ])
        else
            return nothing
        end

    end
    
      
end



function  industrial_model_get_nodes_x0(
    collection,
    dict_init,
    netd::NetworkData;
    pure = :pure, no_control_device = false )

    if no_control_device == false

        return industrial_model_get_components_x0(collection, dict_init, netd ; pure = pure )
    else
        return industrial_model_get_components_x0(collection, dict_init, netd ; pure = pure, no_control_device = true )
    end
    
      
end


function industrial_model_init_operationpoint(
    nd::NetworkData,
    bus_dict_init
    ;pure = :pure, no_control_device = false )

    if no_control_device == false

        return industrial_model_get_nodes_x0(nd.nodes, bus_dict_init, nd ; pure = pure )
        
    else
        
        return industrial_model_get_nodes_x0(nd.nodes, bus_dict_init, nd ;pure = pure, no_control_device = true )
        
    end
    
end



# ------------------------------------------------------
#   im get_state_init
# ------------------------------------------------------


function  im_model_get_components_x0(comps_collection, dict_pf_init, netd::NetworkData )

    tuple_vec_comps_pf_init = external_get_tuple_comps_pf_init(comps_collection, dict_pf_init  )

    vec_comps_x0 = map((arg) -> self_x0(arg[1], arg[2]), tuple_vec_comps_pf_init)

    x0 = vcat(vec_comps_x0...)

    return  vcat(x0[[ get_gens_im_vars_indices_in_system( netd )...; ] ], x0[ [ get_nodes_ur_ui_indices_in_system(netd )...; ]  ])
      
end



function  im_model_get_nodes_x0(
    collection,
    dict_init,
    netd::NetworkData )

    return im_model_get_components_x0(collection, dict_init, netd )
      
end


function im_model_init_operationpoint(
    nd::NetworkData,
    bus_dict_init )

return im_model_get_nodes_x0(nd.nodes, bus_dict_init, nd  )    
end



# ------------------------------------------------------
#  get_state_init
# ------------------------------------------------------


function init_operationpoint(nd::NetworkData, pf_init_dict)

    nodes_state_init = get_nodes_x0(nd.nodes, pf_init_dict; key = "bus_dict_init")

    edges_state_init = get_edges_x0(nd.edges, pf_init_dict; key = "branch_dict_init")

    state_init = vcat(nodes_state_init, edges_state_init)
    
    return state_init 
end


# ------------------------------------------------------
# Init for Nodes and Edges alternative
# ------------------------------------------------------


function pf_f_t_self_PQ_Const_P(uh, src_i, dst_i)

    vh, θh = abs(uh), angle(uh)    

    v   = uh      = vh * exp(im * θh)
    
   i  = 1.0 * pf_dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    i_r      = real(i)

    i_i      = imag(i)
        
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )

    return [i_r, i_i, vh, θh]
    
end


function pf_f_t_self(component::PQ_Const_P,  uh, src_i, dst_i)
    
    # uh, src_i, dst_i

    # p = component.param_values

    # P   = p[1]
    # Q   = p[2]
    # kPL = p[3]
    # kQL = p[4]

    # Y_n = p[5]

    # y_shunt = im * Y_n

    y_shunt = im * component.Y_n
    

    vh, θh = abs(uh), angle(uh)

    v = uh   = vh * exp(im * θh)
    
   i  = 1.0 * pf_dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    i_r      = real(i)

    i_i      = imag(i)

    
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )

    return [i_r, i_i, vh, θh], f_t_self_PQ_Const_P
    
end


function pf_self_init_state(component::PQ_Const_P, uh, src_i, dst_i)

    # p = component.param_values

    # P   = p[1]
    # Q   = p[2]
    # kPL = p[3]
    # kQL = p[4]

    # Y_n = p[5]

    y_shunt = im * Y_n

    y_shunt = im * component.Y_n
    
    vh, θh = abs(uh), angle(uh)
    
    v = uh   = vh * exp(im * θh)
    
    i  = 1.0 * pf_dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    i_r      = real(i)

    i_i      = imag(i)
    
    # p = get_component_params_value(component)

    
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )
        
    state_var     = [real(v), imag(v) ]
    algebraic_var = []

    return  (x0 = [state_var...;algebraic_var...], f_t = [i_r, i_i, vh, θh],  aux = [kPL, kQL] )
    
    # return  (state_var, algebraic_var) , (ph, qh, kPL, kQL)
end


function pf_f_t_self(component::plant_PQ_Const_P,  uh, src_i, dst_i )

    # p = component.param_values

    # P   = p[1]
    # Q   = p[2]
    # kPL = p[3]
    # kQL = p[4]

    # Y_n = p[5]

    # y_shunt = im * Y_n

    y_shunt = im * component.Y_n
    
    vh, θh = abs(uh), angle(uh)    
    
    i  = 1.0 * pf_dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    i_r      = real(i)

    i_i      = imag(i)

    v        = vh * exp(im * θh)
    
    kPL       = (1.0, 0.0, 0.0 )
    kQL       = (1.0, 0.0, 0.0 )

    return [i_r, i_i, vh, θh], f_t_self_PQ_Const_P
    
end


function pf_self_init_state(node::Union{ plant_PQ_Const_P}, uh, src_i, dst_i )

    load_dict_state_syms = node.Load.dict_state_syms
    
    # Load
        
    load_state_algb_var, load_f_t, load_aux = pf_self_init_state( node.Load, uh, src_i, dst_i )

    return ( x0 = [load_state_algb_var...], f_t = (; load_f_t,), aux = load_aux )
    
end


#-----------------------------------------------------


function pf_f_t_self_SM( uh, src_i, dst_i )
        
   i  = 1.0 * pf_dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    i_r      = real(i)

    i_i      = imag(i)

    return [i_r, i_i ]
    
end


function pf_f_t_self(component::Union{pf_SM_2axis_v6},  uh, src_i, dst_i )

    y_shunt = im * component.Y_n

   i0  = 1.0 * pf_dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    i_r      = real(i0)

    i_i      = imag(i0)
    

    return  [i_r, i_i ], f_t_self_SM 
end


function pf_self_init_state(component::Union{pf_SM_2axis_v6}, uh, src_i, dst_i )

    
    p = get_component_params_value(component)

    # Gen 20
    
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

    uh_mag    = p[22]

    y_shunt = im * Y_n

    vh, θh   = abs(uh), angle(uh)
    
    v        = vh * exp(im * θh)

    u_r      = real(v)
    
    u_i      = imag(v)
    
    i0  = 1.0 * pf_dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    i_r      = real(i0)

    i_i      = imag(i0)

    ih0      = abs(i0)
    
    ϕ0       = angle(i0)    
    
    E0       = v + (ra + im * X_q) * i0
    
    δ0       = angle(E0)
    
    v_dq     = v * exp(-im * (δ0 - pi/2))
    vd0      = real(v_dq)
    vq0      = imag(v_dq)

    i_dq     = i0 * exp(-im * (δ0 - pi/2))
    id0      = real(i_dq)
    iq0      = imag(i_dq)
    
    # vd0      = vh * sin(δ0 - θh)
    # vq0      = vh * cos(δ0 - θh)
    
    # id0      = ih0 * sin(δ0 - ϕ0)
    # iq0      = ih0 * cos(δ0 - ϕ0)

    ed0_dash0 =  (X_q - X_q_dash) * iq0
    
    # ed0_dash0 = vd0 + ra * id0 - X_q_dash  * iq0
    
    eq0_dash0 = vq0 + ra * iq0 + X_d_dash  * id0   

    # avr
    
    vf0 = eq0_dash0 + (X_d - X_d_dash) * id0
    
    # Gov
    
    # τe0 = (vq0 + ra * iq0) * iq0 + (vd0 + ra * id0) * id0
    
    # τe0 = ph0 + ra * (id0^2 + iq0^2)

    τe0 = ed0_dash0 * id0 + eq0_dash0 * iq0 + (X_q_dash - X_d_dash) * id0 * iq0

    ph = (ed0_dash0 - (ra * id0 - X_q_dash * iq0) ) * id0 + (eq0_dash0 - (X_d_dash * id0  + ra * iq0)) * iq0  
    
    qh = (eq0_dash0 - (X_d_dash * id0  + ra * iq0)) * id0 - (ed0_dash0 - (ra * id0 - X_q_dash * iq0) ) * iq0  
    
    ω0  = 1.0  # ω0 = ω_ref0 = ωs  

    pe0 = τe0 * ω0
    
    pm0 = pe0
        
    dim            = component.dim
    
    dict_state_syms = component.dict_state_syms

    state_algebraic_var = zeros(dim)    

    list_state_and_algb_val  = [ δ0, ω0, ed0_dash0, eq0_dash0, u_r, u_i, ph, qh, uh_mag ]
    
    list_state_and_algb_sym  = [ :δ, :ω, :ed_dash, :eq_dash, :u_r, :u_i, :ph, :qh, :uh_mag ]
    

    for (sym, val) in zip( list_state_and_algb_sym, list_state_and_algb_val)
        
        state_algebraic_var[ dict_state_syms[ sym ]] = val
        
    end

    return ( x0 = state_algebraic_var, f_t = [i_r, i_i], aux = [vf0, pm0])
    
end


function pf_f_t_self(component::Union{ pf_plant_SM_v6 },  uh, src_i, dst_i )

    Y_n = component.Gen.param_values[20]

    y_shunt = im * Y_n
    
    i  = 1.0 * pf_dynamic_nodal_current_balance(src_i, dst_i) + uh * y_shunt

    i_r      = real(i)

    i_i      = imag(i)
    

    return  [i_r, i_i ], f_t_self_SM 

    
end


function pf_self_init_state(node::Union{pf_plant_SM_v6 }, uh, src_i, dst_i )

    gen_dict_state_syms = node.Gen.dict_state_syms

    plant_dict_state_syms = node.dict_state_syms
    
    # Gen
    gen_state_algb_var,  gen_f_t, vf0_pm0 = pf_self_init_state( node.Gen, uh, src_i, dst_i)

    vf0, pm0 = vf0_pm0

    return ( x0 = [gen_state_algb_var...;], f_t = (; gen_f_t, vf0, pm0  ), aux = [vf0, pm0])
    
end


#----------------------------------------------------------


function pf_get_tuple_nodes_pf_init(nodes_collection, edges_collection, nodes_u )

    Cnb = get_Cnb( edges_collection )
    
    nodes_keys  = collect(keys(nodes_collection))

    nodes = collect(values( nodes_collection ))

    nodes_src_edges = [get_node_src_edges(node_idx, Cnb ) for node_idx in collect(1:length( nodes_keys))]
    
    nodes_dst_edges = [get_node_dst_edges(node_idx, Cnb ) for node_idx in collect(1:length( nodes_keys))]

    edges_src_dst_node_u = get_edges_src_dst_node_u(edges_collection, nodes_u)

    Ybr_cal = [ Symbol(typeof(a_branch )) == :PiModelLine ? calc_branch_Ybr(a_branch.param_values[1:3]..., 1, 1 ) : calc_branch_Ybr(a_branch.param_values[1:4]..., 1 ) for a_branch in collect(values(edges_collection)) ]

   edges_from_and_to_current_Ybr = [ Ybr * src_dst_node_u  for (Ybr, src_dst_node_u) in zip(Ybr_cal, edges_src_dst_node_u)]

    edges_from_current_Ybr = first.(edges_from_and_to_current_Ybr)
    
    edges_to_current_Ybr   = last.(edges_from_and_to_current_Ybr)

    # nodes_src_i = [ src_idx != [] ? edges_from_current_Ybr[src_idx] : [] for src_idx in nodes_src_edges ]

    nodes_src_i = [ src_idx != [] ? edges_from_current_Ybr[src_idx...] : [] for src_idx in nodes_src_edges ]

    # nodes_dst_i = [ dst_idx != [] ? edges_to_current_Ybr[dst_idx] : [] for dst_idx in nodes_dst_edges ]

    nodes_dst_i = [ dst_idx != [] ? edges_to_current_Ybr[dst_idx...] : [] for dst_idx in nodes_dst_edges ]

    
    tuple_vec_comps_pf_init = [ ( nodes[idx], nodes_u[idx], nodes_src_i[idx], nodes_dst_i[idx] ) for idx in 1:length(nodes_keys) ]    
    
    return  tuple_vec_comps_pf_init
      
end

function pf_self_x0(component,  uh, src_i, dst_i)

    x0_f_t_aux = pf_self_init_state(component, uh, src_i, dst_i)
       
    return x0_f_t_aux.x0 
    
end

function pf_self_f_t(component,   uh, src_i, dst_i )

    x0_f_t_aux = pf_self_init_state(component,  uh, src_i, dst_i)
       
    return x0_f_t_aux.f_t 
    
end


function pf_self_aux(component,  uh, src_i, dst_i)

    x0_f_t_aux = pf_self_init_state(component, uh, src_i, dst_i)
       
    return x0_f_t_aux.aux 
    
end

function pf_self_x0_f_t_aux(component,  uh, src_i, dst_i)
       
    return pf_self_init_state(component, uh, src_i, dst_i) 
    
end


function  pf_get_components_x0(nodes_collection, edges_collection, nodes_u )

    tuple_vec_comps_pf_init = pf_get_tuple_nodes_pf_init(nodes_collection, edges_collection, nodes_u )

    vec_comps_x0 = map((arg) -> pf_self_x0( arg[1], arg[2], arg[3], arg[4] ), tuple_vec_comps_pf_init)
    
    return  vcat(vec_comps_x0...)
      
end


function  pf_get_components_f_t( nodes_collection, edges_collection, nodes_u  )

    tuple_vec_comps_pf_init = pf_get_tuple_nodes_pf_init( nodes_collection, edges_collection, nodes_u  )

    vec_comps_f_t = map((arg) -> pf_self_f_t( arg[1], arg[2], arg[3], arg[4]  ) , tuple_vec_comps_pf_init)
    
    return  vec_comps_f_t
      
end


function  pf_get_components_aux( nodes_collection, edges_collection, nodes_u  )

    tuple_vec_comps_pf_init = pf_get_tuple_nodes_pf_init( nodes_collection, edges_collection, nodes_u  )

    vec_comps_aux = map((arg) -> pf_self_aux( arg[1], arg[2], arg[3], arg[4]  ), tuple_vec_comps_pf_init)
    
    return  vec_comps_aux
      
end


function  pf_get_components_x0_f_t_aux( nodes_collection, edges_collection, nodes_u  )

    tuple_vec_comps_pf_init = pf_get_tuple_nodes_pf_init( nodes_collection, edges_collection, nodes_u  )

    vec_comps_x0_f_t_aux = map((arg) -> pf_self_x0_f_t_aux( arg[1], arg[2], arg[3], arg[4] ), tuple_vec_comps_pf_init)
    
    return  vec_comps_x0_f_t_aux
      
end


function  pf_get_nodes_x0( nodes_collection, edges_collection, nodes_u  )
    
    return pf_get_components_x0( nodes_collection, edges_collection, nodes_u )
      
end

function  pf_get_nodes_f_t( nodes_collection, edges_collection, nodes_u )

    tuple_vec_nodes_pf_init = pf_get_tuple_nodes_pf_init( nodes_collection, edges_collection, nodes_u  )

    vec_nodes_f_t = map((arg) -> pf_self_f_t( arg[1], arg[2], arg[3], arg[4] ), tuple_vec_nodes_pf_init )
    
    return  vec_nodes_f_t 
      
end

function  pf_get_nodes_aux(nodes_collection, edges_collection, nodes_u )

    tuple_vec_nodes_pf_init = pf_get_tuple_nodes_pf_init(nodes_collection, edges_collection, nodes_u)

    vec_nodes_aux = map((arg) -> pf_self_aux( arg[1], arg[2], arg[3], arg[4] ), tuple_vec_nodes_pf_init )
    
    return  vec_nodes_aux 
      
end

#-------------------------------------------------------

function pf_self_init_state(branch::PiModelLine, edge_src_dst_node_u)

    Ybr = calc_branch_Ybr(branch.param_values[1:3]..., 1,1)
    
    edges_from_and_to_current = Ybr * edge_src_dst_node_u

    from_current = edges_from_and_to_current[1]
    
    to_current   = edges_from_and_to_current[2]
    
    if_r = real(from_current)
    if_i = imag(from_current)
    it_r = real(to_current)
    it_i = imag(to_current)

    
    return (x0 = [if_r, if_i, it_r, it_i], f_t =[0.0, 0.0])   
    
end



function pf_self_init_state(branch::Transformer, edge_src_dst_node_u )

    Ybr = calc_branch_Ybr(branch.param_values[1:4]..., 1)
    
    edges_from_and_to_current = Ybr * edge_src_dst_node_u

    from_current = edges_from_and_to_current[1]
    
    to_current   = edges_from_and_to_current[2]
    
    if_r = real(from_current)
    if_i = imag(from_current)
    it_r = real(to_current)
    it_i = imag(to_current)

    
    return (x0 = [if_r, if_i, it_r, it_i],
            f_t =[0.0, 0.0])   
    
end

function pf_self_x0(component, edge_src_dst_node_u)

    x0_f_t = pf_self_init_state(component, edge_src_dst_node_u)
       
    return x0_f_t.x0 
    
end


function pf_self_f_t(component,  edge_src_dst_node_u)

    x0_f_t = pf_self_init_state(component, edge_src_dst_node_u)
       
    return x0_f_t.f_t 
    
end

function get_edges_src_dst_node_u(edges_collection, nodes_u)
    edges_orientation = [an_edge.orientation for an_edge in collect(values(edges_collection)) ]

    # edges_src_dst_node_u = [[nodes_u[orient[1]], nodes_u[orient[2]]] for orient in  edges_orientation]

    return [[nodes_u[orient[1]], nodes_u[orient[2]]] for orient in  edges_orientation]
    
    
end


function  pf_get_tuple_edges_pf_init(edges_collection, nodes_u)

    # edges_orientation = [an_edge.orientation for an_edge in collect(values(comps_collection))]

    # edges_src_dst_node_u = [[nodes_u[orient[1]], nodes_u[orient[2]]] for orient in  edges_orientation]

    edges_src_dst_node_u = get_edges_src_dst_node_u(edges_collection, nodes_u)
    
    tuple_vec_comps_pf_init = [(an_edge, an_edge_src_dst_node_u) for (an_edge, an_edge_src_dst_node_u) in zip(collect(values(edges_collection)) , edges_src_dst_node_u ) ]    
    
    return  tuple_vec_comps_pf_init
      
end

function  pf_get_edges_x0(edges_collection,  nodes_u )

    tuple_vec_comps_pf_init = pf_get_tuple_edges_pf_init(edges_collection,  nodes_u )

    vec_comps_x0 = map((arg) -> pf_self_x0(arg[1], arg[2]), tuple_vec_comps_pf_init)
    
    return  vcat(vec_comps_x0...)
      
end


function  pf_get_edges_f_t(edges_collection,   nodes_u )

    tuple_vec_edges_pf_init = pf_get_tuple_edges_pf_init(edges_collection,  nodes_u )

    vec_edges_f_t = map((arg) -> pf_self_f_t(arg[1], arg[2]), tuple_vec_edges_pf_init)
    
    return  vec_edges_f_t
      
end


function pf_init_operationpoint(nd::NetworkData, nodes_u)

    nodes_state_init = pf_get_nodes_x0(nd.nodes, nd.edges, nodes_u)

    edges_state_init = pf_get_edges_x0(nd.edges, nodes_u)

    state_init = vcat(nodes_state_init, edges_state_init)
    
    return state_init 
end




########################################################
# Post Power flow Parameters init
########################################################


function self_post_pf_init_parameters(node::Union{ plant_no_gov_no_avr_wt_loc_load }, key, nd, power_flow_data)

   
    vh  = power_flow_data[1]
    θh  = power_flow_data[2]
    ph  = power_flow_data[3]
    qh  = power_flow_data[4]
    # i_r = power_flow_data[5]
    # i_i = power_flow_data[6]
    pg  = power_flow_data[7]
    qg  = power_flow_data[8]
    
    vf0_pm0  = self_aux(node,  power_flow_data)

    vf0, pm0 =  vf0_pm0

    ir_ii, _ =  f_t_self(node,  power_flow_data)
    
    i_r, i_i = ir_ii
    
    p_order0 = pm0
    
    v_ref = vf0

    ω_ref0 = 1.0
    
    vec_param_symb = Symbol[:P, :Q, :p_order, :v_ref, :ω_ref]
    
    vec_param_pf = Float64[ pg, qg, p_order0, v_ref, ω_ref0]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end




function self_post_pf_init_parameters(node::Union{pf_plant_SM_v6, plant_SM_idq, plant_SM_v6, plant_SM_system_matrix  }, key, nd, power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]
    # i_r = power_flow_data[5]
    # i_i = power_flow_data[6]
    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    
    vf0_pm0  = self_aux(node,  power_flow_data)

    vf0, pm0 =  vf0_pm0

    ir_ii, _ =  f_t_self(node,  power_flow_data)
    
    i_r, i_i = ir_ii
    
    p_order0 = pm0
    
    v_ref = vf0

    ω_ref0 = 1.0
    
    vec_param_symb = Symbol[:P, :Q, :p_order, :v_ref, :ω_ref]
    
    vec_param_pf = Float64[pg, qh, p_order0, v_ref, ω_ref0]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end


function self_post_pf_init_parameters(node::Union{ plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_system_matrix }, key, nd, power_flow_data)
    
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
    
    f_t = self_f_t(node,  power_flow_data)
    
    ω_ref0 = 1.0
    
    v_ref = f_t.avr_f_t[1]
    
    vec_param_symb = Symbol[:P, :Q, :v_ref, :ω_ref]
    
    vec_param_pf = Float64[0.0, qg, v_ref, ω_ref0]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end



function self_post_pf_init_parameters(node::Union{ plant_no_gov_wt_loc_load, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6 }, key, nd, power_flow_data)
    
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
    
    f_t = self_f_t(node,  power_flow_data)
    
    ω_ref0 = 1.0
    
    v_ref = f_t.avr_f_t[1]
    
    vec_param_symb = Symbol[:P, :Q, :v_ref, :ω_ref]
    
    # vec_param_pf = Float64[ph, qh, v_ref, ω_ref0]
    
    vec_param_pf = Float64[0.0, qg, v_ref, ω_ref0]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end



function self_post_pf_init_parameters(node::Union{ plant_wt_loc_load, plant_wt_loc_load_idq, plant_wt_loc_load_v6 }, key, nd, power_flow_data)
    
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
    
    f_t = self_f_t(node,  power_flow_data)
    
    p_order0, ω_ref0 = f_t.gov_f_t
    
    v_ref = f_t.avr_f_t[1]
    
    vec_param_symb = Symbol[:P, :Q, :p_order, :v_ref, :ω_ref]
    
    # vec_param_pf = Float64[ph, qh, p_order0, v_ref, ω_ref0]

    vec_param_pf = Float64[pg, qg, p_order0, v_ref, ω_ref0]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end



function self_post_pf_init_parameters(node::Union{ plant_millano, plant_rscad_idq, plant_rscad_v6, plant_rscad, plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb, plant_SM_2axis_cb_system_matrix }, key, nd, power_flow_data)

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
    
    f_t = self_f_t(node,  power_flow_data)
    
    p_order0, ω_ref0 = f_t.gov_f_t
    
    v_ref = f_t.avr_f_t[1]
    
    vec_param_symb = Symbol[:P, :Q, :p_order, :v_ref, :ω_ref]
    
    vec_param_pf = Float64[ pg, qg, p_order0, v_ref, ω_ref0]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end


function self_post_pf_init_parameters(node::Union{plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb  }, key, nd, power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    f_t = self_f_t(node,  power_flow_data)
    
    p_order0, ω_ref0 = f_t.gov_f_t
    
    v_ref = f_t.avr_f_t[1]
    
    vec_param_symb = Symbol[:P, :Q, :p_order, :v_ref, :ω_ref]
    
    vec_param_pf = Float64[ pg, qg, p_order0, v_ref, ω_ref0]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end

 


function self_post_pf_init_parameters(node::Union{ plant_cb_inf, plant_cb_inf_bus  }, key, nd, power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]
    
    f_t = self_f_t(node,  power_flow_data)

    _, _, vh, θh = f_t.gen_f_t
    
    p_order0, ω_ref0 = f_t.gov_f_t
    
    v_ref = f_t.avr_f_t[1]

    uh  = vh * exp(im * θh)
           
    vec_param_symb = Symbol[:P, :Q, :p_order, :v_ref, :ω_ref, :U]
    
    vec_param_pf = Union{Float64,ComplexF64}[ pg, qg, p_order0, v_ref, ω_ref0, uh]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end
 

function self_post_pf_init_parameters(node::Union{plant_Infinite_bus, plant_Infinite_cb_bus }, key, nd, power_flow_data)

    vh = power_flow_data[1]
    θh = power_flow_data[2]
    ph = power_flow_data[3]
    qh = power_flow_data[4]

    pg   = power_flow_data[7]
    qg   = power_flow_data[8]
    ig_r = power_flow_data[9]
    ig_i = power_flow_data[10]  
    
    # vh_θh_ph_qh = self_aux(node,  power_flow_data)
    # vh, θh, ph, qh = vh_θh_ph_qh
        
    uh  = vh * exp(im * θh)
    
    vec_param_symb  = Symbol[:U, :P, :Q]
    
    vec_param_pf = Union{Float64,ComplexF64}[ uh, pg, qg ]    
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing    
    
end



function self_post_pf_init_parameters(node::Union{ plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load  }, key, nd, power_flow_data)

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
    
    ph_qh_kPL_kQL = self_aux(node,  power_flow_data)
    
    _, _, kPL, kQL = ph_qh_kPL_kQL
    
    vec_param_symb = Symbol[:P, :Q, :kPL, :kQL]
    
    vec_param_pf = Union{Float64, Tuple{Float64, Float64, Float64}}[ ph, qh, kPL, kQL ]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end



function self_post_pf_init_parameters(node::Union{ plant_Transmission_t1, plant_Transmission_t2 }, key, nd, power_flow_data)

    # vh = power_flow_data[1]
    # θh = power_flow_data[2]
    # ph = power_flow_data[3]
    # qh = power_flow_data[4]
    # i_r = power_flow_data[5]
    # i_i = power_flow_data[6]
    # pg   = power_flow_data[7]
    # qg   = power_flow_data[8]
    # ig_r = power_flow_data[9]
    # ig_i = power_flow_data[10]
    
    ph_qh = self_aux(node,  power_flow_data)

    ph, qh = ph_qh   
    
    vec_param_symb  = Symbol[:P, :Q]
    
    vec_param_pf = Float64[ph, qh ]    
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    
    return nothing    
    
end



#--------------------------------------------------------


"""
# Not needed since SM_2axis machines do not have the parameters that are being set.

The parameters can be set for a plant with stand alone SM_2axis machines.

function self_post_pf_init_parameters(node::Union{ SM_2axis_cb_idq, SM_2axis_cb_v6, SM_2axis_cb_direct, SM_2axis_cb_millano, SM_2axis_cb }, key, nd, power_flow_data)

    vf0_pm0 = self_aux(node,  power_flow_data)
    
    vf0, pm0 = vf0_pm0
    
    vec_param_symb = Symbol[:p_order, :v_ref, :ω_ref]
    
    vec_param_pf = Float64[ pm0, vf0, 1.0]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)
        nd.nodes_param[bus_idx][bus_para] = para_value
    end
    return nothing
end


function self_post_pf_init_parameters(node::Union{ SM_2axis_cb_inf, SM_2axis_cb_inf_bus  }, key, nd, power_flow_data)

    vf0_pm0_vh_θh = self_aux(node,  power_flow_data)
    
    vf0, pm0, vh, θh = vf0_pm0_vh_θh

    uh  = vh * exp(im * θh)

    
    vec_param_symb = Symbol[:p_order, :v_ref, :ω_ref, :U]
    
    vec_param_pf = Union{Float64,ComplexF64}[ pm0, vf0, 1.0, uh]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)
        nd.nodes_param[bus_idx][bus_para] = para_value
    end
    return nothing
end

"""

function self_post_pf_init_parameters(node::Union{ Infinite_cb_bus, Infinite_bus }, key, nd, power_flow_data)

    vh_θh_ph_qh = self_aux(node,  power_flow_data)

    vh, θh, ph, qh = vh_θh_ph_qh   
    
    uh  = vh * exp(im * θh)
    
    vec_param_symb = Symbol[:U, :P, :Q]
    
    vec_param_pf = Union{Float64,ComplexF64}[ uh, ph, qh ]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing    
    
end



function self_post_pf_init_parameters(node::Union{ PQ_Const_P, PQ_Const_I, PQ_Const_Z  }, key, nd, power_flow_data)

    ph_qh_kPL_kQL = self_aux(node,  power_flow_data)
    
    ph, qh, kPL, kQL = ph_qh_kPL_kQL

    
    vec_param_symb = Symbol[:P, :Q, :kPL, :kQL]
    
    vec_param_pf = Union{Float64, Tuple{Float64, Float64, Float64}}[ ph, qh, kPL, kQL ]
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    return nothing
end


function self_post_pf_init_parameters(node::Union{ Trans_t1_Node, Trans_t1_Node }, key, nd, power_flow_data)
    
    ph_qh = self_aux(node,  power_flow_data)

    ph, qh = ph_qh   
    
    vec_param_symb  = Symbol[:P, :Q]
    
    vec_param_pf = Union{Float64}[ph, qh ]    
    
    for (a_para, para_value) in zip(vec_param_symb, vec_param_pf)
        bus_idx, bus_para = node_params_idxs(nd, key, a_para)

        # nd.nodes_param[bus_idx][bus_para] = para_value
        set_component_params_value_in_idxs!(nd.nodes_param, bus_idx, bus_para, para_value)
    end
    
    return nothing    
    
end



#--------------------------------------------------------
#--------------------------------------------------------

function external_pf_init_post_pf_parameters!(nd::NetworkData, bus_dict_init, branch_dict_init) 
        
    bus_keys    = collect(keys(nd.nodes))
    
    branch_keys = collect(keys(nd.edges))

    vec_nodetype_init_post_pf_param_tuple = [(nd.nodes[a_key], a_key, nd, bus_dict_init[a_key])  for a_key in  bus_keys]

    for a_node in vec_nodetype_init_post_pf_param_tuple
        self_post_pf_init_parameters(a_node...)
    end

    # for a_key in bus_keys
    #     self_post_pf_init_parameters(working_nd.nodes[a_key], a_key, working_nd, branch_dict_init[a_key]) 
    # end
    
    
    return nothing
end


function external_pf_init_post_pf_parameters!(nd::NetworkData, dict_pf_init; key_nodes="bus_dict_init", key_edges="branch_dict_init")
        
    bus_keys    = collect(keys(nd.nodes))
    
    branch_keys = collect(keys(nd.edges))

    vec_nodetype_init_post_pf_param_tuple = [(nd.nodes[a_key], a_key, nd, dict_pf_init[key_nodes][a_key])  for a_key in  bus_keys]

    for a_node in vec_nodetype_init_post_pf_param_tuple
        self_post_pf_init_parameters(a_node...)
    end


    # for a_key in bus_keys
    #     self_post_pf_init_parameters(working_nd.nodes[a_key], a_key, working_nd, dict_pf_init[key_nodes][a_key]) 
    # end
    
    
    return nothing
end


function init_post_pf_parameters!(nd::NetworkData, dict_pf_init; key_nodes="bus_dict_init", key_edges="branch_dict_init")
    
    working_nd  = deepcopy(nd)
    
    bus_keys    = collect(keys(working_nd.nodes))
    
    branch_keys = collect(keys(working_nd.edges))

    vec_nodetype_init_post_pf_param_tuple = [(working_nd.nodes[a_key], a_key, working_nd, dict_pf_init[key_nodes][a_key])  for a_key in  bus_keys]

    for a_node in vec_nodetype_init_post_pf_param_tuple
        self_post_pf_init_parameters(a_node...)
    end


    # for a_key in bus_keys
    #     self_post_pf_init_parameters(working_nd.nodes[a_key], a_key, working_nd, dict_pf_init[key_nodes][a_key]) 
    # end
    
    
    return working_nd
end

########################################################

