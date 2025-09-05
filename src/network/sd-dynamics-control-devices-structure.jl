
########################################################
# ------------------------------------------------------
# models for building composite plants
# ------------------------------------------------------
########################################################

#########################################################
# ------------------------------------------------------
#  controllers
# ------------------------------------------------------
#########################################################


"""
```math
Lag

X(s)       K
--   = --------
U(s)   1 + s Ta

dx     1         K 
-- = - - x(t) + -- u(t)
dt     Ta        Ta
            
y(t) = x(t)


lead_lag_parallel :

X(s)   1 + s Tc 
--   = --------
U(s)   1 + s Tb

dx     1         (1 - Tc/Tb)
-- = - - x(t) +  ---------- u(t)
dt     Tb            Tb
              Tc
y(t) = x(t) + -- u(t)
              Tb
```

"""

function gov_ieee_tgov1_cb!(dx, x, p_agg, t)
        
    cb_sw, u_gen_ω, f_t, p = p_agg

    sw = cb_sw[1]
    
    T1       = p[1]
    T2       = p[2]
    T3       = p[3]
    Dt       = p[4]
    p_max    = p[5]
    p_min    = p[6]
    R        = p[7]
     
    # p_ref    = f_t[1]

    p_order0 = f_t[1]

    ω_ref0   = f_t[2]
    
    p_ref    = p_order0 * R
    
    xg1, xg2, τm_tilade, ω_ref = x

    # dx[1] = (1/T1) * ( (1/R) * (P_ref - (u_gen_ω/ω_ref - 1.0 )) - x1 )
    # dx[1] = (1/T1) * ( (p_order0 - (1/R) *(u_gen_ω/ω_ref - 1.0 )) - x1 )
    
    if sw == 0

        dx[1] = (1/T1) * (( p_order0 - (1/R) *
            ( u_gen_ω/ω_ref - 1)) - xg1 )
        
    elseif sw == 1
        
        dx[1] = 0.0
        
    elseif sw == 2
        
        dx[1] = 0.0
    end    

    dx[2] = (1/T3) * (( 1 - T2/T3 ) * xg1 - xg2 )
    
    dx[3] = xg2 + (T2/T3) * xg1 - Dt *
        (u_gen_ω/ω_ref - 1) - τm_tilade

    dx[4] = ω_ref0 - ω_ref

    return nothing
end



@kwdef struct gov_ieee_tgov1_cb <: SdGov
    
    T1::Float64    # Steam bowl time constant
    T2::Float64    # 
    T3::Float64    # 
    Dt::Float64    # Turbine damping coefficient
    p_max::Float64 # Maximum turbine output
    p_min::Float64 # Minimum turbine output
    R::Float64     # Droop

    comp_type::Symbol = :gov_ieee_tgov1_cb
    
    state_vars_syms::Vector{Symbol} =
        Symbol[ :xg1, :xg2]
    
    algebraic_vars_syms::Vector{Symbol} =
        Symbol[ :τm_tilade, :ω_ref]
    
    im_algebraic_vars_syms::Vector{Symbol} =
        Symbol[ :τm_tilade ]

    im_vars_syms::Vector{Symbol} =
        Symbol[ state_vars_syms...;im_algebraic_vars_syms...]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )         
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]

    dict_im_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( state_vars_syms  ) )

    dict_im_algebraic_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_algebraic_vars_syms  ) ) 
    
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector(length(state_vars_syms), length(algebraic_vars_syms) )
    param::Vector{Symbol} = Symbol[:T1, :T2, :T3, :Dt,:p_max,:p_min, :R]
    param_values::Vector{Float64} = Float64[T1, T2, T3, Dt, p_max,p_min, R ]
    func::Vector{Function}           = Function[gov_ieee_tgov1_cb!]
    control_sig_syms::Vector{Symbol} = Symbol[ :p_order0, :ω_ref0]
    control_sig::Vector{Float64}     = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol}  = Symbol[:τm_tilade]
    output_sig::Vector{Float64}      = ones(length(output_sig_syms))

    output_sig_state_coeff::Vector{Float64} = [ T2/T3, 1.0 ]

    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    cb_state_event_func::Vector{Function}      = Function[ ]
    cb_state_affect_func::Vector{Function}     = Function[ ]    
    cb_state_syms::Vector{Symbol}              = Symbol[ ]  
    cb_state_conditions::Vector{Float64}       = Float64[ ]     
    cb_state_values::Vector{Function}          = Function[ ]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)
    cb_dyn_state_event_func::Vector{Function}  = Function[ fun_condition_lower_lim, fun_condition_upper_lim ]
    cb_dyn_state_affect_func::Vector{Function} = Function[ fun_affect_anti_windup_lower_lim!, fun_affect_anti_windup_upper_lim! ]    
    cb_dyn_state_syms::Vector{Symbol}   = Symbol[:xg1, :xg1]
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
    cb_dyn_state_conditions::Vector{Float64}   = Float64[ p_min, p_max ]  
    cb_dyn_state_values::Vector{Function}      = Function[ (u,t) -> p_min, (u,t)-> p_max ]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ 0, 0 ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)


    Ax_S_Sgen::Function = Ax_S_Sgen_gov_ieee_tgov1_cb
    
    Ax_A_Sgen::Function = Ax_A_Sgen_gov_ieee_tgov1_cb
    

    Ax_S_S::Function    = Ax_S_S_gov_ieee_tgov1_cb

    Ax_S_A::Function    = Ax_S_A_gov_ieee_tgov1_cb

    Ax_A_S::Function    = Ax_A_S_gov_ieee_tgov1_cb

    Ax_A_A::Function    = Ax_A_A_gov_ieee_tgov1_cb

    
    Bx_S::Function      = Bx_S_gov_ieee_tgov1_cb

    Cx_S::Function      = Cx_S_gov_ieee_tgov1_cb

    Bx_A::Function      = Bx_A_gov_ieee_tgov1_cb

    Cx_A::Function      = Cx_A_gov_ieee_tgov1_cb
    
    stab_Ac::Function   = stab_Ac_gen_ieee_tgov1_cb
    
end


#--------------------------------------------------------


function gov_t0_cb!(dx, x, p_agg, t)
        
    cb_sw, u_gen_ω, f_t, p = p_agg
    
    T3       = p[1]
    T4       = p[2]
    T5       = p[3]
    Tc       = p[4]
    Ts       = p[5]
    p_max    = p[6]
    p_min    = p[7]
    R        = p[8]
     
    p_order0 = f_t[1]    
    ω_ref0   = f_t[2]    
 
    xg1, xg2, xg3, τm_tilade, ω_ref, phat_in, p_in = x
            
    # dx[1] = (1/Ts)*(p_order0 - (1/R)*(u_gen_ω - ω_ref) - xg1)

    dx[1] = (1/Ts)*(phat_in - xg1)
    
    dx[2] = (1/Tc) * ((1 - T3/Tc) * xg1 - xg2 )

    dx[3] = (1/T5)*((1 - T4/T5) * (xg2 + (T3/Tc) * xg1) - xg3)

    dx[4] = xg3 + (T4/T5) * (xg2 + (T3/Tc) * xg1) - τm_tilade

    dx[5] = ω_ref0 - ω_ref

    dx[6] = p_order0 - (1/R) * (u_gen_ω/ω_ref0 - 1 ) - phat_in
    
    dx[7] = phat_in - p_in
            
    return nothing
end


@kwdef struct gov_t0_cb <: SdGov
    
    T3::Float64    # Transient gain time constant
    T4::Float64    # Power fraction time constant
    T5::Float64    # Reheat time constant
    Tc::Float64    # Servo time constant
    Ts::Float64    # Governor time constant
    p_max::Float64 # Maximum turbine output
    p_min::Float64 # Minimum turbine output
    R::Float64     # Droop

    comp_type::Symbol = :gov_t0_cb
    
    state_vars_syms::Vector{Symbol} = Symbol[ :xg1, :xg2, :xg3 ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :τm_tilade, :ω_ref, :phat_in, :p_in ]

    im_algebraic_vars_syms::Vector{Symbol} = Symbol[ :τm_tilade, :phat_in, :p_in ]

    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )         
    
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]

    dict_im_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( state_vars_syms  ) )

    dict_im_algebraic_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_algebraic_vars_syms  ) ) 
        
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector(length(state_vars_syms), length(algebraic_vars_syms) )
    param::Vector{Symbol} = Symbol[:T3, :T4,:T5, :Tc,:Ts,:p_max,:p_min, :R]
    param_values::Vector{Float64} = Float64[T3, T4, T5, Tc, Ts, p_max,p_min, R ]
    func::Vector{Function}           = Function[gov_t0_cb!]
    control_sig_syms::Vector{Symbol} = Symbol[ :p_order0, :ω_ref0]
    control_sig::Vector{Float64}     = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol}  = Symbol[:τm_tilade]
    output_sig::Vector{Float64}      = ones(length(output_sig_syms))

    output_sig_state_coeff::Vector{Float64} = [(T4/T5)*(T3/Tc), (T4/T5), 1.0 ]
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    cb_state_event_func::Vector{Function}      = Function[fun_condition_lower_lim, fun_condition_upper_lim]
    cb_state_affect_func::Vector{Function}     = Function[fun_affect_windup_lower_lim!, fun_affect_windup_upper_lim!]    
    cb_state_syms::Vector{ Union{Symbol,Tuple{Symbol,Symbol}} }  = Union{Symbol,Tuple{Symbol,Symbol}}[(:phat_in, :p_in), (:phat_in, :p_in) ]  
    cb_state_conditions::Vector{Float64}       = Float64[p_min, p_max]     
    cb_state_values::Vector{Function}          = Function[(u,t) -> p_min, (u,t)-> p_max]
    cb_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_state_syms[state_sym] : ( dict_state_syms[state_sym[1]], dict_state_syms[state_sym[1]])  for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)
    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]  
    cb_dyn_state_values::Vector{Function}      = Function[]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    # Ax::Function = Ax_gov_t0_cb
    # Ax_S_Sgen::Function = Ax_S_Sgen_gov_t0_cb
    # Bx::Function = Bx_gov_t0_cb
    # Cx::Function = Cx_gov_t0_cb
    # Ac_S_with_A::Function = Ac_Sgov_with_Agov_gov_t0_cb
    # Ac_A_with_S::Function = Ac_Agov_with_Sgov_gov_t0_cb
    # Ax_A::Function = Ax_Agov_gov_t0_cb
    # Ac_A_with_Sgen::Function = Ac_Agov_with_Sgen_gov_t0_cb
    # Bx_A::Function = Bx_Agov_gov_t0_cb
    # Cx_A::Function = Cx_Agov_gov_t0_cb

    Ax_S_Sgen::Function = Ax_S_Sgen_gov_t0_cb
    Ax_A_Sgen::Function = Ax_A_Sgen_gov_t0_cb
    Ax_S_S::Function       = Ax_S_S_gov_t0_cb
    Ax_S_A::Function       = Ax_S_A_gov_t0_cb
    Ax_A_S::Function       = Ax_A_S_gov_t0_cb
    Ax_A_A::Function       = Ax_A_A_gov_t0_cb
    Bx_S::Function           = Bx_S_gov_t0_cb
    Cx_S::Function           = Cx_S_gov_t0_cb
    Bx_A::Function           = Bx_A_gov_t0_cb
    Cx_A::Function           = Cx_A_gov_t0_cb
    
    stab_Ac::Function = stab_Ac_gen_t0_cb
    
end



function gov_t1_cb!(dx, x, p_agg, t)
        
    cb_sw, u_gen_ω, f_t, p = p_agg
    
    T3       = p[1]
    T4       = p[2]
    T5       = p[3]
    Tc       = p[4]
    Ts       = p[5]
    p_max    = p[6]
    p_min    = p[7]
    R        = p[8]
     
    p_order0 = f_t[1]    
    ω_ref0   = f_t[2]    
 
    xg1, xg2, xg3, τm_tilade, ω_ref, phat_in = x
            
    # dx[1] = (1/Ts) * (p_order0 - (1/R) * (u_gen_ω/ω_ref - 1) - xg1)

    dx[1] = (1/Ts) * (phat_in - xg1)

    dx[2] = (1/Tc) * ((1 - T3/Tc) * xg1 - xg2 )

    dx[3] = (1/T5) * ((1 - T4/T5) * (xg2 + (T3/Tc) * xg1) - xg3 )

    dx[4] = xg3 + (T4/T5) * (xg2 + (T3/Tc) * xg1) - τm_tilade

    dx[5] = ω_ref0 - ω_ref
    
    dx[6] = p_order0 - (1/R) * (u_gen_ω/ω_ref0 - 1 ) - phat_in
            
    return nothing
end



@kwdef struct gov_t1_cb <: SdGov
    
    T3::Float64    # Transient gain time constant
    T4::Float64    # Power fraction time constant
    T5::Float64    # Reheat time constant
    Tc::Float64    # Servo time constant
    Ts::Float64    # Governor time constant
    p_max::Float64 # Maximum turbine output
    p_min::Float64 # Minimum turbine output
    R::Float64     # Droop

    comp_type::Symbol = :gov_t1_cb

    state_vars_syms::Vector{Symbol} = Symbol[ :xg1, :xg2, :xg3 ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :τm_tilade, :ω_ref, :phat_in ]

    im_algebraic_vars_syms::Vector{Symbol} = Symbol[ :τm_tilade, :phat_in ]

    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )            
    
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]


    dict_im_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( state_vars_syms  ) )

    dict_im_algebraic_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_algebraic_vars_syms  ) ) 
    
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector(length(state_vars_syms), length(algebraic_vars_syms) )
    param::Vector{Symbol} = Symbol[:T3, :T4,:T5, :Tc,:Ts,:p_max,:p_min, :R]
    param_values::Vector{Float64} = Float64[T3, T4, T5, Tc, Ts, p_max,p_min, R ]
    func::Vector{Function}           = Function[gov_t1_cb!]
    control_sig_syms::Vector{Symbol} = Symbol[ :p_order0, :ω_ref0]
    control_sig::Vector{Float64}     = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol}  = Symbol[:τm_tilade]
    output_sig::Vector{Float64}      = ones(length(output_sig_syms))
    
    output_sig_state_coeff::Vector{Float64} = [(T4/T5)*(T3/Tc), (T4/T5), 1.0 ]
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    cb_state_event_func::Vector{Function}      = Function[fun_condition_lower_lim, fun_condition_upper_lim]
    cb_state_affect_func::Vector{Function}     = Function[fun_affect_windup_lower_lim!, fun_affect_windup_upper_lim!]    
    cb_state_syms::Vector{Symbol}              = Symbol[:phat_in, :phat_in]  
    cb_state_conditions::Vector{Float64}       = Float64[p_min, p_max]     
    cb_state_values::Vector{Function}          = Function[(u,t) -> p_min, (u,t)-> p_max]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)
    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]  
    cb_dyn_state_values::Vector{Function}      = Function[]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ ]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    # Ax::Function = Ax_gov_t1_cb

    # Ax_S_Sgen::Function = Ax_S_Sgen_gov_t1_cb

    # Bx::Function = Bx_gov_t1_cb

    # Cx::Function = Cx_gov_t1_cb

    # Ac_S_with_A::Function = Ac_Sgov_with_Agov_gov_t1_cb

    # Ac_A_with_S::Function = Ac_Agov_with_Sgov_gov_t1_cb

    # Ax_A::Function = Ax_Agov_gov_t1_cb

    # Ac_A_with_Sgen::Function = Ac_Agov_with_Sgen_gov_t1_cb

    # Bx_A::Function = Bx_Agov_gov_t1_cb

    # Cx_A::Function = Cx_Agov_gov_t1_cb


    Ax_S_Sgen::Function = Ax_S_Sgen_gov_t1_cb
    Ax_A_Sgen::Function = Ax_A_Sgen_gov_t1_cb
    Ax_S_S::Function       = Ax_S_S_gov_t1_cb
    Ax_S_A::Function       = Ax_S_A_gov_t1_cb
    Ax_A_S::Function       = Ax_A_S_gov_t1_cb
    Ax_A_A::Function       = Ax_A_A_gov_t1_cb
    Bx_S::Function           = Bx_S_gov_t1_cb
    Cx_S::Function           = Cx_S_gov_t1_cb
    Bx_A::Function           = Bx_A_gov_t1_cb
    Cx_A::Function           = Cx_A_gov_t1_cb

    stab_Ac::Function = stab_Ac_gen_t1_cb
    
end


function gov_t1_cb_sauer!(dx, x, p_agg, t)
        
    cb_sw, u_gen_ω, f_t, p = p_agg

    # Sauer section 5.2, pg 91.
    # Also see page 82
    
    # TSV = 0.2, PC = 0.7,  RD = 0.05
    # TCH = 0.4, ωs = 2π*60
    
    # Tch      = p[1]
    # Tsv      = p[2]
    # Ts       = p[3]
    # Psv_max  = p[4]
    # Psv_min  = p[5]    
    # RD       = p[6]

    """
    dτm/dt = (1/Tch) * (Psv - τm )
    dPsv   = (1/Tsv) * (Pc - (1/RD) * (ω/ωs - 1.0 ) - Psv)
    0 ≤ Psv ≤ Psv_max
    
    or
    
    dτm/dt = (1/Tch) * (Psv - τm )
    dPsv  = (1/Tsv) * (Pc - ϵ * ωt/(RD * Ts))
    Ts = √(2*H / ωs)
    ωt = Ts * (ω - ωs)
    ϵ = 1/ωs
    0 ≤ Psv ≤ Psv_max
    H should be large enough that ϵ << Ts 
    """
    
    Tc       = p[1]
    Ts       = p[2]
    p_max    = p[3]
    p_min    = p[4]
    R        = p[5]
     
    p_order0 = f_t[1]
    
    ω_ref0   = f_t[2]
    
    #Psv, TM, ωs      TM
    xg1, xg2, ω_ref, τm_tilade = x

    dx[1] = (1/Ts) * (-xg1 + p_order0 - (1 / R) * (u_gen_ω/ω_ref0 - 1)  )

    dx[2] = (1/Tc) * (-xg2 + xg1  )

    dx[3] = ω_ref0 - ω_ref

    dx[4] = xg2  - τm_tilade


    
    return nothing
end


@kwdef struct gov_t1_cb_sauer <: SdGov
    
    Tc::Float64    # Servo time constant
    Ts::Float64    # Governor time constant
    p_max::Float64 # Maximum turbine output
    p_min::Float64 # Minimum turbine output
    R::Float64     # Droop

    comp_type::Symbol = :gov_t1_cb_sauer
    
    state_vars_syms::Vector{Symbol} = Symbol[ :xg1, :xg2]

    algebraic_vars_syms::Vector{Symbol} = Symbol[:ω_ref, :τm_tilade]
    
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[ :τm_tilade ]
    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )       
    
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]

    dict_im_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( state_vars_syms  ) )

    dict_im_algebraic_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_algebraic_vars_syms ) ) 
    
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    param::Vector{Symbol} = Symbol[:Tc,:Ts,:p_max,:p_min,:R]
    param_values::Vector{Float64} = Float64[Tc,Ts,p_max,p_min, R ]
    func::Vector{Function} = Function[gov_t1_cb_sauer!]
    control_sig_syms::Vector{Symbol} = Symbol[ :p_order0, :ω_ref0]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    # output_sig_syms::Vector{Symbol} = Symbol[:xg2]
    
    output_sig_syms::Vector{Symbol} = Symbol[ :τm_tilade ]
    
    output_sig::Vector{Float64} = ones(length(output_sig_syms))

    output_sig_state_coeff::Vector{Float64} = [0.0, 1.0 ]
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    cb_state_event_func::Vector{Function}      = Function[fun_condition_lower_lim, fun_condition_upper_lim]
    cb_state_affect_func::Vector{Function}     = Function[fun_affect_windup_lower_lim!, fun_affect_windup_upper_lim!]    
    cb_state_syms::Vector{Symbol}              = Symbol[ :xg1, :xg1 ]           
    cb_state_conditions::Vector{Float64}       = Float64[p_min, p_max]     
    cb_state_values::Vector{Function}          = Function[(u,t) -> p_min, (u,t)-> p_max]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)
    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]  
    cb_dyn_state_values::Vector{Function}      = Function[]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[0]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)
    
    # Ax::Function = Ax_gov_sauer

    # Ax_S_Sgen::Function = Ax_S_Sgen_gov_t1_cb_sauer

    # Bx::Function = Bx_gov_sauer

    # Cx::Function = Cx_gov_sauer

    # Ac_S_with_A::Function = Ac_Sgov_with_Agov_gov_t1_cb_sauer

    # Ac_A_with_S::Function = Ac_Agov_with_Sgov_gov_t1_cb_sauer

    # Ax_A::Function = Ax_Agov_gov_t1_cb_sauer

    # Ac_A_with_Sgen::Function = Ac_Agov_with_Sgen_gov_t1_cb_sauer

    # Bx_A::Function = Bx_Agov_gov_t1_cb_sauer

    # Cx_A::Function = Cx_Agov_gov_t1_cb_sauer

    Ax_S_Sgen::Function = Ax_S_Sgen_gov_t1_cb_sauer
    Ax_A_Sgen::Function = Ax_A_Sgen_gov_t1_cb_sauer    
    Ax_S_S::Function       = Ax_S_S_gov_t1_cb_sauer    
    Ax_S_A::Function       = Ax_S_A_gov_t1_cb_sauer
    Ax_A_S::Function       = Ax_A_S_gov_t1_cb_sauer
    Ax_A_A::Function       = Ax_A_A_gov_t1_cb_sauer    
    Bx_S::Function           = Bx_S_gov_t1_cb_sauer
    Cx_A::Function           = Cx_A_gov_t1_cb_sauer
    Bx_A::Function           = Bx_A_gov_t1_cb_sauer
    Cx_S::Function           = Cx_S_gov_t1_cb_sauer    

    stab_Ac::Function = stab_Ac_gov_gen_sauer
    
end

#--------------------------------------------------------
#--------------------------------------------------------


function avr_t0_cb!(dx, x, p_agg, t)
    
    cb_sw, u_gen_ur_ui, f_t, p = p_agg

    sw = cb_sw[1]
    
    u_r        = u_gen_ur_ui[1]
    
    u_i        = u_gen_ur_ui[2]

    vh         = abs(u_r + im * u_i)
    
    v_ref0     = f_t[1]
    
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
       
    # vm = x[1], vr1 = x[2], vr2 = x[3], vf_tilade = x[4], v_ref = x[5], 
     
    vm, vr1, vr2, vf_tilade, v_ref, vr1_hat, vf = x

    dx[1] = (1/Tr)  * (vh - vm )

    dx[2] = (Ka * ((v_ref0 - vm) + vr2 - (Kf/Tf) * vf_tilade ) - vr1) /Ta
    
    dx[3] = -vr2/Tf + (Kf/Tf^2) * vf_tilade

    dx[4] = (-1/Te) * (vf_tilade * (Ke + Sevf(
        Ae, Be, vf_tilade)) - vr1_hat)
    
    dx[5] = v_ref0 - v_ref

    if sw == 0
        
        dx[6] = vr1 - vr1_hat
        
    elseif sw == 1
        
        dx[6] = 0.0
        
    elseif sw == 2
        
        dx[6] = 0.0
    end

    dx[7] = vf_tilade - vf
    
    nothing
end




@kwdef struct avr_t0_cb <: SdAvr
    
    Ta::Float64 # Amplifier time constant
    Tf::Float64 # Stabilizer time constant
    Te::Float64 # Field circuit time constant
    Tr::Float64 # Measurement time constant
    Ka::Float64 # Amplifier gain
    Ke::Float64 # Field circuit integral deviation
    Kf::Float64 # Stabilizer gain
    V_R_max::Float64 # Maximum regulator voltage
    V_R_min::Float64 # Minimum regulator voltage
    # S_E_max::Float64 # 
    # S_E0_75max::Float64
    Ae::Float64 # 1st ceiling coeﬃcient
    Be::Float64 # 2nd ceiling coeﬃcient

    comp_type::Symbol = :avr_t0_cb

    state_vars_syms::Vector{Symbol} = Symbol[:vm, :vr1, :vr2, :vf_tilade]
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[:v_ref, :vr1_hat, :vf ]
    
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[ :vr1_hat, :vf ]
    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )       

    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]


    dict_im_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( state_vars_syms  ) )

    dict_im_algebraic_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_algebraic_vars_syms  ) ) 
    
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    param::Vector{Symbol}             = Symbol[ :Ta,:Te,:Tf,:Tr,:Ka,:Ke,:Kf,:V_R_max,:V_R_min, :Ae,:Be ]
    param_values::Vector{Float64} = Float64[ Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be ]
    func::Vector{Function}            = Function[avr_t0_cb!]
    control_sig_syms::Vector{Symbol}  = Symbol[:v_ref]
    control_sig::Vector{Float64}      = ones(length(control_sig_syms))
    # output_sig_syms::Vector{Symbol}   = Symbol[:vf_tilade]
    output_sig_syms::Vector{Symbol}   = Symbol[:vf]
    output_sig::Vector{Float64}       = ones(length(output_sig_syms))

    output_sig_state_coeff::Vector{Float64} = [0.0, 0.0, 0.0, 1.0 ]
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]     
    cb_state_values::Vector{Function}          = Function[]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)
    cb_dyn_state_event_func::Vector{Function}  = Function[ fun_condition_lower_lim, fun_condition_upper_lim ]
    cb_dyn_state_affect_func::Vector{Function} = Function[ fun_affect_anti_windup_lower_lim!, fun_affect_anti_windup_upper_lim! ]
    cb_dyn_state_syms::Vector{ Union{Symbol, Tuple{Symbol,Symbol}} }  = Union{Symbol, Tuple{Symbol,Symbol}}[ (:vr1, :vr1_hat), (:vr1, :vr1_hat) ]
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{ Union{Symbol,Tuple{Symbol, Symbol}}, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
    cb_dyn_state_conditions::Vector{Float64}   = Float64[ V_R_min, V_R_max ]  
    cb_dyn_state_values::Vector{Function}      = Function[(u,t) -> V_R_min, (u, t) -> V_R_max ]
    cb_dyn_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_state_syms[state_sym] : ( dict_state_syms[state_sym[1]], dict_state_syms[state_sym[1]])  for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}  = Int64[0, 0]
    cb_dyn_state_dim::Int64               = length(cb_dyn_state_conditions)

    # Ax::Function = Ax_avr_t0_cb

    # Ax_S_Sgen::Function = Ax_S_Sgen_avr_t0_cb

    # Bx::Function = Bx_avr_t0_cb

    # Cx::Function = Cx_avr_t0_cb

    # Ac_S_with_A::Function = Ac_Savr_with_Aavr_avr_t0_cb

    # Ac_A_with_S::Function = Ac_Aavr_with_Savr_avr_t0_cb

    # Ax_A::Function = Ax_Aavr_avr_t0_cb

    # Ac_A_with_Sgen::Function = Ac_Aavr_with_Sgen_avr_t0_cb

    # Bx_A::Function = Bx_Aavr_avr_t0_cb

    # Cx_A::Function = Cx_Aavr_avr_t0_cb
    

    Ax_S_Sgen::Function = Ax_S_Sgen_avr_t0_cb
    Ax_A_Sgen::Function = Ax_A_Sgen_avr_t0_cb
    Ax_S_S::Function       = Ax_S_S_avr_t0_cb
    Ax_S_A::Function       = Ax_S_A_avr_t0_cb
    Ax_A_S::Function       = Ax_A_S_avr_t0_cb
    Ax_A_A::Function       = Ax_A_A_avr_t0_cb
    Bx_S::Function           = Bx_S_avr_t0_cb
    Cx_S::Function           = Cx_S_avr_t0_cb
    Bx_A::Function           = Bx_A_avr_t0_cb
    Cx_A::Function           = Cx_A_avr_t0_cb

    stab_Ac::Function = stab_Ac_avr_gen_t0_cb
    
end



function avr_t1_cb!(dx, x, p_agg, t)
    
    cb_sw, u_gen_ur_ui, f_t, p = p_agg

    sw = cb_sw[1]
    
    u_r        = u_gen_ur_ui[1]
    
    u_i        = u_gen_ur_ui[2]

    vh         = abs(u_r + im * u_i)
    
    v_ref0     = f_t[1]
    
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
       
     
    vm, vr1, vr2, vf_tilade, v_ref, vf  = x

    dx[1] = (1/Tr)  * (vh - vm )
    
    if sw == 0
        
        dx[2] = (Ka * ((v_ref0 - vm) + vr2 - (Kf/Tf) * vf_tilade ) - vr1) /Ta
        
    elseif sw == 1
        
        dx[2] = 0.0
        
    elseif sw == 2
        
        dx[2] = 0.0
    end

    dx[3] = -vr2/Tf + (Kf/Tf^2) * vf_tilade

    dx[4] = (1/Te) * (vr1 - vf_tilade * (Ke + Sevf(Ae, Be, vf_tilade)) )    

    dx[5] = v_ref0 - v_ref

    dx[6] = vf_tilade - vf
    
    nothing
end




@kwdef struct avr_t1_cb <: SdAvr
    
    Ta::Float64 # Amplifier time constant
    Tf::Float64 # Stabilizer time constant
    Te::Float64 # Field circuit time constant
    Tr::Float64 # Measurement time constant
    Ka::Float64 # Amplifier gain
    Ke::Float64 # Field circuit integral deviation
    Kf::Float64 # Stabilizer gain
    V_R_max::Float64 # Maximum regulator voltage
    V_R_min::Float64 # Minimum regulator voltage
    # S_E_max::Float64 # 
    # S_E0_75max::Float64
    Ae::Float64 # 1st ceiling coeﬃcient
    Be::Float64 # 2nd ceiling coeﬃcient

    comp_type::Symbol = :avr_t1_cb
    
    state_vars_syms::Vector{Symbol} = Symbol[:vm, :vr1, :vr2, :vf_tilade]
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[:v_ref, :vf]
        
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[ :vf ]
    
    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )     
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]

    dict_im_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( state_vars_syms  ) )

    dict_im_algebraic_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_algebraic_vars_syms  ) ) 
    
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    param::Vector{Symbol}             = Symbol[ :Ta,:Te,:Tf,:Tr,:Ka,:Ke,:Kf,:V_R_max,:V_R_min, :Ae,:Be ]
    param_values::Vector{Float64} = Float64[ Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be ]
    func::Vector{Function}            = Function[avr_t1_cb!]
    control_sig_syms::Vector{Symbol}  = Symbol[:v_ref]
    control_sig::Vector{Float64}      = ones(length(control_sig_syms))

    # output_sig_syms::Vector{Symbol}   = Symbol[:vf_tilade]
    
    output_sig_syms::Vector{Symbol}   = Symbol[ :vf ]
    
    output_sig::Vector{Float64}       = ones(length(output_sig_syms))

    output_sig_state_coeff::Vector{Float64} = [0.0, 0.0, 0.0, 1.0 ]
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]     
    cb_state_values::Vector{Function}          = Function[]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)
    cb_dyn_state_event_func::Vector{Function}  = Function[ fun_condition_lower_lim, fun_condition_upper_lim ]
    cb_dyn_state_affect_func::Vector{Function} = Function[ fun_affect_anti_windup_lower_lim!, fun_affect_anti_windup_upper_lim! ]
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[ :vr1, :vr1 ]
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
    cb_dyn_state_conditions::Vector{Float64}   = Float64[ V_R_min, V_R_max ]  
    cb_dyn_state_values::Vector{Function}      = Function[(u,t) -> V_R_min, (u, t) -> V_R_max ]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[0, 0]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    # Ax::Function = Ax_avr_t1_cb

    # Ax_S_Sgen::Function = Ax_S_Sgen_avr_t1_cb

    # Bx::Function = Bx_avr_t1_cb

    # Cx::Function = Cx_avr_t1_cb

    # Ac_S_with_A::Function = Ac_Savr_with_Aavr_avr_t1_cb

    # Ac_A_with_S::Function = Ac_Aavr_with_Savr_avr_t1_cb

    # Ax_A::Function = Ax_Aavr_avr_t1_cb

    # Ac_A_with_Sgen::Function = Ac_Aavr_with_Sgen_avr_t1_cb

    # Bx_A::Function = Bx_Aavr_avr_t1_cb

    # Cx_A::Function = Cx_Aavr_avr_t1_cb

    Ax_S_Sgen::Function = Ax_S_Sgen_avr_t1_cb
    Ax_A_Sgen::Function = Ax_A_Sgen_avr_t1_cb
    Ax_S_S::Function       = Ax_S_S_avr_t1_cb
    Ax_S_A::Function       = Ax_S_A_avr_t1_cb
    Ax_A_S::Function       = Ax_A_S_avr_t1_cb
    Ax_A_A::Function       = Ax_A_A_avr_t1_cb
    Bx_S::Function           = Bx_S_avr_t1_cb
    Cx_S::Function           = Cx_S_avr_t1_cb
    Bx_A::Function           = Bx_A_avr_t1_cb
    Cx_A::Function           = Cx_A_avr_t1_cb
    
    stab_Ac::Function = stab_Ac_avr_gen_t1_cb
        
end


function avr_t1_cb_sauer!(dx, x, p_agg, t)
    
    cb_sw, u_gen_ur_ui, f_t, p = p_agg

    sw = cb_sw[1]
    
    u_r        = u_gen_ur_ui[1]
    
    u_i        = u_gen_ur_ui[2]

    vh         = abs(u_r + im * u_i)
    
    v_ref0     = f_t[1]
    
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
       
    # vr1 = x[1], vr2 = x[2], vf_tilade = x[3], v_ref = x[4]
    
    # VRi, Rfi,  Efdi,     Vrefi
    
    vr1,   vr2, vf_tilade, v_ref, vf  = x   
    
    if sw == 0

        dx[1] = (-vr1 + Ka * vr2 - (Ka * Kf/Tf) * vf_tilade  + Ka * (v_ref0 - vh)) /Ta
                
    elseif sw == 1
        
        dx[1] = 0.0
        
    elseif sw == 2
        
        dx[1] = 0.0
    end

    dx[2] = (1/Tf) * (-vr2 + (Kf/Tf) * vf_tilade  )

    dx[3] = (1/Te) * (-vf_tilade * (Ke + Sevf(Ae, Be, vf_tilade) ) +  vr1 )    

    dx[4] = v_ref0 - v_ref

    dx[5] = vf_tilade - vf
    
    nothing
end



@kwdef struct avr_t1_cb_sauer <: SdAvr
    
    Ta::Float64 # Amplifier time constant
    Tf::Float64 # Stabilizer time constant
    Te::Float64 # Field circuit time constant
    Tr::Float64 # Measurement time constant
    Ka::Float64 # Amplifier gain
    Ke::Float64 # Field circuit integral deviation
    Kf::Float64 # Stabilizer gain
    V_R_max::Float64 # Maximum regulator voltage
    V_R_min::Float64 # Minimum regulator voltage
    # S_E_max::Float64 # 
    # S_E0_75max::Float64
    Ae::Float64 # 1st ceiling coeﬃcient
    Be::Float64 # 2nd ceiling coeﬃcient

    comp_type::Symbol = :avr_t1_cb_sauer

    state_vars_syms::Vector{Symbol} = Symbol[ :vr1, :vr2, :vf_tilade]
    algebraic_vars_syms::Vector{Symbol} = Symbol[:v_ref, :vf]
        
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[ :vf ]

    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )       
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]


    dict_im_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( state_vars_syms  ) )

    dict_im_algebraic_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_algebraic_vars_syms ) ) 
    
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    param::Vector{Symbol} = Symbol[:Ta,:Te,:Tf,:Tr,:Ka,:Ke,:Kf,:V_R_max,:V_R_min, :Ae,:Be]
    param_values::Vector{Float64} = Float64[ Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be ]
    func::Vector{Function}            = Function[avr_t1_cb_sauer!]
    control_sig_syms::Vector{Symbol}  = Symbol[:v_ref0]
    control_sig::Vector{Float64}      = ones(length(control_sig_syms))
    # output_sig_syms::Vector{Symbol}   = Symbol[:vf_tilade]
    output_sig_syms::Vector{Symbol}   = Symbol[:vf]
    output_sig::Vector{Float64}       = ones(length(output_sig_syms))

    output_sig_state_coeff::Vector{Float64} = [0.0, 0.0, 1.0 ]
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]     
    cb_state_values::Vector{Function}          = Function[]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)
    cb_dyn_state_event_func::Vector{Function}  = Function[ fun_condition_lower_lim, fun_condition_upper_lim ]
    cb_dyn_state_affect_func::Vector{Function} = Function[ fun_affect_anti_windup_lower_lim!, fun_affect_anti_windup_upper_lim!]
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[ :vr1, :vr1 ]
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
    cb_dyn_state_conditions::Vector{Float64}   = Float64[ V_R_min, V_R_max ]  
    cb_dyn_state_values::Vector{Function}      = Function[(u,t) -> V_R_min, (u, t) -> V_R_max ]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[0, 0]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)
    
    # Ax::Function = Ax_avr_t1_cb_sauer

    # Ax_S_Sgen::Function = Ax_S_Sgen_avr_t1_cb_sauer

    # Bx::Function = Bx_avr_t1_cb_sauer

    # Cx::Function = Cx_avr_t1_cb_sauer

    # Ac_S_with_A::Function = Ac_Savr_with_Aavr_avr_t1_cb_sauer

    # Ac_A_with_S::Function = Ac_Aavr_with_Savr_avr_t1_cb_sauer

    # Ax_A::Function = Ax_Aavr_avr_t1_cb_sauer

    # Ac_A_with_Sgen::Function = Ac_Aavr_with_Sgen_avr_t1_cb_sauer

    # Bx_A::Function = Bx_Aavr_avr_t1_cb_sauer

    # Cx_A::Function = Cx_Aavr_avr_t1_cb_sauer


    Ax_S_Sgen::Function = Ax_S_Sgen_avr_t1_cb_sauer
    Ax_A_Sgen::Function = Ax_A_Sgen_avr_t1_cb_sauer
    Ax_S_S::Function       = Ax_S_S_avr_t1_cb_sauer
    Ax_S_A::Function       = Ax_S_A_avr_t1_cb_sauer
    Ax_A_S::Function       = Ax_A_S_avr_t1_cb_sauer
    Ax_A_A::Function       = Ax_A_A_avr_t1_cb_sauer
    Bx_S::Function           = Bx_S_avr_t1_cb_sauer
    Cx_S::Function           = Cx_S_avr_t1_cb_sauer
    Bx_A::Function           = Bx_A_avr_t1_cb_sauer
    Cx_A::Function           = Cx_A_avr_t1_cb_sauer
    
    stab_Ac::Function = stab_Ac_avr_gen_t1_cb_sauer
        
end


function avr_t1_with_pss_cb!(dx, x, p_agg, t)
    
    cb_sw, u_pss_vs, u_gen_ur_ui, f_t, p = p_agg
    
    sw = cb_sw[1]
    
    u_r        = u_gen_ur_ui[1]
    
    u_i        = u_gen_ur_ui[2]

    vh         = abs(u_r + im * u_i)
    
    v_ref0     = f_t[1]
    
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
    
    vm, vr1, vr2, vf_tilade, v_ref, vf = x    

    dx[1] = (1/Tr)  * (vh - vm )

    
    if sw == 0
        
        dx[2] = (Ka * ((v_ref0 + u_pss_vs - vm) + vr2 - (Kf/Tf) * vf_tilade ) - vr1) / Ta
        
        # dx[2] = (Ka * ((v_ref - vm) + vr2 - (Kf/Tf) * vf_tilade ) - vr1) / Ta
        
    elseif sw == 1
        
        dx[2] = 0.0
        
    elseif sw == 2
        
        dx[2] = 0.0
    end    
     
    dx[3] = (1/Tf) * ((Kf/Tf) * vf_tilade - vr2 )

    dx[4] = (1/Te) * ( vr1 - vf_tilade * (Ke + Sevf(Ae, Be, vf_tilade)) )

    dx[5] = v_ref0  + u_pss_vs - v_ref

    dx[6] = vf_tilade - vf
    
    nothing
end


@kwdef struct avr_t1_with_pss_cb <: SdAvr
    
    Ta::Float64 # Amplifier time constant
    Tf::Float64 # Stabilizer time constant
    Te::Float64 # Field circuit time constant
    Tr::Float64 # Measurement time constant
    Ka::Float64 # Amplifier gain
    Ke::Float64 # Field circuit integral deviation
    Kf::Float64 # Stabilizer gain
    V_R_max::Float64 # Maximum regulator voltage
    V_R_min::Float64 # Minimum regulator voltage
    # S_E_max::Float64 # 
    # S_E0_75max::Float64
    Ae::Float64 # 1st ceiling coeﬃcient
    Be::Float64 # 2nd ceiling coeﬃcient

    comp_type::Symbol = :avr_t1_with_pss_cb

    state_vars_syms::Vector{Symbol} = Symbol[ :vm, :vr1, :vr2, :vf_tilade ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :v_ref, :vf ]
        
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[ :vf ]

    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )     

    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]

    dict_im_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( state_vars_syms  ) )

    dict_im_algebraic_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_algebraic_vars_syms  ) ) 
    
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    param::Vector{Symbol}  = Symbol[ :Ta,:Te,:Tf,:Tr,:Ka,:Ke,:Kf,:V_R_max,:V_R_min, :Ae,:Be]
    param_values::Vector{Float64} = Float64[ Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be ]
    func::Vector{Function}            = Function[avr_t1_with_pss_cb!]
    control_sig_syms::Vector{Symbol}  = Symbol[:v_ref0]
    control_sig::Vector{Float64}      = ones(length(control_sig_syms))
    # output_sig_syms::Vector{Symbol}   = Symbol[:vf_tilade]
    output_sig_syms::Vector{Symbol}   = Symbol[:vf]
    output_sig::Vector{Float64}       = ones(length(output_sig_syms))
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    cb_state_event_func::Vector{Function}      = Function[]
    cb_state_affect_func::Vector{Function}     = Function[]    
    cb_state_syms::Vector{Symbol}              = Symbol[]           
    cb_state_conditions::Vector{Float64}       = Float64[]     
    cb_state_values::Vector{Function}          = Function[]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)
    cb_dyn_state_event_func::Vector{Function}  = Function[ fun_condition_lower_lim, fun_condition_upper_lim ]
    cb_dyn_state_affect_func::Vector{Function} = Function[ fun_affect_anti_windup_lower_lim!, fun_affect_anti_windup_upper_lim!]
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[ :vr1, :vr1 ]
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
    cb_dyn_state_conditions::Vector{Float64}   = Float64[ V_R_min, V_R_max ]  
    cb_dyn_state_values::Vector{Function}      = Function[(u,t) -> V_R_min, (u, t) -> V_R_max ]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[0, 0]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    Ax_S_Sgen::Function = anonymus_func
    Ax_A_Sgen::Function = anonymus_func
    Ax_S_S::Function       = anonymus_func
    Ax_S_A::Function       = anonymus_func
    Ax_A_S::Function       = anonymus_func
    Ax_A_A::Function       = anonymus_func
    Bx_S::Function           = anonymus_func
    Cx_S::Function           = anonymus_func
    Bx_A::Function           = anonymus_func
    Cx_A::Function           = anonymus_func
    
    stab_Ac::Function = anonymus_func
    
    
end



function pss_t2_cb!(dx, x, p_agg, t)
    
    cb_sw, u_gen_ω, f_t,  p = p_agg
    
    T1     = p[1]
    T2     = p[2]
    T3     = p[3]
    T4     = p[4]
    Tw     = p[5]
    Kw     = p[6]
    vs_max = p[7]
    vs_min = p[8]

    ω_ref0 = f_t[1]
    
    
    ωs     = 1.0
    
    v_SI   =  (ω_ref0 -  u_gen_ω)

    v1, v2, v3, vs = x

    # v1  = x[1], v2 = x[2], v3 = x[3], vs = x[4]

    dx[1] = (-1/Tw) * (Kw * v_SI + v1)
    
    dx[2] = (1/T2) * ((1 - T1/T2)*(Kw * v_SI + v1) - v2)
    
    dx[3] = (1/T4) * ((1 - T3/T4)*(v2 + ((T1/T2) * (Kw * v_SI + v1))) - v3)

    dx[4] = v3 + (T3/T4) * (v2 + (T1/T2) * (Kw * v_SI + v1)) - vs 
    
    nothing
end

@kwdef struct pss_t2_cb <: SdPss
    
    T1::Float64     # First stabilizer time constant
    T2::Float64     # Second stabilizer time constant
    T3::Float64     # Third stabilizer time constant
    T4::Float64     # Fourth stabilizer time constant    
    Tw::Float64     # Wash-out time constant
    Kw::Float64     # Stabilizer gain
    vs_max::Float64 # Max stabilizer output signal
    vs_min::Float64 # Min stabilizer output signal

    comp_type::Symbol = :pss_t2_cb
    
    state_vars_syms::Vector{Symbol} = Symbol[ :v1, :v2, :v3]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :vs]
    
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]


    dict_im_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( state_vars_syms  ) )

    dict_im_algebraic_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( algebraic_vars_syms  ) ) 
    
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
    param::Vector{Symbol} = Symbol[ :T1,:T2,:T3,:T4,:Tw,:Kw,:vs_max,:vs_min ]
    param_values::Vector{Float64} = Float64[ T1, T2, T3, T4, Tw, Kw, vs_max, vs_min ]
    func::Vector{Function} = Function[pss_t2_cb!]
    control_sig_syms::Vector{Symbol} = Symbol[:ω_ref0]
    control_sig::Vector{Float64}     = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol}  = Symbol[:vs ] 
    output_sig::Vector{Float64}      = ones(length(output_sig_syms))
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    cb_state_event_func::Vector{Function}      = Function[fun_condition_lower_lim, fun_condition_upper_lim]
    cb_state_affect_func::Vector{Function}     = Function[fun_affect_windup_lower_lim!, fun_affect_windup_upper_lim!]    
    cb_state_syms::Vector{Symbol}              = Symbol[ :vs, :vs ]           
    cb_state_conditions::Vector{Float64}       = Float64[ vs_min, vs_max ]     
    cb_state_values::Vector{Function}          = Function[(u,t) -> vs_min, (u,t)-> vs_max]
    cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
    cb_state_dim::Int64                        = length(cb_state_conditions)
    cb_dyn_state_event_func::Vector{Function}  = Function[]
    cb_dyn_state_affect_func::Vector{Function} = Function[]    
    cb_dyn_state_syms::Vector{Symbol}          = Symbol[]
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
    cb_dyn_state_conditions::Vector{Float64}   = Float64[]  
    cb_dyn_state_values::Vector{Function}      = Function[]
    cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[0]
    cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)

    Ax_S_Sgen::Function = anonymus_func
    Ax_A_Sgen::Function = anonymus_func
    Ax_S_S::Function       = anonymus_func
    Ax_S_A::Function       = anonymus_func
    Ax_A_S::Function       = anonymus_func
    Ax_A_A::Function       = anonymus_func
    Bx_S::Function           = anonymus_func
    Cx_S::Function           = anonymus_func
    Bx_A::Function           = anonymus_func
    Cx_A::Function           = anonymus_func
    
    stab_Ac::Function = anonymus_func
    
end


# function avr_ieee_dc1a_cb!(dx, x, p_agg, t)
    
#     cb_sw, u_gen_ur_ui, f_t, p = p_agg

#     sw = cb_sw[1]
    
#     u_r        = u_gen_ur_ui[1]
    
#     u_i        = u_gen_ur_ui[2]

#     vh         = abs(u_r + im * u_i)
    
#     v_ref0     = f_t[1]
    
#     Ta         = p[1]
#     Te         = p[2]
#     Tf         = p[3]
#     Tr         = p[4]
#     Ka         = p[5]
#     Ke         = p[6]
#     Kf         = p[7]
#     V_R_max    = p[8]
#     V_R_min    = p[9]
#     Ae         = p[10]
#     Be         = p[11]
       
#     # vm = x[1], vr1 = x[2], vr2 = x[3], vf_tilade = x[4], v_ref = x[5], 
     
#     vm, vr1, vr2, vf_tilade, v_ref  = x

#     dx[1] = (1/Tr)  * (vh - vm )
    
#     if sw == 0
        
#         dx[2] = (Ka * ((v_ref0 - vm) + vr2 - (Kf/Tf) * vf_tilade ) - vr1) /Ta
        
#     elseif sw == 1
        
#         dx[2] = 0.0
        
#     elseif sw == 2
        
#         dx[2] = 0.0
#     end

#     dx[3] = -vr2/Tf + (Kf/Tf^2) * vf_tilade

#     dx[4] = (-1/Te) * (vf_tilade * (Ke + Sevf(Ae, Be, vf_tilade)) - vr1)    

#     dx[5] = v_ref0 - v_ref
    
#     nothing
# end


# @kwdef struct avr_ieee_dc1a_cb <: SdAvr
    
#     Ta::Float64 # Amplifier time constant
#     Tf::Float64 # Stabilizer time constant
#     Te::Float64 # Field circuit time constant
#     Tr::Float64 # Measurement time constant
#     Ka::Float64 # Amplifier gain
#     Ke::Float64 # Field circuit integral deviation
#     Kf::Float64 # Stabilizer gain
#     V_R_max::Float64 # Maximum regulator voltage
#     V_R_min::Float64 # Minimum regulator voltage
#     # S_E_max::Float64 # 
#     # S_E0_75max::Float64
#     Ae::Float64 # 1st ceiling coeﬃcient
#     Be::Float64 # 2nd ceiling coeﬃcient

#     comp_type::Symbol = :avr_ieee_dc1a_cb

#     state_vars_syms::Vector{Symbol} = Symbol[ :vm, :vr1, :vr2, :vf_tilade]
#     algebraic_vars_syms::Vector{Symbol} = Symbol[ :v_ref ]
#     syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
#     state_dim::Int64 = length(state_vars_syms) 
#     dim::Int64 = length(syms)
#     mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
#     dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms) )
#     param::Vector{Symbol}             = Symbol[ :Ta,:Te,:Tf,:Tr,:Ka,:Ke,:Kf,:V_R_max,:V_R_min, :Ae,:Be ]
#     param_values::Vector{Float64} = Float64[ Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be ]
#     func::Vector{Function}            = Function[avr_ieee_dc1a_cb!]
#     control_sig_syms::Vector{Symbol}  = Symbol[:v_ref]
#     control_sig::Vector{Float64}      = ones(length(control_sig_syms))
#     output_sig_syms::Vector{Symbol}   = Symbol[:vf_tilade]
#     output_sig::Vector{Float64}       = ones(length(output_sig_syms))
#     dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
#     cb_state_event_func::Vector{Function}      = Function[]
#     cb_state_affect_func::Vector{Function}     = Function[]    
#     cb_state_syms::Vector{Symbol}              = Symbol[]           
#     cb_state_conditions::Vector{Float64}       = Float64[]     
#     cb_state_values::Vector{Function}          = Function[]
#     cb_state_sym2Idx::Vector{Int64}            = Int64[dict_state_syms[state_sym] for state_sym in  cb_state_syms ]
#     cb_state_dim::Int64                        = length(cb_state_conditions)
#     cb_dyn_state_event_func::Vector{Function}  = Function[ fun_condition_lower_lim, fun_condition_upper_lim ]
#     cb_dyn_state_affect_func::Vector{Function} = Function[ fun_affect_anti_windup_lower_lim!, fun_affect_anti_windup_upper_lim! ]
#     cb_dyn_state_syms::Vector{Symbol}          = Symbol[ :vr1, :vr1 ]
#     cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
#     cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]
#     cb_dyn_state_conditions::Vector{Float64}   = Float64[ V_R_min, V_R_max ]  
#     cb_dyn_state_values::Vector{Function}      = Function[(u,t) -> V_R_min, (u, t) -> V_R_max ]
#     cb_dyn_state_sym2Idx::Vector{Int64}        = Int64[dict_state_syms[state_sym] for state_sym in  cb_dyn_state_syms ]
#     cb_dyn_param_state_sw::Vector{Int64}       = Int64[0, 0]
#     cb_dyn_state_dim::Int64                    = length(cb_dyn_state_conditions)
        
# end

