
# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


########################################################
# ------------------------------------------------------
#  Static Edges Components Structures
# ------------------------------------------------------
########################################################

"""
Generic nodal struct for power flow
"""


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



########################################################
# ------------------------------------------------------
#  Dynamics Edges Components Structures
# ------------------------------------------------------
########################################################


function PiModelLine!(de, e, p_agg, t)

    v_s, v_d, p = p_agg
    
    y          = p[1]
    y_shunt_km = p[2]
    y_shunt_mk = p[3]
    
    y_shunt_from_bus = p[4][1]
    y_shunt_to_bus   = p[4][2]
    
    y_shunt_km_prime = y_shunt_km + y_shunt_from_bus
    y_shunt_mk_prime = y_shunt_mk + y_shunt_to_bus
    
    uh = v_s[1] + v_s[2] * im
    
    uk = v_d[1] + v_d[2] * im
    
    u  = [uh, uk]
    
    # Y  = PiModel(y, y_shunt_km, y_shunt_mk, 1, 1)    
    # Y = calc_branch_Ybr(y, y_shunt_km, y_shunt_mk, 1, 1)
    # Y  = PiModel(y, y_shunt_km_prime, y_shunt_mk_prime, 1, 1)
    Y = calc_branch_Ybr(y, y_shunt_km_prime, y_shunt_mk_prime, 1.0, 1.0)
        
    i  = Y * u
    
    de[1] = real(i[1]) - e[1]  # ih = i[1] 
    de[2] = imag(i[1]) - e[2]    
    de[3] = real(i[2]) - e[3]  # ik = i[2]
    de[4] = imag(i[2]) - e[4]
    
    return nothing
end


@kwdef struct PiModelLine  <: SdBranchElement
    from::String
    to::String
    y::ComplexF64
    y_shunt_km::Float64
    y_shunt_mk::Float64

    comp_type::Symbol = :PiModelLine

    # y_bus_g_sh::Float64 = 0.0
    # y_bus_b_sh::Float64 = 0.0

    # # y_bus_shunts::Tuple{Float64, Float64} = (0.0, 0.0)
    
    # y_bus_shunts::Tuple{Float64, Float64} = (y_bus_g_sh, y_bus_b_sh)
    
    name::String = from * "-" * to
    orientation::Tuple{Int64, Int64} = (parse(Int, split(lowercase(from),"bus")[2] ), parse(Int, split(lowercase(to),"bus")[2]) )
    state_vars_syms::Vector{Symbol} = Symbol[ ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[:i_s_r, :i_s_i, :i_d_r, :i_d_i ]
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    # param::Vector{Symbol} = Symbol[ :y, :y_shunt_km, :y_shunt_mk, :y_bus_shunts ]
    
    # param_values::Vector{Union{Float64, ComplexF64, Tuple{Float64, Float64} }} = Union{Float64, ComplexF64, Tuple{Float64, Float64}}[ y, y_shunt_km, y_shunt_mk, y_bus_shunts ]
    
    # param_real_sym::Vector{Symbol} = Symbol[ :g, :b, :y_shunt_km, :y_shunt_mk, :y_bus_shunts ]
    
    # param_real_values::Vector{ Float64} = Float64[ real(y), imag(y), y_shunt_km, y_shunt_mk, y_bus_shunts ]
    
    # param_real_values::Vector{ Float64} = Float64[ real(y), imag(y), y_shunt_km, y_shunt_mk, y_bus_shunts ]
    
    param::Vector{Symbol} = Symbol[ :y, :y_shunt_km, :y_shunt_mk ]
    
    param_values::Vector{Union{Float64, ComplexF64 }} = Union{Float64, ComplexF64 }[ y, y_shunt_km, y_shunt_mk ]
    
    param_real_sym::Vector{Symbol} = Symbol[ :g, :b, :y_shunt_km, :y_shunt_mk ]
    
    param_real_values::Vector{ Float64} = Float64[ real(y), imag(y), y_shunt_km, y_shunt_mk ]
    
    func::Vector{Function} = Function[PiModelLine!]
    control_sig_syms::Vector{Symbol} = Symbol[]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[ ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    
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
    
    
end


function StaticLine!(de, e, p_agg, t)

    v_s, v_d, p = p_agg
    
    Y = p[1]
    
    uh = v_s[1] + v_s[2]*im
    uk = v_d[1] + v_d[2]*im
    
    i = Y * (uh - uk)
    
    de[1] = real(i)    - e[1]
    de[2] = imag(i)    - e[2] 
    de[3] = - real(i)  - e[3]
    de[4] = - imag(i)  - e[4]
    
    return nothing
end

@kwdef struct StaticLine  <: SdBranchElement
    from::String
    to::String
    G::Float64
    B::Float64
    Y::ComplexF64 = G + im * B
    
    comp_type::Symbol = :StaticLine
    
    name::String = from * "-" * to
    orientation::Tuple{Int64, Int64} = (parse(Int, split(lowercase(from),"bus")[2] ), parse(Int, split(lowercase(to),"bus")[2]) )
    
    state_vars_syms::Vector{Symbol} = Symbol[ ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[:i_s_r, :i_s_i, :i_d_r, :i_d_i ]
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[ :Y ]
    param_values::Vector{Union{Float64, ComplexF64}} = Union{Float64, ComplexF64}[ Y ]

    param_real_sym::Vector{Symbol} = Symbol[ :g, :b ]
    param_real_values::Vector{ Float64} = Float64[ real(Y), imag(Y) ]
    
    func::Vector{Function} = Function[StaticLine!]
    control_sig_syms::Vector{Symbol} = Symbol[]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[ ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    
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
    
    
end


function RLLine!(de, e, p_agg, t)

    v_s, v_d, p = p_agg
    
    R = p[1]
    L = p[2]
    ω0 = p[3]

    @smart_assert R > 0 "The resistance (R) should be greater than  zero (R > 0). This should fail because R = $(R)"
    
    @smart_assert L > 0  "The inductance (L) should be greater than  zero (L > 0). This should fail because L = $(L)"    

    @smart_assert  ω0 > 0  "The rated frequency (ω0)  should be greater than  zero (ω0 > 0). This should fail because ω0 = $(ω0)"
    
    Z = [R -ω0*L; ω0*L R]

    di = (-(Z * [e[1]; e[2]]) .+ [v_d[1]; v_d[2]] .- [v_s[1]; v_s[2]]) ./ L
    de .= [di; -di]

    return nothing
end


@kwdef struct RLLine <: SdBranchElement
    from::String
    to::String
    R::Float64
    L::Float64
    ω0::Float64

    comp_type::Symbol = :RLLine
    
    name::String = from * "-" * to
    orientation::Tuple{Int64, Int64} = (parse(Int, split(lowercase(from),"bus")[2] ), parse(Int, split(lowercase(to),"bus")[2]) )     
    state_vars_syms::Vector{Symbol} = Symbol[ ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[:i_s_r, :i_s_i, :i_d_r, :i_d_i ]
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    param::Vector{Symbol} = Symbol[ :R, :L, :ω0]
    param_values::Vector{Union{Float64, ComplexF64}} = Union{Float64, ComplexF64}[ R, L, ω0  ]

    param_real_sym::Vector{Symbol} = Symbol[ :R, :L, :ω0 ]
    param_real_values::Vector{ Float64} = Float64[ R, L, ω0  ]
        
    func::Vector{Function} = Function[RLLine!]
    control_sig_syms::Vector{Symbol} = Symbol[]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[ ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    
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
    
end


function Transformer!(de, e, p_agg, t)

    v_s, v_d, p = p_agg
    
    y = p[1]
    y_shunt_km = p[2]
    y_shunt_mk = p[3]
    t_ratio    = p[4]
    
    y_shunt_from_bus = p[5][1]
    y_shunt_to_bus   = p[5][2]
    
    y_shunt_km_prime = y_shunt_km + y_shunt_from_bus
    y_shunt_mk_prime = y_shunt_mk + y_shunt_to_bus 
    
    uh = v_s[1] + v_s[2]*im
    uk = v_d[1] + v_d[2]*im
    
    u = [uh, uk]

    
    # Y = PiModel(y, y_shunt_km, y_shunt_mk, t_ratio, 1)
    # Y = calc_branch_Ybr(y,y_shunt_km,y_shunt_mk,t_ratio,1)
    # Y = PiModel(y, y_shunt_km_prime, y_shunt_mk_prime, t_ratio, 1)

    Y = calc_branch_Ybr(y, y_shunt_km_prime, y_shunt_mk_prime, t_ratio, 1)
 
    
    i = Y * u
    
    de[1] = real(i[1]) - e[1]
    de[2] = imag(i[1]) - e[2]    
    de[3] = real(i[2]) - e[3]
    de[4] = imag(i[2]) - e[4]

    return nothing
end


@kwdef struct Transformer <: SdBranchElement
    from::String
    to::String
    y::ComplexF64 
    # y_shunt_km::Float64 = 0 # additional
    # y_shunt_mk::Float64 = 0 # additional
    y_shunt_km::Float64
    y_shunt_mk::Float64
    t_ratio::Float64

    comp_type::Symbol = :Transformer

    # y_bus_shunts::Tuple{Float64, Float64} = (0.0, 0.0)

    # y_bus_g_sh::Float64 = 0.0
    # y_bus_b_sh::Float64 = 0.0

    # y_bus_shunts::Tuple{Float64, Float64} = (y_bus_g_sh, y_bus_b_sh)
    
    name::String = from * "-" * to
    orientation::Tuple{Int64, Int64} = (parse(Int, split(lowercase(from),"bus")[2] ), parse(Int, split(lowercase(to),"bus")[2]) ) 
    state_vars_syms::Vector{Symbol} = Symbol[ ]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ :i_s_r, :i_s_i, :i_d_r, :i_d_i ]
    state_dim::Int64 = length(state_vars_syms) 
    syms::Vector{Symbol} = Symbol[ state_vars_syms...;algebraic_vars_syms...]
    dim::Int64 = length(syms) 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = DAE_MassMatrix(length(state_vars_syms), length(algebraic_vars_syms))
    dae_var::Vector{Bool} = DAE_BoolVector( length(state_vars_syms), length(algebraic_vars_syms)   )
    
    # param::Vector{Symbol} = Symbol[:y,:y_shunt_km, :y_shunt_mk, :t_ratio, :y_bus_shunts ]
    
    # param_values::Vector{Union{Float64, ComplexF64, Tuple{Float64, Float64}}} = Union{Float64, ComplexF64, Tuple{Float64, Float64}}[ y, y_shunt_km, y_shunt_mk, t_ratio, y_bus_shunts]

    # param_real_sym::Vector{Symbol} = Symbol[:g,:b,:y_shunt_km,:y_shunt_mk,:t_ratio, :y_bus_shunts ]
    
    # param_real_values::Vector{Float64} = Float64[real(y), imag(y), y_shunt_km, y_shunt_mk, t_ratio, y_bus_shunts ]


    param::Vector{Symbol} = Symbol[:y,:y_shunt_km, :y_shunt_mk, :t_ratio ]
    
    param_values::Vector{Union{Float64, ComplexF64 }} = Union{Float64, ComplexF64 }[ y, y_shunt_km, y_shunt_mk, t_ratio]

    param_real_sym::Vector{Symbol} = Symbol[:g,:b,:y_shunt_km,:y_shunt_mk,:t_ratio ]
    
    param_real_values::Vector{Float64} = Float64[real(y), imag(y), y_shunt_km, y_shunt_mk, t_ratio ]
    
    func::Vector{Function} = Function[Transformer!]
    control_sig_syms::Vector{Symbol} = Symbol[]
    control_sig::Vector{Float64} = ones(length(control_sig_syms))
    output_sig_syms::Vector{Symbol} = Symbol[ ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    
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
    
end


