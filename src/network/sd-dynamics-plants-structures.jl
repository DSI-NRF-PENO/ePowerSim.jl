# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


#########################################################
# ------------------------------------------------------
#  composite plants
# ------------------------------------------------------
#########################################################

#--------------------------------------------------------
# load plants
#--------------------------------------------------------


function hybrid_pf_plant_Load!(dx, x, (p_agg, plant), t)
    
    # cb_sw, src_i, dst_i, f_t, p, nodes_idx_and_inc_edges, nodes_inc_edges_Ybr_orient, nodes_u_view =  p_agg

    cb_sw, src_i, dst_i, f_t, p, global_pf_param, node_pf_param =  p_agg

    load = plant.Load

    load_states_syms = load.syms
    dict_load_syms   = load.dict_state_syms
    load_output_sigs = load.output_sig_syms

    load_fun! = load.func[5]

    load_fun!(dx, x,
             (cb_sw, src_i, dst_i, f_t, p, global_pf_param,  node_pf_param ), t)    
        
    return nothing
end


function node_pf_plant_Load!(dx, x, (p_agg, plant), t)
    
    # cb_sw, src_i, dst_i, f_t, p, nodes_idx_and_inc_edges, nodes_inc_edges_Ybr_orient, nodes_u_view =  p_agg

    cb_sw, src_i, dst_i, f_t, p,  node_pf_param =  p_agg

    load = plant.Load

    load_states_syms = load.syms
    dict_load_syms   = load.dict_state_syms
    load_output_sigs = load.output_sig_syms

    load_fun! = load.func[4]

    load_fun!(dx, x,
             (cb_sw, src_i, dst_i, f_t, p, node_pf_param ), t)    
        
    return nothing
end


function global_pf_plant_Load!(dx, x, (p_agg, plant), t)
    
    # cb_sw, src_i, dst_i, f_t, p, pf_U, Inet =  p_agg
    
    cb_sw, src_i, dst_i, f_t, p, global_pf_param =  p_agg

    load = plant.Load

    load_states_syms = load.syms
    dict_load_syms   = load.dict_state_syms
    load_output_sigs = load.output_sig_syms

    load_fun! = load.func[3]

    load_fun!(dx, x,
             (cb_sw, src_i, dst_i, f_t, p, global_pf_param ), t)    
        
    return nothing
end



function network_current_plant_Load!(dx, x, (p_agg, plant), t)
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    load = plant.Load

    load_states_syms = load.syms
    dict_load_syms   = load.dict_state_syms
    load_output_sigs = load.output_sig_syms

    load_fun! = load.func[2]

    load_fun!(dx, x,
             (cb_sw, src_i, dst_i, f_t, p), t)    
        
    return nothing
end



function initial_pf_plant_Load!(dx, x, (p_agg, plant), t)
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    load = plant.Load

    load_states_syms = load.syms
    dict_load_syms   = load.dict_state_syms
    load_output_sigs = load.output_sig_syms

    load_fun! = load.func[1]


    load_fun!(dx, x,
             (cb_sw, src_i, dst_i, f_t, p), t)    
        
    return nothing
end


@kwdef struct plant_PQ_Const_P <: SdNonGenPlant
    Load::PQ_Const_P
    Bus::String      = Load.Bus
    name::String     = Load.name
    Bus_num::Int64   = Load.Bus_num
    Bus_type::Symbol = Load.Bus_type

    comp_type::Symbol = :plant_PQ_Const_P

    with_loc_load::Bool = false

    isa_slack::Bool = false
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[Load.algebraic_vars_syms...;]
    state_vars_syms::Vector{Symbol} = Symbol[Load.state_vars_syms...; ]
    syms::Vector{Symbol} = Symbol[state_vars_syms...; algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]     
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)    
    dims::Vector{Int64} = Int64[ Load.dim ]
    dae_var::Vector{Bool} =  Bool[Load.dae_var...;] 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    param::Vector{ Vector{Symbol} } = Vector{Symbol}[ Load.param ]
    param_values::Vector{ Vector{Union{Int64, ComplexF64, Float64, Tuple{Float64, Float64, Float64}}} } = Vector{Union{Int64, ComplexF64, Float64, Tuple{Float64, Float64, Float64}}}[ Load.param_values ]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [ param...; ] ) )


    func::Vector{Function} = Function[ initial_pf_plant_Load!, network_current_plant_Load!, global_pf_plant_Load!, node_pf_plant_Load!,  hybrid_pf_plant_Load! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[  ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms)) 
    output_sig_syms::Vector{Symbol} = Symbol[ Load.output_sig_syms...; ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))     
    cb_state_event_func::Vector{Function} = Function[ Load.cb_state_event_func...; ]  
    cb_state_affect_func::Vector{Function} = Function[ Load.cb_state_affect_func...;  ]    
    cb_state_syms::Vector{Symbol} = Symbol[ Load.cb_state_syms...;  ]           
    cb_state_conditions::Vector{Float64} = Float64[ Load.cb_state_conditions...; ]     
    cb_state_values::Vector{Function} = Function[ Load.cb_state_values...; ]
    cb_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[ state_sym ] for state_sym in cb_state_syms ]
    
    cb_state_dim::Int64  = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[Load.cb_dyn_state_event_func...; ]
    cb_dyn_state_affect_func::Vector{Function} = Function[Load.cb_dyn_state_affect_func...; ]    
    cb_dyn_state_syms::Vector{Symbol} = Symbol[ Load.cb_dyn_state_syms...; ]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ Load.cb_dyn_state_sw_Idx...; ]
            
    cb_dyn_state_conditions::Vector{Float64} = Float64[ Load.cb_dyn_state_conditions...; ]  
    cb_dyn_state_values::Vector{Function} = Function[ Load.cb_dyn_state_values...; ]
    cb_dyn_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[state_sym] for state_sym in cb_dyn_state_syms]
    cb_dyn_param_state_sw::Vector{ Vector{Int64} } = Vector{Int64}[ Load.cb_dyn_param_state_sw ]
    cb_dyn_state_dim::Int64 = length(cb_dyn_state_conditions) 
    
end


@kwdef struct plant_PQ_Const_I <: SdNonGenPlant
    Load::PQ_Const_I
    Bus::String      = Load.Bus
    name::String     = Load.name

    Bus_num::Int64   = Load.Bus_num
    Bus_type::Symbol = Load.Bus_type

    comp_type::Symbol = :plant_PQ_Const_I

    with_loc_load::Bool = false

    isa_slack::Bool = false
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[Load.algebraic_vars_syms...;]
    state_vars_syms::Vector{Symbol} = Symbol[Load.state_vars_syms...; ]
    syms::Vector{Symbol} = Symbol[state_vars_syms...; algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]     
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)    
    dims::Vector{Int64} = Int64[ Load.dim ]
    dae_var::Vector{Bool} =  Bool[Load.dae_var...;] 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    param::Vector{ Vector{Symbol} } = Vector{Symbol}[ Load.param ]
    param_values::Vector{ Vector{Union{Int64, ComplexF64, Float64, Tuple{Float64, Float64, Float64}}} } = Vector{Union{Int64, ComplexF64, Float64, Tuple{Float64, Float64, Float64}}}[ Load.param_values ]
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [ param...; ] ) )
        
    func::Vector{Function} = Function[ initial_pf_plant_Load!, network_current_plant_Load!, global_pf_plant_Load!, node_pf_plant_Load!, hybrid_pf_plant_Load! ]    
    control_sig_syms::Vector{Symbol} = Symbol[  ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms)) 
    output_sig_syms::Vector{Symbol} = Symbol[ Load.output_sig_syms...; ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))     
    cb_state_event_func::Vector{Function} = Function[ Load.cb_state_event_func...; ]  
    cb_state_affect_func::Vector{Function} = Function[ Load.cb_state_affect_func...;  ]    
    cb_state_syms::Vector{Symbol} = Symbol[ Load.cb_state_syms...;  ]           
    cb_state_conditions::Vector{Float64} = Float64[ Load.cb_state_conditions...; ]     
    cb_state_values::Vector{Function} = Function[ Load.cb_state_values...; ]
    cb_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[ state_sym ] for state_sym in cb_state_syms ]
    
    cb_state_dim::Int64  = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[Load.cb_dyn_state_event_func...; ]
    cb_dyn_state_affect_func::Vector{Function} = Function[Load.cb_dyn_state_affect_func...; ]    
    cb_dyn_state_syms::Vector{Symbol} = Symbol[ Load.cb_dyn_state_syms...; ]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ Load.cb_dyn_state_sw_Idx...; ]
            
    cb_dyn_state_conditions::Vector{Float64} = Float64[ Load.cb_dyn_state_conditions...; ]  
    cb_dyn_state_values::Vector{Function} = Function[ Load.cb_dyn_state_values...; ]
    cb_dyn_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[state_sym] for state_sym in cb_dyn_state_syms]
    cb_dyn_param_state_sw::Vector{ Vector{Int64} } = Vector{Int64}[ Load.cb_dyn_param_state_sw ]
    cb_dyn_state_dim::Int64 = length(cb_dyn_state_conditions) 
end



@kwdef struct plant_PQ_Const_Z <: SdNonGenPlant
    Load::PQ_Const_Z
    Bus::String      = Load.Bus
    name::String     = Load.name
    Bus_num::Int64   = Load.Bus_num
    Bus_type::Symbol = Load.Bus_type

    comp_type::Symbol = :plant_PQ_Const_Z

    with_loc_load::Bool = false

    isa_slack::Bool = false
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[Load.algebraic_vars_syms...;]
    state_vars_syms::Vector{Symbol} = Symbol[Load.state_vars_syms...; ]
    syms::Vector{Symbol} = Symbol[state_vars_syms...; algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]     
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)    
    dims::Vector{Int64} = Int64[ Load.dim ]
    dae_var::Vector{Bool} =  Bool[Load.dae_var...;] 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    param::Vector{ Vector{Symbol} } = Vector{Symbol}[ Load.param ]
    param_values::Vector{ Vector{Union{Int64, ComplexF64, Float64, Tuple{Float64, Float64, Float64}}} } = Vector{Union{Int64, ComplexF64, Float64, Tuple{Float64, Float64, Float64}}}[ Load.param_values ]
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [ param...; ] ) )    
    func::Vector{Function} = Function[ initial_pf_plant_Load!, network_current_plant_Load!, global_pf_plant_Load!, node_pf_plant_Load!,  hybrid_pf_plant_Load! ]    
    control_sig_syms::Vector{Symbol} = Symbol[  ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms)) 
    output_sig_syms::Vector{Symbol} = Symbol[ Load.output_sig_syms...; ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))     
    cb_state_event_func::Vector{Function} = Function[ Load.cb_state_event_func...; ]  
    cb_state_affect_func::Vector{Function} = Function[ Load.cb_state_affect_func...;  ]    
    cb_state_syms::Vector{Symbol} = Symbol[ Load.cb_state_syms...;  ]           
    cb_state_conditions::Vector{Float64} = Float64[ Load.cb_state_conditions...; ]     
    cb_state_values::Vector{Function} = Function[ Load.cb_state_values...; ]
    cb_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[ state_sym ] for state_sym in cb_state_syms ]
    cb_state_dim::Int64  = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[Load.cb_dyn_state_event_func...; ]
    cb_dyn_state_affect_func::Vector{Function} = Function[Load.cb_dyn_state_affect_func...; ]    
    cb_dyn_state_syms::Vector{Symbol} = Symbol[ Load.cb_dyn_state_syms...; ]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ Load.cb_dyn_state_sw_Idx...; ]
            
    cb_dyn_state_conditions::Vector{Float64} = Float64[ Load.cb_dyn_state_conditions...; ]  
    cb_dyn_state_values::Vector{Function} = Function[ Load.cb_dyn_state_values...; ]
    cb_dyn_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[state_sym] for state_sym in cb_dyn_state_syms]
    cb_dyn_param_state_sw::Vector{ Vector{Int64} } = Vector{Int64}[ Load.cb_dyn_param_state_sw ]
    cb_dyn_state_dim::Int64 = length(cb_dyn_state_conditions) 
    
end


@kwdef struct plant_PQ_dyn_load <: SdNonGenPlant
    Load::PQ_dyn_load
    Bus::String      = Load.Bus
    name::String     = Load.name

    Bus_num::Int64   = Load.Bus_num
    Bus_type::Symbol = Load.Bus_type

    comp_type::Symbol = :plant_PQ_dyn_load

    with_loc_load::Bool = false

    isa_slack::Bool = false
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[Load.algebraic_vars_syms...;]
    state_vars_syms::Vector{Symbol} = Symbol[Load.state_vars_syms...; ]
    syms::Vector{Symbol} = Symbol[state_vars_syms...; algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]     
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)    
    dims::Vector{Int64} = Int64[ Load.dim ]
    dae_var::Vector{Bool} =  Bool[Load.dae_var...;] 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    param::Vector{ Vector{Symbol} } = Vector{Symbol}[ Load.param ]
    param_values::Vector{ Vector{Union{Int64, ComplexF64, Float64, Tuple{Float64, Float64, Float64}}} } = Vector{Union{Int64, ComplexF64, Float64, Tuple{Float64, Float64, Float64}}}[ Load.param_values ]
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [ param...; ]) )
    
    func::Vector{Function} = Function[ initial_pf_plant_Load!, network_current_plant_Load!, global_pf_plant_Load!, node_pf_plant_Load!,  hybrid_pf_plant_Load! ]    
    control_sig_syms::Vector{Symbol} = Symbol[  ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms)) 
    output_sig_syms::Vector{Symbol} = Symbol[ Load.output_sig_syms...; ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))     
    cb_state_event_func::Vector{Function} = Function[ Load.cb_state_event_func...; ]  
    cb_state_affect_func::Vector{Function} = Function[ Load.cb_state_affect_func...;  ]    
    cb_state_syms::Vector{Symbol} = Symbol[ Load.cb_state_syms...;  ]           
    cb_state_conditions::Vector{Float64} = Float64[ Load.cb_state_conditions...; ]     
    cb_state_values::Vector{Function} = Function[ Load.cb_state_values...; ]
    cb_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[ state_sym ] for state_sym in cb_state_syms ]
    cb_state_dim::Int64  = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[Load.cb_dyn_state_event_func...; ]
    cb_dyn_state_affect_func::Vector{Function} = Function[Load.cb_dyn_state_affect_func...; ]    
    cb_dyn_state_syms::Vector{Symbol} = Symbol[ Load.cb_dyn_state_syms...; ]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ Load.cb_dyn_state_sw_Idx...; ]
            
    cb_dyn_state_conditions::Vector{Float64} = Float64[ Load.cb_dyn_state_conditions...; ]  
    cb_dyn_state_values::Vector{Function} = Function[ Load.cb_dyn_state_values...; ]
    cb_dyn_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[state_sym] for state_sym in cb_dyn_state_syms]
    cb_dyn_param_state_sw::Vector{ Vector{Int64} } = Vector{Int64}[ Load.cb_dyn_param_state_sw ]
    cb_dyn_state_dim::Int64 = length(cb_dyn_state_conditions) 
    
end


#--------------------------------------------------------
# Transmission plants
#--------------------------------------------------------

function hybrid_pf_plant_Transmission!(dx, x, (p_agg, plant), t)
        
    # cb_sw, src_i, dst_i, f_t, p, nodes_idx_and_inc_edges, nodes_inc_edges_Ybr_orient, nodes_u_view =  p_agg

    cb_sw, src_i, dst_i, f_t, p, dyn_global_pf_param, dyn_node_pf_param  =  p_agg

    trans = plant.Trans
    
    trans_f_t = f_t[1] 

    trans_fun! = trans.func[5]

    trans_fun!(dx, x,
             (cb_sw, src_i, dst_i, trans_f_t, p, dyn_global_pf_param, dyn_node_pf_param  ), t)    
        
    return nothing
end


function node_pf_plant_Transmission!(dx, x, (p_agg, plant), t)
        
    # cb_sw, src_i, dst_i, f_t, p, nodes_idx_and_inc_edges, nodes_inc_edges_Ybr_orient, nodes_u_view =  p_agg

    cb_sw, src_i, dst_i, f_t, p, node_pf_param  =  p_agg

    trans = plant.Trans
    
    trans_f_t = f_t[1] 

    trans_fun! = trans.func[4]

    trans_fun!(dx, x,
             (cb_sw, src_i, dst_i, trans_f_t, p, node_pf_param  ), t)    
        
    return nothing
end



function global_pf_plant_Transmission!(dx, x, (p_agg, plant), t)
        
    # cb_sw, src_i, dst_i, f_t, p, pf_U, Inet =  p_agg

    cb_sw, src_i, dst_i, f_t, p, global_pf_param =  p_agg

    trans = plant.Trans
    
    trans_f_t = f_t[1] 


    trans_fun! = trans.func[3]

    trans_fun!(dx, x,
             (cb_sw, src_i, dst_i, trans_f_t, p, global_pf_param  ), t)    
        
    return nothing
end


function network_current_plant_Transmission!(dx, x, (p_agg, plant), t)
        
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    trans = plant.Trans
    
    trans_f_t = f_t[1] 

    trans_fun! = trans.func[2]

    trans_fun!(dx, x,
             (cb_sw, src_i, dst_i, trans_f_t, p), t)    
        
    return nothing
end


function initial_pf_plant_Transmission!(dx, x, (p_agg, plant), t)
        
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    trans = plant.Trans
    
    trans_f_t = f_t[1] 

    trans_fun! = trans.func[1]

    trans_fun!(dx, x,
             (cb_sw, src_i, dst_i, trans_f_t, p), t)    
        
    return nothing
end


@kwdef struct plant_Transmission_t1 <: SdNonGenPlant
    Trans::Trans_t1_Node
    
    Bus::String      = Trans.Bus
    name::String     = Trans.name
    Bus_type::Symbol = Trans.Bus_type
    Bus_num::Int64   = Trans.Bus_num

    comp_type::Symbol =:plant_Transmission_t1

    with_loc_load::Bool = false

    isa_slack::Bool = false
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[Trans.algebraic_vars_syms...;]
    state_vars_syms::Vector{Symbol} = Symbol[Trans.state_vars_syms...; ]
    syms::Vector{Symbol} = Symbol[state_vars_syms...; algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]     
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)    
    dims::Vector{Int64} = Int64[ Trans.dim ]
    dae_var::Vector{Bool} =  Bool[Trans.dae_var...;] 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    param::Vector{Symbol} = Symbol[ Trans.param...;]
    param_values::Vector{ Float64} = Float64[ Trans.param_values...; ]

    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [ param...; ] ) )
    
    func::Vector{Function} = Function[initial_pf_plant_Transmission!, network_current_plant_Transmission!, global_pf_plant_Transmission!, node_pf_plant_Transmission!, hybrid_pf_plant_Transmission! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[  ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms)) 
    output_sig_syms::Vector{Symbol} = Symbol[ Trans.output_sig_syms...; ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))     
    cb_state_event_func::Vector{Function} = Function[ Trans.cb_state_event_func...; ]  
    cb_state_affect_func::Vector{Function} = Function[ Trans.cb_state_affect_func...;  ]    
    cb_state_syms::Vector{Symbol} = Symbol[ Trans.cb_state_syms...;  ]           
    cb_state_conditions::Vector{Float64} = Float64[ Trans.cb_state_conditions...; ]     
    cb_state_values::Vector{Function} = Function[ Trans.cb_state_values...; ]
    cb_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[ state_sym ] for state_sym in cb_state_syms ]
    cb_state_dim::Int64  = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[Trans.cb_dyn_state_event_func...; ]
    cb_dyn_state_affect_func::Vector{Function} = Function[Trans.cb_dyn_state_affect_func...; ]    
    cb_dyn_state_syms::Vector{Symbol} = Symbol[ Trans.cb_dyn_state_syms...; ]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ Trans.cb_dyn_state_sw_Idx...; ]
            
    cb_dyn_state_conditions::Vector{Float64} = Float64[ Trans.cb_dyn_state_conditions...; ]  
    cb_dyn_state_values::Vector{Function} = Function[ Trans.cb_dyn_state_values...; ]
    cb_dyn_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[state_sym] for state_sym in cb_dyn_state_syms]
    cb_dyn_param_state_sw::Vector{ Vector{Int64} } = Vector{Int64}[ Trans.cb_dyn_param_state_sw ]
    cb_dyn_state_dim::Int64 = length(cb_dyn_state_conditions) 
    
end


@kwdef struct plant_Transmission_t2 <: SdNonGenPlant
    Trans::Trans_t2_Node
    Bus::String      = Trans.Bus
    name::String     = Trans.name
    Bus_type::Symbol = Trans.Bus_type
    Bus_num::Int64   = Trans.Bus_num

    comp_type::Symbol = :plant_Transmission_t2

    with_loc_load::Bool = false

    isa_slack::Bool = false
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[Trans.algebraic_vars_syms...;]
    state_vars_syms::Vector{Symbol} = Symbol[Trans.state_vars_syms...; ]
    syms::Vector{Symbol} = Symbol[state_vars_syms...; algebraic_vars_syms...]
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]     
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64 = length(syms)    
    dims::Vector{Int64} = Int64[ Trans.dim ]
    dae_var::Vector{Bool} =  Bool[Trans.dae_var...;] 
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    param::Vector{Symbol} = Symbol[ Trans.param...;]
    param_values::Vector{ Float64} = Float64[ Trans.param_values...; ]
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [ param...; ] ) )
    func::Vector{Function} = Function[ initial_pf_plant_Transmission!, network_current_plant_Transmission!, global_pf_plant_Transmission!, node_pf_plant_Transmission!, hybrid_pf_plant_Transmission! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[  ]
    control_sig::Vector{Float64} = ones(length(control_sig_syms)) 
    output_sig_syms::Vector{Symbol} = Symbol[ Trans.output_sig_syms...; ]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))     
    cb_state_event_func::Vector{Function} = Function[ Trans.cb_state_event_func...; ]  
    cb_state_affect_func::Vector{Function} = Function[ Trans.cb_state_affect_func...;  ]    
    cb_state_syms::Vector{Symbol} = Symbol[ Trans.cb_state_syms...;  ]           
    cb_state_conditions::Vector{Float64} = Float64[ Trans.cb_state_conditions...; ]     
    cb_state_values::Vector{Function} = Function[ Trans.cb_state_values...; ]
    cb_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[ state_sym ] for state_sym in cb_state_syms ]
    cb_state_dim::Int64  = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[Trans.cb_dyn_state_event_func...; ]
    cb_dyn_state_affect_func::Vector{Function} = Function[Trans.cb_dyn_state_affect_func...; ]    
    cb_dyn_state_syms::Vector{Symbol} = Symbol[ Trans.cb_dyn_state_syms...; ]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ Trans.cb_dyn_state_sw_Idx...; ]
            
    cb_dyn_state_conditions::Vector{Float64} = Float64[ Trans.cb_dyn_state_conditions...; ]  
    cb_dyn_state_values::Vector{Function} = Function[ Trans.cb_dyn_state_values...; ]
    cb_dyn_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[state_sym] for state_sym in cb_dyn_state_syms]
    cb_dyn_param_state_sw::Vector{ Vector{Int64} } = Vector{Int64}[ Trans.cb_dyn_param_state_sw ]
    cb_dyn_state_dim::Int64 = length(cb_dyn_state_conditions) 
    
end

#--------------------------------------------------------
# Generation plants
#--------------------------------------------------------



function hybrid_pf_plant_cb!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p, global_pf_param, node_pf_param =  p_agg

    gen_ft, gov_f_t, avr_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[5]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]

    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
    end

            
    p_gen, p_gov, p_avr, plant_control_sig = p
    

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_gov_τm = u_gov[ dict_gov_syms[ gov_output_sigs[1] ] ]

    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, src_i, dst_i, gen_ft, p_gen, global_pf_param, node_pf_param ), t)
    
    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end


function node_pf_plant_cb!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p, node_pf_param =  p_agg

    gen_ft, gov_f_t, avr_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[4]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]

    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
    end

            
    p_gen, p_gov, p_avr, plant_control_sig = p
    

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_gov_τm = u_gov[ dict_gov_syms[ gov_output_sigs[1] ] ]

    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, src_i, dst_i, gen_ft, p_gen,  node_pf_param ), t)
    
    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end



function global_pf_plant_cb!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p, global_pf_param =  p_agg

    gen_ft, gov_f_t, avr_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[3]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]

    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
    end

            
    p_gen, p_gov, p_avr, plant_control_sig = p
    

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_gov_τm = u_gov[ dict_gov_syms[ gov_output_sigs[1] ] ]

    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, src_i, dst_i, gen_ft, p_gen,  global_pf_param ), t)
    
    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end


function network_current_plant_cb!(
    dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    gen_ft, gov_f_t, avr_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[2]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]


    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
    end

            
    p_gen, p_gov, p_avr, plant_control_sig = p
    

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_gov_τm = u_gov[ dict_gov_syms[ gov_output_sigs[1] ] ]

    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, src_i, dst_i, gen_ft, p_gen), t)
    
    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end


function initial_pf_plant_cb!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    gen_ft, gov_f_t, avr_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[1]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]


    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
    end

        
    p_gen, p_gov, p_avr, plant_control_sig = p


    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_gov_τm = u_gov[ dict_gov_syms[ gov_output_sigs[1] ] ]

    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ] 

    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, src_i, dst_i, gen_ft, p_gen), t)
    
    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end



function plant_cb_flattened!( dx, x, (p_agg, plant), t )
        
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    gen_ft, gov_f_t, avr_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[1]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]
    
    states_dims  = [gen.dim, gov.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx    = create_idxs(state_offset, states_dims)
    
    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]

    # dim_gen_cb_sw = gen.cb_dyn_state_dim
    # dim_gov_cb_sw = gov.cb_dyn_state_dim
    # dim_avr_cb_sw = avr.cb_dyn_state_dim

    # dims_cb_sw = [dim_gen_cb_sw, dim_gov_cb_sw, dim_avr_cb_sw]

    # cb_sw_offset = create_offsets(dims_cb_sw)
    
    # cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

    # cb_sw_gen = @view cb_sw[cb_sw_Idx[1]]
    # cb_sw_gov = @view cb_sw[cb_sw_Idx[2]]
    # cb_sw_avr = @view cb_sw[cb_sw_Idx[3]]

    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
    end


    
    # cb_sw_gen, cb_sw_gov, cb_sw_avr = cb_sw

    dims_params   = length.([gen.param, gov.param, avr.param])
    params_offset = create_offsets(dims_params)
    params_Idx    = create_idxs(params_offset, dims_params)
        
    p_gen = @view p[params_Idx[1]]
    p_gov = @view p[params_Idx[2]]
    p_avr = @view p[params_Idx[3]] 

    u_gen_ω     = u_gen[dict_gen_syms[:ω]]

    # u_gen_ur_ui = [u_gen[dict_gen_syms[:u_r]],
    #                u_gen[dict_gen_syms[:u_i]]
    #                ]
    
    # u_gov_τm = u_gov[dict_gov_syms[:τm_tilade]]
    # u_avr_vf = u_avr[dict_avr_syms[:vf_tilade]]
    
    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
     
    u_gov_τm = u_gov[ dict_gov_syms[ gov_output_sigs[1] ] ]
    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, src_i, dst_i, gen_ft, p_gen), t)
    
    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end

@kwdef struct plant_cb_v6 <: SdGenPlant
    Gen::SM_2axis_cb_v6
     Gov::Union{gov_ieee_tgov1_cb, gov_t0_cb, gov_t1_cb, gov_t1_cb_sauer}
    Exc::Union{avr_t0_cb, avr_t1_cb, avr_t1_cb_sauer}
    p_order::Float64 = 1.0
    v_ref::Float64 = 1.0
    ω_ref::Float64 = 1.0

    comp_type::Symbol = :plant_cb_v6               

    Bus::String      = Gen.Bus
    name::String     = Gen.name

    Bus_type::Symbol = Gen.Bus_type
    Bus_num::Int64   = Gen.Bus_num

    with_loc_load::Bool = false

    isa_slack::Bool = false

    isa_condenser::Bool = false
    
    plant_control_sig_syms::Vector{Symbol} = Symbol[ :p_order, :v_ref, :ω_ref]
    plant_control_sig::Vector{Float64} = Float64[ p_order, v_ref, ω_ref]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ Gen.algebraic_vars_syms...; Gov.algebraic_vars_syms...; Exc.algebraic_vars_syms... ]
    state_vars_syms::Vector{Symbol} = Symbol[ Gen.state_vars_syms...; Gov.state_vars_syms...; Exc.state_vars_syms...]
    syms::Vector{Symbol}  = Symbol[ Gen.syms...; Gov.syms...; Exc.syms... ]
        
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[Gen.im_algebraic_vars_syms...; Gov.im_algebraic_vars_syms...; Exc.im_algebraic_vars_syms... ]

    # im_algebraic_vars_syms::Vector{Symbol} = Symbol[ Gov.algebraic_vars_syms...; Exc.algebraic_vars_syms... ]

    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )    
    
    stab_state_vars_syms::Vector{Symbol} = Symbol[ Gen.stab_state_vars_syms...; Gov.state_vars_syms...; Exc.state_vars_syms... ]
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]        
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64       = length(syms)    
    dims::Vector{Int64} = Int64[ Gen.dim...; Gov.dim...; Exc.dim... ]
    dae_var::Vector{Bool} = Bool[ Gen.dae_var...; Gov.dae_var...; Exc.dae_var...]
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    param::Vector{Vector{Symbol}} = Vector{Symbol}[ Gen.param, Gov.param,  Exc.param, plant_control_sig_syms]
    param_values::Vector{Vector{Float64} } = Vector{Float64}[Gen.param_values, Gov.param_values, Exc.param_values, plant_control_sig]

    
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [param...;] ) )    
    func::Vector{Function} = Function[ initial_pf_plant_cb!, network_current_plant_cb!, global_pf_plant_cb!, node_pf_plant_cb!,  hybrid_pf_plant_cb! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[ plant_control_sig_syms...; Gen.control_sig_syms...; Gov.control_sig_syms...; Exc.control_sig_syms...]
    control_sig::Vector{Float64} = [
        plant_control_sig...;ones(length(control_sig_syms) - length(plant_control_sig_syms))...]
    output_sig_syms::Vector{Symbol} = Symbol[ Gen.output_sig_syms...; Gov.output_sig_syms...; Exc.output_sig_syms...]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))         
    cb_state_event_func::Vector{Function}      = Function[ Gen.cb_state_event_func...; Gov.cb_state_event_func...; Exc.cb_state_event_func... ]  
    cb_state_affect_func::Vector{Function}     = Function[ Gen.cb_state_affect_func...; Gov.cb_state_affect_func...; Exc.cb_state_affect_func... ]    

    #
    cb_state_syms::Vector{ Union{Symbol, Tuple{Symbol,Symbol}} } = Union{Symbol, Tuple{Symbol,Symbol}}[ Gen.cb_state_syms...; Gov.cb_state_syms...; Exc.cb_state_syms... ]
    #    

    cb_state_conditions::Vector{Float64}       = Float64[ Gen.cb_state_conditions...; Gov.cb_state_conditions...;  Exc.cb_state_conditions... ]     
    cb_state_values::Vector{Function}           = Function[ Gen.cb_state_values...; Gov.cb_state_values...; Exc.cb_state_values... ]

    #
    cb_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_state_syms[state_sym] : ( dict_state_syms[state_sym[1]], dict_state_syms[state_sym[1]])  for state_sym in  cb_state_syms ]
    #

    im_cb_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_im_vars_syms[ state_sym ] : ( dict_im_vars_syms[ state_sym[1] ], dict_im_vars_syms[ state_sym[1] ])  for state_sym in  cb_state_syms ]
    
    cb_state_dim::Int64 = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[ Gen.cb_dyn_state_event_func...; Gov.cb_dyn_state_event_func...; Exc.cb_dyn_state_event_func... ]
    cb_dyn_state_affect_func::Vector{Function} = Function[ Gen.cb_dyn_state_affect_func...; Gov.cb_dyn_state_affect_func...; Exc.cb_dyn_state_affect_func... ]    

    #
    cb_dyn_state_syms::Vector{ Union{Symbol, Tuple{Symbol,Symbol}} }  = Union{Symbol, Tuple{Symbol,Symbol}}[ Gen.cb_dyn_state_syms...; Gov.cb_dyn_state_syms...; Exc.cb_dyn_state_syms... ]
    #
    
    #
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{ Union{Symbol,Tuple{Symbol, Symbol}}, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    #
    
    #
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[  cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]  
    #    
    cb_dyn_state_conditions::Vector{Float64}   = Float64[ Gen.cb_dyn_state_conditions...; Gov.cb_dyn_state_conditions...; Exc.cb_dyn_state_conditions... ]  
    cb_dyn_state_values::Vector{Function}       = Function[ Gen.cb_dyn_state_values...; Gov.cb_dyn_state_values...; Exc.cb_dyn_state_values... ]

    #
    cb_dyn_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} } = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_state_syms[state_sym] : ( dict_state_syms[state_sym[1]], dict_state_syms[state_sym[1]])  for state_sym in  cb_dyn_state_syms ]
    #    
    
    cb_dyn_param_state_sw::Vector{Int64}       = Int64[ Gen.cb_dyn_param_state_sw...; Gov.cb_dyn_param_state_sw...; Exc.cb_dyn_param_state_sw... ]
    cb_dyn_state_dim::Int64 = length(cb_dyn_param_state_sw )


    dummy_vr1::Float64       = Exc.Ka * (v_ref - Gen.vh)
    dummy_vf_tilade::Float64 = dummy_vr1 * 0.97918


    func_system_matrices::Vector{Function} = Function[ get_a_plant_system_matrices!, get_a_plant_system_matrices ]

    func_update_system_matrices::Function = update_a_plant_system_matrices!

    # system_matrices = get_gen_gov_avr_plant_system_matrices(Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade )

    # stability

    func_system_stab_matrices::Vector{Function} = Function[ get_a_plant_system_stability_matrices!, get_a_plant_system_stability_matrices ]

    func_update_system_stab_matrices::Function = update_a_plant_system_stability_matrices!

    # system_stab_matrices = get_gen_gov_avr_plant_system_stability_matrices(Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade )
    
end



#---------------------------------------------------------------



function hybrid_pf_plant_no_gov!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p, global_pf_param, node_pf_param =  p_agg

    gen_ft, avr_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[5]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]

    if gen.cb_dyn_state_dim == 0 
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen  = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr  = @view cb_sw[cb_sw_Idx[2]]
    end
    
       
    p_gen, p_avr, plant_control_sig = p

    
    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    u_gov_τm = 0.0
    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_avr_vf, src_i, dst_i, gen_ft, p_gen, global_pf_param,  node_pf_param ), t)
    
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end


function node_pf_plant_no_gov!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p, node_pf_param =  p_agg

    gen_ft, avr_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[4]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]

    if gen.cb_dyn_state_dim == 0 
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen  = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr  = @view cb_sw[cb_sw_Idx[2]]
    end
    
       
    p_gen, p_avr, plant_control_sig = p

    
    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    u_gov_τm = 0.0
    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_avr_vf, src_i, dst_i, gen_ft, p_gen,  node_pf_param ), t)
    
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end


function global_pf_plant_no_gov!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p,  global_pf_param =  p_agg

    gen_ft, avr_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[3]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]

    if gen.cb_dyn_state_dim == 0 
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen  = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr  = @view cb_sw[cb_sw_Idx[2]]
    end
    
       
    p_gen, p_avr, plant_control_sig = p

    
    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    u_gov_τm = 0.0
    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen,  u_avr_vf, src_i, dst_i, gen_ft, p_gen,  global_pf_param ), t)
    
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end


function network_current_plant_no_gov!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    gen_ft, avr_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[2]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]


    if gen.cb_dyn_state_dim == 0 
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen  = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr  = @view cb_sw[cb_sw_Idx[2]]
    end
    
       
    p_gen, p_avr, plant_control_sig = p

    
    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    u_gov_τm = 0.0
    
    # gen_fun!(du_gen, u_gen,
    #          (cb_sw_gen, u_gov_τm, u_avr_vf, src_i, dst_i, gen_ft, p_gen), t)
    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_avr_vf, src_i, dst_i, gen_ft, p_gen), t)
    
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end



function initial_pf_plant_no_gov!(
    dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    gen_ft, avr_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[1]
    avr_fun! = avr.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   
    
    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]

    if gen.cb_dyn_state_dim == 0 
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen  = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr  = @view cb_sw[cb_sw_Idx[2]]
    end
    
       
    p_gen, p_avr, plant_control_sig = p

    
    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[ dict_gen_syms[ gen_output_sigs[ 1 ] ] ],
                   u_gen[ dict_gen_syms[ gen_output_sigs[ 2 ] ] ]
                   ]
  
    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    u_gov_τm = 0.0
    
    # gen_fun!(du_gen, u_gen,
    #          (cb_sw_gen, u_gov_τm, u_avr_vf, src_i, dst_i, gen_ft, p_gen), t)

gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_avr_vf, src_i, dst_i, gen_ft, p_gen), t)    
    
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end



function plant_no_gov_flattened!( dx, x, (p_agg, plant), t )
        
    cb_sw, src_i, dst_i, f_t, p =  p_agg

    gen_ft, avr_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    
    gen_fun! = gen.func[1]
    avr_fun! = avr.func[1]
    
    states_dims  = [gen.dim, avr.dim]
    state_offset = create_offsets(states_dims)
    state_Idx    = create_idxs(state_offset, states_dims)
    
    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]

    # dim_gen_cb_sw = gen.cb_dyn_state_dim
    # dim_avr_cb_sw = avr.cb_dyn_state_dim

    # dims_cb_sw = [dim_gen_cb_sw, dim_avr_cb_sw]

    # cb_sw_offset = create_offsets(dims_cb_sw)
    
    # cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

    # cb_sw_gen = @view cb_sw[cb_sw_Idx[1]]
    # cb_sw_avr = @view cb_sw[cb_sw_Idx[2]]

    if gen.cb_dyn_state_dim == 0 
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
    end
        
    dims_params   = length.([gen.param, avr.param])
    params_offset = create_offsets(dims_params)
    params_Idx    = create_idxs(params_offset, dims_params)
        
    p_gen = @view p[params_Idx[1]]
    p_avr = @view p[params_Idx[2]] 

    u_gen_ω     = u_gen[dict_gen_syms[:ω]]

    # u_gen_ur_ui = [u_gen[dict_gen_syms[:u_r]],
    #                u_gen[dict_gen_syms[:u_i]]
    #                ]
    
    # u_gov_τm = u_gov[dict_gov_syms[:τm_tilade]]
    # u_avr_vf = u_avr[dict_avr_syms[:vf_tilade]]
    
    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                   u_gen[dict_gen_syms[gen_output_sigs[2]]]
                   ]
     
    u_avr_vf = u_avr[ dict_avr_syms[ avr_output_sigs[1] ] ]    

    u_gov_τm = 0.0
    
    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_avr_vf, src_i, dst_i, gen_ft, p_gen), t)
    
    
    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t)    
    
    return nothing
        
end


@kwdef struct plant_no_gov_v6  <: SdGenPlant
    Gen::SC_2axis_cb_v6
    Exc::avr_t1_cb
    v_ref::Float64 = 1.0
    ω_ref::Float64 = 1.0

    comp_type::Symbol = :plant_no_gov_v6                

    Bus::String      = Gen.Bus
    name::String     = Gen.name

    Bus_type::Symbol = Gen.Bus_type
    Bus_num::Int64   = Gen.Bus_num

    with_loc_load::Bool = false

    isa_slack::Bool = false

    isa_condenser::Bool = true
    
    plant_control_sig_syms::Vector{Symbol} = Symbol[ :v_ref, :ω_ref]
    plant_control_sig::Vector{Float64} = Float64[ v_ref, ω_ref]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ Gen.algebraic_vars_syms...; Exc.algebraic_vars_syms... ]

    state_vars_syms::Vector{Symbol} = Symbol[ Gen.state_vars_syms...; Exc.state_vars_syms...]

    syms::Vector{Symbol} = Symbol[ Gen.syms...; Exc.syms... ]
    
    stab_state_vars_syms::Vector{Symbol} = Symbol[ Gen.stab_state_vars_syms...; Exc.state_vars_syms... ]
    
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[Gen.im_algebraic_vars_syms...;  Exc.im_algebraic_vars_syms... ]
    
    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]        
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64       = length(syms)    
    dims::Vector{Int64} = Int64[ Gen.dim...; Exc.dim... ]
    dae_var::Vector{Bool} = Bool[ Gen.dae_var...; Exc.dae_var...]
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    param::Vector{Vector{Symbol}} = Vector{Symbol}[ Gen.param, Exc.param, plant_control_sig_syms]
    param_values::Vector{Vector{Float64} } = Vector{Float64}[Gen.param_values, Exc.param_values, plant_control_sig]

    
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [param...;] ) )
    func::Vector{Function} = Function[ initial_pf_plant_no_gov!, network_current_plant_no_gov!,  global_pf_plant_no_gov!, node_pf_plant_no_gov!, hybrid_pf_plant_no_gov! ]    
    control_sig_syms::Vector{Symbol} = Symbol[ plant_control_sig_syms...; Exc.control_sig_syms...]
    control_sig::Vector{Float64} = [
        plant_control_sig...;ones(length(control_sig_syms) - length(plant_control_sig_syms))...]
    output_sig_syms::Vector{Symbol} = Symbol[ Gen.output_sig_syms...; Exc.output_sig_syms...]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))         
    cb_state_event_func::Vector{Function} = Function[ Gen.cb_state_event_func...; Exc.cb_state_event_func... ]  
    cb_state_affect_func::Vector{Function}     = Function[ Gen.cb_state_affect_func...; Exc.cb_state_affect_func... ]    
    cb_state_syms::Vector{Symbol} = Symbol[ Gen.cb_state_syms...; Exc.cb_state_syms... ]           
    cb_state_conditions::Vector{Float64} = Float64[ Gen.cb_state_conditions...;  Exc.cb_state_conditions... ]     
    cb_state_values::Vector{Function} = Function[ Gen.cb_state_values...; Exc.cb_state_values... ]
    cb_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[ state_sym ] for state_sym in cb_state_syms ]

    im_cb_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_im_vars_syms[ state_sym ] : ( dict_im_vars_syms[ state_sym[1] ], dict_im_vars_syms[ state_sym[1] ])  for state_sym in  cb_state_syms ]
    
    cb_state_dim::Int64 = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[ Gen.cb_dyn_state_event_func...; Exc.cb_dyn_state_event_func... ]
    cb_dyn_state_affect_func::Vector{Function} = Function[ Gen.cb_dyn_state_affect_func...; Exc.cb_dyn_state_affect_func... ]    
    cb_dyn_state_syms::Vector{Symbol} = Symbol[ Gen.cb_dyn_state_syms...; Exc.cb_dyn_state_syms... ]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))

    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[Gen.cb_dyn_state_sw_Idx...; Exc.cb_dyn_state_sw_Idx... ]
                
    cb_dyn_state_conditions::Vector{Float64} = Float64[ Gen.cb_dyn_state_conditions...; Exc.cb_dyn_state_conditions... ]
    
    cb_dyn_state_values::Vector{Function} = Function[ Gen.cb_dyn_state_values...; Exc.cb_dyn_state_values... ]
    
    cb_dyn_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[state_sym] for state_sym in cb_dyn_state_syms]
    
    cb_dyn_param_state_sw::Vector{Int64} = Int64[ Gen.cb_dyn_param_state_sw...; Exc.cb_dyn_param_state_sw... ]
    
    cb_dyn_state_dim::Int64 = length(cb_dyn_param_state_sw )

    # stability

    dummy_vr1::Float64       = Exc.Ka * (v_ref - Gen.vh)
    dummy_vf_tilade::Float64 = dummy_vr1 * 0.97918
    
    func_system_matrices::Vector{Function} = Function[get_a_plant_system_matrices!, get_a_plant_system_matrices]
    
    func_update_system_matrices::Function = update_a_plant_system_matrices!
    
    # system_matrices = get_only_gen_avr_plant_system_matrices( Gen, Exc, dummy_vr1, dummy_vf_tilade  )   

    func_system_stab_matrices::Vector{Function} = Function[get_a_plant_system_stability_matrices!, get_a_plant_system_stability_matrices]
    
    func_update_system_stab_matrices::Function = update_a_plant_system_stability_matrices!
    
    # system_stab_matrices = get_only_gen_avr_plant_system_stability_matrices( Gen, Exc, dummy_vr1, dummy_vf_tilade  )
    
end



@kwdef struct plant_no_gov <: SdGenPlant
    Gen::Union{SC_2axis_cb_idq, SC_2axis_cb_v6}
    Exc::Union{avr_t0_cb, avr_t1_cb, avr_t1_cb_sauer}
    v_ref::Float64 = 1.0
    ω_ref::Float64 = 1.0

    comp_type::Symbol = :plant_no_gov                    

    Bus::String      = Gen.Bus
    name::String     = Gen.name

    Bus_type::Symbol = Gen.Bus_type
    Bus_num::Int64   = Gen.Bus_num

    with_loc_load::Bool = false

    isa_slack::Bool = false

    isa_condenser::Bool = true
    
    plant_control_sig_syms::Vector{Symbol} = Symbol[ :v_ref, :ω_ref]
    plant_control_sig::Vector{Float64} = Float64[ v_ref, ω_ref]
    algebraic_vars_syms::Vector{Symbol} = Symbol[ Gen.algebraic_vars_syms...; Exc.algebraic_vars_syms... ]

    state_vars_syms::Vector{Symbol} = Symbol[ Gen.state_vars_syms...; Exc.state_vars_syms...]

    syms::Vector{Symbol}  = Symbol[ Gen.syms...; Exc.syms... ]

    stab_state_vars_syms::Vector{Symbol} = Symbol[ Gen.stab_state_vars_syms...;  Exc.state_vars_syms... ]
    
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[Gen.im_algebraic_vars_syms...;  Exc.im_algebraic_vars_syms... ]

    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]        
    state_dim::Int64 = length(state_vars_syms) 
    dim::Int64       = length(syms)    
    dims::Vector{Int64} = Int64[ Gen.dim...; Exc.dim... ]
    dae_var::Vector{Bool} = Bool[ Gen.dae_var...; Exc.dae_var...]
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    param::Vector{Vector{Symbol}} = Vector{Symbol}[ Gen.param, Exc.param, plant_control_sig_syms]
    param_values::Vector{Vector{Float64} } = Vector{Float64}[Gen.param_values, Exc.param_values, plant_control_sig]

    
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [param...;] ) )
    func::Vector{Function} = Function[ initial_pf_plant_no_gov!, network_current_plant_no_gov!,  global_pf_plant_no_gov!, node_pf_plant_no_gov!, hybrid_pf_plant_no_gov! ]    
    control_sig_syms::Vector{Symbol} = Symbol[ plant_control_sig_syms...; Exc.control_sig_syms...]
    control_sig::Vector{Float64} = [
        plant_control_sig...;ones(length(control_sig_syms) - length(plant_control_sig_syms))...]
    output_sig_syms::Vector{Symbol} = Symbol[ Gen.output_sig_syms...; Exc.output_sig_syms...]
    output_sig::Vector{Float64} = ones(length(output_sig_syms))         
    cb_state_event_func::Vector{Function} = Function[ Gen.cb_state_event_func...; Exc.cb_state_event_func... ]  
    cb_state_affect_func::Vector{Function}     = Function[ Gen.cb_state_affect_func...; Exc.cb_state_affect_func... ]    
    cb_state_syms::Vector{Symbol} = Symbol[ Gen.cb_state_syms...; Exc.cb_state_syms... ]           
    cb_state_conditions::Vector{Float64} = Float64[ Gen.cb_state_conditions...;  Exc.cb_state_conditions... ]     
    cb_state_values::Vector{Function} = Function[ Gen.cb_state_values...; Exc.cb_state_values... ]

    cb_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[ state_sym ] for state_sym in cb_state_syms ]
    cb_state_dim::Int64 = length(cb_state_conditions)

    im_cb_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_im_vars_syms[ state_sym ] : ( dict_im_vars_syms[ state_sym[1] ], dict_im_vars_syms[ state_sym[1] ])  for state_sym in  cb_state_syms ]
    
    cb_dyn_state_event_func::Vector{Function}  = Function[ Gen.cb_dyn_state_event_func...; Exc.cb_dyn_state_event_func... ]
    cb_dyn_state_affect_func::Vector{Function} = Function[ Gen.cb_dyn_state_affect_func...; Exc.cb_dyn_state_affect_func... ]    
    cb_dyn_state_syms::Vector{Symbol} = Symbol[ Gen.cb_dyn_state_syms...; Exc.cb_dyn_state_syms... ]
    
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))

    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[Gen.cb_dyn_state_sw_Idx...; Exc.cb_dyn_state_sw_Idx... ]
                
    cb_dyn_state_conditions::Vector{Float64} = Float64[ Gen.cb_dyn_state_conditions...; Exc.cb_dyn_state_conditions... ]
    
    cb_dyn_state_values::Vector{Function} = Function[ Gen.cb_dyn_state_values...; Exc.cb_dyn_state_values... ]
    
    cb_dyn_state_sym2Idx::Vector{Int64} = Int64[dict_state_syms[state_sym] for state_sym in cb_dyn_state_syms]
    
    cb_dyn_param_state_sw::Vector{Int64} = Int64[ Gen.cb_dyn_param_state_sw...; Exc.cb_dyn_param_state_sw... ]
    
    cb_dyn_state_dim::Int64 = length(cb_dyn_param_state_sw )


    func_system_matrices::Vector{Function} = Function[get_a_plant_system_matrices!, get_a_plant_system_matrices]
    
    func_update_system_matrices::Function = update_a_plant_system_matrices!
    
    # system_matrices = get_only_gen_avr_plant_system_matrices( Gen, Exc, dummy_vr1, dummy_vf_tilade  )

    # stability

    dummy_vr1::Float64       = Exc.Ka * (v_ref - Gen.vh)
    dummy_vf_tilade::Float64 = dummy_vr1 * 0.97918    

    func_system_stab_matrices::Vector{Function} = Function[get_a_plant_system_stability_matrices!, get_a_plant_system_stability_matrices]
    
    func_update_system_stab_matrices::Function = update_a_plant_system_stability_matrices!
    
    # system_stab_matrices = get_only_gen_avr_plant_system_stability_matrices( Gen, Exc, dummy_vr1, dummy_vf_tilade  )
    
end



#--------------------------------------------------------




function hybrid_pf_plant_wt_loc_load!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p, global_pf_param, node_pf_param =  p_agg
  
    gen_ft, gov_f_t, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[5]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim,
                   avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    u_loc_load  = @view x[state_Idx[4]]

    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]
    du_loc_load = @view dx[state_Idx[4]]

    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        cb_sw_loc_load = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0  && loc_load.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[4]]        
    end

    # cb_sw_gen, cb_sw_gov, cb_sw_avr, cb_sw_loc_load = cb_sw
    
    p_gen, p_gov, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_gov_τm = u_gov[dict_gov_syms[gov_output_sigs[1]]]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]]
    ]    

    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, u_loc_load_ir_ii, src_i, dst_i, gen_ft, p_gen,  global_pf_param, node_pf_param ), t)

    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load, u_gen_ur_ui, loc_load_f_t, p_loc_load), t)
        
    return nothing
        
end




function node_pf_plant_wt_loc_load!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p, node_pf_param =  p_agg
  
    gen_ft, gov_f_t, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[4]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim,
                   avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    u_loc_load  = @view x[state_Idx[4]]

    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]
    du_loc_load = @view dx[state_Idx[4]]

    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        cb_sw_loc_load = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0  && loc_load.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[4]]        
    end

    # cb_sw_gen, cb_sw_gov, cb_sw_avr, cb_sw_loc_load = cb_sw
    
    p_gen, p_gov, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_gov_τm = u_gov[dict_gov_syms[gov_output_sigs[1]]]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]]
    ]    

    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, u_loc_load_ir_ii, src_i, dst_i, gen_ft, p_gen, node_pf_param ), t)

    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load, u_gen_ur_ui, loc_load_f_t, p_loc_load), t)
        
    return nothing
        
end


function global_pf_plant_wt_loc_load!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p, global_pf_param =  p_agg
  
    gen_ft, gov_f_t, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[3]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim,
                   avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    u_loc_load  = @view x[state_Idx[4]]

    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]
    du_loc_load = @view dx[state_Idx[4]]

    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        cb_sw_loc_load = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0  && loc_load.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[4]]        
    end

    # cb_sw_gen, cb_sw_gov, cb_sw_avr, cb_sw_loc_load = cb_sw
    
    p_gen, p_gov, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_gov_τm = u_gov[dict_gov_syms[gov_output_sigs[1]]]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]]
    ]    

    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, u_loc_load_ir_ii, src_i, dst_i, gen_ft, p_gen, global_pf_param ), t)

    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load, u_gen_ur_ui, loc_load_f_t, p_loc_load), t)
        
    return nothing
        
end


function network_current_plant_wt_loc_load!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg
  
    gen_ft, gov_f_t, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[2]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim,
                   avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    u_loc_load  = @view x[state_Idx[4]]

    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]
    du_loc_load = @view dx[state_Idx[4]]

    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        cb_sw_loc_load = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0  && loc_load.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[4]]        
    end

    # cb_sw_gen, cb_sw_gov, cb_sw_avr, cb_sw_loc_load = cb_sw
    
    p_gen, p_gov, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_gov_τm = u_gov[dict_gov_syms[gov_output_sigs[1]]]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]]
    ]    

    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, u_loc_load_ir_ii, src_i, dst_i, gen_ft, p_gen), t)

    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load, u_gen_ur_ui, loc_load_f_t, p_loc_load), t)
        
    return nothing
        
end

function initial_pf_plant_wt_loc_load!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg
  
    gen_ft, gov_f_t, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    gov = plant.Gov
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    gov_states_syms = gov.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_gov_syms = gov.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    gov_output_sigs = gov.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[1]
    gov_fun! = gov.func[1]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    # states_dims    = plant.dims

    states_dims = [gen.dim, gov.dim,
                   avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_gov  = @view x[state_Idx[2]]    
    u_avr  = @view x[state_Idx[3]]
    u_loc_load  = @view x[state_Idx[4]]

    du_gen = @view dx[state_Idx[1]] 
    du_gov = @view dx[state_Idx[2]]
    du_avr = @view dx[state_Idx[3]]
    du_loc_load = @view dx[state_Idx[4]]

    if gen.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0

        cb_sw_gen      = Int64[]
        cb_sw_gov      = Int64[]
        cb_sw_loc_load = Int64[]
        
        cb_sw_avr      = cb_sw
        

    elseif gen.cb_dyn_state_dim == 0  && loc_load.cb_dyn_state_dim == 0 && gov.cb_dyn_state_dim != 0

        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]
        
        dims_cb_sw = [gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gov      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
    else
        
        dims_cb_sw = [gen.cb_dyn_state_dim,
                      gov.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_gov      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[3]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[4]]        
    end

    # cb_sw_gen, cb_sw_gov, cb_sw_avr, cb_sw_loc_load = cb_sw
    
    p_gen, p_gov, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_gov_τm = u_gov[dict_gov_syms[gov_output_sigs[1]]]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]]
    ]    

    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_gov_τm, u_avr_vf, u_loc_load_ir_ii, src_i, dst_i, gen_ft, p_gen), t)

    gov_fun!(du_gov, u_gov,
             (cb_sw_gov, u_gen_ω, gov_f_t, p_gov), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr, u_gen_ur_ui, avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load, u_gen_ur_ui, loc_load_f_t, p_loc_load), t)
    
    return nothing
        
end


@kwdef struct plant_wt_loc_load_v6 <: SdGenPlant
    Gen::SM_2axis_wt_loc_load_cb_v6
    Gov::Union{gov_ieee_tgov1_cb, gov_t0_cb, gov_t1_cb, gov_t1_cb_sauer }
    Exc::Union{avr_t1_cb, avr_t0_cb, avr_t1_cb_sauer }
    # Gov::gov_ieee_tgov1_cb    
    # Exc::avr_t1_cb
    Loc_load::loc_Load_t1
    p_order::Float64 = 1.0
    v_ref::Float64 = 1.0
    ω_ref::Float64 = 1.0

    comp_type::Symbol = :plant_wt_loc_load_v6           

    Bus::String      = Gen.Bus
    name::String     = Gen.name

    Bus_type::Symbol = Gen.Bus_type
    Bus_num::Int64   = Gen.Bus_num

    with_loc_load::Bool = true

    isa_slack::Bool = false

    isa_condenser::Bool = false
    
    plant_control_sig_syms::Vector{Symbol} = Symbol[ :p_order, :v_ref, :ω_ref]
    
    plant_control_sig::Vector{Float64} = Float64[ p_order, v_ref, ω_ref]
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[ Gen.algebraic_vars_syms...; Gov.algebraic_vars_syms...; Exc.algebraic_vars_syms...; Loc_load.algebraic_vars_syms...]
    
    state_vars_syms::Vector{Symbol} = Symbol[ Gen.state_vars_syms...; Gov.state_vars_syms...; Exc.state_vars_syms...; Loc_load.state_vars_syms... ]
    
    syms::Vector{Symbol}  = Symbol[ Gen.syms...; Gov.syms...; Exc.syms...; Loc_load.syms... ]

    stab_state_vars_syms::Vector{Symbol} = Symbol[ Gen.stab_state_vars_syms...; Gov.state_vars_syms...; Exc.state_vars_syms...; Loc_load.state_vars_syms... ]
    
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[Gen.im_algebraic_vars_syms...; Gov.im_algebraic_vars_syms...; Exc.im_algebraic_vars_syms... ]
    
    # im_algebraic_vars_syms::Vector{Symbol} = Symbol[  Gov.algebraic_vars_syms...; Exc.algebraic_vars_syms... ]

    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )
    
    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
        
    
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]
    
    state_dim::Int64 = length(state_vars_syms)
    
    dim::Int64 = length(syms)
    
    dims::Vector{Int64} = Int64[ Gen.dim...; Gov.dim...; Exc.dim...; Loc_load.dim ]
    
    dae_var::Vector{Bool} = Bool[ Gen.dae_var...; Gov.dae_var...; Exc.dae_var...; Loc_load.dae_var... ]
    
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    
    param::Vector{Vector{Symbol}} = Vector{Symbol}[ Gen.param, Gov.param,  Exc.param, Loc_load.param, plant_control_sig_syms ]
    param_values::Vector{Vector{Float64} } = Vector{Float64}[Gen.param_values, Gov.param_values, Exc.param_values, Loc_load.param_values, plant_control_sig]

    loc_load_param_idx::Int64 = 4

    
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [param...;] ) )
    func::Vector{Function} = Function[ initial_pf_plant_wt_loc_load!, network_current_plant_wt_loc_load!, global_pf_plant_wt_loc_load!, node_pf_plant_wt_loc_load!, hybrid_pf_plant_wt_loc_load! ]
    
    control_sig_syms::Vector{Symbol} = Symbol[ plant_control_sig_syms...; Gen.control_sig_syms...; Gov.control_sig_syms...; Exc.control_sig_syms...; Loc_load.control_sig_syms... ]
    
    control_sig::Vector{Float64} = [ plant_control_sig...; ones(length(control_sig_syms) - length(plant_control_sig_syms))...]
    
    output_sig_syms::Vector{Symbol} = Symbol[ Gen.output_sig_syms...; Gov.output_sig_syms...; Exc.output_sig_syms...; Loc_load.output_sig_syms... ]
    
    output_sig::Vector{Float64} = ones(length(output_sig_syms))         
    cb_state_event_func::Vector{Function} = Function[ Gen.cb_state_event_func...; Gov.cb_state_event_func...; Exc.cb_state_event_func...; Loc_load.cb_state_event_func... ]
    
    cb_state_affect_func::Vector{Function} = Function[ Gen.cb_state_affect_func...; Gov.cb_state_affect_func...; Exc.cb_state_affect_func...; Loc_load.cb_state_affect_func... ]
    
    #
    cb_state_syms::Vector{ Union{Symbol, Tuple{Symbol,Symbol}} } = Union{Symbol, Tuple{Symbol,Symbol}}[ Gen.cb_state_syms...; Gov.cb_state_syms...; Exc.cb_state_syms...; Loc_load.cb_state_syms... ]
    #
    
    cb_state_conditions::Vector{Float64} = Float64[ Gen.cb_state_conditions...; Gov.cb_state_conditions...; Exc.cb_state_conditions...; Loc_load.cb_state_conditions... ]
    
    cb_state_values::Vector{Function} = Function[ Gen.cb_state_values...; Gov.cb_state_values...; Exc.cb_state_values...; Loc_load.cb_state_values... ]
    
    #
    cb_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_state_syms[state_sym] : ( dict_state_syms[state_sym[1]], dict_state_syms[state_sym[1]])  for state_sym in  cb_state_syms ]

    im_cb_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_im_vars_syms[ state_sym ] : ( dict_im_vars_syms[ state_sym[1] ], dict_im_vars_syms[ state_sym[1] ])  for state_sym in  cb_state_syms ]

    #
    
    cb_state_dim::Int64  = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[ Gen.cb_dyn_state_event_func...; Gov.cb_dyn_state_event_func...; Exc.cb_dyn_state_event_func...; Loc_load.cb_dyn_state_event_func... ]

    cb_dyn_state_affect_func::Vector{Function} = Function[ Gen.cb_dyn_state_affect_func...; Gov.cb_dyn_state_affect_func...; Exc.cb_dyn_state_affect_func...; Loc_load.cb_dyn_state_affect_func... ]

    #
    cb_dyn_state_syms::Vector{ Union{Symbol, Tuple{Symbol,Symbol}} }  = Union{Symbol, Tuple{Symbol,Symbol}}[ Gen.cb_dyn_state_syms...; Gov.cb_dyn_state_syms...; Exc.cb_dyn_state_syms...; Loc_load.cb_dyn_state_syms... ]
    #
    
    #
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Union{Symbol,Tuple{Symbol, Symbol}}, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    #
    
    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[  cb_dict_dyn_state_syms2sw_Idx[sym]  for sym in cb_dyn_state_syms ]    
               
    cb_dyn_state_conditions::Vector{Float64}   = Float64[ Gen.cb_dyn_state_conditions...; Gov.cb_dyn_state_conditions...; Exc.cb_dyn_state_conditions...; Loc_load.cb_dyn_state_conditions... ]
    
    cb_dyn_state_values::Vector{Function} = Function[ Gen.cb_dyn_state_values...; Gov.cb_dyn_state_values...; Exc.cb_dyn_state_values...; Loc_load.cb_dyn_state_values... ]

    #
    cb_dyn_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} } = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_state_syms[state_sym] : ( dict_state_syms[state_sym[1]], dict_state_syms[state_sym[1]])  for state_sym in  cb_dyn_state_syms ]
    #
    
    cb_dyn_param_state_sw::Vector{Int64} = Int64[ Gen.cb_dyn_param_state_sw...; Gov.cb_dyn_param_state_sw...; Exc.cb_dyn_param_state_sw...; Loc_load.cb_dyn_param_state_sw... ]
    
    cb_dyn_state_dim::Int64 = length(cb_dyn_param_state_sw )

    # cb_dyn_param_state_sw::Vector{Vector{Int64}} = Vector{Int64}[ Gen.cb_dyn_param_state_sw, Gov.cb_dyn_param_state_sw, Exc.cb_dyn_param_state_sw, Loc_load.cb_dyn_param_state_sw ]

    # cb_dyn_state_dim::Int64 = length(get_non_null_list(cb_dyn_param_state_sw))


    dummy_vr1::Float64       = Exc.Ka * (v_ref - Gen.vh)
    
    dummy_vf_tilade::Float64 = dummy_vr1 * 0.97918


    func_system_matrices::Vector{Function} = Function[ get_a_plant_system_matrices!, get_a_plant_system_matrices ]

    func_update_system_matrices::Function = update_a_plant_system_matrices!

    # system_matrices = get_gen_gov_avr_plant_system_matrices(Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade )

    # stability

    func_system_stab_matrices::Vector{Function} = Function[ get_a_plant_system_stability_matrices!, get_a_plant_system_stability_matrices ]

    func_update_system_stab_matrices::Function = update_a_plant_system_stability_matrices!

    # system_stab_matrices = get_gen_gov_avr_plant_system_stability_matrices(Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade )
    
    
end

#----------------------------------------------------------


function hybrid_pf_plant_no_gov_wt_loc_load!( dx, x, (p_agg, plant ), t )
    
    cb_sw, src_i, dst_i, f_t, p, global_pf_param, node_pf_param =  p_agg
    
    gen_ft, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[5]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    states_dims = [gen.dim, avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    u_loc_load  = @view x[state_Idx[3]]

    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    du_loc_load = @view dx[state_Idx[3]]

    if gen.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]

    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[3]]
    end
    
    p_gen, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]] ]    

    gen_fun!( du_gen, u_gen,
             ( cb_sw_gen, u_avr_vf,
              u_loc_load_ir_ii, src_i, dst_i,
              gen_ft, p_gen, global_pf_param, node_pf_param ), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr,
              u_gen_ur_ui,
              avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load,
                   u_gen_ur_ui,
                   loc_load_f_t, p_loc_load), t)

    return nothing
        
end



function node_pf_plant_no_gov_wt_loc_load!( dx, x, (p_agg, plant ), t )
    
    cb_sw, src_i, dst_i, f_t, p, node_pf_param =  p_agg
    
    gen_ft, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[4]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    states_dims = [gen.dim, avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    u_loc_load  = @view x[state_Idx[3]]

    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    du_loc_load = @view dx[state_Idx[3]]

    if gen.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]

    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[3]]
    end
    
    p_gen, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]] ]    

    gen_fun!( du_gen, u_gen,
             ( cb_sw_gen, u_avr_vf,
              u_loc_load_ir_ii, src_i, dst_i,
              gen_ft, p_gen, node_pf_param ), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr,
              u_gen_ur_ui,
              avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load,
                   u_gen_ur_ui,
                   loc_load_f_t, p_loc_load), t)

    return nothing
        
end



function global_pf_plant_no_gov_wt_loc_load!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p, global_pf_param =  p_agg
    
    gen_ft, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[3]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    states_dims = [gen.dim, avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    u_loc_load  = @view x[state_Idx[3]]

    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    du_loc_load = @view dx[state_Idx[3]]

    if gen.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]

    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[3]]
    end
    
    p_gen, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]] ]    

    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_avr_vf,
              u_loc_load_ir_ii, src_i, dst_i,
              gen_ft, p_gen, global_pf_param ), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr,
              u_gen_ur_ui,
              avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load,
                   u_gen_ur_ui,
                   loc_load_f_t, p_loc_load), t)

    return nothing
        
end


function network_current_plant_no_gov_wt_loc_load!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg
    
    gen_ft, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[2]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    states_dims = [gen.dim, avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    u_loc_load  = @view x[state_Idx[3]]

    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    du_loc_load = @view dx[state_Idx[3]]

    if gen.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]

    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[3]]
    end
    
    p_gen, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]] ]    

    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_avr_vf,
              u_loc_load_ir_ii, src_i, dst_i,
              gen_ft, p_gen), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr,
              u_gen_ur_ui,
              avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load,
                   u_gen_ur_ui,
                   loc_load_f_t, p_loc_load), t)

    return nothing
        
end


function initial_pf_plant_no_gov_wt_loc_load!( dx, x, (p_agg, plant), t )
    
    cb_sw, src_i, dst_i, f_t, p =  p_agg
    
    gen_ft, avr_f_t, loc_load_f_t = f_t

    gen = plant.Gen
    avr = plant.Exc
    loc_load = plant.Loc_load

    gen_states_syms = gen.syms
    avr_states_syms = avr.syms
    loc_load_states_syms = loc_load.syms

    dict_gen_syms = gen.dict_state_syms
    dict_avr_syms = avr.dict_state_syms
    dict_loc_load_syms = loc_load.dict_state_syms

    gen_output_sigs = gen.output_sig_syms
    avr_output_sigs = avr.output_sig_syms
    loc_load_output_sigs = loc_load.output_sig_syms

    gen_fun! = gen.func[1]
    avr_fun! = avr.func[1]
    loc_load_fun! = loc_load.func[1]

    states_dims = [gen.dim, avr.dim, loc_load.dim]

    state_offset = create_offsets(states_dims)
    state_Idx   = create_idxs(state_offset, states_dims)   

    u_gen  = @view x[state_Idx[1]]    
    u_avr  = @view x[state_Idx[2]]
    u_loc_load  = @view x[state_Idx[3]]

    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    du_loc_load = @view dx[state_Idx[3]]

    if gen.cb_dyn_state_dim == 0 && loc_load.cb_dyn_state_dim == 0
        cb_sw_avr      = cb_sw
        cb_sw_gen      = Int64[]
        cb_sw_loc_load = Int64[]

    else

        dims_cb_sw = [gen.cb_dyn_state_dim,
                      avr.cb_dyn_state_dim,
                      loc_load.cb_dyn_state_dim]

        cb_sw_offset = create_offsets(dims_cb_sw)
        cb_sw_Idx = create_idxs(cb_sw_offset, dims_cb_sw)

        cb_sw_gen      = @view cb_sw[cb_sw_Idx[1]]
        cb_sw_avr      = @view cb_sw[cb_sw_Idx[2]]
        cb_sw_loc_load = @view cb_sw[cb_sw_Idx[3]]
    end
    
    p_gen, p_avr, p_loc_load, plant_control_sig = p

    u_gen_ω  = u_gen[dict_gen_syms[:ω]]

    u_gen_ur_ui = [u_gen[dict_gen_syms[gen_output_sigs[1]]],
                  u_gen[ dict_gen_syms[gen_output_sigs[2]]]
                   ]

    u_avr_vf = u_avr[dict_avr_syms[avr_output_sigs[1]]]    

    u_loc_load_ir_ii = [
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[1]]],
        u_loc_load[
            dict_loc_load_syms[
                loc_load_output_sigs[2]]] ]    

    gen_fun!(du_gen, u_gen,
             (cb_sw_gen, u_avr_vf,
              u_loc_load_ir_ii, src_i, dst_i,
              gen_ft, p_gen), t)

    avr_fun!(du_avr, u_avr,
             (cb_sw_avr,
              u_gen_ur_ui,
              avr_f_t, p_avr), t) 

    loc_load_fun!(du_loc_load, u_loc_load,
                  (cb_sw_loc_load,
                   u_gen_ur_ui,
                   loc_load_f_t, p_loc_load), t)

    return nothing
        
end


@kwdef struct plant_no_gov_wt_loc_load_v6 <: SdGenPlant
    Gen::SC_2axis_wt_loc_load_cb_v6
    Exc::Union{avr_t0_cb, avr_t1_cb, avr_t1_cb_sauer }
    Loc_load::loc_Load_t1
    v_ref::Float64 = 1.0
    ω_ref::Float64 = 1.0

    comp_type::Symbol = :plant_no_gov_wt_loc_load_v6    

    Bus::String  = Gen.Bus    
    name::String  = Gen.name

    Bus_type::Symbol = Gen.Bus_type
    Bus_num::Int64   = Gen.Bus_num

    with_loc_load::Bool = true

    isa_slack::Bool = false

    isa_condenser::Bool = true

    plant_control_sig_syms::Vector{Symbol} = Symbol[ :v_ref, :ω_ref]
    
    plant_control_sig::Vector{Float64} = Float64[ v_ref, ω_ref]
    
    algebraic_vars_syms::Vector{Symbol} = Symbol[ Gen.algebraic_vars_syms...; Exc.algebraic_vars_syms...; Loc_load.algebraic_vars_syms...]
    
    state_vars_syms::Vector{Symbol} = Symbol[ Gen.state_vars_syms...; Exc.state_vars_syms...; Loc_load.state_vars_syms... ]
    
    syms::Vector{Symbol}  = Symbol[ Gen.syms...; Exc.syms...; Loc_load.syms... ]
    
    stab_state_vars_syms::Vector{Symbol} = Symbol[ Gen.stab_state_vars_syms...;  Exc.state_vars_syms...; Loc_load.state_vars_syms... ]
    
    im_algebraic_vars_syms::Vector{Symbol} = Symbol[Gen.im_algebraic_vars_syms...;  Exc.im_algebraic_vars_syms... ]
    
    # im_algebraic_vars_syms::Vector{Symbol} = Symbol[ Exc.algebraic_vars_syms...; ]

    im_vars_syms::Vector{Symbol} = Symbol[ state_vars_syms...; im_algebraic_vars_syms... ]

    dict_im_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( im_vars_syms  ) )
    
    dict_state_vars_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(state_vars_syms))
    
    dict_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(syms))

    dict_stab_state_syms::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(stab_state_vars_syms))
            
    ur_ui_idx::Vector{Int64} = Int64[dict_state_syms[:u_r],dict_state_syms[:u_i]]
    
    state_dim::Int64 = length(state_vars_syms)
    
    dim::Int64 = length(syms)
    
    dims::Vector{Int64} = Int64[ Gen.dim...; Exc.dim...; Loc_load.dim ]
    
    dae_var::Vector{Bool} = Bool[ Gen.dae_var...; Exc.dae_var...; Loc_load.dae_var... ]
    
    mass_matrix::Diagonal{Int64, Vector{Int64}} = Diagonal(Int.(dae_var))
    
    param::Vector{Vector{Symbol}} = Vector{Symbol}[ Gen.param, Exc.param, Loc_load.param, plant_control_sig_syms ]
    param_values::Vector{Vector{Float64} } = Vector{Float64}[Gen.param_values, Exc.param_values, Loc_load.param_values, plant_control_sig]

    loc_load_param_idx::Int64 = 3

    
    
    dict_param_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( [param...;] ) )
    func::Vector{Function} = Function[ initial_pf_plant_no_gov_wt_loc_load!, network_current_plant_no_gov_wt_loc_load!, global_pf_plant_no_gov_wt_loc_load!, node_pf_plant_no_gov_wt_loc_load!, hybrid_pf_plant_no_gov_wt_loc_load!  ]
    
    control_sig_syms::Vector{Symbol} = Symbol[ plant_control_sig_syms...; Gen.control_sig_syms...; Exc.control_sig_syms...; Loc_load.control_sig_syms... ]
    
    control_sig::Vector{Float64} = [ plant_control_sig...; ones(length(control_sig_syms) - length(plant_control_sig_syms))...]
    
    output_sig_syms::Vector{Symbol} = Symbol[ Gen.output_sig_syms...; Exc.output_sig_syms...; Loc_load.output_sig_syms... ]
    
    output_sig::Vector{Float64} = ones(length(output_sig_syms))         
    cb_state_event_func::Vector{Function} = Function[ Gen.cb_state_event_func...; Exc.cb_state_event_func...; Loc_load.cb_state_event_func... ]
    
    cb_state_affect_func::Vector{Function} = Function[ Gen.cb_state_affect_func...; Exc.cb_state_affect_func...; Loc_load.cb_state_affect_func... ]
    
    #
    cb_state_syms::Vector{ Union{Symbol, Tuple{Symbol,Symbol}} } = Union{Symbol, Tuple{Symbol,Symbol}}[ Gen.cb_state_syms...; Exc.cb_state_syms...; Loc_load.cb_state_syms... ]
    #
    
    cb_state_conditions::Vector{Float64} = Float64[ Gen.cb_state_conditions...; Exc.cb_state_conditions...; Loc_load.cb_state_conditions... ]
    
    cb_state_values::Vector{Function} = Function[ Gen.cb_state_values...; Exc.cb_state_values...; Loc_load.cb_state_values... ]
    
    #
    cb_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_state_syms[state_sym] : ( dict_state_syms[state_sym[1]], dict_state_syms[state_sym[1]])  for state_sym in  cb_state_syms ]
    #

    im_cb_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} }  = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_im_vars_syms[ state_sym ] : ( dict_im_vars_syms[ state_sym[1] ], dict_im_vars_syms[ state_sym[1] ])  for state_sym in  cb_state_syms ]
    
    cb_state_dim::Int64  = length(cb_state_conditions)

    cb_dyn_state_event_func::Vector{Function}  = Function[ Gen.cb_dyn_state_event_func...; Exc.cb_dyn_state_event_func...; Loc_load.cb_dyn_state_event_func... ]
    cb_dyn_state_affect_func::Vector{Function} = Function[ Gen.cb_dyn_state_affect_func...; Exc.cb_dyn_state_affect_func...; Loc_load.cb_dyn_state_affect_func... ]    

    #
    cb_dyn_state_syms::Vector{ Union{Symbol, Tuple{Symbol,Symbol}} }  = Union{Symbol, Tuple{Symbol,Symbol}}[ Gen.cb_dyn_state_syms...; Exc.cb_dyn_state_syms...; Loc_load.cb_dyn_state_syms... ]
    #
    
    #
    cb_dict_dyn_state_syms2sw_Idx::OrderedDict{Union{Symbol,Tuple{Symbol, Symbol}}, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate(unique( cb_dyn_state_syms ) ))
    #

    cb_dyn_state_sw_Idx::Vector{Int64 } = Int64[ cb_dict_dyn_state_syms2sw_Idx[sym] for sym in cb_dyn_state_syms]
               
    cb_dyn_state_conditions::Vector{Float64}   = Float64[ Gen.cb_dyn_state_conditions...; Exc.cb_dyn_state_conditions...; Loc_load.cb_dyn_state_conditions... ]
    
    cb_dyn_state_values::Vector{Function} = Function[ Gen.cb_dyn_state_values...; Exc.cb_dyn_state_values...; Loc_load.cb_dyn_state_values... ]
    
    #
    cb_dyn_state_sym2Idx::Vector{ Union{Int64, Tuple{Int64, Int64}} } = Union{Int64, Tuple{Int64, Int64}}[Symbol(typeof(state_sym)) == :Symbol ?  dict_state_syms[state_sym] : ( dict_state_syms[state_sym[1]], dict_state_syms[state_sym[1]])  for state_sym in  cb_dyn_state_syms ]
    #    
    
    cb_dyn_param_state_sw::Vector{Int64} = Int64[ Gen.cb_dyn_param_state_sw...; Exc.cb_dyn_param_state_sw...; Loc_load.cb_dyn_param_state_sw... ]
    
    cb_dyn_state_dim::Int64 = length(cb_dyn_param_state_sw )

    # cb_dyn_param_state_sw::Vector{Vector{Int64}} = Vector{Int64}[ Gen.cb_dyn_param_state_sw, Exc.cb_dyn_param_state_sw, Loc_load.cb_dyn_param_state_sw ]
    
    # cb_dyn_state_dim::Int64 = length(get_non_null_list(cb_dyn_param_state_sw))


    dummy_vr1::Float64       = Exc.Ka * (v_ref - Gen.vh)
    dummy_vf_tilade::Float64 = dummy_vr1 * 0.97918

    func_system_matrices::Vector{Function} = Function[ get_a_plant_system_matrices!, get_a_plant_system_matrices ]
    
    func_update_system_matrices::Function = update_a_plant_system_matrices!
    
    # system_matrices = get_only_gen_avr_plant_system_matrices( Gen, Exc, dummy_vr1, dummy_vf_tilade  )

    # stability

    func_system_stab_matrices::Vector{Function} = Function[ get_a_plant_system_stability_matrices!, get_a_plant_system_stability_matrices ]
    
    func_update_system_stab_matrices::Function = update_a_plant_system_stability_matrices!
    
    # system_stab_matrices = get_only_gen_avr_plant_system_stability_matrices( Gen, Exc, dummy_vr1, dummy_vf_tilade  )
    
end

#--------------------------------------------------------
#--------------------------------------------------------

