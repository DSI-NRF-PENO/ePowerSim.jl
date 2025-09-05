# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

#####################################################
# ---------------------------------------------------
# componets callback
# ---------------------------------------------------
#####################################################


"""
    The callback unit is based on the infrastucture of callbacks in DifferentialEquations.jl. Components that have states with callbacks can define customised callbacks necessary for states. In addition, other set of callbacks can be added externally to states event list and states affect list which are  `list_selected_plants_state_event_cb_paras` and `list_selected_plants_state_affect_cb_paras` resepectively. Addition of callbacks externally is a simple procedure with specific formalism.

A typical event callback parameter requires specification of event function, the index of the state in which the event will occur in the system states indices and the value of the state that constiture and event. Their parameters are represented by `cb_state_event_func`
, `cb_event_in_state_idx_in_states` and `cb_event_trigger_value` respectively. Similarly, affect callback parameters requires specification of affect function, the index of the state that would be affected by the event in the system states indices, a function or value of the state in post event. It should be noted that in the design of callback unit for the package, an event can affect the state it occured in or another state. In addition, some events will require change of parameters in addition to a change in state. Such events will require additional information on the parameter that would be changed.

A illustration of data structures required for an anti-windup process, which is a typical operation in automatic voltage control in synchronous machine are  given as two namedtuple for event parameters and affect parameters. These are


`cb_event_parameters = (cb_state_event_func = cb_fun_condition_lower_lim, cb_event_in_state_idx_in_states = 5, cb_event_trigger_value = -4.16)`

`cb_(with_switch = true, idx = 1, cb_para = (cb_state_affect_func = cb_fun_affect_anti_windup_lower_lim!, cb_event_in_state_idx_in_states = 5, cb_affected_in_state_idx_in_states = 5, cb_post_event_state_value = -4.16))`


"""

function cb_fun_condition_lower_lim(
    u,
    t,
    integrator,
    event_func_param)
    
    (state_Idx,
     state_condition) =
         NamedTupleTools.select(
             event_func_param,
             (:cb_event_in_state_idx_in_states,
             :cb_event_trigger_value))

    u[ state_Idx] - state_condition

end


function cb_fun_condition_upper_lim(
    u,
    t,
    integrator,
    event_func_param)
    
    # (state_Idx,
    #  state_condition) =
    #      param_state_condition
    
    (state_Idx,
     state_condition) =
         NamedTupleTools.select(
             event_func_param,
             (:cb_event_in_state_idx_in_states,
             :cb_event_trigger_value))

    state_condition  - u[ state_Idx]

end


function cb_fun_affect_windup_lower_lim!(
    integrator,
    affect_cb_paras)
    
    # (state_Idx,
    #  state_value) =
    #      param_state_affect

     (state_Idx,
      state_value) =
          NamedTupleTools.select(
              affect_cb_paras,
              (:cb_affected_in_state_idx_in_states,
               :cb_post_event_state_value))
    
    
    integrator.u[ state_Idx] =  isa(state_value, Number) ?
        state_value :
        state_value(
            integrator.u[ state_Idx],
            integrator.t)

end


function cb_fun_affect_windup_upper_lim!(
    integrator,
    affect_cb_paras)
    
    # (state_Idx,
    #  state_value) =
    #      param_state_affect

     (state_Idx,
      state_value) =
          NamedTupleTools.select(
              affect_cb_paras,
              (:cb_affected_in_state_idx_in_states,
               :cb_post_event_state_value))

    " state_value can either be a function or a number"
    
    integrator.u[ state_Idx] =  isa(state_value, Number) ?
        state_value :
        state_value(
            integrator.u[ state_Idx],
            integrator.t)
    
end


"""
Note that the plants_switches is a vector of vectors
plants_switches = [[0], [0], [0]]

Each vector in `plants_switches`
"""
function cb_fun_affect_anti_windup_lower_lim!(
    integrator,
    affect_cb_paras)

    # (state_Idx,
    #  state_value,
    #  state_sw_Idx) =
    #      param_state_affect
    
    (cb_para,
     state_sw_Idx) =
         affect_cb_paras

     (state_Idx,
      state_value) =
          NamedTupleTools.select(
              cb_para,
              (:cb_affected_in_state_idx_in_states,
               
               :cb_post_event_state_value))
    

    integrator.u[ state_Idx] =  isa(state_value, Number) ?
        state_value :
        state_value(
            integrator.u[ state_Idx],
            integrator.t)

     @show "in affect_anti_windup_lower_lim!" "before setting x" state_sw_Idx integrator.u[state_Idx] 

    if integrator(integrator.t, Val{1})[ state_Idx] < 0.0

        integrator.p.plants_cb_paras_switches[
            state_sw_Idx][1] = 2

        @show "in affect_anti_windup_lower_lim!" "after setting dx" integrator.p.plants_cb_paras_switches[state_sw_Idx][1] state_sw_Idx integrator(integrator.t, Val{1})[state_Idx] integrator.u[state_Idx]
        
    end
    
            
end


function cb_fun_affect_anti_windup_upper_lim!(
    integrator,
    affect_cb_paras )
    
    # (state_Idx,
    #  state_value,
    #  state_sw_Idx) =
    #      param_state_affect

    
    (cb_para,
     state_sw_Idx) =
         affect_cb_paras

     (state_Idx,
      state_value) =
          NamedTupleTools.select(
              cb_para,
              (:cb_affected_in_state_idx_in_states,
               
               :cb_post_event_state_value))

    integrator.u[ state_Idx] =  isa(state_value, Number) ?
        state_value :
        state_value(
            integrator.u[ state_Idx],
            integrator.t)

    @show "in affect_anti_windup_upper_lim!" "before setting x" state_sw_Idx integrator.u[state_Idx]
    
    if integrator(integrator.t, Val{1})[
        state_Idx  ] > 0.0

        integrator.p.plants_cb_paras_switches[
            state_sw_Idx][1] = 1

         @show "in affect_anti_windup_upper_lim!" "after setting dx" integrator.p.plants_cb_paras_switches[state_sw_Idx][1] state_sw_Idx integrator(integrator.t, Val{1})[state_Idx] integrator.u[state_Idx]
        
    end
    
    
end


function cb_fun_state_events(
    out,
    u,
    t,
    integrator,
    list_selected_plants_state_event_cb_paras)

    for (idx, a_plant_state_event_cb_paras) in
        enumerate(
            list_selected_plants_state_event_cb_paras)
        
        # event_func =
        #     a_plant_state_event_cb_paras.cb_state_event_func

        event_func =
            getproperty(
                a_plant_state_event_cb_paras,
                :cb_state_event_func )

        event_func_param =
            NamedTupleTools.select(
        a_plant_state_event_cb_paras,
                (:cb_event_in_state_idx_in_states,
                 :cb_event_trigger_value))
        
         # event_func_param =
         #    (getproperty(a_plant_state_event_cb_paras,
         #                 :cb_event_in_state_idx_in_states ),
         #     getproperty(a_plant_state_event_cb_paras,
         #                 :cb_event_trigger_value ))
        
        out[idx] =
            event_func(
                u,
                t,
                integrator,
                event_func_param )
    end
end


function cb_fun_state_affect!(
    integrator,
    idx,
    list_selected_plants_state_affect_cb_paras )
    
    plant_state_affect_cb_paras =
        list_selected_plants_state_affect_cb_paras[idx]
    
    cb_para =
        getproperty(
            plant_state_affect_cb_paras,
            :cb_para )

    affect_func! = getproperty(
        getproperty(
            plant_state_affect_cb_paras,
            :cb_para ),
        :cb_state_affect_func)

    
    # with_switch =
    #     getproperty(
    #         plant_state_affect_cb_paras,
    #         :cb_with_switch )

    
    with_switch =
        getproperty(
            plant_state_affect_cb_paras,
            :with_switch )
    
    
    # state_sw_Idx =
    #     getproperty(
    #         plant_state_affect_cb_paras,
    #         :idx )

    
    state_sw_Idx = with_switch ?
        getproperty(
            plant_state_affect_cb_paras,
            :idx ) : []
    
    # (affect_func!,
    #  affected_in_state_idx_in_states,             
    #  state_value ) =
    #      NamedTupleTools.select(
    #         cb_para,
    #         (:cb_state_affect_func,
    #          :cb_affected_in_state_idx_in_states,
    
    #          :cb_post_event_state_value,
    #          ))

    
    affect_cb_paras =
        with_switch == true ?
        (cb_para,
         state_sw_Idx ) :
             (cb_para )
    
    affect_func!(
        integrator,
        affect_cb_paras)

end


function cb_fun_make_state_callbacks(
    list_selected_plants_state_event_cb_paras,
    list_selected_plants_state_affect_cb_paras )
    
    # Aggregation

    numbers_state_events =
        length(list_selected_plants_state_event_cb_paras )

    save_positions = (true, true)

    cb = VectorContinuousCallback(
        (out,
         u,
         t,
         integrator) -> cb_fun_state_events(
             out,
             u,
             t,
             integrator,
             list_selected_plants_state_event_cb_paras ),
        
        (integrator, idx) -> cb_fun_state_affect!(
            integrator,
            idx,
            list_selected_plants_state_affect_cb_paras),
        numbers_state_events;
        rootfind = SciMLBase.LeftRootFind,
        save_positions = save_positions )
    
    return cb
    
end



function cb_fun_make_state_callbacks(
    generic_model_callback_paras)
    # Aggregation

    (list_selected_plants_state_event_cb_paras,
     list_selected_plants_state_affect_cb_paras) =
         NamedTupleTools.select(
             generic_model_callback_paras,
             (:list_selected_plants_state_event_cb_paras,
              :list_selected_plants_state_affect_cb_paras))
    
    # list_state_event_func =
    #     [ cb_state_event_func...;
    #       cb_dyn_state_event_func...]
    
    # list_state_affect_func =
    #     [ cb_state_affect_func...;
    #       cb_dyn_state_affect_func...]
    
    # list_state_event_func_param =
    #     [ cb_state_event_param...;
    #       cb_dyn_state_event_param...]
    
    # list_state_affect_func_param =
    #     [ cb_state_affect_param...;
    #       cb_dyn_state_affect_param...] 

    numbers_state_events =
        length(
            list_selected_plants_state_event_cb_paras)

    save_positions = (true, true)

    cb = VectorContinuousCallback(
        (out,
         u,
         t,
         integrator) ->
             cb_fun_state_events(
                 out,
                 u,
                 t,
                 integrator,
                 list_selected_plants_state_event_cb_paras ),
        
        (integrator, idx) ->
            cb_fun_state_affect!(
                integrator,
                idx,
                list_selected_plants_state_affect_cb_paras),
        numbers_state_events;
        rootfind = SciMLBase.LeftRootFind,
        save_positions = save_positions )
    
    return cb
    
end


#---------------------------------------------------

"""

An event is associated with a state. An event in a
state (x_1) can affect the state (x_1) it occured in
or another state (x_2).

A state can have multiple events associated with it.

"""
function make_a_component_state_no_callback_paras(
    event_type_sym)
           
        return ( cb_event_type_sym = event_type_sym, )

end

function make_a_component_state_callback_paras(        
    event_type_sym,

    plant_state_idx_in_states,
    plant_states_syms,
    
    event_trigger_value_or_sym,
    post_event_state_value,
    
    event_in_state_sym,
    affected_in_state_sym,

    state_event_func,
    state_affect_func;
    
    comp_paras =
        nothing,
   
    component_parameter_switch =
        nothing )

    #----------------------------------------------
    
    cb_event_type_sym =
        event_type_sym
    
    cb_event_trigger_value =
        typeof(event_trigger_value_or_sym) == Symbol ?
        getproperty(
            comp_paras,
            event_trigger_value_or_sym) :
        event_trigger_value_or_sym
    
    cb_post_event_state_value =
        typeof(post_event_state_value) == Symbol ?
        (u,t) -> getproperty(
            comp_paras,
            event_trigger_value_or_sym) :
        post_event_state_value

    cb_event_in_state_sym =
        event_in_state_sym
    
    cb_affected_in_state_sym =
        affected_in_state_sym
             
    cb_state_event_func =
             state_event_func
         
    cb_state_affect_func =
             state_affect_func

    cb_event_in_state_idx_in_plant_states =
        findfirst((x) -> x == cb_event_in_state_sym,
         plant_states_syms)


    cb_affected_in_state_idx_in_plant_states =
        findfirst((x) -> x == cb_affected_in_state_sym,
         plant_states_syms)

    cb_event_in_state_idx_in_states =
        plant_state_idx_in_states[
            cb_event_in_state_idx_in_plant_states ]
    
    cb_affected_in_state_idx_in_states =
        plant_state_idx_in_states[
            cb_affected_in_state_idx_in_plant_states]

    if event_type_sym == :nothing
        
        return ( cb_event_type_sym = :nothing, )
    
    elseif event_type_sym ∈ [:lower_lim_anti_windup,
                             :upper_lim_anti_windup]
        
        cb_parameter_switch =
            component_parameter_switch

        cb_with_switch = true

        return nt_state_event_callback =
            (;cb_event_type_sym,

             cb_event_trigger_value,

             cb_post_event_state_value,

             cb_event_in_state_sym,
             cb_affected_in_state_sym,

             cb_state_event_func,
             cb_state_affect_func,

             cb_event_in_state_idx_in_states,
             cb_affected_in_state_idx_in_states,

             cb_parameter_switch,
             cb_with_switch)
    else
        
        cb_with_switch = false
        
        return nt_state_event_callback =
            (;cb_event_type_sym,

             cb_event_trigger_value,

             cb_post_event_state_value,

             cb_event_in_state_sym,
             cb_affected_in_state_sym,

             cb_state_event_func,
             cb_state_affect_func,

             cb_event_in_state_idx_in_states,
             cb_affected_in_state_idx_in_states,
             cb_with_switch)        
        
    end
    
end

#----------------------------------------
# callback Two-axis Synchronous generator
#----------------------------------------

function a_SM_plant_generic_model_callback_paras_func(
    plant_state_idx_in_states,
    plant_states_syms;
    kwd_para =
        plant_cb_paras_kwd_para )

    (gen_para,
     avr_para,
     gov_para,
     comps_callback_paras_fun,
     ) =
         kwd_para
             
    
    (gen_cb_paras_fun,
     avr_cb_paras_fun,
     gov_cb_paras_fun) =
         comps_callback_paras_fun
    
    #----------------------------------------

    # gen_fun__Sx_2axis_cb_v6__callback_paras
    
    gen_cb_paras =
        gen_cb_paras_fun(
            plant_state_idx_in_states,
            plant_states_syms;
            comp_paras =
                gen_para )
    
    #----------------------------------------

    # avr_fun__avr_x_cb__callback_paras
    
    avr_cb_paras =
        avr_cb_paras_fun(
            plant_state_idx_in_states,
            plant_states_syms;
            comp_paras =
                avr_para )
    
    avr_cb_parameter_switches =
        [ a_avr_cb_paras.cb_parameter_switch 
          for a_avr_cb_paras in
              avr_cb_paras if
                  a_avr_cb_paras.cb_with_switch ==
                      true ]
    
    #----------------------------------------

    # gov_fun__gov_x_cb__callback_paras
    
    gov_cb_paras =
        gov_cb_paras_fun(
            plant_state_idx_in_states,
            plant_states_syms;
            comp_paras =
                gov_para )
    
    gov_cb_parameter_switches =
        [ a_gov_cb_paras.cb_parameter_switch 
          for a_gov_cb_paras in
              gov_cb_paras if
                  a_gov_cb_paras.cb_with_switch ==
                      true ]
    
    #----------------------------------------

    # plant_cb_parameter_switches
    
     avr_gov_cb_para_sw_in_plant =
        [ avr_cb_parameter_switches;
          gov_cb_parameter_switches ]

    avr_gov_cb_para_sw_idx_in_plant =
        get_vars_or_paras_Idxs_in_flattend(
            [length(avr_cb_parameter_switches),
             length(gov_cb_parameter_switches)];
            dims_given = true )

    
    (avr_cb_para_sw_idx_in_plant,
     gov_cb_para_sw_idx_in_plant) =
         avr_gov_cb_para_sw_idx_in_plant
    
    avr_gov_cb_para_sw_idx_in_plant =
        (; avr_cb_para_sw_idx_in_plant,
         gov_cb_para_sw_idx_in_plant)
    
    #----------------------------------------

    plant_cb_paras = [
        gen_cb_paras;
        avr_cb_paras;
        gov_cb_paras ]
    
    #----------------------------------------

    return (;plant_cb_paras,
             avr_gov_cb_para_sw_in_plant,
             avr_gov_cb_para_sw_idx_in_plant )
end


#----------------------------------------
# callback Two-axis Synchronous condenser
#----------------------------------------


function a_SC_plant_generic_model_callback_paras_func(
    plant_state_idx_in_states,
    plant_states_syms;
    kwd_para =
        plant_cb_paras_kwd_para )

    (gen_para,
     avr_para,
     comps_callback_paras_fun,
     ) =
         kwd_para
             
    
    (gen_cb_paras_fun,
     avr_cb_paras_fun) =
         comps_callback_paras_fun
    
    #----------------------------------------

    gen_cb_paras =
        gen_cb_paras_fun(
            plant_state_idx_in_states,
            plant_states_syms;
            comp_paras =
                gen_para )
    
    #----------------------------------------

    avr_cb_paras =
        avr_cb_paras_fun(
            plant_state_idx_in_states,
            plant_states_syms;
            comp_paras =
                avr_para )
    

    avr_cb_parameter_switches =
        [ a_avr_cb_paras.cb_parameter_switch
          for a_avr_cb_paras in
              avr_cb_paras if
                  a_avr_cb_paras.cb_with_switch ==
                      true]

    
    avr_gov_cb_para_sw_in_plant =
        avr_cb_parameter_switches


    avr_gov_cb_para_sw_idx_in_plant =
        get_vars_or_paras_Idxs_in_flattend(
            [length(avr_cb_parameter_switches),
             0];
            dims_given = true )
    
    (avr_cb_para_sw_idx_in_plant,
     gov_cb_para_sw_idx_in_plant) =
         avr_gov_cb_para_sw_idx_in_plant
    
    avr_gov_cb_para_sw_idx_in_plant =
        (; avr_cb_para_sw_idx_in_plant,
         gov_cb_para_sw_idx_in_plant)
    
    #----------------------------------------

    plant_cb_paras = [
        gen_cb_paras;
        avr_cb_paras]
    
    #----------------------------------------

    return (;plant_cb_paras,
            avr_gov_cb_para_sw_in_plant,
            avr_gov_cb_para_sw_idx_in_plant)
end

# cb_parameter_switch


#-----------------------------------------------------


function plants_generic_model_callback_paras_func(
    state_vars_idx,
    plants_states_syms ;
    kwd_para =
        plants_cb_paras_kwd_para )

    (gens_para,
     avrs_para,
     govs_para,
     comps_callback_paras_funs
     ) =
         kwd_para
             
    gens_cb_paras_fun =
        [an_item.gen_cb_paras_fun for an_item in
             comps_callback_paras_funs]

    avrs_cb_paras_fun = 
        [an_item.avr_cb_paras_fun for an_item in
             comps_callback_paras_funs]

    govs_cb_paras_fun = 
        [an_item.gov_cb_paras_fun for an_item in
             comps_callback_paras_funs]

    plants_cb_paras_wt_switches_paras =
        [ Symbol(split(String(
            nameof(gov_fun)),"__")[2]) ==
                :nothing ?
            a_SC_plant_generic_model_callback_paras_func(
                plant_state_idx_in_states,
                plant_states_syms;
                    kwd_para =
                        (gen_para,
                         avr_para,
                         (gen_fun,
                          avr_fun ))) :
           a_SM_plant_generic_model_callback_paras_func(
               plant_state_idx_in_states,
               plant_states_syms;
               kwd_para =
                   (gen_para,
                    avr_para,
                    gov_para,
                    (gen_fun,
                     avr_fun,
                     gov_fun) ))
          for (plant_state_idx_in_states,
               plant_states_syms,
               gen_para,
               avr_para,
               gov_para,
               gen_fun,
               avr_fun,
               gov_fun ) in
              zip( state_vars_idx,
                   plants_states_syms,
                   
                   gens_para,
                   avrs_para,
                   govs_para,
                   
                   gens_cb_paras_fun,
                   avrs_cb_paras_fun,
                   govs_cb_paras_fun) ]

    #----------------------------------------

    "First flatten plants_cb_paras, and then filter
     out SM gens and SC gens. Their cb paras are
     :nothing since they have no real callbacks"
    
    plants_cb_paras =
        [[ a_plant_cb_para.plant_cb_paras
         for a_plant_cb_para in
             plants_cb_paras_wt_switches_paras]...;]
    
    plants_cb_paras =
        [ a_plant_cb_para
         for a_plant_cb_para in
             plants_cb_paras
             if a_plant_cb_para.cb_event_type_sym !=
                 :nothing ]
        
    plants_avr_gov_cb_para_sw_in_plant =
        [ a_plant_cb_para.avr_gov_cb_para_sw_in_plant
         for a_plant_cb_para in
             plants_cb_paras_wt_switches_paras]
        
    plants_avr_gov_cb_para_sw_idx_in_plant =
        [ a_plant_cb_para.avr_gov_cb_para_sw_idx_in_plant
         for a_plant_cb_para in
             plants_cb_paras_wt_switches_paras]

    #----------------------------------------


    plants_cb_paras_no_switches =     
        [ a_plant_cb_para
         for a_plant_cb_para in
             plants_cb_paras
             if a_plant_cb_para.cb_with_switch ==
                 false ]
    
    plants_cb_paras_wt_switches =     
        [ a_plant_cb_para
         for a_plant_cb_para in
             plants_cb_paras
             if a_plant_cb_para.cb_with_switch ==
                 true ]

    plants_cb_paras_wt_switches_idx =     
        [ idx
         for (idx, a_plant_cb_para) in
             enumerate(plants_cb_paras)
             if a_plant_cb_para.cb_with_switch ==
                 true ]

    "need to properly index plants_cb_paras_wt_switches"
    
    wt_switches_idx2idx = OrderedDict{Int64, Int64}(
        a_wt_switches_idx => a_idx
        for ( a_idx, a_wt_switches_idx) in
            enumerate(plants_cb_paras_wt_switches_idx) )
    
    #----------------------------------------

    plants_cb_paras_switches =
        [plants_avr_gov_cb_para_sw_in_plant...;]        

    plants_cb_paras_switches_idx =
        [ idx for idx in
             1:length(plants_cb_paras_switches)]

    #----------------------------------------

    # list_state_event_func
    
    list_selected_plants_state_event_cb_paras =
        [ NamedTupleTools.select(
            a_cb_paras,
            (:cb_state_event_func,
             :cb_event_in_state_idx_in_states,
             
             :cb_event_trigger_value,
             )) 
          for a_cb_paras in plants_cb_paras ]
    
    # list_state_affect_func
    
    # list_selected_plants_state_affect_cb_paras =
    #     [ a_cb_paras.cb_with_switch == true ?
    #     (with_switch = a_cb_paras.cb_with_switch,
    #      idx = wt_switches_idx2idx[ a_idx],
    #      cb_para = NamedTupleTools.select(
    #         a_cb_paras,
    #         (:cb_state_affect_func,
    #          :cb_event_in_state_idx_in_states,
    #          :cb_affected_in_state_idx_in_states,
             
    #          :cb_post_event_state_value,
    #          )) ) :
    #              (with_switch = a_cb_paras.cb_with_switch,
    #               cb_para = NamedTupleTools.select(
    #                  a_cb_paras,
    #                  (:cb_state_event_func,
    #                   :cb_event_in_state_idx_in_states,
    #                   :cb_affected_in_state_idx_in_states,
    #                   :cb_event_trigger_value,
    #                   :cb_post_event_state_value,
    #                   )) )
    #       for (a_idx, a_cb_paras) in
    #           enumerate(plants_cb_paras) ]

    
    list_selected_plants_state_affect_cb_paras =
        [ a_cb_paras.cb_with_switch == true ?
        (with_switch = a_cb_paras.cb_with_switch,
         idx = wt_switches_idx2idx[ a_idx],
         cb_para = NamedTupleTools.select(
            a_cb_paras,
            (:cb_state_affect_func,
             :cb_event_in_state_idx_in_states,
             :cb_affected_in_state_idx_in_states,
             
             :cb_post_event_state_value,
             )) ) :
                 (with_switch = a_cb_paras.cb_with_switch,
                  cb_para = NamedTupleTools.select(
                     a_cb_paras,
                      (:cb_state_affect_func,
                       
                       :cb_state_event_func,
                       :cb_event_in_state_idx_in_states,
                       :cb_affected_in_state_idx_in_states,
                       :cb_event_trigger_value,
                       :cb_post_event_state_value,
                      )) )
          for (a_idx, a_cb_paras) in
              enumerate(plants_cb_paras) ]
    
    #----------------------------------------
    
    return (;list_selected_plants_state_event_cb_paras,
            list_selected_plants_state_affect_cb_paras,
            plants_cb_paras,
            plants_avr_gov_cb_para_sw_in_plant,
            plants_avr_gov_cb_para_sw_idx_in_plant, 

            plants_cb_paras_wt_switches,
            plants_cb_paras_no_switches,
            
            plants_cb_paras_switches,
            plants_cb_paras_switches_idx)
    
end

#---------------------------------------------------
#---------------------------------------------------

