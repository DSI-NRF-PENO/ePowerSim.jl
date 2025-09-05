# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

#########################################################
# ------------------------------------------------------
# callbacks with values or condition based on functions 
# ------------------------------------------------------
#########################################################

function fun_condition_lower_lim(
    u,
    t,
    integrator,
    param_state_condition)
    
    states_Idx, state_condition = param_state_condition
    
    if Symbol(typeof(states_Idx)) == :Int64
        u[ states_Idx ] - state_condition
    else
        u[ states_Idx[1] ] - state_condition
    end
    

end


function fun_condition_upper_lim(
    u,
    t,
    integrator,
    param_state_condition)
    
    states_Idx, state_condition =  param_state_condition

    if Symbol(typeof(states_Idx)) == :Int64
        
        state_condition  - u[ states_Idx ]
        
    else
        
        state_condition  - u[ states_Idx[1] ]
        
    end
    

end


function fun_affect_windup_lower_lim!(
    integrator,
    param_state_affect)

    states_Idx, state_value =  param_state_affect

    if Symbol(typeof(states_Idx)) == :Int64
        
        integrator.u[states_Idx] =
            state_value(integrator.u[states_Idx],
                        integrator.t)

        # integrator.u[ states_Idx ]  = state_value 

        @show  states_Idx  "in affect_windup_lower_lim!"
        
    else
        integrator.u[states_Idx[2]] =
            state_value(integrator.u[states_Idx[2]],
                        integrator.t)

        # integrator.u[ states_Idx ]  = state_value 

        @show  states_Idx  "in affect_windup_lower_lim!"
    end

end


function fun_affect_windup_upper_lim!(
    integrator,
    param_state_affect)

    states_Idx, state_value =  param_state_affect
    
    if Symbol(typeof(states_Idx)) == :Int64
        
        integrator.u[states_Idx] =
            state_value(integrator.u[states_Idx],
                        integrator.t)

        # integrator.u[ states_Idx ]  = state_value 

        @show  states_Idx  "in affect_windup_lower_lim!"
        
    else
        
        integrator.u[states_Idx[2]] =
            state_value(integrator.u[states_Idx[2]],
                        integrator.t)

        # integrator.u[ states_Idx ]  = state_value 

        @show  states_Idx  "in affect_windup_lower_lim!"
    end    

end


function fun_affect_anti_windup_lower_lim!(
    integrator,
    param_state_affect)

    (states_Idx,
     state_value,
     comp_idx,
     state_sw_Idx) =
        param_state_affect

    if Symbol(typeof(states_Idx)) == :Int64
        
        @show "in affect_anti_windup_lower_lim!" "before setting x" states_Idx integrator.u[states_Idx] 

        integrator.u[states_Idx] =
            state_value(
                integrator.u[states_Idx],
                integrator.t)

        @show "after setting x" integrator.u[states_Idx] integrator.p[2][1][2][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx]


        if integrator(integrator.t, Val{1})[
            states_Idx] < 0.0

            # integrator.p[1][ state_sw_Idx ] = 2

            """
               p[2][1][2] is based on the order
               nodes_cb_sw in parameters passed to
               nodes_edges_dynamics.
               e.g.

                   [2]
               nd, nodes_and_edges_p = p
                [1]
               sd_nodes_param, sd_edges_param =
                   nodes_and_edges_p
                               [2]
               nodes_dyn_func, nodes_cb_sw, ... =
                   sd_nodes_param

            """
            integrator.p[2][1][2][comp_idx][
                state_sw_Idx] = 2 
        end

         @show "after setting dx" integrator.p[2][1][2][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx]

    else

        """ In this case states_Idx is a tuple =
                (states_Idx[1], states_Idx[2]).

           We are using an event in
           states_Idx[1] to change states_Idx[2]
        """
        
        @show "in affect_anti_windup_lower_lim!" "before setting x" states_Idx integrator.u[states_Idx[2]] 

        integrator.u[states_Idx[2]] =
            state_value(integrator.u[states_Idx[2]],
                        integrator.t)

        @show "after setting x" integrator.u[states_Idx[2]] integrator.p[2][1][2][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx[2]]


        if integrator(integrator.t, Val{1})[
            states_Idx[2]  ] < 0.0

            # integrator.p[1][ state_sw_Idx ] = 2

            """
               p[2][1][2] is based on the order
               nodes_cb_sw in parameters passed to
               nodes_edges_dynamics.
               e.g.

                   [2]
               nd, nodes_and_edges_p = p
                [1]
               sd_nodes_param, sd_edges_param =
                   nodes_and_edges_p
                               [2]
               nodes_dyn_func, nodes_cb_sw, ... =
                   sd_nodes_param

            """
            integrator.p[2][1][2][comp_idx][
                state_sw_Idx] = 2 
        end

         @show "after setting dx" integrator.p[2][1][2][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx[2]]
        
    end
            
end


function fun_affect_anti_windup_upper_lim!(
    integrator,
    param_state_affect)

    (states_Idx,
     state_value,
     comp_idx,
     state_sw_Idx) =
        param_state_affect

    if Symbol(typeof(states_Idx)) == :Int64
        
        @show "in affect_anti_windup_lower_lim!" "before setting x" states_Idx integrator.u[states_Idx] 

        integrator.u[states_Idx]  =
            state_value(integrator.u[states_Idx],
                        integrator.t)

        @show "after setting x" integrator.u[states_Idx] integrator.p[2][1][2][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx]


        if integrator(integrator.t, Val{1})[
            states_Idx  ] > 0.0

            # integrator.p[1][ state_sw_Idx ] = 2

            """
               p[2][1][2] is based on the order
               nodes_cb_sw in parameters passed to
               nodes_edges_dynamics.
               e.g.

                   [2]
               nd, nodes_and_edges_p = p
                [1]
               sd_nodes_param, sd_edges_param =
                   nodes_and_edges_p
                               [2]
               nodes_dyn_func, nodes_cb_sw, ... =
                   sd_nodes_param

            """
            integrator.p[2][1][2][comp_idx][
                state_sw_Idx] = 1 
        end

         @show "after setting dx" integrator.p[2][1][2][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx]

    else

        """ In this case states_Idx is a tuple =
            (states_Idx[1], states_Idx[2]).
           We are using an event in
           states_Idx[1] to change states_Idx[2]
        """
        
        @show "in affect_anti_windup_lower_lim!" "before setting x" states_Idx integrator.u[states_Idx[2]] 

        integrator.u[states_Idx[2]]  =
            state_value(integrator.u[states_Idx[2]],
                        integrator.t)

        @show "after setting x" integrator.u[states_Idx[2]] integrator.p[2][1][2][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx[2]]


        if integrator(integrator.t, Val{1})[
            states_Idx[2]  ] > 0.0

            # integrator.p[1][ state_sw_Idx ] = 2

            """
               p[2][1][2] is based on the order
               nodes_cb_sw in parameters passed to
               nodes_edges_dynamics.
               e.g.

                   [2]
               nd, nodes_and_edges_p = p
                [1]
               sd_nodes_param, sd_edges_param =
                   nodes_and_edges_p
                               [2]
               nodes_dyn_func, nodes_cb_sw, ... =
                   sd_nodes_param

            """
            integrator.p[2][1][2][comp_idx][
                state_sw_Idx] = 1 
        end

         @show "after setting dx" integrator.p[2][1][2][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx[2]]
        
    end
    
    
end


function fun_state_events(
    out,
    u,
    t,
    integrator,
    param_state_condition_events)

    (list_state_event_func,
     list_state_event_func_param)  =
         param_state_condition_events  

    for (idx, an_event_func) in
        enumerate(list_state_event_func)
        
        out[idx] =
            an_event_func(
                u,
                t,
                integrator,
                list_state_event_func_param[ idx ] )
    end
end


function fun_state_affect!(
    integrator,
    idx,
    param_state_events_affect)

    (list_state_affect_func,
     list_state_affect_func_param) =
         param_state_events_affect 

    affect_func = list_state_affect_func[idx]

    affect_func(
        integrator,
        list_state_affect_func_param[ idx ] )

end


function fun_make_state_callbacks(
    plants_collection)
    
    dim_nodes_states =
        [a_node.dim
         for a_node in
             plants_collection ]
    
    nodes_states_offset =
        create_offsets( dim_nodes_states )
    
    nodes_states_Idx =
        create_idxs(
            nodes_states_offset,
            dim_nodes_states )

    cb_state_event_func =
        [ a_node.cb_state_event_func
          for a_node in
              plants_collection ]
    
    cb_state_event_func =
        [ cb_state_event_func...;]
    
    cb_state_affect_func =
        [ a_node.cb_state_affect_func
          for a_node in
              plants_collection ]
    
    cb_state_affect_func =
        [ cb_state_affect_func...;]
    
    cb_state_conditions =
        [ a_node.cb_state_conditions
          for a_node in
              plants_collection ]
    
    cb_state_conditions =
        [ cb_state_conditions...; ]
    
    cb_state_values =
        [ a_node.cb_state_values
          for a_node in
              plants_collection ]
    
    cb_state_values =
        [ cb_state_values...;]
    
    cb_state_sym2Idx =
        [ a_node.cb_state_sym2Idx
          for a_node in
              plants_collection ]
    
    cb_state_Idx =
        [ [ Symbol(typeof(state_Idx)) == :Int64 ?
        nodes_states_Idx[ Ind ][ state_Idx ] :
        ( nodes_states_Idx[ Ind ][ state_Idx[1] ],
          nodes_states_Idx[ Ind ][ state_Idx[2] ])
            for state_Idx in
                cb_state_sym2Idx[ Ind ] ]
          for Ind in
              collect( 1:length( nodes_states_Idx ) )
              if cb_state_sym2Idx[ Ind ] != []]
    
    cb_state_Idx =
        [ cb_state_Idx...;]

    cb_state_event_param =
        [ (state_Idx, condition)
          for (state_Idx, condition) in
              zip(cb_state_Idx,
                  cb_state_conditions) ]
    

    cb_state_affect_param =
        [(state_Idx, value)
         for (state_Idx, value) in
             zip(cb_state_Idx,
                 cb_state_values) ]
    
    cb_dyn_state_event_func =
        [ a_node.cb_dyn_state_event_func
          for a_node in
              plants_collection ]
    
    cb_dyn_state_event_func =
        [ cb_dyn_state_event_func...;]
    
    cb_dyn_state_affect_func =
        [ a_node.cb_dyn_state_affect_func
          for a_node in
              plants_collection ]
    
    cb_dyn_state_affect_func =
        [ cb_dyn_state_affect_func...; ]
    
    cb_dyn_state_conditions  =
        [ a_node.cb_dyn_state_conditions
          for a_node in
              plants_collection ]
    
    cb_dyn_state_conditions =
        [ cb_dyn_state_conditions...;]
    
    cb_dyn_state_values =
        [ a_node.cb_dyn_state_values
          for a_node in
              plants_collection ]
    
    cb_dyn_state_values =
        [ cb_dyn_state_values...;]
    
    cb_dyn_state_sym2Idx =
        [ a_node.cb_dyn_state_sym2Idx
          for a_node in
              plants_collection ]
    
    cb_dyn_state_Idx =
        [ [Symbol(typeof(state_Idx)) == :Int64 ?
        nodes_states_Idx[ Ind ][ state_Idx ] :
        (nodes_states_Idx[ Ind ][ state_Idx[1] ],
         nodes_states_Idx[ Ind ][ state_Idx[2] ])
           for state_Idx in cb_dyn_state_sym2Idx[ Ind ]]
          for Ind in collect( 1:length( nodes_states_Idx ))
              if cb_dyn_state_sym2Idx[ Ind ] != [] ]
    
    cb_dyn_state_Idx =
        [ cb_dyn_state_Idx...;]


    cb_dyn_state_sw_Idx =
        [ a_node.cb_dyn_state_sw_Idx
          for a_node in plants_collection ]
    
    # This contains a tuple (comp_idx,  state_sw_Idx) of comp_idx and  state_sw_Idx in
    # `cb_dyn_state_sw_Idx` for each components that have non-null `cb_dyn_state_sw_Idx`

    # comps_cb_dyn_state_sw_para_Idx = [ [(comp_idx, sw_idx) for sw_idx in sw_array] for (comp_idx, sw_array ) in enumerate(cb_dyn_state_sw_Idx) if cb_dyn_state_sw_Idx[comp_idx] != [] ]

    # --------------------------------------
    
    # sw_idx == 1 ? 1 : sw_idx * 2 - 1 is used since all
    # swiches for a component are flattened and each switch
    # contains two items.

    comps_cb_dyn_state_sw_para_Idx =
        [ [(comp_idx, sw_idx == 1 ?
        1 :
        sw_idx * 2 - 1)
           for sw_idx in sw_array]
          for (comp_idx, sw_array ) in
              enumerate(cb_dyn_state_sw_Idx)
              if cb_dyn_state_sw_Idx[comp_idx] != [] ]

    # --------------------------------------
    
    comps_cb_dyn_state_sw_para_Idx =
        [ comps_cb_dyn_state_sw_para_Idx...;]
    
    cb_dyn_state_event_param =
        [ (state_Idx, condition)
          for (state_Idx, condition) in
              zip(cb_dyn_state_Idx,
                  cb_dyn_state_conditions) ]
    
    cb_dyn_state_affect_param =
        [ (state_Idx, value, comp_idx, state_sw_Idx)
          for (state_Idx,value,comp_idx, state_sw_Idx) in
              zip(cb_dyn_state_Idx,
                  cb_dyn_state_values,
                  first.(comps_cb_dyn_state_sw_para_Idx),
                  last.(comps_cb_dyn_state_sw_para_Idx) ) ]

    # Aggregation
    
    list_state_event_func =
        [ cb_state_event_func...;
          cb_dyn_state_event_func...]
    
    list_state_affect_func =
        [ cb_state_affect_func...;
          cb_dyn_state_affect_func...]
    
    list_state_event_func_param =
        [ cb_state_event_param...;
          cb_dyn_state_event_param...]
    
    list_state_affect_func_param =
        [ cb_state_affect_param...;
          cb_dyn_state_affect_param...] 

    numbers_state_events =
        length(list_state_event_func )

    save_positions = (true, true)

    cb = VectorContinuousCallback(
        (out,
         u,
         t,
         integrator) ->
             fun_state_events(
                 out,
                 u,
                 t,
                 integrator,
                 (list_state_event_func,
                  list_state_event_func_param )),
        
        (integrator, idx) ->
            fun_state_affect!(
                integrator,
                idx,
                (list_state_affect_func ,
                 list_state_affect_func_param)),
        numbers_state_events;
        rootfind = SciMLBase.LeftRootFind,
        save_positions = save_positions )
    
    return cb
    
end

#----------------------------------------------
## For a single plant
#----------------------------------------------

function plant_fun_affect_anti_windup_lower_lim!(
    integrator,
    param_state_affect)

    (states_Idx,
     state_value,
     comp_idx,
     state_sw_Idx) =
        param_state_affect

    if Symbol(typeof(states_Idx)) == :Int64

        @show "in affect_anti_windup_lower_lim!" "before setting x" states_Idx integrator.u[states_Idx] comp_idx state_sw_Idx

        integrator.u[states_Idx] =
            state_value(integrator.u[states_Idx],
                        integrator.t)

        @show "after setting x" integrator.u[states_Idx] integrator.p[1][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx]


        if integrator(integrator.t, Val{1})[states_Idx] < 0.0

            # integrator.p[1][ state_sw_Idx ] = 2

            """
               p[2][1][2] is based on the order
               nodes_cb_sw in parameters passed to
               nodes_edges_dynamics.
               e.g.

                   [2]
               nd, nodes_and_edges_p = p
                [1]
               sd_nodes_param, sd_edges_param =
                   nodes_and_edges_p
                               [2]
               nodes_dyn_func, nodes_cb_sw, ... =
                   sd_nodes_param

            """
            integrator.p[1][comp_idx][state_sw_Idx] = 2 
        end

        @show "after setting dx" integrator.p[1][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx]
    else

        @show "in affect_anti_windup_lower_lim!" "before setting x" states_Idx integrator.u[states_Idx[2]] comp_idx state_sw_Idx

        integrator.u[states_Idx[2]]  = state_value(integrator.u[states_Idx[2]], integrator.t)

        @show "after setting x" integrator.u[states_Idx[2]] integrator.p[1][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx[2]]


        if integrator(integrator.t, Val{1})[ states_Idx[2]  ] < 0.0

            # integrator.p[1][ state_sw_Idx ] = 2

            """
               p[2][1][2] is based on the order
               nodes_cb_sw in parameters passed to
               nodes_edges_dynamics.
               e.g.

                   [2]
               nd, nodes_and_edges_p = p
                [1]
               sd_nodes_param, sd_edges_param =
                   nodes_and_edges_p
                               [2]
               nodes_dyn_func, nodes_cb_sw, ... =
                   sd_nodes_param

            """
            integrator.p[1][comp_idx][state_sw_Idx] = 2 
        end

        @show "after setting dx" integrator.p[1][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx[2]]        
        
    end
    

end


function plant_fun_affect_anti_windup_upper_lim!(
    integrator,  param_state_affect)

    (states_Idx,
     state_value,
     comp_idx,
     state_sw_Idx) =
        param_state_affect

    if Symbol(typeof(states_Idx)) == :Int64

        @show "in affect_anti_windup_lower_lim!" "before setting x" states_Idx integrator.u[states_Idx] comp_idx state_sw_Idx

        integrator.u[states_Idx] =
            state_value(integrator.u[states_Idx],
                        integrator.t)

        @show "after setting x" integrator.u[states_Idx] integrator.p[1][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx]


        if integrator(integrator.t, Val{1})[
            states_Idx  ] > 0.0

            # integrator.p[1][ state_sw_Idx ] = 2

            """
               p[2][1][2] is based on the order
               nodes_cb_sw in parameters passed to
               nodes_edges_dynamics.
               e.g.

                   [2]
               nd, nodes_and_edges_p = p
                [1]
               sd_nodes_param, sd_edges_param =
                    nodes_and_edges_p
                               [2]
               nodes_dyn_func, nodes_cb_sw, ... =
                   sd_nodes_param

            """
            integrator.p[1][comp_idx][state_sw_Idx] = 1
        end

        @show "after setting dx" integrator.p[1][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx]
    else

        @show "in affect_anti_windup_lower_lim!" "before setting x" states_Idx integrator.u[states_Idx[2]] comp_idx state_sw_Idx

        integrator.u[states_Idx[2]]  = state_value(integrator.u[states_Idx[2]], integrator.t)

        @show "after setting x" integrator.u[states_Idx[2]] integrator.p[1][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx[2]]


        if integrator(integrator.t, Val{1})[
            states_Idx[2]  ] > 0.0

            # integrator.p[1][ state_sw_Idx ] = 2

            """
               p[2][1][2] is based on the order
               nodes_cb_sw in parameters passed to
               nodes_edges_dynamics.
               e.g.

                   [2]
               nd, nodes_and_edges_p = p
                [1]
               sd_nodes_param, sd_edges_param =
                   nodes_and_edges_p
                               [2]
               nodes_dyn_func, nodes_cb_sw, ... =
                   sd_nodes_param

            """
            integrator.p[1][comp_idx][
                state_sw_Idx] = 1
        end

        @show "after setting dx" integrator.p[1][comp_idx][state_sw_Idx] integrator(integrator.t, Val{1})[states_Idx[2]]        
        
    end
    
end


function plant_fun_state_events(
    out, u, t, integrator,
    param_state_condition_events)

    (list_state_event_func, list_state_event_func_param) =
         param_state_condition_events  

    for (idx, an_event_func) in enumerate(
        list_state_event_func)
        
        out[idx] = an_event_func(
            u, t, integrator,
            list_state_event_func_param[ idx] )
    end
end


function plant_fun_state_affect!(
    integrator, idx,
    param_state_events_affect)

    (list_state_affect_func,
     list_state_affect_func_param) =
         param_state_events_affect 

    affect_func = list_state_affect_func[idx]

    affect_func(integrator,
                list_state_affect_func_param[ idx ] )

end


# plants_collection = [plant]

function plant_fun_make_state_callbacks(plants_collection)
    

    dim_nodes_states =
        [a_node.dim for a_node in plants_collection ]
    
    nodes_states_offset =
        create_offsets( dim_nodes_states )
    
    nodes_states_Idx =
        create_idxs( nodes_states_offset,
                     dim_nodes_states )

    cb_state_event_func =
        [ a_node.cb_state_event_func for a_node in
             plants_collection ]
    
    cb_state_event_func =
        [ cb_state_event_func...;]
    
    cb_state_affect_func =
        [ a_node.cb_state_affect_func for a_node in
             plants_collection ]
    
    cb_state_affect_func =
        [ cb_state_affect_func...;]
    
    cb_state_conditions =
        [ a_node.cb_state_conditions for a_node in
             plants_collection ]
    
    cb_state_conditions =
        [ cb_state_conditions...; ]
    
    cb_state_values =
        [a_node.cb_state_values for a_node in
             plants_collection ]
    
    cb_state_values =
        [ cb_state_values...;]
    
    cb_state_sym2Idx =
        [ a_node.cb_state_sym2Idx for a_node in
             plants_collection ]
    
    # cb_state_Idx             = [ [ nodes_states_Idx[ Ind ][ state_Idx ] for state_Idx in cb_state_sym2Idx[ Ind ] ] for
    #                                 Ind in collect( 1:length( nodes_states_Idx ) )  if cb_state_sym2Idx[ Ind ] != []]

    cb_state_Idx = [ [ Symbol(typeof(state_Idx)) == :Int64 ? nodes_states_Idx[ Ind ][ state_Idx ] : ( nodes_states_Idx[ Ind ][ state_Idx[1] ],  nodes_states_Idx[ Ind ][ state_Idx[2] ])  for state_Idx in cb_state_sym2Idx[ Ind ] ] for
                                     Ind in collect( 1:length( nodes_states_Idx ) )  if cb_state_sym2Idx[ Ind ] != []]
    
    cb_state_Idx =
        [ cb_state_Idx...;]

    cb_state_event_param =
        [(state_Idx,condition) for (state_Idx,condition) in
             zip(cb_state_Idx, cb_state_conditions) ]
    

    cb_state_affect_param =
        [(state_Idx, value) for (state_Idx,  value) in
             zip(cb_state_Idx,  cb_state_values) ]

    
    cb_dyn_state_event_func  =
        [ a_node.cb_dyn_state_event_func for a_node in
             plants_collection ]
    
    cb_dyn_state_event_func =
        [ cb_dyn_state_event_func...;]
    
    # cb_dyn_state_affect_func = [ a_node.cb_dyn_state_affect_func for a_node in plants_collection ]

    cb_dyn_state_affect_func = [[ [plant_fun_affect_anti_windup_lower_lim!, plant_fun_affect_anti_windup_upper_lim!] for k in 1:length(a_node.cb_dyn_state_affect_func)] for a_node in plants_collection  if length(a_node.cb_dyn_state_affect_func) != 0 ]
    
    cb_dyn_state_affect_func =
        [ cb_dyn_state_affect_func...; ]
    
    cb_dyn_state_conditions  =
        [ a_node.cb_dyn_state_conditions for a_node in
             plants_collection ]
    
    cb_dyn_state_conditions  =
        [ cb_dyn_state_conditions...;]
    
    cb_dyn_state_values =
        [ a_node.cb_dyn_state_values for a_node in
             plants_collection ]
    
    cb_dyn_state_values =
        [ cb_dyn_state_values...;]
    
    cb_dyn_state_sym2Idx =
        [ a_node.cb_dyn_state_sym2Idx for a_node in
             plants_collection ]
    
    # cb_dyn_state_Idx         = [ [ nodes_states_Idx[ Ind ][ state_Idx ] for state_Idx in cb_dyn_state_sym2Idx[ Ind ] ] for
    #                                 Ind in collect( 1:length( nodes_states_Idx ) )  if cb_dyn_state_sym2Idx[ Ind ] != [] ]

    cb_dyn_state_Idx = [ [Symbol(typeof(state_Idx)) == :Int64 ? nodes_states_Idx[ Ind ][ state_Idx ] : (nodes_states_Idx[ Ind ][ state_Idx[1] ], nodes_states_Idx[ Ind ][ state_Idx[2] ])  for state_Idx in cb_dyn_state_sym2Idx[ Ind ] ] for
                                    Ind in collect( 1:length( nodes_states_Idx ) )  if cb_dyn_state_sym2Idx[ Ind ] != [] ]
    
    cb_dyn_state_Idx =
        [ cb_dyn_state_Idx...;]


    cb_dyn_state_sw_Idx =
        [ a_node.cb_dyn_state_sw_Idx for a_node in
             plants_collection ]
    
    # This contains a tuple (comp_idx,  state_sw_Idx) of comp_idx and  state_sw_Idx in
    # `cb_dyn_state_sw_Idx` for each components that have non-null `cb_dyn_state_sw_Idx`

    comps_cb_dyn_state_sw_para_Idx = [ [(comp_idx, sw_idx) for sw_idx in sw_array] for ( comp_idx, sw_array ) in enumerate(cb_dyn_state_sw_Idx) if cb_dyn_state_sw_Idx[comp_idx] != [] ]

        
    comps_cb_dyn_state_sw_para_Idx  =
        [ comps_cb_dyn_state_sw_para_Idx...;]
    
    cb_dyn_state_event_param  =
        [(state_Idx,condition) for (state_Idx,condition) in
             zip(cb_dyn_state_Idx,
                 cb_dyn_state_conditions) ]
    
    cb_dyn_state_affect_param = [ (state_Idx, value, comp_idx, state_sw_Idx) for (state_Idx, value, comp_idx, state_sw_Idx) in zip(cb_dyn_state_Idx,  cb_dyn_state_values, first.(comps_cb_dyn_state_sw_para_Idx), last.(comps_cb_dyn_state_sw_para_Idx) ) ]

    # Aggregation
    list_state_event_func =
        [ cb_state_event_func...;
          cb_dyn_state_event_func...]
    
    list_state_affect_func =
        [ cb_state_affect_func...;
          cb_dyn_state_affect_func...]
    
    list_state_event_func_param =
        [ cb_state_event_param...;
          cb_dyn_state_event_param...]
    
    list_state_affect_func_param =
        [ cb_state_affect_param...;
          cb_dyn_state_affect_param...] 

    numbers_state_events =
        length( list_state_event_func )

    save_positions = (true, true)

    cb = VectorContinuousCallback(
        (out, u, t, integrator) -> plant_fun_state_events(
            out, u, t, integrator,
            ( list_state_event_func ,
              list_state_event_func_param )),
        (integrator, idx) -> plant_fun_state_affect!(
            integrator, idx,
            ( list_state_affect_func ,
              list_state_affect_func_param )),
        numbers_state_events;
        rootfind=SciMLBase.LeftRootFind,
        save_positions=save_positions )
    
    return cb
    
end


# #--------------------------------
