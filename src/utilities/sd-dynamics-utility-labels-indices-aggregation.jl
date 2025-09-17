# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247


#####################################################
# ---------------------------------------------------
# Labels, indices and aggregation 
# ---------------------------------------------------
#####################################################

"""
    get_plants_states_syms(
        gens_govs_avrs_states_syms)


Return plants states variables symbols per plant for all generator plants in a vector of vectors.

A generator plant consist of a generator, automatic voltage regulator, governor etc.

"""
function get_plants_states_syms(
    gens_govs_avrs_states_syms)

    return [ :nothing ∈ an_item.gov ?
        [an_item.gen; an_item.avr ] :
        [an_item.gen; an_item.avr; an_item.gov ]
         for an_item in
             gens_govs_avrs_states_syms]


end


"""
    get_generic_state_sym(
        gens_govs_avrs_states_syms,
        gens_nodes_idx;
        label_prefix = "bus")

Return plants states variables symbols or lables for all generator plants  in a flattened
vector.

A generator plant consist of a generator, automatic voltage regulator, governor etc.

"""
function get_generic_state_sym(
    gens_govs_avrs_states_syms,
    gens_nodes_idx;
    label_prefix = "bus")


    plants_states_syms =
        get_plants_states_syms(
            gens_govs_avrs_states_syms)

    # generic_state_sym =
    #     get_labels_by_nodes_idxs_and_vec_vec_syms(
    #         gens_nodes_idx,
    #         plants_states_syms;
    #         label_prefix = "bus" )

    return get_labels_by_nodes_idxs_and_vec_vec_syms(
            gens_nodes_idx,
            plants_states_syms;
            label_prefix = "bus" ) 
    
end


"""
    get_state_labels(
        gens_govs_avrs_states_syms,
        gens_nodes_idx;
        label_prefix = "bus",
        plants_states_by_per_comp = false,
        plants_states_by_per_plant = true )

Return plants states variables symbols or lables for all generator plants  in a flattened
vector.

A generator plant consist of a generator, automatic voltage regulator, governor etc.

"""
function get_state_labels(
    gens_govs_avrs_states_syms,
    gens_nodes_idx;
    label_prefix = "bus",
    plants_states_by_per_comp = false,
    plants_states_by_per_plant = true )


    plants_states_syms =
        get_plants_states_syms(
            gens_govs_avrs_states_syms)

    
    if( plants_states_by_per_comp == true &&
        plants_states_by_per_plant == false)

        plants_states_syms =
            [ :nothing ∈ a_comp.gov ?
            [a_comp.gen; a_comp.avr ] :
            [a_comp.gen; a_comp.avr; a_comp.gov ]
             for a_comp in
                 plants_states_syms ]

        # state_labels =
         return get_labels_by_nodes_idxs_and_vec_vec_syms(
                gens_nodes_idx,
                plants_states_syms;
                label_prefix = label_prefix )
        
    elseif (plants_states_by_per_plant == true &&
        plants_states_by_per_comp == false)
    
        # state_labels =
        return get_labels_by_nodes_idxs_and_vec_vec_syms(
                gens_nodes_idx,
                plants_states_syms;
                label_prefix = label_prefix )

    else
            
        # state_labels =
        return generate_labels_by_nodes_idxs_and_vars(
                gens_nodes_idx,
                plants_states_syms;
                label_prefix = label_prefix )
        
    end
    
    #----------------------------------------
    
end


"""
    get_algebraic_vars_labels(
        dyn_pf_fun_kwd_net_idxs;
        label_prefix = "bus" )

Return plants algebraic variable symbols or lables for all generator plants  in a flattened vector.


# Arguments
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `label_prefix::String="bus"`: the label that should be used as a prefix.


"""
function get_algebraic_vars_labels(
    dyn_pf_fun_kwd_net_idxs;
    label_prefix = "bus" )

    # gens_nodes_with_loc_loads_idx
    
    # :gens_with_loc_loads_idx
    
    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    # algebraic_vars_labels =
    return [generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:vh];
            label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:θh];
             label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:id];
             label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:iq];
             label_prefix = label_prefix ) ]

        
end


"""
    get_network_vars_labels(
        gens_govs_avrs_states_syms,
        dyn_pf_fun_kwd_net_idxs;
        label_prefix = "bus",
        plants_states_by_per_comp = false,
        plants_states_by_per_plant = true )


Return a system network states and algebraic variables label in a flatted.


# Arguments
- `gens_govs_avrs_states_syms`: the namedtuple of states variables symbols  of generator,
   governors and automatic voltage regulator.
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `label_prefix::String="bus"`: the label that should be used as a prefix.


"""
function get_network_vars_labels(
    gens_govs_avrs_states_syms,
    dyn_pf_fun_kwd_net_idxs;
    label_prefix = "bus",
    plants_states_by_per_comp = false,
    plants_states_by_per_plant = true )

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    return [get_state_labels(
        gens_govs_avrs_states_syms,
        gens_nodes_idx);
            
            get_algebraic_vars_labels(
                dyn_pf_fun_kwd_net_idxs;
                label_prefix =
                    label_prefix) ]

    
end

"""
    get_state_vars_idx(
        gens_govs_avrs_states_syms)


Return indices of states  variables in the system in a flatted.


# Arguments
- `gens_govs_avrs_states_syms`: the namedtuple of states variables symbols  of generator,
   governors and automatic voltage regulator.

"""
function get_state_vars_idx(
    gens_govs_avrs_states_syms)

    vec_dims_gens_govs_avrs_states_syms =
        [:nothing ∈ an_item.gov ?
        [length(an_item.gen),
         length(an_item.avr)] :
             [length(an_item.gen),
              length(an_item.avr),
              length(an_item.gov) ]
         for an_item in
             gens_govs_avrs_states_syms]

    dims_gens_govs_avrs_states_syms =
        sum.(vec_dims_gens_govs_avrs_states_syms)
   
    # state_vars_idx =
    #     get_vars_or_paras_Idxs_in_flattend(
    #         dims_gens_govs_avrs_states_syms;
    #         dims_given = true )

    return get_vars_or_paras_Idxs_in_flattend(
            dims_gens_govs_avrs_states_syms;
            dims_given = true )
    
end


"""
    get_plants_states_syms_and_labels(
        gens_govs_avrs_states_syms,
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs)

Return a system namedtuple of `state_vars_idx`, `vec_comp_states_Idx`, `plants_states_syms`,
`generic_state_sym`, `state_labels`, `algebraic_vars_labels`, and `network_vars_labels`.


# Arguments
- `gens_govs_avrs_states_syms`: the namedtuple of states variables symbols  of generator,
   governors and automatic voltage regulator.
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.

"""
function get_plants_states_syms_and_labels(
    gens_govs_avrs_states_syms,
    dyn_pf_fun_kwd_net_idxs,
    dyn_pf_fun_kwd_n2s_idxs)

    #-------------------------------
    
    gens_nodes_idx =
        getproperty(dyn_pf_fun_kwd_net_idxs,
                    :gens_nodes_idx )
    
    #-------------------------------

    vec_vec_gens_govs_avrs_states_syms =
        [ :nothing ∈ an_item.gov ?
        [an_item.gen, an_item.avr ] :
        [an_item.gen, an_item.avr, an_item.gov ]
         for an_item in
             gens_govs_avrs_states_syms]

    vec_dims_gens_govs_avrs_states_syms =
        [:nothing ∈ an_item.gov ?
        [length(an_item.gen),
         length(an_item.avr)] :
             [length(an_item.gen),
              length(an_item.avr),
              length(an_item.gov) ]
         for an_item in
             gens_govs_avrs_states_syms]

    dims_gens_govs_avrs_states_syms =
        sum.(vec_dims_gens_govs_avrs_states_syms)

    state_vars_idx =
        get_vars_or_paras_Idxs_in_flattend(
            dims_gens_govs_avrs_states_syms;
            dims_given = true )

    vec_comp_states_Idx =
        map((x) ->
        get_vars_or_paras_Idxs_in_flattend(
            x;
            dims_given = true ),
            vec_dims_gens_govs_avrs_states_syms )

    plants_states_syms =
        [ :nothing ∈ an_item.gov ?
        [an_item.gen; an_item.avr ] :
        [an_item.gen; an_item.avr; an_item.gov ]
         for an_item in
             gens_govs_avrs_states_syms]

    generic_state_sym =
        get_labels_by_nodes_idxs_and_vec_vec_syms(
            gens_nodes_idx,
            plants_states_syms;
            label_prefix = "bus" ) 
    
    (;state_labels,
     algebraic_vars_labels,
     network_vars_labels) =
        NamedTupleTools.select(
            get_generic_network_vars_labels(
                plants_states_syms,
                dyn_pf_fun_kwd_net_idxs,
                dyn_pf_fun_kwd_n2s_idxs;
                label_prefix =
                    "bus",
                plants_states_by_per_comp =
                    false,
                plants_states_by_per_plant =
                    true,),
            (:state_labels,
             :algebraic_vars_labels,
             :network_vars_labels))
        
    return (;state_vars_idx,
            vec_comp_states_Idx,
            plants_states_syms,
            generic_state_sym,
            state_labels,
            algebraic_vars_labels,
            network_vars_labels)
end


"""
    get_plants_states_syms_wt_labels_wt_names(
        gens_govs_avrs_states_syms,
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs)


Return a system namedtuple of `state_vars_idx`, `vec_comp_states_Idx`, `plants_states_syms`,
`generic_state_sym`, `state_labels`, `algebraic_vars_labels`, `network_vars_labels`,
`model_syms`, `nodes_names`, `gens_nodes_names`, `non_gens_nodes_names`, `SM_gens_nodes_names`,
 and `SC_gens_nodes_names`.


# Arguments
- `gens_govs_avrs_states_syms`: the namedtuple of states variables symbols  of generator,
   governors and automatic voltage regulator.
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.

"""
function get_plants_states_syms_wt_labels_wt_names(
    gens_govs_avrs_states_syms,
    dyn_pf_fun_kwd_net_idxs,
    dyn_pf_fun_kwd_n2s_idxs)

    #-------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx


    (;n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #-------------------------------

    vec_vec_gens_govs_avrs_states_syms =
        [ :nothing ∈ an_item.gov ?
        [an_item.gen, an_item.avr ] :
        [an_item.gen, an_item.avr, an_item.gov ]
         for an_item in
             gens_govs_avrs_states_syms]

    vec_dims_gens_govs_avrs_states_syms =
        [:nothing ∈ an_item.gov ?
        [length(an_item.gen),
         length(an_item.avr)] :
             [length(an_item.gen),
              length(an_item.avr),
              length(an_item.gov) ]
         for an_item in
             gens_govs_avrs_states_syms]

    dims_gens_govs_avrs_states_syms =
        sum.(vec_dims_gens_govs_avrs_states_syms)

    state_vars_idx =
        get_vars_or_paras_Idxs_in_flattend(
            dims_gens_govs_avrs_states_syms;
            dims_given = true )

    vec_comp_states_Idx =
        map((x) ->
        get_vars_or_paras_Idxs_in_flattend(
            x;
            dims_given = true ),
            vec_dims_gens_govs_avrs_states_syms )

    plants_states_syms =
        [ :nothing ∈ an_item.gov ?
        [an_item.gen; an_item.avr ] :
        [an_item.gen; an_item.avr; an_item.gov ]
         for an_item in
             gens_govs_avrs_states_syms]

    generic_state_sym =
        get_labels_by_nodes_idxs_and_vec_vec_syms(
            gens_nodes_idx,
            plants_states_syms;
            label_prefix = "bus" ) 
    
    (;state_labels,
     algebraic_vars_labels,
     network_vars_labels) =
        NamedTupleTools.select(
            get_generic_network_vars_labels(
                plants_states_syms,
                dyn_pf_fun_kwd_net_idxs,
                dyn_pf_fun_kwd_n2s_idxs;
                label_prefix =
                    "bus",
                plants_states_by_per_comp =
                    false,
                plants_states_by_per_plant =
                    true,),
            (:state_labels,
             :algebraic_vars_labels,
             :network_vars_labels))
        
    #----------------------------------------
    # Labels, syms and nodes names
    #----------------------------------------

    model_syms =
        [state_labels;
         algebraic_vars_labels]

    #----------------------------------------
        
    nodes_names =
        ["bus$(n2s_all_nodes_idx[idx])"
         for idx in all_nodes_idx ]

            
    gens_nodes_names = nodes_names[ [n2s_all_nodes_idx[idx]
                     for idx in gens_nodes_idx] ]
            
    non_gens_nodes_names =
        nodes_names[ [n2s_all_nodes_idx[idx]
                     for idx in non_gens_nodes_idx] ]
    
            
    # gens_nodes_names =
    #     nodes_names[ gens_nodes_idx ]

            
    # non_gens_nodes_names =
    #     nodes_names[ non_gens_nodes_idx ]
    
    #-------------------------------
    # I need to get the indices of SM and SC gens

    SM_gens_idx = [ idx
         for (idx, an_item) in
                    enumerate(gens_govs_avrs_states_syms)
                    if  :nothing ∉ an_item.gov ]


    SC_gens_idx = [ idx
         for (idx, an_item) in
                    enumerate(gens_govs_avrs_states_syms)
                    if  :nothing ∈ an_item.gov ]

    SM_gens_nodes_names =
        gens_nodes_names[ SM_gens_idx ]

    SC_gens_nodes_names =
        gens_nodes_names[ SC_gens_idx]



    
    return (;state_vars_idx,
            vec_comp_states_Idx,
            plants_states_syms,
            generic_state_sym,
            state_labels,
            algebraic_vars_labels,
            network_vars_labels,
            
            model_syms,
            nodes_names,
            gens_nodes_names,
            non_gens_nodes_names,
            SM_gens_nodes_names,
            SC_gens_nodes_names )
end


"""
    get_model_syms(
        state_labels,
        dyn_pf_fun_kwd_net_idxs;
        label_prefix = "bus")


Return a system model symbols or labels for state and algebraic variables in a flattened vector.

# Arguments
- `state_labels`: the system state variables labels.
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `label_prefix::String="bus"`: the label that should be used as a prefix.

"""
function get_model_syms(
    state_labels,
    dyn_pf_fun_kwd_net_idxs;
    label_prefix = "bus")

    #----------------------------------------
    # Labels, syms and nodes names
    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    
    algebraic_vars_labels =
        [generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:vh];
            label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:θh];
             label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:id];
             label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:iq];
             label_prefix = label_prefix) ]

    #----------------------------------------
    
    return [state_labels;
         algebraic_vars_labels]
    

end


"""
    get_dynamic_comps_init_out_dyn_callback_funcs(
        gens_govs_avrs_types)

Returns namedtuples of dynamic components functions `comps_callback_paras_funs`, `comps_init_funs`, `comps_output_funs`, `ode_comps_dyn_funs`, `dae_comps_dyn_funs`, `comps_dyn_funs`.


# Arguments
- ``gens_govs_avrs_states_syms`: the namedtuple of states variables symbols of generators, governors and automatic voltage regulators. 
"""
function get_dynamic_comps_init_out_dyn_callback_funcs(
    gens_govs_avrs_types)

    list_gens_type_syms =
        Symbol.([a_type.gen for a_type in
             gens_govs_avrs_types])


    list_govs_type_syms =
        Symbol.([a_type.gov for a_type in
             gens_govs_avrs_types])


    list_avrs_type_syms =
        Symbol.([a_type.avr for a_type in
             gens_govs_avrs_types])

    #----------------------------------------

    gens_cb_paras_func =
        get_gens_callback_paras_func(
         list_gens_type_syms)

    gens_init_func =
        get_gens_init_func(
         list_gens_type_syms)
    
    gens_output_func =
        get_gens_output_func(
            list_gens_type_syms)
    
    ode_gens_dyn_func =
        get_gens_dyn_func(
            list_gens_type_syms,
            :ode)

    
    dae_gens_dyn_func =
        get_gens_dyn_func(
            list_gens_type_syms,
            :dae)

    
    gens_dyn_func =
        get_gens_dyn_func(
            list_gens_type_syms)
    
    
    govs_cb_paras_func =
        get_govs_callback_paras_func(
            list_govs_type_syms)
    
    govs_init_func =
        get_govs_init_func(
         list_govs_type_syms)
    
    govs_output_func =
        get_govs_output_func(
            list_govs_type_syms)
    
    ode_govs_dyn_func =
        get_govs_dyn_func(
            list_govs_type_syms,
            :ode)

    
    dae_govs_dyn_func =
        get_govs_dyn_func(
            list_govs_type_syms,
            :dae)

    
    govs_dyn_func =
        get_govs_dyn_func(
            list_govs_type_syms)
    
    
    ##

    
    avrs_cb_paras_func =
        get_avrs_callback_paras_func(
            list_avrs_type_syms)
    
    avrs_init_func =
        get_avrs_init_func(
         list_avrs_type_syms)
    
    avrs_output_func =
        get_avrs_output_func(
            list_avrs_type_syms)

    
    ode_avrs_dyn_func =
        get_avrs_dyn_func(
            list_avrs_type_syms,
            :ode)

    
    dae_avrs_dyn_func =
        get_avrs_dyn_func(
            list_avrs_type_syms,
            :dae)
    
    avrs_dyn_func =
        get_avrs_dyn_func(
            list_avrs_type_syms)

    
    #----------------------------------------

    return (comps_callback_paras_funs = [
        (gen_cb_paras_fun = a_gen_fun,
         avr_cb_paras_fun = a_avr_fun,
         gov_cb_paras_fun = a_gov_fun) 
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( gens_cb_paras_func,
                    avrs_cb_paras_func,
                    govs_cb_paras_func)],

    comps_init_funs = [
        (gen_init_fun = a_gen_fun,
         avr_init_fun = a_avr_fun,
         gov_init_fun = a_gov_fun)
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( gens_init_func,
                    avrs_init_func,
                    govs_init_func)],

    
    comps_output_funs = [
        (gen_output_fun = a_gen_fun,
         avr_output_fun = a_avr_fun,
         gov_output_fun = a_gov_fun) 
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( gens_output_func,
                    avrs_output_func,
                    govs_output_func)],


    
    ode_comps_dyn_funs = [
        (gen_dyn_fun = a_gen_fun,
         avr_dyn_fun = a_avr_fun,
         gov_dyn_fun = a_gov_fun)
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( ode_gens_dyn_func,
                    ode_avrs_dyn_func,
                    ode_govs_dyn_func)],

    
    dae_comps_dyn_funs = [
        (gen_dyn_fun = a_gen_fun,
         avr_dyn_fun = a_avr_fun,
         gov_dyn_fun = a_gov_fun)
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( dae_gens_dyn_func,
                    dae_avrs_dyn_func,
                    dae_govs_dyn_func)],
    
    comps_dyn_funs = [
        (gen_dyn_fun = a_gen_fun,
         avr_dyn_fun = a_avr_fun,
         gov_dyn_fun = a_gov_fun)
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( gens_dyn_func,
                    avrs_dyn_func,
                    govs_dyn_func)])

end


"""
    get_state_and_algebraic_vars_Idx_in_state(
        state_labels,
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs;
        no_lines_fault = 1 )


Return indices of state and algebraic variables for a normal system and faulted system states in a flattened vector.


 `state_vars_and_i_dq_Idx_in_state`, `state_vars_and_i_dq_wt_fault_Idx_in_state`, `state_algebraic_vars_Idx_in_state`, and `state_algebraic_vars_wt_fault_Idx_in_state`.

# Arguments
- `state_labels`: the system state variables labels.
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.
- `no_lines_fault::Int64 = 1`: the number of lines faults.

"""
function get_state_and_algebraic_vars_Idx_in_state(
    state_labels,
    dyn_pf_fun_kwd_net_idxs,
    dyn_pf_fun_kwd_n2s_idxs;
    no_lines_fault = 1 )

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

     (;n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_n2s_idxs,
             (:n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx))
    
    state_vars_and_i_dq_Idx_in_state =
        get_state_vars_and_i_dq_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx )

    state_vars_and_i_dq_wt_fault_Idx_in_state =
        get_state_vars_and_i_dq_wt_fault_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault =
                no_lines_fault)

    state_algebraic_vars_Idx_in_state =
        get_state_algebraic_vars_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx )

    state_algebraic_vars_wt_fault_Idx_in_state =
        get_state_algebraic_vars_wt_fault_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault =
                no_lines_fault)

    return (;state_vars_and_i_dq_Idx_in_state,
            state_vars_and_i_dq_wt_fault_Idx_in_state,
            state_algebraic_vars_Idx_in_state,
            state_algebraic_vars_wt_fault_Idx_in_state)
    
end


"""
    get_model_nodes_types_names(
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs )


Returns nodes labels for generator buses, non-generator buses and all buses.

# Arguments
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.

"""
function get_model_nodes_types_names(
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs )


    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx


     (;n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_n2s_idxs,
             (:n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx))
    
    nodes_names =
        ["bus$(n2s_all_nodes_idx[idx])"
         for idx in all_nodes_idx ]
            
    gens_nodes_names = nodes_names[ [n2s_all_nodes_idx[idx]
                     for idx in gens_nodes_idx] ]
            
    non_gens_nodes_names =
        nodes_names[ [n2s_all_nodes_idx[idx]
                     for idx in non_gens_nodes_idx] ]

    return (;nodes_names,
            gens_nodes_names,
            non_gens_nodes_names)

end


function get_model_states_comp_idxs_in_Idx(
        network_vars_labels,
        all_nodes_idx,
        n2s_all_nodes_idx;
    vars =
        [:δ, :ω, :ed_dash, :eq_dash] )
    
    nodes_names =
        ["bus$(n2s_all_nodes_idx[idx])"
         for idx in all_nodes_idx ]

    # δ_ω_ed_dash_eq_dash_indices =
    #     get_nodes_state_algb_vars_indices_in_system(
    #         ;network_vars_labels =
    #             network_vars_labels,
    #         nodes_name = nodes_names,
    #         vars = vars )
      
    # δ_idx_in_state =
    #     [idx for idx in
    #          first.(δ_ω_ed_dash_eq_dash_indices)]
    
    # ω_idx_in_state =
    #     [idx for idx in
    #          second.(δ_ω_ed_dash_eq_dash_indices)]
    
    # ed_dash_idx_in_state =
    #     [idx for idx in
    #          third.(δ_ω_ed_dash_eq_dash_indices)]
    
    # eq_dash_idx_in_state =
    #     [idx for idx in
    #          fourth.(δ_ω_ed_dash_eq_dash_indices)]
    
    selected_states_indices =
        get_nodes_state_algb_vars_indices_in_system(
            ;network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_names,
            vars = vars )

    # @show length(selected_states_indices)    
    # @show length(selected_states_indices[1])
    
    dict_nth_funs =
        get_dict_first_to_tenth_funs(
        length(selected_states_indices[1]) )
    
    nt_names = Symbol.([ "$(a_var)_idx_in_state"
                         for a_var in vars ])

    vec_vec_state_idx_in_state =
        [ [idx for idx in (dict_nth_funs[a_var_idx]).(
               selected_states_indices) ]
          for a_var_idx in 1:length(vars) ]

    return namedtuple(
        OrderedDict( a_sym => vec_idx
             for (a_sym, vec_idx) in
                     zip(nt_names,
                         vec_vec_state_idx_in_state) ))


end


"""
    get_mass_matrix_and_bool_dae_vars(
        state_labels,
        algebraic_vars_labels)


Returnsfor a system and generator mass matrices and dae boolean variables `model_mass_matrix`, `model_bool_dae_vars`, `ode_gens_mass_matrix`, `ode_gens_bool_dae_vars`.

"""
function get_mass_matrix_and_bool_dae_vars(
    state_labels,
    algebraic_vars_labels)

    
    model_mass_matrix =
        DAE_MassMatrix(
            length(state_labels),
            length(algebraic_vars_labels) )
    
    model_bool_dae_vars =
        DAE_BoolVector(
            length(state_labels),
            length(algebraic_vars_labels) )
        
    ode_gens_mass_matrix =
        DAE_MassMatrix(
            length(state_labels),
            0 )
    
    ode_gens_bool_dae_vars =
        DAE_BoolVector(
            length(state_labels),
            0 )

    return (;
            model_mass_matrix,
            model_bool_dae_vars,
            ode_gens_mass_matrix,
            ode_gens_bool_dae_vars)

end


#---------------------------------------------------
# indexing of states, parameters and functions
#---------------------------------------------------


function get_generic_red_vh_θh_wt_slack_value_Idx(
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)


     (n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_n2s_idxs,
             (:n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx))
    
    (;slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))

    #--------------------------------------

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_load_idx
    
    dyn_pf_flat_vh_flat_θh_Idx =
        get_generic_flat_vh_flat_θh_Idx(
            all_nodes_idx)
    

    (;dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    
    #--------------------------------------

    non_gens_vh_idx =
         dyn_pf_vh_Idxs[
             non_gens_nodes_idx]

    non_slack_gens_θh_idx =
        dyn_pf_θh_Idxs[
            non_slack_gens_nodes_idx]

    non_gens_θh_idx =
        dyn_pf_θh_Idxs[
            non_gens_nodes_idx]

    red_vh_Idxs = non_gens_vh_idx

    red_θh_Idxs =
        setdiff(dyn_pf_θh_Idxs,
                dyn_pf_θh_Idxs[slack_gens_nodes_idx])

    non_slack_gens_and_non_gens_idx =
        sort([non_slack_gens_nodes_idx;
              non_gens_nodes_idx])

    nodes_with_demands_idx =
        convert(Vector{Int64},
                sort([non_gens_nodes_idx;
              gens_with_loc_load_idx]))

    #----------------------------------------    
    
    red_vh_θh_idx =
        [  non_gens_vh_idx...;
           non_slack_gens_θh_idx...;
           non_gens_θh_idx... ]

    #----------------------------------------

    red_dict_θh_idx2Idx =
        OrderedDict{Int64, Int64}(
            idx => Idx for (idx, Idx) in
                zip(non_slack_gens_and_non_gens_idx,
                     red_θh_Idxs ) )

    red_dict_θh_idx2Idx_in_Idx =
       OrderedDict{Int64, Int64}(
            idx => Idx for (idx, Idx) in
                zip( non_slack_gens_and_non_gens_idx,
                     1:length(red_θh_Idxs) ) )

    red_non_slack_gens_θh_idx2Idx =
        [ red_dict_θh_idx2Idx[idx]
         for idx in
             non_slack_gens_nodes_idx ]

    red_non_slack_gens_θh_idx2Idx_in_Idx =
        [ red_dict_θh_idx2Idx_in_Idx[idx]
          for idx in
             non_slack_gens_nodes_idx ]

    red_non_gens_θh_idx2Idx =
        [ red_dict_θh_idx2Idx[idx]
         for idx in
             non_gens_nodes_idx ]

    red_non_gens_θh_idx2Idx_in_Idx =
        [ red_dict_θh_idx2Idx_in_Idx[idx]
          for idx in
              non_gens_nodes_idx ]
        
    #----------------------------------------        
    
    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_labels =
            [:red_non_gens_vh_Idxs,
             :red_non_slack_gens_θh_Idxs,       
             :red_non_gens_θh_Idxs,
             :red_slack_value_Idxs]

    dim_non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value =
        [ length(non_gens_nodes_idx),
          length(non_slack_gens_nodes_idx),
          length(non_gens_nodes_idx),
          1]

    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value;
            dims_given = true )

    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_Idxs =
        namedtuple(
            OrderedDict{Symbol,UnitRange{Int64}}(
                a_sym => idx_range
                for (a_sym, idx_range) in zip(
                    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_labels,

                    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_Idx )))

    (;red_non_gens_vh_Idxs,
     red_non_slack_gens_θh_Idxs,
     red_non_gens_θh_Idxs,
     red_slack_value_Idxs, ) =
        NamedTupleTools.select(
            non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_Idxs,
            (:red_non_gens_vh_Idxs,
             :red_non_slack_gens_θh_Idxs,       
             :red_non_gens_θh_Idxs,
             :red_slack_value_Idxs, ))
    
    #--------------------------------------------

    return (;dyn_pf_flat_vh_flat_θh_Idx,
            
            non_gens_vh_idx,
            non_slack_gens_θh_idx,
            non_gens_θh_idx,
            
            red_vh_Idxs,
            red_θh_Idxs,
            non_slack_gens_and_non_gens_idx,
            
            red_vh_θh_idx,
            red_dict_θh_idx2Idx,
            red_dict_θh_idx2Idx_in_Idx,
            red_non_slack_gens_θh_idx2Idx,
            red_non_gens_θh_idx2Idx,
            red_non_gens_θh_idx2Idx_in_Idx,

            red_non_gens_vh_Idxs,
            red_non_slack_gens_θh_Idxs,
            red_non_gens_θh_Idxs,

            red_slack_value_Idxs,

            slack_gens_nodes_idx,
            non_slack_gens_nodes_idx,
            gens_nodes_idx,
            non_gens_nodes_idx,
            
            gens_with_loc_load_idx,
            gens_nodes_with_loc_loads_idx,
            nodes_with_demands_idx,
            
            all_nodes_idx,

            n2s_slack_gens_idx,
            n2s_non_slack_gens_idx,
            n2s_gens_idx,
            n2s_non_gens_idx,
            n2s_gens_with_loc_load_idxs,
            n2s_all_nodes_idx,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )
        
end

"""
    get_slack_gens_vh_θh_gens_vh_non_slack_gens_vh(
        vh,
        θh,
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs)


Returns `slack_gens_vh`, `slack_gens_θh`, `gens_vh`, and `non_slack_gens_vh`.

# Arguments
- `vh`: the network nodes voltage magnitudes.
- `θh`: the network nodes voltage angles.
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.

"""
function get_slack_gens_vh_θh_gens_vh_non_slack_gens_vh(
    vh,
    θh,
    dyn_pf_fun_kwd_net_idxs,
    dyn_pf_fun_kwd_n2s_idxs)

    
    (slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_with_loc_load_idx,
    all_nodes_idx) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))    

    (n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))


    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_slack_gens_nodes_idx ]
    
    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]
    
    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]


    transformed_non_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_slack_gens_nodes_idx ]
    
    # non_slack_gens_idx =
    #     setdiff(gens_nodes_idx,
    #             non_slack_gens_nodes_idx)

    slack_gens_vh =
        vh[transformed_slack_gens_nodes_idx]
    
    slack_gens_θh =
        θh[transformed_slack_gens_nodes_idx]
    
    gens_vh =
        vh[transformed_gens_nodes_idx]
    
    non_slack_gens_vh =
        vh[transformed_non_slack_gens_nodes_idx ]
    
    return (;slack_gens_vh,
            slack_gens_θh,

            gens_vh,
            non_slack_gens_vh )
    
end


"""
    get_pf_vh_θh_idx_and_idx2Idx(
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs)


Returns indices of `vh`, `θh`, and indices transformation dictionaries.

# Arguments
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.

"""
function get_pf_vh_θh_idx_and_idx2Idx(
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)


     (n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_n2s_idxs,
             (:n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx))
    
    (slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :gens_with_loc_load_idx,
             :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx
    
    #-------------------------------

    non_slack_gens_and_non_gens_idx =
        sort([non_slack_gens_nodes_idx;
              non_gens_nodes_idx])

    nodes_with_demands_idx =
        convert(Vector{Int64},
                sort([non_gens_nodes_idx;
              gens_with_loc_load_idx]))

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
        
    transformed_non_slack_gens_and_non_gens_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_slack_gens_and_non_gens_idx]
    
    transformed_nodes_with_demands_idx = [
        n2s_all_nodes_idx[idx]
        for idx in nodes_with_demands_idx ]
    
    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]

    #------------------------------------------        
    
    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_load_idx
    
    dyn_pf_flat_vh_flat_θh_Idx =
        get_generic_flat_vh_flat_θh_Idx(
            all_nodes_idx)
    

    (;dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    #--------------------------------------

    dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx =
        get_generic_flat_vh_flat_θh_wt_slack_value_Idx(
            all_nodes_idx)

    (;dyn_slack_value_Idxs,) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
            (:dyn_slack_value_Idxs, ))
    
        
    #--------------------------------------
        
    #  non_gens_vh_idx =
    #      dyn_pf_vh_Idxs[
    #          non_gens_nodes_idx]

    # non_slack_gens_θh_idx =
    #     dyn_pf_θh_Idxs[
    #         non_slack_gens_nodes_idx]

    # non_gens_θh_idx =
    #     dyn_pf_θh_Idxs[
    #         non_gens_nodes_idx]

    # red_vh_Idxs = non_gens_vh_idx

    # red_θh_Idxs =
    #     setdiff(dyn_pf_θh_Idxs,
    #             dyn_pf_θh_Idxs[slack_gens_nodes_idx])
    
     non_gens_vh_idx =
         dyn_pf_vh_Idxs[
             transformed_non_gens_nodes_idx]

    non_slack_gens_θh_idx =
        dyn_pf_θh_Idxs[
            transformed_non_slack_gens_nodes_idx]

    non_gens_θh_idx =
        dyn_pf_θh_Idxs[
            transformed_non_gens_nodes_idx]

    red_vh_Idxs = non_gens_vh_idx

    red_θh_Idxs =
        setdiff(dyn_pf_θh_Idxs,
                dyn_pf_θh_Idxs[
                    transformed_slack_gens_nodes_idx])
    
    #----------------------------------------    
    
    red_vh_θh_idx =
        [  non_gens_vh_idx...;
           non_slack_gens_θh_idx...;
           non_gens_θh_idx... ]

    #----------------------------------------

    red_dict_θh_idx2Idx =
        OrderedDict{Int64, Int64}(
            idx => Idx for (idx, Idx) in
                zip(non_slack_gens_and_non_gens_idx,
                     red_θh_Idxs ) )

    red_dict_θh_idx2Idx_in_Idx =
       OrderedDict{Int64, Int64}(
            idx => Idx for (idx, Idx) in
                zip( non_slack_gens_and_non_gens_idx,
                     1:length(red_θh_Idxs) ) )

    red_non_slack_gens_θh_idx2Idx =
        [ red_dict_θh_idx2Idx[idx]
         for idx in
             non_slack_gens_nodes_idx ]

    red_non_slack_gens_θh_idx2Idx_in_Idx =
        [ red_dict_θh_idx2Idx_in_Idx[idx]
          for idx in
             non_slack_gens_nodes_idx ]

    red_non_gens_θh_idx2Idx =
        [ red_dict_θh_idx2Idx[idx]
         for idx in
             non_gens_nodes_idx ]

    red_non_gens_θh_idx2Idx_in_Idx =
        [ red_dict_θh_idx2Idx_in_Idx[idx]
          for idx in
              non_gens_nodes_idx ]
    
    #-------------------------------

    non_gens_vh_non_slack_gens_θh_non_gens_θh_labels =
            [:red_non_gens_vh_Idxs,
             :red_non_slack_gens_θh_Idxs,       
             :red_non_gens_θh_Idxs]

    dim_non_gens_vh_non_slack_gens_θh_non_gens_θh =
        [ length(non_gens_nodes_idx),
          length(non_slack_gens_nodes_idx),
          length(non_gens_nodes_idx)]

    non_gens_vh_non_slack_gens_θh_non_gens_θh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_non_gens_vh_non_slack_gens_θh_non_gens_θh;
            dims_given = true )

    non_gens_vh_non_slack_gens_θh_non_gens_θh_Idxs =
        namedtuple(
            OrderedDict{Symbol,UnitRange{Int64}}(
                a_sym => idx_range
                for (a_sym, idx_range) in zip(
                    non_gens_vh_non_slack_gens_θh_non_gens_θh_labels,

                    non_gens_vh_non_slack_gens_θh_non_gens_θh_Idx )))

    (;red_non_gens_vh_Idxs,
     red_non_slack_gens_θh_Idxs,
     red_non_gens_θh_Idxs ) =
        NamedTupleTools.select(
            non_gens_vh_non_slack_gens_θh_non_gens_θh_Idxs,
            (:red_non_gens_vh_Idxs,
             :red_non_slack_gens_θh_Idxs,
             :red_non_gens_θh_Idxs))
    

    #--------------------------------------
    
    non_gens_vh_gens_θh_non_gens_θh_slack_value_labels =
            [:non_gens_vh_Idxs,
             :gens_θh_Idxs,       
             :non_gens_θh_Idxs,
             :slack_value_Idxs]


    dim_non_gens_vh_gens_θh_non_gens_θh_slack_value =
        [ length(non_gens_nodes_idx),
          length(gens_nodes_idx),
          length(non_gens_nodes_idx),
          1]

    non_gens_vh_gens_θh_non_gens_θh_slack_value_idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_non_gens_vh_gens_θh_non_gens_θh_slack_value;
            dims_given = true )

    non_gens_vh_gens_θh_non_gens_θh_slack_value_Idxs =
        namedtuple(
            OrderedDict{Symbol,UnitRange{Int64}}(
                a_sym => idx_range
                for (a_sym, idx_range) in zip(
                    non_gens_vh_gens_θh_non_gens_θh_slack_value_labels,
                    non_gens_vh_gens_θh_non_gens_θh_slack_value_idx )))
    
    #--------------------------------------
    
    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_labels =
            [:red_non_gens_vh_Idxs,
             :red_non_slack_gens_θh_Idxs,       
             :red_non_gens_θh_Idxs,
             :red_slack_value_Idxs]

    dim_non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value =
        [ length(non_gens_nodes_idx),
          length(non_slack_gens_nodes_idx),
          length(non_gens_nodes_idx),
          1]

    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value;
            dims_given = true )

    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_Idxs =
        namedtuple(
            OrderedDict{Symbol,UnitRange{Int64}}(
                a_sym => idx_range
                for (a_sym, idx_range) in zip(
                    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_labels,

                    non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_Idx )))

    (;red_slack_value_Idxs, ) =
        NamedTupleTools.select(
            non_gens_vh_non_slack_gens_θh_non_gens_θh_slack_value_Idxs,
            (:red_slack_value_Idxs, ))
    
    #--------------------------------------------
    

    return (;dyn_pf_flat_vh_flat_θh_Idx,
            
            dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
            dyn_slack_value_Idxs,

            non_gens_vh_idx,
            non_slack_gens_θh_idx,
            non_gens_θh_idx,
            
            red_vh_Idxs,
            red_θh_Idxs,
            non_slack_gens_and_non_gens_idx,
            
            red_vh_θh_idx,
            red_dict_θh_idx2Idx,
            red_dict_θh_idx2Idx_in_Idx,
            red_non_slack_gens_θh_idx2Idx,
            red_non_gens_θh_idx2Idx,
            red_non_gens_θh_idx2Idx_in_Idx,

            red_non_gens_vh_Idxs,
            red_non_slack_gens_θh_Idxs,
            red_non_gens_θh_Idxs,

            red_slack_value_Idxs,

            slack_gens_nodes_idx,
            non_slack_gens_nodes_idx,
            gens_nodes_idx,
            non_gens_nodes_idx,
            
            gens_with_loc_load_idx,
            gens_nodes_with_loc_loads_idx,
            nodes_with_demands_idx,
            
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
            
            transformed_non_slack_gens_and_non_gens_idx,
            transformed_nodes_with_demands_idx,
            
            transformed_all_nodes_idx,


            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs,
            
            non_gens_vh_gens_θh_non_gens_θh_slack_value_Idxs)
    
end


"""
    get_labels_by_nodes_idxs_and_vec_vec_syms(
        nodes_idxs,
        plants_states_syms;
        label_prefix = "bus" )


Returns labels for state variables for a list of nodes in `nodes_idxs`.

# Arguments
- `nodes_idxs`: the list of nodes indices.
- `plants_states_syms`: the state variables symbols per plant for all plants.

"""
function get_labels_by_nodes_idxs_and_vec_vec_syms(
    nodes_idxs,
    plants_states_syms;
    label_prefix = "bus" )

    nodes_labels = Symbol[]


    for (node_idx, states_sym) in
        zip(nodes_idxs, plants_states_syms)
        
        append!(nodes_labels,
                [Symbol(label_prefix, node_idx, "_",
                        sym )
                 for sym in states_sym])
    end

    return nodes_labels
    
end

#-----------------------------------------------------



"""
    generate_labels_by_nodes_idxs_and_vars(
        nodes_idxs,
        nodes_vars_syms;
        label_prefix = "bus" )

Return labels for nodes variables.

# Example

t_gens_nodes_idx = [1, 2, 3]

t_state_vars_syms_internal_mode = [:δ, :ω]

generate_labels_by_nodes_idxs_and_vars(
    t_gens_nodes_idx,
    t_state_vars_syms_internal_mode)

"""
function generate_labels_by_nodes_idxs_and_vars(
    nodes_idxs,
    nodes_vars_syms;
    label_prefix = "bus" )


    nodes_labels = Symbol[]


    for node_idx in  nodes_idxs
        
        append!(nodes_labels,
                 [Symbol(label_prefix, node_idx, "_", sym )
                  for sym in nodes_vars_syms])
    end

    return nodes_labels
    
end


"""
    get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
        list_state_vars_syms,
        gens_nodes_idx )


Returns Idxs of variables per node based on a list of variables.

# Example

t_gens_nodes_idx = [1, 2, 3]

t_list_state_vars_syms = [:δ, :ω, :eq_dash, :E_fd]

 get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
    t_list_state_vars_syms,
    t_gens_nodes_idx )

>>> Vector{UnitRange{Int64}}:
 1:4
 5:8
 9:12

"""
function get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
    list_state_vars_syms,
    gens_nodes_idx )

    dims_state_vars =
        length.([list_state_vars_syms
                 for idx in gens_nodes_idx])

    return  get_vars_or_paras_Idxs_in_flattend(
        dims_state_vars;
        dims_given = true )

end


"""
    get_states_idx_by_nodes_idx_wt_vars_syms(
        vec_nodes_states_vars_syms )


Returns indices of state variables provided in `vec_nodes_states_vars_syms` in in the network. 
"""
function get_states_idx_by_nodes_idx_wt_vars_syms(
    vec_nodes_states_vars_syms )

    dims_state_vars = length.(
        vec_nodes_states_vars_syms)

    return  get_vars_or_paras_Idxs_in_flattend(
            dims_state_vars;
            dims_given = true )

end

#------------------------------------------
# From sd-dynamics-dae-generic-formulation
#-------------------------------------------


"""
    get_generic_flat_vh_flat_θh_Idx(
        gens_nodes_idx,
        all_nodes_idx)


Returns a flattend vector of indices of nodes voltages magnitude and nodes voltages angle.

"""
function get_generic_flat_vh_flat_θh_Idx(
    gens_nodes_idx,
    all_nodes_idx)

    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    flat_vh_flat_θh_idx_labels =
        [:dyn_pf_vh_Idxs,
         :dyn_pf_θh_Idxs]

    #----------------------------------------
    
    dim_flat_vh_flat_θh =
        [ dim_all_nodes,
          dim_all_nodes]
    
    #--------------------------------------
    
    flat_vh_flat_θh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_flat_vh_flat_θh;
            dims_given = true )

    return namedtuple(
        OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(flat_vh_flat_θh_idx_labels,
                flat_vh_flat_θh_Idx )))

end



"""
    get_generic_flat_vh_flat_θh_Idx(
        all_nodes_idx)


Returns a flattend vector of indices of nodes voltages magnitude and nodes voltages angle.

"""
function get_generic_flat_vh_flat_θh_Idx(
    all_nodes_idx)

    #----------------------------------------

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    flat_vh_flat_θh_idx_labels =
        [:dyn_pf_vh_Idxs,
         :dyn_pf_θh_Idxs]

    #----------------------------------------
    
    dim_flat_vh_flat_θh =
        [ dim_all_nodes,
          dim_all_nodes]
    
    #--------------------------------------
    
    flat_vh_flat_θh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_flat_vh_flat_θh;
            dims_given = true )

    return namedtuple(
        OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(flat_vh_flat_θh_idx_labels,
                flat_vh_flat_θh_Idx )))

end


"""
    get_generic_flat_vh_flat_θh_id_iq_Idx(
        gens_nodes_idx,
        all_nodes_idx)


Returns a flattend vector of indices of nodes voltages magnitude `vh`, nodes voltages angle `θh`, generators direct-axis currents `id`, generators quadrature axix currents `iq`.

"""
function get_generic_flat_vh_flat_θh_id_iq_Idx(
    gens_nodes_idx,
    all_nodes_idx)

    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    flat_vh_flat_θh_id_iq_idx_labels =
        [:dyn_pf_vh_Idxs,
         :dyn_pf_θh_Idxs,       
         :dyn_pf_id_Idxs,
         :dyn_pf_iq_Idxs]

    #----------------------------------------
    
    dim_flat_vh_flat_θh_id_iq =
        [ dim_all_nodes,
          dim_all_nodes,
          dim_gens,
          dim_gens]
    
    #--------------------------------------
    
    flat_vh_flat_θh_id_iq_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_flat_vh_flat_θh_id_iq;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(flat_vh_flat_θh_id_iq_idx_labels,
                flat_vh_flat_θh_id_iq_Idx )))

end


"""
    get_generic_flat_vh_flat_θh_wt_slack_value_Idx(
        all_nodes_idx)

Returns a flattend vector of indices of nodes voltages magnitude `vh`, nodes voltages angle `θh`, and `slack` variable`. 
"""
function get_generic_flat_vh_flat_θh_wt_slack_value_Idx(
    all_nodes_idx)

    #----------------------------------------

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    flat_vh_flat_θh_slack_value_idx_labels =
        [:dyn_pf_vh_Idxs,
         :dyn_pf_θh_Idxs,       
         :dyn_slack_value_Idxs]

    #----------------------------------------
    
    dim_flat_vh_flat_θh_slack_value =
        [ dim_all_nodes,
          dim_all_nodes,
          1]
    
    #--------------------------------------
    
    flat_vh_flat_θh_slack_value_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_flat_vh_flat_θh_slack_value;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(flat_vh_flat_θh_slack_value_idx_labels,
                flat_vh_flat_θh_slack_value_Idx )))

end


"""
    get_generic_vh_θh_id_iq_vhf_θhf_Idx(
        gens_nodes_idx,
        all_nodes_idx;
        no_lines_fault = 1)


Returns a flattend vector of indices of nodes voltages magnitude `vh`, nodes voltages angle `θh`, generators direct-axis currents `id`, generators quadrature axix currents `iq`, fault node voltages magnitude `vhf`, fault nodes voltage angle `θhf`.
"""
function get_generic_vh_θh_id_iq_vhf_θhf_Idx(
    gens_nodes_idx,
    all_nodes_idx;
    no_lines_fault = 1)

    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    vh_θh_id_iq_vhf_θhf_labels =
        [:dyn_pf_vh_Idxs,
         :dyn_pf_θh_Idxs,       
         :dyn_pf_id_Idxs,
         :dyn_pf_iq_Idxs,
         :dyn_pf_vhf_Idxs,
         :dyn_pf_θhf_Idxs]

    #----------------------------------------
    
    dim_vh_θh_id_iq_vhf_θhf =
        [ dim_all_nodes,
          dim_all_nodes,
          dim_gens,
          dim_gens,
          no_lines_fault,
          no_lines_fault]
    
    #--------------------------------------
    
    vh_θh_id_iq_vhf_θhf_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_vh_θh_id_iq_vhf_θhf;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(vh_θh_id_iq_vhf_θhf_labels,
                vh_θh_id_iq_vhf_θhf_Idx )))

end


"""
    get_generic_vh_vhf_θh_θhf_id_iq_Idx(
        gens_nodes_idx,
        all_nodes_idx;
        no_lines_fault = 1)


Returns a flattend vector of indices of nodes voltages magnitude `vh`, nodes voltages angle `θh`, generators direct-axis currents `id`, generators quadrature axix currents `iq`, fault node voltages magnitude `vhf`, fault nodes voltage angle `θhf`.
"""
function get_generic_vh_vhf_θh_θhf_id_iq_Idx(
    gens_nodes_idx,
    all_nodes_idx;
    no_lines_fault = 1)

    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    vh_vhf_θh_θhf_id_iq_labels =
        [:dyn_pf_vh_Idxs,
         :dyn_pf_vhf_Idxs,
         :dyn_pf_θh_Idxs,
         :dyn_pf_θhf_Idxs,
         :dyn_pf_id_Idxs,
         :dyn_pf_iq_Idxs]

    #----------------------------------------
    
    dim_vh_vhf_θh_θhf_id_iq =
        [ dim_all_nodes,
          no_lines_fault,
          dim_all_nodes,
          no_lines_fault,
          dim_gens,
          dim_gens ]
    
    #--------------------------------------
    
    vh_vhf_θh_θhf_id_iq_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_vh_vhf_θh_θhf_id_iq;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(vh_vhf_θh_θhf_id_iq_labels,
                vh_vhf_θh_θhf_id_iq_Idx )))

end


"""
    get_generic_vh_vhf_Idx(
        all_nodes_idx;
        no_lines_fault = 1)


Returns a flattend vector of indices of nodes voltages magnitude `vh`, and fault node voltages magnitude `vhf`.
"""
function get_generic_vh_vhf_Idx(
    all_nodes_idx;
    no_lines_fault = 1)

    #----------------------------------------

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    vh_vhf_labels =
        [:dyn_pf_vh_Idxs,
         :dyn_pf_vhf_Idxs]

    #----------------------------------------
    
    dim_vh_vhf =
        [ dim_all_nodes,
          no_lines_fault]
    
    #--------------------------------------
    
    vh_vhf_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_vh_vhf;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(vh_vhf_labels,
                vh_vhf_Idx )))

end


"""
    get_generic_θh_θhf_Idx(
        all_nodes_idx;
        no_lines_fault = 1)


Returns a flattend vector of indices of nodes voltages angle `θh`, fault nodes voltage angle `θhf`.
"""
function get_generic_θh_θhf_Idx(
    all_nodes_idx;
    no_lines_fault = 1)

    #----------------------------------------

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    θh_θhf_labels =
        [:dyn_pf_θh_Idxs,       
         :dyn_pf_θhf_Idxs]

    #----------------------------------------
    
    dim_θh_θhf =
        [ dim_all_nodes,
          no_lines_fault]
    
    #--------------------------------------
    
    θh_θhf_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_θh_θhf;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip( θh_θhf_labels,
                 θh_θhf_Idx )))

end


"""
    get_generic_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)


Returns indices of non-generator nodes active power demand `Png`, non-generator nodes reactive power demand `Qng`, generator nodes local active power demand `Pll`, generator nodes local reactive power demand `Qll` in a flattened `Png_Qng_Pll_Qll` vector.
"""
function get_generic_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx
    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    Png_Qng_Pll_Qll_idx_labels =
        [:dyn_P_non_gens_Idxs,
         :dyn_Q_non_gens_Idxs,       
         :dyn_P_gens_loc_load_Idxs,
         :dyn_Q_gens_loc_load_Idxs]

    #----------------------------------------
    
    dim_Png_Qng_Pll_Qll =
        [ dim_non_gens,
          dim_non_gens,
          dim_gens_with_loc_loads,
          dim_gens_with_loc_loads]
    
    #--------------------------------------
    
    Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_Png_Qng_Pll_Qll;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol, UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(Png_Qng_Pll_Qll_idx_labels,
                Png_Qng_Pll_Qll_Idx )))

end

"""
    get_generic_Pg_Qg_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)


Returns indices of generator nodes active power generation `Pg`, generator nodes reactive power generation `Qg`, non-generator nodes active power demand `Png`, non-generator nodes reactive power demand `Qng`, generator nodes local active power demand `Pll`, generator nodes local reactive power demand `Qll` in a flattened `Pg_Qg_Png_Qng_Pll_Qll` vector.
"""
function get_generic_Pg_Qg_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx
    
    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    Pg_Qg_Png_Qng_Pll_Qll_idx_labels =
        [:dyn_P_gens_Idxs,
         :dyn_Q_gens_Idxs,
         :dyn_P_non_gens_Idxs,
         :dyn_Q_non_gens_Idxs,       
         :dyn_P_gens_loc_load_Idxs,
         :dyn_Q_gens_loc_load_Idxs]

    #----------------------------------------
    
    dim_Pg_Qg_Png_Qng_Pll_Qll =
        [ dim_gens,
          dim_gens,
          dim_non_gens,
          dim_non_gens,
          dim_gens_with_loc_loads,
          dim_gens_with_loc_loads]
    
    #--------------------------------------
    
    Pg_Qg_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_Pg_Qg_Png_Qng_Pll_Qll;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(Pg_Qg_Png_Qng_Pll_Qll_idx_labels,
                Pg_Qg_Png_Qng_Pll_Qll_Idx )))

end


"""
    get_generic_scale_Pg_Qg_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)


Returns indices of scale variable `scale`, generator nodes active power generation `Pg`, generator nodes reactive power generation `Qg`, non-generator nodes active power demand `Png`, non-generator nodes reactive power demand `Qng`, generator nodes local active power demand `Pll`, generator nodes local reactive power demand `Qll` in a flattened `scale_Pg_Qg_Png_Qng_Pll_Qll` vector.
"""
function get_generic_scale_Pg_Qg_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    scale_Pg_Qg_Png_Qng_Pll_Qll_idx_labels =
        [:dyn_scale_Idxs,
         :dyn_P_gens_Idxs,
         :dyn_Q_gens_Idxs,
         :dyn_P_non_gens_Idxs,
         :dyn_Q_non_gens_Idxs,       
         :dyn_P_gens_loc_load_Idxs,
         :dyn_Q_gens_loc_load_Idxs]

    #----------------------------------------
    
    dim_scale_Pg_Qg_Png_Qng_Pll_Qll =
        [ 1,
          dim_gens,
          dim_gens,
          dim_non_gens,
          dim_non_gens,
          dim_gens_with_loc_loads,
          dim_gens_with_loc_loads]
    
    #--------------------------------------
    
    scale_Pg_Qg_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_scale_Pg_Qg_Png_Qng_Pll_Qll;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(scale_Pg_Qg_Png_Qng_Pll_Qll_idx_labels,
                scale_Pg_Qg_Png_Qng_Pll_Qll_Idx )))

end


"""
    get_generic_Pg_Png_Qng_Idx(
        dyn_pf_fun_kwd_net_idxs)


Returns indices of generator nodes active power generation `Pg`, generator nodes reactive power generation `Qg`, non-generator nodes active power demand `Png`, non-generator nodes reactive power demand `Qng` in a flattened `Pg_Png_Qng` vector.
"""
function get_generic_Pg_Png_Qng_Idx(
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    Pg_Png_Qng_idx_labels =
        [:dyn_P_gens_Idxs,
         :dyn_P_non_gens_Idxs,
         :dyn_Q_non_gens_Idxs]

    #----------------------------------------
    
    dim_Pg_Png_Qng =
        [ dim_gens,
          dim_non_gens,
          dim_non_gens]
    
    #--------------------------------------
    
    Pg_Png_Qng_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_Pg_Png_Qng;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(Pg_Png_Qng_idx_labels,
                Pg_Png_Qng_Idx )))

end


"""
    get_generic_scale_Pg_Png_Qng(
        dyn_pf_fun_kwd_net_idxs)


Returns indices of scale variable `scale`, generator nodes active power generation `Pg`, non-generator nodes active power demand `Png`, non-generator nodes reactive power demand `Qng` in a flattened `scale_Pg_Qg_Png_Qng` vector.
"""
function get_generic_scale_Pg_Png_Qng_Idx(
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    scale_Pg_Png_Qng_idx_labels =
        [:dyn_scale_Idxs,
         :dyn_P_gens_Idxs,
         :dyn_P_non_gens_Idxs,
         :dyn_Q_non_gens_Idxs]

    #----------------------------------------
    
    dim_scale_Pg_Png_Qng =
        [ 1,
          dim_gens,
          dim_non_gens,
          dim_non_gens]
    
    #--------------------------------------
    
    scale_Pg_Png_Qng_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_scale_Pg_Png_Qng;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(scale_Pg_Png_Qng_idx_labels,
                scale_Pg_Png_Qng_Idx )))

end


"""
    get_edges_r_x_b_ratio_angle_idx(
        edges_size)

Returns indices of `r`, `x`, `b`, `ratio` and `angle` in a flattened vector `edges_r_x_b_ratio_angle`.
"""
function get_edges_r_x_b_ratio_angle_idx(
    edges_size)

    edges_r_x_b_ratio_angle_idx_labels =
        [:r_Idxs,
         :x_Idxs,
         :b_Idxs,
         :ratio_Idxs,
         :angle_Idxs]
    
    dim_edges_r_x_b_ratio_angle =
        [ edges_size,
          edges_size,
          edges_size,
          edges_size,
          edges_size]
    
    
    edges_r_x_b_ratio_angle_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_edges_r_x_b_ratio_angle;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(edges_r_x_b_ratio_angle_idx_labels,
                edges_r_x_b_ratio_angle_Idx )))
        
end


"""
    get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
        Ynet)


Returns namedtuples of `Ynet_rows_Idxs_in_flattend`, and  `Ynet_real_imag_Idxs_in_flattend`.
"""
function get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
    Ynet)

    dims_row_Ynet = length.(Ynet)

    Ynet_rows_Idxs_in_flattend =
        get_vars_or_paras_Idxs_in_flattend(
                dims_row_Ynet;
        dims_given = true )
        
    Ynet_real_imag_idx_labels =
        [:Ynet_real_Idxs,
         :Ynet_imag_Idxs]

    dim_Ynet_real_imag =
        [ sum(length.(Ynet)),
          sum(length.(Ynet))]
    
    
    Ynet_real_imag_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_Ynet_real_imag;
            dims_given = true )

    Ynet_real_imag_Idxs_in_flattend =
        namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
            a_sym => idx_range
        for (a_sym, idx_range) in
            zip(Ynet_real_imag_idx_labels,
                Ynet_real_imag_Idx )))
    
    return (;Ynet_rows_Idxs_in_flattend,
            Ynet_real_imag_Idxs_in_flattend )
        
end


"""
    get_dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)


Returns indices of `v_ref`, `p_order`, `Png`, `Qng`, `Pll`, `Qll` in a flattened vector `v_ref_p_order_Png_Qng_Pll_Qll`.
"""
function get_dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    
    dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_non_gens,
             dim_non_gens,
             dim_gens_with_loc_loads,
             dim_gens_with_loc_loads];
            dims_given = true )

    (dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx) =
         dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx
    
    return (;dyn_v_ref_Idx,
            dyn_p_order_Idx,
            dyn_Png_Idx,
            dyn_Qng_Idx,
            dyn_Pll_Idx,
            dyn_Qll_Idx )

end


"""
    get_dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)


Returns indices of `ω_ref`, `v_ref`, `p_order`, `Png`, `Qng`, `Pll`, `Qll` in a flattened vector `v_ref_p_order_Png_Qng_Pll_Qll`.
"""
function get_dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    
    dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_non_gens,
             dim_non_gens,
             dim_gens_with_loc_loads,
             dim_gens_with_loc_loads];
            dims_given = true )

    (dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx) =
         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx

    
    return (;dyn_ω_ref_Idx,
             dyn_v_ref_Idx,
             dyn_p_order_Idx,
             dyn_Png_Idx,
             dyn_Qng_Idx,
             dyn_Pll_Idx,
             dyn_Qll_Idx )

end


"""
    get_dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)


Returns indices of `δ`, `eq_dash`, `Png`, `Qng`, `Pll`, `Qll` in a flattened vector `δ_eq_dash_Png_Qng_Pll_Qll`.
"""
function get_dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    
    dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_non_gens,
             dim_non_gens,
             dim_gens_with_loc_loads,
             dim_gens_with_loc_loads];
            dims_given = true )

    (dyn_δ_Idx,
     dyn_eq_dash_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx) =
         dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx
    
    return (;
            dyn_δ_Idx,
            dyn_eq_dash_Idx,
            dyn_Png_Idx,
            dyn_Qng_Idx,
            dyn_Pll_Idx,
            dyn_Qll_Idx )

end

"""
    get_dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx(
        dyn_pf_fun_kwd_net_idxs)


Returns indices of `δ`, `ed_dash`, `eq_dash`, `Png`, `Qng`, `Pll`, `Qll` in a flattened vector `δ_ed_dash_eq_dash_Png_Qng_Pll_Qll`.
"""
function get_dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    
    dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_non_gens,
             dim_non_gens,
             dim_gens_with_loc_loads,
             dim_gens_with_loc_loads];
            dims_given = true )

    (dyn_δ_Idx,
     dyn_ed_dash_Idx,
     dyn_eq_dash_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx) =
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx
    
    return (;
            dyn_δ_Idx,
            dyn_ed_dash_Idx,
            dyn_eq_dash_Idx,
            dyn_Png_Idx,
            dyn_Qng_Idx,
            dyn_Pll_Idx,
            dyn_Qll_Idx )

end

"""
    get_dyn_vh_id_iq_V_ref_Tm_Idx(
        gens_nodes_idx;
        reverse_idx = false)


Returns indices of `vh`, `id`, `iq`, `V_ref`, `Tm` in a flattened vector `vh_id_iq_V_ref_Tm`.
"""
function get_dyn_vh_id_iq_V_ref_Tm_Idx(
    gens_nodes_idx;
    reverse_idx = false)

    dim_gens =
        length(gens_nodes_idx)

    
    dyn_vh_id_iq_V_ref_Tm_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    if reverse_idx == false

        (dyn_vh_Idx,
         dyn_id_Idx,
         dyn_iq_Idx,
         dyn_V_ref_Idx,
         dyn_Tm_Idx) =
             dyn_vh_id_iq_V_ref_Tm_Idx
        
        return (;
                dyn_vh_Idx,
                dyn_id_Idx,
                dyn_iq_Idx,
                dyn_V_ref_Idx,
                dyn_Tm_Idx )
    else

        (dyn_Tm_Idx,
         dyn_V_ref_Idx,
         dyn_id_Idx,
         dyn_iq_Idx,
         dyn_vh_Idx) =
             dyn_vh_id_iq_V_ref_Tm_Idx
        
        return (;dyn_Tm_Idx, 
                dyn_V_ref_Idx,
                dyn_id_Idx,
                dyn_iq_Idx,
                dyn_vh_Idx)
        
    end
    
end


"""
    get_dyn_vh_id_iq_V_ref_Tm_Idx(
        gens_nodes_idx)


Returns indices of `V_ref`, `Tm`, `vh`, `id`, `iq` in a flattened vector `V_ref_Tm_vh_id_iq`.
"""
function get_dyn_V_ref_Tm_vh_id_iq_Idx(
    gens_nodes_idx)

    dim_gens =
        length(gens_nodes_idx)

    
    dyn_V_ref_Tm_vh_id_iq_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    (dyn_V_ref_Idx,
     dyn_Tm_Idx,               
     dyn_id_Idx,
     dyn_iq_Idx,
     dyn_vh_Idx) =
         dyn_V_ref_Tm_vh_id_iq_Idx

    return (;dyn_V_ref_Idx,
            dyn_Tm_Idx,               
            dyn_id_Idx,
            dyn_iq_Idx,
            dyn_vh_Idx)
    
end


"""
    get_dyn_vh_id_iq_ωref0_vref0_porder0_Idx(
        gens_nodes_idx;
        reverse_idx = false)


Returns indices of `vh`, `id`, `iq`, `ωref0`, `vref0`, `porder0` in a flattened vector `vh_id_iq_ωref0_vref0_porder0`.
"""
function get_dyn_vh_id_iq_ωref0_vref0_porder0_Idx(
    gens_nodes_idx;
    reverse_idx = false )

    dim_gens =
        length(gens_nodes_idx)

    
    dyn_vh_id_iq_ωref0_vref0_porder0_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    if reverse_idx == false

    (dyn_vh_Idx,
     dyn_id_Idx,
     dyn_iq_Idx,
     dyn_ωref0_Idx,
     dyn_vref0_Idx,
     dyn_porder0_Idx) =
         dyn_vh_id_iq_ωref0_vref0_porder0_Idx
        
        return (;
                dyn_vh_Idx,
                dyn_id_Idx,
                dyn_iq_Idx,
                dyn_ωref0_Idx,
                dyn_vref0_Idx,
                dyn_porder0_Idx )
    else

    (dyn_porder0_Idx,
     dyn_vref0_Idx,
     dyn_ωref0_Idx,
     dyn_id_Idx,
     dyn_iq_Idx,
     dyn_vh_Idx) =
         dyn_vh_id_iq_ωref0_vref0_porder0_Idx
        
        return (;
                dyn_porder0_Idx,
                dyn_vref0_Idx,
                dyn_ωref0_Idx,
                dyn_id_Idx,
                dyn_iq_Idx,
                dyn_vh_Idx)
        
    end
    
end


"""
    get_ωref0_vref0_porder0_id_iq_vh_Idx(
        gens_nodes_idx)


Returns indices of `ωref0`, `vref0`, `porder0`, `id`, `iq`, `vh` in a flattened vector `ωref0_vref0_porder0_id_iq_vh`.
"""
function get_ωref0_vref0_porder0_id_iq_vh_Idx(
    gens_nodes_idx)

    dim_gens =
        length(gens_nodes_idx)

    
    ωref0_vref0_porder0_id_iq_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    (ωref0_Idx,
     vref0_Idx,
     porder0_Idx,     
     id_Idx,
     iq_Idx,
     vh_Idx
     ) =
         ωref0_vref0_porder0_id_iq_vh_Idx

    
     return   (;
               ωref0_Idx,
               vref0_Idx,
               porder0_Idx,     
               id_Idx,
               iq_Idx,
               vh_Idx
               )
    
end


"""
    get_dyn_ωref0_vref0_porder0_id_iq_vh_Idx(
        gens_nodes_idx)


Returns indices of `ωref0`, `vref0`, `porder0`, `id`, `iq`, `vh` in a flattened vector `ωref0_vref0_porder0_id_iq_vh`.
"""
function get_dyn_ωref0_vref0_porder0_id_iq_vh_Idx(
    gens_nodes_idx)

    dim_gens =
        length(gens_nodes_idx)

    
    dyn_ωref0_vref0_porder0_id_iq_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    (dyn_ωref0_Idx,
     dyn_vref0_Idx,
     dyn_porder0_Idx,     
     dyn_id_Idx,
     dyn_iq_Idx,
     dyn_vh_Idx
     ) =
         dyn_ωref0_vref0_porder0_id_iq_vh_Idx

    
     return   (;
               dyn_ωref0_Idx,
               dyn_vref0_Idx,
               dyn_porder0_Idx,     
               dyn_id_Idx,
               dyn_iq_Idx,
               dyn_vh_Idx
               )
    
end


"""
    get_id_iq_pg_vh_Idx(
        gens_nodes_idx)


Returns indices of `id`, `iq`, `pg`, `vh` in a flattened vector `id_iq_pg_vh`.
"""
function get_id_iq_pg_vh_Idx(
    gens_nodes_idx)

    dim_gens =
        length(gens_nodes_idx)

    #--------------------------------------

    id_iq_pg_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    (id_Idx,
     iq_Idx,
     pg_Idx,
     vh_Idx)=
         id_iq_pg_vh_Idx

    
    return (;id_Idx,
            iq_Idx,
            pg_Idx,
            vh_Idx )

end



function get_generic_n_net_comp_type_paras_Idx(
    a_net_comp_type_idx;
    no_paras = 5,
    comp_type = "gens")

    #----------------------------------------

    dim_net_comp_type =
        length(a_net_comp_type_idx)

    #----------------------------------------

    n_paras_idx_labels =
        Symbol.(["$(comp_type)_para_$(a_para)_Idx"
                 for a_para in 1:no_paras])
    
    #----------------------------------------

    dim_n_paras =
        [dim_net_comp_type
         for a_para in
             1:no_paras ]
    
    #--------------------------------------
    
    n_net_comp_type_paras_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_n_paras;
            dims_given = true )

    return namedtuple(OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(n_paras_idx_labels,
                n_net_comp_type_paras_Idx )))

end



function get_generic_n_gens_paras_wt_Png_Qng_Pll_Qll_Idx(
    dyn_pf_fun_kwd_net_idxs;
    no_gens_paras = 3)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    #----------------------------------------

    n_gens_paras_idx_labels =
        Symbol.(["gens_para_$(a_para)_Idx"
                 for a_para in 1:no_gens_paras])
    
    Png_Qng_Pll_Qll_idx_labels =
        [:dyn_P_non_gens_Idxs,
         :dyn_Q_non_gens_Idxs,       
         :dyn_P_gens_loc_load_Idxs,
         :dyn_Q_gens_loc_load_Idxs]

    n_gens_paras_wt_Png_Qng_Pll_Qll_idx_labels =
        [n_gens_paras_idx_labels;
         Png_Qng_Pll_Qll_idx_labels ]
    #----------------------------------------

    dim_n_gens_paras =
        [dim_gens
         for a_para in
             1:no_gens_paras ]
    
    dim_Png_Qng_Pll_Qll =
        [ dim_non_gens,
          dim_non_gens,
          dim_gens_with_loc_loads,
          dim_gens_with_loc_loads]
    
    dim_n_gens_paras_wt_Png_Qng_Pll_Qll =
        [dim_n_gens_paras;
         dim_Png_Qng_Pll_Qll]
    
    #--------------------------------------
    
    n_gens_paras_wt_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            dim_n_gens_paras_wt_Png_Qng_Pll_Qll;
            dims_given = true )

    return namedtuple(
        OrderedDict{Symbol,UnitRange{Int64}}(
        a_sym => idx_range
        for (a_sym, idx_range) in
            zip(n_gens_paras_wt_Png_Qng_Pll_Qll_idx_labels,
                n_gens_paras_wt_Png_Qng_Pll_Qll_Idx )))

end

#-----------------------------------------------------
#-----------------------------------------------------

"""
    get_generic_algebraic_state_sym(
        gens_nodes_idx,
        all_nodes_idx)



Returns algebraic variables symbols.
"""
function get_generic_algebraic_state_sym(
    gens_nodes_idx,
    all_nodes_idx;
    label_prefix = "bus")

    return  [generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:vh];
            label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:θh];
             label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:id];
             label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:iq];
             label_prefix = label_prefix) ]


end

"""
    get_generic_nodes_names(
        dyn_pf_fun_kwd_net_idxs,
        n2s_all_nodes_idx)


Returns namedtuples of nodes types names `all_nodes_names`, `gens_nodes_names`, `non_gens_nodes_names`.


# Arguments
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.

"""
function get_generic_nodes_names(
    dyn_pf_fun_kwd_net_idxs,
    n2s_all_nodes_idx)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx
    
    all_nodes_names =
        ["bus$(n2s_all_nodes_idx[idx])"
         for idx in all_nodes_idx ]

            
    gens_nodes_names =
        ["bus$(idx)"
         for idx in gens_nodes_idx ]

            
    non_gens_nodes_names =
        ["bus$(idx)"
         for idx in non_gens_nodes_idx ]

    return (;all_nodes_names,
            gens_nodes_names,
            non_gens_nodes_names)

end


#-----------------------------------------------------


"""
    get_generic_network_vars_labels(
        plants_states_syms,
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs
        ;label_prefix = "bus",
        plants_states_by_per_comp = false,
        plants_states_by_per_plant = false
        )


Returns namedtuples of `state_labels`, `algebraic_vars_labels`, `network_vars_labels`.


# Arguments
- `plants_states_syms`: the plants states variables symbols per plant for all generator plants. 
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.

"""
function get_generic_network_vars_labels(
    plants_states_syms,
    dyn_pf_fun_kwd_net_idxs,
    dyn_pf_fun_kwd_n2s_idxs
    ;label_prefix = "bus",
    plants_states_by_per_comp = false,
    plants_states_by_per_plant = false
    )

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx


   (;n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    
    if( plants_states_by_per_comp == true &&
        plants_states_by_per_plant == false)

        plants_states_syms =
            [ :nothing ∈ a_comp.gov ?
            [a_comp.gen; a_comp.avr ] :
            [a_comp.gen; a_comp.avr; a_comp.gov ]
             for a_comp in
                 plants_states_syms ]

        state_labels =
            get_labels_by_nodes_idxs_and_vec_vec_syms(
                gens_nodes_idx,
                plants_states_syms;
                label_prefix = label_prefix )
        
    elseif (plants_states_by_per_plant == true &&
        plants_states_by_per_comp == false)
    
        state_labels =
            get_labels_by_nodes_idxs_and_vec_vec_syms(
                gens_nodes_idx,
                plants_states_syms;
                label_prefix = label_prefix )

    else
            
        state_labels =
            generate_labels_by_nodes_idxs_and_vars(
                gens_nodes_idx,
                plants_states_syms;
                label_prefix = label_prefix )
        
    end
    
    #----------------------------------------

    algebraic_vars_labels =
        [generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:vh];
            label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:θh];
             label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:id];
             label_prefix = label_prefix);
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:iq];
             label_prefix = label_prefix) ]

    #----------------------------------------

    network_vars_labels =
        [state_labels;
         algebraic_vars_labels]

    return (;state_labels,
            algebraic_vars_labels,
            network_vars_labels)
    

end

#-----------------------------------------------------


"""
    get_gens_state_vars_idx_in_state(
        network_vars_labels,
        # all_nodes_idx,
        dyn_pf_fun_kwd_net_idxs,
        n2s_all_nodes_idx;
        selected_gens_state_vars_syms =
            (:δ, :ed_dash, :eq_dash) )


Returns state variables indices in the state for a selected set of variables in `selected_gens_state_vars_syms`.


# Arguments
- `network_vars_labels`: the system variables labels.
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `selected_gens_state_vars_syms::Tuple{Symbol}=(:δ, :ed_dash, :eq_dash)`: the selected
   state variables. 

"""
function get_gens_state_vars_idx_in_state(
    network_vars_labels,
    # all_nodes_idx,
    dyn_pf_fun_kwd_net_idxs,
    n2s_all_nodes_idx;
    selected_gens_state_vars_syms =
        (:δ, :ed_dash, :eq_dash) )

    dim_selections =
        length( selected_gens_state_vars_syms )

        
    vec_vec = Vector{Int64}[
        [] for a_para in 1:dim_selections ]


    (;gens_nodes_names,) =
        NamedTupleTools.select(
            get_generic_nodes_names(
    dyn_pf_fun_kwd_net_idxs,
    n2s_all_nodes_idx),
             # get_generic_nodes_names(
             #     all_nodes_idx,
             #     n2s_all_nodes_idx) ,
             (:gens_nodes_names,))
    
    # δ_ω_ed_dash_eq_dash_indices
    
    tup_gens_state_vars_idx_in_state =
        get_nodes_state_algb_vars_indices_in_system(
            ;network_vars_labels =
                network_vars_labels,
            nodes_name = gens_nodes_names,
            vars = selected_gens_state_vars_syms )

    vec_gens_state_vars_idx =
        [[idx for idx in a_state_idx_in_state]
         for a_state_idx_in_state in
             tup_gens_state_vars_idx_in_state]

        
    for (idx, a_sym) in enumerate(
        selected_gens_state_vars_syms)
        for a_vec in vec_gens_state_vars_idx
            push!(vec_vec[idx],
                  a_vec[idx] )
        end
    end


    gens_state_vars_idx_in_state_syms =
        Symbol.(["$(a_sym)_idx_in_state"
                 for a_sym in
                     selected_gens_state_vars_syms])
    
    return namedtuple(OrderedDict{Symbol, Vector{Int64}}(
        a_sym => a_state_idx
        for (a_sym, a_state_idx) in
            zip(gens_state_vars_idx_in_state_syms,
                vec_vec )))

end


"""
    get_state_vars_and_i_dq_Idx_in_state(
        generic_state_labels,
        gens_nodes_idx,
        all_nodes_idx )


Returns state variables, nodes voltage magnitudes, nodes voltage angles,  generators direct axix current id, and generators quadrature axis current indices in the state.


# Arguments
- `network_vars_labels`: the system variables labels.
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `selected_gens_state_vars_syms::Tuple{Symbol}=(:δ, :ed_dash, :eq_dash)`:  

"""
function get_state_vars_and_i_dq_Idx_in_state(
    generic_state_labels,
    gens_nodes_idx,
    all_nodes_idx )

    dim_pure_states =
        length(generic_state_labels)
    
    dim_gens =
        length(gens_nodes_idx)
    
    dim_all_nodes =
        length(all_nodes_idx)

    state_vars_and_i_dq_Idx_in_state = 
        get_vars_or_paras_Idxs_in_flattend(
            [dim_pure_states,
             dim_all_nodes,
             dim_all_nodes,
             dim_gens,
             dim_gens
             ];
            dims_given = true )
    
    (state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state) =
         state_vars_and_i_dq_Idx_in_state
    
    return (;state_var_Idx_in_state,
            vh_Idx_in_state,
            θh_Idx_in_state,
            id_Idx_in_state,
            iq_Idx_in_state)
    
end


"""
    get_state_vars_and_i_dq_wt_fault_Idx_in_state(
        generic_state_labels,
        gens_nodes_idx,
        all_nodes_idx;
        no_lines_fault = 1)


Returns state variables, nodes voltage magnitudes, nodes voltage angles,  generators direct axix current id, and generators quadrature axis current, fault node voltage magnitudes, fault node voltage angles,  indices in the state.


# Arguments
- `generic_state_labels`: the states variables labels.  

"""
function get_state_vars_and_i_dq_wt_fault_Idx_in_state(
    generic_state_labels,
    gens_nodes_idx,
    all_nodes_idx;
    no_lines_fault = 1)

    dim_pure_states =
        length(generic_state_labels)
    
    dim_gens =
        length(gens_nodes_idx)
    
    dim_all_nodes =
        length(all_nodes_idx)

    state_vars_and_i_dq_wt_fault_Idx_in_state = 
        get_vars_or_paras_Idxs_in_flattend(
            [dim_pure_states,
             dim_all_nodes,
             dim_all_nodes,
             dim_gens,
             dim_gens,
             no_lines_fault,
             no_lines_fault
             ];
            dims_given = true )
    
    (state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state,
     vhf_Idx_in_state,
     θhf_Idx_in_state) =
         state_vars_and_i_dq_wt_fault_Idx_in_state
    
    return (;state_var_Idx_in_state,
            vh_Idx_in_state,
            θh_Idx_in_state,
            id_Idx_in_state,
            iq_Idx_in_state,
            vhf_Idx_in_state,
            θhf_Idx_in_state)
    
end



"""
    get_state_algebraic_vars_Idx_in_state(
        generic_state_labels,
        gens_nodes_idx,
        all_nodes_idx)


Returns state variables, and algebraic variables indices in the state.


# Arguments
- `generic_state_labels`: the states variables labels.  

"""
function get_state_algebraic_vars_Idx_in_state(
    generic_state_labels,
    gens_nodes_idx,
    all_nodes_idx )

    dim_pure_states =
        length(generic_state_labels)
    
    dim_gens =
        length(gens_nodes_idx)
    
    dim_all_nodes =
        length(all_nodes_idx)

    dim_algebraic_vars =
        2 * dim_all_nodes + 2 * dim_gens

    state_vars_and_algebraic_Idx_in_state = 
        get_vars_or_paras_Idxs_in_flattend(
            [dim_pure_states,
             dim_algebraic_vars ];
            dims_given = true )
    
    (state_var_Idx_in_state,
     algebraic_var_Idx_in_state) =
         state_vars_and_algebraic_Idx_in_state
    
    return (;state_var_Idx_in_state,
            algebraic_var_Idx_in_state)
    
end


"""
    get_state_algebraic_vars_wt_fault_Idx_in_state(
        generic_state_labels,
        gens_nodes_idx,
        all_nodes_idx;
        no_lines_fault = 1)


Returns state variables, and algebraic variables with fault nodes algebraic variables indices in the state.


# Arguments
- `generic_state_labels`: the states variables labels.  

"""
function get_state_algebraic_vars_wt_fault_Idx_in_state(
    generic_state_labels,
    gens_nodes_idx,
    all_nodes_idx;
    no_lines_fault = 1)

    dim_pure_states =
        length(generic_state_labels)
    
    dim_gens =
        length(gens_nodes_idx)
    
    dim_all_nodes =
        length(all_nodes_idx)

    dim_algebraic_vars_wt_fault =
        2 * dim_all_nodes + 2 * dim_gens + 2 * no_lines_fault 
    

    state_vars_and_algebraic_wt_fault_Idx_in_state = 
        get_vars_or_paras_Idxs_in_flattend(
            [dim_pure_states,
             dim_algebraic_vars_wt_fault ];
            dims_given = true )
    
    (state_var_Idx_in_state,
     algebraic_wt_fault_Idx_in_state) =
         state_vars_and_algebraic_wt_fault_Idx_in_state
    
    return (;state_var_Idx_in_state,
            algebraic_wt_fault_Idx_in_state)
    
end


# ------------------------------------------------------
# ------------------------------------------------------
# Static algebraic, parameter and  indices aggregation
# ------------------------------------------------------
# ------------------------------------------------------

"""
    get_static_Idx_and_syms(
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs;
        no_lines_fault = 1)


Returns  states variables indices, algebraic variables indices, states variables labels, algebraic variables labels and most parameters indices.


# Arguments
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.
"""
function get_static_Idx_and_syms(
    dyn_pf_fun_kwd_net_idxs,
    dyn_pf_fun_kwd_n2s_idxs;
    no_lines_fault = 1)
    
    #-------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

   (;n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))

    #----------------------------------------

    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]
    
    #----------------------------------------
    # Dimensions 
    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

            
    #----------------------------------------
    # Labels, syms and nodes names
    #----------------------------------------

    algebraic_state_sym =
        [generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:vh];
            label_prefix = "bus");
         generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:θh];
             label_prefix = "bus");
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:id];
             label_prefix = "bus");
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:iq];
             label_prefix = "bus") ]

        
    nodes_names =
        ["bus$(n2s_all_nodes_idx[idx])"
         for idx in all_nodes_idx ]

            
    # gens_nodes_names =
    #     nodes_names[ gens_nodes_idx ]

            
    # non_gens_nodes_names =
    #     nodes_names[ non_gens_nodes_idx ]
            
    gens_nodes_names =
        nodes_names[ transformed_gens_nodes_idx ]

            
    non_gens_nodes_names =
        nodes_names[ transformed_non_gens_nodes_idx ]
    
    #----------------------------------------
    # Parameters indices 
    #----------------------------------------
    
    Pg_Qg_Png_Qng_Pll_Qll_Idx =
        get_generic_Pg_Qg_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)

    #--------------------------------------
    
    scale_Pg_Qg_Png_Qng_Pll_Qll_Idx =
        get_generic_scale_Pg_Qg_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)
    
    #--------------------------------------
    
    Png_Qng_Pll_Qll_Idx =
        get_generic_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)

    #--------------------------------------
    
    Pg_Png_Qng_Idx =
        get_generic_Pg_Png_Qng_Idx(
            dyn_pf_fun_kwd_net_idxs)

    #--------------------------------------

    scale_Pg_Png_Qng_Idx =
        get_generic_scale_Pg_Png_Qng_Idx(
    dyn_pf_fun_kwd_net_idxs)


    #--------------------------------------
    
    pf_vh_θh_idx_and_idx2Idx =
        get_pf_vh_θh_idx_and_idx2Idx(
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs)
    
    #--------------------------------------

    dyn_pf_flat_vh_flat_θh_Idx =
        get_generic_flat_vh_flat_θh_Idx(
            gens_nodes_idx,
            all_nodes_idx)

    #--------------------------------------

    dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx =
        get_generic_flat_vh_flat_θh_wt_slack_value_Idx(
            all_nodes_idx)
    
    #--------------------------------------
    
    dyn_pf_vh_vhf_Idx =
        get_generic_vh_vhf_Idx(
            all_nodes_idx;
            no_lines_fault = no_lines_fault)
    
    #--------------------------------------

    dyn_pf_θh_θhf_Idx =
        get_generic_θh_θhf_Idx(
            all_nodes_idx;
            no_lines_fault = no_lines_fault)
    
    #----------------------------------------
    #----------------------------------------

    return (;algebraic_state_sym,
            nodes_names,
            gens_nodes_names,
            non_gens_nodes_names,
            Pg_Qg_Png_Qng_Pll_Qll_Idx,
            scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
            Png_Qng_Pll_Qll_Idx,
            Pg_Png_Qng_Idx,
            
            scale_Pg_Png_Qng_Idx,
            pf_vh_θh_idx_and_idx2Idx,
            
            dyn_pf_flat_vh_flat_θh_Idx,
            dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
            
            dyn_pf_vh_vhf_Idx,
            dyn_pf_θh_θhf_Idx)
    
end



# ------------------------------------------------------
# ------------------------------------------------------
# States, algebraic, parameter and function aggregation
# ------------------------------------------------------
# ----------------------------------------------------

"""
    get_states_Idx_syms_wt_functions(
        net_data_by_components_file,    
        dyn_pf_fun_kwd_net_idxs,
        dyn_pf_fun_kwd_n2s_idxs;
        components_libs_dir =
            nothing,
        no_lines_fault = 1 )


Returns dynamic functions of components, states variables indices, algebraic variables indices, states variables labels, algebraic variables labels and most parameters indices.


# Arguments
- `net_data_by_components_file`: the network data file.
- `dyn_pf_fun_kwd_net_idxs`: the namedtuple of node's type indices.
- `dyn_pf_fun_kwd_n2s_idxs`: the namedtuple of node's type indices translation dictionaries.
- `components_libs_dir`: the components library folder.

"""
function get_states_Idx_syms_wt_functions(
    net_data_by_components_file,    
    dyn_pf_fun_kwd_net_idxs,
    dyn_pf_fun_kwd_n2s_idxs;
    components_libs_dir =
        nothing,
    no_lines_fault = 1 )

    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir,
                     "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end
    
    #-------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))

    # gens_nodes_with_loc_loads_idx =
    #     gens_with_loc_loads_idx

   (;n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------
    # Dimensions 
    #----------------------------------------

    dim_gens =
        length(gens_nodes_idx)

    dim_non_gens =
        length(non_gens_nodes_idx)

    dim_gens_with_loc_loads =
        length(gens_nodes_with_loc_loads_idx)

    dim_all_nodes = length(all_nodes_idx)

    
    #-------------------------------

    (gens_govs_avrs_states_syms,
     gens_govs_avrs_types) =
        get_gens_govs_avrs_states_syms_by_json(
            net_data_by_components_file;
            components_libs_dir =
                components_libs_dir )
    
    #-------------------------------

    vec_vec_gens_govs_avrs_states_syms =
        [ :nothing ∈ an_item.gov ?
        [an_item.gen, an_item.avr ] :
        [an_item.gen, an_item.avr, an_item.gov ]
         for an_item in
             gens_govs_avrs_states_syms]

    vec_dims_gens_govs_avrs_states_syms =
        [:nothing ∈ an_item.gov ?
        [length(an_item.gen),
         length(an_item.avr)] :
             [length(an_item.gen),
              length(an_item.avr),
              length(an_item.gov) ]
         for an_item in
             gens_govs_avrs_states_syms]

    dims_gens_govs_avrs_states_syms =
        sum.(vec_dims_gens_govs_avrs_states_syms)

    state_vars_idx =
        get_vars_or_paras_Idxs_in_flattend(
            dims_gens_govs_avrs_states_syms;
            dims_given = true )

    vec_comp_states_Idx =
        map((x) ->
        get_vars_or_paras_Idxs_in_flattend(
            x;
            dims_given = true ),
            vec_dims_gens_govs_avrs_states_syms )

    plants_states_syms =
        [ :nothing ∈ an_item.gov ?
        [an_item.gen; an_item.avr ] :
        [an_item.gen; an_item.avr; an_item.gov ]
         for an_item in
             gens_govs_avrs_states_syms]

    generic_state_sym =
        get_labels_by_nodes_idxs_and_vec_vec_syms(
            gens_nodes_idx,
            plants_states_syms;
            label_prefix = "bus" )

    
    (;state_labels,
     algebraic_vars_labels,
     network_vars_labels) =
        NamedTupleTools.select(
            get_generic_network_vars_labels(
                plants_states_syms,
                dyn_pf_fun_kwd_net_idxs,
                dyn_pf_fun_kwd_n2s_idxs;
                label_prefix =
                    "bus",
                plants_states_by_per_comp =
                    false,
                plants_states_by_per_plant =
                    true,),
            (:state_labels,
             :algebraic_vars_labels,
             :network_vars_labels))
    

    #----------------------------------------
    #----------------------------------------

    list_gens_type_syms =
        Symbol.([a_type.gen for a_type in
             gens_govs_avrs_types])


    list_govs_type_syms =
        Symbol.([a_type.gov for a_type in
             gens_govs_avrs_types])


    list_avrs_type_syms =
        Symbol.([a_type.avr for a_type in
             gens_govs_avrs_types])

    #----------------------------------------

    gens_cb_paras_func =
        get_gens_callback_paras_func(
         list_gens_type_syms)

    gens_init_func =
        get_gens_init_func(
         list_gens_type_syms)
    
    gens_output_func =
        get_gens_output_func(
            list_gens_type_syms)
    
    ode_gens_dyn_func =
        get_gens_dyn_func(
            list_gens_type_syms,
            :ode)

    
    dae_gens_dyn_func =
        get_gens_dyn_func(
            list_gens_type_syms,
            :dae)

    
    gens_dyn_func =
        get_gens_dyn_func(
            list_gens_type_syms)
    
    

    ##

    
    govs_cb_paras_func =
        get_govs_callback_paras_func(
            list_govs_type_syms)
    
    govs_init_func =
        get_govs_init_func(
         list_govs_type_syms)
    
    govs_output_func =
        get_govs_output_func(
            list_govs_type_syms)
    
    ode_govs_dyn_func =
        get_govs_dyn_func(
            list_govs_type_syms,
            :ode)

    
    dae_govs_dyn_func =
        get_govs_dyn_func(
            list_govs_type_syms,
            :dae)

    
    govs_dyn_func =
        get_govs_dyn_func(
            list_govs_type_syms)
    
    
    ##

    
    avrs_cb_paras_func =
        get_avrs_callback_paras_func(
            list_avrs_type_syms)
    
    avrs_init_func =
        get_avrs_init_func(
         list_avrs_type_syms)
    
    avrs_output_func =
        get_avrs_output_func(
            list_avrs_type_syms)

    
    ode_avrs_dyn_func =
        get_avrs_dyn_func(
            list_avrs_type_syms,
            :ode)

    
    dae_avrs_dyn_func =
        get_avrs_dyn_func(
            list_avrs_type_syms,
            :dae)
    
    avrs_dyn_func =
        get_avrs_dyn_func(
            list_avrs_type_syms)

    
    #----------------------------------------

    comps_callback_paras_funs = [
        (gen_cb_paras_fun = a_gen_fun,
         avr_cb_paras_fun = a_avr_fun,
         gov_cb_paras_fun = a_gov_fun) 
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( gens_cb_paras_func,
                    avrs_cb_paras_func,
                    govs_cb_paras_func)]

    comps_init_funs = [
        (gen_init_fun = a_gen_fun,
         avr_init_fun = a_avr_fun,
         gov_init_fun = a_gov_fun)
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( gens_init_func,
                    avrs_init_func,
                    govs_init_func)]

    
    comps_output_funs = [
        (gen_output_fun = a_gen_fun,
         avr_output_fun = a_avr_fun,
         gov_output_fun = a_gov_fun) 
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( gens_output_func,
                    avrs_output_func,
                    govs_output_func)]


    
    ode_comps_dyn_funs = [
        (gen_dyn_fun = a_gen_fun,
         avr_dyn_fun = a_avr_fun,
         gov_dyn_fun = a_gov_fun)
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( ode_gens_dyn_func,
                    ode_avrs_dyn_func,
                    ode_govs_dyn_func)]

    
    dae_comps_dyn_funs = [
        (gen_dyn_fun = a_gen_fun,
         avr_dyn_fun = a_avr_fun,
         gov_dyn_fun = a_gov_fun)
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( dae_gens_dyn_func,
                    dae_avrs_dyn_func,
                    dae_govs_dyn_func)]
    
    comps_dyn_funs = [
        (gen_dyn_fun = a_gen_fun,
         avr_dyn_fun = a_avr_fun,
         gov_dyn_fun = a_gov_fun)
        for (a_gen_fun,
             a_avr_fun,
             a_gov_fun)
            in zip( gens_dyn_func,
                    avrs_dyn_func,
                    govs_dyn_func)]
        
    #----------------------------------------
    # Labels, syms and nodes names
    #----------------------------------------

    algebraic_state_sym =
        [generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:vh];
            label_prefix = "bus");
         generate_labels_by_nodes_idxs_and_vars(
            all_nodes_idx,
            [:θh];
             label_prefix = "bus");
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:id];
             label_prefix = "bus");
         generate_labels_by_nodes_idxs_and_vars(
            gens_nodes_idx,
            [:iq];
             label_prefix = "bus") ]

    #----------------------------------------

    model_syms =
        [state_labels;
         algebraic_vars_labels]

    #----------------------------------------
        
    nodes_names =
        ["bus$(n2s_all_nodes_idx[idx])"
         for idx in all_nodes_idx ]

            
    gens_nodes_names = nodes_names[ [n2s_all_nodes_idx[idx]
                     for idx in gens_nodes_idx] ]
            
    non_gens_nodes_names =
        nodes_names[ [n2s_all_nodes_idx[idx]
                     for idx in non_gens_nodes_idx] ]
            
    # gens_nodes_names =
    #     nodes_names[ gens_nodes_idx ]

            
    # non_gens_nodes_names =
    #     nodes_names[ non_gens_nodes_idx ]


    #-------------------------------
    # I need to get the indices of SM and SC gens

    SM_gens_idx = [ idx
         for (idx, an_item) in
                    enumerate(gens_govs_avrs_states_syms)
                    if  :nothing ∉ an_item.gov ]


    SC_gens_idx = [ idx
         for (idx, an_item) in
                    enumerate(gens_govs_avrs_states_syms)
                    if  :nothing ∈ an_item.gov ]

    SM_gens_nodes_names =
        gens_nodes_names[ SM_gens_idx ]

    SC_gens_nodes_names =
        gens_nodes_names[ SC_gens_idx]
    
    #----------------------------------------
    # States indices 
    #----------------------------------------
    
    δ_ω_ed_dash_eq_dash_indices =
        get_nodes_state_algb_vars_indices_in_system(
            ;network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_names,
            vars = [:δ, :ω, :ed_dash, :eq_dash] )

    δ_idx_in_state =
        [idx for idx in
             first.(δ_ω_ed_dash_eq_dash_indices)]
    
    ω_idx_in_state =
        [idx for idx in
             second.(δ_ω_ed_dash_eq_dash_indices)]
    
    ed_dash_idx_in_state =
        [idx for idx in
             third.(δ_ω_ed_dash_eq_dash_indices)]
    
    eq_dash_idx_in_state =
        [idx for idx in
             fourth.(δ_ω_ed_dash_eq_dash_indices)]

    generic_model_states_comp_idxs_in_Idx =
        (;δ_idx_in_state,
          ω_idx_in_state,
          ed_dash_idx_in_state,
          eq_dash_idx_in_state )
        
    #----------------------------------------
    
   generic_model_vars_wt_i_dq_Idx_in_state =
        get_vars_or_paras_Idxs_in_flattend(
            [length( generic_state_sym ),
             length( all_nodes_idx ),
             length( all_nodes_idx ),
             length( gens_nodes_idx),
             length( gens_nodes_idx) ];
            dims_given = true )
    
    
    (state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state) =
         generic_model_vars_wt_i_dq_Idx_in_state 

    
    generic_model_vars_wt_i_dq_Idx_in_state =
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state)
    
    generic_model_vars_no_i_dq_Idx_in_state = 
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state )

    #--------------------------------------
    
    gens_state_vars_idx_in_state =
        get_gens_state_vars_idx_in_state(
        state_labels,
        # all_nodes_idx,
        dyn_pf_fun_kwd_net_idxs,
        n2s_all_nodes_idx;
        selected_gens_state_vars_syms =
            (:δ, :ed_dash, :eq_dash) )
    
    # gens_state_vars_idx_in_state =
    #     get_gens_state_vars_idx_in_state(
    #         state_labels, 
    #         all_nodes_idx,
    #         n2s_all_nodes_idx;
    #         selected_gens_state_vars_syms =
    #             (:δ, :ed_dash, :eq_dash) )
    
    state_vars_and_i_dq_Idx_in_state =
        get_state_vars_and_i_dq_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx )

    state_vars_and_i_dq_wt_fault_Idx_in_state =
        get_state_vars_and_i_dq_wt_fault_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault =
                no_lines_fault)

    state_algebraic_vars_Idx_in_state =
        get_state_algebraic_vars_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx )

    state_algebraic_vars_wt_fault_Idx_in_state =
        get_state_algebraic_vars_wt_fault_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault =
                no_lines_fault)
    
    #----------------------------------------
    # DAE matrix or DAE vars
    #----------------------------------------
    
    model_mass_matrix =
        DAE_MassMatrix(
            length(state_labels),
            length(algebraic_vars_labels) )
    
    model_bool_dae_vars =
        DAE_BoolVector(
            length(state_labels),
            length(algebraic_vars_labels) )
        
    ode_gens_mass_matrix =
        DAE_MassMatrix(
            length(state_labels),
            0 )
    
    ode_gens_bool_dae_vars =
        DAE_BoolVector(
            length(state_labels),
            0 )
    
    #----------------------------------------
    # Parameters indices 
    #----------------------------------------
    
    Pg_Qg_Png_Qng_Pll_Qll_Idx =
        get_generic_Pg_Qg_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)

    #--------------------------------------
    
    scale_Pg_Qg_Png_Qng_Pll_Qll_Idx =
        get_generic_scale_Pg_Qg_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)
    
    #--------------------------------------
    
    Png_Qng_Pll_Qll_Idx =
        get_generic_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)

    #--------------------------------------
    
    Pg_Png_Qng_Idx =
        get_generic_Pg_Png_Qng_Idx(
            dyn_pf_fun_kwd_net_idxs)

    #--------------------------------------

    scale_Pg_Png_Qng_Idx =
        get_generic_scale_Pg_Png_Qng_Idx(
    dyn_pf_fun_kwd_net_idxs)

    #--------------------------------------

    dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx =
        get_dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)

    #--------------------------------------

    dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx =
        get_dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)
    
    #--------------------------------------

    dyn_vh_id_iq_V_ref_Tm_Idx =
        get_dyn_vh_id_iq_V_ref_Tm_Idx(
            gens_nodes_idx;
            reverse_idx = false)
    
    dyn_Tm_V_ref_id_iq_vh_Idx =
        get_dyn_vh_id_iq_V_ref_Tm_Idx(
            gens_nodes_idx;
            reverse_idx = true)


    dyn_V_ref_Tm_id_iq_vh_Idx =
        get_dyn_V_ref_Tm_vh_id_iq_Idx(
            gens_nodes_idx)

    #--------------------------------------

    dyn_vh_id_iq_ωref0_vref0_porder0_Idx =
        get_dyn_vh_id_iq_ωref0_vref0_porder0_Idx(
            gens_nodes_idx;
            reverse_idx = false)


    dyn_porder0_vref0_ωref0_id_iq_vh_Idx =
        get_dyn_vh_id_iq_ωref0_vref0_porder0_Idx(
            gens_nodes_idx;
            reverse_idx = true)
    
    dyn_ωref0_vref0_porder0_id_iq_vh_Idx =
        get_dyn_ωref0_vref0_porder0_id_iq_vh_Idx(
            gens_nodes_idx)

    ωref0_vref0_porder0_id_iq_vh_Idx =
        get_ωref0_vref0_porder0_id_iq_vh_Idx(
            gens_nodes_idx)

    #--------------------------------------
    
    id_iq_pg_vh_Idx =
        get_id_iq_pg_vh_Idx(
            gens_nodes_idx)
    
    #--------------------------------------
    
    dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx =
        get_dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)
    
    #--------------------------------------

    dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx =
        get_dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs)

    #--------------------------------------
    
    pf_vh_θh_idx_and_idx2Idx =
        get_pf_vh_θh_idx_and_idx2Idx(
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs)
    
    #--------------------------------------

    dyn_pf_flat_vh_flat_θh_Idx =
        get_generic_flat_vh_flat_θh_Idx(
            gens_nodes_idx,
            all_nodes_idx)

    #--------------------------------------
    
    dyn_pf_flat_vh_flat_θh_id_iq_Idx =
        get_generic_flat_vh_flat_θh_id_iq_Idx(
            gens_nodes_idx,
            all_nodes_idx)

    #--------------------------------------

    dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx =
        get_generic_flat_vh_flat_θh_wt_slack_value_Idx(
            all_nodes_idx)
    
    #--------------------------------------

    dyn_pf_vh_θh_id_iq_vhf_θhf_Idx =
        get_generic_vh_θh_id_iq_vhf_θhf_Idx(
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault = no_lines_fault)
       
    #--------------------------------------

    dyn_pf_vh_vhf_θh_θhf_id_iq_Idx =
         get_generic_vh_vhf_θh_θhf_id_iq_Idx(
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault = no_lines_fault)
    
    #--------------------------------------
    
    dyn_pf_vh_vhf_Idx =
        get_generic_vh_vhf_Idx(
            all_nodes_idx;
            no_lines_fault = no_lines_fault)
    
    #--------------------------------------

    dyn_pf_θh_θhf_Idx =
        get_generic_θh_θhf_Idx(
            all_nodes_idx;
            no_lines_fault = no_lines_fault)
    
    #----------------------------------------
    #----------------------------------------
    
    dyn_3_gens_type_paras_Idx =
        get_generic_n_net_comp_type_paras_Idx(
            gens_nodes_idx;
            no_paras = 3,
            comp_type = "gens")
    
    #--------------------------------------

    dyn_5_gens_type_paras_Idx =
        get_generic_n_net_comp_type_paras_Idx(
            gens_nodes_idx;
            no_paras = 5,
            comp_type = "gens")
    
    #--------------------------------------

    dyn_10_gens_type_paras_Idx =
        get_generic_n_net_comp_type_paras_Idx(
            gens_nodes_idx;
            no_paras = 10,
            comp_type = "gens")
    
    #--------------------------------------

    dyn_10_non_gens_type_paras_Idx =
        get_generic_n_net_comp_type_paras_Idx(
            non_gens_nodes_idx;
            no_paras = 10,
            comp_type = "non_gens")

    #--------------------------------------

    dyn_10_all_nodes_type_paras_Idx =
        get_generic_n_net_comp_type_paras_Idx(
            all_nodes_idx;
            no_paras = 10,
            comp_type = "all_nodes")
    
    #--------------------------------------
    
    dyn_2_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx =
        get_generic_n_gens_paras_wt_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs;
            no_gens_paras = 2)
    
    #--------------------------------------
    
    dyn_3_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx =
        get_generic_n_gens_paras_wt_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs;
            no_gens_paras = 3)

    #--------------------------------------
    
    dyn_4_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx =
        get_generic_n_gens_paras_wt_Png_Qng_Pll_Qll_Idx(
            dyn_pf_fun_kwd_net_idxs;
            no_gens_paras = 4)

    #----------------------------------------
    #----------------------------------------
    
    system_states_idx_kwd_para =
        (;
         state_vars_idx,
         vec_comp_states_Idx, 
         state_vars_and_i_dq_Idx_in_state,
         state_algebraic_vars_Idx_in_state,
         gens_state_vars_idx_in_state)
    
    #--------------------------------------

    
    system_paras_idx_kwd_para =
        (;
         dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx)
        
    #--------------------------------------

    plants_dyn_fun_idx_kwd_para =
        (;
         state_vars_idx,
         vec_comp_states_Idx,

         dyn_vh_id_iq_V_ref_Tm_Idx,
         dyn_V_ref_Tm_id_iq_vh_Idx,    
         dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
         dyn_vh_id_iq_ωref0_vref0_porder0_Idx,
         dyn_5_gens_type_paras_Idx) 

    #--------------------------------------

    plants_algeb_fun_idx_kwd_para =
        (;
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         dyn_2_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
         dyn_3_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx
         ) 
    
    #--------------------------------------
    
    
    return (;
            gens_govs_avrs_states_syms,
            gens_govs_avrs_types,
            
            vec_vec_gens_govs_avrs_states_syms,
            state_vars_idx,
            vec_comp_states_Idx,                        

            plants_states_syms,
            generic_state_sym,
            
            state_labels,
            algebraic_vars_labels,
            network_vars_labels,

            list_gens_type_syms,
            list_govs_type_syms,
            list_avrs_type_syms,
            
            gens_cb_paras_func,            
            gens_init_func,
            gens_output_func,
            gens_dyn_func,
            ode_gens_dyn_func,
            dae_gens_dyn_func,

            govs_cb_paras_func,            
            govs_init_func,
            govs_output_func,
            govs_dyn_func,
            ode_govs_dyn_func,
            dae_govs_dyn_func,

            avrs_cb_paras_func,            
            avrs_init_func,
            avrs_output_func,
            avrs_dyn_func,
            ode_avrs_dyn_func,
            dae_avrs_dyn_func,
            
            comps_callback_paras_funs,
            comps_init_funs,
            comps_output_funs,
            comps_dyn_funs,
            ode_comps_dyn_funs,
            dae_comps_dyn_funs,

            algebraic_state_sym,
            model_syms,
            
            nodes_names,            
            gens_nodes_names,
            non_gens_nodes_names,

            SM_gens_nodes_names,
            SC_gens_nodes_names,

            generic_model_states_comp_idxs_in_Idx,
            
            generic_model_vars_wt_i_dq_Idx_in_state,
            generic_model_vars_no_i_dq_Idx_in_state,
            
            gens_state_vars_idx_in_state,
            state_vars_and_i_dq_Idx_in_state,

            state_vars_and_i_dq_wt_fault_Idx_in_state,
            
            state_algebraic_vars_Idx_in_state,
            state_algebraic_vars_wt_fault_Idx_in_state,
            
            model_mass_matrix,
            model_bool_dae_vars,
            
            ode_gens_mass_matrix,
            ode_gens_bool_dae_vars,

            Pg_Qg_Png_Qng_Pll_Qll_Idx,
            Png_Qng_Pll_Qll_Idx,

            Pg_Png_Qng_Idx,

            scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
            scale_Pg_Png_Qng_Idx,

            dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            
            dyn_vh_id_iq_V_ref_Tm_Idx,
            dyn_V_ref_Tm_id_iq_vh_Idx,

            id_iq_pg_vh_Idx,
            
            dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
            dyn_vh_id_iq_ωref0_vref0_porder0_Idx,

            ωref0_vref0_porder0_id_iq_vh_Idx,
            
            dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
            dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

            pf_vh_θh_idx_and_idx2Idx,
            
            dyn_pf_flat_vh_flat_θh_Idx,
            
            dyn_pf_flat_vh_flat_θh_id_iq_Idx,

            dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,
            
            dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
            dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
            
            dyn_pf_vh_vhf_Idx,
            dyn_pf_θh_θhf_Idx,

            dyn_3_gens_type_paras_Idx,
            dyn_5_gens_type_paras_Idx,

            dyn_10_gens_type_paras_Idx,
            dyn_10_non_gens_type_paras_Idx,
            dyn_10_all_nodes_type_paras_Idx,

            dyn_2_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
            dyn_3_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
            dyn_4_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
            
            system_states_idx_kwd_para,
            system_paras_idx_kwd_para,
            plants_dyn_fun_idx_kwd_para,
            plants_algeb_fun_idx_kwd_para )
    
end


#-----------------------------------------------------
#-----------------------------------------------------

# comment
