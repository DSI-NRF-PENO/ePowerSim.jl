# (c) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.4:  pg 135 - 136, 6.121 - 6.123

#########################################################


# ---------------------------------------------------
# ---------------------------------------------------
# unimplemented
# ---------------------------------------------------
# ---------------------------------------------------

function add_new_node_k_to_neighbourhood_of_node_j()
    
end

function remove_node_k_from_neighbourhood_of_node_j()
    
end


function insert_node_k_in_neighbourhood_list_of_node_j()
    
end

function delete_node_k_from_neighbourhood_list_of_node_j()
    
end


function insert_k_to_j_link_admittance_in_node_j_row_in_Ynet()
    
end

function delete_k_to_j_link_admittance_in_node_j_row_in_Ynet()
    
end



function record_a_faulty_line_nodes_association_idx!(
    dict_nodes_neighbourhood_tuple_list,
    nodes_idx_with_adjacent_nodes_idx,
    faulted_line_k_j,
    node_j,
    node_k)

    @assert node_k ∈ faulted_line_k_j
    
    @assert node_j ∈ faulted_line_k_j

    if node_j ∉ keys( dict_nodes_neighbourhood_tuple_list )

        dict_nodes_neighbourhood_tuple_list[ node_j  ] = []
        
    end

    if node_k ∉ keys( dict_nodes_neighbourhood_tuple_list )

        dict_nodes_neighbourhood_tuple_list[ node_k  ] = []
        
    end
    

    push!( dict_nodes_neighbourhood_tuple_list[ node_k ],
        get_node_a_idx_in_node_b_row_in_Ynet(
            nodes_idx_with_adjacent_nodes_idx,
            node_j,
            node_k ))

    push!(dict_nodes_neighbourhood_tuple_list[ node_j ],
        get_node_a_idx_in_node_b_row_in_Ynet(
            nodes_idx_with_adjacent_nodes_idx,
            node_k,
            node_j ))

    return nothing    
    
end


# ---------------------------------------------------
# ---------------------------------------------------

function get_node_a_idx_in_node_b_row_in_Ynet(
    nodes_idx_with_adjacent_nodes_idx,
    node_a,
    node_b )

   return (node_a,           
           findfirst(
               x -> x == node_a,
               nodes_idx_with_adjacent_nodes_idx[
                   node_b ]))

end


# ---------------------------------------------------
# ---------------------------------------------------
# Faults
# ---------------------------------------------------
# ---------------------------------------------------

# https://discourse.julialang.org/t/changing-size-of-ode-system-resize-works-but-deleteat-doesnt/45900

# https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/

# https://discourse.julialang.org/t/differentialequations-jl-vectorcontinuouscallback-vs-callbackset/88536/4 

# sol = solve(prob, Tsit5(), callback = cb, tstops = [4.0])

# ---------------------------------------------------
# callbacks
# ---------------------------------------------------


"""

no_fault = 0

on_fault = 1

clear_fault = 2

system_fault_status[1] = no_fault

on_fault_time = 4

clear_fault_time = 6

"""

#----------------------------------------

function on_fault_affect!(
    integrator,
    no_lines_fault )

    integrator.p.system_fault_status[1] = 1
    
    # The rationale for twice: vh and θh
    
    resize!(integrator,
            length(integrator.u) +
                2 * no_lines_fault )
    nothing

end


function clear_fault_affect!(
    integrator,
    no_cleared_lines_fault )

    integrator.p.system_fault_status[1] = 2
    
    # The rationale for twice: vh and θh
    
    resize!(integrator,
            length( integrator.u) -
                2 * no_cleared_lines_fault)

    nothing
    

end


function partial_clear_fault_condition(
    u, t, integrator,
    partial_clear_fault_time)
    
    integrator.t == partial_clear_fault_time

end



function partial_clear_fault_affect!(
    integrator,
    no_cleared_lines_fault )

    integrator.p.system_fault_status[1] = 3
    
    # The rationale for twice: vh and θh
    
    resize!(integrator,
            length( integrator.u) -
                2 * no_cleared_lines_fault)

    nothing 
    

end

#----------------------------------------


function on_fault_affect!(
    integrator)

    integrator.p.system_fault_status[1] = 1

end


function clear_fault_affect!(
    integrator )

    integrator.p.system_fault_status[1] = 2
    

end


function partial_clear_fault_affect!(
    integrator )

    integrator.p.system_fault_status[1] = 3
    

end

#----------------------------------------
#----------------------------------------

function clear_fault_wt_state_and_model_para_affect!(
    integrator,
    model_para,
    states )

    integrator.p .=
        model_para

    integrator.u .=
        states

    @show "clear affect wt state affect triggered at time $(integrator.t)"
    
    nothing
    

end


function clear_fault_wt_state_and_model_para_affect!(
    integrator,
    model_para,
    states,
    para_idx,
    state_idx )

    integrator.p[ para_idx ] .=
        model_para[ para_idx ]

    integrator.u[ state_idx ] .=
        states[ state_idx ]

        @show "clear wt state idx triggered at time $(integrator.t)"
    nothing
    
end


function on_pre_fault_condition(
    u, t, integrator,
    on_fault_time)
    
    integrator.t ==  on_fault_time - 0.02 
    
end


function pre_fault_affect!(
    integrator)

    # add_tstop!(integrator, on_fault_time)

    proposed_dt = get_proposed_dt(integrator)
    @show proposed_dt, "$(proposed_dt) "

    # set_proposed_dt!(integrator, 0.001)

    nothing
end


function on_fault_condition(
    u, t, integrator,
    on_fault_time)
    
    integrator.t == on_fault_time
    
end


function on_fault_wt_model_dynamics_para_affect!(
    integrator,
    pe_order_para,
    idx )

    integrator.p[idx] .=
        pe_order_para

    nothing

end



function on_fault_wt_model_dynamics_para_affect!(
    integrator,
    model_para )

    integrator.p .=
        model_para
    
    @show "on fault affect triggered at time $(integrator.t)"

    nothing

end


function clear_fault_condition(
    u, t, integrator,
    clear_fault_time)
    
    integrator.t ==
        clear_fault_time

end


function clear_fault_wt_model_dynamics_para_affect!(
    integrator,
    pe_order_para,
    idx )

    integrator.p[idx] .=
        pe_order_para

    nothing
    

end


function clear_fault_wt_model_dynamics_para_affect!(
    integrator,
    model_para )


    @show "clear affect triggered at time $(integrator.t)"
    
    integrator.p .=
        model_para

    nothing
    

end



function on_fault_Ynet_affect!(
    integrator)

    integrator.p.system_fault_status[1] = 1

end


function clear_fault_Ynet_affect!(
    integrator )

    integrator.p.system_fault_status[1] = 2
    
    nothing
end


#----------------------------------------
# fault on a line
#----------------------------------------

#----------------------------------------
# load drop
#----------------------------------------


#----------------------------------------
# line outage
#----------------------------------------


function on_line_outage_condition(
    u, t, integrator,
    on_line_outage_time)
    
    integrator.t == on_line_outage_time
    
end


function on_line_outage_affect!(
    integrator )

    integrator.p.system_fault_status[1] = 2
    
    nothing
end

function on_generation_adjustment_condition(
    u, t, integrator,
    on_generation_adjustment_time)
    
    integrator.t == on_generation_adjustment_time
    
end


function on_generation_adjustment_affect!(
    integrator,
    ref_adj,
    idx)
    
    integrator.p.generic_model_dynamics_para[idx] .=
        ref_adj

    nothing    
end




#----------------------------------------
#----------------------------------------

function cb_fun_make_fault_callbacks(
    on_fault_time,
    clear_fault_time,
    on_fault_model_para,
    clear_fault_model_para,
    clear_fault_states;
    use_state_in_on_clear_fault =
        false,
    para_idx =
        nothing,
    state_idx =
        nothing)
    

    # DiscreteCallback
    
    if use_state_in_on_clear_fault == false
        
       return  CallbackSet(
    DiscreteCallback(
        (u, t, integrator) ->
            on_fault_condition(
                u, t, integrator,
                on_fault_time ),
        (integrator) ->
            on_fault_wt_model_dynamics_para_affect!(
                integrator,
                on_fault_model_para ); 
        save_positions=(true, true) ),
    
    DiscreteCallback(
        (u, t, integrator) ->
            clear_fault_condition(
                u, t, integrator,
                clear_fault_time ),

        (integrator) ->
            clear_fault_wt_model_dynamics_para_affect!(
                integrator,
                clear_fault_model_para);
        save_positions=(false, true) ))
        
    else
            return CallbackSet(
                DiscreteCallback(
                    (u, t, integrator) ->
            on_fault_condition(
                u, t, integrator,
                on_fault_time ),

                    (integrator) ->
                        on_fault_wt_model_dynamics_para_affect!(
                            integrator, on_fault_model_para); 
                    save_positions=(true,true) ),

                DiscreteCallback(
                    (u, t, integrator) ->
                        clear_fault_condition(
                            u, t,
                            integrator,
                            clear_fault_time ),

                    (integrator) ->
                        clear_fault_wt_state_and_model_para_affect!(
                            integrator, clear_fault_model_para,
            clear_fault_states);
                    save_positions=(false,true) ))
    end
    
        
end

# ---------------------------------------------------

function make_lines_faults_data_set(
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    
    all_nodes_idx,
    n2s_all_nodes_idx,
    
    list_faulted_line_a_b_orientation,
    list_fault_point_from_node_a,
    list_fault_resistance,
    list_no_line_circuit  )



    """
    a_node_fault_event: A node short circuit fault event
    is like adding an infinite shunt admittance to the node.

    a_line_fault_event: Requires an introduction of a new
    node called fault node. The admittance if the node in
    relation to the two nodes at the end of the line will
    depend on there the fault occurs.

    consider an edge list_faulted_line_a_b_orientation
    with nodes (a, b), where `a` is the from node and `b`
    is the to node.

    If we denote the fault node by nf, (a, nf) and (nf, b)

    Let us assign an index after the highest index of the
    network to the fault node. i.e

    """

    
    faulty_Ynet = deepcopy(Ynet)
    
    faulty_nodes_idx_with_adjacent_nodes_idx =
        deepcopy(nodes_idx_with_adjacent_nodes_idx)

    faulty_all_nodes_idx =
        deepcopy(all_nodes_idx )
    
    n2s_faulty_all_nodes_idx =
        deepcopy(n2s_all_nodes_idx)

    list_faulted_line_a_b_orientation =
        deepcopy(list_faulted_line_a_b_orientation)

    list_fault_point_from_node_a =
        deepcopy(list_fault_point_from_node_a)
    
    list_fault_resistance =
        deepcopy(list_fault_resistance)
    
    list_no_line_circuit =
        deepcopy(list_no_line_circuit)

    # -----------------------------
    
    fault_nodes_idx = []

    n2s_fault_nodes_idx =
        OrderedDict{Union{Int64, String, Symbol}, Int64}()

    list_faulty_line_Yl = []

    list_healthy_lines_Yl = []

    list_Ya_nkf = []

    list_Ynkf_b = []

    list_node_b_idx_in_a_node_row = []

    list_node_a_idx_in_b_node_row = []
    
    #----------------------------------------
    
    for (faulted_line_a_b,
         fault_point_from_node_a,
         fault_resistance,
         no_line_circuit) in
        zip(list_faulted_line_a_b_orientation,
            list_fault_point_from_node_a,
            list_fault_resistance,
            list_no_line_circuit)

        # faulted_line_a_b = (a, b) = 5, 4

        # @show faulted_line_a_b

        node_a_idx =
            first(faulted_line_a_b)

        node_b_idx =
            last(faulted_line_a_b)

        nf = n2s_faulty_all_nodes_idx[
            max(faulty_all_nodes_idx...) ] + 1

        push!( faulty_all_nodes_idx, nf)

        n2s_faulty_all_nodes_idx[ nf] = nf

        push!(fault_nodes_idx, nf)

        n2s_fault_nodes_idx[ nf] =
            values(n2s_fault_nodes_idx) == nothing ?
            1 : max(0,
                    collect(values(
                        n2s_fault_nodes_idx))...) + 1

       #----------------------------------------

        # gens_nodes_idx, n2s_gens_idx

        #----------------------------------------
        
        push!(faulty_nodes_idx_with_adjacent_nodes_idx,
                 [nf, node_a_idx, node_b_idx] )

       #----------------------------------------


        node_b_idx_in_a_node_row =
           findfirst(
               x -> x == node_b_idx,
               faulty_nodes_idx_with_adjacent_nodes_idx[
                   node_a_idx ])

        node_a_idx_in_b_node_row  =
           findfirst(
               x -> x == node_a_idx,
               faulty_nodes_idx_with_adjacent_nodes_idx[
                   node_b_idx ])

        push!( list_node_b_idx_in_a_node_row,
               node_b_idx_in_a_node_row )

        push!( list_node_a_idx_in_b_node_row,
               node_a_idx_in_b_node_row )

        #----------------------------------------
        
        # Note, node_a or node_b be could be used

        Yl = faulty_Ynet[
            node_a_idx ][node_b_idx_in_a_node_row]

        #----------------------------------------

       # fault_point_from_node_a = 0.3

        faulty_line_Yl = Yl / no_line_circuit

        healthy_lines_Yl = Yl - Yl / no_line_circuit

        faulty_line_Zl = 1/faulty_line_Yl
        
        Za_nkf = faulty_line_Zl * fault_point_from_node_a

        Znkf_b = faulty_line_Zl * (1 -
            fault_point_from_node_a)
        
        Ya_nkf = 1/Za_nkf

        Ynkf_b = 1/Znkf_b
            
       # Ya_nkf =
       #     faulty_line_Yl/fault_point_from_node_a

       # Ynkf_b =
       #     faulty_line_Yl/(1 - fault_point_from_node_a)

        push!(list_Ya_nkf, Ya_nkf)
        
        push!(list_Ynkf_b, Ynkf_b)

       #----------------------------------------

       fault_admittance = 1/fault_resistance

       # nkf_row_in_Ynet =
       #     [-faulty_line_Yl - fault_admittance,
       #      Ya_nkf, Ynkf_b]

       nkf_row_in_Ynet =
           [-(Ya_nkf + Ynkf_b) - fault_admittance,
            Ya_nkf, Ynkf_b]
        
       push!(faulty_Ynet, nkf_row_in_Ynet)

        #----------------------------------------


        push!( list_faulty_line_Yl,
               faulty_line_Yl )

        push!( list_healthy_lines_Yl,
               healthy_lines_Yl )


       if no_line_circuit == 1

           #----------------------------------------

           # Replace  node_b in the neighbouhood of node_a
           # with nf, similarly replace node_a  in the
           # neighbouhood of node_b with nf

            faulty_nodes_idx_with_adjacent_nodes_idx[
                node_a_idx][node_b_idx_in_a_node_row] = nf

            faulty_nodes_idx_with_adjacent_nodes_idx[
                node_b_idx][node_a_idx_in_b_node_row] = nf

           #----------------------------------------

           # Since the sum of a row should be zero,
           # remove the contribution of Yl to the first
           # element and add the contribution of the
           # new segment

            faulty_Ynet[ node_a_idx][ 1] =
                    faulty_Ynet[ node_a_idx][ 1] +
                    faulty_line_Yl - Ya_nkf

            faulty_Ynet[ node_a_idx][
                    node_b_idx_in_a_node_row] = Ya_nkf

           # Since the sum of a row should be zero,
           # remove the contribution of Yl to the
           # first element and add the contribution of
           # the new segment

           faulty_Ynet[ node_b_idx][ 1] =
               faulty_Ynet[ node_b_idx][ 1] +
               faulty_line_Yl - Ynkf_b


            faulty_Ynet[ node_b_idx][
                    node_a_idx_in_b_node_row] = Ynkf_b

       else

           #----------------------------------------

           # Add the fault node nf to the neighbouhood of
           # node_a and node_b

           push!(
               faulty_nodes_idx_with_adjacent_nodes_idx[
                node_a_idx], nf)

           push!(
               faulty_nodes_idx_with_adjacent_nodes_idx[
                node_b_idx], nf)

           #----------------------------------------

           # Since the sum of a row should be zero,
           # remove the contribution of Yl to the first
           # element and add the contribution of the
           # new segment

           faulty_Ynet[ node_a_idx][ 1] =
               faulty_Ynet[ node_a_idx][ 1] +
               faulty_line_Yl - Ya_nkf

           # since nf is expected to be the last node
           # push Ya_nkf into the list
           
            push!(faulty_Ynet[ node_a_idx], Ya_nkf)

           # Since the sum of a row should be zero,
           # remove the contribution of Yl to the first
           # element and add the contribution of the
           # new segment

           faulty_Ynet[ node_b_idx][ 1] =
               faulty_Ynet[ node_b_idx][ 1] +
               faulty_line_Yl - Ynkf_b

           # since nf is expected to be the last node
           # push Ynkf_b into the list
           
           push!(faulty_Ynet[ node_b_idx], Ynkf_b)

       end

        
    end

    return (;faulty_Ynet,
            faulty_nodes_idx_with_adjacent_nodes_idx,
            
            faulty_all_nodes_idx,
            n2s_faulty_all_nodes_idx,

            fault_nodes_idx,
            n2s_fault_nodes_idx,            

            list_Ya_nkf,
            list_Ynkf_b,
            
            list_faulty_line_Yl,
            list_healthy_lines_Yl,

            list_node_b_idx_in_a_node_row,
            list_node_a_idx_in_b_node_row,

            list_faulted_line_a_b_orientation,
            list_fault_point_from_node_a,
            list_fault_resistance,
            list_no_line_circuit)

    
end


function add_lines_faults_data_set(
    faulty_Ynet,
    faulty_nodes_idx_with_adjacent_nodes_idx,
    
    faulty_all_nodes_idx,
    n2s_faulty_all_nodes_idx,
    
    fault_nodes_idx,
    n2s_fault_nodes_idx,
    
    list_Ya_nkf,
    list_Ynkf_b,
    
    list_faulty_line_Yl,
    list_healthy_lines_Yl,
    
    list_node_b_idx_in_a_node_row,
    list_node_a_idx_in_b_node_row,
    
    list_faulted_line_a_b_orientation,
    list_fault_point_from_node_a,
    list_fault_resistance,
    list_no_line_circuit )
    
    #----------------------------------------
    
    for (faulted_line_a_b,
         fault_point_from_node_a,
         fault_resistance,
         no_line_circuit) in
        zip(list_faulted_line_a_b_orientation,
            list_fault_point_from_node_a,
            list_fault_resistance,
            list_no_line_circuit)

        # faulted_line_a_b = (a, b) = 5, 4

        node_a_idx = first( faulted_line_a_b)

        node_b_idx = last( faulted_line_a_b)

        #--------------------------------------
        # Check for an existing fault entry
        #--------------------------------------

        fault_idx_with_adjacent_nodes_idx =
            faulty_nodes_idx_with_adjacent_nodes_idx[
                collect( keys( n2s_fault_nodes_idx)) ]

        # Each fault_idx_with_adjacent_nodes_idx has only
        # three entries, fault node, node a and node b
        
        fault_idx_with_adjacent_nodes_a_idx =
            second.(fault_idx_with_adjacent_nodes_idx)
        
        fault_idx_with_adjacent_nodes_b_idx =
            third.(fault_idx_with_adjacent_nodes_idx)

        a_idx = findfirst(
            x -> x == node_a_idx,
            fault_idx_with_adjacent_nodes_a_idx )


        b_idx = findfirst(
            x -> x == node_b_idx,
            fault_idx_with_adjacent_nodes_b_idx )

        #--------------------------------------

        # a_idx == b_idx means there is already an entry
        
        if a_idx != b_idx

            nf = n2s_faulty_all_nodes_idx[
                max( faulty_all_nodes_idx...) ] + 1

            push!( faulty_all_nodes_idx, nf)

            n2s_faulty_all_nodes_idx[ nf] = nf

            push!( fault_nodes_idx, nf)

            n2s_fault_nodes_idx[ nf] =
                values( n2s_fault_nodes_idx) == nothing ?
                1 :
                max(0,
                    collect(values(
                        n2s_fault_nodes_idx))...) + 1

           #----------------------------------------

            push!( faulty_nodes_idx_with_adjacent_nodes_idx,
                     [nf, node_a_idx, node_b_idx] )

           #----------------------------------------

           node_b_idx_in_a_node_row =
               findfirst(
                   x -> x == node_b_idx,
                   faulty_nodes_idx_with_adjacent_nodes_idx[
                       node_a_idx ])

            node_a_idx_in_b_node_row  =
               findfirst(
                   x -> x == node_a_idx,
                   faulty_nodes_idx_with_adjacent_nodes_idx[
                       node_b_idx ])

            push!( list_node_b_idx_in_a_node_row,
                   node_b_idx_in_a_node_row )

            push!( list_node_a_idx_in_b_node_row,
                   node_a_idx_in_b_node_row )

            #----------------------------------------

            # Note, node_a or node_b be could be used

            Yl = faulty_Ynet[
                node_a_idx ][ node_b_idx_in_a_node_row]

            #----------------------------------------

           # fault_point_from_node_a = 0.3

           faulty_line_Yl = Yl / no_line_circuit

           healthy_lines_Yl = Yl - Yl / no_line_circuit

           Ya_nkf =
               faulty_line_Yl/fault_point_from_node_a

           Ynkf_b =
               faulty_line_Yl/(1 - fault_point_from_node_a)

            push!( list_Ya_nkf, Ya_nkf)

            push!( list_Ynkf_b, Ynkf_b)

           #----------------------------------------

           fault_admittance = 1/fault_resistance

           nkf_row_in_Ynet =
               [ -faulty_line_Yl - fault_admittance,
                Ya_nkf, Ynkf_b]

           push!( faulty_Ynet, nkf_row_in_Ynet)

            #----------------------------------------

            push!( list_faulty_line_Yl,
                   faulty_line_Yl )

            push!( list_healthy_lines_Yl,
                   healthy_lines_Yl )

            #----------------------------------------


           if no_line_circuit == 1

               #----------------------------------------

               # Replace  node_b in the neighbouhood of
               # node_a with nf, similarly replace node_a
               # in the neighbouhood of node_b with nf

                faulty_nodes_idx_with_adjacent_nodes_idx[
                    node_a_idx][
                        node_b_idx_in_a_node_row] = nf

                faulty_nodes_idx_with_adjacent_nodes_idx[
                    node_b_idx][
                        node_a_idx_in_b_node_row] = nf

               #----------------------------------------

               # Since the sum of a row should be zero,
               # remove the contribution of Yl to the first
               # element and add the contribution of the
               # new segment

                faulty_Ynet[ node_a_idx][ 1] =
                        faulty_Ynet[ node_a_idx][ 1] +
                        faulty_line_Yl - Ya_nkf

                faulty_Ynet[
                    node_a_idx][
                        node_b_idx_in_a_node_row] = Ya_nkf

               # Since the sum of a row should be zero,
               # remove the contribution of Yl to the
               # first element and add the contribution of
               # the new segment

               faulty_Ynet[ node_b_idx][ 1] =
                   faulty_Ynet[ node_b_idx][ 1] +
                   faulty_line_Yl - Ynkf_b


                faulty_Ynet[ node_b_idx][
                        node_a_idx_in_b_node_row] = Ynkf_b

           else

               #----------------------------------------

               # Add the fault node nf to the neighbouhood
               # of node_a and node_b

               push!(
                   faulty_nodes_idx_with_adjacent_nodes_idx[
                       node_a_idx], nf)

               push!(
                   faulty_nodes_idx_with_adjacent_nodes_idx[
                       node_b_idx], nf)

               #----------------------------------------

               # Since the sum of a row should be zero,
               # remove the contribution of Yl to the first
               # element and add the contribution of the
               # new segment

               faulty_Ynet[ node_a_idx][ 1] =
                   faulty_Ynet[ node_a_idx][ 1] +
                   faulty_line_Yl - Ya_nkf

                push!(faulty_Ynet[ node_a_idx], Ya_nkf)

               # Since the sum of a row should be zero,
               # remove the contribution of Yl to the first
               # element and add the contribution of the
               # new segment

               faulty_Ynet[ node_b_idx][ 1] =
                   faulty_Ynet[ node_b_idx][ 1] +
                   faulty_line_Yl - Ynkf_b

                push!(faulty_Ynet[ node_b_idx], Ynkf_b)

           end

        end
        
    end

    return (;faulty_Ynet,
            faulty_nodes_idx_with_adjacent_nodes_idx,
            
            faulty_all_nodes_idx,
            n2s_faulty_all_nodes_idx,

            fault_nodes_idx,
            n2s_fault_nodes_idx,
            
            list_Ya_nkf,
            list_Ynkf_b,
            
            list_faulty_line_Yl,
            list_healthy_lines_Yl,

            list_node_b_idx_in_a_node_row,
            list_node_a_idx_in_b_node_row)
    
end


function remove_lines_faults_data_set(
    clear_fault_selection_list;
    faulty_Ynet =
        faulty_Ynet,
    faulty_nodes_idx_with_adjacent_nodes_idx =
        faulty_nodes_idx_with_adjacent_nodes_idx,
    
    faulty_all_nodes_idx =
        faulty_all_nodes_idx,
    
    n2s_faulty_all_nodes_idx =
        n2s_faulty_all_nodes_idx,
    
    fault_nodes_idx =
        fault_nodes_idx,
    
    n2s_fault_nodes_idx =
        n2s_fault_nodes_idx,
    
    list_Ya_nkf =
        list_Ya_nkf,
    
    list_Ynkf_b =
        list_Ynkf_b,
    
    list_faulty_line_Yl =
        list_faulty_line_Yl,
    
    list_healthy_lines_Yl =
        list_healthy_lines_Yl,
    
    list_node_b_idx_in_a_node_row =
        list_node_b_idx_in_a_node_row,
    
    list_node_a_idx_in_b_node_row =
        list_node_a_idx_in_b_node_row,

    list_faulted_line_a_b_orientation =
        list_faulted_line_a_b_orientation,
    
    list_fault_point_from_node_a =
        list_fault_point_from_node_a,
    
    list_fault_resistance =
        list_fault_resistance,
    
    list_no_line_circuit =
        list_no_line_circuit )

    sorted_selection_list = sort(
        clear_fault_selection_list)

    # # This is used to collect only exixsting faults
    # # in the selection list
    
    # exixting_fault_entries = []
    
    #----------------------------------------
    if fault_nodes_idx != []

        for (faulted_line_a_b,
             fault_point_from_node_a,
             fault_resistance,
             no_line_circuit,

             nf,

             Ya_nkf,
             Ynkf_b,

             faulty_line_Yl,
             healthy_lines_Yl,
             node_b_idx_in_a_node_row,
             node_a_idx_in_b_node_row) in
            zip(list_faulted_line_a_b_orientation[
                sorted_selection_list],

                list_fault_point_from_node_a[
                    sorted_selection_list],

                list_fault_resistance[
                    sorted_selection_list],

                list_no_line_circuit[
                    sorted_selection_list],

                fault_nodes_idx[
                    sorted_selection_list],

                list_Ya_nkf[
                    sorted_selection_list],

                list_Ynkf_b[
                    sorted_selection_list],

                list_faulty_line_Yl[
                    sorted_selection_list],

                list_healthy_lines_Yl[
                    sorted_selection_list],

                list_node_b_idx_in_a_node_row[
                    sorted_selection_list],

                list_node_a_idx_in_b_node_row[
                    sorted_selection_list])

            #--------------------------------------

            node_a_idx = first(faulted_line_a_b)

            node_b_idx = last(faulted_line_a_b)

            #--------------------------------------
            # Check for an existing fault entry
            #--------------------------------------

            fault_idx_with_adjacent_nodes_idx =
                faulty_nodes_idx_with_adjacent_nodes_idx[
                    collect(keys(n2s_fault_nodes_idx)) ]

            # Each fault_idx_with_adjacent_nodes_idx has only
            # three entries, fault node, node a and node b

            fault_idx_with_adjacent_nodes_a_idx =
                second.(fault_idx_with_adjacent_nodes_idx)

            fault_idx_with_adjacent_nodes_b_idx =
                third.(fault_idx_with_adjacent_nodes_idx)

            a_idx = findfirst(
                x -> x == node_a_idx,
                fault_idx_with_adjacent_nodes_a_idx )


            b_idx = findfirst(
                x -> x == node_b_idx,
                fault_idx_with_adjacent_nodes_b_idx )

            #--------------------------------------

            # a_idx == b_idx means there is an entry

            if a_idx == b_idx

                # push!(exixting_fault_entries, nf)

               if no_line_circuit == 1

                   #----------------------------------------

                   deleteat!(
                       faulty_nodes_idx_with_adjacent_nodes_idx[
                        node_a_idx], node_b_idx_in_a_node_row)

                   deleteat!(
                       faulty_nodes_idx_with_adjacent_nodes_idx[
                        node_b_idx], node_a_idx_in_b_node_row)

                   #----------------------------------------

                   # Since the sum of a row should be zero,
                   # remove the contribution of Yl to the first
                   # element and add the contribution of the
                   # new segment

                    faulty_Ynet[ node_a_idx][ 1] =
                            faulty_Ynet[ node_a_idx][ 1] +
                             Ya_nkf

                    deleteat!(faulty_Ynet[ node_a_idx],
                            node_b_idx_in_a_node_row)

                   # Since the sum of a row should be zero,
                   # remove the contribution of Yl to the
                   # first element and add the contribution of
                   # the new segment

                   faulty_Ynet[ node_b_idx][ 1] =
                       faulty_Ynet[ node_b_idx][ 1] + Ynkf_b

                    deleteat!(faulty_Ynet[ node_b_idx],
                            node_a_idx_in_b_node_row )
               else

                   #----------------------------------------

                   # Get the idx of fault node nf in a_node_row
                   # and b_node_row

                   node_nf_idx_in_a_node_row =
                       findfirst(
                           x -> x == nf,
                           faulty_nodes_idx_with_adjacent_nodes_idx[
                               node_a_idx ])

                    node_nf_idx_in_b_node_row  =
                       findfirst(
                           x -> x == nf,
                           faulty_nodes_idx_with_adjacent_nodes_idx[
                               node_b_idx ])

                   # Remove the fault node nf from the
                   # neighbouhood of node_a and node_b

                   deleteat!(
                       faulty_nodes_idx_with_adjacent_nodes_idx[
                        node_a_idx], node_nf_idx_in_a_node_row)

                   deleteat!(
                       faulty_nodes_idx_with_adjacent_nodes_idx[
                        node_b_idx], node_nf_idx_in_b_node_row)

                   #----------------------------------------

                   # Since the sum of a row should be zero,
                   # remove the contribution of the
                   # new segment

                   faulty_Ynet[ node_a_idx][ 1] =
                       faulty_Ynet[ node_a_idx][ 1] + Ya_nkf

                   deleteat!(faulty_Ynet[ node_a_idx],
                             node_nf_idx_in_a_node_row)

                   # Since the sum of a row should be zero,
                   # remove the contribution of the
                   # new segment

                   faulty_Ynet[ node_b_idx][ 1] =
                       faulty_Ynet[ node_b_idx][ 1] + Ynkf_b

                   deleteat!(faulty_Ynet[ node_b_idx],
                             node_nf_idx_in_b_node_row)

               end

            end

        end

    end
    
    if fault_nodes_idx != []


        #----------------------------------------
        # exixting_fault_entries

        deleteat!(faulty_Ynet,
                  fault_nodes_idx[
                      sorted_selection_list])

        deleteat!(faulty_nodes_idx_with_adjacent_nodes_idx,
                  fault_nodes_idx[
                      sorted_selection_list])

        #----------------------------------------

        # faulty_all_nodes_idx =
        #     setdiff(faulty_all_nodes_idx,
        #             fault_nodes_idx[
        #                 sorted_selection_list])


        for nf in fault_nodes_idx[sorted_selection_list]

            delete!(n2s_faulty_all_nodes_idx, nf)

            delete!(n2s_fault_nodes_idx, nf)


            a_fault_node_idx =
                findfirst(
                    x -> x == nf, faulty_all_nodes_idx )

            deleteat!(faulty_all_nodes_idx,
                      a_fault_node_idx)

        end

        # remake n2s_faulty_all_nodes_idx and
        # n2s_fault_nodes_idx

        keys_n2s_faulty_all_nodes =
            sort(collect(keys(n2s_faulty_all_nodes_idx)))

        n2s_faulty_all_nodes_idx =
            OrderedDict{Union{Int64, String, Symbol}, Int64}(
                a_key => a_value for (a_key, a_value) in
                    zip(keys_n2s_faulty_all_nodes,
                        1:length(keys_n2s_faulty_all_nodes)) )

        keys_n2s_fault_nodes =
            sort(collect(keys(n2s_fault_nodes_idx)))

        n2s_fault_nodes_idx =
            OrderedDict{Union{Int64, String, Symbol}, Int64}(
                a_key => a_value for (a_key, a_value) in
                    zip( keys_n2s_fault_nodes,
                        1:length( keys_n2s_fault_nodes)) )

        #----------------------------------------

        deleteat!(list_faulted_line_a_b_orientation,
                sorted_selection_list)

        deleteat!(list_fault_point_from_node_a,
                  sorted_selection_list)

        deleteat!(list_fault_resistance,
                  sorted_selection_list)

        deleteat!(list_no_line_circuit,
                  sorted_selection_list)

        deleteat!(fault_nodes_idx,
                  sorted_selection_list)

        deleteat!(list_Ya_nkf,
                  sorted_selection_list)

        deleteat!(list_Ynkf_b,
                  sorted_selection_list)

        deleteat!(list_faulty_line_Yl,
                  sorted_selection_list)

        deleteat!(list_healthy_lines_Yl,
                  sorted_selection_list)

        deleteat!(list_node_b_idx_in_a_node_row,
                  sorted_selection_list)

        deleteat!(list_node_a_idx_in_b_node_row,
                  sorted_selection_list)
        
    end

    return (;faulty_Ynet,
            faulty_nodes_idx_with_adjacent_nodes_idx,
            
            faulty_all_nodes_idx,
            n2s_faulty_all_nodes_idx,

            fault_nodes_idx,
            n2s_fault_nodes_idx,

            list_faulted_line_a_b_orientation,
            list_fault_point_from_node_a,
            list_fault_resistance,
            list_no_line_circuit,

            list_Ya_nkf,
            list_Ynkf_b,

            list_faulty_line_Yl,
            list_healthy_lines_Yl,

            list_node_b_idx_in_a_node_row,
            list_node_a_idx_in_b_node_row)
    
end


function get_cleared_selected_lines_faults_data_set(
    clear_fault_selection_list;
    faulty_Ynet =
        faulty_Ynet,
    faulty_nodes_idx_with_adjacent_nodes_idx =
        faulty_nodes_idx_with_adjacent_nodes_idx,
    
    faulty_all_nodes_idx =
        faulty_all_nodes_idx,
    
    n2s_faulty_all_nodes_idx =
        n2s_faulty_all_nodes_idx,
    
    fault_nodes_idx =
        fault_nodes_idx,
    
    n2s_fault_nodes_idx =
        n2s_fault_nodes_idx,
    
    list_Ya_nkf =
        list_Ya_nkf,
    
    list_Ynkf_b =
        list_Ynkf_b,
    
    list_faulty_line_Yl =
        list_faulty_line_Yl,
    
    list_healthy_lines_Yl =
        list_healthy_lines_Yl,
    
    list_node_b_idx_in_a_node_row =
        list_node_b_idx_in_a_node_row,
    
    list_node_a_idx_in_b_node_row =
        list_node_a_idx_in_b_node_row,

    list_faulted_line_a_b_orientation =
        list_faulted_line_a_b_orientation,
    
    list_fault_point_from_node_a =
        list_fault_point_from_node_a,
    
    list_fault_resistance =
        list_fault_resistance,
    
    list_no_line_circuit =
        list_no_line_circuit )

    #----------------------------------------
    
    faulty_Ynet =
        deepcopy(faulty_Ynet)
    
    faulty_nodes_idx_with_adjacent_nodes_idx =
        deepcopy(faulty_nodes_idx_with_adjacent_nodes_idx)
    
    pre_clear_fault_Ynet =
        deepcopy(faulty_Ynet)
    
    pre_clear_fault_nodes_idx_with_adjacent_nodes_idx =
        deepcopy(faulty_nodes_idx_with_adjacent_nodes_idx)
    
    pre_clear_fault_all_nodes_idx =
        deepcopy(faulty_all_nodes_idx)
    
    n2s_pre_clear_fault_all_nodes_idx =
        deepcopy(n2s_faulty_all_nodes_idx)
    
    pre_clear_fault_nodes_idx =
        deepcopy(fault_nodes_idx)
    
    n2s_pre_clear_fault_nodes_idx =
        deepcopy(n2s_fault_nodes_idx)
    #
    
    pre_clear_list_Ya_nkf =
        deepcopy(list_Ya_nkf)
    
    pre_clear_list_Ynkf_b =
        deepcopy(list_Ynkf_b)
    
    pre_clear_list_faulty_line_Yl =
        deepcopy(list_faulty_line_Yl)
    
    pre_clear_list_healthy_lines_Yl =
        deepcopy(list_healthy_lines_Yl)
    
    pre_clear_list_node_b_idx_in_a_node_row =
        deepcopy(list_node_b_idx_in_a_node_row)
    
    pre_clear_list_node_a_idx_in_b_node_row =
        deepcopy(list_node_a_idx_in_b_node_row)

    pre_clear_list_faulted_line_a_b_orientation =
        deepcopy(list_faulted_line_a_b_orientation)
    
    pre_clear_list_fault_point_from_node_a =
        deepcopy(list_fault_point_from_node_a)
    
    #
    
    pre_clear_list_fault_resistance =
        deepcopy(list_fault_resistance)

    pre_clear_list_no_line_circuit =
        deepcopy(list_no_line_circuit)
    
    #----------------------------------------
    
    sorted_selection_list = sort(
        clear_fault_selection_list)

    # # This is used to collect only exixsting faults
    # # in the selection list
    
    # exixting_fault_entries = []
    
    #----------------------------------------
    if fault_nodes_idx != []

        for (faulted_line_a_b,
             fault_point_from_node_a,
             fault_resistance,
             no_line_circuit,

             nf,

             Ya_nkf,
             Ynkf_b,

             faulty_line_Yl,
             healthy_lines_Yl,
             node_b_idx_in_a_node_row,
             node_a_idx_in_b_node_row) in
            zip(list_faulted_line_a_b_orientation[
                sorted_selection_list],

                list_fault_point_from_node_a[
                    sorted_selection_list],

                list_fault_resistance[
                    sorted_selection_list],

                list_no_line_circuit[
                    sorted_selection_list],

                fault_nodes_idx[
                    sorted_selection_list],

                list_Ya_nkf[
                    sorted_selection_list],

                list_Ynkf_b[
                    sorted_selection_list],

                list_faulty_line_Yl[
                    sorted_selection_list],

                list_healthy_lines_Yl[
                    sorted_selection_list],

                list_node_b_idx_in_a_node_row[
                    sorted_selection_list],

                list_node_a_idx_in_b_node_row[
                    sorted_selection_list])

            #--------------------------------------

            node_a_idx = first(faulted_line_a_b)

            node_b_idx = last(faulted_line_a_b)

            #--------------------------------------
            # Check for an existing fault entry
            #--------------------------------------

            fault_idx_with_adjacent_nodes_idx =
                faulty_nodes_idx_with_adjacent_nodes_idx[
                    collect(keys(n2s_fault_nodes_idx)) ]

            # Each fault_idx_with_adjacent_nodes_idx has only
            # three entries, fault node, node a and node b

            fault_idx_with_adjacent_nodes_a_idx =
                second.(fault_idx_with_adjacent_nodes_idx)

            fault_idx_with_adjacent_nodes_b_idx =
                third.(fault_idx_with_adjacent_nodes_idx)

            a_idx = findfirst(
                x -> x == node_a_idx,
                fault_idx_with_adjacent_nodes_a_idx )


            b_idx = findfirst(
                x -> x == node_b_idx,
                fault_idx_with_adjacent_nodes_b_idx )

            
            #--------------------------------------

            # a_idx == b_idx means there is an entry

            if a_idx == b_idx

                # push!(exixting_fault_entries, nf)

               if no_line_circuit == 1

            #----------------------------------------

                   deleteat!(
                       faulty_nodes_idx_with_adjacent_nodes_idx[
                        node_a_idx], node_b_idx_in_a_node_row)

                   deleteat!(
                       faulty_nodes_idx_with_adjacent_nodes_idx[
                        node_b_idx], node_a_idx_in_b_node_row)

                   #----------------------------------------

                   # Since the sum of a row should be zero,
                   # remove the contribution of Yl to the first
                   # element and add the contribution of the
                   # new segment

                    faulty_Ynet[ node_a_idx][ 1] =
                            faulty_Ynet[ node_a_idx][ 1] +
                             Ya_nkf

                    deleteat!(faulty_Ynet[ node_a_idx],
                            node_b_idx_in_a_node_row)

                   # Since the sum of a row should be zero,
                   # remove the contribution of Yl to the
                   # first element and add the contribution of
                   # the new segment

                   faulty_Ynet[ node_b_idx][ 1] =
                       faulty_Ynet[ node_b_idx][ 1] + Ynkf_b

                    deleteat!(faulty_Ynet[ node_b_idx],
                            node_a_idx_in_b_node_row )
               else

                   #----------------------------------------

                   # Get the idx of fault node nf in a_node_row
                   # and b_node_row

                   node_nf_idx_in_a_node_row =
                       findfirst(
                           x -> x == nf,
                           faulty_nodes_idx_with_adjacent_nodes_idx[
                               node_a_idx ])

                    node_nf_idx_in_b_node_row  =
                       findfirst(
                           x -> x == nf,
                           faulty_nodes_idx_with_adjacent_nodes_idx[
                               node_b_idx ])

                   # Remove the fault node nf from the
                   # neighbouhood of node_a and node_b

                   deleteat!(
                       faulty_nodes_idx_with_adjacent_nodes_idx[
                        node_a_idx], node_nf_idx_in_a_node_row)

                   deleteat!(
                       faulty_nodes_idx_with_adjacent_nodes_idx[
                        node_b_idx], node_nf_idx_in_b_node_row)

                   #----------------------------------------

                   # Since the sum of a row should be zero,
                   # remove the contribution of the
                   # new segment

                   faulty_Ynet[ node_a_idx][ 1] =
                       faulty_Ynet[ node_a_idx][ 1] + Ya_nkf

                   deleteat!(faulty_Ynet[ node_a_idx],
                             node_nf_idx_in_a_node_row )

                   # Since the sum of a row should be zero,
                   # remove the contribution of the
                   # new segment

                   faulty_Ynet[ node_b_idx][ 1] =
                       faulty_Ynet[ node_b_idx][ 1] + Ynkf_b

                   deleteat!(faulty_Ynet[ node_b_idx],
                             node_nf_idx_in_b_node_row )

               end

            end

        end

    end
    
    if fault_nodes_idx != []


        #----------------------------------------
        # exixting_fault_entries

        deleteat!(faulty_Ynet,
                  fault_nodes_idx[
                      sorted_selection_list])

        deleteat!(
            faulty_nodes_idx_with_adjacent_nodes_idx,
            fault_nodes_idx[
                  sorted_selection_list])

        # if no_line_circuit == 1

        # else
        # end

        #----------------------------------------

        # faulty_all_nodes_idx =
        #     setdiff(faulty_all_nodes_idx,
        #             fault_nodes_idx[
        #                 sorted_selection_list])


        for nf in fault_nodes_idx[sorted_selection_list]

            delete!(n2s_faulty_all_nodes_idx, nf)

            delete!(n2s_fault_nodes_idx, nf)


            a_fault_node_idx =
                findfirst(
                    x -> x == nf, faulty_all_nodes_idx )

            deleteat!(faulty_all_nodes_idx,
                      a_fault_node_idx)

        end

        # remake n2s_faulty_all_nodes_idx and
        # n2s_fault_nodes_idx

        keys_n2s_faulty_all_nodes =
            sort(collect(keys(n2s_faulty_all_nodes_idx)))

        n2s_faulty_all_nodes_idx =
            OrderedDict{Union{Int64,String,Symbol},Int64}(
                a_key => a_value for (a_key, a_value) in
                    zip(keys_n2s_faulty_all_nodes,
                        1:length(
                            keys_n2s_faulty_all_nodes)) )

        keys_n2s_fault_nodes =
            sort(collect(keys(n2s_fault_nodes_idx)))

        n2s_fault_nodes_idx =
            OrderedDict{Union{Int64,String,Symbol},Int64}(
                a_key => a_value for (a_key, a_value) in
                    zip( keys_n2s_fault_nodes,
                        1:length( keys_n2s_fault_nodes)) )

        #----------------------------------------

        deleteat!(list_faulted_line_a_b_orientation,
                sorted_selection_list)

        deleteat!(list_fault_point_from_node_a,
                  sorted_selection_list)

        deleteat!(list_fault_resistance,
                  sorted_selection_list)

        deleteat!(list_no_line_circuit,
                  sorted_selection_list)

        deleteat!(fault_nodes_idx,
                  sorted_selection_list)

        deleteat!(list_Ya_nkf,
                  sorted_selection_list)

        deleteat!(list_Ynkf_b,
                  sorted_selection_list)

        deleteat!(list_faulty_line_Yl,
                  sorted_selection_list)

        deleteat!(list_healthy_lines_Yl,
                  sorted_selection_list)

        deleteat!(list_node_b_idx_in_a_node_row,
                  sorted_selection_list)

        deleteat!(list_node_a_idx_in_b_node_row,
                  sorted_selection_list)
        
    end


    post_clear_fault_Ynet =
        deepcopy(faulty_Ynet)
    
    post_clear_fault_nodes_idx_with_adjacent_nodes_idx =
        deepcopy(faulty_nodes_idx_with_adjacent_nodes_idx)

    post_clear_fault_all_nodes_idx =
        deepcopy(faulty_all_nodes_idx)
    
    n2s_post_clear_fault_all_nodes_idx =
        deepcopy(n2s_faulty_all_nodes_idx)

    post_clear_fault_nodes_idx =
        deepcopy(fault_nodes_idx)
    
    n2s_post_clear_fault_nodes_idx =
        deepcopy(n2s_fault_nodes_idx)

    
    return (;
            
            pre_clear_fault_Ynet,
            pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
            
            pre_clear_fault_all_nodes_idx,
            n2s_pre_clear_fault_all_nodes_idx,
            
            pre_clear_fault_nodes_idx,
            n2s_pre_clear_fault_nodes_idx,
    
            pre_clear_list_Ya_nkf, 
            pre_clear_list_Ynkf_b, 

            pre_clear_list_faulty_line_Yl, 
            pre_clear_list_healthy_lines_Yl, 

            pre_clear_list_node_b_idx_in_a_node_row, 
            pre_clear_list_node_a_idx_in_b_node_row, 

            pre_clear_list_faulted_line_a_b_orientation, 
            pre_clear_list_fault_point_from_node_a,

            pre_clear_list_fault_resistance,
            pre_clear_list_no_line_circuit,
            
            post_clear_fault_Ynet,    
            post_clear_fault_nodes_idx_with_adjacent_nodes_idx,
            
            post_clear_fault_all_nodes_idx,    
            n2s_post_clear_fault_all_nodes_idx,
            
            post_clear_fault_nodes_idx,    
            n2s_post_clear_fault_nodes_idx,
            
            faulty_Ynet,
            faulty_nodes_idx_with_adjacent_nodes_idx,
            
            faulty_all_nodes_idx,
            n2s_faulty_all_nodes_idx,

            fault_nodes_idx,
            n2s_fault_nodes_idx,

            list_faulted_line_a_b_orientation,
            list_fault_point_from_node_a,
            list_fault_resistance,
            list_no_line_circuit,

            list_Ya_nkf,
            list_Ynkf_b,

            list_faulty_line_Yl,
            list_healthy_lines_Yl,

            list_node_b_idx_in_a_node_row,
            list_node_a_idx_in_b_node_row)
    
end




function get_system_simulation_parameters_wt_faults(
    net_data_by_components_file;
    components_libs_dir = nothing,
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,    
    use_nlsolve = false,
    
    pf_alg = NewtonRaphson(),
    
    abstol = 1e-12,    
    reltol = 1e-12,
    
    on_fault_time = 9.0,
    clear_fault_time = 9.02,
       
    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [4],

    list_edges_to_have_fault = [ 2 ], 
    clear_fault_selection_list = [ 1 ],

    with_faults = false)
    
    # """

    # #--------------------------------------    

    # case_name = "case9"

    # json_net_data_file =
    #     "net-static-data-avr-sauer-gov-sauer.json"

    # json_net_data_by_components_file =
    #     json_net_data_file
    
    # components_libs_dir =
    #     joinpath(@__DIR__,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )
    
    # json_case_dir =
    #     joinpath( data_dir, case_name, "json")

    # if json_net_data_by_components_file == nothing

    #     net_data_by_components_file =
    #         joinpath(json_case_dir,
    #                  "net_data_by_components_file.json")
    # else

    #     net_data_by_components_file =
    #         joinpath(json_case_dir,
    #                  json_net_data_by_components_file)
        
    # end

        
    # basekV = 1.0
    
    # use_pu_in_PQ = true
    
    # line_data_in_pu = true

    # use_init_u0 = false
    
    # use_nlsolve = false
    
    # pf_alg        = NewtonRaphson()
    
    # abstol        = 1e-12
    # reltol        = 1e-12
    
    # """

    #--------------------------------------
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src")
        
        components_libs_dir =
            joinpath(
                src_dir,
                "components-lib")

    end


    (;plant_generators_data_from_json,
     plant_loads_data_from_json,
     plant_transmission_data_from_json,
     edge_data_from_json,
     shunt_data_from_json,
     baseMVA_data_from_json) =
        NamedTupleTools.select(
            get_net_data_by_components_from_json_file(
                net_data_by_components_file;
                in_components_type_sym =
                    false ),
            (:plant_generators_data_from_json,
             :plant_loads_data_from_json,
             :plant_transmission_data_from_json,
             :edge_data_from_json,
             :shunt_data_from_json,
             :baseMVA_data_from_json))

    #----------------------------------------

    baseMVA = baseMVA_data_from_json
    
    #------------------------------------------------
    
    (;
     P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist ) =
         get_pf_PQ_param_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json;
             baseMVA =
                 baseMVA,
             use_pu_in_PQ =
                 use_pu_in_PQ )

    
    #----------------------------------------
    #----------------------------------------

    ode_gens_para_selections  =
        (:H, :D,
         :X_d, :X_q,                  
         :X_d_dash, :X_q_dash,
         :T_d_dash, :T_q_dash )

    ode_gens_para_sequence_order =
        (:components_data,
         :gen)
        
    ode_gens_generic_selections =
        (:H, :D,
         :ra, :xℓ,
         :X_d, :X_q,
         :X_d_dash,  :X_q_dash,
         :X_d_2dash, :X_q_2dash,
         :T_d_dash,  :T_q_dash)

    ode_gens_generic_sequence_order =
        (:components_data, :gen)
        
    govs_and_avrs_sequence_order =
        ( :components_data,)
    
    govs_and_avrs_selections =
        ( :gov, :avr )
    
    #----------------------------------------
    
    ode_gens_generic_para =
         get_ode_gens_generic_para(
             plant_generators_data_from_json;
             sequence_order =
                 ode_gens_generic_sequence_order,
             selections =
                 ode_gens_generic_selections)
    
    #------------------------------------------------
    #------------------------------------------------

    (;generic_gens_para,
     generic_govs_para,
     generic_avrs_para) =
         get_generic_gens_avr_gov_para(
             plant_generators_data_from_json;
             gens_sequence_order =
                 ode_gens_generic_sequence_order,
            gens_selections =
                ode_gens_generic_selections,
            govs_and_avrs_sequence_order =
                govs_and_avrs_sequence_order,
            govs_and_avrs_selections =
                govs_and_avrs_selections)
    
    #------------------------------------------------
    #------------------------------------------------
        
    (branches_fbus, branches_tbus,
     branches_r, branches_x, branches_b,
     branches_ratio, branches_angle, branches_type) =
          get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    (branches_fbus,
     branches_tbus) =
         get_edges_fbus_tbus_by_json(
             edge_data_from_json)
    
    (edges_r, edges_x, edges_b,
     edges_ratio, edges_angle, edges_type) =
         get_edges_generic_data_by_json(
             edge_data_from_json )
    
    (Gs, Bs) =
        get_nodes_shunts_Gs_and_Bs_by_json(
            shunt_data_from_json)
    
    #------------------------------------------------
    #------------------------------------------------
    
    # edges_orientation =
    #     get_edges_orientation_by_generic(
    #         branches_fbus, branches_tbus )


    # edges_Ybr_cal =
    #     get_edges_Ybr_by_generic(
    #         edges_r, edges_x, edges_b,
    #         edges_ratio, edges_angle,
    #         edges_type, Gs, Bs;
    #         baseMVA = baseMVA, basekV =  basekV )

    # Ybr_cal_and_edges_orientation =
    #     (;edges_Ybr_cal,
    #      edges_orientation)
    
    # Ynet_wt_nodes_idx_wt_adjacent_nodes =
    #     get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic( get_edges_ftbus_and_generic_data_by_json( edge_data_from_json )..., get_nodes_shunts_idx_and_Gs_and_Bs_by_json( shunt_data_from_json)...; baseMVA = baseMVA, basekV = basekV, baseShunt = baseMVA )

    
    #------------------------------------------------
    #------------------------------------------------

    (;edges_orientation,
     edges_Ybr_cal,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json,
                 edge_data_from_json,
                 shunt_data_from_json;
                 baseMVA =
                     baseMVA,
                 basekV =
                     basekV,
                 use_pu_in_PQ =
                     use_pu_in_PQ,
                 line_data_in_pu =
                     line_data_in_pu ),
             (:edges_orientation,
              :edges_Ybr_cal,
              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))
    
    #------------------------------------------------
    #------------------------------------------------
    
    ode_gens_para =
        NamedTupleTools.select(
            ode_gens_generic_para,
            (:H,
             :D,
             :X_d,
             :X_q,
             :X_d_dash,
             :X_q_dash,
             :T_d_dash,
             :T_q_dash))
    
    #------------------------------------------------
    #------------------------------------------------
        
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )
    
    dyn_pf_fun_kwd_net_idxs =
        NamedTupleTools.select(
            net_nodes_type_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))
    
    dyn_pf_fun_kwd_n2s_idxs =
        NamedTupleTools.select(
            get_dict_net_streamlined_idx_by_nodes_type_idxs(
                net_nodes_type_idxs ),
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))

    nodes_with_demands_idx =
        getproperty(net_nodes_type_idxs,
                    :nodes_with_demands_idx )
    
    #------------------------------------------------
    #------------------------------------------------
    
    dyn_pf_mismatch_vars_kwd_para =
        (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para )

    sta_pf_PQ_para =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA =
                baseMVA,
            use_pu_in_PQ =
                use_pu_in_PQ)

    
    gens_vh_slack_θh_para =
        get_gens_vh_slack_θh_para_by_json(
            plant_generators_data_from_json )

    
    sta_pf_vars_and_paras_idx =
        get_sta_pf_vars_and_paras_idx_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    
    pf_sta_ΔPQ_mismatch_parameters =
        get_pf_sta_ΔPQ_mismatch_parameters_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json;
            baseMVA = baseMVA,
            basekV = basekV,
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu)

    #------------------------------------------------
    #------------------------------------------------
    
    # on_fault_time = 9.0
    # clear_fault_time = 9.02
    
    # list_fault_point_from_node_a = [0.3]
    # list_fault_resistance = [0.001]
    # list_no_line_circuit =  [4]

    # list_edges_to_have_fault = [ 2 ]
    # clear_fault_selection_list = [1]
    
    # no_fault            = 0
    # on_fault            = 1
    # clear_fault         = 2
    # partial_clear_fault = 3

    #------------------------------------------------

    list_faulted_line_a_b_orientation =
        edges_orientation[
            list_edges_to_have_fault  ] 
    
    #--------------------------------------

    no_lines_fault =
        length(list_faulted_line_a_b_orientation)

    #--------------------------------------

    no_current_lines_fault =
        no_lines_fault - length(
            clear_fault_selection_list)
    
    #----------------------------------------    
    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))


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
    #----------------------------------------
    
    (;
     state_vars_idx,
     vec_comp_states_Idx,

     plants_states_syms ,
     generic_state_sym, 

     state_labels,
     algebraic_vars_labels,
     network_vars_labels,
     
     generic_model_states_comp_idxs_in_Idx ,
     generic_model_vars_wt_i_dq_Idx_in_state,

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

     dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     dyn_vh_id_iq_V_ref_Tm_Idx,
     dyn_V_ref_Tm_id_iq_vh_Idx,

     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
     dyn_vh_id_iq_ωref0_vref0_porder0_Idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,

     dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_3_gens_type_paras_Idx,
     dyn_5_gens_type_paras_Idx,

     dyn_2_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
     dyn_3_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
                 
     system_states_idx_kwd_para,
     system_paras_idx_kwd_para,
     plants_dyn_fun_idx_kwd_para,
     plants_algeb_fun_idx_kwd_para ) =
         NamedTupleTools.select(
             get_states_Idx_syms_wt_functions(
                 net_data_by_components_file,
                 dyn_pf_fun_kwd_net_idxs,
                 dyn_pf_fun_kwd_n2s_idxs;
                 components_libs_dir =
                     components_libs_dir),
             (:state_vars_idx,
              :vec_comp_states_Idx,

              :plants_states_syms,
              :generic_state_sym, 

              :state_labels,
              :algebraic_vars_labels,
              :network_vars_labels,

              :generic_model_states_comp_idxs_in_Idx ,
              :generic_model_vars_wt_i_dq_Idx_in_state ,

              :comps_callback_paras_funs,     
              :comps_init_funs,
              :comps_output_funs,
              :comps_dyn_funs,
              :ode_comps_dyn_funs,
              :dae_comps_dyn_funs,

              :algebraic_state_sym,
              :model_syms,

              :nodes_names,
              :gens_nodes_names,
              :non_gens_nodes_names,

              :SM_gens_nodes_names,
              :SC_gens_nodes_names,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,

              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :model_mass_matrix,     
              :model_bool_dae_vars,

              :ode_gens_mass_matrix,
              :ode_gens_bool_dae_vars,

              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Png_Qng_Pll_Qll_Idx,

              :dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :dyn_vh_id_iq_V_ref_Tm_Idx,
              :dyn_V_ref_Tm_id_iq_vh_Idx,

              :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,
              :dyn_vh_id_iq_ωref0_vref0_porder0_Idx,

              :ωref0_vref0_porder0_id_iq_vh_Idx,

              :dyn_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_3_gens_type_paras_Idx,
              :dyn_5_gens_type_paras_Idx,

              :dyn_2_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,
              :dyn_3_gens_type_paras_wt_Png_Qng_Pll_Qll_Idx,

              :system_states_idx_kwd_para,
              :system_paras_idx_kwd_para,
              :plants_dyn_fun_idx_kwd_para,
              :plants_algeb_fun_idx_kwd_para))
    
        
    #----------------------------------------    
    # Static power flow
    #----------------------------------------

    (pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
         NamedTupleTools.select(
             pf_sta_ΔPQ_mismatch_parameters,
             (:pf_kw_para,
              :red_types_Idxs_etc,
              :pf_PQ_param) )

    (red_vh_Idxs,
     red_θh_Idxs) =
         NamedTupleTools.select(
             red_types_Idxs_etc,
             (:red_vh_Idxs,
              :red_θh_Idxs) )
    
    #----------------------------------------
    
    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs))]

    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0) 
          
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------
    
    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (;Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para ) 
    
    generic_results_pf_sta_red_sol =
        get_generic_results_pf_sta_red_sol_u(
            pf_sol;
            generic_red_sol_kwd_para =
                generic_red_sol_kwd_para,
            baseMVA =
                baseMVA,
            basekV =
                1.0 )

    #----------------------------------------    
    #----------------------------------------    
    
    (pf_P_gens,
     pf_Q_gens,
     vh,
     θh,
     gens_vh,
     gens_θh) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens,
             :pf_Q_gens,
             :vh,
             :θh,
             :gens_vh,
             :gens_θh ) )
        
    #----------------------------------------
    #----------------------------------------
    # Dynamic simulation
    #----------------------------------------
    #----------------------------------------

    #----------------------------------------
    # Init
    #----------------------------------------
    
    plants_init_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para ,
         comps_init_funs )
    
    plants_generic_states_init_wt_ref =
        plants_generic_model_init_func(
            gens_vh,
            gens_θh,
            pf_P_gens,
            pf_Q_gens,
            ωs;
            kwd_para =
                plants_init_kwd_para )

    (plants_states_init,
     plants_refs ) =
         NamedTupleTools.select(
             plants_generic_states_init_wt_ref,
             (:plants_states_init,
              :plants_refs))


    ( nt_vec_per_paras,
      vec_vec_per_paras ) =
          get_nt_vec_wt_vec_vec_per_paras(
        plants_refs ;
        nt_syms =
            (:ωs,
             :ω_ref,
             :v_ref,
             :p_order,
             :i_d,
             :i_q,
             :gen_vh ) )

    (ω_ref,
     v_ref,
     p_order,
     gens_i_d,
     gens_i_q ) =
         NamedTupleTools.select(
             nt_vec_per_paras, (
                 :ω_ref,
                 :v_ref,
                 :p_order,
                 :i_d,
                 :i_q))

    #----------------------------------------
    # System model init
    #----------------------------------------

    u0_model_states_init =
        Float64[plants_states_init;
                vh;
                θh;
                gens_i_d;
                gens_i_q]
    
    #----------------------------------------
    # callbacks
    #----------------------------------------

    plants_cb_paras_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para,
         comps_callback_paras_funs )

    #----------------------------------------

     (plants_cb_paras_switches,
      list_selected_plants_state_event_cb_paras,
      list_selected_plants_state_affect_cb_paras,

      avrs_govs_cb_sw,
      avrs_govs_cb_sw_Idx ) =
         NamedTupleTools.select(
             plants_generic_model_callback_paras_func(
                 state_vars_idx,
                 plants_states_syms;
                 kwd_para =
                     plants_cb_paras_kwd_para ) ,
             (:plants_cb_paras_switches,
              :list_selected_plants_state_event_cb_paras,
              :list_selected_plants_state_affect_cb_paras,
              
              :plants_avr_gov_cb_para_sw_in_plant,
              :plants_avr_gov_cb_para_sw_idx_in_plant))
    
    cb_states = cb_fun_make_state_callbacks(
        list_selected_plants_state_event_cb_paras,
        list_selected_plants_state_affect_cb_paras )

    #----------------------------------------
    # System model para
    #----------------------------------------

    (P_non_gens,
     Q_non_gens, 
     P_g_loc_load,
     Q_g_loc_load) =
        NamedTupleTools.select(
            sta_pf_PQ_para,
            (:P_non_gens,
             :Q_non_gens,
             :P_g_loc_load,
             :Q_g_loc_load ) )

    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))

    #----------------------------------------

    gens_δ = u0_model_states_init[
        δ_idx_in_state]

    gens_ed_dash = u0_model_states_init[
        ed_dash_idx_in_state]

    gens_eq_dash = u0_model_states_init[
        eq_dash_idx_in_state]

    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    

    post_sta_PQ =
        (;pf_P_gens,
         pf_Q_gens,
         P_non_gens,
         Q_non_gens,
         P_g_loc_load,
         Q_g_loc_load)
    
    #----------------------------------------
    # model para with callback
    #----------------------------------------

    """

    generic_model_dynamics_para =
        Float64[ω_ref;
                v_ref;
                p_order;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]

    model_dynamics_para =
        (;generic_model_dynamics_para,
          plants_cb_paras_switches)
    """
    
    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh]
    
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll = 
        Float64[ω_ref;
                v_ref;
                p_order;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll
    
    model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll
    
    # model_dynamics_para =
    #     (;generic_model_dynamics_para,
    #      plants_cb_paras_switches )          

    #----------------------------------------    
    # algebraic equation kwd_para
    #----------------------------------------    

     pf_generic_gens_para =
        NamedTupleTools.select(
            get_selected_vec_nt_to_vec_vec(
                generic_gens_para,
                nothing;
                selections =
                    (:ra,
                     :X_d,
                     :X_q,     
                     :X_d_dash,
                     :X_q_dash ),
                vec_datatype = Float64 ),
            (:ra,
              :X_d,
              :X_q,     
              :X_d_dash,
             :X_q_dash ) )

    flat_vh_flat_θh_id_iq_u0 =
        [vh;
         θh;
         gens_i_d;
         gens_i_q ]
    
    dyn_pf_solver =
        (;use_init_u0,
         use_nlsolve,
         pf_alg,
         abstol,
         reltol )

         
     # dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     # dyn_pf_vh_vhf_Idx,
     # dyn_pf_θh_θhf_Idx,
    
    #--------------------------------------
    #--------------------------------------

    # list_faulted_line_a_b_orientation =
    #     edges_orientation[[2,5]] 
    
    # list_fault_point_from_node_a = [0.3, 0.4]
    # list_fault_resistance = [0.001, 0.001]
    # list_no_line_circuit =  [4, 4]

    # clear_fault_selection_list = [1, 2]
    
    (Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         Ynet_wt_nodes_idx_wt_adjacent_nodes
    
    edges_orientation =
        getproperty(
            Ybr_cal_and_edges_orientation,
            :edges_orientation)

    #--------------------------------------
    
    on_fault_net_para =
         make_lines_faults_data_set(
             Ynet,
             nodes_idx_with_adjacent_nodes_idx,

             all_nodes_idx,
             n2s_all_nodes_idx,

             list_faulted_line_a_b_orientation ,
             list_fault_point_from_node_a,
             list_fault_resistance,
             list_no_line_circuit )

    (;faulty_Ynet,
     faulty_nodes_idx_with_adjacent_nodes_idx,

     faulty_all_nodes_idx,
     n2s_faulty_all_nodes_idx,

     fault_nodes_idx,
     n2s_fault_nodes_idx,            

     list_Ya_nkf,
     list_Ynkf_b,

     list_faulty_line_Yl,
     list_healthy_lines_Yl,

     list_node_b_idx_in_a_node_row,
     list_node_a_idx_in_b_node_row) =
         NamedTupleTools.select(
             on_fault_net_para,
             (:faulty_Ynet,
              :faulty_nodes_idx_with_adjacent_nodes_idx,

              :faulty_all_nodes_idx,
              :n2s_faulty_all_nodes_idx,

              :fault_nodes_idx,
              :n2s_fault_nodes_idx,            

              :list_Ya_nkf,
              :list_Ynkf_b,

              :list_faulty_line_Yl,
              :list_healthy_lines_Yl,

              :list_node_b_idx_in_a_node_row,
              :list_node_a_idx_in_b_node_row)) 

    #--------------------------------------

    cleared_selected_lines_faults_net_para =
        get_cleared_selected_lines_faults_data_set(
            clear_fault_selection_list;
            deepcopy(on_fault_net_para)...)

    #--------------------------------------

    # clear_fault_net_para =
    #      remove_lines_faults_data_set(
    #          clear_fault_selection_list;
    #          on_fault_net_para...)

    #--------------------------------------
    
    dyn_pf_vh_θh_id_iq_vhf_θhf_Idx =
        get_generic_vh_θh_id_iq_vhf_θhf_Idx(
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault =
                no_lines_fault)


    dyn_pf_vh_vhf_θh_θhf_id_iq_Idx =
         get_generic_vh_vhf_θh_θhf_id_iq_Idx(
            gens_nodes_idx,
            all_nodes_idx;
             no_lines_fault =
                 no_lines_fault)


    state_algebraic_vars_wt_fault_Idx_in_state =
        get_state_algebraic_vars_wt_fault_Idx_in_state(
            state_labels,
            gens_nodes_idx,
            all_nodes_idx;
            no_lines_fault =
                no_lines_fault)

    #--------------------------------------

    no_cleared_lines_fault =
        length(clear_fault_selection_list)

        # no_lines_fault - no_current_lines_fault
    
    #--------------------------------------

    # https://docs.sciml.ai/DiffEqDocs/latest/solvers/dae_solve/
  
    # initializealg : CheckInit, BrownFullBasicInit, 
    #                 ShampineCollocationInit, NoInit

    #--------------------------------------
        
    cb_on_fault = DiscreteCallback(
        (u, t, integrator) ->
            on_fault_condition(
                u, t, integrator, on_fault_time),
        
        (integrator) ->
            on_fault_affect!(
                integrator, no_lines_fault ); 
        save_positions=(true,true),
        initializealg = ShampineCollocationInit() )


    cb_clear_fault = DiscreteCallback(
        (u, t, integrator) ->
            clear_fault_condition(
                u, t, integrator, clear_fault_time),
        
        (integrator) ->
            clear_fault_affect!(
                integrator, no_cleared_lines_fault);
        save_positions=(true,true),
        initializealg = ShampineCollocationInit())
    
    #--------------------------------------

    cb_faults = CallbackSet(
        cb_on_fault, cb_clear_fault)
        
    #--------------------------------------
    # no resize!
    #--------------------------------------


    cb_on_fault_no_resize = DiscreteCallback(
        (u, t, integrator) ->
            on_fault_condition(
                u, t, integrator, on_fault_time),
        
        (integrator) ->
            on_fault_affect!(
                integrator); 
        save_positions=(true,true),
        initializealg = ShampineCollocationInit() )


    cb_clear_fault_no_resize = DiscreteCallback(
        (u, t, integrator) ->
            clear_fault_condition(
                u, t, integrator, clear_fault_time),
        
        (integrator) ->
            clear_fault_affect!(
                integrator);
        save_positions=(true,true),
        initializealg = ShampineCollocationInit() )
    
    #--------------------------------------

    cb_faults_no_resize = CallbackSet(
        cb_on_fault_no_resize, cb_clear_fault_no_resize)

    #--------------------------------------

    algebraic_generic_model_kwd_para =
        (; loc_load_exist,

         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,

         Ynet_wt_nodes_idx_wt_adjacent_nodes
         )


    algebraic_generic_model_sol_kwd_para =
        (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         flat_vh_flat_θh_id_iq_u0,
         dyn_pf_solver,
         algebraic_generic_model_kwd_para)


    algebraic_generic_model_wt_fault_kwd_para =
        (;loc_load_exist,
                  
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
         dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
         
         on_fault_net_para,

         clear_fault_selection_list,

         no_lines_fault,
         no_cleared_lines_fault,

         list_fault_resistance,

         with_faults,

         cleared_selected_lines_faults_net_para )
    
    algebraic_generic_model_wt_fault_sol_kwd_para =
        (;dyn_pf_flat_vh_flat_θh_id_iq_Idx,
         flat_vh_flat_θh_id_iq_u0,
         dyn_pf_solver,
         algebraic_generic_model_wt_fault_kwd_para,
         cleared_selected_lines_faults_net_para
         )

    #----------------------------------------
    #----------------------------------------
    # ODE system dyamanics simulation
    #----------------------------------------    
    #----------------------------------------
    
    #----------------------------------------
    # generic model_dynamics_kwd_para
    #----------------------------------------

    # ωref0_vref0_porder0_id_iq_vh_Idx =
    #     dyn_ωref0_vref0_porder0_id_iq_vh_Idx

    ode_plants_kwd_para =
        (;state_vars_idx,

         ωref0_vref0_porder0_id_iq_vh_Idx,
         
         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         vec_comp_states_Idx,
         
         avrs_govs_cb_sw,
         avrs_govs_cb_sw_Idx,
         
         ode_comps_dyn_funs,
         comps_output_funs,
         ωs)


    dae_plants_kwd_para =
        (;state_vars_idx,

         ωref0_vref0_porder0_id_iq_vh_Idx,
         
         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         vec_comp_states_Idx,
         
         avrs_govs_cb_sw,
         avrs_govs_cb_sw_Idx,
         
         dae_comps_dyn_funs,
         comps_output_funs,
         ωs)

    
    generic_system_dynamics_wt_fault_kwd_para =
        (;
         gens_nodes_idx,
         ωs,
         loc_load_exist,
         state_vars_idx,

         # id_iq_pg_vh_Idx,

         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
         dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

         gens_state_vars_idx_in_state,
         state_vars_and_i_dq_Idx_in_state,
         
         state_vars_and_i_dq_wt_fault_Idx_in_state,
         
         state_algebraic_vars_Idx_in_state,
         state_algebraic_vars_wt_fault_Idx_in_state,

         dyn_pf_flat_vh_flat_θh_id_iq_Idx,

         dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
         dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

         dyn_pf_vh_vhf_Idx,
         dyn_pf_θh_θhf_Idx,
         
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         pf_generic_gens_para,
         
         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         ode_plants_kwd_para,
         dae_plants_kwd_para,

         algebraic_generic_model_wt_fault_kwd_para,
         algebraic_generic_model_wt_fault_sol_kwd_para,

         algebraic_generic_model_kwd_para,
         algebraic_generic_model_sol_kwd_para,

         no_lines_fault,
         no_current_lines_fault,

         cleared_selected_lines_faults_net_para,

         with_faults,
         generic_results_pf_sta_red_sol )

    return (;
            baseMVA,
            basekV,
            gens_nodes_idx,
            u0_model_states_init,
            model_syms,
            model_mass_matrix,
            model_bool_dae_vars,
            
            algebraic_state_sym,
            ode_gens_mass_matrix,
            ode_gens_bool_dae_vars,

            state_labels,
            algebraic_vars_labels,

            nodes_names,
            gens_nodes_names,
            non_gens_nodes_names,

            SM_gens_nodes_names,
            SC_gens_nodes_names,

            plants_states_init,

            plants_cb_paras_switches,
            avrs_govs_cb_sw,
            avrs_govs_cb_sw_Idx,
            cb_states,

            ωref0_vref0_porder0_id_iq_vh,
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

            fault_nodes_idx,
            
            no_cleared_lines_fault,
            cb_faults,
            cb_faults_no_resize,

            cleared_selected_lines_faults_net_para,

            algebraic_generic_model_wt_fault_sol_kwd_para,
            ode_plants_kwd_para,
            generic_system_dynamics_wt_fault_kwd_para,

            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,
            dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,

            loc_load_exist,

            Ybr_cal_and_edges_orientation,
            sta_pf_PQ_para,
            
            ode_gens_generic_para,
            
            pf_kw_para,

            generic_gens_para,
            generic_avrs_para,
            generic_govs_para,
            comps_init_funs,

            gens_state_vars_idx_in_state,

            generic_results_pf_sta_red_sol,

            post_sta_PQ,            
            vh,
            θh,
            gens_i_d,
            gens_i_q)
    
end


function get_system_pre_or_post_fault_parameters(
    system_status,
    Ynet_wt_nodes_idx_wt_adjacent_nodes,
    cleared_selected_lines_faults_net_para;
    pf_sta_ΔPQ_mismatch_parameters,    
    Ybr_cal_and_edges_orientation,
    sta_pf_PQ_para,        
    ode_gens_generic_para,
    generic_gens_para ,
    generic_avrs_para,
    generic_govs_para ,
    comps_init_funs,
    gens_state_vars_idx_in_state )
    
    (pf_kw_para,
     red_types_Idxs_etc,
     pf_PQ_param) =
         NamedTupleTools.select(
             pf_sta_ΔPQ_mismatch_parameters,
             (:pf_kw_para,
              :red_types_Idxs_etc,
              :pf_PQ_param) )
    

    (red_vh_Idxs,
     red_θh_Idxs) =
         NamedTupleTools.select(
             red_types_Idxs_etc,
             (:red_vh_Idxs,
              :red_θh_Idxs) )
    
    #----------------------------------------
    
    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs))]

    
    (;loc_load_exist,
     pf_kw_gens_vh_slack_θh_para,
     # pf_kw_net_para,
     pf_kw_var_idxs,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs ) =
         NamedTupleTools.select(
             pf_kw_para,
             (:loc_load_exist,
              :pf_kw_gens_vh_slack_θh_para,
              # :pf_kw_net_para,
              :pf_kw_var_idxs,
              :pf_kw_PQ_para_idxs,
              :pf_kw_nodes_types_idxs,
              :pf_kw_n2s_idxs))
    
    if system_status == :pre_fault_state
                    
        pf_kw_net_para =
        (;Ynet,
        nodes_idx_with_adjacent_nodes_idx )
        
    elseif system_status == :post_fault_state

        ( Ynet,
          nodes_idx_with_adjacent_nodes_idx,

          all_nodes_idx,
          n2s_all_nodes_idx,

          fault_nodes_idx,
          n2s_fault_nodes_idx) =
              NamedTupleTools.select(
                  cleared_selected_lines_faults_net_para,
                  (:post_clear_fault_Ynet,
          :post_clear_fault_nodes_idx_with_adjacent_nodes_idx,

          :post_clear_fault_all_nodes_idx,
          :n2s_post_clear_fault_all_nodes_idx,

          :post_clear_fault_nodes_idx,
          :post_pre_clear_fault_nodes_idx))
                            
        Ynet_wt_nodes_idx_wt_adjacent_nodes =
            (;Ynet,
             nodes_idx_with_adjacent_nodes_idx )
        
        pf_kw_net_para =
                (;Ynet,
                 nodes_idx_with_adjacent_nodes_idx )
          
    else
        
        throw("A steady state system status can only be " *
            " pre_fault_state or post_fault_state ")
        
    end

    #--------------------------------------
    
    pf_kw_para =
        (;loc_load_exist,
         pf_kw_gens_vh_slack_θh_para,
         pf_kw_net_para,
         pf_kw_var_idxs,
         pf_kw_PQ_para_idxs,
         pf_kw_nodes_types_idxs,
         pf_kw_n2s_idxs)
    
    #----------------------------------------
    
     kwd_sta_sta_ΔPQ_sol_by_json =
        (;
         pf_alg,
         pf_kw_para,
         red_vh_Idxs,
         red_θh_Idxs,
         sta_red_vh_θh_0 ) 
          
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------
    
    pf_sol =
        get_pf_sta_ΔPQ_mismatch_sol_by_generic(
            pf_PQ_param;
            kwd_para =
                kwd_sta_sta_ΔPQ_sol_by_json )
    
    #----------------------------------------
    # Results    
    #----------------------------------------
    
    generic_red_sol_kwd_para =
        (;Ybr_cal_and_edges_orientation,
          sta_pf_PQ_para,
          ode_gens_generic_para,
          pf_kw_para) 

    generic_results_pf_sta_red_sol =
        get_generic_results_pf_sta_red_sol_u(
            pf_sol,
            Ynet_wt_nodes_idx_wt_adjacent_nodes;
            generic_red_sol_kwd_para =
                generic_red_sol_kwd_para,
            baseMVA = 1.0,
            basekV = 1.0,
            wt_branch_current = false )

    #----------------------------------------    
    #----------------------------------------    
    
    (pf_P_gens,
     pf_Q_gens,
     vh,
     θh,
     gens_vh,
     gens_θh) =
        NamedTupleTools.select(
            generic_results_pf_sta_red_sol,
            (:pf_P_gens,
             :pf_Q_gens,
             :vh,
             :θh,
             :gens_vh,
             :gens_θh) )
            
    #----------------------------------------
    #----------------------------------------
    # Dynamic simulation
    #----------------------------------------
    #----------------------------------------    
    
    plants_init_kwd_para =
        (generic_gens_para ,
         generic_avrs_para,
         generic_govs_para ,
         comps_init_funs )
    
    plants_generic_states_init_wt_ref =
        plants_generic_model_init_func(
            gens_vh,
            gens_θh,
            pf_P_gens,
            pf_Q_gens,
            ωs;
            kwd_para =
                plants_init_kwd_para )

    (plants_states_init,
     plants_refs ) =
         NamedTupleTools.select(
             plants_generic_states_init_wt_ref,
             (:plants_states_init,
              :plants_refs))

    ( nt_vec_per_paras,
      vec_vec_per_paras ) =
          get_nt_vec_wt_vec_vec_per_paras(
        plants_refs ;
        nt_syms =
            (:ωs,
             :ω_ref,
             :v_ref,
             :p_order,
             :i_d,
             :i_q,
             :gen_vh ) )

    (ω_ref,
     v_ref,
     p_order,
     gens_i_d,
     gens_i_q ) =
         NamedTupleTools.select(
             nt_vec_per_paras, (
                 :ω_ref,
                 :v_ref,
                 :p_order,
                 :i_d,
                 :i_q))

    #----------------------------------------
    # System model init
    #----------------------------------------

    u0_model_states_init =
        Float64[plants_states_init;
                vh;
                θh;
                gens_i_d;
                gens_i_q]

    #----------------------------------------
    # System model para
    #----------------------------------------

    (P_non_gens,
     Q_non_gens, 
     P_g_loc_load,
     Q_g_loc_load) =
        NamedTupleTools.select(
            sta_pf_PQ_para,
            (:P_non_gens,
             :Q_non_gens,
             :P_g_loc_load,
             :Q_g_loc_load ) )

    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))

    #----------------------------------------

    gens_δ = u0_model_states_init[
        δ_idx_in_state]

    gens_ed_dash = u0_model_states_init[
        ed_dash_idx_in_state]

    gens_eq_dash = u0_model_states_init[
        eq_dash_idx_in_state]

    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    post_sta_PQ =
        (;pf_P_gens,
         pf_Q_gens,
         P_non_gens,
         Q_non_gens,
         P_g_loc_load,
         Q_g_loc_load)
    
    #----------------------------------------
    #----------------------------------------
    
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh]
    
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll = 
        Float64[ω_ref;
                v_ref;
                p_order;
                P_non_gens;
                Q_non_gens;
                P_g_loc_load;
                Q_g_loc_load]

    #----------------------------------------    
    # algebraic 
    #----------------------------------------    

    flat_vh_flat_θh_id_iq_u0 =
        [vh;
         θh;
         gens_i_d;
         gens_i_q ]
    
    return (;
            vh,
            θh,
            gens_vh,
            gens_θh,
            pf_P_gens,
            pf_Q_gens,
            
            ω_ref,
            v_ref,
            p_order,
            gens_i_d,
            gens_i_q,

            gens_δ,
            gens_ed_dash,
            gens_eq_dash,
            
            u0_model_states_init,
            post_sta_PQ,

            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,
            ωref0_vref0_porder0_id_iq_vh,
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
            flat_vh_flat_θh_id_iq_u0,
            

            generic_results_pf_sta_red_sol,
            plants_generic_states_init_wt_ref,
            plants_refs
            ) 

    
end


#--------------------------------------------------------
#--------------------------------------------------------

function algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
    dx,
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para,
    system_fault_status =
        system_fault_status )

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,

     list_fault_resistance,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,

              :list_fault_resistance,

              :with_faults))

    #----------------------------------------

    fault_resistance = list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2     
    
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
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
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    if system_fault_status[1] == 0

        (Ynet,
         nodes_idx_with_adjacent_nodes_idx) =
             Ynet_wt_nodes_idx_wt_adjacent_nodes
    
        #----------------------------------------    

        (dyn_pf_vh_Idxs,
         dyn_pf_θh_Idxs,       
         dyn_pf_id_Idxs,
         dyn_pf_iq_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_flat_vh_flat_θh_id_iq_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,       
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs))

        #----------------------------------------    
        
        fault_nodes_idx = []

        vh  = x[dyn_pf_vh_Idxs]

        θh  = x[dyn_pf_θh_Idxs]

        gens_i_d  = x[dyn_pf_id_Idxs]

        gens_i_q  = x[dyn_pf_iq_Idxs]

        #----------------------------------------

        gens_vh = vh[gens_nodes_idx]

        gens_θh = θh[gens_nodes_idx]
        
    elseif system_fault_status[1] == 1
        
        (Ynet,
         nodes_idx_with_adjacent_nodes_idx,

         all_nodes_idx,
         n2s_all_nodes_idx,

         fault_nodes_idx,
         n2s_fault_nodes_idx) =
            NamedTupleTools.select(
                on_fault_net_para,
                (:faulty_Ynet,
                 :faulty_nodes_idx_with_adjacent_nodes_idx,

                 :faulty_all_nodes_idx,
                 :n2s_faulty_all_nodes_idx,

                 :fault_nodes_idx,
                 :n2s_fault_nodes_idx))

        #--------------------------------------

        # (dyn_pf_vh_Idxs,
        #  dyn_pf_θh_Idxs,       
        #  dyn_pf_id_Idxs,
        #  dyn_pf_iq_Idxs,
        #  dyn_pf_vhf_Idxs,
        #  dyn_pf_θhf_Idxs) =
        #      NamedTupleTools.select(
        #          dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
        #          (:dyn_pf_vh_Idxs,
        #           :dyn_pf_θh_Idxs,       
        #           :dyn_pf_id_Idxs,
        #           :dyn_pf_iq_Idxs,
        #           :dyn_pf_vhf_Idxs,
        #           :dyn_pf_θhf_Idxs))
        
        (dyn_pf_vh_Idxs,
         dyn_pf_vhf_Idxs,
         dyn_pf_θh_Idxs,
         dyn_pf_θhf_Idxs,
         dyn_pf_id_Idxs,
         dyn_pf_iq_Idxs ) =
             NamedTupleTools.select(
                 dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θh_Idxs,
                  :dyn_pf_θhf_Idxs,
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs ))
        
        
        vhn  = x[dyn_pf_vh_Idxs]

        vhf = x[dyn_pf_vhf_Idxs]

        θhn  = x[dyn_pf_θh_Idxs]
        
        θhf = x[dyn_pf_θhf_Idxs]

        gens_i_d  = x[dyn_pf_id_Idxs]

        gens_i_q  = x[dyn_pf_iq_Idxs]

        vh = [vhn; vhf]

        θh = [θhn; θhf]
        
        #----------------------------------------

        gens_vh = vhn[gens_nodes_idx]

        gens_θh = θhn[gens_nodes_idx]

        
    elseif system_fault_status[1] == 2

        (Ynet,
         nodes_idx_with_adjacent_nodes_idx,

         all_nodes_idx,
         n2s_all_nodes_idx,

         fault_nodes_idx,
         n2s_fault_nodes_idx) =
            NamedTupleTools.select(
                remove_lines_faults_data_set(
                    clear_fault_selection_list;
                    deepcopy(on_fault_net_para)...),
                (:faulty_Ynet,
                 :faulty_nodes_idx_with_adjacent_nodes_idx,

                 :faulty_all_nodes_idx,
                 :n2s_faulty_all_nodes_idx,

                 :fault_nodes_idx,
                 :n2s_fault_nodes_idx))
            
        #----------------------------------------    

        (dyn_pf_vh_Idxs,
         dyn_pf_θh_Idxs,       
         dyn_pf_id_Idxs,
         dyn_pf_iq_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_flat_vh_flat_θh_id_iq_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,       
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs))

        #----------------------------------------    

        # fault_nodes_idx = []

        vh  = x[dyn_pf_vh_Idxs]

        θh  = x[dyn_pf_θh_Idxs]

        gens_i_d = x[dyn_pf_id_Idxs]

        gens_i_q = x[dyn_pf_iq_Idxs]

        #----------------------------------------

        gens_vh = vh[gens_nodes_idx]

        gens_θh = θh[gens_nodes_idx]
        
    else

        throw("system_fault_status " * "
               $(system_fault_status[1]) unknown")
        nothing

    end
    
        
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    # -----------------------------------
    # active power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]
    
    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
                               nth_idx ∈ gens_nodes_idx ?
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ])) : 
        ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 + vh[ n2s_all_nodes_idx[nth_idx]] * sum(
            [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) 
                        for nth_idx in all_nodes_idx ] 


    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------

    
    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                             nth_idx ∈ gens_nodes_idx ?
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
             (fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 + vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx ]

    
    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    
    # ------------------------------------
        
    dx .=
        vcat(P_mismatch,
             Q_mismatch,
             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag)
    
    # ------------------------------------

    return nothing

end



function alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
    dx,
    x, δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para,
    system_fault_status =
        system_fault_status)

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,

     list_fault_resistance) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,

              :list_fault_resistance ))

    #----------------------------------------

    fault_resistance = list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2 
    
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
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
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    if system_fault_status[1] == 0

        (Ynet,
         nodes_idx_with_adjacent_nodes_idx) =
             Ynet_wt_nodes_idx_wt_adjacent_nodes
    
        #----------------------------------------    

        (dyn_pf_vh_Idxs,
         dyn_pf_θh_Idxs,       
         dyn_pf_id_Idxs,
         dyn_pf_iq_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_flat_vh_flat_θh_id_iq_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,       
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs))

        #----------------------------------------    
        
        fault_nodes_idx = []

        vh  = x[dyn_pf_vh_Idxs]

        θh  = x[dyn_pf_θh_Idxs]

        gens_i_d  = x[dyn_pf_id_Idxs]

        gens_i_q  = x[dyn_pf_iq_Idxs]

        #----------------------------------------

        gens_vh = vh[gens_nodes_idx]

        gens_θh = θh[gens_nodes_idx]
        
    elseif system_fault_status[1] == 1
        
        (Ynet,
         nodes_idx_with_adjacent_nodes_idx,

         all_nodes_idx,
         n2s_all_nodes_idx,

         fault_nodes_idx,
         n2s_fault_nodes_idx) =
            NamedTupleTools.select(
                on_fault_net_para,
                (:faulty_Ynet,
                 :faulty_nodes_idx_with_adjacent_nodes_idx,

                 :faulty_all_nodes_idx,
                 :n2s_faulty_all_nodes_idx,

                 :fault_nodes_idx,
                 :n2s_fault_nodes_idx))

        # nornmal_all_nodes_idx =
        #     setdiff(all_nodes_idx, fault_nodes_idx)

        
        (
         nornmal_all_nodes_idx, ) =
             NamedTupleTools.select(
                 dyn_pf_fun_kwd_net_idxs,
                 (
                  :all_nodes_idx,))
        
        #--------------------------------------
        
        (dyn_pf_vh_Idxs,
         dyn_pf_θh_Idxs,
         dyn_pf_id_Idxs,
         dyn_pf_iq_Idxs,
         dyn_pf_vhf_Idxs,
         dyn_pf_θhf_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs,
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θhf_Idxs))
        
        
        vhn  = x[ dyn_pf_vh_Idxs]

        vhf  = x[ dyn_pf_vhf_Idxs]

        θhn  = x[ dyn_pf_θh_Idxs]
        
        θhf  = x[ dyn_pf_θhf_Idxs]

        gens_i_d  = x[ dyn_pf_id_Idxs]

        gens_i_q  = x[ dyn_pf_iq_Idxs]

        vh = [vhn; vhf]

        θh = [θhn; θhf]
        
        #----------------------------------------

        gens_vh = vh[gens_nodes_idx]

        gens_θh = θh[gens_nodes_idx]

        
    elseif system_fault_status[1] == 2

        (Ynet,
         nodes_idx_with_adjacent_nodes_idx,

         all_nodes_idx,
         n2s_all_nodes_idx,

         fault_nodes_idx,
         n2s_fault_nodes_idx) =
            NamedTupleTools.select(
                remove_lines_faults_data_set(
                    clear_fault_selection_list;
                    deepcopy(on_fault_net_para)...),
                (:faulty_Ynet,
                 :faulty_nodes_idx_with_adjacent_nodes_idx,

                 :faulty_all_nodes_idx,
                 :n2s_faulty_all_nodes_idx,

                 :fault_nodes_idx,
                 :n2s_fault_nodes_idx))
            
        #----------------------------------------    

        (dyn_pf_vh_Idxs,
         dyn_pf_θh_Idxs,       
         dyn_pf_id_Idxs,
         dyn_pf_iq_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_flat_vh_flat_θh_id_iq_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,       
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs))

        #----------------------------------------    

        # fault_nodes_idx = []

        vh  = x[dyn_pf_vh_Idxs]

        θh  = x[dyn_pf_θh_Idxs]

        gens_i_d = x[dyn_pf_id_Idxs]

        gens_i_q = x[dyn_pf_iq_Idxs]

        #----------------------------------------

        gens_vh = vh[gens_nodes_idx]

        gens_θh = θh[gens_nodes_idx]
        
    else
 
        throw("system_fault_status " * "
               $(system_fault_status[1]) unknown")

    end
    
        
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    # -----------------------------------
    # active power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]
    
    P_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
                        nth_idx ∈ gens_nodes_idx ?
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ])) : 
        (fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 + vh[ n2s_all_nodes_idx[nth_idx]] * sum(
            [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )  
                        for nth_idx in all_nodes_idx ] 


    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------

    
    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                   nth_idx ∈ gens_nodes_idx ?
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
            ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx ]

    
    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    
    # ------------------------------------

    if system_fault_status[1] == 1

        # nornmal_all_nodes_idx = setdiff(all_nodes_idx, fault_nodes_idx)

        P_mismatch_n = P_mismatch[nornmal_all_nodes_idx]
        Q_mismatch_n = Q_mismatch[nornmal_all_nodes_idx]

        P_mismatch_f = P_mismatch[fault_nodes_idx]
        Q_mismatch_f = Q_mismatch[fault_nodes_idx]
        
        dx .=
            vcat(P_mismatch_n,
                 Q_mismatch_n,
                 gens_stator_equations_mismatch_real,
                 gens_stator_equations_mismatch_imag,
                 P_mismatch_f,
                 Q_mismatch_f)
        
    elseif (system_fault_status[1] == 0) || (
        system_fault_status[1] == 2)

        dx .=
            vcat(P_mismatch,
                 Q_mismatch,
                 gens_stator_equations_mismatch_real,
                 gens_stator_equations_mismatch_imag)
    else
        
        throw("system_fault_status " * "
               $(system_fault_status[1]) unknown")
        
    end
        
    # ------------------------------------

    return nothing

end



function alter2_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
    dx,
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para,
    system_fault_status =
        system_fault_status)

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,
     no_lines_fault,
     no_cleared_lines_fault,

     list_fault_resistance,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,
              
              :no_lines_fault,
              :no_cleared_lines_fault,

              :list_fault_resistance,

              :with_faults))

    #----------------------------------------

    fault_resistance = list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2 
    
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
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
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    if system_fault_status[1] == 0

        (Ynet,
         nodes_idx_with_adjacent_nodes_idx) =
             Ynet_wt_nodes_idx_wt_adjacent_nodes

        #----------------------------------------    
        
        ( n2s_all_nodes_idx, ) =
            NamedTupleTools.select(
                on_fault_net_para,
                ( :n2s_faulty_all_nodes_idx, ))
        
        #----------------------------------------    
        if with_faults == true

            (dyn_pf_vh_Idxs,
             dyn_pf_θh_Idxs,
             dyn_pf_id_Idxs,
             dyn_pf_iq_Idxs,
             dyn_pf_vhf_Idxs,
             dyn_pf_θhf_Idxs) =
                 NamedTupleTools.select(
                     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                     (:dyn_pf_vh_Idxs,
                      :dyn_pf_θh_Idxs,
                      :dyn_pf_id_Idxs,
                      :dyn_pf_iq_Idxs,
                      :dyn_pf_vhf_Idxs,
                      :dyn_pf_θhf_Idxs))


            vhn  = x[dyn_pf_vh_Idxs]

            θhn  = x[dyn_pf_θh_Idxs]

            gens_i_d  = x[dyn_pf_id_Idxs]

            gens_i_q  = x[dyn_pf_iq_Idxs]

            vhf = x[dyn_pf_vhf_Idxs]

            θhf = x[dyn_pf_θhf_Idxs]

            vh = [vhn; vhf]

            θh = [θhn; θhf]
            
        else

            (dyn_pf_vh_Idxs,
             dyn_pf_θh_Idxs,
             dyn_pf_id_Idxs,
             dyn_pf_iq_Idxs,
             dyn_pf_vhf_Idxs,
             dyn_pf_θhf_Idxs) =
                 NamedTupleTools.select(
                     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                     (:dyn_pf_vh_Idxs,
                      :dyn_pf_θh_Idxs,
                      :dyn_pf_id_Idxs,
                      :dyn_pf_iq_Idxs,
                      :dyn_pf_vhf_Idxs,
                      :dyn_pf_θhf_Idxs))


            vh  = x[dyn_pf_vh_Idxs]

            θh  = x[dyn_pf_θh_Idxs]

            gens_i_d  = x[dyn_pf_id_Idxs]

            gens_i_q  = x[dyn_pf_iq_Idxs]
            
        end
        
        #----------------------------------------
        
        gens_vh = vh[ gens_nodes_idx]

        gens_θh = θh[ gens_nodes_idx]

        #----------------------------------------

        # fault_nodes_idx = []

        #----------------------------------------

        # fault_P_mismatch = [ 0.0  for a_vhf in vhf]
        
        # fault_Q_mismatch = [ 0.0  for a_θhf in θhf]
        
        #----------------------------------------
                
    elseif system_fault_status[1] == 1
        
        (Ynet,
         nodes_idx_with_adjacent_nodes_idx,

         all_nodes_idx,
         n2s_all_nodes_idx,

         fault_nodes_idx,
         n2s_fault_nodes_idx) =
            NamedTupleTools.select(
                on_fault_net_para,
                (:faulty_Ynet,
                 :faulty_nodes_idx_with_adjacent_nodes_idx,

                 :faulty_all_nodes_idx,
                 :n2s_faulty_all_nodes_idx,

                 :fault_nodes_idx,
                 :n2s_fault_nodes_idx))

        # nornmal_all_nodes_idx =
        #     setdiff(all_nodes_idx, fault_nodes_idx)
        
        (
         nornmal_all_nodes_idx, ) =
             NamedTupleTools.select(
                 dyn_pf_fun_kwd_net_idxs,
                 (
                  :all_nodes_idx,))
        
        #--------------------------------------
        
        (dyn_pf_vh_Idxs,
         dyn_pf_θh_Idxs,
         dyn_pf_id_Idxs,
         dyn_pf_iq_Idxs,
         dyn_pf_vhf_Idxs,
         dyn_pf_θhf_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs,
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θhf_Idxs))
        
        
        vhn  = x[dyn_pf_vh_Idxs]

        θhn  = x[dyn_pf_θh_Idxs]
        
        gens_i_d  = x[dyn_pf_id_Idxs]

        gens_i_q  = x[dyn_pf_iq_Idxs]

        vhf = x[dyn_pf_vhf_Idxs]

        θhf = x[dyn_pf_θhf_Idxs]

        vh = [vhn; vhf]

        θh = [θhn; θhf]
        
        #----------------------------------------

        gens_vh = vh[gens_nodes_idx]

        gens_θh = θh[gens_nodes_idx]

        # -----------------------------------

        # fault_P_mismatch = [ ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[ n2s_all_nodes_idx[nth_idx]] * sum(
        #         [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
        #             cos( θh[ n2s_all_nodes_idx[nth_idx]] -
        #             θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
        #                for (ynj, idx) in
        #                    zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
        #                         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
        #                     for nth_idx in all_nodes_idx if nth_idx ∈ fault_nodes_idx ] 

        # #----------------------------------------    

        # fault_Q_mismatch = [ nth_idx ∈ fault_nodes_idx ? 
        #     ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
        #         sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
        #                for (ynj, idx) in
        #                    zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
        #                         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
        #                           for nth_idx in all_nodes_idx if nth_idx ∈ fault_nodes_idx ? ]

        #----------------------------------------    
        
    elseif system_fault_status[1] == 2

        (Ynet,
         nodes_idx_with_adjacent_nodes_idx,

         all_nodes_idx,
         n2s_all_nodes_idx,

         fault_nodes_idx,
         n2s_fault_nodes_idx) =
            NamedTupleTools.select(
                remove_lines_faults_data_set(
                    clear_fault_selection_list;
                    deepcopy(on_fault_net_para)...),
                (:faulty_Ynet,
                 :faulty_nodes_idx_with_adjacent_nodes_idx,

                 :faulty_all_nodes_idx,
                 :n2s_faulty_all_nodes_idx,

                 :fault_nodes_idx,
                 :n2s_fault_nodes_idx))

        #----------------------------------------    
        
        ( n2s_all_nodes_idx, ) =
            NamedTupleTools.select(
                on_fault_net_para,
                ( :n2s_faulty_all_nodes_idx, ))
        
        #----------------------------------------    
        
        (dyn_pf_vh_Idxs,
         dyn_pf_θh_Idxs,
         dyn_pf_id_Idxs,
         dyn_pf_iq_Idxs,
         dyn_pf_vhf_Idxs,
         dyn_pf_θhf_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs,
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θhf_Idxs))
        
        
        vhn  = x[dyn_pf_vh_Idxs]

        θhn  = x[dyn_pf_θh_Idxs]
        
        gens_i_d  = x[dyn_pf_id_Idxs]

        gens_i_q  = x[dyn_pf_iq_Idxs]

        vhf = x[dyn_pf_vhf_Idxs]
        
        θhf = x[dyn_pf_θhf_Idxs]

        vh = [vhn; vhf]

        θh = [θhn; θhf]        
        
        #----------------------------------------
        
        gens_vh = vh[gens_nodes_idx]

        gens_θh = θh[gens_nodes_idx]

        #----------------------------------------

        # fault_nodes_idx = []

        #----------------------------------------

        # fault_P_mismatch = [ 0.0  for a_vhf in vhf]
        
        # fault_Q_mismatch = [ 0.0  for a_θhf in θhf]
        
        #----------------------------------------    

    else

        throw("system_fault_status " * "
               $(system_fault_status[1]) unknown")
        nothing

    end
    

    # -----------------------------------
    # active power mismatch
    # -----------------------------------

    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]

    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
             nth_idx ∈ gens_nodes_idx ?
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ])) : ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[ n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                        for nth_idx in all_nodes_idx ] 


    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------


    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) : nth_idx ∈ gens_nodes_idx ? 
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) : ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx ]

    
    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

        
    # ------------------------------------

    if system_fault_status[1] == 1

        # P_mismatch_n = P_mismatch[nornmal_all_nodes_idx]
        # Q_mismatch_n = Q_mismatch[nornmal_all_nodes_idx]

        # fault_P_mismatch = P_mismatch[fault_nodes_idx]
        # fault_Q_mismatch = Q_mismatch[fault_nodes_idx]
        
        dx .=
            vcat(P_mismatch[nornmal_all_nodes_idx],
                 Q_mismatch[nornmal_all_nodes_idx],
                 
                 gens_stator_equations_mismatch_real,
                 gens_stator_equations_mismatch_imag,
                 
                 P_mismatch[fault_nodes_idx],
                 Q_mismatch[fault_nodes_idx])
        
    elseif (system_fault_status[1] == 0) || (
        system_fault_status[1] == 2)

        if with_faults == true
        
            # P_mismatch_n = P_mismatch[all_nodes_idx]
            # Q_mismatch_n = Q_mismatch[all_nodes_idx]

            fault_P_mismatch = [ 0.0  for a_vhf in vhf]

            fault_Q_mismatch = [ 0.0  for a_θhf in θhf]

            dx .=
                vcat(P_mismatch[all_nodes_idx],
                     Q_mismatch[all_nodes_idx],
                     
                     gens_stator_equations_mismatch_real,
                     gens_stator_equations_mismatch_imag,
                     
                     fault_P_mismatch,
                     fault_Q_mismatch)

        else

            # P_mismatch_n =
            #     P_mismatch[all_nodes_idx]
            
            # Q_mismatch_n =
            #     Q_mismatch[all_nodes_idx]

            dx .=
                vcat(P_mismatch[all_nodes_idx],
                     Q_mismatch[all_nodes_idx],
                     
                     gens_stator_equations_mismatch_real,
                     gens_stator_equations_mismatch_imag)
            
            
        end
    else
        
        throw("system_fault_status " * "
               $(system_fault_status[1]) unknown")
        
    end
        
    
    # ------------------------------------

    return nothing

end


function alter3_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
    dx,
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para,
    system_fault_status =
        system_fault_status)

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,
     no_lines_fault,
     no_cleared_lines_fault,

     list_fault_resistance,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,
              
              :no_lines_fault,
              :no_cleared_lines_fault,

              :list_fault_resistance,

              :with_faults))

    #----------------------------------------

    fault_resistance = list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2 
    
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
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
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    (dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs,
     dyn_pf_vhf_Idxs,
     dyn_pf_θhf_Idxs) =
         NamedTupleTools.select(
             dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs,
              :dyn_pf_vhf_Idxs,
              :dyn_pf_θhf_Idxs))

    #----------------------------------------     

    gens_i_d  = x[dyn_pf_id_Idxs]

    gens_i_q  = x[dyn_pf_iq_Idxs]

    gens_vh = (x[dyn_pf_vh_Idxs])[gens_nodes_idx]
    
    gens_θh = (x[dyn_pf_θh_Idxs])[gens_nodes_idx]

    #----------------------------------------
    #----------------------------------------
        
    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        (Ynet,
         nodes_idx_with_adjacent_nodes_idx) =
             Ynet_wt_nodes_idx_wt_adjacent_nodes

        #----------------------------------------    
        
        ( n2s_all_nodes_idx, ) =
            NamedTupleTools.select(
                on_fault_net_para,
                ( :n2s_faulty_all_nodes_idx, ))
        
        #----------------------------------------
        
        if with_faults == true

            vhn  = x[dyn_pf_vh_Idxs]

            θhn  = x[dyn_pf_θh_Idxs]

            vhf = x[dyn_pf_vhf_Idxs]

            θhf = x[dyn_pf_θhf_Idxs]

            vh = [vhn; vhf]

            θh = [θhn; θhf]
            
        else

            vh  = x[dyn_pf_vh_Idxs]

            θh  = x[dyn_pf_θh_Idxs]
            
        end
        
        #----------------------------------------
        
        # gens_vh = vh[ gens_nodes_idx]

        # gens_θh = θh[ gens_nodes_idx]

        #----------------------------------------
        # active power mismatch
        #----------------------------------------
        
        """
        gens:

        id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
                vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

        loads:

        Pl_i =
             -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

        """

        # # P_mismatch
        # dx[dyn_pf_vh_Idxs]

        P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
             (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
              vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
              cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                    nth_idx ∈ gens_nodes_with_loc_loads_idx ?
             (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
              sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
              vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
              cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
              P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
              vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                    cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                          for (ynj, idx) in
                              zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
                 nth_idx ∈ gens_nodes_idx ?
            (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
             sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
             vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
             cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
             vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
             cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                    nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ])) :
                                        ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                            n2s_all_nodes_idx[nth_idx]] * sum(
                    [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                        cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                        θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                    nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                            for nth_idx in all_nodes_idx ] 


        # -----------------------------------
        # reactive power mismatch
        # -----------------------------------

        """
        gens:

        id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
                vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

        loads:

        Ql_i =
             -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

        """

        # # Q_mismatch

        # dx[dyn_pf_θh_Idxs]

        Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
            ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
                vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                    nth_idx ∈ gens_nodes_with_loc_loads_idx ?
              (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
               cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
               vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
               sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
               Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
               vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
               sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                               zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                   nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                       nth_idx ∈ gens_nodes_idx ? 
               (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
                cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
                vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
                sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
                vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                   nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
                                       ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                           n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] *
                                               abs(ynj) *
                    sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                    nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                                  for nth_idx in all_nodes_idx ]

        # ------------------------------------
        
        if with_faults == true
        
            # P_mismatch_n = P_mismatch[all_nodes_idx]
            # Q_mismatch_n = Q_mismatch[all_nodes_idx]

            fault_P_mismatch = [ 0.0  for a_vhf in vhf]

            fault_Q_mismatch = [ 0.0  for a_θhf in θhf]

            dx .=
                vcat(P_mismatch[all_nodes_idx],
                     Q_mismatch[all_nodes_idx],
                     
                     gens_stator_equations_mismatch_real,
                     gens_stator_equations_mismatch_imag,
                     
                     fault_P_mismatch,
                     fault_Q_mismatch)

        else

            # P_mismatch_n =
            #     P_mismatch[all_nodes_idx]
            
            # Q_mismatch_n =
            #     Q_mismatch[all_nodes_idx]

            dx .=
                vcat(P_mismatch[all_nodes_idx],
                     Q_mismatch[all_nodes_idx],
                     
                     gens_stator_equations_mismatch_real,
                     gens_stator_equations_mismatch_imag)
            
            
        end

        # ------------------------------------
        
    elseif system_fault_status[1] == 1
        
        (Ynet,
         nodes_idx_with_adjacent_nodes_idx,

         all_nodes_idx,
         n2s_all_nodes_idx,

         fault_nodes_idx,
         n2s_fault_nodes_idx) =
            NamedTupleTools.select(
                on_fault_net_para,
                (:faulty_Ynet,
                 :faulty_nodes_idx_with_adjacent_nodes_idx,

                 :faulty_all_nodes_idx,
                 :n2s_faulty_all_nodes_idx,

                 :fault_nodes_idx,
                 :n2s_fault_nodes_idx))

        # nornmal_all_nodes_idx =
        #     setdiff(all_nodes_idx, fault_nodes_idx)
        
        (
         nornmal_all_nodes_idx, ) =
             NamedTupleTools.select(
                 dyn_pf_fun_kwd_net_idxs,
                 (
                  :all_nodes_idx,))
        
        #--------------------------------------
        
        vhn  = x[dyn_pf_vh_Idxs]

        θhn  = x[dyn_pf_θh_Idxs]

        vhf = x[dyn_pf_vhf_Idxs]

        θhf = x[dyn_pf_θhf_Idxs]

        vh = [vhn; vhf]

        θh = [θhn; θhf]
        
        #----------------------------------------

        # gens_vh = vh[gens_nodes_idx]

        # gens_θh = θh[gens_nodes_idx]
        
        #----------------------------------------


        fault_P_mismatch = [ ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                            for nth_idx in all_nodes_idx if nth_idx ∈ fault_nodes_idx ] 

        #----------------------------------------    

        fault_Q_mismatch = [
            ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                                  for nth_idx in all_nodes_idx if nth_idx ∈ fault_nodes_idx ]

        
        #----------------------------------------    
        # active power mismatch
        #----------------------------------------    
        
        """
        gens:

        id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
                vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

        loads:

        Pl_i =
             -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

        """

        # # P_mismatch
        # dx[dyn_pf_vh_Idxs]

        P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
             (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
              vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
              cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                    nth_idx ∈ gens_nodes_with_loc_loads_idx ?
             (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
              sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
              vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
              cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
              P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
              vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                    cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                          for (ynj, idx) in
                              zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
                 nth_idx ∈ gens_nodes_idx ?
            (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
             sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
             vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
             cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
             vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
             cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                    nodes_idx_with_adjacent_nodes_idx[n2s_all_nodes_idx[ nth_idx] ])])) :
                                        ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                            n2s_all_nodes_idx[nth_idx]] * sum(
                    [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                        cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                        θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                    nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ]) )
                            for nth_idx in all_nodes_idx ] 

        #----------------------------------------    
        # reactive power mismatch
        #---------------------------------------- 

        """
        gens:

        id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
                vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

        loads:

        Ql_i =
             -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

        """

        # # Q_mismatch

        # dx[dyn_pf_θh_Idxs]

        Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
            ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
                vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                    nth_idx ∈ gens_nodes_with_loc_loads_idx ?
              (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
               cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
               vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
               sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
               Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
               vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
               sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                               zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                   nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                       nth_idx ∈ gens_nodes_idx ? 
               (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
                cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
                vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
                sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
                vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                   nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
                                       ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                           n2s_all_nodes_idx[nth_idx]] * sum([vh[n2s_all_nodes_idx[idx]] *
                                               abs(ynj) *
                    sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                    nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                                  for nth_idx in all_nodes_idx ]

        # ------------------------------------
        
        # P_mismatch_n = P_mismatch[nornmal_all_nodes_idx]
        # Q_mismatch_n = Q_mismatch[nornmal_all_nodes_idx]

        # fault_P_mismatch = P_mismatch[fault_nodes_idx]
        # fault_Q_mismatch = Q_mismatch[fault_nodes_idx]
        
        dx .=
            vcat(P_mismatch[nornmal_all_nodes_idx],
                 Q_mismatch[nornmal_all_nodes_idx],
                 
                 gens_stator_equations_mismatch_real,
                 gens_stator_equations_mismatch_imag,
                 
                 # P_mismatch[fault_nodes_idx],
                 # Q_mismatch[fault_nodes_idx]
                 
                 fault_P_mismatch,
                 fault_Q_mismatch )
        
        # ------------------------------------
        
    elseif system_fault_status[1] == 2

        (Ynet,
         nodes_idx_with_adjacent_nodes_idx,

         all_nodes_idx,
         n2s_all_nodes_idx,

         fault_nodes_idx,
         n2s_fault_nodes_idx) =
            NamedTupleTools.select(
                remove_lines_faults_data_set(
                    clear_fault_selection_list;
                    deepcopy(on_fault_net_para)...),
                (:faulty_Ynet,
                 :faulty_nodes_idx_with_adjacent_nodes_idx,

                 :faulty_all_nodes_idx,
                 :n2s_faulty_all_nodes_idx,

                 :fault_nodes_idx,
                 :n2s_fault_nodes_idx))

        #----------------------------------------    
        
        ( n2s_all_nodes_idx, ) =
            NamedTupleTools.select(
                on_fault_net_para,
                ( :n2s_faulty_all_nodes_idx, ))
        
        #----------------------------------------    
        
        vhn  = x[dyn_pf_vh_Idxs]

        θhn  = x[dyn_pf_θh_Idxs]

        vhf = x[dyn_pf_vhf_Idxs]
        
        θhf = x[dyn_pf_θhf_Idxs]

        vh = [vhn; vhf]

        θh = [θhn; θhf]        
        
        #----------------------------------------
        
        # gens_vh = vh[gens_nodes_idx]

        # gens_θh = θh[gens_nodes_idx]
        
        #----------------------------------------   
        # active power mismatch
        #----------------------------------------    
        
        """
        gens:

        id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
                vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

        loads:

        Pl_i =
             -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

        """

        # # P_mismatch
        # dx[dyn_pf_vh_Idxs]

        P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
             (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
              vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
              cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                    nth_idx ∈ gens_nodes_with_loc_loads_idx ?
             (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
              sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
              vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
              cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
              P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
              vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                    cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                          for (ynj, idx) in
                              zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
                 nth_idx ∈ gens_nodes_idx ?
            (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
             sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
             vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
             cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
             vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
             cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                    nodes_idx_with_adjacent_nodes_idx[n2s_all_nodes_idx[ nth_idx]]) ])) :
                                        ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                            n2s_all_nodes_idx[nth_idx]] * sum(
                    [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                        cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                        θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                    nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                            for nth_idx in all_nodes_idx ] 

        # -----------------------------------
        # reactive power mismatch
        # -----------------------------------

        """
        gens:

        id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
                vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

        loads:

        Ql_i =
             -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

        """

        # # Q_mismatch

        # dx[dyn_pf_θh_Idxs]

        Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
            ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
                vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                    nth_idx ∈ gens_nodes_with_loc_loads_idx ?
              (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
               cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
               vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
               sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
               Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
               vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
               sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                               zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                   nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                       nth_idx ∈ gens_nodes_idx ? 
               (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
                cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
                vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
                sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
                vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                                   nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
                                       ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                           n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] *
                                               abs(ynj) *
                    sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                           for (ynj, idx) in
                               zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                    nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                                  for nth_idx in all_nodes_idx ]

        # ------------------------------------

        fault_P_mismatch = [ 0.0  for a_vhf in vhf]

        fault_Q_mismatch = [ 0.0  for a_θhf in θhf]

        dx .=
            vcat(P_mismatch[all_nodes_idx],
                 Q_mismatch[all_nodes_idx],

                 gens_stator_equations_mismatch_real,
                 gens_stator_equations_mismatch_imag,

                 fault_P_mismatch,
                 fault_Q_mismatch)

        # ------------------------------------

    else

        throw("system_fault_status " * "
               $(system_fault_status[1]) unknown")
        nothing

    end

    # # ------------------------------------

    # if system_fault_status[1] == 1

    #     # P_mismatch_n = P_mismatch[nornmal_all_nodes_idx]
    #     # Q_mismatch_n = Q_mismatch[nornmal_all_nodes_idx]

    #     # fault_P_mismatch = P_mismatch[fault_nodes_idx]
    #     # fault_Q_mismatch = Q_mismatch[fault_nodes_idx]
        
    #     dx .=
    #         vcat(P_mismatch[nornmal_all_nodes_idx],
    #              Q_mismatch[nornmal_all_nodes_idx],
                 
    #              gens_stator_equations_mismatch_real,
    #              gens_stator_equations_mismatch_imag,
                 
    #              # P_mismatch[fault_nodes_idx],
    #              # Q_mismatch[fault_nodes_idx]
                 
    #              fault_P_mismatch,
    #              fault_Q_mismatch )
        
    # elseif (system_fault_status[1] == 0) || (
    #     system_fault_status[1] == 2)

    #     if with_faults == true
        
    #         # P_mismatch_n = P_mismatch[all_nodes_idx]
    #         # Q_mismatch_n = Q_mismatch[all_nodes_idx]

    #         fault_P_mismatch = [ 0.0  for a_vhf in vhf]

    #         fault_Q_mismatch = [ 0.0  for a_θhf in θhf]

    #         dx .=
    #             vcat(P_mismatch[all_nodes_idx],
    #                  Q_mismatch[all_nodes_idx],
                     
    #                  gens_stator_equations_mismatch_real,
    #                  gens_stator_equations_mismatch_imag,
                     
    #                  fault_P_mismatch,
    #                  fault_Q_mismatch)

    #     else

    #         # P_mismatch_n =
    #         #     P_mismatch[all_nodes_idx]
            
    #         # Q_mismatch_n =
    #         #     Q_mismatch[all_nodes_idx]

    #         dx .=
    #             vcat(P_mismatch[all_nodes_idx],
    #                  Q_mismatch[all_nodes_idx],
                     
    #                  gens_stator_equations_mismatch_real,
    #                  gens_stator_equations_mismatch_imag)
            
            
    #     end
    # else
        
    #     throw("system_fault_status " * "
    #            $(system_fault_status[1]) unknown")
        
    # end
        
    
    # ------------------------------------

    return nothing

end



function pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
    dx,
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para )

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,
     no_lines_fault,
     no_cleared_lines_fault,

     # list_fault_resistance,

     with_faults,

     cleared_selected_lines_faults_net_para) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,
              
              :no_lines_fault,
              :no_cleared_lines_fault,

              # :list_fault_resistance,

              :with_faults,

              :cleared_selected_lines_faults_net_para))

    #----------------------------------------

    (list_fault_resistance,) =
         NamedTupleTools.select(
             cleared_selected_lines_faults_net_para,
             (:pre_clear_list_fault_resistance,))
    
    #----------------------------------------
    
    fault_resistance =
        list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2 
    
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
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
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    (dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs,
     dyn_pf_vhf_Idxs,
     dyn_pf_θhf_Idxs) =
         NamedTupleTools.select(
             dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs,
              :dyn_pf_vhf_Idxs,
              :dyn_pf_θhf_Idxs))

    #----------------------------------------     

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]

    #----------------------------------------

    gens_i_d  = x[dyn_pf_id_Idxs]

    gens_i_q  = x[dyn_pf_iq_Idxs]

    gens_vh = (x[dyn_pf_vh_Idxs])[gens_nodes_idx]
    
    gens_θh = (x[dyn_pf_θh_Idxs])[gens_nodes_idx]

    #--------------------------------------

    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]
        
    #----------------------------------------
    
    #  gens real part stator equation mismatch
    #----------------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    #----------------------------------------
    #----------------------------------------

    (Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         Ynet_wt_nodes_idx_wt_adjacent_nodes

    #----------------------------------------    

    ( n2s_all_nodes_idx, ) =
        NamedTupleTools.select(
            on_fault_net_para,
            ( :n2s_faulty_all_nodes_idx, ))

    #----------------------------------------

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]

    #----------------------------------------

    # gens_vh = vh[ gens_nodes_idx]

    # gens_θh = θh[ gens_nodes_idx]

    #----------------------------------------
    # active power mismatch
    #----------------------------------------

    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]

    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
             nth_idx ∈ gens_nodes_idx ?
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ])) :
                                    ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                        n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                        for nth_idx in all_nodes_idx ] 


    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------

    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                   nth_idx ∈ gens_nodes_idx ? 
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
                                   ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                       n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] *
                                           abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx ]

    # ------------------------------------

    fault_P_mismatch = [ 1.0 - a_vhf  for a_vhf in vhf]

    fault_Q_mismatch = [ 0.0 - a_θhf  for a_θhf in θhf]

    dx .=
        vcat(P_mismatch[all_nodes_idx],
             Q_mismatch[all_nodes_idx],

             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag,

             fault_P_mismatch,
             fault_Q_mismatch)

    # ------------------------------------
    
    return nothing

end



function Ynet_pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
    dx,
    x,
    (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
     Ynet);
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para )

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,
     no_lines_fault,
     no_cleared_lines_fault,

     # list_fault_resistance,

     with_faults,

     cleared_selected_lines_faults_net_para) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,
              
              :no_lines_fault,
              :no_cleared_lines_fault,

              # :list_fault_resistance,

              :with_faults,

              :cleared_selected_lines_faults_net_para) )

    #----------------------------------------

    (list_fault_resistance,) =
         NamedTupleTools.select(
             cleared_selected_lines_faults_net_para,
             (:pre_clear_list_fault_resistance,))
    
    #----------------------------------------
    
    fault_resistance =
        list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2 
    
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
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
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    (dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs,
     dyn_pf_vhf_Idxs,
     dyn_pf_θhf_Idxs) =
         NamedTupleTools.select(
             dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs,
              :dyn_pf_vhf_Idxs,
              :dyn_pf_θhf_Idxs))

    #----------------------------------------     

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]

    #----------------------------------------

    gens_i_d  = x[dyn_pf_id_Idxs]

    gens_i_q  = x[dyn_pf_iq_Idxs]

    gens_vh = (x[dyn_pf_vh_Idxs])[gens_nodes_idx]
    
    gens_θh = (x[dyn_pf_θh_Idxs])[gens_nodes_idx]

    #--------------------------------------

    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]
        
    #----------------------------------------
    
    #  gens real part stator equation mismatch
    #----------------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    #----------------------------------------
    #----------------------------------------

    # (Ynet,
    #  nodes_idx_with_adjacent_nodes_idx) =
    #      Ynet_wt_nodes_idx_wt_adjacent_nodes

    (nodes_idx_with_adjacent_nodes_idx,) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:nodes_idx_with_adjacent_nodes_idx,))
    
    #----------------------------------------    

    ( n2s_all_nodes_idx, ) =
        NamedTupleTools.select(
            on_fault_net_para,
            ( :n2s_faulty_all_nodes_idx, ))

    #----------------------------------------

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]

    #----------------------------------------

    # gens_vh = vh[ gens_nodes_idx]

    # gens_θh = θh[ gens_nodes_idx]

    #----------------------------------------
    # active power mismatch
    #----------------------------------------

    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]

    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
             nth_idx ∈ gens_nodes_idx ?
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ])) :
                                    ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                        n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                        for nth_idx in all_nodes_idx ] 


    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------

    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                   nth_idx ∈ gens_nodes_idx ? 
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
                                   ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                       n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] *
                                           abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx ]

    # ------------------------------------

    fault_P_mismatch = [ 1.0 - a_vhf  for a_vhf in vhf]

    fault_Q_mismatch = [ 0.0 - a_θhf  for a_θhf in θhf]

    dx .=
        vcat(P_mismatch[all_nodes_idx],
             Q_mismatch[all_nodes_idx],

             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag,

             fault_P_mismatch,
             fault_Q_mismatch)

    # ------------------------------------
    
    return nothing

end



function fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
    dx,
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para)

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,
     no_lines_fault,
     no_cleared_lines_fault,

     # list_fault_resistance,

     with_faults,

     cleared_selected_lines_faults_net_para) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,
              
              :no_lines_fault,
              :no_cleared_lines_fault,

              # :list_fault_resistance,

              :with_faults,

              :cleared_selected_lines_faults_net_para ))

    #----------------------------------------

    (list_fault_resistance,) =
         NamedTupleTools.select(
             cleared_selected_lines_faults_net_para,
             (:pre_clear_list_fault_resistance,))
    
    #----------------------------------------

    fault_resistance = list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2 
    
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs
    # ,
    # n2s_all_nodes_idx
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs
             # ,
             # :n2s_all_nodes_idx
             ))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx
    # ,
    # all_nodes_idx
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx
             # ,
             # :all_nodes_idx
             ))

    #----------------------------------------

     (
      Ynet,
      nodes_idx_with_adjacent_nodes_idx,

      all_nodes_idx,
      n2s_all_nodes_idx,

      fault_nodes_idx,
      n2s_fault_nodes_idx) =
          NamedTupleTools.select(
              cleared_selected_lines_faults_net_para,
              (:pre_clear_fault_Ynet,
      :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,

      :pre_clear_fault_all_nodes_idx,
      :n2s_pre_clear_fault_all_nodes_idx,

      :pre_clear_fault_nodes_idx,
      :n2s_pre_clear_fault_nodes_idx))
    
    # (Ynet,
    #  nodes_idx_with_adjacent_nodes_idx,

    #  all_nodes_idx,
    #  n2s_all_nodes_idx,

    #  fault_nodes_idx,
    #  n2s_fault_nodes_idx) =
    #     NamedTupleTools.select(
    #         on_fault_net_para,
    #         (:faulty_Ynet,
    #          :faulty_nodes_idx_with_adjacent_nodes_idx,

    #          :faulty_all_nodes_idx,
    #          :n2s_faulty_all_nodes_idx,

    #          :fault_nodes_idx,
    #          :n2s_fault_nodes_idx))


    (nornmal_all_nodes_idx, ) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:all_nodes_idx,))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    (dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs,
     dyn_pf_vhf_Idxs,
     dyn_pf_θhf_Idxs) =
         NamedTupleTools.select(
             dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs,
              :dyn_pf_vhf_Idxs,
              :dyn_pf_θhf_Idxs))

    #----------------------------------------

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]

    #----------------------------------------

    gens_i_d  = x[dyn_pf_id_Idxs]

    gens_i_q  = x[dyn_pf_iq_Idxs]

    gens_vh = (x[dyn_pf_vh_Idxs])[gens_nodes_idx]
    
    gens_θh = (x[dyn_pf_θh_Idxs])[gens_nodes_idx]

    #--------------------------------------

    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]
        
    #----------------------------------------
    
    #  gens real part stator equation mismatch
    #----------------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    #----------------------------------------


    fault_P_mismatch = [ ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
        vh[ n2s_all_nodes_idx[nth_idx]] * sum(
            [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                        for nth_idx in all_nodes_idx if nth_idx ∈ fault_nodes_idx ] 

    #----------------------------------------    

    fault_Q_mismatch = [
        ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
        vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx if nth_idx ∈ fault_nodes_idx  ]


    #----------------------------------------    
    # active power mismatch
    #----------------------------------------    

    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]

    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
             nth_idx ∈ gens_nodes_idx ?
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[n2s_all_nodes_idx[ nth_idx] ])])) :
                                    ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                        n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ]) )
                        for nth_idx in all_nodes_idx ] 

    #----------------------------------------    
    # reactive power mismatch
    #---------------------------------------- 

    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                   nth_idx ∈ gens_nodes_idx ? 
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
                                   ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                       n2s_all_nodes_idx[nth_idx]] * sum([vh[n2s_all_nodes_idx[idx]] *
                                           abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx ]

    # ------------------------------------

    # P_mismatch_n = P_mismatch[nornmal_all_nodes_idx]
    # Q_mismatch_n = Q_mismatch[nornmal_all_nodes_idx]

    # fault_P_mismatch = P_mismatch[fault_nodes_idx]
    # fault_Q_mismatch = Q_mismatch[fault_nodes_idx]

    dx .=
        vcat(P_mismatch[nornmal_all_nodes_idx],
             Q_mismatch[nornmal_all_nodes_idx],

             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag,

             # P_mismatch[fault_nodes_idx],
             # Q_mismatch[fault_nodes_idx]

             fault_P_mismatch,
             fault_Q_mismatch )

    # ------------------------------------
    
    return nothing

end


function Ynet_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
    dx,
    x,
    (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
     Ynet);
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para)

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,
     no_lines_fault,
     no_cleared_lines_fault,

     # list_fault_resistance,

     with_faults,

     cleared_selected_lines_faults_net_para) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,
              
              :no_lines_fault,
              :no_cleared_lines_fault,

              # :list_fault_resistance,

              :with_faults,

              :cleared_selected_lines_faults_net_para ))

    #----------------------------------------

    (list_fault_resistance,) =
         NamedTupleTools.select(
             cleared_selected_lines_faults_net_para,
             (:pre_clear_list_fault_resistance,))
    
    #----------------------------------------

    fault_resistance = list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2 
    
    #----------------------------------------    

    (
     dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs
    # ,
    # n2s_all_nodes_idx
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs
             # ,
             # :n2s_all_nodes_idx
             ))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx
    # ,
    # all_nodes_idx
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx
             # ,
             # :all_nodes_idx
             ))

    #----------------------------------------

     (
      # Ynet,
      nodes_idx_with_adjacent_nodes_idx,

      all_nodes_idx,
      n2s_all_nodes_idx,

      fault_nodes_idx,
      n2s_fault_nodes_idx) =
          NamedTupleTools.select(
              cleared_selected_lines_faults_net_para,
              ( # :pre_clear_fault_Ynet,
      :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,

      :pre_clear_fault_all_nodes_idx,
      :n2s_pre_clear_fault_all_nodes_idx,

      :pre_clear_fault_nodes_idx,
      :n2s_pre_clear_fault_nodes_idx))
    
    # (Ynet,
    #  nodes_idx_with_adjacent_nodes_idx,

    #  all_nodes_idx,
    #  n2s_all_nodes_idx,

    #  fault_nodes_idx,
    #  n2s_fault_nodes_idx) =
    #     NamedTupleTools.select(
    #         on_fault_net_para,
    #         (:faulty_Ynet,
    #          :faulty_nodes_idx_with_adjacent_nodes_idx,

    #          :faulty_all_nodes_idx,
    #          :n2s_faulty_all_nodes_idx,

    #          :fault_nodes_idx,
    #          :n2s_fault_nodes_idx))


    (nornmal_all_nodes_idx, ) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:all_nodes_idx,))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    (dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs,
     dyn_pf_vhf_Idxs,
     dyn_pf_θhf_Idxs) =
         NamedTupleTools.select(
             dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs,
              :dyn_pf_vhf_Idxs,
              :dyn_pf_θhf_Idxs))

    #----------------------------------------

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]

    #----------------------------------------

    gens_i_d  = x[dyn_pf_id_Idxs]

    gens_i_q  = x[dyn_pf_iq_Idxs]

    gens_vh = (x[dyn_pf_vh_Idxs])[gens_nodes_idx]
    
    gens_θh = (x[dyn_pf_θh_Idxs])[gens_nodes_idx]

    #--------------------------------------

    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]
        
    #----------------------------------------
    
    #  gens real part stator equation mismatch
    #----------------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    #----------------------------------------


    fault_P_mismatch = [ ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
        vh[ n2s_all_nodes_idx[nth_idx]] * sum(
            [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                        for nth_idx in all_nodes_idx if nth_idx ∈ fault_nodes_idx ] 

    #----------------------------------------    

    fault_Q_mismatch = [
        ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
        vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx if nth_idx ∈ fault_nodes_idx  ]


    #----------------------------------------    
    # active power mismatch
    #----------------------------------------    

    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]

    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
             nth_idx ∈ gens_nodes_idx ?
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[n2s_all_nodes_idx[ nth_idx] ])])) :
                                    ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                        n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ]) )
                        for nth_idx in all_nodes_idx ] 

    #----------------------------------------    
    # reactive power mismatch
    #---------------------------------------- 

    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                   nth_idx ∈ gens_nodes_idx ? 
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
                                   ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                       n2s_all_nodes_idx[nth_idx]] * sum([vh[n2s_all_nodes_idx[idx]] *
                                           abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx ]

    # ------------------------------------

    # P_mismatch_n = P_mismatch[nornmal_all_nodes_idx]
    # Q_mismatch_n = Q_mismatch[nornmal_all_nodes_idx]

    # fault_P_mismatch = P_mismatch[fault_nodes_idx]
    # fault_Q_mismatch = Q_mismatch[fault_nodes_idx]

    dx .=
        vcat(P_mismatch[nornmal_all_nodes_idx],
             Q_mismatch[nornmal_all_nodes_idx],

             gens_stator_equations_mismatch_real,
             gens_stator_equations_mismatch_imag,

             # P_mismatch[fault_nodes_idx],
             # Q_mismatch[fault_nodes_idx]

             fault_P_mismatch,
             fault_Q_mismatch )

    # ------------------------------------
    
    return nothing

end



function post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
    dx,
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para)

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,
     no_lines_fault,
     no_cleared_lines_fault,

     # list_fault_resistance,

     with_faults,

     cleared_selected_lines_faults_net_para) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,
              
              :no_lines_fault,
              :no_cleared_lines_fault,

              # :list_fault_resistance,

              :with_faults,

              :cleared_selected_lines_faults_net_para))

    (list_fault_resistance,) =
         NamedTupleTools.select(
             cleared_selected_lines_faults_net_para,
             (:pre_clear_list_fault_resistance,))

    #----------------------------------------    
    
    fault_resistance = list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2 
    
    #----------------------------------------    

    (dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs
    # ,
    # n2s_all_nodes_idx
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs
             # ,
             # :n2s_all_nodes_idx
             ))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx
    # ,
    # all_nodes_idx
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx
             # ,
             # :all_nodes_idx
             ))

    #----------------------------------------

     (
      Ynet,
      nodes_idx_with_adjacent_nodes_idx,

      all_nodes_idx,
      n2s_all_nodes_idx,

      fault_nodes_idx,
      n2s_fault_nodes_idx) =
          NamedTupleTools.select(
              cleared_selected_lines_faults_net_para,
              (:post_clear_fault_Ynet,    
               :post_clear_fault_nodes_idx_with_adjacent_nodes_idx,

               :post_clear_fault_all_nodes_idx,    
               :n2s_post_clear_fault_all_nodes_idx,

               :post_clear_fault_nodes_idx,    
               :n2s_post_clear_fault_nodes_idx))    


    (nornmal_all_nodes_idx, ) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:all_nodes_idx,))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    (dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs,
     dyn_pf_vhf_Idxs,
     dyn_pf_θhf_Idxs) =
         NamedTupleTools.select(
             dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs,
              :dyn_pf_vhf_Idxs,
              :dyn_pf_θhf_Idxs))

    #----------------------------------------     

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]

    #----------------------------------------

    gens_i_d  = x[dyn_pf_id_Idxs]

    gens_i_q  = x[dyn_pf_iq_Idxs]

    gens_vh = (x[dyn_pf_vh_Idxs])[gens_nodes_idx]
    
    gens_θh = (x[dyn_pf_θh_Idxs])[gens_nodes_idx]

    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]
        
    #----------------------------------------
    
    #  gens real part stator equation mismatch
    #----------------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    #----------------------------------------
    #----------------------------------------

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]        

    #----------------------------------------

    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]

    #----------------------------------------   
    # active power mismatch
    #----------------------------------------    

    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]

    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
             nth_idx ∈ gens_nodes_idx ?
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[n2s_all_nodes_idx[ nth_idx]]) ])) :
                                    ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                        n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                        for nth_idx in all_nodes_idx ] 

    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------

    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                   nth_idx ∈ gens_nodes_idx ? 
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
                                   ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                       n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] *
                                           abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx ]

    # ------------------------------------

    if length(fault_nodes_idx) == 0

        fault_P_mismatch = [ 1.0 - a_vhf  for a_vhf in vhf]

        fault_Q_mismatch = [ 0.0 - a_θhf  for a_θhf in θhf]

        dx .=
            vcat(P_mismatch[all_nodes_idx],
                 Q_mismatch[all_nodes_idx],

                 gens_stator_equations_mismatch_real,
                 gens_stator_equations_mismatch_imag,

                 fault_P_mismatch,
                 fault_Q_mismatch)        
    else


        fault_P_mismatch =
            [ ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[
                                    n2s_all_nodes_idx[nth_idx]]) ] ) )
              for nth_idx in all_nodes_idx
                  if nth_idx ∈ fault_nodes_idx ] 

        #----------------------------------------    

        fault_Q_mismatch =
            [( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                    sin(θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[
                                    n2s_all_nodes_idx[nth_idx]])])) 
             for nth_idx in all_nodes_idx
                 if nth_idx ∈ fault_nodes_idx  ]


        dx .=
            vcat(P_mismatch[nornmal_all_nodes_idx],
                 Q_mismatch[nornmal_all_nodes_idx],

                 gens_stator_equations_mismatch_real,
                 gens_stator_equations_mismatch_imag,

                 fault_P_mismatch,
                 fault_Q_mismatch )
        
    end

    # ------------------------------------
    
    return nothing

end



function Ynet_post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
    dx,
    x,
    (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
     Ynet);
    kwd_para =
        algebraic_generic_model_wt_fault_kwd_para)

    (loc_load_exist,
              
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     
     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_generic_gens_para,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     on_fault_net_para,
     
     clear_fault_selection_list,
     no_lines_fault,
     no_cleared_lines_fault,

     # list_fault_resistance,

     with_faults,

     cleared_selected_lines_faults_net_para) =
         NamedTupleTools.select(
             kwd_para,
             (:loc_load_exist,

              :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
              
              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :on_fault_net_para,
              
              :clear_fault_selection_list,
              
              :no_lines_fault,
              :no_cleared_lines_fault,

              # :list_fault_resistance,

              :with_faults,

              :cleared_selected_lines_faults_net_para))

    (list_fault_resistance,) =
         NamedTupleTools.select(
             cleared_selected_lines_faults_net_para,
             (:pre_clear_list_fault_resistance,))

    #----------------------------------------    
    
    fault_resistance = list_fault_resistance[1]

    fault_factor_P =
        real(fault_resistance)/(abs(fault_resistance))^2
    
    fault_factor_Q =
        imag(fault_resistance)/(abs(fault_resistance))^2 
    
    #----------------------------------------    

    (dyn_δ_Idxs,
     dyn_ed_dash_Idxs,       
     dyn_eq_dash_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (:dyn_δ_Idx,
              :dyn_ed_dash_Idx,       
              :dyn_eq_dash_Idx,
              
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
    
    #----------------------------------------    

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs
    # ,
    # n2s_all_nodes_idx
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs
             # ,
             # :n2s_all_nodes_idx
             ))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx
    # ,
    # all_nodes_idx
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx
             # ,
             # :all_nodes_idx
             ))

    #----------------------------------------

     (
      # Ynet,
      nodes_idx_with_adjacent_nodes_idx,

      all_nodes_idx,
      n2s_all_nodes_idx,

      fault_nodes_idx,
      n2s_fault_nodes_idx) =
          NamedTupleTools.select(
              cleared_selected_lines_faults_net_para,
              ( # :post_clear_fault_Ynet,    
               :post_clear_fault_nodes_idx_with_adjacent_nodes_idx,

               :post_clear_fault_all_nodes_idx,    
               :n2s_post_clear_fault_all_nodes_idx,

               :post_clear_fault_nodes_idx,    
               :n2s_post_clear_fault_nodes_idx))    


    (nornmal_all_nodes_idx, ) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:all_nodes_idx,))
    
    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_generic_gens_para,
                (:ra, :X_d, :X_q, :X_d_dash,
                 :X_q_dash ) )

    #----------------------------------------

    gens_δ =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_δ_Idxs]
    
    gens_ed_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_ed_dash_Idxs]
    
    gens_eq_dash =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_eq_dash_Idxs]
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
            dyn_Q_non_gens_Idxs]

    if loc_load_exist == true

        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_P_gens_loc_load_Idxs]

        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para[
                dyn_Q_gens_loc_load_Idxs]
    end

    #----------------------------------------    
    #----------------------------------------

    (dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs,
     dyn_pf_vhf_Idxs,
     dyn_pf_θhf_Idxs) =
         NamedTupleTools.select(
             dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs,
              :dyn_pf_vhf_Idxs,
              :dyn_pf_θhf_Idxs))

    #----------------------------------------     

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]

    #----------------------------------------

    gens_i_d  = x[dyn_pf_id_Idxs]

    gens_i_q  = x[dyn_pf_iq_Idxs]

    gens_vh = (x[dyn_pf_vh_Idxs])[gens_nodes_idx]
    
    gens_θh = (x[dyn_pf_θh_Idxs])[gens_nodes_idx]

    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]
        
    #----------------------------------------
    
    #  gens real part stator equation mismatch
    #----------------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """

    # dx[dyn_pf_id_Idxs]
    
    gens_stator_equations_mismatch_real = [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin(
                gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """

    # dx[dyn_pf_iq_Idxs]
    
    gens_stator_equations_mismatch_imag = [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[
                n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  

    #----------------------------------------
    #----------------------------------------

    vhn  = x[dyn_pf_vh_Idxs]

    θhn  = x[dyn_pf_θh_Idxs]

    vhf = x[dyn_pf_vhf_Idxs]

    θhf = x[dyn_pf_θhf_Idxs]

    vh = [vhn; vhf]

    θh = [θhn; θhf]        

    #----------------------------------------

    # gens_vh = vh[gens_nodes_idx]

    # gens_θh = θh[gens_nodes_idx]

    #----------------------------------------   
    # active power mismatch
    #----------------------------------------    

    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    # dx[dyn_pf_vh_Idxs]

    P_mismatch = [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
         nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
             nth_idx ∈ gens_nodes_idx ?
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[n2s_all_nodes_idx[ nth_idx]]) ])) :
                                    ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                        n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) )
                        for nth_idx in all_nodes_idx ] 

    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------

    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    # dx[dyn_pf_θh_Idxs]

    Q_mismatch = [  nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                   nth_idx ∈ gens_nodes_idx ? 
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) :
                                   ( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +  vh[
                                       n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] *
                                           abs(ynj) *
                sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) 
                              for nth_idx in all_nodes_idx ]

    # ------------------------------------

    if length(fault_nodes_idx) == 0

        fault_P_mismatch = [ 1.0 - a_vhf  for a_vhf in vhf]

        fault_Q_mismatch = [ 0.0 - a_θhf  for a_θhf in θhf]

        dx .=
            vcat(P_mismatch[all_nodes_idx],
                 Q_mismatch[all_nodes_idx],

                 gens_stator_equations_mismatch_real,
                 gens_stator_equations_mismatch_imag,

                 fault_P_mismatch,
                 fault_Q_mismatch)        
    else


        fault_P_mismatch =
            [ ( fault_factor_P * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
                    cos( θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[
                                    n2s_all_nodes_idx[nth_idx]]) ] ) )
              for nth_idx in all_nodes_idx
                  if nth_idx ∈ fault_nodes_idx ] 

        #----------------------------------------    

        fault_Q_mismatch =
            [( fault_factor_Q * (vh[ n2s_all_nodes_idx[nth_idx]] )^2 +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum(
                [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                    sin(θh[ n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[
                                    n2s_all_nodes_idx[nth_idx]])])) 
             for nth_idx in all_nodes_idx
                 if nth_idx ∈ fault_nodes_idx  ]


        dx .=
            vcat(P_mismatch[nornmal_all_nodes_idx],
                 Q_mismatch[nornmal_all_nodes_idx],

                 gens_stator_equations_mismatch_real,
                 gens_stator_equations_mismatch_imag,

                 fault_P_mismatch,
                 fault_Q_mismatch )
        
    end

    # ------------------------------------
    
    return nothing

end



# ------------------------------------

function algebraic_generic_pf_ΔPQ_pre_fault_mismatch_sol(
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_sol_kwd_para,
    system_fault_status =
        system_fault_status)

    
    (flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_wt_fault_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:flat_vh_flat_θh_id_iq_u0,
              :dyn_pf_solver,
              :algebraic_generic_model_wt_fault_kwd_para ))

    #----------------------------------------    
    
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol)  =
         NamedTupleTools.select(
             dyn_pf_solver,
             (:use_init_u0,
              :use_nlsolve,
              :pf_alg,
              :abstol,
              :reltol))

    #----------------------------------------        


    if use_init_u0 == true
        
        flat_vh_flat_θh_id_iq =
            flat_vh_flat_θh_id_iq_u0
    else

        flat_vh_flat_θh_id_iq = x
        
    end

    pf_fun_mismatch =
        pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!

    #----------------------------------------    
        
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status),
                      flat_vh_flat_θh_id_iq,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    kwd_para =
                        algebraic_generic_model_wt_fault_kwd_para,
                    system_fault_status =
                        system_fault_status)),
            flat_vh_flat_θh_id_iq,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para ),
                                    pf_alg )

        
    end
    

end



function algebraic_generic_pf_ΔPQ_fault_mismatch_sol(
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_sol_kwd_para )

    
    (flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_wt_fault_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:flat_vh_flat_θh_id_iq_u0,
              :dyn_pf_solver,
              :algebraic_generic_model_wt_fault_kwd_para ))

    #----------------------------------------    
    
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol)  =
         NamedTupleTools.select(
             dyn_pf_solver,
             (:use_init_u0,
              :use_nlsolve,
              :pf_alg,
              :abstol,
              :reltol))

    #----------------------------------------        


    if use_init_u0 == true
        
        flat_vh_flat_θh_id_iq =
            flat_vh_flat_θh_id_iq_u0
    else

        flat_vh_flat_θh_id_iq = x
        
    end

    pf_fun_mismatch =
        fault_state_algebraic_generic_pf_ΔPQ_mismatch!

    #----------------------------------------    
        
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status),
                      flat_vh_flat_θh_id_iq,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    kwd_para =
                        algebraic_generic_model_wt_fault_kwd_para,
                    system_fault_status =
                        system_fault_status)),
            flat_vh_flat_θh_id_iq,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para ),
                                    pf_alg )

        
    end
    
end


function algebraic_generic_pf_ΔPQ_post_fault_mismatch_sol(
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_sol_kwd_para,
    system_fault_status =
        system_fault_status)

    
    (flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_wt_fault_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:flat_vh_flat_θh_id_iq_u0,
              :dyn_pf_solver,
              :algebraic_generic_model_wt_fault_kwd_para ))

    #----------------------------------------    
    
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol)  =
         NamedTupleTools.select(
             dyn_pf_solver,
             (:use_init_u0,
              :use_nlsolve,
              :pf_alg,
              :abstol,
              :reltol))

    #----------------------------------------        


    if use_init_u0 == true
        
        flat_vh_flat_θh_id_iq =
            flat_vh_flat_θh_id_iq_u0
    else

        flat_vh_flat_θh_id_iq = x
        
    end

    pf_fun_mismatch =
        post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!

    #----------------------------------------    
        
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status),
                      flat_vh_flat_θh_id_iq,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    kwd_para =
                        algebraic_generic_model_wt_fault_kwd_para,
                    system_fault_status =
                        system_fault_status)),
            flat_vh_flat_θh_id_iq,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para ),
                                    pf_alg )
        
    end
    
end



# algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!


function algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_sol_kwd_para,
    system_fault_status =
        system_fault_status)

    
    (dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_wt_fault_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:dyn_pf_flat_vh_flat_θh_id_iq_Idx,
              :flat_vh_flat_θh_id_iq_u0,
              :dyn_pf_solver,
              :algebraic_generic_model_wt_fault_kwd_para ))

    #----------------------------------------    
    
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol)  =
         NamedTupleTools.select(
             dyn_pf_solver,
             (:use_init_u0,
              :use_nlsolve,
              :pf_alg,
              :abstol,
              :reltol))

    #----------------------------------------        


    if use_init_u0 == true
        
        flat_vh_flat_θh_id_iq =
            flat_vh_flat_θh_id_iq_u0
    else

        flat_vh_flat_θh_id_iq = x
        
    end

    pf_fun_mismatch =
        algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!
    
        # algebraic_generic_pf_ΔPQ_mismatch!
        # ode_algebraic_generic_model_func!

    #----------------------------------------    
        
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status),
                      flat_vh_flat_θh_id_iq,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    kwd_para =
                        algebraic_generic_model_wt_fault_kwd_para,
                    system_fault_status =
                        system_fault_status) ),
            flat_vh_flat_θh_id_iq,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para),
                                    pf_alg )

        
    end
    

end


function alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
    x, δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_sol_kwd_para,
    system_fault_status =
        system_fault_status)

    
    (
     flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_wt_fault_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (
              :flat_vh_flat_θh_id_iq_u0,
              :dyn_pf_solver,
              :algebraic_generic_model_wt_fault_kwd_para ))

    #----------------------------------------    
    
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol)  =
         NamedTupleTools.select(
             dyn_pf_solver,
             (:use_init_u0,
              :use_nlsolve,
              :pf_alg,
              :abstol,
              :reltol))

    #----------------------------------------        


    if use_init_u0 == true
        
        flat_vh_flat_θh_id_iq =
            flat_vh_flat_θh_id_iq_u0
    else

        flat_vh_flat_θh_id_iq = x
        
    end

    pf_fun_mismatch =
        alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!
    
        # algebraic_generic_pf_ΔPQ_mismatch!
        # ode_algebraic_generic_model_func!

    #----------------------------------------    
        
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status),
                      flat_vh_flat_θh_id_iq,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    kwd_para =
                        algebraic_generic_model_wt_fault_kwd_para,
                    system_fault_status =
                        system_fault_status)),
            flat_vh_flat_θh_id_iq,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para),
                                    pf_alg )

        
    end
    

end


function alter2_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
    x,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    kwd_para =
        algebraic_generic_model_wt_fault_sol_kwd_para,
    system_fault_status =
        system_fault_status)

    
    (
     flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_wt_fault_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (
              :flat_vh_flat_θh_id_iq_u0,
              :dyn_pf_solver,
              :algebraic_generic_model_wt_fault_kwd_para ))

    #----------------------------------------    
    
    (;use_init_u0,
     use_nlsolve,
     pf_alg,
     abstol,
     reltol)  =
         NamedTupleTools.select(
             dyn_pf_solver,
             (:use_init_u0,
              :use_nlsolve,
              :pf_alg,
              :abstol,
              :reltol))

    #----------------------------------------        


    if use_init_u0 == true
        
        flat_vh_flat_θh_id_iq =
            flat_vh_flat_θh_id_iq_u0
    else

        flat_vh_flat_θh_id_iq = x
        
    end

    pf_fun_mismatch =
        alter2_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!

    #----------------------------------------    
        
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status),
                      flat_vh_flat_θh_id_iq,
                      method =
                          :trust_region,
                      autodiff =
                          :forward )
        
    else

        sol = NonlinearSolve.solve( NonlinearProblem(
            NonlinearFunction( (g, x, p) ->
                pf_fun_mismatch(
                    g, x, p;
                    kwd_para =
                        algebraic_generic_model_wt_fault_kwd_para,
                    system_fault_status =
                        system_fault_status)),
            flat_vh_flat_θh_id_iq,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para ),
                                    pf_alg )

        
    end
    

end



function generic_dynamics_wt_fault_by_ode_pf_funcs!(
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_wt_fault_sol_kwd_para,
     algebraic_generic_model_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_wt_fault_sol_kwd_para,
              :algebraic_generic_model_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------

    (pf_vh_Idxs,
     pf_θh_Idxs,       
     pf_id_Idxs,
     pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))

    idx_range =
        first(pf_vh_Idxs):last(pf_iq_Idxs)
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))

    # -------------------------------------
    
    gens_δ = x[
        δ_idx_in_state]
        
    gens_ed_dash = x[
        ed_dash_idx_in_state]
    
    gens_eq_dash = x[
        eq_dash_idx_in_state]
    
    # -------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

        (;state_var_Idx_in_state,
         algebraic_var_Idx_in_state ) =
             NamedTupleTools.select(
                 state_algebraic_vars_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :algebraic_var_Idx_in_state))

        #----------------------------------------

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        vh =  x[vh_Idx_in_state]

        θh =  x[θh_Idx_in_state]

        gens_i_d =  x[id_Idx_in_state]

        gens_i_q =  x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]
        
                
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
   
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state,
         vhf_Idx_in_state,
         θhf_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_wt_fault_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state,
                  :vhf_Idx_in_state,
                  :θhf_Idx_in_state))
    
        (;
         algebraic_wt_fault_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_wt_fault_Idx_in_state,
                 (
                  :algebraic_wt_fault_Idx_in_state,))

        #----------------------------------------

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        du_algebraic_vars_view =
            @view dx[algebraic_wt_fault_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_wt_fault_Idx_in_state]

        #----------------------------------------

        vhn = x[vh_Idx_in_state]

        θhn = x[θh_Idx_in_state]

        gens_i_d = x[id_Idx_in_state]

        gens_i_q = x[ iq_Idx_in_state]

        vhf = x[vhf_Idx_in_state]

        θhf = x[θhf_Idx_in_state]

        #----------------------------------------
        
        vh = [vhn; vhf]
        
        θh = [θhn; θhf]

        #----------------------------------------

        gens_vh =
            vhn[ gens_nodes_idx ]

        gens_θh =
            θhn[ gens_nodes_idx ]
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))


        (;
         algebraic_var_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_Idx_in_state,
                 (
                  :algebraic_var_Idx_in_state,))

        #----------------------------------------

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vh = x[vh_Idx_in_state]

        θh = x[θh_Idx_in_state]

        gens_i_d = x[id_Idx_in_state]

        gens_i_q = x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
    
    
    #----------------------------------------
    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    #----------------------------------------
    #----------------------------------------
    # using plant dynamics 
    #----------------------------------------
    #----------------------------------------
    
    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            ode_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    vh_θh_id_iq =
        [vh;
         θh;
         gens_i_d;
         gens_i_q]

    
    pf_sol =
        algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
            vh_θh_id_iq,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                algebraic_generic_model_wt_fault_sol_kwd_para,
            system_fault_status =
                system_fault_status)

    
    # pf_sol =
    #     algebraic_generic_pf_ΔPQ_mismatch_sol(
    #         vh_θh_id_iq,
    #         δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    #         kwd_para =
    #             algebraic_generic_model_sol_kwd_para )
       
    if system_fault_status[1] == 0
        
        du_algebraic_vars_view .=
            pf_sol[ idx_range ] -  u_algebraic_vars_view
        
    elseif system_fault_status[1] == 1

        (pf_vh_Idxs,
         pf_vhf_Idxs,
         pf_θh_Idxs,
         pf_θhf_Idxs,
         pf_id_Idxs,
         pf_iq_Idxs ) =
             NamedTupleTools.select(
                 dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θh_Idxs,
                  :dyn_pf_θhf_Idxs,
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs ))

        post_pf_sol =
            [pf_sol[pf_vh_Idxs];
             pf_sol[pf_θh_Idxs];
             pf_sol[pf_id_Idxs];
             pf_sol[pf_iq_Idxs];
             pf_sol[pf_vhf_Idxs];
             pf_sol[pf_θhf_Idxs]]
        
        du_algebraic_vars_view .=
            post_pf_sol -  u_algebraic_vars_view
        
    elseif system_fault_status[1] == 2
        
        du_algebraic_vars_view .=
            pf_sol[ idx_range ] -  u_algebraic_vars_view
        
    else
        
        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
    end
    
    #----------------------------------------
        
    return nothing

end

 
function alter_generic_dynamics_wt_fault_by_ode_pf_funcs!(
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

        (;state_var_Idx_in_state,
         algebraic_var_Idx_in_state ) =
             NamedTupleTools.select(
                 state_algebraic_vars_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :algebraic_var_Idx_in_state))

        #----------------------------------------

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vh = x[vh_Idx_in_state]

        θh = x[θh_Idx_in_state]

        gens_i_d = x[id_Idx_in_state]

        gens_i_q = x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        ode_gens_plants_generic_model_func!(
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                ode_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        vh_θh_id_iq =
            [vh;
             θh;
             gens_i_d;
             gens_i_q]


        # vh_θh_id_iq =
        #     x[algebraic_var_Idx_in_state]
        

        pf_sol =
            alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
                vh_θh_id_iq,
                δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_sol_kwd_para,
                system_fault_status =
                    system_fault_status)

        #----------------------------------------

        (pf_vh_Idxs,
         pf_θh_Idxs,       
         pf_id_Idxs,
         pf_iq_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_flat_vh_flat_θh_id_iq_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,       
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs ))
        
        idx_range =
            first(pf_vh_Idxs):last(pf_iq_Idxs)
        
        du_algebraic_vars_view .=
            pf_sol[ idx_range ] -  u_algebraic_vars_view
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
        
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state,
         vhf_Idx_in_state,
         θhf_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_wt_fault_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state,
                  :vhf_Idx_in_state,
                  :θhf_Idx_in_state))
    
        (;
         algebraic_wt_fault_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_wt_fault_Idx_in_state,
                 (
                  :algebraic_wt_fault_Idx_in_state,))

        #----------------------------------------

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        du_algebraic_vars_view =
            @view dx[algebraic_wt_fault_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_wt_fault_Idx_in_state]

        #----------------------------------------

        vhn = x[vh_Idx_in_state] # n : normal

        θhn = x[θh_Idx_in_state]

        gens_i_d = x[id_Idx_in_state]

        gens_i_q = x[ iq_Idx_in_state]

        vhf = x[vhf_Idx_in_state] # f : fault nodes

        θhf = x[θhf_Idx_in_state]

        #----------------------------------------
        
        # vh = [vhn; vhf]
        
        θh = [θhn; θhf]

        #----------------------------------------

        gens_vh =
            vhn[ gens_nodes_idx ]

        gens_θh =
            θhn[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        ode_gens_plants_generic_model_func!(
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                ode_plants_kwd_para )


        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        vh_θh_id_iq =
            [vhn;
             θhn;
             gens_i_d;
             gens_i_q;
             vhf;
             θhf]
        
        # vh_θh_id_iq =
        #     x[algebraic_var_Idx_in_state]

        pf_sol =
            alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
                vh_θh_id_iq,
                δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_sol_kwd_para,
                system_fault_status =
                    system_fault_status)

        #----------------------------------------
        
        (pf_vh_Idxs,
         pf_θh_Idxs,
         pf_id_Idxs,
         pf_iq_Idxs,
         pf_vhf_Idxs,
         pf_θhf_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs,
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θhf_Idxs))
        
        idx_range =
            first(pf_vh_Idxs):last(pf_θhf_Idxs)
        
        du_algebraic_vars_view .=
            pf_sol[ idx_range ] -  u_algebraic_vars_view
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))


        (;state_var_Idx_in_state,
         algebraic_var_Idx_in_state ) =
             NamedTupleTools.select(
                 state_algebraic_vars_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :algebraic_var_Idx_in_state))

        #----------------------------------------

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vh = x[vh_Idx_in_state]

        θh = x[θh_Idx_in_state]

        gens_i_d = x[id_Idx_in_state]

        gens_i_q = x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]
        

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        ode_gens_plants_generic_model_func!(
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                ode_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        vh_θh_id_iq =
            [vh;
             θh;
             gens_i_d;
             gens_i_q]

        pf_sol =
            algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
                vh_θh_id_iq,
                δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_sol_kwd_para,
                system_fault_status =
                    system_fault_status)

        #----------------------------------------

        (pf_vh_Idxs,
         pf_θh_Idxs,       
         pf_id_Idxs,
         pf_iq_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_flat_vh_flat_θh_id_iq_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,       
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs ))
        
        idx_range =
            first(pf_vh_Idxs):last(pf_iq_Idxs)
        
        du_algebraic_vars_view .=
            pf_sol[ idx_range ] -  u_algebraic_vars_view
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


function alter2_generic_dynamics_wt_fault_by_ode_pf_funcs!(
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))
         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        if with_faults == true

            (;state_var_Idx_in_state,
             vh_Idx_in_state,
             θh_Idx_in_state,
             id_Idx_in_state,
             iq_Idx_in_state,
             vhf_Idx_in_state,
             θhf_Idx_in_state) = NamedTupleTools.select(
                 state_vars_and_i_dq_wt_fault_Idx_in_state,
                     (:state_var_Idx_in_state,
                      :vh_Idx_in_state,
                      :θh_Idx_in_state,
                      :id_Idx_in_state,
                      :iq_Idx_in_state,
                      :vhf_Idx_in_state,
                      :θhf_Idx_in_state))

            (
             algebraic_var_Idx_in_state, ) =
                 NamedTupleTools.select(
                     state_algebraic_vars_wt_fault_Idx_in_state,
                     (
                      :algebraic_wt_fault_Idx_in_state,))

            (pf_vh_Idxs,
             pf_θh_Idxs,
             pf_id_Idxs,
             pf_iq_Idxs,
             pf_vhf_Idxs,
             pf_θhf_Idxs) =
                 NamedTupleTools.select(
                     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                     (:dyn_pf_vh_Idxs,
                      :dyn_pf_θh_Idxs,
                      :dyn_pf_id_Idxs,
                      :dyn_pf_iq_Idxs,
                      :dyn_pf_vhf_Idxs,
                      :dyn_pf_θhf_Idxs))

            #----------------------------------------
        
            idx_range =
                first(pf_vh_Idxs):last(pf_θhf_Idxs)

            #----------------------------------------

            vhn = x[vh_Idx_in_state]

            θhn = x[θh_Idx_in_state]

            gens_i_d = x[id_Idx_in_state]

            gens_i_q = x[ iq_Idx_in_state]

            vhf = x[ vhf_Idx_in_state] # f : fault nodes

            θhf = x[ θhf_Idx_in_state]
            

            vh_θh_id_iq_vhf_θhf =
                [vhn;
                 θhn;
                 gens_i_d;
                 gens_i_q;
                 vhf;
                 θhf]
            
        else

            (;state_var_Idx_in_state,
             vh_Idx_in_state,
             θh_Idx_in_state,
             id_Idx_in_state,
             iq_Idx_in_state) =
                 NamedTupleTools.select(
                     state_vars_and_i_dq_Idx_in_state,
                     (:state_var_Idx_in_state,
                      :vh_Idx_in_state,
                      :θh_Idx_in_state,
                      :id_Idx_in_state,
                      :iq_Idx_in_state))

            (
             algebraic_var_Idx_in_state, ) =
                 NamedTupleTools.select(
                      state_algebraic_vars_Idx_in_state,
                     (
                      :algebraic_var_Idx_in_state,))

            (pf_vh_Idxs,
             pf_θh_Idxs,
             pf_id_Idxs,
             pf_iq_Idxs) =
                 NamedTupleTools.select(
                     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                     (:dyn_pf_vh_Idxs,
                      :dyn_pf_θh_Idxs,
                      :dyn_pf_id_Idxs,
                      :dyn_pf_iq_Idxs))

            #----------------------------------------
        
            idx_range =
                first(pf_vh_Idxs):last(pf_iq_Idxs)

            #----------------------------------------
            
            vhn = x[vh_Idx_in_state]

            θhn = x[θh_Idx_in_state]

            gens_i_d = x[id_Idx_in_state]

            gens_i_q = x[ iq_Idx_in_state]
            

            vh_θh_id_iq_vhf_θhf =
                [vhn;
                 θhn;
                 gens_i_d;
                 gens_i_q]
            
        end

        #----------------------------------------

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------
        
        #----------------------------------------
        
        gens_vh =
            vhn[ gens_nodes_idx ]

        gens_θh =
            θhn[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        ode_gens_plants_generic_model_func!(
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                ode_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        # vh_θh_id_iq_vhf_θhf =
        #     x[algebraic_var_Idx_in_state]
        

        pf_sol =
            alter2_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
                vh_θh_id_iq_vhf_θhf,
                δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_sol_kwd_para,
                system_fault_status =
                    system_fault_status)

        #----------------------------------------
        
        du_algebraic_vars_view .=
            pf_sol[ idx_range ] -  u_algebraic_vars_view
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
        
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state,
         vhf_Idx_in_state,
         θhf_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_wt_fault_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state,
                  :vhf_Idx_in_state,
                  :θhf_Idx_in_state))
    
        (
         algebraic_wt_fault_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_wt_fault_Idx_in_state,
                 (
                  :algebraic_wt_fault_Idx_in_state,))
        
        (pf_vh_Idxs,
         pf_θh_Idxs,
         pf_id_Idxs,
         pf_iq_Idxs,
         pf_vhf_Idxs,
         pf_θhf_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs,
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θhf_Idxs))

        #----------------------------------------

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        du_algebraic_vars_view =
            @view dx[algebraic_wt_fault_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_wt_fault_Idx_in_state]

        #----------------------------------------

        vhn = x[vh_Idx_in_state] # n : normal

        θhn = x[θh_Idx_in_state]

        gens_i_d = x[id_Idx_in_state]

        gens_i_q = x[ iq_Idx_in_state]

        vhf = x[vhf_Idx_in_state] # f : fault nodes

        θhf = x[θhf_Idx_in_state]

        #----------------------------------------

        gens_vh =
            vhn[ gens_nodes_idx ]

        gens_θh =
            θhn[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        ode_gens_plants_generic_model_func!(
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                ode_plants_kwd_para )


        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]
        
        # vh_θh_id_iq_vhf_θhf =
        #     x[algebraic_var_Idx_in_state]


        vh_θh_id_iq_vhf_θhf =
            [vhn;
             θhn;
             gens_i_d;
             gens_i_q;
             vhf;
             θhf]
        
        
        pf_sol =
            alter2_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
                vh_θh_id_iq_vhf_θhf,
                δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_sol_kwd_para,
                system_fault_status =
                    system_fault_status)

        #----------------------------------------
        
        idx_range =
            first(pf_vh_Idxs):last(pf_θhf_Idxs)
        
        du_algebraic_vars_view .=
            pf_sol[ idx_range ] -  u_algebraic_vars_view
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]
        
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state,
         vhf_Idx_in_state,
         θhf_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_wt_fault_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state,
                  :vhf_Idx_in_state,
                  :θhf_Idx_in_state))

        (
         algebraic_var_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_wt_fault_Idx_in_state,
                 (
                  :algebraic_wt_fault_Idx_in_state,))
        
        (pf_vh_Idxs,
         pf_θh_Idxs,
         pf_id_Idxs,
         pf_iq_Idxs,
         pf_vhf_Idxs,
         pf_θhf_Idxs) =
             NamedTupleTools.select(
                 dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
                 (:dyn_pf_vh_Idxs,
                  :dyn_pf_θh_Idxs,
                  :dyn_pf_id_Idxs,
                  :dyn_pf_iq_Idxs,
                  :dyn_pf_vhf_Idxs,
                  :dyn_pf_θhf_Idxs))

        #----------------------------------------
        
        idx_range =
            first(pf_vh_Idxs):last(pf_θhf_Idxs)

        #----------------------------------------

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vhn = x[vh_Idx_in_state] # n : normal

        θhn = x[θh_Idx_in_state]

        gens_i_d = x[id_Idx_in_state]

        gens_i_q = x[ iq_Idx_in_state]

        vhf = x[vhf_Idx_in_state] # f : fault nodes

        θhf = x[θhf_Idx_in_state]

        #----------------------------------------
                
        gens_vh =
            vhn[ gens_nodes_idx ]

        gens_θh =
            θhn[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]
        
        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        ode_gens_plants_generic_model_func!(
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                ode_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]


        # vh_θh_id_iq_vhf_θhf =
        #     x[algebraic_var_Idx_in_state]


        vh_θh_id_iq_vhf_θhf =
            [vhn;
             θhn;
             gens_i_d;
             gens_i_q;
             vhf;
             θhf]
        
        pf_sol =
            alter2_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch_sol(
                vh_θh_id_iq_vhf_θhf,
                δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
                kwd_para =
                    algebraic_generic_model_wt_fault_sol_kwd_para,
                system_fault_status =
                    system_fault_status)

        #----------------------------------------
        
        du_algebraic_vars_view .=
            pf_sol[ idx_range ] -  u_algebraic_vars_view
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end



function generic_dynamics_wt_fault_by_dae_pf_funcs!(
    res,
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     dae_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,     

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_wt_fault_kwd_para,
              
              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status
        
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

        (;state_var_Idx_in_state,
         algebraic_var_Idx_in_state ) =
             NamedTupleTools.select(
                 state_algebraic_vars_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :algebraic_var_Idx_in_state))

        #----------------------------------------

        res_states_vars_views =
            [ view(res, idx) for idx in state_vars_idx]

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        res_state_vars_view =
            @view res[state_var_Idx_in_state]

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------
    
        res_algebraic_vars_view =
            @view res[algebraic_var_Idx_in_state]

        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vh = @view x[vh_Idx_in_state]

        θh = @view x[θh_Idx_in_state]

        gens_i_d = @view x[id_Idx_in_state]

        gens_i_q = @view x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        dae_gens_plants_generic_model_func!(
            res_state_vars_view,
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                dae_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (; δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             system_fault_status);
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para  )
        
    elseif system_fault_status[1] == 1

        @show system_fault_status
   
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state,
         vhf_Idx_in_state,
         θhf_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_wt_fault_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state,
                  :vhf_Idx_in_state,
                  :θhf_Idx_in_state))
    
        (;state_var_Idx_in_state,
         algebraic_wt_fault_Idx_in_state ) =
             NamedTupleTools.select(
                 state_algebraic_vars_wt_fault_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :algebraic_wt_fault_Idx_in_state))

        #----------------------------------------

        res_states_vars_views =
            [ view(res, idx) for idx in state_vars_idx]

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        res_state_vars_view =
            @view res[state_var_Idx_in_state]

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        res_algebraic_vars_view =
            @view res[algebraic_wt_fault_Idx_in_state]

        du_algebraic_vars_view =
            @view dx[algebraic_wt_fault_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_wt_fault_Idx_in_state]

        #----------------------------------------

        vhn = x[vh_Idx_in_state] # n : normal

        θhn = x[θh_Idx_in_state]

        gens_i_d = x[id_Idx_in_state]

        gens_i_q = x[ iq_Idx_in_state]

        vhf = x[vhf_Idx_in_state] # f : fault nodes

        θhf = x[θhf_Idx_in_state]

        #----------------------------------------
        
        vh = [vhn; vhf]
        
        θh = [θhn; θhf]

        #----------------------------------------

        gens_vh =
            vhn[ gens_nodes_idx ]

        gens_θh =
            θhn[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        dae_gens_plants_generic_model_func!(
            res_state_vars_view,
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                dae_plants_kwd_para )


        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))


        (;state_var_Idx_in_state,
         algebraic_var_Idx_in_state ) =
             NamedTupleTools.select(
                 state_algebraic_vars_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :algebraic_var_Idx_in_state))

        #----------------------------------------

        res_states_vars_views =
            [ view(res, idx) for idx in state_vars_idx]
        
        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        res_state_vars_view =
            @view res[state_var_Idx_in_state]

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        res_algebraic_vars_view =
            @view res[algebraic_var_Idx_in_state]
        
        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vh = x[vh_Idx_in_state]

        θh = x[θh_Idx_in_state]

        gens_i_d = x[id_Idx_in_state]

        gens_i_q = x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]
        

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        dae_gens_plants_generic_model_func!(
            res_state_vars_view,
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                dae_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status)
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


function alter_generic_dynamics_wt_fault_by_dae_pf_funcs!(
    res,
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     dae_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,     

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_wt_fault_kwd_para,
              
              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]
        
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

        (
         algebraic_var_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_Idx_in_state,
                 (
                  :algebraic_var_Idx_in_state,))

        #----------------------------------------

        res_states_vars_views =
            [ view(res, idx) for idx in state_vars_idx]

        du_states_vars_views =
            [ view(dx, idx)  for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx)   for idx in state_vars_idx]

        #----------------------------------------

        res_state_vars_view =
            @view res[state_var_Idx_in_state]

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------
    
        res_algebraic_vars_view =
            @view res[algebraic_var_Idx_in_state]

        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vh = @view x[vh_Idx_in_state]

        θh = @view x[θh_Idx_in_state]

        gens_i_d = @view x[id_Idx_in_state]

        gens_i_q = @view x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        dae_gens_plants_generic_model_func!(
            res_state_vars_view,
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                dae_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status)
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
   
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state,
         vhf_Idx_in_state,
         θhf_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_wt_fault_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state,
                  :vhf_Idx_in_state,
                  :θhf_Idx_in_state))
    
        (;
         algebraic_wt_fault_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_wt_fault_Idx_in_state,
                 (
                  :algebraic_wt_fault_Idx_in_state,))

        #----------------------------------------

        res_states_vars_views =
            [ view(res, idx) for idx in state_vars_idx]

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        res_state_vars_view =
            @view res[state_var_Idx_in_state]

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        res_algebraic_vars_view =
            @view res[algebraic_wt_fault_Idx_in_state]

        du_algebraic_vars_view =
            @view dx[algebraic_wt_fault_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_wt_fault_Idx_in_state]

        #----------------------------------------

        vh =  @view x[vh_Idx_in_state] # n : normal

        θh =  @view x[θh_Idx_in_state]

        gens_i_d =  @view x[id_Idx_in_state]

        gens_i_q =  @view x[ iq_Idx_in_state]

        #----------------------------------------

        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        dae_gens_plants_generic_model_func!(
            res_state_vars_view,
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                dae_plants_kwd_para )


        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para ;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para ,
            system_fault_status =
                system_fault_status)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))


        # (;state_var_Idx_in_state,
        #  algebraic_var_Idx_in_state ) =
        #      NamedTupleTools.select(
        #          state_algebraic_vars_Idx_in_state,
        #          (:state_var_Idx_in_state,
        #           :algebraic_var_Idx_in_state))

        (
         algebraic_var_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_Idx_in_state,
                 (
                  :algebraic_var_Idx_in_state,))

        #----------------------------------------

        res_states_vars_views =
            [ view(res, idx) for idx in state_vars_idx]
        
        du_states_vars_views =
            [ view(dx, idx)  for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx)   for idx in state_vars_idx]

        #----------------------------------------

        res_state_vars_view =
            @view res[state_var_Idx_in_state]

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        res_algebraic_vars_view =
            @view res[algebraic_var_Idx_in_state]
        
        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vh =  @view x[vh_Idx_in_state]

        θh =  @view x[θh_Idx_in_state]

        gens_i_d =  @view x[id_Idx_in_state]

        gens_i_q =  @view x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]
        

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        dae_gens_plants_generic_model_func!(
            res_state_vars_view,
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                dae_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        alter_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para ;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status)
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


function alter2_generic_dynamics_wt_fault_by_dae_pf_funcs!(
    res,
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     dae_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,     

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_wt_fault_kwd_para,
              
              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]
        
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

        if with_faults == true
            
            (
             algebraic_var_Idx_in_state, ) =
                 NamedTupleTools.select(
                     state_algebraic_vars_wt_fault_Idx_in_state,
                     (
                      :algebraic_wt_fault_Idx_in_state,))
            
        else
            
            (
             algebraic_var_Idx_in_state, ) =
                 NamedTupleTools.select(
                     state_algebraic_vars_Idx_in_state,
                     (
                      :algebraic_var_Idx_in_state,))
            
        end
        

        #----------------------------------------

        res_states_vars_views =
            [ view(res, idx) for idx in state_vars_idx]

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        res_state_vars_view =
            @view res[state_var_Idx_in_state]

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------
    
        res_algebraic_vars_view =
            @view res[algebraic_var_Idx_in_state]

        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vh = @view x[vh_Idx_in_state]

        θh = @view x[θh_Idx_in_state]

        gens_i_d = @view x[id_Idx_in_state]

        gens_i_q = @view x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        dae_gens_plants_generic_model_func!(
            res_state_vars_view,
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                dae_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        alter2_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
             kwd_para =
                 algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status)
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
   
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state,
         vhf_Idx_in_state,
         θhf_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_wt_fault_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state,
                  :vhf_Idx_in_state,
                  :θhf_Idx_in_state))
    
        (
         algebraic_wt_fault_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_wt_fault_Idx_in_state,
                 (
                  :algebraic_wt_fault_Idx_in_state,))

        #----------------------------------------

        res_states_vars_views =
            [ view(res, idx) for idx in state_vars_idx]

        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        res_state_vars_view =
            @view res[state_var_Idx_in_state]

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        res_algebraic_vars_view =
            @view res[algebraic_wt_fault_Idx_in_state]

        du_algebraic_vars_view =
            @view dx[algebraic_wt_fault_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_wt_fault_Idx_in_state]

        #----------------------------------------

        vh =  @view x[vh_Idx_in_state] # n : normal

        θh =  @view x[θh_Idx_in_state]

        gens_i_d =  @view x[id_Idx_in_state]

        gens_i_q =  @view x[ iq_Idx_in_state]

        #----------------------------------------

        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        dae_gens_plants_generic_model_func!(
            res_state_vars_view,
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                dae_plants_kwd_para )


        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        alter2_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) =
             NamedTupleTools.select(
                 state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))


        # (;state_var_Idx_in_state,
        #  algebraic_var_Idx_in_state ) =
        #      NamedTupleTools.select(
        #          state_algebraic_vars_Idx_in_state,
        #          (:state_var_Idx_in_state,
        #           :algebraic_var_Idx_in_state))

        (
         algebraic_var_Idx_in_state, ) =
             NamedTupleTools.select(
                 state_algebraic_vars_wt_fault_Idx_in_state,
                 (
                  :algebraic_wt_fault_Idx_in_state,))

        #----------------------------------------

        res_states_vars_views =
            [ view(res, idx) for idx in state_vars_idx]
        
        du_states_vars_views =
            [ view(dx, idx) for idx in state_vars_idx]

        u_states_vars_views =
            [ view(x, idx) for idx in state_vars_idx]

        #----------------------------------------

        res_state_vars_view =
            @view res[state_var_Idx_in_state]

        du_state_vars_view =
            @view dx[state_var_Idx_in_state]

        u_state_vars_view =
            @view x[state_var_Idx_in_state]

        #----------------------------------------

        res_algebraic_vars_view =
            @view res[algebraic_var_Idx_in_state]
        
        du_algebraic_vars_view =
            @view dx[algebraic_var_Idx_in_state]

        u_algebraic_vars_view =
            @view x[algebraic_var_Idx_in_state]

        #----------------------------------------

        vh =  @view x[vh_Idx_in_state]

        θh =  @view x[θh_Idx_in_state]

        gens_i_d =  @view x[id_Idx_in_state]

        gens_i_q =  @view x[ iq_Idx_in_state]

        #----------------------------------------
        
        gens_vh =
            vh[ gens_nodes_idx ]

        gens_θh =
            θh[ gens_nodes_idx ]

        ωref0_vref0_porder0_id_iq_vh =
            [ω_ref;
             v_ref;
             p_order;
             gens_i_d;
             gens_i_q;
             gens_vh ]
        

        #----------------------------------------
        # using plant dynamics 
        #----------------------------------------

        dae_gens_plants_generic_model_func!(
            res_state_vars_view,
            du_state_vars_view,
            u_state_vars_view,
            ωref0_vref0_porder0_id_iq_vh,
            t;
            kwd_para =
                dae_plants_kwd_para )

        #----------------------------------------
        # power flow equations
        #----------------------------------------

        gens_δ = x[
            δ_idx_in_state]

        gens_ed_dash = x[
            ed_dash_idx_in_state]

        gens_eq_dash = x[
            eq_dash_idx_in_state]

        #----------------------------------------

        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
            [gens_δ;
             gens_ed_dash;
             gens_eq_dash;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

        alter2_algebraic_generic_pf_ΔPQ_wt_fault_and_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para ;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para,
            system_fault_status =
                system_fault_status)
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end



# ------------------------------------------------------
# Fault DAE
# ------------------------------------------------------

function generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!(
    res,
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     Ynet,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     dae_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,     

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_wt_fault_kwd_para,
              
              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    (;state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state,
     vhf_Idx_in_state,
     θhf_Idx_in_state) = NamedTupleTools.select(
         state_vars_and_i_dq_wt_fault_Idx_in_state,
             (:state_var_Idx_in_state,
              :vh_Idx_in_state,
              :θh_Idx_in_state,
              :id_Idx_in_state,
              :iq_Idx_in_state,
              :vhf_Idx_in_state,
              :θhf_Idx_in_state))

    (
     algebraic_wt_fault_Idx_in_state, ) =
         NamedTupleTools.select(
             state_algebraic_vars_wt_fault_Idx_in_state,
             (
              :algebraic_wt_fault_Idx_in_state,))

    #----------------------------------------

    res_states_vars_views =
        [ view(res, idx) for idx in state_vars_idx]

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------

    res_state_vars_view =
        @view res[state_var_Idx_in_state]

    du_state_vars_view =
        @view dx[state_var_Idx_in_state]

    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------

    res_algebraic_vars_view =
        @view res[algebraic_wt_fault_Idx_in_state]

    du_algebraic_vars_view =
        @view dx[algebraic_wt_fault_Idx_in_state]

    u_algebraic_vars_view =
        @view x[algebraic_wt_fault_Idx_in_state]

    #----------------------------------------

    vh =  @view x[vh_Idx_in_state] # n : normal

    θh =  @view x[θh_Idx_in_state]

    gens_i_d =  @view x[id_Idx_in_state]

    gens_i_q =  @view x[iq_Idx_in_state]

    vhf = @view x[vhf_Idx_in_state]
    
    θhf = @view x[θhf_Idx_in_state]
    
    #----------------------------------------

    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]

    #----------------------------------------
    # using plant dynamics 
    #----------------------------------------

    dae_gens_plants_generic_model_func!(
        res_state_vars_view,
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            dae_plants_kwd_para )


    #----------------------------------------
    # power flow equations
    #----------------------------------------

    gens_δ = x[
        δ_idx_in_state]

    gens_ed_dash = x[
        ed_dash_idx_in_state]

    gens_eq_dash = x[
        eq_dash_idx_in_state]

    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
             kwd_para =
                 algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
   

        fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para ;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


         
function Ynet_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!(
    res,
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     Ynet,
     fault_Ynet,
     post_fault_Ynet,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     dae_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,     

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_wt_fault_kwd_para,
              
              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    (;state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state,
     vhf_Idx_in_state,
     θhf_Idx_in_state) = NamedTupleTools.select(
         state_vars_and_i_dq_wt_fault_Idx_in_state,
             (:state_var_Idx_in_state,
              :vh_Idx_in_state,
              :θh_Idx_in_state,
              :id_Idx_in_state,
              :iq_Idx_in_state,
              :vhf_Idx_in_state,
              :θhf_Idx_in_state))

    (
     algebraic_wt_fault_Idx_in_state, ) =
         NamedTupleTools.select(
             state_algebraic_vars_wt_fault_Idx_in_state,
             (
              :algebraic_wt_fault_Idx_in_state,))

    #----------------------------------------

    res_states_vars_views =
        [ view(res, idx) for idx in state_vars_idx]

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------

    res_state_vars_view =
        @view res[state_var_Idx_in_state]

    du_state_vars_view =
        @view dx[state_var_Idx_in_state]

    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------

    res_algebraic_vars_view =
        @view res[algebraic_wt_fault_Idx_in_state]

    du_algebraic_vars_view =
        @view dx[algebraic_wt_fault_Idx_in_state]

    u_algebraic_vars_view =
        @view x[algebraic_wt_fault_Idx_in_state]

    #----------------------------------------

    vh =  @view x[vh_Idx_in_state] # n : normal

    θh =  @view x[θh_Idx_in_state]

    gens_i_d =  @view x[id_Idx_in_state]

    gens_i_q =  @view x[iq_Idx_in_state]

    vhf = @view x[vhf_Idx_in_state]
    
    θhf = @view x[θhf_Idx_in_state]
    
    #----------------------------------------

    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]

    #----------------------------------------
    # using plant dynamics 
    #----------------------------------------

    dae_gens_plants_generic_model_func!(
        res_state_vars_view,
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            dae_plants_kwd_para )


    #----------------------------------------
    # power flow equations
    #----------------------------------------

    gens_δ = x[
        δ_idx_in_state]

    gens_ed_dash = x[
        ed_dash_idx_in_state]

    gens_eq_dash = x[
        eq_dash_idx_in_state]

    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        Ynet_pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             Ynet);
             kwd_para =
                 algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
        
        
        Ynet_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             fault_Ynet);
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        Ynet_post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             post_fault_Ynet) ;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


         
function line_loss_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!(
    res,
    dx,
    x,
    (;generic_model_dynamics_para,
     Ynet,
     fault_Ynet,
     post_fault_Ynet,
     system_fault_status,
     plants_cb_paras_switches),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     dae_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,     

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_wt_fault_kwd_para,
              
              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref = generic_model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = generic_model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    (;state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state,
     vhf_Idx_in_state,
     θhf_Idx_in_state) = NamedTupleTools.select(
         state_vars_and_i_dq_wt_fault_Idx_in_state,
             (:state_var_Idx_in_state,
              :vh_Idx_in_state,
              :θh_Idx_in_state,
              :id_Idx_in_state,
              :iq_Idx_in_state,
              :vhf_Idx_in_state,
              :θhf_Idx_in_state))

    (
     algebraic_wt_fault_Idx_in_state, ) =
         NamedTupleTools.select(
             state_algebraic_vars_wt_fault_Idx_in_state,
             (
              :algebraic_wt_fault_Idx_in_state,))

    #----------------------------------------

    res_states_vars_views =
        [ view(res, idx) for idx in state_vars_idx]

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------

    res_state_vars_view =
        @view res[state_var_Idx_in_state]

    du_state_vars_view =
        @view dx[state_var_Idx_in_state]

    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------

    res_algebraic_vars_view =
        @view res[algebraic_wt_fault_Idx_in_state]

    du_algebraic_vars_view =
        @view dx[algebraic_wt_fault_Idx_in_state]

    u_algebraic_vars_view =
        @view x[algebraic_wt_fault_Idx_in_state]

    #----------------------------------------

    vh =  @view x[vh_Idx_in_state] # n : normal

    θh =  @view x[θh_Idx_in_state]

    gens_i_d =  @view x[id_Idx_in_state]

    gens_i_q =  @view x[iq_Idx_in_state]

    vhf = @view x[vhf_Idx_in_state]
    
    θhf = @view x[θhf_Idx_in_state]
    
    #----------------------------------------

    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]

    #----------------------------------------
    # using plant dynamics 
    #----------------------------------------

    dae_gens_plants_generic_model_func!(
        res_state_vars_view,
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            dae_plants_kwd_para )


    #----------------------------------------
    # power flow equations
    #----------------------------------------

    gens_δ = x[
        δ_idx_in_state]

    gens_ed_dash = x[
        ed_dash_idx_in_state]

    gens_eq_dash = x[
        eq_dash_idx_in_state]

    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        Ynet_pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             Ynet);
             kwd_para =
                 algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
        
        
        Ynet_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             fault_Ynet);
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        Ynet_post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             post_fault_Ynet) ;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


function line_outage_generic_dynamics_wt_pre_post_fault_by_dae_pf_funcs!(
    res,
    dx,
    x,
    (;generic_model_dynamics_para,
     Ynet,
     post_fault_Ynet,
     system_fault_status,
     plants_cb_paras_switches),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, 

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     # dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 
     # pf_generic_gens_para, 
     # Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     dae_plants_kwd_para,

     algebraic_generic_model_kwd_para,

     # cleared_selected_lines_faults_net_para,
     nodes_idx_with_adjacent_nodes_idx,
     pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
     post_clear_fault_nodes_idx_with_adjacent_nodes_idx) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              # :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              # :pf_generic_gens_para,              
              # :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_kwd_para,

              # :cleared_selected_lines_faults_net_para,
              :nodes_idx_with_adjacent_nodes_idx,
              :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
              :post_clear_fault_nodes_idx_with_adjacent_nodes_idx))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))
    
    #----------------------------------------
    
   (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------

    res_states_vars_views =
        [ view(res, idx) for idx in state_vars_idx]
    
    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    res_state_vars_view =
        @view res[state_var_Idx_in_state]
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------
    
    gens_δ = u_state_vars_view[
        δ_idx_in_state]
        
    gens_ed_dash = u_state_vars_view[
        ed_dash_idx_in_state]
    
    gens_eq_dash = u_state_vars_view[
        eq_dash_idx_in_state]    
    
    #----------------------------------------
    
    res_algebraic_vars_view =
        @view res[algebraic_var_Idx_in_state]

    # note du for algebraic_vars = 0
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------

    vh  =  @view x[ vh_Idx_in_state]
    
    θh  =  @view x[ θh_Idx_in_state]
    
    gens_i_d =  @view x[id_Idx_in_state]
    
    gens_i_q =  @view x[iq_Idx_in_state]

    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        θh[
            non_gens_nodes_idx ]
    
    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref =
        generic_model_dynamics_para[
            dyn_v_ref_Idx]
    
    p_order =
        generic_model_dynamics_para[
            dyn_p_order_Idx]
    
    P_non_gens =
        generic_model_dynamics_para[
            dyn_Png_Idx]
    
    Q_non_gens =
        generic_model_dynamics_para[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end


    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    # #----------------------------------------
    # #----------------------------------------
    # # using aggregate plants dynamics 
    # #----------------------------------------
    # #----------------------------------------

    dae_gens_plants_generic_model_func!(
        res_state_vars_view,
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            dae_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        pertubation_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para, Ynet);
            kwd_para =
                algebraic_generic_model_kwd_para,
        nodes_idx_with_adjacent_nodes_idx =
            nodes_idx_with_adjacent_nodes_idx)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        pertubation_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             post_fault_Ynet);
            kwd_para =
                algebraic_generic_model_kwd_para,
        nodes_idx_with_adjacent_nodes_idx =
            post_clear_fault_nodes_idx_with_adjacent_nodes_idx)
           
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end

         
function Ynet_generic_dynamics_wt_pre_post_fault_by_dae_pf_funcs!(
    res,
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     Ynet,
     post_fault_Ynet,
     system_fault_status
     # ,
     # plants_cb_paras_switches
     ),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, 

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     # dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 
     # pf_generic_gens_para, 
     # Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     dae_plants_kwd_para,

     algebraic_generic_model_kwd_para,

     # cleared_selected_lines_faults_net_para,
     nodes_idx_with_adjacent_nodes_idx,
     pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
     post_clear_fault_nodes_idx_with_adjacent_nodes_idx) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              # :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              # :pf_generic_gens_para,              
              # :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_kwd_para,

              # :cleared_selected_lines_faults_net_para,
              :nodes_idx_with_adjacent_nodes_idx,
              :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
              :post_clear_fault_nodes_idx_with_adjacent_nodes_idx))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))
    
    #----------------------------------------
    
   (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------

    res_states_vars_views =
        [ view(res, idx) for idx in state_vars_idx]
    
    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    res_state_vars_view =
        @view res[state_var_Idx_in_state]
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------
    
    gens_δ = u_state_vars_view[
        δ_idx_in_state]
        
    gens_ed_dash = u_state_vars_view[
        ed_dash_idx_in_state]
    
    gens_eq_dash = u_state_vars_view[
        eq_dash_idx_in_state]    
    
    #----------------------------------------
    
    res_algebraic_vars_view =
        @view res[algebraic_var_Idx_in_state]

    # note du for algebraic_vars = 0
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------

    vh  =  @view x[ vh_Idx_in_state]
    
    θh  =  @view x[ θh_Idx_in_state]
    
    gens_i_d =  @view x[id_Idx_in_state]
    
    gens_i_q =  @view x[iq_Idx_in_state]

    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        θh[
            non_gens_nodes_idx ]
    
    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_v_ref_Idx]
    
    p_order =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_p_order_Idx]
    
    P_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Png_Idx]
    
    Q_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end


    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    # #----------------------------------------
    # #----------------------------------------
    # # using aggregate plants dynamics 
    # #----------------------------------------
    # #----------------------------------------

    dae_gens_plants_generic_model_func!(
        res_state_vars_view,
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            dae_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        pertubation_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para, Ynet);
            kwd_para =
                algebraic_generic_model_kwd_para,
        nodes_idx_with_adjacent_nodes_idx =
            nodes_idx_with_adjacent_nodes_idx)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        pertubation_algebraic_generic_pf_ΔPQ_mismatch!(
            res_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             post_fault_Ynet);
            kwd_para =
                algebraic_generic_model_kwd_para,
        nodes_idx_with_adjacent_nodes_idx =
            post_clear_fault_nodes_idx_with_adjacent_nodes_idx)
           
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


# ------------------------------------------------------
# Fault ODE
# ------------------------------------------------------

function mm_generic_dynamics_wt_pre_fault_post_by_ode_pf_funcs!(
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     Ynet,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,     

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_wt_fault_kwd_para,
              
              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    (;state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state,
     vhf_Idx_in_state,
     θhf_Idx_in_state) = NamedTupleTools.select(
         state_vars_and_i_dq_wt_fault_Idx_in_state,
             (:state_var_Idx_in_state,
              :vh_Idx_in_state,
              :θh_Idx_in_state,
              :id_Idx_in_state,
              :iq_Idx_in_state,
              :vhf_Idx_in_state,
              :θhf_Idx_in_state))

    (
     algebraic_wt_fault_Idx_in_state, ) =
         NamedTupleTools.select(
             state_algebraic_vars_wt_fault_Idx_in_state,
             (
              :algebraic_wt_fault_Idx_in_state,))

    #----------------------------------------

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------

    du_state_vars_view =
        @view dx[state_var_Idx_in_state]

    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------

    du_algebraic_vars_view =
        @view dx[algebraic_wt_fault_Idx_in_state]

    u_algebraic_vars_view =
        @view x[algebraic_wt_fault_Idx_in_state]

    #----------------------------------------

    vh =  @view x[vh_Idx_in_state] # n : normal

    θh =  @view x[θh_Idx_in_state]

    gens_i_d =  @view x[id_Idx_in_state]

    gens_i_q =  @view x[iq_Idx_in_state]

    vhf = @view x[vhf_Idx_in_state]
    
    θhf = @view x[θhf_Idx_in_state]
    
    #----------------------------------------

    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]

    #----------------------------------------
    # using plant dynamics 
    #----------------------------------------

    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            dae_plants_kwd_para )


    #----------------------------------------
    # power flow equations
    #----------------------------------------

    gens_δ = x[
        δ_idx_in_state]

    gens_ed_dash = x[
        ed_dash_idx_in_state]

    gens_eq_dash = x[
        eq_dash_idx_in_state]

    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
             kwd_para =
                 algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
   

        fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para ;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


         
function mm_Ynet_generic_dynamics_wt_pre_fault_post_by_ode_pf_funcs!(
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     Ynet,
     fault_Ynet,
     post_fault_Ynet,
     system_fault_status),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,     

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_wt_fault_kwd_para,
              
              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    (;state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state,
     vhf_Idx_in_state,
     θhf_Idx_in_state) = NamedTupleTools.select(
         state_vars_and_i_dq_wt_fault_Idx_in_state,
             (:state_var_Idx_in_state,
              :vh_Idx_in_state,
              :θh_Idx_in_state,
              :id_Idx_in_state,
              :iq_Idx_in_state,
              :vhf_Idx_in_state,
              :θhf_Idx_in_state))

    (
     algebraic_wt_fault_Idx_in_state, ) =
         NamedTupleTools.select(
             state_algebraic_vars_wt_fault_Idx_in_state,
             (
              :algebraic_wt_fault_Idx_in_state,))

    #----------------------------------------

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------

    du_state_vars_view =
        @view dx[state_var_Idx_in_state]

    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------

    du_algebraic_vars_view =
        @view dx[algebraic_wt_fault_Idx_in_state]

    u_algebraic_vars_view =
        @view x[algebraic_wt_fault_Idx_in_state]

    #----------------------------------------

    vh =  @view x[vh_Idx_in_state] # n : normal

    θh =  @view x[θh_Idx_in_state]

    gens_i_d =  @view x[id_Idx_in_state]

    gens_i_q =  @view x[iq_Idx_in_state]

    vhf = @view x[vhf_Idx_in_state]
    
    θhf = @view x[θhf_Idx_in_state]
    
    #----------------------------------------

    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]

    #----------------------------------------
    # using plant dynamics 
    #----------------------------------------

    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            ode_plants_kwd_para )


    #----------------------------------------
    # power flow equations
    #----------------------------------------

    gens_δ = x[
        δ_idx_in_state]

    gens_ed_dash = x[
        ed_dash_idx_in_state]

    gens_eq_dash = x[
        eq_dash_idx_in_state]

    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        Ynet_pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             Ynet);
             kwd_para =
                 algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
        
        
        Ynet_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             fault_Ynet);
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        Ynet_post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             post_fault_Ynet) ;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


         
function mm_line_loss_generic_dynamics_wt_pre_fault_post_by_ode_pf_funcs!(
    dx,
    x,
    (;generic_model_dynamics_para,
     Ynet,
     fault_Ynet,
     post_fault_Ynet,
     system_fault_status,
     plants_cb_paras_switches),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_vars_and_i_dq_wt_fault_Idx_in_state,
     
     state_algebraic_vars_Idx_in_state,
     state_algebraic_vars_wt_fault_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,
     
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
     
     dyn_pf_vh_vhf_Idx,
     dyn_pf_θh_θhf_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_wt_fault_kwd_para,     

     algebraic_generic_model_wt_fault_sol_kwd_para,

     no_lines_fault,
     no_current_lines_fault,

     with_faults) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              
              :state_vars_and_i_dq_wt_fault_Idx_in_state,
              
              :state_algebraic_vars_Idx_in_state,
              :state_algebraic_vars_wt_fault_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
              :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,
              
              :dyn_pf_vh_vhf_Idx,
              :dyn_pf_θh_θhf_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_wt_fault_kwd_para,
              
              :algebraic_generic_model_wt_fault_sol_kwd_para,
              
              :no_lines_fault,
              :no_current_lines_fault,

              :with_faults))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref = generic_model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = generic_model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    #----------------------------------------

    (;state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state,
     vhf_Idx_in_state,
     θhf_Idx_in_state) = NamedTupleTools.select(
         state_vars_and_i_dq_wt_fault_Idx_in_state,
             (:state_var_Idx_in_state,
              :vh_Idx_in_state,
              :θh_Idx_in_state,
              :id_Idx_in_state,
              :iq_Idx_in_state,
              :vhf_Idx_in_state,
              :θhf_Idx_in_state))

    (
     algebraic_wt_fault_Idx_in_state, ) =
         NamedTupleTools.select(
             state_algebraic_vars_wt_fault_Idx_in_state,
             (
              :algebraic_wt_fault_Idx_in_state,))

    #----------------------------------------

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------

    du_state_vars_view =
        @view dx[state_var_Idx_in_state]

    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------

    du_algebraic_vars_view =
        @view dx[algebraic_wt_fault_Idx_in_state]

    u_algebraic_vars_view =
        @view x[algebraic_wt_fault_Idx_in_state]

    #----------------------------------------

    vh =  @view x[vh_Idx_in_state] # n : normal

    θh =  @view x[θh_Idx_in_state]

    gens_i_d =  @view x[id_Idx_in_state]

    gens_i_q =  @view x[iq_Idx_in_state]

    vhf = @view x[vhf_Idx_in_state]
    
    θhf = @view x[θhf_Idx_in_state]
    
    #----------------------------------------

    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]

    #----------------------------------------
    # using plant dynamics 
    #----------------------------------------

    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            ode_plants_kwd_para )


    #----------------------------------------
    # power flow equations
    #----------------------------------------

    gens_δ = x[
        δ_idx_in_state]

    gens_ed_dash = x[
        ed_dash_idx_in_state]

    gens_eq_dash = x[
        eq_dash_idx_in_state]

    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        Ynet_pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             Ynet);
             kwd_para =
                 algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 1

        @show system_fault_status[1]
        
        
        Ynet_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             fault_Ynet);
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        Ynet_post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             post_fault_Ynet) ;
            kwd_para =
                algebraic_generic_model_wt_fault_kwd_para)
        
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end

#---------------------------------------------------

function mm_line_outage_generic_dynamics_wt_pre_post_fault_by_ode_pf_funcs!(
    dx,
    x,
    (;generic_model_dynamics_para,
     Ynet,
     post_fault_Ynet,
     system_fault_status,
     plants_cb_paras_switches),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, 

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     # dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 
     # pf_generic_gens_para, 
     # Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     ode_plants_kwd_para,

     algebraic_generic_model_kwd_para,

     # cleared_selected_lines_faults_net_para,
     nodes_idx_with_adjacent_nodes_idx,
     pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
     post_clear_fault_nodes_idx_with_adjacent_nodes_idx) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              # :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              # :pf_generic_gens_para,              
              # :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_kwd_para,

              # :cleared_selected_lines_faults_net_para,
              :nodes_idx_with_adjacent_nodes_idx,
              :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
              :post_clear_fault_nodes_idx_with_adjacent_nodes_idx))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))
    
    #----------------------------------------
    
   (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    
    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------
    
    gens_δ = u_state_vars_view[
        δ_idx_in_state]
        
    gens_ed_dash = u_state_vars_view[
        ed_dash_idx_in_state]
    
    gens_eq_dash = u_state_vars_view[
        eq_dash_idx_in_state]    
    
    #----------------------------------------

    # note du for algebraic_vars = 0
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------

    vh  =  @view x[ vh_Idx_in_state]
    
    θh  =  @view x[ θh_Idx_in_state]
    
    gens_i_d =  @view x[id_Idx_in_state]
    
    gens_i_q =  @view x[iq_Idx_in_state]

    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        θh[
            non_gens_nodes_idx ]
    
    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref =
        generic_model_dynamics_para[
            dyn_v_ref_Idx]
    
    p_order =
        generic_model_dynamics_para[
            dyn_p_order_Idx]
    
    P_non_gens =
        generic_model_dynamics_para[
            dyn_Png_Idx]
    
    Q_non_gens =
        generic_model_dynamics_para[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    # #----------------------------------------
    # #----------------------------------------
    # # using aggregate plants dynamics 
    # #----------------------------------------
    # #----------------------------------------

    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            ode_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        pertubation_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para, Ynet);
            kwd_para =
                algebraic_generic_model_kwd_para,
        nodes_idx_with_adjacent_nodes_idx =
            nodes_idx_with_adjacent_nodes_idx)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        pertubation_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             post_fault_Ynet);
            kwd_para =
                algebraic_generic_model_kwd_para,
        nodes_idx_with_adjacent_nodes_idx =
            post_clear_fault_nodes_idx_with_adjacent_nodes_idx)
           
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end

         
function mm_Ynet_generic_dynamics_wt_pre_post_fault_by_ode_pf_funcs!(
    dx,
    x,
    (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     Ynet,
     post_fault_Ynet,
     system_fault_status
     # ,
     # plants_cb_paras_switches
     ),
    t;
    kwd_para =
        generic_system_dynamics_wt_fault_kwd_para  )

    #----------------------------------------


    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, 

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     # dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 
     # pf_generic_gens_para, 
     # Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     ode_plants_kwd_para,

     algebraic_generic_model_kwd_para,

     # cleared_selected_lines_faults_net_para,
     nodes_idx_with_adjacent_nodes_idx,
     pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
     post_clear_fault_nodes_idx_with_adjacent_nodes_idx) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              # :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              # :pf_generic_gens_para,              
              # :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_kwd_para,

              # :cleared_selected_lines_faults_net_para,
              :nodes_idx_with_adjacent_nodes_idx,
              :pre_clear_fault_nodes_idx_with_adjacent_nodes_idx,
              :post_clear_fault_nodes_idx_with_adjacent_nodes_idx))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))
    
    #----------------------------------------
    
   (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------
    
    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------
    
    gens_δ = u_state_vars_view[
        δ_idx_in_state]
        
    gens_ed_dash = u_state_vars_view[
        ed_dash_idx_in_state]
    
    gens_eq_dash = u_state_vars_view[
        eq_dash_idx_in_state]    
    
    #----------------------------------------

    # note du for algebraic_vars = 0
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------

    vh  =  @view x[ vh_Idx_in_state]
    
    θh  =  @view x[ θh_Idx_in_state]
    
    gens_i_d =  @view x[id_Idx_in_state]
    
    gens_i_q =  @view x[iq_Idx_in_state]

    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        θh[
            non_gens_nodes_idx ]
    
    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_v_ref_Idx]
    
    p_order =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_p_order_Idx]
    
    P_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Png_Idx]
    
    Q_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end


    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    # #----------------------------------------
    # #----------------------------------------
    # # using aggregate plants dynamics 
    # #----------------------------------------
    # #----------------------------------------

    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            ode_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    #----------------------------------------
    #----------------------------------------

    if system_fault_status[1] == 0

        @show system_fault_status[1]

        pertubation_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para, Ynet);
            kwd_para =
                algebraic_generic_model_kwd_para,
        nodes_idx_with_adjacent_nodes_idx =
            nodes_idx_with_adjacent_nodes_idx)
        
    elseif system_fault_status[1] == 2

        @show system_fault_status[1]

        pertubation_algebraic_generic_pf_ΔPQ_mismatch!(
            du_algebraic_vars_view,
            u_algebraic_vars_view,
            (δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para,
             post_fault_Ynet);
            kwd_para =
                algebraic_generic_model_kwd_para,
        nodes_idx_with_adjacent_nodes_idx =
            post_clear_fault_nodes_idx_with_adjacent_nodes_idx)
           
    else

        throw("system_fault_status " *
            "$(system_fault_status[1]) unknown")
        
    end
        
    return nothing

end


# ------------------------------------------------------
# comments
# ------------------------------------------------------
