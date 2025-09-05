# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


# https://stackoverflow.com/questions/58033504/julia-delete-rows-and-columns-from-an-array-or-matix


########################################################
#-------------------------------------------------------
# Begining of functions defiitions
#-------------------------------------------------------
########################################################


#-------------------------------------------------------
#-------------------------------------------------------


function get_Inet_inj_by_Bnj(
    uh,
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx )

    return [ sum(
        [ im * imag(ynj) * vj for (ynj, vj) in
             zip(Y_bus_vec,
                 uh[nth_node_idx_and_adj_nodes_idx] ) ])
             for (Y_bus_vec,
                  nth_node_idx_and_adj_nodes_idx ) in
                 zip( Ynet,
                      nodes_node_idx_and_incident_edges_other_node_idx ) ]
    
end


function get_Inet_inj_by_Bnj(
    uh,
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx )

    t_nodes_idx_with_adjacent_nodes_idx =
        [[n2s_all_nodes_idx[idx] for idx in
              nth_node_idx_and_adj_nodes_idx]
         for nth_node_idx_and_adj_nodes_idx in
             nodes_idx_with_adjacent_nodes_idx  ]

    return [ sum(
        [ im * imag(ynj) * vj for (ynj, vj) in
             zip(Y_bus_vec,
                 uh[nth_node_idx_and_adj_nodes_idx ] ) ])
             for (Y_bus_vec,
                  nth_node_idx_and_adj_nodes_idx ) in
                 zip( Ynet,
                      t_nodes_idx_with_adjacent_nodes_idx)]
    

    # return [ sum(
    #     [ im * imag(ynj) * vj for (ynj, vj) in
    #          zip(Y_bus_vec,
    #              uh[ n2s_all_nodes_idx[idx]
    #                  for idx in
    #                      nth_node_idx_and_adj_nodes_idx])])
    #          for (Y_bus_vec,
    #               nth_node_idx_and_adj_nodes_idx ) in
    #              zip( Ynet,
    #                   nodes_idx_with_adjacent_nodes_idx ) ]
    
end



function get_Inet_inj(
    uh,
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx )

    return [ sum(
        [ ynj * vj for (ynj, vj) in
             zip(Y_bus_vec,
                 uh[nth_node_idx_and_adj_nodes_idx] ) ])
             for (Y_bus_vec,
                  nth_node_idx_and_adj_nodes_idx ) in
                 zip( Ynet,
                      nodes_node_idx_and_incident_edges_other_node_idx ) ]
    
end



function get_Inet_inj(
    uh,
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx )

    t_nodes_idx_with_adjacent_nodes_idx =
        [[n2s_all_nodes_idx[idx] for idx in
              nth_node_idx_and_adj_nodes_idx]
         for nth_node_idx_and_adj_nodes_idx in
             nodes_idx_with_adjacent_nodes_idx  ]

    return [ sum(
        [ ynj * vj for (ynj, vj) in
             zip(Y_bus_vec,
                 uh[nth_node_idx_and_adj_nodes_idx ] ) ])
             for (Y_bus_vec,
                  nth_node_idx_and_adj_nodes_idx ) in
                 zip( Ynet,
                      t_nodes_idx_with_adjacent_nodes_idx)]
    

    # return [ sum(
    #     [ ynj * vj for (ynj, vj) in
    #          zip(Y_bus_vec,
    #              uh[ n2s_all_nodes_idx[idx]
    #                  for idx in
    #                      nth_node_idx_and_adj_nodes_idx])])
    #          for (Y_bus_vec,
    #               nth_node_idx_and_adj_nodes_idx ) in
    #              zip( Ynet,
    #                   nodes_idx_with_adjacent_nodes_idx ) ]
    
end



function get_Iinj(
    uh,
    P_Q_non_gens_view,
    P_Q_gens_loc_load_view,
    Inet_inj  )


    S_gens_loc_load = x_from_xr_xi.(
        P_Q_gens_loc_load_view )

    S_non_gens = x_from_xr_xi.( P_Q_non_gens_view )
    
    return Inet_inj + (conj.(S_non_gens ))./ (
        conj.( uh )) + (conj.(S_gens_loc_load ))./ (
            conj.( uh ))

end



function get_Iinj(
    uh,
    
    P_non_gens,
    Q_non_gens,
    
    P_g_loc_load,
    Q_g_loc_load,
    
    pf_kw_nodes_types_idxs,
    pf_kw_n2s_idxs,
    
    loc_load_exist,
    Inet_inj;
    system_status =
        nothing)
    
    P_non_gens =
        deepcopy(P_non_gens)
    
    Q_non_gens =
        deepcopy(Q_non_gens)

    (;gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx 
     ) =
         NamedTupleTools.select(
             pf_kw_nodes_types_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))
    
    (;n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             pf_kw_n2s_idxs,
             (:n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx ))

    #----------------------------------------
    
    all_nodes_idx =
        sort([gens_nodes_idx;
              non_gens_nodes_idx])


    if system_status == :fault_state

        (; 
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_with_loc_load_idx,
         all_nodes_idx,

         fault_nodes_idx
         ) =
             NamedTupleTools.select(
                 pf_kw_nodes_types_idxs,
                 (
                  :gens_nodes_idx,
                  :non_gens_nodes_idx,
                  :gens_with_loc_load_idx,
                   :all_nodes_idx,

                   :fault_nodes_idx))

        (; 
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx,

         n2s_fault_nodes_idx) =
             NamedTupleTools.select(
                 pf_kw_n2s_idxs,
                 (
                  :n2s_gens_idx,
                  :n2s_non_gens_idx,
                  :n2s_gens_with_loc_load_idxs,
                  :n2s_all_nodes_idx,

                     :n2s_fault_nodes_idx))

        P_non_gens =
            [P_non_gens;
             zeros(length(fault_nodes_idx))]
        
        Q_non_gens =
            [Q_non_gens;
             zeros(length(fault_nodes_idx))]

    end
    
        
    #----------------------------------------
        
    S_non_gens =
        P_non_gens +
        im * Q_non_gens
    
    if loc_load_exist == true

        S_gens_loc_load =
            P_g_loc_load +
            im * Q_g_loc_load
    end

    return [ idx ∈ gens_with_loc_load_idx ?
        Inet_inj[ n2s_all_nodes_idx[idx]] +
        conj(S_gens_loc_load[
            n2s_gens_with_loc_load_idxs[idx] ] ) /
                conj(uh[ n2s_all_nodes_idx[idx]]) :
                idx ∈  non_gens_nodes_idx ?
                Inet_inj[ n2s_all_nodes_idx[idx]] +
                conj(S_non_gens[
                    n2s_non_gens_idx[idx] ] ) /
                        conj(uh[ n2s_all_nodes_idx[idx]]) :
                        Inet_inj[ n2s_all_nodes_idx[ idx ] ]
             for idx in all_nodes_idx ]    
end



function get_GenSinj(
    Sbus_n,
    
    P_non_gens,
    Q_non_gens,
    
    P_g_loc_load,
    Q_g_loc_load,
    
    pf_kw_nodes_types_idxs,
    pf_kw_n2s_idxs,    
    loc_load_exist;
    system_status =
        nothing)
    
    P_non_gens =
        deepcopy(P_non_gens)
    
    Q_non_gens =
        deepcopy(Q_non_gens)
    
    # (; slack_gens_nodes_idx,
    #  non_slack_gens_nodes_idx,
    #  gens_nodes_idx,
    #  non_gens_nodes_idx,
    #  gens_with_loc_load_idx,
    #  all_nodes_idx #
    #  ) =
    #      pf_kw_nodes_types_idxs
    
    # (; n2s_slack_gens_idx,
    #  n2s_non_slack_gens_idx,
    #  n2s_gens_idx,
    #  n2s_non_gens_idx,
    #  n2s_gens_with_loc_load_idxs,
    #  n2s_all_nodes_idx ) =
    #      pf_kw_n2s_idxs


    (; #slack_gens_nodes_idx,
     #non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     all_nodes_idx 
     ) =
         NamedTupleTools.select(
             pf_kw_nodes_types_idxs,
             ( # :slack_gens_nodes_idx,
              # :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_with_loc_load_idx,
              :all_nodes_idx))
    
    (; # n2s_slack_gens_idx,
     # n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             pf_kw_n2s_idxs,
             ( # :n2s_slack_gens_idx,
               # :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_all_nodes_idx ))
    
    #----------------------------------------


    if system_status == :fault_state

        (; 
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_with_loc_load_idx,
         all_nodes_idx,

         fault_nodes_idx
         ) =
             NamedTupleTools.select(
                 pf_kw_nodes_types_idxs,
                 (
                  :gens_nodes_idx,
                  :non_gens_nodes_idx,
                  :gens_with_loc_load_idx,
                   :all_nodes_idx,

                   :fault_nodes_idx))

        (; 
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx,

         n2s_fault_nodes_idx) =
             NamedTupleTools.select(
                 pf_kw_n2s_idxs,
                 (
                  :n2s_gens_idx,
                  :n2s_non_gens_idx,
                  :n2s_gens_with_loc_load_idxs,
                  :n2s_all_nodes_idx,

                     :n2s_fault_nodes_idx))

        P_non_gens =
            [P_non_gens;
             zeros(length(fault_nodes_idx))]
        
        Q_non_gens =
            [Q_non_gens;
             zeros(length(fault_nodes_idx))]

    end
    
        
    #----------------------------------------
    
    # all_nodes_idx = sort(
    #     [gens_nodes_idx;
    #      non_gens_nodes_idx])
    
    S_non_gens =
        P_non_gens +
        im * Q_non_gens
    
    if loc_load_exist == true

        S_gens_loc_load =
            P_g_loc_load +
            im * Q_g_loc_load
    end
    

    # return [ idx ∈ gens_with_loc_load_idx ?
    #     Sbus_n[idx] +
    #     S_gens_loc_load[
    #         n2s_gens_with_loc_load_idxs[idx]] :
                
    #             idx ∈  non_gens_nodes_idx ?
    #             Sbus_n[idx] +
    #             S_non_gens[
    #                 n2s_non_gens_idx[idx]] : Sbus_n[idx]
    #          for idx in
    #              1:length(Sbus_n) ]


    return [ idx ∈ gens_with_loc_load_idx ?
        Sbus_n[ n2s_all_nodes_idx[idx]] +
        S_gens_loc_load[
            n2s_gens_with_loc_load_idxs[idx]] :
                
                idx ∈  non_gens_nodes_idx ?
                Sbus_n[ n2s_all_nodes_idx[idx]] +
                S_non_gens[
                    n2s_non_gens_idx[idx]] :
                        Sbus_n[ n2s_all_nodes_idx[idx]]
             for idx in
                 all_nodes_idx ]

end


function get_Iinj(
    uh,
    P_Q_non_gens_view,
    P_Q_gens_loc_load_view,
    ( gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx ),
    Inet_inj;
    loc_load_exist = false  )


    if loc_load_exist == true

        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )

        S_gens_loc_load =
            x_from_xr_xi.( P_Q_gens_loc_load_view )

        for (idx, S_node_type) in zip(            
            [ non_gens_nodes_idx,
             gens_with_loc_load_idx ],[S_non_gens,
             S_gens_loc_load ] )

            Inet_inj[idx] .=
                Inet_inj[idx] +
               ((conj.(S_node_type)) ./(conj.(uh[idx])) )
        end
    else
        
        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )

        Inet_inj[ non_gens_nodes_idx ] .=
                Inet_inj[ non_gens_nodes_idx ] +
                ((conj.(S_non_gens)) ./ (
                    conj.(uh[non_gens_nodes_idx]))) 
        
    end

    return Inet_inj
    
end

#-------------------------------------------------------
#-------------------------------------------------------


function get_Ifrom_Ito_by_Bnj(
    nodes_pf_U_view,
    edges_Ybr_cal,
    edges_orientation)

    # the matrix in view y_π, is extracted with [1]
    
    return [ im * imag.(y_π) * [ x_from_xr_xi(
        nodes_pf_U_view[
            orient[1] ]),
                     x_from_xr_xi( nodes_pf_U_view[
                         orient[2] ])]
             for (y_π, orient ) in
                 zip(edges_Ybr_cal,
                     edges_orientation)]

end


function get_Ifrom_Ito_by_Bnj(
    nodes_pf_U_view,
    edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx)

    # the matrix in view y_π, is extracted with [1]
    
    return [ im * imag.(y_π) * [ x_from_xr_xi(
        nodes_pf_U_view[
            n2s_all_nodes_idx[ orient[1]] ]),
                     x_from_xr_xi( nodes_pf_U_view[
                         n2s_all_nodes_idx[ orient[2]]])]
             for (y_π, orient ) in
                 zip(edges_Ybr_cal,
                     edges_orientation)]

end



function get_I_from_I_to_by_Bnj(
    uh,
    edges_Ybr_cal,
    edges_orientation )
    
    return [ im * imag.(y_π) * [uh[orient[1]], uh[orient[2]] ]
             for (y_π, orient ) in
                 zip(edges_Ybr_cal,
                     edges_orientation)]
end


function get_I_from_I_to_by_Bnj(
    uh,
    edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx )
    
    return [im * imag.(y_π) * [uh[ n2s_all_nodes_idx[orient[1]]],
                     uh[ n2s_all_nodes_idx[orient[2]]] ]
             for (y_π, orient ) in
                 zip(edges_Ybr_cal,
                     edges_orientation)]
end


function get_Ifrom_Ito(
    nodes_pf_U_view,
    edges_Ybr_cal,
    edges_orientation)

    # the matrix in view y_π, is extracted with [1]
    
    return [ y_π * [ x_from_xr_xi(
        nodes_pf_U_view[
            orient[1] ]),
                     x_from_xr_xi( nodes_pf_U_view[
                         orient[2] ])]
             for (y_π, orient ) in
                 zip(edges_Ybr_cal,
                     edges_orientation)]

end


function get_Ifrom_Ito(
    nodes_pf_U_view,
    edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx)

    # the matrix in view y_π, is extracted with [1]
    
    return [ y_π * [ x_from_xr_xi(
        nodes_pf_U_view[
            n2s_all_nodes_idx[ orient[1]] ]),
                     x_from_xr_xi( nodes_pf_U_view[
                         n2s_all_nodes_idx[ orient[2]]])]
             for (y_π, orient ) in
                 zip(edges_Ybr_cal,
                     edges_orientation)]

end



function get_I_from_I_to(
    uh,
    edges_Ybr_cal,
    edges_orientation )
    
    return [ y_π * [ uh[orient[1]], uh[orient[2]] ]
             for (y_π, orient ) in
                 zip(edges_Ybr_cal,
                     edges_orientation)]
end


function get_I_from_I_to(
    uh,
    edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx )
    
    return [ y_π * [ uh[ n2s_all_nodes_idx[orient[1]]],
                     uh[ n2s_all_nodes_idx[orient[2]]] ]
             for (y_π, orient ) in
                 zip(edges_Ybr_cal,
                     edges_orientation)]
end


function get_I_from_I_to(
    vh, θh,
    edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx )
    
    uh = vh .* exp.(im * θh)
    
    return [ y_π * [ uh[ n2s_all_nodes_idx[orient[1]]],
                     uh[ n2s_all_nodes_idx[orient[2]]] ]
             for (y_π, orient ) in
                 zip(edges_Ybr_cal,
                     edges_orientation)]
end


function get_a_I_from_I_to(
    a_from_vh, a_to_vh,
    a_from_θh, a_to_θh,
    y_π )
    
    return y_π * [a_from_vh * exp(im * a_from_θh),
                    a_to_vh * exp(im * a_to_θh) ]

end


function get_a_I_from_I_to(
    a_from_vh_a_to_vh,
    a_from_θh_a_to_θh,
    y_π )
    
    return y_π * (a_from_vh_a_to_vh .*
        exp.(im * a_from_θh_a_to_θh))
end


function get_a_branch_current(
    a_from_vh, a_to_vh,
    a_from_θh, a_to_θh,
    y_π)

    I_from_I_to =
        get_a_I_from_I_to(
            a_from_vh, a_to_vh,
            a_from_θh, a_to_θh,
            y_π )

    return first( I_from_I_to ) +
        second( I_from_I_to )

end


function get_a_branch_current(
    a_from_vh_a_to_vh,
    a_from_θh_a_to_θh,
    y_π )

    I_from_I_to =
        get_a_I_from_I_to(
            a_from_vh_a_to_vh,
            a_from_θh_a_to_θh,
            y_π )
    
    return first( I_from_I_to ) +
        second( I_from_I_to )
end


function get_I_from_I_to_by_vh_θh(
    vh, θh ;
    edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx )
    
    uh = vh .* exp.(im * θh)

    from_idxs = [ n2s_all_nodes_idx[idx]
                  for idx in
                      first.(edges_orientation)]
    
    to_idxs = [ n2s_all_nodes_idx[idx]
                for idx in
                    second.(edges_orientation)]

    from_uh = uh[from_idxs]
    
    to_uh   = uh[to_idxs]

    return map((y, f_u, t_u) -> y * [f_u, t_u],
               edges_Ybr_cal, from_uh, to_uh )

    # @show uh
    
    # from_idxs = first.(edges_orientation)
    
    # to_idxs   = second.(edges_orientation)
    
    # from_to_uh = vcat(from_uh, to_uh )

    # return map((y,u) -> y * u, edges_Ybr_cal, from_to_uh)
    
    # return [ y_π * [ uh[ n2s_all_nodes_idx[orient[1]]],
    #                  uh[ n2s_all_nodes_idx[orient[2]]] ]
    #          for (y_π, orient ) in
    #              zip(edges_Ybr_cal,
    #                  edges_orientation)]
end


function get_branches_currents_by_vh_θh(
    vh, θh; edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx )

    I_from_I_to =
        get_I_from_I_to_by_vh_θh(
            vh, θh;
            edges_Ybr_cal,
            edges_orientation,
            n2s_all_nodes_idx )
    
    
    return first.(I_from_I_to) + second.(I_from_I_to)
end


function get_vec_I_from_I_to_by_vh_θh(
    vh, θh ;
    edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx )
    
    uh = vh .* exp.( im * θh )

    from_idxs = [ n2s_all_nodes_idx[idx]
                  for idx in
                      first.(edges_orientation)]
    
    to_idxs = [ n2s_all_nodes_idx[idx]
                for idx in
                    second.(edges_orientation)]

    from_uh = uh[from_idxs, :]
    
    to_uh = uh[to_idxs, :]

     return [ map( ( y, f_u, t_u) -> y * [f_u, t_u ],
                   edges_Ybr_cal,
                   a_from_uh, a_to_uh) 
             for (a_from_uh, a_to_uh) in
                 zip(from_uh, to_uh)  ]

    # @show uh
    
    # from_idxs = first.(edges_orientation)
    
    # to_idxs = second.(edges_orientation)
    

end



function get_vec_branches_currents_by_vh_θh(
    vh, θh; edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx )

    I_from_I_to =
        get_vec_I_from_I_to_by_vh_θh(
            vh, θh;
            edges_Ybr_cal,
            edges_orientation,
            n2s_all_nodes_idx )
    
    
    return [ map( (a_if, a_it) -> a_if .+ a_it,
                 a_if_vec, a_it_vec )
            for (a_if_vec, a_it_vec ) in
                zip(first.(first.(I_from_I_to)),
                    second.(first.(I_from_I_to)))]
    # _, ncol_If_It = size(I_from_I_to)
    
end


function get_vec_branches_from_currents_by_vh_θh(
    vh, θh; edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx )

    I_from_I_to =
        get_vec_I_from_I_to_by_vh_θh(
            vh, θh;
            edges_Ybr_cal,
            edges_orientation,
            n2s_all_nodes_idx )
    
    
    return [ map( (a_if ) -> a_if,
                 a_if_vec )
            for a_if_vec  in
                first.(first.(I_from_I_to)) ]
    # _, ncol_If_It = size(I_from_I_to)
    
end


function get_vec_branches_to_currents_by_vh_θh(
    vh, θh; edges_Ybr_cal,
    edges_orientation,
    n2s_all_nodes_idx )

    I_from_I_to =
        get_vec_I_from_I_to_by_vh_θh(
            vh, θh;
            edges_Ybr_cal,
            edges_orientation,
            n2s_all_nodes_idx )
    
    
    return [ map( ( a_it ) -> a_it,
                 a_it_vec )
            for a_it_vec in
                second.(first.(I_from_I_to)) ]
    # _, ncol_If_It = size(I_from_I_to)
    
end


# -------------------------------------------
# gens local load current
# -------------------------------------------

function get_a_gen_loc_load_current(
    gen_vh,
    gen_θh,           
    gen_P_loc_load,
    gen_Q_loc_load )

    return  conj( gen_P_loc_load +
        im * gen_Q_loc_load ) / (
            gen_vh  * exp( -im * gen_θh ) ) 
end



function get_a_gen_loc_load_current(
    gen_vh,
    gen_θh,
    
    gen_P_loc_load,
    gen_Q_loc_load,
    
    gens_with_loc_load_idx,
    
    gen_idx )

    return gen_idx ∈ gens_with_loc_load_idx ?
        get_a_gen_loc_load_current(
            gen_vh,
            gen_θh,
            gen_P_loc_load,
            gen_Q_loc_load ) :
                0.0 + im *0.0
end


function get_a_gen_loc_load_current(
    gen_vh,
    gen_θh,           
    P_g_loc_loads,
    Q_g_loc_loads,
    
    n2s_gens_with_loc_load_idxs,

    gens_with_loc_load_idx,
    gen_idx )
    
    return gen_idx ∈ gens_with_loc_load_idx ?
        conj( P_g_loc_loads[
            n2s_gens_with_loc_load_idxs[
                gen_idx ] ] +
                    im * Q_g_loc_loads[
                        n2s_gens_with_loc_load_idxs[
                            gen_idx ] ]) /
                                (gen_vh * exp(-im *
                                gen_θh)) :
                                0.0 + im * 0.0 
end


# function get_gens_loc_load_current(
#     nodes_vh,
#     nodes_θh,           
#     P_g_loc_loads,
#     Q_g_loc_loads,

#     n2s_gens_with_loc_load_idxs,
#     n2s_all_nodes_idx, #
    
#     gens_nodes_idx,
#     gens_with_loc_load_idx )
    
#     return [get_a_gen_loc_load_current(
#         nodes_vh[ n2s_all_nodes_idx[gen_idx] ],
#         nodes_θh[ n2s_all_nodes_idx[gen_idx] ],
#         P_g_loc_loads,
#         Q_g_loc_loads,
                
#         n2s_gens_with_loc_load_idxs,

#         gens_with_loc_load_idx,
#         gen_idx )
#             for gen_idx in
#                 gens_nodes_idx] 

# end



function get_gens_loc_load_current(
    nodes_vh,
    nodes_θh,
    
    P_g_loc_loads,
    Q_g_loc_loads,

    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx, #
    
    gens_nodes_idx,
    gens_with_loc_load_idx )
    
    return [ gen_idx ∈ gens_with_loc_load_idx ?
        get_a_gen_loc_load_current(
            nodes_vh[ n2s_all_nodes_idx[gen_idx] ],
            nodes_θh[ n2s_all_nodes_idx[gen_idx] ],
            
            P_g_loc_loads[
                n2s_gens_with_loc_load_idxs[gen_idx]],
            Q_g_loc_loads[
                n2s_gens_with_loc_load_idxs[gen_idx]]) :
                    0.0 + im * 0.0 
             for gen_idx in
                gens_nodes_idx] 

end


function get_gens_local_load_current(
    uh,
    P_Q_gens_loc_load_view,
    gens_with_loc_load_idx;
    loc_load_exist = false ) 

    if loc_load_exist == true

        return ( conj.( x_from_xr_xi.(
            P_Q_gens_loc_load_view ) )) ./ (
                conj.( uh[gens_with_loc_load_idx] ))
    else
        return 0
    end
    
end


function get_gens_local_load_current(
    uh,
    P_Q_gens_loc_load_view ) 

    S_gens_loc_load = x_from_xr_xi.(
        P_Q_gens_loc_load_view )


    return (conj.( S_gens_loc_load )) ./ ( conj.( uh ))
end


# -------------------------------------------
# gens nodes current injection
# -------------------------------------------


function get_a_gen_current_injection(
    nodes_vh,
    nodes_θh,
    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
        
    P_g_loc_load,
    Q_g_loc_load,
    
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx, #
    
    gens_with_loc_load_idx,
    
    gen_idx)

    return gen_idx ∉ gens_with_loc_load_idx ?
        get_a_node_∑_ynj_x_vj(
            nodes_vh,
            nodes_θh,
            Ynet_nth_row,
            nth_node_idx_and_adj_nodes_idx,
            n2s_all_nodes_idx ) :
                get_a_node_∑_ynj_x_vj(
                    nodes_vh,
                    nodes_θh,
                    Ynet_nth_row,
                    nth_node_idx_and_adj_nodes_idx,
                    n2s_all_nodes_idx ) +
                        get_a_gen_loc_load_current(
                            nodes_vh[
                                n2s_all_nodes_idx[
                                    gen_idx]],
                            nodes_θh[
                                n2s_all_nodes_idx[
                                    gen_idx]],
                            
                            P_g_loc_load,
                            Q_g_loc_load,
                            
                            n2s_gens_with_loc_load_idxs,
                            gens_with_loc_load_idx,
                            
                            gen_idx )
                    

end


# function get_a_gen_current_injection(
#     nodes_vh,
#     nodes_θh,           
#     Ynet_nth_row,
#     nth_node_idx_and_adj_nodes_idx,
        
#     P_g_loc_load,
#     Q_g_loc_load,
    
#     n2s_gens_with_loc_load_idxs,
#     n2s_all_nodes_idx, #
#     gens_with_loc_load_idx,

#     gen_idx)

#     return gen_idx ∉ gens_with_loc_load_idx ?
#         get_a_node_∑_ynj_x_vj(
#             nodes_vh,
#             nodes_θh,
#             Ynet_nth_row,
#             nth_node_idx_and_adj_nodes_idx,
#             n2s_all_nodes_idx ) :
#                 get_a_node_∑_ynj_x_vj(
#                     nodes_vh,
#                     nodes_θh,
#                     Ynet_nth_row,
#                     nth_node_idx_and_adj_nodes_idx,
#                     n2s_all_nodes_idx ) +
#                         get_a_gen_loc_load_current(
#                             nodes_vh[
#                                 n2s_all_nodes_idx[gen_idx]],
#                             nodes_θh[
#                                 n2s_all_nodes_idx[gen_idx]], 
#                             P_g_loc_load,
#                             Q_g_loc_load,
#                             n2s_gens_with_loc_load_idxs,
#                             gens_with_loc_load_idx,
#                             gen_idx )                    
# end


function get_gens_current_injection(
    nodes_vh,
    nodes_θh,
    
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    
    P_g_loc_load,
    Q_g_loc_load,
        
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx, #
    
    gens_with_loc_load_idx,    
    gens_nodes_idx )

    return [ get_a_gen_current_injection(
        nodes_vh,
        nodes_θh,           
        Ynet[ n2s_all_nodes_idx[gen_idx]],
        nodes_idx_with_adjacent_nodes_idx[
            n2s_all_nodes_idx[gen_idx]],
        
        P_g_loc_load,
        Q_g_loc_load,
        
        n2s_gens_with_loc_load_idxs,
        n2s_all_nodes_idx, #
        
        gens_with_loc_load_idx,
        
        gen_idx)
             for gen_idx in
                 gens_nodes_idx ]

end



# -------------------------------------------
# gens nodes current injection by B_nj
# -------------------------------------------


function get_a_gen_current_injection_by_Bnj(
    nodes_vh,
    nodes_θh,
    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
        
    P_g_loc_load,
    Q_g_loc_load,
    
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx, #
    
    gens_with_loc_load_idx,
    
    gen_idx)

    return gen_idx ∉ gens_with_loc_load_idx ?
        get_a_node_∑_Bnj_x_vj(
            nodes_vh,
            nodes_θh,
            Ynet_nth_row,
            nth_node_idx_and_adj_nodes_idx,
            n2s_all_nodes_idx ) :
                get_a_node_∑_Bnj_x_vj(
                    nodes_vh,
                    nodes_θh,
                    Ynet_nth_row,
                    nth_node_idx_and_adj_nodes_idx,
                    n2s_all_nodes_idx ) +
                        get_a_gen_loc_load_current(
                            nodes_vh[
                                n2s_all_nodes_idx[
                                    gen_idx]],
                            nodes_θh[
                                n2s_all_nodes_idx[
                                    gen_idx]],
                            
                            P_g_loc_load,
                            Q_g_loc_load,
                            
                            n2s_gens_with_loc_load_idxs,
                            gens_with_loc_load_idx,
                            
                            gen_idx )
                    

end


# function get_a_gen_current_injection(
#     nodes_vh,
#     nodes_θh,           
#     Ynet_nth_row,
#     nth_node_idx_and_adj_nodes_idx,
        
#     P_g_loc_load,
#     Q_g_loc_load,
    
#     n2s_gens_with_loc_load_idxs,
#     n2s_all_nodes_idx, #
#     gens_with_loc_load_idx,

#     gen_idx)

#     return gen_idx ∉ gens_with_loc_load_idx ?
#         get_a_node_∑_ynj_x_vj(
#             nodes_vh,
#             nodes_θh,
#             Ynet_nth_row,
#             nth_node_idx_and_adj_nodes_idx,
#             n2s_all_nodes_idx ) :
#                 get_a_node_∑_ynj_x_vj(
#                     nodes_vh,
#                     nodes_θh,
#                     Ynet_nth_row,
#                     nth_node_idx_and_adj_nodes_idx,
#                     n2s_all_nodes_idx ) +
#                         get_a_gen_loc_load_current(
#                             nodes_vh[
#                                 n2s_all_nodes_idx[gen_idx]],
#                             nodes_θh[
#                                 n2s_all_nodes_idx[gen_idx]], 
#                             P_g_loc_load,
#                             Q_g_loc_load,
#                             n2s_gens_with_loc_load_idxs,
#                             gens_with_loc_load_idx,
#                             gen_idx )                    
# end


function get_gens_current_injection_by_Bnj(
    nodes_vh,
    nodes_θh,
    
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    
    P_g_loc_load,
    Q_g_loc_load,
        
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx, #
    
    gens_with_loc_load_idx,    
    gens_nodes_idx )

    return [ get_a_gen_current_injection_by_Bnj(
        nodes_vh,
        nodes_θh,           
        Ynet[ n2s_all_nodes_idx[gen_idx]],
        nodes_idx_with_adjacent_nodes_idx[
            n2s_all_nodes_idx[gen_idx]],
        
        P_g_loc_load,
        Q_g_loc_load,
        
        n2s_gens_with_loc_load_idxs,
        n2s_all_nodes_idx, #
        
        gens_with_loc_load_idx,
        
        gen_idx)
             for gen_idx in
                 gens_nodes_idx ]

end



# -------------------------------------------
# node net power injection from lines based on Yπ_net
# -------------------------------------------


function get_a_node_∑_S_injection_by_Yπ_net(
    vh,
    θh,
    Yπ_net_nth_row,
    Yshunt,
    nth_node_idx,
    nth_node_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in nth_node_adjacent_nodes_idx]
    
    vh_nth_node =
        vh[ transformed_nth_node_idx ]
    
    θh_nth_node =
        θh[ transformed_nth_node_idx ]

    uh_nth_node =
        vh_nth_node * exp(im * θh_nth_node )

    I_shunt_current_nth_node =
        Yshunt[transformed_nth_node_idx] * uh_nth_node

    vk_adj_nodes =
        vh[ transformed_adj_nodes_idx ]
    
    θk_adj_nodes =
        θh[ transformed_adj_nodes_idx ]

    uk_adj_nodes = vk_adj_nodes .* exp.(
        im * θk_adj_nodes)

    ih_ik = [ (yπ * [uh_nth_node, uk_adj_node ])
               for (yπ, uk_adj_node) in
                   zip(Yπ_net_nth_row,
                       uk_adj_nodes)]
    
    uh_uk = [[uh_nth_node, uk_adj_node ]
             for uk_adj_node in
                 uk_adj_nodes ]
    
    sh_sk = [a_uh_uk .* conj.(a_ih_ik)
             for (a_uh_uk, a_ih_ik) in
                 zip(uh_uk, ih_ik)]

    return  sum(first.((sh_sk)) - second.(sh_sk)) +
        uh_nth_node * conj(I_shunt_current_nth_node)
      
end



function get_nodes_∑_S_injection_by_Yπ_net(
    vh,
    θh,
    Yπ_net,
    Yshunt,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx )

    nodes_node_idx =
        first.(nodes_idx_with_adjacent_nodes_idx)
    
    adjacent_nodes_idx =
        get_only_adjacent_nodes_idx(
            nodes_idx_with_adjacent_nodes_idx)
    
    return [ get_a_node_∑_S_injection_by_Yπ_net(
        vh,
        θh,
        Yπ_net_nth_row,
        Yshunt,
        nth_node_idx,
        nth_node_adjacent_nodes_idx,
        n2s_all_nodes_idx )
             for (Yπ_net_nth_row,
                  nth_node_idx,
                  nth_node_adjacent_nodes_idx) in
                 zip(Yπ_net,
                     nodes_node_idx,
                     adjacent_nodes_idx )]
    
end


function get_a_node_∑_Sh_injection_by_Yπ_net(
    vh,
    θh;
    Yπ_net_nth_row,
    Yshunt,
    nth_node_idx,
    nth_node_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in nth_node_adjacent_nodes_idx]
    
    vh_nth_node =
        vh[ transformed_nth_node_idx ]
    
    θh_nth_node =
        θh[ transformed_nth_node_idx ]

    uh_nth_node =
        vh_nth_node * exp(im * θh_nth_node )

    I_shunt_current_nth_node =
        Yshunt[transformed_nth_node_idx] *
        uh_nth_node

    vk_adj_nodes =
        vh[ transformed_adj_nodes_idx ]
    
    θk_adj_nodes =
        θh[ transformed_adj_nodes_idx ]

    uk_adj_nodes = vk_adj_nodes .* exp.(
        im * θk_adj_nodes)

    # ih_ik = [ (yπ[1] * [uh_nth_node, uk_adj_node ])
    #            for (yπ, uk_adj_node) in
    #                zip(Yπ_net_nth_row,
    #                    uk_adj_nodes)]

    ih_ik = [ (yπ * [uh_nth_node, uk_adj_node ])
               for (yπ, uk_adj_node) in
                   zip(Yπ_net_nth_row,
                       uk_adj_nodes)]
    
    uh_uk = [[uh_nth_node, uk_adj_node ]
             for uk_adj_node in
                 uk_adj_nodes ]
    
    sh_sk = [a_uh_uk .* conj.(a_ih_ik)
             for (a_uh_uk, a_ih_ik) in
                 zip(uh_uk, ih_ik)]

    return  sum( first.((sh_sk))  ) +
        uh_nth_node * conj(I_shunt_current_nth_node)


    # return  sum(first.((sh_sk)) + second.(sh_sk)) +
    #     uh_nth_node * conj(I_shunt_current_nth_node)
    
end


function get_nodes_∑_Sh_injection_by_Yπ_net(
    vh,
    θh;
    Yπ_net,
    Yshunt,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx )

    nodes_node_idx =
        first.(nodes_idx_with_adjacent_nodes_idx)
    
    adjacent_nodes_idx =
        get_only_adjacent_nodes_idx(
            nodes_idx_with_adjacent_nodes_idx)
    
    return [ get_a_node_∑_Sh_injection_by_Yπ_net(
        vh,
        θh;
        Yπ_net_nth_row,
        Yshunt,
        nth_node_idx,
        nth_node_adjacent_nodes_idx,
        n2s_all_nodes_idx )
             for (Yπ_net_nth_row,
                  nth_node_idx,
                  nth_node_adjacent_nodes_idx) in
                 zip(Yπ_net,
                     nodes_node_idx,
                     adjacent_nodes_idx )]
    
end

# ------------------------------------------
# ------------------------------------------

function get_a_node_ph_injection_by_Yπ_net(
    v_h,
    θ_h,
    v_k,
    θ_k;
    y_π)

    i_h, i_k =
        y_π * [v_h * exp(im * θ_h), v_k * exp(
            im * θ_k ) ]

    v_h * exp(im * θ_h) * conj(i_h)
    
    return  real(v_h * exp(im * θ_h) * conj(i_h))
end


function get_a_node_pk_injection_by_Yπ_net(
    v_h,
    θ_h,
    v_k,
    θ_k;
    y_π)

    i_h, i_k =
        y_π * [v_h * exp(im * θ_h), v_k * exp(
            im * θ_k ) ]
        
    return  real(v_k * exp(im * θ_k) * conj(i_k))
end


function get_a_node_qh_injection_by_Yπ_net(
    v_h,
    θ_h,
    v_k,
    θ_k;
    y_π)

    i_h, i_k =
        y_π * [v_h * exp(im * θ_h), v_k * exp(
            im * θ_k ) ]

    v_h * exp(im * θ_h) * conj(i_h)
    
    return  imag(v_h * exp(im * θ_h) * conj(i_h))
end


function get_a_node_qk_injection_by_Yπ_net(
    v_h,
    θ_h,
    v_k,
    θ_k;
    y_π)

    i_h, i_k =
        y_π * [v_h * exp(im * θ_h), v_k * exp(
            im * θ_k ) ]
        
    return  imag(v_k * exp(im * θ_k) * conj(i_k))
end




function get_an_edge_ph_injection(
    v_h,
    θ_h,
    v_k,
    θ_k,
    y_π)

    return (real(y_π[1,1]) * v_h^2 +
        v_h * v_k * (real(y_π[1,2]) * cos(θ_h - θ_k) +
        imag(y_π[1,2]) * sin(θ_h - θ_k))) - (
            real(y_π[2,2]) * v_k^2 +
                v_h * v_k * (real(y_π[1,2]) * cos(
                    θ_h - θ_k) -
                imag(y_π[1,2]) * sin(θ_h - θ_k)))

    
    
end


#

function get_an_edge_ph_injection(
    v_h,
    θ_h,
    v_k,
    θ_k;
    y_π)

    return real(y_π[1,1]) * v_h^2 +
            v_h * v_k * (real(y_π[1,2]) * cos(θ_h - θ_k) +
            imag(y_π[1,2]) * sin(θ_h - θ_k))

end


function get_an_edge_pk_injection(
    v_h,
    θ_h,
    v_k,
    θ_k,
    y_π)

    return real(y_π[2,2]) * v_k^2 +
            v_h * v_k * (real(y_π[1,2]) * cos(θ_h - θ_k) -
            imag(y_π[1,2]) * sin(θ_h - θ_k))
    
    
end


function get_an_edge_pk_injection(
    v_h,
    θ_h,
    v_k,
    θ_k;
    y_π)

    return  real(y_π[2,2]) * v_k^2 +
        v_h * v_k * (real(y_π[1,2]) * cos(
            θ_h - θ_k) - imag(y_π[1,2]) *
                sin(θ_h - θ_k))
    

end


function get_an_edge_qh_injection(
    v_h,
    θ_h,
    v_k,
    θ_k,
    y_π)

    return (-imag(y_π[1,1]) * v_h^2 +
        v_h * v_k * (real(y_π[1,2]) * sin(θ_h - θ_k) -
        imag(y_π[1,2]) * cos(θ_h - θ_k))) - (
            -imag(y_π[2,2]) * v_k^2 -
                v_h * v_k * (real(y_π[1,2]) *
                sin(θ_h - θ_k) +
                imag(y_π[1,2]) * cos(θ_h - θ_k)))
        
end



function get_an_edge_qh_injection(
    v_h,
    θ_h,
    v_k,
    θ_k;
    y_π)

    return -imag(y_π[1,1]) * v_h^2 +
        v_h * v_k * (real(y_π[1,2]) * sin(
            θ_h - θ_k) - imag(y_π[1,2]) *
                cos(θ_h - θ_k))
    


end


function get_an_edge_qk_injection(
    v_h,
    θ_h,
    v_k,
    θ_k,
    y_π)

    return  -imag(y_π[2,2]) * v_k^2 -
        v_h * v_k * (real(y_π[1,2]) * sin(
            θ_h - θ_k) +
            imag(y_π[1,2]) * cos(θ_h - θ_k))
        
end



function get_an_edge_qk_injection(
    v_h,
    θ_h,
    v_k,
    θ_k;
    y_π)
    
    return  -imag(y_π[2,2]) * v_k^2 -
        v_h * v_k * (real(y_π[1,2]) * sin(
            θ_h - θ_k) + imag(y_π[1,2]) *
                cos(θ_h - θ_k))        
end


# ------------------------------------------
# ------------------------------------------


function get_a_node_sh_injection_by_Yπ_net(
    v_h,
    θ_h,
    v_k,
    θ_k;
    Yπ)

    i_h, i_k =
        yπ * [v_h * exp(im * θ_h), v_k * exp(
            im * θ_k ) ]

    v_h * exp(im * θ_h) * conj(i_h)
    
    return  v_h * exp(im * θ_h) * conj(i_h)
end


function get_a_node_sk_injection_by_Yπ_net(
    v_h,
    θ_h,
    v_k,
    θ_k;
    Yπ)

    i_h, i_k =
        yπ * [v_h * exp(im * θ_h), v_k * exp(
            im * θ_k ) ]
        
    return  v_k * exp(im * θ_k) * conj(i_k)
end


# ------------------------------------------
# ------------------------------------------


# function get_nodes_∑_S_injection_by_Yπ_net(
#     vh,
#     θh;
#     Yπ_net,
#     Yshunt,
#     nodes_idx_with_adjacent_nodes_idx,
#     n2s_all_nodes_idx )

#     nodes_node_idx =
#         first.(nodes_idx_with_adjacent_nodes_idx)
    
#     adjacent_nodes_idx =
#         get_only_adjacent_nodes_idx(
#             nodes_idx_with_adjacent_nodes_idx)
    
#     return [ get_a_node_∑_S_injection_by_Yπ_net(
#         vh,
#         θh;
#         Yπ_net_nth_row,
#         Yshunt,
#         nth_node_idx,
#         nth_node_adjacent_nodes_idx,
#         n2s_all_nodes_idx )
#              for (Yπ_net_nth_row,
#                   nth_node_idx,
#                   nth_node_adjacent_nodes_idx) in
#                  zip(Yπ_net,
#                      nodes_node_idx,
#                      adjacent_nodes_idx )]
    
# end


#--------------------------------------
#--------------------------------------

function get_a_node_neighbourhood_edges_ih_ik_by_Yπ_net(
    vh,
    θh;
    Yπ_net_nth_row,
    Yshunt,
    nth_node_idx,
    nth_node_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in nth_node_adjacent_nodes_idx]
    
    vh_nth_node =
        vh[ transformed_nth_node_idx ]
    
    θh_nth_node =
        θh[ transformed_nth_node_idx ]

    uh_nth_node =
        vh_nth_node * exp(im * θh_nth_node )

    I_shunt_current_nth_node =
        Yshunt[transformed_nth_node_idx] *
        uh_nth_node

    vk_adj_nodes =
        vh[ transformed_adj_nodes_idx ]
    
    θk_adj_nodes =
        θh[ transformed_adj_nodes_idx ]

    uk_adj_nodes = vk_adj_nodes .* exp.(
        im * θk_adj_nodes)

    ih_ik = [ (yπ * [uh_nth_node, uk_adj_node ])
               for (yπ, uk_adj_node) in
                   zip(Yπ_net_nth_row,
                       uk_adj_nodes)]
    ih = first.( ih_ik )

    ik = last.( ih_ik )

    return (; ih, ik)
end


function get_a_node_neighbourhood_edges_uh_uk_by_Yπ_net(
    vh,
    θh;
    nth_node_idx,
    nth_node_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in nth_node_adjacent_nodes_idx]
    
    vh_nth_node =
        vh[ transformed_nth_node_idx ]
    
    θh_nth_node =
        θh[ transformed_nth_node_idx ]

    uh_nth_node =
        vh_nth_node * exp(im * θh_nth_node )

    vk_adj_nodes =
        vh[ transformed_adj_nodes_idx ]
    
    θk_adj_nodes =
        θh[ transformed_adj_nodes_idx ]

    uk_adj_nodes = vk_adj_nodes .* exp.(
        im * θk_adj_nodes)
    
    uh_uk = [[uh_nth_node, uk_adj_node ]
             for uk_adj_node in
                 uk_adj_nodes ]

    uh = first.(uh_uk)
    
    uk = last.(uh_uk)

    return (; uh, uk)
end


function get_a_node_neighbourhood_edges_vh_θh_vk_θh_Yπ_net(
    vh,
    θh;
    nth_node_idx,
    nth_node_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in nth_node_adjacent_nodes_idx]
    
    vh_nth_node =
        vh[ transformed_nth_node_idx ]
    
    θh_nth_node =
        θh[ transformed_nth_node_idx ]

    vk_adj_nodes =
        vh[ transformed_adj_nodes_idx ]
    
    θk_adj_nodes =
        θh[ transformed_adj_nodes_idx ]
    
    return [(vh_nth_node,
              θh_nth_node,
              vk_adj_node,
              θk_adj_node)
             for (vk_adj_node, θk_adj_node) in
                 zip(vk_adj_nodes,
                 θk_adj_nodes) ]

end

# 


function get_a_node_neighbourhood_edges_ph_pk_by_Yπ_net(
    vh,
    θh;
    Yπ_net_nth_row,
    Yshunt,
    nth_node_idx,
    nth_node_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in nth_node_adjacent_nodes_idx]
    
    vh_nth_node =
        vh[ transformed_nth_node_idx ]
    
    θh_nth_node =
        θh[ transformed_nth_node_idx ]

    uh_nth_node =
        vh_nth_node * exp(im * θh_nth_node )

    vk_adj_nodes =
        vh[ transformed_adj_nodes_idx ]
    
    θk_adj_nodes =
        θh[ transformed_adj_nodes_idx ]

    ph_pk = [(get_an_edge_ph_injection(
        vh_nth_node, θh_nth_node, vk, θk; y_π),
      get_an_edge_pk_injection(
          vh_nth_node, θh_nth_node, vk, θk; y_π))
     for (vk, θk, y_π) in
         zip(vk_adj_nodes,
             θk_adj_nodes,
             Yπ_net_nth_row) ]
    
    ph = first.(ph_pk)

    pk = last.(ph_pk)

    return (; ph, pk)
    
end


function get_a_node_neighbourhood_edges_qh_qk_by_Yπ_net( vh, θh;
    Yπ_net_nth_row,
    Yshunt,
    nth_node_idx,
    nth_node_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in nth_node_adjacent_nodes_idx]
    
    vh_nth_node =
        vh[ transformed_nth_node_idx ]
    
    θh_nth_node =
        θh[ transformed_nth_node_idx ]

    uh_nth_node =
        vh_nth_node * exp(im * θh_nth_node )

    vk_adj_nodes =
        vh[ transformed_adj_nodes_idx ]
    
    θk_adj_nodes =
        θh[ transformed_adj_nodes_idx ]

    qh_qk = [(get_an_edge_qh_injection(
        vh_nth_node, θh_nth_node, vk, θk; y_π),
      get_an_edge_qk_injection(
          vh_nth_node, θh_nth_node, vk, θk; y_π))
     for (vk, θk, y_π) in
         zip(vk_adj_nodes,
             θk_adj_nodes,
             Yπ_net_nth_row) ]
    
    qh = first.(qh_qk)

    qk = last.(qh_qk)

    return (; qh, qk)
    
end


# -------------------------------------------
# node net current injection based on Ybus
# -------------------------------------------

"""
Columns are used because Julia access matrix via cols. It does not make any diffence, because Ybus is symetric.

`Ybus_nth_col = Ybus[:,nth] == Ybus_nth_row = Ybus[nth,:]`

`Ybus_nth_col' * (vh .* exp.( im * θh ))`

get_a_node_∑_ynj_x_vj_by_Ybus(
    vh,
    θh,
    Ybus[:,1] )

"""
function get_a_node_∑_ynj_x_vj_by_Ybus(
    vh,
    θh,
    Ybus_nth_col )
    
    return  sum(Ybus_nth_col .* (vh .* exp.( im * θh )))
end

function get_a_node_∑_ynj_x_vj_by_Ybus(
    uh,
    Ybus_nth_col )
    
    return  sum(Ybus_nth_col .* uh)
end


function get_nodes_∑_ynj_x_vj_by_Ybus_per_col(
    vh,
    θh,
    Ybus)

    # no_rows could also be used
    # no_row = size(Ybus)[1]

    no_cols = size(Ybus)[2]
    
    
    return [ get_a_node_∑_ynj_x_vj_by_Ybus(
        vh, θh, Ybus[:,a_col] )
             for a_col in 1:no_cols ]
    
end


function get_nodes_∑_ynj_x_vj_by_Ybus_per_col(
    uh,
    Ybus)

    # no_rows could also be used
    # no_row = size(Ybus)[1]

    no_cols = size(Ybus)[2]
    
    
    return [ get_a_node_∑_ynj_x_vj_by_Ybus(
        uh, Ybus[:,a_col] )
             for a_col in 1:no_cols ]
    
end


function get_nodes_∑_ynj_x_vj_by_Ybus(
    vh,
    θh,
    Ybus)
    
    return  Ybus * (vh .* exp.(im * θh))
end


function get_nodes_∑_ynj_x_vj_by_Ybus(
    vh,
    θh;
    Ybus)
    
    return  Matrix(Ybus) * (vh .* exp.(im * θh))
end


# function get_nodes_∑_ynj_x_vj_by_Ybus(
#     uh,
#     Ybus)
    
#     return  Ybus * uh
# end

# -------------------------------------------
# node net current injection based on Yπ_net
# -------------------------------------------

function get_a_node_∑_ynj_x_vj_by_Yπ_net(
    uh,
    Yπ_net_nth_row,
    Yshunt_nth_node,
    nth_node_idx_and_adj_nodes_idx,
    n2s_all_nodes_idx)

    nth_node_idx  =
        nth_node_idx_and_adj_nodes_idx[1]
        
    adj_nodes_idx =
        nth_node_idx_and_adj_nodes_idx[2:end]
    

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [ n2s_all_nodes_idx[idx] for idx in
              adj_nodes_idx ]

    uh_nth_node =
        uh[ transformed_nth_node_idx]

    I_shunt_current_nth_node =
        Yshunt_nth_node * uh_nth_node

    uk_adj_nodes = uh[transformed_adj_nodes_idx]

    return sum( [ (yπ * [uh_nth_node, uk_adj_node ])[1]
               for (yπ, uk_adj_node) in
                   zip(Yπ_net_nth_row,
                       uk_adj_nodes)]) +
                           I_shunt_current_nth_node

end



function get_nodes_∑_ynj_x_vj_by_Yπ_net(
    uh,
    Yπ_net,
    Yshunt,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    # nodes_idx_with_adjacent_nodes_idx,

    return [get_a_node_∑_ynj_x_vj_by_Yπ_net(
        uh,
        Yπ_net_nth_row,
        Yshunt_nth_node,
        nth_node_idx_and_adj_nodes_idx,
        n2s_all_nodes_idx)
             for (Yπ_net_nth_row,
                  Yshunt_nth_node,
                  nth_node_idx_and_adj_nodes_idx ) in
                 zip(Yπ_net,
                     Yshunt,
                     nodes_idx_with_adjacent_nodes_idx )  ]
end



function get_a_node_∑_ynj_x_vj_by_Yπ_net(
    uh;
    Yπ_net_nth_row,
    Yshunt_nth_node,
    nth_node_idx_and_adj_nodes_idx,
    n2s_all_nodes_idx)

    nth_node_idx  =
        nth_node_idx_and_adj_nodes_idx[1]
        
    adj_nodes_idx =
        nth_node_idx_and_adj_nodes_idx[2:end]

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]
    
    transformed_adj_nodes_idx =
        [ n2s_all_nodes_idx[idx] for idx in adj_nodes_idx]
    
    uh_nth_node =
        uh[transformed_nth_node_idx]

    I_shunt_current_nth_node =
        Yshunt_nth_node * uh_nth_node

    uk_adj_nodes = uh[transformed_adj_nodes_idx]

    return sum( [ (yπ * [uh_nth_node,uk_adj_node ])[1]
               for (yπ, uk_adj_node) in
                   zip(Yπ_net_nth_row,
                       uk_adj_nodes)]) +
                           I_shunt_current_nth_node

end


function get_nodes_∑_ynj_x_vj_by_Yπ_net(
    uh;
    Yπ_net,
    Yshunt,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    # nodes_idx_with_adjacent_nodes_idx,

    return [get_a_node_∑_ynj_x_vj_by_Yπ_net(
        uh;
        Yπ_net_nth_row,
        Yshunt_nth_node,
        nth_node_idx_and_adj_nodes_idx,
        n2s_all_nodes_idx)
             for (Yπ_net_nth_row,
                  Yshunt_nth_node,
                  nth_node_idx_and_adj_nodes_idx ) in
                 zip(Yπ_net,
                     Yshunt,
                     nodes_idx_with_adjacent_nodes_idx )  ]
end


function get_a_node_∑_ynj_x_vj_by_Yπ_net(
    vh,
    θh,
    Yπ_net_nth_row,
    Yshunt,
    nth_node_idx,
    nth_node_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in nth_node_adjacent_nodes_idx]
    
    vh_nth_node =
        vh[ transformed_nth_node_idx ]
    
    θh_nth_node =
        θh[ transformed_nth_node_idx ]

    uh_nth_node =
        vh_nth_node * exp(im * θh_nth_node )

    I_shunt_current_nth_node =
        Yshunt[transformed_nth_node_idx] * uh_nth_node

    vk_adj_nodes =
        vh[ transformed_adj_nodes_idx ]
    
    θk_adj_nodes =
        θh[ transformed_adj_nodes_idx ]

    uk_adj_nodes = vk_adj_nodes .* exp.(im * θk_adj_nodes)
    
    return  sum( [ (yπ * [uh_nth_node,uk_adj_node ])[1]
               for (yπ, uk_adj_node) in
                   zip(Yπ_net_nth_row,
                       uk_adj_nodes)]) +
                           I_shunt_current_nth_node
end


function get_nodes_∑_ynj_x_vj_by_Yπ_net(
    vh,
    θh,
    Yπ_net,
    Yshunt,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx )

    nodes_node_idx =
        first.(nodes_idx_with_adjacent_nodes_idx)
    
    adjacent_nodes_idx =
        get_only_adjacent_nodes_idx(
            nodes_idx_with_adjacent_nodes_idx)
    
    return [ get_a_node_∑_ynj_x_vj_by_Yπ_net(
        vh,
        θh,
        Yπ_net_nth_row,
        Yshunt,
        nth_node_idx,
        nth_node_adjacent_nodes_idx,
        n2s_all_nodes_idx )
             for (Yπ_net_nth_row,
                  nth_node_idx,
                  nth_node_adjacent_nodes_idx) in
                 zip(Yπ_net,
                     nodes_node_idx,
                     adjacent_nodes_idx )]
    
end


function get_a_node_∑_ynj_x_vj_by_Yπ_net(
    vh,
    θh;
    Yπ_net_nth_row,
    Yshunt,
    nth_node_idx,
    nth_node_adjacent_nodes_idx,
    n2s_all_nodes_idx)

    transformed_nth_node_idx =
        n2s_all_nodes_idx[nth_node_idx]

    transformed_adj_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in nth_node_adjacent_nodes_idx]
    
    vh_nth_node =
        vh[ transformed_nth_node_idx ]
    
    θh_nth_node =
        θh[ transformed_nth_node_idx ]

    uh_nth_node =
        vh_nth_node * exp(im * θh_nth_node )

    I_shunt_current_nth_node =
        Yshunt[transformed_nth_node_idx] * uh_nth_node

    vk_adj_nodes =
        vh[ transformed_adj_nodes_idx ]
    
    θk_adj_nodes =
        θh[ transformed_adj_nodes_idx ]

    uk_adj_nodes = vk_adj_nodes .* exp.(im * θk_adj_nodes)
    
    return  sum( [ (yπ * [uh_nth_node,uk_adj_node ])[1]
               for (yπ, uk_adj_node) in
                   zip(Yπ_net_nth_row,
                       uk_adj_nodes)]) +
                           I_shunt_current_nth_node
end


function get_nodes_∑_ynj_x_vj_by_Yπ_net(
    vh,
    θh;
    Yπ_net,
    Yshunt,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx )

    nodes_node_idx =
        first.(nodes_idx_with_adjacent_nodes_idx)
    
    adjacent_nodes_idx =
        get_only_adjacent_nodes_idx(
            nodes_idx_with_adjacent_nodes_idx)
    
    return [ get_a_node_∑_ynj_x_vj_by_Yπ_net(
        vh,
        θh;
        Yπ_net_nth_row,
        Yshunt,
        nth_node_idx,
        nth_node_adjacent_nodes_idx,
        n2s_all_nodes_idx )
             for (Yπ_net_nth_row,
                  nth_node_idx,
                  nth_node_adjacent_nodes_idx) in
                 zip(Yπ_net,
                     nodes_node_idx,
                     adjacent_nodes_idx )]
    
end

# -------------------------------------------
# node net current injection based on Ynet
# -------------------------------------------


function get_a_node_∑_ynj_x_vj(
    vh,
    θh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    n2s_all_nodes_idx )
    
    return  sum( [ ynj * vh[ n2s_all_nodes_idx[node_idx] ] *
        exp( im * θh[ n2s_all_nodes_idx[node_idx]] )
               for (ynj, node_idx) in
                   zip(Ynet_nth_row,
                       nth_node_idx_and_adj_nodes_idx)])

end


function get_a_node_∑_ynj_x_vj(
    uh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx )

    # node_idx = first(nth_node_idx_and_adj_nodes_idx )

    return sum([ynj * uh[node_idx]
                for (ynj, node_idx) in
                    zip(Ynet_nth_row,
                        nth_node_idx_and_adj_nodes_idx)])
end


function get_nodes_∑_ynj_x_vj(
    vh,
    θh,
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx,
    n2s_all_nodes_idx )

    return [ get_a_node_∑_ynj_x_vj(
        vh,
        θh,
        Ynet_nth_row,
        nth_node_idx_and_adj_nodes_idx,
        n2s_all_nodes_idx )
             for (Ynet_nth_row,
                  nth_node_idx_and_adj_nodes_idx) in
                 zip(
                     Ynet,
                     nodes_node_idx_and_incident_edges_other_node_idx)]
    
end


function get_nodes_∑_ynj_x_vj(
    uh,
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx )

    return [ get_a_node_∑_ynj_x_vj(
        uh, Ynet_nth_row, nth_node_idx_and_adj_nodes_idx )
             for (Ynet_nth_row, nth_node_idx_and_adj_nodes_idx ) in
                 zip(Ynet,
                     nodes_node_idx_and_incident_edges_other_node_idx )  ]
end



# -------------------------------------------
# node net current injection with B_nj
# -------------------------------------------


function get_a_node_∑_Bnj_x_vj(
    vh,
    θh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    n2s_all_nodes_idx )
    
    return  sum( [ im * imag(ynj) * vh[ n2s_all_nodes_idx[node_idx] ] *
        exp( im * θh[ n2s_all_nodes_idx[node_idx]] )
               for (ynj, node_idx) in
                   zip(Ynet_nth_row,
                       nth_node_idx_and_adj_nodes_idx)])

end


function get_a_node_∑_Bnj_x_vj(
    uh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx )

    # node_idx = first(nth_node_idx_and_adj_nodes_idx )

    return sum([ im * imag(ynj) * uh[node_idx]
                for (ynj, node_idx) in
                    zip(Ynet_nth_row,
                        nth_node_idx_and_adj_nodes_idx)])
end


function get_nodes_∑_Bnj_x_vj(
    vh,
    θh,
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx,
    n2s_all_nodes_idx )

    return [ get_a_node_∑_Bnj_x_vj(
        vh,
        θh,
        Ynet_nth_row,
        nth_node_idx_and_adj_nodes_idx,
        n2s_all_nodes_idx )
             for (Ynet_nth_row,
                  nth_node_idx_and_adj_nodes_idx) in
                 zip(
                     Ynet,
                     nodes_node_idx_and_incident_edges_other_node_idx)]
    
end


function get_nodes_∑_Bnj_x_vj(
    uh,
    Ynet,
    nodes_node_idx_and_incident_edges_other_node_idx )

    return [ get_a_node_∑_Bnj_x_vj(
        uh, Ynet_nth_row, nth_node_idx_and_adj_nodes_idx )
             for (Ynet_nth_row, nth_node_idx_and_adj_nodes_idx ) in
                 zip(Ynet,
                     nodes_node_idx_and_incident_edges_other_node_idx )  ]
end


# -------------------------------------------
#  power injection and current into gens nodes 
# -------------------------------------------

function get_gens_nodes_Pg_net_inj(
    gens_δ,
    gens_vh,
    gens_θh,
    gens_i_d,
    gens_i_q,
    P_g_loc_load,
    Q_g_loc_load,

    n2s_gens_idx,
    gens_nodes_with_loc_loads_idx,

    gens_nodes_idx)

    
    gens_Pg = get_gens_Pg_from_δ_i_dq_vh_θh(
        gens_vh, gens_θh,
        gens_δ,
        gens_i_d, gens_i_q )

    return [idx ∈ gens_nodes_with_loc_loads_idx ?
        gens_Pg[ n2s_gens_idx[ idx]] -
        P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx] ] :
        gens_Pg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]

end


function get_gens_nodes_Qg_net_inj(
    gens_δ,
    gens_vh,
    gens_θh,
    gens_i_d,
    gens_i_q,
    P_g_loc_load,
    Q_g_loc_load,

    n2s_gens_idx,
    gens_nodes_with_loc_loads_idx,

    gens_nodes_idx)

    gens_Qg = get_gens_Qg_from_δ_i_dq_vh_θh(
        gens_vh, gens_θh,
        gens_δ,
        gens_i_d, gens_i_q )


    return [idx ∈ gens_nodes_with_loc_loads_idx ?
        gens_Qg[ n2s_gens_idx[ idx]] -
        Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx] ] :
        gens_Qg[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]
    
    
end


function get_gens_nodes_I_D_net_inj(
    gens_δ,
    gens_vh,
    gens_θh,
    gens_i_d,
    gens_i_q,
    P_g_loc_load,
    Q_g_loc_load,

    n2s_gens_idx,
    gens_nodes_with_loc_loads_idx,

    gens_nodes_idx)

    return [idx ∈ gens_nodes_with_loc_loads_idx ?
        (gens_i_q[ n2s_gens_idx[idx]] *
        cos( gens_δ[ n2s_gens_idx[idx]]) +
        gens_i_d[ n2s_gens_idx[idx]] *
        sin( gens_δ[ n2s_gens_idx[idx]]) -
        (P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] *
        cos(gens_θh[ n2s_gens_idx[ idx]]) +
        Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] *
        sin(gens_θh[ n2s_gens_idx[ idx]]))/
        gens_vh[ n2s_gens_idx[idx]]) :
        (gens_i_q[ n2s_gens_idx[idx]] *
        cos( gens_δ[ n2s_gens_idx[idx]]) +
        gens_i_d[ n2s_gens_idx[idx]] *
        sin( gens_δ[ n2s_gens_idx[idx]])) 
            for idx in gens_nodes_idx ]
    
end


function get_gens_nodes_I_Q_net_inj(
    gens_δ,
    gens_vh,
    gens_θh,
    gens_i_d,
    gens_i_q,
    P_g_loc_load,
    Q_g_loc_load,

    n2s_gens_idx,
    gens_nodes_with_loc_loads_idx,

    gens_nodes_idx)

    return [idx ∈ gens_nodes_with_loc_loads_idx ?
        (gens_i_q[ n2s_gens_idx[idx]] *
        sin( gens_δ[ n2s_gens_idx[idx]]) -
        gens_i_d[ n2s_gens_idx[idx]] *
        cos( gens_δ[ n2s_gens_idx[idx]]) -
        (P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] *
        sin(gens_θh[ n2s_gens_idx[ idx]]) -
        Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] *
        cos(gens_θh[ n2s_gens_idx[ idx]]))/
        gens_vh[ n2s_gens_idx[idx]]) :
        (gens_i_q[ n2s_gens_idx[idx]] *
        sin( gens_δ[ n2s_gens_idx[idx]]) -
        gens_i_d[ n2s_gens_idx[idx]] *
        cos( gens_δ[ n2s_gens_idx[idx]])) 
                  for idx in gens_nodes_idx ]
end

# -------------------------------------------
#  power injection from net to a node 
# -------------------------------------------

function get_nth_node_real_uh_x_∑_ynj_x_vj(
    vh,
    θh,           
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_idx,
    n2s_all_nodes_idx)
        
    return  vh[ n2s_all_nodes_idx[ nth_idx]] * sum(
        [ vh[ n2s_all_nodes_idx[ node_idx]] * abs(ynj) * cos(
            θh[ n2s_all_nodes_idx[ nth_idx]] -
                θh[ n2s_all_nodes_idx[ node_idx]] - angle( ynj ))
             for (ynj, node_idx) in
                 zip(
                     Ynet_nth_row ,
                     nth_node_idx_and_adj_nodes_idx)])
end


function get_nth_node_imag_uh_x_∑_ynj_x_vj(
    vh,
    θh,           
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_idx,
    n2s_all_nodes_idx )
        
    return  vh[ n2s_all_nodes_idx[ nth_idx]] * sum(
        [ vh[ n2s_all_nodes_idx[ node_idx]] * abs(ynj) * sin(
            θh[ n2s_all_nodes_idx[ nth_idx]] -
                θh[ n2s_all_nodes_idx[ node_idx]] - angle(ynj))
          for (ynj, node_idx) in
              zip(Ynet_nth_row,
                  nth_node_idx_and_adj_nodes_idx)])
end


function get_nth_node_real_uh_x_∑_ynj_x_vj(
    vh,
    θh,           
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_idx )
        
    return  vh[nth_idx] * sum(
        [ vh[node_idx] * abs(ynj) * cos(
            θh[nth_idx] - θh[node_idx] - angle( ynj ))
             for (ynj, node_idx) in
                 zip(
                     Ynet_nth_row ,
                     nth_node_idx_and_adj_nodes_idx)])
end


function get_nth_node_imag_uh_x_∑_ynj_x_vj(
    vh,
    θh,           
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_idx )
        
    return  vh[nth_idx] * sum(
        [ vh[node_idx] * abs(ynj) * sin(
            θh[nth_idx] - θh[node_idx] - angle(ynj))
          for (ynj, node_idx) in
              zip(Ynet_nth_row,
                  nth_node_idx_and_adj_nodes_idx)])
end


function get_nodes_real_uh_x_∑_ynj_x_vj(
    vh,
    θh,
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    all_nodes_idx,
    n2s_all_nodes_idx)

    return [
        get_nth_node_real_uh_x_∑_ynj_x_vj(
            vh,
            θh,
            Ynet[ n2s_all_nodes_idx[nth_idx] ],
            nodes_idx_with_adjacent_nodes_idx[
                n2s_all_nodes_idx[nth_idx] ],
            nth_idx,
            n2s_all_nodes_idx)
             for nth_idx in all_nodes_idx ]
    

end


function get_nodes_imag_uh_x_∑_ynj_x_vj(
    vh,
    θh,
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    all_nodes_idx,
    n2s_all_nodes_idx)


    return [
        get_nth_node_imag_uh_x_∑_ynj_x_vj(
            vh,
            θh,
            Ynet[ n2s_all_nodes_idx[nth_idx] ],
            nodes_idx_with_adjacent_nodes_idx[
                n2s_all_nodes_idx[nth_idx] ],
            nth_idx,
            n2s_all_nodes_idx)
             for nth_idx in all_nodes_idx ]
    

end


#----------------------------------------
# non consequtive indices
#----------------------------------------

function get_nth_node_real_uh_x_∑_ynj_x_vj(
    vh,
    θh;
    n2s_all_nodes_idx,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_idx )
        
    return  vh[ n2s_all_nodes_idx[ nth_idx ]] *
        sum( [ vh[ n2s_all_nodes_idx[ node_idx ] ] *
        abs(ynj) * cos( θh[ n2s_all_nodes_idx[ nth_idx ]] -
        θh[n2s_all_nodes_idx[ node_idx] ] - angle( ynj ))
             for (ynj, node_idx) in
                 zip(
                     Ynet_nth_row ,
                     nth_node_idx_and_adj_nodes_idx)])
end


function get_nth_node_imag_uh_x_∑_ynj_x_vj(
    vh,
    θh;
    n2s_all_nodes_idx,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_idx )
        
    return  vh[ n2s_all_nodes_idx[ nth_idx ]] *
        sum([ vh[ n2s_all_nodes_idx[ node_idx ]] *
        abs(ynj) * sin( θh[ n2s_all_nodes_idx[ nth_idx ]] -
                    θh[ n2s_all_nodes_idx[ node_idx ]] -
                    angle(ynj))
          for (ynj, node_idx) in
              zip(Ynet_nth_row,
                  nth_node_idx_and_adj_nodes_idx)])
end


function get_nodes_real_uh_x_∑_ynj_x_vj(
    vh,
    θh;
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    all_nodes_idx,
    n2s_all_nodes_idx)

    return [
        get_nth_node_real_uh_x_∑_ynj_x_vj(
            vh,
            θh;
            n2s_all_nodes_idx,
            Ynet_nth_row =
                Ynet[ n2s_all_nodes_idx[nth_idx] ],
            nth_node_idx_and_adj_nodes_idx =
                nodes_idx_with_adjacent_nodes_idx[
                    n2s_all_nodes_idx[nth_idx] ],
            nth_idx )
             for nth_idx in all_nodes_idx ]
    

end


function get_nodes_imag_uh_x_∑_ynj_x_vj(
    vh,
    θh;
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    all_nodes_idx,
    n2s_all_nodes_idx)


    return [
        get_nth_node_imag_uh_x_∑_ynj_x_vj(
            vh,
            θh;
            n2s_all_nodes_idx,
            Ynet_nth_row =
                Ynet[ n2s_all_nodes_idx[nth_idx] ],
            nth_node_idx_and_adj_nodes_idx =
                nodes_idx_with_adjacent_nodes_idx[
                    n2s_all_nodes_idx[nth_idx] ],
            nth_idx )
             for nth_idx in all_nodes_idx ]
    

end



#----------------------------------------
#----------------------------------------


# function get_nth_node_real_uh_x_∑_ynj_x_vj(
#     vh,
#     θh;           
#     Ynet_nth_row,
#     nth_node_idx_and_adj_nodes_idx,
#     nth_idx )
        
#     return  vh[nth_idx] * sum(
#         [ vh[node_idx] * abs(ynj) * cos(
#             θh[nth_idx] - θh[node_idx] - angle( ynj ))
#              for (ynj, node_idx) in
#                  zip(
#                      Ynet_nth_row ,
#                      nth_node_idx_and_adj_nodes_idx)])
# end


# function get_nth_node_imag_uh_x_∑_ynj_x_vj(
#     vh,
#     θh;           
#     Ynet_nth_row,
#     nth_node_idx_and_adj_nodes_idx,
#     nth_idx )
        
#     return  vh[nth_idx] * sum(
#         [ vh[node_idx] * abs(ynj) * sin(
#             θh[nth_idx] - θh[node_idx] - angle(ynj))
#           for (ynj, node_idx) in
#               zip(Ynet_nth_row,
#                   nth_node_idx_and_adj_nodes_idx)])
# end


# function get_nodes_real_uh_x_∑_ynj_x_vj(
#     vh,
#     θh;
#     Ynet,
#     nodes_idx_with_adjacent_nodes_idx,
#     all_nodes_idx,
#     n2s_all_nodes_idx)

#     return [
#         get_nth_node_real_uh_x_∑_ynj_x_vj(
#             vh,
#             θh;
#             Ynet_nth_row =
#                 Ynet[ n2s_all_nodes_idx[nth_idx] ],
#             nth_node_idx_and_adj_nodes_idx =
#                 nodes_idx_with_adjacent_nodes_idx[
#                     n2s_all_nodes_idx[nth_idx] ],
#             nth_idx )
#              for nth_idx in all_nodes_idx ]
    

# end


# function get_nodes_imag_uh_x_∑_ynj_x_vj(
#     vh,
#     θh;
#     Ynet,
#     nodes_idx_with_adjacent_nodes_idx,
#     all_nodes_idx,
#     n2s_all_nodes_idx)


#     return [
#         get_nth_node_imag_uh_x_∑_ynj_x_vj(
#             vh,
#             θh;
#             Ynet_nth_row =
#                 Ynet[ n2s_all_nodes_idx[nth_idx] ],
#             nth_node_idx_and_adj_nodes_idx =
#                 nodes_idx_with_adjacent_nodes_idx[
#                     n2s_all_nodes_idx[nth_idx] ],
#             nth_idx )
#              for nth_idx in all_nodes_idx ]
    

# end


#----------------------------------------
#----------------------------------------


function get_nth_node_real_uh_x_∑_Bnj_x_vj(
    vh,
    θh,           
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_idx )
        
    return  vh[nth_idx] * sum(
        [ vh[node_idx] * imag(ynj) * sin(
            θh[nth_idx] - θh[node_idx] )
             for (ynj, node_idx) in
                 zip(
                     Ynet_nth_row ,
                     nth_node_idx_and_adj_nodes_idx)])
end


function get_nth_node_imag_uh_x_∑_Bnj_x_vj(
    vh,
    θh,           
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_idx )
        
    return  vh[nth_idx] * sum(
        [ vh[node_idx] * imag(ynj) * cos(
            θh[nth_idx] - θh[node_idx])
          for (ynj, node_idx) in
              zip(Ynet_nth_row,
                  nth_node_idx_and_adj_nodes_idx)])
end


function get_nodes_real_uh_x_∑_Bnj_x_vj(
    vh,
    θh,
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    all_nodes_idx,
    n2s_all_nodes_idx)

    return [
        get_nth_node_real_uh_x_∑_Bnj_x_vj(
            vh,
            θh,
            Ynet[ n2s_all_nodes_idx[nth_idx] ],
            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx] ],
            nth_idx )
             for nth_idx in all_nodes_idx ]
    

end


function get_nodes_imag_uh_x_∑_Bnj_x_vj(
    vh,
    θh,
    Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    all_nodes_idx,
    n2s_all_nodes_idx)


    return [
        get_nth_node_imag_uh_x_∑_Bnj_x_vj(
            vh,
            θh,
            Ynet[ n2s_all_nodes_idx[nth_idx] ],
            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx] ],
            nth_idx )
             for nth_idx in all_nodes_idx ]
    

end


#----------------------------------------
#----------------------------------------


"""
The first element of
nth_node_idx_and_adj_nodes_idx is the idx
of the node
"""
function get_a_node_uh_x_conj_∑_ynj_x_vj(
    vh,
    θh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    n2s_all_nodes_idx)

    node_idx = first(nth_node_idx_and_adj_nodes_idx )
    
    return vh[ n2s_all_nodes_idx[node_idx] ] * exp(
        im * θh[ n2s_all_nodes_idx[node_idx] ] ) *
            conj( sum([ynj * vh[ n2s_all_nodes_idx[idx]] *
            exp(im * θh[ n2s_all_nodes_idx[idx]] )
              for (ynj, idx) in
                   zip(Ynet_nth_row,
                       nth_node_idx_and_adj_nodes_idx)]) )
end

function get_a_node_uh_x_conj_∑_ynj_x_vj_real(
    vh,
    θh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    n2s_all_nodes_idx)
    
    return real(
        get_a_node_uh_x_conj_∑_ynj_x_vj(
            vh,
            θh,
            Ynet_nth_row,
            nth_node_idx_and_adj_nodes_idx,
            n2s_all_nodes_idx))
    
end


function get_a_node_uh_x_conj_∑_ynj_x_vj_imag(
    vh,
    θh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    n2s_all_nodes_idx)
    
    return imag(
        get_a_node_uh_x_conj_∑_ynj_x_vj(
            vh, θh,
            Ynet_nth_row,
            nth_node_idx_and_adj_nodes_idx,
            n2s_all_nodes_idx))
    
end

#----------------------------------------

function get_a_node_uh_x_conj_∑_ynj_x_vj(
    uh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx )

    node_idx = first(nth_node_idx_and_adj_nodes_idx )
    
    return uh[ node_idx ] *
        conj( sum([ ynj * uh[idx] 
              for (ynj, idx) in
                   zip(Ynet_nth_row,
                       nth_node_idx_and_adj_nodes_idx)]))
end


function get_a_node_uh_x_conj_∑_ynj_x_vj_real(
    uh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx )
    
    return real(
        get_a_node_uh_x_conj_∑_ynj_x_vj(
            uh,
            Ynet_nth_row,
            nth_node_idx_and_adj_nodes_idx))
    
end


function get_a_node_uh_x_conj_∑_ynj_x_vj_imag(
    uh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx )
    
    return imag(
        get_a_node_uh_x_conj_∑_ynj_x_vj(
            uh,
            Ynet_nth_row,
            nth_node_idx_and_adj_nodes_idx ))
    
end


# -------------------------------------------
# A gen complex power injection
# -------------------------------------------


function get_a_gen_power_injection_by_vh_θh_id_iq(
    vh_θh_δ_id_iq )

    return  get_a_gen_power_injection_by_vh_θh_id_iq(
    vh_θh_δ_id_iq... )
end

function get_a_gen_power_injection_by_vh_θh_id_iq(
    vh_θh, δ, id_iq)
    vh_θh_id_iq
    return  get_a_gen_power_injection_by_vh_θh_id_iq(
    vh_θh..., δ, id_iq... )
end


function get_a_gen_power_injection_by_vh_θh_id_iq(
    vh_θh, δ, id, iq)
    
    return  get_a_gen_power_injection_by_vh_θh_id_iq(
    vh_θh..., δ, id, iq)
end

function get_a_gen_power_injection_by_vh_θh_id_iq(
    vh, θh, δ, id, iq)
    
    return vh * exp(im * θh ) *
        (id - im * iq) *
        exp(-im * ( δ - π/2 ))
end



# -------------------------------------------


function get_a_gen_active_power_injection_by_vh_θh_id_iq(
    vh, θh, δ, id, iq)
    
    return vh * iq * cos(θh -  δ) -
        vh * id * sin(θh -  δ)
    
end


function get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
    vh, θh, δ, id, iq)
    
    return vh * id * cos(θh -  δ) +
        vh * iq * sin(θh -  δ)
    
end

# -------------------------------------------

function get_a_gen_power_injection(
    vh_θh, δ_ω_ed_dash_eq_dash,
    ra_X_d_dash_X_q_dash )
    
    return get_a_gen_power_injection(
    vh_θh..., δ_ω_ed_dash_eq_dash...,
    ra_X_d_dash_X_q_dash... )
end

function get_a_gen_power_injection(
    vh, θh, δ_ω_ed_dash_eq_dash,
    ra_X_d_dash_X_q_dash )
    
    return get_a_gen_power_injection(
    vh, θh, δ_ω_ed_dash_eq_dash...,
    ra_X_d_dash_X_q_dash... )
end


function get_a_gen_power_injection(
    vh, θh, δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash )

    id_iq = get_a_gen_dyn_idq(
        vh, θh,
        δ, ω, ed_dash, eq_dash,
        ra, X_d_dash, X_q_dash )

    id, iq = id_iq
    
    return vh * exp(im * θh ) *
        (id - im * iq) * exp(-im * ( δ - π/2 ))
end


# -------------------------------------------
# Complex gen stator equation
# -------------------------------------------


function get_a_gen_complex_stator_equations(
    vh_θh, δ_ω_ed_eq, ra_X_d_dash_X_q_dash, id_iq )

    return  get_a_gen_complex_stator_equations(
        vh_θh...,
        δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash...,
        id_iq... )

end

function get_a_gen_complex_stator_equations(
    vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash, id_iq )

    return get_a_gen_complex_stator_equations(
        vh, θh, δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash...,
        id_iq... )

end


function get_a_gen_complex_stator_equations(
    vh, θh, δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash, id, iq )

    return vh * exp(im * θh ) +
        ( ra + im * X_d_dash ) *
        (id + im * iq) * exp( im * ( δ - pi/2)) -
        ( ed_dash + ( X_q_dash -
        X_d_dash) * iq + im * eq_dash  ) *
        exp( im * ( δ - pi/2))

end


# -------------------------------------------


function get_a_gen_stator_equations(
    vh_θh, δ_ω_ed_eq, ra_X_d_dash_X_q_dash, id_iq )

    return  get_a_gen_stator_equations(
        vh_θh..., δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash...,
        id_iq... )

end

function get_a_gen_stator_equations(
    vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash, id_iq )

    return get_a_gen_stator_equations(
        vh, θh, δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash...,
        id_iq... )

end


function get_a_gen_stator_equations(
    vh, θh, δ_ω_ed_eq, ra_Xd_dash_Xq_dash, id, iq )

    return get_a_gen_stator_equations(
        vh, θh, δ_ω_ed_eq...,
        ra_Xd_dash_Xq_dash...,
        id, iq )

end


function get_a_gen_stator_equations(
    vh, θh, δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash, id, iq )

    # [ed_dash - vh * sin(δ - θh) - ra * id + X_q_dash * iq,
    #  eq_dash - vh * cos(δ - θh) - X_d_dash * id - ra * iq ]

    return   ([ed_dash - vh * sin(δ - θh),
               eq_dash - vh * cos(δ - θh)]  -
                  Z_dq(
                      ra, X_d_dash,
                      X_q_dash) * [id, iq])

end


function get_flat_gens_stator_equations_and_Idxs(
    gens_stator_equations )


    """ gse : gens_stator_equations """
    
    gse_eq_1 = first.(
        gens_stator_equations )

    gse_eq_2 = second.(
        gens_stator_equations )
    
    gse_flat = [ gse_eq_1;
                 gse_eq_2 ]
    
    dims_gse =
        length.([ gse_eq_1,
                  gse_eq_2 ])

     _,_,  gse_eq_Idx =
         create_size_offset_Idx(
             dims_gse;
             counter = 0)

    gse_eq_1_Idx = gse_eq_Idx[ 1 ]
    
    gse_eq_2_Idx = gse_eq_Idx[ 2 ]

    gse_Idxs = (; gse_eq_1_Idx,
                gse_eq_2_Idx  )
    
    return (; gse_flat,
            gse_Idxs   )
end
 


########################################################
# ------------------------------------------------------
#  mismatch utitlity functions
# ------------------------------------------------------
########################################################


#-------------------------------------------------------
#-------------------------------------------------------


function update_nodes_current_mismatch!(
    current_mismatch,
    uh,
    (P_Q_gens_view,
     P_Q_non_gens_view,
     P_Q_gens_loc_load_view),
    I_sum_ynj_vj )

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    S_gens_loc_load =
        x_from_xr_xi.( P_Q_gens_loc_load_view )

    S_non_gens =
        x_from_xr_xi.( P_Q_non_gens_view )    

    current_mismatch .=
        I_sum_ynj_vj +
        ( conj.( S_non_gens )) ./ ( conj.( uh )) +
        (conj.( S_gens_loc_load )) ./ ( conj.( uh )) -
        (conj.(S_gens)) ./ ( conj.( uh ))
    
    return nothing
end



function update_nodes_current_mismatch(
    current_mismatch,
    uh,
    (P_Q_non_gens_view,
     P_Q_gens_loc_load_view),
    idq_θ_π_vhθh,
    I_sum_ynj_vj )


    S_gens_loc_load = x_from_xr_xi.(
        P_Q_gens_loc_load_view )

    S_non_gens = x_from_xr_xi.( P_Q_non_gens_view )    

    current_mismatch .=
        I_sum_ynj_vj +
        ( conj.( S_non_gens )) ./ ( conj.( uh )) +
        (conj.( S_gens_loc_load )) ./ ( conj.( uh )) -
        idq_θ_π_vhθh
    
    return nothing
end


function update_nodes_current_mismatch!(
    current_mismatch,
    uh,
    (P_Q_gens_view,
     P_Q_non_gens_view,
     P_Q_gens_loc_load_view),    
    ( gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx ),
    I_sum_ynj_vj;
    loc_load_exist = false )

    if loc_load_exist == true
        
        S_gens = -1 * x_from_xr_xi.( P_Q_gens_view )

        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )

        S_gens_loc_load =
            x_from_xr_xi.( P_Q_gens_loc_load_view )

        for (idx, S_node_type) in zip(            
            [gens_nodes_idx, non_gens_nodes_idx,
             gens_with_loc_load_idx ],
            [ S_gens, S_non_gens, S_gens_loc_load ])

            current_mismatch[idx] .=
                I_sum_ynj_vj[idx] +
                ( (conj.(S_node_type)) ./ (conj.(
                    uh[idx] )) )
        end
        
    else
        
        S_gens = -1 * x_from_xr_xi.( P_Q_gens_view )
        
        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )
        
        for (idx, S_node_type) in zip(            
            [gens_nodes_idx, non_gens_nodes_idx],
            [S_gens, S_non_gens ] )

            current_mismatch[idx] .=
                I_sum_ynj_vj[idx] +
                ((conj.(S_node_type)) ./ (conj.(
                    uh[idx]))) 
        end

        # I_gens =
        #     conj.(S_gens) ./ conj.(uh[gens_nodes_idx])

        # current_mismatch[gens_nodes_idx] .=
        #     I_sum_ynj_vj[gens_nodes_idx] - I_gens

        # I_non_gens =
        #     conj.(S_non_gens ) ./ conj.(uh[non_gens_nodes_idx])

        # current_mismatch[non_gens_nodes_idx] .=
        #     I_sum_ynj_vj[non_gens_nodes_idx] + I_non_gens
        
        
    end
    
    return nothing
end


function update_nodes_current_mismatch(
    current_mismatch,
    uh,
    (P_Q_non_gens_view,
     P_Q_gens_loc_load_view),
    ( gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx ),
    idq_θ_π_vhθh,
    I_sum_ynj_vj;
    loc_load_exist = false )


    if loc_load_exist == true

        current_mismatch[gens_nodes_idx] .=
            I_sum_ynj_vj[gens_nodes_idx] -
            idq_θ_π_vhθh[gens_nodes_idx]

        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )

        S_gens_loc_load =
            x_from_xr_xi.( P_Q_gens_loc_load_view )

        for (idx, S_node_type) in zip(            
            [ non_gens_nodes_idx,
              gens_with_loc_load_idx ],
            [S_non_gens, S_gens_loc_load ] )

            current_mismatch[idx] .=
                I_sum_ynj_vj[idx] +
               ((conj.(S_node_type)) ./(conj.(uh[idx])) )
        end
    else
        
        current_mismatch[gens_nodes_idx] .=
            I_sum_ynj_vj[gens_nodes_idx] -
            idq_θ_π_vhθh[gens_nodes_idx]
        
        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )

        current_mismatch[ non_gens_nodes_idx ] .=
                I_sum_ynj_vj[ non_gens_nodes_idx ] +
                ((conj.(S_non_gens)) ./ (conj.(uh[
                    non_gens_nodes_idx]))) 
        
    end
    
    return nothing
end



function update_nodes_current_mismatch_idq_θπ(
    current_mismatch,
    uh,
    (P_Q_non_gens_view,
     P_Q_gens_loc_load_view),
    ( gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx ),
    idq_θ_π_vhθh,
    I_sum_ynj_vj;
    loc_load_exist = false )


    if loc_load_exist == true

        current_mismatch[gens_nodes_idx] .=
            I_sum_ynj_vj[gens_nodes_idx] -
            idq_θ_π_vhθh[gens_nodes_idx]

        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )

        S_gens_loc_load =
            x_from_xr_xi.( P_Q_gens_loc_load_view )

        for (idx, S_node_type) in zip(            
            [ non_gens_nodes_idx,
              gens_with_loc_load_idx ],
            [S_non_gens, S_gens_loc_load ] )

            current_mismatch[idx] .=
                I_sum_ynj_vj[idx] +
               ((conj.(S_node_type)) ./(conj.(uh[idx])) )
        end
    else
        
        current_mismatch[gens_nodes_idx] .=
            I_sum_ynj_vj[gens_nodes_idx] -
            idq_θ_π_vhθh[gens_nodes_idx]
        
        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )

        current_mismatch[ non_gens_nodes_idx ] .=
                I_sum_ynj_vj[ non_gens_nodes_idx ] +
                ((conj.(S_non_gens)) ./ (conj.(uh[
                    non_gens_nodes_idx]))) 
        
    end
    
    return nothing
end


#-------------------------------------------------------
#-------------------------------------------------------

# i:m where m is the number of gens

function get_dynamic_current_mismatch_gen(
    i, vh, θh,
    id, iq, δ,
    Pl, Ql,
    Ybus,  no_nodes)
    
    return  sum( [Ybus[i][k] * uh[k] for k in 1:no_nodes] ) - ( (id[i] + im *iq[i] ) * exp(im * δ[i] - π/2) + (Pl[i] - im * Ql[i]) / (vh[i] * exp(-im * θh[i])) ) 

end


function get_dynamic_current_mismatch_gen(
    i, vh, θh,
    δ,  ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash,
    Pl, Ql, Ybus, no_nodes)

    idq = get_dynamic_idq_vhθh(
        vh[i], θh[i],
        δ[i], ω[i], ed_dash[i], eq_dash[i],
        ra[i], X_d_dash[i], X_q_dash[i] )
    
    return  sum(
        [Ybus[i][k] * uh[k]
         for k in 1:no_nodes]) - (
             (idq[1] + im *idq[2]) * exp(im * δ[i] - π/2) +
                 (Pl[i] - im * Ql[i]) / (
                     vh[i] * exp( -im * θh[i])) ) 

end


# m+1:no_nodes

function get_dynamic_current_mismatch_load(
    i, vh, θh,
    Pl, Ql,
    Ybus, no_nodes)
    
    return sum(
        [Ybus[i][k] * uh[k]
         for k in 1:no_nodes] ) - (
             Pl[i] - im * Ql[i]) / (
                 vh[i] * exp(-im * θh[i]))  

end



function get_dynamic_current_mismatch_trans()
    
    return 0.0 + im * 0.0

end


function get_dynamic_current_mismatch(
    vh, θh,
    id, iq,
    δ, Pl, Ql,
    Ybus, Dyn_Nodes)
    
    no_nodes = length(collect(values(Dyn_Nodes)))
    
    current_mismatch =  [
        node.Bus_type == :Transmission ?
            get_dynamic_current_mismatch_trans() :
            node.Bus_type == :Load ?
            get_dynamic_current_mismatch_load(
                idx, vh, θh,
                Pl, Ql, Ybus, no_nodes) :
                    get_dynamic_current_mismatch_gen(
                        idx, vh, θh,
                        id, iq, δ, Pl, Ql,
                        Ybus, no_nodes)
        for (idx, node) in
            enumerate( collect(values(Dyn_Nodes)) ) ]

    return vcat( real.( current_mismatch ),
                 imag.( current_mismatch ) )
end

function get_dynamic_current_mismatch(
    vh, θh,
    δ, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash,
    Pl, Ql,
    Ybus, Dyn_Nodes)
    
    no_nodes = length(collect(values(Dyn_Nodes)))
    
    current_mismatch =  [
        node.Bus_type == :Transmission ?
            get_dynamic_current_mismatch_trans() :
            node.Bus_type == :Load ?
            get_dynamic_current_mismatch_load(
                idx, vh, θh,
                Pl, Ql,
                Ybus, no_nodes) :
                    get_dynamic_current_mismatch_gen(
                        idx, vh, θh,
                        δ, ed_dash, eq_dash,
                        ra, X_d_dash, X_q_dash,
                        Pl, Ql, Ybus, no_nodes)
        for (idx, node) in
            enumerate( collect(values(Dyn_Nodes)) ) ]

    return vcat( real.( current_mismatch ),
                 imag.( current_mismatch ) )
end


# -------------------------------------------
# Active gen node power mismatch
# -------------------------------------------


function get_nth_gen_P_pf_mismatch_wt_loc_load(    
    nodes_vh,
    nodes_θh,
    
    δ, id, iq,
    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_P_loc_load,
    nth_idx )
    
    return (id * nodes_vh[nth_idx] *
        sin(δ-nodes_θh[nth_idx])) +
        (iq * nodes_vh[nth_idx] *
        cos(δ-nodes_θh[nth_idx])) -
        nth_P_loc_load -
        sum(
            [nodes_vh[nth_idx] *
                nodes_vh[node_idx] * abs(ynj) *
                cos(
                    nodes_θh[nth_idx] -
                        nodes_θh[node_idx] -
                        angle( ynj ))
             for (ynj, node_idx) in
                 zip( Ynet_nth_row ,
                      nth_node_idx_and_adj_nodes_idx)])
    
end




function get_nth_gen_P_pf_mismatch_no_loc_load(    
    nodes_vh,
    nodes_θh,
    
    δ, id, iq,
    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,    
    nth_idx )

    
    return (id * nodes_vh[nth_idx] *
        sin(δ-nodes_θh[nth_idx])) +
        (iq * nodes_vh[nth_idx] *
        cos(δ - nodes_θh[nth_idx])) -
        sum(
            [nodes_vh[nth_idx] *
                nodes_vh[node_idx] * abs(ynj) *
                cos(
                    nodes_θh[nth_idx] -
                        nodes_θh[node_idx] -
                        angle( ynj ))
             for (ynj, node_idx) in
                 zip( Ynet_nth_row ,
                      nth_node_idx_and_adj_nodes_idx)])
    
end



function get_nth_gen_P_pf_mismatch_wt_loc_load(    
    nodes_vh,
    nodes_θh,    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_gen_P,
    nth_P_loc_load,
    nth_idx )


    # I_net = sum(
    #     [ ynj * nodes_vh[node_idx] *
    #         exp(im * nodes_θh[node_idx])
    #       for (ynj, node_idx) in
    #           zip( Ynet_nth_row ,
    #                nth_node_idx_and_adj_nodes_idx)])
    

    # return nth_gen_P - nth_P_loc_load - real( nodes_vh[nth_idx] *
    #     exp(im * nodes_θh[nth_idx] ) * conj( I_net ))
    
    return nth_gen_P - nth_P_loc_load -
        sum(
            [nodes_vh[nth_idx] *
                nodes_vh[node_idx] * abs(ynj) *
                cos(
                    nodes_θh[nth_idx] -
                        nodes_θh[node_idx] -
                        angle( ynj ))
             for (ynj, node_idx) in
                 zip( Ynet_nth_row ,
                      nth_node_idx_and_adj_nodes_idx)])
    
end



function get_nth_gen_P_pf_mismatch_no_loc_load(    
    nodes_vh,
    nodes_θh,    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_gen_P,
    nth_idx )


    # I_net = sum(
    #     [ ynj * nodes_vh[node_idx] *
    #         exp(im * nodes_θh[node_idx])
    #       for (ynj, node_idx) in
    #           zip( Ynet_nth_row ,
    #                nth_node_idx_and_adj_nodes_idx)])
    

    # return nth_gen_P - real( nodes_vh[nth_idx] *
    #     exp(im * nodes_θh[nth_idx] ) * conj( I_net ))
    
    
    return nth_gen_P -
        sum(
            [nodes_vh[nth_idx] *
                nodes_vh[node_idx] * abs(ynj) *
                cos(
                    nodes_θh[nth_idx] -
                        nodes_θh[node_idx] -
                        angle( ynj ))
             for (ynj, node_idx) in
                 zip( Ynet_nth_row ,
                      nth_node_idx_and_adj_nodes_idx)])
    
end


# -------------------------------------------
# Reactive gen node power mismatch
# -------------------------------------------


function get_nth_gen_Q_pf_mismatch_wt_loc_load(
    nodes_vh,
    nodes_θh,
    
    δ, id, iq,
    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,    
    nth_Q_loc_load,
    nth_idx )
        
    return (id * nodes_vh[nth_idx] *
        cos(δ-nodes_θh[nth_idx])) -
        (iq * nodes_vh[nth_idx] *
        sin(δ-nodes_θh[nth_idx])) -
        nth_Q_loc_load -
        sum(
            [nodes_vh[nth_idx] *
                nodes_vh[node_idx] * abs(ynj) *
                sin( nodes_θh[nth_idx] -
                nodes_θh[node_idx] -
                angle(ynj))
             
             for (ynj, node_idx) in
                 zip( Ynet_nth_row ,
                      nth_node_idx_and_adj_nodes_idx)])     
end



function get_nth_gen_Q_pf_mismatch_no_loc_load(
    nodes_vh,
    nodes_θh,
    
    δ, id, iq,
    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,     
    nth_idx )

    
    return (id * nodes_vh[nth_idx] *
        cos(δ-nodes_θh[nth_idx])) -
        (iq * nodes_vh[nth_idx] *
        sin(δ - nodes_θh[nth_idx])) -
        sum(
            [nodes_vh[nth_idx] *
                nodes_vh[node_idx] * abs(ynj) *
                sin(
                    nodes_θh[nth_idx] -
                        nodes_θh[node_idx] -
                        angle( ynj ))
             for (ynj, node_idx) in
                 zip( Ynet_nth_row ,
                      nth_node_idx_and_adj_nodes_idx)])
     
end



function get_nth_gen_Q_pf_mismatch_wt_loc_load(
    nodes_vh,
    nodes_θh,    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_gen_Q,
    nth_Q_loc_load,
    nth_idx )
        
    return nth_gen_Q - nth_Q_loc_load -
        sum(
            [nodes_vh[nth_idx] *
                nodes_vh[node_idx] * abs(ynj) *
                sin( nodes_θh[nth_idx] -
                nodes_θh[node_idx] -
                angle(ynj))
             
             for (ynj, node_idx) in
                 zip( Ynet_nth_row ,
                      nth_node_idx_and_adj_nodes_idx)])     
end


function get_nth_gen_Q_pf_mismatch_no_loc_load(
    nodes_vh,
    nodes_θh,    
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,
    nth_gen_Q,
    nth_idx )
    
    return nth_gen_Q -
        sum(
            [nodes_vh[nth_idx] *
                nodes_vh[node_idx] * abs(ynj) *
                sin(
                    nodes_θh[nth_idx] -
                        nodes_θh[node_idx] -
                        angle( ynj ))
             for (ynj, node_idx) in
                 zip( Ynet_nth_row ,
                      nth_node_idx_and_adj_nodes_idx)])
     
end

# -------------------------------------------
# Real non gen node power mismatch
# -------------------------------------------


function get_nth_non_gen_P_pf_mismatch( 
    nodes_vh,
    nodes_θh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx,      
    nth_non_gen_P,
    nth_idx )

    # I_net = sum(
    #     [ ynj * nodes_vh[node_idx] *
    #         exp(im * nodes_θh[node_idx])
    #       for (ynj, node_idx) in
    #           zip( Ynet_nth_row ,
    #                nth_node_idx_and_adj_nodes_idx)])
    

    # return nth_non_gen_P - real( nodes_vh[nth_idx] *
    #     exp(im * nodes_θh[nth_idx] ) * conj( I_net ))

    
    return nth_non_gen_P -
        sum([nodes_vh[nth_idx] * nodes_vh[node_idx] * abs(
            ynj) * cos(nodes_θh[nth_idx] -
                nodes_θh[node_idx] - angle( ynj ))
             for (ynj, node_idx) in
                 zip(Ynet_nth_row ,
                     nth_node_idx_and_adj_nodes_idx)])
    
    
end



function get_nth_non_gen_Q_pf_mismatch(
    nodes_vh,
    nodes_θh,
    Ynet_nth_row,
    nth_node_idx_and_adj_nodes_idx, 
    nth_non_gen_Q,
    nth_idx )

    # I_net = sum(
    #     [ ynj * nodes_vh[node_idx] *
    #         exp(im * nodes_θh[node_idx])
    #       for (ynj, node_idx) in
    #           zip( Ynet_nth_row ,
    #                nth_node_idx_and_adj_nodes_idx)])
    

    # return nth_non_gen_Q - imag( nodes_vh[nth_idx] *
    #             exp(im * nodes_θh[nth_idx]) * conj(I_net )) 

    return nth_non_gen_Q -
        sum([nodes_vh[nth_idx] * nodes_vh[node_idx] * abs(
            ynj) * sin(nodes_θh[nth_idx]-
                nodes_θh[node_idx] - angle( ynj ))
         for (ynj, node_idx) in
             zip( Ynet_nth_row ,
                  nth_node_idx_and_adj_nodes_idx)] )

end

# -------------------------------------------
# Complex gen node power mismatch
# -------------------------------------------


function get_a_gen_complex_pf_mismatch(
    
    vh_θh,
    δ_ω_ed_dash_eq_dash,
    ra_X_d_dash_X_q_dash,
    id_iq,    
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    gens_nodes_idx;
    P_loc_load = nothing,
    Q_loc_load = nothing,
    loc_load_exist = false )

    return get_a_gen_complex_pf_mismatch(
        vh_θh...,
        δ_ω_ed_dash_eq_dash...,
        ra_X_d_dash_X_q_dash...,
        id_iq...,    
        uh,
        Y_bus_vec,
        nth_node_idx_and_adj_nodes_idx,
        gens_nodes_idx;
        P_loc_load = P_loc_load,
        Q_loc_load = Q_loc_load,
        loc_load_exist = loc_load_exist )
    
end


function get_a_gen_complex_pf_mismatch(
    
    vh, θh,
    δ_ω_ed_dash_eq_dash,
    ra_X_d_dash_X_q_dash,
    id_iq,    
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    gens_nodes_idx;
    P_loc_load = nothing,
    Q_loc_load = nothing,
    loc_load_exist = false )

    return get_a_gen_complex_pf_mismatch(
        vh, θh,
        δ_ω_ed_dash_eq_dash...,
        ra_X_d_dash_X_q_dash...,
        id_iq...,    
        uh,
        Y_bus_vec,
        nth_node_idx_and_adj_nodes_idx,
        gens_nodes_idx;
        P_loc_load = P_loc_load,
        Q_loc_load = Q_loc_load,
        loc_load_exist = loc_load_exist )
    
end


function get_a_gen_complex_pf_mismatch(
    
    vh, θh,
    δ_ω_ed_dash_eq_dash,
    ra_X_d_dash_X_q_dash,
    id, iq,    
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    gens_nodes_idx;
    P_loc_load = nothing,
    Q_loc_load = nothing,
    loc_load_exist = false )

    return get_a_gen_complex_pf_mismatch(
        vh, θh,
        δ_ω_ed_dash_eq_dash...,
        ra_X_d_dash_X_q_dash...,
        id, iq,    
        uh,
        Y_bus_vec,
        nth_node_idx_and_adj_nodes_idx,
        gens_nodes_idx;
        P_loc_load = P_loc_load,
        Q_loc_load = Q_loc_load,
        loc_load_exist = loc_load_exist )
    
end



function get_a_gen_complex_pf_mismatch(
    
    vh, θh,
    δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash,
    id, iq,    
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    gens_nodes_idx;
    P_loc_load = nothing,
    Q_loc_load = nothing,
    loc_load_exist = false )

    first_node_idx = first(
        nth_node_idx_and_adj_nodes_idx)

    # gen_uh = uh[ first_node_idx ]

    gen_uh = vh * exp( im * θh )

    if first_node_idx ∈ gens_nodes_idx
        if loc_load_exist == true

            loc_load_S =
                P_loc_load + im * Q_loc_load

            a_gen_inj =
                vh * exp(im * θh) * (id -im * iq ) *
                exp(-im * (δ - π/2)) 

            a_gen_I_sum_ynj_vj = sum(
                [ ynj * uh[node_idx]
                  for (ynj, node_idx) in
                      zip( Y_bus_vec,
                           nth_node_idx_and_adj_nodes_idx)])

            return  loc_load_S +
                gen_uh * conj(a_gen_I_sum_ynj_vj) -
                a_gen_inj

        else

            a_gen_inj =
                vh * exp(im * θh) * (id -im * iq ) *
                exp(-im * (δ - π/2)) 

            a_gen_I_sum_ynj_vj = sum(
                [ ynj * uh[node_idx]
                  for (ynj, node_idx) in
                      zip( Y_bus_vec,
                           nth_node_idx_and_adj_nodes_idx)])

            return   gen_uh  * conj(a_gen_I_sum_ynj_vj) -
                a_gen_inj

        end
    else
        nothing
    end
    
end

# -------------------------------------------


function get_a_gen_complex_pf_mismatch_by_id_iq(
    id_iq,    
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    gens_nodes_idx;
    P_loc_load = nothing,
    Q_loc_load = nothing,
    loc_load_exist = false )

    return get_a_gen_complex_pf_mismatch_by_id_iq(
        id_iq,    
        uh,
        Y_bus_vec,
        nth_node_idx_and_adj_nodes_idx,
        gens_nodes_idx;
        P_loc_load = P_loc_load,
        Q_loc_load = Q_loc_load,
        loc_load_exist = loc_load_exist )
end


function get_a_gen_complex_pf_mismatch_by_id_iq(
    id, iq,    
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    gens_nodes_idx;
    P_loc_load = nothing,
    Q_loc_load = nothing,
    loc_load_exist = false )

    first_node_idx = first(
        nth_node_idx_and_adj_nodes_idx)

    gen_uh = uh[ first_node_idx  ]

    if first_node_idx  ∈ gens_nodes_idx

        if loc_load_exist == true

            loc_load_S =
                P_loc_load + im * Q_loc_load

            a_gen_inj =
                gen_uh * (id -im * iq ) *
                exp(-im * (δ - π/2)) 

            a_gen_I_sum_ynj_vj = sum(
                [ ynj * uh[node_idx]
                  for (ynj, node_idx) in
                      zip( Y_bus_vec,
                           nth_node_idx_and_adj_nodes_idx)])

            return  loc_load_S +
                gen_uh *  conj(a_gen_I_sum_ynj_vj) -
                a_gen_inj
        else

            a_gen_inj =
                gen_uh * (id -im * iq ) *
                exp(-im * (δ - π/2)) 

            a_gen_I_sum_ynj_vj = sum(
                [ ynj * uh[node_idx]
                  for (ynj, node_idx) in
                      zip( Y_bus_vec,
                           nth_node_idx_and_adj_nodes_idx)])
            return   gen_uh * conj(a_gen_I_sum_ynj_vj) -
                a_gen_inj

        end        
    else
        nothing
    end

end

# -------------------------------------------

    
function get_a_gen_complex_pf_mismatch_by_δ_ed_eq(
    δ_ω_ed_dash_eq_dash,
    ra_X_d_dash_X_q_dash,   
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    gens_nodes_idx;
    P_loc_load = nothing,
    Q_loc_load = nothing,
    loc_load_exist = false )

    return get_a_gen_complex_pf_mismatch_by_δ_ed_eq(
        δ_ω_ed_dash_eq_dash...,
        ra_X_d_dash_X_q_dash...,   
        uh,
        Y_bus_vec,
        nth_node_idx_and_adj_nodes_idx,
        gens_nodes_idx;
        P_loc_load = P_loc_load,
        Q_loc_load = Q_loc_load,
        loc_load_exist = loc_load_exist )

end



function get_a_gen_complex_pf_mismatch_by_δ_ed_eq(
    δ, ω, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash,   
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    gens_nodes_idx;
    P_loc_load = nothing,
    Q_loc_load = nothing,
    loc_load_exist = false )

    first_node_idx = first(
        nth_node_idx_and_adj_nodes_idx)
    
    gen_uh = uh[ first_node_idx ]
    
    if first_node_idx ∈ gens_nodes_idx

        vh = abs( gen_uh )

        θh = angle( gen_uh )

        id_iq = get_a_gen_dyn_idq(
            vh, θh,
            δ, ω, ed_dash, eq_dash,
            ra, X_d_dash, X_q_dash )

        id, iq = id_iq

        if loc_load_exist == true

            loc_load_S =
                P_loc_load + im * Q_loc_load

            a_gen_inj =
                gen_uh * (id -im * iq ) *
                exp(-im * (δ - π/2)) 

            a_gen_I_sum_ynj_vj = sum(
                [ ynj * uh[node_idx]
                  for (ynj, node_idx) in
                      zip( Y_bus_vec,
                           nth_node_idx_and_adj_nodes_idx)])

            return  loc_load_S +
                gen_uh * conj(a_gen_I_sum_ynj_vj) -
                a_gen_inj
        else

            a_gen_inj =
                gen_uh * (id -im * iq ) *
                exp(-im * (δ - π/2)) 

            a_gen_I_sum_ynj_vj = sum(
                [ ynj * uh[node_idx]
                  for (ynj, node_idx) in
                      zip( Y_bus_vec,
                           nth_node_idx_and_adj_nodes_idx)])
            return   gen_uh * conj(a_gen_I_sum_ynj_vj) -
                a_gen_inj

        end        
        
    else
        nothing
    end

end



# -------------------------------------------
# Complex load node power mismatch
# -------------------------------------------

function get_a_load_complex_pf_mismatch(   
    P_load_Q_load,
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    non_gens_nodes_idx )
    
    return get_a_load_complex_pf_mismatch(   
        P_load_Q_load...,
        uh,
        Y_bus_vec,
        nth_node_idx_and_adj_nodes_idx,
        non_gens_nodes_idx )    

end



function get_a_load_complex_pf_mismatch(   
    P_load,
    Q_load,
    uh,
    Y_bus_vec,
    nth_node_idx_and_adj_nodes_idx,
    non_gens_nodes_idx )

    first_node_idx = first(
        nth_node_idx_and_adj_nodes_idx)

    load_uh = uh[ first_node_idx ]

    if first_node_idx ∈ non_gens_nodes_idx

        load_S =
            P_load + im * Q_load

        a_load_I_sum_ynj_vj = sum(
            [ ynj * uh[node_idx]
              for (ynj, node_idx) in
                  zip( Y_bus_vec,
                       nth_node_idx_and_adj_nodes_idx ) ])

        return  load_S + load_uh  *
            conj(a_load_I_sum_ynj_vj)
        
    else
        
        return nothing
        
    end
    

end


# -------------------------------------------



function get_a_gen_real_stator_equations_mismatch(
    vh, θh,
    δ, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash,
    id, iq )

    return ed_dash - vh * sin(δ - θh) -
        ra * id + X_q_dash * iq

end



function get_a_gen_imag_stator_equations_mismatch(
    vh, θh,
    δ, ed_dash, eq_dash,
    ra, X_d_dash, X_q_dash,
    id, iq )

    return eq_dash - vh * cos(δ - θh) -
        ra * iq - X_d_dash * id 

end

# -------------------------------------------
# -------------------------------------------

function get_current_mismatch_nodes_no_loc_load(
    uh,
    (PQ_gens_pad,
     PQ_non_gens_pad,
    I_sum_ynj_vj ))
    
    return  I_sum_ynj_vj +
        (conj.(x_from_xr_xi.(PQ_non_gens_pad))) ./ (
            conj.(uh)) - (
                conj.(x_from_xr_xi.(PQ_gens_pad))) ./ (
                    conj.(uh))
    
end


function get_current_mismatch_nodes_with_loc_load(
    uh,
    (PQ_gens_pad,
     PQ_non_gens_pad,
     PQ_loc_loads_pad,
    I_sum_ynj_vj ))
    
    return  I_sum_ynj_vj +
        (conj.(x_from_xr_xi.(PQ_non_gens_pad ))) ./ (
            conj.(uh)) + (
                conj.(x_from_xr_xi.(PQ_loc_loads_pad))) ./ (
                    conj.(uh)) - (conj.(
                        x_from_xr_xi.(PQ_gens_pad))) ./ (
                            conj.( uh ))
    
end



function get_current_mismatch_by_idq_no_loc_load(
    uh,
    ( PQ_non_gens_pad,
    idq_wt_pad_view,
    I_sum_ynj_vj ))
    
    return  I_sum_ynj_vj +
        (conj.( x_from_xr_xi.(PQ_loc_loads_pad))) ./ (
            conj.( uh ))  - idq_wt_pad_view
end


function get_current_mismatch_by_idq_with_loc_load(
    uh,
    ( PQ_non_gens_pad,
      PQ_loc_loads_pad,
    idq_wt_pad_view,
    I_sum_ynj_vj ))
    
    return  I_sum_ynj_vj +
        (conj.( x_from_xr_xi.(PQ_non_gens_pad))) ./ (
            conj.(uh)) + (
                conj.(x_from_xr_xi.(PQ_loc_loads_pad))) ./ (
                    conj.( uh )) - idq_wt_pad_view
end



function get_nodes_current_mismatch_idq_θπ(
    uh,
    ( P_Q_non_gens_view,
      P_Q_gens_loc_load_view ),
    idq_θ_π_vhθh,
    I_sum_ynj_vj )

    S_gens_loc_load =
        x_from_xr_xi.( P_Q_gens_loc_load_view )

    S_non_gens =
        x_from_xr_xi.( P_Q_non_gens_view )    

    # I_sum_ynj_vj + ( conj.( S_non_gens )) ./ ( conj.( uh )) + (conj.( S_gens_loc_load )) ./ ( conj.( uh )) - idq_θ_π_vhθh
    
    return  I_sum_ynj_vj +
        (conj.(x_from_xr_xi.(P_Q_non_gens_view))) ./ (
            conj.(uh)) + (conj.(
                x_from_xr_xi.(P_Q_gens_loc_load_view))) ./ (
                    conj.( uh )) - idq_θ_π_vhθh
end


function get_nodes_current_mismatch(
    uh,
    ( P_Q_non_gens_view,
      P_Q_gens_loc_load_view ),
    idq_θ_π_vhθh,
    I_sum_ynj_vj )

    S_gens_loc_load =
        x_from_xr_xi.( P_Q_gens_loc_load_view )

    S_non_gens =
        x_from_xr_xi.( P_Q_non_gens_view )    

    return  I_sum_ynj_vj +
        ( conj.( S_non_gens )) ./ (
            conj.(uh)) + (conj.( S_gens_loc_load )) ./ (
                conj.(uh)) - idq_θ_π_vhθh
end


function get_nodes_current_mismatch(
    uh,
    (P_Q_gens_view,
     P_Q_non_gens_view,
     P_Q_gens_loc_load_view),
     I_sum_ynj_vj )

    S_gens = x_from_xr_xi.( P_Q_gens_view )

    S_gens_loc_load =
        x_from_xr_xi.( P_Q_gens_loc_load_view )

    S_non_gens =
        x_from_xr_xi.( P_Q_non_gens_view )    

    return  I_sum_ynj_vj +
        (conj.( S_non_gens )) ./ (
            conj.(uh)) + (
                conj.(S_gens_loc_load)) ./ (
                    conj.( uh )) - (
                        conj.(S_gens)) ./ ( conj.( uh ))
    
end


# function get_nodes_current_mismatch(
#     uh,
#     (P_Q_gens_view,
#      P_Q_non_gens_view,
#      P_Q_gens_loc_load_view),
#     ( gens_nodes_idx,
#       non_gens_nodes_idx,
#       gens_with_loc_load_idx ),
#     I_sum_ynj_vj;
#     loc_load_exist = false)

#     nodes_size = sum( [length(gens_nodes_idx),
#                        length(non_gens_nodes_idx) ])

#     # This should not be zero initially
    
#     current_mismatch = ones(ComplexF64, nodes_size )

#     if loc_load_exist == true
        
#         S_gens = -1 * x_from_xr_xi.( P_Q_gens_view )

#         S_non_gens =
#             x_from_xr_xi.( P_Q_non_gens_view )

#         S_gens_loc_load =
#             x_from_xr_xi.( P_Q_gens_loc_load_view )

#         for (idx, S_node_type) in zip(
#             [gens_nodes_idx, non_gens_nodes_idx,
#              gens_with_loc_load_idx ],
#             [S_gens, S_non_gens,
#              S_gens_loc_load ] )

#             current_mismatch[idx] .= 
#                 I_sum_ynj_vj[idx] .+
#                ((conj.(S_node_type)) ./ (conj.( uh[idx])))
#         end
#     else
        
#         # S_gens = -1 * x_from_xr_xi.( P_Q_gens_view )

#         # S_non_gens =
#         #     x_from_xr_xi.( P_Q_non_gens_view )

#         # for (idx, S_node_type) in zip(
#         #     [gens_nodes_idx, non_gens_nodes_idx],
#         #     [S_gens, S_non_gens ] )

#         #     current_mismatch[idx] .= I_sum_ynj_vj[idx] +
#         #        ((conj.(S_node_type))./(conj.(uh[idx]))) 
#         # end

        
#         S_gens =
#             x_from_xr_xi.(P_Q_gens_view )

#         I_gens =
#             conj.(S_gens) ./ conj.(uh[gens_nodes_idx])

#         current_mismatch[gens_nodes_idx] .=
#             I_sum_ynj_vj[gens_nodes_idx] - I_gens
        
#         S_non_gens =
#             x_from_xr_xi.(P_Q_non_gens_view )

#         I_non_gens =
#             conj.(S_non_gens) ./ conj.(
#                 uh[non_gens_nodes_idx])

#         current_mismatch[non_gens_nodes_idx] .=
#             I_sum_ynj_vj[non_gens_nodes_idx] + I_non_gens
        
#     end

#     return current_mismatch
    
    
# end



function get_nodes_current_mismatch(
    uh,
    P_gens,
    Q_gens,
    P_non_gens,
    Q_non_gens,
    P_gens_loc_load,
    Q_gens_loc_load,    
    ( gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx),
    ( n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs ),
    I_sum_ynj_vj;
    loc_load_exist = false)
        
    # This should not be zero initially
        
    S_gens = P_gens + im * Q_gens

    S_non_gens = P_non_gens + im * Q_non_gens

    if loc_load_exist == true
            
        S_gens_loc_load =
            P_gens_loc_load + im * Q_gens_loc_load

    end
    
    return [ idx ∈ gens_with_loc_load_idx ?
        I_sum_ynj_vj[idx] + ((conj(S_gens_loc_load[
            n2s_gens_with_loc_load_idxs[idx]]))/(
                conj(uh[idx]))) - ((conj(S_gens[
                    n2s_gens_idx[idx]]))/(
                        conj(uh[idx]))) :
                            idx ∈ non_gens_nodes_idx ?
                            I_sum_ynj_vj[idx] + ((
                                conj(S_non_gens[
                            n2s_non_gens_idx[idx]]))/(
                                conj(uh[idx]))) :
                                    I_sum_ynj_vj[idx] - ((
                                        conj(S_gens[
                                            n2s_gens_idx[
                                                idx]]))/(
                                                    conj(uh[
                                                        idx])))
             for idx in all_nodes_idx]


    # nodes_size = length(all_nodes_idx)
    
    # current_mismatch = ones(ComplexF64, nodes_size )
    
    # if loc_load_exist == true
    #     S_gens = P_gens + im * Q_gens_view 
    #     S_non_gens = P_non_gens + im * Q_non_gens
    #     S_gens_loc_load =
    #         P_gens_loc_load + im * Q_gens_loc_load
    #     for idx in all_nodes_idx
    #         if idx ∈ non_gens_nodes_idx
    #             current_mismatch[idx] = 
    #             I_sum_ynj_vj[idx] +
    #             ((conj(S_non_gens[
    #                 n2s_non_gens_idx[idx]])) / (
    #                     conj( uh[idx])))
    #         elseif idx ∈ gens_with_loc_load_idx
    #             current_mismatch[idx] =
    #                 I_sum_ynj_vj[idx] +
    #                 ((conj(S_gens_loc_load[
    #                     n2s_gens_with_loc_load_idxs[
    #                         idx]]))/(conj( uh[idx]))) -
    #                             ((conj(S_gens[
    #                                 n2s_gens_idx[
    #                                     idx]])) / (
    #                                         conj(uh[
    #                                             idx])))  
    #         else
    #             current_mismatch[idx] =
    #                 I_sum_ynj_vj[idx] - ((conj(S_gens[
    #                     n2s_gens_idx[idx]]))/(
    #                         conj(uh[ idx])))
    #         end
    #     end
    # else
    #     S_gens = P_gens + im * Q_gens_view 
    #     S_non_gens = P_non_gens + im * Q_non_gens        
    #     for idx in all_nodes_idx
    #         if idx ∈ non_gens_nodes_idx
    #             current_mismatch[idx] = 
    #             I_sum_ynj_vj[idx] +
    #             ((conj(S_non_gens[
    #                 n2s_non_gens_idx[idx]])) / (
    #                     conj( uh[idx])))
    #         else
    #             current_mismatch[idx] =
    #                 I_sum_ynj_vj[idx] - ((conj(S_gens[
    #                     n2s_gens_idx[idx]])) / (
    #                         conj(uh[ idx])))                
    #         end
    #     end
    # end
    
    # return current_mismatch
        
end


function get_nodes_current_mismatch(
    vh,
    θh,
    P_gens,
    Q_gens,
    P_non_gens,
    Q_non_gens,
    P_gens_loc_load,
    Q_gens_loc_load,    
    ( gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx),
    ( n2s_gens_idx,
      n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx #
      ),
    I_sum_ynj_vj;
    loc_load_exist = false)

    nodes_size = length(all_nodes_idx)
        
    S_gens = P_gens + im * Q_gens 

    S_non_gens = P_non_gens + im * Q_non_gens

    if loc_load_exist == true
            
        S_gens_loc_load =
            P_gens_loc_load + im * Q_gens_loc_load

    end
    
    return [ idx ∈ gens_with_loc_load_idx ?
        I_sum_ynj_vj[ n2s_all_nodes_idx[idx]] + (
            (conj(S_gens_loc_load[
                n2s_gens_with_loc_load_idxs[idx]]))/(
                    vh[ n2s_all_nodes_idx[idx] ] * exp(
                        -im*θh[ n2s_all_nodes_idx[idx]]))) -
                            ((conj(S_gens[
                    n2s_gens_idx[idx]]))/(
                        vh[ n2s_all_nodes_idx[idx]] * exp(
                            -im*θh[
                                n2s_all_nodes_idx[idx]]))) :
                        idx ∈ non_gens_nodes_idx ?
                        I_sum_ynj_vj[
                            n2s_all_nodes_idx[idx]] + (
                            (conj(S_non_gens[
                            n2s_non_gens_idx[idx]]))/(
                                vh[n2s_all_nodes_idx[idx]] *
                                    exp(-im*θh[ n2s_all_nodes_idx[idx]]))) :
                                    I_sum_ynj_vj[ n2s_all_nodes_idx[idx]] -
                                    ((conj( S_gens[n2s_gens_idx[ idx]]))/(
                                        vh[ n2s_all_nodes_idx[idx]] *
                                            exp(-im*θh[
                                                n2s_all_nodes_idx[idx]])))
             for idx in all_nodes_idx ]
        
end



# function get_nodes_current_mismatch_idq_θπ(
#     uh,
#     ( P_Q_non_gens_view,
#       P_Q_gens_loc_load_view ),
#     ( gens_nodes_idx,
#         non_gens_nodes_idx,
#       gens_with_loc_load_idx ),
#     idq_θ_π_vhθh,
#     I_sum_ynj_vj;
#     loc_load_exist = false)

#     nodes_size =
#         sum(
#             [length(gens_nodes_idx),
#              length(non_gens_nodes_idx) ])

#     # This should not be zero initially
    
#     current_mismatch = ones(ComplexF64, nodes_size )


#     if loc_load_exist == true

#         current_mismatch[gens_nodes_idx] .=
#             I_sum_ynj_vj[gens_nodes_idx] -
#             idq_θ_π_vhθh[gens_nodes_idx]

#         S_non_gens =
#             x_from_xr_xi.( P_Q_non_gens_view )

#         S_gens_loc_load =
#             x_from_xr_xi.( P_Q_gens_loc_load_view )

#         for (idx, S_node_type) in zip(
#             [ non_gens_nodes_idx,
#              gens_with_loc_load_idx ],
#             [ S_non_gens,  S_gens_loc_load ] )

#             current_mismatch[idx] .=
#                 I_sum_ynj_vj[idx] .+
#                ((conj.(S_node_type)) ./(conj.(uh[idx])) )
#         end
#     else
        
#         current_mismatch[gens_nodes_idx] .=
#             I_sum_ynj_vj[gens_nodes_idx] -
#             idq_θ_π_vhθh[gens_nodes_idx]
        
#         S_non_gens =
#             x_from_xr_xi.( P_Q_non_gens_view )

#         current_mismatch[ non_gens_nodes_idx ] .=
#                 I_sum_ynj_vj[ non_gens_nodes_idx ] .+
#                 ((conj.(S_non_gens)) ./ (conj.(uh[non_gens_nodes_idx]))) 
        
#     end

#     return current_mismatch
# end



function get_nodes_current_mismatch_idq_θπ(
    uh,
    ( P_Q_non_gens_view,
      P_Q_gens_loc_load_view ),
    ( gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx ),
    idq_θ_π_vhθh,
    I_sum_ynj_vj;
    loc_load_exist = false)

    nodes_size =
        sum(
            [length(gens_nodes_idx),
             length(non_gens_nodes_idx) ])

    # This should not be zero initially
    
    current_mismatch = ones(ComplexF64, nodes_size )

    if loc_load_exist == true

        current_mismatch[gens_nodes_idx] .=
            I_sum_ynj_vj[gens_nodes_idx] -
            idq_θ_π_vhθh[gens_nodes_idx]

        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )

        S_gens_loc_load =
            x_from_xr_xi.( P_Q_gens_loc_load_view )

        for (idx, S_node_type) in zip(            
            [ non_gens_nodes_idx, gens_with_loc_load_idx ],
            [S_non_gens, S_gens_loc_load ] )

            current_mismatch[idx] .=
                I_sum_ynj_vj[idx] .+
                ((conj.(S_node_type)) ./ (conj.(uh[idx])) )
            
        end
    else
        
        current_mismatch[gens_nodes_idx] .=
            I_sum_ynj_vj[gens_nodes_idx] -
            idq_θ_π_vhθh[gens_nodes_idx]
        
        S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )

        current_mismatch[ non_gens_nodes_idx ] .=
                I_sum_ynj_vj[ non_gens_nodes_idx ] +
                ((conj.(S_non_gens[non_gens_nodes_idx])) ./(conj.(uh[non_gens_nodes_idx]))) 
        
    end

    return current_mismatch
end



function get_nodes_current_mismatch_idq_θπ(
    vh,
    θh,
    ( P_Q_non_gens_view,
      P_Q_gens_loc_load_view ),
    ( gens_nodes_idx,
      non_gens_nodes_idx,
      gens_with_loc_load_idx,
      all_nodes_idx),
     ( n2s_gens_idx,
      n2s_non_gens_idx,
      n2s_gens_with_loc_load_idxs,
      n2s_all_nodes_idx 
      ),
    idq_θ_π_vhθh,
    I_sum_ynj_vj;
    loc_load_exist = false)

    S_non_gens =
            x_from_xr_xi.( P_Q_non_gens_view )
    
    if loc_load_exist == true

        S_gens_loc_load =
            x_from_xr_xi.( P_Q_gens_loc_load_view )
    end
    

    return [ idx ∈ non_gens_nodes_idx ?
        I_sum_ynj_vj[n2s_all_nodes_idx[idx]] +
        ((conj(S_non_gens[n2s_non_gens_idx[idx]]))/(
            vh[n2s_all_nodes_idx[idx]] * exp(
                -im * θh[n2s_all_nodes_idx[idx]]))) :
                    idx ∈ gens_with_loc_load_idx ?
                    I_sum_ynj_vj[n2s_all_nodes_idx[idx]] +
                    ((conj(S_gens_loc_load[
                     n2s_gens_with_loc_load_idxs[idx]]))/
                   (vh[n2s_all_nodes_idx[idx]] *
                   exp(-im * θh[n2s_all_nodes_idx[idx]]))) -
                   idq_θ_π_vhθh[ n2s_all_nodes_idx[idx]] :
                   I_sum_ynj_vj[ n2s_all_nodes_idx[idx]] -
                   idq_θ_π_vhθh[ n2s_all_nodes_idx[idx] ]
      for idx in all_nodes_idx]

end


#-----------------------------------------------------
#-----------------------------------------------------

function diag_id_mismatch(
    gens_vh,
    gens_θh,
    gens_δ,
    gens_ed_dash,
    gens_i_d,
    gens_i_q;
    kwd_para = (ra,
    X_d,
    X_q,
    X_d_dash,
    X_q_dash,
    
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)  )


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

    
    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """
    
    return [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin( gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    

end



function diag_iq_mismatch(
    gens_vh,
    gens_θh,
    gens_δ,
    gens_eq_dash,
    gens_i_d,
    gens_i_q;
    
    kwd_para = (ra,
    X_d,
    X_q,
    X_d_dash,
    X_q_dash,
    
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)  )


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

    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """
    
    return [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[ n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  
    
    

end


function diag_P_mismatch(
    vh,
    θh,
    gens_δ,
    gens_i_d,
    gens_i_q;
    
    kwd_para =
    (Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)
    
    )

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
    
    return [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (gens_i_d[ n2s_gens_idx[ nth_idx]] * vh[ n2s_all_nodes_idx[ nth_idx]] *
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          gens_i_q[ n2s_gens_idx[ nth_idx]] * vh[ n2s_all_nodes_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
        (gens_i_d[ n2s_gens_idx[nth_idx]] * vh[ n2s_all_nodes_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         gens_i_q[ n2s_gens_idx[nth_idx]] * vh[ n2s_all_nodes_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ]))
                        for nth_idx in all_nodes_idx ] 

    
    

end



function diag_Q_mismatch(
    vh,
    θh,
    gens_δ,
    gens_i_d,
    gens_i_q;
    
    kwd_para =
    (Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)
    
    )


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
    
    return [ nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (gens_i_d[ n2s_gens_idx[nth_idx]] * vh[ n2s_all_nodes_idx[nth_idx]] *
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ n2s_all_nodes_idx[nth_idx]] *
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
           (gens_i_d[ n2s_gens_idx[nth_idx]] * vh[ n2s_all_nodes_idx[nth_idx]] *
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            gens_i_q[ n2s_gens_idx[nth_idx] ] * vh[ n2s_all_nodes_idx[nth_idx]] *
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) 
                              for nth_idx in all_nodes_idx ]


    
end



function diag_P_mismatch(
    vh,
    θh,
    gens_δ,
    gens_i_d,
    gens_i_q,
    nothing;
    
    kwd_para =
    (Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)
    
    )


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
    
    return [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * (gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -          
          sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]))  -
                                   P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]]) :
        (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]) ))
                        for nth_idx in all_nodes_idx ] 

    
end


function diag_Q_mismatch(
    vh,
    θh,
    gens_δ,
    gens_i_d,
    gens_i_q,
    nothing;
    
    kwd_para =
    (Ynet,
    nodes_idx_with_adjacent_nodes_idx,
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs)
    
    )
    
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

    return [ nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          ( vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -           
           sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]]) :
           (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])]))) 
                              for nth_idx in all_nodes_idx ]



end


# ------------------------------------------------------
# From sd-dynamics-reduced-order-models
# ---------------------------------------------------

function get_a_gen_Pg_from_δ_i_dq_vh_θh(
    vh, θh,
    δ,
    i_d, i_q )

    "Pg_i = i_d * vh * sin(δ - θh) +
            i_q * vh * cos(δ - θh) "

    return i_d * vh * sin(δ - θh) +
        i_q * vh * cos(δ - θh) 

end


function get_a_gen_Qg_from_δ_i_dq_vh_θh(
    vh, θh,
    δ,
    i_d, i_q )

    
    "Qg_i = i_d * vh * cos(δ - θh) 
            i_q * vh * sin(δ - θh) "


    return i_d * vh * cos(δ - θh) -
        i_q * vh * sin(δ - θh)
    
end


function get_gens_Pg_from_δ_i_dq_vh_θh(
    gens_vh, gens_θh,
    gens_δ,
    gens_id, gens_iq )

    return [ get_a_gen_Pg_from_δ_i_dq_vh_θh(
        vh, θh, δ, i_d, i_q )
           for (vh, θh, δ, i_d, i_q) in
               zip(gens_vh, gens_θh,
                   gens_δ,
                   gens_id, gens_iq ) ]

end


function get_gens_Qg_from_δ_i_dq_vh_θh(
    gens_vh, gens_θh,
    gens_δ,
    gens_id, gens_iq )

    return [ get_a_gen_Qg_from_δ_i_dq_vh_θh(
        vh, θh, δ, i_d, i_q )
           for (vh, θh, δ, i_d, i_q) in
               zip(gens_vh, gens_θh,
                   gens_δ,
                   gens_id, gens_iq ) ]
end


function get_id_iq_vec_by_gens_state_vθ_X_para(
    δ, ed_dash, eq_dash,
    gens_vh, gens_θh,
    ra, X_d_dash, X_q_dash;
    kwd_para = id_iq_Idx )

    (id_Idx,
     iq_Idx) = id_iq_Idx

    diag_Z_d_q =
        vcat(
            hcat( (Diagonal(ra), Diagonal(-X_q_dash)) ),
            hcat( (Diagonal(X_d_dash), Diagonal(ra)) ) )

    b = vcat(ed_dash - vh .* sin.(δ - gens_θh),
             eq_dash - vh .* cos.(δ - gens_θh) )

    vec_id_iq = diag_Z_d_q \  b

    gens_i_d = vec_id_iq[id_Idx]
    
    gens_i_q = vec_id_iq[iq_Idx]

    return (;gens_id, gens_iq )

    
end

function get_dynamic_model_Pq_Qg_id_iq(
    gens_δ, gens_ed_dash, gens_eq_dash,
    gens_vh, gens_θh,
    gens_ra, gens_X_d_dash, gens_X_q_dash;
    kwd_para = id_iq_Idx )

    (;gens_id, gens_iq ) =
        get_id_iq_vec_by_gens_state_vθ_X_para(
            gens_δ, gens_ed_dash, gens_eq_dash,
            gens_vh, gens_θh,
            gens_ra, gens_X_d_dash, gens_X_q_dash;
            kwd_para = id_iq_Idx )
    
    P_gens = get_gens_Pg_from_δ_i_dq_vh_θh(
        gens_vh, gens_θh, gens_δ,
        gens_id, gens_iq )
    
     Q_gens = get_gens_Qg_from_δ_i_dq_vh_θh(
         gens_vh, gens_θh, gens_δ,
         gens_id, gens_iq )

    return (;P_gens, Q_gens, gens_id, gens_iq)
end


function get_dynamic_model_Pq_Qg(
    gens_δ, gens_ed_dash, gens_eq_dash,
    gens_vh, gens_θh,
    gens_ra, gens_X_d_dash, gens_X_q_dash,
    gens_id, gens_iq )
    
    P_gens = get_gens_Pg_from_δ_i_dq_vh_θh(
        gens_vh, gens_θh, gens_δ,
        gens_id, gens_iq )
    
     Q_gens = get_gens_Qg_from_δ_i_dq_vh_θh(
         gens_vh, gens_θh, gens_δ,
         gens_id, gens_iq )

    return (;P_gens, Q_gens)
end

# ------------------------------------------------------

function normalise_angle_θh(
    θh, slack_gens_nodes_idx )

    return θh .- minimum(θh[slack_gens_nodes_idx] )

end


function get_Pg_inj_Qg_inj_Png_Qng(
    Pg_Qg_Png_Qng_Pll_Qll;
    loc_load_exist,
    Pg_Qg_Png_Qng_Pll_Qll_Idx,    
    gens_nodes_idx,
    n2s_gens_idx,
    gens_nodes_with_loc_loads_idx,
    n2s_gens_with_loc_load_idxs)

    
     (Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
          NamedTupleTools.select(
              Pg_Qg_Png_Qng_Pll_Qll_Idx,
              (:dyn_P_gens_Idxs,
               :dyn_Q_gens_Idxs,
               :dyn_P_non_gens_Idxs,
               :dyn_Q_non_gens_Idxs,
               :dyn_P_gens_loc_load_Idxs,
               :dyn_Q_gens_loc_load_Idxs))
    
    
    P_gens =
        Pg_Qg_Png_Qng_Pll_Qll[ Pg_Idx ]

    Q_gens =
        Pg_Qg_Png_Qng_Pll_Qll[ Qg_Idxs ]

    P_non_gens  =
        Pg_Qg_Png_Qng_Pll_Qll[ Png_Idxs ]

    Q_non_gens = 
        Pg_Qg_Png_Qng_Pll_Qll[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            Pg_Qg_Png_Qng_Pll_Qll[ Pgll_Idxs ]

        Q_g_loc_load =
            Pg_Qg_Png_Qng_Pll_Qll[ Qgll_Idxs ]
        
    else

        P_g_loc_load = []
            
        Q_g_loc_load = []
    end

    gens_Pg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        P_gens[ n2s_gens_idx[ idx]] -
        P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] :
        P_gens[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]

    gens_Qg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        Q_gens[ n2s_gens_idx[ idx]] -
        Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] :
        Q_gens[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]


    Pg_inj_Qg_inj_Png_Qng =
        [gens_Pg_inj;
         gens_Qg_inj;
         P_non_gens;
         Q_non_gens]


    Pg_inj_Qg_inj_Png_Qng_Pll_Qll =
        [gens_Pg_inj;
         gens_Qg_inj;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    

    Pg_inj_Png_Qng =
        [gens_Pg_inj;
         P_non_gens;
         Q_non_gens]    
    
    return (;gens_Pg_inj,
            gens_Qg_inj,
            P_gens,
            Q_gens,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,
            Pg_inj_Png_Qng,
            Pg_inj_Qg_inj_Png_Qng,
            Pg_inj_Qg_inj_Png_Qng_Pll_Qll )
end



function get_Pg_inj_Qg_inj_Png_Qng_Pll_Qll(
    Pg_Qg_Png_Qng_Pll_Qll;
    loc_load_exist,
    Pg_Qg_Png_Qng_Pll_Qll_Idx,    
    gens_nodes_idx,
    n2s_gens_idx,
    gens_nodes_with_loc_loads_idx,
    n2s_gens_with_loc_load_idxs)

    
     (Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
          NamedTupleTools.select(
              Pg_Qg_Png_Qng_Pll_Qll_Idx,
              (:dyn_P_gens_Idxs,
               :dyn_Q_gens_Idxs,
               :dyn_P_non_gens_Idxs,
               :dyn_Q_non_gens_Idxs,
               :dyn_P_gens_loc_load_Idxs,
               :dyn_Q_gens_loc_load_Idxs))
    
    
    P_gens =
        Pg_Qg_Png_Qng_Pll_Qll[ Pg_Idx ]

    Q_gens =
        Pg_Qg_Png_Qng_Pll_Qll[ Qg_Idxs ]

    P_non_gens  =
        Pg_Qg_Png_Qng_Pll_Qll[ Png_Idxs ]

    Q_non_gens = 
        Pg_Qg_Png_Qng_Pll_Qll[ Qng_Idxs ]
    
    if loc_load_exist == true

        P_g_loc_load =
            Pg_Qg_Png_Qng_Pll_Qll[ Pgll_Idxs ]

        Q_g_loc_load =
            Pg_Qg_Png_Qng_Pll_Qll[ Qgll_Idxs ]
        
    else

        P_g_loc_load = []
            
        Q_g_loc_load = []
    end

    gens_Pg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        P_gens[ n2s_gens_idx[ idx]] -
        P_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] :
        P_gens[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]

    gens_Qg_inj = [idx ∈ gens_nodes_with_loc_loads_idx ?
        Q_gens[ n2s_gens_idx[ idx]] -
        Q_g_loc_load[ n2s_gens_with_loc_load_idxs[ idx]] :
        Q_gens[ n2s_gens_idx[ idx]]
                  for idx in gens_nodes_idx ]


    Pg_inj_Qg_inj_Png_Qng =
        [gens_Pg_inj;
         gens_Qg_inj;
         P_non_gens;
         Q_non_gens]


    Pg_inj_Qg_inj_Png_Qng_Pll_Qll =
        [gens_Pg_inj;
         gens_Qg_inj;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    

    Pg_inj_Png_Qng =
        [gens_Pg_inj;
         P_non_gens;
         Q_non_gens]    
    
    return (;gens_Pg_inj,
            gens_Qg_inj,
            P_gens,
            Q_gens,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,
            Pg_inj_Png_Qng,
            Pg_inj_Qg_inj_Png_Qng,
            Pg_inj_Qg_inj_Png_Qng_Pll_Qll )
end



# page millano 264


function get_ph_pk_qh_qk_by_Ybr_cal(
    vh,
    θh,
    edges_Ybr_cal;
    edges_orientation,
    n2s_all_nodes_idx)
    
    uh = vh .* exp.( im * θh )

    from_idxs = [n2s_all_nodes_idx[ idx ]
                  for idx in
                      first.( edges_orientation ) ]
    
    to_idxs = [n2s_all_nodes_idx[ idx ]
                for idx in
                    second.( edges_orientation ) ]

    from_vh = vh[ from_idxs ]
    
    to_vh   = vh[ to_idxs ]

    from_θh = θh[ from_idxs ]
    
    to_θh   = θh[ to_idxs ]

    p_h = map(
        (y_π, v_h, v_k, θ_h, θ_k) ->
            real(y_π[1,1]) * v_h^2 +
            v_h * v_k * (real(y_π[1,2]) * cos(θ_h - θ_k) +
                         imag(y_π[1,2]) * sin(θ_h - θ_k)),
            edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )

    p_k = map(
        (y_π, v_h, v_k, θ_h, θ_k) ->
            real(y_π[2,2]) * v_k^2 +
            v_h * v_k * (real(y_π[1,2]) * cos(θ_h - θ_k) -
                         imag(y_π[1,2]) * sin(θ_h - θ_k)),
            edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )

    q_h = map(
        (y_π, v_h, v_k, θ_h, θ_k) ->
            -imag(y_π[1,1]) * v_h^2 +
            v_h * v_k * (real(y_π[1,2]) * sin(θ_h - θ_k) -
                         imag(y_π[1,2]) * cos(θ_h - θ_k)),
            edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )
    
    q_k = map(
        (y_π, v_h, v_k, θ_h, θ_k) ->
            -imag(y_π[2,2]) * v_k^2 -
            v_h * v_k * (real(y_π[1,2]) * sin(θ_h - θ_k) +
                         imag(y_π[1,2]) * cos(θ_h - θ_k)),
            edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )


    # q_k = map(
    #     (y_π, v_h, v_k, θ_h, θ_k) ->
    #         -imag(y_π[2,2]) * v_k^2 -
    #         v_h * v_k * (real(y_π[2,1]) * sin(θ_k - θ_h) -
    #                      imag(y_π[2,1]) * cos(θ_k - θ_h)),
    #      edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )
    

    # p_k = map(
    #     (y_π, v_h, v_k, θ_h, θ_k) ->
    #         real(y_π[2,2]) * v_k^2 +
    #         v_h * v_k * (real(y_π[2,1]) * cos(θ_k - θ_h) +
    #                      imag(y_π[2,1]) * sin(θ_k - θ_h)),
    #     edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )

    
    return (; p_h, p_k, q_h, q_k)

end


function get_ph_pk_qh_qk(
    vh,
    θh;
    Ybr_cal_and_edges_orientation,
    n2s_all_nodes_idx)
    
    (;edges_Ybr_cal,
     edges_orientation) =
         NamedTupleTools.select(
             Ybr_cal_and_edges_orientation,
             (:edges_Ybr_cal,
              :edges_orientation ))
    
    uh = vh .* exp.( im * θh )

    from_idxs = [n2s_all_nodes_idx[ idx ]
                  for idx in
                      first.( edges_orientation ) ]
    
    to_idxs = [n2s_all_nodes_idx[ idx ]
                for idx in
                    second.( edges_orientation ) ]

    from_vh = vh[ from_idxs ]
    
    to_vh   = vh[ to_idxs ]

    from_θh = θh[ from_idxs ]
    
    to_θh   = θh[ to_idxs ]

    p_h = map(
        (y_π, v_h, v_k, θ_h, θ_k) ->
            real(y_π[1,1]) * v_h^2 +
            v_h * v_k * (real(y_π[1,2]) * cos(θ_h - θ_k) +
                         imag(y_π[1,2]) * sin(θ_h - θ_k)),
            edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )

    p_k = map(
        (y_π, v_h, v_k, θ_h, θ_k) ->
            real(y_π[2,2]) * v_k^2 +
            v_h * v_k * (real(y_π[1,2]) * cos(θ_h - θ_k) -
                         imag(y_π[1,2]) * sin(θ_h - θ_k)),
            edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )

    q_h = map(
        (y_π, v_h, v_k, θ_h, θ_k) ->
            -imag(y_π[1,1]) * v_h^2 +
            v_h * v_k * (real(y_π[1,2]) * sin(θ_h - θ_k) -
                         imag(y_π[1,2]) * cos(θ_h - θ_k)),
            edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )
    
    q_k = map(
        (y_π, v_h, v_k, θ_h, θ_k) ->
            -imag(y_π[2,2]) * v_k^2 -
            v_h * v_k * (real(y_π[1,2]) * sin(θ_h - θ_k) +
                         imag(y_π[1,2]) * cos(θ_h - θ_k)),
            edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )


    # q_k = map(
    #     (y_π, v_h, v_k, θ_h, θ_k) ->
    #         -imag(y_π[2,2]) * v_k^2 -
    #         v_h * v_k * (real(y_π[2,1]) * sin(θ_k - θ_h) -
    #                      imag(y_π[2,1]) * cos(θ_k - θ_h)),
    #      edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )
    

    # p_k = map(
    #     (y_π, v_h, v_k, θ_h, θ_k) ->
    #         real(y_π[2,2]) * v_k^2 +
    #         v_h * v_k * (real(y_π[2,1]) * cos(θ_k - θ_h) +
    #                      imag(y_π[2,1]) * sin(θ_k - θ_h)),
    #     edges_Ybr_cal, from_vh, to_vh, from_θh, to_θh )

    
    return (; p_h, p_k, q_h, q_k)

end



function get_an_edge_p_net_q_net_injection(
    vh,
    θh,
    vk,
    θk,
    y_π)

    p_h = real(y_π[1,1]) * vh^2 +
            vh * vk * (real(y_π[1,2]) * cos(θh - θk) +
            imag(y_π[1,2]) * sin(θh - θk))

    p_k = real(y_π[2,2]) * vk^2 +
            vh * vk * (real(y_π[1,2]) * cos(θh - θk) -
                         imag(y_π[1,2]) * sin(θh - θk))

    q_h = -imag(y_π[1,1]) * vh^2 +
            vh * vk * (real(y_π[1,2]) * sin(θh - θk) -
            imag(y_π[1,2]) * cos(θh - θk))
    
    

    q_k = -imag(y_π[2,2]) * vk^2 -
            vh * vk * (real(y_π[1,2]) * sin(θh - θk) +
            imag(y_π[1,2]) * cos(θh - θk))
    
    
    return (; p_net = ( p_h - p_k),  q_net = (q_h- q_k))

end



function get_an_edge_injection_ph_pk_qh_qk(
    vh,
    θh,
    vk,
    θk,
    y_π)

    p_h = real(y_π[1,1]) * vh^2 +
            vh * vk * (real(y_π[1,2]) * cos(θh - θk) +
            imag(y_π[1,2]) * sin(θh - θk))

    p_k = real(y_π[2,2]) * vk^2 +
            vh * vk * (real(y_π[1,2]) * cos(θh - θk) -
                         imag(y_π[1,2]) * sin(θh - θk))

    q_h = -imag(y_π[1,1]) * vh^2 +
            vh * vk * (real(y_π[1,2]) * sin(θh - θk) -
                         imag(y_π[1,2]) * cos(θh - θk))
    

    q_k = -imag(y_π[2,2]) * vk^2 -
            vh * vk * (real(y_π[1,2]) * sin(θh - θk) +
                         imag(y_π[1,2]) * cos(θh - θk))
    
    return (; p_h, p_k, q_h, q_k)

end


function get_an_edge_injection_ph_pk_qh_qk(
    vh,
    θh,
    vk,
    θk;
    y_π)

    p_h = real(y_π[1,1]) * vh^2 +
            vh * vk * (real(y_π[1,2]) * cos(θh - θk) +
            imag(y_π[1,2]) * sin(θh - θk))

    p_k = real(y_π[2,2]) * vk^2 +
            vh * vk * (real(y_π[1,2]) * cos(θh - θk) -
                         imag(y_π[1,2]) * sin(θh - θk))

    q_h = -imag(y_π[1,1]) * vh^2 +
            vh * vk * (real(y_π[1,2]) * sin(θh - θk) -
                         imag(y_π[1,2]) * cos(θh - θk))
    

    q_k = -imag(y_π[2,2]) * vk^2 -
            vh * vk * (real(y_π[1,2]) * sin(θh - θk) +
                         imag(y_π[1,2]) * cos(θh - θk))
    
    return (; p_h, p_k, q_h, q_k)

end



function get_a_transformer_injection_ph_pk_qh_qk(
    vh,
    θh,
    vk,
    θk,
    y_π,
    t_ratio,
    ϕ;
    gFe=0,
    bmu=0)

    p_h = (gFe + real(y_π[1,1])/t_ratio^2)*v_h^2 +
            vh * vk * (real(y_π[1,2])*cos(θh-θk-ϕ) +
            imag(y_π[1,2])*sin(θh-θk))/t_ratio

    p_k = real(y_π[2,2]) * vk^2 +
            vh * vk * (real(y_π[1,2])*cos(θh-θk-ϕ) -
            imag(y_π[1,2])*sin(θh-θk))/t_ratio

    q_h = -(bmu + imag(y_π[1,1])/t_ratio^2)*v_h^2 +
            vh*vk*(real(y_π[1,2])*sin(θh-θk-ϕ) -
                         imag(y_π[1,2])*cos(θh-θk-ϕ))/t_ratio
    
    q_k = -imag(y_π[2,2])*v_k^2 -
            vh*vk*(real(y_π[1,2])*sin(θh-θk-ϕ) +
                         imag(y_π[1,2])*cos(θh-θk-ϕ))/t_ratio
    
    return (; p_h, p_k, q_h, q_k)

end


# ---------------------------------------------------
# ---------------------------------------------------

function get_loss_particpation_by_gens_rating(
    gens_Sn, Pg)

    participating_gens =
        [ abs(P) < 0.0000001 ? 0 : 1  for P in Pg]

    gens_loss_participation_factor =
        (gens_Sn .* participating_gens) ./
        sum(gens_Sn .* participating_gens)

    return (; participating_gens,
            gens_loss_participation_factor)

end


function get_loss_particpation_by_gens_loading( Pg)

    participating_gens =
        [ abs(P) < 0.0000001 ? 0 : 1  for P in Pg]

    gens_loss_participation_factor =
        (Pg .* participating_gens) ./
        sum(Pg .* participating_gens)

    return (; participating_gens,
            gens_loss_participation_factor)

end

# ---------------------------------------------------


function get_gens_active_power_particpation_by_rating(
    gens_Sn, Pg)

    participating_gens =
        [ abs(P) < 0.0000001 ? 0 : 1  for P in Pg]

    gens_active_power_particpation_factor =
        (gens_Sn .* participating_gens) ./
        sum(gens_Sn .* participating_gens)

    return (; participating_gens,
            gens_active_power_particpation_factor)

end


function get_gens_active_power_particpation_by_loading( Pg )

    participating_gens =
        [ abs(P) < 0.0000001 ? 0 : 1  for P in Pg]

    gens_active_power_particpation_factor =
        (Pg .* participating_gens) ./
        sum(Pg .* participating_gens)

    return (; participating_gens,
            gens_active_power_particpation_factor)

end


function get_gens_reactive_power_particpation_by_rating(
    gens_Sn, Qg)

    participating_gens =
        [ abs(Q) < 0.0000001 ? 0 : 1  for Q in Qg]

     gens_reactive_power_particpation_factor =
        (gens_Sn .* participating_gens) ./
        sum(gens_Sn .* participating_gens)

    return (; participating_gens,
            gens_reactive_power_particpation_factor)

end


function get_gens_reactive_power_particpation_by_loading( Qg)

    participating_gens =
        [ abs(Q) < 0.0000001 ? 0 : 1  for Q in Qg]

    gens_reactive_power_particpation_factor =
        (Qg .* participating_gens) ./
        sum(Qg .* participating_gens)

    return (; participating_gens,
            gens_reactive_power_particpation_factor)

end


# ---------------------------------------------------

function get_losses_per_line(
    vh_θh,
    edges_r_x_b_ratio_angle;    
    edges_type,
    edges_orientation,
    n2s_all_nodes_idx,
    dyn_pf_flat_vh_flat_θh_Idx,
    edges_r_x_b_ratio_angle_idx)

    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    
    (; r_Idxs, x_Idxs, b_Idxs,
     ratio_Idxs, angle_Idxs) =
         NamedTupleTools.select(
             edges_r_x_b_ratio_angle_idx,
             (:r_Idxs, :x_Idxs, :b_Idxs,
              :ratio_Idxs, :angle_Idxs))


    edges_r = edges_r_x_b_ratio_angle[r_Idxs]
    edges_x = edges_r_x_b_ratio_angle[x_Idxs]
    edges_b = edges_r_x_b_ratio_angle[b_Idxs]
    edges_ratio = edges_r_x_b_ratio_angle[ratio_Idxs]
    edges_angle = edges_r_x_b_ratio_angle[angle_Idxs]
    
    vh = vh_θh[vh_Idxs]

    θh = vh_θh[θh_Idxs]
    
    from_idxs = [n2s_all_nodes_idx[ idx ]
                  for idx in
                      first.(edges_orientation) ]
    
    to_idxs = [n2s_all_nodes_idx[ idx ]
                for idx in
                    second.(edges_orientation) ]
    
    uh = vh .* exp.( im * θh )

    from_uh = uh[ from_idxs ]
    
    to_uh = uh[ to_idxs ]
    
    return (abs.(inv.( edges_r + im * edges_x ) .* (
            from_uh - to_uh))).^2 .* edges_r
end


function get_sum_line_losses(
    vh_θh,
    edges_r_x_b_ratio_angle;    
    edges_type,
    edges_orientation,
    n2s_all_nodes_idx,
    dyn_pf_flat_vh_flat_θh_Idx,
    edges_r_x_b_ratio_angle_idx)

    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    
    (; r_Idxs, x_Idxs, b_Idxs,
     ratio_Idxs, angle_Idxs) =
         NamedTupleTools.select(
             edges_r_x_b_ratio_angle_idx,
             (:r_Idxs, :x_Idxs, :b_Idxs,
              :ratio_Idxs, :angle_Idxs))


    edges_r = edges_r_x_b_ratio_angle[r_Idxs]
    edges_x = edges_r_x_b_ratio_angle[x_Idxs]
    edges_b = edges_r_x_b_ratio_angle[b_Idxs]
    edges_ratio = edges_r_x_b_ratio_angle[ratio_Idxs]
    edges_angle = edges_r_x_b_ratio_angle[angle_Idxs]
    
    vh = vh_θh[vh_Idxs]

    θh = vh_θh[θh_Idxs]
    
    from_idxs = [n2s_all_nodes_idx[ idx ]
                  for idx in
                      first.(edges_orientation) ]
    
    to_idxs = [n2s_all_nodes_idx[ idx ]
                for idx in
                    second.(edges_orientation) ]
    
    uh = vh .* exp.( im * θh )

    from_uh = uh[ from_idxs ]
    
    to_uh = uh[ to_idxs ]
    
    return sum((abs.(inv.( edges_r + im * edges_x ) .* (
            from_uh - to_uh))).^2 .* edges_r)
end



function get_losses_per_node_by_flattened_Ynet(
    vh_θh, Ynet_real_imag_flattend;
    dyn_pf_flat_vh_flat_θh_Idx,
    nodes_idx_with_adjacent_nodes_idx,
    Ynet_rows_Idxs_in_flattend,
    Ynet_real_imag_Idxs_in_flattend,
    n2s_all_nodes_idx,
    all_nodes_idx )

    (Ynet_real_Idxs,
     Ynet_imag_Idxs ) =
         Ynet_real_imag_Idxs_in_flattend


    Ynet_real =
        Ynet_real_imag_flattend[ Ynet_real_Idxs ]

    Ynet_imag =
        Ynet_real_imag_flattend[ Ynet_imag_Idxs ]


    Ynet = [Ynet_real[idx] + im * Ynet_imag[idx]  
              for idx in Ynet_rows_Idxs_in_flattend ]
    

    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    
    vh = vh_θh[vh_Idxs]

    θh = vh_θh[θh_Idxs]
    
    return  [ (vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ]))
                        for nth_idx in all_nodes_idx ]     
end


function get_total_P_network_loss_by_flattened_Ynet(
    vh_θh,
    Ynet_real_imag_flattend;
    dyn_pf_flat_vh_flat_θh_Idx,
    nodes_idx_with_adjacent_nodes_idx,
    Ynet_rows_Idxs_in_flattend,
    Ynet_real_imag_Idxs_in_flattend,
    n2s_all_nodes_idx,
    all_nodes_idx )

    (Ynet_real_Idxs,
     Ynet_imag_Idxs ) =
         Ynet_real_imag_Idxs_in_flattend

    Ynet_real =
        Ynet_real_imag_flattend[Ynet_real_Idxs]

    Ynet_imag =
        Ynet_real_imag_flattend[Ynet_imag_Idxs]

    Ynet = [Ynet_real[idx] + im * Ynet_imag[idx]  
              for idx in Ynet_rows_Idxs_in_flattend ]
    
    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    
    vh = vh_θh[vh_Idxs]

    θh = vh_θh[θh_Idxs]
    
    return sum([(vh[ n2s_all_nodes_idx[nth_idx]] * sum(
        [vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            cos(
                θh[n2s_all_nodes_idx[nth_idx]] -
                   θh[ n2s_all_nodes_idx[idx]] -
                   angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[
                               nth_idx]],
                         nodes_idx_with_adjacent_nodes_idx[
                           n2s_all_nodes_idx[ nth_idx]])]))
                for nth_idx in all_nodes_idx ] )    
end

function get_total_Q_network_loss_by_flattened_Ynet(
    vh_θh, Ynet_real_imag_flattend;
    dyn_pf_flat_vh_flat_θh_Idx,
    nodes_idx_with_adjacent_nodes_idx,
    Ynet_rows_Idxs_in_flattend,
    Ynet_real_imag_Idxs_in_flattend,
    n2s_all_nodes_idx,
    all_nodes_idx )

    (Ynet_real_Idxs,
     Ynet_imag_Idxs ) =
         Ynet_real_imag_Idxs_in_flattend

    Ynet_real =
        Ynet_real_imag_flattend[Ynet_real_Idxs]

    Ynet_imag =
        Ynet_real_imag_flattend[Ynet_imag_Idxs]

    Ynet = [Ynet_real[idx] + im * Ynet_imag[idx]  
              for idx in Ynet_rows_Idxs_in_flattend ]
    
    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    
    vh = vh_θh[vh_Idxs]

    θh = vh_θh[θh_Idxs]
    
    return  sum([
        (vh[ n2s_all_nodes_idx[ nth_idx ] ] *
            sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin( θh[ n2s_all_nodes_idx[nth_idx]] -
            θh[n2s_all_nodes_idx[idx]] - angle(ynj))
                   for (ynj, idx) in
                       zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                           nodes_idx_with_adjacent_nodes_idx[
                               n2s_all_nodes_idx[nth_idx]])]))
        for nth_idx in all_nodes_idx ])
end






function get_total_P_network_loss_by_sta_pf_PQ_para(
    sta_pf_PQ_para,
    loc_load_exist )

    
    (Pg,
     Png,
     Pll) =
        NamedTupleTools.select(
            sta_pf_PQ_para,
            (:P_gens,
             :P_non_gens,
             :P_g_loc_load ) )
    
    return loc_load_exist == true ?
        sum(Pg) - sum([ Png; Pll]) : sum(Pg) - sum(Png ) 
end



function get_total_Q_network_loss_by_sta_pf_PQ_para(
    sta_pf_PQ_para,
    loc_load_exist )

    
    (
     Qg,
     Qng, 
     Qll) =
        NamedTupleTools.select(
            sta_pf_PQ_para,
            (
             :Q_gens,
             :Q_non_gens,
             :Q_g_loc_load ) )
    
    return loc_load_exist == true ?
        sum(Qg) - sum([ Qng; Qll]) : sum(Qg ) - sum(Qng) 
end




function get_total_P_network_loss(
    vh_θh,
    Ynet;
    dyn_pf_flat_vh_flat_θh_Idx,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx,
    all_nodes_idx )

    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    
    vh = vh_θh[vh_Idxs]

    θh = vh_θh[θh_Idxs]
    
    return sum([ (
        vh[ n2s_all_nodes_idx[nth_idx]] * sum(
            [ vh[n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(
                    θh[n2s_all_nodes_idx[nth_idx]] -
                    θh[ n2s_all_nodes_idx[idx]] -
                    angle(ynj) )
              for (ynj, idx) in
                  zip(Ynet[n2s_all_nodes_idx[
                      nth_idx]],
                      nodes_idx_with_adjacent_nodes_idx[
                          n2s_all_nodes_idx[
                          nth_idx]]) ]))
                        for nth_idx in all_nodes_idx ] )    
end


function get_total_Q_network_loss(
    vh_θh, Ynet;
    dyn_pf_flat_vh_flat_θh_Idx,
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx,
    all_nodes_idx )

    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))
    
    vh = vh_θh[vh_Idxs]

    θh = vh_θh[θh_Idxs]
    
    return  sum([
        (vh[ n2s_all_nodes_idx[ nth_idx ] ] *
            sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin( θh[ n2s_all_nodes_idx[nth_idx]] -
            θh[n2s_all_nodes_idx[idx]] - angle(ynj))
                   for (ynj, idx) in
                       zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                           nodes_idx_with_adjacent_nodes_idx[
                               n2s_all_nodes_idx[nth_idx]])]))
        for nth_idx in all_nodes_idx ])
end


function get_total_P_network_loss(
    vh,
    θh,
    Ynet;
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx,
    all_nodes_idx )
    
    return sum([
        (vh[ n2s_all_nodes_idx[nth_idx]] * sum(
            [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(
                    θh[n2s_all_nodes_idx[nth_idx]] -
                        θh[ n2s_all_nodes_idx[idx]] -
                        angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[
                               nth_idx]],
                          nodes_idx_with_adjacent_nodes_idx[
                         n2s_all_nodes_idx[nth_idx]])]))
                        for nth_idx in all_nodes_idx])    
end


function get_total_Q_network_loss(
    vh,
    θh,
    Ynet;
    nodes_idx_with_adjacent_nodes_idx,
    n2s_all_nodes_idx,
    all_nodes_idx )
    
    return  sum([
        (vh[ n2s_all_nodes_idx[ nth_idx ] ] *
            sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin( θh[ n2s_all_nodes_idx[nth_idx]] -
            θh[n2s_all_nodes_idx[idx]] - angle(ynj))
                   for (ynj, idx) in
                       zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                           nodes_idx_with_adjacent_nodes_idx[
                               n2s_all_nodes_idx[nth_idx]])]))
        for nth_idx in all_nodes_idx ])
end



function get_a_line_loss_by_vh_θh(
    edge_r,
    edge_x,
    f_vh,
    f_θh,
    t_vh,
    t_θh )

    # |z^-1 * (f_uh - t_uh)|^2 * r
    
    return (abs( inv(edge_r + im * edge_x) * (
        (f_vh * exp(im * f_θh)) - (
            t_vh * exp(im * t_θh))) ))^2 * edge_r

end


function get_a_line_loss_by_ur_ui(    
    edge_r,
    edge_x,
    f_uh_real,
    f_uh_imag,
    t_uh_real,
    t_uh_imag )

    # |z^-1 * (f_uh - t_uh)|^2 * r
    
    return abs(inv(edge_r + im * edge_x) * (
        (f_uh_real + im * f_uh_imag) - (
            t_uh_real + im * t_uh_imag)))^2 * edge_r


end


function get_a_line_loss_by_uh(
    edge_r,
    edge_x,
    f_uh,
    t_uh )

    # |z^-1 * (f_uh - t_uh)|^2 * r
    
    return (abs( inv(edge_r + im * edge_x) * (
        f_uh - t_uh) ))^2 * edge_r

end



function get_lines_losses_cal(
    vh,
    θh;
    Ybr_cal_and_edges_orientation,
    n2s_all_nodes_idx)
    
    (;edges_Ybr_cal,
     edges_orientation) =
         NamedTupleTools.select(
             Ybr_cal_and_edges_orientation,
             (:edges_Ybr_cal,
              :edges_orientation ))

    y_lines = -1 .* [ y_π[1,2] for y_π in edges_Ybr_cal]

    r_lines = real.( inv.( y_lines ))
    
    uh = vh .* exp.( im * θh )

    from_idxs = [n2s_all_nodes_idx[ idx ]
                  for idx in
                      first.( edges_orientation ) ]
    
    to_idxs = [n2s_all_nodes_idx[ idx ]
                for idx in
                    second.( edges_orientation ) ]

    from_uh = uh[ from_idxs ]
    
    to_uh = uh[ to_idxs ]

    current_in_lines =
        map((y, f_u, t_u) -> y * (f_u - t_u),
            y_lines, from_uh, to_uh )

    loss_in_lines =
        map((r, i_line) -> r * (abs(i_line))^2,
            r_lines, current_in_lines )

    total_line_loss = sum(loss_in_lines)
    
    return (;total_line_loss,
            loss_in_lines,
            current_in_lines,
            y_lines,
            r_lines,
            from_idxs,
            to_idxs,
            from_uh,
            to_uh)
end



function get_lines_losses_cal(
    vh,
    θh,
    edges_r,
    edges_x,
    edges_b,
    edges_ratio,
    edges_angle;    
    edges_type,
    edges_orientation,
    n2s_all_nodes_idx)

    from_idxs = [n2s_all_nodes_idx[ idx ]
                  for idx in
                      first.(edges_orientation) ]
    
    to_idxs = [n2s_all_nodes_idx[ idx ]
                for idx in
                    second.(edges_orientation) ]

    
    uh = vh .* exp.( im * θh )

    from_uh = uh[ from_idxs ]
    
    to_uh = uh[ to_idxs ]

    z_lines = edges_r + im * edges_x
    
    r_lines = edges_r

    current_in_lines =
        inv.(z_lines) .* (from_uh - to_uh)

    loss_in_lines =
        (abs.(current_in_lines)).^2 .* r_lines
    
    total_line_loss = sum(loss_in_lines)
    
    return (;total_line_loss,
            loss_in_lines,
            current_in_lines,
            z_lines,
            r_lines,
            from_idxs,
            to_idxs,
            from_uh,
            to_uh)
end



function get_lines_total_losses(
    vh,
    θh;
    Ybr_cal_and_edges_orientation,
    n2s_all_nodes_idx)
    
    (;edges_Ybr_cal,
     edges_orientation) =
         NamedTupleTools.select(
             Ybr_cal_and_edges_orientation,
             (:edges_Ybr_cal,
              :edges_orientation ))

    y_lines = -1 .* [ y_π[1,2] for y_π in edges_Ybr_cal]

    r_lines = real.( inv.( y_lines ))
    
    uh = vh .* exp.( im * θh )

    from_idxs = [n2s_all_nodes_idx[ idx ]
                  for idx in
                      first.( edges_orientation ) ]
    
    to_idxs = [n2s_all_nodes_idx[ idx ]
                for idx in
                    second.( edges_orientation ) ]

    from_uh = uh[ from_idxs ]
    
    to_uh = uh[ to_idxs ]

    current_in_lines =
        map((y, f_u, t_u) -> y * (f_u - t_u),
            y_lines, from_uh, to_uh )

    return sum(map((r, i_line) -> r * (abs(i_line))^2,
            r_lines, current_in_lines ))

end

function get_lines_total_losses(
    vh,
    θh,
    edges_r,
    edges_x,
    edges_b,
    edges_ratio,
    edges_angle;    
    edges_type,
    edges_orientation,
    n2s_all_nodes_idx)

    from_idxs = [n2s_all_nodes_idx[ idx ]
                  for idx in
                      first.(edges_orientation) ]
    
    to_idxs = [n2s_all_nodes_idx[ idx ]
                for idx in
                    second.(edges_orientation) ]

    
    uh = vh .* exp.( im * θh )

    from_uh = uh[ from_idxs ]
    
    to_uh = uh[ to_idxs ]

    z_lines = edges_r + im * edges_x
    
    r_lines = edges_r

    current_in_lines =
        inv.(z_lines) .* (from_uh - to_uh)

    return sum((abs.(current_in_lines)).^2 .* r_lines)

end

# ------------------------------------------------------
# comment
# ------------------------------------------------------
