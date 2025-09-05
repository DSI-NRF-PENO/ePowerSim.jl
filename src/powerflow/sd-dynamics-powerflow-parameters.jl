# (c) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.4:  pg 135 - 136, 6.121 - 6.123



########################################################
# ------------------------------------------------------
#  Powerflow parameters and indices
# ------------------------------------------------------
########################################################



function loc_load_exist_bool(netd)

    return get_gens_nodes_with_loc_loads_idx(
            netd.nodes ) == [] ? false : true
end



#-------------------------------------------------------
# functions with netd as an argument
#-------------------------------------------------------

function get_Ynet_no_shunt( netd )

    #----------------------------------------------
    #----------------------------------------------

    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt( netd.nodes )
    
    nodes_incident_edges =
        get_nodes_incident_edges( netd.edges )
    
    edges_Ybr_cal =
        [ Symbol(typeof(a_branch )) == :PiModelLine ?
        calc_branch_Ybr(a_branch.param_values[1:3]...,1,1) :
        calc_branch_Ybr(a_branch.param_values[1:4]..., 1 )
          for a_branch in collect(values(netd.edges)) ]

    edges_orientation =
        [ a_branch.orientation
          for a_branch in collect(values(netd.edges)) ]

    edges_orientation_and_edges_Ybr_cal =
        [ (orient, a_Ybr)
          for (orient, a_Ybr) in
              zip( edges_orientation, edges_Ybr_cal)] 

    #-----------------------------------------
  
    nodes_incident_edges_and_orientation =
        [ edges_orientation[ a_node_incident_edges]
          for a_node_incident_edges in
              nodes_incident_edges ]

    nodes_node_idx_and_incident_edges_other_node_idx =
        [[ [[idx], [
            idx == orient[1] ? orient[2] : orient[1]
                     for orient in
                         a_node_edges_other_nodes_idx]]...;]
         for (idx, a_node_edges_other_nodes_idx ) in
             enumerate(
                 nodes_incident_edges_and_orientation)]    
          
    #------------------------------------------
    
    nodes_incident_edges_orientation_and_Ybr_cal =
        [ edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
          for a_node_incident_edges in
              nodes_incident_edges]
    
    #----------------------------------------------
    #----------------------------------------------    

    Ynet_no_shunt = [ ]

    for (node_idx, orientations_and_edges_Ybr ) in
        enumerate(
            nodes_incident_edges_orientation_and_Ybr_cal)
        
        Ynet_no_shunt_nth_row = []
        
        for k in 0:length( orientations_and_edges_Ybr )
            if  k == 0
                
                push!(
                    Ynet_no_shunt_nth_row,
                    sum([ node_idx == first(
                        first(orient_and_Ybr)) ?
                        last(orient_and_Ybr)[1] :
                        last(orient_and_Ybr)[4]
                          for orient_and_Ybr in
                              orientations_and_edges_Ybr]))

            else
                if  node_idx == first(first(
                    orientations_and_edges_Ybr[k] ))
                    
                    push!(
                        Ynet_no_shunt_nth_row,
                        last(orientations_and_edges_Ybr[
                            k])[3])
                    
                else
                    
                    push!(
                        Ynet_no_shunt_nth_row,
                        last(orientations_and_edges_Ybr[
                            k])[2])
                    
                end
            end
        end
        
        push!(Ynet_no_shunt, Ynet_no_shunt_nth_row) 
    end

    return Ynet_no_shunt

end


function get_Ynet( netd )

    #-------------------------------------------
    #-------------------------------------------

    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt( netd.nodes )

    Ynet_no_shunt =
        get_Ynet_no_shunt( netd )
    
    #-------------------------------------------
    #-------------------------------------------


    Ynet = []

    for (shunt_idx, node_k_idx_and_shunt) in
        enumerate( nodes_idx_and_Yshunt )
        
        Ynet_nth_row = []
        
        for (idx, Ynet_k_element ) in
            enumerate( Ynet_no_shunt[shunt_idx] )
            
            if idx == 1
                
                push!(
                    Ynet_nth_row,
                    last(node_k_idx_and_shunt) +
                        Ynet_k_element)
            else
                push!(Ynet_nth_row, Ynet_k_element)
            end
        end
        push!(Ynet, Ynet_nth_row)
    end

    return Ynet
end


function get_Ybus_no_shunt( netd )

    #-----------------------------------------

    nodes_incident_edges =
        get_nodes_incident_edges( netd.edges )
    
    edges_orientation =
        [ a_branch.orientation
          for a_branch in collect(values(netd.edges)) ]
    
    nodes_incident_edges_and_orientation =
        [ edges_orientation[ a_node_incident_edges ]
          for a_node_incident_edges in
              nodes_incident_edges ]

    nodes_node_idx_and_incident_edges_other_node_idx =
        [[ [[idx], [
            idx == orient[1] ? orient[2] : orient[1]
            for orient in
                a_node_edges_other_nodes_idx ] ]...; ]
         for (idx, a_node_edges_other_nodes_idx ) in
             enumerate(
                 nodes_incident_edges_and_orientation)]    
          
    #------------------------------------------
       
    Ynet_no_shunt = get_Ynet_no_shunt( netd )

    Ybus_row_size = Ybus_col_size = length(Ynet_no_shunt )

    return reshape(
        vcat(
            [[ a_col in nth_row_idx ?
                Ybus_nth_row[
                    findall(x-> x==a_col,
                            nth_row_idx)[1]] : 0.0
               for a_col in 1:Ybus_col_size]
             for (nth_row_idx, Ybus_nth_row) in
                 zip(nodes_node_idx_and_incident_edges_other_node_idx,
                     Ynet_no_shunt )]...;),
        Ybus_row_size, Ybus_col_size )    
end


function get_Ybus( netd  )

    #-----------------------------------------


    nodes_incident_edges =
        get_nodes_incident_edges( netd.edges )
    
    edges_orientation =
        [ a_branch.orientation
          for a_branch in collect(values(netd.edges)) ]
    
    nodes_incident_edges_and_orientation =
        [ edges_orientation[ a_node_incident_edges ]
          for a_node_incident_edges in nodes_incident_edges ]

    nodes_node_idx_and_incident_edges_other_node_idx =
        [[ [[idx], [ idx == orient[1] ? orient[2] : orient[1]
                     for orient in a_node_edges_other_nodes_idx]]...;]
         for (idx, a_node_edges_other_nodes_idx ) in
             enumerate( nodes_incident_edges_and_orientation )  ]    
     
     
    #------------------------------------------

    Ynet = get_Ynet( netd )

    Ybus_row_size = Ybus_col_size = length( Ynet )

    return reshape(
        vcat([[ a_col in nth_row_idx ?
            Ybus_nth_row[ findall(x-> x==a_col, nth_row_idx)[1] ]  :
            0.0
                for a_col in 1:Ybus_col_size]
              for (nth_row_idx, Ybus_nth_row) in
                  zip( nodes_node_idx_and_incident_edges_other_node_idx,
                       Ynet )]...;),
        Ybus_row_size, Ybus_col_size )
    
end

function get_sparseYbus( netd  )

    #-----------------------------------------


    nodes_incident_edges = get_nodes_incident_edges( netd.edges )
    
    edges_orientation =
        [ a_branch.orientation
          for a_branch in collect(values(netd.edges)) ]
    
    nodes_incident_edges_and_orientation =
        [ edges_orientation[ a_node_incident_edges ]
          for a_node_incident_edges in nodes_incident_edges ]

    nodes_node_idx_and_incident_edges_other_node_idx =
        [[ [[idx], [ idx == orient[1] ? orient[2] : orient[1]
                     for orient in a_node_edges_other_nodes_idx]]...;]
         for (idx, a_node_edges_other_nodes_idx ) in
             enumerate( nodes_incident_edges_and_orientation )  ]    
     
     
    #------------------------------------------

    Ynet = get_Ynet( netd )

    Ybus_row_size = Ybus_col_size = length( Ynet )

    return sparse(
        reshape(vcat([[ a_col in nth_row_idx ?
            Ybus_nth_row[ findall(x-> x==a_col,
                                  nth_row_idx)[1] ] :
            0.0
                        for a_col in 1:Ybus_col_size]
                      for  (nth_row_idx, Ybus_nth_row) in
                          zip( nodes_node_idx_and_incident_edges_other_node_idx,
                               Ynet )]...;),
                Ybus_row_size, Ybus_col_size ))
    
end



function get_net_Ynet_etc_parameters( netd )

    # dynamics_case = new_case_IEEE_5_Bus_dynamic_plant_SM_v6_P

    # netd  = NetworkData( dynamics_case()... )

    # ----------------------------------------------------
    # ----------------------------------------------------

    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt(netd.nodes )
    

    nodes_incident_edges =
        get_nodes_incident_edges( netd.edges )

    
    edges_Ybr_cal =
        [ Symbol(typeof(a_branch )) == :PiModelLine ?
        calc_branch_Ybr(
            a_branch.param_values[1:3]..., 1, 1) :
        calc_branch_Ybr(a_branch.param_values[1:4]..., 1 )
          for a_branch in collect(values(netd.edges)) ]

    edges_orientation =
        [ a_branch.orientation
          for a_branch in collect(values(netd.edges)) ]

    # edges_Ybr_cal_view = [view(edges_Ybr_cal, idx:idx)  for idx in 1:length(edges_Ybr_cal) ]


    # edges_orientation_and_edges_Ybr_cal_view = [ (orient, a_view) for (orient, a_view) in zip( edges_orientation, edges_Ybr_cal_view )]

    edges_orientation_and_edges_Ybr_cal =
        [ (orient, a_Ybr)
          for (orient, a_Ybr) in
              zip( edges_orientation, edges_Ybr_cal )] 

    #-----------------------------------------
  
    nodes_incident_edges_and_orientation =
        [ edges_orientation[ a_node_incident_edges ]
          for a_node_incident_edges in
              nodes_incident_edges ]

    nodes_node_idx_and_incident_edges_other_node_idx =
        [[ [[idx], [idx == orient[1] ?
        orient[2] : orient[1]
                    for orient in
                        a_node_edges_other_nodes_idx]]...;]
         for (idx, a_node_edges_other_nodes_idx ) in
             enumerate(
                 nodes_incident_edges_and_orientation)]
    
          
    #------------------------------------------
    
    nodes_incident_edges_orientation_and_Ybr_cal =
        [ edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
          for a_node_incident_edges in
              nodes_incident_edges]

    Ynet_no_shunt = [
        [ k == 0 ?
            sum([node_idx ==
            first(first(orient_and_Ybr)) ?
        last(orient_and_Ybr)[1] :
        last(orient_and_Ybr)[4]
                   for orient_and_Ybr in
                       orientations_and_edges_Ybr]) :
                           node_idx == first(first( orientations_and_edges_Ybr[k] )) ?
                           last(orientations_and_edges_Ybr[k])[3] :
                           last(orientations_and_edges_Ybr[k])[2] for k in
                               0:length( orientations_and_edges_Ybr ) ]
                      for (node_idx, orientations_and_edges_Ybr ) in
                          enumerate(
                              nodes_incident_edges_orientation_and_Ybr_cal ) ]

    Ynet = [
        [  idx == 1 ?
            last(node_k_idx_and_shunt) + Ynet_k_element :
            Ynet_k_element for (idx, Ynet_k_element) in
                enumerate(Ynet_no_shunt[shunt_idx])]
        for (shunt_idx, node_k_idx_and_shunt) in
            enumerate( nodes_idx_and_Yshunt ) ]

    Ybus_row_size = Ybus_col_size = length(Ynet )

    Ybus = reshape(
        vcat([[ a_col in nth_row_idx ?
            Ybus_nth_row[
                findall(x-> x==a_col, nth_row_idx)[1]] :
            0.0
                for a_col in 1:Ybus_col_size]

              for (nth_row_idx, Ybus_nth_row) in
                  zip( nodes_node_idx_and_incident_edges_other_node_idx,
                       Ynet )]...;),
        Ybus_row_size, Ybus_col_size )

    Ybus_no_shunt = reshape(
        vcat([[ a_col in nth_row_idx ?
            Ybus_nth_row[
                findall(x-> x==a_col, nth_row_idx)[1]] :
            0.0
                for a_col in 1:Ybus_col_size]
              for (nth_row_idx, Ybus_nth_row) in
                  zip( nodes_node_idx_and_incident_edges_other_node_idx,
                       Ynet_no_shunt)]...;),
        Ybus_row_size, Ybus_col_size )


    return (; nodes_idx_and_Yshunt,
            nodes_incident_edges,edges_Ybr_cal,
            edges_orientation,edges_orientation_and_edges_Ybr_cal,
            nodes_incident_edges_and_orientation,
            nodes_node_idx_and_incident_edges_other_node_idx,
            nodes_incident_edges_orientation_and_Ybr_cal,
            Ynet_no_shunt, Ynet )
    
    
    # show(stdout, "text/plain", Ybus )    
    # show(stdout, "text/plain", Ybus_no_shunt )
      
           
end

#-------------------------------------------------------


function get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx(
    netd )

    # ----------------------------------------------------
    # ----------------------------------------------------

    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt( netd.nodes )
    

    nodes_incident_edges =
        get_nodes_incident_edges( netd.edges )

    
    edges_Ybr_cal = [
        Symbol(typeof(a_branch )) == :PiModelLine ?
            calc_branch_Ybr(
                a_branch.param_values[1:3]..., 1, 1 ) :
                    calc_branch_Ybr(
                        a_branch.param_values[1:4]..., 1 )
        for a_branch in
            collect(values(netd.edges)) ]

    edges_orientation = [
        a_branch.orientation
        for a_branch in collect(values(netd.edges)) ]

    edges_orientation_and_edges_Ybr_cal = [
        (orient, a_Ybr)
        for (orient, a_Ybr) in
            zip( edges_orientation, edges_Ybr_cal )] 

    #-----------------------------------------
  
    nodes_incident_edges_and_orientation = [
        edges_orientation[ a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges ]

    nodes_node_idx_and_incident_edges_other_node_idx = [
        [ [[idx], [ idx == orient[1] ?
            orient[2] : orient[1]
                    for orient in
                        a_node_edges_other_nodes_idx]]...;]
        for (idx, a_node_edges_other_nodes_idx ) in
            enumerate(
                nodes_incident_edges_and_orientation)]    
          
    #------------------------------------------
    
    nodes_incident_edges_orientation_and_Ybr_cal = [
        edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges]

    Ynet_no_shunt = [
        [ k == 0 ? sum(
            [ node_idx == first( first(orient_and_Ybr) ) ?
        last(orient_and_Ybr)[1] :
        last(orient_and_Ybr)[4]
                          for orient_and_Ybr in
                              orientations_and_edges_Ybr]) :
                                  node_idx == first(first(
                                      orientations_and_edges_Ybr[k] )) ?
                                  last(orientations_and_edges_Ybr[k])[3] :
                                  last(orientations_and_edges_Ybr[k])[2]
          for k in 0:length( orientations_and_edges_Ybr ) ]
                      for (node_idx, orientations_and_edges_Ybr ) in
                          enumerate(
                              nodes_incident_edges_orientation_and_Ybr_cal ) ]


    Ynet = Vector{ComplexF64}[
        [  idx == 1 ?
            last(node_k_idx_and_shunt) + Ynet_k_element :
            Ynet_k_element
           for (idx, Ynet_k_element ) in
               enumerate( Ynet_no_shunt[shunt_idx] ) ]
        for (shunt_idx, node_k_idx_and_shunt) in
            enumerate( nodes_idx_and_Yshunt) ]

    return (; edges_Ybr_cal,
            Ynet, edges_orientation,
            nodes_node_idx_and_incident_edges_other_node_idx )
      
           
end


function get_Ynet_and_nodes_idx_wt_adjacent_nodes_idx(
    netd )

    # -----------------------------------------------

    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt(
            netd.nodes )
    
    nodes_incident_edges =
        get_nodes_incident_edges(
            netd.edges )

    
    edges_Ybr_cal =
        [ Symbol(typeof(a_branch )) == :PiModelLine ?
        calc_branch_Ybr(
            a_branch.param_values[1:3]..., 1, 1 ) :
                calc_branch_Ybr(
                    a_branch.param_values[1:4]..., 1 )
          for a_branch in
              collect(values(netd.edges)) ]

    edges_orientation =
        [ a_branch.orientation
          for a_branch in
              collect(values(netd.edges)) ]

    edges_orientation_and_edges_Ybr_cal =
        [ (orient, a_Ybr)
          for (orient, a_Ybr) in
              zip( edges_orientation,
                   edges_Ybr_cal )] 

    #-----------------------------------------
 
    nodes_incident_edges_and_orientation =
        [ edges_orientation[
            a_node_incident_edges ]
          for a_node_incident_edges in
              nodes_incident_edges ]


    # nodes_idx_with_adjacent_nodes_idx
    # nodes_node_idx_and_incident_edges_other_node_idx
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx]]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            enumerate(
                nodes_incident_edges_and_orientation)]
          
    #------------------------------------------

    nodes_incident_edges_orientation_and_Ybr_cal =[
        edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges]

    Ynet_no_shunt = [
        [ k == 0 ?
            sum( [ node_idx == first( first(
                orient_and_Ybr) ) ?
                    last(orient_and_Ybr)[1] :
                    last(orient_and_Ybr)[4]
                   for orient_and_Ybr in
                       orientations_and_edges_Ybr]) :
                           node_idx == first(first(
                               orientations_and_edges_Ybr[k] )) ?
                                   last(orientations_and_edges_Ybr[k])[3] :
                                   last(orientations_and_edges_Ybr[k])[2]
          for k in 0:length(
              orientations_and_edges_Ybr ) ]
        for (node_idx,orientations_and_edges_Ybr) in 
            enumerate(
                nodes_incident_edges_orientation_and_Ybr_cal ) ]


    Ynet =
        Vector{ComplexF64}[
            [  idx == 1 ?
                last(node_k_idx_and_shunt) +
                Ynet_k_element :
                Ynet_k_element
               for (idx, Ynet_k_element) in
                   enumerate(
                       Ynet_no_shunt[shunt_idx] ) ]
            for (shunt_idx, node_k_idx_and_shunt) in
                enumerate( nodes_idx_and_Yshunt) ]

    return (; Ynet,
            nodes_idx_with_adjacent_nodes_idx )
      
           
end


function get_edges_Ybr_cal_and_edges_orientation(
    netd )

    # -----------------------------------------------

    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt(
            netd.nodes )
    
    nodes_incident_edges =
        get_nodes_incident_edges(
            netd.edges )

    
    edges_Ybr_cal =
        [ Symbol(typeof(a_branch )) == :PiModelLine ?
        calc_branch_Ybr(
            a_branch.param_values[1:3]..., 1, 1 ) :
                calc_branch_Ybr(
                    a_branch.param_values[1:4]..., 1 )
          for a_branch in
              collect(values(netd.edges)) ]

    edges_orientation =
        [ a_branch.orientation
          for a_branch in
              collect(values(netd.edges)) ]

    edges_orientation_and_edges_Ybr_cal =
        [ (orient, a_Ybr)
          for (orient, a_Ybr) in
              zip( edges_orientation,
                   edges_Ybr_cal )] 

    #-----------------------------------------
 
    nodes_incident_edges_and_orientation =
        [ edges_orientation[
            a_node_incident_edges ]
          for a_node_incident_edges in
              nodes_incident_edges ]


    # nodes_idx_with_adjacent_nodes_idx
    # nodes_node_idx_and_incident_edges_other_node_idx
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx]]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            enumerate(
                nodes_incident_edges_and_orientation)]
          
    #------------------------------------------

    nodes_incident_edges_orientation_and_Ybr_cal =[
        edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges]
    
    return (; 
            edges_Ybr_cal,
            edges_orientation )
      
           
end

#-----------------------------------------------------
#-----------------------------------------------------

function get_dynamics_node_powerflow_net_parameters(
    netd )

    # --------------------------------------------------
    # NetworkData
    # ---------------------------------------------------

    nodes_incident_edges = get_nodes_incident_edges(
        netd.edges )

    edges_Ybr_cal_and_orientation =
        [ Symbol(typeof(a_branch )) == :PiModelLine ?
        (calc_branch_Ybr(a_branch.param_values[1:3]..., 1,1),
         a_branch.orientation ) :
             (calc_branch_Ybr(a_branch.param_values[1:4]...,1),
              a_branch.orientation )
          for a_branch in collect(values( netd.edges )) ]

    edges_Ybr_cal_and_orientation_view =
        [view(edges_Ybr_cal_and_orientation, idx:idx)
         for idx in 1:length(edges_Ybr_cal_and_orientation)]

    nodes_incident_edges_Ybr_cal_and_orientation_view =
        [edges_Ybr_cal_and_orientation_view[
            nodes_incident_edges[idx]]
         for idx in  1:length(nodes_incident_edges) ]

    nodes_idx_and_incident_edges =
        [(idx, incident_edges)
         for (idx, incident_edges) in
             enumerate(nodes_incident_edges)]
    
    return (nodes_incident_edges,
            edges_Ybr_cal_and_orientation,
            nodes_incident_edges_Ybr_cal_and_orientation_view,
            nodes_idx_and_incident_edges)
    
end


function get_dynamics_node_powerflow_net_parameters(
    ; dynamics_case = dynamics_case )

    # --------------------------------------------------
    # NetworkData
    # --------------------------------------------------

    netd  = NetworkData(
        dynamics_case()... )

    return get_dynamics_node_powerflow_net_parameters(
        netd  )

    
end

#-------------------------------------------------------
#-------------------------------------------------------

function get_powerflow_net_variables_view(
    state, pf_net_param, netd; init_pf = true )

    (pf_net, pf_idx_and_state,
     pf_param_views, pf_limits,
     pf_Idx, pf_net_misc) =
         pf_net_param
    
    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net
    
    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state 

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

   load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
     nodes_u_Idx,
     gens_idx,
     ur_IDX,
     ui_IDX,
     vh_IDX,
     θh_IDX,
     red_vh_θh_idx,
     ur_idx,
     ui_idx,
     ur_ui_idx) =
         pf_Idx
    
    # ----------------------------------------------------
    # ----------------------------------------------------

    state_view = view(state, 1:length(state))

    # ----------------------------------------------------
    
    # for pf flat start, ur = vh = 1,  ui = θh = 0
    
    if init_pf == true

        state_view[ur_idx] .= ones(  length( vh_IDX ))
        state_view[ui_idx] .= zeros( length( θh_IDX ))

    end
    
    # ----------------------------------------------------

    # x  = similar( state_x0 )

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]
    

    nodes_idx_and_δ_ed_dash_eq_dash_view =
        get_nodes_idx_and_δ_ed_dash_eq_dash_view(
            state, netd )

    nodes_idx_and_δ_ω_ed_dash_eq_dash_view =
        get_nodes_idx_and_δ_ω_ed_dash_eq_dash_view(
            state, netd )
    
    
    nodes_u_view  =
        [ view(state, nodes_u_Idx[Ind])
          for Ind in collect(1:length(nodes_u_Idx)) ]

    # ----------------------------------------------------    
    
    nodes_pf_U_view  =
        [ view(pf_state , nodes_u_Idx[Ind])
          for Ind in collect(1:length(nodes_u_Idx)) ] 
    
    # ----------------------------------------------------
    
    uh_state_view =
        state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]
    
    x0_ur_ui = [
        state_view[ur_ui_idx][ ur_IDX ]...;
        state_view[ur_ui_idx][ ui_IDX ]...]
        
    x0_vh_θh = [
        abs.(uh_state_view)...;
        angle.(uh_state_view)...]

    working_vh_θh_view =
        view(x0_vh_θh, 1:length(x0_vh_θh))

    # ----------------------------------------------------
    
    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet = zeros(ComplexF64, length(uh_state_view))

    # [ Inet ] # [ view(Inet, 1:length(Inet)) ]
    
    Inet_view =  view( Inet, 1:length(Inet) )


    Iinj = zeros(ComplexF64, length(uh_state_view))

    Iinj_view =  view(Iinj, 1:length(Iinj))
    
    
    # ----------------------------------------------------
    # ----------------------------------------------------

    pf_net_vars_views = (
        ; nodes_pf_U_view,
        Inet_view,
        Iinj_view,
        state_view,
        nodes_u_view,
        nodes_idx_and_δ_ω_ed_dash_eq_dash_view,
        x0_vh_θh,
        working_vh_θh_view )
    
    return pf_net_vars_views 
    
end

# ------------------------------------------------------
# ------------------------------------------------------


function get_powerflow_net_parameters( netd )

    # edges_Ybr_cal, edges_orientation, Ynet, nodes_node_idx_and_incident_edges_other_node_idx =  get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )

    (edges_Ybr_cal,
     Ynet,
     edges_orientation,
     nodes_node_idx_and_incident_edges_other_node_idx) =
         get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------

    Ybus =
        get_Ybus_from_Ynet(
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx)
      
    #------------------------------------------        
    
    slack_bus_idx =
        get_components_slack_Idx(
            netd.nodes )

    slack_vh  =
        get_components_slack_vh(
            netd.nodes )

    gens_Idx_and_vh =
        get_generators_Idx_and_vh(
            netd.nodes )

    non_slack_gens_Idx_and_vh =
        get_non_slack_generators_Idx_and_vh(
            netd.nodes )

    slack_ur_ui_Idx_in_state =
        get_components_slack_ur_ui_Idx_in_state(
            netd.nodes )

    non_slack_ur_ui_Idx_in_state =
        get_components_no_slack_ur_ui_Idx_in_state(
            netd.nodes )
    
    ur_ui_Idx_in_state =
        get_components_ur_ui_Idx_in_state(
            netd.nodes )
    
    nodes_u_Idx =
        get_components_ur_ui_Idx_in_state(
            netd.nodes )

    #------------------------------------------        
    
    gens_idx = first.(gens_Idx_and_vh)    
    gens_vh  = last.(gens_Idx_and_vh)

    #------------------------------------------        
    
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   = [ length(ur_idx), length(ui_idx) ]
    ur_ui_offset = create_offsets( ur_ui_dims )
    ur_ui_IDX    = create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    red_vh_θh_idx = [ setdiff(vh_IDX, gens_idx)...;
               setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ----------------------------------------------------
    
    load_trans_nodes_Idx_and_vlimits =
        get_load_trans_nodes_Idx_and_vlimits( netd.nodes )

    # ----------------------------------------------------

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash, :X_q_dash ] )

    ra_Xd_dash_Xq_dash_view  =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d_dash, :X_q_dash  ] )
    
    # ----------------------------------------------------
    
    ra_Xd_Xq_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d, :X_q ] )
        
    P_Q_nodes_view =
        get_components_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )
    
    P_Q_gens_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )

    P_Q_gens_loc_load_view =
        get_gens_loc_load_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :loc_P, :loc_Q ] )

    P_Q_non_gens_view =
        get_non_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )


    #------------------------------------------        

    pf_and_dyn_idx_and_Idx =
        get_pf_and_dyn_idx_and_Idx(netd)
    
    #------------------------------------------
    
    P_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :Q ] )...;]
    
    P_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes ; some_param = [ :Q ] )...;]
    
    # δ_ed_eq_pf_NL_para =
    #     [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;]
    
    P_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_P ] )...;]
    
    Q_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_Q  ] )...;]

    pf_NL_P_Q_para = (; P_gens_NL_para,
                   Q_gens_NL_para,
                   P_load_NL_para,
                   Q_load_NL_para,
                   P_g_loc_load_NL_para,
                   Q_g_loc_load_NL_para)
  
    # pf_NL_P_Q_para 
    #------------------------------------------
    
    pf_net_misc = (; pf_NL_P_Q_para,  )
        
    #------------------------------------------
    
    pf_net =
        (; Ybus,
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx,
         edges_Ybr_cal,
         edges_orientation )
    
    pf_idx_and_state =
        (; slack_vh,
         gens_vh,
         gens_Idx_and_vh,
         non_slack_gens_Idx_and_vh,
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state )    

    pf_param_views =
        (; ra_Xd_Xq_view,
         ra_Xd_dash_Xq_dash_view,
         ra_Xd_Xq_Xd_dash_Xq_dash_view,
         P_Q_nodes_view,
         P_Q_gens_view,
         P_Q_gens_loc_load_view,
         P_Q_non_gens_view )

    pf_limits =
        (; load_trans_nodes_Idx_and_vlimits, )

    pf_Idx =
        (; slack_bus_idx,
         nodes_u_Idx,
         gens_idx,
         ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX,
         red_vh_θh_idx,
         ur_idx,
         ui_idx,
         ur_ui_idx )
    
    pf_net_param =
        (; pf_net,
         pf_idx_and_state,
         pf_param_views,
         pf_limits,
         pf_Idx,
         pf_and_dyn_idx_and_Idx,
         pf_net_misc)


    return pf_net_param
end



function get_powerflow_net_parameters(
    ; dynamics_case = file )

    # dynamics_case = new_case_IEEE_9_Bus_dynamic_plant_SM_v6_P

    netd  = NetworkData( dynamics_case()... )

    # ----------------------------------------------------
    # ----------------------------------------------------

    return get_powerflow_net_parameters( netd )
end

#-----------------------------------------------------------


function get_streamlined_powerflow_net_parameters( netd )

    # edges_Ybr_cal, edges_orientation, Ynet, nodes_node_idx_and_incident_edges_other_node_idx =  get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )

    (edges_Ybr_cal,
     Ynet,
     edges_orientation,
     nodes_node_idx_and_incident_edges_other_node_idx) =
         get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------

    Ybus =
        get_Ybus_from_Ynet(
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx)
      
    #------------------------------------------        
    
    slack_bus_idx =
        get_components_slack_Idx(
            netd.nodes )

    slack_vh  =
        get_components_slack_vh(
            netd.nodes )

    gens_Idx_and_vh =
        get_generators_Idx_and_vh(
            netd.nodes )

    non_slack_gens_Idx_and_vh =
        get_non_slack_generators_Idx_and_vh(
            netd.nodes )

    slack_ur_ui_Idx_in_state =
        get_components_slack_ur_ui_Idx_in_state(
            netd.nodes )

    non_slack_ur_ui_Idx_in_state =
        get_components_no_slack_ur_ui_Idx_in_state(
            netd.nodes )
    
    ur_ui_Idx_in_state =
        get_components_ur_ui_Idx_in_state(
            netd.nodes )
    
    nodes_u_Idx =
        get_components_ur_ui_Idx_in_state(
            netd.nodes )

    #------------------------------------------        
    
    gens_idx = first.(gens_Idx_and_vh)
    
    gens_vh  = last.(gens_Idx_and_vh)

    #------------------------------------------        
    
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   = [ length(ur_idx), length(ui_idx) ]
    ur_ui_offset = create_offsets( ur_ui_dims )
    ur_ui_IDX    = create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    red_vh_θh_idx = [ setdiff(vh_IDX, gens_idx)...;
               setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ----------------------------------------------------
    
    load_trans_nodes_Idx_and_vlimits =
        get_load_trans_nodes_Idx_and_vlimits( netd.nodes )

    # ----------------------------------------------------

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes; 
            param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash, :X_q_dash ],
            gens_view_only = true  )

    ra_Xd_dash_Xq_dash_view  =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d_dash, :X_q_dash  ],
        gens_view_only = true)
    
    # ----------------------------------------------------
    
    ra_Xd_Xq_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d, :X_q ],
        gens_view_only = true)
        
    P_Q_nodes_view =
        get_components_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )
    
    P_Q_gens_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ],
        gens_view_only = true )

    # P_Q_non_gens_view =
    #     get_non_gens_params_view_in_param_values(
    #         netd.nodes_param, netd.nodes;
    #         param_list = [ :P, :Q ] )


    P_Q_non_gens_view =
        get_streamlined_non_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ],
        non_gens_view_only = true )
    

    # P_Q_gens_loc_load_view =
    #     get_gens_loc_load_params_view_in_param_values(
    #         netd.nodes_param, netd.nodes;
    #         param_list = [ :loc_P, :loc_Q ] )


    P_Q_gens_loc_load_view =
        get_streamlined_gens_loc_load_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :loc_P, :loc_Q ],
            gens_view_only = true )
    
    #------------------------------------------        
    #------------------------------------------        

    pf_and_dyn_idx_and_Idx =
        get_pf_and_dyn_idx_and_Idx(netd)
    
    #------------------------------------------      
    #------------------------------------------

    P_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :Q ] )...;]
    
    P_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes ; some_param = [ :Q ] )...;]
    
    # δ_ed_eq_pf_NL_para =
    #     [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;]
    
    P_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_P ] )...;]
    
    Q_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_Q  ] )...;]

    pf_NL_P_Q_para = (; P_gens_NL_para,
                   Q_gens_NL_para,
                   P_load_NL_para,
                   Q_load_NL_para,
                   P_g_loc_load_NL_para,
                   Q_g_loc_load_NL_para)
  
    # pf_NL_P_Q_para 
    #------------------------------------------
    
    pf_net_misc = (; pf_NL_P_Q_para,  )
    
    #------------------------------------------        
    
    pf_net =
        (; Ybus,
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx,
         edges_Ybr_cal,
         edges_orientation )
    
    pf_idx_and_state =
        (; slack_vh,
         gens_vh,
         gens_Idx_and_vh,
         non_slack_gens_Idx_and_vh,
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state )    

    pf_param_views =
        (; ra_Xd_Xq_view,
         ra_Xd_dash_Xq_dash_view,
         ra_Xd_Xq_Xd_dash_Xq_dash_view,
         P_Q_nodes_view,
         P_Q_gens_view,
         P_Q_gens_loc_load_view,
         P_Q_non_gens_view )

    pf_limits =
        (; load_trans_nodes_Idx_and_vlimits, )

    pf_Idx =
        (; slack_bus_idx,
         nodes_u_Idx,
         gens_idx,
         ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX,
         red_vh_θh_idx,
         ur_idx,
         ui_idx,
         ur_ui_idx )
    
    pf_net_param =
        (; pf_net,
         pf_idx_and_state,
         pf_param_views,
         pf_limits,
         pf_Idx,
         pf_and_dyn_idx_and_Idx,
         pf_net_misc )


    return pf_net_param
end


function get_streamlined_powerflow_net_parameters(
    ; dynamics_case = file )

    # dynamics_case = new_case_IEEE_9_Bus_dynamic_plant_SM_v6_P

    netd  = NetworkData( dynamics_case()... )

    # ----------------------------------------------------
    # ----------------------------------------------------

    return get_streamlined_powerflow_net_parameters( netd )
end

# ----------------------------------------------------

"""
Idx_converstion_fun:

`net_to_im_model_indices(idx, dict_conv)`

`net_to_industrial_model_indices(idx, dict_conv)`

dict_conv:

dict_sys_to_im =
    get_net_to_im_indices_dict( netd  )

dict_sys_to_industry =
        get_net_to_industrial_model_indices_dict( netd; no_control_device = false   )

get_a_model_streamlined_powerflow_net_parameters(
    netd;
    dict_sys_to_model_Idx = dict_sys_to_model_Idx,
    Idx_converstion_fun = Idx_converstion_fun )


"""
function get_a_model_streamlined_powerflow_net_parameters(
    netd;
    dict_sys_to_model_Idx = dict_sys_to_model_Idx,
    Idx_converstion_fun = Idx_converstion_fun )

    
    (;edges_Ybr_cal,
     Ynet,
     edges_orientation,
     nodes_node_idx_and_incident_edges_other_node_idx) =
         get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------

    Ybus =
        get_Ybus_from_Ynet(
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx)
      
    #------------------------------------------        
    
    slack_bus_idx =
        get_components_slack_Idx(
            netd.nodes )

    slack_vh  =
        get_components_slack_vh(
            netd.nodes )

    gens_Idx_and_vh =
        get_generators_Idx_and_vh(
            netd.nodes )

    non_slack_gens_Idx_and_vh =
        get_non_slack_generators_Idx_and_vh(
            netd.nodes )

    # Idx_converstion_fun dict_sys_to_model_Idx,


    #------------------------------------------        
    
    gens_idx = first.(gens_Idx_and_vh)    
    gens_vh  = last.(gens_Idx_and_vh)

    #------------------------------------------            
    #------------------------------------------            
    
    slack_ur_ui_Idx_in_state =
        Idx_converstion_fun(
            get_components_slack_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_model_Idx)

    # 

    non_slack_ur_ui_Idx_in_state =
        Idx_converstion_fun(
            get_components_no_slack_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_model_Idx)

    # 
    
    ur_ui_Idx_in_state =
        Idx_converstion_fun(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_model_Idx)
    #  
    
    nodes_u_Idx =
        Idx_converstion_fun(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_model_Idx)

    #------------------------------------------        
    #------------------------------------------            

    
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   = [length(ur_idx),
                    length(ui_idx) ]
    
    ur_ui_offset = create_offsets(ur_ui_dims )
    
    ur_ui_IDX    = create_idxs(
        ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    red_vh_θh_idx = [setdiff(vh_IDX, gens_idx)...;
               setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ----------------------------------------------------
    
    load_trans_nodes_Idx_and_vlimits =
        get_load_trans_nodes_Idx_and_vlimits(
            netd.nodes )

    # ----------------------------------------------------

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes; 
            param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash, :X_q_dash ],
            gens_view_only = true )

    ra_Xd_dash_Xq_dash_view  =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d_dash, :X_q_dash ],
        gens_view_only = true)
    
    # ----------------------------------------------------
    
    ra_Xd_Xq_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d, :X_q ],
        gens_view_only = true)
        
    P_Q_nodes_view =
        get_components_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )
    
    P_Q_gens_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ],
        gens_view_only = true )

    # P_Q_non_gens_view =
    #     get_non_gens_params_view_in_param_values(
    #         netd.nodes_param, netd.nodes;
    #         param_list = [ :P, :Q ] )


    P_Q_non_gens_view =
        get_streamlined_non_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ],
        non_gens_view_only = true )
    

    # P_Q_gens_loc_load_view =
    #     get_gens_loc_load_params_view_in_param_values(
    #         netd.nodes_param, netd.nodes;
    #         param_list = [ :loc_P, :loc_Q ] )


    P_Q_gens_loc_load_view =
        get_streamlined_gens_loc_load_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :loc_P, :loc_Q ],
            gens_view_only = true )
    
    #------------------------------------------        
    #------------------------------------------        

    P_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :Q ] )...;]
    
    P_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes ; some_param = [ :Q ] )...;]
    
    # δ_ed_eq_pf_NL_para =
    #     [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;]
    
    P_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ;
            some_param = [ :loc_P ] )...;]
    
    Q_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_Q  ] )...;]

    pf_NL_P_Q_para = (; P_gens_NL_para,
                   Q_gens_NL_para,
                   P_load_NL_para,
                   Q_load_NL_para,
                   P_g_loc_load_NL_para,
                   Q_g_loc_load_NL_para)
  
    # pf_NL_P_Q_para 
    #------------------------------------------        
    #------------------------------------------        
    
    # pf_and_dyn_idx_and_Idx =
    #     get_pf_and_dyn_idx_and_Idx(netd)
    
    #------------------------------------------
    #------------------------------------------        

    pf_and_dyn_idx_and_Idx =
        get_pf_and_dyn_idx_and_Idx(netd)
    
    #------------------------------------------
    #------------------------------------------

    P_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :Q ] )...;]
    
    P_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes ;
            some_param = [ :Q ] )...;]
    
    # δ_ed_eq_pf_NL_para =
    #     [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;]
    
    P_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ;
            some_param = [ :loc_P ] )...;]
    
    Q_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_Q  ] )...;]

    pf_NL_P_Q_para = (; P_gens_NL_para,
                   Q_gens_NL_para,
                   P_load_NL_para,
                   Q_load_NL_para,
                   P_g_loc_load_NL_para,
                   Q_g_loc_load_NL_para)
  
    # pf_NL_P_Q_para 
    #------------------------------------------
    
    pf_net_misc = (; pf_NL_P_Q_para,  ) 
    
    #------------------------------------------
    
    pf_net =
        (; Ybus,
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx,
         edges_Ybr_cal,
         edges_orientation )
    
    pf_idx_and_state =
        (; slack_vh,
         gens_vh,
         gens_Idx_and_vh,
         non_slack_gens_Idx_and_vh,
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state )    

    pf_param_views =
        (; ra_Xd_Xq_view,
         ra_Xd_dash_Xq_dash_view,
         ra_Xd_Xq_Xd_dash_Xq_dash_view,
         P_Q_nodes_view,
         P_Q_gens_view,
         P_Q_gens_loc_load_view,
         P_Q_non_gens_view )

    pf_limits =
        (; load_trans_nodes_Idx_and_vlimits, )

    pf_Idx =
        (; slack_bus_idx,
         nodes_u_Idx,
         gens_idx,
         ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX,
         red_vh_θh_idx,
         ur_idx,
         ui_idx,
         ur_ui_idx )
    
    pf_net_param =
        (; pf_net,
         pf_idx_and_state,
         pf_param_views,
         pf_limits,
         pf_Idx,
         pf_and_dyn_idx_and_Idx,
         pf_net_misc)


    return pf_net_param
end


function get_a_model_streamlined_powerflow_net_parameters(
    ; dynamics_case = file,
    dict_sys_to_model_Idx = dict_sys_to_model_Idx,
    Idx_converstion_fun = Idx_converstion_fun )

    # dynamics_case = new_case_IEEE_9_Bus_dynamic_plant_SM_v6_P

    netd  = NetworkData( dynamics_case()... )

    # ----------------------------------------------------
    # ----------------------------------------------------

    return get_a_model_streamlined_powerflow_net_parameters(
        netd;
        dict_sys_to_model_Idx = dict_sys_to_model_Idx,
        Idx_converstion_fun = Idx_converstion_fun )
end

# ------------------------------------------------------
# im pf
# ------------------------------------------------------

function get_im_model_powerflow_net_parameters( netd )

 

    dict_sys_to_industry =
        get_net_to_im_indices_dict( netd  )

    (edges_Ybr_cal,
     Ynet,
     edges_orientation,
     nodes_node_idx_and_incident_edges_other_node_idx) =
         get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------

    Ybus =
        get_Ybus_from_Ynet(
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx)
      
    #------------------------------------------        

    slack_bus_idx =
        get_components_slack_Idx(
            netd.nodes )

    slack_vh  =
        get_components_slack_vh(
            netd.nodes )

    gens_Idx_and_vh =
        get_generators_Idx_and_vh(
            netd.nodes )

    non_slack_gens_Idx_and_vh =
        get_non_slack_generators_Idx_and_vh(
            netd.nodes )

    #------------------------------------------
    
    gens_idx = first.(gens_Idx_and_vh)    
    gens_vh  = last.(gens_Idx_and_vh)
    
    #------------------------------------------        

    slack_ur_ui_Idx_in_state =
        net_to_im_model_indices(
            get_components_slack_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )

    non_slack_ur_ui_Idx_in_state =
        net_to_im_model_indices(
            get_components_no_slack_ur_ui_Idx_in_state(
                netd.nodes),
            dict_sys_to_industry  )

    ur_ui_Idx_in_state =
        net_to_im_model_indices(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )

    nodes_u_Idx =
        net_to_im_model_indices(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )
    
    #------------------------------------------       
    #------------------------------------------       
        
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   =
        [ length(ur_idx), length(ui_idx) ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    red_vh_θh_idx =
        [ setdiff(vh_IDX, gens_idx)...;
          setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # -----------------------------------------------
    
    load_trans_nodes_Idx_and_vlimits =
        get_load_trans_nodes_Idx_and_vlimits(
            netd.nodes )

    # -----------------------------------------------

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash,
                :X_q_dash ] )

    ra_Xd_dash_Xq_dash_view  =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d_dash,
                           :X_q_dash  ] )
    
    # ---------------------------------------------
    
    ra_Xd_Xq_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d, :X_q ] )
        
    P_Q_nodes_view =
        get_components_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )
    
    P_Q_gens_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )

    P_Q_gens_loc_load_view =
        get_gens_loc_load_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :loc_P, :loc_Q ] )

    P_Q_non_gens_view =
        get_non_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )

    #------------------------------------------
    
    pf_and_dyn_idx_and_Idx =
        get_pf_and_dyn_idx_and_Idx(netd)
    
    #------------------------------------------

    P_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :Q ] )...;]
    
    P_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes ;
            some_param = [ :Q ] )...;]
    
    # δ_ed_eq_pf_NL_para =
    #     [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;]
    
    P_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ;
            some_param = [ :loc_P ] )...;]
    
    Q_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_Q  ] )...;]

    pf_NL_P_Q_para = (; P_gens_NL_para,
                   Q_gens_NL_para,
                   P_load_NL_para,
                   Q_load_NL_para,
                   P_g_loc_load_NL_para,
                   Q_g_loc_load_NL_para)
  
    # pf_NL_P_Q_para 
    #------------------------------------------
    
    pf_net_misc = (; pf_NL_P_Q_para,  ) 
        
    #------------------------------------------
    
    pf_net =
        (; Ybus,
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx,
         edges_Ybr_cal,
         edges_orientation )
    
    pf_idx_and_state =
        (; slack_vh,
         gens_vh,
         gens_Idx_and_vh,
         non_slack_gens_Idx_and_vh,
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state )    

    pf_param_views =
        (; ra_Xd_Xq_view,
         ra_Xd_dash_Xq_dash_view,
         ra_Xd_Xq_Xd_dash_Xq_dash_view,
         P_Q_nodes_view,
         P_Q_gens_view,
         P_Q_gens_loc_load_view,
         P_Q_non_gens_view )

    pf_limits =
        (; load_trans_nodes_Idx_and_vlimits, )

    pf_Idx =
        (; slack_bus_idx,
         nodes_u_Idx,
         gens_idx,
         ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX,
         red_vh_θh_idx,
         ur_idx,
         ui_idx,
         ur_ui_idx )
    
    pf_net_param =
        (; pf_net,
         pf_idx_and_state,
         pf_param_views,
         pf_limits,
         pf_Idx,
         pf_and_dyn_idx_and_Idx,
         pf_net_misc)


    return pf_net_param
end

#--------------------------------------------------


function get_im_model_powerflow_net_parameters(
    ; dynamics_case = file )

    # dynamics_case =
    #     new_case_IEEE_9_Bus_dynamic_plant_SM_v6_P

    netd  = NetworkData( dynamics_case()... )

    # ----------------------------------------------
    # ----------------------------------------------

    
    return get_im_model_powerflow_net_parameters(
        netd )
end



# --------------------------------------------------
# Industrial pf
# --------------------------------------------------

function get_industrial_model_powerflow_net_parameters(
    netd; no_control_device = false )

 

    if no_control_device == false

        dict_sys_to_industry =
            get_net_to_industrial_model_indices_dict( netd  )
        
    else
        dict_sys_to_industry =
            get_net_to_industrial_model_indices_dict(
                netd;  no_control_device = true  )
    end
    

    (edges_Ybr_cal,
     Ynet,
     edges_orientation,
     nodes_node_idx_and_incident_edges_other_node_idx) =
         get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------

    Ybus =
        get_Ybus_from_Ynet(
            Ynet,
            nodes_node_idx_and_incident_edges_other_node_idx)
      
    #------------------------------------------        

    slack_bus_idx =
        get_components_slack_Idx( netd.nodes )

    slack_vh  =
        get_components_slack_vh( netd.nodes )

    gens_Idx_and_vh =
        get_generators_Idx_and_vh( netd.nodes )

    non_slack_gens_Idx_and_vh =
        get_non_slack_generators_Idx_and_vh(
            netd.nodes )

    slack_ur_ui_Idx_in_state =
        net_to_industrial_model_indices(
            get_components_slack_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )

    non_slack_ur_ui_Idx_in_state =
        net_to_industrial_model_indices(
            get_components_no_slack_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )

    ur_ui_Idx_in_state =
        net_to_industrial_model_indices(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )

    nodes_u_Idx =
        net_to_industrial_model_indices(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )

    #------------------------------------------        
    
    gens_idx = first.(gens_Idx_and_vh)    
    gens_vh  = last.(gens_Idx_and_vh)

    #------------------------------------------        
        
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   = [ length(ur_idx), length(ui_idx) ]
    
    ur_ui_offset = create_offsets( ur_ui_dims )
    
    ur_ui_IDX    = create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    red_vh_θh_idx =
        [ setdiff(vh_IDX, gens_idx)...;
          setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ----------------------------------------------------
    
    load_trans_nodes_Idx_and_vlimits =
        get_load_trans_nodes_Idx_and_vlimits( netd.nodes )

    # ----------------------------------------------------

    ra_Xd_Xq_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d, :X_q,
                           :X_d_dash,
                           :X_q_dash ] )

    ra_Xd_dash_Xq_dash_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d_dash,
                           :X_q_dash  ] )
    
    # ----------------------------------------------------
    
    ra_Xd_Xq_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d, :X_q ] )
        
    P_Q_nodes_view =
        get_components_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )
    
    P_Q_gens_view =
        get_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )

    P_Q_gens_loc_load_view =
        get_gens_loc_load_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :loc_P, :loc_Q ] )

    P_Q_non_gens_view =
        get_non_gens_params_view_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )

    #------------------------------------------

    pf_and_dyn_idx_and_Idx =
        get_pf_and_dyn_idx_and_Idx(netd)

    #------------------------------------------

    P_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :Q ] )...;]
    
    P_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_load_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes ;
            some_param = [ :Q ] )...;]
    
    # δ_ed_eq_pf_NL_para =
    #     [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;]
    
    P_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_P ] )...;]
    
    Q_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_Q  ] )...;]

    pf_NL_P_Q_para = (; P_gens_NL_para,
                   Q_gens_NL_para,
                   P_load_NL_para,
                   Q_load_NL_para,
                   P_g_loc_load_NL_para,
                   Q_g_loc_load_NL_para)
  
    # pf_NL_P_Q_para 
    #------------------------------------------
    
    pf_net_misc = (; pf_NL_P_Q_para,  ) 
    
    #------------------------------------------
    
    
    pf_net =
        (; Ybus,
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx,
         edges_Ybr_cal,
         edges_orientation )
    
    pf_idx_and_state =
        (; slack_vh,
         gens_vh,
         gens_Idx_and_vh,
         non_slack_gens_Idx_and_vh,
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state )    

    pf_param_views =
        (; ra_Xd_Xq_view,
         ra_Xd_dash_Xq_dash_view,
         ra_Xd_Xq_Xd_dash_Xq_dash_view,
         P_Q_nodes_view, P_Q_gens_view,
         P_Q_gens_loc_load_view,
         P_Q_non_gens_view )

    pf_limits =
        (; load_trans_nodes_Idx_and_vlimits, )

    pf_Idx =
        (; slack_bus_idx,
         nodes_u_Idx,
         gens_idx,
         ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX,
         red_vh_θh_idx,
         ur_idx,
         ui_idx,
         ur_ui_idx )
    
    pf_net_param =
        (; pf_net,
         pf_idx_and_state,
         pf_param_views,
         pf_limits,
         pf_Idx,
         pf_and_dyn_idx_and_Idx,
         pf_net_misc)


    return pf_net_param
end

#-----------------------------------------------------------


function get_industrial_model_powerflow_net_parameters(
    ; dynamics_case = file )

    # dynamics_case =
    #     new_case_IEEE_9_Bus_dynamic_plant_SM_v6_P

    netd  = NetworkData( dynamics_case()... )

    # ------------------------------------------------
    # -----------------------------------------------

    return get_industrial_model_powerflow_net_parameters(
        netd )
end


# -------------------------------------------------------
# get_pf_param_and_sys_init
# ------------------------------------------------------

function get_pf_sys_param_sys_views_sys_init(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false )

    # -----------------------------------------------------

    pf_net_param =
        get_powerflow_net_parameters(
            netd )
    
    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits, pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) =
         pf_net_param

    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) =
         pf_net

    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
     nodes_u_Idx,
     gens_idx,
     ur_IDX,
     ui_IDX,
     vh_IDX,
     θh_IDX,
     red_vh_θh_idx,
     ur_idx,
     ui_idx,
     ur_ui_idx) =
         pf_Idx
    
    #----------------------------------------------------

    syms = get_network_vars_labels( netd )
    
    #----------------------------------------------------

    state = zeros(length(
        get_network_vars_labels( netd ) ))

    #----------------------------------------------------

    δ_ed_dash_eq_dash_view =
        get_nodes_δ_ed_dash_eq_dash_view(
            state, netd )

    δ_ω_ed_dash_eq_dash_view =
        get_nodes_δ_ω_ed_dash_eq_dash_view(
            state, netd )
    
    #----------------------------------------------------    
    # ----------------------------------------------------

    state_view = view(state, 1:length(state))

    # ----------------------------------------------------
    
    # for pf flat start, ur = vh = 1,  ui = θh = 0

    # remove init_pf = true
    
    
    
    if init_pf == true 

        state_view[ur_idx] .= ones(  length( ur_IDX ))
        state_view[ui_idx] .= zeros( length( ui_IDX ))

    end
    
    # ----------------------------------------------------

    # x  = similar( state_x0 )

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]
        
    #----------------------------------------------------
    
    nodes_u_view  = [
        view(state, nodes_u_Idx[Ind])
        for Ind in
            collect(1:length(nodes_u_Idx)) ]

    # ----------------------------------------------------  
    
    nodes_pf_U_view  = [
        view(pf_state , nodes_u_Idx[Ind])
        for Ind in
            collect(1:length(nodes_u_Idx)) ] 
    
    # ----------------------------------------------------
    
    uh_state = state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]
    
    x0_ur_ui = [state_view[ur_ui_idx][ ur_IDX ]...;
                state_view[ur_ui_idx][ ui_IDX ]...]
        
    x0_vh_θh = [abs.(uh_state)...; angle.(uh_state)...]

    working_vh_θh_view = view(x0_vh_θh, 1:length(x0_vh_θh))

    x0_vh_view = @view x0_vh_θh[vh_IDX]

    x0_θh_view = @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0 = [ red_vh_θh_0_view; ]

    mismatch = similar( red_vh_θh_0 )
    
    # ----------------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))
    
    # ----------------------------------------------------
    
    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet = zeros(ComplexF64, length( uh_state ))
    
    Inet_view = view( Inet, 1:length( Inet ) )

    Iinj = zeros(ComplexF64, length( uh_state ))

    Iinj_view = view(Iinj, 1:length( Iinj ))
        
    #----------------------------------------------------

    global_pf_views = (
        working_vh_θh_view,
        nodes_pf_U_view,
        Inet_view,
        Iinj_view )

    sd_pf_views = (
        working_vh_θh_view,
        nodes_pf_U_view,
        Inet_view,
        Iinj_view,
        δ_ω_ed_dash_eq_dash_view )

    # ---------------------------------------------------
    # global_pf param
    # ---------------------------------------------------
    
    global_pf_param = (
        pf_net_param,
        sd_pf_views,
        mismatch )

    
    branches_name =
        collect(keys( netd.edges ))
    
    nodes_name =
        collect(keys( netd.nodes ))    

    #----------------------------------------------------
    #----------------------------------------------------
    
    named_tup_pf_result =
        power_balance_powerflow(
            x0_vh_θh, mismatch, sd_pf_views,
            (nodes_name, branches_name) ,
            pf_net_param ;
            maxiter=maxiter,
            ftol=ftol ,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )

    #----------------------------------------------------

    bus_dict_init    = named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init
    
    #---------------------------------------------------
    
    external_pf_init_post_pf_parameters!(
        netd, bus_dict_init, branch_dict_init )

    #---------------------------------------------------

    state .= external_init_operationpoint(
        netd, bus_dict_init, branch_dict_init )

    #----------------------------------------------------
        
    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)

    #---------------------------------------------------

    nodes_f_t  = external_get_nodes_or_edges_f_t(
        netd.nodes, bus_dict_init )
    
    #---------------------------------------------------

    return (; nodes_f_t, nodes_cb_sw, state,
            global_pf_param, named_tup_pf_result )
    
end

# -------------------------------------------------------
# get_pf_in_named_tuple_with_state
# ------------------------------------------------------

function get_pf_in_named_tuple_with_state(
    netd;
    maxiter=40,
    ftol=1000*eps() ,
    xtol=1000*eps(),
    init_pf = true,
    with_δ_ed_eq = false )

    # -----------------------------------------------------

    pf_net_param = get_powerflow_net_parameters( netd )
    
    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx ,
     pf_net_misc) = pf_net_param

    (Ybus,
     Ynet,
     nodes_node_idx_and_incident_edges_other_node_idx,
     edges_Ybr_cal,
     edges_orientation) = pf_net

    (slack_vh,
     gens_vh,
     gens_Idx_and_vh,
     non_slack_gens_Idx_and_vh,
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state) =
         pf_idx_and_state

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    (slack_bus_idx,
     nodes_u_Idx,
     gens_idx,
     ur_IDX,
     ui_IDX,
     vh_IDX,
     θh_IDX,
     red_vh_θh_idx,
     ur_idx,
     ui_idx,
     ur_ui_idx) = pf_Idx
    
    #----------------------------------------------------

    syms = get_network_vars_labels( netd )
    
    #----------------------------------------------------

    state = zeros(length(
        get_network_vars_labels( netd ) ))

    #----------------------------------------------------

    δ_ed_dash_eq_dash_view =
        get_nodes_δ_ed_dash_eq_dash_view(
            state, netd )

    δ_ω_ed_dash_eq_dash_view =
        get_nodes_δ_ω_ed_dash_eq_dash_view(
            state, netd )
    
    #----------------------------------------------------    
    # ----------------------------------------------------

    state_view = view(state, 1:length(state))

    # ----------------------------------------------------
    
    # for pf flat start, ur = vh = 1,  ui = θh = 0

    # remove init_pf = true
    
    
    
    if init_pf == true 

        state_view[ur_idx] .= ones(  length( ur_IDX ))
        state_view[ui_idx] .= zeros( length( ui_IDX ))

    end
    
    # ----------------------------------------------------

    # x  = similar( state_x0 )

    # Decoupling pf_state from state_x0

    pf_state = [ state; ]
        
    #----------------------------------------------------
    
    nodes_u_view  = [ view(state, nodes_u_Idx[Ind])
                      for Ind in
                          collect(1:length(nodes_u_Idx)) ]

    # ----------------------------------------------------  
    
    nodes_pf_U_view  = [
        view(pf_state , nodes_u_Idx[Ind])
        for Ind in collect(1:length(nodes_u_Idx)) ] 
    
    # ----------------------------------------------------
    
    uh_state = state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]
    
    x0_ur_ui = [state_view[ur_ui_idx][ ur_IDX ]...;
                state_view[ur_ui_idx][ ui_IDX ]...]
        
    x0_vh_θh = [abs.(uh_state)...; angle.(uh_state)...]

    working_vh_θh_view = view(x0_vh_θh, 1:length(x0_vh_θh))

    x0_vh_view       = @view x0_vh_θh[vh_IDX]

    x0_θh_view       = @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view = @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0      = [ red_vh_θh_0_view; ]

    mismatch         = similar( red_vh_θh_0 )
    
    # ----------------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))
    
    # ----------------------------------------------------
    
    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet = zeros(ComplexF64, length( uh_state ))
    
    Inet_view =  view( Inet, 1:length( Inet ) )

    Iinj = zeros(ComplexF64, length( uh_state ))

    Iinj_view =  view(Iinj, 1:length( Iinj ))
        
    #----------------------------------------------------

    global_pf_views =
        (working_vh_θh_view,
         nodes_pf_U_view,
         Inet_view,
         Iinj_view )

    sd_pf_views =
        ( working_vh_θh_view,
          nodes_pf_U_view,
          Inet_view,
          Iinj_view,
          δ_ω_ed_dash_eq_dash_view )

    # ---------------------------------------------------
    # global_pf param
    # ---------------------------------------------------
    
    global_pf_param =
        ( pf_net_param,
          sd_pf_views,
          mismatch )

    
    branches_name  = collect(keys( netd.edges ))
    
    nodes_name = collect(keys( netd.nodes ))    

    #----------------------------------------------------
    #----------------------------------------------------
    
    named_tup_pf_result =
        power_balance_powerflow(
            x0_vh_θh, mismatch,
            sd_pf_views,
            (nodes_name, branches_name) ,
            pf_net_param ;
            maxiter=maxiter,
            ftol=ftol ,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )

    #----------------------------------------------------

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init =
        named_tup_pf_result.branch_dict_init
    
    #---------------------------------------------------
    
    external_pf_init_post_pf_parameters!(
        netd, bus_dict_init, branch_dict_init )

    #---------------------------------------------------

    state .= external_init_operationpoint(
        netd, bus_dict_init, branch_dict_init )
   
    #---------------------------------------------------

    return (; named_tup_pf_result, state )
    
end

#-------------------------------------------------------

function get_netd_pf_dict_net_vars_param_x0_cb_sw_f_t(
    netd
    ; init_pf =
        true,
    with_δ_ed_eq =
        false )

    # ------------------------------------------------------
      
    branches_name  = collect(keys( netd.edges ))
    
    nodes_name = collect(keys( netd.nodes ))    
    
    # ------------------------------------------------------
    # ------------------------------------------------------

    pf_net_param =
        get_powerflow_net_parameters(
            netd )

    # ------------------------------------------------------

    (pf_net,
     pf_idx_and_state,
     pf_param_views,
     pf_limits,
     pf_Idx,
     pf_and_dyn_idx_and_Idx,
     pf_net_misc) = pf_net_param

    (ra_Xd_Xq_view,
     ra_Xd_dash_Xq_dash_view,
     ra_Xd_Xq_Xd_dash_Xq_dash_view,
     P_Q_nodes_view,
     P_Q_gens_view,
     P_Q_gens_loc_load_view,
     P_Q_non_gens_view) =
         pf_param_views

    (slack_bus_idx,
     nodes_u_Idx,
     gens_idx,
     ur_IDX,
     ui_IDX,
     vh_IDX,
     θh_IDX,
     red_vh_θh_idx,
     ur_idx,
     ui_idx,
     ur_ui_idx) = pf_Idx
        
    # ------------------------------------------------------

    state =
        zeros( length(
            get_network_vars_labels(
                netd ) ))

    # --------------------------------------------------    

    pf_net_vars_views =
        get_powerflow_net_variables_view(
            state, pf_net_param, netd;
            init_pf = init_pf )

    (nodes_pf_U_view,
     Inet_view,
     Iinj_view,
     state_view,
     nodes_u_view,
     nodes_idx_and_δ_ω_ed_dash_eq_dash_view,
     x0_vh_θh,
     working_vh_θh_view) =
         pf_net_vars_views

    # ------------------------------------------------------

    pf_views = (; working_vh_θh_view,
                nodes_u_view,
                nodes_pf_U_view,
                nodes_idx_and_δ_ω_ed_dash_eq_dash_view,
                ra_Xd_Xq_view,
                ra_Xd_dash_Xq_dash_view,
                ra_Xd_Xq_Xd_dash_Xq_dash_view,
                P_Q_nodes_view,
                P_Q_gens_view,
                P_Q_gens_loc_load_view,
                P_Q_non_gens_view,
                Inet_view,
                Iinj_view )

    # constitution of pf_param
    
    pf_param = (; pf_net,
                pf_idx_and_state,
                pf_views,
                pf_limits,
                pf_Idx  )

    dyn_pf_param = (; nodes_pf_U_view,
                    Inet_view,
                    Iinj_view )

    # ---------------------------------------------------
    # ---------------------------------------------------
    # Results 
    # ---------------------------------------------------
    # ---------------------------------------------------  

    pf_result_dict =
        get_power_balance_powerflow_results(
            x0_vh_θh,
            pf_param,
            nodes_name,
            branches_name;
            with_δ_ed_eq =
                with_δ_ed_eq )

    bus_dict_init =
        pf_result_dict["dict_init"]["bus_dict_init"]

    branch_dict_init =
        pf_result_dict["dict_init"]["branch_dict_init"]
    
    #---------------------------------------------------
    
    external_pf_init_post_pf_parameters!(
        netd, bus_dict_init, branch_dict_init )

    state .= external_init_operationpoint(
        netd, bus_dict_init, branch_dict_init )
        
    nodes_cb_sw = get_nodes_cb_sw(netd.nodes)

    nodes_f_t  = external_get_nodes_or_edges_f_t(
        netd.nodes, bus_dict_init )
    
    #---------------------------------------------------------

    return (; dyn_pf_param,
            pf_param,
            pf_result_dict,
            netd,
            state,
            nodes_cb_sw,
            nodes_f_t )

    # return ( netd, pf_result_dict, pf_net_vars_views, pf_net_param, state,  nodes_cb_sw, nodes_f_t )

    
end

# ------------------------------------------------------

function get_netd_pf_dict_net_vars_param_x0_cb_sw_f_t(
    ; dynamics_case = file,
    init_pf = true,
    with_δ_ed_eq = false )
    
    return get_netd_pf_dict_net_vars_param_x0_cb_sw_f_t(
        NetworkData( dynamics_case()... );
        init_pf = init_pf,
        with_δ_ed_eq = with_δ_ed_eq )
    
end



########################################################

function get_a_model_integrated_u_ur_ui_Idx_in_state(
    netd ;
    Idxs_type =
        :Idxs_im,
    no_control_device =
        false )

    if Idxs_type == :Idxs_hybrid
        
        #------------------------------------------

        slack_ur_ui_Idx_in_state =
            get_components_slack_ur_ui_Idx_in_state(
                    netd.nodes )

        non_slack_ur_ui_Idx_in_state =
            get_components_no_slack_ur_ui_Idx_in_state(
                netd.nodes)

        ur_ui_Idx_in_state =
            get_components_ur_ui_Idx_in_state(
                netd.nodes )

        nodes_u_Idx =
            get_components_ur_ui_Idx_in_state(
                netd.nodes )

       nodes_u_Idx_in_ranges =
           get_nodes_u_Idx_in_ranges(

               nodes_u_Idx )
        
        #------------------------------------------
        
    elseif Idxs_type == :Idxs_im

        dict_sys_to_model_Idx =
            get_net_to_im_indices_dict
        
        Idx_converstion_fun =
            net_to_im_model_indices

        dict_conv =
            dict_sys_to_model_Idx( netd )

        #------------------------------------------

        slack_ur_ui_Idx_in_state =
            Idx_converstion_fun(
            get_components_slack_ur_ui_Idx_in_state(
                netd.nodes ), dict_conv  )

        non_slack_ur_ui_Idx_in_state =
            Idx_converstion_fun(

            get_components_no_slack_ur_ui_Idx_in_state(
                netd.nodes), dict_conv  )

        ur_ui_Idx_in_state = Idx_converstion_fun(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ), dict_conv  )

        nodes_u_Idx = Idx_converstion_fun(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ), dict_conv  )
        
       nodes_u_Idx_in_ranges =
           get_nodes_u_Idx_in_ranges(

               nodes_u_Idx )

        
    elseif Idxs_type == :Idxs_industrial

        if no_control_device == false

            dict_sys_to_model_Idx =
                get_net_to_industrial_model_indices_dict
            
            Idx_converstion_fun =
                net_to_industrial_model_indices

            dict_conv =
                dict_sys_to_model_Idx(
                    netd;
                    no_control_device = false )

            #------------------------------------------

            slack_ur_ui_Idx_in_state =
                Idx_converstion_fun(
                get_components_slack_ur_ui_Idx_in_state(
                    netd.nodes ), dict_conv  )

            non_slack_ur_ui_Idx_in_state =
                Idx_converstion_fun(
                
                get_components_no_slack_ur_ui_Idx_in_state(
                    netd.nodes), dict_conv  )

            ur_ui_Idx_in_state = Idx_converstion_fun(
                get_components_ur_ui_Idx_in_state(
                    netd.nodes ), dict_conv  )

            nodes_u_Idx = Idx_converstion_fun(
                get_components_ur_ui_Idx_in_state(
                    netd.nodes ), dict_conv  )

           nodes_u_Idx_in_ranges =
               get_nodes_u_Idx_in_ranges(

                   nodes_u_Idx )
            
            #----------------------------------------
            
        else

            dict_sys_to_model_Idx =
                get_net_to_industrial_model_indices_dict
            
            Idx_converstion_fun =
                net_to_industrial_model_indices

            dict_conv =
                dict_sys_to_model_Idx(
                    netd;
                    no_control_device = true )

            #----------------------------------------

            slack_ur_ui_Idx_in_state =
                Idx_converstion_fun(
                get_components_slack_ur_ui_Idx_in_state(
                    netd.nodes ), dict_conv )

            non_slack_ur_ui_Idx_in_state =
                Idx_converstion_fun(
                
                get_components_no_slack_ur_ui_Idx_in_state(
                    netd.nodes), dict_conv  )

            ur_ui_Idx_in_state = Idx_converstion_fun(
                get_components_ur_ui_Idx_in_state(
                    netd.nodes ), dict_conv  )

            nodes_u_Idx = Idx_converstion_fun(
                get_components_ur_ui_Idx_in_state(
                    netd.nodes ), dict_conv  )

           nodes_u_Idx_in_ranges =
               get_nodes_u_Idx_in_ranges(

                   nodes_u_Idx )
            
        end
        
    else
        nothing
    end
    
    return (; slack_ur_ui_Idx_in_state, 
            non_slack_ur_ui_Idx_in_state, 
            ur_ui_Idx_in_state, 
            nodes_u_Idx,
            nodes_u_Idx_in_ranges)

end



function get_gens_vh_slack_θh_para(
    nodes_collection )
    
    if isa(nodes_collection,  OrderedDict)
        comp_collection_values =
            collect(values(nodes_collection))

        slack_gens_vh =
            [ comp.Gen.vh for comp in
                 comp_collection_values
                 if comp.Bus_type == :Generator &&
                     comp.isa_slack == true ]


        non_slack_gens_vh =
            [ comp.Gen.vh for comp in
                 comp_collection_values
                 if comp.Bus_type == :Generator &&
                     comp.isa_slack != true ]
        
        slack_gens_θh =
            [ comp.Gen.θh for comp in
                 comp_collection_values
                 if comp.Bus_type == :Generator &&
                     comp.isa_slack == true ]


        gens_vh =
            [ comp.Gen.vh for comp in
                 comp_collection_values
                 if comp.Bus_type == :Generator ]
        
        return (; slack_gens_vh,
                slack_gens_θh,
                gens_vh,
                non_slack_gens_vh ) 
        
    elseif isa(nodes_collection, Union{Array, Vector})


        slack_gens_vh =
            [ comp.Gen.vh for comp in
                 nodes_collection
                 if comp.Bus_type == :Generator &&
                     comp.isa_slack == true ]


        non_slack_gens_vh =
            [ comp.Gen.vh for comp in
                 nodes_collection
                 if comp.Bus_type == :Generator &&
                     comp.isa_slack != true ]
        
        slack_gens_θh =
            [ comp.Gen.θh for comp in
                 nodes_collection
                 if comp.Bus_type == :Generator &&
                     comp.isa_slack == true ]


        gens_vh =
            [ comp.Gen.vh for comp in
                 nodes_collection
                 if comp.Bus_type == :Generator ]
        
        return (; slack_gens_vh,
                slack_gens_θh,
                gens_vh,
                non_slack_gens_vh ) 
    else
        
        comp_collection_values =
            collect(values(nodes_collection))

        slack_gens_vh =
            [ comp.Gen.vh for comp in
                 comp_collection_values
                 if comp.Bus_type == :Generator &&
                     comp.isa_slack == true ]


        non_slack_gens_vh =
            [ comp.Gen.vh for comp in
                 comp_collection_values
                 if comp.Bus_type == :Generator &&
                     comp.isa_slack != true ]
        
        slack_gens_θh =
            [ comp.Gen.θh for comp in
                 comp_collection_values
                 if comp.Bus_type == :Generator &&
                     comp.isa_slack == true ]


        gens_vh =
            [ comp.Gen.vh for comp in
                 comp_collection_values
                 if comp.Bus_type == :Generator ]
        
        return (; slack_gens_vh,
                slack_gens_θh,
                gens_vh,
                non_slack_gens_vh ) 
        
    end

end



function get_net_nodes_type_idxs(
    netd)

    slack_bus_idx =
        get_components_slack_Idx( netd.nodes )
    
    gens_idx = first.(
        get_generators_Idx_and_vh(
            netd.nodes ))        

    # ---------------------------------------------
    
     slack_gens_nodes_idx =
         get_slack_gens_nodes_idx(
             netd.nodes )

     non_slack_gens_nodes_idx =
         get_non_slack_gens_nodes_idx(
             netd.nodes )

     gens_nodes_idx =
         get_gens_nodes_idx( netd.nodes )

     gens_nodes_with_loc_loads_idx =
         get_gens_nodes_with_loc_loads_idx(
             netd.nodes )
    
    loc_load_exist =
        gens_nodes_with_loc_loads_idx == [] ?
        false : true

     load_nodes_idx =
         get_load_nodes_idx( netd.nodes )

     transmission_nodes_idx =
         get_transmission_nodes_idx( netd.nodes )

     non_gens_nodes_idx =
         get_non_gens_nodes_idx( netd.nodes )

     load_trans_nodes_idx =
         get_load_trans_nodes_Idx( netd.nodes )

     all_nodes_idx =
         get_all_nodes_idx( netd.nodes  )

     non_slack_gens_and_non_gens_idx =
         setdiff(all_nodes_idx, slack_gens_nodes_idx)

    # ------------------------------------------------

    return (; slack_bus_idx,
            gens_idx,
            slack_gens_nodes_idx,
            non_slack_gens_nodes_idx,
            gens_nodes_idx,
            gens_nodes_with_loc_loads_idx,
            loc_load_exist,
            load_nodes_idx,
            transmission_nodes_idx,
            non_gens_nodes_idx,
            load_trans_nodes_idx,
            all_nodes_idx,
            non_slack_gens_and_non_gens_idx)
        
end



function get_dict_net_to_streamlined_idx(netd)


    (; slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx) =
         get_net_nodes_type_idxs( netd)

    
    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx

    #------------------------------------------

    streamlined_idx = [ slack_gens_nodes_idx,
                        non_slack_gens_nodes_idx,
                        gens_nodes_idx,
                        non_gens_nodes_idx,
                        load_nodes_idx,
                        all_nodes_idx ]

    vec_dict_net_to_streamlined_idx =
        [ OrderedDict{Int64,Int64}(
            net_idx =>idx
            for (net_idx, idx) in
                zip( a_net_group_idxs,
                     collect(1:length(
                         a_net_group_idxs ))))
      for a_net_group_idxs in
          streamlined_idx ]

    (dict_n2s_slack_gens_idx,
     dict_n2s_non_slack_gens_idx,
     dict_n2s_gens_idx,
     dict_n2s_non_gens_idx,
     dict_n2s_load_idx,
     dict_n2s_all_nodes_idx) =
         vec_dict_net_to_streamlined_idx

    """
    gens_with_loc_load and transmission nodes
    need special attention in cases
    they do not exist

    """        
    if loc_load_exist == true
        dict_n2s_gens_with_loc_load_idxs =
        OrderedDict{Int64, Int64}(
            net_idx =>idx
            for (net_idx, idx) in zip(
                gens_with_loc_load_idx,
                collect(1:length(
                    gens_with_loc_load_idx ))))
    else
        dict_n2s_gens_with_loc_load_idxs = nothing
    end

    
    
    if transmission_nodes_idx == []
        
        dict_n2s_transmission_idxs =
        OrderedDict{Int64, Int64}(
            net_idx =>idx
            for (net_idx, idx) in zip(
                transmission_nodes_idx,
                collect(1:length(
                    transmission_nodes_idx ))))
    else
        dict_n2s_transmission_idxs = nothing
    end

    dicts_net_to_streamlined_idx =
        (; dict_n2s_slack_gens_idx,
         dict_n2s_non_slack_gens_idx,
         dict_n2s_gens_idx,
         dict_n2s_non_gens_idx,
         dict_n2s_load_idx,
         dict_n2s_gens_with_loc_load_idxs,
         dict_n2s_transmission_idxs,
         dict_n2s_all_nodes_idx )
    
    
    return dicts_net_to_streamlined_idx

end


function get_dict_n2s_streamlined_idx(netd)

    (;
     dict_n2s_slack_gens_idx,
     dict_n2s_non_slack_gens_idx,
     dict_n2s_gens_idx,
     dict_n2s_non_gens_idx,
     dict_n2s_load_idx,
     dict_n2s_gens_with_loc_load_idxs,
     dict_n2s_transmission_idxs,
     dict_n2s_all_nodes_idx ) =
         get_dict_net_to_streamlined_idx(netd)


    n2s_slack_gens_idx =
        dict_n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx  =
        dict_n2s_gens_idx
    
    n2s_non_gens_idx  =
        dict_n2s_non_gens_idx
    
    n2s_load_idx =
        dict_n2s_load_idx
    
    n2s_gens_with_loc_load_idxs  =
        dict_n2s_gens_with_loc_load_idxs
    
    n2s_transmission_idxs =
        dict_n2s_transmission_idxs
    
    n2s_all_nodes_idx =
        dict_n2s_all_nodes_idx           
    
    return (;
            n2s_slack_gens_idx,         
            n2s_non_slack_gens_idx,     
            n2s_gens_idx,               
            n2s_non_gens_idx,           
            n2s_load_idx,               
            n2s_gens_with_loc_load_idxs,
            n2s_transmission_idxs,      
            n2s_all_nodes_idx )          
    

end


function get_a_model_integrated_combined_dyn_init(
    netd;
    pf_alg  =
        NewtonRaphson(),
    Idxs_type =
        :Idxs_im  )
    
    """
    #-----------------------------------------------

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer


    pf_alg  = NewtonRaphson()

    Idxs_type = :Idxs_im # :Idxs_hybrid, :Idxs_im, :Idxs_industrial

    #-----------------------------------------------
    
    netd =
        NetworkData( dynamics_case()... )

    """
    
    #-----------------------------------------------
        
    loc_load_exist =
        loc_load_exist_bool( netd )

    #-----------------------------------------------
    
    gens_nodes_collection =
        get_gens_nodes(
            netd.nodes )

    #-----------------------------------------------
        
    nodes_cb_sw =
        get_nodes_cb_sw(netd.nodes)

    #-----------------------------------------------

    gens_nodes_ra_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            netd.nodes_param,
            netd.nodes;
            param_list =
                [ :ra, :X_d_dash, :X_q_dash ],
            gens_view_only = true )
    
    #-----------------------------------------------
    
   im_indices_and_conversion_dict =
         get_im_indices_and_conversion_dict(
             netd  )

    indices_and_conversion_dict =
        im_indices_and_conversion_dict

    im_vars_Idx_in_state =
        indices_and_conversion_dict.im_vars_Idx_in_state

    # (;im_vars_indices_in_system,
    #   pure_states_Idx_in_system,
    #   im_algebraic_vars_Idx_in_system,
    #   ur_ui_Idx_in_system,
    #   im_vars_and_ur_ui_Idx_in_system,
    #   im_vars_Idx_in_state,
    #   nodes_ur_ui_Idx_in_state,
    #   im_state_Idx,
    #   each_gens_im_vars_Idx_in_state,
    #   net_to_im_idx_conversion_dict)  =
    #       im_indices_and_conversion_dict

    #------------------------------------------

    # get_im_model_pf_idx_and_Idx(netd)

    #----------------------------------------
    """
    
    nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
        get_im_δ_ω_ed_dash_eq_dash_Idx(netd )

    """
    
    #----------------------------------------

    # models_gens_nodes_some_vars_Idxs
    
    models_types_gens_δ_ω_ed_dash_eq_dash_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )

    nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_hybrid

    nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_industrial

    nodes_δ_ω_ed_dash_eq_dash_Idxs_im =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_im
    
    nodes_δ_ω_ed_dash_eq_dash_Idxs =
        nodes_δ_ω_ed_dash_eq_dash_Idxs_im

    #------------------------------------------
    #------------------------------------------
    
    integrated_pf_and_init =
         get_a_model_integrated_sta_powerflow_and_init(
             netd;
             pf_alg  =
                 pf_alg)

    #------------------------------------------

    named_tup_pf_result =
        integrated_pf_and_init.named_tup_pf_result

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init

    pf_init_dict = named_tup_pf_result.pf_init_dict
    
    #------------------------------------------
    
    im_state =
        im_model_init_operationpoint(
            netd, bus_dict_init  )

    only_gen  =false
    
    industrial_state =
        industrial_model_init_operationpoint(
            netd, bus_dict_init; pure = :pure,
            no_control_device = only_gen )

    ext_state =
        external_init_operationpoint(
            netd, bus_dict_init,
            branch_dict_init )

    hybrid_state =
        init_operationpoint(netd, pf_init_dict)
    
    #------------------------------------------
    #------------------------------------------
    
    state = im_state

    #----------------------------------------

    # nodes_state_Idx and im_vars_Idx_in_state 

    sim_state_x0 =
        state[ im_vars_Idx_in_state ]

    #----------------------------------------
    
    im_vars_view_in_state =
        get_im_vars_view_in_state(
            state,
            im_vars_Idx_in_state )
    
    #----------------------------------------
    #----------------------------------------
    
    # im_state =
    #     integrated_pf_and_init.im_state


    # industrial_state =
    #     integrated_pf_and_init.industrial_state


    # ext_state =
    #     integrated_pf_and_init.ext_state


    # hybrid_state =
    #     integrated_pf_and_init.hybrid_state
    
    # state = im_state
        
    #----------------------------------------
    
    nodes_name =
        integrated_pf_and_init.nodes_name

    branches_name =
        integrated_pf_and_init.branches_name

    gens_nodes_idx =
        integrated_pf_and_init.gens_nodes_idx
    
    #----------------------------------------
    
    vh = named_tup_pf_result.Vm
    
    θh = named_tup_pf_result.Vθ
    
    #----------------------------------------

    full_vh_θh = [vh; θh]
    
    #----------------------------------------

    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]

    gens_vh_θh =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(gens_vh , gens_θh)]
    
    # gens_vh_θh = get_gens_vh_θh_post_pf(
    #     gens_nodes_collection, bus_dict_init )
    
    #----------------------------------------

    # gens_vh = first.( gens_vh_θh )
    
    # gens_θh = last.( gens_vh_θh ) 
    
    #----------------------------------------    
    
    gens_nodes_ωs_ωref0_vref0_porder0 =
        get_gens_nodes_ωs_ωref0_vref0_porder0(
        gens_nodes_collection,
        bus_dict_init )


    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )
    
    #----------------------------------------

    gen_nodes_δ_ω_ed_dash_eq_dash =
        get_gen_im_nodes_ω_ed_dash_eq_dash(
            state,
            nodes_δ_ω_ed_dash_eq_dash_Idxs )
        
    #----------------------------------------
    
    gens_dynamic_id_iq_pg_vh =
        get_gens_dynamic_id_iq_pg_vh_by_vhθh(
        gens_vh_θh,
        gen_nodes_δ_ω_ed_dash_eq_dash,
            gens_nodes_ra_Xd_dash_Xq_dash )
        
    #------------------------------------------
    #  label syms and mass_matrix
    #------------------------------------------
    
    net_class_names =
           make_case_buses_names(
               netd.nodes )

    # (; network_bus_names,
    #  non_gens_bus_names,
    #  gens_bus_names,
    #  Loads_bus_names,
    #  Trans_bus_names) =
    #      net_class_names

    #------------------------------------------

    im_net_states_and_var_labels =
         generate_im_model_labels(
             ;nodes =
                 netd.nodes )

   # (; net_bus_volts_labels,
   #  gens_nodes_pure_states_labels,
   #  gens_nodes_im_algebraic_vars_labels,
   #  gens_nodes_im_vars_labels) =
   #      im_net_states_and_var_labels
    
 
    #------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes =
                netd.nodes )

    (; im_model_ode_mass_matrix,
     im_model_pf_mass_matrix ) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes =
                 netd.nodes )

    (; gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes =
                 netd.nodes  )

    im_ode_sym =
        gens_nodes_im_vars_labels

    im_ode_mass_matrix =
        im_model_ode_mass_matrix
        
    #-------------------------------------------
    #-------------------------------------------

    para_net_names_labels_syms =
        (;
         net_class_names,
         im_net_states_and_var_labels,
         im_ode_sym )
        
    #----------------------------------------    
    #----------------------------------------
    
    (Ax_Bx_Cx_views,
     Ax_Bx_Cx_matrix,
     states_and_mat_Idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(
            gens_nodes_collection )

    (;vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views )  =
        Ax_Bx_Cx_views
    
    (;Ax_matrix,
     Bx_matrix,
     Cx_matrix ) =
         Ax_Bx_Cx_matrix

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ωs_ω_ref_v_ref_p_order_idxs)  =
         states_and_mat_Idxs
    
    # sys_states_Idxs_and_mat_Idxs =
    #     (; nodes_state_Idx,
    #      Bx_idxs,
    #      Cx_idxs,
    #      id_iq_ph_vh_idxs,
    #      ωs_ω_ref_v_ref_p_order_idxs )
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views ),
        gens_nodes_collection )

    im_plants_system_matrices =
        (;
         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         
         Ax_matrix,
         Bx_matrix,
         Cx_matrix )
    
    #----------------------------------------        

    gens_δ =
        first.( gen_nodes_δ_ω_ed_dash_eq_dash )

    gens_ed_dash =
        third.( gen_nodes_δ_ω_ed_dash_eq_dash )

    gens_eq_dash =
        fourth.( gen_nodes_δ_ω_ed_dash_eq_dash )
    
    # -------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )        

    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
        for (a_vh,
             a_θh,
             a_δ,
             a_ed_dash,
             a_eq_dash,
             a_ra,
             a_X_d_dash,
             a_X_q_dash ) in
            zip( gens_vh,
                 gens_θh,
                 gens_δ,
                 gens_ed_dash,
                 gens_eq_dash,
                 gens_ra,
                 gens_Xd_dash,
                 gens_Xq_dash ) ]

    gens_i_d_0 =
        first.( gens_id_iq )
    
    gens_i_q_0 =
        last.( gens_id_iq )
    
    #----------------------------------------

    gens_vd = [
        a_ed_dash +
            a_Xq_dash * a_iq -
            a_ra * a_id
        for ( a_ed_dash,
              a_Xq_dash,
              a_ra,
              a_id,
              a_iq ) in
            zip(gens_ed_dash,
                gens_Xq_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]

    gens_vq = [
        a_eq_dash - a_Xd_dash * a_id - a_ra * a_id
        for ( a_eq_dash, a_Xd_dash, a_ra, a_id, a_iq ) in
            zip(gens_eq_dash,
                gens_Xd_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]

    #----------------------------------------

    gens_ph = [ a_vd * a_id + a_vq * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    gens_qh = [ a_vq * a_id - a_vd * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]
    
    #----------------------------------------

    # a_id * a_vh * sin(a_δ - a_θh) + a_iq * a_vh * cos(a_δ - a_θh)
    
    gens_ph_by_vh_θh_δ_id_iq = [
        get_a_gen_active_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    # a_id * a_vh * cos(a_δ - a_θh) - a_iq * a_vh * sin(a_δ - a_θh)
    
    gens_qh_by_vh_θh_δ_id_iq = [
        get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]

    #----------------------------------------
    
    S_gens = [
        a_vh * exp(im * (a_θh - a_δ + π/2)) * ( a_id - im * a_iq )
        for ( a_vh, a_θh, a_δ, a_id, a_iq) in
            zip(gens_vh,
                gens_θh,
                gens_δ,
                gens_i_d_0,
                gens_i_q_0) ]
    
    #----------------------------------------
    
    gens_Pei = [
        ed_dash * i_d_0 + eq_dash * i_q_0 +
            (Xq_dash - Xd_dash ) * i_d_0 *  i_q_0
        for (ed_dash,eq_dash,Xd_dash,Xq_dash,i_d_0,i_q_0 ) in
            zip(gens_ed_dash, gens_eq_dash,
                gens_Xd_dash,
                gens_Xq_dash,
                gens_i_d_0, gens_i_q_0 ) ]
    
    #----------------------------------------
    
    intg_vh_θh_id_iq =
        [full_vh_θh; gens_i_d_0; gens_i_q_0]
        
    #----------------------------------------        
    
    (; Ynet,
     nodes_idx_with_adjacent_nodes_idx ) =
         get_Ynet_and_nodes_idx_wt_adjacent_nodes_idx(
             netd )

    nodes_node_idx_and_incident_edges_other_node_idx =
        nodes_idx_with_adjacent_nodes_idx

    dyn_pf_fun_kwd_net_para = 
        (; 
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx )

    #----------------------------------------    
    #----------------------------------------
    
    integrated_pf_vars_and_para_idx =
        get_a_model_integrated_pf_vars_and_para_idx(
            netd;
            Idxs_type =
                Idxs_type,
            no_control_device =
                false )

    #----------------------------------------

    nodes_types_idxs =
        integrated_pf_vars_and_para_idx.nodes_types_idxs

    n2s_idxs =
        integrated_pf_vars_and_para_idx.n2s_idxs

    full_types_Idxs_etc =
        integrated_pf_vars_and_para_idx.full_types_Idxs_etc

    intg_types_Idxs_etc =
        integrated_pf_vars_and_para_idx.intg_types_Idxs_etc

    u_ur_ui_Idx_in_state =
        integrated_pf_vars_and_para_idx.u_ur_ui_Idx_in_state

   (; slack_ur_ui_Idx_in_state, 
    non_slack_ur_ui_Idx_in_state, 
    ur_ui_Idx_in_state, 
    nodes_u_Idx,
    nodes_u_Idx_in_ranges ) =
        u_ur_ui_Idx_in_state

    #----------------------------------------

    # integrated_u_ur_ui_Idx_in_state =
    #     get_a_model_integrated_u_ur_ui_Idx_in_state(
    #         netd;
    #         Idxs_type =
    #             Idxs_type,
    #         no_control_device =
    #             false)

    # (; slack_ur_ui_Idx_in_state, 
    #  non_slack_ur_ui_Idx_in_state, 
    #  ur_ui_Idx_in_state, 
    #  nodes_u_Idx, nodes_u_Idx_in_ranges ) =
    #      integrated_u_ur_ui_Idx_in_state

    # nodes_u_Idx_in_ranges =
    #     get_nodes_u_Idx_in_ranges(

    #         nodes_u_Idx )

    # u_ur_ui_Idx_in_state =
    #     (; slack_ur_ui_Idx_in_state, 
    #      non_slack_ur_ui_Idx_in_state, 
    #      ur_ui_Idx_in_state, 
    #      nodes_u_Idx,
    #      nodes_u_Idx_in_ranges)
    
    #----------------------------------------

    slack_gens_nodes_idx =
        nodes_types_idxs.slack_gens_nodes_idx
    
    non_slack_gens_nodes_idx =
        nodes_types_idxs.non_slack_gens_nodes_idx
    
    gens_nodes_idx =
        nodes_types_idxs.gens_nodes_idx
    
    non_gens_nodes_idx =
        nodes_types_idxs.non_gens_nodes_idx
    
    gens_with_loc_load_idx =
        nodes_types_idxs.gens_with_loc_load_idx
    
    all_nodes_idx = nodes_types_idxs.all_nodes_idx
    
    #----------------------------------------

    gens_nodes_with_loc_loads_idx =
        gens_with_loc_load_idx
    
    dyn_pf_fun_kwd_net_idxs  =  
        (;
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx
         )

    #----------------------------------------
    
    n2s_slack_gens_idx =
        n2s_idxs.n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        n2s_idxs.n2s_non_slack_gens_idx
    
    n2s_gens_idx = n2s_idxs.n2s_gens_idx
    
    n2s_non_gens_idx =
        n2s_idxs.n2s_non_gens_idx
    
    n2s_gens_with_loc_load_idxs =
        n2s_idxs.n2s_gens_with_loc_load_idxs
    
    n2s_all_nodes_idx =
        n2s_idxs.n2s_all_nodes_idx
    
    #----------------------------------------

    dyn_pf_fun_kwd_n2s_idxs =
        (;
         n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx )
    
    #----------------------------------------

    full_gens_vh_Idxs =
        full_types_Idxs_etc.full_gens_vh_Idxs
    
    full_gens_θh_Idxs  =
        full_types_Idxs_etc.full_gens_θh_Idxs
    
    full_non_gens_nodes_vh_Idxs  =
        full_types_Idxs_etc.full_non_gens_nodes_vh_Idxs
    
    full_non_gens_nodes_θh_Idxs  =
        full_types_Idxs_etc.full_non_gens_nodes_θh_Idxs

    # ----
   
    full_vars_Idxs = 
        (;full_gens_vh_Idxs,
         full_gens_θh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_non_gens_nodes_θh_Idxs)
    
    #----------------------------------------
    
    intg_gens_vh_Idxs =
        intg_types_Idxs_etc.intg_gens_vh_Idxs
    
    intg_gens_θh_Idxs =
        intg_types_Idxs_etc.intg_gens_θh_Idxs
    
    intg_non_gens_nodes_vh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_vh_Idxs
    
    intg_non_gens_nodes_θh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_θh_Idxs

    intg_gen_id_Idxs =
        intg_types_Idxs_etc.intg_gen_id_Idxs
    
    intg_gen_iq_Idxs =
        intg_types_Idxs_etc.intg_gen_iq_Idxs
    
    # ----
   
    intg_vars_Idxs = 
        (;
         intg_gens_vh_Idxs,
         intg_gens_θh_Idxs,
         intg_non_gens_nodes_vh_Idxs,
         intg_non_gens_nodes_θh_Idxs,
         intg_gen_id_Idxs,
         intg_gen_iq_Idxs )

    intg_kwd_para =
        (; intg_vars_Idxs, )
    
    #----------------------------------------

    """
    no_ll: no local load
    wt_ll: with local load

    The PQ_δ_ed_eq_para_Idxs is designed in such a way
    that the first five indices of no_ll concides with
    the first five indices of wt_ll

    """

    dyn_PQ_δ_ω_ed_dash_eq_Idxs  =
         get_a_model_integrated_pf_PQ_δ_ω_ed_dash_eq_Idxs(
             netd)

    dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs =
        dyn_PQ_δ_ω_ed_dash_eq_Idxs.dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs
    
    (;
     P_gens_dyn_pf_no_ll_para_Idxs,
     Q_gens_dyn_pf_wt_ll_para_Idxs,
     P_non_gens_dyn_pf_wt_ll_para_Idxs,
     Q_non_gens_dyn_pf_wt_ll_para_Idxs,
     δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs,
     P_g_loc_load_dyn_pf_wt_ll_para_Idxs,
     Q_g_loc_load_dyn_pf_wt_ll_para_Idxs) =
         dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs

    P_gens_dyn_para_Idxs =
        P_gens_dyn_pf_no_ll_para_Idxs
    
    Q_gens_dyn_para_Idxs =
        Q_gens_dyn_pf_wt_ll_para_Idxs
    
    P_non_gens_dyn_para_Idxs =
        P_non_gens_dyn_pf_wt_ll_para_Idxs
    
    Q_non_gens_dyn_para_Idxs =
        Q_non_gens_dyn_pf_wt_ll_para_Idxs
    
    δ_ed_eq_pf_dyn_para_Idxs =
        δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs
    
    P_g_loc_load_dyn_para_Idxs =
        P_g_loc_load_dyn_pf_wt_ll_para_Idxs
    
    Q_g_loc_load_dyn_para_Idxs =
        P_g_loc_load_dyn_pf_wt_ll_para_Idxs

    dyn_pf_fun_kwd_nll_para_vars_Idxs = 
        (;
         P_gens_dyn_para_Idxs,
         Q_gens_dyn_para_Idxs,
         P_non_gens_dyn_para_Idxs,
         Q_non_gens_dyn_para_Idxs,
         δ_ed_eq_pf_dyn_para_Idxs
         ) 

    dyn_pf_fun_kwd_wll_para_vars_Idxs = 
        (;
         P_gens_dyn_para_Idxs,
         Q_gens_dyn_para_Idxs,
         P_non_gens_dyn_para_Idxs,
         Q_non_gens_dyn_para_Idxs,
         δ_ed_eq_pf_dyn_para_Idxs,
         P_g_loc_load_dyn_para_Idxs,
         Q_g_loc_load_dyn_para_Idxs                   
         ) 
    
    #----------------------------------------

    δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        get_gens_nodes_some_state_vars_Idxs_in_flattend(
            gens_nodes_collection
            ; some_state_vars =
                [ :δ, :ω, :ed_dash, :eq_dash ] )
    
    #------------------------------------------

    integrated_pf_PQ_param =
        get_a_model_integrated_pf_PQ_param(
            netd )

    (;
     P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist,
     gens_nodes_ra_Xd_dash_Xq_dash ) =
         integrated_pf_PQ_param

    #----------------------------------------   
    #----------------------------------------
    # get_net_nodes_type_idxs( netd)
    # get_dict_net_to_streamlined_idx( netd)
    # get_a_model_integrated_pf_PQ_δ_ω_ed_dash_eq_Idxs
    #----------------------------------------
    #----------------------------------------

    pf_fun_kwd_para =
        (;
         loc_load_exist,
         dyn_pf_fun_kwd_nll_para_vars_Idxs,
         dyn_pf_fun_kwd_wll_para_vars_Idxs,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_net_para,
         δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
         gens_nodes_ra_Xd_dash_Xq_dash )
    
    #----------------------------------------
    #----------------------------------------
    
    # gens_idq_0 =
    #     [  get_pf_dyn_idq(
    #         vh, θh, gen_δ_ω_ed_eq...,
    #         gen_ra_Xd_dash_Xq_dash... )
    #       for ( vh, θh, gen_δ_ω_ed_eq,
    #             gen_ra_Xd_dash_Xq_dash ) in
    #           zip( gens_vh,
    #                gens_θh,
    #                gen_nodes_δ_ω_ed_dash_eq_dash,
    #                gens_nodes_ra_Xd_dash_Xq_dash ) ]

    # gens_i_d_0 = first.( gens_idq_0 )
        
    # gens_i_q_0 = second.( gens_idq_0 )
    
    #----------------------------------------
    
    full_gens_id_iq =
        (;
         gens_i_d_0,
         gens_i_q_0 )

    full_kwd_para =
        (;
         full_vars_Idxs,
         full_gens_id_iq )

    #----------------------------------------
    
    full_dyn_pf_fun_kwd_para =
        (;
         full_kwd_para,
         pf_fun_kwd_para)
    
    #----------------------------------------

    intg_dyn_pf_fun_kwd_para =
        (;
         intg_vars_Idxs,
         pf_fun_kwd_para )

    #----------------------------------------    
    #----------------------------------------
    # dyn_pf_para
    #----------------------------------------
    #----------------------------------------    

    dyn_pf_wt_ll_para =
        [gens_ph...;
         gens_qh...;
         P_non_gens...;
         Q_non_gens...;
         gen_nodes_δ_ω_ed_dash_eq_dash...;
         P_g_loc_load...;
         Q_g_loc_load...]


    dyn_pf_no_ll_para =
        [gens_ph...;
         gens_qh...;
         P_non_gens...;
         Q_non_gens...;
         gen_nodes_δ_ω_ed_dash_eq_dash... ]

    if loc_load_exist == true

        dyn_pf_fun_flat_para =
            dyn_pf_wt_ll_para
        
    else

        dyn_pf_fun_flat_para =
            dyn_pf_no_ll_para
        
    end

    #-----------------------------------------------
    #-----------------------------------------------

    flat_para_idxs =
        (;
         gens_nodes_vh_θh_idx_in_Idx,
         gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         gens_nodes_id_iq_pg_vh_idx_in_Idx
         
         )

    #-----------------------------------------------

    ode_dyn_pf_flat_idx_in_Idx_idxs =
        (;
         gens_nodes_vh_θh_idx_in_Idx,
         gens_nodes_ωs_ωref0_vref0_porder0_idx_in_Idx,
         gens_nodes_id_iq_pg_vh_idx_in_Idx )

    #-----------------------------------------------

    ode_dyn_pf_para_values =
        (;
         gens_vh_θh,
         gens_nodes_ωs_ωref0_vref0_porder0,
         gens_dynamic_id_iq_pg_vh,
         dyn_pf_fun_flat_para )


    #-----------------------------------------------
    
    # gens_sim_fun_gen_para_Idxs =
    #     gens_vh_θh_ωs_ωref0_vref0_porder0_Idxs

    #-----------------------------------------------
    #-----------------------------------------------

    # ode_para_dyn_pf_para

    #----------------------------------------

    # chunk_size = length(sim_state_x0 )

    # stateDiffCache = similar(sim_state_x0)
    
    # stateDiffCache =
    #     DiffCache(
    #         similar( sim_state_x0 ),
    #         chunk_size )
    
    #----------------------------------------

    # Sys_stateDiffCache = DiffCache(
    #     similar( state ) )

    #----------------------------------------

    gens_sim_state_x0 =
        [ sim_state_x0[ idx ]
         for idx in
             nodes_state_Idx ]

    #----------------------------------------

    # gens_stateDiffCache = [
    #     DiffCache( similar( sim_state_x0[ idx]  ) )
    #     for idx in
    #         nodes_state_Idx ]   

    #-----------------------------------------------

    ode_fun_kwd_para_values =
        (;
         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         gens_nodes_ra_Xd_dash_Xq_dash,
         
         # gens_stateDiffCache,
         
         gens_sim_state_x0,
         
         # im_stateDiffCache,
         
         gens_nodes_collection )

    #----------------------------------------

    dyn_pf_fun_var_values =
        (;
         full_vh_θh,
         intg_vh_θh_id_iq) 

    dyn_pf_fun_para_values =
        (;
         dyn_pf_wt_ll_para,
         dyn_pf_no_ll_para,
         dyn_pf_fun_flat_para )

    dyn_pf_fun_kwd_para =
        (;
         full_dyn_pf_fun_kwd_para,
         intg_dyn_pf_fun_kwd_para )

    #----------------------------------------

    # """ 2 and 4 are dimensions of a gen_vh_θh and 
    # gen_δ_ω_ed_dash_eq_dash respectively in flattend
    # gen_vh_θh_δ_ω_ed_dash_eq_dash for
    # a_gen_voltage_terminal_func """

    # a_gen_vtf_vh_θh_δ_ω_ed_dash_eq_dash_Idx =
    #     get_vars_or_paras_Idxs_in_flattend(
    #         [2, 4]
    #         ; dims_given = true )

    # a_gen_vtf_vh_θh_Idx,  a_gen_vtf_δ_ω_ed_dash_eq_dash_Idx =
    #     a_gen_vtf_vh_θh_vh_θh_δ_ω_ed_dash_eq_dash_Idx


    #----------------------------------------
      
    system_paras_and_kwd_and_Idxs =
        (;
         loc_load_exist,
         
         gens_nodes_collection,
         
         gens_nodes_ra_Xd_dash_Xq_dash,

         im_indices_and_conversion_dict,

         im_vars_Idx_in_state,

         models_types_gens_δ_ω_ed_dash_eq_dash_Idxs,            
         nodes_δ_ω_ed_dash_eq_dash_Idxs,                  

         u_ur_ui_Idx_in_state,
         state,         
         sim_state_x0,

         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,

         Ax_matrix,
         Bx_matrix,
         Cx_matrix,

         states_and_mat_Idxs,

         full_vh_θh,
         intg_vh_θh_id_iq,
         
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_net_para,

         dyn_pf_fun_kwd_nll_para_vars_Idxs,
         dyn_pf_fun_kwd_wll_para_vars_Idxs,

         δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,

         intg_vars_Idxs,

         full_vars_Idxs,
         full_gens_id_iq,
         pf_fun_kwd_para,

         gens_id_iq,

         dyn_pf_wt_ll_para,
         dyn_pf_no_ll_para,
         dyn_pf_fun_flat_para,

         flat_para_idxs,
         
         ode_dyn_pf_flat_idx_in_Idx_idxs,
         
         ode_dyn_pf_para_values,
         
         gens_sim_state_x0,


         )
         
    #----------------------------------------
    #----------------------------------------

    return (;
            loc_load_exist,
            nodes_cb_sw,

            im_indices_and_conversion_dict,
            
            models_types_gens_δ_ω_ed_dash_eq_dash_Idxs,
            
            nodes_δ_ω_ed_dash_eq_dash_Idxs,

            integrated_pf_and_init,

            im_state,
            industrial_state,
            ext_state,
            hybrid_state,
            
            state,

            im_vars_view_in_state,

            full_vh_θh,
            gens_vh_θh,
            
            gens_nodes_ωs_ωref0_vref0_porder0,
            gens_nodes_τm_vf,

            gen_nodes_δ_ω_ed_dash_eq_dash,
            gens_dynamic_id_iq_pg_vh,

            net_class_names,
            im_net_states_and_var_labels,

            im_sym,
            im_mass_matrix,
            im_model_ode_mass_matrix,
            im_model_pf_mass_matrix,

            gens_nodes_im_vars_labels,
            net_bus_volts_labels,

            im_ode_sym,
            im_ode_mass_matrix,

            nodes_state_Idx,
            Bx_idxs,
            Cx_idxs,
            id_iq_ph_vh_idxs,
            ωs_ω_ref_v_ref_p_order_idxs,

            vec_Ax_views,
            vec_Bx_views,
            vec_Cx_views,

            Ax_matrix,
            Bx_matrix,
            Cx_matrix,

            gens_δ,
            gens_ed_dash,
            gens_eq_dash,

            gens_ra,
            gens_Xd_dash,
            gens_Xq_dash,

            gens_id_iq,

            gens_vd, gens_vq, gens_ph, gens_qh,

            gens_ph_by_vh_θh_δ_id_iq,
            gens_qh_by_vh_θh_δ_id_iq,

            S_gens, gens_Pei,

            intg_vh_θh_id_iq,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,

            integrated_pf_vars_and_para_idx,

            u_ur_ui_Idx_in_state,
            
            gens_nodes_with_loc_loads_idx,
            gens_with_loc_load_idx,
            
            dyn_pf_fun_kwd_net_idxs,
            dyn_pf_fun_kwd_n2s_idxs,
            
            full_vars_Idxs,
            intg_vars_Idxs,
            
            intg_types_Idxs_etc,

            intg_kwd_para,
            
            dyn_PQ_δ_ω_ed_dash_eq_Idxs,
            dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs,
            dyn_pf_fun_kwd_nll_para_vars_Idxs,
            dyn_pf_fun_kwd_wll_para_vars_Idxs,

            δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed,
            integrated_pf_PQ_param,

            pf_fun_kwd_para,
            full_gens_id_iq,

            full_kwd_para,
            full_dyn_pf_fun_kwd_para,
            intg_dyn_pf_fun_kwd_para,
            
            dyn_pf_wt_ll_para,
            dyn_pf_no_ll_para,
            dyn_pf_fun_flat_para,

            alt_ode_fun_flat_para,
            type2_alt_ode_fun_flat_para,
            ode_fun_flat_para,
            type2_ode_fun_flat_para,
            alt_ode_dyn_pf_flat_sys_para,
            dyn_type2_alt_ode_dyn_pf_flat_sys_para,
            ode_dyn_pf_flat_sys_para,
            type2_ode_dyn_pf_flat_sys_para,

            flat_para_idxs,
            ode_dyn_pf_flat_idx_in_Idx_idxs,
            ode_dyn_pf_para_values,

            gens_sim_state_x0,
            
            # Sys_stateDiffCache,
            # gens_stateDiffCache,
            
            ode_fun_kwd_para_values,
            dyn_pf_fun_var_values,
            dyn_pf_fun_para_values,
            dyn_pf_fun_kwd_para,

            system_paras_and_kwd_and_Idxs )
        
end


"""
nll : non-local load
dyn : dynamic
iip : inplace

"""
function get_nll_dyn_iip_pf_param(netd)

        
    # ----------------------------------------------

    _, state =
        get_pf_in_named_tuple_with_state( netd)

    # ----------------------------------------------

    """
    hybrid model is used to get gens δ_ω_ed_dash_eq_dash
    from  a standard powerflow for the network
        """            
          

    #  models_gens_nodes_some_vars_Idxs
    
    models_types_gens_δ_ω_ed_dash_eq_dash_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )       

    nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid =
        models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_hybrid
    
    #-----------------------------------------------


    # hybrid_δ_ω_ed_dash_eq_dash_view =
    #     get_gen_nodes_δ_ω_ed_dash_eq_dash_views(
    #         state, nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid )
    
    # δ_ω_ed_dash_eq_dash_view = hybrid_δ_ω_ed_dash_eq_dash_view    

    hybrid_δ_ω_ed_dash_eq_dash =
        get_gen_nodes_δ_ω_ed_dash_eq_dash(
            state, nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid )
        
    δ_ω_ed_dash_eq_dash_view = hybrid_δ_ω_ed_dash_eq_dash    
    δ_ω_ed_dash_eq_dash = δ_ω_ed_dash_eq_dash_view
    
    #-----------------------------------------------
    
    loc_load_exist =
        loc_load_exist_bool(netd)
    
    # ---------------------------------------------------

    pf_net_param =
        get_streamlined_powerflow_net_parameters( netd )
    
    pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx

    # 
    
    # ----------------------------------------------------

    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist

    # ----------------------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx

    gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    gens_with_loc_load_idx =
        net_comp_type_idx.gens_with_loc_load_idx

    # ---------------------------------------------

    gens_with_loc_load_idx =
        loc_load_exist == false ?
        [] : 
        net_comp_type_idx.gens_with_loc_load_idx

    all_nodes_idx =
        net_comp_type_idx.all_nodes_idx

    # ---------------------------------------------
        
    nodes_size = sum(
        [length(gens_nodes_idx),
         length(non_gens_nodes_idx) ])
        
    # ---------------------------------------------

    dict_n2s =
        net_comp_type_idx.dicts_net_to_streamlined_idx

    n2s_gens_idx =
        dict_n2s.dict_n2s_gens_idx
    
    n2s_non_gens_idx =
        dict_n2s.dict_n2s_non_gens_idx
    
    n2s_gens_with_loc_load_idxs =
        dict_n2s.dict_n2s_gens_with_loc_load_idxs

    
    n2s_all_nodes_idx =
        dict_n2s.dict_n2s_all_nodes_idx
    
    # --------------------------------------------------
    
    idx_and_Idx =
        pf_and_dyn_idx_and_Idx.hybrid_pf_etc_idx_and_Idx

    red_vh_Idxs =
        idx_and_Idx.red_vh_Idxs

    red_θh_Idxs =
        idx_and_Idx.red_θh_Idxs

    non_slack_gens_θh_idx2Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx
    
    non_slack_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx_in_Idx
    
    non_gens_θh_idx2Idx =
        idx_and_Idx.non_gens_θh_idx2Idx
    
    non_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_gens_θh_idx2Idx_in_Idx
            
    # ---------------------------------------------

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx
    
    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash
    
    # ----------------------------------------------

    NL_pf_para_Idxs = vec_Idx.NL_pf_para_Idxs
    
    #----------------------------------------------------

    NL_pf_Idxs =
        NL_pf_para_Idxs.no_loc_load_with_δ_ed_NL_pf_para_Idxs

    P_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_gens_NL_para_Idxs

    Q_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_gens_NL_para_Idxs

    P_load_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_load_NL_para_Idxs

    Q_load_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_load_NL_para_Idxs

    δ_ed_eq_pf_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_δ_ed_eq_pf_NL_para_Idxs        
    
    # ---------------------------------------------------
        
    δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        NL_pf_para_Idxs.δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq
    
    # ---------------------------------------------------

    pf_NL_P_Q_para =
        pf_net_misc.pf_NL_P_Q_para
    
    P_gens =
        pf_NL_P_Q_para.P_gens_NL_para
    
    Q_gens =
        pf_NL_P_Q_para.Q_gens_NL_para
        
    P_load =
        pf_NL_P_Q_para.P_load_NL_para
    
    Q_load =
        pf_NL_P_Q_para.Q_load_NL_para
    
    # ---------------------------------------------------

    no_loc_load_with_δ_ed_eq_pf_param =
        ComponentArray(
            P_gens = P_gens,
            Q_gens = Q_gens,
            P_load = P_load,
            Q_load = Q_load,       
            δ_ed_eq_pf =
                [δ_ω_ed_dash_eq_dash_view...;])
        
    # ---------------------------------------------------

    PQ_gens_pad =
        [ idx ∈ gens_nodes_idx ?
        [ P_gens[n2s_gens_idx[idx]],
          Q_gens[n2s_gens_idx[idx]] ] :
              [0.0, 0.0]
          for idx in 1:nodes_size ]

    PQ_non_gens_pad =
        [ idx ∈ non_gens_nodes_idx ?
        [ P_load[n2s_non_gens_idx[idx]],
          Q_load[n2s_non_gens_idx[idx]] ] :
              [0.0, 0.0]
          for idx in 1:nodes_size ]

    P_gens_pad = first.(PQ_gens_pad)
    
    Q_gens_pad = second.(PQ_gens_pad)

    P_non_gens_pad = first.(PQ_non_gens_pad)
    
    Q_non_gens_pad = second.(PQ_non_gens_pad)

    #-------------------------------------

    dims_no_loc_load_with_δ_ed_eq_padded_PQ =
        [length(P_gens_pad),
         length(Q_gens_pad),
         length(P_non_gens_pad),
         length(Q_non_gens_pad),
         length( [δ_ω_ed_dash_eq_dash_view...;]  ) ]

    _,_, no_loc_load_with_δ_ed_eq_padded_PQ_Idx =
        create_size_offset_Idx(
            dims_no_loc_load_with_δ_ed_eq_padded_PQ ;
            counter = 0)

    P_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[1]

    Q_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[2]

    P_non_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[3]

    Q_non_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[4]

    δ_ω_ed_dash_eq_Idx_in_flat_pad =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[5]

    #----------------------------------------------------

    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist
    
    # ---------------------------------------------------
    
    x0_vh = [ idx ∈ gens_nodes_idx ?
        gens_vh[n2s_gens_idx[idx]] : 1.0
              for idx in 1:nodes_size ]
    
    x0_θh = [0.0 for idx in 1:nodes_size ]

    uh_0  = x0_vh .* exp.(im * x0_θh)

    vh_θh_0 = vcat( x0_vh, x0_θh )
     
    red_vh_θh_0 = vh_θh_0[ red_vh_θh_idx  ]

    # -----------------------------------------------

    #  named_tuple_with_state = get_pf_in_named_tuple_with_state( netd)
     
    gens_idq_0 =
        [  get_pf_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( uh_0),
                   angle.( uh_0 ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    gens_i_d_0 = first.( gens_idq_0 )
        
    gens_i_q_0 = second.( gens_idq_0 )
    
    dims_gens_idq =
        length.( [ gens_i_d_0, gens_i_q_0 ]  )

     _,_,  gens_idq_Idx =
         create_size_offset_Idx(
             dims_gens_idq ;
             counter = 0)

    gens_id_Idx = gens_idq_Idx[1]
        
    gens_iq_Idx = gens_idq_Idx[2]


    gens_idq_0_flat = [ gens_i_d_0; gens_i_q_0 ]
        
    # -----------------------------------------------

    dims_red_vh_θh_0_idq_0 =
        length.([ red_vh_θh_0, gens_idq_0_flat ])

     _,_, red_vh_θh_0_idq_0_Idx =
        create_size_offset_Idx(
            dims_red_vh_θh_0_idq_0;
            counter = 0)

    red_vh_θh_0_Idx =
        red_vh_θh_0_idq_0_Idx[1]

    idq_0_Idx =
        red_vh_θh_0_idq_0_Idx[2]
    
    # -----------------------------------------------

    red_vh_θh_0_idq = [ red_vh_θh_0; gens_idq_0_flat ]
    
    # ------------------------------------------------
        
    nodes_types_vh_θh_Idxs =
        pf_and_dyn_idx_and_Idx.hybrid_pf_etc_idx_and_Idx
    
    vh_θh_Idxs_set =
        nodes_types_vh_θh_Idxs.full_nodes_types_Idxs_idx2Idx_etc

    full_slack_gens_vh_Idxs =
        vh_θh_Idxs_set.full_slack_gens_vh_Idxs

    full_non_slack_gens_vh_Idxs =
        vh_θh_Idxs_set.full_non_slack_gens_vh_Idxs

    full_non_gens_nodes_vh_Idxs =
        vh_θh_Idxs_set.full_non_gens_nodes_vh_Idxs

    full_slack_gens_θh_Idxs =
        vh_θh_Idxs_set.full_slack_gens_θh_Idxs

    full_non_slack_gens_θh_Idxs =
        vh_θh_Idxs_set.full_non_slack_gens_θh_Idxs
    
    full_non_gens_nodes_θh_Idxs =
        vh_θh_Idxs_set.full_non_gens_nodes_θh_Idxs

    # ------------------------------------------------

    dims_full_vh_θh_0_idq_0 =
        length.([ vh_θh_0, gens_idq_0_flat ])

     _,_, full_vh_θh_0_idq_Idx =
        create_size_offset_Idx(
            dims_full_vh_θh_0_idq_0;
            counter = 0)

    full_vh_θh_Idx = full_vh_θh_0_idq_Idx[1]

    full_idq_Idx   = full_vh_θh_0_idq_Idx[2]
    
    # ------------------------------------------------

    vh_θh_0_idq = [ vh_θh_0; gens_idq_0_flat ]
    
    # ------------------------------------------------

    red_ΔPQ  = similar( red_vh_θh_0 )
    
    full_ΔPQ = similar( vh_θh_0 )

    red_ΔPQ_Δidq  = similar( red_vh_θh_0_idq )
    
    full_ΔPQ_Δidq = similar( vh_θh_0_idq  )
    
    # ----------------------------------------------- 
    # some selected param
    # -----------------------------------------------

    #

    pois =
        (; P_gens,
         Q_gens,
         P_load,
         Q_load,
         δ_ω_ed_dash_eq_dash_view)

    #

    pois_Idxs =
        (; P_gens_NL_para_Idxs,
         Q_gens_NL_para_Idxs,
         P_load_NL_para_Idxs,
         Q_load_NL_para_Idxs,
         δ_ed_eq_pf_NL_para_Idxs )

    #
    gens_Q_from_red_pf_para =
        (; slack_vh,
         non_slack_gens_Idx_and_vh,
         red_vh_Idxs,
         red_vh_θh_0_Idx,
         idq_0_Idx,
         non_gens_θh_idx2Idx,
         gens_idx,
         n2s_gens_idx,
         non_slack_gens_θh_idx2Idx,
         δ_ω_ed_dash_eq_dash,
         ra_Xd_dash_Xq_dash_view,
         gens_id_Idx,
         gens_iq_Idx )

    #
    
    kwd_net_param =
        (; pf_net,
         pf_idx_and_state,
         pf_param_views,
         pf_and_dyn_idx_and_Idx,
         pf_net_misc,         
         no_loc_load_with_δ_ed_eq_padded_PQ_Idx,
         gens_idq_Idx,
         red_vh_θh_0_idq_0_Idx,
         full_vh_θh_0_idq_Idx,
         δ_ω_ed_dash_eq_dash_view )

    #
    
    kwd_net_full_pf_param =
        (; pf_net,
         pf_idx_and_state,
         pf_param_views,
         pf_and_dyn_idx_and_Idx,
         pf_net_misc,         
         no_loc_load_with_δ_ed_eq_padded_PQ_Idx,
         gens_idq_Idx,
         red_vh_θh_0_idq_0_Idx,
         full_vh_θh_0_idq_Idx,
         δ_ω_ed_dash_eq_dash_view,
         Q_gens)

    nll_dyn_iip_pf_param =
        (;red_ΔPQ,
         red_vh_θh_0,
         full_ΔPQ,
         vh_θh_0,
         red_ΔPQ_Δidq,
         red_vh_θh_0_idq,
         full_ΔPQ_Δidq,
         vh_θh_0_idq,
         gens_Q_from_red_pf_para,
         kwd_net_param,
         kwd_net_full_pf_param )
    
    return (; nll_dyn_iip_pf_param, pois, pois_Idxs) 
end



########################################################
# ------------------------------------------------------
# ------------------------------------------------------
########################################################


function get_a_model_integrated_sta_powerflow_and_init(
    netd;
    pf_alg  =
        NewtonRaphson()  )
    
    #-----------------------------------------------
    
    # ftol=1000*eps()
    
    # xtol=1000*eps()
            
    # pf_alg  = NewtonRaphson()

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
    #     case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    # dynamics_case =
    #     case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    
    # netd =
    #     NetworkData( dynamics_case()... )
    
    #-----------------------------------------------
    
    loc_load_exist =
        loc_load_exist_bool(netd)
    
    # ----------------------------------------------

    (; red_types_Idxs_etc,
     PQ_sta_para_Idxs,
     nodes_types_idxs,
     n2s_idxs ) =
         get_a_model_integrated_sta_pf_vars_and_paras_idx(
             netd )


    (;
     red_non_gens_vh_Idxs,
     red_non_slack_gens_θh_Idxs,
     red_non_gens_θh_Idxs,

     red_vh_Idxs,
     red_θh_Idxs,

     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,
     red_non_slack_gens_θh_idx2Idx_in_Idx,            
     red_non_gens_θh_idx2Idx_in_Idx,

     red_dict_θh_idx2Idx,
     red_dict_θh_idx2Idx_in_Idx,

     red_vh_θh_IDX
     ) =
         red_types_Idxs_etc

    
    (; P_gens_sta_para_Idxs,
    Q_gens_sta_para_Idxs,
    P_non_gens_sta_para_Idxs,
    Q_non_gens_sta_para_Idxs,
    P_g_loc_load_sta_para_Idxs,
     Q_g_loc_load_sta_para_Idxs ) =
         PQ_sta_para_Idxs

    
   (;slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_with_loc_load_idx,
    all_nodes_idx) =
        nodes_types_idxs

    
    (;
     n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx) =
        n2s_idxs

    # ----------------------------------------------
    
    (;
     P_gens,
    Q_gens,
    P_non_gens,
    Q_non_gens,
    P_g_loc_load,
    Q_g_loc_load,
     loc_load_exist
     ) = get_a_model_integrated_sta_pf_PQ_param(
         netd)
    
    # ------------------------------------------------

    (; slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh ) =
         get_gens_vh_slack_θh_para(
             netd.nodes)
    
    # -----------------------------------------------

    (; Ynet,
     nodes_idx_with_adjacent_nodes_idx ) =
         get_Ynet_and_nodes_idx_wt_adjacent_nodes_idx(
             netd )
    
    # ------------------------------------------------

    (; edges_Ybr_cal,
     edges_orientation ) =
         get_edges_Ybr_cal_and_edges_orientation(
             netd )

    # ------------------------------------------------
    
    pf_kw_gens_vh_slack_θh_para =
        (; slack_gens_vh,
         slack_gens_θh,

         gens_vh,
         non_slack_gens_vh )
         

    pf_kw_net_para =
        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx )

    
    pf_kw_var_idxs =
        (; red_vh_Idxs,
         red_non_slack_gens_θh_idx2Idx,
         red_non_gens_θh_idx2Idx )

    
    pf_kw_PQ_para_idxs =
        (;
         P_gens_sta_para_Idxs,
         Q_gens_sta_para_Idxs,
         P_non_gens_sta_para_Idxs,
         Q_non_gens_sta_para_Idxs,
         P_g_loc_load_sta_para_Idxs,
         Q_g_loc_load_sta_para_Idxs ) 
          
    #----------------------------------------
    
    pf_kw_nodes_types_idxs =
        (;slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_with_loc_load_idx ,
         all_nodes_idx) 
              
    #----------------------------------------

    
    pf_kw_n2s_idxs =
        (; n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx)

    #----------------------------------------
    
    rev_n2s_non_slack_gens_idx =
        dict_reverse_keys_values_pair(
            n2s_non_slack_gens_idx )
    
    #----------------------------------------

    pf_kw_para = (
        ;loc_load_exist,
        pf_kw_gens_vh_slack_θh_para,
        pf_kw_net_para,
        pf_kw_var_idxs,
        pf_kw_PQ_para_idxs,
        pf_kw_nodes_types_idxs,
        pf_kw_n2s_idxs
                  )
    
    #----------------------------------------

    if loc_load_exist == true

        pf_PQ_param =
            [P_gens;
             Q_gens;
             P_non_gens;
             Q_non_gens;
             P_g_loc_load;
             Q_g_loc_load]

    else

        pf_PQ_param =
            [P_gens;
             Q_gens;
             P_non_gens;
             Q_non_gens]
        
    end
    
    #----------------------------------------

    sta_red_vh_θh_0 =
        [ ones(length(red_vh_Idxs));
          zeros(length(red_θh_Idxs)) ]

    red_ΔPQ = similar(sta_red_vh_θh_0)
    
    #----------------------------------------
    # Powerflow func and prob
    #----------------------------------------
    
    pf_fun_mismatch =
        get_a_model_integrated_pf_sta_ΔPQ_mismatch
    
    pf_sol =  NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( ( g, x, p ) ->
                pf_fun_mismatch(
                    g, x, p;
                    pf_kw_para =
                        pf_kw_para) ),
            sta_red_vh_θh_0,
            pf_PQ_param ),
        pf_alg )


    #----------------------------------------
    # Results    
    #----------------------------------------

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #----------------------------------------
    
    nodes_name    = collect(keys( netd.nodes )) 

    branches_name = collect(keys( netd.edges ) )
    
    #----------------------------------------

    # all_nodes_idx =
    #     get_all_nodes_idx( netd.nodes  )
    
    red_sol_para =
        (; edges_Ybr_cal,
         edges_orientation,
         nodes_name,
         branches_name,
         pf_kw_para,
         P_g_loc_load,
         Q_g_loc_load,
         P_non_gens,
         Q_non_gens)
    
    # ------------------------------------------------

    named_tup_pf_result =
        get_results_from_pf_sta_red_sol_u(
            pf_sol,
            red_sol_para )
    
    # ------------------------------------------------

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init =
        named_tup_pf_result.branch_dict_init

    pf_init_dict =
        named_tup_pf_result.pf_init_dict
    
    # #----------------------------------------
    # # industrial
    # #----------------------------------------
    
    # indust_state =
    #     industrial_model_init_operationpoint(
    #         netd, bus_dict_init; pure = :pure )
    
    # indust_state =
    #     industrial_model_init_operationpoint(
    #     netd, bus_dict_init; pure = :pure,
    #     no_control_device = false )
    
    # #----------------------------------------
    # # Hybrid, etc
    # #----------------------------------------
    
    # hybrid_state = init_operationpoint(
    #     netd, pf_init_dict)
    
    # #----------------------------------------
    
    # external_init_operationpoint(
    #      netd, bus_dict_init,
    #     branch_dict_init )


    # nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid =
    #     models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_hybrid

    # nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial =
    #     models_types_gens_δ_ω_ed_dash_eq_dash_Idxs.Idxs_industrial


    # gen_nodes_δ_ω_ed_dash_eq_dash =
    #     get_gen_nodes_δ_ω_ed_dash_eq_dash(
    #         state,
    #         nodes_δ_ω_ed_dash_eq_dash_Idxs )

    #----------------------------------------
    # im
    #----------------------------------------
    
    # im_state =
    #     im_model_init_operationpoint(
    #         netd, bus_dict_init  )

    # only_gen  =false
    
    # industrial_state =
    #     industrial_model_init_operationpoint(
    #         netd, bus_dict_init; pure = :pure,
    #         no_control_device = only_gen )

    # ext_state =
    #     external_init_operationpoint(
    #         netd, bus_dict_init,
    #         branch_dict_init )

    hybrid_state = state  =
        init_operationpoint(netd, pf_init_dict)
    
    #----------------------------------------        
    #----------------------------------------
    
    pf_fun_mismatch =
        get_a_model_integrated_pf_sta_ΔPQ_mismatch

    pf_func = NonlinearFunction( ( g, x, p ) ->
        pf_fun_mismatch(
            g, x, p;
            pf_kw_para =
                pf_kw_para ) )

    pf_prob = NonlinearProblem(
        pf_func,
        sta_red_vh_θh_0,
        pf_PQ_param ) 
    
    #----------------------------------------

    pf_sol =
        NonlinearSolve.solve(
            pf_prob, pf_alg )

    #----------------------------------------

    sta_df_dx = ForwardDiff.jacobian(
        ( red_ΔPQ, x ) ->
            pf_fun_mismatch(
                red_ΔPQ, x, pf_PQ_param
                ; pf_kw_para =
                    pf_kw_para),
        red_ΔPQ, pf_sol.u )

    
    sta_df_dp = ForwardDiff.jacobian(
        ( red_ΔPQ, p ) ->
            pf_fun_mismatch(
                red_ΔPQ, pf_sol.u, p
                ; pf_kw_para = pf_kw_para),
        red_ΔPQ, pf_PQ_param )


    # sta_dx_dp  = -(svd( sta_df_dx )) \ sta_df_dp

    #----------------------------------------
    
    
    return (;
            state,
            
            # im_state,
            # industrial_state,
            # ext_state,
            # hybrid_state,
            
            named_tup_pf_result,
            
            sta_df_dx,
            sta_df_dp, 
            
            edges_Ybr_cal,
            edges_orientation,

            nodes_name,
            branches_name,
            
            pf_kw_para,

            get_a_model_integrated_pf_sta_ΔPQ_mismatch,
            red_ΔPQ,
            sta_red_vh_θh_0,
            pf_PQ_param,
            gens_nodes_idx)

end


#-----------------------------------------------------
# integrated
#-----------------------------------------------------


function get_integrated_streamlined_powerflow_net_parameters(
    netd;
    Idxs_type =
        :Idxs_hybrid,
    no_control_device =
        false )

    edges_Ybr_cal, Ynet, edges_orientation,  nodes_node_idx_and_incident_edges_other_node_idx = get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------

    Ybus =
        get_Ybus_from_Ynet(
            Ynet, nodes_node_idx_and_incident_edges_other_node_idx)
      
    #------------------------------------------        
    
    slack_bus_idx =
        get_components_slack_Idx(
            netd.nodes )

    slack_vh  =
        get_components_slack_vh(
            netd.nodes )

    gens_Idx_and_vh =
        get_generators_Idx_and_vh(
            netd.nodes )

    non_slack_gens_Idx_and_vh =
        get_non_slack_generators_Idx_and_vh(
            netd.nodes )

    #------------------------------------------        
    
    gens_idx = first.(gens_Idx_and_vh)
    
    gens_vh  = last.(gens_Idx_and_vh)
    
    # ----------------------------------------------------
    
    load_trans_nodes_Idx_and_vlimits =
        get_load_trans_nodes_Idx_and_vlimits( netd.nodes )

    # ----------------------------------------------------

    ra_Xd_Xq_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            netd.nodes_param, netd.nodes; 
            param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash, :X_q_dash ],
            gens_view_only = true  )

    ra_Xd_dash_Xq_dash  =
        get_gens_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d_dash, :X_q_dash  ],
        gens_view_only = true)
    
    # ----------------------------------------------------
    
    ra_Xd_Xq =
        get_gens_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d, :X_q ],
        gens_view_only = true)
        
    P_Q_nodes =
        get_components_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )
    
    P_Q_gens =
        get_gens_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ],
        gens_view_only = true )

    P_Q_non_gens =
        get_streamlined_non_gens_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ],
        non_gens_view_only = true )
    
    P_Q_gens_loc_load =
        get_streamlined_gens_loc_load_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :loc_P, :loc_Q ],
            gens_view_only = true )
    
    #------------------------------------------        

    if Idxs_type == :Idxs_hybrid

        model_pf_idx_and_Idx =
            get_integrated_hybrid_pf_etc_idx_and_Idx(netd)
        
        (;
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state,
         nodes_u_Idx,         
         ur_idx,
         ui_idx,
         ur_ui_idx,
         vh_θh_idx,
         ur_ui_IDX,
         ur_IDX,
         ui_IDX,
         ir_IDX,
         ii_IDX,
         vh_IDX,
         θh_IDX,         
         slack_bus_idx,
         gens_idx,
         red_vh_θh_idx,
         red_vh_idx,
         red_θh_idx,
         red_vh_Idxs,
         red_θh_Idxs,         
         non_slack_gens_θh_idx2Idx,
         non_slack_gens_θh_idx2Idx_in_Idx,
         non_gens_θh_idx2Idx,
         non_gens_θh_idx2Idx_in_Idx,
         loc_load_exist,
         full_nodes_types_Idxs_idx2Idx_etc,
         intg_nodes_types_Idxs_idx2Idx_etc,
         semi_nodes_types_Idxs_idx2Idx_etc ) =
             model_pf_idx_and_Idx             
        
    elseif Idxs_type == :Idxs_industrial

        model_pf_idx_and_Idx  =
            get_integrated_industrial_model_pf_idx_and_Idx(
                netd;
                no_control_device =
                    no_control_device )
        
        (;
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state,
         nodes_u_Idx,
         ur_idx,
         ui_idx,
         ur_ui_idx,
         vh_θh_idx,
         ur_ui_IDX,
         ur_IDX,
         ui_IDX,
         ir_IDX,
         ii_IDX,
         vh_IDX,
         θh_IDX,
         slack_bus_idx,
         gens_idx,
         red_vh_θh_idx,
         red_vh_idx,
         red_θh_idx,
         red_vh_Idxs,
         red_θh_Idxs,
         non_slack_gens_θh_idx2Idx,
         non_slack_gens_θh_idx2Idx_in_Idx,
         non_gens_θh_idx2Idx,
         non_gens_θh_idx2Idx_in_Idx,
         loc_load_exist,
         full_nodes_types_Idxs_idx2Idx_etc,
         intg_nodes_types_Idxs_idx2Idx_etc,
         semi_nodes_types_Idxs_idx2Idx_etc) =
             model_pf_idx_and_Idx 
             
    elseif Idxs_type == :Idxs_im

        model_pf_idx_and_Idx =
            get_integrated_im_model_pf_idx_and_Idx(netd)
        (;
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state,
         nodes_u_Idx,
         ur_idx,
         ui_idx,
         ur_ui_idx,
         vh_θh_idx,
         ur_ui_IDX,
         ur_IDX,
         ui_IDX,
         ir_IDX,
         ii_IDX,
         vh_IDX,
         θh_IDX,
         slack_bus_idx,
         gens_idx,
         red_vh_θh_idx,
         red_vh_idx,
         red_θh_idx,
         red_vh_Idxs,
         red_θh_Idxs,
         non_slack_gens_θh_idx2Idx,
         non_slack_gens_θh_idx2Idx_in_Idx,
         non_gens_θh_idx2Idx,
         non_gens_θh_idx2Idx_in_Idx,
         loc_load_exist,
         full_nodes_types_Idxs_idx2Idx_etc,
         intg_nodes_types_Idxs_idx2Idx_etc,
         semi_nodes_types_Idxs_idx2Idx_etc ) =
             model_pf_idx_and_Idx
             
    else
         nothing
    end


    #------------------------------------------      

    pf_and_dyn_idx_and_Idx =
        get_integrated_pf_and_dyn_idx_and_Idx(netd)
    
    #------------------------------------------      
    #------------------------------------------

    P_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :Q ] )...;]
    
    # P_load_NL_para =
    #     [get_non_gens_nodes_some_param(
    #         netd.nodes
    #         ; some_param = [ :P] )...;]
    
    # Q_load_NL_para =
    #     [get_non_gens_nodes_some_param(
    #         netd.nodes ; some_param = [ :Q ] )...;]

    
    P_non_gens_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_non_gens_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes ; some_param = [ :Q ] )...;]
    
    # δ_ed_eq_pf_NL_para =
    #     [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;]
    
    P_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_P ] )...;]
    
    Q_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_Q  ] )...;]

    pf_NL_P_Q_para = (; P_gens_NL_para,
                   Q_gens_NL_para,
                   P_non_gens_NL_para,
                   Q_non_gens_NL_para,
                   P_g_loc_load_NL_para,
                   Q_g_loc_load_NL_para)
  
    # pf_NL_P_Q_para 
    #------------------------------------------
    
    pf_net_misc = (; pf_NL_P_Q_para,
                   model_pf_idx_and_Idx  )
    
    #------------------------------------------        
    
    pf_net =
        (; Ybus,
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx,
         edges_Ybr_cal,
         edges_orientation )
    
    pf_idx_and_state =
        (; slack_vh,
         gens_vh,
         gens_Idx_and_vh,
         non_slack_gens_Idx_and_vh,
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state )    

    pf_param_views =
        (; ra_Xd_Xq,
         ra_Xd_dash_Xq_dash,
         ra_Xd_Xq_Xd_dash_Xq_dash,
         P_Q_nodes,
         P_Q_gens,
         P_Q_gens_loc_load,
         P_Q_non_gens )

    pf_limits =
        (; load_trans_nodes_Idx_and_vlimits, )

    pf_Idx =
        (; slack_bus_idx,
         nodes_u_Idx,
         gens_idx,
         ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX,
         red_vh_θh_idx,
         ur_idx,
         ui_idx,
         ur_ui_idx )
    
    pf_net_param =
        (; pf_net,
         pf_idx_and_state,
         pf_param_views,
         pf_limits,
         pf_Idx,
         pf_and_dyn_idx_and_Idx,
         pf_net_misc )

    return pf_net_param
end

         
function get_integrated_streamlined_powerflow_net_parameters(
    ; dynamics_case = file )

    # dynamics_case = new_case_IEEE_9_Bus_dynamic_plant_SM_v6_P

    netd  = NetworkData( dynamics_case()... )

    # ----------------------------------------------------
    # ----------------------------------------------------

    return get_integrated_streamlined_powerflow_net_parameters( netd; Idxs_type = :Idxs_hybrid, no_control_device = false )
end

#-----------------------------------------------



function get_a_model_integrated_streamlined_pf_net_para(
    netd;
    Idxs_type =
        :Idxs_hybrid,
    no_control_device =
        false )

    edges_Ybr_cal, Ynet, edges_orientation,  nodes_node_idx_and_incident_edges_other_node_idx = get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------

    Ybus =
        get_Ybus_from_Ynet(
            Ynet, nodes_node_idx_and_incident_edges_other_node_idx)
      
    #------------------------------------------        
    
    slack_bus_idx =
        get_components_slack_Idx(
            netd.nodes )

    slack_vh  =
        get_components_slack_vh(
            netd.nodes )

    gens_Idx_and_vh =
        get_generators_Idx_and_vh(
            netd.nodes )

    non_slack_gens_Idx_and_vh =
        get_non_slack_generators_Idx_and_vh(
            netd.nodes )

    #------------------------------------------        
    
    gens_idx = first.(gens_Idx_and_vh)
    
    gens_vh  = last.(gens_Idx_and_vh)
    
    # ----------------------------------------------------
    
    load_trans_nodes_Idx_and_vlimits =
        get_load_trans_nodes_Idx_and_vlimits( netd.nodes )

    # ----------------------------------------------------

    ra_Xd_Xq_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            netd.nodes_param, netd.nodes; 
            param_list = [
                :ra, :X_d, :X_q,
                :X_d_dash, :X_q_dash ],
            gens_view_only = true  )

    ra_Xd_dash_Xq_dash  =
        get_gens_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d_dash, :X_q_dash  ],
        gens_view_only = true)
    
    # ----------------------------------------------------
    
    ra_Xd_Xq =
        get_gens_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :ra, :X_d, :X_q ],
        gens_view_only = true)
        
    P_Q_nodes =
        get_components_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ] )
    
    P_Q_gens =
        get_gens_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ],
        gens_view_only = true )

    P_Q_non_gens =
        get_streamlined_non_gens_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :P, :Q ],
        non_gens_view_only = true )
    
    P_Q_gens_loc_load =
        get_streamlined_gens_loc_load_params_in_param_values(
            netd.nodes_param, netd.nodes;
            param_list = [ :loc_P, :loc_Q ],
            gens_view_only = true )
    
    #------------------------------------------        

    
    model_pf_idx_and_Idx =
        get_a_model_integrated_pf_idx_and_Idx(
            netd;
            Idxs_type =
                Idxs_type,
            no_control_device =
                no_control_device)



     (; slack_ur_ui_Idx_in_state,
      non_slack_ur_ui_Idx_in_state,
      ur_ui_Idx_in_state,
      nodes_u_Idx,
      ur_idx,
      ui_idx,
      ur_ui_idx,
      vh_θh_idx,
      ur_ui_IDX,
      ur_IDX,
      ui_IDX,
      ir_IDX,
      ii_IDX,
      vh_IDX,
      θh_IDX,
      slack_bus_idx,
      gens_idx,
      red_vh_θh_idx,
      red_vh_idx,
      red_θh_idx,
      red_vh_Idxs,
      red_θh_Idxs,
      non_slack_gens_θh_idx2Idx,
      non_slack_gens_θh_idx2Idx_in_Idx,
      non_gens_θh_idx2Idx,
      non_gens_θh_idx2Idx_in_Idx,
      loc_load_exist,
      full_nodes_types_Idxs_idx2Idx_etc,
      intg_nodes_types_Idxs_idx2Idx_etc,
      semi_nodes_types_Idxs_idx2Idx_etc ) =
          model_pf_idx_and_Idx
    
    #------------------------------------------      

    pf_and_dyn_idx_and_Idx =
        get_integrated_pf_and_dyn_idx_and_Idx(
            netd)
    
    #------------------------------------------      
    #------------------------------------------

    P_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_gens_NL_para =
        [get_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :Q ] )...;]
    
    # P_load_NL_para =
    #     [get_non_gens_nodes_some_param(
    #         netd.nodes
    #         ; some_param = [ :P] )...;]
    
    # Q_load_NL_para =
    #     [get_non_gens_nodes_some_param(
    #         netd.nodes ; some_param = [ :Q ] )...;]

    
    P_non_gens_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes
            ; some_param = [ :P] )...;]
    
    Q_non_gens_NL_para =
        [get_non_gens_nodes_some_param(
            netd.nodes ; some_param = [ :Q ] )...;]
    
    # δ_ed_eq_pf_NL_para =
    #     [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;]
    
    P_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_P ] )...;]
    
    Q_g_loc_load_NL_para =
        [get_gens_nodes_with_loc_loads_some_param(
            netd.nodes ; some_param = [ :loc_Q  ] )...;]

    pf_NL_P_Q_para = (; P_gens_NL_para,
                   Q_gens_NL_para,
                   P_non_gens_NL_para,
                   Q_non_gens_NL_para,
                   P_g_loc_load_NL_para,
                   Q_g_loc_load_NL_para)
  
    # pf_NL_P_Q_para 
    #------------------------------------------
    
    pf_net_misc = (; pf_NL_P_Q_para,
                   model_pf_idx_and_Idx  )
    
    #------------------------------------------        
    
    pf_net =
        (; Ybus,
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx,
         edges_Ybr_cal,
         edges_orientation )
    
    pf_idx_and_state =
        (; slack_vh,
         gens_vh,
         gens_Idx_and_vh,
         non_slack_gens_Idx_and_vh,
         slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state )    

    pf_param_views =
        (; ra_Xd_Xq,
         ra_Xd_dash_Xq_dash,
         ra_Xd_Xq_Xd_dash_Xq_dash,
         P_Q_nodes,
         P_Q_gens,
         P_Q_gens_loc_load,
         P_Q_non_gens )

    pf_limits =
        (; load_trans_nodes_Idx_and_vlimits, )

    pf_Idx =
        (; slack_bus_idx,
         nodes_u_Idx,
         gens_idx,
         ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX,
         red_vh_θh_idx,
         ur_idx,
         ui_idx,
         ur_ui_idx )
    
    pf_net_param =
        (; pf_net,
         pf_idx_and_state,
         pf_param_views,
         pf_limits,
         pf_Idx,
         pf_and_dyn_idx_and_Idx,
         pf_net_misc )

    return pf_net_param
end

         
function get_integrated_streamlined_powerflow_net_parameters(
    ; dynamics_case = file )

    # dynamics_case = new_case_IEEE_9_Bus_dynamic_plant_SM_v6_P

    netd  = NetworkData( dynamics_case()... )

    # ----------------------------------------------------
    # ----------------------------------------------------

    return get_integrated_streamlined_powerflow_net_parameters( netd; Idxs_type = :Idxs_hybrid, no_control_device = false )
end



#-----------------------------------------------

""" wll: with local loads """
function get_c_integrated_wll_dyn_iip_pf_param(
    netd;
    Idxs_type = :Idxs_hybrid ,
    no_control_device = false )

    """

    no_control_device = false

    Idxs_type = :Idxs_im

    Idxs_type = :Idxs_industrial

    Idxs_type = :Idxs_hybrid

    #-----------------------------------------------
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix

    #-----------------------------------------------
    
    netd =
        NetworkData( dynamics_case()... )

    #-----------------------------------------------

    """
    #-----------------------------------------------
    
    loc_load_exist =
        loc_load_exist_bool(netd)

    # ----------------------------------------------

    """
    hybrid model is used to get gens δ_ω_ed_dash_eq_dash
    from  a standard powerflow for the network

    """

    named_tup_pf_result, state =
        get_pf_in_named_tuple_with_state(
            netd)

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    # ----------------------------------------------

    pf_net_param = get_integrated_streamlined_powerflow_net_parameters( netd; Idxs_type = Idxs_type, no_control_device = no_control_device )

    # ----------------------------------------------
    
    models_gens_nodes_some_vars_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )       

    if Idxs_type == :Idxs_hybrid
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs =
            models_gens_nodes_some_vars_Idxs.Idxs_hybrid

        state = state

        δ_ω_ed_dash_eq_dash = 
            get_gen_nodes_δ_ω_ed_dash_eq_dash(
                state,
                nodes_δ_ω_ed_dash_eq_dash_Idxs )
                
        idx_and_Idx =
            get_hybrid_pf_etc_idx_and_Idx(netd)

        
    elseif Idxs_type == :Idxs_industrial
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs =
            models_gens_nodes_some_vars_Idxs.Idxs_industrial

        state = industrial_model_init_operationpoint(
            netd,
            bus_dict_init
            ;pure =
                :pure,
            no_control_device =
                false )
        
        δ_ω_ed_dash_eq_dash = 
            get_gen_nodes_δ_ω_ed_dash_eq_dash(
                state,
                nodes_δ_ω_ed_dash_eq_dash_Idxs )

        
        idx_and_Idx = get_industrial_model_pf_idx_and_Idx(netd; no_control_device = no_control_device)

    elseif Idxs_type == :Idxs_im
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs =
            models_gens_nodes_some_vars_Idxs.Idxs_im

        state = im_model_init_operationpoint(
            netd,
            bus_dict_init )

        δ_ω_ed_dash_eq_dash = 
            get_gen_nodes_δ_ω_ed_dash_eq_dash(
                state,
                nodes_δ_ω_ed_dash_eq_dash_Idxs )

        
        idx_and_Idx =
            get_im_model_pf_idx_and_Idx(netd)
        
    else
        nothing
    end
    
    #--------------------------------------------   
    #--------------------------------------------
    
    _, pf_idx_and_state, pf_param_views,_, pf_Idx, pf_and_dyn_idx_and_Idx, pf_net_misc = pf_net_param

    # Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, _, _, _ = pf_idx_and_state

    ra_Xd_Xq, ra_Xd_dash_Xq_dash, ra_Xd_Xq_Xd_dash_Xq_dash, _, _, _, _ = pf_param_views

    # load_trans_nodes_Idx_and_vlimits = pf_limits

    _,nodes_u_Idx,_, _, _, _, _, red_vh_θh_idx,_, _,_ = pf_Idx
    
    # ---------------------------------------------

    # loc_load_exist =
    #     pf_and_dyn_idx_and_Idx.loc_load_exist

    # ---------------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx


    slack_gens_nodes_idx =
        net_comp_type_idx.slack_gens_nodes_idx

    non_slack_gens_nodes_idx =
        net_comp_type_idx.non_slack_gens_nodes_idx 
    
    gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    gens_with_loc_load_idx =
        net_comp_type_idx.gens_with_loc_load_idx

    load_nodes_idx =
        net_comp_type_idx.load_nodes_idx

    transmission_nodes_idx =
        net_comp_type_idx.transmission_nodes_idx

    load_trans_nodes_idx =
        net_comp_type_idx.load_trans_nodes_idx
        
    # -----------------------------------------

    gens_with_loc_load_idx =
        loc_load_exist == false ?
        [] : 
        net_comp_type_idx.gens_with_loc_load_idx

    all_nodes_idx =
        net_comp_type_idx.all_nodes_idx

    # ----------------------------------------
        
    nodes_size = length(gens_nodes_idx) +
        length(non_gens_nodes_idx)
    
        
    # ----------------------------------------

    dict_n2s =
        net_comp_type_idx.dicts_net_to_streamlined_idx

    n2s_slack_gens_idx =
        dict_n2s.dict_n2s_slack_gens_idx

    n2s_non_slack_gens_idx =
        dict_n2s.dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx =
        dict_n2s.dict_n2s_gens_idx
    
    n2s_non_gens_idx =
        dict_n2s.dict_n2s_non_gens_idx

    n2s_load_idx =
        dict_n2s.dict_n2s_load_idx
    
    n2s_gens_with_loc_load_idxs =
        dict_n2s.dict_n2s_gens_with_loc_load_idxs

    n2s_transmission_idxs =
        dict_n2s.dict_n2s_transmission_idxs
    
    n2s_all_nodes_idx =
        dict_n2s.dict_n2s_all_nodes_idx
    
    # -------------------------------------------

    red_vh_Idxs =
        idx_and_Idx.red_vh_Idxs

    red_θh_Idxs =
        idx_and_Idx.red_θh_Idxs

    non_slack_gens_θh_idx2Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx
    
    non_slack_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx_in_Idx
    
    non_gens_θh_idx2Idx =
        idx_and_Idx.non_gens_θh_idx2Idx
    
    non_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_gens_θh_idx2Idx_in_Idx
            
    # -----------------------------------------

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx
    
    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash
    
    # ----------------------------------------

    NL_pf_para_Idxs = vec_Idx.NL_pf_para_Idxs
    
    #--------------------------------------------

    NL_pf_Idxs =
        NL_pf_para_Idxs.no_loc_load_with_δ_ed_NL_pf_para_Idxs

    P_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_gens_NL_para_Idxs

    Q_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_gens_NL_para_Idxs

    P_non_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_non_gens_NL_para_Idxs

    Q_non_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_non_gens_NL_para_Idxs

    δ_ed_eq_pf_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_δ_ed_eq_pf_NL_para_Idxs        
    
    # ------------------------------------------
        
    δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        NL_pf_para_Idxs.δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq
    
    # ------------------------------------------

    pf_NL_P_Q_para =
        pf_net_misc.pf_NL_P_Q_para
    
    P_gens =
        pf_NL_P_Q_para.P_gens_NL_para
    
    Q_gens =
        pf_NL_P_Q_para.Q_gens_NL_para
        
    P_non_gens =
        pf_NL_P_Q_para.P_non_gens_NL_para
    
    Q_non_gens =
        pf_NL_P_Q_para.Q_non_gens_NL_para
    
    # --------------------------------------------

    no_loc_load_with_δ_ed_eq_pf_param =
        ComponentArray(
            P_gens = P_gens,
            Q_gens = Q_gens,
            P_non_gens = P_non_gens,
            Q_non_gens = Q_non_gens,       
            δ_ed_eq_pf =
                [δ_ω_ed_dash_eq_dash...;])
        
    # -------------------------------------------

    PQ_gens_pad =
        [ idx ∈ gens_nodes_idx ?
        [ P_gens[n2s_gens_idx[idx]],
          Q_gens[n2s_gens_idx[idx]] ] :
              [0.0, 0.0]
          for idx in 1:nodes_size ]

    PQ_non_gens_pad =
        [ idx ∈ non_gens_nodes_idx ?
        [ P_non_gens[n2s_non_gens_idx[idx]],
          Q_non_gens[n2s_non_gens_idx[idx]] ] :
              [0.0, 0.0]
          for idx in 1:nodes_size ]

    P_gens_pad = first.(PQ_gens_pad)
    
    Q_gens_pad = second.(PQ_gens_pad)

    P_non_gens_pad = first.(PQ_non_gens_pad)
    
    Q_non_gens_pad = second.(PQ_non_gens_pad)

    #-------------------------------------

    dims_no_loc_load_with_δ_ed_eq_padded_PQ =
        [length(P_gens_pad),
         length(Q_gens_pad),
         length(P_non_gens_pad),
         length(Q_non_gens_pad),
         length( [δ_ω_ed_dash_eq_dash...;]  ) ]

    _,_, no_loc_load_with_δ_ed_eq_padded_PQ_Idx =
        create_size_offset_Idx(
            dims_no_loc_load_with_δ_ed_eq_padded_PQ ;
            counter = 0)

    P_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[1]

    Q_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[2]

    P_non_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[3]

    Q_non_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[4]

    δ_ω_ed_dash_eq_Idx_in_flat_pad =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[5]

    #---------------------------------------------

    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist
    
    # ---------------------------------------------
    
    # x0_vh = [ idx ∈ gens_nodes_idx ?
    #     gens_vh[n2s_gens_idx[idx]] : 1.0
    #           for idx in 1:nodes_size ]
    
    # x0_θh =
    #     [0.0
    #      for idx in
    #          1:nodes_size ]
    
    # uh_0  =
    #     x0_vh .* exp.(im * x0_θh)
    
    # vh_θh_0 =
    #     vcat( x0_vh, x0_θh )
     
    # red_vh_θh_0 =
    #     vh_θh_0[ red_vh_θh_idx  ]


    uh_0  = x_from_xr_xi.(
            [state[idx]
             for idx in
                 nodes_u_Idx])

    vh_θh_0 = [ abs.( uh_0  ); angle.( uh_0  )]
     
    red_vh_θh = vh_θh_0[ red_vh_θh_idx  ]

    
    # gens_uh_0 = [ uh_0[  idx  ]
    #             for idx  in 1:nodes_size
    #                 if idx ∈ gens_nodes_idx]

    
    gens_uh_0 =  uh_0[ gens_nodes_idx ]
                
    # --------------------------------------------
    
    #  named_tuple_with_state = get_pf_in_named_tuple_with_state( netd)
     
    gens_idq_0 =
        [  get_a_gen_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( gens_uh_0  ),
                   angle.( gens_uh_0  ),
                   δ_ω_ed_dash_eq_dash,
                   ra_Xd_dash_Xq_dash ) ]

    gens_i_d_0 = first.( gens_idq_0 )
        
    gens_i_q_0 = second.( gens_idq_0 )
    
    dims_gens_idq =
        length.( [ gens_i_d_0, gens_i_q_0 ]  )

     _,_,  gens_idq_Idx =
         create_size_offset_Idx(
             dims_gens_idq ;
             counter = 0)

    gens_id_Idx = gens_idq_Idx[1]
        
    gens_iq_Idx = gens_idq_Idx[2]

    gens_idq_0_flat = [ gens_i_d_0 ; gens_i_q_0 ]

    # -------------------------------------------
    
    dims_red_vh_θh_id_iq_CA =
        length.([ red_vh_θh,
                  gens_i_d_0,
                  gens_i_q_0 ])

     _,_, red_vh_θh_id_iq_CA_Idx =
        create_size_offset_Idx(
            dims_red_vh_θh_id_iq_CA;
            counter = 0)

    red_vh_θh_CA_Idx =
        red_vh_θh_id_iq_CA_Idx[1]

    gen_id_CA_Idx =
        red_vh_θh_id_iq_CA_Idx[2]

    gen_iq_CA_Idx =
        red_vh_θh_id_iq_CA_Idx[3]

    
    # -------------------------------------------
    

   """
    t_δ = first.( δ_ω_ed_dash_eq_dash )

    t_gens_idq_0 = x_from_xr_xi.( gens_idq_0 )

    t_gens_idq_net_0 = t_gens_idq_0 .*
        exp.(im * (t_δ  .- pi/2) )

    t_gens_S = gens_uh_0 .* conj.(t_gens_idq_net_0 )

     """
    
    # -----------------------------------------

    dims_red_vh_θh_0_idq_0 =
        length.([ red_vh_θh, gens_idq_0_flat ])

     _,_, red_vh_θh_0_idq_0_Idx =
        create_size_offset_Idx(
            dims_red_vh_θh_0_idq_0;
            counter = 0)

    flat_red_vh_θh_0_Idx =
        red_vh_θh_0_idq_0_Idx[1]

    flat_idq_0_Idx =
        red_vh_θh_0_idq_0_Idx[2]
    
    # -------------------------------------------

    red_vh_θh_0_idq = [ red_vh_θh; gens_idq_0_flat ]
    
    # --------------------------------------------
    
    model_pf_idx_and_Idx =
        pf_net_misc.model_pf_idx_and_Idx
    
    vh_θh_Idxs_set =
        model_pf_idx_and_Idx.full_nodes_types_Idxs_idx2Idx_etc

    full_slack_gens_vh_Idxs =
        vh_θh_Idxs_set.full_slack_gens_vh_Idxs

    full_non_slack_gens_vh_Idxs =
        vh_θh_Idxs_set.full_non_slack_gens_vh_Idxs

    full_non_gens_nodes_vh_Idxs =
        vh_θh_Idxs_set.full_non_gens_nodes_vh_Idxs

    full_slack_gens_θh_Idxs =
        vh_θh_Idxs_set.full_slack_gens_θh_Idxs

    full_non_slack_gens_θh_Idxs =
        vh_θh_Idxs_set.full_non_slack_gens_θh_Idxs
    
    full_non_gens_nodes_θh_Idxs =
        vh_θh_Idxs_set.full_non_gens_nodes_θh_Idxs

    # -------------------------------------------

    dims_full_vh_θh_0_idq_0 =
        length.([ vh_θh_0, gens_idq_0_flat ])

     _,_, full_vh_θh_0_idq_Idx =
        create_size_offset_Idx(
            dims_full_vh_θh_0_idq_0;
            counter = 0)

    full_vh_θh_Idx = full_vh_θh_0_idq_Idx[1]

    full_idq_Idx   = full_vh_θh_0_idq_Idx[2]
    
    # --------------------------------------------

    # vh_θh_0_idq = [ vh_θh_0; gens_idq_0_flat ]
    
    # --------------------------------------------

    red_ΔPQ  = similar( red_vh_θh )
    
    red_ΔPQ_Δidq  = similar( red_vh_θh_0_idq )
    
    # ------------------------------------------- 
    # some selected param
    # -------------------------------------------
    
    net_idxs =
        (; slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         load_nodes_idx,
         transmission_nodes_idx,
         non_gens_nodes_idx,
         load_trans_nodes_idx,
         gens_with_loc_load_idx)
    
    n2s_idxs =
        (; n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_load_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_transmission_idxs,
         n2s_all_nodes_idx)  

    gens_uh_Q_from_red_sol_para =
        (; slack_vh,
         gens_vh,
         non_slack_gens_Idx_and_vh,
         ra_Xd_dash_Xq_dash,         
         red_vh_Idxs,
         red_θh_Idxs,
         red_vh_θh_idx,
         flat_red_vh_θh_0_Idx,
         flat_idq_0_Idx,
         gens_id_Idx,
         gens_iq_Idx,         
         non_gens_θh_idx2Idx,
         non_slack_gens_θh_idx2Idx,         
         net_idxs,
         n2s_idxs  )
    
    idx2Idx =
        (; non_slack_gens_θh_idx2Idx,
          non_slack_gens_θh_idx2Idx_in_Idx,
          non_gens_θh_idx2Idx,
          non_gens_θh_idx2Idx_in_Idx)

    red_vh_θh_0 = deepcopy(red_vh_θh)
    
    red_vh_θh_id_iq_CA_comp =
        (;red_vh_θh_0,
         gens_i_d_0,
         gens_i_q_0 )

    red_vh_θh_id_iq_CA_Idxs =
        (;red_vh_θh_CA_Idx,
         gen_id_CA_Idx,
         gen_iq_CA_Idx,
         red_vh_θh_id_iq_CA_Idx )    
    pois =
        (; P_gens,
         Q_gens,
         P_non_gens,
         Q_non_gens,
         δ_ω_ed_dash_eq_dash)

    pois_Idxs =
        (; P_gens_NL_para_Idxs,
         Q_gens_NL_para_Idxs,
         P_non_gens_NL_para_Idxs,
         Q_non_gens_NL_para_Idxs,
         δ_ed_eq_pf_NL_para_Idxs )
    
    kwd_net_param =
        (; pf_net_param,
         gens_idq_Idx,
         red_vh_θh_0_idq_0_Idx,
         full_vh_θh_0_idq_Idx )
    
    nll_dyn_iip_pf_param =
        (; 
         red_ΔPQ_Δidq,
         red_vh_θh_0_idq,
         gens_uh_Q_from_red_sol_para,
         kwd_net_param
         )
    
    # return (; nll_dyn_iip_pf_param, pois, pois_Idxs)
    
    flat_δ_ω_ed_dash_eq_dash =
        [δ_ω_ed_dash_eq_dash...;]
    
    return (; nll_dyn_iip_pf_param,
            flat_δ_ω_ed_dash_eq_dash,
            pois,
            pois_Idxs,
            idx2Idx,
            red_vh_θh_id_iq_CA_comp,
            red_vh_θh_id_iq_CA_Idxs) 
end



function get_integrated_dyn_pf_param(
    netd;
    Idxs_type =
        :Idxs_hybrid,
    no_control_device =
        false)

    #------------------------------------
    
    """
    
    abstol=1e-12
    
    reltol=1e-12
    
    init_pf = true
    
    pf_alg = NewtonRaphson()
    
    with_δ_ed_eq = false

    :Idxs_industrial, :Idxs_im
    
    Idxs_type = :Idxs_hybrid

    no_control_device = false
    
    #-----------------------------------
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix

    #-----------------------------------
    
    netd =
        NetworkData( dynamics_case()... )

    #-----------------------------------

    get_c_integrated_nll_dyn_iip_pf_param(
    netd;
    Idxs_type = :Idxs_hybrid ,
        no_control_device = false )

    """
    
    #-----------------------------------
    
    (; nll_dyn_iip_pf_param,
     flat_δ_ω_ed_dash_eq_dash,
     pois,
     pois_Idxs,
     idx2Idx,
     red_vh_θh_id_iq_CA_comp,
     red_vh_θh_id_iq_CA_Idxs) =
        get_c_integrated_nll_dyn_iip_pf_param(
            netd;
            Idxs_type = Idxs_type,
            no_control_device = no_control_device)
    
    (; red_ΔPQ_Δidq,
     red_vh_θh_0_idq,
     gens_uh_Q_from_red_sol_para,
     kwd_net_param
     ) =
         nll_dyn_iip_pf_param 
    
    (; pf_net_param,
    gens_idq_Idx,
    red_vh_θh_0_idq_0_Idx,
    full_vh_θh_0_idq_Idx ) =
         kwd_net_param 

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) =
         gens_uh_Q_from_red_sol_para
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx) =
         n2s_idxs

    (; non_slack_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx_in_Idx,
     non_gens_θh_idx2Idx,
     non_gens_θh_idx2Idx_in_Idx) =
         idx2Idx

    (;red_vh_θh_0,
     gens_i_d_0,
     gens_i_q_0 ) =
         red_vh_θh_id_iq_CA_comp

    (;red_vh_θh_CA_Idx,
     gen_id_CA_Idx,
     gen_iq_CA_Idx,
     red_vh_θh_id_iq_CA_Idx ) =
         red_vh_θh_id_iq_CA_Idxs

    
    #------------------------------------

    pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx, pf_net_misc =
        pf_net_param

    #------------------------------------

    _, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, _ = pf_net

    #---------------------------------------

    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist

    #---------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx

    slack_gens_nodes_idx =
        net_comp_type_idx.slack_gens_nodes_idx

    slack_bus_idx =
        slack_gens_nodes_idx

    non_slack_gens_nodes_idx =
        net_comp_type_idx.non_slack_gens_nodes_idx

    gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    gens_with_loc_load_idx =
        loc_load_exist == false ?
        [] : 
        net_comp_type_idx.gens_with_loc_load_idx

    load_nodes_idx =
        net_comp_type_idx.load_nodes_idx

    transmission_nodes_idx =
        net_comp_type_idx.transmission_nodes_idx

    load_trans_nodes_idx =
        net_comp_type_idx.load_trans_nodes_idx

    all_nodes_idx =
        net_comp_type_idx.all_nodes_idx
    
    #-------------------------------------

    nodes_size = sum(
        [length(gens_nodes_idx),
         length(non_gens_nodes_idx) ])

    #------------------------------------

    dict_n2s =
        net_comp_type_idx.dicts_net_to_streamlined_idx

    n2s_slack_gens_idx =
        dict_n2s.dict_n2s_slack_gens_idx

    n2s_non_slack_gens_idx =
        dict_n2s.dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx =
        dict_n2s.dict_n2s_gens_idx

    n2s_non_gens_idx =
        dict_n2s.dict_n2s_non_gens_idx

    n2s_load_idx =
        dict_n2s.dict_n2s_load_idx
    
    n2s_gens_with_loc_load_idxs =
        dict_n2s.dict_n2s_gens_with_loc_load_idxs

    n2s_transmission_idxs =
        dict_n2s.dict_n2s_transmission_idxs
    
    n2s_all_nodes_idx =
        dict_n2s.dict_n2s_all_nodes_idx
    
    #-------------------------------------
            
    idx_and_Idx =
        pf_net_misc.model_pf_idx_and_Idx

    #---------------------------------------
    
    red_vh_Idxs =
        idx_and_Idx.red_vh_Idxs

    red_θh_Idxs =
        idx_and_Idx.red_θh_Idxs

    non_slack_gens_θh_idx2Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx

    non_slack_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx_in_Idx

    non_gens_θh_idx2Idx =
        idx_and_Idx.non_gens_θh_idx2Idx

    non_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_gens_θh_idx2Idx_in_Idx

    #---------------------------------------

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx

    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash

    #--------------------------------------

    NL_pf_para_Idxs = vec_Idx.NL_pf_para_Idxs

    #--------------------------------------

    NL_pf_Idxs =
        NL_pf_para_Idxs.no_loc_load_with_δ_ed_NL_pf_para_Idxs

    P_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_gens_NL_para_Idxs

    Q_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_gens_NL_para_Idxs

    P_non_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_non_gens_NL_para_Idxs

    Q_non_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_non_gens_NL_para_Idxs

    δ_ed_eq_pf_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_δ_ed_eq_pf_NL_para_Idxs        

    #--------------------------------------
    
    δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        NL_pf_para_Idxs.δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq

    #------------------------------------

    id_Idx, iq_Idx = gens_idq_Idx

    flat_red_vh_θh_0_Idx, flat_idq_0_Idx =
        red_vh_θh_0_idq_0_Idx

    #--------------------------------------

    model_idxs =
        pf_net_misc.model_pf_idx_and_Idx
    
    semi_idxs =
        model_idxs.semi_nodes_types_Idxs_idx2Idx_etc

    intg_idxs =
        model_idxs.intg_nodes_types_Idxs_idx2Idx_etc


    full_idxs =
        model_idxs.full_nodes_types_Idxs_idx2Idx_etc
    
   (;
    intg_slack_gens_vh_Idxs,
    intg_non_slack_gens_vh_Idxs,
    intg_non_gens_nodes_vh_Idxs,

    intg_slack_gens_θh_Idxs,
    intg_non_slack_gens_θh_Idxs,
    intg_non_gens_nodes_θh_Idxs,

    intg_gen_id_Idxs,
    intg_gen_iq_Idxs,

    intg_slack_gens_vh_idx2Idx,
    intg_non_slack_gens_vh_idx2Idx,
    intg_non_gens_vh_idx2Idx,

    intg_slack_gens_θh_idx2Idx,
    intg_non_slack_gens_θh_idx2Idx,
    intg_non_gens_θh_idx2Idx,

    intg_dict_vh_idx2Idx,
    intg_dict_θh_idx2Idx,

    intg_slack_gens_vh_idx2Idx_in_Idx,
    intg_non_slack_gens_vh_idx2Idx_in_Idx,
    intg_non_gens_vh_idx2Idx_in_Idx,

    intg_slack_gens_θh_idx2Idx_in_Idx,
    intg_non_slack_gens_θh_idx2Idx_in_Idx,
    intg_non_gens_θh_idx2Idx_in_Idx,

    intg_dict_vh_idx2Idx_in_Idx,
    intg_dict_θh_idx2Idx_in_Idx
    ) = intg_idxs

   (;

    semi_non_slack_gens_vh_Idxs,
    semi_non_gens_nodes_vh_Idxs,

    semi_non_slack_gens_θh_Idxs,
    semi_non_gens_nodes_θh_Idxs,

    semi_gen_id_Idxs,
    semi_gen_iq_Idxs,

    semi_non_slack_gens_vh_idx2Idx,
    semi_non_gens_vh_idx2Idx,

    semi_non_slack_gens_θh_idx2Idx,
    semi_non_gens_θh_idx2Idx,

    semi_dict_vh_idx2Idx,
    semi_dict_θh_idx2Idx,

    semi_non_slack_gens_vh_idx2Idx_in_Idx,
    semi_non_gens_vh_idx2Idx_in_Idx,

    semi_non_slack_gens_θh_idx2Idx_in_Idx,
    semi_non_gens_θh_idx2Idx_in_Idx,

    semi_dict_vh_idx2Idx_in_Idx,
    semi_dict_θh_idx2Idx_in_Idx
    ) =
        semi_idxs

    (;
     full_slack_gens_vh_Idxs,
    full_non_slack_gens_vh_Idxs,
    full_non_gens_nodes_vh_Idxs,
    
    full_slack_gens_θh_Idxs,
    full_non_slack_gens_θh_Idxs,
    full_non_gens_nodes_θh_Idxs,
    
    full_slack_gens_vh_idx2Idx,
    full_non_slack_gens_vh_idx2Idx,
    full_non_gens_vh_idx2Idx,    
    full_slack_gens_θh_idx2Idx,
    full_non_slack_gens_θh_idx2Idx,
    full_non_gens_θh_idx2Idx,
    full_dict_vh_idx2Idx,
    full_dict_θh_idx2Idx,
    full_slack_gens_vh_idx2Idx_in_Idx,
    full_non_slack_gens_vh_idx2Idx_in_Idx,
    full_non_gens_vh_idx2Idx_in_Idx,
    full_slack_gens_θh_idx2Idx_in_Idx,
    full_non_slack_gens_θh_idx2Idx_in_Idx,
    full_non_gens_θh_idx2Idx_in_Idx,    
    full_dict_vh_idx2Idx_in_Idx,
    full_dict_θh_idx2Idx_in_Idx
    ) =
        full_idxs


    #--------------------------------------

    non_slack_gens_vh =
        last.(non_slack_gens_Idx_and_vh)
    
    non_gen_vh =
        red_vh_θh_0[ red_vh_Idxs ]

    non_slack_and_non_gen_θh =
        red_vh_θh_0[ red_θh_Idxs ]
    
    #--------------------------------------    
    
    intg_gens_vh_Idxs =
        first(intg_slack_gens_vh_Idxs):last(intg_non_slack_gens_vh_Idxs)
    
    intg_gens_θh_Idxs =
        first(intg_slack_gens_θh_Idxs):last(intg_non_slack_gens_θh_Idxs)

    
    intg_vh_θh_id_iq =
        [ gens_vh;
          non_gen_vh;
          [0];
          non_slack_and_non_gen_θh;
          gens_i_d_0;
          gens_i_q_0 ]
    

    intg_ΔPQ_id_iq =
        similar( intg_vh_θh_id_iq )
    
    intg_idxs_and_intg_vh_θh_id_iq =
        (;
         intg_vh_θh_id_iq,
         
         intg_gens_vh_Idxs,
         intg_gens_θh_Idxs,
         
         intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs,

         intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs,

         intg_gen_id_Idxs,
         intg_gen_iq_Idxs,
         intg_idxs)

    #---------------------------------------

    
    semi_vh_θh_id_iq =
        [ non_slack_gens_vh;
          non_gen_vh;          
          non_slack_and_non_gen_θh;
          gens_i_d_0;
          gens_i_q_0 ]

    semi_ΔPQ_id_iq =
        similar( semi_vh_θh_id_iq )
    
    semi_idxs_and_semi_vh_θh_id_iq =
        (;
         semi_vh_θh_id_iq,
         
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs,
         
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs,
         
         semi_gen_id_Idxs,
         semi_gen_iq_Idxs,
         semi_idxs)

    #---------------------------------------

    full_vh_θh =
        [gens_vh;
         non_gen_vh;
         [0];
         non_slack_and_non_gen_θh ]

    full_ΔPQ = similar(full_vh_θh)


    full_gens_vh_Idxs =
        first(full_slack_gens_vh_Idxs):last(full_non_slack_gens_vh_Idxs)
    
    full_gens_θh_Idxs =
        first(full_slack_gens_θh_Idxs):last(full_non_slack_gens_θh_Idxs)
    
    full_idxs_and_full_vh_θh =
        (; full_vh_θh,
         full_gens_vh_Idxs,
         full_gens_θh_Idxs,
         
         full_non_gens_nodes_vh_Idxs,
         full_non_gens_nodes_θh_Idxs,
         full_idxs,
         nodes_size)

    #--------------------------------------


end

# ----------------------------------------------------

function get_integrated_nll_dyn_iip_pf_param(
    netd;
    Idxs_type = :Idxs_hybrid )

    """

    Idxs_type = :Idxs_im

    Idxs_type = :Idxs_industrial

    #-----------------------------------------------
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix

    #-----------------------------------------------
    
    netd =
        NetworkData( dynamics_case()... )

    #-----------------------------------------------

    """
    #-----------------------------------------------
    
    loc_load_exist =
        loc_load_exist_bool(netd)

    # ----------------------------------------------

    """
    hybrid model is used to get gens δ_ω_ed_dash_eq_dash
    from  a standard powerflow for the network

    """

    named_tup_pf_result, state =
        get_pf_in_named_tuple_with_state(
            netd)

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    # --------------------------------------


    """
    Idx_converstion_fun:

       `net_to_im_model_indices(idx, dict_conv)`

       `net_to_industrial_model_indices(idx, dict_conv)`

    dict_conv:

    dict_sys_to_im =
        get_net_to_im_indices_dict( netd  )

    dict_sys_to_industry =
        get_net_to_industrial_model_indices_dict(
            netd; no_control_device = false   )

    get_a_model_streamlined_powerflow_net_parameters(
        netd;
        dict_sys_to_model_Idx = dict_sys_to_model_Idx,
        Idx_converstion_fun = Idx_converstion_fun )


        """
    
    models_gens_nodes_some_vars_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )       

    if Idxs_type == :Idxs_hybrid
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs =
            models_gens_nodes_some_vars_Idxs.Idxs_hybrid

        state = state

        δ_ω_ed_dash_eq_dash = δ_ω_ed_dash_eq_dash_view =
            get_gen_nodes_δ_ω_ed_dash_eq_dash(
                state, nodes_δ_ω_ed_dash_eq_dash_Idxs )
        
        pf_net_param =
            get_streamlined_powerflow_net_parameters(
                netd )
                 
        pf_net_dyn_idx_and_Idx =
            pf_net_param.pf_and_dyn_idx_and_Idx
        
        idx_and_Idx =
            pf_net_dyn_idx_and_Idx.hybrid_pf_etc_idx_and_Idx

        
    elseif Idxs_type == :Idxs_industrial
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs =
            models_gens_nodes_some_vars_Idxs.Idxs_industrial

        state = industrial_model_init_operationpoint(
            netd, bus_dict_init
            ;pure = :pure,
            no_control_device = false )
        
        δ_ω_ed_dash_eq_dash = δ_ω_ed_dash_eq_dash_view =
            get_gen_nodes_δ_ω_ed_dash_eq_dash(
                state, nodes_δ_ω_ed_dash_eq_dash_Idxs )

        
        dict_sys_to_model_Idx =
            get_net_to_industrial_model_indices_dict(
                netd;
                no_control_device = false   )
        
        Idx_converstion_fun =
            net_to_industrial_model_indices

        pf_net_param =
            get_a_model_streamlined_powerflow_net_parameters(
                netd;
                dict_sys_to_model_Idx =
                    dict_sys_to_model_Idx,
                Idx_converstion_fun =
                    Idx_converstion_fun )
        
                         
        pf_net_dyn_idx_and_Idx =
            pf_net_param.pf_and_dyn_idx_and_Idx
        
        idx_and_Idx =
            pf_net_dyn_idx_and_Idx.industrial_model_pf_idx_and_Idx


    elseif Idxs_type == :Idxs_im
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs =
            models_gens_nodes_some_vars_Idxs.Idxs_im

        state = im_model_init_operationpoint(
            netd, bus_dict_init )

        δ_ω_ed_dash_eq_dash = δ_ω_ed_dash_eq_dash_view =
            get_gen_nodes_δ_ω_ed_dash_eq_dash(
                state, nodes_δ_ω_ed_dash_eq_dash_Idxs )

        
        dict_sys_to_model_Idx =
            get_net_to_im_indices_dict(
                netd  )
        
        Idx_converstion_fun =
            net_to_im_model_indices


        pf_net_param =
            get_a_model_streamlined_powerflow_net_parameters(
                netd;
                dict_sys_to_model_Idx =
                    dict_sys_to_model_Idx,
                Idx_converstion_fun =
                    Idx_converstion_fun )
                                 
        pf_net_dyn_idx_and_Idx =
            pf_net_param.pf_and_dyn_idx_and_Idx
        
        idx_and_Idx =
            pf_net_dyn_idx_and_Idx.im_model_pf_idx_and_Idx

    else
        nothing
    end
    
    #-----------------------------------------------        
    #-----------------------------------------------
    
    _, pf_idx_and_state, pf_param_views,_, pf_Idx, pf_and_dyn_idx_and_Idx, pf_net_misc = pf_net_param

    # Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, _, _, _ = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, _, _, _, _ = pf_param_views

    # load_trans_nodes_Idx_and_vlimits = pf_limits

    _,nodes_u_Idx,_, _, _, _, _, red_vh_θh_idx,_, _,_ =
        pf_Idx
    
    # ----------------------------------------------------

    # loc_load_exist =
    #     pf_and_dyn_idx_and_Idx.loc_load_exist

    # ----------------------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx


    slack_gens_nodes_idx =
        net_comp_type_idx.slack_gens_nodes_idx

    non_slack_gens_nodes_idx =
        net_comp_type_idx.non_slack_gens_nodes_idx 
    
    gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    gens_with_loc_load_idx =
        net_comp_type_idx.gens_with_loc_load_idx

    load_nodes_idx =
        net_comp_type_idx.load_nodes_idx

    transmission_nodes_idx =
        net_comp_type_idx.transmission_nodes_idx

    load_trans_nodes_idx =
        net_comp_type_idx.load_trans_nodes_idx
        
    # ---------------------------------------------

    gens_with_loc_load_idx =
        loc_load_exist == false ?
        [] : 
        net_comp_type_idx.gens_with_loc_load_idx

    all_nodes_idx =
        net_comp_type_idx.all_nodes_idx

    # ---------------------------------------------
        
    nodes_size = sum(
        [length(gens_nodes_idx),
         length(non_gens_nodes_idx) ])
        
    # ---------------------------------------------

    dict_n2s =
        net_comp_type_idx.dicts_net_to_streamlined_idx

    n2s_slack_gens_idx =
        dict_n2s.dict_n2s_slack_gens_idx

    n2s_non_slack_gens_idx =
        dict_n2s.dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx =
        dict_n2s.dict_n2s_gens_idx
    
    n2s_non_gens_idx =
        dict_n2s.dict_n2s_non_gens_idx

    n2s_load_idx =
        dict_n2s.dict_n2s_load_idx
    
    n2s_gens_with_loc_load_idxs =
        dict_n2s.dict_n2s_gens_with_loc_load_idxs

    n2s_transmission_idxs =
        dict_n2s.dict_n2s_transmission_idxs
    
    n2s_all_nodes_idx =
        dict_n2s.dict_n2s_all_nodes_idx
    
    # --------------------------------------------------

    red_vh_Idxs =
        idx_and_Idx.red_vh_Idxs

    red_θh_Idxs =
        idx_and_Idx.red_θh_Idxs

    non_slack_gens_θh_idx2Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx
    
    non_slack_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx_in_Idx
    
    non_gens_θh_idx2Idx =
        idx_and_Idx.non_gens_θh_idx2Idx
    
    non_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_gens_θh_idx2Idx_in_Idx
            
    # ---------------------------------------------

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx
    
    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash
    
    # ----------------------------------------------

    NL_pf_para_Idxs = vec_Idx.NL_pf_para_Idxs
    
    #----------------------------------------------------

    NL_pf_Idxs =
        NL_pf_para_Idxs.no_loc_load_with_δ_ed_NL_pf_para_Idxs

    P_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_gens_NL_para_Idxs

    Q_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_gens_NL_para_Idxs

    P_load_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_load_NL_para_Idxs

    Q_load_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_load_NL_para_Idxs

    δ_ed_eq_pf_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_δ_ed_eq_pf_NL_para_Idxs        
    
    # ---------------------------------------------------
        
    δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        NL_pf_para_Idxs.δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq
    
    # ---------------------------------------------------

    pf_NL_P_Q_para =
        pf_net_misc.pf_NL_P_Q_para
    
    P_gens =
        pf_NL_P_Q_para.P_gens_NL_para
    
    Q_gens =
        pf_NL_P_Q_para.Q_gens_NL_para
        
    P_load =
        pf_NL_P_Q_para.P_load_NL_para
    
    Q_load =
        pf_NL_P_Q_para.Q_load_NL_para
    
    # ---------------------------------------------------

    no_loc_load_with_δ_ed_eq_pf_param =
        ComponentArray(
            P_gens = P_gens,
            Q_gens = Q_gens,
            P_load = P_load,
            Q_load = Q_load,       
            δ_ed_eq_pf =
                [δ_ω_ed_dash_eq_dash_view...;])
        
    # ---------------------------------------------------

    PQ_gens_pad =
        [ idx ∈ gens_nodes_idx ?
        [ P_gens[n2s_gens_idx[idx]],
          Q_gens[n2s_gens_idx[idx]] ] :
              [0.0, 0.0]
          for idx in 1:nodes_size ]

    PQ_non_gens_pad =
        [ idx ∈ non_gens_nodes_idx ?
        [ P_load[n2s_non_gens_idx[idx]],
          Q_load[n2s_non_gens_idx[idx]] ] :
              [0.0, 0.0]
          for idx in 1:nodes_size ]

    P_gens_pad = first.(PQ_gens_pad)
    
    Q_gens_pad = second.(PQ_gens_pad)

    P_non_gens_pad = first.(PQ_non_gens_pad)
    
    Q_non_gens_pad = second.(PQ_non_gens_pad)

    #-------------------------------------

    dims_no_loc_load_with_δ_ed_eq_padded_PQ =
        [length(P_gens_pad),
         length(Q_gens_pad),
         length(P_non_gens_pad),
         length(Q_non_gens_pad),
         length( [δ_ω_ed_dash_eq_dash_view...;]  ) ]

    _,_, no_loc_load_with_δ_ed_eq_padded_PQ_Idx =
        create_size_offset_Idx(
            dims_no_loc_load_with_δ_ed_eq_padded_PQ ;
            counter = 0)

    P_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[1]

    Q_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[2]

    P_non_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[3]

    Q_non_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[4]

    δ_ω_ed_dash_eq_Idx_in_flat_pad =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[5]

    #----------------------------------------------------

    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist
    
    # ---------------------------------------------------
    
    # x0_vh = [ idx ∈ gens_nodes_idx ?
    #     gens_vh[n2s_gens_idx[idx]] : 1.0
    #           for idx in 1:nodes_size ]
    
    # x0_θh =
    #     [0.0
    #      for idx in
    #          1:nodes_size ]
    
    # uh_0  =
    #     x0_vh .* exp.(im * x0_θh)
    
    # vh_θh_0 =
    #     vcat( x0_vh, x0_θh )
     
    # red_vh_θh_0 =
    #     vh_θh_0[ red_vh_θh_idx  ]


    uh_0  = x_from_xr_xi.(
            [state[idx]
             for idx in
                 nodes_u_Idx])

    vh_θh_0 = [ abs.( uh_0  ); angle.( uh_0  )]
     
    red_vh_θh_0 = vh_θh_0[ red_vh_θh_idx  ]

    
    gens_uh_0 = [ uh_0[  idx  ]
                for idx  in 1:nodes_size
                    if idx ∈ gens_nodes_idx]

    # ---------------------------------------------------
    
    #  named_tuple_with_state = get_pf_in_named_tuple_with_state( netd)
     
    gens_idq_0 =
        [  get_a_gen_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( gens_uh_0  ),
                   angle.( gens_uh_0  ),
                   δ_ω_ed_dash_eq_dash_view,
                   ra_Xd_dash_Xq_dash_view ) ]

    gens_i_d_0 = first.( gens_idq_0 )
        
    gens_i_q_0 = second.( gens_idq_0 )
    
    dims_gens_idq =
        length.( [ gens_i_d_0, gens_i_q_0 ]  )

     _,_,  gens_idq_Idx =
         create_size_offset_Idx(
             dims_gens_idq ;
             counter = 0)

    gens_id_Idx = gens_idq_Idx[1]
        
    gens_iq_Idx = gens_idq_Idx[2]


    gens_idq_0_flat = [ gens_i_d_0; gens_i_q_0 ]

    # ---------------------------------------------------


   """
    t_δ = first.( δ_ω_ed_dash_eq_dash_view )

    t_gens_idq_0 = x_from_xr_xi.( gens_idq_0 )

    t_gens_idq_net_0 = t_gens_idq_0 .*
        exp.(im * (t_δ  .- pi/2) )

    t_gens_S = gens_uh_0 .* conj.(t_gens_idq_net_0 )

     """
    
    # ---------------------------------------------------

    dims_red_vh_θh_0_idq_0 =
        length.([ red_vh_θh_0, gens_idq_0_flat ])

     _,_, red_vh_θh_0_idq_0_Idx =
        create_size_offset_Idx(
            dims_red_vh_θh_0_idq_0;
            counter = 0)

    flat_red_vh_θh_0_Idx =
        red_vh_θh_0_idq_0_Idx[1]

    flat_idq_0_Idx =
        red_vh_θh_0_idq_0_Idx[2]
    
    # ---------------------------------------------------

    red_vh_θh_0_idq = [ red_vh_θh_0; gens_idq_0_flat ]
    
    # ---------------------------------------------------
        
    nodes_types_vh_θh_Idxs =
        pf_and_dyn_idx_and_Idx.hybrid_pf_etc_idx_and_Idx
    
    vh_θh_Idxs_set =
        nodes_types_vh_θh_Idxs.full_nodes_types_Idxs_idx2Idx_etc

    full_slack_gens_vh_Idxs =
        vh_θh_Idxs_set.full_slack_gens_vh_Idxs

    full_non_slack_gens_vh_Idxs =
        vh_θh_Idxs_set.full_non_slack_gens_vh_Idxs

    full_non_gens_nodes_vh_Idxs =
        vh_θh_Idxs_set.full_non_gens_nodes_vh_Idxs

    full_slack_gens_θh_Idxs =
        vh_θh_Idxs_set.full_slack_gens_θh_Idxs

    full_non_slack_gens_θh_Idxs =
        vh_θh_Idxs_set.full_non_slack_gens_θh_Idxs
    
    full_non_gens_nodes_θh_Idxs =
        vh_θh_Idxs_set.full_non_gens_nodes_θh_Idxs

    # ---------------------------------------------------

    dims_full_vh_θh_0_idq_0 =
        length.([ vh_θh_0, gens_idq_0_flat ])

     _,_, full_vh_θh_0_idq_Idx =
        create_size_offset_Idx(
            dims_full_vh_θh_0_idq_0;
            counter = 0)

    full_vh_θh_Idx = full_vh_θh_0_idq_Idx[1]

    full_idq_Idx   = full_vh_θh_0_idq_Idx[2]
    
    # ---------------------------------------------------

    vh_θh_0_idq = [ vh_θh_0; gens_idq_0_flat ]
    
    # ---------------------------------------------------

    red_ΔPQ  = similar( red_vh_θh_0 )
    
    full_ΔPQ = similar( vh_θh_0 )

    red_ΔPQ_Δidq  = similar( red_vh_θh_0_idq )
    
    full_ΔPQ_Δidq = similar( vh_θh_0_idq  )
    
    # -------------------------------------------------- 
    # some selected param
    # --------------------------------------------------
    
    net_idxs =
        (; slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         load_nodes_idx,
         transmission_nodes_idx,
         non_gens_nodes_idx,
         load_trans_nodes_idx,
         gens_with_loc_load_idx)
    
    n2s_idxs =
        (; n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_load_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_transmission_idxs,
         n2s_all_nodes_idx)  

    gens_uh_Q_from_red_sol_para =
        (; slack_vh,
         gens_vh,
         non_slack_gens_Idx_and_vh,
         ra_Xd_dash_Xq_dash_view,         
         red_vh_Idxs,
         red_θh_Idxs,
         red_vh_θh_idx,
         flat_red_vh_θh_0_Idx,
         flat_idq_0_Idx,
         gens_id_Idx,
         gens_iq_Idx,         
         non_gens_θh_idx2Idx,
         non_slack_gens_θh_idx2Idx,         
         net_idxs,
         n2s_idxs  )

    #
    
    # kwd_net_param =
    #     (; pf_net,
    #      pf_idx_and_state,
    #      pf_param_views,
    #      pf_and_dyn_idx_and_Idx,
    #      pf_net_misc,
    #      #pf_Idx,
    #      no_loc_load_with_δ_ed_eq_padded_PQ_Idx,
    #      gens_idq_Idx,
    #      red_vh_θh_0_idq_0_Idx,
    #      full_vh_θh_0_idq_Idx,
    #      δ_ω_ed_dash_eq_dash_view )

    
    kwd_net_param =
        (; pf_net_param,
         gens_idq_Idx,
         red_vh_θh_0_idq_0_Idx,
         full_vh_θh_0_idq_Idx )
    
    nll_dyn_iip_pf_param =
        (; 
         red_ΔPQ_Δidq,
         red_vh_θh_0_idq,
         gens_uh_Q_from_red_sol_para,
         kwd_net_param
         )
    
    # return (; nll_dyn_iip_pf_param, pois, pois_Idxs)
    flat_δ_ω_ed_dash_eq_dash =
        [δ_ω_ed_dash_eq_dash_view...;]
    
    return (; nll_dyn_iip_pf_param,
            flat_δ_ω_ed_dash_eq_dash ) 
end



function driver_get_integrated_nll_dyn_iip_pf_param()

    abstol=1e-12
    
    reltol=1e-12
    
    init_pf = true
    
    pf_alg = NewtonRaphson()
    
    with_δ_ed_eq = false

    # :Idxs_industrial, :Idxs_im
    
    Idxs_type = :Idxs_hybrid  
    
    #-----------------------------------------------
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix

    #-----------------------------------------------
    
    netd =
        NetworkData( dynamics_case()... )

    #-----------------------------------------------

    (; nll_dyn_iip_pf_param,
     flat_δ_ω_ed_dash_eq_dash ) =
         get_integrated_nll_dyn_iip_pf_param(
             netd;
             Idxs_type = Idxs_type )
    
    (; 
     red_ΔPQ_Δidq,
     red_vh_θh_0_idq,
     gens_uh_Q_from_red_sol_para,
     kwd_net_param
     ) =
         nll_dyn_iip_pf_param
        
    
    (; pf_net_param,
     gens_idq_Idx,
     red_vh_θh_0_idq_0_Idx,
     full_vh_θh_0_idq_Idx ) =
         kwd_net_param 
        
    
    (; slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     gens_with_loc_load_idx) =
         net_idxs
        
    
    (; n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx) =
         n2s_idxs 
          

    (; slack_vh,
     gens_vh,
     non_slack_gens_Idx_and_vh,
     ra_Xd_dash_Xq_dash_view,         
     red_vh_Idxs,
     red_θh_Idxs,     
     red_vh_θh_idx,
     flat_red_vh_θh_0_Idx,
     flat_idq_0_Idx,
     gens_id_Idx,
     gens_iq_Idx,         
     non_gens_θh_idx2Idx,
     non_slack_gens_θh_idx2Idx,         
     net_idxs,
     n2s_idxs  ) =
         gens_uh_Q_from_red_sol_para
        
    return nothing
end


function get_integrated_hybrid_pf_etc_idx_and_Idx(netd)

    #------------------------------------------
    slack_gens_nodes_idx =
         get_slack_gens_nodes_idx(
             netd.nodes )

     non_slack_gens_nodes_idx =
         get_non_slack_gens_nodes_idx(
             netd.nodes )

     gens_nodes_idx =
         get_gens_nodes_idx( netd.nodes )

     gens_nodes_with_loc_loads_idx =
         get_gens_nodes_with_loc_loads_idx(
             netd.nodes )
    
    loc_load_exist =
        gens_nodes_with_loc_loads_idx == [] ?
        false : true
    
     load_nodes_idx =
         get_load_nodes_idx( netd.nodes )

     transmission_nodes_idx =
         get_transmission_nodes_idx( netd.nodes )

     non_gens_nodes_idx =
         get_non_gens_nodes_idx( netd.nodes )

     load_trans_nodes_idx =
         get_load_trans_nodes_Idx( netd.nodes )

     all_nodes_idx =
         get_all_nodes_idx( netd.nodes  )

     non_slack_gens_and_non_gens_idx =
         setdiff(all_nodes_idx,
                 slack_gens_nodes_idx)
        
    # -------------------------------------------------
    
    slack_ur_ui_Idx_in_state =
        get_components_slack_ur_ui_Idx_in_state(
            netd.nodes )

    non_slack_ur_ui_Idx_in_state =
        get_components_no_slack_ur_ui_Idx_in_state(
            netd.nodes )
    
    ur_ui_Idx_in_state =
        get_components_ur_ui_Idx_in_state(
            netd.nodes )
    
    nodes_u_Idx =
        get_components_ur_ui_Idx_in_state(
            netd.nodes )
    
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    dims_ur_ui   = [ length(ur_idx), length(ui_idx) ]
    
    _,_, ur_ui_IDX =
        create_size_offset_Idx(
            dims_ur_ui ;
            counter = 0)

    ur_IDX,  ui_IDX = ur_ui_IDX

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX
    
    # ----------------------------------------------------

    slack_bus_idx =
        get_components_slack_Idx( netd.nodes )
    
    gens_idx = first.(
        get_generators_Idx_and_vh(
            netd.nodes ))        

     # ----------------------------------------------------

     red_vh_θh_idx = [ setdiff(vh_IDX, gens_idx)...;
                setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ---------------------------------------------------   
    
     red_vh_idx = setdiff(vh_IDX, gens_idx)

     red_θh_idx =
         setdiff( θh_IDX, θh_IDX[slack_bus_idx] )

    # ----------------------------------------------------

    red_vh_θh_dims =
        [ length( red_vh_idx ), length( red_θh_idx ) ]    

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(red_vh_θh_dims )

    red_vh_Idxs = red_vh_θh_IDX[1]
    
    red_θh_Idxs = red_vh_θh_IDX[2]
    
    # -------------------------------------------------    

    full_vh_θh_dims =
        [ length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx ),
          length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx )
          ]
    
     _, _, full_vh_θh_IDX =
        create_size_offset_Idx( full_vh_θh_dims )

    full_slack_gens_vh_Idxs =
        full_vh_θh_IDX[1]
    
    full_non_slack_gens_vh_Idxs =
        full_vh_θh_IDX[2]

    full_non_gens_nodes_vh_Idxs =
        full_vh_θh_IDX[3]
    
    full_slack_gens_θh_Idxs =
        full_vh_θh_IDX[4]

    full_non_slack_gens_θh_Idxs =
        full_vh_θh_IDX[5]
    
    full_non_gens_nodes_θh_Idxs =
        full_vh_θh_IDX[6]

    full_nodes_types_vh_Idxs =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs)

    full_nodes_types_θh_Idxs =
        (;full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs)

    full_nodes_types_vh_and_θh_Idxs =
        (; full_nodes_types_vh_Idxs,
         full_nodes_types_θh_Idxs)
    
     # --------------------------------------------------  
    
    full_nodes_type_idxs =
        [slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    full_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...],
                      vh_IDX) )
    
    full_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...],
                      θh_IDX ) )
    
    vec_types_full_vh_idx2Idx = [
        [ full_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx, full_non_slack_gens_vh_idx2Idx, full_non_gens_vh_idx2Idx = vec_types_full_vh_idx2Idx
        

    vec_types_full_θh_idx2Idx = [
        [ full_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx, full_non_slack_gens_θh_idx2Idx, full_non_gens_θh_idx2Idx = vec_types_full_θh_idx2Idx


    full_nodes_types_vh_idx2Idx =
        (; full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx)
    
    full_nodes_types_θh_idx2Idx =
        (; full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx)

    full_nodes_types_vh_and_θh_idx2Idx =
        (; full_nodes_types_vh_idx2Idx,
         full_nodes_types_θh_idx2Idx)

    full_nodes_types_dict_vh_and_θh_idx2Idx =
        (; full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx )

    
     # -----------------------------------------------   

    
    full_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...],
                      1:length(vh_IDX) ) )

    
    full_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...],
                      1:length( θh_IDX ) ) )

    
    vec_types_full_vh_idx2Idx_in_Idx = [
        [ full_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx_in_Idx, full_non_slack_gens_vh_idx2Idx_in_Idx, full_non_gens_vh_idx2Idx_in_Idx = vec_types_full_vh_idx2Idx_in_Idx


    vec_types_full_θh_idx2Idx_in_Idx = [
        [ full_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx_in_Idx, full_non_slack_gens_θh_idx2Idx_in_Idx, full_non_gens_θh_idx2Idx_in_Idx = vec_types_full_θh_idx2Idx_in_Idx


    full_nodes_types_vh_idx2Idx_in_Idx =
        (; full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx   )
    
    full_nodes_types_θh_idx2Idx_in_Idx =
        (; full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx )

    full_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; full_nodes_types_vh_idx2Idx_in_Idx,
         full_nodes_types_θh_idx2Idx_in_Idx )
    
    full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )

     # ------------------------------------------------   

    full_nodes_types_Idxs_idx2Idx_etc =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs,
         full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx,    
         full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx,
         full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx,
         full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx,
         full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx,    
         full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx
         )
    
    
    # --------------------------------------------------   

     dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      1:length(red_θh_Idxs) ) )

    vec_red_nodes_types_idxs =
        [ non_slack_gens_nodes_idx,
          non_gens_nodes_idx ]

    vec_types_red_idx2Idx = [
        [ dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            vec_red_nodes_types_idxs  ]

    non_slack_gens_θh_idx2Idx, non_gens_θh_idx2Idx =
        vec_types_red_idx2Idx

    # -----------------------------------


    vec_types_red_idx2Idx_in_Idx = [
        [ dict_θh_idx2Idx_in_Idx[ idx ]
          for idx in types_idxs ]
        for types_idxs in
            vec_red_nodes_types_idxs  ]

    non_slack_gens_θh_idx2Idx_in_Idx, non_gens_θh_idx2Idx_in_Idx = vec_types_red_idx2Idx_in_Idx

    # ----------------------------------------------------
    # red_vh_idx, red_θh_idx, red_vh_Idxs, red_θh_Idxs
    # non_slack_gens_θh_idx2Idx,
    # non_slack_gens_θh_idx2Idx_in_Idx,
    # non_gens_θh_idx2Idx, non_gens_θh_idx2Idx_in_Idx
    # ----------------------------------------------------

    # ----------------------------------------------------
    # Additional Idx
    # ----------------------------------------------------

    dim_intg_vh_θh_id_iq =
        [
            length( slack_gens_nodes_idx ),            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( slack_gens_nodes_idx ),               
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, intg_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_intg_vh_θh_id_iq  )

    intg_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[1]
    
    intg_non_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[2]

    intg_non_gens_nodes_vh_Idxs =
        intg_vh_θh_id_iq_IDX[3]
    
    intg_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[4]

    intg_non_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[5]
    
    intg_non_gens_nodes_θh_Idxs =
        intg_vh_θh_id_iq_IDX[6]

    intg_gen_id_Idxs =
        intg_vh_θh_id_iq_IDX[7]

    intg_gen_iq_Idxs =
        intg_vh_θh_id_iq_IDX[8]

    intg_nodes_types_vh_Idxs =
        (; intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs)

    intg_nodes_types_θh_Idxs =
        (;intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs)

    intg_nodes_idq_Idxs =
        (; intg_gen_id_Idxs,
           intg_gen_iq_Idxs )

    intg_nodes_types_vh_θh_id_iq_Idxs =
        (; intg_nodes_types_vh_Idxs,
         intg_nodes_types_θh_Idxs,
         intg_nodes_idq_Idxs)
    
     # --------------------------------------------------  
    
    intg_nodes_type_idxs =
        [slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    intg_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      vh_IDX) )
    
    intg_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_intg_vh_idx2Idx = [
        [ intg_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx, intg_non_slack_gens_vh_idx2Idx, intg_non_gens_vh_idx2Idx = vec_types_intg_vh_idx2Idx
        

    vec_types_intg_θh_idx2Idx = [
        [ intg_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx, intg_non_slack_gens_θh_idx2Idx, intg_non_gens_θh_idx2Idx = vec_types_intg_θh_idx2Idx


    intg_nodes_types_vh_idx2Idx =
        (; intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx)
    
    intg_nodes_types_θh_idx2Idx =
        (; intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx)

    intg_nodes_types_vh_and_θh_idx2Idx =
        (; intg_nodes_types_vh_idx2Idx,
         intg_nodes_types_θh_idx2Idx)

    intg_nodes_types_dict_vh_and_θh_idx2Idx =
        (; intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx )

    
     # --------------------------------------------------   
    
    intg_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    intg_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_intg_vh_idx2Idx_in_Idx = [
        [ intg_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx_in_Idx, intg_non_slack_gens_vh_idx2Idx_in_Idx, intg_non_gens_vh_idx2Idx_in_Idx = vec_types_intg_vh_idx2Idx_in_Idx


    vec_types_intg_θh_idx2Idx_in_Idx = [
        [ intg_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx_in_Idx, intg_non_slack_gens_θh_idx2Idx_in_Idx, intg_non_gens_θh_idx2Idx_in_Idx = vec_types_intg_θh_idx2Idx_in_Idx


    intg_nodes_types_vh_idx2Idx_in_Idx =
        (; intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx   )
    
    intg_nodes_types_θh_idx2Idx_in_Idx =
        (; intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx )

    intg_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; intg_nodes_types_vh_idx2Idx_in_Idx,
         intg_nodes_types_θh_idx2Idx_in_Idx )
    
    intg_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx )

     # ------------------------------------------------
    
    intg_nodes_types_Idxs_idx2Idx_etc =
        (;
         intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs,
         
         intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs,
         
         intg_gen_id_Idxs,
         intg_gen_iq_Idxs,
         
         intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx,
         
         intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx,
         
         intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx,
         
         intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx,
         
         intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx,
         
         intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx
         )

     # ----------------------------------------------

    dim_semi_vh_θh_id_iq =
        [
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, semi_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_semi_vh_θh_id_iq  )
    
    semi_non_slack_gens_vh_Idxs =
        semi_vh_θh_id_iq_IDX[1]

    semi_non_gens_nodes_vh_Idxs =
        semi_vh_θh_id_iq_IDX[2]
    
    semi_non_slack_gens_θh_Idxs =
        semi_vh_θh_id_iq_IDX[3]
    
    semi_non_gens_nodes_θh_Idxs =
        semi_vh_θh_id_iq_IDX[4]

    semi_gen_id_Idxs =
        semi_vh_θh_id_iq_IDX[5]

    semi_gen_iq_Idxs =
        semi_vh_θh_id_iq_IDX[6]

    semi_nodes_types_vh_Idxs =
        (;
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs)

    semi_nodes_types_θh_Idxs =
        (;
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs)

    semi_nodes_idq_Idxs =
        (; semi_gen_id_Idxs,
           semi_gen_iq_Idxs )

    semi_nodes_types_vh_θh_id_iq_Idxs =
        (; semi_nodes_types_vh_Idxs,
         semi_nodes_types_θh_Idxs,
         semi_nodes_idq_Idxs)
    
     # -----------------------------------------------  
    
    semi_nodes_type_idxs =
        [
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    semi_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      vh_IDX) )
    
    semi_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_semi_vh_idx2Idx = [
        [ semi_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx, semi_non_gens_vh_idx2Idx = vec_types_semi_vh_idx2Idx
        

    vec_types_semi_θh_idx2Idx = [
        [ semi_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx, semi_non_gens_θh_idx2Idx = vec_types_semi_θh_idx2Idx


    semi_nodes_types_vh_idx2Idx =
        (; 
         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx)
    
    semi_nodes_types_θh_idx2Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx)

    semi_nodes_types_vh_and_θh_idx2Idx =
        (; semi_nodes_types_vh_idx2Idx,
         semi_nodes_types_θh_idx2Idx)

    semi_nodes_types_dict_vh_and_θh_idx2Idx =
        (; semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx )

    
     # ----------------------------------------------
    
    semi_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    semi_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_semi_vh_idx2Idx_in_Idx = [
        [ semi_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx_in_Idx, semi_non_gens_vh_idx2Idx_in_Idx = vec_types_semi_vh_idx2Idx_in_Idx


    vec_types_semi_θh_idx2Idx_in_Idx = [
        [ semi_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx_in_Idx, semi_non_gens_θh_idx2Idx_in_Idx = vec_types_semi_θh_idx2Idx_in_Idx


    semi_nodes_types_vh_idx2Idx_in_Idx =
        (;
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx   )
    
    semi_nodes_types_θh_idx2Idx_in_Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx )

    semi_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; semi_nodes_types_vh_idx2Idx_in_Idx,
         semi_nodes_types_θh_idx2Idx_in_Idx )
    
    semi_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx )

     # -----------------------------------------------

    semi_nodes_types_Idxs_idx2Idx_etc =
        (;
         
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs,
                  
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs,
         
         semi_gen_id_Idxs,
         semi_gen_iq_Idxs,

         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx,

         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx,
         
         semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx,
         
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx,
         
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx,
         
         semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx
         )
        
    # -------------------------------------------------- 

    return (; slack_ur_ui_Idx_in_state,
            non_slack_ur_ui_Idx_in_state,
            ur_ui_Idx_in_state,
            nodes_u_Idx,
            ur_idx,
            ui_idx,
            ur_ui_idx,
            vh_θh_idx,
            ur_ui_IDX,
            ur_IDX,
            ui_IDX,
            ir_IDX,
            ii_IDX,
            vh_IDX,
            θh_IDX,
            slack_bus_idx,
            gens_idx,
            red_vh_θh_idx,
            red_vh_idx,
            red_θh_idx,
            red_vh_Idxs,
            red_θh_Idxs,
            non_slack_gens_θh_idx2Idx,
            non_slack_gens_θh_idx2Idx_in_Idx,
            non_gens_θh_idx2Idx,
            non_gens_θh_idx2Idx_in_Idx,
            loc_load_exist,
            full_nodes_types_Idxs_idx2Idx_etc,
            intg_nodes_types_Idxs_idx2Idx_etc,
            semi_nodes_types_Idxs_idx2Idx_etc)
    
    
end


#------------------------------------------------


function get_integrated_im_model_pf_idx_and_Idx(netd)

    #------------------------------------------        

    slack_bus_idx =
        get_components_slack_Idx( netd.nodes )
    
    gens_idx = first.(
        get_generators_Idx_and_vh(
            netd.nodes ))        

    #------------------------------------------        
    
     slack_gens_nodes_idx =
         get_slack_gens_nodes_idx(
             netd.nodes )

     non_slack_gens_nodes_idx =
         get_non_slack_gens_nodes_idx(
             netd.nodes )

     gens_nodes_idx =
         get_gens_nodes_idx( netd.nodes )

     gens_nodes_with_loc_loads_idx =
         get_gens_nodes_with_loc_loads_idx(
             netd.nodes )
    
    loc_load_exist =
        gens_nodes_with_loc_loads_idx == [] ?
        false : true

     load_nodes_idx =
         get_load_nodes_idx( netd.nodes )

     transmission_nodes_idx =
         get_transmission_nodes_idx( netd.nodes )

     non_gens_nodes_idx =
         get_non_gens_nodes_idx( netd.nodes )

     load_trans_nodes_idx =
         get_load_trans_nodes_Idx( netd.nodes )

     all_nodes_idx =
         get_all_nodes_idx( netd.nodes  )

     non_slack_gens_and_non_gens_idx =
         setdiff(all_nodes_idx, slack_gens_nodes_idx)

    #------------------------------------------        
    
    dict_sys_to_industry =
        get_net_to_im_indices_dict( netd  )

    #------------------------------------------        

    slack_ur_ui_Idx_in_state =
        net_to_im_model_indices(
            get_components_slack_ur_ui_Idx_in_state(
                netd.nodes ), dict_sys_to_industry )

    non_slack_ur_ui_Idx_in_state =
        net_to_im_model_indices(
            get_components_no_slack_ur_ui_Idx_in_state(
                netd.nodes),
            dict_sys_to_industry  )

    ur_ui_Idx_in_state =
        net_to_im_model_indices(
            get_components_ur_ui_Idx_in_state(
                netd.nodes),
            dict_sys_to_industry  )

    nodes_u_Idx =
        net_to_im_model_indices(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )
    
    #------------------------------------------        
        
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   =
        [ length(ur_idx), length(ui_idx) ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    # ----------------------------------------------------

    red_vh_θh_idx =
        [ setdiff(vh_IDX, gens_idx)...;
          setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]

    # ----------------------------------------------------

    red_vh_idx = setdiff(vh_IDX, gens_idx)
    
    red_θh_idx =
        setdiff(θh_IDX, θh_IDX[slack_bus_idx])

    # ----------------------------------------------------

    red_vh_θh_dims =
        [ length( red_vh_idx ), length( red_θh_idx ) ]    

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(red_vh_θh_dims )

    red_vh_Idxs = red_vh_θh_IDX[1]
    
    red_θh_Idxs = red_vh_θh_IDX[2]

     # ------------------------------------------------   

     dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      1:length(red_θh_Idxs) ) )

     non_slack_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_slack_gens_nodes_idx ]

     non_slack_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_slack_gens_nodes_idx ]


     non_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_gens_nodes_idx ]

     non_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_gens_nodes_idx ]
    
     # --------------------------------------------------  

    full_vh_θh_dims =
        [ length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx ),
          length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx )
          ]
    
     _, _, full_vh_θh_IDX =
        create_size_offset_Idx( full_vh_θh_dims )

    full_slack_gens_vh_Idxs =
        full_vh_θh_IDX[1]
    
    full_non_slack_gens_vh_Idxs =
        full_vh_θh_IDX[2]

    full_non_gens_nodes_vh_Idxs =
        full_vh_θh_IDX[3]
    
    full_slack_gens_θh_Idxs =
        full_vh_θh_IDX[4]

    full_non_slack_gens_θh_Idxs =
        full_vh_θh_IDX[5]
    
    full_non_gens_nodes_θh_Idxs =
        full_vh_θh_IDX[6]

    full_nodes_types_vh_Idxs =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs)

    full_nodes_types_θh_Idxs =
        (;full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs)

    full_nodes_types_vh_and_θh_Idxs =
        (; full_nodes_types_vh_Idxs,
         full_nodes_types_θh_Idxs)
    
     # --------------------------------------------------  
    
    full_nodes_type_idxs =
        [slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    full_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      vh_IDX) )
    
    full_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_full_vh_idx2Idx = [
        [ full_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx, full_non_slack_gens_vh_idx2Idx, full_non_gens_vh_idx2Idx = vec_types_full_vh_idx2Idx
        

    vec_types_full_θh_idx2Idx = [
        [ full_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx, full_non_slack_gens_θh_idx2Idx, full_non_gens_θh_idx2Idx = vec_types_full_θh_idx2Idx


    full_nodes_types_vh_idx2Idx =
        (; full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx)
    
    full_nodes_types_θh_idx2Idx =
        (; full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx)

    full_nodes_types_vh_and_θh_idx2Idx =
        (; full_nodes_types_vh_idx2Idx,
         full_nodes_types_θh_idx2Idx)

    full_nodes_types_dict_vh_and_θh_idx2Idx =
        (; full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx )
    
     # -----------------------------------------------   
    
    full_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    full_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_full_vh_idx2Idx_in_Idx = [
        [ full_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx_in_Idx, full_non_slack_gens_vh_idx2Idx_in_Idx, full_non_gens_vh_idx2Idx_in_Idx = vec_types_full_vh_idx2Idx_in_Idx


    vec_types_full_θh_idx2Idx_in_Idx = [
        [ full_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx_in_Idx, full_non_slack_gens_θh_idx2Idx_in_Idx, full_non_gens_θh_idx2Idx_in_Idx = vec_types_full_θh_idx2Idx_in_Idx


    full_nodes_types_vh_idx2Idx_in_Idx =
        (; full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx   )
    
    full_nodes_types_θh_idx2Idx_in_Idx =
        (; full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx )

    full_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; full_nodes_types_vh_idx2Idx_in_Idx,
         full_nodes_types_θh_idx2Idx_in_Idx )
    
    full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )

     # ----------------------------------------------   

    full_nodes_types_Idxs_idx2Idx_etc =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs,
         full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx,    
         full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx,
         full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx,
         full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx,
         full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx,    
         full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx
         )
        

    # ----------------------------------------------------

    # red_vh_idx, red_θh_idx, red_vh_Idxs, red_θh_Idxs
    # non_slack_gens_θh_idx2Idx,
    # non_slack_gens_θh_idx2Idx_in_Idx,
    # non_gens_θh_idx2Idx, non_gens_θh_idx2Idx_in_Idx
    # ----------------------------------------------------

    # ----------------------------------------------------
    # Additional Idx
    # ----------------------------------------------------

    dim_intg_vh_θh_id_iq =
        [
            length( slack_gens_nodes_idx ),            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( slack_gens_nodes_idx ),               
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, intg_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_intg_vh_θh_id_iq  )

    intg_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[1]
    
    intg_non_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[2]

    intg_non_gens_nodes_vh_Idxs =
        intg_vh_θh_id_iq_IDX[3]
    
    intg_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[4]

    intg_non_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[5]
    
    intg_non_gens_nodes_θh_Idxs =
        intg_vh_θh_id_iq_IDX[6]

    intg_gen_id_Idxs =
        intg_vh_θh_id_iq_IDX[7]

    intg_gen_iq_Idxs =
        intg_vh_θh_id_iq_IDX[8]

    intg_nodes_types_vh_Idxs =
        (; intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs)

    intg_nodes_types_θh_Idxs =
        (;intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs)

    intg_nodes_idq_Idxs =
        (; intg_gen_id_Idxs,
           intg_gen_iq_Idxs )

    intg_nodes_types_vh_θh_id_iq_Idxs =
        (; intg_nodes_types_vh_Idxs,
         intg_nodes_types_θh_Idxs,
         intg_nodes_idq_Idxs)
    
     # --------------------------------------------------  
    
    intg_nodes_type_idxs =
        [slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    intg_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      vh_IDX) )
    
    intg_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_intg_vh_idx2Idx = [
        [ intg_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx, intg_non_slack_gens_vh_idx2Idx, intg_non_gens_vh_idx2Idx = vec_types_intg_vh_idx2Idx
        

    vec_types_intg_θh_idx2Idx = [
        [ intg_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx, intg_non_slack_gens_θh_idx2Idx, intg_non_gens_θh_idx2Idx = vec_types_intg_θh_idx2Idx


    intg_nodes_types_vh_idx2Idx =
        (; intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx)
    
    intg_nodes_types_θh_idx2Idx =
        (; intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx)

    intg_nodes_types_vh_and_θh_idx2Idx =
        (; intg_nodes_types_vh_idx2Idx,
         intg_nodes_types_θh_idx2Idx)

    intg_nodes_types_dict_vh_and_θh_idx2Idx =
        (; intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx )

    
     # ----------------------------------------------   
    
    intg_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    intg_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_intg_vh_idx2Idx_in_Idx = [
        [ intg_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx_in_Idx, intg_non_slack_gens_vh_idx2Idx_in_Idx, intg_non_gens_vh_idx2Idx_in_Idx = vec_types_intg_vh_idx2Idx_in_Idx


    vec_types_intg_θh_idx2Idx_in_Idx = [
        [ intg_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx_in_Idx, intg_non_slack_gens_θh_idx2Idx_in_Idx, intg_non_gens_θh_idx2Idx_in_Idx = vec_types_intg_θh_idx2Idx_in_Idx


    intg_nodes_types_vh_idx2Idx_in_Idx =
        (; intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx   )
    
    intg_nodes_types_θh_idx2Idx_in_Idx =
        (; intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx )

    intg_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; intg_nodes_types_vh_idx2Idx_in_Idx,
         intg_nodes_types_θh_idx2Idx_in_Idx )
    
    intg_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx )

     # -----------------------------------------------   

    intg_nodes_types_Idxs_idx2Idx_etc =
        (;
         intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs,
         
         intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs,
         
         intg_gen_id_Idxs,
         intg_gen_iq_Idxs,
         
         intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx,
         
         intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx,
         
         intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx,
         
         intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx,
         
         intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx,
         
         intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx
         )

     # ----------------------------------------------   
     # ----------------------------------------------   

    dim_semi_vh_θh_id_iq =
        [
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, semi_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_semi_vh_θh_id_iq  )
    
    semi_non_slack_gens_vh_Idxs =
        semi_vh_θh_id_iq_IDX[1]

    semi_non_gens_nodes_vh_Idxs =
        semi_vh_θh_id_iq_IDX[2]
    
    semi_non_slack_gens_θh_Idxs =
        semi_vh_θh_id_iq_IDX[3]
    
    semi_non_gens_nodes_θh_Idxs =
        semi_vh_θh_id_iq_IDX[4]

    semi_gen_id_Idxs =
        semi_vh_θh_id_iq_IDX[5]

    semi_gen_iq_Idxs =
        semi_vh_θh_id_iq_IDX[6]

    semi_nodes_types_vh_Idxs =
        (;
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs)

    semi_nodes_types_θh_Idxs =
        (;
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs)

    semi_nodes_idq_Idxs =
        (; semi_gen_id_Idxs,
           semi_gen_iq_Idxs )

    semi_nodes_types_vh_θh_id_iq_Idxs =
        (; semi_nodes_types_vh_Idxs,
         semi_nodes_types_θh_Idxs,
         semi_nodes_idq_Idxs)
    
     # -------------------------------------------------  
    
    semi_nodes_type_idxs =
        [
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    semi_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      vh_IDX) )
    
    semi_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_semi_vh_idx2Idx = [
        [ semi_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx, semi_non_gens_vh_idx2Idx = vec_types_semi_vh_idx2Idx
        

    vec_types_semi_θh_idx2Idx = [
        [ semi_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx, semi_non_gens_θh_idx2Idx = vec_types_semi_θh_idx2Idx


    semi_nodes_types_vh_idx2Idx =
        (; 
         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx)
    
    semi_nodes_types_θh_idx2Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx)

    semi_nodes_types_vh_and_θh_idx2Idx =
        (; semi_nodes_types_vh_idx2Idx,
         semi_nodes_types_θh_idx2Idx)

    semi_nodes_types_dict_vh_and_θh_idx2Idx =
        (; semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx )

    
     # ------------------------------------------------   
    
    semi_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    semi_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_semi_vh_idx2Idx_in_Idx = [
        [ semi_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx_in_Idx, semi_non_gens_vh_idx2Idx_in_Idx = vec_types_semi_vh_idx2Idx_in_Idx


    vec_types_semi_θh_idx2Idx_in_Idx = [
        [ semi_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx_in_Idx, semi_non_gens_θh_idx2Idx_in_Idx = vec_types_semi_θh_idx2Idx_in_Idx


    semi_nodes_types_vh_idx2Idx_in_Idx =
        (;
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx   )
    
    semi_nodes_types_θh_idx2Idx_in_Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx )

    semi_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; semi_nodes_types_vh_idx2Idx_in_Idx,
         semi_nodes_types_θh_idx2Idx_in_Idx )
    
    semi_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx )

     # ------------------------------------------------   

    semi_nodes_types_Idxs_idx2Idx_etc =
        (;
         
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs,
                  
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs,
         
         semi_gen_id_Idxs,
         semi_gen_iq_Idxs,

         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx,

         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx,
         
         semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx,
         
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx,
         
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx,
         
         semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx
         )
        
    # -----------------------------------------------  

    return (; slack_ur_ui_Idx_in_state,
            non_slack_ur_ui_Idx_in_state,
            ur_ui_Idx_in_state,
            nodes_u_Idx,
            ur_idx,
            ui_idx,
            ur_ui_idx,
            vh_θh_idx,
            ur_ui_IDX,
            ur_IDX,
            ui_IDX,
            ir_IDX,
            ii_IDX,
            vh_IDX,
            θh_IDX,
            slack_bus_idx,
            gens_idx,
            red_vh_θh_idx,
            red_vh_idx,
            red_θh_idx,
            red_vh_Idxs,
            red_θh_Idxs,
            non_slack_gens_θh_idx2Idx,
            non_slack_gens_θh_idx2Idx_in_Idx,
            non_gens_θh_idx2Idx,
            non_gens_θh_idx2Idx_in_Idx,
            loc_load_exist,
            full_nodes_types_Idxs_idx2Idx_etc,
            intg_nodes_types_Idxs_idx2Idx_etc,
            semi_nodes_types_Idxs_idx2Idx_etc)


end



function get_integrated_industrial_model_pf_idx_and_Idx(
    netd; no_control_device = false )

    #------------------------------------------        

    slack_bus_idx =
        get_components_slack_Idx( netd.nodes )
    
    gens_idx = first.(
        get_generators_Idx_and_vh(
            netd.nodes ))        

    # ----------------------------------------------------
    
     slack_gens_nodes_idx =
         get_slack_gens_nodes_idx(
             netd.nodes )

     non_slack_gens_nodes_idx =
         get_non_slack_gens_nodes_idx(
             netd.nodes )

     gens_nodes_idx =
         get_gens_nodes_idx( netd.nodes )

     gens_nodes_with_loc_loads_idx =
         get_gens_nodes_with_loc_loads_idx(
             netd.nodes )
    
    loc_load_exist =
        gens_nodes_with_loc_loads_idx == [] ?
        false : true

     load_nodes_idx =
         get_load_nodes_idx( netd.nodes )

     transmission_nodes_idx =
         get_transmission_nodes_idx( netd.nodes )

     non_gens_nodes_idx =
         get_non_gens_nodes_idx( netd.nodes )

     load_trans_nodes_idx =
         get_load_trans_nodes_Idx( netd.nodes )

     all_nodes_idx =
         get_all_nodes_idx( netd.nodes  )

     non_slack_gens_and_non_gens_idx =
         setdiff(all_nodes_idx, slack_gens_nodes_idx)

     # ------------------------------------------------

    if no_control_device == false

        dict_sys_to_industry =
            get_net_to_industrial_model_indices_dict( netd  )
        
    else
        dict_sys_to_industry =
            get_net_to_industrial_model_indices_dict(
                netd;  no_control_device = true  )
    end
    

    #------------------------------------------        

    slack_ur_ui_Idx_in_state =
        net_to_industrial_model_indices(
            get_components_slack_ur_ui_Idx_in_state(
                netd.nodes ), dict_sys_to_industry  )

    non_slack_ur_ui_Idx_in_state =
        net_to_industrial_model_indices(
            get_components_no_slack_ur_ui_Idx_in_state(netd.nodes),
            dict_sys_to_industry  )

    ur_ui_Idx_in_state =
        net_to_industrial_model_indices(
            get_components_ur_ui_Idx_in_state( netd.nodes ),
            dict_sys_to_industry  )

    nodes_u_Idx =
        net_to_industrial_model_indices(
            get_components_ur_ui_Idx_in_state(
                netd.nodes ),
            dict_sys_to_industry  )
    
    #------------------------------------------        
        
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   =
        [ length(ur_idx), length(ui_idx) ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    # ----------------------------------------------------
    
    red_vh_θh_idx =
        [ setdiff(vh_IDX, gens_idx)...;
          setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ----------------------------------------------------

    red_vh_idx = setdiff(vh_IDX, gens_idx)
    
    red_θh_idx =
        setdiff(θh_IDX, θh_IDX[slack_bus_idx])

    # ----------------------------------------------------

    red_vh_θh_dims =
        [ length( red_vh_idx ), length( red_θh_idx ) ]    

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(red_vh_θh_dims )

    red_vh_Idxs = red_vh_θh_IDX[1]
    
    red_θh_Idxs = red_vh_θh_IDX[2]
    
     # -------------------------------------------------- 
    
     dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      1:length(red_θh_Idxs) ) )

     non_slack_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_slack_gens_nodes_idx ]

     non_slack_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_slack_gens_nodes_idx ]


     non_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_gens_nodes_idx ]

     non_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_gens_nodes_idx ]
    
    # ----------------------------------------------------

    full_vh_θh_dims =
        [ length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx ),
          length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx )
          ]
    
     _, _, full_vh_θh_IDX =
        create_size_offset_Idx( full_vh_θh_dims )

    full_slack_gens_vh_Idxs =
        full_vh_θh_IDX[1]
    
    full_non_slack_gens_vh_Idxs =
        full_vh_θh_IDX[2]

    full_non_gens_nodes_vh_Idxs =
        full_vh_θh_IDX[3]
    
    full_slack_gens_θh_Idxs =
        full_vh_θh_IDX[4]

    full_non_slack_gens_θh_Idxs =
        full_vh_θh_IDX[5]
    
    full_non_gens_nodes_θh_Idxs =
        full_vh_θh_IDX[6]

    full_nodes_types_vh_Idxs =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs )

    full_nodes_types_θh_Idxs =
        (;full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs)

    full_nodes_types_vh_and_θh_Idxs =
        (; full_nodes_types_vh_Idxs,
         full_nodes_types_θh_Idxs)
    
     # --------------------------------------------------  
    
    full_nodes_type_idxs =
        [ slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    full_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      vh_IDX) )
    
    full_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_full_vh_idx2Idx = [
        [ full_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx, full_non_slack_gens_vh_idx2Idx, full_non_gens_vh_idx2Idx = vec_types_full_vh_idx2Idx
        

    vec_types_full_θh_idx2Idx = [
        [ full_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx, full_non_slack_gens_θh_idx2Idx, full_non_gens_θh_idx2Idx = vec_types_full_θh_idx2Idx


    full_nodes_types_vh_idx2Idx =
        (; full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx )
    
    full_nodes_types_θh_idx2Idx =
        (; full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx )

    full_nodes_types_vh_and_θh_idx2Idx =
        (; full_nodes_types_vh_idx2Idx,
         full_nodes_types_θh_idx2Idx )

    full_nodes_types_dict_vh_and_θh_idx2Idx =
        (; full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx )
   
    # --------------------------------------------------
    
    full_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    full_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_full_vh_idx2Idx_in_Idx = [
        [ full_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx_in_Idx, full_non_slack_gens_vh_idx2Idx_in_Idx, full_non_gens_vh_idx2Idx_in_Idx = vec_types_full_vh_idx2Idx_in_Idx


    vec_types_full_θh_idx2Idx_in_Idx = [
        [ full_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx_in_Idx, full_non_slack_gens_θh_idx2Idx_in_Idx, full_non_gens_θh_idx2Idx_in_Idx = vec_types_full_θh_idx2Idx_in_Idx


    full_nodes_types_vh_idx2Idx_in_Idx =
        (; full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx   )
    
    full_nodes_types_θh_idx2Idx_in_Idx =
        (; full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx )

    full_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; full_nodes_types_vh_idx2Idx_in_Idx,
         full_nodes_types_θh_idx2Idx_in_Idx )
    
    full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )

     # ------------------------------------------------   

    full_nodes_types_Idxs_idx2Idx_etc =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs,
         full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx,    
         full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx,
         full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx,
         full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx,
         full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx,    
         full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )

    # ----------------------------------------------------

    # red_vh_idx, red_θh_idx, red_vh_Idxs, red_θh_Idxs
    # non_slack_gens_θh_idx2Idx,
    # non_slack_gens_θh_idx2Idx_in_Idx,
    # non_gens_θh_idx2Idx, non_gens_θh_idx2Idx_in_Idx
    # ----------------------------------------------------


    # ----------------------------------------------------
    # Additional Idx
    # ----------------------------------------------------

    dim_intg_vh_θh_id_iq =
        [
            length( slack_gens_nodes_idx ),            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( slack_gens_nodes_idx ),               
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, intg_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_intg_vh_θh_id_iq  )

    intg_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[1]
    
    intg_non_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[2]

    intg_non_gens_nodes_vh_Idxs =
        intg_vh_θh_id_iq_IDX[3]
    
    intg_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[4]

    intg_non_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[5]
    
    intg_non_gens_nodes_θh_Idxs =
        intg_vh_θh_id_iq_IDX[6]

    intg_gen_id_Idxs =
        intg_vh_θh_id_iq_IDX[7]

    intg_gen_iq_Idxs =
        intg_vh_θh_id_iq_IDX[8]

    intg_nodes_types_vh_Idxs =
        (; intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs)

    intg_nodes_types_θh_Idxs =
        (;intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs)

    intg_nodes_idq_Idxs =
        (; intg_gen_id_Idxs,
           intg_gen_iq_Idxs )

    intg_nodes_types_vh_θh_id_iq_Idxs =
        (; intg_nodes_types_vh_Idxs,
         intg_nodes_types_θh_Idxs,
         intg_nodes_idq_Idxs)
    
     # --------------------------------------------------  
    
    intg_nodes_type_idxs =
        [slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    intg_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      vh_IDX) )
    
    intg_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_intg_vh_idx2Idx = [
        [ intg_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx, intg_non_slack_gens_vh_idx2Idx, intg_non_gens_vh_idx2Idx = vec_types_intg_vh_idx2Idx
        

    vec_types_intg_θh_idx2Idx = [
        [ intg_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx, intg_non_slack_gens_θh_idx2Idx, intg_non_gens_θh_idx2Idx = vec_types_intg_θh_idx2Idx


    intg_nodes_types_vh_idx2Idx =
        (; intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx)
    
    intg_nodes_types_θh_idx2Idx =
        (; intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx)

    intg_nodes_types_vh_and_θh_idx2Idx =
        (; intg_nodes_types_vh_idx2Idx,
         intg_nodes_types_θh_idx2Idx)

    intg_nodes_types_dict_vh_and_θh_idx2Idx =
        (; intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx )

    
     # --------------------------------------------------   
    
    intg_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    intg_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_intg_vh_idx2Idx_in_Idx = [
        [ intg_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx_in_Idx, intg_non_slack_gens_vh_idx2Idx_in_Idx, intg_non_gens_vh_idx2Idx_in_Idx = vec_types_intg_vh_idx2Idx_in_Idx


    vec_types_intg_θh_idx2Idx_in_Idx = [
        [ intg_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx_in_Idx, intg_non_slack_gens_θh_idx2Idx_in_Idx, intg_non_gens_θh_idx2Idx_in_Idx = vec_types_intg_θh_idx2Idx_in_Idx


    intg_nodes_types_vh_idx2Idx_in_Idx =
        (; intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx   )
    
    intg_nodes_types_θh_idx2Idx_in_Idx =
        (; intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx )

    intg_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; intg_nodes_types_vh_idx2Idx_in_Idx,
         intg_nodes_types_θh_idx2Idx_in_Idx )
    
    intg_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx )

     # -------------------------------------------------

    intg_nodes_types_Idxs_idx2Idx_etc =
        (;
         intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs,
         
         intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs,
         
         intg_gen_id_Idxs,
         intg_gen_iq_Idxs,
         
         intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx,
         
         intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx,
         
         intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx,
         
         intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx,
         
         intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx,
         
         intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx
         )

     # -------------------------------------------------- 

    dim_semi_vh_θh_id_iq =
        [
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, semi_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_semi_vh_θh_id_iq  )
    
    semi_non_slack_gens_vh_Idxs =
        semi_vh_θh_id_iq_IDX[1]

    semi_non_gens_nodes_vh_Idxs =
        semi_vh_θh_id_iq_IDX[2]
    
    semi_non_slack_gens_θh_Idxs =
        semi_vh_θh_id_iq_IDX[3]
    
    semi_non_gens_nodes_θh_Idxs =
        semi_vh_θh_id_iq_IDX[4]

    semi_gen_id_Idxs =
        semi_vh_θh_id_iq_IDX[5]

    semi_gen_iq_Idxs =
        semi_vh_θh_id_iq_IDX[6]

    semi_nodes_types_vh_Idxs =
        (;
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs)

    semi_nodes_types_θh_Idxs =
        (;
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs)

    semi_nodes_idq_Idxs =
        (; semi_gen_id_Idxs,
           semi_gen_iq_Idxs )

    semi_nodes_types_vh_θh_id_iq_Idxs =
        (; semi_nodes_types_vh_Idxs,
         semi_nodes_types_θh_Idxs,
         semi_nodes_idq_Idxs)
    
     # --------------------------------------------------  
    
    semi_nodes_type_idxs =
        [
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    semi_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      vh_IDX) )
    
    semi_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_semi_vh_idx2Idx = [
        [ semi_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx, semi_non_gens_vh_idx2Idx = vec_types_semi_vh_idx2Idx
        

    vec_types_semi_θh_idx2Idx = [
        [ semi_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx, semi_non_gens_θh_idx2Idx = vec_types_semi_θh_idx2Idx


    semi_nodes_types_vh_idx2Idx =
        (; 
         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx)
    
    semi_nodes_types_θh_idx2Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx)

    semi_nodes_types_vh_and_θh_idx2Idx =
        (; semi_nodes_types_vh_idx2Idx,
         semi_nodes_types_θh_idx2Idx)

    semi_nodes_types_dict_vh_and_θh_idx2Idx =
        (; semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx )

    
     # --------------------------------------------------
    
    semi_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    semi_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )
    
    vec_types_semi_vh_idx2Idx_in_Idx = [
        [ semi_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx_in_Idx, semi_non_gens_vh_idx2Idx_in_Idx = vec_types_semi_vh_idx2Idx_in_Idx


    vec_types_semi_θh_idx2Idx_in_Idx = [
        [ semi_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx_in_Idx, semi_non_gens_θh_idx2Idx_in_Idx = vec_types_semi_θh_idx2Idx_in_Idx


    semi_nodes_types_vh_idx2Idx_in_Idx =
        (;
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx   )
    
    semi_nodes_types_θh_idx2Idx_in_Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx )

    semi_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; semi_nodes_types_vh_idx2Idx_in_Idx,
         semi_nodes_types_θh_idx2Idx_in_Idx )
    
    semi_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx )

     # -------------------------------------------------

    semi_nodes_types_Idxs_idx2Idx_etc =
        (;
         
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs,
                  
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs,
         
         semi_gen_id_Idxs,
         semi_gen_iq_Idxs,

         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx,

         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx,
         
         semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx,
         
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx,
         
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx,
         
         semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx
         )
        
    # --------------------------------------------------

    return (; slack_ur_ui_Idx_in_state,
            non_slack_ur_ui_Idx_in_state,
            ur_ui_Idx_in_state,
            nodes_u_Idx,
            ur_idx,
            ui_idx,
            ur_ui_idx,
            vh_θh_idx,
            ur_ui_IDX,
            ur_IDX,
            ui_IDX,
            ir_IDX,
            ii_IDX,
            vh_IDX,
            θh_IDX,
            slack_bus_idx,
            gens_idx,
            red_vh_θh_idx,
            red_vh_idx,
            red_θh_idx,
            red_vh_Idxs,
            red_θh_Idxs,
            non_slack_gens_θh_idx2Idx,
            non_slack_gens_θh_idx2Idx_in_Idx,
            non_gens_θh_idx2Idx,
            non_gens_θh_idx2Idx_in_Idx,
            loc_load_exist,
            full_nodes_types_Idxs_idx2Idx_etc,
            intg_nodes_types_Idxs_idx2Idx_etc,
            semi_nodes_types_Idxs_idx2Idx_etc )
    
end



function get_integrated_pf_and_dyn_idx_and_Idx(netd)
    
    #------------------------------------------            
    _, _, edges_orientation,  nodes_node_idx_and_incident_edges_other_node_idx = get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------    

    (; slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx) =
         get_net_nodes_type_idxs( netd)
    
    # ----------------------------------------------------
    
    vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
                netd.nodes) )

    
    vec_Idx_gens_nodes_ω_ed_dash_eq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_ω_ed_dash_eq_dash_Idxs_in_plants(
                netd.nodes) )
    
    
    # ----------------------------------------------------

    vec_Idx_ra_Xd_Xq_Xd_dash_Xq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param_dims(
                netd.nodes
                ;some_param =
                    [:ra, :X_d, :X_q,
                     :X_d_dash, :X_q_dash]  ))


    vec_Idx_ra_Xd_dash_Xq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param_dims(
                netd.nodes
                ;some_param =
                    [ :ra, :X_d_dash, :X_q_dash ] ))


    vec_Idx_ra_Xd_Xq =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param_dims(
                netd.nodes
                ;some_param =
                    [ :ra, :X_d, :X_q ] ))
    
    # ------------------------------------------------
    
    vec_Idx_P_Q_nodes =
        get_flattened_to_components_vector_var_Idx(
            get_nodes_some_param_dims(
                netd.nodes
                ; some_param = [:P, :Q] ) )

    
    vec_Idx_P_Q_gens =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P, :Q ] ) )

    
    vec_Idx_P_Q_non_gens =
        get_flattened_to_components_vector_var_Idx(
            get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P, :Q ] ) )

    
    #------------------------------------------
    
    gens_nodes_with_loc_loads_some_param_dims =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param =
                [ :loc_P, :loc_Q  ] )
    
    vec_Idx_P_Q_gens_loc_load =
        gens_nodes_with_loc_loads_some_param_dims == 0 ?
        nothing :
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_with_loc_loads_some_param_dims(
                netd.nodes
                ; some_param = [ :loc_P, :loc_Q  ] ) )

    #------------------------------------------
    #------------------------------------------

    dim_P_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :Q ] ) )
    
    dim_P_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :Q ] ) )       
    
    dim_δ_ed_eq_pf =
        length( [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes)...;] )

    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_P ] ) == 0 ? 0 :
        length( get_gens_nodes_with_loc_loads_some_param_dims( netd.nodes; some_param = [ :loc_P ] )  )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
        length( get_gens_nodes_with_loc_loads_some_param_dims( netd.nodes; some_param = [ :loc_Q  ] )  )
    
    #----------------------------------------
    #----------------------------------------
    
    dim_with_loc_load_with_δ_ed_NL_pf_para =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_δ_ed_eq_pf,
          dim_P_g_loc_load,
          dim_Q_g_loc_load  ]

    _,_, with_loc_load_with_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_with_loc_load_with_δ_ed_NL_pf_para )

    with_loc_load_with_δ_ed_P_gens_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[1]

    with_loc_load_with_δ_ed_Q_gens_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[2]

    with_loc_load_with_δ_ed_P_non_gens_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[3]

    with_loc_load_with_δ_ed_Q_non_gens_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[4]

    with_loc_load_δ_ed_eq_pf_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[5]

    with_loc_load_with_δ_ed_P_g_loc_load_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[6]

    with_loc_load_with_δ_ed_Q_g_loc_load_NL_para_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[7]

    with_loc_load_with_δ_ed_NL_pf_para_Idxs = (
        ; with_loc_load_with_δ_ed_P_gens_NL_para_Idxs,
        with_loc_load_with_δ_ed_Q_gens_NL_para_Idxs,
        with_loc_load_with_δ_ed_P_non_gens_NL_para_Idxs,
        with_loc_load_with_δ_ed_Q_non_gens_NL_para_Idxs,
        with_loc_load_δ_ed_eq_pf_NL_para_Idxs,
        with_loc_load_with_δ_ed_P_g_loc_load_NL_para_Idxs,
        with_loc_load_with_δ_ed_Q_g_loc_load_NL_para_Idxs)

    #----------------------------------------


    dim_with_loc_load_no_δ_ed_NL_pf_para =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens, 
          dim_P_g_loc_load,
          dim_Q_g_loc_load  ]

    _,_, with_loc_load_no_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_with_loc_load_no_δ_ed_NL_pf_para )

    with_loc_load_no_δ_ed_P_gens_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[1]

    with_loc_load_no_δ_ed_Q_gens_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[2]

    with_loc_load_no_δ_ed_P_non_gens_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[3]

    with_loc_load_no_δ_ed_Q_non_gens_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[4]

    with_loc_load_no_δ_ed_P_g_loc_load_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[5]

    with_loc_load_no_δ_ed_Q_g_loc_load_NL_para_Idxs =
        with_loc_load_no_δ_ed_CA_pf_para_Idxs[6]

    with_loc_load_no_δ_ed_NL_pf_para_Idxs = (
        ; with_loc_load_no_δ_ed_P_gens_NL_para_Idxs,
        with_loc_load_no_δ_ed_Q_gens_NL_para_Idxs,
        with_loc_load_no_δ_ed_P_non_gens_NL_para_Idxs,
        with_loc_load_no_δ_ed_Q_non_gens_NL_para_Idxs,
        with_loc_load_no_δ_ed_P_g_loc_load_NL_para_Idxs,
        with_loc_load_no_δ_ed_Q_g_loc_load_NL_para_Idxs)

    #----------------------------------------


    dim_no_loc_load_with_δ_ed_NL_pf_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_δ_ed_eq_pf ]

    _,_, no_loc_load_with_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_no_loc_load_with_δ_ed_NL_pf_para  )

    no_loc_load_with_δ_ed_P_gens_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[1]

    no_loc_load_with_δ_ed_Q_gens_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[2]

    no_loc_load_with_δ_ed_P_non_gens_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[3]

    no_loc_load_with_δ_ed_Q_non_gens_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[4]

    no_loc_load_with_δ_ed_δ_ed_eq_pf_NL_para_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[5]

    no_loc_load_with_δ_ed_NL_pf_para_Idxs = (
        ; no_loc_load_with_δ_ed_P_gens_NL_para_Idxs,
        no_loc_load_with_δ_ed_Q_gens_NL_para_Idxs,
        no_loc_load_with_δ_ed_P_non_gens_NL_para_Idxs,
        no_loc_load_with_δ_ed_Q_non_gens_NL_para_Idxs,
        no_loc_load_with_δ_ed_δ_ed_eq_pf_NL_para_Idxs )

    
    #----------------------------------------


    dim_no_loc_load_no_δ_ed_NL_pf_para =
        [ dim_P_gens, dim_Q_gens, dim_P_non_gens,
          dim_Q_non_gens  ]

    _,_, no_loc_load_no_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_no_loc_load_no_δ_ed_NL_pf_para )

    no_loc_load_no_δ_ed_P_gens_NL_para_Idxs =
        no_loc_load_no_δ_ed_CA_pf_para_Idxs[1]

    no_loc_load_no_δ_ed_Q_gens_NL_para_Idxs =
        no_loc_load_no_δ_ed_CA_pf_para_Idxs[2]

    no_loc_load_no_δ_ed_P_non_gens_NL_para_Idxs =
        no_loc_load_no_δ_ed_CA_pf_para_Idxs[3]

    no_loc_load_no_δ_ed_Q_non_gens_NL_para_Idxs =
        no_loc_load_no_δ_ed_CA_pf_para_Idxs[4]

    no_loc_load_no_δ_ed_NL_pf_para_Idxs = (
        ; no_loc_load_no_δ_ed_P_gens_NL_para_Idxs,
        no_loc_load_no_δ_ed_Q_gens_NL_para_Idxs,
        no_loc_load_no_δ_ed_P_non_gens_NL_para_Idxs,
        no_loc_load_no_δ_ed_Q_non_gens_NL_para_Idxs)
    
    #------------------------------------------

    dims_δ_ω_ed_dash_eq = length.( get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants( netd.nodes))

    _,_, δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq =
        create_size_offset_Idx( dims_δ_ω_ed_dash_eq )

    #------------------------------------------
    
    dicts_net_to_streamlined_idx =
        get_dict_net_to_streamlined_idx(netd)
    
    #------------------------------------------

    NL_pf_para_Idxs = (
        ;  with_loc_load_with_δ_ed_NL_pf_para_Idxs,
        with_loc_load_no_δ_ed_NL_pf_para_Idxs,
        no_loc_load_with_δ_ed_NL_pf_para_Idxs,
        no_loc_load_no_δ_ed_NL_pf_para_Idxs,
        δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq )
    
    #------------------------------------------

    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx
    
    net_comp_type_idx =
        (; slack_gens_nodes_idx,
           non_slack_gens_nodes_idx,
           gens_nodes_idx,
           load_nodes_idx,
           transmission_nodes_idx,
           non_gens_nodes_idx,
           load_trans_nodes_idx,
           gens_with_loc_load_idx,
         all_nodes_idx,
         loc_load_exist,
         NL_pf_para_Idxs,
         dicts_net_to_streamlined_idx )
    
    #------------------------------------------
    
    vec_Idx =
        (; vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash,
         vec_Idx_gens_nodes_ω_ed_dash_eq_dash,
         vec_Idx_ra_Xd_Xq_Xd_dash_Xq_dash,
         vec_Idx_ra_Xd_dash_Xq_dash,
         vec_Idx_ra_Xd_Xq,
         vec_Idx_P_Q_nodes,
         vec_Idx_P_Q_gens,
         vec_Idx_P_Q_non_gens,
         vec_Idx_P_Q_gens_loc_load,
         loc_load_exist,
         NL_pf_para_Idxs )
    
    #------------------------------------------    

    return (;
            net_comp_type_idx,
            vec_Idx,
            loc_load_exist,
            NL_pf_para_Idxs,
            dicts_net_to_streamlined_idx)        
end


#-----------------------------------------------------
# a model a_model
#-----------------------------------------------------

function get_a_model_integrated_pf_var_idx_and_Idx(
    netd;
    Idxs_type =
        Idxs_type,
    no_control_device =
        false )

    #------------------------------------------        

    (; slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx) =
         get_net_nodes_type_idxs( netd)
    
    #------------------------------------------        

    (; slack_ur_ui_Idx_in_state, 
     non_slack_ur_ui_Idx_in_state, 
     ur_ui_Idx_in_state, 
     nodes_u_Idx,
     nodes_u_Idx_in_ranges) =
         get_a_model_integrated_u_ur_ui_Idx_in_state(
             netd ;
             Idxs_type =
                 Idxs_type,
             no_control_device =
                 no_control_device )

    #------------------------------------------        
            
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------    
    
    vh_θh_idx  = ur_ui_idx
   
    #------------------------------------------    

    ur_ui_dims   =
        [ length(ur_idx), length(ui_idx) ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    # ----------------------------------------------------
    
    red_vh_θh_idx =
        [ setdiff(vh_IDX, gens_idx)...;
          setdiff(θh_IDX, θh_IDX[slack_bus_idx])... ]
    
    # ----------------------------------------------------

    red_vh_idx = setdiff(vh_IDX, gens_idx)
    
    red_θh_idx =
        setdiff(θh_IDX, θh_IDX[slack_bus_idx])

    # ----------------------------------------------------

    red_vh_θh_dims =
        [ length( red_vh_idx ), length( red_θh_idx ) ]    

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(red_vh_θh_dims )

    red_vh_Idxs = red_vh_θh_IDX[1]
    
    red_θh_Idxs = red_vh_θh_IDX[2]
    
     # -------------------------------------------------- 
    
     dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      1:length(red_θh_Idxs) ) )

     non_slack_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_slack_gens_nodes_idx ]

     non_slack_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_slack_gens_nodes_idx ]


     non_gens_θh_idx2Idx =
         [ dict_θh_idx2Idx[idx]
          for idx in
              non_gens_nodes_idx ]

     non_gens_θh_idx2Idx_in_Idx =
         [ dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              non_gens_nodes_idx ]
    
    # ----------------------------------------------------

    full_vh_θh_dims =
        [ length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx ),
          length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx )
          ]
    
     _, _, full_vh_θh_IDX =
        create_size_offset_Idx( full_vh_θh_dims )

    full_slack_gens_vh_Idxs =
        full_vh_θh_IDX[1]
    
    full_non_slack_gens_vh_Idxs =
        full_vh_θh_IDX[2]

    full_non_gens_nodes_vh_Idxs =
        full_vh_θh_IDX[3]
    
    full_slack_gens_θh_Idxs =
        full_vh_θh_IDX[4]

    full_non_slack_gens_θh_Idxs =
        full_vh_θh_IDX[5]
    
    full_non_gens_nodes_θh_Idxs =
        full_vh_θh_IDX[6]

    full_nodes_types_vh_Idxs =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs )

    full_nodes_types_θh_Idxs =
        (;full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs)

    full_nodes_types_vh_and_θh_Idxs =
        (; full_nodes_types_vh_Idxs,
         full_nodes_types_θh_Idxs)
    
     # --------------------------------------------------  
    
    full_nodes_type_idxs =
        [ slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    full_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      vh_IDX) )
    
    full_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_full_vh_idx2Idx = [
        [ full_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx, full_non_slack_gens_vh_idx2Idx, full_non_gens_vh_idx2Idx = vec_types_full_vh_idx2Idx
        

    vec_types_full_θh_idx2Idx = [
        [ full_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx, full_non_slack_gens_θh_idx2Idx, full_non_gens_θh_idx2Idx = vec_types_full_θh_idx2Idx


    full_nodes_types_vh_idx2Idx =
        (; full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx )
    
    full_nodes_types_θh_idx2Idx =
        (; full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx )

    full_nodes_types_vh_and_θh_idx2Idx =
        (; full_nodes_types_vh_idx2Idx,
         full_nodes_types_θh_idx2Idx )

    full_nodes_types_dict_vh_and_θh_idx2Idx =
        (; full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx )
   
    # --------------------------------------------------
    
    full_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    full_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_full_vh_idx2Idx_in_Idx = [
        [ full_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx_in_Idx, full_non_slack_gens_vh_idx2Idx_in_Idx, full_non_gens_vh_idx2Idx_in_Idx = vec_types_full_vh_idx2Idx_in_Idx


    vec_types_full_θh_idx2Idx_in_Idx = [
        [ full_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx_in_Idx, full_non_slack_gens_θh_idx2Idx_in_Idx, full_non_gens_θh_idx2Idx_in_Idx = vec_types_full_θh_idx2Idx_in_Idx


    full_nodes_types_vh_idx2Idx_in_Idx =
        (; full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx   )
    
    full_nodes_types_θh_idx2Idx_in_Idx =
        (; full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx )

    full_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; full_nodes_types_vh_idx2Idx_in_Idx,
         full_nodes_types_θh_idx2Idx_in_Idx )
    
    full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )


    # ------------------------------------------------   

    full_gens_vh_Idxs =
        first(full_slack_gens_vh_Idxs):last(full_non_slack_gens_vh_Idxs)
    
    full_gens_θh_Idxs =
        first(full_slack_gens_θh_Idxs):last(full_non_slack_gens_θh_Idxs)

    # ------------------------------------------------   

    full_nodes_types_Idxs_idx2Idx_etc =
        (;full_gens_vh_Idxs,
         full_gens_θh_Idxs,
         
         full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         
         full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs,
         
         full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx,
         
         full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx,
         
         full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx,
         
         full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx,
         
         full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx,
         
         full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )

    # ----------------------------------------------------

    dim_intg_vh_θh_id_iq =
        [
            length( slack_gens_nodes_idx ),            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( slack_gens_nodes_idx ),               
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, intg_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_intg_vh_θh_id_iq  )

    intg_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[1]
    
    intg_non_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[2]

    intg_non_gens_nodes_vh_Idxs =
        intg_vh_θh_id_iq_IDX[3]
    
    intg_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[4]

    intg_non_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[5]
    
    intg_non_gens_nodes_θh_Idxs =
        intg_vh_θh_id_iq_IDX[6]

    intg_gen_id_Idxs =
        intg_vh_θh_id_iq_IDX[7]

    intg_gen_iq_Idxs =
        intg_vh_θh_id_iq_IDX[8]

    intg_nodes_types_vh_Idxs =
        (; intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs)

    intg_nodes_types_θh_Idxs =
        (;intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs)

    intg_nodes_idq_Idxs =
        (; intg_gen_id_Idxs,
           intg_gen_iq_Idxs )

    intg_nodes_types_vh_θh_id_iq_Idxs =
        (; intg_nodes_types_vh_Idxs,
         intg_nodes_types_θh_Idxs,
         intg_nodes_idq_Idxs)
    
     # ------------------------------------------------  
    
    intg_nodes_type_idxs =
        [slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    intg_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      vh_IDX) )
    
    intg_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_intg_vh_idx2Idx = [
        [ intg_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx, intg_non_slack_gens_vh_idx2Idx, intg_non_gens_vh_idx2Idx = vec_types_intg_vh_idx2Idx
        

    vec_types_intg_θh_idx2Idx = [
        [ intg_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx, intg_non_slack_gens_θh_idx2Idx, intg_non_gens_θh_idx2Idx = vec_types_intg_θh_idx2Idx


    intg_nodes_types_vh_idx2Idx =
        (; intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx)
    
    intg_nodes_types_θh_idx2Idx =
        (; intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx)

    intg_nodes_types_vh_and_θh_idx2Idx =
        (; intg_nodes_types_vh_idx2Idx,
         intg_nodes_types_θh_idx2Idx)

    intg_nodes_types_dict_vh_and_θh_idx2Idx =
        (; intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx )

    
     # ------------------------------------------------   
    
    intg_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    intg_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_intg_vh_idx2Idx_in_Idx = [
        [ intg_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx_in_Idx, intg_non_slack_gens_vh_idx2Idx_in_Idx, intg_non_gens_vh_idx2Idx_in_Idx = vec_types_intg_vh_idx2Idx_in_Idx


    vec_types_intg_θh_idx2Idx_in_Idx = [
        [ intg_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx_in_Idx, intg_non_slack_gens_θh_idx2Idx_in_Idx, intg_non_gens_θh_idx2Idx_in_Idx = vec_types_intg_θh_idx2Idx_in_Idx


    intg_nodes_types_vh_idx2Idx_in_Idx =
        (; intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx   )
    
    intg_nodes_types_θh_idx2Idx_in_Idx =
        (; intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx )

    intg_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; intg_nodes_types_vh_idx2Idx_in_Idx,
         intg_nodes_types_θh_idx2Idx_in_Idx )
    
    intg_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx )

     # -------------------------------------------------

    # This assumption is wrong when slack node is not
    # the first node
    
    intg_gens_vh_Idxs =
        first(intg_slack_gens_vh_Idxs):last(intg_non_slack_gens_vh_Idxs)
    
    intg_gens_θh_Idxs =
        first(intg_slack_gens_θh_Idxs):last(intg_non_slack_gens_θh_Idxs)

     # -------------------------------------------------
    
    intg_nodes_types_Idxs_idx2Idx_etc =
        (;intg_gens_vh_Idxs,
         intg_gens_θh_Idxs,
         
         intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs,
         
         intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs,
         
         intg_gen_id_Idxs,
         intg_gen_iq_Idxs,
         
         intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx,
         
         intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx,
         
         intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx,
         
         intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx,
         
         intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx,
         
         intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx
         )

     # -------------------------------------------------- 

    dim_semi_vh_θh_id_iq =
        [
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, semi_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_semi_vh_θh_id_iq  )
    
    semi_non_slack_gens_vh_Idxs =
        semi_vh_θh_id_iq_IDX[1]

    semi_non_gens_nodes_vh_Idxs =
        semi_vh_θh_id_iq_IDX[2]
    
    semi_non_slack_gens_θh_Idxs =
        semi_vh_θh_id_iq_IDX[3]
    
    semi_non_gens_nodes_θh_Idxs =
        semi_vh_θh_id_iq_IDX[4]

    semi_gen_id_Idxs =
        semi_vh_θh_id_iq_IDX[5]

    semi_gen_iq_Idxs =
        semi_vh_θh_id_iq_IDX[6]

    semi_nodes_types_vh_Idxs =
        (;
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs)

    semi_nodes_types_θh_Idxs =
        (;
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs)

    semi_nodes_idq_Idxs =
        (; semi_gen_id_Idxs,
           semi_gen_iq_Idxs )

    semi_nodes_types_vh_θh_id_iq_Idxs =
        (; semi_nodes_types_vh_Idxs,
         semi_nodes_types_θh_Idxs,
         semi_nodes_idq_Idxs)
    
     # --------------------------------------------------  
    
    semi_nodes_type_idxs =
        [
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    semi_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      vh_IDX) )
    
    semi_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_semi_vh_idx2Idx = [
        [ semi_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx, semi_non_gens_vh_idx2Idx = vec_types_semi_vh_idx2Idx
        

    vec_types_semi_θh_idx2Idx = [
        [ semi_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx, semi_non_gens_θh_idx2Idx = vec_types_semi_θh_idx2Idx


    semi_nodes_types_vh_idx2Idx =
        (; 
         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx)
    
    semi_nodes_types_θh_idx2Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx)

    semi_nodes_types_vh_and_θh_idx2Idx =
        (; semi_nodes_types_vh_idx2Idx,
         semi_nodes_types_θh_idx2Idx)

    semi_nodes_types_dict_vh_and_θh_idx2Idx =
        (; semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx )

    
     # --------------------------------------------------
    
    semi_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    semi_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )
    
    vec_types_semi_vh_idx2Idx_in_Idx = [
        [ semi_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx_in_Idx, semi_non_gens_vh_idx2Idx_in_Idx = vec_types_semi_vh_idx2Idx_in_Idx


    vec_types_semi_θh_idx2Idx_in_Idx = [
        [ semi_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx_in_Idx, semi_non_gens_θh_idx2Idx_in_Idx = vec_types_semi_θh_idx2Idx_in_Idx


    semi_nodes_types_vh_idx2Idx_in_Idx =
        (;
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx   )
    
    semi_nodes_types_θh_idx2Idx_in_Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx )

    semi_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; semi_nodes_types_vh_idx2Idx_in_Idx,
         semi_nodes_types_θh_idx2Idx_in_Idx )
    
    semi_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx )

     # -------------------------------------------------

    semi_nodes_types_Idxs_idx2Idx_etc =
        (;
         
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs,
                  
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs,
         
         semi_gen_id_Idxs,
         semi_gen_iq_Idxs,

         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx,

         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx,
         
         semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx,
         
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx,
         
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx,
         
         semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx
         )
        
    # --------------------------------------------------

    return (; slack_ur_ui_Idx_in_state,
            non_slack_ur_ui_Idx_in_state,
            ur_ui_Idx_in_state,
            nodes_u_Idx,
            
            ur_idx,
            ui_idx,
            ur_ui_idx,
            vh_θh_idx,
            
            ur_ui_IDX,
            ur_IDX,
            ui_IDX,
            ir_IDX,
            ii_IDX,
            vh_IDX,
            θh_IDX,
            
            slack_bus_idx,
            gens_idx,
            
            red_vh_θh_idx,
            red_vh_idx,
            red_θh_idx,
            red_vh_Idxs,
            red_θh_Idxs,
            
            non_slack_gens_θh_idx2Idx,
            non_slack_gens_θh_idx2Idx_in_Idx,
            
            non_gens_θh_idx2Idx,
            non_gens_θh_idx2Idx_in_Idx,
            
            loc_load_exist,
            
            full_nodes_types_Idxs_idx2Idx_etc,
            intg_nodes_types_Idxs_idx2Idx_etc,
            semi_nodes_types_Idxs_idx2Idx_etc )
    
end


function get_a_model_integrated_pf_red_full_semi_intg_var_Idx(
    netd;
    Idxs_type =
        Idxs_type,
    no_control_device =
        false )

    #------------------------------------------        

    (; slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx) =
         get_net_nodes_type_idxs( netd)

    net_nodes_type_idxs =
        (; slack_bus_idx,
         gens_idx,
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         loc_load_exist,
         load_nodes_idx,
         transmission_nodes_idx,
         non_gens_nodes_idx,
         load_trans_nodes_idx,
         all_nodes_idx,
         non_slack_gens_and_non_gens_idx )
    
    #------------------------------------------        

    (; slack_ur_ui_Idx_in_state, 
     non_slack_ur_ui_Idx_in_state, 
     ur_ui_Idx_in_state, 
     nodes_u_Idx,
     nodes_u_Idx_in_ranges) =
         get_a_model_integrated_u_ur_ui_Idx_in_state(
             netd ;
             Idxs_type =
                 Idxs_type,
             no_control_device =
                 no_control_device )

    #------------------------------------------        
            
    ur_idx     = first.( ur_ui_Idx_in_state )
    
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
    
    #------------------------------------------
    
    vh_idx = ur_idx

    θh_idx = ui_idx
    
    vh_θh_idx  = ur_ui_idx

    ur_ui_vh_θh_idxs =
        (; ur_idx,
         ui_idx,
         vh_idx,
         θh_idx)

    #------------------------------------------    

    ur_ui_dims   =
        [ length(ur_idx), length(ui_idx) ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    ur_ui_vh_θh_Idxs =
        (;ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX )
    
    # ----------------------------------------------------

    non_gens_vh_idx = setdiff(vh_IDX, gens_idx)

    non_slack_gens_θh_idx =
        θh_IDX[ non_slack_gens_nodes_idx ]

    non_gens_θh_idx =
        θh_IDX[non_gens_nodes_idx]

    red_θh_idx =
        setdiff(θh_IDX, θh_IDX[slack_bus_idx])

    # ----------------------------------------------------
    
    red_vh_θh_idx =
        [  non_gens_vh_idx...;
           non_slack_gens_θh_idx...;
           non_gens_θh_idx... ]

    # ----------------------------------------------------

    red_var_comps_idxs =
        (; non_gens_vh_idx,
         non_slack_gens_θh_idx,
         non_gens_θh_idx )

    # ----------------------------------------------------

    red_vh_θh_dims =
        length.([ non_gens_vh_idx,
                   non_slack_gens_θh_idx,
                   non_gens_θh_idx  ] )  

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(
            red_vh_θh_dims )

    red_vh_Idxs = red_non_gens_vh_Idxs =
        red_vh_θh_IDX[1]
    
    red_non_slack_gens_θh_Idxs =
        red_vh_θh_IDX[2]
    
    red_non_gens_θh_Idxs   =
        red_vh_θh_IDX[3]

    red_θh_Idxs =
        first(red_non_slack_gens_θh_Idxs):last(red_non_gens_θh_Idxs)
    
     # -------------------------------------------------- 
    
     red_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     red_dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
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

    red_types_Idxs_etc =
        (;
         red_non_gens_vh_Idxs,
         red_non_slack_gens_θh_Idxs,
         red_non_gens_θh_Idxs,
         
         red_vh_Idxs,
         red_θh_Idxs,

         red_non_slack_gens_θh_idx2Idx,
         red_non_gens_θh_idx2Idx,
         red_non_slack_gens_θh_idx2Idx_in_Idx,            
         red_non_gens_θh_idx2Idx_in_Idx,

         red_dict_θh_idx2Idx,
         red_dict_θh_idx2Idx_in_Idx,

         red_vh_θh_IDX
         )    
    
    # ----------------------------------------------------
    
    full_vh_θh_dims =
        [ length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx ),
          length( slack_gens_nodes_idx ),
          length( non_slack_gens_nodes_idx ),
          length( non_gens_nodes_idx )
          ]
    
     _, _, full_vh_θh_IDX =
        create_size_offset_Idx( full_vh_θh_dims )

    full_slack_gens_vh_Idxs =
        full_vh_θh_IDX[1]
    
    full_non_slack_gens_vh_Idxs =
        full_vh_θh_IDX[2]

    full_non_gens_nodes_vh_Idxs =
        full_vh_θh_IDX[3]
    
    full_slack_gens_θh_Idxs =
        full_vh_θh_IDX[4]

    full_non_slack_gens_θh_Idxs =
        full_vh_θh_IDX[5]
    
    full_non_gens_nodes_θh_Idxs =
        full_vh_θh_IDX[6]

    full_nodes_types_vh_Idxs =
        (; full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs )

    full_nodes_types_θh_Idxs =
        (;full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs)

    full_nodes_types_vh_and_θh_Idxs =
        (; full_nodes_types_vh_Idxs,
         full_nodes_types_θh_Idxs)
    
     # -------------------------------------------------  
    
    full_nodes_type_idxs =
        [ slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    full_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      vh_IDX) )
    
    full_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_full_vh_idx2Idx = [
        [ full_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx, full_non_slack_gens_vh_idx2Idx, full_non_gens_vh_idx2Idx = vec_types_full_vh_idx2Idx
        

    vec_types_full_θh_idx2Idx = [
        [ full_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx, full_non_slack_gens_θh_idx2Idx, full_non_gens_θh_idx2Idx = vec_types_full_θh_idx2Idx


    full_nodes_types_vh_idx2Idx =
        (; full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx )
    
    full_nodes_types_θh_idx2Idx =
        (; full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx )

    full_nodes_types_vh_and_θh_idx2Idx =
        (; full_nodes_types_vh_idx2Idx,
         full_nodes_types_θh_idx2Idx )

    full_nodes_types_dict_vh_and_θh_idx2Idx =
        (; full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx )
   
    # --------------------------------------------------
    
    full_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    full_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [full_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_full_vh_idx2Idx_in_Idx = [
        [ full_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_vh_idx2Idx_in_Idx, full_non_slack_gens_vh_idx2Idx_in_Idx, full_non_gens_vh_idx2Idx_in_Idx = vec_types_full_vh_idx2Idx_in_Idx


    vec_types_full_θh_idx2Idx_in_Idx = [
        [ full_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            full_nodes_type_idxs ]

    full_slack_gens_θh_idx2Idx_in_Idx, full_non_slack_gens_θh_idx2Idx_in_Idx, full_non_gens_θh_idx2Idx_in_Idx = vec_types_full_θh_idx2Idx_in_Idx


    full_nodes_types_vh_idx2Idx_in_Idx =
        (; full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx   )
    
    full_nodes_types_θh_idx2Idx_in_Idx =
        (; full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx )

    full_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; full_nodes_types_vh_idx2Idx_in_Idx,
         full_nodes_types_θh_idx2Idx_in_Idx )
    
    full_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx )


    # ------------------------------------------------   

    # This assumption is wrong when slack node is not
    # the first node

    full_gens_vh_Idxs =
        first(full_slack_gens_vh_Idxs):last(full_non_slack_gens_vh_Idxs)
    
    full_gens_θh_Idxs =
        first(full_slack_gens_θh_Idxs):last(full_non_slack_gens_θh_Idxs)

    # ------------------------------------------------   

    full_types_Idxs_etc =
        (;full_gens_vh_Idxs,
         full_gens_θh_Idxs,
         
         full_slack_gens_vh_Idxs,
         full_non_slack_gens_vh_Idxs,
         full_non_gens_nodes_vh_Idxs,
         
         full_slack_gens_θh_Idxs,
         full_non_slack_gens_θh_Idxs,
         full_non_gens_nodes_θh_Idxs,
         
         full_slack_gens_vh_idx2Idx,
         full_non_slack_gens_vh_idx2Idx,
         full_non_gens_vh_idx2Idx,
         
         full_slack_gens_θh_idx2Idx,
         full_non_slack_gens_θh_idx2Idx,
         full_non_gens_θh_idx2Idx,
         
         full_dict_vh_idx2Idx,
         full_dict_θh_idx2Idx,
         
         full_slack_gens_vh_idx2Idx_in_Idx,
         full_non_slack_gens_vh_idx2Idx_in_Idx,
         full_non_gens_vh_idx2Idx_in_Idx,
         
         full_slack_gens_θh_idx2Idx_in_Idx,
         full_non_slack_gens_θh_idx2Idx_in_Idx,
         full_non_gens_θh_idx2Idx_in_Idx,
         
         full_dict_vh_idx2Idx_in_Idx,
         full_dict_θh_idx2Idx_in_Idx,

         full_vh_θh_IDX
         )

    # ----------------------------------------------------

    dim_intg_vh_θh_id_iq =
        [
            length( slack_gens_nodes_idx ),            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( slack_gens_nodes_idx ),               
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, intg_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_intg_vh_θh_id_iq  )

    intg_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[1]
    
    intg_non_slack_gens_vh_Idxs =
        intg_vh_θh_id_iq_IDX[2]

    intg_non_gens_nodes_vh_Idxs =
        intg_vh_θh_id_iq_IDX[3]
    
    intg_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[4]

    intg_non_slack_gens_θh_Idxs =
        intg_vh_θh_id_iq_IDX[5]
    
    intg_non_gens_nodes_θh_Idxs =
        intg_vh_θh_id_iq_IDX[6]

    intg_gen_id_Idxs =
        intg_vh_θh_id_iq_IDX[7]

    intg_gen_iq_Idxs =
        intg_vh_θh_id_iq_IDX[8]

    intg_nodes_types_vh_Idxs =
        (; intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs)

    intg_nodes_types_θh_Idxs =
        (;intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs)

    intg_nodes_idq_Idxs =
        (; intg_gen_id_Idxs,
           intg_gen_iq_Idxs )

    intg_nodes_types_vh_θh_id_iq_Idxs =
        (; intg_nodes_types_vh_Idxs,
         intg_nodes_types_θh_Idxs,
         intg_nodes_idq_Idxs)
    
     # -------------------------------------------------  
    
    intg_nodes_type_idxs =
        [slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    intg_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      vh_IDX) )
    
    intg_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_intg_vh_idx2Idx = [
        [ intg_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx, intg_non_slack_gens_vh_idx2Idx, intg_non_gens_vh_idx2Idx = vec_types_intg_vh_idx2Idx
        

    vec_types_intg_θh_idx2Idx = [
        [ intg_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx, intg_non_slack_gens_θh_idx2Idx, intg_non_gens_θh_idx2Idx = vec_types_intg_θh_idx2Idx


    intg_nodes_types_vh_idx2Idx =
        (; intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx)
    
    intg_nodes_types_θh_idx2Idx =
        (; intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx)

    intg_nodes_types_vh_and_θh_idx2Idx =
        (; intg_nodes_types_vh_idx2Idx,
         intg_nodes_types_θh_idx2Idx)

    intg_nodes_types_dict_vh_and_θh_idx2Idx =
        (; intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx )

    
     # ------------------------------------------------   
    
    intg_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    intg_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [intg_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )

    
    vec_types_intg_vh_idx2Idx_in_Idx = [
        [ intg_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_vh_idx2Idx_in_Idx, intg_non_slack_gens_vh_idx2Idx_in_Idx, intg_non_gens_vh_idx2Idx_in_Idx = vec_types_intg_vh_idx2Idx_in_Idx


    vec_types_intg_θh_idx2Idx_in_Idx = [
        [ intg_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    intg_slack_gens_θh_idx2Idx_in_Idx, intg_non_slack_gens_θh_idx2Idx_in_Idx, intg_non_gens_θh_idx2Idx_in_Idx = vec_types_intg_θh_idx2Idx_in_Idx


    intg_nodes_types_vh_idx2Idx_in_Idx =
        (; intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx   )
    
    intg_nodes_types_θh_idx2Idx_in_Idx =
        (; intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx )

    intg_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; intg_nodes_types_vh_idx2Idx_in_Idx,
         intg_nodes_types_θh_idx2Idx_in_Idx )
    
    intg_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx )

     # -------------------------------------------------

    # This assumption is wrong when slack node is not
    # the first node
    
    intg_gens_vh_Idxs =
        first(intg_slack_gens_vh_Idxs):last(intg_non_slack_gens_vh_Idxs)
    
    intg_gens_θh_Idxs =
        first(intg_slack_gens_θh_Idxs):last(intg_non_slack_gens_θh_Idxs)

     # -------------------------------------------------
    
    intg_types_Idxs_etc =
        (;intg_gens_vh_Idxs,
         intg_gens_θh_Idxs,
         
         intg_slack_gens_vh_Idxs,
         intg_non_slack_gens_vh_Idxs,
         intg_non_gens_nodes_vh_Idxs,
         
         intg_slack_gens_θh_Idxs,
         intg_non_slack_gens_θh_Idxs,
         intg_non_gens_nodes_θh_Idxs,
         
         intg_gen_id_Idxs,
         intg_gen_iq_Idxs,
         
         intg_slack_gens_vh_idx2Idx,
         intg_non_slack_gens_vh_idx2Idx,
         intg_non_gens_vh_idx2Idx,
         
         intg_slack_gens_θh_idx2Idx,
         intg_non_slack_gens_θh_idx2Idx,
         intg_non_gens_θh_idx2Idx,
         
         intg_dict_vh_idx2Idx,
         intg_dict_θh_idx2Idx,
         
         intg_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_slack_gens_vh_idx2Idx_in_Idx,
         intg_non_gens_vh_idx2Idx_in_Idx,
         
         intg_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_slack_gens_θh_idx2Idx_in_Idx,
         intg_non_gens_θh_idx2Idx_in_Idx,
         
         intg_dict_vh_idx2Idx_in_Idx,
         intg_dict_θh_idx2Idx_in_Idx,

         intg_vh_θh_id_iq_IDX
         )

     # -------------------------------------------------- 

    dim_semi_vh_θh_id_iq =
        [
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( non_slack_gens_nodes_idx ),
            length( non_gens_nodes_idx ),
            
            length( gens_nodes_idx ),
            length( gens_nodes_idx ) ]
    
     _, _, semi_vh_θh_id_iq_IDX =
        create_size_offset_Idx( dim_semi_vh_θh_id_iq  )
    
    semi_non_slack_gens_vh_Idxs =
        semi_vh_θh_id_iq_IDX[1]

    semi_non_gens_nodes_vh_Idxs =
        semi_vh_θh_id_iq_IDX[2]
    
    semi_non_slack_gens_θh_Idxs =
        semi_vh_θh_id_iq_IDX[3]
    
    semi_non_gens_nodes_θh_Idxs =
        semi_vh_θh_id_iq_IDX[4]

    semi_gen_id_Idxs =
        semi_vh_θh_id_iq_IDX[5]

    semi_gen_iq_Idxs =
        semi_vh_θh_id_iq_IDX[6]

    semi_nodes_types_vh_Idxs =
        (;
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs)

    semi_nodes_types_θh_Idxs =
        (;
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs)

    semi_nodes_idq_Idxs =
        (; semi_gen_id_Idxs,
           semi_gen_iq_Idxs )

    semi_nodes_types_vh_θh_id_iq_Idxs =
        (; semi_nodes_types_vh_Idxs,
         semi_nodes_types_θh_Idxs,
         semi_nodes_idq_Idxs)
    
     # --------------------------------------------------  
    
    semi_nodes_type_idxs =
        [
         non_slack_gens_nodes_idx,
         non_gens_nodes_idx ]


    semi_dict_vh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      vh_IDX) )
    
    semi_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      θh_IDX ) )
    
    vec_types_semi_vh_idx2Idx = [
        [ semi_dict_vh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx, semi_non_gens_vh_idx2Idx = vec_types_semi_vh_idx2Idx
        

    vec_types_semi_θh_idx2Idx = [
        [ semi_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx, semi_non_gens_θh_idx2Idx = vec_types_semi_θh_idx2Idx


    semi_nodes_types_vh_idx2Idx =
        (; 
         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx)
    
    semi_nodes_types_θh_idx2Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx)

    semi_nodes_types_vh_and_θh_idx2Idx =
        (; semi_nodes_types_vh_idx2Idx,
         semi_nodes_types_θh_idx2Idx)

    semi_nodes_types_dict_vh_and_θh_idx2Idx =
        (; semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx )

    
     # --------------------------------------------------
    
    semi_dict_vh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length(vh_IDX) ) )

    
    semi_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( [semi_nodes_type_idxs...;],
                      1:length( θh_IDX ) ) )
    
    vec_types_semi_vh_idx2Idx_in_Idx = [
        [ semi_dict_vh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_vh_idx2Idx_in_Idx, semi_non_gens_vh_idx2Idx_in_Idx = vec_types_semi_vh_idx2Idx_in_Idx


    vec_types_semi_θh_idx2Idx_in_Idx = [
        [ semi_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            semi_nodes_type_idxs ]

     semi_non_slack_gens_θh_idx2Idx_in_Idx, semi_non_gens_θh_idx2Idx_in_Idx = vec_types_semi_θh_idx2Idx_in_Idx


    semi_nodes_types_vh_idx2Idx_in_Idx =
        (;
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx   )
    
    semi_nodes_types_θh_idx2Idx_in_Idx =
        (; 
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx )

    semi_nodes_types_vh_and_θh_idx2Idx_in_Idx =
        (; semi_nodes_types_vh_idx2Idx_in_Idx,
         semi_nodes_types_θh_idx2Idx_in_Idx )
    
    semi_nodes_types_dict_vh_and_θh_idx2Idx_in_Idx =
        (; semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx )

     # -------------------------------------------------

    semi_types_Idxs_etc =
        (;
         
         semi_non_slack_gens_vh_Idxs,
         semi_non_gens_nodes_vh_Idxs,
                  
         semi_non_slack_gens_θh_Idxs,
         semi_non_gens_nodes_θh_Idxs,
         
         semi_gen_id_Idxs,
         semi_gen_iq_Idxs,

         semi_non_slack_gens_vh_idx2Idx,
         semi_non_gens_vh_idx2Idx,

         semi_non_slack_gens_θh_idx2Idx,
         semi_non_gens_θh_idx2Idx,
         
         semi_dict_vh_idx2Idx,
         semi_dict_θh_idx2Idx,
         
         semi_non_slack_gens_vh_idx2Idx_in_Idx,
         semi_non_gens_vh_idx2Idx_in_Idx,
         
         semi_non_slack_gens_θh_idx2Idx_in_Idx,
         semi_non_gens_θh_idx2Idx_in_Idx,
         
         semi_dict_vh_idx2Idx_in_Idx,
         semi_dict_θh_idx2Idx_in_Idx,

         semi_vh_θh_id_iq_IDX
         )
        
    # --------------------------------------------------

    return (
        ; slack_ur_ui_Idx_in_state,
        non_slack_ur_ui_Idx_in_state,
        ur_ui_Idx_in_state,
        nodes_u_Idx,
        nodes_u_Idx_in_ranges,
        

        ur_ui_vh_θh_idxs,
        ur_ui_vh_θh_Idxs,
        
        net_nodes_type_idxs,

        red_types_Idxs_etc,
        full_types_Idxs_etc,
        intg_types_Idxs_etc,
        semi_types_Idxs_etc )
    
end



function get_a_model_integrated_pf_sta_dyn_para_idx(
    netd;
    Idxs_type =
        Idxs_type,
    no_control_device =
        false )

    #------------------------------------------        

    (; slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx) =
         get_net_nodes_type_idxs( netd)

    #------------------------------------------        

    nodes_size =
        length(gens_nodes_idx) +
        length(non_gens_nodes_idx)
    
    #------------------------------------------        

    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx
    
    nodes_types_idxs =
        (;slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_with_loc_load_idx )
    
    #------------------------------------------        

    (; dict_n2s_slack_gens_idx,
     dict_n2s_non_slack_gens_idx,
     dict_n2s_gens_idx,
     dict_n2s_non_gens_idx,
     dict_n2s_load_idx,
     dict_n2s_gens_with_loc_load_idxs,
     dict_n2s_transmission_idxs,
     dict_n2s_all_nodes_idx )  =
         get_dict_net_to_streamlined_idx(
             netd)

    n2s_slack_gens_idx =
        dict_n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx =
        dict_n2s_gens_idx
    
    n2s_non_gens_idx =
        dict_n2s_non_gens_idx
    
    n2s_load_idx =
        dict_n2s_load_idx
    
    n2s_gens_with_loc_load_idxs =
        dict_n2s_gens_with_loc_load_idxs
    
    n2s_transmission_idxs =
        dict_n2s_transmission_idxs

    n2s_all_nodes_idx =
        dict_n2s_all_nodes_idx
    
    n2s_idxs =
        (; n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,         
         n2s_all_nodes_idx )
    
    #------------------------------------------        
    
    (; slack_ur_ui_Idx_in_state, 
     non_slack_ur_ui_Idx_in_state, 
     ur_ui_Idx_in_state, 
     nodes_u_Idx,
     nodes_u_Idx_in_ranges) =
         get_a_model_integrated_u_ur_ui_Idx_in_state(
             netd ;
             Idxs_type =
                 Idxs_type,
             no_control_device =
                 no_control_device )
            
    u_ur_ui_Idx_in_state =
        (;slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges)
    
    #------------------------------------------        
            
    ur_idx     = first.( ur_ui_Idx_in_state )
    ui_idx     = last.(  ur_ui_Idx_in_state )

    ur_ui_idx  = [ ur_idx...; ui_idx... ]
   
    #------------------------------------------
    
    vh_idx = ur_idx

    θh_idx = ui_idx
    
    vh_θh_idx  = ur_ui_idx

    ur_ui_vh_θh_idxs =
        (; ur_idx,
         ui_idx,
         vh_idx,
         θh_idx)

    #------------------------------------------    

    ur_ui_dims   =
        [ nodes_size, nodes_size ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    ur_ui_vh_θh_Idxs =
        (;ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX )
    
    # ----------------------------------------------------

    non_gens_vh_idx = setdiff(vh_IDX, gens_idx)

    non_slack_gens_θh_idx =
        θh_IDX[ non_slack_gens_nodes_idx ]

    non_gens_θh_idx =
        θh_IDX[non_gens_nodes_idx]

    red_θh_idx =
        setdiff(θh_IDX, θh_IDX[slack_bus_idx])

    # ----------------------------------------------------
    # ----------------------------------------------------
    
    red_vh_θh_idx =
        [  non_gens_vh_idx...;
           non_slack_gens_θh_idx...;
           non_gens_θh_idx... ]

    # ----------------------------------------------------

    red_var_comps_idxs =
        (; non_gens_vh_idx,
         non_slack_gens_θh_idx,
         non_gens_θh_idx )

    # ----------------------------------------------------

    red_vh_θh_dims =
        length.([ non_gens_vh_idx,
                   non_slack_gens_θh_idx,
                   non_gens_θh_idx  ] )  

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(
            red_vh_θh_dims )

    red_vh_Idxs = red_non_gens_vh_Idxs =
        red_vh_θh_IDX[1]
    
    red_non_slack_gens_θh_Idxs =
        red_vh_θh_IDX[2]
    
    red_non_gens_θh_Idxs   =
        red_vh_θh_IDX[3]

    red_θh_Idxs =
        first(red_non_slack_gens_θh_Idxs):last(red_non_gens_θh_Idxs)
    
     # -------------------------------------------------- 
    
     red_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     red_dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
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

    red_types_Idxs_etc =
        (;
         red_non_gens_vh_Idxs,
         red_non_slack_gens_θh_Idxs,
         red_non_gens_θh_Idxs,
         
         red_vh_Idxs,
         red_θh_Idxs,

         red_non_slack_gens_θh_idx2Idx,
         red_non_gens_θh_idx2Idx,
         red_non_slack_gens_θh_idx2Idx_in_Idx,            
         red_non_gens_θh_idx2Idx_in_Idx,

         red_dict_θh_idx2Idx,
         red_dict_θh_idx2Idx_in_Idx,

         red_vh_θh_IDX
         )    
    
    # ----------------------------------------------------
    # ----------------------------------------------------

    dim_P_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :Q ] ) )
    
    dim_P_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
            ; some_param = [ :Q ] ) )

    
    dim_δ_ed_eq_pf = length(
        [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
            netd.nodes)...;] )
    
    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes;
            some_param = [ :loc_P ] ) == 0 ? 0 :
                length(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_P ] ) )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
                length(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_Q  ] )  )
    
    #----------------------------------------

    dim_pf_sta_PQ_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_P_g_loc_load,
          dim_Q_g_loc_load ]

    _,_, pf_sta_PQ_para_IDX =
        create_size_offset_Idx(
            dim_pf_sta_PQ_para )

    P_gens_sta_para_Idxs =
         pf_sta_PQ_para_IDX[1]

    Q_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[2]

    P_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[3]

    Q_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[4]

    P_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[5]

    Q_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[6]

    PQ_sta_para_Idxs = (
        ; P_gens_sta_para_Idxs,
        Q_gens_sta_para_Idxs,
        P_non_gens_sta_para_Idxs,
        Q_non_gens_sta_para_Idxs,
        P_g_loc_load_sta_para_Idxs,
        Q_g_loc_load_sta_para_Idxs )

    #----------------------------------------

    dim_pf_dyn_PQ_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_δ_ed_eq_pf,
          dim_P_g_loc_load,
          dim_Q_g_loc_load ]

    _,_, pf_dyn_PQ_para_IDX =
        create_size_offset_Idx(
            dim_pf_dyn_PQ_para )

    P_gens_dyn_para_Idxs =
         pf_dyn_PQ_para_IDX[1]

    Q_gens_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[2]

    P_non_gens_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[3]

    Q_non_gens_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[4]

    δ_ed_eq_pf_para_Idxs =
        pf_dyn_PQ_para_IDX[5]
    
    P_g_loc_load_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[6]

    Q_g_loc_load_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[7]

    PQ_dyn_para_Idxs = (
        ; P_gens_dyn_para_Idxs,
        Q_gens_dyn_para_Idxs,
        P_non_gens_dyn_para_Idxs,
        Q_non_gens_dyn_para_Idxs,
        δ_ed_eq_pf_para_Idxs,
        P_g_loc_load_dyn_para_Idxs,
        Q_g_loc_load_dyn_para_Idxs,
        pf_dyn_PQ_para_IDX)

    #----------------------------------------

    # PQ_dyn_para_Idxs = (
    #     ; P_gens_dyn_para_Idxs,
    #     Q_gens_dyn_para_Idxs,
    #     P_non_gens_dyn_para_Idxs,
    #     Q_non_gens_dyn_para_Idxs,
    #     δ_ed_eq_pf_para_Idxs,
    #     P_g_loc_load_dyn_para_Idxs,
    #     Q_g_loc_load_dyn_para_Idxs )
    
    # PQ_sta_para_Idxs = (
    #     ; P_gens_sta_para_Idxs,
    #     Q_gens_sta_para_Idxs,
    #     P_non_gens_sta_para_Idxs,
    #     Q_non_gens_sta_para_Idxs,
    #     P_g_loc_load_sta_para_Idxs,
    #     Q_g_loc_load_sta_para_Idxs )
    
    # nodes_types_idxs =
    #     (;slack_gens_nodes_idx,
    #      non_slack_gens_nodes_idx,
    #      gens_nodes_idx,
    #      non_gens_nodes_idx,
    #      gens_with_loc_load_idx)
             

    # n2s_idxs =
    #     (; n2s_slack_gens_idx,
    #      n2s_non_slack_gens_idx,
    #      n2s_gens_idx,
    #      n2s_non_gens_idx,
    #      n2s_gens_with_loc_load_idxs )
            
    # u_ur_ui_Idx_in_state =
    #     (;slack_ur_ui_Idx_in_state,
    #      non_slack_ur_ui_Idx_in_state,
    #      ur_ui_Idx_in_state,
    #      nodes_u_Idx, nodes_u_Idx_in_ranges)
    
    # ur_ui_vh_θh_idxs =
    #     (; ur_idx,
    #      ui_idx,
    #      vh_idx,
    #      θh_idx)
    
    
    # ur_ui_vh_θh_Idxs =
    #     (;ur_IDX,
    #      ui_IDX,
    #      vh_IDX,
    #      θh_IDX )

    # red_idxs =
    #     (;non_gens_vh_idx,
    #      non_slack_gens_θh_idx,
    #      non_gens_θh_idx)
    
    
    # red_Idxs =
    #     (; non_gens_vh_Idxs,
    #      non_slack_gens_θh_Idxs,
    #      non_gens_θh_Idxs,
    #      red_θh_Idxs )
    
    # idx2Idx = (; non_slack_gens_θh_idx2Idx,
    #         non_gens_θh_idx2Idx,
    #         non_slack_gens_θh_idx2Idx_in_Idx,            
    #         non_gens_θh_idx2Idx_in_Idx)
    
    # --------------------------------------------------

    return (; red_types_Idxs_etc,
            PQ_sta_para_Idxs,
            PQ_dyn_para_Idxs,
            nodes_types_idxs,
            n2s_idxs,
            u_ur_ui_Idx_in_state,            
            ur_ui_vh_θh_idxs,
            ur_ui_vh_θh_Idxs)
    
end



function get_a_model_integrated_sta_pf_vars_and_paras_idx(
    netd )

    #------------------------------------------        

    (; slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx) =
         get_net_nodes_type_idxs( netd)

    #------------------------------------------        

    nodes_size =
        length(gens_nodes_idx) +
        length(non_gens_nodes_idx)
    
    #------------------------------------------
    
    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx

    nodes_types_idxs =
        (;slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,         
         gens_with_loc_load_idx,
         all_nodes_idx)
    
    #------------------------------------------        

    (; dict_n2s_slack_gens_idx,
     dict_n2s_non_slack_gens_idx,
     dict_n2s_gens_idx,
     dict_n2s_non_gens_idx,
     dict_n2s_load_idx,
     dict_n2s_gens_with_loc_load_idxs,
     dict_n2s_transmission_idxs,
     dict_n2s_all_nodes_idx )  =
         get_dict_net_to_streamlined_idx(
             netd)

    n2s_slack_gens_idx =
        dict_n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx =
        dict_n2s_gens_idx
    
    n2s_non_gens_idx =
        dict_n2s_non_gens_idx
    
    n2s_load_idx =
        dict_n2s_load_idx
    
    n2s_gens_with_loc_load_idxs =
        dict_n2s_gens_with_loc_load_idxs
    
    n2s_transmission_idxs =
        dict_n2s_transmission_idxs

    n2s_all_nodes_idx =
        dict_n2s_all_nodes_idx

    n2s_idxs =
        (; n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx )
    
    #------------------------------------------        

    ur_ui_dims   =
        [ nodes_size, nodes_size ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs( ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

    ur_ui_vh_θh_Idxs =
        (;ur_IDX,
         ui_IDX,
         vh_IDX,
         θh_IDX )
    
    # ----------------------------------------------------

    non_gens_vh_idx = setdiff(vh_IDX, gens_idx)

    non_slack_gens_θh_idx =
        θh_IDX[ non_slack_gens_nodes_idx ]

    non_gens_θh_idx =
        θh_IDX[non_gens_nodes_idx]

    red_θh_idx =
        setdiff(θh_IDX, θh_IDX[slack_bus_idx])

    # ----------------------------------------------------
    
    red_vh_θh_idx =
        [  non_gens_vh_idx...;
           non_slack_gens_θh_idx...;
           non_gens_θh_idx... ]

    # ----------------------------------------------------

    red_var_comps_idxs =
        (; non_gens_vh_idx,
         non_slack_gens_θh_idx,
         non_gens_θh_idx )

    # ----------------------------------------------------

    red_vh_θh_dims =
        length.([ non_gens_vh_idx,
                   non_slack_gens_θh_idx,
                   non_gens_θh_idx  ] )  

    _, _, red_vh_θh_IDX =
        create_size_offset_Idx(
            red_vh_θh_dims )

    red_vh_Idxs = red_non_gens_vh_Idxs =
        red_vh_θh_IDX[1]
    
    red_non_slack_gens_θh_Idxs =
        red_vh_θh_IDX[2]
    
    red_non_gens_θh_Idxs   =
        red_vh_θh_IDX[3]

    red_θh_Idxs =
        first(red_non_slack_gens_θh_Idxs):last(red_non_gens_θh_Idxs)
    
     # -------------------------------------------------- 
    
     red_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

     red_dict_θh_idx2Idx_in_Idx = OrderedDict{Int64, Int64}(
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

    red_types_Idxs_etc =
        (;
         red_non_gens_vh_Idxs,
         red_non_slack_gens_θh_Idxs,
         red_non_gens_θh_Idxs,
         
         red_vh_Idxs,
         red_θh_Idxs,

         red_non_slack_gens_θh_idx2Idx,
         red_non_gens_θh_idx2Idx,
         red_non_slack_gens_θh_idx2Idx_in_Idx,            
         red_non_gens_θh_idx2Idx_in_Idx,

         red_dict_θh_idx2Idx,
         red_dict_θh_idx2Idx_in_Idx,

         red_vh_θh_IDX
         )    
    
    # ----------------------------------------------------

    dim_P_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :Q ] ) )
    
    dim_P_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
            ; some_param = [ :Q ] ) )

    
    dim_δ_ed_eq_pf = length(
        [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
            netd.nodes)...;] )
    
    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes;
            some_param = [ :loc_P ] ) == 0 ? 0 :
                length(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_P ] ) )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
                length(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_Q  ] )  )
    
    #----------------------------------------

    dim_pf_sta_PQ_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_P_g_loc_load,
          dim_Q_g_loc_load ]

    _,_, pf_sta_PQ_para_IDX =
        create_size_offset_Idx(
            dim_pf_sta_PQ_para )

    P_gens_sta_para_Idxs =
         pf_sta_PQ_para_IDX[1]

    Q_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[2]

    P_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[3]

    Q_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[4]

    P_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[5]

    Q_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[6]

    PQ_sta_para_Idxs = (
        ; P_gens_sta_para_Idxs,
        Q_gens_sta_para_Idxs,
        P_non_gens_sta_para_Idxs,
        Q_non_gens_sta_para_Idxs,
        P_g_loc_load_sta_para_Idxs,
        Q_g_loc_load_sta_para_Idxs )

    #----------------------------------------
    #----------------------------------------
    
    # --------------------------------------------------

    return (;red_types_Idxs_etc,
            PQ_sta_para_Idxs,
            nodes_types_idxs,
            n2s_idxs )
    
end




function get_a_model_integrated_dyn_pf_PQ_δ_etc_Idxs(
    netd;
    δ_etc = [:δ, :ω, :ed_dash, :eq_dash],
    δ_etc_first = false )


    dim_P_gens =
        sum(get_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :P] ))

    dim_Q_gens =
        sum(get_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :Q] ))

    dim_P_non_gens =
        sum(get_non_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :P ] ))

    dim_Q_non_gens =
        sum(get_non_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :Q ] ))
    
    dim_δ_etc_pf = sum(get_gens_nodes_some_state_vars_dims(
        netd.nodes;
        some_state_vars =
            δ_etc ))
    
    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes;
            some_param = [ :loc_P ] ) == 0 ? 0 :
                sum(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_P ] ) )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
                sum(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_Q  ] )  )
    
    #----------------------------------------
    #----------------------------------------


    if δ_etc_first == false

        dim_pf_PQ_δ_etc_para  =
            [ dim_P_gens,
              dim_Q_gens,
              dim_P_non_gens,
              dim_Q_non_gens,
              dim_δ_etc_pf,
              dim_P_g_loc_load,
              dim_Q_g_loc_load]

        _,_, pf_PQ_δ_etc_para_IDX =
            create_size_offset_Idx(
                dim_pf_PQ_δ_etc_para )

        dyn_P_gens_Idxs =
             pf_PQ_δ_etc_para_IDX[1]

        dyn_Q_gens_Idxs =
            pf_PQ_δ_etc_para_IDX[2]

        dyn_P_non_gens_Idxs =
            pf_PQ_δ_etc_para_IDX[3]

        dyn_Q_non_gens_Idxs =
            pf_PQ_δ_etc_para_IDX[4]

        dyn_δ_ed_eq_pf_Idxs =
            pf_PQ_δ_etc_para_IDX[5]

        dyn_P_g_loc_load_Idxs =
            pf_PQ_δ_etc_para_IDX[6]

        dyn_Q_g_loc_load_Idxs =
            pf_PQ_δ_etc_para_IDX[7]

        dyn_pf_PQ_δ_etc_Idxs = (
            ;dyn_P_gens_Idxs,
            dyn_Q_gens_Idxs,
            dyn_P_non_gens_Idxs,
            dyn_Q_non_gens_Idxs,
            dyn_δ_ed_eq_pf_Idxs,
            dyn_P_g_loc_load_Idxs,
            dyn_Q_g_loc_load_Idxs )
        
        return (; dyn_pf_PQ_δ_etc_Idxs, )
        
    else

        dim_pf_PQ_δ_etc_para  =
            [ dim_δ_etc_pf,
              dim_P_gens,
              dim_Q_gens,
              dim_P_non_gens,
              dim_Q_non_gens,
              dim_P_g_loc_load,
              dim_Q_g_loc_load]

        _,_, pf_PQ_δ_etc_para_IDX =
            create_size_offset_Idx(
                dim_pf_PQ_δ_etc_para )

        dyn_δ_ed_eq_pf_Idxs =
            pf_PQ_δ_etc_para_IDX[1]

        dyn_P_gens_Idxs =
             pf_PQ_δ_etc_para_IDX[2]

        dyn_Q_gens_Idxs =
            pf_PQ_δ_etc_para_IDX[3]

        dyn_P_non_gens_Idxs =
            pf_PQ_δ_etc_para_IDX[4]

        dyn_Q_non_gens_Idxs =
            pf_PQ_δ_etc_para_IDX[5]

        dyn_P_g_loc_load_Idxs =
            pf_PQ_δ_etc_para_IDX[6]

        dyn_Q_g_loc_load_Idxs =
            pf_PQ_δ_etc_para_IDX[7]

        dyn_δ_etc_pf_PQ_Idxs = (
            ;
            dyn_δ_ed_eq_pf_Idxs,
            dyn_P_gens_Idxs,
            dyn_Q_gens_Idxs,
            dyn_P_non_gens_Idxs,
            dyn_Q_non_gens_Idxs,     
            dyn_P_g_loc_load_Idxs,
            dyn_Q_g_loc_load_Idxs )
        
        return (; dyn_δ_etc_pf_PQ_Idxs, )
        
    end

    #----------------------------------------

    
end


function get_intra_pf_flat_vh_θh_wt_dyn_pf_flat_para_Idx(
    netd;
    δ_etc = [:δ, :ω, :ed_dash, :eq_dash],
    flat_vh_θh_first = true )


    ur_ui_Idx_in_system = [
        get_nodes_ur_ui_indices_in_system(netd )...; ]
    
    dim_flat_vh_θh =
        length( ur_ui_Idx_in_system  )
    
    dim_P_gens =
        sum(get_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :P] ))

    dim_Q_gens =
        sum(get_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :Q] ))

    dim_P_non_gens =
        sum(get_non_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :P ] ))

    dim_Q_non_gens =
        sum(get_non_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :Q ] ))
    
    dim_δ_etc_pf =
        sum(get_gens_nodes_some_state_vars_dims(
        netd.nodes;
        some_state_vars =
            δ_etc ))
    
    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes;
            some_param = [ :loc_P ] ) == 0 ? 0 :
                sum(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_P ] ) )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
                sum(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_Q  ] )  )

   dim_pf_PQ_δ_etc_para  =
       sum([ dim_P_gens,
         dim_Q_gens,
         dim_P_non_gens,
         dim_Q_non_gens,
         dim_δ_etc_pf,
         dim_P_g_loc_load,
             dim_Q_g_loc_load])
    
    #----------------------------------------
    #----------------------------------------

    if flat_vh_θh_first == true

        dims_flat_vh_θh_wt_pf_PQ_δ_etc_para =
            [dim_flat_vh_θh,
             dim_pf_PQ_δ_etc_para ]

        _,_, flat_vh_θh_wt_pf_PQ_δ_etc_para_IDX =
            create_size_offset_Idx(
                dims_flat_vh_θh_wt_pf_PQ_δ_etc_para  )
        
        intra_flat_vh_θh_Idx, intra_dyn_pf_flat_para_Idx =
            flatted_vh_θh_wt_pf_PQ_δ_etc_para_IDX
    else
        dims_pf_PQ_δ_etc_para_wt_flat_vh_θh =
            [ dim_pf_PQ_δ_etc_para,
              dim_flat_vh_θh ]

        _,_, pf_PQ_δ_etc_para_wt_flat_vh_θh_IDX =
            create_size_offset_Idx(
                dims_pf_PQ_δ_etc_para_wt_flat_vh_θh  )
        
        intra_dyn_pf_flat_para_Idx, intra_flat_vh_θh_Idx =
            pf_PQ_δ_etc_para_wt_flat_vh_θh_IDX 
        
    end
    

    return (; intra_flat_vh_θh_Idx,
            intra_dyn_pf_flat_para_Idx)

    
end



function get_intra_pf_flat_ur_ui_wt_dyn_pf_flat_para_Idx(
    netd;
    δ_etc = [:δ, :ω, :ed_dash, :eq_dash],
    flat_ur_ui_first = true )

    intra_flat_vh_θh_Idx, intra_dyn_pf_flat_para_Idx =
        get_intra_pf_flat_vh_θh_wt_dyn_pf_flat_para_Idx(
            netd; δ_etc = δ_etc,
            flat_vh_θh_first = flat_ur_ui_first )

    intra_flat_ur_ui_Idx = intra_flat_vh_θh_Idx
    
    return (; intra_flat_ur_ui_Idx, intra_dyn_pf_flat_para_Idx )

    
end



function get_a_model_integrated_PQ_sta_para_Idxs(
    netd)


    dim_P_gens =
        sum(get_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :P] ))

    dim_Q_gens =
        sum(get_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :Q] ))

    dim_P_non_gens =
        sum(get_non_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :P ] ))

    dim_Q_non_gens =
        sum(get_non_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :Q ] ))
    
    dim_δ_ed_eq_pf =
        sum(get_gens_nodes_some_state_vars_dims(
        netd.nodes;
        some_state_vars =
            [:δ, :ω, :ed_dash, :eq_dash] ))
    
    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes;
            some_param = [ :loc_P ] ) == 0 ? 0 :
                sum(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_P ] ) )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
                sum(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_Q  ] )  )
    
    #----------------------------------------

    dim_pf_sta_PQ_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_P_g_loc_load,
          dim_Q_g_loc_load ]

    _,_, pf_sta_PQ_para_IDX =
        create_size_offset_Idx(
            dim_pf_sta_PQ_para )

    P_gens_sta_para_Idxs =
         pf_sta_PQ_para_IDX[1]

    Q_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[2]

    P_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[3]

    Q_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[4]

    P_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[5]

    Q_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[6]

    PQ_sta_para_Idxs = (
        ; P_gens_sta_para_Idxs,
        Q_gens_sta_para_Idxs,
        P_non_gens_sta_para_Idxs,
        Q_non_gens_sta_para_Idxs,
        P_g_loc_load_sta_para_Idxs,
        Q_g_loc_load_sta_para_Idxs )
    
    #----------------------------------------
    #----------------------------------------

    return PQ_sta_para_Idxs

    
end


function get_a_model_integrated_pf_PQ_δ_ω_ed_dash_eq_dash_Idxs(
    netd)


    dim_P_gens =
        sum(get_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :P] ))

    dim_Q_gens =
        sum(get_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :Q] ))

    dim_P_non_gens =
        sum(get_non_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :P ] ))

    dim_Q_non_gens =
        sum(get_non_gens_nodes_some_param_dims(
            netd.nodes; some_param = [ :Q ] ))
    
    dim_δ_ed_eq_pf = sum(get_gens_nodes_some_state_vars_dims(
        netd.nodes;
        some_state_vars =
            [:δ, :ω, :ed_dash, :eq_dash] ))
    
    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes;
            some_param = [ :loc_P ] ) == 0 ? 0 :
                sum(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_P ] ) )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
                sum(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_Q  ] )  )
    
    #----------------------------------------

    dim_pf_sta_PQ_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_P_g_loc_load,
          dim_Q_g_loc_load ]

    _,_, pf_sta_PQ_para_IDX =
        create_size_offset_Idx(
            dim_pf_sta_PQ_para )

    P_gens_sta_para_Idxs =
         pf_sta_PQ_para_IDX[1]

    Q_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[2]

    P_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[3]

    Q_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[4]

    P_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[5]

    Q_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[6]

    PQ_sta_para_Idxs = (
        ; P_gens_sta_para_Idxs,
        Q_gens_sta_para_Idxs,
        P_non_gens_sta_para_Idxs,
        Q_non_gens_sta_para_Idxs,
        P_g_loc_load_sta_para_Idxs,
        Q_g_loc_load_sta_para_Idxs )
    
    #----------------------------------------

    dim_dyn_pf_no_ll_PQ_δ_ed_eq_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_δ_ed_eq_pf]

    _,_, dyn_pf_no_ll_PQ_δ_ed_eq_para_IDX =
        create_size_offset_Idx(
            dim_dyn_pf_no_ll_PQ_δ_ed_eq_para )

    P_gens_dyn_pf_no_ll_para_Idxs =
         dyn_pf_no_ll_PQ_δ_ed_eq_para_IDX[1]

    Q_gens_dyn_pf_no_ll_para_Idxs =
        dyn_pf_no_ll_PQ_δ_ed_eq_para_IDX[2]

    P_non_gens_dyn_pf_no_ll_para_Idxs =
        dyn_pf_no_ll_PQ_δ_ed_eq_para_IDX[3]

    Q_non_gens_dyn_pf_no_ll_para_Idxs =
        dyn_pf_no_ll_PQ_δ_ed_eq_para_IDX[4]

    δ_ω_ed_eq_ed_dyn_pf_no_ll_para_Idxs =
        dyn_pf_no_ll_PQ_δ_ed_eq_para_IDX[5]

    # Q_g_loc_load_sta_para_Idxs =
    #     pf_sta_PQ_para_IDX[6]

    dyn_pf_no_ll_PQ_δ_ed_eq_para_Idxs = (
        ; P_gens_dyn_pf_no_ll_para_Idxs,
        Q_gens_dyn_pf_no_ll_para_Idxs,
        P_non_gens_dyn_pf_no_ll_para_Idxs,
        Q_non_gens_dyn_pf_no_ll_para_Idxs,
        δ_ω_ed_eq_ed_dyn_pf_no_ll_para_Idxs )

    #----------------------------------------


    dim_dyn_pf_wt_ll_PQ_δ_ed_eq_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_δ_ed_eq_pf,
          dim_P_g_loc_load,
          dim_Q_g_loc_load]

    _,_, dyn_pf_wt_ll_PQ_δ_ed_eq_para_IDX =
        create_size_offset_Idx(
            dim_dyn_pf_wt_ll_PQ_δ_ed_eq_para )

    P_gens_dyn_pf_wt_ll_para_Idxs =
         dyn_pf_wt_ll_PQ_δ_ed_eq_para_IDX[1]

    Q_gens_dyn_pf_wt_ll_para_Idxs =
        dyn_pf_wt_ll_PQ_δ_ed_eq_para_IDX[2]

    P_non_gens_dyn_pf_wt_ll_para_Idxs =
        dyn_pf_wt_ll_PQ_δ_ed_eq_para_IDX[3]

    Q_non_gens_dyn_pf_wt_ll_para_Idxs =
        dyn_pf_wt_ll_PQ_δ_ed_eq_para_IDX[4]

    δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs =
        dyn_pf_wt_ll_PQ_δ_ed_eq_para_IDX[5]

    P_g_loc_load_dyn_pf_wt_ll_para_Idxs =
        dyn_pf_wt_ll_PQ_δ_ed_eq_para_IDX[6]

    Q_g_loc_load_dyn_pf_wt_ll_para_Idxs =
        dyn_pf_wt_ll_PQ_δ_ed_eq_para_IDX[7]
    
    
    dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs = (
        ; P_gens_dyn_pf_no_ll_para_Idxs,
        Q_gens_dyn_pf_wt_ll_para_Idxs,
        P_non_gens_dyn_pf_wt_ll_para_Idxs,
        Q_non_gens_dyn_pf_wt_ll_para_Idxs,
        δ_ω_ed_eq_ed_dyn_pf_wt_ll_para_Idxs,
        P_g_loc_load_dyn_pf_wt_ll_para_Idxs,
        Q_g_loc_load_dyn_pf_wt_ll_para_Idxs)
        
        
    #----------------------------------------

    return (; PQ_sta_para_Idxs,
            dyn_pf_no_ll_PQ_δ_ed_eq_para_Idxs,
            dyn_pf_wt_ll_PQ_δ_ed_eq_para_Idxs)

    
end

#------------------------------------------

function get_a_model_integrated_pf_vars_and_para_idx(
    netd;
    Idxs_type =
        Idxs_type,
    no_control_device =
        false )

    #------------------------------------------        

    (; slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx) =
         get_net_nodes_type_idxs( netd)


    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx
    
    nodes_types_idxs =
        (;slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_with_loc_load_idx,
         all_nodes_idx)
    
    #------------------------------------------        

    (; dict_n2s_slack_gens_idx,
     dict_n2s_non_slack_gens_idx,
     dict_n2s_gens_idx,
     dict_n2s_non_gens_idx,
     dict_n2s_load_idx,
     dict_n2s_gens_with_loc_load_idxs,
     dict_n2s_transmission_idxs,
     dict_n2s_all_nodes_idx )  =
         get_dict_net_to_streamlined_idx(
             netd)

    n2s_slack_gens_idx =
        dict_n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx =
        dict_n2s_gens_idx
    
    n2s_non_gens_idx =
        dict_n2s_non_gens_idx
    
    n2s_load_idx =
        dict_n2s_load_idx
    
    n2s_gens_with_loc_load_idxs =
        dict_n2s_gens_with_loc_load_idxs
    
    n2s_transmission_idxs =
        dict_n2s_transmission_idxs

    n2s_all_nodes_idx =
        dict_n2s_all_nodes_idx

    n2s_idxs =
        (; n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx )
    
    #------------------------------------------        
    
    (;
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges,

     ur_ui_vh_θh_idxs,
     ur_ui_vh_θh_Idxs,

     net_nodes_type_idxs,     

     red_types_Idxs_etc,
     full_types_Idxs_etc,
     intg_types_Idxs_etc,
     semi_types_Idxs_etc ) =
         get_a_model_integrated_pf_red_full_semi_intg_var_Idx(
             netd;
             Idxs_type =
                 Idxs_type,
             no_control_device =
                 no_control_device )

            
    u_ur_ui_Idx_in_state =
        (;slack_ur_ui_Idx_in_state,
         non_slack_ur_ui_Idx_in_state,
         ur_ui_Idx_in_state,
         nodes_u_Idx,
         nodes_u_Idx_in_ranges)
    
    #------------------------------------------        

   (;
    red_non_gens_vh_Idxs,
    red_non_slack_gens_θh_Idxs,
    red_non_gens_θh_Idxs,

    red_vh_Idxs,
    red_θh_Idxs,

    red_non_slack_gens_θh_idx2Idx,
    red_non_gens_θh_idx2Idx,
    red_non_slack_gens_θh_idx2Idx_in_Idx,            
    red_non_gens_θh_idx2Idx_in_Idx,

    red_dict_θh_idx2Idx,
    red_dict_θh_idx2Idx_in_Idx,

    red_vh_θh_IDX
    ) =
        red_types_Idxs_etc   
    
    #------------------------------------------        

   (;full_gens_vh_Idxs,
    full_gens_θh_Idxs,

    full_slack_gens_vh_Idxs,
    full_non_slack_gens_vh_Idxs,
    full_non_gens_nodes_vh_Idxs,

    full_slack_gens_θh_Idxs,
    full_non_slack_gens_θh_Idxs,
    full_non_gens_nodes_θh_Idxs,

    full_slack_gens_vh_idx2Idx,
    full_non_slack_gens_vh_idx2Idx,
    full_non_gens_vh_idx2Idx,

    full_slack_gens_θh_idx2Idx,
    full_non_slack_gens_θh_idx2Idx,
    full_non_gens_θh_idx2Idx,

    full_dict_vh_idx2Idx,
    full_dict_θh_idx2Idx,

    full_slack_gens_vh_idx2Idx_in_Idx,
    full_non_slack_gens_vh_idx2Idx_in_Idx,
    full_non_gens_vh_idx2Idx_in_Idx,

    full_slack_gens_θh_idx2Idx_in_Idx,
    full_non_slack_gens_θh_idx2Idx_in_Idx,
    full_non_gens_θh_idx2Idx_in_Idx,

    full_dict_vh_idx2Idx_in_Idx,
    full_dict_θh_idx2Idx_in_Idx,

    full_vh_θh_IDX
    ) =
        full_types_Idxs_etc
    
    #------------------------------------------        

   (;intg_gens_vh_Idxs,
    intg_gens_θh_Idxs,

    intg_slack_gens_vh_Idxs,
    intg_non_slack_gens_vh_Idxs,
    intg_non_gens_nodes_vh_Idxs,

    intg_slack_gens_θh_Idxs,
    intg_non_slack_gens_θh_Idxs,
    intg_non_gens_nodes_θh_Idxs,

    intg_gen_id_Idxs,
    intg_gen_iq_Idxs,

    intg_slack_gens_vh_idx2Idx,
    intg_non_slack_gens_vh_idx2Idx,
    intg_non_gens_vh_idx2Idx,

    intg_slack_gens_θh_idx2Idx,
    intg_non_slack_gens_θh_idx2Idx,
    intg_non_gens_θh_idx2Idx,

    intg_dict_vh_idx2Idx,
    intg_dict_θh_idx2Idx,

    intg_slack_gens_vh_idx2Idx_in_Idx,
    intg_non_slack_gens_vh_idx2Idx_in_Idx,
    intg_non_gens_vh_idx2Idx_in_Idx,

    intg_slack_gens_θh_idx2Idx_in_Idx,
    intg_non_slack_gens_θh_idx2Idx_in_Idx,
    intg_non_gens_θh_idx2Idx_in_Idx,

    intg_dict_vh_idx2Idx_in_Idx,
    intg_dict_θh_idx2Idx_in_Idx,

    intg_vh_θh_id_iq_IDX
    ) =
        intg_types_Idxs_etc
    
    #------------------------------------------        

   (;
    semi_non_slack_gens_vh_Idxs,
    semi_non_gens_nodes_vh_Idxs,

    semi_non_slack_gens_θh_Idxs,
    semi_non_gens_nodes_θh_Idxs,

    semi_gen_id_Idxs,
    semi_gen_iq_Idxs,

    semi_non_slack_gens_vh_idx2Idx,
    semi_non_gens_vh_idx2Idx,

    semi_non_slack_gens_θh_idx2Idx,
    semi_non_gens_θh_idx2Idx,

    semi_dict_vh_idx2Idx,
    semi_dict_θh_idx2Idx,

    semi_non_slack_gens_vh_idx2Idx_in_Idx,
    semi_non_gens_vh_idx2Idx_in_Idx,

    semi_non_slack_gens_θh_idx2Idx_in_Idx,
    semi_non_gens_θh_idx2Idx_in_Idx,

    semi_dict_vh_idx2Idx_in_Idx,
    semi_dict_θh_idx2Idx_in_Idx,

    semi_vh_θh_id_iq_IDX
    ) =
        semi_types_Idxs_etc
    
    #------------------------------------------        

    dim_P_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :Q ] ) )
    
    dim_P_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
            ; some_param = [ :Q ] ) )

    
    dim_δ_ed_eq_pf = length(
        [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
            netd.nodes)...;] )
    
    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes;
            some_param = [ :loc_P ] ) == 0 ? 0 :
                length(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_P ] ) )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
                length(
                    get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_Q  ] )  )
    
    #----------------------------------------
    #----------------------------------------

    dim_pf_sta_PQ_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_P_g_loc_load,
          dim_Q_g_loc_load ]

    _,_, pf_sta_PQ_para_IDX =
        create_size_offset_Idx( dim_pf_sta_PQ_para )

    P_gens_sta_para_Idxs =
         pf_sta_PQ_para_IDX[1]

    Q_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[2]

    P_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[3]

    Q_non_gens_sta_para_Idxs =
        pf_sta_PQ_para_IDX[4]

    P_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[5]

    Q_g_loc_load_sta_para_Idxs =
        pf_sta_PQ_para_IDX[6]

    PQ_sta_para_Idxs = (
        ; P_gens_sta_para_Idxs,
        Q_gens_sta_para_Idxs,
        P_non_gens_sta_para_Idxs,
        Q_non_gens_sta_para_Idxs,
        P_g_loc_load_sta_para_Idxs,
        Q_g_loc_load_sta_para_Idxs )

    #----------------------------------------

    dim_pf_dyn_PQ_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_δ_ed_eq_pf,
          dim_P_g_loc_load,
          dim_Q_g_loc_load ]

    _,_, pf_dyn_PQ_para_IDX =
        create_size_offset_Idx(
            dim_pf_dyn_PQ_para )

    P_gens_dyn_para_Idxs =
         pf_dyn_PQ_para_IDX[1]

    Q_gens_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[2]

    P_non_gens_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[3]

    Q_non_gens_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[4]

    δ_ed_eq_pf_para_Idxs =
        pf_dyn_PQ_para_IDX[5]
    
    P_g_loc_load_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[6]

    Q_g_loc_load_dyn_para_Idxs =
        pf_dyn_PQ_para_IDX[7]

    PQ_dyn_para_Idxs = (
        ; P_gens_dyn_para_Idxs,
        Q_gens_dyn_para_Idxs,
        P_non_gens_dyn_para_Idxs,
        Q_non_gens_dyn_para_Idxs,
        δ_ed_eq_pf_para_Idxs,
        P_g_loc_load_dyn_para_Idxs,
        Q_g_loc_load_dyn_para_Idxs,
        pf_dyn_PQ_para_IDX )

    #----------------------------------------
    #----------------------------------------

    intg_vh_θh_id_iq_size =
        last(intg_vh_θh_id_iq_IDX[end])

    full_vh_θh_size =
        last(full_vh_θh_IDX[end])
    
    pf_dyn_PQ_para_size =
        sum( dim_pf_dyn_PQ_para )
        
    
    #----------------------------------------
    #----------------------------------------

    dims_intg_vh_θh_id_iq_and_pf_dyn_PQ_para =
        [intg_vh_θh_id_iq_size,
         pf_dyn_PQ_para_size]
    
    _,_, intg_vh_θh_id_iq_and_pf_dyn_PQ_para_IDX =
        create_size_offset_Idx(
            dims_intg_vh_θh_id_iq_and_pf_dyn_PQ_para )

    intg_vh_θh_id_iq_Idxs =
        intg_vh_θh_id_iq_and_pf_dyn_PQ_para_IDX[1]


    intg_pf_dyn_PQ_para_Idxs =
        intg_vh_θh_id_iq_and_pf_dyn_PQ_para_IDX[2]

    #----------------------------------------

    dims_full_vh_θh_and_pf_dyn_PQ_para =
        [full_vh_θh_size,
         pf_dyn_PQ_para_size]
    
    _,_, full_vh_θh_and_pf_dyn_PQ_para_IDX =
        create_size_offset_Idx(
            dims_full_vh_θh_and_pf_dyn_PQ_para )

    full_vh_θh_Idxs =
        full_vh_θh_and_pf_dyn_PQ_para_IDX[1]


    full_pf_dyn_PQ_para_Idxs =
        full_vh_θh_and_pf_dyn_PQ_para_IDX[2]
    
    #----------------------------------------
    #----------------------------------------

    return (;
            PQ_sta_para_Idxs,
            PQ_dyn_para_Idxs,
            
            nodes_types_idxs,
            n2s_idxs,
            
            u_ur_ui_Idx_in_state,            
            ur_ui_vh_θh_idxs,
            ur_ui_vh_θh_Idxs,

            red_types_Idxs_etc,
            full_types_Idxs_etc,
            intg_types_Idxs_etc,
            semi_types_Idxs_etc,            
            
            intg_vh_θh_id_iq_Idxs,
            intg_pf_dyn_PQ_para_Idxs,

            full_vh_θh_Idxs,
            full_pf_dyn_PQ_para_Idxs )
    
end


#------------------------------------------

function get_a_model_integrated_pf_para_idx_and_Idx(
    netd;
    Idxs_type =
        Idxs_type,
    no_control_device =
        false)
    
    #------------------------------------------            
    _, _, edges_orientation, nodes_node_idx_and_incident_edges_other_node_idx = get_Ynet_and_nodes_node_idx_and_incident_edges_other_node_idx( netd )
    
    #------------------------------------------    

    (; slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     load_trans_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx) =
         get_net_nodes_type_idxs( netd)
    
    # ----------------------------------------------------
    
    vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
                netd.nodes) )

    
    vec_Idx_gens_nodes_ω_ed_dash_eq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_ω_ed_dash_eq_dash_Idxs_in_plants(
                netd.nodes) )
        
    # ----------------------------------------------------

    vec_Idx_ra_Xd_Xq_Xd_dash_Xq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param_dims(
                netd.nodes
                ;some_param =
                    [:ra, :X_d, :X_q,
                     :X_d_dash, :X_q_dash]  ))

    vec_Idx_ra_Xd_dash_Xq_dash =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param_dims(
                netd.nodes
                ;some_param =
                    [ :ra, :X_d_dash, :X_q_dash ] ))


    vec_Idx_ra_Xd_Xq =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param_dims(
                netd.nodes
                ;some_param =
                    [ :ra, :X_d, :X_q ] ))
    
    # ------------------------------------------------
    
    vec_Idx_P_Q_nodes =
        get_flattened_to_components_vector_var_Idx(
            get_nodes_some_param_dims(
                netd.nodes
                ; some_param = [:P, :Q] ) )

    
    vec_Idx_P_Q_gens =
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P, :Q ] ) )

    
    vec_Idx_P_Q_non_gens =
        get_flattened_to_components_vector_var_Idx(
            get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P, :Q ] ) )
    
    #------------------------------------------
    
    gens_nodes_with_loc_loads_some_param_dims =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param =
                [ :loc_P, :loc_Q  ] )
    
    vec_Idx_P_Q_gens_loc_load =
        gens_nodes_with_loc_loads_some_param_dims == 0 ?
        nothing :
        get_flattened_to_components_vector_var_Idx(
            get_gens_nodes_with_loc_loads_some_param_dims(
                netd.nodes
                ; some_param = [ :loc_P, :loc_Q  ] ) )

    #------------------------------------------

    dim_P_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_gens =
        length( get_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :Q ] ) )
    
    dim_P_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :P] ) )
    
    dim_Q_non_gens =
        length( get_non_gens_nodes_some_param(
                netd.nodes
                ; some_param = [ :Q ] ) )       
    
    dim_δ_ed_eq_pf = length(
        [get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
            netd.nodes)...;] )

    dim_P_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes;
            some_param = [ :loc_P ] ) == 0 ? 0 :
                length(get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_P ] ) )
    
    dim_Q_g_loc_load =
        get_gens_nodes_with_loc_loads_some_param_dims(
            netd.nodes
            ; some_param = [ :loc_Q  ] ) == 0 ? 0 : 
                length( get_gens_nodes_with_loc_loads_some_param_dims(
                    netd.nodes; some_param = [ :loc_Q  ] )  )
    
    #----------------------------------------

    dim_no_loc_load_with_δ_ed_NL_pf_para  =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_δ_ed_eq_pf ]

    _,_, no_loc_load_with_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_no_loc_load_with_δ_ed_NL_pf_para  )

    no_loc_load_with_δ_ed_P_gens_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[1]

    no_loc_load_with_δ_ed_Q_gens_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[2]

    no_loc_load_with_δ_ed_P_non_gens_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[3]

    no_loc_load_with_δ_ed_Q_non_gens_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[4]

    no_loc_load_with_δ_ed_δ_ed_eq_Idxs =
        no_loc_load_with_δ_ed_CA_pf_para_Idxs[5]

    pf_no_loc_load_with_δ_ed_Idxs = (
        ; no_loc_load_with_δ_ed_P_gens_Idxs,
        no_loc_load_with_δ_ed_Q_gens_Idxs,
        no_loc_load_with_δ_ed_P_non_gens_Idxs,
        no_loc_load_with_δ_ed_Q_non_gens_Idxs,
        no_loc_load_with_δ_ed_δ_ed_eq_Idxs )

    #----------------------------------------
    
    dim_with_loc_load_with_δ_ed_NL_pf_para =
        [ dim_P_gens,
          dim_Q_gens,
          dim_P_non_gens,
          dim_Q_non_gens,
          dim_δ_ed_eq_pf,
          dim_P_g_loc_load,
          dim_Q_g_loc_load  ]

    _,_, with_loc_load_with_δ_ed_CA_pf_para_Idxs =
        create_size_offset_Idx(
            dim_with_loc_load_with_δ_ed_NL_pf_para )

    with_loc_load_with_δ_ed_P_gens_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[1]

    with_loc_load_with_δ_ed_Q_gens_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[2]

    with_loc_load_with_δ_ed_P_non_gens_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[3]

    with_loc_load_with_δ_ed_Q_non_gens_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[4]

    with_loc_load_δ_ed_eq_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[5]

    with_loc_load_with_δ_ed_P_g_loc_load_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[6]

    with_loc_load_with_δ_ed_Q_g_loc_load_Idxs =
        with_loc_load_with_δ_ed_CA_pf_para_Idxs[7]

    pf_with_loc_load_with_δ_ed_Idxs = (
        ; with_loc_load_with_δ_ed_P_gens_Idxs,
        with_loc_load_with_δ_ed_Q_gens_Idxs,
        with_loc_load_with_δ_ed_P_non_gens_Idxs,
        with_loc_load_with_δ_ed_Q_non_gens_Idxs,
        with_loc_load_δ_ed_eq_Idxs,
        with_loc_load_with_δ_ed_P_g_loc_load_Idxs,
        with_loc_load_with_δ_ed_Q_g_loc_load_Idxs)

    #----------------------------------------

    dims_δ_ω_ed_dash_eq = length.(
        get_gens_nodes_δ_ω_ed_dash_eq_dash_Idxs_in_plants(
            netd.nodes))

    _,_, δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq =
        create_size_offset_Idx( dims_δ_ω_ed_dash_eq )

    #------------------------------------------
    
    dicts_net_to_streamlined_idx =
        get_dict_net_to_streamlined_idx(
            netd)
    
    #------------------------------------------

    pf_para_Idxs = (
        ; pf_no_loc_load_with_δ_ed_Idxs,
        pf_with_loc_load_with_δ_ed_Idxs,
        δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq )
    
    #------------------------------------------

    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx
    
    net_comp_type_idx =
        (; slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         load_nodes_idx,
         transmission_nodes_idx,
         non_gens_nodes_idx,
         load_trans_nodes_idx,
         gens_with_loc_load_idx,
         all_nodes_idx )
    
    #------------------------------------------
    
    vec_Idx =
        (; vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash,
         vec_Idx_gens_nodes_ω_ed_dash_eq_dash,
         vec_Idx_ra_Xd_Xq_Xd_dash_Xq_dash,
         vec_Idx_ra_Xd_dash_Xq_dash,
         vec_Idx_ra_Xd_Xq,
         vec_Idx_P_Q_nodes,
         vec_Idx_P_Q_gens,
         vec_Idx_P_Q_non_gens,
         vec_Idx_P_Q_gens_loc_load,
         loc_load_exist )
    
    #------------------------------------------
    
    return (;loc_load_exist,
            dicts_net_to_streamlined_idx,
            pf_para_Idxs,
            net_comp_type_idx,
            vec_Idx )
    
end


#-----------------------------------------------------
#-----------------------------------------------------

function get_a_model_integrated_pf_PQ_param(
    netd; P_Q_only = false )

    #------------------------------------------
                
    loc_load_exist =
        loc_load_exist_bool(netd)

    #------------------------------------------

    # P_gens = get_gens_nodes_some_param(
    #     netd.nodes; some_param = [ :P] )
    
    # Q_gens = get_gens_nodes_some_param(
    #     netd.nodes;some_param = [ :Q ] )    
    
    # P_non_gens = convert.( Float64, get_non_gens_nodes_some_param(
    #             netd.nodes ; some_param = [ :P] ) )
    
    # Q_non_gens =  convert.( Float64, get_non_gens_nodes_some_param(
    #             netd.nodes ; some_param = [ :Q ] ) )

    # P_g_loc_load = get_gens_nodes_with_loc_loads_some_param(
    #     netd.nodes; some_param = [ :loc_P ] )
    
    # Q_g_loc_load = get_gens_nodes_with_loc_loads_some_param(
    #     netd.nodes; some_param = [ :loc_Q  ] )

    #------------------------------------------
    #------------------------------------------

    P_gens = get_gens_nodes_some_param(
        netd.nodes; some_param = [ :P] )
    P_gens = [P_gens...;]
    
    Q_gens = get_gens_nodes_some_param(
        netd.nodes;some_param = [ :Q ])    
    Q_gens = [Q_gens...;]
    
    P_non_gens = get_non_gens_nodes_some_param(
        netd.nodes ; some_param = [ :P] )
    P_non_gens =  convert.( Float64, [P_non_gens...;])
   
    
    Q_non_gens = get_non_gens_nodes_some_param(
        netd.nodes ; some_param = [ :Q ] )
    Q_non_gens = convert.(Float64, [Q_non_gens...;]) 

    # 

    P_g_loc_load = get_gens_nodes_with_loc_loads_some_param(
        netd.nodes; some_param = [:loc_P ])
    P_g_loc_load =  P_g_loc_load == [] ? [] :
        [P_g_loc_load...;]
    
    Q_g_loc_load = get_gens_nodes_with_loc_loads_some_param(
        netd.nodes; some_param = [:loc_Q ])
    Q_g_loc_load = Q_g_loc_load == [] ? [] :
        [Q_g_loc_load...;]

    #----------------------------------------    
    #----------------------------------------

    gens_nodes_ra_Xd_Xq_Xd_dash_Xq_dash =
        get_gens_nodes_some_param(
            netd.nodes; some_param =
                [:ra, :X_d, :X_q, :X_d_dash, :X_q_dash ])
    
    gens_nodes_ra_Xd_Xq = 
        get_gens_nodes_some_param(
            netd.nodes; some_param =
                [ :ra, :X_d, :X_q  ] )
    
    gens_nodes_ra_Xd_dash_Xq_dash =
        get_gens_nodes_some_param(
            netd.nodes;
            some_param =
                [ :ra, :X_d_dash, :X_q_dash ])

    #----------------------------------------

    if  P_Q_only == false

        return (; P_gens,
                Q_gens,
                P_non_gens,
                Q_non_gens,
                P_g_loc_load,
                Q_g_loc_load,
                loc_load_exist,
                gens_nodes_ra_Xd_dash_Xq_dash )
    else

        return (; P_gens,
                Q_gens,
                P_non_gens,
                Q_non_gens,
                P_g_loc_load,
                Q_g_loc_load )        
    end
    
    
end



function get_a_model_integrated_sta_pf_PQ_param(netd)

    #------------------------------------------
                
    loc_load_exist =
        loc_load_exist_bool(netd)

    #------------------------------------------

    P_gens = get_gens_nodes_some_param(
        netd.nodes; some_param = [ :P] )
    P_gens = [P_gens...;]
    
    Q_gens = get_gens_nodes_some_param(
        netd.nodes;some_param = [ :Q ] )    
    Q_gens = [Q_gens...;]
    
    P_non_gens = get_non_gens_nodes_some_param(
        netd.nodes ; some_param = [ :P] )
    P_non_gens =  convert.( Float64, [P_non_gens...;]  )
   
    
    Q_non_gens = get_non_gens_nodes_some_param(
        netd.nodes ; some_param = [ :Q ] )
    Q_non_gens = convert.(Float64, [Q_non_gens...;]  ) 

    # 

    P_g_loc_load = get_gens_nodes_with_loc_loads_some_param(
        netd.nodes; some_param = [ :loc_P ] )
    P_g_loc_load =  P_g_loc_load == [] ? [] : [P_g_loc_load...;]
    
    Q_g_loc_load = get_gens_nodes_with_loc_loads_some_param(
        netd.nodes; some_param = [ :loc_Q  ] )
    Q_g_loc_load = Q_g_loc_load == [] ? [] : [Q_g_loc_load...;]

    #----------------------------------------
    #----------------------------------------

    return (; P_gens,
            Q_gens,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,
            loc_load_exist )
    
end


function get_a_model_integrated_dyn_vars_idx_and_Idx(
    netd;
    Idxs_type =
        :Idxs_hybrid,
    no_control_device =
        false,
    only_gen  =
        false )


    if Idxs_type == :Idxs_hybrid

        # ----------------------------------------------    
        
    elseif Idxs_type == :Idxs_im

        # ----------------------------------------------

        dict_sys_to_im = get_net_to_im_indices_dict(
            netd  )

        (;pure_states_Idx_in_system,
         im_algebraic_vars_Idx_in_system,
         ur_ui_Idx_in_system,
         im_Idx,
         pure_states_Idx,
         im_algebraic_vars_Idx,
         nodes_ur_ui_Idx,
         pure_states_im_algebraic_vars_ur_ui_Idx_in_system,
         im_state_Idx,
         net_to_im_idx_conversion_dict) =
            get_im_indices_and_conversion_dict(
                netd  )

        # ----------------------------------------------
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
            get_im_δ_ω_ed_dash_eq_dash_Idx(netd )

        # ----------------------------------------------

        im_model_each_gen_nodes_pure_states_idx_in_state =
            get_gens_im_pure_states_indices_in_state(netd )

        # ----------------------------------------------

        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
            get_im_δ_ω_ed_dash_eq_dash_Idx(netd )
        
        # ----------------------------------------------

        (; Ax_Bx_Cx_views,
         Ax_Bx_Cx_matrix,
         states_and_mat_Idxs) =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(
            gens_nodes_collection )

        vec_Ax_views, vec_Bx_views, vec_Cx_views  =
            Ax_Bx_Cx_views

        Ax_matrix, Bx_matrix, Cx_matrix =
            Ax_Bx_Cx_matrix


        (; nodes_state_Idx,
         Bx_idxs, Cx_idxs,
         id_iq_ph_vh_idxs,
         ω_ref_ωref_v_ref_idxs)  =
            states_and_mat_Idxs

        init_im_plants_system_matrices_views!(
            (vec_Ax_views, vec_Bx_views, vec_Cx_views),
        gens_nodes_collection )

    sys_states_Idxs_and_mat_Idxs =
        (; nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         id_iq_ph_vh_idxs,
         ω_ref_ωref_v_ref_idxs )
        
        # ----------------------------------------------
        
        (;
         pure_states_Idx_in_system,
        ur_ui_Idx_in_system,
        im_Idx,
        pure_states_Idx,
        nodes_ur_ui_Idx,
        pure_states_im_algebraic_vars_ur_ui_Idx_in_system,
        im_state_Idx,
        net_to_im_idx_conversion_dict,
        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_ω_ed_dash_eq_dash_Idxs_in_state,
        nodes_δ_ed_dash_eq_dash_Idxs_in_state,
         im_model_each_gen_nodes_pure_states_idx_in_state
         )        

        #---------------------------------------------------    
        
    elseif Idxs_type == :Idxs_industrial

        # -------------------------------------------------   

        dict_sys_to_industry =
            get_net_to_industrial_model_indices_dict(
                netd;
                no_control_device =
                    only_gen  )

        (;pure_states_Idx_in_system,
         ur_ui_Idx_in_system, industrial_Idx,
         industrial_model_pure_states_Idx,
         industrial_model_ur_ui_Idx,
         pure_states_and_ur_ui_Idx_in_system,
         industrial_model_state_Idx,
         net_to_industrial_idx_conversion_dict) =
            get_industrial_model_indices_and_conversion_dict(
                netd; no_control_device =
                    only_gen  )

        # -----------------------------------------------

        nodes_δ_ω_ed_dash_eq_dash_Idxs_in_state =
            get_industrial_δ_ω_ed_dash_eq_dash_Idx(
                netd; no_control_device = only_gen )

        nodes_ω_ed_dash_eq_dash_Idxs_in_state   =
            get_industrial_ω_ed_dash_eq_dash_Idx(
                netd; no_control_device = only_gen )

        nodes_δ_ed_dash_eq_dash_Idxs_in_state   =
            get_industrial_δ_ed_dash_eq_dash_Idx(
                netd; no_control_device = only_gen )

        #-----------------------------------------------

        industrial_model_each_gen_nodes_pure_states_idx_in_state =
            get_industrial_gens_pure_states_indices_in_state(
                netd;
                no_control_device =
                    only_gen )

        industrial_model_each_gen_nodes_stab_states_idx_in_state =
            get_industrial_gens_stab_states_indices_in_state(
                netd;
                no_control_device =
                    only_gen )

        #-----------------------------------------------  

        nodes_u_Idx_in_ranges =
            get_nodes_u_Idx_in_ranges(
                nodes_u_Idx )    
        #-----------------------------------------------

        non_gens_idx =
            get_load_trans_nodes_Idx(
                netd.nodes )

        #------------------------------------------------

        industrial_model_nodes_voltage_Idx =
            (; gens_idx,
             non_gens_idx,
             nodes_u_Idx,
             nodes_u_Idx_in_ranges )

        industrial_model_misc_Idx =
            (;
             nodes_u_Idx,
             nodes_u_Idx_in_ranges,
             industrial_model_pure_states_Idx,
             industrial_model_each_gen_nodes_pure_states_idx_in_state, gens_nodes_collection )

        
        #----------------------------------------------
        
        
    else
        nothing
    end
        
end



function get_a_model_integrated_dyn_param(
    netd;
    Idxs_type =
        :Idxs_hybrid,
    no_control_device =
        false,
    pf_alg =
        pf_alg  )
    
    #--------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )

    #----------------------------------------

    if Idxs_type ==  :Idxs_hybrid

        
    elseif Idxs_type ==  :Idxs_im

    #----------------------------------------

        
    #----------------------------------------

        get_a_model_integrated_pf_sta_powerflow(
    netd;
    pf_alg =
        pf_alg )
        
        
    #----------------------------------------

        named_tup_pf_result = power_balance_powerflow(
            x0_vh_θh,
            mismatch,
            sd_pf_views,
            (nodes_name, branches_name) ,
            pf_net_param;
            maxiter=maxiter,
            ftol=ftol,
            xtol=xtol,
            with_δ_ed_eq = with_δ_ed_eq )


        bus_dict_init = named_tup_pf_result.bus_dict_init

        branch_dict_init = named_tup_pf_result.branch_dict_init

    #----------------------------------------

        state .= im_model_init_operationpoint(
            netd, bus_dict_init  )
        
    #----------------------------------------



        gens_vh_θh =
            get_gens_vh_θh( nodes_pf_U_view, gens_idx )
        

        gens_vh_θh_view =
            @view gens_vh_θh[ 1:length(gens_vh_θh ) ]


        gens_nodes_ωs_ωref0_vref0_porder0 =
            get_gens_nodes_ωs_ωref0_vref0_porder0(
            gens_nodes_collection,
            bus_dict_init )


        gens_nodes_ωs_ωref0_vref0_porder0_view = view(
            gens_nodes_ωs_ωref0_vref0_porder0,
            1:length( gens_nodes_ωs_ωref0_vref0_porder0 ) )
        

        gens_nodes_τm_vf =
            get_gens_τm_vf( gens_nodes_collection,
                            bus_dict_init )


        vec_τm_vf_views =
            view( gens_nodes_τm_vf,
                  1:length(gens_nodes_τm_vf) )

        gens_dynamic_id_iq_pg_vh_by_vhθh =
            get_gens_dynamic_id_iq_pg_vh_by_vhθh(
            gens_vh_θh_view,
            gen_nodes_δ_ω_ed_dash_eq_dash_views,
            gen_nodes_ra_Xd_dash_Xq_dash_view )


        im_ωs_ωref_vref_vhθh_idq_etc =  (
            ; gens_dynamic_id_iq_pg_vh_by_vhθh,
            gens_nodes_ωs_ωref0_vref0_porder0,
            gens_nodes_τm_vf,
            gens_vh_θh )
        
    #----------------------------------------
    #----------------------------------------
        
        Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, _ =
            create_gens_nodes_im_aggregate_system_matrices_Idx_and_views( gens_nodes_collection )

        vec_Ax_views, vec_Bx_views, vec_Cx_views  =
            Ax_Bx_Cx_views

        Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix

        init_im_plants_system_matrices_views!(
            (vec_Ax_views,
             vec_Bx_views,
             vec_Cx_views),
            gens_nodes_collection )

        
    elseif Idxs_type ==  :Idxs_industrial

        
    else
        
        nothing
        
    end

    #----------------------------------------    
    #----------------------------------------

    return (; )
    
end


function get_a_model_integrated_syms_mass_matrix_(netd;
    Idxs_type =
        :Idxs_hybrid,
    no_control_device =
        false)
    
    #--------------------------------------

    non_gens_nodes =
        get_non_gens_nodes( netd.nodes )

    gens_nodes_collection =
        get_gens_nodes( netd.nodes )


    #----------------------------------------
    
    (; network_bus_names,
       non_gens_bus_names,
       gens_bus_names,
       Loads_bus_names,
       Trans_bus_names) =
           make_case_buses_names(
               netd.nodes )
    
    net_class_names =
        (;network_bus_names,
         non_gens_bus_names,
         gens_bus_names,
         Loads_bus_names,
         Trans_bus_names )
       
    #----------------------------------------

    if Idxs_type ==  :Idxs_hybrid
        
    elseif Idxs_type ==  :Idxs_im


    (; net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_im_algebraic_vars_labels,
     gens_nodes_im_vars_labels) =
         generate_im_model_labels(
             ; nodes =  netd.nodes )
        

    im_net_states_and_var_labels =
        (; net_bus_volts_labels,
         gens_nodes_pure_states_labels,
         gens_nodes_im_algebraic_vars_labels,
         gens_nodes_im_vars_labels )


    (; im_model_ode_mass_matrix,
     im_model_pf_mass_matrix) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes = netd.nodes )

    (; gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes = netd.nodes  )

    im_sym =
        gens_nodes_im_vars_labels

    im_mass_matrix =
        im_model_ode_mass_matrix
        
    #-------------------------------------------

    para_net_names_labels_syms =
        (; net_class_names,
         im_net_states_and_var_labels,
         im_sym )
    




        
    elseif Idxs_type ==  :Idxs_industrial
        
    else
        
        nothing
        
    end
    
    return (; )
    
end



function get_a_model_integrated_pf_and_dyn_idx_and_Idx(netd;
    Idxs_type =
        :Idxs_hybrid,
    no_control_device =
        false)
    
end



#-------------------------------
# powerflow addendum
#-------------------------------


function auto_diff_get_closure_pf_param(netd)


    #-----------------------------------------------

    pf_net_param =
        get_powerflow_net_parameters( netd )
    
    pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx


    # ----------------------------------------------------

    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist

    # ----------------------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx

    gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    gens_with_loc_load_idx =
        net_comp_type_idx.gens_with_loc_load_idx

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx
    
    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash
    
    #----------------------------------------------------

    syms =
        get_network_vars_labels( netd )
    
    #----------------------------------------------------

    state =
        zeros(length( get_network_vars_labels( netd ) ))

    #----------------------------------------------------
    
    models_networks_labels =
        get_models_networks_labels(netd )
    
    #----------------------------------------------------

    models_gens_nodes_some_vars_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )       

    nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid =
        models_gens_nodes_some_vars_Idxs.Idxs_hybrid

    nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial =
        models_gens_nodes_some_vars_Idxs.Idxs_industrial
    
    nodes_δ_ω_ed_dash_eq_dash_Idxs_im =
        models_gens_nodes_some_vars_Idxs.Idxs_im
    
    #----------------------------------------------------

    δ_ω_ed_dash_eq_dash_view =
        get_nodes_δ_ω_ed_dash_eq_dash_view(
            state, netd )
    
    #----------------------------------------------------    
    
    # ---------------------------------------------------

    state_view =
        view(state, 1:length(state))

    # ----------------------------------------------------    
    # for pf flat start, ur = vh = 1,  ui = θh = 0        
    # ----------------------------------------------------    
    
    if init_pf == true 

        state_view[ur_idx] .= ones(  length( ur_IDX ))
        state_view[ui_idx] .= zeros( length( ui_IDX ))

    end
    
    # ----------------------------------------------------
    # x  = similar( state_x0 )
    #
    # Decoupling pf_state from state_x0
    # ----------------------------------------------------    

    pf_state = [ state; ]
        
    #----------------------------------------------------
    
    nodes_u_view  =
        [ view(state, nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ]

    # ----------------------------------------------------  
    
    nodes_pf_U_view  =
        [ view(pf_state , nodes_u_Idx[Ind])
          for Ind in collect(1:length(nodes_u_Idx)) ] 
    
    # ----------------------------------------------------
    
    uh_state =
        state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]
    
    x0_ur_ui =
        [state_view[ur_ui_idx][ ur_IDX ]...;
         state_view[ur_ui_idx][ ui_IDX ]...]
        
    x0_vh_θh =
        [abs.(uh_state)...;
         angle.(uh_state)...]

    working_vh_θh_view =
        view(x0_vh_θh, 1:length(x0_vh_θh))

    working_vh_θh_DiffCache = similar( x0_vh_θh )
    
    
    x0_vh_view =
        @view x0_vh_θh[vh_IDX]

    x0_θh_view  =
        @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view =
        @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0  =
        [ red_vh_θh_0_view; ]

    red_vh_θh_DiffCache =
        similar( red_vh_θh_0 )

    mismatch =
        similar( red_vh_θh_0 )
    
    # ----------------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))
    
    # ----------------------------------------------------
    
    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet =
        zeros(ComplexF64,
              length( uh_state ))
    
    Inet_view =
        view( Inet,
              1:length( Inet ) )

    Iinj =
        zeros(ComplexF64,
              length( uh_state ))

    Iinj_view =
        view(Iinj,
             1:length( Iinj ))


    idq_wt_pad =
        zeros(ComplexF64, length( uh_state ))

    idq_wt_pad_view =
        view(idq_wt_pad, 1:length( uh_state ) )
    
        
    #----------------------------------------------------

    global_pf_views =
        (;
          nodes_pf_U_view,
          Inet_view,
         Iinj_view,
         idq_wt_pad_view )

    sd_pf_views =
        (;
         nodes_pf_U_view,
         Inet_view,
         Iinj_view,
         idq_wt_pad_view, 
         δ_ω_ed_dash_eq_dash_view )

    # ---------------------------------------------------
    # global_pf param
    # ---------------------------------------------------
    
    global_pf_param =
        (; pf_net_param,
         sd_pf_views,
         mismatch )
    
    # # ---------------------------------------------------
    # # ---------------------------------------------------

    pf_NL_P_Q_para =
        pf_net_misc.pf_NL_P_Q_para
    
    P_gens_NL_para =
        pf_NL_P_Q_para.P_gens_NL_para
    
    Q_gens_NL_para =
        pf_NL_P_Q_para.Q_gens_NL_para
        
    P_load_NL_para =
        pf_NL_P_Q_para.P_load_NL_para
    
    Q_load_NL_para =
        pf_NL_P_Q_para.Q_load_NL_para
    
    P_g_loc_load_NL_para =
        pf_NL_P_Q_para.P_g_loc_load_NL_para
    
    Q_g_loc_load_NL_para =
        pf_NL_P_Q_para.Q_g_loc_load_NL_para
    
    # ---------------------------------------------------
    # ---------------------------------------------------
    
    no_loc_load_no_δ_ed_eq_pf_param = ComponentArray(
        P_gens = pf_NL_P_Q_para.P_gens_NL_para,
        Q_gens = pf_NL_P_Q_para.Q_gens_NL_para,
        P_load = pf_NL_P_Q_para.P_load_NL_para,
        Q_load = pf_NL_P_Q_para.Q_load_NL_para )
    
    # ---------------------------------------------------

    no_loc_load_with_δ_ed_eq_pf_param = ComponentArray(
        P_gens = pf_NL_P_Q_para.P_gens_NL_para,
        Q_gens = pf_NL_P_Q_para.Q_gens_NL_para,
        P_load = pf_NL_P_Q_para.P_load_NL_para,
        Q_load = pf_NL_P_Q_para.Q_load_NL_para,       
        δ_ed_eq_pf = [δ_ω_ed_dash_eq_dash_view...;]  )

    # ---------------------------------------------------

    with_loc_load_no_δ_ed_eq_pf_param = ComponentArray(
        no_loc_load_no_δ_ed_eq_pf_param;
        P_g_loc_load = pf_NL_P_Q_para.P_g_loc_load_NL_para,
        Q_g_loc_load = pf_NL_P_Q_para.Q_g_loc_load_NL_para )

    with_loc_load_with_δ_ed_eq_pf_param = ComponentArray(
        no_loc_load_with_δ_ed_eq_pf_param;
        P_g_loc_load = pf_NL_P_Q_para.P_g_loc_load_NL_para,
        Q_g_loc_load = pf_NL_P_Q_para.Q_g_loc_load_NL_para )

    # ---------------------------------------------------
    # ---------------------------------------------------

    # closure_pf_param 
    
    return  (; pf_net_param, net_comp_type_idx, gens_nodes_idx, non_gens_nodes_idx, gens_with_loc_load_idx, vec_Idx, vec_Idx_δ_ω_ed_dash_eq_dash, syms, state, models_networks_labels , models_gens_nodes_some_vars_Idxs, δ_ω_ed_dash_eq_dash_view, state_view, pf_state, nodes_u_view, nodes_pf_U_view, uh_state, x0_ur_ui, x0_vh_θh, working_vh_θh_view, working_vh_θh_DiffCache, x0_vh_view, x0_θh_view, red_vh_θh_0_view, red_vh_θh_0, red_vh_θh_DiffCache,  mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, no_loc_load_no_δ_ed_eq_pf_param, no_loc_load_with_δ_ed_eq_pf_param, with_loc_load_no_δ_ed_eq_pf_param, with_loc_load_with_δ_ed_eq_pf_param )
    
end



function auto_diff_get_streamlined_closure_pf_param(netd; init_pf = true, with_δ_ed_eq = false )


    #-----------------------------------------------

    pf_net_param =
        get_streamlined_powerflow_net_parameters( netd )
    
    pf_net, pf_idx_and_state, pf_param_views, pf_limits, pf_Idx, pf_and_dyn_idx_and_Idx , pf_net_misc = pf_net_param

    Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, slack_ur_ui_Idx_in_state, non_slack_ur_ui_Idx_in_state, ur_ui_Idx_in_state = pf_idx_and_state

    ra_Xd_Xq_view, ra_Xd_dash_Xq_dash_view, ra_Xd_Xq_Xd_dash_Xq_dash_view, P_Q_nodes_view, P_Q_gens_view, P_Q_gens_loc_load_view, P_Q_non_gens_view = pf_param_views

    load_trans_nodes_Idx_and_vlimits = pf_limits

    slack_bus_idx, nodes_u_Idx, gens_idx, ur_IDX, ui_IDX, vh_IDX, θh_IDX, red_vh_θh_idx, ur_idx, ui_idx, ur_ui_idx = pf_Idx

    # ----------------------------------------------------

    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist

    # ----------------------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx

    gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    gens_with_loc_load_idx =
        net_comp_type_idx.gens_with_loc_load_idx

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx
    
    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash
    
    #----------------------------------------------------

    syms =
        get_network_vars_labels( netd )
    
    #----------------------------------------------------

    state =
        zeros(length( get_network_vars_labels( netd ) ))

    #----------------------------------------------------
    
    models_networks_labels =
        get_models_networks_labels(netd )
    
    #----------------------------------------------------

    models_gens_nodes_some_vars_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )       

    nodes_δ_ω_ed_dash_eq_dash_Idxs_hybrid =
        models_gens_nodes_some_vars_Idxs.Idxs_hybrid

    nodes_δ_ω_ed_dash_eq_dash_Idxs_industrial =
        models_gens_nodes_some_vars_Idxs.Idxs_industrial
    
    nodes_δ_ω_ed_dash_eq_dash_Idxs_im =
        models_gens_nodes_some_vars_Idxs.Idxs_im
    
    #----------------------------------------------------

    δ_ω_ed_dash_eq_dash_view =
        get_nodes_δ_ω_ed_dash_eq_dash_view(
            state, netd )
    
    #----------------------------------------------------    
    
    # ---------------------------------------------------

    state_view =
        view(state, 1:length(state))

    # ----------------------------------------------------   
    # for pf flat start, ur = vh = 1,  ui = θh = 0        
    # ----------------------------------------------------   
    
    if init_pf == true 

        state_view[ur_idx] .= ones(  length( ur_IDX ))
        state_view[ui_idx] .= zeros( length( ui_IDX ))

    end
    
    # ----------------------------------------------------
    # x  = similar( state_x0 )
    #
    # Decoupling pf_state from state_x0
    # ----------------------------------------------------    

    pf_state = [ state; ]
        
    #----------------------------------------------------
    
    nodes_u_view  =
        [ view(state, nodes_u_Idx[Ind])
          for Ind in
              collect(1:length(nodes_u_Idx)) ]

    # ----------------------------------------------------  
    
    nodes_pf_U_view  =
        [ view(pf_state , nodes_u_Idx[Ind])
          for Ind in collect(1:length(nodes_u_Idx)) ] 
    
    # ----------------------------------------------------
    
    uh_state =
        state_view[ur_ui_idx][ ur_IDX ] +
        im * state_view[ur_ui_idx][ ui_IDX ]
    
    x0_ur_ui =
        [state_view[ur_ui_idx][ ur_IDX ]...;
         state_view[ur_ui_idx][ ui_IDX ]...]
        
    x0_vh_θh =
        [abs.(uh_state)...;
         angle.(uh_state)...]

    working_vh_θh_view =
        view(x0_vh_θh, 1:length(x0_vh_θh))

    working_vh_θh_DiffCache = similar( x0_vh_θh )
    
    
    x0_vh_view =
        @view x0_vh_θh[vh_IDX]

    x0_θh_view  =
        @view x0_vh_θh[θh_IDX]

    red_vh_θh_0_view =
        @view x0_vh_θh[ red_vh_θh_idx  ]

    red_vh_θh_0  =
        [ red_vh_θh_0_view; ]

    red_vh_θh_DiffCache =
        similar( red_vh_θh_0 )

    mismatch =
        similar( red_vh_θh_0 )
    
    # ----------------------------------------------------

    Jac_vh_θh = zeros(length( red_vh_θh_idx ),
                      length( red_vh_θh_idx ))
    
    # ----------------------------------------------------
    
    # For storing current. The dims of
    # uh_state_x0 = x0_vh_θh = ir_ii

    Inet =
        zeros(ComplexF64,
              length( uh_state ))
    
    Inet_view =
        view( Inet,
              1:length( Inet ) )

    Iinj =
        zeros(ComplexF64,
              length( uh_state ))

    Iinj_view =
        view(Iinj,
             1:length( Iinj ))


    # idq_wt_pad =
    #     zeros(ComplexF64, length( uh_state ))

    # idq_wt_pad_view =
    #     view(idq_wt_pad, 1:length( uh_state ) )
    
    
    idq_wt_pad =
        zeros(ComplexF64, length( gens_nodes_idx ))
    

    idq_wt_pad_view =
        view(idq_wt_pad, 1:length( gens_nodes_idx ) )
        
    #----------------------------------------------------

    global_pf_views =
        (;
          nodes_pf_U_view,
          Inet_view,
         Iinj_view,
         idq_wt_pad_view )

    sd_pf_views =
        (;
         nodes_pf_U_view,
         Inet_view,
         Iinj_view,
         idq_wt_pad_view, 
         δ_ω_ed_dash_eq_dash_view )

    # ---------------------------------------------------
    # global_pf param
    # ---------------------------------------------------
    
    global_pf_param =
        (; pf_net_param,
         sd_pf_views,
         mismatch )

    # ---------------------------------------------------

    pf_NL_P_Q_para =
        pf_net_misc.pf_NL_P_Q_para
    
    P_gens_NL_para =
        pf_NL_P_Q_para.P_gens_NL_para
    
    Q_gens_NL_para =
        pf_NL_P_Q_para.Q_gens_NL_para
        
    P_load_NL_para =
        pf_NL_P_Q_para.P_load_NL_para
    
    Q_load_NL_para =
        pf_NL_P_Q_para.Q_load_NL_para
    
    P_g_loc_load_NL_para =
        pf_NL_P_Q_para.P_g_loc_load_NL_para
    
    Q_g_loc_load_NL_para =
        pf_NL_P_Q_para.Q_g_loc_load_NL_para
    
    # ---------------------------------------------------
    # ---------------------------------------------------
    
    no_loc_load_no_δ_ed_eq_pf_param = ComponentArray(
        P_gens = pf_NL_P_Q_para.P_gens_NL_para,
        Q_gens = pf_NL_P_Q_para.Q_gens_NL_para,
        P_load = pf_NL_P_Q_para.P_load_NL_para,
        Q_load = pf_NL_P_Q_para.Q_load_NL_para )
    
    # ---------------------------------------------------

    no_loc_load_with_δ_ed_eq_pf_param = ComponentArray(
        P_gens = pf_NL_P_Q_para.P_gens_NL_para,
        Q_gens = pf_NL_P_Q_para.Q_gens_NL_para,
        P_load = pf_NL_P_Q_para.P_load_NL_para,
        Q_load = pf_NL_P_Q_para.Q_load_NL_para,       
        δ_ed_eq_pf = [δ_ω_ed_dash_eq_dash_view...;]  )

    # ---------------------------------------------------

    with_loc_load_no_δ_ed_eq_pf_param = ComponentArray(
        no_loc_load_no_δ_ed_eq_pf_param;
        P_g_loc_load = pf_NL_P_Q_para.P_g_loc_load_NL_para,
        Q_g_loc_load = pf_NL_P_Q_para.Q_g_loc_load_NL_para )

    with_loc_load_with_δ_ed_eq_pf_param = ComponentArray(
        no_loc_load_with_δ_ed_eq_pf_param;
        P_g_loc_load = pf_NL_P_Q_para.P_g_loc_load_NL_para,
        Q_g_loc_load = pf_NL_P_Q_para.Q_g_loc_load_NL_para )

    # ---------------------------------------------------
    # ---------------------------------------------------
    
    # ---------------------------------------------------    
    # ---------------------------------------------------

    # closure_pf_param 
    
    return  (; pf_net_param, net_comp_type_idx, gens_nodes_idx, non_gens_nodes_idx, gens_with_loc_load_idx, vec_Idx, vec_Idx_δ_ω_ed_dash_eq_dash, syms, state, models_networks_labels , models_gens_nodes_some_vars_Idxs, δ_ω_ed_dash_eq_dash_view, state_view, pf_state, nodes_u_view, nodes_pf_U_view, uh_state, x0_ur_ui, x0_vh_θh, working_vh_θh_view, working_vh_θh_DiffCache, x0_vh_view, x0_θh_view, red_vh_θh_0_view, red_vh_θh_0, red_vh_θh_DiffCache,  mismatch, Jac_vh_θh, Inet, Inet_view, Iinj, Iinj_view, idq_wt_pad, idq_wt_pad_view, global_pf_views, sd_pf_views, global_pf_param, no_loc_load_no_δ_ed_eq_pf_param, no_loc_load_with_δ_ed_eq_pf_param, with_loc_load_no_δ_ed_eq_pf_param, with_loc_load_with_δ_ed_eq_pf_param )
    
end

#-----------------------------------------------------
#-----------------------------------------------------


""" nll: no local loads """
function get_c_integrated_nll_dyn_iip_pf_param(
    netd;
    Idxs_type = :Idxs_hybrid ,
    no_control_device = false )

    """

    no_control_device = false

    Idxs_type = :Idxs_im

    Idxs_type = :Idxs_industrial

    Idxs_type = :Idxs_hybrid

    #-----------------------------------------------
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_t2_sauer_system_matrix

    #-----------------------------------------------
    
    netd =
        NetworkData( dynamics_case()... )

    #-----------------------------------------------

    """
    #-----------------------------------------------
    
    loc_load_exist =
        loc_load_exist_bool(netd)

    # ----------------------------------------------

    """
    hybrid model is used to get gens δ_ω_ed_dash_eq_dash
    from  a standard powerflow for the network

    """

    named_tup_pf_result, state =
        get_pf_in_named_tuple_with_state(
            netd)

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    # ----------------------------------------------

    pf_net_param = get_integrated_streamlined_powerflow_net_parameters( netd; Idxs_type = Idxs_type, no_control_device = no_control_device )

    # ----------------------------------------------
    
    models_gens_nodes_some_vars_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ω, :ed_dash, :eq_dash] )       

    if Idxs_type == :Idxs_hybrid
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs =
            models_gens_nodes_some_vars_Idxs.Idxs_hybrid

        state = state

        δ_ω_ed_dash_eq_dash = 
            get_gen_nodes_δ_ω_ed_dash_eq_dash(
                state,
                nodes_δ_ω_ed_dash_eq_dash_Idxs )
                
        idx_and_Idx =
            get_hybrid_pf_etc_idx_and_Idx(netd)

        
    elseif Idxs_type == :Idxs_industrial
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs =
            models_gens_nodes_some_vars_Idxs.Idxs_industrial

        state = industrial_model_init_operationpoint(
            netd,
            bus_dict_init
            ;pure =
                :pure,
            no_control_device =
                false )
        
        δ_ω_ed_dash_eq_dash = 
            get_gen_nodes_δ_ω_ed_dash_eq_dash(
                state,
                nodes_δ_ω_ed_dash_eq_dash_Idxs )

        
        idx_and_Idx = get_industrial_model_pf_idx_and_Idx(netd; no_control_device = no_control_device)

    elseif Idxs_type == :Idxs_im
        
        nodes_δ_ω_ed_dash_eq_dash_Idxs =
            models_gens_nodes_some_vars_Idxs.Idxs_im

        state = im_model_init_operationpoint(
            netd,
            bus_dict_init )

        δ_ω_ed_dash_eq_dash = 
            get_gen_nodes_δ_ω_ed_dash_eq_dash(
                state,
                nodes_δ_ω_ed_dash_eq_dash_Idxs )

        
        idx_and_Idx =
            get_im_model_pf_idx_and_Idx(netd)
        
    else
        nothing
    end
    
    #--------------------------------------------   
    #--------------------------------------------
    
    _, pf_idx_and_state, pf_param_views,_, pf_Idx, pf_and_dyn_idx_and_Idx, pf_net_misc = pf_net_param

    # Ybus, Ynet, nodes_node_idx_and_incident_edges_other_node_idx, edges_Ybr_cal, edges_orientation = pf_net

    slack_vh, gens_vh, gens_Idx_and_vh, non_slack_gens_Idx_and_vh, _, _, _ = pf_idx_and_state

    ra_Xd_Xq, ra_Xd_dash_Xq_dash, ra_Xd_Xq_Xd_dash_Xq_dash, _, _, _, _ = pf_param_views

    # load_trans_nodes_Idx_and_vlimits = pf_limits

    _,nodes_u_Idx,_, _, _, _, _, red_vh_θh_idx,_, _,_ = pf_Idx
    
    # ---------------------------------------------

    # loc_load_exist =
    #     pf_and_dyn_idx_and_Idx.loc_load_exist

    # ---------------------------------------------

    net_comp_type_idx =
        pf_and_dyn_idx_and_Idx.net_comp_type_idx


    slack_gens_nodes_idx =
        net_comp_type_idx.slack_gens_nodes_idx

    non_slack_gens_nodes_idx =
        net_comp_type_idx.non_slack_gens_nodes_idx 
    
    gens_nodes_idx =
        net_comp_type_idx.gens_nodes_idx

    non_gens_nodes_idx =
        net_comp_type_idx.non_gens_nodes_idx

    gens_with_loc_load_idx =
        net_comp_type_idx.gens_with_loc_load_idx

    load_nodes_idx =
        net_comp_type_idx.load_nodes_idx

    transmission_nodes_idx =
        net_comp_type_idx.transmission_nodes_idx

    load_trans_nodes_idx =
        net_comp_type_idx.load_trans_nodes_idx
        
    # -----------------------------------------

    gens_with_loc_load_idx =
        loc_load_exist == false ?
        [] : 
        net_comp_type_idx.gens_with_loc_load_idx

    all_nodes_idx =
        net_comp_type_idx.all_nodes_idx

    # ----------------------------------------
        
    nodes_size = sum(
        [length(gens_nodes_idx),
         length(non_gens_nodes_idx) ])
        
    # ----------------------------------------

    dict_n2s =
        net_comp_type_idx.dicts_net_to_streamlined_idx

    n2s_slack_gens_idx =
        dict_n2s.dict_n2s_slack_gens_idx

    n2s_non_slack_gens_idx =
        dict_n2s.dict_n2s_non_slack_gens_idx
    
    n2s_gens_idx =
        dict_n2s.dict_n2s_gens_idx
    
    n2s_non_gens_idx =
        dict_n2s.dict_n2s_non_gens_idx

    n2s_load_idx =
        dict_n2s.dict_n2s_load_idx
    
    n2s_gens_with_loc_load_idxs =
        dict_n2s.dict_n2s_gens_with_loc_load_idxs

    n2s_transmission_idxs =
        dict_n2s.dict_n2s_transmission_idxs
    
    n2s_all_nodes_idx =
        dict_n2s.dict_n2s_all_nodes_idx
    
    # -------------------------------------------

    red_vh_Idxs =
        idx_and_Idx.red_vh_Idxs

    red_θh_Idxs =
        idx_and_Idx.red_θh_Idxs

    non_slack_gens_θh_idx2Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx
    
    non_slack_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_slack_gens_θh_idx2Idx_in_Idx
    
    non_gens_θh_idx2Idx =
        idx_and_Idx.non_gens_θh_idx2Idx
    
    non_gens_θh_idx2Idx_in_Idx =
        idx_and_Idx.non_gens_θh_idx2Idx_in_Idx
            
    # -----------------------------------------

    vec_Idx = pf_and_dyn_idx_and_Idx.vec_Idx
    
    vec_Idx_δ_ω_ed_dash_eq_dash =
        vec_Idx.vec_Idx_gens_nodes_δ_ω_ed_dash_eq_dash
    
    # ----------------------------------------

    NL_pf_para_Idxs = vec_Idx.NL_pf_para_Idxs
    
    #--------------------------------------------

    NL_pf_Idxs =
        NL_pf_para_Idxs.no_loc_load_with_δ_ed_NL_pf_para_Idxs

    P_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_gens_NL_para_Idxs

    Q_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_gens_NL_para_Idxs

    P_non_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_P_non_gens_NL_para_Idxs

    Q_non_gens_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_Q_non_gens_NL_para_Idxs

    δ_ed_eq_pf_NL_para_Idxs =
        NL_pf_Idxs.no_loc_load_with_δ_ed_δ_ed_eq_pf_NL_para_Idxs        
    
    # ------------------------------------------
        
    δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed =
        NL_pf_para_Idxs.δ_ω_ed_dash_eq_Idxs_in_flattend_δ_ω_ed_dash_eq
    
    # ------------------------------------------

    pf_NL_P_Q_para =
        pf_net_misc.pf_NL_P_Q_para
    
    P_gens =
        pf_NL_P_Q_para.P_gens_NL_para
    
    Q_gens =
        pf_NL_P_Q_para.Q_gens_NL_para
        
    P_non_gens =
        pf_NL_P_Q_para.P_non_gens_NL_para
    
    Q_non_gens =
        pf_NL_P_Q_para.Q_non_gens_NL_para
    
    # --------------------------------------------

    no_loc_load_with_δ_ed_eq_pf_param =
        ComponentArray(
            P_gens = P_gens,
            Q_gens = Q_gens,
            P_non_gens = P_non_gens,
            Q_non_gens = Q_non_gens,       
            δ_ed_eq_pf =
                [δ_ω_ed_dash_eq_dash...;])
        
    # -------------------------------------------

    PQ_gens_pad =
        [ idx ∈ gens_nodes_idx ?
        [ P_gens[n2s_gens_idx[idx]],
          Q_gens[n2s_gens_idx[idx]] ] :
              [0.0, 0.0]
          for idx in 1:nodes_size ]

    PQ_non_gens_pad =
        [ idx ∈ non_gens_nodes_idx ?
        [ P_non_gens[n2s_non_gens_idx[idx]],
          Q_non_gens[n2s_non_gens_idx[idx]] ] :
              [0.0, 0.0]
          for idx in 1:nodes_size ]

    P_gens_pad = first.(PQ_gens_pad)
    
    Q_gens_pad = second.(PQ_gens_pad)

    P_non_gens_pad = first.(PQ_non_gens_pad)
    
    Q_non_gens_pad = second.(PQ_non_gens_pad)

    #-------------------------------------

    dims_no_loc_load_with_δ_ed_eq_padded_PQ =
        [length(P_gens_pad),
         length(Q_gens_pad),
         length(P_non_gens_pad),
         length(Q_non_gens_pad),
         length( [δ_ω_ed_dash_eq_dash...;]  ) ]

    _,_, no_loc_load_with_δ_ed_eq_padded_PQ_Idx =
        create_size_offset_Idx(
            dims_no_loc_load_with_δ_ed_eq_padded_PQ ;
            counter = 0)

    P_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[1]

    Q_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[2]

    P_non_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[3]

    Q_non_gens_pad_Idx =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[4]

    δ_ω_ed_dash_eq_Idx_in_flat_pad =
        no_loc_load_with_δ_ed_eq_padded_PQ_Idx[5]

    #---------------------------------------------

    loc_load_exist =
        pf_and_dyn_idx_and_Idx.loc_load_exist
    
    # ---------------------------------------------
    
    # x0_vh = [ idx ∈ gens_nodes_idx ?
    #     gens_vh[n2s_gens_idx[idx]] : 1.0
    #           for idx in 1:nodes_size ]
    
    # x0_θh =
    #     [0.0
    #      for idx in
    #          1:nodes_size ]
    
    # uh_0  =
    #     x0_vh .* exp.(im * x0_θh)
    
    # vh_θh_0 =
    #     vcat( x0_vh, x0_θh )
     
    # red_vh_θh_0 =
    #     vh_θh_0[ red_vh_θh_idx  ]


    uh_0  = x_from_xr_xi.(
            [state[idx]
             for idx in
                 nodes_u_Idx])

    vh_θh_0 = [ abs.( uh_0  ); angle.( uh_0  )]
     
    red_vh_θh = vh_θh_0[ red_vh_θh_idx  ]

    
    # gens_uh_0 = [ uh_0[  idx  ]
    #             for idx  in 1:nodes_size
    #                 if idx ∈ gens_nodes_idx]

    
    gens_uh_0 =  uh_0[ gens_nodes_idx ]
                
    # --------------------------------------------
    
    #  named_tuple_with_state = get_pf_in_named_tuple_with_state( netd)
     
    gens_idq_0 =
        [  get_a_gen_dyn_idq(
            vh, θh, δ_ω_ed_eq...,
            ra_Xd_dash_Xq_dash... )
          for ( vh, θh, δ_ω_ed_eq,
                ra_Xd_dash_Xq_dash ) in
              zip( abs.( gens_uh_0  ),
                   angle.( gens_uh_0  ),
                   δ_ω_ed_dash_eq_dash,
                   ra_Xd_dash_Xq_dash ) ]

    gens_i_d_0 = first.( gens_idq_0 )
        
    gens_i_q_0 = second.( gens_idq_0 )
    
    dims_gens_idq =
        length.( [ gens_i_d_0, gens_i_q_0 ]  )

     _,_,  gens_idq_Idx =
         create_size_offset_Idx(
             dims_gens_idq ;
             counter = 0)

    gens_id_Idx = gens_idq_Idx[1]
        
    gens_iq_Idx = gens_idq_Idx[2]

    gens_idq_0_flat = [ gens_i_d_0 ; gens_i_q_0 ]

    # -------------------------------------------
    
    dims_red_vh_θh_id_iq_CA =
        length.([ red_vh_θh,
                  gens_i_d_0,
                  gens_i_q_0 ])

     _,_, red_vh_θh_id_iq_CA_Idx =
        create_size_offset_Idx(
            dims_red_vh_θh_id_iq_CA;
            counter = 0)

    red_vh_θh_CA_Idx =
        red_vh_θh_id_iq_CA_Idx[1]

    gen_id_CA_Idx =
        red_vh_θh_id_iq_CA_Idx[2]

    gen_iq_CA_Idx =
        red_vh_θh_id_iq_CA_Idx[3]

    
    # -------------------------------------------
    

   """
    t_δ = first.( δ_ω_ed_dash_eq_dash )

    t_gens_idq_0 = x_from_xr_xi.( gens_idq_0 )

    t_gens_idq_net_0 = t_gens_idq_0 .*
        exp.(im * (t_δ  .- pi/2) )

    t_gens_S = gens_uh_0 .* conj.(t_gens_idq_net_0 )

     """
    
    # -----------------------------------------

    dims_red_vh_θh_0_idq_0 =
        length.([ red_vh_θh, gens_idq_0_flat ])

     _,_, red_vh_θh_0_idq_0_Idx =
        create_size_offset_Idx(
            dims_red_vh_θh_0_idq_0;
            counter = 0)

    flat_red_vh_θh_0_Idx =
        red_vh_θh_0_idq_0_Idx[1]

    flat_idq_0_Idx =
        red_vh_θh_0_idq_0_Idx[2]
    
    # -------------------------------------------

    red_vh_θh_0_idq = [ red_vh_θh; gens_idq_0_flat ]
    
    # --------------------------------------------
    
    model_pf_idx_and_Idx =
        pf_net_misc.model_pf_idx_and_Idx
    
    vh_θh_Idxs_set =
        model_pf_idx_and_Idx.full_nodes_types_Idxs_idx2Idx_etc

    full_slack_gens_vh_Idxs =
        vh_θh_Idxs_set.full_slack_gens_vh_Idxs

    full_non_slack_gens_vh_Idxs =
        vh_θh_Idxs_set.full_non_slack_gens_vh_Idxs

    full_non_gens_nodes_vh_Idxs =
        vh_θh_Idxs_set.full_non_gens_nodes_vh_Idxs

    full_slack_gens_θh_Idxs =
        vh_θh_Idxs_set.full_slack_gens_θh_Idxs

    full_non_slack_gens_θh_Idxs =
        vh_θh_Idxs_set.full_non_slack_gens_θh_Idxs
    
    full_non_gens_nodes_θh_Idxs =
        vh_θh_Idxs_set.full_non_gens_nodes_θh_Idxs

    # -------------------------------------------

    dims_full_vh_θh_0_idq_0 =
        length.([ vh_θh_0, gens_idq_0_flat ])

     _,_, full_vh_θh_0_idq_Idx =
        create_size_offset_Idx(
            dims_full_vh_θh_0_idq_0;
            counter = 0)

    full_vh_θh_Idx = full_vh_θh_0_idq_Idx[1]

    full_idq_Idx   = full_vh_θh_0_idq_Idx[2]
    
    # --------------------------------------------

    # vh_θh_0_idq = [ vh_θh_0; gens_idq_0_flat ]
    
    # --------------------------------------------

    red_ΔPQ  = similar( red_vh_θh )
    
    red_ΔPQ_Δidq  = similar( red_vh_θh_0_idq )
    
    # ------------------------------------------- 
    # some selected param
    # -------------------------------------------
    
    net_idxs =
        (; slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         load_nodes_idx,
         transmission_nodes_idx,
         non_gens_nodes_idx,
         load_trans_nodes_idx,
         gens_with_loc_load_idx)
    
    n2s_idxs =
        (; n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_load_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_transmission_idxs,
         n2s_all_nodes_idx )  

    gens_uh_Q_from_red_sol_para =
        (; slack_vh,
         gens_vh,
         non_slack_gens_Idx_and_vh,
         ra_Xd_dash_Xq_dash,         
         red_vh_Idxs,
         red_θh_Idxs,
         red_vh_θh_idx,
         flat_red_vh_θh_0_Idx,
         flat_idq_0_Idx,
         gens_id_Idx,
         gens_iq_Idx,         
         non_gens_θh_idx2Idx,
         non_slack_gens_θh_idx2Idx,         
         net_idxs,
         n2s_idxs  )
    
    idx2Idx =
        (; non_slack_gens_θh_idx2Idx,
          non_slack_gens_θh_idx2Idx_in_Idx,
          non_gens_θh_idx2Idx,
          non_gens_θh_idx2Idx_in_Idx)

    red_vh_θh_0 = deepcopy(red_vh_θh)
    
    red_vh_θh_id_iq_CA_comp =
        (;red_vh_θh_0,
         gens_i_d_0,
         gens_i_q_0 )

    red_vh_θh_id_iq_CA_Idxs =
        (;red_vh_θh_CA_Idx,
         gen_id_CA_Idx,
         gen_iq_CA_Idx,
         red_vh_θh_id_iq_CA_Idx )    
    pois =
        (; P_gens,
         Q_gens,
         P_non_gens,
         Q_non_gens,
         δ_ω_ed_dash_eq_dash)

    pois_Idxs =
        (; P_gens_NL_para_Idxs,
         Q_gens_NL_para_Idxs,
         P_non_gens_NL_para_Idxs,
         Q_non_gens_NL_para_Idxs,
         δ_ed_eq_pf_NL_para_Idxs )
    
    kwd_net_param =
        (; pf_net_param,
         gens_idq_Idx,
         red_vh_θh_0_idq_0_Idx,
         full_vh_θh_0_idq_Idx )
    
    nll_dyn_iip_pf_param =
        (; 
         red_ΔPQ_Δidq,
         red_vh_θh_0_idq,
         gens_uh_Q_from_red_sol_para,
         kwd_net_param
         )
    
    # return (; nll_dyn_iip_pf_param, pois, pois_Idxs)
    
    flat_δ_ω_ed_dash_eq_dash =
        [δ_ω_ed_dash_eq_dash...;]
    
    return (; nll_dyn_iip_pf_param,
            flat_δ_ω_ed_dash_eq_dash,
            pois,
            pois_Idxs,
            idx2Idx,
            red_vh_θh_id_iq_CA_comp,
            red_vh_θh_id_iq_CA_Idxs) 
end


function get_dynamics_ode_model_and_powerflow_Idxs_set(
    netd;
    Idxs_type =
        :Idxs_im,
    no_control_device =
        false,
    all_wt_size_and_dims =
        false )
    
    """
    # dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    # dynamics_case =
        case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    
    netd =
        NetworkData( dynamics_case()... )
    
    pf_alg          = NewtonRaphson()
    
    # ode_alg       = Rodas4()    
    # algr_name     = "rodas4" 
    
    ode_alg         = ImplicitMidpoint()    
    algr_name       = "ImplicitMidpoint"

    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial

    Idxs_type = :Idxs_im 
    
    no_control_device = false
        
    only_gen        = false

    dt              = 0.01
   
    abstol          = 1e-14
    reltol          = 1e-14
    sim_timespan    = (0.0, 10.0)

    """
    
    #-----------------------------------------------
    
    
    n2s_idxs = n2s_streamlined_idx =
        get_dict_n2s_streamlined_idx(netd)
    
    n2s_slack_gens_idx =
        n2s_idxs.n2s_slack_gens_idx
    
    n2s_non_slack_gens_idx =
        n2s_idxs.n2s_non_slack_gens_idx
    
    n2s_gens_idx = n2s_idxs.n2s_gens_idx
    
    n2s_non_gens_idx =
        n2s_idxs.n2s_non_gens_idx
    
    n2s_gens_with_loc_load_idxs =
        n2s_idxs.n2s_gens_with_loc_load_idxs
    
    n2s_all_nodes_idx =
        n2s_idxs.n2s_all_nodes_idx

    #-----------------------------------------------
    
    nodes_types_idxs = net_nodes_types_idxs =
        get_net_nodes_type_idxs(
            netd)
    
    slack_gens_nodes_idx =
        nodes_types_idxs.slack_gens_nodes_idx
    
    non_slack_gens_nodes_idx =
        nodes_types_idxs.non_slack_gens_nodes_idx
    
    gens_nodes_idx =
        nodes_types_idxs.gens_nodes_idx
    
    non_gens_nodes_idx =
        nodes_types_idxs.non_gens_nodes_idx
    
    gens_with_loc_load_idx =
        nodes_types_idxs.gens_nodes_with_loc_loads_idx
    
    all_nodes_idx = nodes_types_idxs.all_nodes_idx


    gens_nodes_with_loc_loads_idx =
        gens_with_loc_load_idx

    ##

    load_nodes_idx =
        nodes_types_idxs.load_nodes_idx
    
    transmission_nodes_idx =
        nodes_types_idxs.transmission_nodes_idx
    
    #-----------------------------------------------
    #-----------------------------------------------

    node_size = length( all_nodes_idx )
    
    no_of_gens = length(gens_nodes_idx  )
        
    no_of_non_gens = length( non_gens_nodes_idx  )

    no_gens_with_loc_load_idx = length(gens_with_loc_load_idx)

    ##

    no_of_load_nodes = length(load_nodes_idx)
    
    no_of_transmission_nodes = length(transmission_nodes_idx) 

    #-----------------------------------------------    
    
    syms_id =
        [ :id ]
    
    syms_iq =
        [ :iq  ]

    syms_id_iq = [ :id, :iq ]
    
    syms_id_iq_pg_vh =
        [:id, :iq, :pg, :vh ]
    
    syms_ωs_ωref0_vref0_porder0 =
        [ :ωs, :ωref0, :vref0, :porder0 ]
    
    syms_a_gen_vh_θh =
        [:vh :θh ]

    syms_a_gen_δ_ed_dash_eq_dash =
        [ :δ, :ed_dash, :eq_dash]
    
    #-----------------------------------------------    
    
    dim_vh =  dim_θh = dim_ur =  dim_ui = node_size
    
    dim_Pg =  dim_Qg = no_of_gens

    dim_Png =  dim_Qng = no_of_non_gens

    dim_Pgll = dim_Qgll = no_gens_with_loc_load_idx

    dim_gens_id =  dim_gens_iq = no_of_gens
     
    #-----------------------------------------------    

    dim_a_node_ur_ui =
        length([:ur :ui])

    dim_a_node_vh_θh =
        length([:vh :θh])

    dim_a_gen_id_iq =
        length( syms_id_iq )
    
    dim_id_iq_pg_vh =
        length(syms_id_iq_pg_vh)
    
    dim_ωs_ωref0_vref0_porder0 =
        length(syms_ωs_ωref0_vref0_porder0)
    
    dim_a_gen_vh_θh =
        length(syms_a_gen_vh_θh)

    dim_a_gen_δ_ed_dash_eq_dash =
        length(syms_a_gen_δ_ed_dash_eq_dash)
    
    #-----------------------------------------------


    dim_flat_net_ur = node_size
    
    dim_flat_net_ui = node_size
    
    
    dim_flat_net_vh = node_size
    
    dim_flat_net_θh = node_size


    dim_flat_net_ur_ui = node_size * dim_a_node_ur_ui
    
    
    dim_flat_net_vh_θh = node_size * dim_a_node_vh_θh
    
    dim_flat_gens_id_iq_pg_vh =
        no_of_gens * dim_id_iq_pg_vh
        
    
    dim_flat_gens_ωs_ωref0_vref0_porder0 =
        no_of_gens * dim_ωs_ωref0_vref0_porder0
        
    
    dim_flat_gens_vh_θh =
        no_of_gens * dim_a_gen_vh_θh
        

    dim_flat_gens_δ_ed_dash_eq_dash =
        no_of_gens * dim_a_gen_δ_ed_dash_eq_dash

    ###############################################
    
    dim_flat_non_gens_vh_θh =
        no_of_non_gens * dim_a_node_vh_θh

    dim_flat_gens_id_iq =
        dim_gens_id + dim_gens_iq

    dim_flat_vh_θh_id_iq =
        dim_flat_net_vh_θh + dim_flat_gens_id_iq


    #-----------------------------------------------
    
    dim_dyn_pf_fun_flat_para =
        sum([ dim_Pg , dim_Qg,
              dim_Png, dim_Qng,
              dim_flat_gens_δ_ed_dash_eq_dash,
              dim_Pgll, dim_Qgll  ])

    dim_intra_dyn_pf_flat_para =
        dim_dyn_pf_fun_flat_para
    
    #-----------------------------------------------

    
    dim_dyn_ode_pf_flat_para =
        sum([ dim_Pg , dim_Qg,
              dim_Png, dim_Qng,
              dim_Pgll, dim_Qgll  ])

    dim_init_dyn_pf_flat_para =
        dim_dyn_ode_pf_flat_para

    
    
    #-----------------------------------------------
    
    dim_a_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh =
        dim_ωs_ωref0_vref0_porder0 + dim_id_iq_pg_vh

    #-----------------------------------------------
    #-----------------------------------------------

    Pg_Qg_Png_Qng_Pgll_Qgll_Idxs=
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_Pg, dim_Qg,
              dim_Png, dim_Qng,
              dim_Pgll, dim_Qgll ]
            ; dims_given = true )

    (Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs 

    #-----------------------------------------------

    
    pf_PQ_δ_etc_Idxs =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_Pg, dim_Qg,
              dim_Png, dim_Qng,
              dim_flat_gens_δ_ed_dash_eq_dash,
              dim_Pgll, dim_Qgll ]
            ; dims_given = true )

    ( dyn_P_gens_Idxs,
      dyn_Q_gens_Idxs,
      dyn_P_non_gens_Idxs,
      dyn_Q_non_gens_Idxs,
      dyn_δ_ed_eq_pf_Idxs,
      dyn_P_g_loc_load_Idxs,
      dyn_Q_g_loc_load_Idxs) =
          pf_PQ_δ_etc_Idxs 
        
                
    #-----------------------------------------------

    vh_θh_id_iq_Idx = get_vars_or_paras_Idxs_in_flattend(
        [ dim_vh,
          dim_θh,
          dim_gens_id,
          dim_gens_iq  ]
            ; dims_given = true )

    (vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx) =
         vh_θh_id_iq_Idx

    #-----------------------------------------------

    gens_vh_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_a_gen_vh_θh
              for a_gen in
                  gens_nodes_idx  ]
            ; dims_given = true )
    
    #-----------------------------------------------


    gens_ur_ui_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_a_node_ur_ui
              for a_gen in
                  gens_nodes_idx  ]
            ; dims_given = true )
    
    #-----------------------------------------------


    non_gens_vh_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_a_node_vh_θh
              for a_non_gen in
                  non_gens_nodes_idx  ]
            ; dims_given = true )

    #-----------------------------------------------


    non_gens_ur_ui_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_a_node_ur_ui
              for a_non_gen in
                  non_gens_nodes_idx  ]
            ; dims_given = true )

    #-----------------------------------------------

    gens_ωs_ωref0_vref0_porder0_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_ωs_ωref0_vref0_porder0
              for a_gen in
                  gens_nodes_idx  ]
            ; dims_given = true )
    
    #-----------------------------------------------

    gens_dyn_id_iq_pg_vh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_id_iq_pg_vh
              for a_gen in
                  gens_nodes_idx  ]
            ; dims_given = true )
    
    #-----------------------------------------------

    gens_δ_ed_dash_eq_dash_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [  dim_a_gen_δ_ed_dash_eq_dash
              for a_gen in
                  gens_nodes_idx  ]
            ; dims_given = true )

    δ_ed_dash_eq_dash_Idxs_in_flattend =
        gens_δ_ed_dash_eq_dash_idx_in_Idx
    
    #-----------------------------------------------

    vh_θh_idx_in_Idx  =
        get_vars_or_paras_Idxs_in_flattend(
            [  dim_a_node_vh_θh
              for a_node in
                  all_nodes_idx  ]
            ; dims_given = true )

    ur_ui_idx_in_Idx = vh_θh_idx_in_Idx
    
    #-----------------------------------------------

    a_gen_vtf_vh_θh_δ_ed_dash_eq_dash_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_a_gen_vh_θh,  dim_a_gen_δ_ed_dash_eq_dash ]
            ; dims_given = true )
    
    ( a_gen_vtf_vh_θh_Idx,
      a_gen_vtf_δ_ed_dash_eq_dash_Idx ) =
          a_gen_vtf_vh_θh_δ_ed_dash_eq_dash_Idx
    
    #-----------------------------------------------

    a_gen_per_var_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_ωs_ωref0_vref0_porder0,
              dim_id_iq_pg_vh ]
            ; dims_given = true )
    
    per_gen_ωs_ωref0_vref0_porder0_Idx, per_gen_id_iq_pg_vh_Idx =
        a_gen_per_var_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx
    
    #-----------------------------------------------

    per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_flat_gens_ωs_ωref0_vref0_porder0,
              dim_flat_gens_id_iq_pg_vh ]
            ; dims_given = true )

    per_para_ωs_ωref0_vref0_porder0_Idx, per_para_id_iq_pg_vh_Idx =
        per_para_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx
    
    #-----------------------------------------------

    per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [  dim_a_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh
              for a_gen in
                  gens_nodes_idx  ]
            ; dims_given = true )
    
    #-----------------------------------------------

    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_flat_gens_ωs_ωref0_vref0_porder0,
              dim_dyn_pf_fun_flat_para ]
            ; dims_given = true )

    f_ωs_ωref0_vref0_porder0_Idx, f_dyn_pf_para_Idx =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para_Idx

    
    #-----------------------------------------------

    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_ode_pf_para_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_flat_gens_ωs_ωref0_vref0_porder0,
              dim_dyn_ode_pf_flat_para ]
            ; dims_given = true )

    _, f_dyn_ode_pf_para_Idx =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_ode_pf_para_Idx 
        
    #-----------------------------------------------

    ####
    
    intra_pf_flat_vh_θh_wt_dyn_pf_flat_para_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_flat_net_vh_θh,
              dim_dyn_pf_fun_flat_para ]
            ; dims_given = true )

    intra_flat_vh_θh_Idx, intra_dyn_pf_flat_para_Idx =
        intra_pf_flat_vh_θh_wt_dyn_pf_flat_para_Idx

    intra_flat_ur_ui_Idx = intra_flat_vh_θh_Idx

    
    flat_ur_flat_ui_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_flat_net_ur,
              dim_flat_net_ui ]
            ; dims_given = true )

    flat_ur_Idx, flat_ui_Idx = flat_ur_flat_ui_Idx

    flat_vh_Idx = flat_ur_Idx

    flat_θh_Idx = flat_ui_Idx

    flat_vh_flat_θh_Idx = flat_ur_flat_ui_Idx

    #-----------------------------------------------


    #---------------------------------------------
    # im pure states idxs
    #---------------------------------------------

    im_pure_states_syms =
        [ :δ, :ω, :ed_dash, :eq_dash, :xg1, :xg2,
          :vr1, :vr2, :vf_tilade ]
    
    models_types_gens_pure_states_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                im_pure_states_syms )

    vec_im_pure_states_Idxs =
        models_types_gens_pure_states_Idxs.Idxs_im

    im_pure_states_idxs =
        [ vec_im_pure_states_Idxs...; ]

    #---------------------------------------------
    # im algebraic states idxs
    #---------------------------------------------

    im_τm_tilade_vf_states_syms =
        [ :τm_tilade, :vf ]
    
    models_types_gens_τm_tilade_vf_states_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                im_τm_tilade_vf_states_syms  )

    vec_im_τm_tilade_vf_states_Idxs =
        models_types_gens_τm_tilade_vf_states_Idxs.Idxs_im

    im_τm_tilade_vf_states_idxs =
        [ vec_im_τm_tilade_vf_states_Idxs...; ]

    #---------------------------------------------
    # sauer
    #---------------------------------------------

    sauer_states_syms =
        [ :δ, :ω, :eq_dash, :ed_dash,
          :vf_tilade, :vr1, :vr2 ]

    vec_im_sauer_states_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                sauer_states_syms ).Idxs_im

    flat_im_sauer_states_Idxs =
        [vec_im_sauer_states_Idxs...;]

    sauer_δ_eq_dash_ed_dash_syms =
        [ :δ, :eq_dash, :ed_dash ]

    vec_im_sauer_δ_eq_dash_ed_dash_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                sauer_δ_eq_dash_ed_dash_syms ).Idxs_im

    flat_im_sauer_δ_eq_dash_ed_dash_Idxs =
        [vec_im_sauer_δ_eq_dash_ed_dash_Idxs...;]

    #-----------------------------------------------

    
    #-----------------------------------------------
    ################################################
    #-----------------------------------------------    
    
    integrated_pf_vars_and_para_idx =
        get_a_model_integrated_pf_vars_and_para_idx(
            netd;
            Idxs_type =
                Idxs_type,
            no_control_device =
                false )

    #-----------------------------------------------
    
    intg_types_Idxs_etc =
        integrated_pf_vars_and_para_idx.intg_types_Idxs_etc


    gens_vh_Idxs =
        intg_types_Idxs_etc.intg_gens_vh_Idxs
    
    gens_θh_Idxs =
        intg_types_Idxs_etc.intg_gens_θh_Idxs
    
    non_gens_nodes_vh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_vh_Idxs
    
    non_gens_nodes_θh_Idxs =
        intg_types_Idxs_etc.intg_non_gens_nodes_θh_Idxs

    gen_id_Idxs =
        intg_types_Idxs_etc.intg_gen_id_Idxs
    
    gen_iq_Idxs =
        intg_types_Idxs_etc.intg_gen_iq_Idxs
    
    # ----
   
    dyn_pf_vars_Idxs = 
        (;
         gens_vh_Idxs,
         gens_θh_Idxs,
         non_gens_nodes_vh_Idxs,
         non_gens_nodes_θh_Idxs,
         gen_id_Idxs,
         gen_iq_Idxs )

    dyn_pf_vars_Idxs_kwd_para =
        (; dyn_pf_vars_Idxs, )
    
    #-----------------------------------------------

     u_ur_ui_Idx_in_state =
         integrated_pf_vars_and_para_idx.u_ur_ui_Idx_in_state


   (; slack_ur_ui_Idx_in_state, 
    non_slack_ur_ui_Idx_in_state, 
    ur_ui_Idx_in_state, 
    nodes_u_Idx,
    nodes_u_Idx_in_ranges ) =
        u_ur_ui_Idx_in_state
    
    # nodes_u_Idx_in_ranges =
    #     u_ur_ui_Idx_in_state.nodes_u_Idx_in_ranges


    (;
     gens_nodes_u_Idx_in_ranges,
     non_gens_nodes_u_Idx_in_ranges ) =
         gens_and_non_gens_u_Idx_in_ranges(
             all_nodes_idx,
             gens_nodes_idx,
             non_gens_nodes_idx,
                 nodes_u_Idx_in_ranges)
    
    #-----------------------------------------------


    pf_PQ_δ_etc_Idxs = 
        get_a_model_integrated_dyn_pf_PQ_δ_etc_Idxs(
            netd;
            δ_etc = [:δ, :ed_dash, :eq_dash],
            δ_etc_first = false )

    dyn_pf_P_Q_δ_etc_kwd_para_Idxs =
        pf_PQ_δ_etc_Idxs.dyn_pf_PQ_δ_etc_Idxs


    (;     
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
    ) = dyn_pf_P_Q_δ_etc_kwd_para_Idxs
    
    #-----------------------------------------------
    
   im_indices_and_conversion_dict =
         get_im_indices_and_conversion_dict(
             netd  )

    indices_and_conversion_dict =
        im_indices_and_conversion_dict

    im_vars_Idx_in_state =
        indices_and_conversion_dict.im_vars_Idx_in_state

    each_gens_im_vars_Idx_in_state =
        indices_and_conversion_dict.each_gens_im_vars_Idx_in_state

    nodes_ur_ui_Idx_in_state =
        indices_and_conversion_dict.nodes_ur_ui_Idx_in_state
        
    #-----------------------------------------------

    
    models_types_gens_δ_ed_dash_eq_dash_Idxs =
        get_models_gens_nodes_some_vars_Idxs(
            netd; some_state_vars =
                [:δ, :ed_dash, :eq_dash] )

    nodes_δ_ed_dash_eq_dash_Idxs_im =
        models_types_gens_δ_ed_dash_eq_dash_Idxs.Idxs_im
    
    nodes_δ_ed_dash_eq_dash_Idxs =
        nodes_δ_ed_dash_eq_dash_Idxs_im
    
    #-----------------------------------------------    
    
    flat_δ_ed_dash_eq_dash_Idxs_in_state =
        Int64[ nodes_δ_ed_dash_eq_dash_Idxs...; ]

    #-----------------------------------------------
    #-----------------------------------------------
    
    gens_nodes_collection =
        get_gens_nodes(
            netd.nodes )
    
    #-----------------------------------------------
    #-----------------------------------------------    
    
   # δ_ω_ed_dash_eq_dash_Idxs_in_flattend =
   #      get_gens_nodes_some_state_vars_Idxs_in_flattend(
   #          gens_nodes_collection
   #          ; some_state_vars =
   #              [ :δ, :ed_dash, :eq_dash ] )


    _, _, states_and_mat_Idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(
            gens_nodes_collection )

    (;
     nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ωs_ω_ref_v_ref_p_order_idxs )  =
         states_and_mat_Idxs
    
    #-----------------------------------------------
    ################################################
    #-----------------------------------------------    
    
     dyn_pf_fun_kwd_net_idxs  =  
        (;
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx
         )
    
    
    dyn_pf_fun_kwd_n2s_idxs =
        (;
         n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx )    


    Pg_Qg_Png_Qng_Pgll_Qgll_Idxs  =
        (;
         Pg_Idx,
         Qg_Idxs,
         Png_Idxs,
         Qng_Idxs,
         Pgll_Idxs,
         Qgll_Idxs )
    
    pf_PQ_δ_etc_Idxs  = dyn_pf_P_Q_δ_etc_kwd_para_Idxs =
        (;
         dyn_P_gens_Idxs,
         dyn_Q_gens_Idxs,
         dyn_P_non_gens_Idxs,
         dyn_Q_non_gens_Idxs,
         dyn_δ_ed_eq_pf_Idxs,
         dyn_P_g_loc_load_Idxs,
         dyn_Q_g_loc_load_Idxs                   
        ) 


    vh_θh_id_iq_Idx =
        (;
         vh_Idx,
         θh_Idx,
         id_Idx,
         iq_Idx )

   dyn_pf_vars_Idxs = 
        (;
         gens_vh_Idxs,
         gens_θh_Idxs,
         non_gens_nodes_vh_Idxs,
         non_gens_nodes_θh_Idxs,
         gen_id_Idxs,
         gen_iq_Idxs )

    dyn_pf_vars_Idxs_kwd_para =
        (; dyn_pf_vars_Idxs, )
    
    
    states_and_mat_Idxs =
        (;
         nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         id_iq_ph_vh_idxs,
         ωs_ω_ref_v_ref_p_order_idxs ) 
    
    #-----------------------------------------------
    #-----------------------------------------------
    
    # dis

    Idx_set_A = (;
                 gens_vh_θh_idx_in_Idx,
                 non_gens_vh_θh_idx_in_Idx,
                 gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
                 gens_dyn_id_iq_pg_vh_idx_in_Idx,
                 gens_δ_ed_dash_eq_dash_idx_in_Idx,
                 ur_ui_idx_in_Idx,
                 vh_θh_idx_in_Idx, a_gen_vtf_vh_θh_Idx,
                 a_gen_vtf_δ_ed_dash_eq_dash_Idx,
                 per_gen_ωs_ωref0_vref0_porder0_Idx,
                 per_gen_id_iq_pg_vh_Idx,
                 per_para_ωs_ωref0_vref0_porder0_Idx,
                 per_para_id_iq_pg_vh_Idx,
                 per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,
                 f_ωs_ωref0_vref0_porder0_Idx,
                 f_dyn_pf_para_Idx,
                 f_dyn_ode_pf_para_Idx,
                 δ_ed_dash_eq_dash_Idxs_in_flattend,
                 intra_dyn_pf_flat_para_Idx,
                 intra_flat_ur_ui_Idx,
                 intra_flat_vh_θh_Idx )

    
    Idx_set_B = (; slack_ur_ui_Idx_in_state,
                 non_slack_ur_ui_Idx_in_state,
                 ur_ui_Idx_in_state,
                 nodes_u_Idx,
                 nodes_u_Idx_in_ranges,
                 gens_nodes_u_Idx_in_ranges,
                 non_gens_nodes_u_Idx_in_ranges,
                 im_vars_Idx_in_state,
                 each_gens_im_vars_Idx_in_state,
                 nodes_ur_ui_Idx_in_state,
                 nodes_δ_ed_dash_eq_dash_Idxs )

    
    Idx_set_C = (;
                 dyn_pf_fun_kwd_net_idxs,
                 dyn_pf_fun_kwd_n2s_idxs,
                 Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
                 pf_PQ_δ_etc_Idxs,
                 dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
                 vh_θh_id_iq_Idx,
                 dyn_pf_vars_Idxs,
                 dyn_pf_vars_Idxs_kwd_para,
                 states_and_mat_Idxs  )

    Idx_set_misc = (;
                    flat_vh_Idx,
                    flat_θh_Idx,
                    flat_ur_Idx,
                    flat_ui_Idx,
                    flat_ur_flat_ui_Idx,
                    flat_vh_flat_θh_Idx,
                    δ_ed_dash_eq_dash_Idxs_in_flattend,
                    flat_δ_ed_dash_eq_dash_Idxs_in_state,

                    vec_im_pure_states_Idxs,
                    im_pure_states_idxs,
                    vec_im_τm_tilade_vf_states_Idxs,
                    im_τm_tilade_vf_states_idxs,
                    vec_im_sauer_states_Idxs,
                    flat_im_sauer_states_Idxs,
                    vec_im_sauer_δ_eq_dash_ed_dash_Idxs,
                    flat_im_sauer_δ_eq_dash_ed_dash_Idxs
                    )

    #------------------------------------------

        
        sizes_and_dims =
            (;
             node_size,
             no_of_gens,
             no_of_non_gens,
             no_gens_with_loc_load_idx,

             no_of_load_nodes,
             no_of_transmission_nodes,
             
             dim_vh,
             dim_θh,
             
             dim_ur,
             dim_ui,
             
             dim_Pg,
             dim_Qg,
             
             dim_Png,
             dim_Qng,
             
             dim_Pgll,
             dim_Qgll,
             
             dim_gens_id,
             dim_gens_iq,
             
             dim_a_node_ur_ui,
             dim_a_node_vh_θh,

             dim_a_gen_id_iq,
             
             dim_id_iq_pg_vh,
             
             dim_ωs_ωref0_vref0_porder0,
             
             dim_a_gen_vh_θh, 
             dim_a_gen_δ_ed_dash_eq_dash,

             dim_flat_net_ur,
             dim_flat_net_ui,

             dim_flat_net_vh,
             dim_flat_net_θh,

             dim_flat_net_ur_ui,
             dim_flat_net_vh_θh,
             dim_flat_vh_θh_id_iq,
             
             dim_flat_gens_id_iq_pg_vh, 
             dim_flat_gens_ωs_ωref0_vref0_porder0,
                          
             dim_intra_dyn_pf_flat_para,
             dim_init_dyn_pf_flat_para,
             dim_dyn_ode_pf_flat_para,
             
             dim_flat_gens_vh_θh,
             dim_flat_non_gens_vh_θh,
             
             dim_flat_gens_id_iq,
             dim_flat_gens_δ_ed_dash_eq_dash,
             
             dim_dyn_pf_fun_flat_para, 
             dim_a_gen_ωs_ωref0_vref0_porder0_id_iq_pg_vh
             )
       
    #------------------------------------------
    #------------------------------------------
    
    if all_wt_size_and_dims == false
        
        return (;
                Idx_set_A,
                Idx_set_B,
                Idx_set_C,

                Idx_set_misc,
                sizes_and_dims )
        
    else

        all_idxs_and_Idxs =
            (;
             n2s_idxs,
             nodes_types_idxs,

             slack_gens_nodes_idx,
             non_slack_gens_nodes_idx,
             gens_nodes_idx,
             non_gens_nodes_idx,
             gens_nodes_with_loc_loads_idx,
             all_nodes_idx,

             n2s_slack_gens_idx,
             n2s_non_slack_gens_idx,
             n2s_gens_idx,
             n2s_non_gens_idx,
             n2s_gens_with_loc_load_idxs,
             n2s_all_nodes_idx,

             Pg_Idx,
             Qg_Idxs,
             Png_Idxs,
             Qng_Idxs,
             Pgll_Idxs,
             Qgll_Idxs,
             dyn_P_gens_Idxs,
             dyn_Q_gens_Idxs,
             dyn_P_non_gens_Idxs,
             dyn_Q_non_gens_Idxs,
             dyn_δ_ed_eq_pf_Idxs,
             dyn_P_g_loc_load_Idxs,
             dyn_Q_g_loc_load_Idxs,

             vh_Idx,
             θh_Idx,
             id_Idx,
             iq_Idx,

             gens_vh_Idxs,
             gens_θh_Idxs,
             non_gens_nodes_vh_Idxs,
             non_gens_nodes_θh_Idxs,
             gen_id_Idxs,
             gen_iq_Idxs,


             nodes_state_Idx,
             Bx_idxs,
             Cx_idxs,
             id_iq_ph_vh_idxs,
             ωs_ω_ref_v_ref_p_order_idxs,    

             dyn_pf_fun_kwd_net_idxs,
             dyn_pf_fun_kwd_n2s_idxs,
             Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
             pf_PQ_δ_etc_Idxs,
             dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
             vh_θh_id_iq_Idx,

             dyn_pf_vars_Idxs,

             dyn_pf_vars_Idxs_kwd_para,
             states_and_mat_Idxs,

             gens_vh_θh_idx_in_Idx,
             non_gens_vh_θh_idx_in_Idx,
             gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
             gens_dyn_id_iq_pg_vh_idx_in_Idx,
             gens_δ_ed_dash_eq_dash_idx_in_Idx,
             ur_ui_idx_in_Idx,
             vh_θh_idx_in_Idx, a_gen_vtf_vh_θh_Idx,
             a_gen_vtf_δ_ed_dash_eq_dash_Idx,
             per_gen_ωs_ωref0_vref0_porder0_Idx,
             per_gen_id_iq_pg_vh_Idx,
             per_para_ωs_ωref0_vref0_porder0_Idx,
             per_para_id_iq_pg_vh_Idx,
             per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,
             f_ωs_ωref0_vref0_porder0_Idx,
             f_dyn_pf_para_Idx,
             f_dyn_ode_pf_para_Idx,
             δ_ed_dash_eq_dash_Idxs_in_flattend,
             intra_dyn_pf_flat_para_Idx,
             intra_flat_ur_ui_Idx,
             intra_flat_vh_θh_Idx,

             slack_ur_ui_Idx_in_state,
             non_slack_ur_ui_Idx_in_state,
             ur_ui_Idx_in_state,
             nodes_u_Idx,
             nodes_u_Idx_in_ranges,
             gens_nodes_u_Idx_in_ranges,
             non_gens_nodes_u_Idx_in_ranges,
             im_vars_Idx_in_state,
             each_gens_im_vars_Idx_in_state,
             nodes_ur_ui_Idx_in_state,
             nodes_δ_ed_dash_eq_dash_Idxs,
             flat_vh_Idx,
             flat_θh_Idx,
             flat_ur_Idx,
             flat_ui_Idx,
             flat_ur_flat_ui_Idx,
             flat_vh_flat_θh_Idx )

        
        return (; Idx_set_A,
                Idx_set_B,
                Idx_set_C,
                Idx_set_misc,
                all_idxs_and_Idxs,
                sizes_and_dims )
        
    end
    
end


function driver_get_dynamics_ode_model_and_powerflow_Idxs_set()
    
    pf_alg  = NewtonRaphson()
    
    # ode_alg       = Rodas4()    
    # algr_name     = "rodas4" 
    
    ode_alg         = ImplicitMidpoint()    
    algr_name       = "ImplicitMidpoint"


    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial

    Idxs_type = :Idxs_im 
    
    no_control_device = false
        
    only_gen  = false

    all_wt_size_and_dims = false

    dt              = 0.01

    sim_timespan    = (0.0, 10.0)
    
    abstol          = 1e-14

    reltol          = 1e-14

        
    """
    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plants_v6_P_rscad

    """
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    
    netd =
        NetworkData( dynamics_case()... )
    
    #-----------------------------------------------    
            
    loc_load_exist =
        loc_load_exist_bool( netd )

    #-----------------------------------------------
        
    dynamic_model_and_powerflow_Idxs_set =
        get_dynamics_ode_model_and_powerflow_Idxs_set(
            netd;
            Idxs_type =
                Idxs_type,
            no_control_device =
                no_control_device,
            all_wt_size_and_dims =
                all_wt_size_and_dims )
    
end


function get_dynamics_paras_and_Idxs_set(
    netd;
    Idxs_type = :Idxs_im, 
    
    no_control_device = false,
        
    only_gen  = false,

    all_wt_size_and_dims = false,
    
    pf_alg  = NewtonRaphson(),    
    streamlined =
        true )
    
        
    """
    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plants_v6_P_rscad
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    
    netd =
        NetworkData( dynamics_case()... )
    
    pf_alg  = NewtonRaphson()

    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial

    Idxs_type = :Idxs_im 
    
    no_control_device = false
        
    only_gen  = false

    all_wt_size_and_dims = false
    
    abstol          = 1e-14

    reltol          = 1e-14

    """

    #-----------------------------------------------    
            
    loc_load_exist =
        loc_load_exist_bool( netd )

    #-----------------------------------------------
    
    dynamic_model_and_powerflow_Idxs_set =
        get_dynamics_ode_model_and_powerflow_Idxs_set(
            netd;
            Idxs_type =
                Idxs_type,
            no_control_device =
                no_control_device,
            all_wt_size_and_dims =
                all_wt_size_and_dims )
    
    #-----------------------------------------------
    #-----------------------------------------------

    sizes_and_dims =
        dynamic_model_and_powerflow_Idxs_set.sizes_and_dims

    #-----------------------------------------------

    Idx_set_A =
        dynamic_model_and_powerflow_Idxs_set.Idx_set_A

    #-----------------------------------------------

    (;
    gens_vh_θh_idx_in_Idx,
    non_gens_vh_θh_idx_in_Idx,
    gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
    gens_dyn_id_iq_pg_vh_idx_in_Idx,
    gens_δ_ed_dash_eq_dash_idx_in_Idx,
    ur_ui_idx_in_Idx,
    vh_θh_idx_in_Idx, a_gen_vtf_vh_θh_Idx,
    a_gen_vtf_δ_ed_dash_eq_dash_Idx,
    per_gen_ωs_ωref0_vref0_porder0_Idx,
    per_gen_id_iq_pg_vh_Idx,
    per_para_ωs_ωref0_vref0_porder0_Idx,
    per_para_id_iq_pg_vh_Idx,
    per_node_ωs_ωref0_vref0_porder0_id_iq_pg_vh_Idx,
    f_ωs_ωref0_vref0_porder0_Idx,
    f_dyn_pf_para_Idx,
    f_dyn_ode_pf_para_Idx,
    δ_ed_dash_eq_dash_Idxs_in_flattend,
    intra_dyn_pf_flat_para_Idx,
    intra_flat_ur_ui_Idx,
     intra_flat_vh_θh_Idx ) =
         Idx_set_A

    #-----------------------------------------------
    
    Idx_set_B =
        dynamic_model_and_powerflow_Idxs_set.Idx_set_B

    #-----------------------------------------------
    
    (;
     slack_ur_ui_Idx_in_state,
     non_slack_ur_ui_Idx_in_state,
     ur_ui_Idx_in_state,
     nodes_u_Idx,
     nodes_u_Idx_in_ranges,
     gens_nodes_u_Idx_in_ranges,
     non_gens_nodes_u_Idx_in_ranges,
     im_vars_Idx_in_state,
     each_gens_im_vars_Idx_in_state,
     nodes_ur_ui_Idx_in_state,
     nodes_δ_ed_dash_eq_dash_Idxs ) =
         Idx_set_B
    
    #-----------------------------------------------
    #-----------------------------------------------

    Idx_set_C =
        dynamic_model_and_powerflow_Idxs_set.Idx_set_C

    #-----------------------------------------------

    (;
     dyn_pf_fun_kwd_net_idxs,
     dyn_pf_fun_kwd_n2s_idxs,
     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
     pf_PQ_δ_etc_Idxs,
     dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
     vh_θh_id_iq_Idx,
     dyn_pf_vars_Idxs,
     dyn_pf_vars_Idxs_kwd_para,
     states_and_mat_Idxs  ) =
         Idx_set_C

    #----------------------------------------------- 

    Idx_set_misc =
        dynamic_model_and_powerflow_Idxs_set.Idx_set_misc

    
    flat_δ_ed_dash_eq_dash_Idxs_in_state =
        Idx_set_misc.flat_δ_ed_dash_eq_dash_Idxs_in_state

    
    flat_ur_flat_ui_Idx =
        Idx_set_misc.flat_ur_flat_ui_Idx
    
    flat_vh_flat_θh_Idx =
        Idx_set_misc.flat_vh_flat_θh_Idx


    #
    vec_im_pure_states_Idxs =
        Idx_set_misc.vec_im_pure_states_Idxs
    
    im_pure_states_idxs =
        Idx_set_misc.im_pure_states_idxs
    
    vec_im_τm_tilade_vf_states_Idxs =
        Idx_set_misc.vec_im_τm_tilade_vf_states_Idxs
    
    im_τm_tilade_vf_states_idxs =
        Idx_set_misc.im_τm_tilade_vf_states_idxs
    
    vec_im_sauer_states_Idxs =
        Idx_set_misc.vec_im_sauer_states_Idxs
    
    flat_im_sauer_states_Idxs =
        Idx_set_misc.flat_im_sauer_states_Idxs
    
    vec_im_sauer_δ_eq_dash_ed_dash_Idxs =
        Idx_set_misc.vec_im_sauer_δ_eq_dash_ed_dash_Idxs
    
    flat_im_sauer_δ_eq_dash_ed_dash_Idxs =
        Idx_set_misc.flat_im_sauer_δ_eq_dash_ed_dash_Idxs
    
    #-----------------------------------------------
    
    (;
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx
     ) =
         dyn_pf_fun_kwd_net_idxs


    (;
     n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_all_nodes_idx ) =
         dyn_pf_fun_kwd_n2s_idxs   


    (; Pg_Idx,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs) =
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs


    (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_δ_ed_eq_pf_Idxs,
     dyn_P_g_loc_load_Idxs,
     dyn_Q_g_loc_load_Idxs                   
     ) =
         dyn_pf_P_Q_δ_etc_kwd_para_Idxs


    (;vh_Idx,
     θh_Idx,
     id_Idx,
     iq_Idx ) =
         vh_θh_id_iq_Idx


    (;
     gens_vh_Idxs,
     gens_θh_Idxs,
     
     non_gens_nodes_vh_Idxs,
     non_gens_nodes_θh_Idxs,
     
     gen_id_Idxs,
     gen_iq_Idxs ) =
         dyn_pf_vars_Idxs


    (; dyn_pf_vars_Idxs, ) =
        dyn_pf_vars_Idxs_kwd_para 

    (;
     nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ωs_ω_ref_v_ref_p_order_idxs ) =
         states_and_mat_Idxs
    
    #-----------------------------------------------

    Idx_set_misc =
        dynamic_model_and_powerflow_Idxs_set.Idx_set_misc
    
    #-----------------------------------------------    
    #-----------------------------------------------

    
    (; Ynet,
     nodes_idx_with_adjacent_nodes_idx ) =
         get_Ynet_and_nodes_idx_wt_adjacent_nodes_idx(
             netd )

    nodes_node_idx_and_incident_edges_other_node_idx =
        nodes_idx_with_adjacent_nodes_idx
    
    #-----------------------------------------------   
    #  label syms and mass_matrix
    #-----------------------------------------------
    
    net_class_names =
           make_case_buses_names(
               netd.nodes )

    #-----------------------------------------------

    im_net_states_and_var_labels =
         generate_im_model_labels(
             ;nodes =
                 netd.nodes )

    #-----------------------------------------------
    
    im_sym, im_mass_matrix =
        generate_im_sym_and_mass_matrix(
            ; nodes =
                netd.nodes )

    (;
     im_model_ode_mass_matrix,
     im_model_pf_mass_matrix ) =
         generate_im_ode_and_pf_mass_matrix(
             ; nodes =
                 netd.nodes )

    (;
     gens_nodes_im_vars_labels,
     net_bus_volts_labels ) =
         generate_im_ode_and_pf_sym(
             ; nodes =
                 netd.nodes  )

    im_ode_sym =
        gens_nodes_im_vars_labels

    im_ode_mass_matrix =
        im_model_ode_mass_matrix
    
    #-----------------------------------------------

    para_net_names_labels_syms =
        (;
         net_class_names,
         im_net_states_and_var_labels,
         im_ode_sym )

    #-----------------------------------------------    

    (;
     P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load ) =
         get_a_model_integrated_pf_PQ_param(
             netd; P_Q_only = true )
    
    #-----------------------------------------------

    gens_nodes_ra_Xd_dash_Xq_dash =
        get_gens_params_in_param_values(
            netd.nodes_param,
            netd.nodes;
            param_list =
                [ :ra, :X_d_dash, :X_q_dash ],
            gens_view_only = true )        
    
    #-----------------------------------------------    
    #-----------------------------------------------    
    
    gens_nodes_collection =
        get_gens_nodes(
            netd.nodes )
    
    #-----------------------------------------------    
    #-----------------------------------------------
    
    gens_Ax_update_parameters =
        map( get_a_im_Ax_update_parameters,
             gens_nodes_collection )
    
    #-----------------------------------------------    
    #-----------------------------------------------
    
    Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, states_and_mat_Idxs =
        create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(
            gens_nodes_collection )

    (;vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views )  =
        Ax_Bx_Cx_views
    
    (;Ax_matrix,
     Bx_matrix,
     Cx_matrix ) =
         Ax_Bx_Cx_matrix

    (; nodes_state_Idx,
     Bx_idxs,
     Cx_idxs,
     id_iq_ph_vh_idxs,
     ωs_ω_ref_v_ref_p_order_idxs )  =
         states_and_mat_Idxs
    
    init_im_plants_system_matrices_views!(
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views ),
        gens_nodes_collection )

    im_plants_system_matrices =
        (;
         vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         
         Ax_matrix,
         Bx_matrix,
         Cx_matrix )
        
    #-----------------------------------------------    
    #-----------------------------------------------
    
    integrated_pf_and_init =
         get_a_model_integrated_sta_powerflow_and_init(
             netd;
             pf_alg  =
                 pf_alg)
        
    #-----------------------------------------------
    #-----------------------------------------------
    
    nodes_name =
        integrated_pf_and_init.nodes_name

    branches_name =
        integrated_pf_and_init.branches_name

    gens_nodes_idx =
        integrated_pf_and_init.gens_nodes_idx

    #-----------------------------------------------

    named_tup_pf_result =
        integrated_pf_and_init.named_tup_pf_result

    bus_dict_init =
        named_tup_pf_result.bus_dict_init

    branch_dict_init = named_tup_pf_result.branch_dict_init

    pf_init_dict = named_tup_pf_result.pf_init_dict
    
    #-----------------------------------------------
    
    im_state =
        im_model_init_operationpoint(
            netd, bus_dict_init  )

    #-----------------------------------------------
    
    state = im_state
    
    #-----------------------------------------------
    
    im_vars_in_state =
        state[im_vars_Idx_in_state ]
    
    #-----------------------------------------------

    sim_state_x0 = state
    
    # state[ im_vars_Idx_in_state ]
    
    gens_sim_state_x0 =
        [ sim_state_x0[ idx ]
         for idx in
             nodes_state_Idx ]

    #-----------------------------------------------
    
    vh = named_tup_pf_result.Vm
    
    θh = named_tup_pf_result.Vθ

    #-----------------------------------------------

    gens_vh = vh[gens_nodes_idx]

    gens_θh = θh[gens_nodes_idx]

    gens_vh_θh =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(gens_vh , gens_θh)]
    
    #-----------------------------------------------

    non_gens_vh = vh[non_gens_nodes_idx]

    non_gens_θh = θh[non_gens_nodes_idx]

    non_gens_vh_θh =
        [[a_gen_vh, a_gen_θh]
         for (a_gen_vh, a_gen_θh) in
             zip(non_gens_vh , non_gens_θh)]

    #-----------------------------------------------
    
    gens_nodes_ωs_ωref0_vref0_porder0 =
        get_gens_nodes_ωs_ωref0_vref0_porder0(
        gens_nodes_collection,
        bus_dict_init )


    gens_nodes_τm_vf =
        get_gens_τm_vf(
            gens_nodes_collection,
            bus_dict_init )
    
    #-----------------------------------------------

    gens_nodes_δ_ed_dash_eq_dash =
        [state[idx] for idx in
             nodes_δ_ed_dash_eq_dash_Idxs ]
    
    #-----------------------------------------------

    gens_δ =
        first.( gens_nodes_δ_ed_dash_eq_dash )

    gens_ed_dash =
        second.( gens_nodes_δ_ed_dash_eq_dash )

    gens_eq_dash =
        third.( gens_nodes_δ_ed_dash_eq_dash )    
    
    #-----------------------------------------------

    gens_ra =
        first.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xd_dash =
        second.( gens_nodes_ra_Xd_dash_Xq_dash )

    gens_Xq_dash =
        third.( gens_nodes_ra_Xd_dash_Xq_dash )        
    
    #-----------------------------------------------  

    gens_dynamic_id_iq_pg_vh =
        [get_dynamic_id_iq_pg_vh_by_vhθh(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
         for (a_vh,
              a_θh,
              a_δ,
              a_ed_dash,
              a_eq_dash,
              a_ra,
              a_X_d_dash,
              a_X_q_dash ) in
                zip( gens_vh,
                     gens_θh,
                     gens_δ,
                     gens_ed_dash,
                     gens_eq_dash,
                     gens_ra,
                     gens_Xd_dash,
                     gens_Xq_dash ) ]
    
    gens_id_iq = [
        get_a_gen_dyn_idq(
            a_vh,
            a_θh,
            a_δ,
            a_ed_dash,
            a_eq_dash,
            a_ra,
            a_X_d_dash,
            a_X_q_dash )
        for (a_vh,
             a_θh,
             a_δ,
             a_ed_dash,
             a_eq_dash,
             a_ra,
             a_X_d_dash,
             a_X_q_dash ) in
            zip( gens_vh,
                 gens_θh,
                 gens_δ,
                 gens_ed_dash,
                 gens_eq_dash,
                 gens_ra,
                 gens_Xd_dash,
                 gens_Xq_dash ) ]

    gens_i_d_0 =
        first.( gens_id_iq )
    
    gens_i_q_0 =
        last.( gens_id_iq )

    
    #-----------------------------------------------    
    #-----------------------------------------------   

    gens_vd = [
        a_ed_dash +
            a_Xq_dash * a_iq -
            a_ra * a_id
        for ( a_ed_dash,
              a_Xq_dash,
              a_ra,
              a_id,
              a_iq ) in
            zip(gens_ed_dash,
                gens_Xq_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]

    gens_vq = [
        a_eq_dash - a_Xd_dash * a_id - a_ra * a_id
        for ( a_eq_dash, a_Xd_dash, a_ra, a_id, a_iq ) in
            zip(gens_eq_dash,
                gens_Xd_dash,
                gens_ra,
                gens_i_d_0,
                gens_i_q_0 ) ]    
    
    #-----------------------------------------------

    gens_ph = [ a_vd * a_id + a_vq * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    gens_qh = [ a_vq * a_id - a_vd * a_iq
                for ( a_vd, a_vq, a_id, a_iq ) in
                    zip(gens_vd,
                        gens_vq,
                        gens_i_d_0,
                        gens_i_q_0 ) ]

    #-----------------------------------------------

    # a_id * a_vh * sin(a_δ - a_θh) + a_iq * a_vh * cos(a_δ - a_θh)
    
    gens_ph_by_vh_θh_δ_id_iq = [
        get_a_gen_active_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]


    # a_id * a_vh * cos(a_δ - a_θh) - a_iq * a_vh * sin(a_δ - a_θh)
    
    gens_qh_by_vh_θh_δ_id_iq = [
        get_a_gen_reactive_power_injection_by_vh_θh_id_iq(
            a_vh, a_θh, a_δ, a_id, a_iq )
                for ( a_vh, a_θh, a_δ, a_id, a_iq ) in
                    zip(gens_vh,
                        gens_θh,
                        gens_δ,
                        gens_i_d_0,
                        gens_i_q_0 ) ]
    
    #-----------------------------------------------
    
    S_gens = [
        a_vh * exp(im * (a_θh - a_δ + π/2)) *
            ( a_id - im * a_iq )
        for ( a_vh, a_θh, a_δ, a_id, a_iq) in
            zip(gens_vh,
                gens_θh,
                gens_δ,
                gens_i_d_0,
                gens_i_q_0) ]
    
    #-----------------------------------------------
    
    gens_Pei = [
        ed_dash * i_d_0 + eq_dash * i_q_0 +
            (Xq_dash - Xd_dash ) * i_d_0 *  i_q_0
        for (ed_dash,eq_dash,Xd_dash,Xq_dash,i_d_0,i_q_0 ) in
            zip(gens_ed_dash, gens_eq_dash,
                gens_Xd_dash,
                gens_Xq_dash,
                gens_i_d_0, gens_i_q_0 ) ]

    #----------------------------------------    
    #----------------------------------------
    # init dyn_pf_para
    #----------------------------------------
    #----------------------------------------    

    
    if loc_load_exist == true

        dyn_pf_fun_flat_para =
            [gens_ph...;
             gens_qh...;
             P_non_gens...;
             Q_non_gens...;
             P_g_loc_load...;
             Q_g_loc_load...]
        
        init_dyn_pf_flat_para  =
            dyn_pf_fun_flat_para
        
    else

        dyn_pf_fun_flat_para =
            [gens_ph...;
             gens_qh...;
             P_non_gens...;
             Q_non_gens... ]
        
        init_dyn_pf_flat_para  =
            dyn_pf_fun_flat_para
        
    end
    
    #-----------------------------------------------

    dyn_pf_fun_kwd_net_para = 
        (; 
         Ynet,
         nodes_node_idx_and_incident_edges_other_node_idx )

    #-----------------------------------------------    
    #-----------------------------------------------
    #-----------------------------------------------

    dyn_pf_Idxs_kwd_para =
        (;
         dyn_pf_P_Q_δ_etc_kwd_para_Idxs,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_net_para,
         δ_ed_dash_eq_dash_Idxs_in_flattend )
    
    #-----------------------------------------------    
    
    vh_θh_id_iq =
        [ vh; θh; gens_i_d_0; gens_i_q_0 ]

    
    #-----------------------------------------------

    flat_ωs_ωref0_vref0_porder0 =
        f_gens_ωs_ωref0_vref0_porder0 =
        [ gens_nodes_ωs_ωref0_vref0_porder0...; ]
    
    #-----------------------------------------------

    flat_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para =
        f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para =
        [f_gens_ωs_ωref0_vref0_porder0;
         dyn_pf_fun_flat_para ]
    
    #-----------------------------------------------    
    #-----------------------------------------------

    intra_dyn_pf_mismatch_kwd_para =
        (;
         loc_load_exist,
         dyn_pf_vars_Idxs,
         dyn_pf_Idxs_kwd_para,
         gens_nodes_ra_Xd_dash_Xq_dash)
    
    #-----------------------------------------------

    intra_pf_kwd_para =
        (;          
         intra_flat_ur_ui_Idx,

         intra_flat_vh_θh_Idx,
     
         intra_dyn_pf_flat_para_Idx,

         intra_dyn_pf_mismatch_kwd_para,
         
         ur_ui_idx_in_Idx,
          
         nodes_u_Idx_in_ranges,

         nodes_δ_ed_dash_eq_dash_Idxs )    
    
    #-----------------------------------------------

    ode_per_gen_models_func_kwd_paras =  [
        (;
         Ax_view,
         Bx_view,
         Cx_view,
         per_gen_ωs_ωref0_vref0_porder0_Idx,
         per_gen_id_iq_pg_vh_Idx )
        for (Ax_view, Bx_view, Cx_view) in
            zip(
              vec_Ax_views,
              vec_Bx_views,
              vec_Cx_views, ) ]
        
    #-----------------------------------------------


    ode_per_para_model_func_kwd_para =
        (;
         Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         per_para_ωs_ωref0_vref0_porder0_Idx,
         per_para_id_iq_pg_vh_Idx )
         
    #-----------------------------------------------
    
    vtf_gens_fun_kwd_para = [
        (;
         gen_ra_Xd_dash_Xq_dash,
         a_gen_vtf_vh_θh_Idx,
         a_gen_vtf_δ_ed_dash_eq_dash_Idx )

        for gen_ra_Xd_dash_Xq_dash in
            gens_nodes_ra_Xd_dash_Xq_dash ]


    vtf_para_and_idxs =
        (;
         vtf_gens_fun_kwd_para,
         
         gens_nodes_u_Idx_in_ranges,
         
         non_gens_nodes_u_Idx_in_ranges,
         )

    #-----------------------------------------------

    disaggretation_idxs =
        (; 
         gens_ωs_ωref0_vref0_porder0_idx_in_Idx,
         
         
         a_gen_vtf_vh_θh_Idx,
         
         a_gen_vtf_δ_ed_dash_eq_dash_Idx,
         
         
         per_gen_ωs_ωref0_vref0_porder0_Idx,
         
         per_gen_id_iq_pg_vh_Idx,
         
         
         per_para_ωs_ωref0_vref0_porder0_Idx,
         
         per_para_id_iq_pg_vh_Idx,
         
         
         f_ωs_ωref0_vref0_porder0_Idx,
         
         f_dyn_pf_para_Idx,

         f_dyn_ode_pf_para_Idx,
         
         
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs
         )

    #-----------------------------------------------

    post_pf_idxs =
        (;
         vh_Idx,
         θh_Idx,
         id_Idx,
         iq_Idx )

    #-----------------------------------------------


    states_and_matrices =
        (;
         im_state,

         state,

         im_mass_matrix,
         im_ode_mass_matrix,
         im_model_ode_mass_matrix,
         im_model_pf_mass_matrix  )


    labels_and_symbols =
        (;
         net_class_names,
         im_net_states_and_var_labels,
         gens_nodes_im_vars_labels,
         net_bus_volts_labels,
         im_sym,
         im_ode_sym )

    
    #-----------------------------------------------
        
    dynamics_by_gens_kwd_para =
        (;

         loc_load_exist,

         gens_nodes_idx,

         non_gens_nodes_idx,

         each_gens_im_vars_Idx_in_state,

         nodes_state_Idx,

         im_vars_Idx_in_state,

         nodes_ur_ui_Idx_in_state,

         nodes_u_Idx_in_ranges,

         gens_nodes_ra_Xd_dash_Xq_dash,

         nodes_δ_ed_dash_eq_dash_Idxs,

         post_pf_idxs,

         disaggretation_idxs,

         intra_pf_kwd_para,

         vtf_para_and_idxs,

         gens_Ax_update_parameters,
         gens_nodes_collection,

         vec_Ax_views,

         ode_per_para_model_func_kwd_para,

         )

    #-----------------------------------------------

    dynamics_by_gens_para =
        (;
         dynamics_by_gens_kwd_para,
         f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
         states_and_matrices,
         labels_and_symbols )

    #-----------------------------------------------
        
    dynamics_by_per_gen_kwd_para  =
        (;
         loc_load_exist,

         gens_nodes_idx,

         non_gens_nodes_idx,

         each_gens_im_vars_Idx_in_state,

         nodes_state_Idx,

         im_vars_Idx_in_state,

         nodes_ur_ui_Idx_in_state,

         nodes_u_Idx_in_ranges,

         gens_nodes_ra_Xd_dash_Xq_dash,

         nodes_δ_ed_dash_eq_dash_Idxs,

         post_pf_idxs,

         disaggretation_idxs,

         intra_pf_kwd_para,

         vtf_para_and_idxs,

         gens_Ax_update_parameters,
         gens_nodes_collection,

         vec_Ax_views,

         ode_per_gen_models_func_kwd_paras     

         )

    #-----------------------------------------------

    dynamics_by_per_gen_para =
        (;
         dynamics_by_per_gen_kwd_para,
         f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,
         states_and_matrices,
         labels_and_symbols )

    #-----------------------------------------------
    #-----------------------------------------------

    if all_wt_size_and_dims == false

        Idx_sets =
            (;
             Idx_set_A,
             Idx_set_B,
             Idx_set_C,
             Idx_set_misc,
             sizes_and_dims )
        
        para_and_cal =
            (;
             
             gens_Ax_update_parameters,
             im_plants_system_matrices,
             integrated_pf_and_init,
             named_tup_pf_result,
             state,
             sim_state_x0,
             gens_sim_state_x0,
             
             vh, θh,             
             
             gens_nodes_ωs_ωref0_vref0_porder0,
             gens_nodes_τm_vf,
             gens_dynamic_id_iq_pg_vh,             
             gens_nodes_δ_ed_dash_eq_dash,
             gens_nodes_ra_Xd_dash_Xq_dash,

             gens_i_d_0 ,
             gens_i_q_0,
             
             gens_vd, gens_vq,
             gens_ph, gens_qh,
             
             gens_ph_by_vh_θh_δ_id_iq,
             gens_qh_by_vh_θh_δ_id_iq,
             
             S_gens, gens_Pei,

             loc_load_exist,

             dyn_pf_fun_flat_para,
             dyn_pf_fun_kwd_net_para,

             Ynet,
             nodes_node_idx_and_incident_edges_other_node_idx,

             dyn_pf_Idxs_kwd_para,
             vh_θh_id_iq,

             f_gens_ωs_ωref0_vref0_porder0,
             f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,

             intra_dyn_pf_mismatch_kwd_para,
             intra_pf_kwd_para,

             dynamics_by_gens_kwd_para,
             dynamics_by_per_gen_kwd_para,
             dynamics_by_gens_para,
             dynamics_by_per_gen_para )


        if streamlined == true

            return (;
                    loc_load_exist,
                    gens_Ax_update_parameters,
                    im_vars_Idx_in_state,
                    nodes_ur_ui_Idx_in_state,
                    nodes_u_Idx,
                    ur_ui_Idx_in_state,
                    
                    gens_nodes_u_Idx_in_ranges,
                    non_gens_nodes_u_Idx_in_ranges,

                    flat_ur_flat_ui_Idx,

                    state,

                    vec_Ax_views,
                    
                    each_gens_im_vars_Idx_in_state,
                    nodes_state_Idx,
                    nodes_δ_ed_dash_eq_dash_Idxs,
                    
                    gens_nodes_collection,
                    
                    f_dyn_ode_pf_para_Idx,

                    dyn_pf_P_Q_δ_etc_kwd_para_Idxs,

                    Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
                    f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,

                    ur_ui_idx_in_Idx,
                    vh_θh_idx_in_Idx,
                    
                    nodes_u_Idx_in_ranges,
                    
                    intra_flat_ur_ui_Idx,
                    intra_flat_vh_θh_Idx, 
                    intra_dyn_pf_flat_para_Idx,

                    im_mass_matrix,
                    im_ode_mass_matrix,
                    im_model_pf_mass_matrix,
                    labels_and_symbols,
                    gens_nodes_im_vars_labels, 
                    net_bus_volts_labels,
                    im_net_states_and_var_labels,
                    im_sym,
                    im_ode_sym,

                    #
                    post_pf_idxs,
                    intra_pf_kwd_para,
                    intra_dyn_pf_mismatch_kwd_para,

                    dynamics_by_gens_kwd_para,
                    dynamics_by_per_gen_kwd_para,
                    states_and_matrices,
                    δ_ed_dash_eq_dash_Idxs_in_flattend,

                    
                    flat_δ_ed_dash_eq_dash_Idxs_in_state,
                    #

                    dyn_pf_vars_Idxs,
                    dyn_pf_Idxs_kwd_para,
                    gens_nodes_ra_Xd_dash_Xq_dash,

                    gens_nodes_idx,
                    non_gens_nodes_idx,

                    vtf_para_and_idxs,
                    vtf_gens_fun_kwd_para,
                    ode_per_gen_models_func_kwd_paras,

                    init_dyn_pf_flat_para,
                    flat_ωs_ωref0_vref0_porder0,

                    #
                    gens_nodes_τm_vf,
                    gens_dynamic_id_iq_pg_vh,
                    gens_i_d_0, gens_i_q_0,     
                    gens_vd, gens_vq,
                    gens_ph, gens_qh,
                    #
                    
                    gens_ωs_ωref0_vref0_porder0_idx_in_Idx,

                    #

                    vec_im_pure_states_Idxs,
                    im_pure_states_idxs,
                    vec_im_τm_tilade_vf_states_Idxs,
                    im_τm_tilade_vf_states_idxs,
                    vec_im_sauer_states_Idxs,
                    flat_im_sauer_states_Idxs,
                    vec_im_sauer_δ_eq_dash_ed_dash_Idxs,
                    flat_im_sauer_δ_eq_dash_ed_dash_Idxs,
                    
                    #

                    Idx_sets
                    # para_and_cal
                    )
        else

            return (; Idx_sets, para_and_cal )

        end

        
    else
       
        (;
         Idx_set_A,
         Idx_set_B,
         Idx_set_C,
         Idx_set_misc,
         all_idxs_and_Idxs,
         sizes_and_dims ) =
            dynamic_model_and_powerflow_Idxs_set
        
        Idx_sets = (; Idx_set_A, Idx_set_B, Idx_set_C, Idx_set_misc, sizes_and_dims )

        
        para_and_cal =
            (;             
             gens_Ax_update_parameters,
             im_plants_system_matrices,
             integrated_pf_and_init,
             named_tup_pf_result,
             state,
             sim_state_x0,
             gens_sim_state_x0,
             
             vh, θh,             
             
             gens_nodes_ωs_ωref0_vref0_porder0,
             gens_nodes_τm_vf,
             gens_dynamic_id_iq_pg_vh,             
             gens_nodes_δ_ed_dash_eq_dash,
             gens_nodes_ra_Xd_dash_Xq_dash,

             gens_i_d_0,
             gens_i_q_0,
             
             gens_vd, gens_vq,
             gens_ph, gens_qh,
             
             gens_ph_by_vh_θh_δ_id_iq,
             gens_qh_by_vh_θh_δ_id_iq,
             
             S_gens, gens_Pei,

             loc_load_exist,

             dyn_pf_fun_flat_para,
             dyn_pf_fun_kwd_net_para,

             Ynet,
             nodes_node_idx_and_incident_edges_other_node_idx,

             dyn_pf_Idxs_kwd_para,
             vh_θh_id_iq,

             f_gens_ωs_ωref0_vref0_porder0,
             f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,

             intra_dyn_pf_mismatch_kwd_para,
             intra_pf_kwd_para,

             dynamics_by_gens_kwd_para,
             dynamics_by_per_gen_kwd_para,
             dynamics_by_gens_para,
             dynamics_by_per_gen_para )
        
        all_para_and_cal =
            (;
             gens_Ax_update_parameters,
             im_plants_system_matrices,
             integrated_pf_and_init,
             named_tup_pf_result,
             state,
             sim_state_x0,
             gens_sim_state_x0,
             vh, θh,
             gens_vh, gens_θh, gens_vh_θh,
             
             non_gens_vh, non_gens_θh, non_gens_vh_θh,
             gens_nodes_ωs_ωref0_vref0_porder0,
             gens_nodes_τm_vf,
             
             gens_nodes_δ_ed_dash_eq_dash,

             gens_nodes_ra_Xd_dash_Xq_dash,

             gens_δ, gens_ed_dash, gens_eq_dash,

             gens_ra, gens_Xd_dash, gens_Xq_dash,

             gens_dynamic_id_iq_pg_vh,

             gens_id_iq, gens_i_d_0 ,  gens_i_q_0,

             gens_vd, gens_vq, gens_ph, gens_qh,

             gens_ph_by_vh_θh_δ_id_iq,
             gens_qh_by_vh_θh_δ_id_iq,

             S_gens, gens_Pei,

             loc_load_exist,

             dyn_pf_fun_flat_para,
             dyn_pf_fun_kwd_net_para,

             Ynet,
             nodes_node_idx_and_incident_edges_other_node_idx,

             dyn_pf_Idxs_kwd_para,
             vh_θh_id_iq,

             f_gens_ωs_ωref0_vref0_porder0 ,
             f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,

             intra_dyn_pf_mismatch_kwd_para,
             intra_pf_kwd_para,
             
             ode_per_gen_models_func_kwd_paras,
             ode_per_para_model_func_kwd_para,

             vtf_gens_fun_kwd_para,
             vtf_para_and_idxs,

             disaggretation_idxs,

             post_pf_idxs,

             states_and_matrices,

             labels_and_symbols,

             dynamics_by_gens_kwd_para,
             dynamics_by_per_gen_kwd_para,
             dynamics_by_gens_para,
             dynamics_by_per_gen_para )
        
        return  (; Idx_sets,
                 para_and_cal,
                 all_idxs_and_Idxs,
                 all_para_and_cal,
                 sizes_and_dims )
    end



end



function driver_get_dynamics_paras_and_Idxs_set()
    
    pf_alg  = NewtonRaphson()
    
    # ode_alg       = Rodas4()    
    # algr_name     = "rodas4" 
    
    ode_alg         = ImplicitMidpoint()    
    algr_name       = "ImplicitMidpoint"

    # :Idxs_hybrid, :Idxs_im, :Idxs_industrial

    Idxs_type = :Idxs_im 
    
    no_control_device = false
        
    only_gen  = false

    all_wt_size_and_dims = false

    dt              = 0.01

    sim_timespan    = (0.0, 10.0)
    
    abstol          = 1e-14

    reltol          = 1e-14

        
    """
    dynamics_case =
        case_IEEE_14_Bus_dynamic_plant_wt_loc_load_v6_P

    dynamics_case =
        case_IEEE_14_Bus_dynamic_plants_v6_P_rscad

    """
    
    dynamics_case =
        case_IEEE_9_Bus_dynamic_plants_SM_cb_v6_P_sauer
    
    netd =
        NetworkData( dynamics_case()... )
    
    #-----------------------------------------------    
            
    loc_load_exist =
        loc_load_exist_bool( netd )

    #-----------------------------------------------
    
    dynamics_paras_and_Idxs_set =
        get_dynamics_paras_and_Idxs_set(
            netd; Idxs_type = Idxs_type,
            no_control_device = no_control_device,
            only_gen  = only_gen,
            all_wt_size_and_dims = all_wt_size_and_dims,
            pf_alg  = pf_alg )

    
end


# function get_a_im_sim_paras_and_Idxs_set(
#     netd;
#     Idxs_type =
#                 :Idxs_im,
#             no_control_device =
#                 false,
#             only_gen  =
#                 false,
#             all_wt_size_and_dims =
#                 false,
#             pf_alg  =
#                 NewtonRaphson())
#     #-----------------------------------------------
    
#     netd =
#         NetworkData( dynamics_case()... )
    
#     #-----------------------------------------------

#     dynamics_paras_and_Idxs_set =
#         get_dynamics_paras_and_Idxs_set(
#             netd; Idxs_type =
#                 Idxs_type,
#             no_control_device =
#                 no_control_device,
#             only_gen  =
#                 only_gen,
#             all_wt_size_and_dims =
#                 all_wt_size_and_dims,
#             pf_alg  =
#                 pf_alg )

#     (; Idx_sets, para_and_cal ) =
#         dynamics_paras_and_Idxs_set
    
#     #-----------------------------------------------
#     #-----------------------------------------------    
        
#     Idx_set_A =
#         Idx_sets.Idx_set_A

#     Idx_set_B =
#         Idx_sets.Idx_set_B

#     Idx_set_C =
#         Idx_sets.Idx_set_C

#     Idx_set_misc =
#         Idx_sets.Idx_set_misc

#     sizes_and_dims =
#         Idx_sets.sizes_and_dims
#     #-----------------------------------------------
    
#     δ_ed_dash_eq_dash_Idxs_in_flattend =
#         Idx_set_A.δ_ed_dash_eq_dash_Idxs_in_flattend
    
#     #-----------------------------------------------
    
#     Idx_set_B =
#         Idx_sets.Idx_set_B

#     im_vars_Idx_in_state =
#         Idx_set_B.im_vars_Idx_in_state

#     nodes_ur_ui_Idx_in_state =
#         Idx_set_B.nodes_ur_ui_Idx_in_state

#     nodes_u_Idx =
#         Idx_set_B.nodes_u_Idx
   
#     nodes_u_Idx_in_ranges =
#         Idx_set_B.nodes_u_Idx_in_ranges

#     ur_ui_Idx_in_state =
#         Idx_set_B.ur_ui_Idx_in_state

#     gens_nodes_u_Idx_in_ranges =
#         Idx_set_B.gens_nodes_u_Idx_in_ranges

#     non_gens_nodes_u_Idx_in_ranges =
#         Idx_set_B.non_gens_nodes_u_Idx_in_ranges

#     #-----------------------------------------------
    
#     flat_ur_flat_ui_Idx  =
#         Idx_set_misc.flat_ur_flat_ui_Idx
    
#     flat_ur_flat_ui_Idx  =
#         Idx_set_misc.flat_vh_flat_θh_Idx
    
#     #-----------------------------------------------

#     flat_δ_ed_dash_eq_dash_Idxs_in_state =
#         Idx_set_misc.flat_δ_ed_dash_eq_dash_Idxs_in_state
    
#     #-----------------------------------------------    
    
#     para_and_cal =
#         dynamics_paras_and_Idxs_set.para_and_cal

#     #--------------------------------------------------
    
#     loc_load_exist = para_and_cal.loc_load_exist
    
#     #--------------------------------------------------
    
#     gens_Ax_update_parameters =
#         para_and_cal.gens_Ax_update_parameters
    
#     #--------------------------------------------------
    
#     dynamics_by_gens_kwd_para =
#         para_and_cal.dynamics_by_gens_kwd_para
    
#     #--------------------------------------------------

#     state = para_and_cal.state
    
#     #-----------------------------------------------    
    
#     dynamics_by_per_gen_para =
#         para_and_cal.dynamics_by_per_gen_para

    
#     dynamics_by_per_gen_kwd_para =
#         dynamics_by_per_gen_para.dynamics_by_per_gen_kwd_para

#     vec_Ax_views =
#         dynamics_by_per_gen_kwd_para.vec_Ax_views
    
#     each_gens_im_vars_Idx_in_state =
#         dynamics_by_per_gen_kwd_para.each_gens_im_vars_Idx_in_state

#     nodes_state_Idx =
#         dynamics_by_per_gen_kwd_para.nodes_state_Idx


#     nodes_δ_ed_dash_eq_dash_Idxs =
#         dynamics_by_per_gen_kwd_para.nodes_δ_ed_dash_eq_dash_Idxs

    
#     gens_nodes_collection =
#         dynamics_by_per_gen_kwd_para.gens_nodes_collection
    
#     #-----------------------------------------------    
    
#     disaggretation_idxs =
#         dynamics_by_per_gen_kwd_para.disaggretation_idxs

#     f_dyn_ode_pf_para_Idx =
#         disaggretation_idxs.f_dyn_ode_pf_para_Idx

#     Pg_Qg_Png_Qng_Pgll_Qgll_Idxs =
#         disaggretation_idxs.Pg_Qg_Png_Qng_Pgll_Qgll_Idxs

#     f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para =
#         dynamics_by_per_gen_para.f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para


#     ###
#     post_pf_idxs =
#         dynamics_by_per_gen_kwd_para.post_pf_idxs
    
#     #----------------------------------------        
    
#     intra_pf_kwd_para =
#         dynamics_by_per_gen_kwd_para.intra_pf_kwd_para
    
#     ur_ui_idx_in_Idx =
#         intra_pf_kwd_para.ur_ui_idx_in_Idx

#     vh_θh_idx_in_Idx = ur_ui_idx_in_Idx
    
#     nodes_u_Idx_in_ranges =
#         intra_pf_kwd_para.nodes_u_Idx_in_ranges
    
#     intra_flat_ur_ui_Idx =
#         intra_pf_kwd_para.intra_flat_ur_ui_Idx
    
#     intra_flat_vh_θh_Idx =
#         intra_pf_kwd_para.intra_flat_vh_θh_Idx
    
#     intra_dyn_pf_flat_para_Idx =
#         intra_pf_kwd_para.intra_dyn_pf_flat_para_Idx

#     # 

#     intra_dyn_pf_mismatch_kwd_para =
#         intra_pf_kwd_para.intra_dyn_pf_mismatch_kwd_para
    
#     #----------------------------------------    
    
#     states_and_matrices =
#         dynamics_by_per_gen_para.states_and_matrices

#     im_mass_matrix =
#         states_and_matrices.im_mass_matrix
    
#     im_ode_mass_matrix =
#         states_and_matrices.im_ode_mass_matrix

#     im_model_pf_mass_matrix =
#         states_and_matrices.im_model_pf_mass_matrix
    
#     #----------------------------------------    

#     labels_and_symbols =
#         dynamics_by_per_gen_para.labels_and_symbols
    
#     gens_nodes_im_vars_labels =
#         labels_and_symbols.gens_nodes_im_vars_labels

#     net_bus_volts_labels =
#         labels_and_symbols.net_bus_volts_labels

#     im_net_states_and_var_labels =
#         labels_and_symbols.im_net_states_and_var_labels
    
#     im_sym =
#         labels_and_symbols.im_sym
    
#     im_ode_sym =
#         labels_and_symbols.im_ode_sym


#     return (;
#             loc_load_exist,
#             gens_Ax_update_parameters,
#             im_vars_Idx_in_state,
#             nodes_ur_ui_Idx_in_state,
#             nodes_u_Idx,
#             ur_ui_Idx_in_state,
#             gens_nodes_u_Idx_in_ranges,
#             non_gens_nodes_u_Idx_in_ranges,

#             flat_ur_flat_ui_Idx,

#             state,

#             vec_Ax_views,
#             each_gens_im_vars_Idx_in_state,
#             nodes_state_Idx,
#             nodes_δ_ed_dash_eq_dash_Idxs,

#             gens_nodes_collection,
#             f_dyn_ode_pf_para_Idx,

#             Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
#             f_gens_ωs_ωref0_vref0_porder0_wt_dyn_pf_para,

#             ur_ui_idx_in_Idx,
            
#             nodes_u_Idx_in_ranges,
#             intra_flat_ur_ui_Idx,
#             intra_flat_vh_θh_Idx, 
#             intra_dyn_pf_flat_para_Idx,

#             im_mass_matrix,
#             im_ode_mass_matrix,
#             im_model_pf_mass_matrix,
#             labels_and_symbols,
#             gens_nodes_im_vars_labels, 
#             net_bus_volts_labels,
#             im_net_states_and_var_labels,
#             im_sym,
#             im_ode_sym,

#             #
#             post_pf_idxs,
#             intra_pf_kwd_para,
#             intra_dyn_pf_mismatch_kwd_para,

#             dynamics_by_gens_kwd_para,
#             dynamics_by_per_gen_kwd_para,
#             states_and_matrices,
#             δ_ed_dash_eq_dash_Idxs_in_flattend,

#             flat_δ_ed_dash_eq_dash_Idxs_in_state,
#             #

#             Idx_sets,
#             para_and_cal )   

# end

