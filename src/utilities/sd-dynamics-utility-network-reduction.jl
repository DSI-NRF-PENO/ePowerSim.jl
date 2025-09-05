# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247


#####################################################
# ---------------------------------------------------
# Network reduction utility functions
# ---------------------------------------------------
#####################################################

# ---------------------------------------------------
# ---------------------------------------------------
# Kron reduction
# ---------------------------------------------------
# ---------------------------------------------------

# https://delabays.xyz

function kron(A, eliminated)
    
    row_size, col_size = size(A)

    retained = collect(1:row_size)[eliminated]

    A_rr = A[retained, retained]
    A_re = A[retained, eliminated]
    A_er = A[eliminated, retained]
    A_ee = A[eliminated, eliminated]    

    return A_rr - A_re * (A_ee \ A_er)
    
end


function kron(A, eliminated, retained)
    
    row_size, col_size = size(A)

    # retained = collect(1:row_size)[eliminated]

    A_rr = A[retained, retained]
    A_re = A[retained, eliminated]
    A_er = A[eliminated, retained]
    A_ee = A[eliminated, eliminated]    

    return A_rr - A_re * (A_ee \ A_er)
    
end



function kron_by_idx(A,k)


    row_size, col_size = size(A)
    
    A_new = spzeros(row_size, col_size)

    for i in row_size
        if j != k
            
            for j in col_size
                if i != k
                    A_new[k,j] =
                        A[k,j] - (A[k,n] * A[n,j])/A[n,n]

                end            

            end
        end
        
    end

    return A_new
    
end



function kron_symmetric_by_idx(A,k)

    row_size, col_size = size(A)
    
    A_new = spzeros(row_size, col_size)

    for i in row_size
        if j != k
            
            for j in col_size
                if j != k

                    if i != k
                        A_new[k,j] =
                            xA[k,j] - (A[k,n] * A[n,j])/A[n,n]
                        
                        A_new[j,k] = A_new[k,j]

                    end                     
                end


            end
        end
        
    end

    return A_new
    
end


function get_Ybus_from_Ynet(
    Ynet, nodes_node_idx_and_incident_edges_other_node_idx
    )::Matrix{Number} 

    Ybus_row_size = Ybus_col_size = length( Ynet )

    return reshape(
        vcat([[ a_col in nth_row_idx ?
            Ybus_nth_row[
                findall(x-> x==a_col, nth_row_idx)[1]] :
            0.0  for a_col in 1:Ybus_col_size]
              for (nth_row_idx, Ybus_nth_row) in
                  zip( nodes_node_idx_and_incident_edges_other_node_idx,
                       Ynet )]...;),
        Ybus_row_size,
        Ybus_col_size )
        
    
end


#-----------------------------------------------------


function get_sparseYbus_from_Ynet(
    Ynet, nodes_node_idx_and_incident_edges_other_node_idx
    )::SparseMatrixCSC{Number, Int64}

    Ybus_row_size = Ybus_col_size = length(Ynet )
    
    # https://docs.julialang.org/en/v1/stdlib/SparseArrays/#SparseArrays.sparse
    
    return sparse(reshape(
        vcat([[a_col in nth_row_idx ? Ybus_nth_row[
            findall(x-> x==a_col, nth_row_idx)[1]] : 0.0
               for a_col in 1:Ybus_col_size]
              for (nth_row_idx, Ybus_nth_row) in
                  zip(nodes_node_idx_and_incident_edges_other_node_idx,
                      Ynet )]...;),
        Ybus_row_size,
        Ybus_col_size ))
    
    
end


#-------------------------------------------------------

"""
reduced order matrices,

"""

#-------------------------------------------------------

function get_Y_aug_matrices(
    pf_PQ_param,
    vh;
    y_aug_kw_para =
        y_aug_kw_para  )

    #-------------------------------
    
    X_d_dash = y_aug_kw_para.X_d_dash
    
    pf_kw_para = y_aug_kw_para.pf_kw_para

    #-------------------------------
    
    # loc_load_exist =
    #     pf_kw_para.loc_load_exist
    
    # pf_kw_net_para =
    #     pf_kw_para.pf_kw_net_para
    
    # pf_kw_PQ_para_idxs =
    #     pf_kw_para.pf_kw_PQ_para_idxs
    
    # nodes_types_idxs =
    #     pf_kw_para.pf_kw_nodes_types_idxs
    
    # n2s_idxs =
    #     pf_kw_para.pf_kw_n2s_idxs

    (;loc_load_exist,
     pf_kw_net_para,
     pf_kw_PQ_para_idxs,
     pf_kw_nodes_types_idxs,
     pf_kw_n2s_idxs) =
        NamedTupleTools.select(
            pf_kw_para,
            (:loc_load_exist,
             :pf_kw_net_para,
             :pf_kw_PQ_para_idxs,
             :pf_kw_nodes_types_idxs,
             :pf_kw_n2s_idxs) )


    (;Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            pf_kw_net_para,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx))
    
     (;P_gens_sta_para_Idxs,
      Q_gens_sta_para_Idxs,
      P_non_gens_sta_para_Idxs,
      Q_non_gens_sta_para_Idxs,
      P_g_loc_load_sta_para_Idxs,
      Q_g_loc_load_sta_para_Idxs ) =
          NamedTupleTools.select(
              pf_kw_PQ_para_idxs,
              (:P_gens_sta_para_Idxs,
               :Q_gens_sta_para_Idxs,
               :P_non_gens_sta_para_Idxs,
               :Q_non_gens_sta_para_Idxs,
               :P_g_loc_load_sta_para_Idxs,
               :Q_g_loc_load_sta_para_Idxs ) )
    
    
   (;slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_with_loc_load_idx,
    all_nodes_idx) =
        NamedTupleTools.select(
            nodes_types_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_with_loc_load_idx,
             :all_nodes_idx ) ) 

   (;n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx ))
    
    #-------------------------------
    
    dim_gens_nodes =
        length( gens_nodes_idx )

    dim_non_gens_nodes =
        length( non_gens_nodes_idx )
    
    #-------------------------------

    P_non_gens  =
        pf_PQ_param[ P_non_gens_sta_para_Idxs ]

    Q_non_gens = 
        pf_PQ_param[ Q_non_gens_sta_para_Idxs ] 
    
    if loc_load_exist == true

        P_g_loc_load =
            pf_PQ_param[
                P_g_loc_load_sta_para_Idxs ] 

        Q_g_loc_load =
            pf_PQ_param[
                Q_g_loc_load_sta_para_Idxs ] 
        
    else

        P_g_loc_load = [ 0.0]
            
        Q_g_loc_load = [ 0.0]
    end

    #-------------------------------

    # # y_diag_x_d_dash = Diag(1/(im * X_d_dash )j)
    
    # y_x_d_dash  = [ 1/(im * x_d_dash)
    #                    for x_d_dash in X_d_dash]

    # # Y_N2 = Ynet + Diag( y_diag_x_d_dash )
    
    # Ynet_aug_yd_dash = [ idx ∈ gens_nodes_idx ?
    #     [ row_item_idx == 1 ?
    #     Ynet[idx][row_item_idx] +
    #     y_x_d_dash[n2s_gens_idx[idx]] :
    #     Ynet[idx][row_item_idx]
    #       for row_item_idx in 1:length(Ynet[idx]) ] :
    #           Ynet[idx] for idx in all_nodes_idx  ]

    # #-------------------------------
    
    # # y_Li eq 7.205 (sauer)
    
    # y_i_PQ = [idx ∈ non_gens_nodes_idx ?
    #     -((P_non_gens[n2s_non_gens_idx[idx]] -
    #     im * Q_non_gens[
    #         n2s_non_gens_idx[idx]])/vh[
    #             n2s_all_nodes_idx[idx]]^2) :
    #             idx ∈ gens_with_loc_load_idx ? -(
    #                 (P_g_loc_load[
    #                     n2s_gens_with_loc_load_idxs[idx]] -
    #                         im * Q_g_loc_load[
    #                             n2s_gens_with_loc_load_idxs[
    #                                 idx]])/ vh[
    #                                     n2s_all_nodes_idx[
    #                                         idx]]^2) :
    #                                     0.0 + im * 0.0
    #         for idx in all_nodes_idx ]

    # # Y_N2 = Y_N1 + Diag(y_Li)
    
    # Ynet_aug_yd_dash_wt_y_i_PQ = [idx ∈ non_gens_nodes_idx ?
    #     [row_item_idx == 1 ?
    #     Ynet_aug_yd_dash[idx][row_item_idx] + y_i_PQ[idx] :
    #     Ynet_aug_yd_dash[idx][row_item_idx]
    #      for row_item_idx in
    #          1:length(Ynet_aug_yd_dash[idx])] :
    #              idx ∈ gens_with_loc_load_idx ?
    #              [ row_item_idx == 1 ?
    #              Ynet_aug_yd_dash[idx][row_item_idx] +
    #              y_i_PQ[idx] :
    #              Ynet_aug_yd_dash[idx][row_item_idx]
    #                for row_item_idx in
    #                    1:length(Ynet_aug_yd_dash[idx])] :
    #                        Ynet_aug_yd_dash[idx]
    #                               for idx in all_nodes_idx  ]
    
    # #-------------------------------
    
    # y_x_d_dash_diag = Diagonal( y_x_d_dash  )

    # # Ybus = get_sparseYbus_from_Ynet(
    # #     Ynet,
    # #     nodes_idx_with_adjacent_nodes_idx)

    # # Y_aug = Yii_diag + Ybus

    # # get_Ybus_from_Ynet
    
    # # get_sparseYbus_from_Ynet

    # YBus_net_aug_yd_dash_wt_y_i_PQ = get_Ybus_from_Ynet(
    #         Ynet_aug_yd_dash_wt_y_i_PQ,
    #         nodes_idx_with_adjacent_nodes_idx)

    # Y_A = y_x_d_dash_diag
    
    # Y_B = hcat(-y_x_d_dash_diag,
    #        zeros(dim_gens_nodes,
    #              dim_non_gens_nodes))
        
    # Y_C = vcat(-y_x_d_dash_diag,
    #            zeros( dim_non_gens_nodes,
    #                   dim_gens_nodes))

    # Y_D = YBus_net_aug_yd_dash_wt_y_i_PQ 

    # # Y_D = im * imag.( YBus_net_aug_yd_dash_wt_y_i_PQ )

    # Y_internal_nodes = Y_A - Y_B * inv(Y_D) * Y_C

    # return (; Ynet,
    #         y_x_d_dash,
    #         Ynet_aug_yd_dash,
    #         y_i_PQ,
    #         Ynet_aug_yd_dash_wt_y_i_PQ,
    #         y_x_d_dash_diag,
    #         Y_internal_nodes,
    #         Y_A,
    #         Y_B,
    #         Y_C,
    #         Y_D )

    # --------------------------------
    
    # Alternative

    YBus =  get_Ybus_from_Ynet(
        Ynet,
        nodes_idx_with_adjacent_nodes_idx) # .* 1/baseMVA 
    
    # y_Li eq 7.205 (sauer)
    
    y_i_PQ = [idx ∈ non_gens_nodes_idx ?
        -((P_non_gens[n2s_non_gens_idx[idx]] -
        im * Q_non_gens[
            n2s_non_gens_idx[idx]])/vh[
                n2s_all_nodes_idx[idx]]^2) :
                idx ∈ gens_with_loc_load_idx ? -(
                    (P_g_loc_load[
                        n2s_gens_with_loc_load_idxs[idx]] -
                            im * Q_g_loc_load[
                                n2s_gens_with_loc_load_idxs[
                                    idx]])/ vh[
                                        n2s_all_nodes_idx[
                                            idx]]^2) :
                                        0.0 + im * 0.0
            for idx in all_nodes_idx ]

    # --------------------------------
    
    Y_dash = YBus +  Diagonal( -y_i_PQ )

    Y_1 = Y_dash[
        gens_nodes_idx,
        gens_nodes_idx]

    Y_2 = Y_dash[
        gens_nodes_idx,
        non_gens_nodes_idx]

    Y_3 = Y_dash[
        non_gens_nodes_idx,
        gens_nodes_idx]

    Y_4 = Y_dash[
        non_gens_nodes_idx,
        non_gens_nodes_idx]

    # Yred = Y_1 - Y_2 * inv(Y_4) * Y_3

    Yred = Y_1 - Y_2 * (Y_4 \ Y_3)

    # --------------------------------
    # --------------------------------
    
    # y_diag_x_d_dash = Diag(1/(im * X_d_dash )j)
    
    y_x_d_dash  = [ 1/(im * x_d_dash)
                    for x_d_dash in X_d_dash]
    
    diag_y_x_d_dash_wt_aug_zero  = Diagonal(
        [ idx ∈ gens_nodes_idx ?
            1/(im * X_d_dash[n2s_gens_idx[idx] ]) :
            0.0 - im * 0.0 for idx in all_nodes_idx ])

    # --------------------------------
    
    Y_N2 = YBus + diag_y_x_d_dash_wt_aug_zero +
        Diagonal(-y_i_PQ )

    # --------------------------------

    y_x_d_dash_diag = Diagonal(y_x_d_dash)

    Y_A = y_x_d_dash_diag
    
    Y_B = hcat(-y_x_d_dash_diag,
           zeros(dim_gens_nodes,
                 dim_non_gens_nodes))
        
    Y_C = vcat(-y_x_d_dash_diag,
               zeros( dim_non_gens_nodes,
                      dim_gens_nodes))

    Y_D =  vcat(hcat(Y_N2[gens_nodes_idx, gens_nodes_idx],
                    Y_N2[gens_nodes_idx, non_gens_nodes_idx]),
               hcat(Y_N2[non_gens_nodes_idx, gens_nodes_idx],
                    Y_N2[non_gens_nodes_idx,
                         non_gens_nodes_idx]))

    # --------------------------------
    
    # Y_D = im * imag.( YBus_net_aug_yd_dash_wt_y_i_PQ )

    # Y_internal_nodes = Y_A - Y_B * inv(Y_D) * Y_C

    Y_internal_nodes = Y_A - Y_B * (Y_D \ Y_C)
    

    return (; Ynet,
            YBus,
            Yred,
            y_i_PQ,
            Y_N2,
            y_x_d_dash,                        
            Y_A,
            Y_B,
            Y_C,
            Y_D,
            Y_internal_nodes)
        
end



function get_Y_aug_matrices(
    P_non_gens,
    Q_non_gens,
    P_g_loc_load,
    Q_g_loc_load,
    vh,
    Ynet_wt_nodes_idx_wt_adjacent_nodes;
    y_aug_kw_para =
        y_aug_kw_para  )
        
    #-------------------------------

    (;loc_load_exist,
     X_d_dash,    
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs) =
         NamedTupleTools.select(
             y_aug_kw_para,
             (:loc_load_exist,
              :X_d_dash,    
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs)) 

    
    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     # gens_nodes_with_loc_loads_idx,
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
    

    (Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx)) 
    
    YBus =  get_Ybus_from_Ynet(
        Ynet,
        nodes_idx_with_adjacent_nodes_idx) # .* 1/baseMVA 
    
    #-------------------------------
    
    dim_gens_nodes =
        length( gens_nodes_idx )

    dim_non_gens_nodes =
        length( non_gens_nodes_idx )
    
    #-------------------------------
    
    # y_Li eq 7.205 (sauer)
    
    y_i_PQ = [idx ∈ non_gens_nodes_idx ?
        -((P_non_gens[n2s_non_gens_idx[idx]] -
        im * Q_non_gens[
            n2s_non_gens_idx[idx]])/vh[
                n2s_all_nodes_idx[idx]]^2) :
                idx ∈ gens_with_loc_load_idx ? -(
                    (P_g_loc_load[
                        n2s_gens_with_loc_load_idxs[idx]] -
                            im * Q_g_loc_load[
                                n2s_gens_with_loc_load_idxs[
                                    idx]])/ vh[
                                        n2s_all_nodes_idx[
                                            idx]]^2) :
                                        0.0 + im * 0.0
            for idx in all_nodes_idx ]

    # --------------------------------
    
    Y_dash = YBus +  Diagonal( -y_i_PQ )

    Y_1 = Y_dash[
        gens_nodes_idx,
        gens_nodes_idx]

    Y_2 = Y_dash[
        gens_nodes_idx,
        non_gens_nodes_idx]

    Y_3 = Y_dash[
        non_gens_nodes_idx,
        gens_nodes_idx]

    Y_4 = Y_dash[
        non_gens_nodes_idx,
        non_gens_nodes_idx]

    # Yred = Y_1 - Y_2 * inv(Y_4) * Y_3

    Yred = Y_1 - Y_2 * (Y_4 \ Y_3)

    # --------------------------------
    # --------------------------------
    
    # y_diag_x_d_dash = Diag(1/(im * X_d_dash )j)
    
    y_x_d_dash  = [ 1/(im * x_d_dash)
                    for x_d_dash in X_d_dash]
    
    diag_y_x_d_dash_wt_aug_zero  = Diagonal(
        [ idx ∈ gens_nodes_idx ?
            1/(im * X_d_dash[n2s_gens_idx[idx] ]) :
            0.0 - im * 0.0 for idx in all_nodes_idx ])

    # --------------------------------
    
    Y_N2 = YBus + diag_y_x_d_dash_wt_aug_zero +
        Diagonal(-y_i_PQ )

    # --------------------------------

    y_x_d_dash_diag = Diagonal(y_x_d_dash)

    Y_A = y_x_d_dash_diag
    
    Y_B = hcat(-y_x_d_dash_diag,
           zeros(dim_gens_nodes,
                 dim_non_gens_nodes))
        
    Y_C = vcat(-y_x_d_dash_diag,
               zeros( dim_non_gens_nodes,
                      dim_gens_nodes))

    Y_D =  vcat(hcat(Y_N2[gens_nodes_idx, gens_nodes_idx],
                    Y_N2[gens_nodes_idx, non_gens_nodes_idx]),
               hcat(Y_N2[non_gens_nodes_idx, gens_nodes_idx],
                    Y_N2[non_gens_nodes_idx,
                         non_gens_nodes_idx]))

    # --------------------------------
    
    # Y_D = im * imag.( YBus_net_aug_yd_dash_wt_y_i_PQ )

    # Y_internal_nodes = Y_A - Y_B * inv(Y_D) * Y_C

    Y_internal_nodes = Y_A - Y_B * (Y_D \ Y_C)
    

    return (;Ynet,
            YBus,
            Yred,
            y_i_PQ,
            Y_N2,
            y_x_d_dash,                        
            Y_A,
            Y_B,
            Y_C,
            Y_D,
            Y_internal_nodes)
        
end

# ------------------------------------------------------

function get_Yint_and_Yred_matrices(
    vh,
    Ynet_wt_nodes_idx_wt_adjacent_nodes,
    P_non_gens,
    Q_non_gens,
    P_g_loc_load,
    Q_g_loc_load;
    y_aug_kw_para =
        y_aug_kw_para,
    fault_status =
        :pre_fault_state,
    cleared_selected_lines_faults_net_para =
        nothing )
    
    #-------------------------------

    P_non_gens =
        deepcopy(P_non_gens)
    
    Q_non_gens =
        deepcopy(Q_non_gens)
    
    #-------------------------------

    (;loc_load_exist,
     X_d_dash,    
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs) =
         NamedTupleTools.select(
             y_aug_kw_para,
             (:loc_load_exist,
              :X_d_dash,    
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs)) 

    dyn_pf_fun_kwd_n2s_idxs =
        deepcopy(
            dyn_pf_fun_kwd_n2s_idxs)
    
    dyn_pf_fun_kwd_net_idxs =
        deepcopy(
            dyn_pf_fun_kwd_net_idxs)
    
    
    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_load_idx,
     # gens_nodes_with_loc_loads_idx,
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

    if fault_status == :pre_fault_state

        (Ynet,
        nodes_idx_with_adjacent_nodes_idx) =
            NamedTupleTools.select(
                Ynet_wt_nodes_idx_wt_adjacent_nodes,
                (:Ynet,
                 :nodes_idx_with_adjacent_nodes_idx)) 

        YBus =  get_Ybus_from_Ynet(
            Ynet,
            nodes_idx_with_adjacent_nodes_idx)

        #-------------------------------

        dim_gens_nodes =
            length( gens_nodes_idx )

        dim_non_gens_nodes =
            length( non_gens_nodes_idx )

        non_gens_nodes_wt_f_nodes_idx =
            non_gens_nodes_idx

        #-------------------------------

        # y_Li eq 7.205 (sauer)

        y_i_PQ = [idx ∈ non_gens_nodes_idx ?
            -((P_non_gens[n2s_non_gens_idx[idx]] -
            im * Q_non_gens[
                n2s_non_gens_idx[idx]])/vh[
                    n2s_all_nodes_idx[idx]]^2) :
                    idx ∈ gens_with_loc_load_idx ? -(
                        (P_g_loc_load[
                           n2s_gens_with_loc_load_idxs[idx]] -
                                im * Q_g_loc_load[
                                  n2s_gens_with_loc_load_idxs[
                                        idx]])/ vh[
                                            n2s_all_nodes_idx[
                                                idx]]^2) :
                                            0.0 + im * 0.0
                  for idx in all_nodes_idx ]

    # @show size_y_i_PQ = size(y_i_PQ)
        
    elseif fault_status == :fault_state

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

        #-------------------------------
        
        YBus =  get_Ybus_from_Ynet(
            Ynet,
            nodes_idx_with_adjacent_nodes_idx)

        #-------------------------------

        # non_gens_nodes_idx =
        #     [non_gens_nodes_idx;
        #      fault_nodes_idx]
        
        # n2s_non_gens_idx =
        #     OrderedDict(
        #         node_idx => idx
        #         for (idx, node_idx) in
        #             enumerate(
        #                 non_gens_nodes_idx) )

        # P_non_gens =
        #     [ P_non_gens;
        #       zeros(length(fault_nodes_idx))]
        
        # Q_non_gens =
        #     [Q_non_gens;
        #      zeros(length(fault_nodes_idx))]


        # dim_non_gens_nodes =
        #     length( non_gens_nodes_idx ) 
        
        #-------------------------------
        
        dim_gens_nodes =
            length( gens_nodes_idx )

        dim_non_gens_nodes =
            length( non_gens_nodes_idx ) + 
            length(fault_nodes_idx)

        non_gens_nodes_wt_f_nodes_idx =
            [non_gens_nodes_idx;
             fault_nodes_idx]
        
        #-------------------------------

        # y_Li eq 7.205 (sauer)

        y_i_PQ = [idx ∈ non_gens_nodes_idx ?
            -((P_non_gens[n2s_non_gens_idx[idx]] -
            im * Q_non_gens[
                n2s_non_gens_idx[idx]])/vh[
                    n2s_all_nodes_idx[idx]]^2) :
                    idx ∈ gens_with_loc_load_idx ? -(
                        (P_g_loc_load[
                           n2s_gens_with_loc_load_idxs[idx]] -
                                im * Q_g_loc_load[
                                  n2s_gens_with_loc_load_idxs[
                                        idx]])/ vh[
                                            n2s_all_nodes_idx[
                                                idx]]^2) :
                                                    idx ∈ fault_nodes_idx ?
                                                    0.0 + im * 0.0 :
                                            0.0 + im * 0.0
                for idx in all_nodes_idx ]

        # @show size_y_i_PQ = size(y_i_PQ)
        
        
    elseif fault_status == :post_fault_state


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


        YBus =  get_Ybus_from_Ynet(
            Ynet,
            nodes_idx_with_adjacent_nodes_idx)

        #-------------------------------

        dim_gens_nodes =
            length( gens_nodes_idx )

        dim_non_gens_nodes =
            length( non_gens_nodes_idx ) +
            length(fault_nodes_idx)

        non_gens_nodes_wt_f_nodes_idx =
            [non_gens_nodes_idx;
             fault_nodes_idx]

        #-------------------------------


        # y_Li eq 7.205 (sauer)

        y_i_PQ = [idx ∈ non_gens_nodes_idx ?
            -((P_non_gens[n2s_non_gens_idx[idx]] -
            im * Q_non_gens[
                n2s_non_gens_idx[idx]])/vh[
                    n2s_all_nodes_idx[idx]]^2) :
                    idx ∈ gens_with_loc_load_idx ? -(
                        (P_g_loc_load[
                            n2s_gens_with_loc_load_idxs[idx]] -
                                im * Q_g_loc_load[
                                    n2s_gens_with_loc_load_idxs[
                                        idx]])/ vh[
                                            n2s_all_nodes_idx[
                                                idx]]^2) : idx ∈ fault_nodes_idx ? 0.0 + im * 0.0 :
                                            0.0 + im * 0.0
                for idx in all_nodes_idx ]

        # @show size_y_i_PQ = size(y_i_PQ)
        
        
    else
        throw("fault state not known ")
    end
    
    #-------------------------------
    
    # --------------------------------
    
    Y_dash = YBus +  Diagonal( -y_i_PQ )

    # @show size_Y_dash =  size(Y_dash)
    
    Y_1 = Y_dash[
        gens_nodes_idx,
        gens_nodes_idx]

    Y_2 = Y_dash[
        gens_nodes_idx,
        non_gens_nodes_wt_f_nodes_idx]

    Y_3 = Y_dash[
        non_gens_nodes_wt_f_nodes_idx,
        gens_nodes_idx]

    Y_4 = Y_dash[
        non_gens_nodes_wt_f_nodes_idx,
        non_gens_nodes_wt_f_nodes_idx]

    # Yred = Y_1 - Y_2 * inv(Y_4) * Y_3

    Yred = Y_1 - Y_2 * (Y_4 \ Y_3)

    # @show size_Yred =  size(Yred)
    
    # --------------------------------
        
    y_x_d_dash  = [ 1/(im * x_d_dash)
                    for x_d_dash in X_d_dash]

    # @show size_y_x_d_dash =  size(y_x_d_dash)
               
    diag_y_x_d_dash_wt_aug_zero  = Diagonal(
        [ idx ∈ gens_nodes_idx ?
            1/(im * X_d_dash[n2s_gens_idx[idx] ]) :
            0.0 - im * 0.0 for idx in all_nodes_idx ])

    # @show size_diag_y_x_d_dash_wt_aug_zero =
    #     size(diag_y_x_d_dash_wt_aug_zero)
    
    # --------------------------------
    
    Y_N2 = YBus + diag_y_x_d_dash_wt_aug_zero +
        Diagonal(-y_i_PQ )

    # @show size_Y_N2 =  size(Y_N2)
    # --------------------------------

    y_x_d_dash_diag = Diagonal(y_x_d_dash)

    # @show size_y_x_d_dash_diag = size(y_x_d_dash_diag)
    
    Y_A = y_x_d_dash_diag

    # @show size_Y_A = size(Y_A)
    
    Y_B = hcat(-y_x_d_dash_diag,
           zeros(dim_gens_nodes,
                 dim_non_gens_nodes))

    # @show size_Y_B = size(Y_B)
    
    Y_C = vcat(-y_x_d_dash_diag,
               zeros( dim_non_gens_nodes,
                      dim_gens_nodes))

    # @show size_Y_C = size(Y_C)
    
    Y_D =  vcat(hcat(Y_N2[gens_nodes_idx, gens_nodes_idx],
                     Y_N2[gens_nodes_idx,
                          non_gens_nodes_wt_f_nodes_idx]),
                hcat(Y_N2[non_gens_nodes_wt_f_nodes_idx,
                          gens_nodes_idx],
                     Y_N2[non_gens_nodes_wt_f_nodes_idx,
                          non_gens_nodes_wt_f_nodes_idx]))

    # @show size_Y_D = size(Y_D)

    # --------------------------------
    
    # Y_internal_nodes = Y_A - Y_B * inv(Y_D) * Y_C

    Y_internal_nodes = Y_A - Y_B * (Y_D \ Y_C)

    # @show size(Y_internal_nodes)
    
    return (;Yred,
            Y_internal_nodes)
        
end


function get_net_status_Yint_and_Yred_matrices(
    list_edges_to_have_fault,
    fault_nodes_idx,
    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,    
    generic_system_dynamics_wt_fault_kwd_para )
    
    #----------------------------------------
    # power flow equations
    #----------------------------------------

    (;dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     algebraic_generic_model_wt_fault_sol_kwd_para) =
        NamedTupleTools.select(
            generic_system_dynamics_wt_fault_kwd_para,
            (:dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
             :algebraic_generic_model_wt_fault_sol_kwd_para))
    
    (flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_wt_fault_kwd_para) =
         NamedTupleTools.select(
             algebraic_generic_model_wt_fault_sol_kwd_para,
             (:flat_vh_flat_θh_id_iq_u0,
              :dyn_pf_solver,
              :algebraic_generic_model_wt_fault_kwd_para))

    #----------------------------------------    

    (;loc_load_exist,
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     pf_generic_gens_para,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     cleared_selected_lines_faults_net_para) =
        NamedTupleTools.select(
            algebraic_generic_model_wt_fault_kwd_para,
            (:loc_load_exist,
             :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             :dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx,
             :dyn_pf_fun_kwd_n2s_idxs,
             :dyn_pf_fun_kwd_net_idxs,
             :pf_generic_gens_para,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes,
             :cleared_selected_lines_faults_net_para))


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
    
    
    (dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs) =
         NamedTupleTools.select(
             dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs))

    (
     # dyn_δ_Idxs,
     # dyn_ed_dash_Idxs,       
     # dyn_eq_dash_Idxs,        
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             dyn_δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_Idx ,
             (
              # :dyn_δ_Idx,
              # :dyn_ed_dash_Idx,       
              # :dyn_eq_dash_Idx,                 
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
                 :dyn_Qll_Idx))

    #----------------------------------------    
    
    X_d_dash =
        getproperty(
            pf_generic_gens_para,
            :X_d_dash)
    
    P_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll[
            dyn_P_non_gens_Idxs]
    
    Q_non_gens =
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll[
            dyn_Q_non_gens_Idxs]
    
    if loc_load_exist == true
        P_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll[
                dyn_P_gens_loc_load_Idxs]
        
        Q_g_loc_load =
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll[
                dyn_Q_gens_loc_load_Idxs]
    else
        P_g_loc_load = 0.0
        
        Q_g_loc_load = 0.0
    end
    
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
    #----------------------------------------
    
    list_network_status =
        [:pre_fault_state,
         :fault_state,
         :post_fault_state ]

    dict_status_Yint_Yred =
        Dict{Symbol, NamedTuple }()

    #----------------------------------------
    
    for a_system_status in list_network_status

        if a_system_status == :pre_fault_state

            pf_fun_mismatch =
                pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!

        elseif a_system_status ==  :fault_state 

            pf_fun_mismatch =
                fault_state_algebraic_generic_pf_ΔPQ_mismatch!

        elseif a_system_status ==  :post_fault_state

            pf_fun_mismatch =
                post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!
        else
            throw(" system status not known ")
        end

        #----------------------------------------

        number_of_faults =
            length(list_edges_to_have_fault)
        
        # number_of_faults =
        #     length(fault_nodes_idx) 

        twice_number_of_faults =
            2 * number_of_faults

        fault_vh_sym =
            generate_labels_by_nodes_idxs_and_vars(
                fault_nodes_idx,
                [:vh];
                label_prefix = "bus")

        fault_θh_sym =
            generate_labels_by_nodes_idxs_and_vars(
                fault_nodes_idx,
                [:θh];
                label_prefix = "bus")

        fault_nodes_sym =
            [fault_vh_sym;
             fault_θh_sym ]  

        fault_algeb_init =[
            ones(number_of_faults);
            zeros(number_of_faults)
            ]

        #----------------------------------------

        flat_vh_flat_θh_id_iq_vfh_θfh =
            [flat_vh_flat_θh_id_iq_u0;
             fault_algeb_init]

        #----------------------------------------    
        #----------------------------------------    

        if use_nlsolve == true

            sol = nlsolve((g, x) -> pf_fun_mismatch(
                g, x,
                δ_ed_dash_eq_dash_Png_Qng_Pll_Qll;
                    kwd_para =
                        algebraic_generic_model_wt_fault_kwd_para),
                          flat_vh_flat_θh_id_iq_vfh_θfh,
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
                            algebraic_generic_model_wt_fault_kwd_para)),
                flat_vh_flat_θh_id_iq_vfh_θfh,
                δ_ed_dash_eq_dash_Png_Qng_Pll_Qll ),
                                        pf_alg )


        end

        y_aug_kw_para =
            (;loc_load_exist,
             X_d_dash,    
             dyn_pf_fun_kwd_n2s_idxs,
             dyn_pf_fun_kwd_net_idxs)

        
        (;Yred,  Y_internal_nodes) =
            get_Yint_and_Yred_matrices(
                sol[dyn_pf_vh_Idxs],
                Ynet_wt_nodes_idx_wt_adjacent_nodes,
                P_non_gens,
                Q_non_gens,
                P_g_loc_load,
                Q_g_loc_load;
                y_aug_kw_para =
                    y_aug_kw_para,
                fault_status =
                    a_system_status ,
                cleared_selected_lines_faults_net_para =
                    cleared_selected_lines_faults_net_para)
                        
        dict_status_Yint_Yred[a_system_status] =
            (;Yred = Yred,

             Y_internal_nodes = Y_internal_nodes,

             gens_vh =
                 (sol[dyn_pf_vh_Idxs])[gens_nodes_idx],

             gens_θh =
            (sol[dyn_pf_θh_Idxs])[gens_nodes_idx])
    end
    
    return NamedTupleTools.namedtuple(
        dict_status_Yint_Yred)
    

end


function get_a_status_Yint_Yred(
    ; system_status =
        :pre_fault_state,
    
    case_name = "case9",
    timespan      = 10.0,
    on_fault_time = 9.0,
    clear_fault_time = 9.001,

    with_faults = false,
    
    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [4],

    list_edges_to_have_fault = [ 2 ],
    clear_fault_selection_list = [1],
    
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    json_net_data_by_components_file =
        nothing)

    #--------------------------------------

    """
    case_name = "case9"

    case_name = "case14"

    json_net_data_by_components_file =
        "net-static-data-avr-sauer-gov-sauer.json"

    timespan         = 10.0
    on_fault_time    = 9.0
    clear_fault_time = 9.001

    with_faults      = false
    
    list_fault_point_from_node_a = [0.3]
    list_fault_resistance        = [0.001]
    list_no_line_circuit         =  [4]

    list_edges_to_have_fault   = [ 2 ]
    clear_fault_selection_list = [1]
    
    basekV          = 1.0    
    use_pu_in_PQ    = true
    line_data_in_pu = true
  

    """
    #--------------------------------------
    
    basekV = 1.0

    use_pu_in_PQ = true
    
    line_data_in_pu = true 
    
    #--------------------------------------    

    use_init_u0 = false
    
    use_nlsolve = false
    
    pf_alg        = NewtonRaphson()

    #--------------------------------------    
    
    ode_alg       = Rodas4()
    # ode_alg       = FBDF()
    # ode_alg       = QNDF()
    # ode_alg       = radau()
    # ode_alg       = RadauIIA5()
    # ode_alg       = DFBDF()

    dae_alg       = IDA()
    
    abstol        = 1e-12
    reltol        = 1e-12
    
    dt            = 0.01
    # timespan      = 10.0
    
    tspan         = (0.0, timespan)
    
    sim_timespan  = (0.0, timespan)
    
    plot_timespan = (0.0, timespan)
    
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
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end

    case_data_dir =
        joinpath(data_dir,
                 "converted-data",
                 case_name)
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")

    if (json_net_data_by_components_file == "") || (
        json_net_data_by_components_file == nothing )

        net_data_by_components_file =
            joinpath(
                json_case_dir,
                "net_data_by_components_file.json")
    else

        net_data_by_components_file =
            joinpath(
                json_case_dir,
                json_net_data_by_components_file)

    end

    #----------------------------------------    
    
    fault_status =
        (no_fault = 0,
         on_fault = 1,
         clear_fault = 2,
         partial_clear_fault = 3)

    system_fault_status =
        [ fault_status.no_fault]
    
    #----------------------------------------
    
    (;u0_model_states_init,
     model_syms,
     model_mass_matrix,
     
     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     generic_system_dynamics_wt_fault_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states,
     cb_faults,
     cb_faults_no_resize,
     
     fault_nodes_idx,

     state_labels,
     algebraic_vars_labels,

     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,
     post_sta_PQ,

     vh,
     θh,
     gens_i_d,
     gens_i_q,

     algebraic_generic_model_wt_fault_sol_kwd_para,
     
     loc_load_exist,     

     Ybr_cal_and_edges_orientation) =
         NamedTupleTools.select(
             get_system_simulation_parameters_wt_faults(
                 net_data_by_components_file;
                 components_libs_dir =
                     components_libs_dir,
                basekV = basekV,    
                use_pu_in_PQ = use_pu_in_PQ,
                line_data_in_pu = line_data_in_pu,
                use_init_u0 = use_init_u0,
                use_nlsolve = use_nlsolve,

                pf_alg = pf_alg,

                abstol = abstol,
                reltol = reltol,

                on_fault_time = on_fault_time,
                clear_fault_time = clear_fault_time,

                list_fault_point_from_node_a =
                    list_fault_point_from_node_a,

                list_fault_resistance =
                    list_fault_resistance,

                list_no_line_circuit =
                    list_no_line_circuit,

                list_edges_to_have_fault =
                    list_edges_to_have_fault,

                clear_fault_selection_list =
                    clear_fault_selection_list,
             with_faults =
                     with_faults),
             (:u0_model_states_init,
              :model_syms,
              :model_mass_matrix,
              
              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :generic_system_dynamics_wt_fault_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states,
              :cb_faults,
              :cb_faults_no_resize,
              
              :fault_nodes_idx,

              :state_labels,
              :algebraic_vars_labels,

              :δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,
              :post_sta_PQ,

              :vh,
              :θh,
              :gens_i_d,
              :gens_i_q,

              :algebraic_generic_model_wt_fault_sol_kwd_para,
              :loc_load_exist,
              
              :Ybr_cal_and_edges_orientation))

    (;P_non_gens,
     Q_non_gens, 
     P_g_loc_load,
     Q_g_loc_load,) =
         NamedTupleTools.select(
             post_sta_PQ ,
             (:P_non_gens,
              :Q_non_gens, 
              :P_g_loc_load,
              :Q_g_loc_load) )

    (;loc_load_exist,
     dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
     dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     cleared_selected_lines_faults_net_para,

     pf_generic_gens_para) =
        NamedTupleTools.select(
        generic_system_dynamics_wt_fault_kwd_para,
        (:loc_load_exist,
         :dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
         :dyn_pf_vh_vhf_θh_θhf_id_iq_Idx,

         :dyn_pf_fun_kwd_n2s_idxs,
         :dyn_pf_fun_kwd_net_idxs,

         :Ynet_wt_nodes_idx_wt_adjacent_nodes,
         :cleared_selected_lines_faults_net_para,

         :pf_generic_gens_para))
    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    
    (flat_vh_flat_θh_id_iq_u0,
     dyn_pf_solver,
     algebraic_generic_model_wt_fault_kwd_para) =
         NamedTupleTools.select(
             algebraic_generic_model_wt_fault_sol_kwd_para,
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
    
    if system_status == :pre_fault_state

        pf_fun_mismatch =
           pre_fault_state_algebraic_generic_pf_ΔPQ_mismatch!

    elseif system_status == :fault_state 

        pf_fun_mismatch =
            fault_state_algebraic_generic_pf_ΔPQ_mismatch!

    elseif system_status == :post_fault_state

        pf_fun_mismatch =
         post_fault_state_algebraic_generic_pf_ΔPQ_mismatch!

    else
        
        throw(" system status not known ")
        
    end

    #----------------------------------------

    number_of_faults =
        length(list_edges_to_have_fault)

    twice_number_of_faults =
        2 * number_of_faults

    fault_vh_sym =
        generate_labels_by_nodes_idxs_and_vars(
            fault_nodes_idx,
            [:vh];
            label_prefix = "bus")

    fault_θh_sym =
        generate_labels_by_nodes_idxs_and_vars(
            fault_nodes_idx,
            [:θh];
            label_prefix = "bus")

    fault_nodes_sym =
        [fault_vh_sym;
         fault_θh_sym ]  

    fault_algeb_init =[
        ones(number_of_faults);
        zeros(number_of_faults)
        ]

    #----------------------------------------
    
    flat_vh_flat_θh_id_iq_vfh_θfh =
        [vh;
         θh;
         gens_i_d;
         gens_i_q;
         fault_algeb_init ]

    #----------------------------------------    
        
    if use_nlsolve == true

        sol = nlsolve((g, x) -> pf_fun_mismatch(
            g, x,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll;
                kwd_para =
                 algebraic_generic_model_wt_fault_kwd_para),
                      flat_vh_flat_θh_id_iq_vfh_θfh,
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
                algebraic_generic_model_wt_fault_kwd_para)),
            flat_vh_flat_θh_id_iq_vfh_θfh,
            δ_ed_dash_eq_dash_Png_Qng_Pll_Qll),
                                    pf_alg )
        
    end

    (dyn_pf_vh_Idxs, ) =
         NamedTupleTools.select(
             dyn_pf_vh_θh_id_iq_vhf_θhf_Idx,
             (:dyn_pf_vh_Idxs, ))

    pf_generic_gens_para =
        getproperty(
            generic_system_dynamics_wt_fault_kwd_para,
        :pf_generic_gens_para )
    
    X_d_dash = getproperty(
        getproperty(
            generic_system_dynamics_wt_fault_kwd_para,
            :pf_generic_gens_para), :X_d_dash)
    
    # (;X_d_dash,) =
    #     NamedTupleTools.select(
    #         pf_generic_gens_para,
    #         (:X_d_dash))
    
    y_aug_kw_para =
        (;loc_load_exist,
         X_d_dash,    
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs)
    
    return get_Yint_and_Yred_matrices(
        # vh,
        sol[dyn_pf_vh_Idxs],
        Ynet_wt_nodes_idx_wt_adjacent_nodes,
        P_non_gens,
        Q_non_gens,
        P_g_loc_load,
        Q_g_loc_load;
        y_aug_kw_para =
            y_aug_kw_para,
        fault_status =
            system_status ,
        cleared_selected_lines_faults_net_para =
            cleared_selected_lines_faults_net_para)
    
    
end

#-------------------------------------------------------
#-------------------------------------------------------
