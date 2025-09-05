# ####################################################
# using Pkg
# PowerSystemsCoSim_folder = joinpath(@__DIR__,"../..")
# cd(PowerSystemsCoSim_folder)
# Pkg.activate(PowerSystemsCoSim_folder)

# get_system_simulation_parameters

"""

TexTables.jl

LaTeXTabulars.jl

"""

#---------------------------------------------------
# benchmark powerflow functions
#---------------------------------------------------

function plot_benchmark_results(
    results;
    bm_group_sym =
        :ΔI_mismatch,
    
    fun_string_Ybus =
        "get_ΔI_mismatch_by_Ybus",
    fun_string_Ynet =
        "get_ΔI_mismatch_by_Ynet",
    
    network_cases =
        ["case4gs", "case5", "case9",
         "case14",  "case30", "case39",
         "case57", "case118", "case300",
         "case1354pegase", "case2869pegase"],
    
    net_cases_abr =
        ["c4", "c5", "c9", "c14", "c30",
         "c39", "c57", "c118", "c300", "c1354", "c2869"],
    x_coords =
        [4, 5, 9, 14, 30, 39, 57, 118, 300, 1354, 2869],
    
    ctg_string_by_Ybus =
        "ΔI mismatch by Ybus",
    ctg_string_by_Ynet =
        "ΔI mismatch by Ynet",

    benchmark_Δ_mismatch_csv_filename =
        benchmark_ΔI_mismatch_csv_filename,
    # benchmark_ΔPQ_mismatch_csv_filename
    
    benchmark_time_Δ_mismatch_csv_filename =
        benchmark_time_ΔI_mismatch_csv_filename,
        # benchmark_time_ΔPQ_mismatch_csv_filename,

    groupedbar_filename_string =
        "ΔI-mismatch-Ynet-Ybus-groupedbar.pdf",
    scatter_filename_string =
        "ΔI-mismatch-Ynet-Ybus-scatter.pdf",
    
    lable_scatter_mismatch_by_Ybus =
        "ΔI mismatch by Ybus",
    lable_scatter_mismatch_by_Ynet =
        "ΔI mismatch by Ynet" )
    
    df_mismatch = DataFrame( results[bm_group_sym] )


    tran_df_mismatch_time_in_ns =
        transform(
            df_mismatch,
            :first => ByRow(identity) =>
                [bm_group_sym, :case],
            :second => (t -> time.(t) ) => :Time )

    sel_tran_df_mismatch_time_in_ns =
        DataFrames.select(
            tran_df_mismatch_time_in_ns,
            bm_group_sym,
            :case,
            :Time )
        
    CSV.write(benchmark_time_Δ_mismatch_csv_filename,
              sel_tran_df_mismatch_time_in_ns)

    tran_df_mismatch =
        transform(
            df_mismatch,
            :first => ByRow(identity) =>
                [bm_group_sym, :case],
            :second => (t -> log10.(time.(t))) => :Time )

    sel_tran_df_mismatch =
        DataFrames.select(
            tran_df_mismatch,
            bm_group_sym,
            :case,
            :Time)
    
    CSV.write(benchmark_Δ_mismatch_csv_filename,
              sel_tran_df_mismatch)

    # Select rows where column
    # :ΔI_mismatch == "get_ΔI_mismatch_by_Ybus"

    Ybus_sub_sel_tran_df_mismatch =
        subset(sel_tran_df_mismatch,
               bm_group_sym => ByRow(x -> x ==
                   fun_string_Ybus ))


    Ynet_sub_sel_tran_df_mismatch =
        subset(sel_tran_df_mismatch,
               bm_group_sym => ByRow(x -> x ==
                   fun_string_Ynet ))

    #--------------------

    cases_sort_order =
        network_cases

    Ybus_sort_df =
        Ybus_sub_sel_tran_df_mismatch[
            indexin(cases_sort_order,
                    Ybus_sub_sel_tran_df_mismatch.case),:]

    Ynet_sort_df =
        Ynet_sub_sel_tran_df_mismatch[
            indexin(cases_sort_order,
                    Ynet_sub_sel_tran_df_mismatch.case),:]
    
    # net_cases =
    #     ["c4", "c5", "c9", "c14", "c30",
    #      "c39", "c57", "c118", "c300", "c1354", "c2869"]

    # x_coords =
    #     [4, 5, 9, 14, 30, 39, 57, 118, 300, 1354, 2869]
    
    net_cases =
        net_cases_abr

    Ybus_t_values = Ybus_sort_df.Time

    Ynet_t_values = Ynet_sort_df.Time

    labels = net_cases

    # https://github.com/JuliaPlots/StatsPlots.jl#grouped-bar-plots

    # ctg_string_by_Ybus = "ΔI mismatch by Ybus",
    # ctg_string_by_Ynet = "ΔI mismatch by Ynet",

    ctg = repeat([ctg_string_by_Ybus,
                  ctg_string_by_Ynet],
                 inner = 11)

    nam = repeat(net_cases,
                 outer = 2)


    # groupedbar_filename_string =
    #     "ΔI-mismatch-Ynet-Ybus-groupedbar.pdf"
    
    # scatter_filename_string =
    #     "ΔI-mismatch-Ynet-Ybus-scatter.pdf"

    plt_mismatch_groupedbar =
        StatsPlots.groupedbar(
            nam, [Ybus_t_values Ynet_t_values],
            group = ctg,
            xlabel = "Cases",
            ylabel = "Time log10([ns])",
            # title = "ΔI mismatch by admittances",
            bar_width = 0.67,
            lw = 0,
            framestyle = :box)

    plt_scatter_ΔI_mismatch =
        scatter(x_coords, Ybus_t_values, 
            xlabel="cases [number of buses]", 
            ylabel="Time log10([ns])", 
            label = lable_scatter_mismatch_by_Ybus)

    scatter!(plt_scatter_ΔI_mismatch,
             x_coords, Ynet_t_values,
             label= lable_scatter_mismatch_by_Ynet )

        
    figure_groupedbar_filename =
        joinpath(figure_dir,groupedbar_filename_string)

    figure_scatter_filename =
        joinpath(figure_dir,scatter_filename_string)

    savefig(plt_mismatch_groupedbar,
            figure_groupedbar_filename)

    
    savefig(plt_scatter_ΔI_mismatch,
            figure_scatter_filename)

    
end


function create_benchmarkgroups_admittance_vector_funcs(
    SUITE,
    network_cases,
    nt_admittance_vector_funcs;
    data_dir = "",    
    static_net_data_by_components_file =  "",
    components_libs_dir = "" )

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

    #--------------------------------------    
    
    if ( data_dir == "" ) || ( data_dir == nothing )

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end
    
    #--------------------------------------
    
    for case_name in network_cases

        json_case_dir =
            joinpath(
                data_dir,
                "converted-data",
                case_name,
                "json")

        @show case_name


        if (static_net_data_by_components_file == "" ||
            static_net_data_by_components_file == nothing)

            json_net_data_by_components_file =
                joinpath(json_case_dir,
                      "opf-pf-net-default-static-data.json")
        else

            json_net_data_by_components_file =
                joinpath(json_case_dir,
                         static_net_data_by_components_file )

        end

        #--------------------------------------
        
        # static_net_data_by_components_file =  ""
        #-------------------------------

        system_net_static_data =
            get_system_net_static_data(
                case_name;
                script_dir=
                    script_dir,
                data_dir = data_dir,
                json_net_data_by_components_file =
                    json_net_data_by_components_file,
                components_libs_dir = components_libs_dir,
                basekV              = 1.0,    
                use_pu_in_PQ        = true,
                line_data_in_pu     = true,
                pf_alg              =
                    NewtonRaphson(),
                no_lines_fault = 1)
        
        #-------------------------------

        (;edge_data_from_json,
         shunt_data_from_json,
         baseMVA,
         all_nodes_idx,
         n2s_all_nodes_idx) =
             NamedTupleTools.select(
                 system_net_static_data,
                 (:edge_data_from_json,
                  :shunt_data_from_json,
                  :baseMVA,
                  :all_nodes_idx,
                  :n2s_all_nodes_idx ) )
        
        #-------------------------------
        
        for a_admitt_Y_func in values(
            nt_admittance_vector_funcs)

            if nameof(a_admitt_Y_func) == :get_Ynet

                SUITE[:network_admittance][
                    string(a_admitt_Y_func),case_name ] =
                        @benchmarkable(
                            $(a_admitt_Y_func)(
                                $edge_data_from_json,
                                $shunt_data_from_json;
                                baseMVA = $baseMVA,
                                baseShunt = $baseMVA,
                                basekV = 1.0,
                               line_data_in_pu = true ),
                            evals = 1,
                            samples = 1000)

            elseif nameof(a_admitt_Y_func) == :get_Ynet_sp_sh

                SUITE[:network_admittance][
                    string(a_admitt_Y_func),case_name ] =
                        @benchmarkable(
                            $(a_admitt_Y_func)(
                                $edge_data_from_json,
                                $shunt_data_from_json;
                                baseMVA = $baseMVA,
                                baseShunt = $baseMVA,
                                basekV = 1.0,
                               line_data_in_pu = true ),
                            evals = 1,
                            samples = 1000)
                
            elseif nameof(a_admitt_Y_func) == :get_Yπ_net

                SUITE[:network_admittance][
                    string(a_admitt_Y_func),case_name ] =
                        @benchmarkable(
                            $(a_admitt_Y_func)(
                                $edge_data_from_json,
                                $shunt_data_from_json;
                                baseMVA = $baseMVA,
                                baseShunt = $baseMVA,
                                basekV = 1.0,
                                line_data_in_pu = true,
                                orientated_bool = false),
                            evals = 1,
                            samples = 1000)
                
            elseif nameof(a_admitt_Y_func) == :get_Ybus

                SUITE[:network_admittance][
                    string(a_admitt_Y_func),case_name ] =
                        @benchmarkable(
                            $(a_admitt_Y_func)(
                                $edge_data_from_json,
                                $shunt_data_from_json;
                                baseMVA = $baseMVA,
                                basekV = 1.0,
                                line_data_in_pu = true ),
                            evals = 1,
                            samples = 1000)
                
            else
                nothing
            end
        end
    end

end



function create_benchmarkgroups_pf_mismatch(
    SUITE,
    network_cases,
    mismatch_funcs_types;
    
    data_dir = "",    
    static_net_data_by_components_file = "",
    components_libs_dir = "" )

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

    #--------------------------------------    
    
    if ( data_dir == "" ) || ( data_dir == nothing )

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end
    
    #--------------------------------------

    for case_name in network_cases

        json_case_dir =
            joinpath(
                data_dir,
                "converted-data",
                case_name,
                "json")

        @show case_name


        if (static_net_data_by_components_file == "" ||
            static_net_data_by_components_file == nothing)

            json_net_data_by_components_file =
                joinpath(json_case_dir,
                      "opf-pf-net-default-static-data.json")
        else

            json_net_data_by_components_file =
                joinpath(json_case_dir,
                         static_net_data_by_components_file )

        end
        
        #-------------------------------

        system_net_static_data =
            get_system_net_static_data(
                case_name;
                script_dir=
                    script_dir,
                data_dir = data_dir,
                json_net_data_by_components_file =
                   json_net_data_by_components_file,
                components_libs_dir = components_libs_dir,
                basekV              = 1.0,    
                use_pu_in_PQ        = true,
                line_data_in_pu     = true,
                pf_alg              =
                    NewtonRaphson(),
                no_lines_fault = 1)
        
        #-------------------------------

        (sta_red_vh_θh_0,
         pf_PQ_param,
         pf_kw_para) =
             NamedTupleTools.select(
                 system_net_static_data,
                 (:sta_red_vh_θh_0,
                  :pf_PQ_param,
                  :pf_kw_para))

        #-------------------------------

        red_ΔPQ_x = similar(sta_red_vh_θh_0)

        Jac_row_size =
            Jac_col_size = length( sta_red_vh_θh_0 )

        Jac_vh_θh =
            spzeros( Jac_row_size, Jac_col_size )
        
        #-------------------------------

        (;
         edge_data_from_json,
         shunt_data_from_json,
         baseMVA,
         all_nodes_idx,
         n2s_all_nodes_idx) =
             NamedTupleTools.select(
                 system_net_static_data,
                 (
                  :edge_data_from_json,
                  :shunt_data_from_json,
                  :baseMVA,
                  :all_nodes_idx,
                  :n2s_all_nodes_idx ) )

        #-------------------------------
        
        for a_mismatch_func in values(mismatch_funcs_types)

            if nameof(a_mismatch_func) ==
                :get_ΔI_mismatch_by_Ynet 

                Ynet_wt_nodes_idx_wt_adjacent_nodes =
                    get_Ynet( edge_data_from_json,
                          shunt_data_from_json;
                          baseMVA = baseMVA,
                          baseShunt = baseMVA,
                          basekV = 1.0,
                          line_data_in_pu = true )
                
                pf_prob = NonlinearProblem(
                    NonlinearFunction( ( g, x, p ) ->
                                          a_mismatch_func(
                                              g, x, p;
                            Ynet_wt_nodes_idx_wt_adjacent_nodes =
                                Ynet_wt_nodes_idx_wt_adjacent_nodes,
                                              pf_kw_para =
                                     disaggregate_sta_pf_keywords_parameter(
                                       pf_kw_para)),
                                                  jac = (Jac_vh_θh, x, p) ->
                                      sta_pf_Jac!_df_dx_by_Ynet_or_Yπ_net!(
                                           Jac_vh_θh, x, p;
                                          pf_kw_para =
                                              pf_kw_para,
                                               func =
                                                  a_mismatch_func,
                                               net_addmitance_tuple =
                             Ynet_wt_nodes_idx_wt_adjacent_nodes,
                                           by_Ynet_or_Yπ_net =
                                               :Ynet) ),
                                      sta_red_vh_θh_0,
                    pf_PQ_param)

                                
                SUITE[:ΔI_mismatch][
                    string(a_mismatch_func), case_name]  =
                        @benchmarkable(
                    NonlinearSolve.solve(
                        $pf_prob,
                        NewtonRaphson();
                        abstol = abstol,
                        reltol = reltol ),
                    evals = 1,
                    samples = 1000)
                
                # @benchmarkable NonlinearSolve.solve(
                #   $pf_prob,
                #    NewtonRaphson();
                #    abstol = abstol,
                #    reltol = reltol )


            elseif  nameof(a_mismatch_func) ==
                :get_ΔPQ_mismatch_by_Ynet 

                Ynet_wt_nodes_idx_wt_adjacent_nodes =
                    get_Ynet( edge_data_from_json,
                          shunt_data_from_json;
                          baseMVA = baseMVA,
                          baseShunt = baseMVA,
                          basekV = 1.0,
                          line_data_in_pu = true )
                
                pf_prob = NonlinearProblem(
                    NonlinearFunction( ( g, x, p ) ->
                                          a_mismatch_func(
                                              g, x, p;
                            Ynet_wt_nodes_idx_wt_adjacent_nodes =
                                Ynet_wt_nodes_idx_wt_adjacent_nodes,
                                              pf_kw_para =
                                     disaggregate_sta_pf_keywords_parameter(
                                       pf_kw_para)),
                                                  jac = (Jac_vh_θh, x, p) ->
                                      sta_pf_Jac!_df_dx_by_Ynet_or_Yπ_net!(
                                           Jac_vh_θh, x, p;
                                          pf_kw_para =
                                              pf_kw_para,
                                               func =
                                                  a_mismatch_func,
                                               net_addmitance_tuple =
                             Ynet_wt_nodes_idx_wt_adjacent_nodes,
                                           by_Ynet_or_Yπ_net =
                                               :Ynet) ),
                                      sta_red_vh_θh_0,
                    pf_PQ_param)

                                
                SUITE[:ΔPQ_mismatch][
                    string(a_mismatch_func), case_name]  =
                        @benchmarkable(
                    NonlinearSolve.solve(
                        $pf_prob,
                        NewtonRaphson();
                        abstol = abstol,
                        reltol = reltol ),
                    evals = 1,
                    samples = 1000)
                
                # @benchmarkable NonlinearSolve.solve(
                #   $pf_prob,
                #    NewtonRaphson();
                #    abstol = abstol,
                #    reltol = reltol )

                
            elseif nameof(a_mismatch_func) ==
                :get_ΔI_mismatch_by_Yπ_net 

                Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
                    get_Yπ_net(edge_data_from_json,
                             shunt_data_from_json;
                             baseMVA = baseMVA,
                             baseShunt = baseMVA,
                             basekV = 1.0,
                             line_data_in_pu = true,
                             orientated_bool = false)
                
                   pf_prob = NonlinearProblem(
                       NonlinearFunction(( g, x, p ) ->
                           a_mismatch_func(
                               g, x, p;
             Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
                  Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
                               pf_kw_para =
                      disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para)),
                                   jac = (Jac_vh_θh, x, p) ->
                       sta_pf_Jac!_df_dx_by_Ynet_or_Yπ_net!(
                            Jac_vh_θh, x, p;
                                pf_kw_para = pf_kw_para,
                                func =
                                   a_mismatch_func,
                                net_addmitance_tuple =
              Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
                            by_Ynet_or_Yπ_net =
                                :Yπ_net)) ,
                       sta_red_vh_θh_0,
                       pf_PQ_param)

                                
                SUITE[:ΔI_mismatch][
                    string(a_mismatch_func), case_name]  = 
                @benchmarkable(
                    NonlinearSolve.solve(
                        $pf_prob,
                        NewtonRaphson();
                        abstol = abstol,
                        reltol = reltol ),
                    evals = 1,
                    samples = 1000)

                
            elseif  nameof(a_mismatch_func) ==
                :get_ΔPQ_mismatch_by_Yπ_net 

                Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
                    get_Yπ_net(edge_data_from_json,
                             shunt_data_from_json;
                             baseMVA = baseMVA,
                             baseShunt = baseMVA,
                             basekV = 1.0,
                             line_data_in_pu = true,
                             orientated_bool = false)
                
                   pf_prob = NonlinearProblem(
                       NonlinearFunction(( g, x, p ) ->
                           a_mismatch_func(
                               g, x, p;
             Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
                  Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
                               pf_kw_para =
                      disaggregate_sta_pf_keywords_parameter(
                        pf_kw_para)),
                                   jac = (Jac_vh_θh, x, p) ->
                       sta_pf_Jac!_df_dx_by_Ynet_or_Yπ_net!(
                            Jac_vh_θh, x, p;
                                pf_kw_para = pf_kw_para,
                                func =
                                   a_mismatch_func,
                                net_addmitance_tuple =
              Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes,
                            by_Ynet_or_Yπ_net =
                                :Yπ_net)) ,
                       sta_red_vh_θh_0,
                       pf_PQ_param)

                                
                SUITE[:ΔPQ_mismatch][
                    string(a_mismatch_func), case_name]  = 
                @benchmarkable(
                    NonlinearSolve.solve(
                        $pf_prob,
                        NewtonRaphson();
                        abstol = abstol,
                        reltol = reltol ),
                    evals = 1,
                    samples = 1000)
                
            elseif nameof(a_mismatch_func) ==
                :get_ΔI_mismatch_by_Ybus
                
                Ybus = getproperty(get_Ybus(
                    edge_data_from_json,
                         shunt_data_from_json;
                         baseMVA = baseMVA,
                         basekV = 1.0,
                         line_data_in_pu = true ), :Ybus)

               pf_prob = NonlinearProblem(
                   NonlinearFunction( ( g, x, p ) ->
                       a_mismatch_func(
                           g, x, p;
                           Ybus =
                               Ybus,
                           use_autodiff =
                               false,
                           pf_kw_para =
                  disaggregate_sta_pf_keywords_parameter(
                    pf_kw_para)),
                               jac = (Jac_vh_θh, x, p) ->
                   sta_pf_Jac!(
                        Jac_vh_θh, x, p;
                            pf_kw_para =
                  disaggregate_sta_pf_keywords_parameter(
                      pf_kw_para),
                       Ybus =
                           Ybus ) ),
                   sta_red_vh_θh_0,
                   pf_PQ_param)                
                                
                SUITE[:ΔI_mismatch][
                    string(a_mismatch_func), case_name]  =
                @benchmarkable(
                    NonlinearSolve.solve(
                        $pf_prob,
                        NewtonRaphson();
                        abstol = abstol,
                        reltol = reltol ),
                    evals = 1,
                    samples = 1000)
                
            elseif nameof(a_mismatch_func) ==
                :get_ΔPQ_mismatch_by_Ybus 

                
                Ybus = getproperty(get_Ybus(
                    edge_data_from_json,
                         shunt_data_from_json;
                         baseMVA = baseMVA,
                         basekV = 1.0,
                         line_data_in_pu = true ), :Ybus)

               pf_prob = NonlinearProblem(
                   NonlinearFunction( ( g, x, p ) ->
                       a_mismatch_func(
                           g, x, p;
                           Ybus =
                               Ybus,
                           use_autodiff =
                               false,
                           pf_kw_para =
                  disaggregate_sta_pf_keywords_parameter(
                    pf_kw_para)),
                               jac = (Jac_vh_θh, x, p) ->
                   sta_pf_Jac!(
                        Jac_vh_θh, x, p;
                            pf_kw_para =
                  disaggregate_sta_pf_keywords_parameter(
                      pf_kw_para),
                       Ybus =
                           Ybus ) ),
                   sta_red_vh_θh_0,
                   pf_PQ_param)                
                                
                SUITE[:ΔPQ_mismatch][
                    string(a_mismatch_func), case_name]  =
                @benchmarkable(
                    NonlinearSolve.solve(
                        $pf_prob,
                        NewtonRaphson();
                        abstol = abstol,
                        reltol = reltol ),
                    evals = 1,
                    samples = 1000)
                
            else
                nothing
            end
                                        
        end
    end
            
end


#---------------------------------------------------
# sim Functions
#---------------------------------------------------

function sim_full_model_distributed_slack_pf(
    net_data_by_components_file;
    components_libs_dir = "",
    basekV            = 1.0,    
    use_pu_in_PQ      = true,
    line_data_in_pu   = true,
    with_faults       = false,
    wt_network_loss_by_sta_pf_PQ_para_bool = true,
    fractional_digits = 6,
    pf_alg            = NewtonRaphson(),
    abstol            = 1e-12,

    reltol            = 1e-12 )

    
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

    #--------------------------------------    
    
    (;dynamic_parameters,
     pf_parameters,
     opf_parameters,
     states_Idx_syms_wt_functions) =
         NamedTupleTools.select(
             get_generic_system_simulation_parameters(
                 net_data_by_components_file;
                 components_libs_dir,
                 basekV          = basekV,    
                 use_pu_in_PQ    = use_pu_in_PQ,
                 line_data_in_pu = line_data_in_pu,
                 with_faults     = with_faults,
                 pf_alg          =  pf_alg ),
             (:dynamic_parameters,
              :pf_parameters,
              :opf_parameters,
              :states_Idx_syms_wt_functions ))

    
        (;baseMVA,
         basekV,
         loc_load_exist,
         
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx,
            
         nodes_with_demands_idx,
            
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_n2s_idxs,

         
         generic_gens_para,
         ode_gens_para,
         gens_vh_slack_θh_para,

         pf_sta_ΔPQ_mismatch_parameters,
         sta_pf_PQ_para,
         pf_PQ_param,
         pf_generic_gens_para,
         opf_generic_gens_para) =
             NamedTupleTools.select(
                 pf_parameters,
                 (:baseMVA,
                  :basekV,
                  :loc_load_exist,
                  
                  :slack_gens_nodes_idx,
                  :non_slack_gens_nodes_idx,
                  :gens_nodes_idx,
                  :non_gens_nodes_idx,
                  :gens_nodes_with_loc_loads_idx,
                  :all_nodes_idx,
                  
                  :nodes_with_demands_idx,
                  
                  :dyn_pf_fun_kwd_net_idxs,
                  :dyn_pf_fun_kwd_n2s_idxs,

                                    
                  :generic_gens_para,
                  :ode_gens_para,
                  :gens_vh_slack_θh_para,

                  :pf_sta_ΔPQ_mismatch_parameters,
                  :sta_pf_PQ_para,
                  :pf_PQ_param,
                  :pf_generic_gens_para,
                  :opf_generic_gens_para) ) 

        
        (;states_Idx_syms_wt_functions,
         
         plants_cb_paras_switches,
         cb_states,
         
         baseMVA,
         basekV,

         loc_load_exist,

         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx,            
         nodes_with_demands_idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         ode_gens_para,

         sta_pf_PQ_para,
         pf_PQ_param,
         pf_generic_gens_para,

         edges_fbus,
         edges_tbus,
         edges_r,
         edges_x,
         edges_b,
         edges_ratio,
         edges_angle,
         edges_type,
         shunt_idx,
         Gs,
         Bs,

         Ynet_wt_nodes_idx_wt_adjacent_nodes,
         Ybr_cal_and_edges_orientation,

         kwd_sta_sta_ΔPQ_sol_by_json,
         generic_red_sol_kwd_para,
         generic_dyn_sol_kwd_para,

         plants_init_kwd_para,
         plants_cb_paras_kwd_para,

         y_aug_kw_para,

         with_faults,

         edges_r_x_b_ratio_angle_idx,

         Ynet_rows_Idxs_in_flattend,
         Ynet_real_imag_Idxs_in_flattend ) =
             NamedTupleTools.select(
                 dynamic_parameters,
                 (:states_Idx_syms_wt_functions,
         
                  :plants_cb_paras_switches,
                  :cb_states,

                  :baseMVA,
                  :basekV,

                  :loc_load_exist,

                  :slack_gens_nodes_idx,
                  :non_slack_gens_nodes_idx,
                  :gens_nodes_idx,
                  :non_gens_nodes_idx,
                  :gens_nodes_with_loc_loads_idx,
                  :all_nodes_idx,            
                  :nodes_with_demands_idx,

                  :dyn_pf_fun_kwd_n2s_idxs,
                  :dyn_pf_fun_kwd_net_idxs,

                  :generic_gens_para,
                  :generic_avrs_para,
                  :generic_govs_para,

                  :ode_gens_para,

                  :sta_pf_PQ_para,
                  :pf_PQ_param,
                  :pf_generic_gens_para,

                  :edges_fbus,
                  :edges_tbus,
                  :edges_r,
                  :edges_x,
                  :edges_b,
                  :edges_ratio,
                  :edges_angle,
                  :edges_type,
                  :shunt_idx,
                  :Gs,
                  :Bs,

                  :Ynet_wt_nodes_idx_wt_adjacent_nodes,
                  :Ybr_cal_and_edges_orientation,

                  :kwd_sta_sta_ΔPQ_sol_by_json,
                  :generic_red_sol_kwd_para,
                  :generic_dyn_sol_kwd_para,

                  :plants_init_kwd_para,
                  :plants_cb_paras_kwd_para,

                  :y_aug_kw_para,

                  :with_faults,

                  :edges_r_x_b_ratio_angle_idx,

                  :Ynet_rows_Idxs_in_flattend,
                  :Ynet_real_imag_Idxs_in_flattend ))

    
    (;

     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     pf_vh_θh_idx_and_idx2Idx,

     dyn_pf_flat_vh_flat_θh_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     scale_Pg_Png_Qng_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx) = 
         NamedTupleTools.select(
             states_Idx_syms_wt_functions,

             (:Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,

              :pf_vh_θh_idx_and_idx2Idx,
              
              :dyn_pf_flat_vh_flat_θh_Idx,
              
              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :scale_Pg_Png_Qng_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx))
    
    #############################

    #############################

    (non_gens_vh_idx,
     non_slack_gens_θh_idx,
     non_gens_θh_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     red_non_gens_vh_Idxs,
     red_non_slack_gens_θh_Idxs,
     red_non_gens_θh_Idxs,

     red_slack_value_Idxs,

     gens_nodes_idx,

     non_gens_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,

     n2s_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_non_gens_idx,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             pf_vh_θh_idx_and_idx2Idx,
             (:non_gens_vh_idx,
              :non_slack_gens_θh_idx,
              :non_gens_θh_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :red_non_gens_vh_Idxs,
              :red_non_slack_gens_θh_Idxs,
              :red_non_gens_θh_Idxs,

              :red_slack_value_Idxs,

              :gens_nodes_idx,

              :non_gens_nodes_idx,
              :non_slack_gens_and_non_gens_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,

              :n2s_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_non_gens_idx,
              :n2s_all_nodes_idx))

    (;Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         NamedTupleTools.select(
             Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx))

    (;edges_Ybr_cal,
     edges_orientation) =
         NamedTupleTools.select(
             Ybr_cal_and_edges_orientation,
             (:edges_Ybr_cal,
              :edges_orientation ))

    #--------------------------------------------
    #--------------------------------------------

    (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             Pg_Qg_Png_Qng_Pll_Qll_Idx,
             (:dyn_P_gens_Idxs,
              :dyn_Q_gens_Idxs,
              :dyn_P_non_gens_Idxs,
              :dyn_Q_non_gens_Idxs,
              :dyn_P_gens_loc_load_Idxs,
              :dyn_Q_gens_loc_load_Idxs))

    #--------------------------------------------

    red_vh_θh_idx =
        getproperty(
            pf_vh_θh_idx_and_idx2Idx,
            :red_vh_θh_idx)

    #--------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend ) =
         NamedTupleTools.select(
             get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                 Ynet),
             (:Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend))

    (;Ynet_rows_Idxs_in_flattend,
    Ynet_real_imag_Idxs_in_flattend ) =
        NamedTupleTools.select(
            get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                getproperty(
                    Ynet_wt_nodes_idx_wt_adjacent_nodes,
                    :Ynet)),
            (:Ynet_rows_Idxs_in_flattend,
             :Ynet_real_imag_Idxs_in_flattend) )

    (;Ynet_real_Idxs,
     Ynet_imag_Idxs) =
        NamedTupleTools.select(
            Ynet_real_imag_Idxs_in_flattend,
            (:Ynet_real_Idxs,
             :Ynet_imag_Idxs))

    #--------------------------------------------


    (;r_Idxs, x_Idxs, b_Idxs,
     ratio_Idxs, angle_Idxs) =
         NamedTupleTools.select(
             edges_r_x_b_ratio_angle_idx,
             (:r_Idxs, :x_Idxs, :b_Idxs,
              :ratio_Idxs, :angle_Idxs))

    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))

    #--------------------------------------------
    #--------------------------------------------
    
    (Pg,
     Qg,
     Png,
     Qng, 
     Pll,
     Qll) =
        NamedTupleTools.select(
            sta_pf_PQ_para,
            (:P_gens,
             :Q_gens,
             :P_non_gens,
             :Q_non_gens,
             :P_g_loc_load,
             :Q_g_loc_load ) )

    
    Pg_Qg_Png_Qng_Pll_Qll =
        [Pg;
         Qg;
         Png;
         Qng; 
         Pll;
         Qll ]
    
    #--------------------------------------------
    #--------------------------------------------

    Pg_inj_Qg_inj_Png_Qng_Pll_Qll =
        get_Pg_inj_Qg_inj_Png_Qng_Pll_Qll(
            Pg_Qg_Png_Qng_Pll_Qll;
            loc_load_exist,
            Pg_Qg_Png_Qng_Pll_Qll_Idx,    
            gens_nodes_idx,
            n2s_gens_idx,
            gens_nodes_with_loc_loads_idx,
            n2s_gens_with_loc_load_idxs)

    Pg_inj_Qg_inj_Png_Qng =
        getproperty(
            Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
            :Pg_inj_Qg_inj_Png_Qng)

    (;gens_Pg_inj,
     gens_Qg_inj,
     P_non_gens,
     Q_non_gens,
     Pg_inj_Qg_inj_Png_Qng) =
         NamedTupleTools.select(
             Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
             (:gens_Pg_inj,
              :gens_Qg_inj,
              :P_non_gens,
              :Q_non_gens,
              :Pg_inj_Qg_inj_Png_Qng))

    #--------------------------------------------
    #--------------------------------------------

    #--------------------------------------------
    # Loss participation factors
    #--------------------------------------------

    gens_Sn =
        getproperty(ode_gens_para, :Sn)

    gens_loss_participation =
        get_loss_particpation_by_gens_rating(
            gens_Sn, Pg)

    gens_loss_participation_by_loading =
        get_loss_particpation_by_gens_loading( Pg)

    #--------------------------------------------
    # power disturbance resolution participation
    #--------------------------------------------

    gens_active_power_particpation_by_rating =
        get_gens_active_power_particpation_by_rating(
            gens_Sn, Pg)

    gens_active_power_particpation_by_loading =
        get_gens_active_power_particpation_by_loading(
            Pg )

    gens_reactive_power_particpation_by_rating =
        get_gens_reactive_power_particpation_by_rating(
            gens_Sn, Qg)

    gens_reactive_power_particpation_by_loading =
        get_gens_reactive_power_particpation_by_loading(
            Qg)


    active_power_disturbance_resolution_participation =
        gens_active_power_particpation_by_rating

    reactive_power_disturbance_resolution_participation =
        gens_reactive_power_particpation_by_loading

    #--------------------------------------------

    # slack_gens_vh_θh_gens_vh_non_slack_gens_vh =
    #     get_slack_gens_vh_θh_gens_vh_non_slack_gens_vh(
    #         vh,
    #         θh,
    #         dyn_pf_fun_kwd_net_idxs,
    #         dyn_pf_fun_kwd_n2s_idxs)

    slack_gens_vh_θh_gens_vh_non_slack_gens_vh  =
        gens_vh_slack_θh_para
    
    #--------------------------------------------
    # kwd parameters
    #--------------------------------------------

    pf_model_kwd_para =
        (;loc_load_exist,
         slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         gens_loss_participation,
         active_power_disturbance_resolution_participation,
         reactive_power_disturbance_resolution_participation,

         sta_pf_PQ_para,     

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         Pg_Qg_Png_Qng_Pll_Qll_Idx,
         Pg_Png_Qng_Idx,
         pf_vh_θh_idx_and_idx2Idx,

         scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
         scale_Pg_Png_Qng_Idx,
         dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

         shunt_idx,
         edges_fbus,
         edges_tbus,

         edges_type,
         edges_orientation,

         edges_r_x_b_ratio_angle_idx,
         Ynet_rows_Idxs_in_flattend,
         Ynet_real_imag_Idxs_in_flattend,

         baseMVA = baseMVA,
         basekV = basekV )

    generic_sol_kwd_para =
        (;
         Ybr_cal_and_edges_orientation,
         ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para )


    vh_θh =
        [ones(length( all_nodes_idx ));
         zeros(length(all_nodes_idx))]


    #--------------------------------------------

    red_vh_θh_slack_value =
        [vh_θh[non_gens_vh_idx];
         vh_θh[non_slack_gens_θh_idx];
         vh_θh[non_gens_θh_idx];
         [0.0]]

    non_gens_vh_all_θh_slack_value =
        [ones(length(all_nodes_idx ))[non_gens_nodes_idx];
         zeros(length(all_nodes_idx));
         [0.0]]    
    

    vh_θh_slack_value =
        [ones(length(all_nodes_idx ));
         zeros(length(all_nodes_idx));
         [0.0]]   
    #--------------------------------------------


    slack_gens_nodes_idx =
        getproperty(
            dyn_pf_fun_kwd_net_idxs,
            :slack_gens_nodes_idx)[1]

    # Subtract net loss from slack node Pg,
    # loss will subsequently be distributed

    if wt_network_loss_by_sta_pf_PQ_para_bool == true

        ds_gens_Pg_inj = @set gens_Pg_inj[slack_gens_nodes_idx]=
            gens_Pg_inj[slack_gens_nodes_idx] -
            get_total_P_network_loss_by_sta_pf_PQ_para(
                sta_pf_PQ_para,
                loc_load_exist )
        
    else
        
        ds_gens_Pg_inj = @set gens_Pg_inj[slack_gens_nodes_idx]=
            gens_Pg_inj[slack_gens_nodes_idx] -
            get_total_P_network_loss(
                vh,
                θh,
                Ynet;
                nodes_idx_with_adjacent_nodes_idx,
                n2s_all_nodes_idx,
                all_nodes_idx )
    end
    


    ds_Pg_inj_Qg_inj_Png_Qng =
        [ds_gens_Pg_inj;
         gens_Qg_inj;
         P_non_gens;
         Q_non_gens ]



    ds_Pg_Qg_Png_Qng =
        [Pg;
         Qg;
         P_non_gens;
         Q_non_gens ]


    ds_Pg_Qg_Png_Qng_Pll_Qll =
        [Pg;
         Qg;
         P_non_gens;
         Q_non_gens;
         Pll;
         Qll ]

    ds_Pg_inj_Png_Qng =
        [ds_gens_Pg_inj;
         P_non_gens;
         Q_non_gens]

    #--------------------------------------------

    # pf_fun_mismatch =
    #     get_red_model_distributed_slack_pf_ΔPQ_mismatch!


    # red_sol = NonlinearSolve.solve(
    #     NonlinearProblem(
    #         NonlinearFunction( (x, p) ->
    #             pf_fun_mismatch(
    #                 x, p;
    #                 pf_model_kwd_para =
    #                     pf_model_kwd_para )),
    #         red_vh_θh_slack_value,
    #         ds_Pg_inj_Png_Qng ), pf_alg,
    #     abstol = abstol,
    #     reltol = reltol)

    
    # results_pf_red_sol_u =
    #     get_generic_results_ds_pf_red_sol_u(
    #         red_sol;
    #         generic_sol_kwd_para =
    #             generic_sol_kwd_para )

    # results_conti_or_ds_pf_red_sol_u =
    #     get_generic_results_conti_or_ds_pf_red_sol_u(
    #         red_sol,
    #         ds_Pg_inj_Png_Qng;
    #         generic_sol_kwd_para =
    #             generic_sol_kwd_para)

    # -------------------------------------

    pf_fun_mismatch =
        get_model_distributed_slack_pf_ΔPQ_mismatch!


    pf_sol = NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( ( x, p ) ->
                pf_fun_mismatch(
                    x, p;
                    pf_model_kwd_para =
                        pf_model_kwd_para )),
            vh_θh_slack_value,
            ds_Pg_inj_Qg_inj_Png_Qng),
        pf_alg,
        abstol = abstol,
        reltol = reltol)

        
    # -------------------------------------
    
    retcode = pf_sol.retcode
    
    stats = pf_sol.stats


    results_ds_pf_sol_u =
        get_generic_results_ds_pf_sol_u(
            pf_sol;
            generic_sol_kwd_para =
                generic_sol_kwd_para)


    results_conti_or_ds_pf_sol_u =
        get_generic_results_conti_or_ds_pf_sol_u(
            pf_sol,
            ds_Pg_inj_Qg_inj_Png_Qng;
            generic_sol_kwd_para =
            generic_sol_kwd_para)


    ds_slack_value =
        getproperty(
            results_conti_or_ds_pf_sol_u,
            :slack_value)

    ds_vh =
        getproperty(
            results_pf_red_sol_u,
            :vh)

    ds_θh =
        getproperty(
            results_conti_or_ds_pf_sol_u,
            :θh)

    ds_θh_deg =
        getproperty(
            results_conti_or_ds_pf_sol_u,
            :θh_deg)

    ds_pf_P_gens =
        getproperty(
            results_conti_or_ds_pf_sol_u,
            :pf_P_gens)

    ds_pf_Q_gens =
        getproperty(
            results_conti_or_ds_pf_sol_u,
            :pf_Q_gens)

    ds_pf_P_g_gens =
        getproperty(
            results_conti_or_ds_pf_sol_u,
            :pf_P_g_gens)

    ds_pf_Q_g_gens =
        getproperty(
            results_conti_or_ds_pf_sol_u,
            :pf_Q_g_gens)

    transformed_gens_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in gens_nodes_idx ]


    transformed_gens_nodes_with_loc_loads_idx =
        [n2s_all_nodes_idx[idx]
         for idx in gens_nodes_with_loc_loads_idx ]
    

    transformed_non_gens_nodes_idx =
        [ n2s_all_nodes_idx[idx]
         for idx in non_gens_nodes_idx]
    
    df_gens  = DataFrames.DataFrame(
        ;Bus = gens_nodes_idx,
        vh = round.(ds_vh[transformed_gens_nodes_idx];
                    digits = fractional_digits),
        θh = round.(ds_θh[transformed_gens_nodes_idx];
                    digits = fractional_digits),
        θh_deg =
            round.(ds_θh_deg[transformed_gens_nodes_idx];
                   digits = fractional_digits),
        Pg = round.(ds_pf_P_gens;
                   digits = fractional_digits),
        Qg = round.(ds_pf_Q_gens;
                   digits = fractional_digits),)

    
    df_gens_with_loc_load  = DataFrames.DataFrame(
        ;Bus = gens_nodes_with_loc_loads_idx,
        vh = round.(ds_vh[transformed_gens_nodes_with_loc_loads_idx];
                    digits = fractional_digits),
        θh = round.(ds_θh[transformed_gens_nodes_with_loc_loads_idx];
                    digits = fractional_digits),
        θh_deg =
            round.(ds_θh_deg[transformed_gens_nodes_with_loc_loads_idx];
                   digits = fractional_digits),
        Pll = round.(Pll;
                   digits = fractional_digits),
        Qll = round.(Qll;
                   digits = fractional_digits),)

    
    df_non_gens = DataFrames.DataFrame(
        ;Bus = non_gens_nodes_idx,
        vh =round.(ds_vh[transformed_non_gens_nodes_idx];
                   digits = fractional_digits),
        θh = round.(ds_θh[transformed_non_gens_nodes_idx];
                    digits = fractional_digits),
        θh_deg = round.(ds_θh_deg[
            transformed_non_gens_nodes_idx];
                        digits = fractional_digits),
        Png = round.(P_non_gens;
                   digits = fractional_digits),
        Qng = round.(Q_non_gens;
                   digits = fractional_digits),)
    
    df_all_nodes_result = DataFrames.DataFrame(
        ;Bus = all_nodes_idx,
        vh = round.(ds_vh;digits=fractional_digits),
        θh = round.(ds_θh;digits=fractional_digits),
        θh_deg = round.(ds_θh_deg;digits=fractional_digits),
        P = round.( [idx ∈ gens_nodes_idx ? ds_pf_P_gens[ n2s_gens_idx[idx] ] : P_non_gens[ n2s_non_gens_idx[ idx ] ] for idx in all_nodes_idx ]; digits = fractional_digits),
        Q = round.( [idx ∈ gens_nodes_idx ? ds_pf_Q_gens[ n2s_gens_idx[idx] ] : Q_non_gens[ n2s_non_gens_idx[ idx ] ] for idx in all_nodes_idx ]; digits = fractional_digits),)


    return loc_load_exist == true ?
        (;ds_slack_value,
         ds_vh,
         ds_θh,
         ds_θh_deg,
         ds_pf_P_gens,
         ds_pf_Q_gens,
         ds_pf_P_g_gens,
         ds_pf_Q_g_gens,

         df_gens,
         df_gens_with_loc_load,
         df_non_gens,
         df_all_nodes_result,

         results_ds_pf_sol_u,
         results_conti_or_ds_pf_sol_u,
         retcode,

         stats,
         gens_loss_participation,
         gens_loss_participation_by_loading,

         gens_active_power_particpation_by_rating,
         gens_active_power_particpation_by_loading,
         gens_reactive_power_particpation_by_rating,
         gens_reactive_power_particpation_by_loading,

         active_power_disturbance_resolution_participation,
         reactive_power_disturbance_resolution_participation) :
             (;ds_slack_value,
              ds_vh,
              ds_θh,
              ds_θh_deg,
              ds_pf_P_gens,
              ds_pf_Q_gens,
              ds_pf_P_g_gens,
              ds_pf_Q_g_gens,

              df_gens,
              df_non_gens,
              df_all_nodes_result,

              results_ds_pf_sol_u,
              results_conti_or_ds_pf_sol_u,
              retcode,
              stats,
              gens_loss_participation,
              gens_loss_participation_by_loading,

              gens_active_power_particpation_by_rating,
              gens_active_power_particpation_by_loading,
              gens_reactive_power_particpation_by_rating,
              gens_reactive_power_particpation_by_loading,

              active_power_disturbance_resolution_participation,
              reactive_power_disturbance_resolution_participation )
    
end



function sim_red_model_distributed_slack_pf(
    net_data_by_components_file;
    components_libs_dir = "",
    basekV            = 1.0,    
    use_pu_in_PQ      = true,
    line_data_in_pu   = true,
    with_faults       = false,
    wt_network_loss_by_sta_pf_PQ_para_bool = true,
    fractional_digits = 6,
    pf_alg            = NewtonRaphson(),
    abstol      = 1e-12,

    reltol      = 1e-12)

    
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

    #--------------------------------------    
    
    (;dynamic_parameters,
     pf_parameters,
     opf_parameters,
     states_Idx_syms_wt_functions) =
         NamedTupleTools.select(
             get_generic_system_simulation_parameters(
                 net_data_by_components_file;
                 components_libs_dir,
                 basekV          = basekV ,    
                 use_pu_in_PQ    = use_pu_in_PQ,
                 line_data_in_pu = line_data_in_pu,
                 with_faults     = with_faults,
                 pf_alg          =  pf_alg ),
             (:dynamic_parameters,
              :pf_parameters,
              :opf_parameters,
              :states_Idx_syms_wt_functions ))

    
        (;baseMVA,
         basekV,
         loc_load_exist,
         
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx,
            
         nodes_with_demands_idx,
            
         dyn_pf_fun_kwd_net_idxs,
         dyn_pf_fun_kwd_n2s_idxs,

         
         generic_gens_para,
         ode_gens_para,
         gens_vh_slack_θh_para,

         pf_sta_ΔPQ_mismatch_parameters,
         sta_pf_PQ_para,
         pf_PQ_param,
         pf_generic_gens_para,
         opf_generic_gens_para) =
             NamedTupleTools.select(
                 pf_parameters,
                 (:baseMVA,
                  :basekV,
                  :loc_load_exist,
                  
                  :slack_gens_nodes_idx,
                  :non_slack_gens_nodes_idx,
                  :gens_nodes_idx,
                  :non_gens_nodes_idx,
                  :gens_nodes_with_loc_loads_idx,
                  :all_nodes_idx,
                  
                  :nodes_with_demands_idx,
                  
                  :dyn_pf_fun_kwd_net_idxs,
                  :dyn_pf_fun_kwd_n2s_idxs,

                                    
                  :generic_gens_para,
                  :ode_gens_para,
                  :gens_vh_slack_θh_para,

                  :pf_sta_ΔPQ_mismatch_parameters,
                  :sta_pf_PQ_para,
                  :pf_PQ_param,
                  :pf_generic_gens_para,
                  :opf_generic_gens_para) ) 

        
        (;states_Idx_syms_wt_functions,
         
         plants_cb_paras_switches,
         cb_states,
         
         baseMVA,
         basekV,

         loc_load_exist,

         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         all_nodes_idx,            
         nodes_with_demands_idx,

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         generic_gens_para,
         generic_avrs_para,
         generic_govs_para,

         ode_gens_para,

         sta_pf_PQ_para,
         pf_PQ_param,
         pf_generic_gens_para,

         edges_fbus,
         edges_tbus,
         edges_r,
         edges_x,
         edges_b,
         edges_ratio,
         edges_angle,
         edges_type,
         shunt_idx,
         Gs,
         Bs,

         Ynet_wt_nodes_idx_wt_adjacent_nodes,
         Ybr_cal_and_edges_orientation,

         kwd_sta_sta_ΔPQ_sol_by_json,
         generic_red_sol_kwd_para,
         generic_dyn_sol_kwd_para,

         plants_init_kwd_para,
         plants_cb_paras_kwd_para,

         y_aug_kw_para,

         with_faults,

         edges_r_x_b_ratio_angle_idx,

         Ynet_rows_Idxs_in_flattend,
         Ynet_real_imag_Idxs_in_flattend ) =
             NamedTupleTools.select(
                 dynamic_parameters,
                 (:states_Idx_syms_wt_functions,
         
                  :plants_cb_paras_switches,
                  :cb_states,

                  :baseMVA,
                  :basekV,

                  :loc_load_exist,

                  :slack_gens_nodes_idx,
                  :non_slack_gens_nodes_idx,
                  :gens_nodes_idx,
                  :non_gens_nodes_idx,
                  :gens_nodes_with_loc_loads_idx,
                  :all_nodes_idx,            
                  :nodes_with_demands_idx,

                  :dyn_pf_fun_kwd_n2s_idxs,
                  :dyn_pf_fun_kwd_net_idxs,

                  :generic_gens_para,
                  :generic_avrs_para,
                  :generic_govs_para,

                  :ode_gens_para,

                  :sta_pf_PQ_para,
                  :pf_PQ_param,
                  :pf_generic_gens_para,

                  :edges_fbus,
                  :edges_tbus,
                  :edges_r,
                  :edges_x,
                  :edges_b,
                  :edges_ratio,
                  :edges_angle,
                  :edges_type,
                  :shunt_idx,
                  :Gs,
                  :Bs,

                  :Ynet_wt_nodes_idx_wt_adjacent_nodes,
                  :Ybr_cal_and_edges_orientation,

                  :kwd_sta_sta_ΔPQ_sol_by_json,
                  :generic_red_sol_kwd_para,
                  :generic_dyn_sol_kwd_para,

                  :plants_init_kwd_para,
                  :plants_cb_paras_kwd_para,

                  :y_aug_kw_para,

                  :with_faults,

                  :edges_r_x_b_ratio_angle_idx,

                  :Ynet_rows_Idxs_in_flattend,
                  :Ynet_real_imag_Idxs_in_flattend ))

    
    (;

     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     pf_vh_θh_idx_and_idx2Idx,

     dyn_pf_flat_vh_flat_θh_Idx,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     scale_Pg_Png_Qng_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx) = 
         NamedTupleTools.select(
             states_Idx_syms_wt_functions,

             (:Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,

              :pf_vh_θh_idx_and_idx2Idx,
              
              :dyn_pf_flat_vh_flat_θh_Idx,
              
              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :scale_Pg_Png_Qng_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx))
    
    #############################

    #############################

    (non_gens_vh_idx,
     non_slack_gens_θh_idx,
     non_gens_θh_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     red_non_gens_vh_Idxs,
     red_non_slack_gens_θh_Idxs,
     red_non_gens_θh_Idxs,

     red_slack_value_Idxs,

     gens_nodes_idx,

     non_gens_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,

     n2s_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_non_gens_idx,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             pf_vh_θh_idx_and_idx2Idx,
             (:non_gens_vh_idx,
              :non_slack_gens_θh_idx,
              :non_gens_θh_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :red_non_gens_vh_Idxs,
              :red_non_slack_gens_θh_Idxs,
              :red_non_gens_θh_Idxs,

              :red_slack_value_Idxs,

              :gens_nodes_idx,

              :non_gens_nodes_idx,
              :non_slack_gens_and_non_gens_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,

              :n2s_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_non_gens_idx,
              :n2s_all_nodes_idx))

    (;Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         NamedTupleTools.select(
             Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx))

    (;edges_Ybr_cal,
     edges_orientation) =
         NamedTupleTools.select(
             Ybr_cal_and_edges_orientation,
             (:edges_Ybr_cal,
              :edges_orientation ))

    #--------------------------------------------
    #--------------------------------------------

    (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             Pg_Qg_Png_Qng_Pll_Qll_Idx,
             (:dyn_P_gens_Idxs,
              :dyn_Q_gens_Idxs,
              :dyn_P_non_gens_Idxs,
              :dyn_Q_non_gens_Idxs,
              :dyn_P_gens_loc_load_Idxs,
              :dyn_Q_gens_loc_load_Idxs))

    #--------------------------------------------

    red_vh_θh_idx =
        getproperty(
            pf_vh_θh_idx_and_idx2Idx,
            :red_vh_θh_idx)

    #--------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend ) =
         NamedTupleTools.select(
             get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                 Ynet),
             (:Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend))

    (;Ynet_rows_Idxs_in_flattend,
    Ynet_real_imag_Idxs_in_flattend ) =
        NamedTupleTools.select(
            get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                getproperty(
                    Ynet_wt_nodes_idx_wt_adjacent_nodes,
                    :Ynet)),
            (:Ynet_rows_Idxs_in_flattend,
             :Ynet_real_imag_Idxs_in_flattend) )

    (;Ynet_real_Idxs,
     Ynet_imag_Idxs) =
        NamedTupleTools.select(
            Ynet_real_imag_Idxs_in_flattend,
            (:Ynet_real_Idxs,
             :Ynet_imag_Idxs))

    #--------------------------------------------


    (;r_Idxs, x_Idxs, b_Idxs,
     ratio_Idxs, angle_Idxs) =
         NamedTupleTools.select(
             edges_r_x_b_ratio_angle_idx,
             (:r_Idxs, :x_Idxs, :b_Idxs,
              :ratio_Idxs, :angle_Idxs))

    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))

    #--------------------------------------------
    #--------------------------------------------
    
    (Pg,
     Qg,
     Png,
     Qng, 
     Pll,
     Qll) =
        NamedTupleTools.select(
            sta_pf_PQ_para,
            (:P_gens,
             :Q_gens,
             :P_non_gens,
             :Q_non_gens,
             :P_g_loc_load,
             :Q_g_loc_load ) )

    
    Pg_Qg_Png_Qng_Pll_Qll =
        [Pg;
         Qg;
         Png;
         Qng; 
         Pll;
         Qll ]
    
    #--------------------------------------------
    #--------------------------------------------

    Pg_inj_Qg_inj_Png_Qng_Pll_Qll =
        get_Pg_inj_Qg_inj_Png_Qng_Pll_Qll(
            Pg_Qg_Png_Qng_Pll_Qll;
            loc_load_exist,
            Pg_Qg_Png_Qng_Pll_Qll_Idx,    
            gens_nodes_idx,
            n2s_gens_idx,
            gens_nodes_with_loc_loads_idx,
            n2s_gens_with_loc_load_idxs)

    Pg_inj_Qg_inj_Png_Qng =
        getproperty(
            Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
            :Pg_inj_Qg_inj_Png_Qng)

    (;gens_Pg_inj,
     gens_Qg_inj,
     P_non_gens,
     Q_non_gens,
     Pg_inj_Qg_inj_Png_Qng) =
         NamedTupleTools.select(
             Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
             (:gens_Pg_inj,
              :gens_Qg_inj,
              :P_non_gens,
              :Q_non_gens,
              :Pg_inj_Qg_inj_Png_Qng))

    #--------------------------------------------
    #--------------------------------------------

    # Pg =
    #     Pg_Qg_Png_Qng_Pll_Qll[dyn_P_gens_Idxs]

    # Qg =
    #     Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_gens_Idxs]

    # Png =
    #     Pg_Qg_Png_Qng_Pll_Qll[dyn_P_non_gens_Idxs]

    # Qng =
    #     Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_non_gens_Idxs]

    # Pll =
    #     Pg_Qg_Png_Qng_Pll_Qll[dyn_P_gens_loc_load_Idxs]

    # Qll =
    #     Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_gens_loc_load_Idxs]


    #--------------------------------------------
    # Loss participation factors
    #--------------------------------------------

    gens_Sn =
        getproperty(ode_gens_para, :Sn)

    gens_loss_participation =
        get_loss_particpation_by_gens_rating(
            gens_Sn, Pg)

    gens_loss_participation_by_loading =
        get_loss_particpation_by_gens_loading( Pg)

    #--------------------------------------------
    # power disturbance resolution participation
    #--------------------------------------------

    gens_active_power_particpation_by_rating =
        get_gens_active_power_particpation_by_rating(
            gens_Sn, Pg)

    gens_active_power_particpation_by_loading =
        get_gens_active_power_particpation_by_loading(
            Pg )

    gens_reactive_power_particpation_by_rating =
        get_gens_reactive_power_particpation_by_rating(
            gens_Sn, Qg)

    gens_reactive_power_particpation_by_loading =
        get_gens_reactive_power_particpation_by_loading(
            Qg)


    active_power_disturbance_resolution_participation =
        gens_active_power_particpation_by_rating

    reactive_power_disturbance_resolution_participation =
        gens_reactive_power_particpation_by_loading

    #--------------------------------------------

    # slack_gens_vh_θh_gens_vh_non_slack_gens_vh =
    #     get_slack_gens_vh_θh_gens_vh_non_slack_gens_vh(
    #         vh,
    #         θh,
    #         dyn_pf_fun_kwd_net_idxs,
    #         dyn_pf_fun_kwd_n2s_idxs)

    slack_gens_vh_θh_gens_vh_non_slack_gens_vh  =
        gens_vh_slack_θh_para
    
    #--------------------------------------------
    # kwd parameters
    #--------------------------------------------

    pf_model_kwd_para =
        (;loc_load_exist,
         slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         gens_loss_participation,
         active_power_disturbance_resolution_participation,
         reactive_power_disturbance_resolution_participation,

         sta_pf_PQ_para,     

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         Pg_Qg_Png_Qng_Pll_Qll_Idx,
         Pg_Png_Qng_Idx,
         pf_vh_θh_idx_and_idx2Idx,

         scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
         scale_Pg_Png_Qng_Idx,
         dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

         shunt_idx,
         edges_fbus,
         edges_tbus,

         edges_type,
         edges_orientation,

         edges_r_x_b_ratio_angle_idx,
         Ynet_rows_Idxs_in_flattend,
         Ynet_real_imag_Idxs_in_flattend,

         baseMVA = baseMVA,
         basekV = basekV )

    generic_sol_kwd_para =
        (;
         Ybr_cal_and_edges_orientation,
         ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para )


    vh_θh =
        [ones(length( all_nodes_idx ));
         zeros(length(all_nodes_idx))]


    red_vh_θh_slack_value =
        [vh_θh[non_gens_vh_idx];
         vh_θh[non_slack_gens_θh_idx];
         vh_θh[non_gens_θh_idx];
         [0.0]]


    slack_gens_nodes_idx =
        getproperty(
            dyn_pf_fun_kwd_net_idxs,
            :slack_gens_nodes_idx)[1]

    # Subtract net loss from slack node Pg,
    # loss will subsequently be distributed

    if wt_network_loss_by_sta_pf_PQ_para_bool == true

        ds_gens_Pg_inj = @set gens_Pg_inj[
            slack_gens_nodes_idx]=
            gens_Pg_inj[slack_gens_nodes_idx] -
            get_total_P_network_loss_by_sta_pf_PQ_para(
                sta_pf_PQ_para,
                loc_load_exist )
        
    else
        
        ds_gens_Pg_inj = @set gens_Pg_inj[
            slack_gens_nodes_idx]=
            gens_Pg_inj[slack_gens_nodes_idx] -
            get_total_P_network_loss(
                vh,
                θh,
                Ynet;
                nodes_idx_with_adjacent_nodes_idx,
                n2s_all_nodes_idx,
                all_nodes_idx )
    end
    


    ds_Pg_inj_Png_Qng =
        [ds_gens_Pg_inj;
         P_non_gens;
         Q_non_gens]

    pf_fun_mismatch =
        get_red_model_distributed_slack_pf_ΔPQ_mismatch!


    red_sol = NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( (x, p) ->
                pf_fun_mismatch(
                    x, p;
                    pf_model_kwd_para =
                        pf_model_kwd_para )),
            red_vh_θh_slack_value,
            ds_Pg_inj_Png_Qng ), pf_alg,
        abstol = abstol,
        reltol = reltol)

    retcode = red_sol.retcode
    stats = red_sol.stats

    results_pf_red_sol_u =
        get_generic_results_ds_pf_red_sol_u(
            red_sol;
            generic_sol_kwd_para =
                generic_sol_kwd_para )

    results_conti_or_ds_pf_red_sol_u =
        get_generic_results_conti_or_ds_pf_red_sol_u(
            red_sol,
            ds_Pg_inj_Png_Qng;
            generic_sol_kwd_para =
                generic_sol_kwd_para)

    ds_slack_value =
        getproperty(
            results_pf_red_sol_u,
            :slack_value)

    ds_vh =
        getproperty(
            results_pf_red_sol_u,
            :vh)

    ds_θh =
        getproperty(
            results_pf_red_sol_u,
            :θh)

    ds_θh_deg =
        getproperty(
            results_pf_red_sol_u,
            :θh_deg)

    ds_pf_P_gens =
        getproperty(
            results_pf_red_sol_u,
            :pf_P_gens)

    ds_pf_Q_gens =
        getproperty(
            results_pf_red_sol_u,
            :pf_Q_gens)

    ds_pf_P_g_gens =
        getproperty(
            results_pf_red_sol_u,
            :pf_P_g_gens)

    ds_pf_Q_g_gens =
        getproperty(
            results_pf_red_sol_u,
            :pf_Q_g_gens)

    transformed_gens_nodes_idx =
        [n2s_all_nodes_idx[idx]
         for idx in gens_nodes_idx ]


    transformed_gens_nodes_with_loc_loads_idx =
        [n2s_all_nodes_idx[idx]
         for idx in gens_nodes_with_loc_loads_idx ]
    

    transformed_non_gens_nodes_idx =
        [ n2s_all_nodes_idx[idx]
         for idx in non_gens_nodes_idx]
    
    df_gens  = DataFrames.DataFrame(
        ;Bus = gens_nodes_idx,
        vh = round.(ds_vh[transformed_gens_nodes_idx];
                    digits = fractional_digits),
        θh = round.(ds_θh[transformed_gens_nodes_idx];
                    digits = fractional_digits),
        θh_deg =
            round.(ds_θh_deg[transformed_gens_nodes_idx];
                   digits = fractional_digits),
        Pg = round.(ds_pf_P_gens;
                   digits = fractional_digits),
        Qg = round.(ds_pf_Q_gens;
                   digits = fractional_digits),)

    
    df_gens_with_loc_load  = DataFrames.DataFrame(
        ;Bus = gens_nodes_with_loc_loads_idx,
        vh = round.(ds_vh[transformed_gens_nodes_with_loc_loads_idx];
                    digits = fractional_digits),
        θh = round.(ds_θh[transformed_gens_nodes_with_loc_loads_idx];
                    digits = fractional_digits),
        θh_deg =
            round.(ds_θh_deg[transformed_gens_nodes_with_loc_loads_idx];
                   digits = fractional_digits),
        Pll = round.(Pll;
                   digits = fractional_digits),
        Qll = round.(Qll;
                   digits = fractional_digits),)

    
    df_non_gens = DataFrames.DataFrame(
        ;Bus = non_gens_nodes_idx,
        vh =round.(ds_vh[transformed_non_gens_nodes_idx];
                   digits = fractional_digits),
        θh = round.(ds_θh[transformed_non_gens_nodes_idx];
                    digits = fractional_digits),
        θh_deg = round.(ds_θh_deg[
            transformed_non_gens_nodes_idx];
                        digits = fractional_digits),
        Png = round.(P_non_gens;
                   digits = fractional_digits),
        Qng = round.(Q_non_gens;
                   digits = fractional_digits),)
    
    df_all_nodes_result = DataFrames.DataFrame(
        ;Bus = all_nodes_idx,
        vh = round.(ds_vh;digits=fractional_digits),
        θh = round.(ds_θh;digits=fractional_digits),
        θh_deg = round.(ds_θh_deg;digits=fractional_digits),
        P = round.( [idx ∈ gens_nodes_idx ? ds_pf_P_gens[ n2s_gens_idx[idx] ] : P_non_gens[ n2s_non_gens_idx[ idx ] ] for idx in all_nodes_idx ]; digits = fractional_digits),
        Q = round.( [idx ∈ gens_nodes_idx ? ds_pf_Q_gens[ n2s_gens_idx[idx] ] : Q_non_gens[ n2s_non_gens_idx[ idx ] ] for idx in all_nodes_idx ]; digits = fractional_digits),)


    return loc_load_exist == true ?
        (;ds_slack_value,
         ds_vh,
         ds_θh,
         ds_θh_deg,
         ds_pf_P_gens,
         ds_pf_Q_gens,
         ds_pf_P_g_gens,
         ds_pf_Q_g_gens,

         df_gens,
         df_gens_with_loc_load,
         df_non_gens,
         df_all_nodes_result,

         results_pf_red_sol_u,
         results_conti_or_ds_pf_red_sol_u,
         retcode,

         stats,
         gens_loss_participation,
         gens_loss_participation_by_loading,

         gens_active_power_particpation_by_rating,
         gens_active_power_particpation_by_loading,
         gens_reactive_power_particpation_by_rating,
         gens_reactive_power_particpation_by_loading,

         active_power_disturbance_resolution_participation,
         reactive_power_disturbance_resolution_participation) :
             (;ds_slack_value,
              ds_vh,
              ds_θh,
              ds_θh_deg,
              ds_pf_P_gens,
              ds_pf_Q_gens,
              ds_pf_P_g_gens,
              ds_pf_Q_g_gens,

              df_gens,
              df_non_gens,
              df_all_nodes_result,

              results_pf_red_sol_u,
              results_conti_or_ds_pf_red_sol_u,
              retcode,
              stats,
              gens_loss_participation,
              gens_loss_participation_by_loading,

              gens_active_power_particpation_by_rating,
              gens_active_power_particpation_by_loading,
              gens_reactive_power_particpation_by_rating,
              gens_reactive_power_particpation_by_loading,

              active_power_disturbance_resolution_participation,
              reactive_power_disturbance_resolution_participation )
    
end



function sim_red_model_distributed_slack_pf(
    ;system_status=
        :pre_fault_state,
    case_name=
       "",

    with_faults = false,

    json_net_data_by_components_file =
        nothing,

    data_dir=
       "",

    components_libs_dir =
        "",        

    timespan   = timespan,

    on_fault_time = 5.0,
    clear_fault_time = 7.0,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [1],

    list_edges_to_have_fault = [ 8 ],
    clear_fault_selection_list = [1],

    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    fractional_digits=4,)

    
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

    #--------------------------------------    

    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end
    
    #--------------------------------------

    json_case_dir =
        joinpath(
            data_dir,
            "converted-data",
            case_name,
            "json")

    if  (json_net_data_by_components_file == "" ||
        json_net_data_by_components_file == nothing) 

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

    #--------------------------------------

    (baseMVA,
     basekV,

     vh,
     θh,
     u0_model_states_init,
     model_syms,
     model_mass_matrix,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     Pg_Qg_Png_Qng_Pll_Qll_Idx,
     Png_Qng_Pll_Qll_Idx,
     Pg_Png_Qng_Idx,

     pf_vh_θh_idx_and_idx2Idx,

     dyn_pf_flat_vh_flat_θh_Idx,

     edges_fbus,
     edges_tbus,

     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type,
     shunt_idx,
     Gs,
     Bs,

     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     Ybr_cal_and_edges_orientation,

     Pg_Qg_Png_Qng_Pll_Qll,
     loc_load_exist,
     slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

     edges_r_x_b_ratio_angle_idx,
     Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend,

     scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
     scale_Pg_Png_Qng_Idx,
     dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

     sta_pf_PQ_para,
     generic_gens_para,
     ode_gens_para ) =
         NamedTupleTools.select(
             getproperty(
            get_a_status_steady_state_data(
            system_status;
            with_faults = false,
            net_data_by_components_file =
                net_data_by_components_file,

                components_libs_dir =
                    components_libs_dir,        

                timespan   =
                    timespan,

                on_fault_time =
                    on_fault_time,
                clear_fault_time =
                    clear_fault_time,

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

            basekV = basekV ,    
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu),
            :static_prefault_paras),             
             (:baseMVA,
              :basekV,

              :vh,
              :θh,
              :u0_model_states_init,
              :model_syms,
              :model_mass_matrix,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_wt_fault_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :Png_Qng_Pll_Qll_Idx,
              :Pg_Png_Qng_Idx,

              :pf_vh_θh_idx_and_idx2Idx,

              :dyn_pf_flat_vh_flat_θh_Idx,

              :edges_fbus,
              :edges_tbus,

              :edges_r,
              :edges_x,
              :edges_b,
              :edges_ratio,
              :edges_angle,
              :edges_type,
              :shunt_idx,
              :Gs,
              :Bs,

              :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :Ybr_cal_and_edges_orientation,

              :Pg_Qg_Png_Qng_Pll_Qll,
              :loc_load_exist,

              :slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

              :edges_r_x_b_ratio_angle_idx,
              :Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend,

              :scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
              :scale_Pg_Png_Qng_Idx,
              :dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

              :sta_pf_PQ_para,
              :generic_gens_para,
              :ode_gens_para))

    (non_gens_vh_idx,
     non_slack_gens_θh_idx,
     non_gens_θh_idx,

     red_vh_Idxs,
     red_non_slack_gens_θh_idx2Idx,
     red_non_gens_θh_idx2Idx,

     red_non_gens_vh_Idxs,
     red_non_slack_gens_θh_Idxs,
     red_non_gens_θh_Idxs,

     red_slack_value_Idxs,

     gens_nodes_idx,

     non_gens_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx,

     n2s_gens_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_non_gens_idx,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             pf_vh_θh_idx_and_idx2Idx,
             (:non_gens_vh_idx,
              :non_slack_gens_θh_idx,
              :non_gens_θh_idx,

              :red_vh_Idxs,
              :red_non_slack_gens_θh_idx2Idx,
              :red_non_gens_θh_idx2Idx,

              :red_non_gens_vh_Idxs,
              :red_non_slack_gens_θh_Idxs,
              :red_non_gens_θh_Idxs,

              :red_slack_value_Idxs,

              :gens_nodes_idx,

              :non_gens_nodes_idx,
              :non_slack_gens_and_non_gens_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx,

              :n2s_gens_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_non_gens_idx,
              :n2s_all_nodes_idx))

    (;Ynet,
     nodes_idx_with_adjacent_nodes_idx) =
         NamedTupleTools.select(
             Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx))

    (;edges_Ybr_cal,
     edges_orientation) =
         NamedTupleTools.select(
             Ybr_cal_and_edges_orientation,
             (:edges_Ybr_cal,
              :edges_orientation ))

    #--------------------------------------------
    #--------------------------------------------

    (;dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs) =
         NamedTupleTools.select(
             Pg_Qg_Png_Qng_Pll_Qll_Idx,
             (:dyn_P_gens_Idxs,
              :dyn_Q_gens_Idxs,
              :dyn_P_non_gens_Idxs,
              :dyn_Q_non_gens_Idxs,
              :dyn_P_gens_loc_load_Idxs,
              :dyn_Q_gens_loc_load_Idxs))

    #--------------------------------------------

    red_vh_θh_idx =
        getproperty(
            pf_vh_θh_idx_and_idx2Idx,
            :red_vh_θh_idx)

    #--------------------------------------------

    (;Ynet_rows_Idxs_in_flattend,
     Ynet_real_imag_Idxs_in_flattend ) =
         NamedTupleTools.select(
             get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                 Ynet),
             (:Ynet_rows_Idxs_in_flattend,
              :Ynet_real_imag_Idxs_in_flattend))

    (;Ynet_rows_Idxs_in_flattend,
    Ynet_real_imag_Idxs_in_flattend ) =
        NamedTupleTools.select(
            get_Ynet_real_imag_Idxs_wt_rows_Idxs_in_flattend(
                getproperty(
                    Ynet_wt_nodes_idx_wt_adjacent_nodes,
                    :Ynet)),
            (:Ynet_rows_Idxs_in_flattend,
             :Ynet_real_imag_Idxs_in_flattend) )

    (;Ynet_real_Idxs,
     Ynet_imag_Idxs) =
        NamedTupleTools.select(
            Ynet_real_imag_Idxs_in_flattend,
            (:Ynet_real_Idxs,
             :Ynet_imag_Idxs))

    #--------------------------------------------


    (;r_Idxs, x_Idxs, b_Idxs,
     ratio_Idxs, angle_Idxs) =
         NamedTupleTools.select(
             edges_r_x_b_ratio_angle_idx,
             (:r_Idxs, :x_Idxs, :b_Idxs,
              :ratio_Idxs, :angle_Idxs))

    (vh_Idxs, θh_Idxs) =
        NamedTupleTools.select(
            dyn_pf_flat_vh_flat_θh_Idx,
            (:dyn_pf_vh_Idxs,
             :dyn_pf_θh_Idxs))

    #--------------------------------------------
    #--------------------------------------------

    Pg_inj_Qg_inj_Png_Qng_Pll_Qll =
        get_Pg_inj_Qg_inj_Png_Qng_Pll_Qll(
            Pg_Qg_Png_Qng_Pll_Qll;
            loc_load_exist,
            Pg_Qg_Png_Qng_Pll_Qll_Idx,    
            gens_nodes_idx,
            n2s_gens_idx,
            n2s_gens_with_loc_load_idxs)

    Pg_inj_Qg_inj_Png_Qng =
        getproperty(
            Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
            :Pg_inj_Qg_inj_Png_Qng)

    (;gens_Pg_inj,
     gens_Qg_inj,
     P_non_gens,
     Q_non_gens,
     Pg_inj_Qg_inj_Png_Qng) =
         NamedTupleTools.select(
             Pg_inj_Qg_inj_Png_Qng_Pll_Qll,
             (:gens_Pg_inj,
              :gens_Qg_inj,
              :P_non_gens,
              :Q_non_gens,
              :Pg_inj_Qg_inj_Png_Qng))

    #--------------------------------------------
    #--------------------------------------------

    Pg =
        Pg_Qg_Png_Qng_Pll_Qll[dyn_P_gens_Idxs]

    Qg =
        Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_gens_Idxs]

    Png =
        Pg_Qg_Png_Qng_Pll_Qll[dyn_P_non_gens_Idxs]

    Qng =
        Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_non_gens_Idxs]

    Pll =
        Pg_Qg_Png_Qng_Pll_Qll[dyn_P_gens_loc_load_Idxs]

    Qll =
        Pg_Qg_Png_Qng_Pll_Qll[dyn_Q_gens_loc_load_Idxs]


    #--------------------------------------------
    # Loss participation factors
    #--------------------------------------------

    gens_Sn =
        getproperty(ode_gens_para, :Sn)

    gens_loss_participation =
        get_loss_particpation_by_gens_rating(
            gens_Sn, Pg)

    gens_loss_participation_by_loading =
        get_loss_particpation_by_gens_loading( Pg)

    #--------------------------------------------
    # power disturbance resolution participation
    #--------------------------------------------

    gens_active_power_particpation_by_rating =
        get_gens_active_power_particpation_by_rating(
            gens_Sn, Pg)

    gens_active_power_particpation_by_loading =
        get_gens_active_power_particpation_by_loading(
            Pg )

    gens_reactive_power_particpation_by_rating =
        get_gens_reactive_power_particpation_by_rating(
            gens_Sn, Qg)

    gens_reactive_power_particpation_by_loading =
        get_gens_reactive_power_particpation_by_loading(
            Qg)


    active_power_disturbance_resolution_participation =
        gens_active_power_particpation_by_rating

    reactive_power_disturbance_resolution_participation =
        gens_reactive_power_particpation_by_loading

    #--------------------------------------------
    # kwd parameters
    #--------------------------------------------

    pf_model_kwd_para =
        (;loc_load_exist,
         slack_gens_vh_θh_gens_vh_non_slack_gens_vh,

         Ynet_wt_nodes_idx_wt_adjacent_nodes,

         gens_loss_participation,
         active_power_disturbance_resolution_participation,
         reactive_power_disturbance_resolution_participation,

         sta_pf_PQ_para,     

         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         Pg_Qg_Png_Qng_Pll_Qll_Idx,
         Pg_Png_Qng_Idx,
         pf_vh_θh_idx_and_idx2Idx,

         scale_Pg_Qg_Png_Qng_Pll_Qll_Idx,
         scale_Pg_Png_Qng_Idx,
         dyn_pf_flat_vh_flat_θh_wt_slack_value_Idx,

         shunt_idx,
         edges_fbus,
         edges_tbus,

         edges_type,
         edges_orientation,

         edges_r_x_b_ratio_angle_idx,
         Ynet_rows_Idxs_in_flattend,
         Ynet_real_imag_Idxs_in_flattend,

         baseMVA = baseMVA,
         basekV = basekV )

    generic_sol_kwd_para =
        (;
         Ybr_cal_and_edges_orientation,
         ode_gens_para,
         sta_pf_PQ_para,
         pf_model_kwd_para )


    vh_θh =
        [ones(length( vh ));
         zeros(length(θh))]


    red_vh_θh_slack_value =
        [vh_θh[non_gens_vh_idx];
         vh_θh[non_slack_gens_θh_idx];
         vh_θh[non_gens_θh_idx];
         [0.0]]


    slack_gens_nodes_idx =
        getproperty(
            dyn_pf_fun_kwd_net_idxs,
            :slack_gens_nodes_idx)[1]

    # Subtract net loss from slack node Pg,
    # loss will subsequently be distributed


    ds_gens_Pg_inj = @set gens_Pg_inj[slack_gens_nodes_idx]=
        gens_Pg_inj[slack_gens_nodes_idx] -
        get_total_P_network_loss(
            vh,
            θh,
            Ynet;
            nodes_idx_with_adjacent_nodes_idx,
            n2s_all_nodes_idx,
            all_nodes_idx )


    ds_Pg_inj_Png_Qng =
        [ds_gens_Pg_inj;
         P_non_gens;
         Q_non_gens]

    pf_fun_mismatch =
        get_red_model_distributed_slack_pf_ΔPQ_mismatch!


    red_sol = NonlinearSolve.solve(
        NonlinearProblem(
            NonlinearFunction( (x, p) ->
                pf_fun_mismatch(
                    x, p;
                    pf_model_kwd_para =
                        pf_model_kwd_para )),
            red_vh_θh_slack_value,
            ds_Pg_inj_Png_Qng ), pf_alg )

    results_pf_red_sol_u =
        get_generic_results_ds_pf_red_sol_u(
            red_sol;
            generic_sol_kwd_para =
                generic_sol_kwd_para )

    results_conti_or_ds_pf_red_sol_u =
        get_generic_results_conti_or_ds_pf_red_sol_u(
            red_sol,
            ds_Pg_inj_Png_Qng;
            generic_sol_kwd_para =
                generic_sol_kwd_para)

    ds_slack_value =
        getproperty(
            results_pf_red_sol_u,
            :slack_value)

    ds_vh =
        getproperty(
            results_pf_red_sol_u,
            :vh)

    ds_θh =
        getproperty(
            results_pf_red_sol_u,
            :θh)

    ds_θh_deg =
        getproperty(
            results_pf_red_sol_u,
            :θh_deg)

    ds_pf_P_gens =
        getproperty(
            results_pf_red_sol_u,
            :pf_P_gens)

    ds_pf_Q_gens =
        getproperty(
            results_pf_red_sol_u,
            :pf_Q_gens)

    ds_pf_P_g_gens =
        getproperty(
            results_pf_red_sol_u,
            :pf_P_g_gens)

    ds_pf_Q_g_gens =
        getproperty(
            results_pf_red_sol_u,
            :pf_Q_g_gens)


    transformed_gens_nodes_idx =
        [n2s_gens_idx[idx]
         for idx in gens_nodes_idx ]


    transformed_non_gens_nodes_idx =
        [ n2s_non_gens_idx[idx]
         for idx in non_gens_nodes_idx]
    
    df_gens  = DataFrames.DataFrame(
        ;Bus = gens_nodes_idx,
        vh = round.(ds_vh[transformed_gens_nodes_idx];
                    digits = fractional_digits),
        θh = round.(ds_θh[transformed_gens_nodes_idx];
                    digits = fractional_digits),
        θh_deg =
            round.(ds_θh_deg[transformed_gens_nodes_idx];
                   digits = fractional_digits),
        Pg = round.(ds_pf_P_gens;
                   digits = fractional_digits),
        Qg = round.(ds_pf_Q_gens;
                   digits = fractional_digits),)

    df_non_gens = DataFrames.DataFrame(
        ;Bus = non_gens_nodes_idx,
        vh =round.(ds_vh[transformed_non_gens_nodes_idx];
                   digits = fractional_digits),
        θh = round.(ds_θh[transformed_non_gens_nodes_idx];
                    digits = fractional_digits),
        θh_deg = round.(ds_θh_deg[
            transformed_non_gens_nodes_idx];
                        digits = fractional_digits),
        Png = round.(P_non_gens;
                   digits = fractional_digits),
        Qng = round.(Q_non_gens;
                   digits = fractional_digits),)
    
    df_all_nodes_result = DataFrames.DataFrame(
        ;Bus = all_nodes_idx,
        vh = round.(ds_vh;digits=fractional_digits),
        θh = round.(ds_θh;digits=fractional_digits),
        θh_deg = round.(ds_θh_deg;digits=fractional_digits),
        P = round.( [idx ∈ gens_nodes_idx ? ds_pf_P_gens[ n2s_gens_idx[idx] ] : P_non_gens[ n2s_non_gens_idx[ idx ] ] for idx in all_nodes_idx ]; digits = fractional_digits),
        Q = round.( [idx ∈ gens_nodes_idx ? ds_pf_Q_gens[ n2s_gens_idx[idx] ] : Q_non_gens[ n2s_non_gens_idx[ idx ] ] for idx in all_nodes_idx ]; digits = fractional_digits),)

    return (;ds_slack_value,
            ds_vh,
            ds_θh,
            ds_θh_deg,
            ds_pf_P_gens,
            ds_pf_Q_gens,
            ds_pf_P_g_gens,
            ds_pf_Q_g_gens,
            
            df_gens,
            df_non_gens,
            df_all_nodes_result,
            
            results_pf_red_sol_u,
            results_conti_or_ds_pf_red_sol_u)
    
end



function sim_model_by_mass_matrix_ode_by_model_dynamics!(
    # system_state = :pre_fault_state
    system_state ;
    with_faults = false,
    net_data_by_components_file =
        nothing,
    components_libs_dir =
        nothing,        

    timespan         = 20,    
    on_fault_time    = 5.0,
    clear_fault_time = 7.0,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [1],

    list_edges_to_have_fault = [ 8 ],
    clear_fault_selection_list = [1],

    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    
    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    dt            = 0.01,
    Δt            = 1.0 / 2^(4) )

    #---------------------------------------------------

    sim_timespan = (0.0, timespan)

    #---------------------------------------------------
    
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

    # #--------------------------------------    

    # if (data_dir == "") || (data_dir == nothing)

    #     package_dir = pkgdir(ePowerSim)

    #     data_dir = joinpath(package_dir,
    #                  "data")
        
    # end
    
    # #--------------------------------------

    # json_case_dir =
    #     joinpath(
    #         data_dir,
    #         "converted-data",
    #         case_name,
    #         "json")

    # if  (json_net_data_by_components_file == "" ||
    #     json_net_data_by_components_file == nothing) 

    #     net_data_by_components_file =
    #         joinpath(
    #             json_case_dir,
    #             "net_data_by_components_file.json")
    # else

    #     net_data_by_components_file =
    #         joinpath(
    #             json_case_dir,
    #             json_net_data_by_components_file)

    # end

    #--------------------------------------

    (u0_model_states_init,
     model_syms,
     model_mass_matrix,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states) =
         NamedTupleTools.select(
                     getproperty(
                         get_a_status_steady_state_data(
                             system_state;
                              with_faults,
                              net_data_by_components_file,
                              components_libs_dir,        

                              timespan,

                              on_fault_time,
                              clear_fault_time,

                              list_fault_point_from_node_a,
                              list_fault_resistance,
                              list_no_line_circuit,

                              list_edges_to_have_fault,
                              clear_fault_selection_list,

                              basekV,    
                              use_pu_in_PQ,
                              line_data_in_pu),
                         :static_prefault_paras) ,             
             (:u0_model_states_init,
              :model_syms,
              :model_mass_matrix,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_wt_fault_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states ))

    #----------------------------------------
    # ODE system dyamanics simulation
    #----------------------------------------    

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll
    
    model_dynamics_para =
        (;generic_model_dynamics_para,
         plants_cb_paras_switches )
        

    system_dynamics_fun! =
        mm_ode_generic_system_model_dynamics!

    #----------------------------------------

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para

    #----------------------------------------

    system_sol =
        DifferentialEquations.solve(
            ODEProblem(
        ODEFunction(
        (dx,x,p,t) ->
            system_dynamics_fun!(
                dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms,
    mass_matrix = model_mass_matrix ) ,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para ),
            ode_alg,
            callback = cb_states,
            abstol = abstol,
            reltol = reltol )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan)    
    
end


function sim_model_by_mass_matrix_ode_by_funcs_dynamics!(
    # system_state = :pre_fault_state
    system_state ;
    with_faults = false,
    net_data_by_components_file =
        nothing,
    components_libs_dir =
        nothing,        

    timespan         = 20,    
    on_fault_time    = 5.0,
    clear_fault_time = 7.0,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [1],

    list_edges_to_have_fault = [ 8 ],
    clear_fault_selection_list = [1],

    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    
    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    dt            = 0.01,
    Δt            = 1.0 / 2^(4) )

    #---------------------------------------------------

    sim_timespan = (0.0, timespan)

    #---------------------------------------------------
    
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

    #--------------------------------------    

    (u0_model_states_init,
     model_syms,
     model_mass_matrix,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states) =
         NamedTupleTools.select(
                     getproperty(
                         get_a_status_steady_state_data(
                             system_state;
                              with_faults,
                              net_data_by_components_file,
                              components_libs_dir,        

                              timespan,

                              on_fault_time,
                              clear_fault_time,

                              list_fault_point_from_node_a,
                              list_fault_resistance,
                              list_no_line_circuit,

                              list_edges_to_have_fault,
                              clear_fault_selection_list,

                              basekV,    
                              use_pu_in_PQ,
                              line_data_in_pu),
                         :static_prefault_paras) ,             
             (:u0_model_states_init,
              :model_syms,
              :model_mass_matrix,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_wt_fault_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states ))

    #----------------------------------------
    # ODE system dyamanics simulation
    #----------------------------------------    

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll
    
    model_dynamics_para =
        (;generic_model_dynamics_para,
      plants_cb_paras_switches )
        

    system_dynamics_fun! =
        mm_ode_generic_system_model_by_funcs_dynamics!

    #----------------------------------------

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para

    #----------------------------------------

    system_sol =
        DifferentialEquations.solve(
            ODEProblem(
        ODEFunction(
        (dx,x,p,t) ->
            system_dynamics_fun!(
                dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms,
    mass_matrix = model_mass_matrix ) ,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para ),
            ode_alg,
            callback = cb_states,
            abstol = abstol,
            reltol = reltol )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan)    
    
end



function sim_model_by_mass_matrix_by_ode_pf_funcs!(
    # system_state = :pre_fault_state
    system_state ;
    with_faults = false,
    net_data_by_components_file =
        nothing,
    components_libs_dir =
        nothing,        

    timespan         = 20,    
    on_fault_time    = 5.0,
    clear_fault_time = 7.0,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [1],

    list_edges_to_have_fault = [ 8 ],
    clear_fault_selection_list = [1],

    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    
    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    dt            = 0.01,
    Δt            = 1.0 / 2^(4) )

    #---------------------------------------------------

    sim_timespan = (0.0, timespan)

    #---------------------------------------------------
        
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

    #--------------------------------------    

    (u0_model_states_init,
     model_syms,
     model_mass_matrix,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states) =
         NamedTupleTools.select(
                     getproperty(
                         get_a_status_steady_state_data(
                             system_state;
                              with_faults,
                              net_data_by_components_file,
                              components_libs_dir,        

                              timespan,

                              on_fault_time,
                              clear_fault_time,

                              list_fault_point_from_node_a,
                              list_fault_resistance,
                              list_no_line_circuit,

                              list_edges_to_have_fault,
                              clear_fault_selection_list,

                              basekV,    
                              use_pu_in_PQ,
                              line_data_in_pu),
                         :static_prefault_paras) ,             
             (:u0_model_states_init,
              :model_syms,
              :model_mass_matrix,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_wt_fault_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states ))

    #----------------------------------------
    # ODE system dyamanics simulation
    #----------------------------------------    

    model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    system_dynamics_fun! =
        mm_ode_generic_system_dynamics_by_ode_pf_funcs!

    #----------------------------------------

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para

    #----------------------------------------

    system_sol =
        DifferentialEquations.solve(
            ODEProblem(
        ODEFunction(
        (dx,x,p,t) ->
            system_dynamics_fun!(
                dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms,
    mass_matrix = model_mass_matrix ) ,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para ),
            ode_alg,
            abstol = abstol,
            reltol = reltol )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan)    
    
end

#---------------------------------------------------

function sim_model_by_dae_pf_funcs!(
    # system_state = :pre_fault_state
    system_state ;
    with_faults = false,
    net_data_by_components_file =
        nothing,
    components_libs_dir =
        nothing,        

    timespan         = 20,    
    on_fault_time    = 5.0,
    clear_fault_time = 7.0,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [1],

    list_edges_to_have_fault = [ 8 ],
    clear_fault_selection_list = [1],

    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    
    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    dt            = 0.01,
    Δt            = 1.0 / 2^(4) )

    #---------------------------------------------------

    sim_timespan = (0.0, timespan)

    #---------------------------------------------------
        
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

    #--------------------------------------    

    (u0_model_states_init,
     model_syms,
     model_bool_dae_vars,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states) =
         NamedTupleTools.select(
                     getproperty(
                         get_a_status_steady_state_data(
                             system_state;
                              with_faults,
                              net_data_by_components_file,
                              components_libs_dir,        

                              timespan,

                              on_fault_time,
                              clear_fault_time,

                              list_fault_point_from_node_a,
                              list_fault_resistance,
                              list_no_line_circuit,

                              list_edges_to_have_fault,
                              clear_fault_selection_list,

                              basekV,    
                              use_pu_in_PQ,
                              line_data_in_pu),
                         :static_prefault_paras),
             
             (:u0_model_states_init,
              :model_syms,
              :model_bool_dae_vars,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_wt_fault_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states))

    #----------------------------------------    

    du0_model_states_init =
        zeros( length( u0_model_states_init ))

    res = similar( u0_model_states_init )

    #----------------------------------------
    # DAE system dyamanics simulation
    #----------------------------------------    

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para

    model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    system_dynamics_fun! =
        dae_generic_system_dynamics_by_dae_pf_funcs!

    
    #----------------------------------------

    # generic_model_dynamics_para =
    #     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    # model_dynamics_para =
    #     (; generic_model_dynamics_para,
    #      plants_cb_paras_switches )        
    
    # dae_generic_system_model_by_funcs_dynamics!
    # dae_generic_system_model_dynamics!

    #----------------------------------------
    
    system_sol =
        DifferentialEquations.solve(
            DAEProblem(
        DAEFunction(
            (res, dx, x, p, t) ->
                system_dynamics_fun!(
                    res, dx, x,
                    model_dynamics_para,
                    t;
                    kwd_para =
                        model_dynamics_kwd_para);
            syms =
                model_syms),
                du0_model_states_init,
                u0_model_states_init,
                sim_timespan,
                model_dynamics_para,
                differential_vars =
                    model_bool_dae_vars),
            dae_alg,            
            callback = cb_states,            
            abstol = abstol,
            reltol = reltol )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan)
end

#---------------------------------------------------

function sim_model_by_dae_funcs_dynamics!(
    # system_state = :pre_fault_state
    system_state ;
    with_faults = false,
    net_data_by_components_file =
        nothing,
    components_libs_dir =
        nothing,        

    timespan         = 20,    
    on_fault_time    = 5.0,
    clear_fault_time = 7.0,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [1],

    list_edges_to_have_fault = [ 8 ],
    clear_fault_selection_list = [1],

    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    
    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    dt            = 0.01,
    Δt            = 1.0 / 2^(4) )

    #---------------------------------------------------

    sim_timespan = (0.0, timespan)

    #---------------------------------------------------
        
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

    #--------------------------------------    

    (u0_model_states_init,
     model_syms,
     model_bool_dae_vars,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states) =
         NamedTupleTools.select(
                     getproperty(
                         get_a_status_steady_state_data(
                             system_state;
                              with_faults,
                              net_data_by_components_file,
                              components_libs_dir,        

                              timespan,

                              on_fault_time,
                              clear_fault_time,

                              list_fault_point_from_node_a,
                              list_fault_resistance,
                              list_no_line_circuit,

                              list_edges_to_have_fault,
                              clear_fault_selection_list,

                              basekV,    
                              use_pu_in_PQ,
                              line_data_in_pu),
                         :static_prefault_paras) ,   
             (:u0_model_states_init,
              :model_syms,
              :model_bool_dae_vars,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_wt_fault_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states))

    #----------------------------------------    

    du0_model_states_init =
        zeros( length( u0_model_states_init ))

    res = similar( u0_model_states_init )

    #----------------------------------------
    # DAE system dyamanics simulation
    #----------------------------------------    

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para

    # model_dynamics_para =
    #     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    # system_dynamics_fun! =
    #     dae_generic_system_dynamics_by_dae_pf_funcs!

    
    #----------------------------------------

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    model_dynamics_para =
        (; generic_model_dynamics_para,
         plants_cb_paras_switches )        
    
    system_dynamics_fun! =
        dae_generic_system_model_by_funcs_dynamics!
    
    # dae_generic_system_model_dynamics!

    #----------------------------------------
    
    system_sol =
        DifferentialEquations.solve(
            DAEProblem(
        DAEFunction(
            (res, dx, x, p, t) ->
                system_dynamics_fun!(
                    res, dx, x,
                    model_dynamics_para,
                    t;
                    kwd_para =
                        model_dynamics_kwd_para);
            syms =
                model_syms),
                du0_model_states_init,
                u0_model_states_init,
                sim_timespan,
                model_dynamics_para,
                differential_vars =
                    model_bool_dae_vars),
            dae_alg,            
            callback = cb_states,            
            abstol = abstol,
            reltol = reltol )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan)
end


#---------------------------------------------------

function sim_model_by_dae_model_dynamics!(
    # system_state = :pre_fault_state
    system_state ;
    with_faults = false,
    net_data_by_components_file =
        nothing,
    components_libs_dir =
        nothing,        

    timespan         = 20,    
    on_fault_time    = 5.0,
    clear_fault_time = 7.0,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit =  [1],

    list_edges_to_have_fault = [ 8 ],
    clear_fault_selection_list = [1],

    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    
    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    dt            = 0.01,
    Δt            = 1.0 / 2^(4) )

    #---------------------------------------------------

    sim_timespan = (0.0, timespan)

    #---------------------------------------------------
        
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

    #--------------------------------------    

    (u0_model_states_init,
     model_syms,
     model_bool_dae_vars,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states) =
         NamedTupleTools.select(
                     getproperty(
                         get_a_status_steady_state_data(
                             system_state;
                              with_faults,
                              net_data_by_components_file,
                              components_libs_dir,        

                              timespan,

                              on_fault_time,
                              clear_fault_time,

                              list_fault_point_from_node_a,
                              list_fault_resistance,
                              list_no_line_circuit,

                              list_edges_to_have_fault,
                              clear_fault_selection_list,

                              basekV,    
                              use_pu_in_PQ,
                              line_data_in_pu),
                         :static_prefault_paras) ,      
             (:u0_model_states_init,
              :model_syms,
              :model_bool_dae_vars,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_wt_fault_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states))

    #----------------------------------------    

    du0_model_states_init =
        zeros( length( u0_model_states_init ))

    res = similar( u0_model_states_init )

    #----------------------------------------
    # DAE system dyamanics simulation
    #----------------------------------------    

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para
    
    #----------------------------------------

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    model_dynamics_para =
        (; generic_model_dynamics_para,
         plants_cb_paras_switches )        
    
    system_dynamics_fun! =
        dae_generic_system_model_dynamics!

    #----------------------------------------
    
    system_sol =
        DifferentialEquations.solve(
            DAEProblem(
        DAEFunction(
            (res, dx, x, p, t) ->
                system_dynamics_fun!(
                    res, dx, x,
                    model_dynamics_para,
                    t;
                    kwd_para =
                        model_dynamics_kwd_para);
            syms =
                model_syms),
                du0_model_states_init,
                u0_model_states_init,
                sim_timespan,
                model_dynamics_para,
                differential_vars =
                    model_bool_dae_vars),
            dae_alg,            
            callback = cb_states,            
            abstol = abstol,
            reltol = reltol )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan)
end

#---------------------------------------------------
 
function sim_sudden_load_change(
    net_data_by_components_file;

    target_bus_name = "bus4",
    pertubation_factor = 1.10,
    restoration_factor = 1.0,

    pertubation_time = 10.0,
    restoration_time = 10.5,
    
    Δt1 = 1.5,
    Δt2 = 1.5,
    
    
    sim_type,

    case_name,
    components_libs_dir="",
    data_dir="",
    json_case_dir="",
    
    results_dir,
    figure_dir,
    sd_dynamics_sim_csv_filename,

    # # base setting and some booleans
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    timespan      = 20,    
    dt            = 0.01,
    Δt            = 1.0 / 2^(4),

    # # digits rounding
    fractional_digits = 6 )

    #---------------------------------------------------
    # Simulation Period
    #---------------------------------------------------

    time_start    = 0.0

    time_final    = timespan

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

    #--------------------------------------    

    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end


    if (json_case_dir == "") || (json_case_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
    
        json_case_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name,
                     "json")
        
    end

    #---------------------------------------------------
    ## variables, parameters, indices, intialisation
    #---------------------------------------------------

    (;u0_model_states_init,
     model_syms,
     model_bool_dae_vars,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states,

     plants_states_syms,

     gens_nodes_idx,

     state_labels,
     algebraic_vars_labels,
     network_vars_labels,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_system_simulation_parameters(
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
                 reltol = reltol),

             (:u0_model_states_init,
              :model_syms,
              :model_bool_dae_vars,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states,

              :plants_states_syms,
              :gens_nodes_idx,

              :state_labels,
              :algebraic_vars_labels,
              :network_vars_labels,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #----------------------------------------


    (;loc_load_exist,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_labels,
     algebraic_vars_labels) =
         NamedTupleTools.select(
             generic_system_dynamics_kwd_para ,
             (:loc_load_exist,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,

              :state_labels,
              :algebraic_vars_labels))

    #---------------------------------------------------
    ## ODE system dyamanics simulation
    #---------------------------------------------------

    system_dynamics_fun! =
        dae_generic_system_model_by_funcs_dynamics!

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll
        # deepcopy(ω_ref_v_ref_p_order_Png_Qng_Pll_Qll)

    model_dynamics_para =
        (;generic_model_dynamics_para,
          plants_cb_paras_switches )          

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para

    #----------------------------------------
    #----------------------------------------

    du0_model_states_init =
        zeros(length(u0_model_states_init))

    res = similar(u0_model_states_init)

    #----------------------------------------
    # integrator
    #----------------------------------------    

    system_integrator =
        DifferentialEquations.init(
            DAEProblem(
        DAEFunction(
        (res, dx, x, p, t) ->
            system_dynamics_fun!(
                res, dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms),
        du0_model_states_init,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para,
        differential_vars =
            model_bool_dae_vars ),
            dae_alg,
            dt = dt,
            callback = cb_states,
            tstops = [time_final],
            advance_to_tstop = true )

    #---------------------------------------------------
    # parameters df header
    #---------------------------------------------------

    generic_model_dynamics_para_df_header_sym =
        get_make_df_header_generic_model_dynamics_para(
            loc_load_exist,
            dyn_pf_fun_kwd_net_idxs)

    parameter_df =
        DataFrame(
            OrderedDict(a_header => Float64[]
                for a_header in
                    generic_model_dynamics_para_df_header_sym))


    #---------------------------------------------------
    # parameters to be perturbed
    #---------------------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_loads_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx))

    bus_name = target_bus_name
    
    bus_no_or_bus_name = bus_name

    P_or_Q_or_Pll_or_Qll_sym = :P

    var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            P_or_Q_or_Pll_or_Qll_sym;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    P_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :P;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    Q_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :Q;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    vh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :vh ] )[1]

    θh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ;network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :θh ] )[1]

    #---------------------------------------------------

    target_parameter_sym =
        Symbol(
            "$(bus_name)_$(P_or_Q_or_Pll_or_Qll_sym)")

    change_in_parameter_dict =
        OrderedDict(
            :t => Float64[],
            target_parameter_sym =>
                    Float64[] )

    #---------------------------------------------------
    # Save parameters base value
    #---------------------------------------------------

    # Get the Power at bus 4

    var_normal_value =
        system_integrator.p.generic_model_dynamics_para[
            var_idx]
    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    #----------------------------------------
    # simulation steps
    #----------------------------------------    

    tstop1 = pertubation_time             # time_final/10
    tstop2 = restoration_time             # time_final/8
    tstop3 = restoration_time + Δt1       # time_final/6
    tstop4 = restoration_time + Δt1 + Δt2 # time_final/2
    tstop5 = time_final/1.0    

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 1
    # -------------------------------------

    add_tstop!(system_integrator,
               tstop1 )

    step!(system_integrator)

    push!(parameter_df,
        tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    system_sol =
        system_integrator.sol

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 2
    # -------------------------------------

    # A pertubation of Active power at bus 4

    pertubation_factor =
        pertubation_factor

    pertubation_tstop = tstop2

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 3
    # -------------------------------------

    # Bring back to normal the Power at bus 4

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop3

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 4
    # -------------------------------------

    # Continute simulation 

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop4

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym  )

    system_sol = system_integrator.sol


    push!( parameter_df, tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    #---------------------------------------------------
    # simulate tstop4 till the end
    #---------------------------------------------------

    DifferentialEquations.solve!(system_integrator)

    system_sol = system_integrator.sol

    #---------------------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)


    final_intg_model_dynamics_para =
        getproperty(
            system_integrator.p,
            :generic_model_dynamics_para)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          final_intg_model_dynamics_para[var_idx] )

    #---------------------------------------------------
    # Save results to files
    #---------------------------------------------------
    

    sd_dynamics_sim_df = DataFrame(system_sol)

    sd_dynamics_sim_df[!, :] =
        round.(
            sd_dynamics_sim_df[:, :],
            digits=fractional_digits)

    sd_dynamics_sim_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
                     "$(bus_name)-" *
                     "$(sim_type)-states.csv")

    CSV.write(sd_dynamics_sim_csv_filename,
              sd_dynamics_sim_df )

    #---------------------------------------------------

    parameter_df[!, :] =
        round.(
            parameter_df[:, :],
            digits=fractional_digits)

    sd_dynamics_para_sim_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
                     "$(bus_no_or_bus_name)-" *
                     "$(sim_type)-parameters.csv")

    CSV.write(sd_dynamics_para_sim_csv_filename,
              parameter_df )

    #---------------------------------------------------

    generic_model_dynamics_para_df_header_idx =
        1:length(generic_model_dynamics_para_df_header_sym)

    paras_df_header_idx_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(bus_no_or_bus_name)-" *
                     "$(sim_type)-para-df-header-idx.csv")

    CSV.write(paras_df_header_idx_csv_filename,
              DataFrame(OrderedDict(
                  :idx =>
                      generic_model_dynamics_para_df_header_idx,
                  :parameters =>
                      generic_model_dynamics_para_df_header_sym)) )

    #---------------------------------------------------

    generic_model_dynamics_state_df_header_sym =
        [[:t];model_syms]

    generic_model_dynamics_state_df_header_idx =
        1:length( generic_model_dynamics_state_df_header_sym )

    states_df_header_idx_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(bus_no_or_bus_name)-" *
                     "$(sim_type)-states-df-header-idx.csv")

    CSV.write(states_df_header_idx_csv_filename,
              DataFrame(OrderedDict(
                  :idx =>
                      generic_model_dynamics_state_df_header_idx,
                  :states =>
                      generic_model_dynamics_state_df_header_sym) ) )

    #---------------------------------------------------

    save_pertubation_stage_plot(
        case_name,
        system_sol;
        model_syms,
        gens_nodes_names,
        SM_gens_nodes_names,
        non_gens_nodes_names,
        sim_timespan,
        figure_dir,
        P_or_Q_or_Pll_or_Qll_sym,
        bus_idx = bus_no_or_bus_name)

    #---------------------------------------------------
    # extract solution auxilliary results
    #---------------------------------------------------

    sol_auxilliary_results =
        get_sol_auxilliary_results(
            system_sol;
            state_labels,
            algebraic_vars_labels,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    auxilliary_results_plot =
        make_plot_gens_streamedlined_auxilliary_results(;
            sol_auxilliary_results... )


    names_vars_plots =
        propertynames(auxilliary_results_plot)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                auxilliary_results_plot,
                a_vars_plots)

        local filename =
            "$(case_name)-" *
            "$(bus_no_or_bus_name)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "$(String(a_vars_plots)).pdf"

        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))

    end

    #---------------------------------------------------

    target_parameter_plot =
        plot(change_in_parameter_dict[:t],
             change_in_parameter_dict[target_parameter_sym],
             linetype=:steppre,
             yminorticks = 10,
             xminorticks = 10,
             fmt = :pdf,
             lw = 1,
             xlabel = "t [s]",
             ylabel = "$(target_parameter_sym) [p.u]",
             labels = "$(target_parameter_sym)",
             bottom_margin=2Plots.mm,  # Adjust bottom margin
             left_margin=5Plots.mm,   # Adjust left margin
             right_margin=2Plots.mm,  # Adjust right margin
             top_margin=2Plots.mm )

    target_parameter_plot_filename =
            "$(case_name)-" *
            "$(bus_no_or_bus_name)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "target-parameter-" *
            "$(String(target_parameter_sym)).pdf"

    savefig(target_parameter_plot,
                joinpath(
                    figure_dir,
                    target_parameter_plot_filename))


    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan,

            sol_auxilliary_results,
            target_parameter_plot)

end



function sim_sudden_load_change(    
    pertubation_factor,
    restoration_factor,
    pertubation_time,
    restoration_time;
    
    Δt1 = 1.5,
    Δt2 = 1.5,
    
    net_data_by_components_file,

    target_bus_name = "bus4",
    # pertubation_factor = 1.10,
    # restoration_factor = 1.0,
    
    sim_type,

    components_libs_dir="",

    # # base setting and some booleans
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    timespan      = 20,    
    dt            = 0.01,
    Δt            = 1.0 / 2^(4),

    # # digits rounding
    fractional_digits = 6 )

    #---------------------------------------------------
    # Simulation Period
    #---------------------------------------------------

    time_start    = 0.0

    time_final    = timespan

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

    #---------------------------------------------------
    ## variables, parameters, indices, intialisation
    #---------------------------------------------------

    (;u0_model_states_init,
     model_syms,
     model_bool_dae_vars,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states,

     plants_states_syms,

     gens_nodes_idx,

     state_labels,
     algebraic_vars_labels,
     network_vars_labels,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_system_simulation_parameters(
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
                 reltol = reltol),

             (:u0_model_states_init,
              :model_syms,
              :model_bool_dae_vars,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states,

              :plants_states_syms,
              :gens_nodes_idx,

              :state_labels,
              :algebraic_vars_labels,
              :network_vars_labels,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #----------------------------------------


    (;loc_load_exist,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_labels,
     algebraic_vars_labels) =
         NamedTupleTools.select(
             generic_system_dynamics_kwd_para ,
             (:loc_load_exist,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,

              :state_labels,
              :algebraic_vars_labels))

    #---------------------------------------------------
    ## ODE system dyamanics simulation
    #---------------------------------------------------

    system_dynamics_fun! =
        dae_generic_system_model_by_funcs_dynamics!

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll
        # deepcopy(ω_ref_v_ref_p_order_Png_Qng_Pll_Qll)

    model_dynamics_para =
        (;generic_model_dynamics_para,
          plants_cb_paras_switches )          

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para

    #----------------------------------------
    #----------------------------------------

    du0_model_states_init =
        zeros(length(u0_model_states_init))

    res = similar(u0_model_states_init)

    #----------------------------------------
    # integrator
    #----------------------------------------    

    system_integrator =
        DifferentialEquations.init(
            DAEProblem(
        DAEFunction(
        (res, dx, x, p, t) ->
            system_dynamics_fun!(
                res, dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms),
        du0_model_states_init,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para,
        differential_vars =
            model_bool_dae_vars ),
            dae_alg,
            dt = dt,
            callback = cb_states,
            tstops = [time_final],
            advance_to_tstop = true )

    #---------------------------------------------------
    # parameters df header
    #---------------------------------------------------

    generic_model_dynamics_para_df_header_sym =
        get_make_df_header_generic_model_dynamics_para(
            loc_load_exist,
            dyn_pf_fun_kwd_net_idxs)

    parameter_df =
        DataFrame(
            OrderedDict(a_header => Float64[]
                for a_header in
                    generic_model_dynamics_para_df_header_sym  ))


    #---------------------------------------------------
    # parameters to be perturbed
    #---------------------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_loads_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx))

    bus_name = target_bus_name
    
    bus_no_or_bus_name = bus_name

    P_or_Q_or_Pll_or_Qll_sym = :P

    var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            P_or_Q_or_Pll_or_Qll_sym;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    P_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :P;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    Q_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :Q;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    vh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :vh ] )[1]

    θh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :θh ] )[1]

    #---------------------------------------------------

    target_parameter_sym =
        Symbol(
            "$(bus_name)_$(P_or_Q_or_Pll_or_Qll_sym)")

    change_in_parameter_dict =
        OrderedDict(
            :t => Float64[],
            target_parameter_sym =>
                    Float64[] )

    #---------------------------------------------------
    # Save parameters base value
    #---------------------------------------------------

    # Get the Power at bus 4

    var_normal_value =
        system_integrator.p.generic_model_dynamics_para[
            var_idx]
    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    #----------------------------------------
    # simulation steps
    #----------------------------------------    

    tstop1 = pertubation_time             # time_final/10
    tstop2 = restoration_time             # time_final/8
    tstop3 = restoration_time + Δt1       # time_final/6
    tstop4 = restoration_time + Δt1 + Δt2 # time_final/2
    tstop5 = time_final/1.0    

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 1
    # -------------------------------------

    add_tstop!(system_integrator,
               tstop1 )

    step!(system_integrator)

    push!(parameter_df,
        tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    system_sol =
        system_integrator.sol

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 2
    # -------------------------------------

    # A pertubation of Active power at bus 4

    pertubation_factor =
        pertubation_factor

    pertubation_tstop = tstop2

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 3
    # -------------------------------------

    # Bring back to normal the Power at bus 4

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop3

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 4
    # -------------------------------------

    # Continute simulation 

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop4

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym  )

    system_sol = system_integrator.sol


    push!( parameter_df, tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    #---------------------------------------------------
    # simulate tstop4 till the end
    #---------------------------------------------------

    DifferentialEquations.solve!(system_integrator)

    system_sol = system_integrator.sol

    #---------------------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    final_intg_model_dynamics_para =
        getproperty(
            system_integrator.p,
            :generic_model_dynamics_para)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          final_intg_model_dynamics_para[var_idx] )

    #---------------------------------------------------
    # extract solution auxilliary results
    #---------------------------------------------------

    sol_auxilliary_results =
        get_sol_auxilliary_results(
            system_sol;
            state_labels,
            algebraic_vars_labels,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    auxilliary_results_plot =
        make_plot_gens_streamedlined_auxilliary_results(
            ; sol_auxilliary_results... )
    
    #---------------------------------------------------

    target_parameter_plot =
        plot(change_in_parameter_dict[:t],
             change_in_parameter_dict[target_parameter_sym],
             linetype=:steppre,
             yminorticks = 10,
             xminorticks = 10,
             fmt = :pdf,
             lw = 1,
             xlabel = "t [s]",
             ylabel = "$(target_parameter_sym) [p.u]",
             labels = "$(target_parameter_sym)",
             bottom_margin=2Plots.mm,  # Adjust bottom margin
             left_margin=5Plots.mm,   # Adjust left margin
             right_margin=2Plots.mm,  # Adjust right margin
             top_margin=2Plots.mm )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan,

            sol_auxilliary_results,
            auxilliary_results_plot,
            target_parameter_plot)

end

#---------------------------------------------------
 
function sim_mm_sudden_load_change(
    net_data_by_components_file;

    target_bus_name = "bus4",
    pertubation_factor = 1.10,
    restoration_factor = 1.0,

    pertubation_time = 10.0,
    restoration_time = 10.5,
    
    Δt1 = 1.5,
    Δt2 = 1.5,
    
    
    sim_type,

    case_name,
    components_libs_dir="",
    data_dir="",
    json_case_dir="",
    
    results_dir,
    figure_dir,
    sd_dynamics_sim_csv_filename,

    # # base setting and some booleans
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    timespan      = 20,    
    dt            = 0.01,
    Δt            = 1.0 / 2^(4),

    # # digits rounding
    fractional_digits = 6 )
    
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

    #--------------------------------------    

    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end


    if (json_case_dir == "") || (json_case_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
    
        json_case_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name,
                     "json")
        
    end

    #---------------------------------------------------
    # Simulation Period
    #---------------------------------------------------

    time_start    = 0.0

    time_final    = timespan

    tspan         = (0.0, timespan)

    sim_timespan  = (0.0, timespan)

    plot_timespan = (0.0, timespan)

    #---------------------------------------------------
    ## variables, parameters, indices, intialisation
    #---------------------------------------------------

    (;u0_model_states_init,
     model_syms,
     model_bool_dae_vars,
     model_mass_matrix,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states,

     plants_states_syms,

     gens_nodes_idx,

     state_labels,
     algebraic_vars_labels,
     network_vars_labels,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_system_simulation_parameters(
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
                 reltol = reltol),

             (:u0_model_states_init,
              :model_syms,
              :model_bool_dae_vars,
              :model_mass_matrix,
              

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states,

              :plants_states_syms,
              :gens_nodes_idx,

              :state_labels,
              :algebraic_vars_labels,
              :network_vars_labels,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #----------------------------------------


    (;loc_load_exist,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_labels,
     algebraic_vars_labels) =
         NamedTupleTools.select(
             generic_system_dynamics_kwd_para ,
             (:loc_load_exist,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,

              :state_labels,
              :algebraic_vars_labels))
    
    #----------------------------------------
    # ODE system dyamanics simulation
    #----------------------------------------    

    system_dynamics_fun! =
        mm_ode_generic_system_model_by_funcs_dynamics!
        # mm_ode_generic_system_dynamics_by_ode_pf_funcs!

    #----------------------------------------

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para


    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    model_dynamics_para =
        (;generic_model_dynamics_para,
          plants_cb_paras_switches )          

    # model_dynamics_para =
    #     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    #----------------------------------------
    # integrator
    #----------------------------------------    

    system_integrator =
        DifferentialEquations.init(
            ODEProblem(
        ODEFunction(
        (dx,x,p,t) ->
            system_dynamics_fun!(
                dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms,
    mass_matrix = model_mass_matrix ) ,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para ),
            ode_alg,
            # dt = dt,
            
            callback = cb_states,
            tstops = [time_final],
            advance_to_tstop = true,
            
            abstol = abstol,
            reltol = reltol )
    
    #---------------------------------------------------
    # parameters df header
    #---------------------------------------------------

    generic_model_dynamics_para_df_header_sym =
        get_make_df_header_generic_model_dynamics_para(
            loc_load_exist,
            dyn_pf_fun_kwd_net_idxs)

    parameter_df =
        DataFrame(
            OrderedDict(a_header => Float64[]
                for a_header in
                    generic_model_dynamics_para_df_header_sym))


    #---------------------------------------------------
    # parameters to be perturbed
    #---------------------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_loads_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx))

    bus_name = target_bus_name
    
    bus_no_or_bus_name = bus_name

    P_or_Q_or_Pll_or_Qll_sym = :P

    var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            P_or_Q_or_Pll_or_Qll_sym;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    P_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :P;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    Q_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :Q;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    vh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :vh ] )[1]

    θh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :θh ] )[1]

    #---------------------------------------------------

    target_parameter_sym =
        Symbol(
            "$(bus_name)_$(P_or_Q_or_Pll_or_Qll_sym)")

    change_in_parameter_dict =
        OrderedDict(
            :t => Float64[],
            target_parameter_sym =>
                    Float64[] )

    #---------------------------------------------------
    # Save parameters base value
    #---------------------------------------------------

    # Get the Power at bus 4

    var_normal_value =
        system_integrator.p.generic_model_dynamics_para[
            var_idx]
    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    #----------------------------------------
    # simulation steps
    #----------------------------------------    

    tstop1 = pertubation_time             # time_final/10
    tstop2 = restoration_time             # time_final/8
    tstop3 = restoration_time + Δt1       # time_final/6
    tstop4 = restoration_time + Δt1 + Δt2 # time_final/2
    tstop5 = time_final/1.0    

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 1
    # -------------------------------------

    add_tstop!(system_integrator,
               tstop1 )

    step!(system_integrator)

    push!(parameter_df,
        tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    system_sol =
        system_integrator.sol

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 2
    # -------------------------------------

    # A pertubation of Active power at bus 4

    pertubation_factor =
        pertubation_factor

    pertubation_tstop = tstop2

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 3
    # -------------------------------------

    # Bring back to normal the Power at bus 4

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop3

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 4
    # -------------------------------------

    # Continute simulation 

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop4

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym  )

    system_sol = system_integrator.sol


    push!( parameter_df, tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    #---------------------------------------------------
    # simulate tstop4 till the end
    #---------------------------------------------------

    DifferentialEquations.solve!(system_integrator)

    system_sol = system_integrator.sol

    #---------------------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)


    final_intg_model_dynamics_para =
        getproperty(
            system_integrator.p,
            :generic_model_dynamics_para)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          final_intg_model_dynamics_para[var_idx] )

    #---------------------------------------------------
    # Save results to files
    #---------------------------------------------------
    

    sd_dynamics_sim_df = DataFrame(system_sol)

    sd_dynamics_sim_df[!, :] =
        round.(
            sd_dynamics_sim_df[:, :],
            digits=fractional_digits)

    sd_dynamics_sim_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
                     "$(bus_name)-" *
                     "$(sim_type)-states.csv")

    CSV.write(sd_dynamics_sim_csv_filename,
              sd_dynamics_sim_df )

    #---------------------------------------------------

    parameter_df[!, :] =
        round.(
            parameter_df[:, :],
            digits=fractional_digits)

    sd_dynamics_para_sim_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
                     "$(bus_no_or_bus_name)-" *
                     "$(sim_type)-parameters.csv")

    CSV.write(sd_dynamics_para_sim_csv_filename,
              parameter_df )

    #---------------------------------------------------

    generic_model_dynamics_para_df_header_idx =
        1:length(generic_model_dynamics_para_df_header_sym)

    paras_df_header_idx_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(bus_no_or_bus_name)-" *
                     "$(sim_type)-para-df-header-idx.csv")

    CSV.write(paras_df_header_idx_csv_filename,
              DataFrame(OrderedDict(
                  :idx =>
                      generic_model_dynamics_para_df_header_idx,
                  :parameters =>
                      generic_model_dynamics_para_df_header_sym)) )

    #---------------------------------------------------

    generic_model_dynamics_state_df_header_sym =
        [[:t];model_syms]

    generic_model_dynamics_state_df_header_idx =
        1:length( generic_model_dynamics_state_df_header_sym )

    states_df_header_idx_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(bus_no_or_bus_name)-" *
                     "$(sim_type)-states-df-header-idx.csv")

    CSV.write(states_df_header_idx_csv_filename,
              DataFrame(OrderedDict(
                  :idx =>
                      generic_model_dynamics_state_df_header_idx,
                  :states =>
                      generic_model_dynamics_state_df_header_sym) ) )


    #---------------------------------------------------

    save_pertubation_stage_plot(
        case_name,
        system_sol;
        model_syms,
        gens_nodes_names,
        SM_gens_nodes_names,
        non_gens_nodes_names,
        sim_timespan,
        figure_dir,
        P_or_Q_or_Pll_or_Qll_sym,
        bus_idx = bus_no_or_bus_name)

    #---------------------------------------------------
    # extract solution auxilliary results
    #---------------------------------------------------

    sol_auxilliary_results =
        get_sol_auxilliary_results(
            system_sol;
            state_labels,
            algebraic_vars_labels,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    auxilliary_results_plot =
        make_plot_gens_streamedlined_auxilliary_results(;
            sol_auxilliary_results... )


    names_vars_plots =
        propertynames(auxilliary_results_plot)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                auxilliary_results_plot,
                a_vars_plots)

        local filename =
            "$(case_name)-" *
            "$(bus_no_or_bus_name)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "$(String(a_vars_plots)).pdf"

        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))

    end

    #---------------------------------------------------

    target_parameter_plot =
        plot(change_in_parameter_dict[:t],
             change_in_parameter_dict[target_parameter_sym],
             linetype=:steppre,
             yminorticks = 10,
             xminorticks = 10,
             fmt = :pdf,
             lw = 1,
             xlabel = "t [s]",
             ylabel = "$(target_parameter_sym) [p.u]",
             labels = "$(target_parameter_sym)",
             bottom_margin=2Plots.mm,  # Adjust bottom margin
             left_margin=5Plots.mm,   # Adjust left margin
             right_margin=2Plots.mm,  # Adjust right margin
             top_margin=2Plots.mm )

    target_parameter_plot_filename =
            "$(case_name)-" *
            "$(bus_no_or_bus_name)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "target-parameter-" *
            "$(String(target_parameter_sym)).pdf"

    savefig(target_parameter_plot,
                joinpath(
                    figure_dir,
                    target_parameter_plot_filename))


    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan,

            sol_auxilliary_results,
            target_parameter_plot)

end



function sim_mm_sudden_load_change(    
    pertubation_factor,
    restoration_factor,
    pertubation_time,
    restoration_time;
    
    Δt1 = 1.5,
    Δt2 = 1.5,
    
    net_data_by_components_file,

    target_bus_name = "bus4",
    # pertubation_factor = 1.10,
    # restoration_factor = 1.0,
    
    sim_type,

    components_libs_dir="",

    # # base setting and some booleans
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    timespan      = 20,    
    dt            = 0.01,
    Δt            = 1.0 / 2^(4),

    # # digits rounding
    fractional_digits = 6 )

    #---------------------------------------------------
    # Simulation Period
    #---------------------------------------------------

    time_start    = 0.0

    time_final    = timespan

    tspan         = (0.0, timespan)

    sim_timespan  = (0.0, timespan)

    plot_timespan = (0.0, timespan)

    #---------------------------------------------------
    ## variables, parameters, indices, intialisation
    #---------------------------------------------------
    
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

    #--------------------------------------    

    (;u0_model_states_init,
     model_syms,
     model_bool_dae_vars,
     model_mass_matrix,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states,

     plants_states_syms,

     gens_nodes_idx,

     state_labels,
     algebraic_vars_labels,
     network_vars_labels,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_system_simulation_parameters(
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
                 reltol = reltol),

             (:u0_model_states_init,
              :model_syms,
              :model_bool_dae_vars,
              :model_mass_matrix,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states,

              :plants_states_syms,
              :gens_nodes_idx,

              :state_labels,
              :algebraic_vars_labels,
              :network_vars_labels,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #----------------------------------------


    (;loc_load_exist,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_labels,
     algebraic_vars_labels) =
         NamedTupleTools.select(
             generic_system_dynamics_kwd_para ,
             (:loc_load_exist,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,

              :state_labels,
              :algebraic_vars_labels))

    #---------------------------------------------------
    ## ODE system dyamanics simulation
    #---------------------------------------------------

    system_dynamics_fun! =
        mm_ode_generic_system_model_by_funcs_dynamics!
        # mm_ode_generic_system_dynamics_by_ode_pf_funcs!

    #----------------------------------------

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para


    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    model_dynamics_para =
        (;generic_model_dynamics_para,
          plants_cb_paras_switches )          

    # model_dynamics_para =
    #     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    #----------------------------------------
    # integrator
    #----------------------------------------    

    system_integrator =
        DifferentialEquations.init(
            ODEProblem(
        ODEFunction(
        (dx,x,p,t) ->
            system_dynamics_fun!(
                dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms,
    mass_matrix = model_mass_matrix ) ,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para ),
            ode_alg,
            # dt = dt,
            
            callback = cb_states,
            tstops = [time_final],
            advance_to_tstop = true,
            
            abstol = abstol,
            reltol = reltol )
    
    #---------------------------------------------------
    # parameters df header
    #---------------------------------------------------

    generic_model_dynamics_para_df_header_sym =
        get_make_df_header_generic_model_dynamics_para(
            loc_load_exist,
            dyn_pf_fun_kwd_net_idxs)

    parameter_df =
        DataFrame(
            OrderedDict(a_header => Float64[]
                for a_header in
                    generic_model_dynamics_para_df_header_sym  ))


    #---------------------------------------------------
    # parameters to be perturbed
    #---------------------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_loads_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx))

    bus_name = target_bus_name
    
    bus_no_or_bus_name = bus_name

    P_or_Q_or_Pll_or_Qll_sym = :P

    var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            P_or_Q_or_Pll_or_Qll_sym;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    P_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :P;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    Q_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :Q;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    vh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :vh ] )[1]

    θh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :θh ] )[1]

    #---------------------------------------------------

    target_parameter_sym =
        Symbol(
            "$(bus_name)_$(P_or_Q_or_Pll_or_Qll_sym)")

    change_in_parameter_dict =
        OrderedDict(
            :t => Float64[],
            target_parameter_sym =>
                    Float64[] )

    #---------------------------------------------------
    # Save parameters base value
    #---------------------------------------------------

    # Get the Power at bus 4

    var_normal_value =
        system_integrator.p.generic_model_dynamics_para[
            var_idx]
    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    #----------------------------------------
    # simulation steps
    #----------------------------------------    

    tstop1 = pertubation_time             # time_final/10
    tstop2 = restoration_time             # time_final/8
    tstop3 = restoration_time + Δt1       # time_final/6
    tstop4 = restoration_time + Δt1 + Δt2 # time_final/2
    tstop5 = time_final/1.0    

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 1
    # -------------------------------------

    add_tstop!(system_integrator,
               tstop1 )

    step!(system_integrator)

    push!(parameter_df,
        tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    system_sol =
        system_integrator.sol

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 2
    # -------------------------------------

    # A pertubation of Active power at bus 4

    pertubation_factor =
        pertubation_factor

    pertubation_tstop = tstop2

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 3
    # -------------------------------------

    # Bring back to normal the Power at bus 4

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop3

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 4
    # -------------------------------------

    # Continute simulation 

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop4

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym  )

    system_sol = system_integrator.sol


    push!( parameter_df, tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    #---------------------------------------------------
    # simulate tstop4 till the end
    #---------------------------------------------------

    DifferentialEquations.solve!(system_integrator)

    system_sol = system_integrator.sol

    #---------------------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    final_intg_model_dynamics_para =
        getproperty(
            system_integrator.p,
            :generic_model_dynamics_para)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          final_intg_model_dynamics_para[var_idx] )

    #---------------------------------------------------
    # extract solution auxilliary results
    #---------------------------------------------------

    sol_auxilliary_results =
        get_sol_auxilliary_results(
            system_sol;
            state_labels,
            algebraic_vars_labels,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    auxilliary_results_plot =
        make_plot_gens_streamedlined_auxilliary_results(
            ; sol_auxilliary_results... )
    
    #---------------------------------------------------

    target_parameter_plot =
        plot(change_in_parameter_dict[:t],
             change_in_parameter_dict[target_parameter_sym],
             linetype=:steppre,
             yminorticks = 10,
             xminorticks = 10,
             fmt = :pdf,
             lw = 1,
             xlabel = "t [s]",
             ylabel = "$(target_parameter_sym) [p.u]",
             labels = "$(target_parameter_sym)",
             bottom_margin=2Plots.mm,  # Adjust bottom margin
             left_margin=5Plots.mm,   # Adjust left margin
             right_margin=2Plots.mm,  # Adjust right margin
             top_margin=2Plots.mm )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan,

            sol_auxilliary_results,
            auxilliary_results_plot,
            target_parameter_plot)

end

#---------------------------------------------------


function get_sim_mm_sudden_load_or_line_outage_pertubation(
    pertubation_factor,
    restoration_factor,
    pertubation_time,
    restoration_time;
    
    Δt1 = 1.5,
    Δt2 = 1.5,
    
    net_data_by_components_file,

    target_bus_name = "bus4",
    # pertubation_factor = 1.10,
    # restoration_factor = 1.0,
    
    sim_type,

    outage_type =
        :line_outage,
    
    on_fault_time,
    clear_fault_time,
    
    line_outage_time,
    generation_adjustment_time,
    
    components_libs_dir="",

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit  = [1],

    list_edges_to_have_fault   = [ 8 ],
    clear_fault_selection_list = [1],
        
    list_network_status = 
            [:pre_fault_state,
             :post_fault_state],
    
    with_faults =
        false,    
    
    # # base setting and some booleans
    
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    
    timespan      = 20,    
    dt            = 0.01,
    Δt            = 1.0 / 2^(4),

    # # digits rounding
    
    fractional_digits = 6 )
    
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

    #---------------------------------------------------
    # Simulation Period
    #---------------------------------------------------

    time_start    = 0.0

    time_final    = timespan

    tspan         = (0.0, timespan)

    sim_timespan  = (0.0, timespan)

    plot_timespan = (0.0, timespan)

    #---------------------------------------------------
    ## variables, parameters, indices, intialisation
    #---------------------------------------------------

    ntuple_status_steady_state_data =
        get_ntuple_status_steady_state_data(
            ;with_faults =
                with_faults,
            net_data_by_components_file =
                net_data_by_components_file,
            components_libs_dir =
                components_libs_dir,

            timespan =
                timespan,
            on_fault_time =
                on_fault_time,
            clear_fault_time =
                clear_fault_time,

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

            basekV =
                basekV,    
            use_pu_in_PQ =
                use_pu_in_PQ,
            line_data_in_pu =
                line_data_in_pu,
            list_network_status =
                list_network_status )
   
    #---------------------------------------------------

    (;loc_load_exist,
     state_labels,
     algebraic_vars_labels,
     
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     system_fault_status,
     generic_system_dynamics_wt_fault_kwd_para,
     on_fault_net_para,
     cleared_selected_lines_faults_net_para,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

     # model_bool_dae_vars_wt_fault,
     # model_syms_wt_fault,         
     # u0_model_states_init_wt_fault,

     model_bool_dae_vars,     
     model_syms,
     u0_model_states_init,
     model_mass_matrix,     
     
     cb_states,
     plants_cb_paras_switches,

     nodes_names,
     gens_nodes_names,
     non_gens_nodes_names,
     SM_gens_nodes_names,
     SC_gens_nodes_names,

     ωref0_vref0_porder0_id_iq_vh_Idx,
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type,
     Gs,
     Bs,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    ntuple_status_steady_state_data,
                    :pre_fault_state),
                :static_prefault_paras),        
            (:loc_load_exist,
             :state_labels,
             :algebraic_vars_labels,
             :dyn_pf_fun_kwd_n2s_idxs,
             :dyn_pf_fun_kwd_net_idxs,
             
             :system_fault_status,
             :generic_system_dynamics_wt_fault_kwd_para,
             :on_fault_net_para,
             :cleared_selected_lines_faults_net_para,

             :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

             # :model_bool_dae_vars_wt_fault,
             # :model_syms_wt_fault,         
             # :u0_model_states_init_wt_fault,

             :model_bool_dae_vars,     
             :model_syms,
             :u0_model_states_init,
             :model_mass_matrix,             

             :cb_states,
             :plants_cb_paras_switches,

             :nodes_names,
             :gens_nodes_names,
             :non_gens_nodes_names,
             :SM_gens_nodes_names,
             :SC_gens_nodes_names,

             :ωref0_vref0_porder0_id_iq_vh_Idx,
             :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

             :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

             :edges_r,
             :edges_x,
             :edges_b,
             :edges_ratio,
             :edges_angle,
             :edges_type,
             :Gs,
             :Bs,
             :Ybr_cal_and_edges_orientation,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes))


    #----------------------------------------

    # po := post_outage
    
    (po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    ntuple_status_steady_state_data,
                    :post_fault_state),
                :dynamic_status_paras),
            (:ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,))
    
    #----------------------------------------
    #----------------------------------------


    (Ynet, ) =
         NamedTupleTools.select(
        Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:Ynet, ) )

    #----------------------------------------

    (fault_Ynet,
     post_fault_Ynet) =
        NamedTupleTools.select(
        cleared_selected_lines_faults_net_para,
            (:pre_clear_fault_Ynet,
             :post_clear_fault_Ynet, ))

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

    #---------------------------------------

    @show system_fault_status
    
    if system_fault_status[1] != 0

        system_fault_status[1] = 0
        
    end
    
    #---------------------------------------
    #---------------------------------------
    # line loss with only porder_adj
    #---------------------------------------
    #---------------------------------------
        
    cb_line_outage = DiscreteCallback(
        (u, t, integrator) ->
            on_line_outage_condition(
                u, t, integrator,
                line_outage_time),

       on_line_outage_affect!;
        save_positions=(true, true),
        initializealg =
            ShampineCollocationInit() )

    cb_outage_set =
        CallbackSet(cb_line_outage,)

    tstop_outage =
        [line_outage_time]

    
    if outage_type == :line_outage_wt_pref_adjs

        gens_porder_adj =
            po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_p_order_Idx]

        cb_gens_porder_adjustment = DiscreteCallback(
            (u, t, integrator) ->
                on_generation_adjustment_condition(
                    u, t, integrator,
                    generation_adjustment_time),

            (integrator) ->
                on_generation_adjustment_affect!(
                    integrator,
                    gens_porder_adj,
                    dyn_p_order_Idx );
            save_positions=(true, true),
            initializealg =
                ShampineCollocationInit() )


        cb_outage_set =
            CallbackSet(cb_line_outage,
            cb_gens_porder_adjustment)

        tstop_outage =
            [line_outage_time,
             generation_adjustment_time]

        
    elseif outage_type == :line_outage_wt_vpref_adjs

        vref_and_porder_Idx =
            [ dyn_v_ref_Idx; dyn_p_order_Idx]

        vref_and_porder_adj =
            po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                vref_and_porder_Idx]
        
        cb_vref_and_porder_adj = DiscreteCallback(
            (u, t, integrator) ->
                on_generation_adjustment_condition(
                    u, t, integrator,
                    generation_adjustment_time),

            (integrator) ->
                on_generation_adjustment_affect!(
                integrator,
                vref_and_porder_adj,
                vref_and_porder_Idx);
            save_positions=(true, true),
            initializealg =
                ShampineCollocationInit() )

        cb_outage_set =
            CallbackSet(cb_line_outage,
            cb_vref_and_porder_adj)

        tstop_outage =
            [line_outage_time,
             generation_adjustment_time]
        
    else # :line_outage

        nothing
    end

    #----------------------------------------
    #---------------------------------------

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    model_dynamics_para =
        (;generic_model_dynamics_para,
         Ynet,
         post_fault_Ynet,
         system_fault_status,
         plants_cb_paras_switches )

    model_dynamics_kwd_para =
        generic_system_dynamics_wt_fault_kwd_para
    
    #----------------------------------------

    
    system_dynamics_fun! =
        mm_line_outage_generic_dynamics_wt_pre_post_fault_by_ode_pf_funcs!

    #----------------------------------------
    # integrator
    #----------------------------------------    

    # system_integrator =
    #     DifferentialEquations.solve(
    #         ODEProblem(
    #     ODEFunction(
    #     (dx,x,p,t) ->
    #         system_dynamics_fun!(
    #             dx, x,
    #             model_dynamics_para,
    #             t;
    #             kwd_para =
    #                 model_dynamics_kwd_para);
    #         syms =
    #             model_syms,
    #         mass_matrix =
    #             model_mass_matrix ) ,
    #             u0_model_states_init,
    #             sim_timespan,
    #             model_dynamics_para,
    #             callback =
    #                 cb_states),
    #             ode_alg,
    #             callback =
    #                 cb_outage_set,
    #             tstops =
    #                 tstop_outage,
    #             abstol = abstol,
    #             reltol = reltol )
    
    #----------------------------------------
    # integrator
    #----------------------------------------    

    system_integrator =
        DifferentialEquations.init(
            ODEProblem(
        ODEFunction(
        (dx,x,p,t) ->
            system_dynamics_fun!(
                dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms,
            mass_matrix =
                model_mass_matrix ) ,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para ),
            ode_alg,
            # dt = dt,
            
            callback = cb_states,
            tstops = [time_final],
            advance_to_tstop = true,
            
            abstol = abstol,
            reltol = reltol )
    
    #---------------------------------------------------
    # parameters df header
    #---------------------------------------------------

    generic_model_dynamics_para_df_header_sym =
        get_make_df_header_generic_model_dynamics_para(
            loc_load_exist,
            dyn_pf_fun_kwd_net_idxs)

    parameter_df =
        DataFrame(
            OrderedDict(a_header => Float64[]
                for a_header in
              generic_model_dynamics_para_df_header_sym))

    #---------------------------------------------------
    # parameters to be perturbed
    #---------------------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_loads_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx))

    bus_name =
        target_bus_name
    
    bus_no_or_bus_name =
        bus_name

    P_or_Q_or_Pll_or_Qll_sym =
        :P

    var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            P_or_Q_or_Pll_or_Qll_sym;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs)

    P_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :P;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs)

    Q_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :Q;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    vh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :vh ] )[1]

    θh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :θh ] )[1]

    #---------------------------------------------------

    target_parameter_sym =
        Symbol(
            "$(bus_name)_$(P_or_Q_or_Pll_or_Qll_sym)")

    change_in_parameter_dict =
        OrderedDict(
            :t => Float64[],
            target_parameter_sym =>
                    Float64[] )

    #---------------------------------------------------
    # Save parameters base value
    #---------------------------------------------------

    # Get the Power at bus 4

    var_normal_value =
        system_integrator.p.generic_model_dynamics_para[
            var_idx]
    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
       system_integrator.p.generic_model_dynamics_para]...))

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    #----------------------------------------
    # simulation steps
    #----------------------------------------    

    tstop1 = pertubation_time             # time_final/10
    tstop2 = restoration_time             # time_final/8
    tstop3 = restoration_time + Δt1       # time_final/6
    tstop4 = restoration_time + Δt1 + Δt2 # time_final/2
    tstop5 = time_final/1.0    

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 1
    # -------------------------------------

    add_tstop!(
        system_integrator, tstop1 )

    step!(system_integrator)

    push!(parameter_df,
        tuple(
            [[system_integrator.t];
       system_integrator.p.generic_model_dynamics_para]...))

    system_sol =
        system_integrator.sol

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    @show "step 1"
    
    # -------------------------------------
    # Steps 2
    # -------------------------------------

    # A pertubation of Active power at bus 4

    pertubation_factor =
        pertubation_factor

    pertubation_tstop =
        tstop2

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    @show "step 2"
    
    # -------------------------------------
    # Step 3
    # -------------------------------------

    # Bring back to normal the Power at bus 4

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop3

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    @show "step 3"
    
    # -------------------------------------
    # Step 4
    # -------------------------------------

    # Continute simulation 

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop4

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    system_sol = system_integrator.sol

    push!( parameter_df, tuple(
            [[system_integrator.t];
      system_integrator.p.generic_model_dynamics_para]...))

    @show "step 4"
    
    #---------------------------------------------------
    # simulate tstop4 till the end
    #---------------------------------------------------

    DifferentialEquations.solve!(system_integrator)

    system_sol = system_integrator.sol

    #---------------------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    final_intg_model_dynamics_para =
        getproperty(
            system_integrator.p,
            :generic_model_dynamics_para)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          final_intg_model_dynamics_para[var_idx] )

    @show "step 5"
    
    #---------------------------------------------------
    # extract solution auxilliary results
    #---------------------------------------------------

    sol_auxilliary_results =
        get_sol_auxilliary_results(
            system_sol;
            state_labels,
            algebraic_vars_labels,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs)

    auxilliary_results_plot =
        make_plot_gens_streamedlined_auxilliary_results(
            ;sol_auxilliary_results...)
    
    #---------------------------------------------------

    target_parameter_plot =
        plot(change_in_parameter_dict[:t],
             change_in_parameter_dict[target_parameter_sym],
             linetype=:steppre,
             yminorticks = 10,
             xminorticks = 10,
             fmt = :pdf,
             lw = 1,
             xlabel = "t [s]",
             ylabel = "$(target_parameter_sym) [p.u]",
             labels = "$(target_parameter_sym)",
             bottom_margin=2Plots.mm, # Adjust bottom margin
             left_margin=5Plots.mm,   # Adjust left margin
             right_margin=2Plots.mm,  # Adjust right margin
             top_margin=2Plots.mm )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan,

            sol_auxilliary_results,
            auxilliary_results_plot,
            target_parameter_plot )

end



function get_sim_mm_sudden_load_or_line_outage_pertubation(
    ;case_name,
    json_net_data_by_components_file =
        nothing,
    script_dir,
    data_dir  = "",
    components_lib = "",
    
    target_bus_name = "bus4",
    
    sim_type,

    outage_type =
        :line_outage,

    pertubation_factor,
    restoration_factor,
    pertubation_time,
    restoration_time,
    
    Δt1 = 1.5,
    Δt2 = 1.5,
    
    on_fault_time,
    clear_fault_time,
    
    line_outage_time,
    generation_adjustment_time,
    
    components_libs_dir,

    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit  = [1],

    list_edges_to_have_fault   = [ 8 ],
    clear_fault_selection_list = [1],
        
    list_network_status = 
            [:pre_fault_state,
             :post_fault_state],
    
    with_faults =
        false,    
    
    # # base setting and some booleans
    
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,

    # # solvers and settings
    
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    
    timespan      = 20,    
    dt            = 0.01,
    Δt            = 1.0 / 2^(4),

    # # digits rounding
    
    fractional_digits = 6 )
    
    #---------------------------------------------------

    # """

    # case_name =
    #     "case9"
    
    # json_net_data_by_components_file =
    #     json_net_data_by_components_file
    
    # script_dir =
    #     script_dir
    
    # target_bus_name =
    #     "bus5"
    
    # sim_type =
    #     sim_type

    # outage_type =
    #     :line_outage

    # # # Simulation timespan
    
    # # timespan    = 12.0
  
    # dt            = 0.01
    # Δt            = 1.0 / 2^(4)

    # Δt_clear_time = 0.5
  
    # Δt1 = 1.5
    # Δt2 = 1.5

    # on_fault_time =
    #     10.0
    
    # clear_fault_time =
    #     on_fault_time + Δt_clear_time
    
    # line_outage_time =
    #     on_fault_time
    
    # generation_adjustment_time =
    #     line_outage_time + Δt_generation_adjustment_time


    # pertubation_factor =
    #     pertubation_factor
    
    # restoration_factor =
    #     restoration_factor
    
    # pertubation_time =
    #     on_fault_time
    
    # restoration_time =
    #     clear_fault_time

    # timespan =
    #     clear_fault_time + Δt1 + Δt2 + 2.0 # 12.0

    # components_libs_dir =
    #     components_libs_dir

    # list_fault_point_from_node_a = [0.3]
    # list_fault_resistance = [0.001]
    # list_no_line_circuit  = [1]

    # list_edges_to_have_fault   = [ 8 ]
    # clear_fault_selection_list = [1]
        
    # list_network_status = 
    #         [:pre_fault_state,
    #          :post_fault_state]
    
    # with_faults =
    #     false   
    
    # # # base setting and some booleans
    
    # basekV = 1.0
    # use_pu_in_PQ = true
    # line_data_in_pu = true

    # use_init_u0 = false
    # use_nlsolve = false

    # # # solvers and settings
    
    # pf_alg        = NewtonRaphson()    
    # ode_alg       = Rodas4()
    # dae_alg       = IDA()
    
    # abstol        = 1e-12
    # reltol        = 1e-12

    # # # digits rounding
    
    # fractional_digits = 6

    # """
    
    #---------------------------------------------------

    cd(script_dir)

    results_dir =
        joinpath(
            script_dir,
            "results",
            sim_type)

    if !(isdir(results_dir))

        mkpath(results_dir)

    end


    figure_dir =
        joinpath(
            script_dir,
            "figure",
            sim_type)

    if !(isdir(figure_dir))

        mkpath(figure_dir)

    end


    #---------------------------------------------------

    # json_net_data_by_components_file =
    #     "net-static-data-avr-sauer-gov-sauer.json"
    # components_libs_dir =
    #     joinpath(script_dir,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(
    #         script_dir,
    #         "..","..","src",
    #         "data-dir",
    #         "converted_data" )


    # if json_net_data_by_components_file == nothing

    #     net_data_by_components_file =
    #         joinpath(
    #             json_case_dir,
    #             "net_data_by_components_file.json")
    # else

    #     net_data_by_components_file =
    #         joinpath(
    #             json_case_dir,
    #             json_net_data_by_components_file)

    # end


    
    if (components_lib == "") || (
        components_lib == nothing)

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
        joinpath( data_dir,
                  "converted-data",
                  case_name)


    json_case_dir =
        joinpath(case_data_dir,
            "json")

    if (json_net_data_by_components_file == "" ||
        json_net_data_by_components_file == nothing)
    
        
        net_data_by_components_file =
            joinpath(json_case_dir,
                     "net_data_by_components_file.json")
    else
    
        
        net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    sd_dynamics_sim_csv_filename =
        joinpath(
            results_dir,
            "$(case_name)-" *
                "$(sim_type)-" *
                "sim-sd-dynamics.csv")

    #---------------------------------------------------
    # Simulation Period
    #---------------------------------------------------

    time_start    = 0.0

    time_final    = timespan

    tspan         = (0.0, timespan)

    sim_timespan  = (0.0, timespan)

    plot_timespan = (0.0, timespan)

    #---------------------------------------------------

    @assert time_final >= restoration_time + Δt1 + Δt2

    @assert time_final >= clear_fault_time
    
    @assert time_final >= generation_adjustment_time
    
    
    #---------------------------------------------------
    ## variables, parameters, indices, intialisation
    #---------------------------------------------------

    ntuple_status_steady_state_data =
        get_ntuple_status_steady_state_data(
            ;with_faults =
                with_faults,
            net_data_by_components_file =
                net_data_by_components_file,
            components_libs_dir =
                components_libs_dir,

            timespan =
                timespan,
            on_fault_time =
                on_fault_time,
            clear_fault_time =
                clear_fault_time,

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

            basekV =
                basekV,    
            use_pu_in_PQ =
                use_pu_in_PQ,
            line_data_in_pu =
                line_data_in_pu,
            list_network_status =
                list_network_status )
   
    #---------------------------------------------------

    (;loc_load_exist,
     state_labels,
     algebraic_vars_labels,
     
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     system_fault_status,
     generic_system_dynamics_wt_fault_kwd_para,
     on_fault_net_para,
     cleared_selected_lines_faults_net_para,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

     # model_bool_dae_vars_wt_fault,
     # model_syms_wt_fault,         
     # u0_model_states_init_wt_fault,

     model_bool_dae_vars,     
     model_syms,
     u0_model_states_init,
     model_mass_matrix,     
     
     cb_states,
     plants_cb_paras_switches,

     nodes_names,
     gens_nodes_names,
     non_gens_nodes_names,
     SM_gens_nodes_names,
     SC_gens_nodes_names,

     ωref0_vref0_porder0_id_iq_vh_Idx,
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type,
     Gs,
     Bs,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    ntuple_status_steady_state_data,
                    :pre_fault_state),
                :static_prefault_paras),        
            (:loc_load_exist,
             :state_labels,
             :algebraic_vars_labels,
             :dyn_pf_fun_kwd_n2s_idxs,
             :dyn_pf_fun_kwd_net_idxs,
             
             :system_fault_status,
             :generic_system_dynamics_wt_fault_kwd_para,
             :on_fault_net_para,
             :cleared_selected_lines_faults_net_para,

             :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

             # :model_bool_dae_vars_wt_fault,
             # :model_syms_wt_fault,         
             # :u0_model_states_init_wt_fault,

             :model_bool_dae_vars,     
             :model_syms,
             :u0_model_states_init,
             :model_mass_matrix,             

             :cb_states,
             :plants_cb_paras_switches,

             :nodes_names,
             :gens_nodes_names,
             :non_gens_nodes_names,
             :SM_gens_nodes_names,
             :SC_gens_nodes_names,

             :ωref0_vref0_porder0_id_iq_vh_Idx,
             :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

             :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

             :edges_r,
             :edges_x,
             :edges_b,
             :edges_ratio,
             :edges_angle,
             :edges_type,
             :Gs,
             :Bs,
             :Ybr_cal_and_edges_orientation,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes))
    
    #---------------------------------------------------
    ## ODE system dyamanics simulation
    #---------------------------------------------------

    system_dynamics_fun! =
        mm_ode_generic_system_model_by_funcs_dynamics!
        # mm_ode_generic_system_dynamics_by_ode_pf_funcs!

    #----------------------------------------

    model_dynamics_kwd_para =
        generic_system_dynamics_wt_fault_kwd_para
        # generic_system_dynamics_kwd_para
    
    model_dynamics_para =
        (;generic_model_dynamics_para =
         ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
         plants_cb_paras_switches )          

    #----------------------------------------
    # integrator
    #----------------------------------------    

    system_integrator =
        DifferentialEquations.init(
            ODEProblem(
        ODEFunction(
        (dx,x,p,t) ->
            system_dynamics_fun!(
                dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms,
    mass_matrix = model_mass_matrix ) ,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para ),
            ode_alg,
            # dt = dt,
            
            callback = cb_states,
            tstops = [time_final],
            advance_to_tstop = true )

    #---------------------------------------------------
    # parameters df header
    #---------------------------------------------------

    generic_model_dynamics_para_df_header_sym =
        get_make_df_header_generic_model_dynamics_para(
            loc_load_exist,
            dyn_pf_fun_kwd_net_idxs)

    parameter_df =
        DataFrame(
            OrderedDict(a_header => Float64[]
                for a_header in
               generic_model_dynamics_para_df_header_sym))


    #---------------------------------------------------
    # parameters to be perturbed
    #---------------------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_loads_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx))
    
    bus_name =
        target_bus_name
    
    bus_no_or_bus_name =
        bus_name

    # bus_name = "bus4"

    P_or_Q_or_Pll_or_Qll_sym = :P


    var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            P_or_Q_or_Pll_or_Qll_sym;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    P_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :P;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    Q_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :Q;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    vh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :vh ] )[1]

    θh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :θh ] )[1]

    #---------------------------------------------------

    target_parameter_sym =
        Symbol(
            "$(bus_name)_$(P_or_Q_or_Pll_or_Qll_sym)")

    change_in_parameter_dict =
        OrderedDict(
            :t => Float64[],
            target_parameter_sym =>
                    Float64[] )

    #---------------------------------------------------
    # Save parameters base value
    #---------------------------------------------------

    # Get the Power at bus 4

    var_normal_value =
        system_integrator.p.generic_model_dynamics_para[
            var_idx]

    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
       system_integrator.p.generic_model_dynamics_para]...))

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    #----------------------------------------
    # simulation steps
    #----------------------------------------    
    
    tstop1 = pertubation_time             # time_final/10
    tstop2 = restoration_time             # time_final/8
    tstop3 = restoration_time + Δt1       # time_final/6
    tstop4 = restoration_time + Δt1 + Δt2 # time_final/2
    tstop5 = time_final/1.0    

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 1
    # -------------------------------------

    add_tstop!(system_integrator,
               tstop1 )

    step!(system_integrator)

    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
       system_integrator.p.generic_model_dynamics_para]...))

    system_sol =
        system_integrator.sol

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    @show "step 1"
    
    # -------------------------------------
    # Steps 2
    # -------------------------------------

    # A pertubation of Active power at bus 4

    pertubation_factor =
        pertubation_factor

    pertubation_tstop =
        tstop2

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    @show "step 2"
    
    # -------------------------------------
    # Step 3
    # -------------------------------------

    # Bring back to normal the Power at bus 4
    
    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop3

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    @show "step 3"
    
    # -------------------------------------
    # Step 4
    # -------------------------------------

    # Continute simulation 

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop4

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym  )

    system_sol = system_integrator.sol

    push!( parameter_df, tuple(
            [[system_integrator.t];
      system_integrator.p.generic_model_dynamics_para]...))

    @show "step 4"
    
    #---------------------------------------------------
    # simulate tstop4 till the end
    #---------------------------------------------------

    DifferentialEquations.solve!(system_integrator)

    system_sol = system_integrator.sol

    #---------------------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)


    final_intg_model_dynamics_para =
        getproperty(
            system_integrator.p,
            :generic_model_dynamics_para)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          final_intg_model_dynamics_para[var_idx] )

    @show "step 5"

    #---------------------------------------------------

    fractional_digits = 6

    sd_dynamics_sim_df = DataFrame(system_sol)

    sd_dynamics_sim_df[!, :] =
        round.(
            sd_dynamics_sim_df[:, :],
            digits=fractional_digits)

    sd_dynamics_sim_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                    "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
                     "$(bus_name)-" *
                     "$(sim_type)-states.csv")

    CSV.write(sd_dynamics_sim_csv_filename,
              sd_dynamics_sim_df )

    #---------------------------------------------------

    parameter_df[!, :] =
        round.(
            parameter_df[:, :],
            digits=fractional_digits)

    sd_dynamics_para_sim_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                    "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
                     "$(bus_no_or_bus_name)-" *
                     "$(sim_type)-parameters.csv")

    CSV.write(sd_dynamics_para_sim_csv_filename,
              parameter_df )

    #---------------------------------------------------

    generic_model_dynamics_para_df_header_idx =
        1:length(generic_model_dynamics_para_df_header_sym)

    paras_df_header_idx_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(bus_no_or_bus_name)-" *
                     "$(sim_type)-para-df-header-idx.csv")

    CSV.write(paras_df_header_idx_csv_filename,
              DataFrame(OrderedDict(
                  :idx =>
                generic_model_dynamics_para_df_header_idx,
                  :parameters =>
              generic_model_dynamics_para_df_header_sym)))

    #---------------------------------------------------

    generic_model_dynamics_state_df_header_sym =
        [[:t];model_syms]

    generic_model_dynamics_state_df_header_idx =
      1:length( generic_model_dynamics_state_df_header_sym )

    states_df_header_idx_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(bus_no_or_bus_name)-" *
                     "$(sim_type)-states-df-header-idx.csv")

    CSV.write(states_df_header_idx_csv_filename,
              DataFrame(OrderedDict(
                  :idx =>
                generic_model_dynamics_state_df_header_idx,
                  :states =>
           generic_model_dynamics_state_df_header_sym)))

    #---------------------------------------------------

    save_pertubation_stage_plot(
        case_name,
        system_sol;
        model_syms,
        gens_nodes_names,
        SM_gens_nodes_names,
        non_gens_nodes_names,
        sim_timespan,
        figure_dir,
        P_or_Q_or_Pll_or_Qll_sym,
        bus_idx = bus_no_or_bus_name)

    #---------------------------------------------------
    # extract solution auxilliary results
    #---------------------------------------------------

    sol_auxilliary_results =
        get_sol_auxilliary_results(
            system_sol;
            state_labels,
            algebraic_vars_labels,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    auxilliary_results_plot =
        make_plot_gens_streamedlined_auxilliary_results(;
            sol_auxilliary_results... )


    names_vars_plots =
        propertynames(auxilliary_results_plot)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                auxilliary_results_plot,
                a_vars_plots)

        local filename =
            "$(case_name)-" *
            "$(bus_no_or_bus_name)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "$(String(a_vars_plots)).pdf"

        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))

    end

    #---------------------------------------------------

    target_parameter_plot =
        plot(change_in_parameter_dict[:t],
             change_in_parameter_dict[target_parameter_sym],
             linetype=:steppre,
             yminorticks = 10,
             xminorticks = 10,
             fmt = :pdf,
             lw = 1,
             xlabel = "t [s]",
             ylabel = "$(target_parameter_sym) [p.u]",
             labels = "$(target_parameter_sym)",
             bottom_margin=2Plots.mm,  # Adjust bottom margin
             left_margin=5Plots.mm,   # Adjust left margin
             right_margin=2Plots.mm,  # Adjust right margin
             top_margin=2Plots.mm )

    target_parameter_plot_filename =
            "$(case_name)-" *
            "$(bus_no_or_bus_name)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "target-parameter-" *
            "$(String(target_parameter_sym)).pdf"

    savefig(target_parameter_plot,
                joinpath(
                    figure_dir,
                    target_parameter_plot_filename))    

    #---------------------------------------------------
    
    # guick_single_vars_plots_dae_or_ode_sol =
    #     get_guick_single_vars_plots_dae_or_ode_sol(
    #         ;system_sol =
    #         system_sol,
    #         model_syms =
    #             model_syms,
    #         gens_nodes_names =
    #             gens_nodes_names,
    #         SM_gens_nodes_names =
    #             SM_gens_nodes_names,
    #         non_gens_nodes_names =
    #             non_gens_nodes_names,
    #         sim_timespan = sim_timespan )

    #---------------------------------------------------
    # Save results to files
    #---------------------------------------------------

     return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan )

    # fractional_digits = 6

    # sd_dynamics_sim_df = DataFrame(system_sol)

    # sd_dynamics_sim_df[!, :] =
    #     round.(
    #         sd_dynamics_sim_df[:, :],
    #         digits=fractional_digits)

    # sd_dynamics_sim_csv_filename =
    #     joinpath(results_dir,
    #              "$(case_name)-" *
    #                 "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
    #                  "$(bus_name)-" *
    #                  "$(sim_type)-states.csv")

    # CSV.write(sd_dynamics_sim_csv_filename,
    #           sd_dynamics_sim_df )

    # #---------------------------------------------------

    # parameter_df[!, :] =
    #     round.(
    #         parameter_df[:, :],
    #         digits=fractional_digits)

    # sd_dynamics_para_sim_csv_filename =
    #     joinpath(results_dir,
    #              "$(case_name)-" *
    #                 "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
    #                  "$(bus_no_or_bus_name)-" *
    #                  "$(sim_type)-parameters.csv")

    # CSV.write(sd_dynamics_para_sim_csv_filename,
    #           parameter_df )

    # #---------------------------------------------------

    # generic_model_dynamics_para_df_header_idx =
    #     1:length(generic_model_dynamics_para_df_header_sym)

    # paras_df_header_idx_csv_filename =
    #     joinpath(results_dir,
    #              "$(case_name)-" *
    #                  "$(bus_no_or_bus_name)-" *
    #                  "$(sim_type)-para-df-header-idx.csv")

    # CSV.write(paras_df_header_idx_csv_filename,
    #           DataFrame(OrderedDict(
    #               :idx =>
    #             generic_model_dynamics_para_df_header_idx,
    #               :parameters =>
    #           generic_model_dynamics_para_df_header_sym)))

    # #---------------------------------------------------

    # generic_model_dynamics_state_df_header_sym =
    #     [[:t];model_syms]

    # generic_model_dynamics_state_df_header_idx =
    #   1:length( generic_model_dynamics_state_df_header_sym )

    # states_df_header_idx_csv_filename =
    #     joinpath(results_dir,
    #              "$(case_name)-" *
    #                  "$(bus_no_or_bus_name)-" *
    #                  "$(sim_type)-states-df-header-idx.csv")

    # CSV.write(states_df_header_idx_csv_filename,
    #           DataFrame(OrderedDict(
    #               :idx =>
    #             generic_model_dynamics_state_df_header_idx,
    #               :states =>
    #        generic_model_dynamics_state_df_header_sym)))


    # #---------------------------------------------------

    # save_pertubation_stage_plot(
    #     case_name,
    #     system_sol;
    #     model_syms,
    #     gens_nodes_names,
    #     SM_gens_nodes_names,
    #     non_gens_nodes_names,
    #     sim_timespan,
    #     figure_dir,
    #     P_or_Q_or_Pll_or_Qll_sym,
    #     bus_idx = bus_no_or_bus_name)

    # #---------------------------------------------------
    # # extract solution auxilliary results
    # #---------------------------------------------------

    # sol_auxilliary_results =
    #     get_sol_auxilliary_results(
    #         system_sol;
    #         state_labels,
    #         algebraic_vars_labels,

    #         dyn_pf_fun_kwd_n2s_idxs,
    #         dyn_pf_fun_kwd_net_idxs )

    # auxilliary_results_plot =
    #     make_plot_gens_streamedlined_auxilliary_results(;
    #         sol_auxilliary_results... )


    # names_vars_plots =
    #     propertynames(auxilliary_results_plot)

    # for a_vars_plots in names_vars_plots

    #     plots_fig =
    #         getproperty(
    #             auxilliary_results_plot,
    #             a_vars_plots)

    #     local filename =
    #         "$(case_name)-" *
    #         "$(bus_no_or_bus_name)-" *
    #         "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
    #         "pertubation-" *
    #         "$(String(a_vars_plots)).pdf"

    #     savefig(plots_fig,
    #             joinpath(
    #                 figure_dir,
    #                 filename))

    # end

    # #---------------------------------------------------

    # target_parameter_plot =
    #     plot(change_in_parameter_dict[:t],
    #          change_in_parameter_dict[target_parameter_sym],
    #          linetype=:steppre,
    #          yminorticks = 10,
    #          xminorticks = 10,
    #          fmt = :pdf,
    #          lw = 1,
    #          xlabel = "t [s]",
    #          ylabel = "$(target_parameter_sym) [p.u]",
    #          labels = "$(target_parameter_sym)",
    #          bottom_margin=2Plots.mm,  # Adjust bottom margin
    #          left_margin=5Plots.mm,   # Adjust left margin
    #          right_margin=2Plots.mm,  # Adjust right margin
    #          top_margin=2Plots.mm )

    # target_parameter_plot_filename =
    #         "$(case_name)-" *
    #         "$(bus_no_or_bus_name)-" *
    #         "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
    #         "pertubation-" *
    #         "target-parameter-" *
    #         "$(String(target_parameter_sym)).pdf"

    # savefig(target_parameter_plot,
    #             joinpath(
    #                 figure_dir,
    #                 target_parameter_plot_filename))


    # return get_guick_single_vars_plots_dae_or_ode_sol(
    #     ;system_sol =
    #     system_sol,
    #     model_syms =
    #         model_syms,
    #     gens_nodes_names =
    #         gens_nodes_names,
    #     SM_gens_nodes_names =
    #         SM_gens_nodes_names,
    #     non_gens_nodes_names =
    #         non_gens_nodes_names,
    #     sim_timespan = sim_timespan )

end

#---------------------------------------------------

 
function sim_sudden_load_change(
    ;case_name,
    script_dir,
    data_dir       = "",
    components_lib = "",
    
    json_net_data_by_components_file = "",
    sim_type,
    timespan,

    target_bus_name,

    pertubation_factor = 1.10,
    restoration_factor = 1.0,
    
    pertubation_time,
    restoration_time,
    Δt1,
    Δt2,

    fractional_digits,
    
    use_saveat      = false,
    ts              = 0.001,
    
    ode_alg         = Rodas4(),
    dae_alg         = IDA(),
    pf_alg          = NewtonRaphson(),

    dt              = 0.01,

    Δt              = 1.0 / 2^(4),

    basekV          = 1.0,

    use_pu_in_PQ    = true,

    line_data_in_pu = true,

    use_init_u0     = false,

    use_nlsolve     = false,

    abstol          = 1e-12,

    reltol          = 1e-12 )

    #---------------------------------------------------

    cd(script_dir)

    results_dir =
        joinpath(
            script_dir,
            "results",
            sim_type)

    if !(isdir(results_dir))

        mkpath(results_dir)

    end

    figure_dir =
        joinpath(
            script_dir,
            "figure",
            sim_type)

    if !(isdir(figure_dir))

        mkpath(figure_dir)

    end

    #---------------------------------------------------

    # json_net_data_by_components_file =
    #     "net-static-data-avr-sauer-gov-sauer.json"


    # components_libs_dir =
    #     joinpath(script_dir,"..","..","src",
    #              "components-lib")

    # data_dir =
    #     joinpath(
    #         script_dir,
    #         "..","..","src",
    #         "data-dir",
    #         "converted_data" )


    # json_case_dir =
    #     joinpath(
    #         data_dir,
    #         case_name,
    #         "json")

    # if json_net_data_by_components_file == nothing

    #     net_data_by_components_file =
    #         joinpath(
    #             json_case_dir,
    #             "net_data_by_components_file.json")
    # else

    #     net_data_by_components_file =
    #         joinpath(
    #             json_case_dir,
    #             json_net_data_by_components_file)

    # end

    #--------------------------------------
    
    if (components_lib == "") || (
        components_lib == nothing)

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
                 case_name,)


    json_case_dir =
        joinpath(
            case_data_dir,
            "json")

    if (json_net_data_by_components_file == "" ||
        json_net_data_by_components_file == nothing)
            
        net_data_by_components_file =
            joinpath(json_case_dir,
                     "net_data_by_components_file.json" )
    else
            
        net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------

    sd_dynamics_sim_csv_filename =
        joinpath(
            results_dir,
            "$(case_name)-" *
                "$(sim_type)-" *
                "sim-sd-dynamics.csv")

    #---------------------------------------------------
    # Simulation Period
    #---------------------------------------------------

    time_start    = 0.0

    time_final    = timespan

    # dt            = 0.01

    # Δt            = 1.0 / 2^(4)

    tspan         = (0.0, timespan)

    sim_timespan  = (0.0, timespan)

    plot_timespan = (0.0, timespan)


    #---------------------------------------------------
    # base setting and some booleans 
    #---------------------------------------------------

    # basekV = 1.0

    # use_pu_in_PQ = true

    # line_data_in_pu = true

    #---------------------------------------------------
    ## solvers and settings
    #---------------------------------------------------

    # use_init_u0 = false

    # use_nlsolve = false

    # pf_alg        = NewtonRaphson()

    # #--------------------   

    # # ode_alg       = Rodas4()

    # # dae_alg       = IDA()

    # abstol        = 1e-12

    # reltol        = 1e-12

    #---------------------------------------------------
    ## variables, parameters, indices, intialisation
    #---------------------------------------------------

    (;u0_model_states_init,
     model_syms,
     model_bool_dae_vars,

     model_mass_matrix,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
     plants_cb_paras_switches,

     generic_system_dynamics_kwd_para,

     gens_nodes_names,
     SM_gens_nodes_names,
     non_gens_nodes_names,

     cb_states,

     plants_states_syms,

     gens_nodes_idx,

     state_labels,
     algebraic_vars_labels,
     network_vars_labels,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_system_simulation_parameters(
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
                 reltol = reltol),

             (:u0_model_states_init,
              :model_syms,
              :model_bool_dae_vars,

              :model_mass_matrix,

              :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
              :plants_cb_paras_switches,

              :generic_system_dynamics_kwd_para,

              :gens_nodes_names,
              :SM_gens_nodes_names,
              :non_gens_nodes_names,

              :cb_states,

              :plants_states_syms,
              :gens_nodes_idx,

              :state_labels,
              :algebraic_vars_labels,
              :network_vars_labels,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #----------------------------------------


    (;loc_load_exist,
     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,

     state_labels,
     algebraic_vars_labels) =
         NamedTupleTools.select(
             generic_system_dynamics_kwd_para ,
             (:loc_load_exist,
              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,

              :state_labels,
              :algebraic_vars_labels))

    #---------------------------------------------------
    ## ODE system dyamanics simulation
    #---------------------------------------------------

    system_dynamics_fun! =
        mm_ode_generic_system_model_by_funcs_dynamics!
        # dae_generic_system_model_by_funcs_dynamics!

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll
        # deepcopy(ω_ref_v_ref_p_order_Png_Qng_Pll_Qll)

    model_dynamics_para =
        (;generic_model_dynamics_para,
          plants_cb_paras_switches )          

    model_dynamics_kwd_para =
        generic_system_dynamics_kwd_para

    #----------------------------------------
    # integrator
    #----------------------------------------    
    

    if  use_saveat == true

        system_integrator =
            DifferentialEquations.init(
                ODEProblem(
            ODEFunction(
            (dx,x,p,t) ->
                system_dynamics_fun!(
                    dx, x,
                    model_dynamics_para,
                    t;
                    kwd_para =
                        model_dynamics_kwd_para);
            syms =
                model_syms,
        mass_matrix = model_mass_matrix ) ,
            u0_model_states_init,
            sim_timespan,
            model_dynamics_para ),
                ode_alg,
                
                dt = dt,
                saveat = ts,

                callback = cb_states,
                tstops = [time_final],
                advance_to_tstop = true )
        
    else

        system_integrator =
            DifferentialEquations.init(
                ODEProblem(
            ODEFunction(
            (dx,x,p,t) ->
                system_dynamics_fun!(
                    dx, x,
                    model_dynamics_para,
                    t;
                    kwd_para =
                        model_dynamics_kwd_para);
            syms =
                model_syms,
        mass_matrix = model_mass_matrix ) ,
            u0_model_states_init,
            sim_timespan,
            model_dynamics_para ),
                ode_alg,
                dt = dt,

                callback = cb_states,
                tstops = [time_final],
                advance_to_tstop = true )
        
    end
    
    
    #---------------------------------------------------
    # parameters df header
    #---------------------------------------------------

    generic_model_dynamics_para_df_header_sym =
        get_make_df_header_generic_model_dynamics_para(
            loc_load_exist,
            dyn_pf_fun_kwd_net_idxs)

    parameter_df =
        DataFrame(
            OrderedDict(a_header => Float64[]
                for a_header in
           generic_model_dynamics_para_df_header_sym))

    #---------------------------------------------------
    # parameters to be perturbed
    #---------------------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_loads_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx))

    # bus_no_or_bus_name = 4

    bus_no_or_bus_name = target_bus_name

    bus_name = target_bus_name

    # bus_name = "bus4"

    P_or_Q_or_Pll_or_Qll_sym = :P


    var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            P_or_Q_or_Pll_or_Qll_sym;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    P_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :P;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    Q_var_idx =
        get_P_or_Q_idx_in_generic_model_dynamics_para(
            bus_no_or_bus_name,
            :Q;
            loc_load_exist,
            dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    vh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :vh ] )[1]

    θh_var_idx_in_state =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                model_syms,
            bus_name = bus_name,
            vars = [ :θh ] )[1]

    #---------------------------------------------------

    target_parameter_sym =
        Symbol(
            "$(bus_name)_$(P_or_Q_or_Pll_or_Qll_sym)")

    change_in_parameter_dict =
        OrderedDict(
            :t => Float64[],
            target_parameter_sym =>
                    Float64[] )

    #---------------------------------------------------
    # Save parameters base value
    #---------------------------------------------------

    # Get the Power at bus 4

    var_normal_value =
        system_integrator.p.generic_model_dynamics_para[
            var_idx]

    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
      system_integrator.p.generic_model_dynamics_para]...))

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    #----------------------------------------
    # simulation steps
    #----------------------------------------    
    
    tstop1 = pertubation_time             # time_final/10
    tstop2 = restoration_time             # time_final/8
    tstop3 = restoration_time + Δt1       # time_final/6
    tstop4 = restoration_time + Δt1 + Δt2 # time_final/2
    tstop5 = time_final/1.0    
    
    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 1
    # -------------------------------------

    pertubation_count = 1

    add_tstop!(system_integrator,
               tstop1 )

    step!(system_integrator)

    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
     system_integrator.p.generic_model_dynamics_para]...))

    system_sol =
        system_integrator.sol

    #----------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          var_normal_value)

    # -------------------------------------
    # Steps 2
    # -------------------------------------

    pertubation_count = 2

    pertubation_factor =
        pertubation_factor

    pertubation_tstop = tstop2

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 3
    # -------------------------------------

    # Bring back to normal the Power at bus 4

    pertubation_count = 3

    # pertubation_factor = 1.0

    pertubation_factor =
        restoration_factor

    pertubation_tstop = tstop3

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym )

    # -------------------------------------
    # Step 4
    # -------------------------------------

    # maintaint the Power at bus 4

    pertubation_count = 4

    # pertubation_factor = 1.0
    
    pertubation_factor =
        pertubation_factor

    pertubation_tstop = tstop4

    pertubation_by_itegrator(
        var_normal_value,
        pertubation_factor,
        pertubation_tstop,
        var_idx,
        system_integrator;
        parameter_df,
        change_in_parameter_dict,
        target_parameter_sym  )

    system_sol = system_integrator.sol


    push!( parameter_df, tuple(
            [[system_integrator.t];
     system_integrator.p.generic_model_dynamics_para]...
                 ))

    #---------------------------------------------------
    # simulate tstop4 till the end
    #---------------------------------------------------

    DifferentialEquations.solve!(system_integrator)

    system_sol = system_integrator.sol

    #---------------------------------------------------

    push!(change_in_parameter_dict[:t],
                  system_integrator.t)


    final_intg_model_dynamics_para =
        getproperty(
            system_integrator.p,
            :generic_model_dynamics_para)

    push!(change_in_parameter_dict[
        target_parameter_sym],
          final_intg_model_dynamics_para[var_idx] )

    #---------------------------------------------------
    # Save results to files
    #---------------------------------------------------

    fractional_digits = fractional_digits

    sd_dynamics_sim_df = DataFrame(system_sol)

    sd_dynamics_sim_df[!, :] =
        round.(
            sd_dynamics_sim_df[:, :],
            digits=fractional_digits)

    # sd_dynamics_sim_csv_filename =
    #         joinpath(results_dir,
    #                  "$(case_name)-" *
    #                      "$(sim_type)-" *
    #                      "sim-sd-dynamics.csv")


    sd_dynamics_sim_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-"*
                  "$(String(P_or_Q_or_Pll_or_Qll_sym))-"*
                  "$(bus_name)-"*
                  "$(sim_type)-states.csv")

    CSV.write(sd_dynamics_sim_csv_filename,
              sd_dynamics_sim_df )

    #---------------------------------------------------

    parameter_df[!, :] =
        round.(
            parameter_df[:, :],
            digits=fractional_digits)

    sd_dynamics_para_sim_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-"*
                 "$(String(P_or_Q_or_Pll_or_Qll_sym))-"*
                 "$(bus_no_or_bus_name)-"*
                 "$(sim_type)-parameters.csv")

    CSV.write(sd_dynamics_para_sim_csv_filename,
              parameter_df )

    #---------------------------------------------------

    generic_model_dynamics_para_df_header_idx =
        1:length(generic_model_dynamics_para_df_header_sym)

    paras_df_header_idx_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-"*
                 "$(bus_no_or_bus_name)-"*
                 "$(sim_type)-para-df-header-idx.csv")

    CSV.write(paras_df_header_idx_csv_filename,
              DataFrame(OrderedDict(
                  :idx =>
              generic_model_dynamics_para_df_header_idx,
                  :parameters =>
              generic_model_dynamics_para_df_header_sym)))

    #---------------------------------------------------

    generic_model_dynamics_state_df_header_sym =
        [[:t];model_syms]

    generic_model_dynamics_state_df_header_idx =
        1:length(
            generic_model_dynamics_state_df_header_sym)

    states_df_header_idx_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                 "$(bus_no_or_bus_name)-" *
                 "$(sim_type)-states-df-header-idx.csv")

    CSV.write(states_df_header_idx_csv_filename,
              DataFrame(OrderedDict(
                  :idx =>
             generic_model_dynamics_state_df_header_idx,
                  :states =>
             generic_model_dynamics_state_df_header_sym)))

    #---------------------------------------------------

    save_pertubation_stage_plot(
        case_name,
        system_sol;
        model_syms,
        gens_nodes_names,
        SM_gens_nodes_names,
        non_gens_nodes_names,
        sim_timespan,
        figure_dir,
        P_or_Q_or_Pll_or_Qll_sym,
        bus_idx = bus_no_or_bus_name)

    #---------------------------------------------------
    # extract solution auxilliary results
    #---------------------------------------------------

    sol_auxilliary_results =
        get_sol_auxilliary_results(
            system_sol;
            state_labels,
            algebraic_vars_labels,

            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs )

    auxilliary_results_plot =
        make_plot_gens_streamedlined_auxilliary_results(;
            sol_auxilliary_results... )


    names_vars_plots =
        propertynames(auxilliary_results_plot)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                auxilliary_results_plot,
                a_vars_plots)

        local filename =
            "$(case_name)-" *
            "$(bus_no_or_bus_name)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "$(String(a_vars_plots)).pdf"

        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))

    end

    #---------------------------------------------------

    system_sol_jld2 =
        joinpath("$(results_dir)",
                 "system-sol.jld2")
    
    JLD2.@save "$(system_sol_jld2)" system_sol 

    #---------------------------------------------------

    target_parameter_plot =
        plot(change_in_parameter_dict[:t],
             change_in_parameter_dict[
                 target_parameter_sym],
             linetype=:steppre,
             yminorticks = 10,
             xminorticks = 10,
             fmt = :pdf,
             lw = 1,
             xlabel = "t [s]",
             ylabel = "$(target_parameter_sym) [p.u]",
             labels = "$(target_parameter_sym)",
             bottom_margin=2Plots.mm,# Adjust bottom margin
             left_margin=5Plots.mm,  # Adjust left margin
             right_margin=2Plots.mm, # Adjust right margin
             top_margin=2Plots.mm )

    target_parameter_plot_filename =
            "$(case_name)-" *
            "$(bus_no_or_bus_name)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "target-parameter-" *
            "$(String(target_parameter_sym)).pdf"

    savefig(target_parameter_plot,
                joinpath(
                    figure_dir,
                    target_parameter_plot_filename))

    result_sol_plot =
        get_guick_single_vars_plots_dae_or_ode_sol(
            ;system_sol =
            system_sol,
            model_syms =
                model_syms,
            gens_nodes_names =
                gens_nodes_names,
            SM_gens_nodes_names =
                SM_gens_nodes_names,
            non_gens_nodes_names =
                non_gens_nodes_names,
            sim_timespan = sim_timespan )

    return (; result_sol_plot,
            target_parameter_plot)

end


#---------------------------------------------------
#---------------------------------------------------


function sim_line_outage_pertubation_by_dae(
    outage_type,
    # outage_type = :line_outage,    
    # :line_outage_wt_pref_adjs
    # :line_outage_wt_vpref_adjs

    on_fault_time,
    clear_fault_time,
    
    line_outage_time,
    generation_adjustment_time;
    
    net_data_by_components_file,
    timespan = 20,
    
    components_libs_dir =
        nothing,
        
    list_fault_point_from_node_a = [0.3],
    list_fault_resistance = [0.001],
    list_no_line_circuit  = [1],

    list_edges_to_have_fault   = [ 8 ],
    clear_fault_selection_list = [1],
        
    list_network_status = 
            [:pre_fault_state,
             :post_fault_state],
    
    with_faults =
        false,    
    
    basekV          = 1.0,    
    use_pu_in_PQ    = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,
    
    pf_alg  = NewtonRaphson(),
    ode_alg = Rodas4(),    
    dae_alg = IDA(),
    
    abstol  = 1e-12,
    reltol  = 1e-12,
    
    dt = 0.0001,
    Δt = 1.0 / 2^(4) )

    
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

    #---------------------------------------------------
    # Simulation Period
    #---------------------------------------------------

    time_final    = timespan

    tspan         = (0.0, timespan)

    sim_timespan  = (0.0, timespan)

    plot_timespan = (0.0, timespan)
    
    #---------------------------------------------------

    ntuple_status_steady_state_data =
        get_ntuple_status_steady_state_data(
            ;with_faults =
                with_faults,
            net_data_by_components_file =
                net_data_by_components_file,
            components_libs_dir =
                components_libs_dir,

            timespan =
                timespan,
            on_fault_time =
                on_fault_time,
            clear_fault_time =
                clear_fault_time,

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

            basekV =
                basekV,    
            use_pu_in_PQ =
                use_pu_in_PQ,
            line_data_in_pu =
                line_data_in_pu,
            list_network_status =
                list_network_status )
   
    #---------------------------------------------------

    (;state_labels,
     algebraic_vars_labels,
     
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     system_fault_status,
     generic_system_dynamics_wt_fault_kwd_para,
     on_fault_net_para,
     cleared_selected_lines_faults_net_para,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

     # model_bool_dae_vars_wt_fault,
     # model_syms_wt_fault,         
     # u0_model_states_init_wt_fault,

     model_bool_dae_vars,     
     model_syms,
     u0_model_states_init,
     model_mass_matrix,     
     
     cb_states,
     plants_cb_paras_switches,

     nodes_names,
     gens_nodes_names,
     non_gens_nodes_names,
     SM_gens_nodes_names,
     SC_gens_nodes_names,

     ωref0_vref0_porder0_id_iq_vh_Idx,
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type,
     Gs,
     Bs,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    ntuple_status_steady_state_data,
                    :pre_fault_state),
                :static_prefault_paras),        
            (:state_labels,
             :algebraic_vars_labels,
             :dyn_pf_fun_kwd_n2s_idxs,
             :dyn_pf_fun_kwd_net_idxs,
             
             :system_fault_status,
             :generic_system_dynamics_wt_fault_kwd_para,
             :on_fault_net_para,
             :cleared_selected_lines_faults_net_para,

             :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

             # :model_bool_dae_vars_wt_fault,
             # :model_syms_wt_fault,         
             # :u0_model_states_init_wt_fault,

             :model_bool_dae_vars,     
             :model_syms,
             :u0_model_states_init,
             :model_mass_matrix,             

             :cb_states,
             :plants_cb_paras_switches,

             :nodes_names,
             :gens_nodes_names,
             :non_gens_nodes_names,
             :SM_gens_nodes_names,
             :SC_gens_nodes_names,

             :ωref0_vref0_porder0_id_iq_vh_Idx,
             :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

             :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

             :edges_r,
             :edges_x,
             :edges_b,
             :edges_ratio,
             :edges_angle,
             :edges_type,
             :Gs,
             :Bs,
             :Ybr_cal_and_edges_orientation,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #----------------------------------------

    # po := post_outage
    
    (po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    ntuple_status_steady_state_data,
                    :post_fault_state),
                :dynamic_status_paras),
            (:ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,))

    #----------------------------------------
    #----------------------------------------

    (Ynet, ) =
         NamedTupleTools.select(
        Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:Ynet, ) )

    #----------------------------------------

    (fault_Ynet,
     post_fault_Ynet) =
        NamedTupleTools.select(
        cleared_selected_lines_faults_net_para,
            (:pre_clear_fault_Ynet,
             :post_clear_fault_Ynet, ))

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

    #---------------------------------------

    @show system_fault_status
    
    if system_fault_status[1] != 0

        system_fault_status[1] = 0
        
    end
    
    #---------------------------------------
    #---------------------------------------
    # line loss with only porder_adj
    #---------------------------------------
    #---------------------------------------
        
    cb_line_outage = DiscreteCallback(
        (u, t, integrator) ->
            on_line_outage_condition(
                u, t, integrator,
                line_outage_time),

       on_line_outage_affect!;
        save_positions=(true, true),
        initializealg =
            ShampineCollocationInit() )

    cb_outage_set =
        CallbackSet(cb_line_outage,)

    tstop_outage =
        [line_outage_time]

    
    if outage_type == :line_outage_wt_pref_adjs

        gens_porder_adj =
            po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_p_order_Idx]

        cb_gens_porder_adjustment = DiscreteCallback(
            (u, t, integrator) ->
                on_generation_adjustment_condition(
                    u, t, integrator,
                    generation_adjustment_time),

            (integrator) ->
                on_generation_adjustment_affect!(
                    integrator,
                    gens_porder_adj,
                    dyn_p_order_Idx );
            save_positions=(true, true),
            initializealg =
                ShampineCollocationInit() )


        cb_outage_set =
            CallbackSet(cb_line_outage,
            cb_gens_porder_adjustment)

        tstop_outage =
            [line_outage_time,
             generation_adjustment_time]

        
    elseif outage_type == :line_outage_wt_vpref_adjs

        vref_and_porder_Idx =
            [ dyn_v_ref_Idx; dyn_p_order_Idx]

        vref_and_porder_adj =
            po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                vref_and_porder_Idx]
        
        cb_vref_and_porder_adj = DiscreteCallback(
            (u, t, integrator) ->
                on_generation_adjustment_condition(
                    u, t, integrator,
                    generation_adjustment_time),

            (integrator) ->
                on_generation_adjustment_affect!(
                integrator,
                vref_and_porder_adj,
                vref_and_porder_Idx);
            save_positions=(true, true),
            initializealg =
                ShampineCollocationInit() )

        cb_outage_set =
            CallbackSet(cb_line_outage,
            cb_vref_and_porder_adj)

        tstop_outage =
            [line_outage_time,
             generation_adjustment_time]
        
    else # :line_outage

        nothing
    end

    #----------------------------------------
    #----------------------------------------

    du0_model_states_init =
        zeros( length( u0_model_states_init ))

    res = similar( u0_model_states_init )

    
    #---------------------------------------

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    model_dynamics_para =
        (;generic_model_dynamics_para,
         Ynet,
         post_fault_Ynet,
         system_fault_status,
         plants_cb_paras_switches )

    model_dynamics_kwd_para =
        generic_system_dynamics_wt_fault_kwd_para

    #----------------------------------------
    
    system_dynamics_fun! =
        line_outage_generic_dynamics_wt_pre_post_fault_by_dae_pf_funcs!

    # line_loss_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!

    #----------------------------------------
    
    system_sol =
        DifferentialEquations.solve(
            DAEProblem(
        DAEFunction(
            (res, dx, x, p, t) ->
                system_dynamics_fun!(
                    res, dx, x,
                    model_dynamics_para,
                    t;
                    kwd_para =
                        model_dynamics_kwd_para);
            syms =
                model_syms),
                du0_model_states_init,
                u0_model_states_init,
                sim_timespan,
                model_dynamics_para,
                differential_vars =
                    model_bool_dae_vars,
                callback =
                    cb_states),
            dae_alg,
            callback =
                cb_outage_set,
            tstops =
                tstop_outage,
            abstol = abstol,
            reltol = reltol )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan,

            state_labels,
            algebraic_vars_labels,
     
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs,

            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle,
            edges_type,
            Gs,
            Bs,
            Ybr_cal_and_edges_orientation,
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            ntuple_status_steady_state_data )    
end



function sim_line_outage_pertubation_by_mm_ode(
    outage_type,
    # outage_type = :line_outage,    
    # :line_outage_wt_pref_adjs
    # :line_outage_wt_vpref_adjs

    on_fault_time,
    clear_fault_time,
    
    line_outage_time,
    generation_adjustment_time;
    
    net_data_by_components_file,
    timespan = 20,
    
    components_libs_dir =
        nothing,
        
    list_fault_point_from_node_a = [0.3],
    list_fault_resistance        = [0.001],
    list_no_line_circuit         = [1],

    list_edges_to_have_fault     = [ 8 ],
    clear_fault_selection_list   = [1],
        
    list_network_status = 
            [:pre_fault_state,
             :post_fault_state],
    
    with_faults =
        false,    
    
    basekV          = 1.0,    
    use_pu_in_PQ    = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,
    
    pf_alg  = NewtonRaphson(),
    ode_alg = Rodas4(),    
    dae_alg = IDA(),
    
    abstol  = 1e-12,
    reltol  = 1e-12,
    
    dt = 0.0001,
    Δt = 1.0 / 2^(4) )

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

    
    #---------------------------------------------------
    # Simulation Period
    #---------------------------------------------------

    time_final    = timespan

    tspan         = (0.0, timespan)

    sim_timespan  = (0.0, timespan)

    plot_timespan = (0.0, timespan)
    
    #---------------------------------------------------

    ntuple_status_steady_state_data =
        get_ntuple_status_steady_state_data(
            ;with_faults =
                with_faults,
            
            net_data_by_components_file =
                net_data_by_components_file,
            
            components_libs_dir =
                components_libs_dir,

            timespan =
                timespan,
            
            on_fault_time =
                on_fault_time,
            
            clear_fault_time =
                clear_fault_time,

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

            basekV =
                basekV,    
            use_pu_in_PQ =
                use_pu_in_PQ,
            line_data_in_pu =
                line_data_in_pu,
            list_network_status =
                list_network_status )
   
    #---------------------------------------------------

    (;state_labels,
     algebraic_vars_labels,
     
     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     system_fault_status,
     generic_system_dynamics_wt_fault_kwd_para,
     on_fault_net_para,
     cleared_selected_lines_faults_net_para,

     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

     # model_bool_dae_vars_wt_fault,
     # model_syms_wt_fault,         
     # u0_model_states_init_wt_fault,

     model_bool_dae_vars,     
     model_syms,
     u0_model_states_init,
     model_mass_matrix,     
     
     cb_states,
     plants_cb_paras_switches,

     nodes_names,
     gens_nodes_names,
     non_gens_nodes_names,
     SM_gens_nodes_names,
     SC_gens_nodes_names,

     ωref0_vref0_porder0_id_iq_vh_Idx,
     dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type,
     Gs,
     Bs,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    ntuple_status_steady_state_data,
                    :pre_fault_state),
                :static_prefault_paras),        
            (:state_labels,
             :algebraic_vars_labels,
             :dyn_pf_fun_kwd_n2s_idxs,
             :dyn_pf_fun_kwd_net_idxs,
             
             :system_fault_status,
             :generic_system_dynamics_wt_fault_kwd_para,
             :on_fault_net_para,
             :cleared_selected_lines_faults_net_para,

             :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,

             # :model_bool_dae_vars_wt_fault,
             # :model_syms_wt_fault,         
             # :u0_model_states_init_wt_fault,

             :model_bool_dae_vars,     
             :model_syms,
             :u0_model_states_init,
             :model_mass_matrix,             

             :cb_states,
             :plants_cb_paras_switches,

             :nodes_names,
             :gens_nodes_names,
             :non_gens_nodes_names,
             :SM_gens_nodes_names,
             :SC_gens_nodes_names,

             :ωref0_vref0_porder0_id_iq_vh_Idx,
             :dyn_ωref0_vref0_porder0_id_iq_vh_Idx,

             :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

             :edges_r,
             :edges_x,
             :edges_b,
             :edges_ratio,
             :edges_angle,
             :edges_type,
             :Gs,
             :Bs,
             :Ybr_cal_and_edges_orientation,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes))

    #----------------------------------------

    # po := post_outage
    
    (po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    ntuple_status_steady_state_data,
                    :post_fault_state),
                :dynamic_status_paras),
            (:ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,))

    #----------------------------------------
    #----------------------------------------

    (Ynet, ) =
         NamedTupleTools.select(
        Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:Ynet, ) )

    #----------------------------------------

    (fault_Ynet,
     post_fault_Ynet) =
        NamedTupleTools.select(
        cleared_selected_lines_faults_net_para,
            (:pre_clear_fault_Ynet,
             :post_clear_fault_Ynet, ))

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

    #---------------------------------------

    @show system_fault_status
    
    if system_fault_status[1] != 0

        system_fault_status[1] = 0
        
    end
    
    #---------------------------------------
    #---------------------------------------
    # line loss with only porder_adj
    #---------------------------------------
    #---------------------------------------
        
    cb_line_outage = DiscreteCallback(
        (u, t, integrator) ->
            on_line_outage_condition(
                u, t, integrator,
                line_outage_time),

       on_line_outage_affect!;
        save_positions=(true, true),
        initializealg =
            ShampineCollocationInit() )

    cb_outage_set =
        CallbackSet(cb_line_outage,)

    tstop_outage =
        [line_outage_time]

    
    if outage_type == :line_outage_wt_pref_adjs

        gens_porder_adj =
            po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_p_order_Idx]

        cb_gens_porder_adjustment = DiscreteCallback(
            (u, t, integrator) ->
                on_generation_adjustment_condition(
                    u, t, integrator,
                    generation_adjustment_time),

            (integrator) ->
                on_generation_adjustment_affect!(
                    integrator,
                    gens_porder_adj,
                    dyn_p_order_Idx );
            save_positions=(true, true),
            initializealg =
                ShampineCollocationInit() )


        cb_outage_set =
            CallbackSet(cb_line_outage,
            cb_gens_porder_adjustment)

        tstop_outage =
            [line_outage_time,
             generation_adjustment_time]

        
    elseif outage_type == :line_outage_wt_vpref_adjs

        vref_and_porder_Idx =
            [ dyn_v_ref_Idx; dyn_p_order_Idx]

        vref_and_porder_adj =
            po_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                vref_and_porder_Idx]
        
        cb_vref_and_porder_adj = DiscreteCallback(
            (u, t, integrator) ->
                on_generation_adjustment_condition(
                    u, t, integrator,
                    generation_adjustment_time),

            (integrator) ->
                on_generation_adjustment_affect!(
                integrator,
                vref_and_porder_adj,
                vref_and_porder_Idx);
            save_positions=(true, true),
            initializealg =
                ShampineCollocationInit() )

        cb_outage_set =
            CallbackSet(cb_line_outage,
            cb_vref_and_porder_adj)

        tstop_outage =
            [line_outage_time,
             generation_adjustment_time]
        
    else # :line_outage

        nothing
    end

    #----------------------------------------
    #---------------------------------------

    generic_model_dynamics_para =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll

    model_dynamics_para =
        (;generic_model_dynamics_para,
         Ynet,
         post_fault_Ynet,
         system_fault_status,
         plants_cb_paras_switches )

    model_dynamics_kwd_para =
        generic_system_dynamics_wt_fault_kwd_para

    #----------------------------------------
    
    system_dynamics_fun! =
        mm_line_outage_generic_dynamics_wt_pre_post_fault_by_ode_pf_funcs!

    #----------------------------------------

    system_sol =
        DifferentialEquations.solve(
            ODEProblem(
        ODEFunction(
        (dx,x,p,t) ->
            system_dynamics_fun!(
                dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
            syms =
                model_syms,
            mass_matrix =
                model_mass_matrix ) ,
                u0_model_states_init,
                sim_timespan,
                model_dynamics_para,
                callback =
                    cb_states),
                ode_alg,
                callback =
                    cb_outage_set,
                tstops =
                    tstop_outage,
                abstol = abstol,
                reltol = reltol )
    
    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan,

            state_labels,
            algebraic_vars_labels,
     
            dyn_pf_fun_kwd_n2s_idxs,
            dyn_pf_fun_kwd_net_idxs,

            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle,
            edges_type,
            Gs,
            Bs,
            Ybr_cal_and_edges_orientation,
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            ntuple_status_steady_state_data )    
end

#---------------------------------------------------

# sim_line_loss_pertubation
function sim_line_outage_pertubation(
    net_data_by_components_file;

    with_faults,
    sim_type,

    case_name="",
    components_libs_dir="",
    data_dir="",

    # # fault or pertubation
    on_fault_time,
    clear_fault_time,
    line_outage_time,
    generation_adjustment_time,
    
    list_fault_point_from_node_a = [0.01],
    list_fault_resistance = [0.001],
    list_no_line_circuit  =  [1],
    list_edges_to_have_fault   = [ 4 ],
    clear_fault_selection_list = [1],

    # # solvers and settings
    pf_alg        = NewtonRaphson(),    
    ode_alg       = Rodas4(),
    dae_alg       = IDA(),
    abstol        = 1e-12,
    reltol        = 1e-12,

    # # Simulation timespan
    timespan      = 20,    
    dt            = 0.01,
    Δt            = 1.0 / 2^(4),

    # # base setting and some booleans
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    use_init_u0 = false,
    use_nlsolve = false,

    # # digits rounding
    fractional_digits = 6,

    list_network_status = 
                [:pre_fault_state,
                 :post_fault_state],

    list_outage_type =
        [:line_outage,
         :line_outage_wt_pref_adjs,
         :line_outage_wt_vpref_adjs],

    # results folders
    results_dir=nothing,
    figure_dir=nothing,
    sd_dynamics_sim_csv_filename=nothing )
    
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

    #--------------------------------------    

    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    #--------------------------------------
    
    case_data_dir =
        joinpath(data_dir,
                 "converted-data",
                 case_name )

    #--------------------------------------
    
    # Dict{Symbol, NTuple}()
    
    dict_outage_type_sol = Dict() 

    for an_outage_type in list_outage_type

        an_outage_type_sol =
            get_line_loss_outage_wt_or_no_ref_adjs(
                an_outage_type,
                
                on_fault_time,
                clear_fault_time,

                line_outage_time,
                generation_adjustment_time,
                
                net_data_by_components_file;

                timespan,
                
                with_faults,
                components_libs_dir,
                data_dir,    

                list_fault_point_from_node_a,
                list_fault_resistance,
                list_no_line_circuit,

                list_edges_to_have_fault,
                clear_fault_selection_list,

                basekV,    
                use_pu_in_PQ,
                line_data_in_pu,

                use_init_u0, use_nlsolve,

                list_network_status,

                pf_alg, ode_alg, dae_alg,
                abstol, reltol,
                dt, Δt)

        dict_outage_type_sol[an_outage_type] =
            an_outage_type_sol

        save_network_pertubation_sim_plot(
            case_name;
            figure_dir,
            sim_type =
                "$(sim_type)",
            line_in_fault_name =
                "line-8",
            include_v_θ_plot =
                true,
            sub_folder =
                "$(an_outage_type)",
            NamedTupleTools.select(
                an_outage_type_sol,
                (:system_sol,
                :model_syms,
                :gens_nodes_names,
                :SM_gens_nodes_names,
                :non_gens_nodes_names,
                :sim_timespan))...)

        save_sol_auxilliary_results_plot(
            case_name;
             NamedTupleTools.select(
                an_outage_type_sol,
                 (:system_sol,
                  :state_labels,
                  :algebraic_vars_labels,
                  :dyn_pf_fun_kwd_n2s_idxs,
                  :dyn_pf_fun_kwd_net_idxs))...,
            figure_dir,
            sim_type =
                "$(sim_type)",
            line_in_fault_name =
                "line-8",
            sub_folder =
                "$(an_outage_type)" )

        sd_dynamics_sim_df =
            DataFrame(an_outage_type_sol.system_sol)

        sd_dynamics_sim_df[!, :] =
            round.(
                sd_dynamics_sim_df[:, :],
                digits=fractional_digits)

        sd_dynamics_sim_csv_filename =
                joinpath(results_dir,
                         "$(case_name)-" *
                             "$(sim_type)-" *
                             "$(an_outage_type).csv")

        CSV.write(sd_dynamics_sim_csv_filename,
                  sd_dynamics_sim_df )

    end

    return (; dict_outage_type_sol)

end


#---------------------------------------------------

function sim_Ynet_pertubation(
    on_fault_time,
    clear_fault_time,    
    ntuple_status_steady_state_data)

    #----------------------------------------
    # Dynamics
    #----------------------------------------

    # """

    # Se(efd) = Ax * exp( Bx * efd)

    # S_e_max = Ax * exp(Bx * efd_max )

    # S_e_075 = Ax * exp(Bx * 3/4 * efd_max )

    # (S_e_max/S_e_075) = exp(Bx * efd_max * 1/4)

    # Bx = (4/efd_max) * log(S_e_max/S_e_075)

    # Ax = S_e_max/exp(Bx * efd_max )


    # ex1 = 3.3

    # ex2 = 4.5

    # S_ex1 = 0.6602

    # S_ex2 = 4.2662

    # tBx = (4/ex2) * log(S_ex2/S_ex1)

    # # 1.65860

    # tAx = S_ex2/exp(tBx * ex2 )

    # # 0.00244

    # """

    #----------------------------------------

    (;system_fault_status,
     generic_system_dynamics_wt_fault_kwd_para,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     on_fault_net_para,
     cleared_selected_lines_faults_net_para,

     model_bool_dae_vars_wt_fault,
     model_syms_wt_fault,         
     u0_model_states_init_wt_fault,

     cb_states,

     nodes_names,
     gens_nodes_names,
     non_gens_nodes_names,
     SM_gens_nodes_names,
     SC_gens_nodes_names,

     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,
     ωref0_vref0_porder0_id_iq_vh,
     ω_ref_v_ref_p_order_Png_Qng_Pll_Qll) =
        NamedTupleTools.select(
            getproperty(
                getproperty(
                    ntuple_status_steady_state_data,
                    :pre_fault_state),
                :static_prefault_paras),        
            (:system_fault_status,
             :generic_system_dynamics_wt_fault_kwd_para,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes,
             :on_fault_net_para,
             :cleared_selected_lines_faults_net_para,

             :model_bool_dae_vars_wt_fault,
             :model_syms_wt_fault,         
             :u0_model_states_init_wt_fault,

             :cb_states,

             :nodes_names,
             :gens_nodes_names,
             :non_gens_nodes_names,
             :SM_gens_nodes_names,
             :SC_gens_nodes_names,

             :δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,
             :ωref0_vref0_porder0_id_iq_vh,
             :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll))

    #----------------------------------------

    (Ynet, ) =
         NamedTupleTools.select(
        Ynet_wt_nodes_idx_wt_adjacent_nodes,
             (:Ynet, ) )

    #----------------------------------------

    (fault_Ynet,
     post_fault_Ynet) =
        NamedTupleTools.select(
        cleared_selected_lines_faults_net_para,
            (:pre_clear_fault_Ynet,
             :post_clear_fault_Ynet,))

    #----------------------------------------
    
    post_fault_paras =
        getproperty(
            getproperty(
                ntuple_status_steady_state_data,
                :post_fault_state),
            :dynamic_status_paras)
    
    #----------------------------------------

    # """
    
    # (faulty_Ynet,
    #  faulty_nodes_idx_with_adjacent_nodes_idx) =
    #     NamedTupleTools.select(
    #     on_fault_net_para,
    #         (:faulty_Ynet,
    #          :faulty_nodes_idx_with_adjacent_nodes_idx))
    
    # mismatch_faulty_fault_Ynet =
    #     map(x -> round.(x; digits=4),
    #         [t1_array - t2_array
    #          for (t1_array, t2_array) in
    #              zip(faulty_Ynet,
    #                  fault_Ynet)])
    
    # """

    #----------------------------------------
    # DAE system dyamanics simulation
    #----------------------------------------

    cb_on_fault = DiscreteCallback(
        (u, t, integrator) ->
            on_fault_condition(
                u, t, integrator,
                on_fault_time),

        on_fault_Ynet_affect!; 
        save_positions=(true,true),
        initializealg =
            ShampineCollocationInit() )

    cb_clear_fault = DiscreteCallback(
        (u, t, integrator) ->
            clear_fault_condition(
                u, t, integrator,
                clear_fault_time),

       clear_fault_Ynet_affect!;
        save_positions=(true,true),
        initializealg =
            ShampineCollocationInit())

    #--------------------------------------

    cb_faults =
        CallbackSet(
            cb_on_fault,
            cb_clear_fault )

    #----------------------------------------

    model_dynamics_para =
        (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
         Ynet,
         fault_Ynet,
         post_fault_Ynet,
         system_fault_status )          

    model_dynamics_kwd_para =
        generic_system_dynamics_wt_fault_kwd_para

    system_dynamics_fun! =
        Ynet_generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!

    #----------------------------------------

    model_bool_dae_vars =
        model_bool_dae_vars_wt_fault

    model_syms =
        model_syms_wt_fault

    u0_model_states_init =
        u0_model_states_init_wt_fault

    #----------------------------------------

    du0_model_states_init =
        zeros(length(u0_model_states_init))

    res = similar(u0_model_states_init)

    #----------------------------------------

    faults_and_clear_times =
        [on_fault_time,
         clear_fault_time]

    #----------------------------------------

    system_sol =
        DifferentialEquations.solve(
            DAEProblem(
        DAEFunction(
        (res, dx, x, p, t) ->
            system_dynamics_fun!(
                res, dx, x,
                model_dynamics_para,
                t;
                kwd_para =
                    model_dynamics_kwd_para);
        syms =
            model_syms),
        du0_model_states_init,
        u0_model_states_init,
        sim_timespan,
        model_dynamics_para,
        differential_vars =
            model_bool_dae_vars,
                callback = cb_states),
            dae_alg,
            callback = cb_faults,
            tstops = faults_and_clear_times,
            abstol = abstol,
            reltol = reltol )

    return (;system_sol,
            model_syms,
            gens_nodes_names,
            SM_gens_nodes_names,
            non_gens_nodes_names,
            sim_timespan)
    
end

#---------------------------------------------------
#---------------------------------------------------

