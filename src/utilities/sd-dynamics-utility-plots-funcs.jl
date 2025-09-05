# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za



#---------------------------------------------------
# Results utitlities
#---------------------------------------------------


function make_results_dir(
    ;base_dir=@__DIR__,
    res_dir="results_dir")
    
    res_dir = joinpath(base_dir,res_dir)

    if !(isdir(res_dir))
        mkpath(res_dir)
    end

    csv_dir = joinpath(res_dir,"csv")
    
        if !(isdir(csv_dir))
        mkpath(csv_dir)
        end
    
    fig_dir = joinpath(res_dir,"fig")
    
        if !(isdir(fig_dir))
        mkpath(fig_dir)
        end

    sol_dir = joinpath(res_dir,"sol")
    
        if !(isdir(sol_dir))
        mkpath(sol_dir)
        end
    
    return (; res_dir, csv_dir,
            fig_dir, sol_dir)
end



function get_sol_net_pertubation_results(
    system_sol;
    gens_nodes_names,
    plants_states_syms,
    state_labels)

    merged_plants_states_syms =
        [ a_sym for a_sym in
             Set([plants_states_syms...; ]) ]

    dict_state_var_idxs =
        Dict( a_sym =>
        [idx for idx in first.(
             get_nodes_state_algb_vars_indices_in_system(
                ; network_vars_labels =
                    state_labels,
                 nodes_name =
                     gens_nodes_names,
                 vars = [a_sym] ))]
          for a_sym in merged_plants_states_syms )
    

    
    dict_state_var_values =
        Dict( a_sym => system_sol[
            dict_state_var_idxs[a_sym], :]
              for a_sym in merged_plants_states_syms  )

     return NamedTupleTools.namedtuple(
        dict_state_var_values)
    
end


function get_sol_auxilliary_results(
    system_sol;
    state_labels,
    algebraic_vars_labels,
    
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs,

    Ybr_cal_and_edges_orientation =
        nothing,
    selected_gens_state_vars_syms =
        (:δ, :ω, :ed_dash, :eq_dash),
    streamed_lined = true )

    (gens_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :all_nodes_idx))


   (; 
    n2s_all_nodes_idx,
    ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (
             :n2s_all_nodes_idx,
             ))
    
    network_vars_labels =
        [state_labels;
         algebraic_vars_labels]
    
    (δ_idx_in_state,
     ω_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             get_gens_state_vars_idx_in_state(
                 network_vars_labels,
                 dyn_pf_fun_kwd_net_idxs,
                 n2s_all_nodes_idx;
                 selected_gens_state_vars_syms =
                     selected_gens_state_vars_syms),
             (:δ_idx_in_state,
              :ω_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state ))


    (; 
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state) =
         NamedTupleTools.select(
             get_state_vars_and_i_dq_Idx_in_state(
                 state_labels,
                 gens_nodes_idx,
                 all_nodes_idx ),
             ( 
              :vh_Idx_in_state,
              :θh_Idx_in_state,
              :id_Idx_in_state,
                 :iq_Idx_in_state))
        
    gens_δ  =
        system_sol[δ_idx_in_state,  :]

    gens_ω  =
        system_sol[ω_idx_in_state, :]

    gens_ed_dash =
        system_sol[ed_dash_idx_in_state, :]

    gens_eq_dash =
        system_sol[eq_dash_idx_in_state, :]

    gens_id =
        system_sol[id_Idx_in_state, :]

    gens_iq =
        system_sol[iq_Idx_in_state, :]

    gens_vh =
        (system_sol[vh_Idx_in_state, :])[
            gens_nodes_idx,:]

    gens_θh =
        (system_sol[θh_Idx_in_state, :])[
            gens_nodes_idx,:]


    non_gens_vh =
        (system_sol[vh_Idx_in_state, :])[
            non_gens_nodes_idx, :]

    non_gens_θh =
        (system_sol[θh_Idx_in_state, :])[
            non_gens_nodes_idx, :]

    gens_vd = vd = get_a_gen_vd.(
        gens_δ, gens_vh, gens_θh )

    gens_vq = vq = get_a_gen_vq.(
        gens_δ, gens_vh, gens_θh )

    gens_ph = ph = get_a_gen_ph.(
        vd, vq, gens_id, gens_iq)


    gens_qh = qh = get_a_gen_qh.(
        vd, vq, gens_id, gens_iq)    


    if Ybr_cal_and_edges_orientation != nothing

        (; edges_Ybr_cal,
          edges_orientation) =
            NamedTupleTools.select(
                Ybr_cal_and_edges_orientation,
                (:edges_Ybr_cal,
                 :edges_orientation))
        
        vh =
            system_sol[vh_Idx_in_state, :]

        θh =
            system_sol[θh_Idx_in_state, :]

    
       from_idxs = first.(edges_orientation)

       to_idxs = second.(edges_orientation)
        
       I_from_I_to =
            get_vec_I_from_I_to_by_vh_θh(
                vh, θh;
                edges_Ybr_cal,
                edges_orientation,
                n2s_all_nodes_idx )
         
        branches_currents =
            get_vec_branches_currents_by_vh_θh(
                vh, θh;
                edges_Ybr_cal,
                edges_orientation,
                n2s_all_nodes_idx )

        i_b = abs.(branches_currents)
        α   = angle.(branches_currents)
         
        branches_from_currents =
            get_vec_branches_from_currents_by_vh_θh(
                vh,
                θh;
                edges_Ybr_cal,
                edges_orientation,
                n2s_all_nodes_idx)

        f_i_b = abs.( branches_from_currents)
        f_α  = angle.( branches_from_currents)

         
        branches_to_currents =
            get_vec_branches_to_currents_by_vh_θh(
                vh,
                θh;
                edges_Ybr_cal,
                edges_orientation,
                n2s_all_nodes_idx)

        t_i_b = abs.( branches_to_currents)
        t_α   = angle.( branches_to_currents)
        
    end
    
    return (streamed_lined == false) && (
        Ybr_cal_and_edges_orientation != nothing) ? (
        ;t = system_sol.t,
        gens_nodes_idx,
        gens_δ, gens_ω,
        gens_ed_dash, gens_eq_dash,
        gens_id, gens_iq,
        gens_vh, gens_θh,
        non_gens_vh, non_gens_θh,
        gens_vd, gens_vq,
        gens_ph, gens_qh, vd, vq,
        ph, qh,
        vh, θh,
        from_idxs, to_idxs,
        i_b, α,
        I_from_I_to,
        f_i_b, f_α,
        t_i_b, t_α ) :
            Ybr_cal_and_edges_orientation !=
            nothing ? (
                ;t = system_sol.t,
                gens_nodes_idx,
                vd,
                vq,
                ph,
                qh,
                vh,
                θh,
                from_idxs,
                to_idxs,
                i_b,
                α,
                I_from_I_to,
                f_i_b,
                f_α,
                t_i_b,
                t_α ) : (
                ;t = system_sol.t,
                gens_nodes_idx,
                vd,
                vq,
                ph,
                qh )
    
end

# ---------------------------------------------------

function get_make_df_header_generic_model_dynamics_para(
    loc_load_exist,
    dyn_pf_fun_kwd_net_idxs)

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_with_loc_loads_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx))

    
    df_header = loc_load_exist == true ?
        ([[:t]...;
         [ Symbol("ω_ref_gen$(idx)")
                   for idx in gens_nodes_idx ]...;
         [ Symbol("v_ref_gen$(idx)")
           for idx in gens_nodes_idx  ]...;
         [ Symbol("p_order_gen$(idx)")
           for idx in gens_nodes_idx  ]...;

         [ Symbol("P_load_node$(idx)")
           for idx in non_gens_nodes_idx  ]...;
         [ Symbol("Q_load_node$(idx)")
           for idx in non_gens_nodes_idx  ]...;

         [ Symbol("Pg_load_node$(idx)")
           for idx in gens_with_loc_loads_idx ]...;
         [ Symbol("Qg_load_node$(idx)")
           for idx in gens_with_loc_loads_idx ]...;]) :
                       ([[:t]...;
                         [ Symbol("ω_ref_gen$(idx)")
                           for idx in gens_nodes_idx ]...;
                        [ Symbol("v_ref_gen$(idx)")
                          for idx in gens_nodes_idx  ]...;
                        [ Symbol("p_order_gen$(idx)")
                          for idx in gens_nodes_idx  ]...;

                        [ Symbol("P_load_node$(idx)")
                          for idx in non_gens_nodes_idx]...;
                        [ Symbol("Q_load_node$(idx)")
                          for idx in
                              non_gens_nodes_idx]...;])
    return df_header

end

# ------------------------------------------------------
#  Plot functions
# ------------------------------------------------------


function make_plot_helics_federate(
    ntuple_subscriptions;
    fedrate_name="",
    fmt = :pdf,
    lw = 1,
    bottom_margin=2Plots.mm,  # Adjust bottom margin
    left_margin=5Plots.mm,   # Adjust left margin
    right_margin=2Plots.mm,  # Adjust right margin
    top_margin=2Plots.mm     # Adjust top margin
    )

    subscription_propty =
        propertynames(ntuple_subscriptions)

    t = getproperty(ntuple_subscriptions, :t)
    
    dict_plot = Dict{Symbol, Plots.Plot{
        Plots.GRBackend}}()
    
    for a_propty  in subscription_propty

        if a_propty != :t

            legend_labels =
                "$(fedrate_name)_$(a_propty)"

            x_axis_label = "t [s]"
            
            y_axis_label = a_propty == :θ ? "$(a_propty) [rad]" : "$(a_propty) [p.u]"
            
            a_propty_value =
                getproperty(
                    ntuple_subscriptions,
                    a_propty)

            dict_plot[a_propty] =
                plot(t,
                     a_propty_value,
                     yminorticks = 10,
                     fmt = fmt,
                     lw = lw,
                     # ylims = ylims,
                     # tspan = tspan,
                     xlabel = x_axis_label,
                     ylabel = y_axis_label,
                     labels = legend_labels,
                     bottom_margin=bottom_margin, 
                     left_margin=left_margin,
                     right_margin=right_margin,
                     top_margin=top_margin )

        end
        
    end

    return NamedTupleTools.namedtuple(
        dict_plot)
    
end



function make_plot_gens_streamedlined_auxilliary_results(
    ;t,
    gens_nodes_idx,
    vd,
    vq,
    ph,
    qh,
    fmt = :pdf,
    lw = 1,
    bottom_margin=2Plots.mm,  # Adjust bottom margin
    left_margin=5Plots.mm,   # Adjust left margin
    right_margin=2Plots.mm,  # Adjust right margin
    top_margin=2Plots.mm     # Adjust top margin
    )

    streamedlined_gens_propty =
        [:vd,
         :vq,
         :ph,
         :qh]

    streamedlined_gens_value =
        [vd,
         vq,
         ph,
         qh]
        
    dict_plot = Dict{Symbol, Plots.Plot{
        Plots.GRBackend}}()
    
    for (a_gen_propty, a_gen_propty_value) in
        zip(streamedlined_gens_propty,
            streamedlined_gens_value)

        legend_labels =
            ["bus$(idx)_$(a_gen_propty)"
             for idx in
                 gens_nodes_idx]
    
        legend_labels =
            reshape(
                legend_labels, 1,
                length(legend_labels))
        x_axis_label = "t [s]"
        y_axis_label = "$(a_gen_propty) [p.u]"
        
        dict_plot[a_gen_propty] =
            plot(t,
                 a_gen_propty_value',
                yminorticks = 10,
                fmt = fmt,
                lw = lw,
                # ylims = ylims,
                # tspan = tspan,
                xlabel = x_axis_label,
                ylabel = y_axis_label,
                labels = legend_labels,
                bottom_margin=bottom_margin, 
                left_margin=left_margin,
                right_margin=right_margin,
                top_margin=top_margin )
        
    end

    return NamedTupleTools.namedtuple(
        dict_plot)
    
end



function make_plot_gens_streamedlined_auxilliary_results(
    gens_nodes_idx,
    strd_gens_aux_results;
    fmt = :pdf,
    lw = 1,
    bottom_margin=2Plots.mm,  # Adjust bottom margin
    left_margin=5Plots.mm,   # Adjust left margin
    right_margin=2Plots.mm,  # Adjust right margin
    top_margin=2Plots.mm     # Adjust top margin
    )

    # gens_propty_sym =
    #     propertynames(strd_gens_aux_results)
    
    # streamedlined_gens_propty =
    #     NamedTupleTools.select(
    #        gens_propty_sym,
    #         (:vd,
    #          :vq,
    #          :ph,
    #          :qh) )

    streamedlined_gens_propty =
        [:vd,
         :vq,
         :ph,
         :qh]
    
    dict_plot = Dict{Symbol, Plots.Plot{
        Plots.GRBackend}}()
    
    for a_gen_propty in streamedlined_gens_propty

        legend_labels =
            ["bus$(idx)_$(a_gen_propty)"
             for idx in
                 gens_nodes_idx]
    
        legend_labels =
            reshape(
                legend_labels, 1,
                length(legend_labels))
        x_axis_label = "t [s]"
        y_axis_label = "$(a_gen_propty) [p.u]"

        a_gen_propty_value =
            getproperty(
                strd_gens_aux_results,
                a_gen_propty)
        
        dict_plot[a_gen_propty] =
            plot(strd_gens_aux_results.t,
                 a_gen_propty_value',
                yminorticks = 10,
                fmt = fmt,
                lw = lw,
                # ylims = ylims,
                # tspan = tspan,
                xlabel = x_axis_label,
                ylabel = y_axis_label,
                labels = legend_labels,
                bottom_margin=bottom_margin, 
                left_margin=left_margin,
                right_margin=right_margin,
                top_margin=top_margin )
        
    end

    return NamedTupleTools.namedtuple(
        dict_plot)
    

end



function make_plot_branches_current_type(
    t,
    i_b,
    α;
    ib_type = :i_b,
    α_type = :α,
    fmt = :pdf,
    lw = 1,
    bottom_margin=2Plots.mm,  # Adjust bottom margin
    left_margin=5Plots.mm,   # Adjust left margin
    right_margin=2Plots.mm,  # Adjust right margin
    top_margin=2Plots.mm     # Adjust top margin
    )

    streamedlined_propty =
        [ib_type, α_type]
    
    streamedlined_value =
        [ i_b,
          α]
    
    nrow_i_b, _ = size(i_b)
    
    dict_plot = Dict{Symbol, Plots.Plot{
        Plots.GRBackend}}()
    
    for (a_propty, a_propty_value) in
        zip(streamedlined_propty,
            streamedlined_value)

        legend_labels =
            ["branch$(idx)_$(a_propty)"
             for idx in
                 1:nrow_i_b]

        legend_labels =
            reshape(
                legend_labels, 1,
                length(legend_labels))
        

        if last(split(string(:f_i_b),"_")) == "b"

            y_axis_label = " |ib| [p.u]"
            
        else
            
            y_axis_label = " α [rad]"
            
        end
        
        x_axis_label = "t [s]"
        
        dict_plot[a_propty] =
            plot(t,
                 a_propty_value',
                yminorticks = 10,
                fmt = fmt,
                lw = lw,
                # ylims = ylims,
                # tspan = tspan,
                xlabel = x_axis_label,
                ylabel = y_axis_label,
                labels = legend_labels,
                bottom_margin=bottom_margin, 
                left_margin=left_margin,
                right_margin=right_margin,
                top_margin=top_margin )
        
    end

    return NamedTupleTools.namedtuple(
        dict_plot)
    
end



function make_plot_streamedlined_auxilliary_results(
    ;t,
    gens_nodes_idx,
    vd,
    vq,
    ph,
    qh,
    i_b,
    α,    
    fmt = :pdf,
    lw = 1,
    bottom_margin=2Plots.mm,  # Adjust bottom margin
    left_margin=5Plots.mm,   # Adjust left margin
    right_margin=2Plots.mm,  # Adjust right margin
    top_margin=2Plots.mm     # Adjust top margin
    )

    streamedlined_propty =
        [:vd,
         :vq,
         :ph,
         :qh,
         :i_b,
         :α]

    streamedlined_value =
        [vd,
         vq,
         ph,
         qh,
         i_b,
         α]

    nrow_i_b, _ = size(i_b)
    
    dict_plot = Dict{Symbol, Plots.Plot{
        Plots.GRBackend}}()
    
    for (a_propty, a_propty_value) in
        zip(streamedlined_propty,
            streamedlined_value)

        if a_propty == :i_b

            legend_labels =
                ["branch$(idx)_$(a_propty)"
                 for idx in
                     1:nrow_i_b]

            legend_labels =
                reshape(
                    legend_labels, 1,
                    length(legend_labels))
            
        elseif a_propty == :α

            legend_labels =
                ["branch$(idx)_$(a_propty)"
                 for idx in
                     1:nrow_i_b]

            legend_labels =
                reshape(
                    legend_labels, 1,
                    length(legend_labels))            
        else
            
            legend_labels =
                ["bus$(idx)_$(a_propty)"
                 for idx in
                     gens_nodes_idx]

            legend_labels =
                reshape(
                    legend_labels, 1,
                    length(legend_labels))
        end
        
        x_axis_label = "t [s]"
        y_axis_label = "$(a_propty) [p.u]"
        
        dict_plot[a_propty] =
            plot(t,
                 a_propty_value',
                yminorticks = 10,
                fmt = fmt,
                lw = lw,
                # ylims = ylims,
                # tspan = tspan,
                xlabel = x_axis_label,
                ylabel = y_axis_label,
                labels = legend_labels,
                bottom_margin=bottom_margin, 
                left_margin=left_margin,
                right_margin=right_margin,
                top_margin=top_margin )
        
    end

    return NamedTupleTools.namedtuple(
        dict_plot)
    
end


#---------------------------------------------------
#---------------------------------------------------

function save_pertubation_stage_plot(
    case_name,
    system_sol;
    model_syms,
    gens_nodes_names,
    SM_gens_nodes_names,
    non_gens_nodes_names,
    sim_timespan,
    figure_dir,
    P_or_Q_or_Pll_or_Qll_sym,
    
    bus_idx="",

    sub_folder =
        nothing)

   if sub_folder == nothing
        nothing
    else
        figure_dir =
            joinpath(
                figure_dir,
                sub_folder)

         if !(isdir(figure_dir))

             mkpath(figure_dir)

         end
        
    end
    
    group_vars_plots_dae_or_ode =
        get_guick_group_vars_plots_dae_or_ode_sol(
            ; system_sol =
                system_sol,
            model_syms =
                model_syms,
            gens_nodes_names =
                gens_nodes_names,
            SM_gens_nodes_names =
                SM_gens_nodes_names,
            non_gens_nodes_names =
                non_gens_nodes_names,
            sim_timespan =
                sim_timespan )

    names_vars_plots =
        propertynames(group_vars_plots_dae_or_ode)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                group_vars_plots_dae_or_ode,
                a_vars_plots)
        
        filename =
            "$(case_name)-" *
            "$(bus_idx)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "$(String(a_vars_plots)).pdf"
        
        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))
        
    end


    guick_single_vars_plots_dae_or_ode =
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
            sim_timespan =
                sim_timespan )


    single_names_vars_plots =
        propertynames(
            guick_single_vars_plots_dae_or_ode)
    

    for a_vars_plots in single_names_vars_plots

        plots_fig =
            getproperty(
                guick_single_vars_plots_dae_or_ode,
                a_vars_plots)
        
        filename =
            "$(case_name)-" *
            "$(bus_idx)-" *
            "$(String(P_or_Q_or_Pll_or_Qll_sym))-" *
            "pertubation-" *
            "$(String(a_vars_plots)).pdf"
        
        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))        
    end
    
end



function save_co_sim_stage_plot(
    case_name,
    system_sol;
    model_syms,
    gens_nodes_names,
    SM_gens_nodes_names,
    non_gens_nodes_names,
    sim_timespan,
    co_figure_dir,
    fedrate_name,
    sim_type,
    bus_idx="")

    group_vars_plots_dae_or_ode =
        get_guick_group_vars_plots_dae_or_ode_sol(
            ; system_sol =
                system_sol,
            model_syms =
                model_syms,
            gens_nodes_names =
                gens_nodes_names,
            SM_gens_nodes_names =
                SM_gens_nodes_names,
            non_gens_nodes_names =
                non_gens_nodes_names,
            sim_timespan =
                sim_timespan )

    names_vars_plots =
        propertynames(group_vars_plots_dae_or_ode)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                group_vars_plots_dae_or_ode,
                a_vars_plots)
        
        filename =
            "$(case_name)-" *
            "$(fedrate_name)-" *
            "$(sim_type)-" *
            "$(String(a_vars_plots)).pdf"
        
        savefig(plots_fig,
                joinpath(
                    co_figure_dir,
                    filename))
        
    end

    guick_single_vars_plots_dae_or_ode =
        get_guick_single_vars_plots_dae_or_ode_sol(
            ; system_sol =
                system_sol,
            model_syms =
                model_syms,
            gens_nodes_names =
                gens_nodes_names,
            SM_gens_nodes_names =
                SM_gens_nodes_names,
            non_gens_nodes_names =
                non_gens_nodes_names,
            sim_timespan =
                sim_timespan )


    single_names_vars_plots =
        propertynames(
            guick_single_vars_plots_dae_or_ode)
    

    for a_vars_plots in single_names_vars_plots

        plots_fig =
            getproperty(
                guick_single_vars_plots_dae_or_ode,
                a_vars_plots)
        
        filename =
            "$(case_name)-" *
            "$(fedrate_name)-" *
            "$(sim_type)-" *
            "$(String(a_vars_plots)).pdf"
        
        savefig(plots_fig,
                joinpath(
                    co_figure_dir,
                    filename))
        
    end
    
end



function save_co_sim_helics_federate_plot(
    ntuple_subscriptions;
    case_name="",
    co_figure_dir,
    fedrate_name,
    sim_type,
    sub_folder =
        nothing,
    comment = nothing )

   if sub_folder == nothing
        nothing
    else
        co_figure_dir =
            joinpath(
                co_figure_dir,
                sub_folder)

         if !(isdir(co_figure_dir))

             mkpath(co_figure_dir)

         end
        
    end

    
    helics_federate_plots =
        make_plot_helics_federate(
            ntuple_subscriptions;
            fedrate_name=
                fedrate_name)

    
    helics_federate_plots_syms =
        propertynames(
            helics_federate_plots)
    

    for a_vars_plots in helics_federate_plots_syms

        plots_fig =
            getproperty(
                helics_federate_plots,
                a_vars_plots )

        if comment != nothing
        
            filename =
                "$(case_name)-" *
                "$(fedrate_name)-" *
                "$(sim_type)-" *
                "$(comment)-" *
                "$(String(a_vars_plots)).pdf"
        else
            

            filename =
                "$(case_name)-" *
                "$(fedrate_name)-" *
                "$(sim_type)-" *
                "$(String(a_vars_plots)).pdf"
            
        end
        
        savefig(plots_fig,
                joinpath(
                    co_figure_dir,
                    filename))
        
    end
    
end



function save_sol_auxilliary_results_plot(
    case_name,
    ntuple_auxilliary_plot;
    figure_dir,
    sim_type,
    line_in_fault_name,
    sub_folder =
        nothing,
    bus_idx="",
    fedrate_name=
        nothing )

   if sub_folder == nothing
        nothing
    else
        figure_dir =
            joinpath(
                figure_dir,
                sub_folder)

         if !(isdir(figure_dir))

             mkpath(figure_dir)

         end
        
    end

    names_vars_plots =
        propertynames(ntuple_auxilliary_plot)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                ntuple_auxilliary_plot,
                a_vars_plots)
        
        filename =
            "$(case_name)-" *
            "$(sim_type)-" *
            "$(line_in_fault_name)-" *
            "$(String(a_vars_plots)).pdf"

        if fedrate_name != nothing
        
            filename =
                "$(case_name)-" *
                "$(fedrate_name)-" *
                "$(sim_type)-" *
                "$(line_in_fault_name)-" *
                "$(String(a_vars_plots)).pdf"
        end
        
        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))

    end

    
end



function save_sol_auxilliary_results_plot(
    case_name;
    system_sol,    
    state_labels,
    algebraic_vars_labels,
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs,
    
    figure_dir,
    sim_type,
    
    line_in_fault_name,
    bus_idx="",
    
    sub_folder =
        nothing,
    
    fedrate_name=
        nothing,

    selected_gens_state_vars_syms =
        (:δ, :ω, :ed_dash, :eq_dash),
    streamed_lined = true )


    ntuple_auxilliary_plot =
        make_plot_gens_streamedlined_auxilliary_results(
            ;get_sol_auxilliary_results(
                system_sol;
                state_labels,
                algebraic_vars_labels,

                dyn_pf_fun_kwd_n2s_idxs,
                dyn_pf_fun_kwd_net_idxs,

                selected_gens_state_vars_syms =
                    selected_gens_state_vars_syms,
                streamed_lined =
                    streamed_lined )...)


    if sub_folder == nothing
        nothing
    else
        figure_dir =
            joinpath(
                figure_dir,
                sub_folder)

         if !(isdir(figure_dir))

             mkpath(figure_dir)

         end
        
    end

    names_vars_plots =
        propertynames(ntuple_auxilliary_plot)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                ntuple_auxilliary_plot,
                a_vars_plots)
        
        filename =
            "$(case_name)-" *
            "$(sim_type)-" *
            "$(line_in_fault_name)-" *
            "$(String(a_vars_plots)).pdf"

        if fedrate_name != nothing
        
            filename =
                "$(case_name)-" *
                "$(fedrate_name)-" *
                "$(sim_type)-" *
                "$(line_in_fault_name)-" *
                "$(String(a_vars_plots)).pdf"
        end

        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))

    end

    
end



function save_generic_sol_auxilliary_results_plot(
    case_name;
    system_sol,    
    state_labels,
    algebraic_vars_labels,
    dyn_pf_fun_kwd_n2s_idxs,
    dyn_pf_fun_kwd_net_idxs,
    
    figure_dir,
    sim_type,
    line_in_fault_name,
    sub_folder =
        nothing,

    selected_gens_state_vars_syms =
        (:δ, :ω, :ed_dash, :eq_dash),
    streamed_lined = true,
    bus_idx="")


    ntuple_auxilliary_plot =
        make_plot_gens_streamedlined_auxilliary_results(
            ;get_sol_auxilliary_results(
                system_sol;
                state_labels,
                algebraic_vars_labels,

                dyn_pf_fun_kwd_n2s_idxs,
                dyn_pf_fun_kwd_net_idxs,

                selected_gens_state_vars_syms =
                    selected_gens_state_vars_syms,
                streamed_lined =
                    streamed_lined )...)


    if sub_folder == nothing
        nothing
    else
        figure_dir =
            joinpath(
                figure_dir,
                sub_folder)

         if !(isdir(figure_dir))

             mkpath(figure_dir)

         end
        
    end

    names_vars_plots =
        propertynames(ntuple_auxilliary_plot)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                ntuple_auxilliary_plot,
                a_vars_plots)
        
        filename =
            "$(case_name)-" *
            "$(sim_type)-" *
            "$(line_in_fault_name)-" *
            "$(String(a_vars_plots)).pdf"

        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))

    end

    
end



function save_network_pertubation_sim_plot(
    case_name;
    system_sol,
    model_syms,
    gens_nodes_names,
    SM_gens_nodes_names,
    non_gens_nodes_names,
    sim_timespan,
    figure_dir,
    sim_type,
    line_in_fault_name,
    include_v_θ_plot =
        false,
    sub_folder =
        nothing,
    bus_idx="")

    group_vars_plots_dae_or_ode =
        get_guick_group_vars_plots_dae_or_ode_sol(
            ; system_sol =
                system_sol,
            model_syms =
                model_syms,
            gens_nodes_names =
                gens_nodes_names,
            SM_gens_nodes_names =
                SM_gens_nodes_names,
            non_gens_nodes_names =
                non_gens_nodes_names,
            sim_timespan =
                sim_timespan,
            include_v_θ_plot =
                include_v_θ_plot)

    if sub_folder == nothing
        nothing
    else
        figure_dir =
            joinpath(
                figure_dir,
                sub_folder)

         if !(isdir(figure_dir))

             mkpath(figure_dir)

         end
        
    end
    
    names_vars_plots =
        propertynames(group_vars_plots_dae_or_ode)

    for a_vars_plots in names_vars_plots

        plots_fig =
            getproperty(
                group_vars_plots_dae_or_ode,
                a_vars_plots)
        
        filename =
            "$(case_name)-" *
            "$(sim_type)-" *
            "$(line_in_fault_name)-" *
            "$(String(a_vars_plots)).pdf"
        
        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))
        
    end

    guick_single_vars_plots_dae_or_ode =
        get_guick_single_vars_plots_dae_or_ode_sol(
            ; system_sol =
                system_sol,
            model_syms =
                model_syms,
            gens_nodes_names =
                gens_nodes_names,
            SM_gens_nodes_names =
                SM_gens_nodes_names,
            non_gens_nodes_names =
                non_gens_nodes_names,
            sim_timespan =
                sim_timespan,
            include_v_θ_plot =
                include_v_θ_plot )


    single_names_vars_plots =
        propertynames(
            guick_single_vars_plots_dae_or_ode)
    

    for a_vars_plots in single_names_vars_plots

        plots_fig =
            getproperty(
                guick_single_vars_plots_dae_or_ode,
                a_vars_plots)
        
        filename =
            "$(case_name)-" *
            "$(sim_type)-" *
            "$(line_in_fault_name)-" *
            "$(String(a_vars_plots)).pdf"
        
        savefig(plots_fig,
                joinpath(
                    figure_dir,
                    filename))
        
    end
    
end


#---------------------------------------------------
#---------------------------------------------------


function get_guick_group_vars_plots_dae_or_ode_sol(
    ; system_sol = system_sol,
    model_syms = model_syms,
    gens_nodes_names = gens_nodes_names,
    SM_gens_nodes_names = SM_gens_nodes_names,
    non_gens_nodes_names = non_gens_nodes_names,
    sim_timespan = (0, 1.0),
    fmt = :pdf,
    lw = 1,
    ylims = (0.0, 1.2),
    include_v_θ_plot =
        true )


    list_plot_δ_ω_ed_dash_eq_dash =
        [ make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:δ ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "δ [rad]"),
         make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:ω ],
                tspan = sim_timespan,
             fmt = fmt,
             lw = lw,
         ylabel = "ω [p.u]"),
         make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:ed_dash ],
                tspan = sim_timespan ,
             fmt = fmt,
             lw = lw,
         ylabel = "ed' [p.u]"),
         make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:eq_dash ],
                tspan = sim_timespan,
             fmt = fmt,
             lw = lw,
         ylabel = "eq' [p.u]") ]

    plot_gens_states = plot(
            list_plot_δ_ω_ed_dash_eq_dash...;
            layout = (2,2),
            size = (1000, 500),
            lw = lw ,
        xlabel = "t[s]")


    list_plot_vr1_vr2_vf_tilade =
        [make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vr1 ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "vr1 [p.u]"),
         make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vr2 ],
                tspan = sim_timespan,
             fmt = fmt,
             lw = lw,
         ylabel = "vr2 [p.u]"),
         make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vf_tilade ],
                tspan = sim_timespan,
             fmt = fmt,
             lw = lw,
         ylabel = "vf [p.u]") ]

    plot_avrs_states = plot(
            list_plot_vr1_vr2_vf_tilade...;
            layout = (1,3),
            size = (1000, 500),
            lw = lw,
        xlabel = "t[s]")


    list_plot_xg1_xg2 =
        [make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = SM_gens_nodes_names,
                vars = [:xg1 ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "xg1 [p.u]"),
         make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = SM_gens_nodes_names,
                vars = [:xg2 ],
                tspan = sim_timespan,
             fmt = fmt,
             lw = lw,
         ylabel = "xg2 [p.u]") ]

    
    plot_govs_states = plot(
            list_plot_xg1_xg2...;
            layout = (1,2),
            size = (1000, 500),
            lw = lw,
        xlabel = "t[s]")

    if include_v_θ_plot == true

        list_plot_vh_θh =
            [make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vh ],
                tspan = sim_timespan,
                fmt = fmt,
                lw = lw,
            ylabel = "V [p.u]"),
             make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = non_gens_nodes_names,
                vars = [:vh ],
                tspan = sim_timespan,
                 fmt = fmt,
                 lw = lw,
             ylabel = "V [p.u]"),                   
             make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:θh ],
                tspan = sim_timespan,
                 fmt = fmt,
                 lw = lw,
             ylabel = "θ [rad]"),
             make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = non_gens_nodes_names,
                vars = [:θh ],
                tspan = sim_timespan,
                 fmt = fmt,
                 lw = lw,
             ylabel = "θ [rad]") ]

        plot_net_v_θ = plot(
                list_plot_vh_θh...;
                layout = (2,2),
                size = (1000, 500),
                lw = lw,
            xlabel = "t[s]")

        #
        
        list_plot_vh_ylims =
            [make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vh ],
                tspan = sim_timespan,
                fmt = fmt,
                lw = lw,
                ylims = ylims,
            ylabel = "V [p.u]"),
             make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = non_gens_nodes_names,
                vars = [:vh ],
                tspan = sim_timespan,
                 fmt = fmt,
                 lw = lw,
                 ylims = ylims,
             ylabel = "V [p.u]") ]

        plot_net_vh_ylims = plot(
                list_plot_vh_ylims...;
                layout = (1,2),
                size = (1000, 500),
                lw = lw,
            xlabel = "t[s]")
        
        
    end

        

    return include_v_θ_plot == true ?
        (; plot_gens_states,  plot_avrs_states,
         plot_govs_states, plot_net_v_θ, plot_net_vh_ylims
         ) : (; plot_gens_states,  plot_avrs_states,
              plot_govs_states)

end



function get_guick_single_vars_plots_dae_or_ode_sol(
    ; system_sol =
        system_sol,
    model_syms =
        model_syms,
    gens_nodes_names =
        gens_nodes_names,
    SM_gens_nodes_names =
        SM_gens_nodes_names,
    non_gens_nodes_names =
        non_gens_nodes_names,
    sim_timespan = (0 , 1.0),
    fmt = :pdf,
    lw = 1,
    ylims = (0.0, 1.2),
    include_v_θ_plot =
        true )


    δ_a_plot =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:δ ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "δ [rad]")

    ω_a_plot =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:ω ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "ω [p.u]")


    ed_dash_a_plot =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:ed_dash ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "ed' [p.u]")


    eq_dash_a_plot =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:eq_dash ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "eq' [p.u]")


    vr1_a_plot =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vr1 ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "vr1 [p.u]")


    vr2_a_plot =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vr2 ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "vr2 [p.u]")

    vf_tilade_a_plot =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vf_tilade ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "vf [p.u]")

    xg1_a_plot =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = SM_gens_nodes_names,
                vars = [:xg1 ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "xg1 [p.u]")


    xg2_a_plot =
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = SM_gens_nodes_names,
                vars = [:xg2 ],
                tspan = sim_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "xg2 [p.u]")


    if include_v_θ_plot == true

        plot_gens_vh =
            make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vh ],
                tspan = sim_timespan,
                fmt = fmt,
                lw = lw,
            ylabel = "V [p.u]")

        plot_gens_θh =
            make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:θh ],
                tspan = sim_timespan,
                fmt = fmt,
                lw = lw,
            ylabel = "θ [rad]")

        plot_non_gens_vh =
            make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = non_gens_nodes_names,
                vars = [:vh ],
                tspan = sim_timespan,
                fmt = fmt,
                lw = lw,
            ylabel = "V [p.u]")

        plot_non_gens_θh =
            make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = non_gens_nodes_names,
                vars = [:θh ],
                tspan = sim_timespan,
                fmt = fmt,
                lw = lw,
                ylabel = "θ [rad]")

        
        plot_gens_vh_ylims =
             make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = gens_nodes_names,
                vars = [:vh ],
                tspan = sim_timespan,
                fmt = fmt,
                lw = lw,
                ylims = ylims,
                 ylabel = "V [p.u]",
             xlabel = "t[s]")
        
        plot_non_gens_vh_ylims =
            make_plot_of_buses_vars_and_norm_ω_if_in_vars(
                ;sol =
                    system_sol,
                network_vars_labels =
                    model_syms,
                nodes_name = non_gens_nodes_names,
                vars = [:vh ],
                tspan = sim_timespan,
                 fmt = fmt,
                 lw = lw,
                 ylims = ylims,
                ylabel = "V [p.u]",
            xlabel = "t[s]") 
        
    end

    return include_v_θ_plot == true ?
        (; δ_a_plot, ω_a_plot, ed_dash_a_plot, eq_dash_a_plot, vr1_a_plot, vr2_a_plot, vf_tilade_a_plot, xg1_a_plot, xg2_a_plot, plot_gens_vh, plot_gens_θh, plot_non_gens_vh, plot_non_gens_θh ) : (; δ_a_plot, ω_a_plot, ed_dash_a_plot, eq_dash_a_plot, vr1_a_plot, vr2_a_plot, vf_tilade_a_plot, xg1_a_plot, xg2_a_plot, plot_gens_vh_ylims, plot_non_gens_vh_ylims)

end



function make_ode_quick_plot(
    ;sol = ode_sol,
    network_vars_labels = network_vars_labels,
    nodes_name = nodes_name,
    plot_timespan = plot_timespan,
    fmt = :pdf,
    plt_layout = (2,2),
    size = (1000, 500),
    lw = 1,
    xlabel = "t[s]")
    
    vars_plots_flux_decay_by_ext_idq = [
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ;sol =
                sol,
            network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_name,
            vars = [:δ ],
            tspan = plot_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "δ [rad]"),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                sol,
            network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_name,
            vars = [:ω ],
            tspan = plot_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "ω [p.u]"),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                sol,
            network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_name,
            vars = [:eq_dash ],
            tspan = plot_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "eq' [p.u]"),
        make_plot_of_buses_vars_and_norm_ω_if_in_vars(
            ; sol =
                sol,
            network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_name,
            vars = [:E_fd ],
            tspan = plot_timespan,
            fmt = fmt,
            lw = lw,
        ylabel = "Efd [p.u]")]
    
    return plot(
        vars_plots_flux_decay_by_ext_idq...;
        layout = plt_layout,
        size = size,
        lw = lw,
        xlabel = xlabel)


end



"""
Functions for ploting mag and angle for u_r and u_i
in sol
"""

u_mag(t, u_r, u_i)   = (t, abs(u_r + u_i *im ))

u_angle(t, u_r, u_i) = (t, angle(u_r + u_i *im ))

norm_freq(t, ω ) = (t, ω/ωs)

function make_plot_of_a_bus_volt_mag(
    ; sol = sol,
    network_vars_labels = network_vars_labels,
    bus_name = bus_name,
    vars = [:u_r, :u_i],
    tspan = (0.0, 50.0),
    lw = 1,
    fmt = :pdf)

    ur_ui_indices =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                network_vars_labels,
            bus_name = bus_name,
            vars = vars)

    plt_u_mag =  plot(sol,
                      idxs = (u_mag, 0, ur_ui_indices... ),
                      ylabel="|u|",
                      yminorticks=10,
                      lw = lw,
                      fmt=fmt ,
                      tspan = tspan,
                      labels="$(bus_name) |V|" )

end


function make_plot_of_a_bus_volt_angle(
    ; sol = sol,
    network_vars_labels = network_vars_labels,
    bus_name = bus_name,
    vars = [:u_r, :u_i],
    tspan = (0.0, 50.0),
    lw = 1,
    fmt = :pdf)

    ur_ui_indices =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                network_vars_labels,
            bus_name = bus_name,
            vars = vars)


    plt_u_angle =  plot(sol,
                        idxs =
                            (u_angle, 0, ur_ui_indices... ),
                      ylabel="θ",
                        yminorticks=10,
                        lw = lw,
                        fmt=fmt ,
                        tspan = tspan,
                        labels="$(bus_name) θ" )
    

end



function make_plot_of_buses_volt_mag(
    ; sol = sol,
    network_vars_labels =
        network_vars_labels,
    nodes_name = ["bus1", "bus2"],
    vars = [:u_r, :u_i],
    tspan = (0.0, 50.0),
    lw = 1,
    fmt = :pdf)

    tup_vars_indices =
        get_nodes_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_name,
            vars = vars )

    tup_vars_indices = [
        (u_mag, 0, idx...)
        for idx in tup_vars_indices] 

    legend_labels = ["$(bus_name) |V|"
                     for bus_name in nodes_name]
    
    legend_labels = reshape(
        legend_labels, 1,
        length(legend_labels))
    
    plt =  plot(sol, idxs = tup_vars_indices,
                ylabel="|V|",
                yminorticks=10,
                lw = lw,
                fmt=fmt ,
                tspan = tspan,                
                labels = legend_labels )

end

function make_plot_of_buses_volt_angle(
    ; sol = sol,
    network_vars_labels = network_vars_labels,
    nodes_name = ["bus1", "bus2"],
    vars = [:u_r, :u_i],
    tspan = (0.0, 50.0),
    lw =1,
    fmt = :pdf)

    tup_vars_indices =
        get_nodes_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_name,
            vars = vars )

    tup_vars_indices = [ (u_angle, 0, idx...)
                         for idx in tup_vars_indices] 

    legend_labels = ["$(bus_name) θ"
                     for bus_name in nodes_name]
    
    legend_labels = reshape(
        legend_labels, 1, length(legend_labels))
    
    plt =  plot(sol, idxs = tup_vars_indices,
                ylabel = "θ",
                yminorticks=10,
                lw = lw,
                fmt=fmt,
                tspan = tspan,                
                labels = legend_labels )

end


function make_plots_of_buses_volt_mag_with_names(
    fig_dir,
    list_sol,
    list_prefix,
    nodes_name,
    bus_type,                  
    network_vars_labels;
    tspan  = (0.0, 50.0),
    lw = 1,
    suffix = :pdf)
    
    for (a_sol, prefix) in zip(list_sol, list_prefix)
        plt = make_plot_of_buses_volt_mag(
            ; sol = a_sol,
            network_vars_labels = network_vars_labels,
            nodes_name = nodes_name,
            vars = [:u_r, :u_i],
            tspan = tspan,
            lw = lw,
            fmt = suffix)

        savefig(
            plt,
            joinpath(
                fig_dir,
                "$(prefix)-$(bus_type)-volt-mag.$(suffix)"))
        
    end
    
    return nothing
    
end



function make_plots_of_buses_volt_angle_with_names(
    fig_dir,
    list_sol,
    list_prefix,
    nodes_name,
    bus_type,                  
    network_vars_labels;
    tspan  = (0.0, 50.0),
    lw = 1,
    suffix = :pdf)
    
    for (a_sol, prefix) in
        zip(list_sol, list_prefix)
        plt = make_plot_of_buses_volt_angle(
            ; sol = a_sol,
            network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_name,
            vars = [:u_r, :u_i],
            tspan = tspan,
            lw = lw,
            fmt = suffix)

        savefig(plt, joinpath(
            fig_dir,
            "$(prefix)-$(bus_type)-volt-angle.$(suffix)"))  
    end
    
    return nothing
    
end


function get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(; sol = sol, node_syms_labels = node_syms_labels, bus_name = bus_name, vars = [:ω])

    if :ω ∈ vars
        ω_set = Set( [  :ω ] )
        
        vars_set = Set( vars )
        
        vars_ω_setdiff = setdiff( vars_set, ω_set)

        tup_vars_indices = []
        
        legend_labels = []

        if  vars_ω_setdiff != Set{Any}()
        
            tup_no_ω_vars_syms_and_indices =
                get_a_node_state_algb_vars_indices_in_syms(
                    ; node_syms_labels = node_syms_labels,
                    bus_name = bus_name,
                    vars = collect( vars_ω_setdiff ) )

            no_ω_vars_syms =
                first.( tup_no_ω_vars_syms_and_indices )

            no_ω_vars_indices =
                last.( tup_no_ω_vars_syms_and_indices )
            
            tup_no_ω_vars_indices =
                [ (0, idx[1] )
                  for idx in no_ω_vars_indices ]
            
            legend_no_ω_labels =
                ["$(bus_name) $(a_var)"
                 for a_var in no_ω_vars_syms ]
            
            push!(tup_vars_indices, tup_no_ω_vars_indices) 
        
            push!(legend_labels, legend_no_ω_labels)
            
        end
        
        ω_vars_indices =
            get_a_node_state_algb_vars_indices_in_system(
                ; network_vars_labels =
                    node_syms_labels,
                bus_name = bus_name,
                vars = [ :ω ] )

        tup_ω_vars_indices = [
            (norm_freq, 0, idx[1])
              for idx in ω_vars_indices ]
        
        legend_ω_labels =
                ["$(bus_name) $(a_var)"
                 for a_var in [ :ω ] ]
                    
        push!(tup_vars_indices, tup_ω_vars_indices)

        tup_vars_indices = [tup_vars_indices...;]
        
        push!(legend_labels, legend_ω_labels)

        legend_labels = [legend_labels...;]
        legend_labels =
            reshape(legend_labels,
                    1, length(legend_labels ))

        return (; tup_vars_indices, legend_labels )

    else

        vars_indices =
                get_a_node_state_algb_vars_indices_in_system(
                    ; network_vars_labels =
                        node_syms_labels,
                    bus_name = bus_name,
                    vars = vars )
            
        tup_vars_indices =
            [ (0, idx)
              for idx in vars_indices ]
            
            legend_labels =
                ["$(bus_name) $(a_var)"
                 for a_var in vars ]
            
            legend_labels =
                reshape(legend_labels ,
                1, length(legend_labels ))

        return (; tup_vars_indices, legend_labels )
        
    end


end



function make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(
    ; sol = sol,
    node_syms_labels = node_syms_labels,
    bus_name = bus_name,
    vars = [:ω],
    tspan = (0.0, 50.0),
    fmt = :pdf,
    lw = 1,
    ylabel = :nothing,
    bottom_margin=2Plots.mm,  # Adjust bottom margin
    left_margin=5Plots.mm,   # Adjust left margin
    right_margin=2Plots.mm,  # Adjust right margin
    top_margin=2Plots.mm     # Adjust top margin
    )

    (; tup_vars_indices,
     legend_labels ) =
        get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(; sol = sol, node_syms_labels = node_syms_labels, bus_name = bus_name, vars = vars )
    
    return ylabel == :nothing ?  plot(
        sol,
        idxs = tup_vars_indices,
        # idxs = (0, vars_indices... ),
        # ylabel ="$(vars[1])",
        yminorticks = 10,
        fmt = fmt ,
        lw = lw,
        tspan = tspan,
        labels = legend_labels,
        bottom_margin=bottom_margin, 
        left_margin=left_margin,
        right_margin=right_margin,
        top_margin=top_margin ) : plot(
            sol,
            idxs = tup_vars_indices,
            # idxs = (0, vars_indices... ),
            # ylabel ="$(vars[1])",
            yminorticks = 10,
            fmt = fmt ,
            lw = lw,
            tspan = tspan,
            labels = legend_labels,
            bottom_margin=bottom_margin, 
            left_margin=left_margin,
            right_margin=right_margin,
            top_margin=top_margin )
end



function make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    ; sol = sol,
    network_vars_labels =
        network_vars_labels,
    nodes_name = ["bus1", "bus2"],
    vars = [:ω],
    tspan = (0.0, 50.0),
    fmt = :pdf,
    xlabel = :nothing,
    ylabel = :nothing,
    lw = 1,
    ylims = :auto,
    bottom_margin=2Plots.mm,  # Adjust bottom margin
    left_margin=5Plots.mm,   # Adjust left margin
    right_margin=2Plots.mm,  # Adjust right margin
    top_margin=2Plots.mm     # Adjust top margin
    )

 
    tup_vars_indices_legend_labels =  map( (bus_name)-> get_plot_idxs_and_legend_of_a_bus_vars_and_norm_ω_if_in_vars(; sol = sol, node_syms_labels = network_vars_labels, bus_name = bus_name, vars = vars), nodes_name )

    tup_vars_indices =
        first.( tup_vars_indices_legend_labels )

    tup_vars_indices = [tup_vars_indices...;]

    legend_labels =
        last.(tup_vars_indices_legend_labels )

    legend_labels = [ legend_labels...; ]

    legend_labels =
                reshape(legend_labels ,
                1, length(legend_labels ))
    
    return ylabel == :nothing ? plot(
        sol,
        idxs = tup_vars_indices,
        yminorticks = 10,
        fmt = fmt ,
        lw = lw,
        ylims = ylims,
        tspan = tspan,
        labels = legend_labels,
        bottom_margin=bottom_margin, 
        left_margin=left_margin,
        right_margin=right_margin,
        top_margin=top_margin 
    ) : plot(
        sol,
        idxs = tup_vars_indices,
        yminorticks = 10,
        fmt = fmt,
        lw = lw,
        ylims = ylims,
        tspan = tspan,
        ylabel = ylabel,
        labels = legend_labels,
        bottom_margin=bottom_margin, 
        left_margin=left_margin,
        right_margin=right_margin,
        top_margin=top_margin ) 
end



# function make_plot_of_a_bus_vars_and_norm_ω_if_in_vars(; sol = sol, node_syms_labels = node_syms_labels, bus_name = bus_name, vars = [:ω], tspan = (0.0, 50.0), fmt = :pdf)

#     if :ω ∈ vars
#         ω_set = Set( [  :ω ] )
        
#         vars_set = Set( vars )
        
#         vars_ω_setdiff = setdiff( vars_set, ω_set)

#         if  vars_ω_setdiff != Set{Any}()
        
#             tup_no_ω_vars_syms_and_indices =
#                 get_a_node_state_algb_vars_indices_in_syms(
#                     ; node_syms_labels = node_syms_labels,
#                     bus_name = bus_name,
#                     vars = collect( vars_ω_setdiff ) )

#             no_ω_vars_syms =
#                 first.( tup_no_ω_vars_syms_and_indices )

#             no_ω_vars_indices =
#                 last.( tup_no_ω_vars_syms_and_indices )
            
#             tup_no_ω_vars_indices =
#                 [ (0, idx[1] ) for idx in no_ω_vars_indices ]
            
#             legend_no_ω_labels =
#                 ["$(bus_name) $(a_var)"
#                  for a_var in no_ω_vars_syms ]
            
#             legend_no_ω_labels =
#         reshape(legend_no_ω_labels ,
#                 1, length(legend_no_ω_labels ))

#             plt_obj =  plot(sol,
#                       idxs = tup_no_ω_vars_indices,
#                       # idxs = (0, vars_indices... ),
#                       # ylabel="$(vars[1])",
#                       yminorticks=10,
#                       fmt=fmt ,
#                       tspan = tspan,
#                       labels= legend_no_ω_labels  )
#         end
        

#         ω_vars_indices = get_a_node_state_algb_vars_indices_in_system(; network_vars_labels = node_syms_labels, bus_name = bus_name, vars = [ :ω ] )

#         tup_ω_vars_indices = [
#             (norm_freq, 0, idx[1])
#               for idx in ω_vars_indices ]
        
#         legend_ω_labels =
#                 ["$(bus_name) $(a_var)"
#                  for a_var in [ :ω ] ]
            
#             legend_ω_labels =
#         reshape(legend_ω_labels ,
#                 1, length(legend_ω_labels ))
        
#         plt =  plot!(plt_obj, sol,
#                       idxs = tup_ω_vars_indices,
#                       # idxs = (0, vars_indices... ),
#                       # ylabel="$(vars[1])",
#                       #yminorticks=10,
#                       #fmt=fmt ,
#                       #tspan = tspan,
#                          labels= legend_ω_labels  )

#     else


#         vars_indices =
#                 get_a_node_state_algb_vars_indices_in_system(
#                     ; network_vars_labels = node_syms_labels,
#                     bus_name = bus_name,
#                     vars = vars )
            

#         tup_vars_indices =
#                 [ (0, idx) for idx in vars_indices ]
            
#             legend_labels =
#                 ["$(bus_name) $(a_var)"
#                  for a_var in vars ]
            
#             legend_labels =
#         reshape(legend_labels ,
#                 1, length(legend_labels ))

#             plt_obj =  plot(sol,
#                        idxs = tup_vars_indices,
#                        # idxs = (0, vars_indices... ),
#                        # ylabel ="$(vars[1])",
#                        yminorticks = 10,
#                        fmt = fmt ,
#                        tspan = tspan,
#                        labels = legend_labels  )        
        
#     end


# end



function make_plot_of_a_bus_var(
    ; sol = sol,
    network_vars_labels =
        network_vars_labels,
    bus_name = bus_name,
    vars = [:ω],
    tspan = (0.0, 50.0),
    lw = 1,
    fmt = :pdf)

    vars_indices =
        get_a_node_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                network_vars_labels,
            bus_name = bus_name,
            vars = vars)

    plt_u_mag =  plot(
        sol,
        idxs = (0, vars_indices... ),
        ylabel="$(vars[1])",
        yminorticks=10,
        lw = lw,
        fmt=fmt ,
        tspan = tspan,
        labels="$(bus_name) $(vars[1])" )

end



function make_plot_of_buses_var(
    ; sol = sol,
    network_vars_labels =
        network_vars_labels,
    nodes_name = ["bus1", "bus2"],
    vars = [:ω],
    tspan = (0.0, 50.0),
    lw = 1,
    fmt = :pdf)

    tup_vars_indices =
        get_nodes_state_algb_vars_indices_in_system(
            ; network_vars_labels =
                network_vars_labels,
            nodes_name = nodes_name,
            vars = vars )

    tup_vars_indices = [
        (0, idx[1])
        for idx in
            tup_vars_indices] 

    legend_labels = [
        "$(bus_name) $(vars[1])"
        for bus_name in nodes_name]
    
    legend_labels = reshape(
        legend_labels, 1, length(legend_labels))
    
    plt =  plot(
        sol, idxs = tup_vars_indices,
        ylabel="$(vars[1])",
        yminorticks=10,
        tspan = tspan,
        lw = lw,
        fmt = fmt,
        labels = legend_labels )

end


function make_a_plot_for_syms(
    sol, syms_list,
    network_vars_labels ;
    tspan = (0.0, 50.0),
    lw = 1,
    fmt = :pdf )

    plt_plots = []

    for a_state_sym in syms_list

        comp_plot = plot(
            sol;
            idxs = syms_containing(
                network_vars_labels,
                a_state_sym;
                syms = true ),
            ylabel=String(a_state_sym),
            yminorticks=10,
            lw = lw,
            fmt=fmt ,
            tspan = tspan  )
        
        push!(plt_plots, comp_plot )
    end


    if length(plt_plots )     == 1
        
        plt_layout = (1,1)
        
    elseif length(plt_plots ) == 2
        
        plt_layout = (1,2)
        
    elseif length(plt_plots ) == 3
        
        plt_layout = (1,3)
        
    elseif length(plt_plots ) == 4
        
        plt_layout = (2,2)
        
    else
        
        plt_layout = (2,3)
        
    end            

    return plot(
        plt_plots...;
        layout = plt_layout,
        size = (1000, 500),
        lw = lw,
        fmt = fmt,
        xlabel = "t[s]")    

    # return plt    

end



function plot_comp_gen_plant(
    sol_cb_comp_gen_plant;
    tspan =  (0.0, 20.0),
    lw = 1,
    fmt=:pdf)

    (test_plant_sol,
     test_plant_syms,
     gen_states_syms,
     comp_states_syms) =
         sol_cb_comp_gen_plant

    # gen_states_syms
    plt_gen = plot(
        test_plant_sol;
        idxs=[ ind for (ind, sym) in
                  enumerate(test_plant_syms)
                  if sym ∈ gen_states_syms ],
        lw = lw,
        fmt=fmt,
        tspan = tspan  )

    
    plt_comp = plot(
        test_plant_sol;
        idxs=[ ind for (ind, sym) in
                  enumerate(test_plant_syms)
                  if sym ∈ comp_states_syms ],
        lw = lw,
        fmt=fmt,        
        tspan = tspan  )

    return plot(plt_gen,
               plt_comp;
            layout   = (1,2),
            size     = (1000, 500),
            lw       = lw,
             fmt=fmt,
              xlabel = "t[s]")    

    # return plt    

end


function plot_comps_gen_plant(
    sol_comps_gen_plant;
    lw = 1,
    fmt=:pdf,
    tspan =  (0.0, 20.0) )

    # Note: comps_states_syms is a list

    (test_plant_sol,
     test_plant_syms,
     gen_states_syms,
     comps_states_syms) =
         sol_comps_gen_plant

    # gen_states_syms
    plt_gen = plot(
        test_plant_sol;
        idxs=[ ind for (ind, sym) in
                  enumerate(test_plant_syms)
                  if sym ∈ gen_states_syms ],
        yminorticks=10,
        lw = lw,
        fmt=fmt,
        tspan = tspan  )

    plt_plots = [plt_gen]

    for comp_states_syms in comps_states_syms

        comp_plot = plot(
            test_plant_sol; idxs=[
                ind for
                    (ind, sym) in
                    enumerate(test_plant_syms)
                    if sym ∈ comp_states_syms ],
            yminorticks=10, fmt=fmt,lw=lw,tspan = tspan  )
        
        push!(plt_plots, comp_plot )
    end


    if length(plt_plots )     == 1
        
        plt_layout = (1,1)
        
    elseif length(plt_plots ) == 2
        
        plt_layout = (1,2)
        
    elseif length(plt_plots ) == 3
        
        plt_layout = (1,3)
        
    elseif length(plt_plots ) == 4
        
        plt_layout = (2,2)
        
    else
        
        plt_layout = (2,3)
        
    end            

    return plot(
        plt_plots...;
        layout = plt_layout,
        size = (1000, 500),lw=lw,fmt=fmt,xlabel = "t[s]") 

    # return plt    

end



function selected_plot_comps_gen_plant(
    sol_comps_gen_plant;
    tspan =  (0.0, 20.0),
    lw = 1,
    fmt = :pdf,
    selected = (
        gen = [:nothing],
        gov = [:nothing],
        avr = [:nothing],
        pss = [:nothing]) )

    # Note: comps_states_syms is a list

    (test_plant_sol,
     test_plant_syms,
     gen_states_syms,
     comps_states_syms) =
         sol_comps_gen_plant

    # gen_states_syms

    if selected.gen[1] == :nothing
        
        plt_gen = plot(
            test_plant_sol;
            idxs=[ ind for (ind, sym) in
                      enumerate(test_plant_syms)
                      if sym ∈ gen_states_syms ],
            yminorticks=10,
            lw = lw,
            fmt=fmt,
            tspan = tspan  )
    else
        
        plt_gen = plot(
            test_plant_sol;
            idxs=[ ind for (ind, sym) in
                      enumerate(test_plant_syms)
                      if sym ∈ selected.gen ],
            yminorticks=10,
            lw=lw,
            fmt=fmt,
            tspan = tspan  )
    end

    plt_plots = [plt_gen]

    list_selected_tuple = [
        ( getproperty(selected, pn),
          getproperty( comps_states_syms, pn) )
        for pn in propertynames(selected)
            if  pn ∈  propertynames(comps_states_syms) ]
    
        
    for (list_sel_syms, comp_states_syms)  in
        list_selected_tuple
        if list_sel_syms[1] == :nothing
            
            push!(plt_plots, plot(
                test_plant_sol; idxs=[
                    ind for (ind, sym) in
                        enumerate(test_plant_syms)
                        if sym ∈ comp_states_syms ],
                lw=lw,
                fmt=fmt,
                tspan = tspan  ) )
        else
            
            push!(plt_plots, plot(
                test_plant_sol;
                idxs=[ ind
                       for (ind, sym) in
                           enumerate(test_plant_syms)
                           if sym ∈ list_sel_syms ],
                lw=lw,
                fmt=fmt,
                tspan = tspan  ) )            
        end
        
    end

    if length(plt_plots )     == 1
        
        plt_layout = (1,1)    
    
    elseif length(plt_plots )     == 2
        
        plt_layout = (1,2)
        
    elseif length(plt_plots ) == 3
        
        plt_layout = (1,3)
        
    elseif length(plt_plots ) == 4
        
        plt_layout = (2,2)
        
    else
        
        plt_layout = (2,3)
        
    end            

    return plot(
        plt_plots...;
        layout = plt_layout,
        size = (1000, 500),
        lw = lw,
        fmt=fmt,
        xlabel = "t[s]")    

    # return plt    

end


#------------------------------------------------
#------------------------------------------------


function get_vars_to_plot_by_syms(
    vec_selected_syms::Vector{Symbol} )

    return join(String.(vec_selected_syms ),"_")
    
end


"""
This function returns a vector of strings from a list
of list of symbols of variables to be plotted.

list_selected_syms_list = [[:ω, :δ, :ed_dash, :eq_dash],
                           [:xg1, :xg2, :xg3, :τm_tilade],
                           [:vm, :vr1, :vr2, :vf_tilade],
                           [:phat_in, :v_ref]]


list_of_plot_names = get_vars_to_plot_by_syms(
    list_of_lists_selected_syms)

"""
function get_list_vars_to_plot_by_syms(
    vec_of_vec_selected_syms::Vector{Vector{Symbol}} )


    # return map(get_vars_to_plot_by_syms, vec_of_vec_selected_syms)
    
    return map(
        arg -> join(String.( arg ),"_"),
        vec_of_vec_selected_syms )
    
end


"""
`simulate_cases(sd_dyn_sim_func, list_cases; static_case = case_IEEE_5_Bus_static )`

This function simulates several cases in of dynamic power models
in list_cases. The dynamic function for simulation is
`sd_dyn_sim_func` .

The static_case is used to obtain power flow for intialisation
of the states of the dynamic model.

It returns solution for each of the cases in a list.

list_cases = [case_IEEE_5_Bus_dynamic_plant_SM_cb_idq_PQ_Const_I,
              case_IEEE_5_Bus_dynamic_plant_SM_cb_v6_PQ_Const_I,
              case_IEEE_5_Bus_dynamic_plant_SM_cb_idq_PQ_Const_P,
              case_IEEE_5_Bus_dynamic_plant_SM_cb_v6_PQ_Const_P]

list_sol = simulate_cases(sd_dynamics_simulation_external_f_t, list_cases; static_case = case_IEEE_5_Bus_static )

"""
function simulate_cases(
    sd_dyn_sim_func,
    list_cases;
    static_case = case_IEEE_5_Bus_static,
    list_sol = [],
    timespan = 10 )

    # return map((arg) -> sd_dyn_sim_func(; static_case = static_case, dynamics_case = arg ), list_cases)

    for a_case in list_cases

        case_name = String(nameof(a_case))
        println("Simulating $(case_name)")
        
        push!(list_sol, sd_dyn_sim_func(; static_case = static_case, dynamics_case = a_case, timespan = timespan ) )
    end
    
    
end


"""
This function makes plots for selected symbols for a list of solutions:

network_vars_labels = get_network_vars_labels(NetworkData(case_IEEE_5_Bus_dynamic_plant_SM_cb_idq_PQ_Const_I()...))


list_sol = [sol_IEEE_5_Bus_plant_SM_cb_idq_cntrl_and_PQ_Const_I,
            sol_IEEE_5_Bus_plant_SM_cb_v6_cntrl_and_PQ_Const_I,
            sol_IEEE_5_Bus_plant_SM_cb_idq_cntrl_and_PQ_Const_P,
            sol_IEEE_5_Bus_plant_SM_cb_v6_cntrl_and_PQ_Const_P]

list_prefix = String.(nameof.(list_cases))

or

list_prefix = ["IEEE_5_Bus_plant_SM_cb_idq_cntrl_and_PQ_Const_I",
               "IEEE_5_Bus_plant_SM_cb_v6_cntrl_and_PQ_Const_I",
               "IEEE_5_Bus_plant_SM_cb_idq_cntrl_and_PQ_Const_P",
               "IEEE_5_Bus_plant_SM_cb_v6_cntrl_and_PQ_Const_P"]

list_selected_syms_list = [[:ω, :δ, :ed_dash, :eq_dash],
                           [:xg1, :xg2, :xg3, :τm_tilade],
                           [:vm, :vr1, :vr2, :vf_tilade],
                           [:phat_in, :v_ref]]

make_plots_with_names(fig_folder, list_sol, list_prefix,
                      list_selected_syms_list,
                      network_vars_labels; tspan = (0.0, 50.0),
                      suffix = :pdf)
"""
function make_plots_with_names(fig_dir,
                               list_sol,
                               list_prefix,
                               list_selected_syms_list,
                               network_vars_labels;
                               tspan  = (0.0, 50.0),
                               suffix = :pdf,
                               lw=1)

    plotted_vars_names = get_list_vars_to_plot_by_syms(
        list_selected_syms_list)
    
    for (a_sol, prefix) in zip(list_sol, list_prefix)
        for (syms_list, vars_names) in zip(list_selected_syms_list, plotted_vars_names)
            
            plt = make_a_plot_for_syms(a_sol, syms_list,
                                       network_vars_labels;
                                       tspan = tspan,
                                       fmt = suffix,
                                       lw=lw)            
            savefig(plt,
                    joinpath(fig_dir,
                             "$(prefix)-$(vars_names).$(suffix)"))
        end
    end
    
    return nothing
    
end


function make_a_plot_with_names(
    sol,
    prefix,
    selected_syms_list,
    ;fig_dir = fig_dir,
    network_vars_labels = network_vars_labels,
    tspan  = (0.0, 50.0),
    suffix = :pdf,
    lw=1)

    vars_names =
        get_vars_to_plot_by_syms(selected_syms_list)
    
    plt = make_a_plot_for_syms(sol,
                               selected_syms_list,
                               network_vars_labels;
                               tspan = tspan,
                               lw=lw,
                               fmt = suffix)
    
    savefig(plt,
            joinpath(fig_dir,
                     "$(prefix)-$(vars_names).$(suffix)"))
    
    return nothing
    
end

#-----------------------------------------------------


function make_plots_for_industrial_model(
    ; case_fun = dynamics_case,
    sim_model_type = sim_model_type,
    sim_sol = sim_sol,
    para_net_names_labels_syms =
        para_net_names_labels_syms,
    tspan = plot_tspan,
    base_dir = nothing,
    algr_name = algr_name,
    lw=1,
    fmt=:pdf)

    (;net_class_names,
     net_states_and_var_labels,
     industrial_model_sym) =
        para_net_names_labels_syms

    (;network_bus_names,
     non_gens_bus_names,
     gens_bus_names,
     Loads_bus_names,
     Trans_bus_names) =
         net_class_names

    (;net_bus_volts_labels,
     gens_nodes_pure_states_labels,
     gens_nodes_stab_states_label,
     gens_nodes_algebraic_and_states_labels) =
        net_states_and_var_labels


    sim_plt_ω_δ_ed_dash_eq_dash = make_a_plot_for_syms(
        sim_sol ,
        [:ω, :δ, :ed_dash, :eq_dash],
        industrial_model_sym;
        lw=lw,
        fmt=fmt,        
        tspan = tspan )

    
    sim_plt_gov = make_a_plot_for_syms(
        sim_sol,
        [:xg1, :xg2, :τm_tilade],
        industrial_model_sym;
        lw=lw,
        fmt=fmt,
        tspan =  tspan )

    
    sim_plt_avr = make_a_plot_for_syms(
        sim_sol,
        [:vm, :vr1, :vr2, :vf_tilade],
        industrial_model_sym;
        lw=lw,
        fmt=fmt,
        tspan =  tspan )
    
      
    sim_gens_volt_mag = make_plot_of_buses_volt_mag(
        ; sol = sim_sol,
        network_vars_labels = industrial_model_sym,
        nodes_name = gens_bus_names,
        vars = [:u_r, :u_i],
        lw=lw,
        fmt=fmt,
        tspan = tspan )

    
    sim_gens_volt_angle =  make_plot_of_buses_volt_angle(
        ; sol = sim_sol,
        network_vars_labels = industrial_model_sym,
        nodes_name = gens_bus_names,
        vars = [:u_r, :u_i],
        lw=lw,
        fmt=fmt,
        tspan = tspan )

    
    sim_non_gens_volt_mag = make_plot_of_buses_volt_mag(
        ; sol = sim_sol,
        network_vars_labels = industrial_model_sym,
        nodes_name = non_gens_bus_names,
        vars = [:u_r, :u_i],
        lw=lw,
        fmt=fmt,
        tspan = tspan )
    
    sim_non_gens_volt_angle =
        make_plot_of_buses_volt_angle(
            ;sol = sim_sol,
            network_vars_labels =
                industrial_model_sym,
            nodes_name = non_gens_bus_names,
            vars = [:u_r, :u_i],
            lw=lw,
            fmt=fmt,
            tspan = tspan )

    list_plts_volt_mag_ang = [
        sim_gens_volt_mag,
        sim_gens_volt_angle,
        sim_non_gens_volt_mag,
        sim_non_gens_volt_angle ]

    list_plts_gov_avr = [
        sim_plt_gov,
        sim_plt_avr]

    #-----------------------------------------------

        plts_volt_mag_ang = plot(
            list_plts_volt_mag_ang... ;
            layout = (2,2),
            size = (1000, 500),
            lw=lw,
            fmt=fmt,
            xlabel = "t[s]" )

        plts_gov_avr = plot(
            list_plts_gov_avr... ;
            layout = (2,1),
            size = (1000, 500),
            lw=lw,
            fmt=fmt,
            xlabel = "t[s]" )
    
    #-----------------------------------------------

    # folder

    sim_model_type = sim_model_type 
    
    case_name = String(nameof(case_fun))
    
    results_dir = joinpath(sim_model_type, case_name )

    if base_dir==nothing
        
        res_folder, csv_folder, fig_folder, sol_folder =
            make_results_dir(
                ; base_dir=@__DIR__,
                res_dir = results_dir )
        
    else
        
        res_folder, csv_folder, fig_folder, sol_folder =
            make_results_dir(
                ; base_dir=base_dir,
                res_dir = results_dir )
        
    end
    
    #-----------------------------------------------
    
    suffix = :pdf
        
    list_plots_sym = [:plt_ω_δ_ed_dash_eq_dash,
                      :volt_mag_ang,
                      :gov_avr]
    
    list_sim_plots = [
    sim_plt_ω_δ_ed_dash_eq_dash,
    plts_volt_mag_ang,
        plts_gov_avr]

    for ( a_sim_plot, a_plot_sym ) in
        zip( list_sim_plots, list_plots_sym )
        savefig(
            a_sim_plot,
            joinpath(
                fig_folder,
                "$(algr_name)-$(sim_model_type)-$(a_plot_sym).$(suffix)") )
    end
            
    return list_sim_plots
    
    
end

#-----------------------------------------------------
#-----------------------------------------------------


"""
This function saves solution of cases in csv files

list_sol = [sol_IEEE_5_Bus_plant_SM_cb_idq_cntrl_and_PQ_Const_I,
            sol_IEEE_5_Bus_plant_SM_cb_v6_cntrl_and_PQ_Const_I,
            sol_IEEE_5_Bus_plant_SM_cb_idq_cntrl_and_PQ_Const_P,
            sol_IEEE_5_Bus_plant_SM_cb_v6_cntrl_and_PQ_Const_P]

list_prefix = ["IEEE_5_Bus_plant_SM_cb_idq_cntrl_and_PQ_Const_I",
               "IEEE_5_Bus_plant_SM_cb_v6_cntrl_and_PQ_Const_I",
               "IEEE_5_Bus_plant_SM_cb_idq_cntrl_and_PQ_Const_P",
               "IEEE_5_Bus_plant_SM_cb_v6_cntrl_and_PQ_Const_P"]

make_csv_from_sols(csv_folder,
                   list_sol,
                   list_prefix;
                   suffix = :csv)

"""
function make_csv_from_sols(csv_dir,
                               list_sol,
                               list_prefix;
                               suffix = :csv)
    
    for (a_sol, prefix) in zip(list_sol, list_prefix)
        CSV.write(joinpath(csv_dir,
                           "$(prefix).$(suffix)"),
                  DataFrame(a_sol))
    end
    
    return nothing
    
end

"""
https://docs.sciml.ai/DiffEqDocs/stable/features/io/

"""
function save_sols(
    sol_dir, list_sol,
    list_prefix;
    suffix = :jld2 )

    for ( case_name, sol) in zip(list_prefix, list_sol)
        
        sol_filename = joinpath(
            sol_dir,"$(case_name).$(suffix)")
        
        @save  "$(sol_filename)" sol
        
        
    end
    
end


#########################################################
# ------------------------------------------------------
#  Plot functions
# ------------------------------------------------------
#########################################################


function make_a_plot_for_syms(
    sol; netd = netd, syms_list =
        [:ω, :δ, :ed_dash, :eq_dash],
    tspan = (0.0, 50.0),
    lw=1,
    fmt=:pdf)

    plt_plots = []

    for a_state_sym in syms_list

        comp_plot = plot(
            sol; idxs=syms_containing(
                netd, a_state_sym ),
            ylabel=String(a_state_sym),
            yminorticks=10,
            lw=lw,
            fmt=fmt,
            tspan = tspan )
        
        push!(plt_plots, comp_plot )
    end


    if length(plt_plots )     == 1
        
        plt_layout = (1,1)
        
    elseif length(plt_plots ) == 2
        
        plt_layout = (1,2)
        
    elseif length(plt_plots ) == 3
        
        plt_layout = (1,3)
        
    elseif length(plt_plots ) == 4
        
        plt_layout = (2,2)
        
    else
        
        plt_layout = (2,3)
        
    end            

    return plot(plt_plots...;
               layout = plt_layout,
               size = (1000, 500),
               lw = lw,
               xlabel = "t[s]")
        

end



