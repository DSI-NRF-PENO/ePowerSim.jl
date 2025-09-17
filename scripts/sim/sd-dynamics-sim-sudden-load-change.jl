# ####################################################
# using Pkg
# PowerSystemsCoSim_folder = joinpath(@__DIR__,"../..")
# cd(PowerSystemsCoSim_folder)
# Pkg.activate(PowerSystemsCoSim_folder)

#---------------------------------------------------
# global settings
#---------------------------------------------------

freq = 60

Ωb = 2 * pi * freq

ωs = Ωb 

ω_ref0 = ωs

#---------------------------------------------------

using Pkg

using ePowerSim

#---------------------------------------------------

package_dir = pkgdir(ePowerSim)


data =
    joinpath( package_dir, "data")

data_dir = joinpath(data,
             "converted-data")
src_dir =
    joinpath(package_dir, "src")

components_libs_dir =
    joinpath(src_dir,
             "components-lib")

script_dir = @__DIR__



#---------------------------------------------------
## simulation case
#---------------------------------------------------

# case_name = "case14"

case_name = "case9"

timespan  = 20.0

sim_type = "sudden-load-pertubation"

json_net_data_by_components_file =
    "net-static-data-avr-sauer-gov-sauer.json"

#---------------------------------------------------

function sim_sudden_load_change(
    ;case_name,
    json_net_data_by_components_file,
    sim_type,
    timespan,
    script_dir,
    components_libs_dir,
    data_dir,
    json_case_dir)

    # case_name = "case14"

    # case_name = "case9"

    # timespan  = 20.0

    # sim_type = "sudden-load-pertubation"

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

    json_net_data_by_components_file =
        "net-static-data-avr-sauer-gov-sauer.json"

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

    dt            = 0.01

    Δt            = 1.0 / 2^(4)

    tspan         = (0.0, timespan)

    sim_timespan  = (0.0, timespan)

    plot_timespan = (0.0, timespan)


    #---------------------------------------------------
    # base setting and some booleans 
    #---------------------------------------------------

    basekV = 1.0

    use_pu_in_PQ = true

    line_data_in_pu = true

    #---------------------------------------------------
    ## solvers and settings
    #---------------------------------------------------

    use_init_u0 = false

    use_nlsolve = false

    pf_alg        = NewtonRaphson()

    #--------------------   

    ode_alg       = Rodas4()

    dae_alg       = IDA()

    abstol        = 1e-12

    reltol        = 1e-12

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
    ## DAE system dyamanics simulation
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

    bus_no_or_bus_name = 4

    bus_no_or_bus_name = bus_name = "bus4"

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
        # deepcopy(
        #     system_integrator.p.generic_model_dynamics_para[
        #     var_idx] )

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

    tstop1 = time_final/10
    tstop2 = time_final/8
    tstop3 = time_final/6
    tstop4 = time_final/2
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


    # pertubation_factor = 1.0

    # pertubation_tstop = tstop1

    # pertubation_by_itegrator(
    #     var_normal_value,
    #     pertubation_factor,
    #     pertubation_tstop,
    #     var_idx,
    #     system_integrator;
    #     parameter_df,
    #     change_in_parameter_dict,
    #     target_parameter_sym )


    add_tstop!(system_integrator,
               tstop1 )

    step!(system_integrator)

    push!(
        parameter_df,
        tuple(
            [[system_integrator.t];
             system_integrator.p.generic_model_dynamics_para]...
                 ))

    system_sol =
        system_integrator.sol


    # sim_timespan = (0 , tstop1)

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

    pertubation_factor = 1.10

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

    pertubation_factor = 1.0

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

    # Perturb the Power at bus 4

    pertubation_count = 4

    pertubation_factor = 1.0

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


    fractional_digits = 6

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


    return get_guick_single_vars_plots_dae_or_ode_sol(
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

end


#---------------------------------------------------
#---------------------------------------------------


# auxilliary_results_plot =
#     (; plt_gens_ph =
#     plot(sol_auxilliary_results.t,
#          sol_auxilliary_results.gens_ph' ),
     
#      plt_gens_qh =
#     plot(sol_auxilliary_results.t,
#          sol_auxilliary_results.gens_qh' ),

#      plt_gens_ed_dash =
#     plot(sol_auxilliary_results.t,
#          sol_auxilliary_results.gens_ed_dash' ) )
