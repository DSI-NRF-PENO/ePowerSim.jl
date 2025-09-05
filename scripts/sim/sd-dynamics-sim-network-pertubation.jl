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

timespan  = 50.0

sim_type = "network-pertubation"

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

json_net_data_by_components_file =
    "net-static-data-avr-sauer-gov-sauer.json"



json_case_dir =
    joinpath(
        data_dir,
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


sd_dynamics_sim_csv_filename =
    joinpath(
        results_dir,
        "$(case_name)-" *
            "$(sim_type)-" *
            "sim-sd-dynamics.csv")


"""

#---------------------------------------------------
# Simulation Period
#---------------------------------------------------

timespan      = 50.0

time_start    = 0.0

time_final    = timespan

dt            = 0.0001

Δt            = 1.0 / 2^(4)

tspan         = (0.0, timespan)

sim_timespan  = (0.0, timespan)

plot_timespan = (0.0, timespan)

#---------------------------------------------------
## solvers and settings
#---------------------------------------------------

use_init_u0 = false

use_nlsolve = false

pf_alg        = NewtonRaphson()

#---------------------------------------------------

ode_alg       = Rodas4()

# ode_alg       = ImplicitMidpoint()

dae_alg       = IDA()

abstol        = 1e-12

reltol        = 1e-12

#---------------------------------------------------
# base setting and some booleans 
#---------------------------------------------------

basekV = 1.0

use_pu_in_PQ = true

line_data_in_pu = true

#---------------------------------------------------
## fault data
#---------------------------------------------------

on_fault_time = 9.0

clear_fault_time = 9.02

with_faults = false

list_fault_point_from_node_a = [0.01]

list_fault_resistance = [0.001]

list_no_line_circuit =  [1]

list_edges_to_have_fault = [ 2 ]

clear_fault_selection_list = [1]

fault_status =
    (no_fault = 0,
     on_fault = 1,
     clear_fault = 2,
     partial_clear_fault = 3)

system_fault_status =
    [ fault_status.no_fault]


list_system_status =
    [:pre_fault_state,
     :fault_state,
     :post_fault_state]

system_status =
    :pre_fault_state

#---------------------------------------------------
## system fault states
#---------------------------------------------------

system_status = :pre_fault_state,
system_status = :fault_state
system_status = :post_fault_state

case_name = "case9"

"""

#---------------------------------------------------
# case model 
#---------------------------------------------------


timespan   = 50

#-------------------------------------------
# network pertubation 
#-------------------------------------------

network_fault_pertubation_sol =
    get_generic_network_fault_pertubation_by_cb(
        ;case_name = "case9",
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

        with_faults =
            false,
        use_state_in_on_clear_fault =
            false,
        return_extended_results =
            false,

        json_net_data_by_components_file =
            json_net_data_by_components_file,
        components_libs_dir =
            components_libs_dir,
        data_dir =
            data_dir )

save_network_pertubation_sim_plot(
    "case9";
    sim_timespan = (0, timespan),
    figure_dir,
    sim_type = "network-pertubation",
    line_in_fault_name = "line-8",
    network_fault_pertubation_sol...)

#---------------------------------------------------
# Save results to files
#---------------------------------------------------

fractional_digits = 6

sd_dynamics_sim_df =
    DataFrame(network_fault_pertubation_sol)

sd_dynamics_sim_df[!, :] =
    round.(
        sd_dynamics_sim_df[:, :],
        digits=fractional_digits)


sd_dynamics_sim_csv_filename =
        joinpath(results_dir,
                 "$(case_name)-" *
                     "$(sim_type)-" *
                     "sim-sd-dynamics.csv")

CSV.write(sd_dynamics_sim_csv_filename,
          sd_dynamics_sim_df )

#---------------------------------------------------
## test functions
#---------------------------------------------------

"""
generic_network_fault_pertubation_wt_state_plot =
    get_guick_group_vars_plots_dae_or_ode_sol(
        ;include_v_θ_plot =
            false,
         sim_timespan =
             (0, timespan),
         get_generic_network_fault_pertubation(
             ;case_name = "case9",
             timespan   = timespan,

             on_fault_time = 5.0,
             clear_fault_time = 6.0,

             list_fault_point_from_node_a = [0.3],
             list_fault_resistance = [0.001],
             list_no_line_circuit =  [1],

             list_edges_to_have_fault = [ 8 ],
             clear_fault_selection_list = [1],

             basekV = 1.0,    
             use_pu_in_PQ = true,
             line_data_in_pu = true,

             with_faults =
                 false,
             use_state_in_on_clear_fault =
                 true,
             return_extended_results =
                 false,

             json_net_data_by_components_file =
                 json_net_data_by_components_file,
             components_libs_dir =
                 components_libs_dir,
             data_dir =
                 data_dir )...)

generic_network_fault_pertubation_plot =
    get_guick_group_vars_plots_dae_or_ode_sol(
        ;include_v_θ_plot =
            false,
         sim_timespan =
             (0, timespan),
         get_generic_network_fault_pertubation(
             ;case_name = "case9",
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

             with_faults =
                 false,
             use_state_in_on_clear_fault =
                 false,
             return_extended_results =
                 false,

             json_net_data_by_components_file =
                 json_net_data_by_components_file,
             components_libs_dir =
                 components_libs_dir,
             data_dir =
                 data_dir )...)

status_steady_state_parameters =
    get_status_steady_state_parameters(
        net_data_by_components_file;
        components_libs_dir =
            components_libs_dir,
        basekV = 1.0,    
        use_pu_in_PQ = true,
        line_data_in_pu = true,

        use_init_u0 = false,    
        use_nlsolve = false,

        pf_alg = NewtonRaphson(),

        abstol = 1e-12,    
        reltol = 1e-12,

        on_fault_time = 5.0,
        clear_fault_time = 5.2,

        list_fault_point_from_node_a = [0.3],
        list_fault_resistance = [0.001],
        list_no_line_circuit =  [1],

        list_edges_to_have_fault = [ 2 ], 
        clear_fault_selection_list = [ 1 ],

        with_faults = false)


fault_a_status_steady_state_data =
    get_a_status_steady_state_data(
        :pre_fault_state ; # system_status
        with_faults = true,

        net_data_by_components_file =
            net_data_by_components_file,
        components_libs_dir =
            components_libs_dir,

        timespan      = 10.0,
        on_fault_time = 5.0,
        clear_fault_time = 5.2,    

        list_fault_point_from_node_a = [0.3],
        list_fault_resistance = [0.001],
        list_no_line_circuit =  [4],

        list_edges_to_have_fault = [ 2 ],
        clear_fault_selection_list = [1],

        basekV = 1.0,    
        use_pu_in_PQ = true,
        line_data_in_pu = true)

#-------------------------------

system_simulation_parameters_wt_faults =
    get_system_simulation_parameters_wt_faults(
        net_data_by_components_file,
        # system_status
        ;
        components_libs_dir =
            components_libs_dir,
        basekV =
            basekV,    
        use_pu_in_PQ =
            use_pu_in_PQ,
        line_data_in_pu =
            line_data_in_pu,

        use_init_u0 =
            use_init_u0,

        use_nlsolve =
            use_nlsolve,

        pf_alg = pf_alg,

        abstol = 1e-12,
        reltol = 1e-12,

        on_fault_time = 5.0,
        clear_fault_time = 5.2,

        list_fault_point_from_node_a = [ 0.3 ],
        list_fault_resistance = [ 0.001 ],
        list_no_line_circuit =  [ 1 ],

        list_edges_to_have_fault = [ 2 ], 
        clear_fault_selection_list = [ 1 ],

        with_faults = false )

generic_network_fault_pertubation_wt_state =
    get_generic_network_fault_pertubation(
        ;case_name = "case9",
        timespan   = 10.0,

        on_fault_time = 5.0,
        clear_fault_time = 5.2,

        list_fault_point_from_node_a = [0.3],
        list_fault_resistance = [0.001],
        list_no_line_circuit =  [1],

        list_edges_to_have_fault = [ 8 ],
        clear_fault_selection_list = [1],

        basekV = 1.0,    
        use_pu_in_PQ = true,
        line_data_in_pu = true,

        with_faults =
            false,
        use_state_in_on_clear_fault =
            true,
        return_extended_results =
            false,

        json_net_data_by_components_file =
            json_net_data_by_components_file,
        components_libs_dir =
            components_libs_dir,
        data_dir =
            data_dir )


generic_network_fault_pertubation =
    get_generic_network_fault_pertubation(
        ;case_name = "case9",
        timespan   = 10.0,

        on_fault_time = 5.0,
        clear_fault_time = 5.1,

        list_fault_point_from_node_a = [0.3],
        list_fault_resistance = [0.001],
        list_no_line_circuit =  [1],

        list_edges_to_have_fault = [ 8 ],
        clear_fault_selection_list = [1],

        basekV = 1.0,    
        use_pu_in_PQ = true,
        line_data_in_pu = true,

        with_faults =
            false,
        use_state_in_on_clear_fault =
            true,
        return_extended_results =
            true,

        json_net_data_by_components_file =
            json_net_data_by_components_file,
        components_libs_dir =
            components_libs_dir,
        data_dir =
            data_dir )



status_Yint_and_Yred_matrices =
    get_net_status_Yint_and_Yred_matrices(
        list_edges_to_have_fault,
        fault_nodes_idx,
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll,    
        generic_system_dynamics_wt_fault_kwd_para )

"""

#---------------------------------------------------
## comments
#---------------------------------------------------

    
    # #--------------------------------------
    # # Integrator Model simulation
    # #--------------------------------------
    
    # system_integrator =
    #     DifferentialEquations.init(
    #         system_prob,            
    #         ode_alg,
    #         # ImplicitMidpoint(),
    #         # dt = dt,            
    #         tstops = [time_final ],
    #         advance_to_tstop = true,
    #         abstol = abstol,
    #         reltol = reltol )
    
    # proposed_dt =
    #     get_proposed_dt(
    #         system_integrator)

    # @show proposed_dt

    # set_proposed_dt!(
    #         system_integrator, dt)
    
    # # if proposed_dt > on_fault_time
        
    # #     set_proposed_dt!(
    # #         system_integrator, dt)
    # # end

    # current_time = system_integrator.t

    # @show current_time
    
    # while current_time < on_fault_time
            
    #     # DifferentialEquations.step!(
    #     #     system_integrator, dt, true)

    #     Δt = on_fault_time - current_time
        
    #     DifferentialEquations.step!(
    #         system_integrator, Δt , true)
        
    #     current_time = system_integrator.t + dt
        
    #     if current_time >= on_fault_time  - dt
            
    #         system_integrator.p .=
    #             fault_model_para
            
    #         add_tstop!(system_integrator,
    #                    clear_fault_time)
            
    #     end

        
    # end
    
    # post_fault1_ntuple_t_δ_ω =
    #     (t = system_integrator.t,
         
    #      δ_in_state =
    #          system_integrator(
    #              system_integrator.t)[ δ_idx_in_state],
         
    #      ω_in_state =
    #          system_integrator(
    #              system_integrator.t)[ ω_idx_in_state],

    #     porder0 =
    #     system_integrator.p[porder0_Idx] )

    # @show post_fault1_ntuple_t_δ_ω

    # system_integrator.p .=
    #     fault_model_para

    # @show fault_model_para
    
    # set_proposed_dt!(
    #     system_integrator,
    #     dt)
    
    # DifferentialEquations.step!(
    #     system_integrator,
    #     dt,
    #     true)


    # post_fault2_ntuple_t_δ_ω =
    #     (t = system_integrator.t,
         
    #      δ_in_state =
    #          system_integrator(
    #              system_integrator.t)[ δ_idx_in_state],
         
    #      ω_in_state =
    #          system_integrator(
    #              system_integrator.t)[ ω_idx_in_state],

    #      porder0 =
    #     system_integrator.p[porder0_Idx] )

    # @show post_fault2_ntuple_t_δ_ω
    
    # # step : fault
    
    # # add_tstop!(system_integrator,
    # #           on_fault_time)

    # # set_proposed_dt!(
    # #     system_integrator, dt)


    # set_proposed_dt!(
    #     system_integrator, dt)
    
    # DifferentialEquations.step!(
    #     system_integrator,
    #     dt,
    #     true)
    
    # post_fault3_ntuple_t_δ_ω =
    #     (t = system_integrator.t,
         
    #      δ_in_state =
    #          system_integrator(
    #              system_integrator.t)[ δ_idx_in_state],
         
    #      ω_in_state =
    #          system_integrator(
    #              system_integrator.t)[ ω_idx_in_state],

    #      porder0 =
    #     system_integrator.p[porder0_Idx] )

    # @show post_fault3_ntuple_t_δ_ω
    
    # post_fault_current_time = system_integrator.t

    # @show post_fault_current_time
    
    # while post_fault_current_time <= clear_fault_time
            
    #     # DifferentialEquations.step!(
    #     #     system_integrator, dt, true)

    #     # Δt = clear_fault_time - post_fault_current_time

    #     set_proposed_dt!(
    #         system_integrator, dt)
        
    #     Δt =  dt
        
    #     DifferentialEquations.step!(
    #         system_integrator, Δt ,
    #         true)
        
    #     post_fault_current_time =
    #         system_integrator.t

    #     future_time_step =
    #         post_fault_current_time + Δt
        
    #     if  future_time_step >=
    #         clear_fault_time - 10 * Δt
            
    #         system_integrator.p .=
    #             clear_fault_model_para
            
            
    #     end

        
    #     if  future_time_step >= clear_fault_time
            
    #         @show future_time_step

    #         future_time_tstop =
    #             clear_fault_time + 30.0 * 1/60.0

    #         @show  future_time_step
            
    #         add_tstop!(
    #             system_integrator,
    #             future_time_tstop )
                        
    #     end
        
    # end

    # post_clear_ntuple_t = system_integrator.t

    # @show post_clear_ntuple_t


    # system_integrator.p .=
    #             clear_fault_model_para
    
    # set_proposed_dt!(
    #     system_integrator, dt)
    
    # DifferentialEquations.step!(
    #     system_integrator,
    #     dt,
    #     true)
    
    # post_clear_ntuple_t_δ_ω =
    #     (t = system_integrator.t,
         
    #      δ_in_state =
    #          system_integrator(
    #              system_integrator.t)[ δ_idx_in_state],
         
    #      ω_in_state =
    #          system_integrator(
    #              system_integrator.t)[ ω_idx_in_state],

    #      porder0 =
    #     system_integrator.p[porder0_Idx] )

    # @show post_clear_ntuple_t_δ_ω

    # set_proposed_dt!(
    #     system_integrator, dt)
    
    # DifferentialEquations.step!(
    #         system_integrator, dt ,
    #         true)
    
    # # system_integrator.p .=
    # #     clear_fault_model_para
    # #     # fault_model_para

    # # set_proposed_dt!(
    # #     system_integrator, dt)
    
    # # DifferentialEquations.step!(
    # #         system_integrator, dt, true) 

    # # post_clear_t = system_integrator.t

    # # @show post_clear_t
    
    # # proposed_dt =
    # #     get_proposed_dt(
    # #         system_integrator)


    # # current_time = system_integrator.t

    # # if proposed_dt > clear_fault_time - current_time
        
    # #     set_proposed_dt!(
    # #         system_integrator, dt)
    # # end

    # # while system_integrator.t < clear_fault_time
            
    # #     DifferentialEquations.step!(
    # #         system_integrator, dt, true)

    # #     if system_integrator.t > clear_fault_time - 0.02

    # #         system_integrator.p .=
    # #             clear_fault_model_para
    # #     end
        

    # # end
 
    # # pre_clear_ntuple_t_δ_ω =
    # #     (t = system_integrator.t,
         
    # #      δ_in_state =
    # #          system_integrator(
    # #              system_integrator.t)[ δ_idx_in_state],
         
    # #      ω_in_state =
    # #          system_integrator(
    # #     system_integrator.t)[ ω_idx_in_state] )

    # # @show pre_clear_ntuple_t_δ_ω
    
    # # # step : clear
    
    # # # add_tstop!(system_integrator,
    # # #            clear_fault_time)
    
    # # system_integrator.p .=
    # #     clear_fault_model_para

    
    # # DifferentialEquations.step!(
    # #         system_integrator, dt, true) 

    
    # # # system_integrator.u .=
    # # #     clear_fault_states

    # # add_tstop!(system_integrator,
    # #            clear_fault_time + 0.5 )
    
    # # proposed_dt =
    # #     get_proposed_dt(
    # #         system_integrator)


    # # current_time = system_integrator.t

    # # if proposed_dt > clear_fault_time + 0.5 - current_time
        
    # #     set_proposed_dt!(
    # #         system_integrator, dt)
    # # end

    # # while system_integrator.t < clear_fault_time + 0.5 
            
    # #     DifferentialEquations.step!(
    # #         system_integrator, dt, true) 

    # # end
    
    # # DifferentialEquations.step!(
    # #     system_integrator)

    # # step : 30 cycles post clear

    
    # # post_clear_30_cyles_ntuple_t_δ_ω =
    # #     (t = system_integrator.t,
         
    # #      δ_in_state =
    # #          system_integrator(
    # #              system_integrator.t)[ δ_idx_in_state],
         
    # #      ω_in_state =
    # #          system_integrator(
    # #     system_integrator.t)[ ω_idx_in_state] )

    # # @show post_clear_30_cyles_ntuple_t_δ_ω

    # system_sol = DifferentialEquations.solve!(
    #     system_integrator)
    

    #--------------------------------------
    # faults callbacks
    #--------------------------------------

    # cb_on_fault = DiscreteCallback(
    #     (u, t, integrator) ->
    #         on_fault_condition(
    #             u, t, integrator,
    #             on_fault_time ),

    #     (integrator) ->
    #         on_fault_wt_model_dynamics_para_affect!(
    #             integrator,
    #             fault_model_dynamics_para ); 
    #     save_positions=(true,true) )

    # cb_clear_fault_no_state = DiscreteCallback(
    #     (u, t, integrator) ->
    #         clear_fault_condition(
    #             u, t, integrator,
    #             clear_fault_time ),

    #     (integrator) ->
    #         clear_fault_wt_model_dynamics_para_affect!(
    #             integrator,
    #             post_fault_model_dynamics_para);
    #     save_positions=(true,true) )

    # cb_clear_fault_wt_state = DiscreteCallback(
    #     (u, t, integrator) ->
    #         clear_fault_condition(
    #             u, t,
    #             integrator,
    #             clear_fault_time ),

    #     (integrator) ->
    #         clear_fault_wt_state_and_model_para_affect!(
    #             integrator,
    #             post_fault_model_dynamics_para,
    #         post_fault_u0_model_states_init);
    #     save_positions=(true,true) )

    #--------------------------------------

    # cb_faults_no_state = CallbackSet(
    #     cb_on_fault,
    #     cb_clear_fault_no_state)


    # cb_faults_wt_state = CallbackSet(
    #     cb_on_fault,
    #     cb_clear_fault_wt_state)

    # fault_callbacks =
    #     cb_fun_make_fault_callbacks(
    #         on_fault_time,
    #         clear_fault_time,
    #         fault_model_dynamics_para,
    #         post_fault_model_dynamics_para,
    #         post_fault_u0_model_states_init;
    #         use_state_in_on_clear_fault =
    #             use_state_in_on_clear_fault )

    # on_fault_model_para =
    #     fault_model_dynamics_para

    # @show on_fault_time
    
    # @show clear_fault_time


    # pre_fault_cb =
    #     DiscreteCallback(
    #         (u, t, integrator) ->
    #             on_pre_fault_condition(
    #                 u, t, integrator ),
    #         (integrator) ->
    #             pre_fault_affect!(
    #                 integrator,
    #                 on_fault_time)
    #         ; 
    #         save_positions=(true, true)
    #     )
    
    # on_fault_cb =
    #     DiscreteCallback(
    #         (u, t, integrator) ->
    #             on_fault_condition(
    #                 u, t, integrator,
    #                 on_fault_time ),
    #         (integrator) ->
    #             on_fault_wt_model_dynamics_para_affect!(
    #                 integrator,
    #                 fault_model_dynamics_para); 
    #         save_positions=(true, true) )

    # if use_state_in_on_clear_fault == true

    #     on_clear_fault_cb =
    #         DiscreteCallback(
    #             (u, t, integrator) ->
    #                 clear_fault_condition(
    #                     u, t, integrator,
    #                     clear_fault_time ),
    #             (integrator) ->
    #                 clear_fault_wt_state_and_model_para_affect!(
    #                     integrator,
    #                     clear_fault_model_para,
    #                     clear_fault_states);
    #             save_positions=(true,true) )

    # else

    #     on_clear_fault_cb =
    #         DiscreteCallback(
    #             (u, t, integrator) ->
    #                 clear_fault_condition(
    #                     u, t, integrator,
    #                     clear_fault_time ),
    #             (integrator) ->
    #                 clear_fault_wt_model_dynamics_para_affect!(
    #                     integrator,
    #                     clear_fault_model_para);
    #             save_positions=(true,true) )
        
    # end


    #  fault_callbacks =
    #      CallbackSet(pre_fault_cb,
    #                  on_fault_cb ,
    #                  on_clear_fault_cb)
    
    # fault_clear_tstops =
    #     [on_fault_time, clear_fault_time]
    
    # system_sol =
    #     DifferentialEquations.solve(
    #         system_prob,
    #         # ImplicitMidpoint(),
    #         # dt = 0.001,
    #         ode_alg,
    #         callback = fault_callbacks,
    #         tstops = [4.0, on_fault_time,  clear_fault_time],
    #         # advance_to_tstop = true,
    #         abstol = abstol,
    #         reltol = reltol )

    
    # system_function = ODEFunction(
    #     (dx,x,p,t) ->
    #         model_dynamics_fun!(
    #             dx, x,
    #             model_dynamics_para,
    #             t;
    #             kwd_para =
    #                 model_dynamics_kwd_para);
    #     syms =
    #         model_syms,
    #         mass_matrix =
    #             model_mass_matrix )

    #  system_sol = ODEProblem(
    #     system_function,
    #     u0_model_states_init,
    #     sim_timespan,
    #      model_dynamics_para )

    
    # system_sol =
    #     DifferentialEquations.solve(
    #        system_prob,
    #         ode_alg,
    #         callback = fault_callbacks,
    #         tstops = fault_clear_tstops,
    #         advance_to_tstop = true,
    #         abstol = abstol,
    #         reltol = reltol )
    
    # plot_gens_δ =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ;sol =
    #             system_sol,
    #         network_vars_labels =
    #             model_syms,
    #         nodes_name = gens_nodes_names,
    #         vars = [:δ ],
    #         tspan = sim_timespan,
    #         fmt = :png)

    # plot_gens_ω =
    #     make_plot_of_buses_vars_and_norm_ω_if_in_vars(
    #         ;sol =
    #             system_sol,
    #         network_vars_labels =
    #             model_syms,
    #         nodes_name = gens_nodes_names,
    #         vars = [:ω ],
    #         tspan = sim_timespan,
    #         fmt = :png)

# (;u0_model_states_init,
#  model_syms,
#  model_mass_matrix,
#  model_bool_dae_vars,

#  ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
#  generic_system_dynamics_wt_fault_kwd_para,

#  gens_nodes_names,
#  SM_gens_nodes_names,
#  non_gens_nodes_names,

#  cb_states,
#  # cb_faults,
#  cb_faults_no_resize,

#  fault_nodes_idx,
#  state_labels,
#  algebraic_vars_labels,

#  cleared_selected_lines_faults_net_para) =
#      NamedTupleTools.select(
#          get_system_simulation_parameters_wt_faults(
#              net_data_by_components_file;
#              components_libs_dir =
#                  components_libs_dir,
#             basekV = basekV,    
#             use_pu_in_PQ = use_pu_in_PQ,
#             line_data_in_pu = line_data_in_pu,
#             use_init_u0 = use_init_u0,
#             use_nlsolve = use_nlsolve,

#             pf_alg = pf_alg,

#             abstol = abstol,
#             reltol = reltol,

#             on_fault_time = on_fault_time,
#             clear_fault_time = clear_fault_time,

#             list_fault_point_from_node_a =
#                 list_fault_point_from_node_a,

#             list_fault_resistance =
#                 list_fault_resistance,

#             list_no_line_circuit =
#                 list_no_line_circuit,

#             list_edges_to_have_fault =
#                 list_edges_to_have_fault,

#             clear_fault_selection_list =
#                 clear_fault_selection_list,
#          with_faults =
#                  with_faults),
#          (:u0_model_states_init,
#           :model_syms,
#           :model_mass_matrix,
#           :model_bool_dae_vars,

#           :ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
#           :generic_system_dynamics_wt_fault_kwd_para,

#           :gens_nodes_names,
#           :SM_gens_nodes_names,
#           :non_gens_nodes_names,

#           :cb_states,
#           # :cb_faults,
#           :cb_faults_no_resize,

#           :fault_nodes_idx,

#           :state_labels,
#           :algebraic_vars_labels,

#           :cleared_selected_lines_faults_net_para))

# #----------------------------------------
# #----------------------------------------
# # DAE system dyamanics simulation
# #----------------------------------------

# cb_on_fault = DiscreteCallback(
#     (u, t, integrator) ->
#         on_fault_condition(
#             u, t, integrator, on_fault_time),

#     (integrator) ->
#         on_fault_affect!(
#             integrator, no_lines_fault ); 
#     save_positions=(true,true),
#     initializealg = ShampineCollocationInit() )

# cb_clear_fault = DiscreteCallback(
#     (u, t, integrator) ->
#         clear_fault_condition(
#             u, t, integrator, clear_fault_time),

#     (integrator) ->
#         clear_fault_affect!(
#             integrator, no_cleared_lines_fault);
#     save_positions=(true,true),
#     initializealg = ShampineCollocationInit())

# #--------------------------------------

# cb_faults = CallbackSet(
#     cb_on_fault, cb_clear_fault)

# #----------------------------------------

# model_dynamics_para =
#     (;ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
#      system_fault_status )          

# model_dynamics_kwd_para =
#     generic_system_dynamics_wt_fault_kwd_para

# system_dynamics_fun! =
#     generic_dynamics_wt_pre_fault_post_by_dae_pf_funcs!

# #----------------------------------------

# (fault_nodes_idx,) =
#     NamedTupleTools.select(
#     cleared_selected_lines_faults_net_para,
#     (:pre_clear_fault_nodes_idx,))

# #----------------------------------------

# number_of_faults =
#     length(list_edges_to_have_fault)

# twice_number_of_faults =
#     2 * number_of_faults

# fault_vh_sym =
#     generate_labels_by_nodes_idxs_and_vars(
#         fault_nodes_idx,
#         [:vh];
#         label_prefix = "bus")

# fault_θh_sym =
#     generate_labels_by_nodes_idxs_and_vars(
#         fault_nodes_idx,
#         [:θh];
#         label_prefix = "bus")

# fault_nodes_sym =
#     [fault_vh_sym;
#      fault_θh_sym ]  

# fault_algeb_init = [
#     ones(number_of_faults);
#     zeros(number_of_faults) ]

# fault_model_bool_dae_vars =
#     DAE_BoolVector(
#         0,
#         twice_number_of_faults )

# #----------------------------------------

# model_bool_dae_vars =
#     [model_bool_dae_vars;
#      fault_model_bool_dae_vars]

# model_syms =
#     [model_syms; fault_nodes_sym]

# u0_model_states_init =
#     [u0_model_states_init;
#      fault_algeb_init ]

# #----------------------------------------

# du0_model_states_init =
#     zeros(length(u0_model_states_init))

# res = similar(u0_model_states_init)

# #----------------------------------------

# faults_and_clear_times =
#     [on_fault_time,
#      clear_fault_time]

# #----------------------------------------
# #----------------------------------------

# system_integrator =
#     DifferentialEquations.init(
#         DAEProblem(
#     DAEFunction(
#     (res, dx, x, p, t) ->
#         system_dynamics_fun!(
#             res, dx, x,
#             model_dynamics_para,
#             t;
#             kwd_para =
#                 model_dynamics_kwd_para);
#     syms =
#         model_syms),
#     du0_model_states_init,
#     u0_model_states_init,
#     sim_timespan,
#     model_dynamics_para,
#     differential_vars =
#         model_bool_dae_vars,
#         callback = cb_states),
#         dae_alg,
#         dt = dt,
#         callback = cb_faults_no_resize,
#         tstops =
#             [ time_final],
#         abstol = abstol,
#         reltol = reltol,
#     advance_to_tstop = true )

            
# # -------------------------------------
# # Step 1 : start
# # -------------------------------------


# add_tstop!(system_integrator,
#            on_fault_time)

# DifferentialEquations.step!(system_integrator)

# system_sol =
#     system_integrator.sol


# sim_timespan = (0 , on_fault_time)

      
# # -------------------------------------
# # Step 2 : 
# # -------------------------------------

# add_tstop!(system_integrator,
#            clear_fault_time)


# DifferentialEquations.step!(system_integrator)

# system_sol =
#     system_integrator.sol


# sim_timespan = (0 , clear_fault_time)

      
# # -------------------------------------
# # Step 3 : 
# # -------------------------------------

# DifferentialEquations.solve!(
#     system_integrator)

# system_sol =
#     system_integrator.sol
