#########################################################

"""
Since this file is meant for co-sim, it may be used
by several simulations hence, there is a need to
pass argument to the script

"""

##########################################################

# Make a result directory for a co-simulator

if length(ARGS) != 0
    
    Co_simulator_number_in_string = ARGS[1]
    
else
    
    Co_simulator_number_in_string = "0"
    
end

script_dir = @__DIR__

results_dir =
    joinpath(script_dir, "results")

if !(isdir(results_dir))
    
    mkpath(results_dir)
    
end


figure_dir =
    joinpath(script_dir, "figure")

if !(isdir(figure_dir))
    
    mkpath(figure_dir)
    
end

############################################################

co_results_dir =
    joinpath(
        results_dir,
        "co-simulation-$(Co_simulator_number_in_string)")

if !(isdir(co_results_dir))
    
    mkpath(co_results_dir)
    
end

co_figure_dir =
    joinpath(
        figure_dir ,
        "co-simulation-$(Co_simulator_number_in_string)")

if !(isdir(co_figure_dir))
    
    mkpath(co_figure_dir)
    
end

#############################################################
