# # https://github.com/JinraeKim/FSimBase.jl

# import HELICS
# const h = HELICS

#---------------------------------------------------

# using Logging

#---------------------------------------------------
# Additional utility functions
#---------------------------------------------------

function get_base_loads(
    network_data,
    sorted_loads_list)
    
    element_variable_key_list = ["pd", "qd"]
    
    loads_pg_base_dict =  OrderedDict{String, Float64}()
    
    loads_qg_base_dict =  OrderedDict{String, Float64}()
    
    loads_base_dict =
        Dict{String, OrderedDict{String, Float64}}(
        "pd" => loads_pg_base_dict,
            "qd" => loads_qg_base_dict )
    
    for a_variable in element_variable_key_list
        
        for a_load in sorted_loads_list
            loads_base_dict[a_variable][a_load] =
                network_data["load"][a_load][a_variable]
        end
    end

    return loads_base_dict 
end


function get_base_gens(
    network_data,
    sorted_gens_list)

    element_variable_key_list = ["pg", "qg"]
    
    gens_pg_base_dict =  OrderedDict{String, Float64}()
    
    gens_qg_base_dict =  OrderedDict{String, Float64}()
    
    gens_base_dict =
        Dict{String, OrderedDict{String, Float64}}(
        "pg" => gens_pg_base_dict,
            "qg" => gens_qg_base_dict )
    
    for a_variable in element_variable_key_list
        
        for a_gen in sorted_gens_list
            gens_base_dict[a_variable][a_gen] =
                network_data["gen"][a_gen][a_variable]
        end
    end
    
    return gens_base_dict 
end


function perturb_loads!(
    network_data,
    loads_base_dict,
    sorted_loads_list)

    d = Normal()
    
    element_variable_key_list = ["pd", "qd"]

    for a_variable in element_variable_key_list
        for a_load in sorted_loads_list
            network_data["load"][a_load][a_variable] =
                loads_base_dict[a_variable][a_load] + rand(d)/10.0
        end        
    end
end

function perturb_gens!(
    network_data,
    gens_base_dict,
    sorted_gens_list)

     d = Normal()
    element_variable_key_list = ["pg", "qg"]
    
    for a_variable in element_variable_key_list
        for a_gen in sorted_gens_list
            network_data["gen"][a_gen][a_variable] =
                gens_base_dict[a_variable][a_gen] + rand(d)
        end        
    end
end

################################################################
################################################################
# Helics functions 
################################################################
################################################################

function destroy_federate(fed)
    # Finalize the Federate this function halts all
    # communication in the Federate and disconnects
    # it from the core
    
    h.helicsFederateFinalize(fed)
    
    # Disconnect and free a broker    
    h.helicsFederateDestroy(fed)
    
    # This runs some cleanup routines and tries to
    # close out any residual
    # thread that haven't been shutdown yet.
    
    h.helicsCleanupLibrary()
    @info "Federate destroyed"
end



function get_subscriptions(fed)
    
    subscriptions_count =
        h.helicsFederateGetInputCount(fed)
    
    subscriptions_ID = Dict()
    
    if subscriptions_count != 0
        
        for sub_ith in 0:subscriptions_count-1
            
            index = sub_ith + 1
            
            subscriptions_ID[index] =
                h.helicsFederateGetInputByIndex(fed, index)
        end
    end
    return subscriptions_ID
end
    
function get_publications(fed)
    
    publications_count =
        h.helicsFederateGetPublicationCount(fed)
    
    publicationss_ID = Dict()
    
    if publications_count != 0
        
        for pub_ith in 0:publications_count-1
            
            index = pub_ith + 1
            
            publications_ID[index] =
                h.helicsFederateGetPublicationByIndex(fed, index)
        end
    end
    
    return publicationss_ID 
end

function get_filters(fed)
    
    filters_count =
        h.helicsFederateGetFilterCount(fed)
    
    filters_ID = Dict()
    
    if filters_count != 0
        
        for filter_ith in 0:filters_count - 1
            
            index = filter_ith + 1
            
            filters_ID[index] =
                h.helicsFederateGetFilterByIndex(fed, index)
        end
        
    end
    
    return filters_ID
    
end

function get_endpoints(fed)
    
    endpoints_count =
        h.helicsFederateGetEndpointCount(fed)
    
    endpoints_ID = Dict()
    
    if endpoints_count != 0
        
        for endpoint_ith in 1:endpoints_count
            
            index = endpoint_ith + 1
            
            endpoints_ID[index] =
                h.helicsFederateGetEndpointByIndex(
                    fed, index)
         end
    end
    
    return endpoints_ID
    
end


function ConvertComplexVectorToDoubleVector(
    vectorInput::Vector{ComplexF64})
    
    vectorLength = length(vectorInput)
    
    doubleVectorInput =
        Vector{Float64}(undef,0)
    for cVal in vectorInput
        
	push!(doubleVectorInput, cVal.re)
        
        push!(doubleVectorInput, cVal.im)
        
    end
    
   return doubleVectorInput
end

function ConvertDoubleVectorToComplexVector(
    vectorInput::Vector{Float64})
    
    complexVector = Vector{ComplexF64}(undef, 0)
    
    for i in 1:length(vectorInput)
        if 2*i <= length(vectorInput)
	    push!(complexVector, vectorInput[2*i - 1] + 1im * vectorInput[2*i])
        end
    end
    return complexVector    
end


function a_results_to_file(
    results_dir, element_name,
    a_result_dict,
    sim_time_array,
    sorted_elements_list)

    # I need to check if the element dict is empty
    key_test = collect(keys(a_result_dict))[1]
    
    if length(collect(keys(a_result_dict[key_test]))) !=0
        
        element_results_dict_key_list =
            keys(a_result_dict)

        for a_variable in element_results_dict_key_list
            
            time_sim_and_element_results_dict =
                OrderedDict{String, Vector{Float64}}()
            
            time_sim_and_element_results_dict["time_sim"] =
                sim_time_array

            for key in sorted_elements_list
                
                time_sim_and_element_results_dict[
                    "$(element_name)_$(key)"] =
                        a_result_dict[a_variable][key]
            end

            element_df =
                DataFrame(time_sim_and_element_results_dict)

            element_filename =
                "$(element_name)-$(a_variable).csv"
            
            result_csv_file =
                joinpath(results_dir,
                         element_filename)
            
            CSV.write(result_csv_file, element_df)
        end
    end
end


function a_results_to_file(
    results_dir,
    filename,
    results_df)

    result_csv_file =
        joinpath(results_dir,
                 filename)
    
    CSV.write(result_csv_file, results_df)
    
end


