# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

"""
, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6
"""

#########################################################
# ------------------------------------------------------
#  Plant parameters accessors functions
# ------------------------------------------------------
#########################################################

function get_plant_control_sig_syms(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
    
    lens_plant_control_sig_syms=@optic _.plant_control_sig_syms
    
    return getall(plant, lens_plant_control_sig_syms[1])
end


function get_plant_control_sig(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus , plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
    
    return plant.plant_control_sig
end


function get_component_params_sym(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load , plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
    
    return plant.param
end


function get_component_params_value(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.param_values
end



function get_component_dict_state_syms(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.dict_state_syms
end


function get_component_cb_state_event_func(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_state_event_func
end


function get_component_cb_state_affect_func(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_state_affect_func
end


function get_component_cb_state_syms(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_state_syms
end


function get_component_cb_state_conditions(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_state_conditions
end



function get_component_cb_state_values(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_state_values
end


function get_component_cb_state_sym2Idx(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_state_sym2Idx
end



function get_component_cb_state_dim(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_state_dim
end


function get_component_cb_dyn_state_event_func(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_dyn_state_event_func
end



function get_component_cb_dyn_state_affect_func(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_dyn_state_affect_func
end


function get_component_cb_dyn_state_syms(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_dyn_state_syms
end



function get_component_cb_dyn_state_conditions(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_dyn_state_conditions
end



function get_component_cb_dyn_state_values(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_dyn_state_values
end


function get_component_cb_dict_dyn_state_syms2sw_Idx(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_dict_dyn_state_syms2sw_Idx
end


function get_component_cb_dyn_state_sw_Idx(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_dyn_state_sw_Idx
end



function get_component_cb_dyn_state_sym2Idx(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_dyn_state_sym2Idx
end



function get_component_cb_dyn_param_state_sw(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load  } )
   
    return plant.cb_dyn_param_state_sw
end


function get_component_cb_dyn_state_dim(plant::Union{ plant_cb_idq, plant_cb_v6, plant_cb_direct, plant_cb_millano, plant_cb,  plant_pss_cb_idq, plant_pss_cb_v6, plant_pss_cb_direct, plant_pss_cb_sauer,  plant_pss_cb, plant_cb_inf, plant_cb_inf_bus, plant_pss_cb_inf, plant_pss_cb_inf_bus, plant_Infinite_bus, plant_Infinite_cb_bus, plant_PQ_Const_P, plant_PQ_Const_I, plant_PQ_Const_Z, plant_PQ_dyn_load, plant_SM_idq, plant_SM_v6, plant_rscad, plant_rscad_idq, plant_rscad_v6, plant_Transmission_t1, plant_Transmission_t2, plant_millano, plant_millano_idq, plant_millano_v6, plant_no_gov_idq, plant_no_gov_v6, plant_no_gov_wt_loc_load_idq, plant_no_gov_wt_loc_load_v6, plant_wt_loc_load_idq, plant_wt_loc_load_v6, plant_millano, plant_no_gov_wt_loc_load, plant_wt_loc_load, plant_no_gov_no_avr_wt_loc_load } )
   
    return plant.cb_dyn_state_dim
end
