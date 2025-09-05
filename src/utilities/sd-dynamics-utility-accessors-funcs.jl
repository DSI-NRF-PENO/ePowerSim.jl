# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


####################################################
# Accessors Functions
####################################################

#---------------------------------------------------
#---------------------------------------------------
# Definition of lenses or optics
#---------------------------------------------------
#---------------------------------------------------

#---------------------------------------------------
#---------------------------------------------------
# Accessors functions
#---------------------------------------------------
#---------------------------------------------------

#---------------------------------------------------
# plants nodes
#---------------------------------------------------


function get_plant_generators_in_net(comp_data)
    
    lens_plant_generators_in_net=
         @optic _.plant_generators

    return getall(comp_data,
                  lens_plant_generators_in_net)[1] 
end



function get_plant_loads_in_net(comp_data)
    
    lens_plant_loads_in_net=
         @optic _.plant_loads

    return getall(comp_data,
                  lens_plant_loads_in_net)[1] 
end


#---------------------------------------------------
# general components
#---------------------------------------------------


function get_comp_idx_in_net(comp_data)
    
    lens_comp_idx_in_net=
         @optic _.idx

    return getall(comp_data,
                  lens_comp_idx_in_net)[1] 
end


function get_comp_plant_type_in_net(comp_data)

    lens_comp_plant_type_in_net=
        @optic _.plant_type

    return getall(comp_data,
                  lens_comp_plant_type_in_net)[1] 
end


#---------------------------------------------------
# Generation nodes
#---------------------------------------------------


function get_gen_plant_isa_slack_in_net(comp_data)

    lens_gen_plant_isa_slack_in_net=
        @optic _.additional_data.isa_slack

    return getall(comp_data,
                  lens_gen_plant_isa_slack_in_net)[1] 
end


function get_gen_plant_p_order_in_net(comp_data)

    lens_gen_plant_p_order_in_net =
        @optic _.additional_data.p_order

    return getall(comp_data,
                  lens_gen_plant_p_order_in_net)[1] 
end

function get_gen_plant_v_ref_in_net(comp_data)

    lens_gen_plant_v_ref_in_net=
        @optic _.additional_data.v_ref

    return getall(comp_data,
                  lens_gen_plant_v_ref_in_net)[1] 
end

function get_gen_plant_ω_ref_in_net(comp_data)

    lens_gen_plant_ω_ref_in_net =
        @optic _.additional_data.ω_ref

    return getall(comp_data,
                  lens_gen_plant_ω_ref_in_net)[1] 
end


function get_gen_dyn_paras_in_net(comp_data)

    lens_gen_dyn_paras_in_net=
        @optic _.components_data.gen

    return getall(comp_data,
                  lens_gen_dyn_paras_in_net)[1] 
end

function get_gen_avr_paras_in_net(comp_data)

    lens_gen_avr_paras_in_net=
        @optic _.components_data.avr

    return getall(comp_data,
                  lens_gen_avr_paras_in_net)[1] 
end


function get_gen_gov_paras_in_net(comp_data)

    lens_gen_gov_paras_in_net=
        @optic _.components_data.gov

    return getall(comp_data,
                  lens_gen_gov_paras_in_net)[1] 
end


function get_gen_type_in_net(comp_data)

    lens_gen_type_in_net=
        @optic _.components_type.gen

    return getall(comp_data,
                  lens_gen_type_in_net)[1] 
end

function get_avr_type_in_net(comp_data)

    lens_avr_type_in_net=
        @optic _.components_type.avr

    return getall(comp_data,
                  lens_avr_type_in_net)[1] 
end

function get_gov_type_in_net(comp_data)

    lens_gov_type_in_net=
        @optic _.components_type.gov

    return getall(comp_data,
                  lens_gov_type_in_net)[1] 
end

#---------------------------------------------------
# load nodes
#---------------------------------------------------

function get_load_type_in_net(comp_data)

   lens_load_type_in_net=
       @optic _.components_type.load

    return getall(comp_data,
                  lens_load_type_in_net)[1] 
end


function get_load_data_in_net(comp_data)

    lens_load_data_in_net=
        @optic _.components_data.load

    return getall(comp_data,
                  lens_load_data_in_net)[1] 
end


function get_load_P_in_net(comp_data)

    lens_load_P_in_net=
        @optic _.components_data.load.P

    return getall(comp_data,
                  lens_load_P_in_net)[1] 
end


function get_load_Q_in_net(comp_data)

    lens_load_Q_in_net=
        @optic _.components_data.load.Q

    return getall(comp_data,
                  lens_load_Q_in_net)[1] 
end



function get_get_load_PQ_tuple_in_net(comp_data)
    # A lens to get power eld 

    lens_load_P_in_net=
        @optic _.components_data.load.P

    lens_load_Q_in_net=
        @optic _.components_data.load.Q
 
    return (getall(comp_data, lens_load_P_in_net)[1],
            getall(comp_data, lens_load_Q_in_net)[1])
end

#---------------------------------------------------
# transmission nodes
#---------------------------------------------------

function get_transmission_type_in_net(comp_data)

    lens_transmission_type_in_net=
        @optic _.components_type.transmission

    return getall(comp_data,
                  lens_transmission_type_in_net)[1] 
end

function get_transmission_data_in_net(comp_data)

    lens_transmission_data_in_net=
        @optic _.components_data.transmission

    return getall(comp_data,
                  lens_transmission_data_in_net)[1] 
end


function get_transmission_P_in_net(comp_data)

    lens_transmission_P_in_net=
        @optic _.components_data.transmission.P

    return getall(comp_data,
                  lens_transmission_P_in_net)[1] 
end


function get_transmission_Q_in_net(comp_data)

    lens_transmission_Q_in_net=
        @optic _.components_data.transmission.Q

    return getall(comp_data,
                  lens_transmission_Q_in_net)[1] 
end

#---------------------------------------------------
# branches 
#---------------------------------------------------


function get_edge_type_in_net(comp_data)

    lens_edge_type_in_net=
        @optic _.components_type.edge_type

    return getall(comp_data,
                  lens_edge_type_in_net)[1] 
end


#---------------------------------------------------
# selection
#---------------------------------------------------

"""
`get_selected_comp_properties_by_lens`

This is an utility function that return a subset of
properties determined by  `selections` as a vector of
namedtuple.

It first assets that the selection is a subset of properties
of `lensed_comp_data`.

`lensed_comp_data` is the result of application of lens_func
to comp_data.

"""
function get_selected_comp_properties_by_lens(
    lens_func,
    comp_data;
    selections = (:nothing,) )

    lensed_comp_data =
        getall(comp_data, lens_func)[1]

    @assert Set(selections) ⊆ Set(propertynames(
        lensed_comp_data))
            
    return  NamedTupleTools.select(
        lensed_comp_data,
        selections)
        
end



"""
`get_selected_comps_properties_by_lens`

This is an utility function that return a subset of
properties determined by  `selections` as a vector of
namedtuple.


selected_gens_paras =
    get_selected_comps_properties_by_lens(
        get_gen_dyn_paras_in_net,
        plant_generators;
        selections =
            gen_static_para_selections )

selected_gens_paras = [
    (P = 72.3, Q = 27.03, vh = 1.04,
     vmax = 1.1, vmin = 0.9,
     Qmax = 300, Qmin = -300,
     Pmax = 250, Pmin = 10,
     Sn = 390.51248379533274),
    (P = 72.3, Q = 27.03, vh = 1.04,
     vmax = 1.1, vmin = 0.9,
     Qmax = 300, Qmin = -300,
     Pmax = 250, Pmin = 10,
     Sn = 390.51248379533274),
    (P = 72.3, Q = 27.03, vh = 1.04,
     vmax = 1.1, vmin = 0.9,
     Qmax = 300, Qmin = -300,
     Pmax = 250, Pmin = 10,
     Sn = 390.51248379533274)]


It is like using a map on

`get_selected_comp_properties_by_lens`

selected_gens_paras =
    get_selected_comp_properties_by_lens.(
        get_gen_dyn_paras_in_net,
        plant_generators;
        selections =
            gen_static_para_selections )

selected_gens_paras =
    map((y) ->
    get_selected_comp_properties_by_lens(
        get_gen_dyn_paras_in_net,
        y;
        selections =
          gen_static_para_selections),    
    plant_generators )

"""
function get_selected_comps_properties_by_lens(
    lens_func,
    comps_data;
    selections = (:nothing,) )
    
    return  get_selected_comp_properties_by_lens.(
        lens_func,
        comps_data;
        selections =
            selections )        
end



"""
# `get_selected_comp_as_nt_vec_by_lens`

This is an utility function that return a subset of
properties determined by  `selections` as a namedtuple of
vectors.

It first assets that the selection is a subset of properties
of `lensed_comp_data`.

`lensed_comp_data` is the result of application of lens_func
to comp_data.


# get_selected_comp_as_nt_vec_by_lens(
#     get_gen_dyn_paras_in_net,
#     plant_generators;
#     selections =
#         gen_static_para_selections,
#     vec_datatype = Float64 )


selected_comp_as_nt_vec =
    (P = [72.3, 72.3, 72.3],
     Q = [27.03, 27.03, 27.03],
     vh = [1.04, 1.04, 1.04],
     vmax = [1.1, 1.1, 1.1],
     vmin = [0.9, 0.9, 0.9],
     Qmax = [300.0, 300.0, 300.0],
     Qmin = [-300.0, -300.0, -300.0],
     Pmax = [250.0, 250.0, 250.0],
     Pmin = [10.0, 10.0, 10.0],
     Sn = [390.512, 390.512, 390.512])

"""
function get_selected_vec_nt_as_nt_vec_by_lens(
    lens_func,
    comps_data;
    selections = (:nothing,),
    vec_datatype = Float64 )

    dim_selections =
        length( selections )
    
    vec_selected_nt =
        get_selected_comp_properties_by_lens.(
            lens_func,
            comps_data;
            selections =
                selections )
    
    vec_vec = Vector{vec_datatype}[
        [] for a_para in 1:dim_selections ]

    for (idx, a_property) in enumerate(selections)
        for a_namedtuple in vec_selected_nt
            push!(vec_vec[idx],
                  getproperty(a_namedtuple,
                              a_property))
        end
    end
    
    return namedtuple(
        OrderedDict(a_sym => a_value
             for (a_sym, a_value) in
                 zip(selections, vec_vec)))
        
end


#---------------------------------------------------
#---------------------------------------------------

"""
# `get_comp_vector_of_namedtuple_by_lens`

This functions uses a lens function for filtering.


The inputs are a `lens_func` and a collection of namedtuple data

"""
function get_comp_vector_of_nt_by_lens(
    lens_func,
    list_namedtuple_data )

    # if typeof(list_namedtuple_data) ∈ Union{OrderedDict, Dict}

    #     return map(lens_func,
    #             values(list_namedtuple_data))        
    # else

    #     return map(lens_func,
    #             list_namedtuple_data )        
    # end


    return map(lens_func,
            list_namedtuple_data )        
    

end

#---------------------------------------------------
#---------------------------------------------------

"""

obj = (a=1, b=(2, 3, "4"))

# multi-valued optics

lense_multi = @optic₊ (_.a, _.b[2])

lense_multi( obj )

@set lense_multi(obj) = (:x, :y)


obj = (a = (b = (c = 1,),),);


la = @optic _.a
lb = @optic _.b
lc = @optic _.c



"""

