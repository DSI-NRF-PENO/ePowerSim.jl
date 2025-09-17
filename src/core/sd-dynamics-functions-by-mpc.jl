# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

####################################################

#---------------------------------------------------

"""
    get_components_properties_by_json(
        plant_generators_data_from_json;
        <keyword arguments>)

Returns a vector of properties of components selected by the variables `sequence_order` and `selection`.

# Arguments
- `sequence_order::Tuple{Symbol}=(:components_data, :gen)`: the selection order of components data.
- `selections::Tuple{Symbol}=(:P, :Q)`: the tuple of parameters data to be selected.


"""
function get_components_properties_by_json(
    plant_generators_data_from_json ;
    sequence_order =
        (:components_data, :gen),
    selections =
        (:P, :Q))

    return namedtuple_nested_selection(
        plant_generators_data_from_json;
        sequence_order =
            sequence_order,
        selections =
            selections )

end

"""
    get_selected_edges_data_by_json(
        edge_data_from_json;
        sequence_order =
            (:components_data, ) ,
        selections =
            (:r, :x, :b, :ratio,
             :angle) )


Returns a vector of properties of edges selected by the variables `sequence_order` and `selection`.

"""
function get_selected_edges_data_by_json(
    edge_data_from_json;
    sequence_order =
        (:components_data, ) ,
    selections =
        (:r, :x, :b, :ratio,
         :angle) )

    # selected_edges_data_ =
    #     get_components_properties_by_json(
    # edge_data_from_json ;
    #     sequence_order =
    #         sequence_order,
    #     selections =
    #         selections )
    
    return namedtuple_nested_selection(
        edge_data_from_json ;
        sequence_order =
            sequence_order,
        selections =
            selections )
    
end


#---------------------------------------------------
#---------------------------------------------------
# components libs and case data related functions
#---------------------------------------------------
#---------------------------------------------------

"""
    get_components_libs_and_case_data(
        case_name;
        <keyword arguments> )


Returns selected static parameters and dynamic types of components from network model csv files.
"""
function get_components_libs_and_case_data(
    case_name;
    case_data_dir  = "",
    components_lib = "",
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"],

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"],

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"],

    mpc_scalar_column_select =
        ["mpc_baseMVA" ],

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"],

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"])

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
    
    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")

       case_data_dir =
           joinpath( data_dir,
                     "converted-data",
                     case_name,)        
    end

    #--------------------------------------
    
    mpc_data_dir  =
        case_data_dir

    dyn_data_dir  =
        joinpath(case_data_dir,
                 "dyn")

    json_data_dir =
        joinpath(case_data_dir,
                 "json")

    #--------------------------------------

    avrs_libs_dir =
        joinpath( components_libs_dir,
                 "avrs")

    govs_libs_dir =
        joinpath( components_libs_dir,
                 "govs")

    pss_libs_dir =
        joinpath( components_libs_dir,
                 "pss")

    gens_libs_dir =
        joinpath( components_libs_dir,
                 "gens")

    loads_libs_dir =
        joinpath( components_libs_dir,
                 "loads")

    lines_libs_dir =
        joinpath( components_libs_dir,
                 "lines")

    non_gen_plants_libs_dir =
        joinpath( components_libs_dir,
                 "non-gen-plants")

    plants_libs_dir =
        joinpath( components_libs_dir,
                 "plants")

    nodes_libs_dir =
        joinpath( components_libs_dir,
                 "nodes")

    #--------------------------------------

    avrs_type_libs_file_json =
        joinpath( avrs_libs_dir,
                 "avrs-libs.json")

    govs_type_libs_file_json =
        joinpath( govs_libs_dir,
                 "govs-libs.json")

    pss_type_libs_file_json =
        joinpath( pss_libs_dir,
                 "pss-libs.json")

    gens_type_libs_file_json =
        joinpath( gens_libs_dir,
                 "gens-libs.json")

    loads_type_libs_file_json =
        joinpath( loads_libs_dir,
                 "loads-libs.json")

    lines_type_libs_file_json =
        joinpath( lines_libs_dir,
                 "lines-libs.json")

    non_gen_plants_type_libs_file_json =
        joinpath( non_gen_plants_libs_dir,
                 "non-gen-plants-libs.json")

    plants_type_libs_file_json =
        joinpath( plants_libs_dir,
                 "plants-libs.json")

    nodes_type_libs_file_json =
        joinpath( nodes_libs_dir,
                 "nodes-libs.json")

    #--------------------------------------

    avrs_parameters_libs_file_json =
        joinpath( avrs_libs_dir,
                 "avrs-parameters-libs.json")

    govs_parameters_libs_file_json =
        joinpath( govs_libs_dir,
                 "govs-parameters-libs.json")

    pss_parameters_libs_file_json =
        joinpath( pss_libs_dir,
                 "pss-parameters-libs.json")

    gens_dyn_parameters_libs_file_json =
        joinpath( gens_libs_dir,
                 "gens-dyn-parameters-libs.json")

    #--------------------------------------

    gov_json_data_file =
        joinpath( json_data_dir,
                 "dict_gov_nt_params.json") 

    avr_json_data_file =
        joinpath( json_data_dir,
                 "dict_avr_nt_params.json") 

    # dict_gen_sym_type_json_data_file =
    #     joinpath(json_data_dir,"dict_gen_sym_type.json")

    gens_nt_dynamic_params_json_data_file =
        joinpath( json_data_dir,
                 "gens_nt_dynamic_params.json") 

    sym_gens_dynamic_params_json_data_file =
        joinpath( json_data_dir,
                 "sym_gens_dynamic_params.json") 

    dict_gens_dyn_nt_params_json_data_file =
        joinpath( json_data_dir,
                 "dict_gens_dyn_nt_params.json") 

    #--------------------------------------
    #--------------------------------------

    mpc_branch_file =
        joinpath( mpc_data_dir,
                 "mpc_branch.csv")

    mpc_gen_file =
        joinpath( mpc_data_dir,
                 "mpc_gen.csv")

    mpc_bus_file =
        joinpath( mpc_data_dir,
                 "mpc_bus.csv")

    mpc_scalar_file =
        joinpath( mpc_data_dir,
                 "mpc_scalar.csv")

    #--------------------------------------

    dyn_gens_file =
        joinpath( dyn_data_dir,
                 "dyn_gen.csv")

    dyn_plants_file =
        joinpath( dyn_data_dir,
                 "dyn_plant.csv")

    #--------------------------------------

    # mpc_branch_column_select =
    #     ["fbus", "tbus", "r", "x", "b",
    #      "ratio", "angle", "status"]

    # mpc_gen_column_select =
    #     ["bus", "Pg", "Qg", "Qmax", "Qmin",
    #      "Vg", "mBase", "status", "Pmax","Pmin"]

    # mpc_bus_column_select =
    #     ["bus_i", "type", "Pd",
    #      "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    # mpc_scalar_column_select =
    #     ["mpc_baseMVA" ]

    # dyn_gens_column_select =
    #     ["bus","sym_gen_type",
    #      "sym_gen_dynamic_para"]

    # dyn_plants_column_select =
    #     ["bus","Plant_type","Gen",
    #      "isa_slack","Gov","Exc"]

    #--------------------------------------

    dyn_gens_data_types =
        [Int, Symbol, Symbol]

    dyn_plants_data_types =
        [ Int, Symbol, Symbol, Bool, Symbol, Symbol]

    #--------------------------------------

    mpc_branch_selected_data =
        CSV.File( mpc_branch_file;
                 select=mpc_branch_column_select )

    mpc_gen_selected_data =
        CSV.File( mpc_gen_file;
                 select=mpc_gen_column_select )

    mpc_bus_selected_data =
        CSV.File( mpc_bus_file;
                 select=mpc_bus_column_select )

    mpc_scalar_selected_data =
        CSV.File( mpc_scalar_file;
                 select=mpc_scalar_column_select )
    
    #--------------------------------------

    dyn_gens =
        CSV.File(
            dyn_gens_file;
            select = dyn_gens_column_select,
            types  = dyn_gens_data_types )

    dyn_plants =
        CSV.File(
            dyn_plants_file;
            select = dyn_plants_column_select,
            types  = dyn_plants_data_types )

    #--------------------------------------

    mpc_baseMVA =
        mpc_scalar_selected_data.mpc_baseMVA[1]
    
    return (;mpc_branch_selected_data,
            mpc_gen_selected_data,
            mpc_bus_selected_data,
            dyn_gens,
            dyn_plants,
            mpc_baseMVA)
end

#---------------------------------------------------
#---------------------------------------------------
# csv and xlsx network static data related functions
#---------------------------------------------------
#---------------------------------------------------

"""
    get_case_data_by_csv(
        case_name;
        <keyword arguments> )


Returns selected static parameters and dynamic types of components from network model csv files.
"""
function get_case_data_by_csv(
    case_name ;
    case_data_dir = "",
    mpc_branch_column_select = "",
    mpc_gen_column_select = "",
    mpc_bus_column_select = "",
    mpc_scalar_column_select = "",
    dyn_gens_column_select = "",
    dyn_plants_column_select = "",
    dyn_gens_data_types = "",
    dyn_plants_data_types = "",
    wt_plants_data_types_bool =
        true )

    #--------------------------------------    

    if case_data_dir == ""
        
        package_dir = pkgdir(ePowerSim)

        case_data_dir =
            joinpath(package_dir,
                     "data",
                     "converted-data",
                     case_name )
        
    end

    #--------------------------------------
    
    dyn_data_dir  =
        joinpath(case_data_dir,
                 "dyn")

    # added
    if !(isdir( dyn_data_dir ))

        mkpath( dyn_data_dir )

    end

    #--------------------------------------

    csv_branch_file =
        joinpath(case_data_dir,
                 "mpc_branch.csv")

    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")

    csv_gencost_file =
        joinpath(case_data_dir,
                 "mpc_gencost.csv")
    
    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")

    csv_scalar_file =
        joinpath(case_data_dir,
                 "mpc_scalar.csv")

    #--------------------------------------

    csv_load_type_file  = 
        joinpath(case_data_dir,
                 "mpc_load_type.csv")

    # added

    if !( isfile( csv_load_type_file ) )
        
        create_a_default_case_mpc_load_type(
            case_name;
            
            data_dir      = data_dir,
            
            case_data_dir = case_data_dir )
    end

    # added
    
    csv_branch_type_file =
        joinpath(case_data_dir,
                 "mpc_branch_type.csv")


    if !( isfile( csv_branch_type_file ) )
        
        create_a_default_case_mpc_branch_type(
            case_name;
            
            data_dir      = data_dir,
            
            case_data_dir = case_data_dir )
    end
    
    
    #--------------------------------------

    dyn_gens_file =
        joinpath(dyn_data_dir,
                 "dyn_gen.csv")

    # added


    if !( isfile( dyn_gens_file ) )
        
        create_a_default_case_dyn_gens_file(
            case_name;
            data_dir,
            case_data_dir,
            dyn_gens_file)
    end
    
    dyn_plants_file =
        joinpath(dyn_data_dir,
                 "dyn_plant.csv")


    # added
    
    if !( isfile( dyn_plants_file ) )
        
        create_a_default_case_dyn_plants_file(
            case_name;
            data_dir,
            case_data_dir,

            dyn_plants_file )
        
    end
    
    #--------------------------------------

    if mpc_branch_column_select == ""
        
        mpc_branch_column_select =
            ["fbus", "tbus", "r", "x", "b",
             "ratio", "angle", "status"]
        
    end    

    
    if mpc_gen_column_select == ""
        
        mpc_gen_column_select =
            ["bus", "Pg", "Qg", "Qmax", "Qmin",
             "Vg", "mBase", "status", "Pmax","Pmin"]
    end
    

    if mpc_bus_column_select == ""

         mpc_bus_column_select =
             ["bus_i", "type", "Pd",
              "Qd", "Gs", "Bs", "Vmax", "Vmin"]
    end
    

    if mpc_scalar_column_select == ""
        
        mpc_scalar_column_select =
            ["mpc_baseMVA" ]
        
    end
    

    if dyn_gens_column_select == ""
        
        dyn_gens_column_select =
            ["bus","sym_gen_type",
             "sym_gen_dynamic_para"]
    end
    

    if dyn_plants_column_select == ""
        
        dyn_plants_column_select =
            ["bus","Plant_type","Gen",
             "isa_slack","Gov","Exc"]
        
    end
    
    #--------------------------------------

    
    if dyn_gens_data_types == "" &&
        wt_plants_data_types_bool == true
    
        dyn_gens_data_types =
            [Int, Symbol, Symbol]
        
    elseif dyn_gens_data_types != "" &&
        wt_plants_data_types_bool == true
    
        dyn_gens_data_types =
            [Int, Symbol, Symbol]

    elseif dyn_gens_data_types != "" &&
        wt_plants_data_types_bool == false
    
        dyn_gens_data_types =
            [String, Symbol, Symbol]
    else
            
        dyn_gens_data_types =
            [String, Symbol, Symbol]        
    end

    
    if dyn_plants_data_types == "" &&
        wt_plants_data_types_bool == true
    
        dyn_plants_data_types =
            [Int, Symbol, Symbol, Bool, Symbol, Symbol]
    
    elseif dyn_plants_data_types != "" &&
        wt_plants_data_types_bool == true
        
        dyn_plants_data_types =
            [Int, Symbol, Symbol, Bool, Symbol, Symbol]

    elseif dyn_plants_data_types != "" &&
        wt_plants_data_types_bool == false

        dyn_plants_data_types =
            [String, Symbol, Symbol, Bool, Symbol, Symbol]

    else

        dyn_plants_data_types =
            [String, Symbol, Symbol, Bool, Symbol, Symbol]
    end


    #--------------------------------------
    #--------------------------------------

    mpc_branch_selected_data =
        CSV.File(csv_branch_file;
                 select=mpc_branch_column_select )

    mpc_gen_selected_data =
        CSV.File(csv_gen_file;
                 select=mpc_gen_column_select )

    mpc_gencost_data =
        CSV.File(csv_gencost_file)
    
    mpc_bus_selected_data =
        CSV.File(csv_bus_file;
                 select=mpc_bus_column_select )

    mpc_scalar_selected_data =
        CSV.File(csv_scalar_file;
                 select=mpc_scalar_column_select )
    
    #--------------------------------------

    mpc_load_type_data =
        CSV.File(csv_load_type_file )


    mpc_branch_type_data =
        CSV.File(csv_branch_type_file )

    
    #--------------------------------------

    if wt_plants_data_types_bool == true

        dyn_gens =
            CSV.File(
                dyn_gens_file;
                select = dyn_gens_column_select,
                types  = dyn_gens_data_types )

        dyn_plants =
            CSV.File(
                dyn_plants_file;
                select = dyn_plants_column_select,
                types  = dyn_plants_data_types )
        
    else

        dyn_gens =
            CSV.File(
                dyn_gens_file;
                select = dyn_gens_column_select,
                types  = dyn_plants_data_types)

        dyn_plants =
            CSV.File(
                dyn_plants_file;
                select = dyn_plants_column_select,
                types  = dyn_plants_data_types )
        
    end

    #--------------------------------------

    mpc_baseMVA =
        mpc_scalar_selected_data.mpc_baseMVA[1]
    
    return (;mpc_branch_selected_data,
            mpc_gen_selected_data,
            mpc_gencost_data,
            mpc_bus_selected_data,
            dyn_gens,
            dyn_plants,
            mpc_baseMVA,

            mpc_load_type_data,
            mpc_branch_type_data )
end


"""
    get_net_static_data_by_components_by_xlsx(
        ; case_name
        data_dir,
        by_components )


Returns selected static parameters and dynamic types of components from network model xlsx file.
"""
function get_net_static_data_by_components_by_xlsx(
    ;case_name = "case9",        
    data_dir = "",
    by_components = true )

    #--------------------------------------

    if data_dir == ""

        package_dir = pkgdir(ePowerSim)
        
        data_dir = joinpath(package_dir,
                     "data" )
        
    end

    #--------------------------------------

    case_data_dir =
        joinpath(data_dir,
                 "converted-data",
                 case_name )
    
    mpc_data_dir  = case_data_dir

    #--------------------------------------

    xlsx_file =
        joinpath(case_data_dir, "xlsx",
                 "net-static-data.xlsx")

    #--------------------------------------
    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]
    
    #--------------------------------------
    # xlsx    
    #--------------------------------------

    mpc_branch_type_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_branch_type";
                       header = true,
                       infer_eltypes = true))


    mpc_gencost_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_gencost";
                       header = true,
                       infer_eltypes = true))
    
    mpc_branch_selected_data = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_branch";
                       header = true,
                       infer_eltypes = true)),
        mpc_branch_column_select)

    mpc_gen_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_gen";
                       header = true,
                       infer_eltypes = true)),
        mpc_gen_column_select)

    mpc_bus_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_bus";
                       header = true,
                       infer_eltypes = true)),
        mpc_bus_column_select)

    mpc_scalar_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_scalar";
                       header = true,
                       infer_eltypes = true)),
        mpc_scalar_column_select)

    mpc_baseMVA = mpc_scalar_selected_data.mpc_baseMVA[1]

    #--------------------------------------
    
    return  (;mpc_branch_type_data,             
             mpc_branch_selected_data,
             mpc_gen_selected_data,
             mpc_bus_selected_data,
             mpc_baseMVA,
             mpc_gencost_data )
    
end


#---------------------------------------------------
#---------------------------------------------------
# Network branches and nodes  functions
#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
# generic functions 
#---------------------------------------------------

"""
    get_edge_y_data_by_generic(
        fbus, tbus,
        r, x, b,
        ratio, angle,
        status,
        baseMVA, basekV)

Returns per unit π parameters for branches admittances.
"""
function get_edge_y_data_by_generic(
    fbus, tbus,
    r, x, b,
    ratio, angle,
    status,
    baseMVA, basekV)
    
    # -------------------------------------

    baseZ = basekV^2/baseMVA

    baseY = 1/baseZ
    
    # -------------------------------------    

    ys =
        baseY ./ ( r + im * x ) 
       
    y_c = 1 / 2 * ( im *  b ) * baseY

    inv_τ = (ratio == 0.0 || ratio == 0) ?
        1.0 : 1/ratio

    θ_shift = angle

    return (;from = "bus$(fbus)",
            to = "bus$(tbus)",
            y = ys,
            y_shunt_km = y_c,
            y_shunt_mk = y_c,
            t_ratio = inv_τ)
    
end

"""
    get_edge_y_line_data_by_generic(
        fbus, tbus,
        r, x, b,
        ratio, angle,
        status,
        baseMVA, basekV)

Returns per unit π parameters for branches admittances.
"""
function get_edge_y_line_data_by_generic(
    fbus, tbus,
    r, x, b,
    ratio, angle,
    status; baseMVA=1.0, basekV=1.0 )
    
    # -------------------------------------

    baseZ = basekV^2/baseMVA

    baseY = 1/baseZ
    
    # -------------------------------------

    ys = baseY ./ (r + im * x) 
       
    y_c =  1 / 2 * b * baseY

    inv_τ = (ratio == 0.0 || ratio == 0) ?
        1.0 : 1/ratio

    θ_shift = angle

    return ( from = "bus$(fbus)",
             to = "bus$(tbus)",
             y = ys,
             y_shunt_km = y_c,
             y_shunt_mk = y_c )
    
end


"""
    get_edge_y_transformer_data_by_generic(
        fbus, tbus,
        r, x, b,
        ratio, angle,
        status,
        baseMVA, basekV)

Returns per unit π parameters for transformers admittances.
"""
function get_edge_y_transformer_data_by_generic(
    fbus, tbus,
    r, x, b,
    ratio, angle,
    status; baseMVA=1.0, basekV=1.0 )

    # -------------------------------------

    baseZ = basekV^2/baseMVA

    baseY = 1/baseZ
    
    # -------------------------------------
    
    ys = baseY ./ (r + im * x) 
       
    y_c = 1 / 2 * b * baseY

    inv_τ = (ratio == 0.0 || ratio == 0) ?
        1.0 : 1/ratio

    θ_shift = angle

    return (from = "bus$(fbus)",
            to = "bus$(tbus)",
            y = ys,
            y_shunt_km = y_c,
            y_shunt_mk = y_c,
            t_ratio = inv_τ )
    
end

"""
    get_edges_orientation_by_generic(
        branches_fbus,
        branches_tbus )


Returns branches orientations in form of a list of tuples `(src, dst)`.
"""
function  get_edges_orientation_by_generic(
    branches_fbus,
    branches_tbus )

    # edges_orientation
    
    return [(fbus, tbus)
            for (fbus, tbus) in
                zip( branches_fbus, branches_tbus)]
end


"""
    get_edges_Ybr_by_generic(
        r, x, b,
        ratio, angle,
        edge_type,
        Gs, Bs;
        baseMVA = 1.0,
        basekV = 1.0 )


Returns a list of per unit π admittance matrices for branches.
"""
function get_edges_Ybr_by_generic(
    r, x, b,
    ratio, angle,
    edge_type,
    Gs, Bs;
    baseMVA = 1.0,
    basekV = 1.0 )
    
    # -------------------------------------

    baseZ = basekV^2/baseMVA

    baseY = 1/baseZ
    
    # -------------------------------------

    ys = baseY ./ (r + im * x) 
       

    y_c = 1 / 2 * (im *  b) *  baseY

    # y_sh =  (Gs .+ im * Bs)/baseMVA

    inv_τ =
        [ ((a_ratio == 0.0 || a_ratio == 0) &&
        a_edge_type != "Transformer") ? 1.0 :
        ((a_ratio == 0.0 || a_ratio == 0) &&
        a_edge_type == "Transformer") ? 1.0 :
        1/a_ratio
          for (a_ratio, a_edge_type) in
              zip(ratio, edge_type)  ]
    
    θ_shift = angle

    #---------------------------------------------------
    
    Yff = (ys +  y_c) .* (inv_τ).^2

    Ytf = -ys .* (inv_τ ./ exp.(im * θ_shift))

    Yft = -ys .* (inv_τ ./ exp.(-im * θ_shift))

    Ytt = (ys +  y_c)

    #---------------------------------------------------
    
    # edges_Ybr_cal
    
    return [ [yff ytf; yft ytt ]
              for (yff, ytf, yft, ytt) in
                  zip( Yff, Ytf, Yft, Ytt ) ]
    

end


"""
    get_nodes_Yshunt_by_generic(
        buses_Gs,
        buses_Bs;
        baseMVA = 1.0 )


Returns a list of per unit nodes shunt admittance.
"""
function get_nodes_Yshunt_by_generic(
    buses_Gs,
    buses_Bs;
    baseMVA = 1.0 )
             
    # y_sh =
    #     (bus_Gs .+ im *  bus_Bs) ./ mpc_baseMVA
    
    return (buses_Gs .+ im *  buses_Bs) ./ baseMVA

end


"""
    get_nodes_idx_and_Yshunt_non_zero(
        buses_Gs,
        buses_Bs;
        baseMVA = 1.0 )


Returns namedtuples `y_sh_shunt_exist`, `nz_y_sh_idxs`, `nz_y_sh`  for nodes non-zero shunt admittances.
"""
function get_nodes_idx_and_Yshunt_non_zero(
    buses_Gs,
    buses_Bs;
    baseMVA = 1.0 )
    
    y_sh = get_nodes_Yshunt_by_generic(
        buses_Gs, buses_Bs;
        baseMVA =
            baseMVA )

    nz_y_sh_idxs =
        findall(x -> (real(x) != 0 || real(x) != 0.0) || (
            imag(x) != 0 || imag(x) != 0.0), y_sh )

    y_sh_shunt_exist =
        length(nz_y_sh_idxs) != 0 ? true : false

   
    return (;y_sh_shunt_exist,
            nz_y_sh_idxs,
            nz_y_sh = y_sh[nz_y_sh_idxs] )

end


"""
    get_nodes_idx_and_Yshunt_by_generic(
        buses_idx,
        buses_Gs,
        buses_Bs;
        baseMVA = 1.0 )


Returns a list of tuples of nodes indices and nodes shunt admittances.
"""
function get_nodes_idx_and_Yshunt_by_generic(
    buses_idx,
    buses_Gs,
    buses_Bs;
    baseMVA = 1.0 )
             
    y_sh = get_nodes_Yshunt_by_generic(
        buses_Gs, buses_Bs;
        baseMVA =
            baseMVA )
    
    return [(node_idx, y_shunt)
            for (node_idx, y_shunt) in
                zip(buses_idx, y_sh)]

end


"""
    get_edges_Ybr_cal_and_edges_orientation_by_generic(
        branches_fbus,
        branches_tbus,
        r, x, b,
        ratio, angle,
        edge_type,
        Gs, Bs;
        baseMVA = 1.0, basekV = 1.0 )

Returns namedtuple of branches π admittance matrices and branches orientations.
"""
function get_edges_Ybr_cal_and_edges_orientation_by_generic(
    branches_fbus,
    branches_tbus,
    r, x, b,
    ratio, angle,
    edge_type,
    Gs, Bs;
    baseMVA = 1.0, basekV = 1.0 )

    
    edges_Ybr_cal =
        get_edges_Ybr_by_generic(
            r, x, b,
            ratio, angle,
            edge_type, Gs, Bs;
            baseMVA = baseMVA,
            basekV =  basekV )

    edges_orientation =
        get_edges_orientation_by_generic(
            branches_fbus, branches_tbus )
    
    return (;edges_Ybr_cal,
            edges_orientation )
    
end
    
"""
    get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data(
        edges_r, edges_x, edges_b,
        edges_ratio, edges_angle,
        Gs, Bs;
        edges_fbus,
        edges_tbus,
        edges_type,    
        all_nodes_idx,
        n2s_all_nodes_idx,
        baseMVA=1.0,
        basekV=1.0,
        line_data_in_pu = true)


Returns a namedtuple of network admitance vectors and network nodes neigbouhood vectors.
"""
function get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data(
    edges_r, edges_x, edges_b,
    edges_ratio, edges_angle,
    Gs, Bs;
    edges_fbus,
    edges_tbus,
    edges_type,    
    all_nodes_idx,
    n2s_all_nodes_idx,
    baseMVA=1.0,
    basekV=1.0,
    line_data_in_pu = true)

    if line_data_in_pu == true

        return get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
            edges_fbus,
            edges_tbus,
            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle,
            Gs,
            Bs;
            edges_type,
            all_nodes_idx,
            n2s_all_nodes_idx,
            baseMVA = 1.0,
            basekV = 1.0,
            baseShunt = baseMVA )        
    else

        return get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
            edges_fbus,
            edges_tbus,
            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle,
            Gs,
            Bs;
            edges_type,
            all_nodes_idx,
            n2s_all_nodes_idx,
            baseMVA = baseMVA,
            basekV = basekV,
           baseShunt = baseMVA )        
    end

end

# """
# See [`get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data`](@ref)

# """

# @doc (@doc get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data)
"""
    get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data(
        ;edges_fbus, edges_tbus,
        edges_type,
        edges_r, edges_x, edges_b,
        edges_ratio, edges_angle,
        Gs, Bs,
        all_nodes_idx,
        n2s_all_nodes_idx,
        baseMVA=1.0,
        basekV=1.0,
        line_data_in_pu = true)


Returns a namedtuple of network admitance vectors and network nodes neigbouhood vectors.
"""
function get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data(
    ;edges_fbus, edges_tbus,
    edges_type,
    edges_r, edges_x, edges_b,
    edges_ratio, edges_angle,
    Gs, Bs,
    all_nodes_idx,
    n2s_all_nodes_idx,
    baseMVA=1.0,
    basekV=1.0,
    line_data_in_pu = true)

    if line_data_in_pu == true

        return get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
            edges_fbus,
            edges_tbus,
            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle,
            Gs,
            Bs;
            edges_type,
            all_nodes_idx,
            n2s_all_nodes_idx,
            baseMVA = 1.0,
            basekV = 1.0,
            baseShunt = baseMVA )        
    else

        return get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
            edges_fbus,
            edges_tbus,
            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle,
            Gs,
            Bs;
            edges_type,
            all_nodes_idx,
            n2s_all_nodes_idx,
            baseMVA = baseMVA,
            basekV = basekV,
           baseShunt = baseMVA )        
    end

end


# """
# See [`get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data`](@ref)

# """

# @doc (@doc get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data)
"""

    get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data(
        ;edges_fbus, edges_tbus,
        edges_type,
        edges_r, edges_x, edges_b,
        edges_ratio, edges_angle,
        Gs, Bs,
        all_nodes_idx,
        n2s_all_nodes_idx,
        baseMVA=1.0,
        basekV=1.0,
        line_data_in_pu = true)


Returns a namedtuple of network admitance vectors and network nodes neigbouhood vectors.
"""
function get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data(
    ;edges_fbus, edges_tbus,
    edges_type,
    edges_r, edges_x, edges_b,
    edges_ratio, edges_angle,
    Gs, Bs,
    all_nodes_idx,
    n2s_all_nodes_idx,
    baseMVA=1.0,
    basekV=1.0,
    line_data_in_pu = true)

    if line_data_in_pu == true

        return get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
            edges_fbus,
            edges_tbus,
            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle,
            Gs,
            Bs;
            edges_type,
            all_nodes_idx,
            n2s_all_nodes_idx,
            baseMVA = 1.0,
            basekV = 1.0,
            baseShunt = baseMVA )        
    else

        return get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
            edges_fbus,
            edges_tbus,
            edges_r,
            edges_x,
            edges_b,
            edges_ratio,
            edges_angle,
            Gs,
            Bs;
            edges_type,
            all_nodes_idx,
            n2s_all_nodes_idx,
            baseMVA = baseMVA,
            basekV = basekV,
           baseShunt = baseMVA )        
    end

end


# @doc (@doc get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data)
"""
    get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
        branches_fbus,
        branches_tbus,
        r,
        x,
        b,
        ratio,
        angle,        
        Gs,
        Bs;
        edges_type,
        all_nodes_idx,
        n2s_all_nodes_idx,
        baseMVA = 1.0,
        basekV = 1.0,
        baseShunt = 1.0 )


Returns a namedtuple of network admitance vectors and network nodes neigbouhood vectors.
"""
function get_sense_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
    branches_fbus,
    branches_tbus,
    r,
    x,
    b,
    ratio,
    angle,        
    Gs,
    Bs;
    edges_type,
    all_nodes_idx,
    n2s_all_nodes_idx,
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0 )

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt_by_generic(
            transformed_all_nodes_idx,
            Gs, Bs;
            baseMVA =
                baseShunt )
    
    edges_orientation =
        get_edges_orientation_by_generic(
            branches_fbus,
            branches_tbus )

    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )
    
    edges_Ybr_cal =
        get_edges_Ybr_by_generic(
            r, x, b,
            ratio, angle,
            edges_type,
            Gs, Bs;
            baseMVA = baseMVA, basekV = basekV )

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
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx]]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(transformed_all_nodes_idx,
                nodes_incident_edges_and_orientation)]
    
    #------------------------------------------

    nodes_incident_edges_orientation_and_Ybr_cal =[
        edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges ]

    Ynet_no_shunt =
        [ [ k == 0 ?
        sum([ node_idx == first( first(
            orient_and_Ybr) ) ?
                last(orient_and_Ybr)[1] :
                last(orient_and_Ybr)[4]
              for orient_and_Ybr in
                  orientations_and_edges_Ybr]) :
                      node_idx == first(first(
                          orientations_and_edges_Ybr[k] )) ?
                             last(orientations_and_edges_Ybr[
                                  k])[3] :
                     last(orientations_and_edges_Ybr[k])[2]
            for k in 0:length(
                orientations_and_edges_Ybr ) ]
          for (node_idx, orientations_and_edges_Ybr) in 
            zip(transformed_all_nodes_idx,
             nodes_incident_edges_orientation_and_Ybr_cal)]


    Ynet =
        [
          [  idx == 1 ?
                last(node_k_idx_and_shunt) +
                Ynet_k_element :
                Ynet_k_element
               for (idx, Ynet_k_element) in
                   enumerate(
                       Ynet_no_shunt[shunt_idx] ) ]
            for (shunt_idx, node_k_idx_and_shunt) in
                enumerate( nodes_idx_and_Yshunt) ]

    return (; 
            Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
end


# @doc (@doc get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data)
"""
    get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
        branches_fbus,
        branches_tbus,
        r,
        x,
        b,
        ratio,
        angle,        
        Gs,
        Bs;
        edges_type,
        all_nodes_idx,
        n2s_all_nodes_idx,
        baseMVA = 1.0,
        basekV = 1.0,
        baseShunt = 1.0 )


Returns a namedtuple of network admitance vectors and network nodes neigbouhood vectors.
"""
function get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
    branches_fbus,
    branches_tbus,
    r,
    x,
    b,
    ratio,
    angle,        
    Gs,
    Bs;
    edges_type,
    all_nodes_idx,
    n2s_all_nodes_idx,
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0 )

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt_by_generic(
            transformed_all_nodes_idx,
            Gs, Bs;
            baseMVA =
                baseShunt )
    
    edges_orientation =
        get_edges_orientation_by_generic(
            branches_fbus,
            branches_tbus )

    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )
    
    edges_Ybr_cal =
        get_edges_Ybr_by_generic(
            r, x, b,
            ratio, angle,
            edges_type,
            Gs, Bs;
            baseMVA = baseMVA, basekV = basekV )

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
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx]]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(transformed_all_nodes_idx,
                nodes_incident_edges_and_orientation)]
    
    #------------------------------------------

    nodes_incident_edges_orientation_and_Ybr_cal =[
        edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges ]

    Ynet_no_shunt =
        [ [ k == 0 ?
        sum([ node_idx == first( first(
            orient_and_Ybr) ) ?
                last(orient_and_Ybr)[1] :
                last(orient_and_Ybr)[4]
              for orient_and_Ybr in
                  orientations_and_edges_Ybr]) :
                      node_idx == first(first(
                          orientations_and_edges_Ybr[k] )) ?
                             last(orientations_and_edges_Ybr[
                                  k])[3] :
                     last(orientations_and_edges_Ybr[k])[2]
            for k in 0:length(
                orientations_and_edges_Ybr ) ]
          for (node_idx, orientations_and_edges_Ybr) in 
            zip(transformed_all_nodes_idx,
             nodes_incident_edges_orientation_and_Ybr_cal)]


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

    return (; 
            Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
end


# @doc (@doc get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data)
"""
    get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
        branches_fbus,
        branches_tbus,
        r,
        x,
        b,
        ratio,
        angle,
        edge_type,    
        Gs,
        Bs;
        all_nodes_idx,
        n2s_all_nodes_idx,
        baseMVA = 1.0,
        basekV = 1.0,
        baseShunt = 1.0 )


Returns a namedtuple of network admitance vectors and network nodes neigbouhood vectors.
"""
function get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
    branches_fbus,
    branches_tbus,
    r,
    x,
    b,
    ratio,
    angle,
    edge_type,    
    Gs,
    Bs;
    all_nodes_idx,
    n2s_all_nodes_idx,
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0 )

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt_by_generic(
            transformed_all_nodes_idx,
            Gs, Bs;
            baseMVA =
                baseShunt )
    
    edges_orientation =
        get_edges_orientation_by_generic(
            branches_fbus,
            branches_tbus )

    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )
    
    edges_Ybr_cal =
        get_edges_Ybr_by_generic(
            r, x, b,
            ratio, angle,
            edge_type,
            Gs, Bs;
            baseMVA = baseMVA, basekV = basekV )

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
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx]]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(transformed_all_nodes_idx,
                nodes_incident_edges_and_orientation)]
    
    #------------------------------------------

    nodes_incident_edges_orientation_and_Ybr_cal =[
        edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges ]

    Ynet_no_shunt =
        [ [ k == 0 ?
        sum([ node_idx == first( first(
            orient_and_Ybr) ) ?
                last(orient_and_Ybr)[1] :
                last(orient_and_Ybr)[4]
              for orient_and_Ybr in
                  orientations_and_edges_Ybr]) :
                      node_idx == first(first(
                          orientations_and_edges_Ybr[k] )) ?
                             last(orientations_and_edges_Ybr[
                                  k])[3] :
                     last(orientations_and_edges_Ybr[k])[2]
            for k in 0:length(
                orientations_and_edges_Ybr ) ]
          for (node_idx, orientations_and_edges_Ybr) in 
            zip(transformed_all_nodes_idx,
             nodes_incident_edges_orientation_and_Ybr_cal)]


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

    return (; 
            Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
end

# @doc (@doc get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data)
"""
    get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
        branches_fbus,
        branches_tbus,
        r, x, b,
        ratio, angle,
        edge_type,
        buses_idx,
        Gs, Bs;
        baseMVA = 1.0,
        basekV = 1.0,
        baseShunt = 1.0 )


Returns a namedtuple of network admitance vectors and network nodes neigbouhood vectors.
"""
function get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
    branches_fbus,
    branches_tbus,
    r, x, b,
    ratio, angle,
    edge_type,
    buses_idx,
    Gs, Bs;
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0 )

    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt_by_generic(
            buses_idx,
            Gs, Bs;
            baseMVA =
                baseShunt )
    
    edges_orientation =
        get_edges_orientation_by_generic(
            branches_fbus,
            branches_tbus )

    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )
    
    edges_Ybr_cal =
        get_edges_Ybr_by_generic(
            r, x, b,
            ratio, angle,
            edge_type,
            Gs, Bs;
            baseMVA = baseMVA, basekV = basekV )

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

    
    # nodes_idx_with_adjacent_nodes_idx = [
    #     [ [[idx],
    #        [ idx == orient[1] ?
    #            orient[2] : orient[1]
    #          for orient in
    #              a_node_edges_other_nodes_idx]]...; ]
    #     for (idx, a_node_edges_other_nodes_idx) in
    #         enumerate(
    #             nodes_incident_edges_and_orientation)]


    # this address non continous nodes labeling
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx]]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(buses_idx,
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
        for (node_idx, orientations_and_edges_Ybr) in 
            zip(buses_idx,
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

    return (; 
            Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
end


#--------------------------------------

# @doc (@doc get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data)
"""
    get_Ynet(
        edge_data_from_json,
        shunt_data_from_json;
        baseMVA = 1.0,
        basekV = 1.0,
        baseShunt = 1.0,
        line_data_in_pu = true)


Returns a namedtuple of network admitance vectors and network nodes neigbouhood vectors.
"""
function get_Ynet(
    edge_data_from_json,
    shunt_data_from_json;
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0,
    line_data_in_pu = true)

    #--------------------------------------

    (edges_fbus,
     edges_tbus,
     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type) =
         get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    (shunt_idx,
     Gs,
     Bs) =
        get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
            shunt_data_from_json)

    #--------------------------------------

    all_nodes_idx =
        sort(unique([edges_fbus;edges_tbus]))

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx  )
    
    #--------------------------------------

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    transformed_edges_fbus =
        [ n2s_all_nodes_idx[idx] for idx in edges_fbus  ]

    transformed_edges_tbus =
        [  n2s_all_nodes_idx[idx] for idx in edges_tbus ]
    
    #--------------------------------------
    
    if line_data_in_pu == true

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                all_nodes_idx,
                Gs, Bs;
                baseMVA =
                    baseShunt )
        
        edges_orientation =
            get_edges_orientation_by_generic(
                edges_fbus,
                edges_tbus )

        nodes_incident_edges =
            get_nodes_incident_edges_by_orientations(
                edges_orientation )

        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = 1.0,
                basekV = basekV )

    else

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                all_nodes_idx,
                Gs, Bs;
                baseMVA =
                    baseShunt )

        edges_orientation =
            get_edges_orientation_by_generic(
                edges_fbus,
                edges_tbus )

        nodes_incident_edges =
            get_nodes_incident_edges_by_orientations(
                edges_orientation )

        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = baseMVA,
                basekV = basekV )
    end
    

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
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx] ]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(all_nodes_idx,
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
        for (node_idx, orientations_and_edges_Ybr) in 
            zip(all_nodes_idx,
                nodes_incident_edges_orientation_and_Ybr_cal)]
    

    Ynet =
        Vector{ComplexF64}[
            [  idx == 1 ?
                last(node_k_idx_and_shunt) +
                Ynet_k_element :
                Ynet_k_element
               for (idx, Ynet_k_element) in
                   enumerate( Ynet_no_shunt[
                       n2s_all_nodes_idx[
                           first(node_k_idx_and_shunt)]] )
                   ]
            for node_k_idx_and_shunt in
                nodes_idx_and_Yshunt ]
    
    return (;Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
end

# @doc (@doc get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_edges_data)
"""
    get_Ynet(
        edge_data_from_json,
        shunt_data_from_json;
        baseMVA = 1.0,
        basekV = 1.0,
        baseShunt = 1.0,
        line_data_in_pu = true)


Returns a namedtuple of network admitance vectors and network nodes neigbouhood vectors.
"""
function get_Ynet_sp_sh(
    edge_data_from_json,
    shunt_data_from_json;
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0,
    line_data_in_pu = true )

    #--------------------------------------

    (edges_fbus,
     edges_tbus,
     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type) =
         get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    (shunt_idx,
     Gs,
     Bs) =
        get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
            shunt_data_from_json)

    #--------------------------------------

    all_nodes_idx =
        sort(unique([edges_fbus;edges_tbus]))

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx  )
    
    #--------------------------------------

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    transformed_edges_fbus =
        [ n2s_all_nodes_idx[idx] for idx in edges_fbus  ]

    transformed_edges_tbus =
        [  n2s_all_nodes_idx[idx] for idx in edges_tbus ]
    
    #--------------------------------------
    
    if line_data_in_pu == true
        
        (;y_sh_shunt_exist,
         nz_y_sh_idxs,
         nz_y_sh ) =
             get_nodes_idx_and_Yshunt_non_zero(
                 Gs,
                 Bs;
                 baseMVA =
                     baseShunt )
        
        
        edges_orientation =
            get_edges_orientation_by_generic(
                edges_fbus,
                edges_tbus )

        nodes_incident_edges =
            get_nodes_incident_edges_by_orientations(
                edges_orientation )

        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = 1.0,
                basekV = basekV )

    else

        (;y_sh_shunt_exist,
         nz_y_sh_idxs,
         nz_y_sh ) =
             get_nodes_idx_and_Yshunt_non_zero(
                 Gs,
                 Bs;
                 baseMVA =
                     baseShunt )

        edges_orientation =
            get_edges_orientation_by_generic(
                edges_fbus,
                edges_tbus )

        nodes_incident_edges =
            get_nodes_incident_edges_by_orientations(
                edges_orientation )

        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = baseMVA,
                basekV = basekV )
    end
    
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
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx] ]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(all_nodes_idx,
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
        for (node_idx, orientations_and_edges_Ybr) in 
            zip(all_nodes_idx,
                nodes_incident_edges_orientation_and_Ybr_cal)]
    
    #-----------------------------------------

    n2s_nz_y_sh_idxs =
        get_n2s_any( nz_y_sh_idxs )

    #-----------------------------------------

    Ynet = y_sh_shunt_exist != true ? Ynet_no_shunt :
        Vector{ComplexF64}[ node_idx ∈ nz_y_sh_idxs ? 
        [ idx == 1 ?
        Ynet_k_element + nz_y_sh[
            n2s_nz_y_sh_idxs[node_idx]] : Ynet_k_element
               for (idx, Ynet_k_element) in
                   enumerate( Ynet_no_shunt[
                       n2s_all_nodes_idx[ node_idx ] ] )
                   ] : Ynet_no_shunt[ n2s_all_nodes_idx[
                       node_idx]]
                            for node_idx in all_nodes_idx  ]
    
    return (;Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
end

#--------------------------------------


"""
    get_Yπ_net(
        edge_data_from_json,
        shunt_data_from_json;
        # all_nodes_idx,
        # n2s_all_nodes_idx,
        baseMVA = 1.0,
        basekV = 1.0,
        baseShunt = 1.0,
        line_data_in_pu = true,
        orientated_bool = false )


Returns nodes incident edges elementary admittance matrices `yπ`.


swap_yπ_diagonal_elements is used to ensure a
proper orientation of elementary `yπ` based on an
edge orientation `(i,j)` i.e (from, to).
Since an edge is connected to two nodes; node i and
node k, a primary node checks if the first index
in `(i,j)` is the same as its index. if from idx
is the same as the index of the primary node, `yπ`
is accepted, otherwise, values in the diagonal
of `yπ` are swapped. This is important for branches
 element whose elementaray admittance `yπ` matrix
are not symetrix

A primary node is the first node in each row of
`nodes_idx_with_adjacent_nodes_idx`

"""
function get_Yπ_net(
    edge_data_from_json,
    shunt_data_from_json;
    # all_nodes_idx,
    # n2s_all_nodes_idx,
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0,
    line_data_in_pu = true,
    orientated_bool = false )

    #--------------------------------------

    (edges_fbus,
     edges_tbus,
     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type) =
         get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    #--------------------------------------

    (shunt_idx,
     Gs,
     Bs) =
        get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
            shunt_data_from_json)

    #--------------------------------------


    all_nodes_idx =
        sort(unique([edges_fbus;edges_tbus]))

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx  )
    
    #--------------------------------------

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    transformed_edges_fbus =
        [ n2s_all_nodes_idx[idx] for idx in edges_fbus  ]

    transformed_edges_tbus =
        [  n2s_all_nodes_idx[idx] for idx in edges_tbus ]
    
    #--------------------------------------

    if line_data_in_pu == true

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                transformed_all_nodes_idx,
                Gs, Bs;
                baseMVA =
                    baseShunt )

        Yshunt = last.(nodes_idx_and_Yshunt)

        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = 1.0,
                basekV = basekV )        
    else

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                transformed_all_nodes_idx,
                Gs, Bs;
                baseMVA =
                    baseShunt )
        
        Yshunt = last.(nodes_idx_and_Yshunt)

        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = baseMVA,
                basekV = basekV )
        
    end
    
    edges_orientation =
        get_edges_orientation_by_generic(
            edges_fbus,
            edges_tbus )

    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )

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
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx] ]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(all_nodes_idx,
                nodes_incident_edges_and_orientation)]
    
    #------------------------------------------

    nodes_incident_edges_orientation_and_Ybr_cal = [
        edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges]

    if orientated_bool == true

        Yπ_net = [
            [ node_idx == first(first(
                orientations_and_edges_Ybr[k] )) ?
                    last(orientations_and_edges_Ybr[k]) :
                    swap_yπ_diagonal_elements(
                    last(orientations_and_edges_Ybr[k])) 
              for k in 1:length(
                  orientations_and_edges_Ybr ) ]
            for (node_idx, orientations_and_edges_Ybr) in 
                zip(all_nodes_idx,
                    nodes_incident_edges_orientation_and_Ybr_cal ) ]

        
    else

        Yπ_net = [
            [ last.(orientations_and_edges_Ybr) ]
            for orientations_and_edges_Ybr in 
                    nodes_incident_edges_orientation_and_Ybr_cal  ]
        
    end
    
    return (;Yπ_net,
         Yshunt,
         nodes_idx_with_adjacent_nodes_idx )
    
end

#--------------------------------------

"""
    get_Ybus(
        edge_data_from_json,
        shunt_data_from_json;
        basekV = 1.0,
        baseMVA = 1.0,
        line_data_in_pu = true )


Returns network sparse admittance matrix.
"""
function get_Ybus(
    edge_data_from_json,
    shunt_data_from_json;
    basekV = 1.0,
    baseMVA = 1.0,
    line_data_in_pu = true )
    
    (edges_fbus,
     edges_tbus,
     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type) =
         get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    # (Gs, Bs) =
    #     get_nodes_shunts_Gs_and_Bs_by_json(
    #         shunt_data_from_json)

    (shunt_idx,
     Gs,
     Bs) =
        get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
            shunt_data_from_json)

    #--------------------------------------

    all_nodes_idx =
        sort(unique([edges_fbus;edges_tbus]))

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx  )
    
    #--------------------------------------

    transformed_edges_fbus =
        [ n2s_all_nodes_idx[idx] for idx in edges_fbus  ]

    transformed_edges_tbus =
        [  n2s_all_nodes_idx[idx] for idx in edges_tbus ]
    
    #--------------------------------------

    no_nodes = length(Gs)

    no_edges = length(edges_r)

    incidence_matrix  =
        SparseArrays.sparse(
            transformed_edges_fbus,
            1:no_edges, 1, no_nodes, no_edges) +

                SparseArrays.sparse(
                    transformed_edges_tbus,
                    1:no_edges, -1, no_nodes, no_edges)

    #------------------------------------------------
    
    ys =
        1.0 ./ (edges_r + im * edges_x) 
       

    y_c =
        1 / 2 * (im *  edges_b)

    y_sh =  (Gs .+ im *  Bs )

    if line_data_in_pu == true

        baseZ = basekV^2/baseMVA

        baseY = 1/baseZ


        Gs = Gs ./  baseY
        Bs = Bs ./  baseY        

        ys   = baseY ./ (edges_r + im * edges_x) 
        
        y_c  = 1 / 2 * (im *  edges_b) * baseY 

        y_sh =  (Gs .+ im *  Bs )

        
    end

    inv_τ =
        [ (a_ratio == 0.0 || a_ratio == 0) ? 1.0 : 1/a_ratio
              for a_ratio in
                  edges_ratio  ]
    
    θ_shift = edges_angle
    
    Yff = (ys +  y_c) .* (inv_τ).^2

    Ytf = -ys .* (inv_τ ./ exp.(im * θ_shift))

    Yft = -ys .* (inv_τ ./ exp.(-im * θ_shift))

    Ytt = (ys +  y_c)

        
    Y_sh = SparseArrays.spdiagm(y_sh)

    
    Cf = SparseArrays.sparse(
        1:no_edges,
        transformed_edges_fbus,
             1,  no_edges, no_nodes) 

    Ct = SparseArrays.sparse(
        1:no_edges,
        transformed_edges_tbus,
             1,  no_edges, no_nodes) 
    
    Yf = SparseArrays.spdiagm( Yff ) * Cf +
        SparseArrays.spdiagm( Yft )  * Ct

    Yt = SparseArrays.spdiagm( Ytf ) * Cf +
        SparseArrays.spdiagm( Ytt )  * Ct

    Ybus = Cf' * Yf + Ct' * Yt + Y_sh

    return (; Ybus,)
    
    
end


"""
    get_Ynet_and_related_vectors(
        edge_data_from_json,
        shunt_data_from_json;
        basekV = 1.0,
        baseMVA = 1.0,
        line_data_in_pu = true )


Returns constituent objects for building `Ynet`.
"""
function get_Ynet_and_related_vectors(
    edge_data_from_json,
    shunt_data_from_json;
    # all_nodes_idx,
    # n2s_all_nodes_idx,
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0,
    line_data_in_pu = true)

    #--------------------------------------

    (edges_fbus,
     edges_tbus,
     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type) =
         get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    #--------------------------------------
    
    # (Gs, Bs) =
    #     get_nodes_shunts_Gs_and_Bs_by_json(
    #         shunt_data_from_json)

    (shunt_idx,
     Gs,
     Bs) =
        get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
            shunt_data_from_json)


    #--------------------------------------

    all_nodes_idx =
        sort(unique([edges_fbus;edges_tbus]))

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx  )
    
    #--------------------------------------

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
    transformed_edges_fbus =
        [ n2s_all_nodes_idx[idx] for idx in edges_fbus  ]

    transformed_edges_tbus =
        [  n2s_all_nodes_idx[idx] for idx in edges_tbus ]
    
    #--------------------------------------
    
    if line_data_in_pu == true

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                transformed_all_nodes_idx,
                Gs, Bs;
                baseMVA =
                    baseShunt )

        edges_orientation =
            get_edges_orientation_by_generic(
                edges_fbus,
                edges_tbus )

        nodes_incident_edges =
            get_nodes_incident_edges_by_orientations(
                edges_orientation )

        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = 1.0,
                basekV = basekV )

    else

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                transformed_all_nodes_idx,
                Gs, Bs;
                baseMVA =
                    baseShunt )

        edges_orientation =
            get_edges_orientation_by_generic(
                edges_fbus,
                edges_tbus )

        nodes_incident_edges =
            get_nodes_incident_edges_by_orientations(
                edges_orientation )

        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = baseMVA,
                basekV = basekV )
    end
    

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


    # this address non continous nodes labeling
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx] ]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(transformed_all_nodes_idx,
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
        for (node_idx, orientations_and_edges_Ybr) in 
            zip(transformed_all_nodes_idx,
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

    Ynet_wt_nodes_idx_wt_adjacent_nodes =
        (;Ynet,
     nodes_idx_with_adjacent_nodes_idx )
    
    return (;nodes_idx_and_Yshunt ,
            edges_orientation,
            nodes_incident_edges,
            edges_Ybr_cal,
            edges_orientation_and_edges_Ybr_cal,
            nodes_incident_edges_and_orientation,
            nodes_idx_with_adjacent_nodes_idx,
            nodes_incident_edges_orientation_and_Ybr_cal,
            Ynet_no_shunt,
            Ynet_wt_nodes_idx_wt_adjacent_nodes )
    
end

#--------------------------------------


"""
    get_Ybus_and_related_matrices(
        edge_data_from_json,
        shunt_data_from_json;
        basekV = 1.0,
        baseMVA = 1.0,
        line_data_in_pu = true )


Returns constituent objects for building `Ybus`.
"""
function get_Ybus_and_related_matrices(
    edge_data_from_json,
    shunt_data_from_json;
    basekV = 1.0,
    baseMVA = 1.0,
    line_data_in_pu = true )

    #--------------------------------------

    (edges_fbus,
     edges_tbus,
     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type) =
         get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    (shunt_idx,
     Gs,
     Bs) =
        get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
            shunt_data_from_json)


    #--------------------------------------

    all_nodes_idx =
        sort(unique([edges_fbus;edges_tbus]))

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx  )
    
    #--------------------------------------

    transformed_edges_fbus =
        [ n2s_all_nodes_idx[idx] for idx in edges_fbus  ]

    transformed_edges_tbus =
        [  n2s_all_nodes_idx[idx] for idx in edges_tbus ]
    

    #--------------------------------------

    no_nodes = length(Gs)

    no_edges = length(edges_r)

    incidence_matrix  =
        SparseArrays.sparse(
            transformed_edges_fbus,
            1:no_edges, 1, no_nodes, no_edges) +

                SparseArrays.sparse(
                    transformed_edges_tbus,
                    1:no_edges, -1, no_nodes, no_edges)

    #------------------------------------------------
    
    ys =
        1.0 ./ (edges_r + im * edges_x) 
       

    y_c =
        1 / 2 * (im *  edges_b)

    y_sh =  (Gs .+ im *  Bs )

    if line_data_in_pu == true

        baseZ = basekV^2/baseMVA

        baseY = 1/baseZ


        Gs = Gs ./  baseY
        Bs = Bs ./  baseY        

        ys   = baseY ./ (edges_r + im * edges_x) 
        
        y_c  = 1 / 2 * (im *  edges_b) * baseY 

        y_sh =  (Gs .+ im *  Bs )

        
    end

    inv_τ =
        [(a_ratio==0.0 || a_ratio == 0) ? 1.0 : 1/a_ratio
              for a_ratio in
                  edges_ratio  ]
    
    θ_shift = edges_angle
    
    Yff = (ys +  y_c) .* (inv_τ).^2

    Ytf = -ys .* (inv_τ ./ exp.(im * θ_shift))

    Yft = -ys .* (inv_τ ./ exp.(-im * θ_shift))

    Ytt = (ys +  y_c)

        
    Y_sh = SparseArrays.spdiagm(y_sh)

    
    Cf = SparseArrays.sparse(
        1:no_edges,
        transformed_edges_fbus,
             1,  no_edges, no_nodes) 

    Ct = SparseArrays.sparse(
        1:no_edges,
        transformed_edges_tbus,
             1,  no_edges, no_nodes) 
    
    Yf = SparseArrays.spdiagm( Yff ) * Cf +
        SparseArrays.spdiagm( Yft )  * Ct

    Yt = SparseArrays.spdiagm( Ytf ) * Cf +
        SparseArrays.spdiagm( Ytt )  * Ct

    Ybus = Cf' * Yf + Ct' * Yt + Y_sh

    return (;Yff,
            Ytf,
            Yft,
            Ytt,
            Y_sh,
            Cf,
            Ct,
            Yf,
            Yt,
            Ybus)
    
    
end

#--------------------------------------

"""
    get_Yπ_net_and_related_vectors(
        edge_data_from_json,
        shunt_data_from_json;
        all_nodes_idx,
        n2s_all_nodes_idx,
        basekV = 1.0,
        baseMVA = 1.0,
        baseShunt = 1.0,
        line_data_in_pu = true,
        orientated_bool = true )


Returns constituent objects for building `Yπ_net`.
"""
function get_Yπ_net_and_related_vectors(
    edge_data_from_json,
    shunt_data_from_json;
    all_nodes_idx,
    n2s_all_nodes_idx,
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0,
    line_data_in_pu = true,
    orientated_bool = true )

    #--------------------------------------

    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]

    #--------------------------------------

    (edges_fbus,
     edges_tbus,
     edges_r,
     edges_x,
     edges_b,
     edges_ratio,
     edges_angle,
     edges_type) =
         get_edges_ftbus_and_generic_data_by_json(
             edge_data_from_json )

    # (Gs, Bs) =
    #     get_nodes_shunts_Gs_and_Bs_by_json(
    #         shunt_data_from_json)

    (shunt_idx,
     Gs,
     Bs) =
        get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
            shunt_data_from_json)

    #--------------------------------------

    if line_data_in_pu == true

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                transformed_all_nodes_idx,
                Gs, Bs;
                baseMVA =
                    baseShunt )

        Yshunt = last.(nodes_idx_and_Yshunt)


        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = 1.0,
                basekV = basekV )        
    else

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                transformed_all_nodes_idx,
                Gs, Bs;
                baseMVA =
                    baseShunt )

        Yshunt = last.(nodes_idx_and_Yshunt)


        edges_Ybr_cal =
            get_edges_Ybr_by_generic(
                edges_r, edges_x, edges_b,
                edges_ratio, edges_angle,
                edges_type,
                Gs, Bs;
                baseMVA = baseMVA,
                basekV = basekV )
        
    end
        
    edges_orientation =
        get_edges_orientation_by_generic(
            edges_fbus,
            edges_tbus )

    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )

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

    # this address non continous nodes labeling
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx] ]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(transformed_all_nodes_idx,
                nodes_incident_edges_and_orientation)]
    
    #------------------------------------------

    nodes_incident_edges_orientation_and_Ybr_cal =[
        edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges]

    if orientated_bool == true

        Yπ_net = [
            [ node_idx == first(first(
                orientations_and_edges_Ybr[k] )) ?
                    last(orientations_and_edges_Ybr[k]) :
                    swap_yπ_diagonal_elements(
                    last(orientations_and_edges_Ybr[k])) 
              for k in 1:length(
                  orientations_and_edges_Ybr ) ]
            for (node_idx, orientations_and_edges_Ybr) in 
                zip(transformed_all_nodes_idx,
                    nodes_incident_edges_orientation_and_Ybr_cal ) ]

        
    else

        Yπ_net = [
            [ last(orientations_and_edges_Ybr[k])
              for k in 1:length(
                  orientations_and_edges_Ybr ) ]
            for (node_idx, orientations_and_edges_Ybr) in 
                zip(transformed_all_nodes_idx,
                    nodes_incident_edges_orientation_and_Ybr_cal ) ]
        
    end

    Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes =
        (;Yπ_net,
         Yshunt,
         nodes_idx_with_adjacent_nodes_idx )
    
    return (;nodes_idx_and_Yshunt ,
            edges_orientation,
            nodes_incident_edges,
            edges_Ybr_cal,
            edges_orientation_and_edges_Ybr_cal,
            nodes_incident_edges_and_orientation,
            nodes_idx_with_adjacent_nodes_idx,
            nodes_incident_edges_orientation_and_Ybr_cal,
            Yπ_net,
            Yπ_net_Yshunt_wt_nodes_idx_wt_adjacent_nodes )
    
end


#---------------------------------------------------

"""
    get_a_node_type_tup_idx_PQ_data_by_generic(
        idx, Pd, Qd;
        baseMVA = 1.0,
        node_type = :load )


Returns a tuple of `(idx, Bus, P, Q, node_type)` for non-generation node.
"""
function get_a_node_type_tup_idx_PQ_data_by_generic(
    idx, Pd, Qd;
    baseMVA = 1.0,
    node_type = :load )

    return (idx,
            (Bus = "bus$(idx)",
             P   = Pd/baseMVA,
             Q   = Qd/baseMVA),
            node_type)
end

#---------------------------------------------------
# data by mpc
#---------------------------------------------------

"""
    get_shunt_data_by_mpc(
        mpc_bus;
        mpc_baseMVA =
            1.0)


Returns a dict `Dict{:shunt_idx,:shunt_Gs,:shunt_Bs }` for network lines shunt parameters.
"""
function get_shunt_data_by_mpc(
    mpc_bus;
    mpc_baseMVA =
        1.0)

    return Dict(:shunt_idx => mpc_bus.bus_i,
         :shunt_Gs => mpc_bus.Gs,
         :shunt_Bs => mpc_bus.Bs)

end


"""
    get_gencost_data_by_mpc(
        mpc_gencost;
        mpc_baseMVA =
            1.0)


Returns a dict `Dict{:cost_type,:startup,:shutdown,:n,:c_n_1,:c_1,:c_0 }` for generators operation cost.
"""
function get_gencost_data_by_mpc(
    mpc_gencost;
    mpc_baseMVA =
        1.0)

    return Dict(:cost_type => mpc_gencost.var"2",
         :startup => mpc_gencost.startup,
                :shutdown => mpc_gencost.shutdown,
                :n => mpc_gencost.n,
                :c_n_1 => mpc_gencost.var"c(n-1)",
                :c_1 => mpc_gencost.var"...",
                :c_0 => mpc_gencost.c0)

end


function get_branches_impedance_data_by_components_by_mpc(
    mpc_branch;
    vec_edge_type = [],
    use_inductance = false,
    ωs = 376.991  )

    if length(vec_edge_type) == 0
        vec_comp_type =
            [(ratio == 0.0 || ratio == 0) ?
            :PiModelLine : :Transformer
             for ratio in mpc_branch.ratio ]
    else
        vec_comp_type = vec_edge_type
    end
    

    if use_inductance == false
        
       return [(ratio == 0.0 || ratio == 0) ?
           (idx = idx,
            components_data = (;fbus, tbus, r, x, b, ratio,
                         angle, status ),
            components_type = (; edge_type) ) :
                (idx = idx,
                 components_data =
                     (;fbus,tbus,r,x,b,ratio,angle,status),
                 components_type =
                     (; edge_type) )
               for (idx,fbus,tbus,r,x,b,ratio,angle,status,
                    edge_type) in
                 zip(1:length(mpc_branch.fbus),
                     mpc_branch.fbus,
                     mpc_branch.tbus,
                     mpc_branch.r,
                     mpc_branch.x,
                     mpc_branch.b,
                     mpc_branch.ratio,
                     mpc_branch.angle,
                     mpc_branch.status,
                     vec_comp_type)]
    else

        branch_inductance = (1/ωs) .* mpc_branch.x

        branch_capacitance = (1/ωs) .* mpc_branch.b
        
       return [(ratio == 0.0 || ratio == 0) ?
           (idx = idx,
            components_data =
                (;fbus, tbus,
                 r, L, C,
                 ratio,  angle,
                 status ),
            components_type = (; edge_type) ) :
                (idx = idx,
                 components_data =
                     (;fbus, tbus,
                      r, L, C,
                      ratio, angle, status),
                 components_type =
                     (; edge_type) )
               for (idx,
                    fbus, tbus,
                    r, x, b,
                    ratio, angle,
                    status, edge_type) in
                 zip(1:length(mpc_branch.fbus),
                     mpc_branch.fbus,
                     mpc_branch.tbus,
                     mpc_branch.r,
                     branch_inductance,
                     branch_capacitance,
                     mpc_branch.ratio,
                     mpc_branch.angle,
                     mpc_branch.status,
                     vec_comp_type)]
    end
        
end


function get_branches_data_and_types_by_mpc(
    mpc_branch;
    mpc_baseMVA=1.0, basekV = 1.0 )

  return  [ (ratio == 0.0 || ratio == 0) ?
      (idx, get_edge_y_line_data_by_generic(
          fbus, tbus, r, x, b, ratio, angle,
          status; baseMVA=mpc_baseMVA,basekV = basekV), :line) :
              (idx, get_edge_y_transformer_data_by_generic(
                  fbus, tbus,
                  r, x, b,
                  ratio, angle,
                  status;
                  baseMVA=mpc_baseMVA, basekV = basekV ),
               :transformer )
            for (idx,
                 fbus, tbus,
                 r, x, b,
                 ratio, angle,
                 status) in
            zip(1:length(mpc_branch.fbus),
                mpc_branch.fbus,
                mpc_branch.tbus,
                mpc_branch.r,
                mpc_branch.x,
                mpc_branch.b,
                mpc_branch.ratio,
                mpc_branch.angle,
                mpc_branch.status )]
    
    
end

"""
    get_loc_load_tup_idx_PQ_data_by_mpc(
        idx, Pd, Qd;
        mpc_baseMVA = 1.0)


Returns a tuple of `(idx, Bus, loc_P, loc_Q, node_type)` for generation node local load.
"""
function get_loc_load_tup_idx_PQ_data_by_mpc(
    idx,
    Pd,
    Qd;
    mpc_baseMVA = 1.0 )

    return (idx,
            (Bus   = "bus$(idx)",
             loc_P = Pd/mpc_baseMVA,
             loc_Q = Qd/mpc_baseMVA),
            :loc_load)
end

"""
    get_non_gen_node_tup_idx_PQ_data_by_mpc(
        idx, Pd, Qd;
        mpc_baseMVA = 1.0)


Returns a tuple of `(idx, Bus, P, Q)` for non-generation node.
"""
function get_non_gen_node_tup_idx_PQ_data_by_mpc(
    idx,
    Pd,
    Qd;
    mpc_baseMVA = 1.0 )

    return (idx,
            (Bus = "bus$(idx)",
             P   = Pd/mpc_baseMVA,
             Q   = Qd/mpc_baseMVA))
end


"""
    get_gen_node_static_data_tup_by_mpc(
        idx,
        vmax,vmin,
        Pg,Qg,
        Vg,
        Qmax, Qmin,
        Pmax, Pmin;
        mpc_baseMVA = 1.0 )


Returns a generator's index and static parameters as a tuple of index and namedtuple of static parameters.

"""
function get_gen_node_static_data_tup_by_mpc(
    idx,
    vmax,vmin,
    Pg,Qg,
    Vg,
    Qmax, Qmin,
    Pmax, Pmin;
    mpc_baseMVA = 1.0 )

    return (idx,
            (Bus  = "bus$(idx)",
             P    = Pg/mpc_baseMVA,
             Q    = Qg/mpc_baseMVA,             
             vh   = Vg,
             vmax = vmax,
             vmin = vmin,
             Qmax = Qmax/mpc_baseMVA,
             Qmin = Qmin/mpc_baseMVA,
             Pmax = Pmax/mpc_baseMVA,
             Pmin = Pmin/mpc_baseMVA,
             Sn   = (sqrt(Pmax^2 + Qmax^2))/mpc_baseMVA
             ),
            :generator )
end


"""
    get_gen_sub_static_data_tup_by_mpc(
        idx,
        n2s_gens_idx,
        mpc_gen )


Returns a tuple of static parameters for a generator.

"""
function get_gen_sub_static_data_tup_by_mpc(
    idx,
    n2s_gens_idx,
    mpc_gen )

    return (mpc_gen.Pg[n2s_gens_idx[idx]],
            mpc_gen.Qg[n2s_gens_idx[idx]],
            mpc_gen.Vg[n2s_gens_idx[idx]],
            mpc_gen.Qmax[n2s_gens_idx[idx]],
            mpc_gen.Qmin[n2s_gens_idx[idx]],
            mpc_gen.Pmax[n2s_gens_idx[idx]],
            mpc_gen.Pmin[n2s_gens_idx[idx]] )
end


"""
    get_gen_nodes_static_tup_data_by_mpc(
        mpc_bus,
        mpc_gen;
        mpc_baseMVA=1.0 )

Returns a list of static parameters for single or multi generators per node for all generation nodes.
"""
function get_gen_nodes_static_tup_data_by_mpc(
    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0 )

    #------------------------------------------

    (;multi_gens_idx,
     multi_gen_bool,
     multi_gens_nodes) =
         NamedTupleTools.select(
             get_multi_gens_idx_wt_multi_gen_bool(
                 mpc_gen; sorted_bool = false ),
             (:multi_gens_idx,
              :multi_gen_bool,
              :multi_gens_nodes) )

    n2s_multi_gens_idx =
        get_n2s_any( multi_gens_idx )
    
    #------------------------------------------
     
 
    """
    The right idxs to use is the unsorted idx, since they
    are the one aligned with gens static parameters in mpc.
    `gen_node_static_data_tup_by_mpc` can then be sorted
    latter by idx
    """

    gens_nodes_idx =mpc_gen.bus

    n2s_gens_idx =
        get_n2s_any( gens_nodes_idx )

        
    all_nodes_idx =
        get_all_nodes_idx_by_mpc(
            mpc_bus;
            sorted_bool = false)

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx)

    # @show  gens_nodes_idx

    # @show n2s_gens_idx
    
    if multi_gen_bool == true

        
        return [get_gen_node_static_data_tup_by_mpc(
            idx,
            mpc_bus.Vmax[
                n2s_all_nodes_idx[get_node_idx_multi(idx)]],
            mpc_bus.Vmin[
                n2s_all_nodes_idx[get_node_idx_multi(idx)]],
            get_gen_sub_static_data_tup_by_mpc(
                idx, n2s_multi_gens_idx, mpc_gen )...
                    ;mpc_baseMVA = mpc_baseMVA)
                 for idx in multi_gens_idx ]            

    else

        gen_node_static_data_tup_by_mpc =
            [get_gen_node_static_data_tup_by_mpc(
                idx,
                mpc_bus.Vmax[ n2s_all_nodes_idx[idx]],
                mpc_bus.Vmin[ n2s_all_nodes_idx[idx]],
                get_gen_sub_static_data_tup_by_mpc(
                    idx, n2s_gens_idx, mpc_gen )...
                        ;mpc_baseMVA = mpc_baseMVA)
             for idx in gens_nodes_idx ]

        return sort(gen_node_static_data_tup_by_mpc,
                    by = x -> x[1] ) 
    end
    

end

"""
    get_gen_nodes_dict_static_data_by_mpc(
        mpc_bus,
        mpc_gen;
        mpc_baseMVA = 1.0 )

Returns an OrderedDict of static data indexed by generators indices.
"""
function get_gen_nodes_dict_static_data_by_mpc(
    mpc_bus,
    mpc_gen;
    mpc_baseMVA = 1.0 )

    gen_nodes_static_tup =
        get_gen_nodes_static_tup_data_by_mpc(
            mpc_bus,
            mpc_gen;
            mpc_baseMVA=mpc_baseMVA)
     
    tup_type = NamedTuple
    
    return OrderedDict{
        Union{Int64,String,Symbol},
        tup_type}(
            idx => node_data_tup
            for (idx, node_data_tup) in
                zip(first.(gen_nodes_static_tup),
                    second.(gen_nodes_static_tup)) )
    
end


"""
    get_nodes_static_tup_data_by_mpc(
        mpc_bus,
        mpc_gen;
        mpc_baseMVA = 1.0 )

Returns tuples of static data indexed by generating and non-generating plants indices.
"""
function get_nodes_static_tup_data_by_mpc(
    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0 )

    #------------------------------------------

    load_nodes_idx =
        get_load_nodes_idx_by_mpc(
            mpc_bus ;
            sorted_bool = false )

    transmission_nodes_idx =
        get_transmission_nodes_idx_by_mpc(
            mpc_bus;
            sorted_bool = false)

    gens_nodes_idx =  mpc_gen.bus
    # gens_nodes_idx =
    #     get_gens_nodes_idx_by_mpc( mpc_bus)

    all_nodes_idx =
        get_all_nodes_idx_by_mpc(
            mpc_bus;
            sorted_bool = false)

    n2s_gens_idx =
        get_n2s_any(gens_nodes_idx)

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx)

    #------------------------------------------

    (;multi_gens_idx,
     multi_gen_bool,
     multi_gens_nodes,
     gens_countmap) =
         NamedTupleTools.select(
             get_multi_gens_idx_wt_multi_gen_bool(
                 mpc_gen),
             (:multi_gens_idx,
              :multi_gen_bool,
              :multi_gens_nodes,
              :gens_countmap) )

    #------------------------------------------

    if multi_gen_bool == true

        multi_gens_all_nodes_idx =
            get_multi_gens_all_nodes_idx(
                multi_gens_idx,
                multi_gens_nodes,
                gens_countmap,
                all_nodes_idx)
    
        n2s_multi_gens_idx =
            get_n2s_any( multi_gens_idx )

        n2s_multi_gens_all_nodes_idx =
            get_n2s_any( multi_gens_all_nodes_idx )

        return [ idx ∈ String.(load_nodes_idx) ?
            get_a_node_type_tup_idx_PQ_data_by_generic(
                get_node_idx_multi(idx),
                mpc_bus.Pd[
                    n2s_all_nodes_idx[
                    get_node_idx_multi(idx)] ],
                mpc_bus.Qd[
                    n2s_all_nodes_idx[
                    get_node_idx_multi(idx)] ];
                baseMVA =
                    mpc_baseMVA,
                node_type = :load ) :
                    idx ∈ String.(transmission_nodes_idx) ?
                    get_a_node_type_tup_idx_PQ_data_by_generic(
                        get_node_idx_multi(idx),
                        mpc_bus.Pd[
                            n2s_all_nodes_idx[
                                get_node_idx_multi(idx)]],
                        mpc_bus.Qd[
                            n2s_all_nodes_idx[
                                get_node_idx_multi(idx)]];
                        baseMVA =
                            mpc_baseMVA,
                        node_type = :transmission ) :
                            get_gen_node_static_data_tup_by_mpc(
                                idx,
                                mpc_bus.Vmax[
                                    n2s_all_nodes_idx[
                                       get_node_idx_multi(idx)]],
                                mpc_bus.Vmin[
                                    n2s_all_nodes_idx[
                                       get_node_idx_multi(idx)]],
                            get_gen_sub_static_data_tup_by_mpc(
                                idx,
                                n2s_multi_gens_idx,
                                mpc_gen )...;
                                mpc_baseMVA =
                                    mpc_baseMVA )
                 for idx in multi_gens_all_nodes_idx ]
        
    else

        return [ idx ∈ load_nodes_idx ?
            get_a_node_type_tup_idx_PQ_data_by_generic(
                idx, mpc_bus.Pd[ n2s_all_nodes_idx[ idx] ],
                mpc_bus.Qd[ n2s_all_nodes_idx[idx] ];
                baseMVA = mpc_baseMVA, node_type = :load ) :
                    idx ∈ transmission_nodes_idx ?
                    get_a_node_type_tup_idx_PQ_data_by_generic(
                        idx, mpc_bus.Pd[ n2s_all_nodes_idx[idx]],
                        mpc_bus.Qd[ n2s_all_nodes_idx[idx]];
                        baseMVA = mpc_baseMVA,
                        node_type = :transmission ) :
                            get_gen_node_static_data_tup_by_mpc(
                                idx,
                                mpc_bus.Vmax[
                                    n2s_all_nodes_idx[idx]],
                                mpc_bus.Vmin[
                                    n2s_all_nodes_idx[idx]],
                            get_gen_sub_static_data_tup_by_mpc(
                                idx, n2s_gens_idx, mpc_gen )...;
                                mpc_baseMVA = mpc_baseMVA )
                 for idx in all_nodes_idx ]
        
    end
    
end


"""
    get_nodes_dict_static_data_by_mpc(
        mpc_bus,
        mpc_gen;
        mpc_baseMVA = 1.0 )

Returns an OrderedDict of static data indexed by generators indices.
"""
function get_nodes_dict_static_data_by_mpc(
    mpc_bus,
    mpc_gen;
    mpc_baseMVA = 1.0 )

    nodes_static_tup =
        sort(get_nodes_static_tup_data_by_mpc(
            mpc_bus,
            mpc_gen;
            mpc_baseMVA=mpc_baseMVA) , by = x-> [x1])

    # nodes_static_tup =
    #     get_nodes_static_tup_data_by_mpc(
    #         mpc_bus,
    #         mpc_gen;
    #         mpc_baseMVA=mpc_baseMVA)
    

    tup_type = NamedTuple
    
    return OrderedDict{
        Union{Int64,String,Symbol}, tup_type}(
            idx => node_data_tup
            for (idx, node_data_tup) in
                zip(first.(nodes_static_tup),
                    second.(nodes_static_tup)) )
    
end


"""
     get_load_nodes_idx_wt_type_tup_by_mpc(
        mpc_load_type_data )


Returns tuple of load nodes indices, load nodes types and static data.
"""
function get_load_nodes_idx_wt_type_tup_by_mpc(
    mpc_load_type_data )

    # @show mpc_load_type_data

    # @show Tables.columntable(mpc_load_type_data)

    # @show Tables.rowtable(mpc_load_type_data)
    
    if typeof(mpc_load_type_data) == DataFrame

        mpc_load_type_data =
            Tables.rowtable(mpc_load_type_data)

        return [ (idx = a_node.idx,
                  load_type = a_node.load_type)
            for a_node in mpc_load_type_data ]  
    else

        return [
            (idx = a_node.idx,
             load_type = a_node.load_type)
            for a_node in mpc_load_type_data ]        
    end
    
end

"""
    get_load_nodes_static_tup_data_by_mpc(
        mpc_bus;
        mpc_baseMVA = 1.0)


Returns tuple of load nodes indices and static data.
"""
function get_load_nodes_static_tup_data_by_mpc(
    mpc_bus;
    mpc_baseMVA = 1.0)

    all_nodes_idx =
        get_all_nodes_idx_by_mpc( mpc_bus)

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx)
    
    load_nodes_idx =
        get_load_nodes_idx_by_mpc( mpc_bus )
    
    return [(idx,
         get_a_node_type_tup_idx_PQ_data_by_generic(
             idx, mpc_bus.Pd[ n2s_all_nodes_idx[idx]],
             mpc_bus.Qd[ n2s_all_nodes_idx[idx]];
             baseMVA = mpc_baseMVA,
             node_type = :load )) 
            for idx in load_nodes_idx ]
    
end

"""
    get_transmission_nodes_static_tup_data_by_mpc(
        mpc_bus;
        mpc_baseMVA = 1.0)


Returns tuple of transmission nodes indices and static data.
"""
function get_transmission_nodes_static_tup_data_by_mpc(
    mpc_bus;
    mpc_baseMVA=1.0)


    all_nodes_idx =
        get_all_nodes_idx_by_mpc( mpc_bus)

    n2s_all_nodes_idx =
        get_n2s_any( all_nodes_idx)
    
    transmission_nodes_idx =
        get_transmission_nodes_idx_by_mpc(
            mpc_bus)

    if length(transmission_nodes_idx) != 0
        return [(idx,
                 get_a_node_type_tup_idx_PQ_data_by_generic(
                    idx, mpc_bus.Pd[ n2s_all_nodes_idx[idx] ],
                     mpc_bus.Qd[ n2s_all_nodes_idx[ idx ] ];
                     baseMVA = mpc_baseMVA,
                     node_type = :transmission )) 
                for idx in transmission_nodes_idx ]
    else

        return nothing
    end
    
    
end


#---------------------------------------------------


"""
    get_static_loc_loads_json_data_storage_format(
        loc_loads_idx_and_locP_locQ_data;
        type_loc_load = "loc_Load_t1")


Returns OrderedDict of namedtuples for generators local loads indexed by generators indices.
"""
function get_static_loc_loads_json_data_storage_format(
    loc_loads_idx_and_locP_locQ_data;
    type_loc_load = "loc_Load_t1")

    loc_loads_exist =
        length(loc_loads_idx_and_locP_locQ_data) == 0 ?
        false : true

    return loc_loads_exist == false ?  [] : 
        OrderedDict{Union{Int64,String,Symbol},
                    NamedTuple}(
            String(split(lowercase(first(a_tup)),"bus")[2]) =>
                (loc_P = second(a_tup),
                 loc_Q = third( a_tup) )
            for a_tup in
                loc_loads_idx_and_locP_locQ_data )

end


"""
    get_static_load_nodes_json_data_storage_format(
        load_nodes_static_tup_data,
        load_nodes_idx_wt_type_tup;
        plant_type = "plant_PQ_Const_I")

Returns namedtuples of static load nodes data
"""
function get_static_load_nodes_json_data_storage_format(
    load_nodes_static_tup_data,
    load_nodes_idx_wt_type_tup;
    plant_type = "plant_PQ_Const_I")


    return [
        (idx = first(idx_wt_type_tup),
         plant_type = plant_type,
         components_type =
             (load = second(idx_wt_type_tup),),
         components_data =
             (load = second(idx_wt_data_tup),))
        for (idx_wt_type_tup, idx_wt_data_tup) in
            zip(load_nodes_idx_wt_type_tup,
                load_nodes_static_tup_data)]

end


"""
    get_static_transmission_nodes_json_data_storage_format(
        transmission_nodes_static_tup_data,
        load_nodes_idx_wt_type_tup;
        plant_type = "plant_Transmission_t2")

Returns namedtuples of static transmission nodes data.
"""
function get_static_transmission_nodes_json_data_storage_format(
    transmission_nodes_static_tup_data,
    load_nodes_idx_wt_type_tup;
    plant_type = "plant_Transmission_t2")

    transmission_nodes_exist =
        transmission_nodes_static_tup_data == nothing ?
         false : true    

    return transmission_nodes_exist == false ? [] :
        [(idx = first(idx_wt_type_tup),
          plant_type = plant_type,
          components_type =
              (transmission = second(idx_wt_type_tup),),
          components_data =
              (transmission = second(idx_wt_data_tup),))
         for (idx_wt_data_tup, idx_wt_data_tup) in
             zip(load_nodes_idx_wt_type_tup,
                 transmission_nodes_static_tup_data) ]
    
end


"""
    loc_load_exist_bool_by_mpc( mpc_bus )


Returns `true` or `false` if any of the generators has a local load.
"""
function loc_load_exist_bool_by_mpc( mpc_bus )

    gens_nodes_with_loc_loads_idx =
        [a_node
         for (a_node,node_type,node_Pd,node_Qd) in
             zip(mpc_bus.bus_i,
                 mpc_bus.type,
                 mpc_bus.Pd,
                 mpc_bus.Qd)
             
             if (node_type == 3 || node_type == 2) &&
                 ( (node_Pd != 0.0 || node_Pd != 0) ||
                 (node_Qd != 0.0 || node_Qd != 0) ) ]
       
        return length(gens_nodes_with_loc_loads_idx) != 0 ?
            true : false

end


function get_loc_loads_idx_and_locP_locQ_data_by_mpc(
    mpc_bus;
    mpc_baseMVA = 1.0 )

    loc_load_exist =
        loc_load_exist_bool_by_mpc(
            mpc_bus )
    
    if loc_load_exist == true

        return [(Bus = "bus$(idx)",
                  loc_P = node_Pd/mpc_baseMVA,
                  loc_Q = node_Qd/mpc_baseMVA)
                for (idx, node_type, node_Pd, node_Qd) in
                    zip( mpc_bus.bus_i, mpc_bus.type,
                         mpc_bus.Pd, mpc_bus.Qd )
                    if (node_type == 3 || node_type == 2) &&
                        ((node_Pd != 0.0 || node_Pd != 0) ||
                        (node_Qd != 0.0 || node_Qd != 0))]
    else
        return nothing
    end
    
end

#---------------------------------------------------
# sorted indices
#---------------------------------------------------


function get_all_nodes_idx_by_mpc(
    mpc_bus; sorted_bool = true )

    return sorted_bool == false ? copy(mpc_bus.bus_i) : sort(
        copy(mpc_bus.bus_i))
end


# function get_all_nodes_idx_by_mpc(
#     mpc_bus)

#     return sort([a_node
#             for a_node in mpc_bus.bus_i ])
# end


function get_slack_gens_nodes_idx_by_mpc(
    mpc_bus;
    sorted_bool = true )

    return sorted_bool == false ? [a_node
            for (a_node, node_type) in
                zip(mpc_bus.bus_i,
                    mpc_bus.type)
             if node_type == 3 ] :  sort( [a_node
            for (a_node, node_type) in
                zip(mpc_bus.bus_i,
                    mpc_bus.type)
             if node_type == 3 ] )
end


function get_non_slack_gens_nodes_idx_by_mpc(
    mpc_bus;sorted_bool = true  )

    return sorted_bool == false ? [a_node
            for (a_node, node_type) in
                zip( mpc_bus.bus_i,
                     mpc_bus.type)
                if node_type == 2 ] : sort( [a_node
            for (a_node, node_type) in
                zip( mpc_bus.bus_i,
                     mpc_bus.type)
                if node_type == 2 ] )
end


function get_gens_nodes_idx_by_mpc(
    mpc_bus; sorted_bool = true )

    return sorted_bool == false ? [a_node for (a_node, node_type) in
                zip( mpc_bus.bus_i,
                     mpc_bus.type)
                if node_type == 3 || node_type == 2 ] :  sort( [a_node for (a_node, node_type) in
                zip( mpc_bus.bus_i,
                     mpc_bus.type)
                if node_type == 3 || node_type == 2 ] )
end


function get_non_gens_nodes_idx_by_mpc(
    mpc_bus; sorted_bool = true )

    return sorted_bool == false ? [a_node
            for (a_node, node_type) in
                zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 1 ] :  sort([a_node
            for (a_node, node_type) in
                zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 1 ] )
end


function get_non_slack_gens_and_non_gens_idx_by_mpc(
    mpc_bus; sorted_bool = true )

    return sorted_bool == false ? [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 2 || node_type == 1 ] :  sort( [a_node for (a_node, node_type) in
             zip( mpc_bus.bus_i, mpc_bus.type)
             if node_type == 2 || node_type == 1 ] )
end


function get_gens_nodes_with_loc_loads_idx_by_mpc(
    mpc_bus; sorted_bool = true)

    gens_nodes_with_loc_loads_idx =
        [a_node for (a_node,node_type,node_Pd,node_Qd) in
             zip( mpc_bus.bus_i, mpc_bus.type,
                  mpc_bus.Pd, mpc_bus.Qd )
             if (node_type == 3 || node_type == 2) &&
                 ( (node_Pd != 0.0 || node_Pd != 0) ||
                 (node_Qd != 0.0 || node_Qd != 0) ) ]
    
    return length(gens_nodes_with_loc_loads_idx) == 0 ? [] :
        sorted_bool == false ? gens_nodes_with_loc_loads_idx : sort(gens_nodes_with_loc_loads_idx)
    
end


function get_nodes_with_demands_idx_by_mpc(
    mpc_bus; sorted_bool = true )

    return sorted_bool == false ?
        [a_node for (a_node, node_Pd, node_Qd) in
             zip( mpc_bus.bus_i, mpc_bus.Pd, mpc_bus.Qd)
             if  (node_Pd != 0.0 || node_Pd != 0) ||
                 (node_Qd != 0.0 || node_Qd != 0)] : sort(
                     [a_node for (a_node, node_Pd, node_Qd) in
             zip( mpc_bus.bus_i, mpc_bus.Pd, mpc_bus.Qd)
             if  (node_Pd != 0.0 || node_Pd != 0) ||
                 (node_Qd != 0.0 || node_Qd != 0)] )
end


function get_load_nodes_idx_by_mpc(
    mpc_bus; sorted_bool = true )

    return sorted_bool == false ?
        [a_node for (a_node,node_type,node_Pd,node_Qd) in
             zip(mpc_bus.bus_i, mpc_bus.type, mpc_bus.Pd,
                 mpc_bus.Qd)
             if (node_type == 1) && (
                 ((node_Pd != 0.0 || node_Pd != 0) || (
                     node_Qd !=0.0 || node_Qd != 0)))] : sort(
                         [a_node
                          for (a_node,
                               node_type,
                               node_Pd,node_Qd) in
                              zip(mpc_bus.bus_i,
                                  mpc_bus.type,
                                  mpc_bus.Pd,
                                  mpc_bus.Qd)
                if (node_type == 1) && (
                    ((node_Pd != 0.0 || node_Pd != 0) || (
                        node_Qd != 0.0 || node_Qd != 0))) ])
end


function get_transmission_nodes_idx_by_mpc(
    mpc_bus; sorted_bool = true )

    transmission_nodes_idx =
        [a_node for (a_node,a_type, node_Pd, node_Qd) in
             zip(mpc_bus.bus_i, mpc_bus.type, mpc_bus.Pd,
                 mpc_bus.Qd)
             if (a_type == 1) &&
                 (( (node_Pd == 0.0 || node_Pd == 0) &&
                 (node_Qd == 0.0 || node_Qd == 0))) ]

    return  length(transmission_nodes_idx) == 0 ?  [] :
        sorted_bool == false ? transmission_nodes_idx : sort(
            transmission_nodes_idx)
end

#---------------------------------------------------
#---------------------------------------------------


function  get_edges_orientation_by_mpc(
    mpc_branch )

    return get_edges_orientation_by_generic(
        mpc_branch.fbus,
        mpc_branch.tbus  )
end

function get_edges_Ybr_by_mpc(
    mpc_bus,
    mpc_branch,
    mpc_baseMVA;
    basekV = 1.0 )


    edge_type = [ (a_ratio == 0.0 || a_ratio == 0) ?
        "PiModelLine" : "Transformer"
              for a_ratio in
                  mpc_branch.ratio  ]
    
    return get_edges_Ybr_by_generic(
        mpc_branch.r, mpc_branch.x,
        mpc_branch.b, mpc_branch.ratio,
        mpc_branch.angle, edge_type,
        mpc_bus.Gs, mpc_bus.Bs;
        baseMVA = mpc_baseMVA,
        basekV = basekV)    
end



function get_nodes_idx_and_Yshunt_by_mpc(
    mpc_bus,
    mpc_baseMVA )

    y_sh =
        get_nodes_Yshunt_by_generic(
            mpc_bus.Gs, mpc_bus.Bs;
            baseMVA =
                mpc_baseMVA )
    
    return [(node_idx, y_shunt)
            for ( node_idx, y_shunt) in
                zip(mpc_bus.bus_i, y_sh ) ]

end



function get_edges_Ybr_cal_and_edges_orientation_by_mpc(
    mpc_bus,
    mpc_branch,
    baseMVA;
    basekV = 1.0 )

    edges_Ybr_cal =
        get_edges_Ybr_by_mpc(
            mpc_bus,
            mpc_branch,
            mpc_baseMVA;
            basekV = basekV )

    edges_orientation =
        get_edges_orientation_by_mpc(
            mpc_branch )

    return (; edges_Ybr_cal,
            edges_orientation )

    
end

function get_transmission_network_parameters_by_mpc(
    mpc_bus,
    mpc_branch,
    mpc_baseMVA;
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true )

    buses_idx = mpc_bus.bus_i
    
    edges_orientation =
        get_edges_orientation_by_mpc(
            mpc_branch )

    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )

    #-----------------------------------------

    "This gets the orientation of nodes incident edges"
    
    nodes_incident_edges_and_orientation =
        [ edges_orientation[
            a_node_incident_edges ]
          for a_node_incident_edges in
              nodes_incident_edges ]

    #-----------------------------------------
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx]]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(buses_idx,
                nodes_incident_edges_and_orientation)]

    #-----------------------------------------

    if line_data_in_pu == true

        already_in_pu = 1.0
        
    else

        pu_not_needed = 1.0
        baseMVA = mpc_baseMVA
        basekV  = basekV
        
    end
    
    if line_data_in_pu == true && use_pu_in_PQ == true

        edges_Ybr_cal =
            get_edges_Ybr_by_mpc(
                mpc_bus,
                mpc_branch,
                already_in_pu;
                basekV = already_in_pu )

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_mpc(
                mpc_bus,
                mpc_baseMVA )
        

    elseif line_data_in_pu == false && use_pu_in_PQ == true

        edges_Ybr_cal =
            get_edges_Ybr_by_mpc(
                mpc_bus,
                mpc_branch,
                mpc_baseMVA;
                basekV = basekV )

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_mpc(
                mpc_bus,
                mpc_baseMVA )

        
    elseif line_data_in_pu == false && use_pu_in_PQ == false

        edges_Ybr_cal =
            get_edges_Ybr_by_mpc(
                mpc_bus, mpc_branch, pu_not_needed;
                basekV = pu_not_needed )

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_mpc(
                mpc_bus,
                pu_not_needed )
        
    else

        throw("line_data_in_pu = $(line_data_in_pu) and use_pu_in_PQ = $(use_pu_in_PQ) option is not consistent. ")

    end
        
    nodes_Yshunt = second.(nodes_idx_and_Yshunt)
 
    #------------------------------------------

    edges_orientation_and_edges_Ybr_cal =
        [ (orient, a_Ybr)
          for (orient, a_Ybr) in
              zip( edges_orientation,
                   edges_Ybr_cal )] 

    nodes_incident_edges_orientation_and_Ybr_cal = [
        edges_orientation_and_edges_Ybr_cal[
            a_node_incident_edges ]
        for a_node_incident_edges in
            nodes_incident_edges]
    
    #------------------------------------------
    #------------------------------------------

    Ynet_no_shunt = [
        [ k == 0 ?
            sum( [ node_idx == first( first(
                orient_and_Ybr) ) ?
                    last(orient_and_Ybr)[1] :
                    last(orient_and_Ybr)[4]
                for orient_and_Ybr in
                   orientations_and_edges_Ybr]) :
                    node_idx == first(first(
                       orientations_and_edges_Ybr[ k])) ?
                     last(orientations_and_edges_Ybr[k])[3] :
                   last(orientations_and_edges_Ybr[k])[2]
          for k in 0:length(
              orientations_and_edges_Ybr ) ]
        for (node_idx,orientations_and_edges_Ybr) in 
            zip(buses_idx,
           nodes_incident_edges_orientation_and_Ybr_cal)]

    #------------------------------------------

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
    
    #------------------------------------------

    edges_Ybr = edges_Ybr_cal
    
    Ybr_cal_and_edges_orientation =
        (;edges_Ybr_cal,
         edges_orientation )
    
    Ynet_wt_nodes_idx_wt_adjacent_nodes =
        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx)
    
    return (; buses_idx,
            nodes_idx_and_Yshunt,
            nodes_Yshunt,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            edges_Ybr_cal,
            edges_orientation,
            edges_Ybr,
            Ybr_cal_and_edges_orientation,
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            nodes_incident_edges,
            nodes_incident_edges_and_orientation )
    
end


function get_Ynet_and_nodes_idx_wt_adjacent_nodes_idx_by_mpc(
    mpc_bus,
    mpc_branch,
    mpc_baseMVA;
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true )

    (; Ynet,
     nodes_idx_with_adjacent_nodes_idx ) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_mpc(
                 mpc_bus,
                 mpc_branch,
                 mpc_baseMVA;
                 basekV = basekV,
                 use_pu_in_PQ = use_pu_in_PQ,
                 line_data_in_pu = line_data_in_pu ),
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx ))

    return (; 
            Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
end


function get_Ynet_and_nodes_idx_wt_adjacent_nodes_idx_by_mpc2(
    mpc_bus,
    mpc_branch,
    mpc_baseMVA )

    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt_by_mpc(
            mpc_bus, mpc_baseMVA )
    
    edges_orientation =
        get_edges_orientation_by_mpc(
            mpc_branch )

    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )
    
    edges_Ybr_cal =
        get_edges_Ybr_by_mpc(
            mpc_bus, mpc_branch, mpc_baseMVA )

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
                               orientations_and_edges_Ybr[k])) ?
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

    return (; 
            Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
end

#---------------------------------------------------
#---------------------------------------------------
# functions by json
#---------------------------------------------------
#---------------------------------------------------


function get_edges_orientation_by_json(
    edge_data_from_json)

    branches_fbus, branches_tbus =
        get_edges_fbus_tbus_by_json(
            edge_data_from_json)
    
    return get_edges_orientation_by_generic(
        branches_fbus, branches_tbus )
end


function get_sta_pf_PQ_param_by_json(
    plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json;
    baseMVA = 1.0 )

    (;loc_load_exist,            
     gen_nodes_PQ,
     local_load) =
         get_gens_plants_PQ_and_loc_loads_by_json(
             plant_generators_data_from_json )

    (; non_gen_nodes_PQ, ) =
        get_non_gen_plants_PQ_by_json(
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    #---------------------------------


    P_gens = first.(gen_nodes_PQ) ./baseMVA

    Q_gens = second.(gen_nodes_PQ) ./baseMVA
    
    P_non_gens = first.(non_gen_nodes_PQ) ./baseMVA

    Q_non_gens = second.(non_gen_nodes_PQ) ./baseMVA

    if loc_load_exist == true

        P_g_loc_load = first.(local_load) ./baseMVA

        Q_g_loc_load = second.(local_load) ./baseMVA
    else
        
        P_g_loc_load = []

        Q_g_loc_load = []
    end
    
    return (; P_gens,
            Q_gens,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,
            loc_load_exist)
    
end


function get_pf_PQ_param_by_json(
    plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json;
    baseMVA = 1.0,
    use_pu_in_PQ = true )

    if use_pu_in_PQ == true
        
        return get_sta_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA = baseMVA )
    else
        
        return get_sta_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA = 1.0 )
    end
    
end



function get_gens_vh_slack_θh_para_by_json(
    plant_generators_data_from_json )


    # sequence_order =
    #     (:components_data, :gen)
    
    # selections =
    #     (:vh)

    
    gens_vh =
        [ (NamedTupleTools.select(
            nt_plant_gen,
            (:components_data,))).:components_data.gen.vh
          for nt_plant_gen in
              plant_generators_data_from_json]

    
    gens_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_generators_data_from_json,
                (:idx,) ))

    n2s_gens_idx = OrderedDict{Int64,Int64}( net_idx =>idx
            for (net_idx, idx) in
                zip( gens_nodes_idx,
                     collect(1:length(
                         gens_nodes_idx ))
                     ))

    gen_is_slack =
        get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json;
            sequence_order =
                (:additional_data, ),
            selections =
                (:isa_slack,) )

    slack_gens_nodes_idx =
        [ idx
          for (idx, a_is_slack) in
              zip(gens_nodes_idx,
                  gen_is_slack )
              if a_is_slack.isa_slack == true ]

    non_slack_gens_nodes_idx =
        sort(collect(
            setdiff(Set(gens_nodes_idx),
                    Set(slack_gens_nodes_idx))))
    

    # -----------------------------------------------

    slack_gens_vh =
        [ gens_vh[n2s_gens_idx[idx]]
          for idx in gens_nodes_idx if idx ∈
              slack_gens_nodes_idx ]
    
    slack_gens_θh = zeros(length(slack_gens_vh))
    
    non_slack_gens_vh =
        [ gens_vh[n2s_gens_idx[a_non_slack_idx]]
          for a_non_slack_idx in
              non_slack_gens_nodes_idx  ]
    
    # -----------------------------------------------

    return (;
            slack_gens_vh,
            slack_gens_θh,
            gens_vh,
            non_slack_gens_vh )    

end


function get_sta_pf_vars_and_paras_idx_by_json(
    plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json )

    (;slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     nodes_with_demands_idx)  =
         get_net_nodes_type_idxs_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json)
    
    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx
    
    #------------------------------------------        

    nodes_size =
        length(gens_nodes_idx) +
        length(non_gens_nodes_idx)
    
    #------------------------------------------

    nodes_types_idxs =
        (;slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,         
         gens_with_loc_load_idx,
         all_nodes_idx)
    
    #------------------------------------------        

    (;
     n2s_slack_gens_idx,         
     n2s_non_slack_gens_idx,     
     n2s_gens_idx,               
     n2s_non_gens_idx,           
     n2s_load_idx,               
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,      
     n2s_all_nodes_idx ) =
         get_dict_net_streamlined_idx_by_nodes_type_idxs(
             get_net_nodes_type_idxs_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json))


    n2s_idxs =
        (;n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx )


    #-------------------------------

    transformed_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in slack_gens_nodes_idx ]

    transformed_non_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_slack_gens_nodes_idx ]

    transformed_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_nodes_idx ]

    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_gens_nodes_idx ]

    transformed_gens_with_loc_load_idx = [
        n2s_all_nodes_idx[idx]
        for idx in gens_with_loc_load_idx ]
        
    transformed_non_slack_gens_and_non_gens_idx = [
        n2s_all_nodes_idx[idx]
        for idx in non_slack_gens_and_non_gens_idx ]
    
    transformed_all_nodes_idx = [
        n2s_all_nodes_idx[idx]
        for idx in all_nodes_idx ]
    
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
    
    # ------------------------------------------------

    # non_gens_vh_idx = setdiff(vh_IDX, gens_idx)

    # non_slack_gens_θh_idx =
    #     θh_IDX[ non_slack_gens_nodes_idx ]

    # non_gens_θh_idx =
    #     θh_IDX[ non_gens_nodes_idx ]

    # red_θh_idx =
    #     setdiff(θh_IDX, θh_IDX[ slack_bus_idx ])


    non_gens_vh_idx =
        setdiff(vh_IDX,
                transformed_gens_nodes_idx )

    non_slack_gens_θh_idx =
        θh_IDX[ transformed_non_slack_gens_nodes_idx ]

    non_gens_θh_idx =
        θh_IDX[ transformed_non_gens_nodes_idx ]

    red_θh_idx =
        setdiff(θh_IDX,
                θh_IDX[transformed_slack_gens_nodes_idx ])
    
    # -----------------------------------------------
    
    red_vh_θh_idx =
        [  non_gens_vh_idx...;
           non_slack_gens_θh_idx...;
           non_gens_θh_idx... ]

    # -----------------------------------------------

    red_var_comps_idxs =
        (; non_gens_vh_idx,
         non_slack_gens_θh_idx,
         non_gens_θh_idx )

    # ------------------------------------------------

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
    
    red_non_gens_θh_Idxs  =
        red_vh_θh_IDX[3]

    # # This is based on erronous assumtion that gen idxs
    # # are first
    # red_θh_Idxs =
    #     first(red_non_slack_gens_θh_Idxs):last(
    #         red_non_gens_θh_Idxs)

    red_θh_non_slack_gens_θh_and_non_gens_θh_Idxs =
        [red_non_slack_gens_θh_Idxs;
         red_non_gens_θh_Idxs]
    
    red_θh_Idxs =
        minimum(red_θh_non_slack_gens_θh_and_non_gens_θh_Idxs
            ):maximum(
                red_θh_non_slack_gens_θh_and_non_gens_θh_Idxs)
    
     # ----------------------------------------------- 
    
    #  red_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
    #          idx => Idx for (idx, Idx) in
    #              zip( non_slack_gens_and_non_gens_idx,
    #                   red_θh_Idxs ) )

    # red_dict_θh_idx2Idx_in_Idx =
    #     OrderedDict{Int64, Int64}(
    #          idx => Idx for (idx, Idx) in
    #              zip( non_slack_gens_and_non_gens_idx,
    #                   1:length(red_θh_Idxs) ) )

    #  red_non_slack_gens_θh_idx2Idx =
    #      [ red_dict_θh_idx2Idx[idx]
    #       for idx in
    #           non_slack_gens_nodes_idx ]

    #  red_non_slack_gens_θh_idx2Idx_in_Idx =
    #      [ red_dict_θh_idx2Idx_in_Idx[idx]
    #        for idx in
    #           non_slack_gens_nodes_idx ]

    #  red_non_gens_θh_idx2Idx =
    #      [ red_dict_θh_idx2Idx[idx]
    #       for idx in
    #           non_gens_nodes_idx ]

    #  red_non_gens_θh_idx2Idx_in_Idx =
    #      [ red_dict_θh_idx2Idx_in_Idx[idx]
    #        for idx in
    #           non_gens_nodes_idx ]



    
     red_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
             idx => Idx for (idx, Idx) in
                 zip( non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

    red_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
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

    
    # ----------------------------------------------- 

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
    
    # ------------------------------------------------

    dim_P_gens = dim_Q_gens =
        length(gens_nodes_idx)
    
    dim_P_non_gens = dim_Q_non_gens =
        length(non_gens_nodes_idx)

    dim_P_g_loc_load = dim_Q_g_loc_load =
        loc_load_exist == true ?
        length(gens_nodes_with_loc_loads_idx ) : nothing
   
    # ----------------------------------------------------

    if loc_load_exist == true

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

        PQ_sta_para_Idxs =
            (; P_gens_sta_para_Idxs,
             Q_gens_sta_para_Idxs,
             P_non_gens_sta_para_Idxs,
             Q_non_gens_sta_para_Idxs,
             P_g_loc_load_sta_para_Idxs,
             Q_g_loc_load_sta_para_Idxs )

    else

        dim_pf_sta_PQ_para  =
                [ dim_P_gens,
                  dim_Q_gens,
                  dim_P_non_gens,
                  dim_Q_non_gens ]

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

        P_g_loc_load_sta_para_Idxs = nothing

        Q_g_loc_load_sta_para_Idxs = nothing

        PQ_sta_para_Idxs =
            (; P_gens_sta_para_Idxs,
             Q_gens_sta_para_Idxs,
             P_non_gens_sta_para_Idxs,
             Q_non_gens_sta_para_Idxs,
             P_g_loc_load_sta_para_Idxs,
             Q_g_loc_load_sta_para_Idxs )


    end
    
    # --------------------------------------------------

    return (;
            red_types_Idxs_etc,
            PQ_sta_para_Idxs,
            nodes_types_idxs,
            n2s_idxs )    
end



function get_transmission_network_parameters_by_json(
    plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json,
    edge_data_from_json,
    shunt_data_from_json;
    baseMVA = 1.0,
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true )

    # -------------------------------------

    if use_pu_in_PQ == true

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
                    shunt_data_from_json)...;
                baseMVA = baseMVA )        
    else

        nodes_idx_and_Yshunt =
            get_nodes_idx_and_Yshunt_by_generic(
                get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
                    shunt_data_from_json)...;
                baseMVA = 1.0 )        
    end

    nodes_Yshunt = second.( nodes_idx_and_Yshunt )
    
    # -------------------------------------

    if line_data_in_pu == true && use_pu_in_PQ == true
        
        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx ) =
         get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
                get_edges_ftbus_and_generic_data_by_json(
                     edge_data_from_json )...,
                get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
                     shunt_data_from_json )...;
                 baseMVA = 1.0,
                 basekV = 1.0,
                 baseShunt = baseMVA )

        # -----------------------------------------------

        (;edges_Ybr_cal,
         edges_orientation ) =
         get_edges_Ybr_cal_and_edges_orientation_by_generic(
                 get_edges_fbus_tbus_by_json(
                     edge_data_from_json)...,
                 get_edges_generic_data_by_json(
                     edge_data_from_json )...,
                 get_nodes_shunts_Gs_and_Bs_by_json(
                     shunt_data_from_json)...;
                 baseMVA = 1.0,
                 basekV = 1.0)
        
    elseif line_data_in_pu == false && use_pu_in_PQ == true
        
        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx ) =
         get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
                 get_edges_ftbus_and_generic_data_by_json(
                     edge_data_from_json )...,
                 get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
                     shunt_data_from_json)...;
                 baseMVA = baseMVA,
                 basekV = basekV,
                 baseShunt = baseMVA  )

        # -----------------------------------------------

        (;edges_Ybr_cal,
         edges_orientation ) =
         get_edges_Ybr_cal_and_edges_orientation_by_generic(
                 get_edges_fbus_tbus_by_json(
                     edge_data_from_json)...,
                 get_edges_generic_data_by_json(
                     edge_data_from_json )...,
                 get_nodes_shunts_Gs_and_Bs_by_json(
                     shunt_data_from_json)...;
                 baseMVA =  baseMVA,
                 basekV = basekV)
        
    elseif line_data_in_pu == false && use_pu_in_PQ == false
        
        (;Ynet,
         nodes_idx_with_adjacent_nodes_idx ) =
         get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
                 get_edges_ftbus_and_generic_data_by_json(
                     edge_data_from_json )...,
                 get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
                     shunt_data_from_json)...;
                 baseMVA = 1.0,
                 basekV = 1.0,
                 baseShunt = 1.0  )

        # -------------------------------------

        (;edges_Ybr_cal,
         edges_orientation) =
         get_edges_Ybr_cal_and_edges_orientation_by_generic(
                 get_edges_fbus_tbus_by_json(
                     edge_data_from_json)...,
                 get_edges_generic_data_by_json(
                     edge_data_from_json )...,
                 get_nodes_shunts_Gs_and_Bs_by_json(
                     shunt_data_from_json)...;
                 baseMVA =  1.0,
                 basekV = 1.0)
        
    else

        throw("line_data_in_pu = $(line_data_in_pu) and use_pu_in_PQ = $(use_pu_in_PQ) option is not consistent. ")
        
    end
    
    # -------------------------------------

    edges_Ybr = edges_Ybr_cal

    Ybr_cal_and_edges_orientation =
        (;edges_Ybr_cal,
         edges_orientation )    

    Ynet_wt_nodes_idx_wt_adjacent_nodes =
        (;Ynet,
        nodes_idx_with_adjacent_nodes_idx)


    return (;nodes_idx_and_Yshunt,
            nodes_Yshunt,
            Ynet,
            nodes_idx_with_adjacent_nodes_idx,
            edges_Ybr_cal,
            edges_orientation,
            edges_Ybr,
            Ybr_cal_and_edges_orientation,
            Ynet_wt_nodes_idx_wt_adjacent_nodes) 

    
end


function get_pf_sta_ΔPQ_mismatch_parameters_by_json(
    plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json,
    edge_data_from_json,
    shunt_data_from_json;
    baseMVA = 1.0,
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true  )
    
    #-----------------------------------------------
    
    (;
     red_types_Idxs_etc,
     PQ_sta_para_Idxs,
     nodes_types_idxs,
     n2s_idxs ) =
         get_sta_pf_vars_and_paras_idx_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json )

    #-----------------------------------------------

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

    
    (;
     P_gens_sta_para_Idxs,
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

    # -----------------------------------------------

    (; slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh ) =
         get_gens_vh_slack_θh_para_by_json(
             plant_generators_data_from_json )
        
    # ------------------------------------------------

    (;
     P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist ) =
         get_pf_PQ_param_by_json(
             plant_generators_data_from_json,
             plant_loads_data_from_json,
             plant_transmission_data_from_json;
             baseMVA = baseMVA,
             use_pu_in_PQ = use_pu_in_PQ )

    (;
     Ynet,
     nodes_idx_with_adjacent_nodes_idx,
     edges_Ybr_cal,
     edges_orientation) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json,
                 edge_data_from_json,
                 shunt_data_from_json;
                 baseMVA = baseMVA,
                 basekV  = basekV,
                 use_pu_in_PQ = use_pu_in_PQ,
                 line_data_in_pu = line_data_in_pu ),
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx,
              :edges_Ybr_cal,
              :edges_orientation) )

    # -------------------------------------------------

    net_para = (;
                Ynet,
                nodes_idx_with_adjacent_nodes_idx,
                edges_Ybr_cal,
                edges_orientation)
    
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
        (;n2s_slack_gens_idx,
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

    return (;
            pf_kw_para,
            pf_PQ_param,
            red_types_Idxs_etc,
            net_para  )
    
end


#---------------------------------------------------
#---------------------------------------------------
# components data related functions 
#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
# static plant data
#---------------------------------------------------

function get_static_gens_instance_data_by_mpc(
    plants_idx_and_type,
    gens_nodes_static_tup_data)

    @assert first.(plants_idx_and_type) ==
        first.(gens_nodes_static_tup_data)

    return [
        (a_gen_type,
         a_gen_static_tup_data )

        for (a_gen_type,
             a_gen_static_tup_data) in
            zip(second.(plants_idx_and_type),
                second.(gens_nodes_static_tup_data))]
    
end

function get_a_static_gen_plant_wt_loc_load_by_mpc(
    idx,
    plant_type,
    gen_type,
    isa_slack,
    gen_data,
    loc_load_instance_data)

    return (idx = idx,
         plant_type = plant_type,
         additional_data = 
             (isa_slack = isa_slack, ),
         components_type =
             (gen = gen_type,
              loc_load =
                  first(loc_load_instance_data)),
         components_data =
             (gen = second(gen_data),
              loc_load = second(
                  loc_load_instance_data))) 
end


function get_a_static_gen_plant_wt_no_loc_load_by_mpc(
    idx,
    plant_type,
    gen_type,
    isa_slack,
    gen_data)

    return (idx = idx,
         plant_type = plant_type,
         additional_data = 
             (isa_slack = isa_slack, ),
         components_type =
             ( gen = gen_type, ),
         components_data =
             ( gen = second(gen_data), ) ) 
end


function get_static_gens_plant_instances_data_by_mpc(
    dict_plants_gen_sym_type,
    dict_gen_sym_type,
    
    dyn_gens,
    dyn_plants,

    mpc_bus,
    mpc_gen;
    mpc_baseMVA =
        1.0,

    by_components =
        false)

    #------------------------------------------

    loc_load_exist =
        loc_load_exist_bool_by_mpc(
            mpc_bus )

    loc_loads =
        get_loc_loads_idx_and_locP_locQ_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA )

    
    gens_nodes_idx =
        get_gens_nodes_idx_by_mpc(
            mpc_bus)
    
    gens_with_loc_load_idx =
        get_gens_nodes_with_loc_loads_idx_by_mpc(
            mpc_bus)
    
    #--------------------------------------
    
    n2s_gens_idx =
        get_n2s_any(gens_nodes_idx)

    n2s_gens_with_loc_load_idxs =
        loc_load_exist == true ?  get_n2s_any(
            gens_with_loc_load_idx) : get_n2s_any(
                gens_with_loc_load_idx;
                nothing_bool = true)
    
    #------------------------------------------    
    
    plants_idx =
        dyn_plants.bus

    sym_plants_types =
        dyn_plants.Plant_type
    
    sym_gens_types =
        dyn_plants.Gen
    
    bool_isa_slack =
        dyn_plants.isa_slack
    
    plants_types =
        [dict_plants_gen_sym_type[sym_plant_type]
         for sym_plant_type in
             sym_plants_types]

    gens_idx_and_type =
        get_static_gens_idx_and_type(
            plants_idx,
            sym_gens_types)

    plants_idx_and_type =
        get_static_gens_idx_and_type(
            plants_idx,
            sym_plants_types)

    gens_nodes_static_tup_data =
        get_gen_nodes_static_tup_data_by_mpc(
            mpc_bus, mpc_gen;
            mpc_baseMVA=
                mpc_baseMVA )

    #------------------------------------------

    # gens_idx_and_type can also be used
    # as the first arguments
    
    static_gens_instance_data =
        get_static_gens_instance_data_by_mpc(
            plants_idx_and_type,
            gens_nodes_static_tup_data)

    loc_load_instances_data =
        loc_load_exist != false ? [
            (loc_Load_t1, a_loc_load ) 
            for a_loc_load in loc_loads  ] : []

    # @show loc_load_instances_data
    # @show gens_with_loc_load_idx
    
    if by_components == false
        return [ plant_idx ∈ gens_with_loc_load_idx ?
            (plant_idx,(
                plant_type = plant_type,
                gen_type_and_data = gen_data,
                isa_slack = isa_slack,
                loc_load_type_and_data=
                    loc_load_instances_data[
                        n2s_gens_with_loc_load_idxs[
                            plant_idx]] )) :
                    (plant_idx,
                     ( plant_type = plant_type,
                       gen_type_and_data = gen_data,
                       isa_slack = isa_slack) )
                 for (plant_idx,
                      plant_type,
                      gen_type,
                      gen_data,
                      isa_slack ) in
                     zip(gens_nodes_idx,
                         plants_types,
                         sym_gens_types,
                         static_gens_instance_data,
                         bool_isa_slack ) ]
    else
    
        return [ idx ∈ gens_with_loc_load_idx  ?
               get_a_static_gen_plant_wt_loc_load_by_mpc(
                   idx,
                   plant_type,
                   gen_type,
                   isa_slack,
                   gen_data,
                   loc_load_instances_data[
                       n2s_gens_with_loc_load_idxs[idx]]) :
         get_a_static_gen_plant_wt_no_loc_load_by_mpc(
             idx,
             plant_type,
             gen_type,
             isa_slack,
             gen_data)

                 for (idx,
                      plant_type,
                      gen_type,
                       isa_slack,
                      gen_data ) in
                     zip(gens_nodes_idx,
                         plants_types,
                         sym_gens_types,
                         bool_isa_slack,
                         static_gens_instance_data) ]
    end
        
end

#---------------------------------------------------
#---------------------------------------------------

# function get_multi_gens_plant_instances_data_by_mpc(
#     dict_plants_gen_sym_type,
#     dict_gen_sym_type,
#     dict_gens_dyn_nt_params,

#     dict_gov_sym_type,
#     dict_gov_nt_params,

#     dict_avr_sym_type,
#     dict_avr_nt_params,
    
#     dyn_gens,
#     dyn_plants,

#     mpc_bus,
#     mpc_gen;
#     mpc_baseMVA=1.0,

#     p_order = 1.0,
#     v_ref = 1.0,
#     ω_ref = ωs,
#     by_components = false )

#     #------------------------------------------

#     (;multi_gens_idx,
#      multi_gen_bool,
#      multi_gens_nodes) =
#          NamedTupleTools.select(
#              get_multi_gens_idx_wt_multi_gen_bool(
#                  mpc_gen),
#              (:multi_gens_idx,
#               :multi_gen_bool,
#               :multi_gens_nodes) )

#     #------------------------------------------

#     loc_load_exist =
#         loc_load_exist_bool_by_mpc(
#             mpc_bus )

#     loc_loads =
#         get_loc_loads_idx_and_locP_locQ_data_by_mpc(
#             mpc_bus;
#             mpc_baseMVA =
#                 mpc_baseMVA )

    
#     gens_nodes_idx =
#         get_gens_nodes_idx_by_mpc(
#             mpc_bus)
    
#     gens_with_loc_load_idx =
#         get_gens_nodes_with_loc_loads_idx_by_mpc(
#             mpc_bus)

#     #--------------------------------------
    
#     n2s_gens_idx =
#         get_n2s_any(gens_nodes_idx)

#     n2s_gens_with_loc_load_idxs =
#         loc_load_exist == true ?  get_n2s_any(
#             gens_with_loc_load_idx) : get_n2s_any(
#                 gens_with_loc_load_idx;
#                 nothing_bool = true)


#     #------------------------------------------    
#     #------------------------------------------    
    
#     plants_idx =
#         dyn_plants.bus

#     sym_plants_types =
#         dyn_plants.Plant_type
    
#     sym_gens_types =
#         dyn_plants.Gen
    
#     bool_isa_slack =
#         dyn_plants.isa_slack
    
#     sym_govs_nt_params =
#         dyn_plants.Gov
    
#     sym_avrs_nt_params =
#         dyn_plants.Exc


#     plants_types =
#         [dict_plants_gen_sym_type[sym_plant_type]
#          for sym_plant_type in
#              sym_plants_types]

#     gens_type_and_dym_data =
#         get_gens_type_and_dym_data(
#             dict_gen_sym_type,
#             dict_gens_dyn_nt_params,
#             dyn_gens)

#     gens_nodes_static_tup_data =
#         get_gen_nodes_static_tup_data_by_mpc(
#             mpc_bus,
#             mpc_gen;
#             mpc_baseMVA =
#                 mpc_baseMVA)

#     #------------------------------------------

#     gens_instances_data =
#         get_gens_instance_data_by_mpc(
#             gens_type_and_dym_data,
#             gens_nodes_static_tup_data)

#     govs_instances_data =
#         get_govs_instances_data_by_mpc(
#             sym_govs_nt_params,
#             dict_gov_nt_params,
#             dict_gov_sym_type)

#     avrs_instances_data =
#         get_avrs_instances_data_by_mpc(
#             sym_avrs_nt_params,
#             dict_avr_nt_params,
#             dict_avr_sym_type)

#     loc_load_instances_data =
#         loc_load_exist != false ? [
#             (loc_Load_t1, a_loc_load ) 
#             for a_loc_load in loc_loads  ] : []

#     if multi_gen_bool == true

#         if by_components == false
#             return [
#                 get_node_idx_multi(plant_idx) ∈ gens_with_loc_load_idx ?
#                 (plant_idx,(
#                     # plant_idx = plant_idx,
#                     plant_type =
#                         plant_type,
#                     gen_type_and_data =
#                         gen_data,
#                     isa_slack =
#                         isa_slack,
#                     gov_type_and_data =
#                         gov_data,
#                     avr_type_and_data =
#                         avr_data,
#                     loc_load_type_and_data=
#                         loc_load_instances_data[
#                         n2s_gens_with_loc_load_idxs[
#                             get_node_idx_multi(plant_idx) ] ],
#                     p_order  = p_order,
#                     v_ref = v_ref,
#                     ω_ref = ω_ref )) :
#                         (plant_idx,
#                          (# plant_idx = plant_idx,
#                            plant_type =
#                                plant_type,
#                            gen_type_and_data =
#                                gen_data,
#                            isa_slack =
#                                isa_slack,
#                            gov_type_and_data =
#                                gov_data,
#                            avr_type_and_data =
#                                avr_data,
#                            p_order =
#                                p_order,
#                            v_ref =
#                                v_ref,
#                            ω_ref =
#                                ω_ref ) )
#                 for (plant_idx,
#                      plant_type,
#                      gen_data,
#                      isa_slack,
#                      gov_data,
#                      avr_data) in zip(multi_gens_idx,
#                                       plants_types,
#                                       gens_instances_data,
#                                       bool_isa_slack,
#                                       govs_instances_data,
#                                       avrs_instances_data) ]
#         else

#             return [ (get_node_idx_multi(idx) ∈ gens_with_loc_load_idx && first( gov_data != nothing) ) ? get_a_gen_plant_wt_gov_wt_loc_load_by_mpc( idx, plant_type, isa_slack, p_order, v_ref, ω_ref, gen_data, gov_data, avr_data, loc_load_instances_data, n2s_gens_with_loc_load_idxs) :
#                            ( get_node_idx_multi(idx) ∈ gens_with_loc_load_idx && first( gov_data == nothing) ) ? get_a_gen_plant_no_gov_wt_loc_load_by_mpc( idx, plant_type, isa_slack, p_order, v_ref, ω_ref, gen_data, gov_data, avr_data, loc_load_instances_data, n2s_gens_with_loc_load_idxs) :
#                                    ( get_node_idx_multi(idx) ∉ gens_with_loc_load_idx && first( gov_data != nothing) ) ? get_a_gen_plant_wt_gov_no_loc_load_by_mpc(idx, plant_type, isa_slack, p_order, v_ref, ω_ref, gen_data, gov_data, avr_data, loc_load_instances_data, n2s_gens_with_loc_load_idxs) : get_a_gen_plant_no_gov_no_loc_load_by_mpc(idx, plant_type, isa_slack, p_order, v_ref, ω_ref, gen_data, gov_data, avr_data)

#                     for (idx,
#                          plant_type,
#                          gen_data,
#                          isa_slack,
#                          gov_data,
#                          avr_data) in
#                         zip(multi_gens_idx,
#                             plants_types,
#                             gens_instances_data,
#                             bool_isa_slack,
#                             govs_instances_data,
#                             avrs_instances_data) ]


#         end
        
#     else


#         if by_components == false
#             return [ plant_idx ∈ gens_with_loc_load_idx ?
#                 (plant_idx,(
#                     # plant_idx = plant_idx,
#                     plant_type = plant_type,
#                     gen_type_and_data = gen_data,
#                     isa_slack = isa_slack,
#                     gov_type_and_data = gov_data,
#                     avr_type_and_data = avr_data,
#                     loc_load_type_and_data=loc_load_instances_data[
#                         n2s_gens_with_loc_load_idxs[plant_idx]],
#                     p_order  = p_order,
#                     v_ref = v_ref,
#                     ω_ref = ω_ref )) :
#                         (plant_idx,
#                          ( # plant_idx = plant_idx,
#                            plant_type = plant_type,
#                            gen_type_and_data = gen_data,
#                            isa_slack = isa_slack,
#                            gov_type_and_data = gov_data,
#                            avr_type_and_data = avr_data,
#                            p_order = p_order,
#                            v_ref = v_ref,
#                            ω_ref = ω_ref) )
#                      for (plant_idx, plant_type,
#                           gen_data, isa_slack,
#                           gov_data, avr_data) in
#                          zip(gens_nodes_idx,
#                             plants_types,
#                             gens_instances_data,
#                             bool_isa_slack,
#                             govs_instances_data,
#                              avrs_instances_data) ]
#         else

#             return [ ( idx ∈ gens_with_loc_load_idx && first(
#                 gov_data != nothing )) ?
#                    get_a_gen_plant_wt_gov_wt_loc_load_by_mpc(
#                        idx, plant_type, isa_slack,
#                        p_order, v_ref, ω_ref,
#                        gen_data, gov_data, avr_data,
#                        loc_load_instances_data,
#                        n2s_gens_with_loc_load_idxs) :
#                            (idx ∈ gens_with_loc_load_idx &&
#                            first( gov_data == nothing )) ?
#                            get_a_gen_plant_no_gov_wt_loc_load_by_mpc(
#                                idx, plant_type, isa_slack,
#                                p_order, v_ref, ω_ref,
#                                gen_data, gov_data, avr_data,
#                                loc_load_instances_data,
#                                n2s_gens_with_loc_load_idxs) :
#                                    (idx ∉ gens_with_loc_load_idx &&
#                                    first( gov_data != nothing )) ?
#                                    get_a_gen_plant_wt_gov_no_loc_load_by_mpc(
#                                        idx, plant_type,
#                                        isa_slack, p_order,
#                                        v_ref, ω_ref, gen_data,
#                                        gov_data, avr_data,
#                                        loc_load_instances_data,
#                                        n2s_gens_with_loc_load_idxs) :
#                                            get_a_gen_plant_no_gov_no_loc_load_by_mpc(
#                                                idx, plant_type,
#                                                isa_slack, p_order,
#                                                v_ref, ω_ref, gen_data,
#                                                gov_data, avr_data )

#                      for (idx, plant_type, gen_data,
#                           isa_slack, gov_data, avr_data) in
#                          zip(gens_nodes_idx, plants_types,
#                             gens_instances_data, bool_isa_slack,
#                              govs_instances_data,
#                              avrs_instances_data) ]


#         end

        
#     end
    
# end



#---------------------------------------------------
# Gens
#---------------------------------------------------

function get_a_gen_type_and_dym_data(
    idx, n2s_gens_idx,
    dict_gen_sym_type,
    dict_gens_dyn_nt_params,
    dyn_gens)

    # gen_idx =
    #     dyn_gens.bus[n2s_gens_idx[idx]]

    sym_gen_type =
        dyn_gens.sym_gen_type[n2s_gens_idx[idx]]

    sym_gen_dym_para =
        dyn_gens.sym_gen_dynamic_para[n2s_gens_idx[idx]]
    
    return (idx,
             (dict_gen_sym_type[sym_gen_type],
              dict_gens_dyn_nt_params[sym_gen_dym_para]))

end


function get_gens_type_and_dym_data(
    dict_gen_sym_type,
    dict_gens_dyn_nt_params,
    dyn_gens)

    gens_type_and_dym_data = [(idx,
             (dict_gen_sym_type[sym_gen_type],
              dict_gens_dyn_nt_params[sym_gen_dym_para]))
            for (idx, sym_gen_type, sym_gen_dym_para) in
                zip(dyn_gens.bus,
                    dyn_gens.sym_gen_type,
                    dyn_gens.sym_gen_dynamic_para)]
    
     return sort(gens_type_and_dym_data, by = x -> x[1]) 

end



function get_static_gens_idx_and_type(
    plants_idx,
    sym_gens_types)

    static_gens_idx_and_type = [ (idx, sym_gen_type)
            for (idx, sym_gen_type) in
                zip(plants_idx,
                    sym_gens_types)]

    return sort(static_gens_idx_and_type, by = x-> x[1] )

end

#---------------------------------------------------
# Govs
#---------------------------------------------------

function get_govs_instances_data_by_mpc(
    sym_govs_nt_params,
    dict_gov_nt_params,
    dict_gov_sym_type)

    sym_govs_types =
        Symbol.(
            first.(split.(
                String.(sym_govs_nt_params),"__") ) )

    
    govs_types =
        [ sym_gov_type == :nothing ? :nothing :
        dict_gov_sym_type[sym_gov_type]
          for sym_gov_type in
              sym_govs_types  ]

    govs_nt_params =
        [ sym_gov_nt_params == :nothing ? :nothing :
        dict_gov_nt_params[sym_gov_nt_params]
          for sym_gov_nt_params in
              sym_govs_nt_params  ]

    return [
        sym_gov_nt_param == nothing ?
            (nothing, nothing) : (gov_type, gov_nt_params)
            for (gov_type,
                 gov_nt_params,
                 sym_gov_nt_param) in
                zip(govs_types,
                    govs_nt_params,
                    sym_govs_nt_params )]
    
end


#---------------------------------------------------
# Avrs
#---------------------------------------------------

function get_avrs_instances_data_by_mpc(
    sym_avrs_nt_params,
    dict_avr_nt_params,
    dict_avr_sym_type)

    sym_avrs_types =
        Symbol.(
            first.(split.(
                String.(sym_avrs_nt_params),"__") ) )

    avrs_types =
        [ dict_avr_sym_type[sym_avr_type]
          for sym_avr_type in
              sym_avrs_types  ]

    avrs_nt_params =
        [ dict_avr_nt_params[sym_avr_nt_params]
          for sym_avr_nt_params in
              sym_avrs_nt_params  ]

    return [(avr_type, avr_nt_params)
            for (avr_type, avr_nt_params) in
                zip(avrs_types, avrs_nt_params )]
    
end


#---------------------------------------------------
# Dynamic plant data
#---------------------------------------------------

function get_gens_instance_data_by_mpc(
    gens_type_and_dym_data,
    gens_nodes_static_tup_data)

    @assert first.(gens_type_and_dym_data) ==
        first.(gens_nodes_static_tup_data)

    return [(first( a_gen_type_and_dym_data),
            merge( second(a_gen_type_and_dym_data) ,
            a_gen_static_tup_data ) )

        for (a_gen_type_and_dym_data,
             a_gen_static_tup_data) in
            zip(second.(gens_type_and_dym_data),
                second.(gens_nodes_static_tup_data))]
    
end

function get_a_gen_plant_wt_gov_wt_loc_load_by_mpc(
    idx, plant_type, isa_slack,
    p_order, v_ref, ω_ref,
    gen_data, gov_data, avr_data,
    loc_load_instances_data,
    n2s_gens_with_loc_load_idxs)

    return (idx = idx,
         plant_type = plant_type,
         additional_data = 
             (isa_slack = isa_slack,
              p_order  = p_order,
              v_ref = v_ref,
              ω_ref = ω_ref ),
         components_type =
             ( gen = first(gen_data),
               gov = first(gov_data),
               avr = first(avr_data),
               loc_load = first(loc_load_instances_data[
                   n2s_gens_with_loc_load_idxs[idx]]) ),
         components_data =
             ( gen = second(gen_data),
               gov = second(gov_data),
               avr = second(avr_data),
               loc_load = second(loc_load_instances_data[
                   n2s_gens_with_loc_load_idxs[idx]])) ) 
end


function get_a_gen_plant_wt_gov_no_loc_load_by_mpc(
    idx, plant_type, isa_slack,
    p_order, v_ref, ω_ref,
    gen_data, gov_data, avr_data,
    loc_load_instances_data,
    n2s_gens_with_loc_load_idxs)

    return (idx = idx,
         plant_type = plant_type,
         additional_data = 
             (isa_slack = isa_slack,
              p_order  = p_order,
              v_ref = v_ref,
              ω_ref = ω_ref ),
         components_type =
             ( gen = first(gen_data),
               gov = first(gov_data),
               avr = first(avr_data), ),
         components_data =
             ( gen = second(gen_data),
               gov = second(gov_data),
               avr = second(avr_data) ) ) 
end


function get_a_gen_plant_no_gov_wt_loc_load_by_mpc(
    idx, plant_type, isa_slack,
    p_order, v_ref, ω_ref,
    gen_data, gov_data, avr_data,
    loc_load_instances_data,
    n2s_gens_with_loc_load_idxs)

    @assert first(gov_data) == nothing

    return (idx = idx,
         plant_type = plant_type,
         additional_data = 
             (isa_slack = isa_slack,
              v_ref = v_ref,
              ω_ref = ω_ref ),
         components_type =
             ( gen = first(gen_data),
               avr = first(avr_data),
               loc_load = first(loc_load_instances_data[
                   n2s_gens_with_loc_load_idxs[idx]]) ),
         components_data =
             ( gen = second(gen_data),
               avr = second(avr_data),
               loc_load = second(loc_load_instances_data[
                   n2s_gens_with_loc_load_idxs[idx]])) ) 
end


function get_a_gen_plant_no_gov_no_loc_load_by_mpc(
    idx, plant_type, isa_slack,
    p_order, v_ref, ω_ref,
    gen_data, gov_data, avr_data )

    @assert first(gov_data) == nothing

    return (idx = idx,
         plant_type = plant_type,
         additional_data = 
             (isa_slack = isa_slack,
              v_ref = v_ref,
              ω_ref = ω_ref ),
         components_type =
             ( gen = first(gen_data),
               avr = first(avr_data) ),
         components_data =
             ( gen = second(gen_data),
               avr = second(avr_data) )) 
end


function get_gens_plant_instances_data_by_mpc(
    dict_plants_gen_sym_type,
    dict_gen_sym_type,
    dict_gens_dyn_nt_params,

    dict_gov_sym_type,
    dict_gov_nt_params,

    dict_avr_sym_type,
    dict_avr_nt_params,
    
    dyn_gens,
    dyn_plants,

    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0,

    p_order = 1.0,
    v_ref = 1.0,
    ω_ref = ωs,
    by_components = false)

    #------------------------------------------

    loc_load_exist =
        loc_load_exist_bool_by_mpc(
            mpc_bus )

    loc_loads =
        get_loc_loads_idx_and_locP_locQ_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA )

    
    gens_nodes_idx =
        get_gens_nodes_idx_by_mpc(
            mpc_bus ; sorted_bool = false )

    
    gens_with_loc_load_idx =
        get_gens_nodes_with_loc_loads_idx_by_mpc(
            mpc_bus ; sorted_bool = false )
    
    # n2s_gens_idx =
    #     get_a_n2s_net_group(gens_nodes_idx)

    # n2s_gens_with_loc_load_idxs =
    #     get_a_n2s_net_group(
    #         gens_with_loc_load_idx;
    #         loc_load_exist =
    #             loc_load_exist)

    #--------------------------------------

    
    n2s_gens_idx =
        get_n2s_any(gens_nodes_idx)

    n2s_gens_with_loc_load_idxs =
        loc_load_exist == true ?  get_n2s_any(
            gens_with_loc_load_idx) : get_n2s_any(
                gens_with_loc_load_idx;
                nothing_bool = true)
    
    #------------------------------------------    
    
    plants_idx =
        dyn_plants.bus

    sym_plants_types =
        dyn_plants.Plant_type
    
    sym_gens_types =
        dyn_plants.Gen
    
    bool_isa_slack =
        dyn_plants.isa_slack
    
    sym_govs_nt_params =
        dyn_plants.Gov
    
    sym_avrs_nt_params =
        dyn_plants.Exc

    
    # gens_plants_type_and_components_type =
    #     get_plants_type_and_components_type(
    #         dyn_plants)

    
    # sym_govs_types =
    #     Symbol.(
    #         first.(split.(
    #         String.(sym_govs_nt_params),"__") ) )

    # sym_avrs_types =
    #     Symbol.(
    #         first.(split.(
    #             String.(sym_avrs_nt_params),"__") ) )

    plants_types =
        [dict_plants_gen_sym_type[sym_plant_type]
         for sym_plant_type in
             sym_plants_types]

    # nodes_static_tup_data =
    #     get_nodes_static_tup_data_by_mpc(
    #         mpc_bus, mpc_gen;
    #         mpc_baseMVA=mpc_baseMVA)
    
    # nodes_dict_static_data =
    #     get_nodes_dict_static_data_by_mpc(
    #         mpc_bus,
    #         mpc_gen;
    #         mpc_baseMVA=mpc_baseMVA)

    gens_type_and_dym_data =
        get_gens_type_and_dym_data(
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,
            dyn_gens)

    gens_nodes_static_tup_data =
        get_gen_nodes_static_tup_data_by_mpc(
            mpc_bus, mpc_gen;
            mpc_baseMVA=
                mpc_baseMVA)

    #------------------------------------------

    gens_instances_data =
        get_gens_instance_data_by_mpc(
            gens_type_and_dym_data,
            gens_nodes_static_tup_data)

    govs_instances_data =
        get_govs_instances_data_by_mpc(
            sym_govs_nt_params,
            dict_gov_nt_params,
            dict_gov_sym_type)

    avrs_instances_data =
        get_avrs_instances_data_by_mpc(
            sym_avrs_nt_params,
            dict_avr_nt_params,
            dict_avr_sym_type)

    loc_load_instances_data =
        loc_load_exist != false ? [
            (loc_Load_t1, a_loc_load ) 
            for a_loc_load in loc_loads  ] : []

    # @show loc_load_instances_data
    # @show gens_with_loc_load_idx
    
    if by_components == false
        return [ plant_idx ∈ gens_with_loc_load_idx ?
            (plant_idx,(
                # plant_idx = plant_idx,
                plant_type = plant_type,
                gen_type_and_data = gen_data,
                isa_slack = isa_slack,
                gov_type_and_data = gov_data,
                avr_type_and_data = avr_data,
                loc_load_type_and_data=
                    loc_load_instances_data[
                        n2s_gens_with_loc_load_idxs[
                            plant_idx]],
                p_order  = p_order,
                v_ref = v_ref,
                ω_ref = ω_ref )) :
                    (plant_idx,
                     ( # plant_idx = plant_idx,
                       plant_type = plant_type,
                       gen_type_and_data = gen_data,
                       isa_slack = isa_slack,
                       gov_type_and_data = gov_data,
                       avr_type_and_data = avr_data,
                       p_order = p_order,
                       v_ref = v_ref,
                       ω_ref = ω_ref) )
                 for (plant_idx, plant_type,
                      gen_data, isa_slack,
                      gov_data, avr_data) in
                     zip(gens_nodes_idx,
                        plants_types,
                        gens_instances_data,
                        bool_isa_slack,
                        govs_instances_data,
                         avrs_instances_data) ]
    else
    
        return [( idx ∈ gens_with_loc_load_idx && first(
            gov_data != nothing )) ?
               get_a_gen_plant_wt_gov_wt_loc_load_by_mpc(
                   idx, plant_type, isa_slack,
                   p_order, v_ref, ω_ref,
                   gen_data, gov_data, avr_data,
                   loc_load_instances_data,
                   n2s_gens_with_loc_load_idxs) :
                       (idx ∈ gens_with_loc_load_idx &&
                       first( gov_data == nothing )) ?
                       get_a_gen_plant_no_gov_wt_loc_load_by_mpc(
                           idx, plant_type, isa_slack,
                           p_order, v_ref, ω_ref,
                           gen_data, gov_data, avr_data,
                           loc_load_instances_data,
                           n2s_gens_with_loc_load_idxs) :
                               (idx ∉ gens_with_loc_load_idx &&
                               first( gov_data != nothing )) ?
                               get_a_gen_plant_wt_gov_no_loc_load_by_mpc(
                                   idx, plant_type,
                                   isa_slack, p_order,
                                   v_ref, ω_ref, gen_data,
                                   gov_data, avr_data,
                                   loc_load_instances_data,
                                   n2s_gens_with_loc_load_idxs) :
                                       get_a_gen_plant_no_gov_no_loc_load_by_mpc(
                                           idx, plant_type,
                                           isa_slack, p_order,
                                           v_ref, ω_ref, gen_data,
                                           gov_data, avr_data )

                 for (idx, plant_type, gen_data,
                      isa_slack, gov_data, avr_data) in
                     zip(gens_nodes_idx,
                         plants_types,
                         gens_instances_data,
                         bool_isa_slack,
                         govs_instances_data,
                         avrs_instances_data) ]

        
    end
    
    
end


function get_multi_gens_plant_instances_data_by_mpc(
    dict_plants_gen_sym_type,
    dict_gen_sym_type,
    dict_gens_dyn_nt_params,

    dict_gov_sym_type,
    dict_gov_nt_params,

    dict_avr_sym_type,
    dict_avr_nt_params,
    
    dyn_gens,
    dyn_plants,

    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0,

    p_order = 1.0,
    v_ref = 1.0,
    ω_ref = ωs,
    by_components = false )

    #------------------------------------------

    (;multi_gens_idx,
     multi_gen_bool,
     multi_gens_nodes) =
         NamedTupleTools.select(
             get_multi_gens_idx_wt_multi_gen_bool(
                 mpc_gen),
             (:multi_gens_idx,
              :multi_gen_bool,
              :multi_gens_nodes) )

    #------------------------------------------

    loc_load_exist =
        loc_load_exist_bool_by_mpc(
            mpc_bus )

    loc_loads =
        get_loc_loads_idx_and_locP_locQ_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA )

    
    gens_nodes_idx =
        get_gens_nodes_idx_by_mpc(
            mpc_bus)
    
    gens_with_loc_load_idx =
        get_gens_nodes_with_loc_loads_idx_by_mpc(
            mpc_bus)

    #--------------------------------------
    
    n2s_gens_idx =
        get_n2s_any(gens_nodes_idx)

    n2s_gens_with_loc_load_idxs =
        loc_load_exist == true ?  get_n2s_any(
            gens_with_loc_load_idx) : get_n2s_any(
                gens_with_loc_load_idx;
                nothing_bool = true)


    #------------------------------------------    
    #------------------------------------------    
    
    plants_idx =
        dyn_plants.bus

    sym_plants_types =
        dyn_plants.Plant_type
    
    sym_gens_types =
        dyn_plants.Gen
    
    bool_isa_slack =
        dyn_plants.isa_slack
    
    sym_govs_nt_params =
        dyn_plants.Gov
    
    sym_avrs_nt_params =
        dyn_plants.Exc


    plants_types =
        [dict_plants_gen_sym_type[sym_plant_type]
         for sym_plant_type in
             sym_plants_types]

    gens_type_and_dym_data =
        get_gens_type_and_dym_data(
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,
            dyn_gens)

    gens_nodes_static_tup_data =
        get_gen_nodes_static_tup_data_by_mpc(
            mpc_bus,
            mpc_gen;
            mpc_baseMVA =
                mpc_baseMVA)

    #------------------------------------------

    gens_instances_data =
        get_gens_instance_data_by_mpc(
            gens_type_and_dym_data,
            gens_nodes_static_tup_data)

    govs_instances_data =
        get_govs_instances_data_by_mpc(
            sym_govs_nt_params,
            dict_gov_nt_params,
            dict_gov_sym_type)

    avrs_instances_data =
        get_avrs_instances_data_by_mpc(
            sym_avrs_nt_params,
            dict_avr_nt_params,
            dict_avr_sym_type)

    loc_load_instances_data =
        loc_load_exist != false ? [
            (loc_Load_t1, a_loc_load ) 
            for a_loc_load in loc_loads  ] : []

    if multi_gen_bool == true

        if by_components == false
            return [
                get_node_idx_multi(plant_idx) ∈ gens_with_loc_load_idx ?
                (plant_idx,(
                    # plant_idx = plant_idx,
                    plant_type =
                        plant_type,
                    gen_type_and_data =
                        gen_data,
                    isa_slack =
                        isa_slack,
                    gov_type_and_data =
                        gov_data,
                    avr_type_and_data =
                        avr_data,
                    loc_load_type_and_data=
                        loc_load_instances_data[
                        n2s_gens_with_loc_load_idxs[
                            get_node_idx_multi(plant_idx) ] ],
                    p_order  = p_order,
                    v_ref = v_ref,
                    ω_ref = ω_ref )) :
                        (plant_idx,
                         (# plant_idx = plant_idx,
                           plant_type =
                               plant_type,
                           gen_type_and_data =
                               gen_data,
                           isa_slack =
                               isa_slack,
                           gov_type_and_data =
                               gov_data,
                           avr_type_and_data =
                               avr_data,
                           p_order =
                               p_order,
                           v_ref =
                               v_ref,
                           ω_ref =
                               ω_ref ) )
                for (plant_idx,
                     plant_type,
                     gen_data,
                     isa_slack,
                     gov_data,
                     avr_data) in zip(multi_gens_idx,
                                      plants_types,
                                      gens_instances_data,
                                      bool_isa_slack,
                                      govs_instances_data,
                                      avrs_instances_data) ]
        else

            return [ (get_node_idx_multi(idx) ∈ gens_with_loc_load_idx && first( gov_data != nothing) ) ? get_a_gen_plant_wt_gov_wt_loc_load_by_mpc( idx, plant_type, isa_slack, p_order, v_ref, ω_ref, gen_data, gov_data, avr_data, loc_load_instances_data, n2s_gens_with_loc_load_idxs) :
                           ( get_node_idx_multi(idx) ∈ gens_with_loc_load_idx && first( gov_data == nothing) ) ? get_a_gen_plant_no_gov_wt_loc_load_by_mpc( idx, plant_type, isa_slack, p_order, v_ref, ω_ref, gen_data, gov_data, avr_data, loc_load_instances_data, n2s_gens_with_loc_load_idxs) :
                                   ( get_node_idx_multi(idx) ∉ gens_with_loc_load_idx && first( gov_data != nothing) ) ? get_a_gen_plant_wt_gov_no_loc_load_by_mpc(idx, plant_type, isa_slack, p_order, v_ref, ω_ref, gen_data, gov_data, avr_data, loc_load_instances_data, n2s_gens_with_loc_load_idxs) : get_a_gen_plant_no_gov_no_loc_load_by_mpc(idx, plant_type, isa_slack, p_order, v_ref, ω_ref, gen_data, gov_data, avr_data)

                    for (idx,
                         plant_type,
                         gen_data,
                         isa_slack,
                         gov_data,
                         avr_data) in
                        zip(multi_gens_idx,
                            plants_types,
                            gens_instances_data,
                            bool_isa_slack,
                            govs_instances_data,
                            avrs_instances_data) ]


        end
        
    else


        if by_components == false
            return [ plant_idx ∈ gens_with_loc_load_idx ?
                (plant_idx,(
                    # plant_idx = plant_idx,
                    plant_type = plant_type,
                    gen_type_and_data = gen_data,
                    isa_slack = isa_slack,
                    gov_type_and_data = gov_data,
                    avr_type_and_data = avr_data,
                    loc_load_type_and_data=loc_load_instances_data[
                        n2s_gens_with_loc_load_idxs[plant_idx]],
                    p_order  = p_order,
                    v_ref = v_ref,
                    ω_ref = ω_ref )) :
                        (plant_idx,
                         ( # plant_idx = plant_idx,
                           plant_type = plant_type,
                           gen_type_and_data = gen_data,
                           isa_slack = isa_slack,
                           gov_type_and_data = gov_data,
                           avr_type_and_data = avr_data,
                           p_order = p_order,
                           v_ref = v_ref,
                           ω_ref = ω_ref) )
                     for (plant_idx, plant_type,
                          gen_data, isa_slack,
                          gov_data, avr_data) in
                         zip(gens_nodes_idx,
                            plants_types,
                            gens_instances_data,
                            bool_isa_slack,
                            govs_instances_data,
                             avrs_instances_data) ]
        else

            return [ ( idx ∈ gens_with_loc_load_idx && first(
                gov_data != nothing )) ?
                   get_a_gen_plant_wt_gov_wt_loc_load_by_mpc(
                       idx, plant_type, isa_slack,
                       p_order, v_ref, ω_ref,
                       gen_data, gov_data, avr_data,
                       loc_load_instances_data,
                       n2s_gens_with_loc_load_idxs) :
                           (idx ∈ gens_with_loc_load_idx &&
                           first( gov_data == nothing )) ?
                           get_a_gen_plant_no_gov_wt_loc_load_by_mpc(
                               idx, plant_type, isa_slack,
                               p_order, v_ref, ω_ref,
                               gen_data, gov_data, avr_data,
                               loc_load_instances_data,
                               n2s_gens_with_loc_load_idxs) :
                                   (idx ∉ gens_with_loc_load_idx &&
                                   first( gov_data != nothing )) ?
                                   get_a_gen_plant_wt_gov_no_loc_load_by_mpc(
                                       idx, plant_type,
                                       isa_slack, p_order,
                                       v_ref, ω_ref, gen_data,
                                       gov_data, avr_data,
                                       loc_load_instances_data,
                                       n2s_gens_with_loc_load_idxs) :
                                           get_a_gen_plant_no_gov_no_loc_load_by_mpc(
                                               idx, plant_type,
                                               isa_slack, p_order,
                                               v_ref, ω_ref, gen_data,
                                               gov_data, avr_data )

                     for (idx, plant_type, gen_data,
                          isa_slack, gov_data, avr_data) in
                         zip(gens_nodes_idx, plants_types,
                            gens_instances_data, bool_isa_slack,
                             govs_instances_data,
                             avrs_instances_data) ]


        end

        
    end
    
end



#---------------------------------------------------


function get_load_nodes_plant_instances_data_by_mpc(
    mpc_bus;
    mpc_baseMVA =
        1.0,
    load_type =
        PQ_Const_I,
    load_node_type =
        plant_PQ_Const_I,
    by_components =
        false )

    load_nodes_static_tup_data =
        get_load_nodes_static_tup_data_by_mpc(
            mpc_bus;
            mpc_baseMVA = mpc_baseMVA )

    if by_components == false

        return [(idx,
                 (plant_type = load_node_type,
                  load_type_and_data =
                      (load_type,
                       load_node_static_tup_data))) 
         for (idx, load_node_static_tup_data) in
                second.(
                    load_nodes_static_tup_data)]
    else

        return [ (idx = idx,
                 plant_type = load_node_type,
                  components_type =
                      ( load = load_type, ),
                  components_data =
                      ( load = load_node_static_tup_data,) ) 
             for (idx, load_node_static_tup_data) in
                    second.(
                        load_nodes_static_tup_data) ]

    end
    
    
end


function get_transmission_nodes_plant_instances_data_by_mpc(
    mpc_bus; mpc_baseMVA=1.0,
    transmission_type = Trans_t2_Node,
    transmission_node_type = plant_Transmission_t2,
    by_components = false)

    transmission_nodes_idx =
        get_transmission_nodes_idx_by_mpc(
            mpc_bus )

    if length(transmission_nodes_idx) != 0

        transmission_nodes_static_tup_data =
            get_transmission_nodes_static_tup_data_by_mpc(
                mpc_bus;
                mpc_baseMVA=mpc_baseMVA)

        if by_components == false
            return [
                (idx,
                 (plant_type =
                     transmission_node_type,
                  transmission_type_and_data =
                      (transmission_type,
                       transmission_node_static_tup_data)))

                for (idx, transmission_node_static_tup_data) in
                    second.(
                        transmission_nodes_static_tup_data)]
        else

            return [
                (idx = idx,
                 plant_type =
                     transmission_node_type,
                 components_type =
                     ( transmission = transmission_type,),
                 components_data =
                     (transmission =
                     transmission_node_static_tup_data,) )

                for (idx, transmission_node_static_tup_data) in
                    second.(
                        transmission_nodes_static_tup_data)]
            
        end
        
    else

        return []

    end
        
end

#---------------------------------------------------
#---------------------------------------------------


function get_Dyn_Branches_data_list_by_mpc(
    mpc_branch;
    mpc_baseMVA=1.0 )

    branches_data =
        get_branches_data_and_types_by_mpc(
            mpc_branch,
            mpc_baseMVA=mpc_baseMVA )
    
        return [ branch_type == :line ?
            (idx=idx,
             edge_type = PiModelLine,
             edge_data =  branch_data ) :
                 (idx=idx,
                  edge_type = Transformer,
                  edge_data = branch_data )
                 for (idx, branch_data, branch_type) in
                     branches_data ]

    
     
end


function get_Dyn_Branches_data_by_mpc(
    mpc_branch;
    mpc_baseMVA=1.0 )

    branches_data =
        get_branches_data_and_types_by_mpc(
            mpc_branch,
            mpc_baseMVA=mpc_baseMVA )

    Branches =
        [ branch_type == :line ?
        (idx, (edge_type = PiModelLine,
               edge_data =  branch_data )) :
                   (idx, (edge_type = Transformer,
                          edge_data = branch_data ))
          for (idx, branch_data, branch_type) in
              branches_data ]

    
    return OrderedDict("branch$(no)" => branch
                       for (no, branch) in
                           Branches)
end



function get_Dyn_Nodes_data_by_mpc(
    dict_plants_gen_sym_type,
    dict_gen_sym_type,
    dict_gens_dyn_nt_params,

    dict_gov_sym_type,
    dict_gov_nt_params,

    dict_avr_sym_type,
    dict_avr_nt_params,
    
    dyn_gens,
    dyn_plants,

    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0,

    p_order = 1.0,
    v_ref = 1.0,
    ω_ref = ωs,

    load_type =
        PQ_Const_I,
    load_node_type =
        plant_PQ_Const_I,

    transmission_type =
        Trans_t2_Node,
    transmission_node_type =
        plant_Transmission_t2,
    by_components = false)

    ##
    
    gens_nodes_idx =
        get_gens_nodes_idx_by_mpc(
            mpc_bus)
    
    load_nodes_idx =
        get_load_nodes_idx_by_mpc(
            mpc_bus)
    
    transmission_nodes_idx =
        get_transmission_nodes_idx_by_mpc(
            mpc_bus)
    
    all_nodes_idx =
        get_all_nodes_idx_by_mpc(
            mpc_bus)
    
    # n2s_gens_idx =
    #     get_a_n2s_net_group(
    #         gens_nodes_idx)

    # n2s_load_nodes_idx =
    #     get_a_n2s_net_group(
    #         load_nodes_idx)
    
    # n2s_transmission_idxs =
    #     length(transmission_nodes_idx) != 0 ?
    #     get_a_n2s_net_group(
    #         transmission_nodes_idx;
    #         transmission_group = true) : []

    
    n2s_gens_idx =
        get_n2s_any(
            gens_nodes_idx)

    n2s_load_nodes_idx =
        get_n2s_any(
            load_nodes_idx)
    
    n2s_transmission_idxs =
        length(transmission_nodes_idx) != 0 ?
        get_n2s_any(
            transmission_nodes_idx) : get_n2s_any(
                transmission_nodes_idx;
                nothing_bool= true)
    
    gens_plant_instances_data =
        get_gens_plant_instances_data_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus,
            mpc_gen;
            mpc_baseMVA = mpc_baseMVA,

            p_order = p_order,
            v_ref = v_ref,
            ω_ref = ω_ref,
            by_components = by_components )
    
    load_nodes_plant_instances_data =
        get_load_nodes_plant_instances_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA,
            load_type =
                load_type,
            load_node_type =
                load_node_type,
            by_components =
                by_components)

    transmission_nodes_plant_instances_data =
        get_transmission_nodes_plant_instances_data_by_mpc(
            mpc_bus;
            mpc_baseMVA=mpc_baseMVA,
            transmission_type = transmission_type,
            transmission_node_type =
                transmission_node_type,
            by_components = by_components)

    if length(transmission_nodes_idx) != 0        
        tup_idx_and_plants_data =
            [idx ∈ gens_nodes_idx ?
            gens_plant_instances_data[n2s_gens_idx[idx]] :
            idx ∈ load_nodes_idx ?
            load_nodes_plant_instances_data[
                n2s_load_nodes_idx[idx]] :
                    transmission_nodes_plant_instances_data[
                        n2s_transmission_idxs[idx]]
              for idx in all_nodes_idx ]
        
    else
        
        tup_idx_and_plants_data =
            [idx ∈ gens_nodes_idx ?
            gens_plant_instances_data[n2s_gens_idx[idx]] :
            load_nodes_plant_instances_data[
                n2s_load_nodes_idx[idx]]
             for idx in all_nodes_idx]
    end
        
    return OrderedDict( "bus$(idx)" => plant_data
                        for (idx, plant_data) in
                            tup_idx_and_plants_data )    
    ##
    
end


#---------------------------------------------------
#---------------------------------------------------

function get_Dyn_Branches_data_by_components_by_mpc(
    mpc_branch;
    vec_edge_type = [] )
                
        return get_branches_impedance_data_by_components_by_mpc(
            mpc_branch;
            vec_edge_type = vec_edge_type )

    
     
end


function get_Dyn_Nodes_data_by_components_by_mpc(
    dict_plants_gen_sym_type,
    dict_gen_sym_type,
    dict_gens_dyn_nt_params,

    dict_gov_sym_type,
    dict_gov_nt_params,

    dict_avr_sym_type,
    dict_avr_nt_params,
    
    dyn_gens,
    dyn_plants,

    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0,

    p_order = 1.0,
    v_ref = 1.0,
    ω_ref = ωs,

    load_type =
        PQ_Const_I,
    load_node_type =
        plant_PQ_Const_I,

    transmission_type =
        Trans_t2_Node,
    transmission_node_type =
        plant_Transmission_t2,
    by_components = true )

    ## Added by AAY
    
    # gens_nodes_idx =
    #     get_gens_nodes_idx_by_mpc(
    #         mpc_bus)
    
    # load_nodes_idx =
    #     get_load_nodes_idx_by_mpc(
    #         mpc_bus)
    
    # transmission_nodes_idx =
    #     get_transmission_nodes_idx_by_mpc(
    #         mpc_bus)
    
    # all_nodes_idx =
    #     get_all_nodes_idx_by_mpc(
    #         mpc_bus)


    #----------------------
    
    # gens_nodes_idx =
    #     get_unsorted_gens_nodes_idx_by_mpc(
    #         mpc_bus)
    
    # load_nodes_idx =
    #     get_unsorted_load_nodes_idx_by_mpc(
    #         mpc_bus)
    
    # transmission_nodes_idx =
    #     get_transmission_nodes_idx_by_mpc(
    #         mpc_bus ;sorted_bool = false)
    
    # all_nodes_idx =
    #     get_unsorted_all_nodes_idx_by_mpc(
    #         mpc_bus)
    
    
    # n2s_gens_idx =
    #     get_n2s_any(
    #         gens_nodes_idx)

    # n2s_load_nodes_idx =
    #     get_n2s_any(
    #         load_nodes_idx)
    
    # n2s_transmission_idxs =
    #     length(transmission_nodes_idx) != 0 ?
    #     get_n2s_any(
    #         transmission_nodes_idx) : get_n2s_any(
    #             transmission_nodes_idx;
    #             nothing_bool= true)

    gens_plant_instances_data =
        get_gens_plant_instances_data_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus,
            mpc_gen;
            mpc_baseMVA = mpc_baseMVA,

            p_order = p_order,
            v_ref = v_ref,
            ω_ref = ω_ref,
            by_components = by_components)
    
    load_nodes_plant_instances_data =
        get_load_nodes_plant_instances_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA,
            load_type =
                load_type,
            load_node_type =
                load_node_type,
            by_components =
                by_components)

    transmission_nodes_plant_instances_data =
        get_transmission_nodes_plant_instances_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA,
            transmission_type =
                transmission_type,
            transmission_node_type =
                transmission_node_type,
            by_components = by_components)

        
    return Dict(
        :plant_generators =>
            gens_plant_instances_data,
        :plant_loads =>
            load_nodes_plant_instances_data,
        :plant_transmissions =>
            transmission_nodes_plant_instances_data )
       
end


function get_static_Nodes_data_by_components_by_mpc(
    dict_plants_gen_sym_type,
    dict_gen_sym_type,
    
    dyn_gens,
    dyn_plants,

    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0,

    load_type =
        PQ_Const_I,
    load_node_type =
        plant_PQ_Const_I,

    transmission_type =
        Trans_t2_Node,
    transmission_node_type =
        plant_Transmission_t2,
    by_components = true )

    ##
    
    gens_nodes_idx =
        get_gens_nodes_idx_by_mpc(
            mpc_bus)
    
    load_nodes_idx =
        get_load_nodes_idx_by_mpc(
            mpc_bus)
    
    transmission_nodes_idx =
        get_transmission_nodes_idx_by_mpc(
            mpc_bus)
    
    all_nodes_idx =
        get_all_nodes_idx_by_mpc(
            mpc_bus)
    
    n2s_gens_idx =
        get_n2s_any(
            gens_nodes_idx)

    n2s_load_nodes_idx =
        get_n2s_any(
            load_nodes_idx)
    
    n2s_transmission_idxs =
        length(transmission_nodes_idx) != 0 ?
        get_n2s_any(
            transmission_nodes_idx) : get_n2s_any(
                transmission_nodes_idx;
                nothing_bool= true)

    gens_plant_instances_data =
        get_static_gens_plant_instances_data_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,

            dyn_gens,
            dyn_plants,

            mpc_bus,
            mpc_gen;
            
            mpc_baseMVA =
                mpc_baseMVA,

            by_components =
                by_components)
    
    load_nodes_plant_instances_data =
        get_load_nodes_plant_instances_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA,
            load_type =
                load_type,
            load_node_type =
                load_node_type,
            by_components =
                by_components)

    transmission_nodes_plant_instances_data =
        get_transmission_nodes_plant_instances_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA,
            transmission_type =
                transmission_type,
            transmission_node_type =
                transmission_node_type,
            by_components = by_components)

        
    return Dict(
        :plant_generators =>
            gens_plant_instances_data,
        :plant_loads =>
            load_nodes_plant_instances_data,
        :plant_transmissions =>
            transmission_nodes_plant_instances_data )
       
end


function get_multi_gens_Dyn_Nodes_data_by_components_by_mpc(
    dict_plants_gen_sym_type,
    dict_gen_sym_type,
    dict_gens_dyn_nt_params,

    dict_gov_sym_type,
    dict_gov_nt_params,

    dict_avr_sym_type,
    dict_avr_nt_params,
    
    dyn_gens,
    dyn_plants,

    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0,

    p_order = 1.0,
    v_ref = 1.0,
    ω_ref = ωs,

    load_type =
        PQ_Const_I,
    load_node_type =
        plant_PQ_Const_I,

    transmission_type =
        Trans_t2_Node,
    transmission_node_type =
        plant_Transmission_t2,
    by_components = true )

    ##

    #--------------------------------------

    (;multi_gens_idx,
     multi_gen_bool,
     multi_gens_nodes) =
         NamedTupleTools.select(
             get_multi_gens_idx_wt_multi_gen_bool(
                 mpc_gen),
             (:multi_gens_idx,
              :multi_gen_bool,
              :multi_gens_nodes) )

    #--------------------------------------

    
    gens_nodes_idx =
        get_gens_nodes_idx_by_mpc(
            mpc_bus)
    
    load_nodes_idx =
        get_load_nodes_idx_by_mpc(
            mpc_bus)
    
    transmission_nodes_idx =
        get_transmission_nodes_idx_by_mpc(
            mpc_bus)
    
    all_nodes_idx =
        get_all_nodes_idx_by_mpc(
            mpc_bus)

    #--------------------------------------
    
    n2s_gens_idx =
        get_n2s_any(
            gens_nodes_idx)

    n2s_load_nodes_idx =
        get_n2s_any(
            load_nodes_idx)
    
    n2s_transmission_idxs =
        length(transmission_nodes_idx) != 0 ?
        get_n2s_any(
            transmission_nodes_idx) : get_n2s_any(
                transmission_nodes_idx;
                nothing_bool= true)

    #--------------------------------------

    gens_plant_instances_data =
        get_multi_gens_plant_instances_data_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus,
            mpc_gen;
            mpc_baseMVA = mpc_baseMVA,

            p_order = p_order,
            v_ref = v_ref,
            ω_ref = ω_ref,
            by_components =
                by_components)
    
    load_nodes_plant_instances_data =
        get_load_nodes_plant_instances_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA,
            load_type =
                load_type,
            load_node_type =
                load_node_type,
            by_components =
                by_components)

    transmission_nodes_plant_instances_data =
        get_transmission_nodes_plant_instances_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA,
            transmission_type =
                transmission_type,
            transmission_node_type =
                transmission_node_type,
            by_components = by_components)

        
    return Dict(
        :plant_generators =>
            gens_plant_instances_data,
        :plant_loads =>
            load_nodes_plant_instances_data,
        :plant_transmissions =>
            transmission_nodes_plant_instances_data )
       
end


#---------------------------------------------------
#---------------------------------------------------
#  Data conversion related functions
#---------------------------------------------------
#---------------------------------------------------


function create_a_default_case_mpc_branch_type(
    case_name;
    data_dir      = "",
    case_data_dir = "",

    mpc_branch_column_select = "",
    mpc_branch_line_type     =  "PiModelLine",
    mpc_branch_transformer_type =  "Transformer"  )

    #--------------------------------------

    mpc_branch_line_type == "" ?
        "PiModelLine" : mpc_branch_line_type
    
    mpc_branch_transformer_type == "" ?
        "Transformer" : mpc_branch_transformer_type
    
    #--------------------------------------
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end
    
    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
        
    end
    
    #--------------------------------------

    csv_branch_file =
        joinpath(case_data_dir,
                 "mpc_branch.csv")


    csv_branch_type_file =
        joinpath(case_data_dir,
                 "mpc_branch_type.csv")
    
    #--------------------------------------    

    if mpc_branch_column_select == ""
        
        mpc_branch_column_select =
            ["ratio",
             "angle",
             "status"]
        
    end    
    
    #--------------------------------------

    mpc_branch_selected_data =
        CSV.File( csv_branch_file;
                  select = mpc_branch_column_select )

    branch_type =
        [  a_tup.ratio == 0.0 ?
        mpc_branch_line_type : mpc_branch_transformer_type
           for a_tup in
               CSV.File(csv_branch_file;
                        select =
                            mpc_branch_column_select ) ]

    idx = 1:length(branch_type)


    # DataFrame( idx = idx , branch_type = branch_type )
    
    CSV.write( csv_branch_type_file,
               DataFrame( idx = idx ,
                          branch_type = branch_type ) )
    
    return nothing
    
end



function create_a_default_case_mpc_load_type(
    case_name;
    data_dir      = "",
    case_data_dir = "",

    mpc_bus_column_select      = "",
    mpc_load_node_type         = "PQ_Const_P",
    mpc_transmission_node_type =  "Trans_t2_Node" )

    #--------------------------------------

    mpc_load_node_type == "" ?
        "PQ_Const_P" : mpc_load_node_type
    
    mpc_transmission_node_type == "" ?
        "Trans_t2_Node" : mpc_transmission_node_type
    
    #--------------------------------------
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    
    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
        
    end
        
    #--------------------------------------
    
    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")

    
    csv_load_type_file =
        joinpath(case_data_dir,
                 "mpc_load_type.csv")
    
    #--------------------------------------    

    if mpc_bus_column_select == ""

         mpc_bus_column_select =
             ["bus_i", "type",
              "Pd", "Qd"]
    end

    #--------------------------------------
    
    mpc_bus_selected_data =
        CSV.File(csv_bus_file;
                 select=mpc_bus_column_select )

    #--------------------------------------    

    #     mpc_load_node_type, mpc_transmission_node_type

    idx_wt_load_type =
        [  (a_tup.Pd == 0 || a_tup.Pd == 0.0) && (
            a_tup.Qd == 0 || a_tup.Qd == 0.0) ?
                (a_tup.bus_i,
                             mpc_transmission_node_type) : (
            a_tup.bus_i, mpc_load_node_type) 
           for a_tup in mpc_bus_selected_data
               if a_tup.type == 1 ]

    idx    = first.( idx_wt_load_type )
    
    load_type = second.( idx_wt_load_type )

    CSV.write( csv_load_type_file,
               DataFrame( idx = idx ,
                          load_type = load_type ) )
    
    return nothing
    
   
end



function create_a_default_case_dyn_gens_file(
    case_name;
    data_dir = "",
    case_data_dir = "",

    dyn_gens_file = "",
    
    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6", 
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine = "SM_2axis_cb_v6",
    synchronous_condenser = "SC_2axis_cb_v6",

    synchronous_machine_dynamic_param =
        "gen_dynamic_paras_ieee_14_b" )

    #--------------------------------------    
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    
    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
        
    end
        
    #--------------------------------------
    
    dyn_data_dir  =
        joinpath(case_data_dir,
                 "dyn")

    if !(isdir( dyn_data_dir ))

        mkpath( dyn_data_dir )

    end

    #--------------------------------------

    if dyn_gens_file == ""

        dyn_gens_file =
            joinpath(dyn_data_dir,
                     "dyn_gen.csv")
        
    end

    #--------------------------------------

    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")


    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")
    
    #--------------------------------------    

    mpc_bus_column_select =
        ["bus_i",
         "type",
         "Pd",
         "Qd"]
    
    mpc_gen_column_select =
        ["bus",
         "Pg",
         "Qg"]
    
    #--------------------------------------
        
    mpc_bus_selected_data =
        CSV.File( csv_bus_file;
                  select = mpc_bus_column_select )

    gens_nodes_idx_in_bus_cvs =
        [ a_tup.bus_i for a_tup in
             mpc_bus_selected_data if
                 a_tup.type == 3 || a_tup.type == 2 ]

    slack_nodes_idx_in_bus_cvs =
        [a_tup.bus_i for a_tup in
             mpc_bus_selected_data if
                 a_tup.type == 3 ]
    
    mpc_gen_selected_data =
        CSV.File( csv_gen_file;
                  select = mpc_gen_column_select )

    # gens_idx = gens_nodes_idx_in_bus_cvs

    """ This is important, when there are more than one
    generators per bus """

    # [idx for idx in mpc_gen_selected_data.bus]    

    gens_idx = mpc_gen_selected_data.bus
    
    gens_countmap =
        countmap(gens_idx)

    multi_gens_nodes = [k for (k, v) in
                           gens_countmap if v > 1]

    multi_gen_bool =
        length(multi_gens_nodes) == 0 ? false : true

    
    gens_idx = multi_gen_bool == false ? gens_idx : [
        [ idxs ∈ multi_gens_nodes ?
            [ "$(idxs)-$(idx)"
              for idx in 1:gens_countmap[ idxs ]] :
                  "$(idxs)" for idxs in
                      unique(gens_idx)]...;]

    #--------------------------------------
    
    bus_nodes_with_demands_idx =
        [ a_tup.bus_i for a_tup in
             mpc_bus_selected_data if (
                 a_tup.Pd != 0 || a_tup.Pd != 0.0 ) || (
                 a_tup.Qd != 0 || a_tup.Qd != 0.0)]

    #--------------------------------------

    bus_SM =
        [ a_tup.bus for a_tup in
             mpc_gen_selected_data if (
                 (a_tup.bus ∈ slack_nodes_idx_in_bus_cvs) ||
                     ( a_tup.Pg != 0 || a_tup.Pg != 0.0 )) ]

    bus_SC =
        [ a_tup.bus for a_tup in
             mpc_gen_selected_data if (
                 (a_tup.bus ∉ slack_nodes_idx_in_bus_cvs) &&
                     ( a_tup.Pg == 0 || a_tup.Pg == 0.0 )) ]
    
    #--------------------------------------

    """ These are bus_SM and bus_SC with local loads """
    
    bus_SM_in_loc_load =
        [ idx for idx in bus_SM
             if idx ∈  bus_nodes_with_demands_idx ]

    bus_SC_in_loc_load = length(bus_SC) == 0 ? [] :
        [ idx for idx in bus_SC
             if idx ∈  bus_nodes_with_demands_idx ]

    #--------------------------------------
    
    # gens_idx = sort(gens_idx)
    
    sym_gen_type =  multi_gen_bool == false ? [
        idx ∈ bus_SC_in_loc_load ?
            synchronous_condenser_wt_loc_load :
            idx ∈ bus_SM_in_loc_load ?
            synchronous_machine_wt_loc_load :
            idx ∈ bus_SC ? synchronous_condenser :
            synchronous_machine
        for idx in gens_idx
             ] :  [
              get_node_idx_multi(idx) ∈ bus_SC_in_loc_load ?
                     synchronous_condenser_wt_loc_load :
               get_node_idx_multi(idx) ∈ bus_SM_in_loc_load ?
                     synchronous_machine_wt_loc_load :
                     get_node_idx_multi(idx) ∈ bus_SC ?
                     synchronous_condenser :
                     synchronous_machine
                 for idx in gens_idx ] 
    
    sym_gen_dynamic_para =
        [synchronous_machine_dynamic_param
         for idx in gens_idx ]

    #--------------------------------------
    
    CSV.write( dyn_gens_file,
               DataFrame( bus = gens_idx ,
                          sym_gen_type = sym_gen_type,
                          sym_gen_dynamic_para =
                              sym_gen_dynamic_para  ) )

    return nothing

end


function create_a_default_case_dyn_plants_file(
    case_name;
    data_dir      = "",
    case_data_dir = "",

    dyn_plants_file = "",
    
    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6",
    
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    
    synchronous_machine =
        "SM_2axis_cb_v6",
    
    synchronous_condenser =
        "SC_2axis_cb_v6",

    plant_wt_loc_load =
        "plant_wt_loc_load_v6",
    
    plant_no_gov_wt_loc_load =
        "plant_no_gov_wt_loc_load_v6",
    
    plant_no_gov = "plant_no_gov",
    
    plant = "plant_cb_v6",
    
    avr_param = "avr_t1_cb_sauer__1_param",
    
    gov_param = "gov_t1_cb_sauer__1_param" )

    #--------------------------------------    
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    
    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
        
    end
    
    #--------------------------------------
    
    dyn_data_dir  =
        joinpath(case_data_dir,
                 "dyn")

    if !(isdir( dyn_data_dir ))

        mkpath( dyn_data_dir )

    end

    #--------------------------------------

    if dyn_plants_file == ""
        
        dyn_plants_file =
            joinpath(dyn_data_dir,
                     "dyn_plant.csv")
    end
    
    #--------------------------------------
    #--------------------------------------

    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")


    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")
    
    #--------------------------------------    

    mpc_bus_column_select =
        ["bus_i",
         "type",
         "Pd",
         "Qd"]
    
    mpc_gen_column_select =
        ["bus",
         "Pg",
         "Qg"]
    
    #--------------------------------------


    
    mpc_bus_selected_data =
        CSV.File( csv_bus_file;
                  select = mpc_bus_column_select )

    gens_nodes_idx_in_bus_cvs =
        [a_tup.bus_i for a_tup in
             mpc_bus_selected_data if
                 a_tup.type == 3 || a_tup.type == 2 ]

    slack_nodes_idx_in_bus_cvs =
        [a_tup.bus_i for a_tup in
             mpc_bus_selected_data if
                 a_tup.type == 3 ]
    
    mpc_gen_selected_data =
        CSV.File( csv_gen_file;
                  select = mpc_gen_column_select )

    # gens_idx = sort( gens_nodes_idx_in_bus_cvs )

    gens_idx = [ a_tup.bus for a_tup in
                    mpc_gen_selected_data ]
    
    gens_countmap =
        countmap(gens_idx)

    multi_gens_nodes = [k for (k, v) in
                           gens_countmap if v > 1]

    multi_gen_bool =
        length(multi_gens_nodes) == 0 ? false : true
    
    gens_idx = multi_gen_bool == false ? gens_idx : [
        [ idxs ∈ multi_gens_nodes ?
            [ "$(idxs)-$(idx)"
              for idx in 1:gens_countmap[ idxs ]] :
                  "$(idxs)" for idxs in
                      unique(gens_idx)]...; ]

    #--------------------------------------    

    nodes_idx_with_type =
        [ (a_tup.bus_i, a_tup.type )
          for a_tup in
              mpc_bus_selected_data ]

    # nodes_idx_with_type =
    #     sort( nodes_idx_with_type, by = x->x[ 1 ] )

    nodes_idx  = first.(
        nodes_idx_with_type )

    nodes_type = second.(
        nodes_idx_with_type )

    a_slack_bus_idx =
        [ first(a_tup) for a_tup in nodes_idx_with_type
             if second(a_tup) == 3  ]
    
    #--------------------------------------

    bus_nodes_with_demands_idx =
        [ a_tup.bus_i for a_tup in
             mpc_bus_selected_data if
                 ( a_tup.Pd != 0 || a_tup.Pd != 0.0 ) || (
                     a_tup.Qd != 0 || a_tup.Qd != 0.0)]

    #--------------------------------------

    bus_SM =
        [ a_tup.bus for a_tup in
             mpc_gen_selected_data if (
                 (a_tup.bus ∈ slack_nodes_idx_in_bus_cvs) ||
                     ( a_tup.Pg != 0 || a_tup.Pg != 0.0 )) ]

    bus_SC =
        [ a_tup.bus for a_tup in
             mpc_gen_selected_data if (
                 (a_tup.bus ∉ slack_nodes_idx_in_bus_cvs) &&
                     ( a_tup.Pg == 0 || a_tup.Pg == 0.0 )) ]
        
    #--------------------------------------

    bus_SM_in_loc_load =
        [ idx for idx in bus_SM
             if idx ∈  bus_nodes_with_demands_idx ]

    bus_SC_in_loc_load = length(bus_SC) == 0 ? [] :
        [ idx for idx in bus_SC
             if idx ∈  bus_nodes_with_demands_idx ]

    #--------------------------------------
    
    # gens_idx = sort(gens_idx)

    Plant_type = multi_gen_bool == false ?
        [ idx ∈ bus_SC_in_loc_load ?
        plant_no_gov_wt_loc_load :
        idx ∈ bus_SM_in_loc_load ?  plant_wt_loc_load  :
        idx ∈ bus_SC ? plant_no_gov : plant
          for idx in gens_idx  ] :
    [ get_node_idx_multi(idx) ∈ bus_SC_in_loc_load ?
      plant_no_gov_wt_loc_load :
      get_node_idx_multi(idx) ∈ bus_SM_in_loc_load ?
      plant_wt_loc_load :
      get_node_idx_multi(idx) ∈ bus_SC ?
      plant_no_gov : plant
       for idx in gens_idx ]

    
    Gen = multi_gen_bool == false ?
        [ idx ∈ bus_SC_in_loc_load ?
        synchronous_condenser_wt_loc_load :
        idx ∈ bus_SM_in_loc_load ?
        synchronous_machine_wt_loc_load :
        idx ∈ bus_SC ? synchronous_condenser :
        synchronous_machine
          for idx in gens_idx ] :
              [get_node_idx_multi(idx) ∈ bus_SC_in_loc_load ?
        synchronous_condenser_wt_loc_load :
        get_node_idx_multi(idx) ∈ bus_SM_in_loc_load ?
        synchronous_machine_wt_loc_load :
        get_node_idx_multi(idx) ∈ bus_SC ?
        synchronous_condenser : synchronous_machine
          for idx in gens_idx ]

    isa_slack = multi_gen_bool == false ?
        [ idx ∈ a_slack_bus_idx ? true : false
                  for idx in gens_idx ] :
              [ get_node_idx_multi(idx) ∈ a_slack_bus_idx ?
                      true : false
                  for idx in gens_idx ]

    Gov = [ idx ∈ bus_SC ?
        :nothing : gov_param
            for idx in gens_idx ]

    Exc = [ avr_param
            for idx in gens_idx ]
    	
    #--------------------------------------
    
    CSV.write( dyn_plants_file,
               DataFrame(bus = gens_idx,
                         Plant_type = Plant_type,
                         Gen = Gen,
                         isa_slack = isa_slack,
                         Gov = Gov,
                         Exc = Exc ) )

    return nothing

    
end



function create_a_default_case_net_data_xlsx_file(
    case_name;
    data_dir        = "",
    case_data_dir   = "",
    dyn_data_dir    = "",

    xlsx_default_data_file =
        "",

    mpc_branch_column_select    = "",
    mpc_branch_line_type        = "PiModelLine",
    mpc_branch_transformer_type = "Transformer",

    mpc_bus_column_select       = "",
    mpc_load_node_type          = "PQ_Const_P",
    mpc_transmission_node_type  = "Trans_t2_Node",

    synchronous_machine_wt_loc_load   =
        "SM_2axis_wt_loc_load_cb_v6",
    
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine   =
        "SM_2axis_cb_v6",
    
    synchronous_condenser =
        "SC_2axis_cb_v6",

    plant_wt_loc_load           =
        "plant_wt_loc_load_v6",
    plant_no_gov_wt_loc_load    =
        "plant_no_gov_wt_loc_load_v6",
    plant_no_gov        =
        "plant_no_gov",
    plant               =
        "plant_cb_v6",

    synchronous_machine_dynamic_param =
        "gen_dynamic_paras_ieee_14_b",
    avr_param = "avr_t1_cb_sauer__1_param",
    gov_param = "gov_t1_cb_sauer__1_param",
    
    by_components = true )

    #--------------------------------------    
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    
    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
        
    end
        

    #--------------------------------------
    
    xlsx_data_dir  =
        joinpath(case_data_dir,
                 "xlsx")

    if !(isdir( xlsx_data_dir ))

        mkpath( xlsx_data_dir)

    end

    #--------------------------------------

    if xlsx_default_data_file == ""
    
        xlsx_file =
            joinpath(xlsx_data_dir,
                     "net-static-data.xlsx")
        
    else

        xlsx_file =
            joinpath(xlsx_data_dir,
                     xlsx_default_data_file)
        
    end
    
    #--------------------------------------

    csv_branch_file =
        joinpath(case_data_dir,
                 "mpc_branch.csv")

    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")

    csv_gencost_file =
        joinpath(case_data_dir,
                 "mpc_gencost.csv")
    
    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")

    csv_scalar_file =
        joinpath(case_data_dir,
                 "mpc_scalar.csv")

    #--------------------------------------
    #--------------------------------------

    csv_load_type_file  = 
        joinpath(case_data_dir,
                 "mpc_load_type.csv")

    csv_branch_type_file =
        joinpath(case_data_dir,
                 "mpc_branch_type.csv")
    

    if !( isfile( csv_load_type_file ) )
        
        create_a_default_case_mpc_load_type(
            case_name;
            
            data_dir      = data_dir,
            
            case_data_dir = case_data_dir ,

            mpc_bus_column_select =
                mpc_bus_column_select,
            
            mpc_load_node_type =
                mpc_load_node_type,
            
            mpc_transmission_node_type =
                mpc_transmission_node_type )
    end
    

    if !( isfile( csv_branch_type_file ) )
        
        create_a_default_case_mpc_branch_type(
            case_name;
            
            data_dir      = data_dir,
            
            case_data_dir = case_data_dir,

            mpc_branch_column_select =
                mpc_branch_column_select,
            
            mpc_branch_line_type =
                mpc_branch_line_type ,
            
            mpc_branch_transformer_type =
                mpc_branch_transformer_type  )
    end
    
    #--------------------------------------
    
    if dyn_data_dir == ""
    
        dyn_data_dir  =
            joinpath(case_data_dir,
                     "dyn")
        
    end


    if !(isdir( dyn_data_dir ))

        mkpath( dyn_data_dir )

    end

    
    dyn_gens_file =
        joinpath(dyn_data_dir,
                 "dyn_gen.csv")

    dyn_plants_file =
        joinpath(dyn_data_dir,
                 "dyn_plant.csv")
    

    if !( isfile( dyn_gens_file ) )
        
        create_a_default_case_dyn_gens_file(
            case_name;
            data_dir,
            case_data_dir,

            dyn_gens_file,

            synchronous_machine_wt_loc_load, 
            synchronous_condenser_wt_loc_load,
            synchronous_machine,
            synchronous_condenser,

            synchronous_machine_dynamic_param)
    end


    if !( isfile( dyn_plants_file ) )
        
        create_a_default_case_dyn_plants_file(
            case_name;
            data_dir,
            case_data_dir,

            dyn_plants_file,

            synchronous_machine_wt_loc_load, 
            synchronous_condenser_wt_loc_load,
            synchronous_machine,
            synchronous_condenser,

            plant_wt_loc_load,
            plant_no_gov_wt_loc_load,
            plant_no_gov,
            plant,

            avr_param,
            gov_param)
        
    end
    

    #--------------------------------------

    mpc_branch_column_select =
        String["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_bus_column_select =
        String["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_gen_column_select =
        String["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_scalar_column_select =
        String["mpc_baseMVA" ]

    dyn_gens_column_select =
        String["bus","sym_gen_type",
         "sym_gen_dynamic_para"]
    
    dyn_plants_column_select =
       String["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    # dyn_gens_data_types =
    #     [Int, Symbol, Symbol]

    # dyn_plants_data_types =
    #     [ Int, Symbol, Symbol, Bool, Symbol, Symbol]

    (;multi_gen_bool,) =
        NamedTupleTools.select(
            check_multi_gens_bool_by_csv_file(
                csv_gen_file),
            (:multi_gen_bool,) )

    dyn_gens_data_types = multi_gen_bool == false ? [Int, String, String] : [String, String, String]

    dyn_plants_data_types = multi_gen_bool == false ? [ Int, String, String, Bool, String, String] : 
        [ String, String, String, Bool, String, String]
    
    #--------------------------------------

    list_data_types = [
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        dyn_gens_data_types,
        dyn_plants_data_types ]
    
    list_data_select = Union{String,Any}[
        [],
        mpc_branch_column_select,
        mpc_bus_column_select,
        mpc_gen_column_select,
        [],
        [],
        mpc_scalar_column_select,
        dyn_gens_column_select,
        dyn_plants_column_select ]

    list_csv_files =
        [ csv_branch_type_file,
          csv_branch_file,
          csv_bus_file,
          csv_gen_file,
          csv_gencost_file,
          csv_load_type_file,
          csv_scalar_file,
          dyn_gens_file,
          dyn_plants_file ]


    list_df =
        [( a_data_types_list == [] && a_data_select_list == []) ?
        CSV.read(a_csv_file, DataFrame) :
        (a_data_types_list == [] && a_data_select_list != []) ?
        CSV.read(a_csv_file, DataFrame; select = a_data_select_list ) :
        (a_data_types_list != [] && a_data_select_list == []) ?
        CSV.read(a_csv_file, DataFrame; types  = a_data_types_list) :
        CSV.read(a_csv_file, DataFrame; select = a_data_select_list, types = a_data_types_list)
          for (a_data_types_list,a_data_select_list,a_csv_file ) in
              zip(list_data_types,list_data_select,list_csv_files)]
        
    list_sheet_names =
        [ "mpc_branch_type",
          "mpc_branch",
          "mpc_bus",
          "mpc_gen",
          "mpc_gencost",
          "mpc_load_type",
          "mpc_scalar",
          "dyn_gen",
          "dyn_plant"]
    
    @assert length(list_sheet_names) == length(list_df)

    XLSX.openxlsx(xlsx_file, mode="w") do xf
        for i in eachindex( list_sheet_names )
            sheet_name = list_sheet_names[i]
            df = list_df[i]

            if i == firstindex( list_sheet_names )
                sheet = xf[1]
                XLSX.rename!(sheet, sheet_name)
                XLSX.writetable!(sheet, df)
            else
                sheet = XLSX.addsheet!(xf, sheet_name)
                XLSX.writetable!(sheet, df)        
            end
        end
    end    

    return nothing
                
end


function create_a_default_case_net_data_xlsx_file(
    case_name,
    nothing;
    data_dir            = "",
    case_data_dir       = "",
    dyn_data_dir        = "",

    xlsx_default_data_file      = "",
    
    dyn_gens_default_file       = "",
    dyn_plants_default_file     = "",

    mpc_branch_column_select     = "",
    mpc_branch_line_type         =  "PiModelLine",
    mpc_branch_transformer_type  =  "Transformer",

    mpc_bus_column_select = "",
    mpc_load_node_type =
        "PQ_Const_P",
    mpc_transmission_node_type =
        "Trans_t2_Node",

    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6", 
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine =
        "SM_2axis_cb_v6",
    synchronous_condenser =
        "SC_2axis_cb_v6",

    plant_wt_loc_load =
        "plant_wt_loc_load_v6",
    plant_no_gov_wt_loc_load =
        "plant_no_gov_wt_loc_load_v6",
    plant_no_gov =
        "plant_no_gov",
    plant =
        "plant_cb_v6",

    synchronous_machine_dynamic_param =
        "gen_dynamic_paras_ieee_14_b",
    avr_param = "avr_t1_cb_sauer__1_param",
    gov_param = "gov_t1_cb_sauer__1_param",
    
    by_components = true)
 
    #--------------------------------------    
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    
    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
        
    end
        
    #--------------------------------------
    
    xlsx_data_dir  =
        joinpath(case_data_dir,
                 "xlsx")

    if !(isdir( xlsx_data_dir ))

        mkpath( xlsx_data_dir)

    end

    #--------------------------------------
    
    # xlsx_file =
    #     joinpath(xlsx_data_dir,
    #              "net-static-data.xlsx")

    if xlsx_default_data_file == ""
    
        xlsx_file =
            joinpath(xlsx_data_dir,
                     "net-default-static-data.xlsx")
        
    else

        xlsx_file =
            joinpath(xlsx_data_dir,
                     xlsx_default_data_file)
        
    end
    
    #--------------------------------------

    csv_branch_file =
        joinpath(case_data_dir,
                 "mpc_branch.csv")

    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")

    csv_gencost_file =
        joinpath(case_data_dir,
                 "mpc_gencost.csv")
    
    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")

    csv_scalar_file =
        joinpath(case_data_dir,
                 "mpc_scalar.csv")

    #--------------------------------------
    #--------------------------------------

    csv_load_type_file  = 
        joinpath(case_data_dir,
                 "mpc_load_type.csv")

    csv_branch_type_file =
        joinpath(case_data_dir,
                 "mpc_branch_type.csv")
    

    if !( isfile( csv_load_type_file ) )
        
        create_a_default_case_mpc_load_type(
            case_name;
            
            data_dir =
                data_dir,
            
            case_data_dir =
                case_data_dir,

            mpc_bus_column_select =
                mpc_bus_column_select,
            
            mpc_load_node_type =
                mpc_load_node_type,
            
            mpc_transmission_node_type =
                mpc_transmission_node_type )
    end
    

    if !( isfile( csv_branch_type_file ) )
        
        create_a_default_case_mpc_branch_type(
            case_name;
            
            data_dir =
                data_dir,
            
            case_data_dir =
                case_data_dir,

            mpc_branch_column_select =
                mpc_branch_column_select,
            
            mpc_branch_line_type =
                mpc_branch_line_type ,
            
            mpc_branch_transformer_type =
                mpc_branch_transformer_type  )
    end
    
    #--------------------------------------
    
    if dyn_data_dir == ""
    
        dyn_data_dir  =
            joinpath(case_data_dir,
                     "dyn")
        
    end


    if !(isdir( dyn_data_dir ))

        mkpath( dyn_data_dir )

    end

    if dyn_gens_default_file == ""

        dyn_gens_file =
            joinpath(dyn_data_dir,
                     "dyn_gen.csv")

        
    else

        dyn_gens_file =
            joinpath(dyn_data_dir,
                     dyn_gens_default_file)
        
    end

    if dyn_plants_default_file == ""
        
        dyn_plants_file =
            joinpath(dyn_data_dir,
                     "dyn_plant.csv")
        
    else
        
        dyn_plants_file =
            joinpath(dyn_data_dir,
                     dyn_plants_default_file)
        
    end
    
    

    if !( isfile( dyn_gens_file ) )
        
        create_a_default_case_dyn_gens_file(
            case_name;
            data_dir,
            case_data_dir,

            dyn_gens_file,

            synchronous_machine_wt_loc_load, 
            synchronous_condenser_wt_loc_load,
            synchronous_machine,
            synchronous_condenser,

            synchronous_machine_dynamic_param)
    end


    if !( isfile( dyn_plants_file ) )
        
        create_a_default_case_dyn_plants_file(
            case_name;
            data_dir,
            case_data_dir,

            dyn_plants_file,
            
            synchronous_machine_wt_loc_load, 
            synchronous_condenser_wt_loc_load,
            synchronous_machine,
            synchronous_condenser,

            plant_wt_loc_load,
            plant_no_gov_wt_loc_load,
            plant_no_gov,
            plant,

            avr_param,
            gov_param)
        
    end
    
    #--------------------------------------
    #--------------------------------------

    mpc_branch_column_select =
        String["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_bus_column_select =
        String["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_gen_column_select =
        String["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_scalar_column_select =
        String["mpc_baseMVA" ]

    dyn_gens_column_select =
        String["bus","sym_gen_type",
         "sym_gen_dynamic_para"]
    
    dyn_plants_column_select =
       String["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    # dyn_gens_data_types =
    #     [Int, Symbol, Symbol]

    # dyn_plants_data_types =
    #     [ Int, Symbol, Symbol, Bool, Symbol, Symbol]


    dyn_gens_data_types =
        [Int, String, String]

    dyn_plants_data_types =
        [ Int, String, String, Bool, String, String]
    
    #--------------------------------------

    list_data_types = [
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        dyn_gens_data_types,
        dyn_plants_data_types ]
    
    list_data_select = Union{String,Any}[
        [],
        mpc_branch_column_select,
        mpc_bus_column_select,
        mpc_gen_column_select,
        [],
        [],
        mpc_scalar_column_select,
        dyn_gens_column_select,
        dyn_plants_column_select ]

    list_csv_files =
        [ csv_branch_type_file,
          csv_branch_file,
          csv_bus_file,
          csv_gen_file,
          csv_gencost_file,
          csv_load_type_file,
          csv_scalar_file,
          dyn_gens_file,
          dyn_plants_file ]


    # list_df =
    #     [ a_data_select_list == [] ? CSV.read(a_csv_file, DataFrame) : CSV.read(a_csv_file, DataFrame; select = a_data_select_list)
    #       for (a_data_select_list, a_csv_file ) in
    #           zip(list_data_select, list_csv_files) ]
    

    # list_df =
    #     [( a_data_types_list == [] && a_data_select_list == []) ?
    #     DataFrame(CSV.File(a_csv_file )) : (a_data_types_list == [] && a_data_select_list != []) ?
    #     DataFrame(CSV.File(a_csv_file; select = a_data_select_list )) : (a_data_types_list != [] && a_data_select_list == []) ? DataFrame(CSV.File(a_csv_file; types  = a_data_types_list)) : DataFrame(CSV.File(a_csv_file; select = a_data_select_list, types = a_data_types_list))
    #       for ( a_data_types_list, a_data_select_list, a_csv_file ) in
    #           zip(list_data_types, list_data_select, list_csv_files)]


    list_df =
        [( a_data_types_list == [] && a_data_select_list == []) ?
        CSV.read(a_csv_file, DataFrame) : (
            a_data_types_list == [] && a_data_select_list != []) ?
        CSV.read(a_csv_file,DataFrame;select=a_data_select_list) :
        (a_data_types_list != [] && a_data_select_list == []) ?
        CSV.read(a_csv_file,
                 DataFrame;
                 types  = a_data_types_list) :
        CSV.read(a_csv_file,DataFrame;
                 select = a_data_select_list,
                 types = a_data_types_list)
          for (a_data_types_list,a_data_select_list,a_csv_file) in
              zip(list_data_types,list_data_select,list_csv_files)]
        
    list_sheet_names =
        [ "mpc_branch_type",
          "mpc_branch",
          "mpc_bus",
          "mpc_gen",
          "mpc_gencost",
          "mpc_load_type",
          "mpc_scalar",
          "dyn_gen",
          "dyn_plant"]
    
    @assert length(list_sheet_names) == length(list_df)

    XLSX.openxlsx(xlsx_file, mode="w") do xf
        for i in eachindex( list_sheet_names )
            sheet_name = list_sheet_names[i]
            df = list_df[i]

            if i == firstindex( list_sheet_names )
                sheet = xf[1]
                XLSX.rename!(sheet, sheet_name)
                XLSX.writetable!(sheet, df)
            else
                sheet = XLSX.addsheet!(xf, sheet_name)
                XLSX.writetable!(sheet, df)        
            end
        end
    end    

    return nothing
                
end



#---------------------------------------------------
#---------------------------------------------------


function create_default_static_net_json_data_by_xlsx(
    case_name;        
    data_dir = "",
    components_libs_dir = "",
    net_data_by_components_file = "",

    xlsx_data_file = "",
    
    by_components = true,
    wt_plants_data_types_bool = true )

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

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath( case_data_dir, "dyn")

    json_case_dir = joinpath( case_data_dir, "json")
        
    #--------------------------------------


    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
    
        
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     "net-default-static-data.json")
    else
    
        
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------


    if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
        
        xlsx_file =
            joinpath(mpc_data_dir,
                     "xlsx",
                     "net-default-static-data.xlsx")
    else
        
        xlsx_file =
            joinpath(mpc_data_dir,
                     "xlsx",
                     xlsx_data_file)
        
    end

    #--------------------------------------    

    if !( isfile( xlsx_file ) )

        create_a_default_case_net_data_xlsx_file(
            case_name,
            nothing;
            by_components =
                by_components)

    end
    
    #--------------------------------------    
    #--------------------------------------
        
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    	
    mpc_load_type_column_select =
        ["idx", "load_type"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------
    # xlsx
    #--------------------------------------

    mpc_branch_type_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_branch_type";
                       header = true,
                       infer_eltypes = true))


    mpc_gencost_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_gencost";
                       header = true,
                       infer_eltypes = true))
    
    mpc_branch_selected_data = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_branch";
                       header = true,
                       infer_eltypes = true)),
        mpc_branch_column_select)

    mpc_gen_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_gen";
                       header = true,
                       infer_eltypes = true)),
        mpc_gen_column_select)

    mpc_bus_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_bus";
                       header = true,
                       infer_eltypes = true)),
        mpc_bus_column_select)
    
    mpc_load_type_select_data = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_load_type";
                       header = true,
                       infer_eltypes = true)),
        mpc_load_type_column_select )

        
    mpc_scalar_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_scalar";
                       header = true,
                       infer_eltypes = true)),
        mpc_scalar_column_select)

    mpc_baseMVA = mpc_scalar_selected_data.mpc_baseMVA[1]


    if wt_plants_data_types_bool == true

        dyn_gens_data_types =
            [Int, Symbol, Symbol]

        dyn_plants_data_types =
            [Int, Symbol, Symbol, Bool, Symbol, Symbol]

        dyn_gens  = convert_dataframe_selected_cols_types(
            DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_gen";
                           header = true,
                           infer_eltypes = true)),
            dyn_gens_column_select) ,
                     dyn_gens_data_types,
                     dyn_gens_column_select)

        dyn_plants  = convert_dataframe_selected_cols_types(
           DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_plant";
                           header = true,
                           infer_eltypes = true)),
            dyn_plants_column_select) ,
                       dyn_plants_data_types,
                       dyn_plants_column_select)

        
    else

        dyn_gens_data_types =
            [String, Symbol, Symbol]

        dyn_plants_data_types =
            [String, Symbol, Symbol, Bool, Symbol, Symbol]
        
    
        dyn_gens  = convert_dataframe_selected_cols_types(
            DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_gen";
                           header = true)),
            dyn_gens_column_select) ,
                     dyn_gens_data_types,
                     dyn_gens_column_select)

        dyn_plants  = convert_dataframe_selected_cols_types(
           DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_plant";
                           header = true )),
            dyn_plants_column_select) ,
                       dyn_plants_data_types,
                       dyn_plants_column_select)

        # dyn_gens  =  DataFrames.select(
        #     DataFrame(
        #         XLSX.readtable(xlsx_file, "dyn_gen";
        #                    header = true ) ),
        #     dyn_gens_column_select...) 

        # dyn_plants  =  DataFrames.select(
        #     DataFrame(
        #         XLSX.readtable(xlsx_file, "dyn_plant";
        #                    header = true,
        #                    infer_eltypes = true)),
        #     dyn_plants_column_select)
        
    end

    #--------------------------------------    
    #--------------------------------------

    vec_edge_type =
        Symbol.(mpc_branch_type_data.branch_type)
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = mpc_baseMVA)

    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = mpc_baseMVA)
    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
            vec_edge_type =
                vec_edge_type )
    
    gen_nodes_dict_static_data =
        get_gen_nodes_dict_static_data_by_mpc(
            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA =
                mpc_baseMVA )

    gen_nodes_static_tup =
        get_gen_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA=mpc_baseMVA)

    load_nodes_static_tup_data =
        get_load_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA =
                mpc_baseMVA)

    load_nodes_idx_wt_type_tup =
        get_load_nodes_idx_wt_type_tup_by_mpc(
            mpc_load_type_select_data )

    transmission_nodes_static_tup_data =
        get_transmission_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA=
                mpc_baseMVA)

    transmission_nodes_exist =
        transmission_nodes_static_tup_data == nothing ?
         false : false

    loc_loads_idx_and_locP_locQ_data =
        get_loc_loads_idx_and_locP_locQ_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = mpc_baseMVA )

    loc_loads_exist =
        length(loc_loads_idx_and_locP_locQ_data) == 0 ?
        false : true

    #--------------------------------------
    # some indices
    #--------------------------------------

   (;multi_gens_dyn_pf_fun_kwd_net_idxs,
    multi_gens_dyn_pf_fun_kwd_n2s_idxs ) =
        NamedTupleTools.select(
            get_multi_net_nodes_idxs_wt_n2s(
                mpc_gen_selected_data,
                mpc_bus_selected_data),
            (:multi_gens_dyn_pf_fun_kwd_net_idxs,
             :multi_gens_dyn_pf_fun_kwd_n2s_idxs ))

    #--------------------------------------
    #--------------------------------------

    default_net_static_data_dict =
        OrderedDict(
            :loc_loads_exist => [loc_loads_exist],
            
            :transmission_nodes_exist =>
                [transmission_nodes_exist],
            
            :plant_generators =>
                gen_nodes_dict_static_data,
            
            :plant_loc_loads => get_static_loc_loads_json_data_storage_format(loc_loads_idx_and_locP_locQ_data),
            
            :plant_loads =>
                get_static_load_nodes_json_data_storage_format(load_nodes_static_tup_data,load_nodes_idx_wt_type_tup),

            :plant_transmission =>
                get_static_transmission_nodes_json_data_storage_format(transmission_nodes_static_tup_data, load_nodes_idx_wt_type_tup),
            
            :branches =>
                OrderedDict{
                    Union{Int64,String,Symbol},
                    NamedTuple}(a_branch.idx => a_branch
                                for a_branch in
                                    Dyn_Branches_data_list),
            :shunt =>
                [dict_shunt],

            :baseMVA =>
                [mpc_baseMVA],

            :gencost =>
                [dict_gencost] )

        json_net_static_data =
            JSON.json( default_net_static_data_dict )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_static_data )
        end

    #--------------------------------------
    
    return nothing
                
end


function create_default_static_net_json_data_by_mpc(
    case_name;        
    data_dir            = "",
    components_libs_dir = "",
    net_data_by_components_file = "",
    by_components       = true,
    wt_plants_data_types_bool = true )

        
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
                 "converted_data",
                 case_name )

    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath( case_data_dir, "dyn")

    json_case_dir = joinpath( case_data_dir, "json")

    #--------------------------------------        

    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
    
        
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     "net-default-static-data.json")
                      
    else
    
        
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------
    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    if  wt_plants_data_types_bool == true

        dyn_gens_data_types =
            [Int, Symbol, Symbol]

        dyn_plants_data_types =
            [Int, Symbol, Symbol, Bool, Symbol, Symbol]
        
    else

        dyn_gens_data_types =
            [String, Symbol, Symbol]

        dyn_plants_data_types =
            [String, Symbol, Symbol, Bool, Symbol, Symbol]
        
    end

    #--------------------------------------
    #--------------------------------------
    
    (;mpc_branch_selected_data,
     mpc_gen_selected_data,
     mpc_gencost_data,
     mpc_bus_selected_data,
     dyn_gens,
     dyn_plants,
     mpc_baseMVA,

     mpc_load_type_data,

     mpc_branch_type_data) =
         NamedTupleTools.select(
             get_case_data_by_csv(
                 case_name ;
                 case_data_dir =
                     case_data_dir,
                 mpc_branch_column_select =
                     mpc_branch_column_select,
                 mpc_gen_column_select =
                     mpc_gen_column_select,
                 mpc_bus_column_select =
                     mpc_bus_column_select,
                 mpc_scalar_column_select =
                     mpc_scalar_column_select,
                 dyn_gens_column_select =
                     dyn_gens_column_select,
                 dyn_plants_column_select =
                     dyn_plants_column_select,
                 dyn_gens_data_types =
                     dyn_gens_data_types,
                 dyn_plants_data_types =
                     dyn_plants_data_types,
                 wt_plants_data_types_bool =
                     wt_plants_data_types_bool ),
             (:mpc_branch_selected_data,
              :mpc_gen_selected_data,
              :mpc_gencost_data,
              :mpc_bus_selected_data,
              :dyn_gens,
              :dyn_plants,
              :mpc_baseMVA,

              :mpc_load_type_data,

              :mpc_branch_type_data))
            
    #--------------------------------------
    #--------------------------------------

    # vec_edge_type = [ ]

    vec_edge_type =
        Symbol.( mpc_branch_type_data.branch_type )
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = mpc_baseMVA)


    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = mpc_baseMVA) 
    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
            vec_edge_type =
                vec_edge_type )
    
    gen_nodes_dict_static_data =
        get_gen_nodes_dict_static_data_by_mpc(
            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA =
                mpc_baseMVA )

    gen_nodes_static_tup =
        get_gen_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA=mpc_baseMVA)

    load_nodes_static_tup_data =
        get_load_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA =
                mpc_baseMVA)

    load_nodes_idx_wt_type_tup =
        get_load_nodes_idx_wt_type_tup_by_mpc(
            mpc_load_type_data  )

    transmission_nodes_static_tup_data =
        get_transmission_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA=
                mpc_baseMVA)

    transmission_nodes_exist =
        transmission_nodes_static_tup_data == nothing ?
         false : true

    loc_loads_idx_and_locP_locQ_data =
        get_loc_loads_idx_and_locP_locQ_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = mpc_baseMVA )

    # @show loc_loads_idx_and_locP_locQ_data

    loc_loads_exist =
        length(loc_loads_idx_and_locP_locQ_data) == 0 ?
        false : true


    #--------------------------------------
    # some indices
    #--------------------------------------

   (;multi_gens_dyn_pf_fun_kwd_net_idxs,
    multi_gens_dyn_pf_fun_kwd_n2s_idxs ) =
        NamedTupleTools.select(
            get_multi_net_nodes_idxs_wt_n2s(
                mpc_gen_selected_data,
                mpc_bus_selected_data),
            (:multi_gens_dyn_pf_fun_kwd_net_idxs,
             :multi_gens_dyn_pf_fun_kwd_n2s_idxs ))

    #--------------------------------------
    #--------------------------------------

    default_net_static_data_dict =
        OrderedDict(
            :loc_loads_exist => [loc_loads_exist],
            
            :transmission_nodes_exist =>
                [transmission_nodes_exist],
            
            :plant_generators =>
                gen_nodes_dict_static_data,
            
            :plant_loc_loads => get_static_loc_loads_json_data_storage_format(loc_loads_idx_and_locP_locQ_data),
            :plant_loads =>
                get_static_load_nodes_json_data_storage_format(load_nodes_static_tup_data,load_nodes_idx_wt_type_tup),

            :plant_transmission =>
                get_static_transmission_nodes_json_data_storage_format(transmission_nodes_static_tup_data, load_nodes_idx_wt_type_tup),
            
            :branches =>
                OrderedDict{
                    Union{Int64,String,Symbol},
                    NamedTuple}(a_branch.idx => a_branch
                                for a_branch in
                                    Dyn_Branches_data_list),
            :shunt =>
                [dict_shunt],

            :baseMVA =>
                [mpc_baseMVA],

            :gencost =>
                [dict_gencost] )

    #--------------------------------------

        json_net_static_data =
            JSON.json( default_net_static_data_dict )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_static_data )

        end

    #--------------------------------------
    
    return nothing
            
    
end

#---------------------------------------------------
#---------------------------------------------------

function create_a_default_case_net_data_json_file(
    case_name;
    data_dir            = "",
    case_data_dir       = "",
    components_libs_dir = "",

    net_data_by_components_file = "",
    xlsx_data_file              = "",

    by_xlsx_bool = false )

    
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

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
            
        
    end
    
    #--------------------------------------
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")

    if !(isdir( json_case_dir))

        mkpath( json_case_dir)

    end
    
    #--------------------------------------

    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(
                json_case_dir,
                "net_data_by_components_file.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------
    
    if by_xlsx_bool == true

        if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
            
            xlsx_data_file =
                joinpath(case_data_dir,
                         "xlsx",
                         "net-static-data.xlsx" )
        else
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx", xlsx_data_file )
        end
        

        dict_net_data_by_components =
            get_net_data_by_components_by_xlsx(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir,
                xlsx_data_file =
                    xlsx_data_file )

        #--------------------------------------

        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components )

        end
        
    else

        #--------------------------------------

        dict_net_data_by_components =
            get_net_data_by_components_by_mpc(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir )

        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components)

        end
        
    end
    
    #--------------------------------------

    return nothing
    

end


# create_a_default_case_net_static_data_json_file

function create_a_default_case_net_data_json_file(
    case_name,
    nothing;
    data_dir = "",
    case_data_dir = "",
    dyn_data_dir = "",
    components_libs_dir = "",

    mpc_branch_column_select = "",
    mpc_branch_line_type =
        "PiModelLine",
    mpc_branch_transformer_type =
        "Transformer",

    mpc_bus_column_select = "",
    mpc_load_node_type =
        "PQ_Const_P",
    mpc_transmission_node_type =
        "Trans_t2_Node",

    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6", 
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine =
        "SM_2axis_cb_v6",
    synchronous_condenser =
        "SC_2axis_cb_v6",

    plant_wt_loc_load =
        "plant_wt_loc_load_v6",
    plant_no_gov_wt_loc_load =
        "plant_no_gov_wt_loc_load_v6",
    plant_no_gov =
        "plant_no_gov",
    plant =
        "plant_cb_v6",

    synchronous_machine_dynamic_param =
        "gen_dynamic_paras_ieee_14_b",
    avr_param = "avr_t1_cb_sauer__1_param",
    gov_param = "gov_t1_cb_sauer__1_param",
    
    net_data_by_components_file =
        nothing,
    xlsx_data_file =
        nothing,

    by_components =
        true,
    by_xlsx_bool =
        false )


    create_a_default_case_net_data_xlsx_file(
        case_name;
        data_dir,
        case_data_dir,
        dyn_data_dir,

        xlsx_default_data_file =
            xlsx_data_file,

        mpc_branch_column_select,
        mpc_branch_line_type,
        mpc_branch_transformer_type,

        mpc_bus_column_select,
        mpc_load_node_type,
        mpc_transmission_node_type,

        synchronous_machine_wt_loc_load,
        synchronous_condenser_wt_loc_load,
        synchronous_machine,
        synchronous_condenser,

        plant_wt_loc_load,
        plant_no_gov_wt_loc_load,
        plant_no_gov,
        plant,

        synchronous_machine_dynamic_param,
        avr_param,
        gov_param,
    
        by_components )
    
    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end
    
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

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
            
        
    end
    
    #--------------------------------------
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")


    if !(isdir( json_case_dir))

        mkpath( json_case_dir)

    end

    
    #--------------------------------------


    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     "net_data_by_components_file.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------
    
    if by_xlsx_bool == true

        if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx",
                    "net-static-data.xlsx" )
        else
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx",
                    xlsx_data_file )
        end
        
        dict_net_data_by_components =
            get_net_data_by_components_by_xlsx(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir,
                xlsx_data_file =
                    xlsx_data_file )


        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components )

        end
        
    else


        dict_net_data_by_components =
            get_net_data_by_components_by_mpc(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir )

        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components) 
            
        end
        
    end
    
    #--------------------------------------

    return nothing
    
    # return dict_net_data_by_components 
    

end

#---------------------------------------------------
#---------------------------------------------------


function create_a_default_static_case_net_data_json_by_xlsx(
    case_name;
    data_dir            = "",
    case_data_dir       = "",
    components_libs_dir = "",

    net_data_by_components_file = "",
    xlsx_data_file = "",

    by_xlsx_bool =
        false )


    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end
    
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

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
            
        
    end
    
    #--------------------------------------
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")

    if !(isdir( json_case_dir))

        mkpath( json_case_dir)

    end
    
    #--------------------------------------

    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(
                json_case_dir,
                "net-default-static-data.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------

    if (xlsx_data_file == nothing ||
        xlsx_data_file == "")

        xlsx_data_file =
            joinpath(case_data_dir,
                     "xlsx",
                     "net-static-data.xlsx" )
    else

        xlsx_data_file =
            joinpath(
                case_data_dir,
                "xlsx", xlsx_data_file )
    end

    dict_net_data_by_components =
        get_net_data_by_static_components_by_xlsx(
            ;case_name =
                case_name,        
            data_dir =
                data_dir,
            components_libs_dir =
                components_libs_dir,
            xlsx_data_file =
                xlsx_data_file )

    #--------------------------------------

    json_net_data_by_components =
        JSON.json( dict_net_data_by_components )

    # write

    open(json_net_data_by_components_file, "w") do io
        JSON3.pretty( io, json_net_data_by_components )

    end
    
    #--------------------------------------

    return nothing
    

end


function create_a_default_static_case_net_data_json_by_mpc(
    case_name;
    data_dir            = "",
    case_data_dir       = "",
    components_libs_dir = "",

    net_data_by_components_file = "" )


    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end
    
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

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
            
        
    end
    

    #--------------------------------------
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")

    if !(isdir( json_case_dir))

        mkpath( json_case_dir)

    end
    
    #--------------------------------------

    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(
                json_case_dir,
                "net-default-static-data.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end


    dict_net_data_by_components =
        get_net_data_by_static_components_by_mpc(
            ;case_name =
                case_name,        
            data_dir =
                data_dir,
            components_libs_dir =
                components_libs_dir )

    json_net_data_by_components =
        JSON.json( dict_net_data_by_components )

    # write

    open(json_net_data_by_components_file, "w") do io
        JSON3.pretty( io, json_net_data_by_components)

    end

    #--------------------------------------

    return nothing
    

end


#---------------------------------------------------
#---------------------------------------------------

function create_a_default_static_case_net_data_json_file(
    case_name;
    data_dir            = "",
    case_data_dir       = "",
    components_libs_dir = "",

    net_data_by_components_file =
        nothing,
    xlsx_data_file =
        nothing,

    by_xlsx_bool =
        false )


    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end
    
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

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
            
        
    end
    
    #--------------------------------------
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")

    if !(isdir( json_case_dir))

        mkpath( json_case_dir)

    end
    
    #--------------------------------------

    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(
                json_case_dir,
                "net-default-static-data.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------
    
    if by_xlsx_bool == true

        if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
            
            xlsx_data_file =
                joinpath(case_data_dir,
                         "xlsx",
                         "net-static-data.xlsx" )
        else
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx", xlsx_data_file )
        end
        

        dict_net_data_by_components =
            get_net_data_by_static_components_by_xlsx(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir,
                xlsx_data_file =
                    xlsx_data_file )

        #--------------------------------------

        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components )

        end
        
    else

        #--------------------------------------

        dict_net_data_by_components =
            get_net_data_by_static_components_by_mpc(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir )

        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components)

        end
        
    end
    
    #--------------------------------------

    return nothing
    

end


# create_a_default_case_net_static_data_json_file

function create_a_default_static_case_net_data_json_file(
    case_name,
    nothing;
    data_dir = "",
    case_data_dir = "",
    dyn_data_dir = "",
    components_libs_dir = "",

    mpc_branch_column_select = "",
    mpc_branch_line_type =
        "PiModelLine",
    mpc_branch_transformer_type =
        "Transformer",

    mpc_bus_column_select = "",
    mpc_load_node_type =
        "PQ_Const_P",
    mpc_transmission_node_type =
        "Trans_t2_Node",

    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6", 
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine =
        "SM_2axis_cb_v6",
    synchronous_condenser =
        "SC_2axis_cb_v6",

    plant_wt_loc_load =
        "plant_wt_loc_load_v6",
    plant_no_gov_wt_loc_load =
        "plant_no_gov_wt_loc_load_v6",
    plant_no_gov =
        "plant_no_gov",
    plant =
        "plant_cb_v6",

    synchronous_machine_dynamic_param =
        "gen_dynamic_paras_ieee_14_b",
    avr_param = "avr_t1_cb_sauer__1_param",
    gov_param = "gov_t1_cb_sauer__1_param",
    
    net_data_by_components_file =
        nothing,
    xlsx_data_file =
        nothing,

    by_components =
        true,
    by_xlsx_bool =
        false )


    create_a_default_case_net_data_xlsx_file(
        case_name;
        data_dir,
        case_data_dir,
        dyn_data_dir,

        xlsx_default_data_file =
            xlsx_data_file,

        mpc_branch_column_select,
        mpc_branch_line_type,
        mpc_branch_transformer_type,

        mpc_bus_column_select,
        mpc_load_node_type,
        mpc_transmission_node_type,

        synchronous_machine_wt_loc_load,
        synchronous_condenser_wt_loc_load,
        synchronous_machine,
        synchronous_condenser,

        plant_wt_loc_load,
        plant_no_gov_wt_loc_load,
        plant_no_gov,
        plant,

        synchronous_machine_dynamic_param,
        avr_param,
        gov_param,
    
        by_components )
    
    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end

    # #--------------------------------------
    
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

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
            
        
    end
    
    #--------------------------------------
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")


    if !(isdir( json_case_dir))

        mkpath( json_case_dir)

    end

    
    #--------------------------------------


    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     "net-default-static-data.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------
    
    if by_xlsx_bool == true

        if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx",
                    "net-static-data.xlsx" )
        else
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx",
                    xlsx_data_file )
        end
        
        dict_net_data_by_components =
            get_net_data_by_static_components_by_xlsx(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir,
                xlsx_data_file =
                    xlsx_data_file )


        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components )

        end
        
    else


        dict_net_data_by_components =
            get_net_data_by_static_components_by_mpc(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir )

        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components) 
            
        end
        
    end
    
    #--------------------------------------

    return nothing
    
    # return dict_net_data_by_components 
    

end

#---------------------------------------------------

function create_a_case_net_data_by_components_file(
    case_name;
    data_dir = "",
    components_libs_dir = "",
    net_data_by_components_file = "")

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
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")


    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")

    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     "net_data_by_components_file.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end
    
    #--------------------------------------
    #--------------------------------------
        
    dict_net_data_by_components =
        get_net_data_by_components_by_mpc(
            ;case_name =
                case_name,        
            data_dir =
                data_dir,
            components_libs_dir =
                components_libs_dir )
    
    json_net_data_by_components =
        JSON.json( dict_net_data_by_components )

    # write
    
    # open( net_data_by_components_file, "w") do io
    #     JSON3.pretty( io, json_net_data_by_components)
        
    # end

    open( json_net_data_by_components_file, "w") do io
        JSON3.pretty( io, json_net_data_by_components)
        
    end    
    
end

#---------------------------------------------------
#---------------------------------------------------

function create_a_case_net_data_by_components_file_by_xlsx(
    case_name;
    data_dir                    = "",
    components_libs_dir         = "",
    net_data_by_components_file = nothing,
    xlsx_data_file              = nothing )

    #--------------------------------------

    # if data_dir == ""
                
    #     data_dir =
    #         joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )

    # end
    

    # if   components_libs_dir == ""
            
    #     components_libs_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #              "components-lib")
    # end
    
    
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
    
    json_case_dir =
        joinpath(case_data_dir,
                 case_name,
                 "json")
        
    #--------------------------------------

    if (net_data_by_components_file == nothing ||
        net_data_by_components_file == "")
        
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     "net_data_by_components_file.json")
    else
        
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end
    
    
    #--------------------------------------
        
    dict_net_data_by_components =
        get_net_data_by_components_by_xlsx(
            ;case_name =
                case_name,        
            data_dir =
                data_dir,
            components_libs_dir =
                components_libs_dir,
            xlsx_data_file =
                xlsx_data_file )

    #--------------------------------------
    
    json_net_data_by_components =
        JSON.json( dict_net_data_by_components )

    # write
    
    open(json_net_data_by_components_file, "w") do io
        JSON3.pretty( io, json_net_data_by_components )
        
    end
    
end


#---------------------------------------------------
#---------------------------------------------------
# multi gens and plants
#---------------------------------------------------
#---------------------------------------------------

function check_multi_gens_bool_by_csv_file(
    csv_gen_file)
    
    mpc_gen_selected_data =
        CSV.File( csv_gen_file;
               select = ["bus"] )

    gens_idx = mpc_gen_selected_data.bus
    
    gens_countmap =
        countmap(gens_idx)

    multi_gens_nodes = [k for (k, v) in
                           gens_countmap if v > 1]

    multi_gen_bool =
        length(multi_gens_nodes) == 0 ? false : true

    return (;multi_gens_nodes,
            multi_gen_bool,
            gens_countmap)
    
end


function check_multi_gens_bool_by_case(
    case_name;
    data_dir = "",
    case_data_dir = "")

    #--------------------------------------
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    
    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name )
        
    end

    
    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")

    return check_multi_gens_bool_by_csv_file(
        csv_gen_file)

end


function get_node_idx_multi(idx)

    return parse(Int64,String(split("$(idx)",'-')[1])) 

end


function get_component_idx_in_multi( idx )

    return length(split("$(idx)",'-')) == 2 ?
        parse(Int64, String(split("$(idx)",'-')[2]) ) :
        parse(Int64, String(split("$(idx)",'-')[1]) )

end


function get_multi_gens_idx_wt_multi_gen_bool(
    mpc_gen;
    sorted_bool = true )

    if typeof(mpc_gen) == DataFrame

        multi_gens_idx_in_bus_cvs =
            sorted_bool == true ? sort(
                mpc_gen.bus ) : mpc_gen.bus
        
    else                

        multi_gens_idx_in_bus_cvs = sorted_bool == true ? 
            sort([ a_tup.bus for a_tup in
                      mpc_gen] ) : [ a_tup.bus for a_tup in
                      mpc_gen] 

    end
    

    gens_countmap =
        countmap(multi_gens_idx_in_bus_cvs)

    multi_gens_nodes = sorted_bool == true ?
        sort([k for (k, v) in
                  gens_countmap if v > 1] ) :
                      [k for (k, v) in
                           gens_countmap if v > 1]

    multi_gen_bool =
        length(multi_gens_nodes) == 0 ? false : true

    multi_gens_idx = [
        [ idxs ∈ multi_gens_nodes ?
            [ "$(idxs)-$(idx)"
              for idx in 1:gens_countmap[ idxs ]] :
                  "$(idxs)" for idxs in
                      unique(multi_gens_idx_in_bus_cvs)]...;]

    return (; multi_gens_idx,
            multi_gen_bool,
            multi_gens_nodes,
            
            multi_gens_idx_in_bus_cvs,
            gens_countmap)


end


function get_multi_gens_all_nodes_idx(
    multi_gens_idx,
    multi_gens_nodes,
    gens_countmap,
    all_nodes_idx)

    return  [[ idxs ∈ multi_gens_nodes ?
        [ "$(idxs)-$(idx)" for idx in
             1:gens_countmap[ idxs ]] :
        "$(idxs)"  for idxs in
            all_nodes_idx]...;]

end


function get_multi_net_nodes_idxs_wt_n2s(
    mpc_gen,
    mpc_bus)


    # if typeof(mpc_gen) == DataFrame

    #     # mpc_gen = Tables.rowtable(mpc_gen)
    #     mpc_gen = Tables.columntable(mpc_gen)

    # end

    # if typeof(mpc_bus) == DataFrame

    #     # mpc_bus = Tables.rowtable(mpc_bus)
    #     mpc_bus = Tables.columntable(mpc_bus)

    # end
    
    (;slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     nodes_with_demands_idx) =
         NamedTupleTools.select(
             get_net_nodes_type_idxs_by_mpc(
                 mpc_bus),
             (:slack_bus_idx,
             :gens_idx,
             :slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :loc_load_exist,
             :load_nodes_idx,
             :transmission_nodes_idx,
             :non_gens_nodes_idx,
             :all_nodes_idx,
             :non_slack_gens_and_non_gens_idx,
             :nodes_with_demands_idx))

    (;
     n2s_slack_gens_idx,         
     n2s_non_slack_gens_idx,     
     n2s_gens_idx,               
     n2s_non_gens_idx,           
     n2s_load_idx,               
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,      
     n2s_all_nodes_idx,
     n2s_nodes_with_demands_idx)  =
         NamedTupleTools.select(
             get_dict_n2s_streamlined_idx_by_mpc(
                 mpc_bus),
             (:n2s_slack_gens_idx,         
            :n2s_non_slack_gens_idx,     
            :n2s_gens_idx,               
            :n2s_non_gens_idx,           
            :n2s_load_idx,               
            :n2s_gens_with_loc_load_idxs,
            :n2s_transmission_idxs,      
            :n2s_all_nodes_idx,
            :n2s_nodes_with_demands_idx) )
    
    #--------------------------------------    
    #--------------------------------------
    
    (;multi_gens_idx,
     multi_gen_bool,
     multi_gens_nodes,

     multi_gens_idx_in_bus_cvs,
     gens_countmap) =
         NamedTupleTools.select(
             get_multi_gens_idx_wt_multi_gen_bool(
                 mpc_gen;
             sorted_bool = false),
             (:multi_gens_idx,
              :multi_gen_bool,
              :multi_gens_nodes,

              :multi_gens_idx_in_bus_cvs,
              :gens_countmap) )

    multi_gens_all_nodes_idx =
        [[ idxs ∈ multi_gens_nodes ?
        [ "$(idxs)-$(idx)" for idx in
             1:gens_countmap[ idxs ]] :
                 "$(idxs)"  for idxs in
                     all_nodes_idx]...;]
    
    n2s_multi_gens_idx =
            get_n2s_any( multi_gens_idx )

    n2s_multi_gens_all_nodes_idx =
        get_n2s_any( multi_gens_all_nodes_idx )


    multi_gens_dyn_pf_fun_kwd_net_idxs =
        (;slack_bus_idx,
         gens_idx,
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         loc_load_exist,
         load_nodes_idx,
         transmission_nodes_idx,
         non_gens_nodes_idx,
         all_nodes_idx,
         non_slack_gens_and_non_gens_idx,
         nodes_with_demands_idx,

         multi_gens_idx,
         multi_gen_bool,
         multi_gens_nodes,

         multi_gens_idx_in_bus_cvs,
         gens_countmap,
         multi_gens_all_nodes_idx)
    
    multi_gens_dyn_pf_fun_kwd_n2s_idxs =
        (;n2s_slack_gens_idx,         
         n2s_non_slack_gens_idx,     
         n2s_gens_idx,               
         n2s_non_gens_idx,           
         n2s_load_idx,               
         n2s_gens_with_loc_load_idxs,
         n2s_transmission_idxs,      
         n2s_all_nodes_idx,
         n2s_nodes_with_demands_idx,

         n2s_multi_gens_idx,
         n2s_multi_gens_all_nodes_idx)


    return (; multi_gens_dyn_pf_fun_kwd_net_idxs,
            multi_gens_dyn_pf_fun_kwd_n2s_idxs )
    

end


function create_a_default_case_dyn_multi_gens_file(
    case_name;
    data_dir = "",
    case_data_dir = "",
    
    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6", 
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine =
        "SM_2axis_cb_v6",
    synchronous_condenser =
        "SC_2axis_cb_v6",

    synchronous_machine_dynamic_param =
        "gen_dynamic_paras_ieee_14_b" )

    #--------------------------------------
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name)
        
    end

    #--------------------------------------
    
    dyn_data_dir  =
        joinpath(case_data_dir,
                 "dyn")

    if !(isdir( dyn_data_dir ))

        mkpath( dyn_data_dir )

    end

    #--------------------------------------
    
    dyn_gens_file =
        joinpath(dyn_data_dir,
                 "dyn_gen.csv")

    #--------------------------------------
    #--------------------------------------

    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")


    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")
    
    #--------------------------------------    

    mpc_bus_column_select =
        ["bus_i",
         "type",
         "Pd",
         "Qd"]
    
    mpc_gen_column_select =
        ["bus",
         "Pg",
         "Qg"]
    
    #--------------------------------------
    
    mpc_gen_selected_data =
        CSV.File( csv_gen_file;
                  select = mpc_gen_column_select )
    
    multi_gens_idx_in_bus_cvs =
        sort([ a_tup.bus for a_tup in
                  mpc_gen_selected_data])
    
    gens_countmap =
        countmap(multi_gens_idx_in_bus_cvs)

    # multi_gens_nodes =
    #     sort([k for (k, v) in
    #               gens_countmap if v > 1])

    multi_gens_nodes =
        [k for (k, v) in
                  gens_countmap if v > 1]

    multi_gen_bool =
        length(multi_gens_nodes) == 0 ?
        false : true

    #--------------------------------------

    mpc_bus_selected_data =
        CSV.File(csv_bus_file;
                 select = mpc_bus_column_select )
    
    nodes_idx_with_type =
        [ ( a_tup.bus_i, a_tup.type )
          for a_tup in
              mpc_bus_selected_data ]

    # nodes_idx_with_type =
    #     sort( nodes_idx_with_type, by = x->x[ 1 ] )

    nodes_idx_with_type =
        nodes_idx_with_type

    nodes_idx  = first.( nodes_idx_with_type  )

    nodes_type = second.( nodes_idx_with_type )
    
    #--------------------------------------
    
    gens_nodes_idx_in_bus_cvs =
        sort([ a_tup.bus_i for a_tup in
             mpc_bus_selected_data if
                 a_tup.type == 3 || a_tup.type == 2 ])

    multi_gens_idx = multi_gen_bool == false ?
        gens_nodes_idx_in_bus_cvs :
        [[ idxs ∈ multi_gens_nodes ?
        [ "$(idxs)-$(idx)" for idx in
             1:gens_countmap[ idxs ]] :
        "$(idxs)"  for idxs in
            unique(multi_gens_idx_in_bus_cvs)]...;]
    
    slack_nodes_idx_in_bus_cvs =
        [a_tup.bus_i for a_tup in
             mpc_bus_selected_data if
                 a_tup.type == 3 ]
    
    #--------------------------------------
    
    bus_nodes_with_demands_idx =
        [ a_tup.bus_i for a_tup in
             mpc_bus_selected_data if (
                 a_tup.Pd != 0 || a_tup.Pd != 0.0 ) || (
                 a_tup.Qd != 0 || a_tup.Qd != 0.0 ) ]

    #--------------------------------------

    bus_SM =
        [ a_tup.bus for a_tup in
             mpc_gen_selected_data if (
                 (a_tup.bus ∈ slack_nodes_idx_in_bus_cvs) ||
                     ( a_tup.Pg != 0 || a_tup.Pg != 0.0))]

    bus_SC =
        [ a_tup.bus for a_tup in
             mpc_gen_selected_data if (
                 (a_tup.bus ∉ slack_nodes_idx_in_bus_cvs) &&
                     ( a_tup.Pg == 0 || a_tup.Pg == 0.0 ))]
    
    #--------------------------------------

    bus_SM_in_loc_load =
        [ idx for idx in bus_SM
             if idx ∈  bus_nodes_with_demands_idx ]

    bus_SC_in_loc_load = length(bus_SC) == 0 ? [] :
        [ idx for idx in bus_SC
             if idx ∈  bus_nodes_with_demands_idx ]

    #--------------------------------------
    
    # gens_idx = sort(gens_idx)
    
    sym_gen_type =
        [get_node_idx_multi(idx) ∈ bus_SC_in_loc_load ?
        synchronous_condenser_wt_loc_load :
        get_node_idx_multi(idx) ∈ bus_SM_in_loc_load ?
        synchronous_machine_wt_loc_load :
        get_node_idx_multi(idx) ∈ bus_SC ?
        synchronous_condenser :
        synchronous_machine
         for idx in multi_gens_idx ]
    
    sym_gen_dynamic_para =
        [synchronous_machine_dynamic_param
         for idx in multi_gens_idx  ]

    #--------------------------------------
    
    CSV.write( dyn_gens_file,
               DataFrame(
                   bus = multi_gens_idx,
                   sym_gen_type = sym_gen_type,
                   sym_gen_dynamic_para =
                       sym_gen_dynamic_para  ) )

    return nothing



end


function create_a_default_case_dyn_multi_plants_file(
    case_name;
    data_dir      = "",
    case_data_dir = "",

    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6", 
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine =
        "SM_2axis_cb_v6",
    synchronous_condenser =
        "SC_2axis_cb_v6",

    plant_wt_loc_load =
        "plant_wt_loc_load_v6",
    plant_no_gov_wt_loc_load =
        "plant_no_gov_wt_loc_load_v6",
    plant_no_gov = "plant_no_gov",
    plant = "plant_cb_v6",
    
    avr_param = "avr_t1_cb_sauer__1_param",
    gov_param = "gov_t1_cb_sauer__1_param" )

    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end

    #--------------------------------------
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name)
        
    end

    #--------------------------------------
    
    dyn_data_dir  =
        joinpath(case_data_dir,
                 "dyn")


    if !(isdir( dyn_data_dir ))

        mkpath( dyn_data_dir )

    end

    #--------------------------------------

    dyn_plants_file =
        joinpath(dyn_data_dir,
                 "dyn_plant.csv")

    #--------------------------------------
    #--------------------------------------

    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")


    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")
    
    #--------------------------------------    

    mpc_bus_column_select =
        ["bus_i",
         "type",
         "Pd",
         "Qd"]
    
    mpc_gen_column_select =
        ["bus",
         "Pg",
         "Qg"]
    
    #--------------------------------------
    
    mpc_gen_selected_data =
        CSV.File( csv_gen_file;
                  select = mpc_gen_column_select )
    
    multi_gens_idx_in_bus_cvs =
        sort([ a_tup.bus for a_tup in
                  mpc_gen_selected_data])
    
    gens_countmap =
        countmap(multi_gens_idx_in_bus_cvs)

    multi_gens_nodes =
        sort([k for (k, v) in
                  gens_countmap if v > 1])

    multi_gen_bool =
        length(multi_gens_nodes) == 0 ? false : true

    #--------------------------------------

    mpc_bus_selected_data =
        CSV.File(csv_bus_file;
                 select = mpc_bus_column_select )
    
    nodes_idx_with_type =
        [ ( a_tup.bus_i, a_tup.type )
          for a_tup in
              mpc_bus_selected_data ]

    nodes_idx_with_type =
        sort( nodes_idx_with_type, by = x->x[ 1 ] )

    nodes_idx  = first.( nodes_idx_with_type  )

    nodes_type = second.( nodes_idx_with_type )
    
    #--------------------------------------
    
    gens_nodes_idx_in_bus_cvs =
        sort([ a_tup.bus_i for a_tup in
             mpc_bus_selected_data if
                 a_tup.type == 3 || a_tup.type == 2 ])

    multi_gens_idx = multi_gen_bool == false ?
        gens_nodes_idx_in_bus_cvs :
        [[ idxs ∈ multi_gens_nodes ?
        [ "$(idxs)-$(idx)" for idx in
             1:gens_countmap[ idxs ]] :
        "$(idxs)"  for idxs in
            unique(multi_gens_idx_in_bus_cvs)]...;]
    
    slack_nodes_idx_in_bus_cvs =
        sort([a_tup.bus_i for a_tup in
             mpc_bus_selected_data if
                 a_tup.type == 3 ])
    
    #--------------------------------------
    
    bus_nodes_with_demands_idx =
        [ a_tup.bus_i for a_tup in
             mpc_bus_selected_data if (
                 a_tup.Pd != 0 || a_tup.Pd != 0.0 ) || (
                 a_tup.Qd != 0 || a_tup.Qd != 0.0 ) ]

    #--------------------------------------

    bus_SM =
        sort([ a_tup.bus for a_tup in
             mpc_gen_selected_data if (
                 (a_tup.bus ∈ slack_nodes_idx_in_bus_cvs) ||
                     ( a_tup.Pg != 0 || a_tup.Pg != 0.0 )) ])

    bus_SC =
        sort([ a_tup.bus for a_tup in
             mpc_gen_selected_data if (
                 (a_tup.bus ∉ slack_nodes_idx_in_bus_cvs) &&
                     ( a_tup.Pg == 0 || a_tup.Pg == 0.0 )) ])
    
    #--------------------------------------

    bus_SM_in_loc_load =
        [ idx for idx in bus_SM
             if idx ∈  bus_nodes_with_demands_idx ]

    bus_SC_in_loc_load = length(bus_SC) == 0 ? [] :
        [ idx for idx in bus_SC
             if idx ∈  bus_nodes_with_demands_idx ]
    
    #--------------------------------------
    
    # gens_idx = sort(gens_idx)

    Plant_type =
        [get_node_idx_multi(idx) ∈ bus_SC_in_loc_load ?
        plant_no_gov_wt_loc_load :
        get_node_idx_multi(idx) ∈ bus_SM_in_loc_load ?
        plant_wt_loc_load  :
        get_node_idx_multi(idx) ∈ bus_SC ?
        plant_no_gov :
        plant
          for idx in multi_gens_idx ]
    
    Gen =
        [get_node_idx_multi(idx) ∈ bus_SC_in_loc_load ?
        synchronous_condenser_wt_loc_load :
        get_node_idx_multi(idx) ∈ bus_SM_in_loc_load ?
        synchronous_machine_wt_loc_load :
        get_node_idx_multi(idx) ∈ bus_SC ?
        synchronous_condenser :
        synchronous_machine
          for idx in multi_gens_idx ]


    # isa_slack = [ a_type == 3 ? true : false
    #               for (idx, a_type) in
    #                   zip( nodes_idx, nodes_type ) if
    #                       idx ∈ gens_idx ]

    isa_slack =
        [get_node_idx_multi(idx) ∈ slack_nodes_idx_in_bus_cvs ?
        true : false  for idx in multi_gens_idx  ]

    Gov = [get_node_idx_multi(idx) ∈ bus_SC ?
        :nothing : gov_param
            for idx in multi_gens_idx ]

    Exc = [avr_param for idx in multi_gens_idx ]
    	
    #--------------------------------------
    
    CSV.write( dyn_plants_file,
               DataFrame(bus = multi_gens_idx,
                         Plant_type = Plant_type,
                         Gen = Gen,
                         isa_slack = isa_slack,
                         Gov = Gov,
                         Exc = Exc ) )

    return nothing

    
end



function create_a_default_multi_gens_case_net_data_xlsx_file(
    case_name;
    data_dir = "",
    case_data_dir = "",
    dyn_data_dir = "",

    mpc_branch_column_select = "",
    mpc_branch_line_type =
        "PiModelLine",
    mpc_branch_transformer_type =
        "Transformer",

    mpc_bus_column_select = "",
    mpc_load_node_type =
        "PQ_Const_P",
    mpc_transmission_node_type =
        "Trans_t2_Node",

    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6", 
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine =
        "SM_2axis_cb_v6",
    synchronous_condenser =
        "SC_2axis_cb_v6",

    plant_wt_loc_load =
        "plant_wt_loc_load_v6",
    plant_no_gov_wt_loc_load =
        "plant_no_gov_wt_loc_load_v6",
    plant_no_gov =
        "plant_no_gov",
    plant =
        "plant_cb_v6",

    synchronous_machine_dynamic_param =
        "gen_dynamic_paras_ieee_14_b",
    avr_param = "avr_t1_cb_sauer__1_param",
    gov_param = "gov_t1_cb_sauer__1_param",
    
    by_components = true )

    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end
    
    #--------------------------------------
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name)
        
    end

    #--------------------------------------

    
    xlsx_data_dir  =
        joinpath(case_data_dir,
                 "xlsx")

    if !(isdir( xlsx_data_dir ))

        mkpath( xlsx_data_dir)

    end

    #--------------------------------------
    
    xlsx_file =
        joinpath(xlsx_data_dir,
                 "net-static-data.xlsx")
    
    #--------------------------------------

    csv_branch_file =
        joinpath(case_data_dir,
                 "mpc_branch.csv")

    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")

    csv_gencost_file =
        joinpath(case_data_dir,
                 "mpc_gencost.csv")
    
    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")

    csv_scalar_file =
        joinpath(case_data_dir,
                 "mpc_scalar.csv")

    #--------------------------------------
    #--------------------------------------

    csv_load_type_file  = 
        joinpath(case_data_dir,
                 "mpc_load_type.csv")

    csv_branch_type_file =
        joinpath(case_data_dir,
                 "mpc_branch_type.csv")
    

    if !( isfile( csv_load_type_file ) )
        
        create_a_default_case_mpc_load_type(
            case_name;
            data_dir      = data_dir,
            case_data_dir = case_data_dir ,

            mpc_bus_column_select =
                mpc_bus_column_select,
            mpc_load_node_type =
                mpc_load_node_type,
            mpc_transmission_node_type =
                mpc_transmission_node_type )
    end
    

    if !( isfile( csv_branch_type_file ) )
        
        create_a_default_case_mpc_branch_type(
            case_name;
            data_dir      = data_dir,
            case_data_dir = case_data_dir,

            mpc_branch_column_select =
                mpc_branch_column_select,
            mpc_branch_line_type =
                mpc_branch_line_type ,
            mpc_branch_transformer_type =
                mpc_branch_transformer_type  )
    end
    
    #--------------------------------------
    
    if dyn_data_dir == ""
    
        dyn_data_dir  =
            joinpath(case_data_dir,
                     "dyn")
        
    end


    if !(isdir( dyn_data_dir ))

        mkpath( dyn_data_dir )

    end

    
    dyn_gens_file =
        joinpath(dyn_data_dir,
                 "dyn_gen.csv")

    dyn_plants_file =
        joinpath(dyn_data_dir,
                 "dyn_plant.csv")
    

    if !( isfile( dyn_gens_file ) )
        
        create_a_default_case_dyn_multi_gens_file(
            case_name;
            data_dir,
            case_data_dir,

            synchronous_machine_wt_loc_load, 
            synchronous_condenser_wt_loc_load,
            synchronous_machine,
            synchronous_condenser,

            synchronous_machine_dynamic_param)
    end


    if !( isfile( dyn_plants_file ) )
        
        create_a_default_case_dyn_multi_plants_file(
            case_name;
            data_dir,
            case_data_dir,

            synchronous_machine_wt_loc_load, 
            synchronous_condenser_wt_loc_load,
            synchronous_machine,
            synchronous_condenser,

            plant_wt_loc_load,
            plant_no_gov_wt_loc_load,
            plant_no_gov,
            plant,

            avr_param,
            gov_param)
        
    end
    

    #--------------------------------------
    #--------------------------------------


    mpc_branch_column_select =
        String["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_bus_column_select =
        String["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_gen_column_select =
        String["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_scalar_column_select =
        String["mpc_baseMVA" ]

    dyn_gens_column_select =
        String["bus","sym_gen_type",
         "sym_gen_dynamic_para"]
    
    dyn_plants_column_select =
       String["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    # dyn_gens_data_types =
    #     [Int, Symbol, Symbol]

    # dyn_plants_data_types =
    #     [ Int, Symbol, Symbol, Bool, Symbol, Symbol]


    dyn_gens_data_types =
        [String, String, String]

    dyn_plants_data_types =
        [String, String, String, Bool, String, String]
    
    #--------------------------------------

    list_data_types = [
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        dyn_gens_data_types,
        dyn_plants_data_types ]
    
    list_data_select = Union{String,Any}[
        [],
        mpc_branch_column_select,
        mpc_bus_column_select,
        mpc_gen_column_select,
        [],
        [],
        mpc_scalar_column_select,
        dyn_gens_column_select,
        dyn_plants_column_select ]

    list_csv_files =
        [ csv_branch_type_file,
          csv_branch_file,
          csv_bus_file,
          csv_gen_file,
          csv_gencost_file,
          csv_load_type_file,
          csv_scalar_file,
          dyn_gens_file,
          dyn_plants_file ]


    list_df =
        [( a_data_types_list == [] && a_data_select_list == []) ?
        CSV.read(a_csv_file, DataFrame) :
        (a_data_types_list == [] && a_data_select_list != []) ?
        CSV.read(a_csv_file, DataFrame;
                 select = a_data_select_list ) :
        (a_data_types_list != [] && a_data_select_list == []) ?
        CSV.read(a_csv_file, DataFrame;
                 types  = a_data_types_list) :
                     CSV.read(a_csv_file, DataFrame;
                              select = a_data_select_list,
                              types = a_data_types_list)
         for (a_data_types_list,a_data_select_list,a_csv_file) in
             
             zip(list_data_types,
                 list_data_select,list_csv_files)]
        
    list_sheet_names =
        [ "mpc_branch_type",
          "mpc_branch",
          "mpc_bus",
          "mpc_gen",
          "mpc_gencost",
          "mpc_load_type",
          "mpc_scalar",
          "dyn_gen",
          "dyn_plant"]
    
    @assert length(list_sheet_names) == length(list_df)

    XLSX.openxlsx(xlsx_file, mode="w") do xf
        for i in eachindex( list_sheet_names )
            sheet_name = list_sheet_names[i]
            df = list_df[i]

            if i == firstindex( list_sheet_names )
                sheet = xf[1]
                XLSX.rename!(sheet, sheet_name)
                XLSX.writetable!(sheet, df)
            else
                sheet = XLSX.addsheet!(xf, sheet_name)
                XLSX.writetable!(sheet, df)        
            end
        end
    end    

    return nothing
                
end


function create_a_default_multi_gens_case_net_data_xlsx_file(
    case_name,
    nothing;
    data_dir = "",
    case_data_dir = "",
    dyn_data_dir = "",

    mpc_branch_column_select = "",
    mpc_branch_line_type =
        "PiModelLine",
    mpc_branch_transformer_type =
        "Transformer",

    mpc_bus_column_select = "",
    mpc_load_node_type =
        "PQ_Const_P",
    mpc_transmission_node_type =
        "Trans_t2_Node",

    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6", 
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine =
        "SM_2axis_cb_v6",
    synchronous_condenser =
        "SC_2axis_cb_v6",

    plant_wt_loc_load =
        "plant_wt_loc_load_v6",
    plant_no_gov_wt_loc_load =
        "plant_no_gov_wt_loc_load_v6",
    plant_no_gov =
        "plant_no_gov",
    plant =
        "plant_cb_v6",

    synchronous_machine_dynamic_param =
        "gen_dynamic_paras_ieee_14_b",
    avr_param =
        "avr_t1_cb_sauer__1_param",
    gov_param =
        "gov_t1_cb_sauer__1_param",
    
    by_components = true)

    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                     "data")
        
    end

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name)
        
    end

    
    #--------------------------------------
    
    xlsx_data_dir  =
        joinpath(case_data_dir,
                 "xlsx")

    if !(isdir( xlsx_data_dir ))

        mkpath( xlsx_data_dir)

    end

    #--------------------------------------
    
    xlsx_file =
        joinpath(xlsx_data_dir,
                 "net-static-data.xlsx")
    
    #--------------------------------------

    csv_branch_file =
        joinpath(case_data_dir,
                 "mpc_branch.csv")

    csv_gen_file =
        joinpath(case_data_dir,
                 "mpc_gen.csv")

    csv_gencost_file =
        joinpath(case_data_dir,
                 "mpc_gencost.csv")
    
    csv_bus_file =
        joinpath(case_data_dir,
                 "mpc_bus.csv")

    csv_scalar_file =
        joinpath(case_data_dir,
                 "mpc_scalar.csv")

    #--------------------------------------
    #--------------------------------------

    csv_load_type_file  = 
        joinpath(case_data_dir,
                 "mpc_load_type.csv")

    csv_branch_type_file =
        joinpath(case_data_dir,
                 "mpc_branch_type.csv")
    

    if !( isfile( csv_load_type_file ) )
        
        create_a_default_case_mpc_load_type(
            case_name;
            data_dir =
                data_dir,
            case_data_dir =
                case_data_dir,

            mpc_bus_column_select =
                mpc_bus_column_select,
            mpc_load_node_type =
                mpc_load_node_type,
            mpc_transmission_node_type =
                mpc_transmission_node_type )
    end
    

    if !( isfile( csv_branch_type_file ) )
        
        create_a_default_case_mpc_branch_type(
            case_name;
            data_dir =
                data_dir,
            case_data_dir =
                case_data_dir,

            mpc_branch_column_select =
                mpc_branch_column_select,
            mpc_branch_line_type =
                mpc_branch_line_type ,
            mpc_branch_transformer_type =
                mpc_branch_transformer_type  )
    end
    
    #--------------------------------------
    
    if dyn_data_dir == ""
    
        dyn_data_dir  =
            joinpath(case_data_dir,
                     "dyn")
        
    end


    if !(isdir( dyn_data_dir ))

        mkpath( dyn_data_dir )

    end

    
    dyn_gens_file =
        joinpath(dyn_data_dir,
                 "dyn_gen.csv")

    dyn_plants_file =
        joinpath(dyn_data_dir,
                 "dyn_plant.csv")
    

    if !( isfile( dyn_gens_file ) )
        
        create_a_default_case_dyn_multi_gens_file(
            case_name;
            data_dir,
            case_data_dir,

            synchronous_machine_wt_loc_load, 
            synchronous_condenser_wt_loc_load,
            synchronous_machine,
            synchronous_condenser,

            synchronous_machine_dynamic_param)
    end


    if !( isfile( dyn_plants_file ) )
        
        create_a_default_case_dyn_multi_plants_file(
            case_name;
            data_dir,
            case_data_dir,

            synchronous_machine_wt_loc_load, 
            synchronous_condenser_wt_loc_load,
            synchronous_machine,
            synchronous_condenser,

            plant_wt_loc_load,
            plant_no_gov_wt_loc_load,
            plant_no_gov,
            plant,

            avr_param,
            gov_param)
        
    end
    
    #--------------------------------------
    #--------------------------------------

    mpc_branch_column_select =
        String["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_bus_column_select =
        String["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_gen_column_select =
        String["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_scalar_column_select =
        String["mpc_baseMVA" ]

    dyn_gens_column_select =
        String["bus","sym_gen_type",
         "sym_gen_dynamic_para"]
    
    dyn_plants_column_select =
       String["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    # dyn_gens_data_types =
    #     [Int, Symbol, Symbol]

    # dyn_plants_data_types =
    #     [ Int, Symbol, Symbol, Bool, Symbol, Symbol]


    dyn_gens_data_types =
        [String, String, String]

    dyn_plants_data_types =
        [ String, String, String, Bool, String, String]
    
    #--------------------------------------

    list_data_types = [
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        dyn_gens_data_types,
        dyn_plants_data_types ]
    
    list_data_select = Union{String,Any}[
        [],
        mpc_branch_column_select,
        mpc_bus_column_select,
        mpc_gen_column_select,
        [],
        [],
        mpc_scalar_column_select,
        dyn_gens_column_select,
        dyn_plants_column_select ]

    list_csv_files =
        [ csv_branch_type_file,
          csv_branch_file,
          csv_bus_file,
          csv_gen_file,
          csv_gencost_file,
          csv_load_type_file,
          csv_scalar_file,
          dyn_gens_file,
          dyn_plants_file ]

    list_df =
        [( a_data_types_list == [] && a_data_select_list == []) ?
        CSV.read(a_csv_file, DataFrame) : (
            a_data_types_list == [] &&
                a_data_select_list != []) ?
                CSV.read(a_csv_file,
                         DataFrame;
                         select=a_data_select_list) :
        (a_data_types_list != [] && a_data_select_list == []) ?
        CSV.read(a_csv_file,
                 DataFrame;
                 types  = a_data_types_list) :
        CSV.read(a_csv_file,DataFrame;
                 select = a_data_select_list,
                 types = a_data_types_list)
         for (a_data_types_list,a_data_select_list,a_csv_file) in
             zip(list_data_types,
                 list_data_select,list_csv_files)]
        
    list_sheet_names =
        [ "mpc_branch_type",
          "mpc_branch",
          "mpc_bus",
          "mpc_gen",
          "mpc_gencost",
          "mpc_load_type",
          "mpc_scalar",
          "dyn_gen",
          "dyn_plant"]
    
    @assert length(list_sheet_names) == length(list_df)

    XLSX.openxlsx(xlsx_file, mode="w") do xf
        for i in eachindex( list_sheet_names )
            sheet_name = list_sheet_names[i]
            df = list_df[i]

            if i == firstindex( list_sheet_names )
                sheet = xf[1]
                XLSX.rename!(sheet, sheet_name)
                XLSX.writetable!(sheet, df)
            else
                sheet = XLSX.addsheet!(xf, sheet_name)
                XLSX.writetable!(sheet, df)        
            end
        end
    end    

    return nothing
                
end



function create_a_default_multi_gens_case_net_data_json_file(
    case_name;
    data_dir            = "",
    case_data_dir       = "",
    components_libs_dir = "",

    net_data_by_components_file = nothing,
    xlsx_data_file              = nothing,

    by_xlsx_bool = false )

    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end

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

    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name)
        
    end

    #--------------------------------------
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")

    if !(isdir( json_case_dir))

        mkpath( json_case_dir)

    end
    
    #--------------------------------------

    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)    
        
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                  "multi-net_data_by_components_file.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------
    
    if by_xlsx_bool == true

        if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx",
                    "net-static-data.xlsx" )
        else
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx", xlsx_data_file )
        end
        

        if !( isfile( xlsx_data_file ) )
            
            create_a_default_multi_gens_case_net_data_xlsx_file(
                case_name;
                data_dir,
                case_data_dir,
                dyn_data_dir,

                mpc_branch_column_select,
                mpc_branch_line_type,
                mpc_branch_transformer_type,

                mpc_bus_column_select,
                mpc_load_node_type,
                mpc_transmission_node_type,

                synchronous_machine_wt_loc_load,
                synchronous_condenser_wt_loc_load,
                synchronous_machine,
                synchronous_condenser,

                plant_wt_loc_load,
                plant_no_gov_wt_loc_load,
                plant_no_gov,
                plant,

                synchronous_machine_dynamic_param,
                avr_param,
                gov_param,

                by_components )
            
            
        end

        dict_net_data_by_components =
            get_multi_gens_net_data_by_components_by_xlsx(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir,
                xlsx_data_file =
                    xlsx_data_file )

        #--------------------------------------

        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components )

        end
        
    else

        dict_net_data_by_components =
            get_multi_gens_net_data_by_components_by_mpc(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir )

        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components)

        end
        
    end
    
    #--------------------------------------

    return nothing
    

end


function create_a_default_multi_gens_case_net_data_json_file(
    case_name,
    nothing;
    data_dir = "",
    case_data_dir = "",
    dyn_data_dir = "",
    components_libs_dir = "",

    mpc_branch_column_select = "",
    mpc_branch_line_type =
        "PiModelLine",
    mpc_branch_transformer_type =
        "Transformer",

    mpc_bus_column_select = "",
    mpc_load_node_type =
        "PQ_Const_P",
    mpc_transmission_node_type =
        "Trans_t2_Node",

    synchronous_machine_wt_loc_load =
        "SM_2axis_wt_loc_load_cb_v6", 
    synchronous_condenser_wt_loc_load =
        "SC_2axis_wt_loc_load_cb_v6",
    synchronous_machine =
        "SM_2axis_cb_v6",
    synchronous_condenser =
        "SC_2axis_cb_v6",

    plant_wt_loc_load =
        "plant_wt_loc_load_v6",
    plant_no_gov_wt_loc_load =
        "plant_no_gov_wt_loc_load_v6",
    plant_no_gov =
        "plant_no_gov",
    plant =
        "plant_cb_v6",

    synchronous_machine_dynamic_param =
        "gen_dynamic_paras_ieee_14_b",
    avr_param = "avr_t1_cb_sauer__1_param",
    gov_param = "gov_t1_cb_sauer__1_param",
    
    net_data_by_components_file = nothing,
    xlsx_data_file =
        nothing,

    by_components = true,
    by_xlsx_bool = false )

    
    #--------------------------------------
    
    # if case_data_dir == "" || data_dir == ""
        
    #     case_data_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "data-dir",
    #                  "converted_data",
    #                  case_name )
        
    # end

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
    
    if (case_data_dir == "") || (case_data_dir == nothing)

        package_dir = pkgdir(ePowerSim)

        data_dir = joinpath(package_dir,
                            "data")
        
        case_data_dir =
            joinpath(data_dir,
                     "converted-data",
                     case_name)
        
    end

    #--------------------------------------
    
    json_case_dir =
        joinpath(case_data_dir,
                 "json")


    if !(isdir( json_case_dir))

        mkpath( json_case_dir)

    end

    
    #--------------------------------------


    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                  "multi-net-data-by-components-file.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end

    #--------------------------------------
    
    if by_xlsx_bool == true

        if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx",
                    "net-static-data.xlsx" )
        else
            
            xlsx_data_file =
                joinpath(
                    case_data_dir,
                    "xlsx",
                    xlsx_data_file )
        end


        if !( isfile( xlsx_data_file ) )

            create_a_default_multi_gens_case_net_data_xlsx_file(
                case_name,
                nothing;
                data_dir,
                case_data_dir,
                dyn_data_dir,

                mpc_branch_column_select,
                mpc_branch_line_type,
                mpc_branch_transformer_type,

                mpc_bus_column_select,
                mpc_load_node_type,
                mpc_transmission_node_type,

                synchronous_machine_wt_loc_load,
                synchronous_condenser_wt_loc_load,
                synchronous_machine,
                synchronous_condenser,

                plant_wt_loc_load,
                plant_no_gov_wt_loc_load,
                plant_no_gov,
                plant,

                synchronous_machine_dynamic_param,
                avr_param,
                gov_param,

                by_components )
            
            
        end
        
        dict_net_data_by_components =
            get_multi_gens_net_data_by_components_by_xlsx(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir,
                xlsx_data_file =
                    xlsx_data_file )


        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty( io, json_net_data_by_components )

        end
        
    else

        dict_net_data_by_components =
            get_multi_gens_net_data_by_components_by_mpc(
                ;case_name =
                    case_name,        
                data_dir =
                    data_dir,
                components_libs_dir =
                    components_libs_dir )

        json_net_data_by_components =
            JSON.json( dict_net_data_by_components )

        # write

        open(json_net_data_by_components_file, "w") do io
            JSON3.pretty(io, json_net_data_by_components) 
            
        end
        
    end
    
    #--------------------------------------

    return nothing
    
    # return dict_net_data_by_components 
    

end


#---------------------------------------------------
#---------------------------------------------------
#  Reading model data from json net file
#---------------------------------------------------
#---------------------------------------------------

"""
    get_net_data_by_components_from_json_file(
        net_data_by_components_file;
        in_components_type_sym = false )

Returns namedtuples of dynamic network data `plant_generators_data_from_json`, `plant_loads_data_from_json`, `plant_transmission_data_from_json`, `edge_data_from_json`, `shunt_data_from_json`, `baseMVA_data_from_json`, `gencost_data_from_json`.


# Arguments
- `net_data_by_components_file::String`: the path to network json file.
- `in_components_type_sym::Bool=false`: the variable that determines how plants data are stored in json network file.

"""
function get_net_data_by_components_from_json_file(
    net_data_by_components_file;
    in_components_type_sym = false )

    Nodes_Branches_data_by_components = convert(
        OrderedDict{Symbol, Vector},
        copy(JSON3.read(
            net_data_by_components_file)))
     
    plant_generators_data_from_json =
        get_a_gen_plant_data_json_to_nt.(
            [ namedtuple(a_plant_data)
              for a_plant_data in
                  convert(Vector{Dict},
                Nodes_Branches_data_by_components[
                    :plant_generators])] ;
            in_components_type_sym =
                in_components_type_sym )

    plant_loads_data_from_json =
        get_a_non_gen_plant_data_json_to_nt.(
            [ namedtuple(a_plant_data)
              for a_plant_data in
                  convert(Vector{Dict},
                Nodes_Branches_data_by_components[
                    :plant_loads])] ;
            in_components_type_sym =
                in_components_type_sym )

    plant_transmission_data_from_json =
        get_a_non_gen_plant_data_json_to_nt.(
            [ namedtuple(a_plant_data)
              for a_plant_data in
                  convert(Vector{Dict},
                Nodes_Branches_data_by_components[
                    :plant_transmissions])] ;
            in_components_type_sym =
                in_components_type_sym )

    edge_data_from_json =
        get_an_edge_data_json_to_nt.(
            [ namedtuple(an_edge_data)
              for an_edge_data in
                  convert(Vector{Dict},
                Nodes_Branches_data_by_components[
                    :branches ])] ;
            in_components_type_sym =
                in_components_type_sym )

    shunt_data_from_json =
        namedtuple(convert(Vector{Dict},
                Nodes_Branches_data_by_components[
                    :shunt ])[1] )

    baseMVA_data_from_json =
        convert(Vector{Int64},
                Nodes_Branches_data_by_components[
                    :baseMVA])[1]

    gencost_data_from_json =
        namedtuple(convert(Vector{Dict},
                Nodes_Branches_data_by_components[
                    :gencost ])[1] )
    
    return (;
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json,
            baseMVA_data_from_json,
            gencost_data_from_json)
    
end

#---------------------------------------------------

"""
This function `get_Dyn_Nodes_Branches_data_from_json_file`
is not fully implented and tested.
 
"""
function get_Dyn_Nodes_Branches_data_from_json_file(
    json_Dyn_Nodes_Branches_data_file)

    #--------------------------------------
    
    Nodes_Branches_data_from_json_file =
        JSON3.read(
            json_Dyn_Nodes_Branches_data_file)
        
    Dyn_Nodes_data_from_json =
        Nodes_Branches_data_from_json_file[:Nodes]

    Dyn_Branches_data_from_json =
        Nodes_Branches_data_from_json_file[:Branches]

    dyn_branches_data = OrderedDict{
        Symbol, NamedTuple}(
            a_key => (edge_type =
                Symbol(namedtuple(data).edge_type),
                      edge_data = namedtuple(
                          namedtuple(data).edge_data))
            for (a_key, data) in
                pairs(Dyn_Branches_data_from_json))


    dyn_nodes_data = OrderedDict{
        Symbol, NamedTuple}(
            a_key =>
                (edge_type = Symbol(
                    namedtuple(data).edge_type),
                 edge_data = namedtuple(
                     namedtuple(data).edge_data) )
            for (a_key, data) in
                pairs(Dyn_Nodes_data_from_json))
    
    Nodes_Branches_data_from_json_file =
        copy(JSON3.read(
            json_Dyn_Nodes_Branches_data_file))
    
    convert(OrderedDict{Symbol, OrderedDict{
        Symbol, NamedTuple} },
            Nodes_Branches_data_from_json_file)

    OrderedDict{
        Symbol, NamedTuple}(
            a_key => namedtuple(data)
            for (a_key, data) in
                pairs(Dyn_Nodes_data_from_json))

    
    convert(OrderedDict{
        Symbol, NamedTuple}, Dyn_Nodes_data_from_json)
    
    #--------------------------------------

    return Dict{Symbol, NamedTuple}(
        a_key => namedtuple(a_dict)
        for (a_key,a_dict)  in
            convert(Dict{Symbol, Dict{Symbol,Float64}},
                    nt_params_from_file))

    #--------------------------------------
    
end


#---------------------------------------------------
#---------------------------------------------------
#  Reading model dynamic data from csv or xlsx
#---------------------------------------------------
#---------------------------------------------------


function get_net_data_by_components_by_xlsx(
    ;case_name          = "case9",        
    data_dir            = "",
    components_libs_dir = "",
    by_components       = true,
    xlsx_data_file      = nothing )

    #--------------------------------------

    # if data_dir == ""
        
    #     data_dir = joinpath(@__DIR__,"..","..","src",
    #                  "data-dir","converted_data" )
        
    # end

    # #--------------------------------------

    # if components_libs_dir == ""

    #     components_libs_dir =
    #         joinpath(@__DIR__,"..","..","src",
    #                  "components-lib" )
    # end
    
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
    
    # case_data_dir =
    #     joinpath(data_dir,
    #              case_name )
        
    case_data_dir =
        joinpath(data_dir,
                 "converted-data",
                 case_name)

    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath( case_data_dir, "dyn")

    json_data_dir = joinpath( case_data_dir, "json")


    #--------------------------------------
    
    if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
        
        xlsx_file =
            joinpath(case_data_dir, "xlsx",
                     "net-static-data.xlsx")
    else
        
        xlsx_file =
            joinpath(case_data_dir, "xlsx",
                     xlsx_data_file)
        
    end

    #--------------------------------------
    #--------------------------------------

    (;dict_plants_gen_sym_type,

     dict_gen_sym_type,
     dict_gens_dyn_nt_params,

     dict_gov_sym_type,
     dict_gov_nt_params,

     dict_avr_sym_type,
     dict_avr_nt_params) =
         NamedTupleTools.select(
             get_dynamic_components_parameters_libs_from_json(
                 ;components_libs_dir =
                     components_libs_dir),
             (:dict_plants_gen_sym_type,

              :dict_gen_sym_type,
              :dict_gens_dyn_nt_params,

              :dict_gov_sym_type,
              :dict_gov_nt_params,

              :dict_avr_sym_type,
              :dict_avr_nt_params ) )

    #--------------------------------------
    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    dyn_gens_data_types =
        [Int, Symbol, Symbol]

    dyn_plants_data_types =
        [Int, Symbol, Symbol, Bool, Symbol, Symbol]

    #--------------------------------------

    #--------------------------------------
    # xlsx
    #--------------------------------------

    mpc_branch_type_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_branch_type";
                       header = true,
                       infer_eltypes = true))


    mpc_gencost_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_gencost";
                       header = true,
                       infer_eltypes = true))
    
    mpc_branch_selected_data = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_branch";
                       header = true,
                       infer_eltypes = true)),
        mpc_branch_column_select)

    mpc_gen_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_gen";
                       header = true,
                       infer_eltypes = true)),
        mpc_gen_column_select)

    mpc_bus_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_bus";
                       header = true,
                       infer_eltypes = true)),
        mpc_bus_column_select)

    mpc_scalar_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_scalar";
                       header = true,
                       infer_eltypes = true)),
        mpc_scalar_column_select)

    mpc_baseMVA = mpc_scalar_selected_data.mpc_baseMVA[1]

    
    dyn_gens  = convert_dataframe_selected_cols_types(
        DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "dyn_gen";
                       header = true,
                       infer_eltypes = true)),
        dyn_gens_column_select),
                 dyn_gens_data_types,
                 dyn_gens_column_select)

    
    dyn_plants  = convert_dataframe_selected_cols_types(
       DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "dyn_plant";
                       header = true,
                       infer_eltypes = true)),
        dyn_plants_column_select) ,
                   dyn_plants_data_types,
                   dyn_plants_column_select)

    
    #--------------------------------------
    #--------------------------------------

    # vec_edge_type = [ ]

    vec_edge_type =
        Symbol.(mpc_branch_type_data.branch_type)
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = 1.0)

    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = 1.0)
    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
             vec_edge_type =  vec_edge_type )

    Dyn_Nodes_data_by_components_dict =
        get_Dyn_Nodes_data_by_components_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus_selected_data,
            mpc_gen_selected_data ;
            mpc_baseMVA=1.0,

            p_order = 1.0,
            v_ref = 1.0,
            ω_ref = ωs,

            load_type =
                PQ_Const_I,
            load_node_type =
                plant_PQ_Const_I,

            transmission_type =
                Trans_t2_Node,
            transmission_node_type =
                plant_Transmission_t2,
            by_components =
                by_components)
    
    return  OrderedDict(
        :plant_generators =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_generators],
        :plant_loads =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_loads],
        :plant_transmissions =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_transmissions],            

        :branches => Dyn_Branches_data_list,
        
        :shunt => [dict_shunt],
        
        :baseMVA => [mpc_baseMVA],
        
        :gencost => [dict_gencost] )
    
end

#---------------------------------------------------

function get_net_data_by_components_by_mpc(
    ;case_name          = "case9",        
    data_dir            = "",
    components_libs_dir = "",
    by_components       = true )


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

    # case_data_dir =
    #     joinpath(data_dir,
    #              case_name )
    
    case_data_dir =
        joinpath(data_dir,
                 "converted-data",
                 case_name )
    

    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  =
        joinpath( case_data_dir, "dyn")

    json_data_dir =
        joinpath( case_data_dir, "json")
            
    #--------------------------------------
    #--------------------------------------

    (;dict_plants_gen_sym_type,

     dict_gen_sym_type,
     dict_gens_dyn_nt_params,

     dict_gov_sym_type,
     dict_gov_nt_params,

     dict_avr_sym_type,
     dict_avr_nt_params) =
         NamedTupleTools.select(
             get_dynamic_components_parameters_libs_from_json(
                 ;components_libs_dir =
                     components_libs_dir),
             (:dict_plants_gen_sym_type,

              :dict_gen_sym_type,
              :dict_gens_dyn_nt_params,

              :dict_gov_sym_type,
              :dict_gov_nt_params,

              :dict_avr_sym_type,
              :dict_avr_nt_params ) )
    
    #--------------------------------------
    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    dyn_gens_data_types =
        [Int, Symbol, Symbol]

    dyn_plants_data_types =
        [Int, Symbol, Symbol, Bool, Symbol, Symbol]

    #--------------------------------------
    #--------------------------------------
    
    (;mpc_branch_selected_data,
     mpc_gen_selected_data,
     mpc_gencost_data,
     mpc_bus_selected_data,
     dyn_gens,
     dyn_plants,
     mpc_baseMVA,

     mpc_branch_type_data) =
         NamedTupleTools.select(
             get_case_data_by_csv(
                 case_name ;
                 case_data_dir =
                     case_data_dir,
                 mpc_branch_column_select =
                     mpc_branch_column_select,
                 mpc_gen_column_select =
                     mpc_gen_column_select,
                 mpc_bus_column_select =
                     mpc_bus_column_select,
                 mpc_scalar_column_select =
                     mpc_scalar_column_select,
                 dyn_gens_column_select =
                     dyn_gens_column_select,
                 dyn_plants_column_select =
                     dyn_plants_column_select,
                 dyn_gens_data_types =
                     dyn_gens_data_types,
                 dyn_plants_data_types =
                     dyn_plants_data_types ),
             (:mpc_branch_selected_data,
              :mpc_gen_selected_data,
              :mpc_gencost_data,
              :mpc_bus_selected_data,
              :dyn_gens,
              :dyn_plants,
              :mpc_baseMVA,

              :mpc_branch_type_data))

    #--------------------------------------

    # vec_edge_type = [ ]

    vec_edge_type =
        Symbol.( mpc_branch_type_data.branch_type )
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = 1.0)


    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = 1.0) 
    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
             vec_edge_type =  vec_edge_type )

    Dyn_Nodes_data_by_components_dict =
        get_Dyn_Nodes_data_by_components_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus_selected_data,
            mpc_gen_selected_data ;
            mpc_baseMVA=1.0,

            p_order = 1.0,
            v_ref = 1.0,
            ω_ref = ωs,

            load_type =
                PQ_Const_I,
            load_node_type =
                plant_PQ_Const_I,

            transmission_type =
                Trans_t2_Node,
            transmission_node_type =
                plant_Transmission_t2,
            by_components =
                by_components)
    
    return  OrderedDict(
        :plant_generators =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_generators],
        :plant_loads =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_loads],
        :plant_transmissions =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_transmissions],       

        :branches => Dyn_Branches_data_list,
        
        :shunt => [dict_shunt],
        
        :baseMVA => [mpc_baseMVA],

        :gencost => [dict_gencost])
    
end


#---------------------------------------------------
#---------------------------------------------------
#  Reading model static data from csv or xlsx
#---------------------------------------------------
#---------------------------------------------------


function get_net_data_by_static_components_by_xlsx(
    ;case_name          = "case9",        
    data_dir            = "",
    components_libs_dir = "",
    by_components       = true,
    xlsx_data_file      = "" )

    #--------------------------------------
    
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

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath( case_data_dir, "dyn")

    json_data_dir = joinpath( case_data_dir, "json")


    #--------------------------------------
    
    if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
        
        xlsx_file =
            joinpath(case_data_dir, "xlsx",
                     "net-static-data.xlsx")
    else
        
        xlsx_file =
            joinpath(case_data_dir, "xlsx",
                     xlsx_data_file)
        
    end

    #--------------------------------------
    #--------------------------------------
    
    dict_plants_gen_sym_type =
        get_abstract_type_dict_subtypes(
            SdGenPlant )

    dict_gen_sym_type =
        get_abstract_type_dict_subtypes(
            SdGen )

    #--------------------------------------
    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    dyn_gens_data_types =
        [Int, Symbol, Symbol]

    dyn_plants_data_types =
        [Int, Symbol, Symbol, Bool, Symbol, Symbol]

    #--------------------------------------
    # xlsx
    #--------------------------------------

    mpc_branch_type_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_branch_type";
                       header = true,
                       infer_eltypes = true))


    mpc_gencost_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_gencost";
                       header = true,
                       infer_eltypes = true))
    
    mpc_branch_selected_data = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_branch";
                       header = true,
                       infer_eltypes = true)),
        mpc_branch_column_select)

    mpc_gen_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_gen";
                       header = true,
                       infer_eltypes = true)),
        mpc_gen_column_select)

    mpc_bus_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_bus";
                       header = true,
                       infer_eltypes = true)),
        mpc_bus_column_select)

    mpc_scalar_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_scalar";
                       header = true,
                       infer_eltypes = true)),
        mpc_scalar_column_select)

    mpc_baseMVA = mpc_scalar_selected_data.mpc_baseMVA[1]

    
    dyn_gens  = convert_dataframe_selected_cols_types(
        DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "dyn_gen";
                       header = true,
                       infer_eltypes = true)),
        dyn_gens_column_select),
                 dyn_gens_data_types,
                 dyn_gens_column_select)

    
    dyn_plants  = convert_dataframe_selected_cols_types(
       DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "dyn_plant";
                       header = true,
                       infer_eltypes = true)),
        dyn_plants_column_select) ,
                   dyn_plants_data_types,
                   dyn_plants_column_select)

    
    #--------------------------------------
    #--------------------------------------

    vec_edge_type =
        Symbol.(mpc_branch_type_data.branch_type)
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = 1.0)

    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = 1.0)
    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
             vec_edge_type =  vec_edge_type )

    static_Nodes_data_by_components_dict =
        get_static_Nodes_data_by_components_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,

            dyn_gens,
            dyn_plants,

            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA =
                mpc_baseMVA,

            load_type =
                PQ_Const_I,
            load_node_type =
                plant_PQ_Const_I,

            transmission_type =
                Trans_t2_Node,
            transmission_node_type =
                plant_Transmission_t2,
            by_components =
                by_components)
    
    return  OrderedDict(
        :plant_generators =>
                static_Nodes_data_by_components_dict[
                    :plant_generators],
        :plant_loads =>
                static_Nodes_data_by_components_dict[
                    :plant_loads],
        :plant_transmissions =>
                static_Nodes_data_by_components_dict[
                    :plant_transmissions],            

        :branches => Dyn_Branches_data_list,
        
        :shunt => [dict_shunt],
        
        :baseMVA => [mpc_baseMVA],
        
        :gencost => [dict_gencost] )
    
end

#---------------------------------------------------

function get_net_data_by_static_components_by_mpc(
    ;case_name          = "case9",        
    data_dir            = "",
    components_libs_dir = "",
    by_components       = true )


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

    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath( case_data_dir, "dyn")

    json_data_dir = joinpath( case_data_dir, "json")
            
    #--------------------------------------
    #--------------------------------------
    
    dict_plants_gen_sym_type =
        get_abstract_type_dict_subtypes(
            SdGenPlant )

    dict_gen_sym_type =
        get_abstract_type_dict_subtypes(
            SdGen )
    
    #--------------------------------------
    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    dyn_gens_data_types =
        [Int, Symbol, Symbol]

    dyn_plants_data_types =
        [Int, Symbol, Symbol, Bool, Symbol, Symbol]

    #--------------------------------------
    #--------------------------------------
    
    (;mpc_branch_selected_data,
     mpc_gen_selected_data,
     mpc_gencost_data,
     mpc_bus_selected_data,
     dyn_gens,
     dyn_plants,
     mpc_baseMVA,

     mpc_branch_type_data) =
         NamedTupleTools.select(
             get_case_data_by_csv(
                 case_name ;
                 case_data_dir =
                     case_data_dir,
                 mpc_branch_column_select =
                     mpc_branch_column_select,
                 mpc_gen_column_select =
                     mpc_gen_column_select,
                 mpc_bus_column_select =
                     mpc_bus_column_select,
                 mpc_scalar_column_select =
                     mpc_scalar_column_select,
                 dyn_gens_column_select =
                     dyn_gens_column_select,
                 dyn_plants_column_select =
                     dyn_plants_column_select,
                 dyn_gens_data_types =
                     dyn_gens_data_types,
                 dyn_plants_data_types =
                     dyn_plants_data_types ),
             (:mpc_branch_selected_data,
              :mpc_gen_selected_data,
              :mpc_gencost_data,
              :mpc_bus_selected_data,
              :dyn_gens,
              :dyn_plants,
              :mpc_baseMVA,

              :mpc_branch_type_data))

    #--------------------------------------

    vec_edge_type =
        Symbol.( mpc_branch_type_data.branch_type )
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = 1.0)


    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = 1.0) 
    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
             vec_edge_type =  vec_edge_type )

    static_Nodes_data_by_components_dict =
        get_static_Nodes_data_by_components_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,

            dyn_gens,
            dyn_plants,

            mpc_bus_selected_data,
            mpc_gen_selected_data;
            
            mpc_baseMVA =
                mpc_baseMVA,

            load_type =
                PQ_Const_I,
            load_node_type =
                plant_PQ_Const_I,

            transmission_type =
                Trans_t2_Node,
            transmission_node_type =
                plant_Transmission_t2,
            by_components =
                by_components )
    
    return  OrderedDict(
        :plant_generators =>
                static_Nodes_data_by_components_dict[
                    :plant_generators],
        :plant_loads =>
                static_Nodes_data_by_components_dict[
                    :plant_loads],
        :plant_transmissions =>
                static_Nodes_data_by_components_dict[
                    :plant_transmissions],       

        :branches => Dyn_Branches_data_list,
        
        :shunt => [dict_shunt],
        
        :baseMVA => [mpc_baseMVA],

        :gencost => [dict_gencost])
    
end


#---------------------------------------------------
# multi gens and plants
#---------------------------------------------------

function get_multi_gens_net_data_by_components_by_xlsx(
    ;case_name          = "case9",        
    data_dir            = "",
    components_libs_dir = "",
    by_components       = true,
    xlsx_data_file      = nothing,
    wt_plants_data_types_bool = false)


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

    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath( case_data_dir, "dyn")

    json_data_dir = joinpath( case_data_dir, "json")

    #--------------------------------------
    

    if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
        
        xlsx_file =
            joinpath(mpc_data_dir, "xlsx",
                     "net-static-data.xlsx")
    else
        
        xlsx_file =
            joinpath(mpc_data_dir, "xlsx",
                     xlsx_data_file)
        
    end

    #--------------------------------------
    #--------------------------------------

    (;dict_plants_gen_sym_type,

     dict_gen_sym_type,
     dict_gens_dyn_nt_params,

     dict_gov_sym_type,
     dict_gov_nt_params,

     dict_avr_sym_type,
     dict_avr_nt_params) =
         NamedTupleTools.select(
           get_dynamic_components_parameters_libs_from_json(
                 ;components_libs_dir =
                     components_libs_dir),
             (:dict_plants_gen_sym_type,

              :dict_gen_sym_type,
              :dict_gens_dyn_nt_params,

              :dict_gov_sym_type,
              :dict_gov_nt_params,

              :dict_avr_sym_type,
              :dict_avr_nt_params ) )

    #--------------------------------------
    # xlsx
    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]
    
    #--------------------------------------
    #--------------------------------------

    mpc_branch_type_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_branch_type";
                       header = true,
                       infer_eltypes = true))


    mpc_gencost_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_gencost";
                       header = true,
                       infer_eltypes = true))
    
    mpc_branch_selected_data = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_branch";
                       header = true,
                       infer_eltypes = true)),
        mpc_branch_column_select)

    mpc_gen_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_gen";
                       header = true,
                       infer_eltypes = true)),
        mpc_gen_column_select)

    mpc_bus_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_bus";
                       header = true,
                       infer_eltypes = true)),
        mpc_bus_column_select)

    mpc_scalar_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_scalar";
                       header = true,
                       infer_eltypes = true)),
        mpc_scalar_column_select)

    mpc_baseMVA = mpc_scalar_selected_data.mpc_baseMVA[1]
    
    #--------------------------------------

    # vec_edge_type = [ ]

    vec_edge_type =
        Symbol.(mpc_branch_type_data.branch_type)
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = 1.0)

    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = 1.0)

    #--------------------------------------
    #--------------------------------------
    
    if wt_plants_data_types_bool == true

        dyn_gens_data_types =
            [Int, Symbol, Symbol]

        dyn_plants_data_types =
            [Int, Symbol, Symbol, Bool, Symbol, Symbol]
        
        dyn_gens  = convert_dataframe_selected_cols_types(
            DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_gen";
                           header = true,
                           infer_eltypes = true)),
            dyn_gens_column_select),
                     dyn_gens_data_types,
                     dyn_gens_column_select)


        dyn_plants  = convert_dataframe_selected_cols_types(
           DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_plant";
                           header = true,
                           infer_eltypes = true)),
            dyn_plants_column_select) ,
                       dyn_plants_data_types,
                       dyn_plants_column_select)
        
    else

        dyn_gens_data_types =
            [String, Symbol, Symbol]

        dyn_plants_data_types =
            [String, Symbol, Symbol, Bool, Symbol, Symbol]
        
        dyn_gens  = convert_dataframe_selected_cols_types(
            DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_gen";
                           header = true)),
            dyn_gens_column_select),
                     dyn_gens_data_types,
                     dyn_gens_column_select)


        dyn_plants  = convert_dataframe_selected_cols_types(
           DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_plant";
                           header = true)),
            dyn_plants_column_select) ,
                       dyn_plants_data_types,
                       dyn_plants_column_select)
        
    end

    #--------------------------------------    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
             vec_edge_type =  vec_edge_type )

    Dyn_Nodes_data_by_components_dict =
        get_multi_gens_Dyn_Nodes_data_by_components_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA=1.0,

            p_order = 1.0,
            v_ref = 1.0,
            ω_ref = ωs,

            load_type =
                PQ_Const_I,
            load_node_type =
                plant_PQ_Const_I,

            transmission_type =
                Trans_t2_Node,
            transmission_node_type =
                plant_Transmission_t2,
            by_components =
                by_components)
    
    return  OrderedDict(
        :plant_generators =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_generators],
        :plant_loads =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_loads],
        :plant_transmissions =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_transmissions],            

        :branches => Dyn_Branches_data_list,
        
        :shunt => [dict_shunt],
        
        :baseMVA => [mpc_baseMVA],
        
        :gencost => [dict_gencost] )
    
end


function get_multi_gens_net_data_by_components_by_mpc(
    ;case_name          = "case9",        
    data_dir            = "",
    components_libs_dir = "",
    by_components       = true,
    wt_plants_data_types_bool = false)


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
    

    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath( case_data_dir, "dyn")

    json_data_dir = joinpath( case_data_dir, "json")
            
    #--------------------------------------
    #--------------------------------------

    (;dict_plants_gen_sym_type,

     dict_gen_sym_type,
     dict_gens_dyn_nt_params,

     dict_gov_sym_type,
     dict_gov_nt_params,

     dict_avr_sym_type,
     dict_avr_nt_params) =
         NamedTupleTools.select(
           get_dynamic_components_parameters_libs_from_json(
                 ;components_libs_dir =
                     components_libs_dir),
             (:dict_plants_gen_sym_type,

              :dict_gen_sym_type,
              :dict_gens_dyn_nt_params,

              :dict_gov_sym_type,
              :dict_gov_nt_params,

              :dict_avr_sym_type,
              :dict_avr_nt_params ) )
    
    #--------------------------------------
    #--------------------------------------
        
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------
    #--------------------------------------
    
    if wt_plants_data_types_bool == true

        dyn_gens_data_types =
            [Int, Symbol, Symbol]

        dyn_plants_data_types =
            [Int, Symbol, Symbol, Bool, Symbol, Symbol]
        
    else

        dyn_gens_data_types =
            [String, Symbol, Symbol]

        dyn_plants_data_types =
            [String, Symbol, Symbol, Bool, Symbol, Symbol]
        
    end

    #--------------------------------------
    #--------------------------------------
    
    (;mpc_branch_selected_data,
     mpc_gen_selected_data,
     mpc_gencost_data,
     mpc_bus_selected_data,
     dyn_gens,
     dyn_plants,
     mpc_baseMVA,

     mpc_branch_type_data) =
         NamedTupleTools.select(
             get_case_data_by_csv(
                 case_name;
                 case_data_dir =
                     case_data_dir,
                 mpc_branch_column_select =
                     mpc_branch_column_select,
                 mpc_gen_column_select =
                     mpc_gen_column_select,
                 mpc_bus_column_select =
                     mpc_bus_column_select,
                 mpc_scalar_column_select =
                     mpc_scalar_column_select,
                 dyn_gens_column_select =
                     dyn_gens_column_select,
                 dyn_plants_column_select =
                     dyn_plants_column_select,
                 dyn_gens_data_types =
                     dyn_gens_data_types,
                 dyn_plants_data_types =
                     dyn_plants_data_types,
                 wt_plants_data_types_bool =
                     wt_plants_data_types_bool),
             (:mpc_branch_selected_data,
              :mpc_gen_selected_data,
              :mpc_gencost_data,
              :mpc_bus_selected_data,
              :dyn_gens,
              :dyn_plants,
              :mpc_baseMVA,

              :mpc_branch_type_data))

    #--------------------------------------
    
    # vec_edge_type = [ ]

    vec_edge_type =
        Symbol.( mpc_branch_type_data.branch_type )
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = mpc_baseMVA)


    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = mpc_baseMVA) 
    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
             vec_edge_type =  vec_edge_type )

    Dyn_Nodes_data_by_components_dict =
        get_multi_gens_Dyn_Nodes_data_by_components_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus_selected_data,
            mpc_gen_selected_data ;
            mpc_baseMVA=1.0,

            p_order = 1.0,
            v_ref = 1.0,
            ω_ref = ωs,

            load_type =
                PQ_Const_I,
            load_node_type =
                plant_PQ_Const_I,

            transmission_type =
                Trans_t2_Node,
            transmission_node_type =
                plant_Transmission_t2,
            by_components =
                by_components)
    
    return  OrderedDict(
        :plant_generators =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_generators],
        :plant_loads =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_loads],
        :plant_transmissions =>
                Dyn_Nodes_data_by_components_dict[
                    :plant_transmissions],            

        :branches => Dyn_Branches_data_list,
        
        :shunt => [dict_shunt],
        
        :baseMVA => [mpc_baseMVA],

        :gencost => [dict_gencost])
    
end

#---------------------------------------------------
#---------------------------------------------------
#  Reading model static data from csv or xlsx
#---------------------------------------------------
#---------------------------------------------------

function get_default_static_net_json_data_by_xlsx(
    case_name;        
    data_dir            = "",
    components_libs_dir = "",
    by_components       = true,
    xlsx_data_file      = nothing,
    wt_plants_data_types_bool = false )


    #--------------------------------------

    # if data_dir == ""
        
    #     data_dir = joinpath(@__DIR__,"..","..","src",
    #                  "data-dir","converted_data" )
        
    # end

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
    
    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath( case_data_dir, "dyn")

    json_data_dir = joinpath( case_data_dir, "json")

    #--------------------------------------
        
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------
    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath(case_data_dir, "dyn")

    json_data_dir = joinpath(case_data_dir, "json")

    #--------------------------------------
    #--------------------------------------


    if (xlsx_data_file == nothing ||
        xlsx_data_file == "")
        
        xlsx_file =
            joinpath(mpc_data_dir,
                     "xlsx",
                     "net-static-data.xlsx")
    else
        
        xlsx_file =
            joinpath(mpc_data_dir,
                     "xlsx",
                     xlsx_data_file)
        
    end
    
    #--------------------------------------
    # xlsx
    #--------------------------------------

    mpc_branch_type_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_branch_type";
                       header = true,
                       infer_eltypes = true))


    mpc_gencost_data  =  DataFrame(
        XLSX.readtable(xlsx_file, "mpc_gencost";
                       header = true,
                       infer_eltypes = true))
    
    mpc_branch_selected_data = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_branch";
                       header = true,
                       infer_eltypes = true)),
        mpc_branch_column_select)

    mpc_gen_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_gen";
                       header = true,
                       infer_eltypes = true)),
        mpc_gen_column_select)

    mpc_bus_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_bus";
                       header = true,
                       infer_eltypes = true)),
        mpc_bus_column_select)

    mpc_scalar_selected_data  = DataFrames.select(
        DataFrame(
            XLSX.readtable(xlsx_file, "mpc_scalar";
                       header = true,
                       infer_eltypes = true)),
        mpc_scalar_column_select)

    mpc_baseMVA = mpc_scalar_selected_data.mpc_baseMVA[1]


    if wt_plants_data_types_bool == true

        dyn_gens_data_types =
            [Int, Symbol, Symbol]

        dyn_plants_data_types =
            [Int, Symbol, Symbol, Bool, Symbol, Symbol]

        dyn_gens  = convert_dataframe_selected_cols_types(
            DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_gen";
                           header = true,
                           infer_eltypes = true)),
            dyn_gens_column_select) ,
                     dyn_gens_data_types,
                     dyn_gens_column_select)

        dyn_plants  = convert_dataframe_selected_cols_types(
           DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_plant";
                           header = true,
                           infer_eltypes = true)),
            dyn_plants_column_select) ,
                       dyn_plants_data_types,
                       dyn_plants_column_select)

        
    else

        dyn_gens_data_types =
            [String, Symbol, Symbol]

        dyn_plants_data_types =
            [String, Symbol, Symbol, Bool, Symbol, Symbol]
        
    
        dyn_gens  = convert_dataframe_selected_cols_types(
            DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_gen";
                           header = true)),
            dyn_gens_column_select) ,
                     dyn_gens_data_types,
                     dyn_gens_column_select)

        dyn_plants  = convert_dataframe_selected_cols_types(
           DataFrames.select(
            DataFrame(
                XLSX.readtable(xlsx_file, "dyn_plant";
                           header = true )),
            dyn_plants_column_select) ,
                       dyn_plants_data_types,
                       dyn_plants_column_select)

        # dyn_gens  =  DataFrames.select(
        #     DataFrame(
        #         XLSX.readtable(xlsx_file, "dyn_gen";
        #                    header = true ) ),
        #     dyn_gens_column_select...) 

        # dyn_plants  =  DataFrames.select(
        #     DataFrame(
        #         XLSX.readtable(xlsx_file, "dyn_plant";
        #                    header = true,
        #                    infer_eltypes = true)),
        #     dyn_plants_column_select)
        
    end

    #--------------------------------------    
    #--------------------------------------

    vec_edge_type =
        Symbol.(mpc_branch_type_data.branch_type)
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = mpc_baseMVA)

    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = mpc_baseMVA)
    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
            vec_edge_type =
                vec_edge_type )
    
    gen_nodes_dict_static_data =
        get_gen_nodes_dict_static_data_by_mpc(
            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA =
                mpc_baseMVA )

    gen_nodes_static_tup =
        get_gen_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA=mpc_baseMVA)

    load_nodes_static_tup_data =
        get_load_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA =
                mpc_baseMVA)

    transmission_nodes_static_tup_data =
        get_transmission_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA=
                mpc_baseMVA)

    transmission_nodes_exist =
        transmission_nodes_static_tup_data == nothing ?
         false : false

    loc_loads_idx_and_locP_locQ_data =
        get_loc_loads_idx_and_locP_locQ_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = mpc_baseMVA )

    loc_loads_exist =
        length(loc_loads_idx_and_locP_locQ_data) == 0 ?
        false : true


    #--------------------------------------
    # some indices
    #--------------------------------------

   (;multi_gens_dyn_pf_fun_kwd_net_idxs,
    multi_gens_dyn_pf_fun_kwd_n2s_idxs ) =
        NamedTupleTools.select(
            get_multi_net_nodes_idxs_wt_n2s(
                mpc_gen_selected_data,
                mpc_bus_selected_data),
            (:multi_gens_dyn_pf_fun_kwd_net_idxs,
             :multi_gens_dyn_pf_fun_kwd_n2s_idxs ))
    
    #--------------------------------------    
    #--------------------------------------
    
    default_net_static_data_dict =
        OrderedDict(
            :loc_loads_exist => [loc_loads_exist],
            
            :transmission_nodes_exist =>
                [transmission_nodes_exist],
            
            :plant_generators =>
                gen_nodes_dict_static_data,
            
            :plant_loc_loads => loc_loads_exist == true ?
                OrderedDict{Union{Int64,String,Symbol},
                    NamedTuple}(
                        first(a_tup) => second(second(a_tup))
                        for a_tup in
                            loc_loads_idx_and_locP_locQ_data) : [] ,
            :plant_loads =>
                load_nodes_static_tup_data,

            :plant_transmission =>
                transmission_nodes_static_tup_data,
            
            :branches =>
                OrderedDict{
                    Union{Int64,String,Symbol},
                    NamedTuple}(a_branch.idx => a_branch
                                for a_branch in
                                    Dyn_Branches_data_list),
            :shunt =>
                [dict_shunt],

            :baseMVA =>
                [mpc_baseMVA],

            :gencost =>
                [dict_gencost] )

    return ( transmission_nodes_exist == true &&
        loc_loads_exist == true) ?
        (;transmission_nodes_exist,
         loc_loads_exist,
         multi_gens_dyn_pf_fun_kwd_net_idxs,
         multi_gens_dyn_pf_fun_kwd_n2s_idxs,
         dict_shunt,
         dict_gencost,
         Dyn_Branches_data_list,
         gen_nodes_dict_static_data,
         gen_nodes_static_tup,
         load_nodes_static_tup_data,
         transmission_nodes_static_tup_data,
         loc_loads_idx_and_locP_locQ_data,
         default_net_static_data_dict) :
             ( transmission_nodes_exist == false &&
             loc_loads_exist == true) ?
             (;transmission_nodes_exist,
              loc_loads_exist,
              multi_gens_dyn_pf_fun_kwd_net_idxs,
              multi_gens_dyn_pf_fun_kwd_n2s_idxs,
              dict_shunt,
              dict_gencost,
              Dyn_Branches_data_list,
              gen_nodes_dict_static_data,
              gen_nodes_static_tup,
              load_nodes_static_tup_data,
              loc_loads_idx_and_locP_locQ_data,
              default_net_static_data_dict) :
                  (transmission_nodes_exist == true &&
                  loc_loads_exist == false) ?
                  (;transmission_nodes_exist,
                   loc_loads_exist,
                   multi_gens_dyn_pf_fun_kwd_net_idxs,
                   multi_gens_dyn_pf_fun_kwd_n2s_idxs,
                   dict_shunt,
                   dict_gencost,
                   Dyn_Branches_data_list,
                   gen_nodes_dict_static_data,
                   gen_nodes_static_tup,
                   transmission_nodes_static_tup_data,
                   loc_loads_idx_and_locP_locQ_data,
                   default_net_static_data_dict) :
                       (;transmission_nodes_exist,
                        loc_loads_exist,
                        multi_gens_dyn_pf_fun_kwd_net_idxs,
                        multi_gens_dyn_pf_fun_kwd_n2s_idxs,
                        dict_shunt,
                        dict_gencost,
                        Dyn_Branches_data_list,
                        gen_nodes_dict_static_data,
                        gen_nodes_static_tup,
                        default_net_static_data_dict)
        
end



function get_default_static_net_json_data_by_mpc(
    case_name;        
    data_dir            = "",
    components_libs_dir = "",
    by_components       = true,
    wt_plants_data_types_bool = false)


    #--------------------------------------

    # if data_dir == ""
        
    #     data_dir = joinpath(@__DIR__,"..","..","src",
    #                  "data-dir","converted_data" )
        
    # end

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

    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath( case_data_dir, "dyn")

    json_data_dir = joinpath( case_data_dir, "json")

    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    if  wt_plants_data_types_bool == true

        dyn_gens_data_types =
            [Int, Symbol, Symbol]

        dyn_plants_data_types =
            [Int, Symbol, Symbol, Bool, Symbol, Symbol]
        
    else

        dyn_gens_data_types =
            [String, Symbol, Symbol]

        dyn_plants_data_types =
            [String, Symbol, Symbol, Bool, Symbol, Symbol]
        
    end

    #--------------------------------------
    #--------------------------------------
    
    (;mpc_branch_selected_data,
     mpc_gen_selected_data,
     mpc_gencost_data,
     mpc_bus_selected_data,
     dyn_gens,
     dyn_plants,
     mpc_baseMVA,

     mpc_branch_type_data) =
         NamedTupleTools.select(
             get_case_data_by_csv(
                 case_name ;
                 case_data_dir =
                     case_data_dir,
                 mpc_branch_column_select =
                     mpc_branch_column_select,
                 mpc_gen_column_select =
                     mpc_gen_column_select,
                 mpc_bus_column_select =
                     mpc_bus_column_select,
                 mpc_scalar_column_select =
                     mpc_scalar_column_select,
                 dyn_gens_column_select =
                     dyn_gens_column_select,
                 dyn_plants_column_select =
                     dyn_plants_column_select,
                 dyn_gens_data_types =
                     dyn_gens_data_types,
                 dyn_plants_data_types =
                     dyn_plants_data_types,
                 wt_plants_data_types_bool =
                     wt_plants_data_types_bool ),
             (:mpc_branch_selected_data,
              :mpc_gen_selected_data,
              :mpc_gencost_data,
              :mpc_bus_selected_data,
              :dyn_gens,
              :dyn_plants,
              :mpc_baseMVA,

              :mpc_branch_type_data))
            
    #--------------------------------------
    #--------------------------------------

    # vec_edge_type = [ ]

    vec_edge_type =
        Symbol.( mpc_branch_type_data.branch_type )
    
    #--------------------------------------
    
    dict_shunt =
        get_shunt_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = mpc_baseMVA)


    dict_gencost =
        get_gencost_data_by_mpc(
            mpc_gencost_data;
            mpc_baseMVA = mpc_baseMVA) 
    
    #--------------------------------------

    Dyn_Branches_data_list =
        get_Dyn_Branches_data_by_components_by_mpc(
            mpc_branch_selected_data;
            vec_edge_type =
                vec_edge_type )
    
    gen_nodes_dict_static_data =
        get_gen_nodes_dict_static_data_by_mpc(
            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA =
                mpc_baseMVA )

    gen_nodes_static_tup =
        get_gen_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data,
            mpc_gen_selected_data;
            mpc_baseMVA=mpc_baseMVA)

    load_nodes_static_tup_data =
        get_load_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA =
                mpc_baseMVA)

    transmission_nodes_static_tup_data =
        get_transmission_nodes_static_tup_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA=
                mpc_baseMVA)

    transmission_nodes_exist =
        transmission_nodes_static_tup_data == nothing ?
         false : false

    loc_loads_idx_and_locP_locQ_data =
        get_loc_loads_idx_and_locP_locQ_data_by_mpc(
            mpc_bus_selected_data;
            mpc_baseMVA = mpc_baseMVA )

    loc_loads_exist =
        length(loc_loads_idx_and_locP_locQ_data) == 0 ?
        false : true

    #--------------------------------------
    # some indices
    #--------------------------------------

   (;multi_gens_dyn_pf_fun_kwd_net_idxs,
    multi_gens_dyn_pf_fun_kwd_n2s_idxs ) =
        NamedTupleTools.select(
            get_multi_net_nodes_idxs_wt_n2s(
                mpc_gen_selected_data,
                mpc_bus_selected_data),
            (:multi_gens_dyn_pf_fun_kwd_net_idxs,
             :multi_gens_dyn_pf_fun_kwd_n2s_idxs ))

    #--------------------------------------
    #--------------------------------------

    default_net_static_data_dict =
        OrderedDict(
            :loc_loads_exist => [loc_loads_exist],
            
            :transmission_nodes_exist =>
                [transmission_nodes_exist],
            
            :plant_generators =>
                gen_nodes_dict_static_data,
            
            :plant_loc_loads => loc_loads_exist == true ?
                OrderedDict{Union{Int64,String,Symbol},
                    NamedTuple}(
                        first(a_tup) => second(second(a_tup))
                        for a_tup in
                            loc_loads_idx_and_locP_locQ_data) : [] ,
            
            :plant_loads =>
                load_nodes_static_tup_data,

            :plant_transmission =>
                transmission_nodes_static_tup_data,
            
            :branches =>
                OrderedDict{
                    Union{Int64,String,Symbol},
                    NamedTuple}(a_branch.idx => a_branch
                                for a_branch in
                                    Dyn_Branches_data_list),
            :shunt =>
                [dict_shunt],

            :baseMVA =>
                [mpc_baseMVA],

            :gencost =>
                [dict_gencost] )

    return (transmission_nodes_exist == true &&
        loc_loads_exist == true) ?
        (;transmission_nodes_exist,
         loc_loads_exist,
         multi_gens_dyn_pf_fun_kwd_net_idxs,
         multi_gens_dyn_pf_fun_kwd_n2s_idxs,
         dict_shunt,
         dict_gencost,
         Dyn_Branches_data_list,
         gen_nodes_dict_static_data,
         gen_nodes_static_tup,
         load_nodes_static_tup_data,
         transmission_nodes_static_tup_data,
         loc_loads_idx_and_locP_locQ_data,
         default_net_static_data_dict) :
             ( transmission_nodes_exist == false &&
             loc_loads_exist == true) ?
             (;transmission_nodes_exist,
              loc_loads_exist,
              multi_gens_dyn_pf_fun_kwd_net_idxs,
              multi_gens_dyn_pf_fun_kwd_n2s_idxs,
              dict_shunt,
              dict_gencost,
              Dyn_Branches_data_list,
              gen_nodes_dict_static_data,
              gen_nodes_static_tup,
              load_nodes_static_tup_data,
              loc_loads_idx_and_locP_locQ_data,
              default_net_static_data_dict) :
                  (transmission_nodes_exist == true &&
                  loc_loads_exist == false) ?
                  (;transmission_nodes_exist,
                   loc_loads_exist,
                   multi_gens_dyn_pf_fun_kwd_net_idxs,
                   multi_gens_dyn_pf_fun_kwd_n2s_idxs,
                   dict_shunt,
                   dict_gencost,
                   Dyn_Branches_data_list,
                   gen_nodes_dict_static_data,
                   gen_nodes_static_tup,
                   transmission_nodes_static_tup_data,
                   loc_loads_idx_and_locP_locQ_data,
                   default_net_static_data_dict) :
                       (;transmission_nodes_exist,
                        loc_loads_exist,
                        multi_gens_dyn_pf_fun_kwd_net_idxs,
                        multi_gens_dyn_pf_fun_kwd_n2s_idxs,
                        dict_shunt,
                        dict_gencost,
                        Dyn_Branches_data_list,
                        gen_nodes_dict_static_data,
                        gen_nodes_static_tup,
                        default_net_static_data_dict)
    
end


#---------------------------------------------------
#---------------------------------------------------
# Parameters selections  utility functions
#---------------------------------------------------
#---------------------------------------------------


function get_gens_ra_and_reactances_by_json(
    plant_generators_data_from_json ;
    sequence_order =
        (:components_data, :gen) ,
    selections =
        (:ra, :xℓ, :X_d, :X_q,
         :X_d_dash, :X_q_dash,
         :X_d_2dash, :X_q_2dash) )
    
    gens_ra = Float64[]
    gens_xℓ = Float64[]
    gens_X_d = Float64[]
    gens_X_q = Float64[]
    gens_X_d_dash = Float64[]
    gens_X_q_dash = Float64[]
    gens_X_d_2dash = Float64[]
    gens_X_q_2dash = Float64[]
    
    gens_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_generators_data_from_json, (:idx,) )) 

    #  get_components_properties_by_json
    # get_gen_plants_gens_properties_by_json
    
    gens_ra_and_reactances_properties =
        get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json ;
            sequence_order =
                sequence_order,
            selections =
                selections)

    for a_gen_properties in
        gens_ra_and_reactances_properties

        push!(gens_ra, a_gen_properties.ra )
        push!(gens_xℓ, a_gen_properties.xℓ)
        push!(gens_X_d, a_gen_properties.X_d)
        push!(gens_X_q, a_gen_properties.X_q)
        push!(gens_X_d_dash, a_gen_properties.X_d_dash)
        push!(gens_X_q_dash, a_gen_properties.X_q_dash)
        push!(gens_X_d_2dash, a_gen_properties.X_d_2dash)
        push!(gens_X_q_2dash, a_gen_properties.X_q_2dash)
        
    end

    return (;gens_ra,
            gens_xℓ,
            gens_X_d,
            gens_X_q,
            
            gens_X_d_dash,
            gens_X_q_dash,
            
            gens_X_d_2dash,
            gens_X_q_2dash)

end



function get_gens_ode_para_by_json(
    plant_generators_data_from_json;
    gens_sequence_order =
        (:components_data, :gen),
    gens_selections =
        (:D,
         :H,
         :X_d,
         :X_q,
         :X_d_dash,
         :X_q_dash,
         :T_d_dash,
         :T_q_dash ) )
    
    gens_D = Float64[]
    
    gens_H = Float64[]
    
    gens_X_d = Float64[]
    gens_X_q = Float64[]
    
    gens_X_d_dash = Float64[]
    gens_X_q_dash = Float64[]
    
    gens_T_d_dash = Float64[]
    gens_T_q_dash = Float64[]
    
    gens_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_generators_data_from_json, (:idx,))) 

    gens_ode_para =
        get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json ;
            sequence_order =
                gens_sequence_order,
            selections =
                gens_selections )

    
    for a_gen_properties in gens_ode_para

        push!(gens_D, a_gen_properties.D )
        push!(gens_H, a_gen_properties.H )
        
        push!(gens_X_d, a_gen_properties.X_d)
        push!(gens_X_q, a_gen_properties.X_q)
        
        push!(gens_X_d_dash, a_gen_properties.X_d_dash)
        push!(gens_X_q_dash, a_gen_properties.X_q_dash)
        
        push!(gens_T_d_dash, a_gen_properties.T_d_dash)
        push!(gens_T_q_dash, a_gen_properties.T_q_dash)
        
    end

    
    return (;gens_D,
            gens_H,
            
            gens_X_d,
            gens_X_q,
            
            gens_X_d_dash,
            gens_X_q_dash,
            
            gens_T_d_dash,
            gens_T_q_dash )

end


#---------------------------------------------------

function get_avrs_ode_para_by_json_generic(
    plant_generators_data_from_json;
    avrs_sequence_order =
        (:components_data, ), 
    avrs_selections =
        (:avr, ) )

    # get_components_properties_by_json
    
    avrs_ode_para =
        get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json ;
            sequence_order =
                avrs_sequence_order,
            selections =
                avrs_selections )

    
    return [ a_avr.avr for a_avr in avrs_ode_para]

end


function get_avrs_ode_para_by_json(
    plant_generators_data_from_json;
    avrs_sequence_order =
        (:components_data, :avr), 
    avrs_selections =
        (:Ka, :Ta) )
    
    avrs_Ka = Float64[]
    
    avrs_Ta = Float64[]
    

    avrs_ode_para =
        get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json ;
            sequence_order =
                avrs_sequence_order,
            selections =
                avrs_selections )


    for an_avr_properties in avrs_ode_para 

        push!(avrs_Ka, an_avr_properties.Ka )
        
        push!(avrs_Ta, an_avr_properties.Ta )
        
    end

    
    return (; avrs_Ka, avrs_Ta )

end


"""

`get_selected_avrs_ode_para_by_json` is more generic compared to

`get_avrs_ode_para_by_json_generic`,
`get_avrs_ode_para_by_json`

"""
function get_selected_avrs_ode_para_by_json(
    plant_generators_data_from_json;
    avrs_sequence_order =
        (:components_data, :avr), 
    avrs_selections =
        (:Ka, :Ta) )

    dim_avrs_selections =
        length(avrs_selections)
    
    vec_vec_para_values = [
        [] for a_para in
            1:dim_avrs_selections ]
    
    avrs_ode_para =
        get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json ;
            sequence_order =
                avrs_sequence_order,
            selections =
                avrs_selections )


    for an_avr_properties in avrs_ode_para

        for (idx, a_property) in enumerate(avrs_selections)

            push!(vec_vec_para_values[idx],
                  getproperty(
                      an_avr_properties,
                      a_property) )
        end
        
    end
    
    return vec_vec_para_values

end

#---------------------------------------------------

function get_govs_ode_para_by_json(
    plant_generators_data_from_json;
    govs_sequence_order =
        (:components_data, :gov), 
    govs_selections =
        (:R, :Ts, :Tc) )

    govs_R  = Float64[]
    
    govs_Ts = Float64[]
    
    govs_Tc = Float64[]

    govs_ode_para =
        get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json ;
            sequence_order =
                govs_sequence_order,
            selections =
                govs_selections )

   
    for an_gov_properties in govs_ode_para 

        push!(govs_R, an_gov_properties.R )
        
        push!(govs_Ts, an_gov_properties.Ts )

        push!(govs_Tc, an_gov_properties.Tc )
        
    end

    
    return (; govs_R, govs_Ts, govs_Tc )

end

#---------------------------------------------------

"""
`get_selected_comps_ode_para_by_json` can be used for any
component for a plant.

It is more generic compared to:

`get_gens_ra_and_reactances_by_json`,

`get_gens_ode_para_by_json`,

`get_avrs_ode_para_by_json_generic`,

`get_avrs_ode_para_by_json`,

`get_selected_avrs_ode_para_by_json`, and

`get_govs_ode_para_by_json`

"""
function get_selected_comps_ode_para_by_json(
    plant_generators_data_from_json;
    sequence_order =
        ( :components_data, :gen ) ,
    selections =
        (:D,
         :H,
         :X_d,
         :X_q,
         :X_d_dash,
         :X_q_dash,
         :T_d_dash,
         :T_q_dash ) )

    dim_selections =
        length( selections )
    
    comps_ode_para =
        get_components_properties_by_json(
            plant_generators_data_from_json;
            sequence_order =
                sequence_order,
            selections =
                selections )
    

    if intersect((:gen, :avr, :gov, :pss),
                 selections) == Symbol[]
        
        vec_vec_para_values = Vector{Float64}[
            [] for a_para in
                1:dim_selections ]

        for a_comp_properties in comps_ode_para

            for (idx, a_property) in enumerate( selections )

                    
                # push!(vec_vec_para_values[idx], 
                #       getproperty(
                #           a_comp_properties,
                #           a_property) )

                
                if a_property ∈ (:Sn, :vh, :P, :Q,
                                 :Pmin, :Pmax,
                                 :Qmin, :Qmax,
                                 :vmin, :vmax)

                    a_property_value_or_nothing =
                        getproperty(
                              a_comp_properties,
                            a_property)

                    a_property_value =
                    a_property_value_or_nothing == nothing ?
                        99999 : a_property_value_or_nothing

                    push!(vec_vec_para_values[idx],
                          a_property_value )                
                    
                else

                    # a proppery value can be nothing,
                    # for a gov
                    # of synchronosu condenser
                    
                    push!(vec_vec_para_values[idx],
                          getproperty(
                              a_comp_properties,
                              a_property) )
                    
                    
                end

            end

        end

        return vec_vec_para_values
        
    else
                
        vec_vec_para_values = Vector{NamedTuple}[
            [] for a_para in
                1:dim_selections ]
                
        # vec_vec_para_values = [
        #     [] for a_para in
        #         1:dim_selections ]
        
        for a_comp_properties in comps_ode_para
            
            for (idx,a_property) in enumerate( selections)

                a_property_value = getproperty(
                    a_comp_properties,
                    a_property)
                
                if a_property_value != NamedTuple{(
                    :nothing,)}
                    
                    push!(vec_vec_para_values[idx],
                          a_property_value )
                else

                    push!(vec_vec_para_values[idx],
                          (no_gov = :nothing,) )
                                        
                end
                
            end

        end

        return vec_vec_para_values
                
    end
    
end

#---------------------------------------------------
# Gen data selection function by json
#---------------------------------------------------


"""
This function returns a vector of local loads of
generators otherwise, it returns an empty vector.
"""
function get_gen_plants_local_load_by_json(
    plant_generators_data_from_json ;
    loc_load_exist=false,
    sequence_order =
        (:components_data, ),
    selections = (:loc_load,) )

    local_load =  NamedTuple[]
    
    if loc_load_exist == false
        
        return local_load
        
    else
    
        nt_local_load =
            namedtuple_nested_selection(
                plant_generators_data_from_json;
                sequence_order = sequence_order,
                selections = selections )

        for a_local_load in nt_local_load

            push!(local_load,
                  NamedTupleTools.select(
                      a_local_load.loc_load,
                      (:loc_P, :loc_Q )
                  ))
        end
        
        return local_load

    end
    
end


"""
This function returns a vector of properties of
generator nodes, selected by the variables
`sequence_order` and `selection`.

"""
function get_gen_plants_gens_properties_by_json(
    plant_generators_data_from_json ;
    sequence_order =
        (:components_data, :gen),
    selections =
        (:P, :Q))

    return namedtuple_nested_selection(
        plant_generators_data_from_json;
        sequence_order =
            sequence_order,
        selections =
            selections )

end


function get_gens_local_load_and_idx_exist_bool_by_json(
    plant_generators_data_from_json ;
    selections  =
        (:components_data, ) )
    
    local_load = NamedTuple[]
    
    gens_with_loc_load_idx = Int64[]

    nt_local_load = []
    
    for a_gen_data in plant_generators_data_from_json

        nt_gen_components_data =
            NamedTupleTools.select(
                a_gen_data,
                # plant_generators_data_from_json[1],
                selections)
        
        set_property_names =
            Set(collect(propertynames(
                nt_gen_components_data.components_data)))

        issubset(Set([:loc_load]), set_property_names ) ?
            push!(nt_local_load,
                  nt_gen_components_data.components_data.loc_load) :
                      nothing
    end

    if length(nt_local_load) == 0
        
        local_load_exist = false
        
        return (;local_load_exist,
                local_load,
                gens_with_loc_load_idx)
    else

        local_load_exist = true

        for a_nt_local_load in nt_local_load
            
            idx =
                parse(Int, split(
                    lowercase(a_nt_local_load.Bus),
                    "bus")[2])
            
            push!(gens_with_loc_load_idx, idx)

            push!(local_load,
                  (loc_P = a_nt_local_load.loc_P,
                   loc_Q = a_nt_local_load.loc_Q ))
            
        end

        
        return (;local_load_exist,
                local_load,
                gens_with_loc_load_idx)        
        
    end

end


function get_gens_plants_PQ_and_loc_loads_and_idx_by_json(
    plant_generators_data_from_json )

    sequence_order =
        (:components_data, :gen)
    
    selections =
        (:P, :Q)
    
    gens_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_generators_data_from_json,
                (:idx, ))) 

    gen_nodes_PQ =
        convert(Vector{NamedTuple},
                get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json;
            sequence_order =
                sequence_order,
            selections =
                selections))

    gen_is_slack =
        get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json ;
            sequence_order =
                (:additional_data,),
            selections =
                (:isa_slack,) )

    slack_gens_nodes_idx =
        [ idx
          for (idx, a_is_slack) in
              zip(gens_nodes_idx,
                  gen_is_slack )
              if a_is_slack.isa_slack == true ]

    non_slack_gens_nodes_idx =
        sort(collect(
            setdiff(Set(gens_nodes_idx),
                    Set(slack_gens_nodes_idx))))

    (;local_load_exist,
     local_load,
     gens_with_loc_load_idx) =
         get_gens_local_load_and_idx_exist_bool_by_json(
             plant_generators_data_from_json ;
             selections  =
                 (:components_data, ) )

    slack_bus_idx =
        slack_gens_nodes_idx

    gens_idx =
        gens_nodes_idx
    
    gens_nodes_with_loc_loads_idx =
        gens_with_loc_load_idx
    
    return (;slack_bus_idx,
            gens_idx,
            slack_gens_nodes_idx,
            non_slack_gens_nodes_idx,
            gens_nodes_idx,
            gens_with_loc_load_idx,
            gens_nodes_with_loc_loads_idx,
            local_load_exist,
            
            gen_nodes_PQ,
            local_load)

end


function get_gens_plants_and_loc_loads_idx_by_json(
    plant_generators_data_from_json )
    
    gens_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_generators_data_from_json,
                (:idx,) )) 

    gen_is_slack =
        get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json ;
            sequence_order =
                (:additional_data,),
            selections =
                (:isa_slack,) )

    slack_gens_nodes_idx =
        [ idx
          for (idx, a_is_slack) in
              zip(gens_nodes_idx,
                  gen_is_slack )
              if a_is_slack.isa_slack == true ]

    non_slack_gens_nodes_idx =
        sort(collect(
            setdiff(Set(gens_nodes_idx),
                    Set(slack_gens_nodes_idx))))

    (;local_load_exist,
     local_load,
     gens_with_loc_load_idx) =
         get_gens_local_load_and_idx_exist_bool_by_json(
             plant_generators_data_from_json ;
             selections  =
                 (:components_data, ) )

    slack_bus_idx =
        slack_gens_nodes_idx

    gens_idx =
        gens_nodes_idx
    
    gens_nodes_with_loc_loads_idx =
        gens_with_loc_load_idx
    
    return (;slack_bus_idx,
            gens_idx,
            slack_gens_nodes_idx,
            non_slack_gens_nodes_idx,
            gens_nodes_idx,
            gens_with_loc_load_idx,
            gens_nodes_with_loc_loads_idx,
            local_load_exist)

end


function get_gens_plants_PQ_and_loc_loads_by_json(
    plant_generators_data_from_json )

    sequence_order =
        (:components_data, :gen)
    
    selections =
        (:P, :Q)
    
    gen_nodes_PQ =
        convert(Vector{NamedTuple},
                get_gen_plants_gens_properties_by_json(
            plant_generators_data_from_json ;
            sequence_order =
                sequence_order,
            selections =
                selections))

    (loc_load_exist,
     local_load,
     gens_with_loc_load_idx) =
         get_gens_local_load_and_idx_exist_bool_by_json(
             plant_generators_data_from_json ;
             selections  =
                 (:components_data, ) )
    
    return (;
            loc_load_exist,            
            gen_nodes_PQ,
            local_load)

end


#---------------------------------------------------
# Non-gen data selection function by json
#---------------------------------------------------

"""
This function returns a vector of properties of non
generator nodes, selected by the variable `selection`.
"""
function get_non_gen_plants_PQ_and_idx_by_json(
        plant_loads_data_from_json,
        plant_transmission_data_from_json )

    sequence_order =
        (:components_data, )
    
    selections = (:P, :Q)

    non_gen_nodes_PQ = []
    
    if length(plant_transmission_data_from_json) != 0

        load_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_loads_data_from_json, (:idx,) )) 

        transmission_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_transmission_data_from_json,
                (:idx,)))
        
        loads_data =
            namedtuple_nested_selection(
                plant_loads_data_from_json;
                sequence_order =
                    (sequence_order[1],
                     :load),
                selections =
                    selections )
        
        transmission_data =
            namedtuple_nested_selection(
                plant_transmission_data_from_json;
                sequence_order =
                    (sequence_order[1],
                     :transmission),
                selections =
                    selections )

        dict_loads =
            OrderedDict(a_idx => a_nt
                 for (a_idx, a_nt) in
                     zip(load_nodes_idx,
                         loads_data) )

        dict_transmission =
            OrderedDict(a_idx => a_nt
                 for (a_idx, a_nt) in
                     zip(transmission_nodes_idx,
                         transmission_data) )

        non_gen_nodes_idx =
            sort([load_nodes_idx;
                  transmission_nodes_idx])

        for idx in non_gen_nodes_idx

            idx ∈ load_nodes_idx ? push!(
                non_gen_nodes_PQ,
                dict_loads[idx])  : push!(
                    non_gen_nodes_PQ,
                    dict_transmission[idx]) 
        end
        
        return (;load_nodes_idx,
                transmission_nodes_idx,
                non_gen_nodes_idx,
                non_gen_nodes_PQ)
        
    else

        load_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_loads_data_from_json, (:idx,) )) 

        loads_data =
            namedtuple_nested_selection(
                plant_loads_data_from_json;
                sequence_order =
                    (sequence_order[1],
                     :load),
                selections =
                    selections )

        for (idx, a_load) in
            zip(load_nodes_idx, loads_data )

            push!(non_gen_nodes_PQ,
                  a_load)
            
        end

        transmission_nodes_idx = []
        non_gen_nodes_idx = load_nodes_idx
        
        return (;load_nodes_idx,
                transmission_nodes_idx,
                non_gen_nodes_idx,
                non_gen_nodes_PQ )
                
    end

end


function get_non_gen_plants_idx_by_json(
        plant_loads_data_from_json,
        plant_transmission_data_from_json )

    
    if length(plant_transmission_data_from_json) != 0

        load_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_loads_data_from_json,
                (:idx,) )) 

        transmission_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_transmission_data_from_json,
                (:idx,) ) )
        

        non_gen_nodes_idx =
            sort([load_nodes_idx;
                  transmission_nodes_idx])

        
        return (;load_nodes_idx,
                transmission_nodes_idx,
                non_gen_nodes_idx )
        
    else

        load_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_loads_data_from_json,
                (:idx,) )) 

        transmission_nodes_idx = []
        non_gen_nodes_idx = load_nodes_idx
        
        return (;load_nodes_idx,
                transmission_nodes_idx,
                non_gen_nodes_idx )
                
    end

end


function get_non_gen_plants_PQ_by_json(
    plant_loads_data_from_json,
    plant_transmission_data_from_json )

    sequence_order =
        (:components_data, )

    selections = (:P, :Q)
    
    non_gen_nodes_PQ = []
    
    if length(plant_transmission_data_from_json) != 0

        load_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_loads_data_from_json,
                (:idx,))) 

        transmission_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_transmission_data_from_json,
                (:idx,)))
        
        loads_data =
            namedtuple_nested_selection(
                plant_loads_data_from_json;
                sequence_order =
                    (sequence_order[1],
                     :load),
                selections =
                    selections )
        
        transmission_data =
            namedtuple_nested_selection(
                plant_transmission_data_from_json;
                sequence_order =
                    (sequence_order[1],
                     :transmission),
                selections =
                    selections )

        dict_loads =
            OrderedDict(a_idx => a_nt
                 for (a_idx, a_nt) in
                     zip(load_nodes_idx,
                         loads_data) )

        dict_transmission =
            OrderedDict(a_idx => a_nt
                 for (a_idx, a_nt) in
                     zip(transmission_nodes_idx,
                         transmission_data) )

        non_gen_nodes_idx =
            sort([load_nodes_idx;
                  transmission_nodes_idx])

        for idx in non_gen_nodes_idx

            idx ∈ load_nodes_idx ? push!(
                non_gen_nodes_PQ,
                dict_loads[idx])  : push!(
                    non_gen_nodes_PQ,
                    dict_transmission[idx]) 
        end
        
        return (; non_gen_nodes_PQ,)
        
    else

        load_nodes_idx =
            first.(NamedTupleTools.select.(
                plant_loads_data_from_json,
                (:idx,))) 

        loads_data =
            namedtuple_nested_selection(
                plant_loads_data_from_json;
                sequence_order =
                    (sequence_order[1],
                     :load),
                selections =
                    selections )

        for (idx, a_load) in
            zip(load_nodes_idx, loads_data)

            push!(non_gen_nodes_PQ,
                  a_load)
            
        end
        
        return (; non_gen_nodes_PQ,)
                
    end

end

#---------------------------------------------------
# Shunt data selection function by json
#---------------------------------------------------


function get_nodes_shunts_Gs_and_Bs_by_json(
    shunt_data_from_json)

    return (shunt_data_from_json.shunt_Gs,
            shunt_data_from_json.shunt_Bs )

end


function get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
    shunt_data_from_json)

    return (shunt_data_from_json.shunt_idx,
            shunt_data_from_json.shunt_Gs,
            shunt_data_from_json.shunt_Bs )

end


#---------------------------------------------------
# gencost data selection function by json
#---------------------------------------------------


function get_gencost_data_by_json(
    gencost_data_from_json)

    return (;startup = gencost_data_from_json.startup,
            n = gencost_data_from_json.n,
            shutdown = gencost_data_from_json.shutdown,
            c_0 = gencost_data_from_json.c_0,
            c_1 = gencost_data_from_json.c_1,
            c_n_1 = gencost_data_from_json.c_n_1,
            cost_type = gencost_data_from_json.cost_type
            )

end


function get_gens_cost_coeff_in_ascen(
    gencost_data )

    (c_0, c_1, c_n_1) =
        NamedTupleTools.select(
            gencost_data,
            (:c_0, :c_1, :c_n_1))
    
    return [ (coeff_0, coeff_1, coeff_2)
          for (coeff_0, coeff_1, coeff_2) in
              zip(c_0, c_1, c_n_1 )]
    
end


function get_gens_cost_coeff_in_decen(
    gencost_data)

    (c_0, c_1, c_n_1) =
        NamedTupleTools.select(
            gencost_data,
            (:c_0, :c_1, :c_n_1))
    
    [ ( coeff_2, coeff_1,  coeff_0)
          for (coeff_2, coeff_1,  coeff_0) in
              zip( c_n_1, c_1, c_0 )]
    
end

#---------------------------------------------------
# Edges data selection function by json
#---------------------------------------------------


function get_edges_nt_fbus_tbus_by_json(
    edge_data_from_json)


    edges_fbus =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.fbus
          for nt_edge in edge_data_from_json]

    edges_tbus =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.tbus
          for nt_edge in edge_data_from_json]

    return (;edges_fbus, edges_tbus)
    
    
end


function get_edges_fbus_tbus_by_json(
    edge_data_from_json)


    edges_fbus =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.fbus
          for nt_edge in edge_data_from_json]

    edges_tbus =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.tbus
          for nt_edge in edge_data_from_json]

    return (edges_fbus, edges_tbus)
    
    
end


function get_edges_generic_data_by_json(
    edge_data_from_json )
    
    edges_idx =
        [ (NamedTupleTools.select(
            nt_edge, (:idx,))).idx
          for nt_edge in
              edge_data_from_json]
    
    edges_r =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.r
          for nt_edge in edge_data_from_json]

    edges_x =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.x
          for nt_edge in edge_data_from_json]

    edges_b =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.b
          for nt_edge in edge_data_from_json]

    edges_ratio =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.ratio
          for nt_edge in edge_data_from_json]

    edges_angle =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.angle
          for nt_edge in edge_data_from_json]

    edges_type =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_type,))).components_type.edge_type
          for nt_edge in edge_data_from_json]
    

    return ( edges_r, edges_x, edges_b,
             edges_ratio, edges_angle,
             edges_type)

    
end


function get_edges_ftbus_and_generic_data_by_json(
    edge_data_from_json )

    (edges_fbus, edges_tbus) =
        get_edges_fbus_tbus_by_json(
            edge_data_from_json)

    ( edges_r, edges_x,
      edges_b,  edges_ratio,
      edges_angle, edges_type) =
          get_edges_generic_data_by_json(
              edge_data_from_json )

    return (edges_fbus, edges_tbus,
            edges_r, edges_x, edges_b,
            edges_ratio, edges_angle, edges_type)
    
end



function get_edges_idx_and_generic_data_by_json(
    edge_data_from_json;
    sequence_order =
        (:components_data, ),
    selections =
        (:r, :x, :b, :ratio,
         :angle))

    edges_idx =
        [ (NamedTupleTools.select(
            nt_edge, (:idx,))).idx
          for nt_edge in edge_data_from_json]

    edges_fbus =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.fbus
          for nt_edge in edge_data_from_json]

    edges_tbus =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.tbus
          for nt_edge in edge_data_from_json]

    edges_r =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.r
          for nt_edge in edge_data_from_json]


    edges_x =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.x
          for nt_edge in edge_data_from_json]


    edges_b =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.b
          for nt_edge in edge_data_from_json]


    edges_ratio =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.ratio
          for nt_edge in edge_data_from_json]


    edges_angle =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_data,))).components_data.angle
          for nt_edge in edge_data_from_json]

    edges_type =
        [ (NamedTupleTools.select(
            nt_edge,
            (:components_type,))).components_type.edge_type
          for nt_edge in edge_data_from_json]

    return (; edges_idx, edges_fbus, edges_tbus,
            edges_r, edges_x, edges_b,
            edges_ratio, edges_angle, edges_type)
    
end

#---------------------------------------------------
# Migrated from sd-dynamics-fault-events
#---------------------------------------------------


function get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
    r, x, b,
    ratio, angle,    
    Gs, Bs;
    branches_fbus = branches_fbus,
    branches_tbus = branches_tbus,
    edge_type = edge_type,
    buses_idx = buses_idx,
    baseMVA = 1.0,
    basekV = 1.0,
    baseShunt = 1.0 )

    nodes_idx_and_Yshunt =
        get_nodes_idx_and_Yshunt_by_generic(
            buses_idx,
            Gs, Bs;
            baseMVA =
                baseShunt )
    
    edges_orientation =
        get_edges_orientation_by_generic(
            branches_fbus,
            branches_tbus )

    nodes_incident_edges =
        get_nodes_incident_edges_by_orientations(
            edges_orientation )
    
    edges_Ybr_cal =
        get_edges_Ybr_by_generic(
            r, x, b,
            ratio, angle,
            edge_type,
            Gs, Bs;
            baseMVA = baseMVA,
            basekV = basekV )

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
    
    nodes_idx_with_adjacent_nodes_idx = [
        [ [[idx],
           [ idx == orient[1] ?
               orient[2] : orient[1]
             for orient in
                 a_node_edges_other_nodes_idx]]...; ]
        for (idx, a_node_edges_other_nodes_idx) in
            zip(buses_idx,
                nodes_incident_edges_and_orientation)]
    
    #------------------------------------------

    nodes_incident_edges_orientation_and_Ybr_cal = [
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
        for (node_idx, orientations_and_edges_Ybr) in 
            zip(buses_idx,
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

    return (; 
            Ynet,
            nodes_idx_with_adjacent_nodes_idx )
    
end


function get_edges_ftbus_and_generic_data_by_json(
    edge_data_from_json,
    nothing )

    (branches_fbus, branches_tbus) =
        get_edges_fbus_tbus_by_json(
            edge_data_from_json)

    ( r, x,
      b,  ratio,
      angle, edges_type) =
          get_edges_generic_data_by_json(
              edge_data_from_json )

    return (;branches_fbus, branches_tbus,
            r, x, b,
            ratio, angle, edges_type)
    
end


function get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
    shunt_data_from_json,
    nothing)

    buses_idx = shunt_data_from_json.shunt_idx
    Gs = shunt_data_from_json.shunt_Gs
    Bs = shunt_data_from_json.shunt_Bs
    
    return (;buses_idx,
            Gs,
            Bs )

end


#---------------------------------------------------
# input files and namedtuple related functions by json
#---------------------------------------------------

function get_a_gen_plant_data_json_to_nt(
    a_plant_data_json;
    in_components_type_sym = false)

    if in_components_type_sym != false
        return NamedTupleTools.merge(
            NamedTupleTools.select(
                a_plant_data_json,
                (:idx, :plant_type)),
            (components_type =
                get_nt_components_type(
                a_plant_data_json.components_type),),
            (components_data =
                # get_nt_of_dict_to_nt_of_nt(
                # a_plant_data_json.components_data),),
                get_nt_components_data(
                a_plant_data_json.components_data),),
            (additional_data =
                get_nt_additional_data(
                a_plant_data_json.additional_data),) )
        
    else
        return NamedTupleTools.merge(
            NamedTupleTools.select(
                a_plant_data_json,
                (:idx, :plant_type)),
            (components_type = namedtuple(
                a_plant_data_json.components_type),),
            (components_data =
                # get_nt_of_dict_to_nt_of_nt(
                # a_plant_data_json.components_data),),
                get_nt_components_data(
                a_plant_data_json.components_data),),
            (additional_data = get_nt_additional_data(
                a_plant_data_json.additional_data),) )

    end
    

end



function get_a_non_gen_plant_data_json_to_nt(
    a_plant_data_json;
    in_components_type_sym = false)

    if in_components_type_sym != false
        return NamedTupleTools.merge(
            NamedTupleTools.select(
                a_plant_data_json,
                ( :idx, :plant_type )),
            (components_type =
                get_nt_components_type(
                a_plant_data_json.components_type),),
            (components_data =
                # get_nt_of_dict_to_nt_of_nt(
                #     a_plant_data_json.components_data),
                get_nt_components_data(
                    a_plant_data_json.components_data),
             )
        ) 
        
    else
        return NamedTupleTools.merge(
            NamedTupleTools.select(
                a_plant_data_json,
                (:idx, :plant_type)),
            (components_type = namedtuple(
                a_plant_data_json.components_type),),
            (components_data =
                # get_nt_of_dict_to_nt_of_nt(
                #     a_plant_data_json.components_data),
                get_nt_components_data(
                    a_plant_data_json.components_data),
             )) 

    end
    

end


function get_an_edge_data_json_to_nt(
    an_edge_data_json;
    in_components_type_sym = false)

    if in_components_type_sym != false
        return NamedTupleTools.merge(
            NamedTupleTools.select(
                an_edge_data_json,
                ( :idx, )),
            (components_type =
                get_nt_components_type(
                an_edge_data_json.components_type),),
            (components_data =
                namedtuple(
                an_edge_data_json.components_data),) )
        
    else
        return NamedTupleTools.merge(
            NamedTupleTools.select(
                an_edge_data_json,
                ( :idx, )),
            (components_type = namedtuple(
                an_edge_data_json.components_type),),
            (components_data =
                namedtuple(
                an_edge_data_json.components_data),))

    end
    

end


#---------------------------------------------------
#---------------------------------------------------
# components instance related functions by json
#---------------------------------------------------
#---------------------------------------------------

function get_gov_instance_from_json(
    gov_nt_params_from_file,
    dict_gov_sym_type,
    sym_gov_nt_params)

    json_gov_param =
        getproperty(gov_nt_params_from_file,
                    sym_gov_nt_params)

    # gov_param_from_json =
    # namedtuple(
    #     convert(Dict{Symbol,Float64},
    #             json_gov_param ))

    
    gov_param_from_json =
        namedtuple(
            convert(Dict{Symbol,Float64},
                getproperty(gov_nt_params_from_file,
                    sym_gov_nt_params) ))

    sym_gov_type =
        Symbol( split(String(sym_gov_nt_params),"__")[1])

    return dict_gov_sym_type[sym_gov_type](
        ;gov_param_from_json... )
    
end


function get_avr_instance_from_json(
    avr_nt_params_from_file,
    dict_avr_sym_type,
    sym_avr_nt_params)

    json_avr_param =
        getproperty(avr_nt_params_from_file,
                    sym_avr_nt_params)

    # avr_param_from_json =
    # namedtuple(
    #     convert(Dict{Symbol,Float64},
    #             json_avr_param ))

    
    avr_param_from_json =
        namedtuple(
            convert(Dict{Symbol,Float64},
                getproperty(avr_nt_params_from_file,
                    sym_avr_nt_params) ))

    sym_avr_type =
        Symbol( split(String(sym_avr_nt_params),"__")[1])

    return dict_avr_sym_type[sym_avr_type](
        ;avr_param_from_json... )
    
end

#---------------------------------------------------
#---------------------------------------------------

function get_gov_instance(
    gov_type, gov_nt_params)

    return gov_type(;gov_nt_params...)
    
end


function get_govs_instances(
    sym_govs_nt_params,
    dict_gov_nt_params,
    dict_gov_sym_type)

    sym_govs_types =
        Symbol.(
            first.(split.(
                String.(sym_govs_nt_params),"__") ) )

    govs_types =
        [ sym_gov_type == :nothing ? :nothing :
        dict_gov_sym_type[sym_gov_type]
          for sym_gov_type in
              sym_govs_types  ]

    govs_nt_params =
        [ sym_gov_nt_params == :nothing ? :nothing :
        dict_gov_nt_params[sym_gov_nt_params]
          for sym_gov_nt_params in
              sym_govs_nt_params  ]

    return [
        sym_gov_nt_param == nothing ?
            nothing : gov_type(;gov_nt_params...)
            for (gov_type,
                 gov_nt_params,
                 sym_gov_nt_param) in
                zip(govs_types,
                    govs_nt_params,
                    sym_govs_nt_params )]
    
end


function get_avr_instance(
    avr_type, avr_nt_params)

    return avr_type(;avr_nt_params...)
    
end


function get_avrs_instances(
    sym_avrs_nt_params,
    dict_avr_nt_params,
    dict_avr_sym_type)

    sym_avrs_types =
        Symbol.(
            first.(split.(
                String.(sym_avrs_nt_params),"__") ) )

    avrs_types =
        [ dict_avr_sym_type[sym_avr_type]
          for sym_avr_type in
              sym_avrs_types  ]

    avrs_nt_params =
        [ dict_avr_nt_params[sym_avr_nt_params]
          for sym_avr_nt_params in
              sym_avrs_nt_params  ]

    return [avr_type(;avr_nt_params...)
            for (avr_type, avr_nt_params) in
                zip(avrs_types, avrs_nt_params )]
    
end


function get_a_gen_instance_by_static_and_dym_data(
    a_gen_type, a_gen_dym_data, a_gen_a_gen_static_data)

    return  a_gen_type(; a_gen_dym_data...,
                       a_gen_a_gen_static_data...)
    
end


function get_gens_instance_by_static_and_dym_data(
    gens_type_and_dym_data,
    gens_nodes_static_tup_data)

    @assert first.(gens_type_and_dym_data) ==
        first.(gens_nodes_static_tup_data)

    return [
        first(a_gen_type_and_dym_data)(;
            second(a_gen_type_and_dym_data)...,
            a_gen_static_tup_data...)

        for (a_gen_type_and_dym_data,
             a_gen_static_tup_data) in
            zip(second.(gens_type_and_dym_data),
                second.(gens_nodes_static_tup_data))]
    
end


#---------------------------------------------------
#---------------------------------------------------


function get_plants_type_and_components_type(
    dyn_plants)

    return [
        (idx,(plant_type, gen_type,
              isa_slack, gov_type, avr_type)) 
        for (idx, plant_type, gen_type, isa_slack,
             gov_type,avr_type) in
            zip(dyn_plants.bus,
                dyn_plants.Plant_type,
                dyn_plants.Gen,
                dyn_plants.isa_slack,
                dyn_plants.Gov,
                dyn_plants.Exc)]

end


function get_a_gen_plant_instances(
    idx_plant, plant_type, gen, isa_slack,
    gov, avr, loc_load;
    loc_load_exist = false,
    p_order  = 1.0,
    v_ref = 1.0,
    ω_ref = ωs )

    if loc_load_exist != false
        if gov != nothing
            return plant_type(
                Gen = gen,
                Gov = gov,
                Exc = avr,
                isa_slack = isa_slack,
                Loc_load = loc_load,
                p_order = p_order,
                v_ref   = v_ref,
                ω_ref   = ω_ref)
        else            
            return plant_type(
                Gen = gen,
                Exc = avr,
                isa_slack = isa_slack,
                Loc_load = loc_load,
                p_order = p_order,
                v_ref   = v_ref,
                ω_ref   = ω_ref)
        end
    else
        if gov != nothing
            return plant_type(
                Gen = gen,
                Gov = gov,
                Exc = avr,
                isa_slack = isa_slack,
                p_order = p_order,
                v_ref   = v_ref,
                ω_ref   = ω_ref)
        else            
            return plant_type(
                Gen = gen,
                Exc = avr,
                isa_slack = isa_slack,
                p_order = p_order,
                v_ref   = v_ref,
                ω_ref   = ω_ref)
        end
    end

end


function get_a_gen_plant_no_loc_load_instances(
    idx_plant, plant_type, gen, isa_slack,
    gov, avr;    
    p_order  = 1.0,
    v_ref = 1.0,
    ω_ref = ωs )

    if gov != nothing
        return plant_type(
            Gen = gen,
            Gov = gov,
            Exc = avr,
            isa_slack = isa_slack,
            p_order = p_order,
            v_ref   = v_ref,
            ω_ref   = ω_ref)
    else            
        return plant_type(
            Gen = gen,
            Exc = avr,
            isa_slack = isa_slack,
            p_order = p_order,
            v_ref   = v_ref,
            ω_ref   = ω_ref,
            isa_condenser = true)
    end
end


function get_a_gen_plant_wt_loc_load_instances(
    idx_plant, plant_type,
    gen, isa_slack,
    gov, avr, loc_load;
    p_order  = 1.0,
    v_ref = 1.0,
    ω_ref = ωs )

    if gov != nothing
        return plant_type(
            Gen = gen,
            Gov = gov,
            Exc = avr,
            isa_slack = isa_slack,
            Loc_load = loc_load,
            p_order = p_order,
            v_ref   = v_ref,
            ω_ref   = ω_ref)
    else            
        return plant_type(
            Gen = gen,
            Exc = avr,
            isa_slack = isa_slack,
            Loc_load = loc_load,
            p_order = p_order,
            v_ref   = v_ref,
            ω_ref   = ω_ref,
            isa_condenser = true)
    end
        
end



function get_gens_plant_instances(
    dict_plants_gen_sym_type,
    dict_gen_sym_type,
    dict_gens_dyn_nt_params,

    dict_gov_sym_type,
    dict_gov_nt_params,

    dict_avr_sym_type,
    dict_avr_nt_params,
    
    dyn_gens,
    dyn_plants,

    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0,

    p_order = 1.0,
    v_ref = 1.0,
    ω_ref = ωs )

    #------------------------------------------

    loc_load_exist =
        loc_load_exist_bool_by_mpc(
            mpc_bus )

    loc_loads =
        get_loc_loads_idx_and_locP_locQ_data_by_mpc(
            mpc_bus; mpc_baseMVA = mpc_baseMVA )

    
    gens_nodes_idx =
        get_gens_nodes_idx_by_mpc(
            mpc_bus)
    
    gens_with_loc_load_idx =
        get_gens_nodes_with_loc_loads_idx_by_mpc(
            mpc_bus)
    
    # n2s_gens_idx =
    #     get_a_n2s_net_group(gens_nodes_idx)

    # n2s_gens_with_loc_load_idxs =
    #     get_a_n2s_net_group(gens_with_loc_load_idx;
    #                         loc_load_exist = true)

    
    n2s_gens_idx =
        get_n2s_any(gens_nodes_idx)

    n2s_gens_with_loc_load_idxs = loc_load_exist = true ?
        get_n2s_any(gens_with_loc_load_idx) :
        get_n2s_any(gens_with_loc_load_idx;
                    nothing_bool= true )
    
    #------------------------------------------    
    
    plants_idx =
        dyn_plants.bus

    sym_plants_types =
        dyn_plants.Plant_type
    
    sym_gens_types =
        dyn_plants.Gen
    
    bool_isa_slack =
        dyn_plants.isa_slack
    
    sym_govs_nt_params =
        dyn_plants.Gov
    
    sym_avrs_nt_params =
        dyn_plants.Exc

    
    # gens_plants_type_and_components_type =
    #     get_plants_type_and_components_type(
    #         dyn_plants)

    
    # sym_govs_types =
    #     Symbol.(
    #         first.(split.(
    #         String.(sym_govs_nt_params),"__") ) )

    # sym_avrs_types =
    #     Symbol.(
    #         first.(split.(
    #             String.(sym_avrs_nt_params),"__") ) )

    plants_types =
        [dict_plants_gen_sym_type[sym_plant_type]
         for sym_plant_type in
             sym_plants_types]

    # nodes_static_tup_data =
    #     get_nodes_static_tup_data_by_mpc(
    #         mpc_bus, mpc_gen;
    #         mpc_baseMVA=mpc_baseMVA)
    
    # nodes_dict_static_data =
    #     get_nodes_dict_static_data_by_mpc(
    #         mpc_bus,
    #         mpc_gen;
    #         mpc_baseMVA=mpc_baseMVA)

    gens_type_and_dym_data =
        get_gens_type_and_dym_data(
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,
            dyn_gens)

    gens_nodes_static_tup_data =
        get_gen_nodes_static_tup_data_by_mpc(
            mpc_bus,
            mpc_gen;
            mpc_baseMVA=mpc_baseMVA)

    #------------------------------------------

    gens_instances =
        get_gens_instance_by_static_and_dym_data(
            gens_type_and_dym_data,
            gens_nodes_static_tup_data)

    govs_instances =
        get_govs_instances(
            sym_govs_nt_params,
            dict_gov_nt_params,
            dict_gov_sym_type)

    avrs_instances =
        get_avrs_instances(
            sym_avrs_nt_params,
            dict_avr_nt_params,
            dict_avr_sym_type)

    loc_load_instances =
        loc_load_exist != false ? [
            loc_Load_t1(;a_loc_load...) 
            for a_loc_load in loc_loads  ] : []
    
    return [ plant_idx ∈ gens_with_loc_load_idx ?
        (plant_idx,
         get_a_gen_plant_wt_loc_load_instances(
            plant_idx,
            plant_type,
            gen,
            isa_slack,
            gov,
            avr,
            loc_load_instances[
                n2s_gens_with_loc_load_idxs[plant_idx]];
            p_order  = p_order,
            v_ref = v_ref,
            ω_ref = ω_ref )) :
                (plant_idx,
                 get_a_gen_plant_no_loc_load_instances(
                    plant_idx,
                    plant_type,
                    gen,
                    isa_slack,
                    gov,
                    avr;
                    p_order = p_order,
                    v_ref = v_ref,
                    ω_ref = ω_ref) )
             for (plant_idx, plant_type,
                  gen, isa_slack, gov, avr) in
                 zip(gens_nodes_idx,
                    plants_types,
                    gens_instances,
                    bool_isa_slack,
                    govs_instances,
                    avrs_instances)  ]
    
end

#--------------------------------------

function get_load_nodes_plant_instances(
    mpc_bus;
    mpc_baseMVA=1.0,
    load_type =
        PQ_Const_I,
    load_node_type =
        plant_PQ_Const_I)

    load_nodes_static_tup_data =
        get_load_nodes_static_tup_data_by_mpc(
            mpc_bus;
            mpc_baseMVA =
                mpc_baseMVA)

        return [(idx,
          load_node_type(Load =
              load_type(
                  ;load_node_static_tup_data...))) 
         for (idx, load_node_static_tup_data) in
                second.(
                    load_nodes_static_tup_data)]
    
end

function get_transmission_nodes_plant_instances(
    mpc_bus;
    mpc_baseMVA=1.0,
    transmission_type =
        Trans_t2_Node,
    transmission_node_type =
        plant_Transmission_t2 )

    transmission_nodes_idx =
        get_transmission_nodes_idx_by_mpc(
            mpc_bus)

    if length(transmission_nodes_idx) != 0

        transmission_nodes_static_tup_data =
            get_transmission_nodes_static_tup_data_by_mpc(
                mpc_bus;
                mpc_baseMVA=mpc_baseMVA)
        
        return [
            (idx,
             transmission_node_type(Trans =
                 transmission_type(
                     ;transmission_node_static_tup_data...)))

            for (idx,transmission_node_static_tup_data) in
                second.(
                    transmission_nodes_static_tup_data)]
    else

        return []

    end
    
    
end


#---------------------------------------------------
#---------------------------------------------------

function get_Dyn_Branches_by_mpc(
    mpc_branch;
    mpc_baseMVA=1.0,
    basekV = 1.0 )

    branches_data =
        get_branches_data_and_types_by_mpc(
            mpc_branch,
            mpc_baseMVA=mpc_baseMVA, basekV = basekV  )

    Branches =
        [ branch_type == :line ?
        (idx, PiModelLine(;  branch_data... ) ) :
        (idx, Transformer(; branch_data... ))
          for (idx, branch_data, branch_type) in
              branches_data ]
    
    return OrderedDict("branch$(no)" => branch
                       for (no, branch) in
                           Branches)
end




function get_Dyn_Nodes_by_mpc(
    mpc_bus;
    mpc_baseMVA=1.0)

    (;slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     nodes_with_demands_idx) =
         get_net_nodes_type_idxs_by_mpc(
             mpc_bus)

    loads = [(idx, (Bus = "bus$(idx)",
                    P = Pd/mpc_baseMVA,
                    Q = Qd/mpc_baseMVA))
             for (idx, Pd, Qd) in
                 zip(all_nodes_idx,
                     mpc_bus.Pd,
                     mpc_bus.Qd)
                 if idx ∈ load_nodes_idx ]

    transmission = [(idx, (Bus = "bus$(idx)",
                           P = Pd/mpc_baseMVA,
                           Q = Qd/mpc_baseMVA) )
             for (idx, Pd, Qd) in
                 zip(all_nodes_idx,
                     mpc_bus.Pd,
                     mpc_bus.Qd)
                    if idx ∈ transmission_nodes_idx ]

    # get_non_gen_node_tup_idx_PQ_data_by_mpc(
    #     idx, Pd, Qd;mpc_baseMVA = 1.0 )
    # get_loc_load_tup_idx_PQ_data_by_mpc(
    #     idx, Pd, Qd;mpc_baseMVA = 1.0 )
    
    loc_loads = [(idx, (Bus = "bus$(idx)",
                    loc_P = Pd/mpc_baseMVA,
                    loc_Q = Qd/mpc_baseMVA))
             for (idx, Pd, Qd) in
                 zip(all_nodes_idx,
                     mpc_bus.Pd,
                     mpc_bus.Qd)
                 if idx ∈ gens_nodes_with_loc_loads_idx ]
        
    gens_static_data =
        [(idx,
          (Bus  = "bus$(idx)",
           vmax = mpc_bus.Vmax[idx],
           vmin = mpc_bus.Vmin[idx],
           vh   = Vg,
           P    = Pg/mpc_baseMVA,
           Q    = Qg/mpc_baseMVA,
           Qmax = Qmax/mpc_baseMVA,
           Qmin = Qmin/mpc_baseMVA,
           Pmax = Pmax/mpc_baseMVA,
           Pmin = Pmin/mpc_baseMVA
           ) )
         for (idx,Vg,Pg,Qg,Qmax,Qmin,Pmax,Pmin) in
             zip(mpc_gen.bus,  mpc_gen.Vg,
                 mpc_gen.Pg,   mpc_gen.Qg,
                 mpc_gen.Qmax, mpc_gen.Qmin,
                 mpc_gen.Pmax, mpc_gen.Pmin)]
    
end



function get_Dyn_Nodes_by_mpc(
    dict_plants_gen_sym_type,
    dict_gen_sym_type,
    dict_gens_dyn_nt_params,

    dict_gov_sym_type,
    dict_gov_nt_params,

    dict_avr_sym_type,
    dict_avr_nt_params,
    
    dyn_gens,
    dyn_plants,

    mpc_bus,
    mpc_gen;
    mpc_baseMVA=1.0,

    p_order = 1.0,
    v_ref = 1.0,
    ω_ref = ωs,

    load_type =
        PQ_Const_I,
    load_node_type =
        plant_PQ_Const_I,

    transmission_type =
        Trans_t2_Node,
    transmission_node_type =
        plant_Transmission_t2)
    
    gens_nodes_idx =
        get_gens_nodes_idx_by_mpc(
            mpc_bus)
    
    load_nodes_idx =
        get_load_nodes_idx_by_mpc(
            mpc_bus)
    
    transmission_nodes_idx =
        get_transmission_nodes_idx_by_mpc(
            mpc_bus)
    
    all_nodes_idx =
        get_all_nodes_idx_by_mpc(
            mpc_bus)
    
    # n2s_gens_idx =
    #     get_a_n2s_net_group(
    #         gens_nodes_idx)

    # n2s_load_nodes_idx =
    #     get_a_n2s_net_group(
    #         load_nodes_idx)
    
    # n2s_transmission_idxs =
    #     length(transmission_nodes_idx) != 0 ?
    #     get_a_n2s_net_group(
    #         transmission_nodes_idx;
    #         transmission_group = true) : []

    
    n2s_gens_idx =
        get_n2s_any(
            gens_nodes_idx)

    n2s_load_nodes_idx =
        get_n2s_any(
            load_nodes_idx)
    
    n2s_transmission_idxs =
        length(transmission_nodes_idx) != 0 ?
        get_n2s_any(
            transmission_nodes_idx) : get_n2s_any(
                transmission_nodes_idx;
                nothing_bool= true)
    
    
    gens_plant_instances =
        get_gens_plant_instances(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus,
            mpc_gen;
            mpc_baseMVA = mpc_baseMVA,

            p_order = p_order,
            v_ref = v_ref,
            ω_ref = ω_ref )
    
    load_nodes_plant_instances =
        get_load_nodes_plant_instances(
            mpc_bus;
            mpc_baseMVA=mpc_baseMVA,
            load_type = load_type,
            load_node_type = load_node_type)

    transmission_nodes_plant_instances =
        get_transmission_nodes_plant_instances(
            mpc_bus;
            mpc_baseMVA=mpc_baseMVA,
            transmission_type = transmission_type,
            transmission_node_type =
                transmission_node_type )

    if length(transmission_nodes_idx) != 0        
        tup_idx_and_plants =
            [idx ∈ gens_nodes_idx ?
            gens_plant_instances[n2s_gens_idx[idx]] :
            idx ∈ load_nodes_idx ?
            load_nodes_plant_instances[
                n2s_load_nodes_idx[idx]] :
                    transmission_nodes_plant_instances[
                        n2s_transmission_idxs[idx]]
              for idx in all_nodes_idx ]
        
    else
        tup_idx_and_plants =
            [idx ∈ gens_nodes_idx ?
            gens_plant_instances[n2s_gens_idx[idx]] :
            load_nodes_plant_instances[
                n2s_load_nodes_idx[idx]]
             for idx in all_nodes_idx]
    end
        
    return OrderedDict( "bus$(idx)" => plant
                        for (idx, plant) in
                            tup_idx_and_plants )
    
    
end


function get_Dyn_Nodes_Dyn_Branches_by_mpc(
    ;case_name = "case9",        
    data_dir = "",
    components_libs_dir = "" )
    
    #--------------------------------------
    
    if (components_libs_dir == "") || (
        components_libs_dir == nothing )

        package_dir = pkgdir(ePowerSim)

        src_dir =
            joinpath( package_dir, "src" )
        
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
                 case_name )
    
    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath(case_data_dir, "dyn")

    json_data_dir = joinpath(case_data_dir, "json")

    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    dyn_gens_data_types =
        [Int, Symbol, Symbol]

    dyn_plants_data_types =
        [Int, Symbol, Symbol, Bool, Symbol, Symbol]

    #--------------------------------------
    
    (;mpc_branch_selected_data,
     mpc_gen_selected_data,
     mpc_gencost_data,
     mpc_bus_selected_data,
     dyn_gens,
     dyn_plants,
     mpc_baseMVA) =
         NamedTupleTools.select(
             get_case_data_by_csv(
                 case_namex;
                 case_data_dir =
                     case_data_dir,
                 mpc_branch_column_select =
                     mpc_branch_column_select,
                 mpc_gen_column_select =
                     mpc_gen_column_select,
                 mpc_bus_column_select =
                     mpc_bus_column_select,
                 mpc_scalar_column_select =
                     mpc_scalar_column_select,
                 dyn_gens_column_select =
                     dyn_gens_column_select,
                 dyn_plants_column_select =
                     dyn_plants_column_select,
                 dyn_gens_data_types =
                     dyn_gens_data_types,
                 dyn_plants_data_types =
                     dyn_plants_data_types ),
             (:mpc_branch_selected_data,
              :mpc_gen_selected_data,
              :mpc_gencost_data,
              :mpc_bus_selected_data,
              :dyn_gens,
              :dyn_plants,
              :mpc_baseMVA))
    
    #--------------------------------------

    sub_components_strings =
        ["avrs", "govs",  "pss",
         "gens", "loads", "lines"]

    (avrs_libs_dir,
     govs_libs_dir,
     pss_libs_dir,
     gens_libs_dir,
     loads_libs_dir,
     lines_libs_dir) =
         get_sub_components_libs_dir(
             components_libs_dir,
             sub_components_strings )
    
    #--------------------------------------

    components_files_string =
        ["avrs-libs",  "govs-libs",
         "pss-libs",   "gens-libs",
         "loads-libs", "lines-libs" ]

    (avrs_type_libs_file_json,
     govs_type_libs_file_json,
     pss_type_libs_file_json,
     gens_type_libs_file_json,
     loads_type_libs_file_json,
     lines_type_libs_file_json) =
        get_sub_components_libs_files(
            components_libs_dir,
            components_files_string;
            ext = "json" )

    #--------------------------------------

    parameters_libs_files_strings =
        ["avrs-parameters-libs",
         "govs-parameters-libs",
         "pss-parameters-libs",
         "gens-dyn-parameters-libs"]

    (avrs_parameters_libs_file_json,
     govs_parameters_libs_file_json,
     pss_parameters_libs_file_json,
     gens_dyn_parameters_libs_file_json ) =
        get_sub_components_libs_files(
            components_libs_dir,
            parameters_libs_files_strings;
            ext = "json" )

    #--------------------------------------

    gov_json_data_file =
        joinpath(json_data_dir,
                 "dict_gov_nt_params.json") 

    avr_json_data_file =
        joinpath(json_data_dir,
                 "dict_avr_nt_params.json") 

    gens_nt_dynamic_params_json_data_file =
        joinpath(json_data_dir,
                 "gens_nt_dynamic_params.json") 

    sym_gens_dynamic_params_json_data_file =
        joinpath(json_data_dir,
                 "sym_gens_dynamic_params.json") 

    dict_gens_dyn_nt_params_json_data_file =
        joinpath(json_data_dir,
                 "dict_gens_dyn_nt_params.json")
    
    #--------------------------------------

    dict_plants_gen_sym_type =
        get_abstract_type_dict_subtypes(
            SdGenPlant )

    dict_gen_sym_type =
        get_abstract_type_dict_subtypes(
            SdGen )

    dict_gov_sym_type =
        get_abstract_type_dict_subtypes(
            SdGov )

    dict_avr_sym_type =
        get_abstract_type_dict_subtypes(
            SdAvr )

    dict_pss_sym_type =
        get_abstract_type_dict_subtypes(
            SdPss )


    #--------------------------------------
    
    dict_gens_dyn_nt_params =
        get_dict_nt_params_from_json_lib_file(
            gens_dyn_parameters_libs_file_json )    
    
    dict_gov_nt_params =
        get_dict_nt_params_from_json_lib_file(
            govs_parameters_libs_file_json)
    
    dict_avr_nt_params =
        get_dict_nt_params_from_json_lib_file(
            avrs_parameters_libs_file_json)
    
    dict_pss_nt_params =
        get_dict_nt_params_from_json_lib_file(
            pss_parameters_libs_file_json)

    #--------------------------------------
    
    Dyn_Branches =
        get_Dyn_Branches_by_mpc(
            mpc_branch_selected_data;
            mpc_baseMVA=1.0 )

    Dyn_Nodes =
        get_Dyn_Nodes_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus_selected_data,
            mpc_gen_selected_data ;
            mpc_baseMVA=1.0,

            p_order = 1.0,
            v_ref = 1.0,
            ω_ref = ωs,

            load_type =
                PQ_Const_I,
            load_node_type =
                plant_PQ_Const_I,

            transmission_type =
                Trans_t2_Node,
            transmission_node_type =
                plant_Transmission_t2)

    return (; Dyn_Nodes, Dyn_Branches )
end


#---------------------------------------------------


function get_Dyn_Nodes_Dyn_Branches_data_by_mpc(
    ;case_name = "case9",        
    data_dir = "",
    components_libs_dir = "",
    by_components = false )

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
    
    #--------------------------------------
    
    # components_libs_dir =
    #     joinpath(
    #         lib_dir,
    #         components_lib )
    
    #--------------------------------------

    mpc_data_dir  = case_data_dir

    dyn_data_dir  = joinpath(case_data_dir, "dyn")

    json_data_dir = joinpath(case_data_dir, "json")

    #--------------------------------------
    
    mpc_branch_column_select =
        ["fbus", "tbus", "r", "x", "b",
         "ratio", "angle", "status"]

    mpc_gen_column_select =
        ["bus", "Pg", "Qg", "Qmax", "Qmin",
         "Vg", "mBase", "status", "Pmax","Pmin"]

    mpc_bus_column_select =
        ["bus_i", "type", "Pd",
         "Qd", "Gs", "Bs", "Vmax", "Vmin"]

    mpc_scalar_column_select =
        ["mpc_baseMVA" ]

    dyn_gens_column_select =
        ["bus","sym_gen_type",
         "sym_gen_dynamic_para"]

    dyn_plants_column_select =
        ["bus","Plant_type","Gen",
         "isa_slack","Gov","Exc"]

    #--------------------------------------

    dyn_gens_data_types =
        [Int, Symbol, Symbol]

    dyn_plants_data_types =
        [Int, Symbol, Symbol, Bool, Symbol, Symbol]

    #--------------------------------------
    #--------------------------------------
    
    (;mpc_branch_selected_data,
     mpc_gen_selected_data,
     mpc_gencost_data,
     mpc_bus_selected_data,
     dyn_gens,
     dyn_plants,
     mpc_baseMVA) =
         NamedTupleTools.select(
             get_case_data_by_csv(
                 case_name ;
                 case_data_dir =
                     case_data_dir,
                 mpc_branch_column_select =
                     mpc_branch_column_select,
                 mpc_gen_column_select =
                     mpc_gen_column_select,
                 mpc_bus_column_select =
                     mpc_bus_column_select,
                 mpc_scalar_column_select =
                     mpc_scalar_column_select,
                 dyn_gens_column_select =
                     dyn_gens_column_select,
                 dyn_plants_column_select =
                     dyn_plants_column_select,
                 dyn_gens_data_types =
                     dyn_gens_data_types,
                 dyn_plants_data_types =
                     dyn_plants_data_types ),
             (:mpc_branch_selected_data,
              :mpc_gen_selected_data,
              :mpc_gencost_data,
              :mpc_bus_selected_data,
              :dyn_gens,
              :dyn_plants,
              :mpc_baseMVA))
    
    #--------------------------------------    
    #--------------------------------------

    (;dict_plants_gen_sym_type,

     dict_gen_sym_type,
     dict_gens_dyn_nt_params,

     dict_gov_sym_type,
     dict_gov_nt_params,

     dict_avr_sym_type,
     dict_avr_nt_params,

     dict_pss_sym_type,
     dict_pss_nt_params) =
         NamedTupleTools.select(
           get_dynamic_components_parameters_libs_from_json(
                 ;components_libs_dir =
                     components_libs_dir),
             (:dict_plants_gen_sym_type,

              :dict_gen_sym_type,
              :dict_gens_dyn_nt_params,

              :dict_gov_sym_type,
              :dict_gov_nt_params,

              :dict_avr_sym_type, 
              :dict_avr_nt_params,

              :dict_pss_sym_type,
              :dict_pss_nt_params ) )

    #--------------------------------------
    #--------------------------------------
    
    Dyn_Branches_data_dict =
        get_Dyn_Branches_data_by_mpc(
            mpc_branch_selected_data; mpc_baseMVA=1.0 )

    Dyn_Nodes_data_dict =
        get_Dyn_Nodes_data_by_mpc(
            dict_plants_gen_sym_type,
            dict_gen_sym_type,
            dict_gens_dyn_nt_params,

            dict_gov_sym_type,
            dict_gov_nt_params,

            dict_avr_sym_type,
            dict_avr_nt_params,

            dyn_gens,
            dyn_plants,

            mpc_bus_selected_data,
            mpc_gen_selected_data ;
            mpc_baseMVA=1.0,

            p_order = 1.0,
            v_ref = 1.0,
            ω_ref = ωs,

            load_type =
                PQ_Const_I,
            load_node_type =
                plant_PQ_Const_I,

            transmission_type =
                Trans_t2_Node,
            transmission_node_type =
                plant_Transmission_t2,
            by_components = by_components)

    # return (;Dyn_Nodes_data_dict,
    #         Dyn_Branches_data_dict)
    
    return  OrderedDict{Symbol, OrderedDict{
        String, NamedTuple}}(
            :Nodes => Dyn_Nodes_data_dict,
            :Branches => Dyn_Branches_data_dict )
    
end

#---------------------------------------------------
# Case data selection function by json
#---------------------------------------------------

function get_net_parameters_and_Idx_by_components_by_json(
    case_name ;
    data_dir = "",
    net_data_by_components_file = "",    
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true )

    #--------------------------------------

    # case_name = "case14"

    # data_dir =
    #     joinpath(@__DIR__,"..","..","src","data-dir",
    #              "converted_data" )

    # net_data_by_components_file =
    #     joinpath(json_case_dir,
    #              "net_data_by_components_file.json")

    #--------------------------------------
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end
    
    #--------------------------------------

    case_data_dir =
        joinpath(data_dir,
                 "converted-data",
                 case_name)
    
    json_case_dir =
        joinpath( case_data_dir, "json")

    #--------------------------------------

    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     "net_data_by_components_file.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end
    
    #--------------------------------------
    # read    
    #--------------------------------------
 
    (;plant_generators_data_from_json,
     plant_loads_data_from_json,
     plant_transmission_data_from_json,
     edge_data_from_json,
     shunt_data_from_json,
     baseMVA_data_from_json,
     gencost_data_from_json) =
     get_net_data_by_components_from_json_file(
             json_net_data_by_components_file;
             in_components_type_sym = false )

    #--------------------------------------

    baseMVA = baseMVA_data_from_json
    
    #--------------------------------------

    (;slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_with_loc_load_idx,
     gens_nodes_with_loc_loads_idx,
     local_load_exist,
     gen_nodes_PQ,
     local_load) =
         get_gens_plants_PQ_and_loc_loads_and_idx_by_json(
             plant_generators_data_from_json )

    (;load_nodes_idx,
     transmission_nodes_idx,
     non_gen_nodes_idx,
     non_gen_nodes_PQ ) =
         get_non_gen_plants_PQ_and_idx_by_json(             
                 plant_loads_data_from_json,             
                 plant_transmission_data_from_json )

    non_slack_gens_and_non_gens_idx =
        sort([non_slack_gens_nodes_idx;
              non_gen_nodes_idx])
    
    nodes_with_demands_idx =
        convert(Vector{Int64},
                sort([load_nodes_idx;
                      gens_with_loc_load_idx]))
    
    all_nodes_idx =
        sort([gens_nodes_idx;
              non_gen_nodes_idx])

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
         all_nodes_idx,
         non_slack_gens_and_non_gens_idx,
         nodes_with_demands_idx)

    n2s_idx =
        get_dict_net_streamlined_idx_by_nodes_type_idxs(
             net_nodes_type_idxs)

    (;n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         n2s_idx
    
    #------------------------------------------------
    
   gens_ra_and_reactances =
         get_gens_ra_and_reactances_by_json(
             plant_generators_data_from_json ;
             sequence_order =
                 (:components_data, :gen),
             selections =
                 (:ra, :xℓ, :X_d, :X_q,
                  :X_d_dash, :X_q_dash,
                  :X_d_2dash, :X_q_2dash) )


    (;gens_ra,
     gens_xℓ,
     gens_X_d,
     gens_X_q,
     
     gens_X_d_dash,
     gens_X_q_dash,
     
     gens_X_d_2dash,
     gens_X_q_2dash) =
         gens_ra_and_reactances
    
    #------------------------------------------------

    (ra,
     xℓ,
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     X_d_2dash,
     X_q_2dash) =
         get_selected_comps_ode_para_by_json(
             plant_generators_data_from_json;
             sequence_order =
                 (:components_data, :gen),
             selections =
                 (:ra, :xℓ,
                  :X_d, :X_q,
                  
                  :X_d_dash, :X_q_dash,
                  :X_d_2dash, :X_q_2dash) )
    
    #------------------------------------------------

    (H,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash,
     T_d_dash,
     T_q_dash ) = get_selected_comps_ode_para_by_json(
             plant_generators_data_from_json;
             sequence_order =
                 (:components_data, :gen),
             selections =
                 (:H,
                  :X_d, :X_q,                  
                  :X_d_dash, :X_q_dash,
                  :T_d_dash, :T_q_dash ) )
    
    #------------------------------------------------

    ( H,
      Xd_dash ) =
        get_selected_comps_ode_para_by_json(
             plant_generators_data_from_json;
             sequence_order =
                 (:components_data, :gen),
             selections =
                 (:H,
                  :X_d_dash ) )
    
    ode_spcm_gens_para =
        (; H,
          Xd_dash )
    
    #------------------------------------------------

    (H,
     X_d,
     X_q,
     X_d_dash, 
     T_d_dash
      ) = get_selected_comps_ode_para_by_json(
             plant_generators_data_from_json;
             sequence_order =
                 (:components_data, :gen),
             selections =
                 (:H,
                  :X_d,
                  :X_q,                  
                  :X_d_dash
                  :T_d_dash ) )

    ode_flux_decay_gens_para =
        (; H, X_d, X_q, X_d_dash, T_d_dash )
    
    #------------------------------------------------

    (Ka,
     Ta ) =
          get_selected_comps_ode_para_by_json(
              plant_generators_data_from_json;
              sequence_order =
                  ( :components_data, :avr ),
              selections =
                  (:Ka, :Ta ))

    ode_flux_decay_avrs_para =
        (;Ka,
         Ta)
     
    #------------------------------------------------
    
    # get_gens_ode_flux_decay_para_by_json(
    #     plant_generators_data_from_json )
    
    #------------------------------------------------
    
    # edges_orientation =
    #     get_edges_orientation_by_generic(
    #         get_edges_fbus_tbus_by_json(
    #             edge_data_from_json)... )

    # edges_Ybr =
    #     get_edges_Ybr_by_generic(
    #         get_edges_generic_data_by_json(
    #             edge_data_from_json)...,
    #         get_nodes_shunts_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA = baseMVA )

    # nodes_Yshunt =
    #     get_nodes_Yshunt_by_generic(
    #         get_nodes_shunts_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA =
    #             baseMVA )

    # nodes_idx_and_Yshunt =
    #     get_nodes_idx_and_Yshunt_by_generic(
    #         get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA =
    #             baseMVA )
        
    # Ybr_cal_and_edges_orientation =
    #     get_edges_Ybr_cal_and_edges_orientation_by_generic(
    #         get_edges_fbus_tbus_by_json(
    #             edge_data_from_json)...,
    #         get_edges_generic_data_by_json(
    #             edge_data_from_json )...,
    #         get_nodes_shunts_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA = baseMVA)

    # # (; edges_Ybr_cal,
    # #  edges_orientation ) =
    # #      Ybr_cal_and_edges_orientation

    # Ynet_wt_nodes_idx_wt_adjacent_nodes =
    #     get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
    #         get_edges_ftbus_and_generic_data_by_json(
    #             edge_data_from_json )...,
    #         get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA = baseMVA )

    # # (; 
    # #  Ynet,
    # #  nodes_idx_with_adjacent_nodes_idx ) =
    # #      Ynet_wt_nodes_idx_wt_adjacent_nodes

    (;edges_orientation,
     edges_Ybr,
     nodes_Yshunt,
     nodes_idx_and_Yshunt,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json,
                 edge_data_from_json,
                 shunt_data_from_json;
                 baseMVA = baseMVA,
                 basekV  = basekV,
                 use_pu_in_PQ = use_pu_in_PQ,
                 line_data_in_pu = line_data_in_pu ),
             (:edges_orientation,
              :edges_Ybr,
              :nodes_Yshunt,
              :nodes_idx_and_Yshunt,
              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes) )
    

    # sta_pf_PQ_param =
    #     get_sta_pf_PQ_param_by_json(
    #         plant_generators_data_from_json,
    #         plant_loads_data_from_json,
    #         plant_transmission_data_from_json;
    #         baseMVA = baseMVA )

    sta_pf_PQ_param =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA = baseMVA,
            use_pu_in_PQ = use_pu_in_PQ )
    
    # (; P_gens,
    #  Q_gens,
    #  P_non_gens,
    #  Q_non_gens,
    #  P_g_loc_load,
    #  Q_g_loc_load,
    #  loc_load_exist) =
    #      sta_pf_PQ_param

    gens_vh_slack_θh_para =
        get_gens_vh_slack_θh_para_by_json(
            plant_generators_data_from_json )

    # (; slack_gens_vh,
    #  slack_gens_θh,
    #  gens_vh,
    #  non_slack_gens_vh ) =
    #      gens_vh_slack_θh_para

    sta_pf_vars_and_paras_idx =
        get_sta_pf_vars_and_paras_idx_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    # (;
    #  red_types_Idxs_etc,
    #  PQ_sta_para_Idxs,
    #  nodes_types_idxs,
    #  n2s_idxs ) =
    #      sta_pf_vars_and_paras_idx

    pf_sta_ΔPQ_mismatch_parameters =
        get_pf_sta_ΔPQ_mismatch_parameters_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json;
            baseMVA = baseMVA,
            basekV = basekV,
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu )

    # (;
    #  pf_kw_para,
    #  pf_PQ_param,
    #  red_types_Idxs_etc,
    #  net_para  ) =
    #      pf_sta_ΔPQ_mismatch_parameters

    return (;plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json,
            baseMVA_data_from_json,

            gencost_data_from_json,
            
            baseMVA,
            net_nodes_type_idxs,
            n2s_idx,
            edges_orientation,
            edges_Ybr,
            nodes_Yshunt,
            nodes_idx_and_Yshunt,
            Ybr_cal_and_edges_orientation,
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            sta_pf_PQ_param,
            gens_vh_slack_θh_para,
            sta_pf_vars_and_paras_idx,
            pf_sta_ΔPQ_mismatch_parameters)
    
end


function get_net_powerflow_data_by_components_by_json(
    case_name ;
    data_dir = "",
    net_data_by_components_file = "",    
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true )

    #--------------------------------------
    
    if (data_dir == "") || (data_dir == nothing)

        package_dir = pkgdir(ePowerSim)
        
        data_dir =
            joinpath( package_dir, "data")


    end
    
    #--------------------------------------

    case_data_dir =
        joinpath(data_dir,
                 "converted-data",
                 case_name)
    
    json_case_dir =
        joinpath( case_data_dir, "json")

    #--------------------------------------

    if (net_data_by_components_file == "" ||
        net_data_by_components_file == nothing)
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     "net_data_by_components_file.json")
    else
            
        json_net_data_by_components_file =
            joinpath(json_case_dir,
                     net_data_by_components_file)
        
    end
    
    #--------------------------------------
    # read    
    #--------------------------------------
 
    (;plant_generators_data_from_json,
     plant_loads_data_from_json,
     plant_transmission_data_from_json,
     edge_data_from_json,
     shunt_data_from_json,
     baseMVA_data_from_json,
     gencost_data_from_json) =
     get_net_data_by_components_from_json_file(
             json_net_data_by_components_file;
             in_components_type_sym = false )

    #--------------------------------------

    baseMVA = baseMVA_data_from_json
    
    #--------------------------------------

    (;slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_with_loc_load_idx,
     gens_nodes_with_loc_loads_idx,
     local_load_exist,
     gen_nodes_PQ,
     local_load) =
         get_gens_plants_PQ_and_loc_loads_and_idx_by_json(
             plant_generators_data_from_json )

    (;load_nodes_idx,
     transmission_nodes_idx,
     non_gen_nodes_idx,
     non_gen_nodes_PQ ) =
         get_non_gen_plants_PQ_and_idx_by_json(             
                 plant_loads_data_from_json,             
                 plant_transmission_data_from_json )

    non_slack_gens_and_non_gens_idx =
        sort([non_slack_gens_nodes_idx;
              non_gen_nodes_idx])
    
    nodes_with_demands_idx =
        convert(Vector{Int64},
                sort([load_nodes_idx;
                      gens_with_loc_load_idx]))
    
    all_nodes_idx =
        sort([gens_nodes_idx;
              non_gen_nodes_idx])
    
    (;gens_ra,
     gens_xℓ,
     gens_X_d,
     gens_X_q,
     gens_X_d_dash,
     gens_X_q_dash,
     gens_X_d_2dash,
     gens_X_q_2dash) =
         get_gens_ra_and_reactances_by_json(
             plant_generators_data_from_json ;
             sequence_order =
                 (:components_data, :gen),
             selections =
                 (:ra, :xℓ, :X_d, :X_q,
                  :X_d_dash, :X_q_dash,
                  :X_d_2dash, :X_q_2dash) )

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
         all_nodes_idx,
         non_slack_gens_and_non_gens_idx,
         nodes_with_demands_idx)

    n2s_idx =
        get_dict_net_streamlined_idx_by_nodes_type_idxs(
             net_nodes_type_idxs)

    (;n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         n2s_idx
    
    #------------------------------------------------
    
    # edges_orientation =
    #     get_edges_orientation_by_generic(
    #         get_edges_fbus_tbus_by_json(
    #             edge_data_from_json)... )

    # edges_Ybr =
    #     get_edges_Ybr_by_generic(
    #         get_edges_generic_data_by_json(
    #             edge_data_from_json)...,
    #         get_nodes_shunts_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA = baseMVA )

    # nodes_Yshunt =
    #     get_nodes_Yshunt_by_generic(
    #         get_nodes_shunts_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA =
    #             baseMVA )

    # nodes_idx_and_Yshunt =
    #     get_nodes_idx_and_Yshunt_by_generic(
    #         get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA =
    #             baseMVA )
    

    # Ybr_cal_and_edges_orientation =
    #     get_edges_Ybr_cal_and_edges_orientation_by_generic(
    #         get_edges_fbus_tbus_by_json(
    #             edge_data_from_json)...,
    #         get_edges_generic_data_by_json(
    #             edge_data_from_json )...,
    #         get_nodes_shunts_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA = baseMVA)

    # # (; edges_Ybr_cal,
    # #  edges_orientation ) =
    # #      Ybr_cal_and_edges_orientation

    # Ynet_wt_nodes_idx_wt_adjacent_nodes =
    #     get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
    #         get_edges_ftbus_and_generic_data_by_json(
    #             edge_data_from_json )...,
    #         get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA = baseMVA )

    # # (; 
    # #  Ynet,
    # #  nodes_idx_with_adjacent_nodes_idx ) =
    # #      Ynet_wt_nodes_idx_wt_adjacent_nodes

    # -------------------------------------

    (;edges_orientation,
     edges_Ybr,
     nodes_Yshunt,
     nodes_idx_and_Yshunt,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json,
                 edge_data_from_json,
                 shunt_data_from_json;
                 baseMVA = baseMVA,
                 basekV  = basekV,
                 use_pu_in_PQ = use_pu_in_PQ,
                 line_data_in_pu = line_data_in_pu ),
             (:edges_orientation,
              :edges_Ybr,
              :nodes_Yshunt,
              :nodes_idx_and_Yshunt,
              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes) )

    # -------------------------------------
    
    # sta_pf_PQ_param =
    #     get_sta_pf_PQ_param_by_json(
    #         plant_generators_data_from_json,
    #         plant_loads_data_from_json,
    #         plant_transmission_data_from_json;
    #         baseMVA = baseMVA,
    #         use_pu_in_PQ = use_pu_in_PQ )

    sta_pf_PQ_param =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA = baseMVA,
            use_pu_in_PQ = use_pu_in_PQ )
    
    # (; P_gens,
    #  Q_gens,
    #  P_non_gens,
    #  Q_non_gens,
    #  P_g_loc_load,
    #  Q_g_loc_load,
    #  loc_load_exist) =
    #      sta_pf_PQ_param

    gens_vh_slack_θh_para =
        get_gens_vh_slack_θh_para_by_json(
            plant_generators_data_from_json )

    # (; slack_gens_vh,
    #  slack_gens_θh,
    #  gens_vh,
    #  non_slack_gens_vh ) =
    #      gens_vh_slack_θh_para

    sta_pf_vars_and_paras_idx =
        get_sta_pf_vars_and_paras_idx_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    # (;
    #  red_types_Idxs_etc,
    #  PQ_sta_para_Idxs,
    #  nodes_types_idxs,
    #  n2s_idxs ) =
    #      sta_pf_vars_and_paras_idx

    pf_sta_ΔPQ_mismatch_parameters =
        get_pf_sta_ΔPQ_mismatch_parameters_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json;
            baseMVA = baseMVA,
            basekV = basekV,
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu )

    # (;
    #  pf_kw_para,
    #  pf_PQ_param,
    #  red_types_Idxs_etc,
    #  net_para  ) =
    #      pf_sta_ΔPQ_mismatch_parameters

    return (;plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json,
            baseMVA_data_from_json,

            gencost_data_from_json,
            
            baseMVA,
            net_nodes_type_idxs,
            n2s_idx,
            edges_orientation,
            edges_Ybr,
            nodes_Yshunt,
            nodes_idx_and_Yshunt,
            Ybr_cal_and_edges_orientation,
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            sta_pf_PQ_param,
            gens_vh_slack_θh_para,
            sta_pf_vars_and_paras_idx,
            pf_sta_ΔPQ_mismatch_parameters)
    
end


#---------------------------------------------------
#---------------------------------------------------


function get_sta_pf_PQ_param_by_mpc(
    mpc_gen,
    mpc_bus;
    baseMVA = 1.0  )

    loc_load_exist =
        loc_load_exist_bool_by_mpc( mpc_bus )

    P_gens = mpc_gen.Pg ./ baseMVA
    
    Q_gens = mpc_gen.Qg ./ baseMVA
    
    P_non_gens =
        [ node_Pd for (node_type, node_Pd) in
             zip( mpc_bus.type, mpc_bus.Pd)
             if  node_type == 1  ] ./ baseMVA
    
    Q_non_gens =
        [ node_Qd for (node_type, node_Qd ) in
             zip( mpc_bus.type, mpc_bus.Qd )
             if  node_type == 1  ] ./ baseMVA

    
     P_Q_g_loc_load =
         [ ( node_Pd, node_Qd )
          for ( node_type, node_Pd, node_Qd ) in
             zip( mpc_bus.type, mpc_bus.Pd, mpc_bus.Qd )
             if (node_type == 3 || node_type == 2) &&
                 ((node_Pd != 0.0 || node_Pd != 0) || (
                     node_Qd != 0.0 || node_Qd != 0) ) ]
    
    P_g_loc_load = length( P_Q_g_loc_load ) == 0 ?
        [] : first.( P_Q_g_loc_load )

    Q_g_loc_load = length( P_Q_g_loc_load ) == 0 ?
        [] : second.( P_Q_g_loc_load )
    
     # P_g_loc_load =
     #    [node_Pd for ( node_type, node_Pd ) in
     #         zip( mpc_bus.type, mpc_bus.Pd )
     #         if (node_type == 3 || node_type == 2) &&
     #             (node_Pd != 0.0 || node_Pd != 0) ]

     # Q_g_loc_load =
     #    [node_Qd for ( node_type, node_Qd ) in
     #         zip( mpc_bus.type, mpc_bus.Qd )
     #         if (node_type == 3 || node_type == 2) &&
     #             (node_Qd != 0.0 || node_Qd != 0) ]

    P_g_loc_load = length( P_g_loc_load ) == 0 ?
        [] : P_g_loc_load ./ baseMVA

    Q_g_loc_load = length( Q_g_loc_load ) == 0 ?
        [] : Q_g_loc_load ./ baseMVA

    return (;P_gens,
            Q_gens,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,
            loc_load_exist)
    
end


function get_pf_PQ_param_by_mpc(
    mpc_gen,
    mpc_bus;
    baseMVA = 1.0,
    use_pu_in_PQ = true )

    if use_pu_in_PQ == true
        
        (;P_gens,
         Q_gens,
         P_non_gens,
         Q_non_gens,
         P_g_loc_load,
         Q_g_loc_load,
         loc_load_exist) =
             get_sta_pf_PQ_param_by_mpc(
                 mpc_gen,
                 mpc_bus;
                 baseMVA = baseMVA  )
    else
        
        (;P_gens,
         Q_gens,
         P_non_gens,
         Q_non_gens,
         P_g_loc_load,
         Q_g_loc_load,
         loc_load_exist) =
             get_sta_pf_PQ_param_by_mpc(
                 mpc_gen,
                 mpc_bus;
                 baseMVA = 1.0  )        
    end
    
    return (;P_gens,
            Q_gens,
            P_non_gens,
            Q_non_gens,
            P_g_loc_load,
            Q_g_loc_load,
            loc_load_exist )
    
end


function get_gens_vh_slack_θh_para_by_mpc(
    mpc_gen,
    mpc_bus )

    gens_vh = mpc_gen.Vg
    
    # ----------------------------------------------
    
    n2s_streamlined_idx =
        get_dict_n2s_streamlined_idx_by_mpc( mpc_bus )

    n2s_gens_idx =
        n2s_streamlined_idx.n2s_gens_idx
    
    n2s_slack_gens_idx =
        n2s_streamlined_idx.n2s_slack_gens_idx

    n2s_non_slack_gens_idx =
        n2s_streamlined_idx.n2s_non_slack_gens_idx

    # ----------------------------------------------
    
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_mpc( mpc_bus  )

    gens_idx =
        net_nodes_type_idxs.gens_idx
    
    slack_gens_nodes_idx =
        net_nodes_type_idxs.slack_gens_nodes_idx
    
    non_slack_gens_nodes_idx =
        net_nodes_type_idxs.non_slack_gens_nodes_idx

    # -----------------------------------------------

    slack_gens_vh =
        [ gens_vh[n2s_gens_idx[a_slack_idx]]
          for a_slack_idx in
              slack_gens_nodes_idx ]
    
    slack_gens_θh = zeros(length(slack_gens_vh))
    
    non_slack_gens_vh =
        [ gens_vh[n2s_gens_idx[a_non_slack_idx]]
          for a_non_slack_idx in
              non_slack_gens_nodes_idx  ]
    
    # -----------------------------------------------

    return (;slack_gens_vh,
            slack_gens_θh,
            gens_vh,
            non_slack_gens_vh )
end



function get_sta_pf_vars_and_paras_idx_by_mpc(
    mpc_bus)

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
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     nodes_with_demands_idx)  =
         get_net_nodes_type_idxs_by_mpc(
             mpc_bus  )
    
    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx
    
    #------------------------------------------        

    nodes_size =
        length(gens_nodes_idx) +
        length(non_gens_nodes_idx)
    
    #------------------------------------------

    nodes_types_idxs =
        (;slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,         
         gens_with_loc_load_idx,
         all_nodes_idx)
    
    #------------------------------------------        

    (;
     n2s_slack_gens_idx,         
     n2s_non_slack_gens_idx,     
     n2s_gens_idx,               
     n2s_non_gens_idx,           
     n2s_load_idx,               
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,      
     n2s_all_nodes_idx ) =
         get_dict_n2s_streamlined_idx_by_mpc(
             mpc_bus)


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
    
    ur_ui_IDX    = create_idxs(
        ur_ui_offset, ur_ui_dims )
    
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
    
    # ------------------------------------------------

    transformed_gens_idx = [
        n2s_all_nodes_idx[idx] for idx in gens_idx ]
    
    transformed_non_slack_gens_nodes_idx = [
        n2s_all_nodes_idx[idx] for idx in
            non_slack_gens_nodes_idx ]
    
    transformed_non_gens_nodes_idx = [
        n2s_all_nodes_idx[idx] for idx in
            non_gens_nodes_idx ]
    
    transformed_slack_bus_idx = [
        n2s_all_nodes_idx[idx] for idx in
            slack_bus_idx ]

    transformed_non_slack_gens_and_non_gens_idx = [
        n2s_all_nodes_idx[idx] for idx in
            non_slack_gens_and_non_gens_idx ]
    
    non_gens_vh_idx = setdiff(
        vh_IDX,
        transformed_gens_idx)

    non_slack_gens_θh_idx =
        θh_IDX[ transformed_non_slack_gens_nodes_idx ]

    non_gens_θh_idx =
        θh_IDX[ transformed_non_gens_nodes_idx ]

    red_θh_idx =
        setdiff(θh_IDX,
                θh_IDX[ transformed_slack_bus_idx ])
    
    # -----------------------------------------------
    
    red_vh_θh_idx =
        [  non_gens_vh_idx...;
           non_slack_gens_θh_idx...;
           non_gens_θh_idx... ]

    # -----------------------------------------------

    red_var_comps_idxs =
        (; non_gens_vh_idx,
         non_slack_gens_θh_idx,
         non_gens_θh_idx )

    # ------------------------------------------------

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
        first(red_non_slack_gens_θh_Idxs):last(
            red_non_gens_θh_Idxs)
    
     # ----------------------------------------------- 
    
     red_dict_θh_idx2Idx = OrderedDict{Int64, Int64}(
         idx => Idx for (idx, Idx) in
            zip(transformed_non_slack_gens_and_non_gens_idx,
                      red_θh_Idxs ) )

    red_dict_θh_idx2Idx_in_Idx =
        OrderedDict{Int64, Int64}(
           idx => Idx for (idx, Idx) in
            zip( transformed_non_slack_gens_and_non_gens_idx,
                 1:length(red_θh_Idxs) ) )

     red_non_slack_gens_θh_idx2Idx =
         [ red_dict_θh_idx2Idx[idx]
          for idx in
              transformed_non_slack_gens_nodes_idx ]

     red_non_slack_gens_θh_idx2Idx_in_Idx =
         [ red_dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              transformed_non_slack_gens_nodes_idx ]

     red_non_gens_θh_idx2Idx =
         [ red_dict_θh_idx2Idx[idx]
          for idx in
              transformed_non_gens_nodes_idx ]

     red_non_gens_θh_idx2Idx_in_Idx =
         [ red_dict_θh_idx2Idx_in_Idx[idx]
           for idx in
              transformed_non_gens_nodes_idx ]
    
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

         red_vh_θh_IDX )    
    
    # ------------------------------------------------

    dim_P_gens = dim_Q_gens =
        length(gens_nodes_idx)
    
    dim_P_non_gens = dim_Q_non_gens =
        length(non_gens_nodes_idx)

    dim_P_g_loc_load = dim_Q_g_loc_load =
        loc_load_exist == true ?
        length(gens_nodes_with_loc_loads_idx) : nothing
   
    # ----------------------------------------------------

    if loc_load_exist == true

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

        PQ_sta_para_Idxs =
            (; P_gens_sta_para_Idxs,
             Q_gens_sta_para_Idxs,
             P_non_gens_sta_para_Idxs,
             Q_non_gens_sta_para_Idxs,
             P_g_loc_load_sta_para_Idxs,
             Q_g_loc_load_sta_para_Idxs )

    else

        dim_pf_sta_PQ_para  =
                [ dim_P_gens,
                  dim_Q_gens,
                  dim_P_non_gens,
                  dim_Q_non_gens ]

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

        P_g_loc_load_sta_para_Idxs = nothing

        Q_g_loc_load_sta_para_Idxs = nothing

        PQ_sta_para_Idxs =
            (; P_gens_sta_para_Idxs,
             Q_gens_sta_para_Idxs,
             P_non_gens_sta_para_Idxs,
             Q_non_gens_sta_para_Idxs,
             P_g_loc_load_sta_para_Idxs,
             Q_g_loc_load_sta_para_Idxs )


    end
    
    # --------------------------------------------------

    return (;red_types_Idxs_etc,
            PQ_sta_para_Idxs,
            nodes_types_idxs,
            n2s_idxs )    
end


function get_pf_sta_ΔPQ_mismatch_parameters_by_mpc(
    mpc_gen,
    mpc_bus,
    mpc_branch,
    mpc_baseMVA;
    basekV = 1.0,
    use_pu_in_PQ = true,
    line_data_in_pu = true)
    
    #-----------------------------------------------
    
    loc_load_exist =
        loc_load_exist_bool_by_mpc(
            mpc_bus )
    
    #-----------------------------------------------
    
    (;
     red_types_Idxs_etc,
     PQ_sta_para_Idxs,
     nodes_types_idxs,
     n2s_idxs ) =
         get_sta_pf_vars_and_paras_idx_by_mpc( mpc_bus)

    #-----------------------------------------------

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

    
    (;P_gens_sta_para_Idxs,
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

    # -----------------------------------------------

    (;P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist
     ) =
         get_pf_PQ_param_by_mpc(
             mpc_gen,
             mpc_bus;
             baseMVA = mpc_baseMVA,
             use_pu_in_PQ = use_pu_in_PQ )
            
    # ---------------------------------------------------

    (; slack_gens_vh,
     slack_gens_θh,
     gens_vh,
     non_slack_gens_vh ) =
         get_gens_vh_slack_θh_para_by_mpc(
             mpc_gen,
             mpc_bus )

    # ---------------------------------------------------

    (;Ynet,
     nodes_idx_with_adjacent_nodes_idx,
     edges_Ybr_cal,
     edges_orientation) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_mpc(
                 mpc_bus,
                 mpc_branch,
                 mpc_baseMVA;
                 basekV = basekV,
                 use_pu_in_PQ = use_pu_in_PQ,
                 line_data_in_pu = line_data_in_pu ),
             (:Ynet,
              :nodes_idx_with_adjacent_nodes_idx,
              :edges_Ybr_cal,
              :edges_orientation))
    
    net_para = (;
                Ynet,
                nodes_idx_with_adjacent_nodes_idx,
                edges_Ybr_cal,
                edges_orientation)
    
    # -------------------------------------------------
    
    pf_kw_gens_vh_slack_θh_para =
        (;slack_gens_vh,
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
        (;n2s_slack_gens_idx,
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

    return (;
            pf_kw_para,
            pf_PQ_param,
            red_types_Idxs_etc,
            net_para  )
    
end


######################################################

function get_pf_intg_var_Idx_by_mpc(
    mpc_bus )

    # ----------------------------------------------------

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
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     nodes_with_demands_idx)  =
         get_net_nodes_type_idxs_by_mpc(
             mpc_bus  )
    
    gens_with_loc_load_idx =
        gens_nodes_with_loc_loads_idx
    
    # ----------------------------------------------------

    nodes_size =
        length(gens_nodes_idx) +
        length(non_gens_nodes_idx)

    # ----------------------------------------------------

    ur_ui_dims   =
        [ nodes_size, nodes_size ]
    
    ur_ui_offset =
        create_offsets( ur_ui_dims )
    
    ur_ui_IDX    =
        create_idxs(ur_ui_offset, ur_ui_dims )
    
    ur_IDX       = ur_ui_IDX[1]
    ui_IDX       = ur_ui_IDX[2]

    ir_IDX       = ur_IDX
    ii_IDX       = ui_IDX
    
    vh_IDX       = ur_IDX
    θh_IDX       = ui_IDX

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
         create_size_offset_Idx(
             dim_intg_vh_θh_id_iq  )

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

    (intg_slack_gens_vh_idx2Idx,
     intg_non_slack_gens_vh_idx2Idx,
     intg_non_gens_vh_idx2Idx) =
         vec_types_intg_vh_idx2Idx
        

    vec_types_intg_θh_idx2Idx = [
        [ intg_dict_θh_idx2Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    (intg_slack_gens_θh_idx2Idx,
     intg_non_slack_gens_θh_idx2Idx,
     intg_non_gens_θh_idx2Idx) =
         vec_types_intg_θh_idx2Idx


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

    (intg_slack_gens_vh_idx2Idx_in_Idx,
     intg_non_slack_gens_vh_idx2Idx_in_Idx,
     intg_non_gens_vh_idx2Idx_in_Idx) =
         vec_types_intg_vh_idx2Idx_in_Idx


    vec_types_intg_θh_idx2Idx_in_Idx = [
        [ intg_dict_θh_idx2Idx_in_Idx[idx]
          for idx in types_idxs ]
        for types_idxs in
            intg_nodes_type_idxs ]

    (intg_slack_gens_θh_idx2Idx_in_Idx,
     intg_non_slack_gens_θh_idx2Idx_in_Idx,
     intg_non_gens_θh_idx2Idx_in_Idx) =
         vec_types_intg_θh_idx2Idx_in_Idx


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

         intg_vh_θh_id_iq_IDX )
 
          
    # --------------------------------------------------

    return intg_types_Idxs_etc
    
end


function get_dim_dyn_pf_PQ_by_mpc(
    mpc_gen,
    mpc_bus;    
    δ_etc_first =
        false )
    
    #----------------------------------------
    
    (;P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist) =
         get_sta_pf_PQ_param_by_mpc(
             mpc_gen, mpc_bus  )
    
    #----------------------------------------

    dim_P_gens = length(P_gens)
        
    dim_Q_gens = length(P_gens)
        
    dim_P_non_gens = length(P_non_gens)
        
    dim_Q_non_gens = length(Q_non_gens)

    dim_δ = dim_Q_gens

    dim_E_d_dash = dim_Q_gens

    dim_E_q_dash = dim_Q_gens
    
    dim_P_g_loc_load = length(P_g_loc_load)
    
    dim_Q_g_loc_load = length(Q_g_loc_load)


    dim_dyn_pf_PQ_δ_ed_dash_eq_dash =
        (;dim_P_gens,
         dim_Q_gens,
         dim_P_non_gens,
         dim_Q_non_gens,
         dim_δ,
         dim_E_d_dash,
         dim_E_q_dash,
         dim_P_g_loc_load,
         dim_Q_g_loc_load)          

    return(; dim_dyn_pf_PQ_δ_ed_dash_eq_dash )
         
        
end


function get_dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs_by_mpc(
    mpc_gen,
    mpc_bus;    
    δ_etc_first =
        false )
    
    #----------------------------------------
    (;P_gens,
     Q_gens,
     P_non_gens,
     Q_non_gens,
     P_g_loc_load,
     Q_g_loc_load,
     loc_load_exist) =
         get_sta_pf_PQ_param_by_mpc(
             mpc_gen, mpc_bus  )
    
    #----------------------------------------

    dim_P_gens = length(P_gens)
        
    dim_Q_gens = length(P_gens)
        
    dim_P_non_gens = length(P_non_gens)
        
    dim_Q_non_gens = length(Q_non_gens)

    dim_δ = dim_Q_gens

    dim_E_d_dash = dim_Q_gens

    dim_E_q_dash = dim_Q_gens
    
    dim_P_g_loc_load = length(P_g_loc_load)
    
    dim_Q_g_loc_load = length(Q_g_loc_load)
        
    #----------------------------------------
    #----------------------------------------
    
    if δ_etc_first == false

        dim_pf_PQ_δ_etc_para  =
            [dim_P_gens,
             dim_Q_gens,
             dim_P_non_gens,
             dim_Q_non_gens,
             dim_δ,
             dim_E_d_dash,
             dim_E_q_dash,
             dim_P_g_loc_load,
             dim_Q_g_loc_load]

        (dyn_P_gens_Idxs,
         dyn_Q_gens_Idxs,
         dyn_P_non_gens_Idxs,
         dyn_Q_non_gens_Idxs,
         dyn_δ_pf_Idxs,
         dyn_ed_pf_Idxs,
         dyn_eq_pf_Idxs,
         dyn_P_g_loc_load_Idxs,
         dyn_Q_g_loc_load_Idxs) =
             get_vars_or_paras_Idxs_in_flattend(
                 dim_pf_PQ_δ_etc_para;
                 dims_given = true )
        
        dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs =
            ( ;dyn_P_gens_Idxs,
              dyn_Q_gens_Idxs,
              dyn_P_non_gens_Idxs,
              dyn_Q_non_gens_Idxs,
              dyn_δ_pf_Idxs,
              dyn_ed_pf_Idxs,
              dyn_eq_pf_Idxs,
              dyn_P_g_loc_load_Idxs,
              dyn_Q_g_loc_load_Idxs )
        
        return (; dyn_pf_PQ_δ_etc_Idxs, )
        
    else

        dim_pf_PQ_δ_etc_para  =
            [dim_δ,
             dim_E_d_dash,
             dim_E_q_dash,
             dim_P_gens,
             dim_Q_gens,
             dim_P_non_gens,
             dim_Q_non_gens,             
             dim_P_g_loc_load,
             dim_Q_g_loc_load]
        
        (dyn_δ_pf_Idxs,
         dyn_ed_pf_Idxs,
         dyn_eq_pf_Idxs,
         dyn_P_gens_Idxs,
         dyn_Q_gens_Idxs,
         dyn_P_non_gens_Idxs,
         dyn_Q_non_gens_Idxs,
         dyn_P_g_loc_load_Idxs.
         dyn_Q_g_loc_load_Idxs) =
             get_vars_or_paras_Idxs_in_flattend(
                 dim_pf_PQ_δ_etc_para;
                 dims_given = true )

        dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs =
            (; dyn_δ_pf_Idxs,
             dyn_ed_pf_Idxs,
             dyn_eq_pf_Idxs,
             dyn_P_gens_Idxs,
             dyn_Q_gens_Idxs,
             dyn_P_non_gens_Idxs,
             dyn_Q_non_gens_Idxs,
             dyn_P_g_loc_load_Idxs.
             dyn_Q_g_loc_load_Idxs  )
        
        return (; dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs, )
        
    end
    
end


function get_dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs(
    dim_dyn_pf_PQ_δ_ed_dash_eq_dash;    
    δ_etc_first =
        false )

    #----------------------------------------

   (;dim_P_gens,
    dim_Q_gens,
    dim_P_non_gens,
    dim_Q_non_gens,
    dim_δ,
    dim_E_d_dash,
    dim_E_q_dash,
    dim_P_g_loc_load,
    dim_Q_g_loc_load) =
        dim_dyn_pf_PQ_δ_ed_dash_eq_dash           

    #----------------------------------------
    
    if δ_etc_first == false

        dim_pf_PQ_δ_etc_para  =
            [dim_P_gens,
             dim_Q_gens,
             dim_P_non_gens,
             dim_Q_non_gens,
             dim_δ,
             dim_E_d_dash,
             dim_E_q_dash,
             dim_P_g_loc_load,
             dim_Q_g_loc_load]

        (dyn_P_gens_Idxs,
         dyn_Q_gens_Idxs,
         dyn_P_non_gens_Idxs,
         dyn_Q_non_gens_Idxs,
         dyn_δ_pf_Idxs,
         dyn_ed_pf_Idxs,
         dyn_eq_pf_Idxs,
         dyn_P_g_loc_load_Idxs,
         dyn_Q_g_loc_load_Idxs) =
             get_vars_or_paras_Idxs_in_flattend(
                 dim_pf_PQ_δ_etc_para;
                 dims_given = true )
        
        dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs =
            ( ;dyn_P_gens_Idxs,
              dyn_Q_gens_Idxs,
              dyn_P_non_gens_Idxs,
              dyn_Q_non_gens_Idxs,
              dyn_δ_pf_Idxs,
              dyn_ed_pf_Idxs,
              dyn_eq_pf_Idxs,
              dyn_P_g_loc_load_Idxs,
              dyn_Q_g_loc_load_Idxs )
        
        return (; dyn_pf_PQ_δ_etc_Idxs, )
        
    else

        dim_pf_PQ_δ_etc_para  =
            [dim_δ,
             dim_E_d_dash,
             dim_E_q_dash,
             dim_P_gens,
             dim_Q_gens,
             dim_P_non_gens,
             dim_Q_non_gens,             
             dim_P_g_loc_load,
             dim_Q_g_loc_load]
        
        (dyn_δ_pf_Idxs,
         dyn_ed_pf_Idxs,
         dyn_eq_pf_Idxs,
         dyn_P_gens_Idxs,
         dyn_Q_gens_Idxs,
         dyn_P_non_gens_Idxs,
         dyn_Q_non_gens_Idxs,
         dyn_P_g_loc_load_Idxs.
         dyn_Q_g_loc_load_Idxs) =
             get_vars_or_paras_Idxs_in_flattend(
                 dim_pf_PQ_δ_etc_para;
                 dims_given = true )

        dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs =
            (; dyn_δ_pf_Idxs,
             dyn_ed_pf_Idxs,
             dyn_eq_pf_Idxs,
             dyn_P_gens_Idxs,
             dyn_Q_gens_Idxs,
             dyn_P_non_gens_Idxs,
             dyn_Q_non_gens_Idxs,
             dyn_P_g_loc_load_Idxs.
             dyn_Q_g_loc_load_Idxs )
        
        return (; dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs, )
        
    end
    
end



function get_dyn_pf_δ_eq_dash_PQ_Idxs(
    dim_dyn_pf_δ_ed_dash_PQ;    
    δ_etc_first =
        false )

    #----------------------------------------

    (;
     dim_δ,
     dim_E_q_dash,
     
     dim_P_gens,
     dim_Q_gens,
     
     dim_P_non_gens,
     dim_Q_non_gens,

     dim_P_g_loc_load,
     dim_Q_g_loc_load) =
        dim_dyn_pf_δ_ed_dash_PQ           

    #----------------------------------------

    dim_pf_PQ_δ_etc_para  =
        [dim_δ,
         dim_E_q_dash,
         
         dim_P_gens,
         dim_Q_gens,
         
         dim_P_non_gens,
         dim_Q_non_gens,
         
         dim_P_g_loc_load,
         dim_Q_g_loc_load]

    (dyn_δ_pf_Idxs,
     dyn_eq_pf_Idxs,
     
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
     
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
     
     dyn_P_g_loc_load_Idxs.
     dyn_Q_g_loc_load_Idxs) =
         get_vars_or_paras_Idxs_in_flattend(
             dim_dyn_pf_δ_ed_dash_PQ ;
             dims_given = true )

    dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs =
        (; dyn_δ_pf_Idxs,
         dyn_eq_pf_Idxs,
         
         dyn_P_gens_Idxs,
         dyn_Q_gens_Idxs,
         
         dyn_P_non_gens_Idxs,
         dyn_Q_non_gens_Idxs,
         
         dyn_P_g_loc_load_Idxs.
         dyn_Q_g_loc_load_Idxs )

    return (; dyn_pf_PQ_δ_ed_dash_eq_dash_Idxs, )
    
end

#-----------------------------------------------------


function get_static_and_dynamic_pf_var_generic_idx(
    net_nodes_type_idxs )

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
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     nodes_with_demands_idx)  =
         net_nodes_type_idxs


    #-------------------------------

    dim_nodes = length(all_nodes_idx)
    
    dim_vh = dim_θh = dim_nodes
        
    dim_gens = length(gens_nodes_idx)

    dim_P_gens = dim_Q_gens = dim_Pg = dim_Qg = dim_gens

    dim_i_d = dim_i_q = dim_δ = dim_ed = dim_eq =  dim_gens

    dim_non_gens =  length(non_gens_nodes_idx)
    
    dim_P_non_gens = dim_Q_non_gens = dim_non_gens
    
    if loc_load_exist == true

        dim_loc_load = length(gens_nodes_with_loc_loads_idx)
        
        dim_P_g_loc_load =  dim_Q_g_loc_load = dim_loc_load
        
        dim_P_g_loc = dim_Q_g_loc = dim_loc_load
        
    else
        
        dim_loc_load = 0

        dim_P_g_loc_load =  dim_Q_g_loc_load = dim_loc_load
        
        dim_P_g_loc = dim_Q_g_loc = dim_loc_load
    end

    #-------------------------------
        
    δ_ed_dash_eq_dash_Idxs_in_flattend =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_δ,
              dim_ed,
              dim_eq];
            dims_given = true )
    
    #-------------------------------
        
    flat_vh_flat_θh_flat_id_iq_Idx  =
        flat_vh_flat_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_vh,
             dim_θh,
             dim_i_d,
             dim_i_q];
            dims_given = true )

    (nodes_flat_vh_Idx,
     nodes_flat_θh_Idx,
     gens_id_idxs,
     gens_iq_idxs ) =
         flat_vh_flat_θh_flat_id_iq_Idx 

    gens_vh_idxs =
        nodes_flat_vh_Idx[ gens_nodes_idx ]
    
    gens_θh_idxs =
        nodes_flat_θh_Idx[ gens_nodes_idx ]

    non_gens_nodes_vh_idxs =
        nodes_flat_vh_Idx[ non_gens_nodes_idx ]
    
    non_gens_nodes_θh_idxs =
        nodes_flat_vh_Idx[ non_gens_nodes_idx ]
    
    non_pre_ordered_pf_vars_Idxs =
        (;
         gens_vh_idxs,
         gens_θh_idxs,
         
         non_gens_nodes_vh_idxs,
         non_gens_nodes_θh_idxs,

         gens_id_idxs,
         gens_iq_idxs,
         
         flat_vh_flat_θh_flat_id_iq_Idx )
    
    #-------------------------------
        
    flat_vh_flat_θh_flat_Idx  =
        flat_vh_flat_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_vh,
              dim_θh] ;
            dims_given = true )
    

    non_pre_ordered_pf_no_id_iq_vars_Idxs =
        (;
         gens_vh_idxs,
         gens_θh_idxs,
         
         non_gens_nodes_vh_idxs,
         non_gens_nodes_θh_idxs,

         flat_vh_flat_θh_flat_Idx )
    
    #--------------------------------------
        
    gens_vh_idx_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_gens,
              dim_gens ] ;
            dims_given = true )

    
    # flat_Pg_flat_Qg_Idx =
    #     get_vars_or_paras_Idxs_in_flattend(
    #         [ dim_Pg,
    #           dim_Qg ] ;
    #         dims_given = true )
    
    pf_Pg_Qg_Png_Qng_Pll_Qll_Idxs =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_P_gens,
              dim_Q_gens,
              dim_P_non_gens,
              dim_Q_non_gens,
              dim_P_g_loc_load,
              dim_Q_g_loc_load ];
            dims_given = true )


    (Pg_Idxs,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
            pf_Pg_Qg_Png_Qng_Pll_Qll_Idxs

    Pg_Qg_Png_Qng_Pgll_Qgll_Idxs =
        (;Pg_Idxs,
         Qg_Idxs,
         Png_Idxs,
         Qng_Idxs,
         Pgll_Idxs,
         Qgll_Idxs ) 
         
    #--------------------------------------    
            
    dyn_pf_δ_eq_dash_0_P_Q_idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_δ,
             dim_eq,
             dim_P_gens,
             dim_Q_gens,
             dim_P_non_gens,
             dim_Q_non_gens,
             dim_P_g_loc_load,
             dim_Q_g_loc_load ];
             dims_given = true )
    
   ( dyn_δ_Idxs,
     dyn_eq_dash_0_Idxs,
       
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
       
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs                   
    ) =
        dyn_pf_δ_eq_dash_0_P_Q_idx 
 
    dyn_δ_eq_dash_0_pf_P_Q_idx =
        (; dyn_δ_Idxs,
         dyn_eq_dash_0_Idxs,
         
         dyn_P_gens_Idxs,
         dyn_Q_gens_Idxs,
         
         dyn_P_non_gens_Idxs,
         dyn_Q_non_gens_Idxs,
       
         dyn_P_gens_loc_load_Idxs,
         dyn_Q_gens_loc_load_Idxs )

    return (; δ_ed_dash_eq_dash_Idxs_in_flattend,
         non_pre_ordered_pf_vars_Idxs,
         non_pre_ordered_pf_no_id_iq_vars_Idxs,
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         dyn_δ_eq_dash_0_pf_P_Q_idx )
    
end



function get_net_generic_idx(
    plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json )
    
    #--------------------------------------
        
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    #--------------------------------------
    
    (; slack_bus_idx,
     gens_idx,
     slack_gens_nodes_idx,
     non_slack_gens_nodes_idx,
     gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     gens_with_loc_load_idx,
     loc_load_exist,
     load_nodes_idx,
     transmission_nodes_idx,
     non_gens_nodes_idx,
     all_nodes_idx,
     non_slack_gens_and_non_gens_idx,
     nodes_with_demands_idx )  =
         NamedTupleTools.select(
             net_nodes_type_idxs,
             (:slack_bus_idx,
              :gens_idx,
              :slack_gens_nodes_idx,
              :non_slack_gens_nodes_idx,
              :gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :gens_with_loc_load_idx,
              :loc_load_exist,
              :load_nodes_idx,
              :transmission_nodes_idx,
              :non_gens_nodes_idx,
              :all_nodes_idx,
              :non_slack_gens_and_non_gens_idx,
              :nodes_with_demands_idx))

    
    (;n2s_slack_gens_idx,
     n2s_non_slack_gens_idx,
     n2s_gens_idx,
     n2s_non_gens_idx,
     n2s_load_idx,
     n2s_gens_with_loc_load_idxs,
     n2s_transmission_idxs,
     n2s_all_nodes_idx ) =
         NamedTupleTools.select(
             get_dict_net_streamlined_idx_by_nodes_type_idxs(
                 net_nodes_type_idxs ),
             (:n2s_slack_gens_idx,
              :n2s_non_slack_gens_idx,
              :n2s_gens_idx,
              :n2s_non_gens_idx,
              :n2s_load_idx,
              :n2s_gens_with_loc_load_idxs,
              :n2s_transmission_idxs,
              :n2s_all_nodes_idx ))
    
    #-------------------------------


    dim_nodes = length(all_nodes_idx)
    
    dim_vh = dim_θh = dim_nodes
        
    dim_gens = length(gens_nodes_idx)

    dim_P_gens =
        dim_Q_gens = dim_Pg = dim_Qg = dim_gens

    dim_i_d =
        dim_i_q = dim_δ = dim_ed = dim_eq =  dim_gens

    dim_non_gens =
        length(non_gens_nodes_idx)
    
    dim_P_non_gens =
        dim_Q_non_gens = dim_non_gens
    
    if loc_load_exist == true

        dim_loc_load =
            length(gens_nodes_with_loc_loads_idx)
        
        dim_P_g_loc_load =
            dim_Q_g_loc_load = dim_loc_load
        
        dim_P_g_loc =
            dim_Q_g_loc = dim_loc_load
        
    else
        
        dim_loc_load = 0

        dim_P_g_loc_load =
            dim_Q_g_loc_load = dim_loc_load
        
        dim_P_g_loc =
            dim_Q_g_loc = dim_loc_load
    end


    dim_non_gen_PQ = dim_P_non_gens + dim_Q_non_gens +
        dim_P_g_loc_load + dim_Q_g_loc_load
    
    #-------------------------------
        
    δ_ed_dash_eq_dash_Idxs_in_flattend =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_δ,
              dim_ed,
              dim_eq];
            dims_given = true )
    
    #-------------------------------
        
    flat_vh_flat_θh_flat_id_iq_Idx  =
        flat_vh_flat_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_vh,
             dim_θh,
             dim_i_d,
             dim_i_q];
            dims_given = true )
    
    (flat_vh_Idx,
     flat_θh_Idx,
     flat_id_Idx,
     flat_iq_Idx) =
         flat_vh_flat_θh_flat_id_iq_Idx

    flat_vh_flat_θh_flat_id_iq_Idx =
        (;flat_vh_Idx,
         flat_θh_Idx,
         flat_id_Idx,
         flat_iq_Idx)
    

    gens_id_idxs = flat_id_Idx
    
    gens_iq_idxs = flat_iq_Idx 
    
    gens_vh_idxs =
        flat_vh_Idx[ gens_nodes_idx ]
    
    gens_θh_idxs =
        flat_θh_Idx[ gens_nodes_idx ]

    non_gens_nodes_vh_idxs =
        flat_vh_Idx[ non_gens_nodes_idx ]
    
    non_gens_nodes_θh_idxs =
        flat_vh_Idx[ non_gens_nodes_idx ]
    
    non_pre_ordered_pf_vars_Idxs =
        (;
         gens_vh_idxs,
         gens_θh_idxs,
         
         non_gens_nodes_vh_idxs,
         non_gens_nodes_θh_idxs,

         gens_id_idxs,
         gens_iq_idxs,
         
         flat_vh_flat_θh_flat_id_iq_Idx )
    
    #-------------------------------
        
    flat_vh_flat_θh_flat_Idx  =
        flat_vh_flat_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_vh,
              dim_θh] ;
            dims_given = true )
    

    non_pre_ordered_pf_no_id_iq_vars_Idxs =
        (;
         gens_vh_idxs,
         gens_θh_idxs,
         
         non_gens_nodes_vh_idxs,
         non_gens_nodes_θh_idxs,

         flat_vh_flat_θh_flat_Idx )
    
    #--------------------------------------
        
    gens_vh_idx_θh_idx_in_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_gens,
              dim_gens ] ;
            dims_given = true )
    
    pf_Pg_Qg_Png_Qng_Pll_Qll_Idxs =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_P_gens,
              dim_Q_gens,
              dim_P_non_gens,
              dim_Q_non_gens,
              dim_P_g_loc_load,
              dim_Q_g_loc_load ];
            dims_given = true )


    (Pg_Idxs,
     Qg_Idxs,
     Png_Idxs,
     Qng_Idxs,
     Pgll_Idxs,
     Qgll_Idxs ) =
            pf_Pg_Qg_Png_Qng_Pll_Qll_Idxs

    Pg_Qg_Png_Qng_Pgll_Qgll_Idxs =
        (;Pg_Idxs,
         Qg_Idxs,
         Png_Idxs,
         Qng_Idxs,
         Pgll_Idxs,
         Qgll_Idxs ) 
         
    #--------------------------------------
    
    dyn_pf_δ_eq_dash_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_δ,
             dim_eq];
             dims_given = true )
    
   ( dyn_pf_δ_Idxs,
     dyn_pf_eq_dash_Idxs ) =
         dyn_pf_δ_eq_dash_Idx

    dyn_pf_δ_eq_dash_Idx =
        (; dyn_pf_δ_Idxs,
          dyn_pf_eq_dash_Idxs )
             
    #--------------------------------------
    
    dyn_pf_δ_ed_dash_eq_dash_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_δ,
             dim_ed,
             dim_eq];
             dims_given = true )
    
    (dyn_pf_δ_Idxs,
     dyn_pf_ed_dash_Idxs,
     dyn_pf_eq_dash_Idxs ) =
         dyn_pf_δ_ed_dash_eq_dash_Idx

    dyn_pf_δ_ed_dash_eq_dash_Idx =
        (; dyn_pf_δ_Idxs,
         dyn_pf_ed_dash_Idxs,
         dyn_pf_eq_dash_Idxs )

    #--------------------------------------    
    
    dyn_pf_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_P_non_gens,
             dim_Q_non_gens,
             dim_P_g_loc_load,
             dim_Q_g_loc_load ];
             dims_given = true )
    
   ( dyn_pf_P_non_gens_Idxs,
     dyn_pf_Q_non_gens_Idxs,
       
     dyn_pf_P_gens_loc_load_Idxs,
     dyn_pf_Q_gens_loc_load_Idxs) =
         dyn_pf_Png_Qng_Pll_Qll_Idx
         

    dyn_pf_Png_Qng_Pll_Qll_Idx =
        (;dyn_pf_P_non_gens_Idxs,
          dyn_pf_Q_non_gens_Idxs,
          
         dyn_pf_P_gens_loc_load_Idxs,
         dyn_pf_Q_gens_loc_load_Idxs)    
 
    #--------------------------------------
    
    # dyn_pf_δ_eq_dash_0_P_Q_idx
    
    dyn_δ_eq_dash_0_pf_P_Q_idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_δ,
             dim_eq,
             dim_P_gens,
             dim_Q_gens,
             dim_P_non_gens,
             dim_Q_non_gens,
             dim_P_g_loc_load,
             dim_Q_g_loc_load ];
             dims_given = true )
    
   ( dyn_δ_Idxs,
     dyn_eq_dash_0_Idxs,
       
     dyn_P_gens_Idxs,
     dyn_Q_gens_Idxs,
       
     dyn_P_non_gens_Idxs,
     dyn_Q_non_gens_Idxs,
       
     dyn_P_gens_loc_load_Idxs,
     dyn_Q_gens_loc_load_Idxs  ) =
         dyn_δ_eq_dash_0_pf_P_Q_idx
         
    dyn_pf_δ_eq_dash_0_P_Q_idx =
        (; dyn_δ_Idxs,
         dyn_eq_dash_0_Idxs,
         
         dyn_P_gens_Idxs,
         dyn_Q_gens_Idxs,
         
         dyn_P_non_gens_Idxs,
         dyn_Q_non_gens_Idxs,
       
         dyn_P_gens_loc_load_Idxs,
         dyn_Q_gens_loc_load_Idxs )

    #--------------------------------------

    dyn_pf_idq_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_i_d,
             dim_i_q];
            dims_given = true )

    (dyn_pf_id_Idx,
     dyn_pf_iq_Idx) =
         dyn_pf_idq_Idx

    dyn_pf_idq_Idx =
        (;dyn_pf_id_Idx,
         dyn_pf_iq_Idx)

    #--------------------------------------
    
    dyn_pf_idq_δ_ed_dash_non_gen_PQ_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [2 * dim_gens,
             2 * dim_gens,
             dim_non_gen_PQ ];
            dims_given = true )

    ( dyn_pf_flat_idq_Idx,
      dyn_pf_flat_δ_ed_dash_Idx,
      dyn_pf_flat_non_gen_PQ_Idx) =
          dyn_pf_idq_δ_ed_dash_non_gen_PQ_Idx

    dyn_pf_idq_δ_ed_dash_non_gen_PQ_Idx =
        (;dyn_pf_flat_idq_Idx,
         dyn_pf_flat_δ_ed_dash_Idx,
         dyn_pf_flat_non_gen_PQ_Idx)
    
    #--------------------------------------
    
    dyn_pf_idq_δ_ed_dash_ed_dash_non_gen_PQ_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [2 * dim_gens,
             3 * dim_gens,
             dim_non_gen_PQ ];
            dims_given = true )

    (dyn_pf_flat_idq_Idx,
     dyn_pf_flat_δ_ed_dash_eq_dash_Idx,
     dyn_pf_flat_non_gen_PQ_Idx) =
          dyn_pf_idq_δ_ed_dash_ed_dash_non_gen_PQ_Idx

    dyn_pf_idq_δ_ed_dash_ed_dash_non_gen_PQ_Idx =
        (;dyn_pf_flat_idq_Idx,
         dyn_pf_flat_δ_ed_dash_eq_dash_Idx,
         dyn_pf_flat_non_gen_PQ_Idx)    
    
    #--------------------------------------
       
    dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_δ,
             dim_eq,
             dim_P_non_gens,
             dim_Q_non_gens,
             dim_P_g_loc_load,
             dim_Q_g_loc_load ];
             dims_given = true )
    
   ( dyn_pf_δ_Idxs,
     dyn_pf_eq_dash_Idxs,
       
     dyn_pf_P_gens_loc_load_Idxs,
     dyn_pf_Q_gens_loc_load_Idxs  ) =
         dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx
         
    dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx =
        (; dyn_pf_δ_Idxs,
         dyn_pf_eq_dash_Idxs,

         dyn_pf_P_non_gens_Idxs,
         dyn_pf_Q_non_gens_Idxs,
       
         dyn_pf_P_gens_loc_load_Idxs,
         dyn_pf_Q_gens_loc_load_Idxs )

    
    #-------------------------------

    dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_P_gens,
             dim_P_non_gens,
             dim_Q_gens,
             dim_Q_non_gens,
             dim_i_d,
             dim_i_q];
            dims_given = true )
    
    (ΔPg_Idx,
     ΔPng_Idx,
     ΔQg_Idx,
     ΔQng_Idx,
     Δi_d_Idx,
     Δi_q_Idx) =
         dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx

    dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx =
        (; ΔPg_Idx,
         ΔPng_Idx,
         ΔQg_Idx,
         ΔQng_Idx,
         Δi_d_Idx,
         Δi_q_Idx)
    
    #--------------------------------------
    
    dyn_pf_fun_kwd_n2s_idxs =
        (;
         n2s_slack_gens_idx,
         n2s_non_slack_gens_idx,
         n2s_gens_idx,
         n2s_non_gens_idx,
         n2s_gens_with_loc_load_idxs,
         n2s_all_nodes_idx ) 
            
    #-------------------------------

    dyn_pf_fun_kwd_net_idxs =
        (;
         slack_gens_nodes_idx,
         non_slack_gens_nodes_idx,
         gens_nodes_idx,
         non_gens_nodes_idx,
         gens_nodes_with_loc_loads_idx,
         gens_with_loc_load_idx,
         all_nodes_idx,
         
         non_slack_gens_and_non_gens_idx,
         nodes_with_demands_idx)
    
    #-------------------------------
    
    flux_decay_model_state_vars_syms =
        [:δ, :ω, :eq_dash, :E_fd]
    
    flux_decay_model_state_vars_idx_syms =
        [:δ_idx_in_state,
         :ω_idx_in_state,
         :eq_dash_idx_in_state,
         :E_fd_idx_in_state ]
    
    flux_decay_model_states_Idx =
        get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
            flux_decay_model_state_vars_syms,
            gens_nodes_idx )

    flux_decay_model_states_comp_idxs_in_Idx =
         get_vars_idxs_in_range_Idxs(
             flux_decay_model_states_Idx,
             flux_decay_model_state_vars_idx_syms )

    # (;δ_idx_in_state,
    #  ω_idx_in_state,
    #  eq_dash_idx_in_state,
    #  E_fd_idx_in_state ) =
    #      flux_decay_model_states_comp_idxs_in_Idx

    flux_decay_model_vars_Idx_in_state =
        get_vars_or_paras_Idxs_in_flattend(
            [sum(length.(flux_decay_model_states_Idx)),
             dim_nodes,
             dim_nodes];
            dims_given = true )
    
    # flux_decay_model_states_comp_idxs_in_Idx
    
    (state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state ) =
         flux_decay_model_vars_Idx_in_state
    
    flux_decay_model_vars_Idx_in_state =
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state )
    
    #--------------------------------------

   flux_decay_model_vars_wt_i_dq_Idx_in_state =
        get_vars_or_paras_Idxs_in_flattend(
            [sum(length.(flux_decay_model_states_Idx)),
             dim_nodes,
             dim_nodes,
             dim_gens,
             dim_gens];
            dims_given = true )
    
    # flux_decay_model_states_comp_idxs_in_Idx
    
    (state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state) =
         flux_decay_model_vars_wt_i_dq_Idx_in_state 
    
    flux_decay_model_vars_wt_i_dq_Idx_in_state =
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state)
    
    #--------------------------------------
    
    dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_P_non_gens,
             dim_Q_non_gens,
             dim_P_g_loc_load,
              dim_Q_g_loc_load];
            dims_given = true )

    (dyn_V_ref_Idx,
     dyn_Tm_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx)=
         dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx
    
    dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx =
        (;dyn_V_ref_Idx,
         dyn_Tm_Idx,
         dyn_Png_Idx,
         dyn_Qng_Idx,
         dyn_Pll_Idx,
         dyn_Qll_Idx )
    
    #--------------------------------------

    ode_vh_id_iq_V_ref_Tm_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    (ode_vh_Idx,
     ode_id_Idx,
     ode_iq_Idx,
     ode_V_ref_Idx,
     ode_Tm_Idx)=
         ode_vh_id_iq_V_ref_Tm_Idx
    
    ode_vh_id_iq_V_ref_Tm_Idx =
        (;ode_vh_Idx,
         ode_id_Idx,
         ode_iq_Idx,
         ode_V_ref_Idx,
         ode_Tm_Idx )

    
    #--------------------------------------

    id_iq_pg_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    (id_Idx,
     iq_Idx,
     pg_Idx,
     vh_Idx)=
         id_iq_pg_vh_Idx
    
    id_iq_pg_vh_Idx =
        (;id_Idx,
         iq_Idx,
         pg_Idx,
         vh_Idx )
    
    #--------------------------------------

    ωs_ωref0_vref0_porder0_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    (ωs_Idx,
     ωref0_Idx,
     vref0_Idx,
     porder0_Idx)=
         ωs_ωref0_vref0_porder0_Idx
    
    ωs_ωref0_vref0_porder0_Idx =
        (;ωs_Idx,
         ωref0_Idx,
         vref0_Idx,
         porder0_Idx )

    #--------------------------------------

    ωref0_vref0_porder0_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [ dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    ( ωref0_Idx,
     vref0_Idx,
     porder0_Idx)=
         ωref0_vref0_porder0_Idx
    
    ωref0_vref0_porder0_Idx =
        (; ωref0_Idx,
         vref0_Idx,
         porder0_Idx )

    #--------------------------------------

    ωref0_vref0_porder0_id_iq_vh_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens,
             dim_gens];
            dims_given = true )

    (ωref0_Idx,
     vref0_Idx,
      porder0_Idx,
      id_Idx,
      iq_Idx,
      vh_Idx)=
         ωref0_vref0_porder0_id_iq_vh_Idx
    
    ωref0_vref0_porder0_id_iq_vh_Idx =
        (;ωref0_Idx,
         vref0_Idx,
         porder0_Idx,
         id_Idx,
         iq_Idx,
         vh_Idx )

    #--------------------------------------        
    #--------------------------------------
    
    dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx =
        get_vars_or_paras_Idxs_in_flattend(
            [dim_gens,
             dim_gens,
             dim_gens,
             dim_P_non_gens,
             dim_Q_non_gens,
             dim_P_g_loc_load,
             dim_Q_g_loc_load];
            dims_given = true )

    (dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx)=
         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx
    
    dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx =
        (;dyn_ω_ref_Idx,
         dyn_v_ref_Idx,
         dyn_p_order_Idx,
         dyn_Png_Idx,
         dyn_Qng_Idx,
         dyn_Pll_Idx,
         dyn_Qll_Idx )
    
    #--------------------------------------
    # generic sta
    #--------------------------------------
        
    SC_generic_model_state_vars_syms =
        [:δ, :ω, :ed_dash, :eq_dash, :E_fd, :R_f, :V_R]
    
    SC_generic_model_state_vars_idx_syms =
        [:δ_idx_in_state,
         :ω_idx_in_state,
         :ed_dash_idx_in_state,
         :eq_dash_idx_in_state,
         :E_fd_idx_in_state,
         :R_f_idx_in_state,
         :V_R_idx_in_state]
    
    SC_generic_model_states_idx_in_state_Idx =
        get_idxs_in_flattened_by_nodes_idx_wt_vars_syms(
            SC_generic_model_state_vars_syms,
            gens_nodes_idx )

    SC_generic_model_states_comp_idxs_in_Idx =
         get_vars_idxs_in_range_Idxs(
             SC_generic_model_states_idx_in_state_Idx,
             SC_generic_model_state_vars_idx_syms )

    # (δ_idx_in_state,
    #  ω_idx_in_state,
    #  ed_dash_idx_in_state,
    #  eq_dash_idx_in_state,
    #  E_fd_idx_in_state,
    #  R_f_idx_in_state,
    #  V_R_idx_in_state ) =
    #      SC_generic_model_states_comp_idxs_in_Idx
    
    SC_generic_model_vars_no_i_dq_Idx_in_state =
        get_vars_or_paras_Idxs_in_flattend(
            [sum(length.(
                SC_generic_model_states_idx_in_state_Idx)),
             dim_nodes,
             dim_nodes];
            dims_given = true )
    
    # flux_decay_model_states_comp_idxs_in_Idx
    
    (state_var_no_i_dq_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state ) =
         SC_generic_model_vars_no_i_dq_Idx_in_state
    
    SC_generic_model_vars_no_i_dq_Idx_in_state =
        (;state_var_no_i_dq_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state )
    
    #--------------------------------------

   SC_generic_model_vars_wt_i_dq_Idx_in_state =
        get_vars_or_paras_Idxs_in_flattend(
            [sum(length.(
                SC_generic_model_states_idx_in_state_Idx)),
             dim_nodes,
             dim_nodes,
             dim_gens,
             dim_gens];
            dims_given = true )
    
    # flux_decay_model_states_comp_idxs_in_Idx
    
    (state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state) =
         SC_generic_model_vars_wt_i_dq_Idx_in_state 
    
    SC_generic_model_vars_wt_i_dq_Idx_in_state =
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state)
    
    #--------------------------------------
    #--------------------------------------
    
    net_generic_idx =
        (;loc_load_exist,
         δ_ed_dash_eq_dash_Idxs_in_flattend,
         non_pre_ordered_pf_vars_Idxs,
         non_pre_ordered_pf_no_id_iq_vars_Idxs,
         Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
         # dyn_δ_eq_dash_0_pf_P_Q_idx,
         dyn_pf_δ_eq_dash_0_P_Q_idx,
         dyn_pf_fun_kwd_n2s_idxs,
         dyn_pf_fun_kwd_net_idxs,

         dyn_pf_δ_eq_dash_Idx,
         dyn_pf_δ_ed_dash_eq_dash_Idx,
         dyn_pf_Png_Qng_Pll_Qll_Idx,

         
         dyn_pf_idq_Idx,
         dyn_pf_idq_δ_ed_dash_non_gen_PQ_Idx,
         dyn_pf_idq_δ_ed_dash_ed_dash_non_gen_PQ_Idx,


         dyn_pf_δ_eq_dash_Png_Qng_Pll_Qll_Idx,
         
         dyn_pf_Δ_Pg_Png_Qg_Qng_id_iq_Idx,

         flux_decay_model_states_Idx,
         flux_decay_model_states_comp_idxs_in_Idx,

         flux_decay_model_vars_Idx_in_state,
         flux_decay_model_vars_wt_i_dq_Idx_in_state,

         dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,

         ode_vh_id_iq_V_ref_Tm_Idx,
         
         id_iq_pg_vh_Idx,
         ωs_ωref0_vref0_porder0_Idx,

         SC_generic_model_states_idx_in_state_Idx,
         SC_generic_model_states_comp_idxs_in_Idx,
         SC_generic_model_vars_no_i_dq_Idx_in_state,
         SC_generic_model_vars_wt_i_dq_Idx_in_state,

         ωref0_vref0_porder0_Idx,
         ωref0_vref0_porder0_id_iq_vh_Idx,
         dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx)

    return net_generic_idx
    
end




function get_net_generic_parameters(
    plant_generators_data_from_json,
    plant_loads_data_from_json,
    plant_transmission_data_from_json,
    edge_data_from_json,
    shunt_data_from_json,
    baseMVA_data_from_json;
    
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,
    
    in_components_type_sym =
        false,    
    ode_gens_generic_sequence_order =
        (:components_data, :gen),    
    ode_gens_generic_selections =
        (:H, :D,
         :ra, :xℓ,
         :X_d, :X_q,
         :X_d_dash,  :X_q_dash,
         :X_d_2dash, :X_q_2dash,
         :T_d_dash,  :T_q_dash,
         :T_d_2dash, :T_q_2dash,
         :Sn),
    ode_gens_para_sequence_order =
        (:components_data, :gen),
    ode_gens_para_selections  =
        (:H, :D,
         :X_d, :X_q,                  
         :X_d_dash, :X_q_dash,
         :T_d_dash, :T_q_dash, :Sn ),
    ode_spcm_gens_para_sequence_order =
        (:components_data, :gen),    
    ode_spcm_gens_para_selections =
        (:H, :X_d_dash ),    
    ode_flux_decay_gens_para_sequence_order =
        (:components_data, :gen),    
    ode_flux_decay_gens_para_selections =
        (:H,
         :X_d, :X_q,                  
         :X_d_dash, :T_d_dash ), 
    ode_flux_decay_avrs_para_sequence_order =
        ( :components_data, :avr ),
    ode_flux_decay_avrs_para_selections =
        (:Ka, :Ta ),
    govs_and_avrs_sequence_order =
        ( :components_data,),
    govs_and_avrs_selections =
        ( :gov, :avr ) )
    
    #--------------------------------------
 
    # (;plant_generators_data_from_json,
    #  plant_loads_data_from_json,
    #  plant_transmission_data_from_json,
    #  edge_data_from_json,
    #  shunt_data_from_json,
    #  baseMVA_data_from_json,
    #  gencost_data_from_json) =
    #  get_net_data_by_components_from_json_file(
    #          net_data_by_components_file;
    #      in_components_type_sym =
    #          in_components_type_sym )

    #--------------------------------------

    baseMVA = baseMVA_data_from_json
    
    #--------------------------------------
    
    ode_gens_generic_para =
         get_selected_comps_ode_para_by_json(
             plant_generators_data_from_json;
             sequence_order =
                 ode_gens_generic_sequence_order,
             selections =
                 ode_gens_generic_selections )

    # modified by AAY ePowerSim
    
    ode_gens_generic_para = namedtuple(
        OrderedDict(a_sym => a_value
             for (a_sym, a_value) in
                    zip(ode_gens_generic_selections,
                        ode_gens_generic_para)))

    # (H,
    #  D,
    #  ra,
    #  xℓ,
    #  X_d,
    #  X_q,
     
    #  X_d_dash,
    #  X_q_dash,
     
    #  X_d_2dash,
    #  X_q_2dash,
    #  T_d_dash,
    #  T_q_dash, Sn) =
    #      NamedTupleTools.select(
    #          ode_gens_generic_para,
    #          (:H,
    #           :D,
    #           :ra,
    #           :xℓ,
    #           :X_d,
    #           :X_q,

    #           :X_d_dash,
    #           :X_q_dash,

    #           :X_d_2dash,
    #           :X_q_2dash,
    #           :T_d_dash,
    #           :T_q_dash, :Sn))


    # ode_gens_generic_para =
    #     (;
    #      H,
    #      D,
    #      ra,
    #      xℓ,
    #      X_d,
    #      X_q,

    #      X_d_dash,
    #      X_q_dash,

    #      X_d_2dash,
    #      X_q_2dash,
    #      T_d_dash,
    #      T_q_dash, Sn) 


    #------------------------------------------------
    
    generic_gens_para =
        get_components_properties_by_json(
            plant_generators_data_from_json;
            sequence_order =
                 ode_gens_generic_sequence_order ,
             selections =
                 ode_gens_generic_selections )
    

    "To make sure Vector{NamedTuple} is returned
     instead of Vector{Any}"
    generic_gens_para =
        NamedTuple[
            item for item in
                generic_gens_para]

    #------------------------------------------------
    
    ode_gens_para =
        get_selected_comps_ode_para_by_json(
             plant_generators_data_from_json;
             sequence_order =
                 ode_gens_para_sequence_order,
            selections =
                ode_gens_para_selections )


    # modified by AAY ePowerSim
    
    ode_gens_para = namedtuple(
        OrderedDict(a_sym => a_value
             for (a_sym, a_value) in
                    zip(ode_gens_para_selections,
                        ode_gens_para)))
    
    # (H,
    #  D,
    #  X_d,
    #  X_q,
    #  X_d_dash,
    #  X_q_dash,
    #  T_d_dash,
    #  T_q_dash, Sn ) =
    #      NamedTupleTools.select(
    #          ode_gens_para,
    #          (:H,
    #           :D,
    #           :X_d,
    #           :X_q,
    #           :X_d_dash,
    #           :X_q_dash,
    #           :T_d_dash,
    #           :T_q_dash  :Sn))

    # ode_gens_para =
    #     (
    #     ;H,
    #     D,
    #     X_d,
    #     X_q,
    #     X_d_dash,
    #     X_q_dash,
    #     T_d_dash,
    #     T_q_dash, Sn )  
    
    #------------------------------------------------
    
     ode_spcm_gens_para =
        get_selected_comps_ode_para_by_json(
             plant_generators_data_from_json;
            sequence_order =
                ode_spcm_gens_para_sequence_order,
            selections =
                ode_spcm_gens_para_selections
               )


    # modified by AAY ePowerSim

    ode_spcm_gens_para = namedtuple(
        OrderedDict(a_sym => a_value
             for (a_sym, a_value) in
                    zip(ode_spcm_gens_para_selections,
                        ode_spcm_gens_para)))
    
    # ( H,
    #   Xd_dash ) =  ode_spcm_gens_para
    
    # ode_spcm_gens_para =
    #     ( ;H,
    #       Xd_dash )
    
    #------------------------------------------------
    
    ode_flux_decay_gens_para =
        get_selected_comps_ode_para_by_json(
             plant_generators_data_from_json;
          sequence_order =
              ode_flux_decay_gens_para_sequence_order,
          selections =
              ode_flux_decay_gens_para_selections )


    # modified by AAY ePowerSim

     ode_flux_decay_gens_para = namedtuple(
        OrderedDict(a_sym => a_value
             for (a_sym, a_value) in
                    zip(ode_flux_decay_gens_para_selections,
                         ode_flux_decay_gens_para)))
    
    # ( H, X_d, X_q, X_d_dash, T_d_dash ) =
    #     ode_flux_decay_gens_para
    
    # ode_flux_decay_gens_para =
    #     (; H, X_d, X_q, X_d_dash, T_d_dash ) 
            
    #------------------------------------------------
    
    ode_flux_decay_avrs_para =
          get_selected_comps_ode_para_by_json(
              plant_generators_data_from_json;
              sequence_order =
                  ode_flux_decay_avrs_para_sequence_order,
              selections =
                  ode_flux_decay_avrs_para_selections )


    # modified by AAY ePowerSim

     ode_flux_decay_avrs_para = namedtuple(
        OrderedDict(a_sym => a_value
             for (a_sym, a_value) in
                    zip(ode_flux_decay_avrs_para_selections,
                         ode_flux_decay_avrs_para)))
   
   # (Ka,
   #  Ta) = ode_flux_decay_avrs_para

   
   # ode_flux_decay_avrs_para = (; Ka,
   #  Ta)  
    
    generic_govs_para, generic_avrs_para =
        get_selected_comps_ode_para_by_json(
            plant_generators_data_from_json;
            sequence_order =
                govs_and_avrs_sequence_order ,
            selections =
                govs_and_avrs_selections )

    #------------------------------------------------
    
    # edges_orientation =
    #     get_edges_orientation_by_generic(
    #         get_edges_fbus_tbus_by_json(
    #             edge_data_from_json)... )

    # edges_Ybr =
    #     get_edges_Ybr_by_generic(
    #         get_edges_generic_data_by_json(
    #             edge_data_from_json)...,
    #         get_nodes_shunts_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA = baseMVA, basekV = basekV )

    # nodes_Yshunt =
    #     get_nodes_Yshunt_by_generic(
    #         get_nodes_shunts_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA =
    #             baseMVA )

    # nodes_idx_and_Yshunt =
    #     get_nodes_idx_and_Yshunt_by_generic(
    #         get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA =
    #             baseMVA )
    
    # Ybr_cal_and_edges_orientation =
    #     get_edges_Ybr_cal_and_edges_orientation_by_generic(
    #         get_edges_fbus_tbus_by_json(
    #             edge_data_from_json)...,
    #         get_edges_generic_data_by_json(
    #             edge_data_from_json )...,
    #         get_nodes_shunts_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA = baseMVA, basekV = basekV )

    # # (; edges_Ybr_cal,
    # #  edges_orientation ) =
    # #      Ybr_cal_and_edges_orientation

    # Ynet_wt_nodes_idx_wt_adjacent_nodes =
    #   get_Ynet_wt_nodes_idx_wt_adjacent_nodes_by_generic(
    #         get_edges_ftbus_and_generic_data_by_json(
    #             edge_data_from_json )...,
    #         get_nodes_shunts_idx_and_Gs_and_Bs_by_json(
    #             shunt_data_from_json)...;
    #         baseMVA = baseMVA, basekV = basekV )

    # # (; 
    # #  Ynet,
    # #  nodes_idx_with_adjacent_nodes_idx ) =
    # #      Ynet_wt_nodes_idx_wt_adjacent_nodes

    # -------------------------------------

    (;edges_orientation,
     edges_Ybr,
     nodes_Yshunt,
     nodes_idx_and_Yshunt,
     Ybr_cal_and_edges_orientation,
     Ynet_wt_nodes_idx_wt_adjacent_nodes ) =
         NamedTupleTools.select(
             get_transmission_network_parameters_by_json(
                 plant_generators_data_from_json,
                 plant_loads_data_from_json,
                 plant_transmission_data_from_json,
                 edge_data_from_json,
                 shunt_data_from_json;
                 baseMVA = baseMVA,
                 basekV  = basekV,
                 use_pu_in_PQ = use_pu_in_PQ,
                 line_data_in_pu = line_data_in_pu ),
             (:edges_orientation,
              :edges_Ybr,
              :nodes_Yshunt,
              :nodes_idx_and_Yshunt,
              :Ybr_cal_and_edges_orientation,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes) )

    # -------------------------------------
    
    dyn_pf_mismatch_vars_kwd_para =
        (; Ynet_wt_nodes_idx_wt_adjacent_nodes,
          ode_gens_para )
    
    # -------------------------------------

    sta_pf_PQ_para =
        get_pf_PQ_param_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json;
            baseMVA =
                baseMVA,
            use_pu_in_PQ =
                use_pu_in_PQ)

    # (;P_gens,
    #  Q_gens,
    #  P_non_gens,
    #  Q_non_gens,
    #  P_g_loc_load,
    #  Q_g_loc_load,
    #  loc_load_exist) =
    #      sta_pf_PQ_para

    gens_vh_slack_θh_para =
        get_gens_vh_slack_θh_para_by_json(
            plant_generators_data_from_json )

    # (; slack_gens_vh,
    #  slack_gens_θh,
    #  gens_vh,
    #  non_slack_gens_vh ) =
    #      gens_vh_slack_θh_para

    sta_pf_vars_and_paras_idx =
        get_sta_pf_vars_and_paras_idx_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )

    # (;
    #  red_types_Idxs_etc,
    #  PQ_sta_para_Idxs,
    #  nodes_types_idxs,
    #  n2s_idxs ) =
    #      sta_pf_vars_and_paras_idx

    pf_sta_ΔPQ_mismatch_parameters =
        get_pf_sta_ΔPQ_mismatch_parameters_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json;
            baseMVA = baseMVA,
            basekV = basekV,
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu)

    # (;
    #  pf_kw_para,
    #  pf_PQ_param,
    #  red_types_Idxs_etc,
    #  net_para  ) =
    #      pf_sta_ΔPQ_mismatch_parameters

    return (;baseMVA,
            ode_gens_generic_para,
            ode_gens_para,
            ode_spcm_gens_para,
            ode_flux_decay_gens_para,
            
            ode_flux_decay_avrs_para,
            
            edges_orientation,
            edges_Ybr,
            nodes_Yshunt,
            nodes_idx_and_Yshunt,
            Ybr_cal_and_edges_orientation,
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            sta_pf_PQ_para,
            gens_vh_slack_θh_para,
            sta_pf_vars_and_paras_idx,
            pf_sta_ΔPQ_mismatch_parameters,
            dyn_pf_mismatch_vars_kwd_para,

            generic_gens_para,
            generic_avrs_para,
            generic_govs_para )

    #--------------------------------------
    
end


function get_net_generic_parameters_and_idx(
    net_data_by_components_file;
        
    basekV = 1.0,    
    use_pu_in_PQ = true,
    line_data_in_pu = true,

    in_components_type_sym =
        false,
    ode_gens_generic_sequence_order =
        (:components_data, :gen),    
    ode_gens_generic_selections =
        (:H, :D,
         :ra, :xℓ,
         :X_d, :X_q,
         :X_d_dash,  :X_q_dash,
         :X_d_2dash, :X_q_2dash,
         :T_d_dash,  :T_q_dash,
         :T_d_2dash, :T_q_2dash,
         :Sn),
    ode_gens_para_sequence_order =
        (:components_data, :gen),
    ode_gens_para_selections  =
        (:H, :D,
         :X_d, :X_q,                  
         :X_d_dash, :X_q_dash,
         :T_d_dash, :T_q_dash, :Sn ),
    ode_spcm_gens_para_sequence_order =
        (:components_data, :gen),    
    ode_spcm_gens_para_selections =
        (:H, :Xd_dash ),    
    ode_flux_decay_gens_para_sequence_order =
        (:components_data, :gen),    
    ode_flux_decay_gens_para_selections =
        (:H,
         :X_d, :X_q,                  
         :X_d_dash, :T_d_dash ), 
    ode_flux_decay_avrs_para_sequence_order =
        ( :components_data, :avr ),
    ode_flux_decay_avrs_para_selections =
        (:Ka, :Ta ),
    govs_and_avrs_sequence_order =
        ( :components_data,),
    govs_and_avrs_selections =
        ( :gov, :avr ) )

    data_by_components =
        get_net_data_by_components_from_json_file(
             net_data_by_components_file;
            in_components_type_sym =
                in_components_type_sym )
    
    (;
     plant_generators_data_from_json,
     plant_loads_data_from_json,
     plant_transmission_data_from_json,
     edge_data_from_json,
     shunt_data_from_json,
     baseMVA_data_from_json,
     gencost_data_from_json) =
         data_by_components
    
    net_nodes_type_idxs =
        get_net_nodes_type_idxs_by_json(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )
    
    n2s_idxs =
         get_dict_net_streamlined_idx_by_nodes_type_idxs(
             net_nodes_type_idxs )

    net_generic_idx  =
        get_net_generic_idx(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json )
    
    # (;δ_ed_dash_eq_dash_Idxs_in_flattend,
    #  non_pre_ordered_pf_vars_Idxs,
    #  non_pre_ordered_pf_no_id_iq_vars_Idxs,
    #  Pg_Qg_Png_Qng_Pgll_Qgll_Idxs,
    #  dyn_δ_eq_dash_0_pf_P_Q_idx,
    #  dyn_pf_fun_kwd_n2s_idxs,
    #  dyn_pf_fun_kwd_net_idxs) =
    #      net_generic_idx
    
    # net_generic_parameters
    #      get_net_generic_parameters(
    #          net_data_by_components_file;
    #          in_components_type_sym =
    #              in_components_type_sym )

    net_generic_parameters =
        get_net_generic_parameters(
            plant_generators_data_from_json,
            plant_loads_data_from_json,
            plant_transmission_data_from_json,
            edge_data_from_json,
            shunt_data_from_json,
            baseMVA_data_from_json;
                
            basekV = basekV,    
            use_pu_in_PQ = use_pu_in_PQ,
            line_data_in_pu = line_data_in_pu,

            in_components_type_sym =
                in_components_type_sym,

            ode_gens_generic_sequence_order =
                ode_gens_generic_sequence_order,
            ode_gens_generic_selections =
                ode_gens_generic_selections,

            ode_gens_para_sequence_order =
                ode_gens_para_sequence_order,
            ode_gens_para_selections =
                ode_gens_para_selections,

            ode_spcm_gens_para_sequence_order =
                ode_spcm_gens_para_sequence_order,
            ode_spcm_gens_para_selections =
                ode_spcm_gens_para_selections,

            ode_flux_decay_gens_para_sequence_order =
                ode_flux_decay_gens_para_sequence_order,
            ode_flux_decay_gens_para_selections =
                ode_flux_decay_gens_para_selections,

            ode_flux_decay_avrs_para_sequence_order = 
                ode_flux_decay_avrs_para_sequence_order,
            ode_flux_decay_avrs_para_selections =
                ode_flux_decay_avrs_para_selections,

            govs_and_avrs_sequence_order =
                govs_and_avrs_sequence_order,
            govs_and_avrs_selections =
                govs_and_avrs_selections )

    return (;data_by_components,            
            net_nodes_type_idxs,
            n2s_idxs,
            net_generic_idx,
            net_generic_parameters )

end


#-----------------------------------------------------
#-----------------------------------------------------

