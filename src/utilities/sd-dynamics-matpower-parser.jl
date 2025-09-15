# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

#################################################
# -----------------------------------------------
# MATPOWER Parser 
# -----------------------------------------------
#################################################

"""
`extract_matpower_case_name`

Returns the case name of a MATPOWER file.

The first line of every matpower file contains function name
for the file, e.g. "function mpc = case5".
"""
function extract_matpower_case_name(
    string::AbstractString)

    if occursin("function", string)

        case_function_parts =
            split(split(string, ';')[1], '=')
        @assert(length(case_function_parts) == 2)
        

        return strip(case_function_parts[2])

    else
        return nothing

    end
    
end


"""
Test drives `extract_matpower_case_name`
"""
function driver_extract_matpower_case_name(
    data_string::String)
    
    data_lines = split(data_string, '\n')
    
    case_name = nothing
    
    last_index = length(data_lines)
    
    index = 1

    a_fun_found = false
    
    while index <= last_index &&  !a_fun_found
        
        line = strip(data_lines[index])
        
        line = "$(line)"

        if length(line) <= 0 || strip(line)[1] == '%'
            
            index = index + 1
            
            continue
        end

        if occursin("function", line)
            
            case_name = extract_matpower_case_name(line)

            a_fun_found = true
            
            index += 1
        end

        
    end

    return  case_name

end

"""
`parse_matpower_data_header`

Returns `matrix_name`, `columns_header`, `line_count`
for a section in matpower file.
"""
function parse_matpower_data_header(
    lines,
    index;
    start_char = '[',
    end_char= ']' )
    
    last_index = length(lines)
    
    line_count = 0
    
    columns = -1

    @assert(occursin("=", lines[index+line_count]))
    
    matrix_assignment = split(
        lines[index+line_count], '%')[1]
    
    matrix_assignment = strip( matrix_assignment )

    header_line = index - 1
    
    columns_header =
        String.(split(split(
            lines[header_line],"%")[2],
                      r"\s+")[2:end])
    
    
    # @assert(occursin(".",matrix_assignment))
    
    matrix_assignment_parts =
        split(matrix_assignment, '=')
    
    matrix_name = strip(matrix_assignment_parts[1])
    
    line_count = line_count + 1
    
    return matrix_name, columns_header, line_count

end

"""
Test drives `parse_matpower_data_header`
"""
function driver_parse_matpower_data_header(
    data_string::String)
    
    data_lines = split(data_string, '\n')

    matlab_dict = Dict{String,Any}()
    
    case_name = nothing
    
    last_index = length(data_lines)
    
    index = 1
    
    while index <= last_index
        
        line = strip(data_lines[index])
        
        line = "$(line)"

        if length(line) <= 0 || strip(line)[1] == '%'
            
            index = index + 1
            
            continue
        end

        if  occursin("=",line) && occursin("[", line)
                
            matrix_name, columns_header, count_head =
                parse_matpower_data_header(
                    data_lines, index)
                
            matpower_dict["matrix_name"] = matrix_name

            matpower_dict["columns_header"] =
                columns_header
            
            index = index + count_head - 1
             
        end

        index += 1

    end

    return matpower_dict

end


"""
`extract_matpower_scalar`

Breaks up matlab strings of the form 'name = value;'

"""
function extract_matpower_scalar(
    string::AbstractString)

    name_value_pair = split(split(string, ';')[1], '=')
    
    @assert(length(name_value_pair) == 2)
    
    return (strip( replace(name_value_pair[1],
                           "." => "_" )),
            strip( replace(name_value_pair[2],
                           "'" => "" )))

end



"""
Test drives `extract_matpower_scalar`
"""
function driver_extract_matpower_scalar(
    data_string::String)
    
    data_lines = split(data_string, '\n')
    
    case_name = nothing
    
    last_index = length(data_lines)
    
    index = 1

    a_fun_found = false

    list_scalar_name = String[]

    list_scalar_value = []
    
    while index <= last_index 
        
        line = strip(data_lines[index])
        
        line = "$(line)"

        if length(line) <= 0 || strip(line)[1] == '%'
            
            index = index + 1
            
            continue
        end

        if occursin("=",line) && !occursin("{",line) && !occursin("[",line)

            name, value = extract_matpower_scalar(line)
            
            push!( list_scalar_name, name )

            push!( list_scalar_value, value )
                                    
            index += 1
        end
        
        index += 1
    end

    return list_scalar_name, list_scalar_value

end


"""
`parse_matpower_data_as_csv`

Returns a dictionary

 `Dict(
       "name" => matrix_name,
       "columns_header" => columns_header,
       "matrix_body" => matrix_body_lines,
       "line_count" => line_count )`

for a matpower file section.
"""
function parse_matpower_data_as_csv(
    lines,
    index;
    start_char = '[',
    end_char= ']' )
    
    last_index = length(lines)
    
    line_count = 0
    
    columns = -1

    @assert(occursin("=", lines[index+line_count]))
    
    matrix_assignment = split(
        lines[index+line_count], '%')[1]
    
    matrix_assignment = strip( matrix_assignment )

    header_line = index - 1
    
    columns_header =
        String.(split(split(
            lines[header_line],"%")[2],
                      r"\s+")[2:end])
    
    matrix_assignment_parts = split(
        matrix_assignment, '=')
    
    matrix_name = strip(matrix_assignment_parts[1])
    
    matrix_assignment_rhs = ""

    if length(matrix_assignment_parts) > 1

        matrix_assignment_rhs =
            matrix_assignment_parts[2]
    end
        
    line_count = line_count + 1
        
    matrix_body_lines = Vector{String}[  ]
    
    found_close_bracket =
        occursin(string(end_char),
                 matrix_assignment_parts[2] )

    while index + line_count < last_index && !found_close_bracket
        
        line =
            strip(lines[index + line_count] )

        if length(line) == 0 || line[1] == '%'
            
            line_count += 1
            
            continue
        end
        
        line =
            strip(split(line,'%')[1] )
        
        if occursin(string(end_char), line )
            
            found_close_bracket = true
            
        end

        line = strip(split(line, ';')[1])

        line = String.(split(line, r"\s+")[1:end])
        
        push!(matrix_body_lines, line)

        line_count = line_count + 1

    end
    
    return  Dict(
        "name" => matrix_name,
        "columns_header" => columns_header,
        "matrix_body" => matrix_body_lines,
        "line_count" => line_count )

end


"""
Test drives `parse_matpower_data_as_csv`
"""
function driver_parse_matpower_data_as_csv(
    data_string::String)
    
    data_lines = split(data_string, '\n')

    matpower_dict = Dict{String,Any}()
    
    case_name = nothing
    
    last_index = length(data_lines)
    
    index = 1
    
    while index <= last_index
        
        line = strip(data_lines[index])
        
        line = "$(line)"

        if length(line) <= 0 || strip(line)[1] == '%'
            
            index = index + 1
            
            continue
        end

        if  occursin("=",line) && occursin("[", line)
                
            matrix_dict =
                parse_matpower_data_as_csv(
                    data_lines, index)
                
            matpower_dict[matrix_dict["name"]] =
                matrix_dict

            index = index +
                matrix_dict["line_count"] - 1
        end

        index += 1
        
    end

    dict_keys = collect(keys(matpower_dict))

    for a_key in dict_keys

        matpower_dict[a_key]["matrix_body"] =
            matpower_dict[a_key][
                "matrix_body"][1:end-1]
    end

    return matpower_dict

end

"""
`save_matpower_scalar_as_csv`

Saves a matpower scalar section such as `mpc.baseMVA = 100` as a csv
"""
function save_matpower_scalar_as_csv(
    list_scalar_name,
    list_scalar_value,
    case_dir;
    csv_file = nothing )

    mpc_columns_header = list_scalar_name

    csv_body = join(list_scalar_value, ' ')


    data_in_df =
        CSV.File(IOBuffer( csv_body );
                 header = mpc_columns_header,
                          delim=' ' ) |> DataFrame

    if csv_file == nothing
        
        data_file_csv = "mpc_scalar.csv"
        
    else
        
        data_file_csv = csv_file

    end
    
    data_file = joinpath(case_dir, data_file_csv)

    CSV.write(data_file, data_in_df )
    
end


"""
`save_matpower_data_as_csv`

Saves a matpower section as a csv file. The inputs to this function
is a dictionary which is an output of `parse_matpower_data_as_csv`

 `Dict(
       "name" => matrix_name,
       "columns_header" => columns_header,
       "matrix_body" => matrix_body_lines,
       "line_count" => line_count )`

and a directory to store the csv file. The name of the csv file is
taken from `matrix_name` in the dictionary.

"""
function save_matpower_data_as_csv(
    matpower_dict, case_dir )
    
    data_keys = collect(keys(matpower_dict ))

    for a_key in data_keys

        a_mpc_data = matpower_dict[a_key][
            "matrix_body"]

        mpc_columns_header =
            matpower_dict[a_key][
                "columns_header"]

        csv_matrix_body =
            [ join(an_item, ' ')
              for an_item in a_mpc_data ]

       csv_body = join(csv_matrix_body, '\n' )

        data_in_df =
            CSV.File(IOBuffer( csv_body );
                     header = mpc_columns_header,
                              delim=' ' ) |> DataFrame

        a_key_as_name = replace(a_key, "." => "_")
        
        data_file_csv = "$(a_key_as_name).csv"
        
        data_file = joinpath(
            case_dir, data_file_csv)

        CSV.write( data_file, data_in_df )
    end

end


"""
`export_maptpower_data_to_csv`

Exports contents of sections of a matpower file as csv files
to a specified directory.
"""
function export_maptpower_data_to_csv(
    matpower_file,
    csv_output_dir )

    data_string = read(open( matpower_file),String)

    matpower_dict =
        driver_parse_matpower_data_as_csv(data_string)

    save_matpower_data_as_csv(
        matpower_dict, csv_output_dir )

    list_scalar_name, list_scalar_value =
        driver_extract_matpower_scalar(data_string)

    save_matpower_scalar_as_csv(
        list_scalar_name,
        list_scalar_value,
        csv_output_dir )

    
end


"""
`get_matpower_cases_files`

Returns all matpower files in a folder
"""
function get_matpower_cases_files(
    matpower_data_folder )

    return glob(
        "case*.m", matpower_data_folder ) 

end


"""
`created_csv_cases_data_folders`

Creates csv files for matpower cases in associated cases folders.
"""
function created_csv_cases_data_folders(
    matpower_cases_files,
    base_csv_folder  )

    cases_names =
        [ split(splitpath(
            a_case_file)[end], ".")[1]          
          for a_case_file in matpower_cases_files  ]
    
    for a_case in cases_names
        case_dir =
            joinpath(base_csv_folder, a_case)
        
        mkpath(case_dir)

    end
    
end


"""
`export_maptpower_cases_to_csv`

Exports all matpower case files in a given directory to csv
files in associated cases folders.

"""
function export_maptpower_cases_to_csv(
    matpower_data_folder,
    base_csv_folder )

    matpower_cases_files =
        get_matpower_cases_files(
            matpower_data_folder )

    cases_names =
        [ split(splitpath( a_case_file)[end], ".")[1]
          
          for a_case_file in matpower_cases_files  ]

    for (a_matpower_case_file, a_case) in
        zip( matpower_cases_files, cases_names )

        csv_case_dir =
            joinpath( base_csv_folder, a_case )

        mkpath( csv_case_dir )

         export_maptpower_data_to_csv(
             a_matpower_case_file,
             csv_case_dir)
        
    end

end


"""
`get_matpower_scalar_as_iobuffer_by_symbol`

An utility function to get matpower scalar based on a given symbol
`type_key_sym`.
"""
function get_matpower_scalar_as_iobuffer_by_symbol(
    case_file; type_key_sym::Symbol = :mpc_baseMVA )
    
    data_string = read( open( case_file ), String)
    
    list_scalar_name, list_scalar_value =
        driver_extract_matpower_scalar( data_string )
    
    mpc_columns_header = list_scalar_name

    csv_body = join( list_scalar_value, ' ' )

    # data_in_CSV_File =

    return  getproperty(CSV.File(IOBuffer( csv_body );
                 header = mpc_columns_header,
                          delim=' ' ), type_key_sym) 
end

"""
`get_matpower_scalar_as_iobuffer_by_case_file`

An utility function to get matpower scalar based on a given string
`type_key_string`.
"""
function get_matpower_scalar_as_iobuffer_by_case_file(
    case_file; type_key_string::String = "mpc_baseMVA" )

    type_key_sym = Symbol( type_key_string )

    return  get_matpower_scalar_as_iobuffer_by_symbol(
    case_file; type_key_sym = type_key_sym ) 
end

"""
`get_matpower_mpc_type_iobuffer_by_dict`

An utility function to get matpower mpc type (section) based on
a given string `type_key`.
"""
function get_matpower_mpc_type_iobuffer_by_dict(
    matpower_dict; type_key="mpc.gen" )

    mpc_type =
        matpower_dict[type_key ]["matrix_body"] 

    mpc_type_columns_header =
        matpower_dict[type_key ]["columns_header"]


    csv_matrix_body =
        [ join(an_item, ' ')
          for an_item in mpc_type ]

    csv_body = join(csv_matrix_body, '\n' )

    # data_propertynames =
    #     propertynames(data_in_CSV_File)

    # data_in_CSV_File
    
    return  CSV.File(IOBuffer( csv_body );
                 header =
                     mpc_type_columns_header,
                 delim=' ' )


end

"""
`get_matpower_mpc_type_iobuffer_by_case_file`

An utility function to get matpower mpc type (section) based on
a given string `type_key`.
"""
function get_matpower_mpc_type_iobuffer_by_case_file(
    case_file; type_key="mpc.gen" )

    data_string = read(open(case_file),String)

    matpower_dict =
        driver_parse_matpower_data_as_csv(
            data_string )

    return get_matpower_mpc_type_iobuffer_by_dict(
        matpower_dict; type_key=type_key )

end


function get_matpower_dict_by_case_file(
    case_file )

    data_string = read(open(case_file),String)

    # matpower_dict
    
    return driver_parse_matpower_data_as_csv(
        data_string )

end






