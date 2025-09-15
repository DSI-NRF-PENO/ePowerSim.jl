# (c) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# ------------------------------------------------------
#  
# ------------------------------------------------------


#################################################
# -----------------------------------------------
# export header file Parser 
# -----------------------------------------------
#################################################

"""
`extract_function_part_in_export_file`

Returns the function part for a line in
export header file.

"""
function extract_function_part_in_export_file(
    string::AbstractString)

    if occursin("export", string)

        return  strip(split(string, "export")[2] )       

    else
        return nothing

    end
    
end


"""
Test drives `extract_function_part_in_export_file`
"""
function driver_extract_function_part_in_export_file(
    data_string)

    list_exported_funs =
        String[a_string
               for a_string in data_string
                   if occursin("export", a_string) &&  !occursin('#', a_string) ]

    return  (;data_string, list_exported_funs)

end

function extract_multiple_export_funcs(
    export_header_file )

    data_string = open(readlines, export_header_file)
        
    org_data_string, list_exported_funs =

        driver_extract_function_part_in_export_file(
            data_string)

    exported_funs_countmap =
        countmap(list_exported_funs)

    multi_exported_funs =
        sort([(v, k) for (k, v) in
                  exported_funs_countmap if v > 1],
             by = x -> x[1])

    # return (;org_data_string,
    #         list_exported_funs,
    #         multi_exported_funs)


    return (;
            multi_exported_funs, x)
    
    
end

 
# """

# package_dir = pkgdir(ePowerSim)

# src_dir =
#     joinpath(package_dir, "src")

# export_header_file =
#     joinpath(src_dir,
#              "headers",
#              "sd-dynamics-export.jl" )

# multiple_export_funcs =
#     extract_multiple_export_funcs(
#         export_header_file )


# """
x
