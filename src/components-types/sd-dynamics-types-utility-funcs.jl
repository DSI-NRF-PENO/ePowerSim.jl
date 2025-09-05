
########################################################
# ------------------------------------------------------
#  System Dynamics Data Type Utility fubctions
# ------------------------------------------------------
########################################################


"""
```Julia
get_subsubtype(AbstractPowerSystemComponent)
```
It returns the subtypes of subtypes of abstract types of `AbstractPowerSystemComponent`. These are the subtypes of  `SdGen`, `SdNonGen`, `SdGov`, `SdAvr`,  `SdPss`, `SdNonGenPlant`, `SdGenPlant` and  `SdBranchElement`. It is used to get `list_components_types` which is a list of components defined in this package.

"""
function get_subsubtype(AbstractPowerSystemComponent )
    list_two_depth_subtypes =
        [Vector{DataType}[
            subtypes(a_type)
            for a_type in
                subtypes(AbstractPowerSystemComponent)]...;]
    return Tuple(list_two_depth_subtypes)
end


"""
```ComponentsSymbols``` is a struct used to pack all types of components symbols stored in dictionary of symbols dict_syms.

The struct is designed to make a clean implementation of ```symbolof(name)```.

It is a callable struct.
"""
struct ComponentsSymbols
    dict_syms::Dict{Symbol, Dict{String, Vector{Symbol}}}
end


# function (compsyms::ComponentsSymbols)(
#     name::Union{get_subsubtype(AbstractPowerSystemComponent)...})
#     return compsyms.dict_syms[Symbol("$(typeof(name))")]["all_syms"]
# end


# function (compsyms::ComponentsSymbols)(struct_name::Union{Symbol, String})
#     return compsyms.dict_syms[Symbol(struct_name)]["all_syms"]
# end


# function (compsyms::ComponentsSymbols)(name)
#     return compsyms.dict_syms[Symbol("$(name)")]["all_syms"]
# end


# """
# ``` symbolsof =  ComponentsSymbols(dict_syms)``` is a callable struct that can be used insted of

# list_components_types = get_subsubtype(AbstractComponentType)

# symbolsof(name) = dict_syms[Symbol("$(name)")]["all_syms"]

# symbolsof(name::Union{list_components_types...}) = dict_syms[Symbol("$(typeof(name))")]["all_syms"]

# symbolsof(struct_name::Union{Symbol, String}) = dict_syms[Symbol(struct_name)]["all_syms"]

# """
# # symbolsof =  ComponentsSymbols(dict_syms)

# # ------------------------------------------------------
