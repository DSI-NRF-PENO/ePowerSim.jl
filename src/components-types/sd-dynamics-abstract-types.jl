
########################################################
# ------------------------------------------------------
#  System Dynamics Data Type
# ------------------------------------------------------
########################################################


"""
```Julia
abstract type SdGen end

abstract type SdNonGen end

abstract type SdGov end

abstract type SdAvr end

abstract type SdPss end

abstract type SdNonGenPlant end

abstract type SdGenPlant end

# Edges types

abstract type SdBranchElement end

```

Definition of abstract types of structures
"""

abstract type AbstractPowerSystemComponent end

# Nodes types

abstract type SdGen <: AbstractPowerSystemComponent end

abstract type SdNonGen <: AbstractPowerSystemComponent end

abstract type SdGov <: AbstractPowerSystemComponent end

abstract type SdAvr <: AbstractPowerSystemComponent end

abstract type SdPss <: AbstractPowerSystemComponent end

abstract type SdNonGenPlant <: AbstractPowerSystemComponent end

abstract type SdGenPlant <: AbstractPowerSystemComponent end

# Edges types

abstract type SdBranchElement <: AbstractPowerSystemComponent end

####################################################
