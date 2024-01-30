# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Coordinates{N,T}

Parent type of all coordinate types.
"""
abstract type Coordinates{N,T} end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("coordinates/basic.jl")
include("coordinates/gis.jl")
