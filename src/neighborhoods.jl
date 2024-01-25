# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Neighborhood

A neighborhood is a geometry that is not attached to any specific
point in space, and is free to slide over a domain of interest.
"""
abstract type Neighborhood end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("neighborhoods/metricball.jl")
