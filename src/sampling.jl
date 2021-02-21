# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SamplingMethod

A method for sampling from geometries.
"""
abstract type SamplingMethod end

"""
    sample(geometry, method)

Sample elements from `geometry` with `method`.
"""
function sample end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("sampling/regular.jl")
include("sampling/uniform.jl")
