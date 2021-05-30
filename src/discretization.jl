# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DiscretizationMethod

A method for discretizing geometries into meshes.
"""
abstract type DiscretizationMethod end

"""
    discretize(geometry, method)

Discretize `geometry` with discretization `method`.
"""
function discretize end

discretize(multi::Multi, method::DiscretizationMethod) =
  mapreduce(geometry -> discretize(geometry, method), merge, multi)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("discretization/fist.jl")
include("discretization/dehn.jl")
