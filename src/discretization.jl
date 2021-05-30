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

function discretize(polygon::Polygon, method::DiscretizationMethod)
  # build bridges in case the polygon has holes,
  # i.e. reduce to a single outer boundary
  𝒫 = polygon |> unique |> bridge

  # discretize using outer boundary
  discretize(𝒫, method)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("discretization/fist.jl")
include("discretization/dehn.jl")
