# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimplificationMethod

A method for simplifying geometries.
"""
abstract type SimplificationMethod end

"""
    simplify(geometry, method)

Simplify `geometry` with given `method`.
"""
function simplify end

function simplify(domain::Domain, method::SimplificationMethod)
  GeometrySet([simplify(elem, method) for elem in domain])
end

function simplify(multi::Multi, method::SimplificationMethod)
  Multi([simplify(geom, method) for geom in multi])
end

function simplify(polygon::Polygon, method::SimplificationMethod)
  c = [simplify(chain, method) for chain in chains(polygon)]
  PolyArea(c[1], c[2:end])
end

"""
    decimate(geometry, ϵ)

Simplify `geometry` with an appropriate
simplification method and tolerance `ϵ`.
"""
decimate(geometry, ϵ) = simplify(geometry, DouglasPeucker(ϵ))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("simplification/douglaspeucker.jl")
