# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimplificationMethod

A method for simplifying geometric objects.
"""
abstract type SimplificationMethod end

"""
    simplify(object, method)

Simplify `object` with given `method`.
"""
function simplify end

function simplify(box::Box{2}, method::SimplificationMethod)
  c = simplify(boundary(box), method)
  PolyArea(c)
end

function simplify(polygon::Polygon, method::SimplificationMethod)
  c = [simplify(chain, method) for chain in chains(polygon)]
  PolyArea(c[1], c[2:end])
end

function simplify(multi::Multi, method::SimplificationMethod)
  Multi([simplify(geom, method) for geom in multi])
end

function simplify(domain::Domain, method::SimplificationMethod)
  GeometrySet([simplify(elem, method) for elem in domain])
end

"""
    decimate(object, ϵ)

Simplify `object` with an appropriate
simplification method and tolerance `ϵ`.
"""
decimate(object, ϵ) = simplify(object, DouglasPeucker(ϵ))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("simplification/douglaspeucker.jl")
