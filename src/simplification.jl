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

Simplify geometric `object` with given `method`.
"""
function simplify end

simplify(box::Box{𝔼{2}}, method::SimplificationMethod) = Pgon(simplify(boundary(box), method))

simplify(polygon::Polygon, method::SimplificationMethod) = Pgon([simplify(ring, method) for ring in rings(polygon)])

simplify(multi::Multi, method::SimplificationMethod) = Multi([simplify(geom, method) for geom in parent(multi)])

simplify(domain::Domain, method::SimplificationMethod) = GeometrySet([simplify(elem, method) for elem in domain])

# ----------------
# IMPLEMENTATIONS
# ----------------

include("simplification/selinger.jl")
include("simplification/douglaspeucker.jl")
include("simplification/minmax.jl")
