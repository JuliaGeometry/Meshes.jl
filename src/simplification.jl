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

See also [`decimate`](@ref).
"""
function simplify end

simplify(box::Box, method::SimplificationMethod) = _simplify(box, Val(embeddim(box)), method)

_simplify(box::Box, ::Val{2}, method::SimplificationMethod) = PolyArea(simplify(boundary(box), method))

simplify(polygon::Polygon, method::SimplificationMethod) = PolyArea([simplify(ring, method) for ring in rings(polygon)])

simplify(multi::Multi, method::SimplificationMethod) = Multi([simplify(geom, method) for geom in parent(multi)])

simplify(domain::Domain, method::SimplificationMethod) = GeometrySet([simplify(elem, method) for elem in domain])

# ----------------
# IMPLEMENTATIONS
# ----------------

include("simplification/douglaspeucker.jl")
include("simplification/selinger.jl")

# ----------
# UTILITIES
# ----------

"""
    decimate(object, [ϵ]; min=3, max=typemax(Int), maxiter=10)

Simplify `object` with an appropriate simplification method
and deviation tolerance `ϵ`.

If the tolerance `ϵ` is not provided, perform binary search until
the number of vertices is between `min` and `max` or until the
number of iterations reaches a maximum `maxiter`.
"""
decimate(object, ϵ=nothing; min=3, max=typemax(Int), maxiter=10) =
  simplify(object, DouglasPeucker(ϵ, min=min, max=max, maxiter=maxiter))
