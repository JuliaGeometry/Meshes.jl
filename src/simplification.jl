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

function simplify(box::Box{2}, method::SimplificationMethod)
  c = simplify(boundary(box), method)
  PolyArea(c)
end

function simplify(polygon::Polygon, method::SimplificationMethod)
  c = [simplify(chain, method) for chain in chains(polygon)]
  PolyArea(c[1], c[2:end])
end

function simplify(multi::Multi, method::SimplificationMethod)
  Multi([simplify(geom, method) for geom in collect(multi)])
end

function simplify(domain::Domain, method::SimplificationMethod)
  Collection([simplify(elem, method) for elem in domain])
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("simplification/douglaspeucker.jl")

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
