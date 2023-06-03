# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HullMethod

A method for computing hulls of point sets or
collections of geometries.
"""
abstract type HullMethod end

"""
    hull(object, [method])

Compute the hull of the `object` with given `method`.
"""
function hull end

# ----------
# FALLBACKS
# ----------

hull(g::Polytope, method) = hull(vertices(g), method)

hull(p::Primitive, method) = hull(boundary(p), method)

hull(m::Multi, method) = hull(collect(m), method)

hull(d::Domain, method) = hull(collect(d), method)

hull(d::Data, method) = hull(domain(d), method)

# ---------
# DEFAULTS
# ---------

hull(g::Geometry{2}) = hull(g, GrahamScan())

hull(d::Domain{2}) = hull(d, GrahamScan())

# ----------------
# SPECIALIZATIONS
# ----------------

hull(p::Point) = boundingbox(p)

hull(b::Box) = b

hull(b::Ball) = b

hull(s::Sphere) = Ball(center(s), radius(s))

hull(t::Triangle) = t

hull(g::Grid) = boundingbox(g)

# ----------------
# IMPLEMENTATIONS
# ----------------

function hull(geoms::AbstractVector{<:Geometry}, method)
  verts(geom) = vertices(discretize(boundary(geom)))
  hull(mapreduce(verts, vcat, geoms), method)
end

include("hulls/graham.jl")