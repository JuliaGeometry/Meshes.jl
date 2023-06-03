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

hull(d::Domain, method) = hull(collect(d), method)

hull(d::Data, method) = hull(domain(d), method)

hull(g::AbstractVector{<:Geometry}, method) = hull(Multi(g), method)

# ---------
# DEFAULTS
# ---------

hull(g::Geometry{2}) = hull(g, GrahamScan())

hull(d::Domain{2}) = hull(d, GrahamScan())

# ----------------
# SPECIALIZATIONS
# ----------------

hull(p::Point) = p

hull(b::Box) = b

hull(b::Ball) = b

hull(t::Triangle) = t

hull(g::Grid) = boundingbox(g)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("hulls/graham.jl")