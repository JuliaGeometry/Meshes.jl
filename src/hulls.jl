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
    hull(object, method)

Compute the hull of the `object` with given `method`.
"""
function hull end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("hulls/graham.jl")
include("hulls/jarvis.jl")

# ----------
# UTILITIES
# ----------

"""
    convexhull(object)

Compute the convex hull of the `object` with an appropriate method.
"""
function convexhull end

convexhull(p::Point) = p

convexhull(b::Box) = b

convexhull(b::Ball) = b

convexhull(s::Sphere) = Ball(center(s), radius(s))

convexhull(t::Triangle) = t

convexhull(p::Primitive) = convexhull(boundary(p))

convexhull(p::Polytope) = convexhull(vertices(p))

convexhull(m::Multi) = convexhull(collect(m))

convexhull(g::Geometry) = convexhull([g])

convexhull(g::Grid) = boundingbox(g)

convexhull(d::Domain) = convexhull(collect(d))

convexhull(d::Data) = convexhull(domain(d))

convexhull(p::AbstractVector{<:Point{2}}) = hull(p, GrahamScan())

convexhull(g::AbstractVector{<:Geometry{2}}) = mapreduce(pointify, vcat, g) |> unique |> convexhull
