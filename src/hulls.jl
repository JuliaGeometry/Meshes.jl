# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HullMethod

A method for computing hulls of geometries.
"""
abstract type HullMethod end

"""
    hull(points, method)

Compute the hull of `points` with given `method`.
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

Convex hull of `object`.
"""
function convexhull end

# ----------
# FALLBACKS
# ----------

convexhull(p::Polytope) = convexhull(vertices(p))

convexhull(p::Primitive) = convexhull(boundary(p))

convexhull(m::Multi) = convexhull(collect(m))

convexhull(g::Geometry) = convexhull([g])

convexhull(d::Domain) = convexhull(collect(d))

# ----------------
# SPECIALIZATIONS
# ----------------

convexhull(p::Point) = Box(p, p)

convexhull(b::Box) = b

convexhull(b::Ball) = b

convexhull(s::Sphere) = Ball(center(s), radius(s))

convexhull(t::Triangle) = t

convexhull(g::Grid) = Box(extrema(g)...)

convexhull(m::Mesh) = convexhull(vertices(m))

# ----------------
# IMPLEMENTATIONS
# ----------------

convexhull(p::AbstractVector{<:Point{2}}) = hull(p, GrahamScan())

convexhull(g::AbstractVector{<:Geometry{2}}) = mapreduce(pointify, vcat, g) |> unique |> convexhull
