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

convexhull(p::Polytope) = _pconvexhull(eachvertex(p))

convexhull(p::Primitive) = convexhull(boundary(p))

convexhull(m::Multi) = _gconvexhull(parent(m))

convexhull(geoms) = _gconvexhull(geoms)

# ----------------
# SPECIALIZATIONS
# ----------------

convexhull(p::Point) = p

convexhull(b::Box) = b

convexhull(b::Ball) = b

convexhull(s::Sphere) = Ball(center(s), radius(s))

convexhull(t::Triangle) = t

convexhull(g::Grid) = Box(extrema(g)...)

convexhull(m::Mesh) = _pconvexhull(eachvertex(m))

# ----------------
# IMPLEMENTATIONS
# ----------------

_gconvexhull(geoms) = _pconvexhull(p for g in geoms for p in pointify(g))

_pconvexhull(points) = hull(points, GrahamScan())
