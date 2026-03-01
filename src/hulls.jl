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
include("hulls/concave.jl")

# ----------
# UTILITIES
# ----------

"""
    convexhull(object)

Convex hull of `object`.
"""
function convexhull end

"""
    concavehull(object)

Concave hull of `object`.
"""
function concavehull end

# ----------
# FALLBACKS
# ----------

convexhull(p::Polytope) = _pconvexhull(eachvertex(p))

convexhull(p::Primitive) = convexhull(boundary(p))

convexhull(m::Multi) = _gconvexhull(parent(m))

convexhull(geoms) = _gconvexhull(geoms)

concavehull(p::Polytope) = _pconcavehull(eachvertex(p))

concavehull(p::Primitive) = concavehull(boundary(p))

concavehull(m::Multi) = _gconcavehull(parent(m))

concavehull(geoms) = _gconcavehull(geoms)

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

concavehull(p::Point) = p

concavehull(b::Box) = b

concavehull(b::Ball) = b

concavehull(s::Sphere) = Ball(center(s), radius(s))

concavehull(t::Triangle) = t

concavehull(g::Grid) = Box(extrema(g)...)

concavehull(m::Mesh) = _pconcavehull(eachvertex(m))

# ----------------
# IMPLEMENTATIONS
# ----------------

_gconvexhull(geoms) = _pconvexhull(p for g in geoms for p in boundarypoints(g))

_pconvexhull(points) = hull(points, GrahamScan())

_gconcavehull(geoms) = _pconcavehull(parent(geoms))

_pconcavehull(points) = hull(points, Concave())
