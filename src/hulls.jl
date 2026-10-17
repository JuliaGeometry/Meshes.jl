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
include("hulls/moreira.jl")

# ----------
# UTILITIES
# ----------

"""
    convexhull(object)

Convex hull of `object`.
"""
convexhull(object) = _hullwith(object, GrahamScan())

"""
    concavehull(object)

Concave hull of `object`.
"""
concavehull(object) = _hullwith(object, MoreiraMarch())

# ----------
# FALLBACKS
# ----------

_hullwith(object, method) = hull(_hullpoints(object), method)

# ----------------
# SPECIALIZATIONS
# ----------------

_hullwith(p::Point, method) = p

_hullwith(b::Box, method) = b

_hullwith(b::Ball, method) = b

_hullwith(s::Sphere, method) = Ball(center(s), radius(s))

_hullwith(t::Triangle, method) = t

_hullwith(g::Grid, method) = Box(extrema(g)...)

# ----------------
# IMPLEMENTATIONS
# ----------------

_hullpoints(p::Polytope) = eachvertex(p)

_hullpoints(m::Mesh) = eachvertex(m)

_hullpoints(p::Primitive) = _hullpoints(boundary(p))

_hullpoints(m::Multi) = _hullpoints(parent(m))

_hullpoints(geoms) = (p for g in geoms for p in boundarypoints(g))
