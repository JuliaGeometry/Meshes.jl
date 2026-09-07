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
convexhull(object) = _hull(object, GrahamScan())

"""
    concavehull(object)

Concave hull of `object`.
"""
concavehull(object) = _hull(object, MoreiraMarch())

# ----------
# FALLBACKS
# ----------

_hull(object, method) = hull(_hullpoints(object), method)

# ----------------
# SPECIALIZATIONS
# ----------------

_hull(g::Union{Point,Box,Ball,Triangle}, method) = g

_hull(s::Sphere, method) = Ball(center(s), radius(s))

_hull(g::Grid, method) = Box(extrema(g)...)

# ----------------
# IMPLEMENTATIONS
# ----------------

_hullpoints(p::Polytope) = eachvertex(p)

_hullpoints(m::Mesh) = eachvertex(m)

_hullpoints(p::Primitive) = _hullpoints(boundary(p))

_hullpoints(m::Multi) = _hullpoints(parent(m))

_hullpoints(geoms) = (p for g in geoms for p in boundarypoints(g))
