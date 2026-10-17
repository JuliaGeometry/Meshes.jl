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

Convex hull of `object` with [`GrahamScan`](@ref) method.
"""
convexhull(object) = _hullwith(object, GrahamScan())

"""
    concavehull(object)

Concave hull of `object` with [`MoreiraMarch`](@ref) method.
"""
concavehull(object) = _hullwith(object, MoreiraMarch())

# ----------
# FALLBACKS
# ----------

_hullwith(object, method) = hull(_hullpoints(object), method)

_hullwith(p::Primitive, method) = _hullwith(boundary(p), method)

# ----------------
# SPECIALIZATIONS
# ----------------

_hullwith(p::Point, ::HullMethod) = p

_hullwith(b::Box, ::HullMethod) = b

_hullwith(b::Ball, ::HullMethod) = b

_hullwith(s::Sphere, ::HullMethod) = Ball(center(s), radius(s))

_hullwith(t::Triangle, ::HullMethod) = t

_hullwith(g::Grid, ::HullMethod) = Box(extrema(g)...)

# -----------------
# HELPER FUNCTIONS
# -----------------

_hullpoints(p::Polytope) = eachvertex(p)

_hullpoints(m::Mesh) = eachvertex(m)

_hullpoints(m::Multi) = _hullpoints(parent(m))

_hullpoints(geoms) = (p for g in geoms for p in boundarypoints(g))
