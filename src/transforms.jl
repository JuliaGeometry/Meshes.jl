# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometricTransform <: Transform

A method to transform the geometry (e.g. coordinates) of objects.
See https://en.wikipedia.org/wiki/Geometric_transformation.
"""
abstract type GeometricTransform <: TAPI.Transform end

"""
    StatelessGeometricTransform <: StatelessTransform

A stateless [`GeometricTransform`](@ref) as defined in TransformsAPI.jl.
"""
abstract type StatelessGeometricTransform <: TAPI.StatelessTransform end

# helper type for fallback definitions
const StatefulOrStateless = Union{GeometricTransform,StatelessGeometricTransform}

# --------------------
# TRANSFORM FALLBACKS
# --------------------

# convert objects into lists of points
_points(p::Point)      = [p]
_points(g::Geometry)   = vertices(g)
_points(d::Domain)     = mapreduce(_points, vcat, d)
_points(m::Mesh)       = vertices(m)
_points(p::PointSet)   = collect(p)

# convert lists of points into objects
_reconstruct(points, ::Point) = first(points)
_reconstruct(points, ::G) where {G<:Geometry} = G(points)

function TAPI.apply(transform::StatefulOrStateless, object)
  newpoints, pcache = applypoints(transform, _points(object))
  _reconstruct(newpoints, object), pcache
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transforms/taubin.jl")