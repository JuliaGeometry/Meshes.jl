# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometricTransform

A method to transform the geometry (e.g. coordinates) of objects.
See [https://en.wikipedia.org/wiki/Geometric_transformation]
(https://en.wikipedia.org/wiki/Geometric_transformation).
"""
abstract type GeometricTransform <: TAPI.Transform end

"""
    StatelessGeometricTransform

A stateless [`GeometricTransform`](@ref) as defined in
[TransformsAPI.jl](https://github.com/JuliaML/TransformsAPI.jl).
"""
abstract type StatelessGeometricTransform <: TAPI.StatelessTransform end

"""
    newpoints, pcache = applypoint(transform, points, prep)

Implementation of [`apply`](@ref) for points of object.
This function is intended for developers of new transforms.
"""
function applypoint end

"""
    points = revertpoint(transform, newpoints, pcache)

Implementation of [`revert`](@ref) for points of object.
This function is intended for developers of new transforms.
"""
function revertpoint end

# --------------------
# TRANSFORM FALLBACKS
# --------------------

# helper type for fallback definitions
const StatefulOrStateless = Union{GeometricTransform,StatelessGeometricTransform}

# convert objects into lists of points
_points(p::Point)      = [p]
_points(g::Geometry)   = vertices(g)
_points(d::Domain)     = mapreduce(_points, vcat, d)
_points(m::Mesh)       = vertices(m)
_points(p::PointSet)   = collect(p)

# convert lists of points into objects
_reconstruct(points, ::Point) = first(points)
_reconstruct(points, ::G) where {G<:Geometry} = G(points)
_reconstruct(points, mesh::Mesh) = SimpleMesh(points, topology(mesh))

function TAPI.apply(transform::StatefulOrStateless, object)
  prep = TAPI.preprocess(transform, object)
  newpoints, pcache = applypoint(transform, _points(object), prep)
  _reconstruct(newpoints, object), pcache
end

function TAPI.revert(transform::StatefulOrStateless, newobject, cache)
  points = revertpoint(transform, _points(newobject), cache)
  _reconstruct(points, newobject)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transforms/taubin.jl")