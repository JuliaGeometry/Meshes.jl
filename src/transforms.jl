# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometricTransform

A method to transform the geometry (e.g. coordinates) of objects.
See [https://en.wikipedia.org/wiki/Geometric_transformation]
(https://en.wikipedia.org/wiki/Geometric_transformation).
"""
abstract type GeometricTransform <: Transform end

"""
    StatelessGeometricTransform

A stateless [`GeometricTransform`](@ref) as defined in
[TransformsAPI.jl](https://github.com/JuliaML/TransformsAPI.jl).
"""
abstract type StatelessGeometricTransform <: GeometricTransform end

"""
    newpoints, pcache = applypoint(transform, points, prep)

Implementation of [`apply`](@ref) for points of object.
This function is intended for developers of new transform types.
"""
function applypoint end

"""
    points = revertpoint(transform, newpoints, pcache)

Implementation of [`revert`](@ref) for points of object.
This function is intended for developers of new transform types.
"""
function revertpoint end

"""
    newpoints = reapplypoint(transform, points, pcache)

Implementation of [`reapply`](@ref) for points of object.
This function is intended for developers of new transform types.
"""
function reapplypoint end

# -----------------
# HELPER FUNCTIONS
# -----------------

# convert objects into lists of points
_points(p::Point)      = [p]
_points(g::Geometry)   = vertices(g)
_points(d::Domain)     = mapreduce(_points, vcat, d)
_points(m::Mesh)       = vertices(m)
_points(p::PointSet)   = collect(p)

# convert lists of points into objects
_reconstruct(points, ::Point) = first(points)
_reconstruct(points, ::G) where {G<:Geometry} = G(points)
_reconstruct(points, domain::Domain) =
  _reconstruct(points, Collection(collect(domain)))
_reconstruct(points, mesh::Mesh) = SimpleMesh(points, topology(mesh))
_reconstruct(points, ::PointSet) = PointSet(points)

# --------------------
# TRANSFORM FALLBACKS
# --------------------

function apply(transform::GeometricTransform, object)
  prep = preprocess(transform, object)
  newpoints, pcache = applypoint(transform, _points(object), prep)
  _reconstruct(newpoints, object), pcache
end

function revert(transform::GeometricTransform, newobject, cache)
  points = revertpoint(transform, _points(newobject), cache)
  _reconstruct(points, newobject)
end

function reapply(transform::GeometricTransform, object, cache)
  newpoints = reapplypoint(transform, _points(object), cache)
  _reconstruct(newpoints, object)
end

# --------------------
# STATELESS FALLBACKS
# --------------------

reapply(transform::StatelessGeometricTransform, object, cache) =
  apply(transform, object) |> first

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transforms/stdcoords.jl")
include("transforms/scalecoords.jl")
include("transforms/rotatecoords.jl")
include("transforms/translatecoords.jl")
include("transforms/smoothing.jl")
