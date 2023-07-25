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
    PointwiseGeometricTransform

A pointwise [`GeometricTransform`](@ref) defined in terms of a
single point. Pointwise transforms can be easily applied to
all types of geometries and meshes using fallback methods.
"""
abstract type PointwiseGeometricTransform <: GeometricTransform end

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

# --------------------
# TRANSFORM FALLBACKS
# --------------------

function apply(transform::GeometricTransform, object)
  prep = preprocess(transform, object)
  newpoints, pcache = applypoint(transform, pointify(object), prep)
  _reconstruct(newpoints, object), pcache
end

function revert(transform::GeometricTransform, newobject, cache)
  points = revertpoint(transform, pointify(newobject), cache)
  _reconstruct(points, newobject)
end

function reapply(transform::GeometricTransform, object, cache)
  newpoints = reapplypoint(transform, pointify(object), cache)
  _reconstruct(newpoints, object)
end

# convert lists of points into objects
_reconstruct(points, ::Point) = first(points)
_reconstruct(points, ::G) where {G<:Geometry} = G(points)
_reconstruct(points, domain::Domain) = _reconstruct(points, GeometrySet(domain))
_reconstruct(points, mesh::Mesh) = SimpleMesh(points, topology(mesh))
_reconstruct(points, ::PointSet) = PointSet(points)
_reconstruct(points, ::PolyArea) = PolyArea(points)
_reconstruct(points, ::Ring) = Ring(points)
_reconstruct(points, ::Rope) = Rope(points)
_reconstruct(points, ::PL) where {PL<:Polytope} = PL(ntuple(i -> @inbounds(points[i]), nvertices(PL)))

# --------------------
# STATELESS FALLBACKS
# --------------------

reapply(transform::StatelessGeometricTransform, object, cache) = apply(transform, object) |> first

# --------------------
# POINTWISE FALLBACKS
# --------------------

apply(t::PointwiseGeometricTransform, g) = _apply(t, g), nothing

revert(t::PointwiseGeometricTransform, g, cache) = _apply(inv(t), g)

# apply transform recursively
_apply(t::PointwiseGeometricTransform, g::G) where {G<:Geometry} =
  G((_apply(t, getfield(g, n)) for n in fieldnames(G))...)
  
# fallbacks for lists of geometries
_apply(t::PointwiseGeometricTransform, p::NTuple) = ntuple(i -> _apply(t, p[i]), length(p))
_apply(t::PointwiseGeometricTransform, g::AbstractVector{<:Geometry}) = map(gᵢ -> _apply(t, gᵢ), g)
_apply(t::PointwiseGeometricTransform, p::CircularVector{<:Point}) = map(pᵢ -> _apply(t, pᵢ), p)

# stop recursion at specific types
_apply(t::PointwiseGeometricTransform, o::Number) = o

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transforms/rotate.jl")
include("transforms/translate.jl")
include("transforms/stretch.jl")
include("transforms/stdcoords.jl")
include("transforms/repair.jl")
include("transforms/smoothing.jl")
