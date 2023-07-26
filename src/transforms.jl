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
    CoordinateTransform

A method to transform the coordinates of objects.
See [https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations]
(https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations).
"""
abstract type CoordinateTransform <: GeometricTransform end

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

function apply(t::GeometricTransform, o)
  prep = preprocess(t, o)
  p, c = applypoint(t, pointify(o), prep)
  _reconstruct(p, o), c
end

function revert(t::GeometricTransform, o, c)
  p = revertpoint(t, pointify(o), c)
  _reconstruct(p, o)
end

function reapply(t::GeometricTransform, o, c)
  p = reapplypoint(t, pointify(o), c)
  _reconstruct(p, o)
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
# POINTWISE FALLBACKS
# --------------------

apply(t::CoordinateTransform, g) = _apply(t, g), nothing

revert(t::CoordinateTransform, g, c) = _apply(inv(t), g)

# apply transform recursively
_apply(t::CoordinateTransform, g::G) where {G} = G((_apply(t, getfield(g, n)) for n in fieldnames(G))...)

# stop recursion at specific field types
_apply(::CoordinateTransform, x::Number) = x
_apply(::CoordinateTransform, x::Topology) = x

# fallbacks for lists of geometries
_apply(t::CoordinateTransform, g::NTuple{<:Any,<:Geometry}) = map(gᵢ -> _apply(t, gᵢ), g)
_apply(t::CoordinateTransform, g::AbstractVector{<:Geometry}) = map(gᵢ -> _apply(t, gᵢ), g)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transforms/rotate.jl")
include("transforms/translate.jl")
include("transforms/stretch.jl")
include("transforms/stdcoords.jl")
include("transforms/repair.jl")
include("transforms/smoothing.jl")
