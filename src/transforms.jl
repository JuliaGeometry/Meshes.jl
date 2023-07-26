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
    applycoord(transform, object)

Recursively apply coordinate `transform` on `object`.
This function is intended for developers of new
[`CoordinateTransform`](@ref).
"""
function applycoord end

# --------------------
# TRANSFORM FALLBACKS
# --------------------

apply(t::CoordinateTransform, g) = applycoord(t, g), nothing

revert(t::CoordinateTransform, g, c) = applycoord(inv(t), g)

# apply transform recursively
applycoord(t::CoordinateTransform, g::G) where {G<:Union{Geometry,Domain}} =
  G((applycoord(t, getfield(g, n)) for n in fieldnames(G))...)

# stop recursion at non-geometric types
applycoord(::CoordinateTransform, x) = x

# special treatment for lists of geometries
applycoord(t::CoordinateTransform, g::NTuple{<:Any,<:Geometry}) = map(gᵢ -> applycoord(t, gᵢ), g)
applycoord(t::CoordinateTransform, g::AbstractVector{<:Geometry}) = map(gᵢ -> applycoord(t, gᵢ), g)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transforms/rotate.jl")
include("transforms/translate.jl")
include("transforms/stretch.jl")
include("transforms/stdcoords.jl")
include("transforms/repair.jl")
include("transforms/smoothing.jl")
