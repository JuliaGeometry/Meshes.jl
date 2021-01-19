# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Geometry{Dim,T}

A geometry embedded in a `Dim`-dimensional space with coordinates of type `T`.
"""
abstract type Geometry{Dim,T} end

"""
    embeddim(geometry)

Return the number of dimensions of the space where the `geometry` is embedded.
"""
embeddim(::Type{<:Geometry{Dim,T}}) where {Dim,T} = Dim
embeddim(g::Geometry) = embeddim(typeof(g))

"""
    paramdim(geometry)

Return the number of parametric dimensions of the `geometry`. For example, a
sphere embedded in 3D has 2 parametric dimension (polar and azimuthal angles).
"""
paramdim(g::Geometry) = paramdim(typeof(g))

"""
    coordtype(geometry)

Return the machine type of each coordinate used to describe the `geometry`.
"""
coordtype(::Type{<:Geometry{Dim,T}}) where {Dim,T} = T
coordtype(g::Geometry) = coordtype(typeof(g))

"""
    measure(geometry)

Return the measure of the `geometry`, i.e. the length, area, or volume.
"""
function measure end

"""
    boundary(geometry)

Return the boundary of the `geometry`.
"""
function boundary end

# ----------------
# IMPLEMENTATIONS
# ----------------
include("primitives.jl")
include("polytopes.jl")

"""
    GeometrySet(geometries)

A set of `geometries` seen as a single geometry.

## Example

In a geographic map, countries can be described with
multiple polygons (a.k.a. MultiPolygon).
"""
struct GeometrySet{Dim,T,G<:Geometry{Dim,T}} <: Geometry{Dim,T}
  geometries::Vector{G}
end
