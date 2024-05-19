# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Geometry{Dim}

A geometry embedded in a `Dim`-dimensional space.
"""
abstract type Geometry{Dim} end

Broadcast.broadcastable(g::Geometry) = Ref(g)

"""
    embeddim(geometry)

Return the number of dimensions of the space where the `geometry` is embedded.
"""
embeddim(::Type{<:Geometry{Dim}}) where {Dim} = Dim
embeddim(g::Geometry) = embeddim(typeof(g))

"""
    paramdim(geometry)

Return the number of parametric dimensions of the `geometry`. For example, a
sphere embedded in 3D has 2 parametric dimensions (polar and azimuthal angles).

See also [`isparametrized`](@ref).
"""
paramdim(g::Geometry) = paramdim(typeof(g))

"""
    lentype(geometry)

Return the length type of the `geometry`.
"""
lentype(g::Geometry) = lentype(typeof(g))

"""
    centroid(geometry)

Return the centroid of the `geometry`.
"""
centroid(g::Geometry) = center(g)

"""
    extrema(geometry)

Return the top left and bottom right corners of the
bounding box of the `geometry`.
"""
Base.extrema(g::Geometry) = extrema(boundingbox(g))

# -----------
# IO METHODS
# -----------

Base.summary(io::IO, geom::Geometry) = print(io, prettyname(geom))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("primitives.jl")
include("polytopes.jl")
include("multigeoms.jl")

# ------------
# CONVERSIONS
# ------------

Base.convert(::Type{<:Quadrangle}, b::Box{2}) = Quadrangle(vertices(boundary(b))...)

Base.convert(::Type{<:Hexahedron}, b::Box{3}) = Hexahedron(vertices(boundary(b))...)
