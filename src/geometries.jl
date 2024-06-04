# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Geometry{Dim,CRS}

A geometry embedded in a `Dim`-dimensional space with given coordinate reference system `CRS`.
"""
abstract type Geometry{Dim,CRS} end

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
    crs(geometry)

Return the coordinate reference system (CRS) of the `geometry`.
"""
crs(::Type{<:Geometry{Dim,CRS}}) where {Dim,CRS} = CRS
crs(g::Geometry) = crs(typeof(g))

"""
    lentype(geometry)

Return the length type of the `geometry`.
"""
lentype(::Type{<:Geometry{Dim,CRS}}) where {Dim,CRS} = lentype(CRS)
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

# -------
# RANDOM
# -------

Random.rand(::Type{G}) where {G<:Geometry} = rand(Random.default_rng(), G)
Random.rand(::Type{G}, n::Int) where {G<:Geometry} = rand(Random.default_rng(), G, n)
Random.rand(rng::Random.AbstractRNG, ::Type{G}, n::Int) where {G<:Geometry} = [rand(rng, G) for _ in 1:n]

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
