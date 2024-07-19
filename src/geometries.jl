# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Geometry{M,CRS}

A geometry in a given manifold `M` with point coordinates specified
in a coordinate reference system `CRS`.
"""
abstract type Geometry{M,CRS} end

Broadcast.broadcastable(g::Geometry) = Ref(g)

"""
    embeddim(geometry)

Return the number of dimensions of the space where the `geometry` is embedded.
"""
embeddim(::Type{<:Geometry{M,CRS}}) where {M,CRS} = CoordRefSystems.ndims(CRS)
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
crs(::Type{<:Geometry{M,CRS}}) where {M,CRS} = CRS
crs(g::Geometry) = crs(typeof(g))

"""
    manifold(geometry)

Return the manifold where the `geometry` is defined.
"""
manifold(::Type{<:Geometry{M,CRS}}) where {M,CRS} = M
manifold(g::Geometry) = manifold(typeof(g))

"""
    lentype(geometry)

Return the length type of the `geometry`.
"""
lentype(::Type{<:Geometry{CRS}}) where {CRS} = lentype(CRS)
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

function Base.convert(::Type{<:Quadrangle}, b::Box)
  checkdim(b, 2)
  Quadrangle(vertices(boundary(b))...)
end

function Base.convert(::Type{<:Hexahedron}, b::Box)
  checkdim(b, 3)
  Hexahedron(vertices(boundary(b))...)
end
