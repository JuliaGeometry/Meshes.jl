# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Geometry{Dim,T}

A geometry embedded in a `Dim`-dimensional space with coordinates of type `T`.
"""
abstract type Geometry{Dim,T} end

"""
    ndims(geometry)

Return the number of dimensions of the space where the `geometry` is embedded.
"""
Base.ndims(::Geometry{Dim,T}) where {Dim,T} = Dim

"""
    coordtype(geometry)

Return the machine type of each coordinate used to describe the `geometry`.
"""
coordtype(::Geometry{Dim,T}) where {Dim,T}  = T

# -----------
# PRIMITIVES
# -----------

"""
    Primitive{Dim,T}

We say that a geometry is a primitive when it can be expressed as a single
entity with no parts (a.k.a. atomic). For example, a sphere is a primitive
described in terms of a mathematical expression involving a metric and a radius.
See https://en.wikipedia.org/wiki/Geometric_primitive.
"""
abstract type Primitive{Dim,T} <: Geometry{Dim,T} end

include("primitives/box.jl")
include("primitives/sphere.jl")
include("primitives/cylinder.jl")

# ----------
# POLYTOPES
# ----------

"""
    Polytope{Dim,T,N}

We say that a geometry is a polytope when it is made of a collection of "flat" sides.
They are called polygon in 2D and polyhedron in 3D spaces. A polytope can be expressed
by an ordered set of `N` points. These points (a.k.a. vertices) are connected into edges,
faces and cells in 3D. We follow the ordering conventions of the VTK project:
https://lorensen.github.io/VTKExamples/site/VTKBook/05Chapter5/#54-cell-types

Additionally, the following property must hold in order for a geometry to be considered
a polytope: the boundary of a (n+1)-polytope is a collection of n-polytopes, which may
have (n-1)-polytopes in common. See https://en.wikipedia.org/wiki/Polytope.

Meshing algorithms discretize geometries into polytope elements (e.g. triangles,
tetrahedrons, pyramids). Thus, the `Polytope` type can be used for dispatch in
functions that are agnostic to the mesh element type.
"""
abstract type Polytope{Dim,T,N} <: Geometry{Dim,T} end

function (::Type{PT})(vertices...) where {PT<:Polytope}
    P = eltype(vertices)
    return PT{ndims(P),coordtype(P)}(vertices)
end

Base.getindex(p::Polytope, i) = getindex(p.vertices, i)
Base.length(::Type{<:Polytope{Dim,T,N}}) where {Dim,T,N} = N
Base.length(p::Polytope) = length(typeof(p))

Base.iterate(p::Polytope, i) = iterate(p.vertices, i)
Base.iterate(p::Polytope) = iterate(p.vertices)

vertices(p::Polytope) = p.vertices

include("polytopes/line.jl")
include("polytopes/triangle.jl")
include("polytopes/quadrangle.jl")
include("polytopes/pyramid.jl")
include("polytopes/tetrahedron.jl")
include("polytopes/hexahedron.jl")

# -------
# CHAINS
# -------

"""
    Chain(g1, g2, ..., gn)

Construct a chain of geometries `g1`, `g2`, ..., `gn`.
See https://en.wikipedia.org/wiki/Chain_(algebraic_topology).

    Chain(p1, p2, ..., pn)

Alternatively, construct a chain of line geometries
from a sequence of points `p1`, `p2`, ..., `pn`.
"""
struct Chain{Dim,T,N,G<:Geometry{Dim,T}} <: Geometry{Dim,T}
    geometries::NTuple{N,G}
end

Chain(geometries::Vararg{G,N}) where {N,G<:Geometry} = Chain(geometries)

Chain(points::NTuple{N,P}) where {N,P<:Point} =
    Chain(ntuple(i -> Line(points[i], points[i+1]), N-1))

Chain(points::Vararg{P,N}) where {N,P<:Point} = Chain(points)

Base.getindex(c::Chain, i) = getindex(c.geometries, i)
Base.length(::Type{<:Chain{Dim,T,N}}) where {Dim,T,N} = N
Base.length(c::Chain) = length(typeof(c))

# -----------------
# OTHER GEOMETRIES
# -----------------

# TODO: review this
include("faces.jl")

struct Mesh{Dim,T,E<:Polytope{Dim,T},V<:AbstractVector{E}} <: Geometry{Dim,T}
    elements::V # usually a FaceView, to connect a set of points via a set of faces.
end

Mesh(points::AbstractVector{<:Point},
     faces::AbstractVector{<:AbstractFace}) = Mesh(connect(points, faces))

Mesh(points::AbstractVector{<:Point},
     faces::AbstractVector{<:Integer},
     facetype=TriangleFace, skip=1) = Mesh(connect(points, connect(faces, facetype, skip)))

Base.getindex(m::Mesh, i) = getindex(m.elements, i)
Base.length(m::Mesh) = length(m.elements)

Base.iterate(m::Mesh, i) = iterate(m.elements, i)
Base.iterate(m::Mesh) = iterate(m.elements)

elements(m::Mesh) = m.elements

faces(m::Mesh) = faces(m.elements)

coordinates(m::Mesh) = coordinates(m.elements)

volume(m::Mesh) = sum(volume, m)
