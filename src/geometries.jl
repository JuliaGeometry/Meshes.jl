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
include("primitives/ball.jl")
include("primitives/sphere.jl")
include("primitives/cylinder.jl")

# ----------
# POLYTOPES
# ----------

"""
    Polytope{Dim,T,N}

We say that a geometry is a polytope when it is made of a collection of "flat" sides.
They are called polygon in 2D and polyhedron in 3D spaces. A polytope can be expressed
by an ordered set of points. These points (a.k.a. vertices) are connected into edges,
faces and cells in 3D. We follow the ordering conventions of the VTK project:
https://lorensen.github.io/VTKExamples/site/VTKBook/05Chapter5/#54-cell-types

Additionally, the following property must hold in order for a geometry to be considered
a polytope: the boundary of a (n+1)-polytope is a collection of n-polytopes, which may
have (n-1)-polytopes in common. See https://en.wikipedia.org/wiki/Polytope.

Meshing algorithms discretize geometries into polytope elements (e.g. triangles,
tetrahedrons, pyramids). Thus, the `Polytope` type can be used for dispatch in
functions that are agnostic to the mesh element type.
"""
abstract type Polytope{Dim,T} <: Geometry{Dim,T} end

(::Type{PL})(vertices::Vararg{P}) where {PL<:Polytope,P<:Point} = PL(SVector(vertices))
(::Type{PL})(vertices::AbstractVector{TP}) where {PL<:Polytope,TP<:Tuple} = PL(Point.(vertices))
(::Type{PL})(vertices::Vararg{TP}) where {PL<:Polytope,TP<:Tuple} = PL(SVector(vertices))

vertices(p::Polytope) = p.vertices

function Base.show(io::IO, p::Polytope)
  kind = nameof(typeof(p))
  vert = join(p.vertices, ", ")
  print(io, "$kind($vert)")
end

function Base.show(io::IO, ::MIME"text/plain", p::Polytope{Dim,T}) where {Dim,T}
  kind = nameof(typeof(p))
  lines = ["  └─$v" for v in p.vertices]
  println(io, "$kind{$Dim,$T}")
  print(io, join(lines, "\n"))
end

include("polytopes/segment.jl")
include("polytopes/triangle.jl")
include("polytopes/quadrangle.jl")
include("polytopes/pyramid.jl")
include("polytopes/tetrahedron.jl")
include("polytopes/hexahedron.jl")

# -------
# CHAINS
# -------

"""
    Chain(p1, p2, ..., pn)

A polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.
See https://en.wikipedia.org/wiki/Polygonal_chain.
"""
struct Chain{Dim,T,N} <: Geometry{Dim,T}
  vertices::NTuple{N,Point{Dim,T}}
end

Chain(points::Vararg{P,N}) where {N,P<:Point} = Chain(points)
Chain(points::NTuple{N,TP}) where {N,TP<:Tuple} = Chain(Point.(points))
Chain(points::Vararg{TP,N}) where {N,TP<:Tuple} = Chain(points)

vertices(c::Chain) = c.vertices

isclosed(c::Chain) = first(c.vertices) == last(c.vertices)

function Base.show(io::IO, c::Chain{Dim,T,N}) where {Dim,T,N}
  vert = join(c.vertices, ", ")
  print(io, "$N-chain($vert)")
end

function Base.show(io::IO, ::MIME"text/plain", c::Chain{Dim,T,N}) where {Dim,T,N}
  lines = ["  └─$v" for v in c.vertices]
  println(io, "$N-chain{$Dim,$T}")
  print(io, join(lines, "\n"))
end
