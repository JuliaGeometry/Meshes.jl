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
    coordtype(geometry)

Return the machine type of each coordinate used to describe the `geometry`.
"""
coordtype(::Type{<:Geometry{Dim,T}}) where {Dim,T}  = T
coordtype(g::Geometry) = coordtype(typeof(g))

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

# -------
# CHAINS
# -------

"""
    Chain(p1, p2, ..., pn)

A chain from a sequence of points `p1`, `p2`, ..., `pn`.
See https://en.wikipedia.org/wiki/Polygonal_chain.
"""
struct Chain{Dim,T,N} <: Geometry{Dim,T}
  vertices::NTuple{N,Point{Dim,T}}
end

Chain(points::Vararg{P,N}) where {N,P<:Point} = Chain(points)
Chain(points::NTuple{N,TP}) where {N,TP<:Tuple} = Chain(Point.(points))
Chain(points::Vararg{TP,N}) where {N,TP<:Tuple} = Chain(points)

"""
    vertices(chain)

Return the vertices of a chain.
"""
vertices(c::Chain) = c.vertices

"""
    isclosed(chain)

Tells whether or not the chain is closed.
A closed chain is also known as a ring.
"""
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

# ----------
# POLYTOPES
# ----------

"""
    Polytope{Dim,T}

We say that a geometry is a polytope when it is made of a collection of "flat" sides.
They are called polygon in 2D and polyhedron in 3D spaces. A polytope can be expressed
by an ordered set of points. These points (a.k.a. vertices) are connected into edges,
faces and cells in 3D. We follow the ordering conventions of the VTK project:
https://lorensen.github.io/VTKExamples/site/VTKBook/05Chapter5/#54-cell-types

Additionally, the following property must hold in order for a geometry to be considered
a polytope: the boundary of a (n+1)-polytope is a collection of n-polytopes, which may
have (n-1)-polytopes in common. See https://en.wikipedia.org/wiki/Polytope.
"""
abstract type Polytope{Dim,T} <: Geometry{Dim,T} end

"""
    Face{Dim,T,Rank}

We say that a polytope is a face when it can be used as an element in a finite element
mesh (e.g. segments, triangles, tetrahedrons). The `Rank` of the face reflects the actual
parametric dimension of the polytope. For example, a segment is a 1-face, a triangle is a
2-face and a tetrahedron is a 3-face. See https://en.wikipedia.org/wiki/Abstract_polytope.
"""
abstract type Face{Dim,T,Rank} <: Polytope{Dim,T} end

(::Type{F})(vertices::Vararg{P}) where {F<:Face,P<:Point} = F(SVector(vertices))
(::Type{F})(vertices::AbstractVector{TP}) where {F<:Face,TP<:Tuple} = F(Point.(vertices))
(::Type{F})(vertices::Vararg{TP}) where {F<:Face,TP<:Tuple} = F(SVector(vertices))

"""
    rank(face)

Return the rank of the face.
"""
rank(::Type{<:Face{Dim,T,Rank}}) where {Dim,T,Rank} = Rank
rank(f::Face) = rank(typeof(f))

"""
    vertices(face)

Return the vertices of a face.
"""
vertices(f::Face) = f.vertices

==(f1::Face, f2::Face) = f1.vertices == f2.vertices

function Base.show(io::IO, f::Face)
  kind = nameof(typeof(f))
  vert = join(f.vertices, ", ")
  print(io, "$kind($vert)")
end

function Base.show(io::IO, ::MIME"text/plain", f::Face{Dim,T}) where {Dim,T}
  kind = nameof(typeof(f))
  lines = ["  └─$v" for v in f.vertices]
  println(io, "$kind{$Dim,$T}")
  print(io, join(lines, "\n"))
end

include("faces/segment.jl")
include("faces/triangle.jl")
include("faces/quadrangle.jl")
include("faces/pyramid.jl")
include("faces/tetrahedron.jl")
include("faces/hexahedron.jl")

"""
    Polygon(outer, [inner1, inner2, ...])

A polygon with `outer` ring, and optional inner
rings `inner1`, `inner2`, ...
"""
struct Polygon{Dim,T,C<:Chain{Dim,T}} <: Polytope{Dim,T}
  outer::C
  inners::Vector{C}

  function Polygon{Dim,T,C}(outer, inners) where {Dim,T,C}
    @assert isclosed(outer) "invalid outer ring"
    @assert all(isclosed.(inners)) "invalid inner rings"
    new(outer, inners)
  end
end

Polygon(outer::Chain{Dim,T}, inners) where {Dim,T} =
  Polygon{Dim,T,Chain{Dim,T}}(outer, inners)
Polygon(outer::C) where {C<:Chain} = Polygon(outer, C[])

function Base.show(io::IO, p::Polygon)
  outer = p.outer
  inner = isempty(p.inners) ? "" : ", "*join(p.inners, ", ")
  print(io, "Polygon($outer$inner)")
end

function Base.show(io::IO, ::MIME"text/plain", p::Polygon{Dim,T}) where {Dim,T}
  outer = "    └─$(p.outer)"
  inner = ["    └─$v" for v in p.inners]
  println(io, "Polygon{$Dim,$T}")
  println(io, "  outer")
  if isempty(inner)
    print(io, outer)
  else
    println(io, outer)
    println(io, "  inner")
    print(io, join(inner, "\n"))
  end
end
