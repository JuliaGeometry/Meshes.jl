# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Polytope{K,Dim,T}

We say that a geometry is a K-polytope when it is a collection of "flat" sides
that constitue a `K`-dimensional subspace. They are called polygon and polyhedron
respectively for 2D (`K=2`) and 3D (`K=3`) subspaces, embedded in a `Dim`-dimensional
space. The parameter `K` is also known as the rank or parametric dimension of the
polytope: https://en.wikipedia.org/wiki/Abstract_polytope.

The term polytope expresses a particular combinatorial structure. A polyhedron,
for example, can be decomposed into faces. Each face can then be decomposed into
edges, and edges into vertices. Some conventions act as a mapping between vertices
and higher dimensional features (edges, faces, cells...), removing the need to
store all features.

Additionally, the following property must hold in order for a geometry to be considered
a polytope: the boundary of a (K+1)-polytope is a collection of K-polytopes, which may
have (K-1)-polytopes in common. See https://en.wikipedia.org/wiki/Polytope.

### Notes

- Type aliases are `Polygon`, `Polyhedron`.
"""
abstract type Polytope{K,Dim,T} <: Geometry{Dim,T} end

(::Type{PL})(vertices::Vararg{P}) where {PL<:Polytope,P<:Point} = PL(SVector(vertices))
(::Type{PL})(vertices::AbstractVector{TP}) where {PL<:Polytope,TP<:Tuple} = PL(Point.(vertices))
(::Type{PL})(vertices::Vararg{TP}) where {PL<:Polytope,TP<:Tuple} = PL(collect(vertices))

"""
    paramdim(polytope)

Return the parametric dimension or rank of the polytope.
"""
paramdim(::Type{<:Polytope{K}}) where {K} = K

"""
    vertices(polytope)

Return the vertices of a `polytope`.
"""
vertices(p::Polytope) = p.vertices

"""
    nvertices(polytope)

Return the number of vertices in the `polytope`.
"""
nvertices(p::Polytope) = length(vertices(p))

"""
    p1 == p2

Tells whether or not polytopes `p1` and `p2` are equal.
"""
==(p1::Polytope, p2::Polytope) = vertices(p1) == vertices(p2)

"""
    centroid(polytope)

Return the centroid of the `polytope`.
"""
centroid(p::Polytope) = Point(sum(coordinates.(vertices(p))) / length(vertices(p)))

function Base.show(io::IO, p::Polytope)
  kind = prettyname(p)
  vert = join(vertices(p), ", ")
  print(io, "$kind($vert)")
end

function Base.show(io::IO, ::MIME"text/plain", p::Polytope{K,Dim,T}) where {K,Dim,T}
  kind = prettyname(p)
  lines = ["  └─$v" for v in vertices(p)]
  println(io, "$kind{$Dim,$T}")
  print(io, join(lines, "\n"))
end

prettyname(p::Polytope) = prettyname(typeof(p))
function prettyname(PL::Type{<:Polytope})
  n = string(PL)
  i = findfirst('{', n)
  isnothing(i) ? n : n[1:i-1]
end

"""
    Polygon{Dim,T}

A polygon is a 2-polytope, i.e. a polytope with parametric dimension 2.
"""
const Polygon = Polytope{2}

"""
    area(polygon)

Return the area of the `polygon`.
"""
area(p::Polygon) = measure(p)

"""
    Polyhedron{Dim,T}

A polyhedron is a 3-polytope, i.e. a polytope with parametric dimension 3.
"""
const Polyhedron = Polytope{3}

"""
   volume(polyhedron)

Return the volume of the `polyhedron`.
"""
volume(p::Polyhedron) = measure(p)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("polytopes/segment.jl")
include("polytopes/ngon.jl")
include("polytopes/chain.jl")
include("polytopes/polyarea.jl")
include("polytopes/pyramid.jl")
include("polytopes/tetrahedron.jl")
include("polytopes/hexahedron.jl")
