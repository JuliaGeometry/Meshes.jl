# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Polytope{N,Dim,T}

We say that a geometry is a N-polytope when it is a collection of "flat" sides
that constitue a `N`-dimensional subspace. They are called polygon and polyhedron
respectively for 2D (`N=2`) and 3D (`N=3`) subspaces, embedded in a `Dim`-dimensional
space. The parameter `N` is also known as the rank or parametric dimension of the
polytope: https://en.wikipedia.org/wiki/Abstract_polytope.

The term polytope expresses a particular combinatorial structure. A polyhedron,
for example, can be decomposed into faces. Each face can then be decomposed into
edges, and edges into vertices. Some conventions act as a mapping between vertices
and higher dimensional features (edges, faces, cells...), removing the need to
store all features. We follow the ordering conventions of the GMSH project:
https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering

Additionally, the following property must hold in order for a geometry to be considered
a polytope: the boundary of a (N+1)-polytope is a collection of N-polytopes, which may
have (N-1)-polytopes in common. See https://en.wikipedia.org/wiki/Polytope.
"""
abstract type Polytope{N,Dim,T} <: Geometry{Dim,T} end

(::Type{PL})(vertices::Vararg{P}) where {PL<:Polytope,P<:Point} = PL(SVector(vertices))
(::Type{PL})(vertices::AbstractVector{TP}) where {PL<:Polytope,TP<:Tuple} = PL(Point.(vertices))
(::Type{PL})(vertices::Vararg{TP}) where {PL<:Polytope,TP<:Tuple} = PL(collect(vertices))

# type aliases for convenience
const Polygon = Polytope{2}
const Polyhedron = Polytope{3}

"""
    paramdim(polytope)

Return the parametric dimension or rank of the polytope.
"""
paramdim(::Type{<:Polytope{N}}) where {N} = N

"""
    vertices(polytope)

Return the vertices of a `polytope`.
"""
vertices(p::Polytope) = p.vertices

"""
    nvertices(polytope)

Return the number of vertices in the `polytope`.
"""
nvertices(p::Polytope) = length(p.vertices)

"""
    facets(polytope, convention=GMSH)

Return the facets of a `polytope` according to some
ordering `convention`. Default convention is [`GMSH`](@ref).
See https://en.wikipedia.org/wiki/Facet_(geometry)
"""
function facets(p::Polytope{N}, convention=GMSH) where {N}
  faces(p, N-1, convention)
end

"""
    faces(polytope, rank, convention=GMSH)

Return the faces of the `polytope` of given `rank` according
to some ordering `convention`. Default convention is [`GMSH`](@ref).
"""
function faces(p::Polytope, rank, convention=GMSH)
  (materialize(ord, p.vertices) for ord in connectivities(typeof(p), rank, convention))
end

"""
    p1 == p2

Tells whether or not polytopes `p1` and `p2` are equal.
"""
==(p1::Polytope, p2::Polytope) = p1.vertices == p2.vertices

"""
    center(polytope)

Return the center of the `polytope`.
"""
center(p::Polytope) = Point(sum(coordinates.(p.vertices)) / length(p.vertices))

function Base.show(io::IO, p::Polytope)
  kind = nameof(typeof(p))
  vert = join(p.vertices, ", ")
  print(io, "$kind($vert)")
end

function Base.show(io::IO, ::MIME"text/plain", p::Polytope{N,Dim,T}) where {N,Dim,T}
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
include("polytopes/chain.jl")
include("polytopes/polyarea.jl")
