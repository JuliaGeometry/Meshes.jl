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

# -----
# SETS
# -----

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
    facets(polytope)

Return the facets of a `polytope`.
See https://en.wikipedia.org/wiki/Facet_(geometry)
"""
function facets end

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

# -------
# CHAINS
# -------

include("chains.jl")

# -------------
# POLYSURFACES
# -------------

include("polysurfaces.jl")
