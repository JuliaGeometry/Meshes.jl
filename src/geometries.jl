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
    Polytope{Dim,T}

We say that a geometry is a polytope when it is made of a collection of "flat" sides.
They are called polygon in 2D and polyhedron in 3D spaces. A polytope can be expressed
by an ordered set of points. These points (a.k.a. vertices) are connected into edges,
faces and cells in 3D. We follow the ordering conventions of the GMSH project:
https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering

Additionally, the following property must hold in order for a geometry to be considered
a polytope: the boundary of a (n+1)-polytope is a collection of n-polytopes, which may
have (n-1)-polytopes in common. See https://en.wikipedia.org/wiki/Polytope.
"""
abstract type Polytope{Dim,T} <: Geometry{Dim,T} end

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
    Face{Dim,T}

We say that a polytope is a face when it can be used as an element in a finite element
mesh (e.g. segments, triangles, tetrahedrons). The rank of the face reflects the actual
parametric dimension of the polytope. For example, a segment is a 1-face, a triangle is a
2-face and a tetrahedron is a 3-face. See https://en.wikipedia.org/wiki/Abstract_polytope.
"""
abstract type Face{Dim,T} <: Polytope{Dim,T} end

(::Type{F})(vertices::Vararg{P}) where {F<:Face,P<:Point} = F(SVector(vertices))
(::Type{F})(vertices::AbstractVector{TP}) where {F<:Face,TP<:Tuple} = F(Point.(vertices))
(::Type{F})(vertices::Vararg{TP}) where {F<:Face,TP<:Tuple} = F(collect(vertices))

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

# -------
# CHAINS
# -------

"""
    Chain(p1, p2, ..., pn)

A chain from a sequence of points `p1`, `p2`, ..., `pn`.
See https://en.wikipedia.org/wiki/Polygonal_chain.
"""
struct Chain{Dim,T} <: Polytope{Dim,T}
  vertices::Vector{Point{Dim,T}}
end

Chain(points::Vararg{P}) where {P<:Point} = Chain(collect(points))
Chain(points::AbstractVector{TP}) where {TP<:Tuple} = Chain(Point.(points))
Chain(points::Vararg{TP}) where {TP<:Tuple} = Chain(collect(points))

"""
    isclosed(chain)

Tells whether or not the chain is closed.
A closed chain is also known as a ring.
"""
isclosed(c::Chain) = first(c.vertices) == last(c.vertices)

function Base.show(io::IO, c::Chain)
  N = length(c.vertices)
  print(io, "$N-chain")
end

function Base.show(io::IO, ::MIME"text/plain", c::Chain{Dim,T}) where {Dim,T}
  N = length(c.vertices)
  lines = ["  └─$v" for v in c.vertices]
  println(io, "$N-chain{$Dim,$T}")
  print(io, join(lines, "\n"))
end

# ---------
# POLYGONS
# ---------

include("polygons.jl")
