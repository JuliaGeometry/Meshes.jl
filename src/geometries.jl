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

We say that a geometry is a N-polytope when it is a collection of "flat" sides that constitue a N-dimensional subspace.
They are called polygon and polyhedron respectively for 2D (N=2) and 3D (N=3) subspaces, embedded in a `Dim`-dimensional space.
The term polytope expresses a particular combinatorial structure. A polyhedron, for example, can be decomposed into faces. Each face can then be decomposed into edges, and edges into vertices.
Some conventions act as a mapping between vertices and higher dimensional features (edges, faces, cells...), removing the need to store all features. We follow the [ordering conventions](https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering) of the GMSH project.

Additionally, the following property must hold in order for a geometry to be considered
a polytope: the boundary of a (n+1)-polytope is a collection of n-polytopes, which may
have (n-1)-polytopes in common. For more information, see the [Wikipedia](https://en.wikipedia.org/wiki/Polytope) page.
"""
abstract type Polytope{N,Dim,T} <: Geometry{Dim,T} end

const Polygon = Polytope{2}
const Polyhedron = Polytope{3}

paramdim(::Type{<: Polytope{N}}) where {N} = N

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
    center(polytope)

Return the center of the `polytope`.
"""
center(p::Polytope) = Point(sum(coordinates.(p.vertices)) / length(p.vertices))

"""
    Face{N, Dim,T}

We say that a polytope is a face when it can be used as an element in a finite element
mesh (e.g. segments, triangles, tetrahedrons). The rank of the face reflects the actual
parametric dimension of the polytope. For example, a segment is a 1-face, a triangle is a
2-face and a tetrahedron is a 3-face. See https://en.wikipedia.org/wiki/Abstract_polytope.
"""
abstract type Face{N,Dim,T} <: Polytope{N,Dim,T} end

(::Type{F})(vertices::Vararg{P}) where {F<:Face,P<:Point} = F(SVector(vertices))
(::Type{F})(vertices::AbstractVector{TP}) where {F<:Face,TP<:Tuple} = F(Point.(vertices))
(::Type{F})(vertices::Vararg{TP}) where {F<:Face,TP<:Tuple} = F(collect(vertices))

==(f1::Face, f2::Face) = f1.vertices == f2.vertices

function Base.show(io::IO, f::Polytope)
  kind = nameof(typeof(f))
  vert = join(f.vertices, ", ")
  print(io, "$kind($vert)")
end

function Base.show(io::IO, ::MIME"text/plain", f::Polytope{N,Dim,T}) where {N,Dim,T}
  kind = nameof(typeof(f))
  lines = ["  └─$v" for v in f.vertices]
  println(io, "$kind{$Dim,$T}")
  print(io, join(lines, "\n"))
end

nfaces(p, n::Val; ordering=GMSH) = (materialize(c, p.vertices) for c in connectivity(typeof(p), ordering, n))
nfaces(p, ::Val{0}) = p.vertices
facets(p::Polytope{N}; ordering=GMSH) where {N} = nfaces(p, Val(N-1); ordering)


# ----------------
# CONVENTION MODEL
# ----------------

include("conventions.jl")

# ----------------
# COMMON POLYTOPES
# ----------------

include("common_polytopes/segment.jl")
include("common_polytopes/triangle.jl")
include("common_polytopes/quadrangle.jl")
include("common_polytopes/pyramid.jl")
include("common_polytopes/tetrahedron.jl")
include("common_polytopes/hexahedron.jl")

# -----------
# CONVENTIONS
# -----------

include("conventions/gmsh.jl")

# -------
# CHAINS
# -------

include("chains.jl")

# ---------
# POLYSURFACES
# ---------

include("polysurfaces.jl")
