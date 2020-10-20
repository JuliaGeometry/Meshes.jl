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

We say that a geometry is a "primitive" when it can be expressed as a single
entity with no parts. For example, a sphere is a primitive described in terms
of a mathematical expression involving a metric and a radius.

Primitives can be discretized into a collection of finite elements with meshing
algorithms such as triangulation.
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
by an ordered set of points. These points (a.k.a. vertices) are connected into edges,
faces and cells in 3D. The number of points `N` in the polytope is known at compile
time. For example, a triangle has N=3 points. See https://en.wikipedia.org/wiki/Polytope.
"""
abstract type Polytope{Dim,T,N} <: Geometry{Dim,T} end

function (::Type{P})(vertices...) where {P<:Polytope}
    PT  = eltype(vertices)
    Dim = ndims(PT)
    T   = coordtype(PT)
    return P{Dim,T}(vertices)
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
include("polytopes/tetrahedron.jl")

# -----------------
# OTHER GEOMETRIES
# -----------------

# TODO: review this
include("faces.jl")

"""
    AbstractMesh

A mesh is a collection of Polytope elements.
"""
abstract type AbstractMesh{Element<:Polytope} <: AbstractVector{Element} end

"""
    Mesh <: AbstractVector{Element}

The conrecte AbstractMesh implementation
"""
struct Mesh{Dim,T,Element<:Polytope{Dim,T},V<:AbstractVector{Element}} <: AbstractMesh{Element}
    simplices::V # usually a FaceView, to connect a set of points via a set of faces.
end

Base.size(mesh::Mesh) = size(getfield(mesh, :simplices))
Base.getindex(mesh::Mesh, i::Integer) = getfield(mesh, :simplices)[i]

function Mesh(points::AbstractVector{<:Point},
              faces::AbstractVector{<:AbstractFace})
    return Mesh(connect(points, faces))
end

function Mesh(points::AbstractVector{<:Point},
              faces::AbstractVector{<:Integer},
              facetype=TriangleFace, skip=1)
    return Mesh(connect(points, connect(faces, facetype, skip)))
end
