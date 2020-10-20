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
include("polytopes/tetrahedron.jl")

# -----------------
# OTHER GEOMETRIES
# -----------------

# TODO: review this
include("faces.jl")

function coordinates(lines::AbstractArray{<:Line})
    return if lines isa Base.ReinterpretArray
        return coordinates(lines.parent)
    else
        result = PointType[]
        for line in lines
            append!(result, coordinates(line))
        end
        return result
    end
end

"""
    LineString(points::AbstractVector{<:Point})

A LineString is a geometry of connected line segments
"""
struct LineString{Dim,T,V<:AbstractVector{Line{Dim,T}}} <: AbstractVector{Line{Dim,T}}
    points::V
end

coordinates(x::LineString) = coordinates(x.points)
Base.copy(x::LineString) = LineString(copy(x.points))
Base.size(x::LineString) = size(getfield(x, :points))
Base.getindex(x::LineString, i) = getindex(getfield(x, :points), i)

"""
    LineString(points::AbstractVector{<:Point}, skip = 1)

Creates a LineString from a vector of points
With `skip == 1`, the default, it will connect the line like this:
```julia
points = Point[a, b, c, d]
linestring = LineString(points)
@assert linestring == LineString([a => b, b => c, c => d])
```
"""
function LineString(points::AbstractVector{<:Point{Dim,T}}, skip::Int=1) where {Dim,T}
    return LineString(connect(points, Line{Dim,T}, skip))
end

function LineString(points::AbstractVector{<:Pair{Point{N,T},Point{N,T}}}) where {N,T}
    return LineString(reinterpret(Line{N,T}, points))
end

"""
    Polygon(exterior::AbstractVector{<:Point})
    Polygon(exterior::AbstractVector{<:Point}, interiors::Vector{<:AbstractVector{<:Point}})

"""
struct Polygon{Dim,T,L<:AbstractVector{Line{Dim,T}},V<:AbstractVector{L}}
    exterior::L
    interiors::V
end

Base.copy(x::Polygon) = Polygon(copy(x.exterior), copy(x.interiors))

function Base.:(==)(a::Polygon, b::Polygon)
    return (a.exterior == b.exterior) && (a.interiors == b.interiors)
end

Polygon(exterior::L) where {L<:AbstractVector{<:Line}} = Polygon(exterior, L[])

function Polygon(exterior::AbstractVector{P}, skip::Int=1) where {P<:Point{Dim,T}} where {Dim,T}
    return Polygon(LineString(exterior, skip))
end

function Polygon(exterior::AbstractVector{P}, faces::AbstractVector{<:Integer}, skip::Int=1) where {P<:Point{Dim,T}} where {Dim,T}
    return Polygon(LineString(exterior, faces, skip))
end

function Polygon(exterior::AbstractVector{P}, interior::AbstractVector{<:AbstractVector{P}}) where {P<:Point{Dim,T}} where {Dim,T}
    ext = LineString(exterior)
    # We need to take extra care for empty interiors, since
    # if we just map over it, it won't infer the element type correctly!
    int = typeof(ext)[]
    foreach(x -> push!(int, LineString(x)), interior)
    return Polygon(ext, int)
end

Base.ndims(::Polygon{Dim,T}) where {Dim,T} = Dim

function coordinates(polygon::Polygon{N,T}) where {N,T}
    exterior = coordinates(polygon.exterior)
    if isempty(polygon.interiors)
        return exterior
    else
        result = Point{N,T}[]
        append!(result, exterior)
        foreach(x -> append!(result, coordinates(x)), polygon.interiors)
        return result
    end
end

"""
    MultiPolygon(polygons)
"""
struct MultiPolygon{Dim,T,Element<:Polygon{Dim,T},A<:AbstractVector{Element}} <: AbstractVector{Element}
    polygons::A
end

Base.getindex(mp::MultiPolygon, i) = mp.polygons[i]
Base.size(mp::MultiPolygon) = size(mp.polygons)

struct MultiLineString{Dim,T,Element<:LineString{Dim,T},A<:AbstractVector{Element}} <: AbstractVector{Element}
    linestrings::A
end

function MultiLineString(linestrings::AbstractVector{L}) where {L<:AbstractVector{Line{Dim,T}}} where {Dim,T}
    return MultiLineString(linestrings)
end

Base.getindex(ms::MultiLineString, i) = ms.linestrings[i]
Base.size(ms::MultiLineString) = size(ms.linestrings)

"""
    MultiPoint(points::AbstractVector{<:Point})

A collection of points
"""
struct MultiPoint{Dim,T,P<:Point{Dim,T},A<:AbstractVector{P}} <: AbstractVector{P}
    points::A
end

Base.getindex(mpt::MultiPoint, i) = mpt.points[i]
Base.size(mpt::MultiPoint) = size(mpt.points)

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
