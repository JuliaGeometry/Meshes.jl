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
    Polytope{Dim,T}

We say that a geometry is a polytope when it is made of a collection of "flat" sides.
They are called polygon in 2D and polyhedron in 3D spaces. A polytope can be expressed
by an ordered set of points. These points (a.k.a. vertices) are connected into edges,
faces and cells in 3D. See https://en.wikipedia.org/wiki/Polytope.
"""
abstract type Polytope{Dim,T} <: Geometry{Dim,T} end

@propagate_inbounds Base.getindex(x::Polytope, i::Integer) = coordinates(x)[i]
@propagate_inbounds Base.iterate(x::Polytope) = iterate(coordinates(x))
@propagate_inbounds Base.iterate(x::Polytope, i) = iterate(coordinates(x), i)

struct Ngon{Dim,T,N,P<:Point{Dim,T}} <: Polytope{Dim,T}
    points::SVector{N,P}
end

const NNgon{N} = Ngon{Dim,T,N,P} where {Dim,T,P}

function (::Type{<:NNgon{N}})(points::Vararg{P,N}) where {P<:Point{Dim,T},N} where {Dim,T}
    return Ngon{Dim,T,N,P}(SVector(points))
end

coordinates(x::Ngon) = x.points
Base.length(::Type{<:NNgon{N}}) where {N} = N
Base.length(::NNgon{N}) where {N} = N

# TODO: review this
include("faces.jl")

const LineP{Dim,T,P<:Point{Dim,T}} = Ngon{Dim,T,2,P}
const Line{Dim,T} = LineP{Dim,T,Point{Dim,T}}

const TriangleP{Dim,T,P<:Point{Dim,T}} = Ngon{Dim,T,3,P}
const Triangle{Dim,T} = TriangleP{Dim,T,Point{Dim,T}}

function coordinates(lines::AbstractArray{LineP{Dim,T,PointType}}) where {Dim,T,PointType}
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
struct LineString{Dim,T,P<:Point,V<:AbstractVector{<:LineP{Dim,T,P}}} <: AbstractVector{LineP{Dim,T,P}}
    points::V
end

coordinates(x::LineString) = coordinates(x.points)
Base.copy(x::LineString) = LineString(copy(x.points))
Base.size(x::LineString) = size(getfield(x, :points))
Base.getindex(x::LineString, i) = getindex(getfield(x, :points), i)

function LineString(points::AbstractVector{LineP{Dim,T,P}}) where {Dim,T,P}
    return LineString{Dim,T,P,typeof(points)}(points)
end

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
function LineString(points::AbstractVector{<:Point}, skip=1)
    return LineString(connect(points, LineP, skip))
end

function LineString(points::AbstractVector{<:Pair{P,P}}) where {P<:Point{N,T}} where {N,T}
    return LineString(reinterpret(LineP{N,T,P}, points))
end

function LineString(points::AbstractVector{<:Point}, faces::AbstractVector{<:LineFace})
    return LineString(connect(points, faces))
end

"""
    LineString(points::AbstractVector{<:Point}, indices::AbstractVector{<: Integer}, skip = 1)

Creates a LineString from a vector of points and an index list.
With `skip == 1`, the default, it will connect the line like this:


```julia
points = Point[a, b, c, d]; faces = [1, 2, 3, 4]
linestring = LineString(points, faces)
@assert linestring == LineString([a => b, b => c, c => d])
```
To make a segmented line, set skip to 2
```julia
points = Point[a, b, c, d]; faces = [1, 2, 3, 4]
linestring = LineString(points, faces, 2)
@assert linestring == LineString([a => b, c => d])
```
"""
function LineString(points::AbstractVector{<:Point}, indices::AbstractVector{<:Integer}, skip=1)
    faces = connect(indices, LineFace, skip)
    return LineString(points, faces)
end

"""
    Polygon(exterior::AbstractVector{<:Point})
    Polygon(exterior::AbstractVector{<:Point}, interiors::Vector{<:AbstractVector{<:Point}})

"""
struct Polygon{Dim,T,P<:Point{Dim,T},L<:AbstractVector{<:LineP{Dim,T,P}},V<:AbstractVector{L}} <: Polytope{Dim,T}
    exterior::L
    interiors::V
end

Base.copy(x::Polygon) = Polygon(copy(x.exterior), copy(x.interiors))

function Base.:(==)(a::Polygon, b::Polygon)
    return (a.exterior == b.exterior) && (a.interiors == b.interiors)
end

function Polygon(exterior::E, interiors::AbstractVector{E}) where {E<:AbstractVector{LineP{Dim,T,P}}} where {Dim,T,P}
    return Polygon{Dim,T,P,typeof(exterior),typeof(interiors)}(exterior, interiors)
end

Polygon(exterior::L) where {L<:AbstractVector{<:LineP}} = Polygon(exterior, L[])

function Polygon(exterior::AbstractVector{P}, skip::Int=1) where {P<:Point{Dim,T}} where {Dim,T}
    return Polygon(LineString(exterior, skip))
end

function Polygon(exterior::AbstractVector{P}, faces::AbstractVector{<:Integer}, skip::Int=1) where {P<:Point{Dim,T}} where {Dim,T}
    return Polygon(LineString(exterior, faces, skip))
end

function Polygon(exterior::AbstractVector{P}, faces::AbstractVector{<:LineFace}) where {P<:Point{Dim,T}} where {Dim,T}
    return Polygon(LineString(exterior, faces))
end

function Polygon(exterior::AbstractVector{P}, interior::AbstractVector{<:AbstractVector{P}}) where {P<:Point{Dim,T}} where {Dim,T}
    ext = LineString(exterior)
    # We need to take extra care for empty interiors, since
    # if we just map over it, it won't infer the element type correctly!
    int = typeof(ext)[]
    foreach(x -> push!(int, LineString(x)), interior)
    return Polygon(ext, int)
end

function coordinates(polygon::Polygon{N,T,PointType}) where {N,T,PointType}
    exterior = coordinates(polygon.exterior)
    if isempty(polygon.interiors)
        return exterior
    else
        result = PointType[]
        append!(result, exterior)
        foreach(x -> append!(result, coordinates(x)), polygon.interiors)
        return result
    end
end

"""
    MultiPolygon(polygons)
"""
struct MultiPolygon{Dim,T,Element<:Polytope{Dim,T},A<:AbstractVector{Element}} <: AbstractVector{Element}
    polygons::A
end

Base.getindex(mp::MultiPolygon, i) = mp.polygons[i]
Base.size(mp::MultiPolygon) = size(mp.polygons)

struct MultiLineString{Dim,T,Element<:LineString{Dim,T},A<:AbstractVector{Element}} <: AbstractVector{Element}
    linestrings::A
end

function MultiLineString(linestrings::AbstractVector{L}) where {L<:AbstractVector{LineP{Dim,T,P}}} where {Dim,T,P}
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

An abstract mesh is a collection of Polytope elements (Simplices / Ngons).
The connections are defined via faces(mesh), the coordinates of the elements are returned by
coordinates(mesh).
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
