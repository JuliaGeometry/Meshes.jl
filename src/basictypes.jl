"""
Abstract Geometry in R{Dim} with Number type T
"""
abstract type AbstractGeometry{Dim,T<:Number} end
abstract type GeometryPrimitive{Dim,T} <: AbstractGeometry{Dim,T} end
Base.ndims(x::AbstractGeometry{Dim}) where {Dim} = Dim

"""
Geometry made of N connected points. Connected as one flat geometry, it makes a Ngon / Polygon.
Connected as volume it will be a Simplex / Tri / Cube.
Note That `Polytope{N} where N == 3` denotes a Triangle both as a Simplex or Ngon.
"""
abstract type Polytope{Dim,T} <: AbstractGeometry{Dim,T} end
abstract type AbstractPolygon{Dim,T} <: Polytope{Dim,T} end

abstract type AbstractFace{N,T} <: StaticVector{N,T} end
abstract type AbstractSimplexFace{N,T} <: AbstractFace{N,T} end
abstract type AbstractNgonFace{N,T} <: AbstractFace{N,T} end

abstract type AbstractSimplex{Dim,N,T} <: StaticVector{Dim,T} end

"""
Face index, connecting points to form a simplex
"""

@fixed_vector SimplexFace AbstractSimplexFace
const TetrahedronFace{T} = SimplexFace{4,T}
Face(::Type{<:SimplexFace{N}}, ::Type{T}) where {N,T} = SimplexFace{N,T}

"""
Face index, connecting points to form an Ngon
"""

@fixed_vector NgonFace AbstractNgonFace
const LineFace{T} = NgonFace{2,T}
const TriangleFace{T} = NgonFace{3,T}
const QuadFace{T} = NgonFace{4,T}

function Base.show(io::IO, x::TriangleFace{T}) where {T}
    return print(io, "TriangleFace(", join(x, ", "), ")")
end

Face(::Type{<:NgonFace{N}}, ::Type{T}) where {N,T} = NgonFace{N,T}
Face(F::Type{NgonFace{N,FT}}, ::Type{T}) where {FT,N,T} = F

@propagate_inbounds Base.getindex(x::Polytope, i::Integer) = coordinates(x)[i]
@propagate_inbounds Base.iterate(x::Polytope) = iterate(coordinates(x))
@propagate_inbounds Base.iterate(x::Polytope, i) = iterate(coordinates(x), i)

"""
Fixed Size Polygon, e.g.
- N 1-2 : Illegal!
- N = 3 : Triangle
- N = 4 : Quadrilateral (or Quad, Or tetragon)
- N = 5 : Pentagon
- ...
"""
struct Ngon{Dim,T<:Real,N,Point<:AbstractPoint{Dim,T}} <: AbstractPolygon{Dim,T}
    points::SVector{N,Point}
end

const NNgon{N} = Ngon{Dim,T,N,P} where {Dim,T,P}

function (::Type{<:NNgon{N}})(points::Vararg{P,N}) where {P<:AbstractPoint{Dim,T},
                                                          N} where {Dim,T}
    return Ngon{Dim,T,N,P}(SVector(points))
end
Base.show(io::IO, x::NNgon{N}) where {N} = print(io, "Ngon{$N}(", join(x, ", "), ")")

# Interfaces
coordinates(x::Ngon) = x.points
# Base Array interface
Base.length(::Type{<:NNgon{N}}) where {N} = N
Base.length(::NNgon{N}) where {N} = N

"""
The Ngon Polytope element type when indexing an array of points with a SimplexFace
"""
function Polytope(P::Type{<:AbstractPoint{Dim,T}},
                  ::Type{<:AbstractNgonFace{N,IT}}) where {N,Dim,T,IT}
    return Ngon{Dim,T,N,P}
end

"""
The fully concrete Ngon type, when constructed from a point type!
"""
function Polytope(::Type{<:NNgon{N}}, P::Type{<:AbstractPoint{NDim,T}}) where {N,NDim,T}
    return Ngon{NDim,T,N,P}
end

const LineP{Dim,T,P<:AbstractPoint{Dim,T}} = Ngon{Dim,T,2,P}
const Line{Dim,T} = LineP{Dim,T,Point{Dim,T}}

# Simplex{D, T, 3} & Ngon{D, T, 3} are both representing a triangle.
# Since Ngon is supposed to be flat and a triangle is flat, lets prefer Ngon
# for triangle:
const TriangleP{Dim,T,P<:AbstractPoint{Dim,T}} = Ngon{Dim,T,3,P}
const Triangle{Dim,T} = TriangleP{Dim,T,Point{Dim,T}}
const Triangle3d{T} = Triangle{3,T}

Base.show(io::IO, x::TriangleP) = print(io, "Triangle(", join(x, ", "), ")")
Base.summary(io::IO, x::Type{<:TriangleP}) = print(io, "Triangle")

const Quadrilateral{Dim,T} = Ngon{Dim,T,4,P} where {P<:AbstractPoint{Dim,T}}

Base.show(io::IO, x::Quadrilateral) = print(io, "Quad(", join(x, ", "), ")")
Base.summary(io::IO, x::Type{<:Quadrilateral}) = print(io, "Quad")

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
A `Simplex` is a generalization of an N-dimensional tetrahedra and can be thought
of as a minimal convex set containing the specified points.

* A 0-simplex is a point.
* A 1-simplex is a line segment.
* A 2-simplex is a triangle.
* A 3-simplex is a tetrahedron.

Note that this datatype is offset by one compared to the traditional
mathematical terminology. So a one-simplex is represented as `Simplex{2,T}`.
This is for a simpler implementation.

It applies to infinite dimensions. The structure of this type is designed
to allow embedding in higher-order spaces by parameterizing on `T`.
"""
struct Simplex{Dim,T<:Real,N,Point<:AbstractPoint{Dim,T}} <: Polytope{Dim,T}

    points::SVector{N,Point}
end

const NSimplex{N} = Simplex{Dim,T,N,P} where {Dim,T,P}
const TetrahedronP{T,P<:AbstractPoint{3,T}} = Simplex{3,T,4,P}
const Tetrahedron{T} = TetrahedronP{T,Point{3,T}}

Base.show(io::IO, x::TetrahedronP) = print(io, "Tetrahedron(", join(x, ", "), ")")

coordinates(x::Simplex) = x.points

function (::Type{<:NSimplex{N}})(points::Vararg{P,N}) where {P<:AbstractPoint{Dim,T},
                                                             N} where {Dim,T}
    return Simplex{Dim,T,N,P}(SVector(points))
end

# Base Array interface
Base.length(::Type{<:NSimplex{N}}) where {N} = N
Base.length(::NSimplex{N}) where {N} = N

"""
The Simplex Polytope element type when indexing an array of points with a SimplexFace
"""
function Polytope(P::Type{<:AbstractPoint{Dim,T}},
                  ::Type{<:AbstractSimplexFace{N}}) where {N,Dim,T}
    return Simplex{Dim,T,N,P}
end

"""
The fully concrete Simplex type, when constructed from a point type!
"""
function Polytope(::Type{<:NSimplex{N}}, P::Type{<:AbstractPoint{NDim,T}}) where {N,NDim,T}
    return Simplex{NDim,T,N,P}
end
Base.show(io::IO, x::LineP) = print(io, "Line(", x[1], " => ", x[2], ")")

"""
    LineString(points::AbstractVector{<:AbstractPoint})

A LineString is a geometry of connected line segments
"""
struct LineString{Dim,T<:Real,P<:AbstractPoint,V<:AbstractVector{<:LineP{Dim,T,P}}} <:
       AbstractVector{LineP{Dim,T,P}}
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
    LineString(points::AbstractVector{<: AbstractPoint}, skip = 1)

Creates a LineString from a vector of points
With `skip == 1`, the default, it will connect the line like this:
```julia
points = Point[a, b, c, d]
linestring = LineString(points)
@assert linestring == LineString([a => b, b => c, c => d])
```
"""
function LineString(points::AbstractVector{<:AbstractPoint}, skip=1)
    return LineString(connect(points, LineP, skip))
end

function LineString(points::AbstractVector{<:Pair{P,P}}) where {P<:AbstractPoint{N,T}} where {N,
                                                                                              T}
    return LineString(reinterpret(LineP{N,T,P}, points))
end

function LineString(points::AbstractVector{<:AbstractPoint},
                    faces::AbstractVector{<:LineFace})
    return LineString(connect(points, faces))
end

"""
    LineString(points::AbstractVector{<: AbstractPoint}, indices::AbstractVector{<: Integer}, skip = 1)

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
function LineString(points::AbstractVector{<:AbstractPoint},
                    indices::AbstractVector{<:Integer}, skip=1)
    faces = connect(indices, LineFace, skip)
    return LineString(points, faces)
end

"""
    Polygon(exterior::AbstractVector{<:Point})
    Polygon(exterior::AbstractVector{<:Point}, interiors::Vector{<:AbstractVector{<:AbstractPoint}})

"""
struct Polygon{Dim,T<:Real,P<:AbstractPoint{Dim,T},L<:AbstractVector{<:LineP{Dim,T,P}},
               V<:AbstractVector{L}} <: AbstractPolygon{Dim,T}
    exterior::L
    interiors::V
end

Base.copy(x::Polygon) = Polygon(copy(x.exterior), copy(x.interiors))

function Base.:(==)(a::Polygon, b::Polygon)
    return (a.exterior == b.exterior) && (a.interiors == b.interiors)
end

function Polygon(exterior::E,
                 interiors::AbstractVector{E}) where {E<:AbstractVector{LineP{Dim,T,P}}} where {Dim,
                                                                                                T,
                                                                                                P}
    return Polygon{Dim,T,P,typeof(exterior),typeof(interiors)}(exterior, interiors)
end

Polygon(exterior::L) where {L<:AbstractVector{<:LineP}} = Polygon(exterior, L[])

function Polygon(exterior::AbstractVector{P},
                 skip::Int=1) where {P<:AbstractPoint{Dim,T}} where {Dim,T}
    return Polygon(LineString(exterior, skip))
end

function Polygon(exterior::AbstractVector{P}, faces::AbstractVector{<:Integer},
                 skip::Int=1) where {P<:AbstractPoint{Dim,T}} where {Dim,T}
    return Polygon(LineString(exterior, faces, skip))
end

function Polygon(exterior::AbstractVector{P},
                 faces::AbstractVector{<:LineFace}) where {P<:AbstractPoint{Dim,T}} where {Dim,
                                                                                           T}
    return Polygon(LineString(exterior, faces))
end

function Polygon(exterior::AbstractVector{P},
                 interior::AbstractVector{<:AbstractVector{P}}) where {P<:AbstractPoint{Dim,
                                                                                        T}} where {Dim,
                                                                                                   T}
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
    MultiPolygon(polygons::AbstractPolygon)
"""
struct MultiPolygon{Dim,T<:Real,Element<:AbstractPolygon{Dim,T},
                    A<:AbstractVector{Element}} <: AbstractVector{Element}
    polygons::A
end

function MultiPolygon(polygons::AbstractVector{P};
                      kw...) where {P<:AbstractPolygon{Dim,T}} where {Dim,T}
    return MultiPolygon(meta(polygons; kw...))
end

Base.getindex(mp::MultiPolygon, i) = mp.polygons[i]
Base.size(mp::MultiPolygon) = size(mp.polygons)

struct MultiLineString{Dim,T<:Real,Element<:LineString{Dim,T},A<:AbstractVector{Element}} <:
       AbstractVector{Element}

    linestrings::A
end

function MultiLineString(linestrings::AbstractVector{L};
                         kw...) where {L<:AbstractVector{LineP{Dim,T,P}}} where {Dim,T,P}
    return MultiLineString(meta(linestrings; kw...))
end

Base.getindex(ms::MultiLineString, i) = ms.linestrings[i]
Base.size(ms::MultiLineString) = size(ms.linestrings)

"""
    MultiPoint(points::AbstractVector{AbstractPoint})

A collection of points
"""
struct MultiPoint{Dim,T<:Real,P<:AbstractPoint{Dim,T},A<:AbstractVector{P}} <:
       AbstractVector{P}

    points::A
end

function MultiPoint(points::AbstractVector{P};
                    kw...) where {P<:AbstractPoint{Dim,T}} where {Dim,T}
    return MultiPoint(meta(points; kw...))
end

Base.getindex(mpt::MultiPoint, i) = mpt.points[i]
Base.size(mpt::MultiPoint) = size(mpt.points)

"""
    AbstractMesh

An abstract mesh is a collection of Polytope elements (Simplices / Ngons).
The connections are defined via faces(mesh), the coordinates of the elements are returned by
coordinates(mesh). Arbitrary meta information can be attached per point or per face
"""
abstract type AbstractMesh{Element<:Polytope} <: AbstractVector{Element} end

"""
    Mesh <: AbstractVector{Element}
The conrecte AbstractMesh implementation
"""
struct Mesh{Dim,T<:Number,Element<:Polytope{Dim,T},V<:AbstractVector{Element}} <:
       AbstractMesh{Element}
    simplices::V # usually a FaceView, to connect a set of points via a set of faces.
end

Tables.schema(mesh::Mesh) = Tables.schema(getfield(mesh, :simplices))

function Base.getproperty(mesh::Mesh, name::Symbol)
    if name === :position
        # a mesh always has position defined by coordinates...
        return metafree(coordinates(mesh))
    else
        return getproperty(getfield(mesh, :simplices), name)
    end
end

function Base.propertynames(mesh::Mesh)
    names = propertynames(getfield(mesh, :simplices))
    if :position in names
        return names
    else
        # a mesh always has positions!
        return (names..., :position)
    end
end

function Base.summary(io::IO, ::Mesh{Dim,T,Element}) where {Dim,T,Element}
    print(io, "Mesh{$Dim, $T, ")
    summary(io, Element)
    return print(io, "}")
end

Base.size(mesh::Mesh) = size(getfield(mesh, :simplices))
Base.getindex(mesh::Mesh, i::Integer) = getfield(mesh, :simplices)[i]

function Mesh(elements::AbstractVector{<:Polytope{Dim,T}}) where {Dim,T}
    return Mesh{Dim,T,eltype(elements),typeof(elements)}(elements)
end

function Mesh(points::AbstractVector{<:AbstractPoint},
              faces::AbstractVector{<:AbstractFace})
    return Mesh(connect(points, faces))
end

function Mesh(points::AbstractVector{<:AbstractPoint}, faces::AbstractVector{<:Integer},
              facetype=TriangleFace, skip=1)
    return Mesh(connect(points, connect(faces, facetype, skip)))
end
