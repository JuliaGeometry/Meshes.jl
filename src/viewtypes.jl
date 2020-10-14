"""
    TupleView{T, N, Skip, A}

TupleView, groups elements of an array as tuples.
N is the dimension of the tuple, M tells how many elements to skip to the next tuple.
By default TupleView{N} defaults to skip N items.
# a few examples:

```julia

x = [1, 2, 3, 4, 5, 6]
TupleView{2, 1}(x):
> [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]

TupleView{2}(x):
> [(1, 2), (3, 4), (5, 6)]

TupleView{2, 3}(x):
> [(1, 2), (4, 5)]

TupleView{3, 1}(x):
> [(1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, 6)]
```

TupleView can be used together with reinterpret:

```julia
x = [1, 2, 3, 4]
y = reinterpret(Point{2, Int}, TupleView{2, 1}(x))
> [Point(1, 2), Point(2, 3), Point(3, 4)]
```

"""
struct TupleView{T,N,Skip,A} <: AbstractVector{T}
    data::A
    connect::Bool
end

coordinates(tw::TupleView) = coordinates(tw.data)

function Base.size(x::TupleView{T,N,M}) where {T,N,M}
    nitems = length(x.data) ÷ (N - (N - M))
    nitems = nitems - max(N - M, 0)
    return (nitems + x.connect,) # plus one item if we connect
end

function Base.getindex(x::TupleView{T,N,M}, index::Integer) where {T,N,M}
    return ntuple(i -> x.data[mod1(((index - 1) * M) + i, length(x.data))], N)
end

function TupleView{N}(x::AbstractVector; connect=false) where {N}
    return TupleView{N,N}(x, connect=connect)
end

function TupleView{N,M}(x::AbstractVector{T}; connect=false) where {T,N,M}
    return TupleView{NTuple{N,T},N,M,typeof(x)}(x, connect)
end

@inline function connected_line(points::AbstractVector{<:AbstractPoint{N}},
                                skip=N) where {N}
    return connect(points, Line, skip)
end

"""
    connect(points::AbstractVector{<: AbstractPoint}, P::Type{<: Polytype{N}}, skip::Int = N)

Creates a view that connects a number of points to a Polytope `P`.
Between each polytope, `skip` elements are skipped untill the next starts.
Example:
```julia
x = connect(Point[(1, 2), (3, 4), (5, 6), (7, 8)], Line, 2)
x == [Line(Point(1, 2), Point(3, 4)), Line(Point(5, 6), Point(7, 8))]
"""
@inline function connect(points::AbstractVector{Point},
                         P::Type{<:Polytope{N,T} where {N,T}},
                         skip::Int=length(P)) where {Point <: AbstractPoint}
    return reinterpret(Polytope(P, Point), TupleView{length(P),skip}(points))
end

@inline function connect(points::AbstractVector{T}, ::Type{<:AbstractPoint{N}},
                         skip::Int=N) where {T <: Real,N}
    return reinterpret(Point{N,T}, TupleView{N,skip}(points))
end

@inline function connect(points::AbstractVector{T}, P::Type{<:AbstractFace{N}},
                         skip::Int=N) where {T <: Real,N}
    return reinterpret(Face(P, T), TupleView{N,skip}(points))
end

@inline function connect(points::AbstractMatrix{T},
                         P::Type{<:AbstractPoint{N}}) where {T <: Real,N}
    return if size(points, 1) === N
        return reinterpret(Point{N,T}, points)
    elseif size(points, 2) === N
        seglen = size(points, 1)
        columns = ntuple(N) do i
            return view(points, ((i - 1) * seglen + 1):(i * seglen))
        end
        return StructArray{Point{N,T}}((StructArray{NTuple{N,T}}(columns),))
    else
        error("Dim 1 or 2 must be equal to the point dimension!")
    end
end

"""
    FaceView{Elemnt, Point, Face, P, F}

FaceView enables to link one array of points via a face array, to generate one
abstract array of elements.
E.g., this becomes possible:
```
x = FaceView(rand(Point3f, 10), TriangleFace[(1, 2, 3), (2, 4, 5), ...])
x[1] isa Triangle == true
x isa AbstractVector{<: Triangle} == true
# This means we can use it as a mesh:
Mesh(x) # should just work!
Can also be used for Points:

linestring = FaceView(points, LineFace[...])
Polygon(linestring)
```
"""
struct FaceView{Element,Point <: AbstractPoint,Face <: AbstractFace,P <: AbstractVector{Point},F <: AbstractVector{Face}} <: AbstractVector{Element}

    elements::P
    faces::F
end

const SimpleFaceView{Dim,T,NFace,IndexType,PointType <: AbstractPoint{Dim,T},FaceType <: AbstractFace{NFace,IndexType}} = FaceView{Ngon{Dim,T,NFace,PointType},PointType,FaceType,Vector{PointType},Vector{FaceType}}

function Base.getproperty(faceview::FaceView, name::Symbol)
    return getproperty(getfield(faceview, :elements), name)
end

function Base.propertynames(faceview::FaceView)
    return propertynames(getfield(faceview, :elements))
end

Tables.schema(faceview::FaceView) = Tables.schema(getfield(faceview, :elements))

Base.size(faceview::FaceView) = size(getfield(faceview, :faces))

function Base.show(io::IO, ::Type{<:FaceView{Element}}) where {Element}
    if @isdefined Element
        print(io, "FaceView{", Element, "}")
    else
        print(io, "FaceView{T}")
    end
    return
end

@propagate_inbounds function Base.getindex(x::FaceView{Element}, i) where {Element}
    return Element(map(idx -> coordinates(x)[idx], faces(x)[i]))
end

@propagate_inbounds function Base.setindex!(x::FaceView{Element}, element::Element,
                                            i) where {Element}
    face = faces(x)[i]
    for (i, f) in enumerate(face) # TODO unroll!?
        coordinates(x)[face[i]] = element[i]
    end
    return element
end

function connect(points::AbstractVector{P},
                 faces::AbstractVector{F}) where {P <: AbstractPoint,F <: AbstractFace}
    return FaceView{Polytope(P, F),P,F,typeof(points),typeof(faces)}(points, faces)
end

coordinates(mesh::FaceView) = getfield(mesh, :elements)
faces(mesh::FaceView) = getfield(mesh, :faces)
