# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

include("faces.jl")

"""
    TupleView{T,N,M,A}

TupleView, groups elements of an array as tuples.
N is the dimension of the tuple, M tells how many elements to skip to the next tuple.
By default TupleView{N} defaults to skip N items.

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
struct TupleView{T,N,M,A} <: AbstractVector{T}
    data::A
    connect::Bool
end

TupleView{N,M}(x::AbstractVector{T}; connect=false) where {T,N,M} =
    TupleView{NTuple{N,T},N,M,typeof(x)}(x, connect)

TupleView{N}(x::AbstractVector; connect=false) where {N} =
    TupleView{N,N}(x, connect=connect)

function Base.size(x::TupleView{T,N,M}) where {T,N,M}
    nitems = length(x.data) ÷ (N - (N - M))
    nitems = nitems - max(N - M, 0)
    (nitems + x.connect,) # plus one item if we connect
end

Base.getindex(x::TupleView{T,N,M}, index::Integer) where {T,N,M} =
    ntuple(i -> x.data[mod1(((index - 1) * M) + i, length(x.data))], N)

"""
    connect(points::AbstractVector{<:Point}, ::Type{<:Polytope}, skip::Int)

Creates a view that connects a number of points to a Polytope `P`.
Between each polytope, `skip` elements are skipped until the next starts.
Example:
```julia
x = connect(Point[(1, 2), (3, 4), (5, 6), (7, 8)], Segment, 2)
x == [Segment(Point(1, 2), Point(3, 4)), Segment(Point(5, 6), Point(7, 8))]
```
"""
@inline function connect(points::AbstractVector{<:Point},
                         PL::Type{<:Polytope}, skip::Int=length(PL))
    return reinterpret(PL, TupleView{length(PL),skip}(points))
end

@inline function connect(points::AbstractVector{T}, P::Type{<:AbstractFace{N}},
                         skip::Int=N) where {T <: Real,N}
    return reinterpret(Face(P, T), TupleView{N,skip}(points))
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
```
"""
struct FaceView{Element,PT<:Point,Face<:AbstractFace,P<:AbstractVector{PT},F<:AbstractVector{Face}}<:AbstractVector{Element}
    elements::P
    faces::F
end

coordinates(mesh::FaceView) = getfield(mesh, :elements)
faces(mesh::FaceView) = getfield(mesh, :faces)

Base.size(faceview::FaceView) = size(getfield(faceview, :faces))

@propagate_inbounds function Base.getindex(x::FaceView{Element}, i) where {Element}
    return Element(map(idx -> coordinates(x)[idx], faces(x)[i]))
end

@propagate_inbounds function Base.setindex!(x::FaceView{Element}, element::Element, i) where {Element}
    face = faces(x)[i]
    for (i, f) in enumerate(face) # TODO unroll!?
        coordinates(x)[face[i]] = element[i]
    end
    return element
end

function connect(points::AbstractVector{P}, faces::AbstractVector{F}) where {Dim,T,N,V,P<:Point{Dim,T},F<:AbstractFace{N,V}}
    Element = N == 3 ? Triangle{Dim,T} : Tetrahedron{Dim,T}
    return FaceView{Element,P,F,typeof(points),typeof(faces)}(points, faces)
end
