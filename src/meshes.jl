# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct Mesh{Dim,T,E<:Polytope{Dim,T},V<:AbstractVector{E}} <: Geometry{Dim,T}
    elements::V # usually a FaceView, to connect a set of points via a set of faces.
end

Mesh(points::AbstractVector{<:Point},
     faces::AbstractVector{<:AbstractFace}) = Mesh(connect(points, faces))

Mesh(points::AbstractVector{<:Point},
     faces::AbstractVector{<:Integer},
     facetype=TriangleFace, skip=1) = Mesh(connect(points, connect(faces, facetype, skip)))

Base.getindex(m::Mesh, i) = getindex(m.elements, i)
Base.length(m::Mesh) = length(m.elements)

Base.iterate(m::Mesh, i) = iterate(m.elements, i)
Base.iterate(m::Mesh) = iterate(m.elements)

elements(m::Mesh) = m.elements

coordinates(m::Mesh) = coordinates(m.elements)

volume(m::Mesh) = sum(volume, m)

# TOD: review these

convert_simplex(::Type{T}, x::T) where {T} = (x,)

"""
    convert_simplex(::Type{Face{3}}, f::Face{N})

Triangulate an N-Face into a tuple of triangular faces.
"""
@generated function convert_simplex(::Type{TriangleFace}, f::NgonFace{N}) where {N}
    3 <= N || error("decompose not implemented for N <= 3 yet. N: $N")
    v = Expr(:tuple)
    for i in 3:N
        push!(v.args, :(TriangleFace((f[1], f[$(i - 1)], f[$i]))))
    end
    return v
end

"""
    convert_simplex(::Type{Face{2}}, f::Face{N})

Extract all line segments in a Face.
"""
@generated function convert_simplex(::Type{LineFace}, f::NgonFace{N}) where {N}
    2 <= N || error("decompose not implented for N <= 2 yet. N: $N")

    v = Expr(:tuple)
    for i in 1:(N - 1)
        push!(v.args, :(LineFace((f[$i], f[$(i + 1)]))))
    end
    # connect vertices N and 1
    push!(v.args, :(LineFace((f[$N], f[1]))))
    return v
end

function convert_simplex(::Type{Point{N,T}}, p::Point{M,V}) where {N,T,M,V}
    x = coordinates(p)
    return (Point(ntuple(i -> i <= M ? T(x[i]) : T(0), N)),)
end

function convert_simplex(::Type{Vec{N,T}}, v::Vec{M,V}) where {N,T,M,V}
    return (Vec(ntuple(i -> i <= M ? T(v[i]) : T(0), N)),)
end

function decompose(::Type{T}, primitive) where {T}
    return collect_with_eltype(T, primitive)
end

function decompose(::Type{P}, primitive) where {P<:Point}
    return convert.(P, coordinates(primitive))
end

function decompose(::Type{F}, primitive) where {F<:AbstractFace}
    f = faces(primitive)
    f === nothing && return nothing
    return collect_with_eltype(F, f)
end

function collect_with_eltype(::Type{T}, iter) where {T}
    result = T[]
    for element in iter
        for telement in convert_simplex(T, element)
            push!(result, telement)
        end
    end
    return result
end

"""
    mesh(geometry::Meshable{N,T}; facetype=TriangleFace)

Creates a mesh from `geometry`.
"""
function mesh(geometry::Geometry{Dim,T}; facetype=TriangleFace) where {Dim,T}
    points = decompose(Point{Dim,T}, geometry)
    faces  = decompose(facetype, geometry)
    if faces === nothing
        # try to triangulate
        faces = decompose(facetype, points)
    end
    return Mesh(points, faces)
end
