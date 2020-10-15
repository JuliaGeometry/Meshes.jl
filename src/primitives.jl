##
# Generic base overloads
Base.extrema(primitive::GeometryPrimitive) = (minimum(primitive), maximum(primitive))
function widths(x::AbstractRange)
    mini, maxi = Float32.(extrema(x))
    return maxi - mini
end

##
# conversion & decompose
convert_simplex(::Type{T}, x::T) where {T} = (x,)

"""
    convert_simplex(::Type{Face{3}}, f::Face{N})

Triangulate an N-Face into a tuple of triangular faces.
"""
@generated function convert_simplex(::Type{TriangleFace}, f::Union{SimplexFace{N},NgonFace{N}}) where {N}
    3 <= N || error("decompose not implemented for N <= 3 yet. N: $N")
    v = Expr(:tuple)
    for i in 3:N
        push!(v.args, :(TriangleFace(f[1], f[$(i - 1)], f[$i])))
    end
    return v
end

"""
    convert_simplex(::Type{Face{2}}, f::Face{N})

Extract all line segments in a Face.
"""
@generated function convert_simplex(::Type{LineFace}, f::Union{SimplexFace{N},NgonFace{N}}) where {N}
    2 <= N || error("decompose not implented for N <= 2 yet. N: $N")

    v = Expr(:tuple)
    for i in 1:(N - 1)
        push!(v.args, :(LineFace(f[$i], f[$(i + 1)])))
    end
    # connect vertices N and 1
    push!(v.args, :(LineFace(f[$N], f[1])))
    return v
end

# TODO: review these
to_pointn(::Type{T}, x) where {T<:Point} = convert_simplex(T, x)[1]

# TODO: why increase the dimension of the point?
function convert_simplex(::Type{Point{N,T}}, p::Point{M,V}) where {N,T,M,V}
    x = coordinates(p)
    return (Point(ntuple(i -> i <= M ? T(x[i]) : T(0), N)),)
end

# TODO: review these
function convert_simplex(::Type{Vec{N,T}}, v::Vec{M,V}) where {N,T,M,V}
    return (Vec(ntuple(i -> i <= M ? T(v[i]) : T(0), N)),)
end

"""
The unnormalized normal of three vertices.
"""
function orthogonal_vector(v1, v2, v3)
    a = v2 - v1
    b = v3 - v1
    return cross(a, b)
end
