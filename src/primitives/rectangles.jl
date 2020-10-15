"""
    Rectangle{N, T}

A `Rectangle` in `N` dimensions with coordinates of type `T`.
"""
struct Rectangle{N,T} <: GeometryPrimitive{N,T}
    origin::Point{N,T}
    widths::SVector{N,T}
end

origin(r::Rectangle) = r.origin
Base.minimum(r::Rectangle) = coordinates(r.origin)
Base.maximum(r::Rectangle) = coordinates(r.origin) + r.widths
Base.length(::Rectangle{N,T}) where {T,N} = N
widths(r::Rectangle) = r.widths
Base.:(==)(b1::Rectangle, b2::Rectangle) = minimum(b1) == minimum(b2) && widths(b1) == widths(b2)
Base.isequal(b1::Rectangle, b2::Rectangle) = b1 == b2

"""
    split(rectangle, axis, value)

Splits a rectangle into two along an axis at a given location.
"""
split(b::Rectangle, axis, value::Number) = _split(b, axis, value)
function _split(b::H, axis, value) where {H<:Rectangle}
    bmin = minimum(b)
    bmax = maximum(b)
    b1max = setindex(bmax, value, axis)
    b2min = setindex(bmin, value, axis)

    return H(bmin, b1max - bmin), H(b2min, bmax - b2min)
end

function minmax(p::StaticVector, vmin, vmax)
    any(isnan, p) && return (vmin, vmax)
    return min.(p, vmin), max.(p, vmax)
end

# set operations

"""
Perform a union between two Rects.
"""
function Base.union(h1::Rectangle{N,T1}, h2::Rectangle{N,T2}) where {N,T1,T2}
    T = promote_type(T1, T2)
    m = min.(minimum(h1), minimum(h2))
    mm = max.(maximum(h1), maximum(h2))
    return Rectangle{N,T}(Point(m), Vec(mm - m))
end

"""
    intersect(h1::Rectangle, h2::Rectangle)

Perform a intersection between two rectangles.
"""
function intersect(h1::Rectangle{N,T1}, h2::Rectangle{N,T2}) where {N,T1,T2}
    T = promote_type(T1, T2)
    m = max.(minimum(h1), minimum(h2))
    mm = min.(maximum(h1), maximum(h2))
    Rectangle{N,T}(Point(m), Vec(mm - m))
end

# http://en.wikipedia.org/wiki/Allen%27s_interval_algebra
function before(b1::Rectangle{N}, b2::Rectangle{N}) where {N}
    for i in 1:N
        maximum(b1)[i] < minimum(b2)[i] || return false
    end
    return true
end

function overlaps(b1::Rectangle{N}, b2::Rectangle{N}) where {N}
    for i in 1:N
        maximum(b2)[i] > maximum(b1)[i] > minimum(b2)[i] &&
            minimum(b1)[i] < minimum(b2)[i] || return false
    end
    return true
end

function starts(b1::Rectangle{N}, b2::Rectangle{N}) where {N}
    return if minimum(b1) == minimum(b2)
        for i in 1:N
            maximum(b1)[i] < maximum(b2)[i] || return false
        end
        return true
    else
        return false
    end
end

function during(b1::Rectangle{N}, b2::Rectangle{N}) where {N}
    for i in 1:N
        maximum(b1)[i] < maximum(b2)[i] && minimum(b1)[i] > minimum(b2)[i] || return false
    end
    return true
end

function finishes(b1::Rectangle{N}, b2::Rectangle{N}) where {N}
    return if maximum(b1) == maximum(b2)
        for i in 1:N
            minimum(b1)[i] > minimum(b2)[i] || return false
        end
        return true
    else
        return false
    end
end

# containment

"""
    in(b1::Rectangle, b2::Rectangle)

Check if Rectangle `b1` is contained in `b2`. This does not use
strict inequality, so Rects may share faces and this will still
return true.
"""
function Base.in(b1::Rectangle{N}, b2::Rectangle{N}) where {N}
    for i in 1:N
        maximum(b1)[i] <= maximum(b2)[i] && minimum(b1)[i] >= minimum(b2)[i] || return false
    end
    return true
end

"""
    in(pt::Point, b1::Rectangle{N, T})

Check if a point is contained in a Rectangle. This will return true if
the point is on a face of the Rectangle.
"""
function Base.in(pt::Point, b1::Rectangle{N,T}) where {T,N}
    cs = coordinates(pt)
    for i in 1:N
        cs[i] <= maximum(b1)[i] && cs[i] >= minimum(b1)[i] || return false
    end
    return true
end

# decomposition

function faces(::Rectangle{2,T}, nvertices=(2, 2)) where {T}
    w, h = nvertices
    idx = LinearIndices(nvertices)
    quad(i, j) = QuadFace{Int}(idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1])
    return ivec((quad(i, j) for i in 1:(w - 1), j in 1:(h - 1)))
end

function coordinates(rect::Rectangle{2,T}, nvertices=(2, 2)) where {T}
    mini, maxi = extrema(rect)
    xrange, yrange = LinRange.(mini, maxi, nvertices)
    return ivec(Vec(x, y) for x in xrange, y in yrange)
end

function texturecoordinates(::Rectangle{2,T}, nvertices=(2, 2)) where {T}
    xrange, yrange = LinRange.((0, 1), (1, 0), nvertices)
    return ivec(Vec(x, y) for x in xrange, y in yrange)
end

function normals(::Rectangle{2,T}, nvertices=(2, 2)) where {T}
    return Iterators.repeated(Vec(0, 0, 1), prod(nvertices))
end

##
# Rect3D decomposition
function coordinates(rect::Rectangle{3,T}) where {T}
    # TODO use n
    w  = widths(rect)
    o  = coordinates(origin(rect))
    xs = Vec{3,Int}[(0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0), (0, 0, 0), (1, 0, 0),
                    (1, 0, 1), (0, 0, 1), (0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0),
                    (1, 1, 1), (0, 1, 1), (0, 0, 1), (1, 0, 1), (1, 1, 1), (1, 0, 1),
                    (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 1, 0), (0, 1, 0), (0, 1, 1)]
    return ((x .* w .+ o) for x in xs)
end

function texturecoordinates(::Rectangle{3,T}) where {T}
    return coordinates(Rectangle(Point{3,T}(0,0,0), SVector{3,T}(1,1,1)))
end

function faces(::Rectangle{3,T}) where {T}
    return QuadFace{Int}[(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16),
                         (17, 18, 19, 20), (21, 22, 23, 24)]
end
