# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Box(origin, widths)

An axis-aligned box with `origin` and `widths`.
"""
struct Box{N,T} <: Primitive{N,T}
    origin::Point{N,T}
    widths::SVector{N,T}
end

origin(b::Box) = b.origin
Base.minimum(r::Box) = coordinates(r.origin)
Base.maximum(r::Box) = coordinates(r.origin) + r.widths
Base.length(::Box{N,T}) where {T,N} = N
widths(r::Box) = r.widths
Base.:(==)(b1::Box, b2::Box) = minimum(b1) == minimum(b2) && widths(b1) == widths(b2)
Base.isequal(b1::Box, b2::Box) = b1 == b2

"""
    split(Box, axis, value)

Splits a Box into two along an axis at a given location.
"""
split(b::Box, axis, value::Number) = _split(b, axis, value)
function _split(b::H, axis, value) where {H<:Box}
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
function Base.union(h1::Box{N,T1}, h2::Box{N,T2}) where {N,T1,T2}
    T = promote_type(T1, T2)
    m = min.(minimum(h1), minimum(h2))
    mm = max.(maximum(h1), maximum(h2))
    return Box{N,T}(Point(m), Vec(mm - m))
end

"""
    intersect(h1::Box, h2::Box)

Perform a intersection between two rectangles.
"""
function intersect(h1::Box{N,T1}, h2::Box{N,T2}) where {N,T1,T2}
    T = promote_type(T1, T2)
    m = max.(minimum(h1), minimum(h2))
    mm = min.(maximum(h1), maximum(h2))
    Box{N,T}(Point(m), Vec(mm - m))
end

# http://en.wikipedia.org/wiki/Allen%27s_interval_algebra
function before(b1::Box{N}, b2::Box{N}) where {N}
    for i in 1:N
        maximum(b1)[i] < minimum(b2)[i] || return false
    end
    return true
end

function overlaps(b1::Box{N}, b2::Box{N}) where {N}
    for i in 1:N
        maximum(b2)[i] > maximum(b1)[i] > minimum(b2)[i] &&
            minimum(b1)[i] < minimum(b2)[i] || return false
    end
    return true
end

function starts(b1::Box{N}, b2::Box{N}) where {N}
    return if minimum(b1) == minimum(b2)
        for i in 1:N
            maximum(b1)[i] < maximum(b2)[i] || return false
        end
        return true
    else
        return false
    end
end

function during(b1::Box{N}, b2::Box{N}) where {N}
    for i in 1:N
        maximum(b1)[i] < maximum(b2)[i] && minimum(b1)[i] > minimum(b2)[i] || return false
    end
    return true
end

function finishes(b1::Box{N}, b2::Box{N}) where {N}
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
    in(b1::Box, b2::Box)

Check if Box `b1` is contained in `b2`. This does not use
strict inequality, so Rects may share faces and this will still
return true.
"""
function Base.in(b1::Box{N}, b2::Box{N}) where {N}
    for i in 1:N
        maximum(b1)[i] <= maximum(b2)[i] && minimum(b1)[i] >= minimum(b2)[i] || return false
    end
    return true
end

"""
    in(pt::Point, b1::Box{N, T})

Check if a point is contained in a Box. This will return true if
the point is on a face of the Box.
"""
function Base.in(pt::Point, b1::Box{N,T}) where {T,N}
    cs = coordinates(pt)
    for i in 1:N
        cs[i] <= maximum(b1)[i] && cs[i] >= minimum(b1)[i] || return false
    end
    return true
end

# decomposition

function faces(::Box{2,T}, nvertices=(2, 2)) where {T}
    w, h = nvertices
    idx = LinearIndices(nvertices)
    quad(i, j) = QuadFace((idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]))
    return ivec((quad(i, j) for i in 1:(w - 1), j in 1:(h - 1)))
end

function coordinates(rect::Box{2,T}, nvertices=(2, 2)) where {T}
    mini, maxi = extrema(rect)
    xrange, yrange = LinRange.(mini, maxi, nvertices)
    return ivec(Vec(x, y) for x in xrange, y in yrange)
end

##
# Rect3D decomposition
function coordinates(rect::Box{3,T}) where {T}
    # TODO use n
    w  = widths(rect)
    o  = coordinates(origin(rect))
    xs = Vec{3,Int}[(0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0), (0, 0, 0), (1, 0, 0),
                    (1, 0, 1), (0, 0, 1), (0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0),
                    (1, 1, 1), (0, 1, 1), (0, 0, 1), (1, 0, 1), (1, 1, 1), (1, 0, 1),
                    (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 1, 0), (0, 1, 0), (0, 1, 1)]
    return ((x .* w .+ o) for x in xs)
end

function faces(::Box{3,T}) where {T}
    return QuadFace.([(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16),
                      (17, 18, 19, 20), (21, 22, 23, 24)])
end
