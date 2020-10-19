# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Box(min, max)

An axis-aligned box with `min` and `max` corners.
See https://en.wikipedia.org/wiki/Hyperrectangle.

## Example

```julia
Box(Point(0,0,0), Point(1,1,1)) # unit cube
```
"""
struct Box{N,T} <: Primitive{N,T}
    min::Point{N,T}
    max::Point{N,T}
end

Base.minimum(b::Box) = b.min
Base.maximum(b::Box) = b.max
sides(b::Box) = b.max - b.min
volume(b::Box) = prod(b.max - b.min)

function minmax(p::StaticVector, vmin, vmax)
    any(isnan, p) && return (vmin, vmax)
    return min.(p, vmin), max.(p, vmax)
end

# set operations

"""
    union(b1::Box, b2::Box)

Union between boxes.
"""
function Base.union(b1::Box{N,T1}, b2::Box{N,T2}) where {N,T1,T2}
    m1, M1 = minimum(b1), maximum(b1)
    m2, M2 = minimum(b2), maximum(b2)
    T = promote_type(T1, T2)
    m = min.(coordinates(m1), coordinates(m2))
    M = max.(coordinates(M1), coordinates(M2))
    return Box{N,T}(Point(m), Point(M))
end

"""
    intersect(b1::Box, b2::Box)

Intersection between boxes.
"""
function intersect(b1::Box{N,T1}, b2::Box{N,T2}) where {N,T1,T2}
    m1, M1 = minimum(b1), maximum(b1)
    m2, M2 = minimum(b2), maximum(b2)
    T = promote_type(T1, T2)
    m = max.(coordinates(m1), coordinates(m2))
    M = min.(coordinates(M1), coordinates(M2))
    return Box{N,T}(Point(m), Point(M))
end

function before(b1::Box{N}, b2::Box{N}) where {N}
    M1 = coordinates(maximum(b1))
    m2 = coordinates(minimum(b2))
    for i in 1:N
        M1[i] < m2[i] || return false
    end
    return true
end

function overlaps(b1::Box{N}, b2::Box{N}) where {N}
    m1 = coordinates(minimum(b1))
    M1 = coordinates(maximum(b1))
    m2 = coordinates(minimum(b2))
    M2 = coordinates(maximum(b2))
    for i in 1:N
        M2[i] > M1[i] > m2[i] && m1[i] < m2[i] || return false
    end
    return true
end

function starts(b1::Box{N}, b2::Box{N}) where {N}
    m1 = coordinates(minimum(b1))
    M1 = coordinates(maximum(b1))
    m2 = coordinates(minimum(b2))
    M2 = coordinates(maximum(b2))
    return if m1 == m2
        for i in 1:N
            M1[i] < M2[i] || return false
        end
        return true
    else
        return false
    end
end

function during(b1::Box{N}, b2::Box{N}) where {N}
    m1 = coordinates(minimum(b1))
    M1 = coordinates(maximum(b1))
    m2 = coordinates(minimum(b2))
    M2 = coordinates(maximum(b2))
    for i in 1:N
        M1[i] < M2[i] && m1[i] > m2[i] || return false
    end
    return true
end

function finishes(b1::Box{N}, b2::Box{N}) where {N}
    m1 = coordinates(minimum(b1))
    M1 = coordinates(maximum(b1))
    m2 = coordinates(minimum(b2))
    M2 = coordinates(maximum(b2))
    return if M1 == M2
        for i in 1:N
            m1[i] > m2[i] || return false
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
    m1 = coordinates(minimum(b1))
    M1 = coordinates(maximum(b1))
    m2 = coordinates(minimum(b2))
    M2 = coordinates(maximum(b2))
    for i in 1:N
        M1[i] ≤ M2[i] && m1[i] ≥ m2[i] || return false
    end
    return true
end

"""
    in(p::Point, b::Box)

Check if a point is in the box.
"""
function Base.in(p::Point{N,T}, b::Box{N,T}) where {T,N}
    m = coordinates(minimum(b))
    M = coordinates(maximum(b))
    x = coordinates(p)
    for i in 1:N
        m[i] ≤ x[i] && x[i] ≤ M[i] || return false
    end
    return true
end

# decomposition

function coordinates(b::Box{2,T}, nvertices=(2, 2)) where {T}
    m = coordinates(minimum(b))
    M = coordinates(maximum(b))
    xrange, yrange = LinRange.(m, M, nvertices)
    return ivec(Vec(x, y) for x in xrange, y in yrange)
end

function coordinates(b::Box{3,T}) where {T}
    # TODO use n
    w  = sides(b)
    m  = coordinates(minimum(b))
    xs = Vec{3,Int}[(0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0), (0, 0, 0), (1, 0, 0),
                    (1, 0, 1), (0, 0, 1), (0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0),
                    (1, 1, 1), (0, 1, 1), (0, 0, 1), (1, 0, 1), (1, 1, 1), (1, 0, 1),
                    (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 1, 0), (0, 1, 0), (0, 1, 1)]
    return ((x .* w .+ m) for x in xs)
end

function faces(::Box{2,T}, nvertices=(2, 2)) where {T}
    w, h = nvertices
    idx = LinearIndices(nvertices)
    quad(i, j) = QuadFace((idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]))
    return ivec((quad(i, j) for i in 1:(w - 1), j in 1:(h - 1)))
end

function faces(::Box{3,T}) where {T}
    return QuadFace.([(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16),
                      (17, 18, 19, 20), (21, 22, 23, 24)])
end
