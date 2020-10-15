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

"""
    split(rectangle, axis, value)

Splits a rectangle into two along an axis at a given location.
"""
split(b::Rectangle, axis, value::Integer) = _split(b, axis, value)
split(b::Rectangle, axis, value::Number) = _split(b, axis, value)
function _split(b::H, axis, value) where {H<:Rectangle}
    bmin = minimum(b)
    bmax = maximum(b)
    b1max = setindex(bmax, value, axis)
    b2min = setindex(bmin, value, axis)

    return H(bmin, b1max - bmin), H(b2min, bmax - b2min)
end

###
# Algebraic operations

"""
    *(m::Mat, h::Rectangle)

Transform a `Rectangle` using a matrix. Maintains axis-align properties
so a significantly larger Rectangle may be generated.
"""
function Base.:(*)(m::Mat{N1,N1,T1}, h::Rectangle{N2,T2}) where {N1,N2,T1,T2}

    # TypeVar constants
    T = promote_type(T1, T2)
    D = N1 - N2

    # get all points on the Rectangle
    d = decompose(Point, h)
    # make sure our points are sized for the tranform
    pts = (Vec{N1,T}[vcat(pt, ones(Vec{D,T})) for pt in d]...,)::NTuple{2^N2,Vec{N1,T}}

    vmin = Vec{N1,T}(typemax(T))
    vmax = Vec{N1,T}(typemin(T))
    # tranform all points, tracking min and max points
    for pt in pts
        pn = m * pt
        vmin = min.(pn, vmin)
        vmax = max.(pn, vmax)
    end
    return Rectangle{N2,T}(vmin, vmax - vmin)
end

# equal dimension case
function Base.:(*)(m::Mat{N,N,T1}, h::Rectangle{N,T2}) where {N,T1,T2}

    # TypeVar constants
    T = promote_type(T1, T2)

    # get all points on the Rectangle
    pts = decompose(Point, h)

    # make sure our points are sized for the tranform
    vmin = Vec{N,T}(typemax(T))
    vmax = Vec{N,T}(typemin(T))

    # tranform all points, tracking min and max points
    for pt in pts
        pn = m * Vec(pt)
        vmin = min.(pn, vmin)
        vmax = max.(pn, vmax)
    end
    return Rectangle{N,T}(vmin, vmax - vmin)
end

# fast path. TODO make other versions fast without code duplications like now
function Base.:(*)(m::Mat{4,4,T}, h::Rectangle{3,T}) where {T}
    # equal dimension case

    # get all points on the Rectangle
    pts = (Vec{4,T}(0.0, 0.0, 0.0, 1.0), Vec{4,T}(1.0, 0.0, 0.0, 1.0),
           Vec{4,T}(0.0, 1.0, 0.0, 1.0), Vec{4,T}(1.0, 1.0, 0.0, 1.0),
           Vec{4,T}(0.0, 0.0, 1.0, 1.0), Vec{4,T}(1.0, 0.0, 1.0, 1.0),
           Vec{4,T}(0.0, 1.0, 1.0, 1.0), Vec{4,T}(1.0, 1.0, 1.0, 1.0))

    # make sure our points are sized for the tranform
    vmin = Vec{4,T}(typemax(T))
    vmax = Vec{4,T}(typemin(T))
    o, w = origin(h), widths(h)
    _o = Vec{4,T}(o[1], o[2], o[3], T(0))
    _w = Vec{4,T}(w[1], w[2], w[3], T(1))
    # tranform all points, tracking min and max points
    for pt in pts
        pn = m * (_o + (pt .* _w))
        vmin = min.(pn, vmin)
        vmax = max.(pn, vmax)
    end
    _vmin = Vec{3,T}(vmin[1], vmin[2], vmin[3])
    _vmax = Vec{3,T}(vmax[1], vmax[2], vmax[3])
    return Rectangle{3,T}(_vmin, _vmax - _vmin)
end

Base.:(-)(h::Rectangle{N,T}, move::Number) where {N,T} = h - Vec{N,T}(ntuple(i->move,N))
Base.:(+)(h::Rectangle{N,T}, move::Number) where {N,T} = h + Vec{N,T}(ntuple(i->move,N))

function Base.:(-)(h::Rectangle{N,T}, move::StaticVector{N}) where {N,T}
    return Rectangle{N,T}(minimum(h) .- move, widths(h))
end

function Base.:(+)(h::Rectangle{N,T}, move::StaticVector{N}) where {N,T}
    return Rectangle{N,T}(minimum(h) .+ move, widths(h))
end

function Base.:(*)(Rectangle::Rectangle, scaling::Union{Number,StaticVector})
    return Rectangle(minimum(Rectangle) .* scaling, widths(Rectangle) .* scaling)
end

# Enables rectangular indexing into a matrix
function Base.to_indices(A::AbstractMatrix{T}, I::Tuple{Rectangle{2,IT}}) where {T,IT<:Integer}
    rect = I[1]
    mini = minimum(rect)
    wh = widths(rect)
    return ((mini[1] + 1):(mini[1] + wh[1]), (mini[2] + 1):(mini[2] + wh[2]))
end

function minmax(p::StaticVector, vmin, vmax)
    any(isnan, p) && return (vmin, vmax)
    return min.(p, vmin), max.(p, vmax)
end

# Annoying special case for view(Vector{Point}, Vector{Face})
function minmax(tup::Tuple, vmin, vmax)
    for p in tup
        any(isnan, p) && continue
        vmin = min.(p, vmin)
        vmax = max.(p, vmax)
    end
    return vmin, vmax
end

function positive_widths(rect::Rectangle{N,T}) where {N,T}
    mini, maxi = minimum(rect), maximum(rect)
    realmin = min.(mini, maxi)
    realmax = max.(mini, maxi)
    return Rectangle{N,T}(realmin, realmax .- realmin)
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
    diff(h1::Rectangle, h2::Rectangle)

Perform a difference between two Rects.
"""
diff(h1::Rectangle, h2::Rectangle) = h1

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

function update(b::Rectangle{N,T}, v::Vec{N,T2}) where {N,T,T2}
    return update(b, Vec{N,T}(v))
end

function update(b::Rectangle{N,T}, v::Vec{N,T}) where {N,T}
    m = min.(minimum(b), v)
    maxi = maximum(b)
    mm = if any(isnan, maxi)
        v - m
    else
        max.(v, maxi) - m
    end
    return Rectangle{N,T}(m, mm)
end

# Min maximum distance functions between hrectangle and point for a given dimension
function min_dist_dim(rect::Rectangle{N,T}, p::Vec{N,T}, dim::Int) where {N,T}
    return max(zero(T), max(minimum(rect)[dim] - p[dim], p[dim] - maximum(rect)[dim]))
end

function max_dist_dim(rect::Rectangle{N,T}, p::Vec{N,T}, dim::Int) where {N,T}
    return max(maximum(rect)[dim] - p[dim], p[dim] - minimum(rect)[dim])
end

function min_dist_dim(rect1::Rectangle{N,T}, rect2::Rectangle{N,T}, dim::Int) where {N,T}
    return max(zero(T),
               max(minimum(rect1)[dim] - maximum(rect2)[dim],
                   minimum(rect2)[dim] - maximum(rect1)[dim]))
end

function max_dist_dim(rect1::Rectangle{N,T}, rect2::Rectangle{N,T}, dim::Int) where {N,T}
    return max(maximum(rect1)[dim] - minimum(rect2)[dim],
               maximum(rect2)[dim] - minimum(rect1)[dim])
end

# Total minimum maximum distance functions
function min_euclideansq(rect::Rectangle{N,T}, p::Union{Vec{N,T},Rectangle{N,T}}) where {N,T}
    minimum_dist = T(0.0)
    for dim in 1:length(p)
        d = min_dist_dim(rect, p, dim)
        minimum_dist += d * d
    end
    return minimum_dist
end

function max_euclideansq(rect::Rectangle{N,T}, p::Union{Vec{N,T},Rectangle{N,T}}) where {N,T}
    maximum_dist = T(0.0)
    for dim in 1:length(p)
        d = max_dist_dim(rect, p, dim)
        maximum_dist += d * d
    end
    return maximum_dist
end

function min_euclidean(rect::Rectangle{N,T}, p::Union{Vec{N,T},Rectangle{N,T}}) where {N,T}
    return sqrt(min_euclideansq(rect, p))
end

function max_euclidean(rect::Rectangle{N,T}, p::Union{Vec{N,T},Rectangle{N,T}}) where {N,T}
    return sqrt(max_euclideansq(rect, p))
end

# Functions that return both minimum and maximum for convenience
function minmax_dist_dim(rect::Rectangle{N,T}, p::Union{Vec{N,T},Rectangle{N,T}},
                         dim::Int) where {N,T}
    minimum_d = min_dist_dim(rect, p, dim)
    maximum_d = max_dist_dim(rect, p, dim)
    return minimum_d, maximum_d
end

function minmax_euclideansq(rect::Rectangle{N,T}, p::Union{Vec{N,T},Rectangle{N,T}}) where {N,T}
    minimum_dist = min_euclideansq(rect, p)
    maximum_dist = max_euclideansq(rect, p)
    return minimum_dist, maximum_dist
end

function minmax_euclidean(rect::Rectangle{N,T}, p::Union{Vec{N,T},Rectangle{N,T}}) where {N,T}
    minimumsq, maximumsq = minmax_euclideansq(rect, p)
    return sqrt(minimumsq), sqrt(maximumsq)
end

# http://en.wikipedia.org/wiki/Allen%27s_interval_algebra
function before(b1::Rectangle{N}, b2::Rectangle{N}) where {N}
    for i in 1:N
        maximum(b1)[i] < minimum(b2)[i] || return false
    end
    return true
end

meets(b1::Rectangle{N}, b2::Rectangle{N}) where {N} = maximum(b1) == minimum(b2)

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

# equality

Base.:(==)(b1::Rectangle, b2::Rectangle) = minimum(b1) == minimum(b2) && widths(b1) == widths(b2)

Base.isequal(b1::Rectangle, b2::Rectangle) = b1 == b2

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
