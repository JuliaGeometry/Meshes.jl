"""
    Cylinder{N, T}

A `Cylinder` is a 2D rectangle or a 3D cylinder defined by
its `origin` point, its `extremity` and a `radius`.
"""
struct Cylinder{N,T} <: GeometryPrimitive{N,T}
    origin::Point{N,T}
    extremity::Point{N,T}
    radius::T
end

"""
    Cylinder2{T}
    Cylinder3{T}

A `Cylinder2` or `Cylinder3` is a 2D/3D cylinder defined by its origin point,
its extremity and a radius. `origin`, `extremity` and `radius`, must be specified.
"""
const Cylinder2{T} = Cylinder{2,T}
const Cylinder3{T} = Cylinder{3,T}

origin(c::Cylinder{N,T}) where {N,T} = c.origin
extremity(c::Cylinder{N,T}) where {N,T} = c.extremity
radius(c::Cylinder{N,T}) where {N,T} = c.radius
height(c::Cylinder{N,T}) where {N,T} = norm(c.extremity - c.origin)
direction(c::Cylinder{N,T}) where {N,T} = (c.extremity - c.origin) ./ height(c)

function rotation(c::Cylinder{2,T}) where {T}
    d2 = direction(c)
    u = @SVector [d2[1], d2[2], T(0)]
    v = @MVector [u[2], -u[1], T(0)]
    normalize!(v)
    return hcat(v, u, @SVector T[0, 0, 1])
end

function rotation(c::Cylinder{3,T}) where {T}
    d3 = direction(c)
    u = @SVector [d3[1], d3[2], d3[3]]
    if abs(u[1]) > 0 || abs(u[2]) > 0
        v = @MVector [u[2], -u[1], T(0)]
    else
        v = @MVector [T(0), -u[3], u[2]]
    end
    normalize!(v)
    w = @SVector [u[2] * v[3] - u[3] * v[2], -u[1] * v[3] + u[3] * v[1],
                  u[1] * v[2] - u[2] * v[1]]
    return hcat(v, w, u)
end

function coordinates(c::Cylinder{2,T}, nvertices=(2, 2)) where {T}
    o = coordinates(c.origin)
    r = c.radius
    h = height(c)
    rect = Rect(o[1] - r / 2, o[2], r, h)
    M = rotation(c)
    points = Point.(coordinates(rect, nvertices))
    vo = to_pointn(Point{3,T}, origin(c))
    return (M * (to_pointn(Point{3,T}, point) - vo) + coordinates(vo) for point in points)
end

function faces(sphere::Cylinder{2}, nvertices=(2, 2))
    return faces(Rect(0, 0, 1, 1), nvertices)
end

function coordinates(c::Cylinder{3,T}, nvertices=30) where {T}
    if isodd(nvertices)
        nvertices = 2 * (nvertices ÷ 2)
    end
    nvertices = max(8, nvertices)
    nbv = nvertices ÷ 2

    M = rotation(c)
    h = height(c)
    range = 1:(2 * nbv + 2)
    function inner(i)
        return if i == length(range)
            coordinates(c.extremity)
        elseif i == length(range) - 1
            coordinates(c.origin)
        else
            phi = T((2π * (((i + 1) ÷ 2) - 1)) / nbv)
            up = ifelse(isodd(i), T(0), h)
            o  = coordinates(c.origin)
            r  = c.radius
            (M * Vec(r*cos(phi), r*sin(phi), up)) + o
        end
    end

    return (inner(i) for i in range)
end

function faces(c::Cylinder{3}, facets=30)
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing
    nbv = Int(facets / 2)
    indexes = Vector{TriangleFace{Int}}(undef, facets)
    index = 1
    for j in 1:(nbv - 1)
        indexes[index] = (index + 2, index + 1, index)
        indexes[index + 1] = (index + 3, index + 1, index + 2)
        index += 2
    end
    indexes[index] = (1, index + 1, index)
    indexes[index + 1] = (2, index + 1, 1)

    for i in 1:length(indexes)
        i % 2 == 1 ? push!(indexes, (indexes[i][1], indexes[i][3], 2 * nbv + 1)) :
        push!(indexes, (indexes[i][2], indexes[i][1], 2 * nbv + 2))
    end
    return indexes
end
