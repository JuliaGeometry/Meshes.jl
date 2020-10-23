# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cylinder(start, finish, radius)

A right circular cylinder with `start` and `finish` points,
and `radius` of revolution. See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct Cylinder{T} <: Primitive{3,T}
    start::Point{3,T}
    finish::Point{3,T}
    radius::T
end

radius(c::Cylinder) = c.radius
height(c::Cylinder) = norm(c.finish - c.start)
volume(c::Cylinder) = π * radius(c)^2 * height(c)

# TODO: review these
function rotation(c::Cylinder{T}) where T
    d = c.finish - c.start
    u = d ./ norm(d)
    if abs(u[1]) > 0 || abs(u[2]) > 0
        v = MVector(u[2], -u[1], T(0))
    else
        v = MVector(T(0), -u[3], u[2])
    end
    normalize!(v)
    w = SVector(u[2] * v[3] - u[3] * v[2],
               -u[1] * v[3] + u[3] * v[1],
                u[1] * v[2] - u[2] * v[1])
    return hcat(v, w, u)
end

function coordinates(c::Cylinder{T}, nvertices=30) where {T}
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
            coordinates(c.finish)
        elseif i == length(range) - 1
            coordinates(c.start)
        else
            phi = T((2π * (((i + 1) ÷ 2) - 1)) / nbv)
            up = ifelse(isodd(i), T(0), h)
            o  = coordinates(c.start)
            r  = c.radius
            (M * Vec(r*cos(phi), r*sin(phi), up)) + o
        end
    end

    return (inner(i) for i in range)
end

function faces(::Cylinder, facets=30)
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing
    nbv = Int(facets / 2)
    indexes = Vector{TriangleFace}(undef, facets)
    index = 1
    for j in 1:(nbv - 1)
        indexes[index] = TriangleFace((index + 2, index + 1, index))
        indexes[index + 1] = TriangleFace((index + 3, index + 1, index + 2))
        index += 2
    end
    indexes[index] = TriangleFace((1, index + 1, index))
    indexes[index + 1] = TriangleFace((2, index + 1, 1))

    for i in 1:length(indexes)
        i % 2 == 1 ? push!(indexes, TriangleFace((indexes[i][1], indexes[i][3], 2 * nbv + 1))) :
        push!(indexes, TriangleFace((indexes[i][2], indexes[i][1], 2 * nbv + 2)))
    end
    return indexes
end
