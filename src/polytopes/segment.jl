# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

A line segment with points `p1`, `p2`.
"""
struct Segment{Dim,T} <: Polytope{Dim,T,2}
  vertices::NTuple{2,Point{Dim,T}}
end

"""
    intersect(l1, l2)

Intersection of line segmens `l1` and `l2`.
"""
function Base.intersect(a::Segment{2,T1}, b::Segment{2,T2}) where {T1,T2}
    T = promote_type(T1, T2)
    v1, v2 = coordinates.(a)
    v3, v4 = coordinates.(b)
    MT = Mat{2,2,T,4}
    p0 = Point{2,T}(0, 0)

    verticalA = v1[1] == v2[1]
    verticalB = v3[1] == v4[1]

    # if a segment is vertical the linear algebra might have trouble
    # so we will rotate the segments such that neither is vertical
    dorotation = verticalA || verticalB

    if dorotation
        θ = T(0.0)
        if verticalA && verticalB
            θ = T(π / 4)
        elseif verticalA || verticalB # obviously true, but make it clear
            θ34 = -atan(v4[2] - v3[2], v4[1] - v3[1])
            θ12 = -atan(v2[2] - v1[2], v2[1] - v1[1])
            θ = verticalA ? θ34 : θ12
            θ = abs(θ) == T(0) ? (θ12 + θ34) / 2 : θ
            θ = abs(θ) == T(pi) ? (θ12 + θ34) / 2 : θ
        end
        rotation = MT(cos(θ), sin(θ), -sin(θ), cos(θ))
        v1 = rotation * v1
        v2 = rotation * v2
        v3 = rotation * v3
        v4 = rotation * v4
    end

    a = det(MT(v1[1] - v2[1], v1[2] - v2[2], v3[1] - v4[1], v3[2] - v4[2]))

    (abs(a) < eps(T)) && return false, p0 # Lines are parallel

    d1 = det(MT(v1[1], v1[2], v2[1], v2[2]))
    d2 = det(MT(v3[1], v3[2], v4[1], v4[2]))
    x = det(MT(d1, v1[1] - v2[1], d2, v3[1] - v4[1])) / a
    y = det(MT(d1, v1[2] - v2[2], d2, v3[2] - v4[2])) / a

    (x < prevfloat(min(v1[1], v2[1])) || x > nextfloat(max(v1[1], v2[1]))) &&
        return false, p0
    (y < prevfloat(min(v1[2], v2[2])) || y > nextfloat(max(v1[2], v2[2]))) &&
        return false, p0
    (x < prevfloat(min(v3[1], v4[1])) || x > nextfloat(max(v3[1], v4[1]))) &&
        return false, p0
    (y < prevfloat(min(v3[2], v4[2])) || y > nextfloat(max(v3[2], v4[2]))) &&
        return false, p0

    point = Point{2,T}(x, y)

    # don't forget to rotate the answer back
    if dorotation
        point = Point(transpose(rotation) * coordinates(point))
    end

    return true, point
end
