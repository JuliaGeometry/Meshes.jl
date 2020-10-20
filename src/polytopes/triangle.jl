# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Triangle(p1, p2, p3)

A triangle with points `p1`, `p2`, `p3`.
"""
struct Triangle{Dim,T} <: Polytope{Dim,T,3}
    vertices::SVector{3,Point{Dim,T}}
end

function volume(triangle::Triangle)
    v1, v2, v3 = coordinates.(triangle)
    sig = sign(orthogonal_vector(v1, v2, v3) ⋅ v1)
    return sig * abs(v1 ⋅ (v2 × v3)) / 6
end

function Base.in(P::T, triangle::Triangle) where {T<:Point}
    A, B, C = vertices(triangle)
    a = C - B
    b = A - C
    c = B - A

    ap = P - A
    bp = P - B
    cp = P - C

    a_bp = a[1] * bp[2] - a[2] * bp[1]
    c_ap = c[1] * ap[2] - c[2] * ap[1]
    b_cp = b[1] * cp[2] - b[2] * cp[1]

    return ((a_bp ≥ 0) && (b_cp ≥ 0) && (c_ap ≥ 0))
end
