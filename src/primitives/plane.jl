# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Plane(p, u, v)

A plane embedded in R³ passing through point `p`,
defined by non-parallel vectors `u` and `v`.

    Plane(p, n)

Alternatively specify point `p` and a given normal
vector `n` to the plane.
"""
struct Plane{T} <: Primitive{3,T}
  p::Point{3,T}
  u::Vec{3,T}
  v::Vec{3,T}
end

function Plane(p::Point{3,T}, n::Vec{3,T}) where {T}
  # origin of coordinate system
  o = Vec{3,T}(0, 0, 0)

  uv = Vec{3,T}[]
  for i in 1:3
    # subtract projection of Euclidean basis onto normal
    e = Vec{3,T}(ntuple(j -> j == i ? 1 : 0, 3))
    v = e - (e⋅n)/(n⋅n) * n

    # check if it is a valid vector
    isapprox(v, o, atol=atol(T)) || push!(uv, v)

    # we need two vectors
    length(uv) == 2 && break
  end

  Plane(p, uv[1], uv[2])
end

Plane(p::Tuple, u::Tuple, v::Tuple) =
  Plane(Point(p), Vec(u), Vec(v))

Plane(p::Tuple, n::Tuple) =
  Plane(Point(p), Vec(n))

paramdim(::Type{<:Plane}) = 2

isconvex(::Type{<:Plane}) = true

boundary(::Plane) = nothing

(p::Plane)(s, t) = p.p + s*p.u + t*p.v

"""
    origin(plane)

Return the origin of the `plane`.
"""
origin(p::Plane) = p.p

"""
    normal(plane)

Normal vector to the `plane`.
"""
normal(p::Plane) = p.u × p.v
