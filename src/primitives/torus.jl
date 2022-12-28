# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Torus(center, normal, major, minor)

A torus centered at `center` with axis of revolution directed by 
`normal` and with radii `major` and `minor`. 

"""
struct Torus{T} <: Primitive{3,T}
  center::Point{3,T}
  normal::Vec{3,T}
  major::T
  minor::T
end

function Torus(center::Tuple, normal::Tuple, major::T1, minor::T2) where {T1,T2}
  T = promote_type(eltype(center), eltype(normal), T1, T2, Float32)  
  Torus(Point(T.(center)), Vec(T.(normal)), T(major), T(minor))
end

"""
  Torus(p1, p2, p3)

The torus whose equator passes through points `p1`, `p2` and `p3` and with
minor radius `minor`.
"""
function Torus(p1::Point{3}, p2::Point{3}, p3::Point{3}, minor)
  # circumcenter of p1, p2, p3
  n⃗ = normal(Plane(p1, p2, p3))
  v12 = p2 - p1
  v13 = p3 - p1
  p12 = coordinates(p1 + v12/2)
  p13 = coordinates(p1 + v13/2)
  offset = coordinates(p1) ⋅ n⃗
  A = transpose(hcat(n⃗, v12, v13))
  b = [offset, p12 ⋅ v12, p13 ⋅ v13]
  c = Point(inv(A) * b)
  # circumradius
  major = norm(p1 - c)
  #
  T = typeof(major)
  Torus(c, Vec{3,T}(n⃗), major, T(minor))
end

paramdim(::Type{<:Torus}) = 2

isconvex(::Type{<:Torus}) = false

isperiodic(::Type{<:Torus}) = (true, true)

center(t::Torus) = t.center

normal(t::Torus) = t.normal

radii(t::Torus) = (t.major, t.minor)

axis(t::Torus) = Line(t.center, t.center + t.normal)

# https://en.wikipedia.org/wiki/Torus
function measure(t::Torus)
  R, r = t.major, t.minor
  4π^2 * R * r
end

area(t::Torus) = measure(t)

function Base.in(p::Point, t::Torus)
  c, n⃗ = t.center, t.normal
  R, r = t.major, t.minor
  M = uvrotation(Vec(0, 0, 1), n⃗)
  x, y, z = M * (p - c)
  (R - √(x^2 + y^2))^2 + z^2 ≤ r^2
end
