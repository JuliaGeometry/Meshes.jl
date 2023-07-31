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

Torus(center::Point{3,T}, normal::Vec{3,T}, major, minor) where {T} = Torus(center, normal, T(major), T(minor))

Torus(center::Tuple, normal::Tuple, major, minor) = Torus(Point(center), Vec(normal), major, minor)

"""
    Torus(p1, p2, p3, minor)

The torus whose equator passes through points `p1`, `p2` and `p3` and with
minor radius `minor`.
"""
function Torus(p1::Point{3,T}, p2::Point{3,T}, p3::Point{3,T}, minor) where {T}
  c = Circle(p1, p2, p3)
  p = Plane(p1, p2, p3)
  Torus(center(c), normal(p), radius(c), T(minor))
end

Torus(p1::Tuple, p2::Tuple, p3::Tuple, minor) = Torus(Point(p1), Point(p2), Point(p3), minor)

paramdim(::Type{<:Torus}) = 2

isperiodic(::Type{<:Torus}) = (true, true)

center(t::Torus) = t.center

normal(t::Torus) = t.normal

radii(t::Torus) = (t.major, t.minor)

axis(t::Torus) = Line(t.center, t.center + t.normal)

# https://en.wikipedia.org/wiki/Torus
function measure(t::Torus{T}) where {T}
  R, r = t.major, t.minor
  4T(π)^2 * R * r
end

area(t::Torus) = measure(t)

boundary(::Torus) = nothing

function Base.in(p::Point{3,T}, t::Torus{T}) where {T}
  c, n⃗ = t.center, t.normal
  R, r = t.major, t.minor
  Q = rotation_between(n⃗, Vec{3,T}(0, 0, 1))
  x, y, z = Q * (p - c)
  (R - √(x^2 + y^2))^2 + z^2 ≤ r^2
end

function (t::Torus{T})(u, v) where {T}
  if (u < 0 || u > 1) || (v < 0 || v > 1)
    throw(DomainError((u, v), "t(u, v) is not defined for u, v outside [0, 1]²."))
  end
  c, n⃗ = t.center, t.normal
  R, r = t.major, t.minor
  Q = rotation_between(Vec{3,T}(0, 0, 1), n⃗)
  kxy = R^2 - r^2
  kz = √kxy * r
  uₛ = T(π) * (2 * u - 1)
  vₛ = T(π) * (2 * v - 1)
  k = R - r * cos(vₛ)
  x = kxy * cos(uₛ) / k
  y = kxy * sin(uₛ) / k
  z = kz * sin(vₛ) / k
  c + Q * Vec{3,T}(x, y, z)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Torus{T}}) where {T} =
  Torus(rand(rng, Point{3,T}), rand(rng, Vec{3,T}), rand(rng, T), rand(rng, T))
