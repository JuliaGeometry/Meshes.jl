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

The torus whose centerline passes through points `p1`, `p2` and `p3` and with
minor radius `minor`.
"""
function Torus(p1::Point{3,T}, p2::Point{3,T}, p3::Point{3,T}, minor) where {T}
  c = Circle(p1, p2, p3)
  p = Plane(p1, p2, p3)
  Torus(center(c), normal(p), radius(c), T(minor))
end

Torus(p1::Tuple, p2::Tuple, p3::Tuple, minor) = Torus(Point(p1), Point(p2), Point(p3), minor)

paramdim(::Type{<:Torus}) = 2

center(t::Torus) = t.center

normal(t::Torus) = t.normal

radii(t::Torus) = (t.major, t.minor)

axis(t::Torus) = Line(t.center, t.center + t.normal)

function (t::Torus{T})(θ, φ) where {T}
  if (θ < 0 || θ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((θ, φ), "t(θ, φ) is not defined for θ, φ outside [0, 1]²."))
  end
  c, n⃗ = t.center, t.normal
  R, r = t.major, t.minor

  Q = rotation_between(Vec{3,T}(0, 0, 1), n⃗)

  sθ, cθ = sincospi(2 * T(-θ))
  sφ, cφ = sincospi(2 * T(φ))
  x = (R + r * cθ) * cφ
  y = (R + r * cθ) * sφ
  z = r * sθ

  c + Q * Vec{3,T}(x, y, z)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Torus{T}}) where {T} =
  Torus(rand(rng, Point{3,T}), rand(rng, Vec{3,T}), rand(rng, T), rand(rng, T))
