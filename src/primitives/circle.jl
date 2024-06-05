# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Circle(plane, radius)

A circle embedded in 3-dimensional space on a
given `plane` with given `radius`.

See also [`Disk`](@ref).
"""
struct Circle{C<:CRS,P<:Plane{C},ℒ<:Len} <: Primitive{3,C}
  plane::P
  radius::ℒ
  Circle{C,P,ℒ}(plane, radius) where {C<:CRS,P<:Plane{C},ℒ<:Len} = new(plane, radius)
end

Circle(plane::P, radius::ℒ) where {C<:CRS,P<:Plane{C},ℒ<:Len} = Circle{C,P,float(ℒ)}(plane, radius)

Circle(plane::Plane, radius) = Circle(plane, addunit(radius, u"m"))

"""
    Circle(p1, p2, p3)

A circle passing through points `p1`, `p2` and `p3`.
"""
function Circle(p1::Point{3}, p2::Point{3}, p3::Point{3})
  v12 = p2 - p1
  v13 = p3 - p1
  m12 = to(p1 + v12 / 2)
  m13 = to(p1 + v13 / 2)
  n⃗ = normal(Plane(p1, p2, p3))
  F = to(p1) ⋅ n⃗
  M = transpose([n⃗ v12 v13])
  u = [F, m12 ⋅ v12, m13 ⋅ v13]
  O = withdatum(p1, uinv(M) * u)
  r = norm(p1 - O)
  Circle(Plane(O, n⃗), r)
end

Circle(p1::Tuple, p2::Tuple, p3::Tuple) = Circle(Point(p1), Point(p2), Point(p3))

paramdim(::Type{<:Circle}) = 1

plane(c::Circle) = c.plane

center(c::Circle) = c.plane(0, 0)

radius(c::Circle) = c.radius

function (c::Circle)(φ)
  T = numtype(lentype(c))
  if (φ < 0 || φ > 1)
    throw(DomainError(φ, "c(φ) is not defined for φ outside [0, 1]."))
  end
  r = c.radius
  l = r
  sφ, cφ = sincospi(2 * T(φ))
  u = ustrip(l * cφ)
  v = ustrip(l * sφ)
  c.plane(u, v)
end

Random.rand(rng::Random.AbstractRNG, ::Type{Circle}) = Circle(rand(rng, Plane), rand(rng, Met{Float64}))
