# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Circle(plane, radius)

A circle embedded in 3-dimensional space on a
given `plane` with given `radius`.

See also [`Disk`](@ref).
"""
struct Circle{C<:CRS,P<:Plane{C},ℒ<:Len} <: Primitive{C}
  plane::P
  radius::ℒ
  Circle(plane::P, radius::ℒ) where {C<:CRS,P<:Plane{C},ℒ<:Len} = new{C,P,float(ℒ)}(plane, radius)
end

Circle(plane::Plane, radius) = Circle(plane, addunit(radius, u"m"))

"""
    Circle(p1, p2, p3)

A circle passing through points `p1`, `p2` and `p3`.
"""
function Circle(p1::Point, p2::Point, p3::Point)
  v12 = p2 - p1
  v13 = p3 - p1
  m12 = to(p1 + v12 / 2)
  m13 = to(p1 + v13 / 2)
  n⃗ = normal(Plane(p1, p2, p3))
  F = to(p1) ⋅ n⃗
  M = transpose([n⃗ v12 v13])
  u = [F, m12 ⋅ v12, m13 ⋅ v13]
  O = withcrs(p1, uinv(M) * u)
  r = norm(p1 - O)
  Circle(Plane(O, n⃗), r)
end

Circle(p1::Tuple, p2::Tuple, p3::Tuple) = Circle(Point(p1), Point(p2), Point(p3))

paramdim(::Type{<:Circle}) = 1

plane(c::Circle) = c.plane

center(c::Circle) = c.plane(0, 0)

radius(c::Circle) = c.radius

==(c₁::Circle, c₂::Circle) = c₁.plane == c₂.plane && c₁.radius == c₂.radius

Base.isapprox(c₁::Circle, c₂::Circle; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(c₁.plane, c₂.plane; atol, kwargs...) && isapprox(c₁.radius, c₂.radius; atol, kwargs...)

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
