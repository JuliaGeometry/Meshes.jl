# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Circle(plane, radius)

A circle embedded in 3-dimensional space on a
given `plane` with given `radius`.

See also [`Disk`](@ref).
"""
struct Circle{C<:CRS,P<:Plane{C},ℒ<:Len} <: Primitive{𝔼{3},C}
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

radius(c::Circle) = c.radius

center(c::Circle) = plane(c)(0, 0)

==(c₁::Circle, c₂::Circle) = plane(c₁) == plane(c₂) && radius(c₁) == radius(c₂)

Base.isapprox(c₁::Circle, c₂::Circle; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(plane(c₁), plane(c₂); atol, kwargs...) && isapprox(radius(c₁), radius(c₂); atol, kwargs...)

(c::Circle)(φ) = Disk(plane(c), radius(c))(1, φ)
