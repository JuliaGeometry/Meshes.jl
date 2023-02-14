# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Circle(plane, radius)

A circle on a `plane` with given `center` and `radius`.

See also [`Disk`](@ref).
"""
struct Circle{T} <: Primitive{3,T}
  plane::Plane{T}
  center::Point{3,T}
  radius::T
end

Circle(plane::Plane{T}, center::Point{3,T}, radius) where {T} = 
  Circle(plane, center, T(radius))

Circle(plane::Plane{T}, center::Tuple, radius) where {T} = 
  Circle(plane, Point(center), radius)

"""
    Circle(p1, p2, p3)

A 3D circle passing through points `p1`, `p2` and `p3`.
"""
function Circle(p1::Point{3}, p2::Point{3}, p3::Point{3})
  v12 = p2 - p1
  v13 = p3 - p1
  m12 = coordinates(p1 + v12/2)
  m13 = coordinates(p1 + v13/2)
  p = Plane(p1, p2, p3)
  n⃗ = normal(p)
  F = coordinates(p1) ⋅ n⃗
  M = transpose(hcat(n⃗, v12, v13))
  u = [F, m12 ⋅ v12, m13 ⋅ v13]
  O = Point(inv(M) * u)
  r = norm(p1 - O)
  Circle(p, O, r)
end

Circle(p1::Tuple, p2::Tuple, p3::Tuple) =
  Circle(Point(p1), Point(p2), Point(p3))



paramdim(::Type{<:Circle}) = 1

isconvex(::Type{<:Circle}) = false

isperiodic(::Type{<:Circle}) = (true,)

center(c::Circle) = c.center

radius(c::Circle) = c.radius

measure(c::Circle{T}) where {T} = 2 * T(π) * c.radius

Base.length(c::Circle) = measure(c)

function Base.in(p::Point{3,T}, c::Circle{T}) where {T}
  p ∉ c.plane && return false
  s² = sum(abs2, p - center(c))
  r² = radius(c)^2
  isapprox(s², r², atol = atol(T)^2)
end

boundary(::Circle) = nothing