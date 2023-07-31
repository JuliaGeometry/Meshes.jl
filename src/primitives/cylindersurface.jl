# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CylinderSurface(bottom, top, radius)

A circular cylinder surface embedded in R³ with given `radius`,
delimited by `bottom` and `top` planes.

    CylinderSurface(segment, radius)

Alternatively, construct a right circular cylinder surface with given `radius`
and `segment` between `bottom` and `top` planes.

    CylinderSurface(segment)

Or construct a right circular cylinder surface with unit radius along the `segment`.

    CylinderSurface(radius)

Finally, construct a right vertical circular cylinder surface with given `radius`.

See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct CylinderSurface{T} <: Primitive{3,T}
  bot::Plane{T}
  top::Plane{T}
  radius::T
end

function CylinderSurface(start::Point{3,T}, finish::Point{3,T}, radius) where {T}
  dir = finish - start
  bot = Plane(start, dir)
  top = Plane(finish, dir)
  CylinderSurface(bot, top, T(radius))
end

CylinderSurface(start::Tuple, finish::Tuple, radius) = CylinderSurface(Point(start), Point(finish), radius)

CylinderSurface(start::Point{3,T}, finish::Point{3,T}) where {T} = CylinderSurface(start, finish, T(1))

CylinderSurface(start::Tuple, finish::Tuple) = CylinderSurface(Point(start), Point(finish))

CylinderSurface(radius::T) where {T} = CylinderSurface(Point(T(0), T(0), T(0)), Point(T(0), T(0), T(1)), radius)

paramdim(::Type{<:CylinderSurface}) = 2

isparametrized(::Type{<:CylinderSurface}) = true

radius(c::CylinderSurface) = c.radius

bottom(c::CylinderSurface) = c.bot

top(c::CylinderSurface) = c.top

function center(c::CylinderSurface)
  a = coordinates(c.bot(0, 0))
  b = coordinates(c.top(0, 0))
  Point((a .+ b) ./ 2)
end

axis(c::CylinderSurface) = Line(c.bot(0, 0), c.top(0, 0))

function isright(c::CylinderSurface{T}) where {T}
  # cylinder is right if axis
  # is aligned with plane normals
  a = axis(c)
  d = a(T(1)) - a(T(0))
  v = normal(c.bot)
  w = normal(c.top)
  isparallelv = isapprox(norm(d × v), zero(T), atol=atol(T))
  isparallelw = isapprox(norm(d × w), zero(T), atol=atol(T))
  isparallelv && isparallelw
end

measure(c::CylinderSurface{T}) where {T} = (norm(c.bot(0, 0) - c.top(0, 0)) + c.radius) * 2 * c.radius * T(π)

area(c::CylinderSurface) = measure(c)

boundary(::CylinderSurface) = nothing

Base.isapprox(c₁::CylinderSurface{T}, c₂::CylinderSurface{T}) where {T} =
  c₁.bot ≈ c₂.bot && c₁.top ≈ c₂.top && isapprox(c₁.radius, c₂.radius, atol=atol(T))

function (c::CylinderSurface{T})(φ, z) where {T}
  if (φ < 0 || φ > 1) || (z < 0 || z > 1)
    throw(DomainError((φ, z), "c(φ, z) is not defined for φ, z outside [0, 1]²."))
  end
  o = center(c)
  r = radius(c)
  b = bottom(c)
  t = top(c)
  a = axis(c)
  d = a(1) - a(0)
  l = norm(d)

  # rotation to align z axis with cylinder axis
  Q = rotation_between(d, Vec{3,T}(0, 0, 1))

  # new normals of planes in new rotated system
  nᵦ = Q * normal(b)
  nₜ = Q * normal(t)

  # scale coordinates
  φₛ = 2T(π) * φ
  zₛ = z

  # given cylindrical coordinates (r*cos(φ), r*sin(φ), z) and the
  # equation of the plane, we can solve for z and find all points
  # along the ellipse obtained by intersection
  zᵦ = -l / 2 - (r * cos(φₛ) * nᵦ[1] + r * sin(φₛ) * nᵦ[2]) / nᵦ[3]
  zₜ = +l / 2 - (r * cos(φₛ) * nₜ[1] + r * sin(φₛ) * nₜ[2]) / nₜ[3]
  pᵦ = Point(r * cos(φₛ), r * sin(φₛ), zᵦ)
  pₜ = Point(r * cos(φₛ), r * sin(φₛ), zₜ)

  p = pᵦ + zₛ * (pₜ - pᵦ)
  o + Q' * coordinates(p)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{CylinderSurface{T}}) where {T} =
  CylinderSurface(rand(rng, Plane{T}), rand(rng, Plane{T}), rand(rng, T))
