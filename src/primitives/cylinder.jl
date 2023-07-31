# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cylinder(bottom, top, radius)

A solid circular cylinder embedded in R³ with given `radius`,
delimited by `bottom` and `top` planes.

    Cylinder(segment, radius)

Alternatively, construct a right circular cylinder with given `radius`
and `segment` between `bottom` and `top` planes.

    Cylinder(segment)

Or construct a right circular cylinder with unit radius along the `segment`.

    Cylinder(radius)

Finally, construct a right vertical circular cylinder with given `radius`.

See https://en.wikipedia.org/wiki/Cylinder. 
"""
struct Cylinder{T} <: Primitive{3,T}
  bot::Plane{T}
  top::Plane{T}
  radius::T
end

function Cylinder(start::Point{3,T}, finish::Point{3,T}, radius) where {T}
  dir = finish - start
  bot = Plane(start, dir)
  top = Plane(finish, dir)
  Cylinder(bot, top, T(radius))
end

Cylinder(start::Tuple, finish::Tuple, radius) = Cylinder(Point(start), Point(finish), radius)

Cylinder(start::Point{3,T}, finish::Point{3,T}) where {T} = Cylinder(start, finish, T(1))

Cylinder(start::Tuple, finish::Tuple) = Cylinder(Point(start), Point(finish))

Cylinder(radius::T) where {T} = Cylinder(Point(T(0), T(0), T(0)), Point(T(0), T(0), T(1)), radius)

paramdim(::Type{<:Cylinder}) = 3

isparametrized(::Type{<:Cylinder}) = true

radius(c::Cylinder) = c.radius

bottom(c::Cylinder) = c.bot

top(c::Cylinder) = c.top

center(c::Cylinder) = center(boundary(c))

axis(c::Cylinder) = axis(boundary(c))

isright(c::Cylinder) = isright(boundary(c))

measure(c::Cylinder{T}) where {T} = norm(c.bot(0, 0) - c.top(0, 0)) * T(π) * c.radius^2

volume(c::Cylinder) = measure(c)

boundary(c::Cylinder) = CylinderSurface(c.bot, c.top, c.radius)

Base.isapprox(c₁::Cylinder, c₂::Cylinder) = boundary(c₁) ≈ boundary(c₂)

function Base.in(p::Point{3}, c::Cylinder)
  b = c.bot(0, 0)
  t = c.top(0, 0)
  a = t - b
  (p - b) ⋅ a ≥ 0 || return false
  (p - t) ⋅ a ≤ 0 || return false
  norm((p - b) × a) / norm(a) ≤ c.radius
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Cylinder{T}}) where {T} =
  Cylinder(rand(rng, Plane{T}), rand(rng, Plane{T}), rand(rng, T))
