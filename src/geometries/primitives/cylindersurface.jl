# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CylinderSurface(bottom, top, radius)

A circular cylinder surface embedded in R³ with given `radius`,
delimited by `bottom` and `top` planes.

    CylinderSurface(start, finish, radius)

Alternatively, construct a right circular cylinder surface with given `radius`
along the segment with `start` and `finish` end points.

    CylinderSurface(start, finish)

Or construct a right circular cylinder surface with unit radius along the segment
with `start` and `finish` end points.

    CylinderSurface(radius)

Finally, construct a right vertical circular cylinder surface with given `radius`.

See also [`CylinderWall`](@ref) and [`Cylinder`](@ref).
"""
struct CylinderSurface{C<:CRS,P<:Plane{C},ℒ<:Len} <: Primitive{𝔼{3},C}
  bot::P
  top::P
  radius::ℒ
  CylinderSurface(bot::P, top::P, radius::ℒ) where {C<:CRS,P<:Plane{C},ℒ<:Len} = new{C,P,float(ℒ)}(bot, top, radius)
end

CylinderSurface(bot::P, top::P, radius) where {P<:Plane} = CylinderSurface(bot, top, addunit(radius, u"m"))

function CylinderSurface(start::Point, finish::Point, radius)
  dir = finish - start
  bot = Plane(start, dir)
  top = Plane(finish, dir)
  CylinderSurface(bot, top, radius)
end

CylinderSurface(start::Tuple, finish::Tuple, radius) = CylinderSurface(Point(start), Point(finish), radius)

CylinderSurface(start::Point, finish::Point) = CylinderSurface(start, finish, oneunit(lentype(start)))

CylinderSurface(start::Tuple, finish::Tuple) = CylinderSurface(Point(start), Point(finish))

function CylinderSurface(radius)
  z = zero(radius)
  o = oneunit(radius)
  CylinderSurface(Point(z, z, z), Point(z, z, o), radius)
end

paramdim(::Type{<:CylinderSurface}) = 2

bottom(c::CylinderSurface) = c.bot

top(c::CylinderSurface) = c.top

radius(c::CylinderSurface) = c.radius

wall(c::CylinderSurface) = CylinderWall(c.bot, c.top, c.radius)

axis(c::CylinderSurface) = Line(c.bot(0, 0), c.top(0, 0))

function isright(c::CylinderSurface)
  ℒ = lentype(c)
  T = numtype(ℒ)
  # cylinder is right if axis
  # is aligned with plane normals
  a = axis(c)
  d = a(T(1)) - a(T(0))
  v = normal(c.bot)
  w = normal(c.top)
  isparallelv = isapproxzero(norm(d × v))
  isparallelw = isapproxzero(norm(d × w))
  isparallelv && isparallelw
end

==(c₁::CylinderSurface, c₂::CylinderSurface) = c₁.bot == c₂.bot && c₁.top == c₂.top && c₁.radius == c₂.radius

Base.isapprox(c₁::CylinderSurface, c₂::CylinderSurface; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(c₁.bot, c₂.bot; atol, kwargs...) &&
  isapprox(c₁.top, c₂.top; atol, kwargs...) &&
  isapprox(c₁.radius, c₂.radius; atol, kwargs...)

(c::CylinderSurface)(φ, z) = Cylinder(bottom(c), top(c), radius(c))(1, φ, z)

function hasintersectingplanes(c::CylinderSurface)
  x = c.bot ∩ c.top
  !isnothing(x) && evaluate(Euclidean(), axis(c), x) < c.radius
end
