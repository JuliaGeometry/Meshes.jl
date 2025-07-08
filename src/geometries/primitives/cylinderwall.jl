# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CylinderWall(bottom, top, radius)

A circular cylinder wall embedded in R³ with given `radius`,
delimited by `bottom` and `top` planes.

    CylinderWall(start, finish, radius)

Alternatively, construct a right circular cylinder wall with given `radius`
along the segment with `start` and `finish` end points.

    CylinderWall(start, finish)

Or construct a right circular cylinder wall with unit radius along the segment
with `start` and `finish` end points.

    CylinderWall(radius)

Finally, construct a right vertical circular cylinder wall with given `radius`.

See also [`CylinderSurface`](@ref) and [`Cylinder`](@ref).
"""
struct CylinderWall{C<:CRS,P<:Plane{C},ℒ<:Len} <: Primitive{𝔼{3},C}
  bot::P
  top::P
  radius::ℒ
  CylinderWall(bot::P, top::P, radius::ℒ) where {C<:CRS,P<:Plane{C},ℒ<:Len} = new{C,P,float(ℒ)}(bot, top, radius)
end

CylinderWall(bot::P, top::P, radius) where {P<:Plane} = CylinderWall(bot, top, addunit(radius, u"m"))

function CylinderWall(start::Point, finish::Point, radius)
  dir = finish - start
  bot = Plane(start, dir)
  top = Plane(finish, dir)
  CylinderWall(bot, top, radius)
end

CylinderWall(start::Tuple, finish::Tuple, radius) = CylinderWall(Point(start), Point(finish), radius)

CylinderWall(start::Point, finish::Point) = CylinderWall(start, finish, oneunit(lentype(start)))

CylinderWall(start::Tuple, finish::Tuple) = CylinderWall(Point(start), Point(finish))

function CylinderWall(radius)
  z = zero(radius)
  o = oneunit(radius)
  CylinderWall(Point(z, z, z), Point(z, z, o), radius)
end

paramdim(::Type{<:CylinderWall}) = 2

bottom(c::CylinderWall) = c.bot

top(c::CylinderWall) = c.top

radius(c::CylinderWall) = c.radius

axis(c::CylinderWall) = Line(c.bot(0, 0), c.top(0, 0))

function isright(c::CylinderWall)
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

function hasintersectingplanes(c::CylinderWall)
  x = c.bot ∩ c.top
  !isnothing(x) && evaluate(Euclidean(), axis(c), x) < c.radius
end

==(c₁::CylinderWall, c₂::CylinderWall) = c₁.bot == c₂.bot && c₁.top == c₂.top && c₁.radius == c₂.radius

Base.isapprox(c₁::CylinderWall, c₂::CylinderWall; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(c₁.bot, c₂.bot; atol, kwargs...) &&
  isapprox(c₁.top, c₂.top; atol, kwargs...) &&
  isapprox(c₁.radius, c₂.radius; atol, kwargs...)

(c::CylinderWall)(φ, z) = Cylinder(bottom(c), top(c), radius(c))(1, φ, z)
