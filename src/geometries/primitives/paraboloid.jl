# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Paraboloid(apex, radius, focallength)

A paraboloid embedded in R³ and extending up to a distance `radius`
from its focal axis, which is aligned along the z direction and
passes through `apex` (the apex of the paraboloid). The equation
of the paraboloid is the following:

```math
f(x, y) = \\frac{(x - x_0)^2 + (y - y_0)^2}{4f} + z_0\\qquad\\text{for } (x - x_0)^2 + (y - y_0)^2 < r^2,
```

where ``(x_0, y_0, z_0)`` is the apex of the parabola, ``f`` is the
focal length, and ``r`` is the clip radius.

    Paraboloid(apex, radius)

This creates a paraboloid surface with focal length equal to 1.

    Paraboloid(apex)

This creates a paraboloid surface with focal length equal to 1 and a rim with unit
radius.

    Paraboloid()

Same as above, but here the apex is at `(0, 0, 0)`.

See <https://en.wikipedia.org/wiki/Paraboloid>.
"""
struct Paraboloid{C<:CRS,Mₚ<:Manifold,ℒ<:Len} <: Primitive{𝔼{3},C}
  apex::Point{Mₚ,C}
  radius::ℒ
  focallength::ℒ
  Paraboloid(apex::Point{Mₚ,C}, radius::ℒ, focallength::ℒ) where {C<:CRS,Mₚ<:Manifold,ℒ<:Len} =
    new{C,Mₚ,float(ℒ)}(apex, radius, focallength)
end

Paraboloid(apex::Point, radius::Len, focallength::Len) = Paraboloid(apex, promote(radius, focallength)...)

Paraboloid(apex::Point, radius, focallength) = Paraboloid(apex, aslen(radius), aslen(focallength))

Paraboloid(apex::Tuple, radius, focallength) = Paraboloid(Point(apex), radius, focallength)

Paraboloid(apex::Point, radius) = Paraboloid(apex, radius, oneunit(radius))

Paraboloid(apex::Tuple, radius) = Paraboloid(Point(apex), radius)

Paraboloid(apex::Point) = Paraboloid(apex, oneunit(lentype(apex)))

Paraboloid(apex::Tuple) = Paraboloid(Point(apex))

Paraboloid() = Paraboloid(Point(0, 0, 0))

paramdim(::Type{<:Paraboloid}) = 2

"""
    focallength(p::Paraboloid)

Return the focal length of the paraboloid.
"""
focallength(p::Paraboloid) = p.focallength

"""
    radius(p::Paraboloid)

Return the radius of the rim of the paraboloid.
"""
radius(p::Paraboloid) = p.radius

"""
    apex(p::Paraboloid)

Return the apex of the paraboloid.
"""
apex(p::Paraboloid) = p.apex

"""
    axis(p::Paraboloid)

Return the focal axis, connecting the focus with the apex of the paraboloid.
The axis is always aligned with the z direction.
"""
function axis(p::Paraboloid)
  a = apex(p)
  f = focallength(p)
  Line(a, a + Vec(zero(f), zero(f), f))
end

==(p₁::Paraboloid, p₂::Paraboloid) =
  apex(p₁) == apex(p₂) && radius(p₁) == radius(p₂) && focallength(p₁) == focallength(p₂)

Base.isapprox(p₁::Paraboloid, p₂::Paraboloid; atol=atol(lentype(p₁)), kwargs...) =
  isapprox(apex(p₁), apex(p₂); atol, kwargs...) &&
  isapprox(radius(p₁), radius(p₂); atol, kwargs...) &&
  isapprox(focallength(p₁), focallength(p₂); atol, kwargs...)

function (p::Paraboloid)(ρ, θ)
  T = numtype(lentype(p))
  l = T(ρ) * radius(p)
  sθ, cθ = sincospi(2 * T(θ))
  x = l * cθ
  y = l * sθ
  f = focallength(p)
  z = (x^2 + y^2) / 4f
  apex(p) + Vec(x, y, z)
end
