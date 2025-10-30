# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ParaboloidSurface(apex, radius, focallength)

A paraboloid surface embedded in R³ and extending up to a distance
`radius` from its focal axis, which is aligned along the z direction
and passes through `apex` (the apex of the paraboloid). The equation
of the paraboloid is the following:

```math
f(x, y) = \\frac{(x - x_0)^2 + (y - y_0)^2}{4f} + z_0\\qquad\\text{for } x^2 + y^2 < r^2,
```

where ``(x_0, y_0, z_0)`` is the apex of the parabola, ``f`` is the
focal length, and ``r`` is the clip radius.

    ParaboloidSurface(apex, radius)

This creates a paraboloid surface with focal length equal to 1.

    ParaboloidSurface(apex)

This creates a paraboloid surface with focal length equal to 1 and a rim with unit
radius.

    ParaboloidSurface()

Same as above, but here the apex is at `Apex(0, 0, 0)`.

See also <https://en.wikipedia.org/wiki/Paraboloid>.
"""
struct ParaboloidSurface{C<:CRS,Mₚ<:Manifold,ℒ<:Len} <: Primitive{𝔼{3},C}
  apex::Point{Mₚ,C}
  radius::ℒ
  focallength::ℒ
  ParaboloidSurface(apex::Point{Mₚ,C}, radius::ℒ, focallength::ℒ) where {C<:CRS,Mₚ<:Manifold,ℒ<:Len} =
    new{C,Mₚ,float(ℒ)}(apex, radius, focallength)
end

ParaboloidSurface(apex::Point, radius::Len, focallength::Len) = ParaboloidSurface(apex, promote(radius, focallength)...)

ParaboloidSurface(apex::Point, radius, focallength) =
  ParaboloidSurface(apex, aslen(radius), aslen(focallength))

ParaboloidSurface(apex::Tuple, radius, focallength) = ParaboloidSurface(Point(apex), radius, focallength)

ParaboloidSurface(apex::Point, radius) = ParaboloidSurface(apex, radius, oneunit(radius))

ParaboloidSurface(apex::Tuple, radius) = ParaboloidSurface(Point(apex), radius)

ParaboloidSurface(apex::Point) = ParaboloidSurface(apex, oneunit(lentype(apex)))

ParaboloidSurface(apex::Tuple) = ParaboloidSurface(Point(apex))

ParaboloidSurface() = ParaboloidSurface(Point(0, 0, 0))

paramdim(::Type{<:ParaboloidSurface}) = 2

"""
    focallength(p::ParaboloidSurface)

Return the focal length of the paraboloid.
"""
focallength(p::ParaboloidSurface) = p.focallength

"""
    focallength(p::ParaboloidSurface)

Return the radius of the rim of the paraboloid.
"""
radius(p::ParaboloidSurface) = p.radius

"""
    apex(p::ParaboloidSurface)

Return the apex of the paraboloid.
"""
apex(p::ParaboloidSurface) = p.apex

"""
    axis(p::ParaboloidSurface)

Return the focal axis, connecting the focus with the apex of the paraboloid.
The axis is always aligned with the z direction.
"""
function axis(p::ParaboloidSurface)
  a = apex(p)
  f = focallength(p)
  Line(a, a + Vec(zero(f), zero(f), f))
end

==(p₁::ParaboloidSurface, p₂::ParaboloidSurface) =
  apex(p₁) == apex(p₂) && radius(p₁) == radius(p₂) && focallength(p₁) == focallength(p₂)

Base.isapprox(p₁::ParaboloidSurface, p₂::ParaboloidSurface; atol=atol(lentype(p₁)), kwargs...) =
  isapprox(apex(p₁), apex(p₂); atol, kwargs...) &&
  isapprox(radius(p₁), radius(p₂); atol, kwargs...) &&
  isapprox(focallength(p₁), focallength(p₂); atol, kwargs...)

function (p::ParaboloidSurface)(ρ, θ)
  if (ρ < 0 || ρ > 1)
    throw(DomainError((ρ, θ), "p(ρ, θ) is not defined for ρ outside [0, 1]."))
  end
  T = numtype(lentype(p))
  l = T(ρ) * radius(p)
  sθ, cθ = sincospi(2 * T(θ))
  x = l * cθ
  y = l * sθ
  f = focallength(p)
  z = (x^2 + y^2) / 4f
  apex(p) + Vec(x, y, z)
end
