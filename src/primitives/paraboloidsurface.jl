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
struct ParaboloidSurface{T} <: Primitive{3, T}
    apex::Point{3,T}
    radius::T
    focallength::T
end

ParaboloidSurface(apex::Tuple, radius, focallength) =
    ParaboloidSurface(Point(apex), radius, focallength)

ParaboloidSurface(apex::Point{3,T}, radius) where {T} = ParaboloidSurface(apex, T(radius), T(1))

ParaboloidSurface(apex::Point{3,T}) where {T} = ParaboloidSurface(apex, T(1), T(1))

ParaboloidSurface{T}() where {T} = ParaboloidSurface(Point{3, T}(0, 0, 0), T(1), T(1))

ParaboloidSurface() = ParaboloidSurface{Float64}()

paramdim(::Type{<:ParaboloidSurface}) = 2

isparametrized(::Type{<:ParaboloidSurface}) = true

isperiodic(::Type{<:ParaboloidSurface}) = (false, true)

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
axis(p::ParaboloidSurface) = Line(p.apex, p.apex + Vec(0, 0, p.focallength))

Base.isapprox(p₁::ParaboloidSurface{T}, p₂::ParaboloidSurface{T}) where {T} =
  p₁.apex ≈ p₂.apex &&
  isapprox(p₁.focallength, p₂.focallength, atol=atol(T)) &&
  isapprox(p₁.radius, p₂.radius, atol=atol(T))

function (p::ParaboloidSurface{T})(r, θ) where {T}
    # r is the radial coordinate, θ is the angular coordinate
  if (r < 0 || r > 1)
    throw(DomainError((r, θ), "radius r=$r is out of [0, 1]"))
  end
    
    f = p.focallength
    cx, cy, cz = coordinates(p.apex)

    sinθ, cosθ = sincospi(T(2) * θ)
    x = r * p.radius * cosθ
    y = r * p.radius * sinθ
    z = (x^2 + y^2)/4f
    
    Point(cx + x, cy + y, cz + z)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ParaboloidSurface{T}}) where {T} =
  ParaboloidSurface(rand(rng, Point{3, T}), rand(rng, T), rand(rng, T))
