# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ParaboloidSurface(point, radius, focallength)

A paraboloid surface embedded in R³ and extending up to a distance
`radius` from its focal axis, which is aligned along the z direction
and passes through `point` (the apex of the paraboloid). The equation
of the paraboloid is the following:

```math
f(x, y) = \\frac{(x - x_0)^2 + (y - y_0)^2}{4f} + z_0\\qquad\\text{for } x^2 + y^2 < r^2,
```

where ``(x_0, y_0, z_0)`` is the vertex of the parabola, ``f`` is the
focal length, and ``r`` is the clip radius.

    ParaboloidSurface(point, radius)

This creates a paraboloid surface with focal length equal to 1.

    ParaboloidSurface(point)

This creates a paraboloid surface with focal length equal to 1 and a rim with unit
radius.

    ParaboloidSurface()

Same as above, but here the vertex is at `Point(0, 0, 0)`.

See also <https://en.wikipedia.org/wiki/Paraboloid>.
"""
struct ParaboloidSurface{T} <: Primitive{3, T}
    vertex::Point{3,T}
    radius::T
    focallength::T
end

ParaboloidSurface(point::Tuple, radius, focallength) =
    ParaboloidSurface(Point(point), radius, focallength)

ParaboloidSurface(point::Point{3,T}, radius) where {T} = ParaboloidSurface(point, T(radius), T(1))

ParaboloidSurface(vertex::Point{3,T}) where {T} = ParaboloidSurface(vertex, T(1), T(1))

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
    axis(p::ParaboloidSurface)

Return the focal axis, connecting the focus with the vertex of the paraboloid.
The axis is always aligned with the z direction.
"""
axis(p::ParaboloidSurface) = Line(p.vertex, p.vertex + Vec(0, 0, p.focallength))

Base.isapprox(p₁::ParaboloidSurface{T}, p₂::ParaboloidSurface{T}) where {T} =
    (
        p₁.vertex ≈ p₂.vertex &&
        isapprox(p₁.focallength, p₂.focallength, atol=atol(T)) &&
        isapprox(p₁.radius, p₂.radius, atol=atol(T))
    )

function (p::ParaboloidSurface{T})(r, θ) where {T}
    # r is the radial coordinate, θ is the angular coordinate
    (r < 0 || r > 1) && throw(DomainError((r, θ), "radius r=$r is out of [0, 1]"))
    
    f = p.focallength
    cx, cy, cz = coordinates(p.vertex)

    sinθ, cosθ = sincospi(T(2) * θ)
    x = r * p.radius * cosθ
    y = r * p.radius * sinθ
    z = (x^2 + y^2)/4f
    
    Point(cx + x, cy + y, cz + z)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ParaboloidSurface{T}}) where {T} =
  ParaboloidSurface(rand(rng, Point{3, T}), rand(rng, T), rand(rng, T))
