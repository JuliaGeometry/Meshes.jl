# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ParaboloidSurface(point, radius, focallength)

A paraboloid surface embedded in R³ and extending up to a distance `radius` from
its focal axis, which is aligned along the z direction and passes through `point`
(the apex of the paraboloid).

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
    point::Point{3,T}
    radius::T
    focallength::T
end

function ParaboloidSurface(point::Tuple, radius, focallength)
    pt = Point(point)
    ParaboloidSurface(pt, convert(eltype(pt), radius), convert(eltype(pt), focallength))
end

ParaboloidSurface(point::Point{3,T}, radius) where {T} = ParaboloidSurface(point, T(radius), one(T))

ParaboloidSurface(point::Point{3,T}) where {T} = ParaboloidSurface(point, T(1), T(1))

ParaboloidSurface{T}() where {T} = ParaboloidSurface(Point{3, T}(0, 0, 0), T(1), T(1))

ParaboloidSurface() = ParaboloidSurface{Float64}()

@doc raw"""
    point(p::ParaboloidSurface{T})

Return the vertex of the paraboloid as a `Point{3, T}`.
"""
point(p::ParaboloidSurface) = p.point

@doc raw"""
    focallength(p::ParaboloidSurface)

Return the focal length of the paraboloid.
"""
focallength(p::ParaboloidSurface) = p.focallength

@doc raw"""
    focallength(p::ParaboloidSurface)

Return the radius of the rim of the paraboloid.
"""
radius(p::ParaboloidSurface) = p.radius

paramdim(::Type{<:ParaboloidSurface}) = 2

function area(p::ParaboloidSurface{T}) where {T}
    f = p.focallength
    r = p.radius
    8π/3 * f^2 * ((1 + r^2/(2f)^2)^(3/2) - 1)
end

measure(p::ParaboloidSurface) = area(p)

isparametrized(::Type{<:ParaboloidSurface}) = true

isperiodic(::Type{<:ParaboloidSurface}) = (false, true)

@doc raw"""
    axis(p::ParaboloidSurface)

Return the focal axis, connecting the focus with the vertex of the paraboloid.
The axis is always aligned with the z direction.
"""
axis(p::ParaboloidSurface) = Line(p.point, p.point + Vec(0, 0, p.focallength))

Base.isapprox(p₁::ParaboloidSurface{T}, p₂::ParaboloidSurface{T}) where {T} =
    (
        p₁.point ≈ p₂.point &&
        isapprox(p₁.focallength, p₂.focallength, atol=atol(T)) &&
        isapprox(p₁.radius, p₂.radius, atol=atol(T))
    )

function (p::ParaboloidSurface{T})(r, θ) where {T}
    # r is the radial coordinate, θ is the angular coordinate
    (r < 0 || r > 1) && error("radus r=$r is out of [0, 1]")
    
    f = p.focallength
    cx, cy, cz = coordinates(p.point)

    sinθ, cosθ = sincospi(T(2) * θ)
    x = r * p.radius * cosθ
    y = r * p.radius * sinθ
    z = (x^2 + y^2)/4f
    
    Point(cx + x, cy + y, cz + z)
end

discretize(p::ParaboloidSurface) = discretize(p, RegularDiscretization(30))
simplexify(p::ParaboloidSurface) = discretize(p, RegularDiscretization(30))

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ParaboloidSurface{T}}) where {T} =
  ParaboloidSurface(rand(rng, Point{3, T}), rand(rng, T), rand(rng, T))
