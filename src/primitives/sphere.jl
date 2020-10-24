# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Sphere(center, radius)

A sphere with `center` and `radius`.
"""
struct Sphere{Dim,T} <: Primitive{Dim,T}
    center::Point{Dim,T}
    radius::T
end

center(s::Sphere) = s.center
radius(s::Sphere) = s.radius

function Base.in(p, s::Sphere)
    x = coordinates(p)
    c = coordinates(s.center)
    r = s.radius
    sum(abs2, x - c) ≈ r^2
end

function coordinates(s::Sphere{2}, nvertices=64)
    φ = LinRange(0, 2π, nvertices)
    c = coordinates(s.center)
    r = s.radius
    inner(φ) = Vec(r*sin(φ+pi), r*cos(φ+pi)) + c
    ivec((inner(φ) for φ in φ))
end

function coordinates(s::Sphere{3}, nvertices=24)
    φ = LinRange(0, 2π, nvertices)
    θ = LinRange(0, π, nvertices)
    c = coordinates(s.center)
    r = s.radius
    inner(θ, φ) = Vec(r*cos(φ)*sin(θ), r*sin(φ)*sin(θ), r*cos(θ)) + c
    ivec((inner(θ, φ) for θ in θ, φ in φ))
end

function faces(::Sphere{3,T}, nvertices=24) where {T}
    faces(Box(Point{2,T}(0,0), Point{2,T}(1,1)), (nvertices, nvertices))
end
