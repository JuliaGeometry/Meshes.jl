"""
    Sphere(center, radius)

A sphere with `center` and `radius`.
"""
struct Sphere{N,T} <: Primitive{N,T}
    center::Point{N,T}
    radius::T
end

origin(s::Sphere) = s.center
radius(s::Sphere) = s.radius

# TODO: review these
widths(s::Sphere{N,T}) where {N,T} = vfill(Vec{N,T}, 2s.radius)
Base.minimum(s::Sphere) = coordinates(s.center) .- s.radius
Base.maximum(s::Sphere) = coordinates(s.center) .+ s.radius

function Base.in(p, s::Sphere)
    x = coordinates(p)
    c = coordinates(s.center)
    r = s.radius
    sum(abs2, x - c) ≤ r^2
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
    faces(Rectangle(Point{2,T}(0,0), Vec{2,T}(1,1)), (nvertices, nvertices))
end
