"""
    HyperSphere{N, T}

A `HyperSphere` is a generalization of a sphere into N-dimensions.
A `center` and `radius` must be specified.
"""
struct HyperSphere{N,T} <: Primitive{N,T}
    center::Point{N,T}
    radius::T
end

origin(s::HyperSphere) = s.center
radius(s::HyperSphere) = s.radius

# TODO: review these
widths(s::HyperSphere{N,T}) where {N,T} = vfill(Vec{N,T}, 2s.radius)
Base.minimum(s::HyperSphere) = coordinates(s.center) .- s.radius
Base.maximum(s::HyperSphere) = coordinates(s.center) .+ s.radius

"""
    Circle{T}

An alias for a HyperSphere of dimension 2. (i.e. `HyperSphere{2, T}`)
"""
const Circle{T} = HyperSphere{2,T}

"""
    Sphere{T}

An alias for a HyperSphere of dimension 3. (i.e. `HyperSphere{3, T}`)
"""
const Sphere{T} = HyperSphere{3,T}

function Base.in(p::AbstractPoint, s::HyperSphere)
    x = coordinates(p)
    c = coordinates(s.center)
    r = s.radius
    sum(abs2, x - c) ≤ r^2
end

function coordinates(s::Circle, nvertices=64)
    o = coordinates(s.center)
    r = s.radius
    φ = LinRange(0, 2pi, nvertices)
    inner(φ) = Vec(r*sin(φ+pi), r*cos(φ+pi)) + o
    return (inner(φ) for φ in φ)
end

function coordinates(s::Sphere, nvertices=24)
    o = coordinates(s.center)
    r = s.radius
    θ = LinRange(0, pi, nvertices)
    φ = LinRange(0, 2pi, nvertices)
    inner(θ, φ) = Vec(r*cos(φ)*sin(θ), r*sin(φ)*sin(θ), r*cos(θ)) + o
    return ivec((inner(θ, φ) for θ in θ, φ in φ))
end

function faces(::Sphere{T}, nvertices=24) where {T}
    return faces(Rectangle(Point{2,T}(0,0), Vec{2,T}(1,1)), (nvertices, nvertices))
end
