# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Disk(plane, radius)

A disk on a `plane` with given `radius`.

See also [`Circle`](@ref).
"""
struct Disk{T} <: Primitive{3,T}
  plane::Plane{T}
  radius::T
end

Disk(plane::Plane{T}, radius) where {T} = Disk(plane, T(radius))

paramdim(::Type{<:Disk}) = 2

isconvex(::Type{<:Disk}) = true

center(d::Disk) = d.plane(0, 0)

radius(d::Disk) = d.radius

measure(d::Disk{T}) where {T} = T(π) * d.radius^2

area(d::Disk) = measure(d)

function Base.in(p::Point, d::Disk)
  p ∉ d.plane && return false
  s² = sum(abs2, p - center(d))
  r² = radius(d)^2
  s² ≤ r²
end

boundary(d::Disk) = Circle(d.plane, d.radius)

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Disk{T}}) where {T} = Disk(rand(rng, Plane{T}), rand(rng, T))
