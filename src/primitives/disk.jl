# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Disk(plane, radius)

A disk embedded in 3-dimensional space on a
given `plane` with given `radius`.

See also [`Circle`](@ref).
"""
struct Disk{T} <: Primitive{3,T}
  plane::Plane{T}
  radius::T
end

Disk(plane::Plane{T}, radius) where {T} = Disk(plane, T(radius))

paramdim(::Type{<:Disk}) = 2

plane(d::Disk) = d.plane

center(d::Disk) = d.plane(0, 0)

radius(d::Disk) = d.radius

normal(d::Disk) = normal(d.plane)

function (d::Disk{T})(ρ, φ) where {T}
  if (ρ < 0 || ρ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((ρ, φ), "d(ρ, φ) is not defined for ρ, φ outside [0, 1]²."))
  end
  r = d.radius
  l = T(ρ) * r
  sφ, cφ = sincospi(2 * T(φ))
  u = l * cφ
  v = l * sφ
  d.plane(u, v)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Disk{T}}) where {T} = Disk(rand(rng, Plane{T}), rand(rng, T))
