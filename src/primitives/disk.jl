# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Disk(plane, radius)

A disk embedded in 3-dimensional space on a
given `plane` with given `radius`.

See also [`Circle`](@ref).
"""
struct Disk{CRS,P<:Plane{CRS},ℒ<:Len} <: Primitive{3,CRS}
  plane::P
  radius::ℒ
  Disk{CRS,P,ℒ}(plane, radius) where {CRS,P<:Plane{CRS},ℒ<:Len} = new(plane, radius)
end

Disk(plane::P, radius::ℒ) where {CRS,P<:Plane{CRS},ℒ<:Len} = Disk{CRS,P,float(ℒ)}(plane, radius)

Disk(plane::Plane, radius) = Disk(plane, addunit(radius, u"m"))

paramdim(::Type{<:Disk}) = 2

plane(d::Disk) = d.plane

center(d::Disk) = d.plane(0, 0)

radius(d::Disk) = d.radius

normal(d::Disk) = normal(d.plane)

function (d::Disk)(ρ, φ)
  T = numtype(lentype(d))
  if (ρ < 0 || ρ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((ρ, φ), "d(ρ, φ) is not defined for ρ, φ outside [0, 1]²."))
  end
  r = d.radius
  l = T(ρ) * r
  sφ, cφ = sincospi(2 * T(φ))
  u = ustrip(l * cφ)
  v = ustrip(l * sφ)
  d.plane(u, v)
end

Random.rand(rng::Random.AbstractRNG, ::Type{Disk}) = Disk(rand(rng, Plane), rand(rng, Met{Float64}))
