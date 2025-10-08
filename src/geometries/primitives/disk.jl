# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Disk(plane, radius)

A disk embedded in 3-dimensional space on a
given `plane` with given `radius`.

See also [`Circle`](@ref).
"""
struct Disk{C<:CRS,P<:Plane{C},ℒ<:Len} <: Primitive{𝔼{3},C}
  plane::P
  radius::ℒ
  Disk(plane::P, radius::ℒ) where {C<:CRS,P<:Plane{C},ℒ<:Len} = new{C,P,float(ℒ)}(plane, radius)
end

Disk(plane::Plane, radius) = Disk(plane, addunit(radius, u"m"))

paramdim(::Type{<:Disk}) = 2

plane(d::Disk) = d.plane

radius(d::Disk) = d.radius

center(d::Disk) = center(boundary(d))

normal(d::Disk) = normal(plane(d))

==(d₁::Disk, d₂::Disk) = boundary(d₁) == boundary(d₂)

Base.isapprox(d₁::Disk, d₂::Disk; atol=atol(lentype(d₁)), kwargs...) =
  isapprox(boundary(d₁), boundary(d₂); atol, kwargs...)

function (d::Disk)(ρ, φ)
  if (ρ < 0 || ρ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((ρ, φ), "d(ρ, φ) is not defined for ρ, φ outside [0, 1]²."))
  end
  L = lentype(d)
  units = unit(L)
  T = numtype(L)
  l = T(ρ) * radius(d)
  sφ, cφ = sincospi(2 * T(φ))
  u = ustrip(units, l * cφ)
  v = ustrip(units, l * sφ)
  plane(d)(u, v)
end
