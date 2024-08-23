# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ConeSurface(base, apex)

A cone surface with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`Cone`](@ref).
"""
struct ConeSurface{C<:CRS,D<:Disk{C},Mₚ<:Manifold} <: Primitive{𝔼{3},C}
  base::D
  apex::Point{Mₚ,C}
end

function ConeSurface(base::Disk{C}, apex::Tuple) where {C<:Cartesian}
  coords = convert(C, Cartesian{datum(C)}(apex))
  ConeSurface(base, Point(coords))
end

paramdim(::Type{<:ConeSurface}) = 2

base(c::ConeSurface) = c.base

apex(c::ConeSurface) = c.apex

==(c₁::ConeSurface, c₂::ConeSurface) = c₁.base == c₂.base && c₁.apex == c₂.apex

Base.isapprox(c₁::ConeSurface, c₂::ConeSurface; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(c₁.base, c₂.base; atol, kwargs...) && isapprox(c₁.apex, c₂.apex; atol, kwargs...)

function (conesurface::ConeSurface)(uφ, uh)
  T = numtype(lentype(cone))
  if (uφ < 0 || uφ > 1) || (uh < 0 || uh > 1)
    throw(DomainError((uφ, uh), "c(φ, h) is not defined for φ, h outside [0, 1]²."))
  end

  # Aliases
  a = cone.apex
  R = cone.base.radius
  b = cone.base.plane.p
  û = cone.base.plane.u
  v̂ = cone.base.plane.v

  # Scaled parametric coords
  sφ, cφ = sincospi(2uφ)

  # Locate parametric point
  c = b + (R * cφ * û) + (R * sφ * v̂)
  h̄ = uh * (c - a)
  a + h̄
end
