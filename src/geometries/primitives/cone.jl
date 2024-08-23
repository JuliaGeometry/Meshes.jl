# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cone(base, apex)

A cone with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`ConeSurface`](@ref).
"""
struct Cone{C<:CRS,D<:Disk{C},Mₚ<:Manifold} <: Primitive{𝔼{3},C}
  base::D
  apex::Point{Mₚ,C}
end

function Cone(base::Disk{C}, apex::Tuple) where {C<:Cartesian}
  coords = convert(C, Cartesian{datum(C)}(apex))
  Cone(base, Point(coords))
end

paramdim(::Type{<:Cone}) = 3

base(c::Cone) = c.base

apex(c::Cone) = c.apex

height(c::Cone) = norm(center(base(c)) - apex(c))

halfangle(c::Cone) = atan(radius(base(c)), height(c))

==(c₁::Cone, c₂::Cone) = boundary(c₁) == boundary(c₂)

Base.isapprox(c₁::Cone, c₂::Cone; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(boundary(c₁), boundary(c₂); atol, kwargs...)

function (c::Cone)(φ, r, h)
  T = numtype(lentype(c))
  if (φ < 0 || φ > 1) || (r < 0 || r > 1) || (h < 0 || h > 1)
    throw(DomainError((φ, r, h), "c(φ, r, h) is not defined for φ, r, h outside [0, 1]³."))
  end

  a = c.apex
  R = radius(c.base)
  p = c.base.plane.p
  û = normalize(c.base.plane.u)
  v̂ = normalize(c.base.plane.v)
  sφ, cφ = sincospi(2φ)
  b = p + Vec(r * R * cφ * û) + Vec(r * R * sφ * v̂)
  h̄ = h * (b - a)
  a + h̄
end
