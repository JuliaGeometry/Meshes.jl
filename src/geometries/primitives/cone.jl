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

function (cone::Cone)(uφ, ur, uh)
  T = numtype(lentype(cone))
  if (uφ < 0 || uφ > 1) || (ur < 0 || ur > 1) || (uh < 0 || uh > 1)
    throw(DomainError((uφ, ur, uh), "c(φ, r, h) is not defined for φ, r, h outside [0, 1]³."))
  end

  # Aliases
  a = cone.apex
  R = cone.base.radius
  b = cone.base.plane.p
  û = normalize(cone.base.plane.u)
  v̂ = normalize(cone.base.plane.v)

  # Scaled parametric coords
  sφ, cφ = sincospi(2uφ)
  r = R * ur

  # Locate parametric point
  c = b + Vec(r * cφ * û) + Vec(r * sφ * v̂)
  h̄ = uh * (c - a)
  a + h̄
end
