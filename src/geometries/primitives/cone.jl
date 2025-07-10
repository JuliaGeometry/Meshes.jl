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

height(c::Cone) = height(boundary(c))

halfangle(c::Cone) = halfangle(boundary(c))

==(c₁::Cone, c₂::Cone) = boundary(c₁) == boundary(c₂)

Base.isapprox(c₁::Cone, c₂::Cone; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(boundary(c₁), boundary(c₂); atol, kwargs...)

function (c::Cone)(r, φ, h)
  if (r < 0 || r > 1) || (φ < 0 || φ > 1) || (h < 0 || h > 1)
    throw(DomainError((r, φ, h), "c(r, φ, h) is not defined for r, φ, h outside [0, 1]³."))
  end
  b = base(c)(r, φ)
  a = apex(c)
  Segment(b, a)(h)
end
