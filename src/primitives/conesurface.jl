# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ConeSurface(base, apex)

A cone surface with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`Cone`](@ref).
"""
struct ConeSurface{M<:AbstractManifold,C<:CRS,D<:Disk{M,C}} <: Primitive{M,C}
  base::D
  apex::Point{M,C}
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

function (c::ConeSurface)(φ, h)
  T = numtype(lentype(c))
  if (φ < 0 || φ > 1) || (h < 0 || h > 1)
    throw(DomainError((φ, h), "c(φ, h) is not defined for φ, h outside [0, 1]²."))
  end
  n = -normal(c.base)
  v = c.base(T(0), T(0)) - c.apex
  l = norm(v)
  θ = ∠(n, v)
  o = c.apex + T(h) * v
  r = T(h) * l * cos(θ)
  s = Circle(Plane(o, n), r)
  s(T(φ))
end
