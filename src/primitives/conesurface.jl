# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ConeSurface(base, apex)

A cone surface with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`Cone`](@ref).
"""
struct ConeSurface{T} <: Primitive{3,T}
  base::Disk{T}
  apex::Point{3,T}
end

ConeSurface(base::Disk, apex::Tuple) = ConeSurface(base, Point(apex))

paramdim(::Type{<:ConeSurface}) = 2

base(c::ConeSurface) = c.base

apex(c::ConeSurface) = c.apex

function (c::ConeSurface{T})(φ, h) where {T}
  if (φ < 0 || φ > 1) || (h < 0 || h > 1)
    throw(DomainError((φ, h), "c(φ, h) is not defined for φ, h outside [0, 1]²."))
  end
  n = -normal(plane(c.base))
  v = c.base(0, 0) - c.apex
  l = norm(v)
  θ = ∠(n, v)
  o = c.apex + h * v
  r = h * l * cos(θ)
  s = Circle(Plane(o, n), r)
  s(φ)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ConeSurface{T}}) where {T} =
  ConeSurface(rand(rng, Disk{T}), rand(rng, Point{3,T}))
