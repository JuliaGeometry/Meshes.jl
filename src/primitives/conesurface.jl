# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ConeSurface(base, apex)

A cone surface with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`Cone`](@ref).
"""
struct ConeSurface{D<:Disk,P<:Point{3}} <: Primitive{3}
  base::D
  apex::P
end

ConeSurface(base::Disk, apex::Tuple) = ConeSurface(base, Point(apex))

paramdim(::Type{<:ConeSurface}) = 2

lentype(::Type{<:ConeSurface{D}}) where {D} = lentype(D)

base(c::ConeSurface) = c.base

apex(c::ConeSurface) = c.apex

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

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ConeSurface}) =
  ConeSurface(rand(rng, Disk), rand(rng, Point{3}))
