# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ConeSurface(base, apex)

A cone surface with `base` disk and `apex`.
See https://en.wikipedia.org/wiki/Cone.

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

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ConeSurface{T}}) where {T} =
  ConeSurface(rand(rng, Disk{T}), rand(rng, Point{3,T}))
