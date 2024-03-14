# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cone(base, apex)

A cone with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`ConeSurface`](@ref).
"""
struct Cone{T} <: Primitive{3,T}
  base::Disk{T}
  apex::Point{3,T}
end

Cone(base::Disk, apex::Tuple) = Cone(base, Point(apex))

paramdim(::Type{<:Cone}) = 3

base(c::Cone) = c.base

apex(c::Cone) = c.apex

height(c::Cone) = norm(center(base(c)) - apex(c))

halfangle(c::Cone) = atan(radius(base(c)), height(c))

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Cone{T}}) where {T} =
  Cone(rand(rng, Disk{T}), rand(rng, Point{3,T}))
