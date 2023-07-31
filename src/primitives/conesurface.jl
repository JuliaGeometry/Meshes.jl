# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ConeSurface(disk, apex)

A cone surface with base `disk` and `apex`.
See https://en.wikipedia.org/wiki/Cone.

See also [`Cone`](@ref).
"""
struct ConeSurface{T} <: Primitive{3,T}
  disk::Disk{T}
  apex::Point{3,T}
end

ConeSurface(disk::Disk, apex::Tuple) = ConeSurface(disk, Point(apex))

paramdim(::Type{<:ConeSurface}) = 2

boundary(::ConeSurface) = nothing

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ConeSurface{T}}) where {T} =
  ConeSurface(rand(rng, Disk{T}), rand(rng, Point{3,T}))
