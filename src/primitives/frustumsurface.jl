# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  FrustumSurface(bot, top)

A frustum (truncated cone) surface with `bot` and `top` disks.
See https://en.wikipedia.org/wiki/Frustum.

See also [`Frustum`](@ref).
"""
struct FrustumSurface{T} <: Primitive{3,T}
  bot::Disk{T}
  top::Disk{T}

  function FrustumSurface{T}(bot, top) where {T}
    bn = normal(plane(bot))
    tn = normal(plane(top))
    @assert bn ⋅ tn ≈ 1 "Bottom and top plane must be parallel"
    @assert center(bot) ≉  center(top) "Bottom and top centers need to be distinct"
    new(bot, top)
  end
end

FrustumSurface(bot::Disk{T}, top::Disk{T}) where {T} = FrustumSurface{T}(bot, top)

bottom(f::FrustumSurface) = f.bot

top(f::FrustumSurface) = f.top

height(f::FrustumSurface) = norm(center(bottom(f)) - center(top(f)))

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{FrustumSurface{T}}) where {T}
  bottom = rand(rng, Disk{T})
  ax = normal(plane(bottom))
  topplane = Plane{T}(center(bottom)+rand(T)*ax, ax)
  top = Disk{T}(topplane, rand(T))
  FrustumSurface(bottom, top)
end
