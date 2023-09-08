# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TruncatedCone(bot, top)

A truncated cone (frustum) with `bot` and `top` disks.
See https://en.wikipedia.org/wiki/Frustum.

"""
struct TruncatedCone{T} <: Primitive{3,T}
  bot::Disk{T}
  top::Disk{T}

  function TruncatedCone{T}(bot, top) where {T}
    bn = normal(plane(bot))
    tn = normal(plane(top))
    @assert bn ⋅ tn ≈ 1 "Bottom and top plane must be parallel"
    @assert center(bot) ≉  center(top) "Bottom and top centers need to be distinct"
    new(bot, top)
  end
end

TruncatedCone(bot::Disk{T}, top::Disk{T}) where {T} = TruncatedCone{T}(bot, top)

paramdim(::Type{<:TruncatedCone}) = 3

bottom(c::TruncatedCone) = c.bot

top(c::TruncatedCone) = c.top

height(c::TruncatedCone) = norm(center(bottom(c)) - center(top(c)))

halfangle(c::TruncatedCone) = atan(radius(base(c)), height(c))

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{TruncatedCone{T}}) where {T}
  bottom = rand(rng, Disk{T})
  ax = normal(plane(bottom))
  topplane = Plane{T}(center(bottom)+rand(T)*ax, ax)
  top = Disk{T}(topplane, rand(T))
  TruncatedCone(bottom, top)
end
