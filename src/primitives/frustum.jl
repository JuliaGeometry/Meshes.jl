# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Frustum(bot, top)

A frustum (truncated cone) with `bot` and `top` disks.
See <https://en.wikipedia.org/wiki/Frustum>.

See also [`FrustumSurface`](@ref).
"""
struct Frustum{D<:Disk} <: Primitive
  bot::D
  top::D

  function Frustum{D}(bot, top) where {D}
    bn = normal(plane(bot))
    tn = normal(plane(top))
    @assert bn ⋅ tn ≈ 1 "Bottom and top plane must be parallel"
    @assert center(bot) ≉ center(top) "Bottom and top centers need to be distinct"
    new(bot, top)
  end
end

Frustum(bot::D, top::D) where {D<:Disk} = Frustum{D}(bot, top)

paramdim(::Type{<:Frustum}) = 3

bottom(f::Frustum) = f.bot

top(f::Frustum) = f.top

height(f::Frustum) = height(boundary(f))

axis(f::Frustum) = axis(boundary(f))

# TODO
# function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Frustum{T}}) where {T}
#   bottom = rand(rng, Disk{T})
#   ax = normal(plane(bottom))
#   topplane = Plane{T}(center(bottom) + rand(T) * ax, ax)
#   top = Disk{T}(topplane, rand(T))
#   Frustum(bottom, top)
# end
