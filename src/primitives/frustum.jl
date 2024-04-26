# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Frustum(bot, top)

A frustum (truncated cone) with `bot` and `top` disks.
See <https://en.wikipedia.org/wiki/Frustum>.

See also [`FrustumSurface`](@ref).
"""
struct Frustum{D<:Disk} <: Primitive{3}
  bot::D
  top::D

  function Frustum{D}(bot, top) where {D}
    bn = normal(plane(bot))
    tn = normal(plane(top))
    a = bn ⋅ tn
    @assert a ≈ oneunit(a) "Bottom and top plane must be parallel"
    @assert center(bot) ≉ center(top) "Bottom and top centers need to be distinct"
    new(bot, top)
  end
end

Frustum(bot::D, top::D) where {D<:Disk} = Frustum{D}(bot, top)

paramdim(::Type{<:Frustum}) = 3

lentype(::Type{<:Frustum{D}}) where {D} = lentype(D)

bottom(f::Frustum) = f.bot

top(f::Frustum) = f.top

height(f::Frustum) = height(boundary(f))

axis(f::Frustum) = axis(boundary(f))

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Frustum})
  bottom = rand(rng, Disk)
  ax = normal(plane(bottom))
  topplane = Plane(center(bottom) + rand() * ax, ax)
  top = Disk(topplane, rand(Met{Float64}))
  Frustum(bottom, top)
end
