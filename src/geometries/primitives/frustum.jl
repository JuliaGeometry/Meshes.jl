# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Frustum(bot, top)

A frustum (truncated cone) with `bot` and `top` disks.
See <https://en.wikipedia.org/wiki/Frustum>.

See also [`FrustumSurface`](@ref).
"""
struct Frustum{C<:CRS,D<:Disk{C}} <: Primitive{𝔼{3},C}
  bot::D
  top::D

  function Frustum{C,D}(bot, top) where {C<:CRS,D<:Disk{C}}
    bn = normal(plane(bot))
    tn = normal(plane(top))
    a = bn ⋅ tn
    assertion(a ≈ oneunit(a), "Bottom and top plane must be parallel")
    assertion(center(bot) ≉ center(top), "Bottom and top centers need to be distinct")
    new(bot, top)
  end
end

Frustum(bot::D, top::D) where {C<:CRS,D<:Disk{C}} = Frustum{C,D}(bot, top)

paramdim(::Type{<:Frustum}) = 3

bottom(f::Frustum) = f.bot

top(f::Frustum) = f.top

height(f::Frustum) = height(boundary(f))

axis(f::Frustum) = axis(boundary(f))

==(f₁::Frustum, f₂::Frustum) = boundary(f₁) == boundary(f₂)

Base.isapprox(f₁::Frustum, f₂::Frustum; atol=atol(lentype(f₁)), kwargs...) =
  isapprox(boundary(f₁), boundary(f₂); atol, kwargs...)

function (frustum::Frustum)(r, φ, h)
  if (r < 0 || r > 1) || (φ < 0 || φ > 1) || (h < 0 || h > 1)
    throw(DomainError((r, φ, h), "frustum(r, φ, h) is not defined for r, φ, h outside [0, 1]³."))
  end
  a = bottom(frustum)(r, φ)
  b = top(frustum)(r, φ)
  Segment(a, b)(h)
end
