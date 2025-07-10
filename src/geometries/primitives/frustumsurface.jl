# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FrustumSurface(bot, top)

A frustum (truncated cone) surface with `bot` and `top` disks.
See <https://en.wikipedia.org/wiki/Frustum>.

See also [`Frustum`](@ref).
"""
struct FrustumSurface{C<:CRS,D<:Disk{C}} <: Primitive{𝔼{3},C}
  bot::D
  top::D

  function FrustumSurface{C,D}(bot, top) where {C<:CRS,D<:Disk{C}}
    bn = normal(plane(bot))
    tn = normal(plane(top))
    a = bn ⋅ tn
    assertion(a ≈ oneunit(a), "Bottom and top plane must be parallel")
    assertion(center(bot) ≉ center(top), "Bottom and top centers need to be distinct")
    new(bot, top)
  end
end

FrustumSurface(bot::D, top::D) where {C<:CRS,D<:Disk{C}} = FrustumSurface{C,D}(bot, top)

paramdim(::Type{<:FrustumSurface}) = 2

bottom(f::FrustumSurface) = f.bot

top(f::FrustumSurface) = f.top

height(f::FrustumSurface) = norm(center(bottom(f)) - center(top(f)))

axis(f::FrustumSurface) = Line(center(bottom(f)), center(top(f)))

==(f₁::FrustumSurface, f₂::FrustumSurface) = bottom(f₁) == bottom(f₂) && top(f₁) == top(f₂)

Base.isapprox(f₁::FrustumSurface, f₂::FrustumSurface; atol=atol(lentype(f₁)), kwargs...) =
  isapprox(bottom(f₁), bottom(f₂); atol, kwargs...) && isapprox(top(f₁), top(f₂); atol, kwargs...)

(f::FrustumSurface)(φ, z) = Frustum(bottom(f), top(f))(1, φ, z)
