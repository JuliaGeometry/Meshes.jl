# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FrustumSurface(bot, top)

A frustum (truncated cone) surface with `bot` and `top` disks.
See <https://en.wikipedia.org/wiki/Frustum>.

See also [`Frustum`](@ref).
"""
struct FrustumSurface{C<:CRS,D<:Disk{C}} <: Primitive{3,C}
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

==(f₁::FrustumSurface, f₂::FrustumSurface) = f₁.bot == f₂.bot && f₁.top == f₂.top

Base.isapprox(f₁::FrustumSurface, f₂::FrustumSurface; atol=atol(lentype(f₁)), kwargs...) =
  isapprox(f₁.bot, f₂.bot; atol, kwargs...) && isapprox(f₁.top, f₂.top; atol, kwargs...)

function (f::FrustumSurface)(φ, z)
  ℒ = lentype(f)
  T = numtype(ℒ)
  if (φ < 0 || φ > 1) || (z < 0 || z > 1)
    throw(DomainError((φ, z), "f(φ, z) is not defined for φ, z outside [0, 1]²."))
  end
  rb = radius(bottom(f))
  rt = radius(top(f))
  a = axis(f)
  d = a(1) - a(0)
  l = norm(d)

  # rotation to align z axis with cylinder axis
  Q = urotbetween(d, Vec(zero(ℒ), zero(ℒ), oneunit(ℒ)))

  # scale coordinates
  φₛ = 2T(π) * φ
  zₛ = z * l

  # local coordinates, that will be transformed with rotation and position of the FrustumSurface
  x = cos(φₛ) * (rb * (l - zₛ) + rt * zₛ) / l
  y = sin(φₛ) * (rb * (l - zₛ) + rt * zₛ) / l
  z = zₛ
  p = Vec(x, y, z)

  center(bottom(f)) + Q' * p
end
