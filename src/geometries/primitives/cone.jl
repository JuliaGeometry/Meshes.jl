# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cone(base, apex)

A cone with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`ConeSurface`](@ref).
"""
struct Cone{C<:CRS,D<:Disk{C},Mₚ<:Manifold} <: Primitive{𝔼{3},C}
  base::D
  apex::Point{Mₚ,C}
end

function Cone(base::Disk{C}, apex::Tuple) where {C<:Cartesian}
  coords = convert(C, Cartesian{datum(C)}(apex))
  Cone(base, Point(coords))
end

paramdim(::Type{<:Cone}) = 3

base(c::Cone) = c.base

apex(c::Cone) = c.apex

height(c::Cone) = height(boundary(c))

halfangle(c::Cone) = halfangle(boundary(c))

==(c₁::Cone, c₂::Cone) = boundary(c₁) == boundary(c₂)

Base.isapprox(c₁::Cone, c₂::Cone; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(boundary(c₁), boundary(c₂); atol, kwargs...)

function (c::Cone)(r, φ, h)
  b = base(c)(r, φ)
  a = apex(c)
  Segment(promote(b, a)...)(h)
end
