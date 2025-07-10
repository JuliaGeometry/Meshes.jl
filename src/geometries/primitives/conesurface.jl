# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ConeSurface(base, apex)

A cone surface with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`Cone`](@ref).
"""
struct ConeSurface{C<:CRS,D<:Disk{C},Mₚ<:Manifold} <: Primitive{𝔼{3},C}
  base::D
  apex::Point{Mₚ,C}
end

function ConeSurface(base::Disk{C}, apex::Tuple) where {C<:Cartesian}
  coords = convert(C, Cartesian{datum(C)}(apex))
  ConeSurface(base, Point(coords))
end

paramdim(::Type{<:ConeSurface}) = 2

base(c::ConeSurface) = c.base

apex(c::ConeSurface) = c.apex

height(c::ConeSurface) = norm(center(base(c)) - apex(c))

halfangle(c::ConeSurface) = atan(radius(base(c)), height(c))

==(c₁::ConeSurface, c₂::ConeSurface) = base(c₁) == base(c₂) && apex(c₁) == apex(c₂)

Base.isapprox(c₁::ConeSurface, c₂::ConeSurface; atol=atol(lentype(c₁)), kwargs...) =
  isapprox(base(c₁), base(c₂); atol, kwargs...) && isapprox(apex(c₁), apex(c₂); atol, kwargs...)

(c::ConeSurface)(φ, h) = Cone(base(c), apex(c))(1, φ, h)
