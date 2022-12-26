# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  Cone(vertex, center, radius, nsides, closed)

Cone with `vertex` and base with `center` and `radius` made of `nsides` 
segments. If `closed`, the base is included.
"""
struct Cone{T,I<:Integer} <: Primitive{3,T}
  vertex::Point{3,T}
  center::Point{3,T}
  radius::T
  nsides::I
  closed::Bool
end

function Cone(vertex::Tuple, center::Tuple, radius, nsides, closed)
  T = promote_type(eltype(vertex), eltype(center), typeof(radius))
  Cone(Point{3,T}(vertex), Point{3,T}(center), T(radius), nsides, closed)
end

paramdim(::Type{<:Cone}) = 3

isconvex(c::Cone) = c.closed

axis(c::Cone) = Line(c.vertex, c.center)

height(c::Cone) = norm(c.vertex - c.center)