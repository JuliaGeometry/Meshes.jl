# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Torus(center, normal, major, minor)

A torus centered at `center` with axis of revolution directed by 
`normal` and with radii `major` and `minor`. 

"""
struct Torus{T} <: Primitive{3,T}
  center::Point{3,T}
  normal::Vec{3,T}
  major::T
  minor::T
end

function Torus(center::Tuple, normal::Tuple, major::T1, minor::T2) where {T1,T2}
  T = promote_type(eltype(center), eltype(normal), T1, T2, Float32)  
  Torus(Point(T.(center)), Vec(T.(normal)), T(major), T(minor))
end

paramdim(::Type{<:Torus}) = 2

isconvex(::Type{<:Torus}) = false

isperiodic(::Type{<:Torus}) = (true, true)

center(t::Torus) = t.center

radii(t::Torus) = (t.major, t.minor)

axis(t::Torus) = Line(t.center, t.center + t.normal)

# https://en.wikipedia.org/wiki/Torus
function measure(t::Torus)
  R, r = t.major, t.minor
  4π^2 * R * r
end

area(t::Torus) = measure(t)

function Base.in(p::Point, t::Torus)
  R, r = radii(t)
  c = center(t)
  n⃗ = t.normal
  rotation = Rotate(n⃗, Vec(0, 0, 1))
  M = convert(DCM, rotation.rot)
  x, y, z = M * (p - c)
  (R - √(x^2 + y^2))^2 + z^2 ≤ r^2
end
