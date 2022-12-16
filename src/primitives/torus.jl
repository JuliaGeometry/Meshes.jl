# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Torus(majorRadius, minorRadius)

A torus with `majorRadius` and `minorRadius` with the z axis as its axis of 
rotation and the xy plane as its plane of reflection. 

"""
struct Torus{T} <: Primitive{3,T}
  majorRadius::T
  minorRadius::T
end

function Torus(majorRadius::T1, minorRadius::T2) where {T1,T2}
  T = promote(T1, T2)  
  Torus(T(majorRadius), T(minorRadius))
end

paramdim(::Type{<:Torus}) = 2

isconvex(::Type{<:Torus}) = false

isperiodic(::Type{<:Torus}) = (true, true)

minorRadius(t::Torus) = t.minorRadius

majorRadius(t::Torus) = t.majorRadius

# https://en.wikipedia.org/wiki/Torus
function measure(t::Torus)
  R, r = t.majorRadius, t.minorRadius
  2π^2 * R * r^2
end

Base.length(t::Torus) = measure(t)

volume(t::Torus) = measure(t)

boundary(::Torus) = nothing

perimeter(::Torus{T}) where {T} = zero(T)

#= to have such a function for a general torus (I have in mind the torus passing 
  through three points), we need more than the two radii
function Base.in(p::Point, s::Sphere)
  x = coordinates(p)
  c = coordinates(s.center)
  r = s.radius
  sum(abs2, x - c) ≈ r^2
end
 =#