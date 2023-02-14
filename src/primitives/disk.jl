# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Disk(plane, radius)

A disk on a `plane` with given `center` and `radius`.

See also [`Circle`](@ref).
"""
struct Disk{T} <: Primitive{3,T}
  plane::Plane{T}
  center::Point{3,T}
  radius::T
end

Disk(plane::Plane{T}, center::Point{3,T}, radius) where {T} = 
  Disk(plane, center, T(radius))

Disk(plane::Plane{T}, center::Tuple, radius) where {T} = 
  Disk(plane, Point(center), radius)

paramdim(::Type{<:Disk}) = 2

isconvex(::Type{<:Disk}) = true

center(d::Disk) = d.center

radius(d::Disk) = d.radius

measure(d::Disk{T}) where {T} = T(π) * d.radius^2

area(d::Disk) = measure(d)

function Base.in(p::Point, d::Disk)
  p ∉ d.plane && return false
  s² = sum(abs2, p - center(d))
  r² = radius(d)^2
  s² ≤ r²
end

boundary(d::Disk) = Circle(d.plane, d.center, d.radius)