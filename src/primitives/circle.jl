# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Circle(plane, radius)

A circle on a `plane` with given `radius`.

See also [`Disk`](@ref).
"""
struct Circle{T} <: Primitive{3,T}
  plane::Plane{T}
  radius::T
end

Circle(plane::Plane{T}, radius) where {T} = Circle(plane, T(radius))

paramdim(::Type{<:Circle}) = 1

isconvex(::Type{<:Circle}) = false

isperiodic(::Type{<:Circle}) = (true,)

center(c::Circle) = origin(c.plane)

radius(c::Circle) = c.radius

measure(c::Circle{T}) where {T} = 2 * T(π) * c.radius

Base.length(c::Circle) = measure(c)

function Base.in(p::Point{3,T}, c::Circle{T}) where {T}
  p ∉ c.plane && return false
  s² = sum(abs2, p - center(c))
  r² = radius(c)^2
  isapprox(s², r², atol = atol(T)^2)
end

boundary(::Circle) = nothing