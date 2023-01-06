# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cone(disk, apex)

A cone with base `disk` and `apex`.
See https://en.wikipedia.org/wiki/Cone.

See also [`ConeSurface`](@ref).
"""
struct Cone{T} <: Primitive{3,T}
  disk::Disk{T}
  apex::Point{3,T}
end

Cone(disk::Disk, apex::Tuple) = Cone(disk, Point(apex))

paramdim(::Type{<:Cone}) = 3

isconvex(::Type{<:Cone}) = true

boundary(c::Cone) = ConeSurface(c.disk, c.apex)
