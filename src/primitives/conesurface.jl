# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ConeSurface(disk, apex)

A cone surface with base `disk` and `apex`.
See https://en.wikipedia.org/wiki/Cone.

See also [`Cone`](@ref).
"""
struct ConeSurface{T} <: Primitive{3,T}
  disk::Disk{T}
  apex::Point{3,T}
end

ConeSurface(disk::Disk{T}, apex::Point{3,T}) where {T} = ConeSurface{T}(disk, apex)

ConeSurface(disk::Disk, apex::Tuple) = ConeSurface(disk, Point(apex))

paramdim(::Type{<:ConeSurface}) = 2

isconvex(::Type{<:ConeSurface}) = false

boundary(::ConeSurface) = nothing
