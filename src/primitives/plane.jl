# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Plane(o, n)

A plane coincident with point `p`, with normal vector `n`.
"""
struct Plane{Dim,T} <: Primitive{Dim,T}
  p::Point{Dim,T}
  n::Vec{Dim,T}
end

paramdim(::Type{<:Plane}) = 2