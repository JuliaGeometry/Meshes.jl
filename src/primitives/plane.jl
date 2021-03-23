# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Plane(p₀, v, w)

A plane coincident with point `p₀`, defined by non-parallel vectors `v` and `w`.
It can be called as p(s, t) with numeric parameters `s` and `t` to cast it at
`p₀ + s*v + t*w`.
"""
struct Plane{Dim,T} <: Primitive{Dim,T}
  p₀::Point{Dim,T}
  v::Vec{Dim,T}
  w::Vec{Dim,T}
end

function (p::Plane)(s, t)
  p.p₀ + s*p.v + t*p.w
end

paramdim(::Type{<:Plane}) = 2