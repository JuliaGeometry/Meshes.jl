# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ngon(p1, p2, ..., pN)

A N-gon is a polygon with `N` vertices `p1`, `p2`, ..., `pN`
oriented counter-clockwise (CCW). In this case the number of
vertices is fixed and known at compile time. Examples of N-gon
are `Triangle` (N=3), `Quadrangle` (N=4), `Pentagon` (N=5), etc.

### Notes

Although the number of vertices `N` is known at compile time,
we use abstract vectors to store the list of vertices. This
design allows constructing N-gon from views of global vectors
without expensive memory allocations.
"""
struct Ngon{N,Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polygon{Dim,T}
  vertices::V
end

Ngon{N}(vertices::AbstractVector{Point{Dim,T}}) where {N,Dim,T} =
  Ngon{N,Dim,T,typeof(vertices)}(vertices)

Ngon(vertices::AbstractVector{Point{Dim,T}}) where {Dim,T} =
  Ngon{length(vertices)}(vertices)

# type aliases for convenience
const Triangle   = Ngon{3}
const Quadrangle = Ngon{4}
const Pentagon   = Ngon{5}
const Hexagon    = Ngon{6}
const Heptagon   = Ngon{7}
const Octagon    = Ngon{8}
const Nonagon    = Ngon{9}
const Decagon    = Ngon{10}

# number of vertices is known at compile time
nvertices(::Type{<:Ngon{N}}) where {N} = N
nvertices(ngon::Ngon) = nvertices(typeof(ngon))

function signarea(t::Triangle{2})
  vs = t.vertices
  signarea(vs[1], vs[2], vs[3])
end

measure(t::Triangle{2}) = abs(signarea(t))

function measure(q::Quadrangle)
  vs = q.vertices
  Δ₁ = Triangle(view(vs, [1,2,3]))
  Δ₂ = Triangle(view(vs, [3,4,1]))
  measure(Δ₁) + measure(Δ₂)
end

function Base.in(p::Point{2}, t::Triangle{2})
  a, b, c = t.vertices
  abp = signarea(a, b, p)
  bcp = signarea(b, c, p)
  cap = signarea(c, a, p)
  areas = (abp, bcp, cap)
  all(areas .≥ 0) || all(areas .≤ 0)
end

function Base.in(p::Point, q::Quadrangle)
  vs = q.vertices
  Δ₁ = Triangle(view(vs, [1,2,3]))
  Δ₂ = Triangle(view(vs, [3,4,1]))
  p ∈ Δ₁ || p ∈ Δ₂
end
