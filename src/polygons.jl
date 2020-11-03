# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Polygon(outer, [inner1, inner2, ..., innerk])

A polygon with `outer` ring, and optional inner
rings `inner1`, `inner2`, ..., `innerk`.

Rings can be a vector of [`Point`](@ref) or a
vector of tuples with coordinates for convenience.

Most algorithms assume that the outer ring is oriented
counter-clockwise (CCW) and that all inner rings are
oriented clockwise (CW).
"""
struct Polygon{Dim,T,C<:Chain{Dim,T}} <: Polytope{Dim,T}
  outer::C
  inners::Vector{C}

  function Polygon{Dim,T,C}(outer, inners) where {Dim,T,C}
    @assert isclosed(outer) "invalid outer ring"
    @assert all(isclosed.(inners)) "invalid inner rings"
    new(outer, inners)
  end
end

Polygon(outer::C, inners=[]) where {Dim,T,C<:Chain{Dim,T}} =
  Polygon{Dim,T,Chain{Dim,T}}(outer, inners)

Polygon(outer::AbstractVector{P}, inners=[]) where {P<:Point} =
  Polygon(Chain(outer), [Chain(inner) for inner in inners])

Polygon(outer::AbstractVector{TP}, inners=[]) where {TP<:Tuple} =
  Polygon(Point.(outer), [Point.(inner) for inner in inners])

Polygon(outer::Vararg{P}) where {P<:Point} = Polygon(collect(outer))

Polygon(outer::Vararg{TP}) where {TP<:Tuple} = Polygon(collect(Point.(outer)))

"""
    rings(polygon)

Return the outer and inner rings of the polygon.
"""
rings(p::Polygon) = p.outer, p.inners

"""
    hasholes(polygon)

Tells whether or not the `polygon` contains holes.
"""
hasholes(p::Polygon) = !isempty(p.inners)

"""
    windingnumber(point, polygon)

Winding number of `point` with respect to the `polygon`.

## References

* Balbes, R. and Siegel, J. 1990. [A robust method for calculating
  the simplicity and orientation of planar polygons]
  (https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
function windingnumber(point::Point, polygon::Polygon)
  xₒ = point
  xs = vertices(polygon.outer)
  sum(∠(xs[i], xₒ, xs[i+1]) for i in 1:length(xs)-1)
end

function Base.show(io::IO, p::Polygon)
  outer = p.outer
  inner = isempty(p.inners) ? "" : ", "*join(p.inners, ", ")
  print(io, "Polygon($outer$inner)")
end

function Base.show(io::IO, ::MIME"text/plain", p::Polygon{Dim,T}) where {Dim,T}
  outer = "    └─$(p.outer)"
  inner = ["    └─$v" for v in p.inners]
  println(io, "Polygon{$Dim,$T}")
  println(io, "  outer")
  if isempty(inner)
    print(io, outer)
  else
    println(io, outer)
    println(io, "  inner")
    print(io, join(inner, "\n"))
  end
end
