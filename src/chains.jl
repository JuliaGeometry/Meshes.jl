# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Chain(p1, p2, ..., pn)

A polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.
See https://en.wikipedia.org/wiki/Polygonal_chain.
"""
struct Chain{Dim,T} <: Polytope{1,Dim,T}
  vertices::Vector{Point{Dim,T}}
end

Chain(points::Vararg{P}) where {P<:Point} = Chain(collect(points))
Chain(points::AbstractVector{TP}) where {TP<:Tuple} = Chain(Point.(points))
Chain(points::Vararg{TP}) where {TP<:Tuple} = Chain(collect(points))

"""
    isclosed(chain)

Tells whether or not the chain is closed.

A closed chain is also known as a ring.
"""
isclosed(c::Chain) = first(c.vertices) == last(c.vertices)

"""
   issimple(chain)

Tells whether or not the `chain` is simple.

A chain is simple when all its segments only
intersect at end points.
"""
function issimple(c::Chain)
  vs = c.vertices
  ss = [Segment(view(vs, [i,i+1])) for i in 1:length(vs)-1]
  for i in 1:length(ss)
    for j in i+1:length(ss)
      I = intersecttype(ss[i], ss[j])
      if !(I isa CornerTouchingSegments || I isa NonIntersectingSegments)
        return false
      end
    end
  end
  true
end

"""
    windingnumber(point, chain)

Winding number of `point` with respect to the `chain`.

## References

* Balbes, R. and Siegel, J. 1990. [A robust method for calculating
  the simplicity and orientation of planar polygons]
  (https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
function windingnumber(p::Point, c::Chain)
  vₒ, vs = p, c.vertices
  sum(∠(vs[i], vₒ, vs[i+1]) for i in 1:length(vs)-1)
end

"""
    orientation(chain)

Returns the orientation of the `chain` as either
counter-clockwise (CCW) or clockwise (CW).
"""
function orientation(c::Chain{Dim,T}) where {Dim,T}
  # pick any segment
  x1, x2 = c.vertices[1:2]
  x̄ = center(Segment(x1, x2))
  w = windingnumber(x̄, c) - ∠(x1, x̄, x2)
  isapprox(w, π, atol=atol(T)) ? :CCW : :CW
end

function Base.show(io::IO, c::Chain)
  N = length(c.vertices)
  print(io, "$N-chain")
end

function Base.show(io::IO, ::MIME"text/plain", c::Chain{Dim,T}) where {Dim,T}
  N = length(c.vertices)
  lines = ["  └─$v" for v in c.vertices]
  println(io, "$N-chain{$Dim,$T}")
  print(io, join(lines, "\n"))
end
