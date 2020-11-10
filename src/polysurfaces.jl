# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PolySurface(outer, [inner1, inner2, ..., innerk])

A polygonal surface with `outer` ring, and optional inner
rings `inner1`, `inner2`, ..., `innerk`.

Rings can be a vector of [`Point`](@ref) or a
vector of tuples with coordinates for convenience.

Most algorithms assume that the outer ring is oriented
counter-clockwise (CCW) and that all inner rings are
oriented clockwise (CW).
"""
struct PolySurface{Dim,T,C<:Chain{Dim,T}} <: Polytope{2,Dim,T}
  outer::C
  inners::Vector{C}

  function PolySurface{Dim,T,C}(outer, inners) where {Dim,T,C}
    @assert isclosed(outer) "invalid outer ring"
    @assert all(isclosed.(inners)) "invalid inner rings"
    new(outer, inners)
  end
end

PolySurface(outer::C, inners=[]) where {Dim,T,C<:Chain{Dim,T}} =
  PolySurface{Dim,T,Chain{Dim,T}}(outer, inners)

PolySurface(outer::AbstractVector{P}, inners=[]) where {P<:Point} =
  PolySurface(Chain(outer), [Chain(inner) for inner in inners])

PolySurface(outer::AbstractVector{TP}, inners=[]) where {TP<:Tuple} =
  PolySurface(Point.(outer), [Point.(inner) for inner in inners])

PolySurface(outer::Vararg{P}) where {P<:Point} = PolySurface(collect(outer))

PolySurface(outer::Vararg{TP}) where {TP<:Tuple} = PolySurface(collect(Point.(outer)))

"""
    rings(polysurface)

Return the outer and inner rings of the polygon.
"""
rings(p::PolySurface) = p.outer, p.inners

"""
    hasholes(polysurface)

Tells whether or not the `polysurface` contains holes.
"""
hasholes(p::PolySurface) = !isempty(p.inners)

"""
    issimple(polysurface)

Tells whether or not the `polysurface` is a simple polygon.
See https://en.wikipedia.org/wiki/Simple_polygon.
"""
issimple(p::PolySurface) = !hasholes(p) && issimple(p.outer)

"""
    windingnumber(point, polysurface)

Winding number of `point` with respect to the `polysurface`.
"""
windingnumber(point::Point, p::PolySurface) =
  windingnumber(point, p.outer)

"""
    orientation(polysurface)

Returns the orientation of the `polysurface` as either
counter-clockwise (CCW) or clockwise (CW).

For polygons with holes, returns a list of orientations.
"""
function orientation(p::PolySurface)
  if hasholes(p)
    orientation(p.outer), orientation.(p.inners)
  else
    orientation(p.outer)
  end
end

function Base.show(io::IO, p::PolySurface)
  outer = p.outer
  inner = isempty(p.inners) ? "" : ", "*join(p.inners, ", ")
  print(io, "PolySurface($outer$inner)")
end

function Base.show(io::IO, ::MIME"text/plain", p::PolySurface{Dim,T}) where {Dim,T}
  outer = "    └─$(p.outer)"
  inner = ["    └─$v" for v in p.inners]
  println(io, "PolySurface{$Dim,$T}")
  println(io, "  outer")
  if isempty(inner)
    print(io, outer)
  else
    println(io, outer)
    println(io, "  inner")
    print(io, join(inner, "\n"))
  end
end
