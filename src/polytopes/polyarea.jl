# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PolyArea(outer, [inner1, inner2, ..., innerk]; fix=true)

A polygonal area with `outer` ring, and optional inner
rings `inner1`, `inner2`, ..., `innerk`.

Rings can be a vector of [`Point`](@ref) or a
vector of tuples with coordinates for convenience,
in which case the first point should *not* be repeated
at the end of the vector.

The option `fix` tries to correct issues with polygons
in the real world, including issues with:

* `orientation` - Most algorithms assume that the
  outer ring is oriented counter-clockwise (CCW) and
  that all inner rings are oriented clockwise (CW).

* `degeneracy` - Sometimes data is shared with
  degenerate rings (i.e. only 2 vertices).
"""
struct PolyArea{Dim,T,R<:Ring{Dim,T}} <: Polygon{Dim,T}
  outer::R
  inners::Vector{R}

  function PolyArea{Dim,T,R}(outer, inners, fix) where {Dim,T,R}
    if fix
      # fix orientation
      ofix(r, o) = orientation(r) == o ? r : reverse(r)
      outer = ofix(outer, :CCW)
      inners = ofix.(inners, :CW)

      # fix degeneracy
      if nvertices(outer) == 2
        v = vertices(outer)
        A, B = v[1], v[2]
        M = centroid(Segment(A, B))
        outer = Ring(A, M, B)
      end
      inners = filter(c -> nvertices(c) > 2, inners)
    end

    new(outer, inners)
  end
end

PolyArea(outer::R, inners=R[]; fix=true) where {Dim,T,V,R<:Ring{Dim,T,V}} = PolyArea{Dim,T,R}(outer, inners, fix)

PolyArea(outer::AbstractVector{P}, inners=[]; fix=true) where {P<:Point} =
  PolyArea(Ring(outer), [Ring(inner) for inner in inners]; fix=fix)

PolyArea(outer::AbstractVector{TP}, inners=[]; fix=true) where {TP<:Tuple} =
  PolyArea(Point.(outer), [Point.(inner) for inner in inners]; fix=fix)

PolyArea(outer::Vararg{P}; fix=true) where {P<:Point} = PolyArea(collect(outer); fix=fix)

PolyArea(outer::Vararg{TP}; fix=true) where {TP<:Tuple} = PolyArea(collect(Point.(outer)); fix=fix)

==(p1::PolyArea, p2::PolyArea) = p1.outer == p2.outer && p1.inners == p2.inners

function vertices(p::PolyArea{Dim,T}) where {Dim,T}
  vo = vertices(p.outer)
  vi = reduce(vcat, vertices(inner) for inner in p.inners; init=Point{Dim,T}[])
  [vo; vi]
end

nvertices(p::PolyArea) = nvertices(p.outer) + mapreduce(nvertices, +, p.inners, init=0)

centroid(p::PolyArea) = centroid(p.outer)

rings(p::PolyArea) = [p.outer; p.inners]

hasholes(p::PolyArea) = !isempty(p.inners)

issimple(p::PolyArea) = !hasholes(p) && issimple(p.outer)

windingnumber(point::Point, p::PolyArea) = windingnumber(point, p.outer)

function Base.unique!(p::PolyArea)
  unique!(p.outer)
  hasholes(p) && foreach(c -> unique!(c), p.inners)
  p
end

function Base.in(point::Point, polyarea::PolyArea)
  sideof(point, polyarea.outer) == :INSIDE && all(sideof(point, inner) == :OUTSIDE for inner in polyarea.inners)
end

function Base.show(io::IO, p::PolyArea)
  nverts = nvertices.([p.outer; p.inners])
  rings = join(["$n-Ring" for n in nverts], ", ")
  print(io, "PolyArea($rings)")
end

function Base.show(io::IO, ::MIME"text/plain", p::PolyArea{Dim,T}) where {Dim,T}
  nverts = nvertices.([p.outer; p.inners])
  rings = ["$n-Ring" for n in nverts]
  println(io, "PolyArea{$Dim,$T}")
  if length(rings) == 1
    println(io, "  outer")
    print(io, io_lines(rings[1:1], "    "))
  else
    println(io, "  outer")
    println(io, io_lines(rings[1:1], "    "))
    println(io, "  inner")
    print(io, io_lines(rings[2:end], "    "))
  end
end
