# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PolyArea(outer, [inner1, inner2, ..., innerk]; fix=true)

A polygonal area with `outer` chain, and optional inner
chains `inner1`, `inner2`, ..., `innerk`.

Chains can be a vector of [`Point`](@ref) or a
vector of tuples with coordinates for convenience.

The option `fix` tries to correct issues with polygons
in the real world, including issues with:

* `orientation` - Most algorithms assume that the
  outer ring is oriented counter-clockwise (CCW) and
  that all inner rings are oriented clockwise (CW).

* `degeneracy` - Sometimes data is shared with
  degenerate rings (i.e. only 2 vertices).
"""
struct PolyArea{Dim,T,C<:Chain{Dim,T}} <: Polygon{Dim,T}
  outer::C
  inners::Vector{C}

  function PolyArea{Dim,T,C}(outer, inners, fix) where {Dim,T,C}
    @assert isclosed(outer) "invalid outer chain"
    @assert all(isclosed, inners) "invalid inner chains"

    if fix
      # fix orientation
      ofix(c, o) = orientation(c) == o ? c : reverse(c)
      outer = ofix(outer, :CCW)
      inners = ofix.(inners, :CW)

      # fix degeneracy
      if nvertices(outer) == 2
        v = vertices(outer)
        A, B = v[1], v[2]
        s = Segment(A, B)
        M = centroid(s)
        outer = Chain([A, M, B, A])
      end
      inners = filter(c -> nvertices(c) > 2, inners)
    end

    new(outer, inners)
  end
end

PolyArea(outer::C, inners=C[]; fix=true) where {Dim,T,V,C<:Chain{Dim,T,V}} =
  PolyArea{Dim,T,Chain{Dim,T,V}}(outer, inners, fix)

PolyArea(outer::AbstractVector{P}, inners=[]; fix=true) where {P<:Point} =
  PolyArea(Chain(outer), [Chain(inner) for inner in inners]; fix=fix)

PolyArea(outer::AbstractVector{TP}, inners=[]; fix=true) where {TP<:Tuple} =
  PolyArea(Point.(outer), [Point.(inner) for inner in inners]; fix=fix)

PolyArea(outer::Vararg{P}; fix=true) where {P<:Point} =
  PolyArea(collect(outer); fix=fix)

PolyArea(outer::Vararg{TP}; fix=true) where {TP<:Tuple} =
  PolyArea(collect(Point.(outer)); fix=fix)

==(p1::PolyArea, p2::PolyArea) =
  p1.outer == p2.outer && p1.inners == p2.inners

function vertices(p::PolyArea{Dim,T}) where {Dim,T}
  vo = vertices(p.outer)
  vi = reduce(vcat, vertices(inner) for inner in p.inners; init=Point{Dim,T}[])
  [vo; vi]
end

nvertices(p::PolyArea) = nvertices(p.outer) + mapreduce(nvertices, +, p.inners, init=0)

centroid(p::PolyArea) = centroid(p.outer)

chains(p::PolyArea) = [p.outer; p.inners]

hasholes(p::PolyArea) = !isempty(p.inners)

issimple(p::PolyArea) = !hasholes(p) && issimple(p.outer)

windingnumber(point::Point, p::PolyArea) =
  windingnumber(point, p.outer)

function Base.unique!(p::PolyArea)
  close!(unique!(open!(p.outer)))
  hasholes(p) && foreach(c -> close!(unique!(open!(c))), p.inners)
  p
end

function Base.in(point::Point, polyarea::PolyArea)
  sideof(point, polyarea.outer) == :INSIDE &&
    all(sideof(point, inner) == :OUTSIDE for inner in polyarea.inners)
end

function Base.show(io::IO, p::PolyArea)
  nverts = [[npoints(p.outer)]; npoints.(p.inners)]
  chains = join(["$n-chain" for n in nverts], ", ")
  print(io, "PolyArea($chains)")
end

function Base.show(io::IO, ::MIME"text/plain", p::PolyArea{Dim,T}) where {Dim,T}
  nverts = [[npoints(p.outer)]; npoints.(p.inners)]
  chains = ["$n-chain" for n in nverts]
  println(io, "PolyArea{$Dim,$T}")
  if length(chains) == 1
    println(io, "  outer")
    print(io, io_lines(chains[1:1], "    "))
  else
    println(io, "  outer")
    println(io, io_lines(chains[1:1], "    "))
    println(io, "  inner")
    print(io, io_lines(chains[2:end], "    "))
  end
end
