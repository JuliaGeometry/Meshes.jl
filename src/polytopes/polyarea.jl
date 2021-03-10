# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PolyArea(outer, [inner1, inner2, ..., innerk])

A polygonal area with `outer` chain, and optional inner
chains `inner1`, `inner2`, ..., `innerk`.

Chains can be a vector of [`Point`](@ref) or a
vector of tuples with coordinates for convenience.

Most algorithms assume that the outer chain is oriented
counter-clockwise (CCW) and that all inner chains are
oriented clockwise (CW).
"""
struct PolyArea{Dim,T,C<:Chain{Dim,T}} <: Polygon{Dim,T}
  outer::C
  inners::Vector{C}

  function PolyArea{Dim,T,C}(outer, inners) where {Dim,T,C}
    @assert isclosed(outer) "invalid outer chain"
    @assert all(isclosed.(inners)) "invalid inner chains"

    # fix orientation if needed
    fix(c, o) = orientation(c) == o ? c : reverse(c)
    ofixed = fix(outer, :CCW)
    ifixed = map(c -> fix(c, :CW), inners)

    new(ofixed, ifixed)
  end
end

PolyArea(outer::C, inners=[]) where {Dim,T,C<:Chain{Dim,T}} =
  PolyArea{Dim,T,Chain{Dim,T}}(outer, inners)

PolyArea(outer::AbstractVector{P}, inners=[]) where {P<:Point} =
  PolyArea(Chain(outer), [Chain(inner) for inner in inners])

PolyArea(outer::AbstractVector{TP}, inners=[]) where {TP<:Tuple} =
  PolyArea(Point.(outer), [Point.(inner) for inner in inners])

PolyArea(outer::Vararg{P}) where {P<:Point} = PolyArea(collect(outer))

PolyArea(outer::Vararg{TP}) where {TP<:Tuple} = PolyArea(collect(Point.(outer)))

==(p1::PolyArea, p2::PolyArea) =
  p1.outer == p2.outer && p1.inners == p2.inners

nvertices(p::PolyArea) = nvertices(p.outer) + mapreduce(nvertices, +, p.inners, init=0)

centroid(p::PolyArea) = centroid(p.outer)

"""
    chains(polyarea)

Return the outer and inner chains of the polygon.
"""
chains(p::PolyArea) = [p.outer; p.inners]

"""
    hasholes(polyarea)

Tells whether or not the `polyarea` contains holes.
"""
hasholes(p::PolyArea) = !isempty(p.inners)

"""
    issimple(polyarea)

Tells whether or not the `polyarea` is a simple polygon.
See https://en.wikipedia.org/wiki/Simple_polygon.
"""
issimple(p::PolyArea) = !hasholes(p) && issimple(p.outer)

"""
    windingnumber(point, polyarea)

Winding number of `point` with respect to the `polyarea`.
"""
windingnumber(point::Point, p::PolyArea) =
  windingnumber(point, p.outer)

"""
    orientation(polyarea)

Returns the orientation of the `polyarea` as either
counter-clockwise (CCW) or clockwise (CW).

For polygons with holes, returns a list of orientations.
"""
orientation(p::PolyArea, algo=WindingOrientation()) =
  orientation.([p.outer; p.inners], Ref(algo))

"""
    unique(polyarea)

Return a new `polyarea` without duplicate vertices.
"""
Base.unique(p::PolyArea) = unique!(deepcopy(p))

"""
    unique!(polyarea)

Remove duplicate vertices in `polyarea`.
"""
function Base.unique!(p::PolyArea)
  close!(unique!(open!(p.outer)))
  hasholes(p) && foreach(c->close!(unique!(open!(c))), p.inners)
  p
end

"""
    bridge(polyarea)

Transform `polyarea` with holes into a single
outer chain via bridges.

## References

* Held. 1998. [FIST: Fast Industrial-Strength Triangulation of Polygons]
  (https://link.springer.com/article/10.1007/s00453-001-0028-4)
"""
function bridge(p::PolyArea)
  hasholes(p) || return first(chains(p))

  # retrieve chains with coordinates
  pchains = [coordinates.(vertices(c)) for c in chains(p)]

  # sort vertices lexicographically
  coords  = [coord for pchain in pchains for coord in pchain]
  indices = sortperm(sortperm(coords))

  # each chain has its own set of indices
  pinds = Vector{Int}[]; offset = 0
  for nvertex in length.(pchains)
    push!(pinds, indices[offset+1:offset+nvertex])
    offset += nvertex
  end

  # sort chains based on leftmost vertex
  leftmost = argmin.(pinds)
  minimums = getindex.(pinds, leftmost)
  reorder  = sortperm(minimums)
  leftmost = leftmost[reorder]
  minimums = minimums[reorder]
  pchains  = pchains[reorder]
  pinds    = pinds[reorder]

  # initialize outer boundary
  outer = first(pchains)
  oinds = first(pinds)

  # merge islands into outer boundary
  for i in 2:length(pchains)
    l = leftmost[i]
    m = minimums[i]
    c = pchains[i]
    o = pinds[i]

    # find closest vertex in boundary
    dmin, jmin = Inf, 0
    for j in findall(oinds .≤ m)
      d = sum(abs, outer[j] - c[l])
      if d < dmin
        dmin, jmin = d, j
      end
    end

    # insert island at closest vertex
    island = push!(circshift(c, -l+1), c[l])
    iinds  = push!(circshift(o, -l+1), o[l])
    outer = [outer[1:jmin]; island; outer[jmin:end]]
    oinds = [oinds[1:jmin]; iinds;  oinds[jmin:end]]
  end

  # close boundary
  push!(outer, first(outer))

  Chain(Point.(outer))
end

function Base.show(io::IO, p::PolyArea)
  outer = p.outer
  inner = isempty(p.inners) ? "" : ", "*join(p.inners, ", ")
  print(io, "PolyArea($outer$inner)")
end

function Base.show(io::IO, ::MIME"text/plain", p::PolyArea{Dim,T}) where {Dim,T}
  outer = "    └─$(p.outer)"
  inner = ["    └─$v" for v in p.inners]
  println(io, "PolyArea{$Dim,$T}")
  println(io, "  outer")
  if isempty(inner)
    print(io, outer)
  else
    println(io, outer)
    println(io, "  inner")
    print(io, join(inner, "\n"))
  end
end
