# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PolySurface(outer, [inner1, inner2, ..., innerk])

A polygonal surface with `outer` chain, and optional inner
chains `inner1`, `inner2`, ..., `innerk`.

Chains can be a vector of [`Point`](@ref) or a
vector of tuples with coordinates for convenience.

Most algorithms assume that the outer chain is oriented
counter-clockwise (CCW) and that all inner chains are
oriented clockwise (CW).
"""
struct PolySurface{Dim,T,C<:Chain{Dim,T}} <: Polygon{Dim,T}
  outer::C
  inners::Vector{C}

  function PolySurface{Dim,T,C}(outer, inners) where {Dim,T,C}
    @assert isclosed(outer) "invalid outer chain"
    @assert all(isclosed.(inners)) "invalid inner chains"
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

==(p1::PolySurface, p2::PolySurface) =
  p1.outer == p2.outer && p1.inners == p2.inners

"""
    chains(polysurface)

Return the outer and inner chains of the polygon.
"""
chains(p::PolySurface) = [p.outer; p.inners]

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
  orientation.([p.outer; p.inners])
end

"""
    unique(polysurface)

Return a new `polysurface` without duplicate vertices.
"""
Base.unique(p::PolySurface) = unique!(deepcopy(p))

"""
    unique!(polysurface)

Remove duplicate vertices in `polysurface`.
"""
function Base.unique!(p::PolySurface)
  close!(unique!(open!(p.outer)))
  hasholes(p) && foreach(c->close!(unique!(open!(c))), p.inners)
  p
end

"""
    oriented!(polysurface)

Fix orientation of `polysurface` so that outer
chain is counter-clockwise (CCW) and inner chains
are clockwise (CW).
"""
function oriented!(p::PolySurface)
  orients = orientation(p)
  first(orients) == :CCW || reverse!(p.outer)
  for i in 2:length(orients)
    orients[i] == :CW || reverse!(p.inners[i])
  end
  p
end

"""
    bridge(polysurface)

Transform `polysurface` with holes into a single
outer chain via bridges.

## References

* Held. 1998. [FIST: Fast Industrial-Strength Triangulation of Polygons]
  (https://link.springer.com/article/10.1007/s00453-001-0028-4)
"""
function bridge(p::PolySurface)
  # retrieve chains with coordinates
  pchains = [coordinates.(vertices(c)[begin:end-1]) for c in chains(p)]

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

  PolySurface(Point.(outer))
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
