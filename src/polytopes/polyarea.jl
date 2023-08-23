# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PolyArea(outer; fix=true)
    PolyArea([outer, inner₁, inner₂, ..., innerₖ]; fix=true)

A polygonal area with `outer` ring, and optional inner
rings `inner₁`, `inner₂`, ..., `innerₖ`.

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
  degenerate rings (e.g. only 2 vertices).
"""
struct PolyArea{Dim,T,R<:Ring{Dim,T}} <: Polygon{Dim,T}
  rings::Vector{R}

  function PolyArea{Dim,T,R}(rings::Vector{R}) where {Dim,T,R<:Ring{Dim,T}}
    if isempty(rings)
      throw(ArgumentError("cannot create PolyArea without rings"))
    end
    new(rings)
  end
end

function PolyArea(rings::AbstractVector{R}; fix=true) where {Dim,T,R<:Ring{Dim,T}}
  N = length(rings)
  if fix && N > 0
    outer = rings[begin]
    inners = N > 1 ? rings[(begin + 1):end] : R[]

    # fix orientation
    ofix(r, o) = orientation(r) == o ? r : reverse(r)
    outer = ofix(outer, CCW)
    inners = ofix.(inners, CW)

    # fix degeneracy
    if nvertices(outer) == 2
      v = vertices(outer)
      A, B = v[1], v[2]
      M = centroid(Segment(A, B))
      outer = Ring(A, M, B)
    end
    inners = filter(r -> nvertices(r) > 2, inners)

    rings = [outer; inners]
  end

  PolyArea{Dim,T,R}(rings)
end

PolyArea(vertices::AbstractVector{<:AbstractVector}; fix=true) = PolyArea([Ring(v) for v in vertices]; fix)

PolyArea(outer::Ring; fix=true) = PolyArea([outer]; fix)

PolyArea(outer::AbstractVector; fix=true) = PolyArea(Ring(outer); fix)

PolyArea(outer...; fix=true)  = PolyArea(collect(outer); fix)

==(p₁::PolyArea, p₂::PolyArea) = p₁.rings == p₂.rings

function Base.isapprox(p₁::PolyArea, p₂::PolyArea; kwargs...)
  length(p₁.rings) ≠ length(p₂.rings) && return false
  all(isapprox(r₁, r₂; kwargs...) for (r₁, r₂) in zip(p₁.rings, p₂.rings))
end

vertices(p::PolyArea) = [vertices(r) for r in p.rings]

nvertices(p::PolyArea) = mapreduce(nvertices, +, p.rings)

centroid(p::PolyArea) = centroid(first(p.rings))

rings(p::PolyArea) = p.rings

windingnumber(point::Point, p::PolyArea) = windingnumber(point, first(p.rings))

function Base.unique!(p::PolyArea)
  foreach(unique!, p.rings)
  inds = findall(r -> nvertices(r) ≤ 2, p.rings)
  # don't remove first ring (outer)
  setdiff!(inds, first(eachindex(p.rings)))
  isempty(inds) || deleteat!(p.rings, inds)
  p
end

function Base.show(io::IO, p::PolyArea)
  nverts = nvertices.(p.rings)
  rings = join(["$n-Ring" for n in nverts], ", ")
  print(io, "PolyArea($rings)")
end

function Base.show(io::IO, ::MIME"text/plain", p::PolyArea{Dim,T}) where {Dim,T}
  nverts = nvertices.(p.rings)
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

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{<:PolyArea{Dim,T}}) where {Dim,T} =
  PolyArea(rand(rng, Ring{Dim,T}))
