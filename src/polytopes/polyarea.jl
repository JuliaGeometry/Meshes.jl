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
struct PolyArea{M<:AbstractManifold,C<:CRS,R<:Ring{M,C},V<:AbstractVector{R}} <: Polygon{M,C}
  rings::V

  function PolyArea{M,C,R,V}(rings; fix=true) where {M<:AbstractManifold,C<:CRS,R<:Ring{M,C},V<:AbstractVector{R}}
    if isempty(rings)
      throw(ArgumentError("cannot create PolyArea without rings"))
    end

    if fix
      outer = rings[begin]
      inners = length(rings) > 1 ? rings[(begin + 1):end] : R[]

      # fix orientation
      ofix(r, o) = orientation(r) == o ? r : reverse(r)
      outer = ofix(outer, CCW)
      inners = ofix.(inners, CW)

      # fix degeneracy
      if nvertices(outer) == 2
        A, B = vertices(outer)
        P = center(Segment(A, B))
        outer = Ring(A, P, B)
      end
      inners = filter(r -> nvertices(r) > 2, inners)

      rings = [outer; inners]
    end

    new(rings)
  end
end

PolyArea(rings::V; fix=true) where {M<:AbstractManifold,C<:CRS,R<:Ring{M,C},V<:AbstractVector{R}} =
  PolyArea{M,C,R,V}(rings; fix)

PolyArea(vertices::AbstractVector{<:AbstractVector}; fix=true) = PolyArea([Ring(v) for v in vertices]; fix)

PolyArea(outer::Ring; fix=true) = PolyArea([outer]; fix)

PolyArea(outer::AbstractVector; fix=true) = PolyArea(Ring(outer); fix)

PolyArea(outer...; fix=true) = PolyArea(collect(outer); fix)

==(p₁::PolyArea, p₂::PolyArea) = p₁.rings == p₂.rings

Base.isapprox(p₁::PolyArea, p₂::PolyArea; atol=atol(lentype(p₁)), kwargs...) =
  length(p₁.rings) == length(p₂.rings) && all(isapprox(r₁, r₂; atol, kwargs...) for (r₁, r₂) in zip(p₁.rings, p₂.rings))

vertices(p::PolyArea) = mapreduce(vertices, vcat, p.rings)

nvertices(p::PolyArea) = mapreduce(nvertices, +, p.rings)

centroid(p::PolyArea) = centroid(first(p.rings))

rings(p::PolyArea) = p.rings

function Base.unique!(p::PolyArea)
  foreach(unique!, p.rings)
  inds = findall(r -> nvertices(r) ≤ 2, p.rings)
  setdiff!(inds, 1) # don't remove outer ring
  isempty(inds) || deleteat!(p.rings, inds)
  p
end

function Base.show(io::IO, p::PolyArea)
  rings = p.rings
  print(io, "PolyArea(")
  if length(rings) == 1
    r = first(rings)
    printverts(io, vertices(r))
  else
    nverts = nvertices.(rings)
    join(io, ("$n-Ring" for n in nverts), ", ")
  end
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", p::PolyArea)
  rings = p.rings
  summary(io, p)
  println(io)
  println(io, "  outer")
  print(io, "  └─ $(rings[1])")
  if length(rings) > 1
    println(io)
    println(io, "  inner")
    printelms(io, @view(rings[2:end]), "  ")
  end
end
