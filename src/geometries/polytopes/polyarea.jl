# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PolyArea(outer)
    PolyArea([outer, inner₁, inner₂, ..., innerₖ])

A polygonal area with `outer` ring, and optional inner
rings `inner₁`, `inner₂`, ..., `innerₖ`.

Rings can be a vector of [`Point`](@ref) or a
vector of tuples with coordinates for convenience,
in which case the first point should *not* be repeated
at the end of the vector.
"""
struct PolyArea{M<:Manifold,C<:CRS,R<:Ring{M,C},V<:AbstractVector{R}} <: Polygon{M,C}
  rings::V
end

PolyArea(vertices::AbstractVector{<:AbstractVector}) = PolyArea([Ring(v) for v in vertices])

PolyArea(outer::Ring) = PolyArea([outer])

PolyArea(outer::AbstractVector) = PolyArea(Ring(outer))

PolyArea(outer...) = PolyArea(collect(outer))

==(p₁::PolyArea, p₂::PolyArea) = p₁.rings == p₂.rings

Base.isapprox(p₁::PolyArea, p₂::PolyArea; atol=atol(lentype(p₁)), kwargs...) =
  length(p₁.rings) == length(p₂.rings) && all(isapprox(r₁, r₂; atol, kwargs...) for (r₁, r₂) in zip(p₁.rings, p₂.rings))

function vertex(p::PolyArea, ind)
  offset = 0
  for r in p.rings
    nverts = nvertices(r)
    if ind ≤ offset + nverts
      return vertex(r, ind - offset)
    end
    offset += nverts
  end
  throw(BoundsError(p, ind))
end

vertices(p::PolyArea) = collect(eachvertex(p))

nvertices(p::PolyArea) = mapreduce(nvertices, +, p.rings)

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
    printelms(io, @view(rings[2:end]), tab="  ")
  end
end
