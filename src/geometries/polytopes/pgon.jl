# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Pgon(outer)
    Pgon([outer, inner₁, inner₂, ..., innerₖ])

A polygonal area with `outer` ring oriented counter-clockwise,
and optional inner rings `inner₁`, `inner₂`, ..., `innerₖ`
oriented clockwise.

Rings can be specified either as a vector of [`Point`](@ref)s or,
for convenience, as a vector of tuples of coordinates. The first
point should *not* be repeated at the end of the vector. Use the
[`orientation`](@ref) function to verify that the rings satisfy
the orientation requirements described above.

The [`Repair`](@ref) transform can be used to correct the orientation
of rings in polygonal areas that have already been constructed.
"""
struct Pgon{M<:Manifold,C<:CRS,R<:Ring{M,C},V<:AbstractVector{R}} <: Polygon{M,C}
  rings::V
end

Pgon(vertices::AbstractVector{<:AbstractVector}) = Pgon([Ring(v) for v in vertices])

Pgon(outer::Ring) = Pgon([outer])

Pgon(outer::AbstractVector) = Pgon(Ring(outer))

Pgon(outer...) = Pgon(collect(outer))

==(p₁::Pgon, p₂::Pgon) = p₁.rings == p₂.rings

Base.isapprox(p₁::Pgon, p₂::Pgon; atol=atol(lentype(p₁)), kwargs...) =
  length(p₁.rings) == length(p₂.rings) && all(isapprox(r₁, r₂; atol, kwargs...) for (r₁, r₂) in zip(p₁.rings, p₂.rings))

function vertex(p::Pgon, ind)
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

vertices(p::Pgon) = collect(eachvertex(p))

nvertices(p::Pgon) = mapreduce(nvertices, +, p.rings)

rings(p::Pgon) = p.rings

normal(p::Pgon) = newellnormal(vertices(first(p.rings)))

function Base.unique!(p::Pgon)
  foreach(unique!, p.rings)
  inds = findall(r -> nvertices(r) ≤ 2, p.rings)
  setdiff!(inds, 1) # don't remove outer ring
  isempty(inds) || deleteat!(p.rings, inds)
  p
end

function Base.show(io::IO, p::Pgon)
  rings = p.rings
  print(io, "Pgon(")
  if length(rings) == 1
    r = first(rings)
    printverts(io, vertices(r))
  else
    nverts = map(nvertices, rings)
    join(io, ("$n-Ring" for n in nverts), ", ")
  end
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", p::Pgon)
  rings = p.rings
  summary(io, p)
  println(io)
  if length(rings) > 1
    println(io, "  outer ring")
    printelms(io, vertices(rings[1]), tab="    ")
    for i in 2:length(rings)
      println(io)
      println(io, "  inner ring $(i-1)")
      printelms(io, vertices(rings[i]), tab="    ")
    end
  else
    printelms(io, vertices(rings[1]), tab="")
  end
end
