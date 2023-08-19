# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Bridge(δ=0)

Transform polygon with holes into a single outer ring
via bridges of given width `δ` as described in Held 1998.

## References

* Held. 1998. [FIST: Fast Industrial-Strength Triangulation of Polygons]
  (https://link.springer.com/article/10.1007/s00453-001-0028-4)
"""
struct Bridge{T} <: GeometricTransform
  δ::T
end

Bridge() = Bridge(nothing)

function apply(transform::Bridge, poly::Polygon{Dim,T}) where {Dim,T}
  δ = isnothing(transform.δ) ? zero(T) : transform.δ
  if hasholes(poly)
    bridge(rings(poly), δ)
  else
    first(rings(poly)), []
  end
end

function bridge(r::AbstractVector{<:Ring{2,T}}, δ) where {T}
  # sort rings lexicographically
  rings, vinds = repair9(r)
  verts = vertices.(rings)

  # initialize outer boundary
  outer = verts[1]
  oinds = vinds[1]

  # merge holes into outer boundary
  for i in 2:length(verts)
    inner = verts[i]
    iinds = vinds[i]

    # find closest pair of vertices (A, B)
    # connecting outer and inner rings
    omax = 0
    imax = 0
    dmin = typemax(T)
    for jₒ in 1:length(outer), jᵢ in 1:length(inner)
      d = sum(abs, outer[jₒ] - inner[jᵢ])
      if d < dmin
        omax = jₒ
        imax = jᵢ
        dmin = d
      end
    end
    A = outer[omax]
    B = inner[imax]

    # direction and normal to segment A--B
    v = B - A
    u = Vec(-v[2], v[1])
    n = u / norm(u)

    # the point A is split into A′ and A′′ and
    # the point B is split into B′ and B′′ based
    # on a given bridge width δ
    A′ = A + (δ / 2) * n
    A′′ = A - (δ / 2) * n
    B′ = B + (δ / 2) * n
    B′′ = B - (δ / 2) * n

    # insert hole at closest vertex
    outer = [
      outer[begin:(omax - 1)]
      [A′, B′]
      circshift(inner, -imax + 1)[2:end]
      [B′′, A′′]
      outer[(omax + 1):end]
    ]
    oinds = [
      oinds[begin:omax]
      circshift(iinds, -imax + 1)
      [iinds[imax]]
      oinds[omax:end]
    ]
  end

  # find duplicate vertices
  dups = Tuple{Int,Int}[]
  seen = Dict{Int,Int}()
  for (i, ind) in enumerate(oinds)
    if haskey(seen, ind)
      push!(dups, (seen[ind], i))
    else
      seen[ind] = i
    end
  end

  Ring(outer), dups
end
