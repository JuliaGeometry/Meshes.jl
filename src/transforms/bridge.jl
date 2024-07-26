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
struct Bridge{ℒ<:Len} <: GeometricTransform
  δ::ℒ
  Bridge(δ::ℒ) where {ℒ<:Len} = new{float(ℒ)}(δ)
end

Bridge(δ) = Bridge(addunit(δ, u"m"))

Bridge() = Bridge(0.0u"m")

parameters(t::Bridge) = (; δ=t.δ)

function apply(transform::Bridge, poly::PolyArea)
  ℒ = lentype(poly)

  # sort rings lexicographically
  rpoly, rinds = apply(Repair{9}(), poly)

  # retrieve bridge width
  δ = convert(ℒ, transform.δ)

  ring, dups = if hasholes(rpoly)
    bridge(rings(rpoly), rinds, δ)
  else
    first(rings(rpoly)), []
  end

  PolyArea(ring), dups # TODO: check orientation
end

apply(::Bridge, poly::Ngon) = poly, []

function bridge(rings, rinds, δ)
  # extract vertices and indices
  verts = vertices.(rings)
  vinds = rinds

  # retrieve coordinate type
  ℒ = lentype(first(rings))

  # retrieve original CRS
  C = crs(first(rings))

  # initialize outer boundary
  outer = flat.(verts[1])
  oinds = vinds[1]

  # merge holes into outer boundary
  for i in 2:length(verts)
    inner = flat.(verts[i])
    iinds = vinds[i]

    # find closest pair of vertices (A, B)
    # connecting outer and inner rings
    omax = 0
    imax = 0
    dmin = typemax(ℒ)
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
    n = norm(u)

    # the point A is split into A′ and A′′ and
    # the point B is split into B′ and B′′ based
    # on a given bridge width δ
    A′ = A + (δ / 2n) * u
    A′′ = A - (δ / 2n) * u
    B′ = B + (δ / 2n) * u
    B′′ = B - (δ / 2n) * u

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

  points = map(outer) do p
    c = CoordRefSystems.raw(coords(p))
    Point(CoordRefSystems.reconstruct(C, c))
  end

  Ring(points), dups
end
