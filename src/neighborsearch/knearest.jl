# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KNearestSearch(domain, k; distance=EuclideanDistance())

A method for searching `k` nearest neighbors in `domain`
according to `distance`.
"""
struct KNearestSearch{C<:CRS,T} <: BoundedNeighborSearchMethod
  tree::T
  k::Int
end

function KNearestSearch(domain::D, k::Int; distance=EuclideanDistance()) where {D<:Domain}
  C = crs(domain)
  t = _searchtree(domain, distance)
  KNearestSearch{C,typeof(t)}(t, k)
end

KNearestSearch(geoms, k; distance=EuclideanDistance()) = KNearestSearch(GeometrySet(geoms), k; distance)

maxneighbors(method::KNearestSearch) = method.k

function searchdists!(neighbors, distances, pₒ::Point, method::KNearestSearch{C}; mask=nothing) where {C<:CRS}
  tree = method.tree
  k = method.k

  # retrieve unit of length for distances
  u = unit(lentype(C))

  # adjust CRS of query point
  x = svec(convert(C, coords(pₒ)))

  inds, dists = knn(tree, x, k, true)

  if isnothing(mask)
    nneigh = k
    @inbounds for i in 1:k
      neighbors[i] = inds[i]
      distances[i] = dists[i] * u
    end
  else
    nneigh = 0
    @inbounds for i in 1:k
      if mask[inds[i]]
        nneigh += 1
        neighbors[nneigh] = inds[i]
        distances[nneigh] = dists[i] * u
      end
    end
  end

  nneigh
end
