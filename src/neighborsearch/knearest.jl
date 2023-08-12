# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KNearestSearch(domain, k; metric=Euclidean())

A method for searching `k` nearest neighbors in `domain`
according to `metric`.
"""
struct KNearestSearch{D,T} <: BoundedNeighborSearchMethod
  # input fields
  domain::D
  k::Int

  # state fields
  tree::T
end

function KNearestSearch(domain::D, k::Int; metric=Euclidean()) where {D}
  xs = [coordinates(centroid(domain, i)) for i in 1:nelements(domain)]
  tree = metric isa MinkowskiMetric ? KDTree(xs, metric) : BallTree(xs, metric)
  KNearestSearch{D,typeof(tree)}(domain, k, tree)
end

maxneighbors(method::KNearestSearch) = method.k

function searchdists!(neighbors, distances, pₒ::Point, method::KNearestSearch; skip=i->false)
  tree = method.tree
  k = method.k

  inds, dists = knn(tree, coordinates(pₒ), k, true, skip)

  nneigh = length(inds)

  @inbounds for i in 1:nneigh
    neighbors[i] = inds[i]
    distances[i] = dists[i]
  end

  nneigh
end
