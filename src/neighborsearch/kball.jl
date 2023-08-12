# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KBallSearch(domain, k, ball)

A method that searches `k` nearest neighbors and then filters
these neighbors using a norm `ball`.
"""
struct KBallSearch{D,B<:MetricBall,T} <: BoundedNeighborSearchMethod
  # input fields
  domain::D
  k::Int
  ball::B

  # state fields
  tree::T
end

function KBallSearch(domain::D, k::Int, ball::B) where {D,B}
  m = metric(ball)
  xs = [coordinates(centroid(domain, i)) for i in 1:nelements(domain)]
  tree = m isa MinkowskiMetric ? KDTree(xs, m) : BallTree(xs, m)
  KBallSearch{D,B,typeof(tree)}(domain, k, ball, tree)
end

maxneighbors(method::KBallSearch) = method.k

function searchdists!(neighbors, distances, pₒ::Point, method::KBallSearch; skip=i -> false)
  tree = method.tree
  dmax = radius(method.ball)
  k = method.k

  inds, dists = knn(tree, coordinates(pₒ), k, true, skip)

  # keep neighbors inside ball
  keep = dists .≤ dmax

  nneigh = 0
  @inbounds for i in 1:length(keep)
    if keep[i]
      nneigh += 1
      neighbors[nneigh] = inds[i]
      distances[nneigh] = dists[i]
    end
  end

  nneigh
end
