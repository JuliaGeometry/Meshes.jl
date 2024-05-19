# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KBallSearch(domain, k, ball)

A method that searches `k` nearest neighbors and then filters
these neighbors using a norm `ball`.
"""
struct KBallSearch{D<:Domain,B<:MetricBall,T} <: BoundedNeighborSearchMethod
  # input fields
  domain::D
  k::Int
  ball::B

  # state fields
  tree::T
end

function KBallSearch(domain::D, k::Int, ball::B) where {D<:Domain,B<:MetricBall}
  m = metric(ball)
  xs = [ustrip.(coordinates(centroid(domain, i))) for i in 1:nelements(domain)]
  tree = m isa MinkowskiMetric ? KDTree(xs, m) : BallTree(xs, m)
  KBallSearch{D,B,typeof(tree)}(domain, k, ball, tree)
end

KBallSearch(geoms, k, ball) = KBallSearch(GeometrySet(geoms), k, ball)

maxneighbors(method::KBallSearch) = method.k

function searchdists!(neighbors, distances, pₒ::Point, method::KBallSearch; mask=nothing)
  u = unit(lentype(pₒ))
  tree = method.tree
  dmax = radius(method.ball)
  k = method.k

  inds, dists = knn(tree, ustrip.(coordinates(pₒ)), k, true)

  # keep neighbors inside ball
  keep = dists .≤ dmax

  # possibly mask some of the neighbors
  isnothing(mask) || (keep .*= mask[inds])

  nneigh = 0
  @inbounds for i in 1:k
    if keep[i]
      nneigh += 1
      neighbors[nneigh] = inds[i]
      distances[nneigh] = dists[i] * u
    end
  end

  nneigh
end
