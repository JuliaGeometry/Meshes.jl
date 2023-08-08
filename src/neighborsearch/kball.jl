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

search!(neighbors, pₒ::Point, method::KBallSearch; mask=nothing) = first(searchwithdist!(neighbors, pₒ, method; mask))

function searchwithdist!(neighbors, pₒ::Point, method::KBallSearch; mask=nothing)
  k = method.k
  r = radius(method.ball)

  inds, dists = knn(method.tree, coordinates(pₒ), k, true)

  # keep neighbors inside ball
  keep = dists .≤ r

  # possibly mask some of the neighbors
  isnothing(mask) || (keep .*= mask[inds])

  nneigh = 0
  neighdists = empty(dists)
  @inbounds for i in 1:k
    if keep[i]
      nneigh += 1
      neighbors[nneigh] = inds[i]
      push!(neighdists, dists[i])
    end
  end

  nneigh, neighdists
end
