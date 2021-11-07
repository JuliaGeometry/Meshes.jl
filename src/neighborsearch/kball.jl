# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KBallSearch(domain, k, ball)

A method that searches `k` nearest neighbors and then filters
these neighbors using a norm `ball`.
"""
struct KBallSearch{D,B,T} <: BoundedNeighborSearchMethod
  # input fields
  domain::D
  k::Int
  ball::B

  # state fields
  tree::T
end

function KBallSearch(domain::D, k::Int, ball::B) where {D,B}
  m  = metric(ball)
  xs = [coordinates(centroid(domain, i)) for i in 1:nelements(domain)]
  tree = m isa MinkowskiMetric ? KDTree(xs, m) : BallTree(xs, m)
  KBallSearch{D,B,typeof(tree)}(domain, k, ball, tree)
end

maxneighbors(method::KBallSearch) = method.k

function search!(neighbors, pₒ::Point, method::KBallSearch; mask=nothing)
  k = method.k
  r = method.ball isa IsotropicBall ? first(radii(method.ball)) : one(coordtype(pₒ))

  inds, dists = knn(method.tree, coordinates(pₒ), k, true)

  # keep neighbors inside ball
  keep = dists .≤ r

  # possibly mask some of the neighbors
  isnothing(mask) || (keep .*= mask[inds])

  nneigh = 0
  @inbounds for i in 1:k
    if keep[i]
      nneigh += 1
      neighbors[nneigh] = inds[i]
    end
  end

  nneigh
end
