# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KBallSearch(domain, k, ball)

A method that searches `k` nearest neighbors and then filters
these neighbors using a metric `ball`.

See [`MetricBall`](@ref) for additional details.
"""
struct KBallSearch{C<:CRS,T,R} <: BoundedNeighborSearchMethod
  tree::T
  radius::R
  k::Int
end

function KBallSearch(domain::D, k::Int, ball::B) where {D<:Domain,B<:MetricBall}
  C = crs(domain)
  m = metric(ball)
  r = ustrip(unit(lentype(C)), radius(ball))
  X = [svec(centroid(domain, i)) for i in 1:nelements(domain)]
  t = m isa MinkowskiMetric ? KDTree(X, m) : BallTree(X, m)
  KBallSearch{C,typeof(t),typeof(r)}(t, r, k)
end

KBallSearch(geoms, k, ball) = KBallSearch(GeometrySet(geoms), k, ball)

maxneighbors(method::KBallSearch) = method.k

function searchdists!(neighbors, distances, pₒ::Point, method::KBallSearch{C}; mask=nothing) where {C<:CRS}
  t = method.tree
  r = method.radius
  k = method.k

  # retrieve unit of length for distances
  u = unit(lentype(C))

  # adjust CRS of query point
  x = svec(convert(C, coords(pₒ)))

  inds, dists = knn(t, x, k, true)

  # keep neighbors inside ball
  keep = dists .≤ r

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
