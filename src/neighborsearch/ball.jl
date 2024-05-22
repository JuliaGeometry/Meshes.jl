# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallSearch(domain, ball)

A method for searching neighbors in `domain` inside `ball`.
"""
struct BallSearch{D<:Domain,B<:MetricBall,T} <: NeighborSearchMethod
  # input fields
  domain::D
  ball::B

  # state fields
  tree::T
end

function BallSearch(domain::D, ball::B) where {D<:Domain,B<:MetricBall}
  m = metric(ball)
  xs = [ustrip.(to(centroid(domain, i))) for i in 1:nelements(domain)]
  tree = m isa MinkowskiMetric ? KDTree(xs, m) : BallTree(xs, m)
  BallSearch{D,B,typeof(tree)}(domain, ball, tree)
end

BallSearch(geoms, ball) = BallSearch(GeometrySet(geoms), ball)

function search(pₒ::Point, method::BallSearch; mask=nothing)
  tree = method.tree
  dmax = radius(method.ball)

  inds = inrange(tree, ustrip.(to(pₒ)), ustrip(dmax))

  if isnothing(mask)
    inds
  else
    neighbors = Vector{Int}()
    @inbounds for ind in inds
      if mask[ind]
        push!(neighbors, ind)
      end
    end
    neighbors
  end
end
