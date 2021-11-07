# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallSearch(domain, ball)

A method for searching neighbors in `domain` inside `ball`.
"""
struct BallSearch{D,B<:MetricBall,T} <: NeighborSearchMethod
  # input fields
  domain::D
  ball::B

  # state fields
  tree::T
end

function BallSearch(domain::D, ball::B) where {D,B<:MetricBall}
  m  = metric(ball)
  xs = [coordinates(centroid(domain, i)) for i in 1:nelements(domain)]
  tree = m isa MinkowskiMetric ? KDTree(xs, m) : BallTree(xs, m)
  BallSearch{D,B,typeof(tree)}(domain, ball, tree)
end

function search(pₒ::Point, method::BallSearch; mask=nothing)
  inds = searchinds(pₒ, method)
  if mask ≠ nothing
    neighbors = Vector{Int}()
    @inbounds for ind in inds
      if mask[ind]
        push!(neighbors, ind)
      end
    end
    neighbors
  else
    inds
  end
end

searchinds(pₒ, method::BallSearch{D,B,T}) where {D,B<:IsotropicBall,T} =
  inrange(method.tree, coordinates(pₒ), first(radii(method.ball)))

searchinds(pₒ, method::BallSearch{D,B,T}) where {D,B<:AnisotropicBall,T} =
  inrange(method.tree, coordinates(pₒ), one(coordtype(pₒ)))