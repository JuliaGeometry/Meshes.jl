# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NeighborhoodSearch(domain, neighborhood)

A method for searching neighbors in `domain` inside `neighborhood`.
"""
struct NeighborhoodSearch{D,N,T} <: NeighborSearchMethod
  # input fields
  domain::D
  neigh::N

  # state fields
  tree::T
end

function NeighborhoodSearch(domain::D, neigh::N) where {D,N}
  tree = if neigh isa MetricBall
    m = metric(neigh)
    X = coordinates(domain, 1:nelements(domain))
    m isa MinkowskiMetric ? KDTree(X, m) : BallTree(X, m)
  else
    nothing
  end
  NeighborhoodSearch{D,N,typeof(tree)}(domain, neigh, tree)
end

searchinds(pₒ, method::NeighborhoodSearch{D,N,T}) where {D,N<:NormBall,T} =
  inrange(method.tree, coordinates(pₒ), radius(method.neigh))

searchinds(pₒ, method::NeighborhoodSearch{D,N,T}) where {D,N<:Ellipsoid,T} =
  inrange(method.tree, coordinates(pₒ), one(coordtype(pₒ)))

# search method for ball neighborhood
function search(pₒ::Point, method::NeighborhoodSearch; mask=nothing)
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
