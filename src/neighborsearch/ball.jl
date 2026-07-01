# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallSearch(domain, ball)

A method for searching neighbors in `domain` inside metric `ball`.

See [`MetricBall`](@ref) for additional details.
"""
struct BallSearch{C<:CRS,T,R} <: NeighborSearchMethod
  tree::T
  radius::R
end

function BallSearch(domain::D, ball::B) where {D<:Domain,B<:MetricBall}
  C = crs(domain)
  m = metric(ball)
  r = ustrip(unit(lentype(C)), radius(ball))
  X = [svec(centroid(domain, i)) for i in 1:nelements(domain)]
  t = m isa MinkowskiMetric ? KDTree(X, m) : BallTree(X, m)
  BallSearch{C,typeof(t),typeof(r)}(t, r)
end

BallSearch(geoms, ball) = BallSearch(GeometrySet(geoms), ball)

function search(pₒ::Point, method::BallSearch{C}; mask=nothing) where {C<:CRS}
  t = method.tree
  r = method.radius

  # adjust CRS of query point
  x = svec(convert(C, coords(pₒ)))

  inds = inrange(t, x, r)

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
