# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DirectionPartition(direction; tol=1e-6)

A method for partitioning spatial objects along a given `direction`
with bandwidth tolerance `tol`.
"""
struct DirectionPartition{V<:Vec} <: SPredicatePartitionMethod
  direction::V
  tol::Float64

  function DirectionPartition{V}(direction, tol) where {V<:Vec}
    new(normalize(direction), tol)
  end
end

DirectionPartition(direction::V; tol=1e-6) where {V<:Vec} = DirectionPartition{V}(direction, tol)

DirectionPartition(direction::Tuple; tol=1e-6) = DirectionPartition(Vec(direction), tol=tol)

function (p::DirectionPartition)(x, y)
  δ = x - y
  d = p.direction
  norm(δ - (δ ⋅ d) * d) < p.tol
end
