# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DirectionPartition(direction; tol=1e-6)

A method for partitioning spatial objects along a given `direction`
with bandwidth tolerance `tol`.
"""
struct DirectionPartition{Dim,T} <: SPredicatePartitionMethod
  direction::Vec{Dim,T}
  tol::Float64

  function DirectionPartition{Dim,T}(direction, tol) where {Dim,T}
    return new(normalize(direction), tol)
  end
end

function DirectionPartition(direction::Vec{Dim,T}; tol=1e-6) where {Dim,T}
  return DirectionPartition{Dim,T}(direction, tol)
end

function DirectionPartition(direction::NTuple{Dim,T}; tol=1e-6) where {Dim,T}
  return DirectionPartition(Vec(direction); tol=tol)
end

function (p::DirectionPartition)(x, y)
  δ = x - y
  d = p.direction
  return norm(δ - (δ ⋅ d) * d) < p.tol
end
