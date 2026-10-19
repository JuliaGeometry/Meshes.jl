# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MahalanobisDistance(matrix)

A Mahalanobis distance defined by a positive semidefinite symmetric `matrix`.
Given two points `p` and `q`, it is defined as `√((p - q)' * matrix * (p - q))`.

    MahalanobisDistance(radii, rotation=I)

Alternatively, construct the `matrix` based on `radii` for rotated axes
that are defined by a `rotation` specification from the package
[`Rotations.jl`](https://github.com/JuliaGeometry/Rotations.jl).
"""
struct MahalanobisDistance{M<:Symmetric} <: GeometricDistance
  matrix::M
  function MahalanobisDistance(matrix::M) where {M<:Symmetric}
    assertion(eigmin(matrix) ≥ 0, "Mahalanobis matrix must be positive semidefinite")
    new{M}(matrix)
  end
end

function MahalanobisDistance(radii::Tuple, rotation=I)
  r = ustrip.(promote(radii...))
  w = SVector(inv.(r) .^ 2)
  W = Diagonal(w)
  R = rotation
  Q = R * W * transpose(R)
  MahalanobisDistance(Symmetric(Q))
end

function (d::MahalanobisDistance)(p₁::Point, p₂::Point)
  v = p₁ - p₂
  Q = d.matrix
  √(v ⋅ (Q * v))
end
