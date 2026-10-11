# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometricDistance

A distance between points in geometries.
"""
abstract type GeometricDistance end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("distances/euclidean.jl")
include("distances/geodesic.jl")

# -------------
# DISTANCES.JL
# -------------

"""
    evaluate(distance::PreMetric, point₁, point₂)

Evaluate pre-metric `distance` between coordinates of `point₁` and `point₂`.
"""
function evaluate(d::PreMetric, p₁::Point, p₂::Point)
  u₁ = unit(Meshes.lentype(p₁))
  u₂ = unit(Meshes.lentype(p₂))
  u = Unitful.promote_unit(u₁, u₂)
  v₁ = ustrip.(u, to(p₁))
  v₂ = ustrip.(u, to(p₂))
  evaluate(d, v₁, v₂) * u
end
