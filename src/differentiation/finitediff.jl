# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FiniteDifference(ϵ=1e-6)

Finite-difference method with unitless step size `ϵ`.

The step size is applied to parametric coordinates, which
are unitless values in the interval [0, 1].
"""
struct FiniteDifference{T} <: DifferentiationMethod
  ϵ::T
end

FiniteDifference() = FiniteDifference(1e-6)

function derivative(geom::Geometry, uvw, j, method::FiniteDifference)
  # sanity check
  pdim = paramdim(geom)
  if pdim != length(uvw)
    throw(ArgumentError("invalid number of parametric coordinates for geometry"))
  end

  # unitless step size
  ϵ = method.ϵ

  # central difference
  pre = ntuple(i -> i == j ? uvw[i] - ϵ : uvw[i], pdim)
  pos = ntuple(i -> i == j ? uvw[i] + ϵ : uvw[i], pdim)
  (geom(pos...) - geom(pre...)) / 2ϵ
end
