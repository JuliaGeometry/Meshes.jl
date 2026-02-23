# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FiniteDifference(eps=1e-6)

Finite-difference method with unitless step size `eps`.

The step size is applied to parametric coordinates, which
are unitless values in the interval [0, 1].
"""
struct FiniteDifference{T} <: DifferentiationMethod
  eps::T
end

FiniteDifference() = FiniteDifference(1e-6)

function deriv(geom::Geometry, uvw, j, method::FiniteDifference)
  eps = method.eps
  dim = paramdim(geom)
  pre = ntuple(i -> i == j ? uvw[i] - eps : uvw[i], dim)
  pos = ntuple(i -> i == j ? uvw[i] + eps : uvw[i], dim)
  (geom(pos...) - geom(pre...)) / 2eps
end
