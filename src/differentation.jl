# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DifferentiationMethod

A method for calculating derivatives.
"""
abstract type DifferentiationMethod end

"""
    derivative(geom, uvw, j, method)

Calculate the derivative of the `geom`etry's parametric function
at parametric coordinates `uvw` and along `j`-th coordinate using
a differentiation `method`.
"""
function derivative end

# ----------------
# IMPLEMENTATIONS
# ----------------

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
  (geom(pos...) - geom(pre...)) / 2eps
end

# ----------
# FALLBACKS
# ----------

"""
    jacobian(geom, uvw, method=FiniteDifference())

Calculate the Jacobian of the `geom`etry's parametric function
at parametric coordinates `uvw` using a differentiation `method`.
"""
jacobian(geom::Geometry, uvw, method=FiniteDifference()) = ntuple(j -> derivative(geom, uvw, j, method), paramdim(geom))

"""
    differential(geom, uvw, method=FiniteDifference())

Calculate the differential element (length, area, volume, etc.)
of the `geom`etry at parametric coordinates `uvw` using a
differentation `method`.
"""
function differential(geom::Geometry, uvw, method=FiniteDifference())
  J = jacobian(geom, uvw, method)
  if length(J) == 1
    norm(J[1])
  elseif length(J) == 2
    norm(J[1] × J[2])
  elseif length(J) == 3
    abs((J[1] × J[2]) ⋅ J[3])
  end
end
