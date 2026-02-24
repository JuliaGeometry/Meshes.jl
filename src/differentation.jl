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
function derivative(geom::Geometry, uvw, j, method::DifferentiationMethod)
  d = paramdim(geom)
  n = length(uvw)
  d == n || throw(ArgumentError("invalid number of parametric coordinates for geometry"))
  1 ≤ j ≤ n || throw(ArgumentError("attempting to compute derivative along invalid coordinate"))
  deriv(geom, uvw, j, method)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("differentiation/finitediff.jl")
include("differentiation/autodiff.jl")

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
