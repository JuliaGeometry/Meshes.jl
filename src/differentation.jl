# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DifferentiationMethod

A method for calculating derivatives.
"""
abstract type DifferentiationMethod end

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

"""
    jacobian(geom, uvw, method=FiniteDifference())

Calculate the Jacobian of a geometry's parametric function at parametric coordinates `uvw`
using a particular differentiation `method`.
"""
jacobian(geom::Geometry, uvw, method=FiniteDifference()) = jacobianimpl(geom, uvw, method)

function jacobianimpl(geom::Geometry, uvw, method::FiniteDifference)
  # sanity check
  pdim = paramdim(geom)
  if pdim != length(uvw)
    throw(ArgumentError("invalid number of parametric coordinates for geometry"))
  end

  # unitless step size
  ϵ = method.ϵ

  # partial derivatives along dimensions j=1,2,...
  ntuple(pdim) do j
    pre = ntuple(i -> i == j ? uvw[i] - ϵ : uvw[i], pdim)
    pos = ntuple(i -> i == j ? uvw[i] + ϵ : uvw[i], pdim)
    if uvw[j] < 0.01 # right
      (geom(pos...) - geom(uvw...)) / ϵ
    elseif 0.99 < uvw[j] # left
      (geom(uvw...) - geom(pre...)) / ϵ
    else # central
      (geom(pos...) - geom(pre...)) / 2ϵ
    end
  end
end

"""
    differential(geom, uvw, method=FiniteDifference())

Calculate the differential element (length, area, volume, etc.)
of the geometry `geom` at parametric coordinates `uvw` using a
particular differentation `method`.
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
