# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# default differentation backend
const FORWARDDIFF = DI.AutoForwardDiff()

"""
    derivative(geom, uvw, j; dbackend)

Calculate the derivative of the `geom`etry's parametric function
at parametric coordinates `uvw` and along `j`-th coordinate using
a differentiation `dbackend` from DifferentiationInterface.jl.
Default to forward mode automatic differentiation.
"""
function derivative(geom::Geometry, uvw, j; dbackend=FORWARDDIFF)
  # sanity check
  d = paramdim(geom)
  n = length(uvw)
  d == n || throw(ArgumentError("invalid number of parametric coordinates for geometry"))
  1 ≤ j ≤ n || throw(ArgumentError("attempting to compute derivative along invalid coordinate"))

  # strip units to help differentiation backends
  f(t) = ustrip.(to(geom(ntuple(i -> i == j ? t : uvw[i], d)...)))

  # compute derivative and re-add unit
  ∂ = DI.derivative(f, dbackend, uvw[j])
  u = unit(lentype(geom))

  Vec((∂ .* u)...)
end

"""
    jacobian(geom, uvw; dbackend)

Calculate the Jacobian of the `geom`etry's parametric function
at parametric coordinates `uvw` using a differentiation `dbackend`
from DifferentiationInterface.jl. Returns a tuple of vectors, each
corresponding to the derivative along a parametric coordinate.
Default to forward mode automatic differentiation.
"""
jacobian(geom::Geometry, uvw; dbackend=FORWARDDIFF) = ntuple(j -> derivative(geom, uvw, j; dbackend), paramdim(geom))

"""
    differential(geom, uvw; dbackend)

Calculate the differential element (length, area, volume, etc.)
of the `geom`etry at parametric coordinates `uvw` using a
differentiation `dbackend` from DifferentiationInterface.jl.
Default to forward mode automatic differentiation.
"""
function differential(geom::Geometry, uvw; dbackend=FORWARDDIFF)
  J = jacobian(geom, uvw; dbackend)
  if length(J) == 1
    norm(J[1])
  elseif length(J) == 2
    norm(J[1] × J[2])
  elseif length(J) == 3
    abs((J[1] × J[2]) ⋅ J[3])
  end
end
