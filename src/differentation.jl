# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# default differentation method
const FINITEDIFF = DI.AutoFiniteDifferences(fdm=central_fdm(5, 1))

"""
    derivative(geom, uvw, j; ∂backend)

Calculate the derivative of the `geom`etry's parametric function
at parametric coordinates `uvw` and along `j`-th coordinate using
a differentiation `∂backend` from DifferentationInterface.jl.
By default, the `∂backend` is set to finite differences.
"""
function derivative(geom::Geometry, uvw, j; ∂backend=FINITEDIFF)
  # sanity check
  d = paramdim(geom)
  n = length(uvw)
  d == n || throw(ArgumentError("invalid number of parametric coordinates for geometry"))
  1 ≤ j ≤ n || throw(ArgumentError("attempting to compute derivative along invalid coordinate"))

  # strip units to help differentiation backends
  f(t) = ustrip.(to(geom(ntuple(i -> i == j ? t : uvw[i], d)...)))

  # compute derivative and re-add unit
  ∂ = DI.derivative(f, ∂backend, uvw[j])
  u = unit(lentype(geom))

  Vec((∂ .* u)...)
end

"""
    jacobian(geom, uvw; ∂backend)

Calculate the Jacobian of the `geom`etry's parametric function
at parametric coordinates `uvw` using a differentiation `∂backend`
from DifferentationInterface.jl. Returns a tuple of vectors, each
corresponding to the derivative along a parametric coordinate.
By default, `∂backend` is set to finite differences.
"""
jacobian(geom::Geometry, uvw; ∂backend=FINITEDIFF) = ntuple(j -> derivative(geom, uvw, j; ∂backend), paramdim(geom))

"""
    differential(geom, uvw; ∂backend)

Calculate the differential element (length, area, volume, etc.)
of the `geom`etry at parametric coordinates `uvw` using a
differentiation `∂backend` from DifferentiationInterface.jl.
By default, the `∂backend` is set to finite differences.
"""
function differential(geom::Geometry, uvw; ∂backend=FINITEDIFF)
  J = jacobian(geom, uvw; ∂backend)
  if length(J) == 1
    norm(J[1])
  elseif length(J) == 2
    norm(J[1] × J[2])
  elseif length(J) == 3
    abs((J[1] × J[2]) ⋅ J[3])
  end
end
