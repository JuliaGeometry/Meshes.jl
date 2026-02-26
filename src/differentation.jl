# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# default differentation method
const FINITEDIFF = DI.AutoFiniteDifferences(; fdm=central_fdm(5, 1))

"""
    derivative(geom, uvw, j[, method])

Calculate the derivative of the `geom`etry's parametric function
at parametric coordinates `uvw` and along `j`-th coordinate using
an automatic differentiation `method` from DifferentationInterface.jl.
"""
function derivative(geom::Geometry, uvw, j, method::DI.AbstractADType=FINITEDIFF)
  # sanity check
  d = paramdim(geom)
  n = length(uvw)
  d == n || throw(ArgumentError("invalid number of parametric coordinates for geometry"))
  1 ≤ j ≤ n || throw(ArgumentError("attempting to compute derivative along invalid coordinate"))

  # strip units to help autodiff methods
  f(t) = ustrip.(to(geom(ntuple(i -> i == j ? t : uvw[i], d)...)))

  # compute derivative and re-add unit
  ∂ = DI.derivative(f, method, uvw[j])
  u = unit(lentype(geom))

  Vec((∂ .* u)...)
end

# ----------
# FALLBACKS
# ----------

"""
    jacobian(geom, uvw[, method])

Calculate the Jacobian of the `geom`etry's parametric function
at parametric coordinates `uvw` using an automatic differentiation
`method` from DifferentationInterface.jl.
"""
jacobian(geom::Geometry, uvw, method=FINITEDIFF) = ntuple(j -> derivative(geom, uvw, j, method), paramdim(geom))

"""
    differential(geom, uvw[, method])

Calculate the differential element (length, area, volume, etc.)
of the `geom`etry at parametric coordinates `uvw` using an
automatic differentation `method` from DifferentiationInterface.jl.
"""
function differential(geom::Geometry, uvw, method=FINITEDIFF)
  J = jacobian(geom, uvw, method)
  if length(J) == 1
    norm(J[1])
  elseif length(J) == 2
    norm(J[1] × J[2])
  elseif length(J) == 3
    abs((J[1] × J[2]) ⋅ J[3])
  end
end
