# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    integral(fun, geom; n=3)

Calculate the integral over the `geom`etry of the `fun`ction that maps
[`Point`](@ref)s to values in a linear space.

    integral(fun, dom; n=3)

Alternatively, calculate the integral over the `dom`ain (e.g., mesh) by
summing the integrals for each constituent geometry.

Polynomials of degree up to `2n-1` are integrated exactly.

See also [`localintegral`](@ref).
"""
integral(fun, geom::Geometry; n=3) = localintegral(fun ∘ geom, geom; n)

integral(fun, dom::Domain; n=3) = sum(integral(fun, geom; n) for geom in dom)

"""
    localintegral(fun, geom; n=3)

Calculate the integral over the `geom`etry of the `fun`ction that maps
parametric coordinates `uvw` to values in a linear space.

Polynomials of degree up to `2n-1` are integrated exactly.

See also [`integral`](@ref).
"""
localintegral(fun, geom::Geometry; n=3) = _uvwintegral(fun, geom, n)

# --------------
# SPECIAL CASES
# --------------

# the integral of a (bounded) function over a single point is always zero
localintegral(fun, point::Point; n=3) = zero(lentype(point))

# cylinder surface is the union of the lateral surface and the top and bottom disks
localintegral(fun, cylsurf::CylinderSurface; n=3) =
  _uvwintegral(fun, cylsurf, n) + _uvwintegral(fun, top(cylsurf), n) + _uvwintegral(fun, bottom(cylsurf), n)

# cone surface is the union of the lateral surface and the base disk
localintegral(fun, conesurf::ConeSurface; n=3) = _uvwintegral(fun, conesurf, n) + _uvwintegral(fun, base(conesurf), n)

# frustum surface is the union of the lateral surface and the top and bottom disks
localintegral(fun, frustumsurf::FrustumSurface; n=3) =
  _uvwintegral(fun, frustumsurf, n) + _uvwintegral(fun, top(frustumsurf), n) + _uvwintegral(fun, bottom(frustumsurf), n)

# chain is the union of its segments
localintegral(fun, chain::Chain; n=3) = sum(_uvwintegral(fun, seg, n) for seg in segments(chain))

# -----------------
# HELPER FUNCTIONS
# -----------------

function _uvwintegral(fun, geom, n)
  # parametric dimension and number type
  N = paramdim(geom)
  T = numtype(lentype(geom))

  # Gauss-Legendre quadrature points and weights
  ts, ws = gausslegendre(n)
  tgrid = Iterators.product(ntuple(_ -> T.(ts), N)...)
  wgrid = Iterators.product(ntuple(_ -> T.(ws), N)...)

  # compute integral with change of variable and differential element
  Σwᵢfᵢ = sum(zip(tgrid, wgrid)) do (t, w)
    # change of variable [-1, 1] → [0, 1]
    uvw = ntuple(i -> (t[i] + 1) / 2, N)
    prod(w) * fun(uvw...) * differential(geom, uvw)
  end

  # adjust for change of variable
  Σwᵢfᵢ / 2^N
end
