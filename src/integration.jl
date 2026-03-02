# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    integral(fun, geom; n=3)

Calculate the integral over the `geom`etry of the `fun`ction that maps
[`Point`](@ref)s to values in a linear space.

Polynomials of degree up to `2n-1` are integrated exactly.

See also [`localintegral`](@ref).
"""
integral(fun, geom; n=3) = localintegral(fun ∘ geom, geom; n)

"""
    localintegral(fun, geom; n=3)

Calculate the integral over the `geom`etry of the `fun`ction that maps
parametric coordinates `uvw` to values in a linear space.

Polynomials of degree up to `2n-1` are integrated exactly.

See also [`integral`](@ref).
"""
function localintegral(fun, geom; n=3)
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
