# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    integral(fun, geom; order=1)

Calculate the integral over the `geom`etry of the `fun`ction that maps
[`Point`](@ref)s to values in a linear space.

The quadrature `order` can be specified to increase the accuracy.

See also [`localintegral`](@ref).
"""
integral(fun, geom; order=1) = localintegral(fun ∘ geom, geom; order)

"""
    localintegral(fun, geom; order=1)

Calculate the integral over the `geom`etry of the `fun`ction that maps
parametric coordinates `uvw` to values in a linear space.

The quadrature `order` can be specified to increase the accuracy.

See also [`integral`](@ref).
"""
function localintegral(fun, geom; order=1)
  # parametric dimension
  dim = paramdim(geom)

  # Gauss-Legendre quadrature points and weights
  ts, ws = gausslegendre(order)
  tgrid = Iterators.product(ntuple(_ -> ts, dim)...)
  wgrid = Iterators.product(ntuple(_ -> ws, dim)...)

  # compute integral with change of variable and differential element
  Σwᵢfᵢ = sum(zip(tgrid, wgrid)) do (t, w)
    # change of variable [-1, 1] → [0, 1]
    uvw = ntuple(i -> (t[i] + 1) / 2, dim)
    prod(w) * fun(uvw...) * differential(geom, uvw)
  end

  # adjust for change of variable
  Σwᵢfᵢ / 2^dim
end
