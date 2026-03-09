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

integral(fun, cylsurf::CylinderSurface; n=3) =
  localintegral(fun ∘ cylsurf, cylsurf; n) + integral(fun, top(cylsurf); n) + integral(fun, bottom(cylsurf); n)

integral(fun, conesurf::ConeSurface; n=3) =
  localintegral(fun ∘ conesurf, conesurf; n) + integral(fun, base(conesurf); n)

integral(fun, frustumsurf::FrustumSurface; n=3) =
  localintegral(fun ∘ frustumsurf, frustumsurf; n) +
  integral(fun, top(frustumsurf); n) +
  integral(fun, bottom(frustumsurf); n)

integral(fun, multi::Multi; n=3) = sum(integral(fun, geom; n) for geom in parent(multi))

integral(fun, dom::Domain; n=3) = sum(integral(fun, geom; n) for geom in dom)

"""
    localintegral(fun, geom; n=3)

Calculate the integral over the `geom`etry of the `fun`ction that maps
parametric coordinates `uvw` to values in a linear space.

Polynomials of degree up to `2n-1` are integrated exactly.

See also [`integral`](@ref).
"""
localintegral(fun, geom::Geometry; n=3) = _uvwintegral(fun, geom, n)

# ray is parametrized over [0, ∞] interval
localintegral(fun, ray::Ray; n=3) = _uvwintegral(fun, ray, n, trans=t -> @. t / (1 - t))

# line is parametrized over [-∞, ∞] interval
localintegral(fun, line::Line; n=3) = _uvwintegral(fun, line, n, trans=t -> @. log(t) - log(1 - t))

# plane is parametrized over [-∞, ∞] interval
localintegral(fun, plane::Plane; n=3) = _uvwintegral(fun, plane, n, trans=t -> @. log(t) - log(1 - t))

# triangle is parametrized with barycentric coordinates
localintegral(fun, tri::Triangle; n=3) = _uvwintegral(fun, tri, n, trans=t -> (t[1] * t[2], t[2] - t[1] * t[2]))

# specialize quadrangle for performance
localintegral(fun, quad::Quadrangle; n=3) = _uvwintegral(fun, quad, n)

# tetrahedron is parametrized with barycentric coordinates
localintegral(fun, tetra::Tetrahedron; n=3) = _uvwintegral(
  fun,
  tetra,
  n,
  trans=t -> (t[2] * t[3] - t[1] * t[2] * t[3], t[1] * t[2] * t[3], t[3] - t[2] * t[3])
)

# -----------------
# HELPER FUNCTIONS
# -----------------

function _uvwintegral(fun, geom, n; trans=identity)
  # parametric dimension and number type
  N = paramdim(geom)
  T = numtype(lentype(geom))

  # Gauss-Legendre quadrature points and weights
  ts, ws = gausslegendre(n)
  tgrid = Iterators.product(ntuple(_ -> T.(ts), N)...)
  wgrid = Iterators.product(ntuple(_ -> T.(ws), N)...)

  # map quadrature points in [-1, 1] to parametric coordinates in [0, 1],
  # and then map parametric coordinates in [0, 1] to uvw parametrization
  g = trans ∘ (t -> @. (t + 1) / 2)

  # compute integral with change of variable and differential element
  Σwᵢfᵢ = sum(zip(tgrid, wgrid)) do (t, w)
    uvw = g(t)
    prod(w) * fun(uvw...) * differential(geom, uvw)
  end

  # adjust for change of variable
  Σwᵢfᵢ / 2^N
end
