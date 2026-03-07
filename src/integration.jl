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

# cylinder surface is the union of lateral surface and top/bottom disks
integral(fun, cylsurf::CylinderSurface; n=3) =
  localintegral(fun ∘ cylsurf, cylsurf; n) + integral(fun, top(cylsurf); n) + integral(fun, bottom(cylsurf); n)

# cone surface is the union of lateral surface and base disk
integral(fun, conesurf::ConeSurface; n=3) =
  localintegral(fun ∘ conesurf, conesurf; n) + integral(fun, base(conesurf); n)

# frustum surface is the union of lateral surface and top/bottom disks
integral(fun, frustumsurf::FrustumSurface; n=3) =
  localintegral(fun ∘ frustumsurf, frustumsurf; n) +
  integral(fun, top(frustumsurf); n) +
  integral(fun, bottom(frustumsurf); n)

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

# ray is parametrized over [0, ∞] interval
# first map: [-1, 1] → [0, 1] with t -> (t + 1) / 2
# second map: [0, 1] → [0, ∞] with t -> t / (1 - t^2)
localintegral(fun, ray::Ray; n=3) = _uvwintegral(fun, ray, n, trans=(t -> t / (1 - t^2)) ∘ (t -> (t + 1) / 2))

# line is parametrized over [-∞, ∞] interval
localintegral(fun, line::Line; n=3) = _uvwintegral(fun, line, n, trans=t -> t / (1 - t^2))

# plane is parametrized over [-∞, ∞] interval
localintegral(fun, plane::Plane; n=3) = _uvwintegral(fun, plane, n, trans=t -> t / (1 - t^2))

# chain is the union of its segments
localintegral(fun, chain::Chain; n=3) = sum(localintegral(fun, seg; n) for seg in segments(chain))

# polygon is the union of simpler n-gons
localintegral(fun, poly::Polygon; n=3) = sum(localintegral(fun, ngon; n) for ngon in discretize(poly))

# triangle is parametrized with barycentric coordinates
# TODO:

# specialize quadrangle for performance
localintegral(fun, quad::Quadrangle; n=3) = _uvwintegral(fun, quad, n)

# tetrahedron is parametrized with barycentric coordinates
# TODO:

# multi-geometry is the union of simpler geometries
localintegral(fun, multi::Multi; n=3) = sum(localintegral(fun, geom; n) for geom in parent(multi))

# -----------------
# HELPER FUNCTIONS
# -----------------

# we set t -> (t + 1) / 2 by default to map [-1, 1] → [0, 1]
# i.e., quadrature nodes to parametric coordinates in [0, 1]
function _uvwintegral(fun, geom, n; trans=t -> (t + 1) / 2)
  # parametric dimension and number type
  N = paramdim(geom)
  T = numtype(lentype(geom))

  # Gauss-Legendre quadrature points and weights
  ts, ws = gausslegendre(n)
  tgrid = Iterators.product(ntuple(_ -> T.(ts), N)...)
  wgrid = Iterators.product(ntuple(_ -> T.(ws), N)...)

  # compute integral with change of variable and differential element
  Σwᵢfᵢ = sum(zip(tgrid, wgrid)) do (t, w)
    uvw = ntuple(i -> trans(t[i]), N)
    prod(w) * fun(uvw...) * differential(geom, uvw)
  end

  # adjust for change of variable
  Σwᵢfᵢ / 2^N
end
