# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# default integration method
const GAUSSLEGENDRE = II.Backend.Quadrature(gausslegendre(20))

"""
    integral(fun, geom[, method]) 

Calculate the integral over the `geom`etry of the `fun`ction that maps
[`Point`](@ref)s to values in a linear space.

    integral(fun, dom[, method])

Alternatively, calculate the integral over the `dom`ain (e.g., mesh) by
summing the integrals for each constituent geometry.

By default, use Gauss-Legendre quadrature rule with `n` nodes so that
polynomials of degree up to `2n-1` are integrated exactly.

See also [`localintegral`](@ref).
"""
integral(fun, geom::Geometry, method=GAUSSLEGENDRE) = localintegral(fun ∘ geom, geom, method)

# cylinder surface is the union of the curved surface and the top and bottom disks
integral(fun, cylsurf::CylinderSurface, method=GAUSSLEGENDRE) =
  localintegral(fun ∘ cylsurf, cylsurf, method) +
  integral(fun, top(cylsurf), method) +
  integral(fun, bottom(cylsurf), method)

# cone surface is the union of the curved surface and the base disk
integral(fun, conesurf::ConeSurface, method=GAUSSLEGENDRE) =
  localintegral(fun ∘ conesurf, conesurf, method) + integral(fun, base(conesurf), method)

# frustum surface is the union of the curved surface and the top and bottom disks
integral(fun, frustumsurf::FrustumSurface, method=GAUSSLEGENDRE) =
  localintegral(fun ∘ frustumsurf, frustumsurf, method) +
  integral(fun, top(frustumsurf), method) +
  integral(fun, bottom(frustumsurf), method)

# multi-geometry is the union of its constituent geometries
integral(fun, multi::Multi, method=GAUSSLEGENDRE) = sum(integral(fun, geom, method) for geom in parent(multi))

# domain is the union of its constituent geometries
integral(fun, dom::Domain, method=GAUSSLEGENDRE) = sum(integral(fun, geom, method) for geom in dom)

"""
    localintegral(fun, geom[, method])

Calculate the integral over the `geom`etry of the `fun`ction that maps
parametric coordinates `uvw` to values in a linear space.

By default, use Gauss-Legendre quadrature rule with `n` nodes so that
polynomials of degree up to `2n-1` are integrated exactly.

See also [`integral`](@ref).
"""
function localintegral(fun, geom::Geometry, method=GAUSSLEGENDRE)
  # integrand is equal to function times differential element
  integrand(uvw...) = fun(uvw...) * differential(geom, uvw)

  # domain of integration for the given geometry
  domain = ∫domain(geom)

  # perform numerical integration
  II.integral(integrand, domain; backend=method)()
end

function ∫domain(geom::Geometry)
  N = paramdim(geom)
  T = numtype(lentype(geom))
  a = ntuple(_ -> zero(T), N)
  b = ntuple(_ -> one(T), N)
  II.Domain.Box(a, b)
end

function ∫domain(ray::Ray)
  T = numtype(lentype(ray))
  a = (zero(T),)
  b = (II.Infinity(one(T)),)
  II.Domain.Box(a, b)
end

function ∫domain(line::Line)
  T = numtype(lentype(line))
  a = (-II.Infinity(one(T)),)
  b = (II.Infinity(one(T)),)
  II.Domain.Box(a, b)
end

function ∫domain(plane::Plane)
  T = numtype(lentype(plane))
  a = (-II.Infinity(one(T)), -II.Infinity(one(T)))
  b = (II.Infinity(one(T)), II.Infinity(one(T)))
  II.Domain.Box(a, b)
end
