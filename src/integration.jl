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

# cylinder surface is the union of the curved surface and the top and bottom disks
integral(fun, cylsurf::CylinderSurface; n=3) =
  localintegral(fun ∘ cylsurf, cylsurf; n) + integral(fun, top(cylsurf); n) + integral(fun, bottom(cylsurf); n)

# cone surface is the union of the curved surface and the base disk
integral(fun, conesurf::ConeSurface; n=3) =
  localintegral(fun ∘ conesurf, conesurf; n) + integral(fun, base(conesurf); n)

# frustum surface is the union of the curved surface and the top and bottom disks
integral(fun, frustumsurf::FrustumSurface; n=3) =
  localintegral(fun ∘ frustumsurf, frustumsurf; n) +
  integral(fun, top(frustumsurf); n) +
  integral(fun, bottom(frustumsurf); n)

# multi-geometry is the union of its constituent geometries
integral(fun, multi::Multi; n=3) = sum(integral(fun, geom; n) for geom in parent(multi))

# domain is the union of its constituent geometries
integral(fun, dom::Domain; n=3) = sum(integral(fun, geom; n) for geom in dom)

"""
    localintegral(fun, geom; n=3)

Calculate the integral over the `geom`etry of the `fun`ction that maps
parametric coordinates `uvw` to values in a linear space.

Polynomials of degree up to `2n-1` are integrated exactly.

See also [`integral`](@ref).
"""
function localintegral(fun, geom::Geometry; n=3)
  # domain of integration
  domain = ∫dom(geom)

  # Gauss-Legendre quadrature
  backend = II.Backend.Quadrature(gausslegendre(n))

  # integral of function times differential element
  I = II.integral(uvw -> fun(uvw...) * differential(geom, uvw), domain; backend)

  # perform numerical integration
  I()
end

function ∫dom(geom::Geometry)
  N = paramdim(geom)
  T = numtype(lentype(geom))
  a = ntuple(_ -> zero(T), N)
  b = ntuple(_ -> one(T), N)
  II.Domain.Box(a, b)
end

function ∫dom(ray::Ray)
  T = numtype(lentype(ray))
  a = (zero(T),)
  b = (II.Infinity(one(T)),)
  II.Domain.Box(a, b)
end

function ∫dom(line::Line)
  T = numtype(lentype(line))
  a = (-II.Infinity(one(T)),)
  b = (II.Infinity(one(T)),)
  II.Domain.Box(a, b)
end

function ∫dom(plane::Plane)
  T = numtype(lentype(plane))
  a = (-II.Infinity(one(T)), -II.Infinity(one(T)))
  b = (II.Infinity(one(T)), II.Infinity(one(T)))
  II.Domain.Box(a, b)
end
