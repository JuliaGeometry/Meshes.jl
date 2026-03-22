# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# default integration method
const HADAPTIVE = II.Backend.HAdaptiveIntegration()

"""
    integral(fun, geom[, method])

Calculate the integral over the `geom`etry of the `fun`ction that maps
[`Point`](@ref)s to values in a linear space.

    integral(fun, dom[, method])

Alternatively, calculate the integral over the `dom`ain (e.g., mesh) by
summing the integrals for each constituent geometry.

By default, use h-adaptive integration for good accuracy on a wide range of geometries.

See also [`localintegral`](@ref).
"""
integral(fun, geom::Geometry, method=HADAPTIVE) = _integral(fun, geom, method)

# cylinder surface is the union of the curved surface and the top and bottom disks
integral(fun, cylsurf::CylinderSurface, method=HADAPTIVE) =
  localintegral(fun ∘ cylsurf, cylsurf, method) +
  integral(fun, top(cylsurf), method) +
  integral(fun, bottom(cylsurf), method)

# cone surface is the union of the curved surface and the base disk
integral(fun, conesurf::ConeSurface, method=HADAPTIVE) =
  localintegral(fun ∘ conesurf, conesurf, method) + integral(fun, base(conesurf), method)

# frustum surface is the union of the curved surface and the top and bottom disks
integral(fun, frustumsurf::FrustumSurface, method=HADAPTIVE) =
  localintegral(fun ∘ frustumsurf, frustumsurf, method) +
  integral(fun, top(frustumsurf), method) +
  integral(fun, bottom(frustumsurf), method)

# rope is the union of its constituent segments
integral(fun, rope::Rope, method=HADAPTIVE) = sum(integral(fun, seg, method) for seg in segments(rope))

# ring is the union of its constituent segments
integral(fun, ring::Ring, method=HADAPTIVE) = sum(integral(fun, seg, method) for seg in segments(ring))

# polygon is the union of its constituent ngons
integral(fun, poly::Polygon, method=HADAPTIVE) = sum(integral(fun, ngon, method) for ngon in discretize(poly))

# integrate triangles with local integration
integral(fun, tri::Triangle, method=HADAPTIVE) = _integral(fun, tri, method)

# integrate quadrangle with local integration
integral(fun, quad::Quadrangle, method=HADAPTIVE) = _integral(fun, quad, method)

# multi-geometry is the union of its constituent geometries
integral(fun, multi::Multi, method=HADAPTIVE) = sum(integral(fun, geom, method) for geom in parent(multi))

# domain is the union of its constituent geometries
integral(fun, dom::Domain, method=HADAPTIVE) = sum(integral(fun, geom, method) for geom in dom)

# fallback to local integration of fun ∘ geom
_integral(fun, geom, method) = localintegral(fun ∘ geom, geom, method)

"""
    localintegral(fun, geom[, method])

Calculate the integral over the `geom`etry of the `fun`ction that maps
parametric coordinates `uvw` to values in a linear space.

By default, use h-adaptive integration for good accuracy on a wide range of geometries.

See also [`integral`](@ref).
"""
function localintegral(fun, geom::Geometry, method=HADAPTIVE)
  # integrand is equal to function times differential element
  integrand(uvw...) = fun(uvw...) * differential(geom, uvw)

  # domain of integration for the given geometry
  domain = ∫domain(geom)

  # extract units of integral by assuming
  # integrand can be evaluated at zeros
  N = paramdim(geom)
  T = numtype(lentype(geom))
  o = ntuple(_ -> zero(T), N)
  u = unit.(integrand(o...))

  # strip units to help integration backends
  f(uvw...) = ustrip.(integrand(uvw...))

  # perform numerical integration
  II.integral(f, domain; backend=method) .* u
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

function ∫domain(tri::Triangle)
  T = numtype(lentype(tri))
  a = (zero(T), zero(T))
  b = (one(T), zero(T))
  c = (zero(T), one(T))
  II.Domain.Simplex(a, b, c)
end

function ∫domain(tetra::Tetrahedron)
  T = numtype(lentype(tetra))
  a = (zero(T), zero(T), zero(T))
  b = (one(T), zero(T), zero(T))
  c = (zero(T), one(T), zero(T))
  d = (zero(T), zero(T), one(T))
  II.Domain.Simplex(a, b, c, d)
end
