# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# default integration backend
hadaptive(geom) = II.Backend.HAdaptiveIntegration(rtol=rtol(numtype(lentype(geom))))

"""
    integral(fun, geom; ibackend, dbackend)

Calculate the integral over the `geom`etry of the `fun`ction that maps
[`Point`](@ref)s to values in a linear space using an integration `ibackend`
from IntegrationInterface.jl and a differentiation `dbackend` from
DifferentiationInterface.jl.

    integral(fun, dom; ibackend, dbackend)

Alternatively, calculate the integral over the `dom`ain (e.g., mesh) by
summing the integrals for each constituent geometry.

By default, `ibackend` is set to h-adaptive integration for good accuracy
across a wide range of geometries and `dbackend` is set to finite differences.

See also [`localintegral`](@ref).
"""
integral(fun, geom::Geometry; ibackend=hadaptive(geom), dbackend=FINITEDIFF) = _integral(fun, geom, ibackend, dbackend)

# cylinder surface is the union of the curved surface and the top and bottom disks
integral(fun, cylsurf::CylinderSurface; ibackend=hadaptive(cylsurf), dbackend=FINITEDIFF) =
  localintegral(fun ∘ cylsurf, cylsurf; ibackend, dbackend) +
  integral(fun, top(cylsurf); ibackend, dbackend) +
  integral(fun, bottom(cylsurf); ibackend, dbackend)

# cone surface is the union of the curved surface and the base disk
integral(fun, conesurf::ConeSurface; ibackend=hadaptive(conesurf), dbackend=FINITEDIFF) =
  localintegral(fun ∘ conesurf, conesurf; ibackend, dbackend) + integral(fun, base(conesurf); ibackend, dbackend)

# frustum surface is the union of the curved surface and the top and bottom disks
integral(fun, frustumsurf::FrustumSurface; ibackend=hadaptive(frustumsurf), dbackend=FINITEDIFF) =
  localintegral(fun ∘ frustumsurf, frustumsurf; ibackend, dbackend) +
  integral(fun, top(frustumsurf); ibackend, dbackend) +
  integral(fun, bottom(frustumsurf); ibackend, dbackend)

# rope is the union of its constituent segments
integral(fun, rope::Rope; ibackend=hadaptive(rope), dbackend=FINITEDIFF) =
  sum(integral(fun, seg; ibackend, dbackend) for seg in segments(rope))

# ring is the union of its constituent segments
integral(fun, ring::Ring; ibackend=hadaptive(ring), dbackend=FINITEDIFF) =
  sum(integral(fun, seg; ibackend, dbackend) for seg in segments(ring))

# polygon is the union of its constituent ngons
integral(fun, poly::Polygon; ibackend=hadaptive(poly), dbackend=FINITEDIFF) =
  sum(integral(fun, ngon; ibackend, dbackend) for ngon in discretize(poly))

# integrate triangles with local integration
integral(fun, tri::Triangle; ibackend=hadaptive(tri), dbackend=FINITEDIFF) = _integral(fun, tri, ibackend, dbackend)

# integrate quadrangle with local integration
integral(fun, quad::Quadrangle; ibackend=hadaptive(quad), dbackend=FINITEDIFF) = _integral(fun, quad, ibackend, dbackend)

# multi-geometry is the union of its constituent geometries
integral(fun, multi::Multi; ibackend=hadaptive(multi), dbackend=FINITEDIFF) =
  sum(integral(fun, geom; ibackend, dbackend) for geom in parent(multi))

# domain is the union of its constituent geometries
integral(fun, dom::Domain; ibackend=hadaptive(dom), dbackend=FINITEDIFF) =
  sum(integral(fun, geom; ibackend, dbackend) for geom in dom)

# fallback to local integration of fun ∘ geom
_integral(fun, geom, ibackend, dbackend) = localintegral(fun ∘ geom, geom; ibackend, dbackend)

"""
    localintegral(fun, geom; ibackend, dbackend)

Calculate the integral over the `geom`etry of the `fun`ction that maps
parametric coordinates `uvw` to values in a linear space using an integration
`ibackend` from IntegrationInterface.jl and a differentiation `dbackend`
from DifferentiationInterface.jl.

By default, `ibackend` is set to h-adaptive integration for good accuracy
across a wide range of geometries and `dbackend` is set to finite differences.

See also [`integral`](@ref).
"""
function localintegral(fun, geom::Geometry; ibackend=hadaptive(geom), dbackend=FINITEDIFF)
  # integrand is equal to function times differential element
  integrand(uvw...) = fun(uvw...) * differential(geom, uvw; dbackend)

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
  II.integral(f, domain; backend=ibackend) .* u
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
