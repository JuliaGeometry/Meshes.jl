# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# default integration method
const HADAPTIVE = II.Backend.HAdaptiveIntegration(rtol=1e-3)

"""
    integral(fun, geom; ∫backend, ∂backend)

Calculate the integral over the `geom`etry of the `fun`ction that maps
[`Point`](@ref)s to values in a linear space using an integration `∫backend`
from IntegrationInterface.jl and a differentiation `∂backend` from
DifferentiationInterface.jl.

    integral(fun, dom; ∫backend, ∂backend)

Alternatively, calculate the integral over the `dom`ain (e.g., mesh) by
summing the integrals for each constituent geometry.

By default, `∫backend` is set to h-adaptive integration for good accuracy
across a wide range of geometries and `∂backend` is set to finite differences.

See also [`localintegral`](@ref).
"""
integral(fun, geom::Geometry; ∫backend=HADAPTIVE) = _integral(fun, geom, ∫backend)

# cylinder surface is the union of the curved surface and the top and bottom disks
integral(fun, cylsurf::CylinderSurface; ∫backend=HADAPTIVE) =
  localintegral(fun ∘ cylsurf, cylsurf; ∫backend) +
  integral(fun, top(cylsurf); ∫backend) +
  integral(fun, bottom(cylsurf); ∫backend)

# cone surface is the union of the curved surface and the base disk
integral(fun, conesurf::ConeSurface; ∫backend=HADAPTIVE) =
  localintegral(fun ∘ conesurf, conesurf; ∫backend) + integral(fun, base(conesurf); ∫backend)

# frustum surface is the union of the curved surface and the top and bottom disks
integral(fun, frustumsurf::FrustumSurface; ∫backend=HADAPTIVE) =
  localintegral(fun ∘ frustumsurf, frustumsurf; ∫backend) +
  integral(fun, top(frustumsurf); ∫backend) +
  integral(fun, bottom(frustumsurf); ∫backend)

# rope is the union of its constituent segments
integral(fun, rope::Rope; ∫backend=HADAPTIVE) = sum(integral(fun, seg; ∫backend) for seg in segments(rope))

# ring is the union of its constituent segments
integral(fun, ring::Ring; ∫backend=HADAPTIVE) = sum(integral(fun, seg; ∫backend) for seg in segments(ring))

# polygon is the union of its constituent ngons
integral(fun, poly::Polygon; ∫backend=HADAPTIVE) = sum(integral(fun, ngon; ∫backend) for ngon in discretize(poly))

# integrate triangles with local integration
integral(fun, tri::Triangle; ∫backend=HADAPTIVE) = _integral(fun, tri, ∫backend)

# integrate quadrangle with local integration
integral(fun, quad::Quadrangle; ∫backend=HADAPTIVE) = _integral(fun, quad, ∫backend)

# multi-geometry is the union of its constituent geometries
integral(fun, multi::Multi; ∫backend=HADAPTIVE) = sum(integral(fun, geom; ∫backend) for geom in parent(multi))

# domain is the union of its constituent geometries
integral(fun, dom::Domain; ∫backend=HADAPTIVE) = sum(integral(fun, geom; ∫backend) for geom in dom)

# fallback to local integration of fun ∘ geom
_integral(fun, geom, ∫backend) = localintegral(fun ∘ geom, geom; ∫backend)

"""
    localintegral(fun, geom; ∫backend)

Calculate the integral over the `geom`etry of the `fun`ction that maps
parametric coordinates `uvw` to values in a linear space using an
integration `∫backend` from IntegrationInterface.jl.

By default, `∫backend` is set to h-adaptive integration for good accuracy
across a wide range of geometries and `∂backend` is set to finite differences.

See also [`integral`](@ref).
"""
function localintegral(fun, geom::Geometry; ∫backend=HADAPTIVE)
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
  II.integral(f, domain; backend=∫backend) .* u
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
