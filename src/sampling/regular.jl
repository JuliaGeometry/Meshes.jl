# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularSampling(n1, n2, ..., np)

Generate samples regularly using `n1` points along the first
parametric dimension, `n2` points along the second parametric
dimension, ..., `np` points along the last parametric dimension.

## Examples

Sample sphere regularly with 360 longitudes and 180 latitudes:

```julia
sample(Sphere((0,0,0), 1), RegularSampling(360, 180))
```
"""
struct RegularSampling{N} <: ContinuousSamplingMethod
  sizes::Dims{N}
end

RegularSampling(sizes::Vararg{Int,N}) where {N} = RegularSampling(sizes)

function sample(::AbstractRNG, geom::Geometry, method::RegularSampling)
  T = numtype(lentype(geom))
  D = paramdim(geom)
  sz = fitdims(method.sizes, D)
  δₛ = firstoffset(geom)
  δₑ = lastoffset(geom)
  tₛ = ntuple(i -> T(0 + δₛ[i](sz[i])), D)
  tₑ = ntuple(i -> T(1 - δₑ[i](sz[i])), D)
  rs = (range(tₛ[i], stop=tₑ[i], length=sz[i]) for i in 1:D)
  iᵣ = (geom(uv...) for uv in Iterators.product(rs...))
  iₚ = (p for p in extrapoints(geom, sz))
  Iterators.flatmap(identity, (iᵣ, iₚ))
end

firstoffset(g::Geometry) = _firstoffset(g, Val(embeddim(g)))
lastoffset(g::Geometry) = _lastoffset(g, Val(embeddim(g)))
extrapoints(g::Geometry, sz) = _extrapoints(g, Val(embeddim(g)), sz)

_firstoffset(g::Geometry, ::Val) = ntuple(i -> (n -> zero(n)), paramdim(g))
_lastoffset(g::Geometry, ::Val) = ntuple(i -> (n -> isperiodic(g)[i] ? inv(n) : zero(n)), paramdim(g))
_extrapoints(::Geometry, ::Val, sz) = ()

firstoffset(d::Disk) = (n -> inv(n), firstoffset(boundary(d))...)
lastoffset(d::Disk) = (n -> zero(n), lastoffset(boundary(d))...)
extrapoints(d::Disk, sz) = (center(d),)

firstoffset(b::Ball) = (n -> inv(n), firstoffset(boundary(b))...)
lastoffset(b::Ball) = (n -> zero(n), lastoffset(boundary(b))...)
extrapoints(b::Ball, sz) = (center(b),)

_firstoffset(::Sphere, ::Val{3}) = (n -> inv(n + 1), n -> zero(n))
_lastoffset(::Sphere, ::Val{3}) = (n -> inv(n + 1), n -> inv(n))
_extrapoints(s::Sphere, ::Val{3}, sz) = (s(0, 0), s(1, 0))

firstoffset(::Ellipsoid) = (n -> inv(n + 1), n -> zero(n))
lastoffset(::Ellipsoid) = (n -> inv(n + 1), n -> inv(n))
extrapoints(e::Ellipsoid, sz) = (e(0, 0), e(1, 0))

firstoffset(::Cylinder) = (n -> inv(n), n -> zero(n), n -> zero(n))
lastoffset(::Cylinder) = (n -> zero(n), n -> inv(n), n -> zero(n))
function extrapoints(c::Cylinder, sz)
  T = numtype(lentype(c))
  b = bottom(c)(0, 0)
  t = top(c)(0, 0)
  s = Segment(b, t)
  [s(t) for t in range(zero(T), one(T), sz[3])]
end

firstoffset(::CylinderSurface) = (n -> zero(n), n -> zero(n))
lastoffset(::CylinderSurface) = (n -> inv(n), n -> zero(n))
extrapoints(c::CylinderSurface, sz) = (bottom(c)(0, 0), top(c)(0, 0))

firstoffset(::ConeSurface) = (n -> zero(n), n -> inv(n))
lastoffset(::ConeSurface) = (n -> inv(n), n -> zero(n))
extrapoints(c::ConeSurface, sz) = (apex(c), base(c)(0, 0))

firstoffset(::FrustumSurface) = (n -> zero(n), n -> zero(n))
lastoffset(::FrustumSurface) = (n -> inv(n), n -> zero(n))
extrapoints(c::FrustumSurface, sz) = (bottom(c)(0, 0), top(c)(0, 0))

# --------------
# SPECIAL CASES
# --------------

function sample(rng::AbstractRNG, grid::CartesianGrid, method::RegularSampling)
  sample(rng, boundingbox(grid), method)
end
