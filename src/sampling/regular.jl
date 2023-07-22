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

function sample(::AbstractRNG, geom::Geometry{Dim,T}, method::RegularSampling) where {Dim,T}
  V = T <: AbstractFloat ? T : Float64
  D = paramdim(geom)
  sz = fitdims(method.sizes, D)
  δₛ = soffset(geom)
  δₑ = eoffset(geom)
  tₛ = ntuple(i -> V(0 + δₛ[i](sz[i])), D)
  tₑ = ntuple(i -> V(1 - δₑ[i](sz[i])), D)
  rs = (range(tₛ[i], stop=tₑ[i], length=sz[i]) for i in 1:D)
  iᵣ = (geom(uv...) for uv in Iterators.product(rs...))
  iₚ = (p for p in extrapoints(geom))
  Iterators.flatmap(identity, (iᵣ, iₚ))
end

soffset(::Sphere{3}) = (n -> inv(n + 1), n -> zero(n))
eoffset(::Sphere{3}) = (n -> inv(n + 1), n -> inv(n))
extrapoints(s::Sphere{3}) = (s(0, 0), s(1, 0))

soffset(b::Ball) = (n -> inv(n + 1), soffset(boundary(b))...)
eoffset(b::Ball) = (n -> zero(n), eoffset(boundary(b))...)
extrapoints(b::Ball) = (center(b),)

soffset(g::Geometry) = ntuple(i -> (n -> zero(n)), paramdim(g))
eoffset(g::Geometry) = ntuple(i -> (n -> isperiodic(g)[i] ? inv(n) : zero(n)), paramdim(g))
extrapoints(::Geometry) = ()

# --------------
# SPECIAL CASES
# --------------

function sample(::AbstractRNG, cylsurf::CylinderSurface{T}, method::RegularSampling) where {T}
  V = T <: AbstractFloat ? T : Float64
  sz = fitdims(method.sizes, paramdim(cylsurf))
  c = center(cylsurf)
  r = radius(cylsurf)
  b = bottom(cylsurf)
  t = top(cylsurf)
  a = axis(cylsurf)

  φmin, φmax = V(0), V(2π)
  zmin, zmax = V(0), V(1)
  δφ = (φmax - φmin) / sz[1]
  φs = range(φmin, stop=φmax - δφ, length=sz[1])
  zs = range(zmin, stop=zmax, length=sz[2])

  # rotation to align z axis with cylinder axis
  d₃ = a(1) - a(0)
  l = norm(d₃)
  d₃ /= l
  d₁, d₂ = householderbasis(d₃)
  R = transpose([d₁ d₂ d₃])

  # new normals of planes in new rotated system
  nᵦ = R * normal(b)
  nₜ = R * normal(t)

  # given cylindrical coordinates (r*cos(φ), r*sin(φ), z) and the
  # equation of the plane, we can solve for z and find all points
  # along the ellipse obtained by intersection
  zᵦ(φ) = -l / 2 - (r * cos(φ) * nᵦ[1] + r * sin(φ) * nᵦ[2]) / nᵦ[3]
  zₜ(φ) = +l / 2 - (r * cos(φ) * nₜ[1] + r * sin(φ) * nₜ[2]) / nₜ[3]
  cᵦ(φ) = Point(r * cos(φ), r * sin(φ), zᵦ(φ))
  cₜ(φ) = Point(r * cos(φ), r * sin(φ), zₜ(φ))

  function point(φ, z)
    pᵦ, pₜ = cᵦ(φ), cₜ(φ)
    p = pᵦ + z * (pₜ - pᵦ)
    Point(R' * coordinates(p) + coordinates(c))
  end

  iᵣ = (point(φ, z) for z in zs for φ in φs)
  iₚ = (b(0, 0), t(0, 0))
  Iterators.flatmap(identity, (iᵣ, iₚ))
end

function sample(rng::AbstractRNG, grid::CartesianGrid, method::RegularSampling)
  sample(rng, boundingbox(grid), method)
end
