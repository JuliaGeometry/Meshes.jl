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

RegularSampling(sizes::Vararg{Int,N}) where {N} =
  RegularSampling(sizes)

function sample(::AbstractRNG, geom::Geometry{Dim,T},
                method::RegularSampling) where {Dim,T}
  V  = T <: AbstractFloat ? T : Float64
  sz = fitdims(method.sizes, paramdim(geom))
  rs = (range(V(0), V(1), length=s) for s in sz)
  ivec(geom(uv...) for uv in Iterators.product(rs...))
end

function sample(::AbstractRNG, sphere::Sphere{2,T},
                method::RegularSampling) where {T}
  V  = T <: AbstractFloat ? T : Float64
  sz = fitdims(method.sizes, paramdim(sphere))
  c, r = center(sphere), radius(sphere)

  θmin, θmax = V(0), V(2π)
  δθ = (θmax - θmin) / sz[1]
  θrange = range(θmin, stop=θmax-δθ, length=sz[1])

  r⃗(θ) = Vec{2,V}(r*cos(θ), r*sin(θ))

  ivec(c + r⃗(θ) for θ in θrange)
end

# spherical coordinates in ISO 80000-2:2019 convention
function sample(::AbstractRNG, sphere::Sphere{3,T},
                method::RegularSampling) where {T}
  V  = T <: AbstractFloat ? T : Float64
  sz = fitdims(method.sizes, paramdim(sphere))
  c, r = center(sphere), radius(sphere)

  θmin, θmax = V(0), V(π)
  φmin, φmax = V(0), V(2π)
  δθ = (θmax - θmin) / (sz[1] + 1)
  δφ = (φmax - φmin) / (sz[2]    )
  θrange = range(θmin+δθ, stop=θmax-δθ, length=sz[1])
  φrange = range(φmin, stop=φmax-δφ, length=sz[2])

  r⃗(θ, φ) = Vec{3,T}(r*sin(θ)*cos(φ), r*sin(θ)*sin(φ), r*cos(θ))

  ivec(c + r⃗(θ, φ) for θ in θrange, φ in φrange)
end

function sample(::AbstractRNG, ball::Ball{Dim,T},
                method::RegularSampling) where {Dim,T}
  V  = T <: AbstractFloat ? T : Float64
  sz = fitdims(method.sizes, paramdim(ball))
  c, r = center(ball), radius(ball)

  smin, smax = V(0), V(1)
  δs = (smax - smin) / (last(sz) - 1)
  srange = range(smin+δs, stop=smax, length=last(sz))

  # reuse samples on the boundary
  points = sample(Sphere(c, r), RegularSampling(sz[1:Dim-1]))

  scale(p, s) = Point(s * coordinates(p))

  ivec(scale(p, s) for p in points, s in srange)
end

function sample(::AbstractRNG, cylsurf::CylinderSurface{T},
                method::RegularSampling) where {T}
  V  = T <: AbstractFloat ? T : Float64
  sz = fitdims(method.sizes, paramdim(cylsurf))
  r  = radius(cylsurf)
  b  = bottom(cylsurf)
  t  = top(cylsurf)
  a  = axis(cylsurf)

  φmin, φmax = V(0), V(2π)
  zmin, zmax = V(0), V(1)
  δφ = (φmax - φmin) / sz[1]
  φrange = range(φmin, stop=φmax-δφ, length=sz[1])
  zrange = range(zmin, stop=zmax,    length=sz[2])

  # rotation to align z axis with cylinder axis
  d₃  = a(1) - a(0)
  l  = norm(d₃)
  d₃ /= l
  d₁, d₂ = householderbasis(d₃)
  R = transpose([d₁ d₂ d₃])

  # new normals of planes in new rotated system
  nᵦ = R * normal(b)
  nₜ = R * normal(t)

  # given cylindrical coordinates (r*cos(φ), r*sin(φ), z) and the
  # equation of the plane, we can solve for z and find all points
  # along the ellipse obtained by intersection
  zᵦ(φ) = -l/2 - (r*cos(φ)*nᵦ[1] + r*sin(φ)*nᵦ[2]) / nᵦ[3]
  zₜ(φ) = +l/2 - (r*cos(φ)*nₜ[1] + r*sin(φ)*nₜ[2]) / nₜ[3]
  cᵦ(φ) = Point(r*cos(φ), r*sin(φ), zᵦ(φ))
  cₜ(φ) = Point(r*cos(φ), r*sin(φ), zₜ(φ))

  # center of cylinder for final translation
  oᵦ = coordinates(b(0, 0))
  oₜ = coordinates(t(0, 0))
  oₘ = @. (oᵦ + oₜ) / 2

  function point(φ, z)
    pᵦ, pₜ = cᵦ(φ), cₜ(φ)
    p = pᵦ + z*(pₜ - pᵦ)
    Point((R' * coordinates(p)) + oₘ)
  end

  ivec(point(φ, z) for φ in φrange, z in zrange)
end

function sample(rng::AbstractRNG, grid::CartesianGrid,
                method::RegularSampling)
  sample(rng, boundingbox(grid), method)
end

function sample(::AbstractRNG, torus::Torus{T},
                method::RegularSampling) where {T}
  V  = T <: AbstractFloat ? T : Float64
  sz = fitdims(method.sizes, paramdim(torus))
  R, r = radii(torus)

  kxy = R^2 - r^2
  kz = √kxy * r

  umin, umax = V(-π), V(π)
  vmin, vmax = V(-π), V(π)
  δu = (umax - umin) / sz[1]
  δv = (vmax - vmin) / sz[2]
  urange = range(umin, stop=umax-δu, length=sz[1])
  vrange = range(vmin, stop=vmax-δv, length=sz[2])

  c = center(torus)
  n⃗ = normal(torus)
  M = uvrotation(n⃗, Vec{3,T}(0, 0, 1))

  r⃗(u, v) = Vec{3,T}(kxy * cos(u), kxy * sin(u), kz * sin(v)) / (R - r*cos(v))

  ivec(c + M * r⃗(u, v) for u in urange, v in vrange)
end