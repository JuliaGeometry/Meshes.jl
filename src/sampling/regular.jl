# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularSampling(n1, n2, ..., np)

Generate samples regularly using `n1` points along the first
parametric dimension, `n2` points along the second parametric
dimension, ..., `np` poitns along the last parametric dimension.

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

function sample(::AbstractRNG, box::Box,
                method::RegularSampling)
  sz = fitdims(method.sizes, paramdim(box))
  l, u = extrema(box)

  # origin and spacing
  or, sp = l, (u - l) ./ (sz .- 1)

  ivec(or + (ind.I .- 1) .* sp for ind in CartesianIndices(sz))
end

function sample(::AbstractRNG, sphere::Sphere{2,T},
                method::RegularSampling) where {T}
  sz = fitdims(method.sizes, paramdim(sphere))
  c, r = center(sphere), radius(sphere)

  V = T <: AbstractFloat ? T : Float64
  θmin, θmax = V(0), V(2π)
  δθ = (θmax - θmin) / sz[1]
  θrange = range(θmin, stop=θmax-δθ, length=sz[1])

  r⃗(θ) = Vec{2,V}(r*cos(θ), r*sin(θ))

  ivec(c + r⃗(θ) for θ in θrange)
end

# spherical coordinates in ISO 80000-2:2019 convention
function sample(::AbstractRNG, sphere::Sphere{3,T},
                method::RegularSampling) where {T}
  sz = fitdims(method.sizes, paramdim(sphere))
  c, r = center(sphere), radius(sphere)

  V = T <: AbstractFloat ? T : Float64
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
  sz = fitdims(method.sizes, paramdim(ball))
  c, r = center(ball), radius(ball)

  V = T <: AbstractFloat ? T : Float64
  smin, smax = V(0), V(1)
  δs = (smax - smin) / (last(sz) - 1)
  srange = range(smin+δs, stop=smax, length=last(sz))

  # reuse samples on the boundary
  points = sample(Sphere(c, r), RegularSampling(sz[1:Dim-1]))

  scale(p, s) = Point(s * coordinates(p))

  ivec(scale(p, s) for p in points, s in srange)
end

function sample(::AbstractRNG, seg::Segment{Dim,T},
                method::RegularSampling) where {Dim,T}
  sz = fitdims(method.sizes, paramdim(seg))
  trange = range(T(0), T(1), length=sz[1])
  (seg(t) for t in trange)
end

function sample(::AbstractRNG, quad::Quadrangle{Dim,T},
                method::RegularSampling) where {Dim,T}
  sz = fitdims(method.sizes, paramdim(quad))
  urange = range(T(0), T(1), length=sz[1])
  vrange = range(T(0), T(1), length=sz[2])
  ivec(quad(u, v) for u in urange, v in vrange)
end

function sample(::AbstractRNG, hex::Hexahedron{Dim,T},
                method::RegularSampling) where {Dim,T}
  sz = fitdims(method.sizes, paramdim(hex))
  urange = range(T(0), T(1), length=sz[1])
  vrange = range(T(0), T(1), length=sz[2])
  wrange = range(T(0), T(1), length=sz[3])
  ivec(hex(u, v, w) for u in urange, v in vrange, w in wrange)
end
