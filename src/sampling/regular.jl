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

soffset(::CylinderSurface) = (n -> zero(n), n -> zero(n))
eoffset(::CylinderSurface) = (n -> inv(n), n -> zero(n))
extrapoints(c::CylinderSurface) = (bottom(c)(0, 0), top(c)(0, 0))

soffset(g::Geometry) = ntuple(i -> (n -> zero(n)), paramdim(g))
eoffset(g::Geometry) = ntuple(i -> (n -> isperiodic(g)[i] ? inv(n) : zero(n)), paramdim(g))
extrapoints(::Geometry) = ()

# --------------
# SPECIAL CASES
# --------------

function sample(rng::AbstractRNG, grid::CartesianGrid, method::RegularSampling)
  sample(rng, boundingbox(grid), method)
end
