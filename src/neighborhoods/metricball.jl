# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MetricBall(radius, metric=EuclideanDistance())

A metric ball expressed in terms of a single `radius` and a
`metric` (or distance) between geometries.

    MetricBall(radii, rotation=I)

A Mahalanobis ball (a.k.a. ellipsoid) expressed in terms
of multiple `radii` and a `rotation` specification from the
[Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl)
package.

## Examples

N-dimensional Euclidean ball with radius `1.0`:

```julia
julia> MetricBall(1.0)
```

Axis-aligned 3D ellipsoid with radii `(3.0, 2.0, 1.0)`:

```julia
julia> MetricBall((3.0, 2.0, 1.0))
```
"""
struct MetricBall{Dim,ℒ<:Len,R,M} <: Neighborhood
  radii::NTuple{Dim,ℒ}
  rotation::R
  metric::M
end

function MetricBall(radii, rotation, metric)
  r = promote(map(float ∘ aslen, radii)...)
  R = rotation
  m = metric
  MetricBall{length(r),eltype(r),typeof(R),typeof(m)}(r, R, m)
end

MetricBall(radius::Number, metric=EuclideanDistance()) = MetricBall((aslen(radius),), I, metric)

MetricBall(radii::Tuple, rotation=I) = MetricBall(radii, rotation, MahalanobisDistance(radii, rotation))

"""
    hasequalradii(ball)

Tells whether or not the metric `ball` has equal radii.
"""
hasequalradii(ball::MetricBall) = allequal(ball.radii)

"""
    radii(ball)

Return the radii of the metric `ball`.
"""
radii(ball::MetricBall) = ball.radii

"""
    rotation(ball)

Return the rotation of the metric `ball`.
"""
rotation(ball::MetricBall) = ball.rotation

"""
    metric(ball)

Return the metric of the metric `ball`.
"""
metric(ball::MetricBall) = ball.metric

"""
    radius(ball)

Return the effective radius of the metric `ball`,
i.e. the value `r` such that `||v|| ≤ r, ∀ v ∈ ball`
and `||v|| > r, ∀ v ∉ ball``.
"""
function radius(ball::MetricBall)
  r = first(ball.radii)
  ball.metric isa MahalanobisDistance ? oneunit(r) : r
end

function *(α::Real, ball::MetricBall)
  if ball.metric isa MahalanobisDistance
    MetricBall(α .* ball.radii, ball.rotation)
  else
    MetricBall(α .* ball.radii, I, ball.metric)
  end
end

function Base.show(io::IO, ball::MetricBall)
  n = length(ball.radii)
  r = n > 1 ? ball.radii : first(ball.radii)
  m = nameof(typeof(ball.metric))
  print(io, "MetricBall($r, $m)")
end
