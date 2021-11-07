# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MetricBall(radii, [angles]; convention=TaitBryanExtr)
    MetricBall(radius, metric=Euclidean())

A metric ball is a neighborhood that can be expressed in terms
of a metric and a set of `radii`. The two main examples are the
Euclidean ball an the Mahalanobis (ellipsoid) ball.

When multiple `radii` are provided, they can be rotated by `angles` 
according to a given rotation `convention`. Alternatively, a metric
from the [Distances.jl](https://github.com/JuliaStats/Distances.jl)
package can be specified together with a single `radius`.

## Examples

N-dimensional Euclidean ball with radius `1.0`:

```julia
julia> euclidean = MetricBall(1.0)
```

Axis-aligned 3D ellispoid with radii `[3.0, 2.0, 1.0]`:

```julia
julia> mahalanobis = MetricBall([3.0, 2.0, 1.0])
```

See also [`mahalanobis`](@ref).
"""
struct MetricBall{R,M} <: Neighborhood
  radii::R
  metric::M
end

function MetricBall(radii::AbstractVector{T}, angles=nothing;
                    convention=TaitBryanExtr) where {T<:Real}
  ndims  = length(radii)
  angles = isnothing(angles) ? (ndims == 2 ? [zero(T)] : zeros(T, 3)) : angles
  metric = mahalanobis(radii, angles, convention)
  MetricBall{typeof(radii),typeof(metric)}(radii, metric)
end

function MetricBall(radius::T, metric=Euclidean()) where {T<:Real}
  MetricBall{typeof(radius),typeof(metric)}(radius, metric)
end

"""
    radii(ball)

Return the radii of the metric `ball`.
"""
radii(ball::MetricBall) = ball.radii

"""
    metric(ball)

Return the metric of the metric `ball`.
"""
metric(ball::MetricBall) = ball.metric

"""
    range(ball)

Return the range of the metric `ball`, i.e.
the smallest value `r` such that `||v|| ≤ r`
for any `v ∈ ball`.
"""
Base.range(ball::MetricBall{<:Real}) = ball.radii
Base.range(::MetricBall{R,<:Mahalanobis}) where {R} = one(eltype(R))

function Base.show(io::IO, ball::MetricBall)
  r = length(ball.radii) > 1 ? Tuple(ball.radii) : ball.radii
  m = nameof(typeof(ball.metric))
  print(io, "MetricBall($r, $m)")
end