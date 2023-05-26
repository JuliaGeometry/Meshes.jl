# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

default_rotation(::Val{2}, T) = one(Angle2d{T})
default_rotation(::Val{3}, T) = one(QuatRotation{T})

"""
    MetricBall(radii, rotation=nothing)
    MetricBall(radius, metric=Euclidean())

A metric ball is a neighborhood that can be expressed in terms
of a metric and a set of `radii`. The two main examples are the
Euclidean ball an the Mahalanobis (ellipsoid) ball.

When multiple `radii` are provided, they can be rotated by a
`rotation` specification from the [Rotations.jl]
(https://github.com/JuliaGeometry/Rotations.jl)
package. Alternatively, a metric from the [Distances.jl]
(https://github.com/JuliaStats/Distances.jl) package can
be specified together with a single `radius`.

## Examples

N-dimensional Euclidean ball with radius `1.0`:

```julia
julia> euclidean = MetricBall(1.0)
```

Axis-aligned 3D ellispoid with radii `(3.0, 2.0, 1.0)`:

```julia
julia> mahalanobis = MetricBall((3.0, 2.0, 1.0))
```
"""
struct MetricBall{R,M} <: Neighborhood
  radii::R
  metric::M
end

function MetricBall(radii::SVector{Dim,T}, R=default_rotation(Val{Dim}(), T)) where {Dim,T}
  # scaling matrix
  Λ = Diagonal(one(T) ./ radii .^ 2)

  # Mahalanobis metric
  metric = Mahalanobis(Symmetric(R' * Λ * R))

  MetricBall{typeof(radii),typeof(metric)}(radii, metric)
end

MetricBall(radii::NTuple{Dim,T}, rotation=default_rotation(Val{Dim}(), T)) where {Dim,T} =
  MetricBall(SVector(radii), rotation)

# avoid silent calls to inner constructor
MetricBall(radii::AbstractVector{T}, rotation=default_rotation(Val{Dim}(), T)) where {T} =
  MetricBall(SVector{length(radii),T}(radii), rotation)

function MetricBall(radius::T, metric=Euclidean()) where {T<:Real}
  radii = SVector(radius)
  MetricBall{typeof(radii),typeof(metric)}(radii, metric)
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
    radius(ball)

Return the effective radius of the metric `ball`,
i.e. the value `r` such that `||v|| ≤ r, ∀ v ∈ ball`
and `||v|| > r, ∀ v ∉ ball``.
"""
radius(ball::MetricBall) = first(ball.radii)
radius(::MetricBall{R,<:Mahalanobis}) where {R} = one(eltype(R))

"""
    isisotropic(ball)

Tells whether or not the metric `ball` is isotropic,
i.e. if all its radii are equal.
"""
isisotropic(ball::MetricBall) = length(unique(ball.radii)) == 1

function Base.show(io::IO, ball::MetricBall)
  n = length(ball.radii)
  r = n > 1 ? Tuple(ball.radii) : first(ball.radii)
  m = nameof(typeof(ball.metric))
  print(io, "MetricBall($r, $m)")
end
