# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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

Axis-aligned 3D ellipsoid with radii `(3.0, 2.0, 1.0)`:

```julia
julia> mahalanobis = MetricBall((3.0, 2.0, 1.0))
```
"""
struct MetricBall{ℒ,R,M} <: Neighborhood
  radii::ℒ
  rotation::R

  # state fields
  metric::M
end

function MetricBall(radii::SVector{Dim,T}, rotation=default_rotation(Val{Dim}(), T)) where {Dim,T}
  # scaling matrix
  Λ = Diagonal(one(T) ./ radii .^ 2)

  # rotation matrix
  R = rotation

  # anisotropy matrix
  M = Symmetric(R * Λ * R')

  # Mahalanobis metric
  metric = Mahalanobis(M)

  MetricBall(radii, rotation, metric)
end

MetricBall(radii::NTuple{Dim,T}, rotation=default_rotation(Val{Dim}(), T)) where {Dim,T} =
  MetricBall(SVector(radii), rotation)

# avoid silent calls to inner constructor
MetricBall(radii::AbstractVector{T}, rotation=default_rotation(Val{length(radii)}(), T)) where {T} =
  MetricBall(SVector{length(radii),T}(radii), rotation)

MetricBall(radius::T, metric=Euclidean()) where {T<:Number} = MetricBall(SVector(radius), nothing, metric)

default_rotation(::Val{2}, T) = one(Angle2d{T})
default_rotation(::Val{3}, T) = one(QuatRotation{T})

"""
    radii(ball)

Return the radii of the metric `ball`.
"""
radii(ball::MetricBall) = ball.radii

"""
    rotation(ball)

Return the rotation of the metric `ball`.
"""
rotation(ball::MetricBall) = isnothing(ball.rotation) ? I : ball.rotation

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
  ball.metric isa Mahalanobis ? one(r) : r
end

"""
    isisotropic(ball)

Tells whether or not the metric `ball` is isotropic,
i.e. if all its radii are equal.
"""
isisotropic(ball::MetricBall) = length(unique(ball.radii)) == 1

function *(α::Real, ball::MetricBall)
  if ball.metric isa Mahalanobis
    MetricBall(α .* ball.radii, ball.rotation)
  else
    MetricBall(α .* ball.radii, nothing, ball.metric)
  end
end

function Base.show(io::IO, ball::MetricBall)
  n = length(ball.radii)
  r = n > 1 ? Tuple(ball.radii) : first(ball.radii)
  m = nameof(typeof(ball.metric))
  print(io, "MetricBall($r, $m)")
end
