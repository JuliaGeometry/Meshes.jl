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
struct MetricBall{Dim,ℒ<:Len,R,M} <: Neighborhood
  radii::NTuple{Dim,ℒ}
  rotation::R

  # state fields
  metric::M

  MetricBall(radii::NTuple{Dim,ℒ}, rotation::R, metric::M) where {Dim,ℒ<:Len,R,M} =
    new{Dim,float(ℒ),R,M}(radii, rotation, metric)
end

function MetricBall(radii::NTuple{Dim,ℒ}, rotation=nothing) where {Dim,ℒ<:Len}
  # scaling matrix
  Λ = Diagonal(SVector((oneunit(ℒ) ./ radii) .^ 2))

  # rotation matrix
  R = isnothing(rotation) ? default_rotation(Val(Dim), float(numtype(ℒ))) : rotation

  # anisotropy matrix
  M = Symmetric(R * Λ * R')

  # Mahalanobis metric
  metric = Mahalanobis(M)

  MetricBall(radii, R, metric)
end

MetricBall(radii::NTuple{Dim,Len}, rotation=nothing) where {Dim} = MetricBall(promote(radii...), rotation)

MetricBall(radii::Tuple, rotation=nothing) = MetricBall(addunit.(radii, u"m"), rotation)

MetricBall(radius::Len, metric=Euclidean()) = MetricBall((radius,), nothing, metric)

MetricBall(radius::Number, metric=Euclidean()) = MetricBall(addunit(radius, u"m"), metric)

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
  ball.metric isa Mahalanobis ? oneunit(r) : r
end

"""
    isisotropic(ball)

Tells whether or not the metric `ball` is isotropic,
i.e. if all its radii are equal.
"""
isisotropic(ball::MetricBall) = allequal(ball.radii)

function *(α::Real, ball::MetricBall)
  if ball.metric isa Mahalanobis
    MetricBall(α .* ball.radii, ball.rotation)
  else
    MetricBall(α .* ball.radii, nothing, ball.metric)
  end
end

function Base.show(io::IO, ball::MetricBall)
  n = length(ball.radii)
  r = n > 1 ? ball.radii : first(ball.radii)
  m = nameof(typeof(ball.metric))
  print(io, "MetricBall($r, $m)")
end
