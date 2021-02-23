# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NormBall(radius, [metric])

A norm ball with `radius` and `metric`.
Default metric is `Euclidean()`.
"""
struct NormBall{T,M} <: MetricBall
  radius::T
  metric::M
end

function NormBall(radius)
  metric = Euclidean()
  NormBall{typeof(radius),typeof(metric)}(radius, metric)
end

metric(ball::NormBall) = ball.metric

"""
    radius(ball)

Return the radius of the norm `ball`.
"""
radius(ball::NormBall) = ball.radius

function Base.show(io::IO, ball::NormBall{T}) where {T}
  r, d = ball.radius, ball.metric
  print(io, "NormBall{$T}($r, $d)")
end
