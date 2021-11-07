# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    IsotropicBall(radius, [metric])

A ball with single `radius` and `metric`.
Default metric is `Euclidean()`.
"""
struct IsotropicBall{T,M} <: MetricBall
  radius::T
  metric::M
end

function IsotropicBall(radius)
  metric = Euclidean()
  T = typeof(radius)
  M = typeof(metric)
  IsotropicBall{T,M}(radius, metric)
end

metric(ball::IsotropicBall) = ball.metric

radii(ball::IsotropicBall)  = SVector(ball.radius,)

function Base.show(io::IO, ball::IsotropicBall{T}) where {T}
  r, d = ball.radius, ball.metric
  print(io, "IsotropicBall{$T}($r, $d)")
end
