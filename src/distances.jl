# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# flip arguments so that points always come first
evaluate(d::PreMetric, g::Geometry, p::Point) = evaluate(d, p, g)

"""
    evaluate(Euclidean(), point, line)

Evaluate the Euclidean distance between `point` and `line`.
"""
function evaluate(::Euclidean, p::Point, l::Line)
  u = coordinates(p)
  a, b = points(l)
  v = b - a
  α = (u ⋅ v) / (v ⋅ v)
  norm(u - α*v)
end

"""
    evaluate(::PreMetric, point1, point2)

Evaluate pre-metric between coordinates of `point2` and `point2`.
"""
evaluate(d::PreMetric, p1::Point, p2::Point) =
  evaluate(d, coordinates(p1), coordinates(p2))
