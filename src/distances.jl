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
  a, b = points(l)
  u = coordinates(p)
  v = b - a
  α = (u ⋅ v) / (v ⋅ v)
  norm(u - α*v)
end
