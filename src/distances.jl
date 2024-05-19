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
  a, b = l(0), l(1)
  u = p - a
  v = b - a
  α = (u ⋅ v) / (v ⋅ v)
  norm(u - α * v)
end

"""
    evaluate(Euclidean(), line1, line2)
Evaluate the minimum Euclidean distance between `line1` and `line2`.
"""
function evaluate(::Euclidean, line1::Line{Dim}, line2::Line{Dim}) where {Dim}
  λ₁, λ₂, r, rₐ = intersectparameters(line1(0), line1(1), line2(0), line2(1))

  if (r == rₐ == 2) || (r == rₐ == 1)  # lines intersect or are colinear
    return zero(lentype(line1))
  elseif (r == 1) && (rₐ == 2)  # lines are parallel
    return evaluate(Euclidean(), line1(0), line2)
  else  # get distance between closest points on each line
    return evaluate(Euclidean(), line1(λ₁), line2(λ₂))
  end
end

"""
    evaluate(::PreMetric, point1, point2)

Evaluate pre-metric between coordinates of `point1` and `point2`.
"""
evaluate(d::PreMetric, p1::Point, p2::Point) = evaluate(d, coordinates(p1), coordinates(p2))
