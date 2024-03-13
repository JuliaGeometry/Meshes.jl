# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# flip arguments so that points always come first
evaluate(d::PreMetric, g::Geometry, p::Point) = evaluate(d, p, g)

"""
    evaluate(Euclidean(), line1, line2)

Evaluate the minimum Euclidean distance between `line1` and `line2`.
"""
function evaluate(::Euclidean, line1::Line, line2::Line)
  # Inspired by: John Alexiou, Find shortest distance between lines in 3D
  #   URL (version 2022-11-04) https://math.stackexchange.com/q/2217845

  e1 = line1(1) - line1(0)
  e2 = line2(1) - line2(0)
  n = e1 × e2
  nn = n ⋅ n
  r = line2(0) - line1(0)

  # Find the Point's on each line where they are closest
  t1 = ((e2 × n) ⋅ r) / nn
  t2 = ((e1 × n) ⋅ r) / nn
  p1 = line1(t1)
  p2 = line2(t2)

  if iszero(nn)
    # lines are parallel, so pick an arbitrary Point on line1 and
    #   find distance from it to line2
    return evaluate(Euclidean(), line1(0), line2)
  else
    # find distance between closest Point's
    return evaluate(Euclidean(), p1, p2)
  end
end

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
    evaluate(::PreMetric, point1, point2)

Evaluate pre-metric between coordinates of `point1` and `point2`.
"""
evaluate(d::PreMetric, p1::Point, p2::Point) = evaluate(d, coordinates(p1), coordinates(p2))
