# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# Comparisons between unitful and non-unitful quantities

isequalzero(x) = x == zero(x)
isequalone(x) = x == oneunit(x)

isapproxequal(x, y; atol=atol(x), kwargs...) = isapprox(x, y; atol, kwargs...)
isapproxzero(x; atol=atol(x), kwargs...) = isapprox(x, zero(x); atol, kwargs...)
isapproxone(x; atol=atol(x), kwargs...) = isapprox(x, oneunit(x); atol, kwargs...)

ispositive(x) = x > zero(x)
isnegative(x) = x < zero(x)
isnonpositive(x) = x ≤ zero(x)
isnonnegative(x) = x ≥ zero(x)

"""
    mayberound(λ, x, tol)

Round `λ` to `x` if it is within the tolerance `tol`.
"""
function mayberound(λ::T, x::T, atol=atol(T)) where {T}
  isapprox(λ, x, atol=atol) ? x : λ
end

"""
    approxunique(points)

Return a collection of points that are unique up to the default tolerance for the type of the points.
"""
function approxunique(points)
  # spatial index to avoid comparing all pairs of points
  searcher = BallSearch(points, MetricBall(atol(lentype(first(points)))))

  keep = trues(length(points))
  for i in eachindex(points)
    keep[i] || continue
    for j in search(points[i], searcher)
      j > i && points[i] ≈ points[j] && (keep[j] = false)
    end
  end

  points[keep]
end
