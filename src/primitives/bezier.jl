# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BezierCurve(points)

A Bézier curve with control points `points`.
See https://en.wikipedia.org/wiki/Bézier_curve.
A point on the curve `b` can be evaluated by calling
`b(t)` with `t` between `0` and `1`.
The evaluation method defaults to DeCasteljau's algorithm
for accurate evaluation. Horner's method, faster with a
large number of points but less precise, can be used via
`b(t, Horner())`.

## Example

```julia
BezierCurve(Point2[(0.,0.),(1.,-1.)])
```
"""
struct BezierCurve{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Primitive{Dim,T}
  controls::V
end

BezierCurve(points::Vararg) = BezierCurve(collect(points))
BezierCurve(points::AbstractVector{<:Tuple}) = BezierCurve(Point.(points))

ncontrols(b::BezierCurve) = length(b.controls)

"""
Evaluation method used to obtain a point along
a Bézier curve from a parametric expression.
"""
abstract type BezierEvalMethod end

"""
Accurate evaluation using De Casteljau's algorithm.
Recommended for a small number of control points.

See https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm.
"""
struct DeCasteljau <: BezierEvalMethod end

"""
Approximate evaluation using Horner's method.
Recommended for a large number of control points,
if you can afford a precision loss.

See https://en.wikipedia.org/wiki/Horner%27s_method.
"""
struct Horner <: BezierEvalMethod end

(curve::BezierCurve)(t) = curve(t, DeCasteljau())

function (curve::BezierCurve)(t, ::DeCasteljau)
  if t < 0 || t > 1
    throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
  end
  ss = collect(segments(Chain(curve.controls)))
  points = [s(t) for s in ss]
  if length(points) == 1
    points[1]
  else
    BezierCurve(points)(t)
  end
end

function (curve::BezierCurve{<:Any,T})(t, ::Horner) where {T}
  if t < 0 || t > 1
    throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
  end
  tᶜ = one(T) - t
  tᵢ = one(T)
  cᵢ = one(T)
  cs = curve.controls
  bᵢ = coordinates(first(cs)) .* tᶜ
  n = ncontrols(curve) - 1
  for i in 1:n-1
    p = coordinates(cs[i+1])
    tᵢ = tᵢ * t
    cᵢ = (n - i + 1) * cᵢ / i
    bᵢ = (bᵢ + cᵢ * p * tᵢ) * tᶜ
  end
  Point(bᵢ + coordinates(last(cs)) * tᵢ * t)
end
