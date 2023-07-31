# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BezierCurve(points)

A recursive Bézier curve with control points `points`.
See https://en.wikipedia.org/wiki/Bézier_curve.
A point on the curve `b` can be evaluated by calling
`b(t)` with `t` between `0` and `1`.
The evaluation method defaults to DeCasteljau's algorithm
for accurate evaluation. Horner's method, faster with a
large number of points but less precise, can be used via
`b(t, Horner())`.

## Examples

```julia
BezierCurve(Point2[(0.,0.),(1.,-1.)])
```
"""
struct BezierCurve{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Primitive{Dim,T}
  controls::V
end

BezierCurve(points::AbstractVector{<:Tuple}) = BezierCurve(Point.(points))
BezierCurve(points::Vararg) = BezierCurve(collect(points))

paramdim(::Type{<:BezierCurve}) = 1

controls(b::BezierCurve) = b.controls

ncontrols(b::BezierCurve) = length(b.controls)

degree(b::BezierCurve) = ncontrols(b) - 1

function boundary(b::BezierCurve)
  ps = b.controls
  p₁, p₂ = first(ps), last(ps)
  p₁ ≈ p₂ ? nothing : PointSet([p₁, p₂])
end

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

# Apply DeCasteljau's method
function (curve::BezierCurve)(t, ::DeCasteljau)
  if t < 0 || t > 1
    throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
  end
  ss = segments(Rope(curve.controls))
  points = [s(t) for s in ss]
  if length(points) == 1
    points[1]
  else
    BezierCurve(points)(t)
  end
end

# Apply Horner's method on the monomial representation of the
# Bézier curve B = ∑ᵢ aᵢtⁱ with i ∈ [0, n], n the degree of the
# curve, aᵢ = binomial(n, i) * pᵢ * t̄ⁿ⁻ⁱ and t̄ = (1 - t).
# Horner's rule recursively reconstructs B from a sequence bᵢ
# with bₙ = aₙ and bᵢ₋₁ = aᵢ₋₁ + bᵢ * t until b₀ = B.
function (curve::BezierCurve{Dim,T})(t, ::Horner) where {Dim,T}
  if t < 0 || t > 1
    throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
  end
  cs = curve.controls
  t̄ = one(T) - t
  n = degree(curve)
  pₙ = coordinates(last(cs))
  aₙ = pₙ

  # initialization with i = n + 1, so bᵢ₋₁ = bₙ = aₙ
  bᵢ₋₁ = aₙ
  cᵢ₋₁ = one(T)
  t̄ⁿ⁻ⁱ = one(T)
  for i in n:-1:1
    cᵢ₋₁ *= i / (n - i + one(T))
    pᵢ₋₁ = coordinates(cs[i])
    t̄ⁿ⁻ⁱ *= t̄
    aᵢ₋₁ = cᵢ₋₁ * pᵢ₋₁ * t̄ⁿ⁻ⁱ
    bᵢ₋₁ = aᵢ₋₁ + bᵢ₋₁ * t
  end

  b₀ = bᵢ₋₁
  Point(b₀)
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{BezierCurve{Dim,T}}) where {Dim,T} =
  BezierCurve(rand(rng, Point{Dim,T}, 5))
