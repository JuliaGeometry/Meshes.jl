# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BezierCurve(points)

A recursive Bézier curve with control points `points`.
See <https://en.wikipedia.org/wiki/Bézier_curve>.
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
struct BezierCurve{Dim,V<:AbstractVector{<:Point{Dim}}} <: Primitive{Dim}
  controls::V
end

BezierCurve(points::AbstractVector{<:Tuple}) = BezierCurve(Point.(points))
BezierCurve(points...) = BezierCurve(collect(points))

paramdim(::Type{<:BezierCurve}) = 1

lentype(::Type{<:BezierCurve{Dim,V}}) where {Dim,V} = lentype(eltype(V))

controls(b::BezierCurve) = b.controls

ncontrols(b::BezierCurve) = length(b.controls)

degree(b::BezierCurve) = ncontrols(b) - 1

"""
Evaluation method used to obtain a point along
a Bézier curve from a parametric expression.
"""
abstract type BezierEvalMethod end

"""
Accurate evaluation using De Casteljau's algorithm.
Recommended for a small number of control points.

See <https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm>.
"""
struct DeCasteljau <: BezierEvalMethod end

"""
Approximate evaluation using Horner's method.
Recommended for a large number of control points,
if you can afford a precision loss.

See <https://en.wikipedia.org/wiki/Horner%27s_method>.
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
function (curve::BezierCurve)(t, ::Horner)
  T = numtype(lentype(curve))
  if t < 0 || t > 1
    throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
  end
  cs = curve.controls
  t̄ = one(T) - t
  n = degree(curve)
  pₙ = to(last(cs))
  aₙ = pₙ

  # initialization with i = n + 1, so bᵢ₋₁ = bₙ = aₙ
  bᵢ₋₁ = aₙ
  cᵢ₋₁ = one(T)
  t̄ⁿ⁻ⁱ = one(T)
  for i in n:-1:1
    cᵢ₋₁ *= i / (n - i + one(T))
    pᵢ₋₁ = to(cs[i])
    t̄ⁿ⁻ⁱ *= t̄
    aᵢ₋₁ = cᵢ₋₁ * pᵢ₋₁ * t̄ⁿ⁻ⁱ
    bᵢ₋₁ = aᵢ₋₁ + bᵢ₋₁ * t
  end

  b₀ = bᵢ₋₁
  Point(coordinates(b₀))
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{BezierCurve{Dim}}) where {Dim} =
  BezierCurve([rand(rng, Point{Dim}) for _ in 1:5])

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, b::BezierCurve)
  ioctx = IOContext(io, :compact => true)
  print(io, "BezierCurve(controls: [")
  join(ioctx, b.controls, ", ")
  print(io, "])")
end

function Base.show(io::IO, ::MIME"text/plain", b::BezierCurve)
  summary(io, b)
  println(io)
  print(io, "└─ controls: [")
  join(io, b.controls, ", ")
  print(io, "]")
end
