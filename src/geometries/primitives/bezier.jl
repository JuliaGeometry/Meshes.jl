# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BezierCurve(points)

A recursive Bézier curve with control points `points`.
See <https://en.wikipedia.org/wiki/Bézier_curve>.

A point on the Bézier curve `b` can be evaluated by calling
`b(t)` with `t` between `0` and `1`. The evaluation method
defaults to DeCasteljau's algorithm for accurate evaluation.
Horner's method, faster with a large number of points but
less precise, can be used via `b(t, Horner())`.

## Examples

```julia
BezierCurve([(0.,0.),(1.,-1.)])
```
"""
struct BezierCurve{M<:Manifold,C<:CRS,V<:AbstractVector{Point{M,C}}} <: Primitive{M,C}
  controls::V
end

BezierCurve(points::AbstractVector{<:Tuple}) = BezierCurve(Point.(points))
BezierCurve(points...) = BezierCurve(collect(points))

paramdim(::Type{<:BezierCurve}) = 1

controls(b::BezierCurve) = b.controls

ncontrols(b::BezierCurve) = length(b.controls)

degree(b::BezierCurve) = ncontrols(b) - 1

==(b₁::BezierCurve, b₂::BezierCurve) = controls(b₁) == controls(b₂)

Base.isapprox(b₁::BezierCurve, b₂::BezierCurve; atol=atol(lentype(b₁)), kwargs...) =
  ncontrols(b₁) == ncontrols(b₂) &&
  all(isapprox(p₁, p₂; atol, kwargs...) for (p₁, p₂) in zip(controls(b₁), controls(b₂)))

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

(b::BezierCurve)(t) = b(t, DeCasteljau())

# Apply DeCasteljau's method
function (b::BezierCurve)(t, ::DeCasteljau)
  if t < 0 || t > 1
    throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
  end
  rope = Rope(controls(b))
  points = [s(t) for s in segments(rope)]
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
function (b::BezierCurve)(t, ::Horner)
  if t < 0 || t > 1
    throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
  end
  T = numtype(lentype(b))
  cs = controls(b)
  t̄ = one(T) - t
  n = degree(b)
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
  withcrs(b, b₀)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, b::BezierCurve)
  ioctx = IOContext(io, :compact => true)
  print(io, "BezierCurve(controls: [")
  join(ioctx, controls(b), ", ")
  print(io, "])")
end

function Base.show(io::IO, ::MIME"text/plain", b::BezierCurve)
  summary(io, b)
  println(io)
  print(io, "└─ controls: [")
  join(io, controls(b), ", ")
  print(io, "]")
end
