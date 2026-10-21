# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Dilate(δ=0)

Morphological dilation with radius `δ`.
"""
struct Dilate{ℒ<:Len} <: GeometricTransform
  δ::ℒ
  function Dilate(δ::ℒ) where {ℒ<:Len}
    ispositive(δ) || throw(ArgumentError("dilation radius must be positive"))
    new{float(ℒ)}(δ)
  end
end

Dilate(δ) = Dilate(aslen(δ))

Dilate() = Dilate(1.0u"m")

parameters(t::Dilate) = (; δ=t.δ)

function apply(transform::Dilate, point::Point)
  δ = transform.δ
  b = Ball(point, δ)
  aspoly(b), nothing
end

function apply(transform::Dilate, chain::Chain{𝔼{2}})
  δ = transform.δ
  balls = mapreduce(∪, vertices(chain)) do p
    Ball(p, δ) |> aspoly
  end
  quads = mapreduce(∪, segments(chain)) do s
    a, b = extrema(s)
    v⃗ = b - a
    n⃗ = Vec(-v⃗[2], v⃗[1])
    δ⃗ = (δ / norm(n⃗)) * n⃗
    Quadrangle(a - δ⃗, b - δ⃗, b + δ⃗, a + δ⃗)
  end
  balls ∪ quads, nothing
end

function apply(transform::Dilate, poly::Polygon)
  bpoly = transform(boundary(poly))
  poly ∪ bpoly, nothing
end

function apply(transform::Dilate, multi::Multi)
  polys = mapreduce(transform, ∪, parent(multi))
  polys, nothing
end
