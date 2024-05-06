# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DouglasPeucker([ϵ]; min=3, max=typemax(Int), maxiter=10)

Simplify geometries with Douglas-Peucker algorithm. The higher
is the tolerance `ϵ`, the more aggressive is the simplification.

If the tolerance `ϵ` is not provided, perform binary search until
the number of vertices is between `min` and `max` or until the
number of iterations reaches a maximum `maxiter`.

## References

* Douglas, D. and Peucker, T. 1973. [Algorithms for the Reduction of
  the Number of Points Required to Represent a Digitized Line or its
  Caricature](https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
struct DouglasPeucker{T} <: SimplificationMethod
  ϵ::T
  min::Int
  max::Int
  maxiter::Int
end

DouglasPeucker(ϵ=nothing; min=3, max=typemax(Int), maxiter=10) = DouglasPeucker(_ϵ(ϵ), min, max, maxiter)

_ϵ(ϵ::Nothing) = ϵ
_ϵ(ϵ::Number) = _ϵ(addunit(ϵ, u"m"))
_ϵ(ϵ::Len) = float(ϵ)

function simplify(chain::Chain, method::DouglasPeucker)
  v = if isnothing(method.ϵ)
    # perform binary search with other parameters
    βsimplify(vertices(chain), method.min, method.max, method.maxiter)
  else
    # perform Douglas-Peucker ϵ-simplification
    ϵsimplify(vertices(chain), method.ϵ)
  end |> collect
  isclosed(chain) ? Ring(v) : Rope(v)
end

# simplification by means of binary search
function βsimplify(v::AbstractVector{P}, min, max, maxiter) where {P<:Point}
  i = 0
  u = v
  n = length(u)
  a = zero(lentype(P))
  b = initeps(u)
  while !(min ≤ n ≤ max) && i < maxiter
    # midpoint candidate
    ϵ = (a + b) / 2

    # evaluate at midpoint
    u = ϵsimplify(v, ϵ)
    n = length(u)

    # binary search
    n < min && (b = ϵ)
    n > max && (a = ϵ)

    i += 1
  end

  u
end

# initial ϵ guess for a given chain
function initeps(v::AbstractVector{P}) where {P<:Point}
  n = length(v)
  ϵ = typemax(lentype(P))
  l = Line(first(v), last(v))
  d = [evaluate(Euclidean(), v[i], l) for i in 2:(n - 1)]
  ϵ = quantile(d, 0.25)
  2ϵ
end

# simplify chain assuming it is open
function ϵsimplify(v::AbstractVector{P}, ϵ) where {P<:Point}
  # find vertex with maximum distance
  # to reference line
  l = Line(first(v), last(v))
  imax, dmax = 0, zero(lentype(P))
  for i in 2:(length(v) - 1)
    d = evaluate(Euclidean(), v[i], l)
    if d > dmax
      imax = i
      dmax = d
    end
  end

  if dmax < ϵ
    [first(v), last(v)]
  else
    v₁ = ϵsimplify(v[begin:imax], ϵ)
    v₂ = ϵsimplify(v[imax:end], ϵ)
    [v₁[begin:(end - 1)]; v₂]
  end
end
