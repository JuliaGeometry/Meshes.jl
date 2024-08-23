# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MinMaxSimplification(method; min=3, max=typemax(Int), maxiter=10)

Simplify geometries with binary search algorithm and a parent simplification `method`.

The simplification is performed until the number of vertices is in the `[min, max]`
range or until a maximum number of iterations `maxiter` is reached.
"""
struct MinMaxSimplification{M} <: SimplificationMethod
  method::M
  min::Int
  max::Int
  maxiter::Int
end

MinMaxSimplification(method; min=3, max=typemax(Int), maxiter=10)

function simplify(c::Chain, m::MinMaxSimplification)
  i = 0
  s = c
  n = nvertices(c)
  a, b = _initrange(c)
  while !(m.min ≤ n ≤ m.max) && i < m.maxiter
    # midpoint candidate
    ϵ = (a + b) / 2

    # evaluate at midpoint
    s = simplify(c, m.method(ϵ))
    n = nvertices(s)

    # binary search
    n < m.min && (b = ϵ)
    n > m.max && (a = ϵ)

    i += 1
  end

  s
end

# initial range for binary search
function _initrange(c)
  v = vertices(c)
  n = length(v)
  l = Line(first(v), last(v))
  d = [evaluate(Euclidean(), v[i], l) for i in 2:(n - 1)]
  z = zero(lentype(c))
  ϵ = quantile(d, 0.25)
  (z, 2ϵ)
end
