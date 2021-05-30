# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DouglasPeucker(ϵ)

The Douglas-Peucker algorithm for simplifying polygonal chains or
polygonal areas given a deviation tolerance `ϵ`.

## References

* Douglas, D. and Peucker, T. 1973. [Algorithms for the Reduction of
  the Number of Points Required to Represent a Digitized Line or its
  Caricature](https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
struct DouglasPeucker{T} <: SimplificationMethod
  ϵ::T
end

function simplify(chain::Chain, method::DouglasPeucker)
  v = simplify(vertices(chain), method)
  isclosed(chain) ? Chain([v; first(v)]) : Chain(v)
end

# simplify chain assuming it is open
function simplify(v::AbstractVector{Point{Dim,T}},
                  method::DouglasPeucker) where {Dim,T}
  # find vertex with maximum distance
  imax, dmax = 0, zero(T)
  for i in 2:length(v)-1
    d = evaluate(Euclidean(), v[i], Line(first(v), last(v)))
    if d > dmax
      imax = i
      dmax = d
    end
  end

  if dmax < method.ϵ
    [first(v), last(v)]
  else
    v₁ = simplify(v[begin:imax], method)
    v₂ = simplify(v[imax:end],   method)
    [v₁[begin:end-1]; v₂]
  end
end
