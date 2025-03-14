# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DouglasPeuckerSimplification(τ)

Douglas-Peucker's simplification algorithm with tolerance `τ` in length units
(default to meter).

The higher is the tolerance, the more aggressive is the simplification.

## References

* Douglas, D. and Peucker, T. 1973. [Algorithms for the Reduction of
  the Number of Points Required to Represent a Digitized Line or its
  Caricature](https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
struct DouglasPeuckerSimplification{ℒ<:Len} <: SimplificationMethod
  τ::ℒ
  DouglasPeuckerSimplification(τ::ℒ) where {ℒ<:Len} = new{float(ℒ)}(τ)
end

DouglasPeuckerSimplification(τ) = DouglasPeuckerSimplification(addunit(τ, u"m"))

function simplify(chain::Chain, method::DouglasPeuckerSimplification)
  verts = _douglaspeucker(vertices(chain), method.τ)
  isclosed(chain) ? Ring(verts) : Rope(verts)
end

# simplify chain assuming it is open
function _douglaspeucker(v::AbstractVector{P}, τ) where {P<:Point}
  # find vertex with maximum distance to reference line
  l = Line(first(v), last(v))
  imax, dmax = 0, zero(lentype(P))
  for i in 2:(length(v) - 1)
    d = evaluate(Euclidean(), v[i], l)
    if d > dmax
      imax = i
      dmax = d
    end
  end

  if dmax < τ
    [first(v), last(v)]
  else
    v₁ = _douglaspeucker(v[begin:imax], τ)
    v₂ = _douglaspeucker(v[imax:end], τ)
    [v₁[begin:(end - 1)]; v₂]
  end
end
