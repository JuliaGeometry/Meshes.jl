# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DouglasPeuckerSimplification(ϵ)

Simplify geometries with Douglas-Peucker algorithm. The higher
is the tolerance `ϵ`, the more aggressive is the simplification.

## References

* Douglas, D. and Peucker, T. 1973. [Algorithms for the Reduction of
  the Number of Points Required to Represent a Digitized Line or its
  Caricature](https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
struct DouglasPeuckerSimplification{ℒ<:Len} <: SimplificationMethod
  ϵ::ℒ
  DouglasPeuckerSimplification(ϵ::ℒ) where {ℒ<:Len} = new{float(ℒ)}(ϵ)
end

DouglasPeuckerSimplification(ϵ) = DouglasPeuckerSimplification(addunit(ϵ, u"m"))

function simplify(chain::Chain, method::DouglasPeuckerSimplification)
  verts = _douglaspeucker(vertices(chain), method.ϵ) |> collect
  isclosed(chain) ? Ring(verts) : Rope(verts)
end

# simplify chain assuming it is open
function _douglaspeucker(v::AbstractVector{P}, ϵ) where {P<:Point}
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
    v₁ = _douglaspeucker(v[begin:imax], ϵ)
    v₂ = _douglaspeucker(v[imax:end], ϵ)
    [v₁[begin:(end - 1)]; v₂]
  end
end
