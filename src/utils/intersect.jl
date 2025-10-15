# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    pairwiseintersect(segments; [digits])

Compute pairwise intersections between n `segments`
with `digits` precision in O(n⋅log(n+k)) time using
a sweep line algorithm. Similar to an optimal
Bentley-Ottmann algorithm in sparse systems,
and closer to O(n²) in dense systems.

By default, set `digits` based on the absolute
tolerance of the length type of the segments.

## References

* Bentley & Ottmann 1979. [Algorithms for reporting and counting
  geometric intersections](https://ieeexplore.ieee.org/document/1675432)
"""
function pairwiseintersect(segments; digits=_digits(segments))
  # orient segments and round coordinates
  segs = map(segments) do seg
    a, b = coordround.(extrema(seg), digits=digits)
    a > b ? (b, a) : (a, b)
  end

  starts = [CoordRefSystems.values(coords(seg[1]))[1] for seg in segs]
  stops = [CoordRefSystems.values(coords(seg[2]))[1] for seg in segs]

  # sort segments based on start coordinates
  inds = sortperm(starts)
  starts = starts[inds]
  stops = stops[inds]
  segs = segs[inds]

  # sweepline algorithm
  P = eltype(first(segs))
  𝐺 = Dict{P,Vector{Int}}()
  for i in eachindex(segs)
    for j in (i + 1):length(segs)
      _overlaps(starts[i], stops[i], starts[j]) || break

      intersection(Segment(segs[i]), Segment(segs[j])) do I
        if type(I) == Crossing || type(I) == EdgeTouching
          _addintersection!(𝐺, get(I), inds[i], inds[j]; digits=digits)
        end
      end
      continue
    end
  end
  (collect(keys(𝐺)), collect(values(𝐺)))
end

_overlaps(startᵢ, stopᵢ, startₖ) = (startᵢ ≤ startₖ ≤ stopᵢ)

# compute the number of significant digits based on the segment type
# this is used to determine the precision of the points
function _digits(segments)
  seg = first(segments)
  ℒ = lentype(seg)
  τ = ustrip(eps(ℒ))
  round(Int, 0.8 * (-log10(τ))) # 0.8 is a heuristic to avoid numerical issues
end

# add an intersection point to the dictionary with segment indices
function _addintersection!(𝐺, I::Point, index₁::Int, index₂::Int; digits=10)
  p = coordround(I, digits=digits)
  if haskey(𝐺, p)
    append!(𝐺[p], (index₁, index₂))
    unique!(𝐺[p])
  else
    𝐺[p] = Vector{Int}([index₁, index₂])
  end
end
