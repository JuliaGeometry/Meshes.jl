# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    pairwiseintersect(segments; [digits])

Compute pairwise intersections between n `segments`
with `digits` precision in O(n⋅log(n+k)) time
where k is the number of intersections, using
a sweep line algorithm. Similar to an optimal
Bentley-Ottmann algorithm in sparse systems,
and closer to O(n²) in dense systems.

Return intersection points and corresponding indices
of segments involved in the intersection.

By default, set `digits` based on the absolute
tolerance of the length type of the segments.

## Examples

```julia
points, seginds = pairwiseintersect(segments)
points[i] # i-th intersection point
seginds[i] # corresponding segments
```

## References

* Bentley & Ottmann 1979. [Algorithms for reporting and counting
  geometric intersections](https://ieeexplore.ieee.org/document/1675432)
"""
function pairwiseintersect(segments; digits=_digits(segments))
  # orient segments and round coordinates
  segs = map(segments) do s
    a, b = coordround.(extrema(s), digits=digits)
    a > b ? Segment(b, a) : Segment(a, b)
  end

  # extract first (or "x") coordinate from
  # first and last vertices of segments
  x(p) = flat(coords(p)).x
  xₛ = [x(first(vertices(s))) for s in segs]
  xₑ = [x(last(vertices(s))) for s in segs]

  # sort segments based on first coordinates
  inds = sortperm(xₛ)
  segs = segs[inds]
  xₛ = xₛ[inds]
  xₑ = xₑ[inds]

  # sweepline algorithm
  n = length(segs)
  P = eltype(first(segs))
  D = Dict{P,Vector{Int}}()
  for i in 1:n
    for j in (i + 1):n
      # break if segments don't overlap w.r.t. first coordinate
      xₛ[i] ≤ xₛ[j] ≤ xₑ[i] || break

      # perform more expensive intersection algorithm
      intersection(segs[i], segs[j]) do I
        if type(I) == Crossing || type(I) == EdgeTouching
          p = coordround(get(I); digits)
          if haskey(D, p)
            append!(D[p], (inds[i], inds[j]))
          else
            D[p] = [inds[i], inds[j]]
          end
        end
      end
    end
  end

  # remove duplicate indices
  for p in keys(D)
    unique!(D[p])
  end

  collect(keys(D)), collect(values(D))
end

# compute the number of significant digits based on the segment type
# this is used to determine the precision of the points
function _digits(segments)
  s = first(segments)
  ℒ = lentype(s)
  τ = ustrip(eps(ℒ))
  round(Int, 0.8 * (-log10(τ))) # 0.8 is a heuristic to avoid numerical issues
end
