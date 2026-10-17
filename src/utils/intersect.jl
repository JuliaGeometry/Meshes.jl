# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersectparameters(a, b, c, d)

Compute the parameters `λ₁` and `λ₂` of the lines 
`a + λ₁ ⋅ v⃗₁`, with `v⃗₁ = b - a` and
`c + λ₂ ⋅ v⃗₂`, with `v⃗₂ = d - c` spanned by the input
points `a`, `b` resp. `c`, `d` such that to yield line
points with minimal distance or the intersection point
(if lines intersect).

Furthermore, the ranks `r` of the matrix of the linear
system `A ⋅ λ⃗ = y⃗`, with `A = [v⃗₁ -v⃗₂], y⃗ = c - a`
and the rank `rₐ` of the augmented matrix `[A y⃗]` are
calculated in order to identify the intersection type:

- Intersection: r == rₐ == 2
- Collinear: r == rₐ == 1
- No intersection: r ≠ rₐ
  - No intersection and parallel:  r == 1, rₐ == 2
  - No intersection, skew lines: r == 2, rₐ == 3
"""
function intersectparameters(a::Point, b::Point, c::Point, d::Point)
  # augmented linear system
  A = ustrip.([(b - a) (c - d) (c - a)])

  # normalize by maximum absolute coordinate
  A = A / maximum(abs, A)

  # numerical tolerance
  T = eltype(A)
  τ = atol(T)

  # check if a vector is non zero
  isnonzero(v) = !isapprox(v, zero(v), atol=τ)

  # calculate ranks by checking the zero rows of
  # the factor R in the QR matrix factorization
  _, R = qr(A)
  r = sum(isnonzero, eachrow(R[:, SVector(1, 2)]))
  rₐ = sum(isnonzero, eachrow(R))

  # calculate parameters of intersection
  if r ≥ 2
    λ = A[:, SVector(1, 2)] \ A[:, 3]
    λ₁, λ₂ = λ[1], λ[2]
  else # parallel or collinear
    λ₁, λ₂ = zero(T), zero(T)
  end

  λ₁, λ₂, r, rₐ
end

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
  P = eltype(vertices(first(segs)))
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

# compute the intersection of two longitude intervals, taking into account the antimeridian
function lonintersects(minlon₁::T, maxlon₁::T, minlon₂::T, maxlon₂::T) where {T}
  # -180° and 180° represent the same meridian
  πdeg = T(180u"°")

  # convert range to one or two non-decreasing ranges,
  # depending on whether it crosses the antimeridian
  ranges₁ = minlon₁ ≤ maxlon₁ ? ((minlon₁, maxlon₁),) : ((-πdeg, maxlon₁), (minlon₁, πdeg))
  ranges₂ = minlon₂ ≤ maxlon₂ ? ((minlon₂, maxlon₂),) : ((-πdeg, maxlon₂), (minlon₂, πdeg))

  # compute the intersection of the two sets of ranges
  ranges = Tuple{T,T}[]
  for (a₁, b₁) in ranges₁, (a₂, b₂) in ranges₂
    a, b = max(a₁, a₂), min(b₁, b₂)
    a ≤ b && push!(ranges, (a, b))
  end

  # check if the intersection contains the antimeridian
  _hasπdeg(ranges) = any(lon -> lon[1] == -πdeg || lon[2] == πdeg, ranges)
  if _hasπdeg(ranges₁) && _hasπdeg(ranges₂) && !_hasπdeg(ranges)
    push!(ranges, (πdeg, πdeg))
  end

  # ranges meeting at the antimeridian form a single wrapped interval
  if length(ranges) ≥ 2 && first(ranges)[1] == -πdeg && last(ranges)[2] == πdeg
    a, b = last(ranges)[1], first(ranges)[2]
    return a == πdeg && b == -πdeg ? [(πdeg, πdeg)] : [(a, b)]
  end

  ranges
end
