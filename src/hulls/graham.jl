# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GrahamScan()

Compute the convex hull of a set of points or geometries using the
Graham's scan algorithm. See [https://en.wikipedia.org/wiki/Graham_scan]
(https://en.wikipedia.org/wiki/Graham_scan).

The algorithm has complexity `O(n*log(n))` where `n` is the number of points.

## References

* Cormen et al. 2009. [Introduction to Algorithms]
  (https://mitpress.mit.edu/books/introduction-algorithms-third-edition)
"""
struct GrahamScan <: HullMethod end

function hull(points, ::GrahamScan)
  pₒ = first(points)
  Dim = embeddim(pₒ)
  T = coordtype(pₒ)

  @assert Dim == 2 "Graham's scan only defined in 2D"

  # remove duplicates
  p = unique(points)
  n = length(p)

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])

  # sort points lexicographically
  p = p[sortperm(coordinates.(p))]

  # sort points by polar angle
  O = p[1]
  q = p[2:n]
  A = O + Vec{2,T}(0, -1)
  θ = [∠(A, O, B) for B in q]
  q = q[sortperm(θ)]

  # initialise and skip collinear points at beginning 
  r = [O]
  idx = 1
  while idx ≤ length(q) && coordinates(q[idx])[2] == coordinates(O)[2]
    idx += 1
  end
  idx = max(idx, 2)
  push!(r, q[idx - 1])
  if idx > length(q)
    return PolyArea(r)
  else
    push!(r, q[idx])
  end

  # rotational sweep
  for B in q[(idx + 1):end]
    while ∠(r[end - 1], r[end], B) > atol(T) && length(r) ≥ 3
      pop!(r)
    end
    if !iscollinear(r[end - 1], r[end], B)
      push!(r, B)
    end
  end

  PolyArea(r)
end
