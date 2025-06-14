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
  ℒ = lentype(pₒ)

  # sanity check
  assertion(embeddim(pₒ) == 2, "Graham's scan algorithm is only defined in 2D")

  # remove duplicates
  p = unique(points)
  n = length(p)

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])

  # find bottom-left point
  i = argmin(p)
  O = p[i]

  # sort other points by polar angle
  q = p[setdiff(1:n, i)]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  sort!(q, by=B -> ∠(A, O, B))

  # rotational sweep
  r = [O, q[1]]
  for B in q[2:end]
    Δ = signarea(r[end - 1], r[end], B)
    while isnegative(Δ) && length(r) > 2
      pop!(r)
      Δ = signarea(r[end - 1], r[end], B)
    end
    if ispositive(Δ)
      push!(r, B)
    elseif evaluate(Euclidean(), r[end - 1], r[end]) < evaluate(Euclidean(), r[end - 1], B)
      # point is collinear and further away, i.e., r[end - 1] --> r[end] --> B
      pop!(r)
      push!(r, B)
    end
  end

  # return polygonal area
  PolyArea(r)
end
