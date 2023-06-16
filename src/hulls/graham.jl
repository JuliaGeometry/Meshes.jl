# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GrahamScan

Compute the convex hull of a set of points or geometries using the
Graham's scan method. See [https://en.wikipedia.org/wiki/Graham_scan]
(https://en.wikipedia.org/wiki/Graham_scan).

The method has complexity `O(n*log(n))` where `n` is the number of points.

## References

* Cormen et al. 2009. [Introduction to Algorithms]
  (https://mitpress.mit.edu/books/introduction-algorithms-third-edition)
"""
struct GrahamScan <: HullMethod end

function hull(points::AbstractVector{Point{2,T}}, ::GrahamScan) where {T}
  # remove duplicates
  p = unique(points)
  n = length(p)

  # sort by y then by x
  p = p[sortperm(reverse.(coordinates.(p)))]

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])
  if n == 3
    if iscollinear(p...)
      return Segment(first(p), last(p))
    else
      return PolyArea(p)
    end
  end

  # sort by polar angle
  p₀ = p[1]
  q = p[2:n]
  x = p₀ + Vec{2,T}(1, 0)
  θ = [∠(x, p₀, pᵢ) for pᵢ in q]
  q = q[sortperm(θ)]

  # rotational sweep
  c = [p₀, q[1], q[2]]
  for pᵢ in q[3:end]
    while ∠(c[end - 1], c[end], pᵢ) > atol(T)
      pop!(c)
    end
    push!(c, pᵢ)
  end

  PolyArea(c)
end
