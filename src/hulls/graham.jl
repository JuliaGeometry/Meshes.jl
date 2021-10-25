# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GrahamScan

Compute the convex hull of a set of points or geometries using the
Graham's scan method. See [https://en.wikipedia.org/wiki/Graham_scan]
(https://www.pnas.org/content/117/52/33711).

The method has complexity `O(n*log(n))` where `n` is the number of points.

## References

* Cormen et al. 2009. [Introduction to Algorithms]
  (https://mitpress.mit.edu/books/introduction-algorithms-third-edition)
"""
struct GrahamScan <: HullMethod end

function hull(pset::PointSet{2,T}, ::GrahamScan) where {T}
  # remove duplicates
  Q = coordinates.(pset) |> unique

  # sort by y then by x
  sort!(Q, by=reverse)

  # corner cases
  n = length(Q)
  n == 1 && return Point(Q[1])
  n == 2 && return Segment(Point(Q[1]), Point(Q[2]))
  if n == 3
    p₀, p₁, p₂ = Q
    θ = ∠(Point(p₁), Point(p₀), Point(p₂))
    if isapprox(θ, zero(T), atol=atol(T))
      return Segment(Point(p₀), Point(p₂))
    else
      c = Chain(Point(p₀), Point(p₁), Point(p₂), Point(p₀))
      c = orientation(c) == :CCW ? c : reverse(c)
      return PolyArea(c)
    end
  end

  # sort by polar angle
  p₀ = Point(Q[1])
  p  = Point.(Q[2:n])
  x  = p₀ + Vec{2,T}(1, 0)
  θ  = [∠(x, p₀, pᵢ) for pᵢ in p]
  p  = p[sortperm(θ)]

  # rotational sweep
  c = [p₀, p[1], p[2]]
  for pᵢ in p[3:end]
    while ∠(c[end-1], c[end], pᵢ) > atol(T)
      pop!(c)
    end
    push!(c, pᵢ)
  end

  # close chain
  push!(c, c[1])

  PolyArea(c)
end