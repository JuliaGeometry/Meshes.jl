# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    JarvisMarch()

Compute the convex hull of a set of points or geometries using the
Jarvis's march algorithm. See [https://en.wikipedia.org/wiki/Gift_wrapping_algorithm]
(https://en.wikipedia.org/wiki/Gift_wrapping_algorithm).

The algorithm has complexity `O(n*h)` where `n` is the number of points
and `h` is the number of points in the hull.

## References

* Jarvis 1973. [On the identification of the convex hull of a finite set of
  points in the plane](https://www.sciencedirect.com/science/article/abs/pii/0020019073900203)
"""
struct JarvisMarch <: HullMethod end

function hull(points, ::JarvisMarch)
  pₒ = first(points)
  Dim = embeddim(pₒ)
  ℒ = lentype(pₒ)

  assertion(Dim == 2, "Jarvis's march only defined in 2D")

  # remove duplicates
  p = unique(points)
  n = length(p)

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])

  # find bottom-left point
  i = argmin(l -> cartvalues(p[l]), 1:n)

  # candidates for next point
  𝒞 = [1:(i - 1); (i + 1):n]

  # find next point with smallest angle
  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  j = argmin(l -> ∠(A, O, p[l]), 𝒞)

  # initialize ring of indices
  ℐ = [i, j]

  # rotational sweep
  while first(ℐ) != last(ℐ)
    # direction of current segment
    v = p[j] - p[i]

    # update candidates
    𝒞 = setdiff(1:n, [i, j])

    # find next segment
    i = j
    O = p[i]
    A = O + v
    j = argmin(l -> ∠(A, O, p[l]), 𝒞)

    # update ring of indices
    push!(ℐ, j)
  end

  PolyArea(p[ℐ[begin:(end - 1)]])
end
