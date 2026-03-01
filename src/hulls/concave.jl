# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Concave()

Compute the concave hull of a set of points or geometries using a
knearest approach. See [https://en.wikipedia.org/wiki/Gift_wrapping_algorithm]
(https://en.wikipedia.org/wiki/Gift_wrapping_algorithm).

The algorithm has complexity `O(n)` where `n` is the number of points

## References

"""
struct Concave <: HullMethod end

function hull(points, ::Concave; k=3)
  kk = max(k, 3)
  println("Computing concave hull with k=$kk")
  pₒ = first(points)
  ℒ = lentype(pₒ)
  T = numtype(ℒ)

  # sanity check
  ncoords = CoordRefSystems.ncoords(coords(pₒ))
  assertion(ncoords == 2, "This concave hull algorithm is only defined with 2D coordinates")

  # remove duplicates
  p = unique(points)
  n = length(p)

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])
  n == 3 && return PolyArea(p)

  kk = min(kk, n - 1)

  # prevent infinite recursion - if k is too large, fall back to convex hull
  kk >= n - 1 && return convexhull(p)

  # find bottom-left point
  i = argmin(p)
  searcher = KNearestSearch(p, kk)

  # candidates for next point
  𝒞 = [1:(i - 1); (i + 1):n]

  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  j = argmin(l -> ∠(A, O, p[l]), 𝒞)

  # initialize ring of indices
  ℐ = [i, j]

  mask = trues(length(points))
  mask[[i, j]] .= false

  # rotational sweep
  step = 2
  while first(ℐ) != last(ℐ)
    step == 5 && (mask[ℐ[begin]] = 1)
    # direction of current segment
    curr = p[j]
    last = p[i]
    v = curr - last

    # update candidates
    𝒞 = setdiff(1:n, [i, j])

    # find next segment
    i = j
    O = p[i]
    A = O + v

    neighbors = search(curr, searcher; mask=mask)
    k = min(kk, length(neighbors))
    𝒩 = neighbors[1:k]

    sort!(𝒩, by=l -> ∠(A, O, p[l]))

    its = true
    indᵢ = 0
    while its && indᵢ < k
      indᵢ += 1
      cpointᵢ = p[𝒩[indᵢ]]
      lastpoint = cpointᵢ == p[ℐ[begin]] ? 1 : 0
      indⱼ = 2
      its = false
      while !its && indⱼ < length(ℐ) - lastpoint
        its = intersects(Segment(p[ℐ[step]], cpointᵢ), Segment(p[ℐ[step - indⱼ + 1]], p[ℐ[step - indⱼ]]))
        indⱼ += 1
      end
    end

    its && return hull(points, Concave(); k=kk + 1)

    j = 𝒩[indᵢ]
    mask[𝒩[indᵢ]] = false
    step += 1

    push!(ℐ, 𝒩[indᵢ])
  end
  poly = PolyArea(p[ℐ[begin:(end - 1)]])
  # if not all points are in the polygon, increase k and try again
  !all(points .∈ poly) && return hull(points, Concave(); k=kk + 1)
  poly
end
