# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MoreiraSantosMarch(k)

Compute the concave hull of a set of points using k-nearest-neighbours. `k` must
be between 2 and `n - 1`, where `n` is the number of unique points in the input.

The algorithm has complexity `O(k*n*h)` where `n` is the number of points
and `h` is the number of points in the hull, per attempt, and the number of
attempts is bounded above by `n - 1 - k`.

## References

* Moreira, A. & Santos, M. Y. 2007. [Concave hull: a k-nearest-neighbours
  approach for the computation of the region occupied by a set of points]
  (https://www.semanticscholar.org/paper/Concave-hull:-A-k-nearest-neighbours-approach-for-a-Moreira-Santos/319a3450f9909043d46eb7ceb4299efceb984d4f)
"""
struct MoreiraSantosMarch <: HullMethod
  k::Int

  function MoreiraSantosMarch(k)
    assertion(isinteger(k) && k > 2, "k must be greater than 2")
    new(k)
  end
end

function hull(points, method::MoreiraSantosMarch)
  pₒ = first(points)

  # sanity check
  ncoords = CoordRefSystems.ncoords(coords(pₒ))
  assertion(ncoords == 2, "Moreira & Santos's march algorithm is only defined with 2D coordinates")

  # remove duplicates
  p = unique(points)
  n = length(p)

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])
  n == 3 && return PolyArea(p)

  # sanity check
  assertion(method.k < n, "k must be less than the number of unique points")

  # main worker
  _moreirahull(p, method.k)
end

function _moreirahull(p, k)
  n = length(p)
  ℒ = lentype(first(p))

  # clamp the number of neighbours considered at each step
  kk = clamp(k, 3, n - 1)

  # find bottom-left point as the start of the ring
  i = argmin(p)

  # initinitialize k nearest searcher
  searcher = KNearestSearch(p, kk + 1)

  # available points: everything but the start, until the ring has ≥ 4 vertices
  mask = trues(n)
  mask[i] = false

  # initialize ring of indices
  ℐ = [i]

  # find next point with smallest angle
  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  j = moreiranext(searcher, mask, p, ℐ, A, O)

  # if no next point can be found, increase k and try again
  isnothing(j) && return _moreirahull(p, k + 1)

  push!(ℐ, j)
  mask[j] = false

  # rotational sweep
  while first(ℐ) != last(ℐ)
    # start point re-enters candidacy once the ring has enough vertices to close
    length(ℐ) >= 4 && (mask[first(ℐ)] = true)

    v = p[j] - p[i]
    i = j
    O = p[i]
    A = O + v

    j = moreiranext(searcher, mask, p, ℐ, A, O)
    isnothing(j) && return _moreirahull(p, k + 1)

    push!(ℐ, j)
    mask[j] = false
  end

  ring = p[ℐ[begin:(end - 1)]]
  poly = PolyArea(ring)

  # every input point must fall inside or on the resulting hull
  all(q -> q ∈ poly, p) || return _moreirahull(p, k + 1)

  poly
end

# find the nearest candidate by angle whose segment avoids existing hull edges.
function moreiranext(searcher, mask, p, ℐ, A, O)
  𝒞 = search(O, searcher; mask)
  isempty(𝒞) && return nothing

  # iterate over candidates in order of increasing angle, returning the first valid one
  for nᵢ in sort(𝒞, by=l -> ∠(A, O, p[l]))
    cseg = Segment(p[ℐ[end]], p[nᵢ])
    cbox = boundingbox(cseg) # bbox check for efficiency
    tₒ = nᵢ == ℐ[begin] ? 2 : 1
    valid = !any(tₒ:(length(ℐ) - 2)) do t
      eseg = Segment(p[ℐ[t]], p[ℐ[t + 1]])
      intersects(cbox, boundingbox(eseg)) && intersects(cseg, eseg)
    end
    valid && return nᵢ
  end

  nothing
end
