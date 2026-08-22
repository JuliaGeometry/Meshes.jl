# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MoreiraMarch(k)

Compute the concave hull of a set of points or geometries using Moreira's
march algorithm. `k` is the minimum number of nearest neighbors in the
interval `(2, n)` where `n` is the number of unique points in the input.

The algorithm increases `k` until a valid hull is found, or until `k` reaches
`n - 1`. It has complexity `O(k*n*h)` where `h` is the number of points in the
hull, per attempt. The number of attempts is bounded above by `n - 1 - k`.

## References

* Moreira, A. & Santos, M. Y. 2007. [Concave hull: a k-nearest-neighbours
  approach for the computation of the region occupied by a set of points]
  (https://www.semanticscholar.org/paper/Concave-hull:-A-k-nearest-neighbours-approach-for-a-Moreira-Santos/319a3450f9909043d46eb7ceb4299efceb984d4f)
"""
struct MoreiraMarch <: HullMethod
  k::Int

  function MoreiraMarch(k)
    assertion(isinteger(k) && k > 2, "k must be greater than 2")
    new(k)
  end
end

function hull(points, method::MoreiraMarch)
  # sanity check
  ncoords = CoordRefSystems.ncoords(coords(first(points)))
  assertion(ncoords == 2, "Moreira's march algorithm is only defined with 2D coordinates")

  # remove duplicates
  p = unique(points)
  n = length(p)

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])
  n == 3 && return PolyArea(p)

  # sanity check
  assertion(method.k < n, "k must be smaller than the number of unique points")

  # recursively attempt to find a valid hull, increasing k if necessary
  _moreiramarch(p, method.k)
end

function _moreiramarch(p, k)
  # auxiliary variables
  n = length(p)
  ℒ = lentype(first(p))

  # number of neighbors at each step
  m = max(k, 3)

  # k-nearest neighbor search
  searcher = KNearestSearch(p, m)

  # find bottom-left point
  i = argmin(p)

  # initialize ring of indices
  ℐ = [i]

  # initialize search mask
  mask = trues(n)
  mask[i] = false

  # pre-allocate vector of candidate indices
  𝒞 = Vector{Int}(undef, m)

  # find next point with smallest angle
  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  j = _moreiranext!(𝒞, searcher, mask, p, ℐ, A, O)

  # if no next point can be found, increase k and try again
  isnothing(j) && return _moreiramarch(p, k + 1)

  # add point to ring and remove from candidacy
  push!(ℐ, j)
  mask[j] = false

  # rotational sweep
  while first(ℐ) != last(ℐ)
    # start point re-enters candidacy once the ring has enough vertices to close
    length(ℐ) == 4 && (mask[first(ℐ)] = true)

    # direction of current segment
    v = p[j] - p[i]

    # find next point with smallest angle
    i = j
    O = p[i]
    A = O + v
    j = _moreiranext!(𝒞, searcher, mask, p, ℐ, A, O)

    # if no next point can be found, increase k and try again
    isnothing(j) && return _moreiramarch(p, k + 1)

    # add point to ring and remove from candidacy
    push!(ℐ, j)
    mask[j] = false
  end

  # construct polygonal area from ring of indices
  poly = PolyArea(p[ℐ[begin:(end - 1)]])

  # every point must be in the hull, otherwise increase k and try again
  all(∈(poly), p) || return _moreiramarch(p, k + 1)

  poly
end

# find the nearest candidate by angle whose segment avoids existing hull edges
function _moreiranext!(𝒞, searcher, mask, p, ℐ, A, O)
  n = search!(𝒞, O, searcher; mask)
  for i in sort!(view(𝒞, 1:n), by=i -> ∠(A, O, p[i]))
    segᵢ = Segment(p[last(ℐ)], p[i])
    # check if the segment intersects any existing edge
    # skip the last edge if the candidate is the start point
    jₛ = i == first(ℐ) ? 2 : 1
    jₑ = length(ℐ) - 2
    valid = !any(jₛ:jₑ) do j
      segⱼ = Segment(p[ℐ[j]], p[ℐ[j + 1]])
      intersects(segᵢ, segⱼ)
    end
    valid && return i
  end
  nothing
end
