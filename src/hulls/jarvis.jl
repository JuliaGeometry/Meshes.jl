# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    JarvisMarch()

Compute the convex hull of a set of points or geometries using the
Jarvis's march algorithm. See [https://en.wikipedia.org/wiki/Gift_wrapping_algorithm]
(https://en.wikipedia.org/wiki/Gift_wrapping_algorithm).

    JarvisMarch(k)

If `k` is provided, the algorithm will attempt to compute a concave hull using the
k nearest neighbors as proposed by Moreira & Santos 2007. The value of `k` must be
greater than 2 and less than the number of unique points.

The algorithm has complexity `O(n*h)` where `n` is the number of points
and `h` is the number of points in the hull. The concave variant adds a
`k`-nearest-neighbor search and hull-intersection checks per hull vertex.

## References

* Jarvis 1973. [On the identification of the convex hull of a finite set of
  points in the plane](https://www.sciencedirect.com/science/article/abs/pii/0020019073900203)
* Moreira, A. & Santos, M. Y. 2007. [concave hull: a k-nearest neighbours approach for the computation of the region occupied by a set of points](https://www.semanticscholar.org/paper/Concave-hull:-A-k-nearest-neighbours-approach-for-a-Moreira-Santos/319a3450f9909043d46eb7ceb4299efceb984d4f)
"""
struct JarvisMarch <: HullMethod
  k::Union{Nothing,Int}

  function JarvisMarch(k)
    isok = isnothing(k) || (isinteger(k) && k > 2)
    assertion(isok, "k must be greater than 2 or nothing")
    new(k)
  end
end

JarvisMarch() = JarvisMarch(nothing)

function hull(points, method::JarvisMarch)
  pₒ = first(points)
  ℒ = lentype(pₒ)
  k = method.k

  # sanity check
  ncoords = CoordRefSystems.ncoords(coords(pₒ))
  assertion(ncoords == 2, "Jarvis's march algorithm is only defined with 2D coordinates")

  # remove duplicates
  p = unique(points)
  n = length(p)

  # sanity check
  isnothing(k) || assertion(k < n, "k must be less than the number of unique points")

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])

  # initialize searcher and mask of unseen points
  searcher, unseen = jarvissearcher(p, k)

  # find bottom-left point
  i = argmin(p)

  # initialize hull
  ℐ = [i]

  # initialize candidates for next point
  𝒞 = jarviscandidates!(searcher, unseen, p, ℐ)

  # find next point with smallest angle
  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  j = jarvisnext(searcher, 𝒞, p, ℐ, A, O)

  # update hull and mark point as seen
  push!(ℐ, j)
  jarvisupdate!(searcher, unseen, j)

  # rotational sweep
  while first(ℐ) != last(ℐ)
    # direction of current segment
    v = p[j] - p[i]

    # update candidate points
    𝒞 = jarviscandidates!(searcher, unseen, p, ℐ)

    # find next point
    i = j
    O = p[i]
    A = O + v
    j = jarvisnext(searcher, 𝒞, p, ℐ, A, O)

    # no valid next point, should only happen if k is too small
    isnothing(j) && throw(ArgumentError("could not find concave hull with k = $k, try a larger k"))

    # update ring of indices
    push!(ℐ, j)
    jarvisupdate!(searcher, unseen, j)
  end

  PolyArea(p[ℐ[begin:(end - 1)]])
end

# helper to create searcher and mask of unseen points
jarvissearcher(p, k::Nothing) = nothing, nothing
jarvissearcher(p, k::Integer) = KNearestSearch(p, k), trues(length(p))

# helper to get candidate indices for next point
jarviscandidates!(searcher::Nothing, unseen, p, ℐ) = setdiff(1:length(p), last(ℐ, 2))
function jarviscandidates!(searcher::KNearestSearch, unseen, p, ℐ)
  # exclude the current point and its predecessor without allocating a new
  # mask: flip their bits off in place, search, then restore them
  idx = last(ℐ, 2)
  saved = unseen[idx]
  unseen[idx] .= false
  𝒞 = search(p[last(ℐ)], searcher; mask=unseen)
  unseen[idx] .= saved
  𝒞
end

# helper to find next point with smallest angle
jarvisnext(::Nothing, 𝒞, p, ℐ, A, O) = argmin(i -> ∠(A, O, p[i]), 𝒞)
function jarvisnext(::KNearestSearch, 𝒞, p, ℐ, A, O)
  # check candidates in order of increasing angle and accept the first one
  # whose segment does not cross the existing hull, skipping the last edge
  for nᵢ in sort(𝒞, by=i -> ∠(A, O, p[i]))
    cseg = Segment(p[ℐ[end]], p[nᵢ])
    cbox = boundingbox(cseg)
    tₒ = nᵢ == ℐ[begin] ? 2 : 1
    valid = !any(tₒ:(length(ℐ) - 2)) do t
      eseg = Segment(p[ℐ[t]], p[ℐ[t + 1]])
      # quick check to see if segments could intersect before doing more expensive segment intersection check
      intersects(cbox, boundingbox(eseg)) && intersects(cseg, eseg)
    end
    valid && return nᵢ
  end
  nothing
end

# helper to mark point as seen after it is added to the hull
jarvisupdate!(::Nothing, unseen, j) = nothing
jarvisupdate!(::KNearestSearch, unseen, j) = unseen[j] = false
