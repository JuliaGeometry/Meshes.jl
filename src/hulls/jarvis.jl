# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    JarvisMarch()
    JarvisMarch(k)

Compute the convex hull of a set of points or geometries using the
Jarvis's march algorithm. See [https://en.wikipedia.org/wiki/Gift_wrapping_algorithm]
(https://en.wikipedia.org/wiki/Gift_wrapping_algorithm).

If `k` is provided, the algorithm will attempt to compute a concave hull using the
k nearest neighbors as proposed by Moreira & Santos 2007. The value of `k` must be
greater than 2 and less than the number of unique points.

The algorithm has complexity `O(n*h)` where `n` is the number of points
and `h` is the number of points in the hull.

## References

* Jarvis 1973. [On the identification of the convex hull of a finite set of
  points in the plane](https://www.sciencedirect.com/science/article/abs/pii/0020019073900203)
* Moreira, A. & Santos, M. Y. 2007. "concave hull: a k-nearest neighbours approach for the computation of the region occupied by a set of points". In: *Proceedings of the Second International Conference on Computer Graphics Theory and Applications*. SciTePress, pp. 61–68.
"""
struct JarvisMarch{K} <: HullMethod
  k::K
end

JarvisMarch() = JarvisMarch{Nothing}(nothing)

JarvisMarch(k::I) where {I<:Integer} = JarvisMarch{I}(k)

function hull(points, method::JarvisMarch)
  pₒ = first(points)
  ℒ = lentype(pₒ)

  # sanity check
  ncoords = CoordRefSystems.ncoords(coords(pₒ))
  assertion(ncoords == 2, "Jarvis's march algorithm is only defined with 2D coordinates")

  # remove duplicates
  p = unique(points)
  n = length(p)

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])
  k = method.k
  !isnothing(k) && assertion(2 < k < n, "k must be greater than 2 and less than the number of points")

  # find bottom-left point
  i = argmin(p)
  start = i
  # initialize hull with i
  ℐ = [i]

  # initialize searcher and mask of visited points for k-nearest neighbors if needed
  searcher, pointmask = jarvissearcher(k, n, p)

  # find neighbor candidates
  𝒞 = jarviscandidates(searcher, pointmask, n, p, i, i, start)

  # find next point with smallest angle
  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  j = jarvisnext(searcher, 𝒞, p, ℐ, A, O)

  # initialize ring of indices
  push!(ℐ, j)
  jarvisupdate!(searcher, pointmask, j)

  # rotational sweep
  while first(ℐ) != last(ℐ)
    # direction of current segment
    v = p[j] - p[i]

    # find candidates for next point, excluding endpoints of current segment
    𝒞 = jarviscandidates(searcher, pointmask, n, p, j, (i, j), start)
    # no candidates, should only happen if k is too small
    isempty(𝒞) && throw(ArgumentError("could not find concave hull with k = $k, try a larger k"))

    # find next segment
    i = j
    O = p[i]
    A = O + v
    j = jarvisnext(searcher, 𝒞, p, ℐ, A, O)
    # no valid next point, should only happen if k is too small
    isnothing(j) && throw(ArgumentError("could not find concave hull with k = $k, try a larger k"))

    # update ring of indices
    push!(ℐ, j)
    jarvisupdate!(searcher, pointmask, j)
  end

  # return polygonal area
  PolyArea(p[ℐ[begin:(end - 1)]])
end

# helper to find next point for convex hull
jarvisnext(::Nothing, 𝒞, p, ℐ, A, O) = argmin(l -> ∠(A, O, p[l]), 𝒞)

function jarvisnext(::KNearestSearch, 𝒞, p, ℐ, A, O)
  sort!(𝒞, by=l -> ∠(A, O, p[l]))
  # accept first candidate whose segment does not cross the existing hull,
  # skipping the last edge, which shares a vertex with the candidate segment,
  # and the first edge when the candidate is the starting point
  for nᵢ in 𝒞
    cseg = Segment(p[ℐ[end]], p[nᵢ])
    cbox = boundingbox(cseg)
    tₒ = nᵢ == ℐ[begin] ? 2 : 1
    valid = !any(tₒ:(length(ℐ) - 2)) do t
      eseg = Segment(p[ℐ[t]], p[ℐ[t + 1]])
      #check if segments could intersect before doing more expensive segment intersection check
      intersects(cbox, boundingbox(eseg)) && intersects(cseg, eseg)
    end
    valid && return nᵢ
  end
  nothing
end

# helper to get candidate indices for next point
jarviscandidates(::Nothing, pointmask, n, p, current, inds, start) = setdiff(1:n, inds)

function jarviscandidates(searcher::KNearestSearch, pointmask, n, p, current, inds, start)
  # mask out points already in the hull, except for the starting point
  mask = .!pointmask
  mask[start] = true
  for ind in inds
    mask[ind] = false
  end
  search(p[current], searcher; mask=mask)
end

# helper to mark point as visited after it is added to the hull
jarvisupdate!(::Nothing, pointmask, j) = nothing
jarvisupdate!(::KNearestSearch, pointmask, j) = pointmask[j] = true

# helper to create searcher and mask of visited points
jarvissearcher(k::Integer, n, p) = KNearestSearch(p, k), falses(n)
jarvissearcher(k::Nothing, n, p) = nothing, 1:n
