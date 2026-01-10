# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
# script for reconstructing polygons from selected segments
# build polygons from selected segments
function _buildpolys(segments, fills, operation)
  rings, ringinds = _extractrings(segments)
  outer, inner = _classifyrings(rings, ringinds, segments, fills, operation)
  _constructpolys(outer, inner)
end

# extract rings from segments
function _extractrings(segments)
  adj = _buildadjacency(segments)

  P = eltype(eltype(segments))
  rings = Vector{Vector{P}}()
  ringinds = Vector{Vector{Int}}()
  used = falses(length(segments))

  for startind in 1:length(segments)
    used[startind] && continue

    chain, inds = _tracechain(startind, segments, adj, used)

    # only keep rings with more than 2 points
    if length(chain) > 2
      push!(rings, chain)
      push!(ringinds, inds)
    end
  end

  rings, ringinds
end

# trace a chain of segments starting from a starting segment
function _tracechain(startind, segments, adj, used)
  P = eltype(eltype(segments))
  chain = P[]
  inds = Int[]

  curr = startind
  isstart = true

  while !used[curr]
    used[curr] = true
    seg = segments[curr]

    # orient segment in direction of traversal
    prevpt, currpt = isstart ? (seg[1], seg[2]) : (seg[2], seg[1])
    push!(chain, currpt)
    push!(inds, curr)

    # find next segment using smallest CCW angle
    next, nextstart = _findnext(currpt, prevpt, adj, segments, used)
    isnothing(next) && break

    curr = next
    isstart = nextstart
  end

  chain, inds
end

# find next segment to traverse based on smallest CCW angle
function _findnext(currpt, prevpt, adj, segments, used)
  candidates = adj[currpt]
  best = nothing
  beststart = true
  bestangle = Inf * u"rad"

  for (i, isstart) in candidates
    used[i] && continue

    seg = segments[i]
    nextpt = isstart ? seg[2] : seg[1]

    angle = ∠(prevpt, currpt, nextpt)
    # set domain to [0, 2π)
    angle < 0u"rad" && (angle += 2π * u"rad")

    if angle < bestangle
      bestangle = angle
      best = i
      beststart = isstart
    end
  end

  best, beststart
end

# classify rings into outer and inner based on fill information
function _classifyrings(rings, ringinds, segments, fills, operation)
  # build Ring objects
  Grings = [Ring(r) for r in rings]
  ℛ = eltype(Grings)

  outer = Vector{ℛ}()
  inner = Vector{ℛ}()

  for (ring, ind) in zip(Grings, ringinds)
    validsegments = @view segments[ind]
    validfills = @view fills[ind]
    if _ishole(validsegments, validfills, operation)
      push!(inner, ring)
    else
      orientation(ring) == CW && reverse!(ring)
      push!(outer, ring)
    end
  end

  outer, inner
end

# construct polygons from outer rings and their holes
function _constructpolys(outer, inner)
  # map each hole to its containing outer ring
  children = [eltype(inner)[] for _ in outer]

  for hole in inner
    orientation(hole) == CCW && reverse!(hole)

    for (i, out) in enumerate(outer)
      if any(v -> sideof(v, out) ∈ (IN, ON), vertices(hole))
        push!(children[i], hole)
        break
      end
    end
  end

  polys = map(zip(outer, children)) do (out, kids)
    PolyArea([out; kids])
  end

  # return single polygon if only one contained polygon, else Multi
  if length(polys) == 1
    first(polys)
  else
    Multi(polys)
  end
end

# comparison for finding leftmost topmost segment. used in hole detection
function _lefttopmost(s₁, s₂)
  a₁, b₁ = s₁
  a₂, b₂ = s₂

  k₁, p₁ = _toppoint(a₁, b₁)
  k₂, p₂ = _toppoint(a₂, b₂)

  if k₁ != k₂
    k₁ > k₂
  else
    # share top point, use orientation to determine which is leftmost
    other₁ = a₁ == p₁ ? b₁ : a₁
    other₂ = a₂ == p₂ ? b₂ : a₂
    orientation(Ring(p₁, other₁, other₂)) != CW
  end
end

# return topmost point
function _toppoint(a, b)
  ca, cb = coords(flat(a)), coords(flat(b))
  # prefer higher y, then rightmost x
  if ca.y > cb.y || (ca.y == cb.y && ca.x > cb.x)
    (ca.y, ca.x), a
  else
    (cb.y, cb.x), b
  end
end

# build adjacency map from points to segments
function _buildadjacency(segments)
  P = eltype(eltype(segments))
  adj = Dict{P,Vector{Tuple{Int,Bool}}}()

  for (i, seg) in enumerate(segments)
    a, b = seg
    push!(get!(() -> Tuple{Int,Bool}[], adj, a), (i, true))
    push!(get!(() -> Tuple{Int,Bool}[], adj, b), (i, false))
  end

  adj
end

function _ishole(validsegments, validfills, operation)
  # first get top-leftmost segment without sorting (to avoid allocations)
  refsegind = 1
  for i in 2:length(validsegments)
    if _lefttopmost(validsegments[i], validsegments[refsegind])
      refsegind = i
    end
  end

  # get fill bits for this segment
  bits = validfills[refsegind]
  # if filled above for subject or clip, based on operation, it's a hole
  _filled(operation, bits, true)
end
