# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
# -------------------------
# TYPES
# -------------------------

"""
    FillSegments{S}

Container for segments with fill annotations for geometric boolean algorithm.
"""
mutable struct FillSegments{S}
  segments::Vector{S}
  # fill status (subject/clip, filled above/below)
  subjectabove::BitVector
  subjectbelow::BitVector
  clipabove::BitVector
  clipbelow::BitVector
  # segment has been processed?
  computed::BitVector
end

function FillSegments(segments)
  n = length(segments)
  FillSegments(collect(segments), falses(n), falses(n), falses(n), falses(n), falses(n))
end

Base.collect(fs::FillSegments) = fs.segments
Base.length(fs::FillSegments) = length(fs.segments)
Base.getindex(fs::FillSegments, i::Int) = fs.segments[i]
Base.iterate(fs::FillSegments, state...) = iterate(fs.segments, state...)

setfills!(fs::FillSegments, above, below, i, issubject) =
  issubject ? _setsubjectfills!(fs, above, below, i) : _setclipfills!(fs, above, below, i)

_setsubjectfills!(fs::FillSegments, above, below, i) = (fs.subjectabove[i] = above; fs.subjectbelow[i] = below)

_setclipfills!(fs::FillSegments, above, below, i) = (fs.clipabove[i] = above; fs.clipbelow[i] = below)
getfills(fs::FillSegments, i, issubject) =
  issubject ? (fs.subjectabove[i], fs.subjectbelow[i]) : (fs.clipabove[i], fs.clipbelow[i])

_getsubjectfills(fs::FillSegments, i) = (fs.subjectabove[i], fs.subjectbelow[i])
_getclipfills(fs::FillSegments, i) = (fs.clipabove[i], fs.clipbelow[i])

_setcomputed!(fs::FillSegments, i::Int) = (fs.computed[i] = true)

"""
    segmentisless(fs, i, j)

Compare two segments by index for sweep line ordering.
Only works for segments that have been split at intersections.
"""
function segmentisless(fs::FillSegments, i::Int, j::Int)
  segi, segj = fs[i], fs[j]
  a1, a2 = sort(vertices(segi))
  b1, b2 = sort(vertices(segj))

  s1 = sideof(a1, Line(b1, b2))
  s2 = sideof(a2, Line(b1, b2))

  # RIGHT means a is below b
  if s1 == ON
    # if collinear, use endpoints, then index to break tie, otherwise use side of second point and b2 below a2
    s2 == ON ? (a1 < b1 || a2 < b2 || i > j) : s2 == RIGHT
  else
    # if both on same side, return that, otherwise use side of b
    s1 == s2 ? s1 == RIGHT : sideof(b1, Line(a1, a2)) != RIGHT
  end
end

# wrapper for AVL tree key sorting
struct SegmentIndex{S}
  fs::FillSegments{S}
  ind::Int
end

Base.isless(a::SegmentIndex, b::SegmentIndex) = segmentisless(a.fs, a.ind, b.ind)

# helper for storing event points
struct SegmentEvent
  ind::Int
  isstart::Bool
  issubject::Bool
end

# ----------------
# HELPERS
# ----------------

# insert intersection points into rings
function _insertintersections!(intersections, seginds, rings::Vector{Vector{P}}) where {P}
  # group intersections by segment index
  G = Dict{Int,Vector{P}}()
  for (p, segs) in zip(intersections, seginds)
    for s in segs
      push!(get!(G, s, P[]), p)
    end
  end

  sortedseginds = sort(collect(keys(G)))

  # Precompute offsets
  ℒ = length.(rings)
  offsets = [0; cumsum(ℒ)]

  # Track inserted points per ring
  insertcounts = zeros(Int, length(rings))

  for ind in sortedseginds
    # find ring index
    rind = searchsortedfirst(offsets, ind) - 1
    lind = ind - offsets[rind]

    pts = G[ind]
    v = rings[rind]

    n₀ = ℒ[rind]
    startind = lind + insertcounts[rind]

    ps = v[startind]
    pe = lind < n₀ ? v[startind + 1] : v[1]

    filter!(p -> !isapprox(p, ps) && !isapprox(p, pe), pts)
    isempty(pts) && continue

    sort!(pts)
    ps > pe && reverse!(pts)

    if lind < n₀
      for (offset, pt) in enumerate(pts)
        insert!(v, startind + offset, pt)
      end
    else
      append!(v, pts)
    end

    insertcounts[rind] += length(pts)
  end
end

# build event queue from segments
function buildevents(fillsegs::FillSegments{S}, nsegsfirst::Int) where {S}
  P = eltype(vertices(first(fillsegs)))
  events = BinaryTrees.AVLTree{P,Vector{SegmentEvent}}()
  isempty(fillsegs) && return events

  for i in 1:length(fillsegs)
    issubject = i <= nsegsfirst
    s = fillsegs[i]
    a, b = vertices(s)
    a > b && ((a, b) = (b, a))
    _addevent!(events, a, SegmentEvent(i, true, issubject))
    _addevent!(events, b, SegmentEvent(i, false, issubject))
  end
  events
end

function _addevent!(events, p, ev)
  node = BinaryTrees.search(events, p)
  isnothing(node) ? BinaryTrees.insert!(events, p, [ev]) : push!(BinaryTrees.value(node), ev)
end
# select segments based on fill information and operation
function _selectsegments(fillsegs, nsegsfirst, operation)
  P = eltype(vertices(first(fillsegs)))

  selected = Vector{Tuple{P,P}}()
  fills = Vector{UInt8}()
  seen = Set{Tuple{P,P}}()

  for i in 1:length(fillsegs)
    a, b = vertices(fillsegs[i])
    a, b = a < b ? (a, b) : (b, a)
    if (a, b) ∈ seen
      continue
    end

    bits = _filltobits(fillsegs, i)

    # check if this segment is included based on operation
    isfilledabove = _filled(operation, bits, true)
    isfilledbelow = _filled(operation, bits, false)

    if isfilledabove ⊻ isfilledbelow
      push!(seen, (a, b))
      push!(selected, (a, b))
      push!(fills, bits)
    end
  end
  selected, fills
end

function _filled(op::PolygonBoolean, bits, above)
  masksubj = above ? SUBJTOP : SUBJBOTTOM
  maskclip = above ? CLIPTOP : CLIPBOTTOM

  # any subject/clip filled?
  subjfilled = (bits & masksubj) != 0
  clipfilled = (bits & maskclip) != 0

  _finalize(op, subjfilled, clipfilled)
end

_finalize(::PolyIntersection, s, c) = s && c
_finalize(::PolyUnion, s, c) = s || c
_finalize(::PolyDifference, s, c) = s && !c
_finalize(::PolySymDifference, s, c) = s ⊻ c

# annotate fill segments with fill information using a sweep line algorithm
function _annotatefill!(fillsegs::FillSegments{S}, nsegsfirst::Int) where {S}
  isempty(fillsegs) && return
  # build event queue
  events = buildevents(fillsegs, nsegsfirst)
  status = BinaryTrees.AVLTree{SegmentIndex{S},Bool}()

  while !BinaryTrees.isempty(events)
    node = BinaryTrees.minnode(events)
    p, evs = BinaryTrees.key(node), BinaryTrees.value(node)

    # sort: segments ending first, then by index
    sort!(evs, by=ev -> (SegmentIndex{S}(fillsegs, ev.ind), !ev.isstart))

    # remove ending segments
    for ev in evs
      !ev.isstart && BinaryTrees.delete!(status, SegmentIndex{S}(fillsegs, ev.ind))
    end

    # process start events
    startevs = filter(ev -> ev.isstart, evs)

    # pass 1: compute subjectfill (which direction is inside the subject polygon)
    for ev in startevs
      ind, issubject = ev.ind, ev.issubject
      below, _ = BinaryTrees.prevnext(status, SegmentIndex{S}(fillsegs, ind))

      # find closest segment from same polygon below
      belowind = nothing
      curr = below
      while !isnothing(curr)
        if BinaryTrees.value(curr) == issubject
          belowind = BinaryTrees.key(curr).ind
          break
        end
        curr, _ = BinaryTrees.prevnext(status, BinaryTrees.key(curr))
      end

      if isnothing(belowind) # exterior segment, set filled above, below = false
        setfills!(fillsegs, true, false, ind, issubject)
      else
        # otherwise reverse fill info from below segment
        fillabove, fillbelow = getfills(fillsegs, belowind, issubject)
        setfills!(fillsegs, fillbelow, fillabove, ind, issubject)
      end
      BinaryTrees.insert!(status, SegmentIndex{S}(fillsegs, ind), issubject)
    end

    # pass 2: compute fill information relative to the other polygon for remaining segments
    for ev in startevs
      ind, issubject = ev.ind, ev.issubject
      fillsegs.computed[ind] && continue

      below, _ = BinaryTrees.prevnext(status, SegmentIndex{S}(fillsegs, ind))
      if isnothing(below)
        # outside clip polygon
        setfills!(fillsegs, false, false, ind, !issubject)
      else
        bind = BinaryTrees.key(below).ind
        bissubject = BinaryTrees.value(below)
        # inside other polygon?
        if issubject == bissubject
          # below segment is from same polygon, use its status relative to OTHER
          inside = getfills(fillsegs, bind, !issubject)[1] # above value from ext polygon
        else
          # below segment is from other polygon, use its own fill
          inside = getfills(fillsegs, bind, bissubject)[1] # above value from its own polygon
        end
        # mark segment fill accordingly
        setfills!(fillsegs, inside, inside, ind, !issubject)
      end
      # mark as computed
      _setcomputed!(fillsegs, ind)
    end
    # pass 3: exchange fill information for duplicate segments
    for i in 1:length(startevs)
      evi = startevs[i]
      issubjecti = evi.issubject
      !issubjecti && continue # only check from subject polygon

      indi = evi.ind
      ai, bi = vertices(fillsegs[indi])

      for j in 1:length(startevs)
        evj = startevs[j]
        issubjectj = evj.issubject
        issubjecti == issubjectj && continue # skip subject polygon

        indj = evj.ind
        aj, bj = vertices(fillsegs[indj])

        sharestart = ai == aj
        samedir = sharestart && bi == bj
        revdir = ai == bj && bi == aj

        # if segments are collinear and sharing start point, exchange fill info
        if samedir || revdir || (sharestart && sideof(bj, Line(ai, bi)) == ON)
          _exchangefills!(fillsegs, indi, indj)
          break
        end
      end
    end

    BinaryTrees.delete!(events, p)
  end
end

# exchange fill information between two equal segments
function _exchangefills!(fillsegs, i, j)
  # segment i is from subject polygon
  # segment j is from clip polygon
  # they overlap, so they should share fill information

  # get what each segment knows about its own polygon
  sai, sbi = getsubjectfills(fillsegs, i)
  caj, cbj = getclipfills(fillsegs, j)

  # segment i needs to know about clip polygon (copy from j)
  setclipfills!(fillsegs, caj, cbj, i)

  # segment j needs to know about subject polygon (copy from i)
  setsubjectfills!(fillsegs, sai, sbi, j)

  _setcomputed!(fillsegs, i)
  _setcomputed!(fillsegs, j)
end

# comparison for finding leftmost topmost segment
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

# return topmost point and which point it is
function _toppoint(a, b)
  ca, cb = coords(flat(a)), coords(flat(b))
  # prefer higher y, then rightmost x
  if ca.y > cb.y || (ca.y == cb.y && ca.x > cb.x)
    (ca.y, ca.x), a
  else
    (cb.y, cb.x), b
  end
end

# build rings from selected segments
function _buildrings(segments, fills, operation)
  rings, ringinds = _extractrings(segments)
  outer, inner = _classifyrings(rings, ringinds, segments, fills, operation)
  _constructpolys(outer, inner)
end

function _extractrings(segments)
  adj = _buildadjacency(segments)

  P = eltype(eltype(segments))
  rings = Vector{Vector{P}}()
  ringinds = Vector{Vector{Int}}()
  used = falses(length(segments))

  for i in 1:length(segments)
    used[i] && continue

    chain, inds = _tracechain(i, segments, adj, used)

    # only keep rings with more than 2 points
    if length(chain) > 2
      push!(rings, chain)
      push!(ringinds, inds)
    end
  end

  rings, ringinds
end

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
    p₁, p₂ = isstart ? (seg[1], seg[2]) : (seg[2], seg[1])
    push!(chain, p₂)
    push!(inds, curr)

    # find next segment using smallest CCW angle
    next, nextstart = _findnext(p₂, p₁, adj, segments, used)
    isnothing(next) && break

    curr = next
    isstart = nextstart
  end

  chain, inds
end

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

  if length(polys) == 1
    first(polys)
  else
    Multi(polys)
  end
end

function _ishole(segments, fills, operation)
  # first get top-leftmost segment without sorting
  refsegind = 1
  for i in 2:length(segments)
    if _lefttopmost(segments[i], segments[refsegind])
      refsegind = i
    end
  end

  bits = fills[refsegind]
  # if filled above for subject or clip, based on operation, it's a hole
  _filled(operation, bits, true)
end
