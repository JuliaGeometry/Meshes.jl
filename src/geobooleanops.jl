# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
"""
  GeoBoolean

types for geometric boolean operations.

  * GeoIntersection - regions inside both geometries

  * GeoUnion - regions inside at least one geometry

  * GeoDifference - regions inside geometry but not other

  * GeoSymDifference - regions inside exactly one geometry
"""
abstract type GeoBoolean end
struct GeoIntersection <: GeoBoolean end
struct GeoUnion <: GeoBoolean end
struct GeoDifference <: GeoBoolean end
struct GeoSymDifference <: GeoBoolean end

"""
    geobooleanop(geometry, other, operation)

    performs boolean `operation` from `geometry` with `other`

    ## Operands

    - `operation` — the operation to perform; one of:
      - `GeoIntersection` — regions inside both geometries (intersect/∩)
      - `GeoUnion` — regions inside at least one geometry (union/∪)
      - `GeoDifference` — regions inside `geometry` but not `other` (set difference/setdiff)
      - `GeoSymDifference` — regions inside exactly one geometry (symmetric difference/xor/⊻)

    ## Notes

    The algorithm works for both convex and concave polygons using
    Martinez-Rueda clipping.

    ## References
    * Martínez, F., Rueda, A.J., Feito, F.R. 2009. [A new algorithm for computing Boolean operations on
      polygons](https://doi.org/10.1016/j.cag.2009.03.003)

    """
function geobooleanop end
Base.union(a::Polygon, b::Polygon) = geobooleanop(a, b, GeoUnion())
Base.setdiff(a::Polygon, b::Polygon) = geobooleanop(a, b, GeoDifference())
Base.symdiff(a::Polygon, b::Polygon) = geobooleanop(a, b, GeoSymDifference())
Base.xor(a::Polygon, b::Polygon) = geobooleanop(a, b, GeoSymDifference())

function geobooleanop(poly, other, operation)
  poly = deepcopy(poly)
  other = deepcopy(other)
  geobooleanop!(poly, other, operation)
end
function geobooleanop!(ring::Ring, other::Ring, operation::GeoBoolean)
  # make sure other ring is CCW
  occw = orientation(other) == CCW ? other : reverse(other)

  allsegs = Iterators.flatten((segments(ring), segments(occw)))
  intersections, seginds = pairwiseintersect(allsegs)
  isnothing(intersections) && return nothing
  # insert intersections into the two rings
  vᵣ = vertices(ring)
  vₒ = vertices(occw)
  n0::Int = length(vᵣ)
  _insertintersections!(intersections, seginds, vᵣ, vₒ, n0)

  # build annotated segments from split rings
  fsₐ = FillSegments(vcat(collect(segments(ring)), collect(segments(occw))))
  n::Int = length(segments(ring))

  # annotate the polygon relationships of segments
  _annotate!(fsₐ, n)
  selected, fills = _selectsegments(fsₐ, n, operation)

  isempty(selected) && return nothing

  # deduplicate degenerate segments
  P = typeof(first(selected).start)
  dupsegs = similar(selected, 0)
  dupfills = similar(fills, 0)
  dupmap = Dict{Int,Int}()  # old index -> new index

  keymap = Dict{Tuple{P,P},Int}()

  for i in eachindex(selected)
    seg = selected[i]
    f = fills[i]
    s = (seg.start, seg.stop)
    sr = (seg.stop, seg.start)

    if haskey(keymap, s)
      dupmap[i] = keymap[s]
    elseif haskey(keymap, sr)
      dupmap[i] = keymap[sr]
    else
      push!(dupsegs, seg)
      push!(dupfills, f)
      idx = length(dupsegs)
      keymap[s] = idx
      dupmap[i] = idx
    end
  end

  # remove degenerate zero-length segments (exact equality)
  pairs = collect(zip(dupsegs, dupfills))
  filter!(p -> first(p).start != first(p).stop, pairs)

  selected = [seg for (seg, _) in pairs]
  fills = [fill for (_, fill) in pairs]

  # chain segments into closed rings
  function _buildrings(selected, fills, operation)
    P = typeof(first(selected).start)
    used = Set{Int}()
    rings = Vector{Tuple{Vector{P},Bool}}()  # (ring, is_hole)

    # sort so duplicate segments are at the end, can't start chains from them and lead to errors
    indices = collect(eachindex(selected))
    sort!(indices, by=idx -> begin
      fill = fills[idx]
      # segments with mismatched extfill
      fill.extabove == fill.extbelow ? 0 : 1
    end)

    for idxₛ in indices
      idxₛ in used && continue

      segₛ = selected[idxₛ]
      # skip zero-length segments
      segₛ.start == segₛ.stop && continue

      chain = [selected[idxₛ]]
      chainfills = [fills[idxₛ]]
      push!(used, idxₛ)

      while true
        chainₑ = chain[end].stop

        # find candidate segments that can continue the chain
        candidates = Tuple{Int,Bool}[]
        for (idx, seg) in enumerate(selected)
          if seg.start == chainₑ
            push!(candidates, (idx, true))  # segment starts here
          elseif seg.stop == chainₑ
            push!(candidates, (idx, false))  # segment stops here (would be reversed)
          end
        end

        # find best next segment based on fill pattern
        idxₙ = 0
        reversedₙ = false
        currfill = chainfills[end]

        priority = nothing
        for (candidx, isstart) in candidates
          if candidx in used
            continue
          end

          candfill = fills[candidx]
          candseg = selected[candidx]

          # prevent degenerate chains
          lastseg = chain[end]
          degenerate = false

          degenerate = isstart ? (candseg.stop == lastseg.start) : (candseg.start == lastseg.start)
          degenerate && length(chain) == 1 && continue

          # check if extfill matches
          extfillmatch = (candfill.extabove == currfill.extabove && candfill.extbelow == currfill.extbelow)
          if isnothing(priority) || extfillmatch > priority
            priority = extfillmatch
            idxₙ = candidx
            reversedₙ = !isstart  # if segment ends here, we need to reverse it
          end
        end

        idxₙ == 0 && break
        segₙ = selected[idxₙ]
        reversedₙ && (segₙ = (start=segₙ.stop, stop=segₙ.start))

        push!(chain, segₙ)
        push!(chainfills, fills[idxₙ])
        push!(used, idxₙ)

        chain[begin].start == chain[end].stop && break
      end

      # check if chain can close
      if chain[begin].start == chain[end].stop
        ring = [s.start for s in chain]
        # check if this ring represents a hole (for symdiff)
        # rings with extfill=true on both sides are holes
        firstfill = first(chainfills)
        ishole = if operation isa GeoSymDifference
          firstfill.extabove && firstfill.extbelow
        else
          false
        end
        push!(rings, (ring, ishole))
      else
        # chain didn't close naturally - this can happen when FP precision creates micro-segments
        # check if there's another unclosed chain that approximately connects to this one
        # for now: silently skip, printing the warning is slow but can easily happen
        # @warn "Unclosed chain" start = first(chain).start stop = last(chain).stop length = length(chain) note = "potentially due to floating point precision. May lead to open polygons, therefore possibly missing polygons or vertices."
      end
    end

    # separate outer rings from holes
    outers = [r for (r, hole) in rings if !hole]
    holes = [r for (r, hole) in rings if hole]

    # return as vector of rings (outer first, then holes)
    vcat(outers, holes)
  end

  rings = _buildrings(selected, fills, operation)

  (isempty(rings) && return nothing)
  Ring.(rings)
end

function geobooleanop(poly::Polygon, other::Geometry, operation::GeoBoolean)
  # copying variant: deepcopy inputs and delegate to mutating version
  poly = deepcopy(poly)
  other = deepcopy(other)
  geobooleanop!(poly, other, operation)
end

function geobooleanop!(poly::Polygon, other::Geometry, operation::GeoBoolean)
  # perform boolean clipping per ring (mutating variant delegates to ring-level mutator)
  c = [geobooleanop!(ring, boundary(other), operation) for ring in rings(poly)]

  # flatten clipped ring groups into a single vector of Ring
  groups = [g for g in c if !isnothing(g)]
  r = isempty(groups) ? Ring[] : vcat(groups...)
  isempty(r) && return nothing

  # determine which rings overlap
  # build overlap graph
  n = length(r)
  rows = Int[]
  cols = Int[]
  vals = Bool[]
  for i in 1:(n - 1)
    for j in (i + 1):n
      v₁ = vertices(r[i])
      v₂ = vertices(r[j])
      # check if any vertex of one ring is inside the other ring
      if any(p -> sideof(p, r[j]) == IN, v₁) || any(p -> sideof(p, r[i]) == IN, v₂)
        push!(rows, i, j)
        push!(cols, j, i)
        push!(vals, true, true)
      end
    end
  end
  # build a sparse boolean adjacency matrix
  adj = SparseArrays.sparse(rows, cols, vals, n, n)
  # connected components of the overlap graph => groups to merge into PolyArea
  visited = falses(n)
  components = Vector{Vector{Int}}()
  for i in 1:n
    if !visited[i]
      stack = [i]
      component = Int[]
      visited[i] = true
      while !isempty(stack)
        u = pop!(stack)
        push!(component, u)
        for v in 1:n
          if adj[u, v] && !visited[v]
            visited[v] = true
            push!(stack, v)
          end
        end
      end
      push!(components, component)
    end
  end

  # create PolyArea for each group:
  polyareas = PolyArea[]
  for inds in components
    grouprings = r[inds]
    if length(grouprings) == 1
      rr = first(grouprings)
      rr = orientation(rr) == CCW ? rr : reverse(rr)
      push!(polyareas, PolyArea([rr]))
    else
      # this seems to work, so long as outer polygon is legitimately the outside
      outer = first(grouprings)
      outer = orientation(outer) == CCW ? outer : reverse(outer)
      inners = [orientation(rr) == CW ? rr : reverse(rr) for rr in grouprings[2:end]]
      push!(polyareas, PolyArea(vcat([outer], inners)))
    end
  end

  # if multiple polyareas, return MultiPolygon
  if length(polyareas) == 1
    first(polyareas)
  else
    Multi(Tuple(polyareas))
  end
end

"""
    FillSegments{S}

container for segments with fill annotations for clipping algorithm.
  This allows tracking of which regions are filled by each polygon and segment ordering
  during the sweep line algorithm.

Fields:
- `segments`: Vector of Segment geometries
- `ownfillabove`: Vector of Bool/nothing for filled region by itself above
- `ownfillbelow`: Vector of Bool/nothing for filled region by itself below
- `extfillabove`: Vector of Bool/nothing for filled region by other polygon above
- `extfillbelow`: Vector of Bool/nothing for filled region by other polygon below
"""
mutable struct FillSegments{S}
  segments::Vector{S}
  ownfillabove::Vector{Union{Bool,Nothing}}
  ownfillbelow::Vector{Union{Bool,Nothing}}
  extfillabove::Vector{Union{Bool,Nothing}}
  extfillbelow::Vector{Union{Bool,Nothing}}
end

function FillSegments(segs)
  n::Int = length(segs)
  S = typeof(first(segs))
  fillvec = Vector{Union{Bool,Nothing}}(nothing, n)
  FillSegments{S}(collect(segs), copy(fillvec), copy(fillvec), copy(fillvec), copy(fillvec))
end

Base.collect(fs::FillSegments) = fs.segments
Base.length(fs::FillSegments) = length(fs.segments)
Base.iterate(fs::FillSegments) = iterate(fs.segments)

# getters
Base.getindex(fs::FillSegments, i::Int) = fs.segments[i]
ownfillabove(fs::FillSegments, i::Int) = fs.ownfillabove[i]
ownfillbelow(fs::FillSegments, i::Int) = fs.ownfillbelow[i]
extfillabove(fs::FillSegments, i::Int) = fs.extfillabove[i]
extfillbelow(fs::FillSegments, i::Int) = fs.extfillbelow[i]

# setters
setownfillabove!(fs::FillSegments, i::Int, val) = (fs.ownfillabove[i] = val)
setownfillbelow!(fs::FillSegments, i::Int, val) = (fs.ownfillbelow[i] = val)
setextfillabove!(fs::FillSegments, i::Int, val) = (fs.extfillabove[i] = val)
setextfillbelow!(fs::FillSegments, i::Int, val) = (fs.extfillbelow[i] = val)

"""
    segmentisless(fs::FillSegments, i::Int, j::Int)

Compare two segments by index for sweep line ordering, orders by y position.
Returns true if segment i should come before segment j in the status tree.
Note: this assumes segments are split at all intersection points.
"""
function segmentisless(fs::FillSegments, i::Int, j::Int)
  segi = fs[i]
  segj = fs[j]

  a₁, a₂ = sort(vertices(segi))
  b₁, b₂ = sort(vertices(segj))

  sa1 = sideof(a₁, Line(b₁, b₂))
  sa2 = sideof(a₂, Line(b₁, b₂))

  # RIGHT means a is below b
  if sa1 == ON
    if sa2 == ON
      # both collinear
      a₁ < b₁ || a₂ < b₂ || i > j
    else
      sa2 == RIGHT
    end
  else
    sa1 == sa2 ? sa1 == RIGHT : sideof(b₁, Line(a₁, a₂)) != RIGHT
  end
end

# top-level wrapper used by AVL status tree to compare segment ordering
struct SegmentIndex
  fs::FillSegments
  idx::Int
end
Base.isless(a::SegmentIndex, b::SegmentIndex) = segmentisless(a.fs, a.idx, b.idx)

# selection tables for boolean operations
# value: 0=discard, 1=keep filled above, 2=keep filled below
# matching the tables in https://sean.fun/a/polygon-clipping-pt2/#selecting-resulting-segments
const SELECTION = Dict(
  :union => Int8[0, 2, 1, 0, 2, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
  :intersection => Int8[0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 1, 1, 0, 2, 1, 0],
  :difference => Int8[0, 0, 0, 0, 2, 0, 2, 0, 1, 1, 0, 0, 0, 1, 2, 0],
  :xor => Int8[0, 2, 1, 0, 2, 0, 0, 1, 1, 0, 0, 2, 0, 1, 2, 0]
)
_retrieveoperation(::GeoIntersection) = :intersection
_retrieveoperation(::GeoUnion) = :union
_retrieveoperation(::GeoDifference) = :difference
_retrieveoperation(::GeoSymDifference) = :xor
_retrieveoperation(::GeoBoolean) = nothing

# ----------------
# HELPER FUNCTIONS
# ----------------

# add point event to event queue with associated segment event type
function addevent!(events, point, event)
  node = BinaryTrees.search(events, point)
  if isnothing(node)
    BinaryTrees.insert!(events, point, [event])
  else
    push!(BinaryTrees.value(node), event)
  end
end

# build event queue for two polygons' segments.
function buildevents(fs::FillSegments{S}, n::Int) where {S}
  P = eltype(vertices(first(fs)))
  isempty(fs) && return BinaryTrees.AVLTree{P,Vector{Tuple{Int,Bool,Bool}}}()

  events = BinaryTrees.AVLTree{P,Vector{Tuple{Int,Bool,Bool}}}()

  for i in 1:length(fs)
    fromᵣ = i <= n
    seg = fs[i]
    a, b = vertices(seg)
    if a > b
      a, b = b, a
    end

    addevent!(events, a, (i, true, fromᵣ))
    addevent!(events, b, (i, false, fromᵣ))
  end

  events
end

"""
    _annotate!(fs, n)

Annotate fill information for segments in `fs` using a sweep line algorithm.
  This performs one sweep to annotate where polygon regions are filled for
  both the source polygon and the other polygon.

For segments 1:n (from first polygon):
- ownfill computed from segments 1:n below
- extfill computed from segments n+1:end below

For segments n+1:end (from second polygon):
- ownfill computed from segments n+1:end below
- extfill computed from segments 1:n below

"""
function _annotate!(fs::FillSegments{S}, n::Int) where {S}
  isempty(fs) && return

  # build single event queue with origin tracking
  events = buildevents(fs, n)

  # status: AVL tree mapping SegmentIndex -> fromᵣ (which polygon?)
  status = BinaryTrees.AVLTree{SegmentIndex,Bool}()

  # process events in order
  while !BinaryTrees.isempty(events)
    node = BinaryTrees.minnode(events)
    p = BinaryTrees.key(node)
    eventsₐ = BinaryTrees.value(node)

    sort!(eventsₐ, by=ev -> (SegmentIndex(fs, ev[1]), !ev[3]))
    # remove end segments
    eventsₑ = Iterators.filter(ev -> !ev[2], eventsₐ)
    for (segind, _, _) in eventsₑ
      BinaryTrees.delete!(status, SegmentIndex(fs, segind))
    end

    eventsₛ = filter!(ev -> ev[2], eventsₐ)
    # first pass: compute ownfill for all segments
    for (segind, _, fromᵣ) in eventsₛ
      # find what's directly below this segment in the status
      below, _ = BinaryTrees.prevnext(status, SegmentIndex(fs, segind))

      # for ownfill: find the closest segment from SAME polygon below
      ownbelowidx = nothing
      current = below
      while !isnothing(current)
        candidx = BinaryTrees.key(current).idx
        candfromᵣ = BinaryTrees.value(current)
        if candfromᵣ == fromᵣ
          ownbelowidx = candidx
          break
        end
        current, _ = BinaryTrees.prevnext(status, BinaryTrees.key(current))
      end

      # compute ownfill based on same-polygon segment below
      if isnothing(ownbelowidx)
        # no segment from same polygon below: exterior
        setownfillbelow!(fs, segind, false)
        setownfillabove!(fs, segind, true)  # this segment toggles from false to true
      else
        # same polygon segment below: copy and toggle
        setownfillbelow!(fs, segind, ownfillabove(fs, ownbelowidx))
        setownfillabove!(fs, segind, !ownfillbelow(fs, segind))
      end

      # add to status immediately after computing ownfill
      BinaryTrees.insert!(status, SegmentIndex(fs, segind), fromᵣ)
    end

    # second pass: detect duplicate/collinear segments and exchange fill info
    for i in 1:length(eventsₛ)
      segindi, _, fromᵣi = eventsₛ[i]
      segi = fs[segindi]
      ai, bi = vertices(segi)

      # only check segments from first polygon
      fromᵣi || continue

      # look for matching/collinear segment from second polygon
      for j in 1:length(eventsₛ)
        segindj, _, fromᵣj = eventsₛ[j]

        # skip if same polygon
        fromᵣj == fromᵣi && continue

        segj = fs[segindj]
        aj, bj = vertices(segj)

        # check if segments match (same geometry within tolerance)
        # either ai==aj && bi==bj, or ai==bj && bi==aj (reversed)
        sharestart = ai == aj
        samedir = sharestart && bi == bj
        revdir = ai == bj && bi == aj

        # exchange extfill info if exact duplicate found
        if samedir || revdir
          # found exact duplicate - exchange extfill
          setextfillabove!(fs, segindi, ownfillabove(fs, segindj))
          setextfillbelow!(fs, segindi, ownfillbelow(fs, segindj))
          setextfillabove!(fs, segindj, ownfillabove(fs, segindi))
          setextfillbelow!(fs, segindj, ownfillbelow(fs, segindi))
          break
        elseif sharestart
          # exchange extfill for collinear segments sharing a start point
          linei = Line(ai, bi)
          if sideof(bj, linei) == ON
            setextfillabove!(fs, segindi, ownfillabove(fs, segindj))
            setextfillbelow!(fs, segindi, ownfillbelow(fs, segindj))
            setextfillabove!(fs, segindj, ownfillabove(fs, segindi))
            setextfillbelow!(fs, segindj, ownfillbelow(fs, segindi))
            break
          end
        end
      end
    end

    # third pass: compute extfill for non-duplicate segments
    for (segind, _, fromᵣ) in eventsₛ
      # skip if extfill was already set by duplicate exchange
      if isnothing(extfillbelow(fs, segind))
        # find what's directly below this segment in the status
        below, _ = BinaryTrees.prevnext(status, SegmentIndex(fs, segind))

        # for extfill: use the immediately below segment
        if isnothing(below)
          # no segment below: exterior to other polygon
          setextfillbelow!(fs, segind, false)
          setextfillabove!(fs, segind, false)
        else
          belowidx = BinaryTrees.key(below).idx
          belowfromᵣ = BinaryTrees.value(below)

          if fromᵣ == belowfromᵣ
            # same polygon: propagate extfill from below
            inside = extfillabove(fs, belowidx)
          else
            # different polygon: use ownfill from below
            inside = ownfillabove(fs, belowidx)
          end

          setextfillbelow!(fs, segind, inside)
          setextfillabove!(fs, segind, inside)
        end
      end
    end

    BinaryTrees.delete!(events, p)
  end
end

function _selectsegments(fsₐ::FillSegments, n::Int, operation)
  op = _retrieveoperation(operation)
  table = SELECTION[op]
  P = eltype(vertices(first(fsₐ)))
  selected = Vector{NamedTuple{(:start, :stop),Tuple{P,P}}}()
  fills = Vector{NamedTuple{(:ownabove, :ownbelow, :extabove, :extbelow),NTuple{4,Bool}}}()

  for i in 1:length(fsₐ)
    fromᵣ = i <= n

    # get annotations
    ownabove = ownfillabove(fsₐ, i)
    ownbelow = ownfillbelow(fsₐ, i)
    extabove = extfillabove(fsₐ, i)
    extbelow = extfillbelow(fsₐ, i)

    if op == :difference && !fromᵣ
      # for difference, segments from the second polygon should have own<->ext swapped
      ownabove, extabove = extabove, ownabove
      ownbelow, extbelow = extbelow, ownbelow
    end
    # build index from 4 boolean flags (treating as 4-bit number)
    # handle nothing values as false
    ownbitabove = isnothing(ownabove) ? false : ownabove
    ownbitbelow = isnothing(ownbelow) ? false : ownbelow
    extbitabove = isnothing(extabove) ? false : extabove
    extbitbelow = isnothing(extbelow) ? false : extbelow

    idx = (ownbitabove ? 8 : 0) + (ownbitbelow ? 4 : 0) + (extbitabove ? 2 : 0) + (extbitbelow ? 1 : 0) + 1

    # select segments based on table
    if table[idx] != 0
      s = fsₐ[i]
      start, stop = vertices(s)
      # for difference operations, segments from the "subtracted" polygon need to be reversed
      flip = (op == :difference && !fromᵣ)

      if flip
        start, stop = stop, start
      end

      push!(selected, (start=start, stop=stop))
      push!(fills, (ownabove=ownbitabove, ownbelow=ownbitbelow, extabove=extbitabove, extbelow=extbitbelow))
    end
  end

  (selected, fills)
end

function _insertintersections!(intersections, seginds, vᵣ, vₒ, n)
  # build map segidx -> [points]
  P = eltype(intersections)
  G = Dict{Int,Vector{P}}()
  for (p, segs) in zip(intersections, seginds)
    for s in segs
      push!(get!(G, s, P[]), p)
    end
  end

  # insert points into the appropriate vertex arrays in increasing segment order
  idxs = sort(collect(keys(G)))
  offset1 = 1
  offset2 = 1
  for segidx in idxs
    pts = G[segidx]
    isfirst = segidx <= n
    v = isfirst ? vᵣ : vₒ
    i = segidx + (isfirst ? offset1 : offset2 - n)

    ps = v[i - 1]
    pe = v[i]

    # remove any intersection that equals an endpoint
    filter!(p -> p != ps && p != pe, pts)
    isempty(pts) && continue

    sort!(pts)

    if i ≤ length(v)
      # insert points into the vertex array at the correct position
      insert!.(Ref(v), Ref(i), ps < pe ? reverse(pts) : pts)
    else
      # if inserting at end, just append to avoid wrapping issues
      append!(v, ps < pe ? pts : reverse(pts))
    end

    # update offset
    m = length(pts)
    if isfirst
      offset1 += m
    else
      offset2 += m
    end
  end
end
