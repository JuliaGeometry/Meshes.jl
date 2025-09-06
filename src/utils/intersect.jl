# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    bentleyottmann(segments; [digits])

Compute pairwise intersections between n `segments`
with `digits` precision in O((n+k)⋅log(n)) time using
Bentley-Ottmann sweep line algorithm.

By default, set `digits` based on the absolute
tolerance of the length type of the segments.

## References

* Bentley & Ottmann 1979. [Algorithms for reporting and counting
  geometric intersections](https://ieeexplore.ieee.org/document/1675432)
"""
function bentleyottmann(segments; digits=_digits(segments))
  # orient segments and round coordinates
  segs = map(segments) do seg
    a, b = coordround.(extrema(seg), digits=digits)
    a > b ? (b, a) : (a, b)
  end

  # retrieve relevant types
  P = typeof(first(first(segs)))
  S = Tuple{P,P}

  # event queue: stores points with associated sets of starting, ending, and crossing segments
  𝒬 = BinaryTrees.AVLTree{EventPoint{P,S}}()
  # lookup table mapping segments to their linear indices
  lookup = Dict{S,Int}()

  # add points and segments to event queue and lookup table
  buildqueue!(𝒬, lookup, segs)

  # initialize sweepline
  sweep = SweepLine(segs)
  # output dictionary (planar graph 𝐺)
  𝐺 = Dict{P,Vector{Int}}()

  # handle events
  handle!(𝐺, sweep, 𝒬, lookup, digits)

  (collect(keys(𝐺)), collect(values(𝐺)))
end

# -----------------------
# BentleyOttmann Handling
# -----------------------

function handle!(𝐺, sweep, 𝒬, lookup, digits)
  # types
  P = typeof(sweep.p)
  ℒ = lentype(P)
  # vector holding segments intersecting the current event point
  bundle = Vector{SweepSegment{P,ℒ}}()
  # holds segment indices for output
  inds = Set{Int}()

  _innerloops(𝐺, sweep, 𝒬, lookup, bundle, inds, digits)
  𝐺
end

function _innerloops(𝐺, sweep, 𝒬, lookup, bundle, inds, digits)
  P = typeof(sweep.p)
  S = Tuple{P,P}
  ℬₚ, ℰₚ, ℳₚ = (Set{S}() for _ in 1:3)
  while !BinaryTrees.isempty(𝒬)
    # pre allocate dictionary of intersection types
    # pull the next event point and associated data
    p, ℬₚ, ℰₚ, ℳₚ = _pullevent!(𝒬, ℬₚ, ℰₚ, ℳₚ)

    # active segments at current event
    activesegs = activesegments(sweep)

    # update the status structure with the current event point
    _handlestatus!(activesegs, ℬₚ, ℳₚ, ℰₚ, sweep, p)

    # build bundle of segments crossing the current event point
    bundle = _buildbundle!(bundle, ℬₚ, ℳₚ, sweep)
    sort!(bundle)

    _findintersections!(bundle, activesegs, 𝒬, p, ℰₚ, sweep, digits)

    # add intersection points and corresponding segment indices to 𝐺
    _buildoutput!(𝐺, bundle, inds, lookup, ℬₚ, ℰₚ, p)
  end
  𝐺
end

# ------------------------------------
# Status and event handlers
# ------------------------------------

# updates the status structure with the current event point
function _handlestatus!(activesegs, ℬₚ, ℳₚ, ℰₚ, sweep, p)
  for seg in Iterators.flatten((ℰₚ, ℳₚ))
    BinaryTrees.delete!(activesegs, SweepSegment(seg, sweep))
  end

  update!(sweep, p)

  for seg in Iterators.flatten((ℬₚ, ℳₚ))
    BinaryTrees.insert!(activesegs, SweepSegment(seg, sweep))
  end
end

function _handlebottom!(bundle, activesegs, 𝒬, p, digits)
  # check segment below lowest
  s′ = bundle[begin]

  nsₗ, _ = !isnothing(s′) ? BinaryTrees.prevnext(activesegs, s′) : (nothing, nothing)
  if !isnothing(nsₗ)
    _newevent!(𝒬, _sweepline(s′), bundle, p, _segment(s′), _keyseg(nsₗ), digits)
  end
end

function _handletop!(bundle, activesegs, 𝒬, p, digits)
  # check segment above highest
  s″ = bundle[end]

  _, nsᵤ = !isnothing(s″) ? BinaryTrees.prevnext(activesegs, s″) : (nothing, nothing)
  if !isnothing(nsᵤ)
    _newevent!(𝒬, _sweepline(s″), bundle, p, _segment(s″), _keyseg(nsᵤ), digits)
  end
end

# -----------------
# HELPER FUNCTIONS
# -----------------
# builds event queue
function buildqueue!(𝒬, lookup, segs)
  for (i, seg) in enumerate(segs)
    _addpoint!(𝒬, seg[1], [seg], :starts)
    _addpoint!(𝒬, seg[2], [seg], :ends)
    lookup[seg] = i
  end
end

# obtains current event point and segments
function _pullevent!(𝒬, ℬₚ, ℰₚ, ℳₚ)
  node = BinaryTrees.minnode(𝒬)
  ep = BinaryTrees.key(node)
  p = ep.p
  S = typeof(ℬₚ)
  ℬₚ = haskey(ep.types, :starts) ? ep.types[:starts] : S()
  ℰₚ = haskey(ep.types, :ends) ? ep.types[:ends] : S()
  ℳₚ = haskey(ep.types, :crossings) ? ep.types[:crossings] : S()
  # crosses that aren't endpoints (including them can lead to duplicates)
  setdiff!(ℳₚ, ℰₚ)
  BinaryTrees.delete!(𝒬, ep)

  p, ℬₚ, ℰₚ, ℳₚ
end

# builds bundle of segments crossing the current event point
function _buildbundle!(bundle, ℬₚ, ℳₚ, sweep)
  empty!(bundle)
  for seg in Iterators.flatten((ℬₚ, ℳₚ))
    push!(bundle, SweepSegment(seg, sweep))
  end
  bundle
end

# finds new intersections
function _findintersections!(bundle, activesegs, 𝒬, p, ℰₚ, sweep, digits)
  if isempty(bundle) # occurs at endpoints
    endseg = SweepSegment(first(ℰₚ), sweep) # segment at endpoint
    # check newly adjacent segments
    nsₗ, nsᵤ = BinaryTrees.prevnext(activesegs, endseg)
    isnothing(nsₗ) || isnothing(nsᵤ) || _newevent!(𝒬, sweep, bundle, p, _keyseg(nsₗ), _keyseg(nsᵤ), digits)
  else
    # check for intersections with adjacent segments below and above the current event point
    BinaryTrees.isempty(activesegs) || _handlebottom!(bundle, activesegs, 𝒬, p, digits)
    BinaryTrees.isempty(activesegs) || _handletop!(bundle, activesegs, 𝒬, p, digits)
  end
end

# builds output for the current event point
function _buildoutput!(𝐺, bundle, inds, lookup, ℬ, ℰ, p)
  # add intersection points and corresponding segment indices to 𝐺
  if length(bundle) > length(ℬ) # bundle only has ℬ unless p is an intersection
    # add start and crossing segments
    for s in bundle
      push!(inds, lookup[_segment(s)])
    end
    # add ending segments
    for seg in ℰ
      push!(inds, lookup[seg])
    end

    # add indices to output
    indᵥ = collect(inds)
    if haskey(𝐺, p)
      union!(𝐺[p], indᵥ)
    else
      𝐺[p] = indᵥ
    end
  end
  empty!(inds)
  nothing
end

# add new events to the event queue as needed
function _newevent!(𝒬, sweep, bundle, p, seg₁, seg₂, digits)
  intersection(Segment(seg₁), Segment(seg₂)) do I
    t = type(I)
    if t === Crossing || t === EdgeTouching
      i = coordround(get(I), digits=digits)
      if i ≈ p # helps with vertical+horizontal intersections
        push!(bundle, SweepSegment(seg₂, sweep))
        push!(bundle, SweepSegment(seg₁, sweep))
      elseif i > p # add to 𝒬, update existing point if needed
        _addpoint!(𝒬, i, [seg₁, seg₂], :crossings)
      end
    end
  end
end

# add an endpoint and corresponding segment to the event queue at a given position (1=start, 2=end, 3=crossing)
function _addpoint!(𝒬, p, segs, pos)
  P = typeof(p)
  S = Tuple{P,P}
  node = BinaryTrees.search(𝒬, EventPoint(p, S))
  # updates or adds event point based on existing events in 𝒬
  if !isnothing(node)
    update!(node.key, segs, pos)
  else
    # create a new event point for given position
    BinaryTrees.insert!(𝒬, EventPoint(p, segs, pos))
  end
end

# compute the number of significant digits based on the segment type
# this is used to determine the precision of the points
function _digits(segments)
  seg = first(segments)
  ℒ = lentype(seg)
  τ = ustrip(eps(ℒ))
  round(Int, 0.8 * (-log10(τ))) # 0.8 is a heuristic to avoid numerical issues
end

# compute y bounds of the segments
function _ybounds(segs)
  # compute bounding box
  bbox = boundingbox(segs)

  # stretch bounding bbox
  T = numtype(lentype(bbox))
  sbox = bbox |> Stretch(T(1.05))

  # extract y coordinate values
  map(p -> coords(p).y, extrema(sbox))
end

# convenience function to get the segment from a our AVLNode{SweepSegment} structure
_keyseg(ns) = _segment(BinaryTrees.key(ns))

# ----------------
# DATA STRUCTURES
# ----------------

# tracks the event point and constructs the sweepline lazily
mutable struct SweepLine{P<:Point,ℒ<:Number}
  p::P
  const ybounds::Tuple{ℒ,ℒ}
  activesegments::Any # sorted sequence of segments currently intersecting the sweepline
end

# getters
_sweepcoords(sweep::SweepLine) = CoordRefSystems.values(coords(sweep.p))
activesegments(sweep::SweepLine) = sweep.activesegments
update!(sweep::SweepLine, p) = (sweep.p = p; sweep)

#= Sweepline Definition ---

defined as the line segments
|
└——
 o|
  |
where o is the sweepline point and lines are the sweepline

This definition is used to
handle intersections elegantly. It ensures that y coordinates
of non vertical segments are always correct,
vertical segments are always on top of non-vertical segments next to p,
ending segments don't intersect, and all other segments are correctly ordered.

Inspired by LEDA implementation, but modified for Julia
=#
function sweepline(sweep::SweepLine)
  x, y = _sweepcoords(sweep)
  # perturbation to avoid numerical issues
  ϵ = atol(x) + eps(x)
  lower, upper = sweep.ybounds

  p₁ = Point(x + ϵ, lower) # lowest point
  p₂ = Point(x + ϵ, y + 2ϵ) # doubled to avoid precision issues
  p₃ = Point(x - ϵ, y + 2ϵ)
  p₄ = Point(x - ϵ, upper) # highest point
  # tuple of segments representing the sweepline
  # ((lower vertical), (horizontal), (upper vertical))
  ((p₁, p₂), (p₂, p₃), (p₃, p₄))
end

# sweepline intersection with segment
function _sweepintersect(seg::Tuple{P,P}, sweep::SweepLine{P,ℒ}) where {P<:Point,ℒ<:Number}
  rope = sweepline(sweep)
  I = nothing
  # check for intersections between the segment and the sweepline.
  # `Rope()` intersections are not type stable; this loop is a workaround.
  for segᵣ in rope
    I = intersection(Segment(seg), Segment(segᵣ)) do 𝑖
      t = type(𝑖)
      (t === Crossing || t === EdgeTouching || t === CornerTouching) ? get(𝑖) : nothing
    end
    if !isnothing(I)
      break
    end
  end

  # if I isn't found, then we are handling an ending segment, so seg[2] is used instead
  i = isnothing(I) ? seg[2] : I
  _, y = CoordRefSystems.values(coords(i))
  y
end

# tracks a segment and its intersection with the sweepline at the latest calculated event point
mutable struct SweepSegment{P<:Point,ℒ<:Number}
  const seg::Tuple{P,P}
  const sweepline::SweepLine{P,ℒ}
  yintersect::ℒ # information about the intersection with the sweepline
  latestpoint::P # latest point of the sweepline used to calculate the intersection
end

# constructor for SweepSegment using SweepLine
function SweepSegment(seg::Tuple{P,P}, sweep::SweepLine{P,ℒ}) where {P<:Point,ℒ<:Number}
  y = _sweepintersect(seg, sweep)
  SweepSegment{P,ℒ}(seg, sweep, y, sweep.p)
end
# extend SweepLine for SweepSegment
function SweepLine(segs)
  P = typeof(first(first(segs)))
  ℒ = lentype(P)
  ybounds = _ybounds(segs)
  pmin, _ = first(extrema(segs))
  SweepLine{P,ℒ}(pmin, ybounds, BinaryTrees.AVLTree{SweepSegment{P,ℒ}}())
end

# getters for SweepSegment
_segment(s::SweepSegment) = getfield(s, :seg)
_sweepline(s::SweepSegment) = getfield(s, :sweepline)

Base.:(==)(s₁::SweepSegment, s₂::SweepSegment) = s₁.seg == s₂.seg

# compare two segments based on their sweepline intersection relative to the current event point.
function Base.isless(s₁::SweepSegment{P,ℒ}, s₂::SweepSegment{P,ℒ}) where {P<:Point,ℒ}
  # if segments same, return false
  if s₁ == s₂
    return false
  end
  #* calculating the y intersect is the largest performance bottleneck
  ya = s₁.yintersect # s₁ is always up-to-date
  yb = _ycalc!(s₂)

  diff = ustrip(abs(ya - yb))
  tol = eps(ℒ)
  # if segments are separated over y, check ya < yb
  if diff > tol
    ya < yb
  else
    # fallback to lexicographic ordering of segments
    seg₁, seg₂ = s₁.seg, s₂.seg
    seg₁ < seg₂
  end
end

# calculate y-coordinate of intersection with sweepline
function _ycalc!(s::SweepSegment{P,ℒ}) where {P<:Point,ℒ<:Number}
  sweep = s.sweepline
  # if the latest point is the sweepline point, use the precalculated intersection
  if s.latestpoint === sweep.p
    y = s.yintersect
  else
    # otherwise, calculate the intersection with the sweepline
    # and update
    y = convert(ℒ, _sweepintersect(s.seg, sweep))
    s.latestpoint = sweep.p
    s.yintersect = y
  end
  y
end

struct EventPoint{P,S}
  p::P
  types::Dict{Symbol,Set{S}}
end

# empty point constructor
function EventPoint(p::P, ::Type{S}) where {P,S}
  types = Dict{Symbol,Set{S}}()
  EventPoint{P,S}(p, types)
end

function EventPoint(p::P, segs::Vector{S}, pos::Symbol) where {P,S}
  types = Dict{Symbol,Set{S}}()
  # Add the position symbol and set
  types[pos] = Set(segs)
  EventPoint{P,S}(p, types)
end

Base.isless(ep₁::EventPoint, ep₂::EventPoint) = ep₁.p < ep₂.p

function update!(ep::EventPoint{P,S}, segs, pos::Symbol) where {P,S}
  set = get!(ep.types, pos, Set{S}())
  foreach(s -> push!(set, s), segs)
end
