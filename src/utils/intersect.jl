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
  ℒ = lentype(P)
  S = Tuple{P,P}
  U = Set{S}

  # event queue: stores points with associated sets of starting, ending, and crossing segments
  𝒬 = BinaryTrees.AVLTree{P,Tuple{U,U,U}}()

  # status structure: stores segments currently intersecting the sweepline
  ℛ = BinaryTrees.AVLTree{SweepSegment{P,ℒ}}()

  # lookup table mapping segments to their linear indices
  lookup = Dict{S,Int}()

  # add points and segments to event queue and lookup table
  for (i, seg) in enumerate(segs)
    _addpoint!(𝒬, seg[1], [seg], 1)
    _addpoint!(𝒬, seg[2], [seg], 2)
    lookup[seg] = i
  end

  # initialize sweepline
  pmin, _ = first(extrema(segs))
  ybounds = _ybounds(segs)
  sweepline = SweepLine{P,ℒ}(pmin, ybounds)
  # output dictionary (planar graph 𝐺)
  𝐺 = Dict{P,Vector{Int}}()
  # vector holding segments intersecting the current event point
  bundle = Vector{SweepSegment{P,ℒ}}()
  # holds segment indices for output
  inds = Set{Int}()

  # sweep line algorithm
  while !BinaryTrees.isempty(𝒬)
    # current event point
    node = BinaryTrees.minnode(𝒬)
    p = BinaryTrees.key(node)
    BinaryTrees.delete!(𝒬, p)
    # sets of beginning, ending, and crossing segments that include the current event point
    ℬₚ, ℰₚ, ℳₚ = BinaryTrees.value(node)
    # crosses that aren't endpoints (including them can lead to duplicates)
    setdiff!(ℳₚ, ℰₚ)
    # update the status structure with the current event point
    _handlestatus!(ℛ, ℬₚ, ℳₚ, ℰₚ, sweepline, p)

    # build bundle of segments crossing the current event point
    bundle = empty!(bundle)
    for seg in Iterators.flatten((ℬₚ, ℳₚ))
      push!(bundle, SweepSegment(seg, sweepline))
    end
    sort!(bundle)
    # process bundled events
    if isempty(bundle) # occurs at endpoints
      # check newly adjacent segments
      nsₗ, nsᵤ = BinaryTrees.prevnext(ℛ, SweepSegment(first(ℰₚ), sweepline))
      isnothing(nsₗ) || isnothing(nsᵤ) || _newevent!(𝒬, sweepline, bundle, p, _keyseg(nsₗ), _keyseg(nsᵤ), digits)
    else
      # check for intersections with adjacent segments below and above the current event point
      BinaryTrees.isempty(ℛ) || _handlebottom!(bundle, ℛ, 𝒬, p, digits)
      BinaryTrees.isempty(ℛ) || _handletop!(bundle, ℛ, 𝒬, p, digits)
    end

    # add intersection points and corresponding segment indices to 𝐺
    if length(bundle) > length(ℬₚ) # bundle only has ℬₚ unless p is an intersection
      # add start and crossing segments
      for s in bundle
        push!(inds, lookup[_segment(s)])
      end
      # add ending segments
      for seg in ℰₚ
        push!(inds, lookup[seg])
      end

      # add indices to output
      indᵥ = collect(inds)
      if haskey(𝐺, p)
        union!(𝐺[p], indᵥ)
      else
        𝐺[p] = indᵥ
      end
      empty!(inds)
    end
  end
  (collect(keys(𝐺)), collect(values(𝐺)))
end

# ------------------------------------
# Sweep line status and event handling
# ------------------------------------

# updates the status structure with the current event point
function _handlestatus!(ℛ, ℬₚ, ℳₚ, ℰₚ, sweepline, p)
  # remove end segments that are no longer active and crossings to update
  for seg in Iterators.flatten((ℰₚ, ℳₚ))
    BinaryTrees.delete!(ℛ, SweepSegment(seg, sweepline))
  end
  # update sweepline
  sweepline.p = p
  # insert new and crossing segments into the status structure
  for seg in Iterators.flatten((ℬₚ, ℳₚ))
    BinaryTrees.insert!(ℛ, SweepSegment(seg, sweepline))
  end
end

function _handlebottom!(bundle, ℛ, 𝒬, p, digits)
  # bundle is sorted sequence, so the first segment is minimum
  s′ = bundle[begin]
  # element below s′
  nsₗ, _ = !isnothing(s′) ? BinaryTrees.prevnext(ℛ, s′) : (nothing, nothing)
  if !isnothing(nsₗ)
    _newevent!(𝒬, _sweepline(s′), bundle, p, _segment(s′), _keyseg(nsₗ), digits)
  end
end

function _handletop!(bundle, ℛ, 𝒬, p, digits)
  # bundle is sorted sequence, so the last segment is maximum
  s″ = bundle[end]
  # element above s″
  _, nsᵤ = !isnothing(s″) ? BinaryTrees.prevnext(ℛ, s″) : (nothing, nothing)
  if !isnothing(nsᵤ)
    _newevent!(𝒬, _sweepline(s″), bundle, p, _segment(s″), _keyseg(nsᵤ), digits)
  end
end

# -----------------
# HELPER FUNCTIONS
# -----------------

# function to add new events to the event queue as needed
function _newevent!(𝒬, sweepline, bundle, p, seg₁, seg₂, digits)
  intersection(Segment(seg₁), Segment(seg₂)) do I
    t = type(I)
    if t === Crossing || t === EdgeTouching
      i = coordround(get(I), digits=digits)
      if i ≈ p # helps with vertical+horizontal intersections
        push!(bundle, SweepSegment(seg₂, sweepline))
        push!(bundle, SweepSegment(seg₁, sweepline))
      elseif i > p # add to 𝒬, update existing point if needed
        _addpoint!(𝒬, i, [seg₁, seg₂], 3)
      end
    end
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

# add an endpoint and corresponding segment to the event queue at a given position (1=start, 2=end, 3=crossing)
function _addpoint!(𝒬, p, segs, pos)
  U = Set{typeof(first(segs))} # set of segments
  node = BinaryTrees.search(𝒬, p)
  # updates or adds event point based on existing events in 𝒬
  if !isnothing(node)
    union!(BinaryTrees.value(node)[pos], segs)
  else
    # create a new event point with the value of (U(), U(), U()). pos determines whether inserted segments are starts or ends
    vals = ntuple(i -> i == pos ? U(segs) : U(), 3) # (Starts, Ends, Crossings)
    BinaryTrees.insert!(𝒬, p, vals)
  end
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

# handles  the degenerate case of trivial (0-length) segments
_istrivial(seg) = seg[1] == seg[2]

# ----------------
# DATA STRUCTURES
# ----------------

# tracks the event point and constructs the sweepline lazily
mutable struct SweepLine{P<:Point,ℒ<:Number}
  p::P
  ybounds::Tuple{ℒ,ℒ}
end
# getters
_sweeppoint(sweepline::SweepLine) = getfield(sweepline, :p)
_sweepx(sweepline::SweepLine) = CoordRefSystems.values(coords(_sweeppoint(sweepline)))[1]
_sweepy(sweepline::SweepLine) = CoordRefSystems.values(coords(_sweeppoint(sweepline)))[2]
_sweepbounds(sweepline::SweepLine) = getfield(sweepline, :ybounds)

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
function _sweepline(sweepline::SweepLine)
  x = _sweepx(sweepline)
  y = _sweepy(sweepline)
  # perturbation to avoid numerical issues
  ϵ = atol(x) + eps(x)
  lower, upper = _sweepbounds(sweepline)

  p₁ = Point(x + ϵ, lower) # lowest point
  p₂ = Point(x + ϵ, y + 2ϵ) # doubled to avoid precision issues
  p₃ = Point(x - ϵ, y + 2ϵ)
  p₄ = Point(x - ϵ, upper) # highest point
  # tuple of segments representing the sweepline
  # ((lower vertical), (horizontal), (upper vertical))
  ((p₁, p₂), (p₂, p₃), (p₃, p₄))
end

# sweepline intersection with segment
function _sweepintersect(seg::Tuple{P,P}, sweepline::SweepLine{P,ℒ}) where {P<:Point,ℒ<:Number}
  rope = _sweepline(sweepline)
  I = nothing
  # check for intersections between the segment and the sweepline.
  # `Rope()` intersections are not type stable; this loop is a workaround.
  for segᵣ in rope
    I = intersection(Segment(seg), Segment(segᵣ)) do I
      t = type(I)
      if t === Crossing || t === EdgeTouching || t === CornerTouching
        get(I)
      else
        nothing
      end
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

# constructor for SweepSegment using Sweepline
function SweepSegment(seg::Tuple{P,P}, sweepline::SweepLine{P,ℒ}) where {P<:Point,ℒ<:Number}
  y = _sweepintersect(seg, sweepline)
  SweepSegment{P,ℒ}(seg, sweepline, y, _sweeppoint(sweepline))
end

# getters for SweepSegment
_segment(s::SweepSegment) = getfield(s, :seg)
_yintersect(s::SweepSegment) = getfield(s, :yintersect)
_sweepline(s::SweepSegment) = getfield(s, :sweepline)
_sweeppoint(s::SweepSegment) = _sweeppoint(getfield(s, :sweepline))

Base.:(==)(s₁::SweepSegment, s₂::SweepSegment) = _segment(s₁) == _segment(s₂)

# compare two segments based on their sweepline intersection relative to the current event point.
function Base.isless(s₁::SweepSegment{P,ℒ}, s₂::SweepSegment{P,ℒ}) where {P<:Point,ℒ}
  # if segments same, return false
  if s₁ == s₂
    return false
  end
  #* calculating the y intersect is the largest performance bottleneck
  ya = _yintersect(s₁) # s₁ is always up-to-date
  yb = _ycalc!(s₂)

  diff = ustrip(abs(ya - yb))
  tol = eps(ℒ)
  # if segments are separated over y, check ya < yb
  if diff > tol
    ya < yb
  else
    # fallback to lexicographic ordering of segments
    seg₁, seg₂ = _segment(s₁), _segment(s₂)
    seg₁ < seg₂
  end
end

# calculate y-coordinate of intersection with sweepline
function _ycalc!(s::SweepSegment{P,ℒ}) where {P<:Point,ℒ<:Number}
  sweepline = _sweepline(s)
  # if the latest point is the sweepline point, use the precalculated intersection
  if s.latestpoint === _sweeppoint(sweepline)
    y = s.yintersect
  else
    # otherwise, calculate the intersection with the sweepline
    # and update
    y = convert(ℒ, _sweepintersect(_segment(s), sweepline))
    s.latestpoint = _sweeppoint(sweepline)
    s.yintersect = y
  end
  y
end
