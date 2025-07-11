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

### Note

FP32 will likely be incorrect for precision-sensitive tasks. Rounding will help.
"""
function bentleyottmann(segments; digits=_digits(segments))
  # orient segments and round coordinates
  segs = map(segments) do s
    a, b = coordround.(extrema(s), digits=digits)
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
  ℛ = BinaryTrees.AVLTree{_SweepSegment{P,ℒ}}()

  # lookup table for segment indices
  lookup = Dict{S,Int}()

  # initialize event queue and lookup table
  for (i, seg) in enumerate(segs)
    a, b = seg
    _addstartpoint!(𝒬, a, b, U)
    _addendpoint!(𝒬, a, b, U)
    lookup[seg] = i
  end

  # initialize sweepline
  pmin, _ = extrema(first.(segs))
  ybounds = _ybounds(segs)
  sweepline = _SweepLine{P,ℒ}(pmin, ybounds)

  # output dictionary (planar graph 𝐺)
  𝐺 = Dict{P,Vector{Int}}()

  # vector holding segments intersecting the current event point
  bundle = Vector{_SweepSegment{P,ℒ}}()

  # sweep line algorithm
  while !BinaryTrees.isempty(𝒬)
    # current event point
    n = BinaryTrees.minnode(𝒬)
    p = BinaryTrees.key(n)
    BinaryTrees.delete!(𝒬, p)

    # sets of beginning, ending, and crossing segments that include the current event point
    ℬ, ℰ, ℳ = BinaryTrees.value(n)

    # crosses that aren't endpoints (including them can lead to duplicates)
    ℳₚ = setdiff(ℳ, ℰ)

    # handle status structure
    _handlestatus!(ℛ, ℬ, ℳₚ, ℰ, sweepline, p)

    # build bundle of segments crossing the current event point
    bundle = empty!(bundle)
    for s in Iterators.flatten((ℬ, ℳₚ))
      push!(bundle, _SweepSegment(s, sweepline))
    end
    sort!(bundle)

    # process bundled events
    if isempty(bundle) # occurs at endpoints
      # check newly adjacent segments
      sₗ, sᵣ = BinaryTrees.prevnext(ℛ, _SweepSegment(first(ℰ), sweepline))
      isnothing(sₗ) || isnothing(sᵣ) || _newevent!(𝒬, sweepline, bundle, p, _keyseg(sₗ), _keyseg(sᵣ), digits)
    else
      # handle bottom and top events
      BinaryTrees.isempty(ℛ) || _handlebottom!(bundle, ℛ, 𝒬, p, digits)
      BinaryTrees.isempty(ℛ) || _handletop!(bundle, ℛ, 𝒬, p, digits)
    end

    # add intersection points and corresponding segment indices to the output
    if length(bundle) > length(ℬ) # bundle only has ℬ unless p is an intersection
      inds = Set{Int}()

      # start and crossing segments
      for s in bundle
        push!(inds, lookup[_segment(s)])
      end

      # ending segments
      for s in ℰ
        push!(inds, lookup[s])
      end

      # add indices to output
      indᵥ = collect(inds)
      if haskey(𝐺, p)
        union!(𝐺[p], indᵥ)
      else
        𝐺[p] = indᵥ
      end
    end
  end

  (collect(keys(𝐺)), collect(values(𝐺)))
end

# ------------------------------------
# Sweep line status and event handling
# ------------------------------------

function _handlestatus!(ℛ, ℬₚ, ℳₚ, ℰₚ, sweepline, p)
  # remove segments that are no longer active or need to be updated
  for s in ℰₚ ∪ ℳₚ
    sweepseg = _SweepSegment(s, sweepline)
    BinaryTrees.delete!(ℛ, sweepseg)
  end

  # update sweepline
  sweepline.point = p

  # insert segments into the status structure
  for s in ℬₚ ∪ ℳₚ
    BinaryTrees.insert!(ℛ, _SweepSegment(s, sweepline))
  end
end

function _handlebottom!(bundle, ℛ, 𝒬, p, digits)
  # bundle is sorted sequence, so the first segment is minimum
  s′ = bundle[begin]

  sₗ, _ = !isnothing(s′) ? BinaryTrees.prevnext(ℛ, s′) : (nothing, nothing)
  if !isnothing(sₗ)
    _newevent!(𝒬, _sweepline(s′), bundle, p, _segment(s′), _keyseg(sₗ), digits)
  end
end

function _handletop!(bundle, ℛ, 𝒬, p, digits)
  # bundle is sorted sequence, so the last segment is maximum
  s″ = bundle[end]

  _, sᵤ = !isnothing(s″) ? BinaryTrees.prevnext(ℛ, s″) : (nothing, nothing)
  if !isnothing(sᵤ)
    _newevent!(𝒬, _sweepline(s″), bundle, p, _segment(s″), _keyseg(sᵤ), digits)
  end
end

# -----------------
# HELPER FUNCTIONS
# -----------------

# Function to add new events to the event queue as needed
function _newevent!(𝒬, sweepline, bundle, p, s₁, s₂, digits)
  intersection(Segment(s₁), Segment(s₂)) do I
    t = type(I)
    # Only process if the intersection is a Crossing or EdgeTouching
    if t === Crossing || t === EdgeTouching
      i = coordround(get(I), digits=digits)

      # Only consider intersection points at or above/to the right of the current event point
      if i ≥ p
        if i ≈ p
          # If the intersection coincides with the current event point,
          # add both segments to the bundle to avoid duplicate events
          push!(bundle, _SweepSegment(s₂, sweepline))
          push!(bundle, _SweepSegment(s₁, sweepline))
        else
          # Otherwise, insert or update the event queue with the new intersection
          n = BinaryTrees.search(𝒬, i)
          if isnothing(n)
            S = typeof(s₁) # Segment type
            U = Set{S} # Set of segments
            ν = (U(), U(), U([s₁, s₂])) # (start, end, crossing segments)
            BinaryTrees.insert!(𝒬, i, ν)
          else
            # If the node already exists, update the crossing segments
            union!(BinaryTrees.value(n)[3], [s₁, s₂])
          end
        end
      end
    end
  end
end

# Compute the number of significant digits based on the segment type
# This is used to determine the precision of the points
function _digits(segments)
  s = first(segments)
  ℒ = lentype(s)
  τ = ustrip(eps(ℒ))
  round(Int, 0.8 * (-log10(τ)))
end

# Initialize starting and ending points in the event queue

# Ensure each initial event point contains all needed segments
# updates existing events if needed

# Add a segment to the event queue at a given point and position (1=start, 2=end)
function _addinitpoint!(𝒬, p, s, U, pos)
  node = BinaryTrees.search(𝒬, p)
  if !isnothing(node)
    union!(BinaryTrees.value(node)[pos], U([s]))
  else
    vals = ntuple(i -> i == pos ? U([s]) : U(), 3) # (Starts, Ends, Crossings)
    BinaryTrees.insert!(𝒬, p, vals)
  end
end

# Add starting point and segment
_addstartpoint!(𝒬, a, b, U) = _addinitpoint!(𝒬, a, (a, b), U, 1)

# Add ending point and segment
_addendpoint!(𝒬, a, b, U) = _addinitpoint!(𝒬, b, (a, b), U, 2)

# Compute y bounds of the segment domain
function _ybounds(segs)
  T = numtype(lentype(first(first(segs))))
  pmin, pmax = segs |> boundingbox |> Stretch(T(1.05)) |> extrema
  _, ymin = CoordRefSystems.values(coords(pmin))
  _, ymax = CoordRefSystems.values(coords(pmax))
  (ymin, ymax)
end

# Convenience function to get the segment from the node
_keyseg(segment) = _segment(BinaryTrees.key(segment))

# Handles  the degenerate case of trivial (0-length) segments
_istrivial(s) = s[1] == s[2]

# ----------------
# DATA STRUCTURES
# ----------------

# Tracks the event point and constructs the sweepline lazily
mutable struct _SweepLine{P<:Point,ℒ<:Number}
  point::P
  ybounds::Tuple{ℒ,ℒ}
end
# Getters
_sweeppoint(sweepline::_SweepLine) = getfield(sweepline, :point)
_sweepx(sweepline::_SweepLine) = CoordRefSystems.values(coords(_sweeppoint(sweepline)))[1]
_sweepy(sweepline::_SweepLine) = CoordRefSystems.values(coords(_sweeppoint(sweepline)))[2]
_sweepbounds(sweepline::_SweepLine) = getfield(sweepline, :ybounds)

#= Sweepline Definition ---

This odd definition is used to
handle intersections elegantly. It ensures that y coordinates
of non vertical segments are always correct,
vertical segments are always on top of non-vertical segments next to p,
ending segments don't intersect, and all other segments are correctly ordered.

Inspired by LEDA implementation, but modified for Julia
=#
function _sweepline(sweepline::_SweepLine)
  x = _sweepx(sweepline)
  y = _sweepy(sweepline)

  # Perturbation to avoid numerical issues
  ϵ = atol(x) + eps(x)
  lower, upper = _sweepbounds(sweepline)

  p₁ = Point(x + ϵ, lower)
  p₂ = Point(x + ϵ, y + 2ϵ) # doubled to avoid inaccurate ycalc from precision
  p₃ = Point(x - ϵ, y + 2ϵ)
  p₄ = Point(x - ϵ, upper)

  # Tuple of segments representing the sweepline
  ((p₁, p₂), (p₂, p₃), (p₃, p₄))
end

# Sweepline intersection with segment
function _sweepintersect(seg::Tuple{P,P}, sweepline::_SweepLine{P,ℒ}) where {P<:Point,ℒ<:Number}
  rope = _sweepline(sweepline)
  I = nothing
  # Check for intersections between the segment and the sweepline.
  # `Rope()` intersections are not type stable; this loop is a workaround.
  for piece in rope
    I = intersection(Segment(seg), Segment(piece)) do I
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

  # If I isn't set, then we are handling an ending segment, so seg[2] is used
  i = isnothing(I) ? seg[2] : I
  _, y = CoordRefSystems.values(coords(i))
  y
end

# Attaches the intersection point of the input segment with the sweepline
mutable struct _SweepSegment{P<:Point,ℒ<:Number}
  const seg::Tuple{P,P}
  const sweepline::_SweepLine{P,ℒ}
  yintersect::ℒ #information about the intersection with the sweepline
  latestpoint::P # latest point of the sweepline used to calculate the intersection
end

# constructor for _SweepSegment using Sweepline
function _SweepSegment(seg::Tuple{P,P}, sweepline::_SweepLine{P,ℒ}) where {P<:Point,ℒ<:Number}
  y = _sweepintersect(seg, sweepline)
  _SweepSegment{P,ℒ}(seg, sweepline, y, _sweeppoint(sweepline))
end

# getters for _SweepSegment
_segment(sweepsegment::_SweepSegment) = getfield(sweepsegment, :seg)
_yintersect(sweepsegment::_SweepSegment) = getfield(sweepsegment, :yintersect)
_sweepline(sweepsegment::_SweepSegment) = getfield(sweepsegment, :sweepline)
_sweeppoint(sweepsegment::_SweepSegment) = _sweeppoint(getfield(sweepsegment, :sweepline))

Base.:(==)(a::_SweepSegment, b::_SweepSegment) = _segment(a) == _segment(b)

# Compare two segments based on their sweepline intersection relative to the current event point.
function Base.isless(a::_SweepSegment{P,ℒ}, b::_SweepSegment{P,ℒ}) where {P<:Point,ℒ}
  # if segments same, return false
  segₐ, segᵦ = _segment(a), _segment(b)
  if segₐ == segᵦ
    return false
  end

  # Retrieve sweepline point and segment coordinates
  p = _sweeppoint(a)
  a₁, _ = CoordRefSystems.values.(coords.(segₐ))
  b₁, b₂ = CoordRefSystems.values.(coords.(segᵦ))

  # If start of a is on the event point, use orientation to determine order relative to b
  if p == a₁
    s = sign(signarea(b₁, b₂, CoordRefSystems.values(coords(p))))
    # AdaptivePredicates.jl may be more robust (maybe needed in the future)
  else
    s = 0
  end

  if s != 0 || _istrivial(segₐ) || _istrivial(segᵦ)
    return s < 0
  end

  #* Calculating the y intersect is the largest performance bottleneck
  # compute intersection over sweepline
  ya = _yintersect(a) # a is always up-to-date
  yb = _ycalc!(b)

  diff = ustrip(abs(ya - yb))
  tol = eps(ℒ)

  # if segments are separated over y, check ya < yb
  if diff > tol
    ya < yb
  else
    # fallback to lexicographic ordering of segments
    segₐ < segᵦ
  end
end

# calculate y-coordinate of intersection with sweepline
function _ycalc!(a::_SweepSegment{P,ℒ}) where {P<:Point,ℒ<:Number}
  sweepline = _sweepline(a)

  # if the latest point is the sweepline point, use the precalculated intersection
  if a.latestpoint === _sweeppoint(sweepline)
    y = a.yintersect
  else
    # otherwise, calculate the intersection with the sweepline
    # and update
    y = convert(ℒ, _sweepintersect(_segment(a), sweepline))

    a.latestpoint = _sweeppoint(sweepline)
    a.yintersect = y
  end
  y
end
