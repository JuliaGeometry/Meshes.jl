# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
"""
    bentleyottmann(segments; [digits])

Compute pairwise intersections between n `segments`
with `digits` precision in O(n⋅log(n)) time using
Bentley-Ottmann sweep line algorithm.

By default, set `digits` based on the absolute
tolerance of the length type of the segments.

## References

* Bentley & Ottmann 1979. [Algorithms for reporting and counting
  geometric intersections](https://ieeexplore.ieee.org/document/1675432)
"""
function bentleyottmann(segments; digits=_digits(segments))
  TOL = 1 / 10^digits # precomputed tolerance for floating point comparisons
  
  # orient segments
  segs = map(segments) do s
    a, b = coordround.(extrema(s); digits)
    a > b ? (b, a) : (a, b)
  end

  # retrieve types
  T = lentype(first(segs)[1])
  P = typeof(first(segs)[1])
  S = Tuple{P,P}

  # initialization
  𝒬 = BinaryTrees.AVLTree{P}()
  ℛ = BinaryTrees.AVLTree{_SweepSegment{T,P}}()
  ℬ = Dict{P,Vector{S}}()
  ℰ = Dict{P,Vector{S}}()
  lookup = Dict{S,Int}()
  for (i, (a, b)) in enumerate(segs)
    BinaryTrees.insert!(𝒬, a)
    BinaryTrees.insert!(𝒬, b)
    haskey(ℬ, a) ? push!(ℬ[a], (a, b)) : (ℬ[a] = [(a, b)])
    haskey(ℰ, b) ? push!(ℰ[b], (a, b)) : (ℰ[b] = [(a, b)])
    lookup[(a, b)] = i
  end

  # sweep line algorithm
  points = Vector{P}()
  seginds = Vector{Vector{Int}}()
  sweepline = _initsweep(segs)
  while !BinaryTrees.isempty(𝒬)
    # current point (or event)
    p = BinaryTrees.key(BinaryTrees.minnode(𝒬))

    # delete point from event queue
    BinaryTrees.delete!(𝒬, p)
    # handle event, i.e. update 𝒬, ℛ and ℳ
    ℬₚ = get(ℬ, p, S[]) # segments with p at the begin
    ℰₚ = get(ℰ, p, S[]) # segments with p at the end
    P = typeof(p)
    ℳₚ = S[]
    _findintersections!(ℳₚ, ℛ, p, TOL) # segments with p at the middle
    activesegs = Set(ℬₚ ∪ ℳₚ)

    # report intersections
    if !isempty(activesegs) || !isempty(ℰₚ)
      inds = Set{Int}()
      for s in activesegs ∪ ℰₚ
        push!(inds, lookup[s])
      end
      push!(points, p)
      push!(seginds, collect(inds))
    end

    # handle status line
    _handlestatus!(ℛ, ℬₚ, ℳₚ, ℰₚ, sweepline, p, TOL)

    if isempty(activesegs)
      for s in ℰₚ
        sₗ, sᵣ = BinaryTrees.prevnext(ℛ, _SweepSegment(s, sweepline))
        isnothing(sₗ) || isnothing(sᵣ) || _newevent!(𝒬, p, _keyseg(sₗ), _keyseg(sᵣ), digits)
      end
    else
      BinaryTrees.isempty(ℛ) || _handlebottom!(activesegs, ℛ, sweepline, 𝒬, p, digits)

      BinaryTrees.isempty(ℛ) || _handletop!(activesegs, ℛ, sweepline, 𝒬, p, digits)
    end
  end

  (points, seginds)
end

##
## handling functions
##

function _handlestatus!(ℛ, ℬₚ, ℳₚ, ℰₚ, sweepline, p, TOL)
  for s in reverse(ℰₚ ∪ ℳₚ)
    segsweep = _SweepSegment(s, sweepline)
    isnothing(BinaryTrees.search(ℛ, segsweep)) || BinaryTrees.delete!(ℛ, segsweep)
  end

  sweepline.point = _nudge(p, TOL)

  for s in ℬₚ ∪ ℳₚ
    BinaryTrees.insert!(ℛ, _SweepSegment(s, sweepline))
  end
end

function _handlebottom!(activesegs, ℛ, sweepline, 𝒬, p, digits)
  s′ = BinaryTrees.key(_minsearch(ℛ, activesegs, sweepline))

  sₗ, _ = !isnothing(s′) ? BinaryTrees.prevnext(ℛ, s′) : (nothing, nothing)
  if !isnothing(sₗ)
    _newevent!(𝒬, p, _keyseg(sₗ), _segment(s′), digits)
  end
end

function _handletop!(activesegs, ℛ, sweepline, 𝒬, p, digits)
  s″ = BinaryTrees.key(_maxsearch(ℛ, activesegs, sweepline))

  _, sᵤ = !isnothing(s″) ? BinaryTrees.prevnext(ℛ, s″) : (nothing, nothing)
  if !isnothing(sᵤ)
    _newevent!(𝒬, p, _segment(s″), _keyseg(sᵤ), digits)
  end
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _newevent!(𝒬, p, s₁, s₂, digits)
  intersection(Segment(s₁), Segment(s₂)) do I
    if type(I) == Crossing || type(I) == EdgeTouching
      p′ = coordround(get(I); digits)
      if p′ ≥ p && isnothing(BinaryTrees.search(𝒬, p′))
        BinaryTrees.insert!(𝒬, p′)
      end
    end
  end
end

# find segments that intersect with the point p
function _findintersections!(ℳₚ, ℛ, p, TOL)
  x, y = CoordRefSystems.values(coords(p))
  tol = TOL * unit(x) # ensure TOL is in the same unit as x and y
  _search!(BinaryTrees.root(ℛ), ℳₚ, x, y, tol)
  ℳₚ
end

function _search!(node, ℳₚ, x, y, TOL)
  isnothing(node) && return
  seg = _segment(BinaryTrees.key(node))
  x₁, y₁ = CoordRefSystems.values(coords(seg[1]))
  x₂, y₂ = CoordRefSystems.values(coords(seg[2]))

  # Precompute reused values
  dx, dy = x₂ - x₁, y₂ - y₁
  ℒ = hypot(dx, dy) # handling precision issues

  # Ensure the point is not the endpoint (avoids duplicates)
  skip = (x₂ - TOL ≤ x ≤ x₂ + TOL) && (y₂ - TOL ≤ y ≤ y₂ + TOL)
  # if collinear and not an endpoint
  collinear = dy * (x - x₁) - dx * (y - y₁)
  boundcheck = (x₁ - TOL ≤ x ≤ x₂ + TOL) && (y₁ - TOL ≤ y ≤ y₂ + TOL)

  if abs(collinear) ≤ TOL * ℒ #&& boundcheck
    skip || push!(ℳₚ, seg)
  end

  # using difference in y to determine the side of the segment
  # needed to avoid recursion depth and floating point issues
  ŷ = y₁ + (y₂ - y₁) * (x - x₁) / dx # y coordinate of the segment at x

  diff = y - ŷ # difference between the point and the segment

  if diff < -TOL
    _search!(BinaryTrees.left(node), ℳₚ, x, y, TOL)
  elseif diff > TOL
    _search!(BinaryTrees.right(node), ℳₚ, x, y, TOL)
  else
    # if the point is on the segment, check both sides
    _search!(BinaryTrees.left(node), ℳₚ, x, y, TOL)
    _search!(BinaryTrees.right(node), ℳₚ, x, y, TOL)
  end
end

# find the minimum segment among active segments in tree
function _minsearch(ℛ, activesegs, sweepline)
  activeordered = sort([_SweepSegment(s, sweepline) for s in activesegs])
  BinaryTrees.search(ℛ, activeordered[begin])
end

# find the maximum segment among active segments in tree
function _maxsearch(ℛ, activesegs, sweepline)
  activeordered = sort([_SweepSegment(s, sweepline) for s in activesegs])
  BinaryTrees.search(ℛ, activeordered[end])
end

# compute rounding digits
function _digits(segments)
  s = first(segments)
  ℒ = lentype(s)
  τ = ustrip(atol(ℒ))
  round(Int, -log10(τ)) - 1
end

# convenience function to get the segment from the node
function _keyseg(segment)
  _segment(BinaryTrees.key(segment))
end
# nudge the sweepline to get correct ℛ ordering
function _nudge(p, TOL)
  x, y = CoordRefSystems.values(coords(p))
  nudgefactor = unit(x) * TOL * 2
  Point(x + nudgefactor, y + nudgefactor)
end

# ----------------
# DATA STRUCTURES
# ----------------

# tracks sweepline and current y position for searching
mutable struct _SweepLine{P<:Point}
  point::P
end
_sweeppoint(sweepline::_SweepLine) = getfield(sweepline, :point)
_sweepx(sweepline::_SweepLine) = CoordRefSystems.values(coords(_sweeppoint(sweepline)))[1]
_sweepy(sweepline::_SweepLine) = CoordRefSystems.values(coords(_sweeppoint(sweepline)))[2]

# compute the intersection of a segment with the sweepline
function _sweepintersect(seg, sweepline)
  p₁, p₂ = coords.(seg)
  x₁, y₁ = CoordRefSystems.values(p₁)
  x₂, y₂ = CoordRefSystems.values(p₂)
  T = eltype(x₁)

  x = T(_sweepx(sweepline))
  y = T(_sweepy(sweepline))

  if abs(x₁ - x₂) < atol(eltype(x₁))
    return T(min(y, y₂)) # vertical goes at end
  end

  t = (x - x₁) / (x₂ - x₁)
  T(y₁ + t * (y₂ - y₁))
end

# takes input segment and assigns where it intersects sweepline
mutable struct _SweepSegment{T,P<:Point}
  seg::Tuple{P,P}
  xintersect::T
end

# constructor for _SweepSegment
function _SweepSegment(seg::Tuple{P,P}, sweepline::_SweepLine) where {P<:Point}
  xintersect = _sweepintersect(seg, sweepline)
  T = typeof(xintersect)
  _SweepSegment{T,P}(seg, xintersect)
end

_segment(sweepsegment::_SweepSegment) = sweepsegment.seg
_xintersect(sweepsegment::_SweepSegment) = getfield(sweepsegment, :xintersect)

Base.isless(a::_SweepSegment, b::_SweepSegment) = begin
  # if segments same, return false
  if _segment(a) == _segment(b)
    return false
  end
  ya, yb = _xintersect(a), _xintersect(b)

  # if segments are separated over y, check ya < yb
  if abs(ya - yb) > atol(eltype(ya))
    return ya < yb
  end

  # fallback to x-coordinates
  _segment(a) != _segment(b) && _segment(a) < _segment(b)
end

# initialize the sweep line with them minimum
function _initsweep(segs)
  U = lentype(coords(segs[1][1]))
  _SweepLine(Point(U(-Inf), U(-Inf)))
end
