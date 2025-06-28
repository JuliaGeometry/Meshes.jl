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
  refₚ = first(segments)
  T = lentype(refₚ)
  # get the coordinate reference system and manifold for stable point type
  M = manifold(refₚ)
  C = crs(refₚ)
  P = Point{M,C}
  S = Tuple{P,P}
  U = Set{S}

  # orient segments
  segs = map(segments) do s
    a, b = coordround.(extrema(s), digits=digits)
    a > b ? (b, a) : (a, b)
  end

  pmin, pmax = segs |> boundingbox |> Stretch(1.05) |> extrema
  _, ymin = CoordRefSystems.values(coords(pmin))
  _, ymax = CoordRefSystems.values(coords(pmax))
  ybounds = (ymin, ymax)

  # retrieve types
  # initialization
  𝒬 = BinaryTrees.AVLTree{P,Tuple{U,U,U}}()
  ℛ = BinaryTrees.AVLTree{_SweepSegment{P,T}}()
  lookup = Dict{S,Int}()
  # loop through segments and insert them into the event queue 𝒬 with start and ending flags for the segments
  # additionally build the lookup table for segment indices
  for (i, (a, b)) in enumerate(segs)
    # add starting point and segment
    aν = BinaryTrees.search(𝒬, a)
    if !isnothing(aν)
      union!(BinaryTrees.value(aν)[1], U([(a, b)]))
    else
      BinaryTrees.insert!(𝒬, a, (U([(a, b)]), U(), U()))
    end

    # add ending point and segment
    bν = BinaryTrees.search(𝒬, b)
    if !isnothing(bν)
      union!(BinaryTrees.value(bν)[2], U([(a, b)]))
    else
      BinaryTrees.insert!(𝒬, b, (U(), U([(a, b)]), U()))
    end

    # lookup table for segment indices
    lookup[(a, b)] = i
  end

  # initialize sweepline
  sweepline = _SweepLine{P,T}(pmin, ybounds)
  output = Dict{P,Vector{Int}}()
  # sweep line algorithm
  while !BinaryTrees.isempty(𝒬)
    # current point (or event)
    node = BinaryTrees.minnode(𝒬)
    p = BinaryTrees.key(node)
    # delete point from event queue
    BinaryTrees.delete!(𝒬, p)

    # beginning, ending, and crossing segments
    ℬ, ℰ, ℳ = BinaryTrees.value(node)

    # crosses that arent endpoints (including them can lead to duplicates)
    ℳₚ = setdiff(ℳ, ℰ)

    # handle status line
    _handlestatus!(ℛ, ℬ, ℳₚ, ℰ, sweepline, p)

    # build sorted bundle of segments
    bundle = @isdefined(bundle) ? (empty!(bundle); bundle) : Vector{_SweepSegment{P,T}}()
    for s in Iterators.flatten((ℬ, ℳₚ))
      push!(bundle, _SweepSegment(s, sweepline))
    end
    sort!(bundle)

    # process bundled events
    if isempty(bundle)
      # if endpoint, check new adjacent segments
      sₗ, sᵣ = BinaryTrees.prevnext(ℛ, _SweepSegment(first(ℰ), sweepline))
      isnothing(sₗ) || isnothing(sᵣ) || _newevent!(𝒬, sweepline, bundle, p, _keyseg(sₗ), _keyseg(sᵣ), digits)
    else

      # handle bottom and top events
      BinaryTrees.isempty(ℛ) || _handlebottom!(bundle, ℛ, 𝒬, p, digits)

      BinaryTrees.isempty(ℛ) || _handletop!(bundle, ℛ, 𝒬, p, digits)
    end

    # add necessary points and segment indices to output
    if !isempty(bundle) || !isempty(ℰ)
      inds = Set{Int}()
      # for s in bundle# ∪ ℰₚ
      for s in bundle
        push!(inds, lookup[_segment(s)])
      end
      for s in ℰ
        push!(inds, lookup[s])
      end
      indᵥ = collect(inds)
      if haskey(output, p)
        union!(output[p], indᵥ)
      else
        output[p] = indᵥ
      end
    end
  end

  (collect(keys(output)), collect(values(output)))
end

##
## handling functions
##

function _handlestatus!(ℛ, ℬₚ, ℳₚ, ℰₚ, sweepline, p)
  for s in ℰₚ ∪ ℳₚ
    segsweep = _SweepSegment(s, sweepline)
    BinaryTrees.delete!(ℛ, segsweep)
  end

  # update sweepline
  sweepline.point = p

  for s in ℬₚ ∪ ℳₚ
    BinaryTrees.insert!(ℛ, _SweepSegment(s, sweepline))
  end
end

function _handlebottom!(bundle, ℛ, 𝒬, p, digits)
  s′ = bundle[begin]

  sₗ, _ = !isnothing(s′) ? BinaryTrees.prevnext(ℛ, s′) : (nothing, nothing)
  if !isnothing(sₗ)
    _newevent!(𝒬, _sweepline(s′), bundle, p, _segment(s′), _keyseg(sₗ), digits)
  end
end

function _handletop!(bundle, ℛ, 𝒬, p, digits)
  s″ = bundle[end]

  _, sᵤ = !isnothing(s″) ? BinaryTrees.prevnext(ℛ, s″) : (nothing, nothing)
  if !isnothing(sᵤ)
    _newevent!(𝒬, _sweepline(s″), bundle, p, _segment(s″), _keyseg(sᵤ), digits)
  end
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _newevent!(𝒬, sweepline, bundle, p, s₁, s₂, digits)
  intersection(Segment(s₁), Segment(s₂)) do I
    t = type(I)
    if t === Crossing || t === EdgeTouching
      i = coordround(get(I), digits=digits)
      if i ≥ p
        if i ≈ p # if intersection point is the same as the sweepline point
          # this avoids repeatedly inserting the same intersection point
          push!(bundle, _SweepSegment(s₂, sweepline))
          push!(bundle, _SweepSegment(s₁, sweepline))
        else # if intersection point is different from the sweepline point
          node = BinaryTrees.search(𝒬, i)

          if isnothing(node) # insert new event into the event queue
            S = typeof(s₁)
            U = Set{S}
            ν = (U(), U(), U([s₁, s₂]))
            BinaryTrees.insert!(𝒬, i, ν)
          else # union with existing event
            union!(BinaryTrees.value(node)[3], [s₁, s₂])
          end
        end
      end
    end
  end
end

# compute rounding digits for FP precision # current based on reasonable FP64
function _digits(segments)
  s = first(segments)
  ℒ = lentype(s)
  τ = ustrip(eps(ℒ))
  round(Int, 0.8 * (-log10(τ)))
end

# convenience function to get the segment from the node
_keyseg(segment) = _segment(BinaryTrees.key(segment))

# ----------------
# DATA STRUCTURES
# ----------------

# tracks sweepline and current y position for searching
mutable struct _SweepLine{P<:Point,T<:Number}
  point::P
  ybounds::Tuple{T,T}
end
_sweeppoint(sweepline::_SweepLine) = getfield(sweepline, :point)
_sweepx(sweepline::_SweepLine) = CoordRefSystems.values(coords(_sweeppoint(sweepline)))[1]
_sweepy(sweepline::_SweepLine) = CoordRefSystems.values(coords(_sweeppoint(sweepline)))[2]
_sweepbounds(sweepline::_SweepLine) = getfield(sweepline, :ybounds)

#= sweepline definition
This odd definition is used to
handle intersections elegantly. It ensures that y coordinates
of non vertical segments are always correct unless overlapping,
vertical segments are always on top of the segments next to p,
ending segments don't intersect,
and all other segments are correctly ordered.

Inspired by LEDA implementation, but modified for Julia
=#
function _sweepline(sweepline::_SweepLine)
  x = _sweepx(sweepline)
  y = _sweepy(sweepline)
  ϵ = atol(x) + eps(x)
  lower, upper = _sweepbounds(sweepline)

  p₁ = Point(x + ϵ, lower)
  p₂ = Point(x + ϵ, y + 2ϵ)
  p₃ = Point(x - ϵ, y + 2ϵ)
  p₄ = Point(x - ϵ, upper)

  ((p₁, p₂), (p₂, p₃), (p₃, p₄))
end

# sweepline intersection with segment
function _sweepintersect(seg::Tuple{P,P}, sweepline::_SweepLine{P,T}) where {P<:Point,T<:Number}
  rope = _sweepline(sweepline)
  I = nothing
  # this loop checks for intersections between the segment and the sweepline
  # Rope intersections are not type stable, so this is a simple workaround
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
  # happens if ending segment
  i = isnothing(I) ? seg[2] : I
  _, y = CoordRefSystems.values(coords(i))
  y
end

# takes input segment and attaches where it intersects sweepline
mutable struct _SweepSegment{P<:Point,T<:Number}
  const seg::Tuple{P,P}
  const sweepline::_SweepLine{P,T}
  xintersect::T #information about the intersection with the sweepline
  latestpoint::P # latest point of the sweepline used to calculate the intersection
end

# constructor for _SweepSegment
function _SweepSegment(seg::Tuple{P,P}, sweepline::_SweepLine{P,T}) where {P<:Point,T<:Number}
  y = _sweepintersect(seg, sweepline)
  _SweepSegment{P,T}(seg, sweepline, y, _sweeppoint(sweepline))
end

_segment(sweepsegment::_SweepSegment) = getfield(sweepsegment, :seg)
_xintersect(sweepsegment::_SweepSegment) = getfield(sweepsegment, :xintersect)
_sweepline(sweepsegment::_SweepSegment) = getfield(sweepsegment, :sweepline)
_sweeppoint(sweepsegment::_SweepSegment) = _sweeppoint(getfield(sweepsegment, :sweepline))

Base.:(==)(a::_SweepSegment, b::_SweepSegment) = _segment(a) == _segment(b)

function Base.isless(a::_SweepSegment{P,T}, b::_SweepSegment{P,T}) where {P<:Point,T}
  # if segments same, return false
  segₐ, segᵦ = _segment(a), _segment(b)
  if isequal(segₐ, segᵦ)
    return false
  end
  # update to latest sweepline reference
  p = _sweeppoint(a)
  a₁, _ = CoordRefSystems.values.(coords.(segₐ))
  b₁, b₂ = CoordRefSystems.values.(coords.(segᵦ))

  if p == a₁
    s = orient2(b₁, b₂, CoordRefSystems.values(coords(p)))
  else
    s = 0
  end

  if s != 0 || _istrivial(segₐ) || _istrivial(segᵦ)
    return s < 0
  end

  #* This is the largest performance bottleneck
  # compute intersection over sweepline
  ya = _xintersect(a)
  yb = _ycalc!(b)

  diff = ustrip(abs(ya - yb))
  tol = eps(T)
  # if segments are separated over y, check ya < yb
  if diff > tol
    ya < yb
  else
    # fallback to lexicographic ordering of segments
    segₐ < segᵦ
  end
end

# calculate y-coordinate of intersection with sweepline
function _ycalc!(a::_SweepSegment{P,T}) where {P<:Point,T<:Number}
  sweepline = _sweepline(a)
  if a.latestpoint === _sweeppoint(sweepline)
    # if the latest point is the sweepline point, use the intersect
    y = a.xintersect
  else
    # otherwise, calculate the intersection with the sweepline
    # and update the latest point
    y = convert(T, _sweepintersect(_segment(a), sweepline))

    a.latestpoint = _sweeppoint(sweepline)
    a.xintersect = y
  end
  y
end

# handles  the degenerate case of segments that are trivial
_istrivial(s) = s[1] == s[2]
