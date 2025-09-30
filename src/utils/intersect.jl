# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intervalsweep(segments; [digits])

Compute pairwise intersections between n `segments`
with `digits` precision in O(n⋅log(n+k)) time using
an x-interval sweep line algorithm. Similar to an optimal
Bentley-Ottmann algorithm in sparse systems, and closer to O(n²) in dense systems.

By default, set `digits` based on the absolute
tolerance of the length type of the segments.

## References

* Bentley & Ottmann 1979. [Algorithms for reporting and counting
  geometric intersections](https://ieeexplore.ieee.org/document/1675432)
"""
function intervalsweep(segments; digits=_digits(segments))
  # orient segments and round coordinates
  segs = map(segments) do seg
    a, b = coordround.(extrema(seg), digits=digits)
    a > b ? (b, a) : (a, b)
  end
  𝒮 = SweepLineQueue(segs)

  points, seginds = handle!(𝒮; digits=digits)
end

# # compute the number of significant digits based on the segment type
# # this is used to determine the precision of the points
function _digits(segments)
  seg = first(segments)
  ℒ = lentype(seg)
  τ = ustrip(eps(ℒ))
  round(Int, 0.8 * (-log10(τ))) # 0.8 is a heuristic to avoid numerical issues
end

# ----------------
# DATA STRUCTURES
# ----------------

struct SweepLineInterval{T<:Number}
  start::T
  stop::T
  segment::Any
  index::Int
end

struct SweepLineQueue
  intervals::Vector{SweepLineInterval}
end
Base.length(𝒮::SweepLineQueue) = length(𝒮.intervals)
Base.getindex(𝒮::SweepLineQueue, i::Int) = 𝒮.intervals[i]

function SweepLineQueue(segs::Vector{<:Tuple{Point,Point}})
  intervals = Vector{SweepLineInterval}()
  for (i, seg) in enumerate(segs)
    x₁, _ = CoordRefSystems.values(coords(seg[1]))
    x₂, _ = CoordRefSystems.values(coords(seg[2]))
    push!(intervals, SweepLineInterval(min(x₁, x₂), max(x₁, x₂), seg, i))
  end
  SweepLineQueue(sort!(intervals, by=s -> s.start))
end

function overlaps(i₁::SweepLineInterval, i₂::SweepLineInterval)
  i₁.start ≤ i₂.start && i₂.stop ≥ i₁.start
end
# ----------------
# SWEEP LINE HANDLER
# ----------------

"""
  handle!(𝒮::SweepLineQueue)

Iterate through the sweep line queue and compute all intersection points
between overlapping intervals. Returns a tuple of intersection points and
the sets of segment indices that intersect at each point.
"""
function handle!(𝒮::SweepLineQueue; digits=10)
  𝐺 = Dict{Point,Set{Int}}()
  n = length(𝒮)
  for i in 1:n
    current = 𝒮[i]
    for k in (i + 1):n
      candidate = 𝒮[k]
      # If the intervals no longer overlap, break out of the inner loop
      if !overlaps(current, candidate)
        break
      end
      # Check if the segments actually intersect
      I = intersection(Segment(current.segment), Segment(candidate.segment)) do 𝑖
        t = type(𝑖)
        (t === Crossing || t === EdgeTouching) ? get(𝑖) : nothing
      end
      isnothing(I) || _addintersection!(𝐺, I, current.index, candidate.index; digits=digits)
    end
  end
  (collect(keys(𝐺)), collect(values(𝐺)))
end
function _addintersection!(𝐺::Dict{Point,Set{Int}}, I::Point, index₁::Int, index₂::Int; digits=10)
  p = coordround(I, digits=digits)
  if haskey(𝐺, p)
    union!(𝐺[p], index₁)
    union!(𝐺[p], index₂)
  else
    𝐺[p] = Set([index₁, index₂])
  end
end
