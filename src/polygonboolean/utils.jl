# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# insert intersection points into rings
function _insertintersections!(intersections, seginds, allrings)
  # group intersections by segment index
  P = eltype(intersections)
  G = Dict{Int,Vector{P}}()
  for (p, segs) in zip(intersections, seginds)
    for s in segs
      push!(get!(G, s, P[]), p)
    end
  end

  sortedseginds = sort(collect(keys(G)))

  # precompute offsets
  ℒ = size.(allrings, 1)
  offsets = [0; cumsum(ℒ)]

  # track inserted points per ring
  insertcounts = zeros(Int, length(allrings))

  for ind in sortedseginds
    # find ring index
    rind = searchsortedfirst(offsets, ind) - 1
    lind = ind - offsets[rind]

    pts = G[ind]
    v = allrings[rind]

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

    insertcounts[rind] += size(pts, 1)
  end
end

##########
# CONSTANTS for segment fill bitmasks
# segment fills describe whether a segment is filled by a subject or clip polygon
# above or below the segment
##########
const NONE = 0b0000
const SUBJTOP = 0b0001
const SUBJBOTTOM = 0b0010
const CLIPTOP = 0b0100
const CLIPBOTTOM = 0b1000
const BOTHTOP = SUBJTOP | CLIPTOP
const BOTHBOTTOM = SUBJBOTTOM | CLIPBOTTOM

"""
  Segment Fill Rules for Boolean Operations

Internal utility functions for determining segment fill states during polygon boolean operations.

# Fill Rules by Operation

- **Union**: Segment is filled if either subject OR clip polygon is filled
- **Intersection**: Segment is filled if both subject AND clip polygons are filled
- **Difference**: Segment is filled if subject is filled AND clip is NOT filled
- **Symmetric Difference (XOR)**: Segment is filled if exactly one of subject or clip is filled

# Purpose

These rules establish segment-polygon relationships used to:
- Select relevant segments for the output
- Determine hole vs outer polygon relationships
- Build the final boolean operation result
"""
function _filled(operation, bits, above)
  masksubj = above ? SUBJTOP : SUBJBOTTOM
  maskclip = above ? CLIPTOP : CLIPBOTTOM

  # any subject/clip filled?
  subjfilled = (bits & masksubj) != 0
  clipfilled = (bits & maskclip) != 0

  if operation == intersect
    subjfilled && clipfilled
  elseif operation == union
    subjfilled || clipfilled
  elseif operation == setdiff
    subjfilled && !clipfilled
  elseif operation == symdiff || operation == xor
    subjfilled ⊻ clipfilled
  else
    operation(subjfilled, clipfilled)
  end
end
