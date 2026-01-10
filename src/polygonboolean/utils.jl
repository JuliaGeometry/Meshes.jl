# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# insert intersection points into rings
function _insertintersections!(intersections, seginds, allrings::Vector{Vector{P}}) where {P}
  # group intersections by segment index
  G = Dict{Int,Vector{P}}()
  for (p, segs) in zip(intersections, seginds)
    for s in segs
      push!(get!(G, s, P[]), p)
    end
  end

  sortedseginds = sort(collect(keys(G)))

  # precompute offsets
  ℒ = length.(allrings)
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

    insertcounts[rind] += length(pts)
  end
end

# check if segment is filled above/below based on operation and fill bits
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
