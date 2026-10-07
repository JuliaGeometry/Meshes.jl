# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GreinerHormannClipping()

The Greiner-Hormann algorithm for clipping polygons, with the extension
of Foster et al. for degenerate intersections.

## References

* Greiner, G. & Hormann, K. 1998. [Efficient clipping of arbitrary
  polygons](https://dl.acm.org/doi/pdf/10.1145/274363.274364)

* Foster, E. L., Hormann, K. & Popa, R. T. 2019. [Clipping simple polygons
  with degenerate intersections](https://doi.org/10.1016/j.cagx.2019.100007)

### Notes

In contrast to [`SutherlandHodgmanClipping`](@ref), the clipping geometry
can be non-convex, and vertices of one geometry can lie on edges of the
other geometry.
"""
struct GreinerHormannClipping <: ClippingMethod end

function clip(subject::Polygon, other::Polygon, ::GreinerHormannClipping)
  srings = subject |> Repair(11) |> rings
  orings = other |> Repair(11) |> rings

  # intersection phase
  slist, olist = _ghintersect(srings, orings)

  # labeling phase
  _ghmark!(slist, olist)
  _ghflags!(slist, orings)
  _ghflags!(olist, srings)

  # tracing phase
  crings = _ghtrace(slist, olist)
  append!(crings, _ghinner(slist, orings, common=true))
  append!(crings, _ghinner(olist, srings, common=false))

  isempty(crings) ? nothing : _ghpolygons(crings)
end

# vertex of the doubly-linked lists of the algorithm
@kwdef mutable struct GHVertex{P<:Point}
  point::P
  intersect::Bool
  neighbor::Int = 0
  crossing::Bool = false
  entry::Bool = false
  visited::Bool = false
end

# doubly-linked list of vertices split into components (rings)
struct GHList{P<:Point}
  verts::Vector{GHVertex{P}}
  comps::Vector{UnitRange{Int}}
  compof::Vector{Int}
end

# -------------------
# INTERSECTION PHASE
# -------------------

# insert the intersections of the two sets of rings as vertices in both lists
function _ghintersect(rings₁, rings₂)
  vs₁ = [vertices(r) for r in rings₁]
  vs₂ = [vertices(r) for r in rings₂]

  # point type and alpha/beta type
  P = promote_type(eltype(first(vs₁)), eltype(first(vs₂)))
  T = numtype(lentype(P))

  # events attached to original vertices, and inserted along edges
  tag₁ = [zeros(Int, length(v)) for v in vs₁]
  tag₂ = [zeros(Int, length(v)) for v in vs₂]
  ins₁ = [[Tuple{T,P,Int}[] for _ in v] for v in vs₁]
  ins₂ = [[Tuple{T,P,Int}[] for _ in v] for v in vs₂]

  nevents = 0
  for r₁ in eachindex(vs₁), i in eachindex(vs₁[r₁])
    a₁ = vs₁[r₁][i]
    b₁ = vs₁[r₁][i + 1]
    for r₂ in eachindex(vs₂), j in eachindex(vs₂[r₂])
      a₂ = vs₂[r₂][j]
      b₂ = vs₂[r₂][j + 1]

      sa₁ = signarea(a₁, a₂, b₂)
      sb₁ = signarea(b₁, a₂, b₂)
      sa₂ = signarea(a₂, a₁, b₁)
      sb₂ = signarea(b₂, a₁, b₁)

      if !isapproxzero(sa₁ - sb₁) && !isapproxzero(sa₂ - sb₂) # not parallel
        α = sa₁ / (sa₁ - sb₁)
        β = sa₂ / (sa₂ - sb₂)
        (_ghinunit(α) && _ghinunit(β)) || continue
        αzero = isapproxzero(α)
        βzero = isapproxzero(β)
        if !αzero && !βzero
          # X-intersection: new vertex on both edges
          nevents += 1
          _ghins!(ins₁[r₁], i, (α, a₁ + α * (b₁ - a₁), nevents))
          _ghins!(ins₂[r₂], j, (β, a₂ + β * (b₂ - a₂), nevents))
        elseif !αzero
          # T-intersection: vertex of the second ring on an edge of the first
          nevents += 1
          _ghins!(ins₁[r₁], i, (α, a₂, nevents))
          _ghtag!(tag₂[r₂], j, nevents)
        elseif !βzero
          # T-intersection: vertex of the first ring on an edge of the second
          nevents += 1
          _ghins!(ins₂[r₂], j, (β, a₁, nevents))
          _ghtag!(tag₁[r₁], i, nevents)
        else
          # V-intersection: coincident vertices
          nevents += 1
          _ghtag!(tag₁[r₁], i, nevents)
          _ghtag!(tag₂[r₂], j, nevents)
        end
      elseif isapproxzero(sa₁) && isapproxzero(sb₁) && isapproxzero(sa₂) && isapproxzero(sb₂)
        # collinear edges, possibly overlapping
        u₁ = b₁ - a₁
        u₂ = b₂ - a₂
        α = (a₂ - a₁) ⋅ u₁ / (u₁ ⋅ u₁)
        β = (a₁ - a₂) ⋅ u₂ / (u₂ ⋅ u₂)
        αinside = _ghinunit(α) && !isapproxzero(α)
        βinside = _ghinunit(β) && !isapproxzero(β)
        if αinside && βinside
          # X-overlap: both vertices lie inside the other edge
          nevents += 1
          _ghins!(ins₁[r₁], i, (α, a₂, nevents))
          _ghtag!(tag₂[r₂], j, nevents)
          nevents += 1
          _ghins!(ins₂[r₂], j, (β, a₁, nevents))
          _ghtag!(tag₁[r₁], i, nevents)
        elseif βinside
          # T-overlap: vertex of the first ring inside an edge of the second
          nevents += 1
          _ghins!(ins₂[r₂], j, (β, a₁, nevents))
          _ghtag!(tag₁[r₁], i, nevents)
        elseif αinside
          # T-overlap: vertex of the second ring inside an edge of the first
          nevents += 1
          _ghins!(ins₁[r₁], i, (α, a₂, nevents))
          _ghtag!(tag₂[r₂], j, nevents)
        elseif isapproxzero(α) && isapproxzero(β)
          # V-overlap: coincident vertices
          nevents += 1
          _ghtag!(tag₁[r₁], i, nevents)
          _ghtag!(tag₂[r₂], j, nevents)
        end
      end
    end
  end

  list₁, idict₁ = _ghlist(vs₁, tag₁, ins₁)
  list₂, idict₂ = _ghlist(vs₂, tag₂, ins₂)

  # link the two copies of each intersection vertex
  for (event, i) in idict₁
    if haskey(idict₂, event)
      j = idict₂[event]
      list₁.verts[i].neighbor = j
      list₂.verts[j].neighbor = i
    end
  end

  # vertices without a neighbor are not intersections
  for list in (list₁, list₂), v in list.verts
    iszero(v.neighbor) && (v.intersect = false)
  end

  list₁, list₂
end

_ghinunit(λ) = (λ > 0 || isapproxzero(λ)) && (λ < 1 && !isapproxone(λ))

_ghtag!(tag, i, event) = iszero(tag[i]) && (tag[i] = event)

_ghins!(ins, i, event) = push!(ins[i], event)

# build the list of vertices with the intersections inserted along the edges
function _ghlist(vs, tag, ins)
  P = eltype(first(vs))
  verts = GHVertex{P}[]
  comps = UnitRange{Int}[]
  idict = Dict{Int,Int}()
  for r in eachindex(vs)
    start = length(verts) + 1
    for i in eachindex(vs[r])
      thispoint = vs[r][i]
      thisevent = tag[r][i]
      intersect = !iszero(thisevent)
      push!(verts, GHVertex(point=thispoint, intersect=intersect))
      intersect && (idict[thisevent] = length(verts))
      for (_, point, event) in sort(ins[r][i], by=first)
        push!(verts, GHVertex(point=point, intersect=true))
        idict[event] = length(verts)
      end
    end
    push!(comps, start:length(verts))
  end

  compof = zeros(Int, length(verts))
  for (c, range) in enumerate(comps), i in range
    compof[i] = c
  end

  list = GHList(verts, comps, compof)

  list, idict
end

# ---------------
# LABELING PHASE
# ---------------

# mark intersection vertices as crossing or bouncing
function _ghmark!(list₁, list₂)
  types = map(eachindex(list₁.verts)) do i
    list₁.verts[i].intersect ? _ghlocaltype(list₁, list₂, i) : :none
  end

  for i in eachindex(list₁.verts)
    if types[i] == :crossing
      list₁.verts[i].crossing = true
    elseif types[i] ∈ (:lefton, :righton)
      # chain of vertices along a common segment
      x = types[i] == :lefton ? :left : :right
      j = _ghnext(list₁, i)
      while types[j] == :onon
        j = _ghnext(list₁, j)
      end
      y = types[j] == :onleft ? :left : (types[j] == :onright ? :right : :none)
      # delayed crossing if the polygon changes side along the chain
      x ≠ y && y ≠ :none && (list₁.verts[j].crossing = true)
    end
  end

  # the other polygon crosses at the same vertices
  for i in eachindex(list₁.verts)
    v = list₁.verts[i]
    v.intersect && (list₂.verts[v.neighbor].crossing = v.crossing)
  end
end

# local position of the first polygon with respect to the second one
function _ghlocaltype(list₁, list₂, i)
  p₋ = list₁.verts[_ghprev(list₁, i)]
  p₀ = list₁.verts[i]
  p₊ = list₁.verts[_ghnext(list₁, i)]

  k = p₀.neighbor
  q₋ = list₂.verts[_ghprev(list₂, k)]
  q₊ = list₂.verts[_ghnext(list₂, k)]

  # edges of the first polygon that overlap with the second one
  onnext = p₊.intersect && (p₊.neighbor == _ghprev(list₂, k) || p₊.neighbor == _ghnext(list₂, k))
  onprev = p₋.intersect && (p₋.neighbor == _ghprev(list₂, k) || p₋.neighbor == _ghnext(list₂, k))

  if onnext && onprev
    :onon
  elseif onnext
    q = p₊.neighbor == _ghnext(list₂, k) ? q₋ : q₊
    side = _ghsideof(q.point, p₋.point, p₀.point, p₊.point)
    side == :right ? :lefton : :righton
  elseif onprev
    q = p₋.neighbor == _ghprev(list₂, k) ? q₊ : q₋
    side = _ghsideof(q.point, p₋.point, p₀.point, p₊.point)
    side == :right ? :onleft : :onright
  else
    s₋ = _ghsideof(q₋.point, p₋.point, p₀.point, p₊.point)
    s₊ = _ghsideof(q₊.point, p₋.point, p₀.point, p₊.point)
    s₋ ≠ s₊ ? :crossing : :bouncing
  end
end

# wrap around the next component indices
function _ghnext(list, i)
  r = list.comps[list.compof[i]]
  i == last(r) ? first(r) : i + 1
end

# wrap around the previous component indices
function _ghprev(list, i)
  r = list.comps[list.compof[i]]
  i == first(r) ? last(r) : i - 1
end

# side of the point with respect to the chain of two edges
function _ghsideof(point, p₋, p₀, p₊)
  s₁ = signarea(point, p₋, p₀)
  s₂ = signarea(point, p₀, p₊)
  s₃ = signarea(p₋, p₀, p₊)
  if ispositive(s₃)
    ispositive(s₁) && ispositive(s₂) ? :left : :right
  elseif isnegative(s₃)
    ispositive(s₁) || ispositive(s₂) ? :left : :right
  else
    ispositive(s₁) ? :left : :right
  end
end

# label crossing vertices alternately as entry and exit points
function _ghflags!(list, rings)
  for range in list.comps
    start = _ghstart(list, range, rings)
    isnothing(start) && continue
    i, inside = start
    for _ in range
      v = list.verts[i]
      if v.crossing
        v.entry = !inside
        inside = !inside
      end
      i = _ghnext(list, i)
    end
  end
end

# first vertex of the component that is not an intersection, and its status
function _ghstart(list, range, rings)
  for i in range
    list.verts[i].intersect && continue
    point = list.verts[i].point
    any(r -> sideof(point, r) == ON, rings) && continue
    return i, _ghinside(point, rings)
  end
  # all vertices are intersections, use the midpoint of an edge
  for i in range
    j = _ghnext(list, i)
    m = _ghmidpoint(list.verts[i].point, list.verts[j].point)
    any(r -> sideof(m, r) == ON, rings) && continue
    return j, _ghinside(m, rings)
  end
  nothing
end

_ghmidpoint(p₁, p₂) = p₁ + (p₂ - p₁) / 2

# --------------
# TRACING PHASE
# --------------

# trace the components of the clipped polygon
function _ghtrace(list₁, list₂)
  P = typeof(first(list₁.verts).point)
  rings = Vector{P}[]
  for start in eachindex(list₁.verts)
    v = list₁.verts[start]
    # start at entry points so that components keep the orientation of the subject
    (v.intersect && v.crossing && v.entry && !v.visited) || continue

    points = P[]
    i, first₁ = start, true
    rounds = 0
    while (rounds += 1) ≤ length(list₁.verts) + length(list₂.verts)
      list, other = first₁ ? (list₁, list₂) : (list₂, list₁)
      v = list.verts[i]
      v.visited = true
      other.verts[v.neighbor].visited = true
      push!(points, v.point)

      # walk until the next crossing vertex with opposite flag
      forward = v.entry
      while true
        i = forward ? _ghnext(list, i) : _ghprev(list, i)
        w = list.verts[i]
        (w.crossing && w.entry ≠ forward) && break
        push!(points, w.point)
        w.visited = true
        w.intersect && (other.verts[w.neighbor].visited = true)
      end

      w = list.verts[i]
      w.visited = true
      other.verts[w.neighbor].visited = true

      # stop when the component is closed
      (first₁ && i == start) && break
      (!first₁ && w.neighbor == start) && break

      i = w.neighbor
      first₁ = !first₁
    end

    length(points) ≥ 3 && push!(rings, points)
  end
  rings
end

# components without crossing vertices that lie inside the other polygon
function _ghinner(list, rings; common)
  P = typeof(first(list.verts).point)
  inner = Vector{P}[]
  for range in list.comps
    any(i -> list.verts[i].crossing, range) && continue
    start = _ghstart(list, range, rings)
    if isnothing(start)
      # the component encloses the same region as a component of the other polygon
      common && push!(inner, [list.verts[i].point for i in range])
    else
      _, inside = start
      inside && push!(inner, [list.verts[i].point for i in range])
    end
  end
  inner
end

# ----------------
# MERGING RESULTS
# ----------------

_ghinside(point, rings) = isodd(count(r -> sideof(point, r) == IN, rings))

function _ghpolygons(rs)
  # drop collinear vertices introduced by degenerate intersections
  rings = [Ring(_ghsimplify(points)) for points in rs if length(_ghsimplify(points)) ≥ 3]
  isempty(rings) && return nothing

  outers = [r for r in rings if orientation(r) == CCW]
  inners = [r for r in rings if orientation(r) == CW]

  isempty(outers) && return nothing

  polys = map(outers) do outer
    holes = [inner for inner in inners if all(v -> sideof(v, outer) == IN, [first(eachvertex(inner))])]
    PolyArea([outer; holes])
  end

  length(polys) == 1 ? first(polys) : Multi(polys)
end

function _ghsimplify(points)
  n = length(points)
  n < 3 && return points
  keep = [i for i in 1:n if !isapproxzero(signarea(points[mod1(i - 1, n)], points[i], points[mod1(i + 1, n)]))]
  length(keep) < 3 ? points : points[keep]
end
