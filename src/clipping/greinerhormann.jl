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
mutable struct GHVertex{P<:Point}
  point::P
  inter::Bool
  neighbor::Int
  crossing::Bool
  entry::Bool
  visited::Bool
end

# doubly-linked list of vertices split into components (rings)
struct GHList{P<:Point}
  verts::Vector{GHVertex{P}}
  comps::Vector{UnitRange{Int}}
  compof::Vector{Int}
end

function _ghnext(list, i)
  r = list.comps[list.compof[i]]
  i == last(r) ? first(r) : i + 1
end

function _ghprev(list, i)
  r = list.comps[list.compof[i]]
  i == first(r) ? last(r) : i - 1
end

# -------------------
# INTERSECTION PHASE
# -------------------

# insert the intersections of the two sets of rings as vertices in both lists
function _ghintersect(rings₁, rings₂)
  vs₁ = [vertices(r) for r in rings₁]
  vs₂ = [vertices(r) for r in rings₂]

  # events attached to original vertices, and inserted along edges
  tags₁ = [zeros(Int, length(v)) for v in vs₁]
  tags₂ = [zeros(Int, length(v)) for v in vs₂]
  ins₁ = [[Tuple{Float64,eltype(v),Int}[] for _ in v] for v in vs₁]
  ins₂ = [[Tuple{Float64,eltype(v),Int}[] for _ in v] for v in vs₂]

  nevents = 0
  for r₁ in eachindex(vs₁), i in eachindex(vs₁[r₁])
    a₁ = vs₁[r₁][i]
    a₂ = vs₁[r₁][i + 1]
    for r₂ in eachindex(vs₂), j in eachindex(vs₂[r₂])
      b₁ = vs₂[r₂][j]
      b₂ = vs₂[r₂][j + 1]

      sa₁ = signarea(a₁, b₁, b₂)
      sa₂ = signarea(a₂, b₁, b₂)
      sb₁ = signarea(b₁, a₁, a₂)
      sb₂ = signarea(b₂, a₁, a₂)

      if !isapproxzero(sa₁ - sa₂) && !isapproxzero(sb₁ - sb₂) # not parallel
        α = sa₁ / (sa₁ - sa₂)
        β = sb₁ / (sb₁ - sb₂)
        (_ghinunit(α) && _ghinunit(β)) || continue
        αzero = isapproxzero(α)
        βzero = isapproxzero(β)
        if !αzero && !βzero
          # X-intersection: new vertex on both edges
          nevents += 1
          push!(ins₁[r₁][i], (α, a₁ + α * (a₂ - a₁), nevents))
          push!(ins₂[r₂][j], (β, b₁ + β * (b₂ - b₁), nevents))
        elseif !αzero
          # T-intersection: vertex of the second ring on an edge of the first
          nevents += 1
          push!(ins₁[r₁][i], (α, b₁, nevents))
          _ghtag!(tags₂[r₂], j, nevents)
        elseif !βzero
          # T-intersection: vertex of the first ring on an edge of the second
          nevents += 1
          push!(ins₂[r₂][j], (β, a₁, nevents))
          _ghtag!(tags₁[r₁], i, nevents)
        else
          # V-intersection: coincident vertices
          nevents += 1
          _ghtag!(tags₁[r₁], i, nevents)
          _ghtag!(tags₂[r₂], j, nevents)
        end
      elseif isapproxzero(sa₁) && isapproxzero(sa₂) && isapproxzero(sb₁) && isapproxzero(sb₂)
        # collinear edges, possibly overlapping
        u, v = a₂ - a₁, b₂ - b₁
        α = (b₁ - a₁) ⋅ u / (u ⋅ u)
        β = (a₁ - b₁) ⋅ v / (v ⋅ v)
        αin = _ghinunit(α) && !isapproxzero(α)
        βin = _ghinunit(β) && !isapproxzero(β)
        if αin && βin
          # X-overlap: both vertices lie inside the other edge
          nevents += 1
          push!(ins₁[r₁][i], (α, b₁, nevents))
          _ghtag!(tags₂[r₂], j, nevents)
          nevents += 1
          push!(ins₂[r₂][j], (β, a₁, nevents))
          _ghtag!(tags₁[r₁], i, nevents)
        elseif βin
          # T-overlap: vertex of the first ring inside an edge of the second
          nevents += 1
          push!(ins₂[r₂][j], (β, a₁, nevents))
          _ghtag!(tags₁[r₁], i, nevents)
        elseif αin
          # T-overlap: vertex of the second ring inside an edge of the first
          nevents += 1
          push!(ins₁[r₁][i], (α, b₁, nevents))
          _ghtag!(tags₂[r₂], j, nevents)
        elseif isapproxzero(α) && isapproxzero(β)
          # V-overlap: coincident vertices
          nevents += 1
          _ghtag!(tags₁[r₁], i, nevents)
          _ghtag!(tags₂[r₂], j, nevents)
        end
      end
    end
  end

  list₁, map₁ = _ghlist(vs₁, tags₁, ins₁)
  list₂, map₂ = _ghlist(vs₂, tags₂, ins₂)

  # link the two copies of each intersection vertex
  for (event, i) in map₁
    j = get(map₂, event, 0)
    j == 0 && continue
    list₁.verts[i].neighbor = j
    list₂.verts[j].neighbor = i
  end

  # vertices without a neighbor are not intersections
  for list in (list₁, list₂), v in list.verts
    iszero(v.neighbor) && (v.inter = false)
  end

  list₁, list₂
end

_ghtag!(tags, i, event) = iszero(tags[i]) && (tags[i] = event)

_ghinunit(λ) = (λ > 0 || isapproxzero(λ)) && (λ < 1 && !isapproxone(λ))

# build the list of vertices with the intersections inserted along the edges
function _ghlist(vs, tags, ins)
  P = eltype(first(vs))
  verts = GHVertex{P}[]
  comps = UnitRange{Int}[]
  map = Dict{Int,Int}()
  for r in eachindex(vs)
    start = length(verts) + 1
    for i in eachindex(vs[r])
      push!(verts, GHVertex(vs[r][i], !iszero(tags[r][i]), 0, false, false, false))
      iszero(tags[r][i]) || (map[tags[r][i]] = length(verts))
      for (_, point, event) in sort(ins[r][i], by=first)
        push!(verts, GHVertex(point, true, 0, false, false, false))
        map[event] = length(verts)
      end
    end
    push!(comps, start:length(verts))
  end

  compof = zeros(Int, length(verts))
  for (c, range) in enumerate(comps), i in range
    compof[i] = c
  end

  GHList(verts, comps, compof), map
end

# ---------------
# LABELING PHASE
# ---------------

# mark intersection vertices as crossing or bouncing
function _ghmark!(list₁, list₂)
  types = Vector{Symbol}(undef, length(list₁.verts))
  for i in eachindex(list₁.verts)
    types[i] = list₁.verts[i].inter ? _ghlocaltype(list₁, list₂, i) : :none
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
    v.inter && (list₂.verts[v.neighbor].crossing = v.crossing)
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
  onnext = p₊.inter && (p₊.neighbor == _ghprev(list₂, k) || p₊.neighbor == _ghnext(list₂, k))
  onprev = p₋.inter && (p₋.neighbor == _ghprev(list₂, k) || p₋.neighbor == _ghnext(list₂, k))

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
    list.verts[i].inter && continue
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
    (v.inter && v.crossing && v.entry && !v.visited) || continue

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
        w.inter && (other.verts[w.neighbor].visited = true)
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
