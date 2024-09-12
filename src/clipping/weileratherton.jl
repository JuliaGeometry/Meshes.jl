# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    WeilerAthertonClipping()

The Wieler-Atherton algorithm for clipping polygons.

## References

* Weiler, K., & Atherton, P. 1977. [Hidden surface removal using polygon area sorting]
  (https://dl.acm.org/doi/pdf/10.1145/563858.563896)
"""
struct WeilerAthertonClipping <: ClippingMethod end

# VertexType has three different types of vertices used for deciding how to navigate the data
# structure when collecting the polygon rings after clipping.
abstract type VertexType end
abstract type Normal <: VertexType end
abstract type Entering <: VertexType end
abstract type Exiting <: VertexType end

# Data structure for clipping the polygons. Fields left and right are used depending on the
# VertexType. If Normal, left points to the following vertex of the original ring, and right
# to the next intersection vertex on the edge between point and left.point, or the following
# original ring vertex if no intersections on the edge. For Entering and Exiting types, left
# poitns to the following vertex on the clipping ring and right to the following vertex on one
# of the clipped rings.
# TODO Either properly document the usage of left and right, or separate RingVertex into Entering, Exiting and Normal vertices.
mutable struct RingVertex{VT<:VertexType,M<:Manifold,C<:CRS}
  point::Point{M,C}
  left::RingVertex{<:VertexType,M,C}
  right::RingVertex{<:VertexType,M,C}

  function RingVertex{VT,M,C}(point) where {VT<:VertexType,M<:Manifold,C<:CRS}
    v = new(point)
    v.left = v
    v.right = v
    v
  end
end

RingVertex{VT}(point::Point{M,C}) where {VT<:VertexType,M<:Manifold,C<:CRS} = RingVertex{VT,M,C}(point)

isnormal(::RingVertex{Normal}) = true
isnormal(::RingVertex) = false

function appendvertices!(v₁::RingVertex{Normal}, v₂::RingVertex{Normal})
  v₂.left = v₁.left
  v₁.left = v₂
  v₂.right = v₁.right
  v₁.right = v₂
end

# Traversing and selecting following vertices for collecting of the rings after clipping depends
# on the vertex type. Helper nextvertex follows the correct path.
nextvertex(v::RingVertex{Normal}) = v.right
nextvertex(v::RingVertex{Entering}) = v.right
nextvertex(v::RingVertex{Exiting}) = v.left

function clip(poly::Polygon, ring::Ring, ::WeilerAthertonClipping)
  polyrings = rings(poly)

  # If the polygon is contained in the ring, return the polygon right away.
  allcontained = all(v ∈ PolyArea(ring) for v in vertices(polyrings[1]))
  if (orientation(ring) == CCW && allcontained)
    return poly
  end

  # Convert the subject polygon rings and the clipping ring to the RingVertex data structure.
  clippedrings = [gettraversalring(r) for r in polyrings]
  startclipping = gettraversalring(ring)

  # For keeping track of intersected rings, as the non-intersected ones need additional
  # processing at the end.
  intersected = zeros(Bool, length(clippedrings))

  # For collecting all entering vertices, as they are used as starting points for collection
  # of the output rings.
  entering = RingVertex{Entering}[]

  clipping = startclipping
  while true
    # Three consecutive clipping vertices are used to properly identify intersection vertex type
    # for corner cases, so the clipping segment is constructed with the two following vertices
    # after the current one.
    clippingsegment = Segment(clipping.left.point, clipping.left.left.point)

    for (k, startclipped) in enumerate(clippedrings)
      clipped = startclipped
      while true
        # Like for the clipping, the clipped also uses three consecutive vertices.
        clippedsegment = Segment(clipped.left.point, clipped.left.left.point)

        I = intersection(clippingsegment, clippedsegment)
        vertex = vertexfromintersection(I, clipping, clipped)
        success = insertintersections!(vertex, clipping, clipped, entering)
        intersected[k] = intersected[k] || success

        clipped = clipped.left
        clipped.left == startclipped.left && break
      end
    end

    clipping = clipping.left
    clipping.left == startclipping.left && break
  end

  # Handle the case when no interections have been found.
  if !any(intersected)
    if orientation(ring) == CW
      # For an inner clipping ring take all outside rings.
      return PolyArea(collectoutsiderings(ring, polyrings)...)
    else
      # For an outer clipping ring add it to act as the outer ring.
      collectedrings = all(v ∈ PolyArea(polyrings[1]) for v in vertices(ring)) ? [ring] : []
    end
  else
    # Collect rings formed from the intersected rings.
    collectedrings = collectclipped(entering)
  end

  # Complete the collected rings by adding non-intersected rings
  # if they are contained within the collected ones.
  completedrings = [perhapsaddnonintersected!(r, polyrings, intersected) for r in collectedrings]

  # Convert completed ring lists into PolyAreas.
  polys = PolyArea.(filter(!isnothing, completedrings))

  n = length(polys)
  n == 0 ? nothing : (n == 1 ? polys[1] : Multi(polys))
end

clip(poly::Polygon, other::Geometry, method::WeilerAthertonClipping) = clip(poly, boundary(other), method)

function collectoutsiderings(ring, polyrings)
  newpolyrings = Ring[]
  for k in eachindex(polyrings)
    if !all(v ∈ PolyArea(ring) for v in vertices(polyrings[k]))
      push!(newpolyrings, polyrings[k])
    end
  end
  if any(v ∈ PolyArea(polyrings[1]) for v in vertices(ring))
    push!(newpolyrings, ring)
  end
  newpolyrings
end

function perhapsaddnonintersected!(ring, polyrings, intersected)
  # Handle all the non-intersected rings. Add them if they are contained in the clipping ring.
  newpolyrings = [ring]
  for k in eachindex(intersected)
    if !intersected[k]
      ccw = orientation(ring) == CCW

      # Discard if the processed ring is an outer ring but inside another inner ring.
      if ccw && vertices(ring)[1] ∈ PolyArea(polyrings[k]) && k > 1
        return nothing
      end

      vs = vertices(polyrings[k])
      ins = count(v ∈ PolyArea(ring) for v in vs)
      outs = length(vs) - ins
      if ccw == (ins > outs)
        if !ccw
          # An outer ring should be first in the list.
          pushfirst!(newpolyrings, polyrings[k])
        else
          # Inner rings should follow an outer in the list.
          push!(newpolyrings, polyrings[k])
        end
        intersected[k] = true
      end
    end
  end
  newpolyrings
end

# Inserts the intersection in the ring.
function insertintersection!(head::RingVertex, intersection::RingVertex, side::Symbol)
  tail = head.left
  vertex = head

  new = measure(Segment(head.point, intersection.point))
  while true
    os = isnormal(vertex) ? :right : side
    current = measure(Segment(head.point, getfield(vertex, os).point))
    if (new < current) || (getfield(vertex, os) == tail)
      next = getfield(vertex, os)
      if !(new ≈ measure(Segment(head.point, next.point)))
        setfield!(intersection, side, next)
        setfield!(vertex, os, intersection)
      end
      break
    end
    vertex = getfield(vertex, os)
  end
end

# Inserts the intersection into both the clipping and the clipped rings.
function insertintersections!(vertex::Tuple, clipping, clipped, entering)
  (vtype, point) = vertex
  if !isnothing(vtype)
    intersection = RingVertex{vtype}(point)
    insertintersection!(clipping.left, intersection, :left)
    insertintersection!(clipped.left, intersection, :right)

    if vtype == Entering
      push!(entering, intersection)
    end
    return true
  end
  false
end

insertintersections!(vertices::Array, clipping, clipped, entering) =
  any(insertintersections!.(vertices, Ref(clipping), Ref(clipped), Ref(entering)))

# Takes a list of entering vertices and returns all rings that contain those vertices.
function collectclipped(entering::Vector{RingVertex{Entering}})
  rings = Ring[]
  visited = RingVertex[]
  for i in eachindex(entering)
    if entering[i] in visited
      continue
    end

    ring = RingVertex[]
    vertex = entering[i]
    while !(vertex in ring)
      if vertex in visited
        break
      end

      push!(ring, vertex)
      push!(visited, vertex)
      vertex = nextvertex(vertex)
    end

    # Remove duplicates.
    newring = RingVertex[ring[1]]
    for i in 2:length(ring)
      if !(ring[i].point ≈ newring[end].point)
        push!(newring, ring[i])
      end
    end
    ring = newring

    # Polygon might start several vertices after the first collected.
    # This generally happens when there are overlapping edges that lead
    # to several entering vertices without exiting in between. Then, the
    # actual polygon is found by discarding the extra vertices before the
    # proper loop.
    k = findfirst(x -> ring[end].point == x.point, ring[1:(end - 1)])
    if !isnothing(k)
      ring = ring[(k + 1):end]
    end

    if length(ring) > 2
      push!(rings, Ring([r.point for r in ring]))
    end
  end
  rings
end

function vertexfromintersection(I, clipping, clipped)
  type(I) == Crossing && return vertexfromcrossing(get(I), clipping, clipped)
  type(I) == CornerTouching && return vertexfromcornertouching(get(I), clipping, clipped)
  type(I) == EdgeTouching && return vertexfromedgetouching(get(I), clipping, clipped)
  type(I) == Overlapping && return vertexfromoverlapping(get(I), clipping, clipped)
  (nothing, nothing)
end

function vertexfromcrossing(point, clipping, clipped)
  cl = Line(clipping.left.point, clipping.left.left.point)
  vertextype = sideof(clipped.left.left.point, cl) == LEFT ? Entering : Exiting
  (vertextype, point)
end

function vertexfromedgetouching(point, clipping, clipped)
  vertextype = nothing
  if point ≈ clipped.left.point
    # When intersection is at the shared vertex of two edges of the clipped ring,
    # then split the interscting edge of the clipping ring at the intersection point.
    vertextype = decidedirection(
      Segment(clipped.point, clipped.left.point),
      Segment(clipped.left.point, clipped.left.left.point),
      Segment(clipping.left.point, point),
      Segment(point, clipping.left.left.point)
    )
  elseif point ≈ clipping.left.point
    # When intersection is at the shared vertex of two edges of the clipping ring,
    # then split the interscting edge of the clipped ring at the intersection point.
    vertextype = decidedirection(
      Segment(clipped.left.point, point),
      Segment(point, clipped.left.left.point),
      Segment(clipping.point, clipping.left.point),
      Segment(clipping.left.point, clipping.left.left.point)
    )
  end
  (vertextype, point)
end

function vertexfromcornertouching(point, clipping, clipped)
  vertextype = nothing
  # When intersection is at the shared vertices of both the clipping and the clipped rings.
  if (point ≈ clipped.left.point) && (point ≈ clipping.left.point)
    # Only applies if the intersection coincides with the middles of the currently observed
    # vertices for both clipping and clipped rings.
    vertextype = decidedirection(
      Segment(clipped.point, point),
      Segment(point, clipped.left.left.point),
      Segment(clipping.point, point),
      Segment(point, clipping.left.left.point)
    )
  end
  (vertextype, point)
end

function vertexfromoverlapping(segment, clipping, clipped)
  # For both ends of the intersecting segment, check if it coincides with the middle of
  # the obeserved vertices for clipped and clipping rings. If it does, attempt adding a
  # point.

  ret = Tuple[]
  for point in extrema(segment)
    if point ≈ clipped.left.point
      clippingprev = Segment(clipping.left.point, point)
      if measure(clippingprev) ≈ 0.0u"m"
        clippingprev = Segment(clipping.point, point)
      end
      vertextype = decidedirection(
        Segment(clipped.point, point),
        Segment(point, clipped.left.left.point),
        clippingprev,
        Segment(point, clipping.left.left.point)
      )
      push!(ret, (vertextype, point))
    end

    if point ≈ clipping.left.point
      clippedprev = Segment(clipped.left.point, point)
      if measure(clippedprev) ≈ 0.0u"m"
        clippedprev = Segment(clipped.point, point)
      end
      vertextype = decidedirection(
        clippedprev,
        Segment(point, clipped.left.left.point),
        Segment(clipping.point, point),
        Segment(point, clipping.left.left.point)
      )
      push!(ret, (vertextype, point))
    end
  end
  ret
end

# Used to figure out the type of the vertex to add when intersection is other than the crossing
# type. For input it takes four segments which all overlap on a single central vertex.
# Function measures the angles formed between the segments and checks wether the two clipped segments
# start and end in the same region or different regions, as separated by the clipping segments.
function decidedirection(clipped₁, clipped₂, clipping₁, clipping₂)
  T = numtype(lentype(clipped₁))

  # The input segments should all share one common vertex.
  # Form vectors from the common vertex to other vertices.
  a = minimum(clipped₁) - maximum(clipped₁)
  b = maximum(clipped₂) - minimum(clipped₂)
  c = minimum(clipping₁) - maximum(clipping₁)
  d = maximum(clipping₂) - minimum(clipping₂)

  if any(norm(v) ≈ 0.0u"m" for v in (a, b, c, d))
    # An zero length input segment found, no need to calculate.
    return nothing
  end

  tol = atol(0.0) * u"rad"
  twoπ = 2 * T(π) * u"rad"

  β = mod(∠(a, b), twoπ)
  γ = mod(∠(a, c), twoπ)
  δ = mod(∠(a, d), twoπ)

  if isapprox(γ, zero(γ), atol=tol) || isapprox(γ, twoπ, atol=tol)
    if δ < β < twoπ
      return Entering
    else
      return Exiting
    end
  elseif isapprox(δ, zero(δ), atol=tol) || isapprox(δ, twoπ, atol=tol)
    if γ < β < twoπ
      return Exiting
    else
      return Entering
    end
  end

  if γ < δ && (γ < β || isapprox(β, γ, atol=tol)) && (isapprox(β, δ, atol=tol) || β < δ)
    return Exiting
  elseif δ < γ && (δ < β || isapprox(β, δ, atol=tol)) && (isapprox(β, γ, atol=tol) || β < γ)
    return Entering
  end

  nothing
end

# Converts a regular Meshes.Ring into a ring formed with RingVertex data type.
function gettraversalring(ring::Ring)
  vs = vertices(ring)
  start = RingVertex{Normal}(vs[1])

  for v in vs[2:end]
    new = RingVertex{Normal}(v)
    appendvertices!(start, new)
    start = new
  end

  return start.left
end
