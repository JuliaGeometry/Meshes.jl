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
# to the next  intersection vertex on the edge between point and left.point, or the following
# original ring vertex if no intersections on the edge. For Entering and Exiting types, left
# poitns to the following vertex on the clipping ring and right to the following vertex on one
# of the clipped rings.
mutable struct RingVertex{VT<:VertexType,M<:Manifold,C<:CRS}
  point::Point{M,C}
  left::RingVertex{VT,M,C} where {VT<:VertexType}
  right::RingVertex{VT,M,C} where {VT<:VertexType}

  function RingVertex{VT,M,C}(point::Point{M,C}) where {VT<:VertexType,M<:Manifold,C<:CRS}
    v = new(point)
    v.left = v
    v.right = v
    return v
  end
end

RingVertex{VT}(point::Point{M,C}) where {VT<:VertexType,M<:Manifold,C<:CRS} =
  RingVertex{VT,manifold(point),typeof(point.coords)}(point)

isnormal(::RingVertex{Normal}) = true
isnormal(::RingVertex) = false

function appendvertices!(v1::RingVertex{Normal}, v2::RingVertex{Normal})
  v2.left = v1.left
  v1.left = v2
  v2.right = v1.right
  v1.right = v2
end

function clip(poly::Polygon, ring::Ring, ::WeilerAthertonClipping)
  polyrings = rings(poly)

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

        vertextype = nothing
        I = intersection(clippingsegment, clippedsegment)

        # First try to handle Crossing, EdgeTouching and CornerTouching intersections as they might
        # add only a single intersection.

        if type(I) == Crossing
          point = get(I)
          cl = Line(clipping.left.point, clipping.left.left.point)
          vertextype = sideof(clipped.left.left.point, cl) == LEFT ? Entering : Exiting
        elseif type(I) == EdgeTouching
          point = get(I)
          if point ≈ clipped.left.point
            # When intersection is at the shared vertex of two edges of the clipped ring,
            # then split the interscting edge of the clipping ring at the intersection point.
            vertextype = decidedirection(
              Segment(clipped.point, clipped.left.point),
              clippedsegment,
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
              clippingsegment
            )
          end
        elseif type(I) == CornerTouching
          # When intersection is at the shared vertices of both the clipping and the clipped rings.
          point = get(I)
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
        end
        if !isnothing(vertextype)
          insertintersections!(clipping.left, clipped.left, point, vertextype, entering)
          if vertextype == Entering
            intersected[k] = true
          end
        end

        # Overlapping intersections might add up to two intersections, so handle separately.

        if type(I) == Overlapping
          for point in extrema(get(I))
            # For both ends of the intersecting segment, check if it coincides with the middle of
            # the obeserved vertices for clipped and clipping rings. If it does, attempt adding a
            # point.

            if point ≈ clipped.left.point
              clippingprev = Segment(clipping.left.point, point)
              if measure(clippingprev) ≈ 0.0 * u"m"
                clippingprev = Segment(clipping.point, point)
              end

              vertextype = decidedirection(
                Segment(clipped.point, point),
                Segment(point, clipped.left.left.point),
                clippingprev,
                Segment(point, clipping.left.left.point)
              )
              if !isnothing(vertextype)
                insertintersections!(clipping.left, clipped.left, point, vertextype, entering)
                if vertextype == Entering
                  intersected[k] = true
                end
              end
            end

            if point ≈ clipping.left.point
              clippedprev = Segment(clipped.left.point, point)
              if measure(clippedprev) ≈ 0.0 * u"m"
                clippedprev = Segment(clipped.point, point)
              end

              vertextype = decidedirection(
                clippedprev,
                Segment(point, clipped.left.left.point),
                Segment(clipping.point, point),
                Segment(point, clipping.left.left.point)
              )
              if !isnothing(vertextype)
                insertintersections!(clipping.left, clipped.left, point, vertextype, entering)
                if vertextype == Entering
                  intersected[k] = true
                end
              end
            end
          end
        end

        clipped = clipped.left
        if clipped.left == startclipped.left
          break
        end
      end
    end

    clipping = clipping.left
    if clipping.left == startclipping.left
      break
    end
  end

  collectedrings = collectclipped(entering)

  # When no interesction have been registered with any of the rings, the clipping ring either
  # encircles everything or is completely contained in the clipped polygon.
  if all(isequal(false), intersected)
    o = orientation(ring)
    contained = [v ∈ PolyArea(ring) for v in vertices(polyrings[1])]
    if (o == CCW && all(contained))
      return poly
    end
    if (o == CW && !all(contained))
      push!(collectedrings, polyrings[1])
      intersected[1] = true
      for polyring in polyrings[2:end]
        if !all([v ∈ PolyArea(ring) for v in vertices(polyring)])
          push!(collectedrings, polyring)
        end
      end
      if any([v ∈ PolyArea(polyrings[1]) for v in vertices(ring)])
        push!(collectedrings, ring)
      end
      return PolyArea(collectedrings...)
    end
    if all([v ∈ PolyArea(polyrings[1]) for v in vertices(ring)])
      push!(collectedrings, ring)
    end
  end

  # Handle all the non-intersected rings. Add them if they are contained in the clipping ring.
  polys = PolyArea[]
  for r in collectedrings
    newpolyrings = [r]
    valid = true
    for k in eachindex(intersected)
      if !intersected[k]
        if orientation(ring) == CCW
          # Check if majority of vertices are inside the clipping ring.
          ins = count(isequal(true), [v ∈ PolyArea(ring) for v in vertices(polyrings[k])])
          if ins > (length(vertices(polyrings[k])) - ins)
            if orientation(polyrings[k]) == CCW
              pushfirst!(newpolyrings, polyrings[k])
            else
              push!(newpolyrings, polyrings[k])
            end
            intersected[k] = true
          elseif vertices(ring)[1] ∈ PolyArea(polyrings[k]) && k > 1
            # If the ring is contained within one of the inner rings, invalidate it.
            valid = false
          end
        else
          # Check if majority of vertices are outside the clipping ring.
          ins = count(isequal(true), [v ∈ PolyArea(ring) for v in vertices(polyrings[k])])
          if ins < (length(vertices(polyrings[k])) - ins)
            if orientation(polyrings[k]) == CCW
              pushfirst!(newpolyrings, polyrings[k])
            else
              push!(newpolyrings, polyrings[k])
            end
            intersected[k] = true
          end
        end
      end
    end
    if valid
      push!(polys, PolyArea(newpolyrings))
    end
  end

  n = length(polys)
  out = n == 0 ? nothing : (n == 1 ? polys[1] : GeometrySet(polys))
  return out
end

function clip(poly::Polygon, other::Geometry, method::WeilerAthertonClipping)
  return _clip(poly, boundary(other), method)
end

function _clip(poly, multi::Multi, method::WeilerAthertonClipping)
  for r in parent(multi)
    poly = clip(poly, r, method)
    if isnothing(poly)
      return nothing
    end
  end
  return poly
end

_clip(poly, other, method) = clip(poly, other, method::WeilerAthertonClipping)

function clip(dom::Domain, other::Geometry, method::WeilerAthertonClipping)
  clipped = filter(!isnothing, [clip(geom, other, method) for geom in dom])
  return isempty(clipped) ? nothing : (length(clipped) == 1 ? clipped[1] : GeometrySet(clipped))
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
function insertintersections!(clipping, clipped, point, vtype, entering)
  intersection = RingVertex{vtype}(point)
  insertintersection!(clipping, intersection, :left)
  insertintersection!(clipped, intersection, :right)

  if vtype == Entering
    push!(entering, intersection)
  end
end

# Traversing and selecting following vertices for collecting of the rings after clipping depends
# on the vertex type. Helper nextvertex follows the correct path.
nextvertex(v::RingVertex{Normal}) = v.right
nextvertex(v::RingVertex{Entering}) = v.right
nextvertex(v::RingVertex{Exiting}) = v.left

# Takes a list of entering vertices and returns all rings that contain those vertices.
function collectclipped(entering::Vector{RingVertex{Entering}})
  rings::Vector{Ring} = []
  visited::Vector{RingVertex} = []
  for i in eachindex(entering)
    if entering[i] in visited
      continue
    end

    ring::Vector{RingVertex} = []
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
    newring::Vector{RingVertex} = [ring[1]]
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
      push!(rings, Ring([r.point for r in ring]...))
    end
  end
  return rings
end

# Used to figure out the type of the vertex to add when intersection is other than the crossing
# type. For input it takes four segments which all overlap on a single central vertex.
# Function measures the angles formed between the segments and checks wether the two clipped segments
# start and end in the same region or different regions, as separated by the clipping segments.
function decidedirection(clipped₁, clipped₂, clipping₁, clipping₂)
  # The input segments should all share one common vertex.
  # Form vectors from the common vertex to other vertices.
  a = minimum(clipped₁) - maximum(clipped₁)
  b = maximum(clipped₂) - minimum(clipped₂)
  c = minimum(clipping₁) - maximum(clipping₁)
  d = maximum(clipping₂) - minimum(clipping₂)

  if any(norm.([a, b, c, d]) .≈ 0.0 * u"m")
    # An zero length input segment found, no need to calculate.
    return nothing
  end

  β = mod(∠(a, b), 2pi * u"rad")
  γ = mod(∠(a, c), 2pi * u"rad")
  δ = mod(∠(a, d), 2pi * u"rad")

  if isapprox(γ, 0.0, atol=atol(0.0)) || isapprox(γ, 2pi, atol=atol(0.0))
    if δ < β < 2pi
      return Entering
    else
      return Exiting
    end
  elseif isapprox(δ, 0.0, atol=atol(0.0)) || isapprox(δ, 2pi, atol=atol(0.0))
    if γ < β < 2pi
      return Exiting
    else
      return Entering
    end
  end
  if γ < δ && (γ < β || isapprox(β, γ, atol=atol(0.0))) && (isapprox(β, δ, atol=atol(0.0)) || β < δ)
    return Exiting
  elseif δ < γ && (δ < β || isapprox(β, δ, atol=atol(0.0))) && (isapprox(β, γ, atol=atol(0.0)) || β < γ)
    return Entering
  end
  return nothing
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
