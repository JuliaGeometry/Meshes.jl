# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    WeilerAtherton()

The Weiler-Atherton algorithm for clipping polygons.

## References
 
* Weiler, Kevin and Atherton, Peter. 1977. [Hidden surface removal using polygon area sorting]
  (https://dl.acm.org/doi/10.1145/563858.563896)
"""
struct WeilerAtherton <: ClippingMethod end

function clip(poly::Polygon, other::Geometry, method::WeilerAtherton)
  if hasholes(poly)
    @error "not implemented"
  end

  c = [clip(ring, boundary(other), method) for ring in rings(poly)]
  r = filter(!isnothing, vcat(c...))

  if isempty(r)
    nothing
  elseif length(r) > 1
    Multi(PolyArea.(r))
  else
    PolyArea(r[1])
  end
end

function clip(ring::Ring{Dim,T}, other::Ring{Dim,T}, ::WeilerAtherton) where {Dim,T}
  # assumes both ring and other are CCW
  # TODO: Corner case without intersection
  # TODO: Corner case with colinear segments with intersections
  # TODO: Poly point on boundary
  # TODO: PolyArea with holes

  I, pᵣ, pₒ = _pairwiseintersections(ring, other)
  nᵣ = length(pᵣ)
  nₒ = length(pₒ)

  vᵣ = vertices(ring)
  vₒ = vertices(other)

  #TODO: not assume that the first point is not in poly boundary
  isentering = (sideof(vᵣ[pᵣ[1][2]], other) != IN)

  for i in eachindex(pᵣ)
    if pᵣ[i][1] == :INTERSECTION
      kind = isentering ? :INTERSECTION_ENTER : :INTERSECTION_EXIT
      pᵣ[i] = (kind, pᵣ[i][2])
      isentering = !isentering
    end
  end

  isvisitedᵣ = fill(false, length(pᵣ))
  u = Ring{Dim,T}[]

  for i in 1:nᵣ
    if pᵣ[i][1] == :INTERSECTION_ENTER && !isvisitedᵣ[i]
      vᵤ = Point{Dim,T}[]
      j = i
      while(true)
        isvisitedᵣ[j] = true
        push!(vᵤ, I[pᵣ[j][2]])
        j = mod1(j+1, nᵣ)
  
        while(pᵣ[j][1] == :POINT)
          push!(vᵤ, vᵣ[pᵣ[j][2]])
          j = mod1(j+1, nᵣ)
        end

        push!(vᵤ, I[pᵣ[j][2]])
        j = _pointindex(pᵣ[j][2], pₒ)
        j = mod1(j+1, nₒ)

        while(pₒ[j][1] == :POINT)
          push!(vᵤ, vₒ[pₒ[j][2]])
          j = mod1(j+1, nₒ)
        end

        if I[pₒ[j][2]] == I[pᵣ[i][2]]
          break
        end
        j = _pointindex(pₒ[j][2], pᵣ)
      end
      push!(u, Ring(unique(vᵤ)))
    end
  end
  u
end


# ----------------
# HELPER FUNCTIONS
# ----------------

function _pairwiseintersections(ring::Ring{Dim,T}, other::Ring{Dim,T}) where {Dim,T}
  vᵣ = vertices(ring)
  nᵣ = length(vᵣ)
  vₒ = vertices(other)
  nₒ = length(vₒ)

  I = []
  vᵣI = [Int[] for i=1:nᵣ]
  vₒI = [Int[] for i=1:nₒ]
  for i in 1:nᵣ
    for j in 1:nₒ
      sᵣ = Segment(vᵣ[i], vᵣ[i+1])
      sₒ = Segment(vₒ[j], vₒ[j+1])
      sI = sᵣ ∩ sₒ

      if !isnothing(sI)
        push!(I, sI)
        id = length(I)
        push!(vᵣI[i], id)
        push!(vₒI[j], id)
      end
    end
  end

  for i in 1:nᵣ
    sort!(vᵣI[i], by = j -> measure(Segment(vᵣ[i], I[j])))
  end

  for i in 1:nₒ
    sort!(vₒI[i], by = j -> measure(Segment(vₒ[i], I[j])))
  end

  pᵣ = Tuple{Symbol,Int}[]
  for i in 1:nᵣ
    push!(pᵣ, (:POINT, i))
    for j in vᵣI[i]
      push!(pᵣ, (:INTERSECTION, j))
    end
  end

  pₒ = Tuple{Symbol,Int}[]
  for i in 1:nₒ
    push!(pₒ, (:POINT, i))
    for j in vₒI[i]
      push!(pₒ, (:INTERSECTION, j))
    end
  end

  I, pᵣ, pₒ
end

function _pointindex(p, points)
  for i in eachindex(points)
    if points[i][1] != :POINT && points[i][2] == p
      return i
    end
  end
  nothing
end
