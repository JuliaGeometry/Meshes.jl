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

function clip(poly::Polygon{Dim,T}, other::Geometry{Dim,T}, method::WeilerAtherton) where {Dim, T}
  if hasholes(poly)
    outers = clip(rings(poly)[1], boundary(other), method)
    inners = vcat([clip(r, boundary(other), method) for r in rings(poly)[2:end]]...)
    if isempty(outers)
      return nothing
    end

    g = PolyArea{Dim,T}[]
    # build a PolyArea to each outer ring generated
    for o in outers
      ringsₒ = [o]
      # find inner rings inside the outer ring o 
      for i in inners
        if any([sideof(v, o) == IN for v in vertices(i)])
          push!(ringsₒ, i)
        end
      end
      push!(g, PolyArea(ringsₒ))
    end
    Multi(g)
  else
    outers = clip(rings(poly)[1], boundary(other), method)
    if isempty(outers)
      nothing
    else
      Multi(PolyArea.(outers))
    end
  end
end

function clip(ring::Ring{Dim,T}, other::Ring{Dim,T}, ::WeilerAtherton) where {Dim,T}
  vᵣ = orientation(ring) == CCW ? vertices(ring) : reverse(vertices(ring))
  vₒ = orientation(other) == CCW ? vertices(other) : reverse(vertices(other))

  I, pᵣ, pₒ = _pairwiseintersections(vᵣ, vₒ)
  nᵣ = length(pᵣ)
  nₒ = length(pₒ)

  if isempty(I)
    # ring is totally inside or totally outside other
    if sideof(vᵣ[1], other) == IN 
      return [ring]
    else
      return []
    end
  end

  isentering = (sideof(vᵣ[pᵣ[1][2]], other) != IN)

  Iposᵣ = zeros(Int, length(I))
  for i in 1:nᵣ
    if pᵣ[i][1] == :INTERSECTION
      Iposᵣ[pᵣ[i][2]] = i
    end
  end

  Iposₒ = zeros(Int, length(I))
  for i in 1:nₒ
    if pₒ[i][1] == :INTERSECTION
      Iposₒ[pₒ[i][2]] = i
    end
  end

  # specify intersection type
  for i in eachindex(pᵣ)
    if pᵣ[i][1] == :INTERSECTION
      kind = isentering ? :INTERSECTION_ENTER : :INTERSECTION_EXIT
      pᵣ[i] = (kind, pᵣ[i][2])
      isentering = !isentering
    end
  end

  # mark intersection of type INTERSECTION_ENTER
  isvisitedᵣ = fill(false, length(pᵣ))
  u = Ring{Dim,T}[]

  # iterate over all intersections of type INTERSECTION_ENTER
  for i in 1:nᵣ
    if pᵣ[i][1] == :INTERSECTION_ENTER && !isvisitedᵣ[i]
      vᵤ = Point{Dim,T}[]
      j = i
      while(true)
        isvisitedᵣ[j] = true

        # add current ring intersection
        push!(vᵤ, I[pᵣ[j][2]])
        j = mod1(j+1, nᵣ)

        # collect all points of ring until find a intersection
        while(pᵣ[j][1] == :POINT)
          push!(vᵤ, vᵣ[pᵣ[j][2]])
          j = mod1(j+1, nᵣ)
        end

        # add intersection to ring and swap j to other
        push!(vᵤ, I[pᵣ[j][2]])
        j = Iposₒ[pᵣ[j][2]]
        j = mod1(j+1, nₒ)

        # collect all points of other until find a intersection
        while(pₒ[j][1] == :POINT)
          push!(vᵤ, vₒ[pₒ[j][2]])
          j = mod1(j+1, nₒ)
        end

        # stop if reach initial intersection 
        if I[pₒ[j][2]] == I[pᵣ[i][2]]
          break
        end
        
        # swap j back to ring
        j = Iposᵣ[pₒ[j][2]]
      end
      rᵤ = orientation(ring) == CCW ? Ring(unique(vᵤ)) : Ring(reverse(unique(vᵤ)))
      push!(u, rᵤ)
    end
  end
  u
end


# ----------------
# HELPER FUNCTIONS
# ----------------

function _pairwiseintersections(vᵣ, vₒ)
  nᵣ = length(vᵣ)
  nₒ = length(vₒ)

  I = []
  vᵣI = [Int[] for i=1:nᵣ]
  vₒI = [Int[] for i=1:nₒ]
  for i in 1:nᵣ
    p₁, p₂ = vᵣ[i], vᵣ[i+1]
    sᵣ = Segment(p₁, p₂)
    for j in 1:nₒ
      sₒ = Segment(vₒ[j], vₒ[j+1])
      lₒ = Line(vₒ[j], vₒ[j+1])
      sI = sᵣ ∩ sₒ

      if !isnothing(sI) && ((sideof(p₁, lₒ) == RIGHT) ⊻ (sideof(p₂, lₒ) == RIGHT))
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
