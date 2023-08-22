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

function clip(ring::Ring{Dim,T}, other::Ring{Dim,T}, ::WeilerAtherton) where {Dim,T}
  # assumes both ring and other are CCW
  # TODO: Corner case without intersection
  # TODO: Corner case with colinear segments with intersections
  # TODO: Poly point on boundary
  # TODO: PolyArea with holes
  # TODO: Optimize findpointindex

  vᵣ, kindᵣ = verticesswithintersections(ring, other)
  vₒ, kindₒ = verticesswithintersections(other, ring)
  nᵣ = length(vᵣ)
  nₒ = length(vₒ)

  isentering = nothing
  # set specific kinds of intersections to kindᵣ
  for i in eachindex(vᵣ)
    if kindᵣ[i] == :POINT
      isentering = (sideof(vᵣ[i], other) == OUT)
    else
      kindᵣ[i] = isentering ? :INTERSECTION_ENTER : :INTERSECTION_EXIT
      isentering = !isentering
    end
  end

  isvisitedᵣ = fill(false, length(vᵣ))
  u = []

  for i in 1:nᵣ
    if kindᵣ[i] == :INTERSECTION_ENTER && !isvisitedᵣ[i]
      vᵤ = Point{Dim,T}[]
      j = i
      while(true)
        isvisitedᵣ[j] = true
        while(kindᵣ[j] != :INTERSECTION_EXIT)
          push!(vᵤ, vᵣ[j])
          j = mod1(j+1, nᵣ)
        end

        push!(vᵤ, vᵣ[j])
        j = findpointindex(vᵣ[j], vₒ)
        j = mod1(j+1, nₒ)
        
        while(kindₒ[j] != :INTERSECTION)
          push!(vᵤ, vₒ[j])
          j = mod1(j+1, nₒ)
        end

        if vₒ[j] ≈ vᵣ[i]
          break
        end

        j = findpointindex(vₒ[j], vᵣ)
      end
      push!(u, Ring(vᵤ))
    end
  end
  u
end


# ----------------
# HELPER FUNCTIONS
# ----------------

function verticesswithintersections(ring::Ring{Dim,T}, other::Ring{Dim,T}) where {Dim,T}
  vᵣ = vertices(ring)
  nᵣ = length(vᵣ)
  vₒ = vertices(other)
  nₒ = length(vₒ)

  # vector of points and intersection
  vᵢ = Point{Dim,T}[]
  kind = Symbol[]

  # build path to poly (points and intersections)  
  for i in 1:nᵣ
    p₁, p₂ = vᵣ[i], vᵣ[i+1]
    sᵣ = Segment(p₁, p₂)

    # intersections from sᵣ
    I = []

    for j in 1:nₒ
      sₒ = Segment(vₒ[j], vₒ[j+1])

      if !isnothing(sᵣ ∩ sₒ)
        push!(I, sᵣ ∩ sₒ)
      end
    end

    # sort intersections by distance from p₁
    sort!(I, by = p -> measure(Segment(p₁, p)))

    # add points and kind
    push!(vᵢ, p₁)
    push!(kind, :POINT)
    for p in I 
      push!(vᵢ, p)
      push!(kind, :INTERSECTION)
    end
  end

  vᵢ, kind
end

function findpointindex(p, points)
  for i in eachindex(points)
    if points[i] ≈ p
      return i
    end
  end
  nothing
end
