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

  pᵣ, pₒ = _find_and_insert_intersections(vᵣ, vₒ)
  nᵣ= length(pᵣ)

  # corner case without intersections
  if nᵣ == nvertices(ring)
    if sideof(vᵣ[1], other) == IN
      # ring is totally inside other
      return [ring]
    else
      # ring is totally outside other
      return []
    end
  end

  Iposᵣ = _build_dict_intersection_to_position(pᵣ)
  Iposₒ = _build_dict_intersection_to_position(pₒ)

  u = Ring{Dim,T}[]

  isvisited = fill(false, length(pᵣ))
  isentering = (sideof(vᵣ[1], other) == OUT)

  for i in 1:nᵣ
    if pᵣ[i] isa Intersection
      if isentering && !isvisited[i]
        # build new ring starting from the current entering intersection 
        vᵤ = _walkloop(pᵣ, pₒ, i, isvisited, Iposᵣ, Iposₒ)
        rᵤ = orientation(ring) == CCW ? Ring(unique(vᵤ)...) : Ring(reverse(unique(vᵤ)...))

        push!(u, rᵤ)
      end
      isentering = !isentering
    end
  end
  u
end

# ----------------
# HELPER FUNCTIONS
# ----------------

function _find_and_insert_intersections(vᵣ, vₒ)
  nᵣ = length(vᵣ)
  nₒ = length(vₒ)

  vᵣI = [Intersection[] for i=1:nᵣ]
  vₒI = [Intersection[] for i=1:nₒ]

  for i in 1:nᵣ
    for j in 1:nₒ
      sᵣ = Segment(vᵣ[i], vᵣ[i+1])
      sₒ = Segment(vₒ[j], vₒ[j+1])
      lₒ = Line(vₒ[j], vₒ[j+1])
      I = intersection(sᵣ, sₒ)

      if get(I) isa Point && ((sideof(vᵣ[i], lₒ) == RIGHT) ⊻ (sideof(vᵣ[i+1], lₒ) == RIGHT))
        push!(vᵣI[i], I)
        push!(vₒI[j], I)
      end
    end
  end

  pᵣ = []
  for i in 1:nᵣ
    push!(pᵣ, vᵣ[i])
    sort!(vᵣI[i], by = p -> measure(Segment(vᵣ[i], get(p))))
    for I in vᵣI[i]
      push!(pᵣ, I)
    end
  end

  pₒ = []
  for i in 1:nₒ
    push!(pₒ, vₒ[i])
    sort!(vₒI[i], by = p -> measure(Segment(vₒ[i], get(p))))
    for I in vₒI[i]
      push!(pₒ, I)
    end
  end

  pᵣ, pₒ
end

function _build_dict_intersection_to_position(p)
  Ipos = Dict{Intersection, Int}()
  for i in eachindex(p)
    if p[i] isa Intersection
      Ipos[p[i]] = i
    end
  end
  Ipos
end

function _walkloop(pᵣ, pₒ, i, isvisited, Iposᵣ, Iposₒ)
  nᵣ = length(pᵣ)
  nₒ = length(pₒ)
  vᵤ = []

  j = i
  while(true)
    isvisited[j] = true

    # add current intersection
    I = get(pᵣ[j])
    push!(vᵤ, I)
    j = mod1(j+1, nᵣ)

    # walk through points of pᵣ until find a intersection
    while(pᵣ[j] isa Point)
      push!(vᵤ, pᵣ[j])
      j = mod1(j+1, nᵣ)
    end

    # swap to j from pᵣ to pₒ
    j = Iposₒ[pᵣ[j]]

    # add current intersection
    I = get(pₒ[j])
    push!(vᵤ, I)
    j = mod1(j+1, nₒ)

    # walk through points of pₒ until find a intersection
    while(pₒ[j] isa Point)
      push!(vᵤ, pₒ[j])
      j = mod1(j+1, nₒ)
    end

    # swap to j from pₒ to pᵣ
    j = Iposᵣ[pₒ[j]]

    # stop if reach initial intersection 
    if i == j
      break
    end
  end

  vᵤ
end

