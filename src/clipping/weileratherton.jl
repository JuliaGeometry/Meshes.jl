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

  pᵣ = _find_and_insert_intersections(vᵣ, vₒ)
  pₒ = _find_and_insert_intersections(vₒ, vᵣ)
  nᵣ, nₒ = length(pᵣ), length(pₒ)

  if nᵣ == nvertices(ring)
    # ring is totally inside or totally outside other
    if sideof(vᵣ[1], other) == IN 
      return [ring]
    else
      return []
    end
  end

  Iposᵣ = Dict{Intersection, Int}()
  for i in 1:nᵣ
    if pᵣ[i] isa Intersection
      Iposᵣ[pᵣ[i]] = i
    end
  end

  Iposₒ = Dict{Intersection, Int}()
  for i in 1:nₒ
    if pₒ[i] isa Intersection
      Iposₒ[pₒ[i]] = i
    end
  end

  # mark intersections of type INTERSECTION_ENTER
  isvisited = fill(false, length(pᵣ))
  isentering = (sideof(vᵣ[1], other) == OUT)
  u = Ring{Dim,T}[]

  for i in 1:nᵣ
    if pᵣ[i] isa Intersection
      if isentering && !isvisited[i]
        vᵤ = _walkloop(pᵣ, pₒ, i, isvisited)
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
  pᵣ = []
  for i in eachindex(vᵣ)
    push!(pᵣ, vᵣ[i])
    I = []
    for j in eachindex(vₒ)
      sᵣ = Segment(vᵣ[i], vᵣ[i+1])
      sₒ = Segment(vₒ[j], vₒ[j+1])
      lₒ = Line(vₒ[j], vₒ[j+1])
      sI = intersection(sᵣ, sₒ)

      if get(sI) isa Point && ((sideof(vᵣ[i], lₒ) == RIGHT) ⊻ (sideof(vᵣ[i+1], lₒ) == RIGHT))
        push!(I, sI)
      end
    end
    # sort intersections by distance
    sort!(I, by = p -> measure(Segment(vᵣ[i], get(p))))
    for p in I
      push!(pᵣ, p)
    end
  end
  pᵣ
end

function _walkloop(pᵣ, pₒ, i, isvisited)
  nᵣ = length(pᵣ)
  nₒ = length(pₒ)
  vᵤ = []
  j = i
  while(true)
    isvisited[j] = true
    # add current ring intersection
    I = get(pᵣ[j])
    push!(vᵤ, I)
    j = mod1(j+1, nᵣ)
    # collect points of ring until find a intersection
    while(pᵣ[j] isa Point)
      push!(vᵤ, pᵣ[j])
      j = mod1(j+1, nᵣ)
    end

    # add intersection to ring and swap j to other
    I = get(pᵣ[j])
    push!(vᵤ, I)
    #j = Iposₒ[I]
    j = _findIindex(I, pₒ)
    j = mod1(j+1, nₒ)

    # collect points of other until find a intersection
    while(pₒ[j] isa Point)
      push!(vᵤ, pₒ[j])
      j = mod1(j+1, nₒ)
    end
    # swap j back to ring
    I = get(pₒ[j])
    #j = Iposᵣ[I]
    j = _findIindex(I, pᵣ)

    # stop if reach initial intersection 
    if i == j
      break
    end
  end
  vᵤ
end

function _findIindex(I, v)
  for i in eachindex(v)
    if v[i] isa Intersection && get(v[i]) ≈ I
      return i
    end
  end
  nothing
end
