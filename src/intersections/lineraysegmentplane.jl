# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function lineraysegmentdirection(x::Line)
  return x(1) - x(0)
end

function lineraysegmentdirection(x::Ray)
  return direction(x)
end

function lineraysegmentdirection(x::Segment)
  return -reduce(-, coordinates.(vertices(x)))
end

function lineraysegmentorigin(x::Line)
  return x(0)
end

function lineraysegmentorigin(x::Ray)
  return origin(x)
end

function lineraysegmentorigin(x::Segment)
  return vertices(x)[1]
end

function lineraysegmentcrossing(x::Line)
  return CrossingLinePlane
end

function lineraysegmentcrossing(x::Ray)
  return CrossingRayPlane
end

function lineraysegmentcrossing(x::Segment)
  return CrossingSegmentPlane
end

function lineraysegmentoverlapping(x::Line)
  return OverlappingLinePlane
end

function lineraysegmentoverlapping(x::Ray)
  return OverlappingRayPlane
end

function lineraysegmentoverlapping(x::Segment)
  return OverlappingSegmentPlane
end

function lineraysegmenttouching(x::Ray)
  return TouchingRayPlane
end

function lineraysegmenttouching(x::Segment)
  return TouchingSegmentPlane
end

#=
(https://en.wikipedia.org/wiki/Line-plane_intersection)
=#
function Meshes.intersection(f, g::G, p::Plane{T}) where {T, G<:Union{Ray{3,T}, Line{3,T}, Segment{3,T}}}
  p₀ = coordinates(origin(p))
  n  = normal(p)

  # Get the origin and direction of geometry g
  g₀ = lineraysegmentorigin(g)
  gdir = lineraysegmentdirection(g)
  
  # calculate components
  ln = gdir ⋅ n
  
  # if ln is zero, then g is parallel to the plane
  if isapprox(ln, zero(T), atol=atol(T))
    # if the numerator is zero, then g is coincident
    if isapprox(coordinates(p₀ - g₀) ⋅ n, zero(T), atol=atol(T))
      return @IT lineraysegmentoverlapping(g) g f
    else
      return @IT NoIntersection nothing f
    end
  else
    # calculate the length parameter
    λ = -(n ⋅ coordinates(g₀ - p₀)) / ln

    if isa(g, Line)
      return @IT lineraysegmentcrossing(g) g(λ) f
    end

    # if λ is approximately 0, set as so to prevent any domain errors
    if isapprox(λ, zero(T), atol=atol(T))
      return @IT lineraysegmenttouching(g) g(zero(T)) f
    end

    if isa(g, Segment)
      # if λ is approximately 1, set as so to prevent any domain errors
      if isapprox(λ, one(T), atol=atol(T))
        return @IT lineraysegmenttouching(g) g(one(T)) f
      end
    end

    # if λ is out of bounds for the ray/segment, then there is no intersection
    if (λ < zero(T))
      return @IT NoIntersection nothing f
    else
      return @IT lineraysegmentcrossing(g) g(λ) f
    end
  end
end

