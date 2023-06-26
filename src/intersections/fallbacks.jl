# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

intersection(f, r::Ray, p::Polygon) = intersection(f, GeometrySet([r]), discretize(p))

function intersection(f, d₁::Domain{Dim,T}, d₂::Domain{Dim,T}) where {Dim,T}
  # loop over all geometries
  gs = Geometry{Dim,T}[]
  for g₁ in d₁, g₂ in d₂
    g = g₁ ∩ g₂
    isnothing(g) || push!(gs, g)
  end

  # handle intersection at shared facets
  unique!(gs)

  # return intersection
  if isempty(gs)
    return @IT NoIntersection nothing f
  else
    return @IT IntersectingGeometries Multi(gs) f
  end
end
