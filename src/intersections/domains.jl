# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, dom₁::Domain{Dim,T}, dom₂::Domain{Dim,T}) where {Dim,T}
  # loop over all geometries
  gs = Geometry{Dim,T}[]
  for g₁ in dom₁, g₂ in dom₂
    g = g₁ ∩ g₂
    isnothing(g) || push!(gs, g)
  end

  # handle intersection at shared facets
  unique!(gs)

  # return intersection
  if isempty(gs)
    return @IT NotIntersecting nothing f
  else
    return @IT Intersecting Multi(gs) f
  end
end
