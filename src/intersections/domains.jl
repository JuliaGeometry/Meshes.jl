# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, geom::Geometry, pset::PointSet)
  ps = filter(∈(geom), collect(pset))
  if isempty(ps)
    return @IT NotIntersecting nothing f
  else
    return @IT Intersecting PointSet(ps) f
  end
end

intersection(f, dom::Domain, pset::PointSet) = intersection(f, Multi(collect(dom)), pset)

function intersection(f, dom₁::Domain{Dim}, dom₂::Domain{Dim}) where {Dim}
  # loop over all geometries
  gs = Geometry{Dim}[]
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
    return @IT Intersecting GeometrySet(gs) f
  end
end
