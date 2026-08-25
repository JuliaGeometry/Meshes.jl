# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, geom::Geometry, pset::PointSet)
  ps = filter(∈(geom), collect(pset))

  glue(pieces) = length(pieces) == 1 ? only(pieces) : Multi(pieces)

  if isempty(ps)
    return @IT NotIntersecting nothing f
  else
    return @IT Intersecting glue(ps) f
  end
end

intersection(f, dom::Domain, pset::PointSet) = intersection(f, Multi(collect(dom)), pset)

function intersection(f, dom₁::Domain, dom₂::Domain)
  # retrieve manifold and crs
  M = manifold(dom₁)
  C = crs(dom₁)

  # loop over all geometries
  gs = Geometry{M,C}[]
  for g₁ in dom₁, g₂ in dom₂
    g = g₁ ∩ g₂
    isnothing(g) || push!(gs, g)
  end

  # handle intersection at shared facets
  unique!(gs)

  glue(pieces) = length(pieces) == 1 ? only(pieces) : Multi(pieces)

  # return intersection
  if isempty(gs)
    return @IT NotIntersecting nothing f
  else
    return @IT Intersecting glue(gs) f
  end
end
