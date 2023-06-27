# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, geom::Geometry{Dim,T}, pset::PointSet{Dim,T}) where {Dim,T}
  ps = filter(∈(geom), collect(pset))
  if isempty(ps)
    return @IT NotIntersecting nothing f
  else
    return @IT Intersecting PointSet(ps) f
  end
end

function intersection(f, pset::PointSet{Dim,T}, dom::Domain{Dim,T}) where {Dim,T}
  # loop over points once
  ps = Point{Dim,T}[]
  for p in pset
    for g in dom
      if p ∈ g
        push!(ps, p)
        break # interrupt inner loop
      end
    end
  end

  if isempty(ps)
    return @IT NotIntersecting nothing f
  else
    return @IT Intersecting PointSet(ps) f
  end
end

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
    return @IT Intersecting GeometrySet(gs) f
  end
end
