# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

intersection(f, chain₁::Chain, chain₂::Chain) =
  intersection(f, GeometrySet(segments(chain₁)), GeometrySet(segments(chain₂)))

# intersection between a chain and a polygon
# 1. intersect at inner points (Intersecting)
# 2. intersect at boundary points (Touching)
# 3. overlaps at the polygon edges (EdgeTouching)
# 4. no intersection (NotIntersecting)
function intersection(f, chain::Chain, poly::Polygon)
  # retrieve manifold and crs
  M = manifold(chain)
  C = crs(chain)

  # collect intersection pieces and
  # classify final intersection type
  onlytouching = true
  onlyedgetouching = true
  hasintersection = false
  pieces = Geometry{M,C}[]
  for seg in segments(chain)
    I = intersection(seg, poly)
    type(I) == NotIntersecting && continue
    hasintersection = true
    onlytouching &= (type(I) == Touching) || (type(I) == CornerTouching)
    onlyedgetouching &= (type(I) == EdgeTouching)
    geom = get(I)
    if geom isa Multi
      append!(pieces, parent(geom))
    else
      push!(pieces, geom)
    end
  end
  unique!(pieces)

  if !hasintersection
    return @IT NotIntersecting nothing f
  elseif onlytouching
    return @IT Touching maybemulti(pieces) f
  elseif onlyedgetouching
    return @IT EdgeTouching maybemulti(pieces) f
  else
    return @IT Intersecting maybemulti(pieces) f
  end
end
