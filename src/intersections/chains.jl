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
  # assess intersections 
  pieces = Geometry{𝔼{2},crs(chain)}[]
  hasintersection = false
  onlytouching = true
  onlyedgetouching = true
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
  glue(pieces) = length(pieces) == 1 ? only(pieces) : Multi(pieces)

  # classifying intersections
  if !hasintersection
    return @IT NotIntersecting nothing f
  elseif onlytouching
    return @IT Touching glue(pieces) f
  elseif onlyedgetouching
    return @IT EdgeTouching glue(pieces) f
  else
    return @IT Intersecting glue(pieces) f
  end
end
