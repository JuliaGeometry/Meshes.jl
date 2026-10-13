# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, multi::Multi, geom::Geometry)
  # intersect each component geometry
  geoms = parent(multi)
  inter = map(geoms) do geomᵢ
    intersection(f, geomᵢ, geom)
  end

  # filter out non-intersecting results
  valid = filter(i -> type(i) ≠ NotIntersecting, inter)

  # return valid intersection when it exists
  if isempty(valid)
    return @IT NotIntersecting nothing f
  else
    m = maybemulti(map(get, valid))
    return @IT Intersecting m f
  end
end
