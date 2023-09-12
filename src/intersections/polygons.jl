# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, poly::Polygon, geom::Geometry)
  # TODO: use Weiler-Atherton or other more general clipping method
  clipped = clip(poly, geom, SutherlandHodgman())
  if isnothing(clipped)
    @IT NotIntersecting nothing f
  else
    @IT Intersecting clipped f
  end
end
