# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, poly₁::Polygon, poly₂::Polygon)
  println("intersection of polygons")
  println("poly₁: ", poly₁)
  println("poly₂: ", poly₂)
  println("isconvex(poly₁): ", isconvex(poly₁))
  println("isconvex(poly₂): ", isconvex(poly₂))
  # TODO: use Weiler-Atherton or other more general clipping method
  clipped = if isconvex(poly₂)
    clip(poly₁, poly₂, SutherlandHodgmanClipping())
  elseif isconvex(poly₁)
    clip(poly₂, poly₁, SutherlandHodgmanClipping())
  else
    println("test")
    geobooleanop(poly₁, poly₂, GeoIntersection())
    # throw(ErrorException("intersection not implemented between two non-convex polygons"))
  end

  if isnothing(clipped)
    @IT NotIntersecting nothing f
  else
    @IT Intersecting clipped f
  end
end

intersection(f, poly::Polygon, box::Box) = intersection(f, poly, convert(Quadrangle, box))
