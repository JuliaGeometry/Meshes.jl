# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, poly₁::Polygon, poly₂::Polygon)
  # TODO: use Weiler-Atherton or other more general clipping method
  clipped = clip(poly₁, poly₂, SutherlandHodgman())
  if isnothing(clipped)
    @IT NotIntersecting nothing f
  else
    @IT Intersecting clipped f
  end
end

intersection(f, poly::Polygon{2}, box::Box{2}) = intersection(f, poly, convert(Quadrangle, box))
