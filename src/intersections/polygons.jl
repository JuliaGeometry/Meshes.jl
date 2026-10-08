# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, poly₁::Polygon, poly₂::Polygon)
  clipped = clip(poly₁, poly₂, GreinerHormannClipping())

  if isnothing(clipped)
    @IT NotIntersecting nothing f
  else
    @IT Intersecting clipped f
  end
end

intersection(f, poly::Polygon, box::Box) = intersection(f, poly, convert(Quadrangle, box))
