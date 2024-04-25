# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

intersection(f, point₁::Point, point₂::Point) =
  point₁ == point₂ ? (@IT Intersecting point₁ f) : (@IT NotIntersecting nothing f)

@commutativef intersection(f, point::Point, geom::Geometry) =
  point ∈ geom ? (@IT Intersecting point f) : (@IT NotIntersecting nothing f)
