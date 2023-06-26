# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

intersection(f, p::Point, g::Geometry) =
  p ∈ g ? (@IT IntersectingPoint p f) : (@IT NoIntersection nothing f)

intersection(f, p::Point, g::Polygon) =
  p ∈ g ? (@IT IntersectingPoint p f) : (@IT NoIntersection nothing f)
