# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

intersection(f, p₁::Point, p₂::Point) = p₁ == p₂ ? (@IT IntersectingPoint p₁ f) : (@IT NoIntersection nothing f)

intersection(f, p::Point, g::Geometry) = p ∈ g ? (@IT IntersectingPoint p f) : (@IT NoIntersection nothing f)
