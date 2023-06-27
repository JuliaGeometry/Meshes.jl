# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

intersection(f, p₁::Point, p₂::Point) = p₁ == p₂ ? (@IT Intersecting p₁ f) : (@IT NotIntersecting nothing f)

intersection(f, p::Point, g::Geometry) = p ∈ g ? (@IT Intersecting p f) : (@IT NotIntersecting nothing f)
