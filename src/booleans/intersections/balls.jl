# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

intersection(f, ball₁::Ball, ball₂::Ball) = intersection(f, aspoly(ball₁), aspoly(ball₂))

intersection(f, ball::Ball, poly::Polygon) = intersection(f, aspoly(ball), poly)
