# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

getunion(g₁, g₂) = getunion(g₂, g₁) # exploit symmetry

getunion(box₁::Box, box₂::Box) = getunion(aspoly(box₁), aspoly(box₂))

getunion(box::Box, ball::Ball) = getunion(aspoly(box), aspoly(ball))

getunion(box::Box, poly::Union{Polygon,MultiPolygon}) = getunion(aspoly(box), poly)

getunion(ball₁::Ball, ball₂::Ball) = getunion(aspoly(ball₁), aspoly(ball₂))

getunion(ball::Ball, poly::Union{Polygon,MultiPolygon}) = getunion(aspoly(ball), poly)

getunion(poly₁::Union{Polygon,MultiPolygon}, poly₂::Union{Polygon,MultiPolygon}) =
  clip(poly₁, poly₂, GreinerHormannClipping(:∪))
