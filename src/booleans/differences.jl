# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

getdifference(box₁::Box, box₂::Box) = getdifference(aspoly(box₁), aspoly(box₂))

getdifference(box::Box, ball::Ball) = getdifference(aspoly(box), aspoly(ball))

getdifference(ball::Ball, box::Box) = getdifference(aspoly(ball), aspoly(box))

getdifference(box::Box, poly::Polygon) = getdifference(aspoly(box), poly)

getdifference(poly::Polygon, box::Box) = getdifference(poly, aspoly(box))

getdifference(ball₁::Ball, ball₂::Ball) = getdifference(aspoly(ball₁), aspoly(ball₂))

getdifference(ball::Ball, poly::Polygon) = getdifference(aspoly(ball), poly)

getdifference(poly::Polygon, ball::Ball) = getdifference(poly, aspoly(ball))

getdifference(poly₁::Union{Polygon,MultiPolygon}, poly₂::Union{Polygon,MultiPolygon}) =
  clip(poly₁, poly₂, GreinerHormannClipping(:-))
