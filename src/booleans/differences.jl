# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

getdifference(box₁::Box, box₂::Box) = _polygondifference(box₁, box₂)

getdifference(box::Box, ball::Ball) = _polygondifference(box, ball)

getdifference(ball::Ball, box::Box) = _polygondifference(ball, box)

getdifference(box::Box, poly::Polygon) = _polygondifference(box, poly)

getdifference(poly::Polygon, box::Box) = _polygondifference(poly, box)

getdifference(ball₁::Ball, ball₂::Ball) = _polygondifference(ball₁, ball₂)

getdifference(ball::Ball, poly::Polygon) = _polygondifference(ball, poly)

getdifference(poly::Polygon, ball::Ball) = _polygondifference(poly, ball)

getdifference(poly₁::Union{Polygon,MultiPolygon}, poly₂::Union{Polygon,MultiPolygon}) =
  clip(poly₁, poly₂, GreinerHormannClipping(:-))

# -----------------
# HELPER FUNCTIONS
# -----------------

function _polygondifference(geom₁, geom₂)
  mesh₁ = discretize(boundary(geom₁))
  mesh₂ = discretize(boundary(geom₂))
  pgon₁ = Pgon(vertices(mesh₁))
  pgon₂ = Pgon(vertices(mesh₂))
  getdifference(pgon₁, pgon₂)
end
