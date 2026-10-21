# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

getunion(g₁, g₂) = getunion(g₂, g₁) # exploit symmetry

getunion(box₁::Box, box₂::Box) = _polygonunion(box₁, box₂)

getunion(box::Box, ball::Ball) = _polygonunion(box, ball)

getunion(box::Box, poly::Polygon) = _polygonunion(box, poly)

getunion(ball₁::Ball, ball₂::Ball) = _polygonunion(ball₁, ball₂)

getunion(ball::Ball, poly::Polygon) = _polygonunion(ball, poly)

getunion(poly₁::Union{Polygon,MultiPolygon}, poly₂::Union{Polygon,MultiPolygon}) = clip(poly₁, poly₂, GreinerHormannClipping(:∪))

# -----------------
# HELPER FUNCTIONS
# -----------------

function _polygonunion(geom₁, geom₂)
  mesh₁ = discretize(boundary(geom₁))
  mesh₂ = discretize(boundary(geom₂))
  pgon₁ = Pgon(vertices(mesh₁))
  pgon₂ = Pgon(vertices(mesh₂))
  getunion(pgon₁, pgon₂)
end
