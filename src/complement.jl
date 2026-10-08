# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    !(geometry)

Return the complement of the `geometry` with
respect to its bounding box.
"""
!(g::Geometry) = _complement(_boxboundary(g), boundary(g))

function _boxboundary(g)
  ℒ = lentype(g)
  b = boundingbox(g)
  c = to(centroid(b))
  l = sides(b)
  α = (l .+ 2atol(ℒ)) ./ l
  t = Translate(-c...) → Scale(α) → Translate(c...)
  boundary(t(b))
end

_complement(b, r::Ring) = Pgon([b, reverse(r)])

function _complement(b, m::MultiRing)
  rings = parent(m)

  outer = Pgon([b, reverse(first(rings))])
  inners = [Pgon(reverse(rings[i])) for i in 2:length(rings)]

  Multi([[outer]; inners])
end

function _complement(b, s::Sphere{𝔼{2}})
  points = vertices(discretize(s))
  _complement(b, Ring(points))
end
