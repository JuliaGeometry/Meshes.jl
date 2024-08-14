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

_complement(b, r::Ring) = PolyArea([b, reverse(r)])

function _complement(b, m::MultiRing)
  rings = parent(m)

  outer = PolyArea([b, reverse(first(rings))])
  inners = [PolyArea(reverse(rings[i])) for i in 2:length(rings)]

  Multi([[outer]; inners])
end

function _complement(b, p::Primitive)
  ring = Ring(pointify(p))
  PolyArea([b, reverse(ring)])
end
