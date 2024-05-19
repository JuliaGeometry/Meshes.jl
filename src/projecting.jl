# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    proj2D(geometry)

Project 3D `geometry` onto a 2D plane of
maximum variance using singular values.
"""
function proj2D end

proj2D(r::Rope) = Ring(proj2D(vertices(r)))

proj2D(r::Ring) = Ring(proj2D(vertices(r)))

proj2D(p::Ngon) = Ngon(proj2D(collect(vertices(p)))...)

proj2D(p::PolyArea) = PolyArea(proj2D.(rings(p)))

# ---------------
# IMPLEMENTATION
# ---------------

proj2D(points::AbstractVector{<:Point{3}}) = proj(points, svdbasis(points))

function proj(points, basis)
  # retrieve basis
  u, v = basis

  # centroid of projection
  c = centroid(PointSet(points))

  # project points
  map(points) do p
    d = p - c
    x = udot(d, u)
    y = udot(d, v)
    Point(x, y)
  end
end
