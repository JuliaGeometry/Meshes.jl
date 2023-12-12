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

function proj2D(points::AbstractVector{<:Point{3}})
  # retrieve coordinates
  X = reduce(hcat, coordinates.(points))
  μ = colmean(X)

  # compute SVD basis
  u, v = svdbasis(X, μ)

  # project points
  c = Point(μ...)
  map(points) do p
    d = p - c
    Point(d ⋅ u, d ⋅ v)
  end
end
