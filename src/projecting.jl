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

function proj2D(p::PolyArea)
  r = proj2D.(rings(p))
  if hasholes(p)
    PolyArea(first(r), r[2:end])
  else
    PolyArea(first(r))
  end
end

# ---------------
# IMPLEMENTATION
# ---------------

function proj2D(points::AbstractVector{Point{3,T}}) where {T}
  # https://math.stackexchange.com/a/99317
  X = mapreduce(coordinates, hcat, points)
  μ = sum(X, dims=2) / size(X, 2)
  Z = X .- μ
  U = svd(Z).U
  u = U[:, 1]
  v = U[:, 2]
  n = T[0, 0, 1]
  if (u × v) ⋅ n < 0
    u, v = v, u
  end
  [Point(z ⋅ u, z ⋅ v) for z in eachcol(Z)]
end
