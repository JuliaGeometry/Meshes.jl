# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularDiscretization(n1, n2, ..., np)

A method to discretize primitive geometries with
`n1×n2×...×np` elements sampled regularly along
each parametric dimensions. The adequate number
of points is calculated for each type of geometry
and passed to [`RegularSampling`](@ref).
"""
struct RegularDiscretization{N} <: DiscretizationMethod
  sizes::Dims{N}
end

RegularDiscretization(sizes::Vararg{Int,N}) where {N} = RegularDiscretization(sizes)

function discretize(geometry::Geometry, method::RegularDiscretization)
  if isparametrized(geometry)
    verts, tgrid = wrapgrid(geometry, method)
    tmesh = appendtopo(geometry, tgrid)
    SimpleMesh(collect(verts), tmesh)
  else
    box = boundingbox(geometry)
    grid = discretize(box, method)
    view(grid, geometry)
  end
end

# ----------------------
# wrap grid on geometry
# ----------------------

function wrapgrid(g, m)
  sz = fitdims(m.sizes, paramdim(g))
  pd = perdims(g)
  np = @. sz + !pd
  ps = sample(g, RegularSampling(np))
  tg = GridTopology(sz, pd)
  ps, tg
end

perdims(g) = isperiodic(g)
perdims(::Sphere{3}) = (false, true)

# ------------------------
# append to grid topology
# ------------------------

appendtopo(g, tg) = tg

appendtopo(::Ball{2}, tg) = _appendcenter(tg)

appendtopo(::Disk, tg) = _appendcenter(tg)

function appendtopo(::Sphere{3}, tg)
  nx, ny = size(tg)

  # collect quadrangles in the middle
  middle = collect(elements(tg))

  # connect north pole with triangles
  u = nvertices(tg) + 1
  north = map(1:(ny - 1)) do j
    v = cart2corner(tg, 1, j)
    w = cart2corner(tg, 1, j + 1)
    connect((u, v, w))
  end
  v = cart2corner(tg, 1, ny)
  w = cart2corner(tg, 1, 1)
  push!(north, connect((u, v, w)))

  # connect south pole with triangles
  u = nvertices(tg) + 2
  south = map(1:(ny - 1)) do j
    v = cart2corner(tg, nx + 1, j)
    w = cart2corner(tg, nx + 1, j + 1)
    connect((u, w, v))
  end
  v = cart2corner(tg, nx + 1, ny)
  w = cart2corner(tg, nx + 1, 1)
  push!(south, connect((u, w, v)))

  SimpleTopology([middle; north; south])
end

appendtopo(::CylinderSurface, tg) = _appendnorthsouth(tg)

appendtopo(::ConeSurface, tg) = _appendnorthsouth(tg)

function _appendcenter(tg)
  _, ny = size(tg)

  # connect quadrangles in the middle
  quads = collect(elements(tg))

  # connect center with triangles
  u = nvertices(tg) + 1
  tris = map(1:(ny - 1)) do j
    v = cart2corner(tg, 1, j)
    w = cart2corner(tg, 1, j + 1)
    connect((u, v, w))
  end
  v = cart2corner(tg, 1, ny)
  w = cart2corner(tg, 1, 1)
  push!(tris, connect((u, v, w)))

  SimpleTopology([quads; tris])
end

function _appendnorthsouth(tg)
  nx, ny = size(tg)

  # connect quadrangles in the middle
  middle = collect(elements(tg))

  # connect south pole with triangles
  u = nvertices(tg) + 1
  south = map(1:(nx - 1)) do i
    v = cart2corner(tg, i + 1, 1)
    w = cart2corner(tg, i, 1)
    connect((u, v, w))
  end
  v = cart2corner(tg, 1, 1)
  w = cart2corner(tg, nx, 1)
  push!(south, connect((u, v, w)))

  # connect north pole with triangles
  u = nvertices(tg) + 2
  north = map(1:(nx - 1)) do i
    v = cart2corner(tg, i + 1, ny + 1)
    w = cart2corner(tg, i, ny + 1)
    connect((u, w, v))
  end
  v = cart2corner(tg, 1, ny + 1)
  w = cart2corner(tg, nx, ny + 1)
  push!(north, connect((u, w, v)))

  SimpleTopology([middle; north; south])
end

# --------------
# SPECIAL CASES
# --------------

function discretize(box::Box, method::RegularDiscretization)
  sz = fitdims(method.sizes, paramdim(box))
  CartesianGrid(extrema(box)..., dims=sz)
end
