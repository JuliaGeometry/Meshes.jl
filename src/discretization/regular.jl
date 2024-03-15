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
  ip = isperiodic(g)
  np = @. sz + !ip
  ps = sample(g, RegularSampling(np))
  tg = GridTopology(sz, ip)
  ps, tg
end

function wrapgrid(g::Sphere{3}, m)
  sz = fitdims(m.sizes, paramdim(g))
  ip = (false, true)
  np = @. sz + !ip
  ps = sample(g, RegularSampling(np))
  tg = GridTopology(sz, ip)
  ps, tg
end

# ------------------------
# append to grid topology
# ------------------------

appendtopo(g, tg) = tg

appendtopo(::Ball{2}, tg) = _appendcenter(tg)

appendtopo(::Disk, tg) = _appendcenter(tg)

function _appendcenter(tg)
  sz = size(tg)
  ip = isperiodic(tg)
  np = @. sz + !ip
  nx, ny = np

  # connect quadrangles in the middle
  quads = collect(elements(tg))

  # connect center with triangles
  tris = map(1:(ny - 1)) do j
    u = nx * ny + 1
    v = 1 + (j - 1) * nx
    w = 1 + (j) * nx
    connect((u, v, w))
  end
  u = nx * ny + 1
  v = 1 + (ny - 1) * nx
  w = 1
  push!(tris, connect((u, v, w)))

  SimpleTopology([quads; tris])
end

function appendtopo(::Sphere{3}, tg)
  sz = size(tg)
  ip = isperiodic(tg)
  np = @. sz + !ip
  nx, ny = np

  # collect quadrangles in the middle
  middle = collect(elements(tg))

  # connect north pole with triangles
  north = map(1:(ny - 1)) do j
    u = nx * ny + 1
    v = 1 + (j - 1) * nx
    w = 1 + (j) * nx
    connect((u, v, w))
  end
  u = nx * ny + 1
  v = 1 + (ny - 1) * nx
  w = 1
  push!(north, connect((u, v, w)))

  # connect south pole with triangles
  south = map(1:(ny - 1)) do j
    u = nx * ny + 2
    v = (j) * nx
    w = (j + 1) * nx
    connect((u, w, v))
  end
  u = nx * ny + 2
  v = ny * nx
  w = nx
  push!(south, connect((u, w, v)))

  SimpleTopology([middle; north; south])
end

appendtopo(::CylinderSurface, tg) = _appendnorthsouth(tg)

appendtopo(::ConeSurface, tg) = _appendnorthsouth(tg)

function _appendnorthsouth(tg)
  sz = size(tg)
  ip = isperiodic(tg)
  np = @. sz + !ip
  nx, ny = np

  # connect quadrangles in the middle
  middle = collect(elements(tg))

  # connect south pole with triangles
  south = map(1:(nx - 1)) do i
    u = nx * ny + 1
    v = i + 1
    w = i
    connect((u, v, w))
  end
  u = nx * ny + 1
  v = 1
  w = nx
  push!(south, connect((u, v, w)))

  # connect north pole with triangles
  offset = nx * ny - nx
  north = map(1:(nx - 1)) do i
    u = nx * ny + 2
    v = offset + i + 1
    w = offset + i
    connect((u, w, v))
  end
  u = nx * ny + 2
  v = nx * ny - nx + 1
  w = nx * ny
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
