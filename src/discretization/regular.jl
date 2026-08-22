# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularDiscretization(n1, n2, ..., np)

A method to discretize primitive geometries with
`n1×n2×...×np` elements sampled regularly along
each parametric dimension. The adequate number of
points is calculated for each type of geometry and
forwarded to [`RegularSampling`](@ref).
"""
struct RegularDiscretization{N} <: DiscretizationMethod
  sizes::Dims{N}
end

RegularDiscretization(sizes::Vararg{Int,N}) where {N} = RegularDiscretization(sizes)

function discretize(geom::Geometry, method::RegularDiscretization)
  if isparametrized(geom)
    verts, tgrid = wrapgrid(geom, method)
    tmesh = appendtopo(geom, tgrid)
    SimpleMesh(collect(verts), tmesh)
  else
    _discretizewithinbox(geom, method)
  end
end

# box is trivial to discretize into a regular grid
function discretize(box::Box, method::RegularDiscretization)
  sz = fitdims(method.sizes, paramdim(box))
  RegularGrid(extrema(box)..., dims=sz)
end

# triangle is parametrized with barycentric coordinates, so we can't rely on regular sampling
discretize(tri::Triangle, method::RegularDiscretization) = _discretizewithinbox(tri, method)

# tetrahedron is parametrized with barycentric coordinates, so we can't rely on regular sampling
discretize(tetra::Tetrahedron, method::RegularDiscretization) = _discretizewithinbox(tetra, method)

function _discretizewithinbox(geom::Geometry, method::RegularDiscretization)
  box = boundingbox(geom)
  grid = discretize(box, method)
  view(grid, geom)
end

# ----------------------
# wrap grid on geometry
# ----------------------

function wrapgrid(g, m)
  sz = fitdims(m.sizes, paramdim(g))
  pd = isperiodic(g)
  np = @. sz + !pd
  ps = sample(g, RegularSampling(np))
  tg = GridTopology(sz, pd)
  ps, tg
end

# ------------------------
# append to grid topology
# ------------------------

appendtopo(g, tg) = tg

appendtopo(::Ball{𝔼{2}}, tg) = _appendcenter(tg)

appendtopo(::Ball{𝔼{3}}, tg) = _appendcenteraxis(tg)

appendtopo(::Disk, tg) = _appendcenter(tg)

appendtopo(::Sphere{𝔼{3}}, tg) = _appendpoles(tg, 2, true)

appendtopo(::Ellipsoid, tg) = _appendpoles(tg, 2, true)

appendtopo(::Cylinder, tg) = _appendaxis(tg)

appendtopo(::CylinderSurface, tg) = _appendpoles(tg, 1, false)

appendtopo(::ConeSurface, tg) = _appendpoles(tg, 1, false)

appendtopo(::FrustumSurface, tg) = _appendpoles(tg, 1, false)

appendtopo(::ParaboloidSurface, tg) = _appendcenter(tg)

function _appendcenter(tg)
  # auxiliary variables
  _, ny = size(tg)

  # center of disk
  c = nvertices(tg) + 1

  # connect quadrangles in the disk
  quads = collect(elements(tg))

  # connect center with triangles
  tris = map(1:(ny - 1)) do j
    u = cart2corner(tg, 1, j)
    v = cart2corner(tg, 1, j + 1)
    connect((c, u, v))
  end
  u = cart2corner(tg, 1, ny)
  v = cart2corner(tg, 1, 1)
  push!(tris, connect((c, u, v)))

  SimpleTopology([quads; tris])
end

function _appendcenteraxis(tg)
  # auxiliary variables
  nr, nt, np = size(tg)
  it, nvert = nt + 1, nvertices(tg)

  # center and points along the polar axis
  c = nvert + 1
  north(i) = nvert + 1 + i
  south(i) = nvert + 1 + (nr + 1) + i

  # wrap periodic azimuthal index
  wrap(k) = k > np ? 1 : k

  # connect hexahedra in the interior
  hexas = collect(elements(tg))

  # connect polar axis with wedges
  winds = NTuple{6,Int}[]
  for i in 1:nr, k in 1:np
    l = wrap(k + 1)
    u1, v1 = cart2corner(tg, i, 1, k), cart2corner(tg, i, 1, l)
    u2, v2 = cart2corner(tg, i + 1, 1, k), cart2corner(tg, i + 1, 1, l)
    push!(winds, (north(i), u1, v1, north(i + 1), u2, v2))
    u1, v1 = cart2corner(tg, i, it, k), cart2corner(tg, i, it, l)
    u2, v2 = cart2corner(tg, i + 1, it, k), cart2corner(tg, i + 1, it, l)
    push!(winds, (south(i), v1, u1, south(i + 1), v2, u2))
  end
  wedges = [connect(ind, Wedge) for ind in winds]

  # connect center with pyramids
  pinds = NTuple{5,Int}[]
  for j in 1:nt, k in 1:np
    l = wrap(k + 1)
    u1, v1 = cart2corner(tg, 1, j, k), cart2corner(tg, 1, j, l)
    u2, v2 = cart2corner(tg, 1, j + 1, k), cart2corner(tg, 1, j + 1, l)
    push!(pinds, (u1, v1, v2, u2, c))
  end
  pyramids = [connect(ind, Pyramid) for ind in pinds]

  # connect center and poles with tetrahedra
  tinds = NTuple{4,Int}[]
  for k in 1:np
    l = wrap(k + 1)
    push!(tinds, (c, north(1), cart2corner(tg, 1, 1, k), cart2corner(tg, 1, 1, l)))
    push!(tinds, (c, south(1), cart2corner(tg, 1, it, l), cart2corner(tg, 1, it, k)))
  end
  tetras = [connect(ind, Tetrahedron) for ind in tinds]

  SimpleTopology([hexas; wedges; pyramids; tetras])
end

function _appendaxis(tg)
  # auxiliary variables
  _, ny, nz = size(tg)

  # number of grid vertices
  nv = nvertices(tg)

  # connect hexahedra in the volume
  hexas = collect(elements(tg))

  # connect axis with wedges
  wedges = Connectivity{Wedge,6}[]
  for k in 1:nz
    for j in 1:(ny - 1)
      a1 = nv + k
      b1 = cart2corner(tg, 1, j, k)
      c1 = cart2corner(tg, 1, j + 1, k)
      a2 = nv + k + 1
      b2 = cart2corner(tg, 1, j, k + 1)
      c2 = cart2corner(tg, 1, j + 1, k + 1)
      push!(wedges, connect((a1, b1, c1, a2, b2, c2), Wedge))
    end
    a1 = nv + k
    b1 = cart2corner(tg, 1, ny, k)
    c1 = cart2corner(tg, 1, 1, k)
    a2 = nv + k + 1
    b2 = cart2corner(tg, 1, ny, k + 1)
    c2 = cart2corner(tg, 1, 1, k + 1)
    push!(wedges, connect((a1, b1, c1, a2, b2, c2), Wedge))
  end

  SimpleTopology([hexas; wedges])
end

# connect north and south poles to
# grid topology along given dimension
# and counter-clockwise orientation
function _appendpoles(tg, d, ccw)
  # auxiliary variables
  sz = size(tg)
  nd = length(sz)

  # swap indices of poles if necessary
  swap(u, v) = ccw ? (u, v) : (v, u)

  # north and south poles
  n = nvertices(tg) + 1
  s = nvertices(tg) + 2

  # connect quadrangles in the trunk
  trunk = collect(elements(tg))

  # connect north pole with triangles
  north = map(1:(sz[d] - 1)) do j
    iᵤ = ntuple(i -> i == d ? j : 1, nd)
    iᵥ = ntuple(i -> i == d ? j + 1 : 1, nd)
    u = cart2corner(tg, iᵤ...)
    v = cart2corner(tg, iᵥ...)
    connect((n, swap(u, v)...))
  end
  iᵤ = ntuple(i -> i == d ? sz[d] : 1, nd)
  iᵥ = ntuple(i -> 1, nd)
  u = cart2corner(tg, iᵤ...)
  v = cart2corner(tg, iᵥ...)
  push!(north, connect((n, swap(u, v)...)))

  # connect south pole with triangles
  south = map(1:(sz[d] - 1)) do j
    iᵤ = ntuple(i -> i == d ? j : sz[i] + 1, nd)
    iᵥ = ntuple(i -> i == d ? j + 1 : sz[i] + 1, nd)
    u = cart2corner(tg, iᵤ...)
    v = cart2corner(tg, iᵥ...)
    connect((s, swap(v, u)...))
  end
  iᵤ = ntuple(i -> i == d ? sz[d] : sz[i] + 1, nd)
  iᵥ = ntuple(i -> i == d ? 1 : sz[i] + 1, nd)
  u = cart2corner(tg, iᵤ...)
  v = cart2corner(tg, iᵥ...)
  push!(south, connect((s, swap(v, u)...)))

  SimpleTopology([trunk; north; south])
end
