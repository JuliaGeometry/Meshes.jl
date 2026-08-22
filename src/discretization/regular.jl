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

appendtopo(::Sphere{𝔼{3}}, tg) = _appendpoles(tg, 2, true)

appendtopo(::Ellipsoid, tg) = _appendpoles(tg, 2, true)

appendtopo(::Disk, tg) = _appendcenter(tg)

appendtopo(::Cylinder, tg) = _appendaxis(tg)

appendtopo(::CylinderSurface, tg) = _appendpoles(tg, 1, false)

appendtopo(::Cone, tg) = _appendaxisapex(tg)

appendtopo(::ConeSurface, tg) = _appendpoles(tg, 1, false)

appendtopo(::Frustum, tg) = _appendaxis(tg)

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
  nv = nvertices(tg)
  it = nt + 1

  # center and points along the polar axis
  c = nv + 1
  north(i) = nv + 1 + i
  south(i) = nv + 1 + (nr + 1) + i

  # wrap periodic azimuthal index
  wrap(k) = k > np ? 1 : k

  # connect hexahedra in the interior
  hexas = collect(elements(tg))

  # connect polar axis with wedges
  wedges = Connectivity{Wedge,6}[]
  for i in 1:nr, k in 1:np
    l = wrap(k + 1)
    u1, v1 = cart2corner(tg, i, 1, k), cart2corner(tg, i, 1, l)
    u2, v2 = cart2corner(tg, i + 1, 1, k), cart2corner(tg, i + 1, 1, l)
    push!(wedges, connect((north(i), u1, v1, north(i + 1), u2, v2), Wedge))
    u1, v1 = cart2corner(tg, i, it, k), cart2corner(tg, i, it, l)
    u2, v2 = cart2corner(tg, i + 1, it, k), cart2corner(tg, i + 1, it, l)
    push!(wedges, connect((south(i), v1, u1, south(i + 1), v2, u2), Wedge))
  end

  # connect center with pyramids
  pyramids = Connectivity{Pyramid,5}[]
  for j in 1:nt, k in 1:np
    l = wrap(k + 1)
    u1, v1 = cart2corner(tg, 1, j, k), cart2corner(tg, 1, j, l)
    u2, v2 = cart2corner(tg, 1, j + 1, k), cart2corner(tg, 1, j + 1, l)
    push!(pyramids, connect((u1, v1, v2, u2, c), Pyramid))
  end

  # connect center and poles with tetrahedra
  tetras = Connectivity{Tetrahedron,4}[]
  for k in 1:np
    l = wrap(k + 1)
    u1, v1 = cart2corner(tg, 1, 1, k), cart2corner(tg, 1, 1, l)
    u2, v2 = cart2corner(tg, 1, it, k), cart2corner(tg, 1, it, l)
    push!(tetras, connect((c, north(1), u1, v1), Tetrahedron))
    push!(tetras, connect((c, south(1), v2, u2), Tetrahedron))
  end

  SimpleTopology([hexas; wedges; pyramids; tetras])
end

function _appendaxisapex(tg)
  # auxiliary variables
  nr, np, nh = size(tg)
  nv = nvertices(tg)
  ih = nh + 1

  # points along the axis and apex
  axis(k) = nv + k
  apex = nv + ih + 1

  # wrap periodic azimuthal index
  wrap(j) = j > np ? 1 : j

  # connect hexahedra in the volume
  hexas = collect(elements(tg))

  # connect axis with wedges
  wedges = Connectivity{Wedge,6}[]
  for k in 1:nh, j in 1:np
    l = wrap(j + 1)
    u1, v1 = cart2corner(tg, 1, j, k), cart2corner(tg, 1, l, k)
    u2, v2 = cart2corner(tg, 1, j, k + 1), cart2corner(tg, 1, l, k + 1)
    push!(wedges, connect((axis(k), u1, v1, axis(k + 1), u2, v2), Wedge))
  end

  # connect apex with pyramids
  pyramids = Connectivity{Pyramid,5}[]
  for i in 1:nr, j in 1:np
    l = wrap(j + 1)
    u1, v1 = cart2corner(tg, i, j, ih), cart2corner(tg, i, l, ih)
    u2, v2 = cart2corner(tg, i + 1, j, ih), cart2corner(tg, i + 1, l, ih)
    push!(pyramids, connect((u1, u2, v2, v1, apex), Pyramid))
  end

  # connect axis and apex with tetrahedra
  tetras = Connectivity{Tetrahedron,4}[]
  for j in 1:np
    l = wrap(j + 1)
    u, v = cart2corner(tg, 1, j, ih), cart2corner(tg, 1, l, ih)
    push!(tetras, connect((apex, axis(ih), v, u), Tetrahedron))
  end

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
