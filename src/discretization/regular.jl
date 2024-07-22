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
  pd = isperiodic(g)
  np = @. sz + !pd
  ps = sample(g, RegularSampling(np))
  tg = GridTopology(sz, pd)
  ps, tg
end

# ------------------------
# append to grid topology
# ------------------------

appendtopo(g, tg) = _appendtopo(g, Val(embeddim(g)), tg)

_appendtopo(g, ::Val, tp) = tp

_appendtopo(::Ball, ::Val{2}, tg) = _appendcenter(tg)

appendtopo(::Disk, tg) = _appendcenter(tg)

_appendtopo(::Sphere, ::Val{3}, tg) = _appendpoles(tg, 2, true)

appendtopo(::Ellipsoid, tg) = _appendpoles(tg, 2, true)

appendtopo(::Cylinder, tg) = _appendaxis(tg)

appendtopo(::CylinderSurface, tg) = _appendpoles(tg, 1, false)

appendtopo(::ConeSurface, tg) = _appendpoles(tg, 1, false)

appendtopo(::FrustumSurface, tg) = _appendpoles(tg, 1, false)

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

function _appendaxis(tg)
  # auxiliary variables
  _, ny, nz = size(tg)

  # number of grid vertices
  nvert = nvertices(tg)

  # connect hexahedra in the volume
  hexas = collect(elements(tg))

  # connect axis with wedges
  inds = NTuple{6,Int}[]
  for k in 1:nz
    for j in 1:(ny - 1)
      a1 = nvert + k
      b1 = cart2corner(tg, 1, j, k)
      c1 = cart2corner(tg, 1, j + 1, k)
      a2 = nvert + k + 1
      b2 = cart2corner(tg, 1, j, k + 1)
      c2 = cart2corner(tg, 1, j + 1, k + 1)
      push!(inds, (a1, b1, c1, a2, b2, c2))
    end
    a1 = nvert + k
    b1 = cart2corner(tg, 1, ny, k)
    c1 = cart2corner(tg, 1, 1, k)
    a2 = nvert + k + 1
    b2 = cart2corner(tg, 1, ny, k + 1)
    c2 = cart2corner(tg, 1, 1, k + 1)
    push!(inds, (a1, b1, c1, a2, b2, c2))
  end
  wedges = [connect(ind, Wedge) for ind in inds]

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

# --------------
# SPECIAL CASES
# --------------

function discretize(box::Box, method::RegularDiscretization)
  sz = fitdims(method.sizes, paramdim(box))
  CartesianGrid(extrema(box)..., dims=sz)
end
