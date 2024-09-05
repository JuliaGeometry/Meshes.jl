# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    signarea(A, B, C)

Compute signed area of triangle formed by points `A`, `B` and `C`.
"""
function signarea(A::Point, B::Point, C::Point)
  checkdim(A, 2)
  ((B - A) × (C - A)) / 2
end

"""
    householderbasis(n)

Returns a pair of orthonormal tangent vectors `u` and `v` from a normal `n`,
such that `u`, `v`, and `n` form a right-hand orthogonal system.

## References

* D.S. Lopes et al. 2013. ["Tangent vectors to a 3-D surface normal: A geometric tool
  to find orthogonal vectors based on the Householder transformation"]
  (https://doi.org/10.1016/j.cad.2012.11.003)
"""
function householderbasis(n::Vec{3,ℒ}) where {ℒ}
  n̂ = norm(n)
  i = argmax(n .+ n̂)
  n̂ᵢ = Vec(ntuple(j -> j == i ? n̂ : zero(ℒ), 3))
  h = n + n̂ᵢ
  H = (I - 2h * transpose(h) / (transpose(h) * h)) * unit(ℒ)
  u, v = [H[:, j] for j in 1:3 if j != i]
  i == 2 && ((u, v) = (v, u))
  Vec(u), Vec(v)
end

"""
    svdbasis(points)

Returns the 2D basis that retains most of the variance in the list of 3D `points`
using the singular value decomposition (SVD).

See <https://math.stackexchange.com/a/99317>.
"""
function svdbasis(p::AbstractVector{<:Point})
  checkdim(first(p), 3)
  ℒ = lentype(eltype(p))
  X = reduce(hcat, to.(p))
  μ = sum(X, dims=2) / size(X, 2)
  Z = X .- μ
  U = usvd(Z).U
  u = Vec(U[:, 1]...)
  v = Vec(U[:, 2]...)
  n = Vec(zero(ℒ), zero(ℒ), oneunit(ℒ))
  isnegative((u × v) ⋅ n) ? (v, u) : (u, v)
end

"""
    intersectparameters(a, b, c, d)

Compute the parameters `λ₁` and `λ₂` of the lines 
`a + λ₁ ⋅ v⃗₁`, with `v⃗₁ = b - a` and
`c + λ₂ ⋅ v⃗₂`, with `v⃗₂ = d - c` spanned by the input
points `a`, `b` resp. `c`, `d` such that to yield line
points with minimal distance or the intersection point
(if lines intersect).

Furthermore, the ranks `r` of the matrix of the linear
system `A ⋅ λ⃗ = y⃗`, with `A = [v⃗₁ -v⃗₂], y⃗ = c - a`
and the rank `rₐ` of the augmented matrix `[A y⃗]` are
calculated in order to identify the intersection type:

- Intersection: r == rₐ == 2
- Colinear: r == rₐ == 1
- No intersection: r != rₐ
  - No intersection and parallel:  r == 1, rₐ == 2
  - No intersection, skew lines: r == 2, rₐ == 3
"""
function intersectparameters(a::Point, b::Point, c::Point, d::Point)
  A = ustrip.([(b - a) (c - d)])
  y = ustrip.(c - a)
  T = eltype(A)

  # calculate the rank of the augmented matrix by checking
  # the zero entries of the diagonal of R
  _, R = qr([A y])

  # for Dim == 2 one has to check the L1 norm of rows as 
  # there are more columns than rows
  τ = atol(T)
  rₐ = sum(>(τ), sum(abs, R, dims=2))

  # calculate the rank of the rectangular matrix
  r = sum(>(τ), sum(abs, view(R, :, 1:2), dims=2))

  # calculate parameters of intersection or closest point
  if r ≥ 2
    λ = A \ y
    λ₁, λ₂ = λ[1], λ[2]
  else # parallel or collinear
    λ₁, λ₂ = zero(T), zero(T)
  end

  λ₁, λ₂, r, rₐ
end

function cartesianrange(grid::CartesianGrid, limits)
  # grid properties
  or = minimum(grid)
  sp = spacing(grid)
  sz = size(grid)
  nd = length(sz)

  # box points
  bmin = withcrs(grid, ntuple(i -> first(limits[i]), nd))
  bmax = withcrs(grid, ntuple(i -> last(limits[i]), nd))

  # intersection of boxes
  lo, up = extrema(boundingbox(grid) ∩ Box(bmin, bmax))

  # Cartesian indices of new corners
  ilo = max.(ceil.(Int, (lo - or) ./ sp), 1)
  iup = min.(floor.(Int, (up - or) ./ sp) .+ 1, sz)

  # Cartesian range from corner to corner
  CartesianIndex(Tuple(ilo)):CartesianIndex(Tuple(iup))
end

function cartesianrange(grid::RectilinearGrid, limits)
  # grid properties
  sz = size(grid)
  nd = length(sz)

  # box points
  bmin = withcrs(grid, ntuple(i -> first(limits[i]), nd))
  bmax = withcrs(grid, ntuple(i -> last(limits[i]), nd))

  # intersection of boxes
  lo, up = to.(extrema(boundingbox(grid) ∩ Box(bmin, bmax)))

  # integer coordinates of lower point
  ilo = ntuple(nd) do i
    findlast(x -> x ≤ lo[i], xyz(grid)[i])
  end

  # integer coordinates of upper point
  iup = ntuple(nd) do i
    findfirst(x -> x ≥ up[i], xyz(grid)[i])
  end

  # integer coordinates of elements
  CartesianIndex(ilo):CartesianIndex(iup .- 1)
end

function cartesianrange(grid::Grid{𝔼{2}}, limits)
  nx, ny = vsize(grid)
  (lxmin, lxmax), (lymin, lymax) = limits

  a = convert(Cartesian, coords(vertex(grid, (1, 1))))
  b = convert(Cartesian, coords(vertex(grid, (nx, 1))))
  c = convert(Cartesian, coords(vertex(grid, (1, ny))))

  xmin = max(lxmin, a.x)
  ymin = max(lymin, a.y)
  xmax = min(lxmax, b.x)
  ymax = min(lymax, c.y)

  iₛ = findlast(1:nx) do i
    p = vertex(grid, (i, 1))
    c = convert(Cartesian, coords(p))
    c.x ≤ xmin
  end
  iₑ = findfirst(1:nx) do i
    p = vertex(grid, (i, 1))
    c = convert(Cartesian, coords(p))
    c.x ≥ xmax
  end
  jₛ = findlast(1:ny) do i
    p = vertex(grid, (1, i))
    c = convert(Cartesian, coords(p))
    c.y ≤ ymin
  end
  jₑ = findfirst(1:ny) do i
    p = vertex(grid, (1, i))
    c = convert(Cartesian, coords(p))
    c.y ≥ ymax
  end
  if iₛ == iₑ || jₛ == jₑ
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  CartesianIndex(iₛ, jₛ):CartesianIndex(iₑ - 1, jₑ - 1)
end

function cartesianrange(grid::Grid{𝔼{3}}, limits)
  nx, ny, nz = vsize(grid)
  (lxmin, lxmax), (lymin, lymax), (lzmin, lzmax) = limits

  a = convert(Cartesian, coords(vertex(grid, (1, 1, 1))))
  b = convert(Cartesian, coords(vertex(grid, (nx, 1, 1))))
  c = convert(Cartesian, coords(vertex(grid, (1, ny, 1))))
  d = convert(Cartesian, coords(vertex(grid, (1, 1, nz))))

  xmin = max(lxmin, a.x)
  ymin = max(lymin, a.y)
  zmin = max(lzmin, a.z)
  xmax = min(lxmax, b.x)
  ymax = min(lymax, c.y)
  zmax = min(lzmax, d.z)

  iₛ = findlast(1:nx) do i
    p = vertex(grid, (i, 1, 1))
    c = convert(Cartesian, coords(p))
    c.x ≤ xmin
  end
  iₑ = findfirst(1:nx) do i
    p = vertex(grid, (i, 1, 1))
    c = convert(Cartesian, coords(p))
    c.x ≥ xmax
  end
  jₛ = findlast(1:ny) do i
    p = vertex(grid, (1, i, 1))
    c = convert(Cartesian, coords(p))
    c.y ≤ ymin
  end
  jₑ = findfirst(1:ny) do i
    p = vertex(grid, (1, i, 1))
    c = convert(Cartesian, coords(p))
    c.y ≥ ymax
  end
  kₛ = findlast(1:nz) do i
    p = vertex(grid, (1, 1, i))
    c = convert(Cartesian, coords(p))
    c.z ≤ zmin
  end
  kₑ = findfirst(1:nz) do i
    p = vertex(grid, (1, 1, i))
    c = convert(Cartesian, coords(p))
    c.z ≥ zmax
  end

  if iₛ == iₑ || jₛ == jₑ || kₛ == kₑ
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  CartesianIndex(iₛ, jₛ, kₛ):CartesianIndex(iₑ - 1, jₑ - 1, kₑ - 1)
end

function cartesianrange(grid::Grid{🌐}, limits)  
  nlon, nlat = vsize(grid)
  (llonmin, llonmax), (llatmin, llatmax) = limits

  a = convert(Cartesian, coords(vertex(grid, (1, 1))))
  b = convert(Cartesian, coords(vertex(grid, (nlon, 1))))
  c = convert(Cartesian, coords(vertex(grid, (1, nlat))))

  lonmin = max(llonmin, a.lon)
  latmin = max(llatmin, a.lat)
  lonmax = min(llonmax, b.lon)
  latmax = min(llatmax, c.lat)

  iₛ = findlast(1:nlon) do i
    p = vertex(grid, (i, 1))
    c = convert(LatLon, coords(p))
    c.lon ≤ lonmin
  end
  iₑ = findfirst(1:nlon) do i
    p = vertex(grid, (i, 1))
    c = convert(LatLon, coords(p))
    c.lon ≥ lonmax
  end
  jₛ = findlast(1:nlat) do i
    p = vertex(grid, (1, i))
    c = convert(LatLon, coords(p))
    c.lat ≤ latmin
  end
  jₑ = findfirst(1:nlat) do i
    p = vertex(grid, (1, i))
    c = convert(LatLon, coords(p))
    c.lat ≥ latmax
  end

  if iₛ == iₑ || jₛ == jₑ
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  CartesianIndex(iₛ, jₛ):CartesianIndex(iₑ - 1, jₑ - 1)
end
