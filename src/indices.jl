# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ---------------
# LINEAR INDICES
# ---------------

"""
    indices(domain, geometry)

Return the indices of the elements of the `domain` that intersect with the `geometry`.
"""
indices(domain::Domain, geometry::Geometry) = findall(intersects(geometry), domain)

function indices(grid::OrthoRegularGrid, point::Point)
  point ∉ grid && return Int[]

  # grid properties
  orig = minimum(grid)
  spac = spacing(grid)
  dims = size(grid)

  # integer coordinates
  coords = ceil.(Int, (point - orig) ./ spac)

  # fix coordinates that are on the grid border
  coords = clamp.(coords, 1, dims)

  # convert to linear index
  [LinearIndices(dims)[coords...]]
end

function indices(grid::OrthoRegularGrid, chain::Chain)
  dims = size(grid)
  mask = falses(dims)

  for segment in segments(chain)
    p₁, p₂ = vertices(segment)
    _bresenham!(mask, grid, true, p₁, p₂)
  end

  LinearIndices(dims)[mask]
end

function indices(grid::OrthoRegularGrid, poly::Polygon)
  dims = size(grid)
  mask = zeros(Int, dims)
  cpoly = poly ∩ boundingbox(grid)
  isnothing(cpoly) && return Int[]

  for (i, triangle) in enumerate(simplexify(cpoly))
    _fill!(mask, grid, i, triangle)
  end

  LinearIndices(dims)[mask .> 0]
end

function indices(grid::OrthoRegularGrid, box::Box)
  # cartesian range
  range = cartesianrange(grid, box)

  # convert to linear indices
  LinearIndices(size(grid))[range] |> vec
end

indices(grid::OrthoRegularGrid, multi::Multi) = mapreduce(geom -> indices(grid, geom), vcat, parent(multi)) |> unique

function indices(grid::OrthoRectilinearGrid, box::Box)
  # cartesian range
  range = cartesianrange(grid, box)

  # convert to linear indices
  LinearIndices(size(grid))[range] |> vec
end

function indice(grid::TransformedGrid, geometry::Geometry)
  # construct reverse transform from revertible steps
  revtrans = reduce(→, reverse(filter(isrevertible, transform(grid))))

  # find indices in non-transformed space
  indices(parent(grid), revtrans(geometry))
end

# ----------------
# CARTESIAN RANGE
# ----------------

"""
    cartesianrange(grid, box)

Return the Cartesian range of the elements of the `grid` that intersect with the `box`.
"""
cartesianrange(grid::Grid{M}, box::Box{M}) where {M} = _manifoldrange(M, grid, box)

_manifoldrange(::Type{<:𝔼}, grid::Grid, box::Box) = _euclideanrange(grid, box)

_manifoldrange(::Type{<:🌐}, grid::Grid, box::Box) = _geodesicrange(grid, box)

function _euclideanrange(grid::OrthoRegularGrid, box::Box)
  # grid properties
  or = minimum(grid)
  sp = spacing(grid)
  sz = size(grid)

  # intersection of boxes
  lo, up = extrema(boundingbox(grid) ∩ box)

  # Cartesian indices of new corners
  ijkₛ = max.(ceil.(Int, (lo - or) ./ sp), 1)
  ijkₑ = min.(floor.(Int, (up - or) ./ sp) .+ 1, sz)

  # Cartesian range from corner to corner
  CartesianIndex(Tuple(ijkₛ)):CartesianIndex(Tuple(ijkₑ))
end

function _euclideanrange(grid::OrthoRectilinearGrid, box::Box)
  # grid properties
  nd = paramdim(grid)

  # intersection of boxes
  lo, up = to.(extrema(boundingbox(grid) ∩ box))

  # integer coordinates of lower point
  ijkₛ = ntuple(nd) do i
    findlast(x -> x ≤ lo[i], xyz(grid)[i])
  end

  # integer coordinates of upper point
  ijkₑ = ntuple(nd) do i
    findfirst(x -> x ≥ up[i], xyz(grid)[i])
  end

  # integer coordinates of elements
  CartesianIndex(ijkₛ):CartesianIndex(ijkₑ .- 1)
end

function _geodesicrange(grid::Grid, box::Box)
  nlon, nlat = vsize(grid)

  boxmin = convert(LatLon, coords(minimum(box)))
  boxmax = convert(LatLon, coords(maximum(box)))

  a = convert(LatLon, coords(vertex(grid, (1, 1))))
  b = convert(LatLon, coords(vertex(grid, (nlon, 1))))
  c = convert(LatLon, coords(vertex(grid, (1, nlat))))

  swaplon = a.lon > b.lon
  swaplat = a.lat > c.lat

  loninds = swaplon ? (nlon:-1:1) : (1:1:nlon)
  latinds = swaplat ? (nlat:-1:1) : (1:1:nlat)

  gridlonₛ, gridlonₑ = swaplon ? (b.lon, a.lon) : (a.lon, b.lon)
  gridlatₛ, gridlatₑ = swaplat ? (c.lat, a.lat) : (a.lat, c.lat)

  lonmin = max(boxmin.lon, gridlonₛ)
  latmin = max(boxmin.lat, gridlatₛ)
  lonmax = min(boxmax.lon, gridlonₑ)
  latmax = min(boxmax.lat, gridlatₑ)

  iₛ = findlast(loninds) do i
    p = vertex(grid, (i, 1))
    c = convert(LatLon, coords(p))
    c.lon ≤ lonmin
  end
  iₑ = findfirst(loninds) do i
    p = vertex(grid, (i, 1))
    c = convert(LatLon, coords(p))
    c.lon ≥ lonmax
  end

  jₛ = findlast(latinds) do i
    p = vertex(grid, (1, i))
    c = convert(LatLon, coords(p))
    c.lat ≤ latmin
  end
  jₑ = findfirst(latinds) do i
    p = vertex(grid, (1, i))
    c = convert(LatLon, coords(p))
    c.lat ≥ latmax
  end

  if iₛ == iₑ || jₛ == jₑ
    throw(ArgumentError("the passed limits are not valid for the grid"))
  end

  iₛ, iₑ = swaplon ? (iₑ, iₛ) : (iₛ, iₑ)
  jₛ, jₑ = swaplat ? (jₑ, jₛ) : (jₛ, jₑ)

  CartesianIndex(loninds[iₛ], latinds[jₛ]):CartesianIndex(loninds[iₑ] - 1, latinds[jₑ] - 1)
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _fill!(mask, grid, val, triangle)
  v = vertices(triangle)

  # fill edges of triangle
  _bresenham!(mask, grid, val, v[1], v[2])
  _bresenham!(mask, grid, val, v[2], v[3])
  _bresenham!(mask, grid, val, v[3], v[1])

  # fill interior of triangle
  j₁ = findfirst(==(val), mask).I[2]
  j₂ = findlast(==(val), mask).I[2]
  for j in j₁:j₂
    i₁ = findfirst(==(val), @view(mask[:, j]))
    i₂ = findlast(==(val), @view(mask[:, j]))
    mask[i₁:i₂, j] .= val
  end
end

# Bresenham's line algorithm: https://en.wikipedia.org/wiki/Bresenham's_line_algorithm
function _bresenham!(mask, grid, val, p₁, p₂)
  o = minimum(grid)
  s = spacing(grid)

  # integer coordinates
  x₁, y₁ = ceil.(Int, (p₁ - o) ./ s)
  x₂, y₂ = ceil.(Int, (p₂ - o) ./ s)

  # fix coordinates of points that are on the grid border
  xmax, ymax = size(grid)
  x₁ = clamp(x₁, 1, xmax)
  y₁ = clamp(y₁, 1, ymax)
  x₂ = clamp(x₂, 1, xmax)
  y₂ = clamp(y₂, 1, ymax)

  if abs(y₂ - y₁) < abs(x₂ - x₁)
    if x₁ > x₂
      _bresenhamlow!(mask, val, x₂, y₂, x₁, y₁)
    else
      _bresenhamlow!(mask, val, x₁, y₁, x₂, y₂)
    end
  else
    if y₁ > y₂
      _bresenhamhigh!(mask, val, x₂, y₂, x₁, y₁)
    else
      _bresenhamhigh!(mask, val, x₁, y₁, x₂, y₂)
    end
  end
end

function _bresenhamlow!(mask, val, x₁, y₁, x₂, y₂)
  dx = x₂ - x₁
  dy = y₂ - y₁
  yi = 1
  if dy < 0
    yi = -1
    dy = -dy
  end

  D = 2dy - dx
  y = y₁

  for x in x₁:x₂
    mask[x, y] = val

    if D > 0
      y = y + yi
      D = D + 2dy - 2dx
    else
      D = D + 2dy
    end
  end
end

function _bresenhamhigh!(mask, val, x₁, y₁, x₂, y₂)
  dx = x₂ - x₁
  dy = y₂ - y₁
  xi = 1
  if dx < 0
    xi = -1
    dx = -dx
  end

  D = 2dx - dy
  x = x₁

  for y in y₁:y₂
    mask[x, y] = val

    if D > 0
      x = x + xi
      D = D + 2dx - 2dy
    else
      D = D + 2dx
    end
  end
end
