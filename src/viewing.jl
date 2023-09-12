# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# -------------------
# VIEWS WITH INDICES
# -------------------

Base.view(domain::Domain, inds) = DomainView(domain, inds)

# ---------------------
# UNVIEWS WITH INDICES
# ---------------------

"""
    unview(object)

Return the underlying domain/data of the `object` and
the indices of the view. If the `object` is not a view,
then return the `object` with all its indices as a fallback.
"""
function unview end

unview(d::Domain) = d, 1:nelements(d)
unview(v::DomainView) = getfield(v, :domain), getfield(v, :inds)

# ----------------------
# VIEWS WITH GEOMETRIES
# ----------------------

"""
    view(domain, geometry)

Return a view of the `domain` containing all elements that
intersect with the `geometry`.
"""
Base.view(domain::Domain, geometry::Geometry) = view(domain, indices(domain, geometry))

"""
    indices(domain, geometry)

Return the indices of the elements of the `domain`
that intersect with the `geometry`.
"""
indices(domain::Domain, geometry::Geometry) = findall(intersects(geometry), domain)

function indices(grid::Grid{2}, poly::Polygon{2})
  dims = size(grid)
  mask = zeros(Int, dims)
  cpoly = clip(poly, boundingbox(grid), SutherlandHodgman())

  for (i, triangle) in enumerate(simplexify(cpoly))
    _fill!(mask, grid, i, triangle)
  end

  # convert to linear indices
  LinearIndices(dims)[mask .> 0]
end

function indices(grid::Grid{2}, chain::Chain{2})
  dims = size(grid)
  mask = falses(dims)

  for segment in segments(chain)
    p₁, p₂ = vertices(segment)
    _bresenham!(mask, grid, true, p₁, p₂)
  end

  # convert to linear indices
  LinearIndices(dims)[mask]
end

indices(domain::Domain, multi::Multi) = mapreduce(geom -> indices(domain, geom), vcat, parent(multi)) |> unique

function indices(grid::CartesianGrid, box::Box)
  # grid properties
  or = minimum(grid)
  sp = spacing(grid)
  sz = size(grid)

  # intersection of boxes
  lo, up = extrema(boundingbox(grid) ∩ box)

  # Cartesian indices of new corners
  ilo = max.(ceil.(Int, (lo - or) ./ sp), 1)
  iup = min.(floor.(Int, (up - or) ./ sp) .+ 1, sz)

  # Cartesian range from corner to corner
  range = CartesianIndex(Tuple(ilo)):CartesianIndex(Tuple(iup))

  # convert to linear indices
  LinearIndices(sz)[range] |> vec
end

# utils
function _fill!(mask, grid, fval, triangle)
  v = vertices(triangle)

  # fill edges of triangle
  _bresenham!(mask, grid, fval, v[1], v[2])
  _bresenham!(mask, grid, fval, v[2], v[3])
  _bresenham!(mask, grid, fval, v[3], v[1])

  # fill interior of triangle
  j₁ = findfirst(==(fval), mask).I[2]
  j₂ = findlast(==(fval), mask).I[2]
  for j in j₁:j₂
    i₁ = findfirst(==(fval), @view(mask[:, j]))
    i₂ = findlast(==(fval), @view(mask[:, j]))
    mask[i₁:i₂, j] .= fval
  end
end

# Bresenham's line algorithm: https://en.wikipedia.org/wiki/Bresenham's_line_algorithm
function _bresenham!(mask, grid, fval, p₁, p₂)
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
      _bresenhamlow!(mask, fval, x₂, y₂, x₁, y₁)
    else
      _bresenhamlow!(mask, fval, x₁, y₁, x₂, y₂)
    end
  else
    if y₁ > y₂
      _bresenhamhigh!(mask, fval, x₂, y₂, x₁, y₁)
    else
      _bresenhamhigh!(mask, fval, x₁, y₁, x₂, y₂)
    end
  end
end

function _bresenhamlow!(mask, fval, x₁, y₁, x₂, y₂)
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
    mask[x, y] = fval

    if D > 0
      y = y + yi
      D = D + 2dy - 2dx
    else
      D = D + 2dy
    end
  end
end

function _bresenhamhigh!(mask, fval, x₁, y₁, x₂, y₂)
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
    mask[x, y] = fval

    if D > 0
      x = x + xi
      D = D + 2dx - 2dy
    else
      D = D + 2dx
    end
  end
end
