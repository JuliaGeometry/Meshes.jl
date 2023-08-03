# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# -------------------
# VIEWS WITH INDICES
# -------------------

Base.view(domain::Domain, inds) = DomainView(domain, inds)
Base.view(data::Data, inds) = DataView(data, inds)

# specialize view to avoid infinite loops
Base.view(v::DomainView, inds::AbstractVector{Int}) = DomainView(getfield(v, :domain), getfield(v, :inds)[inds])
Base.view(v::DataView, inds::AbstractVector{Int}) = DataView(getfield(v, :data), getfield(v, :inds)[inds])

# specialize view for grids and Cartesian indices
Base.view(grid::Grid, inds::CartesianIndices) = getindex(grid, inds)

# ---------------------
# UNVIEWS WITH INDICES
# ---------------------

"""
    unview(object)

Return the underlying domain/data of the `object` and
the indices of the view. If the `object` is not a view,
then return the `object` with all its indices as a fallback.
"""
unview(object) = object, 1:nitems(object)
unview(v::DomainView) = getfield(v, :domain), getfield(v, :inds)
unview(v::DataView) = getfield(v, :data), getfield(v, :inds)

# ----------------------
# VIEWS WITH GEOMETRIES
# ----------------------

"""
    view(domain, geometry)

Return a view of the `domain` containing all elements that
are inside the `geometry`.
"""
Base.view(domain::Domain, geometry::Geometry) = view(domain, indices(domain, geometry))

function Base.view(data::Data, geometry::Geometry)
  D = typeof(data)
  dom = domain(data)
  tab = values(data)

  # retrieve subdomain
  inds = indices(dom, geometry)
  subdom = view(dom, inds)

  # retrieve subtable
  tinds = _linear(dom, inds)
  subtab = Tables.subset(tab, tinds)

  # data table for elements
  vals = Dict(paramdim(dom) => subtab)

  constructor(D)(subdom, vals)
end

# convert from Cartesian to linear indices if needed
_linear(domain::Domain, inds) = inds
_linear(grid::Grid, inds) = vec(LinearIndices(size(grid))[inds])

"""
    indices(domain, geometry)

Return the indices of the `domain` that are inside the `geometry`.
"""
indices(domain::Domain, geometry::Geometry) = filter(i -> domain[i] ⊆ geometry, 1:nelements(domain))

function indices(grid::Grid{2}, poly::Polygon{2})
  dims = size(grid)
  mask = falses(dims)
  for triangle in simplexify(poly)
    _fill!(mask, grid, triangle)
  end

  # convert to linear indices
  LinearIndices(dims)[mask]
end

indices(domain::Domain, multi::Multi) = mapreduce(geom -> indices(domain, geom), vcat, collect(multi)) |> unique

function indices(grid::CartesianGrid, box::Box)
  # grid properties
  or = minimum(grid)
  sp = spacing(grid)
  sz = size(grid)

  # intersection of boxes
  □ = boundingbox(grid) ∩ box
  lo, up = extrema(□)

  # Cartesian indices of new corners
  ilo = max.(ceil.(Int, (lo - or) ./ sp) .+ 1, 1)
  iup = min.(floor.(Int, (up - or) ./ sp), sz)

  CartesianIndex(Tuple(ilo)):CartesianIndex(Tuple(iup))
end

# ----------
# UTILITIES
# ----------

"""
    slice(object, xmin:xmax, ymin:ymax, ...)

Slice the `object` using real coordinate ranges `xmin:xmax`, `ymin:ymax`, ...

### Notes

This function is equivalent to `view(object, Box(first.(ranges), last.(ranges))`.

In Julia the range `0.5:10.0` is materialized as `[0.5, ..., 9.5]` so it won't
necessarily include the right value. This behavior is different than the more
intuitive behavior of `view(object, Box((0.5,0.5), (10.0,10.0))`.
"""
slice(object, ranges...) = view(object, Box(first.(ranges), last.(ranges)))

function _fill!(mask, grid, triangle)
  v = vertices(triangle)
  
  # fill edges of triangle
  _bresenham!(mask, grid, v[1], v[2])
  _bresenham!(mask, grid, v[2], v[3])
  _bresenham!(mask, grid, v[3], v[1])

  # fill interior of triangle
  j1 = findfirst(mask).I[2]
  j2 = findlast(mask).I[2]
  for j in j1:j2
    i1 = findfirst(mask[:, j])
    i2 = findlast(mask[:, j])
    mask[i1:i2, j] .= true
  end
end

# Bresenham's line algorithm: https://en.wikipedia.org/wiki/Bresenham's_line_algorithm
function _bresenham!(mask, grid, p₁, p₂)
  o = minimum(grid)
  s = spacing(grid)

  # integer coordinates
  x₁, y₁ = ceil.(Int, (p₁ - o) ./ s)
  x₂, y₂ = ceil.(Int, (p₂ - o) ./ s)

  if abs(y₂ - y₁) < abs(x₂ - x₁)
    if x₁ > x₂
      _bresenhamlow!(mask, x₂, y₂, x₁, y₁)
    else
      _bresenhamlow!(mask, x₁, y₁, x₂, y₂)
    end
  else
    if y₁ > y₂
      _bresenhamhigh!(mask, x₂, y₂, x₁, y₁)
    else
      _bresenhamhigh!(mask, x₁, y₁, x₂, y₂)
    end
  end
end

function _bresenhamlow!(mask, x₁, y₁, x₂, y₂)
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
    mask[x, y] = true

    if D > 0
      y = y + yi
      D = D + 2dy - 2dx
    else
      D = D + 2dy
    end
  end
end

function _bresenhamhigh!(mask, x₁, y₁, x₂, y₂)
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
    mask[x, y] = true

    if D > 0
      x = x + xi
      D = D + 2dx - 2dy
    else
      D = D + 2dx
    end
  end
end
