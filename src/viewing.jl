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

function indices(grid::Grid{2}, polygon::Polygon{2})
  mask = zeros(Int, size(grid))
  linds = LinearIndices(size(grid))

  for (n, ring) in enumerate(rings(polygon))
    for seg in segments(ring)
      s = spacing(grid)
      p1, p2 = vertices(seg)
      x1, y1 = ceil.(Int, coordinates(p1) ./ s)
      x2, y2 = ceil.(Int, coordinates(p2) ./ s)
      cind1 = CartesianIndex(x1, y1)
      cind2 = CartesianIndex(x2, y2)
      inds = bresenham(cind1, cind2)
      mask[inds] .= n
    end

    # fill external ring and unfill internal rings
    f = n > 1 : 0 : 1
    for col in eachcol(mask)
      find = findfirst(==(n), col)
      if !isnothing(find)
        lind = findlast(==(n), col)
        if find ≠ lind # skip single vertice
          for i in find:lind
            # don't fill horizontal lines
            col[i] ≠ n && (col[i] = f)
          end
        end
      end
    end
  end

  linds[mask .> 0]
end

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
