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
intersect with the `geometry`.
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
  subtab = Tables.subset(tab, inds)

  # data table for elements
  vals = Dict(paramdim(dom) => subtab)

  constructor(D)(subdom, vals)
end

"""
    indices(domain, geometry)

Return the indices of the elements of the `domain`
that intersect with the `geometry`.
"""
indices(domain::Domain, geometry::Geometry) = filter(i -> intersects(domain[i], geometry), 1:nelements(domain))

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
