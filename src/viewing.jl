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
