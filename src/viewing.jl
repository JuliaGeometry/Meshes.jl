# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    view(domain, geometry)

Return a view of the `domain` containing all elements that
are inside the `geometry`.
"""
@traitfn function Base.view(domain::D, geometry::Geometry) where {D; !IsGrid{D}}
  view(domain, viewindices(domain, geometry))
end

# in the case of grid + box, we can preserve the type of the grid
@traitfn function Base.view(domain::D, box::Box) where {D; IsGrid{D}}
  domain[viewindices(domain, box)]
end

"""
    viewindices(domain, geometry)

Return the indices of the `domain` that are inside the `geometry`.
"""
@traitfn function viewindices(domain::D, geometry::Geometry) where {D; !IsGrid{D}}
  pred(i) = _isinside(domain[i], geometry)
  filter(pred, 1:nelements(domain))
end

_isinside(p::Point, geometry) = p ∈ geometry
_isinside(g::Geometry, geometry) = g ⊆ geometry

@traitfn function viewindices(domain::D, box::Box) where {D; IsGrid{D}}
  # grid properties
  or = coordinates(minimum(domain))
  sp = spacing(domain)
  sz = size(domain)

  # intersection of boxes
  □ = boundingbox(domain) ∩ box
  lo, up = coordinates.(extrema(□))

  # Cartesian indices of new corners
  ilo = @. max(ceil(Int,  (lo - or) / sp) + 1,  1)
  iup = @. min(floor(Int, (up - or) / sp)    , sz)

  CartesianIndex(Tuple(ilo)):CartesianIndex(Tuple(iup))
end
