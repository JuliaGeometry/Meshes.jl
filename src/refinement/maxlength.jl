# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MaxLengthRefinement(length)

Refine mesh into elements with all boundary segments smaller than
or equal to a maximum `length` in length units (default to meters).
"""
struct MaxLengthRefinement{ℒ<:Len} <: RefinementMethod
  length::ℒ
  MaxLengthRefinement(length::ℒ) where {ℒ<:Len} = new{float(ℒ)}(length)
end

MaxLengthRefinement(length) = MaxLengthRefinement(aslen(length))

refine(grid::OrthoRegularGrid, method::MaxLengthRefinement) = _refinesides(grid, method)

refine(grid::RectilinearGrid, method::MaxLengthRefinement) = _refinesides(grid, method)

refine(grid::OrthoStructuredGrid, method::MaxLengthRefinement) = _refinesides(grid, method)

function refine(grid::RegularGrid, method::MaxLengthRefinement)
  iscoarse(e) = perimeter(e) > method.length * nvertices(e)
  while iscoarse(first(grid)) || iscoarse(last(grid))
    grid = refine(grid, RegularRefinement(2))
  end
  grid
end

function refine(mesh::Mesh, method::MaxLengthRefinement)
  iscoarse(e) = perimeter(e) > method.length * nvertices(e)
  while any(iscoarse, mesh)
    mesh = refine(mesh, TriSubdivision())
  end
  mesh
end

#------------------
# HELPER FUNCTIONS
#------------------

function _refinesides(grid, method)
  esize = sides(boundingbox(grid)) ./ size(grid)
  factors = ceil.(Int, esize ./ method.length)
  refine(grid, RegularRefinement(factors))
end
