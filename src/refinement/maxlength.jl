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

refine(grid::OrthoRegularGrid, method::MaxLengthRefinement) = _refineregular(grid, method)

refine(grid::RectilinearGrid, method::MaxLengthRefinement) = _refineregular(grid, method)

refine(grid::OrthoStructuredGrid, method::MaxLengthRefinement) = _refineregular(grid, method)

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

function _refineregular(grid, method)
  esize = sides(boundingbox(grid)) ./ size(grid)
  factors = ceil.(Int, esize ./ method.length)
  refine(grid, RegularRefinement(factors))
end
