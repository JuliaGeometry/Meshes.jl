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

refine(grid::OrthoRegularGrid, method::MaxLengthRefinement) = _refinesides(grid, method.length)

refine(grid::RectilinearGrid, method::MaxLengthRefinement) = _refinesides(grid, method.length)

refine(grid::OrthoStructuredGrid, method::MaxLengthRefinement) = _refinesides(grid, method.length)

function refine(mesh::Mesh, method::MaxLengthRefinement)
  while true
    rmesh = refine(mesh, EdgeRefinement(e -> length(e) > method.length))
    nelements(rmesh) == nelements(mesh) && return mesh
    mesh = rmesh
  end
end

#------------------
# HELPER FUNCTIONS
#------------------

function _refinesides(grid, len)
  esize = sides(boundingbox(grid)) ./ size(grid)
  factors = ceil.(Int, esize ./ len)
  refine(grid, RegularRefinement(factors))
end
