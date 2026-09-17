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

function refine(grid::Grid, method::MaxLengthRefinement)
  while _iscoarse(grid, method.length)
    grid = refine(grid)
  end
  grid
end

function refine(mesh::Mesh, method::MaxLengthRefinement)
  if paramdim(mesh) == 2
    islong(s) = measure(s) > method.length
    while any(islong, segments(mesh))
      mesh = refine(mesh, EdgeRefinement(islong))
    end
  else
    while _iscoarse(mesh, method.length)
      mesh = refine(mesh)
    end
  end
  mesh
end

#------------------
# HELPER FUNCTIONS
#------------------

function _refinesides(grid, len)
  esize = sides(boundingbox(grid)) ./ size(grid)
  factors = ceil.(Int, esize ./ len)
  refine(grid, RegularRefinement(factors))
end

_iscoarse(mesh::Mesh, len) = any(g -> _maxside(g) > len, _elements(mesh))

_elements(mesh::Mesh) = mesh

_elements(grid::Grid) = [grid[begin], grid[(begin + end) ÷ 2], grid[end]]

_maxside(g) = maximum(_sides(g))

_sides(seg::Segment) = (measure(seg),)

function _sides(tri::Triangle)
  A, B, C = vertices(tri)
  AB = Segment(A, B)
  AC = Segment(A, C)
  measure(AB), measure(AC)
end

function _sides(quad::Quadrangle)
  A, B, C, _ = vertices(quad)
  AB = Segment(A, B)
  BC = Segment(B, C)
  measure(AB), measure(BC)
end

function _sides(tetra::Tetrahedron)
  A, B, C, D = vertices(tetra)
  AB = Segment(A, B)
  AC = Segment(A, C)
  AD = Segment(A, D)
  measure(AB), measure(AC), measure(AD)
end

function _sides(hexa::Hexahedron)
  A, B, C, _, E, _, _, _ = vertices(hexa)
  AB = Segment(A, B)
  BC = Segment(B, C)
  AE = Segment(A, E)
  measure(AB), measure(BC), measure(AE)
end
