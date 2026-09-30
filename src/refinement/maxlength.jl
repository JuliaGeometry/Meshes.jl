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

# regular grids with orthogonal sides can be refined in one shot
function refine(grid::OrthoRegularGrid, method::MaxLengthRefinement)
  esize = sides(boundingbox(grid)) ./ size(grid)
  factors = ceil.(Int, esize ./ method.length)
  refine(grid, RegularRefinement(factors))
end

# other grids are iteratively refined with regular refinement
function refine(grid::Grid, method::MaxLengthRefinement)
  while _iscoarse(grid, method.length)
    grid = refine(grid, RegularRefinement(2))
  end
  grid
end

# general meshes are refined with adaptive edge refinement
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

_iscoarse(grid::Grid, len) = any(g -> _maxside(g) > len, _elements(grid))

_elements(grid::Grid) = (grid[begin], grid[(begin + end) ÷ 2], grid[end])

_maxside(g) = maximum(_sides(g))

_sides(seg::Segment) = (length(seg),)

function _sides(quad::Quadrangle)
  A, B, C, _ = vertices(quad)
  AB = Segment(A, B)
  BC = Segment(B, C)
  length(AB), length(BC)
end

function _sides(hexa::Hexahedron)
  A, B, C, _, E, _, _, _ = vertices(hexa)
  AB = Segment(A, B)
  BC = Segment(B, C)
  AE = Segment(A, E)
  length(AB), length(BC), length(AE)
end
