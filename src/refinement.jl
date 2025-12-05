# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RefinementMethod

A method for refining meshes.
"""
abstract type RefinementMethod end

"""
    refine(mesh, [method])

Refine `mesh` with refinement `method`.

If the `method` is omitted, a default is used as a
function of the `mesh`. Grids are refined into finer
grids using regular refinement with a factor of two.
"""
function refine end

refine(mesh::Mesh) = refine(mesh, TriRefinement())

refine(grid::Grid) = refine(grid, RegularRefinement(2))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("refinement/tri.jl")
include("refinement/quad.jl")
include("refinement/trisub.jl")
include("refinement/catmullclark.jl")
include("refinement/regular.jl")
include("refinement/maxlength.jl")
