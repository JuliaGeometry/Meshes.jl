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
"""
function refine end

refine(mesh::Mesh) = refine(mesh, QuadRefinement())

refine(grid::Grid) = refine(grid, RegularRefinement(2))

"""
    refinemaxlen(mesh)

Refine `mesh` to ensure that no element exceeds the
maximum length specified by [`maxlen`](@ref).
"""
refinemaxlen(mesh::Mesh) = refine(mesh, MaxLengthRefinement(maxlen()))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("refinement/tri.jl")
include("refinement/quad.jl")
include("refinement/trisub.jl")
include("refinement/catmullclark.jl")
include("refinement/regular.jl")
include("refinement/maxlength.jl")
