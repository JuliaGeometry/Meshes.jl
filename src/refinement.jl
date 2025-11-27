# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RefinementMethod

A method for refining meshes.
"""
abstract type RefinementMethod end

"""
    refine(mesh, method)

Refine `mesh` with refinement `method`.
"""
function refine end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("refinement/tri.jl")
include("refinement/quad.jl")
include("refinement/regular.jl")
include("refinement/catmullclark.jl")
include("refinement/trisubdivision.jl")
include("refinement/maxlength.jl")
