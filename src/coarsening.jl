# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CoarseningMethod

A method for coarsening meshes.
"""
abstract type CoarseningMethod end

"""
    coarsen(mesh, method)

Coarsen `mesh` with coarsening `method`.
"""
function coarsen end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("coarsening/regular.jl")
