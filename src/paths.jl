# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Path

A path is a sequence of indices that can be
used to traverse a given [`Domain`](@ref).
"""
abstract type Path end

"""
    traverse(domain, path)

Traverse `domain` with `path`.
"""
function traverse end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("paths/linear.jl")
include("paths/random.jl")
include("paths/source.jl")
include("paths/shifted.jl")
include("paths/multigrid.jl")
