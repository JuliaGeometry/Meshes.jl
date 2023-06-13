# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Path

A path is a sequence of indices that can be
used to traverse a given [`Domain`](@ref).

## References

* Nussbaumer et al. 2017. [Which Path to Choose in Sequential Gaussian Simulation]
  (https://link.springer.com/article/10.1007/s11004-017-9699-5)
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

include("traversing/linear.jl")
include("traversing/random.jl")
include("traversing/source.jl")
include("traversing/shifted.jl")
include("traversing/multigrid.jl")
