# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RandomPath()

Traverse a domain with `N` elements in a random
permutation of `1:N`.
"""
struct RandomPath <: Path end

traverse(domain, ::RandomPath) = randperm(nelements(domain))
