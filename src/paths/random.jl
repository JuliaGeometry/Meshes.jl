# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RandomPath

Traverse an object with `N` elements in a random
permutation of `1:N`.
"""
struct RandomPath <: Path end

traverse(object, ::RandomPath) = randperm(nelements(object))
