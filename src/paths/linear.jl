# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LinearPath

Traverse an object with `N` elements in order `1:N`.
"""
struct LinearPath <: Path end

traverse(object, ::LinearPath) = 1:nelements(object)
