# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LinearPath()

Traverse a domain with `n` elements in order `1:n`.
"""
struct LinearPath <: Path end

traverse(domain, ::LinearPath) = 1:nelements(domain)
