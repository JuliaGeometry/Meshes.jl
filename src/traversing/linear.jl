# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LinearPath()

Traverse a domain with `N` elements in order `1:N`.
"""
struct LinearPath <: Path end

traverse(domain, ::LinearPath) = 1:nelements(domain)
