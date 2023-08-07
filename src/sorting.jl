# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SortingMethod

A method for sorting geometric objects.
"""
abstract type SortingMethod end

"""
    sort(domain, method)

Sort the elements of the `domain` with given sorting `method`.
"""
sort(domain::Domain, method::SortingMethod) =
  view(domain, sortinds(domain, method))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("sorting/direction.jl")
