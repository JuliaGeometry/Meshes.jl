# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SortingMethod

A method for sorting objects.
"""
abstract type SortingMethod end

"""
    sort(object, method)

Sort the elements of the `object` with given sorting `method`.
"""
sort(object, method::SortingMethod) = view(object, sortinds(object, method))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("sorting/direction.jl")
