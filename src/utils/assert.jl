# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    assertion(cond, msg)

Throws an `AssertionError(msg)` if `cond` is `false`.
"""
assertion(cond, msg) = cond || throw(AssertionError(msg))

"""
    checkdim(geom, dim)

Throws an `ArgumentError` if the `embeddim` of the geometry `geom`
is different than the specified dimension `dim`. 
"""
checkdim(geom, dim) =
  embeddim(geom) ≠ dim && throw(ArgumentError("geometry must be embedded in $dim-dimensional space"))
