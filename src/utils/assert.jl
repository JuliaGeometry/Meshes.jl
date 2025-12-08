# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    assertion(cond, msg)

Throws an `AssertionError(msg)` if `cond` is `false`.
"""
assertion(cond, msg) = cond || throw(AssertionError(msg))
