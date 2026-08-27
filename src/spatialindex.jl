# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SpatialIndex

Abstract type for data structures that accelerate spatial queries.
"""
abstract type SpatialIndex end

"""
    candidates(query, index)

Return the indices of elements in `index` whose bounding volumes intersect
the bounding box of `query`.

The returned indices are broad-phase candidates. The corresponding elements
are not guaranteed to intersect `query` exactly.
"""
function candidates end

"""
    candidates!(inds, query, index)

Empty `inds`, store the broad-phase candidate indices for `query`, and return
`inds`.
"""
function candidates! end

include("spatialindex/bvh.jl")
