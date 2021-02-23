# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@traitdef IsGrid{D}
@traitimpl IsGrid{D} <- isgrid(D)

"""
    isgrid(D)

Tells whether or not a domain or data of type `D` behaves as a grid
(or array) of elements. This information can be used to optimize
algorithms (e.g. slices, bounding boxes).
"""
isgrid(D) = false
