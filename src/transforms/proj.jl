# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Proj(CRS)

TODO
"""
struct Proj{CRS} <: CoordinateTransform end

Proj(CRS) = Proj{CRS}()

applycoord(::Proj, v::Vec) = v

applycoord(::Proj{CRS}, p::Point) where {CRS} = Point(convert(CRS, coords(p)))
