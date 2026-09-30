# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    union(polygon₁, polygon₂)

Return the union of `polygon₁` and `polygon₂`, or `nothing`
if the polygons are disjoint.
"""
Base.union(poly₁::Polygon, poly₂::Polygon) = _ghclip(poly₁, poly₂, :union)

"""
    setdiff(polygon₁, polygon₂)

Return the part of `polygon₁` that is not covered by `polygon₂`,
or `nothing` if `polygon₂` covers `polygon₁`.
"""
Base.setdiff(poly₁::Polygon, poly₂::Polygon) = _ghclip(poly₁, poly₂, :difference)
