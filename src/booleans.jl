# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    g₁ ∩ g₂

Return the intersection of `g₁` and `g₂`.
"""
Base.intersect(g₁::GeometryOrDomain, g₂::GeometryOrDomain) = get(intersection(g₁, g₂))

"""
    g₁ ∪ g₂

Return the union of `g₁` and `g₂`.
"""
Base.union(g₁::GeometryOrDomain, g₂::GeometryOrDomain) = getunion(g₁, g₂)

"""
    setdiff(g₁, g₂)

Return the set difference of `g₁` and `g₂`.
"""
Base.setdiff(g₁::GeometryOrDomain, g₂::GeometryOrDomain) = getdifference(g₁, g₂)

"""
    !g

Return the complement of `g` with respect to its bounding box.
"""
!(g::Geometry) = getcomplement(g)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("booleans/intersections.jl")
include("booleans/unions.jl")
include("booleans/setdiffs.jl")
include("booleans/complements.jl")
