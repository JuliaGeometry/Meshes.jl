# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DiscretizationMethod

A method for tesselating point sets into meshes.
"""
abstract type TesselationMethod end

"""
    tesselate(pointset, [method])

Tesselate `pointset` with tesselation `method`.

If the `method` is ommitted, a default algorithm is used.
"""
function tesselate end

tesselate(points::AbstractVector{<:Point{2}}, method::TesselationMethod) = tesselate(PointSet(points), method)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("tesselation/delaunay.jl")
