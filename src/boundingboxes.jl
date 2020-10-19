# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    boundingbox(geometry)

Axis-aligned bounding box of the `geometry`.
"""
boundingbox(geometry) = boundingbox(coordinates(geometry))

# fallback implementation treats geometry as a set of points
function boundingbox(geometry::AbstractArray{<:Point{N,T}}) where {N,T}
    vmin = vfill(Vec{N,T}, typemax(T))
    vmax = vfill(Vec{N,T}, typemin(T))
    for p in geometry
        vmin, vmax = minmax(coordinates(p), vmin, vmax)
    end
    Box(Point(vmin), vmax - vmin)
end

# --------------
# SPECIAL CASES
# --------------

boundingbox(b::Box) = b

function boundingbox(s::Sphere)
    mini, maxi = extrema(s)
    Box(Point(mini), maxi .- mini)
end
