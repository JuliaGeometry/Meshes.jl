"""
    boundingbox(geometry)

Axis-aligned bounding box of the `geometry`.
"""
boundingbox(geom) = boundingbox(coordinates(geom))

# fallback implementation treats geometry as a set of points
function boundingbox(geometry::AbstractArray{<:Point{N,T}}) where {N,T}
    vmin = vfill(Vec{N,T}, typemax(T))
    vmax = vfill(Vec{N,T}, typemin(T))
    for p in geometry
        vmin, vmax = minmax(coordinates(p), vmin, vmax)
    end
    Rectangle(Point(vmin), vmax - vmin)
end

# --------------
# SPECIAL CASES
# --------------

boundingbox(r::Rectangle) = r

function boundingbox(s::Sphere)
    mini, maxi = extrema(s)
    Rectangle(Point(mini), maxi .- mini)
end
