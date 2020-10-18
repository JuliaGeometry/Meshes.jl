"""
    boundingbox(geometry)

Axis-aligned bounding box of the `geometry`.
"""
boundingbox(geom) = boundingbox(coordinates(geom))

# fallback implementation treats geometry as a set of points
function boundingbox(geometry::AbstractArray{<:AbstractPoint{N,T}}) where {N,T}
    vmin = vfill(Vec{N,T}, typemax(T))
    vmax = vfill(Vec{N,T}, typemin(T))
    for p in geometry
        vmin, vmax = minmax(coordinates(p), vmin, vmax)
    end
    Rectangle{N,T}(Point(vmin), vmax - vmin)
end

# --------------
# SPECIAL CASES
# --------------

boundingbox(r::Rectangle) = r

function boundingbox(a::Pyramid{T}) where {T}
    w = a.width / 2
    h = a.length
    m = a.middle
    Rectangle{3,T}(m - Point{3,T}(w, w, 0), m + Point{3,T}(w, w, h))
end

function boundingbox(s::Sphere)
    mini, maxi = extrema(s)
    Rectangle(Point(mini), maxi .- mini)
end
