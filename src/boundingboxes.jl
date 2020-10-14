"""
    boundingbox(geometry)

Axis-aligned bounding box of the `geometry`.
"""
boundingbox(geom) = boundingbox(coordinates(geom))

# fallback implementation treats geometry as
# a set of points (i.e. coordinates)
function boundingbox(geometry::AbstractArray{<:AbstractPoint{N,T}}) where {N,T}
    vmin = vfill(Vec{N,T}, typemax(T))
    vmax = vfill(Vec{N,T}, typemin(T))
    for p in geometry
        vmin, vmax = minmax(coordinates(p), vmin, vmax)
    end
    Rect{N,T}(vmin, vmax - vmin)
end

# --------------
# SPECIAL CASES
# --------------

boundingbox(a::Rect) = a

function boundingbox(a::Pyramid{T}) where {T}
    w = a.width / 2
    h = a.length
    m = a.middle
    Rect{3,T}(m - Point{3,T}(w, w, 0), m + Point{3,T}(w, w, h))
end

function boundingbox(a::Sphere{T}) where {T}
    mini, maxi = extrema(a)
    Rect{3,T}(mini, maxi .- mini)
end
