# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    boundingbox(geometry)

Axis-aligned bounding box of the `geometry`.
"""
boundingbox(geometry) = boundingbox(coordinates(geometry))

# fallback implementation treats geometry as a set of points
function boundingbox(points::AbstractArray{<:Point{Dim,T}}) where {Dim,T}
    xmin = MVector(ntuple(i->typemax(T), Dim))
    xmax = MVector(ntuple(i->typemin(T), Dim))
    for p in points
        x = coordinates(p)
        @. xmin = min(x, xmin)
        @. xmax = max(x, xmax)
    end
    Box(Point(xmin), Point(xmax))
end

# -----------
# PRIMITIVES
# -----------

boundingbox(b::Box) = b

function boundingbox(s::Sphere{Dim,T}) where {Dim,T}
    c = center(s)
    r = radius(s)
    r⃗ = vfill(Vec{Dim,T}, r)
    Box(c - r⃗, c + r⃗)
end
