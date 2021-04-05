# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
"""
    mindistance(metric::PreMetric, g::Union{Geometry,Point}, p::Point)

Returns the minimum distance between the the point `p` and the closest point in geometry `g` as
measured by the `metric`.
"""
function mindistance end

# flip arguments to always have geometry be the first argument.
mindistance(metric::PreMetric, p::Point, g::Geometry) = mindistance(metric, g, p)

function mindistance(metric::Union{Euclidean,SqEuclidean}, l::Line, p::Point)
    a, b = l(0), l(1)
    u = p - a
    v = b - a
    α = (u ⋅ v) / (v ⋅ v)
    metric(u, α * v)
end

mindistance(metric::PreMetric, p1::Point, p2::Point) =
    evaluate(metric, coordinates(p1), coordinates(p2))

@deprecate evaluate(d::PreMetric, g::Union{Geometry,Point}, p::Point) mindistance(d, g, p)
@deprecate evaluate(d::PreMetric, p::Point, g::Geometry) mindistance(d, g, p)
