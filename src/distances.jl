# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    closest_point(m::PreMetric, g::Union{Geometry,Point}, p::Point)

Returns the closest point in `g` to the point `p` as measured by the metric `m`.
"""
function closest_point end

closest_point(::PreMetric, p::Point, ::Point) = p

function closest_point(::Union{Euclidean,SqEuclidean}, l::Line, p::Point)
    a, b = l(0), l(1)
    u = p - a
    v = b - a
    α = (u ⋅ v) / (v ⋅ v)
    l(α)
end

"""
    mindistance(metric::PreMetric, g::Union{Geometry,Point}, p::Point)

Returns the minimum distance between the the point `p` and the closest point in geometry `g` as
measured by the `metric`.
"""
function mindistance end

mindistance(metric::PreMetric, p1::Point, p2::Point) =
    evaluate(metric, coordinates(p1), coordinates(p2))
# flip arguments to always have geometry be the first argument.
mindistance(metric::PreMetric, p::Point, g::Geometry) = mindistance(metric, g, p)
mindistance(metric::PreMetric, g::Geometry, p::Point) =
    mindistance(metric, closest_point(metric, g, p), p)

@deprecate evaluate(d::PreMetric, g::Union{Geometry,Point}, p::Point) mindistance(d, g, p)
@deprecate evaluate(d::PreMetric, p::Point, g::Geometry) mindistance(d, g, p)
