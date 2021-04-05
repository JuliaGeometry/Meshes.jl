# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _projected_α(l::Union{Line,Segment}, p::Point)
    a, b = l(0), l(1)
    u = p - a
    v = b - a
    (u ⋅ v) / (v ⋅ v)
end

"""
    closest_point(m::PreMetric, g::Union{Geometry,Point}, p::Point)

Returns the closest point in `g` to the point `p` as measured by the metric `m`.

There may be multiple points in `g` that are closest to `p`. Any of such points may be returned by
this function.
"""
function closest_point end

closest_point(::PreMetric, p::Point, ::Point) = p
closest_point(::Union{Euclidean,SqEuclidean}, l::Line, p::Point) = l(_projected_α(l, p))

function closest_point(::Union{Euclidean,SqEuclidean}, s::Segment, p::Point)
    p1, p2 = vertices(s)
    if p1 == p2
        p1
    else
        α = _projected_α(s, p)
        if α <= 0
            p1
        elseif α >= 1
            p2
        else
            s(α)
        end
    end
end

function closest_point(metric::Union{Euclidean,SqEuclidean}, c::Chain, p::Point)
    reduce(segments(c); init = (; distance = Inf, point = first(vertices(c)))) do best, s
        point = closest_point(metric, s, p)
        distance = mindistance(metric, point, p)
        distance < best.distance ? (; distance, point) : best
    end.point
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
