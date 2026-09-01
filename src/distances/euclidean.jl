# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EuclideanDistance()

Length of the straight line connecting two points in the space
where they are embedded. For points on the ellipsoid (`🌐`) this
is the chord through the interior of the ellipsoid, which is
shorter than any path along the surface.

See also [`GeodesicDistance`](@ref).

## Examples

```julia
d = EuclideanDistance()

d(Point(0, 0), Point(1, 1))

d(Point(LatLon(0, 0)), Point(LatLon(0, 1)))
```
"""
struct EuclideanDistance <: GeometricDistance end

(::EuclideanDistance)(p₁::Point{M}, p₂::Point{M}) where {M} = norm(p₂ - p₁)
