# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeodesicDistance()

Length of the shortest path along the manifold of two points.
In Euclidean space (`𝔼`) this is the straight line connecting
them. On the ellipsoid (`🌐`) this is the geodesic of the
ellipsoid attached to the datum of the coordinate reference
system, and is computed with the series of Karney (2013).

See also [`EuclideanDistance`](@ref).

## Examples

```julia
d = GeodesicDistance()

d(Point(LatLon(0, 0)), Point(LatLon(0, 1)))
```

## References

* Karney, C. F. F. 2013. [Algorithms for geodesics](https://doi.org/10.1007/s00190-012-0578-z)
"""
struct GeodesicDistance <: GeometricDistance end

(::GeodesicDistance)(p₁::Point{𝔼{Dim}}, p₂::Point{𝔼{Dim}}) where {Dim} = norm(p₂ - p₁)

(::GeodesicDistance)(p₁::Point{🌐}, p₂::Point{🌐}) = geodesicdistance(coords(p₁), coords(p₂))
