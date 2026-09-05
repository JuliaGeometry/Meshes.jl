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

function (::GeodesicDistance)(p₁::Point{🌐}, p₂::Point{🌐})
  # convert coordinates to same LatLon CRS
  q₁, q₂ = promote(p₁, p₂)
  c₁ = convert(manifoldcrs(q₁), coords(q₁))
  c₂ = convert(manifoldcrs(q₂), coords(q₂))

  # the manifold only tells us that the points lie on a sphere,
  # the ellipsoid itself comes from the datum of the coordinates
  🌎 = ellipsoid(datum(c₁))

  # the series of Karney need double precision to reach round-off
  T = numtype(lentype(c₁))
  S = promote_type(T, Float64)
  lat₁, lon₁ = S(ustrip(c₁.lat)), S(ustrip(c₁.lon))
  lat₂, lon₂ = S(ustrip(c₂.lat)), S(ustrip(c₂.lon))

  s₁₂, _, _ = _geodesicinverse(🌎, lat₁, lon₁, lat₂, lon₂)

  T(s₁₂) * unit(majoraxis(🌎))
end
