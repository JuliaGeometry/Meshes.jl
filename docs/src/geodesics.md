# Geodesics

```@example geodesics
using Meshes # hide
using CoordRefSystems # hide
```

A geodesic is the shortest path between two points along a manifold.
On the ellipsoid there are two classical problems, and each has its
own function. Both are solved with the series of Karney (2013), on the
ellipsoid attached to the datum of the coordinates.

```@docs
geodesicfwd
geodesicbwd
geodesictangent
geodesicazimuth
```

The direct problem walks a given length along a given azimuth. Plain
numbers are taken to be degrees and meters:

```@example geodesics
p = Point(LatLon(-33.8688, 151.2093))

geodesicfwd(p, 45, 1000000)
```

The inverse problem recovers the azimuth that connects two points.
Its length is [`GeodesicDistance`](@ref), so the two together describe
the geodesic completely:

```@example geodesics
sydney = Point(LatLon(-33.8688, 151.2093))
london = Point(LatLon(51.5074, -0.1278))

geodesicdist = GeodesicDistance()

ϕ = geodesicbwd(sydney, london)
d = geodesicdist(sydney, london)

ϕ, d
```

Walking that azimuth for that length arrives at the second point:

```@example geodesics
geodesicfwd(sydney, geodesicbwd(sydney, london), geodesicdist(sydney, london))
```

The azimuth at the far end points back along the same geodesic, and is
obtained by swapping the arguments:

```@example geodesics
geodesicbwd(london, sydney)
```

## Tangent vectors

An azimuth at a point can also be expressed as a unit vector tangent to
the ellipsoid there, in geocentric Cartesian coordinates:

```@example geodesics
geodesictangent(sydney, geodesicbwd(sydney, london))
```

The conversion goes both ways, and any component of the vector along the
normal of the ellipsoid is ignored:

```@example geodesics
geodesicazimuth(sydney, geodesictangent(sydney, 45))
```
