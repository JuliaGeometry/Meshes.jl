# Distances

```@example distances
using Meshes # hide
using CoordRefSystems # hide
```

Distances are functors that are evaluated with a pair of points. The
definition that is used depends on the manifold of the points, which
is available at compile time.

```@docs
GeometricDistance
EuclideanDistance
GeodesicDistance
```

In Euclidean space the Euclidean distance is the length of the straight
line connecting the points:

```@example distances
d = EuclideanDistance()

d(Point(0, 0), Point(3, 4))
```

On the ellipsoid the same distance is the chord that goes through the
interior of the Earth, whereas the geodesic distance is the length of
the shortest path along the surface:

```@example distances
sydney = Point(LatLon(-33.8688, 151.2093))
london = Point(LatLon(51.5074, -0.1278))

EuclideanDistance()(sydney, london), GeodesicDistance()(sydney, london)
```

The ellipsoid is taken from the datum of the coordinate reference system,
so all parameters are retrieved automatically for the user.
