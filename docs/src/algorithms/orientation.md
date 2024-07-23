# Orientation

```@example orientation
using Meshes # hide
import CairoMakie as Mke # hide
```

Many geometric processing algorithms for 2D geometries
rely on the concept of [`orientation`](@ref), which is
illustrated below.

```@docs
OrientationType
orientation
```

For polygons without holes, the function returns the
orientation of the boundary, which is a [`Ring`](@ref):

```@example orientation
tri = Triangle((0, 0), (1, 0), (0, 1))

orientation(tri)
```

```@example orientation
tri = Triangle((0, 0), (0, 1), (1, 0))

orientation(tri)
```

For polygons with holes, the function returns a
vector with the orientation of all constituent rings:

```@example orientation
outer = [(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)]
hole1 = [(0.2,0.2),(0.4,0.2),(0.4,0.4),(0.2,0.4)]
hole2 = [(0.6,0.2),(0.8,0.2),(0.8,0.4),(0.6,0.4)]
poly  = PolyArea([outer, hole1, hole2])

orientation(poly)
```
