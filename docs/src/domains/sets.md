# Sets

```@example sets
using Meshes # hide
using CoordRefSystems # hide
import CairoMakie as Mke # hide
```

Geometry sets represent a collection of geometries without
any connectivity information (a.k.a., "soup of geometries").

```@docs
GeometrySet
```

```@example sets
GeometrySet(rand(Ball, 3)) |> viz
```

```@docs
PointSet
```

```@example sets
PointSet(rand(Point, 100, crs=Cartesian2D)) |> viz
```