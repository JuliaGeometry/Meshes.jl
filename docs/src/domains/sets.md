# Sets

```@example sets
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example sets
using Meshes, MeshViz # hide
import WGLMakie as Mke # hide
```

Geometry sets represent a collection of geometries without
any connectivity information (a.k.a., "soup of geometries").

```@docs
GeometrySet
```

```@example sets
GeometrySet(rand(Ball{3,Float64}, 3)) |> viz
```

```@docs
PointSet
```

```@example sets
PointSet(rand(Point2, 100)) |> viz
```