# Sets

```@example geomsets
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example geomsets
using Meshes, MeshViz # hide
import WGLMakie as Mke # hide
```

Geometry sets represent a collection of geometries without
any connectivity information (a.k.a., "soup of geometries").

```@docs
GeometrySet
PointSet
```

```@example geomsets
GeometrySet([Point(0, 0), Ball((0, 0), 1), Box((0, 0), (1, 1))])
```

```@example geomsets
GeometrySet(rand(Point2, 10))
```