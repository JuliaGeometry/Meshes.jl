# Points

```@example points
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example points
using Meshes, MeshViz # hide
import WGLMakie as Mke # hide
```

```@docs
Point
```

```@example points
ps = rand(Point3, 100)

viz(ps)
```

```@docs
embeddim(::Point)
coordtype(::Point)
coordinates(::Point)
-(::Point, ::Point)
+(::Point, ::Vec)
-(::Point, ::Vec)
isapprox(::Point, ::Point)
```
