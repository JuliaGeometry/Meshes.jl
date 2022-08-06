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
p1 = Point(0.,0.)
p2 = Point(1.,0.)
p3 = Point(0.,1.)
p4 = Point(1.,1.)

viz([p1, p2, p3, p4], color = [1,2,3,4])
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
