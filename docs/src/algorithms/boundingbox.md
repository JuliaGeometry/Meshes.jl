# Bounding box

```@example boundingbox
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example boundingbox
using Meshes, MeshViz # hide
import WGLMakie as Mke # hide
```

```@docs
boundingbox
```

```@example boundingbox
pset = PointSet(rand(Point2, 100))
bbox = boundingbox(pset)

fig = Mke.Figure(resolution = (800, 400))
viz(fig[1,1], bbox)
viz!(fig[1,1], pset, color = :black)
fig
```
