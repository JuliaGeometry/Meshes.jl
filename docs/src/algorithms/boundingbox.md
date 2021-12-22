# Bounding box

```@docs
boundingbox
```

```@example bbox
using Meshes, MeshViz
import CairoMakie as Mke

pset = PointSet(rand(Point2, 100))
bbox = boundingbox(pset)

fig = Mke.Figure(resolution = (800, 400))
viz(fig[1,1], bbox)
viz!(fig[1,1], pset, color = :black)
fig
```