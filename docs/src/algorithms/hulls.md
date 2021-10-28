# Hulls

```@docs
hull
HullMethod
```

## Graham's scan

```@docs
GrahamScan
```

```@example hull
using Meshes, MeshViz
import CairoMakie

pset = PointSet(rand(Point2, 100))
chul = hull(pset, GrahamScan())

fig = CairoMakie.Figure(resolution = (800, 400))
viz(fig[1,1], chul)
viz!(fig[1,1], pset, color = :black)
fig
```