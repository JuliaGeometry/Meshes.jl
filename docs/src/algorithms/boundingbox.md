# Bounding box

```@example boundingbox
using Meshes # hide
import CairoMakie as Mke # hide
```

```@docs
boundingbox
```

```@example boundingbox
pset = PointSet(rand(Point{2}, 100))
bbox = boundingbox(pset)

fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], bbox)
viz!(fig[1,1], pset, color = :black)
fig
```

```@example boundingbox
box  = Box((-1, -1), (0, 0))
ball = Ball((0, 0), (1))
gset = GeometrySet([box, ball])
bbox = boundingbox(gset)

fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], bbox)
viz!(fig[1,1], boundary(box), color = :gray)
viz!(fig[1,1], boundary(ball), color = :gray)
fig
```
