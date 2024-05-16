# Hulls

```@example hull
using Meshes # hide
import CairoMakie as Mke # hide
```

```@docs
hull
convexhull
HullMethod
GrahamScan
JarvisMarch
```

```@example hull
pset = PointSet([rand(Point{2}) for _ in 1:100])
chul = convexhull(pset)

fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], chul)
viz!(fig[1,1], pset, color = :black)
fig
```

```@example hull
box  = Box((-1, -1), (0, 0))
ball = Ball((0, 0), (1))
gset = GeometrySet([box, ball])
chul = convexhull(gset)

fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], chul)
viz!(fig[1,1], boundary(box), color = :gray)
viz!(fig[1,1], boundary(ball), color = :gray)
fig
```
