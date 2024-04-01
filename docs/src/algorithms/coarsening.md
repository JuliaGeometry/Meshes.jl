# Coarsening

```@example coarsening
using Meshes # hide
import CairoMakie as Mke # hide
```

```@docs
coarsen
CoarseningMethod
```

## RegularCoarsening

```@docs
RegularCoarsening
```

```@example coarsening
grid = CartesianGrid(100, 100)

# refine three times
cor1 = coarsen(grid, RegularCoarsening(2, 2))
cor2 = coarsen(cor1, RegularCoarsening(3, 2))
cor3 = coarsen(cor2, RegularCoarsening(2, 3))

fig = Mke.Figure(size = (800, 800))
viz(fig[1,1], grid, showfacets = true)
viz(fig[1,2], cor1, showfacets = true)
viz(fig[2,1], cor2, showfacets = true)
viz(fig[2,2], cor3, showfacets = true)
fig
```
