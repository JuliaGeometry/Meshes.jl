# Refinement

```@docs
refine
RefinementMethod
```

## TriRefinement

```@docs
TriRefinement
```

```@example
using Meshes, MeshViz
import CairoMakie as Mke

grid = CartesianGrid(10, 10)

# refine three times
ref1 = refine(grid, TriRefinement())
ref2 = refine(ref1, TriRefinement())
ref3 = refine(ref2, TriRefinement())

fig = Mke.Figure(resolution = (800, 800))
viz(fig[1,1], grid, showfacets = true, axis = (title = "original",))
viz(fig[1,2], ref1, showfacets = true, axis = (title = "refine 1",))
viz(fig[2,1], ref2, showfacets = true, axis = (title = "refine 2",))
viz(fig[2,2], ref3, showfacets = true, axis = (title = "refine 3",))
fig
```

## QuadRefinement

```@docs
QuadRefinement
```

```@example
using Meshes, MeshViz
import CairoMakie as Mke

grid = CartesianGrid(10, 10)

# refine three times
ref1 = refine(grid, QuadRefinement())
ref2 = refine(ref1, QuadRefinement())
ref3 = refine(ref2, QuadRefinement())

fig = Mke.Figure(resolution = (800, 800))
viz(fig[1,1], grid, showfacets = true, axis = (title = "original",))
viz(fig[1,2], ref1, showfacets = true, axis = (title = "refine 1",))
viz(fig[2,1], ref2, showfacets = true, axis = (title = "refine 2",))
viz(fig[2,2], ref3, showfacets = true, axis = (title = "refine 3",))
fig
```

## Catmull-Clark

```@docs
CatmullClark
```

```@example
using Meshes, MeshViz
import CairoMakie as Mke

# define a cube in R^3
points = Point3[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
connec = connect.([(1,4,3,2),(5,6,7,8),(1,2,6,5),(3,4,8,7),(1,5,8,4),(2,3,7,6)])
mesh   = SimpleMesh(points, connec)

# refine three times
ref1 = refine(mesh, CatmullClark())
ref2 = refine(ref1, CatmullClark())
ref3 = refine(ref2, CatmullClark())

fig = Mke.Figure(resolution = (800, 800))
viz(fig[1,1], mesh, showfacets = true, axis = (title = "original",))
viz(fig[1,2], ref1, showfacets = true, axis = (title = "refine 1",))
viz(fig[2,1], ref2, showfacets = true, axis = (title = "refine 2",))
viz(fig[2,2], ref3, showfacets = true, axis = (title = "refine 3",))
fig
```
