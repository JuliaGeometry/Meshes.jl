# Refinement

```@example refinement
using Meshes # hide
import CairoMakie as Mke # hide
```

```@docs
refine
RefinementMethod
```

## TriRefinement

```@docs
TriRefinement
```

```@example refinement
grid = CartesianGrid(10, 10)

# refine three times
ref1 = refine(grid, TriRefinement())
ref2 = refine(ref1, TriRefinement())
ref3 = refine(ref2, TriRefinement())

fig = Mke.Figure(size = (800, 800))
viz(fig[1,1], grid, showsegments = true)
viz(fig[1,2], ref1, showsegments = true)
viz(fig[2,1], ref2, showsegments = true)
viz(fig[2,2], ref3, showsegments = true)
fig
```

## QuadRefinement

```@docs
QuadRefinement
```

```@example refinement
grid = CartesianGrid(10, 10)

# refine three times
ref1 = refine(grid, QuadRefinement())
ref2 = refine(ref1, QuadRefinement())
ref3 = refine(ref2, QuadRefinement())

fig = Mke.Figure(size = (800, 800))
viz(fig[1,1], grid, showsegments = true)
viz(fig[1,2], ref1, showsegments = true)
viz(fig[2,1], ref2, showsegments = true)
viz(fig[2,2], ref3, showsegments = true)
fig
```

## TriSubdivision

```@docs
TriSubdivision
```

```@example refinement
grid = CartesianGrid(10, 10)

# refine three times
ref1 = refine(grid, TriSubdivision())
ref2 = refine(ref1, TriSubdivision())
ref3 = refine(ref2, TriSubdivision())

fig = Mke.Figure(size = (800, 800))
viz(fig[1,1], grid, showsegments = true)
viz(fig[1,2], ref1, showsegments = true)
viz(fig[2,1], ref2, showsegments = true)
viz(fig[2,2], ref3, showsegments = true)
fig
```

## Catmull-Clark

```@docs
CatmullClarkRefinement
```

```@example refinement
# define a cube in R^3
points = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
connec = connect.([(1,4,3,2),(5,6,7,8),(1,2,6,5),(3,4,8,7),(1,5,8,4),(2,3,7,6)])
mesh   = SimpleMesh(points, connec)

# refine three times
ref1 = refine(mesh, CatmullClarkRefinement())
ref2 = refine(ref1, CatmullClarkRefinement())
ref3 = refine(ref2, CatmullClarkRefinement())

fig = Mke.Figure(size = (800, 800))
viz(fig[1,1], mesh, showsegments = true)
viz(fig[1,2], ref1, showsegments = true)
viz(fig[2,1], ref2, showsegments = true)
viz(fig[2,2], ref3, showsegments = true)
fig
```

## RegularRefinement

```@docs
RegularRefinement
```

```@example refinement
grid = CartesianGrid(10, 10)

# refine three times
ref1 = refine(grid, RegularRefinement(2, 2))
ref2 = refine(ref1, RegularRefinement(3, 2))
ref3 = refine(ref2, RegularRefinement(2, 3))

fig = Mke.Figure(size = (800, 800))
viz(fig[1,1], grid, showsegments = true)
viz(fig[1,2], ref1, showsegments = true)
viz(fig[2,1], ref2, showsegments = true)
viz(fig[2,2], ref3, showsegments = true)
fig
```

## MaxLengthRefinement

```@docs
MaxLengthRefinement
```

```@example refinement
grid = CartesianGrid(10, 10)

# refine three times
ref1 = refine(grid, MaxLengthRefinement(0.5))
ref2 = refine(ref1, MaxLengthRefinement(0.25))
ref3 = refine(ref2, MaxLengthRefinement(0.1))

fig = Mke.Figure(size = (800, 800))
viz(fig[1,1], grid, showsegments = true)
viz(fig[1,2], ref1, showsegments = true)
viz(fig[2,1], ref2, showsegments = true)
viz(fig[2,2], ref3, showsegments = true)
fig
```
