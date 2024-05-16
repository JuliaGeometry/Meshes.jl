# Visualization

```@example viz
using Meshes # hide
import CairoMakie as Mke # hide
```

The package exports a single [`viz`](@ref) command that
can be used to add objects to the scene with a consistent
set of options.

```@docs
viz
viz!
```

## Geometries

We can visualize a single geometry or multiple geometries in a vector:

```@example viz
triangles = [rand(Triangle{2}) for _ in 1:10]

viz(triangles, color = 1:10)
```

## Domains

Alternatively, we can visualize domains with topological information
such as [`Mesh`](@ref) and show facets efficiently:

```@example viz
grid = CartesianGrid(10, 10, 10)

viz(grid, showsegments = true, segmentcolor = :teal)
```