# Visualization

```@example viz
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example viz
using Meshes # hide
import WGLMakie as Mke # hide
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
triangles = rand(Triangle{2,Float64}, 10)

viz(triangles, color = 1:10)
```

## Domains

Alternatively, we can visualize domains with topological information
such as [`Mesh`](@ref) and show facets efficiently:

```@example viz
grid = CartesianGrid(10, 10, 10)

viz(grid, showfacets = true, facetcolor = :teal)
```