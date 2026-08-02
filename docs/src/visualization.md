# Visualization

```@example viz
using Meshes # hide
using CoordRefSystems # hide
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
triangles = rand(Triangle, 10, crs=Cartesian2D)

viz(triangles, color = 1:10)
```

## Domains

Alternatively, we can visualize domains with topological information
such as [`Mesh`](@ref) and show facets efficiently:

```@example viz
grid = CartesianGrid(10, 10, 10)

viz(grid, showsegments = true, segmentcolor = :teal)
```

### Normals

By default, Makie will compute per-vertex normals, resulting in a smoothed
representation of e.g. a discretized sphere:

```@example viz
sphere = Sphere((0.,0.,0.), 1.)
mesh = discretize(sphere, RegularDiscretization(10,10))
viz(mesh)
```

If this smoothing is not desired, the `normals` keyword can be used to set the
normals that will be used by Makie when rendering. If the number of normals
equals the number of mesh elements then per-face normals are assumed, if it
matches the number of vertices they are used as per-vertex normals.

```@example viz
sphere = Sphere((0.,0.,0.), 1.)
mesh = discretize(sphere, RegularDiscretization(10,10))
viz(mesh; normals=normal.(mesh))
```