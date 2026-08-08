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

### Colors

Colors can be set globaly for the entire geometry as with the triangles above.
If the colors are a vector that is the same size as the number of vertices,
the coloring is per-vertex:

```@example viz
sphere = Sphere((0.,0.,0.), 1.)
mesh = discretize(sphere, RegularDiscretization(10,10))
viz(mesh; color=rand([:red,:green,:blue],nvertices(mesh)))
```

If the number of colors is equal to the number of elements, per-face coloring
is used. In this case, the normals used in the visualization are also
calculated per face, resulting in a faceted appearance, even if the same
color is applied to each face:

```@example viz
sphere = Sphere((0.,0.,0.), 1.)
mesh = discretize(sphere, RegularDiscretization(10,10))
viz(mesh; color=fill(:grey,nelements(mesh)))
```