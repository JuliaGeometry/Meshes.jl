# Merging

Geometries and meshes can be [`merge`](@ref)d into a single
geometric object as illustrated in the following example.
The resulting type depends on the combination of input types,
and can be a [`Mesh`](@ref) or [`Multi`](@ref) geometry.

```@docs
merge(::Mesh, ::Mesh)
```

```@example merge
using Meshes # hide
import CairoMakie as Mke # hide

g = CartesianGrid(2, 2)
t = Triangle((3, 0), (4, 0), (3, 1))

m = merge(g, t)
```
