# Sideof

The [`sideof`](@ref) function can be used to efficiently query the side
of multiple points with respect to a given geometry or mesh.

```@docs
SideType
sideof(::Point, ::Line)
sideof(::Point, ::Ring)
sideof(::Point, ::Mesh)
```

```@example sideof
using Meshes # hide

sideof(Point(0, 0), Line((1, 0), (1, 1)))
```

```@example sideof
points = [Point(0, 0), Point(0.2, 0.2), Point(2, 1)]
polygon = Triangle((0, 0), (1, 0), (0, 1))

sideof(points, boundary(polygon))
```
