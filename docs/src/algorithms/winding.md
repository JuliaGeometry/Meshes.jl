# Winding

The [`winding`](@ref) number is intimately connected to the
[`sideof`](@ref) function, which is used more often in applications.

```@docs
winding
```

```@example winding
using Meshes # hide

points = [Point(0, 0), Point(0.2, 0.2), Point(2, 1)]
polygon = Triangle((0, 0), (1, 0), (0, 1))

winding(points, boundary(polygon))
```
