# Union

```@example union
using Meshes # hide
import CairoMakie as Mke # hide
```

Unions are implemented for some geometries (e.g. polygons)
with the `∪` (`\cup`) operator:

```@docs
union(::Geometry, ::Geometry)
```

```@example union
poly₁ = Pgon([
  [(8, 0), (4, 8), (2, 8), (-2, 0), (0, 0), (1, 2), (5, 2), (6, 0)],
  [(4, 4), (2, 4), (3, 6)]
])
poly₂ = Pgon((0, 1), (6, 1), (6, 7), (4, 7), (4, 3), (2, 3), (2, 7), (0, 7))

poly = poly₁ ∪ poly₂

viz(poly₁)
viz!(poly₂, color = "black", alpha = 0.2)
viz!(boundary(poly), color = "red", segmentsize = 3)
Mke.current_figure()
```
