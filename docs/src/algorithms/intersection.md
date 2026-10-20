# Intersection

```@example intersection
using Meshes # hide
import CairoMakie as Mke # hide
```

Intersections are implemented for various geometries and domains
with the `∩` (`\cap`) operator:

```@docs
intersect(::Geometry, ::Geometry)
```

```@example intersection
poly₁ = Pgon([
  [(8, 0), (4, 8), (2, 8), (-2, 0), (0, 0), (1, 2), (5, 2), (6, 0)],
  [(4, 4), (2, 4), (3, 6)]
])
poly₂ = Pgon((0, 1), (6, 1), (6, 7), (4, 7), (4, 3), (2, 3), (2, 7), (0, 7))

poly = poly₁ ∩ poly₂

viz(poly₁)
viz!(poly₂, color = "black", alpha = 0.2)
viz!(boundary(poly), color = "red", segmentsize = 3)
Mke.current_figure()
```

```@example intersection
seg₁ = Segment((0.0, 0.0), (1.0, 0.0))
seg₂ = Segment((0.5, 0.0), (2.0, 0.0))

seg₁ ∩ seg₂
```

First, the [`intersection`](@ref) function computes the [`Intersection`](@ref)
object, which holds the [`IntersectionType`](@ref) besides the actual geometry:

```@example intersection
I = intersection(seg₁, seg₂)
```

This object supports two methods `type` and `get` to retrieve
the underlying information:

```@example intersection
type(I)
```

```@example intersection
get(I)
```

For performance-sensitive code, it is recommended to use the `intersection`
method with three arguments, including a function to reduce the number of
output types.

In the example below, we use the do syntax to restrict our attention to a
subset of intersection types and to make the return type and `Int` value
in all cases:

```@example intersection
intersection(seg₁, seg₂) do I
  if type(I) == Crossing
    return 1
  elseif type(I) == Overlapping
    return 2
  else
    return 3
  end
end
```

In summary, the following functions and types are related to the `∩` operator:

```@docs
intersection
Intersection
IntersectionType
```
