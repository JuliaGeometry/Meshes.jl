# Intersection

Intersections are implemented for various geometries and domains
with the ∩ (`\cap`) operator:

```@example intersection
using Meshes

s1 = Segment((0.0,0.0), (1.0,0.0))
s2 = Segment((0.5,0.0), (2.0,0.0))

s1 ∩ s2
```

First, the `intersection` function computes the `Intersection` object,
which holds the `IntersectionType` besides the actual geometry:

```@example intersection
I = intersection(s1, s2)
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
intersection(s1, s2) do I
  if type(I) == Crossing
    return 1
  elseif type(I) == Overlapping
    return 2
  else
    return 3
  end
end
```

```@docs
IntersectionType
Intersection
intersection
intersect(::Geometry, ::Geometry)
```

More generally, when the geometries are not convex nor simple,
it is still possible to know whether or not they have an intersection:

```@docs
intersects
supportfun
```

```@example intersection
outer = Point2[(0,0),(1,0),(1,1),(0,1)]
hole1 = Point2[(0.2,0.2),(0.4,0.2),(0.4,0.4),(0.2,0.4)]
hole2 = Point2[(0.6,0.2),(0.8,0.2),(0.8,0.4),(0.6,0.4)]
poly  = PolyArea(outer, [hole1, hole2])
ball1 = Ball((0.5,0.5), 0.05)
ball2 = Ball((0.3,0.3), 0.05)
ball3 = Ball((0.7,0.3), 0.05)
ball4 = Ball((0.3,0.3), 0.15)

intersects(poly, ball1)
```

```@example intersection
intersects(poly, ball2)
```

```@example intersection
intersects(poly, ball3)
```

```@example intersection
intersects(poly, ball4)
```
