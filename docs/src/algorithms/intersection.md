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
