# Intersection

Intersections are implemented for various geometries such
as [`Segment`](@ref), [`Line`](@ref), and [`Box`](@ref):

```@example intersection
using Meshes

s1 = Segment((0.0,0.0), (1.0,0.0))
s2 = Segment((0.5,0.0), (2.0,0.0))

s1 ∩ s2
```

The algorithm first identifies the type of intersection
using `intersecttype`:

```@example intersection
I = intersecttype(s1, s2)
```

and then retrieves the actual value using `get`:

```@example intersection
get(I)
```

More generally, when the geometries are not convex nor simple,
it is still possible to know whether or not they have an intersection:

```@docs
hasintersect
supportfun
```

```@example intersection
outer = Point2[(0,0),(1,0),(1,1),(0,1),(0,0)]
hole1 = Point2[(0.2,0.2),(0.4,0.2),(0.4,0.4),(0.2,0.4),(0.2,0.2)]
hole2 = Point2[(0.6,0.2),(0.8,0.2),(0.8,0.4),(0.6,0.4),(0.6,0.2)]
poly  = PolyArea(outer, [hole1, hole2])
ball1 = Ball((0.5,0.5), 0.05)
ball2 = Ball((0.3,0.3), 0.05)
ball3 = Ball((0.7,0.3), 0.05)
ball4 = Ball((0.3,0.3), 0.15)

hasintersect(poly, ball1)
```

```@example intersection
hasintersect(poly, ball2)
```

```@example intersection
hasintersect(poly, ball3)
```

```@example intersection
hasintersect(poly, ball4)
```
