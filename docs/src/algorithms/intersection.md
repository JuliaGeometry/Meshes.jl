# Intersection

Intersections are implemented for various geometries such
as [`Segment`](@ref), [`Line`](@ref), and [`Box`](@ref):

```@example intersection
using Meshes # hide
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
