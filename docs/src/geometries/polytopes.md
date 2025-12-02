# Polytopes

```@example polytopes
using Meshes # hide
import CairoMakie as Mke # hide
```

## Abstract

```@docs
Polytope
Chain
Polygon
Polyhedron
```

## Concrete

### Segment

```@docs
Segment
```

```@example polytopes
Segment((0, 0), (1, 1)) |> viz
```

### Rope

```@docs
Rope
```

```@example polytopes
Rope((0.0, 0.0), (1.0, 0.5), (1.0, 1.0), (2.0, 0.0)) |> viz
```

### Ring

```@docs
Ring
```

```@example polytopes
Ring((0.0, 0.0), (1.0, 0.5), (1.0, 1.0), (2.0, 0.0)) |> viz
```

### Ngon

```@docs
Ngon
```

```@example polytopes
Triangle((0, 0), (1, 0), (0, 1)) |> viz
```

### PolyArea

```@docs
PolyArea
```

```@example polytopes
outer = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
hole1 = [(0.2, 0.2), (0.2, 0.4), (0.4, 0.4), (0.4, 0.2)]
hole2 = [(0.6, 0.2), (0.6, 0.4), (0.8, 0.4), (0.8, 0.2)]
PolyArea([outer, hole1, hole2]) |> viz
```

### Tetrahedron

```@docs
Tetrahedron
```

```@example polytopes
Tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)) |> viz
```

### Hexahedron

```@docs
Hexahedron
```

```@example polytopes
Hexahedron(
  (0, 0, 0),
  (1, 0, 0),
  (1, 1, 0),
  (0, 1, 0),
  (0, 0, 1),
  (1, 0, 1),
  (1, 1, 1),
  (0, 1, 1)
) |> viz
```

### Pyramid

```@docs
Pyramid
```

```@example polytopes
Pyramid(
  (0, 0, 0),
  (1, 0, 0),
  (1, 1, 0),
  (0, 1, 0),
  (0, 0, 1)
) |> viz
```

### Wedge

```@docs
Wedge
```

```@example polytopes
Wedge(
  (0, 0, 0),
  (1, 0, 0),
  (0, 1, 0),
  (0, 0, 1),
  (1, 0, 1),
  (0, 1, 1)
) |> viz
```
