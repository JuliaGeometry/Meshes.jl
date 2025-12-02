# Primitives

```@example primitives
using Meshes # hide
import CairoMakie as Mke # hide
```

## Abstract

```@docs
Primitive
```

## Concrete

### Point

```@docs
Point
```

```@example primitives
rand(Point, 100) |> viz
```

```@docs
to(::Point)
-(::Point, ::Point)
+(::Point, ::Vec)
-(::Point, ::Vec)
```

### Ray

```@docs
Ray
```

```@example primitives
Ray((1, 2, 3), (1, 0, 1)) |> viz
```

### Line

```@docs
Line
```

```@example primitives
Line((1, 2), (3, 5)) |> viz
```

### BezierCurve

```@docs
BezierCurve
```

```@example primitives
BezierCurve((0, 0), (1, 0), (1, 1)) |> viz
```

### ParametrizedCurve

```@docs
ParametrizedCurve
```

```@example primitives
ParametrizedCurve(t -> Point(cos(t), sin(t), 0.2t), (0, 4π)) |> viz
```

### Plane

```@docs
Plane
```

### Box

```@docs
Box
```

```@example primitives
Box((0, 0, 0), (1, 1, 1)) |> viz
```

### Ball/Sphere

```@docs
Ball
Sphere
```

```@example primitives
Ball((0, 0, 0), 1) |> viz
```

### Ellipsoid

```@docs
Ellipsoid
```

```@example primitives
Ellipsoid((3, 2, 1)) |> viz
```

### Disk/Circle

```@docs
Disk
Circle
```

```@example primitives
Disk(Plane((0, 0, 0), (1, 0, 1)), 1) |> viz
```

### Cylinder/CylinderSurface

```@docs
Cylinder
CylinderSurface
```

```@example primitives
Cylinder(1) |> viz
```

### Cone/ConeSurface

```@docs
Cone
ConeSurface
```

```@example primitives
Cone(
  Disk(Plane((0, 0, 0), (0, 0, 1)), 1),
  (0, 0, 1)
) |> viz
```

### Frustum/FrustumSurface

```@docs
Frustum
FrustumSurface
```

```@example primitives
Frustum(
  Disk(Plane((0, 0, 0), (0, 0, 1)), 2),
  Disk(Plane((0, 0, 10), (0, 0, 1)), 1)
) |> viz
```

### Torus

```@docs
Torus
```

```@example primitives
Torus((0, 0, 0), (1, 0, 0), (0, 1, 0), 0.2) |> viz
```

### ParaboloidSurface

```@docs
ParaboloidSurface
```

```@example primitives
ParaboloidSurface((0, 0, 0), 1, 0.25) |> viz
```
