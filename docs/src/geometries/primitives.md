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
rand(Point3, 100) |> viz
```

```@docs
coordinates(::Point)
-(::Point, ::Point)
+(::Point, ::Vec)
-(::Point, ::Vec)
```

### Ray

```@docs
Ray
```

### Line

```@docs
Line
```

### BezierCurve

```@docs
BezierCurve
```

```@example primitives
BezierCurve((0.,0.), (1.,0.), (1.,1.)) |> viz
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
Box((0.,0.,0.), (1.,1.,1.)) |> viz
```

### Ball

```@docs
Ball
```

```@example primitives
Ball((0.,0.,0.), 1.) |> viz
```

### Sphere

```@docs
Sphere
```

```@example primitives
Sphere((0.,0.,0.), 1.) |> viz
```

### Disk

```@docs
Disk
```

### Circle

```@docs
Circle
```

### Cylinder

```@docs
Cylinder
```

```@example primitives
Cylinder(1.0) |> viz
```

### CylinderSurface

```@docs
CylinderSurface
```

```@example primitives
CylinderSurface(1.0) |> viz
```

### Cone

```@docs
Cone
```

### ConeSurface

```@docs
ConeSurface
```

### Frustum

```@docs
Frustum
```

### FrustumSurface

```@docs
FrustumSurface
```

### Torus

```@docs
Torus
```

```@example primitives
Torus((0.,0.,0.), (1.,0.,0.), (0.,1.,0.), 0.2) |> viz
```
