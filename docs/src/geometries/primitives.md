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

### ParaboloidSurface

Surface of a paraboloid with the focal axis aligned along the ``z`` direction and clipped up to some radius on the ``xy`` direction.

The equation of the paraboloid is the following:

```math
f(x, y) = \frac{(x - x_0)^2 + (y - y_0)^2}{4f} + z_0\qquad\text{for } x^2 + y^2 < r^2,
```
where ``(x_0, y_0, z_0)`` is the vertex of the parabola, ``f`` is the focal length, and ``r`` is the clip radius. See also <https://en.wikipedia.org/wiki/Paraboloid>.

```@docs
ParaboloidSurface
```

The following example shows a paraboloid together with a disk on the ``xy`` plane; the disk has the same radius used for the parabola. 

```@example primitives
v = Point3(5, 2, 4)
r = 1.0
f = 0.25
par = ParaboloidSurface(v, r, f)
disk = Disk(Plane(v, Vec(0, 0, 1)), r)
viz([par, disk], color = [:green, :gray])
```
