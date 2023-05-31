# Primitives

```@example primitives
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example primitives
using Meshes, MeshViz # hide
import WGLMakie as Mke # hide
```

## Abstract

```@docs
Primitive
```

## Concrete

```@docs
Ray
```

```@docs
Line
```

```@docs
BezierCurve
```

```@example primitives
BezierCurve((0.,0.), (1.,0.), (1.,1.)) |> viz
```

```@docs
Plane
```

```@docs
Box
```

```@example primitives
Box((0.,0.,0.), (1.,1.,1.)) |> viz
```

```@docs
Ball
```

```@example primitives
Ball((0.,0.,0.), 1.) |> viz
```

```@docs
Sphere
```

```@example primitives
Sphere((0.,0.,0.), 1.) |> viz
```

```@docs
Disk
```

```@docs
Circle
```

```@docs
Cylinder
```

```@example primitives
Cylinder(1.0) |> viz
```

```@docs
CylinderSurface
```

```@example primitives
CylinderSurface(1.0) |> viz
```

```@docs
Cone
```

```@docs
ConeSurface
```

```@docs
Torus
```

```@example primitives
Torus((0.,0.,0.), (1.,0.,0.), (0.,1.,0.), 0.2) |> viz
```