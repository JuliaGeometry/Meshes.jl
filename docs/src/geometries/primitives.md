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
Line
```

```@docs
Ray
```

```@docs
Plane
```

```@docs
Box
```

```@example primitives
b = Box((0.,0.,0.), (1.,1.,1.))

viz(b)
```

```@docs
Cylinder
```

```@example primitives
c = Cylinder(1.0) # aligned with z axis

viz(c)
```

```@docs
CylinderSurface
```

```@docs
Ball
```

```@docs
Sphere
```

```@example primitives
s = Sphere((0.,0.,0.), 1.)

viz(s)
```

```@docs
Torus
```