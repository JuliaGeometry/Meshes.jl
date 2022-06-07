# Polytopes

```@example polytopes
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example polytopes
using Meshes, MeshViz # hide
import WGLMakie as Mke # hide
```

## Abstract

```@docs
Polytope
Polygon
Polyhedron
```

## Concrete

```@docs
Segment
```

```@example polytopes
s = Segment((0., 0.), (1.,1.))

viz(s)
```

```@docs
Ngon
```

```@example polytopes
t = Triangle((0.,0.), (1.,0.), (0.,1.))

viz(t)
```

```@docs
Chain
```

```@example polytopes
c = Chain((0.,0.), (1.,0.5), (1.,1.), (2.,0.))

viz(c)
```

```@docs
PolyArea
```

```@example polytopes
outer = [(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0),(0.0,0.0)]
hole1 = [(0.2,0.2),(0.4,0.2),(0.4,0.4),(0.2,0.4),(0.2,0.2)]
hole2 = [(0.6,0.2),(0.8,0.2),(0.8,0.4),(0.6,0.4),(0.6,0.2)]
poly  = PolyArea(outer, [hole1, hole2])

viz(poly)
```

```@docs
Tetrahedron
```

```@example polytopes
t = Tetrahedron([(0,0,0),(1,0,0),(0,1,0),(0,0,1)])

viz(t)
```

```@docs
Pyramid
```

```@docs
Hexahedron
```

```@example polytopes
h = Hexahedron([(0,0,0),(1,0,0),(1,1,0),(0,1,0),
                (0,0,1),(1,0,1),(1,1,1),(0,1,1)])

viz(h)
```
