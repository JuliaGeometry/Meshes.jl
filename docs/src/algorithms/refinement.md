# Refinement

```@docs
refine
RefinementMethod
```

## Catmull-Clark

```@docs
CatmullClark
```

```@example
using Meshes # hide
using Plots # hide
gr(fillcolor=false) # hide

# define a cube in R^3
points = Point3[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
connec = connect.([(1,4,3,2),(5,6,7,8),(1,2,6,5),(3,4,8,7),(1,5,8,4),(2,3,7,6)])
mesh   = SimpleMesh(points, connec)

# refine three times
ref1 = refine(mesh, CatmullClark())
ref2 = refine(ref1, CatmullClark())
ref3 = refine(ref2, CatmullClark())

p0 = plot(mesh, title="original")
p1 = plot(ref1, title="refine 1")
p2 = plot(ref2, title="refine 2")
p3 = plot(ref3, title="refine 3")

plot(p0, p1, p2, p3, size=(800,800))
```
