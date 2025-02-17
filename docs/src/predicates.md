# Predicates

```@example intersects
using Meshes # hide
import CairoMakie as Mke # hide
```

This section lists predicates that can be used to check properties of geometric
objects, both of themselves and relative to other geometric objects. 

One important note to make is that these predicates are not necessarily exact.
For example, rather than checking if a point `p` is exactly in a sphere of radius
`r` centered at `c`, we check if `norm(p-c) ≈ r` with an absolute tolerance depending
on the point type, so `p` might be slightly outside the sphere but still be considered
as being inside. This absolute tolerance can be adjusted in specific scopes as discussed
in the [Tolerances](tolerances.md) section.

Robust predicates are often expensive to apply and approximations typically suffice.
If needed, consider [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl) or 
[AdaptivePredicates.jl](https://github.com/JuliaGeometry/AdaptivePredicates.jl).

## isparametrized

```@docs
isparametrized
paramdim
```

## iscurve

```@docs
iscurve
```

## issurface

```@docs
issurface
```

## issolid

```@docs
issolid
```

## isperiodic

```@docs
isperiodic
```

## issimplex

```@docs
issimplex
```

## isclosed

```@docs
isclosed
```

## isconvex

```@docs
isconvex
```

## issimple

```@docs
issimple
```

## hasholes

```@docs
hasholes
```

## point₁ ≤ point₂

```@docs
Base.:<(::Point, ::Point)
Base.:>(::Point, ::Point)
Base.:≤(::Point, ::Point)
Base.:≥(::Point, ::Point)
```

## point₁ ⪯ point₂

```@docs
≺(::Point, ::Point)
≻(::Point, ::Point)
⪯(::Point, ::Point)
⪰(::Point, ::Point)
```

## point ∈ geometry

```@docs
Base.in(::Point, ::Geometry)
```

## geometry₁ ⊆ geometry₂

```@docs
Base.issubset(::Geometry, ::Geometry)
```

## intersects

```@docs
intersects
supportfun
```

```@example intersects
outer = [(0,0),(1,0),(1,1),(0,1)]
hole1 = [(0.2,0.2),(0.4,0.2),(0.4,0.4),(0.2,0.4)]
hole2 = [(0.6,0.2),(0.8,0.2),(0.8,0.4),(0.6,0.4)]
poly  = PolyArea([outer, hole1, hole2])
ball1 = Ball((0.5,0.5), 0.05)
ball2 = Ball((0.3,0.3), 0.05)
ball3 = Ball((0.7,0.3), 0.05)
ball4 = Ball((0.3,0.3), 0.15)

intersects(poly, ball1)
```

```@example intersects
intersects(poly, ball2)
```

```@example intersects
intersects(poly, ball3)
```

```@example intersects
intersects(poly, ball4)
```

## iscollinear 

```@docs 
iscollinear 
``` 

## iscoplanar 

```@docs 
iscoplanar 
```
