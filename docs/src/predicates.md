# Predicates

```@example intersects
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example intersects
using Meshes # hide
import WGLMakie as Mke # hide
```

## isparametrized

```@docs
isparametrized
paramdim
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
outer = Point2[(0,0),(1,0),(1,1),(0,1)]
hole1 = Point2[(0.2,0.2),(0.4,0.2),(0.4,0.4),(0.2,0.4)]
hole2 = Point2[(0.6,0.2),(0.8,0.2),(0.8,0.4),(0.6,0.4)]
poly  = PolyArea(outer, [hole1, hole2])
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
