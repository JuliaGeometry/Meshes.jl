# Calculus

```@example calculus
using Meshes # hide
using Unitful # hide
```

Calculus with geometries is possible thanks to careful
parametrizations and high-level interfaces such as
[DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl)
and [IntegrationInterface.jl](https://github.com/pablosanjose/IntegrationInterface.jl).

Consider the following quadrangle for illustration purposes:

```@example calculus
q = Quadrangle((0, 0, 0), (2, 0, 0), (2, 1, 0), (0, 1, 0))
```

## Differentiation

```@docs
derivative
```

```@example calculus
# derivative at center point along first axis
derivative(q, (0.5, 0.5), 1)
```

```@docs
jacobian
```

```@example calculus
# Jacobian at center (i.e., derivative along both axes)
jacobian(q, (0.5, 0.5))
```

```@docs
differential
```

```@example calculus
# differential element (i.e., infinitesimal area)
differential(q, (0.5, 0.5))
```

## Integration

```@docs
integral
```

```@example calculus
# integral of constant 1 over quadrangle gives area
integral(p -> 1, q)
```

```@example calculus
# unitful integrand is also supported
integral(p -> 1u"W/m^2", q)
```

```@example calculus
# less trivial integrand in terms of Cartesian coordinates
integral(q) do p
  x, y, z = to(p)
  x + 2y + 3z # meter units
end
```

```@docs
localintegral
```

```@example calculus
# local integral in terms of parametric coordinates in [0, 1]²
localintegral((u, v) -> u^2 + 3v, q)
```

All these methods are implemented for curved geometries
(e.g., spherical geometry):

```@example calculus
using CoordRefSystems

# curved quadrangle over the globe
q = Quadrangle(Point(LatLon(0, 0)), Point(LatLon(0, 90)), Point(LatLon(80, 90)), Point(LatLon(80, 0)))

# Jacobian components as rays over the quadrangle
uv = Iterators.product(0:0.1:1, 0:0.1:1)
rs = [Ray.(q(u, v), 0.05 .* jacobian(q, (u, v))) for (u, v) in uv]

# split rays into u and v components
rᵤ = first.(rs) |> vec
rᵥ = last.(rs) |> vec

viz(q, color="teal")
viz!(rᵤ, color="gray80")
viz!(rᵥ, color="gray80")

# Bezier curve over the globe
c = BezierCurve(Point(LatLon(0, 0)), Point(LatLon(40, 90)), Point(LatLon(80, 0)))

# Jacobian components as rays over the curve
ts = 0:0.1:1
rs = [Ray.(c(t), 0.05 .* jacobian(c, (t,))) for t in ts]

# extract t component
rₜ = first.(rs) |> vec

viz!(c, color="yellow", segmentsize=4)
viz!(rₜ, color="magenta")

Mke.current_figure()
```
