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
