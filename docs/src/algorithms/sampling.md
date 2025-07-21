# Sampling

```@example sampling
using Meshes # hide
import CairoMakie as Mke # hide
```

```@docs
sample(::Any, ::SamplingMethod)
SamplingMethod
DiscreteSamplingMethod
ContinuousSamplingMethod
```

## Discrete sampling

### UniformSampling

```@docs
UniformSampling
```

```@example sampling
grid = CartesianGrid(20, 20)

# uniform sampling without replacement
sampler = UniformSampling(100, replace=false)
blocks  = sample(grid, sampler)

viz(blocks)
```

### WeightedSampling

```@docs
WeightedSampling
```

```@example sampling
grid = CartesianGrid(20, 20)

# upper blocks are 10x more likely
weights = [fill(1, 200); fill(10, 200)]

# weighted sampling without replacement
sampler = WeightedSampling(100, weights, replace=false)
blocks  = sample(grid, sampler)

viz(blocks)
```

### BallSampling

```@docs
BallSampling
```

```@example sampling
grid = CartesianGrid(20, 20)

# sample blocks that are apart by a given radius
sampler = BallSampling(5.0)
blocks  = sample(grid, sampler)

viz(blocks)
```

## Continuous sampling

### RegularSampling

```@docs
RegularSampling
```

```@example sampling
grid = CartesianGrid(20, 20)

# sample points regularly
sampler = RegularSampling(20, 30)
points  = sample(grid, sampler) |> collect

viz(points)
```

### HomogeneousSampling

```@docs
HomogeneousSampling
```

```@example sampling
grid = CartesianGrid(20, 20)

# sample points homogeneously
sampler = HomogeneousSampling(100)
points  = sample(grid, sampler) |> collect

viz(points)
```

### MinDistanceSampling

```@docs
MinDistanceSampling
```

```@example sampling
grid = CartesianGrid(20, 20)

# sample points that are apart by a given radius
sampler = MinDistanceSampling(3.0)
points  = sample(grid, sampler) |> collect

viz(points)
```

### FibonacciSampling
```@docs
FibonacciSampling
```

```@example sampling
sphere = Sphere((0.,0.,0.), 1.)

# sample points using the Fibonacci lattice method
sampler  = FibonacciSampling(100)
points  = sample(sphere, sampler) |> collect

viz(points)
```

### AdaptiveSampling

```@docs
AdaptiveSampling
```

```@example sampling
# Sample between 10 and 100 points
sampler = AdaptiveSampling(; min_points=10.0, max_points=100.0)
crv = ParametrizedCurve((t -> Point((sin(t * 2π) * t, cos(t * 2π) * t))), (0, 1))
points  = sample(crv, sampler) |> collect

viz(points)
```

