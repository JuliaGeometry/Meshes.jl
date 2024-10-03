"""
	FibonacciSampling(n)

Generate samples on a sphere using the Fibonacci lattice method.
The parameter `n` specifies the number of points to generate.
In the regular Fibonacci lattice method, the number ϕ is the golden ratio ((1 + √5)/2),
but different numbers can be used, preferably irrational.

## Example

Sample a sphere with 1000 points:

```julia
sample(Sphere((0,0,0), 1), FibonacciSampling(100))

# sample using π instead of the golden ratio
sample(Box((0,0),(1,1)), FibonacciSampling(100,π))
"""
struct FibonacciSampling{T<:Real} <: ContinuousSamplingMethod
  n::Int
  ϕ::T

  function FibonacciSampling(n::Int, ϕ::T) where {T<:Real}
    if n ≤ 0
      throw(ArgumentError("Size must be positive"))
    end
    new{T}(n, ϕ)
  end
end

# default to golden ratio
FibonacciSampling(n::Int) = FibonacciSampling(n, (1 + √5) / 2)

function sample(geom::Geometry, method::FibonacciSampling)
  if paramdim(geom) != 2
    throw(ArgumentError("Fibonacci sampling only defined for 2D geometries"))
  end
  f = _distortion(geom)
  function point(i)
    u, v = mod(i / method.ϕ, 1), i / (method.n - 1)
    geom(f(u, v)...)
  end

  (point(i) for i in 0:(method.n - 1))
end

_distortion(g) = (u, v) -> (u, v)
_distortion(d::Disk) = (u, v) -> (√u, v)
_distortion(b::Ball{𝔼{2}}) = (u, v) -> (√u, v)
_distortion(s::Sphere{𝔼{3}}) = (u, v) -> (acos(1 - 2v) / π, u)
_distortion(b::Ball{🌐}) = (u, v) -> (√u, v)