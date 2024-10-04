# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
	FibonacciSampling(n, ϕ = (1 + √5)/2)

Generate `n` Fibonacci points with parameter `ϕ`.

The golden ratio is used as the default value of `ϕ`,
but other irrational numbers can be used.

See <https://observablehq.com/@meetamit/fibonacci-lattices>
and <https://www.johndcook.com/blog/2023/08/12/fibonacci-lattice>.
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

FibonacciSampling(n::Int) = FibonacciSampling(n, (1 + √5) / 2)

function sample(geom::Geometry, method::FibonacciSampling)
  if paramdim(geom) != 2
    throw(ArgumentError("Fibonacci sampling only defined for 2D geometries"))
  end

  fib = _fibmap(geom)

  function point(i)
    u = mod(i / method.ϕ, 1)
    v = i / (method.n - 1)
    geom(fib(u, v)...)
  end

  (point(i) for i in 0:(method.n - 1))
end

_fibmap(g) = (u, v) -> (u, v)
_fibmap(d::Disk) = (u, v) -> (√u, v)
_fibmap(b::Ball{𝔼{2}}) = (u, v) -> (√u, v)
_fibmap(b::Ball{🌐}) = (u, v) -> (√u, v)
_fibmap(s::Sphere{𝔼{3}}) = (u, v) -> (acos(1 - 2v) / π, u)
