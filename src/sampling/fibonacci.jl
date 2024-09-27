"""
	FibonacciSampling(n)

Generate samples on a sphere using the Fibonacci lattice method.
The parameter `n` specifies the number of points to generate.

## Example

Sample a sphere with 1000 points:

```julia
sample(Sphere((0,0,0), 1), FibonacciSampling(100))
"""
struct FibonacciSampling <: SamplingMethod
  size::Int
  function FibonacciSampling(size::Int)
    if size ≤ 0
      throw(ArgumentError("Size must be positive"))
    end
    new(size)
  end
end

function sample(::AbstractRNG, geom::Union{Sphere,Ball{𝔼{3}}}, method::FibonacciSampling)

  # Compute x,y samples from regular sphere of radius 1
  x, y, z = fibonnaci_sample_sphere(method.size)

  # Translating and scaling points based on sphere properties 
  ## Scale by radius
  scale_radius = Scale(map(x -> ustrip(geom.radius), 1:3)...)
  ## Translate by center
  translate_center = Translate(to(geom.center)...)

  ## Apply scaling and translation to points
  Iterators.map(zip(x, y, z)) do (x, y, z)
    Point(x, y, z) |> translate_center ∘ scale_radius
  end
end

function sample(::AbstractRNG, geom::Disk, method::FibonacciSampling)

  # Compute x,y samples from regular 2D disk of radius 1
  x, y = fibonnaci_sample_disk(method.size)

  # Scaling, rotating and translating points based on plane for disk
  scale_radius = Scale(map(x -> ustrip(geom.radius), 1:3)...)
  translate_center = Translate(to(geom.plane.p)...)
  rotation_plane = Rotate(Vec(0, 0, 1), normal(geom.plane))

  Iterators.map(zip(x, y)) do (x, y)
    Point(x, y, 0) |> translate_center ∘ rotation_plane ∘ scale_radius
  end
end
function sample(::AbstractRNG, geom::Ball{𝔼{2}}, method::FibonacciSampling)

  # Compute x,y samples from regular 2D disk of radius 1
  x, y = fibonnaci_sample_disk(method.size)

  # Scaling, rotating and translating points based on plane for disk
  scale_radius = Scale(map(x -> ustrip(geom.radius), 1:2)...)
  translate_center = Translate(to(geom.center)...)

  Iterators.map(zip(x, y)) do (x, y)
    Point(x, y) |> translate_center ∘ scale_radius
  end
end

"""
    fibonnaci_sample_sphere(n::Int)

Generate `n` sample points uniformly distributed on a sphere using the Fibonacci lattice method.

# Arguments
- `n::Int`: The number of sample points to generate.

# Returns
- `x::Vector{Float64}`: The x-coordinates of the sample points.
- `y::Vector{Float64}`: The y-coordinates of the sample points.
- `z::Vector{Float64}`: The z-coordinates of the sample points.
"""
function fibonnaci_sample_sphere(n::Int)
  # Set the parameters
  goldenRatio = (1 + sqrt(5)) / 2
  i = collect(0:(n - 1))

  # Compute theta and phi
  theta = 2 * pi * mod.(i / goldenRatio, 1)
  phi = acos.(1 .- 2 * i / n)

  # Compute x, y, z coordinates
  x = cos.(theta) .* sin.(phi)
  y = sin.(theta) .* sin.(phi)
  z = cos.(phi)

  return x, y, z
end

"""
    fibonnaci_sample_disk(n::Int)

Generate `n` sample points uniformly distributed on a disk using the Fibonacci lattice method.

# Arguments
- `n::Int`: The number of sample points to generate.

# Returns
- `x::Vector{Float64}`: The x-coordinates of the sample points.
- `y::Vector{Float64}`: The y-coordinates of the sample points.
"""
function fibonnaci_sample_disk(n::Int)
  # Set the parameters
  goldenRatio = (1 + sqrt(5)) / 2
  i = collect(0:(n - 1))

  # Compute theta and phi
  theta = 2 * pi * mod.(i / goldenRatio, 1)
  r = @. √(i / n)
  x = @. cos(theta) * r
  y = @. sin(theta) * r

  return x, y
end