# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    UniformSampling(size, replace=false)

Sample elements uniformly from a given domain/data. Produce a
sample of given `size` with or without replacement depending on
the `replace` option.
"""
struct UniformSampling <: SamplingMethod
  size::Int
  replace::Bool
end

UniformSampling(size::Int) = UniformSampling(size, false)

function sample(object, method::UniformSampling)
  n = nelements(object)
  s = method.size
  r = method.replace
  if s > n && r == false
    @error "invalid sample size for sampling without replacement"
  end
  view(object, sample(1:n, s, replace=r))
end

function sample(triangle::Triangle, method::UniformSampling)
  A, B, C = coordinates.(vertices(triangle))
  function randpoint()
    # sample barycentric coordinates
    u₁, u₂ = rand(2)
    λ₁, λ₂ = 1 - √u₁, u₂ * √u₁
    λ₃     = 1 - λ₁ - λ₂
    Point(λ₁ .* A + λ₂ .* B + λ₃ .* C)
  end
  ivec(randpoint() for _ in 1:method.size)
end
