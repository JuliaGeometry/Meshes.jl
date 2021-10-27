# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HomogeneousSampling(size)

Generate sample of given `size` from geometric object
according to a homogeneous density.
"""
struct HomogeneousSampling <: ContinuousSamplingMethod
  size::Int
end

function sample(rng::AbstractRNG, Ω::DomainOrData,
                method::HomogeneousSampling)
  size    = method.size
  weights = measure.(Ω)

  # sample elements with weights proportial to measure
  w = WeightedSampling(size, weights, replace=true)

  # within each element sample a single point
  h = HomogeneousSampling(1)

  (first(sample(rng, e, h)) for e in sample(rng, Ω, w))
end

function sample(rng::AbstractRNG, triangle::Triangle{Dim,T},
                method::HomogeneousSampling) where {Dim,T}
  A, B, C = coordinates.(vertices(triangle))
  function randpoint()
    # sample barycentric coordinates
    u₁, u₂ = rand(rng, T, 2)
    λ₁, λ₂ = 1 - √u₁, u₂ * √u₁
    λ₃     = 1 - λ₁ - λ₂
    Point(λ₁ .* A + λ₂ .* B + λ₃ .* C)
  end
  (randpoint() for _ in 1:method.size)
end

function sample(rng::AbstractRNG, s::Segment{Dim,T},
                method::HomogeneousSampling) where {Dim,T}
  A, B = vertices(s)
  function randpoint()
    t = rand(rng, T)
    A + t * (B - A)
  end
  (randpoint() for _ in 1:method.size)
end
