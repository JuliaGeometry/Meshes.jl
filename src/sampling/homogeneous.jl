# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HomogeneousSampling(size, [weights])

Generate sample of given `size` from geometric object
according to a homogeneous density. Optionally, provide `weights`
to specify custom sampling weights for the elements of a domain.
"""
struct HomogeneousSampling{W} <: ContinuousSamplingMethod
  size::Int
  weights::W
end

HomogeneousSampling(size::Int) = HomogeneousSampling(size, nothing)

function sample(rng::AbstractRNG, d::Domain, method::HomogeneousSampling)
  size = method.size
  weights = isnothing(method.weights) ? map(measure, d) : method.weights

  # sample elements with weights
  w = WeightedSampling(size, map(ustrip, weights), replace=true)

  # within each element sample a single point
  h = HomogeneousSampling(1)

  (first(sample(rng, e, h)) for e in sample(rng, d, w))
end

function sample(rng::AbstractRNG, g::Geometry, method::HomogeneousSampling)
  if isparametrized(g)
    randpoint() = g(rand(rng, numtype(lentype(g)), paramdim(g))...)
    (randpoint() for _ in 1:(method.size))
  else
    sample(rng, discretize(g), method)
  end
end

# --------------
# SPECIAL CASES
# --------------

function sample(rng::AbstractRNG, t::Triangle, method::HomogeneousSampling)
  function randpoint()
    # sample barycentric coordinates
    u₁, u₂ = rand(rng, numtype(lentype(t)), paramdim(t))
    λ₁, λ₂ = 1 - √u₁, u₂ * √u₁
    t(λ₁, λ₂)
  end
  (randpoint() for _ in 1:(method.size))
end

function sample(rng::AbstractRNG, t::Tetrahedron, method::HomogeneousSampling)
  function randpoint()
    # sample barycentric coordinates
    u₁, u₂, u₃ = rand(rng, numtype(lentype(t)), paramdim(t))
    λ₁ = 1 - ∛u₁
    λ₂ = (1 - λ₁) * (1 - √u₂)
    λ₃ = (1 - λ₁) * √u₂ * u₃
    t(λ₁, λ₂, λ₃)
  end
  (randpoint() for _ in 1:(method.size))
end

sample(rng::AbstractRNG, b::Ball, method::HomogeneousSampling) = _sample(rng, b, Val(paramdim(b)), method)

function _sample(rng::AbstractRNG, b::Ball, ::Val{2}, method::HomogeneousSampling)
  function randpoint()
    u₁, u₂ = rand(rng, numtype(lentype(b)), paramdim(b))
    λ₁, λ₂ = √u₁, u₂
    b(λ₁, λ₂)
  end
  (randpoint() for _ in 1:(method.size))
end

function _sample(rng::AbstractRNG, b::Ball, ::Val{3}, method::HomogeneousSampling)
  function randpoint()
    T = numtype(lentype(b))
    u₁, u₂, u₃ = rand(rng, T, paramdim(b))
    λ₁ = ∛u₁
    λ₂ = acos(1 - 2u₂) / T(π)
    λ₃ = u₃
    b(λ₁, λ₂, λ₃)
  end
  (randpoint() for _ in 1:(method.size))
end
