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
    (_homogeneouspoint(rng, g) for _ in 1:(method.size))
  else
    sample(rng, discretize(g), method)
  end
end

_homogeneouspoint(rng, g) = g(rand(rng, numtype(lentype(g)), paramdim(g))...)

function _homogeneouspoint(rng, t::Triangle)
  u₁, u₂ = rand(rng, numtype(lentype(t)), paramdim(t))
  λ₁, λ₂ = 1 - √u₁, u₂ * √u₁
  t(λ₁, λ₂)
end

function _homogeneouspoint(rng, t::Tetrahedron)
  u₁, u₂, u₃ = rand(rng, numtype(lentype(t)), paramdim(t))
  λ₁ = 1 - ∛u₁
  λ₂ = (1 - λ₁) * (1 - √u₂)
  λ₃ = (1 - λ₁) * √u₂ * u₃
  t(λ₁, λ₂, λ₃)
end

function _homogeneouspoint(rng, b::Ball)
  d = paramdim(b)
  T = numtype(lentype(b))
  u = rand(rng, T, d)
  if d == 1
    b(u...)
  elseif d == 2
    λ₁ = √u[1]
    λ₂ = u[2]
    b(λ₁, λ₂)
  elseif d == 3
    λ₁ = ∛u[1]
    λ₂ = acos(1 - 2u[2]) / T(π)
    λ₃ = u[3]
    b(λ₁, λ₂, λ₃)
  end
end
