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

function sample(rng::AbstractRNG, Ω::DomainOrData, method::HomogeneousSampling)
  size = method.size
  weights = isnothing(method.weights) ? measure.(Ω) : method.weights

  # sample elements with weights
  w = WeightedSampling(size, weights, replace=true)

  # within each element sample a single point
  h = HomogeneousSampling(1)

  (first(sample(rng, e, h)) for e in sample(rng, Ω, w))
end

function sample(rng::AbstractRNG, geom::Geometry{Dim,T}, method::HomogeneousSampling) where {Dim,T}
  if isparametrized(geom)
    randpoint() = geom(rand(rng, T, paramdim(geom))...)
    (randpoint() for _ in 1:(method.size))
  else
    sample(rng, discretize(geom), method)
  end
end

# --------------
# SPECIAL CASES
# --------------

function sample(rng::AbstractRNG, triangle::Triangle{Dim,T}, method::HomogeneousSampling) where {Dim,T}
  function randpoint()
    # sample barycentric coordinates
    u₁, u₂ = rand(rng, T, 2)
    λ₁, λ₂ = 1 - √u₁, u₂ * √u₁
    triangle(λ₁, λ₂)
  end
  (randpoint() for _ in 1:(method.size))
end

function sample(rng::AbstractRNG, tetrahedron::Tetrahedron{Dim,T}, method::HomogeneousSampling) where {Dim,T}
  @error "not implemented"
end

function sample(rng::AbstractRNG, ball::Ball{2,T}, method::HomogeneousSampling) where {T}
  function randpoint()
    u₁, u₂ = rand(rng, T, 2)
    ball(√u₁, u₂)
  end
  (randpoint() for _ in 1:(method.size))
end

function sample(rng::AbstractRNG, ball::Ball{3,T}, method::HomogeneousSampling) where {T}
  function randpoint()
    u₁, u₂, u₃ = rand(rng, T, 3)
    ball(∛u₁, acos(1 - 2u₂) / T(π), u₃)
  end
  (randpoint() for _ in 1:(method.size))
end
