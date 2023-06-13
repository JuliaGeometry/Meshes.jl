# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RandomPath()

Traverse a domain with `N` elements in a random
permutation of `1:N`.
"""
struct RandomPath{R<:AbstractRNG} <: Path
  rng::R
end

RandomPath() = RandomPath(Random.GLOBAL_RNG)

traverse(domain, path::RandomPath) = randperm(path.rng, nelements(domain))
