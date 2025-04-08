# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RandomPath()

Traverse a domain with `n` elements in a random
permutation of `1:n`.
"""
struct RandomPath{R<:AbstractRNG} <: Path
  rng::R
end

RandomPath() = RandomPath(Random.default_rng())

traverse(domain, path::RandomPath) = randperm(path.rng, nelements(domain))
