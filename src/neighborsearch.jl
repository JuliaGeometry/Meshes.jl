# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NeighborSearchMethod

A method for searching neighbors given a reference point.
"""
abstract type NeighborSearchMethod end

"""
    search!(neighbors, pₒ, method; mask=nothing)

Update `neighbors` of point `pₒ` using `method` and return
number of neighbors found. Optionally, specify a `mask` for
all indices of the domain.
"""
function search! end

"""
    search(pₒ, method, mask=nothing)

Return neighbors of point `pₒ` using `method`. Optionally,
specify a `mask` for all indices of the domain.
"""
function search end

"""
    BoundedNeighborSearchMethod

A method for searching neighbors with the property that the number of neighbors
is bounded above by a known constant (e.g. k-nearest neighbors).
"""
abstract type BoundedNeighborSearchMethod <: NeighborSearchMethod end

"""
    maxneighbors(method)

Return the maximum number of neighbors obtained with `method`.
"""
function maxneighbors end

# ----------
# FALLBACKS
# ----------

function search(pₒ::Point, method::BoundedNeighborSearchMethod; mask=nothing)
  neighbors = Vector{Int}(undef, maxneighbors(method))
  nneigh = search!(neighbors, pₒ, method; mask=mask)
  view(neighbors, 1:nneigh)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("neighborsearch/ball.jl")
include("neighborsearch/knearest.jl")
include("neighborsearch/kball.jl")
include("neighborsearch/bounded.jl")
