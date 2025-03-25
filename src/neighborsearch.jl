# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NeighborSearchMethod

A method for searching neighbors given a reference point.
"""
abstract type NeighborSearchMethod end

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

Return the maximum number of neighbors obtained with bounded search `method`.

See [`BoundedNeighborSearchMethod`](@ref) for additional details.
"""
function maxneighbors end

"""
    search!(neighbors, pₒ, method; mask=nothing)

Update `neighbors` of point `pₒ` using bounded search `method` and return
number of neighbors found. Optionally, specify a `mask` for all indices of
the domain.

See [`BoundedNeighborSearchMethod`](@ref) for additional details.
"""
function search! end

"""
    searchdists!(neighbors, distances, pₒ, method; mask=nothing)

Update `neighbors` and `distances` of point `pₒ` using bounded search `method`
and return number of neighbors found. Optionally, specify a `mask` for all
indices of the domain.

See [`BoundedNeighborSearchMethod`](@ref) for additional details.
"""
function searchdists! end

"""
    searchdists(pₒ, method, mask=nothing)

Return neighbors and distances of point `pₒ` using bounded search `method`.
Optionally, specify a `mask` for all indices of the domain.

See [`BoundedNeighborSearchMethod`](@ref) for additional details.
"""
function searchdists end

# ----------
# FALLBACKS
# ----------

function search!(neighbors, pₒ::Point, method::BoundedNeighborSearchMethod; mask=nothing)
  distances = Vector{lentype(pₒ)}(undef, maxneighbors(method))
  searchdists!(neighbors, distances, pₒ, method; mask)
end

function search(pₒ::Point, method::BoundedNeighborSearchMethod; mask=nothing)
  neighbors = Vector{Int}(undef, maxneighbors(method))
  nneigh = search!(neighbors, pₒ, method; mask=mask)
  view(neighbors, 1:nneigh)
end

function searchdists(pₒ::Point, method::BoundedNeighborSearchMethod; mask=nothing)
  neighbors = Vector{Int}(undef, maxneighbors(method))
  distances = Vector{lentype(pₒ)}(undef, maxneighbors(method))
  nneigh = searchdists!(neighbors, distances, pₒ, method; mask)
  view(neighbors, 1:nneigh), view(distances, 1:nneigh)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("neighborsearch/ball.jl")
include("neighborsearch/knearest.jl")
include("neighborsearch/kball.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

# NearestNeighbors.jl only accepts AbstractVector
_svector(c::CRS) = SVector(CoordRefSystems.raw(c))
