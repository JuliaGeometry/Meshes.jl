# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BoundedSearch(method, nmax)

A method for searching at most `nmax` neighbors using `method`.
"""
struct BoundedSearch{M<:NeighborSearchMethod} <: BoundedNeighborSearchMethod
  method::M
  nmax::Int
end

maxneighbors(method::BoundedSearch) = method.nmax

function search!(neighbors, pₒ::Point, method::BoundedSearch; mask=nothing)
  inds = search(pₒ, method.method)
  nmax = method.nmax

  if isnothing(mask)
    nneigh = min(length(inds), nmax)
    @inbounds for i in 1:nneigh
      neighbors[i] = inds[i]
    end
  else
    nneigh = 0
    @inbounds for ind in inds
      if mask[ind]
        nneigh += 1
        neighbors[nneigh] = ind
      end
      nneigh == nmax && break
    end
  end

  return nneigh
end
