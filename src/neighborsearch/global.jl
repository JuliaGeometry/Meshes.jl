# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GlobalSearch(domain)

A method for searching all elements of the `domain`.
"""
struct GlobalSearch{D} <: BoundedNeighborSearchMethod
  domain::D
end

maxneighbors(method::GlobalSearch) = nelements(method.domain)

function search!(neighbors, pₒ::Point, method::GlobalSearch; mask=nothing)
  nelem = nelements(method.domain)

  if isnothing(mask)
    nneigh = nelem
    @inbounds for ind in 1:nelem
      neighbors[ind] = ind
    end
  else
    nneigh = 0
    nmax = sum(mask)
    @inbounds for ind in 1:nelem
      if mask[ind]
        nneigh += 1
        neighbors[nneigh] = ind
      end
      nneigh == nmax && break
    end
  end

  nneigh
end
