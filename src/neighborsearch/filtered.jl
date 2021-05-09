# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FilteredSearch(data, method, maxneighbors; maxpervar=Dict())

A search method for filtering the results of an existing search
`method` that was applied to the `data`. The results are filtered
according to a maximum number of neighbors `maxneighbors`.

Optionally, a dictionary `maxpervar` can be specified to filter
a maximum number per value of categorical variables. The keys
of the dictionary are names of categorical variables and the
values define the maximum number of neighbors per value of the
variable. In the example below, no more than 2 neighbors of same
`rocktype` are selected:

```julia
FilteredSearch(data, method, maxneighbors, maxpervar = (rocktype = 2,))
```

```ascii
     _________
  //     . □   \\           ○  Search point
 //      . ▩    \\        ▩ ▲  Neighbors selected
// △ ▲   . ▩     \\       □ △  Neighbors ignored
 ........○........
\\ △ ▲   .       //
 \\      .      //
  \\ △   . □   //
     ‾‾‾‾‾‾‾‾‾
```
"""
struct FilteredSearch{D,M<:NeighborSearchMethod} <: BoundedNeighborSearchMethod
  data::D
  method::M
  maxneighbors::Int
  maxpervar::Dict{Symbol,Int}
end

FilteredSearch(data::D, method::M, maxneighbors; maxpervar=Dict()) where {D,M} =
  FilteredSearch{D,M}(data, method, maxneighbors, maxpervar)

maxneighbors(method::FilteredSearch) = method.maxneighbors

function search!(neighbors, pₒ::Point, method::FilteredSearch; mask=nothing)
  # retrieve table of values
  table = values(method.data)

  # initialize count per variable value
  maxpervar = method.maxpervar
  count4var = Dict(
    map(keys(maxpervar)) do
      col  = Tables.getcolumn(table, var)
      var => Dict(val => 0 for val in unique(col))
    end
  )

  # initial search, using mask if necessary
  inds = search(pₒ, method.method, mask=mask)

  # maximum per variable
  retained = Int[]
  rowtable = Tables.rowtable(table)
  for ind in inds
    retain = true
    for (var, max) in maxpervar
      for (val, count) in count4var[var]
        if count ≥ max
          retain = false
        end
      end
    end
    if retain
      push!(retained, ind)
      # update counts
      row = rowtable[ind]
      for var in keys(count4var)
        val = row[var]
        count4var[var][val] += 1
      end
    end
  end

  # maximum neighbors
  nmax   = method.maxneighbors
  nneigh = min(nmax, length(retained))
  @inbounds for i in 1:nneigh
    neighbors[i] = retained[i]
  end

  nneigh
end
