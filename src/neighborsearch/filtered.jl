# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FilteredSearch(method, nmax=Inf; maxpercategory, maxpersector, order=Euclidean())

A method for searching at most `nmax` neighbors using `method`. Extra
constraints available: `maxpercategory` and `maxpersector`. The priority is
given to the nearest neighbor using `order`.

## Max per category

It can be a `NamedTuple` or `Dict`, where the first element is the property name
and the second defines the max neighbor per available categories of the property.
In the example below, no more than 2 neighbors from the same hole are selected.

`FilteredSearch(method, maxpercategory = (holeid = 2,))`
[ASCII illustration]

## Max per sector

Only allow a max number of neighbors inside each sector. The sectors are the
cartesian quadrants (2-D) or octants (3-D). If the neighborhood is an `Ellipsoid`,
the sectors will match the rotated quadrants/octants.

`FilteredSearch(method, maxpersector = 2)`
[ASCII illustration]
"""
struct FilteredSearch{M<:NeighborSearchMethod} <: BoundedNeighborSearchMethod
  method::M
  nmax::Int
  maxpercategory
  maxpersector
  order
end

FilteredSearch(method::M, nmax=0; maxpercategory=nothing,
               maxpersector=nothing, order=Euclidean()) where {M} =
  FilteredSearch{M}(method, nmax, maxpercategory, maxpersector, order)

maxneighbors(method::FilteredSearch) = method.nmax

function search!(neighbors, pₒ::Point, method::FilteredSearch; mask=nothing)
  meth = method.method
  obj  = meth.domain
  inds = search(pₒ, meth, mask=mask)
  nmax = method.nmax == 0 ? Inf : method.nmax

  # get distances and give priority to closest neighbors if necessary
  if !(meth isa KNearestSearch)
    indcoords = reduce(hcat, coordinates.([centroid(obj, s) for s in inds]))
    dists     = colwise(method.order, coordinates(pₒ), indcoords)
    sorted    = sortperm(dists)
    inds      = inds[sorted]
  end

  # initialize category and sectors constraints if necessary
  categs  = initcategories(obj, inds, method.maxpercategory)
  sectors = initsectors(meth, method.maxpersector)

  # check about initialization size
  nneigh = 0
  for i in inds
    nneigh == nmax && break

    # check category
    if categs[:use]
      cat, pass = Dict(), false
      for col in keys(categs[:max])
        cat[col] = getproperty(values(obj)[i],col)
        categs[:count][col][cat[col]] >= categs[:max][col] && (pass = true)
      end
      pass && continue
    end

    # check sectors
    if sectors[:use]
      centered  = sectors[:rotmat]' * (centroid(obj,i) - pₒ)
      indsector = getsector(centered)
      sectors[:count][indsector] >= sectors[:max] && continue
    end

    # add neighbor
    nneigh += 1
    if nmax == Inf
      push!(neighbors, i)
    else
      neighbors[nneigh] = i
    end

    # add counters
    sectors[:use] && (sectors[:count][indsector] += 1)
    if categs[:use]
      for col in keys(categs[:max])
        categs[:count][col][cat[col]] += 1
      end
    end

  end

  # slice neighbors if nmax was not reached
  length(neighbors) > nneigh && (neighbors = neighbors[1:nneigh])

  nneigh
end

# initialize categories constraints if necessary
function initcategories(obj, inds, catgs)
  catgs == nothing && return Dict(:use => false)
  tab   = values(obj)
  catvals = Dict(k => unique(getproperty(tab[inds], k)) for k in keys(catgs))
  counter = Dict(k => Dict(zip(v, zeros(Int,length(v)))) for (k, v) in catvals)
  Dict(:use => true, :max => catgs, :count => counter)
end

function initsectors(method, maxpersector)
  # initialize sectors constraints if necessary
  maxpersector == nothing && return Dict(:use => false)
  N = embeddim(method.domain)
  ellips = method isa NeighborhoodSearch && method.neigh isa Ellipsoid

  #add a noise to avoid samples at sectors limits
  rotmatx = N == 2 ? rotmat([0.01], GSLIB) : rotmat([0.01, 0.01, 0.01], GSLIB)

  # add a reverse rotation matrix if it is an ellipsoid
  if ellips
    extrinsic = isextrinsic(method.neigh.convention)
    invertfac = extrinsic ? -1 : 1
    rotmatx  += rotmat(invertfac*method.neigh.angles, method.neigh.convention)'
  end
  Dict(:use => true, :count => zeros(Int, 2^N), :rotmat => rotmatx, :max => maxpersector)
end

function getsector(coords::AbstractVector)
  # get dimensions
  N    = size(coords,1)
  dims = ntuple(i->2, N)

  # get which coord is negative. assign an id to each sector
  signcoord = (coords .< 0) .+ 1
  sectid    = reshape(1:2^N, dims)

  # return sector id to combination of coordinates signs
  sectid[signcoord...]
end
