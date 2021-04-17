# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FilteredSearch(method, nmax=Inf; maxpercategory, maxpersector, metric=Euclidean())

A method for searching at most `nmax` neighbors using `method`. Extra
constraints available: `maxpercategory` and `maxpersector`. The neighbors are
sorted using `metric` distance (priority to the nearest). To preserve the order
of the initial neighbors returned by `method`, set `metric = nothing`.

## Max per category

It can be a `NamedTuple` or `Dict`, where the first element is the property name
and the second defines the max neighbor per available categories of the property.
In the example below, no more than 2 neighbors of same geometry are selected.

`FilteredSearch(method, maxpercategory = (geometry = 2,))`
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

## Max per sector

Only allow a max number of neighbors inside each sector. The sectors are the
cartesian quadrants (2-D) or octants (3-D). If the neighborhood is an `Ellipsoid`,
the sectors will match the rotated quadrants/octants.

`FilteredSearch(method, maxpersector = 2)`
```ascii
     _________
  //     . □   \\           ○  Search point
 //      . ▩    \\        ▩ ▲  Neighbors selected
// ▲ ▲   . ▩     \\       □ △  Neighbors ignored
 ........○........
\\ ▲ ▲   .       //
 \\      .      //
  \\ △   . ▩   //
     ‾‾‾‾‾‾‾‾‾
```
"""
struct FilteredSearch{M<:NeighborSearchMethod} <: BoundedNeighborSearchMethod
  method::M
  nmax
  maxpercategory
  maxpersector
  metric
end

FilteredSearch(method::M, nmax=Inf; maxpercategory=nothing,
               maxpersector=nothing, metric=Euclidean()) where {M} =
  FilteredSearch{M}(method, nmax, maxpercategory, maxpersector, metric)

maxneighbors(method::FilteredSearch) = method.nmax

function search!(neighbors, pₒ::Point, method::FilteredSearch; mask=nothing)
  meth = method.method
  obj  = meth.domain
  nmax = method.nmax

  # initial search, using mask if necessary
  inds = search(pₒ, meth, mask=mask)

  # `inds` is not guranteed to be sorted by distance. also, a different
  # metric might be preferred to give priority to the FilteredSearch.
  # unless `metric = nothing`, the indices are sorted here using metric
  if !isnothing(method.metric)
    dists  = [evaluate(method.metric, coordinates(pₒ), coordinates(centroid(obj, ind))) for ind in inds]
    sorted = sortperm(dists)
    inds   = inds[sorted]
  end

  # initialize category and sectors constraints if necessary. `Dict` containers
  # are here created with the information necessary for further filtering.
  categs  = initcategories(obj, inds, method.maxpercategory)
  sectors = initsectors(meth, method.maxpersector)

  # loop each neighbor candidate
  nneigh = 0
  for i in inds
    # check category. if exceeded the category, the neighbor is ignored;
    # otherwise categs[:count][col][cat[col]] is incremented in the end after
    # nmax and maxpersector checks
    if categs[:use]
      cat, pass = Dict(), false
      tab = Tables.columns(values(obj))
      for col in keys(categs[:max])
        cat[col] = Tables.getcolumn(tab, col)[i]
        categs[:count][col][cat[col]] >= categs[:max][col] && (pass = true)
      end
      pass && continue
    end

    # check sectors; if there are no space for neighbors in the indsector, it's
    # ignored; otherwise sectors[:count][indsector] is incremented in the end
    # after nmax checks
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

    nneigh == nmax && break

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

# initialize categories constraints if necessary. e.g. if `maxpercategory =
# (holeid = 2,)`, dict[:count][:holeid] will have all the possible values of
# holeid, initiliazed to zero for further counting.
function initcategories(obj, inds, catgs)
  catgs == nothing && return Dict(:use => false)
  tab   = Tables.columns(values(obj))
  catvals = Dict(k => unique(Tables.getcolumn(tab, k)[inds]) for k in keys(catgs))
  counter = Dict(k => Dict(zip(v, zeros(Int,length(v)))) for (k, v) in catvals)
  Dict(:use => true, :max => catgs, :count => counter)
end

# initialize sectors constraints if necessary. this will create a counter for
# each sector and a rotation matrix to reverse ellipsoid rotation (if an
# ellipsoid search is used). a rotation noise is also added to simplify filtering.
# it is used here because if data is regular and fall onto axes (e.g. cartesian
# grid data), the sector filtering will return misleading results without it.
function initsectors(method, maxpersector)
  maxpersector == nothing && return Dict(:use => false)
  N = embeddim(method.domain)
  ellips = method isa NeighborhoodSearch && method.neigh isa Ellipsoid

  # add a noise rotation to avoid samples at sectors limits
  rotangs = N == 2 ? [0.001] : [0.001, 0.001, 0.001]
  convention = ellips ? method.neigh.convention : Datamine

  # reverse rotation if it is an ellipsoid
  if ellips
    extrinsic = isextrinsic(method.neigh.convention)
    invertfac = extrinsic ? -1 : 1
    rotangs   = invertfac * (method.neigh.angles .+ rotangs)
  end

  # rotation matrix
  rotmatx = rotmat(rotangs, convention)'

  Dict(:use => true, :count => zeros(Int, 2^N), :rotmat => rotmatx, :max => maxpersector)
end

# this function returns at which quadrant/octant the given centered coordinate
# is. instead of making multiple if statements, a container is used to get the
# combination of positive/negative coordinates and assign a code to it.
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
