# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ComponentsPartition()

A method for partitioning meshes into connected components.
"""
struct ComponentsPartition <: PartitionMethod end

function partitioninds(::AbstractRNG, domain::Domain, ::ComponentsPartition)
  assertion(domain isa Mesh, "connected components only defined for meshes")
  elems = collect(elements(topology(domain)))
  components(elems), Dict()
end

function components(elems::AbstractVector{<:Connectivity})
  # initialize list of connected components
  comps = [[firstindex(elems)]]

  # initialize list of seen vertices
  seen = Set{Int}()
  for v in indices(first(elems))
    push!(seen, v)
  end

  # remaining elements to process
  remaining = collect(eachindex(elems)[2:end])

  added = false
  while !isempty(remaining)
    iter = 1
    while iter ≤ length(remaining)
      elem = elems[remaining[iter]]

      # manually union-split most common polytopes
      # for type stability and maximum performance
      isadjacent = if elem isa Connectivity{Triangle,3}
        adjelem!(seen, elem)
      elseif elem isa Connectivity{Quadrangle,4}
        adjelem!(seen, elem)
      else
        adjelem!(seen, elem)
      end

      if isadjacent
        push!(last(comps), popat!(remaining, iter))
        added = true
      else
        iter += 1
      end
    end

    if added
      # new vertices were "seen" while iterating `remaining`, so
      # we need to iterate again because there may be elements
      # which are now adjacent with the newly "seen" vertices
      added = false
    elseif !isempty(remaining)
      # there are more elements, but none are adjacent to
      # previously seen elements; pop a new element from
      # the original list to start a new connected component
      push!(comps, Int[])
      push!(last(comps), popfirst!(remaining))

      # a disconnected component means that ≥n-1 vertices in
      # the newest element haven't been "seen"; it's possible
      # the new component is connected by a single vertex
      for v in indices(elems[last(last(comps))])
        push!(seen, v)
      end
    end
  end

  comps
end
