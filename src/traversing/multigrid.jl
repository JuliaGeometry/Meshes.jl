# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MultiGridPath()

Traverse a grid using a nested grid system where the path starts at
the coarsest scale and moves to progressively finer scales.
"""
struct MultiGridPath <: Path end

function traverse(grid::Grid, ::MultiGridPath)
  Dim = embeddim(grid)
  dims = size(grid)
  nelems = prod(dims)
  linear = LinearIndices(dims)

  path = Int[]
  steps = dims .- 1
  while length(path) < nelems
    ranges = ntuple(d -> 1:steps[d]:dims[d], Dim)

    for cind in CartesianIndices(ranges)
      lind = linear[cind]
      if lind ∉ path
        push!(path, lind)
      end
    end

    steps = ceil.(Int, steps ./ 2)
  end

  path
end

function traverse(domain::SubGrid, path::MultiGridPath)
  inds = traverse(parent(domain), path)
  pinds = parentindices(domain)
  filter!(∈(pinds), inds)
  inds
end
