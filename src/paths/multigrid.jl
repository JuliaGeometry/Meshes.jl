# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MultiGridPath()

to-do
"""
struct MultiGridPath <: Path end

function traverse(grid::Grid{Dim}, ::MultiGridPath) where {Dim}
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
