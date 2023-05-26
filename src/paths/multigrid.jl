# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MultiGridPath()

to-do
"""
struct MultiGridPath <: Path end

function traverse(grid::Grid{Dim}, ::MultiGridPath) where {Dim}
  path = Int[]
  dims = size(grid)
  steps = dims .- 1
  nelems = prod(dims)
  linear = LinearIndices(dims)
  while length(path) < nelems
    ranges = ntuple(d -> 1:steps[d]:dims[d], Dim)
    cinds = CartesianIndices(ranges)
    linds = [linear[cind] for cind in cinds]
    append!(path, setdiff(linds, path))
    steps = ceil.(Int, steps ./ 2)
  end
  path
end
