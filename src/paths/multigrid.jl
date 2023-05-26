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
  sizes = size(grid)
  steps = sizes .- 1
  nelems = prod(sizes)
  linear = LinearIndices(sizes)
  while length(path) < nelems
    ranges = ntuple(d -> 1:steps[d]:sizes[d], Dim)
    cinds = CartesianIndices(ranges)
    linds = [linear[cind] for cind in cinds]
    append!(path, setdiff(linds, path))
    steps = ceil.(Int, steps ./ 2)
  end
  path
end
