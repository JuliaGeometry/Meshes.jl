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
  nelems = nelements(grid)
  linear = LinearIndices(sizes)
  while length(path) < nelems
    cinds = CartesianIndices(ntuple(d -> 1:steps[d]:sizes[d], Dim))
    inds = [linear[cind] for cind in cinds]
    append!(path, setdiff(inds, path))
    steps = ceil.(Int, steps ./ 2)
  end
  path
end
