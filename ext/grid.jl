# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{Grid}})
  grid = plot[:object][]
  pdim = paramdim(grid)
  edim = embeddim(grid)
  vizgrid!(plot, Val(pdim), Val(edim))
end

vizgrid!(plot, ::Val{2}, ::Val{2}) = vizmesh!(plot, Val(2), Val(2))

function vizgrid!(plot, ::Val{3}, ::Val{3})
  grid = plot[:object]
  color = plot[:color]

  # number of vertices and colors
  nv = Makie.@lift nvertices($grid)
  nc = Makie.@lift $color isa AbstractVector ? length($color) : 1

  if nv[] == nc[]
    error("not implemented")
  else
    vizmesh!(plot, Val(3), Val(3))
  end
end

# ----------------
# SPECIALIZATIONS
# ----------------

include("grid/cartesian.jl")
include("grid/rectilinear.jl")
include("grid/transformed.jl")
include("grid/structured.jl")

# helper functions to create a minimum number
# of line segments within Cartesian/Rectilinear grid
function xysegments(xs, ys)
  coords = []
  for x in xs
    push!(coords, (x, first(ys)))
    push!(coords, (x, last(ys)))
    push!(coords, (NaN, NaN))
  end
  for y in ys
    push!(coords, (first(xs), y))
    push!(coords, (last(xs), y))
    push!(coords, (NaN, NaN))
  end
  x = getindex.(coords, 1)
  y = getindex.(coords, 2)
  x, y
end

function xyzsegments(xs, ys, zs)
  coords = []
  for y in ys, z in zs
    push!(coords, (first(xs), y, z))
    push!(coords, (last(xs), y, z))
    push!(coords, (NaN, NaN, NaN))
  end
  for x in xs, z in zs
    push!(coords, (x, first(ys), z))
    push!(coords, (x, last(ys), z))
    push!(coords, (NaN, NaN, NaN))
  end
  for x in xs, y in ys
    push!(coords, (x, y, first(zs)))
    push!(coords, (x, y, last(zs)))
    push!(coords, (NaN, NaN, NaN))
  end
  x = getindex.(coords, 1)
  y = getindex.(coords, 2)
  z = getindex.(coords, 3)
  x, y, z
end
