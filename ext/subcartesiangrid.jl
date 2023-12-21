# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const SubCartesianGrid{Dim,T} = Meshes.SubDomain{Dim,T,<:CartesianGrid{Dim,T}}

function Makie.plot!(plot::Viz{<:Tuple{SubCartesianGrid}})
  subgrid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # retrieve grid paramaters
  gparams = Makie.@lift let
    grid = parent($subgrid)
    dim = embeddim(grid)
    sp = spacing(grid)

    # coordinates of centroids
    coord(e) = coordinates(centroid(e))
    coords = [coord(e) .+ sp ./ 2 for e in $subgrid]

    # rectangle marker
    marker = Makie.Rect{dim}(-1 .* sp, sp)

    # enable shading in 3D
    shading = dim == 3 ? Makie.FastShading : Makie.NoShading

    coords, marker, shading
  end

  # unpack observable parameters
  coords = Makie.@lift $gparams[1]
  marker = Makie.@lift $gparams[2]
  shading = Makie.@lift $gparams[3]

  # all geometries are equal, use mesh scatter
  Makie.meshscatter!(plot, coords, marker=marker, markersize=1, color=colorant, shading=shading)
end
