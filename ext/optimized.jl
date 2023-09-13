# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ----------------------------------
# recipes optimized for performance
# ----------------------------------

const SubGrid{Dim,T} = Meshes.SubDomain{Dim,T,<:CartesianGrid{Dim,T}}

Makie.plottype(::SubGrid) = Viz{<:Tuple{SubGrid}}

function Makie.plot!(plot::Viz{<:Tuple{SubGrid}})
  gridview = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # retrieve grid paramaters
  gparams = Makie.@lift let
    grid, _ = unview($gridview)
    dims = embeddim(grid)
    sp = spacing(grid)

    # coordinates of centroids
    coord(e) = coordinates(centroid(e))
    coords = [coord(e) .+ sp ./ 2 for e in $gridview]

    # rectangle marker
    marker = Makie.Rect{dims}(-1 .* sp, sp)

    # enable shading in 3D
    shading = dims == 3

    coords, marker, shading
  end

  # unpack observable parameters
  coords = Makie.@lift $gparams[1]
  marker = Makie.@lift $gparams[2]
  shading = Makie.@lift $gparams[3]

  # all geometries are equal, use mesh scatter
  Makie.meshscatter!(plot, coords, marker=marker, markersize=1, color=colorant, shading=shading)
end
