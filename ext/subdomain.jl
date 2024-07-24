# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{SubDomain}})
  subdom = plot[:object]
  M = Makie.@lift manifold($subdom)
  pdim = Makie.@lift paramdim($subdom)
  edim = Makie.@lift embeddim($subdom)
  vizsubdom!(plot, M[], Val(pdim[]), Val(edim[]))
end

function vizsubdom!(plot, ::Type{<:🌐}, pdim::Val, edim::Val)
  @warn "geodesic geometries can't be visualized yet. Visualizing as Euclidean..."
  vizsubdom!(plot, 𝔼, pdim, edim)
end

function vizsubdom!(plot, ::Type{<:𝔼}, ::Val, ::Val)
  subdom = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # collect geometries
  geoms = Makie.@lift collect($subdom)

  vizmany!(plot, geoms, colorant)
end

const SubCartesianGrid{C<:CRS,Mₚ<:AbstractManifold,Dim} = SubDomain{𝔼{Dim},C,<:CartesianGrid{C,Mₚ,Dim}}

function vizsubdom!(plot::Viz{<:Tuple{SubCartesianGrid}}, ::Type{<:𝔼}, ::Val, ::Val)
  subgrid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # retrieve grid paramaters
  gparams = Makie.@lift let
    grid = parent($subgrid)
    dim = embeddim(grid)
    sp = ustrip.(spacing(grid))

    # coordinates of centroids
    coord(e) = ustrip.(to(centroid(e)))
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
