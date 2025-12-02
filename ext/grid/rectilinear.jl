# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{RectilinearGrid}}, M::Type{<:𝔼}, pdim::Val{2}, edim::Val{2})
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  showsegments = plot[:showsegments]

  if crs(grid[]) <: Cartesian
    # process color spec into colorant
    colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

    # number of vertices and colors
    nv = Makie.@lift nvertices($grid)
    nc = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

    # grid coordinates
    xyz = Makie.@lift map(x -> ustrip.(x), Meshes.xyz($grid))
    xs = Makie.@lift $xyz[1]
    ys = Makie.@lift $xyz[2]

    if nc[] == 1
      # visualize bounding box with single color for maximum performance
      # perform simplexification to avoid infinite loops (boxes <-> grids)
      bbox = Makie.@lift boundingbox($grid)
      tmesh = Makie.@lift simplexify($bbox)
      viz!(plot, tmesh, color=colorant)
    else
      if nc[] == nv[]
        # visualize as a simple mesh so that
        # colors can be specified at vertices
        vizmesh!(plot)
      else
        # visualize as built-in heatmap
        sz = Makie.@lift size($grid)
        C = Makie.@lift reshape($colorant, $sz)
        Makie.heatmap!(plot, xs, ys, C)
      end
    end

    if showsegments[]
      vizfacets!(plot)
    end
  else
    vizgridfallback!(plot, M, pdim, edim)
  end
end

function vizgridfacets!(plot::Viz{<:Tuple{RectilinearGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  grid = plot[:object]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  xyz = Makie.@lift map(x -> ustrip.(x), Meshes.xyz($grid))
  tup = Makie.@lift xysegments($xyz...)
  x, y = Makie.@lift($tup[1]), Makie.@lift($tup[2])
  Makie.lines!(plot, x, y, color=segmentcolor, linewidth=segmentsize)
end
