# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{RectilinearGrid}}, M::Type{<:𝔼}, pdim::Val{2}, edim::Val{2})
  showsegments = plot[:showsegments]

  if crs(plot[:object][]) <: Cartesian
    # process color spec into colorant
    Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

    # number of vertices and colors
    Makie.map!(nvertices, plot, [:object], :nv)
    Makie.map!(plot, [:colorant], :nc) do colorant
      colorant isa AbstractVector ? length(colorant) : 1
    end

    # grid coordinates
    Makie.map!(plot, [:object], :xyz) do grid
      map(x -> ustrip.(x), Meshes.xyz(grid))
    end
    Makie.map!(plot, [:xyz], :xs) do xyz
      xyz[1]
    end
    Makie.map!(plot, [:xyz], :ys) do xyz
      xyz[2]
    end

    if plot[:nc][] == plot[:nv][]
      # visualize as a simple mesh so that
      # colors can be specified at vertices
      vizmesh!(plot)
    else
      # visualize as built-in heatmap
      Makie.map!(size, plot, [:object], :sz)
      Makie.map!(plot, [:nc, :colorant, :sz], :C) do nc, colorant, sz
        nc == 1 ? fill(colorant, sz) : reshape(colorant, sz)
      end
      Makie.heatmap!(plot, plot[:xs], plot[:ys], plot[:C])
    end

    if showsegments[]
      vizfacets!(plot)
    end
  else
    vizgridfallback!(plot, M, pdim, edim)
  end
end

function vizgridfacets!(plot::Viz{<:Tuple{RectilinearGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  Makie.map!(plot, [:object], :facets_xy) do grid
    x, y = Meshes.xyz(grid)
    xysegments(ustrip.(x), ustrip.(y))
  end
  Makie.map!(plot, [:facets_xy], :facets_x) do xy
    xy[1]
  end
  Makie.map!(plot, [:facets_xy], :facets_y) do xy
    xy[2]
  end

  Makie.lines!(plot, plot[:facets_x], plot[:facets_y], color=segmentcolor, linewidth=segmentsize)
end
