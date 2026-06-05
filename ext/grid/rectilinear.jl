# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{RectilinearGrid}}, M::Type{<:𝔼}, pdim::Val{2}, edim::Val{2})
  if crs(plot.object[]) <: Cartesian
    # number of vertices and colors
    Makie.map!(plot, [:object, :colorant], [:nv, :nc]) do grid, colorant
      nv = nvertices(grid)
      nc = colorant isa AbstractVector ? length(colorant) : 1
      nv, nc
    end

    if plot.nc[] == plot.nv[]
      # visualize as a simple mesh so that
      # colors can be specified at vertices
      vizmesh!(plot)
    else
      # visualize as built-in heatmap
      Makie.map!(plot, [:object, :colorant], [:x, :y, :C]) do grid, colorant
        sz = size(grid)
        x, y = map(c -> ustrip.(c), Meshes.xyz(grid))
        C = if plot.nc[] == 1
          fill(colorant, sz)
        else
          reshape(colorant, sz)
        end
        x, y, C
      end
      Makie.heatmap!(plot, plot.x, plot.y, plot.C)
    end

    if plot.showsegments[]
      vizfacets!(plot)
    end
  else
    vizgridfallback!(plot, M, pdim, edim)
  end
end

function vizgridfacets!(plot::Viz{<:Tuple{RectilinearGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  Makie.map!(plot, :object, [:xfacets, :yfacets]) do grid
    x, y = map(c -> ustrip.(c), Meshes.xyz(grid))
    xysegments(x, y)
  end
  Makie.lines!(plot, plot.xfacets, plot.yfacets, color=plot.segmentcolor, linewidth=plot.segmentsize)
end
