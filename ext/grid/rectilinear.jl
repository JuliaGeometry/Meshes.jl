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
      # visualize as simple mesh
      vizmesh!(plot)
    else
      # visualize as built-in heatmap
      Makie.map!(plot, [:object, :colorant], [:x, :y, :C, :usetiles]) do grid, colorant
        sz = size(grid)
        ne = nelements(grid)
        x, y = map(c -> ustrip.(c), Meshes.xyz(grid))
        C = colorant isa AbstractVector ? reshape(colorant, sz) : fill(colorant, sz)
        x, y, C, (ne ≥ 100000)
      end

      if plot.usetiles[]
        for (i, tile) in enumerate(TileIterator(axes(plot.C[]), (100, 100)))
          xi = Symbol(:x, i)
          yi = Symbol(:y, i)
          Ci = Symbol(:C, i)
          Makie.map!(plot, [:x, :y, :C], [xi, yi, Ci]) do x, y, C
            xlims, ylims = map((x, y), tile) do c, t
              # note: non-interpolate case requires adding 1 to the end of the tile range
              ctile = first(t):(last(t) + 1)
              extrema(view(c, ctile))
            end
            xlims, ylims, view(C, tile...)
          end
          Makie.heatmap!(plot, plot[xi], plot[yi], plot[Ci])
        end
      else
        Makie.heatmap!(plot, plot.x, plot.y, plot.C)
      end
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
