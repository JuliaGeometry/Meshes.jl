# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

  # visualize as built-in image with or without interpolation
  Makie.map!(plot, [:object, :colorant], [:x, :y, :C, :interpolate]) do grid, colorant
    sz = size(grid)
    nv = nvertices(grid)
    nc = colorant isa AbstractVector ? length(colorant) : 1

    x, y = map(Meshes.xyz(grid)) do c
      extrema(ustrip.(c))
    end

    C = if nc == 1
      fill(colorant, sz)
    elseif nc == nv
      reshape(colorant, sz .+ 1)
    else
      reshape(colorant, sz)
    end

    x, y, C, (nc == nv)
  end
  Makie.image!(plot, plot.x, plot.y, plot.C, interpolate=plot.interpolate)

  if plot.showsegments[]
    vizfacets!(plot)
  end
end

function vizgrid!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{3}, ::Val{3})
  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

  # number of vertices and colors
  Makie.map!(plot, [:object, :colorant], [:nv, :nc]) do grid, colorant
    nv = nvertices(grid)
    nc = colorant isa AbstractVector ? length(colorant) : 1
    nv, nc
  end

  if plot.nc[] == plot.nv[]
    # visualize as a quadrangle mesh so that
    # colors can be specified at vertices
    Makie.map!(plot, :object, :mesh) do grid
      verts = map(asmakie, eachvertex(grid))
      quads = reduce(
        vcat,
        map(elements(topology(grid))) do elem
          i1, i2, i3, i4, i5, i6, i7, i8 = indices(elem)
          [
            # follow the same order of vertices specified
            # in the `boundary` method for `Hexahedron`
            GB.QuadFace(i4, i3, i2, i1),
            GB.QuadFace(i6, i5, i1, i2),
            GB.QuadFace(i3, i7, i6, i2),
            GB.QuadFace(i4, i8, i7, i3),
            GB.QuadFace(i1, i5, i8, i4),
            GB.QuadFace(i6, i7, i8, i5)
          ]
        end
      )
      GB.Mesh(verts, quads)
    end
    Makie.mesh!(plot, plot.mesh, color=plot.colorant, shading=true)
  else
    # visualize as built-in meshscatter
    Makie.map!(plot, :object, [:xyz, :rec]) do grid
      sp = ustrip.(spacing(grid))
      xs, ys, zs = map(Meshes.xyz(grid)) do c
        ustrip.(c[(begin + 1):end])
      end
      xyz = [(x, y, z) for z in zs for y in ys for x in xs]
      rec = Makie.Rect3(-1 .* sp, sp)
      xyz, rec
    end
    Makie.meshscatter!(plot, plot.xyz, marker=plot.rec, markersize=1, color=plot.colorant)
  end

  if plot.showsegments[]
    vizfacets!(plot)
  end
end

function vizgridfacets!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  Makie.map!(plot, :object, [:xfacets, :yfacets]) do grid
    x, y = map(c -> ustrip.(c), Meshes.xyz(grid))
    xysegments(x, y)
  end
  Makie.lines!(plot, plot.xfacets, plot.yfacets, color=plot.segmentcolor, linewidth=plot.segmentsize)
end

function vizgridfacets!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{3}, ::Val{3})
  Makie.map!(plot, :object, [:xfacets, :yfacets, :zfacets]) do grid
    x, y, z = map(c -> ustrip.(c), Meshes.xyz(grid))
    xyzsegments(x, y, z)
  end
  Makie.lines!(plot, plot.xfacets, plot.yfacets, plot.zfacets, color=plot.segmentcolor, linewidth=plot.segmentsize)
end
