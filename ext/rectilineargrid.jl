# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{RectilinearGrid}})
  # retrieve parameters
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  segmentsize = plot[:segmentsize]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # number of dimensions, vertices and colors
  nd = Makie.@lift embeddim($grid)
  nv = Makie.@lift nvertices($grid)
  nc = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

  # size, coordinates and centroid coordinates
  sz = Makie.@lift size($grid)
  xyz = Makie.@lift Meshes.xyz($grid)

  # dispatch different recipes
  if nc[] == 1
    # visualize bounding box with a single
    # color for maximum performance
    bbox = Makie.@lift boundingbox($grid)
    viz!(plot, bbox, color=colorant)
  else
    # search other visualization methods
    if nd[] == 2
      if nc[] == nv[]
        # visualize as a simple mesh so that
        # colors can be specified at vertices
        vizmesh2D!(plot)
      else
        # visualize as built-in heatmap
        x = Makie.@lift $xyz[1]
        y = Makie.@lift $xyz[2]
        C = Makie.@lift reshape($colorant, $sz)
        Makie.heatmap!(plot, x, y, C)
      end
    elseif nd[] == 3
      if nc[] == nv[]
        throw(ErrorException("not implemented"))
      else
        # visualize as a simple mesh
        vizmesh3D!(plot)
      end
    else
      throw(ErrorException("can only visualize 2D and 3D grids"))
    end
  end

  # optimized visualization of facets
  if showfacets[]
    tup = Makie.@lift rectilinearsegments($xyz, $nd)
    args = [Makie.@lift($tup[i]) for i in 1:nd[]]
    Makie.lines!(plot, args..., color=facetcolor, linewidth=segmentsize)
  end
end

# helper function to create a minimum number
# of line segments within Cartesian grid
function rectilinearsegments(xyz, nd)
  if nd == 2
    xs, ys = xyz
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
  elseif nd == 3
    xs, ys, zs = xyz
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
  else
    throw(ErrorException("not implemented"))
  end
end
