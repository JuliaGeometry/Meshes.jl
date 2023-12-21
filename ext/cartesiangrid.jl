# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plottype(::CartesianGrid) = Viz{<:Tuple{CartesianGrid}}

Makie.convert_arguments(::Type{<:Viz}, grid::CartesianGrid) = (grid,)

function Makie.plot!(plot::Viz{<:Tuple{CartesianGrid}})
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

  # origin, spacing and size of grid
  or = Makie.@lift coordinates(minimum($grid))
  sp = Makie.@lift spacing($grid)
  sz = Makie.@lift size($grid)

  # coordinates of centroids
  xyz = [(Makie.@lift range($or[i] + $sp[i] / 2, step=$sp[i], length=$sz[i])) for i in 1:nd[]]

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
        C = Makie.@lift reshape($colorant, $sz)
        Makie.heatmap!(plot, xyz[1], xyz[2], C)
      end
    elseif nd[] == 3
      if nc[] == nv[]
        throw(ErrorException("not implemented"))
      else
        # visualize as built-in meshscatter
        xs = Makie.@lift $(xyz[1]) .+ $sp[1] / 2
        ys = Makie.@lift $(xyz[2]) .+ $sp[2] / 2
        zs = Makie.@lift $(xyz[3]) .+ $sp[3] / 2
        cs = Makie.@lift [(x, y, z) for z in $zs for y in $ys for x in $xs]
        re = Makie.@lift Makie.Rect3(-1 .* $sp, $sp)
        Makie.meshscatter!(plot, cs, marker=re, markersize=1, color=colorant)
      end
    else
      throw(ErrorException("can only visualize 2D and 3D grids"))
    end
  end

  # optimized visualization of facets
  if showfacets[]
    tup = Makie.@lift cartesiansegments($or, $sp, $sz, $nd)
    xyz = [(Makie.@lift $tup[i]) for i in 1:nd[]]
    Makie.lines!(plot, xyz..., color=facetcolor, linewidth=segmentsize)
  end
end

# helper function to create a minimum number
# of line segments within Cartesian grid
function cartesiansegments(or, sp, sz, nd)
  if nd == 2
    xs = range(or[1], step=sp[1], length=sz[1] + 1)
    ys = range(or[2], step=sp[2], length=sz[2] + 1)
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
    xs = range(or[1], step=sp[1], length=sz[1] + 1)
    ys = range(or[2], step=sp[2], length=sz[2] + 1)
    zs = range(or[3], step=sp[3], length=sz[3] + 1)
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
