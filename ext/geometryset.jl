# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plottype(::GeometrySet) = Viz{<:Tuple{GeometrySet}}

function Makie.plot!(plot::Viz{<:Tuple{GeometrySet}})
  collection = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  facetcolor = plot[:facetcolor]
  showfacets = plot[:showfacets]
  pointsize = plot[:pointsize]
  segmentsize = plot[:segmentsize]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # decimate geometries if needed
  geoms = Makie.@lift collect($collection)

  # retrieve parametric dimension
  ranks = Makie.@lift paramdim.($geoms)

  if all(ranks[] .== 0)
    # visualize point set
    coords = Makie.@lift coordinates.($geoms)
    Makie.scatter!(plot, coords, color=colorant, markersize=pointsize)
  elseif all(ranks[] .== 1)
    meshes = Makie.@lift discretize.($geoms)
    vizmany!(plot, meshes)
  elseif all(ranks[] .== 2)
    meshes = Makie.@lift discretize.($geoms)
    vizmany!(plot, meshes)
  elseif all(ranks[] .== 3)
    meshes = Makie.@lift discretize.(boundary.($geoms))
    vizmany!(plot, meshes)
  else # mixed dimension
    # visualize subsets of equal rank
    for rank in (3, 2, 1, 0)
      inds = Makie.@lift findall(g -> paramdim(g) == rank, $geoms)
      if !isempty(inds[])
        gset = Makie.@lift GeometrySet($geoms[$inds])
        if colorant[] isa AbstractVector
          cset = Makie.@lift $colorant[$inds]
        else
          cset = colorant
        end
        viz!(plot, gset, color=cset)
      end
    end
  end

  if showfacets[]
    bounds = Makie.@lift filter(!isnothing, boundary.($geoms))
    if isempty(bounds[])
      # nothing to be done
    elseif all(ranks[] .== 1)
      # all boundaries are point sets
      points = Makie.@lift mapreduce(collect, vcat, $bounds)
      viz!(plot, (Makie.@lift GeometrySet($points)), color=facetcolor, showfacets=false, pointsize=pointsize)
    elseif all(ranks[] .== 2)
      # all boundaries are geometries
      viz!(plot, (Makie.@lift GeometrySet($bounds)), color=facetcolor, showfacets=false, segmentsize=segmentsize)
    elseif all(ranks[] .== 3)
      # we already visualized the boundaries because
      # that is all we can do with 3D geometries
    end
  end
end
