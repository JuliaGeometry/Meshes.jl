# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{GeometrySet}})
  gset = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  pointsize = plot[:pointsize]
  segmentsize = plot[:segmentsize]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # collect geometries in domain
  geoms = Makie.@lift collect($gset)

  # retrieve parametric dimension
  ranks = Makie.@lift paramdim.($geoms)

  if all(ranks[] .== 0)
    points = Makie.@lift pointify.($geoms)
    vizmany!(plot, points)
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
        rset = Makie.@lift GeometrySet($geoms[$inds])
        cset = if colorant[] isa AbstractVector
          Makie.@lift $colorant[$inds]
        else
          colorant
        end
        viz!(plot, rset, color=cset)
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
      viz!(plot, (Makie.@lift GeometrySet($points)), color=facetcolor, pointsize=pointsize)
    elseif all(ranks[] .== 2)
      # all boundaries are geometries
      viz!(plot, (Makie.@lift GeometrySet($bounds)), color=facetcolor, segmentsize=segmentsize)
    elseif all(ranks[] .== 3)
      # we already visualized the boundaries because
      # that is all we can do with 3D geometries
    end
  end
end

function Makie.plot!(plot::Viz{<:Tuple{PointSet}})
  pset = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  pointsize = plot[:pointsize]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # collect geometries in domain
  geoms = Makie.@lift collect($pset)
  coords = Makie.@lift coordinates.($geoms)

  # visualize point set
  Makie.scatter!(plot, coords, color=colorant, markersize=pointsize, overdraw=true)
end
