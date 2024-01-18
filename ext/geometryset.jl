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
      # all boundaries are multipoints
      points = Makie.@lift mapreduce(parent, vcat, $bounds)
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

const PolygonSet{Dim,T} = GeometrySet{Dim,T,<:Polygon{Dim,T}}

function Makie.plot!(plot::Viz{<:Tuple{PolygonSet{2}}})
  pset = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  segmentsize = plot[:segmentsize]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # visualize as built-in poly
  polys = Makie.@lift asgbpoly.(parent($pset))
  if showfacets[]
    Makie.poly!(plot, polys, color=colorant, strokecolor=facetcolor, strokewidth=segmentsize)
  else
    Makie.poly!(plot, polys, color=colorant)
  end
end

# converts Meshes.Polygon to GeometryBasics.Polygon
function asgbpoly(poly)
  rs = rings(poly)
  outer = first(rs)
  opts = [asgbpoint(p) for p in vertices(outer)]
  if hasholes(poly)
    ipts = map(i -> [asgbpoint(p) for p in vertices(rs[i])], 2:length(rs))
    Makie.Polygon(opts, ipts)
  else
    Makie.Polygon(opts)
  end
end

# converts Meshes.Point to GeometryBasics.Point
asgbpoint(p::Point{Dim,T}) where {Dim,T} = Makie.Point{Dim,T}(Tuple(coordinates(p)))
