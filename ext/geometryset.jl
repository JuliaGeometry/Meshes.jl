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

  # get geometries
  geoms = Makie.@lift parent($gset)

  # retrieve parametric dimension
  ranks = Makie.@lift paramdim.($geoms)

  if all(ranks[] .== 0)
    points = Makie.@lift pointify.($geoms)
    vizmany!(plot, points)
  elseif all(ranks[] .== 1)
    vizgset1D!(plot, geoms)
  elseif all(ranks[] .== 2)
    vizgset2D!(plot, geoms)
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

  # get geometries and coordinates
  geoms = Makie.@lift parent($pset)
  coords = Makie.@lift coordinates.($geoms)

  # visualize point set
  Makie.scatter!(plot, coords, color=colorant, markersize=pointsize, overdraw=true)
end

function vizgset1D!(plot, geoms)
  meshes = Makie.@lift discretize.($geoms)
  vizmany!(plot, meshes)
end

const RaySet{Dim,T} = GeometrySet{Dim,T,Ray{Dim,T}}

function vizgset1D!(plot::Viz{<:Tuple{RaySet}}, geoms)
  rset = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]

  if embeddim(rset[]) ∉ (2, 3)
    error("not implemented")
  end

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # visualize as built-in arrows
  origins = Makie.@lift [asmakie(ray(0)) for ray in $geoms]
  directions = Makie.@lift [asmakie(ray(1) - ray(0)) for ray in $geoms]
  Makie.arrows!(plot, origins, directions, color=colorant)
end

function vizgset2D!(plot, geoms)
  meshes = Makie.@lift discretize.($geoms)
  vizmany!(plot, meshes)
end

const PolygonLike{Dim,T} = Union{Polygon{Dim,T},MultiPolygon{Dim,T}}

const PolygonLikeSet{Dim,T} = GeometrySet{Dim,T,<:PolygonLike{Dim,T}}

function vizgset2D!(plot::Viz{<:Tuple{PolygonLikeSet{2}}}, geoms)
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  segmentsize = plot[:segmentsize]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # repeat colors if necessary
  colors = Makie.@lift mayberepeat($colorant, $geoms)

  # visualize as built-in poly
  polys = Makie.@lift asmakie($geoms)
  if showfacets[]
    Makie.poly!(plot, polys, color=colors, strokecolor=facetcolor, strokewidth=segmentsize)
  else
    Makie.poly!(plot, polys, color=colors)
  end
end

asmakie(geoms::AbstractVector{<:Geometry}) = asmakie.(geoms)

asmakie(multis::AbstractVector{<:Multi}) = mapreduce(m -> asmakie.(parent(m)), vcat, multis)

function asmakie(poly::Polygon)
  rs = rings(poly)
  outer = [asmakie(p) for p in vertices(first(rs))]
  if hasholes(poly)
    inners = map(i -> [asmakie(p) for p in vertices(rs[i])], 2:length(rs))
    Makie.Polygon(outer, inners)
  else
    Makie.Polygon(outer)
  end
end

asmakie(p::Point{Dim,T}) where {Dim,T} = Makie.Point{Dim,T}(Tuple(coordinates(p)))

asmakie(v::Vec{Dim,T}) where {Dim,T} = Makie.Vec{Dim,T}(Tuple(v))
