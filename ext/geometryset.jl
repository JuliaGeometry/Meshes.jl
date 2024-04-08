# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{GeometrySet}})
  gset = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # get geometries
  geoms = Makie.@lift parent($gset)

  # get geometry types
  types = Makie.@lift unique(typeof.($geoms))

  for G in types[]
    inds = Makie.@lift findall(g -> g isa G, $geoms)
    gvec = Makie.@lift collect(G, $geoms[$inds])
    colors = Makie.@lift $colorant isa AbstractVector ? $colorant[$inds] : $colorant
    rank = Makie.@lift paramdim(first($gvec))
    if rank[] == 0
      vizgset0D!(plot, gvec, colors)
    elseif rank[] == 1
      vizgset1D!(plot, gvec, colors)
    elseif rank[] == 2
      vizgset2D!(plot, gvec, colors)
    elseif rank[] == 3
      vizgset3D!(plot, gvec, colors)
    end
  end
end

function Makie.plot!(plot::Viz{<:Tuple{PointSet}})
  pset = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  pointsize = plot[:pointsize]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # get geometries and coordinates
  geoms = Makie.@lift parent($pset)
  coords = Makie.@lift coordinates.($geoms)

  # visualize point set
  Makie.scatter!(plot, coords, color=colorant, markersize=pointsize, overdraw=true)
end

const ObservableVector{T} = Makie.Observable{<:AbstractVector{T}}

function vizgset0D!(plot, geoms, colorant)
  points = Makie.@lift pointify.($geoms)
  vizmany!(plot, points, colorant)
end

function vizgset1D!(plot, geoms, colorant)
  meshes = Makie.@lift discretize.($geoms)
  vizmany!(plot, meshes, colorant)
  showfacets1D!(plot, geoms)
end

function vizgset1D!(plot, geoms::ObservableVector{<:Ray}, colorant)
  rset = plot[:object]
  segmentsize = plot[:segmentsize]

  Dim = embeddim(rset[])

  Dim ∈ (2, 3) || error("not implemented")

  # visualize as built-in arrows
  orig = Makie.@lift [asmakie(ray(0)) for ray in $geoms]
  dirs = Makie.@lift [asmakie(ray(1) - ray(0)) for ray in $geoms]
  size = Makie.@lift 0.1 * $segmentsize
  Makie.arrows!(plot, orig, dirs, color=colorant, arrowsize=size)

  showfacets1D!(plot, geoms)
end

function vizgset2D!(plot, geoms, colorant)
  meshes = Makie.@lift discretize.($geoms)
  vizmany!(plot, meshes, colorant)
  showfacets2D!(plot, geoms)
end

const PolygonLike{Dim,T} = Union{Polygon{Dim,T},MultiPolygon{Dim,T}}

function vizgset2D!(plot, geoms::ObservableVector{<:PolygonLike{2}}, colorant)
  showsegments = plot[:showsegments]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  # repeat colors if necessary
  colors = Makie.@lift mayberepeat($colorant, $geoms)

  # visualize as built-in poly
  polys = Makie.@lift asmakie($geoms)
  if showsegments[]
    Makie.poly!(plot, polys, color=colors, strokecolor=segmentcolor, strokewidth=segmentsize)
  else
    Makie.poly!(plot, polys, color=colors)
  end

  showfacets2D!(plot, geoms)
end

function vizgset3D!(plot, geoms, colorant)
  meshes = Makie.@lift discretize.(boundary.($geoms))
  vizmany!(plot, meshes, colorant)
end

function showfacets1D!(plot, geoms)
  showpoints = plot[:showpoints]
  pointcolor = plot[:pointcolor]
  pointsize = plot[:pointsize]

  if showpoints[]
    # all boundaries are points or multipoints
    bounds = Makie.@lift filter(!isnothing, boundary.($geoms))
    bset = Makie.@lift GeometrySet($bounds)
    viz!(plot, bset, color=pointcolor, pointsize=pointsize)
  end
end

function showfacets2D!(plot, geoms)
  showsegments = plot[:showsegments]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  if showsegments[]
    # all boundaries are 1D geometries
    bounds = Makie.@lift filter(!isnothing, boundary.($geoms))
    bset = Makie.@lift GeometrySet($bounds)
    viz!(plot, bset, color=segmentcolor, segmentsize=segmentsize)
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
