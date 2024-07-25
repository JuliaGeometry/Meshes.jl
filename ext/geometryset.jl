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
    M = Makie.@lift manifold(first($gvec))
    pdim = Makie.@lift paramdim(first($gvec))
    edim = Makie.@lift embeddim(first($gvec))
    vizgset!(plot, M[], Val(pdim[]), Val(edim[]), gvec, colors)
  end
end

const ObservableVector{T} = Makie.Observable{<:AbstractVector{T}}

function vizgset!(plot, ::Type{<:🌐}, pdim::Val, edim::Val, geoms, colorant)
  vizgset!(plot, 𝔼, pdim, edim, geoms, colorant)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{0}, ::Val, geoms, colorant)
  points = Makie.@lift pointify.($geoms)
  vizmany!(plot, points, colorant)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{0}, ::Val, geoms::ObservableVector{<:Point}, colorant)
  pointsize = plot[:pointsize]

  # get coordinates
  coords = Makie.@lift map(p -> ustrip.(to(p)), $geoms)

  # visualize points
  Makie.scatter!(plot, coords, color=colorant, markersize=pointsize, overdraw=true)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val, geoms, colorant)
  meshes = Makie.@lift discretize.($geoms)
  vizmany!(plot, meshes, colorant)
  showfacets1D!(plot, geoms)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val, geoms::ObservableVector{<:Ray}, colorant)
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

function vizgset!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val, geoms, colorant)
  meshes = Makie.@lift discretize.($geoms)
  vizmany!(plot, meshes, colorant)
  showfacets2D!(plot, geoms)
end

const PolygonLike = Union{Polygon,MultiPolygon}

function vizgset!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val{2}, geoms::ObservableVector{<:PolygonLike}, colorant)
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
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{3}, ::Val, geoms, colorant)
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

asmakie(p::Point) = Makie.Point(ustrip.(Tuple(to(p))))

asmakie(v::Vec) = Makie.Vec(ustrip.(Tuple(v)))
