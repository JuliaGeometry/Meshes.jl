# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plot!(plot::Viz{<:Tuple{GeometrySet}}) = vizgeoms!(plot)

const ObservableVector{T} = Makie.Observable{<:AbstractVector{T}}

function vizgset!(plot, ::Type{<:🌐}, pdim::Val, edim::Val, geoms, colorant)
  vizgset!(plot, 𝔼, pdim, edim, geoms, colorant)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{0}, ::Val, geoms, colorant)
  points = Makie.@lift pointify.($geoms)
  vizmany!(plot, points, colorant)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{0}, ::Val, geoms::ObservableVector{<:Point}, colorant)
  pointmarker = plot[:pointmarker]
  pointsize = plot[:pointsize]

  # get raw Cartesian coordinates of points
  coords = Makie.@lift map(p -> ustrip.(to(p)), $geoms)

  # visualize points with given marker and size
  Makie.scatter!(plot, coords, color=colorant, marker=pointmarker, markersize=pointsize, overdraw=true)
end

function vizgset!(plot, ::Type{<:🌐}, ::Val{1}, ::Val, geoms, colorant)
  showpoints = plot[:showpoints]

  meshes = Makie.@lift begin
    T = numtype(Meshes.lentype(first($geoms)))
    method = MaxLengthDiscretization(T(1000) * u"km")
    [discretize(g, method) for g in $geoms]
  end
  vizmany!(plot, meshes, colorant)

  if showpoints[]
    vizfacets!(plot, geoms)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val, geoms, colorant)
  showpoints = plot[:showpoints]

  meshes = Makie.@lift discretize.($geoms)
  vizmany!(plot, meshes, colorant)

  if showpoints[]
    vizfacets!(plot, geoms)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val, geoms::ObservableVector{<:Ray}, colorant)
  rset = plot[:object]
  segmentsize = plot[:segmentsize]
  showpoints = plot[:showpoints]

  Dim = embeddim(rset[])

  Dim ∈ (2, 3) || error("not implemented")

  # visualize as built-in arrows
  orig = Makie.@lift [asmakie(ray(0)) for ray in $geoms]
  dirs = Makie.@lift [asmakie(ray(1) - ray(0)) for ray in $geoms]
  size = Makie.@lift 0.1 * $segmentsize
  Makie.arrows!(plot, orig, dirs, color=colorant, arrowsize=size)

  if showpoints[]
    vizfacets!(plot, geoms)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val, geoms, colorant)
  showsegments = plot[:showsegments]

  meshes = Makie.@lift discretize.($geoms)
  vizmany!(plot, meshes, colorant)

  if showsegments[]
    vizfacets!(plot, geoms)
  end
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

vizfacets!(plot::Viz{<:Tuple{GeometrySet}}) = vizgeoms!(plot, facets=false)

function vizfacets!(plot::Viz{<:Tuple{GeometrySet}}, geoms)
  M = Makie.@lift manifold(first($geoms))
  pdim = Makie.@lift paramdim(first($geoms))
  edim = Makie.@lift embeddim(first($geoms))
  vizgsetfacets!(plot, M[], Val(pdim[]), Val(edim[]), geoms)
end

function vizgsetfacets!(plot, ::Type, ::Val{1}, ::Val, geoms)
  pointmarker = plot[:pointmarker]
  pointcolor = plot[:pointcolor]
  pointsize = plot[:pointsize]

  # all boundaries are points or multipoints
  points = Makie.@lift filter(!isnothing, boundary.($geoms))
  pset = Makie.@lift GeometrySet($points)
  viz!(plot, pset, color=pointcolor, pointmarker=pointmarker, pointsize=pointsize)
end

function vizgsetfacets!(plot, ::Type, ::Val{2}, ::Val, geoms)
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  # all boundaries are 1D geometries
  bounds = Makie.@lift filter(!isnothing, boundary.($geoms))
  bset = Makie.@lift GeometrySet($bounds)
  viz!(plot, bset, color=segmentcolor, segmentsize=segmentsize)
end

function vizgeoms!(plot; facets=false)
  gset = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]

  # process color spec into colorant
  colorant = facets ? nothing : Makie.@lift(process($color, $colormap, $colorrange, $alpha))

  # get geometries
  geoms = Makie.@lift parent($gset)

  # get geometry types
  types = Makie.@lift unique(typeof.($geoms))

  for G in types[]
    inds = Makie.@lift findall(g -> g isa G, $geoms)
    gvec = Makie.@lift collect(G, $geoms[$inds])
    M = Makie.@lift manifold(first($gvec))
    pdim = Makie.@lift paramdim(first($gvec))
    edim = Makie.@lift embeddim(first($gvec))
    if facets
      vizgsetfacets!(plot, M[], Val(pdim[]), Val(edim[]), gvec)
    else
      cvec = Makie.@lift $colorant isa AbstractVector ? $colorant[$inds] : $colorant
      vizgset!(plot, M[], Val(pdim[]), Val(edim[]), gvec, cvec)
    end
  end
end
