# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plot!(plot::Viz{<:Tuple{GeometrySet}}) = vizgeoms!(plot)

const ObservableVector{T} = Makie.Observable{<:AbstractVector{T}}

function vizgset!(plot, ::Type{<:🌐}, pdim::Val, edim::Val, geoms, colorant)
  vizgset!(plot, 𝔼, pdim, edim, geoms, colorant)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{0}, ::Val, geoms, colorant)
  points = Makie.@lift aspoints.($geoms)
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

function vizgset!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val, geoms, colorant)
  showpoints = plot[:showpoints]

  meshes = Makie.@lift _discretize.($geoms)
  vizmany!(plot, meshes, colorant)

  if showpoints[]
    vizfacets!(plot, geoms)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val, geoms::ObservableVector{<:Ray}, colorant)
  segmentsize = plot[:segmentsize]
  showpoints = plot[:showpoints]

  Dim = embeddim(first(geoms[]))

  # visualize as built-in arrows
  orig = Makie.@lift [asmakie(ray(0)) for ray in $geoms]
  dirs = Makie.@lift [asmakie(ray(1) - ray(0)) for ray in $geoms]
  if Dim == 2
    tipwidth = Makie.@lift 5 * $segmentsize
    shaftwidth = Makie.@lift 0.2 * $tipwidth
    Makie.arrows2d!(plot, orig, dirs, color=colorant, tipwidth=tipwidth, shaftwidth=shaftwidth)
  elseif Dim == 3
    tipradius = Makie.@lift 0.05 * $segmentsize
    shaftradius = Makie.@lift 0.5 * $tipradius
    Makie.arrows3d!(plot, orig, dirs, color=colorant, tipradius=tipradius, shaftradius=shaftradius)
  else
    error("not implemented")
  end

  if showpoints[]
    vizfacets!(plot, geoms)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val{2}, geoms::ObservableVector{<:Line}, colorant)
  segmentsize = plot[:segmentsize]

  # split vertical and non-vertical lines
  inter = Makie.@lift [line ∩ Line((0, 0), (0, 1)) for line in $geoms]
  vinds = Makie.@lift findall(isnothing, $inter)
  dinds = Makie.@lift setdiff(1:length($geoms), $vinds)

  # split colors accordingly
  if colorant[] isa AbstractVector
    vcolor = Makie.@lift $colorant[$vinds]
    dcolor = Makie.@lift $colorant[$dinds]
  else
    vcolor = colorant
    dcolor = colorant
  end

  # visualize vertical lines
  if !isempty(vinds[])
    vlines = Makie.@lift $geoms[$vinds]
    xcoord = Makie.@lift map($vlines) do vline
      c = coords(vline(0))
      x = convert(Cartesian, c).x
      ustrip(x)
    end
    Makie.vlines!(plot, xcoord, color=vcolor, linewidth=segmentsize)
  end

  # visualize non-vertical lines
  if !isempty(dinds[])
    dlines = Makie.@lift $geoms[$dinds]
    dinter = Makie.@lift $inter[$dinds]
    ycoord = Makie.@lift map($dinter) do point
      c = coords(point)
      y = convert(Cartesian, c).y
      ustrip(y)
    end
    slopes = Makie.@lift map($dlines) do dline
      c1 = convert(Cartesian, coords(dline(0)))
      c2 = convert(Cartesian, coords(dline(1)))
      (c2.y - c1.y) / (c2.x - c1.x)
    end
    Makie.ablines!(plot, ycoord, slopes, color=dcolor, linewidth=segmentsize)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val, geoms, colorant)
  showsegments = plot[:showsegments]

  meshes = Makie.@lift _discretize.($geoms)
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
  meshes = Makie.@lift _discretize.(boundary.($geoms))
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
