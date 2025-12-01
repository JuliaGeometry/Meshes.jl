# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const ObservableVector{T} = Makie.Observable{<:AbstractVector{T}}

Makie.plot!(plot::Viz{<:Tuple{GeometrySet}}) = vizgset!(plot)

# split heterogeneous geometry sets into homogeneous vectors
# of geometries and send these vectors to specialized recipes
function vizgset!(plot; facets=false)
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
      cvec = Makie.@lift if $colorant isa AbstractVector
        $colorant[$inds]
      else
        fill($colorant, length($inds))
      end
      vizgset!(plot, M[], Val(pdim[]), Val(edim[]), gvec, cvec)
    end
  end
end

# ---------------
# IMPLEMENTATION
# ---------------

function vizgset!(plot, ::Type{<:🌐}, pdim::Val, edim::Val, geoms::ObservableVector{<:Geometry}, colors)
  # fallback to Euclidean recipes because Makie doesn't provide
  # more specific recipes for spherical geometries currently
  vizgset!(plot, 𝔼, pdim, edim, geoms, colors)
end

function vizgset!(plot, ::Type{<:𝔼}, pdim::Val, edim::Val, geoms::ObservableVector{<:Geometry}, colors)
  showsegments = plot[:showsegments]
  showpoints = plot[:showpoints]

  if pdim === Val(1)
    meshes = Makie.@lift discretize.($geoms)
    vizmany!(plot, meshes, colors)
    if showpoints[]
      vizfacets!(plot, geoms)
    end
  elseif pdim === Val(2)
    meshes = Makie.@lift discretize.($geoms)
    vizmany!(plot, meshes, colors)
    if showsegments[]
      vizfacets!(plot, geoms)
    end
  elseif pdim == Val(3)
    meshes = Makie.@lift discretize.(boundary.($geoms))
    vizmany!(plot, meshes, colors)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{0}, ::Val, geoms::ObservableVector{<:Point}, colors)
  pointmarker = plot[:pointmarker]
  pointsize = plot[:pointsize]

  # get raw Cartesian coordinates of points
  coords = Makie.@lift map(p -> ustrip.(to(p)), $geoms)

  # visualize points with given marker and size
  Makie.scatter!(plot, coords, color=colors, marker=pointmarker, markersize=pointsize, overdraw=true)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val, geoms::ObservableVector{<:Ray}, colors)
  segmentsize = plot[:segmentsize]
  showpoints = plot[:showpoints]

  edim = embeddim(first(geoms[]))

  # visualize as built-in arrows
  orig = Makie.@lift [asmakie(ray(0)) for ray in $geoms]
  dirs = Makie.@lift [asmakie(ray(1) - ray(0)) for ray in $geoms]
  if edim == 2
    tipwidth = Makie.@lift 5 * $segmentsize
    shaftwidth = Makie.@lift 0.2 * $tipwidth
    Makie.arrows2d!(plot, orig, dirs, color=colors, tipwidth=tipwidth, shaftwidth=shaftwidth)
  elseif edim == 3
    tipradius = Makie.@lift 0.05 * $segmentsize
    shaftradius = Makie.@lift 0.5 * $tipradius
    Makie.arrows3d!(plot, orig, dirs, color=colors, tipradius=tipradius, shaftradius=shaftradius)
  else
    error("not implemented")
  end

  if showpoints[]
    vizfacets!(plot, geoms)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val{2}, geoms::ObservableVector{<:Line}, colors)
  segmentsize = plot[:segmentsize]

  # split vertical and non-vertical lines
  inter = Makie.@lift [line ∩ Line((0, 0), (0, 1)) for line in $geoms]
  vinds = Makie.@lift findall(isnothing, $inter)
  dinds = Makie.@lift setdiff(1:length($geoms), $vinds)

  # split colors accordingly
  vcolor = Makie.@lift $colors[$vinds]
  dcolor = Makie.@lift $colors[$dinds]

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

vizgset!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val{2}, geoms::ObservableVector{<:Box}, colors) =
  vizgsetbox𝔼!(plot, geoms, colors)

vizgset!(plot, ::Type{<:𝔼}, ::Val{3}, ::Val{3}, geoms::ObservableVector{<:Box}, colors) =
  vizgsetbox𝔼!(plot, geoms, colors)

function vizgsetbox𝔼!(plot, geoms::ObservableVector{<:Box}, colors)
  showsegments = plot[:showsegments]

  # visualize as built-in boxes
  boxes = Makie.@lift asmakie.($geoms)
  shading = Makie.@lift embeddim(first($geoms)) == 3
  Makie.mesh!(plot, boxes, color=colors, shading=shading)

  if showsegments[]
    vizfacets!(plot, geoms)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val{2}, geoms::ObservableVector{<:Polygon}, colors)
  showsegments = plot[:showsegments]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  # visualize as built-in polygons
  polys = Makie.@lift asmakie.($geoms)
  if showsegments[]
    Makie.poly!(plot, polys, color=colors, strokecolor=segmentcolor, strokewidth=segmentsize)
  else
    Makie.poly!(plot, polys, color=colors)
  end
end

vizgset!(plot, M::Type{<:🌐}, pdim::Val, edim::Val, geoms::ObservableVector{<:Multi}, colors) =
  vizgsetmulti!(plot, M, pdim, edim, geoms, colors)

vizgset!(plot, M::Type{<:𝔼}, pdim::Val, edim::Val, geoms::ObservableVector{<:Multi}, colors) =
  vizgsetmulti!(plot, M, pdim, edim, geoms, colors)

function vizgsetmulti!(plot, M, pdim, edim, geoms, colors)
  # retrieve parent geometries
  parents = Makie.@lift mapreduce(parent, vcat, $geoms)

  # repeat colors for parents
  pcolors = Makie.@lift [$colors[i] for (i, g) in enumerate($geoms) for _ in 1:length(parent(g))]

  # call recipe for parents
  vizgset!(plot, M, pdim, edim, parents, pcolors)
end

# -------
# FACETS
# -------

vizfacets!(plot::Viz{<:Tuple{GeometrySet}}) = vizgset!(plot, facets=false)

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
