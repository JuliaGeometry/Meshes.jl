# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{GeometrySet}})
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)
  vizgset!(plot)
end

# split heterogeneous geometry sets into homogeneous vectors
# of geometries and send these vectors to specialized recipes
function vizgset!(plot)
  # get geometries
  geoms = parent(plot.object[])

  # get geometry types
  types = unique(map(typeof, geoms))

  for (i, G) in enumerate(types)
    gvec_sym = Symbol(:gvec_, i)
    cvec_sym = Symbol(:cvec_, i)

    Makie.map!(plot, [:object, :colorant], [gvec_sym, cvec_sym]) do g, colorant
      inds = findall(x -> x isa G, parent(g))
      gvec = collect(G, parent(g)[inds])
      cvec = colorant isa AbstractVector ? colorant[inds] : fill(colorant, length(inds))
      gvec, cvec
    end
    gvec = plot[gvec_sym][]
    M = manifold(first(gvec))
    pdim = paramdim(first(gvec))
    edim = embeddim(first(gvec))
    vizgset!(plot, M, Val(pdim), Val(edim), G, gvec_sym, cvec_sym)
  end
end

# ---------------
# IMPLEMENTATION
# ---------------

# fallback to visualization of discretized geometries
function vizgset!(plot, M::Type, pdim::Val, ::Val, G::Type{<:Geometry}, geoms::Symbol, colors::Symbol)
  # refine meshes over the 🌐 manifold until
  # they satisfy the maximum length criterion
  mayberefine = M === 🌐 ? refinemaxlen : identity

  # make sure that the geometries are refined with
  # efficient "grid-like" methods before splitting
  # the resulting quadrangles into triangles
  triangulate = simplexify ∘ mayberefine ∘ discretize

  meshes_sym = Symbol(geoms, :_meshes)

  if pdim === Val(1)
    Makie.map!(plot, [geoms], meshes_sym) do gvec
      map(triangulate, gvec)
    end
    vizmany!(plot, plot[meshes_sym], plot[colors])
    if plot.showpoints[]
      vizfacets!(plot, geoms)
    end
  elseif pdim === Val(2)
    Makie.map!(plot, [geoms], meshes_sym) do gvec
      map(triangulate, gvec)
    end
    vizmany!(plot, plot[meshes_sym], plot[colors])
    if plot.showsegments[]
      vizfacets!(plot, geoms)
    end
  elseif pdim === Val(3)
    Makie.map!(plot, [geoms], meshes_sym) do gvec
      map(triangulate ∘ boundary, gvec)
    end
    vizmany!(plot, plot[meshes_sym], plot[colors])
  end
end

# collect and visualize parents of multi-geometries
function vizgset!(plot, M::Type, pdim::Val, edim::Val, G::Type{<:Multi}, geoms::Symbol, colors::Symbol)
  parents_sym = Symbol(geoms, :_parents)
  pcolors_sym = Symbol(colors, :_pcolors)

  # repeat colors for parents
  Makie.map!(plot, [colors, geoms], [pcolors_sym, parents_sym]) do cvec, gvec
    pcolors = [cvec[i] for (i, g) in enumerate(gvec) for _ in 1:length(parent(g))]
    pgeoms = mapreduce(parent, vcat, gvec)
    pcolors, pgeoms
  end

  P = typeof(first(parent(first(plot[geoms][]))))

  # call recipe for parents
  vizgset!(plot, M, pdim, edim, P, parents_sym, pcolors_sym)
end

function vizgset!(plot, ::Type, ::Val, ::Val, G::Type{<:Point}, geoms::Symbol, colors::Symbol)
  coords_sym = Symbol(geoms, :_coords)

  # get raw Cartesian coordinates of points
  Makie.map!(plot, [geoms], coords_sym) do gvec
    map(p -> ustrip.(to(p)), gvec)
  end

  # visualize points with given marker and size
  Makie.scatter!(
    plot,
    plot[coords_sym],
    color=plot[colors],
    marker=plot.pointmarker,
    markersize=plot.pointsize,
    overdraw=true
  )
end

function vizgset!(plot, ::Type, ::Val, edim::Val, G::Type{<:Ray}, geoms::Symbol, colors::Symbol)
  orig_sym = Symbol(geoms, :_orig)
  dirs_sym = Symbol(geoms, :_dirs)

  # visualize as built-in arrows
  Makie.map!(plot, [geoms], [orig_sym, dirs_sym]) do gvec
    oris = [asmakie(ray(0)) for ray in gvec]
    dirs = [asmakie(ray(1) - ray(0)) for ray in gvec]
    oris, dirs
  end

  if edim === Val(2)
    tipwidth_sym = Symbol(geoms, :_tipwidth)
    shaftwidth_sym = Symbol(geoms, :_shaftwidth)
    Makie.map!(plot, [:segmentsize], [tipwidth_sym, shaftwidth_sym]) do sz
      tw = 5 * sz
      sw = 0.2 * tw
      tw, sw
    end
    Makie.arrows2d!(
      plot,
      plot[orig_sym],
      plot[dirs_sym],
      color=plot[colors],
      tipwidth=plot[tipwidth_sym],
      shaftwidth=plot[shaftwidth_sym]
    )
  elseif edim === Val(3)
    tipradius_sym = Symbol(geoms, :_tipradius)
    shaftradius_sym = Symbol(geoms, :_shaftradius)
    Makie.map!(plot, [:segmentsize], [tipradius_sym, shaftradius_sym]) do sz
      tr = 0.05 * sz
      sr = 0.5 * tr
      tr, sr
    end
    Makie.arrows3d!(
      plot,
      plot[orig_sym],
      plot[dirs_sym],
      color=plot[colors],
      tipradius=plot[tipradius_sym],
      shaftradius=plot[shaftradius_sym]
    )
  else
    error("not implemented")
  end

  if plot.showpoints[]
    vizfacets!(plot, geoms)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val, ::Val{2}, G::Type{<:Line}, geoms::Symbol, colors::Symbol)
  segmentsize = plot.segmentsize

  inter_sym = Symbol(geoms, :_inter)
  vinds_sym = Symbol(geoms, :_vinds)
  dinds_sym = Symbol(geoms, :_dinds)
  vcolor_sym = Symbol(colors, :_vcolor)
  dcolor_sym = Symbol(colors, :_dcolor)

  # split vertical and non-vertical lines
  Makie.map!(plot, [geoms], [inter_sym, vinds_sym, dinds_sym]) do gvec
    inter = [line ∩ Line((0, 0), (0, 1)) for line in gvec]
    vinds = findall(g -> isnothing(g) || g isa Line, inter)
    dinds = setdiff(1:length(gvec), vinds)
    (inter, vinds, dinds)
  end

  # split colors accordingly
  Makie.map!(plot, [colors, vinds_sym], vcolor_sym) do cvec, vinds
    cvec[vinds]
  end
  Makie.map!(plot, [colors, dinds_sym], dcolor_sym) do cvec, dinds
    cvec[dinds]
  end

  xcoord_sym = Symbol(geoms, :_xcoord)
  Makie.map!(plot, [geoms, vinds_sym], xcoord_sym) do gvec, vinds
    vlines = gvec[vinds]
    map(vlines) do vline
      c = coords(vline(0))
      x = convert(Cartesian, c).x
      ustrip(x)
    end
  end
  Makie.vlines!(plot, plot[xcoord_sym], color=plot[vcolor_sym], linewidth=segmentsize)

  ycoord_sym = Symbol(geoms, :_ycoord)
  slopes_sym = Symbol(geoms, :_slopes)
  Makie.map!(plot, [inter_sym, dinds_sym], ycoord_sym) do inter, dinds
    dinter = inter[dinds]
    map(dinter) do I
      y = if I isa Line # horizontal line through origin
        zero(Meshes.lentype(I))
      else # intersection point with vertical axis
        convert(Cartesian, coords(I)).y
      end
      ustrip(y)
    end
  end
  Makie.map!(plot, [geoms, dinds_sym], slopes_sym) do gvec, dinds
    dlines = gvec[dinds]
    map(dlines) do dline
      c1 = convert(Cartesian, coords(dline(0)))
      c2 = convert(Cartesian, coords(dline(1)))
      (c2.y - c1.y) / (c2.x - c1.x)
    end
  end
  Makie.ablines!(plot, plot[ycoord_sym], plot[slopes_sym], color=plot[dcolor_sym], linewidth=segmentsize)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val{2}, G::Type{<:Polygon}, geoms::Symbol, colors::Symbol)
  polys_sym = Symbol(geoms, :_polys)

  # visualize as built-in polygons
  Makie.map!(plot, [geoms], polys_sym) do gvec
    map(asmakie, gvec)
  end

  if plot.showsegments[]
    Makie.poly!(plot, plot[polys_sym], color=plot[colors], strokecolor=plot.segmentcolor, strokewidth=plot.segmentsize)
  else
    Makie.poly!(plot, plot[polys_sym], color=plot[colors])
  end
end

# -------
# FACETS
# -------

function vizfacets!(plot::Viz{<:Tuple{GeometrySet}}, geoms::Symbol)
  gvec = plot[geoms][]
  M = manifold(first(gvec))
  pdim = paramdim(first(gvec))
  edim = embeddim(first(gvec))
  vizgsetfacets!(plot, M, Val(pdim), Val(edim), geoms)
end

function vizgsetfacets!(plot, ::Type, ::Val{1}, ::Val, geoms::Symbol)
  pset_sym = Symbol(geoms, :_pset)

  # all boundaries are points or multipoints
  Makie.map!(plot, [geoms], pset_sym) do gvec
    points = filter(!isnothing, map(boundary, gvec))
    GeometrySet(points)
  end

  viz!(plot, plot[pset_sym], color=plot.pointcolor, pointmarker=plot.pointmarker, pointsize=plot.pointsize)
end

function vizgsetfacets!(plot, ::Type, ::Val{2}, ::Val, geoms::Symbol)
  bset_sym = Symbol(geoms, :_bset)

  # all boundaries are 1D geometries
  Makie.map!(plot, [geoms], bset_sym) do gvec
    bounds = filter(!isnothing, map(boundary, gvec))
    GeometrySet(bounds)
  end

  viz!(plot, plot[bset_sym], color=plot.segmentcolor, segmentsize=plot.segmentsize)
end
