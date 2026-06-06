# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{GeometrySet}})
  colorant!(plot)
  vizgset!(plot)
end

# split heterogeneous geometry sets into homogeneous vectors
# of geometries and send these vectors to specialized recipes
function vizgset!(plot)
  # retrieve geometries and their types
  geoms = parent(plot.object[])
  types = unique(map(typeof, geoms))

  for (i, G) in enumerate(types)
    # add nodes to compute graph for this type of geometry
    gvecid = Symbol(:geoms, i)
    cvecid = Symbol(:colors, i)
    Makie.map!(plot, [:object, :colorant], [gvecid, cvecid]) do gset, colorant
      geos = parent(gset)
      inds = findall(x -> x isa G, geos)
      gvec = collect(G, geos[inds])
      cvec = colorant isa AbstractVector ? colorant[inds] : fill(colorant, length(inds))
      gvec, cvec
    end

    # dispatch specialized recipe for this type of geometry
    g = first(plot[gvecid][])
    M = manifold(g)
    pdim = paramdim(g)
    edim = embeddim(g)
    vizgset!(plot, M, Val(pdim), Val(edim), G, gvecid, cvecid)
  end
end

# ---------------
# IMPLEMENTATION
# ---------------

# fallback to visualization of discretized geometries
function vizgset!(plot, M::Type, pdim::Val, ::Val, ::Type{<:Geometry}, gvecid::Symbol, cvecid::Symbol)
  # visualize 3D geometries as boundary meshes
  maybeboundary = pdim === Val(3) ? boundary : identity

  # refine meshes over the 🌐 manifold until
  # they satisfy the maximum length criterion
  mayberefine = M === 🌐 ? refinemaxlen : identity

  # make sure that the geometries are refined with
  # efficient "grid-like" methods before splitting
  # the resulting quadrangles into triangles
  triangulate = simplexify ∘ mayberefine ∘ discretize ∘ maybeboundary

  # add node to compute graph for triangle meshes
  meshesid = Symbol(gvecid, :_meshes)
  Makie.map!(plot, gvecid, meshesid) do gvec
    map(triangulate, gvec)
  end

  # visualize geometries
  vizmany!(plot, plot[meshesid], plot[cvecid])

  # visualize facets if requested
  pdim === Val(1) && plot.showpoints[] && vizfacets!(plot, gvecid)
  pdim === Val(2) && plot.showsegments[] && vizfacets!(plot, gvecid)
end

# collect and visualize parents of multi-geometries
function vizgset!(plot, M::Type, pdim::Val, edim::Val, ::Type{<:Multi}, gvecid::Symbol, cvecid::Symbol)
  # repeat colors for parents
  parentsid = Symbol(gvecid, :_parents)
  pcolorsid = Symbol(cvecid, :_pcolors)
  Makie.map!(plot, [gvecid, cvecid], [parentsid, pcolorsid]) do gvec, cvec
    parents = mapreduce(parent, vcat, gvec)
    pcolors = [cvec[i] for (i, g) in enumerate(gvec) for _ in 1:length(parent(g))]
    parents, pcolors
  end

  # visualize as set of parents
  P = typeof(first(plot[parentsid][]))
  vizgset!(plot, M, pdim, edim, P, parentsid, pcolorsid)
end

function vizgset!(plot, ::Type, ::Val, ::Val, ::Type{<:Point}, gvecid::Symbol, cvecid::Symbol)
  # get raw Cartesian coordinates of points
  coordsid = Symbol(gvecid, :_coords)
  Makie.map!(plot, gvecid, coordsid) do gvec
    map(p -> ustrip.(to(p)), gvec)
  end

  # visualize points with given marker and size
  Makie.scatter!(
    plot,
    plot[coordsid],
    color=plot[cvecid],
    marker=plot.pointmarker,
    markersize=plot.pointsize,
    overdraw=true
  )
end

function vizgset!(plot, ::Type, ::Val, edim::Val, ::Type{<:Ray}, gvecid::Symbol, cvecid::Symbol)
  # visualize as built-in arrows
  origid = Symbol(gvecid, :_orig)
  dirsid = Symbol(gvecid, :_dirs)
  Makie.map!(plot, gvecid, [origid, dirsid]) do gvec
    orig = [asmakie(ray(0)) for ray in gvec]
    dirs = [asmakie(ray(1) - ray(0)) for ray in gvec]
    orig, dirs
  end

  if edim === Val(2)
    tipwidthid = Symbol(gvecid, :_tipwidth)
    shaftwidthid = Symbol(gvecid, :_shaftwidth)
    Makie.map!(plot, :segmentsize, [tipwidthid, shaftwidthid]) do sz
      tw = 5 * sz
      sw = 0.2 * tw
      tw, sw
    end
    Makie.arrows2d!(
      plot,
      plot[origid],
      plot[dirsid],
      color=plot[cvecid],
      tipwidth=plot[tipwidthid],
      shaftwidth=plot[shaftwidthid]
    )
  elseif edim === Val(3)
    tipradiusid = Symbol(gvecid, :_tipradius)
    shaftradiusid = Symbol(gvecid, :_shaftradius)
    Makie.map!(plot, :segmentsize, [tipradiusid, shaftradiusid]) do sz
      tr = 0.05 * sz
      sr = 0.5 * tr
      tr, sr
    end
    Makie.arrows3d!(
      plot,
      plot[origid],
      plot[dirsid],
      color=plot[cvecid],
      tipradius=plot[tipradiusid],
      shaftradius=plot[shaftradiusid]
    )
  else
    error("not implemented")
  end

  if plot.showpoints[]
    vizfacets!(plot, gvecid)
  end
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val, ::Val{2}, ::Type{<:Line}, gvecid::Symbol, cvecid::Symbol)
  # split vertical and non-vertical lines
  interid = Symbol(gvecid, :_inter)
  vindsid = Symbol(gvecid, :_vinds)
  dindsid = Symbol(gvecid, :_dinds)
  Makie.map!(plot, gvecid, [interid, vindsid, dindsid]) do gvec
    inter = [line ∩ Line((0, 0), (0, 1)) for line in gvec]
    vinds = findall(g -> isnothing(g) || g isa Line, inter)
    dinds = setdiff(1:length(gvec), vinds)
    inter, vinds, dinds
  end

  # split colors accordingly
  vcolorid = Symbol(cvecid, :_vcolor)
  dcolorid = Symbol(cvecid, :_dcolor)
  Makie.map!(plot, [cvecid, vindsid, dindsid], [vcolorid, dcolorid]) do cvec, vinds, dinds
    cvec[vinds], cvec[dinds]
  end

  # compute coordinates and slopes of lines
  xcoordid = Symbol(gvecid, :_xcoord)
  ycoordid = Symbol(gvecid, :_ycoord)
  slopesid = Symbol(gvecid, :_slopes)
  Makie.map!(plot, [gvecid, interid, vindsid, dindsid], [xcoordid, ycoordid, slopesid]) do gvec, inter, vinds, dinds
    # x coordinates of vertical lines
    xcoords = map(gvec[vinds]) do vline
      c = coords(vline(0))
      x = convert(Cartesian, c).x
      ustrip(x)
    end

    # y coordinates of non-vertical lines
    ycoords = map(inter[dinds]) do I
      y = if I isa Line # horizontal line through origin
        zero(Meshes.lentype(I))
      else # intersection point with vertical axis
        convert(Cartesian, coords(I)).y
      end
      ustrip(y)
    end

    # slopes of non-vertical lines
    slopes = map(gvec[dinds]) do dline
      c1 = convert(Cartesian, coords(dline(0)))
      c2 = convert(Cartesian, coords(dline(1)))
      (c2.y - c1.y) / (c2.x - c1.x)
    end

    xcoords, ycoords, slopes
  end

  # visualize vertical and non-vertical lines
  Makie.vlines!(plot, plot[xcoordid], color=plot[vcolorid], linewidth=plot.segmentsize)
  Makie.ablines!(plot, plot[ycoordid], plot[slopesid], color=plot[dcolorid], linewidth=plot.segmentsize)
end

function vizgset!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val{2}, ::Type{<:Polygon}, gvecid::Symbol, cvecid::Symbol)
  # visualize as built-in polygons
  polysid = Symbol(gvecid, :_polys)
  Makie.map!(plot, gvecid, polysid) do gvec
    map(asmakie, gvec)
  end
  if plot.showsegments[]
    Makie.poly!(plot, plot[polysid], color=plot[cvecid], strokecolor=plot.segmentcolor, strokewidth=plot.segmentsize)
  else
    Makie.poly!(plot, plot[polysid], color=plot[cvecid])
  end
end

# -------
# FACETS
# -------

function vizfacets!(plot::Viz{<:Tuple{GeometrySet}}, gvecid::Symbol)
  g = first(plot[gvecid][])
  M = manifold(g)
  pdim = paramdim(g)
  edim = embeddim(g)
  vizgsetfacets!(plot, M, Val(pdim), Val(edim), gvecid)
end

function vizgsetfacets!(plot, ::Type, ::Val{1}, ::Val, gvecid::Symbol)
  # all boundaries are points or multipoints
  psetid = Symbol(gvecid, :_pset)
  Makie.map!(plot, [gvecid], psetid) do gvec
    points = filter(!isnothing, map(boundary, gvec))
    GeometrySet(points)
  end
  viz!(plot, plot[psetid], color=plot.pointcolor, pointmarker=plot.pointmarker, pointsize=plot.pointsize)
end

function vizgsetfacets!(plot, ::Type, ::Val{2}, ::Val, gvecid::Symbol)
  # all boundaries are 1D geometries
  bsetid = Symbol(gvecid, :_bset)
  Makie.map!(plot, [gvecid], bsetid) do gvec
    bounds = filter(!isnothing, map(boundary, gvec))
    GeometrySet(bounds)
  end
  viz!(plot, plot[bsetid], color=plot.segmentcolor, segmentsize=plot.segmentsize)
end
