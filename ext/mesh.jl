# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plot!(plot::Viz{<:Tuple{Mesh}}) = vizmesh!(plot)

function vizmesh!(plot)
  mesh = plot.object[]
  M = manifold(mesh)
  pdim = paramdim(mesh)
  edim = embeddim(mesh)

  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

  vizmesh!(plot, M, Val(pdim), Val(edim))
end

# ---------------
# IMPLEMENTATION
# ---------------

function vizmesh!(plot, ::Type, ::Val{1}, ::Val)
  # compute line segment coordinates and colors
  Makie.map!(plot, [:object, :colorant], [:lcoords, :lcolors]) do mesh, colorant
    # retrieve colors of segments
    colors = colorant isa AbstractVector ? colorant : fill(colorant, nelements(mesh))

    # retrieve coordinates of vertices
    xyzs = map(p -> ustrip.(to(p)), vertices(mesh))

    # retrieve connectivities of segments
    inds = map(collect ∘ indices, elements(topology(mesh)))

    # sentinel coordinates
    xyz = first(xyzs)
    nan = typeof(xyz)(ntuple(i -> NaN, length(xyz)))

    # extract coordinates and colors of segments
    scoords = [xyzs[indsₑ] for indsₑ in inds]
    scolors = [fill(colors[e], length(indsₑ)) for (e, indsₑ) in enumerate(inds)]

    # splice sentinel coordinates to get discrete colors
    lcoords = reduce((xyz₁, xyz₂) -> [xyz₁; [nan]; xyz₂], scoords)
    lcolors = reduce((col₁, col₂) -> [col₁; [first(col₁)]; col₂], scolors)

    lcoords, lcolors
  end

  # visualize as built-in lines
  Makie.lines!(plot, plot.lcoords, color=plot.lcolors, linewidth=plot.segmentsize)
end

function vizmesh!(plot, ::Type, ::Val{2}, ::Val)
  # retrieve triangle mesh parameters
  Makie.map!(plot, [:object, :colorant], :tparams) do mesh, colorant
    # relevant settings
    edim = embeddim(mesh)
    nvert = nvertices(mesh)
    nelem = nelements(mesh)
    verts = map(asmakie, eachvertex(mesh))
    elems = elements(topology(mesh))

    # decide whether or not to reverse connectivity list
    rev = crs(mesh) <: LatLon && orientation(first(mesh)) == CW ? reverse : identity

    # fan triangulation (assume convexity)
    ntri = sum(e -> nvertices(pltype(e)) - 2, elems)
    tris = Vector{GB.TriangleFace{Int}}(undef, ntri)
    tind = 0
    for elem in elems
      I = rev(indices(elem))
      for i in 2:(length(I) - 1)
        tind += 1
        tris[tind] = GB.TriangleFace(I[1], I[i], I[i + 1])
      end
    end

    # element vs. vertex coloring
    if colorant isa AbstractVector
      ncolor = length(colorant)
      if ncolor == nelem # element coloring
        # duplicate vertices and adjust
        # connectivities to avoid linear
        # interpolation of colors
        tind = 0
        elem4tri = Dict{Int,Int}()
        sizehint!(elem4tri, ntri)
        for (eind, e) in enumerate(elems)
          for _ in 1:(nvertices(pltype(e)) - 2)
            tind += 1
            elem4tri[tind] = eind
          end
        end
        nv = 3ntri
        tverts = [verts[i] for tri in tris for i in tri]
        telems = [GB.TriangleFace(i, i + 1, i + 2) for i in range(start=1, step=3, length=ntri)]
        tcolors = map(1:nv) do i
          t = ceil(Int, i / 3)
          e = elem4tri[t]
          colorant[e]
        end
      elseif ncolor == nvert # vertex coloring
        # nothing needs to be done because
        # this is the default in Makie and
        # because the triangulation above
        # does not change the vertices in
        # the original polygonal mesh
        tverts = verts
        telems = tris
        tcolors = colorant
      else
        throw(ArgumentError("Provided $ncolor colors but the mesh has
                            $nvert vertices and $nelem elements."))
      end
    else # single color
      # nothing needs to be done
      tverts = verts
      telems = tris
      tcolors = colorant
    end

    # enable shading in 3D
    tshading = edim == 3

    tverts, telems, tcolors, tshading
  end

  # unpack observable of parameters
  Makie.map!(plot, [:tparams], [:tverts, :telems, :tcolors, :tshading]) do tparams
    (tparams[1], tparams[2], tparams[3], tparams[4])
  end

  # Makie's triangle mesh
  Makie.map!(GB.Mesh, plot, [:tverts, :telems], :mkemesh)

  # main visualization
  Makie.mesh!(plot, plot.mkemesh, color=plot.tcolors, shading=plot.tshading)

  if plot.showsegments[]
    vizfacets!(plot)
  end
end

function vizmesh!(plot, ::Type, ::Val{3}, ::Val)
  Makie.map!(plot, [:object], :meshes) do mesh
    map(discretize ∘ boundary, mesh)
  end
  Makie.map!(plot, [:colorant, :object], :colors) do colorant, mesh
    colorant isa AbstractVector ? colorant : fill(colorant, nelements(mesh))
  end
  vizmany!(plot, plot.meshes, plot.colors)
end

# -------
# FACETS
# -------

function vizfacets!(plot::Viz{<:Tuple{Mesh}})
  mesh = plot.object[]
  M = manifold(mesh)
  pdim = paramdim(mesh)
  edim = embeddim(mesh)
  vizmeshfacets!(plot, M, Val(pdim), Val(edim))
end

function vizmeshfacets!(plot, ::Type, ::Val{2}, ::Val)
  segmentcolor = plot.segmentcolor
  segmentsize = plot.segmentsize

  # retrieve raw coordinates
  Makie.map!(plot, [:object], :coords) do mesh
    # relevant settings
    ℒ = Meshes.lentype(mesh)
    T = Unitful.numtype(ℒ)
    edim = embeddim(mesh)
    topo = topology(mesh)
    nvert = nvertices(mesh)
    verts = vertices(mesh)

    # extract coordinates and insert sentinel
    # vertex with NaN coordinates at the end
    xyz = map(p -> ustrip.(to(p)), verts)
    push!(xyz, SVector(ntuple(i -> T(NaN), edim)))

    # find indices of incident vertices
    inds = Int[]
    try # efficient algorithm with half-edge topology
      t = convert(HalfEdgeTopology, topo)
      ∂ = Boundary{1,0}(t)
      for i in 1:nfacets(t)
        for j in ∂(i)
          push!(inds, j)
        end
        # interleave with a sentinel index
        push!(inds, nvert + 1)
      end
    catch # brute force algorithm with duplicate edges
      t = topo
      ∂ = Boundary{2,0}(t)
      for i in 1:nelements(t)
        for j in ∂(i)
          push!(inds, j)
        end
        # interleave with a sentinel index
        push!(inds, nvert + 1)
      end
    end

    # extract incident vertices
    xyz[inds]
  end

  Makie.lines!(plot, plot.coords, color=segmentcolor, linewidth=segmentsize)
end
