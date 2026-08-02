# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{Mesh}})
  colorant!(plot)
  vizmesh!(plot)
end

function vizmesh!(plot)
  mesh = plot.object[]
  M = manifold(mesh)
  pdim = paramdim(mesh)
  edim = embeddim(mesh)
  vizmesh!(plot, M, Val(pdim), Val(edim))
end

# ---------------
# IMPLEMENTATION
# ---------------

function vizmesh!(plot, ::Type, ::Val{1}, ::Val)
  # retrieve segment coordinates and colors
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

function vizmesh!(plot, ::Type, ::Val{2}, edim::Val)
  # compute triangle mesh and colors
  Makie.map!(plot, [:object, :colorant], [:tmesh]) do mesh, colorant
    # relevant settings
    nvert = nvertices(mesh)
    nelem = nelements(mesh)
    verts = map(asmakie, eachvertex(mesh))
    elems = elements(topology(mesh))

    facecolors = false
    if colorant isa AbstractVector
      colors = Makie.RGBAf.(colorant) # GB.Mesh needs RGBA(f) type
      ncolor = length(colorant)
      if ncolor == nelem # element coloring
        facecolors = true
      elseif ncolor != nvert
        throw(ArgumentError("provided $ncolor colors but the mesh has
                             $nvert vertices and $nelem elements."))
      end
    else
      colors = fill(Makie.RGBAf(colorant), nvert)
    end

    normals_in = plot.normals[]
    normals = nothing
    nnormals = isnothing(normals_in) ? 0 : length(normals_in)
    facenormals = nnormals == nelem
    if nnormals != 0 && nnormals != nelem && nnormals != nvert
        throw(ArgumentError("provided $nnormals normals but the mesh has
                             $nvert vertices and $nelem elements."))
    end

    # decide whether or not to reverse connectivity list and normals
    rev = identity
    normalsign = 1.0
    if crs(mesh) <: LatLon && orientation(first(mesh)) == CW
      rev = reverse
      normalsign = -1.0
    end

    # fan triangulation (assume convexity)
    ntri = sum(e -> nvertices(pltype(e)) - 2, elems)
    tris = Vector{GB.TriangleFace{Int}}(undef, ntri)
    tind = 0
    # initialize normals data
    if nnormals != 0
      if facenormals
        normals = Vector{Makie.Vec{3, Float32}}(undef, ntri)
      else
        # vertex normals can be copied
        normals = map(asmakie, normals_in) .* normalsign
      end
    end
    if facecolors
      tcolors = Vector{typeof(first(colors))}(undef, ntri)
    else # vertex coloring
      # nothing needs to be done because
      # this is the default in Makie and
      # because the triangulation below
      # does not change the vertices in
      # the original polygonal mesh
      tcolors = colors
    end
    for (eind, elem) in enumerate(elems)
      I = rev(indices(elem))
      for i in 2:(length(I) - 1)
        tind += 1
        tris[tind] = GB.TriangleFace(I[1], I[i], I[i + 1])
        if facenormals
          normals[tind] = normalsign * asmakie(normals_in[eind])
        end
        if facecolors
          tcolors[tind] = colors[eind]
        end
      end
    end

    if facenormals
      normals = GB.FaceView(normals, GB.GLTriangleFace.(eachindex(normals)))
    end

    # element vs. vertex coloring
    if facecolors
      tcolors = GB.FaceView(tcolors, GB.GLTriangleFace.(eachindex(tcolors)))
    end

    # triangle mesh
    tmesh = GB.mesh(verts, tris; color=tcolors, normal=normals)
    return (tmesh,)
  end

  # enable shading in 3D
  shading = edim == Val(3)

  # visualize as triangle mesh
  Makie.mesh!(plot, plot.tmesh, shading=shading)

  if plot.showsegments[]
    vizfacets!(plot)
  end
end

function vizmesh!(plot, ::Type, ::Val{3}, ::Val)
  Makie.map!(plot, [:object, :colorant], [:bmesh, :bcolor]) do mesh, colorant
    # discretize boundaries of elements
    meshes = map(discretize ∘ boundary, mesh)
    colors = colorant isa AbstractVector ? colorant : fill(colorant, nelements(mesh))

    # merge into single boundary mesh
    bmesh = reduce(merge, meshes)
    bcolor = [colors[i] for (i, m) in enumerate(meshes) for _ in 1:nelements(m)]

    bmesh, bcolor
  end
  viz!(plot, plot.bmesh, color=plot.bcolor)
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
  # retrieve raw coordinates
  Makie.map!(plot, :object, :coords) do mesh
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

  # visualize as built-in lines
  Makie.lines!(plot, plot.coords, color=plot.segmentcolor, linewidth=plot.segmentsize)
end
