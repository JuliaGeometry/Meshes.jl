# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plot!(plot::Viz{<:Tuple{Mesh}}) = vizmesh!(plot)

function vizmesh!(plot)
  # retrieve mesh and plot attributes
  mesh = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]

  # retrieve manifold and dimensions
  M = Makie.@lift manifold($mesh)
  pdim = Makie.@lift paramdim($mesh)
  edim = Makie.@lift embeddim($mesh)

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  vizmesh!(plot, M[], Val(pdim[]), Val(edim[]), mesh, colorant)
end

# ---------------
# IMPLEMENTATION
# ---------------

function vizmesh!(plot, ::Type{<:🌐}, pdim::Val, edim::Val, mesh, colorant)
  # fallback to Euclidean recipes because Makie doesn't provide
  # more specific recipes for spherical geometries currently
  vizmesh!(plot, 𝔼, pdim, edim, mesh, colorant)
end

function vizmesh!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val, mesh, colorant)
  segmentsize = plot[:segmentsize]

  colors = Makie.@lift $colorant isa AbstractVector ? $colorant : fill($colorant, nelements($mesh))

  # retrieve coordinates of vertices
  xyzs = Makie.@lift map(p -> ustrip.(to(p)), vertices($mesh))

  # retrieve connectivities of segments
  inds = Makie.@lift map(collect ∘ indices, elements(topology($mesh)))

  # extract coordinates and colors of segments
  scoords = Makie.@lift [$xyzs[indsₑ] for indsₑ in $inds]
  scolors = Makie.@lift [fill($colors[e], length(indsₑ)) for (e, indsₑ) in enumerate($inds)]

  # sentinel coordinates
  nan = Makie.@lift let
    v = first($xyzs)
    typeof(v)(ntuple(i -> NaN, length(v)))
  end

  # splice sentinel coordinates to get discrete colors
  lcoords = Makie.@lift reduce((xyz₁, xyz₂) -> [xyz₁; [$nan]; xyz₂], $scoords)
  lcolors = Makie.@lift reduce((col₁, col₂) -> [col₁; [first(col₁)]; col₂], $scolors)

  # visualize segments
  Makie.lines!(plot, lcoords, color=lcolors, linewidth=segmentsize)
end

function vizmesh!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val, mesh, colorant)
  showsegments = plot[:showsegments]

  # retrieve triangle mesh parameters
  tparams = Makie.@lift let
    # relevant settings
    dim = embeddim($mesh)
    nvert = nvertices($mesh)
    nelem = nelements($mesh)
    verts = eachvertex($mesh)
    topo = topology($mesh)
    elems = elements(topo)

    # coordinates of vertices
    coords = map(asmakie, verts)

    # fan triangulation (assume convexity)
    ntris = sum(e -> nvertices(pltype(e)) - 2, elems)
    tris = Vector{GB.TriangleFace{Int}}(undef, ntris)
    tind = 0
    for elem in elems
      I = indices(elem)
      for i in 2:(length(I) - 1)
        tind += 1
        tris[tind] = GB.TriangleFace(I[1], I[i], I[i + 1])
      end
    end

    # element vs. vertex coloring
    if $colorant isa AbstractVector
      ncolor = length($colorant)
      if ncolor == nelem # element coloring
        # duplicate vertices and adjust
        # connectivities to avoid linear
        # interpolation of colors
        tind = 0
        elem4tri = Dict{Int,Int}()
        sizehint!(elem4tri, ntris)
        for (eind, e) in enumerate(elems)
          for _ in 1:(nvertices(pltype(e)) - 2)
            tind += 1
            elem4tri[tind] = eind
          end
        end
        nv = 3ntris
        tcoords = [coords[i] for tri in tris for i in tri]
        tconnec = [GB.TriangleFace(i, i + 1, i + 2) for i in range(start=1, step=3, length=ntris)]
        tcolors = map(1:nv) do i
          t = ceil(Int, i / 3)
          e = elem4tri[t]
          $colorant[e]
        end
      elseif ncolor == nvert # vertex coloring
        # nothing needs to be done because
        # this is the default in Makie and
        # because the triangulation above
        # does not change the vertices in
        # the original polygonal mesh
        tcoords = coords
        tconnec = tris
        tcolors = $colorant
      else
        throw(ArgumentError("Provided $ncolor colors but the mesh has
                            $nvert vertices and $nelem elements."))
      end
    else # single color
      # nothing needs to be done
      tcoords = coords
      tconnec = tris
      tcolors = $colorant
    end

    # enable shading in 3D
    tshading = dim == 3

    tcoords, tconnec, tcolors, tshading
  end

  # unpack observable of parameters
  tcoords = Makie.@lift $tparams[1]
  tconnec = Makie.@lift $tparams[2]
  tcolors = Makie.@lift $tparams[3]
  tshading = Makie.@lift $tparams[4]

  # Makie's triangle mesh
  mkemesh = Makie.@lift GB.Mesh($tcoords, $tconnec)

  # main visualization
  Makie.mesh!(plot, mkemesh, color=tcolors, shading=tshading)

  if showsegments[]
    vizfacets!(plot)
  end
end

function vizmesh!(plot, ::Type{<:𝔼}, ::Val{3}, ::Val, mesh, colorant)
  meshes = Makie.@lift let
    geoms = elements($mesh)
    bounds = boundary.(geoms)
    discretize.(bounds)
  end
  colors = Makie.@lift $colorant isa AbstractVector ? $colorant : fill($colorant, nelements($mesh))
  vizmany!(plot, meshes, colors)
end

# -------
# FACETS
# -------

function vizfacets!(plot::Viz{<:Tuple{Mesh}})
  mesh = plot[:object]
  M = Makie.@lift manifold($mesh)
  pdim = Makie.@lift paramdim($mesh)
  edim = Makie.@lift embeddim($mesh)
  vizmeshfacets!(plot, M[], Val(pdim[]), Val(edim[]))
end

function vizmeshfacets!(plot, ::Type, ::Val{2}, ::Val)
  mesh = plot[:object]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  # retrieve raw coordinates
  coords = Makie.@lift let
    # relevant settings
    ℒ = Meshes.lentype($mesh)
    T = Unitful.numtype(ℒ)
    dim = embeddim($mesh)
    topo = topology($mesh)
    nvert = nvertices($mesh)
    verts = vertices($mesh)

    # extract coordinates and insert sentinel
    # vertex with NaN coordinates at the end
    xyz = map(p -> ustrip.(to(p)), verts)
    push!(xyz, SVector(ntuple(i -> T(NaN), dim)))

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

  Makie.lines!(plot, coords, color=segmentcolor, linewidth=segmentsize)
end
