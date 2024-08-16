# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{Mesh}})
  # retrieve mesh and dimensions
  mesh = plot[:object]
  M = Makie.@lift manifold($mesh)
  pdim = Makie.@lift paramdim($mesh)
  edim = Makie.@lift embeddim($mesh)
  vizmesh!(plot, M[], Val(pdim[]), Val(edim[]))
end

function vizmesh!(plot, ::Type{<:🌐}, pdim::Val, edim::Val)
  vizmesh!(plot, 𝔼, pdim, edim)
end

function vizmesh!(plot, ::Type{<:𝔼}, ::Val{1}, ::Val)
  mesh = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  segmentsize = plot[:segmentsize]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # retrieve coordinates of segments
  coords = Makie.@lift let
    topo = topology($mesh)
    vert = vertices($mesh)
    segmentsof(topo, vert)
  end

  # repeat colors for vertices of segments
  colors = Makie.@lift let
    if $colorant isa AbstractVector
      c = [$colorant[e] for e in 1:nelements($mesh) for _ in 1:3]
      c[begin:(end - 1)]
    else
      $colorant
    end
  end

  # visualize segments
  Makie.lines!(plot, coords, color=colors, linewidth=segmentsize)
end

function vizmesh!(plot, ::Type{<:𝔼}, ::Val{2}, ::Val)
  mesh = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  showsegments = plot[:showsegments]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # retrieve triangle mesh parameters
  tparams = Makie.@lift let
    # relevant settings
    dim = embeddim($mesh)
    nvert = nvertices($mesh)
    nelem = nelements($mesh)
    verts = vertices($mesh)
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
    tshading = dim == 3 ? Makie.FastShading : Makie.NoShading

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
    # retrieve coordinates parameters
    xparams = Makie.@lift let
      # relevant settings
      T = Unitful.numtype(Meshes.lentype($mesh))
      dim = embeddim($mesh)
      topo = topology($mesh)
      nvert = nvertices($mesh)
      verts = vertices($mesh)
      coords = map(p -> ustrip.(to(p)), verts)

      # use a sophisticated data structure
      # to extract the edges from the n-gons
      t = convert(HalfEdgeTopology, topo)
      ∂ = Boundary{1,0}(t)

      # append indices of incident vertices
      # interleaved with a sentinel index
      inds = Int[]
      for i in 1:nfacets(t)
        for j in ∂(i)
          push!(inds, j)
        end
        push!(inds, nvert + 1)
      end

      # fill sentinel index with NaN coordinates
      push!(coords, SVector(ntuple(i -> T(NaN), dim)))

      # extract incident vertices
      coords = coords[inds]

      # split coordinates to match signature
      [getindex.(coords, j) for j in 1:dim]
    end

    # unpack observable of paramaters
    xyz = map(1:embeddim(mesh[])) do i
      Makie.@lift $xparams[i]
    end

    Makie.lines!(plot, xyz..., color=segmentcolor, linewidth=segmentsize)
  end
end

function vizmesh!(plot, ::Type{<:𝔼}, ::Val{3}, ::Val)
  mesh = plot[:object]
  color = plot[:color]
  meshes = Makie.@lift let
    geoms = elements($mesh)
    bounds = boundary.(geoms)
    discretize.(bounds)
  end
  vizmany!(plot, meshes, color)
end

function segmentsof(topo, vert)
  p = first(vert)
  T = Unitful.numtype(Meshes.lentype(p))
  Dim = embeddim(p)
  nan = SVector(ntuple(i -> T(NaN), Dim))
  xs = map(p -> ustrip.(to(p)), vert)

  coords = map(elements(topo)) do e
    inds = indices(e)
    xs[collect(inds)]
  end

  reduce((x, y) -> [x; [nan]; y], coords)
end

function segmentsof(topo::GridTopology, vert)
  xs = map(p -> ustrip.(to(p)), vert)
  ip = first(isperiodic(topo))
  ip ? [xs; [first(xs)]] : xs
end
