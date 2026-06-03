# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plot!(plot::Viz{<:Tuple{SubDomain}}) = vizsubdom!(plot)

function vizsubdom!(plot)
  subdom = plot[:object][]
  M = manifold(subdom)
  pdim = paramdim(subdom)
  edim = embeddim(subdom)
  vizsubdom!(plot, M, Val(pdim), Val(edim))
end

# ---------------
# IMPLEMENTATION
# ---------------

function vizsubdom!(plot, ::Type{<:🌐}, pdim::Val, edim::Val)
  vizsubdom!(plot, 𝔼, pdim, edim)
end

function vizsubdom!(plot, ::Type{<:𝔼}, ::Val, ::Val)
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  showsegments = plot[:showsegments]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]
  showpoints = plot[:showpoints]
  pointmarker = plot[:pointmarker]
  pointcolor = plot[:pointcolor]
  pointsize = plot[:pointsize]

  # construct the geometry set
  Makie.map!(plot, [:object], :gset) do subdom
    GeometrySet(collect(subdom))
  end

  # forward attributes
  viz!(
    plot,
    plot[:gset];
    color,
    alpha,
    colormap,
    colorrange,
    showsegments,
    segmentcolor,
    segmentsize,
    showpoints,
    pointmarker,
    pointcolor,
    pointsize
  )
end

const SubCartesianGrid{M,CRS} = SubDomain{M,CRS,<:CartesianGrid}

function vizsubdom!(plot::Viz{<:Tuple{SubCartesianGrid}}, ::Type{<:𝔼}, ::Val, ::Val)

  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

  # retrieve grid paramaters
  Makie.map!(plot, [:object], :gparams) do subgrid
    grid = parent(subgrid)
    dim = embeddim(grid)
    sp = ustrip.(spacing(grid))

    # coordinates of markers
    coords = map(subgrid) do e
      ustrip.(to(centroid(e))) .+ sp ./ 2
    end

    # rectangle markers
    marker = Makie.Rect{dim}(-1 .* sp, sp)

    # enable shading in 3D
    shading = dim == 3

    coords, marker, shading
  end

  # unpack observable parameters
  Makie.map!(plot, [:gparams], [:coords, :marker, :shading]) do gparams
    (gparams[1], gparams[2], gparams[3])
  end

  # all geometries are equal, use mesh scatter
  Makie.meshscatter!(
    plot,
    plot[:coords],
    marker=plot[:marker],
    markersize=1,
    color=plot[:colorant],
    shading=plot[:shading]
  )
end
