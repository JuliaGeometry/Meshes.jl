# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{Polygon{2}}})
  poly = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  segmentsize = plot[:segmentsize]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # visualize as built-in poly
  gbpoly = Makie.@lift asgbpoly($poly)
  Makie.poly!(plot, gbpoly, color=colorant)

  if showfacets[]
    tup = Makie.@lift polysegments($poly)
    x, y = Makie.@lift($tup[1]), Makie.@lift($tup[2])
    Makie.lines!(plot, x, y, color=facetcolor, linewidth=segmentsize)
  end
end

# converts Meshes.Polygon to GeometryBasics.Polygon
function asgbpoly(poly)
  rs = rings(poly)
  outer = first(rs)
  opts = [asgbpoint(p) for p in vertices(outer)]
  if hasholes(poly)
    ipts = map(i -> [asgbpoint(p) for p in vertices(rs[i])], 2:length(rs))
    Makie.Polygon(opts, ipts)
  else
    Makie.Polygon(opts)
  end
end

# converts Meshes.Point to GeometryBasics.Point
asgbpoint(p::Point{Dim,T}) where {Dim,T} = Makie.Point{Dim,T}(getcoords(p))

# returns the line segments of the 2D Polygon
function polysegments(poly)
  rs = rings(poly)
  outer = first(rs)
  opts = vertices(outer)
  coords = [getcoords(p) for p in vertices(outer)]
  push!(coords, getcoords(first(opts)))
  if hasholes(poly)
    for i in 2:length(rs)
      ipts = vertices(rs[i])
      push!(coords, (NaN, NaN))
      foreach(p -> push!(coords, getcoords(p)), ipts)
      push!(coords, getcoords(first(ipts)))
    end
  end
  x = getindex.(coords, 1)
  y = getindex.(coords, 2)
  x, y
end

# extracts point coordinates as tuple
getcoords(p) = Tuple(coordinates(p))
