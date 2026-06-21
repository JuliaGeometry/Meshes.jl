# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# update plot with processed colors
function colorant!(plot)
  Makie.map!(plot, [:color, :alpha, :colormap, :colorrange], :colorant) do color, alpha, colorscheme, colorrange
    if color isa AbstractVector
      colorfy(color; alpha, colorscheme, colorrange)
    else
      colorfy([color]; alpha, colorscheme, colorrange) |> first
    end
  end
end

asray(v::Vec) = Ray(Point(zero(v)...), v)

function asmakie(v::Vec)
  N = length(v)
  ℒ = eltype(v)
  T = Unitful.numtype(ℒ)
  x = Tuple(ustrip.(v))
  Makie.Vec{N,T}(x)
end

function asmakie(p::Point)
  N = embeddim(p)
  ℒ = Meshes.lentype(p)
  T = Unitful.numtype(ℒ)
  x = Tuple(ustrip.(to(p)))
  Makie.Point{N,T}(x)
end

function asmakie(poly::Polygon)
  rs = rings(poly)
  outer = map(asmakie, eachvertex(rs[1]))
  if hasholes(poly)
    inners = map(i -> map(asmakie, eachvertex(rs[i])), 2:length(rs))
    Makie.Polygon(outer, inners)
  else
    Makie.Polygon(outer)
  end
end
