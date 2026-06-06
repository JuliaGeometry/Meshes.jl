# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# update plot with processed colors
function colorant!(plot)
  Makie.map!(plot, [:color, :alpha, :colormap, :colorrange], :colorant) do color, alphas, colorscheme, colorrange
    if color isa AbstractVector
      colorfy(color; alphas, colorscheme, colorrange)
    else
      colorfy([color]; alphas, colorscheme, colorrange) |> first
    end
  end
end

asray(v::Vec) = Ray(Point(zero(v)...), v)

asmakie(v::Vec) = Makie.Vec{length(v),numtype(eltype(v))}(Tuple(ustrip.(v)))

asmakie(p::Point) = Makie.Point{embeddim(p),numtype(lentype(p))}(Tuple(ustrip.(to(p))))

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
