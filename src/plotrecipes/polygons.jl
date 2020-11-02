# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(polygon::Polygon)
  seriestype --> :scatterpath
  seriescolor --> :black
  label --> "polygon"

  outer, inners = rings(polygon)

  # plot outer ring
  @series begin
    Tuple.(coordinates.(vertices(outer)))
  end

  # plot inner rings
  for inner in inners
    @series begin
      primary --> false
      Tuple.(coordinates.(vertices(inner)))
    end
  end
end
