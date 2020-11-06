# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(polygon::Polygon)
  seriestype --> :path
  seriescolor --> :black
  label --> "polygon"

  outer, inners = rings(polygon)

  # plot outer ring
  @series begin
    outer
  end

  # plot inner rings
  for inner in inners
    @series begin
      primary --> false
      inner
    end
  end
end
