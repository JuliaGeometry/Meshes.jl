# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(polysurface::PolySurface)
  seriestype --> :path
  seriescolor --> :black
  label --> "polysurface"

  outer, inners = chains(polysurface)

  # plot outer chain
  @series begin
    outer
  end

  # plot inner chains
  for inner in inners
    @series begin
      primary --> false
      inner
    end
  end
end
