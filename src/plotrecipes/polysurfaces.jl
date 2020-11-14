# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(polysurface::PolySurface)
  seriestype --> :path
  seriescolor --> :black
  label --> "polysurface"

  pchains = chains(polysurface)

  # plot outer chain
  @series begin
    first(pchains)
  end

  # plot inner chains
  for pchain in pchains[2:end]
    @series begin
      primary --> false
      pchain
    end
  end
end
