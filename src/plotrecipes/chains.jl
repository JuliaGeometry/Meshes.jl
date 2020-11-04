# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(chain::Chain)
  seriestype --> :scatterpath
  seriescolor --> :black
  label --> "chain"

  xs = coordinates.(vertices(chain))

  Tuple.(xs)
end
