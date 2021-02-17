# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# recipe for base case
@recipe function f(segment::Segment)
  seriestype --> :scatterpath
  seriescolor --> :black
  primary --> false

  Tuple.(coordinates.(vertices(segment)))
end

# plot facets in a recursion
@recipe function f(polytope::Polytope)
  for f in facets(polytope)
    @series begin
      primary --> false
      f
    end
  end
end

@recipe function f(chain::Chain)
  seriestype --> :path
  seriescolor --> :black
  label --> "chain"

  xs = coordinates.(vertices(chain))

  isclosed(chain) && push!(xs, xs[begin])

  Tuple.(xs)
end

@recipe function f(polyarea::PolyArea)
  seriestype --> :path
  seriescolor --> :black
  label --> "polyarea"

  pchains = chains(polyarea)

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
