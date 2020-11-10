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
