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
@recipe function f(face::Face)
  for fa in facets(face)
    @series begin
      primary --> false
      fa
    end
  end
end
