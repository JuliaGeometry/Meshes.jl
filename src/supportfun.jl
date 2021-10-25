# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    supportfun(geometry, direction)

Support function of `geometry` for given `direction`.

## References

* Gilbert, E., Johnson, D., Keerthi, S. 1988. [A fast
  Procedure for Computing the Distance Between Complex
  Objects in Three-Dimensional Space]
  (https://ieeexplore.ieee.org/document/2083)
"""
function supportfun(g::Geometry, d::Vec)
  v = vertices(g)
  c = centroid(g)
  _, i = findmax(vᵢ -> (vᵢ - c) ⋅ d, v)
  v[i]
end

function supportfun(b::Ball, d::Vec)
  r = (radius(b) / norm(d)) * d
  center(b) + r
end
