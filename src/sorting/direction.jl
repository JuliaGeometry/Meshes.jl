# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DirectionSort(direction)

Sort geometric objects along a given `direction` vector.
"""
struct DirectionSort{V<:Vec} <: SortingMethod
  direction::V
end

DirectionSort(direction::Tuple) = DirectionSort(Vec(direction))

function sortinds(domain::Domain, method::DirectionSort)
  v = method.direction
  t = map(1:nelements(domain)) do i
    c = centroid(domain, i)
    u = to(c)
    (u ⋅ v) / (v ⋅ v)
  end
  sortperm(t)
end
