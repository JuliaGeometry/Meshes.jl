# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(domain::Domain)
  PointSet(centroid.(domain))
end

@recipe function f(domain::Domain, data::AbstractVector)
  PointSet(centroid.(domain)), data
end
