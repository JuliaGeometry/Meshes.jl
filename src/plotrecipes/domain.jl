# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(domain::Domain)
  linecolor --> :black
  collect(domain)
end

@recipe function f(domain::Domain, data::AbstractVector)
  points = [centroid(domain, i) for i in 1:nelements(domain)]
  PointSet(points), data
end
