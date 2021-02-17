# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(domain::Domain)
  PointSet(domain)
end

@recipe function f(domain::Domain, data::AbstractVector)
  PointSet(domain), data
end
