# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, p::Point{Dim,T}, g::Geometry{Dim,T}) where {Dim,T}
  if p ∈ g
    @IT PertainingPoint p f
  else
    @IT NoIntersection nothing f
  end
end
