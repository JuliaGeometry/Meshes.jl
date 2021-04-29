# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(point::Point{Dim,T}) where {Dim,T}
  seriestype --> :scatter
  seriescolor --> :black
  legend --> false
  [Tuple(coordinates(point))]
end

@recipe function f(points::AbstractVector{Point{Dim,T}}) where {Dim,T}
  seriestype --> :scatter
  seriescolor --> :black
  legend --> false
  Tuple.(coordinates.(points))
end
