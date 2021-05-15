# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(pset::PointSet{Dim,T}) where {Dim,T}
  coords = Tuple.(coordinates.(pset))

  seriestype --> :scatter
  seriescolor --> :black
  legend --> false

  if Dim == 1
    @series begin
      first.(coords), fill(0, length(coords))
    end
  else
    aspect_ratio --> :equal
    @series begin
      coords
    end
  end
end

@recipe function f(pset::PointSet{Dim,T}, data::AbstractVector) where {Dim,T}
  coords = Tuple.(coordinates.(pset))

  seriestype --> :scatter

  if Dim == 1
    first.(coords), data
  else
    aspect_ratio --> :equal
    marker_z --> data
    colorbar --> true
    coords
  end
end
