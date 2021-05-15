# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(pset::PointSet{Dim,T}) where {Dim,T}
  pts = [pset[i] for i in 1:nelements(pset)]

  seriestype --> :scatter
  seriescolor --> :black
  legend --> false

  if Dim == 1
    @series begin
      first.(coordinates.(pts)), fill(0, nelements(pset))
    end
  else
    aspect_ratio --> :equal
    @series begin
      @. Tuple(coordinates(pts))
    end
  end
end

@recipe function f(pset::PointSet{Dim,T}, data::AbstractVector) where {Dim,T}
  pts = [pset[i] for i in 1:nelements(pset)]

  seriestype --> :scatter

  if Dim == 1
    first.(coordinates.(pts)), data
  else
    aspect_ratio --> :equal
    if eltype(data) <: Number
      marker_z --> data
      colorbar --> true
    else # categorical array
      group --> data
    end
    @. Tuple(coordinates(pts))
  end
end
