# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(pset::PointSet{Dim,T}) where {Dim,T}
  X = coordinates(pset, 1:nelements(pset))

  seriestype --> :scatter
  seriescolor --> :black
  legend --> false

  if Dim == 1
    @series begin
      X[1,:], fill(0, nelements(pset))
    end
  elseif Dim == 2
    aspect_ratio --> :equal
    @series begin
      X[1,:], X[2,:]
    end
  elseif Dim == 3
    aspect_ratio --> :equal
    @series begin
      X[1,:], X[2,:], X[3,:]
    end
  else
    @error "cannot plot in more than 3 dimensions"
  end
end

@recipe function f(pset::PointSet{Dim,T}, data::AbstractVector) where {Dim,T}
  X = coordinates(pset, 1:nelements(pset))

  seriestype --> :scatter

  if Dim == 1
    X[1,:], data
  elseif Dim == 2
    aspect_ratio --> :equal
    marker_z --> data
    colorbar --> true
    X[1,:], X[2,:]
  elseif Dim == 3
    aspect_ratio --> :equal
    marker_z --> data
    colorbar --> true
    X[1,:], X[2,:], X[3,:]
  else
    @error "cannot plot in more than 3 dimensions"
  end
end
