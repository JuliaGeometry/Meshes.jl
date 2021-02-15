# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ---------
# FALLBACK
# ---------

@recipe function f(mesh::Mesh)
  seriescolor --> :black
  legend --> false

  for elm in elements(mesh)
    @series begin
      primary --> false
      elm
    end
  end
end

# ---------------
# CARTESIAN GRID
# ---------------

@recipe function f(grid::CartesianGrid{Dim,T}) where {Dim,T}
  or = coordinates(minimum(grid))
  sp = spacing(grid)
  sz = size(grid)

  seriescolor --> :black
  legend --> false

  if Dim == 1
    @series begin
      seriestype --> :scatterpath
      marker --> :vline
      x = range(or[1]-sp[1]/2, step=sp[1], length=sz[1]+1)
      x, fill(zero(T), sz[1]+1)
    end
  elseif Dim == 2
    @series begin
      seriestype --> :path
      aspect_ratio --> :equal
      xs = range(or[1]-sp[1]/2, step=sp[1], length=sz[1]+1)
      ys = range(or[2]-sp[2]/2, step=sp[2], length=sz[2]+1)
      coords = []
      for x in xs
        push!(coords, (x, first(ys)))
        push!(coords, (x, last(ys)))
        push!(coords, (NaN, NaN))
      end
      for y in ys
        push!(coords, (first(xs), y))
        push!(coords, (last(xs), y))
        push!(coords, (NaN, NaN))
      end
      x = getindex.(coords, 1)
      y = getindex.(coords, 2)
      x, y
    end
  elseif Dim == 3
    @series begin
      seriestype --> :path
      aspect_ratio --> :equal
      xs = range(or[1]-sp[1]/2, step=sp[1], length=sz[1]+1)
      ys = range(or[2]-sp[2]/2, step=sp[2], length=sz[2]+1)
      zs = range(or[3]-sp[3]/2, step=sp[3], length=sz[3]+1)
      coords = []
      for y in ys, z in zs
        push!(coords, (first(xs), y, z))
        push!(coords, (last(xs), y, z))
        push!(coords, (NaN, NaN, NaN))
      end
      for x in xs, z in zs
        push!(coords, (x, first(ys), z))
        push!(coords, (x, last(ys), z))
        push!(coords, (NaN, NaN, NaN))
      end
      for x in xs, y in ys
        push!(coords, (x, y, first(zs)))
        push!(coords, (x, y, last(zs)))
        push!(coords, (NaN, NaN, NaN))
      end
      x = getindex.(coords, 1)
      y = getindex.(coords, 2)
      z = getindex.(coords, 3)
      x, y, z
    end
  else
    @error "cannot plot in more than 3 dimensions"
  end
end

@recipe function f(grid::CartesianGrid{Dim,T}, data::AbstractVector) where {Dim,T}
  or = coordinates(minimum(grid))
  sp = spacing(grid)
  sz = size(grid)
  Z  = reshape(data, sz)

  if Dim == 1
    seriestype --> :path
    x = range(or[1], step=sp[1], length=sz[1])
    x, Z
  elseif Dim == 2
    seriestype --> :heatmap
    aspect_ratio --> :equal
    colorbar --> true
    x = range(or[1], step=sp[1], length=sz[1])
    y = range(or[2], step=sp[2], length=sz[2])
    x, y, reverse(rotr90(Z), dims=2)
  elseif Dim == 3
    seriestype --> :heatmap
    aspect_ratio --> :equal
    colorbar --> true
    Z
  else
    @error "cannot plot in more than 3 dimensions"
  end
end
