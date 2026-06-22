# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using PrecompileTools

@setup_workload begin
  # vectors of geometries in GIS
  # TODO: include LatLon after https://github.com/JuliaGeometry/Meshes.jl/issues/1367
  ctype = [Cartesian2D]
  gtype = [Point, Rope, Ring, PolyArea]
  geoms = [rand(G, 3, crs=C) for C in ctype for G in gtype]
  multi = [Multi(geoms[i]) for i in eachindex(geoms)]

  # regular grids (e.g., "raster" images)
  ctype = [Cartesian2D, Cartesian3D, LatLon]
  grids = map(ctype) do C
    N = CoordRefSystems.ncoords(C)
    xyz₁ = ntuple(_ -> 0.0, N)
    xyz₂ = ntuple(_ -> 1.0, N)
    dims = ntuple(_ -> 2, N)
    RegularGrid(Point(C(xyz₁...)), Point(C(xyz₂...)), dims=dims)
  end

  # point clouds in 3D Euclidean space
  ctype = [Cartesian3D]
  clouds = [rand(Point, 3, crs=C) for C in ctype]

  # surface meshes in 3D Euclidean space
  ctype = [Cartesian3D]
  meshes = map(ctype) do C
    points = rand(Point, 3, crs=C)
    connec = connect.([(1, 2, 3)])
    SimpleMesh(points, connec)
  end

  @compile_workload begin
    for i in eachindex(geoms)
      viz(geoms[i], color=1:length(geoms[i]))
      viz(multi[i])
    end

    for i in eachindex(grids)
      viz(grids[i], color=1:nelements(grids[i]))
    end

    for i in eachindex(clouds)
      viz(clouds[i], color=1:length(clouds[i]))
    end

    for i in eachindex(meshes)
      viz(meshes[i], color=1:nelements(meshes[i]))
    end
  end
end
