@testitem "TriRefinement" setup = [Setup] begin
  # CRS propagation
  grid = CartesianGrid((3, 3), merc(0, 0), (T(1), T(1)))
  rgrid = refine(grid, TriRefinement())
  @test crs(rgrid) === crs(grid)

  # quadrangles are divided into 4 triangles
  grid = cartgrid(10, 10)
  rgrid = refine(grid, TriRefinement())
  @test nelements(rgrid) == 4 * nelements(grid)

  # predicate
  points = cart.([(0, 0), (4, 0), (8, 0), (3, 1), (5, 1), (2, 2), (4, 2), (6, 2), (4, 4)])
  connec = connect.([(1, 2, 6), (2, 3, 8), (6, 8, 9), (2, 5, 4), (4, 5, 7), (4, 7, 6), (5, 8, 7)])
  mesh = SimpleMesh(points, connec)
  rmesh = refine(mesh, TriRefinement(e -> measure(e) > T(1) * u"m^2"))
  @test nelements(rmesh) == 13
  @test nvertices(rmesh) == 12
  rmesh = refine(mesh, TriRefinement(e -> measure(e) ≤ T(1) * u"m^2"))
  @test nelements(rmesh) == 15
  @test nvertices(rmesh) == 13

  # latlon
  points = latlon.([(0, 0), (0, 4), (0, 8), (1, 3), (1, 5), (2, 2), (2, 4), (2, 6), (4, 4)])
  connec = connect.([(1, 2, 6), (2, 3, 8), (6, 8, 9), (2, 5, 4), (4, 5, 7), (4, 7, 6), (5, 8, 7)])
  mesh = SimpleMesh(points, connec)
  rmesh = refine(mesh, TriRefinement())
  @test nelements(rmesh) == 21
  @test nvertices(rmesh) == 16
end

@testitem "QuadRefinement" setup = [Setup] begin
  # CRS propagation
  points = merc.([(0, 0), (1, 0), (0, 1), (1, 1), (0.25, 0.25), (0.75, 0.25), (0.5, 0.75)])
  connec = connect.([(1, 2, 6, 5), (1, 5, 7, 3), (2, 4, 7, 6), (3, 7, 4)])
  mesh = SimpleMesh(points, connec)
  ref = refine(mesh, QuadRefinement())
  @test crs(ref) === crs(mesh)

  # latlon
  points = latlon.([(0, 0), (0, 1), (1, 0), (1, 1), (0.25, 0.25), (0.25, 0.75), (0.75, 0.5)])
  connec = connect.([(1, 2, 6, 5), (1, 5, 7, 3), (2, 4, 7, 6), (3, 7, 4)])
  mesh = SimpleMesh(points, connec)
  ref = refine(mesh, QuadRefinement())
  @test nelements(ref) == 15
  @test nvertices(ref) == 22
end

@testitem "TriSubdivision" setup = [Setup] begin
  # CRS propagation
  points = merc.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
  mesh = SimpleMesh(points, connec)
  ref = refine(mesh, TriSubdivision())
  @test crs(ref) === crs(mesh)
end

@testitem "CatmullClark" setup = [Setup] begin
  # CRS propagation
  points = merc.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
  mesh = SimpleMesh(points, connec)
  ref = refine(mesh, CatmullClarkRefinement())
  @test crs(ref) === crs(mesh)
end

@testitem "RegularRefinement" setup = [Setup] begin
  # 2D grids
  grid = CartesianGrid(cart(0, 0), cart(10, 10), dims=(10, 10))
  tgrid = CartesianGrid(cart(0, 0), cart(10, 10), dims=(20, 20))
  @test refine(grid, RegularRefinement(2)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  trgrid = convert(RectilinearGrid, tgrid)
  @test refine(rgrid, RegularRefinement(2)) == trgrid
  sgrid = convert(StructuredGrid, grid)
  tsgrid = convert(StructuredGrid, tgrid)
  @test refine(sgrid, RegularRefinement(2)) == tsgrid
  tfgrid = TransformedGrid(grid, Identity())
  @test refine(tfgrid, RegularRefinement(2)) == refine(grid, RegularRefinement(2))

  # 3D grids
  grid = cartgrid(3, 3, 3)
  tgrid = CartesianGrid(minimum(grid), maximum(grid), dims=(6, 6, 6))
  @test refine(grid, RegularRefinement(2)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  trgrid = convert(RectilinearGrid, tgrid)
  @test refine(rgrid, RegularRefinement(2)) == trgrid
  sgrid = convert(StructuredGrid, grid)
  tsgrid = convert(StructuredGrid, tgrid)
  @test refine(sgrid, RegularRefinement(2)) == tsgrid
  tfgrid = TransformedGrid(grid, Identity())
  @test refine(tfgrid, RegularRefinement(2)) == refine(grid, RegularRefinement(2))
end

@testitem "MaxLengthRefinement" setup = [Setup] begin
  # 2D grids
  grid = CartesianGrid(cart(0, 0), cart(10, 10), dims=(10, 10))
  tgrid = CartesianGrid(cart(0, 0), cart(10, 10), dims=(20, 20))
  method = MaxLengthRefinement(T(0.5) * u"m")
  @test refine(grid, method) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  trgrid = convert(RectilinearGrid, tgrid)
  @test refine(rgrid, method) == trgrid
  sgrid = convert(StructuredGrid, grid)
  tsgrid = convert(StructuredGrid, tgrid)
  @test refine(sgrid, method) == tsgrid

  # general meshes
  mesh = convert(SimpleMesh, cartgrid(2, 2))
  rmesh = refine(mesh, MaxLengthRefinement(T(0.5) * u"m"))
  @test all(e -> perimeter(e) / 3 ≤ T(0.5) * u"m", rmesh)
end

@testitem "Refine" setup = [Setup] begin
  # 2D grids
  grid = cartgrid(10, 10)
  rgrid = refine(grid)
  @test size(rgrid) == (20, 20)

  # general meshes
  grid = cartgrid(10, 10)
  mesh = topoconvert(SimpleTopology, grid)
  rmesh = refine(mesh)
  @test eltype(rmesh) <: Triangle
  @test nelements(rmesh) == 400
end
