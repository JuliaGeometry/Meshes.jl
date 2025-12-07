@testitem "TriRefinement" setup = [Setup] begin
  # CRS propagation
  grid = CartesianGrid(merc(0, 0), merc(3, 3))
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
  grid = cartgrid(10, 10, 10)
  tgrid = CartesianGrid(minimum(grid), maximum(grid), dims=(20, 20, 20))
  @test refine(grid, RegularRefinement(2)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  trgrid = convert(RectilinearGrid, tgrid)
  @test refine(rgrid, RegularRefinement(2)) == trgrid
  sgrid = convert(StructuredGrid, grid)
  tsgrid = convert(StructuredGrid, tgrid)
  @test refine(sgrid, RegularRefinement(2)) == tsgrid
  tfgrid = TransformedGrid(grid, Identity())
  @test refine(tfgrid, RegularRefinement(2)) == refine(grid, RegularRefinement(2))

  # non-multiple dimensions (2D grids)
  grid = CartesianGrid(cart(0, 0), cart(13, 17), dims=(13, 17))
  tgrid = CartesianGrid(cart(0, 0), cart(13, 17), dims=(65, 51))
  @test refine(grid, RegularRefinement(5, 3)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  @test size(refine(rgrid, RegularRefinement(5, 3))) == (65, 51)
  sgrid = convert(StructuredGrid, grid)
  @test size(refine(sgrid, RegularRefinement(5, 3))) == (65, 51)
  tfgrid = TransformedGrid(grid, Identity())
  @test size(refine(tfgrid, RegularRefinement(5, 3))) == (65, 51)

  # non-multiple dimensions (3D grids)
  grid = CartesianGrid(cart(0, 0, 0), cart(13, 17, 23), dims=(13, 17, 23))
  tgrid = CartesianGrid(cart(0, 0, 0), cart(13, 17, 23), dims=(91, 85, 69))
  @test refine(grid, RegularRefinement(7, 5, 3)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  @test size(refine(rgrid, RegularRefinement(7, 5, 3))) == (91, 85, 69)
  sgrid = convert(StructuredGrid, grid)
  @test size(refine(sgrid, RegularRefinement(7, 5, 3))) == (91, 85, 69)
  tfgrid = TransformedGrid(grid, Identity())
  @test size(refine(tfgrid, RegularRefinement(7, 5, 3))) == (91, 85, 69)

  # preserve topology
  topo = GridTopology((50, 50), (true, false))
  ttopo = GridTopology((100, 100), (true, false))
  grid = CartesianGrid(cart(0, 0), T.((1, 1)), topo)
  @test topology(refine(grid, RegularRefinement(2))) == ttopo
  rgrid = convert(RectilinearGrid, grid)
  @test topology(refine(rgrid, RegularRefinement(2))) == ttopo
  sgrid = convert(StructuredGrid, grid)
  @test topology(refine(sgrid, RegularRefinement(2))) == ttopo
  tfgrid = TransformedGrid(grid, Identity())
  @test topology(refine(tfgrid, RegularRefinement(2))) == ttopo

  # large 2D grid
  grid = CartesianGrid(cart(0, 0), cart(16200, 8100), dims=(203, 203))
  tgrid = CartesianGrid(cart(0, 0), cart(16200, 8100), dims=(16240, 8120))
  @test refine(grid, RegularRefinement(80, 40)) == tgrid

  # LatLon grid
  grid = RegularGrid(latlon(0, 0), latlon(45, 45), dims=(5, 5))
  tgrid = RegularGrid(latlon(0, 0), latlon(45, 45), dims=(10, 10))
  @test refine(grid, RegularRefinement(2)) == tgrid
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
