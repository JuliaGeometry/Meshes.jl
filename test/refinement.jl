@testitem "TriRefinement" setup = [Setup] begin
  grid = cartgrid(3, 3)
  ref1 = refine(grid, TriRefinement())
  ref2 = refine(ref1, TriRefinement())

  if visualtests
    fig = Mke.Figure(size=(900, 300))
    viz(fig[1, 1], grid, showsegments=true)
    viz(fig[1, 2], ref1, showsegments=true)
    viz(fig[1, 3], ref2, showsegments=true)
    @test_reference "data/trirefine-$T.png" fig
  end

  # CRS propagation
  grid = CartesianGrid((3, 3), merc(0, 0), (T(1), T(1)))
  ref = refine(grid, TriRefinement())
  @test crs(ref) === crs(grid)

  # predicate
  points = cart.([(0, 0), (4, 0), (8, 0), (3, 1), (5, 1), (2, 2), (4, 2), (6, 2), (4, 4)])
  connec = connect.([(1, 2, 6), (2, 3, 8), (6, 8, 9), (2, 5, 4), (4, 5, 7), (4, 7, 6), (5, 8, 7)])
  mesh = SimpleMesh(points, connec)
  ref = refine(mesh, TriRefinement(e -> measure(e) > T(1) * u"m^2"))
  @test nelements(ref) == 13
  @test nvertices(ref) == 12
  ref = refine(mesh, TriRefinement(e -> measure(e) ≤ T(1) * u"m^2"))
  @test nelements(ref) == 15
  @test nvertices(ref) == 13

  # latlon
  points = latlon.([(0, 0), (0, 4), (0, 8), (1, 3), (1, 5), (2, 2), (2, 4), (2, 6), (4, 4)])
  connec = connect.([(1, 2, 6), (2, 3, 8), (6, 8, 9), (2, 5, 4), (4, 5, 7), (4, 7, 6), (5, 8, 7)])
  mesh = SimpleMesh(points, connec)
  ref = refine(mesh, TriRefinement())
  @test nelements(ref) == 21
  @test nvertices(ref) == 16
end

@testitem "QuadRefinement" setup = [Setup] begin
  points = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.25, 0.25), (0.75, 0.25), (0.5, 0.75)])
  connec = connect.([(1, 2, 6, 5), (1, 5, 7, 3), (2, 4, 7, 6), (3, 7, 4)])
  mesh = SimpleMesh(points, connec)
  ref1 = refine(mesh, QuadRefinement())
  ref2 = refine(ref1, QuadRefinement())
  ref3 = refine(ref2, QuadRefinement())

  if visualtests
    fig = Mke.Figure(size=(900, 300))
    viz(fig[1, 1], ref1, showsegments=true)
    viz(fig[1, 2], ref2, showsegments=true)
    viz(fig[1, 3], ref3, showsegments=true)
    @test_reference "data/quadrefine-$T.png" fig
  end

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

@testitem "CatmullClark" setup = [Setup] begin
  points = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
  mesh = SimpleMesh(points, connec)
  ref1 = refine(mesh, CatmullClarkRefinement())
  ref2 = refine(ref1, CatmullClarkRefinement())
  ref3 = refine(ref2, CatmullClarkRefinement())

  if visualtests
    fig = Mke.Figure(size=(900, 300))
    viz(fig[1, 1], ref1, showsegments=true)
    viz(fig[1, 2], ref2, showsegments=true)
    viz(fig[1, 3], ref3, showsegments=true)
    @test_reference "data/catmullclark-1-$T.png" fig
  end

  points = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.25, 0.25), (0.75, 0.25), (0.5, 0.75)])
  connec = connect.([(1, 2, 6, 5), (1, 5, 7, 3), (2, 4, 7, 6), (3, 7, 4)])
  mesh = SimpleMesh(points, connec)
  ref1 = refine(mesh, CatmullClarkRefinement())
  ref2 = refine(ref1, CatmullClarkRefinement())
  ref3 = refine(ref2, CatmullClarkRefinement())

  if visualtests
    fig = Mke.Figure(size=(900, 300))
    viz(fig[1, 1], ref1, showsegments=true)
    viz(fig[1, 2], ref2, showsegments=true)
    viz(fig[1, 3], ref3, showsegments=true)
    @test_reference "data/catmullclark-2-$T.png" fig
  end

  points = cart.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
  connec = connect.([(1, 4, 3, 2), (5, 6, 7, 8), (1, 2, 6, 5), (3, 4, 8, 7), (1, 5, 8, 4), (2, 3, 7, 6)])
  mesh = SimpleMesh(points, connec)
  ref1 = refine(mesh, CatmullClarkRefinement())
  ref2 = refine(ref1, CatmullClarkRefinement())
  ref3 = refine(ref2, CatmullClarkRefinement())

  if visualtests
    fig = Mke.Figure(size=(900, 300))
    viz(fig[1, 1], ref1, showsegments=true)
    viz(fig[1, 2], ref2, showsegments=true)
    viz(fig[1, 3], ref3, showsegments=true)
    @test_reference "data/catmullclark-3-$T.png" fig
  end

  # CRS propagation
  points = merc.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
  mesh = SimpleMesh(points, connec)
  ref = refine(mesh, CatmullClarkRefinement())
  @test crs(ref) === crs(mesh)
end

@testitem "TriSubdivision" setup = [Setup] begin
  points = cart.([(-1, -1, -1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)])
  connec = connect.([(1, 2, 3), (3, 2, 4), (4, 2, 1), (1, 3, 4)])
  mesh = SimpleMesh(points, connec)
  ref1 = refine(mesh, TriSubdivision())
  ref2 = refine(ref1, TriSubdivision())
  ref3 = refine(ref2, TriSubdivision())

  if visualtests
    fig = Mke.Figure(size=(900, 300))
    viz(fig[1, 1], ref1, showsegments=true)
    viz(fig[1, 2], ref2, showsegments=true)
    viz(fig[1, 3], ref3, showsegments=true)
    @test_reference "data/trisubdivision-$T.png" fig
  end
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
