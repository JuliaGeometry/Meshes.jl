@testset "Refinement" begin
  @testset "TriRefinement" begin
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

    # datum propagation
    c = Cartesian{WGS84Latest}(T(0), T(0))
    grid = CartesianGrid((3, 3), Point(c), (T(1), T(1)))
    ref = refine(grid, TriRefinement())
    @test datum(Meshes.crs(ref)) === WGS84Latest
  end

  @testset "QuadRefinement" begin
    points = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.25, 0.25), (0.75, 0.25), (0.5, 0.75)])
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

    # datum propagation
    tuples = [T.((0, 0)), T.((1, 0)), T.((0, 1)), T.((1, 1)), T.((0.25, 0.25)), T.((0.75, 0.25)), T.((0.5, 0.75))]
    points = Point.(Cartesian{WGS84Latest}.(tuples))
    connec = connect.([(1, 2, 6, 5), (1, 5, 7, 3), (2, 4, 7, 6), (3, 7, 4)])
    mesh = SimpleMesh(points, connec)
    ref = refine(mesh, QuadRefinement())
    @test datum(Meshes.crs(ref)) === WGS84Latest
  end

  @testset "RegularRefinement" begin
    # 2D grids
    grid = CartesianGrid(point(0, 0), point(10, 10), dims=(10, 10))
    tgrid = CartesianGrid(point(0, 0), point(10, 10), dims=(20, 20))
    @test refine(grid, RegularRefinement(2)) == tgrid
    rgrid = convert(RectilinearGrid, grid)
    trgrid = convert(RectilinearGrid, tgrid)
    @test refine(rgrid, RegularRefinement(2)) == trgrid
    sgrid = convert(StructuredGrid, grid)
    tsgrid = convert(StructuredGrid, tgrid)
    @test refine(sgrid, RegularRefinement(2)) == tsgrid

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
  end

  @testset "CatmullClark" begin
    points = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
    mesh = SimpleMesh(points, connec)
    ref1 = refine(mesh, CatmullClark())
    ref2 = refine(ref1, CatmullClark())
    ref3 = refine(ref2, CatmullClark())

    if visualtests
      fig = Mke.Figure(size=(900, 300))
      viz(fig[1, 1], ref1, showsegments=true)
      viz(fig[1, 2], ref2, showsegments=true)
      viz(fig[1, 3], ref3, showsegments=true)
      @test_reference "data/catmullclark-1-$T.png" fig
    end

    points = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.25, 0.25), (0.75, 0.25), (0.5, 0.75)])
    connec = connect.([(1, 2, 6, 5), (1, 5, 7, 3), (2, 4, 7, 6), (3, 7, 4)])
    mesh = SimpleMesh(points, connec)
    ref1 = refine(mesh, CatmullClark())
    ref2 = refine(ref1, CatmullClark())
    ref3 = refine(ref2, CatmullClark())

    if visualtests
      fig = Mke.Figure(size=(900, 300))
      viz(fig[1, 1], ref1, showsegments=true)
      viz(fig[1, 2], ref2, showsegments=true)
      viz(fig[1, 3], ref3, showsegments=true)
      @test_reference "data/catmullclark-2-$T.png" fig
    end

    points = point.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
    connec = connect.([(1, 4, 3, 2), (5, 6, 7, 8), (1, 2, 6, 5), (3, 4, 8, 7), (1, 5, 8, 4), (2, 3, 7, 6)])
    mesh = SimpleMesh(points, connec)
    ref1 = refine(mesh, CatmullClark())
    ref2 = refine(ref1, CatmullClark())
    ref3 = refine(ref2, CatmullClark())

    if visualtests
      fig = Mke.Figure(size=(900, 300))
      viz(fig[1, 1], ref1, showsegments=true)
      viz(fig[1, 2], ref2, showsegments=true)
      viz(fig[1, 3], ref3, showsegments=true)
      @test_reference "data/catmullclark-3-$T.png" fig
    end

    # datum propagation
    tuples = [T.((0, 0)), T.((1, 0)), T.((0, 1)), T.((1, 1)), T.((0.5, 0.5))]
    points = Point.(Cartesian{WGS84Latest}.(tuples))
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
    mesh = SimpleMesh(points, connec)
    ref = refine(mesh, CatmullClark())
    @test datum(Meshes.crs(ref)) === WGS84Latest
  end

  @testset "TriSubdivision" begin
    points = point.([(-1, -1, -1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)])
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
end
