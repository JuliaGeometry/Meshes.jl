@testset "Paths" begin
  grid = CartesianGrid{T}(100, 100)

  for path in [LinearPath(), RandomPath(), ShiftedPath(LinearPath(), 0), SourcePath(1:3)]
    p = traverse(grid, path)
    @test length(p) == 100 * 100
  end

  @testset "LinearPath" begin
    p = traverse(grid, LinearPath())
    @test p == 1:(100 * 100)
  end

  @testset "RandomPath" begin
    p = traverse(grid, RandomPath())
    @test all(1 .≤ collect(p) .≤ 100 * 100)
  end

  @testset "SourcePath" begin
    grid = CartesianGrid{T}(3, 3)
    pset = PointSet(centroid.(grid))

    for sdomain in [grid, pset]
      p = traverse(sdomain, SourcePath([1, 9]))
      @test collect(p) == [1, 9, 2, 4, 6, 8, 5, 3, 7]

      p = traverse(sdomain, SourcePath([1]))
      @test collect(p) == [1, 2, 4, 5, 3, 7, 6, 8, 9]
    end
  end

  @testset "ShiftedPath" begin
    grid = CartesianGrid{T}(3, 3)
    path = LinearPath()
    for offset in [0, 1, -1]
      spath = ShiftedPath(path, offset)
      p = traverse(grid, path)
      sp = traverse(grid, spath)
      @test length(sp) == 9
      @test collect(sp) == circshift(p, -offset)
    end
  end

  @testset "MultiGridPath" begin
    path = MultiGridPath()

    grid = CartesianGrid{T}(3, 3)
    @test traverse(grid, path) == [1, 3, 7, 9, 2, 4, 5, 6, 8]

    grid = CartesianGrid{T}(3, 4)
    @test traverse(grid, path) == [1, 3, 10, 12, 2, 7, 8, 9, 4, 5, 6, 11]

    grid = CartesianGrid(3, 3, 2)
    @test traverse(grid, path) == [1, 3, 7, 9, 10, 12, 16, 18, 2, 4, 5, 6, 8, 11, 13, 14, 15, 17]

    grid = RectilinearGrid(T.(0:3), T.(0:3))
    @test traverse(grid, path) == [1, 3, 7, 9, 2, 4, 5, 6, 8]

    grid = RectilinearGrid(T.(0:0.5:2), T.(0:0.5:2))
    @test traverse(grid, path) == [1, 4, 13, 16, 3, 9, 11, 2, 5, 6, 7, 8, 10, 12, 14, 15]

    cgrid = CartesianGrid{T}(4, 4)
    rgrid = RectilinearGrid(T.(0:4), T.(0:4))
    @test traverse(cgrid, path) == traverse(rgrid, path)

    if visualtests
      grid = CartesianGrid{T}(7, 7)
      elems = [grid[i] for i in traverse(grid, path)]
      fig = viz(elems, color=1:length(elems))
      @test_reference "data/multi-grid-path-7x7.png" fig

      grid = CartesianGrid{T}(6, 6)
      elems = [grid[i] for i in traverse(grid, path)]
      fig = viz(elems, color=1:length(elems))
      @test_reference "data/multi-grid-path-6x6.png" fig
    end
  end
end
