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
end
