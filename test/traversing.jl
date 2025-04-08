@testitem "Traversing" setup = [Setup] begin
  # LinearPath
  grid = cartgrid(100, 100)
  path = LinearPath()
  @test traverse(grid, path) == 1:10000

  # RandomPath
  grid = cartgrid(100, 100)
  path = RandomPath(StableRNG(123))
  @test all(1 .≤ traverse(grid, path) .≤ 10000)
  grid = cartgrid(3, 3)
  path = RandomPath(StableRNG(123))
  @test traverse(grid, path) == [4, 7, 2, 1, 3, 8, 5, 6, 9]

  # SourcePath
  grid = cartgrid(3, 3)
  pset = PointSet(centroid.(grid))
  for sdomain in [grid, pset]
    spath = SourcePath([1, 9])
    @test traverse(sdomain, spath) == [1, 9, 2, 4, 6, 8, 5, 3, 7]
    spath = SourcePath([1])
    @test traverse(sdomain, spath) == [1, 2, 4, 5, 3, 7, 6, 8, 9]
  end

  # ShiftedPath
  grid = cartgrid(3, 3)
  path = LinearPath()
  for offset in [0, 1, -1]
    spath = ShiftedPath(path, offset)
    t = traverse(grid, path)
    st = traverse(grid, spath)
    @test length(st) == 9
    @test collect(st) == circshift(t, -offset)
  end

  # MultiGridPath
  path = MultiGridPath()
  grid = cartgrid(3, 3)
  @test traverse(grid, path) == [1, 3, 7, 9, 2, 4, 5, 6, 8]
  grid = cartgrid(3, 4)
  @test traverse(grid, path) == [1, 3, 10, 12, 2, 7, 8, 9, 4, 5, 6, 11]
  grid = CartesianGrid(3, 3, 2)
  @test traverse(grid, path) == [1, 3, 7, 9, 10, 12, 16, 18, 2, 4, 5, 6, 8, 11, 13, 14, 15, 17]
  grid = RectilinearGrid(T.(0:3), T.(0:3))
  @test traverse(grid, path) == [1, 3, 7, 9, 2, 4, 5, 6, 8]
  grid = RectilinearGrid(T.(0:0.5:2), T.(0:0.5:2))
  @test traverse(grid, path) == [1, 4, 13, 16, 3, 9, 11, 2, 5, 6, 7, 8, 10, 12, 14, 15]
  cgrid = cartgrid(4, 4)
  rgrid = RectilinearGrid(T.(0:4), T.(0:4))
  @test traverse(cgrid, path) == traverse(rgrid, path)
  grid = cartgrid(3, 4)
  vgrid = view(grid, 3:10)
  @test traverse(vgrid, path) == [3, 10, 7, 8, 9, 4, 5, 6]

  # visual tests
  if visualtests
    paths = [LinearPath(), RandomPath(StableRNG(123)), ShiftedPath(LinearPath(), 10), SourcePath(1:3), MultiGridPath()]

    fnames = ["linear-path", "random-path", "shifted-path", "source-path", "multi-grid-path"]

    for (path, fname) in zip(paths, fnames)
      for n in (6, 7)
        agrid = cartgrid(n, n)
        pinds = collect(traverse(agrid, path))
        pgrid = view(agrid, pinds)
        fig = viz(pgrid, color=1:nelements(pgrid))
        @test_reference "data/$fname-$(n)x$(n).png" fig
      end
    end
  end
end
