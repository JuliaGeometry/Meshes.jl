@testitem "Traversing" setup = [Setup] begin
  grid = cartgrid(100, 100)
  for path in [LinearPath(), RandomPath(), ShiftedPath(LinearPath(), 0), SourcePath(1:3)]
    p = traverse(grid, path)
    @test length(p) == 100 * 100
  end

  grid = cartgrid(100, 100)
  p = traverse(grid, LinearPath())
  @test p == 1:(100 * 100)

  grid = cartgrid(100, 100)
  p = traverse(grid, RandomPath())
  @test all(1 .≤ collect(p) .≤ 100 * 100)
  path = RandomPath(StableRNG(123))
  grid = cartgrid(3, 3)
  @test traverse(grid, path) == [4, 7, 2, 1, 3, 8, 5, 6, 9]

  grid = cartgrid(3, 3)
  pset = PointSet(centroid.(grid))
  for sdomain in [grid, pset]
    t = traverse(sdomain, SourcePath([1, 9]))
    @test collect(t) == [1, 9, 2, 4, 6, 8, 5, 3, 7]

    t = traverse(sdomain, SourcePath([1]))
    @test collect(t) == [1, 2, 4, 5, 3, 7, 6, 8, 9]
  end

  grid = cartgrid(3, 3)
  path = LinearPath()
  for offset in [0, 1, -1]
    spath = ShiftedPath(path, offset)
    t = traverse(grid, path)
    st = traverse(grid, spath)
    @test length(st) == 9
    @test collect(st) == circshift(t, -offset)
  end

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

  if visualtests
    paths = [LinearPath(), RandomPath(StableRNG(123)), ShiftedPath(LinearPath(), 10), SourcePath(1:3), MultiGridPath()]

    fnames = ["linear-path", "random-path", "shifted-path", "source-path", "multi-grid-path"]

    for (path, fname) in zip(paths, fnames)
      for d in (6, 7)
        agrid = cartgrid(d, d)
        elems = [agrid[i] for i in traverse(agrid, path)]
        fig = viz(elems, color=1:length(elems))
        @test_reference "data/$fname-$(d)x$(d).png" fig
      end
    end
  end
end
