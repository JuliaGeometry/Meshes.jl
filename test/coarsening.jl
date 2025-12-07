@testitem "RegularCoarsening" setup = [Setup] begin
  # 2D grids
  grid = CartesianGrid(cart(0.0, 0.0), cart(10.0, 10.0), dims=(20, 20))
  tgrid = CartesianGrid(cart(0.0, 0.0), cart(10.0, 10.0), dims=(10, 10))
  @test coarsen(grid, RegularCoarsening(2)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  crgrid = convert(RectilinearGrid, tgrid)
  @test coarsen(rgrid, RegularCoarsening(2)) == crgrid
  sgrid = convert(StructuredGrid, grid)
  csgrid = convert(StructuredGrid, tgrid)
  @test coarsen(sgrid, RegularCoarsening(2)) == csgrid
  tfgrid = TransformedGrid(grid, Identity())
  @test coarsen(tfgrid, RegularCoarsening(2)) == tgrid

  # 3D grids
  grid = cartgrid(100, 100, 100)
  tgrid = CartesianGrid(minimum(grid), maximum(grid), dims=(50, 25, 20))
  @test coarsen(grid, RegularCoarsening(2, 4, 5)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  crgrid = convert(RectilinearGrid, tgrid)
  @test coarsen(rgrid, RegularCoarsening(2, 4, 5)) == crgrid
  sgrid = convert(StructuredGrid, grid)
  csgrid = convert(StructuredGrid, tgrid)
  @test coarsen(sgrid, RegularCoarsening(2, 4, 5)) == csgrid
  tfgrid = TransformedGrid(grid, Identity())
  @test coarsen(tfgrid, RegularCoarsening(2, 4, 5)) == coarsen(grid, RegularCoarsening(2, 4, 5))

  # non-multiple dimensions (2D grids)
  grid = CartesianGrid(cart(0, 0), cart(13, 17), dims=(13, 17))
  tgrid = CartesianGrid(cart(0, 0), cart(13, 17), dims=(3, 6))
  @test coarsen(grid, RegularCoarsening(5, 3)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  @test size(coarsen(rgrid, RegularCoarsening(5, 3))) == (3, 6)
  sgrid = convert(StructuredGrid, grid)
  @test size(coarsen(sgrid, RegularCoarsening(5, 3))) == (3, 6)
  tfgrid = TransformedGrid(grid, Identity())
  @test size(coarsen(tfgrid, RegularCoarsening(5, 3))) == (3, 6)

  # non-multiple dimensions (3D grids)
  grid = CartesianGrid(cart(0, 0, 0), cart(13, 17, 23), dims=(13, 17, 23))
  tgrid = CartesianGrid(cart(0, 0, 0), cart(13, 17, 23), dims=(2, 4, 8))
  @test coarsen(grid, RegularCoarsening(7, 5, 3)) == tgrid
  rgrid = convert(RectilinearGrid, grid)
  @test size(coarsen(rgrid, RegularCoarsening(7, 5, 3))) == (2, 4, 8)
  sgrid = convert(StructuredGrid, grid)
  @test size(coarsen(sgrid, RegularCoarsening(7, 5, 3))) == (2, 4, 8)
  tfgrid = TransformedGrid(grid, Identity())
  @test size(coarsen(tfgrid, RegularCoarsening(7, 5, 3))) == (2, 4, 8)

  # preserve topology
  topo = GridTopology((100, 100), (true, false))
  ttopo = GridTopology((50, 50), (true, false))
  grid = CartesianGrid(cart(0, 0), T.((1, 1)), topo)
  @test topology(coarsen(grid, RegularCoarsening(2))) == ttopo
  rgrid = convert(RectilinearGrid, grid)
  @test topology(coarsen(rgrid, RegularCoarsening(2))) == ttopo
  sgrid = convert(StructuredGrid, grid)
  @test topology(coarsen(sgrid, RegularCoarsening(2))) == ttopo
  tfgrid = TransformedGrid(grid, Identity())
  @test topology(coarsen(tfgrid, RegularCoarsening(2))) == ttopo

  # large 2D grid
  grid = CartesianGrid(cart(0, 0), cart(16200, 8100), dims=(16200, 8100))
  tgrid = CartesianGrid(cart(0, 0), cart(16200, 8100), dims=(203, 203))
  @test coarsen(grid, RegularCoarsening(80, 40)) == tgrid

  # LatLon grid
  grid = RegularGrid(latlon(0, 0), latlon(45, 45), dims=(10, 10))
  tgrid = RegularGrid(latlon(0, 0), latlon(45, 45), dims=(5, 5))
  @test coarsen(grid, RegularCoarsening(2)) == tgrid
end
