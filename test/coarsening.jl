@testitem "RegularCoarsening" setup = [Setup] begin
  # 2D grids
  grid = CartesianGrid(cart(0.0, 0.0), cart(10.0, 10.0), dims=(20, 20))
  cgrid = CartesianGrid(cart(0.0, 0.0), cart(10.0, 10.0), dims=(10, 10))
  @test coarsen(grid, RegularCoarsening(2)) == cgrid
  rgrid = convert(RectilinearGrid, grid)
  crgrid = convert(RectilinearGrid, cgrid)
  @test coarsen(rgrid, RegularCoarsening(2)) == crgrid
  sgrid = convert(StructuredGrid, grid)
  csgrid = convert(StructuredGrid, cgrid)
  @test coarsen(sgrid, RegularCoarsening(2)) == csgrid
  tfgrid = TransformedGrid(grid, Identity())
  @test coarsen(tfgrid, RegularCoarsening(2)) == cgrid

  # 3D grids
  grid = cartgrid(100, 100, 100)
  cgrid = CartesianGrid(minimum(grid), maximum(grid), dims=(50, 25, 20))
  @test coarsen(grid, RegularCoarsening(2, 4, 5)) == cgrid
  rgrid = convert(RectilinearGrid, grid)
  crgrid = convert(RectilinearGrid, cgrid)
  @test coarsen(rgrid, RegularCoarsening(2, 4, 5)) == crgrid
  sgrid = convert(StructuredGrid, grid)
  csgrid = convert(StructuredGrid, cgrid)
  @test coarsen(sgrid, RegularCoarsening(2, 4, 5)) == csgrid
  tfgrid = TransformedGrid(grid, Identity())
  @test coarsen(tfgrid, RegularCoarsening(2, 4, 5)) == coarsen(grid, RegularCoarsening(2, 4, 5))

  # non-multiple dimensions (2D grids)
  grid = CartesianGrid(cart(0, 0), cart(13, 17), dims=(13, 17))
  cgrid = CartesianGrid(cart(0, 0), cart(13, 17), dims=(3, 6))
  @test coarsen(grid, RegularCoarsening(5, 3)) == cgrid
  rgrid = convert(RectilinearGrid, grid)
  @test size(coarsen(rgrid, RegularCoarsening(5, 3))) == (3, 6)
  sgrid = convert(StructuredGrid, grid)
  @test size(coarsen(sgrid, RegularCoarsening(5, 3))) == (3, 6)
  tfgrid = TransformedGrid(grid, Identity())
  @test size(coarsen(tfgrid, RegularCoarsening(5, 3))) == (3, 6)

  # non-multiple dimensions (3D grids)
  grid = CartesianGrid(cart(0, 0, 0), cart(13, 17, 23), dims=(13, 17, 23))
  cgrid = CartesianGrid(cart(0, 0, 0), cart(13, 17, 23), dims=(2, 4, 8))
  @test coarsen(grid, RegularCoarsening(7, 5, 3)) == cgrid
  rgrid = convert(RectilinearGrid, grid)
  @test size(coarsen(rgrid, RegularCoarsening(7, 5, 3))) == (2, 4, 8)
  sgrid = convert(StructuredGrid, grid)
  @test size(coarsen(sgrid, RegularCoarsening(7, 5, 3))) == (2, 4, 8)
  tfgrid = TransformedGrid(grid, Identity())
  @test size(coarsen(tfgrid, RegularCoarsening(7, 5, 3))) == (2, 4, 8)

  # large 2D grid
  grid = CartesianGrid(cart(0, 0), cart(16200, 8100), dims=(16200, 8100))
  cgrid = CartesianGrid(cart(0, 0), cart(16200, 8100), dims=(203, 203))
  @test coarsen(grid, RegularCoarsening(80, 40)) == cgrid

  # preserve topology
  topo = GridTopology((100, 100), (true, false))
  ctopo = GridTopology((50, 50), (true, false))
  grid = CartesianGrid(cart(0, 0), T.((1, 1)), topo)
  @test topology(coarsen(grid, RegularCoarsening(2))) == ctopo
  rgrid = convert(RectilinearGrid, grid)
  @test topology(coarsen(rgrid, RegularCoarsening(2))) == ctopo
  sgrid = convert(StructuredGrid, grid)
  @test topology(coarsen(sgrid, RegularCoarsening(2))) == ctopo
  tfgrid = TransformedGrid(grid, Identity())
  @test topology(coarsen(tfgrid, RegularCoarsening(2))) == ctopo
end
